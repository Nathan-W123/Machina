#include "sparlab/fem/Multigrid.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"
#include "sparlab/core/Timer.hpp"

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>
#include <Eigen/SparseCholesky>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <sstream>
#include <utility>

namespace sparlab {
namespace {

// Inner products sum fixed-size chunks in a fixed order, so their value does
// not depend on the number of threads. Short loops stay serial.
constexpr Index kChunk = 4096;
constexpr Index kParallelRows = 16384;

// ---------------------------------------------------------------------------
// Compressed sparse row storage and kernels
// ---------------------------------------------------------------------------

struct Csr {
  Index rows = 0;
  Index cols = 0;
  std::vector<Index> ptr;
  std::vector<Index> col;
  std::vector<Scalar> val;
  Index nnz() const { return ptr.empty() ? 0 : ptr.back(); }
};

/// Non-owning view: an owned Csr, or a symmetric Eigen matrix read as CSR
/// (the column-major storage of a symmetric matrix is its row-major storage).
struct CsrView {
  Index rows = 0;
  Index cols = 0;
  const Index* ptr = nullptr;
  const Index* col = nullptr;
  const Scalar* val = nullptr;
  Index nnz() const { return rows == 0 ? 0 : ptr[rows]; }
};

CsrView view_of(const Csr& m) {
  return {m.rows, m.cols, m.ptr.data(), m.col.data(), m.val.data()};
}

CsrView view_of_symmetric(const SparseMatrix& a) {
  return {static_cast<Index>(a.rows()), static_cast<Index>(a.cols()), a.outerIndexPtr(),
          a.innerIndexPtr(), a.valuePtr()};
}

SparseMatrix to_eigen(const CsrView& m) {
  // Build the transpose's column-major storage from the CSR arrays, then
  // transpose back: exact, and independent of symmetry.
  TripletList t;
  t.reserve(static_cast<std::size_t>(m.nnz()));
  for (Index i = 0; i < m.rows; ++i) {
    for (Index k = m.ptr[i]; k < m.ptr[i + 1]; ++k) t.emplace_back(i, m.col[k], m.val[k]);
  }
  SparseMatrix out(m.rows, m.cols);
  out.setFromTriplets(t.begin(), t.end());
  out.makeCompressed();
  return out;
}

Scalar chunked_dot(const Scalar* a, const Scalar* b, Index n) {
  const Index chunks = (n + kChunk - 1) / kChunk;
  if (chunks <= 1) {
    Scalar s = 0.0;
    for (Index i = 0; i < n; ++i) s += a[i] * b[i];
    return 0.0 + s;
  }
  std::vector<Scalar> partial(static_cast<std::size_t>(chunks), 0.0);
#if defined(SPARLAB_HAVE_OPENMP)
#pragma omp parallel for schedule(static) if (n >= kParallelRows)
#endif
  for (Index c = 0; c < chunks; ++c) {
    const Index begin = c * kChunk;
    const Index end = std::min(n, begin + kChunk);
    Scalar s = 0.0;
    for (Index i = begin; i < end; ++i) s += a[i] * b[i];
    partial[static_cast<std::size_t>(c)] = s;
  }
  Scalar total = 0.0;
  for (Scalar p : partial) total += p;
  return total;
}

/// y = A x
void spmv(const CsrView& a, const Scalar* x, Scalar* y) {
#if defined(SPARLAB_HAVE_OPENMP)
#pragma omp parallel for schedule(static) if (a.rows >= kParallelRows)
#endif
  for (Index i = 0; i < a.rows; ++i) {
    Scalar s = 0.0;
    for (Index k = a.ptr[i]; k < a.ptr[i + 1]; ++k) s += a.val[k] * x[a.col[k]];
    y[i] = s;
  }
}

/// r = b - A x
void residual(const CsrView& a, const Scalar* x, const Scalar* b, Scalar* r) {
#if defined(SPARLAB_HAVE_OPENMP)
#pragma omp parallel for schedule(static) if (a.rows >= kParallelRows)
#endif
  for (Index i = 0; i < a.rows; ++i) {
    Scalar s = 0.0;
    for (Index k = a.ptr[i]; k < a.ptr[i + 1]; ++k) s += a.val[k] * x[a.col[k]];
    r[i] = b[i] - s;
  }
}

/// C = A B (Gustavson), rows computed independently and sorted by column.
Csr multiply(const CsrView& a, const CsrView& b) {
  if (a.cols != b.rows) throw SolverError("multigrid: inner dimensions of a product differ");
  Csr c;
  c.rows = a.rows;
  c.cols = b.cols;
  c.ptr.assign(static_cast<std::size_t>(a.rows) + 1, 0);
#if defined(SPARLAB_HAVE_OPENMP)
#pragma omp parallel if (a.rows >= kParallelRows)
#endif
  {
    std::vector<Index> marker(static_cast<std::size_t>(b.cols), -1);
#if defined(SPARLAB_HAVE_OPENMP)
#pragma omp for schedule(static)
#endif
    for (Index i = 0; i < a.rows; ++i) {
      Index count = 0;
      for (Index k = a.ptr[i]; k < a.ptr[i + 1]; ++k) {
        const Index j = a.col[k];
        for (Index kk = b.ptr[j]; kk < b.ptr[j + 1]; ++kk) {
          const Index col = b.col[kk];
          if (marker[static_cast<std::size_t>(col)] != i) {
            marker[static_cast<std::size_t>(col)] = i;
            ++count;
          }
        }
      }
      c.ptr[static_cast<std::size_t>(i) + 1] = count;
    }
  }
  for (Index i = 0; i < a.rows; ++i) {
    c.ptr[static_cast<std::size_t>(i) + 1] += c.ptr[static_cast<std::size_t>(i)];
  }
  c.col.resize(static_cast<std::size_t>(c.nnz()));
  c.val.resize(static_cast<std::size_t>(c.nnz()));
#if defined(SPARLAB_HAVE_OPENMP)
#pragma omp parallel if (a.rows >= kParallelRows)
#endif
  {
    // Dense accumulator per thread; each output entry sums its products in
    // the order they are met, which does not depend on the thread count.
    std::vector<Index> marker(static_cast<std::size_t>(b.cols), -1);
    std::vector<Scalar> acc(static_cast<std::size_t>(b.cols), 0.0);
    std::vector<Index> cols;
#if defined(SPARLAB_HAVE_OPENMP)
#pragma omp for schedule(static)
#endif
    for (Index i = 0; i < a.rows; ++i) {
      cols.clear();
      for (Index k = a.ptr[i]; k < a.ptr[i + 1]; ++k) {
        const Index j = a.col[k];
        const Scalar av = a.val[k];
        for (Index kk = b.ptr[j]; kk < b.ptr[j + 1]; ++kk) {
          const std::size_t sc = static_cast<std::size_t>(b.col[kk]);
          if (marker[sc] != i) {
            marker[sc] = i;
            acc[sc] = av * b.val[kk];
            cols.push_back(b.col[kk]);
          } else {
            acc[sc] += av * b.val[kk];
          }
        }
      }
      std::sort(cols.begin(), cols.end());
      Index p = c.ptr[static_cast<std::size_t>(i)];
      for (Index col : cols) {
        c.col[static_cast<std::size_t>(p)] = col;
        c.val[static_cast<std::size_t>(p)] = acc[static_cast<std::size_t>(col)];
        ++p;
      }
    }
  }
  return c;
}

/// C = A B into the existing pattern of C (from an earlier `multiply` of
/// matrices with the same patterns). Products are summed in the same order as
/// `multiply`, so the values are identical to a full recomputation.
void multiply_numeric(const CsrView& a, const CsrView& b, Csr& c) {
#if defined(SPARLAB_HAVE_OPENMP)
#pragma omp parallel if (a.rows >= kParallelRows)
#endif
  {
    std::vector<Index> slot(static_cast<std::size_t>(b.cols), 0);
#if defined(SPARLAB_HAVE_OPENMP)
#pragma omp for schedule(static)
#endif
    for (Index i = 0; i < a.rows; ++i) {
      const Index begin = c.ptr[static_cast<std::size_t>(i)];
      const Index end = c.ptr[static_cast<std::size_t>(i) + 1];
      for (Index p = begin; p < end; ++p) {
        slot[static_cast<std::size_t>(c.col[static_cast<std::size_t>(p)])] = p;
        c.val[static_cast<std::size_t>(p)] = 0.0;
      }
      for (Index k = a.ptr[i]; k < a.ptr[i + 1]; ++k) {
        const Index j = a.col[k];
        const Scalar av = a.val[k];
        for (Index kk = b.ptr[j]; kk < b.ptr[j + 1]; ++kk) {
          c.val[static_cast<std::size_t>(slot[static_cast<std::size_t>(b.col[kk])])] +=
              av * b.val[kk];
        }
      }
    }
  }
}

Csr transpose(const CsrView& a) {
  Csr t;
  t.rows = a.cols;
  t.cols = a.rows;
  t.ptr.assign(static_cast<std::size_t>(a.cols) + 1, 0);
  for (Index k = 0; k < a.nnz(); ++k) ++t.ptr[static_cast<std::size_t>(a.col[k]) + 1];
  for (Index j = 0; j < a.cols; ++j) {
    t.ptr[static_cast<std::size_t>(j) + 1] += t.ptr[static_cast<std::size_t>(j)];
  }
  t.col.resize(static_cast<std::size_t>(a.nnz()));
  t.val.resize(static_cast<std::size_t>(a.nnz()));
  std::vector<Index> next(t.ptr.begin(), t.ptr.end() - 1);
  for (Index i = 0; i < a.rows; ++i) {
    for (Index k = a.ptr[i]; k < a.ptr[i + 1]; ++k) {
      const Index p = next[static_cast<std::size_t>(a.col[k])]++;
      t.col[static_cast<std::size_t>(p)] = i;
      t.val[static_cast<std::size_t>(p)] = a.val[k];
    }
  }
  return t;
}

/// (A + A^T) / 2 for a matrix whose rows are sorted by column.
Csr symmetrize(const Csr& a) {
  const Csr t = transpose(view_of(a));
  Csr s;
  s.rows = a.rows;
  s.cols = a.cols;
  s.ptr.assign(static_cast<std::size_t>(a.rows) + 1, 0);
  s.col.reserve(static_cast<std::size_t>(a.nnz()));
  s.val.reserve(static_cast<std::size_t>(a.nnz()));
  for (Index i = 0; i < a.rows; ++i) {
    Index p = a.ptr[static_cast<std::size_t>(i)];
    Index q = t.ptr[static_cast<std::size_t>(i)];
    const Index pe = a.ptr[static_cast<std::size_t>(i) + 1];
    const Index qe = t.ptr[static_cast<std::size_t>(i) + 1];
    while (p < pe || q < qe) {
      const Index cp = p < pe ? a.col[static_cast<std::size_t>(p)] : std::numeric_limits<Index>::max();
      const Index cq = q < qe ? t.col[static_cast<std::size_t>(q)] : std::numeric_limits<Index>::max();
      if (cp == cq) {
        s.col.push_back(cp);
        s.val.push_back(0.5 * (a.val[static_cast<std::size_t>(p)] + t.val[static_cast<std::size_t>(q)]));
        ++p;
        ++q;
      } else if (cp < cq) {
        s.col.push_back(cp);
        s.val.push_back(0.5 * a.val[static_cast<std::size_t>(p)]);
        ++p;
      } else {
        s.col.push_back(cq);
        s.val.push_back(0.5 * t.val[static_cast<std::size_t>(q)]);
        ++q;
      }
    }
    s.ptr[static_cast<std::size_t>(i) + 1] = static_cast<Index>(s.col.size());
  }
  return s;
}

/// Largest eigenvalue of D^-1 A (= of D^-1/2 A D^-1/2) by Lanczos. The Ritz
/// value is a lower bound that converges from below in a few steps.
Scalar lanczos_lambda_max(const CsrView& a, const std::vector<Scalar>& inv_diag, int steps) {
  const Index n = a.rows;
  if (n == 0) return 0.0;
  std::vector<Scalar> s(static_cast<std::size_t>(n));
  for (Index i = 0; i < n; ++i) s[static_cast<std::size_t>(i)] = std::sqrt(inv_diag[static_cast<std::size_t>(i)]);
  Vector v(n);
  for (Index i = 0; i < n; ++i) {
    v(i) = 1.0 + 0.5 * std::sin(1.7 * static_cast<Scalar>(i) + 0.3);
  }
  v /= std::sqrt(chunked_dot(v.data(), v.data(), n));
  Vector v_prev = Vector::Zero(n);
  Vector w(n);
  Vector tmp(n);
  std::vector<Scalar> alpha;
  std::vector<Scalar> beta;
  const int m = std::max(1, std::min<int>(steps, static_cast<int>(n)));
  for (int j = 0; j < m; ++j) {
    for (Index i = 0; i < n; ++i) tmp(i) = s[static_cast<std::size_t>(i)] * v(i);
    spmv(a, tmp.data(), w.data());
    for (Index i = 0; i < n; ++i) w(i) *= s[static_cast<std::size_t>(i)];
    if (j > 0) w -= beta.back() * v_prev;
    const Scalar aj = chunked_dot(w.data(), v.data(), n);
    alpha.push_back(aj);
    w -= aj * v;
    const Scalar bj = std::sqrt(chunked_dot(w.data(), w.data(), n));
    if (j + 1 == m || !(bj > 1.0e-12 * std::abs(aj))) break;
    beta.push_back(bj);
    v_prev = v;
    v = w / bj;
  }
  const Index k = static_cast<Index>(alpha.size());
  Vector diag(k);
  Vector sub(std::max<Index>(k - 1, 0));
  for (Index i = 0; i < k; ++i) diag(i) = alpha[static_cast<std::size_t>(i)];
  for (Index i = 0; i + 1 < k; ++i) sub(i) = beta[static_cast<std::size_t>(i)];
  if (k == 1) return diag(0);
  Eigen::SelfAdjointEigenSolver<Matrix> es;
  es.computeFromTridiagonal(diag, sub, Eigen::EigenvaluesOnly);
  return es.eigenvalues().maxCoeff();
}

}  // namespace

// ---------------------------------------------------------------------------
// Hierarchy
// ---------------------------------------------------------------------------

struct AmgPreconditioner::Impl {
  struct Level {
    Csr owned;                       // coarse levels own their operator
    CsrView a;                       // the operator (level 0: the caller's matrix)
    std::vector<Scalar> inv_diag;
    Scalar lambda_max = 0.0;
    // Unknowns grouped by node (contiguous ranges).
    Index num_nodes = 0;
    std::vector<Index> node_ptr;
    Matrix nullspace;                // unknowns x m
    // Transfer to the next level.
    std::vector<Index> aggregate;    // node -> aggregate
    Csr p_tent;
    Csr p;
    Csr r;                           // P^T
    Csr ap;                          // A P, kept for its pattern when reusing
    Csr coarse_raw;                  // P^T A P before symmetrisation, likewise
    // Workspace of the cycle.
    mutable Vector x, b, res, d;
  };

  AmgOptions options;
  std::vector<Level> levels;
  AmgStats stats;
  // Coarsest-level direct solve.
  Eigen::LDLT<Matrix> dense;
  Eigen::SimplicialLDLT<SparseMatrix, Eigen::Lower, Eigen::AMDOrdering<StorageIndex>> sparse;
  bool use_dense = true;
  bool ready = false;
  // Fingerprint of the structure the aggregates were built for.
  Index fp_rows = -1;
  Index fp_nnz = -1;
  const std::vector<Index>* fp_unknowns = nullptr;
  Index fp_unknowns_size = -1;
  const SparseMatrix* matrix = nullptr;

  void build_finest_nodes(const DofLayout& layout, Index n);
  void aggregate_level(Level& lvl);
  void tentative_prolongator(Level& lvl, Level& next);
  void numeric_level(std::size_t l);
  void smooth_prolongator(Level& lvl, bool numeric_only);
  void coarse_factorize(Level& lvl);
  void check_coarse_pivots(bool ok, const Vector& d, Index n);

  void smooth(const Level& lvl, bool pre, bool zero_guess) const;
  void cycle(std::size_t l) const;
};

void AmgPreconditioner::Impl::build_finest_nodes(const DofLayout& layout, Index n) {
  Level& lvl = levels.front();
  const std::vector<Index>& unknowns = *layout.unknowns;
  const int dim = layout.dim;
  const int ndpn = layout.nodal_dofs();
  const Matrix& coords = *layout.coordinates;
  lvl.node_ptr.clear();
  lvl.node_ptr.push_back(0);
  Index previous_node = -1;
  for (Index i = 0; i < n; ++i) {
    const Index g = unknowns[static_cast<std::size_t>(i)];
    const Index node = g / ndpn;
    if (node >= coords.cols() || g < 0) {
      std::ostringstream os;
      os << "multigrid: unknown " << i << " maps to DOF " << g
         << ", outside the model's " << coords.cols() << " nodes";
      throw SolverError(os.str());
    }
    if (i > 0 && g <= unknowns[static_cast<std::size_t>(i) - 1]) {
      throw SolverError("multigrid: the DOF layout must list unknowns in ascending order");
    }
    if (node != previous_node) {
      if (i > 0) lvl.node_ptr.push_back(i);
      previous_node = node;
    }
  }
  lvl.node_ptr.push_back(n);
  lvl.num_nodes = static_cast<Index>(lvl.node_ptr.size()) - 1;

  // Rigid-body modes about the centre of the mesh, in coordinates scaled to
  // unit size so translations and rotations are of comparable magnitude.
  Vector3 lo = Vector3::Constant(std::numeric_limits<Scalar>::max());
  Vector3 hi = Vector3::Constant(-std::numeric_limits<Scalar>::max());
  for (Index c = 0; c < coords.cols(); ++c) {
    Vector3 x = Vector3::Zero();
    x.head(dim) = coords.col(c);
    lo = lo.cwiseMin(x);
    hi = hi.cwiseMax(x);
  }
  const Vector3 centre = 0.5 * (lo + hi);
  const Scalar scale = std::max((hi - lo).maxCoeff(), 1.0e-300);
  const int m = dim == 2 ? 3 : 6;
  lvl.nullspace.setZero(n, m);
  for (Index i = 0; i < n; ++i) {
    const Index g = unknowns[static_cast<std::size_t>(i)];
    const Index node = g / ndpn;
    const int k = static_cast<int>(g % ndpn);
    Vector3 x = Vector3::Zero();
    x.head(dim) = coords.col(node);
    const Vector3 r = (x - centre) / scale;
    if (k >= 3) {
      // A nodal rotation of a shell or beam node: the rigid rotation about
      // axis k - 3 turns every node by the same angle. With the translations
      // written in coordinates scaled by 1/scale, the rotation entry of the
      // same mode is 1/scale, so each column stays a rigid motion.
      lvl.nullspace(i, k) = 1.0 / scale;
      continue;
    }
    lvl.nullspace(i, k) = 1.0;
    if (dim == 2) {
      lvl.nullspace(i, 2) = k == 0 ? -r.y() : r.x();
    } else {
      // Columns 3, 4, 5: rotations about x, y and z, u = omega x r.
      const Vector3 rx = Vector3::UnitX().cross(r);
      const Vector3 ry = Vector3::UnitY().cross(r);
      const Vector3 rz = Vector3::UnitZ().cross(r);
      lvl.nullspace(i, 3) = rx(k);
      lvl.nullspace(i, 4) = ry(k);
      lvl.nullspace(i, 5) = rz(k);
    }
  }
  stats.near_null_space_dimension = m;
}

void AmgPreconditioner::Impl::aggregate_level(Level& lvl) {
  const CsrView& a = lvl.a;
  const Index nn = lvl.num_nodes;
  std::vector<Index> node_of(static_cast<std::size_t>(a.rows));
  for (Index node = 0; node < nn; ++node) {
    for (Index i = lvl.node_ptr[static_cast<std::size_t>(node)];
         i < lvl.node_ptr[static_cast<std::size_t>(node) + 1]; ++i) {
      node_of[static_cast<std::size_t>(i)] = node;
    }
  }
  // Frobenius norms of the node blocks: diagonal first, then the graph.
  std::vector<Scalar> block_diag(static_cast<std::size_t>(nn), 0.0);
  for (Index node = 0; node < nn; ++node) {
    Scalar s = 0.0;
    for (Index i = lvl.node_ptr[static_cast<std::size_t>(node)];
         i < lvl.node_ptr[static_cast<std::size_t>(node) + 1]; ++i) {
      for (Index k = a.ptr[i]; k < a.ptr[i + 1]; ++k) {
        if (node_of[static_cast<std::size_t>(a.col[k])] == node) s += a.val[k] * a.val[k];
      }
    }
    block_diag[static_cast<std::size_t>(node)] = std::sqrt(s);
  }
  std::vector<Index> strong_ptr(static_cast<std::size_t>(nn) + 1, 0);
  std::vector<Index> strong;
  std::vector<Scalar> strong_value;
  strong.reserve(static_cast<std::size_t>(a.nnz()) / 4 + 1);
  strong_value.reserve(strong.capacity());
  {
    std::vector<Index> marker(static_cast<std::size_t>(nn), -1);
    std::vector<Scalar> acc(static_cast<std::size_t>(nn), 0.0);
    std::vector<Index> touched;
    const Scalar theta = options.strength_threshold;
    for (Index node = 0; node < nn; ++node) {
      touched.clear();
      for (Index i = lvl.node_ptr[static_cast<std::size_t>(node)];
           i < lvl.node_ptr[static_cast<std::size_t>(node) + 1]; ++i) {
        for (Index k = a.ptr[i]; k < a.ptr[i + 1]; ++k) {
          const Index other = node_of[static_cast<std::size_t>(a.col[k])];
          if (other == node) continue;
          if (marker[static_cast<std::size_t>(other)] != node) {
            marker[static_cast<std::size_t>(other)] = node;
            acc[static_cast<std::size_t>(other)] = 0.0;
            touched.push_back(other);
          }
          acc[static_cast<std::size_t>(other)] += a.val[k] * a.val[k];
        }
      }
      std::sort(touched.begin(), touched.end());
      const Scalar di = block_diag[static_cast<std::size_t>(node)];
      for (Index other : touched) {
        const Scalar s = std::sqrt(acc[static_cast<std::size_t>(other)]);
        const Scalar ref = std::sqrt(di * block_diag[static_cast<std::size_t>(other)]);
        if (s > 0.0 && s >= theta * ref) {
          strong.push_back(other);
          strong_value.push_back(ref > 0.0 ? s / ref : 0.0);
        }
      }
      strong_ptr[static_cast<std::size_t>(node) + 1] = static_cast<Index>(strong.size());
    }
  }

  // Greedy aggregation (Vanek, Mandel, Brezina), three passes in node order.
  std::vector<Index>& agg = lvl.aggregate;
  agg.assign(static_cast<std::size_t>(nn), -1);
  Index count = 0;
  for (Index node = 0; node < nn; ++node) {
    if (agg[static_cast<std::size_t>(node)] >= 0) continue;
    const Index b = strong_ptr[static_cast<std::size_t>(node)];
    const Index e = strong_ptr[static_cast<std::size_t>(node) + 1];
    if (b == e) continue;
    bool free_neighbourhood = true;
    for (Index k = b; k < e; ++k) {
      if (agg[static_cast<std::size_t>(strong[static_cast<std::size_t>(k)])] >= 0) {
        free_neighbourhood = false;
        break;
      }
    }
    if (!free_neighbourhood) continue;
    agg[static_cast<std::size_t>(node)] = count;
    for (Index k = b; k < e; ++k) agg[static_cast<std::size_t>(strong[static_cast<std::size_t>(k)])] = count;
    ++count;
  }
  const std::vector<Index> first_pass = agg;
  for (Index node = 0; node < nn; ++node) {
    if (first_pass[static_cast<std::size_t>(node)] >= 0) continue;
    Index best = -1;
    Scalar best_value = -1.0;
    for (Index k = strong_ptr[static_cast<std::size_t>(node)];
         k < strong_ptr[static_cast<std::size_t>(node) + 1]; ++k) {
      const Index other = strong[static_cast<std::size_t>(k)];
      const Index target = first_pass[static_cast<std::size_t>(other)];
      if (target >= 0 && strong_value[static_cast<std::size_t>(k)] > best_value) {
        best = target;
        best_value = strong_value[static_cast<std::size_t>(k)];
      }
    }
    if (best >= 0) agg[static_cast<std::size_t>(node)] = best;
  }
  for (Index node = 0; node < nn; ++node) {
    if (agg[static_cast<std::size_t>(node)] >= 0) continue;
    agg[static_cast<std::size_t>(node)] = count;
    for (Index k = strong_ptr[static_cast<std::size_t>(node)];
         k < strong_ptr[static_cast<std::size_t>(node) + 1]; ++k) {
      const Index other = strong[static_cast<std::size_t>(k)];
      if (agg[static_cast<std::size_t>(other)] < 0) agg[static_cast<std::size_t>(other)] = count;
    }
    ++count;
  }
}

void AmgPreconditioner::Impl::tentative_prolongator(Level& lvl, Level& next) {
  const Index n = lvl.a.rows;
  const Index nn = lvl.num_nodes;
  const Index m = static_cast<Index>(lvl.nullspace.cols());
  const std::vector<Index>& agg = lvl.aggregate;
  const Index na = agg.empty() ? 0 : *std::max_element(agg.begin(), agg.end()) + 1;

  // Nodes of each aggregate, ascending.
  std::vector<Index> agg_ptr(static_cast<std::size_t>(na) + 1, 0);
  for (Index node = 0; node < nn; ++node) ++agg_ptr[static_cast<std::size_t>(agg[static_cast<std::size_t>(node)]) + 1];
  for (Index a = 0; a < na; ++a) agg_ptr[static_cast<std::size_t>(a) + 1] += agg_ptr[static_cast<std::size_t>(a)];
  std::vector<Index> agg_nodes(static_cast<std::size_t>(nn));
  {
    std::vector<Index> fill(agg_ptr.begin(), agg_ptr.end() - 1);
    for (Index node = 0; node < nn; ++node) {
      agg_nodes[static_cast<std::size_t>(fill[static_cast<std::size_t>(agg[static_cast<std::size_t>(node)])]++)] = node;
    }
  }

  std::vector<Matrix> q(static_cast<std::size_t>(na));
  std::vector<Matrix> r(static_cast<std::size_t>(na));
#if defined(SPARLAB_HAVE_OPENMP)
#pragma omp parallel for schedule(dynamic, 256) if (na >= 4096)
#endif
  for (Index a = 0; a < na; ++a) {
    Index rows = 0;
    for (Index k = agg_ptr[static_cast<std::size_t>(a)]; k < agg_ptr[static_cast<std::size_t>(a) + 1]; ++k) {
      const Index node = agg_nodes[static_cast<std::size_t>(k)];
      rows += lvl.node_ptr[static_cast<std::size_t>(node) + 1] - lvl.node_ptr[static_cast<std::size_t>(node)];
    }
    Matrix block(rows, m);
    Index row = 0;
    for (Index k = agg_ptr[static_cast<std::size_t>(a)]; k < agg_ptr[static_cast<std::size_t>(a) + 1]; ++k) {
      const Index node = agg_nodes[static_cast<std::size_t>(k)];
      for (Index i = lvl.node_ptr[static_cast<std::size_t>(node)];
           i < lvl.node_ptr[static_cast<std::size_t>(node) + 1]; ++i) {
        block.row(row++) = lvl.nullspace.row(i);
      }
    }
    Eigen::ColPivHouseholderQR<Matrix> qr(block);
    const Matrix rfull = qr.matrixR().template triangularView<Eigen::Upper>();
    const Index kmax = std::min(rows, m);
    const Scalar lead = kmax > 0 ? std::abs(rfull(0, 0)) : 0.0;
    Index rank = 0;
    for (Index i = 0; i < kmax; ++i) {
      if (std::abs(rfull(i, i)) > 1.0e-10 * lead) ++rank;
    }
    rank = std::max<Index>(rank, 1);
    q[static_cast<std::size_t>(a)] = qr.householderQ() * Matrix::Identity(rows, rank);
    r[static_cast<std::size_t>(a)] =
        rfull.topRows(rank) * qr.colsPermutation().transpose();
  }

  // Coarse unknowns: aggregate a owns [offset_a, offset_a + rank_a).
  next.node_ptr.assign(static_cast<std::size_t>(na) + 1, 0);
  for (Index a = 0; a < na; ++a) {
    next.node_ptr[static_cast<std::size_t>(a) + 1] =
        next.node_ptr[static_cast<std::size_t>(a)] + static_cast<Index>(q[static_cast<std::size_t>(a)].cols());
  }
  next.num_nodes = na;
  const Index nc = next.node_ptr.back();
  next.nullspace.setZero(nc, m);
  for (Index a = 0; a < na; ++a) {
    next.nullspace.middleRows(next.node_ptr[static_cast<std::size_t>(a)], r[static_cast<std::size_t>(a)].rows()) =
        r[static_cast<std::size_t>(a)];
  }

  // P-hat: row i (in aggregate a, local position p) holds Q_a(p, :).
  Csr& pt = lvl.p_tent;
  pt.rows = n;
  pt.cols = nc;
  pt.ptr.assign(static_cast<std::size_t>(n) + 1, 0);
  std::vector<Index> agg_of_row(static_cast<std::size_t>(n));
  std::vector<Index> local_row(static_cast<std::size_t>(n));
  for (Index a = 0; a < na; ++a) {
    Index local = 0;
    for (Index k = agg_ptr[static_cast<std::size_t>(a)]; k < agg_ptr[static_cast<std::size_t>(a) + 1]; ++k) {
      const Index node = agg_nodes[static_cast<std::size_t>(k)];
      for (Index i = lvl.node_ptr[static_cast<std::size_t>(node)];
           i < lvl.node_ptr[static_cast<std::size_t>(node) + 1]; ++i) {
        agg_of_row[static_cast<std::size_t>(i)] = a;
        local_row[static_cast<std::size_t>(i)] = local++;
      }
    }
  }
  for (Index i = 0; i < n; ++i) {
    pt.ptr[static_cast<std::size_t>(i) + 1] =
        pt.ptr[static_cast<std::size_t>(i)] +
        static_cast<Index>(q[static_cast<std::size_t>(agg_of_row[static_cast<std::size_t>(i)])].cols());
  }
  pt.col.resize(static_cast<std::size_t>(pt.nnz()));
  pt.val.resize(static_cast<std::size_t>(pt.nnz()));
  for (Index i = 0; i < n; ++i) {
    const Index a = agg_of_row[static_cast<std::size_t>(i)];
    const Matrix& qa = q[static_cast<std::size_t>(a)];
    Index p = pt.ptr[static_cast<std::size_t>(i)];
    for (Index j = 0; j < qa.cols(); ++j) {
      pt.col[static_cast<std::size_t>(p)] = next.node_ptr[static_cast<std::size_t>(a)] + j;
      pt.val[static_cast<std::size_t>(p)] = qa(local_row[static_cast<std::size_t>(i)], j);
      ++p;
    }
  }
}

void AmgPreconditioner::Impl::smooth_prolongator(Level& lvl, bool numeric_only) {
  // P = (I - omega D^-1 A) P-hat. Every row of P-hat is contained in the
  // pattern of A P-hat (A has a non-zero diagonal), so P takes that pattern.
  const Scalar omega = options.prolongator_damping / lvl.lambda_max;
  Csr ap;
  if (numeric_only) {
    ap = std::move(lvl.p);
    multiply_numeric(lvl.a, view_of(lvl.p_tent), ap);
  } else {
    ap = multiply(lvl.a, view_of(lvl.p_tent));
  }
  const Csr& pt = lvl.p_tent;
#if defined(SPARLAB_HAVE_OPENMP)
#pragma omp parallel for schedule(static) if (ap.rows >= kParallelRows)
#endif
  for (Index i = 0; i < ap.rows; ++i) {
    const Scalar scale = omega * lvl.inv_diag[static_cast<std::size_t>(i)];
    for (Index k = ap.ptr[static_cast<std::size_t>(i)]; k < ap.ptr[static_cast<std::size_t>(i) + 1]; ++k) {
      ap.val[static_cast<std::size_t>(k)] *= -scale;
    }
    for (Index kk = pt.ptr[static_cast<std::size_t>(i)]; kk < pt.ptr[static_cast<std::size_t>(i) + 1]; ++kk) {
      const Index c = pt.col[static_cast<std::size_t>(kk)];
      const Index* begin = ap.col.data() + ap.ptr[static_cast<std::size_t>(i)];
      const Index* end = ap.col.data() + ap.ptr[static_cast<std::size_t>(i) + 1];
      const Index* hit = std::lower_bound(begin, end, c);
      if (hit == end || *hit != c) {
        throw SolverError("multigrid: internal error, the smoothed prolongator lost an entry");
      }
      ap.val[static_cast<std::size_t>(hit - ap.col.data())] += pt.val[static_cast<std::size_t>(kk)];
    }
  }
  lvl.p = std::move(ap);
  if (numeric_only) {
    // Same pattern: refill the transpose's values in its fixed order.
    Csr& r = lvl.r;
    std::vector<Index> next(r.ptr.begin(), r.ptr.end() - 1);
    for (Index i = 0; i < lvl.p.rows; ++i) {
      for (Index k = lvl.p.ptr[static_cast<std::size_t>(i)]; k < lvl.p.ptr[static_cast<std::size_t>(i) + 1]; ++k) {
        r.val[static_cast<std::size_t>(next[static_cast<std::size_t>(lvl.p.col[static_cast<std::size_t>(k)])]++)] =
            lvl.p.val[static_cast<std::size_t>(k)];
      }
    }
  } else {
    lvl.r = transpose(view_of(lvl.p));
  }
}

void AmgPreconditioner::Impl::numeric_level(std::size_t l) {
  Level& lvl = levels[l];
  const CsrView& a = lvl.a;
  lvl.inv_diag.assign(static_cast<std::size_t>(a.rows), 0.0);
  for (Index i = 0; i < a.rows; ++i) {
    Scalar d = 0.0;
    for (Index k = a.ptr[i]; k < a.ptr[i + 1]; ++k) {
      if (a.col[k] == i) d += a.val[k];
    }
    if (!(d > 0.0)) {
      std::ostringstream os;
      os << "multigrid level " << l << ": diagonal entry " << i << " is " << d
         << "; a stiffness matrix must have a positive diagonal. The model has an "
            "unconnected degree of freedom or a non-positive stiffness";
      throw SolverError(os.str());
    }
    lvl.inv_diag[static_cast<std::size_t>(i)] = 1.0 / d;
  }
  lvl.lambda_max = lanczos_lambda_max(a, lvl.inv_diag, options.lanczos_steps);
  if (!(lvl.lambda_max > 0.0) || !std::isfinite(lvl.lambda_max)) {
    std::ostringstream os;
    os << "multigrid level " << l << ": the spectral estimate of D^-1 A is "
       << lvl.lambda_max << "; the matrix is not positive definite";
    throw SolverError(os.str());
  }
  lvl.x.setZero(a.rows);
  lvl.b.setZero(a.rows);
  lvl.res.setZero(a.rows);
  lvl.d.setZero(a.rows);
}

void AmgPreconditioner::Impl::coarse_factorize(Level& lvl) {
  const Index n = lvl.a.rows;
  use_dense = n <= 4000;
  if (use_dense) {
    Matrix dense_a = Matrix::Zero(n, n);
    for (Index i = 0; i < n; ++i) {
      for (Index k = lvl.a.ptr[i]; k < lvl.a.ptr[i + 1]; ++k) dense_a(i, lvl.a.col[k]) = lvl.a.val[k];
    }
    dense.compute(dense_a);
    check_coarse_pivots(dense.info() == Eigen::Success, dense.vectorD(), n);
    return;
  }
  const SparseMatrix a = to_eigen(lvl.a);
  sparse.compute(a);
  check_coarse_pivots(sparse.info() == Eigen::Success, sparse.vectorD(), n);
}

void AmgPreconditioner::Impl::check_coarse_pivots(bool ok, const Vector& d, Index n) {
  // The coarse space contains every rigid-body motion exactly, so an
  // unresisted motion of the model shows up here as a zero (round-off sized)
  // pivot however large the model is.
  const Scalar max_abs = d.size() > 0 ? d.cwiseAbs().maxCoeff() : 0.0;
  const Scalar ratio = max_abs > 0.0 ? d.minCoeff() / max_abs : 0.0;
  stats.coarse_pivot_ratio = ratio;
  if (!ok || !(ratio > options.coarse_pivot_tolerance)) {
    std::ostringstream os;
    os << "multigrid: the coarsest operator (" << n << " unknowns) is singular or "
       << "indefinite (smallest / largest pivot " << ratio
       << "). A rigid-body mode survived to the coarse level, so the model is "
          "under-constrained: too few supports, or a part connected to nothing";
    throw SolverError(os.str());
  }
}

void AmgPreconditioner::Impl::smooth(const Level& lvl, bool pre, bool zero_guess) const {
  const CsrView& a = lvl.a;
  const Index n = a.rows;
  Vector& x = lvl.x;
  const Vector& b = lvl.b;
  if (options.smoother == AmgSmoother::SymmetricGaussSeidel) {
    if (zero_guess) x.setZero();
    for (int sweep = 0; sweep < options.smoother_degree; ++sweep) {
      if (pre) {
        for (Index i = 0; i < n; ++i) {
          Scalar s = b(i);
          for (Index k = a.ptr[i]; k < a.ptr[i + 1]; ++k) s -= a.val[k] * x(a.col[k]);
          x(i) += s * lvl.inv_diag[static_cast<std::size_t>(i)];
        }
      } else {
        for (Index i = n - 1; i >= 0; --i) {
          Scalar s = b(i);
          for (Index k = a.ptr[i]; k < a.ptr[i + 1]; ++k) s -= a.val[k] * x(a.col[k]);
          x(i) += s * lvl.inv_diag[static_cast<std::size_t>(i)];
        }
      }
    }
    return;
  }
  // Chebyshev (Saad, Algorithm 12.1, preconditioned by D) on
  // [lambda_max / ratio, 1.1 lambda_max] of D^-1 A.
  const Scalar upper = 1.1 * lvl.lambda_max;
  const Scalar lower = upper / options.chebyshev_ratio;
  const Scalar theta = 0.5 * (upper + lower);
  const Scalar delta = 0.5 * (upper - lower);
  const Scalar sigma = theta / delta;
  Scalar rho = 1.0 / sigma;
  Vector& r = lvl.res;
  Vector& d = lvl.d;
  if (zero_guess) {
    x.setZero();
    r = b;
  } else {
    residual(a, x.data(), b.data(), r.data());
  }
  const std::vector<Scalar>& dinv = lvl.inv_diag;
  const bool parallel = n >= kParallelRows;
  (void)parallel;
#if defined(SPARLAB_HAVE_OPENMP)
#pragma omp parallel for schedule(static) if (parallel)
#endif
  for (Index i = 0; i < n; ++i) {
    d(i) = dinv[static_cast<std::size_t>(i)] * r(i) / theta;
    x(i) += d(i);
  }
  for (int k = 1; k < options.smoother_degree; ++k) {
    // r <- r - A d  (the residual of the updated x)
#if defined(SPARLAB_HAVE_OPENMP)
#pragma omp parallel for schedule(static) if (parallel)
#endif
    for (Index i = 0; i < n; ++i) {
      Scalar s = 0.0;
      for (Index kk = a.ptr[i]; kk < a.ptr[i + 1]; ++kk) s += a.val[kk] * d(a.col[kk]);
      r(i) -= s;
    }
    const Scalar rho_new = 1.0 / (2.0 * sigma - rho);
    const Scalar c1 = rho_new * rho;
    const Scalar c2 = 2.0 * rho_new / delta;
#if defined(SPARLAB_HAVE_OPENMP)
#pragma omp parallel for schedule(static) if (parallel)
#endif
    for (Index i = 0; i < n; ++i) {
      d(i) = c1 * d(i) + c2 * dinv[static_cast<std::size_t>(i)] * r(i);
      x(i) += d(i);
    }
    rho = rho_new;
  }
}

void AmgPreconditioner::Impl::cycle(std::size_t l) const {
  const Level& lvl = levels[l];
  if (l + 1 == levels.size()) {
    if (use_dense) {
      lvl.x = dense.solve(lvl.b);
    } else {
      lvl.x = sparse.solve(lvl.b);
    }
    return;
  }
  smooth(lvl, /*pre=*/true, /*zero_guess=*/true);
  residual(lvl.a, lvl.x.data(), lvl.b.data(), lvl.res.data());
  const Level& next = levels[l + 1];
  spmv(view_of(lvl.r), lvl.res.data(), next.b.data());
  cycle(l + 1);
  // x += P x_coarse
  const Csr& p = lvl.p;
#if defined(SPARLAB_HAVE_OPENMP)
#pragma omp parallel for schedule(static) if (p.rows >= kParallelRows)
#endif
  for (Index i = 0; i < p.rows; ++i) {
    Scalar s = 0.0;
    for (Index k = p.ptr[static_cast<std::size_t>(i)]; k < p.ptr[static_cast<std::size_t>(i) + 1]; ++k) {
      s += p.val[static_cast<std::size_t>(k)] * next.x(p.col[static_cast<std::size_t>(k)]);
    }
    lvl.x(i) += s;
  }
  smooth(lvl, /*pre=*/false, /*zero_guess=*/false);
}

// ---------------------------------------------------------------------------
// Public interface
// ---------------------------------------------------------------------------

std::string to_string(AmgSmoother smoother) {
  switch (smoother) {
    case AmgSmoother::Chebyshev: return "chebyshev";
    case AmgSmoother::SymmetricGaussSeidel: return "gauss_seidel";
  }
  return "unknown";
}

AmgSmoother parse_amg_smoother(const std::string& text) {
  if (text == "chebyshev") return AmgSmoother::Chebyshev;
  if (text == "gauss_seidel" || text == "symmetric_gauss_seidel") {
    return AmgSmoother::SymmetricGaussSeidel;
  }
  throw ConfigError("unknown multigrid smoother '" + text +
                    "' (expected chebyshev|gauss_seidel)");
}

AmgPreconditioner::AmgPreconditioner(AmgOptions options) : impl_(std::make_unique<Impl>()) {
  impl_->options = options;
  if (!(options.strength_threshold >= 0.0 && options.strength_threshold < 1.0)) {
    throw ConfigError("multigrid strength_threshold must lie in [0, 1)");
  }
  if (options.max_levels < 1) throw ConfigError("multigrid max_levels must be at least 1");
  if (options.coarse_size < 1) throw ConfigError("multigrid coarse_size must be positive");
  if (options.smoother_degree < 1) {
    throw ConfigError("multigrid smoother_degree must be at least 1");
  }
  if (!(options.chebyshev_ratio > 1.0)) {
    throw ConfigError("multigrid chebyshev_ratio must exceed 1");
  }
  if (!(options.prolongator_damping > 0.0 && options.prolongator_damping < 2.0)) {
    throw ConfigError("multigrid prolongator_damping must lie in (0, 2)");
  }
  if (options.lanczos_steps < 2) throw ConfigError("multigrid lanczos_steps must be at least 2");
}

AmgPreconditioner::~AmgPreconditioner() = default;

void AmgPreconditioner::setup(const SparseMatrix& a, const DofLayout& layout) {
  Impl& im = *impl_;
  Timer timer;
  if (a.rows() != a.cols() || a.rows() == 0) {
    throw SolverError("multigrid: the system matrix must be square and non-empty");
  }
  if (!a.isCompressed()) {
    throw SolverError("multigrid: the system matrix must be in compressed storage");
  }
  if (layout.dim != 2 && layout.dim != 3) {
    throw SolverError("multigrid: the DOF layout must be two- or three-dimensional");
  }
  if (layout.nodal_dofs() != layout.dim &&
      !(layout.dim == 3 && layout.nodal_dofs() == kMaxDofsPerNode)) {
    std::ostringstream os;
    os << "multigrid: a " << layout.dim << "-D layout with " << layout.nodal_dofs()
       << " DOFs per node is not supported (translations only, or three translations "
          "and three rotations in 3-D)";
    throw SolverError(os.str());
  }
  if (layout.coordinates == nullptr || layout.unknowns == nullptr ||
      static_cast<Index>(layout.unknowns->size()) != a.rows()) {
    std::ostringstream os;
    os << "multigrid: the DOF layout lists "
       << (layout.unknowns ? layout.unknowns->size() : 0) << " unknowns for a matrix of "
       << a.rows() << " rows";
    throw SolverError(os.str());
  }
  const Index n = static_cast<Index>(a.rows());
  AmgStats& st = im.stats;
  st.check_seconds = st.spectral_seconds = st.aggregation_seconds = 0.0;
  st.tentative_seconds = st.smoothing_seconds = st.galerkin_seconds = st.coarse_seconds = 0.0;
  Timer phase;
  const bool same_structure = im.ready && im.fp_rows == n &&
                              im.fp_nnz == static_cast<Index>(a.nonZeros()) &&
                              im.fp_unknowns == layout.unknowns &&
                              im.fp_unknowns_size == static_cast<Index>(layout.unknowns->size());
  const bool reuse = im.options.reuse_aggregates && same_structure;

  if (!same_structure) {
    // Symmetry check (once per structure): every stored (i, j) has a (j, i)
    // with the same value up to round-off.
    const CsrView v = view_of_symmetric(a);
    Scalar max_abs = 0.0;
    for (Index k = 0; k < v.nnz(); ++k) max_abs = std::max(max_abs, std::abs(v.val[k]));
    for (Index i = 0; i < n; ++i) {
      for (Index k = v.ptr[i]; k < v.ptr[i + 1]; ++k) {
        const Index j = v.col[k];
        const Index* begin = v.col + v.ptr[j];
        const Index* end = v.col + v.ptr[j + 1];
        const Index* hit = std::lower_bound(begin, end, i);
        const bool found = hit != end && *hit == i;
        const Scalar other = found ? v.val[hit - v.col] : 0.0;
        if (std::abs(other - v.val[k]) > 1.0e-10 * max_abs) {
          std::ostringstream os;
          os << "multigrid: the system matrix is not symmetric (entry (" << i << ", " << j
             << ") = " << v.val[k] << " but (" << j << ", " << i << ") = " << other
             << "); conjugate gradients needs a symmetric positive-definite matrix";
          throw SolverError(os.str());
        }
      }
    }
  }
  for (Index k = 0; k < static_cast<Index>(a.nonZeros()); ++k) {
    if (!std::isfinite(a.valuePtr()[k])) {
      throw SolverError("multigrid: the system matrix has a non-finite entry");
    }
  }

  st.check_seconds = phase.elapsed_seconds();
  if (!reuse) {
    im.levels.clear();
    // Level operators are viewed through raw pointers into their neighbours'
    // storage, so the level array must never reallocate.
    im.levels.reserve(static_cast<std::size_t>(im.options.max_levels) + 1);
    im.levels.emplace_back();
    im.build_finest_nodes(layout, n);
  }
  im.levels.front().a = view_of_symmetric(a);
  im.matrix = &a;

  std::size_t l = 0;
  const std::size_t reused_levels = im.levels.size();
  while (true) {
    Impl::Level& lvl = im.levels[l];
    phase.reset();
    im.numeric_level(l);
    st.spectral_seconds += phase.elapsed_seconds();
    const bool last_reused = reuse && l + 1 == reused_levels;
    const bool stop = reuse ? last_reused
                            : (lvl.a.rows <= im.options.coarse_size ||
                               static_cast<int>(l) + 1 >= im.options.max_levels);
    if (stop) break;
    if (!reuse) {
      phase.reset();
      im.aggregate_level(lvl);
      st.aggregation_seconds += phase.elapsed_seconds();
      phase.reset();
      Impl::Level next;
      im.tentative_prolongator(lvl, next);
      st.tentative_seconds += phase.elapsed_seconds();
      const Index nc = next.node_ptr.back();
      if (nc >= static_cast<Index>(0.85 * static_cast<Scalar>(lvl.a.rows)) || nc == 0) {
        // Aggregation no longer reduces the problem: solve this level directly.
        lvl.aggregate.clear();
        lvl.p_tent = Csr();
        log::debug("multigrid: coarsening stalled at level ", l, " (", lvl.a.rows,
                   " -> ", nc, " unknowns); solving it directly");
        break;
      }
      im.levels.push_back(std::move(next));
    }
    Impl::Level& cur = im.levels[l];
    phase.reset();
    im.smooth_prolongator(cur, reuse);
    st.smoothing_seconds += phase.elapsed_seconds();
    phase.reset();
    Impl::Level& coarse = im.levels[l + 1];
    if (reuse) {
      multiply_numeric(cur.a, view_of(cur.p), cur.ap);
      multiply_numeric(view_of(cur.r), view_of(cur.ap), cur.coarse_raw);
    } else {
      cur.ap = multiply(cur.a, view_of(cur.p));
      cur.coarse_raw = multiply(view_of(cur.r), view_of(cur.ap));
    }
    coarse.owned = symmetrize(cur.coarse_raw);
    coarse.a = view_of(coarse.owned);
    st.galerkin_seconds += phase.elapsed_seconds();
    log::trace("multigrid level ", l, ": A ", cur.a.rows, " rows / ", cur.a.nnz(),
               " nnz, P ", cur.p.nnz(), " nnz, AP ", cur.ap.nnz(), " nnz, coarse ",
               coarse.a.rows, " rows / ", coarse.a.nnz(), " nnz");
    ++l;
  }
  im.levels.resize(l + 1);
  phase.reset();
  im.coarse_factorize(im.levels.back());
  st.coarse_seconds = phase.elapsed_seconds();

  im.fp_rows = n;
  im.fp_nnz = static_cast<Index>(a.nonZeros());
  im.fp_unknowns = layout.unknowns;
  im.fp_unknowns_size = static_cast<Index>(layout.unknowns->size());
  im.ready = true;

  AmgStats& s = st;
  s.levels.clear();
  Scalar nnz_sum = 0.0;
  Scalar rows_sum = 0.0;
  for (std::size_t k = 0; k < im.levels.size(); ++k) {
    const Impl::Level& lvl = im.levels[k];
    AmgLevelStats ls;
    ls.unknowns = lvl.a.rows;
    ls.nonzeros = lvl.a.nnz();
    ls.aggregates = k + 1 < im.levels.size() ? im.levels[k + 1].num_nodes : 0;
    ls.lambda_max = lvl.lambda_max;
    s.levels.push_back(ls);
    nnz_sum += static_cast<Scalar>(ls.nonzeros);
    rows_sum += static_cast<Scalar>(ls.unknowns);
  }
  s.operator_complexity = nnz_sum / static_cast<Scalar>(s.levels.front().nonzeros);
  {
    Scalar stored = 0.0;
    for (std::size_t k = 0; k < im.levels.size(); ++k) {
      const Impl::Level& lvl = im.levels[k];
      if (k > 0) stored += static_cast<Scalar>(lvl.a.nnz());
      stored += static_cast<Scalar>(lvl.p.nnz() + lvl.r.nnz() + lvl.ap.nnz() +
                                    lvl.coarse_raw.nnz() + lvl.p_tent.nnz());
    }
    const Scalar nc = static_cast<Scalar>(im.levels.back().a.rows);
    stored += im.use_dense ? nc * nc : 0.0;
    s.storage_nonzeros = static_cast<Index>(std::min(stored, 2.0e9));
  }
  s.grid_complexity = rows_sum / static_cast<Scalar>(s.levels.front().unknowns);
  s.setup_seconds = timer.elapsed_seconds();
  s.reused_aggregates = reuse;
  log::debug("multigrid: ", s.levels.size(), " levels, operator complexity ",
             s.operator_complexity, ", coarsest ", s.levels.back().unknowns,
             " unknowns, setup ", s.setup_seconds, " s", reuse ? " (aggregates reused)" : "",
             " [checks ", s.check_seconds, ", spectral ", s.spectral_seconds,
             ", aggregation ", s.aggregation_seconds, ", tentative ", s.tentative_seconds,
             ", smoothing ", s.smoothing_seconds, ", Galerkin ", s.galerkin_seconds,
             ", coarse ", s.coarse_seconds, "]");
}

void AmgPreconditioner::apply(const Vector& r, Vector& z) const {
  const Impl& im = *impl_;
  if (!im.ready) throw SolverError("multigrid: apply() called before setup()");
  const Impl::Level& fine = im.levels.front();
  if (r.size() != fine.a.rows) throw SolverError("multigrid: residual has the wrong length");
  fine.b = r;
  im.cycle(0);
  z = fine.x;
}

const AmgStats& AmgPreconditioner::stats() const { return impl_->stats; }
const AmgOptions& AmgPreconditioner::options() const { return impl_->options; }
bool AmgPreconditioner::ready() const { return impl_->ready; }
int AmgPreconditioner::num_levels() const { return static_cast<int>(impl_->levels.size()); }

SparseMatrix AmgPreconditioner::prolongator(int level) const {
  return to_eigen(view_of(impl_->levels.at(static_cast<std::size_t>(level)).p));
}
SparseMatrix AmgPreconditioner::tentative_prolongator(int level) const {
  return to_eigen(view_of(impl_->levels.at(static_cast<std::size_t>(level)).p_tent));
}
SparseMatrix AmgPreconditioner::level_matrix(int level) const {
  return to_eigen(impl_->levels.at(static_cast<std::size_t>(level)).a);
}
Matrix AmgPreconditioner::near_null_space(int level) const {
  return impl_->levels.at(static_cast<std::size_t>(level)).nullspace;
}

Scalar deterministic_dot(const Vector& a, const Vector& b) {
  if (a.size() != b.size()) throw SolverError("dot product of vectors of different lengths");
  return chunked_dot(a.data(), b.data(), static_cast<Index>(a.size()));
}

void symmetric_spmv(const SparseMatrix& a, const Vector& x, Vector& y) {
  if (!a.isCompressed()) throw SolverError("symmetric_spmv needs compressed storage");
  y.resize(a.rows());
  spmv(view_of_symmetric(a), x.data(), y.data());
}

PcgResult pcg_solve(const SparseMatrix& a, const Vector& b, Vector& x,
                    const AmgPreconditioner& m, Scalar tolerance, int max_iterations) {
  const Index n = static_cast<Index>(a.rows());
  if (b.size() != n) throw SolverError("pcg: right-hand side has the wrong length");
  if (x.size() != n) x = Vector::Zero(n);
  const CsrView av = view_of_symmetric(a);
  PcgResult result;
  const Scalar bnorm = std::sqrt(chunked_dot(b.data(), b.data(), n));
  if (bnorm == 0.0) {
    x.setZero();
    result.converged = true;
    return result;
  }
  const Scalar target = tolerance * bnorm;
  Vector r(n), z(n), p(n), q(n);
  int total = 0;
  bool recursive_converged = false;
  Scalar best_true = std::numeric_limits<Scalar>::max();
  // Up to two restarts from the true residual, in case the recursively
  // updated residual has drifted from b - A x.
  for (int restart = 0; restart < 3; ++restart) {
    residual(av, x.data(), b.data(), r.data());
    Scalar rnorm = std::sqrt(chunked_dot(r.data(), r.data(), n));
    result.relative_residual = rnorm / bnorm;
    if (rnorm <= target) {
      result.converged = true;
      break;
    }
    m.apply(r, z);
    p = z;
    Scalar rz = chunked_dot(r.data(), z.data(), n);
    bool converged = false;
    while (total < max_iterations) {
      spmv(av, p.data(), q.data());
      const Scalar pq = chunked_dot(p.data(), q.data(), n);
      if (!(pq > 0.0) || !std::isfinite(pq)) {
        std::ostringstream os;
        os << "preconditioned CG broke down at iteration " << total + 1 << " (p^T A p = " << pq
           << "): the stiffness matrix or the multigrid preconditioner is not positive "
              "definite. Check the supports, or use a direct solver to confirm";
        throw SolverError(os.str());
      }
      const Scalar alpha = rz / pq;
#if defined(SPARLAB_HAVE_OPENMP)
#pragma omp parallel for schedule(static) if (n >= kParallelRows)
#endif
      for (Index i = 0; i < n; ++i) {
        x(i) += alpha * p(i);
        r(i) -= alpha * q(i);
      }
      ++total;
      rnorm = std::sqrt(chunked_dot(r.data(), r.data(), n));
      if (rnorm <= target) {
        converged = true;
        break;
      }
      m.apply(r, z);
      const Scalar rz_new = chunked_dot(r.data(), z.data(), n);
      const Scalar beta = rz_new / rz;
      rz = rz_new;
#if defined(SPARLAB_HAVE_OPENMP)
#pragma omp parallel for schedule(static) if (n >= kParallelRows)
#endif
      for (Index i = 0; i < n; ++i) p(i) = z(i) + beta * p(i);
    }
    residual(av, x.data(), b.data(), r.data());
    result.relative_residual = std::sqrt(chunked_dot(r.data(), r.data(), n)) / bnorm;
    if (result.relative_residual <= tolerance) {
      result.converged = true;
      break;
    }
    if (!converged) break;  // iteration budget spent
    recursive_converged = true;
    // Restarting did not improve the true residual: it has reached the
    // accuracy with which b - A x can be evaluated in floating point, which
    // for a stiff system can sit above a very tight tolerance.
    if (!(result.relative_residual < 0.5 * best_true)) break;
    best_true = result.relative_residual;
  }
  // Attainable accuracy: the recursion met the tolerance and restarts no
  // longer reduce the true residual. That residual is returned for the caller
  // to check (StaticAnalysis enforces its residual_tolerance on it).
  if (!result.converged && recursive_converged) result.converged = true;
  result.iterations = total;
  return result;
}

}  // namespace sparlab
