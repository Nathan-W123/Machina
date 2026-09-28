/// \file Multigrid.hpp
/// \brief Smoothed-aggregation algebraic multigrid (SA-AMG) and the
///        preconditioned conjugate gradient method built on it.
///
/// A sparse Cholesky factorisation of a 3-D stiffness matrix fills in: its
/// cost grows like \f$n^{2}\f$ and its memory like \f$n^{4/3}\f$, which is
/// what capped the solid models of this code base at a few tens of thousands
/// of unknowns. Conjugate gradients needs only matrix-vector products, and
/// with a multigrid preconditioner the iteration count stays nearly constant
/// as the mesh is refined, so the whole solve scales close to linearly.
///
/// ### Hierarchy (Vanek, Mandel and Brezina, 1996)
/// On each level \f$\ell\f$, starting from the reduced stiffness matrix
/// \f$A_0 = K_{ff}\f$:
///
///  1. **Strength of connection** between the unknowns' nodes \f$I, J\f$
///     uses the Frobenius norms of the matrix blocks,
///     \f$ \|A_{IJ}\|_F \ge \theta \sqrt{\|A_{II}\|_F \|A_{JJ}\|_F} \f$.
///     Weak links - between solid and SIMP void material, for instance -
///     are ignored when aggregating.
///  2. **Aggregation** groups each node with its strongly connected
///     neighbours (greedy, in node order, three passes), so an aggregate is a
///     small patch of the mesh.
///  3. **Tentative prolongator** \f$\hat P\f$: the near-null space \f$B\f$ -
///     the rigid-body modes (2 translations and 1 rotation in 2-D, 3 and 3 in
///     3-D), which the stiffness matrix barely resists - is restricted to each
///     aggregate and orthonormalised by a rank-revealing QR,
///     \f$B_a = Q_a R_a\f$. \f$Q_a\f$ is the aggregate's block of \f$\hat P\f$
///     and \f$R_a\f$ its rows of the coarse near-null space, so
///     \f$\hat P B_c = B\f$ exactly: the coarse space represents every
///     rigid-body motion of every aggregate.
///  4. **Prolongator smoothing** \f$P = (I - \omega D^{-1} A)\hat P\f$ with
///     \f$\omega = \tfrac{4}{3} / \lambda_{max}(D^{-1}A)\f$ (Lanczos
///     estimate), which lowers the energy of the coarse basis functions.
///  5. **Galerkin coarse operator** \f$A_{\ell+1} = P^T A_\ell P\f$
///     (symmetrised), which keeps every level symmetric positive definite.
///
/// The recursion stops at `coarse_size` unknowns, where a dense Cholesky
/// factorisation solves exactly. A coarse operator that is not positive
/// definite means a rigid-body mode survived to the coarse level - the model
/// is under-constrained - and is reported as such.
///
/// ### Cycle and Krylov method
/// One V-cycle with a symmetric smoother (Chebyshev polynomial in
/// \f$D^{-1}A\f$, or a forward / backward Gauss-Seidel pair) is a symmetric
/// positive-definite preconditioner, so the outer method is plain
/// preconditioned CG, stopped when \f$\|b - A x\|_2 \le \text{tol}\,
/// \|b\|_2\f$. An initial guess (the previous design iteration's
/// displacement, the previous subspace iterate) can be supplied.
///
/// ### Determinism
/// Every kernel is either row-parallel with rows written independently
/// (matrix-vector products, the sparse matrix-matrix products of the setup,
/// the Chebyshev smoother) or sequential (Gauss-Seidel, aggregation), and
/// every inner product sums fixed-size chunks in a fixed order. Results are
/// therefore bitwise identical for any number of OpenMP threads.
#pragma once

#include "sparlab/core/Types.hpp"

#include <memory>
#include <string>
#include <vector>

namespace sparlab {

enum class AmgSmoother {
  Chebyshev,             ///< polynomial in D^-1 A: row-parallel, symmetric
  SymmetricGaussSeidel   ///< forward sweep before, backward sweep after: sequential
};

std::string to_string(AmgSmoother smoother);
AmgSmoother parse_amg_smoother(const std::string& text);

struct AmgOptions {
  /// Strength-of-connection threshold theta in [0, 1).
  Scalar strength_threshold = 0.02;
  /// Largest number of levels, the finest included.
  int max_levels = 12;
  /// The hierarchy stops once a level has at most this many unknowns; that
  /// level is solved by a dense Cholesky factorisation.
  Index coarse_size = 1500;
  AmgSmoother smoother = AmgSmoother::Chebyshev;
  /// Chebyshev polynomial degree, or the number of Gauss-Seidel sweeps each
  /// way.
  int smoother_degree = 3;
  /// Chebyshev smooths the eigenvalues of D^-1 A in
  /// [lambda_max / ratio, 1.1 lambda_max].
  Scalar chebyshev_ratio = 30.0;
  /// omega * lambda_max(D^-1 A) of the prolongator smoother.
  Scalar prolongator_damping = 4.0 / 3.0;
  /// Lanczos steps for the lambda_max estimates.
  int lanczos_steps = 12;
  /// Smallest accepted pivot ratio of the coarsest factorisation; below it
  /// the model is reported as under-constrained. SIMP void (a 1e-9 stiffness
  /// floor) stays orders of magnitude above it.
  Scalar coarse_pivot_tolerance = 1.0e-13;
  /// Keep the aggregates and tentative prolongators of the first setup when
  /// the next matrix has the same size and sparsity (the design iterations of
  /// a topology optimisation); smoothing, coarse operators and smoothers are
  /// rebuilt from the new values every time.
  bool reuse_aggregates = true;
};

/// Where each unknown of the reduced system lives, from which the rigid-body
/// near-null space is built. The pointers must stay valid during `setup`.
struct DofLayout {
  int dim = 0;                                  ///< spatial dimension (2 or 3)
  /// DOFs per node: `dim` translations for a continuum model (0 means `dim`),
  /// or 6 - three translations, then three rotations - for a shell or beam
  /// model, whose rigid-body modes include the nodal rotations.
  int dofs_per_node = 0;
  const Matrix* coordinates = nullptr;          ///< dim x num_nodes [m]
  const std::vector<Index>* unknowns = nullptr; ///< global DOF (node * ndpn + k) of
                                                ///< each unknown, ascending
  int nodal_dofs() const { return dofs_per_node > 0 ? dofs_per_node : dim; }
};

struct AmgLevelStats {
  Index unknowns = 0;
  Index nonzeros = 0;
  Index aggregates = 0;       ///< nodes of the next level (0 on the coarsest)
  Scalar lambda_max = 0.0;    ///< Lanczos estimate of lambda_max(D^-1 A)
};

struct AmgStats {
  std::vector<AmgLevelStats> levels;
  /// sum_l nnz(A_l) / nnz(A_0): memory and smoothing work relative to A_0.
  Scalar operator_complexity = 0.0;
  /// sum_l n_l / n_0.
  Scalar grid_complexity = 0.0;
  Scalar setup_seconds = 0.0;
  bool reused_aggregates = false;
  int near_null_space_dimension = 0;
  /// Smallest / largest pivot of the coarsest factorisation.
  Scalar coarse_pivot_ratio = 0.0;
  /// Stored numbers beyond the finest matrix: coarse operators, P and P^T,
  /// and the dense coarse factor.
  Index storage_nonzeros = 0;
  /// Setup time by phase [s]: input checks, spectral estimates, strength and
  /// aggregation, tentative prolongators, prolongator smoothing, Galerkin
  /// products and the coarse factorisation.
  Scalar check_seconds = 0.0;
  Scalar spectral_seconds = 0.0;
  Scalar aggregation_seconds = 0.0;
  Scalar tentative_seconds = 0.0;
  Scalar smoothing_seconds = 0.0;
  Scalar galerkin_seconds = 0.0;
  Scalar coarse_seconds = 0.0;
};

/// The multigrid hierarchy of one matrix, applied as a preconditioner.
class AmgPreconditioner {
 public:
  explicit AmgPreconditioner(AmgOptions options = AmgOptions());
  ~AmgPreconditioner();
  AmgPreconditioner(const AmgPreconditioner&) = delete;
  AmgPreconditioner& operator=(const AmgPreconditioner&) = delete;

  /// Build (or, with `reuse_aggregates`, rebuild numerically) the hierarchy of
  /// the symmetric positive-definite matrix `a`. The matrix is referenced, not
  /// copied, and must outlive every `apply` until the next `setup`.
  /// \throws SolverError for a non-symmetric, non-square or non-finite matrix,
  ///         a layout that does not match it, or a coarse operator that is not
  ///         positive definite (an under-constrained model).
  void setup(const SparseMatrix& a, const DofLayout& layout);

  /// z = M^{-1} r: one V-cycle from a zero initial guess.
  void apply(const Vector& r, Vector& z) const;

  const AmgStats& stats() const;
  const AmgOptions& options() const;
  bool ready() const;

  /// Test access: prolongator (smoothed) and coarse operator of level `l`,
  /// and the near-null space of level `l`.
  /// \{
  SparseMatrix prolongator(int level) const;
  SparseMatrix tentative_prolongator(int level) const;
  SparseMatrix level_matrix(int level) const;
  Matrix near_null_space(int level) const;
  int num_levels() const;
  /// \}

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

struct PcgResult {
  int iterations = 0;
  Scalar relative_residual = 0.0;  ///< true residual ||b - A x|| / ||b|| at exit
  bool converged = false;
};

/// Preconditioned conjugate gradients on the matrix of `m`'s last setup.
/// `x` holds the initial guess on entry (zero or a warm start) and the
/// solution on exit.
/// \throws SolverError on a breakdown (p^T A p <= 0: the matrix or the
///         preconditioner is not positive definite).
PcgResult pcg_solve(const SparseMatrix& a, const Vector& b, Vector& x,
                    const AmgPreconditioner& m, Scalar tolerance, int max_iterations);

/// Deterministic parallel kernels, exposed for the tests.
/// \{
Scalar deterministic_dot(const Vector& a, const Vector& b);
/// y = A x for a symmetric matrix stored column-major (read as CSR).
void symmetric_spmv(const SparseMatrix& a, const Vector& x, Vector& y);
/// \}

}  // namespace sparlab
