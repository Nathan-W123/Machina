#include "sparlab/fem/Buckling.hpp"

#include "sparlab/fem/Loads.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"
#include "sparlab/fem/StaticAnalysis.hpp"

#include <Eigen/Dense>
#include <Eigen/SparseCholesky>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <random>
#include <sstream>

namespace sparlab {
namespace {

constexpr Scalar kPi = 3.14159265358979323846264338327950288;

/// Eigenvalues of the pencil below this fraction of the largest |mu| are
/// treated as zero: a load factor 1e10 times the smallest one means nothing.
constexpr Scalar kZeroFraction = 1.0e-10;

/// Systems up to this many free DOFs use the dense generalised eigensolver.
constexpr Eigen::Index kDenseLimit = 400;

void validate(const BucklingOptions& options) {
  if (options.num_modes < 1) throw ConfigError("buckling.num_modes must be at least 1");
  if (!(options.tolerance > 0.0)) throw ConfigError("buckling.tolerance must be positive");
  if (!(options.residual_tolerance > 0.0)) {
    throw ConfigError("buckling.residual_tolerance must be positive");
  }
  if (options.max_iterations < 1) {
    throw ConfigError("buckling.max_iterations must be at least 1");
  }
}

/// Orthonormal basis (Euclidean) of the columns of `z` by column-pivoted QR,
/// dropping numerically dependent columns.
Matrix orthonormal_basis(const Matrix& z) {
  Eigen::ColPivHouseholderQR<Matrix> qr(z);
  qr.setThreshold(1.0e-12);
  const Eigen::Index rank = qr.rank();
  Matrix q = qr.householderQ() * Matrix::Identity(z.rows(), rank);
  return q;
}

}  // namespace

SparseMatrix assemble_geometric_stiffness(const FemModel& model, const Assembler& assembler,
                                          const Vector& displacement,
                                          const Vector* stress_scale,
                                          const Vector* temperature) {
  const Mesh& mesh = model.mesh();
  if (displacement.size() != model.dofs().num_dofs()) {
    std::ostringstream os;
    os << "geometric stiffness: displacement vector has length " << displacement.size()
       << " but the model has " << model.dofs().num_dofs() << " DOFs";
    throw ModelError(os.str());
  }
  if (stress_scale != nullptr && stress_scale->size() != mesh.num_elements()) {
    std::ostringstream os;
    os << "geometric stiffness: stress scale vector has length " << stress_scale->size()
       << " but the mesh has " << mesh.num_elements() << " elements";
    throw ModelError(os.str());
  }
  const int npe = mesh.nodes_per_elem();
  const int edofs = npe * model.dofs_per_node();
  const Element& element = model.element();
  const std::vector<IntegrationPoint> rule = element.integration_rule(model.integration());
  return assembler.assemble_elementwise([&](Index e) -> Matrix {
    const Scalar s = stress_scale != nullptr ? (*stress_scale)(e) : 1.0;
    if (s == 0.0) return Matrix::Zero(edofs, edofs);
    const Vector ue = model.dofs().gather(mesh.element_nodes(e), npe, displacement);
    const Matrix coords = mesh.element_coordinates(e);
    if (temperature == nullptr || model.material_of(e).thermal_expansion() == 0.0) {
      return element.geometric_stiffness(coords, model.constitutive_of(e), ue, s,
                                         model.thickness(), model.integration());
    }
    // Thermal prestress: sigma = s D (B u - eps0) at each point (with the
    // condensed internal modes of an incompatible-mode element).
    const Vector alpha = linear_internal_parameters(model, e, ue, temperature);
    std::vector<Vector> stresses;
    stresses.reserve(rule.size());
    for (const IntegrationPoint& ip : rule) {
      const Vector eps0 = element_thermal_strain(model, e, ip.point, *temperature);
      stresses.push_back(s * (model.constitutive_of(e) *
                              (linear_point_strain(model, e, ip.point, ue, alpha) - eps0)));
    }
    return element.geometric_stiffness_of_stress(coords, stresses, model.thickness(),
                                                 model.integration());
  });
}

BucklingResult solve_buckling(const FemModel& model, const Assembler& assembler,
                              const SparseMatrix& k_full, const SparseMatrix& k_g_full,
                              const BucklingOptions& options, const FreeSolve& solve,
                              const Matrix* initial_subspace,
                              const Vector* stiffness_scale,
                              const std::vector<char>* solid_mask) {
  validate(options);
  const SparseMatrix k = assembler.reduce_free_free(k_full);
  const SparseMatrix g = -assembler.reduce_free_free(k_g_full);  // G = -K_G
  const Eigen::Index n = k.rows();
  if (n == 0) {
    throw ModelError("buckling analysis has no free degrees of freedom; every DOF is "
                     "constrained");
  }
  if (g.rows() != n) throw ModelError("buckling: K and K_G have different sizes");
  const Scalar g_norm = g.nonZeros() > 0 ? g.coeffs().cwiseAbs().maxCoeff() : 0.0;

  BucklingResult result;
  const int m_req = static_cast<int>(std::min<Eigen::Index>(options.num_modes, n));
  const int q = static_cast<int>(std::min<Eigen::Index>(n, std::max(2 * m_req, m_req + 8)));
  result.subspace_size = q;

  // Descending eigenvalues mu = 1 / lambda of the pencil (G, K) and their
  // K-normalised vectors, for the top of the spectrum.
  Vector mu_top;
  Matrix phi_top;

  if (!(g_norm > 0.0)) {
    // No stress at all (an unloaded case): nothing can buckle.
    result.converged = true;
    result.no_positive_load_factor = true;
  } else if (n <= kDenseLimit || q >= n) {
    const Matrix kd(k);
    const Matrix gd(g);
    Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> ges(gd, kd);
    if (ges.info() != Eigen::Success) {
      throw SolverError("the dense buckling eigensolve failed; K_ff is not positive "
                        "definite. Check the boundary conditions");
    }
    const Vector mu = ges.eigenvalues();  // ascending
    mu_top = mu.reverse().head(m_req);
    phi_top = ges.eigenvectors().rowwise().reverse().leftCols(m_req);
    result.subspace = phi_top;
    result.converged = true;
  } else {
    std::unique_ptr<LinearSolver> own;
    FreeSolve kinv = solve;
    if (!kinv) {
      own = make_linear_solver(options.linear);
      DofLayout layout;
      layout.dim = model.dim();
      layout.dofs_per_node = model.dofs_per_node();
      layout.coordinates = &model.mesh().coordinates();
      layout.unknowns = &model.dofs().free_dofs();
      own->set_layout(layout);
      own->factorize(k);
      result.linear_solver = own->name();
      kinv = [&own, &result](const Vector& b) {
        Vector x = own->solve(b);
        result.linear_iterations += own->last_iterations();
        return x;
      };
    }

    std::mt19937 rng(options.seed);
    std::uniform_real_distribution<Scalar> dist(-1.0, 1.0);
    const auto random_column = [&]() {
      Vector v(n);
      for (Eigen::Index i = 0; i < n; ++i) v(i) = dist(rng);
      return v;
    };
    Matrix x(n, q);
    int col = 0;
    if (initial_subspace != nullptr && initial_subspace->rows() == n) {
      // Keep two random directions so modes the previous design did not have
      // can still enter.
      const int keep = static_cast<int>(std::min<Eigen::Index>(initial_subspace->cols(), q - 2));
      x.leftCols(keep) = initial_subspace->leftCols(keep);
      col = keep;
    }
    for (; col < q; ++col) x.col(col) = random_column();

    // The buckling spectral transformation: iterate with (K - sigma G)^{-1} K,
    // whose eigenvalues nu = lambda / (lambda - sigma) exceed 1 for the load
    // factors above sigma and lie in (0, 1) for every negative one. For
    // 0 < sigma < lambda_1 the shifted matrix K - sigma G = K + sigma K_G is
    // positive definite, and by Sylvester's law of inertia the number of
    // negative pivots of its LDL^T factorisation is the number of load
    // factors in (0, sigma): that count places sigma safely.
    Eigen::SimplicialLDLT<SparseMatrix> shifted;
    const auto inertia = [&](Scalar s) -> int {
      shifted.compute(SparseMatrix(k - s * g));
      ++result.factorizations;
      if (shifted.info() != Eigen::Success) return -1;
      const Vector pivots = shifted.vectorD();
      int negative = 0;
      for (Eigen::Index i = 0; i < pivots.size(); ++i) {
        if (!(pivots(i) != 0.0) || !std::isfinite(pivots(i))) return -1;
        negative += pivots(i) < 0.0 ? 1 : 0;
      }
      return negative;
    };
    // Choose sigma below the smallest positive load factor. A guess that is
    // an estimate of lambda_1 is lowered until no load factor lies below it;
    // without one the search grows from 1 (the load case itself) by factors
    // of 8. Returns false when no positive load factor lies below 1e15.
    bool transformed = false;
    const auto choose_sigma = [&](Scalar guess) -> bool {
      Scalar s = guess > 0.0 ? guess : 1.0;
      int negative = inertia(s);
      if (negative != 0) {
        for (int halvings = 1; negative != 0; ++halvings) {
          if (halvings > 120) {
            throw SolverError("the buckling shift search found no shift below the smallest "
                              "load factor; K_ff is not positive definite");
          }
          s *= 0.5;
          negative = inertia(s);
        }
      } else if (!(guess > 0.0)) {
        Scalar clean = s;
        for (;;) {
          if (s > 1.0e15) return false;
          const Scalar next = 8.0 * s;
          if (inertia(next) != 0) break;
          s = next;
          clean = s;
        }
        s = clean;
        if (inertia(s) != 0) throw SolverError("inconsistent inertia in the buckling shift search");
      }
      result.sigma = s;
      transformed = true;
      return true;
    };
    bool refined = false;
    if (options.transform == BucklingTransform::ShiftInvert) {
      if (!choose_sigma(options.shift_hint > 0.0 ? 0.8 * options.shift_hint : 0.0)) {
        result.no_positive_load_factor = true;
        result.converged = true;
      }
    }

    Vector previous = Vector::Constant(m_req, std::numeric_limits<Scalar>::quiet_NaN());
    Vector mu_desc;
    Matrix ritz;
    for (int iter = 1; iter <= options.max_iterations && !result.no_positive_load_factor;
         ++iter) {
      Matrix z(n, q);
      if (transformed) {
        const Matrix y = k * x;
        for (int c = 0; c < q; ++c) z.col(c) = shifted.solve(y.col(c));
      } else {
        const Matrix y = g * x;
        for (int c = 0; c < q; ++c) z.col(c) = kinv(y.col(c));
      }
      if (!z.allFinite()) {
        throw SolverError("buckling subspace iteration produced a non-finite basis; the "
                          "stiffness matrix is singular or severely ill-conditioned");
      }
      const Matrix basis = orthonormal_basis(z);
      const Eigen::Index r = basis.cols();
      if (r < m_req) {
        std::ostringstream os;
        os << "the buckling subspace collapsed to " << r << " independent vectors, fewer "
           << "than the " << m_req << " requested modes: the stress state loads fewer "
           << "independent deformation patterns than that. Request fewer modes";
        throw SolverError(os.str());
      }
      const Matrix gr = basis.transpose() * (g * basis);
      const Matrix kr = basis.transpose() * (k * basis);
      Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> ges(0.5 * (gr + gr.transpose()),
                                                           0.5 * (kr + kr.transpose()));
      if (ges.info() != Eigen::Success) {
        std::ostringstream os;
        os << "the projected " << r << " x " << r << " buckling eigenproblem failed at "
           << "subspace iteration " << iter << "; K_ff is not positive definite on the "
           << "subspace. Check the boundary conditions";
        throw SolverError(os.str());
      }
      mu_desc = ges.eigenvalues().reverse();
      ritz = basis * ges.eigenvectors().rowwise().reverse();
      x.leftCols(r) = ritz;
      for (Eigen::Index c = r; c < q; ++c) x.col(c) = random_column();

      const Scalar scale = mu_desc.cwiseAbs().maxCoeff();
      const Scalar zero = kZeroFraction * scale;
      const Scalar mu_m = mu_desc(m_req - 1) > zero ? mu_desc(m_req - 1) : 0.0;

      // Convergence of the m top Ritz values: relative change for the
      // positive ones (it is the relative change of lambda), absolute
      // against the spectrum's scale for the others, which are not reported.
      Scalar change = 0.0;
      Scalar residual = 0.0;
      for (int i = 0; i < m_req; ++i) {
        const Scalar mu = mu_desc(i);
        const Scalar denom = mu > zero ? mu : scale;
        const Scalar delta = std::isnan(previous(i)) ? std::numeric_limits<Scalar>::infinity()
                                                     : std::abs(mu - previous(i)) / denom;
        change = std::max(change, delta);
        if (mu > zero) {
          const Vector kphi = k * ritz.col(i);
          const Vector res = g * ritz.col(i) - mu * kphi;
          residual = std::max(residual, res.norm() / std::max(mu * kphi.norm(), 1.0e-300));
        }
      }
      previous = mu_desc.head(m_req);
      result.iterations = iter;
      result.final_change = change;
      log::trace("buckling subspace iteration ", iter, ": change ", change, ", residual ",
                 residual, transformed ? ", transformed, sigma " : ", plain",
                 transformed ? result.sigma : 0.0);
      if (change < options.tolerance && residual < options.residual_tolerance) {
        result.converged = true;
        break;
      }

      if (!transformed && options.transform == BucklingTransform::Auto) {
        // Negative eigenvalues larger in magnitude than the m-th wanted one
        // take the slots the wanted ones converge in; once more of them than
        // the spare slots appear, switch to the transformation, which maps
        // all of them below the wanted ones.
        int competing = 0;
        for (Eigen::Index j = 0; j < r; ++j) competing += mu_desc(j) < -mu_m ? 1 : 0;
        if (competing + m_req > q - 2) {
          const Scalar estimate = mu_desc(0) > zero ? 0.8 / mu_desc(0) : 0.0;
          if (!choose_sigma(options.shift_hint > 0.0 ? 0.8 * options.shift_hint : estimate)) {
            result.no_positive_load_factor = true;
            result.converged = true;
            break;
          }
          previous.setConstant(std::numeric_limits<Scalar>::quiet_NaN());
        }
      } else if (transformed && !refined && mu_desc(0) > zero && change < 1.0e-3 &&
                 result.sigma < 0.5 / mu_desc(0)) {
        // A shift found by the coarse search can sit far below lambda_1,
        // where the transformation separates the modes poorly; move it up to
        // 0.8 of the estimate once lambda_1 has settled.
        refined = true;
        const Scalar before = result.sigma;
        if (!choose_sigma(0.8 / mu_desc(0))) {
          throw SolverError("the refined buckling shift found no load factor");
        }
        log::debug("buckling: shift raised from ", before, " to ", result.sigma);
      }
    }
    result.transformed = transformed;
    if (!result.converged) {
      std::ostringstream os;
      os << "buckling subspace iteration did not converge in " << options.max_iterations
         << " iterations: the relative change of the requested load factors stalled at "
         << result.final_change << " against a tolerance of " << options.tolerance
         << " (every eigenpair residual must also fall below "
         << options.residual_tolerance
         << "). Increase buckling.max_iterations, request fewer modes, or loosen "
            "buckling.tolerance";
      throw ConvergenceError(os.str());
    }
    if (!result.no_positive_load_factor) {
      mu_top = mu_desc.head(m_req);
      phi_top = ritz.leftCols(m_req);
      result.subspace = ritz;
    }
  }

  // Keep the positive eigenvalues: the smallest positive load factors.
  std::vector<int> keep;
  if (mu_top.size() > 0) {
    const Scalar zero = kZeroFraction * mu_top.cwiseAbs().maxCoeff();
    for (Eigen::Index i = 0; i < mu_top.size(); ++i) {
      if (mu_top(i) > zero) keep.push_back(static_cast<int>(i));
    }
  }
  const int found = static_cast<int>(keep.size());
  result.load_factors.resize(found);
  result.residuals.resize(found);
  result.solid_energy_fraction.resize(found);
  result.mode_shapes.setZero(model.dofs().num_dofs(), found);
  const auto& free = model.dofs().free_dofs();
  for (int j = 0; j < found; ++j) {
    const int i = keep[static_cast<std::size_t>(j)];
    Vector phi = phi_top.col(i);
    const Scalar knorm = std::sqrt(std::max(phi.dot(k * phi), 0.0));
    if (!(knorm > 0.0)) throw SolverError("a buckling mode has zero strain energy");
    phi /= knorm;
    const Scalar lambda = 1.0 / mu_top(i);
    result.load_factors(j) = lambda;
    const Vector kphi = k * phi;
    result.residuals(j) = (kphi - lambda * (g * phi)).norm() / std::max(kphi.norm(), 1.0e-300);
    for (std::size_t r = 0; r < free.size(); ++r) {
      result.mode_shapes(free[r], j) = phi(static_cast<Eigen::Index>(r));
    }
    if (result.residuals(j) > options.residual_tolerance) {
      std::ostringstream os;
      os << "buckling mode " << j + 1 << " has an eigenpair residual of "
         << result.residuals(j) << ", above the tolerance " << options.residual_tolerance;
      result.warnings.push_back(os.str());
      log::warn(os.str());
    }
  }
  if (found == 0) {
    result.no_positive_load_factor = true;
  } else if (found < m_req) {
    std::ostringstream os;
    os << "only " << found << " of the " << m_req
       << " requested modes have a positive load factor; the others belong to the "
          "reversed load";
    result.warnings.push_back(os.str());
    log::info(os.str());
  }

  // Share of each mode's strain energy in solid elements.
  const Mesh& mesh = model.mesh();
  const int npe = mesh.nodes_per_elem();
  Vector pe;
  for (int j = 0; j < found; ++j) {
    if (solid_mask == nullptr) {
      result.solid_energy_fraction(j) = 1.0;
      continue;
    }
    Scalar total = 0.0;
    Scalar solid = 0.0;
    const Vector shape = result.mode_shapes.col(j);
    for (Index e = 0; e < mesh.num_elements(); ++e) {
      model.dofs().gather(mesh.element_nodes(e), npe, shape, pe);
      const Scalar s = stiffness_scale != nullptr ? (*stiffness_scale)(e) : 1.0;
      const Scalar energy = s * pe.dot(assembler.element_stiffness(e) * pe);
      total += energy;
      if ((*solid_mask)[static_cast<std::size_t>(e)]) solid += energy;
    }
    result.solid_energy_fraction(j) = total > 0.0 ? solid / total : 0.0;
  }

  if (result.no_positive_load_factor) {
    log::info("buckling: no positive load factor; the stress state stiffens the structure "
              "against every computed mode (the reversed load may still buckle it)");
  } else {
    log::info("buckling: ", found, " mode(s), lambda_1 = ", result.load_factors(0), ", ",
              result.iterations, " subspace iteration(s), subspace ", q,
              result.transformed ? ", spectral transformation at sigma = " : "",
              result.transformed ? std::to_string(result.sigma) : std::string());
  }
  return result;
}

BucklingResult analyse_buckling(const FemModel& model, const Assembler& assembler,
                                std::size_t load_case, const BucklingOptions& options) {
  if (load_case >= model.load_vectors().size()) {
    throw ConfigError("buckling analysis asked for a load case the model does not have");
  }
  StaticAnalysisOptions static_options;
  static_options.linear = options.linear;
  StaticAnalysis analysis(model, assembler, static_options);
  analysis.prepare();
  const Vector u = analysis.solve_load_vector(model.load_vectors()[load_case]);
  const Vector& temperature = model.load_case_data(load_case).temperature;
  const SparseMatrix k_g = assemble_geometric_stiffness(
      model, assembler, u, nullptr, temperature.size() > 0 ? &temperature : nullptr);
  const DofManager& dofs = model.dofs();
  const FreeSolve solve = [&](const Vector& b) {
    return dofs.restrict_to_free(analysis.solve_homogeneous(dofs.expand(b)));
  };
  BucklingResult result =
      solve_buckling(model, assembler, analysis.stiffness(), k_g, options, solve);
  result.load_case = model.load_case_specs()[load_case].name;
  result.linear_solver = analysis.solver().name();
  return result;
}

Scalar euler_cantilever_load(Scalar e, Scalar second_moment, Scalar length) {
  if (!(e > 0.0) || !(second_moment > 0.0) || !(length > 0.0)) {
    throw ConfigError("euler_cantilever_load needs positive E, I and L");
  }
  return kPi * kPi * e * second_moment / (4.0 * length * length);
}

Scalar engesser_cantilever_load(Scalar e, Scalar poisson, Scalar second_moment, Scalar area,
                                Scalar length) {
  const Scalar euler = euler_cantilever_load(e, second_moment, length);
  const Scalar shear = (5.0 / 6.0) * e / (2.0 * (1.0 + poisson)) * area;
  return euler / (1.0 + euler / shear);
}

}  // namespace sparlab
