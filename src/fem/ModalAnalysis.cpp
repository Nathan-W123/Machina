#include "sparlab/fem/ModalAnalysis.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"
#include "sparlab/fem/LinearSolver.hpp"
#include "sparlab/fem/ModelDiagnostics.hpp"

#include <Eigen/Dense>
#include <Eigen/SparseCholesky>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <sstream>

namespace sparlab {
namespace {

constexpr Scalar kTwoPi = 6.283185307179586476925286766559;

/// Deterministic starting subspace following Bathe: the first vector is the
/// mass diagonal, the following ones are unit vectors at the DOFs with the
/// largest m_ii / k_ii ratio (the most "flexible and heavy" DOFs), and any
/// remainder is filled from a seeded PRNG.
Matrix starting_subspace(const SparseMatrix& k, const SparseMatrix& m, int q,
                         unsigned int seed) {
  const Eigen::Index n = k.rows();
  Matrix x = Matrix::Zero(n, q);

  const Vector kd = k.diagonal();
  const Vector md = m.diagonal();
  x.col(0) = md;

  std::vector<Eigen::Index> order(static_cast<std::size_t>(n));
  std::iota(order.begin(), order.end(), Eigen::Index{0});
  std::stable_sort(order.begin(), order.end(), [&](Eigen::Index a, Eigen::Index b) {
    const Scalar ra = kd(a) > 0.0 ? md(a) / kd(a) : 0.0;
    const Scalar rb = kd(b) > 0.0 ? md(b) / kd(b) : 0.0;
    return ra > rb;
  });

  int col = 1;
  for (std::size_t i = 0; i < order.size() && col < q; ++i, ++col) {
    x(order[i], col) = 1.0;
  }
  if (col < q) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<Scalar> dist(-1.0, 1.0);
    for (; col < q; ++col) {
      for (Eigen::Index i = 0; i < n; ++i) x(i, col) = dist(rng);
    }
  }
  return x;
}

}  // namespace

Scalar cantilever_bending_frequency(Scalar e, Scalar rho, Scalar length, Scalar height,
                                    Scalar thickness, int mode_number) {
  static const Scalar beta_l[5] = {1.87510407, 4.69409113, 7.85475744, 10.99554073,
                                   14.13716839};
  if (mode_number < 1 || mode_number > 5) {
    throw ConfigError(
        "cantilever_bending_frequency supports mode numbers 1..5 (tabulated beta_n L)");
  }
  if (!(e > 0.0) || !(rho > 0.0) || !(length > 0.0) || !(height > 0.0) ||
      !(thickness > 0.0)) {
    throw ConfigError("cantilever_bending_frequency requires positive E, rho, L, h, t");
  }
  // I = t h^3 / 12, A = h t  =>  EI / (rho A L^4) = E h^2 / (12 rho L^4)
  const Scalar bl = beta_l[mode_number - 1];
  return (bl * bl) / (kTwoPi * length * length) *
         std::sqrt(e * height * height / (12.0 * rho));
}

Scalar rod_axial_frequency(Scalar e, Scalar rho, Scalar length, int mode_number) {
  if (mode_number < 1) {
    throw ConfigError("rod_axial_frequency requires mode_number >= 1");
  }
  if (!(e > 0.0) || !(rho > 0.0) || !(length > 0.0)) {
    throw ConfigError("rod_axial_frequency requires positive E, rho and L");
  }
  const Scalar wave_speed = std::sqrt(e / rho);
  return static_cast<Scalar>(2 * mode_number - 1) * wave_speed / (4.0 * length);
}

ModalResult solve_modal(const FemModel& model, const Assembler& assembler,
                        const ModalAnalysisOptions& options,
                        const Vector* stiffness_scale, const Vector* mass_scale) {
  if (options.num_modes < 1) {
    throw ConfigError("modal analysis requires num_modes >= 1");
  }

  const SparseMatrix k_full = assembler.assemble_stiffness(stiffness_scale);
  const SparseMatrix m_full = assembler.assemble_mass(options.mass_type, mass_scale);
  const SparseMatrix k = assembler.reduce_free_free(k_full);
  const SparseMatrix m = assembler.reduce_free_free(m_full);

  const Eigen::Index n = k.rows();
  if (n == 0) {
    throw ModelError(
        "modal analysis has no free degrees of freedom; every DOF is constrained");
  }

  ModalResult result;
  // Every translational direction carries the full mass once, so the sum of
  // the assembled matrix is dim times the structural mass. A shell or beam
  // model also carries rotary inertia on its rotational DOFs; its mass is the
  // kinetic energy of a unit rigid translation along x, t^T M t.
  if (model.dofs_per_node() == model.dim()) {
    result.total_mass = m_full.sum() / static_cast<Scalar>(model.dim());
  } else {
    Vector t = Vector::Zero(model.dofs().num_dofs());
    for (Index node = 0; node < model.mesh().num_nodes(); ++node) {
      t(model.dofs().dof(node, 0)) = 1.0;
    }
    result.total_mass = t.dot(m_full * t);
  }

  const int m_req = std::min<int>(options.num_modes, static_cast<int>(n));
  if (m_req < options.num_modes) {
    std::ostringstream os;
    os << "requested " << options.num_modes << " modes but the model has only " << n
       << " free DOFs; returning " << m_req;
    result.warnings.push_back(os.str());
    log::warn(os.str());
  }

  const int q = std::min<int>(static_cast<int>(n), std::max(2 * m_req, m_req + 8));
  result.subspace_size = q;

  Vector lambda;
  Matrix phi;

  if (q >= static_cast<int>(n) || n <= 400) {
    // Small system: a dense generalised symmetric eigensolve is both faster and
    // a useful independent reference for the iterative path.
    const Matrix kd(k);
    const Matrix md(m);
    Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> ges(kd, md);
    if (ges.info() != Eigen::Success) {
      throw SolverError(
          "dense generalised eigensolve failed; the mass matrix is not positive "
          "definite (check material.density and the mass interpolation floor)");
    }
    lambda = ges.eigenvalues().head(m_req);
    phi = ges.eigenvectors().leftCols(m_req);
    result.converged = true;
    result.iterations = 0;
    result.final_change = 0.0;
  } else {
    // Shift-invert subspace iteration.
    Scalar sigma = 0.0;
    if (options.shift_factor != 0.0) {
      const Vector kd = k.diagonal();
      const Vector md = m.diagonal();
      Scalar min_ratio = std::numeric_limits<Scalar>::max();
      for (Eigen::Index i = 0; i < n; ++i) {
        if (md(i) > 0.0) min_ratio = std::min(min_ratio, kd(i) / md(i));
      }
      if (min_ratio < std::numeric_limits<Scalar>::max()) {
        sigma = options.shift_factor * min_ratio;
      }
    }

    SparseMatrix a = k;
    if (sigma != 0.0) a = k - sigma * m;

    LinearSolverOptions lin = options.linear;
    const bool iterative_type = lin.type == LinearSolverType::AmgCg ||
                                lin.type == LinearSolverType::ConjugateGradient ||
                                lin.type == LinearSolverType::Auto;
    if (sigma > 0.0 && iterative_type) {
      // A positive shift can make K - sigma M indefinite, which CG cannot
      // handle; the factorisation handles it and its pivots report trouble.
      log::info("modal analysis: shift sigma = ", sigma,
                " > 0, so the shifted matrix may be indefinite; using a sparse "
                "Cholesky factorisation instead of ", to_string(lin.type));
      lin.type = LinearSolverType::SimplicialLdlt;
    }
    auto solver = make_linear_solver(lin);
    DofLayout layout;
    layout.dim = model.dim();
    layout.dofs_per_node = model.dofs_per_node();
    layout.coordinates = &model.mesh().coordinates();
    layout.unknowns = &model.dofs().free_dofs();
    solver->set_layout(layout);
    solver->factorize(a);
    result.linear_solver = solver->name();
    const bool warm = solver->iterative() && lin.warm_start;

    Matrix x = starting_subspace(k, m, q, options.seed);
    Vector prev = Vector::Constant(m_req, std::numeric_limits<Scalar>::max());
    Vector mu_previous;

    for (int iter = 1; iter <= options.max_iterations; ++iter) {
      const Matrix y = m * x;
      Matrix xbar(n, q);
      for (int c = 0; c < q; ++c) {
        // Near convergence (K - sigma M)^-1 M x = x / (lambda - sigma): the
        // Ritz pair of the last iteration is an initial guess close to the
        // answer.
        if (warm && mu_previous.size() == q && mu_previous(c) - sigma > 0.0) {
          const Vector guess = x.col(c) / (mu_previous(c) - sigma);
          xbar.col(c) = solver->solve_from(y.col(c), guess);
        } else {
          xbar.col(c) = solver->solve(y.col(c));
        }
        result.linear_iterations += solver->last_iterations();
      }
      if (!xbar.allFinite()) {
        throw SolverError(
            "subspace iteration produced a non-finite basis; the shifted matrix "
            "K - sigma M is singular. Reduce modal.shift_factor");
      }

      const Matrix kr = (xbar.transpose() * (k * xbar)).eval();
      const Matrix mr = (xbar.transpose() * (m * xbar)).eval();
      const Matrix krs = 0.5 * (kr + kr.transpose());
      const Matrix mrs = 0.5 * (mr + mr.transpose());

      Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> ges(krs, mrs);
      if (ges.info() != Eigen::Success) {
        std::ostringstream os;
        os << "the projected " << q << " x " << q
           << " eigenproblem became numerically singular at subspace iteration " << iter
           << "; the trial vectors lost independence. Reduce modal.num_modes or "
              "tighten modal.tolerance";
        throw SolverError(os.str());
      }

      const Vector mu = ges.eigenvalues();
      x = xbar * ges.eigenvectors();
      mu_previous = mu;

      const Vector current = mu.head(m_req);
      Scalar change = 0.0;
      for (int i = 0; i < m_req; ++i) {
        const Scalar denom = std::max(std::abs(current(i)), 1.0e-300);
        change = std::max(change, std::abs(current(i) - prev(i)) / denom);
      }
      prev = current;
      result.iterations = iter;
      result.final_change = change;

      // The highest requested mode converges last, so the eigenvalue change
      // alone can look settled while its eigenvector is still inaccurate. The
      // eigenpair residual is therefore part of the stopping rule rather than
      // something checked only afterwards. Two sparse matrix-vector products
      // per mode is negligible next to the q back-substitutions above.
      Scalar max_residual = 0.0;
      for (int i = 0; i < m_req; ++i) {
        const Vector mx = m * x.col(i);
        const Vector residual = k * x.col(i) - current(i) * mx;
        const Scalar reference = std::max((current(i) * mx).norm(), 1.0e-300);
        max_residual = std::max(max_residual, residual.norm() / reference);
      }
      log::trace("subspace iteration ", iter, ": relative eigenvalue change ", change,
                 ", max eigenpair residual ", max_residual);

      if (change < options.tolerance && max_residual < options.residual_tolerance) {
        result.converged = true;
        break;
      }
    }

    if (!result.converged) {
      std::ostringstream os;
      os << "subspace iteration did not converge in " << options.max_iterations
         << " iterations; the relative eigenvalue change stalled at "
         << result.final_change << " against a tolerance of " << options.tolerance
         << " (the stopping rule also requires every eigenpair residual below "
         << options.residual_tolerance
         << "). Increase modal.max_iterations, loosen modal.tolerance or "
            "modal.residual_tolerance, or request fewer modes";
      throw ConvergenceError(os.str());
    }

    lambda = prev;
    phi = x.leftCols(m_req);
  }

  // Re-normalise so that phi^T M phi = I exactly (the projection guarantees it
  // only up to round-off).
  for (int i = 0; i < m_req; ++i) {
    const Scalar mnorm = std::sqrt(std::max(phi.col(i).dot(m * phi.col(i)), 0.0));
    if (!(mnorm > 0.0)) {
      throw SolverError(
          "a computed mode shape has zero modal mass; the mass matrix is singular on "
          "the free DOFs");
    }
    phi.col(i) /= mnorm;
  }

  // Validity screening.
  const Scalar stiffness_scale_norm =
      k.coeffs().size() > 0 ? k.coeffs().cwiseAbs().maxCoeff() : 0.0;
  const Scalar mass_scale_norm =
      m.coeffs().size() > 0 ? m.coeffs().cwiseAbs().maxCoeff() : 1.0;
  const Scalar lambda_scale =
      mass_scale_norm > 0.0 ? stiffness_scale_norm / mass_scale_norm : 1.0;

  result.eigenvalues = lambda;
  result.angular_frequencies.setZero(m_req);
  result.frequencies_hz.setZero(m_req);
  result.modal_residuals.setZero(m_req);

  for (int i = 0; i < m_req; ++i) {
    const Scalar lam = lambda(i);
    if (!std::isfinite(lam)) {
      std::ostringstream os;
      os << "eigenvalue " << i << " is not finite; the eigenproblem is ill-posed";
      throw SolverError(os.str());
    }
    if (lam < -options.rigid_body_ratio * lambda_scale) {
      std::ostringstream os;
      os << "eigenvalue " << i << " is negative (" << lam
         << " 1/s^2). K_ff is not positive definite, which is physically impossible "
            "for a properly constrained linear-elastic model: check the boundary "
            "conditions and the SIMP stiffness floor";
      throw SolverError(os.str());
    }
    if (lam < options.rigid_body_ratio * lambda_scale) {
      std::ostringstream os;
      os << "eigenvalue " << i << " is " << lam
         << " 1/s^2, indistinguishable from zero at this stiffness/mass scale ("
         << lambda_scale
         << "): the model retains a rigid-body or mechanism mode. Add constraints.";
      result.warnings.push_back(os.str());
      log::warn(os.str());
    }

    const Scalar omega = std::sqrt(std::max(lam, 0.0));
    result.angular_frequencies(i) = omega;
    result.frequencies_hz(i) = omega / kTwoPi;

    const Vector residual = k * phi.col(i) - lam * (m * phi.col(i));
    const Scalar reference = std::max((lam * (m * phi.col(i))).norm(), 1.0e-300);
    result.modal_residuals(i) = residual.norm() / reference;
    if (result.modal_residuals(i) > options.residual_tolerance) {
      std::ostringstream os;
      os << "mode " << i << " has an eigenpair residual of " << result.modal_residuals(i)
         << ", above the tolerance " << options.residual_tolerance;
      result.warnings.push_back(os.str());
      log::warn(os.str());
    }
  }

  // Expand mode shapes to full length (prescribed DOFs stay zero: homogeneous
  // Dirichlet data is the correct boundary condition for free vibration).
  result.mode_shapes.setZero(model.dofs().num_dofs(), m_req);
  const auto& free = model.dofs().free_dofs();
  for (int i = 0; i < m_req; ++i) {
    for (std::size_t r = 0; r < free.size(); ++r) {
      result.mode_shapes(free[r], i) = phi(static_cast<Eigen::Index>(r), i);
    }
  }

  // Fraction of modal kinetic energy carried by low-density elements. Only
  // meaningful for a SIMP design, where it exposes spurious localised modes.
  result.modal_mass_fraction.setZero(m_req);
  if (mass_scale != nullptr) {
    const Mesh& mesh = model.mesh();
    const int npe = mesh.nodes_per_elem();
    const Scalar threshold = 0.3;
    Vector pe;
    for (int i = 0; i < m_req; ++i) {
      Scalar total = 0.0;
      Scalar low = 0.0;
      const Vector shape = result.mode_shapes.col(i);
      for (Index e = 0; e < mesh.num_elements(); ++e) {
        model.dofs().gather(mesh.element_nodes(e), npe, shape, pe);
        const Scalar s = (*mass_scale)(e);
        const Scalar ke = s * pe.dot(assembler.element_mass(e) * pe);
        total += ke;
        if (s < threshold) low += ke;
      }
      result.modal_mass_fraction(i) = total > 0.0 ? low / total : 0.0;
    }
  }

  log::info("modal analysis: ", m_req, " modes, subspace ", q, ", ", result.iterations,
            " iteration(s), f1 = ", result.frequencies_hz(0), " Hz, total mass ",
            result.total_mass, " kg");
  return result;
}

}  // namespace sparlab
