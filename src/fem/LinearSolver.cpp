#include "sparlab/fem/LinearSolver.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"

#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>

#include <algorithm>
#include <cmath>
#include <sstream>

namespace sparlab {
namespace {

const char* kSingularHint =
    "the reduced stiffness matrix is singular or nearly singular. Typical causes: "
    "(a) too few displacement boundary conditions, leaving rigid-body or mechanism "
    "modes; (b) a disconnected region of elements with no constraint of its own; "
    "(c) a zero or extremely small SIMP stiffness floor (topology.min_density / "
    "simp.emin_ratio)";

template <typename Solver>
void report_info(const Solver& solver, const std::string& who) {
  if (solver.info() == Eigen::Success) return;
  std::ostringstream os;
  os << who << " failed with Eigen status ";
  switch (solver.info()) {
    case Eigen::NumericalIssue: os << "NumericalIssue"; break;
    case Eigen::NoConvergence: os << "NoConvergence"; break;
    case Eigen::InvalidInput: os << "InvalidInput"; break;
    default: os << "Unknown"; break;
  }
  os << "; " << kSingularHint;
  throw SolverError(os.str());
}

void check_square_finite(const SparseMatrix& a, const std::string& who) {
  if (a.rows() != a.cols()) {
    std::ostringstream os;
    os << who << " requires a square matrix (got " << a.rows() << " x " << a.cols()
       << ")";
    throw SolverError(os.str());
  }
  if (a.rows() == 0) {
    throw SolverError(who +
                      " received an empty system: every degree of freedom is "
                      "constrained, so there is nothing to solve");
  }
  for (Eigen::Index k = 0; k < a.outerSize(); ++k) {
    for (SparseMatrix::InnerIterator it(a, k); it; ++it) {
      if (!std::isfinite(it.value())) {
        std::ostringstream os;
        os << who << ": system matrix entry (" << it.row() << ", " << it.col()
           << ") is not finite";
        throw SolverError(os.str());
      }
    }
  }
}

class LdltSolver final : public LinearSolver {
 public:
  explicit LdltSolver(Scalar pivot_tolerance) : pivot_tol_(pivot_tolerance) {}

  void factorize(const SparseMatrix& a) override {
    check_square_finite(a, name());
    solver_.compute(a);
    report_info(solver_, name());
    const Vector d = solver_.vectorD();
    const Scalar max_abs = d.cwiseAbs().maxCoeff();
    if (!(max_abs > 0.0)) {
      throw SolverError(name() + ": all LDL^T pivots vanish; " + kSingularHint);
    }
    const Scalar min_pivot = d.minCoeff();
    if (min_pivot <= 0.0) {
      std::ostringstream os;
      os << name() << ": LDL^T produced a non-positive pivot (" << min_pivot
         << "), so the matrix is not positive definite. " << kSingularHint;
      throw SolverError(os.str());
    }
    const Scalar ratio = min_pivot / max_abs;
    if (ratio < pivot_tol_) {
      std::ostringstream os;
      os << name() << ": smallest/largest LDL^T pivot ratio is " << ratio
         << ", below the tolerance " << pivot_tol_ << ". " << kSingularHint;
      throw SolverError(os.str());
    }
    log::debug(name(), ": factorised ", a.rows(), " unknowns, ", a.nonZeros(),
               " stored entries, pivot ratio ", ratio);
  }

  Vector solve(const Vector& b) override {
    Vector x = solver_.solve(b);
    report_info(solver_, name() + " back-substitution");
    return x;
  }

  std::string name() const override { return "SimplicialLDLT"; }

  Index storage_nonzeros() const override {
    return static_cast<Index>(solver_.matrixL().nestedExpression().nonZeros() +
                              solver_.vectorD().size());
  }

 private:
  Eigen::SimplicialLDLT<SparseMatrix, Eigen::Lower, Eigen::AMDOrdering<StorageIndex>> solver_;
  Scalar pivot_tol_;
};

class LltSolver final : public LinearSolver {
 public:
  void factorize(const SparseMatrix& a) override {
    check_square_finite(a, name());
    solver_.compute(a);
    report_info(solver_, name());
  }
  Vector solve(const Vector& b) override {
    Vector x = solver_.solve(b);
    report_info(solver_, name() + " back-substitution");
    return x;
  }
  std::string name() const override { return "SimplicialLLT"; }

 private:
  Eigen::SimplicialLLT<SparseMatrix, Eigen::Lower, Eigen::AMDOrdering<StorageIndex>> solver_;
};

class SparseLuSolver final : public LinearSolver {
 public:
  void factorize(const SparseMatrix& a) override {
    check_square_finite(a, name());
    solver_.compute(a);
    report_info(solver_, name());
  }
  Vector solve(const Vector& b) override {
    Vector x = solver_.solve(b);
    report_info(solver_, name() + " back-substitution");
    return x;
  }
  std::string name() const override { return "SparseLU"; }

 private:
  Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<StorageIndex>> solver_;
};

class CgSolver final : public LinearSolver {
 public:
  CgSolver(Scalar tolerance, int max_iterations)
      : tolerance_(tolerance), max_iterations_(max_iterations) {}

  void factorize(const SparseMatrix& a) override {
    check_square_finite(a, name());
    solver_.setTolerance(tolerance_);
    if (max_iterations_ > 0) solver_.setMaxIterations(max_iterations_);
    solver_.compute(a);
    report_info(solver_, name() + " setup");
  }

  Vector solve(const Vector& b) override { return check(solver_.solve(b)); }

  Vector solve_from(const Vector& b, const Vector& initial_guess) override {
    if (initial_guess.size() != b.size()) return solve(b);
    return check(solver_.solveWithGuess(b, initial_guess));
  }

  std::string name() const override { return "ConjugateGradient(Jacobi)"; }
  bool iterative() const override { return true; }
  int last_iterations() const override { return iterations_; }
  Scalar last_error() const override { return error_; }

 private:
  Vector check(Vector x) {
    if (solver_.info() != Eigen::Success) {
      std::ostringstream os;
      os << name() << " did not converge: " << solver_.iterations()
         << " iterations reached an error estimate of " << solver_.error()
         << " against a tolerance of " << tolerance_
         << ". Increase max_iterations, loosen the tolerance, or switch to a direct "
            "solver";
      throw ConvergenceError(os.str());
    }
    iterations_ = static_cast<int>(solver_.iterations());
    error_ = solver_.error();
    return x;
  }

  Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower | Eigen::Upper,
                           Eigen::DiagonalPreconditioner<Scalar>>
      solver_;
  Scalar tolerance_;
  int max_iterations_;
  int iterations_ = 0;
  Scalar error_ = 0.0;
};

/// Conjugate gradients preconditioned by one smoothed-aggregation V-cycle.
class AmgCgSolver final : public LinearSolver {
 public:
  explicit AmgCgSolver(const LinearSolverOptions& options)
      : options_(options), amg_(options.amg) {}

  void set_layout(const DofLayout& layout) override {
    layout_ = layout;
    has_layout_ = true;
  }

  void factorize(const SparseMatrix& a) override {
    check_square_finite(a, name());
    if (!has_layout_) {
      throw SolverError(name() +
                        " needs the DOF layout (node coordinates of the unknowns) to "
                        "build its rigid-body near-null space; the caller did not supply "
                        "one. Use a direct solver for a bare matrix");
    }
    matrix_ = &a;
    amg_.setup(a, layout_);
  }

  Vector solve(const Vector& b) override {
    Vector x = Vector::Zero(b.size());
    return run(b, x);
  }

  Vector solve_from(const Vector& b, const Vector& initial_guess) override {
    Vector x = initial_guess.size() == b.size() ? initial_guess : Vector::Zero(b.size());
    return run(b, x);
  }

  std::string name() const override {
    return std::string("AMG-CG(") + to_string(options_.amg.smoother) + ")";
  }
  bool iterative() const override { return true; }
  int last_iterations() const override { return iterations_; }
  Scalar last_error() const override { return error_; }
  const AmgStats* amg_stats() const override { return &amg_.stats(); }
  Index storage_nonzeros() const override { return amg_.stats().storage_nonzeros; }

 private:
  Vector run(const Vector& b, Vector& x) {
    if (matrix_ == nullptr) throw SolverError(name() + ": solve called before factorize");
    const int cap = options_.max_iterations > 0 ? options_.max_iterations : 1000;
    const PcgResult result =
        pcg_solve(*matrix_, b, x, amg_, options_.iterative_tolerance, cap);
    iterations_ = result.iterations;
    error_ = result.relative_residual;
    if (!result.converged) {
      std::ostringstream os;
      os << name() << " did not converge: " << result.iterations
         << " iterations left a relative residual of " << result.relative_residual
         << " against a tolerance of " << options_.iterative_tolerance << " ("
         << amg_.stats().levels.size() << " levels, operator complexity "
         << amg_.stats().operator_complexity
         << "). Raise solver.linear.max_iterations, lower "
            "solver.linear.amg.strength_threshold, raise amg.smoother_degree, or use "
            "\"simplicial_ldlt\"";
      throw ConvergenceError(os.str());
    }
    return x;
  }

  LinearSolverOptions options_;
  AmgPreconditioner amg_;
  DofLayout layout_;
  bool has_layout_ = false;
  const SparseMatrix* matrix_ = nullptr;
  int iterations_ = 0;
  Scalar error_ = 0.0;
};

/// Direct factorisation for small systems, multigrid CG for large ones. The
/// choice is made at every factorisation from the number of unknowns and the
/// dimension of the layout; a solver of the chosen kind is kept, so a
/// multigrid hierarchy is reused across factorisations.
class AutoSolver final : public LinearSolver {
 public:
  explicit AutoSolver(const LinearSolverOptions& options) : options_(options) {}

  void set_layout(const DofLayout& layout) override {
    layout_ = layout;
    has_layout_ = true;
    if (inner_) inner_->set_layout(layout);
  }

  void factorize(const SparseMatrix& a) override {
    const int dim = has_layout_ ? layout_.dim : 2;
    // A shell or beam mesh (six DOFs per node) is a surface or a line: its
    // Cholesky fill grows like a plane problem's, not a solid's.
    const bool solid = dim == 3 && (!has_layout_ || layout_.nodal_dofs() == 3);
    const Index limit = solid ? options_.auto_direct_limit_3d : options_.auto_direct_limit_2d;
    const LinearSolverType wanted =
        (has_layout_ && static_cast<Index>(a.rows()) > limit) ? LinearSolverType::AmgCg
                                                               : LinearSolverType::SimplicialLdlt;
    if (!inner_ || wanted != chosen_) {
      LinearSolverOptions o = options_;
      o.type = wanted;
      inner_ = make_linear_solver(o);
      if (has_layout_) inner_->set_layout(layout_);
      chosen_ = wanted;
      log::debug("linear solver 'auto': ", a.rows(), " unknowns (limit ", limit, " in ",
                 dim, "-D) -> ", inner_->name());
    }
    inner_->factorize(a);
  }

  Vector solve(const Vector& b) override { return inner().solve(b); }
  Vector solve_from(const Vector& b, const Vector& x0) override {
    return inner().solve_from(b, x0);
  }
  std::string name() const override {
    return inner_ ? "auto:" + inner_->name() : std::string("auto");
  }
  bool iterative() const override { return inner_ && inner_->iterative(); }
  int last_iterations() const override { return inner_ ? inner_->last_iterations() : 0; }
  Scalar last_error() const override { return inner_ ? inner_->last_error() : 0.0; }
  const AmgStats* amg_stats() const override {
    return inner_ ? inner_->amg_stats() : nullptr;
  }
  Index storage_nonzeros() const override {
    return inner_ ? inner_->storage_nonzeros() : 0;
  }

 private:
  LinearSolver& inner() {
    if (!inner_) throw SolverError("linear solver 'auto': solve called before factorize");
    return *inner_;
  }

  LinearSolverOptions options_;
  DofLayout layout_;
  bool has_layout_ = false;
  std::unique_ptr<LinearSolver> inner_;
  LinearSolverType chosen_ = LinearSolverType::SimplicialLdlt;
};

class DenseLuSolver final : public LinearSolver {
 public:
  void factorize(const SparseMatrix& a) override {
    check_square_finite(a, name());
    if (a.rows() > 4000) {
      std::ostringstream os;
      os << name() << " refuses a system with " << a.rows()
         << " unknowns; the dense path exists for verifying small models (<= 4000 "
            "DOFs) against the sparse solvers";
      throw SolverError(os.str());
    }
    dense_ = Matrix(a);
    lu_ = dense_.partialPivLu();
    // partialPivLu has no failure flag; detect singularity from the U diagonal.
    const Matrix u = lu_.matrixLU().triangularView<Eigen::Upper>();
    const Vector diag = u.diagonal().cwiseAbs();
    if (diag.size() == 0 || diag.maxCoeff() == 0.0 ||
        diag.minCoeff() / diag.maxCoeff() < 1.0e-14) {
      std::ostringstream os;
      os << name() << ": dense LU pivot ratio "
         << (diag.maxCoeff() > 0 ? diag.minCoeff() / diag.maxCoeff() : 0.0)
         << " indicates a singular system. " << kSingularHint;
      throw SolverError(os.str());
    }
  }

  Vector solve(const Vector& b) override { return lu_.solve(b); }
  std::string name() const override { return "DenseLU"; }

 private:
  Matrix dense_;
  Eigen::PartialPivLU<Matrix> lu_;
};

}  // namespace

std::string to_string(LinearSolverType type) {
  switch (type) {
    case LinearSolverType::SimplicialLdlt: return "simplicial_ldlt";
    case LinearSolverType::SimplicialLlt: return "simplicial_llt";
    case LinearSolverType::SparseLu: return "sparse_lu";
    case LinearSolverType::ConjugateGradient: return "conjugate_gradient";
    case LinearSolverType::DenseLu: return "dense_lu";
    case LinearSolverType::AmgCg: return "amg_cg";
    case LinearSolverType::Auto: return "auto";
  }
  return "unknown";
}

LinearSolverType parse_linear_solver_type(const std::string& text) {
  std::string lower;
  std::transform(text.begin(), text.end(), std::back_inserter(lower),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  if (lower == "simplicial_ldlt" || lower == "ldlt" || lower == "cholesky")
    return LinearSolverType::SimplicialLdlt;
  if (lower == "simplicial_llt" || lower == "llt") return LinearSolverType::SimplicialLlt;
  if (lower == "sparse_lu" || lower == "lu") return LinearSolverType::SparseLu;
  if (lower == "conjugate_gradient" || lower == "cg")
    return LinearSolverType::ConjugateGradient;
  if (lower == "dense_lu" || lower == "dense") return LinearSolverType::DenseLu;
  if (lower == "amg_cg" || lower == "amg" || lower == "multigrid") return LinearSolverType::AmgCg;
  if (lower == "auto") return LinearSolverType::Auto;
  throw ConfigError("unknown linear solver '" + text +
                    "' (expected simplicial_ldlt|simplicial_llt|sparse_lu|"
                    "conjugate_gradient|dense_lu|amg_cg|auto)");
}

std::unique_ptr<LinearSolver> make_linear_solver(const LinearSolverOptions& options) {
  switch (options.type) {
    case LinearSolverType::SimplicialLdlt:
      return std::make_unique<LdltSolver>(options.pivot_tolerance);
    case LinearSolverType::SimplicialLlt:
      return std::make_unique<LltSolver>();
    case LinearSolverType::SparseLu:
      return std::make_unique<SparseLuSolver>();
    case LinearSolverType::ConjugateGradient:
      return std::make_unique<CgSolver>(options.iterative_tolerance,
                                        options.max_iterations);
    case LinearSolverType::DenseLu:
      return std::make_unique<DenseLuSolver>();
    case LinearSolverType::AmgCg:
      return std::make_unique<AmgCgSolver>(options);
    case LinearSolverType::Auto:
      return std::make_unique<AutoSolver>(options);
  }
  throw ConfigError("unhandled linear solver type");
}

Scalar scaled_residual(const SparseMatrix& a, const Vector& x, const Vector& b) {
  const Vector r = a * x - b;
  const Scalar scale = std::max(b.norm(), std::numeric_limits<Scalar>::min());
  return r.norm() / scale;
}

Scalar backward_error(const SparseMatrix& a, const Vector& x, const Vector& b) {
  const Vector r = a * x - b;
  Vector gross = b.cwiseAbs();
  for (Eigen::Index col = 0; col < a.outerSize(); ++col) {
    const Scalar xc = std::abs(x(col));
    for (SparseMatrix::InnerIterator it(a, col); it; ++it) {
      gross(it.row()) += std::abs(it.value()) * xc;
    }
  }
  return r.norm() / std::max(gross.norm(), std::numeric_limits<Scalar>::min());
}

bool residual_accepted(Scalar scaled, Scalar backward, Scalar tolerance) {
  return scaled <= tolerance || backward <= kRoundoffBackwardError;
}

Vector solve_and_verify(const SparseMatrix& a, const Vector& b,
                        const LinearSolverOptions& options, Scalar* out_residual) {
  auto solver = make_linear_solver(options);
  solver->factorize(a);
  Vector x = solver->solve(b);
  if (!x.allFinite()) {
    throw SolverError(solver->name() +
                      " returned a non-finite solution vector; " + kSingularHint);
  }
  const Scalar residual = scaled_residual(a, x, b);
  const Scalar backward = backward_error(a, x, b);
  if (out_residual) *out_residual = residual;
  if (!residual_accepted(residual, backward, options.residual_tolerance)) {
    std::ostringstream os;
    os << solver->name() << " left a scaled residual of " << residual
       << ", above the tolerance " << options.residual_tolerance
       << ", and a backward error of " << backward << ", above round-off ("
       << kRoundoffBackwardError << "). The system is likely ill-conditioned; "
       << kSingularHint;
    throw SolverError(os.str());
  }
  log::debug(solver->name(), ": solved ", a.rows(), " unknowns, scaled residual ",
             residual);
  return x;
}

bool is_positive_definite(const SparseMatrix& a, Scalar min_pivot_ratio,
                          Scalar* out_min_pivot_ratio) {
  if (a.rows() != a.cols() || a.rows() == 0) return false;
  Eigen::SimplicialLDLT<SparseMatrix, Eigen::Lower, Eigen::AMDOrdering<StorageIndex>>
      ldlt;
  ldlt.compute(a);
  if (ldlt.info() != Eigen::Success) {
    if (out_min_pivot_ratio) *out_min_pivot_ratio = 0.0;
    return false;
  }
  const Vector d = ldlt.vectorD();
  const Scalar max_abs = d.cwiseAbs().maxCoeff();
  if (!(max_abs > 0.0)) {
    if (out_min_pivot_ratio) *out_min_pivot_ratio = 0.0;
    return false;
  }
  const Scalar ratio = d.minCoeff() / max_abs;
  if (out_min_pivot_ratio) *out_min_pivot_ratio = ratio;
  return ratio > min_pivot_ratio;
}

}  // namespace sparlab
