/// \file LinearSolver.hpp
/// \brief Sparse and dense linear solvers with explicit failure reporting.
///
/// Every solver reports failure through `SolverError`; none of them returns a
/// silently wrong answer. After each solve the caller can (and by default does)
/// verify the scaled residual
/// \f$ \|A x - b\|_2 / \max(\|b\|_2, \epsilon) \f$ against a recorded
/// tolerance, or - where forming \f$Ax - b\f$ cancels so much that no
/// solver can meet that tolerance in double precision - a backward error at
/// round-off (`residual_accepted`). Recorded tolerances appear in every run
/// summary.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/fem/Multigrid.hpp"

#include <limits>
#include <memory>
#include <string>

namespace sparlab {

enum class LinearSolverType {
  SimplicialLdlt,     ///< sparse Cholesky (LDL^T) with AMD ordering - default
  SimplicialLlt,      ///< sparse Cholesky (LL^T), requires strict SPD
  SparseLu,           ///< sparse LU, works for indefinite systems
  ConjugateGradient,  ///< diagonal-preconditioned CG (matrix-free friendly)
  DenseLu,            ///< dense partial-pivot LU; verification of small models
  AmgCg,              ///< CG preconditioned by smoothed-aggregation multigrid
  Auto                ///< SimplicialLdlt up to a size limit, AmgCg above it
};

std::string to_string(LinearSolverType type);
LinearSolverType parse_linear_solver_type(const std::string& text);

struct LinearSolverOptions {
  /// Auto: an exact sparse Cholesky factorisation up to the size limits
  /// below, multigrid CG above them.
  LinearSolverType type = LinearSolverType::Auto;
  /// Relative residual tolerance requested from iterative solvers.
  Scalar iterative_tolerance = 1.0e-12;
  /// Iteration cap for iterative solvers; 0 selects Eigen's default (2n).
  int max_iterations = 0;
  /// Scaled-residual tolerance enforced after every solve. Exceeding it raises
  /// SolverError rather than returning a dubious displacement field.
  Scalar residual_tolerance = 1.0e-8;
  /// Smallest acceptable |pivot| / max|pivot| ratio in the Cholesky factor.
  /// Below this the system is reported as numerically singular.
  Scalar pivot_tolerance = 1.0e-14;
  /// Multigrid parameters (AmgCg, and Auto above its limit).
  AmgOptions amg;
  /// Start iterative solves from the previous solution of the same problem
  /// (the last design iteration, the last subspace iterate) where a caller
  /// has one. Direct solvers ignore it.
  bool warm_start = true;
  /// Auto: the largest number of unknowns still factorised directly, for
  /// plane (2-D) and solid (3-D) models. Above it, multigrid CG is used.
  Index auto_direct_limit_2d = 50000;
  Index auto_direct_limit_3d = 10000;
};

/// Polymorphic wrapper over the Eigen solvers.
class LinearSolver {
 public:
  virtual ~LinearSolver() = default;

  /// Factorise (or store, for iterative solvers) the system matrix.
  /// \throws SolverError on numerical failure, with a diagnosis of the likely
  ///         modelling cause.
  virtual void factorize(const SparseMatrix& a) = 0;

  /// Solve for one right-hand side. `factorize` must have been called.
  virtual Vector solve(const Vector& b) = 0;

  /// Solve starting from `initial_guess`. Iterative solvers use it; direct
  /// solvers ignore it and give exactly the `solve(b)` answer.
  virtual Vector solve_from(const Vector& b, const Vector& initial_guess) {
    (void)initial_guess;
    return solve(b);
  }

  /// Where the unknowns live (their nodes and coordinates). The multigrid
  /// solver needs it to build the rigid-body near-null space; it must be set
  /// before `factorize` and stay valid while the solver is used. Other
  /// solvers ignore it.
  virtual void set_layout(const DofLayout& layout) { (void)layout; }

  virtual std::string name() const = 0;

  /// True when solutions are approximate to a residual tolerance.
  virtual bool iterative() const { return false; }

  /// Iterations consumed by the last solve (0 for direct solvers).
  virtual int last_iterations() const { return 0; }

  /// Solver-reported error estimate of the last solve (0 for direct solvers).
  virtual Scalar last_error() const { return 0.0; }

  /// Multigrid hierarchy statistics of the last factorisation, or nullptr.
  virtual const AmgStats* amg_stats() const { return nullptr; }

  /// Numbers the solver stores beyond the matrix itself: the factor's
  /// non-zeros for a Cholesky factorisation, the coarse operators,
  /// transfer operators and coarse factor of a multigrid hierarchy. A
  /// measure of solver memory (8 bytes of value plus 4 of index each,
  /// roughly). 0 when not tracked.
  virtual Index storage_nonzeros() const { return 0; }
};

std::unique_ptr<LinearSolver> make_linear_solver(const LinearSolverOptions& options);

/// Scaled residual \f$\|Ax-b\| / \max(\|b\|, \text{tiny})\f$.
Scalar scaled_residual(const SparseMatrix& a, const Vector& x, const Vector& b);

/// Backward error of a solution: the residual over the norm of its gross
/// terms, \f$\|Ax - b\| / \|\,|A||x| + |b|\,\|\f$. The rounding of the sums
/// that form \f$Ax - b\f$ is a small multiple of \f$\epsilon\f$ times that
/// denominator, so no solver - direct or iterative - can leave a smaller
/// backward error than a few \f$\epsilon\f$; a backward-stable solve (the
/// Cholesky factorisation of an SPD matrix) leaves about that. It never
/// exceeds the scaled residual. The two part company in a bending-dominated
/// structure: in a slender beam the entries of \f$Ax\f$ are sums of terms far
/// larger than their result, so the Cholesky solve of a Tet10 cantilever with
/// L/h = 100 leaves a scaled residual of 6.5e-7 and a backward error of
/// 2.4e-16.
Scalar backward_error(const SparseMatrix& a, const Vector& x, const Vector& b);

/// The backward error at or below which a residual is round-off.
inline constexpr Scalar kRoundoffBackwardError = 64.0 * std::numeric_limits<Scalar>::epsilon();

/// A solve is accepted when its scaled residual is within `tolerance` or its
/// backward error is at round-off (`kRoundoffBackwardError`): the solution is
/// then the exact one of a system whose entries differ from the assembled
/// ones in the last bits, the best any floating-point solve can deliver.
bool residual_accepted(Scalar scaled, Scalar backward, Scalar tolerance);

/// Factorise, solve and verify the residual in one call.
/// \throws SolverError when the residual is not accepted (`residual_accepted`)
///         or the solution is not finite.
Vector solve_and_verify(const SparseMatrix& a, const Vector& b,
                        const LinearSolverOptions& options,
                        Scalar* out_residual = nullptr);

/// Test whether a symmetric sparse matrix is numerically positive definite by
/// inspecting the LDL^T pivots.
/// \param min_pivot_ratio smallest accepted pivot / largest pivot.
bool is_positive_definite(const SparseMatrix& a, Scalar min_pivot_ratio = 1.0e-14,
                          Scalar* out_min_pivot_ratio = nullptr);

}  // namespace sparlab
