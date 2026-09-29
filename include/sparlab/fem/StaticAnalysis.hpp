/// \file StaticAnalysis.hpp
/// \brief Linear static solution, reaction recovery and equilibrium checks.
///
/// The reduced system is
/// \f[ K_{ff}\, u_f = f_f - K_{fp}\, u_p, \f]
/// and reactions follow from the *full* residual
/// \f[ r = K u - f, \f]
/// whose entries at prescribed DOFs are the support reactions and whose entries
/// at free DOFs are the (tiny) solver residual. This gives exact reactions with
/// no separate assembly.
///
/// Energy definitions used consistently across the code base:
///   * strain energy \f$ U = \tfrac{1}{2} u^T K u \f$ [J] - with a thermal
///     strain the elastic energy
///     \f$\tfrac12 u^T K u - u^T f_{th} + \tfrac12\int\varepsilon_0^T D\,\varepsilon_0\f$,
///   * compliance \f$ C = f^T u \f$ [J] with \f$f\f$ the *applied* load vector
///     (thermal equivalent loads included).
/// For homogeneous Dirichlet data and no thermal strain \f$C = 2U\f$; the two
/// are reported separately so the identity can be checked (it is, in the test
/// suite).
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/fem/Assembler.hpp"
#include "sparlab/fem/FemModel.hpp"
#include "sparlab/fem/LinearSolver.hpp"

#include <string>
#include <vector>

namespace sparlab {

/// Global force / moment balance of a solved load case. Moments are taken
/// about the origin; on a 2-D model only the z component is non-zero. The
/// relative errors divide by the gross size of the applied nodal loads,
/// \f$\sum_n |f_n|\f$ and \f$\sum_n |x_n|\,|f_n|\f$ (plus nodal moments), the
/// scale of the round-off in the sums, so a self-equilibrated load - a
/// thermal strain, a self-weight carried by a traction - is judged fairly.
struct EquilibriumCheck {
  Vector3 applied_force = Vector3::Zero();    ///< sum of applied nodal forces [N]
  Vector3 reaction_force = Vector3::Zero();   ///< sum of support reactions [N]
  Vector3 force_residual = Vector3::Zero();   ///< applied + reaction [N]
  Vector3 applied_moment = Vector3::Zero();   ///< about the origin [N m]
  Vector3 reaction_moment = Vector3::Zero();  ///< about the origin [N m]
  Vector3 moment_residual = Vector3::Zero();  ///< applied + reaction [N m]
  Scalar relative_force_error = 0.0;          ///< |residual| / sum_n |f_n|
  Scalar relative_moment_error = 0.0;         ///< |residual| / sum_n |x_n| |f_n|
};

/// Result of one solved load case.
struct StaticSolution {
  std::string load_case_name;
  Scalar weight = 1.0;
  Vector displacement;                ///< full length, [m]
  Vector reactions;                   ///< full length, non-zero at prescribed DOFs [N]
  Scalar compliance = 0.0;            ///< f^T u [J]
  Scalar strain_energy = 0.0;         ///< 0.5 u^T K u [J]
  Scalar max_displacement_magnitude = 0.0;  ///< [m]
  Index max_displacement_node = -1;
  Scalar scaled_residual = 0.0;       ///< ||K_ff u_f - rhs|| / ||rhs||
  Scalar backward_error = 0.0;        ///< ||K_ff u_f - rhs|| / || |K_ff||u_f| + |rhs| ||
  int solver_iterations = 0;
  std::string solver_name;            ///< the linear solver that produced it
  EquilibriumCheck equilibrium;
};

struct StaticAnalysisOptions {
  LinearSolverOptions linear;
  /// Run `require_well_posed` before assembling. Disabled internally by the
  /// topology optimizer, whose SIMP stiffness floor keeps K non-singular even
  /// in void regions.
  bool check_model = true;
  /// Relative force-balance tolerance; exceeding it raises SolverError.
  Scalar equilibrium_tolerance = 1.0e-6;
};

/// Linear static solver. One instance reuses a single factorisation of
/// \f$K_{ff}\f$ across all load cases, which is what makes multi-load-case
/// topology optimisation affordable.
class StaticAnalysis {
 public:
  StaticAnalysis(const FemModel& model, const Assembler& assembler,
                 StaticAnalysisOptions options = StaticAnalysisOptions());

  /// Solve with `solver`, owned by the caller, instead of a private solver.
  /// A topology optimisation keeps one solver across design iterations so a
  /// multigrid hierarchy is reused. The solver must outlive this analysis.
  void use_external_solver(LinearSolver* solver);

  /// Factorise \f$K_{ff}\f$ for the given per-element stiffness factors.
  /// \param stiffness_scale optional \f$E(\rho_e)/E_0\f$ (length num_elements).
  void prepare(const Vector* stiffness_scale = nullptr);

  /// Solve every load case of the model. `prepare` is called automatically if
  /// it has not been called since the last change of `stiffness_scale`.
  std::vector<StaticSolution> solve_all(const Vector* stiffness_scale = nullptr);

  /// Solve one right-hand side with the current factorisation, returning the
  /// full displacement vector. `initial_guess` (full length) seeds an
  /// iterative solver and is ignored by a direct one.
  Vector solve_load_vector(const Vector& applied_force,
                           const Vector* initial_guess = nullptr);

  /// Solve \f$K_{ff}\,\lambda_f = r_f\f$ with the current factorisation and
  /// zeros at the prescribed DOFs: the adjoint problem of any functional
  /// evaluated on this model, whatever the prescribed displacements.
  Vector solve_homogeneous(const Vector& rhs, const Vector* initial_guess = nullptr);

  /// Iterations of the last solve (0 for a direct solver).
  int last_iterations() const;

  /// The solver in use (after `prepare`).
  const LinearSolver& solver() const;

  /// Where the unknowns of \f$K_{ff}\f$ live, for the multigrid solver.
  const DofLayout& layout() const { return layout_; }

  /// Weighted compliance \f$\sum_l w_l\, f_l^T u_l\f$ with normalised weights.
  static Scalar weighted_compliance(const std::vector<StaticSolution>& solutions,
                                    const std::vector<Scalar>& normalised_weights);

  /// The assembled global stiffness matrix of the last `prepare` call.
  const SparseMatrix& stiffness() const { return k_full_; }

  const StaticAnalysisOptions& options() const { return options_; }

 private:
  StaticSolution build_solution(const std::string& name, Scalar weight,
                                const Vector& applied_force,
                                const LoadCaseData* data = nullptr);

  LinearSolver& active_solver();
  Vector free_guess(const Vector* initial_guess) const;

  const FemModel& model_;
  const Assembler& assembler_;
  StaticAnalysisOptions options_;
  DofLayout layout_;
  std::unique_ptr<LinearSolver> owned_solver_;
  LinearSolver* external_solver_ = nullptr;
  SparseMatrix k_full_;
  SparseMatrix k_ff_;
  SparseMatrix k_fp_;
  Vector prescribed_;
  bool prepared_ = false;
};

}  // namespace sparlab
