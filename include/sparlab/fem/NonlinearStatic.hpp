/// \file NonlinearStatic.hpp
/// \brief Geometrically non-linear static analysis: large displacement and
///        rotation by the total Lagrangian formulation, solved by Newton's
///        method under load control or along the equilibrium path by the
///        arc-length method.
///
/// **Equilibrium.** With the load factor \f$\lambda\f$ and the free DOFs of
/// the displacement \f$u\f$, the residual is
/// \f[
///   R(u, \lambda) = f_{int}(u, \lambda) - \lambda f_{dead}
///                 - f_{pressure}(u, \lambda) - f_{rotation}(u, \lambda) ,
/// \f]
/// * \f$f_{int}\f$ from TotalLagrangian.hpp, with a temperature change scaled
///   by \f$\lambda\f$;
/// * \f$f_{dead}\f$ the point loads, tractions, self-weight and body forces
///   of the load case - they keep their magnitude and direction;
/// * a pressure **follows** the deforming face by default: it is integrated
///   over the current face, \f$-\lambda p\int N a(x)\f$ with the area vector
///   of the deformed geometry, and contributes the (non-symmetric) load
///   stiffness of FaceGeometry.hpp to the tangent; `follower_pressure = false`
///   keeps it on the undeformed face;
/// * a rotation's centrifugal load acts at the deformed position,
///   \f$\lambda\rho\omega^2 (I - ee^T)(X + u - c)\f$, which is linear in
///   \f$u\f$ and softens the tangent by \f$\lambda\omega^2 M_\perp\f$ (spin
///   softening);
/// * prescribed displacements are scaled by \f$\lambda\f$ (under load
///   control; the arc-length method needs them homogeneous).
///
/// **Load control** steps \f$\lambda\f$ from 0 to 1 and solves each step by
/// Newton's method with the consistent tangent from the tangent predictor,
/// optionally with a line search on the energy; a step that fails to converge
/// (or inverts an element) is cut in half and retried, and a quickly
/// converging run lengthens its steps again. Past a limit point load control
/// either fails or lands silently on a distant branch of equilibria - a
/// snap-through, which is not quasi-static - so two more kinds of step are cut:
/// one whose iterations meet a tangent with a negative pivot although it
/// started from a stable state, and one whose converged state lies farther
/// from the tangent predictor than the predicted increment itself (on a
/// smooth path that distance shrinks with the step; across a limit point it
/// does not). A run that no halving gets past stops and names the cause.
/// **Arc-length** (Crisfield's cylindrical form) constrains the norm of
/// the displacement increment instead, so it follows the path through limit
/// points - snap-through, snap-back of the load - where load control must
/// fail; the arc length adapts to the iteration count, and the root of the
/// constraint closest to the previous direction is taken.
///
/// **Convergence** of an iteration requires both
/// \f$\|R_f\| \le \epsilon_R \max(\|f_{ext}\|, \|r_p\|)\f$ (the external load
/// or, under pure displacement control, the reactions) and
/// \f$\|\delta u\| \le \epsilon_u \|\Delta u\|\f$ for the step increment - or
/// the residual at its round-off floor. Under load control a line search on
/// the energy along the Newton direction guards each iteration. Every
/// converged step records the load factor, the iterations, the residual,
/// the monitored displacements and reactions and - from the \f$LDL^T\f$
/// factor of a symmetric tangent - the number of negative pivots, which by
/// Sylvester's law of inertia counts the negative eigenvalues: a positive
/// count on a load path means the state is unstable (past a limit or
/// bifurcation point). A follower-pressure tangent that is symmetric once
/// assembled (the rim of the loaded surface held) is factorised as such.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/fem/Assembler.hpp"
#include "sparlab/fem/FemModel.hpp"
#include "sparlab/fem/StaticAnalysis.hpp"
#include "sparlab/material/Hyperelastic.hpp"

#include <limits>
#include <string>
#include <vector>

namespace sparlab {

/// A quantity to record along the path: the mean displacement component
/// over a node region, or the sum of the reaction component over it - the
/// force a displacement-controlled path is read by.
struct NonlinearMonitor {
  enum class Quantity { Displacement, Reaction };
  std::string name;
  SelectorGroup region;
  int component = 0;  ///< 0 = x, 1 = y, 2 = z
  Quantity quantity = Quantity::Displacement;
};

struct NonlinearOptions {
  HyperelasticModel law = HyperelasticModel::SaintVenantKirchhoff;
  enum class Method { LoadControl, ArcLength };
  Method method = Method::LoadControl;
  /// Load control: the number of equal steps to lambda = 1 it starts with.
  /// Arc-length: the first step is the first of `steps` equal load
  /// increments, which fixes the initial arc length.
  int steps = 10;
  int max_steps = 500;           ///< converged steps, cuts not counted
  int max_iterations = 25;       ///< Newton iterations per step
  int max_cuts = 12;             ///< successive halvings of a failing step
  Scalar residual_tolerance = 1.0e-8;
  Scalar displacement_tolerance = 1.0e-8;
  bool line_search = true;
  bool follower_pressure = true;
  /// Arc-length: the load factor at which the run stops (lambda reaches it
  /// exactly, the last step switching to load control).
  Scalar target_load_factor = 1.0;
  /// Arc-length: iterations per step the arc length adapts towards.
  int desired_iterations = 5;
  /// Arc-length: bounds on the arc length relative to the first one.
  Scalar min_arc_ratio = 1.0e-6;
  Scalar max_arc_ratio = 10.0;
  /// Load control: load factors in (0, 1), increasing, that the path must
  /// pass through exactly - each becomes a converged step (steps never jump
  /// over one), so results can be read at chosen load levels.
  std::vector<Scalar> load_factors;
  /// Quantities recorded at every converged step.
  std::vector<NonlinearMonitor> monitors;
};

std::string to_string(NonlinearOptions::Method method);
/// "load_control" or "arc_length".
NonlinearOptions::Method parse_nonlinear_method(const std::string& text);
std::string to_string(NonlinearMonitor::Quantity quantity);
/// "displacement" or "reaction".
NonlinearMonitor::Quantity parse_monitor_quantity(const std::string& text);

/// One converged step.
struct NonlinearStep {
  int index = 0;
  Scalar load_factor = 0.0;
  int iterations = 0;
  Scalar residual = 0.0;          ///< relative residual at convergence
  Scalar arc_length = 0.0;        ///< arc-length method only
  int negative_pivots = -1;       ///< -1 when not available (LU)
  int cuts = 0;                   ///< halvings before this step converged
  std::vector<Scalar> monitors;   ///< monitored displacements [m] and reactions [N]
  Scalar max_displacement = 0.0;  ///< largest nodal displacement magnitude [m]
};

struct NonlinearResult {
  std::string load_case_name;
  std::string method;
  std::string law;
  std::vector<std::string> monitor_names;
  std::vector<std::string> monitor_units;  ///< "m" or "N"
  std::vector<NonlinearStep> steps;
  Vector displacement;        ///< final, full length [m]
  Vector reactions;           ///< f_int - f_ext at the prescribed DOFs [N]
  Scalar load_factor = 0.0;   ///< final lambda
  bool completed = false;     ///< reached the target load factor
  std::string termination;    ///< why the run stopped
  /// Load control stopped at a critical point: the lowest load factor a step
  /// reached only by meeting an unstable tangent or by jumping to another
  /// branch. The limit or bifurcation point lies between `load_factor` and
  /// it. NaN otherwise.
  Scalar critical_bound = std::numeric_limits<Scalar>::quiet_NaN();
  int total_iterations = 0;
  int total_cuts = 0;
  Scalar strain_energy = 0.0; ///< int W dV0 at the final state [J]
  Scalar max_green_strain = 0.0;
  Scalar min_jacobian = 1.0;
  EquilibriumCheck equilibrium;  ///< final force balance, deformed loads
  /// The final tangent was symmetric and factorised by LDL^T (which reports
  /// its inertia); false when a follower pressure left it non-symmetric and
  /// LU factorised it.
  bool symmetric_tangent = true;
  std::string linear_solver;
  /// Final element stresses (Voigt): Cauchy and second Piola-Kirchhoff
  /// [Pa], the von Mises stress of the Cauchy stress, and sigma_zz of a
  /// plane-strain model.
  Matrix element_cauchy;
  Matrix element_piola_kirchhoff;
  Vector element_von_mises;
  Vector element_cauchy_zz;
};

/// The residual \f$R(u, \lambda)\f$, the load rate
/// \f$q = -\partial R/\partial\lambda\f$ and the tangent
/// \f$\partial R/\partial u\f$ of load case `load_case` at a state, exactly as
/// the solver evaluates them (full-length vectors, full matrix): for the
/// derivative checks of the tests and for diagnostics.
struct NonlinearState {
  Vector residual;
  Vector external;
  Vector load_rate;
  SparseMatrix tangent;
  Scalar energy = 0.0;
};
NonlinearState evaluate_nonlinear_state(const FemModel& model, const Assembler& assembler,
                                        std::size_t load_case, const NonlinearOptions& options,
                                        const Vector& u, Scalar lambda);

/// Non-linear static analysis of one load case of a finalised model.
class NonlinearStaticAnalysis {
 public:
  NonlinearStaticAnalysis(const FemModel& model, const Assembler& assembler,
                          NonlinearOptions options);

  /// Solve load case `load_case`. A run that cannot continue - a step cut
  /// `max_cuts` times, a singular tangent, the step budget exhausted - stops
  /// and reports why in `termination` with `completed = false`, returning
  /// the last converged state; it never returns an unconverged state.
  NonlinearResult solve(std::size_t load_case);

  const NonlinearOptions& options() const { return options_; }

 private:
  const FemModel& model_;
  const Assembler& assembler_;
  NonlinearOptions options_;
};

}  // namespace sparlab
