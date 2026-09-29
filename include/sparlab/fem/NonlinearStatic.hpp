/// \file NonlinearStatic.hpp
/// \brief Non-linear static analysis: large displacement and rotation by the
///        total Lagrangian formulation and J2 plasticity, solved by Newton's
///        method under load control (along a load path that may unload) or
///        along the equilibrium path by the arc-length method.
///
/// **Kinematics.** `finite` (the default) is the total Lagrangian
/// formulation of TotalLagrangian.hpp: large displacement and rotation.
/// `small_strain` keeps the linear strain \f$\varepsilon = Bu\f$ and the
/// undeformed geometry - no geometric stiffness, pressures on the undeformed
/// faces, the centrifugal load of the undeformed body - so that with elastic
/// materials it is the linear analysis and with plastic ones the classical
/// materially-non-linear-only analysis. The solver reports the rotation and
/// the neglected quadratic part of the Green strain, and warns when they are
/// not small.
///
/// **Plasticity.** A material with a yield stress (IsotropicMaterial.hpp) is
/// elastoplastic: its elements take the J2 return of Plasticity.hpp at every
/// integration point (Elastoplastic.hpp), in the linear strain or, with
/// finite kinematics, in the Green-Lagrange strain. Each point keeps the
/// internal variables of the last converged step; every iteration returns
/// from them (not from the last iterate), and a converged step commits the
/// new ones, so a rejected or halved step leaves no trace. The Q4 and Hex8
/// elements of an elastoplastic material use the mean dilatation (B-bar) by
/// default, which keeps them from locking under the isochoric plastic flow
/// (MeanDilatation).
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
/// **Load control** steps \f$\lambda\f$ from 0 to 1 - or along a `load_path`
/// of turning points, such as \f$0 \to 1 \to 0\f$ to load and unload - and
/// solves each step by
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
/// does not). A run that no halving gets past stops and names the cause. The
/// second test is made only for elastic materials with finite kinematics:
/// the elastic predictor of a step in which points start to yield
/// underestimates the increment by the ratio of elastic to plastic
/// stiffness, and small-strain kinematics has no distant branches to jump
/// to. Beyond a plastic collapse load (no hardening) there is no equilibrium
/// at all, and load control stops with the last converged load below it.
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
#include "sparlab/fem/Contact.hpp"
#include "sparlab/fem/Elastoplastic.hpp"
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

/// Which elements of an elastoplastic material average their dilatation (see
/// Elastoplastic.hpp); elastic materials always keep their plain operator,
/// and plane stress and one-point elements have nothing to average.
enum class MeanDilatation {
  /// Q4 and Hex8, whose full integration locks under isochoric plastic flow
  /// (on the thick-cylinder verification the fully integrated Q4's collapse
  /// load is 1.2 % high at 16 cells through the wall and its plateau keeps
  /// rising). Not Tet10: its four-point element converges to the collapse
  /// load at third order without it, and at second order from below with it,
  /// its constant element pressure oscillating (verify_plasticity.cpp).
  Auto,
  All,   ///< every element with several integration points
  None
};

struct NonlinearOptions {
  /// The elastic law with finite kinematics (small strain is linear
  /// elasticity); an elastoplastic material takes the Saint Venant-Kirchhoff
  /// form.
  HyperelasticModel law = HyperelasticModel::SaintVenantKirchhoff;
  Kinematics kinematics = Kinematics::Finite;
  MeanDilatation mean_dilatation = MeanDilatation::Auto;
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
  /// Load control: the load factors to visit in turn from 0, each reached
  /// exactly - a load path with turning points, such as {1, 0} (load, then
  /// unload) or {1, -1, 1} (a cycle). Each leg starts with `steps` equal
  /// steps. Empty: from 0 to 1. Not with `load_factors` or arc length.
  std::vector<Scalar> load_path;
  /// Quantities recorded at every converged step.
  std::vector<NonlinearMonitor> monitors;
  /// Unilateral contact (Contact.hpp): with `small_strain` kinematics and
  /// load control only.
  ContactOptions contact;
};

std::string to_string(MeanDilatation mode);
/// "auto", "all" or "none".
MeanDilatation parse_mean_dilatation(const std::string& text);
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
  /// Elastoplastic models: the integration points that yielded in this step
  /// and the largest accumulated plastic strain after it; -1 and 0 otherwise.
  int yielding_points = -1;
  Scalar max_plastic_strain = 0.0;
  /// Contact: per pair, the slave nodes in contact and the resultant
  /// contact force on the slave body [N]; empty without contact.
  std::vector<int> contact_active;
  std::vector<Vector3> contact_force;
};

struct NonlinearResult {
  std::string load_case_name;
  std::string method;
  std::string law;
  std::string kinematics;
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
  /// Load control stopped because no step converged beyond this load factor
  /// (in `max_cuts` halvings, successes short of it notwithstanding): a limit
  /// load - a limit point, or a plastic collapse without hardening - may lie
  /// between `load_factor` and it, or Newton's method may merely have failed
  /// there. NaN otherwise.
  Scalar unreached_load_factor = std::numeric_limits<Scalar>::quiet_NaN();
  int total_iterations = 0;
  int total_cuts = 0;
  Scalar strain_energy = 0.0; ///< stored energy at the final state (elastic + hardening) [J]
  /// The largest strain component: Green-Lagrange with finite kinematics,
  /// the linear strain with small strain [-].
  Scalar max_green_strain = 0.0;
  Scalar min_jacobian = 1.0;
  /// Small strain: the largest infinitesimal rotation [rad] and the largest
  /// component of the neglected quadratic strain H^T H / 2 [-].
  Scalar max_rotation = 0.0;
  Scalar max_quadratic_strain = 0.0;
  /// Elastoplastic models: whether any material yields at all, whether mean
  /// dilatation was applied, the largest accumulated plastic strain of each
  /// element's points, and how many of the elastoplastic elements'
  /// integration points have yielded.
  bool plastic = false;
  bool mean_dilatation = false;
  Vector element_plastic_strain;
  Scalar max_plastic_strain = 0.0;
  int plastic_points = 0;
  int total_points = 0;
  /// Conditions the run met that bear on its validity (also logged).
  std::vector<std::string> warnings;
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
  /// Contact at the final state: every slave node that takes part, and the
  /// pair totals; empty without contact. The contact forces of rigid
  /// obstacles count as reactions in `equilibrium`.
  std::vector<ContactNodeResult> contact_nodes;
  std::vector<ContactPairResult> contact_pairs;
};

/// The residual \f$R(u, \lambda)\f$, the load rate
/// \f$q = -\partial R/\partial\lambda\f$ and the tangent
/// \f$\partial R/\partial u\f$ of load case `load_case` at a state, exactly as
/// the solver evaluates them (full-length vectors, full matrix): for the
/// derivative checks of the tests and for diagnostics. Elastoplastic points
/// return from the virgin state.
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
