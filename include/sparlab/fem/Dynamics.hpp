/// \file Dynamics.hpp
/// \brief Transient dynamics by the HHT-alpha method and the steady-state
///        harmonic (frequency) response.
///
/// **Transient.** The semi-discrete equations of motion of a load case,
/// \f[
///   M\,\ddot u + C\,\dot u + K\,u = A(t)\, f, \qquad
///   u_p(t) = A(t)\, g,
/// \f]
/// with the load vector \f$f\f$ of the case (forces, pressures, body and
/// thermal loads) and its prescribed displacements \f$g\f$ both following the
/// amplitude \f$A(t)\f$, and Rayleigh damping \f$C = a M + b K\f$. They are
/// integrated with a constant step \f$\Delta t\f$ by the method of Hilber,
/// Hughes and Taylor (1977): the Newmark relations
/// \f[
///   u_{n+1} = u_n + \Delta t\, v_n
///     + \Delta t^2 \left[(\tfrac12 - \beta)\, a_n + \beta\, a_{n+1}\right],
///   \qquad
///   v_{n+1} = v_n + \Delta t\left[(1-\gamma)\, a_n + \gamma\, a_{n+1}\right]
/// \f]
/// with equilibrium at the weighted point
/// \f[
///   M a_{n+1} + (1+\alpha)(C v_{n+1} + K u_{n+1}) - \alpha (C v_n + K u_n)
///     = (1+\alpha) f_{n+1} - \alpha f_n ,
/// \f]
/// \f$\alpha \in [-1/3, 0]\f$, \f$\beta = (1-\alpha)^2/4\f$,
/// \f$\gamma = 1/2 - \alpha\f$. \f$\alpha = 0\f$ is the trapezoidal rule
/// (average acceleration): unconditionally stable, second order, and for a
/// linear undamped model it conserves the energy
/// \f$\tfrac12 v^T M v + \tfrac12 u^T K u\f$ exactly, the change of the energy
/// over a step being exactly the work of the trapezoidal forces - so the
/// energy balance the solver reports is a round-off check. A negative
/// \f$\alpha\f$ keeps second order and damps the modes whose period is short
/// against the step. The prescribed DOFs follow the same Newmark kinematics,
/// their motion fixed; the effective stiffness
/// \f$c_0 M + (1+\alpha) c_3 C + (1+\alpha) K\f$ is factorised once.
///
/// **Frequency response.** For a harmonic load \f$f e^{i\omega t}\f$ and a
/// harmonic prescribed displacement \f$g e^{i\omega t}\f$ the steady state
/// solves
/// \f[
///   \left[K (1 + i\eta) - \omega^2 M + i\omega C\right] U = f
/// \f]
/// on the free DOFs, a complex sparse LU per frequency, with the structural
/// (hysteretic) loss factor \f$\eta\f$ and the Rayleigh damping above.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/fem/Assembler.hpp"
#include "sparlab/fem/FemModel.hpp"
#include "sparlab/fem/LinearSolver.hpp"
#include "sparlab/fem/NonlinearStatic.hpp"
#include "sparlab/fem/Selector.hpp"

#include <complex>
#include <string>
#include <vector>

namespace sparlab {

using ComplexScalar = std::complex<Scalar>;
using ComplexVector = Eigen::VectorXcd;

/// The time function A(t) that scales a load case - its loads and its
/// prescribed displacements together.
struct Amplitude {
  enum class Kind {
    Step,      ///< A = scale for t >= 0 (the load is applied at t = 0)
    Table,     ///< piecewise linear through (times, values), constant beyond the ends
    Harmonic   ///< A = scale sin(2 pi frequency t + phase)
  };
  Kind kind = Kind::Step;
  std::vector<Scalar> times;   ///< Table: increasing [s]
  std::vector<Scalar> values;  ///< Table: one per time
  Scalar frequency = 0.0;      ///< Harmonic [Hz]
  Scalar phase = 0.0;          ///< Harmonic [rad]
  Scalar scale = 1.0;          ///< multiplies every kind

  Scalar value(Scalar t) const;
  /// dA/dt at t (a Table's slope to the right of t; 0 for a Step).
  Scalar rate(Scalar t) const;
  /// d2A/dt2 at t (0 for a Step and a Table).
  Scalar second_rate(Scalar t) const;
  /// \throws ConfigError for a malformed amplitude.
  void validate() const;
};

std::string to_string(Amplitude::Kind kind);
/// "step", "table" or "harmonic".
Amplitude::Kind parse_amplitude_kind(const std::string& text);

/// A quantity recorded at every step (transient) or frequency (harmonic):
/// the mean displacement, velocity or acceleration component over a node
/// region, or the sum of the reaction component over it.
struct DynamicMonitor {
  enum class Quantity { Displacement, Velocity, Acceleration, Reaction };
  std::string name;
  SelectorGroup region;
  int component = 0;  ///< 0 = x, 1 = y, 2 = z
  Quantity quantity = Quantity::Displacement;
};

std::string to_string(DynamicMonitor::Quantity quantity);
/// "displacement", "velocity", "acceleration" or "reaction".
DynamicMonitor::Quantity parse_dynamic_quantity(const std::string& text);

/// The Newmark parameters of the HHT-alpha method.
struct HhtParameters {
  Scalar alpha = 0.0;
  Scalar beta = 0.25;
  Scalar gamma = 0.5;
  /// \throws ConfigError unless -1/3 <= alpha <= 0.
  static HhtParameters from_alpha(Scalar alpha);
};

struct TransientOptions {
  Scalar time_step = 0.0;  ///< [s], constant
  Scalar end_time = 0.0;   ///< [s], a whole number of steps
  /// HHT-alpha in [-1/3, 0]; 0 is the trapezoidal rule.
  Scalar alpha = 0.0;
  MassType mass_type = MassType::Consistent;
  Scalar mass_damping = 0.0;       ///< Rayleigh a in C = a M + b K [1/s]
  Scalar stiffness_damping = 0.0;  ///< Rayleigh b [s]
  Amplitude amplitude;
  /// The initial state: at rest (u = 0 on the free DOFs, v = 0, the
  /// prescribed DOFs at A(0) g moving at A'(0) g), or the static solution
  /// under A(0) (v = 0) - a structure preloaded and then released.
  enum class Start { Rest, Static };
  Start start = Start::Rest;
  /// Every n-th step's displacement and velocity are kept for output (and
  /// the last step's); 0 keeps only the final state.
  int snapshot_every = 0;
  std::vector<DynamicMonitor> monitors;
  LinearSolverOptions linear;
  /// Non-linear transient: the internal forces, loads and tangent of the
  /// load case's non-linear system (`nonlinear_options`: kinematics, elastic
  /// law, plasticity with its mean dilatation, follower pressure, the Newton
  /// tolerances and iteration limit) - Newton's method at every step on the
  /// HHT-alpha residual, the plastic history committed on convergence; the
  /// Rayleigh damping takes the linear elastic stiffness. Starts at rest.
  bool nonlinear = false;
  NonlinearOptions nonlinear_options;
  /// Non-linear: a step whose Newton iteration fails is repeated in halves
  /// from the last converged state, at most this many halvings deep.
  int max_cuts = 8;
};

std::string to_string(TransientOptions::Start start);
/// "rest" or "static".
TransientOptions::Start parse_transient_start(const std::string& text);

/// One step of the path (step 0 is the initial state).
struct TransientStep {
  int index = 0;
  Scalar time = 0.0;           ///< [s]
  std::vector<Scalar> monitors;
  Scalar max_displacement = 0.0;  ///< largest nodal displacement magnitude [m]
  Scalar kinetic_energy = 0.0;    ///< v^T M v / 2 [J]
  Scalar strain_energy = 0.0;     ///< u^T K u / 2 [J]
  /// Energy dissipated by C since t = 0, the trapezoidal sum of
  /// dt/4 (v_n + v_n+1)^T C (v_n + v_n+1) [J].
  Scalar damping_energy = 0.0;
  /// Work done since t = 0 by the loads and the reactions of the prescribed
  /// motion, the trapezoidal sum of (u_n+1 - u_n)^T (F_n + F_n+1) / 2 [J].
  Scalar external_work = 0.0;
  /// Non-linear: Newton iterations over the step (its sub-steps included),
  /// the halvings it needed, the integration points that yielded in it and
  /// the largest accumulated plastic strain after it (-1 and 0 without
  /// plasticity).
  int iterations = 0;
  int cuts = 0;
  int yielding_points = -1;
  Scalar max_plastic_strain = 0.0;
};

struct TransientSnapshot {
  Scalar time = 0.0;
  Vector displacement;  ///< full length [m]
  Vector velocity;      ///< full length [m/s]
};

struct TransientResult {
  std::string load_case_name;
  std::vector<std::string> monitor_names;
  std::vector<std::string> monitor_units;
  std::vector<std::vector<Index>> monitor_nodes;  ///< the nodes each monitor covers
  std::vector<TransientStep> steps;  ///< the initial state, then every step
  HhtParameters parameters;
  Scalar time_step = 0.0;
  int num_steps = 0;
  Vector displacement;   ///< final, full length [m]
  Vector velocity;       ///< final [m/s]
  Vector acceleration;   ///< final [m/s^2]
  /// Final M a + C v + K u - f at the prescribed DOFs [N], zero elsewhere.
  Vector reactions;
  std::vector<TransientSnapshot> snapshots;
  Scalar total_mass = 0.0;  ///< [kg]
  /// The largest |E_0 + W - T - U - D| over the path, over the largest of
  /// the energies involved: exact to round-off for the trapezoidal rule
  /// (alpha = 0); with alpha < 0 the energy the method dissipates.
  Scalar energy_balance_error = 0.0;
  /// Energy the HHT-alpha method dissipated, E_0 + W - T - U - D at the end
  /// [J] (zero to round-off for alpha = 0 on a linear model; with a
  /// non-linear model it also holds the plastic dissipation and the
  /// trapezoidal rule's own error, which a linear model does not have).
  Scalar numerical_dissipation = 0.0;
  std::string linear_solver;
  std::vector<std::string> warnings;
  /// Non-linear: whether the run reached end_time, and why it stopped if
  /// not (the result then holds the last converged step).
  bool nonlinear = false;
  bool completed = true;
  std::string termination;
  int total_iterations = 0;
  bool plastic = false;
  /// Non-linear: some elastoplastic element averages its dilatation.
  bool mean_dilatation = false;
  Scalar max_plastic_strain = 0.0;
};

/// Integrate load case `load_case` of a finalised model in time.
/// \throws ConfigError for inconsistent options, ModelError for a model
///         without free DOFs, SolverError for a singular effective stiffness.
TransientResult solve_transient(const FemModel& model, const Assembler& assembler,
                                std::size_t load_case, const TransientOptions& options);

struct FrequencyResponseOptions {
  std::vector<Scalar> frequencies;  ///< [Hz], >= 0, in the order solved
  MassType mass_type = MassType::Consistent;
  Scalar mass_damping = 0.0;        ///< Rayleigh a [1/s]
  Scalar stiffness_damping = 0.0;   ///< Rayleigh b [s]
  Scalar structural_damping = 0.0;  ///< loss factor eta, K (1 + i eta)
  /// Frequencies [Hz] whose full complex field is kept for output (the
  /// nearest solved frequency each).
  std::vector<Scalar> snapshot_frequencies;
  std::vector<DynamicMonitor> monitors;
};

struct FrequencyPoint {
  Scalar frequency = 0.0;  ///< [Hz]
  /// Complex amplitudes of the monitors (displacement [m], velocity
  /// i omega U, acceleration -omega^2 U, reaction [N]).
  std::vector<ComplexScalar> monitors;
  /// The largest displacement any node reaches over a cycle [m] (see
  /// harmonic_peak_displacements).
  Scalar max_displacement = 0.0;
};

struct FrequencySnapshot {
  Scalar frequency = 0.0;
  ComplexVector displacement;  ///< full length, complex amplitude [m]
};

struct FrequencyResponseResult {
  std::string load_case_name;
  std::vector<std::string> monitor_names;
  std::vector<std::string> monitor_units;
  std::vector<std::vector<Index>> monitor_nodes;  ///< the nodes each monitor covers
  std::vector<FrequencyPoint> points;
  std::vector<FrequencySnapshot> snapshots;
  std::vector<std::string> warnings;
};

/// The largest displacement every node reaches over a cycle of the harmonic
/// motion Re(U e^{i omega t}), from the complex amplitudes U = a + i b of its
/// translations: the node traces the ellipse a cos(omega t) - b sin(omega t),
/// and this is its semi-major axis. It equals the modulus
/// sqrt(|a|^2 + |b|^2) when the components move in phase and is smaller for
/// an elliptical orbit (a circular one of radius r has modulus sqrt(2) r).
/// One value per node [m].
Vector harmonic_peak_displacements(const FemModel& model, const ComplexVector& u);

/// The steady-state harmonic response of load case `load_case` at each
/// frequency. \throws SolverError at a frequency where the dynamic stiffness
/// is singular (an undamped natural frequency).
FrequencyResponseResult solve_frequency_response(const FemModel& model,
                                                 const Assembler& assembler,
                                                 std::size_t load_case,
                                                 const FrequencyResponseOptions& options);

}  // namespace sparlab
