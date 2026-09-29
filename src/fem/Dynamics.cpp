/// \file Dynamics.cpp
/// \brief HHT-alpha transient integration and the steady-state harmonic
///        response (see Dynamics.hpp for the formulation).
#include "sparlab/fem/Dynamics.hpp"

#include "NonlinearSystem.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"
#include "sparlab/fem/Multigrid.hpp"

#include <Eigen/SparseLU>

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

namespace sparlab {
namespace {

constexpr Scalar kPi = 3.14159265358979323846;

using ComplexSparse = Eigen::SparseMatrix<ComplexScalar, Eigen::ColMajor, StorageIndex>;

Vector restrict(const Vector& full, const std::vector<Index>& dofs) {
  Vector out(static_cast<Eigen::Index>(dofs.size()));
  for (std::size_t i = 0; i < dofs.size(); ++i) out(static_cast<Eigen::Index>(i)) = full(dofs[i]);
  return out;
}

void scatter(Vector& full, const std::vector<Index>& dofs, const Vector& part) {
  for (std::size_t i = 0; i < dofs.size(); ++i) full(dofs[i]) = part(static_cast<Eigen::Index>(i));
}

/// The nodes of every monitor, checked against the model.
std::vector<std::vector<Index>> monitor_nodes(const FemModel& model,
                                              const std::vector<DynamicMonitor>& monitors) {
  std::vector<std::vector<Index>> out;
  for (const DynamicMonitor& m : monitors) {
    if (m.component < 0 || m.component >= model.dim()) {
      throw ConfigError("monitor '" + m.name + "' asks for a component the model lacks");
    }
    std::vector<Index> nodes = m.region.select_nodes(model.mesh());
    if (nodes.empty()) throw ConfigError("monitor '" + m.name + "' selects no node");
    out.push_back(std::move(nodes));
  }
  return out;
}

std::string monitor_unit(DynamicMonitor::Quantity q) {
  switch (q) {
    case DynamicMonitor::Quantity::Displacement: return "m";
    case DynamicMonitor::Quantity::Velocity: return "m/s";
    case DynamicMonitor::Quantity::Acceleration: return "m/s^2";
    case DynamicMonitor::Quantity::Reaction: return "N";
  }
  return "";
}

/// A monitor's value: the mean of a kinematic component over its nodes, or
/// the sum of the reaction component.
template <typename VectorType>
typename VectorType::Scalar monitor_value(const FemModel& model, const DynamicMonitor& m,
                                          const std::vector<Index>& nodes,
                                          const VectorType& u, const VectorType& v,
                                          const VectorType& a, const VectorType& r) {
  const VectorType* source = &u;
  switch (m.quantity) {
    case DynamicMonitor::Quantity::Displacement: source = &u; break;
    case DynamicMonitor::Quantity::Velocity: source = &v; break;
    case DynamicMonitor::Quantity::Acceleration: source = &a; break;
    case DynamicMonitor::Quantity::Reaction: source = &r; break;
  }
  typename VectorType::Scalar sum(0);
  for (Index node : nodes) sum += (*source)(model.dofs().dof(node, m.component));
  if (m.quantity == DynamicMonitor::Quantity::Reaction) return sum;
  return sum / static_cast<Scalar>(nodes.size());
}

Scalar max_nodal_magnitude(const FemModel& model, const Vector& u) {
  const int dim = model.dim();
  Scalar top = 0.0;
  for (Index node = 0; node < model.mesh().num_nodes(); ++node) {
    Scalar s = 0.0;
    for (int k = 0; k < dim; ++k) {
      const Scalar x = u(model.dofs().dof(node, k));
      s += x * x;
    }
    top = std::max(top, std::sqrt(s));
  }
  return top;
}


/// Mass of the model: every translational direction carries it once, so the
/// assembled matrix sums to dim times it.
Scalar model_mass(const FemModel& model, const SparseMatrix& m_full) {
  if (model.dofs_per_node() == model.dim()) {
    return m_full.sum() / static_cast<Scalar>(model.dim());
  }
  Vector t = Vector::Zero(model.dofs().num_dofs());
  for (Index node = 0; node < model.mesh().num_nodes(); ++node) t(model.dofs().dof(node, 0)) = 1.0;
  return t.dot(m_full * t);
}

DofLayout layout_of(const FemModel& model) {
  DofLayout layout;
  layout.dim = model.dim();
  layout.dofs_per_node = model.dofs_per_node();
  layout.coordinates = &model.mesh().coordinates();
  layout.unknowns = &model.dofs().free_dofs();
  return layout;
}

void check_load_case(const FemModel& model, std::size_t load_case) {
  if (!model.finalized()) throw ConfigError("the model must be finalised before a dynamic analysis");
  if (load_case >= model.load_case_specs().size()) {
    throw ConfigError("dynamic analysis: load case index out of range");
  }
  if (model.dofs().free_dofs().empty()) {
    throw ModelError("the dynamic analysis has no free degrees of freedom; every DOF is "
                     "constrained");
  }
}

/// The non-linear transient (see TransientOptions::nonlinear): Newton's
/// method on the HHT-alpha residual
///   R(u1) = M a1 + (1 + alpha) (C v1 + r(u1, A1)) - alpha (C v0 + r(u0, A0)),
/// r = f_int - f_ext of the non-linear system, a1 and v1 the Newmark
/// functions of u1, with the tangent c0 M + (1 + alpha) (c3 C + K_T). A step
/// whose iteration fails is finished in halves from the last converged state
/// (whose plastic history is already committed), and the halves grow back.
TransientResult nonlinear_transient(const FemModel& model, const Assembler& assembler,
                                    std::size_t load_case, const TransientOptions& options,
                                    int steps, const HhtParameters& p, TransientResult result,
                                    const std::vector<std::vector<Index>>& nodes,
                                    const SparseMatrix& m, const SparseMatrix& k0) {
  if (options.start != TransientOptions::Start::Rest) {
    throw ConfigError("the non-linear transient starts at rest; a preloaded state is a "
                      "non-linear static analysis of its own (give the load an amplitude "
                      "that ramps it up instead)");
  }
  if (options.max_cuts < 0) throw ConfigError("transient: max_cuts must be >= 0");
  const NonlinearOptions& nl = options.nonlinear_options;
  if (nl.max_iterations < 1 || !(nl.residual_tolerance > 0.0) ||
      !(nl.displacement_tolerance > 0.0)) {
    throw ConfigError("transient: the non-linear tolerances must be positive and "
                      "max_iterations >= 1");
  }
  detail::NonlinearSystem system(model, assembler, load_case, nl);
  result.nonlinear = true;
  result.plastic = system.plastic();

  const DofManager& dofs = model.dofs();
  const std::vector<Index>& free = dofs.free_dofs();
  const std::vector<Index>& fixed = dofs.constrained_dofs();
  const Index n = dofs.num_dofs();
  const bool damped = options.mass_damping > 0.0 || options.stiffness_damping > 0.0;
  SparseMatrix c(n, n);
  if (damped) c = options.mass_damping * m + options.stiffness_damping * k0;
  const Vector g_p = restrict(dofs.prescribed_vector(), fixed);
  const Amplitude& amp = options.amplitude;
  const Scalar dt = options.time_step;
  const Scalar w = 1.0 + p.alpha;
  const Scalar eps = std::numeric_limits<Scalar>::epsilon();

  // The initial state at rest, the prescribed DOFs moving with the amplitude.
  Scalar lambda = amp.value(0.0);
  Vector u = Vector::Zero(n);
  Vector v = Vector::Zero(n);
  Vector a = Vector::Zero(n);
  scatter(u, fixed, lambda * g_p);
  scatter(v, fixed, amp.rate(0.0) * g_p);
  scatter(a, fixed, amp.second_rate(0.0) * g_p);
  detail::Evaluation ev = system.evaluate(u, lambda, false);
  {
    // M_ff a_f = -(M a_p + C v + f_int - f_ext)_f.
    Vector known = ev.residual + m * a;
    if (damped) known += c * v;
    const SparseMatrix m_ff = assembler.reduce_free_free(m);
    std::unique_ptr<LinearSolver> mass = make_linear_solver(options.linear);
    mass->set_layout(layout_of(model));
    mass->factorize(m_ff);
    scatter(a, free, mass->solve(-restrict(known, free)));
  }

  // The force on the body: the loads, and at the prescribed DOFs the
  // reactions M a + C v + f_int - f_ext on top of them.
  const auto body_force = [&](const detail::Evaluation& e, const Vector& vv, const Vector& aa,
                              Vector& reactions) {
    Vector force = e.external;
    const Vector inertia = m * aa;
    Vector damping;
    if (damped) damping = c * vv;
    reactions = Vector::Zero(n);
    for (Index d : fixed) {
      const Scalar total = inertia(d) + (damped ? damping(d) : 0.0) + e.internal(d);
      reactions(d) = total - e.external(d);
      force(d) = total;
    }
    return force;
  };
  Vector r;
  Vector force = body_force(ev, v, a, r);
  Vector residual = ev.residual;  // f_int - f_ext of the converged state
  Scalar strain_energy = ev.energy;
  const Scalar e0 = 0.5 * v.dot(m * v) + strain_energy;
  Scalar damping_energy = 0.0;
  Scalar work = 0.0;
  Scalar energy_scale = std::abs(e0);
  Scalar worst_balance = 0.0;
  bool energy_warned = false;

  const auto record = [&](int index, Scalar t, int iterations, int cuts, int yielding) {
    TransientStep s;
    s.index = index;
    s.time = t;
    for (std::size_t i = 0; i < options.monitors.size(); ++i) {
      s.monitors.push_back(monitor_value(model, options.monitors[i], nodes[i], u, v, a, r));
    }
    s.max_displacement = max_nodal_magnitude(model, u);
    s.kinetic_energy = 0.5 * v.dot(m * v);
    s.strain_energy = strain_energy;
    s.damping_energy = damping_energy;
    s.external_work = work;
    s.iterations = iterations;
    s.cuts = cuts;
    if (system.plastic()) {
      s.yielding_points = yielding;
      s.max_plastic_strain = system.max_plastic_strain();
    }
    energy_scale = std::max({energy_scale, s.kinetic_energy, std::abs(s.strain_energy),
                             damping_energy, std::abs(work)});
    const Scalar balance = e0 + work - s.kinetic_energy - s.strain_energy - damping_energy;
    worst_balance = std::max(worst_balance, std::abs(balance));
    // The balance holds the plastic dissipation (positive) and the method's
    // own error, small on a resolved path. Energy created - or, without
    // plasticity, lost - beyond a percent of the energies involved means the
    // step is too long for the non-linear response or the state has jumped
    // to a spurious branch (a load applied suddenly at a node can crush the
    // loaded element).
    // (With alpha < 0 the trapezoidal work is not the method's own energy
    // measure, so only the trapezoidal rule is checked.)
    if (!energy_warned && index > 0 && p.alpha == 0.0 &&
        (balance < -1.0e-2 * energy_scale ||
         (!system.plastic() && balance > 1.0e-2 * energy_scale))) {
      std::ostringstream os;
      os << "the energy balance of the non-linear transient is off by " << balance << " J ("
         << balance / energy_scale << " of the energies involved) at t = " << t
         << " s: the step is too long for the non-linear response, or the state has jumped "
            "to a spurious branch - shorten time_step, or apply the load gradually (a table "
            "amplitude) rather than suddenly";
      result.warnings.push_back(os.str());
      log::warn(os.str());
      energy_warned = true;
    }
    result.total_iterations += iterations;
    result.steps.push_back(std::move(s));
  };
  const auto snapshot = [&](Scalar t) {
    TransientSnapshot snap;
    snap.time = t;
    snap.displacement = u;
    snap.velocity = v;
    result.snapshots.push_back(std::move(snap));
  };
  record(0, 0.0, 0, 0, 0);
  if (options.snapshot_every > 0) snapshot(0.0);

  detail::TangentFactor factor;
  // One step of length h from the converged state: true on convergence,
  // with the new state and its evaluation in the out-parameters.
  const auto advance = [&](Scalar t, Scalar h, Vector& u1, Vector& v1, Vector& a1,
                           detail::Evaluation& ev1, int& iterations) -> bool {
    const Scalar c0 = 1.0 / (p.beta * h * h);
    const Scalar c1 = 1.0 / (p.beta * h);
    const Scalar c2 = 1.0 / (2.0 * p.beta) - 1.0;
    const Scalar c3 = p.gamma / (p.beta * h);
    const Scalar c4 = 1.0 - p.gamma / p.beta;
    const Scalar c5 = h * (1.0 - p.gamma / (2.0 * p.beta));
    const Scalar lambda1 = amp.value(t + h);
    const Vector a_hat = -c0 * u - c1 * v - c2 * a;
    const Vector v_hat = -c3 * u + c4 * v + c5 * a;
    Vector previous = p.alpha * residual;
    if (damped) previous += p.alpha * (c * v);
    // Predictor: constant acceleration over the step on the free DOFs.
    u1 = u + h * v + 0.5 * h * h * a;
    scatter(u1, fixed, lambda1 * g_p);
    Scalar last_correction = 0.0;
    for (int it = 1; it <= nl.max_iterations; ++it) {
      try {
        ev1 = system.evaluate(u1, lambda1, true);
      } catch (const SolverError& error) {
        log::debug("  transient step at t = ", t + h, ": ", error.what());
        return false;
      }
      a1 = c0 * u1 + a_hat;
      v1 = c3 * u1 + v_hat;
      const Vector inertia = m * a1;
      Vector total = inertia + w * ev1.residual - previous;
      Vector damping_force;
      if (damped) {
        damping_force = c * v1;
        total += w * damping_force;
      }
      const Vector rf = restrict(total, free);
      if (!rf.allFinite()) return false;
      SparseMatrix jacobian = c0 * m + w * ev1.tangent;
      if (damped) jacobian += (w * c3) * c;
      Scalar reactions = 0.0;
      for (Index d : fixed) reactions += total(d) * total(d);
      const Scalar scale = std::max({restrict(inertia, free).norm(), ev1.external.norm(),
                                     damped ? restrict(damping_force, free).norm() : 0.0,
                                     std::sqrt(reactions), ev1.thermal, 1.0e-300});
      const Scalar floor =
          detail::residual_floor(ev1, detail::stiffness_gross(jacobian, u1, free), scale);
      // The correction is judged against the step's increment or the
      // displacement itself, whichever is larger: at a turning point of the
      // motion a step barely moves, and a limit relative to its increment
      // alone would ask for corrections below the round-off of u.
      const Scalar increment = restrict(u1 - u, free).norm();
      const Scalar limit = std::max(
          nl.displacement_tolerance * std::max({increment, restrict(u1, free).norm(), 1.0e-300}),
          64.0 * eps * u1.norm());
      log::debug("  transient newton ", it, " at t = ", t + h, ": |R| = ", rf.norm(), " (scale ",
                 scale, ", round-off floor ", floor, "), |du| = ", last_correction, " (limit ",
                 limit, ")");
      if ((rf.norm() <= nl.residual_tolerance * scale && last_correction <= limit) ||
          (it > 1 && rf.norm() <= floor)) {
        iterations = it - 1;
        return true;
      }
      if (!factor.factorize(detail::free_block(jacobian, free), system.symmetric())) {
        return false;
      }
      const Vector du = factor.solve(-rf);
      if (!du.allFinite()) return false;
      detail::add_to(u1, free, du, 1.0);
      last_correction = du.norm();
    }
    return false;
  };

  Scalar t = 0.0;
  Scalar h = dt;
  for (int step = 1; step <= steps && result.completed; ++step) {
    const Scalar target = step * dt;
    int iterations = 0;
    int cuts = 0;
    int yielding = 0;
    while (t < target - 1.0e-9 * dt) {
      h = std::min(h, target - t);
      Vector u1, v1, a1;
      detail::Evaluation ev1;
      int used = 0;
      if (!advance(t, h, u1, v1, a1, ev1, used)) {
        if (cuts >= options.max_cuts) {
          std::ostringstream os;
          os << "Newton's method did not converge at t = " << t + h << " s after " << cuts
             << " halving(s) of the step; the last converged state is at t = " << t << " s";
          result.completed = false;
          result.termination = os.str();
          result.warnings.push_back(os.str());
          log::warn(os.str());
          break;
        }
        h *= 0.5;
        ++cuts;
        continue;
      }
      iterations += used;
      yielding = std::max(yielding, ev1.yielding_points);
      system.commit(ev1);
      Vector r1;
      const Vector force1 = body_force(ev1, v1, a1, r1);
      if (damped) {
        const Vector vs = v + v1;
        damping_energy += 0.25 * h * vs.dot(c * vs);
      }
      work += 0.5 * (u1 - u).dot(force + force1);
      u = std::move(u1);
      v = std::move(v1);
      a = std::move(a1);
      r = std::move(r1);
      force = force1;
      residual = ev1.residual;
      strain_energy = ev1.energy;
      t += h;
      if (h < dt) h = std::min(2.0 * h, dt);
    }
    if (!result.completed) break;
    t = target;  // the sub-steps summed to it within round-off
    record(step, target, iterations, cuts, yielding);
    if (options.snapshot_every > 0 && (step % options.snapshot_every == 0 || step == steps)) {
      snapshot(target);
    }
  }

  result.displacement = u;
  result.velocity = v;
  result.acceleration = a;
  result.reactions = r;
  result.linear_solver = factor.name();
  result.max_plastic_strain = system.max_plastic_strain();
  const TransientStep& last = result.steps.back();
  result.numerical_dissipation =
      e0 + last.external_work - last.kinetic_energy - last.strain_energy - last.damping_energy;
  result.energy_balance_error = energy_scale > 0.0 ? worst_balance / energy_scale : 0.0;
  if (options.snapshot_every == 0) snapshot(last.time);
  std::ostringstream os;
  os << "non-linear transient '" << result.load_case_name << "': "
     << static_cast<int>(result.steps.size()) - 1 << " of " << steps << " step(s) of " << dt
     << " s, " << result.total_iterations << " Newton iteration(s)"
     << (result.plastic ? ", largest plastic strain " + std::to_string(result.max_plastic_strain)
                        : std::string());
  log::info(os.str());
  return result;
}

}  // namespace

// ---------------------------------------------------------------------------
// Amplitude
// ---------------------------------------------------------------------------
Scalar Amplitude::value(Scalar t) const {
  switch (kind) {
    case Kind::Step:
      return scale;
    case Kind::Harmonic:
      return scale * std::sin(2.0 * kPi * frequency * t + phase);
    case Kind::Table: {
      if (t <= times.front()) return scale * values.front();
      if (t >= times.back()) return scale * values.back();
      const auto it = std::upper_bound(times.begin(), times.end(), t);
      const std::size_t i = static_cast<std::size_t>(it - times.begin());
      const Scalar s = (t - times[i - 1]) / (times[i] - times[i - 1]);
      return scale * (values[i - 1] + s * (values[i] - values[i - 1]));
    }
  }
  return 0.0;
}

Scalar Amplitude::rate(Scalar t) const {
  switch (kind) {
    case Kind::Step:
      return 0.0;
    case Kind::Harmonic:
      return scale * 2.0 * kPi * frequency * std::cos(2.0 * kPi * frequency * t + phase);
    case Kind::Table: {
      if (t < times.front() || t >= times.back()) return 0.0;
      const auto it = std::upper_bound(times.begin(), times.end(), t);
      const std::size_t i = static_cast<std::size_t>(it - times.begin());
      return scale * (values[i] - values[i - 1]) / (times[i] - times[i - 1]);
    }
  }
  return 0.0;
}

Scalar Amplitude::second_rate(Scalar t) const {
  if (kind != Kind::Harmonic) return 0.0;
  const Scalar w = 2.0 * kPi * frequency;
  return -scale * w * w * std::sin(w * t + phase);
}

void Amplitude::validate() const {
  if (!std::isfinite(scale)) throw ConfigError("amplitude: the scale must be finite");
  if (kind == Kind::Harmonic) {
    if (!(frequency > 0.0) || !std::isfinite(frequency)) {
      throw ConfigError("harmonic amplitude: the frequency must be positive");
    }
    if (!std::isfinite(phase)) throw ConfigError("harmonic amplitude: the phase must be finite");
  }
  if (kind == Kind::Table) {
    if (times.size() < 2 || times.size() != values.size()) {
      throw ConfigError("table amplitude: give at least two times and one value per time");
    }
    for (std::size_t i = 0; i < times.size(); ++i) {
      if (!std::isfinite(times[i]) || !std::isfinite(values[i])) {
        throw ConfigError("table amplitude: times and values must be finite");
      }
      if (i > 0 && !(times[i] > times[i - 1])) {
        throw ConfigError("table amplitude: the times must increase strictly");
      }
    }
  }
}

std::string to_string(Amplitude::Kind kind) {
  switch (kind) {
    case Amplitude::Kind::Step: return "step";
    case Amplitude::Kind::Table: return "table";
    case Amplitude::Kind::Harmonic: return "harmonic";
  }
  return "";
}

Amplitude::Kind parse_amplitude_kind(const std::string& text) {
  if (text == "step") return Amplitude::Kind::Step;
  if (text == "table") return Amplitude::Kind::Table;
  if (text == "harmonic") return Amplitude::Kind::Harmonic;
  throw ConfigError("unknown amplitude '" + text + "'; expected \"step\", \"table\" or \"harmonic\"");
}

std::string to_string(DynamicMonitor::Quantity quantity) {
  switch (quantity) {
    case DynamicMonitor::Quantity::Displacement: return "displacement";
    case DynamicMonitor::Quantity::Velocity: return "velocity";
    case DynamicMonitor::Quantity::Acceleration: return "acceleration";
    case DynamicMonitor::Quantity::Reaction: return "reaction";
  }
  return "";
}

DynamicMonitor::Quantity parse_dynamic_quantity(const std::string& text) {
  if (text == "displacement") return DynamicMonitor::Quantity::Displacement;
  if (text == "velocity") return DynamicMonitor::Quantity::Velocity;
  if (text == "acceleration") return DynamicMonitor::Quantity::Acceleration;
  if (text == "reaction") return DynamicMonitor::Quantity::Reaction;
  throw ConfigError("unknown monitor quantity '" + text +
                    "'; expected \"displacement\", \"velocity\", \"acceleration\" or "
                    "\"reaction\"");
}

Vector harmonic_peak_displacements(const FemModel& model, const ComplexVector& u) {
  const int dim = model.dim();
  Vector out(model.mesh().num_nodes());
  for (Index node = 0; node < model.mesh().num_nodes(); ++node) {
    // x(t) = a cos(omega t) - b sin(omega t) with U = a + i b: an ellipse
    // whose semi-major axis squared is (|a|^2 + |b|^2)/2 plus
    // sqrt(((|a|^2 - |b|^2)/2)^2 + (a.b)^2).
    Scalar aa = 0.0;
    Scalar bb = 0.0;
    Scalar ab = 0.0;
    for (int k = 0; k < dim; ++k) {
      const ComplexScalar z = u(model.dofs().dof(node, k));
      aa += z.real() * z.real();
      bb += z.imag() * z.imag();
      ab += z.real() * z.imag();
    }
    out(node) = std::sqrt(0.5 * (aa + bb) + std::hypot(0.5 * (aa - bb), ab));
  }
  return out;
}

std::string to_string(TransientOptions::Start start) {
  return start == TransientOptions::Start::Rest ? "rest" : "static";
}

TransientOptions::Start parse_transient_start(const std::string& text) {
  if (text == "rest") return TransientOptions::Start::Rest;
  if (text == "static") return TransientOptions::Start::Static;
  throw ConfigError("unknown initial state '" + text + "'; expected \"rest\" or \"static\"");
}

HhtParameters HhtParameters::from_alpha(Scalar alpha) {
  if (!(alpha >= -1.0 / 3.0 - 1.0e-15 && alpha <= 0.0)) {
    std::ostringstream os;
    os << "the HHT-alpha parameter must lie in [-1/3, 0] (got " << alpha
       << "); 0 is the trapezoidal rule, -0.05 to -0.1 a common light numerical damping";
    throw ConfigError(os.str());
  }
  HhtParameters p;
  p.alpha = alpha;
  p.beta = 0.25 * (1.0 - alpha) * (1.0 - alpha);
  p.gamma = 0.5 - alpha;
  return p;
}

// ---------------------------------------------------------------------------
// Transient
// ---------------------------------------------------------------------------
TransientResult solve_transient(const FemModel& model, const Assembler& assembler,
                                std::size_t load_case, const TransientOptions& options) {
  check_load_case(model, load_case);
  const Scalar dt = options.time_step;
  if (!(dt > 0.0) || !std::isfinite(dt)) throw ConfigError("transient: time_step must be positive");
  if (!(options.end_time > 0.0) || !std::isfinite(options.end_time)) {
    throw ConfigError("transient: end_time must be positive");
  }
  const Scalar ratio = options.end_time / dt;
  const Scalar rounded = std::round(ratio);
  if (rounded < 1.0 || std::abs(ratio - rounded) > 1.0e-9 * std::max(1.0, ratio)) {
    std::ostringstream os;
    os << "transient: end_time " << options.end_time << " s is not a whole number of steps of "
       << dt << " s (" << ratio << "); adjust either";
    throw ConfigError(os.str());
  }
  const int steps = static_cast<int>(rounded);
  if (options.mass_damping < 0.0 || options.stiffness_damping < 0.0) {
    throw ConfigError("transient: Rayleigh damping coefficients must be non-negative");
  }
  if (options.snapshot_every < 0) throw ConfigError("transient: snapshot_every must be >= 0");
  options.amplitude.validate();
  const HhtParameters p = HhtParameters::from_alpha(options.alpha);

  const DofManager& dofs = model.dofs();
  const std::vector<Index>& free = dofs.free_dofs();
  const std::vector<Index>& fixed = dofs.constrained_dofs();
  const Index n = dofs.num_dofs();

  TransientResult result;
  result.load_case_name = model.load_case_specs()[load_case].name;
  result.parameters = p;
  result.time_step = dt;
  result.num_steps = steps;
  const std::vector<std::vector<Index>> nodes = monitor_nodes(model, options.monitors);
  for (const DynamicMonitor& m : options.monitors) {
    result.monitor_names.push_back(m.name);
    result.monitor_units.push_back(monitor_unit(m.quantity));
  }

  const SparseMatrix k = assembler.assemble_stiffness();
  const SparseMatrix m = assembler.assemble_mass(options.mass_type);
  result.total_mass = model_mass(model, m);
  if (!(result.total_mass > 0.0)) {
    throw ModelError("the model has no mass; give the material a positive density for a "
                     "dynamic analysis");
  }
  if (options.nonlinear) {
    return nonlinear_transient(model, assembler, load_case, options, steps, p, std::move(result),
                               nodes, m, k);
  }
  const bool damped = options.mass_damping > 0.0 || options.stiffness_damping > 0.0;
  SparseMatrix c(n, n);
  if (damped) c = options.mass_damping * m + options.stiffness_damping * k;

  const Vector f = model.load_vectors()[load_case];
  const Vector g = dofs.prescribed_vector();
  const Vector g_p = restrict(g, fixed);

  // Newmark coefficients: a1 = c0 (u1 - u0) - c1 v0 - c2 a0,
  // v1 = c3 (u1 - u0) + c4 v0 + c5 a0.
  const Scalar c0 = 1.0 / (p.beta * dt * dt);
  const Scalar c1 = 1.0 / (p.beta * dt);
  const Scalar c2 = 1.0 / (2.0 * p.beta) - 1.0;
  const Scalar c3 = p.gamma / (p.beta * dt);
  const Scalar c4 = 1.0 - p.gamma / p.beta;
  const Scalar c5 = dt * (1.0 - p.gamma / (2.0 * p.beta));
  const Scalar w = 1.0 + p.alpha;

  SparseMatrix k_eff = c0 * m + w * k;
  if (damped) k_eff += (w * c3) * c;
  const SparseMatrix k_eff_ff = assembler.reduce_free_free(k_eff);
  const SparseMatrix k_eff_fp = assembler.reduce_free_prescribed(k_eff);

  const DofLayout layout = layout_of(model);
  std::unique_ptr<LinearSolver> solver = make_linear_solver(options.linear);
  solver->set_layout(layout);
  solver->factorize(k_eff_ff);
  result.linear_solver = solver->name();

  // The initial state.
  const Amplitude& amp = options.amplitude;
  Scalar a_now = amp.value(0.0);
  Vector u = Vector::Zero(n);
  Vector v = Vector::Zero(n);
  Vector a = Vector::Zero(n);
  scatter(u, fixed, a_now * g_p);
  scatter(v, fixed, amp.rate(0.0) * g_p);
  scatter(a, fixed, amp.second_rate(0.0) * g_p);
  if (options.start == TransientOptions::Start::Static) {
    const SparseMatrix k_ff = assembler.reduce_free_free(k);
    const SparseMatrix k_fp = assembler.reduce_free_prescribed(k);
    std::unique_ptr<LinearSolver> statics = make_linear_solver(options.linear);
    statics->set_layout(layout);
    statics->factorize(k_ff);
    Vector rhs = a_now * restrict(f, free);
    if (fixed.size() > 0) rhs -= k_fp * restrict(u, fixed);
    scatter(u, free, statics->solve(rhs));
  }
  {
    // M_ff a_f = A(0) f_f - (M a_p + C v + K u)_f.
    Vector known = k * u;
    if (damped) known += c * v;
    Vector m_a = m * a;  // the prescribed accelerations only, so far
    const Vector rhs = restrict(a_now * f - known - m_a, free);
    const SparseMatrix m_ff = assembler.reduce_free_free(m);
    std::unique_ptr<LinearSolver> mass = make_linear_solver(options.linear);
    mass->set_layout(layout);
    mass->factorize(m_ff);
    scatter(a, free, mass->solve(rhs));
  }

  // The reaction at the prescribed DOFs, M a + C v + K u - A f, and the full
  // force F = A f + r whose trapezoidal work the energy balance counts.
  const auto reactions_of = [&](const Vector& uu, const Vector& vv, const Vector& aa,
                                Scalar amplitude) {
    Vector internal = m * aa + k * uu;
    if (damped) internal += c * vv;
    Vector r = Vector::Zero(n);
    for (Index d : fixed) r(d) = internal(d) - amplitude * f(d);
    return r;
  };

  Vector r = reactions_of(u, v, a, a_now);
  Vector force = a_now * f + r;
  const Scalar e0 = 0.5 * v.dot(m * v) + 0.5 * u.dot(k * u);
  Scalar damping_energy = 0.0;
  Scalar work = 0.0;
  Scalar energy_scale = std::abs(e0);
  Scalar worst_balance = 0.0;

  const auto record = [&](int index, Scalar t) {
    TransientStep s;
    s.index = index;
    s.time = t;
    for (std::size_t i = 0; i < options.monitors.size(); ++i) {
      s.monitors.push_back(monitor_value(model, options.monitors[i], nodes[i], u, v, a, r));
    }
    s.max_displacement = max_nodal_magnitude(model, u);
    s.kinetic_energy = 0.5 * v.dot(m * v);
    s.strain_energy = 0.5 * u.dot(k * u);
    s.damping_energy = damping_energy;
    s.external_work = work;
    energy_scale = std::max({energy_scale, s.kinetic_energy, s.strain_energy, damping_energy,
                             std::abs(work)});
    const Scalar balance = e0 + work - s.kinetic_energy - s.strain_energy - damping_energy;
    worst_balance = std::max(worst_balance, std::abs(balance));
    result.steps.push_back(std::move(s));
  };
  const auto snapshot = [&](Scalar t) {
    TransientSnapshot snap;
    snap.time = t;
    snap.displacement = u;
    snap.velocity = v;
    result.snapshots.push_back(std::move(snap));
  };
  record(0, 0.0);
  if (options.snapshot_every > 0) snapshot(0.0);

  // A harmonic amplitude needs its period resolved: the trapezoidal rule
  // lengthens a period T by (omega dt)^2 / 12 of itself.
  if (amp.kind == Amplitude::Kind::Harmonic) {
    const Scalar per_period = 1.0 / (amp.frequency * dt);
    if (per_period < 20.0) {
      std::ostringstream os;
      os << "the harmonic amplitude's period is resolved by only " << per_period
         << " steps; the trapezoidal rule lengthens it by (omega dt)^2 / 12 = "
         << std::pow(2.0 * kPi / per_period, 2) / 12.0 << " - take time_step below T / 20";
      result.warnings.push_back(os.str());
      log::warn(os.str());
    }
  }

  for (int step = 1; step <= steps; ++step) {
    const Scalar t = step * dt;
    const Scalar a_next = amp.value(t);
    const Vector a_hat = -c0 * u - c1 * v - c2 * a;
    const Vector v_hat = -c3 * u + c4 * v + c5 * a;
    Vector rhs = (w * a_next - p.alpha * a_now) * f - m * a_hat;
    if (damped) rhs += -w * (c * v_hat) + p.alpha * (c * v);
    if (p.alpha != 0.0) rhs += p.alpha * (k * u);
    Vector u_next = Vector::Zero(n);
    scatter(u_next, fixed, a_next * g_p);
    Vector rhs_f = restrict(rhs, free);
    if (fixed.size() > 0) rhs_f -= k_eff_fp * restrict(u_next, fixed);
    const Vector u_f = solver->solve(rhs_f);
    if (!u_f.allFinite()) {
      throw SolverError("transient: the effective stiffness solve returned a non-finite "
                        "displacement at t = " + std::to_string(t) + " s");
    }
    scatter(u_next, free, u_f);
    const Vector a_next_vec = c0 * u_next + a_hat;
    const Vector v_next = c3 * u_next + v_hat;
    const Vector r_next = reactions_of(u_next, v_next, a_next_vec, a_next);
    const Vector force_next = a_next * f + r_next;
    if (damped) {
      const Vector vs = v + v_next;
      damping_energy += 0.25 * dt * vs.dot(c * vs);
    }
    work += 0.5 * (u_next - u).dot(force + force_next);
    u = u_next;
    v = v_next;
    a = a_next_vec;
    r = r_next;
    force = force_next;
    a_now = a_next;
    record(step, t);
    if (options.snapshot_every > 0 && (step % options.snapshot_every == 0 || step == steps)) {
      snapshot(t);
    }
  }

  result.displacement = u;
  result.velocity = v;
  result.acceleration = a;
  result.reactions = r;
  const TransientStep& last = result.steps.back();
  result.numerical_dissipation =
      e0 + last.external_work - last.kinetic_energy - last.strain_energy - last.damping_energy;
  result.energy_balance_error = energy_scale > 0.0 ? worst_balance / energy_scale : 0.0;
  if (options.snapshot_every == 0) snapshot(result.steps.back().time);
  std::ostringstream os;
  os << "transient '" << result.load_case_name << "': " << steps << " step(s) of " << dt
     << " s (HHT alpha = " << p.alpha << "), energy balance " << result.energy_balance_error;
  log::info(os.str());
  return result;
}

// ---------------------------------------------------------------------------
// Frequency response
// ---------------------------------------------------------------------------
FrequencyResponseResult solve_frequency_response(const FemModel& model,
                                                 const Assembler& assembler,
                                                 std::size_t load_case,
                                                 const FrequencyResponseOptions& options) {
  check_load_case(model, load_case);
  if (options.frequencies.empty()) throw ConfigError("frequency response: give frequencies");
  for (Scalar freq : options.frequencies) {
    if (!(freq >= 0.0) || !std::isfinite(freq)) {
      throw ConfigError("frequency response: frequencies must be finite and non-negative");
    }
  }
  if (options.mass_damping < 0.0 || options.stiffness_damping < 0.0 ||
      options.structural_damping < 0.0) {
    throw ConfigError("frequency response: damping coefficients must be non-negative");
  }

  const DofManager& dofs = model.dofs();
  const std::vector<Index>& free = dofs.free_dofs();
  const std::vector<Index>& fixed = dofs.constrained_dofs();
  const Index n = dofs.num_dofs();

  FrequencyResponseResult result;
  result.load_case_name = model.load_case_specs()[load_case].name;
  const std::vector<std::vector<Index>> nodes = monitor_nodes(model, options.monitors);
  for (const DynamicMonitor& m : options.monitors) {
    result.monitor_names.push_back(m.name);
    result.monitor_units.push_back(monitor_unit(m.quantity));
  }

  const SparseMatrix k = assembler.assemble_stiffness();
  const SparseMatrix m = assembler.assemble_mass(options.mass_type);
  const ComplexSparse k_ff = assembler.reduce_free_free(k).cast<ComplexScalar>();
  const ComplexSparse m_ff = assembler.reduce_free_free(m).cast<ComplexScalar>();
  const ComplexSparse k_fp = assembler.reduce_free_prescribed(k).cast<ComplexScalar>();
  const ComplexSparse m_fp = assembler.reduce_free_prescribed(m).cast<ComplexScalar>();
  const ComplexVector f = model.load_vectors()[load_case].cast<ComplexScalar>();
  const ComplexVector g = dofs.prescribed_vector().cast<ComplexScalar>();
  ComplexVector g_p(static_cast<Eigen::Index>(fixed.size()));
  for (std::size_t i = 0; i < fixed.size(); ++i) g_p(static_cast<Eigen::Index>(i)) = g(fixed[i]);
  ComplexVector f_f(static_cast<Eigen::Index>(free.size()));
  for (std::size_t i = 0; i < free.size(); ++i) f_f(static_cast<Eigen::Index>(i)) = f(free[i]);

  // The snapshot frequencies, each the nearest solved one.
  std::vector<char> keep(options.frequencies.size(), 0);
  for (Scalar want : options.snapshot_frequencies) {
    std::size_t best = 0;
    for (std::size_t i = 1; i < options.frequencies.size(); ++i) {
      if (std::abs(options.frequencies[i] - want) < std::abs(options.frequencies[best] - want)) {
        best = i;
      }
    }
    keep[best] = 1;
  }

  const ComplexScalar i_unit(0.0, 1.0);
  // The static response to the same loads, the reference of the dynamic
  // amplification: undamped at a natural frequency the dynamic stiffness is
  // singular, and a backward-stable solve still returns a huge answer with
  // a small backward error - one amplified a million times over the static
  // response is round-off.
  Scalar static_peak = 0.0;
  {
    const SparseMatrix k_real = assembler.reduce_free_free(k);
    Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<StorageIndex>> statics;
    statics.compute(k_real);
    if (statics.info() == Eigen::Success) {
      Vector rhs = restrict(model.load_vectors()[load_case], free);
      if (fixed.size() > 0) {
        rhs -= assembler.reduce_free_prescribed(k) * restrict(dofs.prescribed_vector(), fixed);
      }
      const Vector u_static = statics.solve(rhs);
      if (u_static.allFinite()) static_peak = u_static.cwiseAbs().maxCoeff();
    }
  }
  Eigen::SparseLU<ComplexSparse, Eigen::COLAMDOrdering<StorageIndex>> lu;
  bool analysed = false;
  for (std::size_t j = 0; j < options.frequencies.size(); ++j) {
    const Scalar freq = options.frequencies[j];
    const Scalar omega = 2.0 * kPi * freq;
    // K (1 + i eta) - omega^2 M + i omega (a M + b K).
    const ComplexScalar kc = 1.0 + i_unit * (options.structural_damping +
                                             omega * options.stiffness_damping);
    const ComplexScalar mc = -omega * omega + i_unit * omega * options.mass_damping;
    const ComplexSparse a_ff = kc * k_ff + mc * m_ff;
    ComplexVector rhs = f_f;
    if (fixed.size() > 0) rhs -= (kc * k_fp + mc * m_fp) * g_p;
    if (!analysed) {
      lu.analyzePattern(a_ff);
      analysed = true;
    }
    lu.factorize(a_ff);
    std::ostringstream where;
    where << freq << " Hz";
    if (lu.info() != Eigen::Success) {
      throw SolverError("frequency response: the dynamic stiffness is singular at " +
                        where.str() +
                        " - an undamped natural frequency (or a mechanism); give the model "
                        "damping (structural_damping or Rayleigh) or move the frequency");
    }
    const ComplexVector u_f = lu.solve(rhs);
    if (!u_f.allFinite()) {
      throw SolverError("frequency response: a non-finite response at " + where.str() +
                        " - the dynamic stiffness is singular there; add damping");
    }
    // Backward error of the solve: |r| over || |A| |U| + |b| ||.
    const ComplexVector residual = a_ff * u_f - rhs;
    Eigen::VectorXd gross = rhs.cwiseAbs();
    for (Eigen::Index col = 0; col < a_ff.outerSize(); ++col) {
      const Scalar uc = std::abs(u_f(col));
      for (ComplexSparse::InnerIterator it(a_ff, col); it; ++it) {
        gross(it.row()) += std::abs(it.value()) * uc;
      }
    }
    const Scalar backward = residual.norm() / std::max(gross.norm(), 1.0e-300);
    if (backward > 1.0e-10) {
      std::ostringstream os;
      os << "frequency response: the solve at " << where.str() << " has a backward error of "
         << backward << " - the dynamic stiffness is nearly singular (a natural frequency of a "
         << "lightly damped model); add damping or move the frequency";
      throw SolverError(os.str());
    }
    ComplexVector u = ComplexVector::Zero(n);
    for (std::size_t i = 0; i < free.size(); ++i) u(free[i]) = u_f(static_cast<Eigen::Index>(i));
    for (std::size_t i = 0; i < fixed.size(); ++i) u(fixed[i]) = g_p(static_cast<Eigen::Index>(i));
    // Reactions at the prescribed DOFs: (K (1 + i eta) - w^2 M + i w C) U - f.
    const ComplexVector ku = k.cast<ComplexScalar>() * u;
    const ComplexVector mu = m.cast<ComplexScalar>() * u;
    ComplexVector r = ComplexVector::Zero(n);
    for (Index d : fixed) r(d) = kc * ku(d) + mc * mu(d) - f(d);
    const ComplexVector v = (i_unit * omega) * u;
    const ComplexVector acc = (-omega * omega) * u;

    const Scalar peak = u_f.cwiseAbs().maxCoeff();
    if (static_peak > 0.0 && peak > 1.0e6 * static_peak) {
      std::ostringstream os;
      os << "the response at " << freq << " Hz is " << peak / static_peak
         << " times the static one: the frequency lies within about "
         << static_peak / peak << " of an undamped natural frequency, where the answer is "
         << "round-off; give the model damping (structural_damping or Rayleigh) or move the "
         << "frequency";
      result.warnings.push_back(os.str());
      log::warn(os.str());
    }
    FrequencyPoint point;
    point.frequency = freq;
    for (std::size_t i = 0; i < options.monitors.size(); ++i) {
      point.monitors.push_back(monitor_value(model, options.monitors[i], nodes[i], u, v, acc, r));
    }
    point.max_displacement = harmonic_peak_displacements(model, u).maxCoeff();
    result.points.push_back(std::move(point));
    if (keep[j]) {
      FrequencySnapshot snap;
      snap.frequency = freq;
      snap.displacement = u;
      result.snapshots.push_back(std::move(snap));
    }
  }
  std::ostringstream os;
  os << "frequency response '" << result.load_case_name << "': " << options.frequencies.size()
     << " frequenc" << (options.frequencies.size() == 1 ? "y" : "ies");
  log::info(os.str());
  return result;
}

}  // namespace sparlab
