#include "sparlab/fem/NonlinearStatic.hpp"

#include "NonlinearSystem.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"
#include "sparlab/elements/FaceGeometry.hpp"
#include "sparlab/fem/Loads.hpp"
#include "sparlab/fem/StressRecovery.hpp"
#include "sparlab/fem/TotalLagrangian.hpp"

#include <Eigen/Geometry>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <sstream>

namespace sparlab {
namespace {

using namespace detail;

/// Small-strain kinematics neglects the quadratic part of the Green strain,
/// H^T H / 2; a run warns once it exceeds this fraction of the largest strain.
constexpr Scalar kQuadraticWarning = 0.1;
/// The small-strain theory and the elastoplastic law assume small strains; a
/// run warns beyond this strain.
constexpr Scalar kStrainWarning = 0.05;

}  // namespace

NonlinearState evaluate_nonlinear_state(const FemModel& model, const Assembler& assembler,
                                        std::size_t load_case, const NonlinearOptions& options,
                                        const Vector& u, Scalar lambda) {
  const NonlinearSystem system(model, assembler, load_case, options);
  Evaluation ev = system.evaluate(u, lambda, true);
  NonlinearState out;
  out.residual = std::move(ev.residual);
  out.external = std::move(ev.external);
  out.load_rate = std::move(ev.load_rate);
  out.tangent = std::move(ev.tangent);
  out.energy = ev.energy;
  return out;
}

std::string to_string(MeanDilatation mode) {
  switch (mode) {
    case MeanDilatation::Auto: return "auto";
    case MeanDilatation::All: return "all";
    case MeanDilatation::None: return "none";
  }
  return "auto";
}

MeanDilatation parse_mean_dilatation(const std::string& text) {
  if (text == "auto") return MeanDilatation::Auto;
  if (text == "all") return MeanDilatation::All;
  if (text == "none") return MeanDilatation::None;
  throw ConfigError("unknown mean_dilatation '" + text +
                    "'; expected \"auto\" (Q4 and Hex8), \"all\" or \"none\"");
}

std::string to_string(NonlinearOptions::Method method) {
  return method == NonlinearOptions::Method::ArcLength ? "arc_length" : "load_control";
}

NonlinearOptions::Method parse_nonlinear_method(const std::string& text) {
  if (text == "load_control") return NonlinearOptions::Method::LoadControl;
  if (text == "arc_length") return NonlinearOptions::Method::ArcLength;
  throw ConfigError("unknown non-linear method '" + text +
                    "'; expected \"load_control\" or \"arc_length\"");
}

std::string to_string(NonlinearMonitor::Quantity quantity) {
  return quantity == NonlinearMonitor::Quantity::Reaction ? "reaction" : "displacement";
}

NonlinearMonitor::Quantity parse_monitor_quantity(const std::string& text) {
  if (text == "displacement") return NonlinearMonitor::Quantity::Displacement;
  if (text == "reaction") return NonlinearMonitor::Quantity::Reaction;
  throw ConfigError("unknown monitor quantity '" + text +
                    "'; expected \"displacement\" or \"reaction\"");
}

NonlinearStaticAnalysis::NonlinearStaticAnalysis(const FemModel& model,
                                                 const Assembler& assembler,
                                                 NonlinearOptions options)
    : model_(model), assembler_(assembler), options_(std::move(options)) {
  if (!model.finalized()) throw ModelError("the model must be finalised before a solve");
  if (model.dofs_per_node() != model.dim()) {
    throw ConfigError("the non-linear analysis is written for continuum elements");
  }
  if (options_.steps < 1 || options_.max_iterations < 1 || options_.max_steps < 1) {
    throw ConfigError("the non-linear analysis needs steps, max_steps and max_iterations >= 1");
  }
  if (!(options_.residual_tolerance > 0.0) || !(options_.displacement_tolerance > 0.0)) {
    throw ConfigError("the non-linear tolerances must be positive");
  }
  // Refuse what the laws refuse before any step, so that it is reported as
  // what it is and not as a step that failed to converge.
  if (options_.law == HyperelasticModel::NeoHookean &&
      options_.kinematics == Kinematics::SmallStrain) {
    throw ConfigError("the small-strain kinematics is linear elasticity; the neo-Hookean law "
                      "needs \"finite\" kinematics");
  }
  if (options_.law == HyperelasticModel::NeoHookean) {
    for (Index e = 0; e < model.mesh().num_elements(); ++e) {
      if (!model.material_of(e).plasticity().enabled()) continue;
      throw ConfigError("material '" + model.material_of(e).name() +
                        "' is elastoplastic, which takes the Saint Venant-Kirchhoff form; the "
                        "neo-Hookean law cannot be combined with plasticity");
    }
  }
  if (options_.law == HyperelasticModel::NeoHookean &&
      model.stress_state() == StressState::PlaneStress) {
    throw ConfigError("the neo-Hookean law needs plane strain or a solid mesh; in plane "
                      "stress use \"saint_venant_kirchhoff\"");
  }
  if (options_.contact.enabled) {
    if (options_.kinematics != Kinematics::SmallStrain) {
      throw ConfigError("contact is formulated for small displacements - the contact geometry "
                        "of the reference configuration and a gap linear in the displacement - "
                        "so it needs \"small_strain\" kinematics in the non-linear analysis");
    }
    if (options_.method == NonlinearOptions::Method::ArcLength) {
      throw ConfigError("contact is solved under load control; the arc-length method is not "
                        "available with it");
    }
  }
  if (!options_.load_path.empty()) {
    if (options_.method == NonlinearOptions::Method::ArcLength) {
      throw ConfigError("'load_path' is followed by load control; the arc-length method "
                        "follows the equilibrium path and cannot unload");
    }
    if (!options_.load_factors.empty()) {
      throw ConfigError("give either 'load_path' or 'load_factors', not both");
    }
    Scalar previous = 0.0;
    for (Scalar f : options_.load_path) {
      if (!std::isfinite(f) || f == previous) {
        throw ConfigError("'load_path' needs finite load factors, each different from the one "
                          "before it (the path starts at 0)");
      }
      previous = f;
    }
  }
  if (!options_.load_factors.empty()) {
    if (options_.method == NonlinearOptions::Method::ArcLength) {
      throw ConfigError("'load_factors' fixes the load levels of load control; the "
                        "arc-length method chooses its own");
    }
    Scalar previous = 0.0;
    for (Scalar f : options_.load_factors) {
      if (!(f > previous) || !(f < 1.0)) {
        throw ConfigError("'load_factors' must increase strictly within (0, 1)");
      }
      previous = f;
    }
  }
}

NonlinearResult NonlinearStaticAnalysis::solve(std::size_t load_case) {
  const Mesh& mesh = model_.mesh();
  const int dim = mesh.dim();
  const Index n = model_.dofs().num_dofs();
  const std::vector<Index>& free_dofs = model_.dofs().free_dofs();
  const std::vector<Index>& fixed = model_.dofs().constrained_dofs();
  const LoadCaseSpec& spec = model_.load_case_specs().at(load_case);
  if (options_.law == HyperelasticModel::NeoHookean &&
      model_.load_case_data(load_case).temperature.size() > 0) {
    for (Index e = 0; e < mesh.num_elements(); ++e) {
      if (model_.material_of(e).thermal_expansion() == 0.0) continue;
      throw ConfigError("load case '" + spec.name +
                        "' has a temperature change, whose thermal strain the neo-Hookean "
                        "law does not model; use \"saint_venant_kirchhoff\"");
    }
  }
  NonlinearSystem system(model_, assembler_, load_case, options_);
  std::unique_ptr<ContactProblem> contact;
  if (options_.contact.enabled) contact = std::make_unique<ContactProblem>(model_, options_.contact);
  // The contact status of the last converged state (all open at the start).
  std::vector<ContactStatus> contact_status;
  if (contact) contact_status.assign(contact->nodes().size(), ContactStatus::Open);
  std::string contact_failure;  // why the last contact iteration gave up
  const bool arc = options_.method == NonlinearOptions::Method::ArcLength;
  const bool small = options_.kinematics == Kinematics::SmallStrain;
  const bool logarithmic = options_.kinematics == Kinematics::FiniteLogarithmic;
  // The jump test compares the converged state with the elastic tangent's
  // prediction: meaningful with finite kinematics and elastic materials only
  // (see the file comment).
  const bool jump_test = !small && !system.plastic();

  Vector prescribed(static_cast<Eigen::Index>(fixed.size()));
  for (std::size_t i = 0; i < fixed.size(); ++i) {
    prescribed(static_cast<Eigen::Index>(i)) = model_.dofs().prescribed_value(fixed[i]);
  }
  if (arc && prescribed.size() > 0 && prescribed.cwiseAbs().maxCoeff() > 0.0) {
    throw ConfigError("the arc-length method needs homogeneous prescribed displacements; "
                      "drive a prescribed displacement with \"load_control\" instead");
  }

  NonlinearResult result;
  result.load_case_name = spec.name;
  result.method = to_string(options_.method);
  result.law = logarithmic ? "hencky" : to_string(options_.law);
  result.kinematics = to_string(options_.kinematics);
  result.plastic = system.plastic();
  for (Index e = 0; e < mesh.num_elements(); ++e) {
    if (system.averaged(e)) result.mean_dilatation = true;
  }
  std::vector<std::vector<Index>> monitor_nodes;
  for (const NonlinearMonitor& m : options_.monitors) {
    if (m.component < 0 || m.component >= dim) {
      throw ConfigError("monitor '" + m.name + "' asks for a component the model lacks");
    }
    std::vector<Index> nodes = m.region.select_nodes(mesh);
    if (nodes.empty()) throw ConfigError("monitor '" + m.name + "' selects no node");
    monitor_nodes.push_back(std::move(nodes));
    result.monitor_names.push_back(m.name);
    result.monitor_units.push_back(m.quantity == NonlinearMonitor::Quantity::Reaction ? "N"
                                                                                     : "m");
  }

  // A displacement monitor is the mean over its nodes; a reaction monitor
  // the sum of the full residual f_int - f_ext (the reaction at a prescribed
  // DOF, round-off at a free one).
  const auto monitors_of = [&](const Vector& u, const Vector& residual) {
    std::vector<Scalar> values;
    for (std::size_t i = 0; i < monitor_nodes.size(); ++i) {
      const NonlinearMonitor& m = options_.monitors[i];
      const bool reaction = m.quantity == NonlinearMonitor::Quantity::Reaction;
      Scalar sum = 0.0;
      for (Index node : monitor_nodes[i]) {
        sum += (reaction ? residual : u)(node * dim + m.component);
      }
      values.push_back(reaction ? sum : sum / static_cast<Scalar>(monitor_nodes[i].size()));
    }
    return values;
  };
  const auto max_magnitude = [&](const Vector& u) {
    Scalar top = 0.0;
    for (Index node = 0; node < mesh.num_nodes(); ++node) {
      top = std::max(top, u.segment(node * dim, dim).norm());
    }
    return top;
  };
  // Scale of the residual: the external load, the reactions under pure
  // displacement control, or the thermal forces of a temperature change.
  const auto residual_scale = [&](const Evaluation& ev) {
    const Scalar reactions = restrict(ev.residual, fixed).norm();
    return std::max({ev.external.norm(), reactions, ev.thermal, 1.0e-300});
  };

  // A correction below this is converged: the displacement tolerance
  // relative to the step's increment, no finer than the round-off of the
  // displacement itself.
  const auto correction_limit = [&](Scalar increment, const Vector& state) {
    return std::max(options_.displacement_tolerance * std::max(increment, 1.0e-300),
                    64.0 * std::numeric_limits<Scalar>::epsilon() * state.norm());
  };

  TangentFactor factor;
  Vector u = Vector::Zero(n);
  Scalar lambda = 0.0;
  const Scalar target = arc ? options_.target_load_factor : 1.0;
  if (!(target > 0.0)) throw ConfigError("the target load factor must be positive");
  const Scalar initial_step = target / static_cast<Scalar>(options_.steps);
  Scalar step = initial_step;
  Scalar arc_length = 0.0;
  Scalar first_arc = 0.0;
  Vector previous_increment;  // free DOFs, the last converged step's
  int cuts_in_a_row = 0;

  // Record the converged state (u, lambda) with its full residual; with
  // contact, the pair totals of the step (from the state it started at).
  std::vector<ContactPairResult> step_contact;
  const auto record = [&](int iterations, Scalar residual, int cuts, int pivots,
                          const Vector& full_residual, int yielding) {
    NonlinearStep s;
    for (const ContactPairResult& p : step_contact) {
      s.contact_active.push_back(p.active);
      s.contact_force.push_back(p.force);
    }
    s.index = static_cast<int>(result.steps.size()) + 1;
    s.load_factor = lambda;
    s.iterations = iterations;
    s.residual = residual;
    s.arc_length = arc_length;
    s.negative_pivots = pivots;
    s.cuts = cuts;
    s.monitors = monitors_of(u, full_residual);
    s.max_displacement = max_magnitude(u);
    if (system.plastic()) {
      s.yielding_points = yielding;
      s.max_plastic_strain = system.max_plastic_strain();
    }
    result.steps.push_back(s);
    result.total_iterations += iterations;
    std::ostringstream os;
    os << "non-linear step " << s.index << ": lambda = " << lambda << ", " << iterations
       << " iteration(s), residual " << residual;
    for (std::size_t i = 0; i < s.monitors.size(); ++i) {
      os << ", " << result.monitor_names[i] << " = " << s.monitors[i];
    }
    if (system.plastic()) {
      os << ", " << yielding << " point(s) yielding, largest plastic strain "
         << s.max_plastic_strain;
    }
    for (const ContactPairResult& p : step_contact) {
      os << ", contact '" << p.name << "' " << p.active << " node(s), force "
         << p.force.head(dim).transpose();
    }
    log::info(os.str());
  };

  // Line search along a Newton direction du from `state` at load factor lam,
  // on the energy: g(a) = du . R_f(state + a du) is the derivative of the
  // potential along du and vanishes where the energy is least. The full step
  // is kept when |g(1)| <= 0.8 |g(0)| (Crisfield's slack criterion);
  // otherwise regula falsi on g, bracketed by g(0) < 0 < g(1), finds a
  // shorter one. The norm of the residual would be the wrong measure: in a
  // slender structure it is dominated by the stiff axial forces that a large
  // rotation stirs up while the Newton step is still good, and backtracking on
  // it stalls Newton to a linear crawl. g weighs the residual by the step, so
  // it keeps the full step whenever the energy agrees. A direction that is
  // not one of descent (g(0) >= 0: a non-symmetric follower-pressure tangent,
  // or an indefinite one) takes the full step. A step into an inverted
  // element is halved until it is valid.
  const auto line_search = [&](const Vector& state, const Vector& du, const Vector& r,
                               Scalar lam) -> Scalar {
    const Scalar g0 = du.dot(r);
    if (!(g0 < 0.0)) return 1.0;
    const auto g = [&](Scalar a) -> Scalar {
      Vector probe = state;
      add_to(probe, free_dofs, du, a);
      try {
        return du.dot(restrict(system.evaluate(probe, lam, false).residual, free_dofs));
      } catch (const SolverError&) {
        return std::numeric_limits<Scalar>::quiet_NaN();
      }
    };
    constexpr Scalar slack = 0.8;
    constexpr Scalar shortest = 0.1;
    Scalar a_hi = 1.0;
    Scalar g_hi = g(a_hi);
    for (int k = 0; k < 6 && !std::isfinite(g_hi); ++k) {
      a_hi *= 0.5;
      g_hi = g(a_hi);
    }
    if (!std::isfinite(g_hi) || std::abs(g_hi) <= slack * std::abs(g0) || g_hi < 0.0) {
      return a_hi;
    }
    Scalar a_lo = 0.0;
    Scalar g_lo = g0;
    Scalar best = a_hi;  // the last step with a valid state
    for (int it = 0; it < 5; ++it) {
      const Scalar a =
          std::clamp(a_lo - g_lo * (a_hi - a_lo) / (g_hi - g_lo), shortest * a_hi, a_hi);
      const Scalar ga = g(a);
      if (!std::isfinite(ga)) {
        a_hi = a;
        g_hi = std::abs(g0);  // an inverted element: far beyond the minimum
        continue;
      }
      best = a;
      if (std::abs(ga) <= slack * std::abs(g0)) break;
      if (ga < 0.0) {
        a_lo = a;
        g_lo = ga;
      } else {
        a_hi = a;
        g_hi = ga;
      }
    }
    return best;
  };

  // Newton iterations from the converged state `state` (at lambda) under
  // load control to lambda_new, from the tangent predictor
  // K_T du = dl q - R (with the prescribed-displacement increment). Returns
  // true on convergence with `state` updated, `full_residual` its residual
  // (all DOFs) and `pivots` the negative pivots of the last tangent
  // factorised - at a state already within the residual tolerance, so the
  // inertia of the converged state (-1 when LU factorised it).
  //
  // Two kinds of step are abandoned even though Newton might converge, with
  // `failure` saying which, because past a limit point load control either
  // fails or - worse - lands silently on a distant branch of equilibria (a
  // snap-through, which is not quasi-static):
  //   * from a stable state (a positive definite tangent), an iteration that
  //     meets a tangent with a negative pivot: near a stable stretch of the
  //     path the iterates of a short enough step stay near it, and their
  //     tangents positive definite;
  //   * a converged state farther from the tangent predictor than the
  //     predicted increment itself: on a smooth path the corrector shrinks
  //     with the step, O(dl) relative to the predictor (0.01 to 0.8 on the
  //     verification problems), while a jump to another branch leaves it
  //     O(1) or larger (21 on the snapping arch) however short the step.
  // Halving the step tells an iterate that strayed from a limit point: the
  // first passes on a shorter step, the second never does.
  enum class StepFailure { None, Diverged, Unstable, Jumped };
  const auto newton_load_control = [&](Vector& state, Scalar lambda_new, int& iterations,
                                       Scalar& residual_out, int& pivots,
                                       Vector& full_residual, int& yielding,
                                       StepFailure& failure) -> bool {
    failure = StepFailure::Diverged;
    const Scalar dl = lambda_new - lambda;
    Evaluation ev = system.evaluate(state, lambda, true);
    if (!factor.factorize(free_block(ev.tangent, free_dofs), system.symmetric())) {
      return false;
    }
    const bool stable_start = factor.negative_pivots() == 0;
    const Vector dp = dl * prescribed;
    Vector rhs = dl * restrict(ev.load_rate, free_dofs) - restrict(ev.residual, free_dofs);
    if (dp.size() > 0) rhs -= free_prescribed_times(ev.tangent, free_dofs, fixed, dp);
    const Vector predicted = factor.solve(rhs);
    Vector trial = state;
    add_to(trial, free_dofs, predicted, 1.0);
    for (std::size_t i = 0; i < fixed.size(); ++i) {
      trial(fixed[i]) = lambda_new * prescribed(static_cast<Eigen::Index>(i));
    }
    const Vector start = state;
    const Vector predicted_state = trial;
    for (int it = 1; it <= options_.max_iterations; ++it) {
      ev = system.evaluate(trial, lambda_new, true);
      const Vector r = restrict(ev.residual, free_dofs);
      const Scalar scale = residual_scale(ev);
      if (!r.allFinite()) return false;
      if (!factor.factorize(free_block(ev.tangent, free_dofs), system.symmetric())) return false;
      const Scalar k_gross = stiffness_gross(ev.tangent, trial, free_dofs);
      pivots = factor.negative_pivots();
      if (stable_start && pivots > 0) {
        log::debug("  newton ", it, " at lambda = ", lambda_new, " met ", pivots,
                   " negative pivot(s) from a stable state");
        failure = StepFailure::Unstable;
        return false;
      }
      const Vector du = factor.solve(-r);
      if (!du.allFinite()) return false;
      const Scalar alpha = options_.line_search ? line_search(trial, du, r, lambda_new) : 1.0;
      add_to(trial, free_dofs, du, alpha);
      const Scalar increment = restrict(trial - start, free_dofs).norm();
      const Scalar rel = r.norm() / scale;
      log::debug("  newton ", it, " at lambda = ", lambda_new, ": |R| = ", r.norm(), " (scale ",
                 scale, ", round-off floor ", residual_floor(ev, k_gross, scale), "), |du| = ",
                 alpha * du.norm(),
                 " (limit ", correction_limit(increment, trial), "), line search ", alpha);
      // Converged: residual and correction within tolerance - or the
      // residual at its round-off floor, where the correction is round-off
      // amplified by the conditioning and cannot shrink further.
      if ((r.norm() <= options_.residual_tolerance * scale &&
           alpha * du.norm() <= correction_limit(increment, trial)) ||
          r.norm() <= residual_floor(ev, k_gross, scale)) {
        // The correction just applied is below tolerance: accept, and check
        // the residual of the accepted state.
        Evaluation final_ev = system.evaluate(trial, lambda_new, false);
        const Scalar final_norm = restrict(final_ev.residual, free_dofs).norm();
        const Scalar final_scale = residual_scale(final_ev);
        if (final_norm <= std::max(options_.residual_tolerance * final_scale,
                                   residual_floor(final_ev, k_gross, final_scale))) {
          const Scalar corrector = restrict(trial - predicted_state, free_dofs).norm();
          log::debug("  step to lambda = ", lambda_new, ": corrector / predictor = ",
                     corrector / predicted.norm());
          if (jump_test && predicted.norm() > 0.0 && corrector > predicted.norm()) {
            failure = StepFailure::Jumped;
            return false;
          }
          failure = StepFailure::None;
          state = trial;
          iterations = it;
          residual_out = final_norm / final_scale;
          full_residual = final_ev.residual;
          yielding = final_ev.yielding_points;
          system.commit(final_ev);
          return true;
        }
      }
      if (rel > 1.0e8) return false;  // diverging
    }
    return false;
  };

  // Newton iterations with contact from the converged state `state` (at
  // lambda) to lambda_new: a semismooth Newton method on the condensed
  // contact equations (Contact.hpp), from the converged state with the
  // prescribed displacements of lambda_new. Converged when the contact status
  // of every node is that of the previous iteration and the condensed
  // residual - equilibrium, the gap of the nodes in contact, the slip
  // conditions - is within tolerance (or at its round-off floor).
  const auto newton_contact = [&](Vector& state, Scalar lambda_new, int& iterations,
                                  Scalar& residual_out, Vector& full_residual, int& yielding,
                                  std::vector<ContactStatus>& status_out) -> bool {
    Vector trial = state;
    for (std::size_t i = 0; i < fixed.size(); ++i) {
      trial(fixed[i]) = lambda_new * prescribed(static_cast<Eigen::Index>(i));
    }
    std::vector<ContactStatus> previous;
    for (int it = 1; it <= options_.max_iterations; ++it) {
      Evaluation ev = system.evaluate(trial, lambda_new, true);
      const ContactProblem::Linearization lin =
          contact->linearize(trial, lambda_new, ev.residual, ev.tangent, state, lambda);
      const Scalar scale = std::max(residual_scale(ev), lin.force_scale);
      const Scalar norm = lin.rhs.norm();
      if (!std::isfinite(norm)) return false;
      const Scalar k_gross = stiffness_gross(ev.tangent, trial, free_dofs);
      const bool settled = !previous.empty() && lin.status == previous;
      int in_contact = 0;
      for (ContactStatus st : lin.status) in_contact += st != ContactStatus::Open ? 1 : 0;
      log::debug("  contact newton ", it, " at lambda = ", lambda_new, ": |R| = ", norm,
                 " (scale ", scale, "), ", in_contact, " node(s) in contact",
                 settled ? ", status settled" : "");
      if (settled && (norm <= options_.residual_tolerance * scale ||
                      norm <= residual_floor(ev, k_gross, scale))) {
        state = trial;
        iterations = it;
        residual_out = norm / scale;
        full_residual = ev.residual;
        yielding = ev.yielding_points;
        status_out = lin.status;
        system.commit(ev);
        return true;
      }
      if (norm > 1.0e8 * scale) return false;  // diverging
      previous = lin.status;
      Eigen::SparseLU<SparseMatrix> lu;
      lu.analyzePattern(lin.matrix);
      lu.factorize(lin.matrix);
      if (lu.info() != Eigen::Success) {
        contact_failure =
            "the contact system is singular: a body is not restrained against a rigid-body "
            "motion that the contact leaves free (a body held only by frictionless contact "
            "can slide along it) - restrain it with supports or prescribed displacements";
        return false;
      }
      const Vector du = lu.solve(lin.rhs);
      if (!du.allFinite()) return false;
      add_to(trial, free_dofs, du, 1.0);
    }
    contact_failure = "the contact status did not settle within max_iterations";
    return false;
  };

  bool path_completed = false;
  if (!arc) {
    // The load levels to pass through exactly: the turning points of the
    // load path, or the intermediate load factors and then lambda = 1.
    const bool path = !options_.load_path.empty();
    std::vector<Scalar> stations = options_.load_path;
    if (!path) {
      stations = options_.load_factors;
      stations.push_back(target);
    }
    Scalar path_scale = 0.0;
    for (Scalar v : stations) path_scale = std::max(path_scale, std::abs(v));
    std::size_t next_station = 0;
    // The path runs in legs between its turning points; each leg starts
    // with steps of 1/steps of its length, and s = direction * lambda grows
    // along it.
    Scalar direction = stations[0] > 0.0 ? 1.0 : -1.0;
    Scalar leg_step =
        (path ? std::abs(stations[0]) : target) / static_cast<Scalar>(options_.steps);
    step = leg_step;
    // A suspected critical point: the nearest s a step reached only by
    // meeting an unstable tangent or jumping, the kind of the evidence, and
    // the halvings it has cost from any state. Steps keep closing in on it
    // (a success lengthens the next step again); once max_cuts halvings have
    // failed to pass it the run stops there. A step that converges beyond it
    // proves it a stray iterate and clears it.
    Scalar critical_beyond = std::numeric_limits<Scalar>::infinity();
    StepFailure critical_kind = StepFailure::None;
    int critical_cuts = 0;
    // Likewise the nearest s that a step failed to reach by not converging:
    // successes short of it keep resetting the count of halvings in a row,
    // so without this a load beyond a limit - a plastic collapse, where there
    // is no equilibrium at all - would be crept towards until the step budget
    // ran out.
    Scalar stall_beyond = std::numeric_limits<Scalar>::infinity();
    int stall_cuts = 0;
    while (next_station < stations.size()) {
      if (static_cast<int>(result.steps.size()) >= options_.max_steps) {
        result.termination = path ? "the step budget (max_steps) ran out before the end of "
                                    "the load path"
                                  : "the step budget (max_steps) ran out before lambda = 1";
        break;
      }
      const Scalar station = stations[next_station];
      // A step that would stop just short of a station stretches to it.
      const Scalar lambda_new = direction * (station - lambda) <= step + 1.0e-9 * path_scale
                                    ? station
                                    : lambda + direction * step;
      Vector state = u;
      int iterations = 0;
      Scalar residual = 0.0;
      int pivots = -1;
      int yielding = 0;
      Vector full_residual;
      StepFailure failure = StepFailure::Diverged;
      bool ok = false;
      std::vector<ContactStatus> new_status;
      try {
        if (contact) {
          ok = newton_contact(state, lambda_new, iterations, residual, full_residual, yielding,
                              new_status);
          failure = ok ? StepFailure::None : StepFailure::Diverged;
        } else {
          ok = newton_load_control(state, lambda_new, iterations, residual, pivots,
                                   full_residual, yielding, failure);
        }
      } catch (const SolverError& ex) {
        log::debug("non-linear step at lambda = ", lambda_new, " failed: ", ex.what());
        ok = false;
        failure = StepFailure::Diverged;
      }
      if (!ok) {
        ++cuts_in_a_row;
        ++result.total_cuts;
        if (failure == StepFailure::Unstable || failure == StepFailure::Jumped) {
          critical_beyond = std::min(critical_beyond, direction * lambda_new);
          // A jump is the stronger evidence: it names a limit point.
          if (critical_kind != StepFailure::Jumped) critical_kind = failure;
          ++critical_cuts;
        } else {
          stall_beyond = std::min(stall_beyond, direction * lambda_new);
          ++stall_cuts;
        }
        if (critical_cuts > options_.max_cuts) {
          result.critical_bound = direction * critical_beyond;
          std::ostringstream os;
          if (critical_kind == StepFailure::Jumped) {
            os << "between lambda = " << lambda << " and " << result.critical_bound
               << " the path turns at a limit point: every step beyond it lands on a distant "
                  "branch of equilibria, a snap-through that load control cannot follow - the "
                  "arc-length method can";
          } else {
            os << "the tangent stiffness loses positive definiteness between lambda = " << lambda
               << " and " << result.critical_bound << ": a limit or bifurcation point"
               << (system.plastic() ? " or a plastic collapse (whose tangent is singular)" : "")
               << ", which load control cannot pass - the arc-length method follows the path "
                  "through a limit point";
          }
          result.termination = os.str();
          break;
        }
        if (stall_cuts > options_.max_cuts || cuts_in_a_row > options_.max_cuts) {
          std::ostringstream os;
          if (stall_cuts > options_.max_cuts) {
            result.unreached_load_factor = direction * stall_beyond;
            os << "no step converged beyond lambda = " << result.unreached_load_factor
               << " in " << stall_cuts << " halvings; the last converged load factor is "
               << lambda;
          } else {
            os << "a load step from lambda = " << lambda << " failed to converge after "
               << options_.max_cuts << " halvings (step " << step << ")";
          }
          os << ". The load may exceed a limit point - try the arc-length method";
          if (system.plastic()) {
            os << " - or a plastic collapse load, beyond which a material without hardening "
                  "has no equilibrium";
          }
          if (contact && !contact_failure.empty()) os << ". With contact: " << contact_failure;
          result.termination = os.str();
          break;
        }
        step *= 0.5;
        log::info("non-linear step to lambda = ", lambda_new,
                  failure == StepFailure::Unstable ? " met an unstable tangent"
                  : failure == StepFailure::Jumped ? " jumped to another branch"
                                                   : " did not converge",
                  "; halving to ", step);
        continue;
      }
      if (contact) {
        // The step's contact totals (its slip measured from the state it
        // started at), then commit the slip.
        step_contact = contact->pair_results(
            contact->node_results(state, lambda_new, full_residual, new_status, u, lambda));
        contact->commit(state, lambda_new, new_status, u, lambda);
        contact_status = new_status;
      }
      u = state;
      lambda = lambda_new;
      record(iterations, residual, cuts_in_a_row, pivots, full_residual, yielding);
      cuts_in_a_row = 0;
      if (direction * lambda > critical_beyond) {
        log::debug("load factor ", direction * critical_beyond,
                   " was a stray iterate, not a critical point");
        critical_beyond = std::numeric_limits<Scalar>::infinity();
        critical_kind = StepFailure::None;
        critical_cuts = 0;
      }
      if (direction * lambda > stall_beyond) {
        stall_beyond = std::numeric_limits<Scalar>::infinity();
        stall_cuts = 0;
      }
      if (iterations <= 3) step = std::min(leg_step, 1.5 * step);
      if (lambda_new == station) {
        ++next_station;
        if (next_station < stations.size()) {
          const Scalar turn = stations[next_station] > lambda ? 1.0 : -1.0;
          if (turn != direction) {
            // A turning point: a new leg, from its own first step.
            direction = turn;
            leg_step = std::abs(stations[next_station] - lambda) /
                       static_cast<Scalar>(options_.steps);
            step = leg_step;
            critical_beyond = std::numeric_limits<Scalar>::infinity();
            critical_kind = StepFailure::None;
            critical_cuts = 0;
            stall_beyond = std::numeric_limits<Scalar>::infinity();
            stall_cuts = 0;
          }
        }
      }
    }
    path_completed = next_station == stations.size();
  } else {
    // Crisfield's cylindrical arc-length method.
    bool finished = false;
    while (!finished) {
      if (static_cast<int>(result.steps.size()) >= options_.max_steps) {
        result.termination = "the step budget (max_steps) ran out before the target load factor";
        break;
      }
      Evaluation ev = system.evaluate(u, lambda, true);
      if (!factor.factorize(free_block(ev.tangent, free_dofs), system.symmetric())) {
        result.termination = "the tangent stiffness could not be factorised";
        break;
      }
      if (!result.steps.empty()) result.steps.back().negative_pivots = factor.negative_pivots();
      Vector uq = factor.solve(restrict(ev.load_rate, free_dofs));
      if (arc_length == 0.0) {
        arc_length = initial_step * uq.norm();
        first_arc = arc_length;
        if (!(arc_length > 0.0)) {
          result.termination = "the load produces no displacement: nothing to follow";
          break;
        }
      }
      // Predictor along the tangent, in the direction of the last increment.
      Scalar sign = 1.0;
      if (previous_increment.size() > 0 && previous_increment.dot(uq) < 0.0) sign = -1.0;
      Scalar dlambda = sign * arc_length / uq.norm();
      Vector du = dlambda * uq;
      bool converged = false;
      int iterations = 0;
      Scalar residual = 0.0;
      Vector step_residual;  // full residual of the converged state
      Evaluation accepted_ev;  // its evaluation, whose internal variables are committed
      try {
        for (int it = 1; it <= options_.max_iterations; ++it) {
          Vector trial = u;
          add_to(trial, free_dofs, du, 1.0);
          ev = system.evaluate(trial, lambda + dlambda, true);
          const Vector r = restrict(ev.residual, free_dofs);
          const Scalar scale = residual_scale(ev);
          if (!r.allFinite()) break;
          if (!factor.factorize(free_block(ev.tangent, free_dofs), system.symmetric())) break;
          const Scalar k_gross = stiffness_gross(ev.tangent, trial, free_dofs);
          const Vector ur = factor.solve(-r);
          uq = factor.solve(restrict(ev.load_rate, free_dofs));
          const Scalar a = uq.dot(uq);
          const Vector base = du + ur;
          const Scalar b = 2.0 * uq.dot(base);
          const Scalar c = base.dot(base) - arc_length * arc_length;
          const Scalar disc = b * b - 4.0 * a * c;
          if (!(disc >= 0.0) || !(a > 0.0)) break;
          const Scalar root = std::sqrt(disc);
          const Scalar s1 = (-b + root) / (2.0 * a);
          const Scalar s2 = (-b - root) / (2.0 * a);
          const Vector d1 = base + s1 * uq;
          const Vector d2 = base + s2 * uq;
          const bool first = d1.dot(du) >= d2.dot(du);
          const Scalar ds = first ? s1 : s2;
          const Vector correction = (first ? d1 : d2) - du;
          const Scalar rel = r.norm() / scale;
          du = first ? d1 : d2;
          dlambda += ds;
          Vector accepted = u;
          add_to(accepted, free_dofs, du, 1.0);
          if ((r.norm() <= options_.residual_tolerance * scale &&
               correction.norm() <= correction_limit(du.norm(), accepted)) ||
              r.norm() <= residual_floor(ev, k_gross, scale)) {
            Evaluation final_ev = system.evaluate(accepted, lambda + dlambda, false);
            const Scalar final_norm = restrict(final_ev.residual, free_dofs).norm();
            const Scalar final_scale = residual_scale(final_ev);
            if (final_norm <= std::max(options_.residual_tolerance * final_scale,
                                       residual_floor(final_ev, k_gross, final_scale))) {
              converged = true;
              iterations = it;
              residual = final_norm / final_scale;
              step_residual = final_ev.residual;
              accepted_ev = std::move(final_ev);
              break;
            }
          }
          if (rel > 1.0e8) break;
        }
      } catch (const SolverError& ex) {
        log::debug("arc-length step failed: ", ex.what());
        converged = false;
      }
      if (!converged) {
        ++cuts_in_a_row;
        ++result.total_cuts;
        arc_length *= 0.5;
        if (cuts_in_a_row > options_.max_cuts || arc_length < options_.min_arc_ratio * first_arc) {
          std::ostringstream os;
          os << "an arc-length step from lambda = " << lambda << " failed to converge after "
             << cuts_in_a_row << " halvings of the arc length";
          result.termination = os.str();
          break;
        }
        continue;
      }
      if (lambda + dlambda >= target) {
        // Land exactly on the target by load control from the last state.
        // (From the last converged state, whose internal variables are still
        // the committed ones.)
        Vector state = u;
        int its = 0;
        Scalar res = 0.0;
        int pivots = -1;
        int yielding = 0;
        Vector landing_residual;
        StepFailure failure = StepFailure::Diverged;
        bool ok = false;
        try {
          ok = newton_load_control(state, target, its, res, pivots, landing_residual, yielding,
                                   failure);
        } catch (const SolverError&) {
          ok = false;
        }
        if (ok) {
          u = state;
          lambda = target;
          record(its, res, cuts_in_a_row, pivots, landing_residual, yielding);
          finished = true;
          break;
        }
        // The path is not monotone near the target: take the arc-length step.
      }
      const int yielding = accepted_ev.yielding_points;
      system.commit(accepted_ev);
      add_to(u, free_dofs, du, 1.0);
      lambda += dlambda;
      previous_increment = du;
      // The inertia of this state is set by the next step's predictor, which
      // factorises the tangent here (or by the final factorisation).
      record(iterations, residual, cuts_in_a_row, -1, step_residual, yielding);
      cuts_in_a_row = 0;
      const Scalar factor_growth =
          std::sqrt(static_cast<Scalar>(options_.desired_iterations) /
                    static_cast<Scalar>(std::max(iterations, 1)));
      arc_length *= std::clamp(factor_growth, 0.5, 2.0);
      arc_length = std::clamp(arc_length, options_.min_arc_ratio * first_arc,
                              options_.max_arc_ratio * first_arc);
    }
  }

  // Final state: inertia, reactions, balance, stresses.
  const Evaluation final_ev = system.evaluate(u, lambda, true);
  if (factor.factorize(free_block(final_ev.tangent, free_dofs), system.symmetric()) &&
      !result.steps.empty()) {
    result.steps.back().negative_pivots = factor.negative_pivots();
  }
  result.symmetric_tangent = factor.symmetric();
  result.linear_solver = factor.name();
  result.displacement = u;
  result.load_factor = lambda;
  result.completed = arc ? lambda >= target * (1.0 - 1.0e-12) : path_completed;
  if (result.completed && result.termination.empty()) {
    result.termination = arc                               ? "reached the target load factor"
                         : !options_.load_path.empty() ? "completed the load path"
                                                           : "reached lambda = 1";
  }
  result.strain_energy = final_ev.energy;
  result.reactions = Vector::Zero(n);
  for (Index d : fixed) result.reactions(d) = final_ev.residual(d);

  // Contact at the final state (its slip committed). A rigid obstacle is a
  // support: the forces it exerts - the residual at the free DOFs of its
  // slave nodes - join the reactions in the force balance below.
  Vector obstacle = Vector::Zero(n);
  if (contact) {
    result.contact_nodes =
        contact->node_results(u, lambda, final_ev.residual, contact_status, u, lambda);
    result.contact_pairs = contact->pair_results(result.contact_nodes);
    for (const ContactProblem::Node& c : contact->nodes()) {
      if (!options_.contact.pairs[c.pair].rigid) continue;
      for (int comp : c.free) {
        const Index d = c.node * dim + comp;
        obstacle(d) = final_ev.residual(d);
      }
    }
    for (const std::string& note : contact->exclusions()) {
      result.warnings.push_back(note);
      log::warn("load case '", spec.name, "': ", note);
    }
  }

  // Force and moment balance of the deformed body.
  EquilibriumCheck& eq = result.equilibrium;
  Scalar force_scale = 0.0;
  Scalar moment_scale = 0.0;
  for (Index node = 0; node < mesh.num_nodes(); ++node) {
    const Vector3 x = mesh.node(node);
    Vector3 xd = x;
    Vector3 fa = Vector3::Zero();
    Vector3 fr = Vector3::Zero();
    for (int k = 0; k < dim; ++k) {
      xd(k) += u(node * dim + k);
      fa(k) = final_ev.external(node * dim + k);
      fr(k) = result.reactions(node * dim + k) + obstacle(node * dim + k);
    }
    eq.applied_force += fa;
    eq.reaction_force += fr;
    eq.applied_moment += xd.cross(fa);
    eq.reaction_moment += xd.cross(fr);
    force_scale += fa.norm() + fr.norm();
    moment_scale += xd.norm() * (fa.norm() + fr.norm());
  }
  if (dim == 2) {
    eq.applied_moment = Vector3(0.0, 0.0, eq.applied_moment.z());
    eq.reaction_moment = Vector3(0.0, 0.0, eq.reaction_moment.z());
  }
  eq.force_residual = eq.applied_force + eq.reaction_force;
  eq.moment_residual = eq.applied_moment + eq.reaction_moment;
  eq.relative_force_error = eq.force_residual.norm() / std::max(force_scale, 1.0e-300);
  eq.relative_moment_error = eq.moment_residual.norm() / std::max(moment_scale, 1.0e-300);

  const Index ne = mesh.num_elements();
  const int npe = mesh.nodes_per_elem();
  const int nv = voigt_components(dim);
  result.element_cauchy = Matrix::Zero(nv, ne);
  result.element_piola_kirchhoff = Matrix::Zero(nv, ne);
  result.element_von_mises = Vector::Zero(ne);
  result.element_cauchy_zz = Vector::Zero(ne);
  if (system.plastic()) result.element_plastic_strain = Vector::Zero(ne);
  if (logarithmic) {
    result.element_kirchhoff = Matrix::Zero(6, ne);
    result.element_log_strain = Matrix::Zero(6, ne);
  }
  const LoadCaseData& data = model_.load_case_data(load_case);
  const Vector* temperature = data.temperature.size() > 0 ? &data.temperature : nullptr;
  result.min_jacobian = std::numeric_limits<Scalar>::infinity();
  const int points = elastoplastic_points(model_);
  // A 3-D Voigt stress in the model's own Voigt order.
  const auto own_voigt = [&](const Vector6& v) {
    if (dim == 3) return Vector(v);
    Vector out(3);
    out << v(0), v(1), v(3);
    return out;
  };
  for (Index e = 0; e < ne; ++e) {
    const Index* nodes = mesh.element_nodes(e);
    Vector ue(dim * npe);
    for (int a = 0; a < npe; ++a) {
      for (int k = 0; k < dim; ++k) ue(dim * a + k) = u(nodes[a] * dim + k);
    }
    if (system.routed_through_return(e)) {
      const ElastoplasticStress st =
          elastoplastic_stress(model_, e, ue, system.committed(e), system.averaged(e),
                               temperature, lambda, options_.kinematics);
      result.element_cauchy.col(e) = own_voigt(st.cauchy);
      result.element_piola_kirchhoff.col(e) = own_voigt(st.piola_kirchhoff);
      if (dim == 2) result.element_cauchy_zz(e) = st.cauchy(2);
      result.element_von_mises(e) = st.von_mises;
      result.max_green_strain = std::max(result.max_green_strain, st.max_strain);
      result.min_jacobian = std::min(result.min_jacobian, st.min_jacobian);
      result.max_rotation = std::max(result.max_rotation, st.max_rotation);
      result.max_quadratic_strain = std::max(result.max_quadratic_strain, st.max_quadratic_strain);
      if (logarithmic) {
        result.element_kirchhoff.col(e) = st.kirchhoff;
        result.element_log_strain.col(e) = st.logarithmic_strain;
      }
      if (system.elastoplastic(e)) {
        result.element_plastic_strain(e) = st.max_equivalent_plastic_strain;
        result.max_plastic_strain =
            std::max(result.max_plastic_strain, st.max_equivalent_plastic_strain);
        result.plastic_points += st.plastic_points;
        result.total_points += points;
      }
      continue;
    }
    const TotalLagrangianStress st =
        total_lagrangian_stress(model_, e, ue, options_.law, temperature, lambda);
    result.element_cauchy.col(e) = st.cauchy;
    result.element_piola_kirchhoff.col(e) = st.piola_kirchhoff;
    result.element_cauchy_zz(e) = st.cauchy_zz;
    result.element_von_mises(e) =
        dim == 3 ? von_mises(st.cauchy, StressState::ThreeDimensional, 0.0)
                 : von_mises_plane(st.cauchy(0), st.cauchy(1), st.cauchy(2), st.cauchy_zz);
    result.max_green_strain = std::max(result.max_green_strain, st.max_green_strain);
    result.min_jacobian = std::min(result.min_jacobian, st.min_jacobian);
  }

  // What bears on the validity of the run.
  const auto warn = [&](const std::string& text) {
    result.warnings.push_back(text);
    log::warn("load case '", spec.name, "': ", text);
  };
  if (small && result.max_quadratic_strain > kQuadraticWarning * result.max_green_strain &&
      result.max_green_strain > 0.0) {
    std::ostringstream os;
    os << "the small-strain kinematics neglects the quadratic part of the Green strain, "
          "H^T H / 2, which reaches "
       << result.max_quadratic_strain / result.max_green_strain
       << " of the largest strain (largest rotation " << result.max_rotation
       << " rad); \"finite\" kinematics models it";
    warn(os.str());
  }
  // (The logarithmic kinematics is meant for large strain.)
  if ((small || system.plastic()) && !logarithmic && result.max_green_strain > kStrainWarning) {
    std::ostringstream os;
    os << "the largest strain is " << result.max_green_strain << ", beyond the small strains "
       << (system.plastic() ? "the elastoplastic law assumes"
                            : "of the small-strain (linear elastic) theory");
    warn(os.str());
  }
  if (system.plastic() && model_.stress_state() != StressState::PlaneStress) {
    const ElementType type = mesh.element_type();
    if (type == ElementType::Tri3 || type == ElementType::Tet4) {
      warn("constant-strain elements (" + to_string(type) +
           ") can lock under the isochoric flow of a fully plastic state in plane strain "
           "and 3-D, depending on the mesh pattern (a collapse load comes out too high); "
           "check one against Q4 or Hex8 with mean dilatation, or Tet10");
    } else if ((type == ElementType::Quad4 || type == ElementType::Hex8) &&
               !result.mean_dilatation) {
      warn("without mean dilatation the fully integrated " + to_string(type) +
           " locks under isochoric plastic flow in plane strain and 3-D: a collapse load "
           "comes out too high and the plastic plateau keeps rising");
    }
  }
  if (!result.completed) {
    log::warn("non-linear analysis of load case '", spec.name, "' stopped at lambda = ",
              lambda, ": ", result.termination);
  } else {
    log::info("non-linear analysis of load case '", spec.name, "': ", result.steps.size(),
              " step(s), ", result.total_iterations, " iteration(s), ", result.total_cuts,
              " cut(s), lambda = ", lambda, ", largest ",
              small ? "strain " : logarithmic ? "log strain " : "Green strain ",
              result.max_green_strain);
    if (system.plastic()) {
      log::info("  ", result.plastic_points, " of ", result.total_points,
                " elastoplastic integration point(s) have yielded; largest plastic strain ",
                result.max_plastic_strain);
    }
  }
  return result;
}

}  // namespace sparlab
