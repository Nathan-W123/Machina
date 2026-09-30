/// \file verify_dynamics.cpp
/// \brief Verification of the transient and harmonic analyses against exact
///        solutions of the discrete equations and of the continuum.
///
/// Studies:
///   * `transient-modal`      the HHT-alpha transient of Q4 and Hex8 cantilevers
///                            against the exact solution of the same discrete
///                            equations by modal superposition - every mode of
///                            a dense eigensolve, each advanced by its own
///                            scalar HHT-alpha recursion - with consistent and
///                            lumped mass, the trapezoidal rule and alpha < 0,
///                            Rayleigh damping, a sudden, a harmonic and a
///                            released load; and the trapezoidal rule's energy
///                            balance repeated in 80-bit arithmetic, to show
///                            that what remains of it is round-off;
///   * `rod-harmonic`         the steady harmonic response of a fixed-free rod,
///                            driven at its end or at its base, undamped, with
///                            structural and with Rayleigh damping: to
///                            round-off against the exact solution of the
///                            discrete equations (their dispersion relation),
///                            and at second order in the element length
///                            against the exact continuum solution;
///   * `rod-transient`        the same rod under an end force ramped up
///                            smoothly: its end displacement against the
///                            continuum's modal series, converging at second
///                            order as the element length and the time step
///                            shrink together;
///   * `nonlinear-oscillator` one element in uniaxial strain - a single degree
///                            of freedom - under a sudden load, with the
///                            Saint Venant-Kirchhoff law at large strain and
///                            with small-strain J2 plasticity and linear
///                            hardening: to the Newton tolerance against the
///                            scalar HHT-alpha recursion of the same equation,
///                            and at second order in the step against the
///                            exact motion (a Runge-Kutta reference for the
///                            elastic law, the piecewise closed form for the
///                            plastic one);
///   * `explicit-rod`         the rod-transient case integrated by explicit
///                            central differences (ExplicitDynamics.hpp) with
///                            the lumped mass, Q4 through the generic element
///                            dispatch and Hex8 through the dedicated kernel:
///                            second order against the continuum; the element
///                            eigenvalue and power-iteration step estimates
///                            against the exact limit; a slow damped ramp
///                            reaching the static displacement.

#include "VerifySupport.hpp"

#include "AppSupport.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/fem/Assembler.hpp"
#include "sparlab/fem/Dynamics.hpp"
#include "sparlab/fem/ExplicitDynamics.hpp"
#include "sparlab/io/CsvWriter.hpp"
#include "sparlab/mesh/StructuredMesh.hpp"

#include <Eigen/Eigenvalues>
#include <Eigen/SparseCholesky>

#include <algorithm>
#include <cmath>
#include <complex>
#include <functional>
#include <limits>
#include <sstream>
#include <vector>

namespace sparlab {
namespace verify {
namespace {

constexpr Scalar kPi = 3.14159265358979323846;
using Complex = std::complex<Scalar>;

std::string fmt(Scalar v, int digits = 6) {
  return std::isnan(v) ? std::string() : app::format(v, digits);
}

Selector x_box(Scalar xmin, Scalar xmax) {
  Selector s;
  s.kind = SelectorKind::Box;
  s.xmin = xmin;
  s.xmax = xmax;
  return s;
}

std::string element_name(ElementType type) { return type == ElementType::Quad4 ? "Q4" : "Hex8"; }

Vector free_part(const FemModel& model, const Vector& full) {
  const std::vector<Index>& free = model.dofs().free_dofs();
  Vector out(static_cast<Eigen::Index>(free.size()));
  for (std::size_t i = 0; i < free.size(); ++i) out(static_cast<Eigen::Index>(i)) = full(free[i]);
  return out;
}

/// The mean x-displacement monitor over the nodes of the plane x = at, or
/// the sum of the x-reactions over them.
DynamicMonitor x_monitor(const std::string& name, Scalar at, DynamicMonitor::Quantity q) {
  DynamicMonitor m;
  m.name = name;
  m.region.name = name;
  m.region.members.push_back(x_box(at, at));
  m.component = 0;
  m.quantity = q;
  return m;
}

// ---------------------------------------------------------------------------
// Transient against the exact discrete modal solution
// ---------------------------------------------------------------------------

/// A cantilever 1 m long, 0.1 m deep (Q4 in plane stress, 0.01 m thick) or
/// 0.1 m x 0.1 m (Hex8), steel, clamped at x = 0 with a tip load of -100 N
/// in y shared by the end nodes.
FemModel cantilever_model(ElementType type) {
  const bool plane = type == ElementType::Quad4;
  StructuredMeshSpec spec;
  spec.nx = plane ? 20 : 10;
  spec.ny = plane ? 4 : 2;
  spec.nz = 2;
  spec.lx = 1.0;
  spec.ly = 0.1;
  spec.lz = 0.1;
  FemModel model(plane ? make_structured_quad_mesh(spec) : make_structured_hex_mesh(spec),
                 IsotropicMaterial(200.0e9, 0.3, 7850.0, "steel"), plane ? 0.01 : 1.0,
                 plane ? StressState::PlaneStress : StressState::ThreeDimensional,
                 IntegrationOptions());
  DisplacementConstraint root;
  root.region.name = "root";
  root.region.members.push_back(x_box(-std::numeric_limits<Scalar>::infinity(), 0.0));
  root.fix_x = root.fix_y = true;
  root.fix_z = !plane;
  model.constraints().push_back(root);
  LoadCaseSpec lc;
  lc.name = "tip";
  PointLoadSpec tip;
  tip.region.name = "tip";
  tip.region.members.push_back(x_box(1.0, std::numeric_limits<Scalar>::infinity()));
  tip.force = Vector3(0.0, -100.0, 0.0);
  tip.distribute_total = true;
  lc.point_loads.push_back(tip);
  model.load_case_specs().push_back(lc);
  model.finalize();
  return model;
}

using Extended = long double;
using ExtendedMatrix = Eigen::Matrix<Extended, Eigen::Dynamic, Eigen::Dynamic>;
using ExtendedVector = Eigen::Matrix<Extended, Eigen::Dynamic, 1>;

/// Every mode of the free DOFs, from a dense generalized eigensolve of
/// (K, M) in 80-bit extended precision - so that the reference's own
/// round-off lies far below the double-precision solver's.
struct DenseModes {
  ExtendedVector omega2;  ///< ascending [1/s^2]
  ExtendedMatrix phi;     ///< free DOFs x modes, M-orthonormal
};

DenseModes dense_modes(const Assembler& assembler, MassType type) {
  const ExtendedMatrix k =
      Matrix(assembler.reduce_free_free(assembler.assemble_stiffness())).cast<Extended>();
  const ExtendedMatrix m =
      Matrix(assembler.reduce_free_free(assembler.assemble_mass(type))).cast<Extended>();
  Eigen::GeneralizedSelfAdjointEigenSolver<ExtendedMatrix> es(k, m);
  if (es.info() != Eigen::Success) throw SolverError("dense generalized eigensolve failed");
  return {es.eigenvalues(), es.eigenvectors()};
}

/// The HHT-alpha method on one mode, q'' + c q' + w2 q = p(t_n), in its
/// acceleration (predictor-corrector) form - not the displacement form of
/// the solver: Newmark predictors, the equation of motion at the weighted
/// point solved for the new acceleration, then the correctors.
template <typename T>
std::vector<T> hht_mode(T w2, T c, const std::vector<T>& p, T dt, T alpha, T q0, T v0) {
  const T beta = T(0.25) * (T(1) - alpha) * (T(1) - alpha);
  const T gamma = T(0.5) - alpha;
  const T w = T(1) + alpha;
  T q = q0;
  T v = v0;
  T a = p[0] - c * v0 - w2 * q0;
  std::vector<T> out{q};
  for (std::size_t n = 1; n < p.size(); ++n) {
    const T q_pred = q + dt * v + dt * dt * (T(0.5) - beta) * a;
    const T v_pred = v + dt * (T(1) - gamma) * a;
    const T a1 = (w * p[n] - alpha * p[n - 1] - w * (c * v_pred + w2 * q_pred) +
                  alpha * (c * v + w2 * q)) /
                 (T(1) + w * c * gamma * dt + w * w2 * beta * dt * dt);
    q = q_pred + beta * dt * dt * a1;
    v = v_pred + gamma * dt * a1;
    a = a1;
    out.push_back(q);
  }
  return out;
}

/// The undamped trapezoidal recursion from rest under f A(t_n), in the
/// arithmetic T: the largest |E_0 + W - T - U| over the energies involved.
template <typename T>
Scalar trapezoidal_balance(const SparseMatrix& k, const SparseMatrix& m, const Vector& f,
                           const std::vector<Scalar>& amplitude, Scalar dt) {
  using Sparse = Eigen::SparseMatrix<T>;
  using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;
  const Sparse kt = k.cast<T>();
  const Sparse mt = m.cast<T>();
  const Vec ft = f.cast<T>();
  const T h = static_cast<T>(dt);
  const T c0 = T(4) / (h * h);
  const T c1 = T(4) / h;
  const T c3 = T(2) / h;
  const Sparse keff = c0 * mt + kt;
  Eigen::SimplicialLDLT<Sparse> step(keff);
  Eigen::SimplicialLDLT<Sparse> mass(mt);
  if (step.info() != Eigen::Success || mass.info() != Eigen::Success) {
    throw SolverError("trapezoidal_balance: factorisation failed");
  }
  const Eigen::Index n = ft.size();
  Vec u = Vec::Zero(n);
  Vec v = Vec::Zero(n);
  T a0 = static_cast<T>(amplitude[0]);
  Vec a = mass.solve(Vec(a0 * ft));
  T work = 0;
  T worst = 0;
  T scale = 0;
  for (std::size_t s = 1; s < amplitude.size(); ++s) {
    const T a1 = static_cast<T>(amplitude[s]);
    const Vec a_hat = -c0 * u - c1 * v - a;
    const Vec v_hat = -c3 * u - v;
    const Vec u1 = step.solve(Vec(a1 * ft - mt * a_hat));
    work += T(0.5) * (u1 - u).dot(Vec((a0 + a1) * ft));
    a = c0 * u1 + a_hat;
    v = c3 * u1 + v_hat;
    u = u1;
    a0 = a1;
    const T kinetic = T(0.5) * v.dot(mt * v);
    const T strain = T(0.5) * u.dot(kt * u);
    scale = std::max({scale, kinetic, strain, T(std::abs(work))});
    worst = std::max(worst, T(std::abs(work - kinetic - strain)));
  }
  return static_cast<Scalar>(worst / scale);
}

}  // namespace

StudyOutcome study_transient_modal(const std::string& out_dir, json::Value& summary) {
  CsvWriter csv(path_join(out_dir, "transient_modal.csv"),
                {"element", "mass", "case", "alpha", "steps", "free_dofs",
                 "difference_first_half[-]", "max_relative_difference[-]", "energy_balance[-]",
                 "numerical_dissipation_over_work[-]", "kappa_eff_eps[-]"});
  json::Value block = json::Value::make_object();
  json::Value records = json::Value::make_array();
  Scalar worst = 0.0;
  Scalar worst_half = 0.0;  // over the first half of each run
  Scalar worst_balance = 0.0;
  bool dissipation_positive = true;
  const Scalar eps = std::numeric_limits<Scalar>::epsilon();

  struct Case {
    std::string name;
    Scalar alpha;
    bool damped;
    int load;  ///< 0 sudden, 1 harmonic, 2 released from the static state
  };
  const std::vector<Case> cases = {{"sudden load undamped", 0.0, false, 0},
                                   {"harmonic load Rayleigh damping", 0.0, true, 1},
                                   {"harmonic load Rayleigh damping", -0.1, true, 1},
                                   {"sudden load undamped", -0.3, false, 0},
                                   {"released from the static state", -0.05, true, 2}};
  json::Value precision = json::Value::make_object();
  for (const ElementType type : {ElementType::Quad4, ElementType::Hex8}) {
    const FemModel model = cantilever_model(type);
    const Assembler assembler(model);
    for (const MassType mass : {MassType::Consistent, MassType::Lumped}) {
      const DenseModes modes = dense_modes(assembler, mass);
      const Scalar w1 = std::sqrt(static_cast<Scalar>(modes.omega2(0)));
      const Scalar f1 = w1 / (2.0 * kPi);
      const Scalar period = 1.0 / f1;
      const Vector f = free_part(model, model.load_vectors()[0]);
      const ExtendedVector modal_load = modes.phi.transpose() * f.cast<Extended>();
      for (const Case& c : cases) {
        TransientOptions options;
        options.mass_type = mass;
        options.alpha = c.alpha;
        options.time_step = period / 40.0;
        options.end_time = 120 * options.time_step;
        options.snapshot_every = 1;
        if (c.damped) {
          // 2 % of critical in the first mode, from both terms.
          options.mass_damping = 0.02 * w1;
          options.stiffness_damping = 0.02 / w1;
        }
        if (c.load == 1) {
          options.amplitude.kind = Amplitude::Kind::Harmonic;
          options.amplitude.frequency = 1.3 * f1;
          options.amplitude.phase = 0.4;
        } else if (c.load == 2) {
          options.start = TransientOptions::Start::Static;
          options.amplitude.kind = Amplitude::Kind::Table;
          options.amplitude.times = {0.0, 0.5 * options.time_step};
          options.amplitude.values = {1.0, 0.0};
        }
        const TransientResult r = solve_transient(model, assembler, 0, options);
        const int steps = r.num_steps;
        const Scalar dt = options.time_step;
        ExtendedMatrix q(modes.omega2.size(), steps + 1);
        for (Eigen::Index j = 0; j < modes.omega2.size(); ++j) {
          std::vector<Extended> p;
          for (int n = 0; n <= steps; ++n) {
            p.push_back(modal_load(j) * static_cast<Extended>(options.amplitude.value(n * dt)));
          }
          const Extended damping = static_cast<Extended>(options.mass_damping) +
                                   static_cast<Extended>(options.stiffness_damping) *
                                       modes.omega2(j);
          const Extended q0 = c.load == 2 ? p[0] / modes.omega2(j) : Extended(0);
          const std::vector<Extended> history =
              hht_mode<Extended>(modes.omega2(j), damping, p, static_cast<Extended>(dt),
                                 static_cast<Extended>(c.alpha), q0, Extended(0));
          for (int n = 0; n <= steps; ++n) q(j, n) = history[static_cast<std::size_t>(n)];
        }
        // The largest difference over the first half of the run and over
        // all of it: round-off accumulates step by step.
        Scalar difference = 0.0;
        Scalar difference_half = 0.0;
        Scalar scale = 0.0;
        for (int n = 0; n <= steps; ++n) {
          const Vector reference = (modes.phi * q.col(n)).cast<Scalar>();
          const Vector ours =
              free_part(model, r.snapshots[static_cast<std::size_t>(n)].displacement);
          difference = std::max(difference, (ours - reference).cwiseAbs().maxCoeff());
          if (2 * n <= steps) difference_half = difference;
          scale = std::max(scale, reference.cwiseAbs().maxCoeff());
        }
        const Scalar relative = difference / scale;
        const Scalar relative_half = difference_half / scale;
        // The round-off scale of the solver's own solves: the condition
        // number of its effective stiffness c0 M + (1 + alpha)(c3 C + K),
        // the ratio of its extreme modal entries, times eps.
        const HhtParameters hht = HhtParameters::from_alpha(c.alpha);
        const Scalar c0 = 1.0 / (hht.beta * dt * dt);
        const Scalar c3 = hht.gamma / (hht.beta * dt);
        const auto effective = [&](Scalar w2) {
          return c0 + (1.0 + c.alpha) * (c3 * (options.mass_damping +
                                               options.stiffness_damping * w2) + w2);
        };
        const Scalar kappa =
            effective(static_cast<Scalar>(modes.omega2(modes.omega2.size() - 1))) /
            effective(static_cast<Scalar>(modes.omega2(0)));
        const Scalar dissipation =
            r.steps.back().external_work != 0.0
                ? r.numerical_dissipation / std::abs(r.steps.back().external_work)
                : std::nan("");
        worst = std::max(worst, relative);
        worst_half = std::max(worst_half, relative_half);
        if (c.alpha == 0.0) {
          worst_balance = std::max(worst_balance, r.energy_balance_error);
        } else if (c.load == 0) {
          // HHT-alpha dissipates the energy of the modes it damps.
          dissipation_positive = dissipation_positive && r.numerical_dissipation > 0.0;
        }
        csv.raw_row({element_name(type), to_string(mass), c.name, fmt(c.alpha, 3),
                     fmt(static_cast<Scalar>(steps)),
                     fmt(static_cast<Scalar>(model.dofs().num_free())), fmt(relative_half, 4),
                     fmt(relative, 4),
                     fmt(r.energy_balance_error, 4), fmt(dissipation, 4), fmt(kappa * eps, 4)});
        json::Value rec = json::Value::make_object();
        rec.set("element", json::Value::make_string(element_name(type)));
        rec.set("mass", json::Value::make_string(to_string(mass)));
        rec.set("case", json::Value::make_string(c.name));
        rec.set("alpha", json::Value::make_number(c.alpha));
        rec.set("difference_first_half", json::Value::make_number(relative_half));
        rec.set("max_relative_difference", json::Value::make_number(relative));
        rec.set("energy_balance", json::Value::make_number(r.energy_balance_error));
        rec.set("numerical_dissipation_over_work", json::Value::make_number(dissipation));
        rec.set("kappa_eff_eps", json::Value::make_number(kappa * eps));
        records.push_back(rec);
      }
      // The trapezoidal rule's energy balance of the sudden load, repeated in
      // double and 80-bit extended precision: a balance error that falls with
      // the unit round-off is round-off.
      if (type == ElementType::Quad4 && mass == MassType::Consistent) {
        const SparseMatrix k = assembler.reduce_free_free(assembler.assemble_stiffness());
        const SparseMatrix m = assembler.reduce_free_free(assembler.assemble_mass(mass));
        const std::vector<Scalar> steps_amplitude(121, 1.0);
        const Scalar in_double = trapezoidal_balance<double>(k, m, f, steps_amplitude, period / 40.0);
        const Scalar in_extended =
            trapezoidal_balance<long double>(k, m, f, steps_amplitude, period / 40.0);
        const Scalar eps_extended = static_cast<Scalar>(std::numeric_limits<long double>::epsilon());
        precision.set("model", json::Value::make_string("Q4 cantilever, consistent mass, sudden "
                                                        "load, trapezoidal rule, 120 steps"));
        precision.set("balance_double", json::Value::make_number(in_double));
        precision.set("balance_extended", json::Value::make_number(in_extended));
        precision.set("eps_double", json::Value::make_number(eps));
        precision.set("eps_extended", json::Value::make_number(eps_extended));
        precision.set("balance_ratio", json::Value::make_number(in_double / in_extended));
        precision.set("eps_ratio", json::Value::make_number(eps / eps_extended));
      }
    }
  }
  csv.close();
  block.set("kind", json::Value::make_string(
                        "verification (exact solution of the discrete equations by modal "
                        "superposition)"));
  block.set("cases", records);
  block.set("energy_precision", precision);
  block.set("note",
            json::Value::make_string(
                "Cantilever 1 m x 0.1 m, steel (E = 200 GPa, nu = 0.3, rho = 7850 kg/m^3), Q4 "
                "20 x 4 in plane stress 0.01 m thick and Hex8 10 x 2 x 2, clamped, tip load "
                "-100 N; 120 steps of T1/40. The reference: every mode of a dense generalized "
                "eigensolve (K, M) of the free DOFs, each integrated by the scalar HHT-alpha "
                "recursion in acceleration form, summed. Rayleigh damping 2 % of critical in "
                "the first mode from each of a M and b K; the harmonic load at 1.3 f1; the "
                "released case starts from the static solution and the load drops to zero "
                "within the first step. The reference is computed in 80-bit extended "
                "precision (a double-precision modal sum is itself only good to about "
                "omega_max^2/omega_1^2 eps, 1e-10 here), so the difference is the solver's "
                "double-precision round-off, accumulating from step to step (compare the "
                "first half of each run with the whole). kappa_eff_eps: the condition number "
                "of the effective stiffness (the ratio of its extreme modal entries) times "
                "eps, the round-off scale of one step's solve. Energy balance: "
                "the largest |E_0 + W - T - U - D| over the energies involved, exact to "
                "round-off for alpha = 0; the numerical dissipation of alpha < 0 over the "
                "external work."));
  summary.set("transient_modal", block);

  const Scalar ratio = precision.find("balance_ratio")->number_value();
  block.set("largest_difference", json::Value::make_number(worst));
  block.set("largest_difference_first_half", json::Value::make_number(worst_half));
  std::ostringstream note;
  note << "largest difference " << app::format(worst, 3)
       << " (over the first half of the runs " << app::format(worst_half, 3)
       << "); largest trapezoidal energy balance "
       << app::format(worst_balance, 3)
       << "; in 80-bit arithmetic it falls " << app::format(ratio, 3) << " times (eps "
       << app::format(precision.find("eps_ratio")->number_value(), 3) << " times)";
  StudyOutcome outcome;
  outcome.name = "HHT-alpha transient vs exact discrete modal solution (Q4, Hex8)";
  outcome.kind = "verification";
  outcome.metric =
      "largest relative difference over steps, models, masses and cases";
  outcome.value = worst;
  outcome.tolerance = 1.0e-9;
  outcome.passed = worst <= 1.0e-9 && worst_balance <= 1.0e-10 && dissipation_positive &&
                   ratio > 100.0;
  outcome.note = note.str();
  return outcome;
}

// ---------------------------------------------------------------------------
// Fixed-free rod: harmonic response
// ---------------------------------------------------------------------------

namespace {

/// A steel rod 1 m long, 0.05 m square, nu = 0, its lateral displacements
/// held (an exactly one-dimensional model): Q4 in plane stress 0.05 m
/// thick with one element across, or Hex8 with one element across. The
/// end x = 0 is held (or driven, `base`), the end x = L carries `force`.
struct Rod {
  Scalar length = 1.0;
  Scalar side = 0.05;
  Scalar youngs = 200.0e9;
  Scalar density = 7850.0;
  Scalar area() const { return side * side; }
  Scalar wave_speed() const { return std::sqrt(youngs / density); }
  /// The first natural frequency of the fixed-free rod, c / (4 L) [Hz].
  Scalar f1() const { return wave_speed() / (4.0 * length); }
};

FemModel rod_model(const Rod& rod, ElementType type, Index n, Scalar force, Scalar base) {
  const bool plane = type == ElementType::Quad4;
  StructuredMeshSpec spec;
  spec.nx = n;
  spec.ny = 1;
  spec.nz = 1;
  spec.lx = rod.length;
  spec.ly = rod.side;
  spec.lz = rod.side;
  FemModel model(plane ? make_structured_quad_mesh(spec) : make_structured_hex_mesh(spec),
                 IsotropicMaterial(rod.youngs, 0.0, rod.density, "steel"),
                 plane ? rod.side : 1.0,
                 plane ? StressState::PlaneStress : StressState::ThreeDimensional,
                 IntegrationOptions());
  DisplacementConstraint lateral;
  lateral.region.name = "lateral";
  lateral.region.members.push_back(Selector());
  lateral.fix_y = true;
  lateral.fix_z = !plane;
  model.constraints().push_back(lateral);
  DisplacementConstraint root;
  root.region.name = "root";
  root.region.members.push_back(x_box(-std::numeric_limits<Scalar>::infinity(), 0.0));
  root.set(0, true, base);
  model.constraints().push_back(root);
  LoadCaseSpec lc;
  lc.name = "drive";
  if (force != 0.0) {
    PointLoadSpec end;
    end.region.name = "end";
    end.region.members.push_back(x_box(rod.length, std::numeric_limits<Scalar>::infinity()));
    end.force = Vector3(force, 0.0, 0.0);
    end.distribute_total = true;
    lc.point_loads.push_back(end);
  } else {
    lc.prescribed_displacement_only = true;
  }
  model.load_case_specs().push_back(lc);
  model.finalize();
  return model;
}

/// The damping of a harmonic case: the complex modulus and density that
/// K (1 + i eta) - omega^2 M + i omega (a M + b K) amounts to,
/// E* = E (1 + i eta + i omega b), rho* = rho (1 - i a / omega).
struct HarmonicDamping {
  Scalar eta = 0.0;
  Scalar a = 0.0;
  Scalar b = 0.0;
  Complex modulus(Scalar youngs, Scalar omega) const {
    return youngs * Complex(1.0, eta + omega * b);
  }
  Complex density(Scalar rho, Scalar omega) const { return rho * Complex(1.0, -a / omega); }
};

/// The exact solution of the discrete equations of n linear rod elements
/// at omega: u_j = g cos(j theta) + D sin(j theta), theta from the
/// dispersion relation of the interior stencil,
///   consistent mass: cos theta = (6 - 2 W) / (6 + W),
///   lumped mass:     cos theta = 1 - W / 2,        W = omega^2 h^2 rho* / E*,
/// u_0 = g, and D from the equation of the end node under the force F.
/// Returns u_0..u_n and the reaction at node 0.
std::vector<Complex> rod_discrete(const Rod& rod, Index n, bool lumped, Complex e_star,
                                  Complex rho_star, Scalar omega, Scalar force, Scalar base,
                                  Complex& reaction) {
  const Scalar h = rod.length / static_cast<Scalar>(n);
  const Scalar area = rod.area();
  const Complex stiff = e_star * area / h;               // E* A / h
  const Complex mass = omega * omega * rho_star * area * h;  // omega^2 rho* A h
  const Complex w = omega * omega * h * h * rho_star / e_star;
  const Complex cos_theta = lumped ? 1.0 - 0.5 * w : (6.0 - 2.0 * w) / (6.0 + w);
  const Complex theta = std::acos(cos_theta);
  const auto c = [&](Index j) { return std::cos(static_cast<Scalar>(j) * theta); };
  const auto s = [&](Index j) { return std::sin(static_cast<Scalar>(j) * theta); };
  // End node n: stiff (u_n - u_{n-1}) - mass (2 u_n + u_{n-1}) / 6 = F
  // (consistent) or stiff (u_n - u_{n-1}) - mass u_n / 2 = F (lumped).
  const auto end_row = [&](Complex un, Complex un1) {
    return lumped ? stiff * (un - un1) - 0.5 * mass * un
                  : stiff * (un - un1) - mass * (2.0 * un + un1) / 6.0;
  };
  const Complex from_base = end_row(base * c(n), base * c(n - 1));
  const Complex per_d = end_row(s(n), s(n - 1));
  const Complex d = (force - from_base) / per_d;
  std::vector<Complex> u;
  for (Index j = 0; j <= n; ++j) u.push_back(base * c(j) + d * s(j));
  // The reaction at node 0: its row of the dynamic stiffness times u.
  reaction = lumped ? stiff * (u[0] - u[1]) - 0.5 * mass * u[0]
                    : stiff * (u[0] - u[1]) - mass * (2.0 * u[0] + u[1]) / 6.0;
  return u;
}

/// The exact continuum solution at x:
/// u = g cos(k (L - x)) / cos(k L) + F sin(k x) / (E* A k cos(k L)),
/// k = omega sqrt(rho* / E*).
Complex rod_continuum(const Rod& rod, Complex e_star, Complex rho_star, Scalar omega,
                      Scalar force, Scalar base, Scalar x) {
  const Complex k = omega * std::sqrt(rho_star / e_star);
  const Complex ckl = std::cos(k * rod.length);
  return base * std::cos(k * (rod.length - x)) / ckl +
         force * std::sin(k * x) / (e_star * rod.area() * k * ckl);
}

}  // namespace

StudyOutcome study_rod_harmonic(const std::string& out_dir, json::Value& summary) {
  const Rod rod;
  const Scalar f1 = rod.f1();
  const Scalar w1 = 2.0 * kPi * f1;
  struct Case {
    std::string name;
    Scalar ratio;  ///< f / f1
    HarmonicDamping damping;
    bool base;     ///< driven at the base (else by the end force)
  };
  HarmonicDamping none;
  HarmonicDamping structural;
  structural.eta = 0.02;
  HarmonicDamping rayleigh;  // 1 % of critical at f1 from each term
  rayleigh.a = 0.01 * w1;
  rayleigh.b = 0.01 / w1;
  const std::vector<Case> cases = {{"end force undamped", 0.5, none, false},
                                   {"end force undamped", 2.5, none, false},
                                   {"end force eta 0.02 at resonance", 1.0, structural, false},
                                   {"end force eta 0.02 at resonance", 3.0, structural, false},
                                   {"end force Rayleigh at resonance", 1.0, rayleigh, false},
                                   {"base motion undamped", 1.5, none, true},
                                   {"base motion eta 0.02 at resonance", 1.0, structural, true}};
  const std::vector<Index> ladder = {10, 20, 40, 80, 160};
  const Scalar force = 1000.0;
  const Scalar base = 1.0e-5;

  CsvWriter csv(path_join(out_dir, "rod_harmonic.csv"),
                {"element", "mass", "case", "f_over_f1", "n", "discrete_difference[-]",
                 "reaction_difference[-]", "continuum_error[-]", "order[-]", "end_amplitude[m]",
                 "exact_end_amplitude[m]"});
  json::Value records = json::Value::make_array();
  Scalar worst_discrete = 0.0;
  Scalar worst_order = 1.0e300;
  Scalar worst_error = 0.0;
  for (const ElementType type : {ElementType::Quad4, ElementType::Hex8}) {
    for (const MassType mass : {MassType::Consistent, MassType::Lumped}) {
      const bool lumped = mass == MassType::Lumped;
      for (const Case& c : cases) {
        const Scalar freq = c.ratio * f1;
        const Scalar omega = 2.0 * kPi * freq;
        const Complex e_star = c.damping.modulus(rod.youngs, omega);
        const Complex rho_star = c.damping.density(rod.density, omega);
        const Scalar f_load = c.base ? 0.0 : force;
        const Scalar g = c.base ? base : 0.0;
        const Complex exact_end = rod_continuum(rod, e_star, rho_star, omega, f_load, g, rod.length);
        std::vector<Scalar> errors;
        std::vector<Scalar> hs;
        json::Value meshes = json::Value::make_array();
        for (const Index n : ladder) {
          const FemModel model = rod_model(rod, type, n, f_load, g);
          const Assembler assembler(model);
          FrequencyResponseOptions options;
          options.frequencies = {freq};
          options.mass_type = mass;
          options.structural_damping = c.damping.eta;
          options.mass_damping = c.damping.a;
          options.stiffness_damping = c.damping.b;
          options.snapshot_frequencies = {freq};
          options.monitors.push_back(
              x_monitor("end_ux", rod.length, DynamicMonitor::Quantity::Displacement));
          options.monitors.push_back(
              x_monitor("root_rx", 0.0, DynamicMonitor::Quantity::Reaction));
          const FrequencyResponseResult r = solve_frequency_response(model, assembler, 0, options);
          Complex reaction;
          const std::vector<Complex> discrete =
              rod_discrete(rod, n, lumped, e_star, rho_star, omega, f_load, g, reaction);
          // Every node's u_x against its station's discrete value.
          const Scalar h = rod.length / static_cast<Scalar>(n);
          Scalar diff = 0.0;
          Scalar scale = 0.0;
          for (const Complex& z : discrete) scale = std::max(scale, std::abs(z));
          for (Index node = 0; node < model.mesh().num_nodes(); ++node) {
            const Index j = static_cast<Index>(std::lround(model.mesh().node(node).x() / h));
            const Complex ours = r.snapshots.front().displacement(model.dofs().dof(node, 0));
            diff = std::max(diff, std::abs(ours - discrete[static_cast<std::size_t>(j)]));
          }
          const Scalar discrete_difference = diff / scale;
          const Scalar reaction_difference =
              std::abs(r.points.front().monitors[1] - reaction) / std::abs(reaction);
          const Complex end = r.points.front().monitors[0];
          const Scalar error = std::abs(end - exact_end) / std::abs(exact_end);
          errors.push_back(error);
          hs.push_back(h);
          const std::size_t last = errors.size() - 1;
          const Scalar order = last > 0 ? observed_order(hs[last - 1], errors[last - 1], hs[last],
                                                         errors[last])
                                        : std::nan("");
          worst_discrete = std::max({worst_discrete, discrete_difference, reaction_difference});
          csv.raw_row({element_name(type), to_string(mass), c.name, fmt(c.ratio, 3),
                       fmt(static_cast<Scalar>(n)), fmt(discrete_difference, 4),
                       fmt(reaction_difference, 4), fmt(error, 6), fmt(order, 4),
                       fmt(std::abs(end), 10), fmt(std::abs(exact_end), 10)});
          json::Value rec = json::Value::make_object();
          rec.set("n", json::Value::make_number(static_cast<Scalar>(n)));
          rec.set("discrete_difference", json::Value::make_number(discrete_difference));
          rec.set("reaction_difference", json::Value::make_number(reaction_difference));
          rec.set("continuum_error", json::Value::make_number(error));
          rec.set("order", json::Value::make_number(order));
          meshes.push_back(rec);
          if (n == ladder.back()) {
            worst_order = std::min(worst_order, order);
            worst_error = std::max(worst_error, error);
          }
        }
        json::Value entry = json::Value::make_object();
        entry.set("element", json::Value::make_string(element_name(type)));
        entry.set("mass", json::Value::make_string(to_string(mass)));
        entry.set("case", json::Value::make_string(c.name));
        entry.set("f_over_f1", json::Value::make_number(c.ratio));
        entry.set("meshes", meshes);
        records.push_back(entry);
      }
    }
  }
  csv.close();

  // A sweep through the first three resonances with structural damping for
  // the figure: the exact continuum response and two meshes.
  {
    CsvWriter sweep(path_join(out_dir, "rod_harmonic_sweep.csv"),
                    {"f_over_f1", "exact_abs[m]", "exact_phase[deg]", "q4_n10_abs[m]",
                     "q4_n10_phase[deg]", "q4_n40_abs[m]", "q4_n40_phase[deg]"});
    std::vector<Scalar> ratios;
    for (int i = 0; i <= 600; ++i) ratios.push_back(0.05 + 5.45 * i / 600.0);
    std::vector<FrequencyResponseResult> runs;
    for (const Index n : {Index(10), Index(40)}) {
      const FemModel model = rod_model(rod, ElementType::Quad4, n, force, 0.0);
      const Assembler assembler(model);
      FrequencyResponseOptions options;
      for (Scalar ratio : ratios) options.frequencies.push_back(ratio * f1);
      options.structural_damping = structural.eta;
      options.monitors.push_back(
          x_monitor("end_ux", rod.length, DynamicMonitor::Quantity::Displacement));
      runs.push_back(solve_frequency_response(model, assembler, 0, options));
    }
    for (std::size_t i = 0; i < ratios.size(); ++i) {
      const Scalar omega = 2.0 * kPi * ratios[i] * f1;
      const Complex exact = rod_continuum(rod, structural.modulus(rod.youngs, omega),
                                          structural.density(rod.density, omega), omega, force,
                                          0.0, rod.length);
      const Complex a = runs[0].points[i].monitors[0];
      const Complex b = runs[1].points[i].monitors[0];
      sweep.row({ratios[i], std::abs(exact), std::arg(exact) * 180.0 / kPi, std::abs(a),
                 std::arg(a) * 180.0 / kPi, std::abs(b), std::arg(b) * 180.0 / kPi});
    }
    sweep.close();
  }

  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string(
                        "verification (exact discrete solution; exact continuum solution with "
                        "its convergence rate)"));
  block.set("f1_Hz", json::Value::make_number(f1));
  block.set("cases", records);
  block.set("note",
            json::Value::make_string(
                "Fixed-free steel rod L = 1 m, 0.05 m square, nu = 0, lateral displacements held "
                "(one element across: an exactly one-dimensional model), f1 = c/(4L). End force "
                "1000 N, or base motion 1e-5 m. Damping: structural eta = 0.02, or Rayleigh 1 % "
                "of critical at f1 from each of a M and b K. Discrete: every node against "
                "u_j = g cos(j theta) + D sin(j theta) with the dispersion relation of the "
                "stencil (cos theta = (6 - 2W)/(6 + W) consistent, 1 - W/2 lumped, W = omega^2 "
                "h^2 rho*/E*), and the root reaction; continuum: the end amplitude against "
                "u = g cos(k(L - x))/cos(kL) + F sin(kx)/(E* A k cos(kL)), k = omega "
                "sqrt(rho*/E*), E* = E(1 + i eta + i omega b), rho* = rho(1 - i a/omega)."));
  summary.set("rod_harmonic", block);

  std::ostringstream note;
  note << "discrete difference " << app::format(worst_discrete, 3)
       << "; continuum order on the finest pair >= " << app::format(worst_order, 3)
       << ", largest error at n = 160 " << app::format(worst_error, 3);
  StudyOutcome outcome;
  outcome.name = "harmonic response of a rod vs exact discrete and continuum solutions";
  outcome.kind = "verification";
  outcome.metric = "largest relative difference to the exact discrete solution";
  outcome.value = worst_discrete;
  outcome.tolerance = 1.0e-9;
  outcome.passed = worst_discrete <= 1.0e-9 && worst_order >= 1.9;
  outcome.note = note.str();
  return outcome;
}

// ---------------------------------------------------------------------------
// Fixed-free rod: transient under a smoothly ramped end force
// ---------------------------------------------------------------------------

namespace {

/// sin(x T) / x and (1 - cos(x T)) / x = 2 sin^2(x T / 2) / x, with their
/// limits at x = 0.
Scalar sin_over(Scalar x, Scalar t) {
  return std::abs(x * t) < 1.0e-8 ? t : std::sin(x * t) / x;
}
Scalar one_minus_cos_over(Scalar x, Scalar t) {
  if (std::abs(x * t) < 1.0e-8) return 0.5 * x * t * t;
  const Scalar s = std::sin(0.5 * x * t);
  return 2.0 * s * s / x;
}

/// The ramp s(t) = sin^2(pi t / (2 T_r)) up to T_r, then 1.
Scalar ramp(Scalar t, Scalar t_ramp) {
  if (t >= t_ramp) return 1.0;
  const Scalar s = std::sin(0.5 * kPi * t / t_ramp);
  return s * s;
}

/// The end displacement of the fixed-free rod under F0 s(t) from rest:
/// modes sin(k_n x), k_n = (2n - 1) pi / (2L), omega_n = c k_n, modal mass
/// rho A L / 2. Integrating Duhamel's integral by parts separates the
/// static response, F0 L s(t) / (E A), from a dynamic remainder
///   -sum_n 2 F0 / (rho A L omega_n^2) int_0^min(t, T_r) s'(tau) cos(omega_n (t - tau)) dtau,
/// whose terms fall as 1/omega_n^4 (s' vanishes at both ends of the ramp).
Scalar rod_ramp_response(const Rod& rod, Scalar force, Scalar t_ramp, Scalar t, int modes) {
  const Scalar c = rod.wave_speed();
  const Scalar beta = 0.5 * kPi / t_ramp;  // s' = beta sin(2 beta tau)
  const Scalar upper = std::min(t, t_ramp);
  const Scalar two_beta = 2.0 * beta;
  Scalar dynamic = 0.0;
  // Sum from the smallest terms up.
  for (int n = modes; n >= 1; --n) {
    const Scalar omega = c * (2.0 * n - 1.0) * kPi / (2.0 * rod.length);
    // int_0^T sin(2 beta tau) cos(omega (t - tau)) dtau
    //   = cos(omega t) I_c + sin(omega t) I_s
    const Scalar i_c = 0.5 * (one_minus_cos_over(two_beta + omega, upper) +
                              one_minus_cos_over(two_beta - omega, upper));
    const Scalar i_s = 0.5 * (sin_over(two_beta - omega, upper) -
                              sin_over(two_beta + omega, upper));
    const Scalar j = beta * (std::cos(omega * t) * i_c + std::sin(omega * t) * i_s);
    dynamic += 2.0 / (rod.density * rod.area() * rod.length * omega * omega) * j;
  }
  return force * (rod.length * ramp(t, t_ramp) / (rod.youngs * rod.area()) - dynamic);
}

}  // namespace

StudyOutcome study_rod_transient(const std::string& out_dir, json::Value& summary) {
  const Rod rod;
  const Scalar c = rod.wave_speed();
  const Scalar period = 4.0 * rod.length / c;  // the first mode's
  const Scalar t_ramp = 0.6 * period;
  const Scalar t_end = 2.5 * period;
  const Scalar force = 1000.0;
  const std::vector<Index> ladder = {20, 40, 80, 160};
  const int series_modes = 4000;
  // Courant number 1: dt = h / c. The ramp is tabulated at every step of the
  // finest run, so every run sees the smooth ramp at its own steps.
  const Scalar dt_fine = rod.length / static_cast<Scalar>(ladder.back()) / c;
  Amplitude table;
  table.kind = Amplitude::Kind::Table;
  const int ramp_steps = static_cast<int>(std::lround(t_ramp / dt_fine));
  for (int k = 0; k <= ramp_steps; ++k) {
    table.times.push_back(k * dt_fine);
    table.values.push_back(ramp(k * dt_fine, t_ramp));
  }

  // The exact end displacement at every step of the finest run; a run with
  // n elements steps every (finest n / n)-th of them.
  const int fine_steps = static_cast<int>(std::lround(t_end / dt_fine));
  std::vector<Scalar> exact_fine;
  for (int i = 0; i <= fine_steps; ++i) {
    exact_fine.push_back(rod_ramp_response(rod, force, t_ramp, i * dt_fine, series_modes));
  }

  CsvWriter csv(path_join(out_dir, "rod_transient.csv"),
                {"element", "mass", "n", "dt[s]", "steps", "max_error[-]", "order[-]",
                 "energy_balance[-]"});
  CsvWriter history(path_join(out_dir, "rod_transient_history.csv"),
                    {"t[s]", "exact[m]", "q4_n20[m]", "q4_n40[m]", "q4_n160[m]"});
  json::Value records = json::Value::make_array();
  Scalar worst_order = 1.0e300;
  Scalar worst_error = 0.0;
  Scalar worst_balance = 0.0;
  std::vector<std::vector<Scalar>> q4_histories;
  std::vector<Scalar> fine_times;
  for (const ElementType type : {ElementType::Quad4, ElementType::Hex8}) {
    for (const MassType mass : {MassType::Consistent, MassType::Lumped}) {
      std::vector<Scalar> errors;
      std::vector<Scalar> hs;
      json::Value meshes = json::Value::make_array();
      for (const Index n : ladder) {
        const FemModel model = rod_model(rod, type, n, force, 0.0);
        const Assembler assembler(model);
        const Scalar h = rod.length / static_cast<Scalar>(n);
        TransientOptions options;
        options.time_step = h / c;
        const int steps = static_cast<int>(std::lround(t_end / options.time_step));
        options.end_time = steps * options.time_step;
        options.mass_type = mass;
        options.amplitude = table;
        options.monitors.push_back(
            x_monitor("end_ux", rod.length, DynamicMonitor::Quantity::Displacement));
        const TransientResult r = solve_transient(model, assembler, 0, options);
        const std::size_t stride = static_cast<std::size_t>(ladder.back() / n);
        Scalar diff = 0.0;
        Scalar scale = 0.0;
        std::vector<Scalar> ends;
        for (const TransientStep& s : r.steps) {
          const Scalar exact = exact_fine[static_cast<std::size_t>(s.index) * stride];
          diff = std::max(diff, std::abs(s.monitors[0] - exact));
          scale = std::max(scale, std::abs(exact));
          ends.push_back(s.monitors[0]);
        }
        const Scalar error = diff / scale;
        errors.push_back(error);
        hs.push_back(h);
        const std::size_t last = errors.size() - 1;
        const Scalar order =
            last > 0 ? observed_order(hs[last - 1], errors[last - 1], hs[last], errors[last])
                     : std::nan("");
        worst_balance = std::max(worst_balance, r.energy_balance_error);
        csv.raw_row({element_name(type), to_string(mass), fmt(static_cast<Scalar>(n)),
                     fmt(options.time_step, 6), fmt(static_cast<Scalar>(steps)), fmt(error, 6),
                     fmt(order, 4), fmt(r.energy_balance_error, 4)});
        json::Value rec = json::Value::make_object();
        rec.set("n", json::Value::make_number(static_cast<Scalar>(n)));
        rec.set("time_step_s", json::Value::make_number(options.time_step));
        rec.set("max_error", json::Value::make_number(error));
        rec.set("order", json::Value::make_number(order));
        rec.set("energy_balance", json::Value::make_number(r.energy_balance_error));
        meshes.push_back(rec);
        if (n == ladder.back()) {
          worst_order = std::min(worst_order, order);
          worst_error = std::max(worst_error, error);
        }
        if (type == ElementType::Quad4 && mass == MassType::Consistent) {
          q4_histories.push_back(ends);
          if (n == ladder.back()) {
            for (const TransientStep& s : r.steps) fine_times.push_back(s.time);
          }
        }
      }
      json::Value entry = json::Value::make_object();
      entry.set("element", json::Value::make_string(element_name(type)));
      entry.set("mass", json::Value::make_string(to_string(mass)));
      entry.set("meshes", meshes);
      records.push_back(entry);
    }
  }
  csv.close();
  // The end displacement over time: the exact series and Q4 (consistent)
  // at n = 20, 40 and 160, on the finest run's steps (the coarser runs'
  // steps are every 8th and 4th of them).
  for (std::size_t i = 0; i < fine_times.size(); ++i) {
    const auto coarse = [&](std::size_t run, std::size_t stride) {
      return i % stride == 0 ? q4_histories[run][i / stride] : std::nan("");
    };
    history.row({fine_times[i], exact_fine[i], coarse(0, 8), coarse(1, 4), q4_histories[3][i]});
  }
  history.close();

  // The modal series' own truncation: halving the number of modes.
  Scalar truncation = 0.0;
  for (int i = 1; i <= 50; ++i) {
    const Scalar t = t_end * i / 50.0;
    truncation = std::max(truncation,
                          std::abs(rod_ramp_response(rod, force, t_ramp, t, series_modes) -
                                   rod_ramp_response(rod, force, t_ramp, t, series_modes / 2)));
  }
  truncation /= force * rod.length / (rod.youngs * rod.area());

  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string(
                        "verification (convergence to the exact continuum solution)"));
  block.set("cases", records);
  block.set("series_modes", json::Value::make_number(series_modes));
  block.set("series_truncation_estimate", json::Value::make_number(truncation));
  block.set("note",
            json::Value::make_string(
                "The rod of rod-harmonic under an end force 1000 N s(t), s = sin^2(pi t / (2 "
                "T_r)) up to T_r = 0.6 T1, then 1, from rest; T1 = 4L/c. Trapezoidal rule "
                "with Courant number 1 (dt = h/c), 2.5 T1. Error: the largest end-displacement "
                "difference over the steps against the exact modal series (4000 modes, the "
                "static part summed exactly), over the largest exact displacement. Truncation "
                "estimate: the series with 4000 against 2000 modes, over F L/(E A)."));
  summary.set("rod_transient", block);

  std::ostringstream note;
  note << "order on the finest pair >= " << app::format(worst_order, 3)
       << ", largest error at n = 160 " << app::format(worst_error, 3) << ", series truncation "
       << app::format(truncation, 2) << ", energy balance <= " << app::format(worst_balance, 2);
  StudyOutcome outcome;
  outcome.name = "transient of a rod under a ramped end force vs exact continuum solution";
  outcome.kind = "verification";
  outcome.metric = "smallest observed convergence order (h and dt halved together)";
  outcome.value = worst_order;
  outcome.tolerance = 1.9;
  outcome.passed = worst_order >= 1.9 && truncation < 1.0e-3 * worst_error;
  outcome.note = note.str();
  return outcome;
}

// ---------------------------------------------------------------------------
// Non-linear oscillator: one element in uniaxial strain
// ---------------------------------------------------------------------------

namespace {

/// One element L x L (x L), its x = 0 face held in x, every lateral
/// displacement held: uniaxial strain, one degree of freedom (the x = L
/// face's u), loaded by `force` shared by that face's nodes.
FemModel oscillator_model(ElementType type, const IsotropicMaterial& material, Scalar length,
                          Scalar force) {
  const bool plane = type == ElementType::Quad4;
  StructuredMeshSpec spec;
  spec.lx = spec.ly = spec.lz = length;
  FemModel model(plane ? make_structured_quad_mesh(spec) : make_structured_hex_mesh(spec),
                 material, plane ? length : 1.0,
                 plane ? StressState::PlaneStrain : StressState::ThreeDimensional,
                 IntegrationOptions());
  DisplacementConstraint lateral;
  lateral.region.name = "lateral";
  lateral.region.members.push_back(Selector());
  lateral.fix_y = true;
  lateral.fix_z = !plane;
  model.constraints().push_back(lateral);
  DisplacementConstraint root;
  root.region.name = "root";
  root.region.members.push_back(x_box(-std::numeric_limits<Scalar>::infinity(), 0.0));
  root.fix_x = true;
  model.constraints().push_back(root);
  LoadCaseSpec lc;
  lc.name = "pull";
  PointLoadSpec end;
  end.region.name = "end";
  end.region.members.push_back(x_box(length, std::numeric_limits<Scalar>::infinity()));
  end.force = Vector3(force, 0.0, 0.0);
  end.distribute_total = true;
  lc.point_loads.push_back(end);
  model.load_case_specs().push_back(lc);
  model.finalize();
  return model;
}

/// The axial force of the element at end displacement u and its slope, from
/// the state before the step (plastic strain and accumulated plastic strain
/// of the (1, -1/2, -1/2) uniaxial pattern), which `commit` advances.
struct UniaxialLaw {
  Scalar area = 0.0;
  Scalar length = 0.0;
  Scalar lambda = 0.0;
  Scalar mu = 0.0;
  bool finite = false;   ///< Saint Venant-Kirchhoff at large strain, else small strain
  Scalar yield = 0.0;    ///< > 0: J2 with linear isotropic hardening
  Scalar hardening = 0.0;
  Scalar plastic = 0.0;  ///< committed epsilon_p (11 component)
  Scalar accumulated = 0.0;

  /// Force and slope at u from the committed state; `trial` returns the new
  /// plastic state.
  void evaluate(Scalar u, Scalar& force, Scalar& slope, Scalar* new_plastic = nullptr,
                Scalar* new_accumulated = nullptr) const {
    const Scalar m = lambda + 2.0 * mu;
    if (finite) {
      // F11 = 1 + u/L, E11 = (F11^2 - 1)/2, S11 = (lambda + 2 mu) E11,
      // nominal stress P11 = F11 S11.
      const Scalar f11 = 1.0 + u / length;
      const Scalar e11 = 0.5 * (f11 * f11 - 1.0);
      force = area * f11 * m * e11;
      slope = area / length * m * (e11 + f11 * f11);
      return;
    }
    const Scalar eps = u / length;
    const Scalar bulk = lambda + 2.0 * mu / 3.0;
    Scalar s11 = 2.0 * mu * (2.0 * eps / 3.0 - plastic);
    Scalar ds11 = 4.0 * mu / 3.0;
    Scalar p_new = plastic;
    Scalar a_new = accumulated;
    if (yield > 0.0) {
      const Scalar trial = 1.5 * std::abs(s11);
      const Scalar limit = yield + hardening * accumulated;
      if (trial > limit) {
        const Scalar sign = s11 > 0.0 ? 1.0 : -1.0;
        const Scalar gamma = (trial - limit) / (3.0 * mu + hardening);
        p_new = plastic + sign * gamma;
        a_new = accumulated + gamma;
        s11 -= 2.0 * mu * sign * gamma;
        ds11 = 4.0 * mu / 3.0 * hardening / (3.0 * mu + hardening);
      }
    }
    if (new_plastic) *new_plastic = p_new;
    if (new_accumulated) *new_accumulated = a_new;
    force = area * (bulk * eps + s11);
    slope = area / length * (bulk + ds11);
  }
  void commit(Scalar u) {
    Scalar f, k, p, a;
    evaluate(u, f, k, &p, &a);
    plastic = p;
    accumulated = a;
  }
};

/// The HHT-alpha method on m u'' + c u' + N(u) = F (constant), from rest:
/// Newton's method to machine precision on
/// m a1 + (1 + alpha)(c v1 + N(u1)) - alpha (c v0 + N(u0)) = F.
std::vector<Scalar> scalar_hht(UniaxialLaw law, Scalar mass, Scalar damping, Scalar force,
                               Scalar dt, int steps, Scalar alpha) {
  const Scalar beta = 0.25 * (1.0 - alpha) * (1.0 - alpha);
  const Scalar gamma = 0.5 - alpha;
  const Scalar c0 = 1.0 / (beta * dt * dt);
  const Scalar c1 = 1.0 / (beta * dt);
  const Scalar c2 = 1.0 / (2.0 * beta) - 1.0;
  const Scalar c3 = gamma / (beta * dt);
  const Scalar c4 = 1.0 - gamma / beta;
  const Scalar c5 = dt * (1.0 - gamma / (2.0 * beta));
  const Scalar w = 1.0 + alpha;
  Scalar u = 0.0;
  Scalar v = 0.0;
  Scalar n0 = 0.0;
  Scalar k0 = 0.0;
  law.evaluate(0.0, n0, k0);
  Scalar a = (force - n0) / mass;
  std::vector<Scalar> out{u};
  for (int s = 1; s <= steps; ++s) {
    Scalar u1 = u + dt * v + 0.5 * dt * dt * a;
    for (int it = 0; it < 60; ++it) {
      Scalar n1 = 0.0;
      Scalar k1 = 0.0;
      law.evaluate(u1, n1, k1);
      const Scalar a1 = c0 * (u1 - u) - c1 * v - c2 * a;
      const Scalar v1 = c3 * (u1 - u) + c4 * v + c5 * a;
      const Scalar r = mass * a1 + w * (damping * v1 + n1) - alpha * (damping * v + n0) - force;
      const Scalar du = -r / (c0 * mass + w * (c3 * damping + k1));
      u1 += du;
      if (std::abs(du) <= 4.0 * std::numeric_limits<Scalar>::epsilon() * std::abs(u1)) break;
    }
    law.commit(u1);
    Scalar n1 = 0.0;
    Scalar k1 = 0.0;
    law.evaluate(u1, n1, k1);
    const Scalar a1 = c0 * (u1 - u) - c1 * v - c2 * a;
    v = c3 * (u1 - u) + c4 * v + c5 * a;
    a = a1;
    u = u1;
    n0 = n1;
    out.push_back(u);
  }
  return out;
}

/// The classical Runge-Kutta method on the elastic oscillator (no history):
/// u at every multiple of `every` of its step h.
std::vector<Scalar> runge_kutta(const UniaxialLaw& law, Scalar mass, Scalar force, Scalar h,
                                int steps, int every) {
  const auto accel = [&](Scalar u) {
    Scalar n = 0.0;
    Scalar k = 0.0;
    law.evaluate(u, n, k);
    return (force - n) / mass;
  };
  Scalar u = 0.0;
  Scalar v = 0.0;
  std::vector<Scalar> out{u};
  for (int s = 1; s <= steps; ++s) {
    const Scalar k1u = v;
    const Scalar k1v = accel(u);
    const Scalar k2u = v + 0.5 * h * k1v;
    const Scalar k2v = accel(u + 0.5 * h * k1u);
    const Scalar k3u = v + 0.5 * h * k2v;
    const Scalar k3v = accel(u + 0.5 * h * k2u);
    const Scalar k4u = v + h * k3v;
    const Scalar k4v = accel(u + h * k3u);
    u += h / 6.0 * (k1u + 2.0 * k2u + 2.0 * k3u + k4u);
    v += h / 6.0 * (k1v + 2.0 * k2v + 2.0 * k3v + k4v);
    if (s % every == 0) out.push_back(u);
  }
  return out;
}

/// The exact motion of the elastoplastic oscillator (linear isotropic
/// hardening, uniaxial strain) under a sudden force from rest, piecewise in
/// closed form: elastic until u reaches u_y, plastic (slope k_p) until the
/// velocity vanishes at u_m, then elastic about the shifted equilibrium.
/// Valid while the unloading stays elastic, which the caller checks.
struct PlasticMotion {
  Scalar mass, force, k_e, k_p, u_y, n_y;
  Scalar t_y, v_y, u_star, radius, phase, t_m, u_m, n_m;
  PlasticMotion(Scalar m, Scalar f, Scalar ke, Scalar kp, Scalar uy)
      : mass(m), force(f), k_e(ke), k_p(kp), u_y(uy), n_y(ke * uy) {
    const Scalar we = std::sqrt(k_e / mass);
    const Scalar wp = std::sqrt(k_p / mass);
    t_y = std::acos(1.0 - k_e * u_y / force) / we;
    v_y = force / k_e * we * std::sin(we * t_y);
    u_star = u_y + (force - n_y) / k_p;
    radius = std::hypot(u_y - u_star, v_y / wp);
    phase = std::atan2(v_y / wp, u_y - u_star);
    t_m = t_y + phase / wp;
    u_m = u_star + radius;
    n_m = n_y + k_p * (u_m - u_y);
  }
  Scalar operator()(Scalar t) const {
    const Scalar we = std::sqrt(k_e / mass);
    const Scalar wp = std::sqrt(k_p / mass);
    if (t <= t_y) return force / k_e * (1.0 - std::cos(we * t));
    if (t <= t_m) return u_star + radius * std::cos(wp * (t - t_y) - phase);
    const Scalar centre = u_m - (n_m - force) / k_e;
    return centre + (u_m - centre) * std::cos(we * (t - t_m));
  }
};

}  // namespace

StudyOutcome study_nonlinear_oscillator(const std::string& out_dir, json::Value& summary) {
  const Scalar length = 0.1;
  const Scalar youngs = 200.0e9;
  const Scalar poisson = 0.3;
  const Scalar rho = 7850.0;
  const Scalar lambda = youngs * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
  const Scalar mu = youngs / (2.0 * (1.0 + poisson));
  const Scalar area = length * length;
  const Scalar k_e = area * (lambda + 2.0 * mu) / length;

  CsvWriter csv(path_join(out_dir, "nonlinear_oscillator.csv"),
                {"law", "element", "mass", "alpha", "steps_per_period", "discrete_difference[-]",
                 "error[-]", "order[-]", "iterations", "final_balance[J]",
                 "exact_dissipation[J]"});
  CsvWriter history(path_join(out_dir, "nonlinear_oscillator_history.csv"),
                    {"law", "t[s]", "exact[m]", "coarse[m]", "fine[m]"});
  json::Value records = json::Value::make_array();
  Scalar worst_discrete = 0.0;
  Scalar worst_order = 1.0e300;
  Scalar reference_check = std::nan("");  // Runge-Kutta at h against h/2
  std::ostringstream note;

  struct Law {
    std::string name;
    bool finite;
    Scalar yield;
    Scalar hardening;
    Scalar force;  ///< the sudden load [N]
  };
  // Saint Venant-Kirchhoff pulled to a peak strain near 10 %; the plastic
  // law yields on the first swing (the load is 0.8 of the yield force, the
  // sudden application nearly doubles it) and unloads elastically.
  const Scalar sy = 250.0e6;
  const Scalar u_yield = length * sy / (2.0 * mu);
  const std::vector<Law> laws = {{"Saint Venant-Kirchhoff finite strain", true, 0.0, 0.0,
                                  0.05 * area * (lambda + 2.0 * mu)},
                                 {"J2 linear hardening small strain", false, sy, 20.0e9,
                                  0.8 * k_e * u_yield}};
  for (const Law& law_spec : laws) {
    IsotropicMaterial material(youngs, poisson, rho, "steel");
    if (law_spec.yield > 0.0) {
      PlasticityParameters p;
      p.yield_stress = law_spec.yield;
      p.hardening_modulus = law_spec.hardening;
      material.set_plasticity(p);
    }
    UniaxialLaw law;
    law.area = area;
    law.length = length;
    law.lambda = lambda;
    law.mu = mu;
    law.finite = law_spec.finite;
    law.yield = law_spec.yield;
    law.hardening = law_spec.hardening;
    for (const ElementType type : {ElementType::Quad4, ElementType::Hex8}) {
      const FemModel model = oscillator_model(type, material, length, law_spec.force);
      const Assembler assembler(model);
      for (const MassType mass_type : {MassType::Lumped, MassType::Consistent}) {
        // The oscillator's mass: rho A L / 2 lumped, rho A L / 3 consistent
        // (the end face's share of the element's consistent mass).
        const Scalar mass = rho * area * length / (mass_type == MassType::Lumped ? 2.0 : 3.0);
        const Scalar period = 2.0 * kPi * std::sqrt(mass / k_e);
        for (const Scalar alpha : {0.0, -0.1}) {
          // The exact motion over 1.5 periods, and its plastic dissipation.
          const Scalar t_end = 1.5 * period;
          std::function<Scalar(Scalar)> exact;
          Scalar exact_dissipation = 0.0;
          if (law_spec.yield > 0.0) {
          const Scalar bulk = lambda + 2.0 * mu / 3.0;
            const Scalar k_p = area / length *
                               (bulk + 4.0 * mu / 3.0 * law_spec.hardening /
                                           (3.0 * mu + law_spec.hardening));
            const PlasticMotion motion(mass, law_spec.force, k_e, k_p, u_yield);
            // The unloading must stay elastic: the deviatoric stress falls by
            // at most (4/3) mu times the strain range of the swing.
            const Scalar alpha_m =
                (2.0 * mu * motion.u_m / length - sy) / (3.0 * mu + law_spec.hardening);
            const Scalar swing = 2.0 * (motion.n_m - law_spec.force) / k_e / length;
            if (!(4.0 / 3.0 * mu * swing < 4.0 / 3.0 * (sy + law_spec.hardening * alpha_m))) {
              throw SolverError("nonlinear-oscillator: the unloading of the plastic case would "
                                "yield again; its closed form does not apply");
            }
            if (!(motion.t_m < t_end)) {
              throw SolverError("nonlinear-oscillator: the plastic swing outlasts the run");
            }
            exact = motion;
            // sigma_y times the accumulated plastic strain over the volume:
            // the plastic work less the stored hardening energy.
            exact_dissipation = sy * alpha_m * area * length;
          }
          // The step ladder: 20 to 320 steps per period (a whole number of
          // steps in 1.5 periods).
          std::vector<Scalar> errors;
          std::vector<Scalar> dts;
          std::vector<std::vector<Scalar>> runs;
          std::vector<Scalar> run_times;
          for (const int per_period : {20, 40, 80, 160, 320}) {
            TransientOptions options;
            options.time_step = period / per_period;
            const int steps = per_period * 3 / 2;
            options.end_time = steps * options.time_step;
            options.mass_type = mass_type;
            options.alpha = alpha;
            options.nonlinear = true;
            options.nonlinear_options.kinematics =
                law_spec.finite ? Kinematics::Finite : Kinematics::SmallStrain;
            options.nonlinear_options.residual_tolerance = 1.0e-12;
            options.nonlinear_options.displacement_tolerance = 1.0e-12;
            options.monitors.push_back(
                x_monitor("end_ux", length, DynamicMonitor::Quantity::Displacement));
            const TransientResult r = solve_transient(model, assembler, 0, options);
            if (!r.completed) throw SolverError("nonlinear-oscillator: " + r.termination);
            const std::vector<Scalar> reference =
                scalar_hht(law, mass, 0.0, law_spec.force, options.time_step, steps, alpha);
            Scalar diff = 0.0;
            Scalar scale = 0.0;
            Scalar error = 0.0;
            std::vector<Scalar> exact_values;
            if (!exact) {
              // Runge-Kutta at 1/64 of the finest step of the ladder.
              const int every = 64 * 320 / per_period;
              exact_values = runge_kutta(law, mass, law_spec.force, period / (320.0 * 64.0),
                                         steps * every, every);
            }
            std::vector<Scalar> ends;
            for (const TransientStep& s : r.steps) {
              const std::size_t i = static_cast<std::size_t>(s.index);
              diff = std::max(diff, std::abs(s.monitors[0] - reference[i]));
              scale = std::max(scale, std::abs(reference[i]));
              const Scalar truth = exact ? exact(s.time) : exact_values[i];
              error = std::max(error, std::abs(s.monitors[0] - truth));
              ends.push_back(s.monitors[0]);
            }
            const Scalar discrete_difference = diff / scale;
            error /= scale;
            errors.push_back(error);
            dts.push_back(options.time_step);
            const std::size_t last = errors.size() - 1;
            const Scalar order = last > 0 ? observed_order(dts[last - 1], errors[last - 1],
                                                           dts[last], errors[last])
                                          : std::nan("");
            worst_discrete = std::max(worst_discrete, discrete_difference);
            csv.raw_row({law_spec.name, element_name(type), to_string(mass_type), fmt(alpha, 3),
                         fmt(static_cast<Scalar>(per_period)), fmt(discrete_difference, 4),
                         fmt(error, 6), fmt(order, 4),
                         fmt(static_cast<Scalar>(r.total_iterations)),
                         fmt(r.numerical_dissipation, 8), fmt(exact_dissipation, 8)});
            json::Value rec = json::Value::make_object();
            rec.set("law", json::Value::make_string(law_spec.name));
            rec.set("element", json::Value::make_string(element_name(type)));
            rec.set("mass", json::Value::make_string(to_string(mass_type)));
            rec.set("alpha", json::Value::make_number(alpha));
            rec.set("steps_per_period", json::Value::make_number(per_period));
            rec.set("discrete_difference", json::Value::make_number(discrete_difference));
            rec.set("error", json::Value::make_number(error));
            rec.set("order", json::Value::make_number(order));
            rec.set("final_balance_J", json::Value::make_number(r.numerical_dissipation));
            if (law_spec.yield > 0.0) {
              rec.set("exact_dissipation_J", json::Value::make_number(exact_dissipation));
            }
            records.push_back(rec);
            if (per_period == 320) worst_order = std::min(worst_order, order);
            if (type == ElementType::Quad4 && mass_type == MassType::Lumped && alpha == 0.0) {
              runs.push_back(ends);
              if (per_period == 320) {
                for (const TransientStep& s : r.steps) run_times.push_back(s.time);
              }
            }
          }
          if (!exact && type == ElementType::Quad4 && mass_type == MassType::Lumped &&
              alpha == 0.0) {
            // The Runge-Kutta reference against itself at half its step.
            const int steps = 3 * 320 / 2;
            const std::vector<Scalar> once =
                runge_kutta(law, mass, law_spec.force, period / (320.0 * 64.0), steps * 64, 64);
            const std::vector<Scalar> halved = runge_kutta(
                law, mass, law_spec.force, period / (320.0 * 128.0), steps * 128, 128);
            Scalar top = 0.0;
            Scalar gap = 0.0;
            for (std::size_t i = 0; i < once.size(); ++i) {
              gap = std::max(gap, std::abs(once[i] - halved[i]));
              top = std::max(top, std::abs(halved[i]));
            }
            reference_check = gap / top;
          }
          if (!runs.empty()) {
            // The history for the figure: 20 and 320 steps per period, and
            // the exact motion on the finest steps.
            std::vector<Scalar> rk;
            if (!exact) {
              rk = runge_kutta(law, mass, law_spec.force, period / (320.0 * 64.0),
                               static_cast<int>(run_times.size() - 1) * 64, 64);
            }
            for (std::size_t i = 0; i < run_times.size(); ++i) {
              const Scalar truth = exact ? exact(run_times[i]) : rk[i];
              history.raw_row({law_spec.name, fmt(run_times[i], 10), fmt(truth, 12),
                               i % 16 == 0 ? fmt(runs.front()[i / 16], 12) : std::string(),
                               fmt(runs.back()[i], 12)});
            }
          }
        }
      }
    }
    note << law_spec.name << "; ";
  }
  csv.close();
  history.close();

  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string(
                        "verification (exact discrete recursion; exact motion with its "
                        "convergence rate)"));
  block.set("cases", records);
  block.set("runge_kutta_self_difference", json::Value::make_number(reference_check));
  block.set("note",
            json::Value::make_string(
                "One element 0.1 m (Q4 in plane strain 0.1 m thick, Hex8), steel E = 200 GPa, "
                "nu = 0.3, rho = 7850 kg/m^3, its x = 0 face held in x and every lateral "
                "displacement held: uniaxial strain, a single degree of freedom with mass rho A "
                "L/2 (lumped) or rho A L/3 (consistent). Sudden end force from rest. Saint "
                "Venant-Kirchhoff at finite strain: N(u) = A (lambda + 2 mu) F11 (F11^2 - "
                "1)/2, F11 = 1 + u/L, 0.05 A (lambda + 2 mu); J2 small strain with sigma_y = "
                "250 MPa and linear isotropic hardening 20 GPa, 0.8 of the yield force. "
                "Discrete: the scalar HHT-alpha recursion of m u'' + N(u) = F with Newton to "
                "machine precision (the plastic law's radial return in closed form). Exact: "
                "classical Runge-Kutta at 1/64 of the finest step (elastic; its difference "
                "from itself at half the step is runge_kutta_self_difference), the piecewise "
                "closed form (plastic: elastic to u_y, slope K + (4/3) mu H/(3 mu + H) to the "
                "turning point, elastic after). Errors over the largest displacement; "
                "exact_dissipation = sigma_y alpha_p A L, against the run's final energy "
                "balance."));
  summary.set("nonlinear_oscillator", block);

  StudyOutcome outcome;
  outcome.name = "non-linear oscillator (finite-strain elastic, elastoplastic) vs exact motion";
  outcome.kind = "verification";
  outcome.metric = "largest relative difference to the scalar HHT-alpha recursion";
  outcome.value = worst_discrete;
  outcome.tolerance = 1.0e-9;
  outcome.passed = worst_discrete <= 1.0e-9 && worst_order >= 1.8 && reference_check < 1.0e-10;
  std::ostringstream full;
  full << note.str() << "smallest order on the finest pair " << app::format(worst_order, 3)
       << ", Runge-Kutta reference self-difference " << app::format(reference_check, 2);
  outcome.note = full.str();
  return outcome;
}

// ---------------------------------------------------------------------------
// Explicit central differences on the rod (ExplicitDynamics.hpp)
// ---------------------------------------------------------------------------

namespace {

/// The drive of an explicit run of a rod model: its constraints held at
/// zero, the load case's end force scaled by `ramp`.
ExplicitDrive rod_drive(const FemModel& model, Scalar duration, const Amplitude& ramp) {
  ExplicitDrive drive;
  drive.duration = duration;
  drive.fixed = model.dofs().constrained_dofs();
  drive.fixed_start = Vector::Zero(static_cast<Eigen::Index>(drive.fixed.size()));
  drive.fixed_end = drive.fixed_start;
  drive.ramp = ramp;
  return drive;
}

/// The x displacement of the rod's end (the mean over the x = L nodes).
Scalar rod_end(const FemModel& model, const Vector& u, Scalar length) {
  const Mesh& mesh = model.mesh();
  const int dim = mesh.dim();
  Scalar sum = 0.0;
  int count = 0;
  for (Index n = 0; n < mesh.num_nodes(); ++n) {
    if (mesh.node(n).x() < length - 1.0e-9) continue;
    sum += u(n * dim);
    ++count;
  }
  return sum / count;
}

}  // namespace

StudyOutcome study_explicit_rod(const std::string& out_dir, json::Value& summary) {
  const Rod rod;
  const Scalar c = rod.wave_speed();
  const Scalar period = 4.0 * rod.length / c;  // the first mode's
  const Scalar t_ramp = 0.6 * period;
  const Scalar t_end = 2.5 * period;
  const Scalar force = 1000.0;
  const Scalar courant = 0.5;
  const std::vector<Index> ladder = {20, 40, 80, 160};
  const int series_modes = 4000;
  const Scalar static_end = force * rod.length / (rod.youngs * rod.area());
  NonlinearOptions linear;
  linear.kinematics = Kinematics::SmallStrain;

  // (1) The transient of rod-transient, integrated explicitly: lumped mass,
  // Courant number 0.5, the ramp tabulated at the finest run's steps (the
  // coarser runs' steps are among them).
  const Scalar dt_fine = courant * rod.length / static_cast<Scalar>(ladder.back()) / c;
  Amplitude table;
  table.kind = Amplitude::Kind::Table;
  const int ramp_steps = static_cast<int>(std::lround(t_ramp / dt_fine));
  for (int k = 0; k <= ramp_steps; ++k) {
    table.times.push_back(k * dt_fine);
    table.values.push_back(ramp(k * dt_fine, t_ramp));
  }
  table.times.push_back(t_end * 10.0);
  table.values.push_back(1.0);
  CsvWriter csv(path_join(out_dir, "explicit_rod.csv"),
                {"element", "kernel", "n", "dt[s]", "steps", "max_error[-]", "order[-]",
                 "energy_balance[-]"});
  json::Value records = json::Value::make_array();
  Scalar worst_order = 1.0e300;
  Scalar worst_error = 0.0;
  Scalar worst_balance = 0.0;
  bool dedicated_used = false;
  for (const ElementType type : {ElementType::Quad4, ElementType::Hex8}) {
    std::vector<Scalar> errors;
    std::vector<Scalar> hs;
    json::Value meshes = json::Value::make_array();
    std::string kernel;
    for (const Index n : ladder) {
      const FemModel model = rod_model(rod, type, n, force, 0.0);
      const Assembler assembler(model);
      const Scalar h = rod.length / static_cast<Scalar>(n);
      ExplicitOptions options;
      options.time_step = courant * h / c;
      const long steps = std::lround(t_end / options.time_step);
      const auto stride = static_cast<int>(ladder.back() / n);
      options.history_every = 1000000;
      // The end displacement at every 8th step of the finest run.
      options.snapshot_every = 8 / stride;
      ExplicitDynamics integrator(model, assembler, linear, {}, options, 0);
      const ExplicitResult r =
          integrator.run(rod_drive(model, static_cast<Scalar>(steps) * options.time_step, table),
                         ExplicitState{Vector::Zero(model.dofs().num_dofs()), {}, {}, {}});
      if (!r.completed) throw SolverError("explicit-rod: " + r.termination);
      kernel = r.kernel;
      dedicated_used = dedicated_used || r.kernel == "dedicated Hex8";
      Scalar diff = 0.0;
      Scalar scale = 0.0;
      for (const ExplicitSnapshot& snap : r.snapshots) {
        const Scalar exact = rod_ramp_response(rod, force, t_ramp, snap.time, series_modes);
        diff = std::max(diff, std::abs(rod_end(model, snap.displacement, rod.length) - exact));
        scale = std::max(scale, std::abs(exact));
      }
      const Scalar error = diff / scale;
      errors.push_back(error);
      hs.push_back(h);
      const std::size_t last = errors.size() - 1;
      const Scalar order =
          last > 0 ? observed_order(hs[last - 1], errors[last - 1], hs[last], errors[last])
                   : std::nan("");
      worst_balance = std::max(worst_balance, r.max_energy_error);
      csv.raw_row({element_name(type), r.kernel, fmt(static_cast<Scalar>(n)),
                   fmt(options.time_step, 6), fmt(static_cast<Scalar>(r.steps)), fmt(error, 6),
                   fmt(order, 4), fmt(r.max_energy_error, 4)});
      json::Value rec = json::Value::make_object();
      rec.set("n", json::Value::make_number(static_cast<Scalar>(n)));
      rec.set("time_step_s", json::Value::make_number(options.time_step));
      rec.set("steps", json::Value::make_number(static_cast<Scalar>(r.steps)));
      rec.set("max_error", json::Value::make_number(error));
      rec.set("order", json::Value::make_number(order));
      rec.set("energy_balance", json::Value::make_number(r.max_energy_error));
      meshes.push_back(rec);
      if (n == ladder.back()) {
        worst_order = std::min(worst_order, order);
        worst_error = std::max(worst_error, error);
      }
    }
    json::Value entry = json::Value::make_object();
    entry.set("element", json::Value::make_string(element_name(type)));
    entry.set("kernel", json::Value::make_string(kernel));
    entry.set("meshes", meshes);
    records.push_back(entry);
  }
  csv.close();

  // (2) The stable step: the element eigenvalue bound and power iteration
  // against 2 / omega_max of the constrained lumped system (dense), n = 20.
  json::Value stability = json::Value::make_array();
  bool bounds_hold = true;
  Scalar worst_bound_ratio = 0.0;
  for (const ElementType type : {ElementType::Quad4, ElementType::Hex8}) {
    const FemModel model = rod_model(rod, type, 20, force, 0.0);
    const Assembler assembler(model);
    const Matrix k(assembler.reduce_free_free(assembler.assemble_stiffness()));
    const Matrix m(assembler.reduce_free_free(assembler.assemble_mass(MassType::Lumped)));
    Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> es(k, m, Eigen::EigenvaluesOnly);
    const Scalar exact = 2.0 / std::sqrt(es.eigenvalues().maxCoeff());
    ExplicitOptions eo;
    const ExplicitDynamics by_element(model, assembler, linear, {}, eo, 0);
    const Vector zero = Vector::Zero(model.dofs().num_dofs());
    const Scalar element = by_element.element_time_steps(zero).minCoeff();
    eo.stable_step.method = StableStepOptions::Method::PowerIteration;
    eo.stable_step.power_iterations = 400;
    const ExplicitDynamics by_power(model, assembler, linear, {}, eo, 0);
    const Scalar power = by_power.power_iteration_time_step(zero);
    bounds_hold = bounds_hold && element <= exact * (1.0 + 1.0e-12) &&
                  power <= exact * (1.0 + 1.0e-12);
    worst_bound_ratio = std::max(worst_bound_ratio, exact / element);
    json::Value rec = json::Value::make_object();
    rec.set("element", json::Value::make_string(element_name(type)));
    rec.set("exact_constrained_s", json::Value::make_number(exact));
    rec.set("element_eigenvalue_s", json::Value::make_number(element));
    rec.set("power_iteration_s", json::Value::make_number(power));
    rec.set("one_dimensional_s", json::Value::make_number(rod.length / 20.0 / c));
    stability.push_back(rec);
  }

  // (3) Quasi-static: the end force ramped over 20 periods and held for 5
  // with mass-proportional damping (20 % of critical in the first mode)
  // reaches the static displacement F L / (E A).
  Scalar static_error = 0.0;
  Scalar static_kinetic = 0.0;
  {
    const FemModel model = rod_model(rod, ElementType::Hex8, 40, force, 0.0);
    const Assembler assembler(model);
    ExplicitOptions options;
    options.mass_damping = 2.0 * 0.2 * 2.0 * kPi / period;
    options.history_every = 1000;
    ExplicitDynamics integrator(model, assembler, linear, {}, options, 0);
    Amplitude slow;
    slow.kind = Amplitude::Kind::Table;
    const int knots = 400;
    for (int k = 0; k <= knots; ++k) {
      const Scalar t = 20.0 * period * k / knots;
      slow.times.push_back(t);
      slow.values.push_back(ramp(t, 20.0 * period));
    }
    slow.times.push_back(100.0 * period);
    slow.values.push_back(1.0);
    const ExplicitResult r =
        integrator.run(rod_drive(model, 25.0 * period, slow),
                       ExplicitState{Vector::Zero(model.dofs().num_dofs()), {}, {}, {}});
    if (!r.completed) throw SolverError("explicit-rod (quasi-static): " + r.termination);
    static_error =
        std::abs(rod_end(model, r.final_state.displacement, rod.length) - static_end) / static_end;
    static_kinetic = r.records.back().kinetic / r.records.back().internal;
  }

  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string(
                        "verification (convergence to the exact continuum solution; stable "
                        "step bounds; quasi-static limit)"));
  block.set("cases", records);
  block.set("stable_step", stability);
  block.set("quasi_static_error", json::Value::make_number(static_error));
  block.set("quasi_static_kinetic_ratio", json::Value::make_number(static_kinetic));
  block.set("note",
            json::Value::make_string(
                "The rod of rod-transient (1 m, 0.05 m square, E = 200 GPa, nu = 0, lateral "
                "displacements held) under the end force 1000 N s(t), s = sin^2(pi t / (2 T_r)) "
                "up to T_r = 0.6 T1, from rest, 2.5 T1; explicit central differences with the "
                "lumped mass at Courant number 0.5 (dt = h / (2c)), small strain; Q4 through the "
                "generic element dispatch, Hex8 through the dedicated kernel. Error: the largest "
                "end-displacement difference at every 8th step of the finest run against the "
                "exact modal series (4000 modes), over the largest exact displacement. Stable "
                "step: the element eigenvalue bound and power iteration (unconstrained) "
                "against 2 / omega_max of the constrained lumped system; both must lie below "
                "it. Quasi-static: Hex8, n = 40, the force ramped over 20 T1 and held 5 T1 with "
                "mass damping 0.4 omega_1, against F L / (E A)."));
  summary.set("explicit_rod", block);

  std::ostringstream note;
  note << "order on the finest pair >= " << app::format(worst_order, 3)
       << ", largest error at n = 160 " << app::format(worst_error, 3)
       << ", energy balance <= " << app::format(worst_balance, 2)
       << "; stable-step bounds hold: " << (bounds_hold ? "yes" : "NO")
       << " (exact / element bound <= " << app::format(worst_bound_ratio, 3)
       << "); quasi-static end displacement error " << app::format(static_error, 2)
       << " (T / W_int " << app::format(static_kinetic, 2) << ")";
  StudyOutcome outcome;
  outcome.name = "explicit central differences on a rod vs exact continuum solution";
  outcome.kind = "verification";
  outcome.metric = "smallest observed convergence order (h and dt halved together)";
  outcome.value = worst_order;
  outcome.tolerance = 1.9;
  outcome.passed = worst_order >= 1.9 && worst_balance < 1.0e-2 && bounds_hold &&
                   static_error < 1.0e-3 && dedicated_used;
  outcome.note = note.str();
  return outcome;
}

}  // namespace verify
}  // namespace sparlab
