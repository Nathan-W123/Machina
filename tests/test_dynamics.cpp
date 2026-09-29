/// \file test_dynamics.cpp
/// \brief Transient dynamics (HHT-alpha) and the harmonic response, exact
///        against the modal decomposition of the same discrete model.
///
/// A linear model decouples in its M-orthonormal modes (Rayleigh damping is
/// diagonal in them too), and the HHT-alpha method - linear in (u, v, a) -
/// applied to the whole system is the same method applied to every mode. The
/// tests integrate each mode with the scalar recursion written in the
/// acceleration (predictor-corrector) form, not the displacement form the
/// solver uses, and compare the global solution step by step. The
/// trapezoidal rule's exact energy balance, the momentum balance with the
/// reactions of a prescribed motion, and the modal sum of the frequency
/// response are checked the same way.
#include "TestSupport.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/fem/Dynamics.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <Eigen/Eigenvalues>

#include <cmath>
#include <complex>
#include <vector>

using namespace sparlab;
using namespace sparlab::testing;
using Catch::Approx;

namespace {

constexpr Scalar kPi = 3.14159265358979323846;

struct ModalBasis {
  Vector omega2;  ///< ascending [1/s^2]
  Matrix phi;     ///< free DOFs x modes, M-orthonormal
};

ModalBasis modal_basis(const Assembler& assembler, MassType type) {
  const Matrix k(assembler.reduce_free_free(assembler.assemble_stiffness()));
  const Matrix m(assembler.reduce_free_free(assembler.assemble_mass(type)));
  Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> es(k, m);
  REQUIRE(es.info() == Eigen::Success);
  return {es.eigenvalues(), es.eigenvectors()};
}

/// The HHT-alpha method on one mode, q'' + c q' + w2 q = p(t_n), in the
/// acceleration form: the Newmark predictors from the last step, the
/// equation of motion at the alpha-weighted point solved for the new
/// acceleration, then the correctors.
std::vector<Scalar> integrate_mode(Scalar w2, Scalar c, const std::vector<Scalar>& p, Scalar dt,
                                   Scalar alpha, Scalar q0, Scalar v0) {
  const Scalar beta = 0.25 * (1.0 - alpha) * (1.0 - alpha);
  const Scalar gamma = 0.5 - alpha;
  const Scalar w = 1.0 + alpha;
  Scalar q = q0;
  Scalar v = v0;
  Scalar a = p[0] - c * v0 - w2 * q0;
  std::vector<Scalar> out{q};
  for (std::size_t n = 1; n < p.size(); ++n) {
    const Scalar q_pred = q + dt * v + dt * dt * (0.5 - beta) * a;
    const Scalar v_pred = v + dt * (1.0 - gamma) * a;
    // a1 + w (c v1 + w2 q1) - alpha (c v + w2 q) = w p1 - alpha p0.
    const Scalar a1 = (w * p[n] - alpha * p[n - 1] - w * (c * v_pred + w2 * q_pred) +
                       alpha * (c * v + w2 * q)) /
                      (1.0 + w * c * gamma * dt + w * w2 * beta * dt * dt);
    q = q_pred + beta * dt * dt * a1;
    v = v_pred + gamma * dt * a1;
    a = a1;
    out.push_back(q);
  }
  return out;
}

Vector free_part(const FemModel& model, const Vector& full) {
  const std::vector<Index>& free = model.dofs().free_dofs();
  Vector out(static_cast<Eigen::Index>(free.size()));
  for (std::size_t i = 0; i < free.size(); ++i) out(static_cast<Eigen::Index>(i)) = full(free[i]);
  return out;
}

/// The global transient solution against the modal recursion: the largest
/// nodal difference over every step, over the largest displacement.
Scalar modal_error(const FemModel& model, const Assembler& assembler,
                   const TransientOptions& options) {
  const TransientResult r = solve_transient(model, assembler, 0, options);
  REQUIRE(static_cast<int>(r.snapshots.size()) == r.num_steps + 1);
  const ModalBasis basis = modal_basis(assembler, options.mass_type);
  const Vector f = free_part(model, model.load_vectors()[0]);
  const Vector load = basis.phi.transpose() * f;  // modal loads per unit amplitude
  const int steps = r.num_steps;
  const Scalar dt = options.time_step;
  Matrix q(basis.omega2.size(), steps + 1);
  for (Eigen::Index j = 0; j < basis.omega2.size(); ++j) {
    std::vector<Scalar> p;
    for (int n = 0; n <= steps; ++n) p.push_back(load(j) * options.amplitude.value(n * dt));
    const Scalar c = options.mass_damping + options.stiffness_damping * basis.omega2(j);
    // At rest q = 0; from the static state q = p(0) / w2 (the free DOFs'
    // static solution, since this model prescribes no displacement).
    const Scalar q0 = options.start == TransientOptions::Start::Static ? p[0] / basis.omega2(j)
                                                                       : 0.0;
    const std::vector<Scalar> history = integrate_mode(basis.omega2(j), c, p, dt,
                                                       options.alpha, q0, 0.0);
    for (int n = 0; n <= steps; ++n) q(j, n) = history[static_cast<std::size_t>(n)];
  }
  Scalar worst = 0.0;
  Scalar scale = 0.0;
  for (int n = 0; n <= steps; ++n) {
    const Vector reference = basis.phi * q.col(n);
    const Vector ours = free_part(model, r.snapshots[static_cast<std::size_t>(n)].displacement);
    worst = std::max(worst, (ours - reference).cwiseAbs().maxCoeff());
    scale = std::max(scale, reference.cwiseAbs().maxCoeff());
  }
  REQUIRE(scale > 0.0);
  return worst / scale;
}

}  // namespace

TEST_CASE("amplitudes evaluate, differentiate and validate", "[dynamics]") {
  Amplitude step;
  step.scale = 2.5;
  REQUIRE(step.value(0.0) == 2.5);
  REQUIRE(step.value(7.0) == 2.5);
  REQUIRE(step.rate(1.0) == 0.0);

  Amplitude table;
  table.kind = Amplitude::Kind::Table;
  table.times = {0.0, 1.0, 3.0};
  table.values = {0.0, 2.0, -2.0};
  REQUIRE(table.value(-1.0) == 0.0);
  REQUIRE(table.value(0.5) == Approx(1.0).epsilon(1e-15));
  REQUIRE(table.value(2.0) == Approx(0.0).margin(1e-15));
  REQUIRE(table.value(9.0) == -2.0);
  REQUIRE(table.rate(0.5) == Approx(2.0).epsilon(1e-15));
  REQUIRE(table.rate(1.0) == Approx(-2.0).epsilon(1e-15));  // the slope to the right
  REQUIRE(table.rate(4.0) == 0.0);

  Amplitude wave;
  wave.kind = Amplitude::Kind::Harmonic;
  wave.frequency = 5.0;
  wave.phase = 0.3;
  wave.scale = 2.0;
  const Scalar w = 2.0 * kPi * 5.0;
  REQUIRE(wave.value(0.01) == Approx(2.0 * std::sin(w * 0.01 + 0.3)).epsilon(1e-14));
  REQUIRE(wave.rate(0.01) == Approx(2.0 * w * std::cos(w * 0.01 + 0.3)).epsilon(1e-14));
  REQUIRE(wave.second_rate(0.01) ==
          Approx(-2.0 * w * w * std::sin(w * 0.01 + 0.3)).epsilon(1e-14));

  Amplitude bad = table;
  bad.times = {0.0, 1.0, 1.0};
  REQUIRE_THROWS_AS(bad.validate(), ConfigError);
  bad.times = {0.0, 1.0};
  REQUIRE_THROWS_AS(bad.validate(), ConfigError);
  Amplitude silent = wave;
  silent.frequency = 0.0;
  REQUIRE_THROWS_AS(silent.validate(), ConfigError);
  REQUIRE_NOTHROW(table.validate());
  REQUIRE_NOTHROW(wave.validate());

  const HhtParameters trapezoidal = HhtParameters::from_alpha(0.0);
  REQUIRE(trapezoidal.beta == 0.25);
  REQUIRE(trapezoidal.gamma == 0.5);
  const HhtParameters damped = HhtParameters::from_alpha(-0.1);
  REQUIRE(damped.beta == Approx(0.3025).epsilon(1e-15));
  REQUIRE(damped.gamma == Approx(0.6).epsilon(1e-15));
  REQUIRE_THROWS_AS(HhtParameters::from_alpha(0.05), ConfigError);
  REQUIRE_THROWS_AS(HhtParameters::from_alpha(-0.4), ConfigError);
}

TEST_CASE("the scalar trapezoidal recursion has the exact discrete period", "[dynamics]") {
  // q'' + w^2 q = 0 from (q0, 0): the trapezoidal rule advances the phase by
  // theta per step with tan(theta / 2) = w dt / 2 and keeps the amplitude -
  // the reference recursion of the tests below, checked against it.
  const Scalar w = 3.0;
  const Scalar dt = 0.2;
  const std::vector<Scalar> p(200, 0.0);
  const std::vector<Scalar> q = integrate_mode(w * w, 0.0, p, dt, 0.0, 1.0, 0.0);
  const Scalar theta = 2.0 * std::atan(w * dt / 2.0);
  Scalar worst = 0.0;
  for (std::size_t n = 0; n < q.size(); ++n) {
    worst = std::max(worst, std::abs(q[n] - std::cos(static_cast<Scalar>(n) * theta)));
  }
  REQUIRE(worst <= 1.0e-12);
}

TEST_CASE("a transient run is the HHT-alpha recursion of every mode", "[dynamics][transient]") {
  SolidCantileverCase solid;
  solid.length = 0.3;
  solid.height = 0.06;
  solid.width = 0.04;
  solid.poisson = 0.3;
  const FemModel plate = make_small_plate(4, 3);
  const FemModel block = make_cantilever_3d(solid, 3, 2, 1);
  for (const FemModel* model : {&plate, &block}) {
    const Assembler assembler(*model);
    const ModalBasis probe = modal_basis(assembler, MassType::Consistent);
    const Scalar f1 = std::sqrt(probe.omega2(0)) / (2.0 * kPi);
    for (const MassType mass : {MassType::Consistent, MassType::Lumped}) {
      for (const Scalar alpha : {0.0, -0.1}) {
        TransientOptions options;
        options.mass_type = mass;
        options.alpha = alpha;
        options.time_step = 0.05 / f1;
        options.end_time = 40 * options.time_step;
        options.mass_damping = 0.02 * 2.0 * kPi * f1;
        options.stiffness_damping = 1.0e-3 / (2.0 * kPi * f1);
        options.snapshot_every = 1;
        options.amplitude.kind = Amplitude::Kind::Harmonic;
        options.amplitude.frequency = 0.7 * f1;
        options.amplitude.phase = 0.4;
        INFO((model == &plate ? "Q4" : "Hex8") << " " << to_string(mass) << " alpha " << alpha);
        REQUIRE(modal_error(*model, assembler, options) <= 1.0e-10);

        // A preloaded structure released: the static state, then the load
        // drops to zero.
        TransientOptions release = options;
        release.start = TransientOptions::Start::Static;
        release.amplitude = Amplitude();
        release.amplitude.kind = Amplitude::Kind::Table;
        release.amplitude.times = {0.0, 0.5 * options.time_step};
        release.amplitude.values = {1.0, 0.0};
        REQUIRE(modal_error(*model, assembler, release) <= 1.0e-10);
      }
    }
  }
}

TEST_CASE("the trapezoidal rule balances energy exactly and HHT-alpha dissipates it",
          "[dynamics][transient]") {
  const FemModel model = make_small_plate(4, 3);
  const Assembler assembler(model);
  const ModalBasis basis = modal_basis(assembler, MassType::Consistent);
  const Scalar f1 = std::sqrt(basis.omega2(0)) / (2.0 * kPi);
  TransientOptions options;
  options.time_step = 0.03 / f1;
  options.end_time = 300 * options.time_step;
  options.amplitude.kind = Amplitude::Kind::Harmonic;
  options.amplitude.frequency = 1.3 * f1;

  // Undamped: T + U = E_0 + W to round-off at every step.
  const TransientResult free_run = solve_transient(model, assembler, 0, options);
  REQUIRE(free_run.energy_balance_error <= 1.0e-12);
  REQUIRE(std::abs(free_run.numerical_dissipation) <=
          1.0e-12 * free_run.steps.back().external_work + 1.0e-18);
  Scalar peak = 0.0;
  for (const TransientStep& s : free_run.steps) peak = std::max(peak, s.kinetic_energy);
  REQUIRE(peak > 0.0);

  // Damped: the dissipated energy closes the balance just as exactly.
  TransientOptions damped = options;
  damped.mass_damping = 0.05 * 2.0 * kPi * f1;
  damped.stiffness_damping = 2.0e-3 / (2.0 * kPi * f1);
  const TransientResult damped_run = solve_transient(model, assembler, 0, damped);
  REQUIRE(damped_run.energy_balance_error <= 1.0e-12);
  REQUIRE(damped_run.steps.back().damping_energy > 0.0);

  // HHT-alpha: the method itself dissipates energy, E_0 + W - T - U > 0.
  TransientOptions hht = options;
  hht.alpha = -0.1;
  const TransientResult hht_run = solve_transient(model, assembler, 0, hht);
  REQUIRE(hht_run.numerical_dissipation > 1.0e-6 * hht_run.steps.back().external_work);
}

TEST_CASE("a prescribed motion balances energy and momentum with its reactions",
          "[dynamics][transient]") {
  // The left edge of a free plate is shaken in y; the plate has no other
  // support, so its y-momentum changes only by the reactions:
  // sum_y r = sum_y (M a + C v) (K and the stiffness part of C take rigid
  // translations to zero).
  StructuredMeshSpec spec;
  spec.nx = 4;
  spec.ny = 3;
  spec.lx = 0.4;
  spec.ly = 0.3;
  FemModel model(make_structured_quad_mesh(spec), default_material(), 0.005,
                 StressState::PlaneStress, IntegrationOptions());
  DisplacementConstraint shaken;
  Selector left;
  left.kind = SelectorKind::Box;
  left.xmax = 0.0;
  shaken.region.members.push_back(left);
  shaken.set(0, true, 0.0);
  shaken.set(1, true, 1.0e-4);
  model.constraints().push_back(shaken);
  LoadCaseSpec lc;
  lc.name = "shake";
  lc.prescribed_displacement_only = true;
  model.load_case_specs().push_back(lc);
  model.finalize();
  const Assembler assembler(model);
  const ModalBasis basis = modal_basis(assembler, MassType::Consistent);
  const Scalar f1 = std::sqrt(basis.omega2(0)) / (2.0 * kPi);

  TransientOptions options;
  options.time_step = 0.04 / f1;
  options.end_time = 150 * options.time_step;
  options.amplitude.kind = Amplitude::Kind::Harmonic;
  options.amplitude.frequency = 0.8 * f1;
  options.mass_damping = 0.03 * 2.0 * kPi * f1;
  options.stiffness_damping = 1.0e-3 / (2.0 * kPi * f1);
  DynamicMonitor reaction;
  reaction.name = "left_force_y";
  reaction.region.members.push_back(left);
  reaction.component = 1;
  reaction.quantity = DynamicMonitor::Quantity::Reaction;
  options.monitors.push_back(reaction);
  const TransientResult r = solve_transient(model, assembler, 0, options);
  REQUIRE(r.energy_balance_error <= 1.0e-12);

  const SparseMatrix m = assembler.assemble_mass(MassType::Consistent);
  const SparseMatrix k = assembler.assemble_stiffness();
  const Vector inertia = m * r.acceleration +
                         options.mass_damping * (m * r.velocity) +
                         options.stiffness_damping * (k * r.velocity);
  Scalar momentum_rate = 0.0;
  Scalar gross = 0.0;
  for (Index node = 0; node < model.mesh().num_nodes(); ++node) {
    momentum_rate += inertia(node * 2 + 1);
    gross += std::abs(inertia(node * 2 + 1));
  }
  const Scalar reaction_sum = r.steps.back().monitors[0];
  REQUIRE(std::abs(reaction_sum - momentum_rate) <= 1.0e-10 * gross);
  REQUIRE(std::abs(reaction_sum) > 1.0e-6 * gross);
}

TEST_CASE("the harmonic response is the sum over all modes", "[dynamics][frequency]") {
  SolidCantileverCase solid;
  solid.length = 0.3;
  solid.height = 0.06;
  solid.width = 0.04;
  solid.poisson = 0.3;
  const FemModel plate = make_small_plate(4, 3);
  const FemModel block = make_cantilever_3d(solid, 3, 2, 1);
  for (const FemModel* model : {&plate, &block}) {
    const Assembler assembler(*model);
    for (const MassType mass : {MassType::Consistent, MassType::Lumped}) {
      const ModalBasis basis = modal_basis(assembler, mass);
      const Scalar f1 = std::sqrt(basis.omega2(0)) / (2.0 * kPi);
      FrequencyResponseOptions options;
      options.mass_type = mass;
      options.structural_damping = 0.02;
      options.mass_damping = 0.01 * 2.0 * kPi * f1;
      options.stiffness_damping = 1.0e-4 / (2.0 * kPi * f1);
      options.frequencies = {0.0, 0.5 * f1, f1, 1.7 * f1, 3.0 * f1};
      options.snapshot_frequencies = options.frequencies;
      DynamicMonitor tip;
      tip.name = "tip";
      Selector all;
      tip.region.members.push_back(all);
      tip.component = 1;
      for (const auto q : {DynamicMonitor::Quantity::Displacement,
                           DynamicMonitor::Quantity::Velocity,
                           DynamicMonitor::Quantity::Acceleration}) {
        tip.quantity = q;
        options.monitors.push_back(tip);
      }
      const FrequencyResponseResult r = solve_frequency_response(*model, assembler, 0, options);
      REQUIRE(r.snapshots.size() == options.frequencies.size());
      const Vector f = free_part(*model, model->load_vectors()[0]);
      const Vector load = basis.phi.transpose() * f;
      const std::complex<Scalar> i_unit(0.0, 1.0);
      INFO((model == &plate ? "Q4" : "Hex8") << " " << to_string(mass));
      for (std::size_t j = 0; j < options.frequencies.size(); ++j) {
        const Scalar omega = 2.0 * kPi * options.frequencies[j];
        ComplexVector reference = ComplexVector::Zero(basis.phi.rows());
        for (Eigen::Index m = 0; m < basis.omega2.size(); ++m) {
          const Scalar w2 = basis.omega2(m);
          const std::complex<Scalar> d =
              w2 * (1.0 + i_unit * options.structural_damping) - omega * omega +
              i_unit * omega * (options.mass_damping + options.stiffness_damping * w2);
          reference += basis.phi.col(m).cast<std::complex<Scalar>>() * (load(m) / d);
        }
        const std::vector<Index>& free = model->dofs().free_dofs();
        Scalar worst = 0.0;
        for (std::size_t i = 0; i < free.size(); ++i) {
          worst = std::max(worst, std::abs(r.snapshots[j].displacement(free[i]) -
                                           reference(static_cast<Eigen::Index>(i))));
        }
        REQUIRE(worst <= 1.0e-10 * reference.cwiseAbs().maxCoeff());
        // Velocity and acceleration monitors are i omega U and -omega^2 U.
        const std::complex<Scalar> u = r.points[j].monitors[0];
        REQUIRE(std::abs(r.points[j].monitors[1] - i_unit * omega * u) <=
                1.0e-14 * (std::abs(u) * omega + 1e-300));
        REQUIRE(std::abs(r.points[j].monitors[2] + omega * omega * u) <=
                1.0e-14 * (std::abs(u) * omega * omega + 1e-300));
      }
    }
  }
}

TEST_CASE("the harmonic peak of a node is the semi-major axis of its orbit",
          "[dynamics][frequency]") {
  const FemModel model = make_small_plate(1, 1);  // four nodes
  const std::complex<Scalar> i_unit(0.0, 1.0);
  ComplexVector u = ComplexVector::Zero(model.dofs().num_dofs());
  const auto set = [&](Index node, std::complex<Scalar> ux, std::complex<Scalar> uy) {
    u(model.dofs().dof(node, 0)) = ux;
    u(model.dofs().dof(node, 1)) = uy;
  };
  set(0, 3.0, 4.0);                        // in phase: a straight line of amplitude 5
  set(1, 2.0, 2.0 * i_unit);               // a circle of radius 2 (modulus 2 sqrt 2)
  set(2, 3.0, 1.0 * i_unit);               // an ellipse with semi-axes 3 and 1
  set(3, (1.0 + i_unit) * 0.6, (1.0 + i_unit) * 0.8);  // in phase at 45 degrees: 1 sqrt 2
  const Vector peak = harmonic_peak_displacements(model, u);
  REQUIRE(peak(0) == Approx(5.0).epsilon(1e-15));
  REQUIRE(peak(1) == Approx(2.0).epsilon(1e-15));
  REQUIRE(peak(2) == Approx(3.0).epsilon(1e-15));
  REQUIRE(peak(3) == Approx(std::sqrt(2.0)).epsilon(1e-15));
  // Sampled over a cycle, the orbit's largest radius is the peak.
  for (Index node = 0; node < 4; ++node) {
    Scalar sampled = 0.0;
    for (int k = 0; k < 3600; ++k) {
      const std::complex<Scalar> phase = std::exp(i_unit * (2.0 * kPi * k / 3600.0));
      const Scalar x = (u(model.dofs().dof(node, 0)) * phase).real();
      const Scalar y = (u(model.dofs().dof(node, 1)) * phase).real();
      sampled = std::max(sampled, std::hypot(x, y));
    }
    REQUIRE(sampled <= peak(node) * (1.0 + 1e-14));
    REQUIRE(sampled >= peak(node) * (1.0 - 1e-6));
  }
}

TEST_CASE("the dynamic analyses refuse inconsistent input", "[dynamics]") {
  const FemModel model = make_small_plate(4, 3);
  const Assembler assembler(model);
  TransientOptions options;
  options.time_step = 1.0e-4;
  options.end_time = 1.05e-3;  // not a whole number of steps
  REQUIRE_THROWS_AS(solve_transient(model, assembler, 0, options), ConfigError);
  options.end_time = 1.0e-3;
  options.alpha = -0.5;
  REQUIRE_THROWS_AS(solve_transient(model, assembler, 0, options), ConfigError);
  options.alpha = 0.0;
  options.mass_damping = -1.0;
  REQUIRE_THROWS_AS(solve_transient(model, assembler, 0, options), ConfigError);
  options.mass_damping = 0.0;
  REQUIRE_NOTHROW(solve_transient(model, assembler, 0, options));

  FrequencyResponseOptions harmonic;
  REQUIRE_THROWS_AS(solve_frequency_response(model, assembler, 0, harmonic), ConfigError);
  harmonic.frequencies = {-1.0};
  REQUIRE_THROWS_AS(solve_frequency_response(model, assembler, 0, harmonic), ConfigError);

  // A massless material has no dynamics.
  FemModel massless(make_structured_quad_mesh([] {
                      StructuredMeshSpec s;
                      s.nx = 2;
                      s.ny = 2;
                      return s;
                    }()),
                    IsotropicMaterial(70.0e9, 0.3, 0.0, "massless"), 0.01,
                    StressState::PlaneStress, IntegrationOptions());
  DisplacementConstraint left;
  Selector box;
  box.kind = SelectorKind::Box;
  box.xmax = 0.0;
  left.region.members.push_back(box);
  left.fix_x = left.fix_y = true;
  massless.constraints().push_back(left);
  LoadCaseSpec lc;
  lc.name = "pull";
  PointLoadSpec pull;
  Selector right = box;
  right.xmax = std::numeric_limits<Scalar>::infinity();
  right.xmin = 1.0;
  pull.region.members.push_back(right);
  pull.force = Vector3(1.0, 0.0, 0.0);
  lc.point_loads.push_back(pull);
  massless.load_case_specs().push_back(lc);
  massless.finalize();
  const Assembler massless_assembler(massless);
  options.end_time = 1.0e-3;
  REQUIRE_THROWS_AS(solve_transient(massless, massless_assembler, 0, options), ModelError);
}

namespace {

/// The plate of make_small_plate with a plastic steel-like material.
FemModel plastic_plate() {
  StructuredMeshSpec spec;
  spec.nx = 18;  // square cells: an oblong fully integrated Q4 locks in bending
  spec.ny = 3;
  spec.lx = 0.6;
  spec.ly = 0.1;
  IsotropicMaterial steel(200.0e9, 0.3, 7800.0, "steel");
  PlasticityParameters plasticity;
  plasticity.yield_stress = 250.0e6;
  plasticity.hardening_modulus = 2.0e9;
  steel.set_plasticity(plasticity);
  FemModel model(make_structured_quad_mesh(spec), steel, 0.01, StressState::PlaneStrain,
                 IntegrationOptions());
  DisplacementConstraint root;
  Selector left;
  left.kind = SelectorKind::Box;
  left.xmax = 0.0;
  root.region.members.push_back(left);
  root.fix_x = root.fix_y = true;
  model.constraints().push_back(root);
  LoadCaseSpec lc;
  lc.name = "tip";
  PointLoadSpec tip;
  Selector right;
  right.kind = SelectorKind::Box;
  right.xmin = 0.6;
  tip.region.members.push_back(right);
  tip.force = Vector3(0.0, -1.3e4, 0.0);  // plastic strain 2e-3, rotation 0.016 rad
  tip.distribute_total = true;
  lc.point_loads.push_back(tip);
  model.load_case_specs().push_back(lc);
  model.finalize();
  return model;
}

}  // namespace

TEST_CASE("the non-linear transient of a linear model is the linear transient",
          "[dynamics][transient][nonlinear]") {
  // Small-strain kinematics with an elastic material is a linear system:
  // Newton converges in one iteration to the linear step.
  const FemModel model = make_small_plate(4, 3);
  const Assembler assembler(model);
  const ModalBasis basis = modal_basis(assembler, MassType::Consistent);
  const Scalar f1 = std::sqrt(basis.omega2(0)) / (2.0 * kPi);
  for (const Scalar alpha : {0.0, -0.1}) {
    TransientOptions options;
    options.alpha = alpha;
    options.time_step = 0.05 / f1;
    options.end_time = 60 * options.time_step;
    options.mass_damping = 0.02 * 2.0 * kPi * f1;
    options.stiffness_damping = 1.0e-3 / (2.0 * kPi * f1);
    options.amplitude.kind = Amplitude::Kind::Harmonic;
    options.amplitude.frequency = 0.9 * f1;
    options.snapshot_every = 1;
    const TransientResult linear = solve_transient(model, assembler, 0, options);
    options.nonlinear = true;
    options.nonlinear_options.kinematics = Kinematics::SmallStrain;
    options.nonlinear_options.residual_tolerance = 1.0e-12;
    options.nonlinear_options.displacement_tolerance = 1.0e-12;
    const TransientResult nonlinear = solve_transient(model, assembler, 0, options);
    REQUIRE(nonlinear.completed);
    REQUIRE(nonlinear.snapshots.size() == linear.snapshots.size());
    Scalar worst = 0.0;
    Scalar scale = 0.0;
    for (std::size_t i = 0; i < linear.snapshots.size(); ++i) {
      worst = std::max(worst, (nonlinear.snapshots[i].displacement -
                               linear.snapshots[i].displacement).cwiseAbs().maxCoeff());
      scale = std::max(scale, linear.snapshots[i].displacement.cwiseAbs().maxCoeff());
    }
    INFO("alpha " << alpha);
    REQUIRE(worst <= 1.0e-10 * scale);
    for (const TransientStep& s : nonlinear.steps) REQUIRE(s.iterations <= 2);
  }
}

TEST_CASE("the non-linear transient converges at second order in the step",
          "[dynamics][transient][nonlinear]") {
  // A slender elastic strip with finite kinematics, its tip load ramped up
  // smoothly over a quarter of the first period (applied suddenly at a node,
  // it would crush the loaded element - a spurious state the run warns
  // about): the displacement after half a period converges at the
  // trapezoidal rule's second order, and the energy balance closes to its
  // O(dt^2) error.
  CantileverCase c;
  c.length = 1.0;
  c.height = 0.05;
  c.thickness = 0.05;
  c.youngs = 1.0e9;
  c.density = 1000.0;
  c.poisson = 0.3;
  c.tip_load = -300.0;  // a tip deflection of a sixth of the span
  const FemModel model = make_cantilever(c, 20, 2);
  const Assembler assembler(model);
  const ModalBasis basis = modal_basis(assembler, MassType::Consistent);
  const Scalar period = 2.0 * kPi / std::sqrt(basis.omega2(0));
  std::vector<Scalar> balance;
  const auto tip_history_end = [&](int per_period) {
    TransientOptions options;
    options.time_step = period / per_period;
    options.end_time = 0.5 * period;
    options.nonlinear = true;
    options.nonlinear_options.residual_tolerance = 1.0e-11;
    options.nonlinear_options.displacement_tolerance = 1.0e-11;
    options.amplitude.kind = Amplitude::Kind::Table;
    for (int i = 0; i <= 40; ++i) {
      options.amplitude.times.push_back(0.25 * period * i / 40.0);
      const Scalar s = std::sin(0.5 * kPi * i / 40.0);
      options.amplitude.values.push_back(s * s);
    }
    const TransientResult r = solve_transient(model, assembler, 0, options);
    REQUIRE(r.completed);
    REQUIRE(r.warnings.empty());
    balance.push_back(r.energy_balance_error);
    return r.displacement;
  };
  const Vector coarse = tip_history_end(40);
  const Vector medium = tip_history_end(80);
  const Vector fine = tip_history_end(160);
  const Vector finest = tip_history_end(640);
  const Scalar e1 = (coarse - finest).cwiseAbs().maxCoeff();
  const Scalar e2 = (medium - finest).cwiseAbs().maxCoeff();
  const Scalar e3 = (fine - finest).cwiseAbs().maxCoeff();
  const Scalar order_a = std::log2(e1 / e2);
  const Scalar order_b = std::log2(e2 / e3);
  INFO("errors " << e1 << " " << e2 << " " << e3);
  REQUIRE(order_a > 1.8);
  REQUIRE(order_b > 1.8);
  INFO("energy balance " << balance[0] << " " << balance[1] << " " << balance[2]);
  REQUIRE(std::log2(balance[0] / balance[1]) > 1.8);
  REQUIRE(balance[0] < 1.0e-3);
  // The deflection is large: finite and small-strain kinematics differ.
  REQUIRE(finest.cwiseAbs().maxCoeff() > 0.1 * c.length);
}

TEST_CASE("a slow elastoplastic transient approaches the static solution",
          "[dynamics][transient][nonlinear][plasticity]") {
  // Ramped smoothly and damped, the load leaves the plate at rest in the
  // state the non-linear static analysis finds, with its plastic strain.
  // While the ramp runs, the damping force - about 2 zeta / (omega T_ramp)
  // of the load, spread like the mass rather than at the tip - shifts the
  // plastic loading path a little, so the gap to the static state closes as
  // the ramp lengthens: doubling it must shrink the gap. The energy the
  // plastic flow dissipated is positive.
  const FemModel model = plastic_plate();
  const Assembler assembler(model);
  const ModalBasis basis = modal_basis(assembler, MassType::Consistent);
  const Scalar f1 = std::sqrt(basis.omega2(0)) / (2.0 * kPi);
  const Scalar period = 1.0 / f1;

  // The static reference in 320 load steps: the backward-Euler return
  // integrates the non-proportional stress path of bending at first order in
  // the step, and 20 steps would leave it 0.7 % from its limit (measured).
  NonlinearOptions statics;
  statics.kinematics = Kinematics::SmallStrain;
  statics.steps = 320;
  statics.residual_tolerance = 1.0e-11;
  statics.displacement_tolerance = 1.0e-11;
  const NonlinearResult reference = NonlinearStaticAnalysis(model, assembler, statics).solve(0);
  REQUIRE(reference.completed);
  REQUIRE(reference.max_plastic_strain > 1.0e-3);

  const auto gap = [&](Scalar ramp_periods, Scalar& dissipation, Scalar& plastic) {
    TransientOptions options;
    options.time_step = period / 20.0;
    options.end_time = (ramp_periods + 30.0) * period;
    options.mass_damping = 0.3 * 2.0 * kPi * f1;  // 30 % of critical in the first mode
    // A smooth ramp, sin^2: a linear ramp's kinks would each start an
    // oscillation of about 1 / (omega T_ramp) of the deflection, and near
    // collapse the overshoot at its end becomes permanent plastic flow.
    options.amplitude.kind = Amplitude::Kind::Table;
    for (int i = 0; i <= 60; ++i) {
      options.amplitude.times.push_back(ramp_periods * period * i / 60.0);
      const Scalar s = std::sin(0.5 * kPi * i / 60.0);
      options.amplitude.values.push_back(s * s);
    }
    options.nonlinear = true;
    options.nonlinear_options = statics;
    const TransientResult r = solve_transient(model, assembler, 0, options);
    REQUIRE(r.completed);
    REQUIRE(r.plastic);
    dissipation = r.numerical_dissipation;
    plastic = r.max_plastic_strain;
    return (r.displacement - reference.displacement).cwiseAbs().maxCoeff() /
           reference.displacement.cwiseAbs().maxCoeff();
  };
  Scalar dissipation_a = 0.0, plastic_a = 0.0, dissipation_b = 0.0, plastic_b = 0.0;
  const Scalar gap_a = gap(30.0, dissipation_a, plastic_a);
  const Scalar gap_b = gap(60.0, dissipation_b, plastic_b);
  INFO("gap to the static state: ramp 30 T " << gap_a << ", 60 T " << gap_b);
  REQUIRE(gap_b < 0.7 * gap_a);
  REQUIRE(gap_b < 2.0e-3);
  REQUIRE(plastic_b == Approx(reference.max_plastic_strain).epsilon(1.0e-2));
  // E_0 + W - T - U - D: the plastic dissipation (with the method's own
  // small error), positive.
  REQUIRE(dissipation_a > 0.0);
  REQUIRE(dissipation_b > 0.0);
}

TEST_CASE("an undamped frequency response at a natural frequency is flagged",
          "[dynamics][frequency]") {
  const FemModel model = make_small_plate(4, 3);
  const Assembler assembler(model);
  const ModalBasis basis = modal_basis(assembler, MassType::Consistent);
  const Scalar f1 = std::sqrt(basis.omega2(0)) / (2.0 * kPi);
  FrequencyResponseOptions options;
  options.frequencies = {0.5 * f1, f1};
  FrequencyResponseResult r;
  // Exactly singular to working precision, the LU fails outright; otherwise
  // the answer is flagged as round-off.
  try {
    r = solve_frequency_response(model, assembler, 0, options);
    REQUIRE(r.warnings.size() == 1);
    REQUIRE(r.warnings[0].find("undamped natural frequency") != std::string::npos);
  } catch (const SolverError& error) {
    REQUIRE(std::string(error.what()).find("natural frequency") != std::string::npos);
  }
  // With damping the resonant response is finite and unflagged.
  options.structural_damping = 0.02;
  r = solve_frequency_response(model, assembler, 0, options);
  REQUIRE(r.warnings.empty());
  REQUIRE(r.points[1].max_displacement > 10.0 * r.points[0].max_displacement);
}
