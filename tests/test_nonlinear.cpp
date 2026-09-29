/// \file test_nonlinear.cpp
/// \brief The large-deformation analysis: material laws, the total Lagrangian
///        element, the complete tangent of the solver's residual, and exact
///        solutions of homogeneous deformation.
///
/// Every check is exact or a derivative check: the laws' stresses against
/// their energies and their tangents against their stresses, the element
/// tangent against central differences of the internal force, objectivity
/// under a rigid rotation, the small-strain limit, homogeneous deformations
/// that every element reproduces exactly, and the closed-form stretch of a
/// cube under dead and follower pressure.
#include "TestSupport.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/fem/NonlinearStatic.hpp"
#include "sparlab/fem/TotalLagrangian.hpp"
#include "sparlab/io/Config.hpp"
#include "sparlab/material/Hyperelastic.hpp"

#include <Eigen/Eigenvalues>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>

using namespace sparlab;
using namespace sparlab::testing;
using Catch::Approx;

namespace {

constexpr Scalar kInf = std::numeric_limits<Scalar>::infinity();

SelectorGroup box_region(Scalar xmin, Scalar xmax, Scalar ymin = -kInf, Scalar ymax = kInf,
                         Scalar zmin = -kInf, Scalar zmax = kInf) {
  SelectorGroup g;
  Selector s;
  s.kind = SelectorKind::Box;
  s.xmin = xmin;
  s.xmax = xmax;
  s.ymin = ymin;
  s.ymax = ymax;
  s.zmin = zmin;
  s.zmax = zmax;
  g.members.push_back(s);
  g.name = "box";
  return g;
}

SelectorGroup point_region(const Vector3& x) {
  SelectorGroup g;
  Selector s;
  s.kind = SelectorKind::NearestNode;
  s.point = x;
  g.members.push_back(s);
  g.name = "point";
  return g;
}

DisplacementConstraint fix(const SelectorGroup& region, bool x, bool y, bool z) {
  DisplacementConstraint bc;
  bc.region = region;
  bc.fix_x = x;
  bc.fix_y = y;
  bc.fix_z = z;
  return bc;
}

Matrix random_matrix(int rows, int cols, Scalar amplitude, unsigned seed) {
  std::mt19937 gen(seed);
  std::uniform_real_distribution<Scalar> dist(-1.0, 1.0);
  Matrix m(rows, cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) m(i, j) = amplitude * dist(gen);
  }
  return m;
}

Mesh single_cell(ElementType type) {
  StructuredMeshSpec spec;
  spec.nx = spec.ny = spec.nz = 1;
  spec.lx = 1.0;
  spec.ly = 0.8;
  spec.lz = 0.6;
  switch (type) {
    case ElementType::Quad4: return make_perturbed_quad_mesh([&] { auto s = spec; s.nx = 2; s.ny = 2; return s; }(), 0.2);
    case ElementType::Tri3: return make_perturbed_tri_mesh([&] { auto s = spec; s.nx = 2; s.ny = 2; return s; }(), 0.2);
    case ElementType::Hex8: return make_perturbed_hex_mesh([&] { auto s = spec; s.nx = s.ny = s.nz = 2; return s; }(), 0.2);
    case ElementType::Tet4: return make_perturbed_tet_mesh([&] { auto s = spec; s.nx = s.ny = s.nz = 2; return s; }(), 0.2);
    case ElementType::Tet10: return make_perturbed_tet10_mesh([&] { auto s = spec; s.nx = s.ny = s.nz = 2; return s; }(), 0.2);
  }
  throw ConfigError("unhandled element type");
}

}  // namespace

TEST_CASE("the hyperelastic laws are consistent: energy, stress and tangent",
          "[nonlinear][material]") {
  IsotropicMaterial m(200.0e9, 0.3, 7800.0, "steel");
  struct Case {
    HyperelasticModel law;
    StressState state;
    int dim;
  };
  const std::vector<Case> cases = {
      {HyperelasticModel::SaintVenantKirchhoff, StressState::PlaneStress, 2},
      {HyperelasticModel::SaintVenantKirchhoff, StressState::PlaneStrain, 2},
      {HyperelasticModel::SaintVenantKirchhoff, StressState::ThreeDimensional, 3},
      {HyperelasticModel::NeoHookean, StressState::PlaneStrain, 2},
      {HyperelasticModel::NeoHookean, StressState::ThreeDimensional, 3}};
  for (const Case& c : cases) {
    // A large deformation: stretches, shears and a rotation.
    const Matrix grad = random_matrix(c.dim, c.dim, 0.3, 7u);
    const HyperelasticResponse r = evaluate_hyperelastic(c.law, m, c.state, grad);
    const Matrix dg = random_matrix(c.dim, c.dim, 1.0, 11u);
    const Scalar h = 1.0e-6;
    const HyperelasticResponse rp = evaluate_hyperelastic(c.law, m, c.state, grad + h * dg);
    const HyperelasticResponse rm = evaluate_hyperelastic(c.law, m, c.state, grad - h * dg);
    const Vector de = (green_lagrange_voigt(grad + h * dg) - green_lagrange_voigt(grad - h * dg)) /
                      (2.0 * h);
    // dW = S : dE and dS = C dE.
    const Scalar dw = (rp.energy - rm.energy) / (2.0 * h);
    REQUIRE(dw == Approx(r.stress.dot(de)).epsilon(1e-7));
    const Vector ds = (rp.stress - rm.stress) / (2.0 * h);
    const Vector predicted = r.tangent * de;
    REQUIRE((ds - predicted).norm() <= 1e-6 * predicted.norm());
    // The tangent of a hyperelastic law is symmetric.
    REQUIRE((r.tangent - r.tangent.transpose()).norm() <= 1e-12 * r.tangent.norm());
  }
}

TEST_CASE("both laws reduce to linear elasticity for small strain", "[nonlinear][material]") {
  IsotropicMaterial m(70.0e9, 0.33, 2700.0);
  for (const HyperelasticModel law :
       {HyperelasticModel::SaintVenantKirchhoff, HyperelasticModel::NeoHookean}) {
    const Matrix grad = random_matrix(3, 3, 1.0e-7, 3u);
    const HyperelasticResponse r =
        evaluate_hyperelastic(law, m, StressState::ThreeDimensional, grad);
    const Matrix eps = 0.5 * (grad + grad.transpose());
    Vector ev(6);
    ev << eps(0, 0), eps(1, 1), eps(2, 2), 2.0 * eps(0, 1), 2.0 * eps(1, 2), 2.0 * eps(2, 0);
    const Vector linear = m.three_dimensional_matrix() * ev;
    // The difference is second order in the strain, 1e-7 here: no round-off
    // floor above it now that the laws work from the displacement gradient.
    REQUIRE((r.stress - linear).norm() <= 1e-6 * linear.norm());
    REQUIRE((r.tangent - m.three_dimensional_matrix()).norm() <=
            1e-5 * m.three_dimensional_matrix().norm());
  }
}

TEST_CASE("the neo-Hookean law refuses plane stress and a temperature", "[nonlinear][material]") {
  IsotropicMaterial m(1.0e6, 0.45, 1000.0);
  REQUIRE_THROWS_AS(evaluate_hyperelastic(HyperelasticModel::NeoHookean, m,
                                          StressState::PlaneStress, Matrix::Zero(2, 2)),
                    ConfigError);
  m.set_thermal(1.0e-4, 20.0, 0.2);
  REQUIRE_THROWS_AS(evaluate_hyperelastic(HyperelasticModel::NeoHookean, m,
                                          StressState::ThreeDimensional, Matrix::Zero(3, 3),
                                          10.0),
                    ConfigError);
  // An inverted material point is reported, not evaluated.
  m.set_thermal(0.0, 0.0, 0.0);
  REQUIRE_THROWS_AS(evaluate_hyperelastic(HyperelasticModel::NeoHookean, m,
                                          StressState::ThreeDimensional,
                                          -2.0 * Matrix::Identity(3, 3)),
                    SolverError);
  REQUIRE(parse_hyperelastic_model("neo_hookean") == HyperelasticModel::NeoHookean);
  REQUIRE_THROWS_AS(parse_hyperelastic_model("mooney"), ConfigError);
}

TEST_CASE("the total Lagrangian tangent is the derivative of the internal force",
          "[nonlinear][element]") {
  for (const ElementType type : {ElementType::Quad4, ElementType::Tri3, ElementType::Hex8,
                                 ElementType::Tet4, ElementType::Tet10}) {
    for (const HyperelasticModel law :
         {HyperelasticModel::SaintVenantKirchhoff, HyperelasticModel::NeoHookean}) {
      Mesh mesh = single_cell(type);
      const int dim = mesh.dim();
      const StressState state =
          dim == 2 ? StressState::PlaneStrain : StressState::ThreeDimensional;
      IsotropicMaterial mat(200.0e9, 0.3, 7800.0);
      FemModel model(std::move(mesh), mat, 1.0, state, IntegrationOptions());
      model.constraints().push_back(fix(box_region(-kInf, 0.0), true, true, dim == 3));
      LoadCaseSpec lc;
      lc.name = "load";
      PointLoadSpec p;
      p.region = point_region(Vector3(1.0, 0.8, 0.6));
      p.force = Vector3(0.0, -1.0, 0.0);
      lc.point_loads.push_back(p);
      model.load_case_specs().push_back(lc);
      model.finalize();
      const Index e = 0;
      const int nd = dim * model.mesh().nodes_per_elem();
      // Large displacements: ~15 % of the cell plus a rigid turn of 0.4 rad.
      const Matrix x0 = model.mesh().element_coordinates(e);
      Matrix rot = Matrix::Identity(dim, dim);
      rot(0, 0) = std::cos(0.4);
      rot(0, 1) = -std::sin(0.4);
      rot(1, 0) = std::sin(0.4);
      rot(1, 1) = std::cos(0.4);
      // (Strains of up to ~10 %: enough to exercise every non-linear term,
      // short of inverting the cell.)
      const Matrix disp = (rot - Matrix::Identity(dim, dim)) * x0 +
                          random_matrix(dim, static_cast<int>(x0.cols()), 0.02, 5u);
      Vector ue(nd);
      for (int a = 0; a < model.mesh().nodes_per_elem(); ++a) {
        for (int k = 0; k < dim; ++k) ue(dim * a + k) = disp(k, a);
      }
      const TotalLagrangianElement base =
          total_lagrangian_element(model, e, ue, law, nullptr, 1.0, true);
      const Scalar h = 1.0e-7;
      Scalar worst = 0.0;
      for (int j = 0; j < nd; ++j) {
        Vector up = ue;
        Vector um = ue;
        up(j) += h;
        um(j) -= h;
        const Vector fd = (total_lagrangian_element(model, e, up, law, nullptr, 1.0, false)
                               .internal_force -
                           total_lagrangian_element(model, e, um, law, nullptr, 1.0, false)
                               .internal_force) /
                          (2.0 * h);
        worst = std::max(worst, (fd - base.tangent.col(j)).cwiseAbs().maxCoeff());
      }
      INFO(to_string(type) << " " << to_string(law));
      REQUIRE(worst <= 1.0e-6 * base.tangent.cwiseAbs().maxCoeff());
      // And the internal force is the derivative of the energy.
      const Vector dir = random_matrix(nd, 1, 1.0, 13u);
      const Scalar de =
          (total_lagrangian_element(model, e, ue + h * dir, law, nullptr, 1.0, false).energy -
           total_lagrangian_element(model, e, ue - h * dir, law, nullptr, 1.0, false).energy) /
          (2.0 * h);
      REQUIRE(de == Approx(base.internal_force.dot(dir)).epsilon(1e-6));
    }
  }
}

TEST_CASE("a rigid rotation leaves no internal force and the reference tangent is linear",
          "[nonlinear][element]") {
  for (const ElementType type : {ElementType::Quad4, ElementType::Hex8, ElementType::Tet10}) {
    Mesh mesh = single_cell(type);
    const int dim = mesh.dim();
    const StressState state = dim == 2 ? StressState::PlaneStress : StressState::ThreeDimensional;
    IsotropicMaterial mat(200.0e9, 0.3, 7800.0);
    FemModel model(std::move(mesh), mat, dim == 2 ? 0.05 : 1.0, state, IntegrationOptions());
    const Matrix x0 = model.mesh().element_coordinates(0);
    const int nd = dim * model.mesh().nodes_per_elem();
    // A 90-degree turn about z (and a translation).
    Matrix rot = Matrix::Zero(dim, dim);
    rot(0, 1) = -1.0;
    rot(1, 0) = 1.0;
    if (dim == 3) rot(2, 2) = 1.0;
    const Matrix disp = (rot - Matrix::Identity(dim, dim)) * x0 +
                        Vector::Constant(dim, 0.3).replicate(1, x0.cols());
    Vector ue(nd);
    for (int a = 0; a < model.mesh().nodes_per_elem(); ++a) {
      for (int k = 0; k < dim; ++k) ue(dim * a + k) = disp(k, a);
    }
    for (const HyperelasticModel law :
         {HyperelasticModel::SaintVenantKirchhoff, HyperelasticModel::NeoHookean}) {
      if (law == HyperelasticModel::NeoHookean && state == StressState::PlaneStress) continue;
      const TotalLagrangianElement rigid =
          total_lagrangian_element(model, 0, ue, law, nullptr, 1.0, false);
      INFO(to_string(type) << " " << to_string(law));
      // Scale: E times the cell's measure (below 1 here), times round-off -
      // a rotation by 90 degrees leaves E = 0 only to round-off in H^T H.
      REQUIRE(rigid.internal_force.cwiseAbs().maxCoeff() <= 1.0e-12 * 200.0e9);
      REQUIRE(std::abs(rigid.energy) <= 1.0e-13 * 200.0e9);
      // At zero displacement the tangent is the linear stiffness.
      const TotalLagrangianElement zero =
          total_lagrangian_element(model, 0, Vector::Zero(nd), law, nullptr, 1.0, true);
      const Matrix k = model.element().stiffness(x0, model.constitutive(),
                                                 dim == 2 ? model.thickness() : 1.0,
                                                 model.integration());
      REQUIRE((zero.tangent - k).cwiseAbs().maxCoeff() <= 1.0e-10 * k.cwiseAbs().maxCoeff());
    }
  }
}

TEST_CASE("the solver's residual has the tangent and load rate it iterates with",
          "[nonlinear][solver]") {
  // A Hex8 block with every non-linear load: a follower pressure, a rotation
  // at the deformed position (spin softening), self-weight and a
  // temperature, at a large displacement state.
  StructuredMeshSpec spec;
  spec.nx = 3;
  spec.ny = 2;
  spec.nz = 2;
  spec.lx = 0.3;
  spec.ly = 0.1;
  spec.lz = 0.1;
  IsotropicMaterial steel(200.0e9, 0.3, 7800.0, "steel");
  // alpha dT = 0.05 at lambda = 1, so that the thermal stretch's own
  // non-linearity (the 1/theta of the multiplicative split) is exercised.
  steel.set_thermal(5.0e-4, 20.0, 45.0);
  FemModel model(make_perturbed_hex_mesh(spec, 0.2), steel, 1.0,
                 StressState::ThreeDimensional, IntegrationOptions());
  model.constraints().push_back(fix(box_region(-kInf, 0.0), true, true, true));
  LoadCaseSpec lc;
  lc.name = "all";
  PressureLoadSpec pressure;
  pressure.region = box_region(-kInf, kInf, -kInf, kInf, 0.1, kInf);
  pressure.pressure = 5.0e7;
  lc.pressures.push_back(pressure);
  lc.gravity = Vector3(0.0, 0.0, -9.81);
  lc.centrifugal.enabled = true;
  lc.centrifugal.angular_velocity = 300.0;
  lc.centrifugal.axis = Vector3(0.0, 0.0, 1.0);
  lc.temperature.source = TemperatureSpec::Source::Uniform;
  lc.temperature.uniform = 120.0;
  model.load_case_specs().push_back(lc);
  model.finalize();
  Assembler assembler(model);
  NonlinearOptions options;
  const Index n = model.dofs().num_dofs();
  Vector u(n);
  {
    const Matrix r = random_matrix(static_cast<int>(n), 1, 0.01, 17u);
    u = r.col(0);
  }
  for (Index d : model.dofs().constrained_dofs()) u(d) = 0.0;
  const Scalar lambda = 0.7;
  const NonlinearState state = evaluate_nonlinear_state(model, assembler, 0, options, u, lambda);
  const Matrix k(state.tangent);
  const Scalar h = 1.0e-8;
  Scalar worst = 0.0;
  for (Index d : model.dofs().free_dofs()) {
    Vector up = u;
    Vector um = u;
    up(d) += h;
    um(d) -= h;
    const Vector fd = (evaluate_nonlinear_state(model, assembler, 0, options, up, lambda).residual -
                       evaluate_nonlinear_state(model, assembler, 0, options, um, lambda).residual) /
                      (2.0 * h);
    worst = std::max(worst, (fd - k.col(d)).cwiseAbs().maxCoeff());
  }
  REQUIRE(worst <= 1.0e-6 * k.cwiseAbs().maxCoeff());
  // The follower pressure makes the tangent non-symmetric.
  REQUIRE((k - k.transpose()).cwiseAbs().maxCoeff() > 1.0e-6 * k.cwiseAbs().maxCoeff());
  // q = -dR/dlambda.
  const Scalar hl = 1.0e-6;
  const Vector dr = (evaluate_nonlinear_state(model, assembler, 0, options, u, lambda + hl).residual -
                     evaluate_nonlinear_state(model, assembler, 0, options, u, lambda - hl).residual) /
                    (2.0 * hl);
  REQUIRE((dr + state.load_rate).cwiseAbs().maxCoeff() <=
          1.0e-6 * state.load_rate.cwiseAbs().maxCoeff());
}

TEST_CASE("a homogeneous large deformation is reproduced exactly on distorted meshes",
          "[nonlinear][patch]") {
  // Every boundary node carries u = (F - I) X for a deformation gradient
  // with 30 % stretch, 20 % shear and a rotation; the interior nodes must
  // follow exactly, and the Cauchy stress must be the law's.
  for (const ElementType type : {ElementType::Quad4, ElementType::Hex8, ElementType::Tet10}) {
    StructuredMeshSpec spec;
    spec.nx = spec.ny = spec.nz = 3;
    Mesh mesh = type == ElementType::Quad4 ? make_perturbed_quad_mesh(spec, 0.25)
                : type == ElementType::Hex8 ? make_perturbed_hex_mesh(spec, 0.25)
                                            : make_perturbed_tet10_mesh(spec, 0.25);
    const int dim = mesh.dim();
    Matrix f = Matrix::Identity(dim, dim);
    f(0, 0) = 1.3;
    f(0, 1) = 0.2;
    f(1, 0) = -0.1;
    f(1, 1) = 0.9;
    if (dim == 3) {
      f(2, 2) = 1.1;
      f(2, 0) = 0.15;
    }
    const StressState state = dim == 2 ? StressState::PlaneStrain : StressState::ThreeDimensional;
    IsotropicMaterial mat(1.0e7, 0.3, 1000.0);
    for (const HyperelasticModel law :
         {HyperelasticModel::SaintVenantKirchhoff, HyperelasticModel::NeoHookean}) {
      FemModel model(mesh, mat, 1.0, state, IntegrationOptions());
      // Boundary nodes: prescribe u = (F - I) X component by component.
      std::vector<char> on_boundary(static_cast<std::size_t>(model.mesh().num_nodes()), 0);
      for (const Mesh::BoundaryFace& face : model.mesh().boundary_faces()) {
        for (Index node : face.nodes) on_boundary[static_cast<std::size_t>(node)] = 1;
      }
      for (Index node = 0; node < model.mesh().num_nodes(); ++node) {
        if (!on_boundary[static_cast<std::size_t>(node)]) continue;
        const Vector x = model.mesh().coordinates().col(node);
        const Vector target = (f - Matrix::Identity(dim, dim)) * x;
        DisplacementConstraint bc;
        Selector s;
        s.kind = SelectorKind::NodeIds;
        s.ids.assign(1, node);
        bc.region.members.push_back(s);
        for (int k = 0; k < dim; ++k) bc.set(k, true, target(k));
        model.constraints().push_back(bc);
      }
      LoadCaseSpec lc;
      lc.name = "patch";
      lc.prescribed_displacement_only = true;
      model.load_case_specs().push_back(lc);
      model.finalize();
      Assembler assembler(model);
      NonlinearOptions options;
      options.law = law;
      options.steps = 4;
      NonlinearStaticAnalysis analysis(model, assembler, options);
      const NonlinearResult r = analysis.solve(0);
      INFO(to_string(type) << " " << to_string(law));
      REQUIRE(r.completed);
      Scalar err = 0.0;
      for (Index node = 0; node < model.mesh().num_nodes(); ++node) {
        const Vector x = model.mesh().coordinates().col(node);
        const Vector exact = (f - Matrix::Identity(dim, dim)) * x;
        err = std::max(err, (r.displacement.segment(node * dim, dim) - exact).norm());
      }
      REQUIRE(err <= 1.0e-9);
      const Matrix grad = f - Matrix::Identity(dim, dim);
      const HyperelasticResponse law_response = evaluate_hyperelastic(law, mat, state, grad);
      const Vector cauchy = cauchy_from_piola_kirchhoff(grad, law_response, state);
      for (Index e = 0; e < model.mesh().num_elements(); ++e) {
        REQUIRE((r.element_cauchy.col(e) - cauchy).norm() <= 1.0e-7 * cauchy.norm());
      }
      REQUIRE(r.equilibrium.relative_force_error <= 1.0e-10);
    }
  }
}

namespace {

/// A unit cube of 2 x 2 x 2 Hex8 cells on a statically determinate (3-2-1)
/// support, under a pressure p on every face.
FemModel pressurised_cube(const IsotropicMaterial& mat, Scalar p) {
  StructuredMeshSpec spec;
  spec.nx = spec.ny = spec.nz = 2;
  FemModel model(make_structured_hex_mesh(spec), mat, 1.0, StressState::ThreeDimensional,
                 IntegrationOptions());
  model.constraints().push_back(fix(point_region(Vector3(0.0, 0.0, 0.0)), true, true, true));
  model.constraints().push_back(fix(point_region(Vector3(1.0, 0.0, 0.0)), false, true, true));
  model.constraints().push_back(fix(point_region(Vector3(0.0, 1.0, 0.0)), false, false, true));
  LoadCaseSpec lc;
  lc.name = "squeeze";
  PressureLoadSpec all;
  Selector every;
  every.kind = SelectorKind::All;
  all.region.members.push_back(every);
  all.pressure = p;
  lc.pressures.push_back(all);
  model.load_case_specs().push_back(lc);
  model.finalize();
  return model;
}

/// The stretch s of F = s I under a dead pressure p: the first
/// Piola-Kirchhoff stress -p, K3 s (s^2 - 1) / 2 = -p, by Newton from s = 1.
Scalar dead_stretch(Scalar k3, Scalar p) {
  Scalar s = 1.0;
  for (int it = 0; it < 50; ++it) {
    const Scalar g = 0.5 * k3 * s * (s * s - 1.0) + p;
    const Scalar dg = 0.5 * k3 * (3.0 * s * s - 1.0);
    s -= g / dg;
  }
  return s;
}

}  // namespace

TEST_CASE("a free cube under pressure takes the exact stretch, follower and dead",
          "[nonlinear][pressure]") {
  // F = s I under a pressure p on every face. A follower pressure is a
  // Cauchy stress -p: for Saint Venant-Kirchhoff sigma = K3 (s^2 - 1)/(2 s)
  // with K3 = 3 lambda + 2 mu, so s = (-p + sqrt(p^2 + K3^2)) / K3. A dead
  // pressure acts on the undeformed area, a first Piola-Kirchhoff stress -p:
  // K3 s (s^2 - 1) / 2 = -p. A dead compressive pressure destabilises the
  // cube early (the next test), so it is checked below its critical level.
  const Scalar e = 1.0e8;
  const Scalar nu = 0.3;
  const Scalar k3 = e / (1.0 - 2.0 * nu);
  IsotropicMaterial mat(e, nu, 1000.0);
  for (const bool follower : {true, false}) {
    const Scalar p = follower ? 0.1 * k3 : 0.004 * k3;
    FemModel model = pressurised_cube(mat, p);
    Assembler assembler(model);
    NonlinearOptions options;
    options.follower_pressure = follower;
    options.steps = 5;
    NonlinearStaticAnalysis analysis(model, assembler, options);
    const NonlinearResult r = analysis.solve(0);
    REQUIRE(r.completed);
    const Scalar s = follower ? (-p + std::sqrt(p * p + k3 * k3)) / k3 : dead_stretch(k3, p);
    INFO((follower ? "follower" : "dead") << " pressure, stretch " << s);
    for (Index node = 0; node < model.mesh().num_nodes(); ++node) {
      const Vector3 x = model.mesh().node(node);
      REQUIRE((r.displacement.segment(3 * node, 3) - (s - 1.0) * x).norm() <= 1.0e-9);
    }
    // The supports react nothing: the pressure on a closed surface balances.
    REQUIRE(r.reactions.cwiseAbs().maxCoeff() <= 1.0e-6 * p);
    // Closed surface: the follower tangent is symmetric once assembled, and
    // is factorised as such, reporting the stable path's inertia.
    REQUIRE(r.symmetric_tangent);
    for (const NonlinearStep& step : r.steps) REQUIRE(step.negative_pivots == 0);
  }
}

TEST_CASE("load control stops where a dead compressive pressure destabilises the cube",
          "[nonlinear][pressure][stability]") {
  // Under a dead pressure the loads do not turn with the body, and a
  // compressed body loses stability against rotation, which the 3-2-1
  // support resists only by the deformation near its three points. The
  // homogeneous state F = s I stays an equilibrium but turns unstable where
  // the smallest eigenvalue of the free-free tangent at that state changes
  // sign: found here by bisection on the load factor, independently of the
  // path-following, which must stop just below it and say why.
  const Scalar e = 1.0e8;
  const Scalar nu = 0.3;
  const Scalar k3 = e / (1.0 - 2.0 * nu);
  const Scalar p = 0.1 * k3;
  IsotropicMaterial mat(e, nu, 1000.0);
  FemModel model = pressurised_cube(mat, p);
  Assembler assembler(model);
  NonlinearOptions options;
  options.follower_pressure = false;
  options.steps = 5;
  const std::vector<Index>& free_dofs = model.dofs().free_dofs();
  const auto smallest_eigenvalue = [&](Scalar lambda) {
    const Scalar s = dead_stretch(k3, lambda * p);
    Vector u(model.dofs().num_dofs());
    for (Index node = 0; node < model.mesh().num_nodes(); ++node) {
      u.segment(3 * node, 3) = (s - 1.0) * model.mesh().node(node);
    }
    const NonlinearState state = evaluate_nonlinear_state(model, assembler, 0, options, u, lambda);
    const Matrix full = Matrix(state.tangent);
    const auto n = static_cast<Eigen::Index>(free_dofs.size());
    Matrix k(n, n);
    for (Eigen::Index i = 0; i < n; ++i) {
      for (Eigen::Index j = 0; j < n; ++j) {
        k(i, j) = full(free_dofs[static_cast<std::size_t>(i)], free_dofs[static_cast<std::size_t>(j)]);
      }
    }
    return Eigen::SelfAdjointEigenSolver<Matrix>(0.5 * (k + k.transpose())).eigenvalues()(0);
  };
  Scalar lo = 0.0;
  Scalar hi = 1.0;
  REQUIRE(smallest_eigenvalue(lo) > 0.0);
  REQUIRE(smallest_eigenvalue(hi) < 0.0);
  for (int it = 0; it < 60; ++it) {
    const Scalar mid = 0.5 * (lo + hi);
    (smallest_eigenvalue(mid) > 0.0 ? lo : hi) = mid;
  }
  const Scalar critical = 0.5 * (lo + hi);
  INFO("critical load factor " << critical);

  NonlinearStaticAnalysis analysis(model, assembler, options);
  const NonlinearResult r = analysis.solve(0);
  REQUIRE_FALSE(r.completed);
  REQUIRE(r.termination.find("positive definiteness") != std::string::npos);
  // The run brackets the critical point, to within 1 %.
  REQUIRE(r.load_factor <= critical);
  REQUIRE(r.critical_bound >= critical);
  REQUIRE(r.critical_bound - r.load_factor <= 0.01 * critical);
  // Every recorded state is the homogeneous one, and stable.
  for (const NonlinearStep& step : r.steps) REQUIRE(step.negative_pivots == 0);
  const Scalar s = dead_stretch(k3, r.load_factor * p);
  for (Index node = 0; node < model.mesh().num_nodes(); ++node) {
    const Vector3 x = model.mesh().node(node);
    REQUIRE((r.displacement.segment(3 * node, 3) - (s - 1.0) * x).norm() <= 1.0e-9);
  }
}

TEST_CASE("load control and arc-length reach the same large-deflection state",
          "[nonlinear][solver]") {
  StructuredMeshSpec spec;
  spec.nx = 8;
  spec.ny = 1;
  spec.nz = 1;
  spec.lx = 1.0;
  spec.ly = 0.05;
  spec.lz = 0.05;
  IsotropicMaterial mat(200.0e9, 0.3, 7800.0);
  const Scalar inertia = std::pow(0.05, 4) / 12.0;
  const Scalar load = 3.0 * 200.0e9 * inertia;  // P L^2 / EI = 3
  Vector reference;
  for (const NonlinearOptions::Method method :
       {NonlinearOptions::Method::LoadControl, NonlinearOptions::Method::ArcLength}) {
    FemModel model(make_structured_tet10_mesh(spec), mat, 1.0, StressState::ThreeDimensional,
                   IntegrationOptions());
    model.constraints().push_back(fix(box_region(-kInf, 0.0), true, true, true));
    LoadCaseSpec lc;
    lc.name = "tip";
    PointLoadSpec tip;
    tip.region = box_region(1.0, kInf);
    tip.force = Vector3(0.0, -load, 0.0);
    lc.point_loads.push_back(tip);
    model.load_case_specs().push_back(lc);
    model.finalize();
    Assembler assembler(model);
    NonlinearOptions options;
    options.method = method;
    options.steps = 8;
    options.residual_tolerance = 1.0e-10;
    NonlinearStaticAnalysis analysis(model, assembler, options);
    const NonlinearResult r = analysis.solve(0);
    INFO(to_string(method));
    REQUIRE(r.completed);
    REQUIRE(r.load_factor == Approx(1.0).margin(1e-12));
    REQUIRE(r.equilibrium.relative_force_error <= 1.0e-9);
    REQUIRE(r.equilibrium.relative_moment_error <= 1.0e-9);
    // A large deflection: the tip moves by more than a third of the span.
    if (reference.size() == 0) {
      reference = r.displacement;
      REQUIRE(r.displacement.cwiseAbs().maxCoeff() > 0.3);
    } else {
      REQUIRE((r.displacement - reference).cwiseAbs().maxCoeff() <=
              1.0e-7 * reference.cwiseAbs().maxCoeff());
    }
  }
}

TEST_CASE("the non-linear solution tends to the linear one as the load vanishes",
          "[nonlinear][solver]") {
  // A cantilever under a tip load P. Its deviation from the linear solution
  // has two parts with different orders: the axial shortening of the bent
  // beam (int theta^2 / 2) is second order in P but absent from linear
  // theory, so relative to the O(P) displacement it is O(P); the change of
  // the tip's deflection itself is relatively O(P^2) (the elastica's
  // 1 - c k^2, k = P L^2 / EI). Both rates must show, down to a load where
  // the second is 2e-8 of the deflection - which needs the laws' full
  // precision at small strain, a tolerance of 1e-12 and a solver that knows
  // its round-off floor (about 2e-11 of the deflection for this beam).
  std::vector<Scalar> field;
  std::vector<Scalar> tip;
  for (const Scalar load : {-3000.0, -300.0, -30.0}) {
    CantileverCase c;
    c.poisson = 0.3;
    c.tip_load = load;
    FemModel model = make_cantilever(c, 20, 4);
    Assembler assembler(model);
    StaticAnalysis linear(model, assembler);
    const Vector ul = linear.solve_all().front().displacement;
    NonlinearOptions options;
    options.steps = 1;
    options.residual_tolerance = 1.0e-12;
    options.displacement_tolerance = 1.0e-12;
    NonlinearStaticAnalysis analysis(model, assembler, options);
    const NonlinearResult r = analysis.solve(0);
    REQUIRE(r.completed);
    field.push_back((r.displacement - ul).norm() / ul.norm());
    tip.push_back(std::abs(tip_deflection(model, r.displacement) / tip_deflection(model, ul) - 1.0));
  }
  REQUIRE(field[1] / field[0] == Approx(0.1).epsilon(0.02));
  REQUIRE(field[2] / field[1] == Approx(0.1).epsilon(0.02));
  REQUIRE(tip[1] / tip[0] == Approx(0.01).epsilon(0.03));
  REQUIRE(tip[2] / tip[1] == Approx(0.01).epsilon(0.03));
}

TEST_CASE("the arc-length method refuses prescribed displacements; options are validated",
          "[nonlinear][solver]") {
  StructuredMeshSpec spec;
  spec.nx = 4;
  spec.ny = 2;
  FemModel model(make_structured_quad_mesh(spec), default_material(), 0.01,
                 StressState::PlaneStress, IntegrationOptions());
  model.constraints().push_back(fix(box_region(-kInf, 0.0), true, true, false));
  DisplacementConstraint pull;
  pull.region = box_region(1.0, kInf);
  pull.set(0, true, 1.0e-3);
  model.constraints().push_back(pull);
  LoadCaseSpec lc;
  lc.name = "pull";
  lc.prescribed_displacement_only = true;
  model.load_case_specs().push_back(lc);
  model.finalize();
  Assembler assembler(model);
  NonlinearOptions options;
  options.method = NonlinearOptions::Method::ArcLength;
  NonlinearStaticAnalysis analysis(model, assembler, options);
  REQUIRE_THROWS_AS(analysis.solve(0), ConfigError);
  NonlinearOptions bad;
  bad.steps = 0;
  REQUIRE_THROWS_AS(NonlinearStaticAnalysis(model, assembler, bad), ConfigError);
  REQUIRE(parse_nonlinear_method("arc_length") == NonlinearOptions::Method::ArcLength);
  REQUIRE_THROWS_AS(parse_nonlinear_method("riks"), ConfigError);
}

TEST_CASE("the nonlinear block of a deck parses, validates and drives the solver",
          "[nonlinear][config]") {
  // A plane cantilever 1 m x 0.1 m under a tip force of 2 kN.
  const auto deck = [](const std::string& nonlinear) {
    return json::parse(R"({
      "mesh": { "type": "structured_quad", "nx": 20, "ny": 2, "lx": 1.0, "ly": 0.1 },
      "material": { "youngs_modulus": 70e9, "poisson_ratio": 0.3 },
      "model": { "thickness": 0.01, "stress_state": "plane_strain" },
      "boundary_conditions": [ { "fix": ["x", "y"], "region": { "box": { "xmax": 0.0 } } } ],
      "load_cases": [
        { "name": "other", "point_loads": [ { "force": [1.0, 0.0],
            "region": { "box": { "xmin": 1.0 } } } ] },
        { "name": "tip", "point_loads": [ { "force": [0.0, -2000.0],
            "region": { "box": { "xmin": 1.0 } } } ] } ],
      "nonlinear": )" + nonlinear + "}");
  };
  const Configuration config = parse_configuration(deck(R"({
      "enabled": true, "material_model": "neo_hookean", "method": "load_control",
      "steps": 4, "load_factors": [0.1, 0.35],
      "residual_tolerance": 1e-9, "max_cuts": 8, "follower_pressure": false,
      "monitors": [
        { "name": "tip_v", "component": "y", "region": { "box": { "xmin": 1.0 } } },
        { "name": "root_ry", "component": "y", "quantity": "reaction",
          "region": { "box": { "xmax": 0.0 } } } ],
      "load_cases": ["tip"] })"),
                                                   "inline", true);
  REQUIRE(config.nonlinear.enabled);
  const NonlinearOptions& o = config.nonlinear.options;
  REQUIRE(o.law == HyperelasticModel::NeoHookean);
  REQUIRE(o.method == NonlinearOptions::Method::LoadControl);
  REQUIRE(o.steps == 4);
  REQUIRE(o.max_cuts == 8);
  REQUIRE(o.residual_tolerance == 1.0e-9);
  REQUIRE_FALSE(o.follower_pressure);
  REQUIRE(o.load_factors == std::vector<Scalar>{0.1, 0.35});
  REQUIRE(o.monitors.size() == 2);
  REQUIRE(o.monitors[1].quantity == NonlinearMonitor::Quantity::Reaction);
  REQUIRE(config.nonlinear_load_cases() == std::vector<std::size_t>{1});

  // The run passes through the stations and the reaction monitor carries the
  // support force: the root balances the tip force at every step.
  FemModel model = build_model(config);
  Assembler assembler(model);
  const NonlinearResult r = NonlinearStaticAnalysis(model, assembler, o).solve(1);
  REQUIRE(r.completed);
  REQUIRE(r.monitor_units == std::vector<std::string>{"m", "N"});
  std::vector<Scalar> factors;
  for (const NonlinearStep& s : r.steps) {
    factors.push_back(s.load_factor);
    REQUIRE(s.monitors[1] == Approx(2000.0 * s.load_factor).epsilon(1.0e-8));
    REQUIRE(s.monitors[0] < 0.0);
  }
  REQUIRE(std::find(factors.begin(), factors.end(), 0.1) != factors.end());
  REQUIRE(std::find(factors.begin(), factors.end(), 0.35) != factors.end());
  REQUIRE(factors.back() == 1.0);

  // What the block refuses.
  const auto refuses = [&](const std::string& block) {
    REQUIRE_THROWS_AS(parse_configuration(deck(block), "inline", true), ConfigError);
  };
  refuses(R"({ "enabled": true, "method": "riks" })");
  refuses(R"({ "enabled": true, "material_model": "mooney_rivlin" })");
  refuses(R"({ "enabled": true, "load_cases": ["missing"] })");
  refuses(R"({ "enabled": true, "load_factors": [0.5, 0.2] })");
  refuses(R"({ "enabled": true, "load_factors": [0.5, 1.0] })");
  refuses(R"({ "enabled": true, "method": "arc_length", "load_factors": [0.5] })");
  refuses(R"({ "enabled": true, "steps": 0 })");
  refuses(R"({ "enabled": true, "min_arc_ratio": 2.0 })");
  refuses(R"({ "enabled": true, "monitors": [ { "component": "y", "quantity": "stress",
                "region": { "box": { "xmin": 1.0 } } } ] })");
  refuses(R"({ "enabled": true, "monitors": [ { "component": "z",
                "region": { "box": { "xmin": 1.0 } } } ] })");
  refuses(R"({ "enabled": true, "stepz": 4 })");  // unknown key, strict parsing
}

TEST_CASE("neo-Hookean in plane stress is refused by the deck and by the solver",
          "[nonlinear][config]") {
  const json::Value doc = json::parse(R"({
      "mesh": { "type": "structured_quad", "nx": 4, "ny": 2, "lx": 1.0, "ly": 0.1 },
      "material": { "youngs_modulus": 70e9, "poisson_ratio": 0.3 },
      "model": { "thickness": 0.01, "stress_state": "plane_stress" },
      "boundary_conditions": [ { "fix": ["x", "y"], "region": { "box": { "xmax": 0.0 } } } ],
      "load_cases": [ { "name": "tip", "point_loads": [ { "force": [0.0, -1.0],
            "region": { "box": { "xmin": 1.0 } } } ] } ],
      "nonlinear": { "enabled": true, "material_model": "neo_hookean" } })");
  REQUIRE_THROWS_AS(parse_configuration(doc, "inline", true), ConfigError);
  FemModel model = make_small_plate();
  Assembler assembler(model);
  NonlinearOptions o;
  o.law = HyperelasticModel::NeoHookean;
  REQUIRE_THROWS_AS(NonlinearStaticAnalysis(model, assembler, o), ConfigError);
}

TEST_CASE("a freely heated body takes exactly the thermal stretch 1 + alpha dT, stress-free",
          "[nonlinear][thermal]") {
  // The thermal strain of the Saint Venant-Kirchhoff law is the
  // Green-Lagrange strain of the free thermal stretch, so a body free to
  // expand by alpha dT = 5 % reaches F = (1 + alpha dT) I with no stress; an
  // additive alpha dT in the Green strain would stop at sqrt(1 + 2 alpha dT),
  // 1.2 % of the expansion short.
  const Scalar alpha = 5.0e-4;
  const Scalar t_ref = 293.15;
  const Scalar dt = 100.0;
  IsotropicMaterial mat(1.0e9, 0.3, 1000.0);
  mat.set_thermal(alpha, t_ref, 1.0);
  for (const ElementType type : {ElementType::Quad4, ElementType::Hex8, ElementType::Tet10}) {
    StructuredMeshSpec spec;
    spec.nx = spec.ny = spec.nz = 2;
    const bool solid = type != ElementType::Quad4;
    Mesh mesh = type == ElementType::Quad4  ? make_structured_quad_mesh(spec)
                : type == ElementType::Hex8 ? make_structured_hex_mesh(spec)
                                            : make_structured_tet10_mesh(spec);
    FemModel model(std::move(mesh), mat, solid ? 1.0 : 0.01,
                   solid ? StressState::ThreeDimensional : StressState::PlaneStress,
                   IntegrationOptions());
    if (solid) {
      model.constraints().push_back(fix(point_region(Vector3(0.0, 0.0, 0.0)), true, true, true));
      model.constraints().push_back(fix(point_region(Vector3(1.0, 0.0, 0.0)), false, true, true));
      model.constraints().push_back(fix(point_region(Vector3(0.0, 1.0, 0.0)), false, false, true));
    } else {
      model.constraints().push_back(fix(point_region(Vector3(0.0, 0.0, 0.0)), true, true, false));
      model.constraints().push_back(fix(point_region(Vector3(1.0, 0.0, 0.0)), false, true, false));
    }
    LoadCaseSpec lc;
    lc.name = "heat";
    lc.temperature.source = TemperatureSpec::Source::Uniform;
    lc.temperature.uniform = t_ref + dt;
    model.load_case_specs().push_back(lc);
    model.finalize();
    Assembler assembler(model);
    NonlinearOptions options;
    options.steps = 4;
    const NonlinearResult r = NonlinearStaticAnalysis(model, assembler, options).solve(0);
    INFO(to_string(type));
    REQUIRE(r.completed);
    const int dim = model.dim();
    for (Index node = 0; node < model.mesh().num_nodes(); ++node) {
      const Vector3 x = model.mesh().node(node);
      for (int k = 0; k < dim; ++k) {
        REQUIRE(r.displacement(dim * node + k) == Approx(alpha * dt * x(k)).margin(1.0e-12));
      }
    }
    REQUIRE(r.element_von_mises.cwiseAbs().maxCoeff() <= 1.0e-6 * mat.youngs_modulus() * 1.0e-6);
    REQUIRE(r.reactions.cwiseAbs().maxCoeff() <= 1.0e-6);
  }
}

TEST_CASE("the centrifugal load acts at the deformed position with the mass, not rho times it",
          "[nonlinear][rotation]") {
  // Moving a body by d normal to the axis moves every point's distance from
  // it by d, so the total centrifugal force changes by rho omega^2 V d
  // exactly - whatever the element, the tangent's spin-softening term
  // omega^2 M_perp must carry the density once.
  const Scalar rho = 1100.0;
  const Scalar omega = 50.0;
  IsotropicMaterial mat(1.0e9, 0.3, rho);
  for (const ElementType type : {ElementType::Quad4, ElementType::Hex8, ElementType::Tet10}) {
    StructuredMeshSpec spec;
    spec.nx = spec.ny = spec.nz = 2;
    spec.lx = 0.5;
    spec.ly = 0.4;
    spec.lz = 0.3;
    const bool solid = type != ElementType::Quad4;
    Mesh mesh = type == ElementType::Quad4  ? make_structured_quad_mesh(spec)
                : type == ElementType::Hex8 ? make_structured_hex_mesh(spec)
                                            : make_structured_tet10_mesh(spec);
    const Scalar thickness = solid ? 1.0 : 0.02;
    FemModel model(std::move(mesh), mat, thickness,
                   solid ? StressState::ThreeDimensional : StressState::PlaneStress,
                   IntegrationOptions());
    model.constraints().push_back(
        fix(point_region(Vector3(0.0, 0.0, 0.0)), true, true, solid));
    LoadCaseSpec lc;
    lc.name = "spin";
    lc.centrifugal.enabled = true;
    lc.centrifugal.angular_velocity = omega;
    lc.centrifugal.axis = Vector3::UnitZ();
    lc.centrifugal.point = Vector3(-0.2, 0.1, 0.0);
    model.load_case_specs().push_back(lc);
    model.finalize();
    Assembler assembler(model);
    const int dim = model.dim();
    const Scalar volume = solid ? 0.5 * 0.4 * 0.3 : 0.5 * 0.4 * thickness;
    const Scalar d = 1.0e-3;
    Vector u0 = Vector::Zero(model.dofs().num_dofs());
    Vector u1 = u0;
    for (Index node = 0; node < model.mesh().num_nodes(); ++node) u1(dim * node) = d;
    NonlinearOptions options;
    const NonlinearState s0 = evaluate_nonlinear_state(model, assembler, 0, options, u0, 1.0);
    const NonlinearState s1 = evaluate_nonlinear_state(model, assembler, 0, options, u1, 1.0);
    Vector3 change = Vector3::Zero();
    for (Index node = 0; node < model.mesh().num_nodes(); ++node) {
      for (int k = 0; k < dim; ++k) change(k) += s1.external(dim * node + k) - s0.external(dim * node + k);
    }
    INFO(to_string(type));
    REQUIRE(change.x() == Approx(rho * omega * omega * volume * d).epsilon(1.0e-12));
    REQUIRE(std::abs(change.y()) <= 1.0e-12 * rho * omega * omega * volume * d);
    // The tangent's rotation term is the same matrix: K_T u1 - K_T u0 over a
    // rigid translation is -omega^2 M_perp u1 (the material part of a rigid
    // motion vanishes).
    const Vector ku = s1.tangent * u1;
    Vector3 tangent_change = Vector3::Zero();
    for (Index node = 0; node < model.mesh().num_nodes(); ++node) {
      for (int k = 0; k < dim; ++k) tangent_change(k) += ku(dim * node + k);
    }
    REQUIRE(tangent_change.x() == Approx(-rho * omega * omega * volume * d).epsilon(1.0e-10));
  }
}

TEST_CASE("a heated cube held between two walls takes the multiplicative thermoelastic state",
          "[nonlinear][thermal]") {
  // Held at u_x = 0 on x = 0 and x = 1, free sideways, heated by dT: F =
  // diag(1, l, l). With theta = 1 + alpha dT and E_theta = (theta^2 - 1)/2,
  // S = D (E - E_theta)/theta gives S_yy = 0 at E_yy = (1 + nu) E_theta, so
  // l = sqrt(1 + 2 (1 + nu) E_theta), and S_xx = -Young E_theta / theta, a
  // Cauchy stress sigma_xx = S_xx / l^2. The additive split (no 1/theta)
  // would give a stress 5 % higher at alpha dT = 0.05.
  const Scalar young = 1.0e9;
  const Scalar nu = 0.3;
  const Scalar alpha = 5.0e-4;
  const Scalar dt = 100.0;
  IsotropicMaterial mat(young, nu, 1000.0);
  mat.set_thermal(alpha, 0.0, 1.0);
  for (const ElementType type : {ElementType::Hex8, ElementType::Tet10}) {
    StructuredMeshSpec spec;
    spec.nx = spec.ny = spec.nz = 1;
    FemModel model(type == ElementType::Hex8 ? make_structured_hex_mesh(spec)
                                             : make_structured_tet10_mesh(spec),
                   mat, 1.0, StressState::ThreeDimensional, IntegrationOptions());
    model.constraints().push_back(fix(box_region(-kInf, 0.0), true, false, false));
    model.constraints().push_back(fix(box_region(1.0, kInf), true, false, false));
    model.constraints().push_back(fix(point_region(Vector3(0.0, 0.0, 0.0)), false, true, true));
    model.constraints().push_back(fix(point_region(Vector3(0.0, 1.0, 0.0)), false, false, true));
    model.constraints().push_back(fix(point_region(Vector3(0.0, 0.0, 1.0)), false, true, false));
    LoadCaseSpec lc;
    lc.name = "heat";
    lc.temperature.source = TemperatureSpec::Source::Uniform;
    lc.temperature.uniform = dt;
    model.load_case_specs().push_back(lc);
    model.finalize();
    Assembler assembler(model);
    NonlinearOptions options;
    options.steps = 4;
    const NonlinearResult r = NonlinearStaticAnalysis(model, assembler, options).solve(0);
    INFO(to_string(type));
    REQUIRE(r.completed);
    const Scalar theta = 1.0 + alpha * dt;
    const Scalar e_theta = 0.5 * (theta * theta - 1.0);
    const Scalar l = std::sqrt(1.0 + 2.0 * (1.0 + nu) * e_theta);
    const Scalar sigma_xx = -young * e_theta / (theta * l * l);
    for (Index node = 0; node < model.mesh().num_nodes(); ++node) {
      const Vector3 x = model.mesh().node(node);
      REQUIRE(r.displacement(3 * node) == Approx(0.0).margin(1.0e-12));
      REQUIRE(r.displacement(3 * node + 1) == Approx((l - 1.0) * x.y()).margin(1.0e-12));
      REQUIRE(r.displacement(3 * node + 2) == Approx((l - 1.0) * x.z()).margin(1.0e-12));
    }
    for (Index e = 0; e < model.mesh().num_elements(); ++e) {
      REQUIRE(r.element_cauchy(0, e) == Approx(sigma_xx).epsilon(1.0e-10));
      REQUIRE(std::abs(r.element_cauchy(1, e)) <= 1.0e-9 * std::abs(sigma_xx));
      REQUIRE(std::abs(r.element_cauchy(2, e)) <= 1.0e-9 * std::abs(sigma_xx));
    }
  }
}
