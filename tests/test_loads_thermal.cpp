/// \file test_loads_thermal.cpp
/// \brief Self-weight, body forces, rotation, pressure, thermal strain, steady
///        conduction and multi-material models against exact solutions.
///
/// Every reference here is exact, not a convergence target: one-dimensional
/// problems whose finite-element nodal values are exact for consistent loads,
/// fields that an element represents exactly (a quadratic displacement on the
/// Tet10), resultants that equal a mass or a projected area, and states of
/// uniform stress.
#include "TestSupport.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/elements/FaceGeometry.hpp"
#include "sparlab/fem/Buckling.hpp"
#include "sparlab/fem/HeatConduction.hpp"
#include "sparlab/fem/Loads.hpp"
#include "sparlab/fem/StressRecovery.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <limits>

using namespace sparlab;
using namespace sparlab::testing;
using Catch::Approx;

namespace {

constexpr Scalar kInf = std::numeric_limits<Scalar>::infinity();

StructuredMeshSpec spec2(Index nx, Index ny, Scalar lx, Scalar ly) {
  StructuredMeshSpec s;
  s.nx = nx;
  s.ny = ny;
  s.lx = lx;
  s.ly = ly;
  return s;
}

StructuredMeshSpec spec3(Index nx, Index ny, Index nz, Scalar lx, Scalar ly, Scalar lz) {
  StructuredMeshSpec s = spec2(nx, ny, lx, ly);
  s.nz = nz;
  s.lz = lz;
  return s;
}

SelectorGroup box(Scalar xmin, Scalar xmax, Scalar ymin = -kInf, Scalar ymax = kInf,
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

SelectorGroup nearest(const Vector3& x) {
  SelectorGroup g;
  Selector s;
  s.kind = SelectorKind::NearestNode;
  s.point = x;
  g.members.push_back(s);
  g.name = "point";
  return g;
}

DisplacementConstraint fixed(const SelectorGroup& region, bool x, bool y, bool z = false) {
  DisplacementConstraint bc;
  bc.region = region;
  bc.fix_x = x;
  bc.fix_y = y;
  bc.fix_z = z;
  return bc;
}

IsotropicMaterial thermal_steel() {
  IsotropicMaterial m(200.0e9, 0.3, 7800.0, "steel");
  m.set_thermal(1.2e-5, 20.0, 50.0);
  return m;
}

struct Solved {
  std::vector<StaticSolution> solutions;
  std::vector<StressField> stresses;
};

Solved solve(FemModel& model) {
  model.finalize();
  Assembler assembler(model);
  StaticAnalysis analysis(model, assembler);
  Solved out;
  out.solutions = analysis.solve_all();
  for (std::size_t l = 0; l < out.solutions.size(); ++l) {
    const Vector& t = model.load_case_data(l).temperature;
    out.stresses.push_back(recover_stresses(model, assembler, out.solutions[l].displacement,
                                            nullptr, t.size() > 0 ? &t : nullptr));
  }
  return out;
}

LoadCaseSpec uniform_temperature(Scalar t) {
  LoadCaseSpec lc;
  lc.name = "heat";
  lc.temperature.source = TemperatureSpec::Source::Uniform;
  lc.temperature.uniform = t;
  return lc;
}

}  // namespace

// --- body loads -------------------------------------------------------------

TEST_CASE("self-weight sums to the model's mass times g on every element type",
          "[loads][body]") {
  const Vector3 g(0.0, -9.81, 0.0);
  const Vector3 g3(1.5, -9.81, 2.0);
  struct Case {
    const char* name;
    Mesh mesh;
    StressState state;
    Scalar thickness;
    Vector3 gravity;
  };
  std::vector<Case> cases;
  cases.push_back({"Q4", make_perturbed_quad_mesh(spec2(5, 3, 2.0, 1.0), 0.2, 3u),
                   StressState::PlaneStress, 0.02, g});
  cases.push_back({"Tri3", make_perturbed_tri_mesh(spec2(5, 3, 2.0, 1.0), 0.2, 3u),
                   StressState::PlaneStrain, 0.5, g});
  cases.push_back({"Hex8", make_perturbed_hex_mesh(spec3(3, 2, 2, 1.0, 0.6, 0.5), 0.2, 3u),
                   StressState::ThreeDimensional, 1.0, g3});
  cases.push_back({"Tet4", make_perturbed_tet_mesh(spec3(3, 2, 2, 1.0, 0.6, 0.5), 0.2, 3u),
                   StressState::ThreeDimensional, 1.0, g3});
  cases.push_back({"Tet10", make_perturbed_tet10_mesh(spec3(2, 2, 2, 1.0, 0.6, 0.5), 0.2, 3u),
                   StressState::ThreeDimensional, 1.0, g3});
  for (Case& c : cases) {
    INFO(c.name);
    FemModel model(c.mesh, IsotropicMaterial(70.0e9, 0.3, 2700.0), c.thickness, c.state,
                   IntegrationOptions());
    model.constraints().push_back(fixed(box(-kInf, 0.0), true, true, c.mesh.dim() == 3));
    LoadCaseSpec lc;
    lc.name = "weight";
    lc.gravity = c.gravity;
    model.load_case_specs().push_back(lc);
    model.finalize();
    Assembler assembler(model);
    const Scalar mass = assembler.total_mass();
    const Vector3 resultant = model.load_case_data(0).body_resultant;
    for (int k = 0; k < c.mesh.dim(); ++k) {
      REQUIRE(resultant(k) == Approx(mass * c.gravity(k)).epsilon(1e-12));
    }
  }
}

TEST_CASE("a hanging bar under self-weight is exact at the nodes (Q4) and everywhere (Tet10)",
          "[loads][body][verification]") {
  const Scalar length = 2.0;
  const Scalar e = 70.0e9;
  const Scalar rho = 2700.0;
  const Scalar g = 9.81;
  // Bar along x, fixed at x = 0, gravity along +x, Poisson ratio 0 so the
  // problem is one-dimensional: u(x) = rho g / E (L x - x^2 / 2).
  const auto exact = [&](Scalar x) { return rho * g / e * (length * x - 0.5 * x * x); };
  const Scalar scale = exact(length);

  SECTION("Q4, nodally exact") {
    FemModel model(make_structured_quad_mesh(spec2(10, 2, length, 0.1)),
                   IsotropicMaterial(e, 0.0, rho), 0.01, StressState::PlaneStress,
                   IntegrationOptions());
    model.constraints().push_back(fixed(box(-kInf, 0.0), true, false));
    model.constraints().push_back(fixed(nearest(Vector3::Zero()), false, true));
    LoadCaseSpec lc;
    lc.name = "weight";
    lc.gravity = Vector3(g, 0.0, 0.0);
    model.load_case_specs().push_back(lc);
    const Solved s = solve(model);
    const Vector& u = s.solutions[0].displacement;
    for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
      const Scalar x = model.mesh().node(n).x();
      REQUIRE(std::abs(u(2 * n) - exact(x)) < 1e-12 * scale);
      REQUIRE(std::abs(u(2 * n + 1)) < 1e-12 * scale);
    }
    // The reaction carries the whole weight.
    REQUIRE(s.solutions[0].equilibrium.reaction_force.x() ==
            Approx(-rho * g * length * 0.1 * 0.01).epsilon(1e-12));
  }
  SECTION("Tet10, exact everywhere") {
    FemModel model(make_structured_tet10_mesh(spec3(3, 1, 1, length, 0.1, 0.1)),
                   IsotropicMaterial(e, 0.0, rho), 1.0, StressState::ThreeDimensional,
                   IntegrationOptions());
    model.constraints().push_back(fixed(box(-kInf, 0.0), true, true, true));
    LoadCaseSpec lc;
    lc.name = "weight";
    lc.gravity = Vector3(g, 0.0, 0.0);
    model.load_case_specs().push_back(lc);
    const Solved s = solve(model);
    const Vector& u = s.solutions[0].displacement;
    for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
      const Scalar x = model.mesh().node(n).x();
      REQUIRE(std::abs(u(3 * n) - exact(x)) < 1e-11 * scale);
    }
  }
}

TEST_CASE("a rotating bar carries the exact centrifugal stretch", "[loads][body][verification]") {
  // A bar along x in [0, L] x [0, h] turning about the z axis through the
  // origin: b = rho w^2 (x, y). With nu = 0 the solution separates,
  // u_x = rho w^2 / E (L^2 x / 2 - x^3 / 6), u_y = rho w^2 / E (h^2 y / 2 - y^3 / 6),
  // and the bilinear element is exact at its nodes.
  const Scalar length = 1.0;
  const Scalar h = 0.2;
  const Scalar e = 200.0e9;
  const Scalar rho = 7800.0;
  const Scalar w = 300.0;
  FemModel model(make_structured_quad_mesh(spec2(8, 3, length, h)), IsotropicMaterial(e, 0.0, rho),
                 0.01, StressState::PlaneStress, IntegrationOptions());
  model.constraints().push_back(fixed(box(-kInf, 0.0), true, false));
  model.constraints().push_back(fixed(box(-kInf, kInf, -kInf, 0.0), false, true));
  LoadCaseSpec lc;
  lc.name = "spin";
  lc.centrifugal.enabled = true;
  lc.centrifugal.angular_velocity = w;
  model.load_case_specs().push_back(lc);
  const Solved s = solve(model);
  const Vector& u = s.solutions[0].displacement;
  const Scalar c = rho * w * w / e;
  const Scalar scale = c * length * length * length / 3.0;
  for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
    const Vector3 x = model.mesh().node(n);
    REQUIRE(std::abs(u(2 * n) - c * (0.5 * length * length * x.x() - x.x() * x.x() * x.x() / 6.0)) <
            1e-12 * scale);
    REQUIRE(std::abs(u(2 * n + 1) - c * (0.5 * h * h * x.y() - x.y() * x.y() * x.y() / 6.0)) <
            1e-12 * scale);
  }
}

TEST_CASE("a plane model rejects an in-plane rotation axis and massless gravity", "[loads][body]") {
  {
    FemModel model(make_structured_quad_mesh(spec2(2, 2, 1.0, 1.0)), default_material(), 0.01,
                   StressState::PlaneStress, IntegrationOptions());
    model.constraints().push_back(fixed(box(-kInf, 0.0), true, true));
    LoadCaseSpec lc;
    lc.name = "bad_axis";
    lc.centrifugal.enabled = true;
    lc.centrifugal.angular_velocity = 10.0;
    lc.centrifugal.axis = Vector3::UnitX();
    model.load_case_specs().push_back(lc);
    REQUIRE_THROWS_AS(model.finalize(), ConfigError);
  }
  {
    FemModel model(make_structured_quad_mesh(spec2(2, 2, 1.0, 1.0)),
                   IsotropicMaterial(70.0e9, 0.3, 0.0), 0.01, StressState::PlaneStress,
                   IntegrationOptions());
    model.constraints().push_back(fixed(box(-kInf, 0.0), true, true));
    LoadCaseSpec lc;
    lc.name = "weightless";
    lc.gravity = Vector3(0.0, -9.81, 0.0);
    model.load_case_specs().push_back(lc);
    REQUIRE_THROWS_AS(model.finalize(), ConfigError);
  }
}

// --- pressure ---------------------------------------------------------------

TEST_CASE("pressure on a curved Tet10 surface has the projected-area resultant",
          "[loads][pressure]") {
  // A quarter of a thick cylinder, r in [a, b], theta in [0, pi/2], z in [0, L],
  // built by mapping every node of a structured Tet10 box - edge nodes
  // included - through (r, theta, z) -> (r cos theta, r sin theta, z), so the
  // outer face is curved. The pressure on it pushes inward; its resultant is
  // -p times the projected area b L in both x and y, whatever the faceting.
  const Scalar a = 0.5;
  const Scalar b = 1.0;
  const Scalar len = 0.4;
  const Scalar pi = 3.14159265358979323846;
  Mesh box_mesh = make_structured_tet10_mesh(spec3(2, 4, 2, b - a, 0.5 * pi, len));
  Matrix coords = box_mesh.coordinates();
  for (Index n = 0; n < coords.cols(); ++n) {
    const Scalar r = a + coords(0, n);
    const Scalar th = coords(1, n);
    coords(0, n) = r * std::cos(th);
    coords(1, n) = r * std::sin(th);
  }
  Mesh mesh(coords, box_mesh.connectivity(), ElementType::Tet10);
  FemModel model(mesh, default_material(), 1.0, StressState::ThreeDimensional,
                 IntegrationOptions());
  model.constraints().push_back(fixed(box(-kInf, 1e-9, -kInf, kInf), true, false, false));
  model.constraints().push_back(fixed(box(-kInf, kInf, -kInf, 1e-9), false, true, false));
  model.constraints().push_back(fixed(box(-kInf, kInf, -kInf, kInf, -kInf, 0.0), false, false, true));
  LoadCaseSpec lc;
  lc.name = "outer_pressure";
  PressureLoadSpec p;
  Selector outer;
  outer.kind = SelectorKind::Annulus;
  outer.center = Vector3::Zero();
  outer.inner_radius = b - 1e-9;
  outer.radius = b + 1e-9;
  p.region.members.push_back(outer);
  p.region.name = "outer";
  p.pressure = 2.0e6;
  lc.pressures.push_back(p);
  model.load_case_specs().push_back(lc);
  model.finalize();
  const Vector& f = model.load_case_data(0).mechanical;
  Vector3 total = Vector3::Zero();
  for (Index n = 0; n < mesh.num_nodes(); ++n) total += Vector3(f(3 * n), f(3 * n + 1), f(3 * n + 2));
  REQUIRE(total.x() == Approx(-2.0e6 * b * len).epsilon(1e-12));
  REQUIRE(total.y() == Approx(-2.0e6 * b * len).epsilon(1e-12));
  REQUIRE(std::abs(total.z()) < 1e-9 * 2.0e6 * b * len);
}

TEST_CASE("the follower-pressure stiffness is the derivative of the pressure load",
          "[loads][pressure]") {
  // d f_a / d x_b of face_pressure_forces against central differences, on
  // distorted (non-planar, curved) faces of every shape: the stiffness a
  // large-deflection analysis needs when the pressure follows the face.
  struct Case {
    FaceShape shape;
    Matrix x;
  };
  std::vector<Case> cases;
  {
    Matrix x(2, 2);
    x << 0.1, 0.9, -0.2, 0.3;
    cases.push_back({FaceShape::Line2, x});
  }
  {
    Matrix x(3, 3);
    x << 0.0, 1.1, 0.2, 0.1, -0.1, 0.9, 0.05, 0.2, -0.1;
    cases.push_back({FaceShape::Tri3, x});
  }
  {
    Matrix x(3, 4);
    x << 0.0, 1.0, 1.2, -0.1, 0.0, 0.1, 0.9, 1.0, 0.0, 0.15, -0.1, 0.2;
    cases.push_back({FaceShape::Quad4, x});
  }
  {
    // Corners, then edge nodes pushed off the chords: a curved face.
    Matrix x(3, 6);
    x << 0.0, 1.0, 0.0, 0.5, 0.55, 0.0, 0.0, 0.0, 1.0, 0.05, 0.5, 0.45, 0.0, 0.1, -0.05,
        0.12, 0.2, 0.1;
    cases.push_back({FaceShape::Tri6, x});
  }
  const Scalar p = 3.0e5;
  const Scalar thickness = 0.02;
  for (const Case& c : cases) {
    const int nf = face_shape_nodes(c.shape);
    const int dim = static_cast<int>(c.x.rows());
    const Matrix k = face_pressure_stiffness(c.shape, c.x, p, thickness, 4);
    Scalar scale = 0.0;
    Scalar worst = 0.0;
    const Scalar h = 1.0e-6;
    for (int b = 0; b < nf; ++b) {
      for (int j = 0; j < dim; ++j) {
        Matrix plus = c.x;
        Matrix minus = c.x;
        plus(j, b) += h;
        minus(j, b) -= h;
        const Matrix fd = (face_pressure_forces(c.shape, plus, p, thickness, 4) -
                           face_pressure_forces(c.shape, minus, p, thickness, 4)) /
                          (2.0 * h);
        for (int a = 0; a < nf; ++a) {
          for (int i = 0; i < 3; ++i) {
            const Scalar analytic = k(3 * a + i, 3 * b + j);
            scale = std::max(scale, std::abs(analytic));
            worst = std::max(worst, std::abs(analytic - fd(i, a)));
          }
        }
      }
    }
    INFO("face shape " << static_cast<int>(c.shape));
    REQUIRE(scale > 0.0);
    REQUIRE(worst <= 1.0e-7 * scale);
  }
}

// --- thermal strain ---------------------------------------------------------

TEST_CASE("free thermal expansion is stress-free on every element and idealisation",
          "[thermal][verification]") {
  const Scalar dt = 100.0;
  const IsotropicMaterial steel = thermal_steel();
  const Scalar alpha = steel.thermal_expansion();
  const Scalar sigma_scale = steel.youngs_modulus() * alpha * dt;
  struct Case {
    const char* name;
    Mesh mesh;
    StressState state;
  };
  std::vector<Case> cases;
  cases.push_back({"Q4 plane stress", make_perturbed_quad_mesh(spec2(4, 3, 2.0, 1.0), 0.2, 5u),
                   StressState::PlaneStress});
  cases.push_back({"Q4 plane strain", make_perturbed_quad_mesh(spec2(4, 3, 2.0, 1.0), 0.2, 5u),
                   StressState::PlaneStrain});
  cases.push_back({"Tri3 plane stress", make_perturbed_tri_mesh(spec2(4, 3, 2.0, 1.0), 0.2, 5u),
                   StressState::PlaneStress});
  cases.push_back({"Hex8", make_perturbed_hex_mesh(spec3(3, 2, 2, 2.0, 1.0, 0.5), 0.2, 5u),
                   StressState::ThreeDimensional});
  cases.push_back({"Tet4", make_perturbed_tet_mesh(spec3(3, 2, 2, 2.0, 1.0, 0.5), 0.2, 5u),
                   StressState::ThreeDimensional});
  cases.push_back({"Tet10", make_perturbed_tet10_mesh(spec3(2, 2, 2, 2.0, 1.0, 0.5), 0.2, 5u),
                   StressState::ThreeDimensional});
  for (Case& c : cases) {
    INFO(c.name);
    const int dim = c.mesh.dim();
    FemModel model(c.mesh, steel, dim == 2 ? 0.05 : 1.0, c.state, IntegrationOptions());
    // Statically determinate 3-2-1 supports: no thermal stress may arise.
    model.constraints().push_back(fixed(nearest(Vector3::Zero()), true, true, dim == 3));
    model.constraints().push_back(fixed(nearest(Vector3(2.0, 0.0, 0.0)), false, true, dim == 3));
    if (dim == 3) model.constraints().push_back(fixed(nearest(Vector3(0.0, 1.0, 0.0)), false, false, true));
    model.load_case_specs().push_back(uniform_temperature(steel.reference_temperature() + dt));
    const Solved s = solve(model);
    const Scalar stretch =
        (c.state == StressState::PlaneStrain ? 1.0 + steel.poisson_ratio() : 1.0) * alpha * dt;
    const Vector& u = s.solutions[0].displacement;
    for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
      const Vector3 x = model.mesh().node(n);
      for (int k = 0; k < dim; ++k) REQUIRE(std::abs(u(dim * n + k) - stretch * x(k)) < 1e-12 * 2.0);
    }
    REQUIRE(s.stresses[0].element_stress.cwiseAbs().maxCoeff() < 1e-8 * sigma_scale);
    REQUIRE(s.solutions[0].strain_energy < 1e-10 * s.solutions[0].compliance + 1e-30);
    if (c.state == StressState::PlaneStrain) {
      // The restrained out-of-plane expansion leaves sigma_zz = -E alpha dT.
      for (Index e = 0; e < model.mesh().num_elements(); ++e) {
        REQUIRE(s.stresses[0].element_sigma_zz(e) == Approx(-sigma_scale).epsilon(1e-9));
      }
    }
  }
}

TEST_CASE("a fully restrained block carries the hydrostatic thermal stress",
          "[thermal][verification]") {
  const Scalar dt = -60.0;
  const IsotropicMaterial steel = thermal_steel();
  const Scalar nu = steel.poisson_ratio();
  const Scalar e = steel.youngs_modulus();
  const Scalar alpha = steel.thermal_expansion();
  SECTION("Hex8: sigma = -E alpha dT / (1 - 2 nu) in every direction") {
    FemModel model(make_structured_hex_mesh(spec3(2, 2, 2, 1.0, 1.0, 1.0)), steel, 1.0,
                   StressState::ThreeDimensional, IntegrationOptions());
    model.constraints().push_back(fixed(box(-kInf, 0.0), true, false, false));
    model.constraints().push_back(fixed(box(1.0, kInf), true, false, false));
    model.constraints().push_back(fixed(box(-kInf, kInf, -kInf, 0.0), false, true, false));
    model.constraints().push_back(fixed(box(-kInf, kInf, 1.0, kInf), false, true, false));
    model.constraints().push_back(fixed(box(-kInf, kInf, -kInf, kInf, -kInf, 0.0), false, false, true));
    model.constraints().push_back(fixed(box(-kInf, kInf, -kInf, kInf, 1.0, kInf), false, false, true));
    model.load_case_specs().push_back(uniform_temperature(steel.reference_temperature() + dt));
    const Solved s = solve(model);
    const Scalar expected = -e * alpha * dt / (1.0 - 2.0 * nu);
    for (Index el = 0; el < model.mesh().num_elements(); ++el) {
      for (int k = 0; k < 3; ++k) {
        REQUIRE(s.stresses[0].element_stress(k, el) == Approx(expected).epsilon(1e-10));
      }
    }
    // Elastic energy of the restrained state: 3/2 sigma eps_mech V.
    const Scalar eps_mech = -alpha * dt;
    REQUIRE(s.solutions[0].strain_energy == Approx(1.5 * expected * eps_mech).epsilon(1e-9));
  }
  SECTION("Q4 plane stress: sigma = -E alpha dT / (1 - nu)") {
    FemModel model(make_structured_quad_mesh(spec2(3, 2, 1.0, 1.0)), steel, 0.1,
                   StressState::PlaneStress, IntegrationOptions());
    model.constraints().push_back(fixed(box(-kInf, 0.0), true, false));
    model.constraints().push_back(fixed(box(1.0, kInf), true, false));
    model.constraints().push_back(fixed(box(-kInf, kInf, -kInf, 0.0), false, true));
    model.constraints().push_back(fixed(box(-kInf, kInf, 1.0, kInf), false, true));
    model.load_case_specs().push_back(uniform_temperature(steel.reference_temperature() + dt));
    const Solved s = solve(model);
    for (Index el = 0; el < model.mesh().num_elements(); ++el) {
      REQUIRE(s.stresses[0].element_stress(0, el) ==
              Approx(-e * alpha * dt / (1.0 - nu)).epsilon(1e-10));
    }
  }
}

TEST_CASE("a linear temperature gradient bends a free Tet10 bar without stress",
          "[thermal][tet10][verification]") {
  // T = T_ref + a + b y makes eps = alpha (a + b y) I, which is compatible:
  // u = alpha [(a + b y) x, a y + b y^2 / 2 - b (x^2 + z^2) / 2, (a + b y) z],
  // quadratic, so the Tet10 reproduces it with zero stress.
  const IsotropicMaterial steel = thermal_steel();
  const Scalar alpha = steel.thermal_expansion();
  const Scalar a = 40.0;
  const Scalar b = 300.0;  // K/m
  const Scalar length = 1.0;
  const Scalar h = 0.1;
  const Scalar w = 0.08;
  const Mesh mesh = make_structured_tet10_mesh(spec3(4, 2, 2, length, h, w));
  FemModel model2(mesh, steel, 1.0, StressState::ThreeDimensional, IntegrationOptions());
  model2.constraints().push_back(fixed(nearest(Vector3::Zero()), true, true, true));
  model2.constraints().push_back(fixed(nearest(Vector3(0.0, h, 0.0)), true, false, true));
  model2.constraints().push_back(fixed(nearest(Vector3(0.0, 0.0, w)), true, false, false));
  // The linear field enters through a Regions spec with one node per region.
  LoadCaseSpec lc2;
  lc2.name = "gradient";
  lc2.temperature.source = TemperatureSpec::Source::Regions;
  lc2.temperature.uniform = steel.reference_temperature();
  for (Index n = 0; n < mesh.num_nodes(); ++n) {
    RegionValue rv;
    Selector sel;
    sel.kind = SelectorKind::NodeIds;
    sel.ids.assign(1, n);
    rv.region.members.push_back(sel);
    rv.value = steel.reference_temperature() + a + b * mesh.node(n).y();
    lc2.temperature.regions.push_back(rv);
  }
  model2.load_case_specs().push_back(lc2);
  const Solved s = solve(model2);
  const Vector& u = s.solutions[0].displacement;
  Scalar max_error = 0.0;
  for (Index n = 0; n < model2.mesh().num_nodes(); ++n) {
    const Vector3 x = model2.mesh().node(n);
    const Vector3 exact(alpha * (a + b * x.y()) * x.x(),
                        alpha * (a * x.y() + 0.5 * b * x.y() * x.y() -
                                 0.5 * b * (x.x() * x.x() + x.z() * x.z())),
                        alpha * (a + b * x.y()) * x.z());
    for (int k = 0; k < 3; ++k) max_error = std::max(max_error, std::abs(u(3 * n + k) - exact(k)));
  }
  REQUIRE(max_error < 1e-12);
  REQUIRE(s.stresses[0].element_stress.cwiseAbs().maxCoeff() <
          1e-6 * steel.youngs_modulus() * alpha * b * h);
}

TEST_CASE("the thermal prestress of a restrained bar enters the buckling analysis",
          "[thermal][buckling]") {
  // A bar with both end edges held axially and heated: with nu = 0 the
  // prebuckling state is exactly the uniform sigma_xx = -E alpha dT. The same
  // stress is produced mechanically by shortening the bar by alpha dT L, so
  // the two buckling analyses must agree to round-off, and both near the
  // clamped-clamped Euler load 4 pi^2 E I / L^2 of the coarse plane model.
  IsotropicMaterial m(200.0e9, 0.0, 7800.0, "steel");
  m.set_thermal(1.0e-5, 0.0, 50.0);
  const Scalar length = 1.0;
  const Scalar h = 0.02;
  const Scalar dt = 10.0;
  const auto build = [&](bool thermal) {
    FemModel model(make_structured_quad_mesh(spec2(200, 8, length, h)), m, 0.01,
                   StressState::PlaneStress, IntegrationOptions());
    // Ends held axially along their whole edge (no end rotation: clamped in
    // bending) and vertically at mid-height only, so the heated bar expands
    // freely in y and both prestresses are exactly the uniform sigma_xx.
    model.constraints().push_back(fixed(box(-kInf, 0.0), true, false));
    DisplacementConstraint right = fixed(box(length, kInf), true, false);
    if (!thermal) right.value_x = -m.thermal_expansion() * dt * length;
    model.constraints().push_back(right);
    model.constraints().push_back(fixed(nearest(Vector3(0.0, 0.5 * h, 0.0)), false, true));
    model.constraints().push_back(fixed(nearest(Vector3(length, 0.5 * h, 0.0)), false, true));
    LoadCaseSpec lc = thermal ? uniform_temperature(dt) : LoadCaseSpec();
    if (!thermal) {
      lc.name = "shortened";
      lc.prescribed_displacement_only = true;
    }
    model.load_case_specs().push_back(lc);
    model.finalize();
    return model;
  };
  FemModel heated = build(true);
  FemModel shortened = build(false);
  BucklingOptions options;
  options.num_modes = 1;
  Assembler ah(heated);
  Assembler as(shortened);
  const BucklingResult rh = analyse_buckling(heated, ah, 0, options);
  const BucklingResult rs = analyse_buckling(shortened, as, 0, options);
  REQUIRE(rh.load_factors.size() == 1);
  REQUIRE(rs.load_factors.size() == 1);
  REQUIRE(rh.load_factors(0) == Approx(rs.load_factors(0)).epsilon(1e-9));
  // Clamped-clamped Euler: N_cr = 4 pi^2 E I / L^2 against N = E A alpha dT.
  const Scalar pi = 3.14159265358979323846;
  const Scalar euler = 4.0 * pi * pi * (h * h / 12.0) / (length * length) /
                       (m.thermal_expansion() * dt);
  REQUIRE(rh.load_factors(0) == Approx(euler).epsilon(0.05));
}

// --- conduction -------------------------------------------------------------

TEST_CASE("steady conduction reproduces the exact one-dimensional profiles",
          "[thermal][conduction][verification]") {
  const Scalar k = 45.0;
  const Scalar length = 0.5;
  IsotropicMaterial m(200.0e9, 0.3, 7800.0, "steel");
  m.set_thermal(1.2e-5, 20.0, k);
  const auto left = box(-kInf, 0.0);
  const auto right = box(length, kInf);

  SECTION("heat source between two fixed temperatures (Tet10: exact everywhere)") {
    const Scalar q = 2.0e6;  // W/m^3
    FemModel model(make_structured_tet10_mesh(spec3(3, 1, 1, length, 0.1, 0.1)), m, 1.0,
                   StressState::ThreeDimensional, IntegrationOptions());
    ConductionSpec spec;
    spec.prescribed.push_back({left, false, 100.0});
    spec.prescribed.push_back({right, false, 40.0});
    RegionValue source;
    source.whole_model = true;
    source.value = q;
    spec.sources.push_back(source);
    const ConductionResult r = solve_conduction(model, spec, LinearSolverOptions());
    for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
      const Scalar x = model.mesh().node(n).x();
      const Scalar exact = 100.0 + (40.0 - 100.0) * x / length + q * x * (length - x) / (2.0 * k);
      REQUIRE(r.temperature(n) == Approx(exact).epsilon(1e-12));
    }
    // All the generated heat leaves through the two fixed faces.
    REQUIRE(r.summary.applied_heat == Approx(q * length * 0.1 * 0.1).epsilon(1e-12));
    REQUIRE(r.summary.relative_balance_error < 1e-12);
  }
  SECTION("convection at the far end (Q4: exact linear profile)") {
    const Scalar h = 120.0;
    const Scalar ambient = 25.0;
    FemModel model(make_structured_quad_mesh(spec2(6, 2, length, 0.1)), m, 0.02,
                   StressState::PlaneStress, IntegrationOptions());
    ConductionSpec spec;
    spec.prescribed.push_back({left, false, 300.0});
    ConvectionSpec c;
    c.region = right;
    c.film_coefficient = h;
    c.ambient = ambient;
    spec.convection.push_back(c);
    const ConductionResult r = solve_conduction(model, spec, LinearSolverOptions());
    for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
      const Scalar x = model.mesh().node(n).x();
      const Scalar exact = 300.0 + (ambient - 300.0) * h * x / (k + h * length);
      REQUIRE(r.temperature(n) == Approx(exact).epsilon(1e-12));
    }
    REQUIRE(r.summary.relative_balance_error < 1e-12);
  }
  SECTION("a flux into the far end (Hex8: exact linear profile)") {
    const Scalar flux = 5.0e4;  // W/m^2
    FemModel model(make_perturbed_hex_mesh(spec3(4, 2, 2, length, 0.1, 0.1), 0.0, 1u), m, 1.0,
                   StressState::ThreeDimensional, IntegrationOptions());
    ConductionSpec spec;
    spec.prescribed.push_back({left, false, 20.0});
    spec.fluxes.push_back({right, false, flux});
    const ConductionResult r = solve_conduction(model, spec, LinearSolverOptions());
    for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
      const Scalar x = model.mesh().node(n).x();
      REQUIRE(r.temperature(n) == Approx(20.0 + flux * x / k).epsilon(1e-12));
    }
    REQUIRE(r.summary.prescribed_heat == Approx(flux * 0.1 * 0.1).epsilon(1e-12));
  }
  SECTION("fluxes alone are refused") {
    FemModel model(make_structured_quad_mesh(spec2(2, 2, length, 0.1)), m, 0.02,
                   StressState::PlaneStress, IntegrationOptions());
    ConductionSpec spec;
    spec.fluxes.push_back({right, false, 1.0e3});
    REQUIRE_THROWS_AS(solve_conduction(model, spec, LinearSolverOptions()), ConfigError);
  }
}

TEST_CASE("a conducted temperature field drives the thermal strain of a load case",
          "[thermal][conduction]") {
  // A Tet10 bar with T = 20 at x = 0 and 120 at x = L: the conducted field is
  // linear, T - T_ref = b x, and the free expansion is the quadratic
  // u = alpha b [x^2/2 - (y^2 + z^2)/2, x y, x z], which the Tet10 reproduces
  // with zero stress on statically determinate supports.
  IsotropicMaterial m(200.0e9, 0.3, 7800.0, "steel");
  m.set_thermal(1.2e-5, 20.0, 45.0);
  const Scalar length = 1.0;
  const Scalar h = 0.1;
  FemModel model(make_structured_tet10_mesh(spec3(4, 1, 1, length, h, h)), m, 1.0,
                 StressState::ThreeDimensional, IntegrationOptions());
  model.constraints().push_back(fixed(nearest(Vector3::Zero()), true, true, true));
  model.constraints().push_back(fixed(nearest(Vector3(length, 0.0, 0.0)), false, true, true));
  model.constraints().push_back(fixed(nearest(Vector3(0.0, h, 0.0)), false, false, true));
  LoadCaseSpec lc;
  lc.name = "conducted";
  lc.temperature.source = TemperatureSpec::Source::Conduction;
  lc.temperature.conduction.prescribed.push_back({box(-kInf, 0.0), false, 20.0});
  lc.temperature.conduction.prescribed.push_back({box(length, kInf), false, 120.0});
  model.load_case_specs().push_back(lc);
  const Solved s = solve(model);
  REQUIRE(model.load_case_data(0).conduction_solved);
  const Scalar b = 100.0 / length;
  const Scalar alpha = m.thermal_expansion();
  const Vector& u = s.solutions[0].displacement;
  for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
    const Vector3 x = model.mesh().node(n);
    REQUIRE(std::abs(u(3 * n) - alpha * b * (0.5 * x.x() * x.x() -
                                              0.5 * (x.y() * x.y() + x.z() * x.z()))) < 1e-12);
    REQUIRE(std::abs(u(3 * n + 1) - alpha * b * x.x() * x.y()) < 1e-12);
    REQUIRE(std::abs(u(3 * n + 2) - alpha * b * x.x() * x.z()) < 1e-12);
  }
  REQUIRE(s.stresses[0].element_stress.cwiseAbs().maxCoeff() <
          1e-6 * m.youngs_modulus() * alpha * 100.0);
}

// --- materials --------------------------------------------------------------

TEST_CASE("two materials in series give the exact tip displacement", "[materials]") {
  const Scalar length = 2.0;
  const Scalar area = 0.1 * 0.01;
  const Scalar force = 1.0e4;
  const IsotropicMaterial steel(200.0e9, 0.0, 7800.0, "steel");
  const IsotropicMaterial aluminium(70.0e9, 0.0, 2700.0, "aluminium");
  FemModel model(make_structured_quad_mesh(spec2(8, 2, length, 0.1)), steel, 0.01,
                 StressState::PlaneStress, IntegrationOptions());
  std::vector<Index> right_half;
  for (Index e = 0; e < model.mesh().num_elements(); ++e) {
    if (model.mesh().element_centroid(e).x() > 0.5 * length) right_half.push_back(e);
  }
  model.assign_material(aluminium, right_half);
  REQUIRE(model.num_materials() == 2);
  REQUIRE_FALSE(model.single_material());
  model.constraints().push_back(fixed(box(-kInf, 0.0), true, false));
  model.constraints().push_back(fixed(nearest(Vector3::Zero()), false, true));
  LoadCaseSpec lc;
  lc.name = "pull";
  TractionLoadSpec p;
  p.region = box(length, kInf);
  p.traction = Vector3(force / area, 0.0, 0.0);
  lc.tractions.push_back(p);
  model.load_case_specs().push_back(lc);
  const Solved s = solve(model);
  const Scalar expected = force * 0.5 * length / area * (1.0 / 200.0e9 + 1.0 / 70.0e9);
  for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
    if (model.mesh().node(n).x() == length) {
      REQUIRE(s.solutions[0].displacement(2 * n) == Approx(expected).epsilon(1e-12));
    }
  }
  // Both halves carry the same stress, force / area.
  REQUIRE(s.stresses[0].element_stress.row(0).minCoeff() == Approx(force / area).epsilon(1e-10));
  REQUIRE(s.stresses[0].element_stress.row(0).maxCoeff() == Approx(force / area).epsilon(1e-10));
  // The mass of the model uses each element's density.
  Assembler assembler(model);
  REQUIRE(assembler.total_mass() == Approx(0.5 * length * area * (7800.0 + 2700.0)).epsilon(1e-12));
}
