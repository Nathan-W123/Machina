/// \file test_contact.cpp
/// \brief Unilateral contact: rigid obstacles and mortar pairs, frictionless
///        and with Coulomb friction.
///
/// The exact checks are homogeneous states, which the linear elements
/// reproduce on distorted meshes: a block pressed onto a rigid plane (with
/// and without an initial gap), two blocks with non-matching meshes pressed
/// together (the contact patch test, which the dual mortar method passes and
/// node-to-segment contact does not), and a block dragged along a plane with
/// friction, whose every node slips and whose tangential force is then mu
/// times the normal one. The rest are the invariants of any solution -
/// complementarity, the force balance - and the refusals.
#include "TestSupport.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/fem/Contact.hpp"
#include "sparlab/fem/NonlinearStatic.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <limits>

using namespace sparlab;
using namespace sparlab::testing;
using Catch::Approx;

namespace {

constexpr Scalar kInf = std::numeric_limits<Scalar>::infinity();

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
  return g;
}

Mesh block(ElementType type, Index nx, Index ny, Index nz, Scalar lx, Scalar ly, Scalar lz,
           Scalar y0 = 0.0, unsigned int seed = 12345u) {
  StructuredMeshSpec spec;
  spec.nx = nx;
  spec.ny = ny;
  spec.nz = nz;
  spec.lx = lx;
  spec.ly = ly;
  spec.lz = lz;
  spec.y0 = y0;
  switch (type) {
    case ElementType::Quad4: return make_perturbed_quad_mesh(spec, 0.25, seed);
    case ElementType::Tri3: return make_perturbed_tri_mesh(spec, 0.25, seed);
    case ElementType::Hex8: return make_perturbed_hex_mesh(spec, 0.2, seed);
    case ElementType::Tet4: return make_perturbed_tet_mesh(spec, 0.2, seed);
    default: break;
  }
  return make_structured_tet10_mesh(spec);
}

/// Two meshes as one (two bodies).
Mesh merge(const Mesh& a, const Mesh& b) {
  Matrix coords(a.dim(), a.num_nodes() + b.num_nodes());
  coords << a.coordinates(), b.coordinates();
  std::vector<Index> connectivity = a.connectivity();
  for (Index n : b.connectivity()) connectivity.push_back(n + a.num_nodes());
  return Mesh(coords, connectivity, a.element_type());
}

StressState state_of(const Mesh& m) {
  return m.dim() == 2 ? StressState::PlaneStress : StressState::ThreeDimensional;
}

DisplacementConstraint fix(const SelectorGroup& region, int comp, Scalar value = 0.0) {
  DisplacementConstraint c;
  c.region = region;
  c.set(comp, true, value);
  return c;
}

NonlinearOptions contact_options(const ContactPairSpec& pair) {
  NonlinearOptions o;
  o.kinematics = Kinematics::SmallStrain;
  o.steps = 1;
  o.residual_tolerance = 1.0e-12;
  o.displacement_tolerance = 1.0e-12;
  o.contact.enabled = true;
  o.contact.pairs.push_back(pair);
  return o;
}

}  // namespace

TEST_CASE("a block pressed onto a rigid plane takes a uniform pressure, exact on distorted meshes",
          "[contact]") {
  // The top of the block is pushed down by delta; its bottom meets the plane
  // after the initial gap g0. Frictionless and free to expand sideways, the
  // block is in uniaxial stress sigma_yy = -E (delta - g0) / H whatever the
  // mesh, so every bottom node is in contact with pressure E (delta - g0) / H
  // and zero gap, and the displacement is the exact linear field.
  const IsotropicMaterial m = default_material();
  const Scalar e = m.youngs_modulus();
  const Scalar nu = m.poisson_ratio();
  const Scalar lx = 0.4;
  const Scalar h = 0.2;
  const Scalar lz = 0.3;
  const Scalar delta = 1.0e-4;
  for (const ElementType type :
       {ElementType::Quad4, ElementType::Tri3, ElementType::Hex8, ElementType::Tet4}) {
    for (const Scalar g0 : {0.0, 0.25 * delta}) {
      INFO(to_string(type) << ", initial gap " << g0);
      Mesh mesh = block(type, 4, 3, 3, lx, h, lz);
      const int dim = mesh.dim();
      const StressState state = state_of(mesh);
      const Scalar thickness = dim == 2 ? 0.01 : 1.0;
      FemModel model(std::move(mesh), m, thickness, state, IntegrationOptions());
      model.constraints().push_back(fix(box(-kInf, kInf, h, kInf), 1, -delta));
      // The sliding modes the frictionless contact leaves free: x (and z) at
      // one corner, and in 3-D the turn about y at another.
      model.constraints().push_back(fix(box(-kInf, 0.0, -kInf, 0.0, -kInf, 0.0), 0));
      if (dim == 3) {
        model.constraints().push_back(fix(box(-kInf, 0.0, -kInf, 0.0, -kInf, 0.0), 2));
        model.constraints().push_back(fix(box(lx, kInf, -kInf, 0.0, -kInf, 0.0), 2));
      }
      LoadCaseSpec lc;
      lc.name = "press";
      lc.prescribed_displacement_only = true;
      model.load_case_specs().push_back(lc);
      model.finalize();
      Assembler assembler(model);
      ContactPairSpec pair;
      pair.name = "floor";
      pair.slave = box(-kInf, kInf, -kInf, 0.0);
      pair.obstacle.kind = RigidObstacle::Kind::Plane;
      pair.obstacle.point = Vector3(0.0, -g0, 0.0);
      pair.obstacle.direction = Vector3::UnitY();
      NonlinearStaticAnalysis analysis(model, assembler, contact_options(pair));
      const NonlinearResult r = analysis.solve(0);
      REQUIRE(r.completed);
      const Scalar strain = (delta - g0) / h;
      const Scalar p = e * strain;
      REQUIRE(!r.contact_nodes.empty());
      for (const ContactNodeResult& c : r.contact_nodes) {
        REQUIRE(c.status == ContactStatus::Slip);  // frictionless nodes in contact
        REQUIRE(std::abs(c.gap) <= 1.0e-12 * delta);
        REQUIRE(c.pressure == Approx(p).epsilon(1.0e-9));
      }
      // The exact displacement field.
      Scalar worst = 0.0;
      for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
        const Vector3 x = model.mesh().node(n);
        Vector3 exact = Vector3::Zero();
        exact.y() = -g0 - strain * x.y();
        exact.x() = nu * strain * x.x();
        if (dim == 3) exact.z() = nu * strain * x.z();
        for (int k = 0; k < dim; ++k) {
          worst = std::max(worst, std::abs(r.displacement(n * dim + k) - exact(k)));
        }
      }
      REQUIRE(worst <= 1.0e-10 * delta);
      // The plane's force is the pressure over the bottom area, and it
      // balances the reactions of the pushed top.
      const Scalar area = lx * (dim == 2 ? thickness : lz);
      REQUIRE(r.contact_pairs.size() == 1);
      REQUIRE(r.contact_pairs[0].force.y() == Approx(p * area).epsilon(1.0e-9));
      REQUIRE(r.contact_pairs[0].area == Approx(area).epsilon(1.0e-12));
      REQUIRE(r.equilibrium.relative_force_error <= 1.0e-10);
    }
  }
}

TEST_CASE("two blocks with non-matching meshes pass the contact patch test (dual mortar)",
          "[contact]") {
  // An upper block on a lower one, their meshes not matching at the
  // interface, a small initial gap g0 between them; the upper block's top
  // is pushed down by delta, the lower block's bottom is held vertically.
  // With the same material the two blocks are in the same uniaxial stress
  // and expand sideways alike, so the exact state is homogeneous in each:
  // sigma_yy = -E (delta - g0) / (H1 + H2) everywhere, the contact pressure
  // uniform. The dual mortar method reproduces it on distorted meshes to
  // round-off - the contact patch test, which node-to-segment contact fails
  // on non-matching meshes.
  const IsotropicMaterial m = default_material();
  const Scalar e = m.youngs_modulus();
  const Scalar nu = m.poisson_ratio();
  const Scalar lx = 0.4;
  const Scalar lz = 0.3;
  const Scalar h1 = 0.2;
  const Scalar h2 = 0.15;
  const Scalar g0 = 2.0e-5;
  const Scalar delta = 1.0e-4;
  for (const ElementType type :
       {ElementType::Quad4, ElementType::Tri3, ElementType::Hex8, ElementType::Tet4}) {
    INFO(to_string(type));
    const Mesh lower = block(type, 4, 3, 3, lx, h1, lz, 0.0, 11u);
    const Mesh upper = block(type, 5, 2, 4, lx, h2, lz, h1 + g0, 17u);
    Mesh mesh = merge(lower, upper);
    const int dim = mesh.dim();
    const Scalar thickness = dim == 2 ? 0.01 : 1.0;
    FemModel model(std::move(mesh), m, thickness, state_of(lower), IntegrationOptions());
    const Scalar top = h1 + g0 + h2;
    model.constraints().push_back(fix(box(-kInf, kInf, top, kInf), 1, -delta));
    model.constraints().push_back(fix(box(-kInf, kInf, -kInf, 0.0), 1));
    // x (and z) at the corner of each block on x = 0 (z = 0), and in 3-D the
    // turn about y of each block.
    for (const Scalar y : {0.0, h1 + g0}) {
      model.constraints().push_back(fix(box(-kInf, 0.0, y, y, -kInf, 0.0), 0));
      if (dim == 3) {
        model.constraints().push_back(fix(box(-kInf, 0.0, y, y, -kInf, 0.0), 2));
        model.constraints().push_back(fix(box(lx, kInf, y, y, -kInf, 0.0), 2));
      }
    }
    LoadCaseSpec lc;
    lc.name = "press";
    lc.prescribed_displacement_only = true;
    model.load_case_specs().push_back(lc);
    model.finalize();
    Assembler assembler(model);
    ContactPairSpec pair;
    pair.name = "interface";
    pair.rigid = false;
    pair.slave = box(-kInf, kInf, h1 + g0, h1 + g0);   // bottom of the upper block
    pair.master = box(-kInf, kInf, h1, h1);             // top of the lower block
    NonlinearStaticAnalysis analysis(model, assembler, contact_options(pair));
    const NonlinearResult r = analysis.solve(0);
    REQUIRE(r.completed);
    const Scalar strain = (delta - g0) / (h1 + h2);
    const Scalar p = e * strain;
    REQUIRE(!r.contact_nodes.empty());
    Scalar worst_p = 0.0;
    for (const ContactNodeResult& c : r.contact_nodes) {
      REQUIRE(c.status == ContactStatus::Slip);
      REQUIRE(std::abs(c.gap) <= 1.0e-12 * delta);
      worst_p = std::max(worst_p, std::abs(c.pressure - p) / p);
    }
    REQUIRE(worst_p <= 1.0e-9);
    // Homogeneous in each block: the lower one from its held bottom, the
    // upper one from the lower one's top plus the closed gap.
    Scalar worst = 0.0;
    for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
      const Vector3 x = model.mesh().node(n);
      const bool in_upper = n >= lower.num_nodes();
      Vector3 exact = Vector3::Zero();
      exact.y() = in_upper ? -g0 - strain * (x.y() - g0) : -strain * x.y();
      exact.x() = nu * strain * x.x();
      if (dim == 3) exact.z() = nu * strain * x.z();
      for (int k = 0; k < dim; ++k) {
        worst = std::max(worst, std::abs(r.displacement(n * dim + k) - exact(k)));
      }
    }
    REQUIRE(worst <= 1.0e-10 * delta);
    REQUIRE(r.contact_pairs.size() == 1);
    const Scalar area = lx * (dim == 2 ? thickness : lz);
    // The lower block pushes the upper (slave) one up.
    REQUIRE(r.contact_pairs[0].force.y() == Approx(p * area).epsilon(1.0e-9));
    REQUIRE(r.equilibrium.relative_force_error <= 1.0e-10);
  }
}

TEST_CASE("a block dragged along a rigid plane slips at every node with the Coulomb traction",
          "[contact][friction]") {
  // The top is pushed down by delta and sideways by s, in proportion. With
  // s well beyond what the shear stiffness transmits without slip
  // (s > mu 2 (1 + nu) delta for a block in simple shear) every bottom node
  // slips, its friction traction is mu p against the slip, and so the
  // tangential force on the block is mu times the normal one - whatever the
  // distribution of the pressure. A short sideways push sticks instead: the
  // tangential force stays below mu N and no sticking node carries more
  // than mu p.
  const IsotropicMaterial m = default_material();
  const Scalar mu = 0.3;
  const Scalar lx = 0.4;
  const Scalar h = 0.2;
  const Scalar lz = 0.3;
  const Scalar delta = 1.0e-4;
  for (const ElementType type : {ElementType::Quad4, ElementType::Hex8, ElementType::Tet4}) {
    for (const Scalar ratio : {5.0, 0.1}) {
      INFO(to_string(type) << ", sideways / down " << ratio);
      Mesh mesh = block(type, 4, 3, 3, lx, h, lz);
      const int dim = mesh.dim();
      const StressState state = state_of(mesh);
      const Scalar thickness = dim == 2 ? 0.01 : 1.0;
      FemModel model(std::move(mesh), m, thickness, state, IntegrationOptions());
      const SelectorGroup top = box(-kInf, kInf, h, kInf);
      model.constraints().push_back(fix(top, 1, -delta));
      model.constraints().push_back(fix(top, 0, ratio * delta));
      if (dim == 3) model.constraints().push_back(fix(top, 2));
      LoadCaseSpec lc;
      lc.name = "drag";
      lc.prescribed_displacement_only = true;
      model.load_case_specs().push_back(lc);
      model.finalize();
      Assembler assembler(model);
      ContactPairSpec pair;
      pair.name = "floor";
      pair.slave = box(-kInf, kInf, -kInf, 0.0);
      pair.obstacle.kind = RigidObstacle::Kind::Plane;
      pair.obstacle.direction = Vector3::UnitY();
      pair.friction = mu;
      NonlinearOptions options = contact_options(pair);
      options.steps = 4;
      NonlinearMonitor fx;
      fx.name = "fx";
      fx.region = top;
      fx.component = 0;
      fx.quantity = NonlinearMonitor::Quantity::Reaction;
      NonlinearMonitor fy = fx;
      fy.name = "fy";
      fy.component = 1;
      options.monitors = {fx, fy};
      NonlinearStaticAnalysis analysis(model, assembler, options);
      const NonlinearResult r = analysis.solve(0);
      REQUIRE(r.completed);
      const Scalar rx = r.steps.back().monitors[0];
      const Scalar ry = r.steps.back().monitors[1];
      // The reactions are the forces the supports exert on the block: the
      // top holds it down (and drags it along +x).
      REQUIRE(ry < 0.0);
      const ContactPairResult& pr = r.contact_pairs.at(0);
      REQUIRE(pr.active == pr.nodes);
      if (ratio > 1.0) {
        REQUIRE(pr.slipping == pr.nodes);
        for (const ContactNodeResult& c : r.contact_nodes) {
          REQUIRE(c.status == ContactStatus::Slip);
          // The traction is mu p against the slip (mostly +x; in 3-D the
          // sideways expansion adds some z).
          REQUIRE(c.slip.x() > 0.0);
          REQUIRE(c.traction.norm() == Approx(mu * c.pressure).epsilon(1.0e-9));
          REQUIRE(c.traction.dot(c.slip.normalized()) ==
                  Approx(-c.traction.norm()).epsilon(1.0e-9));
        }
        if (dim == 2) {
          // Every traction along -x: the tangential force is mu N exactly.
          REQUIRE(std::abs(rx) == Approx(mu * std::abs(ry)).epsilon(1.0e-9));
          REQUIRE(pr.force.x() == Approx(-mu * pr.force.y()).epsilon(1.0e-9));
        } else {
          REQUIRE(std::abs(rx) <= mu * std::abs(ry));
          REQUIRE(std::abs(rx) >= 0.99 * mu * std::abs(ry));
        }
      } else {
        REQUIRE(pr.sticking > 0);
        REQUIRE(std::abs(rx) < mu * std::abs(ry));
        for (const ContactNodeResult& c : r.contact_nodes) {
          if (c.status != ContactStatus::Stick) continue;
          REQUIRE(c.traction.norm() <= mu * c.pressure * (1.0 + 1.0e-9));
          REQUIRE(c.slip.norm() <= 1.0e-12 * delta);
        }
      }
      // Complementarity and the force balance with the plane as a support.
      for (const ContactNodeResult& c : r.contact_nodes) {
        REQUIRE(c.pressure >= 0.0);
        REQUIRE(std::abs(c.gap) <= 1.0e-12 * delta);
      }
      REQUIRE(r.equilibrium.relative_force_error <= 1.0e-10);
    }
  }
}

TEST_CASE("rigid obstacles measure the signed distance and its normal", "[contact]") {
  RigidObstacle plane;
  plane.kind = RigidObstacle::Kind::Plane;
  plane.point = Vector3(0.0, 1.0, 0.0);
  plane.direction = Vector3(0.0, 2.0, 0.0);  // normalised on use
  Vector3 n;
  REQUIRE(plane.gap(Vector3(3.0, 1.5, 0.0), n, 3) == Approx(0.5));
  REQUIRE((n - Vector3::UnitY()).norm() < 1.0e-15);

  RigidObstacle cyl;
  cyl.kind = RigidObstacle::Kind::Cylinder;
  cyl.point = Vector3(1.0, 0.0, 0.0);
  cyl.direction = Vector3::UnitZ();
  cyl.radius = 0.5;
  REQUIRE(cyl.gap(Vector3(1.0, 2.0, 7.0), n, 3) == Approx(1.5));  // along the axis irrelevant
  REQUIRE((n - Vector3::UnitY()).norm() < 1.0e-15);
  cyl.inside = true;
  REQUIRE(cyl.gap(Vector3(1.3, 0.0, 0.0), n, 2) == Approx(0.2));
  REQUIRE((n + Vector3::UnitX()).norm() < 1.0e-15);

  RigidObstacle sphere;
  sphere.kind = RigidObstacle::Kind::Sphere;
  sphere.radius = 1.0;
  REQUIRE(sphere.gap(Vector3(0.0, 0.0, 3.0), n, 3) == Approx(2.0));
  REQUIRE((n - Vector3::UnitZ()).norm() < 1.0e-15);
  REQUIRE_THROWS_AS(sphere.gap(Vector3::Zero(), n, 3), ModelError);
  REQUIRE(parse_obstacle_kind("cylinder") == RigidObstacle::Kind::Cylinder);
  REQUIRE_THROWS_AS(parse_obstacle_kind("cone"), ConfigError);
}

TEST_CASE("contact refuses what it cannot model and says why", "[contact]") {
  const IsotropicMaterial m = default_material();
  const auto model_of = [&](ElementType type) {
    Mesh mesh = block(type, 2, 2, 2, 0.4, 0.2, 0.3);
    const StressState state = state_of(mesh);
    FemModel model(std::move(mesh), m, mesh.dim() == 2 ? 0.01 : 1.0, state,
                   IntegrationOptions());
    return model;
  };
  ContactPairSpec pair;
  pair.name = "floor";
  pair.slave = box(-kInf, kInf, -kInf, 0.0);
  pair.obstacle.direction = Vector3::UnitY();
  {
    // Finite kinematics and the arc-length method.
    FemModel model = model_of(ElementType::Quad4);
    LoadCaseSpec lc;
    lc.name = "x";
    model.load_case_specs().push_back(lc);
    model.constraints().push_back(fix(box(-kInf, kInf, 0.2, kInf), 1, -1.0e-4));
    model.finalize();
    Assembler assembler(model);
    NonlinearOptions o = contact_options(pair);
    o.kinematics = Kinematics::Finite;
    REQUIRE_THROWS_AS(NonlinearStaticAnalysis(model, assembler, o), ConfigError);
    o = contact_options(pair);
    o.method = NonlinearOptions::Method::ArcLength;
    REQUIRE_THROWS_AS(NonlinearStaticAnalysis(model, assembler, o), ConfigError);
    // A slave surface that selects nothing, and a node on both surfaces.
    ContactOptions c;
    c.enabled = true;
    ContactPairSpec none = pair;
    none.slave = box(10.0, kInf);
    c.pairs = {none};
    REQUIRE_THROWS_AS(ContactProblem(model, c), ConfigError);
    ContactPairSpec both = pair;
    both.rigid = false;
    both.master = box(-kInf, kInf, -kInf, 0.0);
    c.pairs = {both};
    REQUIRE_THROWS_AS(ContactProblem(model, c), ConfigError);
    ContactPairSpec negative = pair;
    negative.friction = -0.1;
    c.pairs = {negative};
    REQUIRE_THROWS_AS(ContactProblem(model, c), ConfigError);
  }
  {
    // Quadratic tetrahedra.
    FemModel model = model_of(ElementType::Tet10);
    LoadCaseSpec lc;
    lc.name = "x";
    model.load_case_specs().push_back(lc);
    model.finalize();
    ContactOptions c;
    c.enabled = true;
    c.pairs = {pair};
    REQUIRE_THROWS_AS(ContactProblem(model, c), ConfigError);
  }
  {
    // A block held only by frictionless contact can slide: the run stops
    // and says why rather than returning a meaningless state.
    FemModel model = model_of(ElementType::Quad4);
    model.constraints().push_back(fix(box(-kInf, kInf, 0.2, kInf), 1, -1.0e-4));
    LoadCaseSpec lc;
    lc.name = "press";
    lc.prescribed_displacement_only = true;
    model.load_case_specs().push_back(lc);
    model.finalize();
    Assembler assembler(model);
    NonlinearOptions o = contact_options(pair);
    o.max_cuts = 2;
    NonlinearStaticAnalysis analysis(model, assembler, o);
    const NonlinearResult r = analysis.solve(0);
    REQUIRE_FALSE(r.completed);
    REQUIRE(r.termination.find("rigid-body motion") != std::string::npos);
  }
}

TEST_CASE("a block pressed past yield onto a rigid plane follows the exact uniaxial curve",
          "[contact][plasticity]") {
  // J2 with linear isotropic hardening in uniaxial stress: sigma = E eps
  // up to yield, then sigma_y + E H (eps - sigma_y / E) / (E + H). The block
  // is pressed down by delta along a path that yields; its bottom is in
  // contact with the plane throughout, and the contact pressure is the
  // uniaxial stress at every step, on distorted meshes. The plastic history
  // is committed only with the contact status settled.
  IsotropicMaterial m(200.0e9, 0.3, 7800.0, "steel");
  PlasticityParameters pp;
  pp.yield_stress = 250.0e6;
  pp.hardening_modulus = 10.0e9;
  m.set_plasticity(pp);
  const Scalar e = m.youngs_modulus();
  const Scalar hh = pp.hardening_modulus;
  const Scalar sy = pp.yield_stress;
  const Scalar lx = 0.4;
  const Scalar h = 0.2;
  const Scalar lz = 0.3;
  const Scalar delta = 3.0 * sy / e * h;  // three times the yield strain
  const auto exact = [&](Scalar strain) {
    return strain * e <= sy ? e * strain : sy + e * hh * (strain - sy / e) / (e + hh);
  };
  for (const ElementType type : {ElementType::Quad4, ElementType::Hex8}) {
    INFO(to_string(type));
    Mesh mesh = block(type, 4, 3, 3, lx, h, lz);
    const int dim = mesh.dim();
    const StressState state = state_of(mesh);
    const Scalar thickness = dim == 2 ? 0.01 : 1.0;
    FemModel model(std::move(mesh), m, thickness, state, IntegrationOptions());
    model.constraints().push_back(fix(box(-kInf, kInf, h, kInf), 1, -delta));
    model.constraints().push_back(fix(box(-kInf, 0.0, -kInf, 0.0, -kInf, 0.0), 0));
    if (dim == 3) {
      model.constraints().push_back(fix(box(-kInf, 0.0, -kInf, 0.0, -kInf, 0.0), 2));
      model.constraints().push_back(fix(box(lx, kInf, -kInf, 0.0, -kInf, 0.0), 2));
    }
    LoadCaseSpec lc;
    lc.name = "press";
    lc.prescribed_displacement_only = true;
    model.load_case_specs().push_back(lc);
    model.finalize();
    Assembler assembler(model);
    ContactPairSpec pair;
    pair.name = "floor";
    pair.slave = box(-kInf, kInf, -kInf, 0.0);
    pair.obstacle.direction = Vector3::UnitY();
    NonlinearOptions options = contact_options(pair);
    options.steps = 6;
    NonlinearMonitor fy;
    fy.name = "fy";
    fy.region = box(-kInf, kInf, h, kInf);
    fy.component = 1;
    fy.quantity = NonlinearMonitor::Quantity::Reaction;
    options.monitors = {fy};
    NonlinearStaticAnalysis analysis(model, assembler, options);
    const NonlinearResult r = analysis.solve(0);
    REQUIRE(r.completed);
    REQUIRE(r.plastic);
    REQUIRE(r.max_plastic_strain > 0.0);
    const Scalar area = lx * (dim == 2 ? thickness : lz);
    for (const NonlinearStep& s : r.steps) {
      const Scalar p = exact(s.load_factor * delta / h);
      REQUIRE(-s.monitors[0] == Approx(p * area).epsilon(1.0e-9));
      REQUIRE(s.contact_force.at(0).y() == Approx(p * area).epsilon(1.0e-9));
    }
    const Scalar p = exact(delta / h);
    for (const ContactNodeResult& c : r.contact_nodes) {
      REQUIRE(c.status != ContactStatus::Open);
      REQUIRE(c.pressure == Approx(p).epsilon(1.0e-9));
    }
  }
}

TEST_CASE("friction on a mortar pair: stacked blocks stick, a dragged block slips everywhere",
          "[contact][friction]") {
  // Stacked blocks of one material in uniaxial stress expand sideways
  // alike: nothing slides, so with friction every node sticks, carries no
  // traction, and the state is the frictionless one. An upper block dragged
  // sideways over a lower one held at its bottom slips at every node, and
  // in 2-D its tangential force is then mu times the normal one exactly.
  const IsotropicMaterial m = default_material();
  const Scalar e = m.youngs_modulus();
  const Scalar mu = 0.3;
  const Scalar lx = 0.4;
  const Scalar lz = 0.3;
  const Scalar h1 = 0.2;
  const Scalar h2 = 0.15;
  const Scalar delta = 1.0e-4;
  for (const ElementType type : {ElementType::Quad4, ElementType::Hex8}) {
    for (const bool drag : {false, true}) {
      INFO(to_string(type) << (drag ? ", dragged" : ", stacked"));
      const Mesh lower = block(type, 4, 3, 3, lx, h1, lz, 0.0, 11u);
      const Mesh upper = block(type, 5, 2, 4, lx, h2, lz, h1, 17u);
      Mesh mesh = merge(lower, upper);
      const int dim = mesh.dim();
      const Scalar thickness = dim == 2 ? 0.01 : 1.0;
      FemModel model(std::move(mesh), m, thickness, state_of(lower), IntegrationOptions());
      const SelectorGroup top = box(-kInf, kInf, h1 + h2, kInf);
      model.constraints().push_back(fix(top, 1, -delta));
      if (drag) {
        model.constraints().push_back(fix(top, 0, 5.0 * delta));
        if (dim == 3) model.constraints().push_back(fix(top, 2));
        for (int k = 0; k < dim; ++k) {
          model.constraints().push_back(fix(box(-kInf, kInf, -kInf, 0.0), k));
        }
      } else {
        model.constraints().push_back(fix(box(-kInf, kInf, -kInf, 0.0), 1));
        for (const Scalar y : {0.0, h1}) {
          model.constraints().push_back(fix(box(-kInf, 0.0, y, y, -kInf, 0.0), 0));
          if (dim == 3) {
            model.constraints().push_back(fix(box(-kInf, 0.0, y, y, -kInf, 0.0), 2));
            model.constraints().push_back(fix(box(lx, kInf, y, y, -kInf, 0.0), 2));
          }
        }
      }
      LoadCaseSpec lc;
      lc.name = "load";
      lc.prescribed_displacement_only = true;
      model.load_case_specs().push_back(lc);
      model.finalize();
      Assembler assembler(model);
      ContactPairSpec pair;
      pair.name = "interface";
      pair.rigid = false;
      pair.friction = mu;
      // The two blocks touch; each surface is selected through its own
      // block's side of the interface plane (the upper block's bottom
      // nodes are numbered after the lower block's).
      pair.slave = box(-kInf, kInf, h1, h1);
      pair.master = pair.slave;
      pair.slave.members.clear();
      pair.master.members.clear();
      SelectorGroup up_nodes;
      Selector ids;
      ids.kind = SelectorKind::NodeIds;
      for (Index n = lower.num_nodes(); n < lower.num_nodes() + upper.num_nodes(); ++n) {
        if (std::abs(model.mesh().node(n).y() - h1) < 1.0e-12) ids.ids.push_back(n);
      }
      pair.slave.members.push_back(ids);
      Selector lids;
      lids.kind = SelectorKind::NodeIds;
      for (Index n = 0; n < lower.num_nodes(); ++n) {
        if (std::abs(model.mesh().node(n).y() - h1) < 1.0e-12) lids.ids.push_back(n);
      }
      pair.master.members.push_back(lids);
      NonlinearOptions options = contact_options(pair);
      options.steps = 3;
      NonlinearMonitor fx;
      fx.name = "fx";
      fx.region = top;
      fx.component = 0;
      fx.quantity = NonlinearMonitor::Quantity::Reaction;
      NonlinearMonitor fy = fx;
      fy.name = "fy";
      fy.component = 1;
      options.monitors = {fx, fy};
      NonlinearStaticAnalysis analysis(model, assembler, options);
      const NonlinearResult r = analysis.solve(0);
      REQUIRE(r.completed);
      const ContactPairResult& pr = r.contact_pairs.at(0);
      REQUIRE(pr.active == pr.nodes);
      if (!drag) {
        const Scalar p = e * delta / (h1 + h2);
        REQUIRE(pr.sticking == pr.nodes);
        for (const ContactNodeResult& c : r.contact_nodes) {
          REQUIRE(c.pressure == Approx(p).epsilon(1.0e-9));
          REQUIRE(c.traction.norm() <= 1.0e-9 * p);
        }
      } else {
        REQUIRE(pr.slipping == pr.nodes);
        for (const ContactNodeResult& c : r.contact_nodes) {
          REQUIRE(c.slip.x() > 0.0);
          REQUIRE(c.traction.norm() == Approx(mu * c.pressure).epsilon(1.0e-9));
        }
        const Scalar rx = r.steps.back().monitors[0];
        const Scalar ry = r.steps.back().monitors[1];
        if (dim == 2) {
          REQUIRE(std::abs(rx) == Approx(mu * std::abs(ry)).epsilon(1.0e-9));
        } else {
          REQUIRE(std::abs(rx) <= mu * std::abs(ry) * (1.0 + 1.0e-12));
        }
      }
      REQUIRE(r.equilibrium.relative_force_error <= 1.0e-10);
    }
  }
}
