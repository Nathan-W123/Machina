/// \file test_forming.cpp
/// \brief Incremental forming: tool trajectories and the penalty contact of
///        rigid tools in the current configuration with Coulomb friction.
///
/// The contact kernel is checked against central finite differences of its
/// own residual (sphere, plane and cylinder; frictionless, sticking and
/// slipping), and the trajectories against their definition.
#include "TestSupport.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"
#include "sparlab/fem/RigidTool.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <random>
#include <string>

using namespace sparlab;
using namespace sparlab::testing;
using Catch::Approx;
using Catch::Matchers::ContainsSubstring;

namespace {

constexpr Scalar kInf = std::numeric_limits<Scalar>::infinity();

SelectorGroup box(Scalar xmin, Scalar xmax, Scalar ymin = -kInf, Scalar ymax = kInf,
                  Scalar zmin = -kInf, Scalar zmax = kInf, const std::string& name = "box") {
  SelectorGroup g;
  g.name = name;
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

ToolTrajectory path(std::vector<Scalar> times, std::vector<Vector3> points) {
  ToolTrajectory t;
  t.times = std::move(times);
  t.points = std::move(points);
  return t;
}

/// A finalised model with one load case that applies nothing.
FemModel finalised(Mesh mesh, const IsotropicMaterial& m, Scalar thickness = 1.0,
                   StressState state = StressState::ThreeDimensional,
                   const std::vector<DisplacementConstraint>& bcs = {}) {
  FemModel model(std::move(mesh), m, thickness, state, IntegrationOptions());
  model.constraints() = bcs;
  LoadCaseSpec lc;
  lc.name = "none";
  lc.prescribed_displacement_only = true;
  model.load_case_specs().push_back(lc);
  model.finalize();
  return model;
}

Mesh hex_block(Index nx, Index ny, Index nz, Scalar lx, Scalar ly, Scalar lz, Scalar x0 = 0.0,
               Scalar y0 = 0.0, Scalar z0 = 0.0) {
  StructuredMeshSpec spec;
  spec.nx = nx;
  spec.ny = ny;
  spec.nz = nz;
  spec.lx = lx;
  spec.ly = ly;
  spec.lz = lz;
  spec.x0 = x0;
  spec.y0 = y0;
  spec.z0 = z0;
  return make_structured_hex_mesh(spec);
}

}  // namespace

// ---------------------------------------------------------------------------
// 1. Trajectories
// ---------------------------------------------------------------------------
TEST_CASE("a tool trajectory interpolates linearly, holds its ends and reads a CSV file",
          "[forming]") {
  const ToolTrajectory t = path({0.0, 1.0, 3.0}, {Vector3(0, 0, 1), Vector3(2, 0, 1),
                                                  Vector3(2, 4, 0)});
  t.validate();
  CHECK((t.position(0.5) - Vector3(1, 0, 1)).norm() < 1e-15);
  CHECK((t.position(2.0) - Vector3(2, 2, 0.5)).norm() < 1e-15);
  CHECK((t.position(-1.0) - Vector3(0, 0, 1)).norm() == 0.0);  // clamped before
  CHECK((t.position(7.0) - Vector3(2, 4, 0)).norm() == 0.0);   // and after
  CHECK((t.position(1.0) - Vector3(2, 0, 1)).norm() < 1e-15);  // a knot exactly
  CHECK((t.velocity(0.25) - Vector3(2, 0, 0)).norm() < 1e-15);
  CHECK((t.velocity(1.0) - Vector3(0, 2, -0.5)).norm() < 1e-15);  // right derivative
  CHECK(t.velocity(3.0).norm() == 0.0);
  CHECK(t.velocity(-0.1).norm() == 0.0);
  // Monotone and random-access lookups agree (the cached segment).
  for (Scalar s : {2.5, 0.1, 2.9, 1.5, 0.0, 2.0}) {
    const Vector3 expected = s <= 1.0 ? Vector3(2 * s, 0, 1)
                                      : Vector3(2, 2 * (s - 1), 1 - 0.5 * (s - 1));
    CHECK((t.position(s) - expected).norm() < 1e-14);
  }
  CHECK(t.knots_in(0.0, 3.0) == std::vector<Scalar>{1.0});
  CHECK(t.knots_in(1.0, 3.0).empty());

  // Refusals.
  CHECK_THROWS_WITH(path({0.0}, {Vector3::Zero()}).validate(),
                    ContainsSubstring("at least two knots"));
  CHECK_THROWS_WITH(path({0.0, 0.0}, {Vector3::Zero(), Vector3::Ones()}).validate(),
                    ContainsSubstring("increase strictly"));
  CHECK_THROWS_WITH(path({0.0, 1.0}, {Vector3::Zero(), Vector3(0, std::nan(""), 0)}).validate(),
                    ContainsSubstring("not finite"));
  CHECK_THROWS_WITH(path({0.0, 1.0}, {Vector3::Zero()}).validate(),
                    ContainsSubstring("point(s)"));

  // CSV: spaces, blank and comment lines tolerated.
  const std::filesystem::path dir =
      std::filesystem::temp_directory_path() / "sparlab_test_forming_csv";
  std::filesystem::create_directories(dir);
  const auto write = [&](const std::string& name, const std::string& text) {
    const std::string p = (dir / name).string();
    std::ofstream(p) << text;
    return p;
  };
  const ToolTrajectory csv = ToolTrajectory::from_csv(
      write("ok.csv", "# a path\n t , x, y ,z\n0, 0,0,0.01\n\n 1.5 ,0.02, 0, 0.01 \n"));
  REQUIRE(csv.times.size() == 2);
  CHECK(csv.times[1] == 1.5);
  CHECK((csv.points[1] - Vector3(0.02, 0, 0.01)).norm() == 0.0);
  CHECK_THROWS_WITH(ToolTrajectory::from_csv(write("h.csv", "time,x,y,z\n0,0,0,0\n1,1,1,1\n")),
                    ContainsSubstring("line 1") && ContainsSubstring("header"));
  CHECK_THROWS_WITH(ToolTrajectory::from_csv(write("f.csv", "t,x,y,z\n0,0,0,0\n1,1,abc,1\n")),
                    ContainsSubstring("line 3") && ContainsSubstring("abc"));
  CHECK_THROWS_WITH(ToolTrajectory::from_csv(write("n.csv", "t,x,y,z\n0,0,0,0\n1,1,1\n")),
                    ContainsSubstring("line 3") && ContainsSubstring("4 values"));
  CHECK_THROWS_WITH(ToolTrajectory::from_csv(write("d.csv", "t,x,y,z\n0,0,0,0\n0,1,1,1\n")),
                    ContainsSubstring("increase strictly"));
  CHECK_THROWS_AS(ToolTrajectory::from_csv((dir / "missing.csv").string()), IoError);
  std::filesystem::remove_all(dir);
}

// ---------------------------------------------------------------------------
// 2. The contact tangent against finite differences
// ---------------------------------------------------------------------------
namespace {

/// The largest error of the contact tangent against central differences of
/// the contact residual, over the DOFs of the nodes in contact, relative to
/// the largest tangent entry.
Scalar tangent_error(const ToolContact& contact, const FemModel& model, const Vector& u,
                     Scalar t, int& active, int& slipping, bool& symmetric) {
  const int dim = model.dim();
  const ToolContactEvaluation ev = contact.evaluate(u, t, true);
  active = 0;
  slipping = 0;
  for (const ToolResultant& r : ev.tools) {
    active += r.active_nodes;
    slipping += r.slipping_nodes;
  }
  symmetric = ev.symmetric;
  Scalar top = 0.0;
  for (const auto& b : ev.tangent) top = std::max(top, b.k.cwiseAbs().maxCoeff());
  const Scalar h = 1.0e-9;
  Scalar worst = 0.0;
  for (const auto& b : ev.tangent) {
    for (int c = 0; c < dim; ++c) {
      const Index d = b.node * dim + c;
      Vector up = u;
      Vector dn = u;
      up(d) += h;
      dn(d) -= h;
      Vector fd = (contact.evaluate(up, t, false).residual -
                   contact.evaluate(dn, t, false).residual) / (2.0 * h);
      // The tangent is node-diagonal: the node's own block, nothing else.
      for (int a = 0; a < dim; ++a) fd(b.node * dim + a) -= b.k(a, c);
      worst = std::max(worst, fd.cwiseAbs().maxCoeff() / top);
    }
  }
  return worst;
}

}  // namespace

TEST_CASE("the penalty contact tangent is the exact derivative of its forces, with friction",
          "[forming]") {
  const IsotropicMaterial m = default_material();
  // A block 10 mm x 10 mm x 2 mm, its top at z = 0, with its nodes moved
  // off the grid so no symmetry helps.
  std::mt19937 rng(7u);
  std::uniform_real_distribution<Scalar> jitter(-1.0e-4, 1.0e-4);
  for (const RigidTool::Shape shape :
       {RigidTool::Shape::Sphere, RigidTool::Shape::Plane, RigidTool::Shape::Cylinder}) {
    for (const std::string regime : {"frictionless", "stick", "slip", "new"}) {
      INFO(to_string(shape) << ", " << regime);
      FemModel model = finalised(hex_block(20, 20, 1, 0.01, 0.01, 0.001, 0.0, 0.0, -0.001), m);
      RigidTool tool;
      tool.name = "tool";
      tool.shape = shape;
      tool.radius = 0.005;
      tool.normal = Vector3(0.1, -0.2, 1.0);          // plane: tilted, facing +z
      tool.axis = Vector3(1.0, 0.3, 0.0);             // cylinder
      tool.surface = box(-kInf, kInf, -kInf, kInf, 0.0, kInf);
      tool.friction = regime == "frictionless" ? 0.0 : regime == "slip" ? 0.1 : 0.8;
      // Centre (or plane point) placed so that several top nodes penetrate,
      // by up to 0.1 mm (the plane, tilted, more).
      const Vector3 centre = shape == RigidTool::Shape::Plane
                                 ? Vector3(0.005, 0.005, -2.0e-5)
                                 : Vector3(0.0052, 0.0047, 0.005 - 1.0e-4);
      tool.trajectory = path({0.0, 1.0}, {centre, centre + Vector3(1e-3, 0, 0)});
      ToolContact contact(model, {tool});
      const Index n = model.dofs().num_dofs();
      Vector u = Vector::Zero(n);
      for (Index d = 0; d < n; ++d) u(d) = 0.02 * jitter(rng);  // up to 2 um
      const Scalar t = 0.0;
      contact.begin_increment(u, t, t, {1}, 1.0e-3);
      // A committed history for every node in contact at (u, t): a force
      // and an anchor so the trial sticks or slips clearly.
      ToolContactEvaluation ev0 = contact.evaluate(u, t, false);
      REQUIRE(ev0.tools[0].active_nodes >= 3);
      if (regime != "frictionless" && regime != "new") {
        ToolHistory h(1);
        const Mesh& mesh = model.mesh();
        for (Index node : contact.slave_nodes(0)) {
          Vector3 x = mesh.node(node) + u.segment(node * 3, 3);
          Vector3 nrm;
          Matrix3 dn;
          Scalar d = 0.0;
          const Vector3 c = tool.trajectory.position(t);
          if (tool.gap(x, c, 3, nrm, dn, d) >= 0.0) continue;
          const Vector3 r = x - c;
          const Vector3 along(1, 0.4, 0);
          const Vector3 tangent = (along - along.dot(nrm) * nrm).normalized();
          ToolNodeHistory entry;
          // Stick: a small force and a short slip; slip: the anchor far away.
          // The committed force leans off the tangent plane (as one of an
          // earlier normal does), which the transport rotates back.
          entry.force = (regime == "stick" ? 5.0 : 50.0) * (tangent + 0.2 * nrm);
          entry.relative = r - (regime == "stick" ? 1.0e-8 : 1.0e-3) * tangent;
          h[0][node] = entry;
        }
        contact.set_history(h);
      }
      int active = 0;
      int slipping = 0;
      bool symmetric = true;
      const Scalar err = tangent_error(contact, model, u, t, active, slipping, symmetric);
      log::info("contact tangent vs central differences, ", to_string(shape), ", ", regime, ": ",
                err, " (", active, " node(s) in contact, ", slipping, " slipping)");
      CHECK(err < 1.0e-6);
      if (regime == "slip") {
        CHECK(slipping == active);
        CHECK_FALSE(symmetric);
      }
      if (regime == "stick" || regime == "new") CHECK(slipping == 0);
      if (regime == "frictionless") CHECK(symmetric);
    }
  }

  // 2-D: a circle on a plane-strain strip, sliding with friction.
  {
    StructuredMeshSpec spec;
    spec.nx = 24;
    spec.ny = 2;
    spec.lx = 0.012;
    spec.ly = 0.002;
    spec.y0 = -0.002;
    FemModel model = finalised(make_structured_quad_mesh(spec), m, 0.01, StressState::PlaneStrain);
    RigidTool tool;
    tool.radius = 0.004;
    tool.surface = box(-kInf, kInf, 0.0, kInf);
    tool.friction = 0.2;
    const Vector3 c(0.0061, 0.004 - 1.0e-4, 0.0);
    tool.trajectory = path({0.0, 1.0}, {c, c + Vector3(1e-3, 0, 0)});
    ToolContact contact(model, {tool});
    Vector u = Vector::Zero(model.dofs().num_dofs());
    contact.begin_increment(u, 0.0, 0.0, {1}, 1.0e-3);
    for (Index d = 0; d < u.size(); d += 2) u(d) = 3.0e-5;  // slides the strip along x
    int active = 0;
    int slipping = 0;
    bool symmetric = true;
    const Scalar err2 = tangent_error(contact, model, u, 0.0, active, slipping, symmetric);
    log::info("contact tangent vs central differences, 2-D circle, slip: ", err2, " (", active,
              " node(s) in contact, ", slipping, " slipping)");
    CHECK(err2 < 1.0e-6);
    CHECK(active >= 2);
    CHECK(slipping == active);
  }
}

TEST_CASE("tool contact refuses what it cannot model", "[forming]") {
  const IsotropicMaterial m = default_material();
  FemModel model = finalised(hex_block(2, 2, 1, 0.01, 0.01, 0.002), m);
  RigidTool tool;
  tool.radius = 0.005;
  tool.surface = box(-kInf, kInf, -kInf, kInf, 0.002, kInf);
  tool.trajectory = path({0.0, 1.0}, {Vector3(0, 0, 0.01), Vector3(0, 0, 0.0)});
  RigidTool bad = tool;
  bad.radius = 0.0;
  CHECK_THROWS_WITH(ToolContact(model, {bad}), ContainsSubstring("positive radius"));
  bad = tool;
  bad.friction = -0.1;
  CHECK_THROWS_WITH(ToolContact(model, {bad}), ContainsSubstring("friction"));
  bad = tool;
  bad.surface = box(-kInf, kInf, -kInf, kInf, 0.5, kInf);
  CHECK_THROWS_WITH(ToolContact(model, {bad}), ContainsSubstring("selects no boundary face"));
  CHECK_THROWS_WITH(ToolContact(model, {tool, tool}), ContainsSubstring("used twice"));
  // A node through the sphere's centre: the increment is cut.
  ToolContact contact(model, {tool});
  Vector u = Vector::Zero(model.dofs().num_dofs());
  CHECK_THROWS_AS(contact.evaluate(u, 1.0, false), SolverError);
}
