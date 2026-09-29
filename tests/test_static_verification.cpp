/// \file test_static_verification.cpp
/// \brief Patch test, stress recovery, beam-theory comparison and mesh
///        convergence.
///
/// The patch test and the stress-recovery test are *verification*: they compare
/// against exact answers for the discrete problem and must hold to round-off.
/// The beam-theory comparison is *validation*: plane-stress elasticity and
/// 1-D beam theory are different models, so the test asserts the size of the
/// expected agreement rather than exactness, and separately checks that the
/// discretisation error itself converges.
#include "TestSupport.hpp"

#include "sparlab/fem/StressRecovery.hpp"
#include "sparlab/mesh/StructuredMesh.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <vector>

using namespace sparlab;
using namespace sparlab::testing;
using Catch::Approx;

namespace {

/// Build a patch-test model on a (possibly distorted) mesh with the exact
/// linear field prescribed on its boundary.
struct PatchResult {
  Scalar displacement_error = 0.0;
  Scalar strain_error = 0.0;
  Scalar stress_error = 0.0;
};

PatchResult run_patch_test(Scalar perturbation) {
  StructuredMeshSpec spec;
  spec.nx = 4;
  spec.ny = 3;
  spec.lx = 2.0;
  spec.ly = 1.5;

  const Vector2 offset(1.0e-4, -2.0e-4);
  Matrix2 gradient;
  gradient << 3.0e-4, 1.0e-4, 1.0e-4, -2.0e-4;
  const Vector3 exact_strain(gradient(0, 0), gradient(1, 1),
                             gradient(0, 1) + gradient(1, 0));

  const IsotropicMaterial material(200.0e9, 0.3, 7850.0, "patch");
  FemModel model(make_perturbed_quad_mesh(spec, perturbation, 7u), material, 0.02,
                 StressState::PlaneStress, IntegrationOptions());

  const std::vector<Mesh::BoundaryFace> edges = model.mesh().boundary_faces();
  std::vector<char> on_boundary(static_cast<std::size_t>(model.mesh().num_nodes()), 0);
  for (const Mesh::BoundaryFace& e : edges) {
    for (Index n : e.nodes) on_boundary[static_cast<std::size_t>(n)] = 1;
  }
  std::vector<Index> boundary_nodes;
  for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
    if (on_boundary[static_cast<std::size_t>(n)]) boundary_nodes.push_back(n);
  }

  DisplacementConstraint bc;
  bc.region.name = "patch_boundary";
  Selector sel;
  sel.kind = SelectorKind::NodeIds;
  sel.ids = boundary_nodes;
  bc.region.members.push_back(sel);
  bc.fix_x = true;
  bc.fix_y = true;
  model.constraints().push_back(bc);

  // The patch test is driven entirely by the prescribed boundary field.
  LoadCaseSpec load;
  load.name = "patch";
  load.prescribed_displacement_only = true;
  model.load_case_specs().push_back(load);
  model.finalize();

  for (Index n : boundary_nodes) {
    const Vector2 x = model.mesh().node(n).head<2>();
    const Vector2 u = offset + gradient * x;
    model.dofs().prescribe(n, 0, u.x());
    model.dofs().prescribe(n, 1, u.y());
  }

  Assembler assembler(model);
  StaticAnalysisOptions options;
  options.linear.residual_tolerance = 1.0e-9;
  StaticAnalysis analysis(model, assembler, options);
  const Vector u = analysis.solve_all().front().displacement;

  PatchResult result;
  Scalar u_scale = 0.0;
  for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
    const Vector2 x = model.mesh().node(n).head<2>();
    const Vector2 expected = offset + gradient * x;
    result.displacement_error =
        std::max(result.displacement_error,
                 std::max(std::abs(u(n * 2 + 0) - expected.x()),
                          std::abs(u(n * 2 + 1) - expected.y())));
    u_scale = std::max(u_scale, expected.cwiseAbs().maxCoeff());
  }
  result.displacement_error /= u_scale;

  const StressField field = recover_stresses(model, assembler, u);
  const Vector3 exact_stress = material.plane_stress_matrix() * exact_strain;
  for (Index e = 0; e < model.mesh().num_elements(); ++e) {
    result.strain_error = std::max(
        result.strain_error,
        (field.element_strain.col(e) - exact_strain).cwiseAbs().maxCoeff());
    result.stress_error = std::max(
        result.stress_error,
        (field.element_stress.col(e) - exact_stress).cwiseAbs().maxCoeff());
  }
  result.strain_error /= exact_strain.cwiseAbs().maxCoeff();
  result.stress_error /= exact_stress.cwiseAbs().maxCoeff();
  return result;
}

}  // namespace

TEST_CASE("constant-strain patch test passes on uniform and distorted meshes",
          "[patch][verification]") {
  for (Scalar perturbation : {0.0, 0.15, 0.30, 0.40}) {
    const PatchResult r = run_patch_test(perturbation);
    INFO("perturbation = " << perturbation);
    REQUIRE(r.displacement_error == Approx(0.0).margin(1.0e-11));
    REQUIRE(r.strain_error == Approx(0.0).margin(1.0e-11));
    REQUIRE(r.stress_error == Approx(0.0).margin(1.0e-11));
  }
}

TEST_CASE("stress recovery reproduces a known displacement field",
          "[stress][verification]") {
  // Impose a known quadratic-free (linear) field directly, bypassing the solve,
  // and check every recovered quantity against the closed form.
  StructuredMeshSpec spec;
  spec.nx = 3;
  spec.ny = 2;
  spec.lx = 0.6;
  spec.ly = 0.4;

  const IsotropicMaterial material(70.0e9, 0.3, 2700.0, "known_field");
  FemModel model(make_perturbed_quad_mesh(spec, 0.2, 3u), material, 0.01,
                 StressState::PlaneStress, IntegrationOptions());

  DisplacementConstraint bc;
  Selector box;
  box.kind = SelectorKind::Box;
  box.xmax = 0.0;
  bc.region.members.push_back(box);
  bc.fix_x = true;
  bc.fix_y = true;
  model.constraints().push_back(bc);
  LoadCaseSpec load;
  load.name = "unused";
  PointLoadSpec p;
  Selector nearest;
  nearest.kind = SelectorKind::NearestNode;
  nearest.point = Vector3(spec.lx, 0.0, 0.0);
  p.region.members.push_back(nearest);
  p.force = Vector3(1.0, 0.0, 0.0);
  load.point_loads.push_back(p);
  model.load_case_specs().push_back(load);
  model.finalize();

  const Scalar exx = 1.5e-4;
  const Scalar eyy = -0.5e-4;
  const Scalar gxy = 2.0e-4;
  Vector u = Vector::Zero(model.dofs().num_dofs());
  for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
    const Vector3 x = model.mesh().node(n);
    u(n * 2 + 0) = exx * x.x() + 0.5 * gxy * x.y();
    u(n * 2 + 1) = 0.5 * gxy * x.x() + eyy * x.y();
  }

  Assembler assembler(model);
  const StressField field = recover_stresses(model, assembler, u);

  const Vector3 exact_strain(exx, eyy, gxy);
  const Vector3 exact_stress = material.plane_stress_matrix() * exact_strain;
  const Scalar exact_vm =
      von_mises(exact_stress, StressState::PlaneStress, material.poisson_ratio());

  for (Index e = 0; e < model.mesh().num_elements(); ++e) {
    REQUIRE((field.element_strain.col(e) - exact_strain).cwiseAbs().maxCoeff() /
                exact_strain.cwiseAbs().maxCoeff() == Approx(0.0).margin(1.0e-12));
    REQUIRE((field.element_stress.col(e) - exact_stress).cwiseAbs().maxCoeff() /
                exact_stress.cwiseAbs().maxCoeff() == Approx(0.0).margin(1.0e-12));
    REQUIRE(field.element_von_mises(e) == Approx(exact_vm).epsilon(1.0e-11));
    // Principal-stress invariants.
    REQUIRE(field.element_principal_max(e) + field.element_principal_min(e) ==
            Approx(exact_stress(0) + exact_stress(1)).epsilon(1.0e-11));
  }

  // Nodal averaging of a constant field reproduces the constant exactly.
  for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
    REQUIRE((field.nodal_stress.col(n) - exact_stress).cwiseAbs().maxCoeff() /
                exact_stress.cwiseAbs().maxCoeff() == Approx(0.0).margin(1.0e-11));
  }

  // Total strain energy equals 1/2 eps^T sigma * volume for a uniform field.
  const Scalar volume = model.domain_volume();
  const Scalar expected_energy = 0.5 * exact_strain.dot(exact_stress) * volume;
  REQUIRE(field.element_strain_energy.sum() ==
          Approx(expected_energy).epsilon(1.0e-10));

  // The scaled variant multiplies the stress but not the strain.
  const Vector scale = Vector::Constant(model.mesh().num_elements(), 0.25);
  const StressField scaled = recover_stresses(model, assembler, u, &scale);
  REQUIRE(scaled.element_stress.col(0).isApprox(0.25 * exact_stress, 1.0e-10));
  REQUIRE(scaled.element_solid_stress.col(0).isApprox(exact_stress, 1.0e-10));
  REQUIRE(scaled.element_strain.col(0).isApprox(exact_strain, 1.0e-10));
}

TEST_CASE("element and nodal stress fields agree with pointwise evaluation",
          "[stress]") {
  FemModel model = make_small_plate(5, 4);
  Assembler assembler(model);
  StaticAnalysis analysis(model, assembler, StaticAnalysisOptions());
  const Vector u = analysis.solve_all().front().displacement;
  const StressField field = recover_stresses(model, assembler, u);

  // The element value is the average over the stiffness quadrature points;
  // reproduce it independently.
  const std::vector<NaturalPoint> points =
      model.element().stress_evaluation_points(model.integration());
  for (Index e = 0; e < std::min<Index>(model.mesh().num_elements(), 8); ++e) {
    Vector3 average = Vector3::Zero();
    for (const NaturalPoint& p : points) {
      average += element_stress_at(model, e, p, u);
    }
    average /= static_cast<Scalar>(points.size());
    REQUIRE((field.element_stress.col(e) - average).cwiseAbs().maxCoeff() /
                std::max(average.cwiseAbs().maxCoeff(), 1.0) ==
            Approx(0.0).margin(1.0e-12));
  }
}

TEST_CASE("cantilever tip deflection matches beam theory",
          "[beam][validation]") {
  // nu = 0 removes the Poisson coupling that 1-D beam theory omits, leaving
  // shear deformation as the only difference; Timoshenko theory accounts for
  // it, so the agreement should be tight.
  CantileverCase c;
  c.poisson = 0.0;
  FemModel model = make_cantilever(c, 160, 16);
  Assembler assembler(model);
  StaticAnalysis analysis(model, assembler, StaticAnalysisOptions());
  const StaticSolution sol = analysis.solve_all().front();
  const Scalar tip = tip_deflection(model, sol.displacement);

  const Scalar eb = c.euler_bernoulli_tip();
  const Scalar ti = c.timoshenko_tip();

  INFO("FEM " << tip << " m, Euler-Bernoulli " << eb << " m, Timoshenko " << ti << " m");
  // Timoshenko is the closer reference for a 2-D model of an L/h = 10 beam.
  REQUIRE(std::abs(tip - ti) / std::abs(ti) < 0.02);
  // Euler-Bernoulli under-predicts the deflection because it omits shear.
  REQUIRE(std::abs(tip) > std::abs(eb));
  REQUIRE(std::abs(tip - eb) / std::abs(eb) < 0.06);
}

TEST_CASE("cantilever deflection converges under mesh refinement",
          "[convergence][verification]") {
  CantileverCase c;
  c.poisson = 0.0;

  const std::vector<std::pair<Index, Index>> meshes = {
      {10, 2}, {20, 4}, {40, 8}, {80, 16}};
  std::vector<Scalar> tips;
  std::vector<Scalar> sizes;
  std::vector<Scalar> compliances;

  for (const auto& m : meshes) {
    FemModel model = make_cantilever(c, m.first, m.second);
    Assembler assembler(model);
    StaticAnalysis analysis(model, assembler, StaticAnalysisOptions());
    const StaticSolution sol = analysis.solve_all().front();
    tips.push_back(tip_deflection(model, sol.displacement));
    compliances.push_back(sol.compliance);
    sizes.push_back(c.length / static_cast<Scalar>(m.first));
  }

  // Monotone convergence towards the finest value, with the error shrinking.
  const Scalar finest = tips.back();
  std::vector<Scalar> errors;
  for (std::size_t i = 0; i + 1 < tips.size(); ++i) {
    errors.push_back(std::abs(tips[i] - finest) / std::abs(finest));
  }
  for (std::size_t i = 0; i + 1 < errors.size(); ++i) {
    INFO("error[" << i << "] = " << errors[i] << ", error[" << i + 1
                  << "] = " << errors[i + 1]);
    REQUIRE(errors[i + 1] < errors[i]);
  }

  // Observed convergence order between the two finest refinements: the Q4
  // displacement error is expected to be second order in h.
  const std::size_t last = errors.size() - 1;
  const Scalar order = std::log(errors[last - 1] / errors[last]) /
                       std::log(sizes[last - 1] / sizes[last]);
  INFO("observed order " << order);
  REQUIRE(order > 1.5);
  REQUIRE(order < 2.6);

  // Compliance converges from below: a coarser mesh is stiffer, and the
  // displacement-based Q4 is a lower bound on the true compliance.
  for (std::size_t i = 0; i + 1 < compliances.size(); ++i) {
    REQUIRE(compliances[i] < compliances[i + 1]);
  }
}

TEST_CASE("plane strain is stiffer than plane stress for the same beam",
          "[validation]") {
  CantileverCase c;
  c.poisson = 0.3;

  const auto solve = [&](StressState state) {
    StructuredMeshSpec spec;
    spec.nx = 40;
    spec.ny = 8;
    spec.lx = c.length;
    spec.ly = c.height;
    FemModel model(make_structured_quad_mesh(spec),
                   IsotropicMaterial(c.youngs, c.poisson, c.density), c.thickness,
                   state, IntegrationOptions());
    DisplacementConstraint root;
    Selector box;
    box.kind = SelectorKind::Box;
    box.xmax = 0.0;
    root.region.members.push_back(box);
    root.fix_x = true;
    root.fix_y = true;
    model.constraints().push_back(root);
    LoadCaseSpec load;
    load.name = "tip";
    PointLoadSpec tip;
    Selector tip_box;
    tip_box.kind = SelectorKind::Box;
    tip_box.xmin = c.length;
    tip.region.members.push_back(tip_box);
    tip.force = Vector3(0.0, c.tip_load, 0.0);
    load.point_loads.push_back(tip);
    model.load_case_specs().push_back(load);
    model.finalize();
    Assembler assembler(model);
    StaticAnalysis analysis(model, assembler, StaticAnalysisOptions());
    return analysis.solve_all().front().compliance;
  };

  const Scalar plane_stress = solve(StressState::PlaneStress);
  const Scalar plane_strain = solve(StressState::PlaneStrain);
  REQUIRE(plane_strain < plane_stress);
}

TEST_CASE("distributed tractions give a mesh-independent resultant",
          "[loads][verification]") {
  const Scalar pressure = -2.0e5;  // Pa
  const Scalar thickness = 0.01;
  const Scalar length = 1.0;

  Scalar previous = 0.0;
  for (Index nx : {8, 16, 32}) {
    StructuredMeshSpec spec;
    spec.nx = nx;
    spec.ny = nx / 4;
    spec.lx = length;
    spec.ly = 0.25;

    FemModel model(make_structured_quad_mesh(spec), default_material(), thickness,
                   StressState::PlaneStress, IntegrationOptions());
    DisplacementConstraint bc;
    Selector box;
    box.kind = SelectorKind::Box;
    box.xmax = 0.0;
    bc.region.members.push_back(box);
    bc.fix_x = true;
    bc.fix_y = true;
    model.constraints().push_back(bc);

    LoadCaseSpec load;
    load.name = "top_pressure";
    TractionLoadSpec traction;
    traction.region.name = "top_edge";
    Selector top;
    top.kind = SelectorKind::Box;
    top.ymin = spec.ly;
    traction.region.members.push_back(top);
    traction.traction = Vector3(0.0, pressure, 0.0);
    load.tractions.push_back(traction);
    model.load_case_specs().push_back(load);
    model.finalize();

    // The consistent nodal forces must sum to pressure * length * thickness.
    const Vector& f = model.load_vectors().front();
    Scalar total = 0.0;
    for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
      total += f(n * 2 + 1);
    }
    REQUIRE(total == Approx(pressure * length * thickness).epsilon(1.0e-12));

    Assembler assembler(model);
    StaticAnalysis analysis(model, assembler, StaticAnalysisOptions());
    const StaticSolution sol = analysis.solve_all().front();
    REQUIRE(sol.equilibrium.relative_force_error == Approx(0.0).margin(1.0e-10));

    if (previous != 0.0) {
      // Refinement changes the answer by only a few percent.
      REQUIRE(std::abs(sol.compliance - previous) / previous < 0.15);
    }
    previous = sol.compliance;
  }
}

TEST_CASE("a prescribed rigid motion passes the equilibrium check", "[static][equilibrium]") {
  // A free plate whose left edge is moved rigidly: nothing is strained, the
  // reactions are round-off, and their balance is judged against the gross
  // force they are formed from rather than against themselves.
  StructuredMeshSpec spec;
  spec.nx = 6;
  spec.ny = 3;
  spec.lx = 0.6;
  spec.ly = 0.3;
  FemModel model(make_structured_quad_mesh(spec), default_material(), 0.01,
                 StressState::PlaneStrain, IntegrationOptions());
  DisplacementConstraint shaken;
  Selector left;
  left.kind = SelectorKind::Box;
  left.xmax = 0.0;
  shaken.region.members.push_back(left);
  shaken.set(0, true, 2.0e-4);
  shaken.set(1, true, -1.0e-3);
  model.constraints().push_back(shaken);
  LoadCaseSpec lc;
  lc.name = "translate";
  lc.prescribed_displacement_only = true;
  model.load_case_specs().push_back(lc);
  model.finalize();
  Assembler assembler(model);
  StaticAnalysis analysis(model, assembler, StaticAnalysisOptions());
  const StaticSolution s = analysis.solve_all().front();
  REQUIRE(s.equilibrium.relative_force_error <= 1.0e-12);
  // The body translates rigidly.
  for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
    REQUIRE(s.displacement(2 * n) == Approx(2.0e-4).epsilon(1e-9));
    REQUIRE(s.displacement(2 * n + 1) == Approx(-1.0e-3).epsilon(1e-9));
  }
}
