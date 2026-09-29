/// \file sparlab_verify.cpp
/// \brief Verification and validation studies.
///
/// **Verification** asks "does the code solve the equations it claims to
/// solve?" and is judged against exact answers for the discrete problem:
/// the patch test, solver agreement, and the finite-difference check of the
/// analytical sensitivities. These must pass to tight tolerances.
///
/// **Validation** asks "are the modelling assumptions appropriate?" and is
/// judged against an independent theory whose assumptions differ from the
/// model's: Euler-Bernoulli and Timoshenko beam deflection, and
/// Euler-Bernoulli natural frequencies. A finite gap there is expected and is
/// reported as a modelling difference, not an error. The studies quantify it
/// rather than asserting agreement.
///
/// Studies implemented:
///   * `sensitivity`      analytical vs central-difference topology gradients
///                        over a range of perturbation sizes;
///   * `mesh-convergence` cantilever tip deflection and compliance vs mesh
///                        size, with observed convergence order, against
///                        Euler-Bernoulli and Timoshenko theory at two
///                        Poisson ratios;
///   * `modal`            cantilever bending frequencies vs mesh size against
///                        Euler-Bernoulli theory;
///   * `solver-agreement` every linear solver (four direct, Jacobi CG,
///                        multigrid CG and the automatic choice) on one
///                        small model;
///   * `patch-test`       constant-strain patch test on a distorted mesh;
///   * `patch-test-3d`    the same on a distorted Hex8 mesh;
///   * `mesh-convergence-3d` solid (Hex8) cantilever vs Timoshenko theory;
///   * `sensitivity-3d`   the finite-difference gradient check on a hex mesh;
///   * `modal-3d`         solid cantilever frequencies about both axes;
///   * `patch-test-simplex` the patch test on distorted Tri3 and Tet4 meshes;
///   * `mesh-convergence-simplex` Tri3 and Tet4 cantilevers vs Timoshenko;
///   * `multigrid`        multigrid CG against Cholesky and Jacobi CG under
///                        refinement (Hex8 and Tet4);
///   * `sensitivity-projection` the gradient check through the Heaviside
///                        projection at three sharpnesses;
///   * `patch-test-quadratic` pure bending on distorted Tet10 meshes (exact)
///                        against Tet4 and Hex8;
///   * `mesh-convergence-tet10` Tet10, Tet4 and Hex8 cantilevers vs Timoshenko;
///   * `buckling-euler`   buckling of a cantilever column (Q4, Hex8, Tet4, Tet10)
///                        vs Euler with Engesser's shear correction;
///   * `sensitivity-buckling` the adjoint gradient of the buckling constraint;
///   * `sensitivity-overhang` the compliance gradient through the overhang filter.
///
/// Studies of the pressure, volume and thermal loads (verify_loads.cpp):
///   * `lame-cylinder`    thick cylinder under internal pressure vs Lame;
///   * `rotating-disk`    spinning disk and cylinder vs the exact solutions;
///   * `thermal-cylinder` conduction with generation and convection, and the
///                        thermal stress it causes, vs the exact solutions;
///   * `bimetal-strip`    curvature of a heated two-material strip vs
///                        Timoshenko;
///   * `self-weight`      a bar hanging under its own weight vs the exact field.
///
/// Studies of the geometrically non-linear analysis (verify_nonlinear.cpp):
///   * `elastica`         a slender Tet10 cantilever under a dead tip force
///                        vs Euler's elastica, up to k = PL^2/EI = 10;
///   * `hyperelastic-cylinder` a thick tube at finite strain vs the exact
///                        solutions: neo-Hookean, inflated by a follower
///                        pressure and spinning; Saint Venant-Kirchhoff,
///                        heated by a conducted temperature;
///   * `arch-snap-through` a shallow arch's snap-through: arc length vs
///                        displacement control, inertia, and load control
///                        stopping at the limit point.

#include "AppSupport.hpp"
#include "VerifySupport.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Timer.hpp"
#include "sparlab/fem/Assembler.hpp"
#include "sparlab/fem/Buckling.hpp"
#include "sparlab/fem/ModalAnalysis.hpp"
#include "sparlab/fem/ModelDiagnostics.hpp"
#include "sparlab/fem/StaticAnalysis.hpp"
#include "sparlab/fem/StressRecovery.hpp"
#include "sparlab/io/Config.hpp"
#include "sparlab/io/CsvWriter.hpp"
#include "sparlab/io/Json.hpp"
#include "sparlab/mesh/StructuredMesh.hpp"
#include "sparlab/topopt/BucklingConstraint.hpp"
#include "sparlab/topopt/OverhangFilter.hpp"
#include "sparlab/topopt/Projection.hpp"
#include "sparlab/topopt/Sensitivity.hpp"
#include "sparlab/topopt/TopologyOptimizer.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>

using namespace sparlab;
using verify::observed_order;
using verify::StudyOutcome;

namespace {

// ---------------------------------------------------------------------------
// Shared cantilever definition used by the validation studies.
//
// L/h = 10 keeps the beam slender enough for beam theory to be a reasonable
// reference while remaining a genuinely two-dimensional model. The tip load is
// specified as a resultant spread over the tip edge nodes so the total force is
// mesh independent.
// ---------------------------------------------------------------------------
struct CantileverSpec {
  Scalar length = 1.0;        ///< L [m]
  Scalar height = 0.1;        ///< h [m]
  Scalar thickness = 0.01;    ///< t [m]
  Scalar youngs = 70.0e9;     ///< E [Pa]
  Scalar poisson = 0.3;       ///< nu [-]
  Scalar density = 2700.0;    ///< rho [kg/m^3]
  Scalar tip_load = -1000.0;  ///< P [N], negative = downward
};

struct CantileverResult {
  Index nx = 0;
  Index ny = 0;
  Index num_elements = 0;
  Index num_dofs = 0;
  Scalar element_size = 0.0;
  Scalar tip_deflection_mean = 0.0;   ///< mean u_y over the tip edge [m]
  Scalar tip_deflection_mid = 0.0;    ///< u_y at the tip mid-height node [m]
  Scalar compliance = 0.0;
  Scalar strain_energy = 0.0;
  Scalar max_von_mises = 0.0;
  Scalar reaction_error = 0.0;
  Scalar solve_seconds = 0.0;
};

/// Euler-Bernoulli tip deflection of a cantilever with an end load:
/// delta = P L^3 / (3 E I), I = t h^3 / 12.
Scalar euler_bernoulli_tip(const CantileverSpec& s) {
  const Scalar i = s.thickness * s.height * s.height * s.height / 12.0;
  return s.tip_load * s.length * s.length * s.length / (3.0 * s.youngs * i);
}

/// Timoshenko tip deflection: adds the shear term P L / (k G A) with the
/// rectangular-section shear coefficient k = 5/6.
Scalar timoshenko_tip(const CantileverSpec& s) {
  const Scalar i = s.thickness * s.height * s.height * s.height / 12.0;
  const Scalar a = s.thickness * s.height;
  const Scalar g = s.youngs / (2.0 * (1.0 + s.poisson));
  const Scalar k = 5.0 / 6.0;
  return s.tip_load * s.length * s.length * s.length / (3.0 * s.youngs * i) +
         s.tip_load * s.length / (k * g * a);
}

FemModel build_cantilever(const CantileverSpec& spec, Index nx, Index ny) {
  StructuredMeshSpec mesh_spec;
  mesh_spec.nx = nx;
  mesh_spec.ny = ny;
  mesh_spec.lx = spec.length;
  mesh_spec.ly = spec.height;

  IsotropicMaterial material(spec.youngs, spec.poisson, spec.density, "verification");
  FemModel model(make_structured_quad_mesh(mesh_spec), material, spec.thickness,
                 StressState::PlaneStress, IntegrationOptions());

  DisplacementConstraint root;
  root.region.name = "root";
  Selector box;
  box.kind = SelectorKind::Box;
  box.xmax = 0.0;
  root.region.members.push_back(box);
  root.fix_x = true;
  root.fix_y = true;
  model.constraints().push_back(root);

  LoadCaseSpec load;
  load.name = "tip_load";
  PointLoadSpec tip;
  tip.region.name = "tip_edge";
  Selector tip_box;
  tip_box.kind = SelectorKind::Box;
  tip_box.xmin = spec.length;
  tip.region.members.push_back(tip_box);
  tip.force = Vector3(0.0, spec.tip_load, 0.0);
  tip.distribute_total = true;
  load.point_loads.push_back(tip);
  model.load_case_specs().push_back(load);

  model.finalize();
  return model;
}

CantileverResult solve_cantilever(const CantileverSpec& spec, Index nx, Index ny) {
  CantileverResult out;
  out.nx = nx;
  out.ny = ny;
  FemModel model = build_cantilever(spec, nx, ny);
  Assembler assembler(model);

  Timer timer;
  StaticAnalysisOptions options;
  StaticAnalysis analysis(model, assembler, options);
  const std::vector<StaticSolution> solutions = analysis.solve_all();
  out.solve_seconds = timer.elapsed_seconds();

  const StaticSolution& sol = solutions.front();
  const StressField field = recover_stresses(model, assembler, sol.displacement);

  out.num_elements = model.mesh().num_elements();
  out.num_dofs = model.dofs().num_dofs();
  out.element_size = spec.length / static_cast<Scalar>(nx);
  out.compliance = sol.compliance;
  out.strain_energy = sol.strain_energy;
  out.max_von_mises = field.element_von_mises.maxCoeff();
  out.reaction_error = sol.equilibrium.relative_force_error;

  // Average u_y over the tip edge, plus the mid-height value.
  const StructuredGridInfo& info = *model.mesh().structured_info();
  Scalar sum = 0.0;
  for (Index j = 0; j <= ny; ++j) {
    const Index node = structured_node_index(info, nx, j);
    sum += sol.displacement(node * 2 + 1);
  }
  out.tip_deflection_mean = sum / static_cast<Scalar>(ny + 1);
  const Index mid = structured_node_index(info, nx, ny / 2);
  out.tip_deflection_mid = sol.displacement(mid * 2 + 1);
  return out;
}


// ---------------------------------------------------------------------------
// Studies
// ---------------------------------------------------------------------------

/// Analytical topology sensitivities vs central differences.
StudyOutcome study_sensitivity(const std::string& out_dir, json::Value& summary,
                               Scalar tolerance) {
  // A deliberately small mesh: the check needs two extra linear solves per
  // element per step size.
  CantileverSpec spec;
  spec.length = 0.6;
  spec.height = 0.3;
  spec.tip_load = -500.0;
  FemModel model = build_cantilever(spec, 12, 6);
  Assembler assembler(model);

  // Radius of 1.5 cells: large enough for the filter to couple neighbours (so
  // the chain rule is genuinely exercised), small enough to stay cheap.
  const Scalar cell = spec.length / 12.0;
  DensityFilter filter(model.mesh(), FilterType::Density, 1.5 * cell);

  std::vector<PassiveRegionSpec> passive;
  {
    // One passive solid and one passive void patch, so the exclusion logic is
    // exercised and reported.
    PassiveRegionSpec solid;
    solid.region.name = "passive_solid_patch";
    Selector s;
    s.kind = SelectorKind::Box;
    s.xmax = 0.05;
    s.ymin = 0.10;
    s.ymax = 0.20;
    solid.region.members.push_back(s);
    solid.solid = true;
    passive.push_back(solid);

    PassiveRegionSpec hole;
    hole.region.name = "passive_void_patch";
    Selector c;
    c.kind = SelectorKind::Circle;
    c.center = Vector3(0.30, 0.15, 0.0);
    c.radius = 0.04;
    hole.region.members.push_back(c);
    hole.solid = false;
    passive.push_back(hole);
  }

  DesignDomain domain(model, /*volume_fraction=*/0.5, /*initial_density=*/0.5, passive);

  SimpOptions simp;
  simp.penalty = 3.0;
  simp.emin_ratio = 1.0e-9;
  StaticAnalysisOptions analysis_options;
  // The FD check differences two compliance values that agree to ~10 digits, so
  // the linear solves must be as accurate as the direct solver allows.
  analysis_options.linear.type = LinearSolverType::SimplicialLdlt;
  analysis_options.linear.residual_tolerance = 1.0e-9;

  ComplianceObjective objective(model, assembler, filter, domain, simp,
                                analysis_options);

  // A non-uniform starting design: a uniform field makes every sensitivity
  // nearly identical and hides indexing mistakes.
  Vector x = domain.initial_design();
  for (Index e = 0; e < domain.num_elements(); ++e) {
    if (!domain.is_free(e)) continue;
    const Vector3 c = model.mesh().element_centroid(e);
    x(e) = 0.35 + 0.30 * std::sin(7.0 * c.x()) * std::cos(5.0 * c.y());
  }
  domain.clamp(x);

  // Every element is tested: the mesh is small enough that two extra linear
  // solves per element per step size is still fast, and full coverage is what
  // makes an indexing mistake impossible to miss.
  std::vector<Index> targets;
  for (Index e = 0; e < domain.num_elements(); ++e) targets.push_back(e);

  const std::vector<Scalar> steps = {1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4,
                                     1.0e-5, 1.0e-6, 1.0e-7, 1.0e-8};

  CsvWriter step_csv(path_join(out_dir, "sensitivity_steps.csv"),
                     {"step[-]", "num_tested", "num_excluded", "max_absolute_error[J]",
                      "max_relative_error[-]", "rms_relative_error[-]",
                      "directional_relative_error[-]", "gradient_inf_norm[J]",
                      "passed[-]"});
  CsvWriter element_csv(
      path_join(out_dir, "sensitivity_elements.csv"),
      {"step[-]", "element", "cx[m]", "cy[m]", "analytical[J]", "finite_difference[J]",
       "absolute_error[J]", "relative_error[-]", "excluded[-]"});

  json::Value step_records = json::Value::make_array();
  Scalar best_error = std::numeric_limits<Scalar>::max();
  Scalar best_step = 0.0;
  Index excluded_count = 0;
  std::vector<std::string> exclusion_reasons;

  for (Scalar step : steps) {
    const SensitivityCheckResult check =
        verify_sensitivities(objective, domain, x, targets, step, tolerance);
    step_csv.row({step, static_cast<Scalar>(check.num_tested),
                  static_cast<Scalar>(check.num_excluded), check.max_absolute_error,
                  check.max_relative_error, check.rms_relative_error,
                  check.directional_relative_error, check.gradient_infinity_norm,
                  check.passed ? 1.0 : 0.0});
    for (const SensitivityCheckEntry& entry : check.entries) {
      const Vector3 c = model.mesh().element_centroid(entry.element);
      element_csv.row({step, static_cast<Scalar>(entry.element), c.x(), c.y(),
                       entry.analytical, entry.finite_difference, entry.absolute_error,
                       entry.relative_error, entry.excluded ? 1.0 : 0.0});
    }
    if (check.max_relative_error < best_error) {
      best_error = check.max_relative_error;
      best_step = step;
    }
    excluded_count = check.num_excluded;
    if (exclusion_reasons.empty()) {
      for (const SensitivityCheckEntry& entry : check.entries) {
        if (!entry.excluded) continue;
        if (std::find(exclusion_reasons.begin(), exclusion_reasons.end(),
                      entry.exclusion_reason) == exclusion_reasons.end()) {
          exclusion_reasons.push_back(entry.exclusion_reason);
        }
      }
    }

    json::Value rec = json::Value::make_object();
    rec.set("step", json::Value::make_number(step));
    rec.set("num_tested", json::Value::make_number(check.num_tested));
    rec.set("num_excluded", json::Value::make_number(check.num_excluded));
    rec.set("max_absolute_error_J",
            json::Value::make_number(check.max_absolute_error));
    rec.set("max_relative_error", json::Value::make_number(check.max_relative_error));
    rec.set("rms_relative_error", json::Value::make_number(check.rms_relative_error));
    rec.set("directional_relative_error",
            json::Value::make_number(check.directional_relative_error));
    rec.set("passed", json::Value::make_bool(check.passed));
    step_records.push_back(rec);
  }
  step_csv.close();
  element_csv.close();

  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("verification"));
  block.set("mesh", json::Value::make_string("12 x 6 Q4, 72 elements"));
  block.set("filter", json::Value::make_string("density, radius = 1.5 cells"));
  block.set("simp_penalty", json::Value::make_number(simp.penalty));
  block.set("design_point",
            json::Value::make_string(
                "non-uniform: x = 0.35 + 0.30 sin(7x) cos(5y) on free elements"));
  block.set("tested_elements", json::Value::make_number(targets.size()));
  block.set("excluded_elements", json::Value::make_number(excluded_count));
  block.set("exclusion_reasons", json::array_of(exclusion_reasons));
  block.set("tolerance", json::Value::make_number(tolerance));
  block.set("best_step", json::Value::make_number(best_step));
  block.set("best_max_relative_error", json::Value::make_number(best_error));
  block.set("steps", step_records);
  summary.set("sensitivity", block);

  StudyOutcome outcome;
  outcome.name = "sensitivity (analytical vs central differences)";
  outcome.kind = "verification";
  outcome.metric = "min over steps of max relative error";
  outcome.value = best_error;
  outcome.tolerance = tolerance;
  outcome.passed = best_error <= tolerance;
  std::ostringstream note;
  note << "best at h = " << best_step << "; " << excluded_count << " of "
       << targets.size() << " candidates excluded at bounds";
  outcome.note = note.str();
  return outcome;
}

/// Cantilever tip deflection and compliance vs mesh size.
StudyOutcome study_mesh_convergence(const std::string& out_dir, json::Value& summary) {
  CsvWriter csv(path_join(out_dir, "mesh_convergence.csv"),
                {"poisson[-]", "nx", "ny", "num_elements", "num_dofs", "h[m]",
                 "tip_mean[m]", "tip_mid[m]", "euler_bernoulli[m]", "timoshenko[m]",
                 "rel_error_eb[-]", "rel_error_timoshenko[-]", "compliance[J]",
                 "strain_energy[J]", "compliance_over_2U[-]", "max_von_mises[Pa]",
                 "reaction_rel_error[-]", "solve_seconds[s]"});

  json::Value records = json::Value::make_array();
  const std::vector<Scalar> poissons = {0.0, 0.3};
  const std::vector<std::pair<Index, Index>> meshes = {
      {10, 2}, {20, 4}, {40, 8}, {80, 16}, {160, 32}, {320, 64}};

  Scalar worst_richardson = 0.0;
  std::map<Scalar, Scalar> order_by_poisson;
  std::map<Scalar, Scalar> finest_error_by_poisson;

  for (Scalar nu : poissons) {
    CantileverSpec spec;
    spec.poisson = nu;
    const Scalar eb = euler_bernoulli_tip(spec);
    const Scalar ti = timoshenko_tip(spec);

    std::vector<CantileverResult> results;
    for (const auto& m : meshes) {
      const CantileverResult r = solve_cantilever(spec, m.first, m.second);
      results.push_back(r);
      csv.row({nu, static_cast<Scalar>(r.nx), static_cast<Scalar>(r.ny),
               static_cast<Scalar>(r.num_elements), static_cast<Scalar>(r.num_dofs),
               r.element_size, r.tip_deflection_mean, r.tip_deflection_mid, eb, ti,
               std::abs(r.tip_deflection_mean - eb) / std::abs(eb),
               std::abs(r.tip_deflection_mean - ti) / std::abs(ti), r.compliance,
               r.strain_energy,
               r.strain_energy != 0.0 ? r.compliance / (2.0 * r.strain_energy) : 0.0,
               r.max_von_mises, r.reaction_error, r.solve_seconds});

      json::Value rec = json::Value::make_object();
      rec.set("poisson", json::Value::make_number(nu));
      rec.set("nx", json::Value::make_number(r.nx));
      rec.set("ny", json::Value::make_number(r.ny));
      rec.set("num_dofs", json::Value::make_number(r.num_dofs));
      rec.set("h_m", json::Value::make_number(r.element_size));
      rec.set("tip_mean_m", json::Value::make_number(r.tip_deflection_mean));
      rec.set("compliance_J", json::Value::make_number(r.compliance));
      rec.set("relative_error_vs_timoshenko",
              json::Value::make_number(std::abs(r.tip_deflection_mean - ti) /
                                       std::abs(ti)));
      records.push_back(rec);
    }

    // Self-convergence: use the finest mesh as the reference for the observed
    // order, which measures discretisation error alone (verification) rather
    // than the beam-theory gap (validation).
    const CantileverResult& finest = results.back();
    std::vector<Scalar> errors;
    std::vector<Scalar> sizes;
    for (std::size_t i = 0; i + 1 < results.size(); ++i) {
      errors.push_back(std::abs(results[i].tip_deflection_mean -
                                finest.tip_deflection_mean) /
                       std::abs(finest.tip_deflection_mean));
      sizes.push_back(results[i].element_size);
    }
    Scalar order = 0.0;
    if (errors.size() >= 2) {
      order = observed_order(sizes[errors.size() - 2], errors[errors.size() - 2],
                             sizes[errors.size() - 1], errors[errors.size() - 1]);
    }
    order_by_poisson[nu] = order;
    finest_error_by_poisson[nu] =
        std::abs(finest.tip_deflection_mean - ti) / std::abs(ti);
    worst_richardson = std::max(worst_richardson, errors.empty() ? 0.0 : errors.back());
  }
  csv.close();

  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("verification (self-convergence) + "
                                             "validation (beam theory)"));
  block.set("geometry", json::Value::make_string(
                            "L = 1 m, h = 0.1 m (L/h = 10), t = 0.01 m, "
                            "E = 70 GPa, tip resultant -1000 N"));
  block.set("records", records);
  json::Value orders = json::Value::make_object();
  for (const auto& kv : order_by_poisson) {
    orders.set(app::format(kv.first, 3), json::Value::make_number(kv.second));
  }
  block.set("observed_convergence_order_tip_deflection", orders);
  json::Value finest = json::Value::make_object();
  for (const auto& kv : finest_error_by_poisson) {
    finest.set(app::format(kv.first, 3), json::Value::make_number(kv.second));
  }
  block.set("finest_mesh_relative_error_vs_timoshenko", finest);
  block.set("note",
            json::Value::make_string(
                "The plane-stress model and beam theory are different models. The "
                "Timoshenko value (bending + shear) is the closer reference; the "
                "residual gap at nu = 0.3 is the anticlastic/Poisson effect that "
                "1-D beam theory omits, and is a modelling difference rather than a "
                "code error. Self-convergence against the finest mesh isolates the "
                "discretisation error."));
  summary.set("mesh_convergence", block);

  StudyOutcome outcome;
  outcome.name = "mesh convergence (cantilever tip deflection)";
  outcome.kind = "verification + validation";
  outcome.metric = "relative error vs Timoshenko on the finest mesh (nu = 0)";
  outcome.value = finest_error_by_poisson[0.0];
  outcome.tolerance = 0.02;
  outcome.passed = outcome.value <= outcome.tolerance;
  std::ostringstream note;
  note << "observed order " << app::format(order_by_poisson[0.0], 3) << " (nu = 0), "
       << app::format(order_by_poisson[0.3], 3) << " (nu = 0.3)";
  outcome.note = note.str();
  return outcome;
}

/// Cantilever bending frequencies vs Euler-Bernoulli theory.
StudyOutcome study_modal(const std::string& out_dir, json::Value& summary) {
  CsvWriter csv(path_join(out_dir, "modal_convergence.csv"),
                {"nx", "ny", "num_dofs", "h[m]", "total_mass[kg]", "expected_mass[kg]",
                 "mass_error[-]", "f1_fem[Hz]", "f1_theory[Hz]", "f1_rel_error[-]",
                 "f2_fem[Hz]", "f2_theory[Hz]", "f2_rel_error[-]",
                 "f_axial_fem[Hz]", "f_axial_theory[Hz]", "f_axial_rel_error[-]",
                 "max_eigenpair_residual[-]", "seconds[s]"});

  CantileverSpec spec;
  spec.poisson = 0.0;  // removes the Poisson coupling beam theory omits
  const std::vector<std::pair<Index, Index>> meshes = {
      {20, 4}, {40, 8}, {80, 16}, {160, 32}};

  json::Value records = json::Value::make_array();
  Scalar finest_error = 0.0;
  Scalar finest_axial_error = 0.0;
  Scalar worst_mass_error = 0.0;
  for (const auto& m : meshes) {
    FemModel model = build_cantilever(spec, m.first, m.second);
    Assembler assembler(model);
    ModalAnalysisOptions options;
    options.num_modes = 8;  // enough to contain the first axial mode
    options.mass_type = MassType::Consistent;

    Timer timer;
    const ModalResult modal = solve_modal(model, assembler, options);
    const Scalar seconds = timer.elapsed_seconds();

    const Scalar f1_theory = cantilever_bending_frequency(
        spec.youngs, spec.density, spec.length, spec.height, spec.thickness, 1);
    const Scalar f2_theory = cantilever_bending_frequency(
        spec.youngs, spec.density, spec.length, spec.height, spec.thickness, 2);
    const Scalar expected_mass =
        spec.density * spec.length * spec.height * spec.thickness;
    const Scalar mass_error =
        std::abs(modal.total_mass - expected_mass) / expected_mass;
    worst_mass_error = std::max(worst_mass_error, mass_error);

    const Scalar f1 = modal.frequencies_hz(0);
    // The second *bending* mode: in a 2-D model the axial mode can appear
    // between the bending modes, so the comparison uses the next frequency
    // whose value is closest to the theoretical second bending frequency.
    Scalar f2 = 0.0;
    Scalar best = std::numeric_limits<Scalar>::max();
    for (Eigen::Index i = 1; i < modal.frequencies_hz.size(); ++i) {
      const Scalar d = std::abs(modal.frequencies_hz(i) - f2_theory);
      if (d < best) {
        best = d;
        f2 = modal.frequencies_hz(i);
      }
    }
    // The first axial mode: identified as the mode whose kinetic energy is
    // overwhelmingly in the x direction, then compared with fixed-free rod
    // theory. This is an independent analytical reference from the bending
    // series and is exact at nu = 0.
    const Scalar f_axial_theory =
        rod_axial_frequency(spec.youngs, spec.density, spec.length, 1);
    Scalar f_axial = 0.0;
    Scalar best_axial_ratio = 0.0;
    for (Eigen::Index mode = 0; mode < modal.mode_shapes.cols(); ++mode) {
      Scalar ex = 0.0;
      Scalar ey = 0.0;
      for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
        const Scalar ux = modal.mode_shapes(n * 2 + 0, mode);
        const Scalar uy = modal.mode_shapes(n * 2 + 1, mode);
        ex += ux * ux;
        ey += uy * uy;
      }
      const Scalar ratio = ex / std::max(ex + ey, 1.0e-300);
      if (ratio > best_axial_ratio) {
        best_axial_ratio = ratio;
        f_axial = modal.frequencies_hz(mode);
      }
    }
    const Scalar axial_error =
        std::abs(f_axial - f_axial_theory) / f_axial_theory;

    csv.row({static_cast<Scalar>(m.first), static_cast<Scalar>(m.second),
             static_cast<Scalar>(model.dofs().num_dofs()),
             spec.length / static_cast<Scalar>(m.first), modal.total_mass,
             expected_mass, mass_error, f1, f1_theory,
             std::abs(f1 - f1_theory) / f1_theory, f2, f2_theory,
             std::abs(f2 - f2_theory) / f2_theory, f_axial, f_axial_theory,
             axial_error, modal.modal_residuals.maxCoeff(), seconds});

    json::Value rec = json::Value::make_object();
    rec.set("nx", json::Value::make_number(m.first));
    rec.set("ny", json::Value::make_number(m.second));
    rec.set("f1_hz", json::Value::make_number(f1));
    rec.set("f1_theory_hz", json::Value::make_number(f1_theory));
    rec.set("f1_relative_error",
            json::Value::make_number(std::abs(f1 - f1_theory) / f1_theory));
    rec.set("f2_hz", json::Value::make_number(f2));
    rec.set("f2_theory_hz", json::Value::make_number(f2_theory));
    rec.set("f_axial_hz", json::Value::make_number(f_axial));
    rec.set("f_axial_theory_hz", json::Value::make_number(f_axial_theory));
    rec.set("f_axial_relative_error", json::Value::make_number(axial_error));
    rec.set("mass_relative_error", json::Value::make_number(mass_error));
    rec.set("max_eigenpair_residual",
            json::Value::make_number(modal.modal_residuals.maxCoeff()));
    records.push_back(rec);
    finest_error = std::abs(f1 - f1_theory) / f1_theory;
    finest_axial_error = axial_error;
  }
  csv.close();

  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("validation (Euler-Bernoulli frequencies) "
                                             "+ verification (mass conservation)"));
  block.set("geometry",
            json::Value::make_string("L = 1 m, h = 0.1 m, t = 0.01 m, E = 70 GPa, "
                                     "rho = 2700 kg/m^3, nu = 0"));
  block.set("bending_reference",
            json::Value::make_string(
                "f_n = (beta_n L)^2 / (2 pi L^2) sqrt(E h^2 / (12 rho)), "
                "beta_1 L = 1.87510407, beta_2 L = 4.69409113"));
  block.set("axial_reference",
            json::Value::make_string(
                "fixed-free rod, f_n = (2n-1)/(4L) sqrt(E/rho); exact for a bar "
                "in uniaxial stress, which plane stress reproduces at nu = 0"));
  block.set("records", records);
  block.set("max_mass_relative_error", json::Value::make_number(worst_mass_error));
  block.set("note",
            json::Value::make_string(
                "The 2-D model includes shear deformation and rotary inertia that "
                "Euler-Bernoulli theory omits, so the computed frequencies are "
                "expected to sit slightly below the theoretical values, and the gap "
                "grows with mode number. Conservation of total mass "
                "(sum(M)/2 = rho V) is a verification check and holds to round-off."));
  summary.set("modal_validation", block);

  StudyOutcome outcome;
  outcome.name = "modal frequencies (cantilever vs Euler-Bernoulli)";
  outcome.kind = "validation";
  outcome.metric = "f1 relative error on the finest mesh";
  outcome.value = finest_error;
  outcome.tolerance = 0.02;
  outcome.passed = finest_error <= outcome.tolerance;
  std::ostringstream note;
  note << "mass conservation error <= " << app::format(worst_mass_error, 3)
       << ", axial-mode error " << app::format(finest_axial_error, 3)
       << " vs fixed-free rod theory";
  outcome.note = note.str();
  return outcome;
}

/// All linear solvers on one small model.
StudyOutcome study_solver_agreement(const std::string& out_dir, json::Value& summary) {
  CantileverSpec spec;
  FemModel model = build_cantilever(spec, 12, 4);
  Assembler assembler(model);

  const std::vector<LinearSolverType> types = {
      LinearSolverType::SimplicialLdlt, LinearSolverType::SimplicialLlt,
      LinearSolverType::SparseLu,       LinearSolverType::ConjugateGradient,
      LinearSolverType::AmgCg,          LinearSolverType::Auto,
      LinearSolverType::DenseLu};

  CsvWriter csv(path_join(out_dir, "solver_agreement.csv"),
                {"solver", "num_dofs", "max_displacement[m]", "compliance[J]",
                 "scaled_residual[-]", "relative_difference_vs_dense[-]",
                 "iterations", "seconds[s]"});

  Vector reference;
  std::map<std::string, Scalar> differences;
  std::vector<std::vector<std::string>> rows;
  json::Value records = json::Value::make_array();

  // Solve with the dense solver first so it can serve as the reference.
  std::vector<LinearSolverType> ordered = {LinearSolverType::DenseLu};
  for (LinearSolverType t : types) {
    if (t != LinearSolverType::DenseLu) ordered.push_back(t);
  }

  for (LinearSolverType type : ordered) {
    StaticAnalysisOptions options;
    options.linear.type = type;
    options.linear.residual_tolerance = 1.0e-7;
    options.linear.iterative_tolerance = 1.0e-14;
    // A small coarse size gives the 130-unknown model a real three-level
    // hierarchy instead of a single direct solve.
    options.linear.amg.coarse_size = 16;
    Timer timer;
    StaticAnalysis analysis(model, assembler, options);
    const std::vector<StaticSolution> solutions = analysis.solve_all();
    const Scalar seconds = timer.elapsed_seconds();
    const StaticSolution& sol = solutions.front();

    if (reference.size() == 0) reference = sol.displacement;
    const Scalar diff = (sol.displacement - reference).cwiseAbs().maxCoeff() /
                        std::max(reference.cwiseAbs().maxCoeff(), 1.0e-300);
    differences[to_string(type)] = diff;

    rows.push_back({to_string(type),
                    app::format(static_cast<Scalar>(model.dofs().num_dofs()), 9),
                    app::format(sol.max_displacement_magnitude, 17),
                    app::format(sol.compliance, 17),
                    app::format(sol.scaled_residual, 6), app::format(diff, 6),
                    app::format(static_cast<Scalar>(sol.solver_iterations), 6),
                    app::format(seconds, 4)});

    json::Value rec = json::Value::make_object();
    rec.set("solver", json::Value::make_string(to_string(type)));
    rec.set("max_displacement_m",
            json::Value::make_number(sol.max_displacement_magnitude));
    rec.set("compliance_J", json::Value::make_number(sol.compliance));
    rec.set("scaled_residual", json::Value::make_number(sol.scaled_residual));
    rec.set("relative_difference_vs_dense", json::Value::make_number(diff));
    rec.set("iterations", json::Value::make_number(sol.solver_iterations));
    records.push_back(rec);
  }
  for (const auto& row : rows) csv.raw_row(row);
  csv.close();

  Scalar worst = 0.0;
  for (const auto& kv : differences) worst = std::max(worst, kv.second);

  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("verification"));
  block.set("mesh", json::Value::make_string("12 x 4 Q4, 130 DOFs"));
  block.set("reference_solver", json::Value::make_string("dense_lu"));
  block.set("records", records);
  block.set("max_relative_difference", json::Value::make_number(worst));
  summary.set("solver_agreement", block);

  StudyOutcome outcome;
  outcome.name = "solver agreement (sparse vs dense)";
  outcome.kind = "verification";
  outcome.metric = "max relative displacement difference vs dense LU";
  outcome.value = worst;
  outcome.tolerance = 1.0e-8;
  outcome.passed = worst <= outcome.tolerance;
  return outcome;
}

/// Constant-strain patch test on a distorted mesh.
StudyOutcome study_patch_test(const std::string& out_dir, json::Value& summary) {
  StructuredMeshSpec spec;
  spec.nx = 4;
  spec.ny = 4;
  spec.lx = 2.0;
  spec.ly = 1.0;

  // The exact linear displacement field: u = a + G x with a constant strain.
  const Vector2 offset(1.0e-4, -2.0e-4);
  Matrix2 gradient;
  gradient << 3.0e-4, 1.0e-4, 1.0e-4, -2.0e-4;
  const Vector3 exact_strain(gradient(0, 0), gradient(1, 1),
                             gradient(0, 1) + gradient(1, 0));

  CsvWriter csv(path_join(out_dir, "patch_test.csv"),
                {"perturbation[-]", "num_elements", "max_displacement_error[m]",
                 "relative_displacement_error[-]", "max_strain_error[-]",
                 "relative_strain_error[-]", "max_stress_error[Pa]",
                 "relative_stress_error[-]"});

  json::Value records = json::Value::make_array();
  Scalar worst_relative = 0.0;
  for (Scalar perturbation : {0.0, 0.15, 0.30}) {
    Mesh mesh = make_perturbed_quad_mesh(spec, perturbation, 20240917u);
    IsotropicMaterial material(200.0e9, 0.3, 7850.0, "patch_steel");
    FemModel model(std::move(mesh), material, 0.02, StressState::PlaneStress,
                   IntegrationOptions());

    // Prescribe the exact field on every boundary node.
    const std::vector<Mesh::BoundaryFace> edges = model.mesh().boundary_faces();
    std::vector<char> on_boundary(
        static_cast<std::size_t>(model.mesh().num_nodes()), 0);
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

    // A zero load case: the patch test is driven purely by the prescribed
    // boundary displacements.
    // The patch test is driven entirely by the prescribed boundary field.
    LoadCaseSpec load;
    load.name = "patch";
    load.prescribed_displacement_only = true;
    model.load_case_specs().push_back(load);
    model.finalize();

    // Overwrite the prescribed values with the exact field.
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
    const std::vector<StaticSolution> solutions = analysis.solve_all();
    const Vector& u = solutions.front().displacement;

    Scalar max_u_error = 0.0;
    Scalar u_scale = 0.0;
    for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
      const Vector2 x = model.mesh().node(n).head<2>();
      const Vector2 expected = offset + gradient * x;
      max_u_error = std::max(
          max_u_error, std::max(std::abs(u(n * 2 + 0) - expected.x()),
                                std::abs(u(n * 2 + 1) - expected.y())));
      u_scale = std::max(u_scale, expected.cwiseAbs().maxCoeff());
    }

    const StressField field = recover_stresses(model, assembler, u);
    const Vector3 exact_stress = material.plane_stress_matrix() * exact_strain;
    Scalar max_strain_error = 0.0;
    Scalar max_stress_error = 0.0;
    for (Index e = 0; e < model.mesh().num_elements(); ++e) {
      const Vector3 strain_e = field.element_strain.col(e);
      const Vector3 stress_e = field.element_stress.col(e);
      max_strain_error =
          std::max(max_strain_error, (strain_e - exact_strain).cwiseAbs().maxCoeff());
      max_stress_error =
          std::max(max_stress_error, (stress_e - exact_stress).cwiseAbs().maxCoeff());
    }

    const Scalar rel_u = max_u_error / std::max(u_scale, 1.0e-300);
    const Scalar rel_e =
        max_strain_error / std::max(exact_strain.cwiseAbs().maxCoeff(), 1.0e-300);
    const Scalar rel_s =
        max_stress_error / std::max(exact_stress.cwiseAbs().maxCoeff(), 1.0e-300);
    worst_relative = std::max({worst_relative, rel_u, rel_e, rel_s});

    csv.row({perturbation, static_cast<Scalar>(model.mesh().num_elements()),
             max_u_error, rel_u, max_strain_error, rel_e, max_stress_error, rel_s});

    json::Value rec = json::Value::make_object();
    rec.set("perturbation", json::Value::make_number(perturbation));
    rec.set("relative_displacement_error", json::Value::make_number(rel_u));
    rec.set("relative_strain_error", json::Value::make_number(rel_e));
    rec.set("relative_stress_error", json::Value::make_number(rel_s));
    records.push_back(rec);
  }
  csv.close();

  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("verification"));
  block.set("field", json::Value::make_string(
                         "u = (1e-4, -2e-4) + [[3e-4, 1e-4], [1e-4, -2e-4]] x, "
                         "constant strain (3e-4, -2e-4, 2e-4)"));
  block.set("records", records);
  block.set("max_relative_error", json::Value::make_number(worst_relative));
  block.set("note", json::Value::make_string(
                        "The interior displacement, the recovered strain and the "
                        "recovered stress must all reproduce the exact linear field "
                        "to round-off, including on distorted meshes."));
  summary.set("patch_test", block);

  StudyOutcome outcome;
  outcome.name = "patch test (constant strain, distorted mesh)";
  outcome.kind = "verification";
  outcome.metric = "max relative error in u, strain and stress";
  outcome.value = worst_relative;
  outcome.tolerance = 1.0e-10;
  outcome.passed = worst_relative <= outcome.tolerance;
  return outcome;
}


// ---------------------------------------------------------------------------
// Three-dimensional studies (Hex8)
// ---------------------------------------------------------------------------

/// Solid cantilever: length along x, height along y, width along z.
struct SolidCantileverSpec {
  Scalar length = 1.0;
  Scalar height = 0.1;
  Scalar width = 0.05;
  Scalar youngs = 70.0e9;
  Scalar poisson = 0.0;
  Scalar density = 2700.0;
  Scalar tip_load = -1000.0;
};

Scalar euler_bernoulli_tip_3d(const SolidCantileverSpec& s) {
  const Scalar i = s.width * s.height * s.height * s.height / 12.0;
  return s.tip_load * s.length * s.length * s.length / (3.0 * s.youngs * i);
}

Scalar timoshenko_tip_3d(const SolidCantileverSpec& s) {
  const Scalar i = s.width * s.height * s.height * s.height / 12.0;
  const Scalar a = s.width * s.height;
  const Scalar g = s.youngs / (2.0 * (1.0 + s.poisson));
  return s.tip_load * s.length * s.length * s.length / (3.0 * s.youngs * i) +
         s.tip_load * s.length / ((5.0 / 6.0) * g * a);
}

FemModel build_cantilever_3d(const SolidCantileverSpec& spec, Index nx, Index ny, Index nz) {
  StructuredMeshSpec mesh_spec;
  mesh_spec.nx = nx;
  mesh_spec.ny = ny;
  mesh_spec.nz = nz;
  mesh_spec.lx = spec.length;
  mesh_spec.ly = spec.height;
  mesh_spec.lz = spec.width;

  IsotropicMaterial material(spec.youngs, spec.poisson, spec.density, "verification_3d");
  FemModel model(make_structured_hex_mesh(mesh_spec), material, 1.0,
                 StressState::ThreeDimensional, IntegrationOptions());

  DisplacementConstraint root;
  root.region.name = "root";
  Selector box;
  box.kind = SelectorKind::Box;
  box.xmax = 0.0;
  root.region.members.push_back(box);
  root.fix_x = root.fix_y = root.fix_z = true;
  model.constraints().push_back(root);

  LoadCaseSpec load;
  load.name = "tip_load";
  PointLoadSpec tip;
  tip.region.name = "tip_face";
  Selector tip_box;
  tip_box.kind = SelectorKind::Box;
  tip_box.xmin = spec.length;
  tip.region.members.push_back(tip_box);
  tip.force = Vector3(0.0, spec.tip_load, 0.0);
  tip.distribute_total = true;
  load.point_loads.push_back(tip);
  model.load_case_specs().push_back(load);
  model.finalize();
  return model;
}

/// Mean u_y over the tip face [m].
Scalar tip_deflection_3d(const FemModel& model, const Vector& u) {
  const StructuredGridInfo& info = *model.mesh().structured_info();
  Scalar sum = 0.0;
  Index count = 0;
  for (Index k = 0; k <= info.nz; ++k) {
    for (Index j = 0; j <= info.ny; ++j) {
      sum += u(structured_node_index(info, info.nx, j, k) * 3 + 1);
      ++count;
    }
  }
  return sum / static_cast<Scalar>(count);
}

/// Constant-strain patch test on distorted Hex8 meshes.
StudyOutcome study_patch_test_3d(const std::string& out_dir, json::Value& summary) {
  StructuredMeshSpec spec;
  spec.nx = 3;
  spec.ny = 3;
  spec.nz = 3;
  spec.lx = 1.5;
  spec.ly = 1.0;
  spec.lz = 1.2;

  const Vector3 offset(1.0e-4, -2.0e-4, 0.5e-4);
  Matrix3 gradient;
  gradient << 3.0e-4, 1.0e-4, -0.5e-4,
              1.0e-4, -2.0e-4, 0.7e-4,
              -0.5e-4, 0.7e-4, 1.5e-4;
  Vector6 exact_strain;
  exact_strain << gradient(0, 0), gradient(1, 1), gradient(2, 2),
      gradient(0, 1) + gradient(1, 0), gradient(1, 2) + gradient(2, 1),
      gradient(2, 0) + gradient(0, 2);

  CsvWriter csv(path_join(out_dir, "patch_test_3d.csv"),
                {"perturbation[-]", "num_elements", "relative_displacement_error[-]",
                 "relative_strain_error[-]", "relative_stress_error[-]"});
  json::Value records = json::Value::make_array();
  Scalar worst_relative = 0.0;
  for (Scalar perturbation : {0.0, 0.15, 0.30}) {
    IsotropicMaterial material(200.0e9, 0.3, 7850.0, "patch_steel");
    FemModel model(make_perturbed_hex_mesh(spec, perturbation, 20240917u), material, 1.0,
                   StressState::ThreeDimensional, IntegrationOptions());

    std::vector<char> on_boundary(static_cast<std::size_t>(model.mesh().num_nodes()), 0);
    for (const Mesh::BoundaryFace& f : model.mesh().boundary_faces()) {
      for (Index n : f.nodes) on_boundary[static_cast<std::size_t>(n)] = 1;
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
    bc.fix_x = bc.fix_y = bc.fix_z = true;
    model.constraints().push_back(bc);
    LoadCaseSpec load;
    load.name = "patch";
    load.prescribed_displacement_only = true;
    model.load_case_specs().push_back(load);
    model.finalize();
    for (Index n : boundary_nodes) {
      const Vector3 u = offset + gradient * model.mesh().node(n);
      for (int k = 0; k < 3; ++k) model.dofs().prescribe(n, k, u(k));
    }

    Assembler assembler(model);
    StaticAnalysisOptions options;
    options.linear.residual_tolerance = 1.0e-9;
    StaticAnalysis analysis(model, assembler, options);
    const std::vector<StaticSolution> solutions = analysis.solve_all();
    const Vector& u = solutions.front().displacement;

    Scalar max_u_error = 0.0;
    Scalar u_scale = 0.0;
    for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
      const Vector3 expected = offset + gradient * model.mesh().node(n);
      for (int k = 0; k < 3; ++k) {
        max_u_error = std::max(max_u_error, std::abs(u(n * 3 + k) - expected(k)));
      }
      u_scale = std::max(u_scale, expected.cwiseAbs().maxCoeff());
    }
    const StressField field = recover_stresses(model, assembler, u);
    const Vector6 exact_stress = material.three_dimensional_matrix() * exact_strain;
    Scalar max_strain_error = 0.0;
    Scalar max_stress_error = 0.0;
    for (Index e = 0; e < model.mesh().num_elements(); ++e) {
      const Vector6 strain_e = field.element_strain.col(e);
      const Vector6 stress_e = field.element_stress.col(e);
      max_strain_error = std::max(max_strain_error,
                                  (strain_e - exact_strain).cwiseAbs().maxCoeff());
      max_stress_error = std::max(max_stress_error,
                                  (stress_e - exact_stress).cwiseAbs().maxCoeff());
    }
    const Scalar rel_u = max_u_error / std::max(u_scale, 1.0e-300);
    const Scalar rel_e = max_strain_error / exact_strain.cwiseAbs().maxCoeff();
    const Scalar rel_s = max_stress_error / exact_stress.cwiseAbs().maxCoeff();
    worst_relative = std::max({worst_relative, rel_u, rel_e, rel_s});
    csv.row({perturbation, static_cast<Scalar>(model.mesh().num_elements()), rel_u, rel_e,
             rel_s});
    json::Value rec = json::Value::make_object();
    rec.set("perturbation", json::Value::make_number(perturbation));
    rec.set("relative_displacement_error", json::Value::make_number(rel_u));
    rec.set("relative_strain_error", json::Value::make_number(rel_e));
    rec.set("relative_stress_error", json::Value::make_number(rel_s));
    records.push_back(rec);
  }
  csv.close();

  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("verification"));
  block.set("element", json::Value::make_string("Hex8, 3 x 3 x 3 cells"));
  block.set("field", json::Value::make_string(
                         "u = (1e-4, -2e-4, 0.5e-4) + G x with symmetric G, six constant "
                         "strain components"));
  block.set("records", records);
  block.set("max_relative_error", json::Value::make_number(worst_relative));
  summary.set("patch_test_3d", block);

  StudyOutcome outcome;
  outcome.name = "patch test 3-D (constant strain, distorted Hex8 mesh)";
  outcome.kind = "verification";
  outcome.metric = "max relative error in u, strain and stress";
  outcome.value = worst_relative;
  outcome.tolerance = 1.0e-10;
  outcome.passed = worst_relative <= outcome.tolerance;
  return outcome;
}

/// Solid cantilever tip deflection vs mesh size against Timoshenko theory.
StudyOutcome study_mesh_convergence_3d(const std::string& out_dir, json::Value& summary) {
  CsvWriter csv(path_join(out_dir, "mesh_convergence_3d.csv"),
                {"nx", "ny", "nz", "num_elements", "num_dofs", "h[m]", "tip_mean[m]",
                 "euler_bernoulli[m]", "timoshenko[m]", "rel_error_eb[-]",
                 "rel_error_timoshenko[-]", "compliance[J]", "strain_energy[J]",
                 "max_von_mises[Pa]", "reaction_rel_error[-]", "solve_seconds[s]"});
  SolidCantileverSpec spec;
  const Scalar eb = euler_bernoulli_tip_3d(spec);
  const Scalar ti = timoshenko_tip_3d(spec);
  const std::vector<std::array<Index, 3>> meshes = {
      {16, 2, 1}, {32, 4, 2}, {64, 8, 4}, {96, 12, 6}};

  json::Value records = json::Value::make_array();
  std::vector<Scalar> tips;
  std::vector<Scalar> sizes;
  for (const auto& m : meshes) {
    FemModel model = build_cantilever_3d(spec, m[0], m[1], m[2]);
    Assembler assembler(model);
    Timer timer;
    StaticAnalysis analysis(model, assembler, StaticAnalysisOptions());
    const StaticSolution sol = analysis.solve_all().front();
    const Scalar seconds = timer.elapsed_seconds();
    const StressField field = recover_stresses(model, assembler, sol.displacement);
    const Scalar tip = tip_deflection_3d(model, sol.displacement);
    const Scalar h = spec.length / static_cast<Scalar>(m[0]);
    tips.push_back(tip);
    sizes.push_back(h);
    csv.row({static_cast<Scalar>(m[0]), static_cast<Scalar>(m[1]),
             static_cast<Scalar>(m[2]), static_cast<Scalar>(model.mesh().num_elements()),
             static_cast<Scalar>(model.dofs().num_dofs()), h, tip, eb, ti,
             std::abs(tip - eb) / std::abs(eb), std::abs(tip - ti) / std::abs(ti),
             sol.compliance, sol.strain_energy, field.element_von_mises.maxCoeff(),
             sol.equilibrium.relative_force_error, seconds});
    json::Value rec = json::Value::make_object();
    rec.set("nx", json::Value::make_number(m[0]));
    rec.set("ny", json::Value::make_number(m[1]));
    rec.set("nz", json::Value::make_number(m[2]));
    rec.set("num_dofs", json::Value::make_number(model.dofs().num_dofs()));
    rec.set("h_m", json::Value::make_number(h));
    rec.set("tip_mean_m", json::Value::make_number(tip));
    rec.set("compliance_J", json::Value::make_number(sol.compliance));
    rec.set("relative_error_vs_timoshenko",
            json::Value::make_number(std::abs(tip - ti) / std::abs(ti)));
    rec.set("solve_seconds", json::Value::make_number(seconds));
    records.push_back(rec);
  }
  csv.close();

  // Self-convergence order against the finest mesh.
  std::vector<Scalar> errors;
  for (std::size_t i = 0; i + 1 < tips.size(); ++i) {
    errors.push_back(std::abs(tips[i] - tips.back()) / std::abs(tips.back()));
  }
  Scalar order = 0.0;
  if (errors.size() >= 2) {
    order = observed_order(sizes[errors.size() - 2], errors[errors.size() - 2],
                           sizes[errors.size() - 1], errors[errors.size() - 1]);
  }
  const Scalar finest_error = std::abs(tips.back() - ti) / std::abs(ti);

  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("verification (self-convergence) + "
                                             "validation (beam theory)"));
  block.set("geometry", json::Value::make_string(
                            "L = 1 m, h = 0.1 m, b = 0.05 m, E = 70 GPa, nu = 0, "
                            "tip resultant -1000 N over the tip face"));
  block.set("records", records);
  block.set("observed_convergence_order_tip_deflection", json::Value::make_number(order));
  block.set("finest_mesh_relative_error_vs_timoshenko",
            json::Value::make_number(finest_error));
  block.set("note", json::Value::make_string(
                        "The fully integrated Hex8 is stiff in bending on coarse meshes "
                        "(shear locking), so the deflection approaches the Timoshenko "
                        "value from below as the section is refined. Self-convergence "
                        "against the finest mesh isolates the discretisation error."));
  summary.set("mesh_convergence_3d", block);

  StudyOutcome outcome;
  outcome.name = "mesh convergence 3-D (Hex8 cantilever tip deflection)";
  outcome.kind = "verification + validation";
  outcome.metric = "relative error vs Timoshenko on the finest mesh";
  outcome.value = finest_error;
  outcome.tolerance = 0.03;
  outcome.passed = finest_error <= outcome.tolerance;
  std::ostringstream note;
  note << "observed order " << app::format(order, 3) << " (self-convergence)";
  outcome.note = note.str();
  return outcome;
}

/// Finite-difference check of the compliance gradient on a hex mesh.
StudyOutcome study_sensitivity_3d(const std::string& out_dir, json::Value& summary,
                                  Scalar tolerance) {
  StructuredMeshSpec spec;
  spec.nx = 6;
  spec.ny = 3;
  spec.nz = 2;
  spec.lx = 0.6;
  spec.ly = 0.3;
  spec.lz = 0.2;
  FemModel model(make_structured_hex_mesh(spec), IsotropicMaterial(70.0e9, 0.3, 2700.0),
                 1.0, StressState::ThreeDimensional, IntegrationOptions());
  DisplacementConstraint root;
  root.region.name = "root";
  Selector box;
  box.kind = SelectorKind::Box;
  box.xmax = 0.0;
  root.region.members.push_back(box);
  root.fix_x = root.fix_y = root.fix_z = true;
  model.constraints().push_back(root);
  LoadCaseSpec load;
  load.name = "tip";
  PointLoadSpec tip;
  tip.region.name = "tip_edge";
  Selector tip_box;
  tip_box.kind = SelectorKind::Box;
  tip_box.xmin = spec.lx;
  tip_box.ymax = 0.0;
  tip.region.members.push_back(tip_box);
  tip.force = Vector3(0.0, -500.0, 100.0);
  load.point_loads.push_back(tip);
  model.load_case_specs().push_back(load);
  model.finalize();
  Assembler assembler(model);

  const Scalar cell = 0.1;
  DensityFilter filter(model.mesh(), FilterType::Density, 1.5 * cell);
  PassiveRegionSpec hole;
  hole.region.name = "passive_void_sphere";
  Selector sphere;
  sphere.kind = SelectorKind::Sphere;
  sphere.center = Vector3(0.35, 0.15, 0.10);
  sphere.radius = 0.06;
  hole.region.members.push_back(sphere);
  hole.solid = false;
  DesignDomain domain(model, 0.5, 0.5, {hole});

  SimpOptions simp;
  StaticAnalysisOptions analysis_options;
  analysis_options.linear.residual_tolerance = 1.0e-9;
  ComplianceObjective objective(model, assembler, filter, domain, simp, analysis_options);

  Vector x = domain.initial_design();
  for (Index e = 0; e < domain.num_elements(); ++e) {
    if (!domain.is_free(e)) continue;
    const Vector3 c = model.mesh().element_centroid(e);
    x(e) = 0.45 + 0.25 * std::sin(9.0 * c.x()) * std::cos(7.0 * c.y()) * std::cos(5.0 * c.z());
  }
  domain.clamp(x);

  const std::vector<Scalar> steps = {1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7};
  CsvWriter csv(path_join(out_dir, "sensitivity_steps_3d.csv"),
                {"step[-]", "num_tested", "num_excluded", "max_absolute_error[J]",
                 "max_relative_error[-]", "rms_relative_error[-]",
                 "directional_relative_error[-]", "passed[-]"});
  json::Value step_records = json::Value::make_array();
  Scalar best_error = std::numeric_limits<Scalar>::max();
  Scalar best_step = 0.0;
  Index excluded = 0;
  for (Scalar step : steps) {
    const SensitivityCheckResult check =
        verify_sensitivities(objective, domain, x, {}, step, tolerance);
    csv.row({step, static_cast<Scalar>(check.num_tested),
             static_cast<Scalar>(check.num_excluded), check.max_absolute_error,
             check.max_relative_error, check.rms_relative_error,
             check.directional_relative_error, check.passed ? 1.0 : 0.0});
    if (check.max_relative_error < best_error) {
      best_error = check.max_relative_error;
      best_step = step;
    }
    excluded = check.num_excluded;
    json::Value rec = json::Value::make_object();
    rec.set("step", json::Value::make_number(step));
    rec.set("num_tested", json::Value::make_number(check.num_tested));
    rec.set("num_excluded", json::Value::make_number(check.num_excluded));
    rec.set("max_relative_error", json::Value::make_number(check.max_relative_error));
    rec.set("directional_relative_error",
            json::Value::make_number(check.directional_relative_error));
    rec.set("passed", json::Value::make_bool(check.passed));
    step_records.push_back(rec);
  }
  csv.close();

  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("verification"));
  block.set("mesh", json::Value::make_string("6 x 3 x 2 Hex8, 36 elements"));
  block.set("filter", json::Value::make_string("density, radius = 1.5 cells"));
  block.set("tested_elements", json::Value::make_number(domain.num_elements()));
  block.set("excluded_elements", json::Value::make_number(excluded));
  block.set("tolerance", json::Value::make_number(tolerance));
  block.set("best_step", json::Value::make_number(best_step));
  block.set("best_max_relative_error", json::Value::make_number(best_error));
  block.set("steps", step_records);
  summary.set("sensitivity_3d", block);

  StudyOutcome outcome;
  outcome.name = "sensitivity 3-D (analytical vs central differences, Hex8)";
  outcome.kind = "verification";
  outcome.metric = "min over steps of max relative error";
  outcome.value = best_error;
  outcome.tolerance = tolerance;
  outcome.passed = best_error <= tolerance;
  std::ostringstream note;
  note << "best at h = " << best_step << "; " << excluded << " of "
       << domain.num_elements() << " candidates excluded at bounds";
  outcome.note = note.str();
  return outcome;
}

/// Solid cantilever frequencies about both bending axes vs beam theory.
StudyOutcome study_modal_3d(const std::string& out_dir, json::Value& summary) {
  CsvWriter csv(path_join(out_dir, "modal_convergence_3d.csv"),
                {"nx", "ny", "nz", "num_dofs", "total_mass[kg]", "expected_mass[kg]",
                 "mass_error[-]", "f1_weak_fem[Hz]", "f1_weak_theory[Hz]",
                 "f1_weak_rel_error[-]", "f1_strong_fem[Hz]", "f1_strong_theory[Hz]",
                 "f1_strong_rel_error[-]", "max_eigenpair_residual[-]", "seconds[s]"});
  SolidCantileverSpec spec;
  const std::vector<std::array<Index, 3>> meshes = {{30, 3, 2}, {60, 4, 4}, {90, 6, 6}};
  json::Value records = json::Value::make_array();
  Scalar finest_weak = 0.0;
  Scalar finest_strong = 0.0;
  Scalar worst_mass = 0.0;
  for (const auto& m : meshes) {
    FemModel model = build_cantilever_3d(spec, m[0], m[1], m[2]);
    Assembler assembler(model);
    ModalAnalysisOptions options;
    options.num_modes = 6;
    Timer timer;
    const ModalResult modal = solve_modal(model, assembler, options);
    const Scalar seconds = timer.elapsed_seconds();

    const Scalar expected_mass = spec.density * spec.length * spec.height * spec.width;
    const Scalar mass_error = std::abs(modal.total_mass - expected_mass) / expected_mass;
    worst_mass = std::max(worst_mass, mass_error);
    // Weak-axis bending (about z, in the width direction) is the first mode;
    // strong-axis bending (about y) is the closest computed frequency to its
    // theoretical value.
    const Scalar f_weak_theory = cantilever_bending_frequency(
        spec.youngs, spec.density, spec.length, spec.width, spec.height, 1);
    const Scalar f_strong_theory = cantilever_bending_frequency(
        spec.youngs, spec.density, spec.length, spec.height, spec.width, 1);
    const Scalar f_weak = modal.frequencies_hz(0);
    Scalar f_strong = 0.0;
    Scalar best = std::numeric_limits<Scalar>::max();
    for (Eigen::Index i = 1; i < modal.frequencies_hz.size(); ++i) {
      const Scalar d = std::abs(modal.frequencies_hz(i) - f_strong_theory);
      if (d < best) {
        best = d;
        f_strong = modal.frequencies_hz(i);
      }
    }
    const Scalar weak_error = std::abs(f_weak - f_weak_theory) / f_weak_theory;
    const Scalar strong_error = std::abs(f_strong - f_strong_theory) / f_strong_theory;
    csv.row({static_cast<Scalar>(m[0]), static_cast<Scalar>(m[1]),
             static_cast<Scalar>(m[2]), static_cast<Scalar>(model.dofs().num_dofs()),
             modal.total_mass, expected_mass, mass_error, f_weak, f_weak_theory,
             weak_error, f_strong, f_strong_theory, strong_error,
             modal.modal_residuals.maxCoeff(), seconds});
    json::Value rec = json::Value::make_object();
    rec.set("nx", json::Value::make_number(m[0]));
    rec.set("ny", json::Value::make_number(m[1]));
    rec.set("nz", json::Value::make_number(m[2]));
    rec.set("f1_weak_hz", json::Value::make_number(f_weak));
    rec.set("f1_weak_theory_hz", json::Value::make_number(f_weak_theory));
    rec.set("f1_weak_relative_error", json::Value::make_number(weak_error));
    rec.set("f1_strong_hz", json::Value::make_number(f_strong));
    rec.set("f1_strong_theory_hz", json::Value::make_number(f_strong_theory));
    rec.set("f1_strong_relative_error", json::Value::make_number(strong_error));
    rec.set("mass_relative_error", json::Value::make_number(mass_error));
    records.push_back(rec);
    finest_weak = weak_error;
    finest_strong = strong_error;
  }
  csv.close();

  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("validation (Euler-Bernoulli frequencies) "
                                             "+ verification (mass conservation)"));
  block.set("geometry", json::Value::make_string(
                            "L = 1 m, h = 0.1 m (y), b = 0.05 m (z), E = 70 GPa, "
                            "rho = 2700 kg/m^3, nu = 0"));
  block.set("records", records);
  block.set("max_mass_relative_error", json::Value::make_number(worst_mass));
  block.set("note", json::Value::make_string(
                        "The first mode bends about the weak axis (through the 0.05 m "
                        "width), which needs several Hex8 through the thickness to "
                        "shed its shear-locking stiffness; the strong-axis mode "
                        "converges faster. Mass conservation (sum(M)/3 = rho V) holds "
                        "to round-off."));
  summary.set("modal_validation_3d", block);

  StudyOutcome outcome;
  outcome.name = "modal frequencies 3-D (Hex8 cantilever vs Euler-Bernoulli)";
  outcome.kind = "validation";
  outcome.metric = "f1 (weak axis) relative error on the finest mesh";
  outcome.value = finest_weak;
  outcome.tolerance = 0.03;
  outcome.passed = finest_weak <= outcome.tolerance;
  std::ostringstream note;
  note << "strong-axis error " << app::format(finest_strong, 3)
       << ", mass conservation error <= " << app::format(worst_mass, 3);
  outcome.note = note.str();
  return outcome;
}

// ---------------------------------------------------------------------------
// Linear simplices, multigrid and projection
// ---------------------------------------------------------------------------

/// Largest relative errors of displacement, strain and stress when a linear
/// field is prescribed on the boundary of `mesh` (any element type).
struct PatchErrors {
  Scalar displacement = 0.0;
  Scalar strain = 0.0;
  Scalar stress = 0.0;
};

PatchErrors run_patch(Mesh mesh) {
  const int dim = mesh.dim();
  const IsotropicMaterial material(200.0e9, 0.3, 7850.0, "patch_steel");
  FemModel model(std::move(mesh), material, dim == 2 ? 0.02 : 1.0,
                 dim == 2 ? StressState::PlaneStress : StressState::ThreeDimensional,
                 IntegrationOptions());
  Matrix3 g;
  g << 3.0e-4, 1.0e-4, -0.5e-4,
       1.0e-4, -2.0e-4, 0.7e-4,
       -0.5e-4, 0.7e-4, 1.5e-4;
  const Vector3 offset(1.0e-4, -2.0e-4, 0.5e-4);
  if (dim == 2) {
    g.row(2).setZero();
    g.col(2).setZero();
  }
  std::vector<char> on_boundary(static_cast<std::size_t>(model.mesh().num_nodes()), 0);
  for (const Mesh::BoundaryFace& f : model.mesh().boundary_faces()) {
    for (Index n : f.nodes) on_boundary[static_cast<std::size_t>(n)] = 1;
  }
  std::vector<Index> boundary;
  for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
    if (on_boundary[static_cast<std::size_t>(n)]) boundary.push_back(n);
  }
  DisplacementConstraint bc;
  Selector sel;
  sel.kind = SelectorKind::NodeIds;
  sel.ids = boundary;
  bc.region.members.push_back(sel);
  bc.fix_x = bc.fix_y = true;
  bc.fix_z = dim == 3;
  model.constraints().push_back(bc);
  LoadCaseSpec load;
  load.name = "patch";
  load.prescribed_displacement_only = true;
  model.load_case_specs().push_back(load);
  model.finalize();
  for (Index n : boundary) {
    const Vector3 u = offset + g * model.mesh().node(n);
    for (int k = 0; k < dim; ++k) model.dofs().prescribe(n, k, u(k));
  }
  Assembler assembler(model);
  StaticAnalysisOptions options;
  options.linear.type = LinearSolverType::SimplicialLdlt;
  options.linear.residual_tolerance = 1.0e-9;
  StaticAnalysis analysis(model, assembler, options);
  const Vector u = analysis.solve_all().front().displacement;

  PatchErrors out;
  Scalar u_scale = 0.0;
  for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
    const Vector3 expected = offset + g * model.mesh().node(n);
    for (int k = 0; k < dim; ++k) {
      out.displacement = std::max(out.displacement, std::abs(u(n * dim + k) - expected(k)));
      u_scale = std::max(u_scale, std::abs(expected(k)));
    }
  }
  out.displacement /= u_scale;
  Vector strain;
  Vector stress;
  if (dim == 2) {
    const Vector3 e(g(0, 0), g(1, 1), g(0, 1) + g(1, 0));
    strain = e;
    stress = material.plane_stress_matrix() * e;
  } else {
    Vector6 e;
    e << g(0, 0), g(1, 1), g(2, 2), g(0, 1) + g(1, 0), g(1, 2) + g(2, 1), g(2, 0) + g(0, 2);
    strain = e;
    stress = material.three_dimensional_matrix() * e;
  }
  const StressField field = recover_stresses(model, assembler, u);
  for (Index e = 0; e < model.mesh().num_elements(); ++e) {
    out.strain = std::max(out.strain, (Vector(field.element_strain.col(e)) - strain)
                                          .cwiseAbs()
                                          .maxCoeff());
    out.stress = std::max(out.stress, (Vector(field.element_stress.col(e)) - stress)
                                          .cwiseAbs()
                                          .maxCoeff());
  }
  out.strain /= strain.cwiseAbs().maxCoeff();
  out.stress /= stress.cwiseAbs().maxCoeff();
  return out;
}

/// Constant-strain patch test on distorted Tri3 and Tet4 meshes.
StudyOutcome study_patch_test_simplex(const std::string& out_dir, json::Value& summary) {
  CsvWriter csv(path_join(out_dir, "patch_test_simplex.csv"),
                {"element", "perturbation[-]", "num_elements",
                 "relative_displacement_error[-]", "relative_strain_error[-]",
                 "relative_stress_error[-]"});
  json::Value records = json::Value::make_array();
  Scalar worst = 0.0;
  for (const ElementType type : {ElementType::Tri3, ElementType::Tet4}) {
    const std::vector<Scalar> perturbations =
        type == ElementType::Tri3 ? std::vector<Scalar>{0.0, 0.2, 0.4}
                                  : std::vector<Scalar>{0.0, 0.15, 0.3};
    for (Scalar perturbation : perturbations) {
      StructuredMeshSpec spec;
      spec.nx = type == ElementType::Tri3 ? 5 : 3;
      spec.ny = type == ElementType::Tri3 ? 4 : 3;
      spec.nz = 3;
      spec.lx = 1.5;
      spec.ly = 1.0;
      spec.lz = 1.2;
      Mesh mesh = type == ElementType::Tri3
                      ? make_perturbed_tri_mesh(spec, perturbation, 20240917u)
                      : make_perturbed_tet_mesh(spec, perturbation, 20240917u);
      const Index ne = mesh.num_elements();
      const PatchErrors e = run_patch(std::move(mesh));
      worst = std::max({worst, e.displacement, e.strain, e.stress});
      csv.raw_row({to_string(type), app::format(perturbation, 3),
                   app::format(static_cast<Scalar>(ne), 9), app::format(e.displacement, 4),
                   app::format(e.strain, 4), app::format(e.stress, 4)});
      json::Value rec = json::Value::make_object();
      rec.set("element", json::Value::make_string(to_string(type)));
      rec.set("perturbation", json::Value::make_number(perturbation));
      rec.set("num_elements", json::Value::make_number(ne));
      rec.set("relative_displacement_error", json::Value::make_number(e.displacement));
      rec.set("relative_strain_error", json::Value::make_number(e.strain));
      rec.set("relative_stress_error", json::Value::make_number(e.stress));
      records.push_back(rec);
    }
  }
  csv.close();
  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("verification"));
  block.set("records", records);
  block.set("max_relative_error", json::Value::make_number(worst));
  summary.set("patch_test_simplex", block);

  StudyOutcome outcome;
  outcome.name = "patch test (constant strain, distorted Tri3 and Tet4 meshes)";
  outcome.kind = "verification";
  outcome.metric = "max relative displacement / strain / stress error";
  outcome.value = worst;
  outcome.tolerance = 1.0e-10;
  outcome.passed = worst <= outcome.tolerance;
  return outcome;
}

/// Mean deflection over the nodes of the tip section x = L.
Scalar tip_mean_deflection(const FemModel& model, const Vector& u, Scalar length) {
  const int dim = model.dim();
  Scalar sum = 0.0;
  Index count = 0;
  for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
    if (model.mesh().node(n).x() >= length - 1.0e-12) {
      sum += u(n * dim + 1);
      ++count;
    }
  }
  return sum / static_cast<Scalar>(count);
}

/// Linear-triangle and linear-tetrahedron cantilevers against beam theory.
StudyOutcome study_mesh_convergence_simplex(const std::string& out_dir, json::Value& summary) {
  CsvWriter csv(path_join(out_dir, "mesh_convergence_simplex.csv"),
                {"element", "nx", "ny", "nz", "num_elements", "num_dofs", "h[m]",
                 "tip_mean[m]", "timoshenko[m]", "rel_error_timoshenko[-]",
                 "compliance[J]", "reaction_rel_error[-]", "solver", "solve_seconds[s]"});
  json::Value blocks = json::Value::make_object();
  Scalar worst_finest = 0.0;
  std::ostringstream note;
  for (const ElementType type : {ElementType::Tri3, ElementType::Tet4}) {
    const bool solid = type == ElementType::Tet4;
    CantileverSpec plane;
    SolidCantileverSpec spec3;
    const Scalar reference = solid ? timoshenko_tip_3d(spec3) : timoshenko_tip(plane);
    const std::vector<std::array<Index, 3>> grids =
        solid ? std::vector<std::array<Index, 3>>{{16, 2, 1}, {32, 4, 2}, {64, 8, 4},
                                                  {96, 12, 6}, {128, 16, 8}}
              : std::vector<std::array<Index, 3>>{{20, 2, 0}, {40, 4, 0}, {80, 8, 0},
                                                  {160, 16, 0}, {320, 32, 0}};
    std::vector<Scalar> tips;
    std::vector<Scalar> sizes;
    json::Value records = json::Value::make_array();
    for (const auto& grid : grids) {
      StructuredMeshSpec ms;
      ms.nx = grid[0];
      ms.ny = grid[1];
      ms.nz = std::max<Index>(grid[2], 1);
      ms.lx = solid ? spec3.length : plane.length;
      ms.ly = solid ? spec3.height : plane.height;
      ms.lz = spec3.width;
      const IsotropicMaterial material =
          solid ? IsotropicMaterial(spec3.youngs, spec3.poisson, spec3.density, "tet")
                : IsotropicMaterial(plane.youngs, plane.poisson, plane.density, "tri");
      FemModel model(solid ? make_structured_tet_mesh(ms) : make_structured_tri_mesh(ms),
                     material, solid ? 1.0 : plane.thickness,
                     solid ? StressState::ThreeDimensional : StressState::PlaneStress,
                     IntegrationOptions());
      DisplacementConstraint root;
      Selector box;
      box.kind = SelectorKind::Box;
      box.xmax = 0.0;
      root.region.members.push_back(box);
      root.fix_x = root.fix_y = true;
      root.fix_z = solid;
      model.constraints().push_back(root);
      LoadCaseSpec load;
      load.name = "tip_load";
      PointLoadSpec tip;
      Selector tip_box;
      tip_box.kind = SelectorKind::Box;
      tip_box.xmin = ms.lx;
      tip.region.members.push_back(tip_box);
      tip.force = Vector3(0.0, solid ? spec3.tip_load : plane.tip_load, 0.0);
      tip.distribute_total = true;
      load.point_loads.push_back(tip);
      model.load_case_specs().push_back(load);
      model.finalize();
      Assembler assembler(model);
      Timer timer;
      StaticAnalysis analysis(model, assembler, StaticAnalysisOptions());
      const StaticSolution sol = analysis.solve_all().front();
      const Scalar seconds = timer.elapsed_seconds();
      const Scalar tip_mean = tip_mean_deflection(model, sol.displacement, ms.lx);
      const Scalar h = ms.lx / static_cast<Scalar>(ms.nx);
      tips.push_back(tip_mean);
      sizes.push_back(h);
      const Scalar error = std::abs(tip_mean - reference) / std::abs(reference);
      csv.raw_row({to_string(type), app::format(static_cast<Scalar>(grid[0]), 6),
                   app::format(static_cast<Scalar>(grid[1]), 6),
                   app::format(static_cast<Scalar>(grid[2]), 6),
                   app::format(static_cast<Scalar>(model.mesh().num_elements()), 9),
                   app::format(static_cast<Scalar>(model.dofs().num_dofs()), 9),
                   app::format(h, 6), app::format(tip_mean, 12), app::format(reference, 12),
                   app::format(error, 6), app::format(sol.compliance, 12),
                   app::format(sol.equilibrium.relative_force_error, 4), sol.solver_name,
                   app::format(seconds, 4)});
      json::Value rec = json::Value::make_object();
      rec.set("nx", json::Value::make_number(grid[0]));
      rec.set("num_dofs", json::Value::make_number(model.dofs().num_dofs()));
      rec.set("h_m", json::Value::make_number(h));
      rec.set("tip_mean_m", json::Value::make_number(tip_mean));
      rec.set("relative_error_vs_timoshenko", json::Value::make_number(error));
      rec.set("solver", json::Value::make_string(sol.solver_name));
      rec.set("solver_iterations", json::Value::make_number(sol.solver_iterations));
      records.push_back(rec);
    }
    std::vector<Scalar> errors;
    for (std::size_t i = 0; i + 1 < tips.size(); ++i) {
      errors.push_back(std::abs(tips[i] - tips.back()) / std::abs(tips.back()));
    }
    const std::size_t k = errors.size();
    const Scalar order = k >= 2 ? observed_order(sizes[k - 2], errors[k - 2], sizes[k - 1],
                                                 errors[k - 1])
                                : 0.0;
    const Scalar finest = std::abs(tips.back() - reference) / std::abs(reference);
    worst_finest = std::max(worst_finest, finest);
    json::Value block = json::Value::make_object();
    block.set("records", records);
    block.set("observed_convergence_order_tip_deflection", json::Value::make_number(order));
    block.set("finest_mesh_relative_error_vs_timoshenko", json::Value::make_number(finest));
    blocks.set(to_string(type), block);
    note << to_string(type) << ": finest error " << app::format(100.0 * finest, 3)
         << " %, observed order " << app::format(order, 3) << "; ";
  }
  csv.close();
  blocks.set("kind", json::Value::make_string("verification (self-convergence) + "
                                              "validation (Timoshenko beam theory)"));
  blocks.set("note", json::Value::make_string(
                         "The same cantilevers as the Q4 and Hex8 studies on meshes made by "
                         "splitting their cells (2 triangles per quadrilateral, 6 Kuhn "
                         "tetrahedra per hexahedron). Constant-strain simplices lock in "
                         "bending, so the deflection approaches the beam value from below "
                         "and more slowly than with the bilinear / trilinear elements."));
  summary.set("mesh_convergence_simplex", blocks);

  StudyOutcome outcome;
  outcome.name = "mesh convergence (Tri3 and Tet4 cantilever tip deflection)";
  outcome.kind = "verification + validation";
  outcome.metric = "larger finest-mesh relative error vs Timoshenko";
  outcome.value = worst_finest;
  outcome.tolerance = 0.03;
  outcome.passed = worst_finest <= outcome.tolerance;
  outcome.note = note.str();
  return outcome;
}

/// Multigrid-preconditioned CG against the Cholesky factorisation and against
/// Jacobi-preconditioned CG, as the mesh is refined.
StudyOutcome study_multigrid(const std::string& out_dir, json::Value& summary) {
  CsvWriter csv(path_join(out_dir, "multigrid_scaling.csv"),
                {"element", "nx", "num_dofs", "amg_iterations", "amg_levels",
                 "operator_complexity", "amg_seconds[s]", "jacobi_iterations",
                 "jacobi_seconds[s]", "ldlt_seconds[s]", "relative_difference_vs_ldlt[-]"});
  json::Value records = json::Value::make_array();
  Scalar worst_difference = 0.0;
  Scalar worst_growth = 0.0;
  std::ostringstream note;
  for (const ElementType type : {ElementType::Hex8, ElementType::Tet4}) {
    const std::vector<Index> sizes = type == ElementType::Hex8
                                         ? std::vector<Index>{8, 16, 32, 48}
                                         : std::vector<Index>{6, 12, 24, 36};
    std::vector<int> amg_counts;
    for (Index nx : sizes) {
      StructuredMeshSpec ms;
      ms.nx = nx;
      ms.ny = nx / 2;
      ms.nz = std::max<Index>(1, nx / 4);
      ms.lx = 2.0;
      ms.ly = 1.0;
      ms.lz = 0.5;
      FemModel model(type == ElementType::Hex8 ? make_structured_hex_mesh(ms)
                                               : make_structured_tet_mesh(ms),
                     IsotropicMaterial(70.0e9, 0.3, 2700.0, "mg"), 1.0,
                     StressState::ThreeDimensional, IntegrationOptions());
      DisplacementConstraint root;
      Selector box;
      box.kind = SelectorKind::Box;
      box.xmax = 0.0;
      root.region.members.push_back(box);
      root.fix_x = root.fix_y = root.fix_z = true;
      model.constraints().push_back(root);
      LoadCaseSpec load;
      load.name = "tip";
      PointLoadSpec p;
      Selector corner;
      corner.kind = SelectorKind::NearestNode;
      corner.point = Vector3(2.0, 0.0, 0.0);
      p.region.members.push_back(corner);
      p.force = Vector3(0.0, -1000.0, 0.0);
      load.point_loads.push_back(p);
      model.load_case_specs().push_back(load);
      model.finalize();
      Assembler assembler(model);
      const auto solve = [&](LinearSolverType t, Scalar* seconds, int* iterations,
                             int* levels, Scalar* complexity) {
        StaticAnalysisOptions o;
        o.linear.type = t;
        o.linear.iterative_tolerance = 1.0e-10;
        o.linear.max_iterations = 20000;
        Timer timer;
        StaticAnalysis analysis(model, assembler, o);
        const StaticSolution sol = analysis.solve_all().front();
        *seconds = timer.elapsed_seconds();
        *iterations = sol.solver_iterations;
        if (levels != nullptr) {
          const AmgStats* st = analysis.solver().amg_stats();
          *levels = st ? static_cast<int>(st->levels.size()) : 0;
          *complexity = st ? st->operator_complexity : 0.0;
        }
        return sol.displacement;
      };
      Scalar amg_seconds = 0.0;
      Scalar jacobi_seconds = 0.0;
      Scalar ldlt_seconds = 0.0;
      int amg_iterations = 0;
      int jacobi_iterations = 0;
      int ldlt_iterations = 0;
      int levels = 0;
      Scalar complexity = 0.0;
      const Vector u_amg =
          solve(LinearSolverType::AmgCg, &amg_seconds, &amg_iterations, &levels, &complexity);
      solve(LinearSolverType::ConjugateGradient, &jacobi_seconds, &jacobi_iterations, nullptr,
            nullptr);
      const Vector u_ldlt = solve(LinearSolverType::SimplicialLdlt, &ldlt_seconds,
                                  &ldlt_iterations, nullptr, nullptr);
      const Scalar difference = (u_amg - u_ldlt).norm() / u_ldlt.norm();
      worst_difference = std::max(worst_difference, difference);
      if (levels >= 2) amg_counts.push_back(amg_iterations);
      csv.raw_row({to_string(type), app::format(static_cast<Scalar>(nx), 6),
                   app::format(static_cast<Scalar>(model.dofs().num_dofs()), 9),
                   app::format(static_cast<Scalar>(amg_iterations), 6),
                   app::format(static_cast<Scalar>(levels), 3), app::format(complexity, 4),
                   app::format(amg_seconds, 4),
                   app::format(static_cast<Scalar>(jacobi_iterations), 6),
                   app::format(jacobi_seconds, 4), app::format(ldlt_seconds, 4),
                   app::format(difference, 4)});
      json::Value rec = json::Value::make_object();
      rec.set("element", json::Value::make_string(to_string(type)));
      rec.set("nx", json::Value::make_number(nx));
      rec.set("num_dofs", json::Value::make_number(model.dofs().num_dofs()));
      rec.set("amg_iterations", json::Value::make_number(amg_iterations));
      rec.set("amg_levels", json::Value::make_number(levels));
      rec.set("operator_complexity", json::Value::make_number(complexity));
      rec.set("jacobi_iterations", json::Value::make_number(jacobi_iterations));
      rec.set("amg_seconds", json::Value::make_number(amg_seconds));
      rec.set("jacobi_seconds", json::Value::make_number(jacobi_seconds));
      rec.set("ldlt_seconds", json::Value::make_number(ldlt_seconds));
      rec.set("relative_difference_vs_ldlt", json::Value::make_number(difference));
      records.push_back(rec);
    }
    // Over the sizes that form a hierarchy; smaller ones are below
    // coarse_size and solved directly in one iteration.
    const Scalar growth = amg_counts.empty()
                              ? 0.0
                              : static_cast<Scalar>(amg_counts.back()) /
                                    static_cast<Scalar>(std::max(amg_counts.front(), 1));
    worst_growth = std::max(worst_growth, growth);
    if (!amg_counts.empty()) {
      note << to_string(type) << ": AMG iterations " << amg_counts.front() << " -> "
           << amg_counts.back() << " over the multi-level sizes; ";
    }
  }
  csv.close();
  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("verification"));
  block.set("tolerance", json::Value::make_number(1.0e-10));
  block.set("records", records);
  block.set("max_relative_difference_vs_ldlt", json::Value::make_number(worst_difference));
  block.set("max_iteration_growth_coarsest_to_finest", json::Value::make_number(worst_growth));
  summary.set("multigrid", block);

  StudyOutcome outcome;
  outcome.name = "multigrid CG (agreement with Cholesky, iterations under refinement)";
  outcome.kind = "verification";
  outcome.metric = "max relative displacement difference vs LDL^T";
  outcome.value = worst_difference;
  outcome.tolerance = 1.0e-8;
  outcome.passed = worst_difference <= outcome.tolerance && worst_growth <= 1.6;
  note << "iteration growth coarsest -> finest " << app::format(worst_growth, 3)
       << " (limit 1.6)";
  outcome.note = note.str();
  return outcome;
}

/// Central-difference check of the compliance gradient through the Heaviside
/// projection, on a quadrilateral and a tetrahedral mesh.
StudyOutcome study_sensitivity_projection(const std::string& out_dir, json::Value& summary,
                                          Scalar tolerance) {
  CsvWriter csv(path_join(out_dir, "sensitivity_projection.csv"),
                {"element", "beta[-]", "eta[-]", "step[-]", "num_tested",
                 "max_relative_error[-]", "max_scaled_error[-]",
                 "directional_relative_error[-]"});
  json::Value records = json::Value::make_array();
  Scalar worst = 0.0;
  for (const ElementType type : {ElementType::Quad4, ElementType::Tet4}) {
    StructuredMeshSpec ms;
    const bool solid = type == ElementType::Tet4;
    ms.nx = solid ? 5 : 12;
    ms.ny = solid ? 2 : 6;
    ms.nz = 2;
    ms.lx = solid ? 0.5 : 0.6;
    ms.ly = solid ? 0.2 : 0.3;
    ms.lz = 0.2;
    FemModel model(solid ? make_structured_tet_mesh(ms) : make_structured_quad_mesh(ms),
                   IsotropicMaterial(70.0e9, 0.3, 2700.0, "proj"), solid ? 1.0 : 0.01,
                   solid ? StressState::ThreeDimensional : StressState::PlaneStress,
                   IntegrationOptions());
    DisplacementConstraint root;
    Selector box;
    box.kind = SelectorKind::Box;
    box.xmax = 0.0;
    root.region.members.push_back(box);
    root.fix_x = root.fix_y = true;
    root.fix_z = solid;
    model.constraints().push_back(root);
    LoadCaseSpec load;
    load.name = "tip";
    PointLoadSpec p;
    Selector corner;
    corner.kind = SelectorKind::NearestNode;
    corner.point = Vector3(ms.lx, 0.0, 0.0);
    p.region.members.push_back(corner);
    p.force = Vector3(0.0, -500.0, solid ? 100.0 : 0.0);
    load.point_loads.push_back(p);
    model.load_case_specs().push_back(load);
    model.finalize();
    Assembler assembler(model);
    const DensityFilter filter(model.mesh(), FilterType::Density,
                               1.5 * model.mesh().mean_element_size());
    DesignDomain domain(model, 0.5, 0.5, {});
    StaticAnalysisOptions options;
    options.linear.type = LinearSolverType::SimplicialLdlt;
    options.linear.residual_tolerance = 1.0e-9;
    ComplianceObjective objective(model, assembler, filter, domain, SimpOptions(), options);
    Vector x = domain.initial_design();
    for (Index e = 0; e < domain.num_elements(); ++e) {
      const Vector3 c = model.mesh().element_centroid(e);
      x(e) = 0.5 + 0.3 * std::sin(9.0 * c.x()) * std::cos(7.0 * c.y() + 3.0 * c.z());
    }
    domain.clamp(x);
    for (Scalar beta : {2.0, 8.0, 32.0}) {
      objective.set_projection(beta, 0.5);
      Scalar best = std::numeric_limits<Scalar>::max();
      for (Scalar step : {1.0e-4, 1.0e-5, 1.0e-6}) {
        const SensitivityCheckResult check =
            verify_sensitivities(objective, domain, x, {}, step, tolerance);
        best = std::min(best, std::max(check.max_scaled_error,
                                       check.directional_relative_error));
        csv.raw_row({to_string(type), app::format(beta, 4), "0.5", app::format(step, 3),
                     app::format(static_cast<Scalar>(check.num_tested), 6),
                     app::format(check.max_relative_error, 4),
                     app::format(check.max_scaled_error, 4),
                     app::format(check.directional_relative_error, 4)});
      }
      worst = std::max(worst, best);
      json::Value rec = json::Value::make_object();
      rec.set("element", json::Value::make_string(to_string(type)));
      rec.set("beta", json::Value::make_number(beta));
      rec.set("best_max_scaled_or_directional_error", json::Value::make_number(best));
      records.push_back(rec);
    }
  }
  csv.close();
  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("verification"));
  block.set("records", records);
  block.set("worst_best_error", json::Value::make_number(worst));
  block.set("note", json::Value::make_string(
                        "Each entry is judged against max(|analytical|, |FD|, 1e-3 "
                        "||gradient||_inf): at large beta the projection's derivative "
                        "vanishes away from the threshold, and those entries are compared "
                        "with the round-off floor of a central difference. The "
                        "directional derivative along the full gradient is checked as "
                        "well."));
  summary.set("sensitivity_projection", block);

  StudyOutcome outcome;
  outcome.name = "sensitivity through the Heaviside projection (Q4 and Tet4)";
  outcome.kind = "verification";
  outcome.metric =
      "worst over beta of the best max(scaled entry error, directional error)";
  outcome.value = worst;
  outcome.tolerance = tolerance;
  outcome.passed = worst <= tolerance;
  return outcome;
}

// ---------------------------------------------------------------------------
// Quadratic tetrahedra, buckling and manufacturing constraints
// ---------------------------------------------------------------------------

/// Clamped-root solid cantilever of `mesh` under a traction on the tip face
/// whose resultant is `force` (consistent nodal loads for every element).
FemModel solid_cantilever(Mesh mesh, const SolidCantileverSpec& spec, const Vector3& force) {
  FemModel model(std::move(mesh), IsotropicMaterial(spec.youngs, spec.poisson, spec.density, "c"),
                 1.0, StressState::ThreeDimensional, IntegrationOptions());
  DisplacementConstraint root;
  Selector box;
  box.kind = SelectorKind::Box;
  box.xmax = 0.0;
  root.region.members.push_back(box);
  root.fix_x = root.fix_y = root.fix_z = true;
  model.constraints().push_back(root);
  LoadCaseSpec load;
  load.name = "tip";
  TractionLoadSpec tip;
  Selector end;
  end.kind = SelectorKind::Box;
  end.xmin = spec.length;
  tip.region.members.push_back(end);
  tip.traction = force / (spec.height * spec.width);
  load.tractions.push_back(tip);
  model.load_case_specs().push_back(load);
  model.finalize();
  return model;
}

/// The pure-bending field u_x = -k x z, u_y = nu k y z,
/// u_z = k/2 (x^2 + nu (z^2 - y^2)): an exact solution with sigma_xx = -E k z.
Vector3 pure_bending(const Vector3& x, Scalar kappa, Scalar nu) {
  return Vector3(-kappa * x.x() * x.z(), nu * kappa * x.y() * x.z(),
                 0.5 * kappa * (x.x() * x.x() + nu * (x.z() * x.z() - x.y() * x.y())));
}

/// Quadratic (pure-bending) patch test: Tet10 must reproduce it exactly on
/// distorted meshes; Tet4 and Hex8 cannot, and their error is reported.
StudyOutcome study_patch_test_quadratic(const std::string& out_dir, json::Value& summary) {
  CsvWriter csv(path_join(out_dir, "patch_test_quadratic.csv"),
                {"element", "perturbation[-]", "num_elements", "num_dofs",
                 "relative_displacement_error[-]", "relative_stress_error[-]"});
  const Scalar kappa = 2.0e-3;
  const IsotropicMaterial material(200.0e9, 0.3, 7850.0, "patch");
  json::Value records = json::Value::make_array();
  Scalar worst_tet10 = 0.0;
  Scalar best_linear = std::numeric_limits<Scalar>::max();
  for (const ElementType type : {ElementType::Tet10, ElementType::Tet4, ElementType::Hex8}) {
    for (const Scalar perturbation : {0.0, 0.15, 0.3}) {
      StructuredMeshSpec ms;
      ms.nx = ms.ny = ms.nz = 3;
      ms.lx = 1.5;
      ms.ly = 1.0;
      ms.lz = 1.2;
      Mesh mesh = type == ElementType::Tet10  ? make_perturbed_tet10_mesh(ms, perturbation, 7u)
                  : type == ElementType::Tet4 ? make_perturbed_tet_mesh(ms, perturbation, 7u)
                                              : make_perturbed_hex_mesh(ms, perturbation, 7u);
      FemModel model(std::move(mesh), material, 1.0, StressState::ThreeDimensional,
                     IntegrationOptions());
      std::vector<char> on(static_cast<std::size_t>(model.mesh().num_nodes()), 0);
      for (const Mesh::BoundaryFace& f : model.mesh().boundary_faces()) {
        for (Index n : f.nodes) on[static_cast<std::size_t>(n)] = 1;
      }
      std::vector<Index> boundary;
      for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
        if (on[static_cast<std::size_t>(n)]) boundary.push_back(n);
      }
      DisplacementConstraint bc;
      Selector sel;
      sel.kind = SelectorKind::NodeIds;
      sel.ids = boundary;
      bc.region.members.push_back(sel);
      bc.fix_x = bc.fix_y = bc.fix_z = true;
      model.constraints().push_back(bc);
      LoadCaseSpec load;
      load.name = "bending";
      load.prescribed_displacement_only = true;
      model.load_case_specs().push_back(load);
      model.finalize();
      for (Index n : boundary) {
        const Vector3 u = pure_bending(model.mesh().node(n), kappa, 0.3);
        for (int k = 0; k < 3; ++k) model.dofs().prescribe(n, k, u(k));
      }
      Assembler assembler(model);
      StaticAnalysisOptions options;
      options.linear.type = LinearSolverType::SimplicialLdlt;
      StaticAnalysis analysis(model, assembler, options);
      const Vector u = analysis.solve_all().front().displacement;
      Scalar error = 0.0;
      Scalar scale = 0.0;
      for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
        const Vector3 exact = pure_bending(model.mesh().node(n), kappa, 0.3);
        error = std::max(error, (u.segment<3>(3 * n) - exact).cwiseAbs().maxCoeff());
        scale = std::max(scale, exact.cwiseAbs().maxCoeff());
      }
      // Element-average stress against the exact linear sigma_xx at the
      // element's volume centroid.
      const StressField field = recover_stresses(model, assembler, u);
      Scalar stress_error = 0.0;
      const Scalar stress_scale = material.youngs_modulus() * kappa * ms.lz;
      for (Index e = 0; e < model.mesh().num_elements(); ++e) {
        const Vector3 c = model.mesh().element_coordinates(e).leftCols(corner_nodes(type))
                              .rowwise().mean();
        Vector6 exact = Vector6::Zero();
        exact(0) = -material.youngs_modulus() * kappa * c.z();
        stress_error = std::max(stress_error,
                                (Vector6(field.element_stress.col(e)) - exact).cwiseAbs().maxCoeff());
      }
      const Scalar rel = error / scale;
      const Scalar rel_stress = stress_error / stress_scale;
      if (type == ElementType::Tet10) {
        worst_tet10 = std::max({worst_tet10, rel, rel_stress});
      } else if (perturbation > 0.0) {
        // On the undistorted grid the linear elements hit the nodal values
        // exactly (nodal superconvergence of a symmetric grid), so only the
        // distorted meshes test them.
        best_linear = std::min(best_linear, rel);
      }
      csv.raw_row({to_string(type), app::format(perturbation, 3),
                   app::format(static_cast<Scalar>(model.mesh().num_elements()), 6),
                   app::format(static_cast<Scalar>(model.dofs().num_dofs()), 6),
                   app::format(rel, 4), app::format(rel_stress, 4)});
      json::Value rec = json::Value::make_object();
      rec.set("element", json::Value::make_string(to_string(type)));
      rec.set("perturbation", json::Value::make_number(perturbation));
      rec.set("relative_displacement_error", json::Value::make_number(rel));
      rec.set("relative_stress_error", json::Value::make_number(rel_stress));
      records.push_back(rec);
    }
  }
  csv.close();
  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("verification"));
  block.set("records", records);
  block.set("note", json::Value::make_string(
                        "Pure bending about y prescribed on the boundary nodes of a "
                        "3 x 3 x 3 box, interior nodes solved. The field is quadratic, so "
                        "straight-sided Tet10 cells (edge nodes at the midpoints of the "
                        "distorted edges) contain it and must reproduce it and its linear "
                        "stress exactly; the Tet4 and the Hex8 cannot, and their errors show "
                        "the size of what the Tet10 removes. On the undistorted grid "
                        "(perturbation 0) the linear elements still hit the nodal values "
                        "exactly - nodal superconvergence of a symmetric grid, not a passed "
                        "patch test: the Tet4 element stresses are wrong there - so the "
                        "linear elements are judged on the distorted meshes."));
  summary.set("patch_test_quadratic", block);
  StudyOutcome outcome;
  outcome.name = "quadratic patch test (pure bending, distorted Tet10 meshes)";
  outcome.kind = "verification";
  outcome.metric = "max relative displacement / stress error of the Tet10";
  outcome.value = worst_tet10;
  outcome.tolerance = 1.0e-9;
  outcome.passed = worst_tet10 <= outcome.tolerance && best_linear > 1.0e-4;
  std::ostringstream note;
  note << "smallest error of the linear elements on the distorted meshes "
       << app::format(best_linear, 3) << " (they cannot pass)";
  outcome.note = note.str();
  return outcome;
}

/// Tet10 against Tet4 and Hex8 on the same solid cantilevers, loaded by a
/// consistent tip traction, against Timoshenko beam theory.
StudyOutcome study_mesh_convergence_tet10(const std::string& out_dir, json::Value& summary) {
  CsvWriter csv(path_join(out_dir, "mesh_convergence_tet10.csv"),
                {"element", "nx", "ny", "nz", "num_elements", "num_dofs", "h[m]", "tip_mean[m]",
                 "timoshenko[m]", "rel_error_timoshenko[-]", "solve_seconds[s]"});
  const SolidCantileverSpec spec;
  const Scalar reference = timoshenko_tip_3d(spec);
  const std::vector<std::array<Index, 3>> grids = {{10, 2, 1}, {20, 4, 2}, {40, 8, 4}};
  json::Value blocks = json::Value::make_object();
  Scalar tet10_finest = 0.0;
  Scalar tet10_coarsest = 0.0;
  std::ostringstream note;
  for (const ElementType type : {ElementType::Hex8, ElementType::Tet4, ElementType::Tet10}) {
    json::Value records = json::Value::make_array();
    Scalar finest = 0.0;
    for (std::size_t g = 0; g < grids.size(); ++g) {
      StructuredMeshSpec ms;
      ms.nx = grids[g][0];
      ms.ny = grids[g][1];
      ms.nz = grids[g][2];
      ms.lx = spec.length;
      ms.ly = spec.height;
      ms.lz = spec.width;
      Mesh mesh = type == ElementType::Hex8   ? make_structured_hex_mesh(ms)
                  : type == ElementType::Tet4 ? make_structured_tet_mesh(ms)
                                              : make_structured_tet10_mesh(ms);
      FemModel model = solid_cantilever(std::move(mesh), spec, Vector3(0.0, spec.tip_load, 0.0));
      Assembler assembler(model);
      StaticAnalysisOptions options;
      options.linear.type = LinearSolverType::SimplicialLdlt;
      Timer timer;
      StaticAnalysis analysis(model, assembler, options);
      const StaticSolution sol = analysis.solve_all().front();
      const Scalar seconds = timer.elapsed_seconds();
      const Scalar tip = tip_mean_deflection(model, sol.displacement, spec.length);
      const Scalar error = std::abs(tip - reference) / std::abs(reference);
      const Scalar h = spec.length / static_cast<Scalar>(ms.nx);
      finest = error;
      if (type == ElementType::Tet10 && g == 0) tet10_coarsest = error;
      csv.raw_row({to_string(type), app::format(static_cast<Scalar>(ms.nx), 6),
                   app::format(static_cast<Scalar>(ms.ny), 6),
                   app::format(static_cast<Scalar>(ms.nz), 6),
                   app::format(static_cast<Scalar>(model.mesh().num_elements()), 9),
                   app::format(static_cast<Scalar>(model.dofs().num_dofs()), 9),
                   app::format(h, 6), app::format(tip, 12), app::format(reference, 12),
                   app::format(error, 6), app::format(seconds, 4)});
      json::Value rec = json::Value::make_object();
      rec.set("nx", json::Value::make_number(ms.nx));
      rec.set("num_dofs", json::Value::make_number(model.dofs().num_dofs()));
      rec.set("tip_mean_m", json::Value::make_number(tip));
      rec.set("relative_error_vs_timoshenko", json::Value::make_number(error));
      records.push_back(rec);
    }
    if (type == ElementType::Tet10) tet10_finest = finest;
    json::Value block = json::Value::make_object();
    block.set("records", records);
    block.set("finest_mesh_relative_error_vs_timoshenko", json::Value::make_number(finest));
    blocks.set(to_string(type), block);
    note << to_string(type) << " finest " << app::format(100.0 * finest, 3) << " %; ";
  }
  csv.close();
  blocks.set("kind", json::Value::make_string("verification + validation (Timoshenko)"));
  blocks.set("note", json::Value::make_string(
                         "The 1 m x 0.1 m x 0.05 m cantilever of the Hex8 study (nu = 0), "
                         "loaded by a tip traction whose consistent nodal loads each "
                         "element computes, on the same grids for the three elements. The "
                         "Tet10 meshes split each cell into six Kuhn tetrahedra with edge "
                         "nodes at the edge midpoints."));
  summary.set("mesh_convergence_tet10", blocks);
  StudyOutcome outcome;
  outcome.name = "mesh convergence (Tet10 vs Tet4 and Hex8 cantilever)";
  outcome.kind = "verification + validation";
  outcome.metric = "Tet10 finest-mesh relative error vs Timoshenko";
  outcome.value = tet10_finest;
  outcome.tolerance = 0.01;
  outcome.passed = tet10_finest <= outcome.tolerance;
  note << "Tet10 on the coarsest (10 x 2 x 1) grid " << app::format(100.0 * tet10_coarsest, 3)
       << " %";
  outcome.note = note.str();
  return outcome;
}

/// Buckling load of a cantilever column (clamped base, axial tip traction)
/// against Euler with Engesser's shear correction.
StudyOutcome study_buckling_euler(const std::string& out_dir, json::Value& summary) {
  CsvWriter csv(path_join(out_dir, "buckling_euler.csv"),
                {"element", "nx", "ny", "nz", "num_dofs", "load_factor_1[N]", "load_factor_2[N]",
                 "engesser[N]", "euler[N]", "rel_error_engesser[-]", "iterations",
                 "transformed"});
  const Scalar length = 1.0;
  const Scalar h = 0.05;
  const Scalar t2d = 0.02;
  const Scalar e = 200.0e9;
  const Scalar nu = 0.3;
  json::Value blocks = json::Value::make_object();
  std::map<std::string, Scalar> finest;
  const auto run = [&](const ElementType type, Mesh mesh) {
    const int dim = mesh.dim();
    const Scalar t = dim == 2 ? t2d : 1.0;
    const Scalar width = dim == 2 ? t2d : h;
    const Scalar inertia = width * h * h * h / 12.0;
    const Scalar area = width * h;
    const Scalar reference = engesser_cantilever_load(e, nu, inertia, area, length);
    FemModel model(std::move(mesh), IsotropicMaterial(e, nu, 7850.0, "steel"), t,
                   dim == 2 ? StressState::PlaneStress : StressState::ThreeDimensional,
                   IntegrationOptions());
    DisplacementConstraint root;
    Selector box;
    box.kind = SelectorKind::Box;
    box.xmax = 0.0;
    root.region.members.push_back(box);
    root.fix_x = root.fix_y = true;
    root.fix_z = dim == 3;
    model.constraints().push_back(root);
    LoadCaseSpec load;
    load.name = "axial";
    TractionLoadSpec top;
    Selector end;
    end.kind = SelectorKind::Box;
    end.xmin = length;
    top.region.members.push_back(end);
    top.traction = Vector3(-1.0 / area, 0.0, 0.0);  // 1 N: lambda is the load in N
    load.tractions.push_back(top);
    model.load_case_specs().push_back(load);
    model.finalize();
    Assembler assembler(model);
    BucklingOptions options;
    options.num_modes = 2;
    options.linear.type = LinearSolverType::SimplicialLdlt;
    const BucklingResult r = analyse_buckling(model, assembler, 0, options);
    const Scalar error = r.load_factors(0) / reference - 1.0;
    const StructuredGridInfo& info = *model.mesh().structured_info();
    (void)info;
    return std::make_tuple(r, reference, error, model.dofs().num_dofs(),
                           euler_cantilever_load(e, inertia, length));
  };
  const auto record = [&](const ElementType type, const StructuredMeshSpec& ms, int dim,
                          json::Value& records) {
    Mesh mesh = dim == 2                    ? make_structured_quad_mesh(ms)
                : type == ElementType::Hex8 ? make_structured_hex_mesh(ms)
                : type == ElementType::Tet4 ? make_structured_tet_mesh(ms)
                                            : make_structured_tet10_mesh(ms);
    const std::string name = dim == 2 ? "Quad4" : to_string(type);
    const auto [r, reference, error, dofs, euler] = run(type, std::move(mesh));
    csv.raw_row({name, app::format(static_cast<Scalar>(ms.nx), 6),
                 app::format(static_cast<Scalar>(ms.ny), 6),
                 app::format(static_cast<Scalar>(dim == 3 ? ms.nz : 0), 6),
                 app::format(static_cast<Scalar>(dofs), 9), app::format(r.load_factors(0), 10),
                 app::format(r.load_factors(1), 10), app::format(reference, 10),
                 app::format(euler, 10), app::format(error, 6),
                 app::format(static_cast<Scalar>(r.iterations), 6), r.transformed ? "1" : "0"});
    json::Value rec = json::Value::make_object();
    rec.set("nx", json::Value::make_number(ms.nx));
    rec.set("num_dofs", json::Value::make_number(dofs));
    rec.set("load_factor_1_N", json::Value::make_number(r.load_factors(0)));
    rec.set("engesser_N", json::Value::make_number(reference));
    rec.set("relative_error_vs_engesser", json::Value::make_number(error));
    records.push_back(rec);
    finest[name] = error;
  };
  {
    json::Value records = json::Value::make_array();
    for (Index nx : {20, 40, 80, 160}) {
      StructuredMeshSpec ms;
      ms.nx = nx;
      ms.ny = nx / 20 * 2;
      ms.lx = length;
      ms.ly = h;
      record(ElementType::Quad4, ms, 2, records);
    }
    blocks.set("Quad4", records);
  }
  for (const ElementType type : {ElementType::Hex8, ElementType::Tet4, ElementType::Tet10}) {
    json::Value records = json::Value::make_array();
    for (const std::array<Index, 3>& g :
         std::vector<std::array<Index, 3>>{{10, 1, 1}, {20, 2, 2}, {40, 4, 4}}) {
      StructuredMeshSpec ms;
      ms.nx = g[0];
      ms.ny = g[1];
      ms.nz = g[2];
      ms.lx = length;
      ms.ly = h;
      ms.lz = h;
      record(type, ms, 3, records);
    }
    blocks.set(to_string(type), records);
  }
  csv.close();
  blocks.set("kind", json::Value::make_string("verification + validation (Euler-Engesser)"));
  blocks.set("note", json::Value::make_string(
                         "A 1 m column of 0.05 m square (3-D) or 0.05 m x 0.02 m (plane "
                         "stress) section, clamped at x = 0, under a unit axial compressive "
                         "tip traction, so the first load factor is the critical load in "
                         "newtons. Reference: Euler P = pi^2 E I / (4 L^2) with Engesser's "
                         "shear correction 1/P_s = 1/P + 1/(5/6 G A); displacement elements "
                         "converge to it from above."));
  summary.set("buckling_euler", blocks);
  StudyOutcome outcome;
  outcome.name = "linear buckling of a cantilever column vs Euler-Engesser";
  outcome.kind = "verification + validation";
  outcome.metric = "larger finest-mesh error of Q4 and Tet10";
  outcome.value = std::max(std::abs(finest["Quad4"]), std::abs(finest["Tet10"]));
  outcome.tolerance = 0.01;
  outcome.passed = outcome.value <= outcome.tolerance;
  std::ostringstream note;
  for (const auto& kv : finest) {
    note << kv.first << " " << app::format(100.0 * kv.second, 3) << " %; ";
  }
  outcome.note = note.str();
  return outcome;
}

/// A clamped-left block with a tip load for the gradient checks below.
FemModel gradient_block(Mesh mesh, Scalar compression) {
  const int dim = mesh.dim();
  const BoundingBox bb = mesh.bounding_box();
  FemModel model(std::move(mesh), IsotropicMaterial(70.0e9, 0.3, 2700.0, "block"),
                 dim == 2 ? 0.01 : 1.0,
                 dim == 2 ? StressState::PlaneStress : StressState::ThreeDimensional,
                 IntegrationOptions());
  DisplacementConstraint root;
  Selector box;
  box.kind = SelectorKind::Box;
  box.xmax = bb.lower.x();
  root.region.members.push_back(box);
  root.fix_x = root.fix_y = true;
  root.fix_z = dim == 3;
  model.constraints().push_back(root);
  LoadCaseSpec load;
  load.name = "tip";
  TractionLoadSpec end;
  Selector face;
  face.kind = SelectorKind::Box;
  face.xmin = bb.upper.x();
  end.region.members.push_back(face);
  end.traction = Vector3(-compression, -0.1 * compression, 0.0);
  load.tractions.push_back(end);
  model.load_case_specs().push_back(load);
  model.finalize();
  return model;
}

/// Central differences of the aggregated buckling constraint and of lambda_1,
/// through the density filter and the projection.
StudyOutcome study_sensitivity_buckling(const std::string& out_dir, json::Value& summary,
                                        Scalar tolerance) {
  CsvWriter csv(path_join(out_dir, "sensitivity_buckling.csv"),
                {"mesh", "beta[-]", "step[-]", "num_tested", "ks_max_scaled_error[-]",
                 "lambda1_max_scaled_error[-]", "lambda_1[-]", "lambda_2[-]"});
  json::Value records = json::Value::make_array();
  Scalar worst = 0.0;
  struct Case {
    std::string name;
    Mesh mesh;
    Scalar beta;
  };
  std::vector<Case> cases;
  {
    StructuredMeshSpec ms;
    ms.nx = 16;
    ms.ny = 8;
    ms.lx = 0.8;
    ms.ly = 0.4;
    cases.push_back({"Q4 16 x 8", make_structured_quad_mesh(ms), 0.0});
    cases.push_back({"Q4 16 x 8", make_structured_quad_mesh(ms), 4.0});
    StructuredMeshSpec ms3;
    ms3.nx = 5;
    ms3.ny = ms3.nz = 2;
    ms3.lx = 0.5;
    ms3.ly = ms3.lz = 0.2;
    cases.push_back({"Hex8 5 x 2 x 2", make_structured_hex_mesh(ms3), 0.0});
  }
  for (Case& c : cases) {
    FemModel model = gradient_block(std::move(c.mesh), 1.0e6);
    Assembler assembler(model);
    const DensityFilter filter(model.mesh(), FilterType::Density,
                               1.5 * model.mesh().mean_element_size());
    DesignDomain domain(model, 0.5, 0.5, {});
    StaticAnalysisOptions options;
    options.linear.type = LinearSolverType::SimplicialLdlt;
    ComplianceObjective objective(model, assembler, filter, domain, SimpOptions(), options);
    if (c.beta > 0.0) objective.set_projection(c.beta, 0.5);
    BucklingConstraintOptions bo;
    bo.enabled = true;
    bo.min_load_factor = 1.0;
    bo.num_modes = 3;
    bo.ks_parameter = 20.0;
    bo.eigen.tolerance = 1.0e-14;
    bo.eigen.residual_tolerance = 1.0e-11;
    bo.eigen.max_iterations = 5000;
    BucklingConstraint buckling(model, assembler, bo);
    Vector x = domain.initial_design();
    for (Index e = 0; e < domain.num_elements(); ++e) {
      const Vector3 p = model.mesh().element_centroid(e);
      x(e) = 0.5 + 0.3 * std::sin(11.0 * p.x()) * std::cos(7.0 * p.y() + 5.0 * p.z());
    }
    domain.clamp(x);
    const ObjectiveEvaluation eval = objective.evaluate(x, true);
    const BucklingEvaluation be = buckling.evaluate(objective, eval, 0, true);
    const Vector dl = objective.chain_to_design(eval, be.dlambda_dphysical.front());
    Scalar best = std::numeric_limits<Scalar>::max();
    for (Scalar step : {1.0e-3, 1.0e-4, 1.0e-5}) {
      Scalar ks_err = 0.0;
      Scalar l_err = 0.0;
      const Scalar ks_floor = 1.0e-3 * be.dg_dx.cwiseAbs().maxCoeff();
      const Scalar l_floor = 1.0e-3 * dl.cwiseAbs().maxCoeff();
      Index tested = 0;
      for (Index e = 0; e < domain.num_elements(); e += 3) {
        Vector xp = x;
        Vector xm = x;
        xp(e) += step;
        xm(e) -= step;
        const BucklingEvaluation bp =
            buckling.evaluate(objective, objective.evaluate(xp, false), 0, false);
        const BucklingEvaluation bm =
            buckling.evaluate(objective, objective.evaluate(xm, false), 0, false);
        const Scalar fd_ks = (bp.constraint - bm.constraint) / (2.0 * step);
        const Scalar fd_l = (bp.load_factors(0) - bm.load_factors(0)) / (2.0 * step);
        ks_err = std::max(ks_err, std::abs(fd_ks - be.dg_dx(e)) /
                                      std::max({std::abs(fd_ks), std::abs(be.dg_dx(e)), ks_floor}));
        l_err = std::max(l_err, std::abs(fd_l - dl(e)) /
                                    std::max({std::abs(fd_l), std::abs(dl(e)), l_floor}));
        ++tested;
      }
      best = std::min(best, std::max(ks_err, l_err));
      csv.raw_row({c.name, app::format(c.beta, 3), app::format(step, 3),
                   app::format(static_cast<Scalar>(tested), 6), app::format(ks_err, 4),
                   app::format(l_err, 4), app::format(be.load_factors(0), 8),
                   app::format(be.load_factors(1), 8)});
    }
    worst = std::max(worst, best);
    json::Value rec = json::Value::make_object();
    rec.set("mesh", json::Value::make_string(c.name));
    rec.set("beta", json::Value::make_number(c.beta));
    rec.set("lambda_1", json::Value::make_number(be.load_factors(0)));
    rec.set("best_max_scaled_error", json::Value::make_number(best));
    records.push_back(rec);
  }
  csv.close();
  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("verification"));
  block.set("records", records);
  block.set("note", json::Value::make_string(
                        "Aggregated buckling constraint KS(1 / lambda_i) - 1 over the lowest "
                        "three load factors and lambda_1 itself, adjoint gradient against "
                        "central differences through the density filter (and the projection "
                        "at beta = 4) on every third element; entries judged against "
                        "max(|analytical|, |FD|, 1e-3 ||gradient||_inf), eigensolves "
                        "converged to 1e-14."));
  summary.set("sensitivity_buckling", block);
  StudyOutcome outcome;
  outcome.name = "buckling-constraint sensitivity (Q4 and Hex8)";
  outcome.kind = "verification";
  outcome.metric = "worst case of the best (over steps) max scaled error";
  outcome.value = worst;
  outcome.tolerance = tolerance;
  outcome.passed = worst <= tolerance;
  return outcome;
}

/// Central differences of the compliance through the overhang filter and the
/// projection.
StudyOutcome study_sensitivity_overhang(const std::string& out_dir, json::Value& summary,
                                        Scalar tolerance) {
  CsvWriter csv(path_join(out_dir, "sensitivity_overhang.csv"),
                {"mesh", "build_direction", "beta[-]", "step[-]", "num_tested",
                 "max_scaled_error[-]", "directional_relative_error[-]"});
  json::Value records = json::Value::make_array();
  Scalar worst = 0.0;
  for (int dim : {2, 3}) {
    StructuredMeshSpec ms;
    ms.nx = dim == 2 ? 16 : 6;
    ms.ny = dim == 2 ? 8 : 4;
    ms.nz = 4;
    ms.lx = dim == 2 ? 0.8 : 0.6;
    ms.ly = dim == 2 ? 0.4 : 0.4;
    ms.lz = 0.4;
    FemModel model = gradient_block(
        dim == 2 ? make_structured_quad_mesh(ms) : make_structured_hex_mesh(ms), 1.0e5);
    Assembler assembler(model);
    const DensityFilter filter(model.mesh(), FilterType::Density,
                               1.5 * model.mesh().mean_element_size());
    DesignDomain domain(model, 0.5, 0.5, {});
    StaticAnalysisOptions options;
    options.linear.type = LinearSolverType::SimplicialLdlt;
    ComplianceObjective objective(model, assembler, filter, domain, SimpOptions(), options);
    OverhangOptions oh;
    oh.filter = true;
    oh.direction = parse_build_direction(dim == 2 ? "+y" : "+z", dim);
    const OverhangFilter am(model.mesh(), oh);
    objective.set_overhang(&am);
    Vector x = domain.initial_design();
    for (Index e = 0; e < domain.num_elements(); ++e) {
      const Vector3 p = model.mesh().element_centroid(e);
      x(e) = 0.5 + 0.3 * std::sin(11.0 * p.x()) * std::cos(7.0 * p.y() + 5.0 * p.z());
    }
    domain.clamp(x);
    for (Scalar beta : {0.0, 4.0, 16.0}) {
      if (beta > 0.0) {
        objective.set_projection(beta, 0.5);
      } else {
        objective.disable_projection();
      }
      Scalar best = std::numeric_limits<Scalar>::max();
      for (Scalar step : {1.0e-3, 1.0e-4, 1.0e-5}) {
        const SensitivityCheckResult check =
            verify_sensitivities(objective, domain, x, {}, step, tolerance);
        best = std::min(best, std::max(check.max_scaled_error, check.directional_relative_error));
        csv.raw_row({dim == 2 ? "Q4" : "Hex8", oh.direction.label(), app::format(beta, 3),
                     app::format(step, 3), app::format(static_cast<Scalar>(check.num_tested), 6),
                     app::format(check.max_scaled_error, 4),
                     app::format(check.directional_relative_error, 4)});
      }
      worst = std::max(worst, best);
      json::Value rec = json::Value::make_object();
      rec.set("mesh", json::Value::make_string(dim == 2 ? "Q4" : "Hex8"));
      rec.set("beta", json::Value::make_number(beta));
      rec.set("best_error", json::Value::make_number(best));
      records.push_back(rec);
    }
  }
  csv.close();
  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("verification"));
  block.set("records", records);
  block.set("note", json::Value::make_string(
                        "Compliance gradient through the density filter, the overhang "
                        "filter (its adjoint recursion) and, at beta > 0, the projection, "
                        "against central differences; entries judged against max(|analytical|, "
                        "|FD|, 1e-3 ||gradient||_inf), plus the directional derivative."));
  summary.set("sensitivity_overhang", block);
  StudyOutcome outcome;
  outcome.name = "sensitivity through the overhang filter (Q4 and Hex8)";
  outcome.kind = "verification";
  outcome.metric = "worst over beta of the best max(scaled entry error, directional error)";
  outcome.value = worst;
  outcome.tolerance = tolerance;
  outcome.passed = worst <= tolerance;
  return outcome;
}

}  // namespace

int main(int argc, char** argv) {
  return app::run_guarded([&]() -> int {
    const std::vector<std::string> known = {"study", "output", "verbosity",
                                            "sensitivity-tolerance", "help"};
    app::CommandLine cli(argc, argv, known);
    if (cli.has("help")) {
      return app::print_usage(
          "sparlab_verify", "[--study <name>] [--output <dir>]",
          {{"--study <name>",
            "all (default) | sensitivity | mesh-convergence | modal | "
            "solver-agreement | patch-test | patch-test-3d | mesh-convergence-3d | "
            "sensitivity-3d | modal-3d | patch-test-simplex | mesh-convergence-simplex | "
            "multigrid | sensitivity-projection | patch-test-quadratic | "
            "mesh-convergence-tet10 | buckling-euler | sensitivity-buckling | "
            "sensitivity-overhang | lame-cylinder | rotating-disk | thermal-cylinder | "
            "bimetal-strip | self-weight | elastica | hyperelastic-cylinder | "
            "arch-snap-through | plastic-cylinder | plastic-bending | plastic-cycle | "
            "chaboche-cycle | hill-directional | logarithmic-uniaxial | logarithmic-tube | "
            "transient-modal | rod-harmonic | rod-transient | nonlinear-oscillator"},
           {"--output <dir>", "output directory (default results/verification)"},
           {"--sensitivity-tolerance <t>",
            "pass threshold on the max relative gradient error (default 1e-5)"},
           {"--verbosity <lvl>", "trace|debug|info|warn|error|silent"},
           {"--help", "show this message"}});
    }
    app::apply_verbosity(cli);

    const std::string study = cli.value("study", "all");
    const std::string out_dir = cli.value("output", "results/verification");
    const Scalar sensitivity_tolerance = cli.number("sensitivity-tolerance", 1.0e-5);
    ensure_directory(out_dir);

    const bool all = study == "all";
    json::Value summary = json::Value::make_object();
    {
      json::Value header = json::Value::make_object();
      header.set("study", json::Value::make_string(study));
      header.set("units", json::Value::make_string("SI: m, N, Pa, kg, Hz, J"));
      header.set("terminology",
                 json::Value::make_string(
                     "verification = the code solves the stated discrete equations "
                     "correctly (judged against exact answers); validation = the "
                     "modelling assumptions are appropriate (judged against an "
                     "independent theory, where a finite gap is expected)"));
      summary.set("about", header);
    }

    Timer wall;
    std::vector<StudyOutcome> outcomes;
    if (all || study == "patch-test") {
      outcomes.push_back(study_patch_test(out_dir, summary));
    }
    if (all || study == "solver-agreement") {
      outcomes.push_back(study_solver_agreement(out_dir, summary));
    }
    if (all || study == "sensitivity") {
      outcomes.push_back(study_sensitivity(out_dir, summary, sensitivity_tolerance));
    }
    if (all || study == "mesh-convergence") {
      outcomes.push_back(study_mesh_convergence(out_dir, summary));
    }
    if (all || study == "modal") {
      outcomes.push_back(study_modal(out_dir, summary));
    }
    if (all || study == "patch-test-3d") {
      outcomes.push_back(study_patch_test_3d(out_dir, summary));
    }
    if (all || study == "sensitivity-3d") {
      outcomes.push_back(study_sensitivity_3d(out_dir, summary, sensitivity_tolerance));
    }
    if (all || study == "mesh-convergence-3d") {
      outcomes.push_back(study_mesh_convergence_3d(out_dir, summary));
    }
    if (all || study == "modal-3d") {
      outcomes.push_back(study_modal_3d(out_dir, summary));
    }
    if (all || study == "patch-test-simplex") {
      outcomes.push_back(study_patch_test_simplex(out_dir, summary));
    }
    if (all || study == "mesh-convergence-simplex") {
      outcomes.push_back(study_mesh_convergence_simplex(out_dir, summary));
    }
    if (all || study == "multigrid") {
      outcomes.push_back(study_multigrid(out_dir, summary));
    }
    if (all || study == "sensitivity-projection") {
      outcomes.push_back(study_sensitivity_projection(out_dir, summary, sensitivity_tolerance));
    }
    if (all || study == "patch-test-quadratic") {
      outcomes.push_back(study_patch_test_quadratic(out_dir, summary));
    }
    if (all || study == "mesh-convergence-tet10") {
      outcomes.push_back(study_mesh_convergence_tet10(out_dir, summary));
    }
    if (all || study == "buckling-euler") {
      outcomes.push_back(study_buckling_euler(out_dir, summary));
    }
    if (all || study == "sensitivity-buckling") {
      outcomes.push_back(study_sensitivity_buckling(out_dir, summary, sensitivity_tolerance));
    }
    if (all || study == "sensitivity-overhang") {
      outcomes.push_back(study_sensitivity_overhang(out_dir, summary, sensitivity_tolerance));
    }
    if (all || study == "lame-cylinder") {
      outcomes.push_back(verify::study_lame_cylinder(out_dir, summary));
    }
    if (all || study == "rotating-disk") {
      outcomes.push_back(verify::study_rotating_disk(out_dir, summary));
    }
    if (all || study == "thermal-cylinder") {
      outcomes.push_back(verify::study_thermal_cylinder(out_dir, summary));
    }
    if (all || study == "bimetal-strip") {
      outcomes.push_back(verify::study_bimetal_strip(out_dir, summary));
    }
    if (all || study == "self-weight") {
      outcomes.push_back(verify::study_self_weight(out_dir, summary));
    }
    if (all || study == "elastica") {
      outcomes.push_back(verify::study_elastica(out_dir, summary));
    }
    if (all || study == "hyperelastic-cylinder") {
      outcomes.push_back(verify::study_hyperelastic_cylinder(out_dir, summary));
    }
    if (all || study == "arch-snap-through") {
      outcomes.push_back(verify::study_arch_snap_through(out_dir, summary));
    }
    if (all || study == "plastic-cylinder") {
      outcomes.push_back(verify::study_plastic_cylinder(out_dir, summary));
    }
    if (all || study == "plastic-bending") {
      outcomes.push_back(verify::study_plastic_bending(out_dir, summary));
    }
    if (all || study == "plastic-cycle") {
      outcomes.push_back(verify::study_plastic_cycle(out_dir, summary));
    }
    if (all || study == "chaboche-cycle") {
      outcomes.push_back(verify::study_chaboche_cycle(out_dir, summary));
    }
    if (all || study == "hill-directional") {
      outcomes.push_back(verify::study_hill_directional(out_dir, summary));
    }
    if (all || study == "logarithmic-uniaxial") {
      outcomes.push_back(verify::study_logarithmic_uniaxial(out_dir, summary));
    }
    if (all || study == "logarithmic-tube") {
      outcomes.push_back(verify::study_logarithmic_tube(out_dir, summary));
    }
    if (all || study == "transient-modal") {
      outcomes.push_back(verify::study_transient_modal(out_dir, summary));
    }
    if (all || study == "rod-harmonic") {
      outcomes.push_back(verify::study_rod_harmonic(out_dir, summary));
    }
    if (all || study == "rod-transient") {
      outcomes.push_back(verify::study_rod_transient(out_dir, summary));
    }
    if (all || study == "nonlinear-oscillator") {
      outcomes.push_back(verify::study_nonlinear_oscillator(out_dir, summary));
    }
    if (outcomes.empty()) {
      throw ConfigError("unknown study '" + study +
                        "'; run with --help to list the available studies");
    }

    json::Value results = json::Value::make_array();
    bool all_passed = true;
    for (const StudyOutcome& o : outcomes) {
      json::Value rec = json::Value::make_object();
      rec.set("study", json::Value::make_string(o.name));
      rec.set("kind", json::Value::make_string(o.kind));
      rec.set("metric", json::Value::make_string(o.metric));
      rec.set("value", json::Value::make_number(o.value));
      rec.set("tolerance", json::Value::make_number(o.tolerance));
      rec.set("passed", json::Value::make_bool(o.passed));
      if (!o.note.empty()) rec.set("note", json::Value::make_string(o.note));
      results.push_back(rec);
      all_passed = all_passed && o.passed;
    }
    summary.set("outcomes", results);
    summary.set("all_passed", json::Value::make_bool(all_passed));
    summary.set("total_seconds", json::Value::make_number(wall.elapsed_seconds()));

    std::ofstream out(path_join(out_dir, "summary.json"));
    if (!out) throw IoError("cannot write the verification summary");
    out << json::dump(summary, 2) << '\n';
    out.close();

    std::cout << "\nverification / validation results\n";
    std::cout << std::string(96, '-') << "\n";
    for (const StudyOutcome& o : outcomes) {
      std::cout << (o.passed ? "  PASS  " : "  FAIL  ") << o.name << "\n";
      std::cout << "          " << o.kind << " | " << o.metric << " = "
                << app::format(o.value, 4) << " (tolerance "
                << app::format(o.tolerance, 3) << ")\n";
      if (!o.note.empty()) std::cout << "          " << o.note << "\n";
    }
    std::cout << std::string(96, '-') << "\n";
    std::cout << (all_passed ? "all studies passed" : "AT LEAST ONE STUDY FAILED")
              << " in " << app::format(wall.elapsed_seconds(), 4) << " s\n";
    std::cout << "results: " << out_dir << "\n";
    return all_passed ? 0 : 1;
  });
}
