/// \file sparlab_topopt.cpp
/// \brief SIMP minimum-compliance topology optimisation of one configuration.
///
/// Pipeline:
///   1. build the FE model and check that the *solid* model is well posed;
///   2. build the density filter and the design domain (passive regions);
///   3. run the optimality-criteria loop;
///   4. recover stresses on the final design;
///   5. interpret the density field as a solid body (threshold + largest
///      connected group) and, where the boundary conditions still apply,
///      analyse that body directly;
///   6. compare natural frequencies of the full solid domain, an equal-mass
///      uniform plate, and the optimised topology;
///   7. when a buckling check is configured, compare the buckling load
///      factors of the full solid domain and the interpreted structure - the
///      part that would be exported;
///   8. write every artefact plus a summary, including the "before" (design
///      domain) and "after" (thresholded structure) geometries as VTK + STL.
///
/// Step 6 uses a property specific to the 2-D idealisation: uniformly scaling
/// the thickness scales K and M equally, so an equal-mass *uniform* plate has
/// exactly the same natural frequencies as the full solid domain while its
/// compliance rises by 1/nu. The run verifies that numerically rather than
/// asserting it, and the equal-mass compliance is the baseline the optimised
/// design is judged against. A 3-D solid has no thickness to scale, so a 3-D
/// run reports the full solid as its only reference and says so.

#include "AppSupport.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Timer.hpp"
#include "sparlab/fem/Assembler.hpp"
#include "sparlab/fem/Buckling.hpp"
#include "sparlab/fem/ModalAnalysis.hpp"
#include "sparlab/fem/ModelDiagnostics.hpp"
#include "sparlab/fem/StaticAnalysis.hpp"
#include "sparlab/fem/StressRecovery.hpp"
#include "sparlab/io/Config.hpp"
#include "sparlab/io/ResultWriter.hpp"
#include "sparlab/mesh/SubMesh.hpp"
#include "sparlab/topopt/LengthScale.hpp"
#include "sparlab/topopt/OverhangFilter.hpp"
#include "sparlab/topopt/TopologyOptimizer.hpp"

#include <iostream>
#include <memory>

using namespace sparlab;

namespace {

/// Build a solid FE model on an extracted sub-mesh, re-applying the geometric
/// boundary conditions and (optionally) the load cases.
struct SubModelBuild {
  std::unique_ptr<FemModel> model;
  bool loads_applied = false;
  std::string note;
};

/// Buckling check of the listed load cases of a plain (solid) model, with one
/// factorisation shared by the static solves and the eigensolves.
std::vector<BucklingResult> check_buckling(const FemModel& model, const Assembler& assembler,
                                           const Configuration& config,
                                           const std::vector<std::size_t>& cases) {
  StaticAnalysis analysis(model, assembler, config.analysis);
  const std::vector<StaticSolution> solutions = analysis.solve_all();
  const DofManager& dofs = model.dofs();
  const FreeSolve solve = [&](const Vector& b) {
    return dofs.restrict_to_free(analysis.solve_homogeneous(dofs.expand(b)));
  };
  std::vector<BucklingResult> out;
  for (std::size_t l : cases) {
    const SparseMatrix k_g =
        assemble_geometric_stiffness(model, assembler, solutions[l].displacement);
    BucklingResult result = solve_buckling(model, assembler, analysis.stiffness(), k_g,
                                           config.buckling.options, solve);
    result.load_case = solutions[l].load_case_name;
    result.linear_solver = analysis.solver().name();
    out.push_back(std::move(result));
  }
  return out;
}

SubModelBuild build_solid_submodel(const Configuration& config, Mesh sub_mesh,
                                   bool with_loads) {
  SubModelBuild out;
  auto model = std::make_unique<FemModel>(std::move(sub_mesh), config.material(),
                                          config.thickness, config.stress_state,
                                          config.integration);
  model->constraints() = config.constraints;
  if (with_loads) {
    model->load_case_specs() = config.load_cases;
    try {
      model->finalize(/*require_load_cases=*/true);
      out.loads_applied = true;
    } catch (const ConfigError& e) {
      // The load region no longer intersects the retained material. Report it
      // instead of quietly skipping the comparison.
      out.note = std::string("load cases could not be applied to the interpreted "
                             "solid structure: ") +
                 e.what();
      model = std::make_unique<FemModel>(Mesh(model->mesh()), config.material(),
                                         config.thickness, config.stress_state,
                                         config.integration);
      model->constraints() = config.constraints;
      model->finalize(/*require_load_cases=*/false);
    }
  } else {
    model->finalize(/*require_load_cases=*/false);
  }
  out.model = std::move(model);
  return out;
}

}  // namespace

int main(int argc, char** argv) {
  return app::run_guarded([&]() -> int {
    const std::vector<std::string> known = {
        "config",        "output",      "verbosity",   "strict-config",
        "volume-fraction", "penalty",   "filter-radius", "filter-radius-elements",
        "filter-type",   "max-iterations", "nx",       "ny",           "nz",
        "youngs-modulus", "load-weights", "modes",     "no-vtk",
        "no-csv",        "tag",         "method",      "stress-limit", "no-stress",
        "solver",        "projection",  "no-projection", "beta-max",   "buckling",
        "min-load-factor", "no-buckling-constraint", "robust", "no-robust",
        "overhang",      "no-overhang-filter", "help"};
    app::CommandLine cli(argc, argv, known);
    if (cli.has("help") || argc == 1) {
      return app::print_usage(
          "sparlab_topopt", "--config <deck.json> [--output <dir>] [overrides]",
          {{"--config <file>", "input deck with a topology section"},
           {"--output <dir>", "output directory (default results/<case>)"},
           {"--volume-fraction <v>", "override topology.volume_fraction"},
           {"--penalty <p>", "override topology.simp.penalty"},
           {"--filter-radius <r>", "override the filter radius in metres"},
           {"--filter-radius-elements <k>", "filter radius as a multiple of the cell size"},
           {"--filter-type <t>", "none|density|sensitivity"},
           {"--max-iterations <n>", "override optimizer.max_iterations"},
           {"--nx <n> --ny <n> [--nz <n>]", "override the mesh resolution"},
           {"--youngs-modulus <E>", "override material.youngs_modulus [Pa]"},
           {"--load-weights <w1,w2,...>", "override the per-load-case weights"},
           {"--modes <n>", "enable modal analysis with n modes"},
           {"--method <oc|mma>", "override optimizer.method"},
           {"--stress-limit <Pa>", "enable the aggregated stress constraint at this "
                                   "allowable von Mises stress (switches to MMA)"},
           {"--no-stress", "disable the deck's stress constraint"},
           {"--buckling <n>", "check n buckling modes of the full solid domain and the "
                              "interpreted structure"},
           {"--min-load-factor <l>", "enable the buckling constraint lambda >= l "
                                     "(switches to MMA)"},
           {"--no-buckling-constraint", "disable the deck's buckling constraint"},
           {"--robust / --no-robust", "switch the robust (eroded/dilated) formulation on "
                                      "(with the projection) or off"},
           {"--overhang <dir>", "apply the overhang filter for build direction dir "
                                "(+x, -y, +z, ...)"},
           {"--no-overhang-filter", "keep the deck's overhang check but not its filter"},
           {"--solver <type>", "override solver.linear.type (simplicial_ldlt, amg_cg, "
                               "auto, ...)"},
           {"--projection / --no-projection", "switch the Heaviside projection on (with "
                                              "the deck's or default settings) or off"},
           {"--beta-max <b>", "override topology.projection.beta_max (enables it)"},
           {"--tag <name>", "suffix appended to the case name in summary.json"},
           {"--no-vtk / --no-csv", "skip the corresponding output"},
           {"--strict-config", "treat unknown configuration keys as errors"},
           {"--verbosity <lvl>", "trace|debug|info|warn|error|silent"},
           {"--help", "show this message"}});
    }
    app::apply_verbosity(cli);

    Configuration config =
        load_configuration(cli.require("config"), cli.has("strict-config"));

    // ---- command-line overrides (used by the parametric study) -------------
    if (cli.has("volume-fraction")) {
      config.topology.volume_fraction = cli.number("volume-fraction", 0.4);
    }
    if (cli.has("penalty")) {
      config.topology.optimizer.simp.penalty = cli.number("penalty", 3.0);
    }
    if (cli.has("filter-radius")) {
      config.topology.filter_radius = cli.number("filter-radius", 0.0);
      config.topology.filter_radius_elements = 0.0;
    }
    if (cli.has("filter-radius-elements")) {
      config.topology.filter_radius_elements =
          cli.number("filter-radius-elements", 1.5);
      config.topology.filter_radius = 0.0;
    }
    if (cli.has("filter-type")) {
      config.topology.filter_type = parse_filter_type(cli.value("filter-type"));
    }
    if (cli.has("max-iterations")) {
      config.topology.optimizer.max_iterations = cli.integer("max-iterations", 200);
    }
    if ((cli.has("nx") || cli.has("ny") || cli.has("nz")) &&
        !is_structured(config.mesh_kind)) {
      throw ConfigError("--nx/--ny/--nz set the resolution of a generated box mesh, but "
                        "this deck reads its mesh from '" + config.mesh_file.path +
                        "'; refine that mesh in the mesher instead");
    }
    if (cli.has("nx")) config.mesh_spec.nx = cli.integer("nx", config.mesh_spec.nx);
    if (cli.has("ny")) config.mesh_spec.ny = cli.integer("ny", config.mesh_spec.ny);
    if (cli.has("nz")) {
      if (config.dim() != 3) {
        throw ConfigError("--nz applies to a solid (structured_hex or structured_tet) "
                          "mesh only");
      }
      config.mesh_spec.nz = cli.integer("nz", config.mesh_spec.nz);
    }
    if (cli.has("youngs-modulus")) {
      config.set_material(config.material().with_youngs_modulus(
          cli.number("youngs-modulus", config.material().youngs_modulus())));
    }
    if (cli.has("load-weights")) {
      const std::vector<Scalar> weights = cli.number_list("load-weights");
      if (weights.size() != config.load_cases.size()) {
        std::ostringstream os;
        os << "--load-weights supplied " << weights.size() << " weight(s) but the deck "
           << "defines " << config.load_cases.size() << " load case(s)";
        throw ConfigError(os.str());
      }
      for (std::size_t l = 0; l < weights.size(); ++l) {
        config.load_cases[l].weight = weights[l];
      }
    }
    if (cli.has("modes")) {
      config.modal.enabled = true;
      config.modal.options.num_modes = cli.integer("modes", 6);
    }
    if (cli.has("method")) {
      config.topology.optimizer.method = parse_optimizer_method(cli.value("method"));
    }
    if (cli.has("stress-limit")) {
      config.topology.optimizer.stress.enabled = true;
      config.topology.optimizer.stress.limit = cli.number("stress-limit", 0.0);
      config.topology.optimizer.method = OptimizerMethod::MMA;
    }
    if (cli.has("no-stress")) config.topology.optimizer.stress.enabled = false;
    if (cli.has("buckling")) {
      config.buckling.enabled = true;
      config.buckling.options.num_modes = cli.integer("buckling", 4);
    }
    if (cli.has("min-load-factor")) {
      config.topology.optimizer.buckling.enabled = true;
      config.topology.optimizer.buckling.min_load_factor = cli.number("min-load-factor", 1.0);
      config.topology.optimizer.method = OptimizerMethod::MMA;
      config.topology.optimizer.buckling.validate();
    }
    if (cli.has("no-buckling-constraint")) config.topology.optimizer.buckling.enabled = false;
    if (cli.has("solver")) {
      config.analysis.linear.type = parse_linear_solver_type(cli.value("solver"));
      config.topology.optimizer.analysis.linear.type = config.analysis.linear.type;
      config.modal.options.linear.type = config.analysis.linear.type;
      config.buckling.options.linear.type = config.analysis.linear.type;
      config.topology.optimizer.buckling.eigen.linear.type = config.analysis.linear.type;
    }
    if (cli.has("projection") && cli.has("no-projection")) {
      throw ConfigError("--projection and --no-projection contradict each other");
    }
    if (cli.has("projection")) config.topology.optimizer.projection.enabled = true;
    if (cli.has("beta-max")) {
      config.topology.optimizer.projection.enabled = true;
      config.topology.optimizer.projection.beta_max = cli.number("beta-max", 32.0);
    }
    if (cli.has("no-projection")) config.topology.optimizer.projection.enabled = false;
    if (cli.has("robust") && cli.has("no-robust")) {
      throw ConfigError("--robust and --no-robust contradict each other");
    }
    if (cli.has("robust")) {
      config.topology.optimizer.projection.enabled = true;
      config.topology.optimizer.projection.robust = true;
      config.topology.length_scale_check = true;
    }
    if (cli.has("no-robust")) config.topology.optimizer.projection.robust = false;
    config.topology.optimizer.projection.validate();
    if (cli.has("overhang")) {
      config.topology.optimizer.overhang.direction =
          parse_build_direction(cli.value("overhang"), config.dim());
      config.topology.optimizer.overhang.filter = true;
      config.topology.overhang_check = true;
    }
    if (cli.has("no-overhang-filter")) config.topology.optimizer.overhang.filter = false;
    if (cli.has("no-vtk")) config.output.write_vtk = false;
    if (cli.has("no-csv")) config.output.write_csv = false;
    if (cli.has("tag")) config.name += "_" + cli.value("tag");

    if (!config.topology.enabled) {
      throw ConfigError(
          "the deck has topology.enabled = false; sparlab_topopt needs a topology "
          "section. Use sparlab_solve for a plain analysis");
    }
    if (config.nonlinear.enabled) {
      throw ConfigError(
          "the deck enables the non-linear analysis, but the optimiser designs for linear "
          "statics; run the non-linear analysis of a design with sparlab_solve");
    }
    {
      bool plastic = config.material().plasticity().enabled();
      for (const MaterialRegion& r : config.material_regions) {
        plastic = plastic || r.material.plasticity().enabled();
      }
      if (plastic) {
        log::warn("the deck's material has a yield stress, but the optimiser designs for "
                  "linear elasticity: the plasticity is ignored (a stress constraint can "
                  "keep the design below yield)");
      }
    }

    const std::string out_dir =
        cli.value("output", app::default_output_directory(config.name));

    TimingLedger timings;
    Timer wall;

    FemModel model = [&] {
      ScopedTimer t(timings, "model_build");
      return build_model(config);
    }();

    // The solid model must be well posed; the SIMP floor would otherwise mask a
    // missing boundary condition behind a very soft but non-singular matrix.
    require_well_posed(model);

    // Self-weight and thermal loads grow and shrink with the material the
    // optimiser places, and a second material needs its own interpolation;
    // treating either as a fixed load would optimise the wrong problem.
    if (!model.single_material()) {
      throw ConfigError("topology optimisation works on one material; remove "
                        "'material_regions' (use passive regions for fixed parts)");
    }
    for (const LoadCaseSpec& spec : model.load_case_specs()) {
      if (spec.has_body_loads() || spec.has_temperature()) {
        throw ConfigError("load case '" + spec.name +
                          "': gravity, body forces, rotation and temperature loads depend "
                          "on the design, and the optimiser does not yet interpolate "
                          "them; analyse such cases with sparlab_solve");
      }
    }

    Assembler assembler(model);
    const Scalar filter_radius = config.resolved_filter_radius(model.mesh());
    DensityFilter filter(model.mesh(), config.topology.filter_type, filter_radius);
    DesignDomain domain = build_design_domain(config, model);

    TopologyOptimizationResult result;
    {
      ScopedTimer t(timings, "optimization");
      TopologyOptimizer optimizer(model, assembler, filter, domain,
                                  config.topology.optimizer);
      result = optimizer.run();
    }

    // ---- reference analyses on the solid domain ---------------------------
    Scalar solid_compliance = 0.0;
    {
      ScopedTimer t(timings, "solid_reference");
      StaticAnalysisOptions opts = config.analysis;
      StaticAnalysis analysis(model, assembler, opts);
      const std::vector<StaticSolution> sols = analysis.solve_all();
      solid_compliance =
          StaticAnalysis::weighted_compliance(sols, model.normalised_weights());
    }

    // ---- stresses of the optimised design ---------------------------------
    std::vector<StressField> stresses;
    {
      ScopedTimer t(timings, "stress_recovery");
      for (const Vector& u : result.displacements) {
        stresses.push_back(
            recover_stresses(model, assembler, u, &result.stiffness_factors));
      }
    }

    // ---- interpret the density field as a solid body -----------------------
    std::unique_ptr<TopologyInterpretation> interpretation;
    {
      ScopedTimer t(timings, "interpretation");
      interpretation = std::make_unique<TopologyInterpretation>(
          interpret_density_as_solid(model.mesh(), result.physical_density,
                                     domain.element_volumes(),
                                     config.topology.optimizer.interpretation_threshold));
    }

    // ---- analysis of the interpreted solid --------------------------------
    // Thresholding the density field is the only way a SIMP result becomes a
    // structure, so what *that* structure delivers is reported on every run,
    // not only when modal analysis happens to be switched on. It is one extra
    // factorisation of a smaller model.
    //
    // The interpreted structure can legitimately be ill-posed - the largest
    // connected group need not reach a support, and on a coarse mesh with an
    // unresolved filter it often does not. That is a fact about the design,
    // not a reason to throw away a completed optimisation, so the failure is
    // recorded in the summary and warned about rather than made fatal.
    std::unique_ptr<SubModelBuild> sub_build;
    std::unique_ptr<Assembler> sub_assembler;
    json::Value interpreted = json::Value::make_object();
    {
      ScopedTimer t(timings, "interpreted_analysis");
      try {
        sub_build = std::make_unique<SubModelBuild>(build_solid_submodel(
            config, interpretation->sub.mesh, /*with_loads=*/true));
        require_well_posed(*sub_build->model);
        sub_assembler = std::make_unique<Assembler>(*sub_build->model);

        const SubModelBuild& sub = *sub_build;
        interpreted.set("loads_applied", json::Value::make_bool(sub.loads_applied));
        if (!sub.note.empty()) {
          interpreted.set("note", json::Value::make_string(sub.note));
        }
        interpreted.set("num_elements",
                        json::Value::make_number(sub.model->mesh().num_elements()));
        interpreted.set("num_nodes",
                        json::Value::make_number(sub.model->mesh().num_nodes()));
        interpreted.set(
            "mass_kg",
            json::Value::make_number(sub.model->domain_volume() *
                                     config.material().density()));
        if (sub.loads_applied) {
          StaticAnalysisOptions opts = config.analysis;
          StaticAnalysis sub_analysis(*sub.model, *sub_assembler, opts);
          const std::vector<StaticSolution> sols = sub_analysis.solve_all();
          const Scalar c = StaticAnalysis::weighted_compliance(
              sols, sub.model->normalised_weights());
          interpreted.set("weighted_compliance_J", json::Value::make_number(c));
          // The gap between the penalised relaxed model and the structure a
          // threshold actually produces. It can fall either side of 1:
          // thresholding promotes every element above the threshold to full
          // material, which the SIMP penalty had been crediting with only
          // rho^p of its stiffness (stiffer, ratio < 1), and deletes
          // everything below it (softer, ratio > 1). Which effect wins is a
          // property of the design, so it is measured rather than assumed.
          interpreted.set("compliance_vs_simp_ratio",
                          json::Value::make_number(
                              result.compliance > 0.0 ? c / result.compliance : 0.0));
          Scalar max_vm = 0.0;
          json::Value per_case = json::Value::make_array();
          for (const StaticSolution& s : sols) {
            const StressField f =
                recover_stresses(*sub.model, *sub_assembler, s.displacement);
            max_vm = std::max(max_vm, f.element_von_mises.maxCoeff());
            per_case.push_back(json::Value::make_number(f.element_von_mises.maxCoeff()));
          }
          interpreted.set("max_von_mises_Pa", json::Value::make_number(max_vm));
          interpreted.set("max_von_mises_per_load_case_Pa", per_case);
          if (result.stress_constrained) {
            // The check that matters: the thresholded structure with full
            // material, against the limit the optimiser was given.
            const Scalar limit = config.topology.optimizer.stress.limit;
            interpreted.set("stress_limit_Pa", json::Value::make_number(limit));
            interpreted.set("max_von_mises_over_limit",
                            json::Value::make_number(max_vm / limit));
            interpreted.set("meets_stress_limit", json::Value::make_bool(max_vm <= limit));
          }
        }
      } catch (const std::exception& error) {
        sub_build.reset();
        sub_assembler.reset();
        interpreted.set("analysis_failed", json::Value::make_bool(true));
        interpreted.set("reason", json::Value::make_string(error.what()));
        log::warn("the density field thresholded at ",
                  config.topology.optimizer.interpretation_threshold,
                  " does not form an analysable structure: ", error.what(),
                  ". The optimisation result itself is unaffected and is still "
                  "reported; treat the interpreted geometry as invalid");
      }
    }

    // ---- modal comparison --------------------------------------------------
    std::unique_ptr<ModalResult> modal_solid;
    std::unique_ptr<ModalResult> modal_thin;
    std::unique_ptr<ModalResult> modal_topology;

    if (config.modal.enabled) {
      ScopedTimer t(timings, "modal_analysis");
      modal_solid = std::make_unique<ModalResult>(
          solve_modal(model, assembler, config.modal.options));

      if (config.modal.compare_mass_matched_baseline && config.dim() == 2) {
        // Equal-mass uniform plate: thickness scaled by the achieved volume
        // fraction. K and M both scale with thickness, so the frequencies must
        // match the full solid; the run checks that instead of assuming it.
        // The argument is a 2-D one, so a solid mesh skips it (see below).
        Configuration thin = config;
        thin.thickness = config.thickness * result.volume_fraction;
        FemModel thin_model = build_model(thin);
        Assembler thin_assembler(thin_model);
        modal_thin = std::make_unique<ModalResult>(
            solve_modal(thin_model, thin_assembler, config.modal.options));
      }

      if (config.modal.analyse_optimised_topology) {
        if (sub_build != nullptr) {
          modal_topology = std::make_unique<ModalResult>(
              solve_modal(*sub_build->model, *sub_assembler, config.modal.options));
        } else {
          log::warn("skipping the modal analysis of the optimised topology: the "
                    "interpreted structure could not be analysed (see the "
                    "interpreted_solid_analysis block of summary.json)");
        }
      }
    }

    // ---- manufacturability checks -------------------------------------------
    // What the final design delivers on the two process limits this code
    // models: unsupported overhangs for a build direction, and the minimum
    // member and gap size. Nothing else about manufacturability is checked.
    json::Value manufacturing = json::Value::make_object();
    std::unique_ptr<OverhangReport> overhang_report;
    std::unique_ptr<LengthScaleScan> length_report;
    {
      const Scalar threshold = config.topology.optimizer.interpretation_threshold;
      if (config.topology.overhang_check) {
        const OverhangFilter stencil(model.mesh(), config.topology.optimizer.overhang);
        overhang_report = std::make_unique<OverhangReport>(check_overhang(
            stencil, result.physical_density, domain.element_volumes(), threshold));
        manufacturing.set("overhang", overhang_json(*overhang_report,
                                                    result.overhang_filtered));
      }
      if (config.topology.length_scale_check) {
        const Scalar cell = model.mesh().mean_element_size();
        length_report = std::make_unique<LengthScaleScan>(scan_length_scale(
            model.mesh(), result.physical_density, domain.element_volumes(), threshold,
            0.5 * cell, config.topology.length_scale_max_radius_elements * cell,
            config.topology.length_scale_tolerance));
        manufacturing.set("length_scale", length_scale_json(*length_report));
      }
    }

    // ---- buckling check ------------------------------------------------------
    // The load factors of the full solid domain and of the interpreted
    // structure - the part the STL describes - for the same load cases. The
    // SIMP design itself is not checked here: its void material would need
    // the pseudo-mode treatment of the buckling constraint, and it is not
    // what would be built.
    std::vector<BucklingResult> buckling_solid;
    std::vector<BucklingResult> buckling_topology;
    json::Value buckling_check = json::Value::make_object();
    if (config.buckling.enabled) {
      ScopedTimer t(timings, "buckling_check");
      const std::vector<std::size_t> cases = config.buckling_load_cases();
      buckling_solid = check_buckling(model, assembler, config, cases);
      buckling_check.set("full_solid", buckling_json(buckling_solid, config.buckling.options,
                                                     "the full solid design domain"));
      if (config.buckling.analyse_optimised_topology) {
        if (sub_build != nullptr && sub_build->loads_applied) {
          buckling_topology = check_buckling(*sub_build->model, *sub_assembler, config, cases);
          buckling_check.set(
              "interpreted_structure",
              buckling_json(buckling_topology, config.buckling.options,
                            "the density field thresholded at "
                            "solid_interpretation.threshold, largest face-connected "
                            "group, full material: the exported structure_after part"));
        } else {
          const std::string reason =
              sub_build == nullptr
                  ? "the interpreted structure could not be analysed (see "
                    "interpreted_solid_analysis)"
                  : "the load cases could not be applied to the interpreted structure";
          buckling_check.set("interpreted_structure_skipped",
                             json::Value::make_string(reason));
          log::warn("skipping the buckling check of the interpreted structure: ", reason);
        }
      }
    }

    // ---- output ------------------------------------------------------------
    ResultWriter writer(out_dir, config);
    {
      ScopedTimer t(timings, "output");
      writer.write_config();
      writer.write_mesh(model);
      writer.write_history(result);
      if (config.output.write_csv) {
        writer.write_density(model.mesh(), domain, result);
        for (std::size_t l = 0; l < result.displacements.size(); ++l) {
          const std::string& name = config.load_cases[l].name;
          writer.write_displacement(model.mesh(), name, result.displacements[l]);
          writer.write_stress(model.mesh(), name, stresses[l],
                              &result.physical_density);
        }
      }
      if (config.output.write_density_history) {
        writer.write_density_history(result);
      }
      if (config.output.write_vtk) {
        for (std::size_t l = 0; l < result.displacements.size(); ++l) {
          writer.write_static_vtk(model.mesh(), config.load_cases[l].name,
                                  result.displacements[l], stresses[l],
                                  &result.physical_density,
                                  &result.stiffness_factors);
        }
      }
      if (modal_solid) writer.write_modal(model.mesh(), *modal_solid, "solid");
      if (modal_topology) {
        writer.write_modal(interpretation->sub.mesh, *modal_topology, "topology");
      }
      if (!buckling_solid.empty()) writer.write_buckling(model.mesh(), buckling_solid, "solid");
      if (!buckling_topology.empty()) {
        writer.write_buckling(interpretation->sub.mesh, buckling_topology, "topology");
      }
    }

    // ---- before / after geometry ------------------------------------------
    // "Before" is the design domain with the uniform starting density;
    // "after" is the thresholded, largest-group structure with the physical
    // density of each retained cell. Both are solids in the STL sense only
    // through the interpretation recorded next to them.
    json::Value geometry = json::Value::make_object();
    {
      ScopedTimer t(timings, "geometry_export");
      geometry.set("before", writer.write_geometry(model.mesh(), domain.initial_design(),
                                                   model.thickness(), "structure_before",
                                                   "design domain before optimisation"));
      Vector retained_density(interpretation->sub.mesh.num_elements());
      for (std::size_t i = 0; i < interpretation->sub.element_map.size(); ++i) {
        retained_density(static_cast<Eigen::Index>(i)) =
            result.physical_density(interpretation->sub.element_map[i]);
      }
      geometry.set("after", writer.write_geometry(interpretation->sub.mesh, retained_density,
                                                  model.thickness(), "structure_after",
                                                  "interpreted structure after optimisation"));
      geometry.set("note", json::Value::make_string(
                               "structure_after is the density field thresholded at "
                               "solid_interpretation.threshold with only the largest "
                               "face-connected group kept, written as the cells' outer "
                               "surface" +
                               std::string(model.dim() == 2
                                               ? " extruded by model.thickness; "
                                               : "; ") +
                               "it is an interpretation of a SIMP result, not a "
                               "checked design. structure_before is the full design "
                               "domain with the uniform starting density."));
    }

    timings.add("total", wall.elapsed_seconds());

    json::Value summary = make_topology_summary(
        config, model, domain, filter, result, interpretation.get(), modal_solid.get(),
        modal_topology.get(), modal_thin.get(), timings);

    // Mass-stiffness comparison at equal mass.
    {
      json::Value comparison = json::Value::make_object();
      comparison.set("full_solid_compliance_J",
                     json::Value::make_number(solid_compliance));
      const Scalar equal_mass = result.volume_fraction > 0.0
                                    ? solid_compliance / result.volume_fraction
                                    : 0.0;
      if (config.dim() == 2) {
        comparison.set("equal_mass_uniform_plate_compliance_J",
                       json::Value::make_number(equal_mass));
      }
      comparison.set("optimised_compliance_J",
                     json::Value::make_number(result.compliance));
      if (config.dim() == 2) {
        comparison.set("stiffness_gain_over_equal_mass_plate",
                       json::Value::make_number(
                           result.compliance > 0.0 ? equal_mass / result.compliance : 0.0));
      }
      comparison.set("compliance_penalty_vs_full_solid",
                     json::Value::make_number(solid_compliance > 0.0
                                                  ? result.compliance / solid_compliance
                                                  : 0.0));
      comparison.set(
          "note",
          json::Value::make_string(
              config.dim() == 2
                  ? "In this 2-D idealisation K and M both scale linearly with "
                    "thickness, so an equal-mass uniform plate has the same natural "
                    "frequencies as the full solid domain and a compliance of C_solid / "
                    "volume_fraction. That plate is the equal-mass baseline for the "
                    "optimised design."
                  : "A 3-D solid has no thickness to scale, so the equal-mass uniform "
                    "plate of the 2-D runs does not exist here. The only reference "
                    "reported is the full solid domain; the interpreted structure "
                    "(interpreted_solid_analysis) is the design's own re-solve."));
      if (modal_thin && modal_solid && modal_solid->frequencies_hz.size() > 0) {
        const Eigen::Index n =
            std::min(modal_thin->frequencies_hz.size(), modal_solid->frequencies_hz.size());
        Scalar max_rel = 0.0;
        for (Eigen::Index i = 0; i < n; ++i) {
          const Scalar ref = std::max(std::abs(modal_solid->frequencies_hz(i)), 1e-30);
          max_rel = std::max(max_rel,
                             std::abs(modal_thin->frequencies_hz(i) -
                                      modal_solid->frequencies_hz(i)) /
                                 ref);
        }
        comparison.set("thickness_scaling_frequency_max_relative_difference",
                       json::Value::make_number(max_rel));
      }
      summary.set("mass_stiffness_comparison", comparison);
    }
    if (!interpreted.members().empty()) {
      summary.set("interpreted_solid_analysis", interpreted);
    }
    summary.set("geometry_export", geometry);
    if (!buckling_check.members().empty()) summary.set("buckling_check", buckling_check);
    if (!manufacturing.members().empty()) summary.set("manufacturing_checks", manufacturing);
    writer.write_json("summary.json", summary);

    // ---- console report ----------------------------------------------------
    std::cout << "case: " << config.name << "\n";
    std::cout << "  design:      " << domain.describe() << "\n";
    std::cout << "  filter:      " << to_string(filter.type()) << ", radius "
              << app::format(filter.radius()) << " m ("
              << app::format(filter.average_support()) << " elements support)\n";
    std::cout << "  converged:   " << (result.converged ? "yes" : "NO") << " after "
              << result.iterations << " iterations (" << result.linear_solves
              << " linear solves)\n";
    std::cout << "  compliance:  " << app::format(result.compliance)
              << " J (full solid " << app::format(solid_compliance) << " J";
    if (config.dim() == 2) {
      std::cout << ", equal-mass uniform plate "
                << app::format(result.volume_fraction > 0.0
                                   ? solid_compliance / result.volume_fraction
                                   : 0.0)
                << " J";
    }
    std::cout << ")\n";
    std::cout << "  volume:      " << app::format(result.volume) << " m^3, fraction "
              << app::format(result.volume_fraction) << " (target "
              << app::format(domain.volume_fraction()) << ", relative violation "
              << app::format(result.volume_constraint_violation) << ")\n";
    std::cout << "  grey level:  " << app::format(result.gray_level);
    if (result.projected) {
      std::cout << " (projected at beta " << app::format(result.final_beta)
                << "; before projection " << app::format(gray_level(result.filtered_density))
                << ")";
    }
    std::cout << "\n";
    std::cout << "  solver:      " << result.linear_solver;
    if (result.linear_iterations > 0) {
      std::cout << ", " << result.linear_iterations << " CG iterations over "
                << result.linear_solves << " solves";
      if (result.has_amg_stats) {
        std::cout << " (" << result.amg_stats.levels.size() << " levels, operator complexity "
                  << app::format(result.amg_stats.operator_complexity, 4) << ")";
      }
    }
    std::cout << "\n";
    std::cout << "  method:      " << to_string(result.method);
    if (result.method == OptimizerMethod::MMA) {
      std::cout << " (largest constraint value "
                << app::format(result.constraint_violation) << ", "
                << (result.feasible ? "feasible" : "INFEASIBLE") << ")";
    }
    std::cout << "\n";
    if (result.stress_constrained) {
      std::cout << "  stress:      max relaxed ratio " << app::format(result.max_stress_ratio)
                << " of the " << app::format(config.topology.optimizer.stress.limit)
                << " Pa limit";
      const json::Value* over = interpreted.find("max_von_mises_over_limit");
      if (over != nullptr) {
        std::cout << "; interpreted structure re-solve " << app::format(over->number_value())
                  << " of the limit";
      }
      std::cout << "\n";
    }
    if (result.robust || result.erosion_checked) {
      const RobustRecord& rr = result.robust_record;
      std::cout << (result.robust ? "  robust:      compliance eroded "
                                  : "  erosion:     compliance eroded ")
                << app::format(rr.compliance_eroded)
                << " J, blueprint " << app::format(rr.compliance_intermediate)
                << " J, dilated " << app::format(rr.compliance_dilated)
                << " J; volume fractions " << app::format(rr.volume_fraction_eroded) << " / "
                << app::format(rr.volume_fraction_intermediate) << " / "
                << app::format(rr.volume_fraction_dilated) << "\n";
    }
    if (overhang_report) {
      std::cout << "  overhang:    build " << overhang_report->build_direction << ", "
                << overhang_report->unsupported_elements << " of "
                << overhang_report->solid_elements << " solid elements unsupported ("
                << app::format(100.0 * overhang_report->unsupported_fraction)
                << " % of the solid volume)"
                << (result.overhang_filtered ? ", filter on" : ", filter off") << "\n";
    }
    if (length_report) {
      const Scalar cell = model.mesh().mean_element_size();
      std::cout << "  length scale: members at least ~" << app::format(length_report->solid_min_size)
                << " m (" << app::format(length_report->solid_min_size / cell)
                << " cells), gaps at least ~" << app::format(length_report->void_min_size)
                << " m (" << app::format(length_report->void_min_size / cell) << " cells)"
                << (length_report->solid_bound_reached_cap ? ", solid bound at the scan cap" : "")
                << "\n";
    }
    std::cout << "  interpreted: threshold "
              << app::format(interpretation->threshold) << " keeps "
              << interpretation->elements_retained << " of "
              << model.mesh().num_elements() << " elements in "
              << interpretation->components_above_threshold << " group(s), discarding "
              << app::format(interpretation->volume_discarded_as_islands)
              << " m^3 as islands\n";
    if (modal_solid) {
      std::cout << "  f1 solid:    " << app::format(modal_solid->frequencies_hz(0))
                << " Hz (mass " << app::format(modal_solid->total_mass) << " kg)\n";
    }
    if (modal_topology) {
      std::cout << "  f1 topology: " << app::format(modal_topology->frequencies_hz(0))
                << " Hz (mass " << app::format(modal_topology->total_mass) << " kg)\n";
    }
    if (result.buckling_constrained) {
      std::cout << "  buckling:    SIMP design lambda_1 " << app::format(result.min_load_factor)
                << " (required " << app::format(config.topology.optimizer.buckling.min_load_factor)
                << ")\n";
    }
    for (std::size_t i = 0; i < buckling_solid.size(); ++i) {
      const auto first = [](const BucklingResult& b) {
        return b.no_positive_load_factor ? std::string("none (stiffening)")
                                         : app::format(b.load_factors(0));
      };
      std::cout << "  lambda_1 '" << buckling_solid[i].load_case << "': full solid "
                << first(buckling_solid[i]);
      if (i < buckling_topology.size()) {
        std::cout << ", interpreted structure " << first(buckling_topology[i]);
      }
      std::cout << "\n";
    }
    std::cout << "  geometry:    structure_before.{vtk,stl} and structure_after.{vtk,stl} ("
              << geometry.find("after")->find("num_triangles")->number_value()
              << " triangles, "
              << app::format(geometry.find("after")->find("enclosed_volume_m3")->number_value())
              << " m^3 enclosed)\n";
    std::cout << "  runtime:     " << app::format(timings.get("total")) << " s\n";
    std::cout << "  results:     " << out_dir << "\n";
    for (const std::string& w : result.warnings) {
      std::cout << "  warning:     " << w << "\n";
    }
    return 0;
  });
}
