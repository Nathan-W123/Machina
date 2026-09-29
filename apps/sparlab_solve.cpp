/// \file sparlab_solve.cpp
/// \brief Linear static (and optional modal, buckling and geometrically
///        non-linear) analysis of one configuration.
///
/// Solves every load case of the deck, recovers stresses and reactions, checks
/// global equilibrium, optionally runs a modal analysis, a linear buckling
/// check of the load cases (reusing the static factorisation) and a
/// large-deflection analysis of the selected load cases, and writes the full
/// result set to the output directory.

#include "AppSupport.hpp"

#include "sparlab/core/Timer.hpp"
#include "sparlab/fem/Assembler.hpp"
#include "sparlab/fem/Buckling.hpp"
#include "sparlab/fem/ModalAnalysis.hpp"
#include "sparlab/fem/ModelDiagnostics.hpp"
#include "sparlab/fem/NonlinearStatic.hpp"
#include "sparlab/fem/StaticAnalysis.hpp"
#include "sparlab/fem/StressRecovery.hpp"
#include "sparlab/io/CalculixWriter.hpp"
#include "sparlab/io/Config.hpp"
#include "sparlab/io/ResultWriter.hpp"

#include <iostream>

using namespace sparlab;

int main(int argc, char** argv) {
  return app::run_guarded([&]() -> int {
    const std::vector<std::string> known = {"config",   "output",        "verbosity",
                                            "strict-config", "modes",   "no-vtk",
                                            "no-csv",   "export-calculix", "solver",
                                            "buckling", "nonlinear", "help"};
    app::CommandLine cli(argc, argv, known);
    if (cli.has("help") || argc == 1) {
      return app::print_usage(
          "sparlab_solve", "--config <deck.json> [--output <dir>]",
          {{"--config <file>", "input deck describing mesh, material, BCs, loads"},
           {"--output <dir>", "output directory (default results/<case>)"},
           {"--modes <n>", "override modal.num_modes and enable modal analysis"},
           {"--buckling <n>", "enable the linear buckling check with n modes per load "
                              "case"},
           {"--nonlinear", "enable the geometrically non-linear analysis (the deck's "
                           "'nonlinear' settings, or their defaults)"},
           {"--solver <type>", "override solver.linear.type (simplicial_ldlt, amg_cg, "
                               "auto, ...)"},
           {"--no-vtk", "skip VTK output"},
           {"--no-csv", "skip per-node/per-element CSV output"},
           {"--export-calculix", "also write one CalculiX .inp per load case, a "
                                 "heat-transfer .inp for a conducted temperature and an "
                                 "NLGEOM .inp for a non-linear case, into the output "
                                 "directory (cross-validation)"},
           {"--strict-config", "treat unknown configuration keys as errors"},
           {"--verbosity <lvl>", "trace|debug|info|warn|error|silent"},
           {"--help", "show this message"}});
    }
    app::apply_verbosity(cli);

    Configuration config =
        load_configuration(cli.require("config"), cli.has("strict-config"));
    if (cli.has("modes")) {
      config.modal.enabled = true;
      config.modal.options.num_modes = cli.integer("modes", 6);
    }
    if (cli.has("no-vtk")) config.output.write_vtk = false;
    if (cli.has("no-csv")) config.output.write_csv = false;
    if (cli.has("buckling")) {
      config.buckling.enabled = true;
      config.buckling.options.num_modes = cli.integer("buckling", 4);
    }
    if (cli.has("nonlinear")) {
      config.nonlinear.enabled = true;
      (void)config.nonlinear_load_cases();  // validates the deck's names
    }
    if (cli.has("solver")) {
      config.analysis.linear.type = parse_linear_solver_type(cli.value("solver"));
      config.modal.options.linear.type = config.analysis.linear.type;
      config.buckling.options.linear.type = config.analysis.linear.type;
    }

    const std::string out_dir =
        cli.value("output", app::default_output_directory(config.name));

    TimingLedger timings;
    Timer wall;

    FemModel model = [&] {
      ScopedTimer t(timings, "model_build");
      return build_model(config);
    }();

    const ModelDiagnostics diagnostics = diagnose_model(model);
    if (!diagnostics.well_posed()) {
      // Report every problem, then fail: an ill-posed static model cannot
      // produce a meaningful answer.
      require_well_posed(model);
    }

    Assembler assembler(model);
    log::info("assembler: element matrix cache ",
              (assembler.uses_element_cache() ? "enabled (uniform structured mesh)"
                                              : "disabled (general mesh)"));

    std::vector<StaticSolution> solutions;
    StaticAnalysis analysis(model, assembler, config.analysis);
    {
      ScopedTimer t(timings, "static_solve");
      solutions = analysis.solve_all();
    }

    // Buckling of each checked load case, with the static factorisation.
    std::vector<BucklingResult> buckling;
    if (config.buckling.enabled) {
      ScopedTimer t(timings, "buckling_analysis");
      const DofManager& dofs = model.dofs();
      const FreeSolve solve = [&](const Vector& b) {
        return dofs.restrict_to_free(analysis.solve_homogeneous(dofs.expand(b)));
      };
      for (std::size_t l : config.buckling_load_cases()) {
        const Vector& temperature = model.load_case_data(l).temperature;
        const SparseMatrix k_g = assemble_geometric_stiffness(
            model, assembler, solutions[l].displacement, nullptr,
            temperature.size() > 0 ? &temperature : nullptr);
        BucklingResult result = solve_buckling(model, assembler, analysis.stiffness(), k_g,
                                               config.buckling.options, solve);
        result.load_case = solutions[l].load_case_name;
        result.linear_solver = analysis.solver().name();
        buckling.push_back(std::move(result));
      }
    }

    // Geometrically non-linear analysis of each selected load case.
    std::vector<NonlinearResult> nonlinear;
    if (config.nonlinear.enabled) {
      ScopedTimer t(timings, "nonlinear_analysis");
      NonlinearStaticAnalysis nl(model, assembler, config.nonlinear.options);
      for (std::size_t l : config.nonlinear_load_cases()) nonlinear.push_back(nl.solve(l));
    }

    std::vector<StressField> stresses;
    {
      ScopedTimer t(timings, "stress_recovery");
      stresses.reserve(solutions.size());
      for (std::size_t l = 0; l < solutions.size(); ++l) {
        const Vector& temperature = model.load_case_data(l).temperature;
        stresses.push_back(recover_stresses(model, assembler, solutions[l].displacement,
                                            nullptr,
                                            temperature.size() > 0 ? &temperature : nullptr));
      }
    }

    std::unique_ptr<ModalResult> modal;
    if (config.modal.enabled) {
      ScopedTimer t(timings, "modal_analysis");
      modal = std::make_unique<ModalResult>(
          solve_modal(model, assembler, config.modal.options));
    }

    ResultWriter writer(out_dir, config);
    {
      ScopedTimer t(timings, "output");
      writer.write_config();
      writer.write_mesh(model);
      for (std::size_t l = 0; l < solutions.size(); ++l) {
        const std::string& name = solutions[l].load_case_name;
        if (config.output.write_csv) {
          writer.write_displacement(model.mesh(), name, solutions[l].displacement);
          writer.write_stress(model.mesh(), name, stresses[l]);
          writer.write_reactions(model.mesh(), model.dofs(), name,
                                 solutions[l].reactions);
        }
        const Vector& temperature = model.load_case_data(l).temperature;
        if (config.output.write_csv && temperature.size() > 0) {
          writer.write_temperature(model.mesh(), name, temperature);
        }
        if (config.output.write_vtk) {
          writer.write_static_vtk(model.mesh(), name, solutions[l].displacement,
                                  stresses[l], nullptr, nullptr,
                                  temperature.size() > 0 ? &temperature : nullptr);
        }
      }
      if (modal) writer.write_modal(model.mesh(), *modal);
      if (!buckling.empty()) writer.write_buckling(model.mesh(), buckling);
      for (const NonlinearResult& r : nonlinear) writer.write_nonlinear(model, r);
      if (cli.has("export-calculix")) {
        // The non-linear cases go out as NLGEOM decks too, when CalculiX has
        // the same material law.
        CalculixNonlinearExport nlgeom;
        const CalculixNonlinearExport* nonlinear_export = nullptr;
        if (config.nonlinear.enabled) {
          if (config.nonlinear.options.law == HyperelasticModel::SaintVenantKirchhoff) {
            nlgeom.load_cases = config.nonlinear_load_cases();
            nlgeom.increments = config.nonlinear.options.steps;
            nlgeom.follower_pressure = config.nonlinear.options.follower_pressure;
            nonlinear_export = &nlgeom;
          } else {
            log::warn("the non-linear cases are not exported to CalculiX: its NEO HOOKE is a "
                      "different strain energy from SparLab's neo-Hookean law");
          }
        }
        const std::vector<std::string> decks = write_calculix_decks(
            model, writer.file("calculix"), config.name, nonlinear_export);
        for (const std::string& deck : decks) log::info("wrote CalculiX deck ", deck);
      }
    }

    timings.add("total", wall.elapsed_seconds());
    json::Value summary = make_static_summary(config, model, diagnostics, solutions,
                                              stresses, modal.get(), timings);
    if (!buckling.empty()) {
      summary.set("buckling", buckling_json(buckling, config.buckling.options,
                                            "the model as meshed"));
    }
    if (!nonlinear.empty()) {
      summary.set("nonlinear",
                  nonlinear_json(nonlinear, config.nonlinear.options, model, solutions));
    }
    writer.write_json("summary.json", summary);

    std::cout << "case: " << config.name << "\n";
    std::cout << "  mesh:        " << model.mesh().num_elements() << " elements, "
              << model.dofs().num_dofs() << " DOFs (" << model.dofs().num_free()
              << " free)\n";
    if (!solutions.empty()) {
      std::cout << "  solver:      " << solutions.front().solver_name;
      int iterations = 0;
      for (const StaticSolution& s : solutions) iterations += s.solver_iterations;
      if (iterations > 0) {
        std::cout << ", " << iterations << " CG iterations over " << solutions.size()
                  << " load case(s)";
      }
      std::cout << "\n";
    }
    for (std::size_t l = 0; l < solutions.size(); ++l) {
      const StaticSolution& s = solutions[l];
      std::cout << "  load case '" << s.load_case_name << "': compliance "
                << app::format(s.compliance) << " J, strain energy "
                << app::format(s.strain_energy) << " J, max |u| "
                << app::format(s.max_displacement_magnitude) << " m, max von Mises "
                << app::format(stresses[l].element_von_mises.maxCoeff())
                << " Pa\n";
      const auto vec = [&](const Vector3& v) {
        std::string text = "(" + app::format(v.x()) + ", " + app::format(v.y());
        if (model.dim() == 3) text += ", " + app::format(v.z());
        return text + ")";
      };
      std::cout << "      equilibrium: applied " << vec(s.equilibrium.applied_force)
                << " N, reactions " << vec(s.equilibrium.reaction_force)
                << " N, relative error "
                << app::format(s.equilibrium.relative_force_error) << "\n";
    }
    if (modal) {
      std::cout << "  modal: ";
      for (Eigen::Index i = 0; i < modal->frequencies_hz.size(); ++i) {
        std::cout << (i ? ", " : "") << "f" << i + 1 << " = "
                  << app::format(modal->frequencies_hz(i)) << " Hz";
      }
      std::cout << "\n";
      std::cout << "      total mass " << app::format(modal->total_mass) << " kg\n";
    }
    for (const BucklingResult& b : buckling) {
      std::cout << "  buckling '" << b.load_case << "': ";
      if (b.no_positive_load_factor) {
        std::cout << "no positive load factor (the load case only stiffens the model)";
      } else {
        for (Eigen::Index i = 0; i < b.load_factors.size(); ++i) {
          std::cout << (i ? ", " : "") << "lambda" << i + 1 << " = "
                    << app::format(b.load_factors(i));
        }
      }
      std::cout << "\n";
    }
    for (const NonlinearResult& r : nonlinear) {
      std::cout << "  non-linear '" << r.load_case_name << "' (" << r.method << ", " << r.law
                << "): " << (r.completed ? "completed" : "STOPPED") << " at lambda = "
                << app::format(r.load_factor) << " after " << r.steps.size() << " step(s), "
                << r.total_iterations << " iteration(s)";
      if (!r.steps.empty()) {
        std::cout << ", max |u| " << app::format(r.steps.back().max_displacement) << " m";
      }
      std::cout << "\n";
      if (!r.completed) std::cout << "      " << r.termination << "\n";
    }
    std::cout << "  runtime:     " << app::format(timings.get("total")) << " s\n";
    std::cout << "  results:     " << out_dir << "\n";
    return 0;
  });
}
