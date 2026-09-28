/// \file test_io.cpp
/// \brief JSON parsing/serialisation, configuration validation and writers.
#include "TestSupport.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Timer.hpp"
#include "sparlab/io/CalculixWriter.hpp"
#include "sparlab/io/Config.hpp"
#include "sparlab/io/CsvWriter.hpp"
#include "sparlab/io/Json.hpp"
#include "sparlab/io/ResultWriter.hpp"
#include "sparlab/io/VtkWriter.hpp"
#include "sparlab/topopt/DensityFilter.hpp"
#include "sparlab/topopt/DesignDomain.hpp"
#include "sparlab/topopt/TopologyOptimizer.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace sparlab;
using namespace sparlab::testing;
using Catch::Approx;

namespace {

/// A minimal but complete deck used by the configuration tests.
std::string minimal_deck() {
  return R"({
    "name": "unit_case",
    "description": "a deck used by the test suite",
    "mesh": { "type": "structured_quad", "nx": 8, "ny": 4, "lx": 0.8, "ly": 0.4 },
    "material": { "name": "al", "youngs_modulus": 70e9, "poisson_ratio": 0.33,
                  "density": 2700 },
    "model": { "thickness": 0.01, "stress_state": "plane_stress" },
    "boundary_conditions": [
      { "name": "root", "fix": ["x", "y"], "region": { "box": { "xmax": 0.0 } } }
    ],
    "load_cases": [
      { "name": "tip", "weight": 1.0,
        "point_loads": [
          { "name": "tip_edge", "force": [0.0, -1000.0], "distribution": "total",
            "region": { "box": { "xmin": 0.8 } } }
        ]
      }
    ]
  })";
}

std::string temp_path(const std::string& name) {
  return path_join("results/_test_tmp", name);
}

}  // namespace

TEST_CASE("JSON parser handles the full value grammar", "[io][json]") {
  const std::string text = R"({
    // a line comment
    "string": "with \"escapes\" and \n newline and é",
    /* a block comment */
    "number": -1.25e-3,
    "integer": 42,
    "yes": true,
    "no": false,
    "nothing": null,
    "array": [1, 2, [3, 4], {"nested": "value"}],
    "object": { "a": 1, "b": 2 },
    "trailing": [1, 2, 3,],
  })";
  const json::Value doc = json::parse(text, "test");
  REQUIRE(doc.is_object());
  REQUIRE(doc.find("string")->string_value() ==
          std::string("with \"escapes\" and \n newline and \xc3\xa9"));
  REQUIRE(doc.find("number")->number_value() == Approx(-1.25e-3));
  REQUIRE(doc.find("integer")->number_value() == Approx(42.0));
  REQUIRE(doc.find("yes")->bool_value());
  REQUIRE_FALSE(doc.find("no")->bool_value());
  REQUIRE(doc.find("nothing")->is_null());
  REQUIRE(doc.find("array")->array_items().size() == 4);
  REQUIRE(doc.find("array")->array_items()[2].array_items().size() == 2);
  REQUIRE(doc.find("object")->find("b")->number_value() == Approx(2.0));
  REQUIRE(doc.find("trailing")->array_items().size() == 3);
  REQUIRE(doc.find("absent") == nullptr);
}

TEST_CASE("JSON syntax errors report the location", "[io][json][diagnostics]") {
  const auto message = [](const std::string& text) {
    try {
      json::parse(text, "deck.json");
    } catch (const ConfigError& e) {
      return std::string(e.what());
    }
    return std::string("no error");
  };

  REQUIRE(message("{\"a\": }").find("deck.json:1:") != std::string::npos);
  REQUIRE_THROWS_AS(json::parse("{\"a\": 1", "x"), ConfigError);
  REQUIRE_THROWS_AS(json::parse("{a: 1}", "x"), ConfigError);
  REQUIRE_THROWS_AS(json::parse("[1, 2", "x"), ConfigError);
  REQUIRE_THROWS_AS(json::parse("{\"a\": tru}", "x"), ConfigError);
  REQUIRE_THROWS_AS(json::parse("{\"a\": \"unterminated}", "x"), ConfigError);
  REQUIRE_THROWS_AS(json::parse("{\"a\": 1} extra", "x"), ConfigError);
  REQUIRE_THROWS_AS(json::parse("{\"a\": 1, \"a\": 2}", "x"), ConfigError);
  REQUIRE_THROWS_AS(json::parse("{\"a\": /* unterminated", "x"), ConfigError);
  REQUIRE_THROWS_AS(json::parse("", "x"), ConfigError);
  REQUIRE_THROWS_AS(json::parse_file("/definitely/not/here.json"), IoError);

  // A multi-line document reports the right line number.
  const std::string multi = "{\n  \"a\": 1,\n  \"b\": @\n}";
  REQUIRE(message(multi).find(":3:") != std::string::npos);
}

TEST_CASE("JSON serialisation round-trips", "[io][json]") {
  json::Value doc = json::Value::make_object();
  doc.set("name", json::Value::make_string("round\"trip\n"));
  doc.set("value", json::Value::make_number(1.0 / 3.0));
  doc.set("flag", json::Value::make_bool(true));
  doc.set("empty_object", json::Value::make_object());
  doc.set("empty_array", json::Value::make_array());
  Vector v(3);
  v << 1.5, -2.5, 1.0e-17;
  doc.set("vector", json::array_of(v));
  doc.set("names", json::array_of(std::vector<std::string>{"a", "b"}));

  const std::string dumped = json::dump(doc, 2);
  const json::Value again = json::parse(dumped, "roundtrip");
  REQUIRE(again.find("name")->string_value() == doc.find("name")->string_value());
  REQUIRE(again.find("value")->number_value() == doc.find("value")->number_value());
  REQUIRE(again.find("vector")->array_items()[2].number_value() == Approx(1.0e-17));
  REQUIRE(again.find("empty_array")->array_items().empty());

  // Compact form parses too.
  REQUIRE_NOTHROW(json::parse(json::dump(doc, 0), "compact"));

  // Non-finite numbers become null rather than invalid JSON.
  json::Value bad = json::Value::make_object();
  bad.set("inf", json::Value::make_number(std::numeric_limits<double>::infinity()));
  REQUIRE(json::dump(bad, 0).find("null") != std::string::npos);
  REQUIRE_NOTHROW(json::parse(json::dump(bad, 0), "inf"));

  // set() overwrites in place and preserves ordering.
  doc.set("flag", json::Value::make_bool(false));
  REQUIRE_FALSE(doc.find("flag")->bool_value());
}

TEST_CASE("ConfigNode reports paths, types and unread keys", "[io][config]") {
  const json::Value doc = json::parse(R"({
    "a": { "b": 3.5, "c": "text", "d": true, "v": [1.0, 2.0] },
    "list": [ {"x": 1}, {"x": 2} ],
    "typo_key": 1
  })",
                                      "t");
  const ConfigNode root(doc);

  REQUIRE(root.child("a").child("b").number() == Approx(3.5));
  REQUIRE(root.child("a").string_or("c", "") == "text");
  REQUIRE(root.child("a").boolean_or("d", false));
  REQUIRE(root.child("a").vector2_or("v", Vector2::Zero()).isApprox(Vector2(1.0, 2.0)));
  REQUIRE(root.child("a").number_or("missing", 7.0) == Approx(7.0));
  REQUIRE(root.array("list").size() == 2);
  REQUIRE(root.array("list")[1].child("x").integer() == 2);
  REQUIRE(root.array("missing").empty());

  // Type errors name the path.
  try {
    root.child("a").child("c").number();
    FAIL("expected a ConfigError");
  } catch (const ConfigError& e) {
    REQUIRE(std::string(e.what()).find("a.c") != std::string::npos);
    REQUIRE(std::string(e.what()).find("number") != std::string::npos);
  }
  REQUIRE_THROWS_AS(root.require("nope"), ConfigError);
  REQUIRE_THROWS_AS(root.child("a").child("b").string(), ConfigError);
  REQUIRE_THROWS_AS(root.child("a").child("b").boolean(), ConfigError);
  REQUIRE_THROWS_AS(root.child("a").child("c").vector2(), ConfigError);
  REQUIRE_THROWS_AS(root.child("a").child("b").integer(), ConfigError);
  REQUIRE_THROWS_AS(root.child("list").child("x"), ConfigError);

  // positive_number and bounded_number validate ranges.
  REQUIRE(root.child("a").positive_number("b") == Approx(3.5));
  REQUIRE_THROWS_AS(root.child("a").bounded_number("b", 0.0, 1.0), ConfigError);

  // "typo_key" was never read.
  const std::vector<std::string> unused = root.unused_keys();
  REQUIRE(std::find(unused.begin(), unused.end(), "typo_key") != unused.end());
  REQUIRE(std::find(unused.begin(), unused.end(), "a.b") == unused.end());
}

TEST_CASE("a minimal deck parses into a solvable model", "[io][config]") {
  const json::Value doc = json::parse(minimal_deck(), "minimal");
  const Configuration config = parse_configuration(doc, "minimal", /*strict=*/true);

  REQUIRE(config.name == "unit_case");
  REQUIRE(config.mesh_spec.nx == 8);
  REQUIRE(config.mesh_spec.ly == Approx(0.4));
  REQUIRE(config.material().youngs_modulus() == Approx(70.0e9));
  REQUIRE(config.material().poisson_ratio() == Approx(0.33));
  REQUIRE(config.thickness == Approx(0.01));
  REQUIRE(config.stress_state == StressState::PlaneStress);
  REQUIRE(config.constraints.size() == 1);
  REQUIRE(config.constraints.front().fix_x);
  REQUIRE(config.constraints.front().fix_y);
  REQUIRE(config.load_cases.size() == 1);
  REQUIRE(config.load_cases.front().point_loads.size() == 1);
  REQUIRE(config.load_cases.front().point_loads.front().force.y() == Approx(-1000.0));
  REQUIRE(config.load_cases.front().point_loads.front().distribute_total);
  // Defaults.
  REQUIRE(config.integration.stiffness_points == 2);
  REQUIRE(config.integration.mass_points == 3);
  // "auto": exact Cholesky below the size limits, multigrid CG above.
  REQUIRE(config.analysis.linear.type == LinearSolverType::Auto);
  REQUIRE(config.modal.options.linear.type == LinearSolverType::Auto);
  REQUIRE_FALSE(config.topology.optimizer.projection.enabled);
  REQUIRE_FALSE(config.modal.enabled);
  REQUIRE_FALSE(config.topology.enabled);

  FemModel model = build_model(config);
  REQUIRE(model.mesh().num_elements() == 32);
  REQUIRE(model.finalized());
  Assembler assembler(model);
  StaticAnalysis analysis(model, assembler, config.analysis);
  const StaticSolution sol = analysis.solve_all().front();
  REQUIRE(sol.compliance > 0.0);
  REQUIRE(sol.equilibrium.relative_force_error == Approx(0.0).margin(1.0e-10));
}

TEST_CASE("configuration errors are specific and actionable",
          "[io][config][diagnostics]") {
  const auto parse_modified = [](const std::string& body) {
    const json::Value doc = json::parse(body, "bad");
    return parse_configuration(doc, "bad");
  };

  REQUIRE_THROWS_AS(parse_modified("[1, 2]"), ConfigError);
  REQUIRE_THROWS_AS(parse_modified("{}"), ConfigError);

  // Missing mesh entries.
  REQUIRE_THROWS_AS(
      parse_modified(R"({"mesh": {"nx": 2, "ny": 2, "ly": 1.0},
                         "material": {"youngs_modulus": 1, "poisson_ratio": 0.3},
                         "boundary_conditions": [], "load_cases": []})"),
      ConfigError);

  // Unknown mesh type.
  REQUIRE_THROWS_AS(
      parse_modified(R"({"mesh": {"type": "tri3", "nx": 2, "ny": 2, "lx": 1, "ly": 1},
                         "material": {"youngs_modulus": 1, "poisson_ratio": 0.3},
                         "boundary_conditions": [], "load_cases": []})"),
      ConfigError);

  // Empty boundary conditions.
  REQUIRE_THROWS_AS(
      parse_modified(R"({"mesh": {"nx": 2, "ny": 2, "lx": 1, "ly": 1},
                         "material": {"youngs_modulus": 1e9, "poisson_ratio": 0.3},
                         "boundary_conditions": [],
                         "load_cases": [{"point_loads": [
                            {"force": [1,0], "region": {"all": true}}]}]})"),
      ConfigError);

  // A load case with no loads.
  REQUIRE_THROWS_AS(
      parse_modified(R"({"mesh": {"nx": 2, "ny": 2, "lx": 1, "ly": 1},
                         "material": {"youngs_modulus": 1e9, "poisson_ratio": 0.3},
                         "boundary_conditions": [
                            {"fix": ["x","y"], "region": {"box": {"xmax": 0}}}],
                         "load_cases": [{"name": "empty"}]})"),
      ConfigError);

  // An invalid fix component.
  REQUIRE_THROWS_AS(
      parse_modified(R"({"mesh": {"nx": 2, "ny": 2, "lx": 1, "ly": 1},
                         "material": {"youngs_modulus": 1e9, "poisson_ratio": 0.3},
                         "boundary_conditions": [
                            {"fix": ["z"], "region": {"box": {"xmax": 0}}}],
                         "load_cases": [{"point_loads": [
                            {"force": [1,0], "region": {"all": true}}]}]})"),
      ConfigError);

  // A region naming two primitives at once.
  REQUIRE_THROWS_AS(
      parse_modified(R"({"mesh": {"nx": 2, "ny": 2, "lx": 1, "ly": 1},
                         "material": {"youngs_modulus": 1e9, "poisson_ratio": 0.3},
                         "boundary_conditions": [
                            {"fix": ["x"], "region": {"box": {"xmax": 0},
                                                      "circle": {"center": [0,0],
                                                                 "radius": 1}}}],
                         "load_cases": [{"point_loads": [
                            {"force": [1,0], "region": {"all": true}}]}]})"),
      ConfigError);

  // A region naming no primitive.
  REQUIRE_THROWS_AS(
      parse_modified(R"({"mesh": {"nx": 2, "ny": 2, "lx": 1, "ly": 1},
                         "material": {"youngs_modulus": 1e9, "poisson_ratio": 0.3},
                         "boundary_conditions": [
                            {"fix": ["x"], "region": {"name": "empty"}}],
                         "load_cases": [{"point_loads": [
                            {"force": [1,0], "region": {"all": true}}]}]})"),
      ConfigError);

  // Modal analysis without a density.
  REQUIRE_THROWS_AS(
      parse_modified(R"({"mesh": {"nx": 2, "ny": 2, "lx": 1, "ly": 1},
                         "material": {"youngs_modulus": 1e9, "poisson_ratio": 0.3},
                         "boundary_conditions": [
                            {"fix": ["x","y"], "region": {"box": {"xmax": 0}}}],
                         "load_cases": [{"point_loads": [
                            {"force": [1,0], "region": {"all": true}}]}],
                         "modal": {"enabled": true}})"),
      ConfigError);

  // An empty BC region selects no nodes: caught at finalize().
  {
    const json::Value doc = json::parse(R"({
      "mesh": {"nx": 2, "ny": 2, "lx": 1, "ly": 1},
      "material": {"youngs_modulus": 1e9, "poisson_ratio": 0.3, "density": 1},
      "boundary_conditions": [
        {"name": "nowhere", "fix": ["x","y"],
         "region": {"box": {"xmin": 100.0}}}],
      "load_cases": [{"point_loads": [
        {"force": [1,0], "region": {"all": true}}]}]
    })",
                                        "empty_bc");
    const Configuration config = parse_configuration(doc, "empty_bc");
    REQUIRE_THROWS_AS(build_model(config), ConfigError);
  }
}

TEST_CASE("unknown configuration keys are reported", "[io][config][diagnostics]") {
  const std::string deck = R"({
    "name": "typo_case",
    "mesh": { "nx": 4, "ny": 2, "lx": 1.0, "ly": 0.5 },
    "material": { "youngs_modulus": 1e9, "poisson_ratio": 0.3, "density": 1000 },
    "model": { "thicknes": 0.01 },
    "boundary_conditions": [
      { "fix": ["x","y"], "region": { "box": { "xmax": 0.0 } } }
    ],
    "load_cases": [
      { "point_loads": [ { "force": [0,-1], "region": { "all": true } } ] }
    ]
  })";
  const json::Value doc = json::parse(deck, "typo");
  // In strict mode the typo is an error; the misspelled key means `thickness`
  // silently takes its default otherwise.
  REQUIRE_THROWS_AS(parse_configuration(doc, "typo", /*strict=*/true), ConfigError);
  // Non-strict mode warns and continues with the default thickness.
  const Configuration config = parse_configuration(doc, "typo", /*strict=*/false);
  REQUIRE(config.thickness == Approx(1.0));
}

TEST_CASE("a topology deck parses with filter, passive regions and overrides",
          "[io][config][topopt]") {
  const std::string deck = R"({
    "name": "topo_case",
    "mesh": { "nx": 20, "ny": 10, "lx": 1.0, "ly": 0.5 },
    "material": { "youngs_modulus": 70e9, "poisson_ratio": 0.3, "density": 2700 },
    "model": { "thickness": 0.01 },
    "boundary_conditions": [
      { "fix": ["x","y"], "region": { "box": { "xmax": 0.0 } } }
    ],
    "load_cases": [
      { "name": "down", "weight": 2.0,
        "point_loads": [ { "force": [0,-1000], "region": { "box": { "xmin": 1.0 } } } ] },
      { "name": "side", "weight": 1.0,
        "tractions": [ { "traction": [1e5, 0], "region": { "box": { "ymin": 0.5 } } } ] }
    ],
    "modal": { "enabled": true, "num_modes": 4, "mass_type": "lumped" },
    "topology": {
      "enabled": true,
      "volume_fraction": 0.35,
      "simp": { "penalty": 3.5, "emin_ratio": 1e-8,
                "mass_interpolation": "linear" },
      "filter": { "type": "density", "radius_elements": 2.0 },
      "optimizer": { "max_iterations": 50, "move_limit": 0.15,
                     "continuation_steps": 2, "penalty_start": 1.5,
                     "continuation_iterations": 10 },
      "passive_regions": [
        { "name": "hole", "type": "void",
          "region": { "circle": { "center": [0.25, 0.25], "radius": 0.06 } } },
        { "name": "collar", "type": "solid",
          "region": { "annulus": { "center": [0.25, 0.25],
                                   "inner_radius": 0.06, "radius": 0.09 } } }
      ]
    },
    "output": { "vtk": false, "csv": true }
  })";
  const json::Value doc = json::parse(deck, "topo");
  const Configuration config = parse_configuration(doc, "topo", /*strict=*/true);

  REQUIRE(config.topology.enabled);
  REQUIRE(config.topology.volume_fraction == Approx(0.35));
  REQUIRE(config.topology.optimizer.simp.penalty == Approx(3.5));
  REQUIRE(config.topology.optimizer.simp.mass_law == MassInterpolation::Linear);
  REQUIRE(config.topology.filter_type == FilterType::Density);
  REQUIRE(config.topology.filter_radius_elements == Approx(2.0));
  REQUIRE(config.topology.optimizer.oc.move_limit == Approx(0.15));
  REQUIRE(config.topology.optimizer.continuation_steps == 2);
  REQUIRE(config.topology.passive_regions.size() == 2);
  REQUIRE_FALSE(config.topology.passive_regions[0].solid);
  REQUIRE(config.topology.passive_regions[1].solid);
  REQUIRE(config.modal.enabled);
  REQUIRE(config.modal.options.mass_type == MassType::Lumped);
  REQUIRE_FALSE(config.output.write_vtk);
  REQUIRE(config.load_cases.size() == 2);
  REQUIRE(config.load_cases[1].tractions.size() == 1);

  FemModel model = build_model(config);
  REQUIRE(model.load_vectors().size() == 2);
  const std::vector<Scalar> weights = model.normalised_weights();
  REQUIRE(weights[0] == Approx(2.0 / 3.0));

  // radius_elements resolves against the mesh: cell size is 0.05 m.
  REQUIRE(config.resolved_filter_radius(model.mesh()) == Approx(0.1));

  const DesignDomain domain = build_design_domain(config, model);
  REQUIRE(domain.num_passive_void() > 0);
  REQUIRE(domain.num_passive_solid() > 0);
}

TEST_CASE("the topology summary records the load cases behind its objective",
          "[io][writers][topopt]") {
  // The design study reads per-load-case compliance out of summary.json and
  // recombines it, so the summary has to name the cases and record the
  // normalised weights actually used. Without the normalised weights the
  // reported objective cannot be reconstructed, because the objective is a
  // weighted *mean*, not a weighted sum.
  const std::string deck = R"({
    "name": "summary_case",
    "mesh": { "nx": 16, "ny": 8, "lx": 0.8, "ly": 0.4 },
    "material": { "youngs_modulus": 70e9, "poisson_ratio": 0.3, "density": 2700 },
    "model": { "thickness": 0.01 },
    "boundary_conditions": [
      { "fix": ["x","y"], "region": { "box": { "xmax": 0.0 } } }
    ],
    "load_cases": [
      { "name": "down", "weight": 3.0,
        "point_loads": [ { "force": [0,-1000], "region": { "box": { "xmin": 0.8 } } } ] },
      { "name": "side", "weight": 1.0,
        "point_loads": [ { "force": [500,0], "region": { "box": { "xmin": 0.8 } } } ] }
    ],
    "topology": {
      "enabled": true,
      "volume_fraction": 0.4,
      "filter": { "type": "density", "radius_elements": 1.5 },
      "optimizer": { "max_iterations": 12, "continuation_steps": 1 }
    }
  })";
  const Configuration config =
      parse_configuration(json::parse(deck, "summary"), "summary", /*strict=*/true);
  FemModel model = build_model(config);
  const DesignDomain domain = build_design_domain(config, model);
  const DensityFilter filter(model.mesh(), config.topology.filter_type,
                             config.resolved_filter_radius(model.mesh()));
  Assembler assembler(model);
  TopologyOptimizer optimizer(model, assembler, filter, domain,
                              config.topology.optimizer);
  const TopologyOptimizationResult result = optimizer.run();

  TimingLedger timings;
  const json::Value summary =
      make_topology_summary(config, model, domain, filter, result, nullptr, nullptr,
                            nullptr, nullptr, timings);
  const json::Value& setup = *summary.find("optimization_setup");
  const json::Value& res = *summary.find("optimization_result");

  const std::vector<json::Value>& names = setup.find("load_case_names")->array_items();
  REQUIRE(names.size() == 2);
  REQUIRE(names[0].string_value() == "down");
  REQUIRE(names[1].string_value() == "side");

  // Deck weights are recorded as written...
  const std::vector<json::Value>& raw =
      setup.find("load_case_weights")->array_items();
  REQUIRE(raw[0].number_value() == Approx(3.0));
  REQUIRE(raw[1].number_value() == Approx(1.0));
  // ...and the normalised weights as used.
  const std::vector<json::Value>& norm =
      setup.find("load_case_weights_normalised")->array_items();
  REQUIRE(norm[0].number_value() == Approx(0.75));
  REQUIRE(norm[1].number_value() == Approx(0.25));

  // The headline objective reconstructs exactly from the recorded parts.
  const std::vector<json::Value>& per_case =
      res.find("load_case_compliance_J")->array_items();
  REQUIRE(per_case.size() == 2);
  Scalar recombined = 0.0;
  for (std::size_t l = 0; l < per_case.size(); ++l) {
    recombined += norm[l].number_value() * per_case[l].number_value();
  }
  REQUIRE(recombined ==
          Approx(res.find("compliance_J")->number_value()).epsilon(1.0e-12));
}

TEST_CASE("CSV and VTK writers produce well-formed files", "[io][writers]") {
  ensure_directory("results/_test_tmp");

  SECTION("CSV enforces the declared column count") {
    const std::string path = temp_path("rows.csv");
    {
      CsvWriter csv(path, {"i", "a[m]", "b[N]"}, 9);
      csv.row(0, {1.5, -2.5});
      csv.row({3.0, 4.0, 5.0});
      csv.raw_row({"x", "y", "z"});
      REQUIRE_THROWS_AS(csv.row({1.0}), IoError);
      REQUIRE_THROWS_AS(csv.row(1, {1.0}), IoError);
      REQUIRE_THROWS_AS(csv.raw_row({"only_one"}), IoError);
      csv.close();
    }
    std::ifstream in(path);
    REQUIRE(in.good());
    std::string line;
    std::getline(in, line);
    REQUIRE(line == "i,a[m],b[N]");
    std::getline(in, line);
    REQUIRE(line == "0,1.5,-2.5");
    std::getline(in, line);
    REQUIRE(line == "3,4,5");
    std::getline(in, line);
    REQUIRE(line == "x,y,z");
    std::remove(path.c_str());
  }

  SECTION("CSV reports an unwritable path") {
    REQUIRE_THROWS_AS(CsvWriter("/definitely/not/a/dir/x.csv", {"a"}), IoError);
  }

  SECTION("VTK writes the expected header and data blocks") {
    FemModel model = make_small_plate(3, 2);
    Assembler assembler(model);
    StaticAnalysis analysis(model, assembler, StaticAnalysisOptions());
    const Vector u = analysis.solve_all().front().displacement;

    const std::string path = temp_path("fields.vtk");
    VtkWriter writer(model.mesh(), "unit test");
    writer.add_point_vectors("displacement", u);
    writer.add_point_scalars("dummy_point", Vector::Ones(model.mesh().num_nodes()));
    writer.add_cell_scalars("density", Vector::Ones(model.mesh().num_elements()));
    writer.write(path);

    std::ifstream in(path);
    std::stringstream buffer;
    buffer << in.rdbuf();
    const std::string content = buffer.str();
    REQUIRE(content.find("# vtk DataFile Version 3.0") == 0);
    REQUIRE(content.find("DATASET UNSTRUCTURED_GRID") != std::string::npos);
    REQUIRE(content.find("POINTS 12 double") != std::string::npos);
    REQUIRE(content.find("CELLS 6 30") != std::string::npos);
    REQUIRE(content.find("CELL_TYPES 6") != std::string::npos);
    REQUIRE(content.find("VECTORS displacement double") != std::string::npos);
    REQUIRE(content.find("SCALARS density double 1") != std::string::npos);
    std::remove(path.c_str());

    // Length mismatches are rejected.
    VtkWriter bad(model.mesh());
    REQUIRE_THROWS_AS(bad.add_point_scalars("x", Vector::Ones(3)), IoError);
    REQUIRE_THROWS_AS(bad.add_cell_scalars("y", Vector::Ones(3)), IoError);
    REQUIRE_THROWS_AS(bad.add_point_vectors("z", Vector::Ones(3)), IoError);
    REQUIRE_THROWS_AS(bad.write("/definitely/not/a/dir/x.vtk"), IoError);
  }

  SECTION("directory creation is recursive and idempotent") {
    const std::string nested = "results/_test_tmp/a/b/c";
    REQUIRE_NOTHROW(ensure_directory(nested));
    REQUIRE_NOTHROW(ensure_directory(nested));
    // A path that exists as a file cannot become a directory.
    const std::string file = temp_path("not_a_dir");
    std::ofstream(file) << "x";
    REQUIRE_THROWS_AS(ensure_directory(file), IoError);
    std::remove(file.c_str());
  }
}

TEST_CASE("CalculiX decks are written per load case with the matching element",
          "[io][writers][cross-validation]") {
  ensure_directory("results/_test_tmp");
  FemModel model = make_small_plate(3, 2);
  const std::vector<std::string> decks =
      write_calculix_decks(model, "results/_test_tmp/ccx", "unit");
  REQUIRE(decks.size() == 1);
  REQUIRE(calculix_element_type(model) == "CPS4");
  std::ifstream in(decks.front());
  std::stringstream buffer;
  buffer << in.rdbuf();
  const std::string text = buffer.str();
  REQUIRE(text.find("*ELEMENT, TYPE=CPS4, ELSET=EALL") != std::string::npos);
  REQUIRE(text.find("*NODE, NSET=NALL\n1, 0, 0, 0\n") != std::string::npos);
  // The thickness follows the section header; it is written with full
  // precision, so parse it back rather than matching its text.
  const std::string section_header = "*SOLID SECTION, ELSET=EALL, MATERIAL=MAT1\n";
  const std::size_t section_at = text.find(section_header);
  REQUIRE(section_at != std::string::npos);
  const std::size_t thickness_at = section_at + section_header.size();
  const std::size_t thickness_end = text.find('\n', thickness_at);
  REQUIRE(thickness_end != std::string::npos);
  const double thickness = std::stod(text.substr(thickness_at, thickness_end - thickness_at));
  REQUIRE(thickness == Catch::Approx(0.005).epsilon(1e-15));
  REQUIRE(text.find("*BOUNDARY") != std::string::npos);
  REQUIRE(text.find("*CLOAD") != std::string::npos);
  REQUIRE(text.find("*END STEP") != std::string::npos);
  // One *BOUNDARY line per prescribed DOF, one *CLOAD line per non-zero force.
  std::size_t boundary_lines = 0;
  std::size_t cload_lines = 0;
  std::istringstream lines(text);
  std::string line;
  int section = 0;
  while (std::getline(lines, line)) {
    if (line.rfind("*BOUNDARY", 0) == 0) { section = 1; continue; }
    if (line.rfind("*CLOAD", 0) == 0) { section = 2; continue; }
    if (line.rfind("*", 0) == 0) { section = 0; continue; }
    if (section == 1) ++boundary_lines;
    if (section == 2) ++cload_lines;
  }
  REQUIRE(boundary_lines == static_cast<std::size_t>(model.dofs().num_constrained()));
  REQUIRE(cload_lines == 2);  // fx and fy at the loaded corner node
  std::remove(decks.front().c_str());

  StructuredMeshSpec spec;
  spec.nx = spec.ny = spec.nz = 2;
  FemModel solid(make_structured_hex_mesh(spec), default_material(), 1.0,
                 StressState::ThreeDimensional, IntegrationOptions());
  REQUIRE(calculix_element_type(solid) == "C3D8");
  REQUIRE_THROWS_AS(write_calculix_decks(solid, "results/_test_tmp/ccx3", "x"), IoError);
}

namespace {

std::string read_text(const std::string& path) {
  std::ifstream in(path);
  std::stringstream buffer;
  buffer << in.rdbuf();
  return buffer.str();
}

/// Lines of `text` after the card `card` up to the next card.
std::vector<std::string> card_lines(const std::string& text, const std::string& card) {
  std::vector<std::string> out;
  std::istringstream lines(text);
  std::string line;
  bool inside = false;
  while (std::getline(lines, line)) {
    if (line.rfind("*", 0) == 0) {
      inside = line == card;
      continue;
    }
    if (inside) out.push_back(line);
  }
  return out;
}

}  // namespace

TEST_CASE("CalculiX decks carry pressure, body, thermal and conduction loads natively",
          "[io][writers][cross-validation]") {
  ensure_directory("results/_test_tmp");
  StructuredMeshSpec spec;
  spec.nx = 2;
  spec.ny = 2;
  spec.nz = 2;
  IsotropicMaterial steel(200.0e9, 0.3, 7850.0, "steel");
  steel.set_thermal(1.2e-5, 293.15, 45.0);
  FemModel model(make_structured_hex_mesh(spec), steel, 1.0, StressState::ThreeDimensional,
                 IntegrationOptions());
  IsotropicMaterial aluminium(70.0e9, 0.33, 2700.0, "aluminium");
  aluminium.set_thermal(2.3e-5, 293.15, 167.0);
  model.assign_material(aluminium, {0, 1});
  DisplacementConstraint root;
  Selector x0;
  x0.kind = SelectorKind::Box;
  x0.xmax = 0.0;
  root.region.members.push_back(x0);
  root.fix_x = root.fix_y = root.fix_z = true;
  model.constraints().push_back(root);

  LoadCaseSpec loads;
  loads.name = "mixed";
  PressureLoadSpec top;
  Selector z1;
  z1.kind = SelectorKind::Box;
  z1.zmin = 1.0;
  top.region.members.push_back(z1);
  top.pressure = 1.0e5;
  loads.pressures.push_back(top);
  PressureLoadSpec bottom = top;
  Selector z0;
  z0.kind = SelectorKind::Box;
  z0.zmax = 0.0;
  bottom.region.members = {z0};
  loads.pressures.push_back(bottom);
  loads.gravity = Vector3(0.0, 0.0, -9.81);
  loads.centrifugal.enabled = true;
  loads.centrifugal.angular_velocity = 10.0;
  loads.centrifugal.axis = Vector3(0.0, 0.0, 2.0);
  loads.temperature.source = TemperatureSpec::Source::Conduction;
  RegionValue hot;
  hot.region.members.push_back(x0);
  hot.value = 373.15;
  loads.temperature.conduction.prescribed.push_back(hot);
  ConvectionSpec film;
  film.region.members.push_back(z1);
  film.film_coefficient = 25.0;
  film.ambient = 293.15;
  loads.temperature.conduction.convection.push_back(film);
  model.load_case_specs().push_back(loads);
  model.finalize();

  const std::vector<std::string> decks =
      write_calculix_decks(model, "results/_test_tmp/native", "unit");
  REQUIRE(decks.size() == 2);
  REQUIRE(decks[1].find("_conduction.inp") != std::string::npos);

  const std::string text = read_text(decks[0]);
  // One material and section per material, on element sets.
  REQUIRE(text.find("*ELSET, ELSET=M1") != std::string::npos);
  REQUIRE(text.find("*ELSET, ELSET=M2\n1, 2\n") != std::string::npos);
  REQUIRE(text.find("*SOLID SECTION, ELSET=M2, MATERIAL=MAT2") != std::string::npos);
  REQUIRE(text.find("*DENSITY\n7850\n") != std::string::npos);
  REQUIRE(text.find("*EXPANSION, ZERO=293.14999999999998\n2.3e-05\n") != std::string::npos);
  REQUIRE(text.find("*INITIAL CONDITIONS, TYPE=TEMPERATURE\nNALL, 293.14999999999998") !=
          std::string::npos);
  // The pressure goes out as face loads: the top faces (z = 1) are CalculiX's
  // face 2 of a C3D8, the bottom faces its face 1, four of each; no nodal
  // force remains.
  const std::vector<std::string> dload = card_lines(text, "*DLOAD");
  int p1 = 0;
  int p2 = 0;
  for (const std::string& line : dload) {
    if (line.find(", P1, ") != std::string::npos) ++p1;
    if (line.find(", P2, ") != std::string::npos) ++p2;
  }
  REQUIRE(p1 == 4);
  REQUIRE(p2 == 4);
  REQUIRE(text.find("*CLOAD") == std::string::npos);
  REQUIRE(text.find("EALL, GRAV, 9.8100000000000005, 0, 0, -1\n") != std::string::npos);
  REQUIRE(text.find("EALL, CENTRIF, 100, 0, 0, 0, 0, 0, 1\n") != std::string::npos);
  REQUIRE(card_lines(text, "*TEMPERATURE").size() ==
          static_cast<std::size_t>(model.mesh().num_nodes()));

  const std::string heat = read_text(decks[1]);
  REQUIRE(heat.find("*ELEMENT, TYPE=DC3D8, ELSET=EALL") != std::string::npos);
  REQUIRE(heat.find("*HEAT TRANSFER, STEADY STATE") != std::string::npos);
  REQUIRE(heat.find("*CONDUCTIVITY\n167\n") != std::string::npos);
  // Nine nodes on x = 0 held on DOF 11; four convecting faces, CalculiX's F2.
  REQUIRE(card_lines(heat, "*BOUNDARY").size() == 9);
  for (const std::string& line : card_lines(heat, "*BOUNDARY")) {
    REQUIRE(line.find(", 11, 11, 373.") != std::string::npos);
  }
  const std::vector<std::string> films = card_lines(heat, "*FILM");
  REQUIRE(films.size() == 4);
  for (const std::string& line : films) {
    REQUIRE(line.find(", F2, 293.14999999999998, 25") != std::string::npos);
  }
  for (const std::string& path : decks) std::remove(path.c_str());
}

TEST_CASE("path_join handles separators", "[io]") {
  REQUIRE(path_join("a", "b") == "a/b");
  REQUIRE(path_join("a/", "b") == "a/b");
  REQUIRE(path_join("", "b") == "b");
}
