/// \file verify_logarithmic.cpp
/// \brief Verification of the logarithmic-strain finite plasticity
///        (`kinematics: finite_logarithmic`) at large strain.
///
/// Studies:
///   * `logarithmic-uniaxial`  a block on symmetry planes, a distorted
///                             2 x 2 x 2 Hex8 patch, pulled along x: to a
///                             stretch of 2 (log strain 0.69) with J2 and
///                             linear or Voce hardening, through a
///                             tension-compression cycle to +-0.3 log strain
///                             with Armstrong-Frederick backstresses, and to a
///                             stretch of 1.6 along the rolling and the
///                             transverse direction of a Hill48 sheet. For a
///                             coaxial homogeneous stretch the log strains of
///                             the increments add, so the Kirchhoff stress
///                             against the axial log strain is the 1-D
///                             small-strain law exactly, and Hill's plastic
///                             width over thickness log strain is r exactly;
///   * `logarithmic-tube`      the thick tube of `plastic-cylinder` driven by a
///                             follower bore pressure along the arc-length
///                             path past its collapse, now at finite strain:
///                             Q4 and Hex8 with the mean dilatation of ln J
///                             converge to the collapse pressure of the
///                             locking-free Tet10, the fully integrated Q4
///                             locks.
#include "VerifySupport.hpp"

#include "AppSupport.hpp"

#include "sparlab/core/Timer.hpp"
#include "sparlab/fem/Assembler.hpp"
#include "sparlab/fem/NonlinearStatic.hpp"
#include "sparlab/io/CsvWriter.hpp"
#include "sparlab/mesh/StructuredMesh.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <vector>

namespace sparlab {
namespace verify {
namespace {

constexpr Scalar kInf = std::numeric_limits<Scalar>::infinity();

std::string fmt(Scalar v, int digits = 6) {
  return std::isnan(v) ? std::string() : app::format(v, digits);
}

/// The nodes on the plane x_axis = value.
SelectorGroup plane_at(int axis, Scalar value) {
  SelectorGroup g;
  g.name = "plane";
  Selector s;
  s.kind = SelectorKind::Box;
  s.xmin = axis == 0 ? value : -kInf;
  s.xmax = axis == 0 ? value : kInf;
  s.ymin = axis == 1 ? value : -kInf;
  s.ymax = axis == 1 ? value : kInf;
  s.zmin = axis == 2 ? value : -kInf;
  s.zmax = axis == 2 ? value : kInf;
  g.members.push_back(s);
  return g;
}

NonlinearMonitor monitor(const std::string& name, int axis, Scalar value, bool reaction) {
  NonlinearMonitor m;
  m.name = name;
  m.region = plane_at(axis, value);
  m.component = axis;
  m.quantity =
      reaction ? NonlinearMonitor::Quantity::Reaction : NonlinearMonitor::Quantity::Displacement;
  return m;
}

/// The 1-D law of the return in uniaxial stress, exact along a path of
/// monotonic legs: J2 with linear and Voce isotropic hardening and
/// Armstrong-Frederick backstresses, each leg's plastic flow found from the
/// branch solution X_i = nu C_i/gamma_i + (X_i0 - nu C_i/gamma_i)
/// exp(-gamma_i D) with the consistency condition (bisection on D, along
/// which the yield function falls monotonically).
struct UniaxialLaw {
  PlasticityParameters p;
  Scalar e = 0.0;
  Scalar plastic_strain = 0.0;
  Scalar accumulated = 0.0;
  std::vector<Scalar> x;

  Scalar back(int i, Scalar nu, Scalar d) const {
    const Backstress t = p.kinematic_term(i);
    const Scalar x0 = x[static_cast<std::size_t>(i)];
    if (t.recovery == 0.0) return x0 + nu * t.modulus * d;
    const Scalar limit = nu * t.modulus / t.recovery;
    return limit + (x0 - limit) * std::exp(-t.recovery * d);
  }
  Scalar back_sum(Scalar nu, Scalar d) const {
    Scalar s = 0.0;
    for (int i = 0; i < p.kinematic_terms(); ++i) s += back(i, nu, d);
    return s;
  }
  Scalar at(Scalar eps) {
    if (x.empty()) x.assign(static_cast<std::size_t>(p.kinematic_terms()), 0.0);
    const Scalar trial = e * (eps - plastic_strain);
    const Scalar relative = trial - back_sum(1.0, 0.0);
    if (std::abs(relative) <= p.yield(accumulated)) return trial;
    const Scalar nu = relative > 0.0 ? 1.0 : -1.0;
    const auto f = [&](Scalar d) {
      return nu * (e * (eps - plastic_strain - nu * d) - back_sum(nu, d)) -
             p.yield(accumulated + d);
    };
    Scalar lo = 0.0;
    Scalar hi = 1.0e-3;
    while (f(hi) > 0.0) hi *= 2.0;
    for (int it = 0; it < 200; ++it) {
      const Scalar mid = 0.5 * (lo + hi);
      (f(mid) > 0.0 ? lo : hi) = mid;
    }
    const Scalar d = 0.5 * (lo + hi);
    for (int i = 0; i < p.kinematic_terms(); ++i) x[static_cast<std::size_t>(i)] = back(i, nu, d);
    plastic_strain += nu * d;
    accumulated += d;
    return e * (eps - plastic_strain);
  }
};

}  // namespace

StudyOutcome study_logarithmic_uniaxial(const std::string& out_dir, json::Value& summary) {
  const Scalar youngs = 200.0e9;
  const Scalar poisson = 0.3;
  const Scalar length = 1.0;
  const Scalar width = 0.5;
  const Scalar height = 0.4;
  const Scalar area = width * height;
  struct Case {
    std::string name;
    PlasticityParameters p;
    Scalar pull = 0.0;               ///< u_x of the end face at lambda = 1 [m]
    std::vector<Scalar> path;        ///< load path (empty: 0 -> 1)
    int steps = 20;
    int hill = 0;                    ///< 0: J2 (1-D law); 1: Hill48 along RD; 2: along TD
  };
  PlasticityParameters linear;
  linear.yield_stress = 250.0e6;
  linear.hardening_modulus = 1.5e9;
  PlasticityParameters voce;
  voce.yield_stress = 250.0e6;
  voce.hardening_modulus = 0.5e9;
  voce.saturation_stress = 150.0e6;
  voce.saturation_rate = 12.0;
  PlasticityParameters chaboche;
  chaboche.yield_stress = 250.0e6;
  chaboche.hardening_modulus = 0.5e9;
  chaboche.saturation_stress = 60.0e6;
  chaboche.saturation_rate = 25.0;
  chaboche.num_backstresses = 2;
  chaboche.backstresses[0] = {40.0e9, 400.0};
  chaboche.backstresses[1] = {5.0e9, 40.0};
  PlasticityParameters hill = linear;
  hill.hardening_modulus = 1.0e9;
  hill.criterion = YieldCriterion::Hill48;
  hill.hill.r0 = 1.9;
  hill.hill.r45 = 1.5;
  hill.hill.r90 = 2.3;
  PlasticityParameters hill_td = hill;
  hill_td.hill.rolling_direction = Vector3::UnitY();
  const Scalar back = std::expm1(-0.3) / std::expm1(0.3);
  const std::vector<Case> cases = {
      {"J2, linear hardening, to a stretch of 2", linear, length, {}, 20, 0},
      {"J2, Voce hardening, to a stretch of 2", voce, length, {}, 20, 0},
      {"J2, Voce and two Armstrong-Frederick backstresses, +-0.3 log strain", chaboche,
       length * std::expm1(0.3), {1.0, back, 1.0}, 12, 0},
      {"Hill48 along RD, to a stretch of 1.6", hill, 0.6 * length, {}, 12, 1},
      {"Hill48 along TD, to a stretch of 1.6", hill_td, 0.6 * length, {}, 12, 2}};

  CsvWriter csv(path_join(out_dir, "logarithmic_uniaxial.csv"),
                {"case", "step", "load_factor[-]", "log_strain[-]", "kirchhoff[Pa]",
                 "kirchhoff_exact[Pa]", "stress_error[-]", "r[-]", "r_exact[-]", "r_error[-]",
                 "lateral_log_strain_y[-]", "lateral_log_strain_z[-]"});
  Scalar worst_stress = 0.0;
  Scalar worst_r = 0.0;
  Scalar largest_plastic = 0.0;
  bool completed = true;
  int most_iterations = 0;
  std::size_t steps = 0;
  json::Value records = json::Value::make_array();
  for (const Case& c : cases) {
    IsotropicMaterial material(youngs, poisson, 7800.0, "steel");
    material.set_plasticity(c.p);
    StructuredMeshSpec spec;
    spec.nx = spec.ny = spec.nz = 2;
    spec.lx = length;
    spec.ly = width;
    spec.lz = height;
    FemModel model(make_perturbed_hex_mesh(spec, 0.3), material, 1.0,
                   StressState::ThreeDimensional, IntegrationOptions());
    for (int axis = 0; axis < 3; ++axis) {
      DisplacementConstraint sym;
      sym.region = plane_at(axis, 0.0);
      sym.set(axis, true, 0.0);
      model.constraints().push_back(sym);
    }
    DisplacementConstraint end;
    end.region = plane_at(0, length);
    end.set(0, true, c.pull);
    model.constraints().push_back(end);
    LoadCaseSpec lc;
    lc.name = "pull";
    lc.prescribed_displacement_only = true;
    model.load_case_specs().push_back(lc);
    model.finalize();
    Assembler assembler(model);
    NonlinearOptions options;
    options.kinematics = Kinematics::FiniteLogarithmic;
    options.steps = c.steps;
    options.load_path = c.path;
    options.residual_tolerance = 1.0e-12;
    options.displacement_tolerance = 1.0e-12;
    options.monitors = {monitor("force", 0, length, true), monitor("uy", 1, width, false),
                        monitor("uz", 2, height, false)};
    const NonlinearResult r = NonlinearStaticAnalysis(model, assembler, options).solve(0);
    completed = completed && r.completed;
    steps += r.steps.size();
    largest_plastic = std::max(largest_plastic, r.max_plastic_strain);
    UniaxialLaw law;
    law.p = material.plasticity();
    law.e = youngs;
    const Hill48Parameters& h = material.plasticity().hill;
    // Hill48 with linear hardening along an axis: E (k sigma_y + k^2 H eps)
    // / (E + k^2 H) with k = sigma_theta / sigma_0 = 1 (RD), (F + H)^-1/2
    // (TD); the plastic width over thickness log strain r0 or r90.
    const Scalar k = c.hill == 2 ? 1.0 / std::sqrt(h.F + h.H) : 1.0;
    const Scalar r_exact = c.hill == 1 ? h.H / h.G : h.H / h.F;
    Scalar case_stress = 0.0;
    Scalar case_r = 0.0;
    for (const NonlinearStep& s : r.steps) {
      const Scalar l1 = 1.0 + s.load_factor * c.pull / length;
      const Scalar eps = std::log(l1);
      const Scalar tau = s.monitors[0] * l1 / area;  // the force times l1 over A0
      const Scalar ey = std::log1p(s.monitors[1] / width);
      const Scalar ez = std::log1p(s.monitors[2] / height);
      Scalar exact = 0.0;
      Scalar error = 0.0;
      Scalar ratio = std::nan("");
      if (c.hill == 0) {
        exact = law.at(eps);
        error = std::abs(tau - exact) / law.p.yield_stress;
      } else if (s.yielding_points > 0) {
        const Scalar hk = law.p.hardening_modulus;
        exact = youngs * (k * law.p.yield_stress + k * k * hk * eps) / (youngs + k * k * hk);
        error = std::abs(tau - exact) / law.p.yield_stress;
        // The plastic lateral log strains: the total less -nu tau / E.
        ratio = (ey + poisson * tau / youngs) / (ez + poisson * tau / youngs);
        case_r = std::max(case_r, std::abs(ratio - r_exact) / r_exact);
      } else {
        exact = youngs * eps;
        error = std::abs(tau - exact) / law.p.yield_stress;
      }
      case_stress = std::max(case_stress, error);
      most_iterations = std::max(most_iterations, s.iterations);
      csv.raw_row({c.name, fmt(static_cast<Scalar>(s.index)), fmt(s.load_factor, 10),
                   fmt(eps, 10), fmt(tau, 12), fmt(exact, 12), fmt(error, 4), fmt(ratio, 12),
                   c.hill > 0 ? fmt(r_exact, 12) : std::string(),
                   c.hill > 0 && !std::isnan(ratio) ? fmt(std::abs(ratio - r_exact) / r_exact, 4)
                                                    : std::string(),
                   fmt(ey, 10), fmt(ez, 10)});
    }
    worst_stress = std::max(worst_stress, case_stress);
    worst_r = std::max(worst_r, case_r);
    json::Value rec = json::Value::make_object();
    rec.set("case", json::Value::make_string(c.name));
    rec.set("completed", json::Value::make_bool(r.completed));
    rec.set("steps", json::Value::make_number(static_cast<Scalar>(r.steps.size())));
    rec.set("max_kirchhoff_error_over_yield", json::Value::make_number(case_stress));
    if (c.hill > 0) rec.set("max_r_error", json::Value::make_number(case_r));
    rec.set("max_log_strain", json::Value::make_number(r.max_green_strain));
    rec.set("max_equivalent_plastic_strain", json::Value::make_number(r.max_plastic_strain));
    records.push_back(rec);
  }
  csv.close();
  const Scalar worst = std::max(worst_stress, worst_r);
  json::Value block = json::Value::make_object();
  block.set("kind",
            json::Value::make_string("verification (exact large-strain uniaxial response)"));
  block.set("cases", records);
  block.set("max_kirchhoff_error_over_yield", json::Value::make_number(worst_stress));
  block.set("max_r_error", json::Value::make_number(worst_r));
  block.set("note",
            json::Value::make_string(
                "A block 1 m x 0.5 m x 0.4 m, a distorted 2 x 2 x 2 Hex8 patch on the "
                "symmetry planes x = 0, y = 0, z = 0, its face x = 1 pulled along x, with "
                "logarithmic-strain kinematics (mean dilatation of ln J): E = 200 GPa, nu = "
                "0.3, sigma_y = 250 MPa. To a stretch of 2 in 20 steps with linear (1.5 GPa) "
                "or Voce (0.5 GPa + 150 MPa at rate 12) hardening; through 0 -> +0.3 -> -0.3 "
                "-> +0.3 log strain in 12 steps a leg with Voce hardening and two "
                "Armstrong-Frederick backstresses (40 GPa, 400), (5 GPa, 40); to a stretch of "
                "1.6 in 12 steps along RD and TD of a Hill48 sheet (r = 1.9, 1.5, 2.3, linear "
                "hardening 1 GPa). The Kirchhoff stress, the end force times the stretch over "
                "the reference area, against the 1-D small-strain law at the axial log "
                "strain, exact for a coaxial homogeneous stretch; for Hill48 the plastic "
                "lateral log strains (total + nu tau / E, from the mean face displacements) "
                "against r0 and r90."));
  summary.set("logarithmic_uniaxial", block);

  StudyOutcome outcome;
  outcome.name = "large-strain uniaxial tension and cycle, logarithmic strain, vs exact";
  outcome.kind = "verification";
  outcome.metric =
      "largest Kirchhoff-stress error over sigma_y and relative r-value error, to a stretch of 2";
  outcome.value = worst;
  outcome.tolerance = 1.0e-12;
  outcome.passed = completed && worst <= 1.0e-12;
  std::ostringstream note;
  note << cases.size() << " runs, " << steps << " steps: stress "
       << app::format(worst_stress, 3) << ", r-value " << app::format(worst_r, 3)
       << "; plastic strain up to " << app::format(largest_plastic, 3)
       << ", at most " << most_iterations << " Newton iterations per step";
  outcome.note = note.str();
  return outcome;
}

StudyOutcome study_logarithmic_tube(const std::string& out_dir, json::Value& summary) {
  const Scalar a = 0.1;
  const Scalar b = 0.2;
  const Scalar youngs = 200.0e9;
  const Scalar poisson = 0.3;
  const Scalar sy = 250.0e6;
  const Scalar p_limit = 2.0 / std::sqrt(3.0) * sy * std::log(b / a);
  PlasticityParameters params;
  params.yield_stress = sy;
  IsotropicMaterial material(youngs, poisson, 7800.0, "steel");
  material.set_plasticity(params);
  struct Variant {
    ElementType type;
    bool mean_dilatation;
    std::vector<Index> ladder;
  };
  // Tet10, locking-free without averaging (plastic-cylinder), gives the
  // reference; the mean dilatation of ln J is judged on Q4 and Hex8, and
  // the fully integrated Q4 shows the locking it removes.
  const std::vector<Variant> variants = {{ElementType::Tet10, false, {4, 8}},
                                         {ElementType::Quad4, true, {4, 8, 16}},
                                         {ElementType::Hex8, true, {4, 8, 16}},
                                         {ElementType::Quad4, false, {4, 8, 16}}};
  CsvWriter csv(path_join(out_dir, "logarithmic_tube.csv"),
                {"element", "mean_dilatation", "n_r", "num_dofs", "peak_pressure_ratio[-]",
                 "difference_to_tet10[-]", "order[-]", "max_log_strain[-]",
                 "max_plastic_strain[-]", "steps", "seconds"});
  Scalar reference = std::nan("");
  Scalar worst = 0.0;
  Scalar worst_order = std::numeric_limits<Scalar>::infinity();
  bool passed = true;
  std::ostringstream note;
  json::Value block = json::Value::make_object();
  for (const Variant& v : variants) {
    std::vector<Scalar> peaks;
    std::vector<Scalar> hs;
    json::Value records = json::Value::make_array();
    Scalar order = std::nan("");
    for (const Index nr : v.ladder) {
      Mesh mesh = sector_mesh(v.type, nr, 2 * nr, a, b);
      const int dim = mesh.dim();
      FemModel model(std::move(mesh), material, 1.0,
                     dim == 2 ? StressState::PlaneStrain : StressState::ThreeDimensional,
                     IntegrationOptions());
      add_quarter_supports(model);
      LoadCaseSpec load;
      load.name = "bore_pressure";
      PressureLoadSpec bore;
      bore.region = cylinder_surface("bore", a, true);
      bore.pressure = p_limit;
      load.pressures.push_back(bore);
      model.load_case_specs().push_back(load);
      model.finalize();
      Assembler assembler(model);
      NonlinearOptions options;
      options.kinematics = Kinematics::FiniteLogarithmic;
      options.mean_dilatation = v.mean_dilatation ? MeanDilatation::All : MeanDilatation::None;
      options.method = NonlinearOptions::Method::ArcLength;
      options.target_load_factor = 1.5;
      options.steps = 20;
      options.max_steps = 60;
      options.max_arc_ratio = 5.0;
      options.residual_tolerance = 1.0e-10;
      options.displacement_tolerance = 1.0e-10;
      Timer timer;
      const NonlinearResult r = NonlinearStaticAnalysis(model, assembler, options).solve(0);
      const Scalar seconds = timer.elapsed_seconds();
      Scalar peak = 0.0;
      for (const NonlinearStep& s : r.steps) peak = std::max(peak, s.load_factor);
      peaks.push_back(peak);
      hs.push_back((b - a) / static_cast<Scalar>(nr));
      if (v.type == ElementType::Tet10) reference = peak;  // the finest Tet10
      const Scalar difference = v.type == ElementType::Tet10 ? 0.0 : peak - reference;
      if (peaks.size() > 1 && v.type != ElementType::Tet10) {
        const std::size_t n = peaks.size();
        order = std::log(std::abs((peaks[n - 2] - reference) / (peaks[n - 1] - reference))) /
                std::log(hs[n - 2] / hs[n - 1]);
      }
      csv.raw_row({to_string(v.type), v.mean_dilatation ? "yes" : "no",
                   fmt(static_cast<Scalar>(nr)),
                   fmt(static_cast<Scalar>(model.dofs().num_dofs()), 9), fmt(peak, 10),
                   fmt(difference, 4), fmt(order, 4), fmt(r.max_green_strain, 4),
                   fmt(r.max_plastic_strain, 4), fmt(static_cast<Scalar>(r.steps.size())),
                   fmt(seconds, 4)});
      json::Value rec = json::Value::make_object();
      rec.set("n_r", json::Value::make_number(static_cast<Scalar>(nr)));
      rec.set("peak_pressure_ratio", json::Value::make_number(peak));
      rec.set("difference_to_tet10", json::Value::make_number(difference));
      rec.set("max_log_strain", json::Value::make_number(r.max_green_strain));
      rec.set("max_plastic_strain", json::Value::make_number(r.max_plastic_strain));
      rec.set("steps", json::Value::make_number(static_cast<Scalar>(r.steps.size())));
      records.push_back(rec);
    }
    json::Value entry = json::Value::make_object();
    entry.set("meshes", records);
    if (v.type != ElementType::Tet10) entry.set("order", json::Value::make_number(order));
    const std::string name =
        to_string(v.type) + (v.mean_dilatation ? " mean dilatation" : " standard");
    block.set(name, entry);
    if (v.type == ElementType::Tet10) {
      note << "Tet10 " << app::format(peaks.front(), 6) << ", " << app::format(peaks.back(), 6)
           << " p_L; ";
      continue;
    }
    note << name << " " << app::format(peaks.back() - reference, 3);
    if (!std::isnan(order)) note << " (order " << app::format(order, 3) << ")";
    note << "; ";
    if (v.mean_dilatation) {
      worst = std::max(worst, std::abs(peaks.back() - reference));
      if (!std::isnan(order)) worst_order = std::min(worst_order, order);
      passed = passed && std::abs(peaks.back() - reference) <= 1.0e-3 &&
               (std::isnan(order) || order >= 1.8);
    }
  }
  csv.close();
  block.set("kind", json::Value::make_string(
                        "verification (locking-free convergence to the Tet10 collapse pressure)"));
  block.set("small_strain_collapse_pressure_Pa", json::Value::make_number(p_limit));
  block.set("note",
            json::Value::make_string(
                "The tube of plastic-cylinder (a = 0.1 m, b = 0.2 m, plane strain, E = 200 "
                "GPa, nu = 0.3, sigma_y = 250 MPa, no hardening) with logarithmic-strain "
                "kinematics and the bore pressure following the expanding bore, along the "
                "arc-length path past the peak (60 steps, target 1.5 p_L). At finite strain "
                "the peak falls below the small-strain limit p_L = (2/sqrt 3) sigma_y ln(b/a) "
                "as the wall thins; the finest Tet10 (locking-free unaveraged, section 24) is "
                "the reference for Q4 and Hex8 with the mean dilatation of ln J (judged: "
                "within 1e-3 p_L on the finest mesh, order >= 1.8) and for the fully "
                "integrated Q4, which locks."));
  summary.set("logarithmic_tube", block);

  StudyOutcome outcome;
  outcome.name = "thick tube to collapse at finite strain, logarithmic strain (mean dilatation)";
  outcome.kind = "verification";
  outcome.metric = "peak-pressure difference of Q4 and Hex8 (mean dilatation of ln J), finest "
                   "mesh, to the finest Tet10";
  outcome.value = worst;
  outcome.tolerance = 1.0e-3;
  outcome.passed = passed;
  note << "smallest order of the averaged elements " << app::format(worst_order, 3);
  outcome.note = note.str();
  return outcome;
}

}  // namespace verify
}  // namespace sparlab
