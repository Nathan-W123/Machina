/// \file verify_plasticity.cpp
/// \brief Verification of the elastoplastic analysis against exact solutions
///        of small-strain J2 plasticity.
///
/// Studies:
///   * `plastic-cylinder`  a thick tube in plane strain, perfectly plastic,
///                         driven by its bore pressure to plastic collapse
///                         along the arc-length path: the collapse pressure
///                         against the exact limit load
///                         p_L = (2/sqrt 3) sigma_y ln(b/a), the flatness of the
///                         collapse plateau, and the stress field on it against
///                         the exact fully plastic field - for Q4, Hex8 and
///                         Tet10 with mean dilatation, and without it (with
///                         Tri3) to show volumetric locking;
///   * `plastic-bending`   pure bending of a plane-stress beam by prescribed end
///                         rotations: the moment-curvature relation of an
///                         elastic-perfectly plastic rectangle, then unloading
///                         to the curvature at which the exact moment
///                         vanishes, with the residual stress it leaves;
///   * `plastic-cycle`     a bar on a distorted Hex8 mesh strained through a
///                         full cycle with linear and Voce isotropic plus
///                         kinematic hardening, against the exact uniaxial
///                         response at every step (the Bauschinger effect).
///
/// The tube and the beam are exact for the continuum the model discretises,
/// so their errors must vanish as the mesh is refined, at the rate the
/// element allows; the uniaxial state is exact on any mesh.

#include "VerifySupport.hpp"

#include "AppSupport.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Timer.hpp"
#include "sparlab/fem/Assembler.hpp"
#include "sparlab/fem/NonlinearStatic.hpp"
#include "sparlab/io/CsvWriter.hpp"
#include "sparlab/mesh/StructuredMesh.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <vector>

namespace sparlab {
namespace verify {
namespace {

std::string fmt(Scalar v, int digits = 6) {
  return std::isnan(v) ? std::string() : app::format(v, digits);
}

IsotropicMaterial plastic_material(Scalar youngs, Scalar poisson, const PlasticityParameters& p,
                                   const std::string& name) {
  IsotropicMaterial m(youngs, poisson, 7800.0, name);
  m.set_plasticity(p);
  return m;
}

SelectorGroup nodes_group(const std::string& name, const std::vector<Index>& nodes) {
  SelectorGroup g;
  g.name = name;
  Selector s;
  s.kind = SelectorKind::NodeIds;
  s.ids = nodes;
  g.members.push_back(s);
  return g;
}

NonlinearMonitor reaction_monitor(const std::string& name, const SelectorGroup& region,
                                  int component) {
  NonlinearMonitor m;
  m.name = name;
  m.region = region;
  m.component = component;
  m.quantity = NonlinearMonitor::Quantity::Reaction;
  return m;
}

// ---------------------------------------------------------------------------
// Thick cylinder to plastic collapse
// ---------------------------------------------------------------------------

struct CylinderVariant {
  ElementType type;
  bool mean_dilatation;
  std::vector<Index> ladder;  ///< cells through the wall
  bool recommended;           ///< the default choice for the element: judged
};

std::string variant_name(const CylinderVariant& v) {
  return to_string(v.type) + (v.mean_dilatation ? " mean dilatation" : " standard");
}

}  // namespace

StudyOutcome study_plastic_cylinder(const std::string& out_dir, json::Value& summary) {
  const Scalar a = 0.1;
  const Scalar b = 0.2;
  const Scalar youngs = 200.0e9;
  const Scalar poisson = 0.3;
  const Scalar sy = 250.0e6;
  const Scalar root3 = std::sqrt(3.0);
  // The exact collapse pressure of a perfectly plastic von Mises tube in
  // plane strain (any Poisson ratio: at collapse the elastic strain rates
  // vanish and the flow is isochoric), and the bore pressure at first yield
  // of Lame's elastic solution, where sigma_z = 2 nu A.
  const Scalar p_limit = 2.0 / root3 * sy * std::log(b / a);
  const Scalar p_yield =
      sy * (b * b - a * a) /
      std::sqrt(3.0 * std::pow(b, 4) + std::pow(1.0 - 2.0 * poisson, 2) * std::pow(a, 4));
  PlasticityParameters params;
  params.yield_stress = sy;
  const IsotropicMaterial material = plastic_material(youngs, poisson, params, "steel");

  // The defaults (MeanDilatation::Auto) are judged: Q4 and Hex8 with mean
  // dilatation, Tet10 without; the other variants show what the choice
  // does - the fully integrated Q4 locks, the mean-dilatation Tet10 converges
  // from below with an oscillating pressure - and constant-strain Tri3.
  const std::vector<CylinderVariant> variants = {
      {ElementType::Quad4, true, {4, 8, 16, 32}, true},
      {ElementType::Hex8, true, {4, 8, 16}, true},
      {ElementType::Tet10, false, {2, 4, 8}, true},
      {ElementType::Quad4, false, {4, 8, 16}, false},
      {ElementType::Tet10, true, {2, 4, 8}, false},
      {ElementType::Tri3, false, {4, 8, 16}, false}};
  // The arc-length path runs past the collapse (target 1.5 p_L is never
  // reached by a locking-free element): 60 steps whose arc length may grow to
  // five times the first reach the plateau and run along it for a stretch.
  const Scalar target = 1.5;

  CsvWriter csv(path_join(out_dir, "plastic_cylinder.csv"),
                {"element", "mean_dilatation", "n_r", "n_theta", "h[m]", "num_dofs",
                 "collapse_pressure_ratio[-]", "collapse_error[-]", "collapse_order[-]",
                 "plateau_rise[-]", "bore_displacement[m]", "stress_rms_error[-]",
                 "stress_order[-]", "yielded_points", "points", "steps", "termination",
                 "seconds"});
  // The load paths of every run, and the stress through the wall on the
  // plateau of each variant's finest mesh.
  CsvWriter paths(path_join(out_dir, "plastic_cylinder_paths.csv"),
                  {"element", "mean_dilatation", "n_r", "step", "pressure_ratio[-]",
                   "bore_displacement[m]", "yielding_points"});
  CsvWriter profile(path_join(out_dir, "plastic_cylinder_stress.csv"),
                    {"element", "mean_dilatation", "n_r", "r[m]", "sigma_r[Pa]",
                     "sigma_theta[Pa]", "sigma_z[Pa]", "sigma_r_exact[Pa]",
                     "sigma_theta_exact[Pa]", "sigma_z_exact[Pa]"});
  json::Value block = json::Value::make_object();
  bool passed = true;
  Scalar worst_error = 0.0;
  Scalar worst_order_shortfall = -1.0e300;
  std::ostringstream note;
  for (const CylinderVariant& v : variants) {
    std::vector<Scalar> hs;
    std::vector<Scalar> errors;
    std::vector<Scalar> stress_errors;
    json::Value records = json::Value::make_array();
    for (const Index nr : v.ladder) {
      const Index nt = 2 * nr;
      Mesh mesh = sector_mesh(v.type, nr, nt, a, b);
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
      options.kinematics = Kinematics::SmallStrain;
      options.mean_dilatation = v.mean_dilatation ? MeanDilatation::All : MeanDilatation::None;
      options.method = NonlinearOptions::Method::ArcLength;
      options.target_load_factor = target;
      options.steps = 20;
      options.max_steps = 60;
      options.max_arc_ratio = 5.0;
      options.residual_tolerance = 1.0e-10;
      options.displacement_tolerance = 1.0e-10;
      NonlinearMonitor bore_u;
      bore_u.name = "bore_ux";
      bore_u.region = nodes_group("bore_x_axis", [&] {
        std::vector<Index> out;
        for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
          const Vector3 x = model.mesh().node(n);
          if (std::abs(x.y()) < 1.0e-12 && x.x() <= a * (1.0 + 1.0e-9)) out.push_back(n);
        }
        return out;
      }());
      bore_u.component = 0;
      options.monitors.push_back(bore_u);
      Timer timer;
      NonlinearStaticAnalysis analysis(model, assembler, options);
      const NonlinearResult r = analysis.solve(0);
      const Scalar seconds = timer.elapsed_seconds();

      // The collapse pressure: the largest load factor of the path (every
      // converged state is an equilibrium within yield, so a lower bound on
      // the discrete collapse load, which the plateau attains).
      Scalar top = 0.0;
      for (const NonlinearStep& s : r.steps) top = std::max(top, s.load_factor);
      const std::size_t n_steps = r.steps.size();
      const Scalar rise = n_steps > 10 ? r.steps.back().load_factor -
                                             r.steps[n_steps - 11].load_factor
                                       : std::nan("");
      const Scalar error = top - 1.0;

      for (const NonlinearStep& s : r.steps) {
        paths.raw_row({to_string(v.type), v.mean_dilatation ? "yes" : "no",
                       fmt(static_cast<Scalar>(nr)), fmt(static_cast<Scalar>(s.index)),
                       fmt(s.load_factor, 12), fmt(s.monitors[0], 8),
                       fmt(static_cast<Scalar>(s.yielding_points))});
      }
      const bool finest = nr == v.ladder.back();

      // The stress field on the plateau against the exact fully plastic one.
      Scalar err_sq = 0.0;
      const Mesh& m = model.mesh();
      for (Index e = 0; e < m.num_elements(); ++e) {
        const Vector3 c = m.element_centroid(e);
        const Scalar radius = std::hypot(c.x(), c.y());
        const Scalar cs = c.x() / radius;
        const Scalar sn = c.y() / radius;
        const Scalar sxx = r.element_cauchy(0, e);
        const Scalar syy = r.element_cauchy(1, e);
        const Scalar sxy = r.element_cauchy(dim == 2 ? 2 : 3, e);
        const Scalar szz = dim == 2 ? r.element_cauchy_zz(e) : r.element_cauchy(2, e);
        const Scalar s_r = sxx * cs * cs + syy * sn * sn + 2.0 * sxy * sn * cs;
        const Scalar s_t = sxx * sn * sn + syy * cs * cs - 2.0 * sxy * sn * cs;
        const Scalar exact_r = 2.0 / root3 * sy * std::log(radius / b);
        const Scalar exact_t = exact_r + 2.0 / root3 * sy;
        const Scalar exact_z = 0.5 * (exact_r + exact_t);
        err_sq += (s_r - exact_r) * (s_r - exact_r) + (s_t - exact_t) * (s_t - exact_t) +
                  (szz - exact_z) * (szz - exact_z);
        if (finest) {
          profile.raw_row({to_string(v.type), v.mean_dilatation ? "yes" : "no",
                           fmt(static_cast<Scalar>(nr)), fmt(radius, 8), fmt(s_r, 10),
                           fmt(s_t, 10), fmt(szz, 10), fmt(exact_r, 10), fmt(exact_t, 10),
                           fmt(exact_z, 10)});
        }
      }
      const Scalar stress_error =
          std::sqrt(err_sq / (3.0 * static_cast<Scalar>(m.num_elements()))) / sy;
      hs.push_back((b - a) / static_cast<Scalar>(nr));
      errors.push_back(std::abs(error));
      stress_errors.push_back(stress_error);
      const std::size_t last = errors.size() - 1;
      const Scalar order =
          last > 0 ? observed_order(hs[last - 1], errors[last - 1], hs[last], errors[last])
                   : std::nan("");
      const Scalar stress_order =
          last > 0 ? observed_order(hs[last - 1], stress_errors[last - 1], hs[last],
                                    stress_errors[last])
                   : std::nan("");
      const Scalar bore_disp = r.steps.empty() ? 0.0 : r.steps.back().monitors[0];
      csv.raw_row({to_string(v.type), v.mean_dilatation ? "yes" : "no",
                   fmt(static_cast<Scalar>(nr)), fmt(static_cast<Scalar>(nt)), fmt(hs.back(), 8),
                   fmt(static_cast<Scalar>(model.dofs().num_dofs()), 9), fmt(top, 12),
                   fmt(error, 6), fmt(order, 4), fmt(rise, 4), fmt(bore_disp, 6),
                   fmt(stress_error, 6), fmt(stress_order, 4),
                   fmt(static_cast<Scalar>(r.plastic_points), 9),
                   fmt(static_cast<Scalar>(r.total_points), 9),
                   fmt(static_cast<Scalar>(n_steps)), r.termination, fmt(seconds, 4)});
      json::Value rec = json::Value::make_object();
      rec.set("n_r", json::Value::make_number(static_cast<Scalar>(nr)));
      rec.set("num_dofs", json::Value::make_number(static_cast<Scalar>(model.dofs().num_dofs())));
      rec.set("collapse_pressure_ratio", json::Value::make_number(top));
      rec.set("collapse_error", json::Value::make_number(error));
      rec.set("plateau_rise_last_10_steps", json::Value::make_number(rise));
      rec.set("bore_displacement_m", json::Value::make_number(bore_disp));
      rec.set("fully_plastic_stress_rms_error", json::Value::make_number(stress_error));
      rec.set("yielded_points", json::Value::make_number(r.plastic_points));
      rec.set("points", json::Value::make_number(r.total_points));
      rec.set("steps", json::Value::make_number(static_cast<Scalar>(n_steps)));
      rec.set("termination", json::Value::make_string(r.termination));
      records.push_back(rec);
    }
    const std::size_t last = errors.size() - 1;
    const Scalar order = observed_order(hs[last - 1], errors[last - 1], hs[last], errors[last]);
    json::Value entry = json::Value::make_object();
    entry.set("meshes", records);
    entry.set("collapse_error_order", json::Value::make_number(order));
    const Scalar stress_order = observed_order(hs[last - 1], stress_errors[last - 1], hs[last],
                                               stress_errors[last]);
    entry.set("stress_error_order", json::Value::make_number(stress_order));
    entry.set("default_choice", json::Value::make_bool(v.recommended));
    note << variant_name(v) << " " << app::format(errors.back(), 3) << " (order "
         << app::format(order, 3) << ", stress " << app::format(stress_errors.back(), 3)
         << "); ";
    if (v.recommended) {
      // The collapse pressure converges to p_L at second order or better,
      // and the stress field on the plateau with it.
      worst_error = std::max(worst_error, errors.back());
      worst_order_shortfall = std::max(worst_order_shortfall, 1.8 - order);
      passed = passed && errors.back() <= 1.0e-3 && order >= 1.8 &&
               stress_errors.back() <= 0.01;
    }
    block.set(variant_name(v), entry);
  }
  csv.close();
  paths.close();
  profile.close();

  block.set("kind", json::Value::make_string(
                        "verification (exact collapse load and fully plastic stress field)"));
  block.set("collapse_pressure_exact_Pa", json::Value::make_number(p_limit));
  block.set("first_yield_pressure_exact_Pa", json::Value::make_number(p_yield));
  block.set("note",
            json::Value::make_string(
                "Quarter section of a tube a = 0.1 m, b = 0.2 m in plane strain (E = 200 "
                "GPa, nu = 0.3, sigma_y = 250 MPa, no hardening), small strain, loaded by a "
                "bore pressure along the arc-length path past its collapse (60 steps, target "
                "1.5 p_L). Collapse pressure: the largest load factor of the path, a lower "
                "bound on the discrete collapse load which the plateau attains, against "
                "p_L = (2/sqrt 3) sigma_y ln(b/a). Plateau rise: the load factor gained over "
                "the last ten steps (zero on a true collapse plateau; a locking element "
                "keeps a residual stiffness). Stress: element-average sigma_r, sigma_theta, "
                "sigma_z at the final state against the fully plastic field sigma_r = "
                "(2/sqrt 3) sigma_y ln(r/b), sigma_theta = sigma_r + 2 sigma_y/sqrt 3, "
                "sigma_z = (sigma_r + sigma_theta)/2, RMS over sigma_y. Symmetry supports; "
                "the 3-D sections one cell deep with u_z = 0."));
  summary.set("plastic_cylinder", block);

  note << "p_L = " << app::format(p_limit, 6) << " Pa, first yield "
       << app::format(p_yield / p_limit, 4) << " p_L";
  StudyOutcome outcome;
  outcome.name = "thick tube to plastic collapse vs exact limit load (mean dilatation, locking)";
  outcome.kind = "verification";
  outcome.metric = "largest collapse-pressure error of the default elements, finest mesh";
  outcome.value = worst_error;
  outcome.tolerance = 1.0e-3;
  outcome.passed = passed;
  outcome.note = note.str();
  return outcome;
}

// ---------------------------------------------------------------------------
// Pure bending of an elastic-perfectly plastic beam
// ---------------------------------------------------------------------------

namespace {

/// The exact moment of a rectangle t x h in uniaxial bending at curvature k,
/// elastic-perfectly plastic.
Scalar bending_moment(Scalar k, Scalar youngs, Scalar sy, Scalar t, Scalar h) {
  const Scalar ky = 2.0 * sy / (youngs * h);
  const Scalar inertia = t * h * h * h / 12.0;
  if (std::abs(k) <= ky) return youngs * inertia * k;
  const Scalar mp = sy * t * h * h / 4.0;
  const Scalar ratio = ky / std::abs(k);
  return std::copysign(mp * (1.0 - ratio * ratio / 3.0), k);
}

/// The exact axial stress at the distance y from the neutral axis, loaded to
/// curvature k1 and unloaded elastically to k (strain -k y, compression on
/// top for k > 0).
Scalar bending_stress(Scalar y, Scalar k1, Scalar k, Scalar youngs, Scalar sy) {
  const Scalar loaded = std::clamp(-youngs * k1 * y, -sy, sy);
  return loaded - youngs * (k - k1) * y;
}

}  // namespace

StudyOutcome study_plastic_bending(const std::string& out_dir, json::Value& summary) {
  const Scalar length = 0.2;
  const Scalar h = 0.05;
  const Scalar t = 0.01;
  const Scalar youngs = 200.0e9;
  const Scalar poisson = 0.3;
  const Scalar sy = 250.0e6;
  const Scalar ky = 2.0 * sy / (youngs * h);
  const Scalar mp = sy * t * h * h / 4.0;
  const Scalar inertia = t * h * h * h / 12.0;
  PlasticityParameters params;
  params.yield_stress = sy;
  const IsotropicMaterial material = plastic_material(youngs, poisson, params, "steel");
  // Loaded to 3 curvatures at first yield, then back along the elastic line
  // to the curvature at which the moment vanishes.
  const Scalar k1 = 3.0 * ky;
  const Scalar m1 = bending_moment(k1, youngs, sy, t, h);
  const Scalar k_res = k1 - m1 / (youngs * inertia);
  // Curvatures at which the loading branch is compared.
  const std::vector<Scalar> ratios = {0.5, 1.0, 1.5, 2.0, 3.0};

  CsvWriter csv(path_join(out_dir, "plastic_bending.csv"),
                {"ny", "nx", "h[m]", "num_dofs", "curvature_ratio[-]", "moment[N m]",
                 "moment_exact[N m]", "moment_error[-]"});
  CsvWriter residual_csv(path_join(out_dir, "plastic_bending_residual.csv"),
                         {"ny", "y[m]", "residual_stress[Pa]", "residual_stress_exact[Pa]"});
  json::Value block = json::Value::make_object();
  json::Value records = json::Value::make_array();
  std::vector<Scalar> hs;
  std::vector<Scalar> loading_errors;
  std::vector<Scalar> residual_moments;
  std::vector<Scalar> residual_stress_errors;
  for (const Index ny : {Index{4}, Index{8}, Index{16}, Index{32}}) {
    const Index nx = ny * static_cast<Index>(std::lround(length / h));
    StructuredMeshSpec spec;
    spec.nx = nx;
    spec.ny = ny;
    spec.lx = length;
    spec.ly = h;
    const Mesh mesh = make_structured_quad_mesh(spec);
    // One run per curvature on the loading branch, and one loaded to k1 and
    // unloaded to k_res.
    Scalar worst_loading = 0.0;
    Scalar residual_moment = 0.0;
    Scalar residual_error = 0.0;
    Index num_dofs = 0;
    // The moment of the right end's reactions about mid-depth after the
    // prescribed curvature k (times the load factors of `path`).
    const auto solve = [&](Scalar k, const std::vector<Scalar>& path, NonlinearResult& out) {
      FemModel model(mesh, material, t, StressState::PlaneStress, IntegrationOptions());
      // u_x = -k (x - L/2)(y - h/2) on both ends; u_y = 0 at the left end's
      // mid-depth node.
      for (Index n = 0; n < mesh.num_nodes(); ++n) {
        const Vector3 x = mesh.node(n);
        const bool left = x.x() <= 1.0e-12;
        const bool right = x.x() >= length - 1.0e-12;
        if (!left && !right) continue;
        DisplacementConstraint bc;
        bc.region = nodes_group("end", {n});
        bc.set(0, true, -k * (x.x() - 0.5 * length) * (x.y() - 0.5 * h));
        if (left && std::abs(x.y() - 0.5 * h) < 1.0e-12) bc.set(1, true, 0.0);
        model.constraints().push_back(bc);
      }
      LoadCaseSpec lc;
      lc.name = "bend";
      lc.prescribed_displacement_only = true;
      model.load_case_specs().push_back(lc);
      model.finalize();
      Assembler assembler(model);
      NonlinearOptions options;
      options.kinematics = Kinematics::SmallStrain;
      options.steps = 10;
      options.load_path = path;
      options.residual_tolerance = 1.0e-11;
      options.displacement_tolerance = 1.0e-11;
      out = NonlinearStaticAnalysis(model, assembler, options).solve(0);
      num_dofs = model.dofs().num_dofs();
      Scalar moment = 0.0;
      for (Index n = 0; n < mesh.num_nodes(); ++n) {
        const Vector3 x = mesh.node(n);
        if (x.x() < length - 1.0e-12) continue;
        moment -= out.reactions(2 * n) * (x.y() - 0.5 * h);
      }
      return moment;
    };
    for (const Scalar ratio : ratios) {
      NonlinearResult r;
      const Scalar k = ratio * ky;
      const Scalar moment = solve(k, {1.0}, r);
      const Scalar exact = bending_moment(k, youngs, sy, t, h);
      const Scalar err = std::abs(moment - exact) / mp;
      worst_loading = std::max(worst_loading, r.completed ? err : 1.0);
      csv.raw_row({fmt(static_cast<Scalar>(ny)), fmt(static_cast<Scalar>(nx)),
                   fmt(h / static_cast<Scalar>(ny), 8), fmt(static_cast<Scalar>(num_dofs), 9),
                   fmt(ratio, 4), fmt(moment, 10), fmt(exact, 10), fmt(err, 6)});
    }
    {
      NonlinearResult r;
      const Scalar moment = solve(k1, {1.0, k_res / k1}, r);
      residual_moment = std::abs(moment) / mp;
      // The residual stress in the column of elements left of mid-span,
      // against the exact one at each element's centroid.
      Scalar err_sq = 0.0;
      Index count = 0;
      const Scalar dx = length / static_cast<Scalar>(nx);
      for (Index e = 0; e < mesh.num_elements(); ++e) {
        const Vector3 c = mesh.element_centroid(e);
        if (std::abs(c.x() - (0.5 * length - 0.5 * dx)) > 1.0e-9) continue;
        const Scalar y = c.y() - 0.5 * h;
        const Scalar exact = bending_stress(y, k1, k_res, youngs, sy);
        err_sq += (r.element_cauchy(0, e) - exact) * (r.element_cauchy(0, e) - exact);
        ++count;
        residual_csv.raw_row({fmt(static_cast<Scalar>(ny)), fmt(c.y(), 8),
                              fmt(r.element_cauchy(0, e), 10), fmt(exact, 10)});
      }
      residual_error = std::sqrt(err_sq / static_cast<Scalar>(std::max<Index>(count, 1))) / sy;
      if (!r.completed) residual_error = 1.0;
    }
    hs.push_back(h / static_cast<Scalar>(ny));
    loading_errors.push_back(worst_loading);
    residual_moments.push_back(residual_moment);
    residual_stress_errors.push_back(residual_error);
    json::Value rec = json::Value::make_object();
    rec.set("ny", json::Value::make_number(static_cast<Scalar>(ny)));
    rec.set("nx", json::Value::make_number(static_cast<Scalar>(nx)));
    rec.set("num_dofs", json::Value::make_number(static_cast<Scalar>(num_dofs)));
    rec.set("max_moment_error_over_mp", json::Value::make_number(worst_loading));
    rec.set("residual_moment_over_mp", json::Value::make_number(residual_moment));
    rec.set("residual_stress_rms_error_over_sy", json::Value::make_number(residual_error));
    records.push_back(rec);
  }
  csv.close();
  residual_csv.close();
  const std::size_t last = hs.size() - 1;
  const Scalar order =
      observed_order(hs[last - 1], loading_errors[last - 1], hs[last], loading_errors[last]);
  const Scalar residual_order = observed_order(hs[last - 1], residual_stress_errors[last - 1],
                                               hs[last], residual_stress_errors[last]);
  block.set("meshes", records);
  block.set("moment_error_order", json::Value::make_number(order));
  block.set("residual_stress_error_order", json::Value::make_number(residual_order));
  block.set("plastic_moment_Nm", json::Value::make_number(mp));
  block.set("unloaded_curvature_ratio", json::Value::make_number(k_res / ky));
  block.set("kind", json::Value::make_string("verification (exact elastoplastic bending)"));
  block.set("note",
            json::Value::make_string(
                "A beam 0.2 m x 0.05 m, 0.01 m thick, in plane stress (E = 200 GPa, "
                "nu = 0.3, sigma_y = 250 MPa, no hardening), small strain, bent by end "
                "displacements u_x = -k (x - L/2)(y - h/2) on square Q4 cells: the "
                "section is in uniaxial stress, compatible in plane stress with the free "
                "transverse strain, so M(k) = E I k up to k_y = 2 sigma_y/(E h) and "
                "M_p (1 - (k_y/k)^2 / 3) beyond, M_p = sigma_y t h^2/4. Moments from the "
                "end reactions at k = 0.5, 1, 1.5, 2, 3 k_y, error over M_p. Unloading: "
                "from 3 k_y back to the curvature at which the exact moment vanishes "
                "(k_res = 1.556 k_y, an elastic unloading that stays short of reverse "
                "yield), where the moment must vanish and the residual stress is the "
                "loaded profile less E (3 k_y - k_res) times the distance from the neutral "
                "axis; element stresses of the column left of mid-span against it, RMS "
                "over sigma_y."));
  summary.set("plastic_bending", block);

  std::ostringstream note;
  note << "moment error " << app::format(loading_errors.back(), 3) << " M_p (order "
       << app::format(order, 3) << "), residual moment " << app::format(residual_moments.back(), 3)
       << " M_p, residual stress " << app::format(residual_stress_errors.back(), 3)
       << " sigma_y (order " << app::format(residual_order, 3) << ")";
  StudyOutcome outcome;
  outcome.name = "pure bending: moment-curvature and residual stress vs exact";
  outcome.kind = "verification";
  outcome.metric = "largest moment error over M_p on the loading branch, finest mesh";
  outcome.value = loading_errors.back();
  outcome.tolerance = 1.0e-3;
  outcome.passed = loading_errors.back() <= 1.0e-3 && order >= 1.5 &&
                   residual_moments.back() <= 1.0e-3 && residual_stress_errors.back() <= 0.02;
  outcome.note = note.str();
  return outcome;
}

// ---------------------------------------------------------------------------
// A uniaxial cycle with combined hardening
// ---------------------------------------------------------------------------

StudyOutcome study_plastic_cycle(const std::string& out_dir, json::Value& summary) {
  const Scalar youngs = 200.0e9;
  const Scalar poisson = 0.3;
  PlasticityParameters params;
  params.yield_stress = 250.0e6;
  params.hardening_modulus = 1.0e9;
  params.saturation_stress = 100.0e6;
  params.saturation_rate = 30.0;
  params.kinematic_hardening_modulus = 4.0e9;
  const IsotropicMaterial material = plastic_material(youngs, poisson, params, "steel");
  const Scalar length = 1.0;
  const Scalar height = 0.1;
  const Scalar width = 0.1;
  const Scalar strain_amplitude = 0.01;

  StructuredMeshSpec spec;
  spec.nx = 4;
  spec.ny = 2;
  spec.nz = 2;
  spec.lx = length;
  spec.ly = height;
  spec.lz = width;
  FemModel model(make_perturbed_hex_mesh(spec, 0.25), material, 1.0,
                 StressState::ThreeDimensional, IntegrationOptions());
  // x held on the left face; rigid motion removed at points the uniform
  // lateral contraction leaves in place; the right face driven along x.
  const auto box_nodes = [&](Scalar xmin, Scalar xmax, Scalar ymax, Scalar zmax) {
    std::vector<Index> out;
    for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
      const Vector3 x = model.mesh().node(n);
      if (x.x() >= xmin - 1.0e-12 && x.x() <= xmax + 1.0e-12 && x.y() <= ymax + 1.0e-12 &&
          x.z() <= zmax + 1.0e-12) {
        out.push_back(n);
      }
    }
    return out;
  };
  DisplacementConstraint left;
  left.region = nodes_group("left", box_nodes(0.0, 0.0, height, width));
  left.fix_x = true;
  model.constraints().push_back(left);
  DisplacementConstraint origin;
  origin.region = nodes_group("origin", box_nodes(0.0, 0.0, 0.0, 0.0));
  origin.fix_y = origin.fix_z = true;
  model.constraints().push_back(origin);
  DisplacementConstraint turn;
  turn.region = nodes_group("turn", [&] {
    std::vector<Index> out;
    for (Index n : box_nodes(0.0, 0.0, height, 0.0)) {
      if (model.mesh().node(n).y() >= height - 1.0e-12) out.push_back(n);
    }
    return out;
  }());
  turn.fix_z = true;
  model.constraints().push_back(turn);
  const std::vector<Index> right_nodes = box_nodes(length, length, height, width);
  DisplacementConstraint pull;
  pull.region = nodes_group("right", right_nodes);
  pull.set(0, true, strain_amplitude * length);
  model.constraints().push_back(pull);
  LoadCaseSpec lc;
  lc.name = "cycle";
  lc.prescribed_displacement_only = true;
  model.load_case_specs().push_back(lc);
  model.finalize();
  Assembler assembler(model);
  NonlinearOptions options;
  options.kinematics = Kinematics::SmallStrain;
  options.steps = 20;
  options.load_path = {1.0, -1.0, 1.0};
  options.residual_tolerance = 1.0e-12;
  options.displacement_tolerance = 1.0e-12;
  options.monitors.push_back(reaction_monitor("force", nodes_group("right", right_nodes), 0));
  const NonlinearResult r = NonlinearStaticAnalysis(model, assembler, options).solve(0);

  // The exact uniaxial response along the same strains: at every level the
  // state satisfies sigma = E (eps - eps_p) and, while flowing,
  // sigma - H_kin eps_p = s (sigma_y(alpha)) with s the flow direction; each
  // step is monotonic in strain, over which the flow keeps its direction.
  const PlasticityParameters& p = params;
  Scalar eps_p = 0.0;
  Scalar alpha = 0.0;
  const Scalar area = height * width;
  CsvWriter csv(path_join(out_dir, "plastic_cycle.csv"),
                {"step", "load_factor[-]", "strain[-]", "stress[Pa]", "stress_exact[Pa]",
                 "error[-]", "yielding_points"});
  Scalar worst = 0.0;
  for (const NonlinearStep& s : r.steps) {
    const Scalar eps = s.load_factor * strain_amplitude;
    const Scalar trial = youngs * (eps - eps_p);
    const Scalar relative = trial - p.kinematic_hardening_modulus * eps_p;
    Scalar exact = trial;
    if (std::abs(relative) > p.yield(alpha)) {
      // Flow in the direction of the relative stress: solve
      // |relative| - (E + H_kin) d - sigma_y(alpha + d) = 0 for d >= 0.
      const Scalar sign = relative > 0.0 ? 1.0 : -1.0;
      Scalar d = 0.0;
      for (int it = 0; it < 100; ++it) {
        const Scalar g = std::abs(relative) - (youngs + p.kinematic_hardening_modulus) * d -
                         p.yield(alpha + d);
        const Scalar dg = -(youngs + p.kinematic_hardening_modulus) - p.yield_slope(alpha + d);
        d -= g / dg;
        if (std::abs(g) <= 1.0e-14 * p.yield_stress) break;
      }
      eps_p += sign * d;
      alpha += d;
      exact = youngs * (eps - eps_p);
    }
    const Scalar stress = s.monitors[0] / area;
    const Scalar err = std::abs(stress - exact) / p.yield_stress;
    worst = std::max(worst, err);
    csv.raw_row({fmt(static_cast<Scalar>(s.index)), fmt(s.load_factor, 8), fmt(eps, 8),
                 fmt(stress, 12), fmt(exact, 12), fmt(err, 4),
                 fmt(static_cast<Scalar>(s.yielding_points))});
  }
  csv.close();
  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("verification (exact uniaxial response)"));
  block.set("completed", json::Value::make_bool(r.completed));
  block.set("steps", json::Value::make_number(static_cast<Scalar>(r.steps.size())));
  block.set("max_stress_error_over_yield", json::Value::make_number(worst));
  block.set("final_equivalent_plastic_strain", json::Value::make_number(alpha));
  block.set("note",
            json::Value::make_string(
                "A bar 1 m x 0.1 m x 0.1 m on a distorted 4 x 2 x 2 Hex8 mesh (mean "
                "dilatation), strained along x through eps = 0 -> 0.01 -> -0.01 -> 0.01 in "
                "20 steps per leg: E = 200 GPa, nu = 0.3, sigma_y = 250 MPa, linear "
                "isotropic hardening 1 GPa plus Voce 100 MPa at rate 30, Prager kinematic "
                "hardening 4 GPa. The end force over the area against the exact uniaxial "
                "stress at every step; a homogeneous state is exact on any mesh."));
  summary.set("plastic_cycle", block);

  StudyOutcome outcome;
  outcome.name = "uniaxial cycle with combined hardening vs exact (distorted Hex8)";
  outcome.kind = "verification";
  outcome.metric = "largest stress error over sigma_y along the cycle";
  outcome.value = worst;
  outcome.tolerance = 1.0e-9;
  outcome.passed = r.completed && worst <= 1.0e-9;
  std::ostringstream note;
  note << r.steps.size() << " steps, final accumulated plastic strain "
       << app::format(alpha, 4);
  outcome.note = note.str();
  return outcome;
}

}  // namespace verify
}  // namespace sparlab
