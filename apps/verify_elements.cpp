/// \file verify_elements.cpp
/// \brief Verification of the incompatible-mode Hex8 (Hex8Incompatible.hpp)
///        and the through-thickness rule of a sheet in bending: the
///        locking of the standard Hex8 that they remove, measured against
///        beam and plate theory and against a converged fine mesh.
///
/// Studies:
///   * `im-cantilever`     the slender cantilever of the Hex8 studies (1 m x
///                         0.1 m x 0.05 m, E = 70 GPa) with ONE element
///                         through its depth, nx x 1 x 1 for nx = 10, 20, 40
///                         and nu = 0 and 0.3, against Timoshenko (and
///                         Euler-Bernoulli), and a sweep of L/h = 10, 100,
///                         1000 on 20 x 1 x 1: the error of the incompatible
///                         modes does not grow with the slenderness, the
///                         standard Hex8's grows as (L/h)^2;
///   * `im-macneal-harder` the straight beam of MacNeal and Harder (1985, *A
///                         proposed standard set of problems to test finite
///                         element accuracy*, Finite Elem. Anal. Des. 1): 6 x
///                         0.2 x 0.1, E = 1e7, nu = 0.3, 6 x 1 x 1 cells,
///                         rectangular, parallelogram (45 deg) and
///                         trapezoidal (+-45 deg) cells, unit tip loads in and
///                         out of plane, against the references 0.1081 and
///                         0.4321. The trapezoidal degradation is reported,
///                         not judged: it is the known limit of the modes on
///                         non-affine cells (MacNeal 1987);
///   * `im-plate`          a square plate, t/a = 1/50, clamped and simply
///                         supported (at the mid-surface), under a uniform
///                         pressure, on a quarter with TWO element layers:
///                         the centre deflection against Kirchhoff's
///                         (Timoshenko and Woinowsky-Krieger 1959);
///   * `im-springback`     a sheet strip bent past yield by prescribed end
///                         rotations and released (the forming driver):
///                         ONE element layer with the incompatible modes and
///                         a 2 x 2 x 5 rule against a converged fine
///                         standard mesh (and 1 layer with 2 x 2 x 2, and
///                         the standard element on 1 layer, both reported).

#include "VerifySupport.hpp"

#include "AppSupport.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Timer.hpp"
#include "sparlab/elements/Quadrature.hpp"
#include "sparlab/fem/Assembler.hpp"
#include "sparlab/fem/Forming.hpp"
#include "sparlab/fem/StaticAnalysis.hpp"
#include "sparlab/io/CsvWriter.hpp"
#include "sparlab/mesh/StructuredMesh.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <sstream>
#include <vector>

namespace sparlab {
namespace verify {
namespace {

std::string fmt(Scalar v, int digits = 6) { return app::format(v, digits); }

json::Value num(Scalar v) { return json::Value::make_number(v); }
json::Value str(const std::string& s) { return json::Value::make_string(s); }

IntegrationOptions options_of(ElementFormulation formulation, int thickness_points = 0,
                              int thickness_axis = 2) {
  IntegrationOptions o;
  o.formulation = formulation;
  o.thickness_points = thickness_points;
  o.thickness_axis = thickness_axis;
  return o;
}

const char* short_name(ElementFormulation f) {
  return f == ElementFormulation::IncompatibleModes ? "incompatible" : "standard";
}

SelectorGroup box_group(const std::string& name, const Vector3& lo, const Vector3& hi) {
  SelectorGroup g;
  g.name = name;
  Selector s;
  s.kind = SelectorKind::Box;
  s.xmin = lo.x();
  s.ymin = lo.y();
  s.zmin = lo.z();
  s.xmax = hi.x();
  s.ymax = hi.y();
  s.zmax = hi.z();
  g.members.push_back(s);
  return g;
}

SelectorGroup node_group(const std::string& name, const std::vector<Index>& nodes) {
  SelectorGroup g;
  g.name = name;
  Selector s;
  s.kind = SelectorKind::NodeIds;
  s.ids = nodes;
  g.members.push_back(s);
  return g;
}

DisplacementConstraint fixed(const SelectorGroup& region, bool x, bool y, bool z) {
  DisplacementConstraint c;
  c.region = region;
  c.fix_x = x;
  c.fix_y = y;
  c.fix_z = z;
  return c;
}

/// The linear static solution of a finalised model's first load case, and
/// its relative force-balance error.
Vector solve_static(const FemModel& model, Scalar* balance = nullptr,
                    Scalar equilibrium_tolerance = 1.0e-6) {
  Assembler assembler(model);
  StaticAnalysisOptions options;
  options.equilibrium_tolerance = equilibrium_tolerance;
  StaticAnalysis analysis(model, assembler, options);
  StaticSolution sol = analysis.solve_all().front();
  if (balance != nullptr) *balance = sol.equilibrium.relative_force_error;
  return std::move(sol.displacement);
}

/// Mean of displacement component `c` over the nodes inside a box.
Scalar mean_displacement(const FemModel& model, const Vector& u, const Vector3& lo,
                         const Vector3& hi, int c) {
  const Mesh& mesh = model.mesh();
  Scalar sum = 0.0;
  Index count = 0;
  for (Index n = 0; n < mesh.num_nodes(); ++n) {
    const Vector3 x = mesh.node(n);
    if ((x.array() < lo.array()).any() || (x.array() > hi.array()).any()) continue;
    sum += u(3 * n + c);
    ++count;
  }
  if (count == 0) throw ModelError("verification: no node in the measuring box");
  return sum / static_cast<Scalar>(count);
}

constexpr Scalar kEps = 1.0e-9;

// ---------------------------------------------------------------------------
// Slender cantilever, one element through the depth
// ---------------------------------------------------------------------------

struct Cantilever {
  Scalar length = 1.0;
  Scalar height = 0.1;  ///< bending direction, y
  Scalar width = 0.05;  ///< z
  Scalar youngs = 70.0e9;
  Scalar poisson = 0.0;
  Scalar load = -1000.0;  ///< tip resultant along y [N]
  Scalar inertia() const { return width * height * height * height / 12.0; }
  Scalar euler_bernoulli() const {
    return load * length * length * length / (3.0 * youngs * inertia());
  }
  Scalar timoshenko() const {
    const Scalar g = youngs / (2.0 * (1.0 + poisson));
    return euler_bernoulli() + load * length / ((5.0 / 6.0) * g * width * height);
  }
};

/// Mean tip deflection of the cantilever on an nx x ny x nz mesh.
/// The force balance of the solve is returned in `balance`: at L/h = 1000
/// the condition number of the stiffness (the axial and shear stiffness of
/// an element against the bending of the beam, (L/h)^2 apart) leaves a
/// relative imbalance of 2e-4 while the deflection is still resolved to
/// 1e-4 (it converges with nx like the others), so the solve accepts 1e-3.
Scalar cantilever_tip(const Cantilever& c, Index nx, Index ny, Index nz,
                      ElementFormulation formulation, Scalar& balance) {
  StructuredMeshSpec ms;
  ms.nx = nx;
  ms.ny = ny;
  ms.nz = nz;
  ms.lx = c.length;
  ms.ly = c.height;
  ms.lz = c.width;
  FemModel model(make_structured_hex_mesh(ms), IsotropicMaterial(c.youngs, c.poisson, 2700.0, "al"),
                 1.0, StressState::ThreeDimensional, options_of(formulation));
  const Scalar big = 1.0e9;
  model.constraints().push_back(
      fixed(box_group("root", Vector3(-big, -big, -big), Vector3(kEps, big, big)), true, true, true));
  LoadCaseSpec lc;
  lc.name = "tip";
  PointLoadSpec tip;
  tip.region = box_group("tip", Vector3(c.length - kEps * c.length, -big, -big),
                         Vector3(big, big, big));
  tip.force = Vector3(0.0, c.load, 0.0);
  tip.distribute_total = true;
  lc.point_loads.push_back(tip);
  model.load_case_specs().push_back(lc);
  model.finalize();
  const Vector u = solve_static(model, &balance, 1.0e-3);
  return mean_displacement(model, u, Vector3(c.length - kEps * c.length, -big, -big),
                           Vector3(big, big, big), 1);
}

}  // namespace

StudyOutcome study_im_cantilever(const std::string& out_dir, json::Value& summary) {
  CsvWriter csv(path_join(out_dir, "im_cantilever.csv"),
                {"formulation", "nu[-]", "L/h[-]", "nx", "tip[m]", "timoshenko[m]",
                 "euler_bernoulli[m]", "ratio_timoshenko[-]", "ratio_euler_bernoulli[-]",
                 "force_balance[-]"});
  json::Value records = json::Value::make_array();
  Scalar judged = 0.0;       // |ratio - 1| of 20 x 1 x 1, nu = 0, incompatible
  Scalar worst_fine = 0.0;   // the same over nx = 20, 40 and both nu
  Scalar worst_sweep = 0.0;  // over L/h = 10, 100, 1000 (20 x 1 x 1, nu = 0)
  Scalar standard_20 = 0.0;
  Scalar standard_1000 = 0.0;
  Scalar balance = 0.0;
  Scalar worst_balance = 0.0;
  const auto record = [&](ElementFormulation f, const Cantilever& c, Index nx, Scalar tip) {
    worst_balance = std::max(worst_balance, balance);
    const Scalar rt = tip / c.timoshenko();
    const Scalar re = tip / c.euler_bernoulli();
    csv.raw_row({short_name(f), fmt(c.poisson, 3), fmt(c.length / c.height, 6),
                 fmt(static_cast<Scalar>(nx), 4), fmt(tip, 10), fmt(c.timoshenko(), 10),
                 fmt(c.euler_bernoulli(), 10), fmt(rt, 8), fmt(re, 8), fmt(balance, 3)});
    json::Value r = json::Value::make_object();
    r.set("formulation", str(short_name(f)));
    r.set("nu", num(c.poisson));
    r.set("slenderness", num(c.length / c.height));
    r.set("nx", num(static_cast<Scalar>(nx)));
    r.set("tip_m", num(tip));
    r.set("ratio_to_timoshenko", num(rt));
    r.set("ratio_to_euler_bernoulli", num(re));
    records.push_back(r);
    return rt;
  };
  for (const ElementFormulation f :
       {ElementFormulation::Standard, ElementFormulation::IncompatibleModes}) {
    for (const Scalar nu : {0.0, 0.3}) {
      Cantilever c;
      c.poisson = nu;
      for (const Index nx : {Index{10}, Index{20}, Index{40}}) {
        const Scalar rt = record(f, c, nx, cantilever_tip(c, nx, 1, 1, f, balance));
        if (f == ElementFormulation::IncompatibleModes) {
          if (nx == 20 && nu == 0.0) judged = std::abs(rt - 1.0);
          if (nx >= 20) worst_fine = std::max(worst_fine, std::abs(rt - 1.0));
        } else if (nx == 20 && nu == 0.0) {
          standard_20 = rt;
        }
      }
    }
    for (const Scalar slenderness : {10.0, 100.0, 1000.0}) {
      Cantilever c;
      c.length = slenderness * c.height;
      const Scalar rt = record(f, c, 20, cantilever_tip(c, 20, 1, 1, f, balance));
      if (f == ElementFormulation::IncompatibleModes) {
        worst_sweep = std::max(worst_sweep, std::abs(rt - 1.0));
      } else if (slenderness == 1000.0) {
        standard_1000 = rt;
      }
    }
  }
  csv.close();
  json::Value block = json::Value::make_object();
  block.set("kind", str("validation (Timoshenko beam theory)"));
  block.set("records", records);
  block.set("incompatible_20x1x1_error", num(judged));
  block.set("incompatible_worst_error_nx_ge_20", num(worst_fine));
  block.set("incompatible_worst_error_slenderness_sweep", num(worst_sweep));
  block.set("standard_20x1x1_ratio", num(standard_20));
  block.set("standard_slenderness_1000_ratio", num(standard_1000));
  block.set("worst_relative_force_balance", num(worst_balance));
  block.set("note",
            str("L = 1 m, h = 0.1 m (y, one element), b = 0.05 m (z, one element), E = 70 "
                "GPa, clamped root face, -1000 N tip resultant on the tip nodes; the mean "
                "tip deflection over Timoshenko's (shear factor 5/6). The standard Hex8 on "
                "one element through the depth locks in shear; the incompatible modes "
                "represent the linear bending stress in each element exactly, so their "
                "error is the discretisation of the moment along x (and the 3-D root) and "
                "does not grow with L/h. With nu = 0.3 the clamped root face also "
                "restrains the Poisson contraction, which beam theory ignores."));
  summary.set("im_cantilever", block);
  StudyOutcome o;
  o.name = "slender cantilever, one incompatible-mode Hex8 through the depth";
  o.kind = "validation";
  o.metric = "|tip / Timoshenko - 1| on 20 x 1 x 1, nu = 0 (and every nx >= 20, nu, L/h <= 5 %)";
  o.value = judged;
  o.tolerance = 0.02;
  o.passed = judged <= o.tolerance && worst_fine <= 0.05 && worst_sweep <= 0.05;
  std::ostringstream note;
  note << "standard Hex8 on 20 x 1 x 1: " << app::format(standard_20, 4)
       << " of Timoshenko (L/h = 10), " << app::format(standard_1000, 4)
       << " at L/h = 1000; incompatible worst over the L/h sweep "
       << app::format(worst_sweep, 3);
  o.note = note.str();
  return o;
}

namespace {

// ---------------------------------------------------------------------------
// MacNeal-Harder straight beam
// ---------------------------------------------------------------------------

enum class CellShape { Rectangular, Parallelogram, Trapezoidal };

const char* shape_name(CellShape s) {
  switch (s) {
    case CellShape::Rectangular: return "rectangular";
    case CellShape::Parallelogram: return "parallelogram";
    case CellShape::Trapezoidal: return "trapezoidal";
  }
  return "?";
}

/// Tip deflection of the MacNeal-Harder straight beam (6 x 0.2 x 0.1, 6
/// cells) under a unit tip load along `component` (1: in plane, along the
/// 0.2 depth; 2: out of plane, along the 0.1 thickness).
Scalar macneal_harder_tip(CellShape shape, ElementFormulation formulation, int component,
                          Index refine = 1) {
  const Scalar length = 6.0;
  const Scalar depth = 0.2;
  const Scalar thickness = 0.1;
  StructuredMeshSpec ms;
  ms.nx = 6 * refine;
  ms.ny = refine;
  ms.nz = refine;
  ms.lx = length;
  ms.ly = depth;
  ms.lz = thickness;
  const Mesh straight = make_structured_hex_mesh(ms);
  Matrix x = straight.coordinates();
  // Interior interfaces inclined at 45 deg in the plane of the depth: the
  // node at y = 0 moves back by depth/2, the one at y = depth forward, all
  // the same way (parallelogram) or alternately (trapezoid).
  for (Index n = 0; n < straight.num_nodes(); ++n) {
    const Scalar xi = x(0, n);
    const Index i = static_cast<Index>(std::lround(xi));
    if (shape == CellShape::Parallelogram) {
      x(0, n) = xi + (x(1, n) / depth - 0.5) * depth;
      continue;
    }
    if (shape == CellShape::Trapezoidal && std::abs(xi - static_cast<Scalar>(i)) > kEps) {
      throw ModelError("verification: the trapezoidal beam is built on 6 cells only");
    }
    // The parallelogram beam has inclined ends too (every cell the same
    // parallelogram); the trapezoidal one keeps its ends square.
    if (shape == CellShape::Rectangular) continue;
    if ((i == 0 || i == 6) && shape == CellShape::Trapezoidal) continue;
    const Scalar sign = shape == CellShape::Parallelogram ? 1.0 : (i % 2 == 0 ? 1.0 : -1.0);
    const Scalar offset = (x(1, n) / depth - 0.5) * depth * sign;
    x(0, n) = xi + offset;
  }
  Mesh mesh(std::move(x), straight.connectivity(), ElementType::Hex8);
  mesh.validate();
  FemModel model(std::move(mesh), IsotropicMaterial(1.0e7, 0.3, 1.0, "macneal"), 1.0,
                 StressState::ThreeDimensional, options_of(formulation));
  // The root and tip faces (inclined for the parallelogram) hold the nodes
  // within depth/2 of x = 0 and x = L.
  const Scalar big = 1.0e9;
  const Scalar reach = 0.5 * depth + kEps;
  model.constraints().push_back(
      fixed(box_group("root", Vector3(-big, -big, -big), Vector3(reach, big, big)), true, true, true));
  LoadCaseSpec lc;
  lc.name = "tip";
  PointLoadSpec tip;
  tip.region = box_group("tip", Vector3(length - reach, -big, -big), Vector3(big, big, big));
  tip.force = Vector3::Zero();
  tip.force(component) = 1.0;
  tip.distribute_total = true;
  lc.point_loads.push_back(tip);
  model.load_case_specs().push_back(lc);
  model.finalize();
  const Vector u = solve_static(model);
  return mean_displacement(model, u, Vector3(length - reach, -big, -big), Vector3(big, big, big),
                           component);
}

}  // namespace

StudyOutcome study_im_macneal_harder(const std::string& out_dir, json::Value& summary) {
  CsvWriter csv(path_join(out_dir, "im_macneal_harder.csv"),
                {"formulation", "cells", "load", "tip[m]", "reference[m]", "normalised[-]"});
  json::Value records = json::Value::make_array();
  const Scalar reference[3] = {0.0, 0.1081, 0.4321};
  Scalar judged = 0.0;  // worst |normalised - 1| of incompatible, rectangular and parallelogram
  Scalar trapezoid_in_plane = 0.0;
  Scalar parallelogram_out_of_plane = 0.0;
  for (const ElementFormulation f :
       {ElementFormulation::Standard, ElementFormulation::IncompatibleModes}) {
    for (const CellShape s :
         {CellShape::Rectangular, CellShape::Parallelogram, CellShape::Trapezoidal}) {
      for (const int component : {1, 2}) {
        const Scalar tip = macneal_harder_tip(s, f, component);
        const Scalar normalised = tip / reference[component];
        csv.raw_row({short_name(f), shape_name(s), component == 1 ? "in-plane" : "out-of-plane",
                     fmt(tip, 8), fmt(reference[component], 6), fmt(normalised, 6)});
        json::Value r = json::Value::make_object();
        r.set("formulation", str(short_name(f)));
        r.set("cells", str(shape_name(s)));
        r.set("load", str(component == 1 ? "in-plane" : "out-of-plane"));
        r.set("normalised_tip", num(normalised));
        records.push_back(r);
        if (f == ElementFormulation::IncompatibleModes) {
          if (s == CellShape::Rectangular || (s == CellShape::Parallelogram && component == 1)) {
            judged = std::max(judged, std::abs(normalised - 1.0));
          } else if (s == CellShape::Trapezoidal && component == 1) {
            trapezoid_in_plane = normalised;
          } else if (s == CellShape::Parallelogram) {
            parallelogram_out_of_plane = normalised;
          }
        }
      }
    }
  }
  // The parallelogram beam's own converged answer: its inclined clamped
  // root and tip make it a different body from the beam of the references.
  Scalar converged[3] = {0.0, 0.0, 0.0};
  for (const int component : {1, 2}) {
    converged[component] = macneal_harder_tip(CellShape::Parallelogram,
                                               ElementFormulation::IncompatibleModes, component,
                                               8) /
                           reference[component];
    csv.raw_row({"incompatible 48 x 8 x 8", "parallelogram",
                 component == 1 ? "in-plane" : "out-of-plane",
                 fmt(converged[component] * reference[component], 8),
                 fmt(reference[component], 6), fmt(converged[component], 6)});
  }
  csv.close();
  json::Value block = json::Value::make_object();
  block.set("kind", str("validation (MacNeal-Harder reference values)"));
  block.set("records", records);
  block.set("incompatible_parallelogram_out_of_plane", num(parallelogram_out_of_plane));
  block.set("parallelogram_fine_mesh_in_plane", num(converged[1]));
  block.set("parallelogram_fine_mesh_out_of_plane", num(converged[2]));
  block.set("incompatible_trapezoidal_in_plane", num(trapezoid_in_plane));
  block.set("note",
            str("MacNeal-Harder straight cantilever, 6 x 0.2 x 0.1, E = 1e7, nu = 0.3, 6 x 1 "
                "x 1 Hex8, clamped root face, unit tip resultant on the tip nodes; mean tip "
                "deflection over 0.1081 (in plane) and 0.4321 (out of plane). The "
                "parallelogram beam is sheared by 45 deg in its depth plane, root and tip "
                "faces included. Judged: the incompatible modes on rectangular cells (both "
                "loads) and on parallelogram cells in plane (affine cells, on which the "
                "Taylor-corrected modes are exact in pure bending). Reported: the "
                "parallelogram out of plane, where its inclined clamped root (bending "
                "coupled with twist) is a boundary layer one element cannot resolve - the "
                "same beam on 48 x 8 x 8 cells is itself stiffer than the reference - and "
                "trapezoidal cells, which are not affine and on which the modes lock "
                "partly (MacNeal 1987)."));
  summary.set("im_macneal_harder", block);
  StudyOutcome o;
  o.name = "MacNeal-Harder straight beam, incompatible-mode Hex8";
  o.kind = "validation";
  o.metric = "worst |tip / reference - 1|: rectangular cells, both loads; parallelogram, in plane";
  o.value = judged;
  o.tolerance = 0.03;
  o.passed = judged <= o.tolerance;
  o.note = "reported only: parallelogram out of plane " +
           app::format(parallelogram_out_of_plane, 4) + " (48 x 8 x 8: " +
           app::format(converged[2], 4) + "), trapezoidal in plane " +
           app::format(trapezoid_in_plane, 4) + " of the reference";
  return o;
}

namespace {

// ---------------------------------------------------------------------------
// Thin square plate, two element layers
// ---------------------------------------------------------------------------

/// Centre deflection of a square plate (side a, thickness t, uniform
/// pressure q on the top face) on the quarter [0, a/2]^2 with n x n x 2
/// cells; clamped (every node of the outer side faces fixed) or simply
/// supported (u_z = 0 on the outer edges' mid-surface nodes).
Scalar plate_centre(Scalar a, Scalar t, Scalar q, Index n, bool clamped,
                    ElementFormulation formulation, Index& dofs) {
  StructuredMeshSpec ms;
  ms.nx = n;
  ms.ny = n;
  ms.nz = 2;
  ms.lx = 0.5 * a;
  ms.ly = 0.5 * a;
  ms.lz = t;
  FemModel model(make_structured_hex_mesh(ms), IsotropicMaterial(200.0e9, 0.3, 7800.0, "steel"),
                 1.0, StressState::ThreeDimensional, options_of(formulation));
  const Scalar big = 1.0e9;
  const Scalar h = 0.5 * a;
  model.constraints().push_back(
      fixed(box_group("sym_x", Vector3(-big, -big, -big), Vector3(kEps, big, big)), true, false,
            false));
  model.constraints().push_back(
      fixed(box_group("sym_y", Vector3(-big, -big, -big), Vector3(big, kEps, big)), false, true,
            false));
  if (clamped) {
    model.constraints().push_back(fixed(
        box_group("edge_x", Vector3(h - kEps, -big, -big), Vector3(big, big, big)), true, true, true));
    model.constraints().push_back(fixed(
        box_group("edge_y", Vector3(-big, h - kEps, -big), Vector3(big, big, big)), true, true, true));
  } else {
    // u_z = 0 over the side faces: they rotate freely about the edge.
    model.constraints().push_back(fixed(
        box_group("edge_x", Vector3(h - kEps, -big, -big), Vector3(big, big, big)), false, false, true));
    model.constraints().push_back(fixed(
        box_group("edge_y", Vector3(-big, h - kEps, -big), Vector3(big, big, big)), false, false, true));
  }
  LoadCaseSpec lc;
  lc.name = "pressure";
  PressureLoadSpec p;
  p.region = box_group("top", Vector3(-big, -big, t - kEps), Vector3(big, big, big));
  p.pressure = q;
  lc.pressures.push_back(p);
  model.load_case_specs().push_back(lc);
  model.finalize();
  dofs = model.dofs().num_dofs();
  const Vector u = solve_static(model);
  return -mean_displacement(model, u, Vector3(-big, -big, -big), Vector3(kEps, kEps, big), 2);
}

}  // namespace

StudyOutcome study_im_plate(const std::string& out_dir, json::Value& summary) {
  const Scalar a = 1.0;
  const Scalar t = a / 50.0;
  const Scalar q = 1.0e4;
  const Scalar rigidity = 200.0e9 * t * t * t / (12.0 * (1.0 - 0.3 * 0.3));
  // Kirchhoff centre deflections of the square plate under a uniform load
  // (Timoshenko and Woinowsky-Krieger 1959, tables 8 and 35).
  const Scalar kirchhoff_ss = 0.00406235 * q * a * a * a * a / rigidity;
  const Scalar kirchhoff_cc = 0.00126532 * q * a * a * a * a / rigidity;
  CsvWriter csv(path_join(out_dir, "im_plate.csv"),
                {"formulation", "support", "n", "num_dofs", "centre[m]", "kirchhoff[m]",
                 "ratio[-]"});
  json::Value records = json::Value::make_array();
  Scalar judged = 0.0;
  Scalar standard_16_cc = 0.0;
  for (const ElementFormulation f :
       {ElementFormulation::Standard, ElementFormulation::IncompatibleModes}) {
    for (const bool clamped : {true, false}) {
      for (const Index n : {Index{4}, Index{8}, Index{16}}) {
        Index dofs = 0;
        const Scalar w = plate_centre(a, t, q, n, clamped, f, dofs);
        const Scalar reference = clamped ? kirchhoff_cc : kirchhoff_ss;
        const Scalar ratio = w / reference;
        csv.raw_row({short_name(f), clamped ? "clamped" : "simply_supported",
                     fmt(static_cast<Scalar>(n), 4), fmt(static_cast<Scalar>(dofs), 8),
                     fmt(w, 10), fmt(reference, 10), fmt(ratio, 6)});
        json::Value r = json::Value::make_object();
        r.set("formulation", str(short_name(f)));
        r.set("support", str(clamped ? "clamped" : "simply_supported"));
        r.set("n", num(static_cast<Scalar>(n)));
        r.set("ratio_to_kirchhoff", num(ratio));
        records.push_back(r);
        if (n == 16) {
          if (f == ElementFormulation::IncompatibleModes) {
            judged = std::max(judged, std::abs(ratio - 1.0));
          } else if (clamped) {
            standard_16_cc = ratio;
          }
        }
      }
    }
  }
  csv.close();
  json::Value block = json::Value::make_object();
  block.set("kind", str("validation (Kirchhoff plate theory)"));
  block.set("records", records);
  block.set("note",
            str("Square plate a = 1 m, t = 0.02 m (t/a = 1/50), E = 200 GPa, nu = 0.3, "
                "uniform pressure 1e4 Pa on the top face, quarter model n x n x 2 Hex8 "
                "(two layers through the thickness), clamped side faces or u_z = 0 on the "
                "mid-surface nodes of the edges; mean centre deflection over Kirchhoff's "
                "0.00126 q a^4/D (clamped) and 0.00406 q a^4/D (simply supported). The "
                "3-D solid differs from Kirchhoff by its transverse shear and the 3-D "
                "support, O((t/a)^2) - a fraction of a percent here."));
  summary.set("im_plate", block);
  StudyOutcome o;
  o.name = "thin square plate (t/a = 1/50), two incompatible-mode Hex8 layers";
  o.kind = "validation";
  o.metric = "worst |centre / Kirchhoff - 1| on 16 x 16 x 2, clamped and simply supported";
  o.value = judged;
  o.tolerance = 0.03;
  o.passed = judged <= o.tolerance;
  o.note = "standard Hex8, clamped, 16 x 16 x 2: " + app::format(standard_16_cc, 4) +
           " of Kirchhoff";
  return o;
}

namespace {

// ---------------------------------------------------------------------------
// Elastoplastic bending and springback of a strip
// ---------------------------------------------------------------------------

struct StripResult {
  Scalar loaded_curvature = 0.0;    ///< [1/m]
  Scalar released_curvature = 0.0;  ///< [1/m]
  Scalar max_plastic_strain = 0.0;
  Index dofs = 0;
  int local_iterations = 0;
  Scalar seconds = 0.0;
  bool completed = false;
};

/// The end-face rotation of a strip from its nodes at x = x_end: the
/// least-squares slope du_x/dz [rad].
Scalar face_rotation(const Mesh& mesh, const Vector& u, Scalar x_end) {
  Scalar sz = 0.0;
  Scalar su = 0.0;
  Scalar szz = 0.0;
  Scalar szu = 0.0;
  Scalar count = 0.0;
  for (Index n = 0; n < mesh.num_nodes(); ++n) {
    const Vector3 x = mesh.node(n);
    if (std::abs(x.x() - x_end) > kEps) continue;
    sz += x.z();
    su += u(3 * n);
    szz += x.z() * x.z();
    szu += x.z() * u(3 * n);
    count += 1.0;
  }
  return (count * szu - sz * su) / (count * szz - sz * sz);
}

/// The strip L x w x t bent to curvature kappa by end rotations (u_x =
/// -kappa (x - L/2)(z - t/2) on both end faces), then released onto
/// statically determinate supports; small strain.
StripResult bend_and_release(const IsotropicMaterial& material, Scalar length, Scalar width,
                             Scalar thickness, Scalar kappa, Index nx, Index ny, Index nz,
                             const IntegrationOptions& integration) {
  StructuredMeshSpec ms;
  ms.nx = nx;
  ms.ny = ny;
  ms.nz = nz;
  ms.lx = length;
  ms.ly = width;
  ms.lz = thickness;
  FemModel model(make_structured_hex_mesh(ms), material, 1.0, StressState::ThreeDimensional,
                 integration);
  LoadCaseSpec lc;
  lc.name = "none";
  lc.prescribed_displacement_only = true;
  model.load_case_specs().push_back(lc);
  model.finalize();
  Assembler assembler(model);
  const Mesh& mesh = model.mesh();
  const auto corner = [&](Scalar x, Scalar y, Scalar z) {
    for (Index n = 0; n < mesh.num_nodes(); ++n) {
      if ((mesh.node(n) - Vector3(x, y, z)).cwiseAbs().maxCoeff() <= kEps) return n;
    }
    throw ModelError("verification: strip corner node not found");
  };
  FormingStep bend;
  bend.name = "bend";
  bend.type = FormingStep::Type::Form;
  bend.increments = 8;
  for (Index n = 0; n < mesh.num_nodes(); ++n) {
    const Vector3 x = mesh.node(n);
    const bool end = x.x() <= kEps || x.x() >= length - kEps;
    if (!end) continue;
    StepConstraint c;
    c.name = "end_" + std::to_string(n);
    c.mode = StepConstraint::Mode::Absolute;
    c.constraint.region = node_group(c.name, {n});
    c.constraint.fix_x = true;
    c.constraint.value_x = -kappa * (x.x() - 0.5 * length) * (x.z() - 0.5 * thickness);
    // The end faces' corners stay at their height (the pure-bending field is
    // even in y - w/2 and z - t/2), and one corner fixes u_y.
    const bool y_edge = x.y() <= kEps || x.y() >= width - kEps;
    const bool z_edge = x.z() <= kEps || x.z() >= thickness - kEps;
    if (y_edge && z_edge) c.constraint.fix_z = true;
    if (x.x() <= kEps && x.y() <= kEps && x.z() <= kEps) c.constraint.fix_y = true;
    bend.constraints.push_back(c);
  }
  FormingStep release;
  release.name = "release";
  release.type = FormingStep::Type::Release;
  release.increments = 8;
  const auto hold = [&](Index n, bool x, bool y, bool z) {
    StepConstraint c;
    c.name = "support_" + std::to_string(n);
    c.mode = StepConstraint::Mode::Hold;
    c.constraint.region = node_group(c.name, {n});
    c.constraint.fix_x = x;
    c.constraint.fix_y = y;
    c.constraint.fix_z = z;
    return c;
  };
  release.constraints.push_back(hold(corner(0.0, 0.0, 0.0), true, true, true));
  release.constraints.push_back(hold(corner(0.0, width, 0.0), true, false, true));
  release.constraints.push_back(hold(corner(length, 0.0, 0.0), false, false, true));
  FormingOptions fo;
  fo.kinematics = Kinematics::SmallStrain;
  fo.steps = {bend, release};
  fo.residual_tolerance = 1.0e-9;
  fo.displacement_tolerance = 1.0e-9;
  StripResult out;
  out.dofs = model.dofs().num_dofs();
  Timer timer;
  const FormingResult r = FormingAnalysis(model, assembler, fo).run();
  out.seconds = timer.elapsed_seconds();
  out.completed = r.completed;
  if (!r.completed || r.steps.size() != 2) return out;
  const auto curvature = [&](const Vector& u) {
    return (face_rotation(mesh, u, 0.0) - face_rotation(mesh, u, length)) / length;
  };
  out.loaded_curvature = curvature(r.steps[0].displacement);
  out.released_curvature = curvature(r.steps[1].displacement);
  out.max_plastic_strain = r.steps[0].max_plastic_strain;
  return out;
}

}  // namespace

StudyOutcome study_im_springback(const std::string& out_dir, json::Value& summary) {
  // A 2 mm steel sheet strip, 2 mm wide, 20 mm long, bent to 4 times the
  // curvature at first yield of its surface.
  const Scalar length = 0.02;
  const Scalar width = 0.002;
  const Scalar thickness = 0.002;
  const Scalar youngs = 200.0e9;
  const Scalar sy = 250.0e6;
  IsotropicMaterial material(youngs, 0.3, 7800.0, "steel");
  PlasticityParameters p;
  p.yield_stress = sy;
  p.hardening_modulus = 1.0e9;
  material.set_plasticity(p);
  const Scalar kappa_y = 2.0 * sy / (youngs * thickness);
  const Scalar kappa = 4.0 * kappa_y;

  struct Variant {
    std::string name;
    ElementFormulation formulation;
    int thickness_points;
    Index nx, ny, nz;
    bool reference;
  };
  const std::vector<Variant> variants = {
      {"incompatible, 1 layer, 2x2x5", ElementFormulation::IncompatibleModes, 5, 10, 1, 1, false},
      {"incompatible, 1 layer, 2x2x7", ElementFormulation::IncompatibleModes, 7, 10, 1, 1, false},
      {"incompatible, 1 layer, 2x2x2", ElementFormulation::IncompatibleModes, 0, 10, 1, 1, false},
      {"incompatible, 2 layers, 2x2x3", ElementFormulation::IncompatibleModes, 3, 10, 1, 2, false},
      {"standard, 1 layer, 2x2x5", ElementFormulation::Standard, 5, 10, 1, 1, false},
      {"standard, 4 layers, 2x2x2", ElementFormulation::Standard, 0, 40, 4, 4, true},
      {"standard, 8 layers, 2x2x2", ElementFormulation::Standard, 0, 80, 8, 8, true},
  };
  CsvWriter csv(path_join(out_dir, "im_springback.csv"),
                {"variant", "num_dofs", "loaded_curvature[1/m]", "released_curvature[1/m]",
                 "springback[1/m]", "springback_ratio[-]", "max_plastic_strain[-]",
                 "seconds[s]"});
  json::Value records = json::Value::make_array();
  std::vector<StripResult> results;
  for (const Variant& v : variants) {
    const StripResult r =
        bend_and_release(material, length, width, thickness, kappa, v.nx, v.ny, v.nz,
                         options_of(v.formulation, v.thickness_points, 2));
    results.push_back(r);
    const Scalar springback = r.loaded_curvature - r.released_curvature;
    csv.raw_row({v.name, fmt(static_cast<Scalar>(r.dofs), 8), fmt(r.loaded_curvature, 10),
                 fmt(r.released_curvature, 10), fmt(springback, 10),
                 fmt(springback / r.loaded_curvature, 8), fmt(r.max_plastic_strain, 6),
                 fmt(r.seconds, 4)});
    json::Value rec = json::Value::make_object();
    rec.set("variant", str(v.name));
    rec.set("completed", json::Value::make_bool(r.completed));
    rec.set("num_dofs", num(static_cast<Scalar>(r.dofs)));
    rec.set("loaded_curvature_per_m", num(r.loaded_curvature));
    rec.set("released_curvature_per_m", num(r.released_curvature));
    rec.set("springback_per_m", num(springback));
    rec.set("seconds", num(r.seconds));
    records.push_back(rec);
  }
  csv.close();
  // The reference: Richardson extrapolation of the 4- and 8-layer standard
  // meshes (second order), and its distance to the 8-layer value as the
  // reference's own uncertainty.
  const auto springback_of = [&](std::size_t i) {
    return results[i].loaded_curvature - results[i].released_curvature;
  };
  const Scalar s4 = springback_of(5);
  const Scalar s8 = springback_of(6);
  const Scalar reference = s8 + (s8 - s4) / 3.0;
  const Scalar uncertainty = std::abs(reference - s8) / std::abs(reference);
  const auto error_of = [&](std::size_t i) {
    return results[i].completed ? std::abs(springback_of(i) - reference) / std::abs(reference)
                                : 1.0;
  };
  // The quadrature error of the moment alone: sigma(z) of a rigid-plastic
  // shell outside the elastic core |z| < t/8 (kappa = 4 kappa_y, the linear
  // hardening neglected), integrated by the n-point Gauss rule against
  // exactly; the springback is M / (E I), so it carries this error.
  const auto moment_quadrature_error = [&](int n) {
    const Scalar core = kappa_y / kappa;  // the elastic core, relative to t/2
    Scalar gauss = 0.0;
    for (const QuadraturePoint1D& p : gauss_legendre_line_extended(n)) {
      gauss += p.weight * p.xi * std::clamp(p.xi / core, -1.0, 1.0);
    }
    const Scalar exact = 2.0 * (core * core / 3.0 + 0.5 * (1.0 - core * core));
    return std::abs(gauss - exact) / exact;
  };
  const Scalar judged = error_of(0);
  const Scalar seven = error_of(1);
  const Scalar standard_one_layer = error_of(4);
  bool completed = true;
  for (const StripResult& r : results) completed = completed && r.completed;
  json::Value block = json::Value::make_object();
  block.set("kind", str("verification (converged fine mesh)"));
  block.set("records", records);
  block.set("reference_springback_per_m", num(reference));
  block.set("reference_uncertainty", num(uncertainty));
  block.set("error_incompatible_1_layer_2x2x5", num(judged));
  block.set("error_incompatible_1_layer_2x2x7", num(seven));
  block.set("error_incompatible_1_layer_2x2x2", num(error_of(2)));
  block.set("error_incompatible_2_layers_2x2x3", num(error_of(3)));
  block.set("error_standard_1_layer_2x2x5", num(standard_one_layer));
  block.set("moment_quadrature_error_5_points", num(moment_quadrature_error(5)));
  block.set("moment_quadrature_error_7_points", num(moment_quadrature_error(7)));
  block.set("note",
            str("A steel strip 20 x 2 x 2 mm (E = 200 GPa, nu = 0.3, sigma_y = 250 MPa, "
                "linear hardening 1 GPa), small strain, bent by end rotations u_x = -k (x - "
                "L/2)(z - t/2) to k = 4 k_y (k_y = 2 sigma_y/(E t)), then released onto "
                "statically determinate supports (the forming driver's release). The "
                "curvature is the relative rotation of the end faces (least-squares du_x/dz) "
                "over L; the springback is its change on release. Reference: the standard "
                "Hex8 on 4 and 8 layers of cubic cells, Richardson-extrapolated. One "
                "incompatible-mode layer represents the linear strain through the "
                "thickness and unloads exactly (M / (E I)); what is left is the Gauss "
                "quadrature of the kinked elastic-plastic stress through the thickness, "
                "estimated by moment_quadrature_error_n (the rigid-plastic shell outside "
                "the elastic core |z| < t/8, integrated by n points against exactly)."));
  summary.set("im_springback", block);
  StudyOutcome o;
  o.name = "elastoplastic bending and springback, one incompatible-mode layer, 2x2x5 rule";
  o.kind = "verification";
  o.metric = "relative springback error vs the extrapolated fine standard mesh (2x2x7 <= 2 %)";
  o.value = judged;
  o.tolerance = 0.05;
  o.passed = completed && judged <= o.tolerance && seven <= 0.02 && uncertainty <= 0.01 &&
             judged < 0.5 * standard_one_layer;
  std::ostringstream note;
  note << "2x2x7: " << app::format(seven, 3) << "; 1 layer 2x2x2: "
       << app::format(error_of(2), 3) << "; standard 1 layer 2x2x5: "
       << app::format(standard_one_layer, 3) << "; quadrature estimate 5 / 7 points "
       << app::format(moment_quadrature_error(5), 3) << " / "
       << app::format(moment_quadrature_error(7), 3) << "; reference uncertainty "
       << app::format(uncertainty, 3);
  o.note = note.str();
  return o;
}

}  // namespace verify
}  // namespace sparlab
