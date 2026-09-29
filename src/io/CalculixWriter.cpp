#include "sparlab/io/CalculixWriter.hpp"

#include "sparlab/core/Exceptions.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

namespace sparlab {
namespace {

std::string sanitise(const std::string& name) {
  std::string out;
  for (char c : name) {
    const bool ok = (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
                    (c >= '0' && c <= '9') || c == '-' || c == '_';
    out.push_back(ok ? c : '_');
  }
  return out.empty() ? "unnamed" : out;
}

/// CalculiX reads at most 20 characters per data field and rejects a longer
/// one ("1.4408030432795649e-18" is 22): the most significant digits of `v`
/// that fit - 14 for a two-digit exponent, 13 at the least for any double.
std::string field(Scalar v) {
  for (int digits = 17; digits >= 6; --digits) {
    std::ostringstream os;
    os << std::setprecision(digits) << v;
    if (os.str().size() <= 20) return os.str();
  }
  std::ostringstream os;
  os << std::scientific << std::setprecision(12) << v;
  return os.str();
}

/// Corner nodes (local, 0-based) of CalculiX's faces 1, 2, ... of an element:
/// edges 1-2, 2-3, ... of the plane elements; for C3D8 the faces 1-2-3-4,
/// 5-8-7-6, 1-5-6-2, 2-6-7-3, 3-7-8-4, 4-8-5-1; for C3D4 and C3D10 the faces
/// 1-2-3, 1-4-2, 2-4-3, 3-4-1 (CalculiX manual, section 6.2).
const std::vector<std::vector<int>>& calculix_faces(ElementType type) {
  static const std::vector<std::vector<int>> quad = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
  static const std::vector<std::vector<int>> tri = {{0, 1}, {1, 2}, {2, 0}};
  static const std::vector<std::vector<int>> hex = {{0, 1, 2, 3}, {4, 7, 6, 5}, {0, 4, 5, 1},
                                                    {1, 5, 6, 2}, {2, 6, 7, 3}, {3, 7, 4, 0}};
  static const std::vector<std::vector<int>> tet = {{0, 1, 2}, {0, 3, 1}, {1, 3, 2}, {2, 3, 0}};
  switch (type) {
    case ElementType::Quad4: return quad;
    case ElementType::Tri3: return tri;
    case ElementType::Hex8: return hex;
    case ElementType::Tet4:
    case ElementType::Tet10: return tet;
  }
  throw IoError("no CalculiX face table for this element type");
}

/// CalculiX's number (1-based) of SparLab's local face `local_face`: the face
/// with the same corner nodes.
int calculix_face_number(ElementType type, int local_face) {
  const std::vector<int>& ours =
      element_local_faces(type)[static_cast<std::size_t>(local_face)];
  const std::vector<std::vector<int>>& theirs = calculix_faces(type);
  const std::size_t corners = theirs.front().size();
  std::vector<int> key(ours.begin(), ours.begin() + static_cast<std::ptrdiff_t>(corners));
  std::sort(key.begin(), key.end());
  for (std::size_t f = 0; f < theirs.size(); ++f) {
    std::vector<int> candidate = theirs[f];
    std::sort(candidate.begin(), candidate.end());
    if (candidate == key) return static_cast<int>(f) + 1;
  }
  throw IoError("a SparLab face has no CalculiX counterpart");
}

/// CalculiX heat-transfer element for the model's element type. CalculiX 2.21
/// ignores the plane DC2D4 / DC2D3 cards (it reads no integration point for
/// them), but conducts heat in its plane elements, expanded through the
/// thickness like in a mechanical step, when the step is *HEAT TRANSFER.
std::string calculix_conduction_type(const FemModel& model) {
  switch (model.mesh().element_type()) {
    case ElementType::Quad4:
    case ElementType::Tri3: return calculix_element_type(model);
    case ElementType::Hex8: return "DC3D8";
    case ElementType::Tet4: return "DC3D4";
    case ElementType::Tet10: return "DC3D10";
  }
  throw IoError("no CalculiX heat-transfer element for this mesh");
}

/// Write an index list as a set card, at most 16 entries per line.
void write_set(std::ostream& out, const std::string& keyword, const std::string& name,
               const std::vector<Index>& ids) {
  out << "*" << keyword << ", " << (keyword == "NSET" ? "NSET=" : "ELSET=") << name << "\n";
  for (std::size_t i = 0; i < ids.size(); ++i) {
    out << ids[i] + 1 << (i + 1 == ids.size() || (i + 1) % 16 == 0 ? "\n" : ", ");
  }
}

void write_nodes_and_elements(std::ostream& out, const Mesh& mesh, const std::string& type) {
  out << "*NODE, NSET=NALL\n";
  for (Index n = 0; n < mesh.num_nodes(); ++n) {
    const Vector3 x = mesh.node(n);
    out << n + 1 << ", " << field(x.x()) << ", " << field(x.y()) << ", " << field(x.z())
        << "\n";
  }
  out << "*ELEMENT, TYPE=" << type << ", ELSET=EALL\n";
  for (Index e = 0; e < mesh.num_elements(); ++e) {
    const Index* nodes = mesh.element_nodes(e);
    out << e + 1;
    for (int a = 0; a < mesh.nodes_per_elem(); ++a) out << ", " << nodes[a] + 1;
    out << "\n";
  }
}

/// Element sets M1, M2, ... of the model's materials (EALL for one material),
/// in material order.
std::vector<std::string> write_material_sets(std::ostream& out, const FemModel& model) {
  if (model.single_material()) return {"EALL"};
  std::vector<std::vector<Index>> members(static_cast<std::size_t>(model.num_materials()));
  for (Index e = 0; e < model.mesh().num_elements(); ++e) {
    members[static_cast<std::size_t>(model.element_material(e))].push_back(e);
  }
  std::vector<std::string> names;
  for (std::size_t m = 0; m < members.size(); ++m) {
    const std::string name = "M" + std::to_string(m + 1);
    names.push_back(members[m].empty() ? std::string() : name);
    if (!members[m].empty()) write_set(out, "ELSET", name, members[m]);
  }
  return names;
}

/// The boundary faces of `faces_in_region` as (element, CalculiX face) pairs.
std::vector<std::pair<Index, int>> calculix_faces_of(const Mesh& mesh,
                                                     const std::vector<Mesh::BoundaryFace>& boundary,
                                                     const SelectorGroup& region) {
  std::vector<std::pair<Index, int>> out;
  for (const Mesh::BoundaryFace& face : faces_in_region(mesh, boundary, region)) {
    out.emplace_back(face.element,
                     calculix_face_number(mesh.element_type(), face.local_face));
  }
  return out;
}

/// The one stress-free temperature of the model's materials: CalculiX measures
/// thermal strain from the initial nodal temperature, which a node shared by
/// materials of different reference temperatures cannot carry.
Scalar common_reference_temperature(const FemModel& model) {
  const Scalar t_ref = model.material().reference_temperature();
  for (const IsotropicMaterial& m : model.materials()) {
    if (m.thermal_expansion() != 0.0 && m.reference_temperature() != t_ref) {
      throw IoError("CalculiX measures thermal strain from the initial nodal temperature, "
                    "so materials with different reference temperatures (" +
                    field(t_ref) + " and " + field(m.reference_temperature()) +
                    " K) cannot be exported");
    }
  }
  return t_ref;
}

/// The (yield stress, equivalent plastic strain) table of an isotropic
/// hardening law, up to a plastic strain of 10 (CalculiX holds the last
/// yield stress beyond it): exact for linear hardening; for Voce saturation
/// chords within 1e-4 of the saturation stress Q, since the chord over h
/// departs from the curve by at most h^2 |sigma_y''| / 8 with
/// |sigma_y''| = Q delta^2 e^(-delta alpha) largest at its start.
std::vector<std::pair<Scalar, Scalar>> hardening_table(const PlasticityParameters& p) {
  constexpr Scalar kEnd = 10.0;
  std::vector<std::pair<Scalar, Scalar>> out;
  out.emplace_back(p.yield(0.0), 0.0);
  if (p.saturation_stress == 0.0) {
    if (p.hardening_modulus > 0.0) out.emplace_back(p.yield(kEnd), kEnd);
    return out;
  }
  const Scalar q = p.saturation_stress;
  const Scalar d = p.saturation_rate;
  Scalar a = 0.0;
  while (a < kEnd) {
    const Scalar curvature = q * d * d * std::exp(-d * a);
    a = std::min(kEnd, a + std::sqrt(8.0 * 1.0e-4 * q / curvature));
    out.emplace_back(p.yield(a), a);
    if (d * a > 40.0) {
      // Saturated to round-off: the rest of the curve is linear.
      if (a < kEnd) out.emplace_back(p.yield(kEnd), kEnd);
      break;
    }
  }
  return out;
}

/// The mechanical deck of load case `l`: linear, or with `nonlinear` a
/// non-linear analysis along its load path.
void write_static_deck(std::ostream& out, const FemModel& model, std::size_t l,
                       const std::string& case_name,
                       const std::vector<Mesh::BoundaryFace>& boundary,
                       const CalculixNonlinearExport* nonlinear = nullptr,
                       const CalculixTransientExport* transient = nullptr) {
  const Mesh& mesh = model.mesh();
  const int dim = mesh.dim();
  // CalculiX numbers the DOFs of a node 1-3 (translations) and 4-6
  // (rotations), the same order as SparLab's.
  const int ndpn = model.dofs_per_node();
  const std::string element = calculix_element_type(model);
  const LoadCaseSpec& spec = model.load_case_specs()[l];
  const LoadCaseData& data = model.load_case_data(l);
  const bool thermal = data.temperature.size() > 0;
  const bool dynamic = transient != nullptr;
  const bool needs_density =
      dynamic || spec.gravity.squaredNorm() > 0.0 || spec.centrifugal.enabled;
  bool plastic = false;
  if (nonlinear != nullptr || (dynamic && transient->options.nonlinear)) {
    for (const IsotropicMaterial& m : model.materials()) plastic = plastic || m.plasticity().enabled();
  }
  const bool damped = dynamic && (transient->options.mass_damping > 0.0 ||
                                  transient->options.stiffness_damping > 0.0);
  const int steps = dynamic ? static_cast<int>(std::lround(transient->options.end_time /
                                                           transient->options.time_step))
                            : 0;

  out << "*HEADING\n";
  out << "SparLab cross-validation export: " << case_name << " / " << spec.name << " ("
      << element << ", " << mesh.num_elements() << " elements)\n";
  write_nodes_and_elements(out, mesh, element);
  const std::vector<std::string> sets = write_material_sets(out, model);
  for (std::size_t m = 0; m < sets.size(); ++m) {
    if (sets[m].empty()) continue;
    const IsotropicMaterial& mat = model.materials()[m];
    out << "*MATERIAL, NAME=MAT" << m + 1 << "\n*ELASTIC\n" << field(mat.youngs_modulus())
        << ", " << field(mat.poisson_ratio()) << "\n";
    if (plastic && mat.plasticity().enabled()) {
      out << "*PLASTIC\n";
      for (const auto& point : hardening_table(mat.plasticity())) {
        out << field(point.first) << ", " << field(point.second) << "\n";
      }
    }
    if (needs_density) out << "*DENSITY\n" << field(mat.density()) << "\n";
    if (damped) {
      // Rayleigh damping of a direct integration, C = a M + b K.
      out << "*DAMPING, ALPHA=" << field(transient->options.mass_damping)
          << ", BETA=" << field(transient->options.stiffness_damping) << "\n";
    }
    if (thermal) {
      out << "*EXPANSION, ZERO=" << field(mat.reference_temperature()) << "\n"
          << field(mat.thermal_expansion()) << "\n";
    }
    out << "*SOLID SECTION, ELSET=" << sets[m] << ", MATERIAL=MAT" << m + 1 << "\n";
    if (dim == 2) out << field(model.thickness()) << "\n";
  }
  // Body-force regions get their own element sets.
  std::vector<std::string> body_sets;
  for (std::size_t b = 0; b < spec.body_forces.size(); ++b) {
    const BodyForceSpec& body = spec.body_forces[b];
    if (body.whole_model) {
      body_sets.push_back("EALL");
      continue;
    }
    const std::string name = "BODY" + std::to_string(b + 1);
    write_set(out, "ELSET", name, body.region.select_elements(mesh));
    body_sets.push_back(name);
  }
  Scalar t_ref = 0.0;
  if (thermal) {
    t_ref = common_reference_temperature(model);
    out << "*INITIAL CONDITIONS, TYPE=TEMPERATURE\nNALL, " << field(t_ref) << "\n";
  }

  // Point loads and tractions as the assembled nodal forces; every other load
  // in CalculiX's own form, so that CalculiX integrates it itself - except a
  // dead pressure under NLGEOM, where CalculiX's face load would follow the
  // face: it stays among the nodal forces of the undeformed faces. (Without
  // NLGEOM a face load acts on the undeformed face.)
  bool finite_step = nonlinear != nullptr && nonlinear->nlgeom;
  bool follower = nonlinear != nullptr && nonlinear->follower_pressure;
  if (dynamic && transient->options.nonlinear) {
    finite_step = transient->options.nonlinear_options.kinematics == Kinematics::Finite;
    follower = transient->options.nonlinear_options.follower_pressure;
  }
  const bool pressure_faces = !spec.pressures.empty() && (!finite_step || follower);
  Vector concentrated = data.mechanical;
  if (pressure_faces) {
    LoadCaseSpec pressures_only;
    pressures_only.name = spec.name;
    pressures_only.pressures = spec.pressures;
    concentrated -= assemble_load_vector(mesh, model.element(), pressures_only,
                                         model.thickness(), model.integration());
  }

  // A transient's amplitude at every step time: the method reads the loads
  // at those times only, so the table is exact where it is used.
  const std::string on_amplitude = dynamic ? ", AMPLITUDE=A1" : "";
  if (dynamic) {
    const TransientOptions& o = transient->options;
    out << "*AMPLITUDE, NAME=A1\n";
    for (int k = 0; k <= steps; ++k) {
      const Scalar t = k * o.time_step;
      out << field(t) << ", " << field(o.amplitude.value(t)) << "\n";
    }
  }

  // One step per leg of the load path, each ramping every load from the last
  // leg's level to its own.
  std::vector<Scalar> levels{1.0};
  if (nonlinear != nullptr && !nonlinear->load_path.empty()) levels = nonlinear->load_path;
  for (const Scalar level : levels) {
    if (dynamic) {
      // CalculiX's HHT-alpha method with fixed steps; a non-linear run
      // converged tightly (1e-6 of the residual and the correction, as the
      // static decks).
      const TransientOptions& o = transient->options;
      const bool finite = o.nonlinear && o.nonlinear_options.kinematics == Kinematics::Finite;
      out << "*STEP" << (finite ? ", NLGEOM" : "") << ", INC=" << 10 * steps << "\n";
      if (o.nonlinear) {
        out << "*CONTROLS, PARAMETERS=FIELD\n1.e-6, 1.e-6\n"
            << "*CONTROLS, PARAMETERS=TIME INCREMENTATION\n20, 30, 200, 200\n";
      }
      out << "*DYNAMIC, DIRECT, ALPHA=" << field(o.alpha) << "\n"
          << field(o.time_step) << ", " << field(steps * o.time_step) << "\n";
    } else if (nonlinear == nullptr) {
      out << "*STEP\n*STATIC\n";
    } else if (plastic || !nonlinear->nlgeom) {
      // The increments of SparLab's run: the return depends on them.
      const int increments = std::max(nonlinear->increments, 1);
      out << "*STEP" << (nonlinear->nlgeom ? ", NLGEOM" : "") << ", INC=" << 10 * increments
          << "\n"
          << "*CONTROLS, PARAMETERS=FIELD\n1.e-6, 1.e-6\n"
          << "*CONTROLS, PARAMETERS=TIME INCREMENTATION\n20, 30, 200, 200\n"
          << "*STATIC, DIRECT\n"
          << field(1.0 / static_cast<Scalar>(increments)) << ", 1.\n";
    } else {
      // Automatic increments from 1 / increments (at most 1/50) down to 1e-6
      // of it; the step time is the load factor. At least 50 increments,
      // because CalculiX lags the deformed-position centrifugal load within
      // an increment (its answer converges to SparLab's as the increments
      // shrink: measured 1e-5 of the displacement at 10, 1.4e-6 at 80).
      // CalculiX's default convergence test (largest residual below 0.005 of
      // the mean force, largest correction below 0.01 of the increment)
      // leaves errors of about 1e-4 of the displacement; a comparison needs
      // it converged: 1e-6 on both, up to 200 iterations, and never the
      // looser residual test it otherwise switches to after the ninth
      // iteration.
      const int increments = std::max(nonlinear->increments, 50);
      const Scalar dt = 1.0 / static_cast<Scalar>(increments);
      out << "*STEP, NLGEOM, INC=" << 100 * increments << "\n"
          << "*CONTROLS, PARAMETERS=FIELD\n1.e-6, 1.e-6\n"
          << "*CONTROLS, PARAMETERS=TIME INCREMENTATION\n20, 30, 200, 200\n"
          << "*STATIC\n"
          << field(dt) << ", 1., " << field(1.0e-6 * dt) << ", " << field(dt) << "\n";
    }
    // Held DOFs on one card; in a transient the prescribed motions follow
    // the amplitude on a second.
    for (const bool moving : {false, true}) {
      bool header = false;
      for (Index d : model.dofs().constrained_dofs()) {
        const Scalar value = model.dofs().prescribed_value(d);
        if (dynamic && moving != (value != 0.0)) continue;
        if (!dynamic && moving) continue;
        if (!header) out << "*BOUNDARY" << (dynamic && moving ? on_amplitude : "") << "\n";
        header = true;
        const Index node = d / ndpn;
        const int component = static_cast<int>(d % ndpn) + 1;
        out << node + 1 << ", " << component << ", " << component << ", "
            << field(level * value) << "\n";
      }
    }
    bool any = false;
    for (Index n = 0; n < mesh.num_nodes(); ++n) {
      for (int k = 0; k < ndpn; ++k) {
        const Scalar f = concentrated(n * ndpn + k);
        if (f == 0.0) continue;
        if (!any) out << "*CLOAD" << on_amplitude << "\n";
        any = true;
        out << n + 1 << ", " << k + 1 << ", " << field(level * f) << "\n";
      }
    }
    const bool distributed = pressure_faces || spec.has_body_loads();
    if (distributed) out << "*DLOAD" << on_amplitude << "\n";
    for (const PressureLoadSpec& p : spec.pressures) {
      if (!pressure_faces) break;
      for (const auto& face : calculix_faces_of(mesh, boundary, p.region)) {
        out << face.first + 1 << ", P" << face.second << ", " << field(level * p.pressure)
            << "\n";
      }
    }
    if (spec.gravity.squaredNorm() > 0.0) {
      const Scalar g = spec.gravity.norm();
      const Vector3 e = spec.gravity / g;
      out << "EALL, GRAV, " << field(level * g) << ", " << field(e.x()) << ", "
          << field(e.y()) << ", " << field(e.z()) << "\n";
    }
    for (std::size_t b = 0; b < spec.body_forces.size(); ++b) {
      const Vector3& density = spec.body_forces[b].force_density;
      static const char* const labels[3] = {"BX", "BY", "BZ"};
      for (int k = 0; k < dim; ++k) {
        if (density(k) != 0.0) {
          out << body_sets[b] << ", " << labels[k] << ", " << field(level * density(k)) << "\n";
        }
      }
    }
    if (spec.centrifugal.enabled) {
      const CentrifugalSpec& c = spec.centrifugal;
      const Vector3 axis = c.axis.normalized();
      out << "EALL, CENTRIF, "
          << field(level * c.angular_velocity * c.angular_velocity) << ", "
          << field(c.point.x()) << ", " << field(c.point.y()) << ", " << field(c.point.z())
          << ", " << field(axis.x()) << ", " << field(axis.y()) << ", " << field(axis.z())
          << "\n";
    }
    if (thermal) {
      out << "*TEMPERATURE\n";
      for (Index n = 0; n < mesh.num_nodes(); ++n) {
        out << n + 1 << ", " << field(t_ref + level * (data.temperature(n) - t_ref)) << "\n";
      }
    }
    if (dynamic) {
      // The displacements at the snapshot increments (only the last without
      // snapshots).
      const int every = transient->options.snapshot_every > 0
                            ? transient->options.snapshot_every
                            : steps;
      out << "*NODE FILE, FREQUENCY=" << every << "\nU\n*END STEP\n";
    } else {
      out << "*NODE FILE\nU\n*EL FILE\nS\n*END STEP\n";
    }
  }
}

/// The steady heat-transfer deck of load case `l`'s conduction problem.
void write_conduction_deck(std::ostream& out, const FemModel& model, std::size_t l,
                           const std::string& case_name,
                           const std::vector<Mesh::BoundaryFace>& boundary) {
  const Mesh& mesh = model.mesh();
  const LoadCaseSpec& spec = model.load_case_specs()[l];
  const LoadCaseData& data = model.load_case_data(l);
  const ConductionSpec& c = spec.temperature.conduction;
  const std::string element = calculix_conduction_type(model);

  out << "*HEADING\n";
  out << "SparLab conduction export: " << case_name << " / " << spec.name << " (" << element
      << ", " << mesh.num_elements() << " elements)\n";
  write_nodes_and_elements(out, mesh, element);
  const std::vector<std::string> sets = write_material_sets(out, model);
  for (std::size_t m = 0; m < sets.size(); ++m) {
    if (sets[m].empty()) continue;
    const IsotropicMaterial& mat = model.materials()[m];
    out << "*MATERIAL, NAME=MAT" << m + 1 << "\n*CONDUCTIVITY\n"
        << field(mat.conductivity()) << "\n";
    out << "*SOLID SECTION, ELSET=" << sets[m] << ", MATERIAL=MAT" << m + 1 << "\n";
    if (mesh.dim() == 2) out << field(model.thickness()) << "\n";
  }
  std::vector<std::string> source_sets;
  for (std::size_t s = 0; s < c.sources.size(); ++s) {
    if (c.sources[s].whole_model) {
      source_sets.push_back("EALL");
      continue;
    }
    const std::string name = "SOURCE" + std::to_string(s + 1);
    write_set(out, "ELSET", name, c.sources[s].region.select_elements(mesh));
    source_sets.push_back(name);
  }
  // The start of CalculiX's (linear, one-increment) iteration: the mean
  // solved temperature, so no value of the answer is handed over.
  out << "*INITIAL CONDITIONS, TYPE=TEMPERATURE\nNALL, " << field(data.temperature.mean())
      << "\n";

  out << "*STEP\n*HEAT TRANSFER, STEADY STATE\n1., 1.\n";
  // Prescribed temperatures: later regions win, as in the solve.
  std::vector<Scalar> prescribed(static_cast<std::size_t>(mesh.num_nodes()),
                                 std::numeric_limits<Scalar>::quiet_NaN());
  for (const RegionValue& p : c.prescribed) {
    for (Index n : p.region.select_nodes(mesh)) prescribed[static_cast<std::size_t>(n)] = p.value;
  }
  if (!c.prescribed.empty()) {
    out << "*BOUNDARY\n";
    for (Index n = 0; n < mesh.num_nodes(); ++n) {
      const Scalar v = prescribed[static_cast<std::size_t>(n)];
      if (!std::isnan(v)) out << n + 1 << ", 11, 11, " << field(v) << "\n";
    }
  }
  if (!c.fluxes.empty() || !c.sources.empty()) out << "*DFLUX\n";
  for (const RegionValue& q : c.fluxes) {
    for (const auto& face : calculix_faces_of(mesh, boundary, q.region)) {
      out << face.first + 1 << ", S" << face.second << ", " << field(q.value) << "\n";
    }
  }
  for (std::size_t s = 0; s < c.sources.size(); ++s) {
    out << source_sets[s] << ", BF, " << field(c.sources[s].value) << "\n";
  }
  if (!c.convection.empty()) out << "*FILM\n";
  for (const ConvectionSpec& h : c.convection) {
    for (const auto& face : calculix_faces_of(mesh, boundary, h.region)) {
      out << face.first + 1 << ", F" << face.second << ", " << field(h.ambient) << ", "
          << field(h.film_coefficient) << "\n";
    }
  }
  out << "*NODE FILE\nNT\n*END STEP\n";
}

}  // namespace

std::string calculix_element_type(const FemModel& model) {
  switch (model.mesh().element_type()) {
    case ElementType::Quad4:
      return model.stress_state() == StressState::PlaneStrain ? "CPE4" : "CPS4";
    case ElementType::Hex8:
      return "C3D8";
    case ElementType::Tri3:
      return model.stress_state() == StressState::PlaneStrain ? "CPE3" : "CPS3";
    case ElementType::Tet4:
      return "C3D4";
    case ElementType::Tet10:
      return "C3D10";  // same node order as the Tet10
  }
  throw IoError("no CalculiX element type for this mesh");
}

std::string calculix_transient_obstacle(const FemModel& model, std::size_t l,
                                        const TransientOptions& options) {
  if (options.mass_type == MassType::Lumped) {
    return "CalculiX's implicit dynamics uses the consistent mass, and the run's is lumped";
  }
  if (options.start == TransientOptions::Start::Static) {
    return "CalculiX's dynamic step starts at rest, and the run starts from the static state";
  }
  if (options.amplitude.value(0.0) != 0.0) {
    return "CalculiX does not start a dynamic step in equilibrium with a load already acting "
           "at t = 0 as SparLab does (measured: 2.6 % apart in the first step of a sudden "
           "load); ramp the amplitude up from zero";
  }
  if (l >= model.load_case_specs().size()) return "no such load case";
  const LoadCaseSpec& spec = model.load_case_specs()[l];
  if (model.load_case_data(l).temperature.size() > 0) return "the load case is thermal";
  if (spec.centrifugal.enabled) return "the load case rotates";
  if (options.nonlinear) {
    const NonlinearOptions& nl = options.nonlinear_options;
    if (nl.kinematics == Kinematics::Finite && nl.law != HyperelasticModel::SaintVenantKirchhoff) {
      return "CalculiX's NEO HOOKE is a different strain energy from SparLab's neo-Hookean law";
    }
    if (nl.kinematics == Kinematics::FiniteLogarithmic) {
      return calculix_logarithmic_obstacle();
    }
    for (const IsotropicMaterial& m : model.materials()) {
      const std::string obstacle = calculix_plasticity_obstacle(m);
      if (!obstacle.empty()) return obstacle;
    }
  }
  return "";
}

std::string calculix_logarithmic_obstacle() {
  return "the logarithmic-strain kinematics (Hencky elasticity, plasticity additive in the log "
         "strain) has no counterpart in CalculiX, whose NLGEOM elasticity is the Saint "
         "Venant-Kirchhoff law and whose finite-strain plasticity is multiplicative";
}

std::string calculix_plasticity_obstacle(const IsotropicMaterial& material) {
  const PlasticityParameters& p = material.plasticity();
  if (!p.enabled()) return "";
  if (p.criterion == YieldCriterion::Hill48) {
    return "material '" + material.name() +
           "' yields by Hill's 1948 anisotropic criterion, and CalculiX's *PLASTIC is von "
           "Mises (its anisotropic plasticity is a user material)";
  }
  for (int i = 0; i < p.num_backstresses; ++i) {
    if (p.backstresses[static_cast<std::size_t>(i)].recovery > 0.0) {
      return "material '" + material.name() +
             "' has Armstrong-Frederick (Chaboche) backstresses with dynamic recovery, which "
             "CalculiX's *PLASTIC does not model";
    }
  }
  if (p.kinematic_hardening_modulus > 0.0 || p.num_backstresses > 0) {
    return "material '" + material.name() +
           "' has linear kinematic hardening, and CalculiX's HARDENING=KINEMATIC does not "
           "reproduce Prager's linear kinematic hardening (a single element in uniaxial "
           "tension softens)";
  }
  return "";
}

std::vector<std::string> write_calculix_decks(const FemModel& model, const std::string& stem,
                                              const std::string& case_name,
                                              const CalculixNonlinearExport* nonlinear,
                                              const CalculixTransientExport* transient) {
  if (!model.finalized()) throw IoError("the model must be finalised before export");
  const std::vector<LoadCaseSpec>& specs = model.load_case_specs();
  bool faces_needed = false;
  for (const LoadCaseSpec& spec : specs) {
    faces_needed = faces_needed || !spec.pressures.empty() ||
                   spec.temperature.source == TemperatureSpec::Source::Conduction;
  }
  const std::vector<Mesh::BoundaryFace> boundary =
      faces_needed ? model.mesh().boundary_faces() : std::vector<Mesh::BoundaryFace>();

  std::vector<std::string> paths;
  const auto write = [&](const std::string& path, const auto& body) {
    std::ofstream out(path, std::ios::out | std::ios::trunc);
    if (!out) throw IoError("cannot open '" + path + "' for writing");
    body(out);
    out.flush();
    if (!out) throw IoError("failed while writing '" + path + "'");
    paths.push_back(path);
  };
  for (std::size_t l = 0; l < specs.size(); ++l) {
    const std::string base = stem + "_" + sanitise(specs[l].name);
    write(base + ".inp",
          [&](std::ostream& out) { write_static_deck(out, model, l, case_name, boundary); });
    if (model.load_case_data(l).conduction_solved) {
      write(base + "_conduction.inp", [&](std::ostream& out) {
        write_conduction_deck(out, model, l, case_name, boundary);
      });
    }
    if (nonlinear != nullptr && std::find(nonlinear->load_cases.begin(),
                                          nonlinear->load_cases.end(),
                                          l) != nonlinear->load_cases.end()) {
      for (const IsotropicMaterial& m : model.materials()) {
        const std::string obstacle = calculix_plasticity_obstacle(m);
        if (!obstacle.empty()) {
          throw IoError("the non-linear analysis of load case '" + specs[l].name +
                        "' cannot be exported to CalculiX: " + obstacle);
        }
      }
      write(base + (nonlinear->nlgeom ? "_nlgeom.inp" : "_small_strain.inp"),
            [&](std::ostream& out) {
              write_static_deck(out, model, l, case_name, boundary, nonlinear);
            });
    }
    if (transient != nullptr && std::find(transient->load_cases.begin(),
                                          transient->load_cases.end(),
                                          l) != transient->load_cases.end()) {
      const std::string obstacle = calculix_transient_obstacle(model, l, transient->options);
      if (!obstacle.empty()) {
        throw IoError("the transient of load case '" + specs[l].name +
                      "' cannot be exported to CalculiX: " + obstacle);
      }
      write(base + "_dynamic.inp", [&](std::ostream& out) {
        write_static_deck(out, model, l, case_name, boundary, nullptr, transient);
      });
    }
  }
  return paths;
}

}  // namespace sparlab
