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

/// The mechanical deck of load case `l`.
void write_static_deck(std::ostream& out, const FemModel& model, std::size_t l,
                       const std::string& case_name,
                       const std::vector<Mesh::BoundaryFace>& boundary) {
  const Mesh& mesh = model.mesh();
  const int dim = mesh.dim();
  // CalculiX numbers the DOFs of a node 1-3 (translations) and 4-6
  // (rotations), the same order as SparLab's.
  const int ndpn = model.dofs_per_node();
  const std::string element = calculix_element_type(model);
  const LoadCaseSpec& spec = model.load_case_specs()[l];
  const LoadCaseData& data = model.load_case_data(l);
  const bool thermal = data.temperature.size() > 0;
  const bool needs_density = spec.gravity.squaredNorm() > 0.0 || spec.centrifugal.enabled;

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
    if (needs_density) out << "*DENSITY\n" << field(mat.density()) << "\n";
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
  if (thermal) {
    out << "*INITIAL CONDITIONS, TYPE=TEMPERATURE\nNALL, "
        << field(common_reference_temperature(model)) << "\n";
  }

  out << "*STEP\n*STATIC\n";
  out << "*BOUNDARY\n";
  for (Index d : model.dofs().constrained_dofs()) {
    const Index node = d / ndpn;
    const int component = static_cast<int>(d % ndpn) + 1;
    out << node + 1 << ", " << component << ", " << component << ", "
        << field(model.dofs().prescribed_value(d)) << "\n";
  }

  // Point loads and tractions as the assembled nodal forces; every other load
  // in CalculiX's own form, so that CalculiX integrates it itself.
  Vector concentrated = data.mechanical;
  if (!spec.pressures.empty()) {
    LoadCaseSpec pressures_only;
    pressures_only.name = spec.name;
    pressures_only.pressures = spec.pressures;
    concentrated -= assemble_load_vector(mesh, model.element(), pressures_only,
                                         model.thickness(), model.integration());
  }
  bool any = false;
  for (Index n = 0; n < mesh.num_nodes(); ++n) {
    for (int k = 0; k < ndpn; ++k) {
      const Scalar f = concentrated(n * ndpn + k);
      if (f == 0.0) continue;
      if (!any) out << "*CLOAD\n";
      any = true;
      out << n + 1 << ", " << k + 1 << ", " << field(f) << "\n";
    }
  }
  const bool distributed = !spec.pressures.empty() || spec.has_body_loads();
  if (distributed) out << "*DLOAD\n";
  for (const PressureLoadSpec& p : spec.pressures) {
    for (const auto& face : calculix_faces_of(mesh, boundary, p.region)) {
      out << face.first + 1 << ", P" << face.second << ", " << field(p.pressure) << "\n";
    }
  }
  if (spec.gravity.squaredNorm() > 0.0) {
    const Scalar g = spec.gravity.norm();
    const Vector3 e = spec.gravity / g;
    out << "EALL, GRAV, " << field(g) << ", " << field(e.x()) << ", " << field(e.y()) << ", "
        << field(e.z()) << "\n";
  }
  for (std::size_t b = 0; b < spec.body_forces.size(); ++b) {
    const Vector3& density = spec.body_forces[b].force_density;
    static const char* const labels[3] = {"BX", "BY", "BZ"};
    for (int k = 0; k < dim; ++k) {
      if (density(k) != 0.0) {
        out << body_sets[b] << ", " << labels[k] << ", " << field(density(k)) << "\n";
      }
    }
  }
  if (spec.centrifugal.enabled) {
    const CentrifugalSpec& c = spec.centrifugal;
    const Vector3 axis = c.axis.normalized();
    out << "EALL, CENTRIF, " << field(c.angular_velocity * c.angular_velocity) << ", "
        << field(c.point.x()) << ", " << field(c.point.y()) << ", " << field(c.point.z())
        << ", " << field(axis.x()) << ", " << field(axis.y()) << ", " << field(axis.z())
        << "\n";
  }
  if (thermal) {
    out << "*TEMPERATURE\n";
    for (Index n = 0; n < mesh.num_nodes(); ++n) {
      out << n + 1 << ", " << field(data.temperature(n)) << "\n";
    }
  }
  out << "*NODE FILE\nU\n*EL FILE\nS\n*END STEP\n";
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

std::vector<std::string> write_calculix_decks(const FemModel& model, const std::string& stem,
                                              const std::string& case_name) {
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
  }
  return paths;
}

}  // namespace sparlab
