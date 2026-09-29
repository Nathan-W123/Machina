#include "sparlab/io/ResultWriter.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"
#include "sparlab/core/Version.hpp"
#include "sparlab/io/CsvWriter.hpp"
#include "sparlab/io/StlWriter.hpp"
#include "sparlab/io/VtkWriter.hpp"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace sparlab {
namespace {

std::string sanitise(const std::string& name) {
  std::string out;
  out.reserve(name.size());
  for (char c : name) {
    const bool ok = (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
                    (c >= '0' && c <= '9') || c == '-' || c == '_';
    out.push_back(ok ? c : '_');
  }
  if (out.empty()) out = "unnamed";
  return out;
}

/// [x, y] on a 2-D model, [x, y, z] on a 3-D one.
json::Value point_json(const Vector3& v, int dim) {
  json::Value out = json::Value::make_array();
  for (int i = 0; i < dim; ++i) out.push_back(json::Value::make_number(v(i)));
  return out;
}

const char* component_name(int k) { return dof_component_name(k); }

/// Vector magnitude over the model's components only, of a vector holding
/// `dim` translations per node.
Scalar magnitude(const Vector& full, Index n, int dim) {
  return dim == 2 ? std::hypot(full(n * dim + 0), full(n * dim + 1))
                  : std::hypot(full(n * dim + 0), full(n * dim + 1), full(n * dim + 2));
}

/// DOFs per node of a full-length nodal vector (a displacement, a mode, a
/// reaction): the translations of a continuum model, six for a shell or beam
/// model.
int nodal_dofs(const Mesh& mesh, const Vector& full) {
  const Index nn = mesh.num_nodes();
  if (nn == 0) return mesh.dim();
  const Index ndpn = static_cast<Index>(full.size()) / nn;
  if (ndpn * nn != full.size() || (ndpn != mesh.dim() && ndpn != kMaxDofsPerNode)) {
    std::ostringstream os;
    os << "a nodal vector of length " << full.size() << " does not fit the mesh's " << nn
       << " nodes with " << mesh.dim() << " or " << kMaxDofsPerNode << " DOFs per node";
    throw IoError(os.str());
  }
  return static_cast<int>(ndpn);
}

/// The translations of a full-length nodal vector, `dim` per node.
Vector translations(const Mesh& mesh, const Vector& full) {
  const int ndpn = nodal_dofs(mesh, full);
  const int dim = mesh.dim();
  if (ndpn == dim) return full;
  Vector out(static_cast<Eigen::Index>(mesh.num_nodes()) * dim);
  for (Index n = 0; n < mesh.num_nodes(); ++n) {
    for (int k = 0; k < dim; ++k) out(n * dim + k) = full(n * ndpn + k);
  }
  return out;
}

/// The rotations (three per node) of a shell or beam nodal vector; empty
/// for a continuum one.
Vector rotations(const Mesh& mesh, const Vector& full) {
  const int ndpn = nodal_dofs(mesh, full);
  if (ndpn != kMaxDofsPerNode) return Vector();
  Vector out(static_cast<Eigen::Index>(mesh.num_nodes()) * 3);
  for (Index n = 0; n < mesh.num_nodes(); ++n) {
    for (int k = 0; k < 3; ++k) out(n * 3 + k) = full(n * ndpn + 3 + k);
  }
  return out;
}

json::Value equilibrium_json(const EquilibriumCheck& eq, int dim) {
  json::Value out = json::Value::make_object();
  out.set("applied_force_N", point_json(eq.applied_force, dim));
  out.set("reaction_force_N", point_json(eq.reaction_force, dim));
  out.set("force_residual_N", point_json(eq.force_residual, dim));
  if (dim == 2) {
    // A plane model has a single moment component, about z.
    out.set("applied_moment_Nm", json::Value::make_number(eq.applied_moment.z()));
    out.set("reaction_moment_Nm", json::Value::make_number(eq.reaction_moment.z()));
    out.set("moment_residual_Nm", json::Value::make_number(eq.moment_residual.z()));
  } else {
    out.set("applied_moment_Nm", point_json(eq.applied_moment, 3));
    out.set("reaction_moment_Nm", point_json(eq.reaction_moment, 3));
    out.set("moment_residual_Nm", point_json(eq.moment_residual, 3));
  }
  out.set("relative_force_error", json::Value::make_number(eq.relative_force_error));
  out.set("relative_moment_error", json::Value::make_number(eq.relative_moment_error));
  return out;
}

json::Value modal_json(const ModalResult& modal) {
  json::Value out = json::Value::make_object();
  out.set("num_modes", json::Value::make_number(modal.eigenvalues.size()));
  out.set("eigenvalues_1_per_s2", json::array_of(modal.eigenvalues));
  out.set("angular_frequencies_rad_per_s", json::array_of(modal.angular_frequencies));
  out.set("frequencies_hz", json::array_of(modal.frequencies_hz));
  out.set("eigenpair_residuals", json::array_of(modal.modal_residuals));
  if (modal.modal_mass_fraction.size() > 0) {
    out.set("low_density_kinetic_energy_fraction",
            json::array_of(modal.modal_mass_fraction));
  }
  out.set("total_mass_kg", json::Value::make_number(modal.total_mass));
  out.set("subspace_size", json::Value::make_number(modal.subspace_size));
  out.set("iterations", json::Value::make_number(modal.iterations));
  out.set("converged", json::Value::make_bool(modal.converged));
  out.set("final_relative_change", json::Value::make_number(modal.final_change));
  if (!modal.linear_solver.empty()) {
    out.set("linear_solver", json::Value::make_string(modal.linear_solver));
    out.set("linear_iterations", json::Value::make_number(modal.linear_iterations));
  }
  out.set("warnings", json::array_of(modal.warnings));
  return out;
}

json::Value timings_json(const TimingLedger& timings) {
  json::Value out = json::Value::make_object();
  for (const auto& kv : timings.totals()) {
    out.set(kv.first, json::Value::make_number(kv.second));
  }
  return out;
}

json::Value mesh_quality_json(const MeshQuality& q) {
  json::Value out = json::Value::make_object();
  out.set("metric", json::Value::make_string(q.metric));
  out.set("min", json::Value::make_number(q.min));
  out.set("mean", json::Value::make_number(q.mean));
  out.set("worst_element", json::Value::make_number(q.worst_element));
  out.set("poor_threshold", json::Value::make_number(q.poor_threshold));
  out.set("poor_elements", json::Value::make_number(q.poor_elements));
  return out;
}

json::Value count_map_json(const std::map<std::string, Index>& counts) {
  json::Value out = json::Value::make_object();
  for (const auto& entry : counts) out.set(entry.first, json::Value::make_number(entry.second));
  return out;
}

/// What the mesh reader found in a mesh file, so a summary records how the
/// file was interpreted (sets, repairs, warnings) and not only its size.
json::Value mesh_file_json(const Configuration& config) {
  const MeshReadReport& r = config.mesh_report;
  json::Value out = json::Value::make_object();
  out.set("path", json::Value::make_string(config.mesh_file.path));
  out.set("format", json::Value::make_string(r.format));
  if (!r.version.empty()) out.set("version", json::Value::make_string(r.version));
  out.set("scale", json::Value::make_number(r.scale));
  out.set("nodes_in_file", json::Value::make_number(r.nodes_in_file));
  out.set("nodes_used", json::Value::make_number(r.nodes_used));
  out.set("unreferenced_nodes_dropped", json::Value::make_number(r.unreferenced_nodes));
  out.set("elements_in_file", json::Value::make_number(r.elements_in_file));
  out.set("cells", json::Value::make_number(r.cells));
  out.set("boundary_elements", json::Value::make_number(r.boundary_elements));
  out.set("cells_reoriented", json::Value::make_number(r.reoriented));
  out.set("duplicate_nodes", json::Value::make_number(r.duplicate_nodes));
  out.set("duplicates_merged", json::Value::make_bool(r.duplicates_merged));
  out.set("node_sets", count_map_json(r.node_sets));
  out.set("element_sets", count_map_json(r.element_sets));
  out.set("ignored", count_map_json(r.ignored));
  json::Value warnings = json::Value::make_array();
  for (const std::string& w : r.warnings) warnings.push_back(json::Value::make_string(w));
  out.set("warnings", warnings);
  return out;
}

json::Value mesh_stats_json(const Configuration& config, const FemModel& model) {
  const Mesh& mesh = model.mesh();
  const int dim = mesh.dim();
  json::Value out = json::Value::make_object();
  out.set("source", json::Value::make_string(to_string(config.mesh_kind)));
  out.set("element_type", json::Value::make_string(to_string(mesh.element_type())));
  out.set("dim", json::Value::make_number(dim));
  out.set("num_nodes", json::Value::make_number(mesh.num_nodes()));
  out.set("num_elements", json::Value::make_number(mesh.num_elements()));
  out.set("num_dofs", json::Value::make_number(model.dofs().num_dofs()));
  out.set("num_free_dofs", json::Value::make_number(model.dofs().num_free()));
  out.set("num_prescribed_dofs",
          json::Value::make_number(model.dofs().num_constrained()));
  const BoundingBox bb = mesh.bounding_box();
  json::Value box = json::Value::make_array();
  for (int i = 0; i < dim; ++i) box.push_back(json::Value::make_number(bb.lower(i)));
  for (int i = 0; i < dim; ++i) box.push_back(json::Value::make_number(bb.upper(i)));
  out.set(dim == 2 ? "bounding_box_xmin_ymin_xmax_ymax_m"
                   : "bounding_box_xmin_ymin_zmin_xmax_ymax_zmax_m",
          box);
  out.set("thickness_m", json::Value::make_number(model.thickness()));
  out.set("domain_volume_m3", json::Value::make_number(model.domain_volume()));
  if (mesh.structured_info().has_value()) {
    const StructuredGridInfo& info = *mesh.structured_info();
    json::Value grid = json::Value::make_object();
    grid.set("nx", json::Value::make_number(info.nx));
    grid.set("ny", json::Value::make_number(info.ny));
    if (dim == 3) grid.set("nz", json::Value::make_number(info.nz));
    grid.set("lx_m", json::Value::make_number(info.lx));
    grid.set("ly_m", json::Value::make_number(info.ly));
    if (dim == 3) grid.set("lz_m", json::Value::make_number(info.lz));
    grid.set("uniform", json::Value::make_bool(info.uniform));
    out.set("structured_grid", grid);
  } else if (config.mesh_kind == MeshKind::StructuredTri ||
             config.mesh_kind == MeshKind::StructuredTet) {
    // The simplex box meshes split the cells of this grid.
    json::Value grid = json::Value::make_object();
    grid.set("nx", json::Value::make_number(config.mesh_spec.nx));
    grid.set("ny", json::Value::make_number(config.mesh_spec.ny));
    if (dim == 3) grid.set("nz", json::Value::make_number(config.mesh_spec.nz));
    grid.set("lx_m", json::Value::make_number(config.mesh_spec.lx));
    grid.set("ly_m", json::Value::make_number(config.mesh_spec.ly));
    if (dim == 3) grid.set("lz_m", json::Value::make_number(config.mesh_spec.lz));
    grid.set("elements_per_cell", json::Value::make_number(dim == 2 ? 2 : 6));
    out.set("split_grid", grid);
  }
  out.set("mean_element_size_m", json::Value::make_number(mesh.mean_element_size()));
  out.set("quality", mesh_quality_json(mesh.quality()));
  if (config.mesh_kind == MeshKind::File) out.set("file", mesh_file_json(config));
  return out;
}

json::Value material_json(const IsotropicMaterial& m, StressState state) {
  json::Value out = json::Value::make_object();
  out.set("name", json::Value::make_string(m.name()));
  out.set("youngs_modulus_Pa", json::Value::make_number(m.youngs_modulus()));
  out.set("poisson_ratio", json::Value::make_number(m.poisson_ratio()));
  out.set("density_kg_per_m3", json::Value::make_number(m.density()));
  out.set("shear_modulus_Pa", json::Value::make_number(m.shear_modulus()));
  out.set("stress_state", json::Value::make_string(to_string(state)));
  if (m.thermal_expansion() != 0.0 || m.conductivity() != 0.0) {
    out.set("thermal_expansion_per_K", json::Value::make_number(m.thermal_expansion()));
    out.set("reference_temperature_K", json::Value::make_number(m.reference_temperature()));
    out.set("conductivity_W_per_mK", json::Value::make_number(m.conductivity()));
  }
  if (m.plasticity().enabled()) {
    const PlasticityParameters& p = m.plasticity();
    json::Value j = json::Value::make_object();
    j.set("yield_stress_Pa", json::Value::make_number(p.yield_stress));
    j.set("hardening_modulus_Pa", json::Value::make_number(p.hardening_modulus));
    j.set("kinematic_hardening_modulus_Pa",
          json::Value::make_number(p.kinematic_hardening_modulus));
    j.set("saturation_stress_Pa", json::Value::make_number(p.saturation_stress));
    j.set("saturation_rate", json::Value::make_number(p.saturation_rate));
    out.set("plasticity", j);
  }
  return out;
}

/// Every material of a multi-material model with the elements that use it.
json::Value materials_json(const FemModel& model) {
  json::Value out = json::Value::make_array();
  std::vector<Index> count(static_cast<std::size_t>(model.num_materials()), 0);
  for (Index e = 0; e < model.mesh().num_elements(); ++e) {
    ++count[static_cast<std::size_t>(model.element_material(e))];
  }
  for (int i = 0; i < model.num_materials(); ++i) {
    json::Value m = material_json(model.materials()[static_cast<std::size_t>(i)],
                                  model.stress_state());
    m.set("elements", json::Value::make_number(count[static_cast<std::size_t>(i)]));
    out.push_back(m);
  }
  return out;
}

/// Load parts of a load case beyond point loads and tractions: the body
/// loads' resultant and the temperature field (with its conduction solve).
json::Value load_parts_json(const FemModel& model, std::size_t l) {
  const LoadCaseData& data = model.load_case_data(l);
  const LoadCaseSpec& spec = model.load_case_specs()[l];
  const int dim = model.dim();
  json::Value out = json::Value::make_object();
  if (data.body.size() > 0) {
    out.set("body_load_resultant_N", point_json(data.body_resultant, dim));
    if (spec.gravity.squaredNorm() > 0.0) {
      out.set("gravity_m_per_s2", point_json(spec.gravity, dim));
    }
    if (spec.centrifugal.enabled) {
      out.set("angular_velocity_rad_per_s",
              json::Value::make_number(spec.centrifugal.angular_velocity));
    }
  }
  if (data.temperature.size() > 0) {
    json::Value t = json::Value::make_object();
    t.set("min_K", json::Value::make_number(data.temperature.minCoeff()));
    t.set("max_K", json::Value::make_number(data.temperature.maxCoeff()));
    t.set("thermal_self_energy_J", json::Value::make_number(data.thermal_self_energy));
    if (data.conduction_solved) {
      const ConductionSummary& c = data.conduction;
      json::Value cs = json::Value::make_object();
      cs.set("applied_heat_W", json::Value::make_number(c.applied_heat));
      cs.set("heat_through_prescribed_temperatures_W",
             json::Value::make_number(c.prescribed_heat));
      cs.set("gross_inflow_at_prescribed_temperatures_W",
             json::Value::make_number(c.prescribed_inflow));
      cs.set("gross_outflow_at_prescribed_temperatures_W",
             json::Value::make_number(c.prescribed_outflow));
      cs.set("relative_heat_balance_error", json::Value::make_number(c.relative_balance_error));
      cs.set("scaled_residual", json::Value::make_number(c.scaled_residual));
      cs.set("prescribed_nodes", json::Value::make_number(c.prescribed_nodes));
      cs.set("flux_faces", json::Value::make_number(c.flux_faces));
      cs.set("convection_faces", json::Value::make_number(c.convection_faces));
      cs.set("source_elements", json::Value::make_number(c.source_elements));
      cs.set("linear_solver", json::Value::make_string(c.solver));
      t.set("conduction", cs);
      t.set("source", json::Value::make_string("steady conduction"));
    } else {
      t.set("source", json::Value::make_string(
                          spec.temperature.source == TemperatureSpec::Source::Uniform
                              ? "uniform"
                              : "regions"));
    }
    out.set("temperature", t);
  }
  return out;
}

json::Value diagnostics_json(const ModelDiagnostics& diag) {
  json::Value out = json::Value::make_object();
  out.set("well_posed", json::Value::make_bool(diag.well_posed()));
  out.set("num_element_groups", json::Value::make_number(diag.components.size()));
  out.set("problems", json::array_of(diag.problems));
  json::Value groups = json::Value::make_array();
  for (const MeshComponent& c : diag.components) {
    json::Value g = json::Value::make_object();
    g.set("num_elements", json::Value::make_number(c.elements.size()));
    g.set("num_nodes", json::Value::make_number(c.nodes.size()));
    g.set("prescribed_dofs", json::Value::make_number(c.prescribed_dofs));
    g.set("rigid_body_null_dimension",
          json::Value::make_number(c.rigid_null_dimension));
    groups.push_back(g);
  }
  out.set("element_groups", groups);
  return out;
}

/// Coordinate column headers for a node or centroid table.
std::vector<std::string> coordinate_headers(int dim, const char* prefix) {
  std::vector<std::string> out;
  for (int k = 0; k < dim; ++k) {
    out.push_back(std::string(prefix) + component_name(k) + "[m]");
  }
  return out;
}

/// Per-component headers such as ux[m], uy[m](, uz[m]).
std::vector<std::string> component_headers(int dim, const char* stem, const char* unit) {
  std::vector<std::string> out;
  for (int k = 0; k < dim; ++k) {
    out.push_back(std::string(stem) + component_name(k) + "[" + unit + "]");
  }
  return out;
}

std::vector<std::string> concat(std::vector<std::string> a,
                                const std::vector<std::string>& b) {
  a.insert(a.end(), b.begin(), b.end());
  return a;
}

/// Digits of the file numbers of a series of `count` files (at least four).
int series_width(std::size_t count) {
  int width = 1;
  for (std::size_t c = count > 0 ? count - 1 : 0; c >= 10; c /= 10) ++width;
  return std::max(width, 4);
}

/// Largest nodal magnitude of every node of a vector of `dim` translations
/// per node.
Vector nodal_magnitudes(const Mesh& mesh, const Vector& translation) {
  Vector out(mesh.num_nodes());
  for (Index n = 0; n < mesh.num_nodes(); ++n) out(n) = magnitude(translation, n, mesh.dim());
  return out;
}

}  // namespace

ResultWriter::ResultWriter(std::string directory, const Configuration& config)
    : directory_(std::move(directory)), config_(config) {
  ensure_directory(directory_);
  log::debug("result directory: ", directory_);
}

std::string ResultWriter::file(const std::string& name) const {
  return path_join(directory_, name);
}

void ResultWriter::write_json(const std::string& file_name,
                              const json::Value& value) const {
  const std::string path = file(file_name);
  std::ofstream out(path, std::ios::out | std::ios::trunc);
  if (!out) throw IoError("cannot open '" + path + "' for writing");
  out << json::dump(value, 2) << '\n';
  out.flush();
  if (!out) throw IoError("failed while writing '" + path + "'");
}

void ResultWriter::write_config() const {
  write_json("config.json", config_.document);
}

void ResultWriter::write_mesh(const FemModel& model) const {
  const Mesh& mesh = model.mesh();
  const int dim = mesh.dim();
  json::Value doc = json::Value::make_object();
  doc.set("case", json::Value::make_string(config_.name));
  doc.set("element_type", json::Value::make_string(to_string(mesh.element_type())));
  doc.set("dim", json::Value::make_number(dim));
  doc.set("nodes_per_element", json::Value::make_number(mesh.nodes_per_elem()));

  json::Value nodes = json::Value::make_array();
  for (Index n = 0; n < mesh.num_nodes(); ++n) {
    nodes.push_back(point_json(mesh.node(n), dim));
  }
  doc.set("nodes_m", nodes);

  json::Value elements = json::Value::make_array();
  for (Index e = 0; e < mesh.num_elements(); ++e) {
    json::Value row = json::Value::make_array();
    const Index* en = mesh.element_nodes(e);
    for (int a = 0; a < mesh.nodes_per_elem(); ++a) {
      row.push_back(json::Value::make_number(en[a]));
    }
    elements.push_back(row);
  }
  doc.set("elements", elements);
  if (!model.single_material()) {
    // Index into the summary's "materials" of every element.
    json::Value materials = json::Value::make_array();
    for (Index e = 0; e < mesh.num_elements(); ++e) {
      materials.push_back(json::Value::make_number(model.element_material(e)));
    }
    doc.set("element_materials", materials);
  }

  // Prescribed DOFs, for the boundary-condition figure.
  const int ndpn = model.dofs_per_node();
  doc.set("dofs_per_node", json::Value::make_number(ndpn));
  json::Value constraints = json::Value::make_array();
  for (Index d : model.dofs().constrained_dofs()) {
    json::Value entry = json::Value::make_object();
    const int k = static_cast<int>(d % ndpn);
    entry.set("node", json::Value::make_number(d / ndpn));
    entry.set("component", json::Value::make_string(component_name(k)));
    entry.set(k < 3 ? "value_m" : "value_rad",
              json::Value::make_number(model.dofs().prescribed_value(d)));
    constraints.push_back(entry);
  }
  doc.set("prescribed_dofs", constraints);

  // Applied nodal forces per load case (non-zero entries only).
  json::Value cases = json::Value::make_array();
  const std::vector<Vector>& loads = model.load_vectors();
  for (std::size_t l = 0; l < loads.size(); ++l) {
    json::Value entry = json::Value::make_object();
    entry.set("name", json::Value::make_string(model.load_case_specs()[l].name));
    entry.set("weight",
              json::Value::make_number(model.load_case_specs()[l].weight));
    json::Value forces = json::Value::make_array();
    for (Index n = 0; n < mesh.num_nodes(); ++n) {
      bool nonzero = false;
      for (int k = 0; k < ndpn; ++k) nonzero = nonzero || loads[l](n * ndpn + k) != 0.0;
      if (!nonzero) continue;
      json::Value f = json::Value::make_object();
      f.set("node", json::Value::make_number(n));
      for (int k = 0; k < dim; ++k) {
        f.set(std::string("f") + component_name(k) + "_N",
              json::Value::make_number(loads[l](n * ndpn + k)));
      }
      if (ndpn == kMaxDofsPerNode) {
        for (int k = 0; k < 3; ++k) {
          f.set(std::string("m") + component_name(k) + "_Nm",
                json::Value::make_number(loads[l](n * ndpn + 3 + k)));
        }
      }
      forces.push_back(f);
    }
    entry.set("nodal_forces", forces);
    cases.push_back(entry);
  }
  doc.set("load_cases", cases);

  write_json("mesh.json", doc);
}

void ResultWriter::write_displacement(const Mesh& mesh, const std::string& load_case,
                                      const Vector& full_displacement,
                                      const std::string& stem) const {
  const int dim = mesh.dim();
  const Vector displacement = translations(mesh, full_displacement);
  const Vector rotation = rotations(mesh, full_displacement);
  std::vector<std::string> header{"node"};
  header = concat(header, coordinate_headers(dim, ""));
  header = concat(header, component_headers(dim, "u", "m"));
  header.push_back("umag[m]");
  if (rotation.size() > 0) header = concat(header, component_headers(3, "r", "rad"));
  CsvWriter csv(file(stem + "_" + sanitise(load_case) + ".csv"), header);
  for (Index n = 0; n < mesh.num_nodes(); ++n) {
    const Vector3 x = mesh.node(n);
    std::vector<Scalar> row;
    for (int k = 0; k < dim; ++k) row.push_back(x(k));
    for (int k = 0; k < dim; ++k) row.push_back(displacement(n * dim + k));
    row.push_back(magnitude(displacement, n, dim));
    if (rotation.size() > 0) {
      for (int k = 0; k < 3; ++k) row.push_back(rotation(n * 3 + k));
    }
    csv.row(n, row);
  }
  csv.close();
}

void ResultWriter::write_stress(const Mesh& mesh, const std::string& load_case,
                                const StressField& field, const Vector* density) const {
  const int dim = mesh.dim();
  std::vector<std::string> header;
  if (dim == 2) {
    header = {"element", "cx[m]", "cy[m]", "area[m2]", "density[-]", "exx[-]",
              "eyy[-]", "gxy[-]", "sxx[Pa]", "syy[Pa]", "sxy[Pa]", "von_mises[Pa]",
              "solid_von_mises[Pa]", "principal_max[Pa]", "principal_min[Pa]",
              "strain_energy[J]"};
  } else {
    header = {"element", "cx[m]", "cy[m]", "cz[m]", "volume[m3]", "density[-]",
              "exx[-]", "eyy[-]", "ezz[-]", "gxy[-]", "gyz[-]", "gzx[-]",
              "sxx[Pa]", "syy[Pa]", "szz[Pa]", "sxy[Pa]", "syz[Pa]", "szx[Pa]",
              "von_mises[Pa]", "solid_von_mises[Pa]", "principal_max[Pa]",
              "principal_mid[Pa]", "principal_min[Pa]", "strain_energy[J]"};
  }
  CsvWriter csv(file("stress_" + sanitise(load_case) + ".csv"), header);
  const Eigen::Index nv = field.element_strain.rows();
  for (Index e = 0; e < mesh.num_elements(); ++e) {
    const Vector3 c = mesh.element_centroid(e);
    std::vector<Scalar> row;
    for (int k = 0; k < dim; ++k) row.push_back(c(k));
    row.push_back(mesh.element_measure(e));
    row.push_back(density ? (*density)(e) : 1.0);
    for (Eigen::Index i = 0; i < nv; ++i) row.push_back(field.element_strain(i, e));
    for (Eigen::Index i = 0; i < nv; ++i) row.push_back(field.element_stress(i, e));
    row.push_back(field.element_von_mises(e));
    row.push_back(field.element_solid_von_mises(e));
    row.push_back(field.element_principal_max(e));
    if (dim == 3) row.push_back(field.element_principal_mid(e));
    row.push_back(field.element_principal_min(e));
    row.push_back(field.element_strain_energy(e));
    csv.row(e, row);
  }
  csv.close();
}

void ResultWriter::write_reactions(const Mesh& mesh, const DofManager& dofs,
                                   const std::string& load_case,
                                   const Vector& full_reactions,
                                   const std::string& stem) const {
  const int dim = mesh.dim();
  const int ndpn = dofs.dofs_per_node();
  const Vector reactions = translations(mesh, full_reactions);
  const Vector moments = rotations(mesh, full_reactions);
  std::vector<std::string> header{"node"};
  header = concat(header, coordinate_headers(dim, ""));
  header = concat(header, component_headers(dim, "r", "N"));
  header.push_back("rmag[N]");
  if (moments.size() > 0) header = concat(header, component_headers(3, "m", "Nm"));
  CsvWriter csv(file(stem + "_" + sanitise(load_case) + ".csv"), header);
  for (Index n = 0; n < mesh.num_nodes(); ++n) {
    bool constrained = false;
    for (int k = 0; k < ndpn; ++k) {
      constrained = constrained || dofs.is_constrained(n * ndpn + k);
    }
    if (!constrained) continue;
    const Vector3 x = mesh.node(n);
    std::vector<Scalar> row;
    for (int k = 0; k < dim; ++k) row.push_back(x(k));
    for (int k = 0; k < dim; ++k) row.push_back(reactions(n * dim + k));
    row.push_back(magnitude(reactions, n, dim));
    if (moments.size() > 0) {
      for (int k = 0; k < 3; ++k) row.push_back(moments(n * 3 + k));
    }
    csv.row(n, row);
  }
  csv.close();
}

void ResultWriter::write_temperature(const Mesh& mesh, const std::string& load_case,
                                     const Vector& temperature) const {
  const int dim = mesh.dim();
  std::vector<std::string> header{"node"};
  header = concat(header, coordinate_headers(dim, ""));
  header.push_back("temperature[K]");
  CsvWriter csv(file("temperature_" + sanitise(load_case) + ".csv"), header);
  for (Index n = 0; n < mesh.num_nodes(); ++n) {
    const Vector3 x = mesh.node(n);
    std::vector<Scalar> row;
    for (int k = 0; k < dim; ++k) row.push_back(x(k));
    row.push_back(temperature(n));
    csv.row(n, row);
  }
  csv.close();
}

void ResultWriter::write_static_vtk(const Mesh& mesh, const std::string& load_case,
                                    const Vector& full_displacement,
                                    const StressField& field, const Vector* density,
                                    const Vector* stiffness_factor,
                                    const Vector* temperature) const {
  const int dim = mesh.dim();
  const Vector displacement = translations(mesh, full_displacement);
  const Vector rotation = rotations(mesh, full_displacement);
  VtkWriter writer(mesh, "SparLab static solution: " + config_.name + " / " + load_case);
  Vector mag(mesh.num_nodes());
  for (Index n = 0; n < mesh.num_nodes(); ++n) mag(n) = magnitude(displacement, n, dim);
  writer.add_point_vectors("displacement", displacement);
  if (rotation.size() > 0) writer.add_point_vectors("rotation", rotation);
  writer.add_point_scalars("displacement_magnitude", mag);
  writer.add_point_scalars("nodal_von_mises", field.nodal_von_mises);
  if (temperature != nullptr) writer.add_point_scalars("temperature", *temperature);
  if (field.element_sigma_zz.size() > 0) {
    writer.add_cell_scalars("sigma_zz", field.element_sigma_zz);
  }
  if (dim == 2) {
    writer.add_cell_scalars("sigma_xx", field.element_stress.row(0).transpose());
    writer.add_cell_scalars("sigma_yy", field.element_stress.row(1).transpose());
    writer.add_cell_scalars("sigma_xy", field.element_stress.row(2).transpose());
  } else {
    writer.add_cell_scalars("sigma_xx", field.element_stress.row(0).transpose());
    writer.add_cell_scalars("sigma_yy", field.element_stress.row(1).transpose());
    writer.add_cell_scalars("sigma_zz", field.element_stress.row(2).transpose());
    writer.add_cell_scalars("sigma_xy", field.element_stress.row(3).transpose());
    writer.add_cell_scalars("sigma_yz", field.element_stress.row(4).transpose());
    writer.add_cell_scalars("sigma_zx", field.element_stress.row(5).transpose());
  }
  writer.add_cell_scalars("von_mises", field.element_von_mises);
  writer.add_cell_scalars("principal_max", field.element_principal_max);
  writer.add_cell_scalars("principal_min", field.element_principal_min);
  writer.add_cell_scalars("strain_energy", field.element_strain_energy);
  if (density) writer.add_cell_scalars("density", *density);
  if (stiffness_factor) writer.add_cell_scalars("stiffness_factor", *stiffness_factor);
  writer.write(file("fields_" + sanitise(load_case) + ".vtk"));
}

void ResultWriter::write_modal(const Mesh& mesh, const ModalResult& modal,
                               const std::string& tag) const {
  const int dim = mesh.dim();
  const std::string suffix = tag.empty() ? "" : "_" + sanitise(tag);
  {
    CsvWriter csv(file("modes" + suffix + ".csv"),
                  {"mode", "eigenvalue[1/s2]", "omega[rad/s]", "frequency[Hz]",
                   "eigenpair_residual[-]", "low_density_ke_fraction[-]"});
    for (Eigen::Index i = 0; i < modal.eigenvalues.size(); ++i) {
      const Scalar frac =
          modal.modal_mass_fraction.size() > i ? modal.modal_mass_fraction(i) : 0.0;
      csv.row(static_cast<Index>(i),
              {modal.eigenvalues(i), modal.angular_frequencies(i),
               modal.frequencies_hz(i), modal.modal_residuals(i), frac});
    }
    csv.close();
  }

  if (!config_.output.write_mode_shapes) return;

  // Mode shapes carry the nodal rotations too on a shell or beam model; the
  // CSV lists them after the translations of each mode.
  const int ndpn = modal.mode_shapes.cols() > 0
                       ? nodal_dofs(mesh, Vector(modal.mode_shapes.col(0)))
                       : dim;
  std::vector<std::string> header{"node"};
  header = concat(header, coordinate_headers(dim, ""));
  for (Eigen::Index i = 0; i < modal.mode_shapes.cols(); ++i) {
    for (int k = 0; k < ndpn; ++k) {
      std::ostringstream name;
      if (k < 3) {
        name << "u" << component_name(k) << "_mode" << i << "[m]";
      } else {
        name << component_name(k) << "_mode" << i << "[rad]";
      }
      header.push_back(name.str());
    }
  }
  CsvWriter csv(file("mode_shapes" + suffix + ".csv"), header);
  for (Index n = 0; n < mesh.num_nodes(); ++n) {
    const Vector3 x = mesh.node(n);
    std::vector<Scalar> row;
    for (int k = 0; k < dim; ++k) row.push_back(x(k));
    for (Eigen::Index i = 0; i < modal.mode_shapes.cols(); ++i) {
      for (int k = 0; k < ndpn; ++k) row.push_back(modal.mode_shapes(n * ndpn + k, i));
    }
    csv.row(n, row);
  }
  csv.close();

  if (config_.output.write_vtk) {
    for (Eigen::Index i = 0; i < modal.mode_shapes.cols(); ++i) {
      std::ostringstream title;
      title << "SparLab mode " << i << " at " << modal.frequencies_hz(i) << " Hz";
      VtkWriter writer(mesh, title.str());
      const Vector shape = translations(mesh, modal.mode_shapes.col(i));
      writer.add_point_vectors("mode_shape", shape);
      Vector mag(mesh.num_nodes());
      for (Index n = 0; n < mesh.num_nodes(); ++n) mag(n) = magnitude(shape, n, dim);
      writer.add_point_scalars("mode_shape_magnitude", mag);
      std::ostringstream name;
      name << "mode" << suffix << "_" << i << ".vtk";
      writer.write(file(name.str()));
    }
  }
}

void ResultWriter::write_buckling(const Mesh& mesh, const std::vector<BucklingResult>& results,
                                  const std::string& tag) const {
  const int dim = mesh.dim();
  const std::string suffix = tag.empty() ? "" : "_" + sanitise(tag);
  {
    CsvWriter csv(file("buckling" + suffix + ".csv"),
                  {"load_case", "mode", "load_factor[-]", "eigenpair_residual[-]",
                   "solid_energy_fraction[-]"});
    for (const BucklingResult& r : results) {
      for (Eigen::Index i = 0; i < r.load_factors.size(); ++i) {
        std::ostringstream lf;
        lf << std::setprecision(17) << r.load_factors(i);
        std::ostringstream res;
        res << std::setprecision(6) << r.residuals(i);
        std::ostringstream solid;
        solid << std::setprecision(6) << r.solid_energy_fraction(i);
        csv.raw_row({r.load_case, std::to_string(i + 1), lf.str(), res.str(), solid.str()});
      }
    }
    csv.close();
  }
  if (!config_.output.write_mode_shapes || !config_.output.write_vtk) return;
  for (const BucklingResult& r : results) {
    for (Eigen::Index i = 0; i < r.mode_shapes.cols(); ++i) {
      Vector shape = translations(mesh, r.mode_shapes.col(i));
      Scalar largest = 0.0;
      for (Index n = 0; n < mesh.num_nodes(); ++n) largest = std::max(largest, magnitude(shape, n, dim));
      if (largest > 0.0) shape /= largest;
      std::ostringstream title;
      title << "SparLab buckling mode " << i + 1 << " of load case " << r.load_case
            << " at load factor " << r.load_factors(i) << " (max |phi| = 1)";
      VtkWriter writer(mesh, title.str());
      writer.add_point_vectors("buckling_mode", shape);
      Vector mag(mesh.num_nodes());
      for (Index n = 0; n < mesh.num_nodes(); ++n) mag(n) = magnitude(shape, n, dim);
      writer.add_point_scalars("buckling_mode_magnitude", mag);
      std::ostringstream name;
      name << "buckling" << suffix << "_" << sanitise(r.load_case) << "_" << i + 1 << ".vtk";
      writer.write(file(name.str()));
    }
  }
}

void ResultWriter::write_nonlinear(const FemModel& model, const NonlinearResult& result) const {
  const Mesh& mesh = model.mesh();
  const int dim = mesh.dim();
  const std::string lc = sanitise(result.load_case_name);
  {
    std::vector<std::string> header{"step",          "load_factor[-]",     "iterations[-]",
                                    "cuts[-]",       "residual[-]",        "arc_length[m]",
                                    "negative_pivots[-]", "max_displacement[m]"};
    if (result.plastic) {
      header.push_back("yielding_points[-]");
      header.push_back("max_plastic_strain[-]");
    }
    for (std::size_t i = 0; i < result.monitor_names.size(); ++i) {
      header.push_back(result.monitor_names[i] + "[" + result.monitor_units[i] + "]");
    }
    CsvWriter csv(file("nonlinear_" + lc + ".csv"), header);
    for (const NonlinearStep& s : result.steps) {
      std::vector<Scalar> row{s.load_factor, static_cast<Scalar>(s.iterations),
                              static_cast<Scalar>(s.cuts), s.residual, s.arc_length,
                              static_cast<Scalar>(s.negative_pivots), s.max_displacement};
      if (result.plastic) {
        row.push_back(static_cast<Scalar>(s.yielding_points));
        row.push_back(s.max_plastic_strain);
      }
      row.insert(row.end(), s.monitors.begin(), s.monitors.end());
      csv.row(s.index, row);
    }
    csv.close();
  }
  const bool plane_strain = dim == 2 && model.stress_state() == StressState::PlaneStrain;
  const Index ne = mesh.num_elements();
  if (config_.output.write_csv) {
    write_displacement(mesh, result.load_case_name, result.displacement,
                       "nonlinear_displacement");
    write_reactions(mesh, model.dofs(), result.load_case_name, result.reactions,
                    "nonlinear_reactions");
    // Cauchy (true) stress on the deformed element and the second
    // Piola-Kirchhoff stress, both averaged over the integration points; the
    // coordinates are the reference centroid.
    std::vector<std::string> header{"element"};
    header = concat(header, coordinate_headers(dim, "c"));
    header.push_back(dim == 2 ? "area[m2]" : "volume[m3]");
    const std::vector<std::string> pairs =
        dim == 2 ? std::vector<std::string>{"xx", "yy", "xy"}
                 : std::vector<std::string>{"xx", "yy", "zz", "xy", "yz", "zx"};
    for (const std::string& p : pairs) header.push_back("s" + p + "[Pa]");
    if (plane_strain) header.push_back("szz[Pa]");
    header.push_back("von_mises[Pa]");
    for (const std::string& p : pairs) header.push_back("pk2_" + p + "[Pa]");
    if (result.plastic) header.push_back("equivalent_plastic_strain[-]");
    CsvWriter csv(file("nonlinear_stress_" + lc + ".csv"), header);
    for (Index e = 0; e < ne; ++e) {
      const Vector3 c = mesh.element_centroid(e);
      std::vector<Scalar> row;
      for (int k = 0; k < dim; ++k) row.push_back(c(k));
      row.push_back(mesh.element_measure(e));
      for (Eigen::Index i = 0; i < result.element_cauchy.rows(); ++i) {
        row.push_back(result.element_cauchy(i, e));
      }
      if (plane_strain) row.push_back(result.element_cauchy_zz(e));
      row.push_back(result.element_von_mises(e));
      for (Eigen::Index i = 0; i < result.element_piola_kirchhoff.rows(); ++i) {
        row.push_back(result.element_piola_kirchhoff(i, e));
      }
      if (result.plastic) row.push_back(result.element_plastic_strain(e));
      csv.row(e, row);
    }
    csv.close();
  }
  if (config_.output.write_vtk) {
    std::ostringstream title;
    title << "SparLab non-linear solution: " << config_.name << " / " << result.load_case_name
          << " at load factor " << result.load_factor;
    VtkWriter writer(mesh, title.str());
    Vector mag(mesh.num_nodes());
    for (Index n = 0; n < mesh.num_nodes(); ++n) mag(n) = magnitude(result.displacement, n, dim);
    writer.add_point_vectors("displacement", result.displacement);
    writer.add_point_scalars("displacement_magnitude", mag);
    const std::vector<std::string> pairs =
        dim == 2 ? std::vector<std::string>{"xx", "yy", "xy"}
                 : std::vector<std::string>{"xx", "yy", "zz", "xy", "yz", "zx"};
    for (std::size_t i = 0; i < pairs.size(); ++i) {
      writer.add_cell_scalars("cauchy_" + pairs[i],
                              result.element_cauchy.row(static_cast<Eigen::Index>(i)).transpose());
    }
    if (plane_strain) writer.add_cell_scalars("cauchy_zz", result.element_cauchy_zz);
    writer.add_cell_scalars("von_mises", result.element_von_mises);
    if (result.plastic) {
      writer.add_cell_scalars("equivalent_plastic_strain", result.element_plastic_strain);
    }
    writer.write(file("nonlinear_" + lc + ".vtk"));
  }
}

void ResultWriter::write_transient(const FemModel& model, const TransientResult& result) const {
  const Mesh& mesh = model.mesh();
  const int dim = mesh.dim();
  const std::string lc = sanitise(result.load_case_name);
  {
    // The energy columns: kinetic T = v'Mv/2, the strain energy U (with
    // plasticity the stored energy, elastic plus hardening), the damping
    // dissipation D and the external work W since t = 0; the balance
    // E_0 + W - T - U - D is the energy the step sequence lost: round-off
    // for the trapezoidal rule on a linear model, the numerical dissipation
    // for alpha < 0, and on top the plastic dissipation of a plastic model.
    std::vector<std::string> header{"step",
                                    "time[s]",
                                    "max_displacement[m]",
                                    "kinetic_energy[J]",
                                    result.plastic ? "stored_energy[J]" : "strain_energy[J]",
                                    "damping_energy[J]",
                                    "external_work[J]",
                                    "energy_balance[J]"};
    if (result.nonlinear) {
      header.push_back("iterations[-]");
      header.push_back("cuts[-]");
    }
    if (result.plastic) {
      header.push_back("yielding_points[-]");
      header.push_back("max_plastic_strain[-]");
    }
    for (std::size_t i = 0; i < result.monitor_names.size(); ++i) {
      header.push_back(result.monitor_names[i] + "[" + result.monitor_units[i] + "]");
    }
    CsvWriter csv(file("transient_" + lc + ".csv"), header);
    const Scalar e0 = result.steps.empty()
                          ? 0.0
                          : result.steps.front().kinetic_energy + result.steps.front().strain_energy;
    for (const TransientStep& s : result.steps) {
      std::vector<Scalar> row{s.time,
                              s.max_displacement,
                              s.kinetic_energy,
                              s.strain_energy,
                              s.damping_energy,
                              s.external_work,
                              e0 + s.external_work - s.kinetic_energy - s.strain_energy -
                                  s.damping_energy};
      if (result.nonlinear) {
        row.push_back(static_cast<Scalar>(s.iterations));
        row.push_back(static_cast<Scalar>(s.cuts));
      }
      if (result.plastic) {
        row.push_back(static_cast<Scalar>(s.yielding_points));
        row.push_back(s.max_plastic_strain);
      }
      row.insert(row.end(), s.monitors.begin(), s.monitors.end());
      csv.row(s.index, row);
    }
    csv.close();
  }
  if (config_.output.write_csv) {
    // The final state: displacement, velocity and acceleration per node.
    const Vector u = translations(mesh, result.displacement);
    const Vector v = translations(mesh, result.velocity);
    const Vector a = translations(mesh, result.acceleration);
    std::vector<std::string> header{"node"};
    header = concat(header, coordinate_headers(dim, ""));
    header = concat(header, component_headers(dim, "u", "m"));
    header = concat(header, component_headers(dim, "v", "m/s"));
    header = concat(header, component_headers(dim, "a", "m/s^2"));
    CsvWriter csv(file("transient_state_" + lc + ".csv"), header);
    for (Index n = 0; n < mesh.num_nodes(); ++n) {
      const Vector3 x = mesh.node(n);
      std::vector<Scalar> row;
      for (int k = 0; k < dim; ++k) row.push_back(x(k));
      for (const Vector* field : {&u, &v, &a}) {
        for (int k = 0; k < dim; ++k) row.push_back((*field)(n * dim + k));
      }
      csv.row(n, row);
    }
    csv.close();
    write_reactions(mesh, model.dofs(), result.load_case_name, result.reactions,
                    "transient_reactions");
  }
  if (config_.output.write_vtk && !result.snapshots.empty()) {
    // A numbered series with ParaView's file-series index, which carries the
    // time of every file.
    const int width = series_width(result.snapshots.size());
    json::Value files = json::Value::make_array();
    for (std::size_t k = 0; k < result.snapshots.size(); ++k) {
      const TransientSnapshot& snap = result.snapshots[k];
      std::ostringstream title;
      title << "SparLab transient: " << config_.name << " / " << result.load_case_name
            << " at t = " << snap.time << " s";
      VtkWriter writer(mesh, title.str());
      const Vector u = translations(mesh, snap.displacement);
      const Vector v = translations(mesh, snap.velocity);
      writer.add_point_vectors("displacement", u);
      writer.add_point_scalars("displacement_magnitude", nodal_magnitudes(mesh, u));
      writer.add_point_vectors("velocity", v);
      writer.add_point_scalars("velocity_magnitude", nodal_magnitudes(mesh, v));
      std::ostringstream name;
      name << "transient_" << lc << "_" << std::setw(width) << std::setfill('0') << k << ".vtk";
      writer.write(file(name.str()));
      json::Value entry = json::Value::make_object();
      entry.set("name", json::Value::make_string(name.str()));
      entry.set("time", json::Value::make_number(snap.time));
      files.push_back(entry);
    }
    json::Value series = json::Value::make_object();
    series.set("file-series-version", json::Value::make_string("1.0"));
    series.set("files", files);
    write_json("transient_" + lc + ".vtk.series", series);
  }
}

void ResultWriter::write_frequency_response(const FemModel& model,
                                            const FrequencyResponseResult& result) const {
  constexpr Scalar kDegrees = 180.0 / 3.14159265358979323846;
  const Mesh& mesh = model.mesh();
  const std::string lc = sanitise(result.load_case_name);
  {
    // Every monitor as its complex amplitude Z (the response is
    // Re(Z e^{i omega t}) under the load f cos(omega t)): real and imaginary
    // parts, modulus and phase arg Z, the lag behind the load being -arg Z.
    std::vector<std::string> header{"point", "frequency[Hz]", "max_displacement[m]"};
    for (std::size_t i = 0; i < result.monitor_names.size(); ++i) {
      const std::string& name = result.monitor_names[i];
      const std::string& unit = result.monitor_units[i];
      header.push_back(name + "_re[" + unit + "]");
      header.push_back(name + "_im[" + unit + "]");
      header.push_back(name + "_abs[" + unit + "]");
      header.push_back(name + "_phase[deg]");
    }
    CsvWriter csv(file("frequency_response_" + lc + ".csv"), header);
    for (std::size_t j = 0; j < result.points.size(); ++j) {
      const FrequencyPoint& p = result.points[j];
      std::vector<Scalar> row{p.frequency, p.max_displacement};
      for (const ComplexScalar& z : p.monitors) {
        row.push_back(z.real());
        row.push_back(z.imag());
        row.push_back(std::abs(z));
        row.push_back(std::arg(z) * kDegrees);
      }
      csv.row(static_cast<Index>(j), row);
    }
    csv.close();
  }
  if (config_.output.write_csv) {
    // The complex nodal amplitudes at every snapshot frequency, one file
    // each (numbered as the VTK series), with the frequency in the name of
    // the first column's header.
    const int dim = mesh.dim();
    const int width = series_width(result.snapshots.size());
    for (std::size_t k = 0; k < result.snapshots.size(); ++k) {
      const FrequencySnapshot& snap = result.snapshots[k];
      const Vector re = translations(mesh, Vector(snap.displacement.real()));
      const Vector im = translations(mesh, Vector(snap.displacement.imag()));
      const Vector peak = harmonic_peak_displacements(model, snap.displacement);
      std::ostringstream first;
      first << "node@" << std::setprecision(17) << snap.frequency << "Hz";
      std::vector<std::string> header{first.str()};
      header = concat(header, coordinate_headers(dim, ""));
      for (int c = 0; c < dim; ++c) {
        header.push_back(std::string("u") + component_name(c) + "_re[m]");
        header.push_back(std::string("u") + component_name(c) + "_im[m]");
      }
      header.push_back("peak[m]");
      std::ostringstream name;
      name << "frequency_response_field_" << lc << "_" << std::setw(width) << std::setfill('0')
           << k << ".csv";
      CsvWriter field(file(name.str()), header);
      for (Index n = 0; n < mesh.num_nodes(); ++n) {
        const Vector3 x = mesh.node(n);
        std::vector<Scalar> row;
        for (int c = 0; c < dim; ++c) row.push_back(x(c));
        for (int c = 0; c < dim; ++c) {
          row.push_back(re(n * dim + c));
          row.push_back(im(n * dim + c));
        }
        row.push_back(peak(n));
        field.row(n, row);
      }
      field.close();
    }
  }
  if (config_.output.write_vtk && !result.snapshots.empty()) {
    // One file per snapshot frequency, with a file-series index whose "time"
    // is the frequency in Hz (so ParaView steps through the frequencies).
    const int width = series_width(result.snapshots.size());
    json::Value files = json::Value::make_array();
    for (std::size_t k = 0; k < result.snapshots.size(); ++k) {
      const FrequencySnapshot& snap = result.snapshots[k];
      std::ostringstream title;
      title << "SparLab frequency response: " << config_.name << " / " << result.load_case_name
            << " at " << snap.frequency << " Hz";
      VtkWriter writer(mesh, title.str());
      const Vector re = translations(mesh, Vector(snap.displacement.real()));
      const Vector im = translations(mesh, Vector(snap.displacement.imag()));
      writer.add_point_vectors("displacement_real", re);
      writer.add_point_vectors("displacement_imag", im);
      // The largest displacement of each node over a cycle.
      writer.add_point_scalars("displacement_peak",
                               harmonic_peak_displacements(model, snap.displacement));
      std::ostringstream name;
      name << "frequency_response_" << lc << "_" << std::setw(width) << std::setfill('0') << k
           << ".vtk";
      writer.write(file(name.str()));
      json::Value entry = json::Value::make_object();
      entry.set("name", json::Value::make_string(name.str()));
      entry.set("time", json::Value::make_number(snap.frequency));
      files.push_back(entry);
    }
    json::Value series = json::Value::make_object();
    series.set("file-series-version", json::Value::make_string("1.0"));
    series.set("files", files);
    write_json("frequency_response_" + lc + ".vtk.series", series);
  }
}

void ResultWriter::write_history(const TopologyOptimizationResult& result) const {
  const bool mma = result.method == OptimizerMethod::MMA;
  std::vector<std::string> header{
      "iteration", "penalty[-]", "compliance[J]", "volume[m3]", "volume_fraction[-]",
      "max_design_change[-]", "lagrange_multiplier[-]", "grey_level[-]",
      mma ? "mma_iterations[-]" : "oc_bisections[-]", "volume_converged[-]",
      "seconds[s]"};
  if (mma) {
    header.insert(header.end(), {"max_stress_ratio[-]", "stress_constraint[-]",
                                 "constraint_violation[-]"});
  }
  header.push_back("linear_iterations[-]");
  if (result.projected) header.push_back("beta[-]");
  if (result.buckling_constrained) {
    header.insert(header.end(), {"min_load_factor[-]", "buckling_constraint[-]",
                                 "buckling_subspace_iterations[-]"});
  }
  if (result.robust) {
    header.insert(header.end(), {"eroded_volume_fraction[-]", "dilated_volume_fraction[-]",
                                 "dilated_target_fraction[-]"});
  }
  CsvWriter csv(file("history.csv"), header);
  for (const TopologyIteration& it : result.history) {
    std::vector<Scalar> row{it.penalty, it.compliance, it.volume, it.volume_fraction,
                            it.max_change, it.lambda, it.gray_level,
                            static_cast<Scalar>(it.bisections),
                            it.volume_converged ? 1.0 : 0.0, it.seconds};
    if (mma) {
      row.insert(row.end(), {it.max_stress_ratio, it.stress_constraint,
                             it.constraint_violation});
    }
    row.push_back(static_cast<Scalar>(it.linear_iterations));
    if (result.projected) row.push_back(it.beta);
    if (result.buckling_constrained) {
      row.insert(row.end(), {it.min_load_factor, it.buckling_constraint,
                             static_cast<Scalar>(it.buckling_iterations)});
    }
    if (result.robust) {
      row.insert(row.end(), {it.eroded_volume_fraction, it.dilated_volume_fraction,
                             it.dilated_target_fraction});
    }
    csv.row(it.iteration, row);
  }
  csv.close();
}

void ResultWriter::write_density(const Mesh& mesh, const DesignDomain& domain,
                                 const TopologyOptimizationResult& result) const {
  const int dim = mesh.dim();
  std::vector<std::string> header{"element"};
  header = concat(header, coordinate_headers(dim, "c"));
  if (dim == 2) header.push_back("area[m2]");
  header = concat(header, {"volume[m3]", "design_x[-]", "physical_density[-]",
                           "stiffness_factor[-]", "passive_tag[-]", "strain_energy[J]"});
  const bool projected =
      result.projected && result.filtered_density.size() == mesh.num_elements();
  if (projected) header.push_back("filtered_density[-]");
  const bool printable =
      result.overhang_filtered && result.printable_density.size() == mesh.num_elements();
  if (printable) header.push_back("printable_density[-]");
  const bool robust = (result.robust || result.erosion_checked) &&
                      result.robust_record.eroded_density.size() == mesh.num_elements();
  if (robust) header.insert(header.end(), {"eroded_density[-]", "dilated_density[-]"});
  CsvWriter csv(file("density_final.csv"), header);
  for (Index e = 0; e < mesh.num_elements(); ++e) {
    const Vector3 c = mesh.element_centroid(e);
    std::vector<Scalar> row;
    for (int k = 0; k < dim; ++k) row.push_back(c(k));
    if (dim == 2) row.push_back(mesh.element_measure(e));
    row.push_back(domain.element_volumes()(e));
    row.push_back(result.design(e));
    row.push_back(result.physical_density(e));
    row.push_back(result.stiffness_factors(e));
    row.push_back(static_cast<Scalar>(
        static_cast<int>(domain.tags()[static_cast<std::size_t>(e)])));
    row.push_back(result.element_strain_energy.size() > e
                      ? result.element_strain_energy(e)
                      : 0.0);
    if (projected) row.push_back(result.filtered_density(e));
    if (printable) row.push_back(result.printable_density(e));
    if (robust) {
      row.push_back(result.robust_record.eroded_density(e));
      row.push_back(result.robust_record.dilated_density(e));
    }
    csv.row(e, row);
  }
  csv.close();
}

void ResultWriter::write_density_history(const TopologyOptimizationResult& result,
                                         int max_frames) const {
  if (result.snapshots.empty()) return;
  const int count = static_cast<int>(result.snapshots.size());
  const int stride = std::max(1, (count + max_frames - 1) / std::max(max_frames, 1));

  std::vector<int> frames;
  for (int i = 0; i < count; i += stride) frames.push_back(i);
  if (frames.back() != count - 1) frames.push_back(count - 1);

  const Eigen::Index ne = result.snapshots.front().size();
  std::vector<std::string> header{"iteration"};
  for (Eigen::Index e = 0; e < ne; ++e) {
    std::ostringstream os;
    os << "rho_" << e;
    header.push_back(os.str());
  }
  // Six significant digits keep the animation file small; the final density is
  // written at full precision in density_final.csv.
  CsvWriter csv(file("density_history.csv"), header, 6);
  for (int f : frames) {
    std::vector<Scalar> row;
    row.reserve(static_cast<std::size_t>(ne));
    const Vector& snapshot = result.snapshots[static_cast<std::size_t>(f)];
    for (Eigen::Index e = 0; e < ne; ++e) row.push_back(snapshot(e));
    csv.row(result.snapshot_iterations[static_cast<std::size_t>(f)], row);
  }
  csv.close();
  log::debug("wrote ", frames.size(), " density snapshots (of ", count,
             " recorded) to density_history.csv");
}

json::Value ResultWriter::write_geometry(const Mesh& mesh, const Vector& density,
                                        Scalar thickness, const std::string& stem,
                                        const std::string& what) const {
  if (density.size() != mesh.num_elements()) {
    std::ostringstream os;
    os << "geometry export '" << stem << "': density has " << density.size()
       << " entries but the mesh has " << mesh.num_elements() << " cells";
    throw IoError(os.str());
  }
  const std::string vtk_name = stem + ".vtk";
  const std::string stl_name = stem + ".stl";

  VtkWriter vtk(mesh, "SparLab " + what + ": " + config_.name);
  vtk.add_cell_scalars("density", density);
  vtk.write(file(vtk_name));

  const Scalar extrusion = mesh.dim() == 2 ? thickness : 1.0;
  const TriangleSurface surface = boundary_surface(mesh, extrusion);
  write_stl(file(stl_name), surface, config_.name + " " + what);
  const SurfaceStats stats = surface_stats(surface);

  // Closure is checked exactly on the edge pairing. The enclosed volume is
  // compared with the cell volume as well: identical for planar faces, and
  // different by the flat-triangle approximation of a bilinear patch where a
  // distorted cell was cut, which is reported so it is never mistaken for
  // exactness.
  Scalar cell_volume = 0.0;
  for (Index e = 0; e < mesh.num_elements(); ++e) {
    cell_volume += mesh.element_measure(e) * extrusion;
  }
  const Scalar mismatch =
      std::abs(stats.enclosed_volume - cell_volume) / std::max(cell_volume, 1.0e-300);
  if (!stats.closed) {
    log::warn("geometry export '", stem, "': the STL surface has ", stats.unmatched_edges,
              " unmatched directed edges, so it is not closed or not consistently "
              "oriented; do not use it as a solid");
  } else if (stats.enclosed_volume <= 0.0) {
    log::warn("geometry export '", stem, "': the STL surface encloses a non-positive "
              "volume (", stats.enclosed_volume, " m^3); its normals point inward");
  } else if (stats.non_manifold_edges > 0) {
    log::warn("geometry export '", stem, "': the STL surface is closed but has ",
              stats.non_manifold_edges, " non-manifold edge(s) where cells touch only "
              "along an edge; a slicer may split it into several shells, and the "
              "interpretation threshold or connectivity rule decides whether such "
              "cells are one part");
  }

  json::Value out = json::Value::make_object();
  out.set("vtk", json::Value::make_string(vtk_name));
  out.set("stl", json::Value::make_string(stl_name));
  out.set("stl_format", json::Value::make_string("binary, outward normals"));
  out.set("num_cells", json::Value::make_number(mesh.num_elements()));
  out.set("num_triangles", json::Value::make_number(stats.num_triangles));
  out.set("closed_surface", json::Value::make_bool(stats.closed));
  out.set("unmatched_edges", json::Value::make_number(stats.unmatched_edges));
  out.set("non_manifold_edges", json::Value::make_number(stats.non_manifold_edges));
  out.set("surface_area_m2", json::Value::make_number(stats.area));
  out.set("enclosed_volume_m3", json::Value::make_number(stats.enclosed_volume));
  out.set("cell_volume_m3", json::Value::make_number(cell_volume));
  out.set("volume_relative_mismatch", json::Value::make_number(mismatch));
  json::Value lower = json::Value::make_array();
  json::Value upper = json::Value::make_array();
  for (int i = 0; i < 3; ++i) {
    lower.push_back(json::Value::make_number(stats.bounds.lower(i)));
    upper.push_back(json::Value::make_number(stats.bounds.upper(i)));
  }
  out.set("bounds_lower_m", lower);
  out.set("bounds_upper_m", upper);
  if (mesh.dim() == 2) {
    out.set("extruded_thickness_m", json::Value::make_number(extrusion));
  }
  log::debug("wrote ", what, ": ", mesh.num_elements(), " cells, ", stats.num_triangles,
             " triangles, ", stats.enclosed_volume, " m^3");
  return out;
}

// ----------------------------------------------------------------------------
// Summaries
// ----------------------------------------------------------------------------

json::Value make_provenance(const Configuration& config) {
  json::Value out = json::Value::make_object();
  out.set("code", json::Value::make_string("SparLab"));
  out.set("version", json::Value::make_string(SPARLAB_VERSION_STRING));
  out.set("build_type", json::Value::make_string(SPARLAB_BUILD_TYPE));
  out.set("compiler", json::Value::make_string(SPARLAB_COMPILER));
  out.set("eigen_version", json::Value::make_string(SPARLAB_EIGEN_VERSION));
  out.set("config_source", json::Value::make_string(config.source_path));
  out.set("case", json::Value::make_string(config.name));
  out.set("description", json::Value::make_string(config.description));
  out.set("units", json::Value::make_string("SI: m, N, Pa, kg, kg/m^3, Hz, J"));

  const std::time_t now = std::time(nullptr);
  char buffer[64] = {0};
  std::tm utc {};
#if defined(_WIN32)
  gmtime_s(&utc, &now);
#else
  gmtime_r(&now, &utc);
#endif
  std::strftime(buffer, sizeof(buffer), "%Y-%m-%dT%H:%M:%SZ", &utc);
  out.set("timestamp_utc", json::Value::make_string(buffer));
  return out;
}

namespace {

json::Value amg_stats_json(const AmgStats& st) {
  json::Value out = json::Value::make_object();
  out.set("levels", json::Value::make_number(static_cast<Scalar>(st.levels.size())));
  json::Value levels = json::Value::make_array();
  for (const AmgLevelStats& l : st.levels) {
    json::Value lv = json::Value::make_object();
    lv.set("unknowns", json::Value::make_number(l.unknowns));
    lv.set("nonzeros", json::Value::make_number(l.nonzeros));
    lv.set("lambda_max_Dinv_A", json::Value::make_number(l.lambda_max));
    levels.push_back(lv);
  }
  out.set("level_sizes", levels);
  out.set("operator_complexity", json::Value::make_number(st.operator_complexity));
  out.set("grid_complexity", json::Value::make_number(st.grid_complexity));
  out.set("near_null_space_dimension",
          json::Value::make_number(st.near_null_space_dimension));
  out.set("coarse_pivot_ratio", json::Value::make_number(st.coarse_pivot_ratio));
  out.set("last_setup_seconds", json::Value::make_number(st.setup_seconds));
  out.set("aggregates_reused_in_last_setup", json::Value::make_bool(st.reused_aggregates));
  return out;
}

json::Value linear_solver_json(const LinearSolverOptions& o) {
  json::Value out = json::Value::make_object();
  out.set("type", json::Value::make_string(to_string(o.type)));
  out.set("warm_start", json::Value::make_bool(o.warm_start));
  if (o.type == LinearSolverType::Auto) {
    out.set("auto_direct_limit_plane", json::Value::make_number(o.auto_direct_limit_2d));
    out.set("auto_direct_limit_solid", json::Value::make_number(o.auto_direct_limit_3d));
  }
  if (o.type == LinearSolverType::AmgCg || o.type == LinearSolverType::Auto) {
    json::Value amg = json::Value::make_object();
    amg.set("strength_threshold", json::Value::make_number(o.amg.strength_threshold));
    amg.set("max_levels", json::Value::make_number(o.amg.max_levels));
    amg.set("coarse_size", json::Value::make_number(o.amg.coarse_size));
    amg.set("smoother", json::Value::make_string(to_string(o.amg.smoother)));
    amg.set("smoother_degree", json::Value::make_number(o.amg.smoother_degree));
    amg.set("chebyshev_ratio", json::Value::make_number(o.amg.chebyshev_ratio));
    amg.set("prolongator_damping", json::Value::make_number(o.amg.prolongator_damping));
    amg.set("lanczos_steps", json::Value::make_number(o.amg.lanczos_steps));
    amg.set("reuse_aggregates", json::Value::make_bool(o.amg.reuse_aggregates));
    amg.set("coarse_pivot_tolerance",
            json::Value::make_number(o.amg.coarse_pivot_tolerance));
    out.set("multigrid", amg);
  }
  return out;
}

json::Value tolerance_json(const Configuration& config) {
  json::Value out = json::Value::make_object();
  out.set("linear_solver", json::Value::make_string(
                               to_string(config.analysis.linear.type)));
  out.set("linear_residual_tolerance",
          json::Value::make_number(config.analysis.linear.residual_tolerance));
  out.set("linear_iterative_tolerance",
          json::Value::make_number(config.analysis.linear.iterative_tolerance));
  out.set("cholesky_pivot_tolerance",
          json::Value::make_number(config.analysis.linear.pivot_tolerance));
  out.set("linear_solver_settings", linear_solver_json(config.analysis.linear));
  out.set("equilibrium_tolerance",
          json::Value::make_number(config.analysis.equilibrium_tolerance));
  out.set("modal_tolerance", json::Value::make_number(config.modal.options.tolerance));
  out.set("modal_residual_tolerance",
          json::Value::make_number(config.modal.options.residual_tolerance));
  out.set("buckling_tolerance", json::Value::make_number(config.buckling.options.tolerance));
  out.set("buckling_residual_tolerance",
          json::Value::make_number(config.buckling.options.residual_tolerance));
  out.set("optimizer_change_tolerance",
          json::Value::make_number(config.topology.optimizer.change_tolerance));
  out.set("optimizer_volume_tolerance",
          json::Value::make_number(config.topology.optimizer.oc.volume_tolerance));
  return out;
}

}  // namespace

json::Value make_static_summary(const Configuration& config, const FemModel& model,
                                const ModelDiagnostics& diagnostics,
                                const std::vector<StaticSolution>& solutions,
                                const std::vector<StressField>& stresses,
                                const ModalResult* modal,
                                const TimingLedger& timings) {
  json::Value out = json::Value::make_object();
  out.set("provenance", make_provenance(config));
  out.set("mesh", mesh_stats_json(config, model));
  out.set("material", material_json(model.material(), model.stress_state()));
  if (!model.single_material()) out.set("materials", materials_json(model));
  out.set("tolerances", tolerance_json(config));
  out.set("model_diagnostics", diagnostics_json(diagnostics));

  json::Value cases = json::Value::make_array();
  for (std::size_t l = 0; l < solutions.size(); ++l) {
    const StaticSolution& sol = solutions[l];
    json::Value entry = json::Value::make_object();
    entry.set("name", json::Value::make_string(sol.load_case_name));
    {
      const json::Value parts = load_parts_json(model, l);
      if (!parts.members().empty()) entry.set("loads", parts);
    }
    entry.set("weight", json::Value::make_number(sol.weight));
    entry.set("compliance_J", json::Value::make_number(sol.compliance));
    entry.set("strain_energy_J", json::Value::make_number(sol.strain_energy));
    entry.set("compliance_over_twice_strain_energy",
              json::Value::make_number(
                  sol.strain_energy != 0.0
                      ? sol.compliance / (2.0 * sol.strain_energy)
                      : 0.0));
    entry.set("max_displacement_magnitude_m",
              json::Value::make_number(sol.max_displacement_magnitude));
    entry.set("max_displacement_node",
              json::Value::make_number(sol.max_displacement_node));
    entry.set("scaled_residual", json::Value::make_number(sol.scaled_residual));
    entry.set("backward_error", json::Value::make_number(sol.backward_error));
    entry.set("solver_iterations", json::Value::make_number(sol.solver_iterations));
    entry.set("linear_solver", json::Value::make_string(sol.solver_name));
    entry.set("equilibrium", equilibrium_json(sol.equilibrium, model.dim()));
    if (l < stresses.size()) {
      entry.set("max_von_mises_Pa",
                json::Value::make_number(stresses[l].element_von_mises.maxCoeff()));
      entry.set("total_strain_energy_J",
                json::Value::make_number(stresses[l].element_strain_energy.sum()));
    }
    cases.push_back(entry);
  }
  out.set("load_cases", cases);

  if (modal != nullptr) out.set("modal", modal_json(*modal));
  out.set("timings_s", timings_json(timings));
  return out;
}

json::Value buckling_json(const std::vector<BucklingResult>& results,
                          const BucklingOptions& options, const std::string& what) {
  json::Value out = json::Value::make_object();
  out.set("structure", json::Value::make_string(what));
  out.set("method", json::Value::make_string(
                        "linear buckling (K_ff + lambda K_G,ff(u)) phi = 0 of the linear "
                        "static stress state of each load case, by subspace iteration on "
                        "(-K_G, K) with a Rayleigh-Ritz projection (dense eigensolve up "
                        "to 400 free DOFs)"));
  out.set("meaning", json::Value::make_string(
                         "the structure is predicted to buckle at load_factor times the "
                         "load case; an upper bound for a real, imperfect structure "
                         "(bifurcation of the ideal geometry, no follower forces, no "
                         "post-buckling)"));
  out.set("tolerance", json::Value::make_number(options.tolerance));
  out.set("residual_tolerance", json::Value::make_number(options.residual_tolerance));
  out.set("requested_modes", json::Value::make_number(options.num_modes));
  json::Value cases = json::Value::make_array();
  for (const BucklingResult& r : results) {
    json::Value c = json::Value::make_object();
    c.set("load_case", json::Value::make_string(r.load_case));
    c.set("load_factors", json::array_of(r.load_factors));
    c.set("eigenpair_residuals", json::array_of(r.residuals));
    c.set("solid_energy_fraction", json::array_of(r.solid_energy_fraction));
    c.set("no_positive_load_factor", json::Value::make_bool(r.no_positive_load_factor));
    c.set("iterations", json::Value::make_number(r.iterations));
    c.set("subspace_size", json::Value::make_number(r.subspace_size));
    c.set("spectral_transformation", json::Value::make_bool(r.transformed));
    if (r.transformed) {
      c.set("transformation_shift", json::Value::make_number(r.sigma));
      c.set("shifted_factorizations", json::Value::make_number(r.factorizations));
    }
    c.set("converged", json::Value::make_bool(r.converged));
    c.set("final_relative_change", json::Value::make_number(r.final_change));
    if (!r.linear_solver.empty()) {
      c.set("linear_solver", json::Value::make_string(r.linear_solver));
    }
    c.set("warnings", json::array_of(r.warnings));
    cases.push_back(c);
  }
  out.set("load_cases", cases);
  return out;
}

json::Value nonlinear_json(const std::vector<NonlinearResult>& results,
                           const NonlinearOptions& options, const FemModel& model,
                           const std::vector<StaticSolution>& linear) {
  const bool arc = options.method == NonlinearOptions::Method::ArcLength;
  const bool small = options.kinematics == Kinematics::SmallStrain;
  const bool svk = options.law == HyperelasticModel::SaintVenantKirchhoff;
  bool plastic = false;
  for (Index e = 0; e < model.mesh().num_elements(); ++e) {
    if (model.material_of(e).plasticity().enabled()) plastic = true;
  }
  json::Value out = json::Value::make_object();
  std::string formulation =
      small ? "small-strain statics: linear strain on the undeformed geometry, Newton's method "
              "with the consistent tangent"
            : "geometrically non-linear statics, total Lagrangian: second Piola-Kirchhoff "
              "stress and Green-Lagrange strain on the reference configuration, consistent "
              "tangent (material plus initial-stress part), Newton's method";
  if (plastic) {
    formulation += small ? "; J2 plasticity by the backward-Euler radial return"
                         : "; J2 plasticity in the Green-Lagrange strain and second "
                           "Piola-Kirchhoff stress (small strain, large rotation) by the "
                           "backward-Euler radial return";
  }
  out.set("formulation", json::Value::make_string(formulation));
  out.set("kinematics", json::Value::make_string(to_string(options.kinematics)));
  out.set("material_model", json::Value::make_string(small ? "linear_elastic"
                                                           : to_string(options.law)));
  out.set("method", json::Value::make_string(to_string(options.method)));
  out.set("follower_pressure", json::Value::make_bool(!small && options.follower_pressure));
  out.set("load_scaling",
          json::Value::make_string(
              small ? "every load of the case - forces, pressures, body forces, the rotation's "
                      "centrifugal load, a temperature change and prescribed displacements - "
                      "scales with the load factor lambda and acts on the undeformed geometry"
                    : "every load of the case - forces, pressures, body forces, the rotation's "
                      "centrifugal load, a temperature change and prescribed displacements - "
                      "scales with the load factor lambda; pressures follow the deformed faces "
                      "when follower_pressure is set, a rotation acts at the deformed "
                      "position"));
  if (plastic) {
    json::Value p = json::Value::make_object();
    p.set("law", json::Value::make_string(
                     "J2 (von Mises) with linear and Voce isotropic hardening and Prager's "
                     "linear kinematic hardening; backward-Euler radial return with the "
                     "consistent tangent"));
    p.set("mean_dilatation", json::Value::make_string(to_string(options.mean_dilatation)));
    out.set("plasticity", p);
  }
  json::Value opts = json::Value::make_object();
  opts.set("steps", json::Value::make_number(options.steps));
  opts.set("max_steps", json::Value::make_number(options.max_steps));
  opts.set("max_iterations", json::Value::make_number(options.max_iterations));
  opts.set("max_cuts", json::Value::make_number(options.max_cuts));
  opts.set("residual_tolerance", json::Value::make_number(options.residual_tolerance));
  opts.set("displacement_tolerance", json::Value::make_number(options.displacement_tolerance));
  opts.set("line_search", json::Value::make_bool(options.line_search));
  if (!options.load_path.empty()) opts.set("load_path", json::array_of(options.load_path));
  if (arc) {
    opts.set("target_load_factor", json::Value::make_number(options.target_load_factor));
    opts.set("desired_iterations", json::Value::make_number(options.desired_iterations));
    opts.set("min_arc_ratio", json::Value::make_number(options.min_arc_ratio));
    opts.set("max_arc_ratio", json::Value::make_number(options.max_arc_ratio));
  }
  out.set("options", opts);

  const int dim = model.dim();
  json::Value cases = json::Value::make_array();
  for (const NonlinearResult& r : results) {
    json::Value c = json::Value::make_object();
    c.set("load_case", json::Value::make_string(r.load_case_name));
    c.set("completed", json::Value::make_bool(r.completed));
    c.set("termination", json::Value::make_string(r.termination));
    c.set("load_factor", json::Value::make_number(r.load_factor));
    if (!std::isnan(r.critical_bound)) {
      json::Value bracket = json::Value::make_array();
      bracket.push_back(json::Value::make_number(r.load_factor));
      bracket.push_back(json::Value::make_number(r.critical_bound));
      c.set("critical_load_factor_bracket", bracket);
    }
    if (!std::isnan(r.unreached_load_factor)) {
      c.set("unreached_load_factor", json::Value::make_number(r.unreached_load_factor));
    }
    c.set("steps", json::Value::make_number(static_cast<Scalar>(r.steps.size())));
    c.set("iterations", json::Value::make_number(r.total_iterations));
    c.set("cuts", json::Value::make_number(r.total_cuts));
    Scalar max_u = 0.0;
    for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
      max_u = std::max(max_u, magnitude(r.displacement, n, dim));
    }
    c.set("max_displacement_m", json::Value::make_number(max_u));
    for (const StaticSolution& s : linear) {
      if (s.load_case_name != r.load_case_name) continue;
      // The linear solution at the same load factor (it is linear in lambda).
      c.set("linear_max_displacement_m",
            json::Value::make_number(r.load_factor * s.max_displacement_magnitude));
    }
    c.set("strain_energy_J", json::Value::make_number(r.strain_energy));
    if (small) {
      c.set("max_strain", json::Value::make_number(r.max_green_strain));
      c.set("max_rotation_rad", json::Value::make_number(r.max_rotation));
      c.set("max_neglected_quadratic_strain", json::Value::make_number(r.max_quadratic_strain));
    } else {
      c.set("max_green_strain", json::Value::make_number(r.max_green_strain));
      c.set("min_jacobian", json::Value::make_number(r.min_jacobian));
    }
    if (r.plastic) {
      json::Value p = json::Value::make_object();
      p.set("mean_dilatation_applied", json::Value::make_bool(r.mean_dilatation));
      p.set("yielded_points", json::Value::make_number(r.plastic_points));
      p.set("elastoplastic_points", json::Value::make_number(r.total_points));
      p.set("max_equivalent_plastic_strain", json::Value::make_number(r.max_plastic_strain));
      int yielded_elements = 0;
      for (Eigen::Index e = 0; e < r.element_plastic_strain.size(); ++e) {
        if (r.element_plastic_strain(e) > 0.0) ++yielded_elements;
      }
      p.set("yielded_elements", json::Value::make_number(yielded_elements));
      int first_yield_step = -1;
      for (const NonlinearStep& s : r.steps) {
        if (s.yielding_points > 0) {
          first_yield_step = s.index;
          break;
        }
      }
      // The step in which a point first yielded, and the load factor it
      // ended at (yield began within it).
      p.set("first_yielding_step", json::Value::make_number(first_yield_step));
      if (first_yield_step > 0) {
        p.set("first_yielding_step_load_factor",
              json::Value::make_number(
                  r.steps[static_cast<std::size_t>(first_yield_step - 1)].load_factor));
      }
      c.set("plasticity", p);
    }
    c.set("max_von_mises_Pa", json::Value::make_number(
                                  r.element_von_mises.size() > 0 ? r.element_von_mises.maxCoeff()
                                                                 : 0.0));
    Scalar min_lambda = 0.0;
    Scalar max_lambda = 0.0;
    int max_pivots = -1;
    for (const NonlinearStep& s : r.steps) {
      min_lambda = std::min(min_lambda, s.load_factor);
      max_lambda = std::max(max_lambda, s.load_factor);
      max_pivots = std::max(max_pivots, s.negative_pivots);
    }
    c.set("max_load_factor_on_path", json::Value::make_number(max_lambda));
    c.set("min_load_factor_on_path", json::Value::make_number(min_lambda));
    const int final_pivots = r.steps.empty() ? -1 : r.steps.back().negative_pivots;
    c.set("final_negative_pivots", json::Value::make_number(final_pivots));
    c.set("max_negative_pivots_on_path", json::Value::make_number(max_pivots));
    json::Value monitors = json::Value::make_object();
    if (!r.steps.empty()) {
      for (std::size_t i = 0; i < r.monitor_names.size(); ++i) {
        monitors.set(r.monitor_names[i] + "_" + r.monitor_units[i],
                     json::Value::make_number(r.steps.back().monitors[i]));
      }
    }
    c.set("final_monitors", monitors);
    c.set("equilibrium", equilibrium_json(r.equilibrium, dim));
    c.set("symmetric_tangent", json::Value::make_bool(r.symmetric_tangent));
    c.set("linear_solver", json::Value::make_string(r.linear_solver));

    std::vector<std::string> warnings;
    if (!r.completed) warnings.push_back("the run stopped early: " + r.termination);
    for (const std::string& w : r.warnings) warnings.push_back(w);
    if (final_pivots > 0) {
      warnings.push_back("the final state is unstable: its tangent has " +
                         std::to_string(final_pivots) +
                         " negative eigenvalue(s), so the path has passed a limit or "
                         "bifurcation point");
    } else if (max_pivots > 0) {
      warnings.push_back("the path passed through unstable states (up to " +
                         std::to_string(max_pivots) +
                         " negative eigenvalue(s) of the tangent) before the final one");
    }
    if (!r.symmetric_tangent) {
      warnings.push_back("the tangent is non-symmetric (follower pressure), factorised by LU, "
                         "which reveals no inertia: stability is not assessed");
    }
    if (!small && !r.plastic && svk && r.max_green_strain > 0.05) {
      std::ostringstream os;
      os << "the largest Green-Lagrange strain is " << r.max_green_strain
         << "; the Saint Venant-Kirchhoff law is meant for small strain (large rotation) - "
            "use \"neo_hookean\" for large strain";
      warnings.push_back(os.str());
    }
    if (!small && svk && r.min_jacobian < 1.0 / std::sqrt(3.0)) {
      std::ostringstream os;
      os << "an integration point is compressed to a volume ratio J = " << r.min_jacobian
         << " < 1/sqrt(3): below a stretch of 1/sqrt(3) the Saint Venant-Kirchhoff law's "
            "compressive force falls again, so it does not model strong compression - use "
            "\"neo_hookean\"";
      warnings.push_back(os.str());
    }
    c.set("warnings", json::array_of(warnings));
    cases.push_back(c);
  }
  out.set("load_cases", cases);
  return out;
}

namespace {

json::Value amplitude_json(const Amplitude& a) {
  json::Value out = json::Value::make_object();
  out.set("type", json::Value::make_string(to_string(a.kind)));
  out.set("scale", json::Value::make_number(a.scale));
  if (a.kind == Amplitude::Kind::Table) {
    out.set("times_s", json::array_of(a.times));
    out.set("values", json::array_of(a.values));
  } else if (a.kind == Amplitude::Kind::Harmonic) {
    out.set("frequency_Hz", json::Value::make_number(a.frequency));
    out.set("phase_rad", json::Value::make_number(a.phase));
  }
  return out;
}

}  // namespace

json::Value transient_json(const std::vector<TransientResult>& results,
                           const TransientOptions& options) {
  const HhtParameters p = HhtParameters::from_alpha(options.alpha);
  const NonlinearOptions& nl = options.nonlinear_options;
  json::Value out = json::Value::make_object();
  out.set("method",
          json::Value::make_string(
              "HHT-alpha (Hilber, Hughes and Taylor 1977) with a constant time step: Newmark "
              "kinematics, equilibrium at the weighted point M a1 + (1 + alpha)(C v1 + K u1) - "
              "alpha (C v0 + K u0) = (1 + alpha) f1 - alpha f0; alpha = 0 is the trapezoidal "
              "rule (average acceleration), which conserves the energy of a linear undamped "
              "model exactly"));
  out.set("equations",
          json::Value::make_string(
              options.nonlinear
                  ? "M a + C v + f_int(u) = f_ext(A(t)): the internal forces, loads and tangent "
                    "of the load case's non-linear system, Newton's method at every step on the "
                    "HHT-alpha residual with the tangent c0 M + (1 + alpha)(c3 C + K_T), the "
                    "plastic history committed on convergence"
                  : "M a + C v + K u = A(t) f, the prescribed displacements A(t) g; the "
                    "effective stiffness c0 M + (1 + alpha)(c3 C + K) factorised once"));
  out.set("alpha", json::Value::make_number(p.alpha));
  out.set("beta", json::Value::make_number(p.beta));
  out.set("gamma", json::Value::make_number(p.gamma));
  out.set("time_step_s", json::Value::make_number(options.time_step));
  out.set("end_time_s", json::Value::make_number(options.end_time));
  out.set("mass", json::Value::make_string(to_string(options.mass_type)));
  json::Value damping = json::Value::make_object();
  damping.set("mass_1_per_s", json::Value::make_number(options.mass_damping));
  damping.set("stiffness_s", json::Value::make_number(options.stiffness_damping));
  damping.set("stiffness_matrix",
              json::Value::make_string(options.nonlinear
                                           ? "the linear elastic stiffness of the undeformed model"
                                           : "the model's stiffness"));
  out.set("rayleigh_damping", damping);
  out.set("amplitude", amplitude_json(options.amplitude));
  out.set("start", json::Value::make_string(to_string(options.start)));
  out.set("snapshot_every", json::Value::make_number(options.snapshot_every));
  out.set("linear_solver_settings", linear_solver_json(options.linear));
  out.set("energy_balance",
          json::Value::make_string(
              "E_0 + W - T - U - D: the initial energy plus the trapezoidal work of the loads "
              "and the reactions of the prescribed motion, less the kinetic, strain (stored) "
              "and damping energies; zero to round-off for the trapezoidal rule on a linear "
              "model, the numerical dissipation of alpha < 0, and it holds the plastic "
              "dissipation of a plastic model"));
  if (options.nonlinear) {
    const bool small = nl.kinematics == Kinematics::SmallStrain;
    json::Value n = json::Value::make_object();
    n.set("kinematics", json::Value::make_string(to_string(nl.kinematics)));
    n.set("material_model",
          json::Value::make_string(small ? "linear_elastic" : to_string(nl.law)));
    n.set("mean_dilatation", json::Value::make_string(to_string(nl.mean_dilatation)));
    n.set("follower_pressure", json::Value::make_bool(!small && nl.follower_pressure));
    n.set("residual_tolerance", json::Value::make_number(nl.residual_tolerance));
    n.set("displacement_tolerance", json::Value::make_number(nl.displacement_tolerance));
    n.set("max_iterations", json::Value::make_number(nl.max_iterations));
    n.set("max_cuts", json::Value::make_number(options.max_cuts));
    out.set("nonlinear", n);
  }
  json::Value cases = json::Value::make_array();
  for (const TransientResult& r : results) {
    json::Value c = json::Value::make_object();
    c.set("load_case", json::Value::make_string(r.load_case_name));
    c.set("completed", json::Value::make_bool(r.completed));
    if (!r.completed) c.set("termination", json::Value::make_string(r.termination));
    c.set("steps", json::Value::make_number(static_cast<Scalar>(r.steps.size()) - 1.0));
    c.set("requested_steps", json::Value::make_number(r.num_steps));
    c.set("total_mass_kg", json::Value::make_number(r.total_mass));
    if (!r.steps.empty()) {
      const TransientStep& first = r.steps.front();
      const TransientStep& last = r.steps.back();
      c.set("final_time_s", json::Value::make_number(last.time));
      std::size_t peak = 0;
      for (std::size_t i = 1; i < r.steps.size(); ++i) {
        if (r.steps[i].max_displacement > r.steps[peak].max_displacement) peak = i;
      }
      c.set("max_displacement_m", json::Value::make_number(r.steps[peak].max_displacement));
      c.set("max_displacement_time_s", json::Value::make_number(r.steps[peak].time));
      c.set("final_max_displacement_m", json::Value::make_number(last.max_displacement));
      json::Value energy = json::Value::make_object();
      energy.set("initial_J", json::Value::make_number(first.kinetic_energy + first.strain_energy));
      energy.set("final_kinetic_J", json::Value::make_number(last.kinetic_energy));
      energy.set(r.plastic ? "final_stored_J" : "final_strain_J",
                 json::Value::make_number(last.strain_energy));
      energy.set("damping_dissipation_J", json::Value::make_number(last.damping_energy));
      energy.set("external_work_J", json::Value::make_number(last.external_work));
      energy.set("final_balance_J", json::Value::make_number(r.numerical_dissipation));
      energy.set("largest_relative_balance", json::Value::make_number(r.energy_balance_error));
      c.set("energy", energy);
      json::Value monitors = json::Value::make_array();
      for (std::size_t i = 0; i < r.monitor_names.size(); ++i) {
        std::size_t lo = 0;
        std::size_t hi = 0;
        for (std::size_t s = 1; s < r.steps.size(); ++s) {
          if (r.steps[s].monitors[i] < r.steps[lo].monitors[i]) lo = s;
          if (r.steps[s].monitors[i] > r.steps[hi].monitors[i]) hi = s;
        }
        json::Value m = json::Value::make_object();
        m.set("name", json::Value::make_string(r.monitor_names[i]));
        m.set("unit", json::Value::make_string(r.monitor_units[i]));
        if (i < options.monitors.size()) {
          m.set("quantity", json::Value::make_string(to_string(options.monitors[i].quantity)));
          m.set("component", json::Value::make_number(options.monitors[i].component));
        }
        if (i < r.monitor_nodes.size()) m.set("nodes", json::array_of(r.monitor_nodes[i]));
        m.set("max", json::Value::make_number(r.steps[hi].monitors[i]));
        m.set("max_time_s", json::Value::make_number(r.steps[hi].time));
        m.set("min", json::Value::make_number(r.steps[lo].monitors[i]));
        m.set("min_time_s", json::Value::make_number(r.steps[lo].time));
        m.set("final", json::Value::make_number(last.monitors[i]));
        monitors.push_back(m);
      }
      c.set("monitors", monitors);
    }
    if (r.nonlinear) {
      int most = 0;
      int cuts = 0;
      for (const TransientStep& s : r.steps) {
        most = std::max(most, s.iterations);
        cuts += s.cuts;
      }
      c.set("newton_iterations", json::Value::make_number(r.total_iterations));
      c.set("max_iterations_in_a_step", json::Value::make_number(most));
      c.set("step_halvings", json::Value::make_number(cuts));
      if (r.plastic) {
        json::Value pl = json::Value::make_object();
        pl.set("mean_dilatation_applied", json::Value::make_bool(r.mean_dilatation));
        pl.set("max_equivalent_plastic_strain", json::Value::make_number(r.max_plastic_strain));
        int first_yield_step = -1;
        for (const TransientStep& s : r.steps) {
          if (s.yielding_points > 0) {
            first_yield_step = s.index;
            break;
          }
        }
        // The step in which a point first yielded, and the interval of time
        // it spans (yield began within it).
        pl.set("first_yielding_step", json::Value::make_number(first_yield_step));
        if (first_yield_step > 0) {
          const auto k = static_cast<std::size_t>(first_yield_step);
          json::Value interval = json::Value::make_array();
          interval.push_back(json::Value::make_number(r.steps[k - 1].time));
          interval.push_back(json::Value::make_number(r.steps[k].time));
          pl.set("first_yielding_interval_s", interval);
        }
        c.set("plasticity", pl);
      }
    }
    c.set("linear_solver", json::Value::make_string(r.linear_solver));
    c.set("warnings", json::array_of(r.warnings));
    cases.push_back(c);
  }
  out.set("load_cases", cases);
  return out;
}

json::Value frequency_response_json(const std::vector<FrequencyResponseResult>& results,
                                    const FrequencyResponseOptions& options) {
  constexpr Scalar kDegrees = 180.0 / 3.14159265358979323846;
  json::Value out = json::Value::make_object();
  out.set("method",
          json::Value::make_string(
              "steady-state harmonic response: [K (1 + i eta) - omega^2 M + i omega C] U = f on "
              "the free DOFs, C = a M + b K, a complex sparse LU per frequency; the load case's "
              "loads and prescribed displacements are the amplitudes of f cos(omega t) and "
              "g cos(omega t), the response is Re(U e^{i omega t}) and a phase is arg U (the "
              "lag behind the load is -arg U)"));
  out.set("frequencies", json::Value::make_number(static_cast<Scalar>(options.frequencies.size())));
  out.set("frequencies_Hz", json::array_of(options.frequencies));
  if (!options.frequencies.empty()) {
    const auto range = std::minmax_element(options.frequencies.begin(), options.frequencies.end());
    out.set("min_frequency_Hz", json::Value::make_number(*range.first));
    out.set("max_frequency_Hz", json::Value::make_number(*range.second));
  }
  out.set("mass", json::Value::make_string(to_string(options.mass_type)));
  json::Value damping = json::Value::make_object();
  damping.set("mass_1_per_s", json::Value::make_number(options.mass_damping));
  damping.set("stiffness_s", json::Value::make_number(options.stiffness_damping));
  damping.set("structural_loss_factor", json::Value::make_number(options.structural_damping));
  out.set("damping", damping);
  json::Value cases = json::Value::make_array();
  for (const FrequencyResponseResult& r : results) {
    json::Value c = json::Value::make_object();
    c.set("load_case", json::Value::make_string(r.load_case_name));
    if (!r.points.empty()) {
      std::size_t peak = 0;
      for (std::size_t j = 1; j < r.points.size(); ++j) {
        if (r.points[j].max_displacement > r.points[peak].max_displacement) peak = j;
      }
      c.set("max_displacement_m", json::Value::make_number(r.points[peak].max_displacement));
      c.set("max_displacement_frequency_Hz", json::Value::make_number(r.points[peak].frequency));
      json::Value monitors = json::Value::make_array();
      for (std::size_t i = 0; i < r.monitor_names.size(); ++i) {
        std::size_t top = 0;
        for (std::size_t j = 1; j < r.points.size(); ++j) {
          if (std::abs(r.points[j].monitors[i]) > std::abs(r.points[top].monitors[i])) top = j;
        }
        const ComplexScalar z = r.points[top].monitors[i];
        json::Value m = json::Value::make_object();
        m.set("name", json::Value::make_string(r.monitor_names[i]));
        m.set("unit", json::Value::make_string(r.monitor_units[i]));
        if (i < options.monitors.size()) {
          m.set("quantity", json::Value::make_string(to_string(options.monitors[i].quantity)));
          m.set("component", json::Value::make_number(options.monitors[i].component));
        }
        if (i < r.monitor_nodes.size()) m.set("nodes", json::array_of(r.monitor_nodes[i]));
        m.set("peak_amplitude", json::Value::make_number(std::abs(z)));
        m.set("peak_frequency_Hz", json::Value::make_number(r.points[top].frequency));
        m.set("phase_at_peak_deg", json::Value::make_number(std::arg(z) * kDegrees));
        monitors.push_back(m);
      }
      c.set("monitors", monitors);
    }
    std::vector<Scalar> snapshots;
    for (const FrequencySnapshot& s : r.snapshots) snapshots.push_back(s.frequency);
    c.set("snapshot_frequencies_Hz", json::array_of(snapshots));
    c.set("warnings", json::array_of(r.warnings));
    cases.push_back(c);
  }
  out.set("load_cases", cases);
  return out;
}

json::Value overhang_json(const OverhangReport& report, bool filtered) {
  json::Value out = json::Value::make_object();
  out.set("build_direction", json::Value::make_string(report.build_direction));
  out.set("overhang_angle_deg", json::Value::make_number(report.overhang_angle_degrees));
  out.set("threshold", json::Value::make_number(report.threshold));
  out.set("filtered_during_optimisation", json::Value::make_bool(filtered));
  out.set("solid_elements", json::Value::make_number(report.solid_elements));
  out.set("unsupported_elements", json::Value::make_number(report.unsupported_elements));
  out.set("solid_volume_m3", json::Value::make_number(report.solid_volume));
  out.set("unsupported_volume_m3", json::Value::make_number(report.unsupported_volume));
  out.set("unsupported_fraction", json::Value::make_number(report.unsupported_fraction));
  out.set("lowest_unsupported_layer", json::Value::make_number(report.lowest_unsupported_layer));
  out.set("note", json::Value::make_string(
                      "a solid element (density >= threshold) off the build plate is "
                      "supported when a solid element lies in its stencil in the layer "
                      "below: directly below or diagonally (3 elements in 2-D, a cross of 5 "
                      "in 3-D). Only this overhang rule is modelled; no other process "
                      "limit (support removal, thermal distortion, minimum wall) is"));
  return out;
}

json::Value length_scale_json(const LengthScaleScan& scan) {
  json::Value out = json::Value::make_object();
  out.set("solid_min_size_m", json::Value::make_number(scan.solid_min_size));
  out.set("void_min_size_m", json::Value::make_number(scan.void_min_size));
  out.set("solid_bound_reached_scan_cap", json::Value::make_bool(scan.solid_bound_reached_cap));
  out.set("void_bound_reached_scan_cap", json::Value::make_bool(scan.void_bound_reached_cap));
  out.set("radius_step_m", json::Value::make_number(scan.step));
  out.set("tolerance", json::Value::make_number(scan.tolerance));
  json::Value probes = json::Value::make_array();
  for (const LengthScaleReport& p : scan.probes) {
    json::Value e = json::Value::make_object();
    e.set("radius_m", json::Value::make_number(p.radius));
    e.set("solid_fraction_removed_by_opening",
          json::Value::make_number(p.solid_violation_fraction));
    e.set("void_fraction_filled_by_closing", json::Value::make_number(p.void_violation_fraction));
    probes.push_back(e);
  }
  out.set("probes", probes);
  out.set("note", json::Value::make_string(
                      "morphological opening (solid) and closing (void) of the thresholded "
                      "design with balls of growing radius r on the element centroids: a "
                      "probe flags members thinner (gaps narrower) than about 2r. The minimum "
                      "sizes are twice the largest radius flagging at most `tolerance` of the "
                      "solid (void) volume, to within one radius step; a measurement of the "
                      "design, not a guarantee"));
  return out;
}

json::Value make_topology_summary(const Configuration& config, const FemModel& model,
                                  const DesignDomain& domain,
                                  const DensityFilter& filter,
                                  const TopologyOptimizationResult& result,
                                  const TopologyInterpretation* interpretation,
                                  const ModalResult* modal_initial,
                                  const ModalResult* modal_optimised,
                                  const ModalResult* modal_mass_matched,
                                  const TimingLedger& timings) {
  json::Value out = json::Value::make_object();
  out.set("provenance", make_provenance(config));
  out.set("mesh", mesh_stats_json(config, model));
  out.set("material", material_json(model.material(), model.stress_state()));
  out.set("tolerances", tolerance_json(config));

  {
    json::Value setup = json::Value::make_object();
    // The objective is a weighted mean over load cases, so the summary has to
    // name the cases and record their weights: `load_case_compliance_J` below
    // is in this order, and without the weights the reported total cannot be
    // reconstructed from the parts. Both the weights as written in the deck and
    // the normalised weights actually used are recorded, because only the
    // latter reproduce the objective:
    //   compliance_J = sum_l load_case_weights_normalised[l] *
    //                        load_case_compliance_J[l].
    {
      std::vector<std::string> names;
      std::vector<Scalar> weights;
      names.reserve(model.load_case_specs().size());
      weights.reserve(model.load_case_specs().size());
      for (const auto& spec : model.load_case_specs()) {
        names.push_back(spec.name);
        weights.push_back(spec.weight);
      }
      setup.set("load_case_names", json::array_of(names));
      setup.set("load_case_weights", json::array_of(weights));
      setup.set("load_case_weights_normalised",
                json::array_of(model.normalised_weights()));
    }
    setup.set("volume_fraction_target",
              json::Value::make_number(domain.volume_fraction()));
    setup.set("volume_target_m3", json::Value::make_number(domain.volume_target()));
    setup.set("initial_density", json::Value::make_number(domain.initial_density()));
    setup.set("domain_volume_m3", json::Value::make_number(domain.domain_volume()));
    setup.set("num_design_variables", json::Value::make_number(domain.num_elements()));
    setup.set("num_free_variables",
              json::Value::make_number(domain.num_free_variables()));
    setup.set("num_passive_solid", json::Value::make_number(domain.num_passive_solid()));
    setup.set("num_passive_void", json::Value::make_number(domain.num_passive_void()));
    setup.set("passive_solid_volume_m3",
              json::Value::make_number(domain.passive_solid_volume()));
    setup.set("passive_void_volume_m3",
              json::Value::make_number(domain.passive_void_volume()));
    setup.set("filter_type", json::Value::make_string(to_string(filter.type())));
    setup.set("filter_radius_m", json::Value::make_number(filter.radius()));
    setup.set("filter_average_support",
              json::Value::make_number(filter.average_support()));
    setup.set("simp_penalty",
              json::Value::make_number(config.topology.optimizer.simp.penalty));
    setup.set("simp_emin_ratio",
              json::Value::make_number(config.topology.optimizer.simp.emin_ratio));
    setup.set("mass_interpolation",
              json::Value::make_string(
                  to_string(config.topology.optimizer.simp.mass_law)));
    setup.set("move_limit",
              json::Value::make_number(config.topology.optimizer.oc.move_limit));
    setup.set("damping", json::Value::make_number(config.topology.optimizer.oc.damping));
    setup.set("continuation_steps",
              json::Value::make_number(config.topology.optimizer.continuation_steps));
    setup.set("penalty_start",
              json::Value::make_number(config.topology.optimizer.penalty_start));
    setup.set("max_iterations",
              json::Value::make_number(config.topology.optimizer.max_iterations));
    setup.set("change_tolerance",
              json::Value::make_number(config.topology.optimizer.change_tolerance));
    setup.set("objective_tolerance",
              json::Value::make_number(config.topology.optimizer.objective_tolerance));
    setup.set("objective_window",
              json::Value::make_number(config.topology.optimizer.objective_window));
    setup.set("method", json::Value::make_string(to_string(result.method)));
    if (result.method == OptimizerMethod::MMA) {
      json::Value mma = json::Value::make_object();
      const MmaOptions& mo = config.topology.optimizer.mma;
      mma.set("move_limit", json::Value::make_number(mo.move_limit));
      mma.set("asymptote_init", json::Value::make_number(mo.asymptote_init));
      mma.set("asymptote_increase", json::Value::make_number(mo.asymptote_increase));
      mma.set("asymptote_decrease", json::Value::make_number(mo.asymptote_decrease));
      mma.set("constraint_penalty", json::Value::make_number(mo.c));
      mma.set("subproblem_tolerance", json::Value::make_number(mo.epsimin));
      mma.set("constraint_tolerance",
              json::Value::make_number(config.topology.optimizer.constraint_tolerance));
      setup.set("mma", mma);
    }
    {
      const ProjectionOptions& po = config.topology.optimizer.projection;
      json::Value pj = json::Value::make_object();
      pj.set("enabled", json::Value::make_bool(po.enabled));
      if (po.enabled) {
        pj.set("eta", json::Value::make_number(po.eta));
        pj.set("beta_start", json::Value::make_number(po.beta_start));
        pj.set("beta_max", json::Value::make_number(po.beta_max));
        pj.set("beta_factor", json::Value::make_number(po.beta_factor));
        pj.set("beta_interval", json::Value::make_number(po.beta_interval));
        pj.set("advance_on_convergence", json::Value::make_bool(po.advance_on_convergence));
        pj.set("formulation",
               json::Value::make_string(
                   "rho_bar = (tanh(beta eta) + tanh(beta (rho_tilde - eta))) / "
                   "(tanh(beta eta) + tanh(beta (1 - eta))) on the filtered density; "
                   "SIMP, volume and stress act on rho_bar"));
        pj.set("robust", json::Value::make_bool(po.robust));
        if (po.robust) {
          pj.set("robust_delta", json::Value::make_number(po.robust_delta));
          pj.set("eta_eroded", json::Value::make_number(po.eroded_eta()));
          pj.set("eta_dilated", json::Value::make_number(po.dilated_eta()));
          pj.set("robust_volume_interval",
                 json::Value::make_number(po.robust_volume_interval));
          pj.set("robust_formulation",
                 json::Value::make_string(
                     "objective and all constraints but the volume on the eroded design "
                     "(eta + delta); volume on the dilated design (eta - delta) against "
                     "V* V_dilated / V_blueprint, rescaled every robust_volume_interval "
                     "iterations; the blueprint (eta) is reported and exported"));
        } else if (po.erosion_check) {
          pj.set("erosion_check", json::Value::make_bool(true));
          pj.set("robust_delta", json::Value::make_number(po.robust_delta));
        }
      }
      setup.set("projection", pj);
    }
    if (result.overhang_filtered) {
      const OverhangOptions& oh = config.topology.optimizer.overhang;
      json::Value am = json::Value::make_object();
      am.set("build_direction", json::Value::make_string(oh.direction.label()));
      am.set("smax_exponent", json::Value::make_number(oh.smax_exponent));
      am.set("smax_reference", json::Value::make_number(oh.smax_reference));
      am.set("smin_epsilon", json::Value::make_number(oh.smin_epsilon));
      am.set("formulation",
             json::Value::make_string(
                 "Langelaar's AM filter on the filtered density, layer by layer from the "
                 "build plate: xi_e = smin(rho_e, smax over the supporting elements of the "
                 "layer below - 3 in 2-D, a cross of 5 in 3-D); passive elements keep "
                 "their density; the projection acts on xi"));
      setup.set("overhang_filter", am);
    }
    if (result.stress_constrained) {
      json::Value sc = json::Value::make_object();
      const StressConstraintOptions& so = config.topology.optimizer.stress;
      sc.set("limit_Pa", json::Value::make_number(so.limit));
      sc.set("p_norm", json::Value::make_number(so.p_norm));
      sc.set("relaxation_q", json::Value::make_number(so.relaxation));
      sc.set("scaling_blend", json::Value::make_number(so.scaling_blend));
      sc.set("feasibility_tolerance", json::Value::make_number(so.feasibility_tolerance));
      sc.set("formulation", json::Value::make_string(
                                "rho^q sigma_vm(solid) at the element centre, p-norm over "
                                "elements with adaptive scale c, one constraint per load "
                                "case: c * g_PN - 1 <= 0"));
      setup.set("stress_constraint", sc);
    }
    if (result.buckling_constrained) {
      json::Value bc = json::Value::make_object();
      const BucklingConstraintOptions& bo = config.topology.optimizer.buckling;
      bc.set("min_load_factor", json::Value::make_number(bo.min_load_factor));
      bc.set("num_modes", json::Value::make_number(bo.num_modes));
      bc.set("ks_parameter", json::Value::make_number(bo.ks_parameter));
      bc.set("ks_overestimate_bound",
             json::Value::make_number(std::log(static_cast<Scalar>(bo.num_modes)) /
                                      bo.ks_parameter));
      bc.set("solid_threshold", json::Value::make_number(bo.solid_threshold));
      bc.set("eigen_tolerance", json::Value::make_number(bo.eigen.tolerance));
      bc.set("eigen_residual_tolerance",
             json::Value::make_number(bo.eigen.residual_tolerance));
      bc.set("formulation",
             json::Value::make_string(
                 "(K(rho) + lambda K_G(rho, u)) phi = 0 with E_K = E_min + rho^p (E_0 - "
                 "E_min) and the stress in K_G from E_G = rho^p E_0 (no E_min floor, "
                 "against pseudo modes in void); one constraint per load case: KS_P("
                 "lambda_req / lambda_i) - 1 <= 0 over the lowest num_modes positive "
                 "load factors, which overestimates max_i lambda_req / lambda_i by at "
                 "most ln(num_modes) / P"));
      setup.set("buckling_constraint", bc);
    }
    out.set("optimization_setup", setup);
  }

  {
    json::Value res = json::Value::make_object();
    res.set("converged", json::Value::make_bool(result.converged));
    res.set("stop_reason", json::Value::make_string(result.stop_reason));
    res.set("final_design_change",
            json::Value::make_number(result.final_design_change));
    res.set("final_relative_objective_change",
            json::Value::make_number(result.final_objective_change));
    res.set("iterations", json::Value::make_number(result.iterations));
    res.set("linear_solves", json::Value::make_number(result.linear_solves));
    {
      json::Value ls = json::Value::make_object();
      ls.set("solver", json::Value::make_string(result.linear_solver));
      ls.set("iterative_iterations_total",
             json::Value::make_number(result.linear_iterations));
      ls.set("iterative_iterations_per_solve",
             json::Value::make_number(result.linear_solves > 0
                                          ? static_cast<Scalar>(result.linear_iterations) /
                                                static_cast<Scalar>(result.linear_solves)
                                          : 0.0));
      if (result.has_amg_stats) ls.set("multigrid", amg_stats_json(result.amg_stats));
      res.set("linear_solver", ls);
    }
    if (result.projected) {
      json::Value pj = json::Value::make_object();
      pj.set("final_beta", json::Value::make_number(result.final_beta));
      pj.set("eta", json::Value::make_number(
                        result.robust ? result.robust_record.eta_intermediate
                                      : result.projection_eta));
      pj.set("filtered_grey_level",
             json::Value::make_number(gray_level(result.filtered_density)));
      pj.set("note", json::Value::make_string(
                         "grey_level is that of the projected (physical) density; "
                         "filtered_grey_level that of the density before projection"));
      res.set("projection", pj);
    }
    if (result.robust || result.erosion_checked) {
      const RobustRecord& rr = result.robust_record;
      json::Value rb = json::Value::make_object();
      json::Value designs = json::Value::make_array();
      const auto design = [&](const char* name, Scalar eta, Scalar c, Scalar vf) {
        json::Value d = json::Value::make_object();
        d.set("design", json::Value::make_string(name));
        d.set("eta", json::Value::make_number(eta));
        d.set("compliance_J", json::Value::make_number(c));
        d.set("volume_fraction", json::Value::make_number(vf));
        designs.push_back(d);
      };
      design("eroded", rr.eta_eroded, rr.compliance_eroded, rr.volume_fraction_eroded);
      design("intermediate (blueprint)", rr.eta_intermediate, rr.compliance_intermediate,
             rr.volume_fraction_intermediate);
      design("dilated", rr.eta_dilated, rr.compliance_dilated, rr.volume_fraction_dilated);
      rb.set("designs", designs);
      if (result.robust) {
        rb.set("dilated_target_fraction",
               json::Value::make_number(rr.dilated_target_fraction));
        rb.set("note", json::Value::make_string(
                           "compliance_J and volume_fraction below are the blueprint's; "
                           "the optimiser minimised the eroded design's compliance"));
        res.set("robust", rb);
      } else {
        rb.set("note", json::Value::make_string(
                           "a design optimised without the robust formulation, evaluated "
                           "once at the eroded and dilated thresholds eta +- robust_delta "
                           "of the final beta: its compliance if the part came out "
                           "uniformly thinner or thicker than drawn"));
        res.set("erosion_check", rb);
      }
    }
    res.set("compliance_J", json::Value::make_number(result.compliance));
    res.set("load_case_compliance_J", json::array_of(result.load_case_compliance));
    res.set("volume_m3", json::Value::make_number(result.volume));
    res.set("volume_fraction", json::Value::make_number(result.volume_fraction));
    res.set("volume_constraint_relative_violation",
            json::Value::make_number(result.volume_constraint_violation));
    res.set("volume_constraint_satisfied",
            json::Value::make_bool(result.volume_constraint_violation <= 1.0e-6));
    res.set("grey_level", json::Value::make_number(result.gray_level));
    res.set("mass_kg",
            json::Value::make_number(result.volume * model.material().density()));
    res.set("total_seconds", json::Value::make_number(result.total_seconds));
    res.set("seconds_per_iteration",
            json::Value::make_number(result.iterations > 0
                                         ? result.total_seconds / result.iterations
                                         : 0.0));
    res.set("warnings", json::array_of(result.warnings));
    if (!result.history.empty()) {
      res.set("initial_compliance_J",
              json::Value::make_number(result.history.front().compliance));
      res.set("compliance_reduction_factor",
              json::Value::make_number(result.compliance > 0.0
                                           ? result.history.front().compliance /
                                                 result.compliance
                                           : 0.0));
    }
    if (result.method == OptimizerMethod::MMA) {
      res.set("constraint_violation", json::Value::make_number(result.constraint_violation));
      res.set("feasible", json::Value::make_bool(result.feasible));
    }
    if (result.stress_constrained) {
      json::Value stress = json::Value::make_object();
      stress.set("max_relaxed_stress_ratio",
                 json::Value::make_number(result.max_stress_ratio));
      stress.set("max_relaxed_stress_Pa",
                 json::Value::make_number(result.max_stress_ratio *
                                          config.topology.optimizer.stress.limit));
      json::Value cases = json::Value::make_array();
      for (const StressConstraintRecord& rec : result.stress) {
        json::Value c = json::Value::make_object();
        c.set("load_case", json::Value::make_string(rec.load_case));
        c.set("max_relaxed_stress_Pa", json::Value::make_number(rec.max_relaxed_stress_Pa));
        c.set("max_relaxed_stress_ratio", json::Value::make_number(rec.max_relaxed_ratio));
        c.set("max_element", json::Value::make_number(rec.max_element));
        c.set("p_norm_ratio", json::Value::make_number(rec.p_norm_ratio));
        c.set("scale", json::Value::make_number(rec.scale));
        c.set("constraint_value", json::Value::make_number(rec.constraint));
        c.set("max_solid_stress_ratio_retained",
              json::Value::make_number(rec.max_solid_ratio_retained));
        cases.push_back(c);
      }
      stress.set("load_cases", cases);
      stress.set("note", json::Value::make_string(
                             "The constraint bounds the relaxed (rho^q) aggregated stress "
                             "of the SIMP model. max_solid_stress_ratio_retained is the "
                             "unrelaxed solid-material stress over the elements at or above "
                             "the interpretation threshold, and "
                             "interpreted_solid_analysis.max_von_mises_over_limit is the "
                             "re-solve of the thresholded structure: those are the numbers "
                             "that say whether the structure meets the limit."));
      res.set("stress", stress);
    }
    if (result.buckling_constrained) {
      json::Value buckling = json::Value::make_object();
      buckling.set("min_load_factor", json::Value::make_number(result.min_load_factor));
      json::Value cases = json::Value::make_array();
      for (const BucklingConstraintRecord& rec : result.buckling) {
        json::Value c = json::Value::make_object();
        c.set("load_case", json::Value::make_string(rec.load_case));
        c.set("load_factors", json::array_of(rec.load_factors));
        c.set("solid_energy_fraction", json::array_of(rec.solid_energy_fraction));
        c.set("ks", json::Value::make_number(rec.ks));
        c.set("constraint_value", json::Value::make_number(rec.constraint));
        c.set("no_positive_load_factor", json::Value::make_bool(rec.no_positive_load_factor));
        cases.push_back(c);
      }
      buckling.set("load_cases", cases);
      buckling.set("note", json::Value::make_string(
                               "Load factors of the SIMP model the constraint acts on "
                               "(void material keeps rho^p of its stress stiffness). The "
                               "buckling_check block re-solves the interpreted structure, "
                               "which is the number that says whether the part is stable."));
      res.set("buckling", buckling);
    }
    out.set("optimization_result", res);
  }

  if (interpretation != nullptr) {
    json::Value interp = json::Value::make_object();
    interp.set("threshold", json::Value::make_number(interpretation->threshold));
    interp.set("elements_above_threshold",
               json::Value::make_number(interpretation->elements_above_threshold));
    interp.set("elements_retained",
               json::Value::make_number(interpretation->elements_retained));
    interp.set("connected_groups_above_threshold",
               json::Value::make_number(interpretation->components_above_threshold));
    interp.set("volume_above_threshold_m3",
               json::Value::make_number(interpretation->volume_above_threshold));
    interp.set("volume_retained_m3",
               json::Value::make_number(interpretation->volume_retained));
    interp.set("volume_discarded_as_islands_m3",
               json::Value::make_number(interpretation->volume_discarded_as_islands));
    interp.set("note", json::Value::make_string(
                           "A SIMP density field is a relative material distribution. "
                           "This block records how it was turned into a solid body for "
                           "the modal comparison: threshold, largest edge-connected "
                           "group, and the material discarded in the process."));
    out.set("solid_interpretation", interp);
  }

  if (modal_initial != nullptr) out.set("modal_initial_solid", modal_json(*modal_initial));
  if (modal_mass_matched != nullptr) {
    out.set("modal_mass_matched_baseline", modal_json(*modal_mass_matched));
  }
  if (modal_optimised != nullptr) {
    out.set("modal_optimised_topology", modal_json(*modal_optimised));
  }

  out.set("timings_s", timings_json(timings));
  return out;
}

}  // namespace sparlab
