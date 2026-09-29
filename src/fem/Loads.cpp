#include "sparlab/fem/Loads.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"

#include <cmath>
#include <sstream>

namespace sparlab {
namespace {

void require_continuum(const FemModel& model, const char* what) {
  if (model.dofs_per_node() != model.dim()) {
    std::ostringstream os;
    os << what << " is formulated for continuum elements here; the "
       << to_string(model.mesh().element_type()) << " element supplies its own";
    throw ConfigError(os.str());
  }
}

void require_in_plane(const Vector3& v, int dim, const std::string& what) {
  if (dim == 2 && v.z() != 0.0) {
    std::ostringstream os;
    os << what << " has a z component of " << v.z()
       << " but the model is two-dimensional and carries no out-of-plane load";
    throw ConfigError(os.str());
  }
}

/// Unit-density consistent mass of element `e`, cached for a uniform grid.
class UnitMass {
 public:
  explicit UnitMass(const FemModel& model) : model_(model) {
    const auto& info = model.mesh().structured_info();
    uniform_ = info.has_value() && info->uniform;
  }
  const Matrix& operator()(Index e) {
    if (uniform_ && cached_.size() > 0) return cached_;
    cached_ = model_.element().consistent_mass(model_.mesh().element_coordinates(e), 1.0,
                                               model_.thickness(), model_.integration());
    return cached_;
  }

 private:
  const FemModel& model_;
  bool uniform_ = false;
  Matrix cached_;
};

}  // namespace

Vector assemble_body_load_vector(const FemModel& model, const LoadCaseSpec& spec) {
  const Mesh& mesh = model.mesh();
  const int dim = mesh.dim();
  const int ndpn = model.dofs_per_node();
  const int npe = mesh.nodes_per_elem();
  const Index ne = mesh.num_elements();
  Vector f = Vector::Zero(model.dofs().num_dofs());
  if (!spec.has_body_loads()) return f;
  require_continuum(model, "a body load");

  require_in_plane(spec.gravity, dim, "gravity of load case '" + spec.name + "'");
  // Per-element body force density per node: b(x_a) [N/m^3], node-major.
  std::vector<Vector> density_field(static_cast<std::size_t>(ne));
  const auto field = [&](Index e) -> Vector& {
    Vector& b = density_field[static_cast<std::size_t>(e)];
    if (b.size() == 0) b = Vector::Zero(npe * ndpn);
    return b;
  };

  if (spec.gravity.squaredNorm() > 0.0) {
    Index massless = 0;
    for (Index e = 0; e < ne; ++e) {
      const Scalar rho = model.material_of(e).density();
      if (rho <= 0.0) {
        ++massless;
        continue;
      }
      Vector& b = field(e);
      for (int a = 0; a < npe; ++a) {
        for (int k = 0; k < dim; ++k) b(ndpn * a + k) += rho * spec.gravity(k);
      }
    }
    if (massless == ne) {
      throw ConfigError("load case '" + spec.name +
                        "' applies gravity, but no material has a density; set "
                        "material.density to give the model a weight");
    }
    if (massless > 0) {
      log::warn("load case '", spec.name, "': gravity leaves ", massless,
                " element(s) of zero-density material weightless");
    }
  }

  for (const BodyForceSpec& body : spec.body_forces) {
    require_in_plane(body.force_density, dim,
                     "body force '" + body.region.name + "' of load case '" + spec.name + "'");
    std::vector<Index> elements;
    if (body.whole_model) {
      elements.resize(static_cast<std::size_t>(ne));
      for (Index e = 0; e < ne; ++e) elements[static_cast<std::size_t>(e)] = e;
    } else {
      elements = body.region.select_elements(mesh);
      if (elements.empty()) {
        throw ConfigError("body force region '" + body.region.name + "' in load case '" +
                          spec.name + "' selected no element");
      }
    }
    for (Index e : elements) {
      Vector& b = field(e);
      for (int a = 0; a < npe; ++a) {
        for (int k = 0; k < dim; ++k) b(ndpn * a + k) += body.force_density(k);
      }
    }
  }

  if (spec.centrifugal.enabled) {
    const CentrifugalSpec& c = spec.centrifugal;
    if (!(c.axis.norm() > 0.0) || !c.axis.allFinite()) {
      throw ConfigError("load case '" + spec.name + "': the centrifugal axis must be a "
                        "non-zero direction");
    }
    const Vector3 e_axis = c.axis.normalized();
    if (dim == 2 && std::abs(std::abs(e_axis.z()) - 1.0) > 1.0e-12) {
      throw ConfigError("load case '" + spec.name +
                        "': a plane model can only rotate about an axis normal to its "
                        "plane (\"axis\": [0, 0, 1]); an in-plane axis would load it "
                        "out of its plane");
    }
    const Matrix3 perp = Matrix3::Identity() - e_axis * e_axis.transpose();
    const Scalar w2 = c.angular_velocity * c.angular_velocity;
    for (Index e = 0; e < ne; ++e) {
      const Scalar rho = model.material_of(e).density();
      if (rho <= 0.0) continue;
      Vector& b = field(e);
      const Index* nodes = mesh.element_nodes(e);
      for (int a = 0; a < npe; ++a) {
        const Vector3 r = perp * (mesh.node(nodes[a]) - c.point);
        for (int k = 0; k < dim; ++k) b(ndpn * a + k) += rho * w2 * r(k);
      }
    }
  }

  UnitMass unit_mass(model);
  for (Index e = 0; e < ne; ++e) {
    const Vector& b = density_field[static_cast<std::size_t>(e)];
    if (b.size() == 0) continue;
    const Vector fe = unit_mass(e) * b;
    model.dofs().scatter_add(mesh.element_nodes(e), npe, fe, f);
  }
  return f;
}

Vector resolve_region_temperatures(const Mesh& mesh, const TemperatureSpec& spec) {
  Vector t = Vector::Constant(mesh.num_nodes(), spec.uniform);
  if (spec.source == TemperatureSpec::Source::Regions) {
    for (const RegionValue& region : spec.regions) {
      const std::vector<Index> nodes = region.region.select_nodes(mesh);
      if (nodes.empty()) {
        throw ConfigError("temperature region '" + region.region.name +
                          "' selected no nodes; check its coordinates");
      }
      for (Index n : nodes) t(n) = region.value;
    }
  }
  return t;
}

Vector linear_internal_parameters(const FemModel& model, Index e, const Vector& ue,
                                  const Vector* temperature) {
  const Element& element = model.element();
  const int ni = element.num_internal_dofs();
  if (ni == 0) return Vector();
  const Mesh& mesh = model.mesh();
  const Matrix coords = mesh.element_coordinates(e);
  const Matrix& d = model.constitutive_of(e);
  const InternalCondensation c =
      element.condense_internal(coords, d, model.thickness(), model.integration());
  Vector alpha = -(c.coupling * ue);
  if (temperature != nullptr && model.material_of(e).thermal_expansion() != 0.0) {
    const Scalar t = mesh.dim() == 2 ? model.thickness() : 1.0;
    Vector fa = Vector::Zero(ni);
    for (const IntegrationPoint& ip : element.integration_rule(model.integration())) {
      const StrainOperator op = element.strain_operator(coords, ip.point);
      const Vector eps0 = element_thermal_strain(model, e, ip.point, *temperature);
      fa.noalias() += (t * ip.weight * op.detJ) *
                      (element.internal_strain_operator(coords, ip.point).transpose() * (d * eps0));
    }
    alpha.noalias() += c.inverse * fa;
  }
  return alpha;
}

Vector linear_point_strain(const FemModel& model, Index e, const NaturalPoint& point,
                           const Vector& ue, const Vector& alpha) {
  const Matrix coords = model.mesh().element_coordinates(e);
  Vector strain = model.element().strain_operator(coords, point).b * ue;
  if (alpha.size() > 0) {
    strain.noalias() += model.element().internal_strain_operator(coords, point) * alpha;
  }
  return strain;
}

Scalar element_temperature_change(const FemModel& model, Index element,
                                  const NaturalPoint& point, const Vector& temperature) {
  const Mesh& mesh = model.mesh();
  const Vector n = model.element().shape_functions(point);
  const Index* nodes = mesh.element_nodes(element);
  Scalar t = 0.0;
  for (int a = 0; a < mesh.nodes_per_elem(); ++a) t += n(a) * temperature(nodes[a]);
  return t - model.material_of(element).reference_temperature();
}

Vector element_thermal_strain(const FemModel& model, Index element, const NaturalPoint& point,
                              const Vector& temperature) {
  const Scalar dt = element_temperature_change(model, element, point, temperature);
  return model.material_of(element).thermal_strain(model.stress_state(), dt);
}

ThermalLoad assemble_thermal_load(const FemModel& model, const Vector& temperature,
                                  const Vector* stiffness_scale) {
  require_continuum(model, "a thermal load");
  const Mesh& mesh = model.mesh();
  if (temperature.size() != mesh.num_nodes()) {
    std::ostringstream os;
    os << "temperature field has " << temperature.size() << " entries for "
       << mesh.num_nodes() << " nodes";
    throw ModelError(os.str());
  }
  const Element& element = model.element();
  const int npe = mesh.nodes_per_elem();
  const Scalar t = mesh.dim() == 2 ? model.thickness() : 1.0;
  const std::vector<IntegrationPoint> rule = element.integration_rule(model.integration());
  ThermalLoad out;
  out.force = Vector::Zero(model.dofs().num_dofs());
  Index expanding = 0;
  for (Index e = 0; e < mesh.num_elements(); ++e) {
    if (model.material_of(e).thermal_expansion() == 0.0) continue;
    ++expanding;
    const Scalar s = stiffness_scale ? (*stiffness_scale)(e) : 1.0;
    const Matrix& d = model.constitutive_of(e);
    const Matrix coords = mesh.element_coordinates(e);
    Vector fe = Vector::Zero(element.num_dofs());
    const int ni = element.num_internal_dofs();
    Vector fa = Vector::Zero(ni);  // int B~^T D eps0 of the internal modes
    for (const IntegrationPoint& ip : rule) {
      const StrainOperator op = element.strain_operator(coords, ip.point);
      const Vector eps0 = element_thermal_strain(model, e, ip.point, temperature);
      const Vector sigma0 = s * (d * eps0);
      const Scalar dv = t * ip.weight * op.detJ;
      fe.noalias() += dv * (op.b.transpose() * sigma0);
      out.self_energy += 0.5 * dv * eps0.dot(sigma0);
      if (ni > 0) fa.noalias() += dv * (element.internal_strain_operator(coords, ip.point).transpose() * sigma0);
    }
    if (ni > 0) {
      // Condensed: f* = int B-hat^T D eps0 = f_u - (K_aa^-1 K_au)^T f_a, and
      // the self energy of the element's own relaxation, less
      // 1/2 f_a^T K_aa^-1 f_a (both with the scaled D, whose factor cancels
      // in the coupling and scales the inverse).
      const InternalCondensation c = element.condense_internal(coords, d, model.thickness(),
                                                               model.integration());
      fe.noalias() -= c.coupling.transpose() * fa;
      if (s > 0.0) out.self_energy -= 0.5 * fa.dot(c.inverse * fa) / s;
    }
    model.dofs().scatter_add(mesh.element_nodes(e), npe, fe, out.force);
  }
  if (expanding == 0) {
    log::warn("a temperature field is applied, but no material has a thermal expansion "
              "coefficient; set material.thermal_expansion for a thermal strain");
  }
  return out;
}

}  // namespace sparlab
