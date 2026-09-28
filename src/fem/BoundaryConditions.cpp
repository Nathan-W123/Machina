#include "sparlab/fem/BoundaryConditions.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"
#include "sparlab/elements/FaceGeometry.hpp"

#include <algorithm>
#include <set>
#include <sstream>

namespace sparlab {
namespace {

void require_in_plane(const Vector3& v, int dim, const std::string& what) {
  if (dim == 2 && v.z() != 0.0) {
    std::ostringstream os;
    os << what << " has a z component of " << v.z()
       << " but the model is two-dimensional; a plane model carries no out-of-plane "
          "component (use a 3-D mesh or drop it)";
    throw ConfigError(os.str());
  }
}

}  // namespace

std::vector<Mesh::BoundaryFace> faces_in_region(const Mesh& mesh,
                                                const std::vector<Mesh::BoundaryFace>& boundary,
                                                const SelectorGroup& region) {
  std::vector<char> in_region(static_cast<std::size_t>(mesh.num_nodes()), 0);
  for (Index n : region.select_nodes(mesh)) in_region[static_cast<std::size_t>(n)] = 1;
  std::vector<Mesh::BoundaryFace> out;
  for (const Mesh::BoundaryFace& face : boundary) {
    const bool inside = std::all_of(face.nodes.begin(), face.nodes.end(), [&](Index n) {
      return in_region[static_cast<std::size_t>(n)] != 0;
    });
    if (inside) out.push_back(face);
  }
  return out;
}

void DisplacementConstraint::set(int component, bool fix, Scalar value) {
  switch (component) {
    case 0: fix_x = fix; value_x = value; return;
    case 1: fix_y = fix; value_y = value; return;
    case 2: fix_z = fix; value_z = value; return;
    case 3: fix_rx = fix; value_rx = value; return;
    case 4: fix_ry = fix; value_ry = value; return;
    case 5: fix_rz = fix; value_rz = value; return;
    default: break;
  }
  std::ostringstream os;
  os << "DOF component " << component << " is outside [0, " << kMaxDofsPerNode - 1 << "]";
  throw ModelError(os.str());
}

Index apply_constraints(const Mesh& mesh,
                        const std::vector<DisplacementConstraint>& constraints,
                        DofManager& dofs) {
  const int dim = mesh.dim();
  const int ndpn = dofs.dofs_per_node();
  std::set<Index> touched;
  for (const DisplacementConstraint& bc : constraints) {
    bool any = false;
    for (int k = 0; k < kMaxDofsPerNode; ++k) any = any || bc.fixes(k);
    if (!any) {
      throw ConfigError("boundary condition '" + bc.region.name +
                        "' constrains no component; remove it or list the "
                        "components to fix");
    }
    if (dim == 2 && bc.fix_z) {
      throw ConfigError("boundary condition '" + bc.region.name +
                        "' fixes z, but the model is two-dimensional and has no z "
                        "degree of freedom");
    }
    if (ndpn < kMaxDofsPerNode && bc.fixes_rotation()) {
      throw ConfigError("boundary condition '" + bc.region.name +
                        "' fixes a rotation, but the model's nodes carry translations "
                        "only; rotations are degrees of freedom of shell and beam "
                        "elements");
    }
    const std::vector<Index> nodes = bc.region.select_nodes(mesh);
    if (nodes.empty()) {
      throw ConfigError("boundary condition '" + bc.region.name +
                        "' selected no nodes; check its coordinates against the mesh "
                        "extents (an unconstrained model yields a singular stiffness "
                        "matrix)");
    }
    for (Index n : nodes) {
      for (int k = 0; k < ndpn; ++k) {
        if (!bc.fixes(k)) continue;
        dofs.prescribe(n, k, bc.value(k));
        touched.insert(dofs.dof(n, k));
      }
    }
    if (static_cast<int>(log::level()) <= static_cast<int>(log::Level::Debug)) {
      std::ostringstream fixed;
      for (int k = 0; k < ndpn; ++k) {
        if (bc.fixes(k)) fixed << ' ' << dof_component_name(k);
      }
      log::debug("boundary condition '", bc.region.name, "': ", nodes.size(),
                 " nodes, fixed:", fixed.str());
    }
  }
  return static_cast<Index>(touched.size());
}

Vector assemble_load_vector(const Mesh& mesh, const Element& element,
                            const LoadCaseSpec& load_case, Scalar thickness,
                            const IntegrationOptions& integration) {
  const int dim = mesh.dim();
  const int ndpn = element.dofs_per_node();
  Vector f = Vector::Zero(mesh.num_nodes() * ndpn);

  for (const PointLoadSpec& load : load_case.point_loads) {
    require_in_plane(load.force, dim,
                     "point load '" + load.region.name + "' in load case '" +
                         load_case.name + "'");
    if (ndpn < kMaxDofsPerNode && load.moment.squaredNorm() != 0.0) {
      throw ConfigError("point load '" + load.region.name + "' in load case '" +
                        load_case.name +
                        "' applies a moment, but the model's nodes carry translations "
                        "only; a moment needs shell or beam elements (on a continuum "
                        "mesh, apply it as a couple of forces or a traction)");
    }
    const std::vector<Index> nodes = load.region.select_nodes(mesh);
    if (nodes.empty()) {
      throw ConfigError("point load region '" + load.region.name +
                        "' in load case '" + load_case.name +
                        "' selected no nodes; check its coordinates");
    }
    const Scalar scale =
        load.distribute_total ? 1.0 / static_cast<Scalar>(nodes.size()) : 1.0;
    for (Index n : nodes) {
      for (int k = 0; k < dim; ++k) f(n * ndpn + k) += scale * load.force(k);
      if (ndpn == kMaxDofsPerNode) {
        for (int k = 0; k < 3; ++k) f(n * ndpn + 3 + k) += scale * load.moment(k);
      }
    }
    log::debug("point load '", load.region.name, "' in case '", load_case.name, "': ",
               nodes.size(), " nodes, resultant (", load.force.x(), ", ",
               load.force.y(), (dim == 3 ? ", " : ""), (dim == 3 ? load.force.z() : 0.0),
               ") N", load.distribute_total ? " distributed" : " per node");
  }

  if (!load_case.tractions.empty()) {
    const std::vector<Mesh::BoundaryFace> faces = mesh.boundary_faces();
    for (const TractionLoadSpec& load : load_case.tractions) {
      require_in_plane(load.traction, dim,
                       "traction '" + load.region.name + "' in load case '" +
                           load_case.name + "'");
      Index matched = 0;
      Scalar total_measure = 0.0;
      for (const Mesh::BoundaryFace& face : faces_in_region(mesh, faces, load.region)) {
        const Matrix coords = mesh.element_coordinates(face.element);
        const Vector fe = element.boundary_traction(coords, face.local_face,
                                                    load.traction, thickness, integration);
        const Index* enodes = mesh.element_nodes(face.element);
        for (int a = 0; a < mesh.nodes_per_elem(); ++a) {
          for (int k = 0; k < ndpn; ++k) f(enodes[a] * ndpn + k) += fe(ndpn * a + k);
        }
        total_measure += mesh.face_measure(face);
        ++matched;
      }
      if (matched == 0) {
        throw ConfigError("traction region '" + load.region.name + "' in load case '" +
                          load_case.name + "' matched no boundary " +
                          (dim == 2 ? "edge" : "face") +
                          "; the region must contain every node of at least one " +
                          (dim == 2 ? "edge" : "face") + " on the mesh boundary");
      }
      log::debug("traction '", load.region.name, "' in case '", load_case.name, "': ",
                 matched, (dim == 2 ? " boundary edges, length " : " boundary faces, area "),
                 total_measure, (dim == 2 ? " m, traction (" : " m^2, traction ("),
                 load.traction.x(), ", ", load.traction.y(),
                 (dim == 3 ? ", " : ""), (dim == 3 ? load.traction.z() : 0.0), ") Pa");
    }
  }

  if (!load_case.pressures.empty()) {
    if (ndpn != dim) {
      throw ConfigError("load case '" + load_case.name +
                        "': a boundary pressure acts on the faces of a continuum mesh; "
                        "shell and beam models take their own surface loads");
    }
    const std::vector<Mesh::BoundaryFace> faces = mesh.boundary_faces();
    const FaceShape shape = face_shape_of(mesh.element_type());
    // A curved six-node face needs three points per direction for N_a times
    // its quadratic area vector.
    const int points = shape == FaceShape::Tri6 ? std::max(integration.edge_points, 3)
                                                : integration.edge_points;
    const std::vector<std::vector<int>>& table = element_local_faces(mesh.element_type());
    for (const PressureLoadSpec& load : load_case.pressures) {
      Index matched = 0;
      Scalar total_measure = 0.0;
      for (const Mesh::BoundaryFace& face : faces_in_region(mesh, faces, load.region)) {
        const Matrix xf = element_face_coordinates(mesh, face.element, face.local_face);
        const Matrix fe = face_pressure_forces(shape, xf, load.pressure, thickness, points);
        const Index* enodes = mesh.element_nodes(face.element);
        const std::vector<int>& fn = table[static_cast<std::size_t>(face.local_face)];
        for (std::size_t a = 0; a < fn.size(); ++a) {
          const Index node = enodes[fn[a]];
          for (int k = 0; k < dim; ++k) f(node * ndpn + k) += fe(k, static_cast<Eigen::Index>(a));
        }
        total_measure += mesh.face_measure(face);
        ++matched;
      }
      if (matched == 0) {
        throw ConfigError("pressure region '" + load.region.name + "' in load case '" +
                          load_case.name + "' matched no boundary " +
                          (dim == 2 ? "edge" : "face") +
                          "; the region must contain every node of at least one " +
                          (dim == 2 ? "edge" : "face") + " on the mesh boundary");
      }
      log::debug("pressure '", load.region.name, "' in case '", load_case.name, "': ",
                 matched, (dim == 2 ? " boundary edges, length " : " boundary faces, area "),
                 total_measure, ", pressure ", load.pressure, " Pa");
    }
  }

  return f;
}

}  // namespace sparlab
