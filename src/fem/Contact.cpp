#include "sparlab/fem/Contact.hpp"

#include "NonlinearSystem.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"
#include "sparlab/elements/FaceGeometry.hpp"
#include "sparlab/elements/Quadrature.hpp"
#include "sparlab/fem/BoundaryConditions.hpp"
#include "sparlab/fem/FemModel.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <sstream>

namespace sparlab {
namespace {

/// A slave node's coverage by the master surface, sum_l M_jl / D_j, below
/// which it is left out of contact: a node the master surface does not
/// cover (or covers only in part, at the master's edge) has no consistent
/// weighted gap - a rigid translation of both bodies would change it.
constexpr Scalar kMinCoverage = 0.99;
/// Above this the master surface overlaps itself in projection (both sides
/// of a thin part selected, a corner folded onto the slave).
constexpr Scalar kMaxCoverage = 1.01;

/// A face of a contact surface in the reference configuration.
struct SurfaceFace {
  Index element = 0;
  std::vector<Index> nodes;  ///< global, in the face table's order
  Matrix coords;             ///< dim x nf
  Matrix dual;               ///< psi_a = sum_b dual(a, b) N_b
  Vector weights;            ///< int N_a dA per face node (times the thickness in 2-D)
  Vector3 lo = Vector3::Zero();  ///< bounding box
  Vector3 hi = Vector3::Zero();
  Scalar size = 0.0;         ///< diagonal of the bounding box [m]
};

Vector3 padded(const Matrix& x, Eigen::Index col) {
  Vector3 p = Vector3::Zero();
  p.head(x.rows()) = x.col(col);
  return p;
}

/// Reference coordinates of the nodes of a face shape.
std::vector<std::pair<Scalar, Scalar>> face_node_coordinates(FaceShape shape) {
  switch (shape) {
    case FaceShape::Line2: return {{-1.0, 0.0}, {1.0, 0.0}};
    case FaceShape::Tri3: return {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}};
    case FaceShape::Quad4: return {{-1.0, -1.0}, {1.0, -1.0}, {1.0, 1.0}, {-1.0, 1.0}};
    case FaceShape::Tri6: break;
  }
  throw ConfigError("contact needs faces of linear elements");
}

/// The unit outward normal of a face at reference point (s, t).
Vector3 face_normal(FaceShape shape, const Matrix& coords, Scalar s, Scalar t) {
  Vector n;
  Matrix dn;
  face_shape_functions(shape, s, t, n, dn);
  const FaceGeometryPoint g = face_geometry(shape, coords, s, t, n, dn);
  if (!(g.jacobian > 0.0)) throw ModelError("a contact face is degenerate (zero area)");
  return g.area / g.jacobian;
}

std::vector<SurfaceFace> build_faces(const FemModel& model,
                                     const std::vector<Mesh::BoundaryFace>& selected) {
  const Mesh& mesh = model.mesh();
  const FaceShape shape = face_shape_of(mesh.element_type());
  const Scalar thickness = mesh.dim() == 2 ? model.thickness() : 1.0;
  const std::vector<std::vector<int>>& table = element_local_faces(mesh.element_type());
  std::vector<SurfaceFace> faces;
  faces.reserve(selected.size());
  for (const Mesh::BoundaryFace& bf : selected) {
    SurfaceFace f;
    f.element = bf.element;
    const Index* en = mesh.element_nodes(bf.element);
    for (int l : table[static_cast<std::size_t>(bf.local_face)]) f.nodes.push_back(en[l]);
    f.coords = element_face_coordinates(mesh, bf.element, bf.local_face);
    f.weights = face_shape_integrals(shape, f.coords, thickness, 3);
    const Matrix products = face_shape_products(shape, f.coords, thickness, 3);
    // Dual basis: psi = A N with A = D M^-1, so that int psi_a N_b = D_a delta_ab.
    f.dual = f.weights.asDiagonal() * products.inverse();
    f.lo = f.hi = padded(f.coords, 0);
    for (Eigen::Index a = 1; a < f.coords.cols(); ++a) {
      const Vector3 p = padded(f.coords, a);
      f.lo = f.lo.cwiseMin(p);
      f.hi = f.hi.cwiseMax(p);
    }
    f.size = (f.hi - f.lo).norm();
    faces.push_back(std::move(f));
  }
  return faces;
}

Scalar cross2(const Vector3& a, const Vector3& b) { return a.x() * b.y() - a.y() * b.x(); }

/// Accumulated mortar integrals of one slave node.
struct MortarRow {
  std::map<Index, Scalar> m;  ///< master node -> M_jl
  Scalar gap = 0.0;           ///< int psi_j g0
};

/// 2-D segment integration of one slave edge against one master edge
/// (Popp, Gee and Wall 2009): the master nodes are projected onto the slave
/// edge along the slave's interpolated normal field, which fixes the
/// overlap; at each Gauss point of the overlap the slave point is projected
/// along that normal onto the master edge.
void integrate_segment_2d(const SurfaceFace& slave, const std::vector<Vector3>& slave_normals,
                          const SurfaceFace& master, Scalar thickness,
                          std::map<Index, MortarRow>& rows) {
  const Vector3 x0 = padded(slave.coords, 0);
  const Vector3 x1 = padded(slave.coords, 1);
  const Vector3 y0 = padded(master.coords, 0);
  const Vector3 y1 = padded(master.coords, 1);
  const Vector3& n0 = slave_normals[0];
  const Vector3& n1 = slave_normals[1];
  // The master edge must face the slave one (its outward normal against the
  // slave's): the far side of a thin master part is not paired.
  const Vector3 ydir = y1 - y0;
  if (Vector3(ydir.y(), -ydir.x(), 0.0).dot(n0 + n1) >= 0.0) return;
  const Vector3 xm = 0.5 * (x0 + x1);
  const Vector3 hs = 0.5 * (x1 - x0);
  const Vector3 nm = 0.5 * (n0 + n1);
  const Vector3 dn = 0.5 * (n1 - n0);
  // (X(xi) - Y) x n(xi) = 0 is quadratic in xi.
  const auto to_slave = [&](const Vector3& y, Scalar& xi) -> bool {
    const Vector3 r = xm - y;
    const Scalar a = cross2(hs, dn);
    const Scalar b = cross2(r, dn) + cross2(hs, nm);
    const Scalar c = cross2(r, nm);
    const Scalar scale = hs.norm() * (nm.norm() + dn.norm());
    if (std::abs(a) <= 1.0e-12 * scale) {
      if (std::abs(b) <= 1.0e-14 * scale) return false;
      xi = -c / b;
      return true;
    }
    const Scalar disc = b * b - 4.0 * a * c;
    if (disc < 0.0) return false;
    const Scalar root = std::sqrt(disc);
    // The numerically stable pair of roots; the one nearer the edge.
    const Scalar q = -0.5 * (b + (b >= 0.0 ? root : -root));
    const Scalar r1 = q / a;
    const Scalar r2 = q != 0.0 ? c / q : r1;
    xi = std::abs(r1) <= std::abs(r2) ? r1 : r2;
    return true;
  };
  Scalar xa = 0.0;
  Scalar xb = 0.0;
  if (!to_slave(y0, xa) || !to_slave(y1, xb)) return;
  const Scalar lo = std::max(-1.0, std::min(xa, xb));
  const Scalar hi = std::min(1.0, std::max(xa, xb));
  if (hi - lo <= 1.0e-12) return;
  const Vector3 ym = 0.5 * (y0 + y1);
  const Vector3 hm = 0.5 * (y1 - y0);
  const Scalar jac = hs.norm();  // |dX/dxi|
  for (const QuadraturePoint1D& g : gauss_legendre_line(4)) {
    const Scalar xi = 0.5 * (lo + hi) + 0.5 * (hi - lo) * g.xi;
    const Scalar w = g.weight * 0.5 * (hi - lo) * jac * thickness;
    const Scalar na = 0.5 * (1.0 - xi);
    const Scalar nb = 0.5 * (1.0 + xi);
    const Vector3 x = na * x0 + nb * x1;
    const Vector3 n = (na * n0 + nb * n1).normalized();
    // (Y(eta) - x) x n = 0, linear in eta.
    const Scalar den = cross2(hm, n);
    if (std::abs(den) <= 1.0e-14 * hm.norm()) continue;
    const Scalar eta = -cross2(ym - x, n) / den;
    if (eta < -1.0 - 1.0e-8 || eta > 1.0 + 1.0e-8) continue;
    const Scalar mc = 0.5 * (1.0 - eta);
    const Scalar md = 0.5 * (1.0 + eta);
    const Vector3 y = mc * y0 + md * y1;
    const Scalar g0 = n.dot(y - x);
    const Eigen::Vector2d shape(na, nb);
    const Eigen::Vector2d psi = slave.dual * shape;
    for (int a = 0; a < 2; ++a) {
      MortarRow& row = rows[slave.nodes[static_cast<std::size_t>(a)]];
      row.m[master.nodes[0]] += w * psi(a) * mc;
      row.m[master.nodes[1]] += w * psi(a) * md;
      row.gap += w * psi(a) * g0;
    }
  }
}

using Polygon = std::vector<Eigen::Vector2d>;

Scalar signed_area(const Polygon& p) {
  Scalar a = 0.0;
  for (std::size_t i = 0; i < p.size(); ++i) {
    const Eigen::Vector2d& u = p[i];
    const Eigen::Vector2d& v = p[(i + 1) % p.size()];
    a += u.x() * v.y() - u.y() * v.x();
  }
  return 0.5 * a;
}

/// Sutherland-Hodgman: `subject` clipped by the convex, counter-clockwise
/// polygon `clip`. Meshes that share lines put vertices on the other
/// polygon's edges and edges along them; a point within round-off of an
/// edge's line counts as inside, and an intersection parameter is kept on
/// its segment, so that a nearly parallel edge cannot throw a point far
/// off (the hull taken afterwards orders what remains).
Polygon clip_polygon(const Polygon& subject, const Polygon& clip) {
  Scalar size = 0.0;
  for (const Eigen::Vector2d& p : clip) size = std::max(size, p.cwiseAbs().maxCoeff());
  for (const Eigen::Vector2d& p : subject) size = std::max(size, p.cwiseAbs().maxCoeff());
  Polygon out = subject;
  for (std::size_t i = 0; i < clip.size() && !out.empty(); ++i) {
    const Eigen::Vector2d a = clip[i];
    const Eigen::Vector2d b = clip[(i + 1) % clip.size()];
    const Eigen::Vector2d e = b - a;
    const Scalar tol = 1.0e-12 * e.norm() * size;
    const auto inside = [&](const Eigen::Vector2d& p) {
      return e.x() * (p.y() - a.y()) - e.y() * (p.x() - a.x()) >= -tol;
    };
    Polygon in = std::move(out);
    out.clear();
    for (std::size_t k = 0; k < in.size(); ++k) {
      const Eigen::Vector2d& p = in[k];
      const Eigen::Vector2d& q = in[(k + 1) % in.size()];
      const bool pin = inside(p);
      const bool qin = inside(q);
      if (pin) out.push_back(p);
      if (pin != qin) {
        const Eigen::Vector2d d = q - p;
        const Scalar den = e.x() * d.y() - e.y() * d.x();
        if (std::abs(den) > 1.0e-14 * e.norm() * d.norm()) {
          const Scalar s =
              std::clamp((e.x() * (a.y() - p.y()) - e.y() * (a.x() - p.x())) / den, 0.0, 1.0);
          out.push_back(p + s * d);
        }
      }
    }
  }
  return out;
}

/// Convex hull (Andrew's monotone chain), counter-clockwise, without
/// repeated or collinear points. The intersection of two convex polygons is
/// convex, but where a vertex of one lies on an edge of the other - as on
/// meshes that share lines - round-off in the clipping's inside test can
/// emit its points out of order; their hull is the intersection.
Polygon convex_hull(Polygon points) {
  if (points.size() < 3) return points;
  std::sort(points.begin(), points.end(), [](const Eigen::Vector2d& a, const Eigen::Vector2d& b) {
    return a.x() < b.x() || (a.x() == b.x() && a.y() < b.y());
  });
  const auto cross = [](const Eigen::Vector2d& o, const Eigen::Vector2d& a,
                        const Eigen::Vector2d& b) {
    return (a.x() - o.x()) * (b.y() - o.y()) - (a.y() - o.y()) * (b.x() - o.x());
  };
  Polygon hull(2 * points.size());
  std::size_t k = 0;
  for (std::size_t i = 0; i < points.size(); ++i) {
    while (k >= 2 && cross(hull[k - 2], hull[k - 1], points[i]) <= 0.0) --k;
    hull[k++] = points[i];
  }
  for (std::size_t i = points.size() - 1, t = k + 1; i > 0; --i) {
    while (k >= t && cross(hull[k - 2], hull[k - 1], points[i - 1]) <= 0.0) --k;
    hull[k++] = points[i - 1];
  }
  hull.resize(k > 0 ? k - 1 : 0);
  return hull;
}

/// The reference coordinates (s, t) of the point of a face that projects,
/// along the auxiliary plane's normal, onto `x` (a point of the plane with
/// in-plane basis e1, e2): Newton on e_i . (X(s, t) - x) = 0.
bool project_to_face(FaceShape shape, const Matrix& coords, const Vector3& x, const Vector3& e1,
                     const Vector3& e2, Scalar& s, Scalar& t) {
  s = shape == FaceShape::Quad4 ? 0.0 : 1.0 / 3.0;
  t = s;
  Vector n;
  Matrix dn;
  for (int it = 0; it < 30; ++it) {
    face_shape_functions(shape, s, t, n, dn);
    Vector3 p = Vector3::Zero();
    Vector3 ps = Vector3::Zero();
    Vector3 pt = Vector3::Zero();
    for (Eigen::Index a = 0; a < coords.cols(); ++a) {
      const Vector3 xa = padded(coords, a);
      p += n(a) * xa;
      ps += dn(a, 0) * xa;
      pt += dn(a, 1) * xa;
    }
    const Eigen::Vector2d f(e1.dot(p - x), e2.dot(p - x));
    Eigen::Matrix2d j;
    j << e1.dot(ps), e1.dot(pt), e2.dot(ps), e2.dot(pt);
    const Scalar det = j.determinant();
    if (std::abs(det) <= 1.0e-300) return false;
    const Eigen::Vector2d step = j.inverse() * f;
    s -= step.x();
    t -= step.y();
    if (step.norm() <= 1.0e-14) return true;
  }
  return false;
}

/// 3-D segment integration of one slave face against one master face on
/// the slave face's auxiliary plane (Puso and Laursen 2004; Popp et al.
/// 2010): both faces are projected onto the plane through the slave face's
/// centre normal to its interpolated nodal normal, the master polygon is
/// clipped by the slave one, and the clip polygon is integrated in
/// triangles, every quadrature point projected back onto both faces.
void integrate_segment_3d(FaceShape shape, const SurfaceFace& slave,
                          const std::vector<Vector3>& slave_normals, const SurfaceFace& master,
                          std::map<Index, MortarRow>& rows) {
  const int nf = static_cast<int>(slave.nodes.size());
  const Scalar sc = shape == FaceShape::Quad4 ? 0.0 : 1.0 / 3.0;
  Vector n;
  Matrix dn;
  face_shape_functions(shape, sc, sc, n, dn);
  Vector3 x0 = Vector3::Zero();
  Vector3 n0 = Vector3::Zero();
  for (int a = 0; a < nf; ++a) {
    x0 += n(a) * padded(slave.coords, a);
    n0 += n(a) * slave_normals[static_cast<std::size_t>(a)];
  }
  n0.normalize();
  // The master face must face the slave one.
  const Scalar mc = master.nodes.size() == 4 ? 0.0 : 1.0 / 3.0;
  if (face_normal(shape, master.coords, mc, mc).dot(n0) >= 0.0) return;
  Vector3 e1 = padded(slave.coords, 1) - padded(slave.coords, 0);
  e1 -= e1.dot(n0) * n0;
  e1.normalize();
  const Vector3 e2 = n0.cross(e1);
  const auto plane = [&](const Vector3& p) {
    return Eigen::Vector2d(e1.dot(p - x0), e2.dot(p - x0));
  };
  Polygon sp;
  for (int a = 0; a < nf; ++a) sp.push_back(plane(padded(slave.coords, a)));
  if (signed_area(sp) < 0.0) std::reverse(sp.begin(), sp.end());
  Polygon mp;
  for (Eigen::Index a = 0; a < master.coords.cols(); ++a) mp.push_back(plane(padded(master.coords, a)));
  if (signed_area(mp) < 0.0) std::reverse(mp.begin(), mp.end());
  const Polygon cp = convex_hull(clip_polygon(mp, sp));
  if (cp.size() < 3) return;
  const Scalar area = signed_area(cp);
  if (!(area > 1.0e-12 * std::abs(signed_area(sp)))) return;
  Eigen::Vector2d centre = Eigen::Vector2d::Zero();
  for (const Eigen::Vector2d& p : cp) centre += p;
  centre /= static_cast<Scalar>(cp.size());
  Vector sn;
  Matrix sdn;
  for (std::size_t k = 0; k < cp.size(); ++k) {
    const Eigen::Vector2d a = cp[k] - centre;
    const Eigen::Vector2d b = cp[(k + 1) % cp.size()] - centre;
    const Scalar tri = 0.5 * (a.x() * b.y() - a.y() * b.x());
    if (!(tri > 0.0)) continue;
    for (const QuadraturePoint2D& q : collapsed_gauss_triangle(4)) {
      const Eigen::Vector2d p = centre + q.xi * a + q.eta * b;
      const Scalar w = q.weight * 2.0 * tri;
      const Vector3 x = x0 + p.x() * e1 + p.y() * e2;
      Scalar s = 0.0;
      Scalar t = 0.0;
      Scalar ms = 0.0;
      Scalar mt = 0.0;
      if (!project_to_face(shape, slave.coords, x, e1, e2, s, t) ||
          !project_to_face(shape, master.coords, x, e1, e2, ms, mt)) {
        continue;
      }
      face_shape_functions(shape, s, t, sn, sdn);
      Vector mn;
      Matrix mdn;
      face_shape_functions(shape, ms, mt, mn, mdn);
      Vector3 xs = Vector3::Zero();
      Vector3 ym = Vector3::Zero();
      for (int a2 = 0; a2 < nf; ++a2) xs += sn(a2) * padded(slave.coords, a2);
      for (Eigen::Index a2 = 0; a2 < master.coords.cols(); ++a2) {
        ym += mn(a2) * padded(master.coords, a2);
      }
      const Scalar g0 = n0.dot(ym - xs);
      const Vector psi = slave.dual * sn;
      for (int a2 = 0; a2 < nf; ++a2) {
        MortarRow& row = rows[slave.nodes[static_cast<std::size_t>(a2)]];
        for (Eigen::Index l = 0; l < mn.size(); ++l) {
          row.m[master.nodes[static_cast<std::size_t>(l)]] += w * psi(a2) * mn(l);
        }
        row.gap += w * psi(a2) * g0;
      }
    }
  }
}

}  // namespace

// ---------------------------------------------------------------------------
// Obstacles and names
// ---------------------------------------------------------------------------
Scalar RigidObstacle::gap(const Vector3& x, Vector3& normal, int dim) const {
  switch (kind) {
    case Kind::Plane: {
      Vector3 n = direction;
      if (dim == 2) n.z() = 0.0;
      const Scalar len = n.norm();
      if (!(len > 0.0)) throw ConfigError("a contact plane needs a non-zero normal");
      normal = n / len;
      return normal.dot(x - point);
    }
    case Kind::Cylinder:
    case Kind::Sphere: {
      Vector3 r = x - point;
      if (kind == Kind::Cylinder) {
        Vector3 axis = dim == 2 ? Vector3::UnitZ() : direction;
        const Scalar len = axis.norm();
        if (!(len > 0.0)) throw ConfigError("a contact cylinder needs a non-zero axis");
        axis /= len;
        r -= r.dot(axis) * axis;
      }
      if (dim == 2) r.z() = 0.0;
      const Scalar d = r.norm();
      if (!(d > 1.0e-12 * std::max(radius, 1.0e-300))) {
        throw ModelError("a contact node lies on the axis (or at the centre) of a rigid "
                         "obstacle, where its normal is undefined");
      }
      const Vector3 radial = r / d;
      if (inside) {
        normal = -radial;
        return radius - d;
      }
      normal = radial;
      return d - radius;
    }
  }
  return 0.0;
}

std::string to_string(RigidObstacle::Kind kind) {
  switch (kind) {
    case RigidObstacle::Kind::Plane: return "plane";
    case RigidObstacle::Kind::Cylinder: return "cylinder";
    case RigidObstacle::Kind::Sphere: return "sphere";
  }
  return "plane";
}

RigidObstacle::Kind parse_obstacle_kind(const std::string& text) {
  if (text == "plane") return RigidObstacle::Kind::Plane;
  if (text == "cylinder") return RigidObstacle::Kind::Cylinder;
  if (text == "sphere") return RigidObstacle::Kind::Sphere;
  throw ConfigError("unknown obstacle type '" + text +
                    "'; expected \"plane\", \"cylinder\" or \"sphere\"");
}

std::string to_string(ContactStatus status) {
  switch (status) {
    case ContactStatus::Open: return "open";
    case ContactStatus::Stick: return "stick";
    case ContactStatus::Slip: return "slip";
  }
  return "open";
}

// ---------------------------------------------------------------------------
// Construction: the geometry of every pair in the reference configuration
// ---------------------------------------------------------------------------
ContactProblem::ContactProblem(const FemModel& model, const ContactOptions& options)
    : model_(model), options_(options) {
  const Mesh& mesh = model.mesh();
  const int dim = mesh.dim();
  if (mesh.element_type() == ElementType::Tet10) {
    throw ConfigError("contact needs linear elements (Q4, Tri3, Hex8, Tet4): the dual basis "
                      "of the contact pressure does not exist on the six-node face of a Tet10, "
                      "whose corner weights int N dA vanish");
  }
  if (options.pairs.empty()) throw ConfigError("'contact' is enabled but defines no pair");
  if (!(options.complementarity > 0.0) || !(options.search_factor > 0.0)) {
    throw ConfigError("contact: 'complementarity' and 'search_factor' must be positive");
  }
  const FaceShape shape = face_shape_of(mesh.element_type());
  const Scalar thickness = dim == 2 ? model.thickness() : 1.0;
  const std::vector<Mesh::BoundaryFace> boundary = mesh.boundary_faces();
  const DofManager& dofs = model.dofs();
  std::map<Index, std::string> slave_owner;
  std::set<Index> all_masters;
  std::vector<std::set<Index>> pair_slaves(options.pairs.size());

  for (std::size_t k = 0; k < options.pairs.size(); ++k) {
    const ContactPairSpec& pair = options.pairs[k];
    const std::string label = "contact pair '" + pair.name + "'";
    if (!(pair.friction >= 0.0) || !std::isfinite(pair.friction)) {
      throw ConfigError(label + ": the friction coefficient must be >= 0");
    }
    if (pair.friction > 0.0) frictional_ = true;
    if (pair.rigid && pair.obstacle.kind != RigidObstacle::Kind::Plane &&
        !(pair.obstacle.radius > 0.0)) {
      throw ConfigError(label + ": the obstacle needs a positive radius");
    }
    const std::vector<SurfaceFace> slave =
        build_faces(model, faces_in_region(mesh, boundary, pair.slave));
    if (slave.empty()) {
      throw ConfigError(label + ": the slave region selects no boundary face (every node of "
                        "a face must lie in it)");
    }
    std::vector<SurfaceFace> master;
    if (!pair.rigid) {
      master = build_faces(model, faces_in_region(mesh, boundary, pair.master));
      if (master.empty()) {
        throw ConfigError(label + ": the master region selects no boundary face");
      }
    }

    // Slave nodes: weights, averaged outward normals, the stiffest adjacent
    // material.
    std::map<Index, Scalar> weight;
    std::map<Index, Vector3> normal_sum;
    std::map<Index, Scalar> modulus;
    const std::vector<std::pair<Scalar, Scalar>> corners = face_node_coordinates(shape);
    std::vector<std::vector<Vector3>> face_normals(slave.size());
    for (std::size_t f = 0; f < slave.size(); ++f) {
      const SurfaceFace& face = slave[f];
      const Scalar e = model.material_of(face.element).youngs_modulus();
      for (std::size_t a = 0; a < face.nodes.size(); ++a) {
        const Index node = face.nodes[a];
        weight[node] += face.weights(static_cast<Eigen::Index>(a));
        // (emplace: a default-constructed Eigen vector is uninitialised)
        normal_sum.emplace(node, Vector3::Zero()).first->second +=
            face_normal(shape, face.coords, corners[a].first, corners[a].second);
        modulus[node] = std::max(modulus[node], e);
      }
    }
    std::map<Index, Vector3> slave_normal;
    for (const auto& [node, sum] : normal_sum) {
      const Scalar len = sum.norm();
      if (!(len > 1.0e-12)) {
        throw ConfigError(label + ": the slave surface folds back on itself at node " +
                          std::to_string(node) + " (its face normals cancel)");
      }
      slave_normal[node] = sum / len;
    }
    for (std::size_t f = 0; f < slave.size(); ++f) {
      for (Index node : slave[f].nodes) face_normals[f].push_back(slave_normal[node]);
    }
    for (const auto& [node, w] : weight) {
      (void)w;
      const auto it = slave_owner.find(node);
      if (it != slave_owner.end()) {
        throw ConfigError(label + ": node " + std::to_string(node) + " is also a slave node of "
                          "contact pair '" + it->second + "'; a node takes one contact "
                          "constraint");
      }
      slave_owner[node] = pair.name;
      pair_slaves[k].insert(node);
    }

    // Mortar integrals of a deformable pair.
    std::map<Index, MortarRow> rows;
    if (!pair.rigid) {
      std::set<Index> masters;
      for (const SurfaceFace& face : master) masters.insert(face.nodes.begin(), face.nodes.end());
      for (Index node : masters) {
        if (weight.count(node) != 0) {
          throw ConfigError(label + ": node " + std::to_string(node) + " lies on both the slave "
                            "and the master surface; the two surfaces must be distinct");
        }
        all_masters.insert(node);
      }
      for (std::size_t f = 0; f < slave.size(); ++f) {
        const SurfaceFace& s = slave[f];
        const Scalar reach = options.search_factor * s.size;
        for (const SurfaceFace& m : master) {
          if ((m.lo.array() > (s.hi.array() + reach)).any() ||
              (m.hi.array() < (s.lo.array() - reach)).any()) {
            continue;
          }
          if (dim == 2) {
            integrate_segment_2d(s, face_normals[f], m, thickness, rows);
          } else {
            integrate_segment_3d(shape, s, face_normals[f], m, rows);
          }
        }
      }
    }

    // The nodes that take part.
    int uncovered = 0;
    int partial = 0;
    int overlapped = 0;
    int held = 0;
    int blocked = 0;
    for (const auto& [node, d] : weight) {
      Node c;
      c.node = node;
      c.pair = k;
      c.weight = d;
      c.friction = pair.friction;
      const Vector3 x = mesh.node(node);
      if (pair.rigid) {
        Vector3 nu;
        const Scalar g = pair.obstacle.gap(x, nu, dim);
        c.normal = nu;
        c.initial_gap = d * g;
        c.motion = pair.obstacle.motion;
        if (dim == 2) c.motion.z() = 0.0;
      } else {
        const auto it = rows.find(node);
        Scalar covered = 0.0;
        if (it != rows.end()) {
          for (const auto& [l, v] : it->second.m) covered += v;
        }
        const Scalar ratio = covered / d;
        if (it == rows.end() || ratio < 1.0e-6) {
          ++uncovered;
          continue;
        }
        if (ratio < kMinCoverage) {
          ++partial;
          continue;
        }
        if (ratio > kMaxCoverage) {
          ++overlapped;
          continue;
        }
        // Scale the row to sum_l M_jl = D_j exactly: a rigid translation of
        // both surfaces then leaves the gap unchanged (on a flat face the
        // factor is 1 to round-off; on a curved one it absorbs the
        // auxiliary plane's area error).
        const Scalar factor = d / covered;
        for (const auto& [l, v] : it->second.m) {
          if (std::abs(v) > 1.0e-14 * d) c.masters.emplace_back(l, v * factor);
        }
        c.initial_gap = it->second.gap * factor;
        c.normal = -slave_normal[node];  // the pressure pushes into the slave body
      }
      for (int comp = 0; comp < dim; ++comp) {
        if (!dofs.is_constrained(dofs.dof(node, comp))) c.free.push_back(comp);
      }
      if (c.free.empty()) {
        ++held;
        continue;
      }
      for (int comp : c.free) c.free_normal(comp) = c.normal(comp);
      if (c.free_normal.norm() < 1.0e-6) {
        ++blocked;
        continue;
      }
      // Tangential directions: the free axes, orthogonalised against the
      // free normal (the one most parallel to it dropped).
      const Vector3 fn = c.free_normal.normalized();
      std::vector<std::pair<Scalar, int>> axes;
      for (int comp : c.free) axes.emplace_back(std::abs(fn(comp)), comp);
      std::sort(axes.begin(), axes.end());
      axes.pop_back();
      for (const auto& [unused, comp] : axes) {
        (void)unused;
        Vector3 t = Vector3::Zero();
        t(comp) = 1.0;
        t -= t.dot(fn) * fn;
        for (const Vector3& prev : c.tangents) t -= t.dot(prev) * prev;
        const Scalar len = t.norm();
        if (len > 1.0e-10) c.tangents.push_back(t / len);
      }
      const Scalar h = dim == 2 ? d / thickness : std::sqrt(d);
      c.complementarity = options.complementarity * modulus[node] / h;
      nodes_.push_back(std::move(c));
    }
    excluded_.push_back(uncovered + partial + overlapped + held + blocked);
    const auto note = [&](int count, const std::string& what) {
      if (count == 0) return;
      std::ostringstream os;
      os << label << ": " << count << " slave node(s) " << what;
      exclusions_.push_back(os.str());
    };
    note(uncovered, "have no master surface opposite them within the search distance and "
                    "cannot come into contact");
    note(partial, "lie at the edge of the master surface, which covers them only in part, "
                  "and are left out of contact");
    note(overlapped, "see the master surface overlap itself in projection and are left out "
                     "of contact");
    note(held, "have every displacement component prescribed and are left out of contact");
    note(blocked, "are prescribed along their contact normal and are left out of contact");
  }
  for (Index node : all_masters) {
    const auto it = slave_owner.find(node);
    if (it != slave_owner.end()) {
      throw ConfigError("node " + std::to_string(node) + " is a slave node of contact pair '" +
                        it->second + "' and a master node of another; a node is one or the "
                        "other");
    }
  }
  if (nodes_.empty()) {
    throw ConfigError("contact: no slave node can come into contact (see the pairs' surfaces "
                      "and supports)");
  }
  slip_.assign(nodes_.size(), Vector3::Zero());
}

// ---------------------------------------------------------------------------
// Gap, tractions and slip
// ---------------------------------------------------------------------------
Scalar ContactProblem::weighted_gap(std::size_t i, const Vector& u, Scalar lambda) const {
  const Node& n = nodes_[i];
  const int dim = model_.dim();
  Scalar g = n.initial_gap;
  for (int k = 0; k < dim; ++k) {
    g += n.weight * n.normal(k) * (u(n.node * dim + k) - lambda * n.motion(k));
    for (const auto& [l, m] : n.masters) g -= m * n.normal(k) * u(l * dim + k);
  }
  return g;
}

void ContactProblem::tractions(const Node& n, const Vector& residual, Scalar& pressure,
                               Eigen::Vector2d& tangential) const {
  const int dim = model_.dim();
  Scalar dot = 0.0;
  for (int comp : n.free) dot += n.free_normal(comp) * residual(n.node * dim + comp);
  pressure = dot / (n.weight * n.free_normal.squaredNorm());
  tangential.setZero();
  for (std::size_t k = 0; k < n.tangents.size(); ++k) {
    Scalar s = 0.0;
    for (int comp : n.free) s += n.tangents[k](comp) * residual(n.node * dim + comp);
    tangential(static_cast<Eigen::Index>(k)) = s / n.weight;
  }
}

Eigen::Vector2d ContactProblem::weighted_slip(const Node& n, const Vector& u, Scalar lambda,
                                              const Vector& u_start, Scalar lambda_start) const {
  const int dim = model_.dim();
  Vector3 rel = Vector3::Zero();
  for (int k = 0; k < dim; ++k) {
    rel(k) = n.weight * (u(n.node * dim + k) - u_start(n.node * dim + k) -
                         (lambda - lambda_start) * n.motion(k));
    for (const auto& [l, m] : n.masters) rel(k) -= m * (u(l * dim + k) - u_start(l * dim + k));
  }
  Eigen::Vector2d out = Eigen::Vector2d::Zero();
  for (std::size_t k = 0; k < n.tangents.size(); ++k) {
    out(static_cast<Eigen::Index>(k)) = n.tangents[k].dot(rel);
  }
  return out;
}

// ---------------------------------------------------------------------------
// The condensed Newton system
// ---------------------------------------------------------------------------
ContactProblem::Linearization ContactProblem::linearize(const Vector& u, Scalar lambda,
                                                        const Vector& residual,
                                                        const SparseMatrix& tangent,
                                                        const Vector& u_start,
                                                        Scalar lambda_start) const {
  const DofManager& dofs = model_.dofs();
  const int dim = model_.dim();
  const std::vector<Index>& free_dofs = dofs.free_dofs();
  const Index nf = static_cast<Index>(free_dofs.size());
  using RowMajor = Eigen::SparseMatrix<Scalar, Eigen::RowMajor>;
  const RowMajor kr(detail::free_block(tangent, free_dofs));

  // Every free row is a combination of rows of K_ff plus explicit entries.
  struct Recipe {
    std::vector<std::pair<Index, Scalar>> rows;      ///< (row of K_ff, coefficient)
    std::vector<std::pair<Index, Scalar>> entries;   ///< (column, value)
    Scalar rhs = 0.0;
  };
  std::vector<Recipe> recipe(static_cast<std::size_t>(nf));
  for (Index r = 0; r < nf; ++r) {
    Recipe& rc = recipe[static_cast<std::size_t>(r)];
    rc.rows.emplace_back(r, 1.0);
    rc.rhs = -residual(free_dofs[static_cast<std::size_t>(r)]);
  }
  const auto red = [&](Index node, int comp) { return dofs.reduced_index(node * dim + comp); };

  Linearization out;
  out.status.assign(nodes_.size(), ContactStatus::Open);
  Scalar force2 = 0.0;
  for (std::size_t i = 0; i < nodes_.size(); ++i) {
    const Node& n = nodes_[i];
    Scalar p = 0.0;
    Eigen::Vector2d tau;
    tractions(n, residual, p, tau);
    const Scalar g = weighted_gap(i, u, lambda);
    const Scalar c = n.complementarity;
    const Scalar trial = p - c * g / n.weight;
    if (!(trial > 0.0)) continue;  // open: its rows stay the equilibrium p = 0
    const Scalar d = n.weight;
    const std::size_t m = n.tangents.size();
    const Scalar nf2 = n.free_normal.squaredNorm();

    // The status: in contact; with friction sticking or slipping.
    ContactStatus status = ContactStatus::Slip;
    Eigen::Vector2d lam_t = -tau;  // tangential traction along the slip
    Eigen::Vector2d slip = Eigen::Vector2d::Zero();
    Eigen::Vector2d v = Eigen::Vector2d::Zero();
    if (n.friction > 0.0) {
      slip = weighted_slip(n, u, lambda, u_start, lambda_start);
      v = lam_t + (c / d) * slip;
      // (Components beyond the node's m tangential directions are zero.)
      const Scalar bound = n.friction * trial;
      status = (m == 0 || v.norm() < bound) ? ContactStatus::Stick : ContactStatus::Slip;
    }
    out.status[i] = status;

    // Contact force of the node from its free residual components:
    // f = G R_F, G = nu nu_F^T / |nu_F|^2 (+ T T^T with friction).
    Eigen::Matrix3d gmat = n.normal * n.free_normal.transpose() / nf2;
    if (n.friction > 0.0) {
      for (const Vector3& t : n.tangents) gmat += t * t.transpose();
    }
    Vector3 rj = Vector3::Zero();
    for (int comp : n.free) rj(comp) = residual(n.node * dim + comp);
    // The scale of the condensed residual: the contact force, and the force
    // c g~ that closes the node's gap (all there is before contact settles).
    force2 += (gmat * rj).squaredNorm() + (c * g) * (c * g);

    // Master rows take the condensed contact force: R_l + sum (M/D) G R_j = 0.
    for (const auto& [l, mjl] : n.masters) {
      for (int a = 0; a < dim; ++a) {
        const Index r = red(l, a);
        if (r < 0) continue;
        Recipe& rc = recipe[static_cast<std::size_t>(r)];
        for (int b : n.free) {
          const Scalar coef = (mjl / d) * gmat(a, b);
          if (coef == 0.0) continue;
          rc.rows.emplace_back(red(n.node, b), coef);
          rc.rhs -= coef * residual(n.node * dim + b);
        }
      }
    }

    // The node's own rows: the gap, then one per tangential direction.
    std::vector<Index> own;
    for (int comp : n.free) own.push_back(red(n.node, comp));
    // Explicit entries of a kinematic row sum_b w_b (D du_jb - sum_l M_jl du_lb).
    const auto kinematic = [&](Recipe& rc, const Vector3& w, Scalar scale) {
      for (int b = 0; b < dim; ++b) {
        if (w(b) == 0.0) continue;
        const Index cj = red(n.node, b);
        if (cj >= 0) rc.entries.emplace_back(cj, scale * d * w(b));
        for (const auto& [l, mjl] : n.masters) {
          const Index cl = red(l, b);
          if (cl >= 0) rc.entries.emplace_back(cl, -scale * mjl * w(b));
        }
      }
    };
    {
      Recipe& rc = recipe[static_cast<std::size_t>(own[0])];
      rc.rows.clear();
      rc.entries.clear();
      kinematic(rc, n.normal, c);
      rc.rhs = -c * g;
    }
    for (std::size_t k = 0; k < m; ++k) {
      Recipe& rc = recipe[static_cast<std::size_t>(own[k + 1])];
      rc.rows.clear();
      rc.entries.clear();
      const Vector3& t = n.tangents[k];
      if (n.friction == 0.0) {
        // Frictionless: no tangential traction, t . R_j = 0.
        Scalar r0 = 0.0;
        for (int b : n.free) {
          rc.rows.emplace_back(red(n.node, b), t(b));
          r0 += t(b) * residual(n.node * dim + b);
        }
        rc.rhs = -r0;
      } else if (status == ContactStatus::Stick) {
        // Sticking: no slip over the step.
        kinematic(rc, t, c);
        rc.rhs = -c * slip(static_cast<Eigen::Index>(k));
      } else {
        // Slipping: lam_t = mu p v / |v|, linearised and multiplied by D.
        const Scalar vn = v.norm();
        const Eigen::Vector2d vh = v / vn;
        const Eigen::Matrix2d q = Eigen::Matrix2d::Identity() - vh * vh.transpose();
        const Scalar mu = n.friction;
        for (int b : n.free) {
          Scalar coef = -mu * vh(static_cast<Eigen::Index>(k)) * n.free_normal(b) / nf2;
          for (std::size_t i2 = 0; i2 < m; ++i2) {
            const Scalar delta = i2 == k ? 1.0 : 0.0;
            coef -= (delta - mu * p * q(static_cast<Eigen::Index>(k), static_cast<Eigen::Index>(i2)) /
                                 vn) *
                    n.tangents[i2](b);
          }
          if (coef != 0.0) rc.rows.emplace_back(red(n.node, b), coef);
        }
        Vector3 w = Vector3::Zero();
        for (std::size_t i2 = 0; i2 < m; ++i2) {
          w += q(static_cast<Eigen::Index>(k), static_cast<Eigen::Index>(i2)) * n.tangents[i2];
        }
        kinematic(rc, w, -mu * p * c / vn);
        rc.rhs = -d * (lam_t(static_cast<Eigen::Index>(k)) -
                       mu * p * vh(static_cast<Eigen::Index>(k)));
      }
    }
  }
  out.force_scale = std::sqrt(force2);

  TripletList triplets;
  out.rhs.resize(nf);
  for (Index r = 0; r < nf; ++r) {
    const Recipe& rc = recipe[static_cast<std::size_t>(r)];
    out.rhs(r) = rc.rhs;
    for (const auto& [src, coef] : rc.rows) {
      for (RowMajor::InnerIterator it(kr, src); it; ++it) {
        triplets.emplace_back(r, it.col(), coef * it.value());
      }
    }
    for (const auto& [col, value] : rc.entries) triplets.emplace_back(r, col, value);
  }
  out.matrix.resize(nf, nf);
  out.matrix.setFromTriplets(triplets.begin(), triplets.end());
  out.matrix.makeCompressed();
  return out;
}

// ---------------------------------------------------------------------------
// Results
// ---------------------------------------------------------------------------
std::vector<ContactNodeResult> ContactProblem::node_results(
    const Vector& u, Scalar lambda, const Vector& residual,
    const std::vector<ContactStatus>& status, const Vector& u_start,
    Scalar lambda_start) const {
  std::vector<ContactNodeResult> out;
  out.reserve(nodes_.size());
  for (std::size_t i = 0; i < nodes_.size(); ++i) {
    const Node& n = nodes_[i];
    ContactNodeResult r;
    r.node = n.node;
    r.pair = n.pair;
    r.weight = n.weight;
    r.normal = n.normal;
    r.gap = weighted_gap(i, u, lambda) / n.weight;
    r.status = status[i];
    r.slip = slip_[i];
    if (status[i] != ContactStatus::Open) {
      Scalar p = 0.0;
      Eigen::Vector2d tau;
      tractions(n, residual, p, tau);
      r.pressure = p;
      if (n.friction > 0.0) {
        for (std::size_t k = 0; k < n.tangents.size(); ++k) {
          r.traction += tau(static_cast<Eigen::Index>(k)) * n.tangents[k];
        }
      }
      const Eigen::Vector2d s = weighted_slip(n, u, lambda, u_start, lambda_start);
      for (std::size_t k = 0; k < n.tangents.size(); ++k) {
        r.slip += (s(static_cast<Eigen::Index>(k)) / n.weight) * n.tangents[k];
      }
    }
    out.push_back(r);
  }
  return out;
}

std::vector<ContactPairResult> ContactProblem::pair_results(
    const std::vector<ContactNodeResult>& nodes) const {
  std::vector<ContactPairResult> out(options_.pairs.size());
  for (std::size_t k = 0; k < out.size(); ++k) {
    out[k].name = options_.pairs[k].name;
    out[k].rigid = options_.pairs[k].rigid;
    out[k].friction = options_.pairs[k].friction;
    out[k].min_gap = std::numeric_limits<Scalar>::infinity();
  }
  for (const ContactNodeResult& r : nodes) {
    ContactPairResult& p = out[r.pair];
    ++p.nodes;
    p.min_gap = std::min(p.min_gap, r.gap);
    p.max_slip = std::max(p.max_slip, r.slip.norm());
    if (r.status == ContactStatus::Open) continue;
    ++p.active;
    if (p.friction > 0.0) {
      if (r.status == ContactStatus::Stick) ++p.sticking;
      if (r.status == ContactStatus::Slip) ++p.slipping;
    }
    p.area += r.weight;
    p.force += r.weight * (r.pressure * r.normal + r.traction);
    p.max_pressure = std::max(p.max_pressure, r.pressure);
  }
  for (std::size_t k = 0; k < out.size(); ++k) {
    if (out[k].nodes == 0) out[k].min_gap = 0.0;
    out[k].excluded = excluded_[k];
  }
  return out;
}

void ContactProblem::commit(const Vector& u, Scalar lambda,
                            const std::vector<ContactStatus>& status, const Vector& u_start,
                            Scalar lambda_start) {
  for (std::size_t i = 0; i < nodes_.size(); ++i) {
    if (status[i] == ContactStatus::Open) continue;
    const Node& n = nodes_[i];
    const Eigen::Vector2d s = weighted_slip(n, u, lambda, u_start, lambda_start);
    for (std::size_t k = 0; k < n.tangents.size(); ++k) {
      slip_[i] += (s(static_cast<Eigen::Index>(k)) / n.weight) * n.tangents[k];
    }
  }
}

}  // namespace sparlab


