#include "sparlab/elements/FaceGeometry.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/elements/Quadrature.hpp"

#include <Eigen/Geometry>

#include <sstream>

namespace sparlab {
namespace {

bool is_line(FaceShape shape) { return shape == FaceShape::Line2; }

Matrix3 skew(const Vector3& v) {
  Matrix3 m;
  m << 0.0, -v.z(), v.y(),
       v.z(), 0.0, -v.x(),
       -v.y(), v.x(), 0.0;
  return m;
}

/// Face coordinates padded to three rows.
Matrix padded(const Matrix& face_coords) {
  if (face_coords.rows() == 3) return face_coords;
  if (face_coords.rows() != 2) {
    throw MeshError("face coordinates must have two or three rows");
  }
  Matrix x = Matrix::Zero(3, face_coords.cols());
  x.topRows(2) = face_coords;
  return x;
}

void check_nodes(FaceShape shape, const Matrix& face_coords) {
  if (face_coords.cols() != face_shape_nodes(shape)) {
    std::ostringstream os;
    os << "a face of this shape has " << face_shape_nodes(shape) << " nodes, received "
       << face_coords.cols() << " coordinate columns";
    throw MeshError(os.str());
  }
}

}  // namespace

FaceShape face_shape_of(ElementType type) {
  switch (type) {
    case ElementType::Quad4:
    case ElementType::Tri3: return FaceShape::Line2;
    case ElementType::Hex8: return FaceShape::Quad4;
    case ElementType::Tet4: return FaceShape::Tri3;
    case ElementType::Tet10: return FaceShape::Tri6;
  }
  throw MeshError("element type has no boundary-face shape registered");
}

int face_shape_nodes(FaceShape shape) {
  switch (shape) {
    case FaceShape::Line2: return 2;
    case FaceShape::Tri3: return 3;
    case FaceShape::Quad4: return 4;
    case FaceShape::Tri6: return 6;
  }
  return 0;
}

std::vector<FacePoint> face_quadrature(FaceShape shape, int points) {
  if (points < 1 || points > 4) {
    std::ostringstream os;
    os << "face quadrature supports 1 to 4 points per direction, got " << points;
    throw ConfigError(os.str());
  }
  std::vector<FacePoint> rule;
  switch (shape) {
    case FaceShape::Line2:
      for (const QuadraturePoint1D& g : gauss_legendre_line(points)) {
        rule.push_back({g.xi, 0.0, g.weight});
      }
      break;
    case FaceShape::Quad4:
      for (const QuadraturePoint2D& g : gauss_legendre_square(points)) {
        rule.push_back({g.xi, g.eta, g.weight});
      }
      break;
    case FaceShape::Tri3:
    case FaceShape::Tri6:
      for (const QuadraturePoint2D& g : collapsed_gauss_triangle(points)) {
        rule.push_back({g.xi, g.eta, g.weight});
      }
      break;
  }
  return rule;
}

void face_shape_functions(FaceShape shape, Scalar s, Scalar t, Vector& n, Matrix& dn) {
  const int nf = face_shape_nodes(shape);
  n.resize(nf);
  dn.setZero(nf, 2);
  switch (shape) {
    case FaceShape::Line2:
      n << 0.5 * (1.0 - s), 0.5 * (1.0 + s);
      dn(0, 0) = -0.5;
      dn(1, 0) = 0.5;
      return;
    case FaceShape::Tri3:
      n << 1.0 - s - t, s, t;
      dn << -1.0, -1.0,
            1.0, 0.0,
            0.0, 1.0;
      return;
    case FaceShape::Quad4: {
      static constexpr Scalar c[4][2] = {{-1.0, -1.0}, {1.0, -1.0}, {1.0, 1.0}, {-1.0, 1.0}};
      for (int a = 0; a < 4; ++a) {
        n(a) = 0.25 * (1.0 + s * c[a][0]) * (1.0 + t * c[a][1]);
        dn(a, 0) = 0.25 * c[a][0] * (1.0 + t * c[a][1]);
        dn(a, 1) = 0.25 * c[a][1] * (1.0 + s * c[a][0]);
      }
      return;
    }
    case FaceShape::Tri6: {
      const Scalar l0 = 1.0 - s - t;
      const Scalar l1 = s;
      const Scalar l2 = t;
      n << l0 * (2.0 * l0 - 1.0), l1 * (2.0 * l1 - 1.0), l2 * (2.0 * l2 - 1.0),
          4.0 * l0 * l1, 4.0 * l1 * l2, 4.0 * l2 * l0;
      dn << -(4.0 * l0 - 1.0), -(4.0 * l0 - 1.0),
            4.0 * l1 - 1.0, 0.0,
            0.0, 4.0 * l2 - 1.0,
            4.0 * (l0 - l1), -4.0 * l1,
            4.0 * l2, 4.0 * l1,
            -4.0 * l2, 4.0 * (l0 - l2);
      return;
    }
  }
}

FaceGeometryPoint face_geometry(FaceShape shape, const Matrix& face_coords, Scalar s, Scalar t,
                                const Vector& n, const Matrix& dn) {
  (void)s;
  (void)t;
  const Matrix x = padded(face_coords);
  FaceGeometryPoint g;
  g.x = x * n;
  g.xs = x * dn.col(0);
  if (is_line(shape)) {
    // Per unit thickness: the right-hand normal of the direction of travel.
    g.area = Vector3(g.xs.y(), -g.xs.x(), 0.0);
  } else {
    g.xt = x * dn.col(1);
    g.area = g.xs.cross(g.xt);
  }
  g.jacobian = g.area.norm();
  return g;
}

Matrix element_face_coordinates(const Mesh& mesh, Index element, int local_face) {
  const std::vector<std::vector<int>>& faces = element_local_faces(mesh.element_type());
  if (local_face < 0 || local_face >= static_cast<int>(faces.size())) {
    std::ostringstream os;
    os << "local face " << local_face << " is outside the element's " << faces.size()
       << " faces";
    throw MeshError(os.str());
  }
  const std::vector<int>& fn = faces[static_cast<std::size_t>(local_face)];
  const Index* nodes = mesh.element_nodes(element);
  Matrix x(mesh.dim(), static_cast<Eigen::Index>(fn.size()));
  for (std::size_t a = 0; a < fn.size(); ++a) {
    x.col(static_cast<Eigen::Index>(a)) = mesh.coordinates().col(nodes[fn[a]]);
  }
  return x;
}

Matrix face_pressure_forces(FaceShape shape, const Matrix& face_coords, Scalar pressure,
                            Scalar thickness, int points) {
  check_nodes(shape, face_coords);
  const int nf = face_shape_nodes(shape);
  const Scalar depth = is_line(shape) ? thickness : 1.0;
  Matrix f = Matrix::Zero(3, nf);
  Vector n;
  Matrix dn;
  Scalar measure = 0.0;
  for (const FacePoint& p : face_quadrature(shape, points)) {
    face_shape_functions(shape, p.s, p.t, n, dn);
    const FaceGeometryPoint g = face_geometry(shape, face_coords, p.s, p.t, n, dn);
    measure += p.weight * g.jacobian;
    for (int a = 0; a < nf; ++a) f.col(a) -= (pressure * depth * p.weight * n(a)) * g.area;
  }
  if (!(measure > 0.0)) throw MeshError("zero-measure face in a pressure load");
  return f;
}

Matrix face_pressure_stiffness(FaceShape shape, const Matrix& face_coords, Scalar pressure,
                               Scalar thickness, int points) {
  check_nodes(shape, face_coords);
  const int nf = face_shape_nodes(shape);
  const Scalar depth = is_line(shape) ? thickness : 1.0;
  Matrix k = Matrix::Zero(3 * nf, 3 * nf);
  Vector n;
  Matrix dn;
  // d(a)/dx_b: for a line a = (y_s, -x_s, 0), so da = dN_b/ds R dx_b with R
  // the in-plane quarter turn; for a surface a = x_s x x_t, so
  // da = [dN_b/ds (-skew(x_t)) + dN_b/dt skew(x_s)] dx_b.
  Matrix3 quarter = Matrix3::Zero();
  quarter(0, 1) = 1.0;
  quarter(1, 0) = -1.0;
  for (const FacePoint& p : face_quadrature(shape, points)) {
    face_shape_functions(shape, p.s, p.t, n, dn);
    const FaceGeometryPoint g = face_geometry(shape, face_coords, p.s, p.t, n, dn);
    for (int b = 0; b < nf; ++b) {
      Matrix3 da;
      if (is_line(shape)) {
        da = dn(b, 0) * quarter;
      } else {
        da = dn(b, 0) * (-skew(g.xt)) + dn(b, 1) * skew(g.xs);
      }
      for (int a = 0; a < nf; ++a) {
        k.block<3, 3>(3 * a, 3 * b) -= (pressure * depth * p.weight * n(a)) * da;
      }
    }
  }
  return k;
}

Vector face_shape_integrals(FaceShape shape, const Matrix& face_coords, Scalar thickness,
                            int points) {
  check_nodes(shape, face_coords);
  const Scalar depth = is_line(shape) ? thickness : 1.0;
  Vector out = Vector::Zero(face_shape_nodes(shape));
  Vector n;
  Matrix dn;
  for (const FacePoint& p : face_quadrature(shape, points)) {
    face_shape_functions(shape, p.s, p.t, n, dn);
    const FaceGeometryPoint g = face_geometry(shape, face_coords, p.s, p.t, n, dn);
    out += (depth * p.weight * g.jacobian) * n;
  }
  return out;
}

Matrix face_shape_products(FaceShape shape, const Matrix& face_coords, Scalar thickness,
                           int points) {
  check_nodes(shape, face_coords);
  const Scalar depth = is_line(shape) ? thickness : 1.0;
  const int nf = face_shape_nodes(shape);
  Matrix out = Matrix::Zero(nf, nf);
  Vector n;
  Matrix dn;
  for (const FacePoint& p : face_quadrature(shape, points)) {
    face_shape_functions(shape, p.s, p.t, n, dn);
    const FaceGeometryPoint g = face_geometry(shape, face_coords, p.s, p.t, n, dn);
    out += (depth * p.weight * g.jacobian) * (n * n.transpose());
  }
  return out;
}

}  // namespace sparlab
