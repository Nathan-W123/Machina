/// \file FaceGeometry.hpp
/// \brief Shape functions, quadrature and outward area vectors of element
///        boundary faces.
///
/// A boundary face of an element is itself an isoparametric entity: the edge
/// of a plane element is a two-node line, a Hex8 face a four-node bilinear
/// quadrilateral, a Tet4 face a flat three-node triangle and a Tet10 face a
/// six-node quadratic triangle. Surface loads (a normal pressure, a heat flux,
/// convection), contact and the load stiffness of a follower pressure are all
/// integrals over such faces, written once here against the face's own nodes.
///
/// **Orientation.** The face tables of Mesh.hpp wind every face so its
/// right-hand normal points out of the element. With face coordinates
/// \f$(s, t)\f$ and \f$x(s,t) = \sum_a N_a(s,t)\,x_a\f$, the outward area
/// vector per unit reference measure is
/// \f[
///   a = x_{,s} \times x_{,t} \quad\text{(3-D faces)},\qquad
///   a = (y_{,s},\, -x_{,s},\, 0) \quad\text{(2-D edges, per unit thickness)},
/// \f]
/// so \f$\int_\Gamma f\,n\,d\Gamma = \int f\,a\,ds\,dt\f$, and \f$|a|\f$ is the
/// area (length) element. For a 2-D edge wound counter-clockwise around its
/// element the right-hand normal of the direction of travel points outward.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/mesh/Mesh.hpp"

#include <vector>

namespace sparlab {

/// Boundary-face topologies.
enum class FaceShape {
  Line2,  ///< edge of a Q4 or Tri3: s in [-1, 1]
  Tri3,   ///< face of a Tet4: unit triangle, (s, t) >= 0, s + t <= 1
  Quad4,  ///< face of a Hex8: s, t in [-1, 1]
  Tri6    ///< face of a Tet10: corners, then edge nodes 0-1, 1-2, 2-0
};

/// The boundary-face topology of an element type.
/// \throws MeshError for an element type without faces.
FaceShape face_shape_of(ElementType type);

/// Nodes of a face topology.
int face_shape_nodes(FaceShape shape);

/// One quadrature point on the reference face (weight without the area
/// element).
struct FacePoint {
  Scalar s = 0.0;
  Scalar t = 0.0;
  Scalar weight = 0.0;
};

/// Quadrature rule on a reference face, exact for polynomials of degree
/// `2 * points - 1` per direction on lines and quadrilaterals, and with
/// `points` x `points` collapsed Gauss points on triangles (1 <= points <= 4).
std::vector<FacePoint> face_quadrature(FaceShape shape, int points);

/// Shape functions and their derivatives with respect to (s, t) at a point:
/// `n` has face_shape_nodes entries, `dn` is nodes x 2 (d/ds, d/dt; the t
/// column is zero on a line).
void face_shape_functions(FaceShape shape, Scalar s, Scalar t, Vector& n, Matrix& dn);

/// Geometry of a face at one reference point.
struct FaceGeometryPoint {
  Vector3 x = Vector3::Zero();      ///< position [m]
  Vector3 xs = Vector3::Zero();     ///< dx/ds [m]
  Vector3 xt = Vector3::Zero();     ///< dx/dt [m] (zero on a line)
  Vector3 area = Vector3::Zero();   ///< outward area vector a (see the file comment)
  Scalar jacobian = 0.0;            ///< |a| [m^2] (3-D) or [m] (2-D)
};

/// Evaluate the face geometry at (s, t) from the face's nodal coordinates
/// (dim x face nodes; 2-D coordinates are padded with z = 0).
FaceGeometryPoint face_geometry(FaceShape shape, const Matrix& face_coords, Scalar s,
                                Scalar t, const Vector& n, const Matrix& dn);

/// Nodal coordinates (dim x face nodes) of local face `local_face` of
/// element `e`, in the face table's order.
Matrix element_face_coordinates(const Mesh& mesh, Index element, int local_face);

/// Consistent nodal forces of a uniform normal pressure `p` [Pa] on one face,
/// positive pushing into the body: \f$f_a = -p\int N_a\,a\,ds\,dt\f$ (times the
/// thickness of a plane model). Returned per face node, 3 components each.
/// Exact on curved faces: the area vector is integrated, not a mean normal.
Matrix face_pressure_forces(FaceShape shape, const Matrix& face_coords, Scalar pressure,
                            Scalar thickness, int points);

/// Load stiffness of a follower pressure on one face in its current
/// configuration: \f$\partial f_a / \partial x_b\f$ (3 x 3 blocks, face nodes x
/// face nodes, as a (3 nf) x (3 nf) matrix), with f from face_pressure_forces.
/// A pressure that follows the deforming surface contributes
/// \f$-\partial f / \partial u\f$ to the tangent; the block is not symmetric
/// on an open surface.
Matrix face_pressure_stiffness(FaceShape shape, const Matrix& face_coords, Scalar pressure,
                               Scalar thickness, int points);

/// \f$\int_\Gamma N_a\,d\Gamma\f$ per face node (times the thickness in 2-D):
/// the consistent nodal weights of a uniform flux.
Vector face_shape_integrals(FaceShape shape, const Matrix& face_coords, Scalar thickness,
                            int points);

/// \f$\int_\Gamma N_a N_b\,d\Gamma\f$ (times the thickness in 2-D): the
/// consistent face "mass" of convection.
Matrix face_shape_products(FaceShape shape, const Matrix& face_coords, Scalar thickness,
                           int points);

}  // namespace sparlab
