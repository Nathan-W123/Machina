/// \file Tet10.hpp
/// \brief Ten-node quadratic tetrahedron (solid element).
///
/// Nodes 0-3 are the corners of the reference tetrahedron (0,0,0), (1,0,0),
/// (0,1,0), (0,0,1); nodes 4-9 sit on the edges 0-1, 1-2, 2-0, 0-3, 1-3,
/// 2-3, which is the VTK (VTK_QUADRATIC_TETRA), Abaqus / CalculiX (C3D10)
/// and scikit-fem order. Gmsh numbers the last two edge nodes the other way
/// round; the mesh reader swaps them. With the barycentric coordinates
/// \f$L_0 = 1-\xi-\eta-\zeta,\ L_1 = \xi,\ L_2 = \eta,\ L_3 = \zeta\f$:
/// \f[
///   N_i = L_i(2L_i - 1)\ \ (i = 0..3), \qquad
///   N_{ab} = 4 L_a L_b\ \ \text{on edge } a\text{-}b .
/// \f]
/// The map is isoparametric, so edge nodes may lie off the straight edge and
/// follow a curved CAD surface (a Gmsh `Mesh.ElementOrder = 2` mesh does
/// this). Strain varies linearly over a straight-sided element, which is what
/// frees the element from the shear locking of the constant-strain Tet4: it
/// reproduces any quadratic displacement field - pure bending among them -
/// exactly.
///
/// Integration:
///   * stiffness   the symmetric 4-point rule (degree 2), exact for a
///                 straight-sided element and the rule of Abaqus' C3D10;
///   * mass        the 64-point collapsed Gauss rule (degree 5), exact for
///                 \f$N^T N\f$ on a straight-sided element;
///   * face loads  a 6-node triangular face integrated with the 16-point
///                 collapsed Gauss rule on \f$N^T \bar t\,|x_{,r}\times x_{,s}|\f$,
///                 which puts \f$S\bar t/3\f$ on each edge node and nothing on
///                 the corners of a flat face of area \f$S\f$.
///
/// The consistent mass matrix of the quadratic tetrahedron has negative row
/// sums at the corners, so a row-sum lumped mass would be indefinite; the
/// assembler uses Hinton-Rock-Zienkiewicz lumping (the scaled diagonal) for
/// this element instead. Voigt order, engineering shear strains, node-major
/// DOFs and the unit-thickness rule are those of the Hex8 and Tet4.
#pragma once

#include "sparlab/elements/Element.hpp"

namespace sparlab {

/// Shape functions (size 10) of the quadratic tetrahedron at a natural point.
Eigen::Matrix<Scalar, 10, 1> tet10_shape_functions(Scalar xi, Scalar eta, Scalar zeta);

/// Natural gradients (row a = dN_a / d(xi, eta, zeta)) at a natural point.
Eigen::Matrix<Scalar, 10, 3> tet10_shape_gradients_natural(Scalar xi, Scalar eta,
                                                           Scalar zeta);

/// Volume of a quadratic tetrahedron from its 3 x 10 nodal coordinates,
/// integrated with the 27-point collapsed rule (exact: det J is a cubic).
/// `min_det`, when given, receives the smallest Jacobian determinant over the
/// stiffness integration points, so a cell folded by a badly curved edge is
/// detectable even when its volume is positive.
Scalar tet10_volume(const Matrix& coords, Scalar* min_det = nullptr);

/// Ratio of the smallest to the largest Jacobian determinant over the ten
/// nodes: 1 for a straight-sided element (affine map), below 1 as curved
/// edges distort it, and <= 0 once the map folds at a node.
Scalar tet10_jacobian_ratio(const Matrix& coords);

/// Area of a 6-node triangular face (corners 0-2, then the edge nodes of
/// 0-1, 1-2, 2-0) from its 3 x 6 nodal coordinates.
Scalar tri6_face_area(const Matrix& face_coords);

class Tet10Element final : public Element {
 public:
  ElementType type() const override { return ElementType::Tet10; }
  int dim() const override { return 3; }
  int num_nodes() const override { return 10; }
  int num_faces() const override { return 4; }
  bool diagonal_scaled_lumping() const override { return true; }
  NaturalPoint reference_centroid() const override;

  Matrix stiffness(const Matrix& coords, const Matrix& d, Scalar thickness,
                   const IntegrationOptions& opts) const override;

  Matrix consistent_mass(const Matrix& coords, Scalar density, Scalar thickness,
                         const IntegrationOptions& opts) const override;

  StrainOperator strain_operator(const Matrix& coords,
                                 const NaturalPoint& point) const override;

  Vector shape_functions(const NaturalPoint& point) const override;

  std::vector<NaturalPoint> stress_evaluation_points(
      const IntegrationOptions& opts) const override;

  std::vector<IntegrationPoint> integration_rule(
      const IntegrationOptions& opts) const override;

  Vector boundary_traction(const Matrix& coords, int local_face, const Vector3& traction,
                           Scalar thickness, const IntegrationOptions& opts) const override;
};

}  // namespace sparlab
