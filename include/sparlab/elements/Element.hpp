/// \file Element.hpp
/// \brief Abstract element interface.
///
/// Every element implementation supplies the kernels the rest of the library
/// needs: a stiffness matrix, a consistent mass matrix, a strain-displacement
/// operator at a requested parametric location, its stiffness quadrature
/// rule, and a consistent nodal load vector for a traction on one boundary
/// face. The interface is dimension-generic: coordinates, constitutive
/// matrices and the strain operator are dynamic Eigen matrices sized by
/// `dim()` and `num_voigt()`, so the assembler, stress recovery, load
/// application and the geometric stiffness of a buckling analysis are written
/// once for every topology: the plane Q4 and Tri3 and the solid Hex8, Tet4
/// and Tet10. Adding a new topology means adding one subclass and one entry
/// to `make_element`, plus its face table and file-format codes.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/mesh/Mesh.hpp"

#include <memory>
#include <string>
#include <vector>

namespace sparlab {

/// Element-local coordinates on the reference domain: the square / cube
/// [-1, 1]^dim for Q4 and Hex8, the unit triangle / tetrahedron for Tri3,
/// Tet4 and Tet10. `zeta` is ignored by plane elements.
struct NaturalPoint {
  Scalar xi = 0.0;
  Scalar eta = 0.0;
  Scalar zeta = 0.0;
};

/// Element formulation of the displacement elements.
enum class ElementFormulation {
  /// The isoparametric element with its compatible displacement field alone.
  Standard,
  /// Hex8 only: Wilson-Taylor incompatible modes, three condensed bubble
  /// fields per displacement component (Hex8Incompatible.hpp).
  IncompatibleModes
};

std::string to_string(ElementFormulation formulation);
/// "standard" or "incompatible_modes".
/// \throws ConfigError for any other text.
ElementFormulation parse_element_formulation(const std::string& text);

/// Quadrature orders used by the Q4 and Hex8 kernels. The linear simplices
/// (Tri3, Tet4) evaluate every kernel in closed form, exactly, and the Tet10
/// uses fixed simplex rules (Tet10.hpp); all three ignore these orders.
///
/// The element formulation travels with the orders because every kernel
/// that integrates an element receives them: `make_element` builds the
/// element of the formulation, and the kernels of the non-linear analysis
/// find it through `Element::num_internal_nodes`.
struct IntegrationOptions {
  int stiffness_points = 2;  ///< points per direction for K_e
  int mass_points = 3;       ///< points per direction for M_e
  int edge_points = 2;       ///< points per direction on a boundary face for tractions
  /// Hex8: points of the stiffness rule through the thickness of a sheet,
  /// along the natural axis `thickness_axis` (1 to kMaxThicknessPoints of
  /// Quadrature.hpp); 0 keeps `stiffness_points`. The rule is then
  /// `stiffness_points` x `stiffness_points` x `thickness_points`
  /// (`gauss_legendre_box`), in the point order of every other rule, and
  /// the per-point history of the non-linear analysis follows it.
  int thickness_points = 0;
  /// Natural axis of the element through the sheet thickness: 0 (xi), 1 (eta)
  /// or 2 (zeta), which the structured hexahedral generator aligns with x, y
  /// and z.
  int thickness_axis = 2;
  ElementFormulation formulation = ElementFormulation::Standard;
};

/// Strain-displacement data evaluated at one parametric point.
struct StrainOperator {
  Matrix b;           ///< num_voigt x num_dofs operator [1/m]
  Scalar detJ = 0.0;  ///< Jacobian determinant [m^dim]
};

/// One point of an element's stiffness quadrature rule.
struct IntegrationPoint {
  NaturalPoint point;
  Scalar weight = 0.0;  ///< weight on the reference domain, without det J
};

/// Static condensation of an element's internal (incompatible) modes for a
/// linear-elastic constitutive matrix D. With the internal parameters
/// \f$\alpha\f$ and their strain operator \f$\tilde B\f$
/// (`Element::internal_strain_operator`),
/// \f[
///   K_{\alpha\alpha} = \int \tilde B^T D \tilde B\,dV ,\qquad
///   K_{\alpha u} = \int \tilde B^T D B\,dV ,\qquad
///   \alpha = K_{\alpha\alpha}^{-1}\Big(\int \tilde B^T D\,\varepsilon_0\,dV
///             - K_{\alpha u} u_e\Big)
/// \f]
/// for an element with the initial (thermal) strain \f$\varepsilon_0\f$,
/// so that the strain \f$Bu_e + \tilde B\alpha\f$ is
/// \f$\hat B u_e + \tilde B K_{\alpha\alpha}^{-1}\int\tilde B^T
/// D\varepsilon_0\f$ with \f$\hat B = B - \tilde B\,K_{\alpha\alpha}^{-1}
/// K_{\alpha u}\f$, and \f$\int\hat B^T D\hat B = K^*\f$, the condensed
/// stiffness. Empty (0 rows) for an element without internal modes.
struct InternalCondensation {
  Matrix coupling;  ///< \f$K_{\alpha\alpha}^{-1}K_{\alpha u}\f$, num_internal_dofs x num_dofs [-]
  Matrix inverse;   ///< \f$K_{\alpha\alpha}^{-1}\f$ [m/N]
};

/// Abstract continuum element.
class Element {
 public:
  virtual ~Element() = default;

  virtual ElementType type() const = 0;
  /// Spatial dimension (2 or 3): the translational DOFs per node.
  virtual int dim() const = 0;
  virtual int num_nodes() const = 0;
  /// Number of boundary faces (edges in 2-D).
  virtual int num_faces() const = 0;

  /// Degrees of freedom per node, in node-major order: the `dim` translations
  /// of a continuum element; the three translations followed by the three
  /// rotations about the global axes for a shell or beam element.
  virtual int dofs_per_node() const { return dim(); }

  /// Lump the mass matrix by scaling its diagonal to the element mass per
  /// DOF component (Hinton, Rock and Zienkiewicz) rather than by row sums:
  /// the quadratic tetrahedron, whose corner row sums are negative, and the
  /// structural elements, whose rotational rows carry rotary inertia.
  virtual bool diagonal_scaled_lumping() const { return false; }

  /// Natural coordinates of the element centroid: the origin of the
  /// reference square / cube, (1/3, 1/3) on the triangle, (1/4, 1/4, 1/4) on
  /// the tetrahedron.
  virtual NaturalPoint reference_centroid() const { return NaturalPoint(); }

  /// Degrees of freedom carried by the element (num_nodes * dofs_per_node).
  int num_dofs() const { return num_nodes() * dofs_per_node(); }

  /// Voigt components of the element's strain operator.
  int num_voigt() const { return voigt_components(dim()); }

  /// Element stiffness matrix
  /// \f$ K_e = \int_{\Omega_e} t\, B^T D B \, d\Omega \f$ [N/m].
  /// \param coords dim x num_nodes nodal coordinates [m].
  /// \param d constitutive matrix, num_voigt x num_voigt [Pa].
  /// \param thickness out-of-plane thickness [m] (must be 1 for 3-D elements).
  virtual Matrix stiffness(const Matrix& coords, const Matrix& d, Scalar thickness,
                           const IntegrationOptions& opts) const = 0;

  /// Consistent element mass matrix
  /// \f$ M_e = \int_{\Omega_e} \rho\, t\, N^T N \, d\Omega \f$ [kg].
  virtual Matrix consistent_mass(const Matrix& coords, Scalar density, Scalar thickness,
                                 const IntegrationOptions& opts) const = 0;

  /// Strain-displacement operator and Jacobian determinant at a natural point.
  /// \throws MeshError when detJ <= 0 (inverted or degenerate element).
  virtual StrainOperator strain_operator(const Matrix& coords,
                                         const NaturalPoint& point) const = 0;

  /// Shape function values at a natural point (size num_nodes).
  virtual Vector shape_functions(const NaturalPoint& point) const = 0;

  /// Parametric locations of the stiffness quadrature points, in the order used
  /// by `stiffness`. Needed for stress recovery at Gauss points.
  virtual std::vector<NaturalPoint> stress_evaluation_points(
      const IntegrationOptions& opts) const = 0;

  /// The stiffness quadrature rule itself: the points of
  /// `stress_evaluation_points` with their reference weights, so that
  /// \f$\sum_g w_g \det J_g\, f(\xi_g)\f$ integrates \f$f\f$ exactly as
  /// `stiffness` does.
  virtual std::vector<IntegrationPoint> integration_rule(
      const IntegrationOptions& opts) const = 0;

  /// Geometric (initial-stress) stiffness of the element in the stress state
  /// of the element displacement `ue`:
  /// \f[
  ///   K_{G,e} = \int_{\Omega_e} t\, G^T \sigma\, G \,d\Omega \otimes I_{dim},
  ///   \qquad \sigma = s\,D B u_e ,
  /// \f]
  /// where \f$G\f$ holds the shape-function gradients (dim x num_nodes) and
  /// \f$\sigma\f$ is the stress tensor at each stiffness integration point,
  /// so \f$\phi^T K_{G,e}\phi = \int t \sum_k \nabla\phi_k^T\sigma\nabla\phi_k\f$
  /// is the second-order work of the stress on the rotations of a mode
  /// \f$\phi\f$. \f$K_{G,e}\f$ is linear in \f$u_e\f$ and indefinite. `t` is
  /// the thickness of a plane element and 1 for a solid.
  /// \param d constitutive matrix the stress is computed with [Pa].
  /// \param stress_scale factor \f$s\f$ on the stress (1 for a plain analysis).
  /// Shell and beam elements override it with their own kinematics.
  virtual Matrix geometric_stiffness(const Matrix& coords, const Matrix& d, const Vector& ue,
                                     Scalar stress_scale, Scalar thickness,
                                     const IntegrationOptions& opts) const;

  /// Geometric stiffness of a given stress field: `stresses` holds the Voigt
  /// stress at each point of `integration_rule(opts)`, in its order. This is
  /// the form a prestress that is not \f$D B u_e\f$ needs - a thermal one,
  /// \f$D(Bu_e - \varepsilon_0)\f$, say. Continuum elements only.
  Matrix geometric_stiffness_of_stress(const Matrix& coords, const std::vector<Vector>& stresses,
                                       Scalar thickness, const IntegrationOptions& opts) const;

  /// Derivative of \f$\phi_e^T K_{G,e}(u_e)\phi_e\f$ with respect to \f$u_e\f$.
  /// Because \f$K_{G,e}\f$ is linear in \f$u_e\f$ this is the vector
  /// \f$g_e\f$ with \f$\phi_e^T K_{G,e}(u_e)\phi_e = g_e^T u_e\f$:
  /// \f[
  ///   g_e = s \int_{\Omega_e} t\, B^T D\, \hat\Phi \, d\Omega,
  ///   \qquad \Phi_{ij} = \sum_k \frac{\partial\phi_k}{\partial x_i}
  ///   \frac{\partial\phi_k}{\partial x_j},
  /// \f]
  /// with \f$\hat\Phi\f$ the Voigt vector of \f$\Phi\f$ with doubled shear
  /// entries. It is the adjoint load of a buckling-load sensitivity.
  virtual Vector geometric_stiffness_derivative(const Matrix& coords, const Matrix& d,
                                                const Vector& phi, Scalar stress_scale,
                                                Scalar thickness,
                                                const IntegrationOptions& opts) const;

  /// Internal "pseudo-nodes" of the element: fields that are condensed
  /// inside the element (static condensation) and carry `dim()` parameters
  /// each, like a node - 3 for the incompatible-mode Hex8, 0 for every
  /// other element.
  virtual int num_internal_nodes() const { return 0; }

  /// Internal parameters condensed in the element (num_internal_nodes * dim).
  int num_internal_dofs() const { return num_internal_nodes() * dim(); }

  /// Gradients of the internal fields with respect to the reference
  /// coordinates at a natural point, dim x num_internal_nodes [1/m] (a
  /// dim x 0 matrix for an element without internal modes): the columns
  /// that extend the shape-function gradients \f$G\f$ of the element's
  /// nodes, so the displacement gradient is \f$H = U G^T + A\tilde G^T\f$
  /// with the internal parameters \f$A\f$ (dim x num_internal_nodes).
  /// \throws MeshError when det J <= 0.
  virtual Matrix internal_mode_gradients(const Matrix& coords, const NaturalPoint& point) const;

  /// The small-strain operator of the internal parameters at a natural
  /// point, num_voigt x num_internal_dofs [1/m], built from
  /// `internal_mode_gradients` as `strain_operator` is from the shape
  /// gradients (parameter k of pseudo-node m is column dim * m + k).
  Matrix internal_strain_operator(const Matrix& coords, const NaturalPoint& point) const;

  /// Condensation of the internal parameters of a linear-elastic element
  /// with constitutive matrix `d` over the stiffness rule (see
  /// InternalCondensation); empty for an element without internal modes.
  /// \throws SolverError when \f$K_{\alpha\alpha}\f$ is singular.
  InternalCondensation condense_internal(const Matrix& coords, const Matrix& d,
                                         Scalar thickness, const IntegrationOptions& opts) const;

  /// The strain operator with the internal parameters condensed,
  /// \f$\hat B = B - \tilde B\,K_{\alpha\alpha}^{-1}K_{\alpha u}\f$, at a
  /// natural point (the plain B of an element without internal modes).
  /// Never pass it where the shape-function gradients are read off B.
  Matrix condensed_strain_operator(const Matrix& coords, const NaturalPoint& point,
                                   const InternalCondensation& condensation) const;

  /// Consistent nodal forces for a constant traction on local face
  /// `local_face` (an edge in 2-D):
  /// \f$ f_e = \int_{\Gamma_e} t\, N^T \bar{t} \, d\Gamma \f$ [N].
  /// \param traction traction vector [Pa]; the z component must be zero for a
  ///        plane element.
  /// \return vector of length num_dofs [N].
  virtual Vector boundary_traction(const Matrix& coords, int local_face,
                                   const Vector3& traction, Scalar thickness,
                                   const IntegrationOptions& opts) const = 0;

  /// Local node indices of boundary face `local_face`, from the topology's
  /// shared face table.
  /// \throws MeshError for an out-of-range face index.
  const std::vector<int>& face_nodes(int local_face) const;
};

/// Factory for the supported element topologies.
/// \throws ConfigError for an unsupported type.
std::unique_ptr<Element> make_element(ElementType type);

/// The element of the formulation `opts.formulation` for a topology.
/// \throws ConfigError for an unsupported type, or a formulation the
///         topology does not have (incompatible modes: Hex8 only).
std::unique_ptr<Element> make_element(ElementType type, const IntegrationOptions& opts);

}  // namespace sparlab
