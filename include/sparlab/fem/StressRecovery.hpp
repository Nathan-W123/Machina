/// \file StressRecovery.hpp
/// \brief Element strain / stress recovery and derived invariants.
///
/// Strain is obtained from the element's strain-displacement operator,
/// \f$ \boldsymbol{\varepsilon} = \mathbf{B}(\xi,\eta[,\zeta])\, \mathbf{u}_e \f$,
/// and stress from the constitutive law, \f$ \boldsymbol{\sigma} = s_e \mathbf{D}
/// \boldsymbol{\varepsilon} \f$, where \f$ s_e \f$ is the element stiffness
/// scale factor (1 for a plain analysis, \f$E(\rho_e)/E_0\f$ for a SIMP
/// design). With \f$ s_e \f$ applied the result is the *macroscopic* stress
/// carried by the element; the unscaled value is the stress in the underlying
/// solid material and is reported alongside it.
///
/// Reported element values are the average over the stiffness quadrature
/// points. Nodal values are measure-weighted averages of the adjacent element
/// values (simple averaging, not superconvergent patch recovery), and are used
/// for smooth contour plots only.
///
/// **Thermal strain.** With a temperature field the stress is
/// \f$\sigma = s_e D (B u_e - \varepsilon_0)\f$ (IsotropicMaterial.hpp); the
/// reported strain stays the total strain \f$B u_e\f$, and the element strain
/// energy is the elastic one, \f$\tfrac12\int(Bu-\varepsilon_0)^T D
/// (Bu-\varepsilon_0)\,dV\f$.
///
/// von Mises stress uses the correct out-of-plane stress for the active
/// idealisation: \f$\sigma_{zz} = 0\f$ for plane stress,
/// \f$\sigma_{zz} = \nu(\sigma_{xx}+\sigma_{yy}) - s_e E\alpha\Delta T\f$
/// for plane strain and the computed \f$\sigma_{zz}\f$ in 3-D, so
/// \f[
///   \sigma_{vm} = \sqrt{\tfrac{1}{2}\left[(\sigma_{xx}-\sigma_{yy})^2 +
///   (\sigma_{yy}-\sigma_{zz})^2 + (\sigma_{zz}-\sigma_{xx})^2\right] +
///   3\left(\sigma_{xy}^2 + \sigma_{yz}^2 + \sigma_{zx}^2\right)}.
/// \f]
/// Principal stresses are the roots of the 2-D Mohr circle or, in 3-D, the
/// eigenvalues of the symmetric stress tensor.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/fem/Assembler.hpp"
#include "sparlab/fem/FemModel.hpp"

namespace sparlab {

/// Voigt stress/strain storage: num_voigt x N columns, one per entity.
using VoigtField = Matrix;

struct StressField {
  VoigtField element_strain;            ///< nv x num_elements [-]
  VoigtField element_stress;            ///< nv x num_elements, macroscopic [Pa]
  VoigtField element_solid_stress;      ///< nv x num_elements, solid material [Pa]
  Vector element_von_mises;             ///< num_elements, from element_stress [Pa]
  Vector element_solid_von_mises;       ///< num_elements, from solid stress [Pa]
  Vector element_principal_max;         ///< num_elements, sigma_1 [Pa]
  Vector element_principal_mid;         ///< num_elements, sigma_2 [Pa] (3-D only, else empty)
  Vector element_principal_min;         ///< num_elements, smallest principal stress [Pa]
  Vector element_strain_energy;         ///< num_elements [J]
  /// Plane strain only: the out-of-plane stress sigma_zz per element [Pa].
  Vector element_sigma_zz;
  /// With a temperature field: the element's mean temperature [K].
  Vector element_temperature;

  VoigtField nodal_stress;              ///< nv x num_nodes, averaged [Pa]
  Vector nodal_von_mises;               ///< num_nodes [Pa]
};

/// von Mises stress from a Voigt stress vector (3 or 6 components).
Scalar von_mises(const Vector& voigt_stress, StressState state, Scalar poisson);

/// von Mises stress of a plane state with an explicit out-of-plane stress.
Scalar von_mises_plane(Scalar sxx, Scalar syy, Scalar sxy, Scalar szz);

/// In-plane principal stresses (sigma_1 >= sigma_2) of a 2-D Voigt vector.
void principal_stresses(const Vector& voigt_stress, Scalar& s1, Scalar& s2);

/// Principal stresses (descending) of a 6-component Voigt vector.
Vector3 principal_stresses_3d(const Vector& voigt_stress);

/// Strain at one parametric point of one element.
Vector element_strain_at(const FemModel& model, Index element, const NaturalPoint& point,
                         const Vector& displacement);

/// Stress at one parametric point of one element (macroscopic, scaled).
/// \param temperature optional nodal temperatures [K] (thermal strain).
Vector element_stress_at(const FemModel& model, Index element, const NaturalPoint& point,
                         const Vector& displacement, Scalar stiffness_scale = 1.0,
                         const Vector* temperature = nullptr);

/// Recover all strain / stress data for a displacement field.
/// \param stiffness_scale optional per-element SIMP factors \f$E(\rho)/E_0\f$.
/// \param temperature optional nodal temperatures [K] of the load case.
StressField recover_stresses(const FemModel& model, const Assembler& assembler,
                             const Vector& displacement,
                             const Vector* stiffness_scale = nullptr,
                             const Vector* temperature = nullptr);

}  // namespace sparlab
