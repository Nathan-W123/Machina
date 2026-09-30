/// \file Loads.hpp
/// \brief Volume loads - self-weight, body forces, rotation - and thermal
///        loads of a continuum model.
///
/// **Body loads from the consistent mass.** A body force density that is
/// affine in position, \f$b(x) = b_0 + B x\f$, is interpolated exactly by the
/// shape functions of any isoparametric element,
/// \f$b(x(\xi)) = \sum_b N_b(\xi)\,b(x_b)\f$, because the element maps
/// \f$x(\xi) = \sum_b N_b(\xi)\,x_b\f$ exactly. Its consistent nodal forces are
/// therefore
/// \f[
///   f_a = \int_{\Omega_e} N_a\,b\,dV = \sum_b \Big(\int_{\Omega_e} N_a N_b\,dV\Big)\,b(x_b)
///       = (M_e^{(1)}\,\hat b)_a ,
/// \f]
/// the unit-density consistent mass matrix times the nodal values of \f$b\f$ -
/// exact to the mass matrix's own quadrature, on curved Tet10 cells too.
/// Self-weight is the constant \f$b = \rho g\f$, a uniform body force density is
/// constant, and the centrifugal load of a rotation \f$\omega\f$ about an axis
/// through \f$c\f$ with direction \f$e\f$ is the affine
/// \f$b = \rho\,\omega^2 (I - e e^T)(x - c)\f$. The resultant of self-weight is
/// the model's mass times \f$g\f$ to round-off, which the tests check.
///
/// **Thermal loads.** A temperature change \f$\Delta T(x) = N^T (T_e - T_{ref})\f$
/// induces the free strain \f$\varepsilon_0 = \alpha\Delta T\f$ (Voigt form per
/// idealisation, IsotropicMaterial.hpp). With \f$\sigma = D(Bu - \varepsilon_0)\f$
/// equilibrium reads \f$K u = f + f_{th}\f$ with
/// \f[
///   f_{th} = \int_{\Omega_e} B^T D\,\varepsilon_0\,t\,dV ,
/// \f]
/// integrated with the stiffness rule of the element, so that a field of
/// free expansion that the element can represent (a uniform or linear
/// \f$\Delta T\f$ on affine cells) is reproduced exactly with zero stress.
/// An element with condensed internal modes (the incompatible-mode Hex8)
/// takes the condensed operator, \f$f_{th} = \int\hat B^TD\varepsilon_0\f$
/// (Element.hpp, InternalCondensation); the modes take no load from a
/// uniform \f$\varepsilon_0\f$, whose \f$\int\tilde B^TD\varepsilon_0\f$
/// vanishes.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/fem/FemModel.hpp"

namespace sparlab {

/// Consistent nodal forces [N] of the body loads of a load case: gravity,
/// uniform body force densities and the centrifugal load. Zero-length
/// output never happens; a case without body loads gets a zero vector.
/// \throws ConfigError for a body-force region that selects no element, or a
///         centrifugal axis that is not a direction.
Vector assemble_body_load_vector(const FemModel& model, const LoadCaseSpec& spec);

/// Nodal temperatures of a `Uniform` or `Regions` temperature field: the
/// uniform value everywhere, then each region's value on its nodes, later
/// regions winning.
/// \throws ConfigError for a region that selects no node.
Vector resolve_region_temperatures(const Mesh& mesh, const TemperatureSpec& spec);

/// Thermal equivalent load of a nodal temperature field.
struct ThermalLoad {
  Vector force;              ///< f_th [N]
  /// 1/2 int eps0^T D eps0 dV [J], less, for an element with internal
  /// modes, the energy their relaxation releases,
  /// \f$\tfrac12 f_\alpha^T K_{\alpha\alpha}^{-1} f_\alpha\f$ with
  /// \f$f_\alpha = \int\tilde B^T D\varepsilon_0\f$: the elastic energy
  /// is then \f$\tfrac12 u^TK^*u - f_{th}^Tu + \f$ it.
  Scalar self_energy = 0.0;
};

/// Assemble \f$f_{th}\f$ and the self energy for nodal temperatures
/// `temperature` [K] (each element's material supplies alpha and T_ref).
/// \param stiffness_scale optional per-element factors on D (a SIMP design).
ThermalLoad assemble_thermal_load(const FemModel& model, const Vector& temperature,
                                  const Vector* stiffness_scale = nullptr);

/// Thermal strain \f$\varepsilon_0\f$ (Voigt, 3 or 6 components) at a natural
/// point of an element for nodal temperatures `temperature` [K].
Vector element_thermal_strain(const FemModel& model, Index element, const NaturalPoint& point,
                              const Vector& temperature);

/// The internal (incompatible-mode) parameters of element e of a linear
/// elastic model at the element displacement `ue`, with the thermal strain
/// of the nodal temperatures `temperature` (nullptr: none):
/// \f$\alpha = K_{\alpha\alpha}^{-1}(\int\tilde B^TD\varepsilon_0\,dV -
/// K_{\alpha u}u_e)\f$ (Element.hpp, InternalCondensation). Empty for an
/// element without internal modes. A stiffness factor on D cancels.
Vector linear_internal_parameters(const FemModel& model, Index e, const Vector& ue,
                                  const Vector* temperature);

/// The small strain \f$Bu_e + \tilde B\alpha\f$ at a natural point of
/// element e, with the internal parameters of `linear_internal_parameters`
/// (empty: the compatible \f$Bu_e\f$ alone).
Vector linear_point_strain(const FemModel& model, Index e, const NaturalPoint& point,
                           const Vector& ue, const Vector& alpha);

/// Temperature change \f$T - T_{ref}\f$ at a natural point of an element.
Scalar element_temperature_change(const FemModel& model, Index element,
                                  const NaturalPoint& point, const Vector& temperature);

}  // namespace sparlab
