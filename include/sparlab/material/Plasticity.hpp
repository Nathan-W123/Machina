/// \file Plasticity.hpp
/// \brief J2 (von Mises) plasticity: the backward-Euler radial return and its
///        consistent tangent (Simo and Hughes 1998, boxes 3.1 and 3.2), in
///        3-D, plane strain and plane stress.
///
/// **The law.** The strain splits additively, \f$\varepsilon = \varepsilon^e
/// + \varepsilon^p + \varepsilon^\theta\f$ (the thermal strain
/// \f$\alpha\Delta T\,I\f$), the stress is \f$\sigma = D\varepsilon^e\f$, and
/// with the relative stress \f$\xi = \mathrm{dev}\,\sigma - \beta\f$ the
/// yield function is
/// \f[
///   f = \|\xi\| - \sqrt{2/3}\,\sigma_y(\bar\alpha) \le 0 ,
/// \f]
/// with associative flow \f$\dot\varepsilon^p = \dot\gamma\,n\f$,
/// \f$n = \xi/\|\xi\|\f$, the accumulated plastic strain
/// \f$\dot{\bar\alpha} = \sqrt{2/3}\,\dot\gamma\f$ and the back stress
/// \f$\dot\beta = \tfrac{2}{3}H_{kin}\,\dot\gamma\,n\f$
/// (IsotropicMaterial.hpp has \f$\sigma_y\f$).
///
/// **Integration.** From the state of the last converged step, the elastic
/// trial stress \f$s^{tr} = 2G\,\mathrm{dev}(\varepsilon - \varepsilon^p_n)\f$
/// either stays inside the yield surface or returns radially onto it:
/// \f$\Delta\gamma\f$ solves
/// \f[
///   \|\xi^{tr}\| - (2G + \tfrac{2}{3}H_{kin})\,\Delta\gamma
///   - \sqrt{2/3}\,\sigma_y(\bar\alpha_n + \sqrt{2/3}\,\Delta\gamma) = 0
/// \f]
/// (in closed form for linear hardening, by Newton with Voce saturation),
/// and \f$\sigma = \sigma^{tr} - 2G\Delta\gamma\,n\f$. The consistent tangent
/// is
/// \f[
///   C = K\,1\otimes 1 + 2G\theta\,(I - \tfrac{1}{3}1\otimes 1)
///       - 2G\bar\theta\,n\otimes n ,\quad
///   \theta = 1 - \frac{2G\Delta\gamma}{\|\xi^{tr}\|} ,\quad
///   \bar\theta = \frac{1}{1 + (\sigma_y' + H_{kin})/(3G)} - (1 - \theta) .
/// \f]
/// The return is exact for linear hardening under a strain increment that
/// keeps its direction, and first-order accurate otherwise.
///
/// **Idealisations.** Every point carries the full 3-D state (strain with
/// engineering shears, stress tensorial, in the Voigt order
/// \f$\{11, 22, 33, 12, 23, 31\}\f$). Plane strain returns with the given
/// \f$\varepsilon_{33}\f$ (0, or the B-bar value of Elastoplastic.hpp).
/// Plane stress finds \f$\varepsilon_{33}\f$ by Newton so that
/// \f$\sigma_{33} = 0\f$ and condenses it out of the tangent,
/// \f$C_{ab} - C_{a3}C_{3b}/C_{33}\f$.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/material/IsotropicMaterial.hpp"

namespace sparlab {

/// Internal variables of one integration point.
struct PlasticState {
  Vector6 plastic_strain = Vector6::Zero();  ///< \f$\varepsilon^p\f$, engineering shears [-]
  Vector6 back_stress = Vector6::Zero();     ///< \f$\beta\f$ [Pa]
  Scalar equivalent_plastic_strain = 0.0;    ///< \f$\bar\alpha\f$ [-]
  Scalar thickness_strain = 0.0;             ///< plane stress: total \f$\varepsilon_{33}\f$ [-]
  /// The return that produced this state was plastic. A point still on its
  /// yield surface is then assumed to go on loading: its tangent at zero
  /// strain increment is the continuum elastoplastic one (\f$\Delta\gamma = 0\f$
  /// in the consistent tangent), not the elastic one - which only changes the
  /// predictor of the next step, never a converged state.
  bool loading = false;
};

/// The result of one return.
struct PlasticResponse {
  Vector6 stress = Vector6::Zero();    ///< \f$\sigma\f$ [Pa]
  Matrix6 tangent = Matrix6::Zero();   ///< \f$d\sigma/d\varepsilon\f$ [Pa]; condensed in plane stress
  PlasticState state;                  ///< the updated internal variables
  bool yielding = false;               ///< the return was plastic
  Scalar energy = 0.0;                 ///< stored energy density, elastic + hardening [J/m^3]
  Scalar strain_33 = 0.0;              ///< the total \f$\varepsilon_{33}\f$ it used [-]
};

/// Stress, consistent tangent and updated state at the total strain
/// `strain` (all six components; in plane stress component 2 is ignored and
/// found) from the converged state `committed`, with the thermal strain
/// \f$\alpha\,\Delta T\f$ of `delta_t`.
/// \throws SolverError when the plane-stress iteration or the Voce return
///         fails to converge.
PlasticResponse j2_return(const IsotropicMaterial& material, StressState state,
                          const Vector6& strain, const PlasticState& committed,
                          Scalar delta_t = 0.0);

/// The von Mises stress of a 3-D Voigt stress [Pa].
Scalar von_mises_stress(const Vector6& stress);

}  // namespace sparlab
