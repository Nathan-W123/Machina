/// \file Elastoplastic.hpp
/// \brief One element of an elastoplastic analysis: internal force,
///        consistent tangent and updated internal variables from the J2
///        return of Plasticity.hpp at each integration point, with small-
///        strain or finite (total Lagrangian) kinematics.
///
/// **Small strain.** With the linear strain operator \f$B\f$ (the element's
/// own, in the 3-D Voigt order - a plane model's out-of-plane row carries the
/// plane-strain or mean-dilatation \f$\varepsilon_{33}\f$), the strain at a
/// point is \f$\varepsilon = Bu_e\f$ and
/// \f[
///   f_e = \int B^T\sigma\,dV ,\qquad K_e = \int B^T C\,B\,dV ,
/// \f]
/// with \f$\sigma\f$ and the consistent tangent \f$C\f$ of the return from
/// the last converged state.
///
/// **Finite kinematics** (large displacement and rotation, small strain).
/// The return takes the Green-Lagrange strain \f$E\f$ for the strain and
/// gives the second Piola-Kirchhoff stress \f$S\f$: the additive split
/// \f$E = E^e + E^p + E^\theta\f$ with \f$S = D E^e\f$, the elastoplastic
/// counterpart of the Saint Venant-Kirchhoff law and, like it, exact under
/// any rigid rotation and meant for strains that stay small (with small
/// strains the materially-non-linear-only relations may be written between
/// \f$S\f$ and \f$E\f$ at large displacement and rotation: Bathe, *Finite
/// Element Procedures*, 1996, ch. 6). Then
/// \f[
///   f_e = \int B_{NL}^T S\,dV_0 ,\qquad
///   K_e = \int B_{NL}^T C\,B_{NL}\,dV_0 + \int (G^T\mathbf{S}\,G)\otimes I\,dV_0
/// \f]
/// with the operators of TotalLagrangian.hpp. The thermal strain
/// \f$E^\theta\f$ is the Green strain of the free thermal stretch, as in
/// Hyperelastic.hpp, so a freely heated body is stress-free, but the elastic
/// stiffness is not rescaled by the thermal stretch as the Saint
/// Venant-Kirchhoff law's is (a relative difference of \f$\alpha\Delta T\f$
/// in a thermal stress, below the model's small-strain premise).
///
/// **Mean dilatation.** Plastic flow is isochoric. On an element with
/// several integration points (Q4, Hex8, Tet10) each point's dilatation is
/// then a constraint, and in plane strain and 3-D the constraints outnumber
/// the element's degrees of freedom: the element locks and overestimates the
/// collapse load (Nagtegaal, Parks and Rice 1974). The mean-dilatation
/// operator (Hughes 1980) replaces every point's dilatation by the element's
/// volume average, \f$\bar B = B + \tfrac{1}{3}\,m\,(\bar b - b)^T\f$ with
/// \f$b = B^T m\f$ the dilatation row, \f$\bar b\f$ its average and
/// \f$m = \{1,1,1,0,0,0\}\f$, which leaves one constraint per element. In
/// plane strain the average dilatation gives \f$\varepsilon_{33} =
/// (\bar\theta - \theta)/3\f$ at a point, zero on average. With finite
/// kinematics the same is done to the Green strain,
/// \f$\bar E = E + \tfrac{1}{3}\,m\,(\overline{\mathrm{tr}E} - \mathrm{tr}E)\f$,
/// which is objective and whose second variation adds
/// \f$\tfrac{1}{3}\,\mathrm{tr}S\,(\overline{G^TG} - G^TG)\otimes I\f$ to the
/// geometric stiffness; being derived from the strain, it keeps the tangent
/// symmetric and consistent. Plane stress has no incompressibility
/// constraint (the thickness strain is free) and keeps \f$B\f$. On a
/// constant-strain element (Tri3, Tet4) the average is the point value and
/// nothing changes - those elements lock.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/fem/FemModel.hpp"
#include "sparlab/material/Plasticity.hpp"

#include <string>
#include <vector>

namespace sparlab {

/// Kinematics of a non-linear analysis.
enum class Kinematics {
  Finite,      ///< total Lagrangian: large displacement and rotation
  SmallStrain  ///< linear strain: small displacement and rotation
};

std::string to_string(Kinematics kinematics);
/// "finite" or "small_strain".
Kinematics parse_kinematics(const std::string& text);

struct ElastoplasticElement {
  Vector internal_force;               ///< \f$f_e\f$ [N]
  Matrix tangent;                      ///< \f$K_e\f$ [N/m], if requested
  /// \f$\partial f_e/\partial\lambda\f$ at fixed displacement for a
  /// temperature change \f$\lambda\Delta T_0\f$ [N]; empty without one.
  Vector thermal_force_rate;
  Scalar energy = 0.0;                 ///< stored energy, elastic + hardening [J]
  std::vector<PlasticState> states;    ///< per integration point, after the return
  int yielding_points = 0;             ///< points whose return was plastic
};

/// The element at displacement `ue` from the converged states `committed`
/// (one per point of the element's stiffness rule), with the nodal
/// temperature changes `temperature` (per unit load factor; may be null)
/// scaled by `temperature_scale`. `mean_dilatation` selects B-bar (ignored
/// in plane stress).
/// \throws SolverError when a return fails to converge.
ElastoplasticElement elastoplastic_element(const FemModel& model, Index e, const Vector& ue,
                                           const std::vector<PlasticState>& committed,
                                           bool mean_dilatation, const Vector* temperature,
                                           Scalar temperature_scale, bool want_tangent,
                                           Kinematics kinematics = Kinematics::SmallStrain);

/// Stress state of an element at a converged state, from the point returns:
/// the point averages of the Cauchy stress (small strain: \f$\sigma\f$) and
/// of the second Piola-Kirchhoff stress (small strain: \f$\sigma\f$ again),
/// 3-D Voigt, and measures of the deformation.
struct ElastoplasticStress {
  Vector6 cauchy = Vector6::Zero();           ///< [Pa]
  Vector6 piola_kirchhoff = Vector6::Zero();  ///< [Pa]
  Scalar von_mises = 0.0;                     ///< of the average Cauchy stress [Pa]
  Scalar max_point_von_mises = 0.0;           ///< the largest of its points [Pa]
  Scalar max_equivalent_plastic_strain = 0.0; ///< \f$\max\bar\alpha\f$ over its points [-]
  int plastic_points = 0;                     ///< points with \f$\bar\alpha > 0\f$
  /// The largest tensor strain component over its points: \f$\varepsilon\f$
  /// (small strain) or \f$E\f$ (finite, with the thickness strain) [-].
  Scalar max_strain = 0.0;
  /// Small strain: the largest infinitesimal rotation \f$|\mathrm{skew}\,H|\f$
  /// and the largest component of the neglected \f$\tfrac12 H^TH\f$ [-].
  Scalar max_rotation = 0.0;
  Scalar max_quadratic_strain = 0.0;
  Scalar min_jacobian = 1.0;                  ///< finite: the smallest det F [-]
};
ElastoplasticStress elastoplastic_stress(const FemModel& model, Index e, const Vector& ue,
                                         const std::vector<PlasticState>& committed,
                                         bool mean_dilatation, const Vector* temperature,
                                         Scalar temperature_scale,
                                         Kinematics kinematics = Kinematics::SmallStrain);

/// Number of integration points of the model's stiffness rule.
int elastoplastic_points(const FemModel& model);

}  // namespace sparlab
