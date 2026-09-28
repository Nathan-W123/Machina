/// \file LoadCaseData.hpp
/// \brief What a finalised model knows about each load case: its load vector
///        split by origin, and its temperature field.
///
/// The parts are kept apart because the analyses treat them differently: a
/// nonlinear or transient analysis scales the mechanical and body loads with
/// the load factor or an amplitude but enters the temperature through the
/// constitutive law, not as an equivalent load; the CalculiX export writes
/// each part in its own native form so the reference code computes the body
/// and thermal loads itself; and the elastic strain energy of a thermal case
/// needs the thermal part separately.
#pragma once

#include "sparlab/core/Types.hpp"

#include <string>

namespace sparlab {

/// Scalars of a steady heat-conduction solve.
struct ConductionSummary {
  /// Heat entering through sources, fluxes and convection towards the
  /// ambient temperature [W] (convection counted at the solved temperature).
  Scalar applied_heat = 0.0;
  /// Net heat leaving through the prescribed-temperature nodes [W]: minus the
  /// sum of their "reactions" K T - F.
  Scalar prescribed_heat = 0.0;
  /// Gross heat entering and leaving through prescribed temperatures [W]
  /// (positive and negative reactions summed separately): the flow through a
  /// wall held at two temperatures, whose net is zero.
  Scalar prescribed_inflow = 0.0;
  Scalar prescribed_outflow = 0.0;
  /// |applied - prescribed| over the largest heat flow of the problem.
  Scalar relative_balance_error = 0.0;
  Scalar scaled_residual = 0.0;    ///< ||K T - Q|| / ||Q|| on the free nodes
  Scalar min_temperature = 0.0;    ///< [K]
  Scalar max_temperature = 0.0;    ///< [K]
  int prescribed_nodes = 0;
  int flux_faces = 0;
  int convection_faces = 0;
  int source_elements = 0;
  std::string solver;
};

/// One load case of a finalised model.
struct LoadCaseData {
  /// Point loads, tractions and pressures [N].
  Vector mechanical;
  /// Gravity, body forces and the centrifugal load [N]; empty when the case
  /// has none.
  Vector body;
  /// Equivalent nodal forces of the thermal strain,
  /// \f$\int B^T D\,\varepsilon_0\,dV\f$ [N]; empty without a temperature.
  Vector thermal;
  /// \f$\tfrac12\int \varepsilon_0^T D\,\varepsilon_0\,dV\f$ [J], so the
  /// elastic strain energy of a solution is
  /// \f$\tfrac12 u^T K u - u^T f_{th} + \tfrac12\int\varepsilon_0^T D\varepsilon_0\f$.
  Scalar thermal_self_energy = 0.0;
  /// Nodal temperatures [K]; empty when the case has no temperature field.
  Vector temperature;
  bool conduction_solved = false;
  ConductionSummary conduction;
  /// Total resultant of the body part (sum of its nodal forces) [N].
  Vector3 body_resultant = Vector3::Zero();
};

}  // namespace sparlab
