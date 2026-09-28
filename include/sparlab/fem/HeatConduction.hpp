/// \file HeatConduction.hpp
/// \brief Steady-state heat conduction: the temperature field of a thermal
///        load case.
///
/// The temperature \f$T\f$ solves
/// \f[
///   -\nabla\cdot(k\,\nabla T) = Q \ \text{in}\ \Omega, \qquad
///   T = \bar T \ \text{on}\ \Gamma_T, \qquad
///   k\,\partial_n T = q \ \text{on}\ \Gamma_q, \qquad
///   -k\,\partial_n T = h\,(T - T_\infty) \ \text{on}\ \Gamma_h ,
/// \f]
/// discretised with the structural mesh's own shape functions,
/// \f[
///   \Big(\int_\Omega k\,\nabla N^T\nabla N\,t\,dV + \int_{\Gamma_h} h\,N N^T t\,d\Gamma\Big) T
///   = \int_\Omega N\,Q\,t\,dV + \int_{\Gamma_q} N\,q\,t\,d\Gamma
///     + \int_{\Gamma_h} h\,T_\infty N\,t\,d\Gamma ,
/// \f]
/// where \f$t\f$ is the thickness of a plane model (a plane model conducts
/// in its plane, insulated on its faces) and 1 for a solid. The prescribed
/// temperatures are imposed by partitioning, like prescribed displacements.
/// The conductivity is each element's material's.
///
/// **Heat balance.** The residual \f$r = K T - F\f$ at the prescribed nodes is
/// the heat entering there; at steady state sources, fluxes, net convective
/// inflow and that prescribed inflow sum to zero. The solve reports both
/// sides and their mismatch, the thermal counterpart of the static solver's
/// force balance.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/fem/FemModel.hpp"
#include "sparlab/fem/LinearSolver.hpp"

namespace sparlab {

struct ConductionResult {
  Vector temperature;         ///< nodal temperatures [K]
  ConductionSummary summary;  ///< heat balance and solver statistics
};

/// Solve steady conduction on `model`'s mesh.
/// \throws ConfigError for a region that selects nothing, a material without
///         conductivity, or a problem with neither a prescribed temperature
///         nor convection (the temperature would be determined only up to a
///         constant); SolverError when the solve misses its residual
///         tolerance or the heat balance.
ConductionResult solve_conduction(const FemModel& model, const ConductionSpec& spec,
                                  const LinearSolverOptions& linear);

/// Element conductivity matrix \f$\int k\,\nabla N^T \nabla N\,t\,dV\f$ of
/// element `e` (num_nodes x num_nodes) [W/K].
Matrix element_conductivity(const FemModel& model, Index e);

}  // namespace sparlab
