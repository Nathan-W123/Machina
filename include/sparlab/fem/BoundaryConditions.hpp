/// \file BoundaryConditions.hpp
/// \brief Displacement constraint and load specifications, and their
///        application to a DofManager / global force vector.
///
/// Sign convention: forces and displacements are positive along the global
/// +x / +y / +z axes. Tractions are given as a stress vector \f$\bar{t}\f$ [Pa]
/// in global components and are integrated over the selected boundary faces
/// (edges in 2-D), producing consistent nodal forces (not lumped ones).
///
/// Vectors are three-component throughout; on a 2-D model the z entry of a
/// force or traction must be zero and `fix_z` is an error, so a deck cannot
/// silently lose a component it thought it had applied.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/elements/Element.hpp"
#include "sparlab/fem/DofManager.hpp"
#include "sparlab/fem/Selector.hpp"

#include <string>
#include <vector>

namespace sparlab {

/// Prescribed displacement on the x, y and/or z component of a node region,
/// and, on a shell or beam model, prescribed rotations about the global axes.
struct DisplacementConstraint {
  SelectorGroup region;
  bool fix_x = false;
  bool fix_y = false;
  bool fix_z = false;
  bool fix_rx = false;   ///< rotation about x (shell and beam nodes only)
  bool fix_ry = false;
  bool fix_rz = false;
  Scalar value_x = 0.0;  ///< prescribed u_x [m]
  Scalar value_y = 0.0;  ///< prescribed u_y [m]
  Scalar value_z = 0.0;  ///< prescribed u_z [m]
  Scalar value_rx = 0.0;  ///< prescribed rotation about x [rad]
  Scalar value_ry = 0.0;
  Scalar value_rz = 0.0;

  /// Component `k` in node-major DOF order: 0-2 the translations, 3-5 the
  /// rotations.
  bool fixes(int component) const {
    switch (component) {
      case 0: return fix_x;
      case 1: return fix_y;
      case 2: return fix_z;
      case 3: return fix_rx;
      case 4: return fix_ry;
      case 5: return fix_rz;
      default: return false;
    }
  }
  Scalar value(int component) const {
    switch (component) {
      case 0: return value_x;
      case 1: return value_y;
      case 2: return value_z;
      case 3: return value_rx;
      case 4: return value_ry;
      case 5: return value_rz;
      default: return 0.0;
    }
  }
  void set(int component, bool fix, Scalar value = 0.0);
  bool fixes_rotation() const { return fix_rx || fix_ry || fix_rz; }
};

/// Concentrated nodal force (and, on a shell or beam model, moment) applied
/// to a node region.
struct PointLoadSpec {
  SelectorGroup region;
  Vector3 force = Vector3::Zero();  ///< [N]; z must be 0 on a 2-D model
  /// Nodal moment about the global axes [N m]; only a model whose nodes carry
  /// rotations (shells, beams) can take one.
  Vector3 moment = Vector3::Zero();
  /// When true, `force` is the *resultant* and is divided equally among the
  /// selected nodes (mesh-independent total load). When false, `force` is
  /// applied to each selected node.
  bool distribute_total = true;
};

/// Constant traction applied to the mesh boundary faces (edges in 2-D) whose
/// nodes all lie inside the region.
struct TractionLoadSpec {
  SelectorGroup region;
  Vector3 traction = Vector3::Zero();  ///< [Pa]; z must be 0 on a 2-D model
};

/// Uniform normal pressure on the boundary faces (edges in 2-D) whose nodes
/// all lie inside the region, positive pushing into the body. It is
/// integrated with the face's own area vector, so it acts normal to a curved
/// face at every point (FaceGeometry.hpp).
struct PressureLoadSpec {
  SelectorGroup region;
  Scalar pressure = 0.0;  ///< [Pa]
};

/// Uniform body force density on the elements of a region, or of the whole
/// model when `whole_model` is set.
struct BodyForceSpec {
  SelectorGroup region;
  bool whole_model = true;
  Vector3 force_density = Vector3::Zero();  ///< [N/m^3]; z must be 0 on a 2-D model
};

/// Steady rotation about an axis: every element carries the centrifugal body
/// force \f$\rho\,\omega^2 r_\perp\f$, with \f$r_\perp\f$ the distance vector
/// from the axis.
struct CentrifugalSpec {
  bool enabled = false;
  Scalar angular_velocity = 0.0;          ///< omega [rad/s]
  Vector3 axis = Vector3::UnitZ();        ///< unit direction of the axis
  Vector3 point = Vector3::Zero();        ///< a point on the axis [m]
};

/// A prescribed nodal value on a region: a temperature [K] on nodes, or a heat
/// input on faces or elements, depending on where it is used.
struct RegionValue {
  SelectorGroup region;
  bool whole_model = false;  ///< elements only: every element
  Scalar value = 0.0;
};

/// Convection from boundary faces to an ambient temperature,
/// \f$q = h\,(T - T_\infty)\f$ leaving the body.
struct ConvectionSpec {
  SelectorGroup region;
  Scalar film_coefficient = 0.0;  ///< h [W/(m^2 K)]
  Scalar ambient = 0.0;           ///< T_inf [K]
};

/// Steady heat conduction \f$-\nabla\cdot(k\nabla T) = Q\f$ whose solution is
/// the temperature field of a load case (HeatConduction.hpp).
struct ConductionSpec {
  std::vector<RegionValue> prescribed;   ///< temperatures on node regions [K]
  std::vector<RegionValue> fluxes;       ///< heat flux into the body on faces [W/m^2]
  std::vector<ConvectionSpec> convection;
  std::vector<RegionValue> sources;      ///< volumetric heat generation [W/m^3]
};

/// The temperature field of a load case.
struct TemperatureSpec {
  enum class Source {
    None,        ///< no temperature change: no thermal strain
    Uniform,     ///< every node at `uniform`
    Regions,     ///< `uniform` everywhere, then each region's value (later wins)
    Conduction   ///< the solution of `conduction`
  };
  Source source = Source::None;
  Scalar uniform = 0.0;              ///< [K]
  std::vector<RegionValue> regions;  ///< node regions and their temperatures [K]
  ConductionSpec conduction;
};

/// One load case. Multiple cases are combined by the topology optimizer with
/// the given weights; the FE solver solves each independently.
struct LoadCaseSpec {
  std::string name = "load_case";
  Scalar weight = 1.0;  ///< weight in the multi-load compliance objective [-]
  std::vector<PointLoadSpec> point_loads;
  std::vector<TractionLoadSpec> tractions;
  std::vector<PressureLoadSpec> pressures;
  /// Gravitational (or any uniform) acceleration of the whole model [m/s^2]:
  /// the body force \f$\rho\,g\f$ of every element, its self-weight.
  Vector3 gravity = Vector3::Zero();
  std::vector<BodyForceSpec> body_forces;
  CentrifugalSpec centrifugal;
  TemperatureSpec temperature;
  /// Set when the case carries no applied force and is driven purely by the
  /// prescribed displacements of the model (enforced-deflection studies, the
  /// constant-strain patch test). Declaring it explicitly keeps an accidentally
  /// empty load case an error rather than a silently trivial solution.
  bool prescribed_displacement_only = false;

  /// Gravity, body forces or a rotation: loads that act on the volume.
  bool has_body_loads() const {
    return gravity.squaredNorm() > 0.0 || !body_forces.empty() || centrifugal.enabled;
  }
  bool has_temperature() const { return temperature.source != TemperatureSpec::Source::None; }
  /// Any applied loading at all, of any kind.
  bool has_loads() const {
    return !point_loads.empty() || !tractions.empty() || !pressures.empty() ||
           has_body_loads() || has_temperature();
  }
};

/// The boundary faces (edges in 2-D) of `boundary` whose nodes all lie in
/// `region`: the faces a traction, pressure, heat flux or convection region
/// acts on, in the order of `boundary`.
std::vector<Mesh::BoundaryFace> faces_in_region(const Mesh& mesh,
                                                const std::vector<Mesh::BoundaryFace>& boundary,
                                                const SelectorGroup& region);

/// Apply every constraint to `dofs`.
/// \return the number of distinct DOFs constrained.
/// \throws ConfigError when a region selects no nodes (almost always an
///         authoring mistake, and silently ignoring it produces a singular
///         system later), when a 2-D model is asked to fix z, or when a
///         model without rotational DOFs is asked to fix a rotation.
Index apply_constraints(const Mesh& mesh,
                        const std::vector<DisplacementConstraint>& constraints,
                        DofManager& dofs);

/// Assemble the global force vector [N] of one load case.
/// \throws ConfigError for empty regions, a traction region that matches no
///         boundary face, or an out-of-plane component on a 2-D model.
Vector assemble_load_vector(const Mesh& mesh, const Element& element,
                            const LoadCaseSpec& load_case, Scalar thickness,
                            const IntegrationOptions& integration);

}  // namespace sparlab
