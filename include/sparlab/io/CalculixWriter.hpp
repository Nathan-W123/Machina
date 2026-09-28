/// \file CalculixWriter.hpp
/// \brief Export of a finalised model as CalculiX (Abaqus-style) input decks.
///
/// One `.inp` per load case is written, each a complete linear static job:
/// nodes (1-based), elements (`CPS4` / `CPS3` for plane stress, `CPE4` /
/// `CPE3` for plane strain, `C3D8` / `C3D4` / `C3D10` for solids, in the same
/// local node order SparLab uses), one isotropic `*ELASTIC` material and
/// `*SOLID SECTION` per material of the model (with the thickness on its data
/// line for plane elements), the prescribed DOFs as `*BOUNDARY` cards, and
/// `*NODE FILE, U` so the solver writes nodal displacements to its `.frd`
/// result file.
///
/// Point loads and tractions go out as their assembled nodal forces
/// (`*CLOAD`). Every other load goes out in CalculiX's own form, so that
/// CalculiX integrates it with its own code and the comparison tests SparLab's
/// load assembly as well as its stiffness and solve:
///   * pressure as `*DLOAD` face loads `P1`...`P6` on the same faces;
///   * self-weight as `GRAV`, with `*DENSITY` in the materials;
///   * body force densities as `BX`, `BY`, `BZ` on an element set;
///   * rotation as `CENTRIF` (omega^2, a point on the axis, its direction);
///   * a temperature field as `*EXPANSION` (with `ZERO` at the reference
///     temperature), `*INITIAL CONDITIONS, TYPE=TEMPERATURE` at the reference
///     temperature and the nodal `*TEMPERATURE`s of the case.
/// A conducted temperature field also gets `<stem>_<case>_conduction.inp`, a
/// steady `*HEAT TRANSFER` job - on `DC3D8` / `DC3D4` / `DC3D10` cells, and on
/// the plane elements themselves in 2-D (CalculiX 2.21 reads no integration
/// point for `DC2D4` / `DC2D3`) - with `*CONDUCTIVITY`, the prescribed temperatures on
/// DOF 11, `*DFLUX` surface fluxes `S1`...`S6` and body fluxes `BF`, `*FILM`
/// convection `F1`...`F6`, and `*NODE FILE, NT`, so CalculiX solves the
/// conduction problem itself.
///
/// The exported problem is the *same discrete problem* SparLab solves
/// (identical mesh, element type, integration order and materials), so
/// agreement is expected to solver precision for the solid elements, up to
/// the quadrature CalculiX applies to a distributed load on a curved Tet10
/// face or cell, and to the level at which CalculiX's plane elements - which
/// it expands into solid elements through the thickness - reproduce a plane
/// Q4 or Tri3. Both are measured, not assumed; see
/// python/scripts/cross_validate.py.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/fem/FemModel.hpp"

#include <string>
#include <vector>

namespace sparlab {

/// Write `<stem>_<load case>.inp` for every load case of `model`, followed by
/// `<stem>_<load case>_conduction.inp` for a case whose temperature is
/// conducted.
/// \return the paths written, in load-case order.
/// \throws IoError when a file cannot be written, or for a thermal case whose
///         materials have different reference temperatures (CalculiX measures
///         thermal strain from the initial nodal temperature).
std::vector<std::string> write_calculix_decks(const FemModel& model, const std::string& stem,
                                              const std::string& case_name);

/// CalculiX element keyword for the model's element type and stress state.
std::string calculix_element_type(const FemModel& model);

}  // namespace sparlab
