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
/// A non-linear analysis of a case also gets `<stem>_<case>_nlgeom.inp`
/// (finite kinematics) or `<stem>_<case>_small_strain.inp`: the same model
/// and loads in a `*STEP` (with `NLGEOM` for finite kinematics) whose
/// `*STATIC` step ramps every load (and prescribed displacement) with the
/// step time from 0 to 1, as SparLab's load factor does - one step per leg of
/// a load path, each ramping from the last leg's level to its own. An
/// elastic NLGEOM step takes automatically sized increments no longer than
/// 1 / `increments`; an elastoplastic one takes exactly `increments` fixed
/// increments per leg (`DIRECT`), the steps of SparLab's run, because the
/// backward-Euler return depends on the increments on a non-proportional
/// path. An elastoplastic material goes out with `*PLASTIC`: its isotropic
/// hardening curve as (yield stress, equivalent plastic strain) pairs, exact
/// for linear hardening (CalculiX holds the last value beyond the table, so
/// the table runs to a plastic strain of 10) and a piecewise-linear table
/// within 1e-4 of the saturation stress for Voce hardening. Kinematic
/// hardening is not exported: CalculiX 2.21's `HARDENING=KINEMATIC`, given
/// the table of a linear rule, softens a single element in uniaxial tension
/// at the rate the rule hardens it (measured), so such a case gets no
/// non-linear deck; nor does a Hill48 material (CalculiX's `*PLASTIC` is von
/// Mises) or one with Armstrong-Frederick backstresses. Under NLGEOM CalculiX's
/// `*ELASTIC` material is Saint Venant-Kirchhoff (second Piola-Kirchhoff
/// stress linear in the Green-Lagrange strain) and its `*DLOAD` pressure
/// follows the deforming face, so a follower pressure stays a face load
/// while a dead one goes out as its nodal forces on the undeformed faces
/// (`*CLOAD`). SparLab's neo-Hookean law has no CalculiX counterpart (its
/// `NEO HOOKE` splits the energy into an isochoric part and (J - 1)^2), so
/// such a run is not exported.
///
/// A transient run of a case also gets `<stem>_<case>_dynamic.inp`: the model
/// with `*DENSITY` and, for Rayleigh damping, `*DAMPING, ALPHA=a, BETA=b` in
/// every material (CalculiX's C = a M + b K of a direct integration), the
/// amplitude sampled at every step time as an `*AMPLITUDE` table (the method
/// reads the loads at the step times only, so the table is exact there),
/// and one `*STEP` (with `NLGEOM` for finite kinematics) of `*DYNAMIC, DIRECT,
/// ALPHA=alpha` - CalculiX's HHT-alpha method, whose sign convention and
/// Newmark parameters are SparLab's (measured: the two agree to the output
/// rounding with and without damping at alpha = 0 and -0.1) - with every
/// load and prescribed displacement on that amplitude and `*NODE FILE` at the
/// snapshot increments. CalculiX's implicit dynamics uses the consistent
/// mass and starts at rest, so a lumped-mass or preloaded run is not
/// exported; nor is a thermal or rotating load case.
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
#include "sparlab/fem/Dynamics.hpp"
#include "sparlab/fem/FemModel.hpp"

#include <string>
#include <vector>

namespace sparlab {

/// The non-linear decks to write beside the linear ones.
struct CalculixNonlinearExport {
  std::vector<std::size_t> load_cases;  ///< indices of the cases to export
  /// Increments per leg of the load path: the largest increment of an
  /// elastic NLGEOM step, the fixed increment of an elastoplastic one.
  int increments = 10;
  bool follower_pressure = true;        ///< as NonlinearOptions::follower_pressure
  /// Finite kinematics (`*STEP, NLGEOM`); false for small strain.
  bool nlgeom = true;
  /// The load factors of the path's turning points, one `*STEP` each; empty
  /// for a single step from 0 to 1.
  std::vector<Scalar> load_path;
};

/// The transient decks to write beside the linear ones (`*DYNAMIC, DIRECT`).
struct CalculixTransientExport {
  std::vector<std::size_t> load_cases;  ///< indices of the cases to export
  TransientOptions options;             ///< the run's options (step, alpha, damping, ...)
};

/// Why the transient options cannot go out as a CalculiX `*DYNAMIC` deck of
/// load case `l` (lumped mass, a preloaded start, a thermal or rotating
/// load), or an empty string when they can.
std::string calculix_transient_obstacle(const FemModel& model, std::size_t l,
                                        const TransientOptions& options);

/// Why a material's plasticity cannot go out as CalculiX's `*PLASTIC` in a
/// non-linear or transient deck - Hill48 anisotropy, Armstrong-Frederick
/// backstresses with recovery, or linear kinematic hardening (Prager's
/// modulus or backstresses without recovery) - or an empty string when it
/// can (elastic, or von Mises with isotropic hardening).
std::string calculix_plasticity_obstacle(const IsotropicMaterial& material);

/// Write `<stem>_<load case>.inp` for every load case of `model`, followed by
/// `<stem>_<load case>_conduction.inp` for a case whose temperature is
/// conducted, `<stem>_<load case>_nlgeom.inp` for each case `nonlinear`
/// names and `<stem>_<load case>_dynamic.inp` for each case `transient`
/// names.
/// \return the paths written, in load-case order.
/// \throws IoError when a file cannot be written, for a thermal case whose
///         materials have different reference temperatures (CalculiX measures
///         thermal strain from the initial nodal temperature), for a
///         transient that calculix_transient_obstacle refuses, or for a
///         non-linear case with a material that calculix_plasticity_obstacle
///         refuses.
std::vector<std::string> write_calculix_decks(const FemModel& model, const std::string& stem,
                                              const std::string& case_name,
                                              const CalculixNonlinearExport* nonlinear = nullptr,
                                              const CalculixTransientExport* transient = nullptr);

/// CalculiX element keyword for the model's element type and stress state.
std::string calculix_element_type(const FemModel& model);

}  // namespace sparlab
