/// \file ResultWriter.hpp
/// \brief Writes every run artefact into one output directory.
///
/// Layout produced for a case named `<case>` in `<out>`:
/// \code
///   <out>/summary.json              scalars, tolerances, timings, provenance
///   <out>/config.json               verbatim echo of the input deck
///   <out>/mesh.json                 nodes, connectivity, BCs, applied loads
///   <out>/displacement_<lc>.csv     nodal displacements per load case
///   <out>/stress_<lc>.csv           element strains/stresses per load case
///   <out>/reactions_<lc>.csv        support reactions per load case
///   <out>/fields_<lc>.vtk           ParaView fields per load case
///   <out>/modes.csv                 eigenvalues, frequencies, residuals
///   <out>/mode_shapes.csv           M-orthonormal mode shapes
///   <out>/mode_<k>.vtk              one file per mode shape
///   <out>/buckling.csv              buckling load factors per checked load case
///   <out>/buckling_<lc>_<k>.vtk     buckling mode shapes (max |phi| = 1)
///   <out>/nonlinear_<lc>.csv        load-displacement path of a non-linear run
///   <out>/nonlinear_displacement_<lc>.csv  its final nodal displacements
///   <out>/nonlinear_stress_<lc>.csv its final Cauchy and 2nd Piola-Kirchhoff stresses
///   <out>/nonlinear_reactions_<lc>.csv     its final support reactions
///   <out>/nonlinear_<lc>.vtk        its final fields
///   <out>/transient_<lc>.csv        time history of a transient: energies, monitors
///   <out>/transient_state_<lc>.csv  its final displacement, velocity, acceleration
///   <out>/transient_reactions_<lc>.csv     its final support reactions
///   <out>/transient_<lc>_<k>.vtk    its snapshots, indexed with their times by
///                                   transient_<lc>.vtk.series (ParaView file series)
///   <out>/frequency_response_<lc>.csv      complex monitor amplitudes per frequency
///   <out>/frequency_response_<lc>_<k>.vtk  complex fields at the snapshot frequencies,
///                                   indexed by frequency_response_<lc>.vtk.series
///   <out>/history.csv               optimisation iteration history
///   <out>/density_final.csv         final design and physical density
///   <out>/density_history.csv       density snapshots for the animation
///   <out>/structure_before.{vtk,stl} the design domain the optimiser started from
///   <out>/structure_after.{vtk,stl}  the thresholded, largest-group structure
/// \endcode
/// Every file is plain text. `summary.json` is the single source of truth for
/// the numbers quoted in the documentation.
#pragma once

#include "sparlab/core/Timer.hpp"
#include "sparlab/core/Types.hpp"
#include "sparlab/fem/Buckling.hpp"
#include "sparlab/fem/Dynamics.hpp"
#include "sparlab/fem/ModalAnalysis.hpp"
#include "sparlab/fem/ModelDiagnostics.hpp"
#include "sparlab/fem/NonlinearStatic.hpp"
#include "sparlab/fem/StaticAnalysis.hpp"
#include "sparlab/fem/StressRecovery.hpp"
#include "sparlab/io/Config.hpp"
#include "sparlab/io/Json.hpp"
#include "sparlab/mesh/SubMesh.hpp"
#include "sparlab/topopt/LengthScale.hpp"
#include "sparlab/topopt/OverhangFilter.hpp"
#include "sparlab/topopt/TopologyOptimizer.hpp"

#include <string>
#include <vector>

namespace sparlab {

class ResultWriter {
 public:
  /// Create (or reuse) the output directory.
  ResultWriter(std::string directory, const Configuration& config);

  const std::string& directory() const { return directory_; }

  /// Echo the input deck so a result directory is self-describing.
  void write_config() const;

  /// Nodes, connectivity, prescribed DOFs and the applied load vectors.
  void write_mesh(const FemModel& model) const;

  /// Nodal displacements of one load case (`<stem>_<lc>.csv`).
  void write_displacement(const Mesh& mesh, const std::string& load_case,
                          const Vector& displacement,
                          const std::string& stem = "displacement") const;

  /// Element strain/stress data of one load case.
  void write_stress(const Mesh& mesh, const std::string& load_case,
                    const StressField& field, const Vector* density = nullptr) const;

  /// Support reactions of one load case (`<stem>_<lc>.csv`).
  void write_reactions(const Mesh& mesh, const DofManager& dofs,
                       const std::string& load_case, const Vector& reactions,
                       const std::string& stem = "reactions") const;

  /// Combined VTK file for one load case, with the nodal temperatures of a
  /// thermal case when given.
  void write_static_vtk(const Mesh& mesh, const std::string& load_case,
                        const Vector& displacement, const StressField& field,
                        const Vector* density, const Vector* stiffness_factor,
                        const Vector* temperature = nullptr) const;

  /// Nodal temperatures of one load case (`temperature_<lc>.csv`).
  void write_temperature(const Mesh& mesh, const std::string& load_case,
                         const Vector& temperature) const;

  /// Eigenvalues, frequencies, residuals and (optionally) mode shapes.
  void write_modal(const Mesh& mesh, const ModalResult& modal,
                   const std::string& tag = "") const;

  /// Buckling load factors of the checked load cases (`buckling<_tag>.csv`)
  /// and, when mode shapes are written, one VTK file per mode with the shape
  /// scaled to a largest nodal displacement of 1 (a buckling mode has no
  /// amplitude of its own).
  void write_buckling(const Mesh& mesh, const std::vector<BucklingResult>& results,
                      const std::string& tag = "") const;

  /// A non-linear run: the load-displacement path (`nonlinear_<lc>.csv`,
  /// always written) and, as the output settings ask, the final nodal
  /// displacements, element stresses and reactions as CSV and the final
  /// fields as VTK.
  void write_nonlinear(const FemModel& model, const NonlinearResult& result) const;

  /// A transient run: the history of monitors, energies and (non-linear)
  /// iterations (`transient_<lc>.csv`, always written); as the output
  /// settings ask, the final displacement, velocity and acceleration with
  /// the reactions as CSV, and the snapshots as a numbered VTK series
  /// (`transient_<lc>_0000.vtk`, ...) with its index of times.
  void write_transient(const FemModel& model, const TransientResult& result) const;

  /// A frequency response: the complex amplitudes of the monitors at every
  /// frequency (`frequency_response_<lc>.csv`, always written) and, as the
  /// output settings ask, the snapshot fields (real part, imaginary part and
  /// amplitude of the displacement) as VTK.
  void write_frequency_response(const FemModel& model,
                                const FrequencyResponseResult& result) const;

  /// Optimisation iteration history.
  void write_history(const TopologyOptimizationResult& result) const;

  /// Final design, physical density and per-element energy.
  void write_density(const Mesh& mesh, const DesignDomain& domain,
                     const TopologyOptimizationResult& result) const;

  /// Density snapshots used by the evolution animation.
  /// \param max_frames snapshots are strided down to at most this many frames.
  void write_density_history(const TopologyOptimizationResult& result,
                             int max_frames = 60) const;

  /// One geometry as a solid: `<stem>.vtk` (the cells, with `density` as cell
  /// data) and `<stem>.stl` (the watertight boundary surface, extruded by
  /// `thickness` for a 2-D mesh). Returns the surface statistics and the
  /// mismatch between the enclosed and the cell volume, which the summary
  /// records.
  /// \param density per-cell density written to the VTK file (length
  ///        num_elements of `mesh`).
  json::Value write_geometry(const Mesh& mesh, const Vector& density, Scalar thickness,
                             const std::string& stem, const std::string& what) const;

  /// Write an arbitrary JSON document into the directory.
  void write_json(const std::string& file_name, const json::Value& value) const;

  /// Path of a file inside the output directory.
  std::string file(const std::string& name) const;

 private:
  std::string directory_;
  const Configuration& config_;
};

/// Build the JSON summary of a static (plus optional modal) analysis.
json::Value make_static_summary(const Configuration& config, const FemModel& model,
                                const ModelDiagnostics& diagnostics,
                                const std::vector<StaticSolution>& solutions,
                                const std::vector<StressField>& stresses,
                                const ModalResult* modal, const TimingLedger& timings);

/// Build the JSON summary of a topology-optimisation run.
json::Value make_topology_summary(const Configuration& config, const FemModel& model,
                                  const DesignDomain& domain,
                                  const DensityFilter& filter,
                                  const TopologyOptimizationResult& result,
                                  const TopologyInterpretation* interpretation,
                                  const ModalResult* modal_initial,
                                  const ModalResult* modal_optimised,
                                  const ModalResult* modal_mass_matched,
                                  const TimingLedger& timings);

/// Common provenance block: version, build type, timestamp, tolerances.
json::Value make_provenance(const Configuration& config);

/// Summary block of a set of buckling checks: per load case the load
/// factors, residuals, solid-energy fractions and solver statistics, with the
/// method, its tolerances and what a load factor means.
json::Value buckling_json(const std::vector<BucklingResult>& results,
                          const BucklingOptions& options, const std::string& what);

/// Summary block of the non-linear runs: the formulation, the options and
/// per load case the outcome, the final state and its checks. `linear` holds
/// the linear solutions of all load cases (indexed like the deck's cases) so
/// each run is set beside its linear counterpart.
json::Value nonlinear_json(const std::vector<NonlinearResult>& results,
                           const NonlinearOptions& options, const FemModel& model,
                           const std::vector<StaticSolution>& linear);

/// The `transient` block of summary.json: the method, its options and one
/// entry per load case (energies, peaks of the monitors, the final state).
json::Value transient_json(const std::vector<TransientResult>& results,
                           const TransientOptions& options);

/// The `frequency_response` block of summary.json: the options and one entry
/// per load case (the peak of every monitor and where it lies).
json::Value frequency_response_json(const std::vector<FrequencyResponseResult>& results,
                                    const FrequencyResponseOptions& options);

/// Summary block of the overhang check of the final design.
json::Value overhang_json(const OverhangReport& report, bool filtered);

/// Summary block of the minimum length-scale scan of the final design.
json::Value length_scale_json(const LengthScaleScan& scan);

}  // namespace sparlab
