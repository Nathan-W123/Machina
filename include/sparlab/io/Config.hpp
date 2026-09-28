/// \file Config.hpp
/// \brief Declarative JSON input deck: schema, parsing and model construction.
///
/// An input deck fully determines a run. The schema is documented in
/// docs/configuration.md; the section names below mirror it one to one:
/// \code
/// {
///   "name": "cantilever_beam",
///   "mesh":     { "type": "structured_quad", "nx":.., "ny":.., "lx":.., "ly":.. },
///               // or { "type": "structured_hex", "nx","ny","nz", "lx","ly","lz" }
///               // or "structured_tri" / "structured_tet" with the same keys
///               // or { "type": "file", "path": "part.msh", "scale": 0.001 }
///               // "order": 2 turns tetrahedra into 10-node Tet10 cells
///   "material": { "youngs_modulus":.., "poisson_ratio":.., "density":.. },
///   "model":    { "thickness":.., "stress_state": "plane_stress" },
///   "boundary_conditions": [ { "fix": ["x","y"], "region": {..} } ],
///   "load_cases":          [ { "name":.., "weight":.., "point_loads": [..] } ],
///   "solver":   { "linear": {..}, "equilibrium_tolerance":.. },
///   "modal":    { "enabled": true, "num_modes":.. },
///   "buckling": { "enabled": true, "num_modes":.., "load_cases": [..] },
///   "topology": { "enabled": true, "volume_fraction":.., "filter": {..} },
///   "output":   { "vtk": true, "csv": true }
/// }
/// \endcode
///
/// The mesh fixes the spatial dimension of the whole deck: a solid mesh
/// (`structured_hex`, `structured_tet`, or a file of tetrahedra / hexahedra)
/// takes three-entry vectors, may fix `z`, and uses the three-dimensional
/// constitutive law; a plane mesh (`structured_quad`, `structured_tri`, or a
/// file of triangles / quadrilaterals) takes two-entry vectors and a plane
/// stress state. Mixing the two is an error rather than a silently truncated
/// component. A mesh file is therefore read while the deck is parsed; its
/// `path` is relative to the deck's own directory, and its physical groups
/// (Gmsh) or *NSET / *ELSET cards (Abaqus) become regions through
/// `{"group": "<name>"}`.
///
/// Every numeric field is in SI units. Unknown keys are reported (not ignored),
/// because a misspelled key that silently takes its default is a direct route
/// to a wrong published number.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/fem/Buckling.hpp"
#include "sparlab/fem/FemModel.hpp"
#include "sparlab/fem/LinearSolver.hpp"
#include "sparlab/fem/ModalAnalysis.hpp"
#include "sparlab/fem/StaticAnalysis.hpp"
#include "sparlab/io/Json.hpp"
#include "sparlab/io/MeshReader.hpp"
#include "sparlab/mesh/StructuredMesh.hpp"
#include "sparlab/topopt/DensityFilter.hpp"
#include "sparlab/topopt/DesignDomain.hpp"
#include "sparlab/topopt/TopologyOptimizer.hpp"

#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace sparlab {

/// Mesh sources a deck can name.
enum class MeshKind {
  StructuredQuad,  ///< "structured_quad": 2-D Q4 grid
  StructuredHex,   ///< "structured_hex": 3-D Hex8 grid
  StructuredTri,   ///< "structured_tri": the Q4 grid, each cell split into 2 Tri3
  StructuredTet,   ///< "structured_tet": the Hex8 grid, each cell split into 6 Tet4
  File             ///< "file": an unstructured Gmsh (.msh) or Abaqus (.inp) mesh
};

std::string to_string(MeshKind kind);

/// True for the generated box meshes, whose resolution a deck (or `--nx`)
/// sets; false for a mesh read from a file, whose resolution the mesher set.
bool is_structured(MeshKind kind);

/// A mesh read from a file.
struct MeshFileConfig {
  std::string path;           ///< as written in the deck
  std::string resolved_path;  ///< resolved against the deck's directory
  MeshReadOptions read;
};

struct ModalConfig {
  bool enabled = false;
  ModalAnalysisOptions options;
  /// Also run modal analysis on the density-thresholded topology (topopt runs).
  bool analyse_optimised_topology = true;
  /// Compare against a uniformly thinned domain of the same mass (2-D only:
  /// the argument rests on thickness scaling and has no 3-D counterpart).
  bool compare_mass_matched_baseline = true;
};

/// Linear buckling check of the load cases (sparlab_solve), and of the full
/// solid domain and the interpreted structure after an optimisation
/// (sparlab_topopt).
struct BucklingConfig {
  bool enabled = false;
  BucklingOptions options;
  /// Load cases to check, by name; empty checks every case.
  std::vector<std::string> load_cases;
  /// Also check the thresholded structure of a topology run.
  bool analyse_optimised_topology = true;
};

struct TopologyConfig {
  bool enabled = false;
  Scalar volume_fraction = 0.4;
  /// Starting density for free variables; negative means "use volume_fraction".
  Scalar initial_density = -1.0;
  FilterType filter_type = FilterType::Density;
  /// Filter radius in metres. When `filter_radius_elements` > 0 the radius is
  /// computed as that multiple of the mean element size instead.
  Scalar filter_radius = 0.0;
  Scalar filter_radius_elements = 1.5;
  TopologyOptimizerOptions optimizer;
  std::vector<PassiveRegionSpec> passive_regions;
  /// A `topology.overhang` section names a build direction: the final design
  /// is checked for unsupported material, and with `filter` the overhang
  /// filter shapes the optimisation (optimizer.overhang holds the settings).
  bool overhang_check = false;
  /// Minimum length-scale check of the final design: probe radii from half a
  /// cell size up to `length_scale_max_radius_elements` cells in half-cell
  /// steps, a probe passing while it flags at most `length_scale_tolerance`
  /// of the volume. On by default with the robust formulation.
  bool length_scale_check = false;
  Scalar length_scale_max_radius_elements = 6.0;
  Scalar length_scale_tolerance = 0.02;
};

/// A material assigned to an element region, overriding the deck's primary
/// `material` there (later regions win).
struct MaterialRegion {
  std::string name;
  SelectorGroup region;
  IsotropicMaterial material;
};

struct OutputConfig {
  bool write_csv = true;
  bool write_vtk = true;
  bool write_density_history = true;
  bool write_mode_shapes = true;
};

/// A complete, validated run specification.
class Configuration {
 public:
  std::string name = "case";
  std::string description;
  std::string source_path;

  MeshKind mesh_kind = MeshKind::StructuredQuad;
  StructuredMeshSpec mesh_spec;   ///< structured kinds
  /// Polynomial order of the tetrahedra: 1 (Tet4) or 2 (Tet10). A
  /// `structured_tet` deck with order 2 splits the grid into Tet10 cells; a
  /// file of linear tetrahedra with order 2 is elevated to straight-sided
  /// Tet10 cells (`mesh_elevated`); a file of 10-node tetrahedra is order 2
  /// by itself.
  int mesh_order = 1;
  bool mesh_elevated = false;
  MeshFileConfig mesh_file;       ///< MeshKind::File
  /// The mesh read from `mesh_file` while parsing (shared by copies of the
  /// configuration, which the thinned-plate comparison makes), and what the
  /// reader found.
  std::shared_ptr<const Mesh> file_mesh;
  MeshReadReport mesh_report;
  Scalar thickness = 1.0;
  StressState stress_state = StressState::PlaneStress;
  IntegrationOptions integration;

  std::vector<DisplacementConstraint> constraints;
  std::vector<LoadCaseSpec> load_cases;
  std::vector<MaterialRegion> material_regions;

  StaticAnalysisOptions analysis;
  ModalConfig modal;
  BucklingConfig buckling;
  TopologyConfig topology;
  OutputConfig output;

  /// The parsed document, retained so run summaries can record exactly what
  /// was requested.
  json::Value document;

  /// Spatial dimension implied by the mesh (2 or 3).
  int dim() const;

  /// One-line description of the mesh source for messages and summaries,
  /// e.g. "structured_tet 40 x 20 x 10" or "file 'bracket.msh' (Tri3)".
  std::string describe_mesh() const;

  const IsotropicMaterial& material() const;
  void set_material(IsotropicMaterial material);

  /// Filter radius resolved against a concrete mesh [m].
  Scalar resolved_filter_radius(const Mesh& mesh) const;

  /// Indices of the load cases the buckling check covers.
  std::vector<std::size_t> buckling_load_cases() const;

 private:
  std::optional<IsotropicMaterial> material_;
};

/// Parse a configuration from an already-parsed document.
/// \param source name used in error messages.
/// \param strict when true, unknown keys raise ConfigError instead of a warning.
/// \param base_directory directory a relative mesh-file path is resolved
///        against; empty means the current working directory.
Configuration parse_configuration(const json::Value& document, const std::string& source,
                                  bool strict = false,
                                  const std::string& base_directory = "");

/// Read and parse a configuration file; a relative mesh-file path is resolved
/// against the file's own directory.
Configuration load_configuration(const std::string& path, bool strict = false);

/// Build the mesh described by a configuration (without the model).
Mesh build_mesh(const Configuration& config);

/// Build the mesh, material, boundary conditions and load cases described by a
/// configuration and finalise the model.
FemModel build_model(const Configuration& config);

/// Build the design domain for a topology optimisation run.
DesignDomain build_design_domain(const Configuration& config, const FemModel& model);

/// Parse a region specification (used by tests and by the config parser).
/// \param dim mesh dimension (2 or 3), which decides how many entries the
///        vector-valued fields must carry; 0 accepts either form.
SelectorGroup parse_region(const ConfigNode& node, const std::string& default_name,
                           int dim = 0);

}  // namespace sparlab
