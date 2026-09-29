# Architecture

## Layering

Each box is one directory. An arrow means "depends on"; there are no cycles and
no upward dependencies, so any layer can be replaced without touching the ones
below it.

```
                          +-------------------------+
                          |  apps/                  |
                          |  sparlab_solve          |   command-line drivers
                          |  sparlab_topopt         |   argument parsing,
                          |  sparlab_verify         |   orchestration, console
                          |  sparlab_bench          |   reports
                          +------------+------------+
                                       |
        +------------------------------+------------------------------+
        |                              |                              |
        v                              v                              v
+---------------+          +-------------------------+       +-----------------+
|  io/          |          |  topopt/                |       |  fem/           |
|  Json         |          |  DesignDomain           |       |  ModalAnalysis  |
|  Config       |--------->|  SimpInterpolation      |------>|  Buckling       |
|  MeshReader   |          |  DensityFilter          |       |  StaticAnalysis |
|  CsvWriter    |          |  OverhangFilter         |       |  StressRecovery |
|  VtkWriter    |          |  Projection             |       |  ModelDiagnostics|
|  StlWriter    |          |  Sensitivity            |       |  LinearSolver   |
|  CalculixWriter|         |  StressConstraint       |       |  Multigrid      |
|  ResultWriter |          |  BucklingConstraint     |       |  Assembler      |
+---------------+          |  LengthScale            |       |  FemModel       |
                           |  OptimalityCriteria     |       |  BoundaryConds  |
                           |  Mma                    |       |  DofManager     |
                           |  TopologyOptimizer      |       |  Selector       |
                           +-------------------------+       |  Loads          |
                                                             |  HeatConduction |
                                                             |  NonlinearStatic|
                                                             |  TotalLagrangian|
                                                             |  Elastoplastic  |
                                                             +--------+--------+
                                                                      |
                                    +---------------------------------+
                                    |                 |               |
                                    v                 v               v
                           +----------------+  +-------------+  +------------+
                           |  elements/     |  |  material/  |  |  mesh/     |
                           |  Element (abc) |  |  Isotropic  |  |  Mesh      |
                           |  Quad4  Tri3   |  |  Material   |  |  Structured|
                           |  Hex8   Tet4   |  |  Hyper-     |  |  SubMesh   |
                           |  Tet10         |  |  elastic    |  |            |
                           |  Quadrature    |  |  Plasticity |  |            |
                           |  FaceGeometry  |  |             |  |            |
                           +-------+--------+  +------+------+  +-----+------+
                                   |                  |               |
                                   +--------+---------+---------------+
                                            v
                                   +--------------------+
                                   |  core/             |
                                   |  Types  Exceptions |
                                   |  Logging  Timer    |
                                   |  Version           |
                                   +--------------------+

   python/sparlab_viz/  reads only the files io/ writes:
       loaders -> style -> fields / solid -> plots / plots3d / studies
   python/scripts/cross_validate.py  drives CalculiX and scikit-fem on the
       exported decks and compares nodal displacements, temperatures,
       buckling load factors and the final states of the non-linear runs
       (CalculiX NLGEOM and *PLASTIC; its own total Lagrangian and J2
       solvers in scikit-fem)
   python/scripts/make_meshes.py  generates the Gmsh meshes the
       real-geometry decks read (the meshes are committed)
   python/scripts/tet10_part_study.py  meshes the engine mount at several
       sizes and orders and compares Tet4 with Tet10
```

The core is dimension-generic at run time: `Mesh::dim()` is 2 or 3, the
element kernels are written for their own dimension (`Quad4` and `Tri3` on a
2-D, `Hex8`, `Tet4` and `Tet10` on a 3-D mesh) and everything above the element
layer - DOF
numbering, assembly, partitioning, stress recovery, modal analysis, the
optimiser, the writers - is written against `dim`, `nodes_per_element` and
the Voigt length rather than against a fixed 2. The Q4 kernels keep
fixed-size internal matrices, so the plane path produces bit-for-bit what it
did before the solid path existed - checked by re-running every benchmark
deck and comparing each summary, CSV and VTK file with the earlier output.

## What each component owns

| Component | Owns | Deliberately does *not* |
|-----------|------|--------------------------|
| `core/` | scalar and matrix aliases, the exception hierarchy, the logger, timing, build provenance | anything numerical |
| `mesh/Mesh` | nodal coordinates (2-D or 3-D), flat connectivity with an explicit stride, validation, boundary edge / face extraction, named node and element sets, cell-quality statistics | materials, DOFs, physics |
| `mesh/StructuredMesh` | the rectangular Q4 and the box Hex8 generators, the Tri3 and Tet4 meshes split from them and the Tet10 meshes elevated from those, each with a node-perturbing variant for patch tests; `elevate_to_tet10` for any Tet4 mesh | reading files |
| `mesh/SubMesh` | extracting an element subset with compact renumbering; edge- (face-) connected components; the density-to-solid interpretation | deciding *whether* to interpret |
| `material/IsotropicMaterial` | `(E, nu, rho)`, the plane-stress, plane-strain and three-dimensional matrices, the thermal properties, the plasticity parameters and the hardening curve, input validation | element integration |
| `material/Hyperelastic` | the Saint Venant-Kirchhoff and compressible neo-Hookean laws at a displacement gradient: energy, second Piola-Kirchhoff stress, tangent, the multiplicative thermal split, the thickness stretch of plane stress, the Cauchy stress | the element loop |
| `material/Plasticity` | the J2 return at one point from its converged internal variables - radial return, Voce Newton, the consistent and continuum tangents, the plane-stress thickness strain, the stored energy | where the strain comes from, and the history |
| `elements/Element` | the abstract kernel interface: stiffness, consistent mass, strain operator, edge or face traction, local face tables, and - from the element's integration rule - the geometric stiffness and its derivative with respect to the element displacements | the element's own geometry |
| `elements/Quad4`, `elements/Hex8` | the bilinear quadrilateral and the trilinear hexahedron | any other topology |
| `elements/Tri3`, `elements/Tet4` | the constant-strain triangle and tetrahedron, in closed form | quadrature |
| `elements/Tet10` | the isoparametric ten-node tetrahedron: shape functions, 4-point stiffness, collapsed-Gauss mass and face loads, volume and Jacobian ratio of curved cells | lumping (the assembler's HRZ) |
| `elements/Quadrature` | Gauss-Legendre rules on the line, the square and the cube; the symmetric 4-point tetrahedron rule and collapsed Gauss rules on the triangle and the tetrahedron | where they are used |
| `elements/FaceGeometry` | the shape functions and area vectors of edges and (curved) faces, the consistent pressure load on a face and its derivative with respect to the face's nodes (the follower-pressure stiffness), face integrals for fluxes and convection | which faces are loaded |
| `fem/DofManager` | DOF numbering, prescribed values, the free/prescribed partition, gather and scatter | assembly |
| `fem/Selector` | region selection by geometry (box, circle, annulus, sphere, ids, nearest node) or by a mesh file's named group, unions and complements | what a region is *for* |
| `fem/BoundaryConditions` | constraint and load specifications, and turning them into prescribed DOFs and a global force vector | solving |
| `fem/FemModel` | the complete discrete model: mesh, material, idealisation, integration orders, DOFs, load cases | any numerics |
| `fem/Assembler` | sparse assembly of `K` and `M`, block extraction, the uniform-mesh element cache, the cached sparsity pattern later assemblies scatter into | choosing scale factors |
| `fem/LinearSolver` | the solver backends behind one interface (direct, multigrid CG, Jacobi CG, the automatic choice), warm starts, pivot inspection, residual verification | model semantics |
| `fem/Multigrid` | the smoothed-aggregation hierarchy, its reuse across refactorisations, the V-cycle, preconditioned CG, the deterministic parallel kernels | which solver a model uses |
| `fem/ModelDiagnostics` | per-component rigid-body and floating-region detection before any factorisation | fixing the model |
| `fem/StaticAnalysis` | the reduced solve, reaction recovery, global equilibrium checks, one factorisation shared across load cases | stress |
| `fem/StressRecovery` | strain, stress, von Mises, principal stresses, element strain energy, nodal averaging | plotting |
| `fem/ModalAnalysis` | the generalised eigenproblem, subspace iteration, validity screening, the analytical references | design variables |
| `fem/Buckling` | the geometric stiffness assembly, the buckling eigenproblem by subspace iteration with the spectral transformation and inertia-placed shift, warm starts, the solid-energy diagnostic, the Euler and Engesser references | the constraint built on it |
| `fem/Loads` | the body loads from the consistent mass (self-weight, force densities, rotation), region temperatures, the thermal load and thermal strain | the conduction solve |
| `fem/HeatConduction` | steady conduction on the structural mesh: conductivity, sources, fluxes, convection, prescribed temperatures, the heat balance | the thermal strain it causes |
| `fem/TotalLagrangian` | one element's internal force, tangent (material and initial-stress parts), energy and thermal load rate at a displacement, and its Cauchy stress, largest Green strain and smallest `J`; the reference gradients and the Green-strain operator it shares | the global system |
| `fem/Elastoplastic` | one elastoplastic element with small-strain or finite kinematics: the point strains (with the mean dilatation), the returns, internal force, tangent (with the geometric part), thermal load rate and the updated internal variables; its stresses and deformation measures | keeping the history |
| `fem/NonlinearStatic` | the non-linear system of a load case (follower pressure, the deformed-position centrifugal load, dead loads, temperature, prescribed displacements, all scaled by one load factor), the committed internal variables of every elastoplastic point and their commit on convergence, Newton with the energy line search, load control along a load path with its critical-point, unreached-load and collapse bracketing, the arc-length method, inertia, monitors, the deformed force balance, the validity warnings | the material laws |
| `topopt/DesignDomain` | design variables, bounds, passive tags, element volumes, feasibility checks | the objective |
| `topopt/SimpInterpolation` | `E(rho)`, its derivative, the mass laws | assembly |
| `topopt/DensityFilter` | the filter operator, its exact adjoint, the Sigmund variant | the objective |
| `topopt/OverhangFilter` | the layer-by-layer AM filter on a structured grid, its adjoint recursion, the overhang check | the build plate's physics |
| `topopt/Projection` | the smoothed Heaviside projection, its derivative, the `beta` schedule's and the robust formulation's options | when to raise `beta` |
| `topopt/Sensitivity` | the weighted-compliance objective, the physical density (filtered, overhang-filtered and projected, at any threshold), the exact gradients through the whole chain, the linear solver an adjoint solve reuses, and the finite-difference verifier | the update rule |
| `topopt/StressConstraint` | the relaxed, p-norm-aggregated von Mises constraint per load case, its adaptive scale and its adjoint gradient through the filter | the optimiser |
| `topopt/BucklingConstraint` | the KS-aggregated lower bound on the load factors per load case, the stress interpolation without the stiffness floor, the adjoint sensitivity, the warm-started eigensolves | the optimiser |
| `topopt/LengthScale` | the morphological opening and closing of a thresholded design and the scan of probe radii | enforcing a length scale |
| `topopt/OptimalityCriteria` | one OC step with multiplier bisection, against a fixed or a moving (robust) volume target | the loop |
| `topopt/Mma` | one MMA step: asymptotes, the convex subproblem and its primal-dual interior-point solve | the objective and constraints |
| `topopt/TopologyOptimizer` | the loop for either method, penalty and `beta` continuation, the robust formulation's three designs and volume rescaling, the erosion check, convergence, history, snapshots | I/O |
| `io/Json` | a self-contained JSON reader/writer and a path-aware, typo-catching config reader | the schema |
| `io/Config` | the input-deck schema, validation, and model construction | numerics |
| `io/MeshReader` | Gmsh (MSH 2.2 / 4.1) and Abaqus / CalculiX `.inp` import: cell types, named sets, orientation repair, unused and duplicate nodes, units, the read report | boundary conditions or materials in the file |
| `io/CsvWriter`, `io/VtkWriter` | plain-text export with explicit precision | what to export |
| `io/StlWriter` | the boundary surface of a mesh (extruded for a plane one) as an indexed triangle surface, its closure and manifold checks, binary STL in and out | choosing what to export |
| `io/CalculixWriter` | one CalculiX input deck per load case for the same discrete problem, every field within CalculiX's 20 characters, and its non-linear counterpart (`*STEP, NLGEOM`, or small strain) with `*PLASTIC` and one step per leg of a load path | running CalculiX |
| `io/ResultWriter` | the result-directory layout, the geometry export and the summary documents | computing anything |
| `apps/` | argument parsing, orchestration, console reports, exit codes | physics |
| `python/sparlab_viz` | reading result files and drawing figures (plane fields in `fields`/`plots`, solid surfaces in `solid`/`plots3d`) | recomputing physics |
| `python/scripts/cross_validate.py` | driving CalculiX (static, `*BUCKLE`, `*HEAT TRANSFER`, `NLGEOM`, `*PLASTIC`) and scikit-fem (including an independent geometric stiffness, total Lagrangian solver and J2 solver) on the exported problems and comparing nodal displacements, temperatures and load factors | judging which code is right |
| `python/scripts/make_meshes.py` | the Gmsh parts of the real-geometry decks in linear and quadratic tetrahedra, their physical groups, deterministic mesher options | the analysis |
| `python/scripts/tet10_part_study.py` | meshing, solving and tabulating the Tet4 / Tet10 comparison on the engine mount | the element formulation |

## Design decisions worth naming

**Dirichlet conditions by partitioning, not penalty.** The reduced system
`K_ff u_f = f_f - K_fp u_p` keeps `K_ff` symmetric positive definite, gives
exact reactions from the full residual `r = K u - f`, and leaves the eigenvalue
spectrum undistorted. A penalty method would have made the modal analysis wrong
in a way that is hard to notice.

**One factorisation per optimiser iteration.** `StaticAnalysis::prepare()`
factorises `K_ff` once and `solve_load_vector()` reuses it, so a three-load-case
design costs one factorisation and three back-substitutions per iteration rather
than three factorisations. Combined with self-adjointness (no adjoint solve for
the compliance gradient), that is what makes multi-load-case optimisation on a
38 000-element mesh a 5-minute job rather than an hour. With the multigrid
solver the same interface means one hierarchy update per iteration and one
warm-started CG solve per load case: the objective keeps its solver, and the
solver keeps its aggregates and each load case's last solution, between
iterations.

**A direct solver as the reference, multigrid for scale.** The committed 2-D
decks and the verification studies run on the sparse Cholesky factorisation,
whose answer does not depend on a tolerance - except the multigrid study
itself, and the two finest Tet4 meshes of the simplex convergence study,
where `auto` selects multigrid (the summary records which ran). The
multigrid CG solver is the
path to large 3-D models; it is verified against the direct solver (to the
iterative tolerance, on every element type) rather than trusted, and `auto`
switches to it only above a size where the factorisation's fill starts to
dominate. Writing the multigrid rather than linking one kept the build at
Eigen and a C++17 compiler, and let it use what the model knows - node
coordinates for the rigid-body near-null space, the DOF blocks for the
strength of connection - which a generic algebraic package would have to be
told.

**Assembly into a cached pattern.** An optimisation re-assembles a matrix
whose pattern never changes, hundreds of times. The first assembly records
where every element block lands; later ones scatter into those slots in
element order - the order the triplet path sums duplicates in - so the
matrix is bitwise identical, which a test asserts, and the per-assembly sort
and triplet storage disappear. On the largest scaling mesh that was the
difference between assembly dominating a multigrid solve and being a small
part of it.

**The projection defines the physical density.** With the Heaviside
projection on, `rho_bar` rather than the filtered density is what the SIMP
law, the volume constraint, the stress constraint, the grey level and the
interpretation see, and the chain rule runs through it for every function.
Nothing downstream needed to know whether the projection is on, and the
finite-difference verification covers the combination.

**A mesh file is checked, not trusted.** The readers treat a file from
another tool the way a reviewer would: unsupported and mixed cell types are
refused with the re-export that fixes them, orientation is repaired and
counted, unused and coincident nodes are reported, the unit scale is checked
for plausibility, and everything - including the keywords that were ignored -
lands in the summary. Named groups become node and element sets, so a deck
addresses `"group": "bolt_holes"` instead of coordinates that would break
the moment the geometry changes.

**Element-matrix cache for uniform grids.** On a uniform structured mesh every
cell is geometrically identical, so `K_e^0` and `M_e^0` are integrated once. The
generic per-element path is always compiled in and the test suite asserts the
two agree, so the optimisation cannot silently diverge from the general case.

**Purely geometric region selectors.** Boundary conditions, loads and passive
regions are all specified by geometry, never by node or element indices (though
indices are available). The same deck therefore applies unchanged at any mesh
resolution, and a topology extracted from a density field can be re-constrained
without re-authoring anything - which is what makes the static and modal
analysis of the interpreted structure possible at all. That re-analysis runs on
every topology job, because the compliance of the *thresholded* structure is
what the design would really deliver; if the extracted structure turns out not
to be analysable, the reason is recorded in the summary and warned about rather
than allowed to destroy a completed optimisation.

**Passive regions as equal bounds.** Pinning the design variable rather than the
physical density keeps the objective an exact function of the free variables,
which is what makes the finite-difference verification of the gradient
meaningful.

**A hand-written JSON parser.** The build needs nothing but Eigen and a C++17
compiler. The parser accepts `//` and `/* */` comments and a trailing comma,
because an annotated input deck is far easier to maintain, and it reports the
line and column of a syntax error. The `ConfigNode` wrapper records every key it
reads so unknown keys can be reported: a misspelled key that silently takes its
default is one of the easiest ways to publish a wrong number.

**The Python layer never computes physics.** It reads `mesh.json`,
`summary.json` and the CSV tables. Derived values in figures are limited to
ratios and least-squares slopes, and each is labelled as derived. That means a
number in a figure can always be traced to a solver run.

**One dimension-generic core, not two solvers.** Rather than a 3-D copy of
the 2-D code, the dimension is a run-time property of the mesh. What that
bought: every verification study, the optimiser, the writers and the
diagnostics run on the solid path unchanged, and the plane path is protected
by the fixed-size Q4 kernels and the byte-for-byte reproduction check. What it
cost: a handful of `dim`-dependent branches in the constitutive law, the
von Mises formula, the rigid-body-mode table and the face tables.

**MMA returns the iterate it judged.** Convergence is decided on the iterate
whose objective and constraints were just evaluated, and that iterate is the
design returned - not the unevaluated update, which under the objective-stall
criterion can still differ by up to the move limit and need not be feasible.
A stalled objective on an infeasible design is reported as non-convergence,
never as a result.

**Geometry export says what it is.** The STL surface is the boundary of the
thresholded cells, so it inherits the mesh's staircase and every
interpretation choice. The writer checks closure (every directed edge has its
reverse) and manifoldness (no edge used twice each way), compares the enclosed
volume with the cell volume, and records all of it in the summary; a surface
that is closed but non-manifold because cells touch along an edge is reported
as such rather than passed off as a solid.

## Extension points

The interfaces were shaped so the following are additions, not rewrites.

**A new element topology.** Implement `Element` (the kernels plus the local
face table), add an `ElementType` enumerator, and register it in
`make_element()`. That is exactly how `Hex8`, and later `Tri3` and `Tet4`,
were added: the layers above them are written against `dim()`,
`num_nodes()`, `num_faces()` and `num_voigt()` and needed no change beyond
the writers' cell-type tables (VTK, CalculiX, STL faces). `Mesh` stores
connectivity as a flat array with an explicit stride, so one mesh holds one
cell type; mixed meshes (a quad-dominant plane mesh, a hex mesh with prisms)
need that stride replaced by a per-element offset table, which is the one
change outside the new subclasses. The readers refuse such files today and
say how to re-export them.

**Plane strain.** Already implemented and unit-tested end to end: set
`model.stress_state` to `plane_strain`. The constitutive matrix, the von Mises
formula's out-of-plane term and the test coverage are all in place. Plane stress
remains the default and the idealisation the verification studies exercise, so
plane strain is documented as available rather than as validated to the same
depth.

**A new material model.** `FemModel` holds one `IsotropicMaterial` and hands its
constitutive matrix to the element kernels. An orthotropic or
temperature-dependent material replaces that with a small interface returning
`D` per element; nothing above the element layer depends on the material being
isotropic.

**Another mesh format.** `io/MeshReader` reads Gmsh and Abaqus / CalculiX
files into a `Mesh` with named sets and a `MeshReadReport`; a third format
is one more function with the same signature, dispatched on `format` or the
file extension. The checks after parsing (orientation, unused and duplicate
nodes, units, quality) are shared, so a new reader only has to produce
coordinates, cells and sets.

**A different optimiser.** `ComplianceObjective::evaluate` returns the objective,
the volume constraint and both gradients; `optimality_criteria_update` and
`MmaOptimizer::update` are the two update rules, selected by
`optimizer.method`. An augmented-Lagrangian or a trust-region step would be
a third, again without touching the objective.

**A new constraint.** `StressConstraint` is the template, and
`BucklingConstraint` the second instance of it: it takes the model,
the assembler and the filter, evaluates its value and gradient from an
`ObjectiveEvaluation` (reusing the objective's solver for its adjoint solve
through `ComplianceObjective::solve_adjoint`, and its chain rule through the
filter and projection through `ComplianceObjective::chain_to_design`) and
hands one row of `fval` / `dfdx` to MMA. A frequency constraint, a
displacement bound or a second volume budget on a region fits the same
shape.

**A new objective.** The self-adjoint shortcut in `Sensitivity.cpp` is specific
to compliance. A non-self-adjoint objective (a displacement at a point, a
frequency) needs an adjoint solve, which the stress constraint already
exercises: the solver is kept and reused for the adjoint right-hand side,
with its own warm start per adjoint.

**Another linear solver.** `LinearSolver` is the interface: `factorize`,
`solve`, `solve_from` (with an initial guess), `set_layout` (the node
coordinates and DOF map an algebraic multigrid needs), and the statistics
the summary records. A new backend - a supernodal Cholesky, a GPU CG - is a
subclass and an entry in `make_linear_solver()`; the verification studies'
solver-agreement check then covers it.

**A new external code for cross-validation.** `CalculixWriter` shows the
shape: write the same nodes, connectivity, supports and nodal loads in the
target's format, then add a `solve_with_<code>` function and a `compare`
entry in `python/scripts/cross_validate.py`, which already handles the
node-by-node comparison, the tolerances and the summary.

## File map

```
include/sparlab/          public headers, one per component, documented
src/                      implementations, mirroring include/
apps/                     four command-line drivers plus shared CLI support
tests/                    eighteen Catch2 translation units plus shared fixtures
configs/benchmarks/       the fifteen benchmark decks (four plane compliance
                          cases and the MBB beam's projected, robust and
                          overhang variants, the stress-constrained L-bracket,
                          the buckling-constrained column, the solid bracket,
                          its projected, overhang and 356k-DOF variants, and
                          the two parts read from mesh files)
configs/meshes/           the Gmsh / Abaqus meshes those parts read (the
                          engine mount also in curved Tet10), each with a
                          .json record of how it was generated
configs/studies/          the design-study baseline deck
configs/verification/     the static decks the cross-validation solves (Q4,
                          Hex8, Tri3, Tet4, the lug at nu = 0, the Tet10
                          engine mount, and the buckling columns)
python/sparlab_viz/       loaders, style, field artists (plane and solid), figures
python/scripts/           figure and table drivers, the cross-validation driver
scripts/                  run scripts (benchmarks, verification, cross-validation,
                          the Tet10 study, the design study, scaling)
docs/                     this documentation
docs/figures/             committed figures and the animations
docs/results/             committed result tables, regenerated from summaries
results/                  generated output (git-ignored)
```
