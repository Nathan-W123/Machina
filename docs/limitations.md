# Assumptions and limitations

This is what SparLab does **not** do, and what its results do **not** claim.
Nothing here is a bug list; every item is a deliberate scope boundary, and each
says what would be needed to lift it.

## Physics and idealisation

**Plane or solid continuum, nothing in between.** Plane stress (default),
plane strain, and three-dimensional solids on hexahedral or tetrahedral
meshes. There are no
plate, shell or beam elements, so a thin-walled part is either a plane
load-path study or a solid mesh with enough elements through the wall to
resolve its bending - which the solid element needs several of (below).
Buckling is seen only as far as the continuum model allows: the linear
buckling analysis finds the in-plane buckling of a plane model's struts and
every mode of a solid mesh fine enough to bend, but a plane model cannot
buckle out of its plane, so plate buckling of a thin lightened web - which
can govern it - stays invisible to the 2-D decks.

**Elastic and linear unless asked otherwise.** The static, modal and
buckling analyses and the topology optimisation are linear and elastic:
small strain, small displacement, one factorisation - a material's
`plasticity` block is ignored there (the solve and optimise apps say so). The
deformation figures of those runs are exaggerated by a stated factor purely
for visibility. `sparlab_solve` adds a non-linear static analysis
(`nonlinear` in the deck; `docs/formulation.md`, sections 7c and 7d): large
displacement and rotation, with the Saint Venant-Kirchhoff law (small strain)
or a compressible neo-Hookean law (large strain), follower pressures, a
centrifugal load at the deformed position, the finite-strain thermal split,
J2 plasticity with isotropic (linear and Voce) and kinematic (Prager)
hardening in small strain or in the Green-Lagrange strain at large rotation,
load paths that unload and reverse, load control that stops at limit and
bifurcation points and at a plastic collapse and says so, and the arc-length
method through limit points. What it does not do:

* plasticity is rate-independent J2 at small strain: no finite-strain
  plasticity (with `finite` kinematics the J2 return in the Green-Lagrange
  strain is sound for large rotation with small strains, and the run warns
  beyond a strain of 0.05), no creep, viscoplasticity, damage, fracture,
  non-associative or pressure-dependent yield, nonlinear (Armstrong-Frederick)
  kinematic hardening, anisotropy or temperature-dependent properties; the
  elastic part of the neo-Hookean law cannot be combined with plasticity;
* constant-strain elements (Tri3, Tet4) can lock under the isochoric flow of
  a fully plastic state in plane strain and 3-D, depending on the mesh
  pattern, and fully integrated Q4 and Hex8 lock without the mean
  dilatation that is their default - a locked collapse load comes out high,
  which is unconservative;
* contact (`contact` block) is small-sliding - the contact geometry of the
  reference configuration, a gap linear in the displacement - with
  `small_strain` kinematics only; tools that travel far over a surface need
  the forming analysis (`sparlab_form`, `docs/forming.md`), whose penalty
  contact is evaluated in the current configuration but is rigid-tool only,
  approximate (a penetration of about `p h / (s E)`), implicit and
  quasi-static, and has no remeshing, trimming, thermal effects or tool spin
  (its own limitations are listed there); the static path is quasi-static,
  and a snap-through that the arc-length method follows is a sequence of
  equilibria, not the dynamic jump a real structure would make - the
  non-linear transient (below) integrates that jump, with its own limits;
* at a bifurcation of a perfect structure there is no branch switching: load
  control stops there, and the arc-length method stays on the fundamental
  path. A post-buckling analysis needs an imperfection built into the mesh
  (a perturbed geometry), which the deck does not generate;
* the Saint Venant-Kirchhoff law is valid for small strain only - the summary
  warns above a Green-Lagrange strain of 0.05, and under strong compression
  (a stretch below `1/sqrt(3)`) it softens unphysically. The neo-Hookean law
  has no plane-stress form and no thermal strain; its constants are the small-
  strain `E` and `nu`, not a fit to test data of a real elastomer;
* a follower pressure over free edges makes the tangent non-symmetric; LU
  factorises it, which reports no inertia, so the stability of such a state
  is not assessed (the summary says so);
* body forces and self-weight stay dead loads per reference volume, a
  prescribed displacement scales with the load factor, and every load of a
  case shares one load factor - a `load_path` scales them all together
  (load, unload, reverse), but there are no load sequences of different
  loads (a preload followed by a service load) within a case;
* the topology optimisation, the modal analysis and the buckling analysis do
  not use the non-linear state (`sparlab_topopt` refuses the block).

**Buckling is linear bifurcation.** The buckling check is the eigenvalue
problem `(K + lambda K_G(u)) phi = 0` of the linear static state: the
bifurcation load of the perfect geometry, which is an upper bound on the
collapse load of a real part with its imperfections, residual stresses and
plasticity. No imperfection sensitivity, no follower loads, no knock-down
factor. A load factor of 6 is not a safety factor of 6 against collapse; for
a shell-like or imperfection-sensitive structure the difference can be
large. The non-linear analysis gives the limit load of the geometry as
meshed - an imperfect mesh included - and brackets the first loss of
stability under load control, but not a post-buckling branch that starts at
a bifurcation of the perfect geometry.

**Dynamics: direct integration with a constant step, and the direct harmonic
response.** The transient analysis integrates `M a + C v + K u = A(t) f` (or
its non-linear counterpart) by HHT-alpha with a constant step, and the
harmonic response solves the complex dynamic stiffness at every frequency.
What they do not do:

* one amplitude scales every load and prescribed displacement of a case: no
  loads that follow different time histories within a case, no initial
  velocity field, and no load applied at a moving position;
* the step is constant and chosen by the user: no automatic step control by
  an error estimate (a failing Newton step is halved, and the halves grow
  back), and no explicit (central-difference) integration for short
  impact-like events, which the implicit method can only follow with a small
  step at a high cost per step;
* damping is Rayleigh's (and, for the harmonic response, a structural loss
  factor): no modal damping ratios, no discrete dampers, no frequency-dependent
  or viscoelastic material, and in a non-linear run the Rayleigh stiffness
  term uses the linear elastic stiffness, not the tangent;
* the non-linear transient starts at rest (a preloaded start is a static
  analysis of its own) and has the non-linear statics' limits (section
  above); the trapezoidal rule conserves the energy of a non-linear system
  only to `O(dt^2)`, and the run warns when the balance shows energy created
  beyond 1 % - it does not bound the error of the motion itself;
* no modal superposition (neither for the transient nor for the harmonic
  response), no response spectrum, random vibration or power spectral
  densities, no rotor dynamics (Coriolis and gyroscopic terms), no
  fluid-structure coupling, no wave-absorbing boundaries, and no fatigue;
* the harmonic response is linear: the load case's linear stiffness, no
  prestress from a static preload (the geometric stiffness of a preloaded
  structure is not added), and no harmonic balance for a non-linear
  structure;
* a natural frequency of the modal analysis is the undamped eigenvalue of the
  constrained model; the dynamic analyses are not coupled to the topology
  optimisation (`sparlab_topopt` refuses their blocks).

**Body loads: self-weight, force densities and steady rotation.** Gravity,
uniform body force densities on element regions and the centrifugal load of a
steady rotation are integrated exactly from the consistent mass. Rotation is
the static centrifugal load only: no Coriolis or gyroscopic terms. In the
linear analyses it acts at the undeformed position; the non-linear analysis
applies it at the deformed position, with the spin-softening stiffness that
follows. There is no inertia relief (a
free-flying body balanced by its own acceleration), so a body load must be
reacted by supports. The wing-rib deck's "fuel inertia" case is still a
*representative edge pressure*, not a body-force calculation.

**Thermal strain is linear thermoelasticity.** A temperature field - uniform,
given on regions, or the solution of steady conduction with fixed
temperatures, surface fluxes, convection and heat generation - enters as the
free strain `alpha (T - T_ref)` with temperature-independent `E`, `nu`,
`alpha` and `k`. The conduction is steady (no transient heat transfer, heat
capacity or time), has no radiation, no contact conductance between parts and
no heat from deformation; the temperature is solved first and does not depend
on the displacement. Materials may differ by element region
(`material_regions`), each with its own reference temperature, but the
CalculiX export needs one common reference temperature (CalculiX measures
thermal strain from the initial nodal temperature). First-order hexahedra and
quadrilaterals integrate the interpolated temperature (the consistent load),
which cannot represent the free expansion of a linear temperature gradient
exactly - that needs a quadratic displacement - so a free Hex8 or Q4 part in
a gradient shows a small spurious stress that vanishes under refinement; the
Tet10 represents it exactly. In the non-linear analysis the thermal stretch
`1 + alpha dT` splits off multiplicatively (Saint Venant-Kirchhoff only),
with the same temperature-independent constants. Topology optimisation does
not take body loads, temperatures or several materials yet (`sparlab_topopt`
refuses them).

**Linear elements, and one quadratic element.** The four-node
quadrilateral, the three-node triangle, the eight-node hexahedron and the
four-node tetrahedron are stiff in bending and need several elements through
a bending depth; the constant-strain simplices are much stiffer than the
bilinear and trilinear elements. The mesh-convergence studies quantify it:
the coarsest solid cantilever is 16 % too stiff on Hex8 and 50 % on Tet4. The
ten-node tetrahedron removes most of this for solids: on a 10 x 2 x 1-cell
cantilever, where the Hex8 is 33 % and the Tet4 66 % too stiff, it is within
0.06 % of beam theory, and on the engine mount meshed from CAD the Tet4 compliance
is 19 % below the finest Tet10 run on the benchmark's own mesh
(`docs/benchmarks.md`, section 11). There is no quadratic plane element
(Q8, Tri6) and no Hex20, so a 2-D model or a hexahedral mesh still needs
refinement to bend. On a Tet10 cell with curved edges the 4-point stiffness
rule is not exact - the integrand is rational - as in every code that uses
C3D10; the study's curved and straight-sided meshes of the same cells agree
to 0.13 % at 3 mm. Enhanced assumed strain or B-bar would help the linear
elements per DOF; `docs/architecture.md` says what adding an element
involves.

**Meshes: one cell type per mesh.** The structured generators make boxes of
quadrilaterals, triangles, hexahedra or tetrahedra, linear or (tetrahedra,
`mesh.order: 2`) quadratic. Real geometry comes in through Gmsh (MSH 2.2 /
4.1, ASCII) and Abaqus / CalculiX `.inp` files, with these limits: one cell
type per mesh (a quad-dominant or a hex-dominant mesh with prisms is
refused, with the re-export that fixes it); the only second-order cell is
the ten-node tetrahedron (C3D10, Gmsh type 11), whose edge nodes may follow
a curved surface - every other boundary is a polygon of cell edges, and
elevating a Tet4 mesh keeps its facets;
Abaqus parts must be flat (no instances with translations, one part);
supports, loads and materials in the file are ignored in favour of the deck;
and SparLab does not mesh - the two committed parts come from
`python/scripts/make_meshes.py`, which needs Gmsh. There is no adaptive
refinement, and on a graded mesh a filter radius given in cells
(`radius_elements`) multiplies the *mean* cell size, so it is a different
physical radius in the fine and the coarse regions than a reader might
assume; a radius in metres avoids that.

**Point loads are singular.** A concentrated force on a node produces a stress
that does not converge under refinement - the peak grows without bound. Reported
peak stresses near a point load are mesh-dependent artefacts. The benchmark
decks spread their resultants over a small region for this reason, and the
distributed-traction path is the honest way to apply a pressure.

**Plane strain is verified on curved sections, not on beams.** The thick
cylinder under pressure, the conducted thermal cylinder, the rotating
cylinder and the finite-strain tube run in plane strain against exact
solutions, the Hex8 and Tet10 sections held at `u_z = 0` reproduce the Q4
answers, and CalculiX's plane-strain expansion (`CPE4`) agrees node by node.
The beam and frequency studies against Timoshenko and Euler-Bernoulli run
in plane stress only.

## Attachment and joint representation

The bolt holes, bearing collars and pin lug in the aerospace bracket are
**representations**, not joint models:

* a bolt shank is idealised as perfectly rigid - every node inside a bolt hole
  is fully fixed. Bearing compliance, fastener flexibility, hole clearance and
  bolt preload are absent;
* the lug load is applied to the ring of nodes around the hole as if bearing
  pressure were uniform. A real pin joint has a cosine-like bearing distribution
  and local yielding;
* passive solid collars stand in for machined bosses. Their size was chosen to
  look like plausible hardware, not derived from a bearing-stress allowable;
* no fretting, no bushing, no through-thickness load transfer.

A real joint analysis needs 3-D contact and a fastener model. Nothing here
should be read as a joint substantiation.

## Loads

**The load cases are illustrative.** The bracket's 9 kN down-load, its 4.5 kN
reversal and its 4.5 kN lateral case, and the wing rib's skin pressures, are
plausible magnitudes chosen to exercise multi-load-case optimisation. They are
**not** derived from a flight-loads analysis, not factored to any airworthiness
requirement, and carry no limit/ultimate distinction.

**The objective is weighted compliance, not a load envelope.** Minimising
`sum_l w_l C_l` makes the design stiff on *average* across the weighted cases.
It does not guarantee any single case meets a deflection limit, and it is not
the same as optimising for the worst case. A worst-case (min-max) formulation
needs a different objective and a bound-constrained optimiser.

## Topology optimisation

**A density field is not geometry.** The SIMP result is a relative material
distribution over a fixed mesh. Elements at intermediate density represent
neither solid nor void; they are an artefact of relaxing a 0/1 problem into a
continuous one. SparLab never presents a density field as a solid body without
saying how it was converted: the "interpretation" block in every topology
summary records the threshold, how many elements survived, how many
disconnected groups they formed and how much material was discarded as islands.
The modal comparison of an "optimised structure" is run on that explicitly
extracted sub-mesh.

**Manufacturability: two rules modelled, nothing more.** The overhang
filter models one additive-manufacturing rule - no material without support
directly or diagonally below it, a 45-degree limit on square or cubic cells -
and the robust formulation models one uniform manufacturing error, the whole
part coming out thinner or thicker by the erosion. A design that passes the
overhang check is printable *under that rule*: support removal, residual
stress and distortion, surface finish, anisotropic material, a minimum wall
and the build plate's own limits are not modelled. Nothing models a draw
direction, milling-tool access, fillet radii or a machining setup. The
filter radius and the robust formulation impose a minimum *length scale*,
which is a geometric property, not a process. Without these options the
designs are load-path guidance for a designer, not parts; with them they
are still not certified parts.

**The optimum is local.** SIMP with `p > 1` is non-convex; a different starting
design, penalty schedule or mesh can converge to a different local optimum. The
design study shows exactly that: the mesh-refinement arm with the filter radius
fixed in cells produces visibly different topologies at different resolutions.
Continuation reduces the dependence but does not remove it.

**The projection sharpens the design; alone it does not give it a length
scale.** A density filter alone leaves an intermediate band about one radius
wide at every material boundary (the MBB benchmark settles at a grey level
of 0.28, the solid bracket at 0.34). The Heaviside projection removes most
of it - the projected solid bracket ends at 0.013 - and with it the gap
between the objective and the thresholded structure. A single projection at
`eta = 0.5` does not secure a minimum member size: features thinner than the
filter radius survive, and on the projected MBB beam two diagonals vanish
under a 0.1 erosion and cost 22 % of the stiffness. The robust formulation
(`topology.projection.robust`) is the tool for one, and the length-scale
check measures what a design delivers. The size follows from the filter
radius and `robust_delta`; it is not a stated number the solver enforces.
At `robust_delta = 0.1` the MBB's smallest member grew only from 3 to 4
cells, and at 0.05 the robust column still has members about one cell
thick. At a large `beta` the optimisation also becomes
strongly non-convex: single elements at the solid-void interface flip by the
full move limit while the compliance stands still, which is why the
projected decks stop on the objective-stall criterion rather than on the
design change. The `beta` schedule is a further setting the result depends
on, and it is recorded with every run.

**Three constraint types.** Optimality criteria handles the volume
constraint alone (of the dilated design in a robust run); MMA handles the
volume plus one aggregated stress constraint and one aggregated buckling
constraint per load case. Per-mode frequency constraints, displacement
bounds and multiple volume budgets on separate regions are not written -
`docs/architecture.md` says what each would take, and `StressConstraint` and
`BucklingConstraint` are the templates - and the OC path remains the default
because it converges with no tuning where MMA needs a move limit chosen per
problem.

**The buckling constraint acts on the SIMP model, not on the part.** The
constrained load factors are those of the density field, in which
intermediate densities keep `rho^p` of their stress stiffness. The exported
part - thresholded, full material, grey members gone - can buckle earlier:
on the column benchmark the constraint holds `lambda >= 6` on the SIMP model
while the part buckles at 3.65 with a plain projection and 5.81 with the
robust formulation. The buckling check of the exported part is therefore
the number to judge, and a design should be re-run with a margin when it
falls short. The KS aggregate overestimates the largest ratio by at most
`ln(m)/P` and is conservative; the sensitivity formula holds for simple
eigenvalues, and at a repeated eigenvalue only the aggregate stays
differentiable. Void material can still produce pseudo modes; the
energy-fraction diagnostic names them but does not remove them. The cost is
one eigensolve and one adjoint per aggregated mode every iteration - 242
linear solves per iteration on the column against one without the
constraint - and the benchmark needed a move limit of 0.05 and 3000 Newton
iterations for MMA's subproblem.

**The overhang filter needs a structured grid, and OC struggles with it.**
The supports of an element are grid cells of the layer below, so the filter
and its check run only on `structured_quad` and `structured_hex` meshes,
with the build direction along a grid axis; the overhang angle is
`atan(layer thickness / cell width)`, 45 degrees on square cells. A
45-degree staircase touches its supports along an edge, which the
interpretation's face connectivity does not count as a joint, so a filtered
3-D design can leave small islands (0.7 % of the bracket's material). With
OC, the filtered MBB and bracket runs locked into period-2 cycles at the
final `beta`; both overhang decks use MMA.

**The robust formulation needs care with OC.** Each iteration projects the
filtered field twice more, for the dilated and the blueprint designs,
without an extra solve: only the eroded design is analysed. The dilated
design's volume target is rescaled from the blueprint's, so the blueprint
meets the volume fraction only as closely as the ratio of the two volumes
holds between rescalings. On the robust MBB, which rescales every 20
iterations, interface elements flip at `beta = 32`, the blueprint's fraction
alternates between about 0.498 and 0.502, and the returned design is 0.41 %
below it. Rescaled at every iteration, OC locked into a period-2 cycle at
`beta = 32` that met no stopping criterion.

**The stress constraint bounds a relaxed aggregate, not a part's stress.**
The constrained quantity is `rho^q sigma_vm` at each element centre,
aggregated with a p-norm whose scale is re-fitted every iteration. Three
consequences: intermediate-density elements carry a stress the material would
not; the aggregate underestimates the true maximum between re-fits, so the
relaxed peak hovers about the limit rather than sitting under it; and a
stress concentration the mesh does not resolve is not bounded at all - a
re-entrant corner in a density design is mesh sensitive whether or not the
constraint is on. Every constrained run therefore re-solves its thresholded
structure with full material and reports that peak against the limit, which
is a consistency check, not a yield or fatigue substantiation. The reported
macroscopic stress in a low-density element is still an average over a
near-void cell and not a material stress, and the figures still mask elements
below the interpretation threshold.

**MMA needs tuning where OC does not.** With a moving constraint surface the
default move limit of 0.2 left the stress-constrained L-bracket oscillating
at its iteration cap. At 0.1 it still spends long stretches in a cycle of
period 5 at the move limit and needs about 500 iterations to converge on the
design change; one of its subproblems needed more Newton iterations than the
original budget of 200 per barrier level. Asymptote parameters, the p-norm
exponent and the subproblem budget are further knobs, each recorded in the
summary. A run that hits a cap says so, but it does not say which knob to
turn.

**The design change stalls on fine meshes.** OC's `max |dx|` measure does not
reach `1e-2` on the finer benchmark meshes because thin members migrate one cell
at a time long after the objective has settled. A measured example: on the
cantilever benchmark the compliance improves by only 0.35% between iteration 250
and iteration 800, while `max |dx|` stays between 0.004 and 0.09 in nine
iterations out of ten and reaches the move limit of 0.2 around iteration 260.
This is why an objective-stall criterion is the backstop - the relative spread
of the compliance over a window, which an oscillating design cannot satisfy
mid-cycle - and why every run records both indicators and which one fired.

**Compliance is not monotone by construction.** OC is a fixed-point update with
a move limit, not a line-search method, and each continuation stage raises the
penalty, which raises the compliance of a fixed density field. Histories are
recorded and plotted unsmoothed.

**The exported geometry is the interpretation, staircase and all.**
`structure_after.stl` is the boundary of the cells at or above the threshold
- a voxel surface, not a smoothed or fitted one - and inherits every choice
the interpretation made. Two retained cells that touch only along an edge
leave a non-manifold edge in it (32 on the 3-D bracket); the writer reports
the count and a slicer may split the surface into shells there. A plane
design is extruded by its thickness. None of this is a printable or
machinable part: no smoothing, no minimum feature size beyond the filter
radius, no fillets, no draft.

## Numerics

**The automatic solver choice is a size rule.** `auto`, the default,
factorises up to 50 000 free unknowns in 2-D and 10 000 in 3-D and uses
multigrid CG above that. The limits sit above the measured single-solve
crossover on purpose (`docs/benchmarks.md`, *Runtime and scaling*): a
factorisation serves every load case and every modal solve, and it has no
iteration count that can go wrong. They are not tuned per problem. A deck
that names a solver gets that solver, and `solver.linear.auto_direct_limit`
moves the limits.

**The direct solver bounds the solid problem size.** `SimplicialLDLT` with
AMD ordering is what makes a small optimisation loop cheap, with one
factorisation shared by every load case and no adjoint solve for
compliance. In 3-D its measured cost grows like `n^2.35`: 43 s per
factorisation at 47 775 DOFs, with 789 factor entries per DOF. No run here
factorises a solid system above that size. A nested-dissection ordering
(METIS) or a supernodal factorisation (CHOLMOD) would move the limit, and
neither is a dependency.

**Multigrid is tested on what the benchmarks need, not in general.** The
smoothed-aggregation hierarchy uses the rigid-body modes as its near-null
space, and its strength threshold keeps solid and SIMP void in separate
aggregates. On the problems tested the CG iteration count barely moves
with the mesh: 14 to 17 on Q4 and Hex8 blocks up to 830 115 unknowns, 16
to 21 on Tet4, and 18 to 24 per solve inside the optimisation runs with
their `1e-9` stiffness floor. Not tested: nearly incompressible material
(`nu` near 0.5), where the rigid-body modes are not a sufficient near-null
space; badly shaped or strongly graded meshes; and plane strain. A solve
that does not converge is reported with the residual reached and the keys
to change, never returned as an answer. The setup is not free either. For
one solve from scratch, Jacobi-preconditioned CG is as fast as multigrid up
to 28 413 unknowns on the Hex8 block, and a modal analysis with CG repeats
the whole solve for every vector of the subspace. The two modal analyses
of the 356 475-DOF bracket take a fifth of its run.

**Element loops are serial; only the iterative solvers are threaded.**
Assembly, sensitivities, filtering and stress recovery loop over elements
on one core, and the sparse Cholesky factorisation is serial too. The
multigrid kernels use OpenMP, and Eigen threads the matrix-vector product
of its Jacobi-preconditioned CG. Both give the same bits on one thread as
on many; a test asserts it for multigrid. Direct-solver runtimes are
therefore single-core numbers and iterative-solver runtimes are
four-thread numbers. Parallel element loops would need per-thread buffers
or a colouring to keep the results deterministic.

**The SIMP stiffness floor limits conditioning.** With `emin_ratio = 1e-9`
the smallest LDL^T pivot ratio sits around `1e-9`, comfortably above the
`1e-14` singularity threshold, but the condition number of `K_ff` is
correspondingly large. The multigrid solver copes because void and solid do
not share aggregates. Jacobi-preconditioned CG is kept as a baseline and is
measured on uniform material only; no benchmark uses it.

**Mechanism modes are not detected analytically.** The pre-solve
diagnostics detect rigid-body under-constraint exactly, per connected
element group, and detect floating regions. An internal mechanism, such as
two blocks joined at a single node, is not caught analytically. With the
direct solver it surfaces as a non-positive Cholesky pivot, which the
solver reports as a singular system with the same hint text. With the
multigrid solver it should surface as a singular coarsest level or a solve
that does not converge. A rigid-body mode reaching the coarsest level is
tested, an internal mechanism is not, so a model suspected of one is better
checked with the direct solver.

**Nodal stress averaging is simple, not superconvergent.** Nodal values are
area-weighted averages of adjacent element values, used for smooth contours
only. A superconvergent patch recovery would give better nodal accuracy; element
values, which are what the CSV tables and the reported peaks use, are unaffected.

## Verification and validation

**Verification is strong; validation is narrow.** The implementation is
verified against exact answers on all five element types. The patch tests
pass to 4e-15 (Q4), 2e-15 (Hex8) and 9e-15 (Tri3 and Tet4), and the Tet10
passes the quadratic (pure-bending) patch test to 1.2e-14. The sparse and
dense solvers agree to 9e-12, and multigrid CG agrees with Cholesky to
4e-12. The compliance gradient matches central differences to 2e-8 on Q4
and Hex8, and to 1e-7 through the Heaviside projection on Q4 and Tet4; the
buckling constraint's gradient to 1.7e-6 and the overhang filter's to
2.8e-6. Mass is conserved to 5e-14. A homogeneous large deformation on
distorted meshes is reproduced to 1e-9 m, and the tube at finite strain
converges at every element's order to its exact solution. A transient run
equals the exact solution of its discrete equations (every mode integrated
exactly, in extended precision) to 3.3e-10, a rod's harmonic response its
exact discrete solution to 2.9e-10, and the rod's harmonic and transient
responses and a one-element oscillator's finite-strain elastic and
elastoplastic motion converge at second order to their exact solutions.
Linear static displacements are cross-validated node by node against two
independent codes on thirty-eight problems with sixty-three load cases,
covering all five element types and both mesh-file formats; so are buckling
load factors on three columns, large-deflection states on four decks,
elastoplastic states on ten, transient histories on five and harmonic
responses on two. scikit-fem agrees to solver round-off, displacements and
load factors alike: every linear difference lies below the round-off scale
of its system (its condition number times eps), 1.5e-10 or better except on
the three worst-conditioned systems - a plane-strain strip and two slender
Tet10 cantilevers, 1.4e-9 to 2.2e-8; scikit-fem's own total Lagrangian
solver agrees with SparLab's non-linear states to 1e-14, and its own J2
solver with the elastoplastic states, small strain and finite, E-bar
included, to 7e-12; an HHT-alpha integration and a direct complex solve on
its matrices agree with the transient histories to 2.2e-8 (the final
acceleration of an elastoplastic run; its displacements to 1.2e-12) and
with the harmonic responses to 1.0e-10. CalculiX's displacements agree to
the rounding of its own result file, 4.3e-6 or better, linear, `NLGEOM`,
`*PLASTIC` and `*DYNAMIC` alike, wherever the two codes solve the same
discrete problem; its plane-stress
comparisons at `nu != 0` and its finite-strain plasticity are recorded but
not judged, its kinematic hardening is not compared (it does not
reproduce Prager's rule), and neither are its plane elements in dynamics
(their `*DYNAMIC` response contradicts CalculiX's own `*FREQUENCY`); its
buckling factors differ by up to 8.3e-5 for a reason not identified.
Validation against independent theory covers exactly five references:
Euler-Bernoulli and Timoshenko cantilever deflection, Euler-Bernoulli
bending frequencies, fixed-free rod axial frequencies, the Euler-Engesser
buckling load of a clamped column, and Euler's elastica. Stresses,
frequencies, non-linear static load paths and optimised designs are not
compared with another code, and there is **no comparison against
experiment**.

**The MBB compliance is not compared to a published number.** SparLab reports
what it computes (218.8 J with the density filter, 203.2 J with the sensitivity
filter, at the 60x20-equivalent radius on a 180x60 mesh). Which value a
published result corresponds to depends on the filter variant, the radius
convention, the convergence criterion and the penalty schedule, so quoting
agreement without re-deriving the reference under the same settings would be
guesswork. The internal consistency checks - the verified gradient, the exact
volume constraint, the KKT spread across interior elements - are what support
the result.

## Reporting

**The mass-stiffness comparison is 2-D specific.** In the plane idealisation
`K` and `M` both scale linearly with thickness, so a uniformly thinned plate
has exactly the full-solid natural frequencies while its compliance rises as
`1/thickness`. That is what makes "equal-mass uniform plate" a well-defined
baseline, and the runs verify the frequency invariance numerically rather
than assuming it. A solid has no thickness to thin, so the solid cases
report only the full solid domain as their reference and no "stiffness
gain"; their `mass_stiffness_comparison.note` says so. The baseline is also
not fair where passive voids remove part of the domain from the design but
not from the plate, which is why the L-bracket quotes its unconstrained run
instead. The lug bracket's holes and solid rims are in both the part and the
plate, so there the comparison is like for like.

**Without a projection, solid density designs are greyer and the re-solve
matters more.** A 1.5-cell filter averages 17 hexahedra rather than 9
quadrilaterals, and the intermediate boundary layer is a surface rather
than a curve. The 3-D bracket without a projection therefore settles at a
grey level of 0.34, and its thresholded structure is 39 % stiffer than the
SIMP field. The Heaviside projection closes most of that gap: the projected
bracket ends at a grey level of 0.013 with its thresholded structure within
2 % of the objective, and the 356 475-DOF bracket within 0.03 %. The engine
mount, whose projection stops at `beta = 16` on a tetrahedral mesh, still
differs by 4.5 %. The re-solved value is the one to plan with, because it
is the analysis of the part that is actually exported.

**A projected run reports the projected density.** With the projection on,
the compliance, the volume, the grey level and every figure refer to the
projected density, the one the finite-element model sees. The grey level
before projection is recorded beside it (`filtered_grey_level`) and is much
higher: 0.12 to 0.35 on the benchmarks. The volume constraint holds at every
design update, but two kinds of iterate are analysed before an update has
put them on the target: the uniform start seen through the projection, and
the first design after each increase of `beta`. Their recorded volume is off
by up to 16 % at the start (the lug bracket) and up to 1.8 % after a `beta`
step (the projected solid bracket). The convergence figure says so beside
those points.

**Runtimes are machine-specific.** Every reported time comes from one machine
and one build, both recorded in the result summaries. Absolute values will
differ elsewhere; the scaling exponents are the portable part, and even those
carry cache effects at these sizes.
