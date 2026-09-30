# Conventions

Everything in SparLab is strict SI. There are no unit-conversion helpers and no
implicit scaling anywhere in the code, so a value read from a result file is
already in the unit its column header states. The one conversion factor is
explicit: a mesh file drawn in millimetres is read with `mesh.scale = 0.001`,
applied once as the file is read and recorded in the summary. A mesh whose
extent comes out above 20 m or below 0.1 mm draws a warning naming that key.

## Units

| Quantity | Unit | Symbol in code / files |
|----------|------|------------------------|
| Length, coordinate, displacement | metre | `m` |
| Out-of-plane thickness | metre | `thickness_m` |
| Force, reaction | newton | `N` |
| Traction, stress, Young's modulus | pascal | `Pa` |
| Mass density | kilogram per cubic metre | `kg/m^3` |
| Mass | kilogram | `kg` |
| Area | square metre | `m^2` |
| Volume (area x thickness) | cubic metre | `m^3` |
| Strain, density, volume fraction | dimensionless | `[-]` |
| Compliance, strain energy | joule | `J` (= N m) |
| Eigenvalue `lambda` | 1 / second squared | `1/s^2` |
| Angular frequency `omega` | radian per second | `rad/s` |
| Frequency `f` | hertz | `Hz` |
| Stiffness entry | newton per metre | `N/m` |
| Time, time step | second | `s` |
| Velocity, acceleration | metre per second (squared) | `m/s`, `m/s^2` |
| Rayleigh damping `a` (of `M`), `b` (of `K`) | 1 / second, second | `1/s`, `s` |
| Phase of a harmonic amplitude | degree (radian in an amplitude's `phase`) | `deg` |

The MBB benchmark deliberately uses `E = 1 Pa`, a 1 m cell and `F = 1 N` so that
its compliance can be read against the dimensionless values usually quoted for
that problem. It is still SI - just a particular choice of scale.

## Coordinate system

Right-handed Cartesian, with the model in the x-y plane and z out of plane:

```
      y
      ^
      |            thickness t extends +- t/2 in z
      |
      +-------> x
     /
    z  (out of the page, towards the reader)
```

* The structured mesh generator places the domain with its lower-left corner at
  `(x0, y0)`, default `(0, 0)`.
* Node numbering on a structured grid is `node(i, j) = j * (nx + 1) + i`, so the
  x index runs fastest.
* Element numbering is `elem(i, j) = j * nx + i`, again x fastest. This is why a
  C-order reshape of an element vector to `(ny, nx)` produces an image whose
  row 0 is the *bottom* of the domain.
* Element nodes are stored counter-clockwise from the lower-left corner:

```
  3 --------- 2          eta
  |           |           ^
  |    (e)    |           |
  |           |           +---> xi
  0 --------- 1
```

Counter-clockwise ordering makes the Jacobian determinant positive; a clockwise
or self-intersecting element is rejected by `Mesh::validate()` with a message
naming the element.

### Solid meshes

A `structured_hex` mesh is a box with its lower corner at `(x0, y0, z0)`.
Node numbering is `node(i, j, k) = k (nx+1)(ny+1) + j (nx+1) + i` and element
numbering `elem(i, j, k) = k nx ny + j nx + i`, x fastest, then y, then z. The
eight nodes of a Hex8 follow the VTK hexahedron convention: the bottom face
(`zeta = -1`) counter-clockwise seen from above, then the top face in the same
order:

```
        7 --------- 6
       /|          /|          zeta
      4 --------- 5 |           ^   eta
      | |         | |           |  /
      | 3 --------|-2           | /
      |/          |/            +-----> xi
      0 --------- 1
```

The six faces are numbered `0: bottom (0 3 2 1)`, `1: top (4 5 6 7)`,
`2: front y = -1 (0 1 5 4)`, `3: right x = +1 (1 2 6 5)`,
`4: back y = +1 (2 3 7 6)`, `5: left x = -1 (3 0 4 7)`, each wound so its
right-hand normal points out of the element. That winding is what the STL
export relies on for outward normals.

### Simplices

A Tri3 lists its three nodes counter-clockwise (positive area), with edges
`0: (0 1)`, `1: (1 2)`, `2: (2 0)`. A Tet4 lists its four nodes so that
`(x1 - x0) . ((x2 - x0) x (x3 - x0)) > 0` - positive volume, the VTK
convention - with faces `0: (0 2 1)`, `1: (0 1 3)`, `2: (1 2 3)`,
`3: (0 3 2)`, opposite nodes 3, 2, 0 and 1 and wound outward like the
hexahedron's. `structured_tri` and `structured_tet` meshes keep the node
numbering of the grid they split: element `2 elem(i, j) + s` is triangle `s`
of cell `(i, j)`, and element `6 elem(i, j, k) + s` tetrahedron `s` of the
cell.

### Quadratic tetrahedra

A Tet10 lists its four corners as a Tet4 does (positive volume), then the
six edge nodes of the edges `(0 1)`, `(1 2)`, `(2 0)`, `(0 3)`, `(1 3)`,
`(2 3)` - the VTK (`VTK_QUADRATIC_TETRA`, cell type 24), Abaqus / CalculiX
(`C3D10`) and scikit-fem order. Gmsh (element type 11) numbers the last two
the other way round, and its reader swaps them. A face lists its three
corners as the Tet4 face does, then the edge nodes of its corner pairs
`(0 1)`, `(1 2)`, `(2 0)`; faces are matched between cells by their corners
only. `mesh.order: 2` elevates a Tet4 mesh by appending one node per edge
after the existing nodes, in order of first appearance, so the corner nodes
keep their numbers.

### Meshes read from a file

Nodes keep the order of the file, renumbered consecutively from 0 after the
nodes no cell uses are dropped; cells keep the order of the file. A cell in
the opposite orientation to the one above (a clockwise triangle, a mirrored
tetrahedron or hexahedron) is re-ordered as it is read, and the number of
such cells is recorded under `mesh.file.cells_reoriented`.

## Degrees of freedom

`dim` translations per node, numbered node-major:

```
  dof(node n, component c) = dim * n + c,     c = 0 -> u_x,  1 -> u_y,  2 -> u_z
```

An element's DOF vector follows the same ordering:
`{u_x1, u_y1, u_x2, u_y2, ...}` on a plane mesh and `{u_x1, u_y1, u_z1, ...}`
on a solid one.

There are no rotational degrees of freedom: these are continuum elements, not
beam or shell elements.

## Signs

* Displacements and forces are positive along `+x` and `+y`. A downward load is
  a negative `y` force.
* Tractions are given as a stress vector in **global** components, not as
  pressure normal to a surface. A downward traction on a horizontal top edge is
  `[0, -p]`; to push *into* a surface you must supply the correct sign yourself.
  This is deliberate: a "pressure" convention needs an outward normal, which is
  ambiguous on a corner.
* Tensile stress is positive. `sigma_xx > 0` means the material is in tension
  along x.
* Reactions are reported as the force the *supports apply to the structure*, so
  applied load plus reaction sums to zero. Every static result records that sum
  and its relative error.
* Positive moment is counter-clockwise about `+z`, i.e. `M = x F_y - y F_x`,
  taken about the origin. On a solid mesh the moment is the full vector
  `M = x cross F` about the origin, reported component by component.
* A non-linear load factor `lambda` multiplies every load of the case at
  once - forces, pressures, body loads, the rotation's centrifugal load (so
  the angular velocity grows as `sqrt(lambda)`), the temperature change and
  prescribed displacements - and `lambda = 1` is the load case as written.
  A `load_path` takes `lambda` through its turning points in order, and a
  negative `lambda` reverses every load. A reaction monitor sums the support
  reactions over its nodes with the sign above; a displacement monitor
  averages the nodal displacements.
* Plasticity: the equivalent (accumulated) plastic strain is
  `a = int sqrt(2/3) ||d eps_p||`, which in uniaxial tension is the plastic
  strain itself; the yield stress `sigma_y(a)` and the hardening moduli are
  uniaxial values. The plastic strain is stored with engineering shears
  (as every strain), the back stress with tensorial components (as every
  stress).
* A buckling load factor `lambda` multiplies the load case: the structure is
  predicted to buckle under `lambda f`. Load factors are reported positive
  and ascending; the negative ones of the pencil (buckling under the
  reversed load) are not reported, and a load case with no positive one is
  flagged `no_positive_load_factor`. Buckling modes are normalised to
  `phi^T K phi = 1` (and in the VTK files to a largest displacement of 1);
  their sign is arbitrary.
* A transient's amplitude `A(t)` multiplies every load and prescribed
  displacement of the case at time `t`; step 0 of its history is the
  initial state. A reaction monitor sums `M a + C v + K u - A f` over the
  prescribed DOFs of its nodes - the force the supports (or the shaking
  device) apply - so the force needed to accelerate the structure shows in
  it.
* A harmonic response is `Re(U e^{i omega t})` under the load `f cos(omega t)`:
  its phase `arg U` is measured against the load, `0` in phase and `180 deg`
  in antiphase, so the displacement under a downward load starts at `180 deg`
  at low frequency. A velocity amplitude is `i omega U`, an acceleration
  `-omega^2 U`; "the largest displacement over a cycle" is the semi-major
  axis of the ellipse a node traces.
* A build direction `+y` means the part grows along `+y` from a plate at the
  low-`y` end of the domain; `-y` from a plate at the high end. Layer 0 is
  the layer on the plate.

## Voigt notation

Stress and strain are stored as vectors with the *engineering* shear strain:

```
  plane:  sigma = { sigma_xx, sigma_yy, sigma_xy }
          eps   = { eps_xx,   eps_yy,   gamma_xy },    gamma_xy = 2 eps_xy

  solid:  sigma = { sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_yz, sigma_zx }
          eps   = { eps_xx,   eps_yy,   eps_zz,   gamma_xy, gamma_yz, gamma_zx }
```

This pairing makes `sigma^T eps` the strain-energy density, which is why the
constitutive matrix carries `G` rather than `2G` in its shear entries. The CSV
columns follow the same order (`sxx, syy, szz, sxy, syz, szx`).

## Energy definitions

| Quantity | Definition | Note |
|----------|------------|------|
| Strain energy | `U = 1/2 u^T K u` | always the elastic energy stored |
| Compliance | `C = f^T u` | `f` is the **applied** load vector, excluding reactions |

For homogeneous Dirichlet data these satisfy `C = 2U` exactly, and every static
result reports the ratio `C / 2U` so the identity can be checked. With a
non-zero prescribed displacement the two differ: the patch test, for example,
has `C = 0` and `U > 0`, because all the work is done by the supports.

A transient reports the kinetic energy `T = 1/2 v^T M v`, the strain energy
`U` above (with plasticity the stored energy, elastic plus hardening), the
energy dissipated by damping `D = sum dt/4 (v_n + v_n+1)^T C (v_n + v_n+1)`
and the work of the loads and reactions
`W = sum 1/2 (u_n+1 - u_n)^T (F_n + F_n+1)`, and the balance
`E_0 + W - T - U - D` - zero to round-off for the trapezoidal rule on a
linear model, the numerical dissipation of HHT-alpha, and with plasticity
the plastic dissipation.

The topology objective is the weight-normalised sum of the per-case
compliances,

```
  c(x) = sum_l w_l f_l^T u_l ,     sum_l w_l = 1
```

so a run's objective value does not change if every weight is scaled by the same
factor.

## Density and volume fraction

* A design variable `x_e` and a physical density `rho_e` both live in `[0, 1]`
  and are dimensionless.
* The volume fraction is `sum_e rho_e v_e / sum_e v_e`, with `v_e` the element
  volume in `m^3`. It is computed on the **physical** density - filtered, and
  projected when the Heaviside projection is on - which is the quantity the
  constraint is written on.
* Passive solid regions count towards the volume budget. A mass constraint on a
  real part includes its attachment bosses, so the budget must too.
* Mass is `rho_material * sum_e rho_e v_e`, i.e. the SIMP density multiplies the
  material density linearly.

## Interpreting a density field

A SIMP density field is a *relative material distribution*, not a solid body.
SparLab never presents one as geometry without saying how it was converted:

* every figure of a density field states that black is `rho = 1`;
* the "interpretation" block in a topology summary records the threshold, how
  many elements survived it, how many disconnected groups they formed, and how
  much material was discarded as islands;
* the modal analysis of an "optimised structure" is run on that explicitly
  extracted sub-mesh, not on the density field.

See `docs/limitations.md` for what this interpretation does and does not claim.

## Tolerances

Every tolerance is configurable and every run records the value it used in
`summary.json` under `tolerances`. The defaults are:

| Tolerance | Default | What it bounds |
|-----------|---------|----------------|
| `solver.linear.residual_tolerance` | `1e-8` | scaled residual `\|\|Ku-f\|\| / \|\|f\|\|` after each solve, whatever the solver; a solve whose backward error `\|\|Ku-f\|\| / \|\| \|K\|\|u\| + \|f\| \|\|` is at most `64 eps` is accepted too (its residual is round-off; `docs/formulation.md`, section 5) |
| `solver.linear.pivot_tolerance` | `1e-14` | smallest / largest LDL^T pivot before the system is called singular |
| `solver.linear.iterative_tolerance` | `1e-12` | relative residual the CG solvers (Jacobi and multigrid) iterate to |
| `solver.linear.amg.coarse_pivot_tolerance` | `1e-13` | smallest / largest pivot of the multigrid's coarsest factorisation before the model is called under-constrained |
| `solver.equilibrium_tolerance` | `1e-6` | relative global force-balance error |
| `modal.tolerance` | `1e-10` | relative eigenvalue change between subspace iterations |
| `modal.residual_tolerance` | `1e-6` | relative eigenpair residual, also part of the stopping rule |
| `topology.optimizer.change_tolerance` | `1e-2` | `max \|dx\|` between iterations |
| `topology.optimizer.objective_tolerance` | `0` (off; `5e-5` in most decks) | relative spread `(max - min) / \|c_k\|` of the last `objective_window + 1` compliances |
| `topology.optimizer.volume_tolerance` | `1e-10` | relative volume error of the multiplier bisection (OC) |
| `topology.optimizer.constraint_tolerance` | `1e-4` | largest MMA constraint value an iterate may have and still count as feasible; convergence requires it |
| `topology.optimizer.mma.subproblem_tolerance` | `1e-7` | final barrier parameter of the MMA subproblem; each level converges to a KKT residual of 0.9 times its barrier |
| `topology.optimizer.mma.max_newton_iterations` | `500` | Newton iterations per barrier level before the subproblem is reported as failed |
| `topology.stress.feasibility_tolerance` | `1e-3` | relative margin used when the summary reports whether the relaxed stress maximum meets the limit |
| `buckling.tolerance` (and `topology.buckling_constraint.tolerance`) | `1e-8` | relative change of the requested load factors between subspace iterations |
| `buckling.residual_tolerance` | `1e-6` | relative eigenpair residual `\|\|K phi + lambda K_G phi\|\| / \|\|K phi\|\|`, also part of the stopping rule |
| `topology.length_scale_check.tolerance` | `0.02` | share of the solid (void) volume a morphological probe may change and still pass |
| `nonlinear.residual_tolerance` | `1e-8` | out-of-balance force of a converged step over the load scale (the largest of the applied loads, the reactions and the thermal forces); a residual at its round-off floor, `max(1024 eps` gross element forces`, 64 eps \|\| \|K_T\|\|u\| \|\|)`, is accepted too while that floor is below `1e-6` of the scale |
| J2 return, Voce hardening | `1e-13` | the consistency residual of the plastic multiplier over the yield radius; 60 Newton iterations, else an error |
| J2 return, plane stress | `1e-12` | `\|sigma_33\|` over the larger of the stress and the yield stress; 30 Newton iterations, else an error |
| elastic return at the yield surface | `1e-12` | a trial overshoot of the yield radius below this is elastic (round-off must not start plastic flow); a point that yielded in its last step and sits within `1e-10` of its surface keeps the continuum tangent |
| `nonlinear.displacement_tolerance` | `1e-8` | last Newton correction over the displacement increment of the step, no finer than `64 eps \|\|u\|\|` |
| `nonlinear.max_cuts` | `12` | successive halvings of a failing step, and of the steps closing in on a critical point, before the run stops |
| `transient.residual_tolerance`, `transient.displacement_tolerance` | `1e-8` | non-linear transient: the HHT-alpha residual over the largest force involved (inertia, loads, damping, reactions), and the last Newton correction over the larger of the step's increment and the displacement, no finer than `64 eps \|\|u\|\|`; the residual's round-off floor is accepted as in the statics |
| `transient.max_cuts` | `8` | successive halvings of a transient step whose Newton iteration fails |
| harmonic solve backward error | `1e-10` | `\|\|A U - f\|\| / \|\| \|A\|\|U\| + \|f\| \|\|` of the complex solve at every frequency; above it the frequency is too close to an undamped natural frequency and the run stops |
| harmonic amplification | `1e6` | a response this many times the static one is flagged as round-off near a natural frequency |
| non-linear transient energy balance | `1e-2` | energy created (or, without plasticity, lost) beyond this share of the energies involved draws a warning (trapezoidal rule only) |
| `forming.steps[].explicit.energy_tolerance` | `0.05` | explicit forming step: the energy balance error, relative to the largest energy of the run so far, above which the run warns (`summary.json`: `analysis.explicit_tolerances`) |
| `forming.steps[].explicit.energy_limit` | `0.5` | explicit forming step: the balance error that stops the run (an unstable integration) |
| `forming.steps[].explicit.kinetic_ratio_warning` | `0.1` | explicit forming step: kinetic over internal energy after the first contact above which the run is warned about as not quasi-static |
| `forming.steps[].explicit.mass_scaling.max_added_mass_fraction` | `0.05` | explicit forming step: added over physical mass above which mass scaling is warned about |
| explicit Hex8 kernel self-check | `1e-10` | relative difference between the dedicated kernel and the element dispatch at the start of a run above which the dispatch is used instead |

The decks that use the multigrid solver (`bracket_3d_projected`,
`bracket_3d_large`, `engine_mount_3d`) set `iterative_tolerance` to `1e-10`.
The multigrid verification study runs at that tolerance and measures how far
the displacement then sits from the direct solution (`docs/verification.md`).

Exceeding a solver tolerance raises an exception with a diagnosis; exceeding an
optimiser tolerance is reported as non-convergence in both the console output
and `summary.json`. Nothing is silently accepted.

## Determinism

The only randomness in the code is the filler part of the subspace-iteration
starting bases (modal and buckling) and the node perturbation of the
`make_perturbed_*_mesh` generators (quadrilateral, hexahedral, and the
triangle and tetrahedron meshes split from them, linear and quadratic). All
take an explicit seed (`modal.seed` and `buckling.seed`, default `20240917`)
and default to deterministic values, so repeated runs of the same deck on
the same build give bit-identical results - the benchmark comparison runs
reproduced their earlier numbers to every printed digit. The MMA subproblem solver
and the stress-constraint scaling are deterministic by construction: they
start from fixed values and involve no random choice.

Threads do not change results either. The multigrid solver is SparLab's
only OpenMP code: its parallel kernels write disjoint rows, and every inner
product sums fixed-size chunks in a fixed order, so a solve gives the same
bits on one thread as on many (a test asserts it). Eigen threads one more
kernel when OpenMP is on, the matrix-vector product inside its
Jacobi-preconditioned CG, and computes each row of it on a single thread,
so that solver is deterministic too. The element loops -
assembly, sensitivities, stress recovery - run sequentially in element
order. The committed meshes come from Gmsh run single-threaded with fixed
options by `python/scripts/make_meshes.py`, which records the Gmsh version
and every option beside each mesh; the Tet10 study regenerates the 4 mm
engine mount identically (39 936 cells, 8 640 nodes).
