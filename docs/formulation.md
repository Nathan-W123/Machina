# Mathematical formulation

This document states exactly what SparLab solves. Notation follows
`docs/conventions.md`; every equation below is implemented in the file named
beside it.

## 1. Continuum problem

Small-strain linear elasticity on a domain `Omega` - a plane of constant
out-of-plane thickness `t` (2-D) or a solid (3-D, where `t = 1` everywhere
below) - with displacement boundary `Gamma_u`, traction boundary `Gamma_t` and
pressure boundary `Gamma_p`:

```
  div sigma + b = 0            in Omega          (equilibrium)
  eps = 1/2 (grad u + grad u^T)                  (kinematics)
  sigma = D : (eps - eps_0)                      (constitutive)
  u = u_bar                    on Gamma_u
  sigma . n = t_bar            on Gamma_t
  sigma . n = -p n             on Gamma_p
```

`b` is the body force density - self-weight `rho g`, a uniform force density,
or the centrifugal load `rho omega^2 r_perp` of a steady rotation - and
`eps_0 = alpha (T - T_ref)` (per idealisation, below) the free thermal strain of
a temperature field. The weak form, for all admissible `v`, is

```
  integral_Omega  t eps(v)^T D eps(u) dOmega
      =  integral_Omega t v^T b dOmega  +  integral_Omega t eps(v)^T D eps_0 dOmega
       + integral_Gamma_t  t v^T t_bar dGamma  -  integral_Gamma_p  t p v^T n dGamma
       + sum_k v(x_k)^T F_k
```

with the last term the point loads. Each element may carry its own material
(`material_regions`), so `D`, `rho`, `alpha` and `T_ref` are per element.
Section 4b gives the discretisation of the load terms.

### Plane idealisations

With Voigt ordering `{sigma_xx, sigma_yy, sigma_xy}` and
`{eps_xx, eps_yy, gamma_xy}` (`gamma_xy = 2 eps_xy`), an isotropic material
gives (`material/IsotropicMaterial.cpp`):

**Plane stress** (`sigma_zz = 0`; a thin sheet loaded in its plane)

```
           E     | 1    nu       0     |
  D  =  -------- | nu    1       0     |
        1 - nu^2 | 0     0  (1-nu)/2   |
```

**Plane strain** (`eps_zz = 0`; a long prismatic body, or a plane of symmetry)

```
                 E          | 1-nu   nu        0      |
  D  =  ------------------   | nu    1-nu       0      |
        (1+nu)(1-2nu)        | 0      0    (1-2nu)/2   |
```

In plane strain the out-of-plane stress is `sigma_zz = nu (sigma_xx + sigma_yy)`
and the von Mises formula must use it; `fem/StressRecovery.cpp` does.

The `(3,3)` entry of both matrices equals the shear modulus
`G = E / (2(1+nu))`, which is a consequence of pairing Voigt stress with
*engineering* shear strain and is checked in the test suite.

Both idealisations are implemented and unit-tested. Plane stress is the default
on a plane mesh and is the path the plane benchmark and verification studies
exercise.

### Three dimensions

On a solid (`structured_hex`) mesh the full isotropic law is used, with Voigt
ordering `{sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_yz, sigma_zx}` and
the three engineering shear strains. In terms of the Lame constants
`lambda = E nu / ((1+nu)(1-2nu))` and `G = E / (2(1+nu))`,

```
        | lambda+2G   lambda     lambda    0  0  0 |
        | lambda    lambda+2G    lambda    0  0  0 |
  D  =  | lambda      lambda   lambda+2G   0  0  0 |          (6 x 6)
        |    0          0         0        G  0  0 |
        |    0          0         0        0  G  0 |
        |    0          0         0        0  0  G |
```

which the test suite checks entry by entry against this Lame form, and whose
von Mises reduction it checks against uniaxial, hydrostatic (zero) and
pure-shear (`sqrt(3) tau`) states. There is no thickness: a solid deck rejects
`model.thickness`.

## 2. Element formulation

### Shape functions

The four-node bilinear isoparametric quadrilateral (`elements/Quad4.cpp`), on
the reference square `(xi, eta) in [-1,1]^2`:

```
  N_1 = 1/4 (1 - xi)(1 - eta)        N_3 = 1/4 (1 + xi)(1 + eta)
  N_2 = 1/4 (1 + xi)(1 - eta)        N_4 = 1/4 (1 - xi)(1 + eta)
```

with natural derivatives

```
  dN_1/dxi  = -1/4 (1 - eta)     dN_1/deta = -1/4 (1 - xi)
  dN_2/dxi  = +1/4 (1 - eta)     dN_2/deta = -1/4 (1 + xi)
  dN_3/dxi  = +1/4 (1 + eta)     dN_3/deta = +1/4 (1 + xi)
  dN_4/dxi  = -1/4 (1 + eta)     dN_4/deta = +1/4 (1 - xi)
```

These satisfy `sum_a N_a = 1` everywhere, `N_a = delta_ab` at node `b`, and
`sum_a grad N_a = 0` - all three are unit-tested.

### Isoparametric mapping

Geometry and displacement use the same basis:

```
  x(xi, eta) = sum_a N_a(xi, eta) x_a        u(xi, eta) = sum_a N_a(xi, eta) u_a
```

The Jacobian and its determinant are

```
  J_ij = dx_i / dxi_j = sum_a x_a,i  dN_a/dxi_j            (a 2 x 2 matrix)
  detJ = J_11 J_22 - J_12 J_21
```

`detJ <= 0` at any quadrature point means the element is inverted, collapsed or
wound clockwise, and raises `MeshError` naming the element and the parametric
point. Physical shape-function gradients follow from
`dN_a/dx = (dN_a/dxi) J^{-1}` (the 2x2 inverse is written out explicitly; the
assembly loop runs it millions of times).

### Strain-displacement operator

```
          | dN_1/dx    0     dN_2/dx    0     ...  |
  B  =    |   0     dN_1/dy     0    dN_2/dy  ...  |      (3 x 8)
          | dN_1/dy dN_1/dx  dN_2/dy dN_2/dx  ...  |
```

so `eps = B u_e` with `u_e` in node-major order. Because the Q4 space contains
all linear fields exactly, `B u_e` reproduces a constant strain field exactly at
*every* parametric point, on distorted meshes as well as uniform ones. That is
what the patch test checks, and it holds to round-off.

### Quadrature

Gauss-Legendre tensor-product rules on the reference square
(`elements/Quadrature.cpp`). An `n`-point rule per direction integrates
polynomials of degree `2n - 1` exactly.

| Kernel | Default rule | Why |
|--------|--------------|-----|
| Stiffness `K_e` | 2 x 2 (2 x 2 x 2 for the Hex8) | exact for the Q4 on a parallelogram and the Hex8 on a parallelepiped; full integration, no reduced-integration hourglass modes |
| Mass `M_e` | 3 x 3 (3 x 3 x 3) | `N^T N` is biquadratic, and on a general cell `detJ` is not constant |
| Edge / face traction | 2 points per direction | exact for a linear traction times a linear shape function |

All are configurable per run via `model.integration`. On a rectangle the 2x2,
3x3 and 4x4 stiffness matrices agree to round-off, which is unit-tested; the
same holds for the Hex8 on a box.

### Element matrices

```
  K_e = integral_Omega_e  t B^T D B  dOmega
      = sum_g  w_g  t  detJ(xi_g)  B(xi_g)^T D B(xi_g)              [N/m]

  M_e = integral_Omega_e  rho t N^T N dOmega
      = sum_g  w_g  rho t detJ(xi_g)  N(xi_g)^T N(xi_g)             [kg]
```

where `N` is the `2 x 8` shape-function matrix. Both are symmetrised after
integration, since any asymmetry is pure round-off.

`K_e` has exactly three zero eigenvalues - the two translations and the
infinitesimal rotation - and five strictly positive ones (six and eighteen
for the Hex8). Both parts of that statement are unit-tested on a uniform and
on a distorted element of each type.

Row-sum lumping of `M_e` gives the optional diagonal mass matrix. For the Q4 the
row sums add up to the exact element mass, so total mass is conserved either
way.

### Consistent edge loads

For an edge with end nodes `a`, `b` and a constant traction `t_bar`,

```
  f_e = integral_Gamma_e  t N^T t_bar ds
```

parametrised by `s in [-1,1]` with the constant Jacobian `|x_b - x_a| / 2`. For
a constant traction on a straight edge this splits the resultant evenly between
the two nodes, and the total is exactly `t_bar * (edge length) * t`. The test
suite checks the total is mesh independent under refinement.

### The trilinear hexahedron

On a solid mesh the element is the eight-node trilinear hexahedron
(`elements/Hex8.cpp`) on the reference cube `(xi, eta, zeta) in [-1,1]^3`,
nodes in the VTK order of `docs/conventions.md`:

```
  N_a = 1/8 (1 + xi_a xi)(1 + eta_a eta)(1 + zeta_a zeta),   (xi_a, eta_a, zeta_a) = +-1
```

The isoparametric mapping, the `3 x 3` Jacobian (inverted explicitly through
its adjugate; `detJ <= 0` raises `MeshError` naming the element and the
point) and the `6 x 24` strain-displacement operator

```
          | dN_a/dx     0        0     |
          |   0      dN_a/dy     0     |
  B_a  =  |   0         0     dN_a/dz  |        columns of node a, rows in the
          | dN_a/dy  dN_a/dx     0     |        order xx, yy, zz, xy, yz, zx
          |   0      dN_a/dz  dN_a/dy  |
          | dN_a/dz     0     dN_a/dx  |
```

follow the Q4 pattern exactly, and so do the checks: `sum_a N_a = 1`,
`sum_a grad N_a = 0`, exact reproduction of a linear field with six
independent constant strains on a distorted cell, `sum_g w_g detJ` against the
volume of a sheared box, `K_e` symmetric with exactly **six** zero eigenvalues
(three translations, three rotations) and the rest positive, and the 2x2x2
and 3x3x3 rules agreeing on a box. Face tractions integrate the
bilinear face `x(s, t)` with the area element `|x_s cross x_t|`, so a
constant traction on a flat face splits its resultant evenly over the four
corners and the total is exactly `t_bar * area` at any resolution.

The fully integrated Hex8 shares the Q4's stiffness in bending: it needs
several elements through a bending depth, which the 3-D mesh-convergence
study quantifies (`docs/verification.md`), or the incompatible modes below.

### The Hex8 with incompatible modes

With `model.element_formulation: "incompatible_modes"` a Hex8 mesh takes the
element of Wilson, Taylor, Doherty and Ghaboussi (1973, *Incompatible
displacement models*, in *Numerical and Computer Methods in Structural
Mechanics*) with the correction of Taylor, Beresford and Wilson (1976, *A
non-conforming element for stress analysis*, IJNME 10)
(`elements/Hex8Incompatible.cpp`, `fem/IncompatibleModes.cpp`). The
trilinear element cannot bend: a pure-bending field puts a parasitic shear
`gamma_xz ~ x` into every cell, because its sides stay straight, so an
element long against its depth is too stiff by a factor that grows as
`(L/h)^2`. That is what spoils a sheet meshed with one or two elements
through its thickness, and with it the springback, which is the elastic
unloading of a bending moment.

**The modes.** Three bubble fields `phi_k = 1 - xi_k^2` per displacement
component (nine parameters `alpha` per element), which vanish at the nodes
and are not continuous between elements. Their gradients are taken with the
Jacobian `J_0` of the element centre, scaled by the ratio of determinants:

```
  g~_k(xi) = (j_0 / j(xi)) (-2 xi_k) J_0^-T e_k            (Taylor's correction)
  H = U G^T + A G~^T                                        (A: 3 x 3, column m = pseudo-node m)
```

so `sum_q w_q j_q g~_k = 0` for every rule symmetric in each direction: a
constant strain puts no load on the modes, `alpha = 0`, and the element
passes the patch test on any geometry (Wilson's original fails it on a
distorted cell). The modes enter every kinematic operator as three extra
"pseudo-nodes": the element is an eleven-node element of 33 degrees of
freedom `(u, alpha)`, and the linear `B`, the Green-Lagrange
`B_NL(F, [G | G~])` and the geometric stiffness `([G|G~]^T S [G|G~]) x I` are
those of the library applied to eleven gradient columns. In small strain the
strain span of the modes is the enhanced field EAS-9 (Andelfinger and Ramm
1993), and the condensed stiffness equals Simo and Rifai's (1990). In finite
deformation this is the enhanced deformation gradient `F = I + H` of Simo
and Armero (1992, *Geometrically non-linear enhanced strain mixed methods
and the method of incompatible modes*, IJNME 33) and Simo, Armero and Taylor
(1993, CMAME 110); it is objective (a rigid turn `R` is carried by `R A`).

**Linear analyses** condense in closed form, with fixed-size 24 x 24,
24 x 9 and 9 x 9 blocks:

```
  K* = K_uu - K_ua K_aa^-1 K_au,     K_aa = int B~^T D B~ dV,   K_au = int B~^T D B dV
  B^ = B - B~ K_aa^-1 K_au          (int B^^T D B^ dV = K*)
```

`strain_operator` stays the compatible `B` (its normal rows are read as
shape-function gradients elsewhere); the linear consumers take the
condensed `B^`: stress recovery (whose modes also take their share of a
non-uniform thermal strain), the thermal load `f* = int B^^T D eps_0 dV` and
its self energy, and
the geometric stiffness of linear buckling, which takes the stress of `B^`
on the compatible gradients (the modes are left out of the initial-stress
term, a linearised-buckling approximation). The mass and the surface loads
are those of the Hex8: the modes carry no inertia and no load.

**Non-linear analyses** (`nonlinear`, `forming`; small strain, `finite`,
`finite_logarithmic`; elastic, neo-Hookean, J2, Hill48, Chaboche). The modes
are internal, `r_alpha(u, alpha) = 0`, and on every evaluation of an element
the parameters are found by Newton's method from the parameters committed
at the last converged step, with every point's return from its committed
state:

```
  alpha_0 = alpha_n,      alpha_{i+1} = alpha_i + s_i d_i,      d_i = -K_aa^-1 r_alpha
  f*  = f_u - K_ua K_aa^-1 r_alpha        (the last correction: the error of f* is quadratic)
  K*  = K_uu - K_ua K_aa^-1 K_au          (K_au, not K_ua^T, for a non-symmetric tangent)
  q*  = q_u - K_ua K_aa^-1 q_a            (the thermal load rate)
```

so an element is a pure function of `(u, lambda)` and the committed
history, which the line search, the step cuts, the arc length and the
bitwise-deterministic parallel evaluation rely on (the Simo-Rifai update of
`alpha` with the global iterate would not be). It stops at
`|r_alpha|_inf <= 1e-8` of the gross scale `max_i sum_q w_q |B~_q|^T |sigma_q|`
plus a round-off floor, and raises `SolverError`, which cuts the step, after
25 iterations or at a singular `K_aa`. The full step (`s_i = 1`) is taken
whenever it reduces `|r_alpha|`; far from the solution - a return far from
its committed state, a load reversal - the points switch between elastic
and plastic along the step, `r_alpha` is only piecewise smooth and Newton's
iterates can cycle, so `s_i` then comes from a regula-falsi search for the
root of `d_i^T r_alpha(alpha_i + s d_i)`, the derivative of the incremental
potential along `d_i` wherever there is one (Crisfield 1991, *Non-linear
finite element analysis of solids and structures*, vol. 1, sec. 9.3); a
trial that inverts a point counts as beyond the root, and where that
derivative brackets no root on the step (`K_aa` not positive definite, a
non-convex potential) the step is halved on `|r_alpha|`, the best trial
kept. Elastic small strain
is linear in `alpha`: one iteration. Plastic steps converge quadratically;
in the unit tests a converged step takes at most 3 local iterations, a
reversal to six times the strain in one increment 3 to 12. Everything runs
on fixed-size Eigen types (6 x 33 operators, the 9 x 9, 24 x 9, 24 x 24
blocks) in a per-thread workspace: the local iteration allocates nothing.

`NonlinearSystem` keeps the committed `alpha` of every element, elastic ones
too, and commits it with the plastic states; the stress output evaluates
the points at the committed parameters (the converged ones), so element
stresses, strains and Jacobians are those of the enhanced `F = I + H`. The
forming driver carries the parameters in its restartable `AnalysisState`
(`internal`, per element); a state of the other formulation is refused. The
committed parameters are only the start of the local iteration, not
history: a restart without them finds the same state to the local
tolerance, from further away.

**Mean dilatation** (E-bar) is not combined with the modes: `all` is
refused and `auto` does not apply it to the incompatible-mode Hex8.
Averaging the dilatation of the augmented element leaves a zero-energy
mode (the field `u = k (x z, y z, (z^2 - x^2 - y^2)/2)`, whose strain `k z I`
averages to nothing), and averaging only the compatible dilatation leaves
a sheet one element thick with a fraction of its bending stiffness. The
modes relax the isochoric constraint of plastic flow themselves: they
absorb the linear part of every point's dilatation, which is the part that
locks.

**Hourglassing in compression.** The enhanced elements of finite
deformation can lose stability under large compression although the
material and the compatible element are stable (Wriggers and Reese 1996,
*A note on enhanced strain methods for large deformations*, CMAME 135). On
a neo-Hookean unit cube (`nu = 0.3`) in homogeneous uniaxial compression,
the condensed tangent of the compression test loses definiteness between
the stretches 0.66 and 0.65 in a spurious hourglass mode, while the
standard Hex8's stays positive definite to 0.3 (unit test). The thinning of
a sheet under equibiaxial stretch - tensile stress - stays stable to a
thickness stretch below 0.5. Springback (small elastic strains) is far from
it; a formed part in large in-plane compression (a flange, a wrinkle) is
not, and should be checked against the standard element there.

**Limits.** The Taylor-corrected modes are exact in pure bending on
affine cells (boxes, parallelepipeds) and lose much of their effect on
trapezoidal ones: the MacNeal-Harder beam on trapezoidal cells gives 5 % of
the reference deflection in plane and 3 % out of plane (MacNeal 1987, *A
theorem regarding the locking of tapered four-noded membrane elements*,
IJNME 24); curved walls of cones and pyramids are of that kind. The export
to CalculiX writes `C3D8I`, which is not cross-validated here (no CalculiX
on the machine).

### A sheet's rule through its thickness

`model.integration.thickness_points` (1 to 7) with
`thickness_direction` (`x`, `y` or `z`, default `z`) gives every Hex8,
standard or with incompatible modes, the `n x n x m` Gauss rule
(`gauss_legendre_box`; `n = stiffness_points`, `m = thickness_points` along
the natural axis that the structured generator aligns with that global
axis; a mesh whose elements have another axis across the sheet is warned
about). The elastic-plastic stress through a bent sheet is kinked where
the plastic zones begin, which two Gauss points a layer integrate badly:
springback is sensitive to it (Wagoner and Li 2007, *Simulation of
springback: through-thickness integration*, Int. J. Plasticity 23). The
rules of 5 to 7 points come from Newton's iteration on the Legendre
recurrence; equal orders give the cube rule point for point and bit for
bit. The history of the non-linear analysis (one plastic state per point),
the stress-evaluation points and every consumer follow the rule, in the
point order of the other rules (`zeta` the slowest index, `xi` the fastest).
The integration of the moment of a rigid-plastic shell outside an elastic
core of a quarter of the half-thickness (bending to `4 k_y`) is short by
3.5 % with 5 points and 0.9 % with 7; the springback study measures 3.2 %
and 0.6 % (`docs/verification.md`, section 27). The CalculiX export refuses
the rule (CalculiX's `C3D8`/`C3D8I` have the fixed 2 x 2 x 2 rule).

### The linear simplices

Meshes read from a mesh generator are mostly triangles and tetrahedra, so
SparLab also has the three-node triangle (`elements/Tri3.cpp`, plane) and
the four-node tetrahedron (`elements/Tet4.cpp`, solid). Their shape functions
are the barycentric coordinates, linear in `x`:

```
  Tri3:  N_a = (alpha_a + beta_a x + gamma_a y) / (2 A)
  Tet4:  N_a = (a_a + b_a x + c_a y + d_a z) / (6 V)
```

so `grad N_a` - and with it `B`, the strain and the stress - is constant over
the cell, and the element matrices have closed forms:

```
  K_e = t A B^T D B            (Tri3, B 3 x 6)
  K_e = V B^T D B              (Tet4, B 6 x 12)

  M_e = rho t A / 12 (1 + delta_ab) I_2      (Tri3, node blocks)
  M_e = rho V / 20 (1 + delta_ab) I_3        (Tet4, node blocks)
```

which is the exact consistent mass. A traction on a straight edge (a flat
triangular face) splits its resultant equally over the two (three) nodes. A
cell with non-positive area or volume raises `MeshError` naming the element
and its measure; the file readers re-order mirrored cells before that point,
so the error means a folded cell, not a node-order convention. Stresses are
reported at the centroid, which is exact for a constant-strain element.

The structured generators split each Q4 into two triangles, with the
diagonal alternating from cell to cell so there is no preferred direction,
and each Hex8 into the six Kuhn tetrahedra around its `0-6` diagonal, which
cuts every face of the grid along the same diagonal from both sides and so
stays conforming. The same box then exists as Q4, Tri3, Hex8 and Tet4 meshes
with identical nodes, which is what the simplex verification studies use.

Constant-strain elements lock in bending: a linear field cannot represent
the linearly varying strain through a beam's depth. They reproduce every
constant-strain state exactly (the patch test) but converge on a bending
problem more slowly per unknown than the bilinear and trilinear elements,
which the simplex mesh-convergence study measures (`docs/verification.md`).

### The quadratic tetrahedron

A part meshed from CAD is almost always a tetrahedral mesh, and on linear
tetrahedra it reads too stiff: the Tet4 above cannot bend. The ten-node
tetrahedron (`elements/Tet10.cpp`) adds a node on every edge. With the
barycentric coordinates `L_0 = 1 - xi - eta - zeta, L_1 = xi, L_2 = eta,
L_3 = zeta` of the reference cell,

```
  corners 0..3:     N_i  = L_i (2 L_i - 1)
  edge nodes 4..9:  N_ab = 4 L_a L_b     on the edges 0-1, 1-2, 2-0, 0-3, 1-3, 2-3
```

which is the VTK, Abaqus/CalculiX (C3D10) and scikit-fem node order; Gmsh
numbers the last two edge nodes the other way round and its reader swaps
them. The map `x(xi) = sum_a N_a x_a` is isoparametric, so an edge node may
sit off the straight edge on a curved CAD surface, as a Gmsh
`Mesh.ElementOrder = 2` mesh places it. On a straight-sided cell `grad N`
is linear, so the strain varies linearly and any quadratic displacement
field - pure bending among them - is reproduced exactly (the quadratic
patch test, `docs/verification.md`).

| Integral | Rule | Exact for |
|----------|------|-----------|
| stiffness `int B^T D B` and geometric stiffness | symmetric 4-point rule (degree 2), the rule of C3D10 | straight-sided cells |
| consistent mass `int rho N^T N` | 64-point collapsed Gauss (degree 5) | straight-sided cells |
| volume `int det J` | 27-point collapsed Gauss | any cell: `det J` is a cubic |
| face traction on a 6-node face | 16-point collapsed Gauss on `N^T t |x_,r x x_,s|` | flat faces: `S t / 3` on each edge node, nothing on the corners |

On a curved cell `B^T D B` is rational and the 4-point rule approximates it,
as every code using C3D10 does. The consistent mass matrix of the quadratic
tetrahedron has negative row sums at the corners, so a row-sum lumped mass
would be indefinite; lumping uses Hinton-Rock-Zienkiewicz scaling of the
diagonal instead, which on a straight cell gives each corner `m/36` and each
edge node `4m/27`. Element quality is the Tet4 measure of the corner
tetrahedron times the Jacobian ratio `min det J / max det J` over the ten
nodes (1 for a straight cell, falling as curved edges distort it); a cell
whose map folds at an integration point raises `MeshError`.

Two ways to a Tet10 mesh: read one (Abaqus/CalculiX `C3D10`, Gmsh element
type 11), or set `mesh.order: 2` on a tetrahedral deck, which elevates the
Tet4 mesh by putting a node at every edge midpoint (straight-sided cells;
the faceted geometry of the linear mesh stays). Boundary faces keep their
corner nodes as the key, so node sets and faces carry over.

## 3. Assembly

`fem/Assembler.cpp` builds a triplet list and compresses it to CSC:

```
  K = sum_e  s_e  A_e^T K_e^0 A_e            M = sum_e  m_e  A_e^T M_e^0 A_e
```

`A_e` is the (implicit) DOF gather operator, and `s_e`, `m_e` are per-element
scale factors: both are 1 for a plain analysis, and the topology optimizer
supplies `s_e = E(rho_e)/E_0` and a mass interpolation factor `m_e`.

On a **uniform structured mesh** all cells are geometrically identical, so
`K_e^0` and `M_e^0` are integrated once and reused; this is the dominant saving
in an optimisation loop. The generic per-element path is always available, and
the test suite asserts the two produce bit-comparable matrices.

Every assembly after the first reuses the **sparsity pattern**. The first
call records, per element, where each of its node blocks lands in the CSC
arrays of `K` (and of the reduced `K_ff`); later calls scatter the scaled
element matrices straight into those slots, in element order, which is the
order the triplet path sums duplicates in - so the result is bitwise
identical to `setFromTriplets`, which a test asserts, while skipping the sort
and the triplet storage. On the 830 115-DOF Hex8 block of the scaling
benchmark this took the assembly from 12.4 s to about 1 s and the peak memory
of the run from 4.8 GB to 3.9 GB. A zero scale factor on some element falls
back to the triplet path, whose pattern would differ.

## 4. Boundary conditions and the reduced system

Dirichlet conditions are applied by **partitioning**, not by penalty terms
(`fem/DofManager.cpp`, `fem/StaticAnalysis.cpp`). Splitting the DOFs into free
`f` and prescribed `p`,

```
  | K_ff  K_fp | | u_f |     | f_f     |
  |            | |     |  =  |         |
  | K_pf  K_pp | | u_p |     | f_p + r |
```

the first block row gives the system actually solved,

```
  K_ff u_f = f_f - K_fp u_p
```

and the reactions come from the *full* residual

```
  r = K u - f
```

whose entries at prescribed DOFs are the support reactions and whose entries at
free DOFs are the (tiny) solver residual. This gives exact reactions with no
second assembly, and keeps `K_ff` symmetric positive definite - no artificial
large diagonal entries to wreck the conditioning or the eigenvalue spectrum.

Because `K t = 0` for a rigid translation `t`, summing the residual over all
DOFs in one direction gives `sum(r) + sum(f) = 0` identically. That is the
global force balance every static result reports, and it is a genuine check on
the assembly and the solve rather than an identity of the post-processing. The
reported relative error divides the residual resultant by the gross size of the
applied nodal forces, `sum_n |f_n|` (and the moment residual by
`sum_n |x_n| |f_n|`), which is the scale of the round-off in the sums. Dividing
by the resultant instead would fail every self-equilibrated load - a thermal
strain, or a self-weight carried by a traction - whose resultant is itself
round-off. A case driven by prescribed displacements alone applies no force;
its reactions balance among themselves and are measured against their own
gross size, `sum_n |r_n|` (and `sum_n |x_n| |r_n|`).

## 4b. Pressure, volume and thermal loads

**Pressure** (`elements/FaceGeometry.cpp`). A boundary face with nodes `x_a` and
face shape functions `N_a(s, t)` carries

```
  f_a = - integral p N_a a(s, t) ds dt ,    a = x_s x x_t  (surface),  a = (y_s, -x_s) t  (edge)
```

with `a` the unnormalised area vector, so `|a| ds dt` is the area element and
the load follows the normal at every point of a curved face. On a Tet10 face
(six nodes, `x_s`, `x_t` linear) the integrand is of degree 4 and is integrated
with three collapsed Gauss points per direction, which is exact; on a flat
face it is of degree 2. The linear analyses load the undeformed face. The
derivative of the load with respect to the face's own nodal positions,
`d f_a / d x_b = -p N_a [N_b,s (-skew(x_t)) + N_b,t skew(x_s)]` integrated over a
surface (`face_pressure_stiffness`), is what a large-deflection analysis needs
for a pressure that follows the deforming face; it matches central differences
of the load to `1e-7` on distorted faces of all four shapes (tests).

**Body loads** (`fem/Loads.cpp`). A force density affine in position,
`b(x) = b_0 + B x`, is interpolated exactly by the shape functions of an
isoparametric element (`x = sum_b N_b x_b` exactly, curved cells included), so

```
  f_a = integral N_a b dV = sum_b ( integral N_a N_b dV ) b(x_b) = ( M_e(rho = 1) b_hat )_a ,
```

the unit-density consistent mass times the nodal values of `b`. Self-weight is
constant, a uniform force density is constant, and the centrifugal load of a
rotation `omega` about an axis `e` through `c` is the affine
`rho omega^2 (I - e e^T)(x - c)`: all three are integrated exactly to the mass
matrix's own quadrature, which is exact for these integrands on straight
cells. The resultant of self-weight is the model's mass times `g` to
round-off (tested).

**Thermal strain** (`material/IsotropicMaterial.cpp`, `fem/Loads.cpp`). With
`dT = N^T (T_e - T_ref)` interpolated from the nodal temperatures,

```
  plane stress   eps_0 = alpha dT {1, 1, 0}               (sigma_zz = 0, eps_zz free)
  plane strain   eps_0 = (1 + nu) alpha dT {1, 1, 0}      (eps_zz = 0;
                 sigma_zz = nu (sigma_xx + sigma_yy) - E alpha dT)
  3-D            eps_0 = alpha dT {1, 1, 1, 0, 0, 0}
  f_th = integral B^T D eps_0 t dV ,
```

integrated with the element's stiffness rule - exact for the linear elements,
the four-point rule of the Tet10 - so that any free expansion the element can
represent (a uniform `dT`; a linear one on affine simplices) is reproduced with
zero stress. The stress is `D (B u - eps_0)`, and the elastic strain energy is

```
  U = 1/2 u^T K u - u^T f_th + 1/2 integral eps_0^T D eps_0 t dV ,
```

not `1/2 u^T K u`: the thermal part of the load works against the free
expansion, not against the stiffness. A linear buckling analysis includes the
thermal prestress in `K_G`.

**Steady conduction** (`fem/HeatConduction.cpp`). The temperature of a
`conduction` load case solves `-div(k grad T) = Q` with prescribed
temperatures, surface fluxes `q` into the body, and convection
`h (T - T_inf)` leaving it:

```
  ( integral k grad N^T grad N t dV + integral_Gamma_h h N N^T t dGamma ) T
      = integral N Q t dV + integral_Gamma_q N q t dGamma + integral_Gamma_h h T_inf N t dGamma ,
```

on the structural mesh and shape functions (a plane model conducts in its
plane, insulated on its faces), prescribed temperatures by partitioning. The
residual `r = K T - F` at the prescribed nodes is the heat entering there; the
solve reports sources, fluxes and net convective inflow against it and
refuses a mismatch above `1e-6` of the largest heat flow of the problem.

## 5. Linear solvers

| Name | Implementation | Use |
|------|----------------|-----|
| `auto` (default) | `simplicial_ldlt` up to 50 000 free unknowns in 2-D and 10 000 in 3-D, `amg_cg` above | every deck that does not name a solver |
| `simplicial_ldlt` | Eigen `SimplicialLDLT` with AMD ordering | the reference; SPD systems; reports non-positive or tiny pivots as a singular system |
| `amg_cg` | CG preconditioned by SparLab's smoothed-aggregation multigrid (`fem/Multigrid.cpp`) | large systems, 3-D above all |
| `conjugate_gradient` | Eigen `ConjugateGradient`, Jacobi preconditioner | the baseline that shows what multigrid buys |
| `simplicial_llt` | Eigen `SimplicialLLT` with AMD | strictly SPD |
| `sparse_lu` | Eigen `SparseLU` with COLAMD | indefinite systems |
| `dense_lu` | Eigen `PartialPivLU` | verification of small models; refuses above 4000 DOFs |

After every solve the scaled residual `||A x - b|| / max(||b||, tiny)` is
compared with `solver.linear.residual_tolerance`, together with the backward
error `||A x - b|| / || |A||x| + |b| ||`, the residual over the norm of its
gross terms. The rounding of the sums that form `A x - b` is a few machine
epsilons times that denominator, so no solver can leave a smaller backward
error than that. A backward-stable solve, such as the Cholesky factorisation
of an SPD matrix, leaves about that much.

The two measures part company in a bending-dominated structure. In a slender
beam the entries of `K u` are sums of terms far larger than their result.
Take a Tet10 cantilever with `L/h = 100`, 50 cells long and 2 x 2 cells in
section. There `|| |K||u| + |f| ||` is `2.7e9` times `||f||`. Its Cholesky
solve leaves a scaled residual of `6.5e-7` and a backward error of
`2.4e-16`, about one `eps`.

A solve is accepted when either the scaled residual is within tolerance or
the backward error is at round-off, which is at most `64 eps`, about
`1.4e-14`. In the second case the residual is no larger than the rounding of
the sums that form it, which no floating-point solve can go below. Only a
backward-stable result reaches the test: the direct solvers check their
pivots, and the CG solvers throw unless they met their own tolerance
(`1e-12` by default). If both measures fail, the solve raises `SolverError`
with both measured values. A non-finite solution raises the same error. The
summary records both measures for each load case.

Neither measure bounds the *forward* error, which can be the condition
number times larger. On a Tet10 cantilever at `L/h = 100` (24 x 2 x 2
cells), scikit-fem's solution of the same system - itself backward stable,
`1.0e-16` against SparLab's `2.1e-16` - differs from SparLab's by `6.6e-8`
of the largest displacement; at `L/h = 50` by `8.0e-9`, the `(L/h)^3` of a
slender beam's conditioning. Nor can the backward error decide between two
codes on such a problem: scaling SparLab's solution by `1 + 1e-6` moves its
backward error only from `2.1e-16` to `3.3e-16`, because the residual that
error leaves, `1e-6 f`, is tiny beside the sums `|K||u|`. The
cross-validation therefore compares displacements, and its elastica deck
runs at `L/h = 50` (`docs/verification.md`, section 14). The
LDL^T path additionally inspects its pivots: a non-positive pivot, or a
`min/max` pivot ratio below `pivot_tolerance`, is reported as a singular system
together with the three modelling causes that usually produce it.

### Why an iterative solver

A sparse Cholesky factorisation of a 3-D stiffness matrix fills in. With a
nested-dissection-quality ordering its cost grows like `n^2` and its memory
like `n^(4/3)`, and the AMD ordering used here does somewhat worse: on the
Hex8 block of the scaling benchmark the factor holds about 560 entries per
unknown at 28 413 unknowns, against about 70 in `K` itself. Conjugate
gradients needs only products with `K`; the question is how many.
Preconditioned only by its diagonal, CG needs `O(sqrt(kappa))` iterations
and `kappa` grows like `h^-2`, so the count doubles with every halving of
the element size. A multigrid preconditioner keeps it nearly constant.

### Smoothed-aggregation multigrid

`fem/Multigrid.cpp` builds the hierarchy of Vanek, Mandel and Brezina (1996)
on `A_0 = K_ff`. On each level:

1. **strength of connection** between two nodes `I`, `J` compares the
   Frobenius norms of the matrix blocks,
   `||A_IJ||_F >= theta sqrt(||A_II||_F ||A_JJ||_F)` with `theta = 0.02`, so
   a weak link - solid next to SIMP void - does not glue two aggregates
   together;
2. **aggregation** groups each node with its strongly connected neighbours,
   greedily in node order in three passes;
3. **tentative prolongator.** The near-null space `B` is the rigid-body
   modes (two translations and a rotation in 2-D; three and three in 3-D),
   built from the node coordinates of the free unknowns. Restricted to each
   aggregate and orthonormalised by a rank-revealing QR, `B_a = Q_a R_a`,
   `Q_a` becomes the aggregate's block of `P_hat` and `R_a` its rows of the
   coarse near-null space, so `P_hat B_c = B` exactly - the coarse space
   represents every rigid motion of every aggregate;
4. **prolongator smoothing** `P = (I - omega D^-1 A) P_hat` with
   `omega = (4/3) / lambda_max(D^-1 A)`, the eigenvalue estimated by twelve
   Lanczos steps;
5. **Galerkin coarse operator** `A_{l+1} = P^T A_l P`, symmetrised, which
   keeps every level symmetric positive definite.

The recursion stops at 1500 unknowns, where a dense LDL^T factorisation
solves exactly. A coarsest pivot below `1e-13` of the largest means a rigid
body mode reached the coarse level - the model is under-constrained - and is
reported as such; the SIMP stiffness floor of `1e-9` stays orders of
magnitude above that.

One V-cycle, with a degree-3 Chebyshev smoother in `D^-1 A` on the interval
`[lambda_max / 30, 1.1 lambda_max]` before and after the coarse correction
(or a forward and a backward Gauss-Seidel sweep), is a symmetric
positive-definite operator, so the outer method is plain preconditioned CG,
stopped at `||b - A x|| <= tol ||b||` (`iterative_tolerance`, `1e-12` by
default). CG restarts from the true residual up to twice. When the
recursively updated residual meets the tolerance but the true residual
`b - A x` does not, and a restart no longer halves it, that is the accuracy
with which `b - A x` can be evaluated in floating point, and the solve is
accepted - the residual check above (`residual_tolerance`, `1e-8`) still
applies to it. A solve that exhausts its iteration budget raises
`ConvergenceError` naming the residual reached, the hierarchy and the keys
to change.

**In an optimisation loop** the matrix changes every iteration but its
pattern does not. The aggregates, the tentative prolongators and the
sparsity of every product are kept (`reuse_aggregates`), and only the
numerical part - smoothing, the Galerkin products, the smoother bounds and
the coarse factorisation - is redone; a test checks the reused hierarchy
equals a fresh one bit for bit. Each state and adjoint solve starts from the
previous iteration's solution of the same load case (`warm_start`).

**Determinism.** Every parallel kernel either writes rows independently
(matrix-vector and matrix-matrix products, the Chebyshev smoother) or runs
sequentially (aggregation, Gauss-Seidel), and every inner product sums
fixed-size chunks in a fixed order, so results are bitwise identical for any
number of OpenMP threads, which a test asserts.

The measured behaviour - iteration counts under refinement, time and memory
against the direct solver, agreement with it - is in `docs/verification.md`
and `docs/benchmarks.md`.

## 6. Stress recovery

```
  eps_e = B(xi_g) u_e                  sigma_e = s_e D eps_e
```

Element values are the average over the stiffness quadrature points; nodal
values are area-weighted averages of the adjacent element values (simple
averaging, used for smooth contours only - not a superconvergent patch
recovery). With `s_e` applied the result is the *macroscopic* stress carried by
the element, and the unscaled solid-material stress is reported alongside it.

Derived quantities:

```
  sigma_vm = sqrt( 1/2 [ (sxx-syy)^2 + (syy-szz)^2 + (szz-sxx)^2 ] + 3 (sxy^2 + syz^2 + szx^2) )
  szz = 0                            (plane stress;  syz = szx = 0 in both plane cases)
  szz = nu (sxx + syy)               (plane strain)

  sigma_{1,2} = (sxx+syy)/2  +-  sqrt( ((sxx-syy)/2)^2 + sxy^2 )          (plane)
  sigma_{1,2,3} = eigenvalues of the symmetric 3 x 3 stress tensor         (solid)

  U_e = 1/2 s_e u_e^T K_e^0 u_e      (element strain energy, exact)
```

On a solid mesh the three principal stresses come from a self-adjoint
eigen-decomposition of the stress tensor, and the test suite checks that
their sum, the sum of their pairwise products and their product reproduce the
three stress invariants.

## 7. Modal analysis

Undamped free vibration on the free DOFs (`fem/ModalAnalysis.cpp`):

```
  K_ff phi = lambda M_ff phi ,    lambda = omega^2 ,    f = omega / (2 pi)
```

Prescribed DOFs are removed by the same partitioning as the static solve, so no
constrained DOF can contribute a spurious mode.

**Algorithm.** Bathe subspace iteration with optional shift-invert. With
`q = min(n, max(2m, m+8))` trial vectors,

```
  Xbar_{k+1} = (K_ff - sigma M_ff)^{-1} M_ff X_k
  K_r = Xbar^T K_ff Xbar ,   M_r = Xbar^T M_ff Xbar          (q x q, dense)
  solve  K_r Q = mu M_r Q ,   X_{k+1} = Xbar Q
```

One sparse Cholesky factorisation serves every iteration and every vector. The
starting basis is deterministic: the mass diagonal, then unit vectors at the
DOFs with the largest `m_ii / k_ii` ratio, then a seeded PRNG if more are
needed. Because `Q^T M_r Q = I`, the rotated basis is automatically
M-orthonormal, and the modes are renormalised at the end so `phi^T M phi = 1`
holds to round-off.

For systems with 400 or fewer free DOFs a dense
`GeneralizedSelfAdjointEigenSolver` is used instead: it is faster at that size
*and* gives the test suite an independent reference for the iterative path. The
two agree to 1e-8 relative on a 40x8 cantilever.

**Stopping rule.** Both the relative eigenvalue change and every eigenpair
residual `||K phi - lambda M phi|| / ||lambda M phi||` must be below their
tolerances. Including the residual matters because the highest requested mode
converges last: its eigenvalue can look settled while its vector is not.

**Validity screening.** Each converged pair is checked for a non-finite
eigenvalue (error), a negative eigenvalue beyond round-off (error - `K_ff` is
not positive definite, which is impossible for a properly constrained
linear-elastic model), and an eigenvalue indistinguishable from zero at the
problem's stiffness/mass scale (warning - an unsuppressed rigid-body or
mechanism mode).

**Mass conservation.** Summing every entry of the assembled mass matrix gives
`sum(M) = dim * rho V`, one factor per translation direction (2 on a plane
mesh, 3 on a solid one). This identity is exact for both the consistent and
the lumped matrix and is used as a verification check (measured error
`<= 5e-14` in both dimensions).

### Analytical references

Two independent closed forms are used as validation references:

```
  bending, fixed-free beam:   f_n = (beta_n L)^2 / (2 pi L^2) sqrt(E h^2 / (12 rho))
                              beta_1 L = 1.87510407, beta_2 L = 4.69409113, ...

  axial, fixed-free rod:      f_n = (2n - 1) / (4 L) sqrt(E / rho)
```

The axial reference is *exact* for a bar in uniaxial stress, which plane stress
at `nu = 0` reproduces exactly; the measured agreement is 4e-6 relative. The
bending reference is not exact for a 2-D model, because Euler-Bernoulli theory
omits the shear deformation and rotary inertia the 2-D model includes - so the
computed frequencies must sit slightly *below* theory, and the gap must grow
with mode number. Both behaviours are asserted in the test suite, which is
stronger than asserting agreement.

## 7b. Linear buckling

`fem/Buckling.cpp`. A load case `f` gives the linear displacement `K u = f`
and with it the stress `sigma = D B u`. Linear (bifurcation) buckling asks
for which multiple `lambda` of that stress state the stiffness loses
definiteness:

```
  (K_ff + lambda K_G,ff(u)) phi = 0 ,
  K_G,e(u_e) = int_e t G^T S(sigma) G dOmega ,   S = blockdiag(sigma, ..., sigma)
```

with `G` the matrix of shape-function gradients (so that `phi^T K_G phi =
int sigma_ij phi_k,i phi_k,j`), integrated with the element's stiffness
rule. The structure is predicted to buckle at the load `lambda f`: `lambda
> 1` is a safety factor on the load case, `lambda < 1` means it buckles
before reaching it. A load that only stretches the structure has no
positive `lambda` and is reported so (the reversed load may buckle it at a
negative `lambda`). Assumptions: linear elasticity up to buckling, the
pre-buckling state is the linear solution scaled, loads keep their
direction, and the geometry is perfect - `lambda` is the bifurcation load
of the ideal structure, an upper bound on the collapse load of a real one.

**Algorithm.** `mu = 1/lambda` solves the symmetric-definite pencil
`(-K_G) phi = mu K phi`, and the smallest positive load factors are its
largest eigenvalues. Subspace iteration on `X <- K^-1 (-K_G) X` with a
Rayleigh-Ritz projection onto `q = min(n, max(2m, m + 8))` vectors reuses
the factorisation of `K_ff` the static solve made. `-K_G` is indefinite,
so load factors of the reversed load converge alongside; when more of them
than the spare slots outrank the wanted ones - a structure mostly in
tension, or SIMP void in tension - the iteration switches to the buckling
spectral transformation (Grimes, Lewis and Simon 1994)

```
  X <- (K + sigma K_G)^-1 K X ,     nu = lambda / (lambda - sigma)
```

which maps every negative `lambda` into `(0, 1)` and the wanted ones above
1. For `0 < sigma < lambda_1` the shifted matrix is positive definite, and
by Sylvester's law of inertia the negative pivots of its `LDL^T`
factorisation count the load factors below `sigma`; that count places
`sigma`. Convergence needs both the relative change of the requested load
factors and every residual `||K phi + lambda K_G phi|| / ||K phi||` below
their tolerances (`1e-8`, `1e-6` by default). Systems with at most 400 free
DOFs use a dense generalised eigensolve. Modes are normalised to `phi^T K
phi = 1` and each carries the fraction of its strain energy in elements at
density `>= 0.5` (1 for a plain analysis), so a mode localised in void
material is visible.

Where it runs: `sparlab_solve --buckling n` (or a `buckling` deck section)
checks the analysed model; in a topology run the same section checks the
full solid domain and the exported part - the density thresholded at 0.5,
largest face-connected group, full material - after the optimisation. The
buckling constraint on the SIMP design is in `docs/topology_optimization.md`
(section 5d).

**Plane models.** A plane-stress model buckles only in its plane: a strut
of a 2-D design can buckle sideways within the plane, but a thin plate
cannot buckle out of it, which needs shell or solid elements.

## 7c. Geometrically non-linear statics

`fem/TotalLagrangian.cpp`, `material/Hyperelastic.cpp`,
`fem/NonlinearStatic.cpp`. Large displacement and rotation, written on the
reference configuration (total Lagrangian). With `H = grad_X u` the
displacement gradient, `F = I + H` the deformation gradient and
`E = (H + H^T + H^T H) / 2` the Green-Lagrange strain, equilibrium is the
stationarity of

```
  Pi(u, lambda) = integral_Omega0 W(E) dV0 - lambda f_ext(u) . u   (conservative loads)
  R(u, lambda)  = f_int(u) - f_ext(u, lambda) = 0 ,
  f_int,e       = integral B_NL^T S dV0 ,       S = dW/dE  (second Piola-Kirchhoff)
  K_T,e         = integral ( B_NL^T D_T B_NL + G^T S G ) dV0 - d f_ext,e / d u ,
```

where `B_NL` maps nodal displacement variations to `delta E` (its rows are
`F_ki N_a,i` on the normal components, `F_ki N_a,j + F_kj N_a,i` on the
engineering shears), `D_T = dS/dE` is the material tangent, and the second
term is the initial-stress stiffness of section 7b with the current `S`. The
element's stiffness quadrature integrates everything; a plane model
multiplies by the thickness and a solid by 1.

**Materials.**

- *Saint Venant-Kirchhoff*: `S = D E`, `W = 1/2 E : D : E`, with the linear
  elasticity matrix `D` of the stress state - large rotation, small strain.
  It is the simplest frame-invariant extension of the linear law and reduces
  to it as the strains vanish, but under compression its force falls again
  below a stretch of `1/sqrt(3)`, so it is not a large-strain law. In plane
  stress `S_33 = 0` fixes the thickness strain
  `E_33 = [(1 + nu) E_theta - nu (E_11 + E_22)] / (1 - nu)` (with `E_theta = 0`
  without temperature), which enters `J` and the Cauchy stress.
- *Compressible neo-Hookean*: `W = mu/2 (tr C - 3) - mu ln J + lambda/2 (ln J)^2`,
  `S = mu (I - C^-1) + lambda ln J C^-1`,
  `D_T,ijkl = lambda C^-1_ij C^-1_kl + (mu - lambda ln J)(C^-1_ik C^-1_jl + C^-1_il C^-1_jk)`,
  with `C = F^T F`, `J = det F` and the Lame constants of `E` and `nu`, so it
  too reduces to linear elasticity at small strain. `J - 1` is formed from the
  invariants of `H` and `ln J` by `log1p`, without the cancellation of
  `det F - 1` at small strain. Plane strain and solids only; an inverted
  point (`J <= 0`) is an error that halves the load step.

**Thermal strain at finite strain.** A temperature change `dT` stretches a
free element by `theta = 1 + alpha dT` in every direction. The thermal
stretch splits off multiplicatively, `F = F_e theta I` (Lu and Pister 1975),
so `E = theta^2 E_e + E_theta I` with `E_theta = alpha dT (1 + alpha dT / 2)`, the
Green-Lagrange strain of the free stretch. The elastic energy per unit volume
of the expanded, stress-free body is `1/2 E_e : D : E_e`, which per unit
reference volume is `theta^3` times that:

```
  W = (E - E_theta I) : D : (E - E_theta I) / (2 theta) ,    S = D (E - E_theta I) / theta ,
  dS/dlambda = -D [ alpha dT0 I + alpha dT0 (E - E_theta I) / theta^2 ]   (dT = lambda dT0),
```

in Voigt form, `E_theta I` and `alpha dT0 I` being the linear model's
thermal-strain vectors of the stress state at those changes (the
plane-strain `(1 + nu)` factor included). A freely heated body takes the
stretch `1 + alpha dT` exactly, stress-free, and a restrained one carries
the stress of the expanded material. The two obvious alternatives fail
here. Subtracting the linear thermal strain, `S = D (E - alpha dT I)`, gives
the free body the stretch `sqrt(1 + 2 alpha dT)`. Subtracting the Green
strain of the free stretch without the `1/theta`, `S = D (E - E_theta I)`,
gets the free stretch right, but measures the elastic energy per reference
instead of per expanded volume and gives a restrained body a stress `theta`
times too high. For a cube held between two walls and free sideways at
`alpha dT = 0.05` (`E = 1 GPa`, `nu = 0.3`), the closed forms give
`sigma_xx = -43.07 MPa` with the split and `-45.22 MPa` (`theta = 1.05` times
it) without the `1/theta`. The unit tests reproduce the split's free state
on Q4, Hex8 and Tet10 cells and its restrained state on Hex8 and Tet10
cells, to `1e-10` or better. The thermal strain is written for Saint
Venant-Kirchhoff only.

**Loads**, all scaled by the one load factor `lambda`:

- point loads, tractions, self-weight and body forces are *dead*: fixed
  vectors of the reference configuration (a body force per unit reference
  volume, which conserves mass);
- a *follower pressure* acts on the deformed face,
  `f_a = -lambda p integral N_a (x_s x x_t) ds dt` at the current nodal
  positions, with the load stiffness of section 4b. Over a free edge that
  stiffness is not symmetric, and the assembled tangent is tested for
  symmetry (`max |A - A^T| <= 1e-12 max |A|`): a symmetric one is factorised
  by `LDL^T`, a non-symmetric one by `LU`. `follower_pressure: false`
  keeps the pressure on the reference faces;
- a *rotation* loads each point at its deformed position,
  `rho omega^2 P (x - c)` with `P = I - e e^T`, so
  `f = lambda omega^2 (M (x) P) (x_hat - c_hat)` with `M` the consistent mass
  matrix (whose entries carry `rho`) and the load stiffness
  `-lambda omega^2 M (x) P`, the spin softening;
- the temperature change is `lambda dT0`, and prescribed displacements are
  `lambda u_0`.

**Newton's method under load control.** From a converged state at
`lambda`, the step to `lambda + dl` starts from the tangent predictor
`K_T du = dl q - R`, `q = -dR/dlambda` (with the increment of the prescribed
displacements), then corrects with `K_T du = -R`. Each direction is scaled
by an energy line search: `g(a) = du . R(u + a du)` is the derivative of the
potential along `du`; the full step stays when `|g(1)| <= 0.8 |g(0)|`
(Crisfield's slack), else regula falsi on `g` finds a shorter one, no
shorter than 0.1. The norm of the residual would be the wrong measure: in a
slender structure a good Newton step stirs up large axial forces, and
backtracking on `||R||` stalls Newton into a linear crawl (on the elastica it
took 158 iterations and a halving where the energy search takes 87). A direction of ascent, as
a non-symmetric or indefinite tangent can give, takes the full step. A step
converges when

```
  ( ||R_f|| <= tol_r s  and  ||a du|| <= max(tol_u ||u - u_0||, 64 eps ||u||) )
  or ||R_f|| <= max(1024 eps g_f, 64 eps || |K_T| |u| ||_f) ,
```

`s` the largest of the applied loads, the reactions and the thermal forces,
`g_f` the gross assembly `sum_e |f_e| + |f_ext|`, `u - u_0` the step's
displacement increment; the residual of the accepted state is checked once
more. The second line is the round-off floor, below which no tolerance can
be met: the residual stagnates at about 150 `eps g_f` on the plane
cantilever, and at 0.3 `eps || |K_T| |u| ||` on the elastica, where each
stored displacement is exact only to its own rounding. A floor counts only
while it is below `1e-6 s`: past a limit or a plastic collapse load the
displacement runs away and `|K_T| |u|` with it, until the "floor" exceeds
the load itself - measured: a perfectly plastic tube driven 5 % past its
collapse pressure "converged" at strains of 1e11 before this bound, and now
stops with the collapse bracketed. A step that fails is halved (at most
`max_cuts` times in a row); three or fewer iterations lengthen the next one
by half, never beyond the first.

**Limit and bifurcation points.** Past a limit point load control has no
nearby solution: Newton fails, or - worse - converges on a distant branch,
a snap-through that is not quasi-static. Two tests reject such a step even
when Newton converges:

- from a stable state (no negative pivot), an iteration that meets a
  tangent with a negative pivot. Near a stable stretch of the path, the
  iterates of a short step stay near it;
- a converged state farther from the tangent predictor than the predicted
  increment. On a smooth path the corrector shrinks with the step - the
  ratio measured 0.01 to 0.78 on the verification problems - while a jump to
  another branch leaves it at order 1 or larger (21 on the snapping arch),
  however short the step.

The lowest load factor rejected this way is a suspected critical point. The
steps close in on it by halving; a step that converges beyond it clears it,
and after `max_cuts` failed halvings the run stops and reports the bracket
`[lambda, lambda_c]` - a limit point if a step jumped, a loss of
definiteness (limit or bifurcation) otherwise. At a bifurcation of a
perfect structure the fundamental path continues, but unstable, which is
why a negative pivot stops load control there rather than letting it go on.

**Arc length.** Crisfield's cylindrical method makes the load factor an
unknown and fixes the step length instead: `||du||^2 = ds^2`, with `du` the
free displacement increment of the step. The predictor follows the tangent
solution `K_T u_q = q` in the direction of the last converged increment.
Each iteration solves `K_T u_r = -R` and `K_T u_q = q` and sets
`du <- du + u_r + d(lambda) u_q`, taking the root of the quadratic in
`d(lambda)` whose `du` turns least from the previous iterate. The first arc
length is that of the first of `steps` load increments; each next one
scales by `sqrt(desired_iterations / iterations)`, clamped to `[0.5, 2]`
per step and to `[min_arc_ratio, max_arc_ratio]` times the first. A failed
step halves the arc length. A step that would pass the target load factor
is replaced by a load-controlled step from the last state onto it.

**Stability.** For a symmetric tangent the `LDL^T` factorisation reports
its inertia: the number of negative pivots equals the number of negative
eigenvalues (Sylvester), so a converged state is stable when there are
none. Holding one displacement instead of its force - displacement control
- removes one row and column from the free system. Haynsworth's inertia
additivity then says the force-controlled tangent has exactly one more
negative pivot than the displacement-controlled one wherever the force
falls as that displacement grows, `dF/d delta < 0`, and the same number
elsewhere. The arch study checks this at every step. A non-symmetric
tangent (a follower pressure over free edges) is factorised by `LU`, which
reports no inertia; the summary then says that stability was not assessed.

**Results.** Each converged step records the load factor, iterations,
halvings, residual, arc length, negative pivots, largest displacement and
the monitors - a mean displacement over a node region, or the sum of its
reactions. The final state adds the Cauchy stress `sigma = F S F^T / J`
(element averages of the integration points, with `sigma_zz` in plane
strain and the thickness stretch in `F` in plane stress), its von Mises
stress, `S` itself, the largest Green-Lagrange strain, the smallest `J`,
the strain energy, and the balance of the applied loads against the
reactions in the deformed configuration: forces, and moments about the
deformed positions.

## 7d. Plasticity

`material.plasticity` makes a material elastoplastic in the non-linear
analysis: J2 (von Mises) plasticity with isotropic and kinematic hardening,
integrated at every integration point by the backward-Euler radial return of
Simo and Hughes (*Computational Inelasticity*, 1998, boxes 3.1 and 3.2),
with its consistent tangent (`src/material/Plasticity.cpp`,
`src/fem/Elastoplastic.cpp`).

**The law.** The strain splits additively,
`eps = eps_e + eps_p + alpha dT I`, the stress is `sigma = K tr(eps_e) I +
2 G dev(eps_e)`, and with the relative stress `xi = dev(sigma) - beta` the
yield function is

```
  f = ||xi|| - sqrt(2/3) sigma_y(a) <= 0 ,
  sigma_y(a) = sigma_y0 + H a + Q (1 - exp(-delta a)) ,
```

`||.||` the tensor norm. The flow is associative, `d eps_p = d gamma n` with
`n = xi / ||xi||`; the accumulated plastic strain grows by
`da = sqrt(2/3) d gamma` (in uniaxial tension it is the plastic strain) and
the back stress by `d beta = (2/3) H_kin d gamma n` (Prager). Every point
carries the full 3-D state; Voigt strains have engineering shears, stresses
tensorial components.

**The return.** From the internal variables of the last converged step, the
trial deviatoric stress `s_tr = 2 G dev(eps - eps_p,n - alpha dT I)` and
`xi_tr = s_tr - beta_n` either stay inside the surface - an elastic step -
or return onto it along `n = xi_tr / ||xi_tr||`, which the return does not
change (radial return). The multiplier solves

```
  g(dg) = ||xi_tr|| - (2 G + 2/3 H_kin) dg - sqrt(2/3) sigma_y(a_n + sqrt(2/3) dg) = 0 ,
```

in closed form for linear hardening and by Newton for Voce saturation (`g`
is convex and decreasing there, so Newton from `dg = 0` rises monotonically
onto the root), to `1e-13` of the radius `sqrt(2/3) sigma_y` - or, far
outside the surface, where one ulp of `||xi_tr||` exceeds that (beyond
about 450 radii: a deviatoric trial increment of order one, which the
Newton iterates of a large-strain run reach), to the round-off of
`||xi_tr||` itself (`8 eps ||xi_tr||`). Then `sigma = sigma_tr - 2 G dg n`,
`eps_p = eps_p,n + dg n`, `beta = beta_n + (2/3) H_kin dg n`. On a path along
which `n` keeps its direction - uniaxial stress, for one - the return is
exact for linear hardening and any step: the unit tests and the cycle study
compare it with the closed-form uniaxial curves to round-off (`1e-12`,
`2e-14`). The consistent tangent, the derivative of the returned stress,

```
  C = K 1 x 1 + 2 G theta (I - 1 x 1 / 3) - 2 G theta_bar n x n ,
  theta = 1 - 2 G dg / ||xi_tr|| ,
  theta_bar = 1 / (1 + (sigma_y'(a) + H_kin) / (3 G)) - (1 - theta) ,
```

keeps Newton's quadratic convergence; it is symmetric (associative flow) and
checked against central differences of the return (1e-6, in 3-D, plane strain
and plane stress). A point that yielded in its last converged step and still
sits on its surface returns elastically at zero increment - but with the
continuum elastoplastic tangent (`dg = 0` in `C`), so that the next step's
tangent predictor already sees the softer response of continued loading.

**Plane states.** Plane strain returns with `eps_33 = 0` (or the
mean-dilatation value below) and reports `sigma_zz`. Plane stress finds
`eps_33` at every point by Newton on `sigma_33(eps_33) = 0`, from the
elastic predictor (exact for an elastic step) with the consistent `C_33,33`,
to `1e-12` of the stress; the tangent is condensed,
`C_ab - C_a3 C_3b / C_33`. A thermal strain enters every normal component;
at fixed displacement its load rate is `-C alpha dT_0 m` with the consistent
(condensed) tangent.

**Hill48 anisotropy and Chaboche kinematic hardening.** With
`yield_criterion: hill48` or `backstresses` the return is a closest-point
projection (Simo and Hughes, box 3.4) of the general law
(`src/material/Plasticity.cpp`, `return_general_3d`; the header
`Plasticity.hpp` carries the full derivation):

```
  f = sqrt(xi^T P xi) - sigma_y(a) <= 0 ,   xi = sigma - sum_i alpha_i ,
  d eps_p = d lambda m ,   m = P xi / sigma_bar ,   da = d lambda ,
  d alpha_i = (2/3) C_i d eps_p - gamma_i alpha_i d lambda ,
```

P the Hill48 matrix in global axes, `T^T P_mat T` with `P_mat` built from
`F, G, H, L, M, N` in the frame RD, TD, ND and `T` the rotation of tensorial
Voigt stresses, found once when the material is set (von Mises is
`F = G = H = 1/2, L = M = N = 3/2`). `m` is the engineering form of the flow
direction (the derivative by Voigt components doubles the shears); its
normal part sums to zero, so the flow is isochoric and mean dilatation stays
valid. `lambda` is the work-conjugate accumulated strain,
`xi : d eps_p = sigma_bar d lambda`, which for von Mises is the radial
return's `sqrt(2/3) dg`. Over a step with the direction of its end each
backstress is `alpha_i = theta_i alpha_i,n + (2/3) C_i chi_i m_t` - the
exponential `theta = exp(-gamma dl)`, `chi = (1 - theta) / gamma`, exact for a
fixed direction, or backward Euler `theta = 1 / (1 + gamma dl)`,
`chi = theta dl` - and, the elasticity being isotropic and `m` deviatoric,
`xi = xi_tilde(dl) - b(dl) m_t` with `xi_tilde = s_tr - sum theta_i alpha_i,n`
and `b = 2 G dl + (2/3) sum C_i chi_i`: the backstresses drop out and Newton
solves seven equations for `(xi, dl)`,

```
  R1 = xi - xi_tilde + b S^-1 P xi / sigma_bar = 0 ,   R2 = sigma_bar - sigma_y(a_n + dl) = 0 ,
```

from the trial point with a backtracking line search, to `1e-13 sigma_y` -
tight enough that stress recovery's return from the committed state is
elastic. Differentiating the converged equations by the strain, which enters
only through `s_tr = 2 G dev(eps)` (the linear map `D_dev`), gives
`[A; a^T] = J^-1 [D_dev; 0]` with J the Newton Jacobian, and the consistent
tangent

```
  C = C_e - 2 G [ m_t a^T + (dl / sigma_bar) S^-1 (P - m m^T) A ] ,
```

symmetric for linear kinematic hardening and non-symmetric as soon as a
backstress recovers; the element then keeps it unsymmetrised and the solver
factorises it by LU. At `dl = 0` it is the continuum tangent the loading flag
uses. Von Mises without backstresses keeps the radial return above, bit for
bit; the general return with `r = 1` reproduces it (stress, tangent, state) to
`1e-12`. Hill48 with isotropic hardening follows the exact directional
hardening curve and r-value in uniaxial stress along any direction of the
sheet, and von Mises with backstresses the exact Armstrong-Frederick branch
solutions - monotonic, reversed and cyclic - at any step size, since both
paths keep their flow direction; backward Euler converges to them at first
order (`docs/verification.md`, section 24).

**The element.** With small strain the point strain is `eps = B u_e`, the
internal force `int B^T sigma dV` and the tangent `int B^T C B dV`. Plastic
flow is isochoric, and on an element with several integration points each
point's dilatation becomes a constraint: in plane strain and 3-D the fully
integrated Q4 and Hex8 then lock - measured on the thick tube below, the
collapse pressure comes out 1.2 % high at 16 cells through the wall and the
collapse plateau keeps rising. The mean-dilatation operator of Hughes (1980)
replaces every point's dilatation by the element's volume average,
`B_bar = B + (1/3) m (b_bar - b)^T` with the dilatation row `b = B^T m` and
`m = {1, 1, 1, 0, 0, 0}`; in plane strain it gives `eps_33 = (theta_bar -
theta) / 3` at a point, zero on average. It is the default for Q4 and Hex8
(`mean_dilatation: auto`) and brings them to the exact collapse load at
second order. The four-point Tet10 needs none: it converges at third order
without it, and with it at second order from below, its constant element
pressure oscillating. Plane stress has no incompressibility constraint, and
one-point elements have nothing to average.

**Finite kinematics.** With `kinematics: finite` the return takes the
Green-Lagrange strain `E` for the strain and gives the second
Piola-Kirchhoff stress `S`: `E = E_e + E_p + E_theta` with `S = D E_e`,
the elastoplastic counterpart of the Saint Venant-Kirchhoff law, exact under
any rigid rotation and meant for strains that stay small (at large
displacement and rotation with small strain the materially non-linear
relations may be written between `S` and `E`; Bathe, *Finite Element
Procedures*, 1996, ch. 6). The element takes `B_NL` of section 7c,
`f = int B_NL^T S dV_0`, `K = int B_NL^T C B_NL dV_0 + int (G^T S G) x I dV_0`.
The thermal strain is the Green strain of the free thermal stretch, so free
heating stays stress-free, but the elastic stiffness is not rescaled by the
thermal stretch as the Saint Venant-Kirchhoff law's is (a relative difference
of `alpha dT` in a thermal stress). Mean dilatation acts on the Green strain,
`E_bar = E + (1/3) m (mean(tr E) - tr E)`; its second variation adds
`(1/3) tr(S) (mean(G^T G) - G^T G) x I` to the geometric stiffness, and the
tangent stays symmetric and consistent (checked by central differences with
a rigid turn of 0.5 rad). Objectivity - a rigidly rotated state returns the
same `S`, plastic state and energy, and its forces turn with it - exact
homogeneous finite deformations on distorted meshes, heated or not, and free
heating are unit-tested, and an independent implementation of the same model
in scikit-fem reproduces four large-deflection decks - a Tet10 cantilever,
and a clamped beam driven into membrane action and back with combined
hardening on Hex8 and Q4 (E-bar) and in plane stress
(`docs/verification.md`, section 24). The model is not finite-strain
plasticity: beyond a Green strain of 0.05 the run warns - large strains
take the logarithmic kinematics below. A large rotation
superposed on plastic flow also sharpens a real effect: under a large
hydrostatic pressure (a 3 % volume change, 5 GPa) the
geometric stiffness outweighs a small plastic tangent and a homogeneous
state loses stability - at strains well beyond the model's range.

**Logarithmic-strain finite plasticity.** With `kinematics:
finite_logarithmic` the strains may be large (Miehe, Apel and Lambrecht,
*Anisotropic additive plasticity in the logarithmic strain space*, CMAME
191, 2002; `src/material/LogarithmicStrain.cpp`). Total Lagrangian as
above, but the return - the same J2, Hill48, Chaboche, plane-stress and
thermal return, unchanged - takes the logarithmic (Hencky) strain

```
  E_log = ln(C) / 2 = sum_a (ln lambda_a / 2) N_a N_a ,   C = F^T F = I + 2 E ,
  tr E_log = ln J ,
```

from the spectral decomposition of E (`lambda_a = 1 + 2 mu_a`, and
`ln(lambda_a) / 2 = log1p(2 mu_a) / 2`, which keeps the relative accuracy of
E at a strain of `1e-6`), and gives the stress T conjugate to it with its
algorithmic tangent `C_alg = dT/dE_log`. The additive split
`E_log = E_log,e + E_log,p` with `T = D E_log,e` is exact for a coaxial
homogeneous stretch - the log strains of successive coaxial stretches add -
so uniaxial J2 at any strain is the 1-D small-strain law in the log strain,
with T the Kirchhoff stress. The flow stays isochoric (`det F_p = 1`), the
Hill48 axes are those of the reference configuration and turn with the
material (no plastic spin of the texture), and the model is objective: C,
and so everything but the push-forward, ignores a superposed rotation.

With `dE_log = P dE` (`P = 2 dE_log/dC`) and `dE = B_NL du`:

```
  S = P^T T ,   f = int B_NL^T S dV_0 ,
  K = int B_NL^T (P^T C_alg P + T : L) B_NL dV_0 + int (G^T S G) x I dV_0 ,
```

`T : L` the Hessian of `T : E_log(E)` at fixed T (`L = 4 d2E_log/dC dC`,
the sixth-order term). By the Daleckii-Krein formulas, in the principal axes
of C,

```
  [dE_log]_ab = l[lambda_a, lambda_b] [dE]_ab ,
  T : d2E_log[dE1, dE2] = 2 sum_abc T_ab l[lambda_a, lambda_c, lambda_b]
                            ([dE1]_ac [dE2]_cb + [dE2]_ac [dE1]_cb) ,
```

with the divided differences of ln, `l[x, y] = (ln x - ln y)/(x - y)` and
`l[x, y, z] = (l[x, y] - l[y, z])/(x - z)`, `l[x, x] = 1/x`,
`l[x, x, x] = -1/(2 x^2)`. They carry the whole treatment of equal and nearly
equal eigenvalues: `l[x, y]` is evaluated as `log1p((x - y)/y)/(x - y)`,
exact to round-off at any gap and at none; `l[x, y, z]` as the difference
quotient when the spread `x - z` exceeds `5e-3 y`, and otherwise as the
Taylor series about the middle eigenvalue to fifth order,
`y^-2 sum_n (-1)^(n+1) h_n(u, w)/(n + 2)` with `u = (x - y)/y`,
`w = (z - y)/y` and `h_n` the complete homogeneous symmetric polynomials -
the threshold balancing the cancellation of the quotient (about
`4 eps / 5e-3 = 2e-13`) against the truncation of the series
(`(5e-3)^6 = 2e-14`), so that distinct, two equal, three equal and nearly
equal eigenvalues are all evaluated to about `1e-13` and continuously.
Checked against an independent matrix logarithm (Schur-Pade, `4.8e-15`) and
against fourth-order central differences of `E_log(E)` at distinct
eigenvalues, two equal, `C = I`, `C = s^2 I`, and gaps of `1e-9`, `1e-7`,
`1e-4` (two and three eigenvalues) and on both sides of the series threshold
(`P` to `4.1e-12`, `T : L` to `4.3e-12`: the differences' own accuracy).

The elasticity is Hencky's, quadratic in `E_log`, for elastic elements too
(the `material_model` is then not used; `neo_hookean` is refused). The
thermal strain is the log strain of the free thermal stretch
`1 + alpha dT` - the stretch of `finite` kinematics - `ln(1 + alpha dT) I`,
which the return takes as the temperature change `ln(1 + alpha dT)/alpha`:
a freely heated body expands stress-free by `1 + alpha dT`. In plane stress
the return finds `E_log,33` with `T_33 = 0`, which is `S_33 = 0` (the
normal is a principal direction of C), and the thickness stretches by
`exp(E_log,33)`.

`E_log` sees only `C = F^T F`, which is positive definite at an inverted
point (`det F < 0`) too, where `E_log`, the energy and `tr E_log =
ln |det F|` are those of a valid state with one axis reflected - an
equilibrium the solver would otherwise accept. Each point therefore checks
the sign of `det F` (the in-plane determinant in a plane model, whose `F_33`
is 1 or `exp(E_log,33)`), and an inverted one fails the step as with the
neo-Hookean law: the line search shortens a step into it and load control
cuts one that ends there, so a block crushed to zero volume stops short of
it instead of converging on the far side.

*Mean dilatation* acts on the logarithmic volumetric strain `ln J`: each
point's `ln J` is replaced by the element's volume average,
`E_log_bar = E_log + (1/3) m (mean(ln J) - ln J)`, exactly the log strain of
the modified deformation `F_bar = (J_bar/J)^(1/3) F`, so it stays objective.
With `b = m^T P B_NL` the variation of `ln J` and `B_bar = P B_NL +
(1/3) m (mean(b) - b)`, the force is `int B_bar^T T` and the consistent
tangent

```
  K = int [ B_bar^T C_alg B_bar + B_NL^T (T' : L) B_NL + (G^T S' G) x I ] dV_0 ,
  T' = dev T + p_bar m ,   S' = P^T T' ,
```

`p_bar` the element's mean Kirchhoff pressure: the second variation of
`int T : E_log_bar` is that of `int (T - p m) : E_log` plus `int p` times the
second variation of `mean(ln J)`, and both are initial-stress terms because
`ln J = m : E_log`. On the thick tube driven past its collapse at finite
strain (`--study logarithmic-tube`), Q4 and Hex8 with this averaging
converge at second order to the collapse pressure of the locking-free Tet10
(`4.7e-4 p_L` on 16 cells through the wall), while the fully integrated Q4
locks (`0.117 p_L` too high on 4 cells).

*Stress recovery* reports the Kirchhoff stress `tau = F S F^T` (with the
modified `F_bar`, `F_bar S_bar F_bar^T`, the same), the Cauchy stress
`tau / J_bar` with `J_bar = exp(tr E_log_bar)`, the second Piola-Kirchhoff
stress `S = P^T T`, and the log strain the return took (with the thickness
strain in plane stress).

Uniaxial tension to a stretch of 2 with linear and Voce hardening, a
tension-compression cycle to `+-0.3` log strain with Armstrong-Frederick
backstresses, and Hill48 pulled along RD and TD to a stretch of 1.6 (whose
plastic lateral log strains keep the ratio `r` exactly) are exact to
`1.2e-13` on a distorted Hex8 patch (`--study logarithmic-uniaxial`); the
element tangent is the central difference of the internal force (`1.1e-9`)
for Hex8, Q4 in plane strain and plane stress and Tet10, elastic and
plastic, J2, Hill48 and Chaboche, with and without mean dilatation at a
stretch of 1.35 and a turn of 0.6 rad; at strains of `1e-6` the element
agrees with the Green-Lagrange one to the order of the strain; every
non-linear verification deck runs with it (`docs/verification.md`,
section 24). On the plastic clamped beam deck it costs 1.45 times the
Green-Lagrange run (3.07 s against 2.12 s single-threaded, the same 81 Newton
iterations): an averaged Hex8 with its tangent takes 35 us against 23 us, for
the 3 x 3 eigen-decomposition (0.5 us a point) and the `T : L` term.

**History.** Each point keeps `(eps_p, beta, a)` of the last converged step.
Every Newton iterate, line-search probe and rejected step returns from
those, never from the last iterate, and only a converged step commits the
new ones - so halving a step leaves no trace, and the path is the sequence
of converged steps. A `load_path` visits turning points in order (unloading
leaves the permanent set and residual stresses). Load control's jump test
(section 7c) is off with plasticity: the elastic predictor of a step in which
points start to yield underestimates the increment by the ratio of elastic
to plastic stiffness. Without hardening a structure has a collapse load
beyond which no equilibrium exists; the largest converged load factor is a
lower bound on it, the arc-length method runs onto its plateau, and load
control stops with it bracketed.

**Small strain.** `kinematics: small_strain` keeps the linear strain on the
undeformed geometry: no geometric stiffness, pressures on the undeformed
faces, the rotation's load at the undeformed positions. With elastic
materials it reproduces the linear analysis (unit-tested to 1e-10, one
Newton iteration per step). The run reports the largest strain, the largest
infinitesimal rotation and the largest component of the quadratic strain
`H^T H / 2` it neglects, and warns when that exceeds a tenth of the largest
strain - where membrane action (a structure restrained against the motion a
rotation causes) makes it matter.

## 7e. Dynamics

The `transient` block integrates a load case in time and the
`frequency_response` block computes its steady harmonic response
(`src/fem/Dynamics.cpp`).

**Equations of motion.** On the free DOFs,

```
  M a + C v + K u = A(t) f ,     u_p = A(t) g ,
```

with the load vector `f` of the case (forces, pressures, body and thermal
loads), its prescribed displacements `g` and one amplitude `A(t)` for both:
a step, a piecewise-linear table (constant beyond its ends) or
`scale sin(2 pi f t + phase)`. `M` is the consistent mass (each element's
`int rho N^T N dV`, integrated exactly: 3 x 3 (x 3) Gauss points for Q4 and
Hex8, a 64-point collapsed rule for the Tet10) or the lumped one - the row
sums for the linear elements, the diagonal scaled to the element's mass for
the Tet10, whose corner row sums are negative (Hinton, Rock and Zienkiewicz,
1976). Damping is Rayleigh's, `C = a M + b K`, which leaves the modes
uncoupled with the damping ratio `zeta_j = a / (2 omega_j) + b omega_j / 2`.

**HHT-alpha.** The step `dt` is constant. Hilber, Hughes and Taylor (1977)
keep the Newmark relations

```
  u1 = u0 + dt v0 + dt^2 [(1/2 - beta) a0 + beta a1] ,
  v1 = v0 + dt [(1 - gamma) a0 + gamma a1]
```

and weight the equilibrium of the step's two ends,

```
  M a1 + (1 + alpha)(C v1 + K u1) - alpha (C v0 + K u0) = (1 + alpha) f1 - alpha f0 ,
  beta = (1 - alpha)^2 / 4 ,   gamma = 1/2 - alpha ,   -1/3 <= alpha <= 0 .
```

The method is unconditionally stable and second-order accurate for every
`alpha` in that range; its amplification of a mode whose period is short
against the step tends to `(1 + alpha) / (1 - alpha)`, so `alpha < 0` damps
the poorly resolved high frequencies while barely touching the resolved low
ones. CalculiX's `*DYNAMIC` uses the same sign convention, and its default is
`alpha = -0.05` (measured: a deck without `ALPHA` reproduces `ALPHA=-0.05`
digit for digit and differs from `ALPHA=0` by 1e-4).
`alpha = 0` is the trapezoidal rule (average acceleration): no numerical
dissipation, and a period lengthened by `(omega dt)^2 / 12` of itself - the
discrete solution of an undamped mode advances its phase by `theta` per step
with `tan(theta / 2) = omega dt / 2`, which the unit tests confirm to
round-off. A harmonic amplitude resolved by fewer than 20 steps per period
draws a warning.

**The displacement form.** SparLab solves for `u1`. With
`c0 = 1 / (beta dt^2)`, `c1 = 1 / (beta dt)`, `c2 = 1 / (2 beta) - 1`,
`c3 = gamma / (beta dt)`, `c4 = 1 - gamma / beta`,
`c5 = dt (1 - gamma / (2 beta))`,

```
  a1 = c0 (u1 - u0) - c1 v0 - c2 a0 ,   v1 = c3 (u1 - u0) + c4 v0 + c5 a0 ,
  [c0 M + (1 + alpha)(c3 C + K)] u1 = (1 + alpha) f1 - alpha f0
      + M (c0 u0 + c1 v0 + c2 a0) + (1 + alpha) C (c3 u0 - c4 v0 - c5 a0)
      + alpha (C v0 + K u0) ,
```

on the free DOFs, the prescribed ones moved to the right-hand side. The
effective stiffness is factorised once for the whole run. A prescribed DOF
follows the same kinematics: its displacement is `A(t_n) g` at every step,
its velocity and acceleration come from the Newmark relations of that
history (starting from `A'(0) g` and `A''(0) g`), so a shaken support moves
exactly as the method moves every other DOF.

**Initial state.** At rest (`u = 0`, `v = 0` on the free DOFs) or, with
`start: static`, in the static equilibrium of `A(0) f` (a preloaded structure
released); the initial acceleration solves
`M_ff a0 = A(0) f_f - (M a_p + C v + K u)_f`, so the run starts in
equilibrium - a load that acts at `t = 0` accelerates the structure at once.

**Energy balance.** With the kinetic energy `T = v^T M v / 2`, the strain
energy `U = u^T K u / 2`, the energy dissipated by damping
`D = sum dt/4 (v_n + v_n+1)^T C (v_n + v_n+1)` and the work of the loads and
of the reactions of the prescribed motion
`W = sum (u_n+1 - u_n)^T (F_n + F_n+1) / 2`, the trapezoidal rule satisfies
`E_0 + W - T - U - D = 0` at every step, exactly: `u1 - u0 = dt/2 (v0 + v1)`
and `v1 - v0 = dt/2 (a0 + a1)` turn the change of `T + U` into the
trapezoidal work less the damping term. The run reports the largest
`|E_0 + W - T - U - D|` over the path, over the largest of the energies:
round-off for `alpha = 0` on a linear model (below 1e-10 in every study; the
same recursion in 80-bit arithmetic brings it down 1.7e3 times, as eps
falls 2.0e3 times), and for `alpha < 0` the energy the method itself
dissipates.

**Non-linear transient.** With `nonlinear: true` the elastic forces become
the internal forces of the load case's non-linear system (section 7c: finite
kinematics with the Saint Venant-Kirchhoff or neo-Hookean law, or small
strain; the J2 plasticity of section 7d with its mean dilatation; follower
pressure), and every step solves the HHT-alpha residual

```
  R(u1) = M a1 + (1 + alpha)(C v1 + f_int(u1) - f_ext(t1)) - alpha (C v0 + f_int(u0) - f_ext(t0)) = 0
```

by Newton's method with the tangent `c0 M + (1 + alpha)(c3 C + K_T)`, from a
constant-acceleration predictor. `C` takes the linear elastic stiffness of
the undeformed model. A step converges when the residual is below the
tolerance times the largest force involved (inertia, loads, damping,
reactions) and the last correction below the tolerance times the larger of
the step's increment and the displacement - an increment alone would ask for
corrections below the rounding of `u` at a turning point of the motion - or
when the residual sits at its round-off floor (section 7c). The plastic
history is committed only on convergence, so a step whose Newton iteration
fails is repeated in halves from the last converged state (at most
`max_cuts` halvings deep), the halves growing back to `dt`. The trapezoidal
rule does not conserve the energy of a non-linear system exactly; its error
is `O(dt^2)`, and the balance also holds the plastic dissipation - the
plastic work less the stored hardening energy - so the run warns only when
energy is created (or, without plasticity, lost) beyond 1 % of the energies
involved: the step is then too long, or a load applied suddenly at a node has
crushed the element under it onto a spurious branch of the Saint
Venant-Kirchhoff law (apply such a load with a table amplitude). The
non-linear transient starts at rest; a preloaded start is a static analysis of
its own.

**Harmonic response.** For a load `f cos(omega t)` and a prescribed motion
`g cos(omega t)`, the steady state is `Re(U e^{i omega t})` with

```
  [K (1 + i eta) - omega^2 M + i omega (a M + b K)] U = f   on the free DOFs, U_p = g ,
```

`eta` the structural (hysteretic) loss factor. Every frequency is a complex
sparse LU (the sparsity pattern analysed once); its backward error
`|| A U - f || / || |A| |U| + |f| ||` must stay below 1e-10, and a response
more than a million times the static one is flagged - at an undamped natural
frequency the dynamic stiffness is singular to working precision and the
answer is round-off. A monitor records the complex amplitude of a
displacement, of a velocity `i omega U` or an acceleration `-omega^2 U`, or of
a reaction; its phase `arg U` is measured against the load (0 in phase,
180 deg in antiphase), and a lightly damped mode turns it by -180 deg through
its resonance, -90 deg of that at the peak. The
largest displacement of a node over a cycle is the semi-major axis of the
ellipse `a cos(omega t) - b sin(omega t)` it traces (`U = a + i b`):
`sqrt((|a|^2 + |b|^2)/2 + sqrt(((|a|^2 - |b|^2)/2)^2 + (a . b)^2))`, the
modulus `sqrt(|a|^2 + |b|^2)` only when its components move in phase.

## 8. Topology optimisation

See `docs/topology_optimization.md` for the SIMP interpolation, the filters, the
sensitivity derivation, the optimality-criteria and MMA updates, and the
aggregated stress constraint with its adjoint.

## 9. Cross-validation

The discrete problem a deck defines is exported verbatim - the same nodes,
connectivity, materials and supports; point loads and tractions as their
consistent nodal forces; pressures, self-weight, body forces, rotation and
temperatures in CalculiX's own form (`P` faces, `GRAV`, `BX/BY/BZ`, `CENTRIF`,
`*EXPANSION` and `*TEMPERATURE`) so that CalculiX integrates them itself, and a
conducted temperature as its own `*HEAT TRANSFER` job - to CalculiX (`*.inp`,
elements CPS4/CPE4, CPS3/CPE3, C3D8, C3D4, C3D10) and rebuilt in scikit-fem
(`ElementQuad1`, `ElementTriP1`, `ElementHex1`, `ElementTetP1`, and
`ElementTetP2` on an isoparametric `MeshTet2` with the 4-point rule, with the
same Lame constants, `lambda* = 2 lambda G / (lambda + 2G)` for plane
stress), and the three nodal displacement fields are compared node by node.
Where the run computed buckling load factors, scikit-fem assembles the
geometric stiffness of its own static solution at the same quadrature
points and solves the pencil densely, and CalculiX runs the exported deck
as a `*BUCKLE` step; the load factors are compared mode by mode. CalculiX's
plane elements are not plane elements internally: it expands them into a
layer of solid elements, which reproduces plane stress only for `nu = 0`, so
a plane-stress comparison at `nu != 0` compares two idealisations and is
recorded without being judged; the plane-strain expansion is exact. scikit-fem
also integrates the pressure, volume and thermal loads itself from the same
deck CalculiX reads. Three of CalculiX's own formulation choices differ from
SparLab's and were identified by reproducing them in scikit-fem: the
element-average temperature for the thermal strain of a first-order
hexahedron, and the four-point (body load) and three-point (face pressure)
rules of the C3D10, which are not exact for a centrifugal load or for a
pressure on a curved face. `docs/verification.md` has the measured differences
and what they mean.

A transient run is integrated again by scikit-fem - its own `K` and `M`, the
same lumping, HHT-alpha in the acceleration (predictor-corrector) form rather
than SparLab's displacement form, and for a non-linear run an independent
J2 / Saint Venant-Kirchhoff system in Newton's method at every step - and,
on solid elements, by CalculiX as a `*DYNAMIC, DIRECT, ALPHA` step with the
amplitude tabulated at every step and Rayleigh `*DAMPING`; the monitors at
every step, the snapshot fields and the final state are compared. CalculiX's
expanded plane elements are not used in dynamics: their step response
contradicts CalculiX's own `*FREQUENCY` result for the same mesh (measured).
A harmonic response is compared with a direct complex solve in scikit-fem;
CalculiX's `*STEADY STATE DYNAMICS` is a modal superposition.
