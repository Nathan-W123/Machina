# Input-deck reference

A run is fully determined by one JSON deck. Every numeric field is SI
(`docs/conventions.md`). Unknown keys are **reported**, and with
`--strict-config` they are an error - a misspelled key that silently takes its
default is one of the easiest ways to publish a wrong number.

Two conveniences beyond strict JSON: `//` and `/* */` comments, and a trailing
comma before `}` or `]`. Everything else follows RFC 8259, and a syntax error
reports its line and column.

## Top level

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `name` | string | `"case"` | case name; also the default output directory |
| `description` | string | `""` | free text, echoed into `summary.json` |
| `units` | string | `"SI"` | informational only |
| `mesh` | object | required | see below |
| `material` | object | required | see below |
| `material_regions` | array | `[]` | other materials on element regions; see below |
| `model` | object | `{}` | idealisation and integration |
| `boundary_conditions` | array | required, non-empty | displacement constraints |
| `load_cases` | array | required, non-empty | loading |
| `solver` | object | `{}` | linear solver and equilibrium tolerances |
| `modal` | object | `{}` | free-vibration analysis |
| `buckling` | object | `{}` | linear buckling check of the load cases |
| `nonlinear` | object | `{}` | non-linear statics: large deflection, finite strain, plasticity |
| `transient` | object | `{}` | transient dynamics (HHT-alpha), linear or non-linear |
| `frequency_response` | object | `{}` | steady harmonic response |
| `topology` | object | `{}` | topology optimisation |
| `output` | object | `{}` | which artefacts to write |

## `mesh`

```json
"mesh": { "type": "structured_quad", "nx": 240, "ny": 160,
          "lx": 0.30, "ly": 0.20, "x0": 0.0, "y0": 0.0 }

"mesh": { "type": "structured_hex", "nx": 32, "ny": 16, "nz": 8,
          "lx": 0.24, "ly": 0.12, "lz": 0.06 }

"mesh": { "type": "file", "path": "../meshes/engine_mount_3d.inp", "scale": 0.001 }

"mesh": { "type": "structured_tet", "nx": 20, "ny": 2, "nz": 2,
          "lx": 1.0, "ly": 0.05, "lz": 0.05, "order": 2 }
```

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `type` | string | `structured_quad` | `structured_quad` (plane, Q4), `structured_tri` (plane, Tri3: each cell split into two triangles), `structured_hex` (solid, Hex8), `structured_tet` (solid, Tet4: each cell split into six Kuhn tetrahedra), or `file` (read from a mesh file) |
| `nx`, `ny` | integer | required for the structured types | elements per direction, `>= 1` |
| `nz` | integer | required for `structured_hex` / `structured_tet` | elements through the third direction |
| `lx`, `ly` | number | required for the structured types | domain extents [m], `> 0` |
| `lz` | number | required for `structured_hex` / `structured_tet` | extent in z [m] |
| `x0`, `y0`, `z0` | number | `0` | lower corner [m] |
| `order` | integer | `1` | `2` makes the tetrahedra quadratic (Tet10): on `structured_tet`, and on a `file` mesh of linear tetrahedra, which is elevated with an edge node at every edge midpoint (straight-sided cells, named sets carried over). A file of 10-node tetrahedra is read as Tet10 without it; `order: 1` on such a file, or `2` on any other cell type, is an error |

The mesh type fixes the dimension of the whole deck: a solid mesh has three
displacement components per node, regions may use `zmin`/`zmax` and
`sphere`, boundary conditions may fix `"z"`, forces and tractions carry three
components, `model.thickness` must be absent (a solid has none) and
`model.stress_state` defaults to `three_dimensional`. Naming `nz`, `lz` or
`z0` on a plane mesh is an error, as is any z component on a plane deck.

**Meshes read from a file** (`"type": "file"`)

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `path` | string | required | the mesh file; a relative path resolves against the deck's directory |
| `format` | string | `auto` | `auto` (from the extension: `.msh` is Gmsh, `.inp` is Abaqus / CalculiX), `gmsh` or `abaqus` |
| `scale` | number | `1` | factor applied to every coordinate: `0.001` for a mesh drawn in millimetres |
| `merge_duplicate_nodes` | bool | `false` | merge coincident nodes instead of only reporting them |
| `duplicate_tolerance` | number [m] | `0` | distance below which two nodes coincide, after scaling; `0` means `1e-9` times the bounding-box diagonal |

Gmsh MSH 2.2 and 4.1 ASCII and Abaqus / CalculiX `.inp` files are read. The
cells of the highest dimension in the file become the mesh: linear triangles
or quadrilaterals give a plane model, linear tetrahedra or hexahedra a solid
one; 10-node tetrahedra (Abaqus `C3D10`, Gmsh type 11, whose last two edge
nodes are swapped into SparLab's order on reading) give a quadratic solid
whose edge nodes may follow curved surfaces. A file that mixes cell types,
or holds any other second-order cell, prisms or pyramids, is refused with a
message saying how to re-export it. Lower
dimensional elements (the boundary curves and facets a mesher writes for its
physical groups) only feed the named sets. Each Gmsh physical group becomes a
node set, and an element set when it holds cells; an Abaqus `*NSET` is a node
set, and an `*ELSET` is an element set when it holds cells and a node set
when it holds boundary elements. Inverted or mirrored cells are re-ordered
and counted, nodes no cell uses are dropped, coincident nodes are reported,
a plane mesh must lie in `z = 0`, and a domain larger than 20 m or smaller
than 0.1 mm draws a warning about units. Boundary conditions, loads and
materials in an `.inp` file are not imported, and the report names every
ignored keyword. What was read and done is recorded under `mesh.file` in
`summary.json`, with the cell-quality statistics under `mesh.quality`.

The resolution of a file mesh comes from the mesher, so `nx` ... `z0` are
errors on a `file` mesh, and so are `--nx`/`--ny`/`--nz` on the command line.
`python/scripts/make_meshes.py` regenerates the meshes the decks use (`make
meshes`, which needs the `gmsh` Python package), the quadratic engine mount
among them.

## `material`

```json
"material": { "name": "Al 7075-T6", "youngs_modulus": 71.7e9,
              "poisson_ratio": 0.33, "density": 2810.0 }
```

| Key | Type | Default | Constraint |
|-----|------|---------|------------|
| `name` | string | `"material"` | carried into result files |
| `youngs_modulus` | number [Pa] | required | `> 0` |
| `poisson_ratio` | number | required | in `(-1, 0.5)` so both idealisations stay positive definite |
| `density` | number [kg/m^3] | `0` | `>= 0`; a positive value is required for modal analysis, self-weight and rotation |
| `thermal_expansion` | number [1/K] | `0` | linear expansion coefficient `alpha`; `0` means a temperature causes no strain |
| `reference_temperature` | number [K] | `0` | the stress-free temperature `T_ref`: the thermal strain is `alpha (T - T_ref)` |
| `conductivity` | number [W/(m K)] | `0` | `k`, needed by a conducted temperature field (`load_cases[].temperature.conduction`) |
| `plasticity` | object | none | J2 plasticity, below; without it the material is elastic |

Temperatures may be given in kelvin or in degrees Celsius as long as every
temperature of the deck - `reference_temperature`, and the load cases' uniform,
regional, prescribed and ambient values - uses the same scale: only differences
enter the thermal strain and the conduction problem. The thermal strain is
`alpha dT {1, 1, 0}` in plane stress (the free out-of-plane expansion leaves the
in-plane law unchanged), `(1 + nu) alpha dT {1, 1, 0}` in plane strain (the
restrained out-of-plane expansion; `sigma_zz = nu (sigma_xx + sigma_yy) - E alpha dT`
is reported as `sigma_zz`), and `alpha dT {1, 1, 1, 0, 0, 0}` in 3-D.

### `material.plasticity`

```json
"plasticity": { "yield_stress": 250e6, "hardening_modulus": 2e9,
                "saturation_stress": 80e6, "saturation_rate": 25,
                "kinematic_hardening_modulus": 3e9 }
```

Rate-independent plasticity, which the non-linear analysis models (`nonlinear`,
below; `docs/formulation.md`, section 7d): von Mises (J2) or Hill's 1948
anisotropic yield criterion, with isotropic and kinematic hardening. The yield
stress after an accumulated plastic strain `alpha` is
`sigma_y(alpha) = yield_stress + hardening_modulus alpha + saturation_stress (1 - exp(-saturation_rate alpha))`
(linear isotropic hardening plus Voce saturation), and the yield surface
moves with the back stress of Prager's linear kinematic hardening,
`d beta = (2/3) kinematic_hardening_modulus d eps_p`, and/or of Chaboche
backstresses (below). In uniaxial tension the
plastic slope of the stress against the plastic strain is
`hardening_modulus + kinematic_hardening_modulus` (plus the Voce part).

| Key | Type | Default | Constraint |
|-----|------|---------|------------|
| `yield_stress` | number [Pa] | required | `> 0`: the initial uniaxial yield stress `sigma_y0` (with Hill48, the equivalent stress at yield: the uniaxial yield stress along the rolling direction for the r-value and stress-ratio calibrations) |
| `hardening_modulus` | number [Pa] | `0` | `>= 0`, `H`: linear isotropic hardening (`0` with the others `0` is perfect plasticity) |
| `saturation_stress` | number [Pa] | `0` | `>= 0`, `Q`: the Voce saturation increment of the yield stress |
| `saturation_rate` | number [-] | `0` | `> 0` when `saturation_stress` is given, `delta` |
| `kinematic_hardening_modulus` | number [Pa] | `0` | `>= 0`, `H_kin`: Prager's linear kinematic hardening |
| `yield_criterion` | string | `von_mises` | `von_mises` or `hill48`; `hill48` needs the `anisotropy` block and the block needs `hill48` |
| `anisotropy` | object | - | the Hill48 calibration and material frame (below) |
| `backstresses` | array | `[]` | at most 4 Armstrong-Frederick backstresses `{ "modulus": C_i, "recovery": gamma_i }` (below) |
| `kinematic_integration` | string | `exponential` | how a backstress is integrated over a step: `exponential` or `backward_euler` |

**Hill48 anisotropy** (`"yield_criterion": "hill48"`). Sheet metal rolled
along RD, with transverse direction TD and normal ND, yields when
`F (s_TT - s_NN)^2 + G (s_NN - s_RR)^2 + H (s_RR - s_TT)^2 + 2 L s_TN^2 + 2 M s_NR^2 + 2 N s_RT^2 = sigma_y(alpha)^2`
of the relative stress in that frame, with associative (isochoric) flow.

```json
"plasticity": { "yield_stress": 170e6, "saturation_stress": 90e6, "saturation_rate": 12,
  "yield_criterion": "hill48",
  "anisotropy": { "r0": 1.9, "r45": 1.5, "r90": 2.3, "out_of_plane_shear": [1.5, 1.5],
                  "rolling_direction": [1, 0, 0], "sheet_normal": [0, 0, 1] },
  "backstresses": [ { "modulus": 60e9, "recovery": 600 }, { "modulus": 8e9, "recovery": 60 },
                    { "modulus": 1e9, "recovery": 0 } ],
  "kinematic_integration": "exponential" }
```

| `anisotropy` key | Type | Default | Meaning |
|------------------|------|---------|---------|
| `r0`, `r45`, `r90` | numbers [-] | - | Lankford coefficients at 0, 45 and 90 deg to RD, `> 0`: `F = r0 / (r90 (1 + r0))`, `G = 1 / (1 + r0)`, `H = r0 / (1 + r0)`, `N = (r0 + r90)(1 + 2 r45) / (2 r90 (1 + r0))`, so that `sigma_y` is the yield stress along RD |
| `stress_ratios` | object | - | `{ "sigma_45": s45, "sigma_90": s90, "sigma_biaxial": sb }`, the yield stresses at 45 and 90 deg and in equibiaxial tension over the one along RD, `> 0`: `G + H = 1`, `F + H = 1 / s90^2`, `F + G = 1 / sb^2`, `2 N = 4 / s45^2 - F - G` |
| `coefficients` | object | - | `{ "F", "G", "H", "L", "M", "N" }` directly (the yield stress along RD is then `yield_stress / sqrt(G + H)`) |
| `out_of_plane_shear` | `[L, M]` | `[1.5, 1.5]` | the out-of-plane shear coefficients of the r-value and stress-ratio calibrations (in-plane tests do not measure them; `1.5` is isotropic) |
| `rolling_direction` | `[x, y, z]` | `[1, 0, 0]` | RD in global axes, projected normal to ND |
| `sheet_normal` | `[x, y, z]` | `[0, 0, 1]` | ND in global axes; TD = ND x RD |
| `rolling_angle` | number [deg] | - | RD in the x-y plane at this angle from x, ND = z (instead of the two vectors) |

Exactly one of the r-values, `stress_ratios` and `coefficients` is given
(`r = 1`, or all ratios `1`, is von Mises). The coefficients must make a
convex surface: `L, M, N > 0`, `F + H > 0` and `FG + GH + HF > 0`. The frame
directions always have three components, in a plane model too, whose `z`
must then be one of RD, TD, ND (a plane model has no out-of-plane shear
strain): for a sheet in the model plane (plane stress) ND = z, for a
section through the thickness (plane strain) for instance ND = y and TD = z.
With finite kinematics the frame is the reference one and turns with the
material. Hill48 in uniaxial tension at `theta` from RD yields at
`sigma_y [F sin^4 + G cos^4 + H cos^2 2theta + 2 N sin^2 cos^2]^(-1/2)` with the
r-value `[H + (2N - F - G - 4H) sin^2 cos^2] / [F sin^2 + G cos^2]` - both
verified to round-off (`--study hill-directional`).

**Chaboche backstresses** (`backstresses`). The back stress is the sum of up
to four Armstrong-Frederick terms,
`d alpha_i = (2/3) C_i d eps_p - gamma_i alpha_i d alpha` (with Prager's
`kinematic_hardening_modulus` as one more, without recovery - at most four in
all), each saturating at `C_i / gamma_i` in uniaxial tension - for von
Mises, `sigma = sigma_y(alpha) + sum_i (C_i / gamma_i)(1 - exp(-gamma_i alpha))`
on a monotonic path (with Hill48 the backstress grows along the Hill flow
direction, not along the stress, and a uniaxial path is not proportional).
`modulus` is required and `> 0`, `recovery` defaults to `0` (Prager) and is
`>= 0`. With any `recovery > 0` the consistent tangent is not symmetric: the
iterations of the non-linear analysis then factorise it by LU, which reveals
no inertia, so the test for a negative pivot within a step is off (the
tangent at a converged state - the continuum one, symmetric - is still
factorised by LDL^T, and its inertia reported). `exponential`
integrates each backstress exactly for a flow direction fixed over the step
(so a uniaxial path is exact at any step size); `backward_euler` is
first-order accurate.

The CalculiX export writes no non-linear or transient deck for a Hill48
material, for backstresses or for Prager's rule, and says why: CalculiX's
`*PLASTIC` is von Mises, has no Armstrong-Frederick rule, and its
`HARDENING=KINEMATIC` does not reproduce Prager's linear rule.

The linear analyses - static, modal, buckling - and the topology
optimisation treat the material as elastic; `sparlab_solve` says so when a
plastic material meets no non-linear analysis.

## `material_regions`

```json
"material_regions": [
  { "name": "aluminium_half", "region": { "box": { "xmin": 0.1 } },
    "material": { "name": "aluminium", "youngs_modulus": 70.0e9, "poisson_ratio": 0.33,
                  "density": 2700.0, "thermal_expansion": 2.3e-5,
                  "reference_temperature": 293.15, "conductivity": 167.0 } }
]
```

Each entry gives the elements of `region` (selected at their centroids, like
every element region) the `material`, which takes the same keys as the top-level
`material`. Later entries override earlier ones on shared elements; every
other element keeps the top-level material. A region that selects no element
is an error. Stresses, strain energies, the self-weight, the thermal strain and
the conduction all use each element's own material, and the CalculiX export
writes one `*MATERIAL` and `*SOLID SECTION` per material.

Topology optimisation (`sparlab_topopt`) rejects `material_regions` for now:
its interpolation is written for one solid material.

## `model`

```json
"model": { "thickness": 0.008, "stress_state": "plane_stress",
           "integration": { "stiffness_points": 2, "mass_points": 3,
                            "edge_points": 2 } }
```

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `thickness` | number [m] | `1.0` | out-of-plane thickness, `> 0`; plane meshes only - a solid mesh rejects the key |
| `stress_state` | string | `plane_stress` (`three_dimensional` on a solid mesh) | `plane_stress`, `plane_strain` or `three_dimensional`; the plane idealisations need a plane mesh (Q4 or Tri3) and `three_dimensional` a solid one (Hex8, Tet4 or Tet10) |
| `integration.stiffness_points` | integer | `2` | Gauss points per direction for `K_e`, 1-4 (2x2 for the Q4, 2x2x2 for the Hex8); the linear simplices have a constant strain and integrate exactly with one point whatever is set here |
| `integration.mass_points` | integer | `3` | Gauss points per direction for `M_e`, 1-4; the simplices use the exact closed-form consistent mass instead |
| `integration.face_points` | integer | `2` | Gauss points per direction on a loaded edge or face, 1-4; `edge_points` is accepted as a synonym |

## Regions

Boundary conditions, loads and passive regions all take a `region` object. A
region is the union of one or more primitives, optionally complemented.

```json
"region": { "name": "bolt_attachments", "invert": false, "tolerance": 0.0,
            "any_of": [
              { "circle": { "center": [0.030, 0.045], "radius": 0.0105 } },
              { "circle": { "center": [0.030, 0.155], "radius": 0.0105 } }
            ] }
```

A single primitive may also be written directly in the `region` object, without
`any_of`.

| Primitive | Form | Selects |
|-----------|------|---------|
| `all` | `"all": true` | everything |
| `box` | `{ "xmin":, "xmax":, "ymin":, "ymax":, "zmin":, "zmax": }`, any subset | points inside the axis-aligned box; omitted sides are unbounded; `zmin`/`zmax` on a solid mesh only |
| `circle` | `{ "center": [x, y(, z)], "radius": r, "axis": "z" }` | the closed disc; on a solid mesh, the infinite cylinder of that radius about the line through `center` along `axis` (`x`, `y` or `z`, default `z`) |
| `annulus` | `{ "center": [...], "inner_radius": r0, "radius": r1, "axis": "z" }` | the closed ring, `r0 < r1`; a tube on a solid mesh |
| `sphere` | `{ "center": [x, y, z], "radius": r }` | the closed ball (solid meshes only) |
| `node_ids` | `[...]` | explicit nodes (node regions only) |
| `element_ids` | `[...]` | explicit elements (element regions only) |
| `nearest_node` | `[x, y(, z)]` | the single node closest to the point (node regions only) |
| `group` | `"name"` | a named set of a mesh read from a file: its node set in a node region, its element set in an element region |

A `center` or `nearest_node` point carries as many components as the mesh has
dimensions. A `group` that the mesh does not have is an error listing the
sets it does have; so is a `group` on a structured mesh, which has none.

| Region key | Type | Default | Meaning |
|------------|------|---------|---------|
| `name` | string | auto | used in diagnostics |
| `invert` | bool | `false` | select the complement of the union |
| `tolerance` | number [m] | `0` | absolute geometric tolerance; `0` means `1e-9` times the mesh diagonal |

A degenerate box such as `{"xmax": 0.0}` reliably captures the line of nodes at
`x = 0` because of that relative tolerance. Element regions are tested at the
element centroid. Naming two primitives in one region object without `any_of` is
an error, as is a region that selects nothing.

## `boundary_conditions`

```json
"boundary_conditions": [
  { "name": "clamped_root", "fix": ["x", "y"], "value": [0.0, 0.0],
    "region": { "box": { "xmax": 0.0 } } }
]
```

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `name` | string | auto | diagnostics |
| `fix` | array of `"x"` / `"y"` (/ `"z"` on a solid mesh) | required, non-empty | which components to prescribe |
| `value` | `[u_x, u_y(, u_z)]` [m] | zeros | prescribed displacement |
| `region` | object | required | node region |

Non-zero values are handled by static condensation, not by modifying the load
vector, so the reactions stay exact.

## `load_cases`

```json
"load_cases": [
  { "name": "down_limit", "weight": 1.0,
    "point_loads": [
      { "name": "lug", "force": [0.0, -9000.0], "distribution": "total",
        "region": { "annulus": { "center": [0.270, 0.100],
                                 "inner_radius": 0.0075, "radius": 0.0105 } } }
    ],
    "tractions": [
      { "name": "upper_skin", "traction": [0.0, -1.2e5],
        "region": { "box": { "ymin": 0.10 } } }
    ]
  }
]
```

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `name` | string | auto | used in file names and reports |
| `weight` | number | `1.0` | weight in the multi-load objective, `>= 0`; weights are normalised to sum to 1 |
| `point_loads` | array | `[]` | concentrated nodal forces |
| `tractions` | array | `[]` | distributed edge loads |
| `pressures` | array | `[]` | normal pressures on boundary edges or faces |
| `gravity` | `[g_x, g_y(, g_z)]` [m/s^2] | none | uniform acceleration of the whole model: the self-weight `rho g` of every element |
| `body_forces` | array | `[]` | uniform force densities on element regions |
| `centrifugal` | object | none | steady rotation about an axis |
| `temperature` | object | none | the case's temperature field, whose thermal strain loads the structure |
| `prescribed_displacement_only` | bool | `false` | declares a case with no applied force, driven by the prescribed displacements (the patch test, enforced-deflection studies) |

**Point load**

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `force` | `[F_x, F_y(, F_z)]` [N] | required | see `distribution` |
| `distribution` | `"total"` / `"per_node"` | `"total"` | `total` divides the resultant among the selected nodes, so the total force is mesh independent; `per_node` applies `force` to each node |
| `region` | object | required | node region |

**Traction**

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `traction` | `[t_x, t_y(, t_z)]` [Pa] | required | stress vector in **global** components, not a normal pressure |
| `region` | object | required | node region; a boundary edge (face, on a solid mesh) is loaded when *all* its nodes lie inside |

Tractions are integrated to consistent nodal forces, so the resultant is exactly
`traction * (loaded length) * thickness` on a plane mesh and
`traction * (loaded area)` on a solid one, at any mesh resolution. A traction
region matching no boundary edge or face is an error.

**Pressure**

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `pressure` | number [Pa] | required | positive pushes *into* the body |
| `region` | object | required | node region; like a traction, a boundary edge or face is loaded when all its nodes lie inside |

A pressure acts along the normal of the face at every point: it is integrated
with the face's own area vector (`x_s x x_t` on a surface, the rotated tangent
on an edge), so on a curved Tet10 face - whose edge nodes lie on the true
geometry - it follows the curvature, and its resultant on a closed surface is
zero to round-off. On a curved six-node face it is integrated with three
points per direction, exact for the degree-4 integrand.

**Body force**

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `force_density` | `[b_x, b_y(, b_z)]` [N/m^3] | required | force per unit volume |
| `region` | object | every element | element region (centroids) |

**Centrifugal**

```json
"centrifugal": { "angular_velocity": 300.0, "axis": [0.0, 0.0, 1.0], "point": [0.0, 0.05, 0.0] }
```

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `angular_velocity` | number [rad/s] | required | `omega` |
| `axis` | `[a_x, a_y, a_z]` | `[0, 0, 1]` | direction of the axis (normalised); a plane model can only spin about `z` |
| `point` | `[x, y(, z)]` [m] | origin | a point on the axis |

Every element carries `rho omega^2 r_perp`, `r_perp` its distance vector from
the axis. The body loads are integrated from the consistent mass,
`f = M_e(rho = 1) b(x_nodes)`, which is exact for any force density affine in
position - self-weight, a uniform body force and the centrifugal load - on
straight and curved cells alike (`include/sparlab/fem/Loads.hpp`).

**Temperature**

Three ways to give a case a temperature field [K]; the thermal strain of each
element uses its own material's `thermal_expansion` and `reference_temperature`.

```json
"temperature": { "uniform": 343.15 }

"temperature": { "uniform": 293.15,
                 "regions": [ { "name": "end", "value": 393.15,
                                "region": { "box": { "xmin": 0.15 } } } ] }

"temperature": { "conduction": {
    "prescribed": [ { "name": "root", "value": 373.15, "region": { "box": { "xmax": 0.0 } } } ],
    "flux":       [ { "name": "end",  "value": 1.0e4,  "region": { "box": { "xmin": 0.2 } } } ],
    "sources":    [ { "name": "heater", "value": 2.0e5, "region": { "box": { "xmin": 0.15 } } } ],
    "convection": [ { "name": "top", "film_coefficient": 25.0, "ambient": 293.15,
                      "region": { "box": { "zmin": 0.02 } } } ] } }
```

| Key | Meaning |
|-----|---------|
| `uniform` [K] | every node at this temperature; with `regions` the base value, which defaults to the primary material's `reference_temperature` |
| `regions` | node regions and their temperatures; later entries win on shared nodes |
| `conduction.prescribed` | fixed temperatures [K] on node regions |
| `conduction.flux` | heat flux into the body [W/m^2] on boundary edges / faces |
| `conduction.sources` | volumetric heat generation [W/m^3]; `region` (element centroids) is optional and defaults to every element |
| `conduction.convection` | `q = film_coefficient (T - ambient)` leaving the body [W/(m^2 K), K] on boundary edges / faces |

`conduction` solves steady heat conduction, `-div(k grad T) = Q`, on the same
mesh with the same shape functions (a plane model conducts in its plane with its
thickness and insulated faces), and the solved temperature drives the thermal
strain. It needs a prescribed temperature or convection somewhere: with fluxes
and sources alone the temperature is fixed only up to a constant, which is
refused. The run logs the heat balance - heat entering through sources, fluxes
and convection against heat leaving through the prescribed temperatures - and
stops with `SolverError` if it misses by more than `1e-6` of the largest heat
flow. `conduction` cannot be combined with `uniform` or `regions`. The solved
temperatures are written to `temperature_<case>.csv` and to the VTK file.

A temperature changes the applied load by `f_th = int B^T D eps_0 dV`,
integrated with the element's stiffness rule, and the elastic strain energy
becomes `1/2 u^T K u - u^T f_th + 1/2 int eps_0^T D eps_0 dV`; the stresses are
`D (B u - eps_0)`. A linear buckling analysis includes the thermal prestress.

A load case with no load at all - no point load, traction, pressure, gravity,
body force, rotation or temperature - is an error unless
`prescribed_displacement_only` is set.

## `solver`

```json
"solver": {
  "linear": { "type": "auto", "residual_tolerance": 1e-8,
              "iterative_tolerance": 1e-12, "max_iterations": 0,
              "pivot_tolerance": 1e-14, "warm_start": true,
              "auto_direct_limit": { "plane": 50000, "solid": 10000 },
              "amg": { "smoother": "chebyshev", "smoother_degree": 3,
                       "strength_threshold": 0.02, "coarse_size": 1500 } },
  "equilibrium_tolerance": 1e-6,
  "check_model": true
}
```

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `linear.type` | string | `auto` | `auto`, `simplicial_ldlt`, `amg_cg` (also `amg`, `multigrid`), `conjugate_gradient`, `simplicial_llt`, `sparse_lu`, `dense_lu` |
| `linear.residual_tolerance` | number | `1e-8` | scaled residual `||K u - f|| / ||f||` accepted after each solve, whatever the solver. A solve whose backward error `||K u - f|| / || |K| |u| + |f| ||` is at round-off (64 machine epsilon) is accepted too: its residual is the rounding of the sums that form `K u`, which in a slender structure can exceed `1e-8 ||f||` (`docs/formulation.md`, section 5) |
| `linear.iterative_tolerance` | number | `1e-12` | relative residual the CG solvers iterate to |
| `linear.max_iterations` | integer | `0` | CG iteration cap; `0` means `2n` for Jacobi CG and 1000 for multigrid CG |
| `linear.pivot_tolerance` | number | `1e-14` | smallest accepted `min/max` LDL^T pivot ratio |
| `linear.warm_start` | bool | `true` | start each iterative solve from the previous solution of the same load case (or adjoint); an optimisation loop then needs a fraction of the iterations |
| `linear.auto_direct_limit.plane` | integer | `50000` | `auto` factorises a plane system with at most this many free unknowns and uses multigrid CG above it |
| `linear.auto_direct_limit.solid` | integer | `10000` | the same limit for a solid system, where the Cholesky fill grows much faster |
| `equilibrium_tolerance` | number | `1e-6` | relative force-balance error before `SolverError` |
| `check_model` | bool | `true` | run the rigid-body and floating-region diagnostics before assembling |

The solvers: `simplicial_ldlt` is Eigen's sparse LDL^T with an AMD ordering,
the reference; `amg_cg` is conjugate gradients preconditioned by SparLab's
smoothed-aggregation algebraic multigrid (`docs/formulation.md` gives the
construction); `conjugate_gradient` is Eigen's Jacobi-preconditioned CG, kept
as the baseline that shows what the multigrid buys; `simplicial_llt`,
`sparse_lu` and `dense_lu` exist for the solver-agreement verification. The
summary records which solver ran (`auto:AMG-CG(chebyshev)`, for instance),
its settings and, for the iterative ones, the iteration counts.

**`linear.amg`** - the multigrid preconditioner. The defaults are the
measured ones; change them to study the solver, not to make a run work.

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `strength_threshold` | number | `0.02` | node pairs whose coupling block is weaker than this (relative, Frobenius) are not aggregated together |
| `max_levels` | integer | `12` | hierarchy depth cap |
| `coarse_size` | integer | `1500` | unknowns at or below which a level is solved directly (dense LDL^T) |
| `smoother` | string | `chebyshev` | `chebyshev` or `gauss_seidel` (symmetric) |
| `smoother_degree` | integer | `3` | Chebyshev degree, or Gauss-Seidel sweeps |
| `chebyshev_ratio` | number | `30` | the smoother targets eigenvalues in `[lambda_max / ratio, lambda_max]` |
| `prolongator_damping` | number | `4/3` | Jacobi damping of the tentative prolongator, divided by `lambda_max(D^-1 A)` |
| `lanczos_steps` | integer | `12` | Lanczos steps of the `lambda_max` estimate |
| `reuse_aggregates` | bool | `true` | keep the aggregation (and the sparsity of every product) across refactorisations of a matrix with the same pattern, as in an optimisation loop |
| `coarse_pivot_tolerance` | number | `1e-13` | a coarsest-level pivot below this, relative to the largest, is reported as an under-constrained model |

## `modal`

```json
"modal": { "enabled": true, "num_modes": 6, "mass_type": "consistent",
           "tolerance": 1e-10, "residual_tolerance": 1e-6,
           "max_iterations": 200, "shift_factor": 0.0,
           "rigid_body_ratio": 1e-10, "seed": 20240917,
           "analyse_optimised_topology": true,
           "compare_mass_matched_baseline": true }
```

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `enabled` | bool | `false` | run free-vibration analysis |
| `num_modes` | integer | `6` | lowest modes requested, `>= 1` |
| `mass_type` | string | `consistent` | or `lumped` (row-sum) |
| `tolerance` | number | `1e-10` | relative eigenvalue change for convergence |
| `residual_tolerance` | number | `1e-6` | eigenpair residual, also part of the stopping rule |
| `max_iterations` | integer | `200` | subspace-iteration cap |
| `shift_factor` | number | `0` | shift as a multiple of `min_i K_ii/M_ii`; `0` disables |
| `rigid_body_ratio` | number | `1e-10` | below this times the stiffness/mass scale, a mode is flagged as rigid-body |
| `seed` | integer | `20240917` | seed for the filler part of the starting basis |
| `analyse_optimised_topology` | bool | `true` | also compute *mode shapes* of the thresholded solid interpretation |
| `compare_mass_matched_baseline` | bool | `true` | topology runs also analyse an equal-mass uniform plate |

## `buckling`

```json
"buckling": { "enabled": true, "num_modes": 4, "tolerance": 1e-8,
              "residual_tolerance": 1e-6, "max_iterations": 400,
              "seed": 20240917, "load_cases": ["axial"],
              "analyse_optimised_topology": true }
```

A linear (bifurcation) buckling analysis of each load case's linear static
state, `(K + lambda K_G(u)) phi = 0` (`docs/formulation.md`, section 7b).
`sparlab_solve` runs it on the analysed model; `sparlab_topopt` runs it
after the optimisation on the full solid domain and on the exported part.

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `enabled` | bool | `false` | run the buckling analysis |
| `num_modes` | integer | `4` | lowest positive load factors requested, `>= 1` |
| `tolerance` | number | `1e-8` | relative change of the requested load factors for convergence |
| `residual_tolerance` | number | `1e-6` | largest accepted eigenpair residual `||K phi + lambda K_G phi|| / ||K phi||` |
| `max_iterations` | integer | `400` | subspace-iteration cap |
| `seed` | integer | `20240917` | seed of the random part of the starting subspace |
| `load_cases` | array of strings | all | load cases to analyse, by name |
| `analyse_optimised_topology` | bool | `true` | in a topology run, also analyse the exported part |

A load factor above 1 is a safety factor on the load case; a load case that
only stretches the structure has no positive load factor and is reported
as such. The results are `buckling_<solid|topology>.csv` with every mode's
load factor, residual and energy share in solid material, and, with
`output.vtk` and `output.mode_shapes`, one VTK file per mode.

## `nonlinear`

```json
"nonlinear": {
  "enabled": true,
  "kinematics": "finite",
  "material_model": "saint_venant_kirchhoff",
  "method": "load_control",
  "steps": 10,
  "residual_tolerance": 1e-8,
  "follower_pressure": true,
  "load_factors": [0.25, 0.5],
  "monitors": [
    { "name": "tip_deflection", "component": "y", "region": { "box": { "xmin": 1.0 } } },
    { "name": "root_force", "component": "y", "quantity": "reaction",
      "region": { "box": { "xmax": 0.0 } } }
  ],
  "load_cases": ["tip_force"]
}
```

A non-linear static analysis of each selected load case: large displacement
and rotation in the total Lagrangian form (`kinematics: finite`), large
strain too in the logarithmic strain (`finite_logarithmic`), or small strain
(`small_strain`), with plasticity for the materials that have a
`plasticity` block, the equilibrium solved by Newton's method along a load
path (`docs/formulation.md`, sections 7c and 7d). One load factor lambda
scales every load of the case together - forces, pressures, self-weight and
body forces, the rotation, the temperature change and prescribed
displacements - from 0 to 1, or along a `load_path` that may unload.
`sparlab_solve` runs it after the linear analysis, which stays in the summary
for comparison; `sparlab_topopt` refuses a deck with it enabled (the
optimisation is linear).

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `enabled` | bool | `false` | run the non-linear analysis |
| `kinematics` | string | `finite` | `finite`: the total Lagrangian formulation in the Green-Lagrange strain, large displacement and rotation. `finite_logarithmic`: the total Lagrangian formulation in the logarithmic (Hencky) strain `ln(C)/2` - large strain: the plasticity is the small-strain return in the log strain (exact for coaxial stretches of any size), the elasticity Hencky's, for elastic materials too; the thermal strain is `ln(1 + alpha dT)`; mean dilatation averages `ln J`; the summary reports the largest log strain and the stress CSV and VTK add the Kirchhoff stress and the log strain (`docs/formulation.md`, section 7d). `small_strain`: the linear strain on the undeformed geometry - no geometric stiffness, pressures on the undeformed faces, the rotation's load at the undeformed positions; with elastic materials it is the linear analysis, with plastic ones the classical small-strain elastoplastic analysis |
| `material_model` | string | `saint_venant_kirchhoff` | the elastic law with `finite` kinematics. `saint_venant_kirchhoff`: the linear law between Green-Lagrange strain and second Piola-Kirchhoff stress - large rotation, small strain, any stress state. `neo_hookean`: the compressible neo-Hookean law - large strain, plane strain or solid meshes only, elastic materials only. A plastic material takes the Saint Venant-Kirchhoff form (its J2 return in the Green strain). `finite_logarithmic` kinematics has its own law (Hencky) and refuses `neo_hookean` |
| `mean_dilatation` | string or bool | `auto` | the elements of a plastic material that average their dilatation over the element (B-bar; its Green-strain form with `finite` kinematics, the average of `ln J` with `finite_logarithmic`), which keeps them from locking under the isochoric plastic flow. `auto`: Q4 and Hex8, which lock without it; `all`: also Tet10; `none`. `true` and `false` stand for `all` and `none`. Plane stress and one-point elements (Tri3, Tet4) have nothing to average |
| `method` | string | `load_control` | `load_control`: Newton at prescribed load factors. `arc_length`: Crisfield's cylindrical arc-length method, which follows the path through limit points |
| `steps` | integer | `10` | load control: the equal steps to lambda = 1 it starts with (halved on failure, lengthened again after easy steps, never beyond this size). Arc length: the first arc length is that of the first of `steps` equal load increments |
| `max_steps` | integer | `500` | converged steps before the run stops |
| `max_iterations` | integer | `25` | Newton iterations per step |
| `max_cuts` | integer | `12` | successive halvings of a failing step before the run stops |
| `residual_tolerance` | number | `1e-8` | out-of-balance force over the load scale (the largest of the applied loads, the reactions and the thermal forces) |
| `displacement_tolerance` | number | `1e-8` | Newton correction over the displacement increment of the step |
| `line_search` | bool | `true` | energy line search along each Newton direction |
| `follower_pressure` | bool | `true` | pressures act on the deformed faces, their direction and area following the deformation; `false` keeps them on the reference faces (dead) |
| `load_factors` | array of numbers | none | load control: load factors in (0, 1), strictly increasing, that the path passes through exactly, so the results are recorded at chosen load levels. Not with `arc_length` |
| `load_path` | array of numbers | none | load control: the load factors to visit in turn from 0, each reached exactly - a path with turning points, such as `[1, 0]` (load, then unload: the permanent set and residual stresses of a plastic body) or `[1, -1, 1]` (a cycle). Each leg starts with `steps` equal steps. Not with `load_factors` or `arc_length` |
| `target_load_factor` | number | `1` | arc length: the load factor where the run stops, reached exactly by a last load-controlled step |
| `desired_iterations` | integer | `5` | arc length: the Newton iterations per step the arc length adapts towards |
| `min_arc_ratio`, `max_arc_ratio` | number | `1e-6`, `10` | arc length: bounds on the arc length relative to the first one |
| `monitors` | array | none | quantities recorded at every converged step. Each has a `name`, a node `region`, a `component` (`x`, `y`, `z`) and a `quantity`: `displacement` (default) is the mean over the region's nodes [m], `reaction` the sum of the support reactions over them [N] |
| `load_cases` | array of strings | all | load cases to analyse, by name |

A load case whose path meets a limit or bifurcation point cannot be followed
by load control past it. The analysis rejects a step whose Newton iterations,
from a stable state, meet a tangent with a negative pivot, or whose converged
state lies farther from the tangent predictor than the predicted increment
itself (a jump to another branch). It closes in on the point by halving and
stops when `max_cuts` halvings have failed to pass it. `summary.json` then
records `completed: false`, the reason and the bracket of the critical load
factor (`critical_load_factor_bracket`), and the console line says `STOPPED`.
The arc-length method follows such a path, and records the negative pivots
of each converged state, which mark the unstable stretches.

With plasticity the second of those tests (the jump to another branch) is
not made: the elastic predictor of a step in which points start to yield
underestimates the increment by the ratio of elastic to plastic stiffness,
and small-strain kinematics has no distant branches. Beyond a plastic
collapse load - a structure of a material without hardening - there is no
equilibrium at all: load control stops at the last converged load factor
with `unreached_load_factor`, the nearest load factor no step reached, or
the bracket above when the singular collapse tangent shows a negative pivot.
The largest converged load factor is a lower bound on the collapse load (an
equilibrium within yield); the arc-length method runs onto the collapse
plateau, whose load factor is the collapse load of the model
(`docs/verification.md`, the plastic cylinder).

Refused with a message before any step: the neo-Hookean law in plane stress
(its plane-stress form, which needs the thickness stretch that makes S_33
vanish, is not implemented), the neo-Hookean law under a temperature change
in a material with thermal expansion (the thermal strain is modelled for
Saint Venant-Kirchhoff only), the neo-Hookean law with `small_strain`
kinematics or a plastic material, `follower_pressure: true` with
`small_strain` kinematics, a `load_path` with `load_factors` or with the
arc-length method, and the arc-length method with a non-zero prescribed
displacement (its constraint measures the free displacements only; drive a
prescribed displacement by load control, which the arch study does to follow
a snap-through in displacement control).

The results, per load case:

| File | Content |
|------|---------|
| `nonlinear_<lc>.csv` | one row per converged step: load factor, Newton iterations, halvings, the relative residual at convergence, the arc length, the negative pivots of the tangent (`-1` when a non-symmetric tangent was factorised by LU, which reports no inertia), the largest displacement, with plasticity the integration points that yielded in the step and the largest accumulated plastic strain after it, and the monitors |
| `nonlinear_displacement_<lc>.csv`, `nonlinear_reactions_<lc>.csv` | nodal displacements and reactions at the final load factor (with `output.csv`) |
| `nonlinear_stress_<lc>.csv` | element Cauchy stress (with `szz` in plane strain), its von Mises stress, the second Piola-Kirchhoff stress (with small strain both are sigma) and, with plasticity, the largest accumulated plastic strain of the element's points (with `output.csv`) |
| `nonlinear_<lc>.vtk` | displacement, Cauchy stress, von Mises stress and, with plasticity, the equivalent plastic strain at the final state (with `output.vtk`) |

The `nonlinear` block of `summary.json` holds, per load case, whether it
reached lambda = 1, the steps, iterations and halvings, the largest
displacement beside the linear analysis's, the strain energy, the largest
Green-Lagrange strain, the smallest volume ratio J, the largest von Mises
stress, the final monitors, the force and moment balance of the deformed
body, which factorisation the tangent took, and warnings - for a run that
stopped, a final state that is not stable, a non-symmetric tangent whose
stability is therefore not assessed, and a Saint Venant-Kirchhoff state
whose strains exceed the law's range (Green strain above 0.05, or a volume
ratio J below 1/sqrt(3): below a stretch of 1/sqrt(3) the law's compressive
force falls again). With `small_strain` kinematics it records the largest
strain, the largest rotation and the largest component of the neglected
quadratic strain `H^T H / 2`, and warns when that exceeds a tenth of the
largest strain (it matters where the structure is restrained against the
motion a rotation causes: membrane action) or a strain exceeds 0.05. With
plasticity its `plasticity` block records the integration points that have
yielded, the largest accumulated plastic strain, the elements that yielded,
whether mean dilatation was applied, and the step in which a point first
yielded; the warnings add strains beyond 0.05 (the elastoplastic law assumes
small strains) and elements that lock under isochoric plastic flow.

## `transient`

```json
"transient": {
  "enabled": true,
  "time_step": 1e-4,
  "end_time": 3e-2,
  "alpha": -0.05,
  "mass": "consistent",
  "damping": { "mass": 20.0, "stiffness": 1e-5 },
  "amplitude": { "type": "table", "times": [0, 0.004], "values": [0, 1] },
  "start": "rest",
  "snapshot_every": 30,
  "monitors": [
    { "name": "tip_uy", "component": "y", "region": { "box": { "xmin": 1.0 } } },
    { "name": "tip_vy", "component": "y", "quantity": "velocity",
      "region": { "box": { "xmin": 1.0 } } },
    { "name": "root_force", "component": "y", "quantity": "reaction",
      "region": { "box": { "xmax": 0.0 } } }
  ],
  "load_cases": ["tip"]
}
```

A transient analysis of each selected load case: the equations of motion
`M a + C v + K u = A(t) f`, the case's prescribed displacements following
`A(t) g`, integrated with a constant step by the HHT-alpha method
(`docs/formulation.md`, section 7e). Every load of the case - forces,
pressures, body forces, self-weight, the rotation's load, a temperature
change - and every prescribed displacement is scaled by the one amplitude.
`sparlab_solve` runs it after the static (and non-linear) analyses;
`sparlab_topopt` refuses a deck with it enabled.

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `enabled` | bool | `false` | run the transient analysis |
| `time_step` | number | required | the constant step `dt` [s] |
| `end_time` | number | required | the duration [s], a whole number of steps |
| `alpha` | number | `0` | HHT-alpha in `[-1/3, 0]`: `0` is the trapezoidal rule (no numerical dissipation, exact energy balance for a linear model), `-0.05` to `-0.1` damps the modes the step does not resolve (CalculiX's default is `-0.05`) |
| `mass` | string | `consistent` | `consistent` or `lumped` (row sums; HRZ diagonal scaling for the Tet10) |
| `damping` | object | none | Rayleigh damping `C = mass M + stiffness K`: `mass` [1/s] and `stiffness` [s], both `>= 0`. The damping ratio of a mode is `mass / (2 omega) + stiffness omega / 2` |
| `amplitude` | object | step of 1 | `type`: `step` (A = `scale` from `t = 0`), `table` (piecewise linear through `times` [s] and `values`, constant beyond the ends), or `harmonic` (`scale sin(2 pi frequency t + phase)`, `frequency` [Hz], `phase` [rad]); `scale` multiplies every kind |
| `start` | string | `rest` | `rest`: `u = 0` and `v = 0` on the free DOFs. `static`: the static equilibrium under `A(0)`, at rest - a preloaded structure released (with a table that drops to zero). Linear transient only |
| `snapshot_every` | integer | `0` | keep every n-th step's displacement and velocity (and the last) for the VTK series; `0` keeps the final state only |
| `monitors` | array | none | quantities recorded at every step. Each has a `name`, a node `region`, a `component` (`x`, `y`, `z`) and a `quantity`: `displacement` (default), `velocity` or `acceleration` - the mean over the region's nodes [m, m/s, m/s^2] - or `reaction`, the sum of the support reactions over them [N] |
| `nonlinear` | bool | `false` | the non-linear transient: the internal forces, loads and tangent of the non-linear system (finite kinematics, the elastic laws, J2 plasticity) in Newton's method at every step. It takes the keys below |
| `kinematics`, `material_model`, `mean_dilatation`, `follower_pressure` | | as `nonlinear` | as in the `nonlinear` block; only with `nonlinear: true` |
| `residual_tolerance`, `displacement_tolerance` | number | `1e-8` | Newton's tolerances: the residual over the largest force involved (inertia, loads, damping, reactions), the correction over the larger of the step's increment and the displacement |
| `max_iterations` | integer | `25` | Newton iterations per step |
| `max_cuts` | integer | `8` | successive halvings of a step whose Newton iteration fails, from the last converged state |
| `load_cases` | array of strings | all | load cases to integrate, by name |

A step load applied at `t = 0` accelerates the structure at once (the run
starts in equilibrium, the initial acceleration solving `M a = A(0) f - K u`),
which suits a linear model. A load applied suddenly at a node of a
non-linear model can crush the element under it onto a spurious branch of the
Saint Venant-Kirchhoff law; ramp it up with a `table`, and keep a harmonic
load's period resolved by at least 20 steps (the trapezoidal rule lengthens a
period by `(omega dt)^2 / 12` of itself; the run warns below 20). Damping
acts through the linear elastic stiffness in a non-linear run. Refused with a
message: an `end_time` that is not a whole number of steps, `alpha` outside
`[-1/3, 0]`, negative damping, a non-linear key without `nonlinear: true`, and
`start: static` with `nonlinear: true`.

The results, per load case:

| File | Content |
|------|---------|
| `transient_<lc>.csv` | one row per step (step 0 the initial state): time, the largest nodal displacement, the kinetic, strain (with plasticity: stored, elastic plus hardening) and damping energies, the external work, the energy balance `E_0 + W - T - U - D`, with a non-linear run the Newton iterations and halvings, with plasticity the points that yielded and the largest plastic strain, and the monitors |
| `transient_state_<lc>.csv`, `transient_reactions_<lc>.csv` | the final displacement, velocity and acceleration per node, and the final support reactions (with `output.csv`) |
| `transient_<lc>_<k>.vtk`, `transient_<lc>.vtk.series` | the snapshots (displacement and velocity with their magnitudes), numbered, and ParaView's file-series index of their times (with `output.vtk`) |

The `transient` block of `summary.json` holds the method (with `beta`,
`gamma`), the options, and per load case whether it completed, the steps,
the model's mass, the largest displacement and when it occurred, the
energies at the end, the largest relative energy balance over the path
(round-off for the trapezoidal rule on a linear model; the numerical
dissipation with `alpha < 0`), the final balance in joules (with plasticity
the plastic dissipation), every monitor's maximum, minimum, their times, its
final value and its nodes, with a non-linear run the Newton iterations and
halvings and the plastic strain, and warnings.

## `frequency_response`

```json
"frequency_response": {
  "enabled": true,
  "frequencies": { "start": 10.0, "end": 3000.0, "count": 120, "spacing": "log" },
  "mass": "consistent",
  "damping": { "structural": 0.02, "mass": 2.0, "stiffness": 1e-6 },
  "snapshot_frequencies": [80.0, 500.0],
  "monitors": [
    { "name": "tip_uy", "component": "y", "region": { "box": { "xmin": 1.0 } } },
    { "name": "tip_ay", "component": "y", "quantity": "acceleration",
      "region": { "box": { "xmin": 1.0 } } }
  ]
}
```

The steady harmonic response of each selected load case: its loads as the
amplitudes of `f cos(omega t)` and its prescribed displacements as those of
`g cos(omega t)` (a shaken support), solved directly at every frequency,
`[K (1 + i eta) - omega^2 M + i omega (a M + b K)] U = f` (section 7e).

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `enabled` | bool | `false` | run the harmonic response |
| `frequencies` | array or object | required | the frequencies [Hz], `>= 0`: a list, or `{ "start", "end", "count", "spacing" }` with `spacing` `linear` (default) or `log` (then `start > 0`) |
| `mass` | string | `consistent` | `consistent` or `lumped` |
| `damping` | object | none | `structural`: the loss factor `eta` of `K (1 + i eta)`; `mass`, `stiffness`: Rayleigh's `a`, `b`. All `>= 0` |
| `snapshot_frequencies` | array | none | frequencies [Hz] whose complex field is written (each the nearest solved one) |
| `monitors` | array | none | as in `transient`; recorded as complex amplitudes: a velocity is `i omega U`, an acceleration `-omega^2 U` |
| `load_cases` | array of strings | all | load cases to analyse, by name |

Without damping the dynamic stiffness is singular at a natural frequency: a
frequency where the LU fails, or whose solve has a backward error above
1e-10, stops the run with a message; a response more than a million times
the static one is flagged as round-off. Give the model damping, or move the
frequency.

| File | Content |
|------|---------|
| `frequency_response_<lc>.csv` | one row per frequency: the largest nodal displacement over a cycle and, per monitor, its real and imaginary parts, modulus and phase `arg U` in degrees (against the load: 0 in phase, 180 in antiphase) |
| `frequency_response_field_<lc>_<k>.csv` | the complex nodal amplitudes at each snapshot frequency (in the first column's header) and each node's largest displacement over a cycle (with `output.csv`) |
| `frequency_response_<lc>_<k>.vtk`, `frequency_response_<lc>.vtk.series` | the real and imaginary parts of the displacement and its peak over a cycle, with a file-series index whose "time" is the frequency (with `output.vtk`) |

The `frequency_response` block of `summary.json` holds the method, the
frequencies, the damping and per load case the largest displacement and its
frequency, every monitor's peak amplitude, the frequency and phase there and
its nodes, the snapshot frequencies and warnings.

## `topology`

```json
"topology": {
  "enabled": true,
  "volume_fraction": 0.35,
  "initial_density": -1,
  "simp": { "penalty": 3.0, "emin_ratio": 1e-9, "mass_floor": 1e-9,
            "mass_interpolation": "penalty_matched" },
  "filter": { "type": "density", "radius_elements": 3.0 },
  "optimizer": { "max_iterations": 600, "change_tolerance": 1e-2,
                 "objective_tolerance": 5e-5, "objective_window": 20,
                 "move_limit": 0.2, "damping": 0.5,
                 "continuation_steps": 3, "penalty_start": 1.5,
                 "continuation_iterations": 30,
                 "volume_tolerance": 1e-10, "max_bisections": 200,
                 "history_stride": 1, "interpretation_threshold": 0.5 },
  "passive_regions": [
    { "name": "bolt_holes", "type": "void",
      "region": { "circle": { "center": [0.030, 0.045], "radius": 0.010 } } }
  ]
}
```

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `enabled` | bool | `false` | `sparlab_topopt` requires `true` |
| `volume_fraction` | number | `0.4` | target `nu`, in `(0, 1]`; measured on the physical density over the **whole** domain, passive regions included |
| `initial_density` | number | `-1` | starting value for free variables; negative means "use `volume_fraction`". The value used is recorded as `optimization_setup.initial_density`. Start below the target when passive solid regions take a large share of the budget, or the first update may not reach it (`engine_mount_3d`) |

**`simp`**

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `penalty` | number | `3.0` | `p >= 1` |
| `emin_ratio` | number | `1e-9` | `E_min/E_0`, in `(0, 1)` |
| `mass_floor` | number | `1e-9` | mass floor for `penalty_matched` |
| `mass_interpolation` | string | `penalty_matched` | or `linear` |

**`filter`**

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `type` | string | `density` | `density`, `sensitivity` or `none` |
| `radius` | number [m] | - | radius in metres; takes precedence when `> 0` |
| `radius_elements` | number | `1.5` | radius as a multiple of the mean element size |

Specify the radius in **metres** when comparing meshes: that holds the minimum
length scale fixed, which is what makes the optimum mesh convergent.

**`optimizer`**

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `max_iterations` | integer | `200` | iteration cap |
| `change_tolerance` | number | `1e-2` | stop when `max |dx|` falls below this |
| `objective_tolerance` | number | `0` | stop when the relative spread `(max - min) / |c_k|` of the last `objective_window + 1` compliances falls below this; `0` disables. For a monotone history that is the change over the window, and an oscillating design cannot satisfy it mid-cycle |
| `objective_window` | integer | `5` | window for the above |
| `move_limit` | number | `0.2` | maximum per-iteration density change, in `(0, 1]` |
| `damping` | number | `0.5` | OC exponent `eta`, in `(0, 1]` |
| `continuation_steps` | integer | `1` | `> 1` enables penalty continuation |
| `penalty_start` | number | `1.0` | starting penalty, `>= 1` and `<= simp.penalty` |
| `continuation_iterations` | integer | `25` | iterations per continuation stage |
| `volume_tolerance` | number | `1e-10` | relative volume error of the multiplier bisection (OC) |
| `max_bisections` | integer | `200` | bisection cap (OC) |
| `history_stride` | integer | `1` | record a density snapshot every N iterations; `0` disables |
| `interpretation_threshold` | number | `0.5` | density threshold used when reading the field as geometry, in `(0, 1)` |
| `method` | string | `oc` | `oc` (optimality criteria, volume constraint only) or `mma` (method of moving asymptotes, any number of constraints) |
| `constraint_tolerance` | number | `1e-4` | MMA: largest constraint value `g_i <= tol` for an iterate to count as feasible; convergence requires feasibility |
| `mma.move_limit` | number | `move_limit` | MMA per-iteration bound on the design change |
| `mma.asymptote_init` | number | `0.5` | initial asymptote distance as a fraction of the variable range |
| `mma.asymptote_increase` | number | `1.2` | asymptote widening when a variable keeps moving the same way |
| `mma.asymptote_decrease` | number | `0.7` | asymptote tightening when it oscillates |
| `mma.constraint_penalty` | number | `1000` | Svanberg's `c_i`: the price of the artificial variables that relax a constraint |
| `mma.subproblem_tolerance` | number | `1e-7` | final barrier parameter of the interior-point subproblem solver; each barrier level is converged when the KKT residual falls below 0.9 times it |
| `mma.max_newton_iterations` | integer | `500` | Newton iterations allowed per barrier level, `>= 10`; running out is a `ConvergenceError` naming the residual block that stalled |

```json
"topology": {
  "optimizer": { "method": "mma", "move_limit": 0.1, "constraint_tolerance": 1e-4,
                 "mma": { "asymptote_init": 0.5, "asymptote_increase": 1.2,
                          "asymptote_decrease": 0.7 } },
  "stress": { "enabled": true, "limit": 9.4e6, "p_norm": 8.0,
              "relaxation": 0.5, "scaling_blend": 0.5, "feasibility_tolerance": 1e-3 }
}
```

**`stress`** - an aggregated von Mises stress constraint, one per load case.
Requires `method: "mma"` and the density filter (the sensitivity filter has no
gradient the adjoint could use).

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `enabled` | bool | `false` | add the constraint |
| `limit` | number [Pa] | required | allowable von Mises stress |
| `p_norm` | number | `8` | exponent `P` of the p-norm aggregate; larger tracks the maximum more closely at the price of a rougher constraint |
| `relaxation` | number | `0.5` | qp-relaxation exponent `q`: the relaxed stress is `rho^q sigma_vm(solid)`, so a void element cannot violate the limit |
| `scaling_blend` | number | `0.5` | weight `alpha` in the adaptive scale `c_k = alpha s_max/g_PN + (1 - alpha) c_{k-1}` that makes the aggregate track the true maximum |
| `feasibility_tolerance` | number | `1e-3` | relative margin the summary uses when it reports whether the *relaxed* maximum meets the limit |

`docs/topology_optimization.md` gives the formulation and says what the
constraint does and does not guarantee.

```json
"topology": {
  "projection": { "enabled": true, "eta": 0.5, "beta_start": 1, "beta_max": 32,
                  "beta_factor": 2, "beta_interval": 40,
                  "advance_on_convergence": true }
}
```

**`projection`** - the smoothed Heaviside projection of the filtered density,
with continuation of its sharpness `beta`. The projected density is the
physical one: SIMP, the volume constraint, the stress constraint and the
grey level all act on it. Needs the density filter (or no filter); with the
sensitivity filter it is a `ConfigError`.

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `enabled` | bool | `false` | project the filtered density |
| `eta` | number | `0.5` | threshold of the step, in `(0, 1)` |
| `beta_start` | number | `1` | first sharpness, `> 0` |
| `beta_max` | number | `32` | final sharpness, `>= beta_start` |
| `beta_factor` | number | `2` | multiplier per continuation stage, `> 1` |
| `beta_interval` | integer | `50` | most iterations spent at one `beta` |
| `advance_on_convergence` | bool | `true` | move to the next `beta` as soon as the design change at the current one falls below `change_tolerance` |
| `robust` | bool | `false` | robust formulation: objective and constraints on the eroded design (`eta + robust_delta`), volume on the dilated one (`eta - robust_delta`), the blueprint (`eta`) reported and exported (`docs/topology_optimization.md`, section 3c) |
| `robust_delta` | number | `0.1` | `delta` of the eroded and dilated thresholds; `eta +- delta` must lie in `(0, 1)` |
| `robust_volume_interval` | integer | `1` | iterations between rescalings of the dilated volume target; use about 20 with OC |
| `erosion_check` | bool | `false` | without `robust`: evaluate the final design at `eta +- robust_delta` once, to report how much it loses if the part comes out thinner or thicker |

Convergence is only declared at `beta_max`, and the objective-stall history
restarts at every `beta` step because the objective changes with it. A run
that stops at the iteration cap before reaching `beta_max` says so in a
warning.

```json
"topology": {
  "optimizer": { "method": "mma", "move_limit": 0.05, "mma": { "max_newton_iterations": 3000 } },
  "buckling_constraint": { "enabled": true, "min_load_factor": 6.0, "num_modes": 6,
                           "ks_parameter": 40.0, "solid_threshold": 0.5,
                           "load_cases": ["axial"] }
}
```

**`buckling_constraint`** - the lowest positive buckling load factors of the
SIMP design at or above a required multiple of each constrained load case,
aggregated by a KS function (`docs/topology_optimization.md`, section 5d).
Requires `method: "mma"` and the density filter (or none).

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `enabled` | bool | `false` | add the constraint |
| `min_load_factor` | number | `1` | required load factor `lambda_req`, `> 0` |
| `num_modes` | integer | `6` | load factors aggregated per load case, in `[1, 50]` |
| `ks_parameter` | number | `40` | KS parameter `P`, in `[1, 500]`; the aggregate overestimates the largest ratio by at most `ln(num_modes)/P` |
| `solid_threshold` | number | `0.5` | density at which an element counts as solid in the pseudo-mode diagnostic |
| `load_cases` | array of strings | all | load cases to constrain, by name |
| `tolerance`, `residual_tolerance`, `max_iterations`, `seed` | | as `buckling` | settings of the eigensolve |

The load factors it holds are the SIMP model's; the `buckling` section
checks the exported part, which is the number to judge.

```json
"topology": {
  "overhang": { "filter": true, "build_direction": "+y",
                "smax_exponent": 40, "smax_reference": 0.5, "smin_epsilon": 1e-4 },
  "length_scale_check": { "enabled": true, "max_radius_elements": 8, "tolerance": 0.02 }
}
```

**`overhang`** - the additive-manufacturing overhang rule
(`docs/topology_optimization.md`, section 3d). Its presence switches on the
overhang check of the final design; `filter` also applies the filter during
the optimisation. Needs a `structured_quad` or `structured_hex` mesh and,
with the filter, the density filter (or none).

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `build_direction` | string | required | `+x`, `-x`, `+y`, `-y` (and `+z`, `-z` in 3-D): the axis the part grows along, from the plate at the low end (`+`) or the high end (`-`) |
| `filter` | bool | `false` | apply the overhang filter between the density filter and the projection |
| `smax_exponent` | number | `40` | `P` of the smooth maximum, in `[2, 200]`, and large enough that `Q = P + ln(n)/ln(xi_0)` stays positive for the supports present |
| `smax_reference` | number | `0.5` | `xi_0`, the value at which the smooth maximum of equal supports is exact, in `(0, 1)` |
| `smin_epsilon` | number | `1e-4` | `epsilon` of the smooth minimum, in `(0, 0.1)` |

**`length_scale_check`** - the morphological measurement of the smallest
member and gap of the final thresholded design. On by default in a robust
run; the section's presence switches it on otherwise.

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `enabled` | bool | on with `robust` or when the section is present | run the check |
| `max_radius_elements` | number | `6` | largest probe radius, in cells (probes step by half a cell), `>= 0.5` |
| `tolerance` | number | `0.02` | share of the solid (void) volume a probe may change and still pass, in `[0, 0.5)`; opening rounds convex corners, so zero is too strict |

**`passive_regions`** - an array of

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `name` | string | auto | diagnostics |
| `type` | string | `solid` | `solid` (`x = 1`) or `void` (`x = density`) |
| `density` | number | `0` | density imposed on a void region, in `[0, 1)` |
| `region` | object | required | element region |

## `output`

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `csv` | bool | `true` | per-node and per-element CSV tables |
| `vtk` | bool | `true` | legacy VTK field files for ParaView |
| `density_history` | bool | `true` | density snapshots for the animation |
| `mode_shapes` | bool | `true` | mode-shape CSV and per-mode VTK |

`summary.json` and `mesh.json` are always written: a result directory should be
self-describing.

## Command-line overrides

`sparlab_topopt` accepts overrides so a parametric study never has to edit or
duplicate a deck. Each corresponds to one deck field:

```
--volume-fraction        topology.volume_fraction
--penalty                topology.simp.penalty
--filter-radius          topology.filter.radius            (metres)
--filter-radius-elements topology.filter.radius_elements   (cells)
--filter-type            topology.filter.type
--max-iterations         topology.optimizer.max_iterations
--method oc|mma          topology.optimizer.method
--stress-limit <Pa>      topology.stress.limit, and enables the constraint and MMA
--no-stress              topology.stress.enabled = false
--nx --ny (--nz)         mesh.nx, mesh.ny (mesh.nz, solid meshes only; an error on a file mesh)
--youngs-modulus         material.youngs_modulus
--load-weights w1,w2,... one weight per load case, in deck order
--modes N                modal.enabled = true, modal.num_modes = N
--solver TYPE            solver.linear.type
--projection             topology.projection.enabled = true (deck or default schedule)
--no-projection          topology.projection.enabled = false
--beta-max B             topology.projection.beta_max, and enables the projection
--buckling N             buckling.enabled = true, buckling.num_modes = N
--min-load-factor L      topology.buckling_constraint.min_load_factor, and enables the constraint and MMA
--no-buckling-constraint topology.buckling_constraint.enabled = false
--robust / --no-robust   topology.projection.robust (--robust also enables the projection)
--overhang DIR           topology.overhang.build_direction = DIR, filter = true
--no-overhang-filter     topology.overhang.filter = false (the check stays)
--tag NAME               appended to the case name in summary.json
```

`sparlab_solve` takes `--solver`, `--modes` and `--buckling` too, and
`--nonlinear`, which enables the non-linear analysis with the deck's
`nonlinear` settings or their defaults. With `--export-calculix` it
additionally writes one CalculiX input deck per load case
(`calculix_<case>.inp`: CPS4 / CPE4, CPS3 / CPE3, C3D8, C3D4 or C3D10
elements, the same nodes, supports and nodal loads, every field within
CalculiX's 20 characters) into the output directory, which is what
`scripts/run_cross_validation.sh` feeds to `ccx`, as a static step and, for a
run with buckling, as a `*BUCKLE` step. A non-linear Saint Venant-Kirchhoff
case also goes out as `calculix_<case>_nlgeom.inp`: a `*STEP, NLGEOM` with
CalculiX's own load cards (a follower `*DLOAD` pressure, `CENTRIF`, `GRAV`;
a dead pressure as the nodal forces of the reference faces), increments of
at most 1/50 of the load (CalculiX lags a centrifugal load at the deformed
position within an increment) and convergence controls of 1e-6 on residual
and correction. A neo-Hookean case is not exported: CalculiX's `NEO HOOKE` is
a different strain energy; nor is a `finite_logarithmic` case (a non-linear
or transient one): CalculiX's `NLGEOM` elasticity is the Saint
Venant-Kirchhoff law and its finite-strain plasticity multiplicative, neither
the Hencky law nor plasticity additive in the log strain.

Run any app with `--help` for its full flag list.
