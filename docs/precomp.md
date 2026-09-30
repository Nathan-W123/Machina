# precomp: the Python layer for springback pre-compensation

`precomp` is the Python package around SparLab's forming solver for robotic
single-point incremental forming (SPIF). It turns a target part into a tool
path, a forming deck and a simulated formed shape; measures real parts from
scans; and closes the loop by displacement adjustment - with the finite-element
model, a learned surrogate (`precomp.ml`, built on the API below and documented
in [precomp_ml.md](precomp_ml.md)),
or the measured part itself.

It adds no physics of its own. Every number it reports is either read from a
solver result, measured from a scan, or computed by a geometric procedure
whose accuracy is checked against an exact answer in `python/tests` (the
table in [Verification](#verification) gives the measured errors). What it
does *not* do is listed in [Limitations](#limitations).

## Install and test

```bash
pip install -e .                 # the package `precomp` from python/precomp, and the `precomp` command
pip install -e '.[dev]'          # + pytest
pip install -e '.[torch]'        # + torch, for the neural models of precomp.ml only
python3 -m pytest python/tests -q    # 153 tests (42 for precomp.ml), ~2.5 min; the 6 integration
                                     # tests (~3 min of it) skip without build/bin/sparlab_form
```

Python 3.10 or newer; numpy, scipy, pandas, scikit-learn, joblib, contourpy
and matplotlib. `python/sparlab_viz` is not part of the package and is
untouched. The tests import the package from the source tree, so they also
run without installing it.

The forming simulation needs the C++ executable `sparlab_form`. It is looked
up as `FormingSetup.executable`, else `$PRECOMP_SPARLAB_FORM`, else
`build/bin/sparlab_form` relative to the working directory or the
repository root; its absence is an error naming every place looked at.

## Conventions

**Units.** Strict SI, as in the rest of SparLab (`docs/conventions.md`):
m, N, Pa, s, kg/m^3. Angles carry their unit in their name (`wall_angle_deg`;
radians otherwise). Millimetres appear only in human-facing output - figure
axes and report tables - and are always labelled `[mm]`. A scan in
millimetres is read with an explicit `scale=0.001`; a cloud whose size comes
out above 20 m or below 0.1 mm draws a `UnitsWarning` naming `scale`.

**The surface.** The blank lies flat in the x-y plane with its tool side (top)
at z = 0 and the material below, z in [-t, 0]. SPIF pushes the tool down, so
a formed part is a height field z = f(x, y) <= 0 *of the tool-side surface*.
Targets, commanded shapes, simulated and scanned parts are all this one
object, a `HeightMap`. Parts are centred on the origin of the blank, which is
also the centre of the mesh.

**Arrays.** A `HeightMap` holds `z[ny, nx]` on a regular `Grid` of spacing
`h`: row j is y = y0 + j h, column i is x = x0 + i h, x fastest (C order) -
the node order of SparLab's structured grids, so row 0 is the *lowest* y.
`z` is finite everywhere; nodes without data (off a scan, outside a
simulated surface) are marked `False` in `mask` and carry the nearest valid
height as a fill value, so derivatives stay defined.

**Signs.** Normals point +z, towards the tool. A signed deviation is positive
when the formed surface lies on the tool side of the target (the part came
out too shallow). The wall angle is arctan |grad z|, 0 on a flat region. The
tool force columns of the solver are the force the sheet exerts on the tool
(fz > 0 while it pushes down; `docs/forming.md`, section 3), which is the
load the robot carries.

## Architecture

Each box is a module; an arrow means "imports". There are no cycles, and
nothing below `api` knows about machine learning.

```
                +--------------------------------------------+
                |  cli   the precomp command (argparse)      |
                |  api   predict / compensate / compare_scan |---> precomp.ml (lazily:
                +------+-------------+-------------+---------+      cli, registry)
                       |             |             |
          +------------v--+   +------v-----+   +---v--------+
          | compensation  |   | report     |   | robot      |
          | DA, predictors|   | figures,md |   | compliance |
          +---+--------+--+   +------+-----+   +---+--------+
              |        |             |             |
     +--------v--+  +--v-------------v--+          |
     | fea       |  | metrology         |          |
     | setup,deck|  | clouds, ICP,      |          |
     | runner,   |  | deviation, metrics|          |
     | results   |  +---------+---------+          |
     +--+-----+--+            |                    |
        |     |               |    +---------------v------------------+
        |     +---------------|--->| toolpath  drop cutter, paths, CSV |
        |                     |    +---------------+------------------+
        |                     |                    |
  +-----v------+      +-------v--------------------v----------------+
  | materials  |      | geometry  HeightMap, Grid, STL, part families |
  +------------+      +----------------------------------------------+
```

| Module | Contents |
|--------|----------|
| `precomp.geometry` | `Grid`, `HeightMap` (interpolation, gradients, normals, wall angle, curvature, smoothing, resampling, points, STL, `.npz`), `read_stl` / `write_stl` / `raycast_top`, the part families |
| `precomp.materials` | `Material` (elasticity, linear + Voce isotropic and Prager kinematic hardening; optional Armstrong-Frederick backstress and Hill48 r-values), Swift / Hollomon conversion, the nominal alloy library, `to_sparlab()` |
| `precomp.toolpath` | `tool_center_surface` (drop cutter), `contour_toolpath`, `spiral_toolpath`, `Toolpath` (trajectory and robot CSV, summary), `dsif_support_points` (the DSIF support tool, synchronised), `rim_pass_path`, `lift_heights` / `lift_cutter` / `reach_from_below` (a ball under the sheet), `signed_outline_distance` |
| `precomp.fea` | `FormingSetup`, `build_deck`, `deck_hash`, `run_deck`, `simulate` (content-addressed cache), `simulate_many` (process pool), `load_result` / `FormingResult`; `precomp.fea.support` - backing plate, DSIF support and rim pass in the deck, `command_upper_bound` ([Rim support](#rim-support)) |
| `precomp.metrology` | `read_point_cloud`, `align` (robust point-to-plane ICP), `signed_deviation`, `vertical_deviation`, `metrics`, region masks |
| `precomp.compensation` | `displacement_adjustment`, `update_from_scan`, `limit_wall_angle`, `FEAPredictor`, `SurrogatePredictor`, `CompositePredictor`, the `FieldModel` protocol |
| `precomp.robot` | `CartesianCompliance` (deflection, pre-compensation, from joint stiffness), `forces_on_path`, `compensate_toolpath` |
| `precomp.api` | `predict`, `compensate` -> `CompensationResult`, `compare_scan` -> `ScanComparison` |
| `precomp.report` | deviation and height maps, histograms, sections, convergence plots, `write_report` (Markdown) |
| `precomp.cli` | the `precomp` command |
| `precomp.ml` | learned deviation models, conformal intervals, the training envelope, surrogate compensation - [precomp_ml.md](precomp_ml.md) |

## Geometry

### HeightMap

| Method | Result | Notes |
|--------|--------|-------|
| `interpolate(x, y, masked=True)` | z [m] | bilinear; NaN off the grid, and (with `masked`) where the stencil touches an invalid node |
| `gradient(sigma=None)` | (dz/dx, dz/dy) | central differences, or Gaussian derivatives at length scale `sigma` [m] |
| `normals(sigma)`, `wall_angle(sigma)` | (ny, nx, 3), rad | unit normal (-z_x, -z_y, 1)/norm |
| `curvature(sigma)` | (H [1/m], K [1/m^2]) | of the graph surface; H > 0 where the surface bends towards +z |
| `smooth(sigma)`, `resample(grid)` | HeightMap | resampled nodes off the source are invalid |
| `to_points()`, `from_points(points, grid, max_gap=None)` | (n, 3), HeightMap | scattered -> grid by linear Delaunay interpolation; outside the hull, or farther than `max_gap` from data, invalid |
| `to_stl(path)`, `from_stl(path, grid, scale)` | | two triangles per valid cell; reading casts a ray along -z per node and keeps the top-most hit |
| `save(path)`, `load(path)` | `.npz` | `z`, `mask`, `grid` (JSON), `metadata` (JSON), `format = "precomp.HeightMap/1"` |

The Gaussian derivatives use kernels normalised by their moments (the
smoothing kernel sums to 1, the second-derivative kernel to 0), so they are
exact for any quadratic whatever its offset from z = 0. `scipy.ndimage`'s own
derivative filters are not: on a surface 15 mm below the sheet plane their
second-derivative kernel leaked the height into the curvature (a saddle's
Gaussian curvature came out +15 instead of -7.1 1/m^2), which the tests now
guard against.

### Part families

Each family is a frozen dataclass with physical parameters, `height(x, y)`,
`heightmap(grid)`, `bounds`, `sample(rng)`, `max_wall_angle_deg()`,
`max_depth()`, `footprint_radius()` and `to_dict()` / `part_from_dict()` (a
`family` key). Every constructor refuses a wall steeper than
`MAX_WALL_ANGLE_DEG = 65` and inconsistent dimensions with a message naming
the parameters; `sample` redraws until the part is valid.

| Family | Parameters | Construction |
|--------|------------|--------------|
| `truncated_cone` | `top_radius`, `wall_angle_deg`, `depth`, `top_fillet`, `bottom_fillet` | profile of the radial distance |
| `pyramid` | `half_width_x`, `half_width_y`, `wall_angle_deg`, `depth`, `corner_radius`, fillets | profile of the distance to a rounded rectangle |
| `dome` | `opening_radius`, `depth`, `rim_fillet`, `aspect` | spherical cap tangent to a rim fillet (closed form); `aspect` != 1 scales y (ellipsoidal cap) |
| `elliptic_cone` | `semi_axis_x`, `semi_axis_y`, `wall_angle_deg`, `depth`, fillets | profile of the exact distance to an ellipse |
| `two_level` | `top_radius`, `wall_angle_1_deg`, `depth_1`, `ledge_width`, `wall_angle_2_deg`, `depth_2`, `fillet_radius` | stepped cone, four filleted corners |
| `freeform` | `dents` ((cx, cy, amplitude, sigma), ...), `inner_radius`, `outer_radius`, `seed` | sum of Gaussian dents under a C2 radial window; `sample` scales the amplitudes to at most 55 deg and 35 mm deep |
| `saddle` | half-widths, `corner_radius`, `wall_angle_deg`, `depth`, `saddle_amplitude`, fillets | rounded-rectangle tray whose floor is z = -(depth + a ((y/hy)^2 - (x/hx)^2)) |

The cone-like families are a C1 profile g(w) of an inward distance w from the
nominal opening: a polyline whose corners are replaced by circular fillets
tangent to both segments (`FilletedProfile`). The distance functions are
exact and of unit gradient wherever the profile slopes, so every straight
wall has exactly the design angle. `top_radius` (and the half-widths and
semi-axes) are the *nominal* opening where the wall, extended, meets the
sheet plane; the rim fillet widens the real opening by its tangent length.
The flange stays at z = 0 everywhere outside `footprint_radius()`.

## Materials

`Material.to_sparlab()` writes SparLab's `material` block
(`docs/configuration.md`, `material.plasticity`):

```json
{"name": "AA5754-O", "youngs_modulus": 7e10, "poisson_ratio": 0.33, "density": 2670,
 "plasticity": {"yield_stress": 1e8, "hardening_modulus": 2.38e8,
                "saturation_stress": 1.24e8, "saturation_rate": 14.4,
                "kinematic_hardening_modulus": 0,
                "yield_criterion": "hill48",
                "anisotropy": {"r0": 0.75, "r45": 0.7, "r90": 0.8,
                               "rolling_direction": [1, 0, 0], "sheet_normal": [0, 0, 1]}}}
```

A material with r-values (`r0`, `r45`, `r90`, all three) yields by Hill's
1948 criterion calibrated by them, with the rolling direction along the
deck's x axis and the sheet normal along z; `yield_stress` is then the
uniaxial yield stress along the rolling direction. A material with an
Armstrong-Frederick term (`af_C` > 0 [Pa], `af_gamma` >= 0 [-]: d beta =
2/3 C d eps_p - gamma beta d a, added to Prager's) writes it as one Chaboche
backstress, `"backstresses": [{"modulus": af_C, "recovery": af_gamma}]`.
Without them the block is von Mises with Prager's rule only. These are the
keys SparLab reads; every deck is run with `--strict-config`.

Power laws are converted, never passed through: `from_swift(K, eps0, n)`
keeps the initial yield stress K eps0^n exact and fits H, Q and delta to
sigma = K (eps0 + a)^n over plastic strains 0-0.6 (non-negative least squares
per delta, delta by a bounded scalar search); `from_hollomon(K, n,
yield_stress)` is Swift through the yield point. The fit is recorded in
`hardening_source`. The library holds nominal handbook-order values, marked
so in `source` - **not certified data**:

| Material | E [GPa] | yield [MPa] | Hollomon K [MPa], n | r0 / r45 / r90 | fit RMS / max error |
|----------|--------:|------------:|--------------------|----------------|--------------------:|
| AA5754-O | 70.0 | 100 | 420, 0.30 | 0.75 / 0.70 / 0.80 | 1.6 % / 6.2 % |
| AA2024-O | 73.1 | 75 | 330, 0.20 | 0.70 / 0.80 / 0.70 | 3.3 % / 18.6 % |
| AA6061-T6 | 68.9 | 276 | 410, 0.06 | 0.60 / 0.70 / 0.65 | 0.9 % / 4.3 % |
| AA7075-O | 71.7 | 103 | 400, 0.17 | 0.70 / 0.90 / 0.80 | 3.2 % / 20.2 % |
| DC04 | 210 | 170 | 530, 0.22 | 1.90 / 1.50 / 2.20 | 1.6 % / 6.5 % |
| Ti-6Al-4V annealed | 113.8 | 880 | 1100, 0.04 | - | 0.4 % / 1.7 % |

The largest errors sit in the first percent of plastic strain, where a
Hollomon curve with a small eps0 rises almost vertically and one Voce term
cannot follow it. Kinematic hardening is zero throughout: its parameters come
from reverse-loading tests that handbooks do not give.

## Tool paths

`tool_center_surface(hm, R)` is the drop cutter of a ball of radius R: the
tool-centre height c(x) = max over |d| <= R of [z(x + d) + sqrt(R^2 - |d|^2)],
a grey-scale dilation with a non-flat spherical structuring element. It is
exact over the grid nodes - the ball touches at least one node and contains
none. Between nodes the surface is known only by interpolation, which can
rise above the node samples by about h^2 / (8 R cos^3 theta) at a wall of
angle theta; that is the gouge a grid of spacing h can hide (4 um measured
for h = 0.25 mm, R = 5 mm on a 50 deg cone).

`contour_toolpath` puts level k at the constant tool-centre height
R - k step_down (the tip at -k step_down) and follows the level set of the
drop-cutter surface there, so every point touches the target. The last level
lies step_down / 1000 above the deepest point: its loop traces the edge of a
flat floor, which is not swept, as usual in SPIF. Loops come from `contourpy`,
are closed, oriented counter-clockwise (or clockwise, or alternating), start
at the point nearest the previous loop's start, and are resampled at equal
arc length. Several loops on one level (two pockets) are visited by nearest
neighbour through a retract; a move between levels descends directly when
nothing is in the way. Points are finally lifted onto the bilinear
drop-cutter surface wherever a marching-squares chord cut beneath it.

`spiral_toolpath` descends continuously: through revolution k the height
falls linearly with arc length from level k-1 to level k, and the x-y
position is found by bisection on the line from the point of loop k-1
towards its nearest point on loop k, where the drop-cutter surface has
exactly that height. The height never increases and the tool touches the
surface after the first revolution, which starts with the tip on the sheet.
It needs one loop per level.

`Toolpath` holds the points, a `level` per point (`AIR = -1` for approach,
retract and traverse), the arc length `s` and the pseudo-time `t = s / L`
(strictly increasing; duplicate points are dropped). `to_sparlab_csv` writes
`t,x,y,z` with 17 significant digits; `to_robot_csv(path, feed)` writes
`x,y,z,i,j,k,feed` in metres, a unit tool axis and m/s; `summary(feed)` gives
lengths, levels, height range and the time at that feed.

## The forming deck

`build_deck(setup, commanded, out_dir, toolpath=None, target=None)` writes:

| File | Content | In the hash |
|------|---------|:-----------:|
| `deck.json` | the SparLab deck: `mesh`, `material`, `model`, `forming` | yes |
| `toolpath.csv` | tool-centre trajectory `t,x,y,z`, t from 0 to 1 | yes |
| `support_path.csv` | with `support: "dsif"` only: the support tool's trajectory, on the same t knots (and, with the rim pass, t from 2 to 3) | yes |
| `commanded.npz` | the commanded surface (provenance) | no |
| `precomp_deck.json` | the setup, grid, tool-path summary, package version, the support's report (provenance) | no |

`target` is the part the fixture is made for; only a support strategy reads
it ([Rim support](#rim-support)).

The mesh is `structured_hex` (or `structured_tet`) of the blank
[-L/2, L/2]^2 x [-t, 0], `round(blank_size / element_size)` elements per side
(L is `blank_size` rounded to whole elements) and `layers` through the
thickness. The `forming` block follows `docs/forming.md` (section 2), and the
deck has nothing else: no top-level `boundary_conditions` (every step lists
its own), no `load_cases` (the analysis applies none) and no `nonlinear`
block (sparlab_form does not read one - kinematics and Newton settings live
in `forming`):

```json
"forming": {
  "kinematics": "finite_logarithmic",
  "tools": [{"name": "tool", "shape": "sphere", "radius": 0.005,
             "surface": {"name": "tool_side", "box": {"zmin": 0.0, "xmin": -0.08, "xmax": 0.08,
                                                     "ymin": -0.08, "ymax": 0.08}},
             "friction": 0.1, "trajectory": {"file": "toolpath.csv"}}],
  "steps": [
    {"name": "form",    "type": "form",    "tools": ["tool"], "time": [0.0, 0.998],
     "max_tool_travel": 0.001,             "boundary_conditions": [<clamp>]},
    {"name": "unload",  "type": "release", "tools": [],       "boundary_conditions": [<clamp>]},
    {"name": "release", "type": "release", "tools": [],       "boundary_conditions": [<3-2-1>]}
  ],
  "output": {"vtk": true, "snapshots": "steps"}
}
```

The clamp holds x, y, z of every node within `clamp_margin` of the blank
edges. "form" runs the trajectory from t = 0 to its last point in contact
with the part (the final retract is not simulated); "unload" removes the
tool and ramps its force out with the clamp still on - the springback in the
fixture (sparlab_form warns that the clamp holds reactions: it is not
statically determinate, which is the point); with `release: "321"`,
"release" replaces the clamp by a statically determinate support on three
top nodes of the former clamp - x, y, z at (-a, -a), y, z at (a, -a), z at
(-a, a) - six constraints that remove the rigid-body modes and restrain no
springback. With `"clamped_only"` there is no third step. Every step
constraint is `"mode": "hold"`: it keeps its DOFs where the step finds them,
so the support nodes stay where the clamp held them and the released part
keeps its place on the fixture.

`kinematics` is `"finite_logarithmic"` by default - the large-strain
formulation, whose plastic return works in the logarithmic strain, since
SPIF reaches plastic strains of order one - or `"finite"` / `"small_strain"`.
`contact` sets the tool's `penalty` (the scale s of kappa = s E / h, default
10) and `tangential_penalty` (default 1); `solver` sets the keys of
`forming.newton` (`max_iterations`, `residual_tolerance`,
`displacement_tolerance`, `line_search`, `max_cuts`, `max_increments`) and
`friction_tangent`, `solver`, `mean_dilatation` of `forming`. Any other key
is refused when the setup is made. `build_deck` refuses a tool path whose
tool would reach the clamped frame.

### Runs and the cache

`simulate(setup, commanded, work_dir, target=None)` builds the deck in
`work_dir/staging/<uuid>`, hashes it - SHA-256 of the canonical deck JSON,
the trajectory files its tools read (`toolpath.csv`, then a support's) and
the output of `sparlab_form --version` - and looks the hash up in
`work_dir/runs/<h[:2]>/<h>/`:

```
runs/ab/ab12.../   deck.json  toolpath.csv  commanded.npz  precomp_deck.json
                   output/          the solver's result directory
                   sparlab_form.log stdout and stderr of the run
                   run.json         command, exit code, runtime [s]
                   COMPLETE | FAILED  written last
```

A `COMPLETE` entry is returned without running: an identical deck (same
physics, same trajectory, same solver build) is never run twice. Execution
settings - executable path, timeout, threads - are not part of the hash. A
run executes in staging and is published by one atomic rename, so a reader
never sees a half-written entry and two processes racing on one deck cannot
corrupt it. A failure is recorded in `FAILED` (reason, exit code, log tail)
and raised as `FormingError`; later calls raise the recorded failure again
without running unless `retry_failed=True`. Each run is a subprocess
`sparlab_form --config deck.json --output output --strict-config` with
`OMP_NUM_THREADS = setup.threads` (1 by default): a deck key the solver does
not read fails the run (exit 2) rather than taking a default silently. Any
exit status but 0 is a `FormingError` whose record names it (2 configuration
error, 3 a step stopped - with the reason from `summary.json` - and 4 I/O).

`simulate_many(jobs, work_dir, max_workers, executor="process")` runs
(setup, commanded) pairs - or (setup, commanded, target) triples - on a
process pool (deck building and tool paths run in parallel too) and returns
one `SimulationOutcome` per job, in order - a failed job carries its reason
and is never dropped.

### Result files

The loader reads the result directory as `docs/forming.md` (section 3)
defines it:

| File | Columns / content | Required |
|------|-------------------|:--------:|
| `summary.json` | `completed`, `termination`, `steps` (per step `name`, `type`, `completed`, `files_stem`, ...), `tools`, `timing`, `warnings`, ... | yes |
| `<files_stem>_nodes.csv` | `node,X,Y,Z,ux,uy,uz` - reference coordinates and displacement [m] | yes, per completed step |
| `<files_stem>_elements.csv` | `element,eq_plastic_strain,von_mises_Pa` | no |
| `tool_forces.csv` | `step,increment,t,tool,cx,cy,cz,fx,fy,fz,active_nodes,max_penetration_m` | no (needed for forces) |
| `mesh.json`, `config.json`, `<files_stem>.vtk` | the mesh, the deck, the fields for ParaView | no |

`files_stem` is `step_<k>_<name>`, k counting from 1; the loader takes it
from the summary, so snapshot files are never mistaken for steps, and the
steps it loads are the completed ones. Every table is validated - columns
present, no missing values, unique node ids, the same nodes and reference
coordinates in every step, a node table for every completed step - and a
violation is an error naming the file. `FormingResult.formed_surface(step,
grid)` takes the nodes with reference Z = 0 (the tool side), moves them to
their deformed positions, triangulates them as the reference grid and
interpolates linearly at the grid nodes; a deformed surface that folds over
in plan is refused. `thickness_map` is the distance between the deformed top
and bottom node of each through-thickness column. `forming_forces()` is the
force table with |f| added: the force the sheet exerts **on the tool**
(fz > 0 while it pushes down), one row per converged increment of a step
the tool is active in (`step` from 1, `active_nodes` 0 where it is in the
air).

## Metrology

`read_point_cloud(path, scale)` reads `.xyz` / `.txt` / `.csv` / `.asc`
(whitespace, comma or semicolon; an optional header naming x, y, z),
ASCII and binary `.ply` vertices, `.stl` corners and `.npy` arrays.

`align(points, target, mode)` registers a scan by point-to-plane ICP. The
residual of a point is its distance along the target normal to the tangent
plane at the target point vertically below it - zero exactly on the surface.
Each iteration re-weights robustly (Tukey's biweight at 4.685 MAD-sigmas by
default, Huber, or none; optional trimming) and solves for a small rotation
about the weighted centroid and a translation. `mode` is `rigid`,
`translation_z` or `none`; `fixture_mask` restricts the fit to the target's
flange, the way the forming cell holds the part. Directions the data cannot
fix get no update and `rank` says how many were fixed: 5 for an axisymmetric
part (its rotation about the axis is free), 3 on a flat flange (only z and
the two tilts). The initial guess `auto` refines identity and the four
principal-axes candidates for ten iterations each and keeps identity unless
another candidate halves its median residual.

`signed_deviation(formed, target)` gives, at each target node, the distance
along the target normal to the formed surface: by Newton's method on the
formed height field, or, for a raw cloud, by intersecting the normal with a
plane fitted to the nearest eight points. `metrics(dev, mask, tolerance)`
returns `rms`, `mae`, `max_abs`, `p95_abs`, `bias`, `std`, `n` and
`pct_within_tol` (percent); `part_mask`, `flange_mask`, `wall_mask` (wall
angle above 5 deg) and `region_masks` select the regions.

## Compensation

`displacement_adjustment(target, predictor)` iterates
c_(k+1) = c_k - alpha (f(c_k) - t) (Gan and Wagoner, 2004) over any predictor,
measuring the error vertically or along the target normal (then moving the
commanded surface along the normals and resampling). Each update may be
smoothed (Gaussian, normalised over the adjusted region so the flange does not
leak in), keeps the hold region - the flange by default - at the target
height, keeps the surface at z <= 0 (or below `upper_bound`, what a support
pushing from below can realise - [Rim support](#rim-support)), and is
projected onto the formable set:
`limit_wall_angle` returns the smallest surface above the commanded one whose
central-difference wall angle nowhere exceeds 65 deg (iterated dilation with a
cone; it only raises, never cuts deeper). The result keeps the best predicted
iterate (`commanded`, `formed`, `best_iteration`), the unverified next update
(`proposed`) and a per-iteration history.

For a linear springback f(c) = c + G c with spectral radius rho(G) < 1 and
alpha = 1, the error on the adjusted region obeys e_(k+1) = -G e_k: geometric
convergence at rate rho(G). `update_from_scan(commanded, scan, target)` is
one step of it with the measured part as the predictor (aligned on the flange
by default).

Predictors: `FEAPredictor(setup, work_dir, target=None)` simulates (through
the cache; `target` is the part a support's fixture is made for);
`SurrogatePredictor(model, setup)` and `CompositePredictor(base, model, setup)`
use any object with `predict_deviation(commanded, setup) -> (mean, std)`.

With `sparlab_form` as the predictor (`benchmarks/fea_da_cone`, made by
`python/scripts/fea_da_demo.py`): a 3 mm deep, 45 deg cone on a 40 x 40 x
1 mm AA5754-O blank (20 x 20 x 2 Hex8, 4 mm tool) formed with a clamp as
close as a backing plate came out 0.465 mm RMS (0.970 mm at most) from the
target over the part, and 0.405 mm (0.833 mm) after one DA step; with a
wider clamp 0.739 mm (1.392 mm) and 0.630 mm (1.196 mm). Each simulation
took 2-4 minutes on one thread. DA converges slowly here: the rim sags
where the command cannot rise above the sheet plane, the unswept floor
follows the command by about half, and the response is not local (the
record's README).

## Rim support

A tool pressing from above cannot stop the sheet between the clamp and the
first contour from bending down: on the springback benchmark
(`benchmarks/springback`) every part came out 0.9-1.1 mm too deep near its
rim, 79-87 % of the squared deviation left after any compensation, since a
command held at z <= 0 cannot ask for the rim to come up. `FormingSetup`'s
`support` adds support from below; `precomp.fea.support` writes it into the
deck with the rigid tools `sparlab_form` already has (docs/forming.md), so
nothing in the solver changed.

| `support` | What the deck gets | Steps |
|-----------|--------------------|-------|
| `"none"` (default) | the forming tool only - the deck, and its hash, as before | form, unload, release |
| `"backing_plate"` | `plate`: a `plane` tool (normal +z) through the sheet's bottom z = -t, held still, on the bottom faces wholly outside the opening (the part's outline grown by `clearance`) | form (tool + plate), unload (both removed, clamp on), release |
| `"dsif"` | `support`: a ball of `radius` on the bottom faces of the free window, trajectory `support_path.csv` on the forming tool's pseudo-time | form (tool + support), unload (both removed), release |
| `"dsif"`, `rim_pass: true` | as dsif, and the support's rim pass on pseudo-time 2-3 | form, unload, rim_pass (support alone, clamp on), rim_unload, release |

`support_settings` (defaults in `precomp.fea.setup.SUPPORT_SETTINGS`;
`FormingSetup.resolved_support()` fills them in):

| Strategy | Key | Default | Meaning |
|----------|-----|---------|---------|
| backing_plate | `clearance` | 1 mm | the opening is the outline grown by this |
| both | `friction`, `penalty` | the tool's | contact of plate or support |
| dsif | `radius` | `tool_radius` | the support ball |
| dsif | `squeeze` | 0 | the balls' gap along the normal is (1 - squeeze) t_n; > 0 squeezes, < 0 leaves a gap |
| dsif | `thickness_law` | `"sine"` | t_n = t cos(wall angle) (the sine law), or `"initial"`: t |
| dsif | `clearance` | 2 mm | how far below the part the support waits in the air |
| dsif | `rim_pass` | false | the support-only pass along the rim after forming |
| dsif | `rim_inside`, `rim_outside` | 2 mm, 1 mm | the band the rim pass sweeps, from the outline |
| dsif | `rim_spacing` | 1 mm | between its loops |
| dsif | `rim_max_raise` | 2 mm | cap on a command above the sheet plane in that band |

**The plate's opening.** `sparlab_form` selects a tool's slave surface as
the boundary faces whose nodes all lie in a region; its region primitives
(boxes, circles, annuli, spheres, unions, complements) describe no arbitrary
outline, but its explicit `node_ids` primitive does, with the structured
mesh's documented numbering (node (i, j, k) is k (n+1)^2 + j (n+1) + i,
include/sparlab/mesh/StructuredMesh.hpp; the Tet4 mesh keeps those nodes).
`plate_nodes` lists the corners of every bottom cell whose four corners lie
farther than `clearance` from the part - so the plate's faces are exactly
those cells. On a mesh the opening is therefore up to one element wider
than asked; the provenance (`precomp_deck.json`, `support.plate`) records
the clearance realised (1.35 mm for 1 mm asked on the 2 mm mesh of the
validation cone). The plane holds the sheet up to the sheet plane outside
the opening and nothing holds it down onto the plate; it is removed with
the tool in "unload", before the 3-2-1 release.

**The fixture is made for the part.** The plate's opening and the rim pass
band follow the outline of `target` - the same for every commanded iterate
of that part, as a plate cut for it would be. `build_deck`, `simulate`,
`simulate_many` (a third job item), `FEAPredictor` and
`precomp.api.compensate` pass it on; without it the outline is the
commanded surface's, which moves as compensation raises the rim.
`precomp.ml` simulators get it with every job of a setup that has a
support (`precomp.ml.generate.sim_job`).

**The DSIF support path** (`precomp.toolpath.dsif_support_points`): one
support-ball centre per point of the forming path, written with the forming
tool's t, so both move together in one step. In contact it sits opposite
the tool along the drop-cutter normal n, `p - (R1 + t_n + R2) n`, moved
towards the tool by squeeze t_n - but never above the lift cutter (the drop
cutter upside down, `lift_heights`) of the *in-process* sheet: the target
above the tip's current height, flat at that height below it. On the first
revolution the normal-offset position lies inside the still flat sheet; the
bound puts the ball under it instead. In the air the support waits below
the whole part. The deck builder refuses a support whose contact point
would leave the free window.

**The rim pass** (`precomp.toolpath.rim_pass_path`): closed loops at signed
distances `rim_outside, ..., -rim_inside` from the outline (outside in),
each at the lift-cutter height of the commanded underside - the highest the
ball can rise there without entering the commanded sheet. Where the formed
sheet sagged below the command, the ball pushes it back up to it.

**What a command above the sheet plane can mean.** The forming tool presses
from above; its path is made for the command's part below the plane
(`precomp.fea.deck.forming_surface`). A plate at the sheet's bottom holds
the sheet up to the plane outside its opening and pushes nothing higher;
the synchronised DSIF support pushes up only opposite the forming tool,
whose path never rises above the plane. Only the rim pass pushes the sheet
up on its own, so `command_upper_bound(setup, target)` - a function of the
command - is 0 everywhere except, with the rim pass, on the band and flange
strip it sweeps, where the command may rise above the plane by at most
`rim_max_raise` and only as far as the support ball can push the commanded
underside from below: `reach_from_below` (the morphological opening of the
underside z - t by the ball), cut and repeated to a fixed point. Seen from
below, a part's rim is a concave corner and a raised band next to a held
flange is a slot: a ball larger than their radius cannot reach into them,
and with the 4 mm support ball on the validation cone none of the raise DA
asked for survives. `displacement_adjustment(..., upper_bound=)` (and
`update_from_scan`) keeps every command below it (None: z <= 0, as
before); the wall-angle limit is the forming tool's and applies to the part
below the plane. `compensation_masks` adjusts the swept flange strip rather
than holding it at the plane (the rim pass realises it from below).
`build_deck` refuses a command above the bound - no tool could produce it -
and a command above the plane without a target. `precomp.api.compensate`
does all of this for a setup with a support.

**A trained model and the support.** `support` and `support_settings` are
physics fields: a model's setup envelope records them, so a model trained
on one strategy refuses a setup with another (`setup_mismatch`), and a
model trained before the fields existed - on single-point forming only -
refuses any support.

**Validation** (`benchmarks/support_cone`, `python/scripts/support_validation.py`;
SparLab simulations, not experiment): the springback benchmark's smallest
held-out cone (15.4 mm across, 3.09 mm deep, 38.6 deg; 40 x 40 x 1 mm
AA5754-O, 20 x 20 x 2 Hex8, 5 mm clamp, 4 mm tool, spiral, 3-2-1 release).
Deviation of the released part from the target over the part [mm]; rim sag
= mean vertical deviation of the part less than 1 mm deep:

| support | RMS | rim sag | RMS after one FE-DA step | rim sag | runtime (one thread) |
|---------|--:|--:|--:|--:|--:|
| none | 0.767 | -1.151 | 0.658 | -0.992 | 144 s |
| backing plate, 1 mm clearance (1.35 mm on the mesh) | 0.370 | -0.541 | 0.297 | -0.399 | 434 s |
| DSIF, sine-law gap | 0.275 | -0.221 | 0.298 | -0.053 | 251 s |
| DSIF, gap t | 0.357 | -0.494 | 0.318 | -0.163 | 248 s |
| DSIF, sine-law gap, rim pass | 0.338 | -0.074 | 0.267 | -0.119 | 290 s |

Support from below removes most of the sag before any compensation. The
sine-law DSIF's lead over the plate is mostly a squeeze: the model's wall
comes out about 20 % thicker than t cos(theta), and with a gap of t DSIF
lands near the plate. What is left moves to the deeper part (DSIF leaves it
0.24 mm shallow, and one DA step with alpha = 1 made it worse); the rim
pass cuts the sag further but lifts the part. DA with the rim pass asked
for up to 0.31 mm above the plane at the rim; the 4 mm ball reaches none of
it (a 2 mm ball would reach 0.11 mm), so every command stayed at z <= 0.
Every run completed, every tool's penetration stayed below 1.3 um, the
plate and support were pushed down and the tool up, and the release left
no reaction; the record's README has the details and the limitations.

## Robot compliance

`CartesianCompliance(C)` with C the 3 x 3 compliance of the tool-centre point
[m/N] (or one matrix per path point): `deflected_path(p, f) = p + C f`;
`precompensate_path(p, f) = p - C f`, exact when f is the force at the target
positions, and with a `force_model` of the commanded path a fixed-point
iteration p_cmd = p - C f(p_cmd), which converges when ||C df/dp|| < 1.
`from_joint_stiffness(J, K)` gives C = J_v K^-1 J_v^T for one pose;
`forces_on_path(path, forces)` interpolates a simulated force history at the
path's pseudo-time.

## High-level API and command line

```python
from precomp import api
from precomp.fea import FormingSetup
from precomp.geometry import Grid, TruncatedCone

target = TruncatedCone(0.05, 50.0, 0.02, 0.005, 0.005).heightmap(Grid.centered(0.2, 5e-4))
setup = FormingSetup("AA5754-O", thickness=1e-3, tool_radius=5e-3, step_down=5e-4)
result = api.compensate(target, setup, method="fea", iterations=3, work_dir="runs")
result.save("out/cone")            # compensated.npz, predicted.npz, toolpath.csv, compensation.json
scan = api.compare_scan("part.ply", target, scale=1e-3, fixture="flange", tolerance=2e-4)
```

```bash
precomp part --list
precomp part --family truncated_cone --param top_radius=0.05 wall_angle_deg=50 depth=0.02 \
             --size 0.2 --spacing 5e-4 --out cone.npz --stl cone.stl
precomp part --family freeform --sample --seed 7 --out freeform.npz
precomp material                              # the library
precomp material AA5754-O --sparlab           # its SparLab block
precomp setup --material DC04 --set thickness=8e-4 tool_radius=0.006 --out setup.json
precomp toolpath --surface cone.npz --setup setup.json --out toolpath.csv --robot robot.csv --feed 0.03
precomp simulate --setup setup.json --commanded cone.npz --work-dir runs --out formed.npz
precomp compensate --target cone.npz --setup setup.json --iterations 3 --work-dir runs --out comp
precomp scan compare --scan part.ply --scale 0.001 --target cone.npz --fixture flange \
                     --tolerance 2e-4 --out cmp
precomp scan update --scan part.ply --scale 0.001 --commanded comp/compensated.npz \
                    --target cone.npz --out next.npz
precomp report --target cone.npz --formed formed.npz --tolerance 2e-4 --out report
precomp compensate --target cone.npz --setup setup.json --method surrogate --model models/gbm \
                   --verify-fea --work-dir runs --out comp     # see precomp_ml.md
```

Exit status 0 on success, 2 for a usage or input error (the message names
the cause), 3 when a simulation fails; `--debug` shows the traceback.

## Extension points for precomp.ml

The ML layer builds on this API and is imported only when used:

* **Models as predictors.** An object with
  `predict_deviation(commanded, setup) -> (mean, std)` - the vertical
  deviation dz = z_formed - z_commanded [m] on the commanded grid, `std` [m]
  or None - is a `FieldModel`: `SurrogatePredictor`, `CompositePredictor`,
  `api.predict` and `api.compensate` accept it. Optionally
  `predict_interval(commanded, setup, level) -> (lower, upper)` (bounds of dz
  [m]) and `assess(commanded, setup) -> dict` (out-of-distribution report);
  `compensate` reports both when present, and refuses a target or a
  compensated shape that `assess` puts outside the model's envelope unless
  `allow_out_of_envelope=True` (`--allow-out-of-envelope`). A model may
  declare `target`: "dz" (the total deviation - every `precomp.ml` bundle)
  or "residual" (a correction of a simulation); method "hybrid" adds the
  model to a simulation and so refuses a "dz" model, which would count the
  springback twice. `api.predict` and `api.compensate` also take the
  directory of a `precomp.ml` bundle as the model (loaded with
  `precomp.ml.registry.load_model`, imported only then).
* **Command line.** `precomp dataset ...`, `precomp train ...`,
  `precomp evaluate ...` and `precomp active ...` hand their full argument
  list to `precomp.ml.cli.main(argv) -> int`; `precomp compensate --method
  surrogate --model DIR` loads the model with
  `precomp.ml.registry.load_model(DIR)` and, with `--verify-fea` (the default
  for surrogate and hybrid; needs `--work-dir`), simulates the compensated
  part once; the report then shows that simulation, not the model's
  prediction, and states the envelope verdict. Without `precomp.ml` those
  commands exit with status 2 and say so.
* **Data.** `Part.sample(rng)` / `to_dict()`, `FormingSetup.to_dict()`,
  `simulate_many` (outcomes with the run hash, cache hit and runtime),
  `FormingResult.formed_surface()` / `forming_forces()` and the `Toolpath`
  pseudo-time are the inputs a data set records.

## Verification

Each row is a test in `python/tests` against an exact answer; the measured
value is from the test's own configuration.

| Check | Exact answer | Measured | Test |
|-------|--------------|----------|------|
| Drop cutter on a plane | c = z + R | bitwise equal | `test_toolpath` |
| Drop cutter on a 50 deg filleted cone, 120 centres, h = 0.25 mm | distance to the analytic surface = R | min - R = -4.2 um (bound ~6 um), exactly R at 10+ node contacts | `test_toolpath` |
| Contour and spiral paths | no point below the drop-cutter surface; spiral height non-increasing | gouge 0, height increase 0, contact to 1e-9 m | `test_toolpath` (30 part/tool combinations checked during development) |
| ICP, elliptic cone, 20 000 points, 1e-5 m noise, 2 % outliers up to 10 mm, (2, -1.5, 3) deg and (3, -2, 1.5) mm | the inverse transform | 1.1 um RMS point error, 0.0012 deg | `test_metrology` |
| Signed deviation of a cone offset by 0.3 mm along its normal | 0.3 mm | straight wall 0.86 um (bound h^2 tan(alpha) / 8 r = 1.5 um), fillets 13.6 um (bound h^2 / 8 rho cos^2 alpha = 20 um) | `test_metrology` |
| DA on formed = c - 0.2 S c (S a Gaussian smoothing) | contraction <= rho(G) = 0.2 per iteration | 0.194 per iteration, 5e-5 of the initial error after 6 steps | `test_compensation` |
| DA through the forming chain (test-double solver, formed = 0.9 c) | contraction 0.1 | 0.1 to 1e-6 relative | `test_api_cli` |
| One DA step from a tilted, lifted scan aligned on the flange | c - (f - t) | 3e-7 m (the in-plane position a flat flange cannot fix) | `test_compensation` |
| Wall-angle projection | angle <= 65 deg, never lower | 65.000 deg, raised only | `test_compensation` |
| Curvature of a spherical bowl (central differences) | H = 1/R, K = 1/R^2 | 1e-4 relative | `test_geometry` |
| Smoothed derivatives of a quadratic 15 mm below the sheet plane | its exact derivatives | 1e-9 relative (scipy's own filters: O(1) errors) | `test_geometry`, `test_parts` |
| Distance to an ellipse | brute force over 400 001 points | 1e-9 m (the sampling) | `test_parts` |
| Formed surface from a synthetic result (stretched, deflected sheet) | the deformed nodes | 1e-15 m | `test_fea_results` |
| STL round trip | the height map | 1e-15 m (ASCII), 2e-9 m (binary, float32) | `test_geometry` |
| DSIF support on a 40 deg cone, 4 mm balls, t = 1 mm | one centre per tool point; gap t cos(theta) along the normal where the sheet is formed; never inside the in-process sheet; below the part in the air | median gap error < 1 um, never squeezed below it; 0 points inside (1 um at the floor, the last level's offset) | `test_toolpath` |
| Reach of a ball from below | the underside where the ball fits; less in a slot narrower than it | 1e-6 m on a wide Gaussian bump; < 0.1 of a 1 mm high, 1 mm wide slot for a 4 mm ball | `test_fea_support` |
| Backing plate node ids | bottom nodes (z = -t) outside the opening, in the numbering of the mesh the solver writes | all, on `mesh.json` of the run | `test_fea_support` |
| DA with a rim pass on a sag of 1 mm | above the plane only in the swept band and strip, within the ball's reach (a fixed point of the bound), the rest of the flange held | raised 0.3+ mm with a 1.5 mm ball, less with a 4 mm one; 0 elsewhere | `test_fea_support` |

With the real `sparlab_form` (`test_integration_sparlab.py`, skipped without
`build/bin/sparlab_form`), against the contract of `docs/forming.md` and
physical sense rather than exact answers:

| Check | Expected | Measured | Test |
|-------|----------|----------|------|
| Every deck variant precomp writes: Hill48, von Mises and Chaboche materials; the three kinematics; contact and every solver key; clamped only; spiral; Tet4; backing plate (Hex8, Tet4); DSIF with squeeze, with the rim pass, clamped only | accepted by `--strict-config`, the analysis built | 14 of 14; a `contact` key in the tool refused (exit 2, naming it) | `test_integration_sparlab` |
| The test double and `sparlab_form` on one deck | the same files, CSV columns, summary and `mesh.json` keys, mesh and step windows | identical | `test_integration_sparlab` |
| Tiny SPIF: 20 x 20 x 1 mm AA5754-O blank, 8 x 8 x 2 Hex8, 4 mm tool, one spiral revolution to 1 mm ending in contact, unload, 3-2-1 release | depth about the tool's 1 mm; the sheet under the removed tool rises; no reaction after the release; the force on the tool upwards, of order 1 kN | depth 0.937 / 0.913 / 0.953 mm after form / unload / release; rise 25-45 um; release reactions 8e-13 N; fz >= 0 at all 51 increments, peak 1 070 N; 71 increments, 328 iterations, 3.1 s | `test_integration_sparlab` |
| The original small cone: 80 mm blank, 20 x 20 x 1 Hex8, two contours to 4 mm | completes, depth below twice the target's | completes, 270 increments, about 100 s | `test_integration_sparlab` |
| Tiny cone (24 mm blank, 12 x 12 x 2 Hex8, 3 mm tool, 1 mm deep) on a backing plate, 0.5 mm clearance | the plate pushed down in "form" only, never through; no reaction after the release | plate in contact 33 of 40 increments, fz -2.2 kN to 0, penetration 0.27 um; release reaction 4e-13 N; 16 s | `test_integration_sparlab` |
| The same cone with a DSIF support and the rim pass | tool pushed up, support down; the rim pass pushes up in its own step; no reaction after the release | tool fz 0 to 589 N, support -505 N to 0 (steps 1 and 3), penetration <= 0.42 um; release reaction 5e-13 N; 24 s | `test_integration_sparlab` |

## Limitations

* **The C++ contract.** The deck, the command line and the result files
  follow `docs/forming.md` (sections 2 and 3) as of this revision; every run
  uses `--strict-config`, so a key the solver stops reading fails loudly. The
  test double (`python/tests/fake_sparlab_form.py`) mirrors the contract and
  `test_integration_sparlab.py` checks it against the real executable, file
  by file and key by key, when `build/bin/sparlab_form` exists; a change of
  the C++ contract is a change of `precomp.fea.deck.forming_block`,
  `precomp.fea.results` and the double.
* **Forming model.** One spherical forming tool, and at most one support
  (a backing plate or a DSIF support ball); the path's final retract is not
  simulated (the "unload" step removes the tools). The node-to-surface
  contact resolves a tool poorly on elements not much smaller than it
  (`docs/forming.md`, section 6), and Hex8 sheets with few layers are stiff
  in bending: the numbers of a coarse deck are indicative. Where the tool
  circles inside its own radius (a small floor), the sheet inside the loop
  can end below the tool tip: one 1 mm contour of a 4 mm tool on a 3.3 mm
  radius over a 20 mm blank left the sheet 1.31-1.35 mm deep, on 8 x 8 x 2,
  16 x 16 x 2 and 16 x 16 x 4 Hex8 alike.
* **Materials.** The library values are nominal, not certified; Hill48 with
  r < 1 (the aluminium alloys) is known to underestimate the equibiaxial
  yield stress, which later criteria (Yld2000) correct and SparLab does not
  have; Swift and Hollomon laws are approximated by linear + Voce hardening
  (errors in the table above).
* **Geometry.** A pyramid's corner radius must be at least the horizontal run
  of its wall plus the bottom fillet's tangent length (so the corners are C1
  cones, not creases): sharp-cornered pyramids are not in the family. The
  numerical wall angle of `freeform` and `saddle` is measured on a refined
  grid to about 0.02 deg. Every part is centred on the blank.
* **Tool paths.** The drop cutter is exact on the grid; the gouge hidden
  between nodes is about h^2 / (8 R cos^3 theta). The flat floor is not
  swept. The spiral needs one loop per level (no islands, one pocket), and
  so does the rim pass.
* **Support.** The DSIF support follows the sine law (or the initial
  thickness), not the sheet's actual thickness: where the model's wall comes
  out thicker than t cos(theta) - about 20 % on the validation cone's
  2-layer mesh - the sine-law gap squeezes it, which the `squeeze` setting
  does not show; `thickness_law: "initial"` leaves the gap t. The support's
  springback and deflection, and a support ball that would touch the sheet
  away from its contact point in the air moves of contour paths, are not
  modelled. The backing plate is flat at the sheet plane and one element
  coarser than its clearance on the mesh; a raised plate or a die is not
  offered. Upward commands exist only where the rim pass reaches them.
* **Metrology.** Scans must be height fields of the tool side (no
  overhangs). A flat flange fixes only z and the two tilts; the rotation of an
  axisymmetric part about its axis is never observable. The point-cloud
  deviation fits local planes and is biased by about kappa r^2 / 4 on curved
  surfaces (r the radius of eight neighbours); gridding the scan first, as
  `compare_scan` does by default, is more accurate for dense scans.
* **Compensation.** The formability projection only raises the surface; a
  target steeper than the limit is formed shallower there, not refused. A
  scan with holes gets no update where it has no data.
* **Robot.** One constant compliance matrix (or one per point) is a
  first-order model; forces are interpolated in the path's pseudo-time, not
  in real time.
* **Thickness.** `thickness_map` is the straight distance between the top
  and bottom node of a mesh column, which is the sheet thickness only while
  the column stays close to the normal.
