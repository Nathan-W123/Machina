# Displacement adjustment with SparLab's forming FEA: a small cone

The first physics numbers of `precomp`: every formed shape here is a
`sparlab_form` simulation (implicit, quasi-static, finite-strain plasticity),
not the proxy of `benchmarks/ml_proxy`. They are **simulations of a coarse
model with nominal material data, not validated against a formed part**.

Reproduce (two simulations per run, about 4.5 and 7 minutes on one thread of
a shared 4-core machine):

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release && cmake --build build -j2
python3 python/scripts/fea_da_demo.py --out benchmarks/fea_da_cone/clamp_9mm --work <scratch> --clamp-margin 0.009
python3 python/scripts/fea_da_demo.py --out benchmarks/fea_da_cone/clamp_5mm --work <scratch2> --clamp-margin 0.005
```

The script runs the `precomp` command as a user would - `precomp part`,
`precomp setup`, `precomp compensate --method fea --iterations 2` - and
writes, per run, what it ran and measured: `run.json` (the three commands,
the precomp and SparLab versions, runtimes, load average), `setup.json`,
`part.json`, `compensation.json` (as `precomp compensate` wrote it),
`history.csv` (per iteration: deviation over the part, the simulation's
cost) and `profiles.csv` (azimuthal means every 0.5 mm of radius). The
simulations are deterministic: the `clamp_5mm` record was written by the
final script from the run cache of its first run (`cache_hit` in
`history.csv`; that run's `precomp compensate` took 427 s), `clamp_9mm`
from scratch (270 s).

## Setup

| | |
|-|-|
| Target | truncated cone, top radius 9 mm, wall 45 deg, depth 3 mm, 2 mm rim and floor fillets (footprint 19.7 mm across); 161 x 161 grid at 0.25 mm |
| Blank | 40 x 40 x 1 mm AA5754-O from the library: E 70 GPa, Hollomon K = 420 MPa, n = 0.30 from 100 MPa, fitted by linear + Voce hardening; Hill48 with r0 / r45 / r90 = 0.75 / 0.70 / 0.80 - nominal handbook values |
| Mesh | 20 x 20 x 2 Hex8 (2 mm in plane, 0.5 mm through the thickness), 3 969 DOFs; `finite_logarithmic` kinematics, mean dilatation |
| Clamp | the frame within 5 mm of the blank edge (a 30 mm square window) or 9 mm (a 22 mm window, 1.2 mm beyond the part at the mid-sides, like a backing plate) |
| Tool and path | sphere of 4 mm radius, friction 0.1, penalty scale 10; a spiral with 1 mm step-down and 1 mm point spacing, 1 mm tool travel per increment |
| Steps | form (clamped) - unload (tool removed, clamped) - release onto 3-2-1 supports at three corners of the former clamp |
| Compensation | displacement adjustment, vertical error, alpha = 1, no smoothing, the flange held at z = 0, commanded surface kept at z <= 0 and within 65 deg; two iterations = two simulations (the target as commanded, then once compensated) |

## Results

Deviation of the released part from the target, over the part (the nodes
where the target is below the sheet plane, 4 809 nodes), in mm:

| Clamp | Iteration | Vertical RMS | Vertical max | Normal RMS | Normal max | Bias | Wall RMS | Flange RMS | Simulation |
|-------|-----------|-------------:|-------------:|-----------:|-----------:|-----:|---------:|-----------:|------------|
| 9 mm | 0 (target) | **0.465** | 0.970 | 0.439 | 0.903 | -0.220 | 0.538 | 0.158 | 143 s, 338 increments, 2 162 iterations, 6 cuts |
| 9 mm | 1 (compensated) | **0.405** | 0.833 | 0.388 | 0.796 | -0.118 | 0.470 | 0.144 | 126 s, 295 increments, 1 895 iterations, 7 cuts |
| 5 mm | 0 (target) | **0.739** | 1.392 | 0.728 | 1.332 | -0.460 | 0.866 | 0.391 | 249 s, 432 increments, 2 852 iterations, 16 cuts |
| 5 mm | 1 (compensated) | **0.630** | 1.196 | 0.629 | 1.163 | -0.287 | 0.733 | 0.347 | 176 s, 286 increments, 1 989 iterations, 8 cuts |

One compensation step lowers the RMS deviation by 13 % (backing-plate
clamp: 0.465 to 0.405 mm) and 15 % (wide clamp: 0.739 to 0.630 mm), and the
largest by 14 %. Largest equivalent plastic strain 0.32-0.46, peak tool
force 1.0-1.3 kN, at most three nodes in contact, penetration below
0.7 um; the release leaves no reaction (below 1e-12 N). Every run warns
that the clamp of the "unload" step holds reactions (1.0-1.3 kN): that is
the springback the clamp prevents, which the release lets go.

Azimuthal means with the backing-plate clamp (mm; formed = after the
release unless marked):

| r | target | formed (form step) | formed | commanded, iteration 1 | formed, iteration 1 |
|--:|-------:|-------------------:|-------:|-----------------------:|--------------------:|
| 0 | -3.000 | -2.480 | -2.588 | -3.412 | -2.778 |
| 2 | -3.000 | -2.622 | -2.731 | -3.269 | -2.908 |
| 4 | -3.000 | -2.894 | -3.001 | -2.999 | -2.981 |
| 6 | -2.817 | -2.464 | -2.570 | -3.064 | -2.417 |
| 7 | -1.992 | -1.968 | -2.072 | -1.913 | -1.856 |
| 8 | -0.997 | -1.382 | -1.483 | -0.510 | -1.295 |
| 9 | -0.176 | -0.901 | -0.999 | 0.000 | -0.882 |
| 10 | 0.000 | -0.528 | -0.622 | 0.000 | -0.552 |
| 12 | 0.000 | -0.108 | -0.190 | 0.000 | -0.178 |
| 16 | 0.000 | 0.000 | -0.049 | 0.000 | -0.053 |

## What the numbers say

* **Displacement adjustment contracts slowly on this process.** A process
  map close to the identity would lose most of its error in one step; this
  one keeps 85-87 %. The profiles show three reasons:
  * *The rim sags.* Between the clamp and the first contour the sheet bends
    down, 0.8 mm below the target at r = 9 mm (1.3 mm with the wide
    clamp). Correcting it would
    command the surface above the sheet plane, which neither DA (it keeps
    z <= 0) nor a tool pressing from above can do: the command is 0 there
    and the sag stays (-0.88 mm after compensation).
  * *The floor is not swept.* The centre came out 0.41 mm shallow;
    commanding it 0.41 mm deeper moved it 0.19 mm - the tool reaches only
    the floor's edge.
  * *The response is not local.* Raising the upper wall moves every
    contour inwards, and at r = 6 mm the compensated part came out 0.15 mm
    shallower although it was commanded 0.25 mm deeper.
  The release onto the corner supports also lowers the part by about
  0.1 mm (the flange by 0.02-0.05 mm), an offset that vertical DA without
  registration keeps chasing.
* **The unsupported flange matters more than the compensation.** Clamping
  close to the part, as a backing plate does, cuts the part's deviation
  before any compensation by 37 % (0.739 to 0.465 mm), and the flange's by
  60 % (more of it is clamped).
* **What would help**, in order: more DA iterations (the error was still
  falling), smoothing of the update (the non-local response),
  registration of the released part on its flange before measuring (the
  release offset), and a finer mesh - 2 mm elements under a 4 mm tool carry
  the load on at most three nodes (`docs/forming.md`, section 6).

## Limitations

Coarse model (two Hex8 layers, stiff in bending; 2 mm elements); nominal,
uncertified material data, and Hill48 with r < 1; one tool path style; a
3 mm deep part on a 40 mm blank; no experiment. The runtimes were measured
on one thread of a shared machine (load average 3-4) and vary by tens of
percent between runs.
