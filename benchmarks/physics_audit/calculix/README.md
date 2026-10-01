# Physics audit - lens: independent code cross-check (CalculiX 2.21)

STATUS: complete (resumed after a machine restart; every finished run of the
first auditor was reused, the interrupted 3D plate runs were redone).
Re-run everything with `run_all.sh`; run log `logs/run_log.txt`. Bulky solver
outputs (CalculiX .frd/.dat, VTK) are git-ignored; every number below is in a
committed `result.json` / `*.csv`.

## Summary (ranked by effect on the predicted released shape)

| # | finding | size on the shape | evidence |
|---|---|---|---|
| 1 | SparLab's forming mechanics (finite-strain J2 + Voce, penalty contact on a rigid tool, Coulomb friction, tool-removal unload, release onto determinate supports) reproduce CalculiX on the identical discrete problem | 3D plate with plunge + drag: released shape <= 2.3 um max, springback field 0.6 um rms; forces <= 0.7 % (plunge), 0.9 % / 3.6 % (drag Fz / Fx) | `punch/plate_compare_plate_h2_*_drag3.csv`, `punch/strip_h2_*/result.json` |
| 2 | Benchmark mesh (2 mm in-plane) discretisation error | 3D test: 0.038 mm rms / 0.089 mm max released shape, springback under-predicted by 16 % (B-bar Hex8); strip: 0.08 mm rms / 0.10 mm max vs C3D20R | `punch/plate_mesh_study.csv`, `strip_reference/compare_ref_h0.25_nz8_s2s_retract.csv` |
| 3 | Full-integration Hex8 without B-bar (CalculiX C3D8, SparLab `mean_dilatation: none`) locks plastically - the benchmark's default `auto` avoids it | 0.5 mm (strip), 0.14 mm rms (plate) | same |
| 4 | Residual difference SparLab (fine linear hexes, node-wise penalty) vs CalculiX C3D20R surface-to-surface reference | 6 % punch force, 0.010-0.013 mm rms / 0.03 mm max released profile | `strip_reference/compare_ref_h0.25_nz8_s2s_retract.csv` |
| 5 | Log-additive vs multiplicative large-strain plasticity | stress <= 0.3 % (uniaxial to 0.5 and back), <= 0.5 % in shear to gamma 0.4, 2.2 % at gamma = 1 | `homogeneous/results_inc400.csv` |
| 6 | Tool travel per increment 1 mm (benchmark) vs 0.1 mm; 5 vs 2 thickness points; 4 vs 2 layers | <= 0.004 mm | plate / strip tables |
| 7 | Not cross-checkable in CalculiX: Hill48 plasticity, Chaboche / Prager kinematic hardening (no built-in CalculiX counterpart; finite-strain *PLASTIC needs isotropic elasticity) | - (material lens) | CalculiX manual *PLASTIC |


## Setup common to every test (`common.py`)

* Same nodes in both codes (structured hex grids built in `common.Grid`);
  SparLab decks in SI, CalculiX decks in mm/N/MPa.
* Material: the benchmark's AA5754-O flow curve (sigma_y0 = 100 MPa, linear
  H = 237.6 MPa + Voce Q = 123.9 MPa, b = 14.38), E = 70 GPa, nu = 0.33,
  **J2** (CalculiX has no built-in Hill48 plasticity; its large-strain
  plasticity needs isotropic elasticity - CalculiX manual, *PLASTIC,
  https://web.mit.edu/calculix_v2.7/CalculiX/ccx_2.7/doc/ccx/node232.html).
  CalculiX gets the flow curve as a dense *PLASTIC table (572 points).
* SparLab: `finite_logarithmic` (Hencky strain, additive split).
  CalculiX: `NLGEOM` + isotropic elasticity -> Simo's multiplicative
  finite-strain J2 (same manual page).
* Contact: SparLab penalty on an analytic sphere/cylinder,
  kappa = 10 E / h_elem; CalculiX NODE TO SURFACE, LINEAR pressure-overclosure
  of the same slope, Coulomb mu = 0.1 with stick slope = kappa; the punch is a
  meshed C3D8 shell whose every node is prescribed (rigid).
* Elements: SparLab `standard` Hex8 has B-bar (mean dilatation) by default
  (`mean_dilatation: auto`); SparLab `std_nobbar` = `mean_dilatation: none`
  corresponds to CalculiX C3D8 (full 2x2x2 integration);
  SparLab `incompatible_modes` corresponds to CalculiX C3D8I.

## Test 3 - homogeneous large-strain plasticity (DONE)

`case_homogeneous.py`, one element, results `homogeneous/results_inc400.csv`
(CalculiX with 400 increments per step; `results_inc20.csv` 20, same answer).

| path | quantity | SparLab vs CalculiX |
|---|---|---|
| uniaxial, log strain 0 -> 0.5 | von Mises (Cauchy) | within 0.29 % at every step (340.84 vs 341.82 MPa at 0.5) |
| uniaxial, reversed 0.5 -> 0 | von Mises | within 0.26 % (458.6 vs 457.4 MPa at the end) |
| uniaxial | PEEQ | within 0.0013 absolute (0.4951 vs 0.4968 at 0.5) |
| simple shear gamma <= 0.4 | von Mises | within 0.5 % |
| simple shear gamma = 0.6 / 0.8 / 1.0 | von Mises | -0.40 % / -1.09 % / -2.19 % (352.1 vs 359.9 MPa) |

Reason for the shear drift: in simple shear the principal stretch axes rotate
(31.7 deg at gamma = 1). SparLab's yield function acts on the stress
conjugate to the Hencky strain; with rotating axes its Cauchy stress is not
coaxial with it, so the Cauchy Mises stress sits 2 % below the flow stress at
SparLab's own PEEQ (352.1 vs 359.1 MPa), while CalculiX's multiplicative
model returns Cauchy Mises = flow stress at its PEEQ exactly. Both are
legitimate finite-strain models; the difference is a property of the
log-additive formulation, not a bug. At the strains of the benchmark parts
(largest PEEQ ~ 0.3, through-thickness shear well below gamma = 0.5) it is
< 0.5 % in stress.

## Test 1 + 2 - punch forming of a clamped strip, unload, release onto determinate supports (h = 2 mm, 2 layers: DONE)

`case_punch.py`, kind=strip: plane-strain strip 40 x 1 mm, clamped |x| >= 15,
cylindrical punch R = 4 mm pushed 4 mm, mu = 0.1; unload (tool off, clamped);
release onto a pin + roller at x = -18 / +18 (top). Benchmark mesh density:
2 mm in-plane, 2 layers of 0.5 mm.

| case (SparLab / CalculiX) | punch force at 4 mm, N/mm | max rel. force diff along plunge | top-surface dz after form / unload / release (max) | springback (unload) |
|---|---|---|---|---|
| `incompatible_modes` / C3D8I | 94.40 / 94.40 | 0.33 % | 0.0004 / 0.0005 / 0.0004 mm | 0.1451 / 0.1450 mm |
| `std_nobbar` / C3D8 | 133.04 / 133.23 | 0.23 % | 0.0011 / 0.0009 / 0.0047 mm | 0.1962 / 0.1960 mm |
| `standard` (B-bar, the benchmark element) / C3D8 | 97.83 / 133.23 | 32.7 % | 0.11 / 0.13 / 0.54 mm | 0.144 / 0.196 mm |

Evidence: `punch/strip_h2_*_mu0.1/result.json`, `force.csv`, `shape_*.csv`.

Reading: where the two codes solve the same discrete problem (same element
technology), SparLab reproduces CalculiX to < 0.4 % in force and < 5 um in the
formed, unloaded and released shapes, i.e. the contact, friction,
large-strain plasticity, unload and statically-determinate release are
implemented consistently. The benchmark's B-bar Hex8 has no CalculiX
counterpart; it lands within 4 % (force) and 6 um (released deepest point,
-3.913 vs -3.907 mm) of the incompatible-modes element, while the plain
full-integration Hex8 (CalculiX C3D8) locks: 41 % too stiff and a released
deepest point 0.54 mm shallower. Which of these is right is answered by the
element-independent reference below.

## Test 1 + 2 in 3D - spherical punch, plunge + drag with friction, unload, 3-2-1 release (DONE)

`case_punch.py kind=plate`: the benchmark's blank (40 x 40 x 1 mm, clamped
outside the 30 x 30 window), its mesh (2 mm, 2 layers), its tool (sphere
R = 4 mm, mu = 0.1); punch plunged 3 mm at the centre, then dragged 3 mm along
+x at depth (friction engaged, the ISF tool motion); unload with the tool
removed (clamped); release onto the benchmark's 3-2-1 supports at the nodes
nearest (+-17.5, +-17.5, 0) (here (+-18, +-18)). Table:
`punch/plate_compare_plate_h2_im_mu0.1_drag3.csv`; per-run
`punch/plate_h2_im_mu0.1_drag3/result.json`.

| quantity | SparLab IM | CalculiX C3D8I | diff |
|---|---|---|---|
| punch force end of plunge (3 mm) | 1216.1 N | 1213.9 N | +0.18 % (max 0.71 % along the plunge) |
| drag: mean vertical force | 764.8 N | 772.1 N | -0.9 % |
| drag: mean horizontal (friction + ploughing) force | 370.0 N | 357.1 N | +3.6 % |
| top-surface z after form / unload / release | | | max 1.7 / 1.9 / 2.3 um, rms 0.5 / 0.5 / 0.6 um |
| springback field (release - form), rms | 0.0669 mm | 0.0666 mm | rms of difference 0.6 um |
| largest PEEQ | 0.284 | 0.283 | 0.5 % |

Full-integration pair (SparLab `std_nobbar` / CalculiX C3D8, same test,
`punch/plate_compare_plate_h2_nobbar_mu0.1_drag3.csv`): plunge force 1858.7 /
1858.8 N, drag Fz 1360 / 1341 N (+1.4 %), drag Fx 617 / 644 N (-4.1 %),
released shape 0.8 um rms / 3.5 um max, springback field 0.8 um rms.
The drag (tangential) force is the least-matched quantity (3.6-4.1 %): both
codes regularise stick by a tangential penalty but CalculiX's node-to-surface
friction acts on the meshed (faceted, 32 x 32 cubed-sphere) punch, SparLab's
on the exact sphere; its effect on the shapes is < 3 um here.
Run times: SparLab 60 s (IM) / 34 s, CalculiX 1073 s / 1147 s (1 thread each).

The 6-mm drag of the first auditor (`plate_h2_im_mu0.1_drag`,
`plate_h2_nobbar_mu0.1_drag`) could not be completed in CalculiX: C3D8I
diverged at 3.8 mm of drag (increment below minimum) and the C3D8 run, after
completing the 6 mm drag, diverged at the first increment of the punch
retraction. Unloading by removing the contact pair (`*MODEL CHANGE, TYPE=CONTACT
PAIR, REMOVE`, the same operation as SparLab's tool-removal unload) is now the
default of `case_punch.py`; on the strip it gives the same unloaded shape as
the retraction to 0.02 um (C3D8, `strip_h2_std_nz2_tp5` vs
`strip_h2_std_mu0.1`) and 0.1 um (C3D8I, `strip_h2_im_mu0.1_remove`).

## Element-independent reference for the strip: CalculiX C3D20R (DONE)

`strip_reference.py`: the same strip, punch, material, steps and supports in
CalculiX C3D20R (quadratic, reduced integration), SURFACE TO SURFACE contact
(node-to-surface with quadratic slaves diverged at the first contact:
`strip_reference/ref_h0.5_nz4`, `ref_h0.25_nz8`). Two references:
0.5 mm x 4 layers (`ref_h0.5_nz4_s2s`, contact removed for the unload) and
0.25 mm x 8 layers (`ref_h0.25_nz8_s2s_retract`, punch retracted - the
contact removal made the fine model's plasticity update fail,
`ref_h0.25_nz8_s2s`). The reference is converged: the two differ by
0.005 mm rms / 0.015 mm max in the released profile and 0.7 % in force
(106.4 vs 105.7 N/mm at 4 mm).
Tables against the finer one: `strip_reference/compare_ref_h0.25_nz8_s2s_retract.csv`
(profile |x| <= 19 mm) and `strip_reference/regions_ref_h0.25_nz8_s2s_retract.csv`
(core |x| < 8, wall 8..15, flange 15..19); the same against the 0.5 mm one
`*_ref_h0.5_nz4_s2s.csv`.

| SparLab case | F at 4 mm vs ref | released profile rms / max vs ref | core / wall rms | deepest released point (ref -3.869 mm) |
|---|---|---|---|---|
| **benchmark: `standard` (B-bar), 2 mm, 2 layers** | -8.1 % | 0.080 / 0.101 mm | 0.084 / 0.094 | -3.913 |
| same, 1 mm max tool travel (benchmark) instead of 0.1 | -8.2 % | 0.079 / 0.099 mm | | -3.911 |
| same, 5 thickness points / 4 layers + 5 points | -8.1 / -7.8 % | 0.080 / 0.101, 0.077 / 0.103 mm (vs 0.5 mm ref) | | -3.913 / -3.911 |
| `incompatible_modes`, 2 mm, 2 / 4 layers | -11.3 / -11.2 % | 0.083 / 0.129, 0.082 / 0.125 mm | 0.065 / 0.112 | -3.907 / -3.903 |
| `standard`, 1 mm, 4 layers | -2.5 % | 0.025 / 0.049 mm | 0.029 / 0.028 | -3.901 |
| `incompatible_modes`, 1 mm, 2 / 4 layers | -4.6 % | 0.019 / 0.044, 0.019 / 0.042 mm | 0.026 / 0.017 | -3.914 / -3.911 |
| `standard`, 0.5 mm, 4 layers | -5.3 % | 0.013 / 0.032 mm | 0.018 / 0.011 | -3.901 |
| `incompatible_modes`, 0.5 mm, 4 / 8 layers | -6.2 / -6.1 % | 0.010 / 0.031 mm | 0.014 / 0.008 | -3.900 |
| `std_nobbar` (= C3D8), 2 mm, 2 layers | +25.0 % | 0.209 / 0.496 mm | | -3.373 |

Reading:
* On the benchmark's mesh the dominant discretisation error is the
  **in-plane** element size (2 mm against a 4 mm tool radius): ~0.08 mm rms /
  0.10-0.13 mm max in the formed and released profile, largest in the wall
  (0.09-0.11 mm rms). 1 mm gives ~0.02 mm, 0.5 mm ~0.01 mm. Layers (2 -> 4)
  and thickness points (2 -> 5) change the released shape by <= 0.004 mm, the
  1 mm tool travel per increment by <= 0.003 mm.
* The benchmark element (B-bar Hex8) and the incompatible-modes Hex8 converge
  to the same answer (at 0.5 mm: deepest released point within 0.001 mm of each other, both 0.010-0.013 mm rms from the reference). The plain
  full-integration Hex8 (C3D8 / `mean_dilatation: none`) locks plastically
  (0.5 mm error, +25 % force) and does not cure with refinement; SparLab's
  default (`auto`) avoids it.
* What remains at 0.5 mm - SparLab 6 % below the reference in punch force and
  0.03 mm deeper at the released bottom (0.010-0.013 mm rms) - is not a
  SparLab defect: CalculiX's own C3D8I on the same 0.5 mm x 8 mesh with
  node-to-surface contact agrees with SparLab to 0.2 % in force and 0.7 um in
  shape (`punch/strip_h0.5_im_nz8/result.json`). It is the difference between
  linear hexes with node-wise penalty contact and quadratic hexes with
  surface-to-surface (segment-integrated) contact around a 4 mm-radius tool -
  a modelling uncertainty of the contact patch, small next to the 2 mm mesh
  error.

## 3D plate: element / mesh / increment study (DONE)

`plate_mesh.py` -> `punch/plate_mesh_study.csv` (log `punch/plate_mesh_study.log`).
Same plate test (plunge 3 mm + drag 3 mm, unload, 3-2-1 release). Both codes
agree on the identical discrete problem (previous section, < 2.3 um), so the
mesh study is run in SparLab only; every case is sampled at the 2 mm grid's
top nodes and compared with the finest run (`standard`, 1 mm, 4 layers).
Regions: core max(|x|,|y|) < 8 mm, wall 8..15, rim (outside the clamp window) >= 15.

| SparLab case | plunge F | drag Fx | springback rms (release - form) | released shape vs finest: rms all / core / wall / rim, max | springback field error rms | runtime |
|---|---|---|---|---|---|---|
| **benchmark: `standard`, 2 mm, 2 layers, 0.1 mm travel** | 1259 N | 467 N | 0.0527 mm | 0.038 / 0.064 / 0.049 / 0.010, 0.089 mm | 0.015 mm | 43 s |
| same, 1 mm travel (benchmark value) | 1259 | 468 | 0.0526 | 0.038 / 0.063 / 0.049 / 0.010, 0.089 | 0.015 | 51 s |
| same, 5 thickness points | 1259 | 470 | 0.0526 | 0.038 / 0.064 / 0.049 / 0.010, 0.090 | 0.015 | 101 s |
| `incompatible_modes`, 2 mm, 2 layers | 1216 | 370 | 0.0669 | 0.042 / 0.068 / 0.055 / 0.013, 0.107 | 0.012 | 60 s |
| `std_nobbar` (C3D8), 2 mm | 1859 | 617 | 0.1449 | 0.140 / 0.207 / 0.179 / 0.063, 0.274 | 0.086 | 34 s |
| `standard`, 1 mm, 2 layers | 1258 | 427 | 0.0587 | 0.0045 / 0.008 / 0.005 / 0.003, 0.014 | 0.005 | 423 s |
| `incompatible_modes`, 1 mm, 2 layers | 1222 | 394 | 0.0647 | 0.009 / 0.015 / 0.009 / 0.008, 0.025 | 0.006 | 332 s |
| `standard`, 1 mm, 4 layers (finest) | 1267 | 437 | 0.0625 | - | - | 1130 s |

Reading: on the benchmark's 2 mm mesh the released shape carries ~0.04 mm rms
(0.09 mm max) of discretisation error in this small test and the springback
field is off by 0.015 mm rms - about 1/4 of the springback itself (the B-bar
Hex8 under-predicts springback by 16 %, the incompatible-modes Hex8
over-predicts it by 7 %); the drag force differs by up to +/-15 %. A 1 mm
in-plane mesh brings these to 0.005-0.009 mm and 0.005-0.006 mm at 7-10x the
runtime. Thickness points and the 1 mm tool travel per increment do not
matter (<= 0.0015 mm).


## What this means for the benchmark (0.73 mm RMS, 0.9-1.1 mm rim sag)

* Verified fact (this lens): the solver's mechanics are not a source of
  error at the 0.01 mm level. Where SparLab and CalculiX solve the same discrete
  problem, the released shapes agree to <= 2.3 um and the springback fields to
  < 1 um rms.
* Verified fact: the benchmark's 2 mm in-plane mesh adds discretisation error:
  0.04 mm rms (3D plunge+drag test) to 0.08 mm rms (plane-strain strip) in the
  released shape, and the B-bar Hex8 under-predicts springback by 16 % in the
  3D test.
* Inference, not measured on the benchmark parts: if the 16 % springback
  error carries over to the benchmark's ~1 mm rim sag, it is worth roughly
  0.1-0.2 mm. The convergence lens should confirm this on the actual parts at
  1 mm.
* Not covered by this lens: Hill48 and kinematic hardening (CalculiX has no
  built-in equivalent), the material data themselves, and the process model
  (one tool against a free sheet, where the real process uses two robots).
  On present evidence these, not the numerics, decide whether the 0.73 mm
  matches a real part.

## Recommendations

1. Make 1 mm in-plane the production / validation mesh (`element_size: 0.001`,
   2 layers). In the 3D test it cut the shape error from 0.038 to 0.0045 mm rms
   and the springback error from 0.015 to 0.005 mm rms. Cost: 7-10x runtime
   (43 s -> 423 s here). At minimum, re-run 2-3 benchmark parts at 1 mm and
   check how much the uncompensated RMS and the rim sag move.
2. Keep `mean_dilatation: auto`, and never use `none` with the standard
   Hex8. Either element is fine at 1 mm. At 2 mm, B-bar and incompatible
   modes bracket the answer (springback -16 % / +7 %).
3. Keep `max_tool_travel` at 1 mm and the default thickness points. Neither
   changes the shape by more than 0.004 mm here, and both cost runtime.
4. Do not use CalculiX as a production cross-solver for full toolpaths. It
   was 18-33x slower and failed to complete a 6 mm drag. Keep the small
   decks in this directory as a regression check (about 20 min per 3D case).
5. Hill48 and Chaboche can only be validated against material tests:
   uniaxial tests at 0/45/90 deg for r-values, and a tension-compression or
   cyclic bending test for the Bauschinger effect. CalculiX cannot provide
   that check.

## Sources

* CalculiX *PLASTIC (large-strain theory for isotropic elasticity only;
  hardening options): https://web.mit.edu/calculix_v2.7/CalculiX/ccx_2.7/doc/ccx/node232.html
* CalculiX C3D8I (incompatible modes remove shear locking, reduce volumetric
  locking): https://www.feacluster.com/CalculiX/ccx_2.18/doc/ccx/node30.html
* CalculiX *MODEL CHANGE, TYPE=CONTACT PAIR, REMOVE (used for the unload):
  https://www.feacluster.com/CalculiX/ccx_2.18/doc/ccx/node299.html
