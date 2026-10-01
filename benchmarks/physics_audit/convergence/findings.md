## Findings, ranked by their effect on the predicted final shape

Numbers are for the representative part; "reference" is 0.625 mm / 1 IM layer
/ 2 x 2 x 7, "limit" its Richardson extrapolation (Table 3). Evidence paths are
relative to this directory.

**1. The in-plane element size is the dominant discretisation error, and the
benchmark's 2 mm is not converged: ~0.10-0.12 mm RMS, 0.24 mm max on the
released shape.** Against the reference the benchmark's released surface
differs by 0.099 mm RMS / 0.242 mm max over the part (`diffs.csv`,
`h2_L2_std` vs `h0.625_L1_im_tp7`); with the reference's own error the bound
is 0.118 mm (Table 4). Its metrics are biased the same way: vertical RMS
deviation 0.767 against a limit of 0.674 mm (+0.093, 14 %), vertical max
1.426 against 1.212 (+0.215), rim sag -1.151 against -1.019 (0.13 mm too
much sag), the floor 0.20 mm too high at the centre (`fig_profile_error.png`,
`profiles.csv`). The IM family converges cleanly at second order over 2 / 1 /
0.625 mm (observed p = 1.98 for the RMS, 2.03 for the max, 1.99 for the rim
sag, 1.78 for the field; `richardson.csv`), so the extrapolation is
trustworthy for those quantities. Refining 2 -> 1 mm moves the released shape
by 0.070 (standard) to 0.107 mm RMS (IM), 1 -> 0.833 mm by 0.017 and
1 -> 0.625 mm by 0.025 mm (Table 2). The error sits where the sheet bends
sharply: the pillow of the untouched floor (r < 3 mm, +0.20 mm), the rim of
the part and the free ring up to the clamp (r = 5-15 mm, -0.04 to -0.12 mm).
Most of it is already in the *formed* shape (formed-shape change 0.099 mm RMS,
springback-field change 0.014 mm; Table 2), i.e. in the forming, not in the
unloading. Cause: a 4 mm tool on 2 mm elements touches 1.3 nodes on average
(at most 3; 3.1 on 1 mm; `cases/*/output/tool_forces.csv`), and a 2 mm
element cannot carry the bending localised under the tool and at the clamp
edge.

**2. At 2 mm the standard element's answer is right for the wrong reason.**
The standard Hex8 (2 x 0.5 mm cells) locks in bending: the incompatible-mode
element changes the 2 mm answer by 0.049-0.058 mm RMS, 0.24-0.28 mm max and
the rim sag by -0.055 to -0.067 mm (Table 2) - *away* from the converged
answer (IM 2 mm: field error bound 0.15 mm, the standard element 0.12 mm).
The stiffening of the locked element partly cancels the coarse mesh's
over-sag. At 1 mm the two technologies agree within 0.015 mm RMS / 0.054 max,
and both converge to the same answer. So the benchmark's 2 mm accuracy is a
cancellation of two errors whose balance will shift with the part (wall angle,
depth), the tool radius and the thickness; it cannot be carried to other
setups. More standard layers do nothing for it (2 -> 4 layers: 0.003 mm RMS,
0.017 max; 2.9 times the cost): as the channel-springback study of
Padmanabhan et al. (2007) found, the in-plane to thickness size ratio matters,
not the number of layers.

**3. The rim sag is physical in this model, not a mesh artefact.** It shrinks
from -1.151 (benchmark) to -1.038 (reference) and -1.019 mm (limit): 10-12 %.
The benchmark's conclusion that the band next to the rim dominates the error
(its `regions.csv`: 79-87 % of the squared error) survives refinement; for
this part the upper-band RMS goes 1.164 -> 1.052 mm (`results.csv`,
`upper_rms_mm`). Compensation strategies should be judged against the model's
converged sag (~1.0 mm here), not the 2 mm model's (1.15 mm).

**4. The springback itself (released minus formed) is small and converges
slowly, but its discretisation error is below 0.03 mm.** Springback RMS over
the part: 0.120 mm (benchmark), 0.135 (1 mm standard), 0.142 / 0.137 / 0.132 /
0.131 mm (IM at 2 / 1 / 0.833 / 0.625 mm) - no clean order (the IM values drift
by 0.005 mm per step), but the springback *field* changes by at most 0.014 mm
RMS / 0.025 mm max between the benchmark and the reference (Table 2). About
half of it is not the part springing back but the clamped frame bowing after
the release: the frame's mid-edges drop 0.06 mm (benchmark) to 0.10 mm
(reference) below the 3-2-1 corners, which lowers the whole window
(`analyze.py`'s surfaces; the springback profile is a near-uniform -0.11 to
-0.15 mm from the centre to the clamp, `profiles.csv`, `sb_tot_mm`). The
unclamping springback alone (release minus unload, `sb_rel`) is 0.151 ->
0.155 mm RMS; its IM series 0.186 / 0.164 / 0.155 converges at p = 0.8
(limit 0.136). For compensation this means: the springback is 0.12-0.14 mm of
a 0.7 mm deviation, and the numerical error of the *formed* shape (finding 1)
is larger than the whole springback's numerical error.

**5. Explicit forming: 0.017-0.019 mm from implicit on the final shape, but 13-20 %
less springback - and it is not cheaper here.** At 0.5 m/s with selective
dynamic mass scaling to 1 us (added mass 230-257 times the physical,
equivalent speed 7.5-8 m/s; kinetic ratio 0.006-0.020, penetration 0.5-0.8 %
of the element, energy balance 3e-4: every validity check of docs/forming.md
7.2 passes) the formed shape agrees with the implicit one to 0.007 mm RMS /
0.023 mm max, the released shape to 0.017 / 0.024 mm, the same on 2 and 1 mm
(Table 2) - but the springback is 0.105 against 0.120 mm (2 mm) and 0.113
against 0.135 mm (1 mm), 13-17 % less, four times the 3.5 % the smoke study of
docs/forming.md 7.3 measured at a similar speed. Cost: 604 s against 215 s
(2 mm), 2 463 s against 1 601 s (1 mm) - explicit is 1.5-2.8 times *more*
expensive on these 4 000-15 000 DOF sheets, and its fast kernel covers only
the standard 2 x 2 x 2 Hex8: the incompatible-mode element and a thickness
rule go through the generic element dispatch, which docs/forming.md 7.5 puts
at about 35 times the kernel's cost. The 0.625 mm explicit run (standard
Hex8, 8 192 elements, 38 025 DOFs, 6 617 s, rerun after the restart) confirms it: its released shape is 0.019 mm RMS / 0.060 mm max from the implicit IM reference (0.625 mm; element technologies differ, which at 1 mm accounts for 0.015 mm), vertical RMS 0.681 against 0.688, rim sag -1.029 against -1.038, but springback 0.105 against 0.131 mm (-20 %; the deficit grows with refinement: -13, -17, -20 % at 2, 1, 0.625 mm, `results.csv`) at 1.5 times the cost of the implicit reference (6 617 s against 4 441 s). Validity checks pass (kinetic ratio 0.0017, penetration 0.9 % of the element, energy error 5e-5, added mass 338 times the physical).

**6. Converged already (<= 0.008 mm RMS, <= 0.027 mm max on the released shape):**
one IM layer against two (2 x 2 x 3) at 1 mm (0.008 / 0.027 mm: one layer is
enough through the thickness),
tool travel per increment 1 -> 0.5 mm (0.002 / 0.012 mm; the 0.5 mm run needed
fewer cuts, 4 against 9, and was not slower), penalty scale 10 -> 30
(0.002 / 0.015 mm; the largest penetration 0.6 -> 0.2 um) and 3 -> 10
(0.004 / 0.023 mm; 2.0 -> 0.6 um), Newton tolerances 1e-6 -> 1e-8
(0.003 / 0.016 mm; 2 cuts against 9, and not slower), release increments
10 -> 40 (1e-5 mm: the unloading is elastic), thickness
points 5 -> 7 in the IM layer (0.002 / 0.004 mm at 2 mm, 0.003 / 0.025 mm at
1 mm), and 5 instead of 2 points through each standard layer (0.002 / 0.008 mm) (Table 2). None of the solver controls is worth tightening;
the tolerance budget belongs to the mesh. The run time at a given mesh is set
by the Newton iterations and cuts, not by these settings' nominal cost: the
softer penalty 3 needed 1 094 iterations and 2 cuts against 1 906 and 9 at
the default 10 (77 s against 215 s), the stiffer 30 needed 2 844 and 10.

**7. The metrology interpolation of a coarse mesh costs 0.04 mm RMS by
itself.** `formed_surface` interpolates the nodes linearly on flat triangles.
Representing the 1 mm solution by its every-second node (a 2 mm node grid)
changes it by 0.042 mm RMS / 0.160 mm max; a bicubic spline through the same
nodes cuts that to 0.015 / 0.063 mm (`diffs.csv`, rows `[every 2th node]` and
`[every 2th node, bicubic]`), at the cost of a second. On the benchmark mesh it
removes only a slice of the error (0.078 -> 0.067 mm against the 1 mm
solution), because the nodal solution itself is off (finding 1); on 1 mm and
finer meshes the flat-triangle error is 4 times smaller and does not matter.

**8. Not converged at 2 mm, outside the shape metrics:** the mean tool force
(566 N benchmark, 528 N reference, IM 416 / 522 / 528 N at 2 / 1 / 0.625 mm;
Richardson limit 529 N), relevant to the robot's force feedback and load
predictions; the peak equivalent plastic strain (0.38 benchmark, 0.46-0.55 at
1 mm, 0.59 at 0.625 mm: localised, still rising), relevant to any thinning or
fracture check; and the deepest point (formed depth 3.082 -> 3.094 -> 3.108
mm, rising at a constant rate: the tool's imprint on the floor is resolved
only as the mesh approaches the contact width).

**9. A second, steeper part converges the same way (added on resume).**
`truncated_cone-s2026-0001` (wall 54.0 deg, depth 2.63 mm, top fillet
2.33 mm, bottom fillet 1.49 mm), IM 2 x 2 x 5 at 2 / 1 / 0.833 / 0.625 mm
(Table 5, `results_c1.csv`, `diffs_c1.csv`, `richardson_fit_c1.csv`). The
benchmark setup (reproduced bit for bit from its run cache) is 0.080 mm RMS /
0.256 mm max from the 0.625 mm solution; RMS deviation 0.703 against a limit
of 0.642-0.646 mm, rim sag -1.122 against -1.026..-1.036 mm. The 0.833 mm
IM setup is 0.012 / 0.044 mm from the finest mesh (field error with the
finest mesh's own 0.009-0.012: 0.021-0.024 mm), its RMS error +0.011..+0.015,
rim sag -0.016..-0.027. The worst-point deviation converges erratically
(p = 1.3-6 by triple; limit 1.107-1.176 mm): 0.833 mm is +0.011..+0.080 from
it. The springback on this part is larger (0.150 benchmark, 0.165 mm
converged): the benchmark under-predicts it by 9 % (part 1: 0.120 against
0.131, -8 %), the 2 mm IM element over-predicts it on both parts.

**10. The Richardson limit itself carries +-0.01-0.03 mm.** The observed
order on part 1 depends on which three meshes are used: p = 1.6 (2/1/0.833),
1.8-2.1 (triples with 2 mm and a fine mesh), 2.5-3.3 (meshes <= 1 mm), up to
5.8 for the max (`richardson_fit.csv`, 10 triples and 2 fits after adding
0.714 mm). The limits spread over 0.659-0.682 (RMS), 1.187-1.238 (max),
-0.997..-1.030 mm (rim sag) (Table 6). The 2 mm mesh is pre-asymptotic, so
the triples without it (the upper ends of the bands) are the more credible;
the conclusions do not change, but the earlier single-triple "errors"
(Table 4) are central estimates, not bounds.

**11. Penalty 3 halves the cost at the price of ~one mesh step.** On the
0.833 mm IM mesh penalty 3 changes the released shape by 0.009 / 0.021 mm
(part 1) and 0.005 / 0.020 mm (part 2), always towards more sag, and runs
in 798 against 1 586 s and 1 172 against 1 965 s (half the Newton
iterations: 1 651 against 3 351). At 0.714 mm penalty 3 (1 448 s) is as
accurate as 0.833 mm with penalty 10 (1 586 s): no net gain on part 1
(Table 1-2). Use it as a cost lever only where the 0.005-0.01 mm is affordable.

