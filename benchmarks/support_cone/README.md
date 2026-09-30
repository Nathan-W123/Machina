# Rim support on one small cone: backing plate and DSIF against single-point forming

Every formed part here is a **SparLab simulation** (`sparlab_form`: implicit,
quasi-static, `finite_logarithmic` kinematics, Hill48 plasticity with
**nominal handbook data** for AA5754-O). Nothing was formed or measured on a
machine. The numbers compare support strategies *on this model* - one
part, one process point, one run each - not on real parts
([Limitations](#limitations)).

The question: the springback benchmark (`benchmarks/springback`) found that
79-87 % of the squared deviation left after any compensation is a rim sag -
near the rim every part came out 0.9-1.1 mm too deep, because a single tool
pushing from above cannot hold the sheet between the clamp and the first
contour, and a command kept at z <= 0 cannot ask for it to come up. Does
support from below - `FormingSetup.support`, `docs/precomp.md`
("Rim support") - remove it?

## Headline

The benchmark's smallest held-out cone (truncated_cone-s2026-0000: 15.4 mm
across, 3.09 mm deep, 38.6 deg wall, 1.08 / 2.08 mm fillets), formed as the
target (1 FE run) and after one step of displacement adjustment with the FE
model in the loop (+DA1, 2 FE runs). Deviation of the released part from
the target over the part [mm]; **rim sag** = the mean vertical deviation of
the upper band (the part less than 1 mm deep; negative = too deep):

| strategy | vertical RMS | +DA1 | max | +DA1 | rim sag | +DA1 | upper band RMS | +DA1 | deeper part RMS (bias) | +DA1 | flange bias | Newton iterations | +DA1 |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| single-point (none) | **0.767** | 0.658 | 1.426 | 1.248 | **-1.151** | -0.992 | 1.164 | 1.011 | 0.367 (-0.19) | 0.295 (+0.01) | -0.252 | 1 906 | 1 037 |
| backing plate, 1 mm clearance | **0.370** | 0.297 | 0.758 | 0.626 | **-0.541** | -0.399 | 0.553 | 0.423 | 0.193 (-0.00) | 0.187 (+0.12) | -0.045 | 5 007 | 3 605 |
| DSIF, sine-law gap | **0.275** | 0.298 | 0.537 | 0.646 | **-0.221** | -0.053 | 0.258 | 0.161 | 0.284 (+0.24) | 0.355 (+0.34) | +0.043 | 3 282 | 3 680 |
| DSIF, gap t (sensitivity) | **0.357** | 0.318 | 0.745 | 0.669 | **-0.494** | -0.163 | 0.511 | 0.229 | 0.222 (+0.09) | 0.360 (+0.33) | -0.035 | 3 255 | 2 928 |
| DSIF, sine-law gap + rim pass | **0.361** | 0.258 | 0.692 | 0.597 | **-0.037** | -0.142 | 0.142 | 0.188 | 0.442 (+0.42) | 0.291 (+0.25) | +0.075 | 3 848 | 3 891 |

(Normal RMS within 0.03 mm of the vertical RMS throughout: `cases.csv`. The
cost is given in Newton iterations, which do not depend on the machine's
load; the single-thread runtimes of the first run, three at a time on an
otherwise idle 4-core machine, were 83-434 s - see
[Reproduce](#reproduce).)

**In one sentence:** on this model, support from below removes most of the
rim sag before any compensation - a backing plate halves the deviation
(0.767 to 0.370 mm RMS) and a DSIF support ball with the sine-law gap cuts
it by two thirds (0.275 mm), both better than single-point forming after a
DA step (0.658 mm) - but DSIF's lead over the plate comes with a squeeze of
the model's too-thick lower wall (with a gap of t it gives 0.357 mm), and
what is left moves from the rim to the deeper part, where one DA step helps
the plate (0.297 mm) and not DSIF.

## What the numbers say

* **The plate.** With the sheet outside the opening held at the plane, the
  rim sag halves (-1.15 to -0.54 mm) and the flange comes out almost flat
  (bias -0.25 to -0.05 mm); the deeper part is right on average. One DA
  step then brings the part to 0.297 mm, the best plate result and 2.2
  times better than single-point forming after the same step (0.658 mm).
  The rest of the sag sits in the 1.22 mm the opening leaves outside the
  outline on the 2 mm mesh (1 mm asked) plus the rim fillet (-0.63 mm at
  r = 7.5 mm, just inside the outline at 7.7 mm; -0.46 mm at r = 8 mm): the
  plate supports nothing inside its opening. It is also the most expensive
  run (5 007 Newton iterations against 1 906): 26-27 nodes are in contact
  with the plate and its reaction reaches 2.1-2.2 kN, more than the tool's
  1.3-1.5 kN - the plate's edge is the fulcrum about which the clamp holds
  the sheet down.
* **DSIF.** The support ball opposite the tool nearly removes the sag
  (-0.22 mm) and the flange (+0.04 mm), but the deeper part comes out
  0.24 mm too shallow (the floor 0.5 mm): the error moved rather than went.
  The sine law puts the balls t cos(theta) = 0.78 mm apart on the wall,
  while the model's wall is 0.94 mm thick (0.96 mm in single-point forming;
  the 2-layer mesh thins less than the sine law), so where the support sits
  opposite the tool it squeezes the wall. It sits there on a minority of
  the path only: wherever the opposite position would put the ball into the
  formed part's underside (seen from below the rim is a corner), the path
  keeps it lower (the lift bound of `dsif_support_points`). Of the 111
  contact knots of the path that is 81 %: all 45 of the first level (the
  tip less than 1.5 mm deep, where the rim forms; planned gap 1.7 mm, the
  ball 1.1 mm below its opposite position) and about two thirds of the
  wall and floor (0.2 mm below on the wall). The other 19 %, every one with
  the tip deeper than 2 mm, have the sine law's 0.78 mm gap, and the tool
  and support forces peak in that part of the path (1.7 and 1.6 kN, with
  the tip 2.2-2.4 mm deep; 1.8 times the single-point tool force). With the
  gap at the initial thickness t (sensitivity row) the support reaches its
  opposite position at 32 % of the knots, with a 1.00 mm gap, and otherwise
  plans nearly the same path (median planned gap 1.09 mm with the sine law,
  1.11 mm with t; 10th percentile 0.78 against 1.00 mm); the forces drop
  to 1.2 / 0.8 kN and the result lands near the plate's (0.357 against
  0.370 mm, sag -0.49 against -0.54 mm). So the sine-law DSIF's lead over
  the plate (0.275 against 0.370 mm) comes with that squeeze of the lower
  wall and floor, not with the support pressing harder at the rim; how a
  squeeze more than 2 mm deep lifts the rim this much is not established by
  one run each. The `squeeze` setting (0 here) does not show it. Without
  it, DSIF still halves single-point forming's deviation. One DA step moves
  the DSIF rim to within 0.05 mm but deepens the command by 0.54 mm and the
  part comes out *shallower* in the deeper region (bias +0.24 to +0.34 mm):
  the non-local response the springback benchmark saw on the deeper part is
  stronger here, and DSIF + DA1 (0.298 mm) is worse than DSIF alone. With
  the initial-thickness gap DA1 helps (0.357 to 0.318 mm).
* **The rim pass.** After forming, the support alone traces four loops from
  1 mm outside to 2 mm inside the outline, 323 N on average where it
  touches (29 % of its 203 increments, 568 N at most). It cuts the rim sag
  to -0.04 mm (upper band RMS 0.142 mm, the lowest) but lifts the whole
  part: the deeper part +0.42 mm too shallow, flange +0.08 mm, overall
  0.361 mm, worse than DSIF without it. With one DA step it gives the best
  part of the study, 0.258 mm (max 0.597 mm). (These two runs changed with
  the review's fixes: the loops now sit at a lift height exact over the
  grid's nodes rather than an interpolated one, slightly higher; the first
  run gave 0.338 / 0.267 mm, sag -0.07 / -0.12 mm.)
* **Commands above the sheet plane did not happen.** DA with the rim pass
  may raise the command above the plane only as far as the rim pass the
  deck writes for that command pushes the sheet (`command_upper_bound`: the
  top its ball sweeps along the loops, plus t, to a fixed point). DA asked
  for up to 0.31 mm on 790 nodes of the rim band and flange strip; seen
  from below the rim is a concave corner and the outermost loop is held
  down by the flange beside it, and the 4 mm ball reaches none of it.
  Computed on the same command, not simulated: a 3 or 2 mm ball reaches
  3-5 um on 2 nodes, a 1.5 mm ball 0.08 mm on 15, a 1 mm ball 0.16 mm on
  444. The DA1 command therefore stayed at z <= 0; the rim pass's gain came
  from pushing the sagged sheet back to the command, not from overbending
  it. (The first version of the bound took the reach of a ball placed
  anywhere under the command, not along the loops, and put a 2 mm ball at
  0.11 mm here: more than the written rim pass produces.)
* **Every run is physically sane** by the checks recorded
  (`cases.csv`, `tools.csv`): every step completed; the sheet pushes the
  plate and the support down and the tool up, never the other way; the
  largest penetration of any tool is 1.27 um; the plate carries forces only
  in "form" and the support only in "form" and "rim_pass"; the clamp holds
  0.24-0.94 kN in the unload steps (the only warning of every run - two with
  the rim pass - the springback the fixture prevents), and the 3-2-1
  release leaves at most 1.4e-12 N; 2-11 cut-backs per run, all recovered.
  In contact the support ball's body stays within 14.5 mm of the centre in
  plan, inside the 15 mm free window (no clash with a lower clamp frame).

Azimuthal means of the released parts [mm] (`profiles.csv`; the outline is
at r = 7.7 mm):

| r [mm] | target | none | plate | DSIF | DSIF gap t | DSIF + rim pass | none +DA1 | plate +DA1 | DSIF +DA1 | DSIF + rim +DA1 |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 0 | -3.09 | -2.88 | -2.71 | -2.55 | -2.65 | -2.40 | -2.80 | -2.65 | -2.53 | -2.54 |
| 2 | -3.09 | -3.09 | -3.05 | -2.82 | -2.91 | -2.66 | -3.00 | -2.99 | -2.75 | -2.74 |
| 4 | -2.65 | -2.60 | -2.48 | -2.31 | -2.42 | -2.13 | -2.36 | -2.39 | -2.26 | -2.35 |
| 5 | -1.85 | -2.22 | -1.97 | -1.71 | -1.90 | -1.52 | -1.98 | -1.81 | -1.58 | -1.72 |
| 6 | -1.05 | -1.82 | -1.37 | -1.03 | -1.30 | -0.83 | -1.62 | -1.16 | -0.84 | -0.97 |
| 7 | -0.25 | -1.49 | -0.86 | -0.54 | -0.81 | -0.36 | -1.34 | -0.72 | -0.38 | -0.46 |
| 7.5 | -0.02 | -1.35 | -0.65 | -0.35 | -0.62 | -0.18 | -1.21 | -0.55 | -0.20 | -0.26 |
| 8 | 0.00 | -1.21 | -0.46 | -0.20 | -0.45 | -0.04 | -1.09 | -0.39 | -0.05 | -0.09 |
| 9 | 0.00 | -0.98 | -0.23 | -0.04 | -0.24 | 0.11 | -0.89 | -0.20 | 0.12 | 0.10 |
| 10 | 0.00 | -0.79 | -0.10 | 0.05 | -0.11 | 0.17 | -0.72 | -0.09 | 0.22 | 0.21 |
| 12 | 0.00 | -0.48 | -0.04 | 0.08 | -0.04 | 0.15 | -0.44 | -0.04 | 0.21 | 0.21 |
| 15 | 0.00 | -0.18 | -0.03 | 0.07 | -0.01 | 0.09 | -0.17 | -0.03 | 0.16 | 0.14 |

## Setup

| | |
|-|-|
| Part | truncated cone, top radius 7.32 mm, wall 38.6 deg, depth 3.09 mm, fillets 1.08 (rim) and 2.08 mm (floor); the springback benchmark's test part truncated_cone-s2026-0000; 161 x 161 grid at 0.25 mm |
| Blank, material, mesh | 40 x 40 x 1 mm AA5754-O (nominal library data), 20 x 20 x 2 Hex8, mean dilatation, `finite_logarithmic` - the springback benchmark's model |
| Fixture, tool, path | clamp 5 mm wide at the blank edge; sphere of 4 mm radius, friction 0.1; spiral, 1 mm step-down, 1 mm point spacing, 1 mm travel per increment; form - unload (clamped) - release onto 3-2-1 supports |
| none | single-point forming: the benchmark's deck, byte for byte (the same RMS, 0.767 / 0.658 mm, as its uncompensated and FE-DA-1 rows) |
| backing plate | `support: "backing_plate"`, clearance 1 mm (1.22 mm to the outline realised on the mesh: the plate carries 320 of the 400 bottom faces), friction 0.1, removed with the tool in "unload" |
| DSIF | `support: "dsif"`: a 4 mm support ball, friction 0.1, squeeze 0, gap t cos(theta) (sine law) or t (`thickness_law: "initial"`), both balls on one pseudo-time, removed together in "unload" |
| rim pass | `rim_pass: true`: loops 1 mm outside to 2 mm inside the outline, 1 mm apart, at the lift-cutter height of the commanded underside (exact over the grid's nodes); the rim pass step and its removal follow "unload", clamp on |
| DA1 | `displacement_adjustment`, vertical error, alpha = 1, no smoothing, 65 deg wall limit, flange held, the fixture made for the target; with the rim pass the flange strip it sweeps is adjusted and the command may rise above the plane as far as the rim pass written for it reaches (`compensation_masks`, `command_upper_bound`) |

## Reproduce

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release && cmake --build build -j3
OPENBLAS_NUM_THREADS=1 python3 python/scripts/support_validation.py \
    --out benchmarks/support_cone --work benchmarks/support_cone/work --workers 3 --da-steps 1
```

The script reads the protocol from `config.json` (written from its
defaults when absent), runs every strategy's DA chain through precomp's
content-addressed cache in `work/runs` (git-ignored; an identical deck is
never run twice, so an interrupted run resumes) and writes `cases.csv` (one
row per case: deviation by region, cost, physics checks, wall thickness,
the plate's realised clearance), `tools.csv` (per case and tool: contact,
penetration, forces), `profiles.csv`, `tables.md` and `run.json`
(versions; the timing there is of the last invocation).

The tables are from the second run, after the review's fixes, on
30 September 2026 between 18:10 and 18:45 UTC with `sparlab 1.0.0
(bcd7e4c5003c)` (the same C++ source as the first run's 8ecc88656a15) and
precomp at 671209c. Eight of its ten decks are byte-identical to the first
run's (16:55-17:18 UTC) and gave identical results - every deviation,
increment, Newton iteration, cut-back, force and penetration; the
simulations are deterministic. The two with the rim pass changed with its
loops. The runtimes in `cases.csv` are of the second run, three at a time
on a 4-core machine that two simulations of another job shared (load
average 4-9): 1.05-2.6 times the first run's for the same decks, and not a
measure of the strategies' cost.

## Limitations

* **Simulation only**, with the springback benchmark's coarse model
  (two Hex8 layers, stiff in bending; 2 mm elements under a 4 mm ball, at
  most 3-6 nodes in contact) and nominal material data. The springback of a
  real part - and the effect of a real support - may differ in size and
  shape.
* **One part, one run per strategy.** Differences of a few hundredths of a
  millimetre are within what another part, mesh or step size could change;
  the ranking of plate against DSIF depends on the DSIF gap (sensitivity
  row) and on the plate's clearance, of which one value was run.
* **The DSIF gap follows the sine law, not the model's sheet.** The model's
  wall is about 20 % thicker than t cos(theta), so where the support
  reaches its opposite position - 19 % of the path's contact knots here,
  all with the tip deeper than 2 mm - the sine-law gap squeezes it; the
  sine-law DSIF's lead over the plate comes with that squeeze. The support's
  own compliance and a real machine's alignment errors are not modelled.
* **The plate's opening is quantised by the mesh** (1 mm asked, 1.22 mm to
  the outline realised) and the plate is flat at the sheet plane; nothing
  holds the sheet down onto it.
* **One DA step, alpha = 1.** DSIF + DA1 got worse: a relaxed or smoothed
  update, or a second step, may do better; not run.
* **Commands above the plane** were allowed only where the rim pass's ball
  reaches from below, which a 4 mm ball does nowhere on this rim; a smaller
  rim-pass ball would reach some (numbers above, computed), not simulated.
