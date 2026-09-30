# Springback pre-compensation: uncompensated vs FE displacement adjustment vs ML first shot

Every formed part in this benchmark is a **SparLab simulation** (`sparlab_form`:
implicit, quasi-static, `finite_logarithmic` kinematics, Hill48 plasticity with
**nominal handbook data** for AA5754-O). Nothing here was formed or measured on
a machine. The numbers say how the three strategies compare *on this model*;
they are not a claim about real parts (see [Limitations](#limitations)).

The question: given a target part, how close does the released part come to
it when you

1. **form the target as it is** (1 FE run),
2. run **displacement adjustment with the FE model in the loop** (FE-DA-k: the
   part formed after k compensation steps, k + 1 FE runs), or
3. let a **surrogate trained on simulations** compensate the target by
   displacement adjustment on its own prediction, then form the result once
   (ML first shot: 1 FE run per part, plus the shared training set)?

## Headline

Final stage (n18): 48 training and 16 calibration parts (144 simulated
samples), the **same 8 held-out test parts** (2 per family) at every stage.
Deviation of the released part from the target over the part, mean over the
8 test parts [mm]:

| Method | FE runs per part | vertical RMS | median | vertical max | normal RMS | RMS / uncompensated | parts better than uncompensated | parts better than FE-DA-1 |
|---|--:|--:|--:|--:|--:|--:|--:|--:|
| uncompensated (target as commanded) | 1 | **0.731** | 0.735 | 1.341 | 0.729 | 1.00 | - | - |
| FE-DA-1 (FE in the loop, one step) | 2 | **0.623** | 0.620 | 1.129 | 0.629 | 0.85 | 8 / 8 | - |
| FE-DA-2 (FE in the loop, two steps) | 3 | **0.602** | 0.618 | 1.156 | 0.608 | 0.82 | 8 / 8 | 7 / 8 |
| ML first shot, GBM ensemble | 1 (+ training set) | **0.634** | 0.634 | 1.147 | 0.640 | 0.87 | 7 / 8 | 4 / 8 |
| ML first shot, MLP ensemble | 1 (+ training set) | **0.626** | 0.631 | 1.165 | 0.631 | 0.86 | 7 / 8 | 5 / 8 |

The learning curve over the stages (same test parts; mean vertical RMS [mm];
dz error = the model's RMS error of the predicted deviation on the 16 test
samples):

| Stage | train / calibration parts | ML-GBM | ML-MLP | GBM dz error | MLP dz error | FE-DA-1 | FE-DA-2 | uncompensated |
|---|---|--:|--:|--:|--:|--:|--:|--:|
| n6 (pilot) | 12 / 4 | 0.627 | - | 0.120 (17.7 %) | - | 0.623 | - | 0.731 |
| n10 | 24 / 8 | 0.615 | 0.613 | 0.098 (14.0 %) | 0.109 (15.7 %) | 0.623 | 0.602 | 0.731 |
| n14 | 36 / 12 | 0.621 | **0.605** | 0.092 (13.2 %) | 0.108 (15.3 %) | 0.623 | 0.602 | 0.731 |
| n18 | 48 / 16 | 0.634 | 0.626 | 0.085 (12.3 %) | 0.106 (14.8 %) | 0.623 | 0.602 | 0.731 |

**In one sentence:** on this model a surrogate trained on 12-48 simulated
parts compensates a new part in one shot about as well as one step of
FE-in-the-loop displacement adjustment (which costs a second simulation of
that part), and a second FE step is better still - but every method removes
only 13-18 % of the deviation, because 79-87 % of the remaining squared
error is a rim sag that no command at z <= 0 can correct.

## Design

| | |
|-|-|
| Parts | four families - truncated cone, pyramid (rounded rectangle), dome (spherical / ellipsoidal cap), elliptic cone - drawn per family from a scrambled Sobol sequence (`precomp.ml.design_points`, seed 2026) over footprint 14-22 mm, depth 2-4 mm, walls 30-55 deg, fillets 1-4 mm (`config.json`, `part_bounds`); 0.25 mm grid |
| Blank, material | 40 x 40 x 1 mm AA5754-O from the precomp library: E 70 GPa, Hollomon K = 420 MPa, n = 0.30 from 100 MPa fitted by linear + Voce hardening, Hill48 r0 / r45 / r90 = 0.75 / 0.70 / 0.80 - nominal, not certified |
| Mesh, solver | 20 x 20 x 2 Hex8 (2 mm in plane, 0.5 mm through the thickness), mean dilatation; `finite_logarithmic`; implicit Newton with cut-backs - the `fea_da_cone` model |
| Fixture, tool | clamp 5 mm wide at the blank edge (a 30 mm window); 4 mm radius sphere, friction 0.1; spiral path, 1 mm step-down, 1 mm point spacing, 1 mm tool travel per increment; form - unload (clamped) - release onto 3-2-1 supports |
| Split | by design-point index within a family, fixed before the first run and stable as the design grows: indices 0-1 **test** (8 parts), then every 4th **calibration**, the rest **train** |
| Data | per train and calibration part two simulated commanded variants: the target, and one FE displacement-adjustment step from its simulated part (the shape the surrogate is queried at during compensation); test parts the same two, for the dz evaluation only |
| FE-DA | `precomp.compensation.displacement_adjustment` with a `SimulatorPredictor` around the same simulator: c_(k+1) = c_k - (f(c_k) - t), vertical error, alpha = 1, no smoothing, flange held, commanded surface kept at z <= 0 and within 65 deg - exactly `precomp compensate --method fea` |
| ML | `precomp.ml.train_surrogate` with `GBMEnsemble` (8 x 300 iterations) or `MLPEnsemble` (5 x (128, 128), 40 epochs), 57 point features, 2000 stratified nodes per sample, calibrated (conformal) on the calibration parts, training envelope enforced; `surrogate_compensate` (at most 8 predictions, stop at < 2 % predicted improvement or when an iterate leaves the envelope, keep the best predicted iterate); then one `sparlab_form` run of the compensated shape - the part the method delivers |
| Metric | vertical deviation of the released part (after the 3-2-1 release) from the target over the part (target below the sheet plane): RMS and max over the nodes, and the RMS of the normal deviation; regions: the upper band of the part (< 1 mm deep), the deeper part, the flange |

The split never changes a part's role: stage n6 (6 design points per family)
trains on 12 parts, n10 on 24, n14 on 36, n18 on 48, always tested on the
same 8 parts,
so the stages form a learning curve on a fixed test set.

## Results by stage

Per test part, vertical RMS [mm] (largest |deviation| in brackets), stage
n18; `*` = the target was outside the model's training envelope and was
compensated with the override:

| part | uncompensated | FE-DA-1 | FE-DA-2 | ML-GBM | ML-MLP |
|---|--:|--:|--:|--:|--:|
| dome-s2026-0000 | 0.654 (1.07) | 0.532 (0.90) | 0.485 (0.83) | 0.654 (1.07) | 0.654 (1.07) |
| dome-s2026-0001 | 0.776 (1.40) | 0.625 (1.25) | 0.555 (1.17) | 0.592 (1.22) | 0.616 (1.25) |
| elliptic_cone-s2026-0000 | 0.621 (1.21) | 0.512 (1.06) | 0.493 (1.01) | 0.513 (1.06) * | 0.495 (1.03) * |
| elliptic_cone-s2026-0001 | 0.841 (1.51) | 0.726 (1.30) | 0.714 (1.29) | 0.733 (1.32) | 0.728 (1.31) |
| pyramid-s2026-0000 | 0.791 (1.45) | 0.615 (1.19) | 0.604 (1.17) | 0.601 (1.16) | 0.593 (1.15) |
| pyramid-s2026-0001 | 0.699 (1.38) | 0.615 (1.09) | 0.648 (1.33) | 0.614 (1.09) | 0.614 (1.09) |
| truncated_cone-s2026-0000 | 0.767 (1.43) | 0.658 (1.25) | 0.632 (1.20) | 0.670 (1.26) | 0.659 (1.25) |
| truncated_cone-s2026-0001 | 0.703 (1.29) | 0.700 (0.98) | 0.686 (1.25) | 0.699 (0.98) | 0.647 (1.19) |

The other stages are in `stage_n6/`, `stage_n10/`, `stage_n14/` (`tables.md`).

Paired differences, ML minus FE-DA on the same part (vertical RMS, mm; mean
+- standard deviation over the 8 parts, and on how many parts ML is lower):

| Stage | GBM - DA-1 | MLP - DA-1 | GBM - DA-2 | MLP - DA-2 |
|---|--:|--:|--:|--:|
| n6 | +0.004 +- 0.051 (5/8) | - | - | - |
| n10 | -0.008 +- 0.018 (4/8) | -0.010 +- 0.018 (5/8) | +0.013 +- 0.027 (3/8) | +0.011 +- 0.021 (2/8) |
| n14 | -0.002 +- 0.017 (3/8) | -0.018 +- 0.018 (7/8) | +0.019 +- 0.026 (1/8) | +0.003 +- 0.028 (4/8) |
| n18 | +0.012 +- 0.047 (4/8) | +0.003 +- 0.051 (5/8) | +0.032 +- 0.060 (2/8) | +0.024 +- 0.067 (3/8) |

With 8 test parts, differences of 0.01-0.02 mm between the ML shot and
FE-DA-1 are within the part-to-part scatter; the n18 figures are dominated
by one part (below).

**The surrogate as a model.** The GBM's dz error on the test parts fell from
0.120 mm (12 training parts) to 0.085 mm (48) against a dz of 0.70 mm RMS; it
is 1.1-1.5 times larger on the compensated variants (0.098-0.127 mm, 15-20 %)
than on the targets (0.073-0.113 mm, 10-15 %) - the shapes the model is
queried at during compensation are the harder ones. The MLP predicts dz less well
(0.106-0.109 mm) yet gave the better compensations at n10-n18. With 12 and 16
calibration parts the 90 % conformal intervals are supported: expected
per-part coverage 0.96-1.00 over the test samples, the worst part 0.74
(GBM, n18) - as documented, the guarantee is on average over parts, not per
part. The surrogate's *predicted* residual after its own compensation was
within 0.045 mm of the verified one on every test part at n14 and n18
(0.056 mm at n10, 0.077 mm at the pilot; e.g. 0.697 predicted, 0.670
verified): the model correctly foresaw that
compensation cannot get much below 0.6 mm here.

**The envelope costs a part at n18.** Enforcing the training envelope is a
safety feature, and at n18 it decided dome-s2026-0000: the first compensated
iterate of both models fell outside the (re-fitted) envelope, so
`surrogate_compensate` kept the best iterate inside - the target itself - and
the "ML" part is the uncompensated one (0.654 mm; its verification run is the
cached uncompensated run). At n10 and n14 the same part was compensated to
0.49-0.54 mm. Without that part the n18 means are GBM 0.632, MLP 0.622,
FE-DA-1 0.636, FE-DA-2 0.619 mm (7 parts). The envelope flagged a quarter of
the test samples at n18 against 0-12 % before: with more parts its 99 %
thresholds, taken over held-out training parts, came out tighter. One target
(elliptic_cone-s2026-0000) was outside the envelope at n14 and n18 and was
compensated with the override (flagged `*`; it was compensated well).

## Where the error is

`regions.csv`, stage n18, mean over the test parts [mm]; the upper band is
the part less than 1 mm deep (the rim fillet and the top of the wall), the
share is that band's part of the squared error over the part:

| Method | upper band RMS | bias | deeper part RMS | bias | flange RMS | share of error in upper band |
|---|--:|--:|--:|--:|--:|--:|
| uncompensated | 1.099 | -1.086 | 0.328 | -0.137 | 0.376 | 0.87 |
| FE-DA-1 | 0.879 | -0.856 | 0.346 | +0.126 | 0.324 | 0.79 |
| FE-DA-2 | 0.905 | -0.884 | 0.272 | +0.021 | 0.329 | 0.87 |
| ML-GBM | 0.896 | -0.874 | 0.353 | +0.115 | 0.328 | 0.79 |
| ML-MLP | 0.908 | -0.886 | 0.330 | +0.054 | 0.329 | 0.82 |

Near the rim every part comes out about 0.9-1.1 mm **too deep**: between the
5 mm clamp and the first contour the sheet bends down (the `fea_da_cone`
record measured 1.3 mm at the rim of a 45 deg cone with this clamp). Vertical
DA would have to command that band *above* the sheet plane; a tool pressing
from above cannot, and the compensation keeps the command at z <= 0 with the
flange held, so the band keeps 80-90 % of its error whatever the method. The
methods differ in the deeper part: one FE-DA step over-corrects it (bias
-0.14 to +0.13 mm, RMS 0.328 to 0.346), a second step brings it to 0.272 mm
with no bias, and the MLP's shot lands in between (0.330, +0.05). The
ranking of the methods is decided there, on 13-21 % of the squared error.

## Cost

| | |
|-|-|
| One simulation | 243 s mean, 227 s median, 103-609 s wall on one thread (166 runs of stage n18); 223 increments, 1 517 Newton iterations and 5 cut-backs on average; largest plastic strain 0.16-0.59. Two ran at a time on a 4-core machine shared with other builds and solver runs (load average 2-9): the first two timing runs, at load 9, took 627 and 637 s, the `fea_da_cone` cone about half that on a quiet machine |
| All simulations | 207 distinct `sparlab_form` runs in the cache, 13.8 h of single-thread wall time, 8.3 h elapsed at two at a time (01:06-09:24 UTC, all four stages) |
| Per stage (elapsed) | n6 2 h 28 min (56 runs), n10 1 h 51 min (32 data + 8 DA-2 + 16 verifying), n14 1 h 44 min (32 + 16), n18 1 h 53 min (32 + 16) |
| Training | GBM 129-324 s, MLP 52-123 s on two threads (12-48 parts, 48 000-192 000 rows); evaluation and surrogate DA of 8 targets 1-6 min |
| Per test part | FE-DA-k: k + 1 simulations of that part (about 4 min each here). ML: one simulation of that part after about 10 s of surrogate DA; the training set (2 simulations per training or calibration part: 128 at n18) is paid once for all parts of the design space |

## Failures

**None of the 207 simulations failed** (`failures.csv` of every stage is
empty; no generation, FE-DA or verification job failed). Within the runs
the solver cut the increment back 1 079 times (5 per run on average), 1 077
of them because an element inverted under the tool ("the load step is too
large or the mesh too coarse") - every one recovered by halving the step.
Every run carries one warning, expected with this deck: the clamp of the
"unload" step holds reactions (0.7-1.1 kN), the springback the clamp prevents,
which the 3-2-1 release lets go. What did not work as intended is the ML
side, not the solver: at n18 the envelope stopped both surrogates on one test
part before any compensation (above), and at the pilot the GBM stopped there
for the same reason; one target (elliptic_cone-s2026-0000) needed the
envelope override at n14 and n18.

## What the numbers say

* **The ML first shot is as good as one FE-in-the-loop step, at half the
  per-part simulation cost.** Averaged over the stages the surrogate's part
  (one simulation) is within +-0.02 mm of FE-DA-1's (two simulations), and
  better than the uncompensated part in 53 of 56 part-stage-model cases (the
  three exceptions: the envelope kept the target, dome-s2026-0000). That is the whole promise
  of a surrogate on this model: it replaces the first trial.
* **It is not better than FE-DA with a second step**, which ends lower on
  average (0.602 mm) and wins most paired comparisons against both models -
  though it made one part worse than the first step (pyramid-s2026-0001,
  0.615 to 0.648 mm, an over-correction of the floor), which the surrogates
  did not.
* **More data improved the prediction, not the compensation.** Four times as
  many training parts cut the GBM's dz error by 30 % (0.120 to 0.085 mm), but
  the verified deviation stayed at 0.60-0.63 mm: the compensation is limited
  by the rim sag, which no command can fix, and by the non-local, not quite
  monotone response of the deeper part, not by the model's accuracy. The
  best stage (n14, MLP 0.605 mm) is within noise of the others.
* **What limits every method is the process, not the algorithm.** 79-87 %
  of the error after compensation is the rim sag between clamp and part. A
  backing plate close to the part (the `fea_da_cone` record: 37 % less
  uncompensated deviation), a registration of the released part on its
  flange, or a tool path that works the rim would each do more than any
  compensation here. Displacement adjustment restricted to z <= 0 is the
  wrong tool for a region that comes out too deep next to the sheet plane.
* **The surrogate knew.** Its predicted residual matched the verified one
  within 0.045 mm on every part from n14 on, so a user would have been told before
  forming that compensation buys little here - arguably the most useful
  output of the ML model on this process.

## Reproduce

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release && cmake --build build -j2
pip install -e '.[torch]'
W=benchmarks/springback/work          # run cache, data set, models (git-ignored, ~300 MB)
python3 python/scripts/springback_benchmark.py --out benchmarks/springback --work $W \
    --n-per-family 6  --models gbm     --workers 2 --da-fe-runs 2    # pilot
python3 python/scripts/springback_benchmark.py --out benchmarks/springback --work $W \
    --n-per-family 10 --models gbm mlp --workers 2 --da-fe-runs 3
python3 python/scripts/springback_benchmark.py --out benchmarks/springback --work $W \
    --n-per-family 14 --models gbm mlp --workers 2 --da-fe-runs 3
python3 python/scripts/springback_benchmark.py --out benchmarks/springback --work $W \
    --n-per-family 18 --models gbm mlp --workers 2 --da-fe-runs 3
```

Each command resumes: an identical deck is fetched from the content-addressed
run cache (`<work>/runs`), the data set skips the samples it has, model
bundles (`<work>/models/<stage>`) are reused when their training and
calibration samples are the same, and a surrogate-compensated shape
(`<work>/ml/<stage>`) is kept. A larger `--n-per-family` adds design points
to the same data set - the first n Sobol points of a family do not depend on
n - so the earlier stages' runs are reused. A rerun of a stage from the cache
regenerates its tables and keeps the first run's timings in `run.json`
(`reruns` lists the others). Simulations are single-threaded
(`threads: 1`); `--workers` of them run at a time; training uses
`--threads` (2).

Versions: precomp 0.1.0, `sparlab 1.0.0 (a014858e7dbc-dirty)` (the
Release build in `build/`; `-dirty` as the executable reports it: built from
a source tree with uncommitted changes), Python 3.11.15, numpy 2.4.6, scipy 1.17.1,
pandas 3.0.6, scikit-learn 1.9.1, torch 2.14.0+cpu; Linux, 4 cores shared
with other jobs (load average 4-9 during the runs). Each stage's `run.json`
records the command, versions, git revision, timings and load averages.

## Files

`config.json` is the protocol (design, solver block, split, model settings).
Per stage, `stage_n<N>/`:

| File | Content |
|------|---------|
| `tables.md` | the tables of this README for that stage |
| `headline.csv` | per test part and method: FE runs, vertical RMS / max / bias, normal RMS / max, notes (envelope override, failure) |
| `summary.csv` | per method: means over the test parts, ratio to uncompensated, parts better than uncompensated and than FE-DA-1 |
| `regions.csv` | per test part and method: RMS and bias in the upper band, the deeper part and the flange |
| `ml_compensation_<model>.csv` | the surrogate compensation: why it stopped, predictions used, predicted residual, envelope verdicts, the verifying run |
| `dz_per_part.csv` | the models' dz prediction error on the test samples (`precomp.ml.evaluate_surrogate`) |
| `design.csv` | every design point: family, split, footprint, depth, wall angle, parameters |
| `simulations.csv` | every run of the stage: runtime, increments, Newton iterations, cut-backs, warnings, largest plastic strain, status and failure reason |
| `failures.csv` | failed generation, FE-DA and verification jobs with their reasons |
| `run.json` | command, versions, timings, load averages, counts |

## Swapping the solver

The script knows the solver only through `config["solver"]`: a simulator name
(`"sparlab_form"` -> `precomp.ml.SparlabSimulator`) and the `FormingSetup`
fields. Data generation, FE-DA (through `SimulatorPredictor`) and the
verification runs all call the `precomp.ml.Simulator` protocol
(`run([(setup, commanded)]) -> [SimOutcome]`). Explicit forming (a
`form_explicit` step type) therefore changes that block - the setup fields
that select the step type and its settings, once `FormingSetup` and
`precomp.fea.deck` write them - and nothing else; the run cache keys on the
deck, so implicit and explicit runs never mix. Use a new `--work` (or at
least a new data set) for another solver: a model trained on one solver's
runs is not a model of the other's, and the training envelope records the
setup fields it saw.

## Limitations

* **Simulation only.** No part was formed or measured; the material data are
  nominal handbook values (Hollomon fitted by linear + Voce, Hill48 with
  r < 1, which underestimates the equibiaxial yield stress); the forming
  model is coarse: two Hex8 layers (stiff in bending), 2 mm elements under a
  4 mm tool (at most a few nodes in contact; frequent inverted-element
  cut-backs), one tool path style, one fixture, no thinning check. The
  springback of a real part may differ in size and shape.
* **Small test set.** 8 test parts, 2 per family; 7 of 8 of every
  comparison is one part away from 6 of 8. The paired differences between
  ML and FE-DA-1 (0.01-0.02 mm) are within the part-to-part scatter. The
  test parts are the first two Sobol points of each family, not a random
  draw.
* **One process point.** One material, thickness, tool, step-down,
  friction, clamp and blank size; the models were trained on exactly this
  setup and their envelope refuses others. Nothing here says how the
  methods compare on deeper parts or steeper walls (these parts are 2-4 mm
  deep, 30-55 deg), where DA contracts differently.
* **Two commanded variants per training part** (the target and one FE-DA
  step) to fit the simulation budget; the perturbed variant of
  `precomp.ml.generate` was not simulated (`train_kinds` in `config.json`
  adds it to the same data set).
* **FE-DA as configured**: vertical error, alpha = 1, no smoothing, no
  registration of the released part, z <= 0, flange held - the defaults of
  `precomp compensate --method fea`. Relaxation, smoothing or normal-direction
  DA may converge faster; three simulations is not convergence.
* **Envelope sensitivity.** The envelope's thresholds are 99 % quantiles over
  a few dozen held-out training parts; re-fitting them per stage changed
  which iterates are accepted, which alone decided one of eight parts at
  n18. The override was used only where the *target* was outside, never to
  rescue a stopped iteration.
* **Machine.** Runtimes were measured on a shared 4-core machine at load
  2-9 and vary by tens of percent between runs; the simulations themselves
  are deterministic (a rerun reproduces the same result from the same deck).
