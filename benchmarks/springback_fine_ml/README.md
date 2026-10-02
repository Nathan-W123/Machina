# ML first shot on the fine mesh with a backing plate

Every formed part here is a **SparLab simulation** (`sparlab_form`: implicit,
quasi-static, `finite_logarithmic` kinematics, Hill48 plasticity with
**nominal handbook data** for AA5754-O). Nothing was formed or measured on a
machine. The numbers compare compensation methods *on this model*. They are
not a claim about real parts (see [Limitations](#limitations)).

This is stage 2 of the plan in `benchmarks/springback_fine/README.md`
(section "Recommendation for the ML dataset stage"): train the surrogates of
`benchmarks/springback` (GBM and MLP ensembles) on fine-mesh backing-plate
simulations, let them compensate the same 8 held-out test parts, and form
each compensated shape **once** on the converged mesh. The baselines are the
fine benchmark's backing-plate runs of the same parts.

## Headline

Stage n18: 48 training and 16 calibration parts, 128 simulated samples at
contact penalty 3, no generation failures. The 8 test parts (2 per family)
are verified at penalty 10, the setting of the baselines. Deviation of the
released part from the target over the part, mean over the 8 test parts
[mm]. **Rim sag** is the mean vertical deviation of the part less than 1 mm
deep (negative means too deep). **Interior** is the part deeper than that.

| method | FE runs / part | vertical RMS | median | max | rim sag | interior RMS (bias) | RMS / uncomp. | parts better than uncomp. | parts better than FE-DA-1 | minus FE-DA-1, mean (sd) |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| uncompensated | 1 | **0.300** | 0.298 | 0.583 | -0.436 | 0.149 (-0.05) | 1.00 | - | 1 / 8 | +0.037 (0.043) |
| FE-DA-1 (FE in the loop, one step) | 2 | **0.263** | 0.252 | 0.535 | -0.334 | 0.162 (+0.08) | 0.87 | 7 / 8 | - | - |
| FE-DA-2 (two steps) | 3 | **0.272** | 0.259 | 0.595 | -0.318 | 0.188 (+0.12) | 0.90 | 6 / 8 | 3 / 8 | +0.009 (0.020) |
| FE-DA-1, fixed decks | 2 | **0.264** | 0.256 | 0.538 | -0.337 | 0.163 (+0.09) | 0.88 | 7 / 8 | 1 / 8 | +0.001 (0.004) |
| ML first shot, GBM ensemble | 1 (+ training set) | **0.253** | 0.242 | 0.521 | -0.348 | 0.126 (+0.06) | 0.84 | 7 / 8, 1 tie | 5 / 8 | -0.010 (0.018) |
| ML first shot, MLP ensemble | 1 (+ training set) | **0.250** | 0.246 | 0.519 | -0.350 | 0.116 (+0.04) | 0.83 | 7 / 8, 1 tie | 6 / 8 | -0.012 (0.017) |

"FE-DA-1, fixed decks" (`fe_da_fixed/`) repeats FE-DA-1 with the two deck
fixes made after the fine benchmark (the spiral's last revolution walked to
its loop, d522474, and the DA command clipped to the tool's reach, 459edeb).
Three of the eight decks changed. The mean moved from 0.263 to 0.264 mm, so
the fixes do not change the baseline. The ML data and the ML shots use the
fixed decks.

**The tie.** On truncated_cone-s2026-0001 both surrogates predicted that no
compensation helps and kept the target as the command. The verifying deck is
then the uncompensated deck (same hash, a cache hit), and the result is the
uncompensated 0.344 mm. `summary.csv` and `tables.md` count it as "better
than uncompensated" (8 / 8) because of a rounding difference in the stored
baseline. It is a tie, and the table above says so.

**In one sentence:** on the plate, a surrogate trained on 128 fine-mesh
simulations gives a first shot of 0.253 (GBM) / 0.250 mm (MLP) with one FE
run per part, as good as one step of FE displacement adjustment (0.263 mm,
two FE runs per part) and slightly better on average, but the 0.010-0.012 mm
margin is within the part-to-part spread of 8 test parts and is not
established.

## What the numbers say

* **ML matches FE-DA-1 with one run fewer per part.** Paired over the 8
  parts, ML minus FE-DA-1 is -0.010 +- 0.018 mm (GBM) and -0.012 +-
  0.017 mm (MLP); paired t-test p = 0.17 and 0.08, Wilcoxon p = 0.25 and
  0.08. Against the fixed-deck FE-DA-1, MLP gives p = 0.04 (t) and 0.02
  (Wilcoxon), GBM 0.10 and 0.20. Those are four of eight uncorrected tests
  on 8 parts: read them as "no worse than FE-DA-1", not as "better". Both
  are clearly better than forming the target as it is (-0.047 and
  -0.050 mm; paired t p <= 0.01, Wilcoxon p = 0.02 over the 7 untied parts).
* **Where the gain comes from.** The interior: 0.126 / 0.116 mm RMS against
  0.162 mm after FE-DA-1, with a smaller positive bias (+0.06 / +0.04
  against +0.08 mm). The rim sag is slightly worse than FE-DA-1's (-0.348 /
  -0.350 against -0.334 mm). As in the fine benchmark, the rim sag dominates
  the remaining deviation, and no method here removes it.
* **The surrogate avoids the DA overshoot on truncated_cone-s2026-0001.**
  This is the part where FE-DA makes the plate result worse (0.386 / 0.426
  against 0.344 mm, see `benchmarks/springback_fine`). Both surrogates
  predicted no improvement from any DA step and kept the target, so they
  deliver the uncompensated 0.344 mm. That alone is 0.005 mm of the
  0.010-0.012 mm mean margin over FE-DA-1. On the other seven parts ML and
  FE-DA-1 are within 0.03 mm of each other, either way, most within
  0.01 mm; the largest ML gains are truncated_cone-s2026-0000 (GBM 0.221
  against 0.251 mm) and pyramid-s2026-0001 (MLP 0.269 against 0.299 mm).
* **Predicted against achieved.** The verified RMS (penalty 10) is within
  -0.08 to +0.04 mm of the surrogate's predicted residual (penalty 3). MLP
  is optimistic on 6 of 8 parts (by 0.004-0.035 mm), GBM on 4. Both are
  pessimistic by about 0.08 mm on elliptic_cone-s2026-0001, where the shot
  stopped on leaving the envelope. The prediction does not rank the parts
  reliably: only the verifying run says how good a shot is.
* **dz accuracy.** On the 16 penalty-10 test samples (uncompensated and
  fixed FE-DA-1 of each test part) the models predict the formed deviation
  with an RMS error of 0.049 mm (GBM) and 0.052 mm (MLP), 17-18 % of the
  0.289 mm deviation itself. The conformal 90 % intervals cover 95 % (GBM)
  and 88 % (MLP) of the test nodes on average; with 8 parts the MLP's
  88 % is not evidence of under-coverage. These errors include the penalty
  mismatch (trained at 3, tested at 10), which the fine benchmark measured
  at +0.000-0.009 mm RMS on the plate.
* **The envelope.** One test target, elliptic_cone-s2026-0000, is outside
  the training envelope (1 of 16 test samples flagged). It was compensated
  with the override and is marked in `tables.md`. Its result (0.263 /
  0.261 mm) is worse than FE-DA-1 (0.252 mm) but equal to the fixed-deck
  FE-DA-1 (0.261 mm).

Per part, vertical RMS [mm] (max |dev| in brackets):

| part | uncompensated | FE-DA-1 | FE-DA-2 | FE-DA-1 fixed | ML-GBM | ML-MLP |
|---|--:|--:|--:|--:|--:|--:|
| dome-s2026-0000 | 0.307 (0.50) | 0.213 (0.41) | 0.203 (0.40) | 0.213 (0.41) | 0.214 (0.41) | 0.217 (0.41) |
| dome-s2026-0001 | 0.290 (0.59) | 0.223 (0.51) | 0.216 (0.50) | 0.223 (0.51) | 0.222 (0.52) | 0.215 (0.52) |
| elliptic_cone-s2026-0000 | 0.277 (0.57) | 0.252 (0.50) | 0.294 (0.69) | 0.261 (0.53) | 0.263 (0.53) * | 0.261 (0.53) * |
| elliptic_cone-s2026-0001 | 0.305 (0.61) | 0.266 (0.56) | 0.260 (0.55) | 0.270 (0.55) | 0.262 (0.54) | 0.265 (0.56) |
| pyramid-s2026-0000 | 0.290 (0.57) | 0.213 (0.48) | 0.217 (0.48) | 0.210 (0.48) | 0.197 (0.48) | 0.202 (0.49) |
| pyramid-s2026-0001 | 0.308 (0.69) | 0.299 (0.61) | 0.304 (0.63) | 0.299 (0.61) | 0.301 (0.62) | 0.269 (0.56) |
| truncated_cone-s2026-0000 | 0.280 (0.55) | 0.251 (0.46) | 0.257 (0.61) | 0.251 (0.46) | 0.221 (0.47) | 0.231 (0.49) |
| truncated_cone-s2026-0001 | 0.344 (0.60) | 0.386 (0.74) | 0.426 (0.89) | 0.386 (0.74) | 0.344 (0.60) = | 0.344 (0.60) = |

`*` target outside the training envelope (compensated with the override);
`=` the surrogate kept the target, so the shot is the uncompensated run.

## Against the 2 mm benchmark

In `benchmarks/springback` (2 mm mesh, single-point forming, same design,
same split, stage n18) the ML first shot was slightly *worse* than FE-DA-1:
0.634 (GBM) / 0.626 (MLP) against 0.623 mm, 4 and 5 of 8 parts better. Here
it is slightly better: 5 and 6 of 8 parts. The dz error is smaller in
absolute terms (0.049 / 0.052 mm against 0.085 / 0.106 mm) but larger
relative to the deviation (17-18 % against 12-15 %), because the plate
halves the deviation. The difference is most likely the setting,
not the model: with the plate the deviation is half as large and FE-DA has less to
gain (0.300 to 0.263 mm), FE-DA-2 is worse than FE-DA-1, and on one
part FE-DA goes the wrong way while the surrogate stops.

## Cost

| | runs | runtime per run, mean | CPU-h |
|---|--:|--:|--:|
| training and calibration data (64 parts x 2, penalty 3) | 128 | 1 222 s | 43.4 |
| training (GBM / MLP, 4 threads) | - | 110 s / 57 s | - |
| surrogate compensation per part (GBM / MLP) | - | 3.6 s / 2.1 s | - |
| ML verifying run per test part (penalty 10) | 1 | 2 670-2 750 s (one cache hit) | 0.75 |
| FE-DA-1 per part (the fine benchmark, penalty 10) | 2 | 3 430 s | 1.94 |

Run 4 at a time on a 4-core machine, so runtimes are under a load average of
about 4. FE-DA-1 costs about one extra penalty-10 run (about 1 CPU-h) per
part. The training set cost about 43 CPU-h, so ML pays for its data only
after about 45 parts of this kind, and only if its accuracy holds on them.

## Failures

None. No data part failed (`data_status.csv`: 64 parts, 128 runs, no
retries), no verifying run failed (`stage_n18/failures.csv` is empty), and
no baseline run failed. The driver was restarted several times after
container restarts. Each restart resumed from the run cache and the data
set, and only the runs in flight were repeated.

## Stage n26

The driver continues with stage n26 (26 design points per family: 32 more
data parts, 64 runs, then both models and their shots again on the same 8
test parts). It started at 23:56 UTC on 2 October. At the stage n18 pace
(about 35 min per part, 4 at a time) the data takes about 5 h and the ML
shots about 3 h. Its results will land in `stage_n26/` and this README
will need updating then.

## Setup

| | |
|-|-|
| Design | `benchmarks/springback`'s: four families (truncated cone, pyramid, dome, elliptic cone), its part bounds, Sobol seed 2026, 0.25 mm grid. Split by part: indices 0-1 of each family are the test parts (the fine benchmark's 8), then every 4th calibration, the rest training. Stage n<N> holds the first N points per family |
| Mesh, solver | 48 x 48 x 1 Hex8 with incompatible modes (0.833 mm in plane), 5 points through the thickness, `finite_logarithmic`; data at `springback_fine_bulk` (penalty 3), test runs at `springback_fine` (penalty 10) |
| Process | 40 x 40 x 1 mm AA5754-O (nominal data); 4 mm sphere, spiral path; backing plate with 1 mm clearance to the target's outline |
| Data per part | the target formed as it is, and one FE-DA step from that result (vertical, alpha 1, no smoothing, 65 deg wall limit, the plate's masks and bound, the command clipped to the tool's reach), both formed |
| Models | `GBMEnsemble` (8 members, 300 iterations) and `MLPEnsemble` (5 members, 128 x 128, 40 epochs), 57 point features, 2000 nodes per sample, conformal calibration on the calibration parts, training envelope enforced |
| ML first shot | `surrogate_compensate` with the training setup (<= 8 predictions, stop at < 2 % predicted improvement, keep the best predicted iterate, reach clip), then one simulation of the compensated shape at penalty 10 |
| Metric | vertical deviation of the released part from the target over the part: RMS and max; upper band (< 1 mm deep), interior, flange |

Files in `stage_n18/`: `headline.csv` (per part and method), `summary.csv`,
`tables.md`, `ml_compensation_<model>.csv` (iterations, predicted and
verified RMS, envelope), `dz_per_part.csv`, `design.csv`,
`simulations.csv` (every run of the stage with its deck hash),
`failures.csv`, `run.json` (versions, timings, dataset counts).

## Reproduce

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release && cmake --build build -j3
benchmarks/springback_fine_ml/run.sh          # stages 18 and 26, background, resumable
# the fixed-deck FE-DA-1 baseline (fe_da_fixed/):
OPENBLAS_NUM_THREADS=1 PYTHONPATH=python python3 python/scripts/springback_fine_benchmark.py \
    --out benchmarks/springback_fine_ml/fe_da_fixed --work benchmarks/springback_fine_ml/work \
    --workers 3 --strategies backing_plate --tool-reach --fe-runs 2
```

The driver is `python/scripts/springback_fine_ml.py`; its log is
`work/ml.log`. Every run goes through the content-addressed run cache in
`work/runs` (git-ignored), the data set skips the samples it has, and model
bundles and compensated shapes are reused, so an interrupted run repeats
only the runs in flight. Stage n18's data ran from 09:37 to 21:08 UTC on
2 October (with restarts); its ML part from 21:08 to 23:56 UTC.

## Limitations

* **Simulation only, nominal material.** The surrogate learns SparLab, and
  SparLab with handbook AA5754-O data. Whether either matches a real plate
  forming cell is not tested here.
* **8 test parts.** Every comparison is over 8 parts from four families.
  A 0.01 mm mean difference with a 0.017 mm standard deviation between
  parts is not resolved; neither is the ranking of GBM and MLP.
* **Penalty mismatch.** The models are trained at penalty 3 and do not see
  the penalty; the verifying runs are at penalty 10. The fine benchmark
  measured the difference on the plate at 0.000-0.009 mm RMS, small against
  the 0.05 mm dz error, but it is a bias the surrogate cannot correct.
* **One step of DA in the data.** Each training part has the uncompensated
  run and one DA step from it. The surrogate compensates by iterating on
  its own prediction, so its later iterates are judged on commands further
  from the data than one DA step. The envelope check catches some of this
  (two shots per model stopped on leaving it).
* **The design's range.** The test parts lie in the training design's range
  by construction, except elliptic_cone-s2026-0000 (outside the envelope,
  overridden). New families, materials, thicknesses or supports are not
  covered.
