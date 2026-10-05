# Combined compensation v2 on the 8 test parts: retrained transfer correction + stronger uncertainty penalty

**SparLab simulation** (sparlab_form), penalty 10, DSIF + rim pass (`smoke_optimizer_edge.edge_setup(True)`). Predictions are the surrogate's, not results. 8 test parts, one verifying run each: a small sample.

**Fix 1 and fix 2 are applied together (no ablation)**: the difference v1 -> v2 cannot be split between them.

* fix 2: optimiser std weight 1.0 (v1: 0.25) on the mean squared model std, chosen a priori (the std then weighs as much as the residual); not tuned on the test parts. The transfer model has no conformal calibration (every rim-pass part is a TRAIN part of the base), so the raw ensemble std is used;
* fix 1: the 8 role=fit TRAIN parts were compensated with the v1 model + the new penalty, simulated, and added to work/data_rim as `compensated` samples; the correction was refitted on 8 uncompensated + 8 compensated fit samples (the 4 holdout parts and all test parts kept out); its form picked by leave-one-part-out on TRAIN data only.

## Test parts: vertical RMS of formed - target [mm]

plate columns: backing plate (stage_n26 headline); v1 = combined_rim (transfer v1, std weight 0.25); v2 = this run.

| part | plate uncompensated | ML-MLP plate | combined v1 | combined v2 | v1 predicted | v2 predicted |
|---|--:|--:|--:|--:|--:|--:|
| dome-s2026-0000 | 0.307 | 0.211 | 0.071 | 0.073 | 0.073 | 0.089 |
| dome-s2026-0001 | 0.290 | 0.217 | 0.145 | 0.138 | 0.095 | 0.193 |
| elliptic_cone-s2026-0000 | 0.278 | 0.258 | 0.202 | 0.205 | 0.186 | 0.200 |
| elliptic_cone-s2026-0001 | 0.305 | 0.258 | 0.205 | 0.224 | 0.137 | 0.156 |
| pyramid-s2026-0000 | 0.290 | 0.196 | 0.103 | 0.138 | 0.067 | 0.096 |
| pyramid-s2026-0001 | 0.308 | 0.268 | 0.216 | 0.255 | 0.185 | 0.216 |
| truncated_cone-s2026-0000 | 0.280 | 0.234 | 0.205 | 0.214 | 0.100 | 0.136 |
| truncated_cone-s2026-0001 | 0.344 | 0.344 | 0.365 | 0.374 | 0.269 | 0.356 |
| **mean** | **0.300** | **0.248** | **0.189** | **0.203** | **0.139** | **0.180** |

v2 simulated on 8 of 8 parts (means over the available ones).

| mean over parts [mm] | plate uncompensated | ML-MLP plate | combined v1 | combined v2 |
|---|--:|--:|--:|--:|
| upper-band RMS | 0.445 | 0.367 | 0.129 | 0.130 |
| deep RMS | 0.149 | 0.115 | 0.213 | 0.231 |

Predicted vs simulated (mean vertical RMS): v1 0.139 predicted / 0.189 simulated; v2 0.180 predicted / 0.203 simulated. Predicted upper-band / deep RMS: v1 0.085 / 0.155, v2 0.101 / 0.205.

## Transfer correction: old vs new (TRAIN rim-pass data only)

Leave-one-part-out over the 8 fit parts: each part's samples (uncompensated and compensated) left out together, the correction refitted on the other 7 parts, dz prediction error RMS [mm] on the left-out samples (mean). 'old' = v1's protocol (linear, uncompensated samples only).

| correction | trained on | all 16 | uncompensated 8 | compensated 8 |
|---|---|--:|--:|--:|
| base MLP (no correction) | - | 0.501 | 0.486 | 0.515 |
| old: linear | uncompensated | 0.108 | 0.100 | 0.116 |
| new: gbm (chosen) | uncompensated + compensated | 0.100 | 0.095 | 0.106 |
| new: linear | uncompensated + compensated | 0.106 | 0.097 | 0.115 |

Chosen: **gbm** (mean LOPO dz RMS over the left-out fit samples (both samples of a part leave together; TRAIN data only); 'gbm' only if more than 5 % below 'linear'); v2 accepted by its own validation: True. v1 (as used) on the 8 compensated fit runs (not in its fit): dz error 0.115 mm.

4 holdout TRAIN parts (uncompensated rim-pass runs, in neither fit), dz error [mm]: base 0.526, v1 0.124, v2 0.118.

## The 8 compensated fit runs (TRAIN, data for fix 1)

v1 model + std weight 1.0, DSIF + rim pass, penalty 10; vertical RMS [mm].

| part | rim pass uncompensated | predicted | simulated | upper band | deep | status |
|---|--:|--:|--:|--:|--:|---|
| dome-s2026-0010 | 0.294 | 0.051 | 0.064 | 0.070 | 0.056 | ok |
| dome-s2026-0011 | 0.389 | 0.128 | 0.173 | 0.121 | 0.194 | ok |
| elliptic_cone-s2026-0008 | 0.237 | 0.100 | 0.170 | 0.113 | 0.203 | ok |
| elliptic_cone-s2026-0015 | 0.324 | 0.119 | 0.255 | 0.094 | 0.312 | ok |
| pyramid-s2026-0006 | 0.184 | 0.075 | 0.139 | 0.062 | 0.179 | ok |
| pyramid-s2026-0023 | 0.336 | 0.118 | 0.141 | 0.114 | 0.155 | ok |
| truncated_cone-s2026-0008 | 0.531 | 0.207 | 0.212 | 0.149 | 0.244 | ok |
| truncated_cone-s2026-0015 | 0.392 | 0.231 | 0.287 | 0.091 | 0.367 | ok |

Made by `python/scripts/combined_rim2.py` (run: `benchmarks/springback_fine_ml/combined_rim2/run.sh`). v1 results: `benchmarks/springback_fine_ml/combined_rim/`. Per-part rows: `results.csv`, fit runs: `fit_runs.csv`, correction fit: `refit.json`.
