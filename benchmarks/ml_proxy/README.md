# precomp.ml on the ProxySimulator - proxy numbers, NOT physics

**Every number in this directory is measured on data from the
`ProxySimulator`, an analytic stand-in for springback (rim under-forming, a
pillow on flat bases, global bending after a 3-2-1 release; formulas in
[docs/precomp_ml.md](../../docs/precomp_ml.md#data)). It is not physics and
has not been validated against anything.** These tables show that the
pipeline runs end to end and what its outputs look like. They say nothing
about how well the models predict SparLab simulations or real parts; that
benchmark follows when `sparlab_form` exists (`precomp dataset generate
--simulator sparlab`).

Reproduce (about 10 + 6 minutes on two threads of a shared 4-core machine;
the runs are deterministic - a second run reproduced every table):

```bash
python3 python/scripts/ml_proxy_demo.py --out benchmarks/ml_proxy --work <scratch> --n-per-family 8
python3 python/scripts/ml_proxy_demo.py --out benchmarks/ml_proxy/holdout_freeform --work <scratch> \
    --n-per-family 8 --holdout-family freeform --skip mlp unet
```

## Protocol (fixed before the first run)

* **Design**: all seven part families x 8 Sobol design points (seed 2026) =
  56 parts; one material (AA5754-O - one material per model); tool radius
  4-8 mm, step-down 0.2-1 mm, thickness 0.6-1.5 mm, friction 0.05-0.2;
  z-level contour tool paths; 2 mm grid on the 0.2 m blank. Three commanded
  variants per part (uncompensated, perturbed, one DA step): 168 samples,
  no generation failure. dz over the parts: 0.17-2.19 mm, 0.65-0.95 mm RMS
  on average per split.
* **Splits** (grouped by part - a part's three variants stay together): the
  family `saddle` held out of training and calibration entirely; of the other
  48 parts, 10 test, 8 calibration, 30 training ([design.csv](design.csv)).
* **Models**: GBMEnsemble (8 members, 300 iterations), MLPEnsemble (5 members,
  128-128), FieldUNet (40 x 40 grid, 12-24-48 channels, 100 epochs); 2000
  stratified points per training sample; the `time_frac` feature from the
  tool path; target divided by the elastic unloading curvature. Each
  calibrated at 90 % on the same 8 calibration parts.
* **Metrics**: dz error on the commanded part's nodes, per part, then
  averaged over parts; relative error = RMS error / RMS dz of the part.
  Coverage: pooled over the test nodes for the single calibration split, and
  the expected coverage over 200 random calibration/test partitions of the
  18 calibration + test parts (`partition_coverage`).
* **Compensation**: GBM surrogate displacement adjustment (at most 8
  iterations, stop at < 2 % improvement) of the 10 uncompensated test targets;
  the compensated shape and the target each formed by the proxy; improvement
  factor = RMS vertical deviation over the part, target formed as is /
  compensated shape formed.

The `holdout_freeform` run repeats the GBM part with `freeform` held out
instead. It was added after the first run showed that `saddle` is not
out-of-distribution for these models (below): freeform, whose dents differ
most from the other families (and which the unit tests hold out), shows the
envelope on a family that is.

## Results (proxy - not physics)

Held-out prediction error of dz and interval coverage ([summary.csv](summary.csv)):

| Model | Split | Parts | RMS dz [mm] | RMS error [mm] | max abs error [mm] | rel. error | coverage 90 % (one split) | expected coverage 90 % | worst part's coverage |
|-------|-------|------:|------:|------:|------:|------:|------:|------:|------:|
| GBMEnsemble | test | 10 | 0.790 | 0.071 | 0.215 | 8.7 % | 93.2 % | 89.5 +- 6.1 % | 3 % |
| GBMEnsemble | held-out saddle | 8 | 0.648 | 0.048 | 0.143 | 7.4 % | 96.5 % | - | 83 % |
| MLPEnsemble | test | 10 | 0.790 | 0.026 | 0.102 | 3.4 % | 97.6 % | 89.0 +- 5.6 % | 92 % |
| MLPEnsemble | held-out saddle | 8 | 0.648 | 0.037 | 0.120 | 6.0 % | 97.4 % | - | 86 % |
| FieldUNet | test | 10 | 0.790 | 0.086 | 0.273 | 11.3 % | 97.5 % | 88.4 +- 7.0 % | 90 % |
| FieldUNet | held-out saddle | 8 | 0.648 | 0.081 | 0.227 | 12.7 % | 97.6 % | - | 84 % |
| GBMEnsemble (freeform run) | test | 10 | 0.735 | 0.058 | 0.180 | 7.7 % | 93.5 % | 87.4 +- 10.2 % | 45 % |
| GBMEnsemble (freeform run) | held-out freeform | 8 | 0.945 | 0.120 | 0.288 | 13.3 % | 99.3 % | - | 95 % |

The expected coverage (mean over calibration draws) is within 1.7 % of the
nominal 90 % for all three models in the main run (2.6 % in the freeform
run); the single calibration split of 8 parts lands 4-9 % above it, and the
per-part coverage ranges down to 3 % for the worst GBM test part (a
compensated freeform part whose error, 0.19 mm RMS, is almost all a
part-wide offset of +0.19 mm). Per region
([per_region.csv](per_region.csv)) the base carries the largest error for
every model (GBM 0.11 mm RMS against 0.07-0.08 mm on rim and wall).

Envelope - share of samples flagged out of distribution by the GBM surrogate
([ood_probes.csv](ood_probes.csv)):

| Case | main run (saddle held out) | freeform run | reasons given |
|------|------:|------:|------|
| test parts (in distribution) | 0 / 30 | 0 / 30 | - |
| held-out family | saddle 0 / 24 | freeform 20 / 24 (83 %) | mean wall angle, max depth, radius of gyration |
| test targets with material DC04 (trained on AA5754-O only) | 10 / 10 | 10 / 10 | Young's modulus, r-values |
| test targets with thickness 2.5 mm (trained 0.6-1.5 mm) | 10 / 10 | 10 / 10 | thickness, elastic curvature |

The saddle - a rounded-rectangle tray like the pyramid family with a
mildly curved floor - lies inside the training envelope and is predicted as
well as the in-distribution parts (7.4 % relative error for the GBM), so not
flagging it is the right answer; freeform parts are predicted worse (13.3 %)
and flagged. No in-distribution test part was flagged in either run.

Surrogate compensation of the 10 uncompensated test targets, verified by the
proxy ([compensation.csv](compensation.csv)):

| | main run | freeform run |
|---|---:|---:|
| RMS deviation over the part, target formed as is [mm] | 0.78 (mean) | 0.72 |
| RMS deviation, compensated shape formed [mm] | 0.073 (mean), 0.067 (median) | 0.061, 0.057 |
| improvement factor, median (min - max) | 10.8 (5.3 - 41.1) | 14.0 (5.0 - 35.5) |
| largest deviation of a compensated part [mm] | 0.54 | 0.47 |
| predicted residual RMS [mm] | < 0.0002 | < 0.0002 |
| stopped by | stagnation 9, iteration limit 1 | stagnation 10 |

The predicted residual after compensation is essentially zero - the
surrogate is inverted exactly - while the achieved deviation is the model's
error at the compensated shape: the prediction is not the result, which is
why `verify_with_fea` exists. Every compensated shape was inside the
envelope.

## Files

| File | Content |
|------|---------|
| `summary.csv` | per model and split: parts, mean RMS dz and error [mm], relative error, coverage (single split, expected, worst part), share flagged |
| `per_part.csv` | per model and sample: error metrics [m], dz RMS, coverage and width at 90 %, mean std, envelope verdict, score, reasons |
| `per_region.csv` | per model, split and region (rim, wall, base, all): pooled error metrics [m], share within 0.2 mm, coverage |
| `calibration.csv` | per model and split: coverage at 50-99 %, per-part minimum and median, mean width [m] |
| `ood_probes.csv` | envelope verdicts of the probe cases |
| `compensation.csv` | per compensated test target: iterations, stop reason, predicted / achieved / uncompensated RMS and max deviation [m], improvement factor |
| `design.csv` | every sample: split, process parameters, depth, dz RMS [m] |
| `run.json` | seed, sizes, timings, command |
| `holdout_freeform/` | the same tables for the freeform hold-out (GBM only) |

Every CSV has a `data_source` column reading `proxy - not physics`.
