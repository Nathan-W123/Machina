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

Reproduce (about 25 + 6 minutes on two threads of a shared 4-core machine;
the runs are deterministic - a second run reproduced every table):

```bash
python3 python/scripts/ml_proxy_demo.py --out benchmarks/ml_proxy --work <scratch> --n-per-family 10
python3 python/scripts/ml_proxy_demo.py --out benchmarks/ml_proxy/holdout_freeform --work <scratch> \
    --n-per-family 10 --holdout-family freeform --skip mlp unet lofo
```

## Protocol (revision 2, fixed before its first run)

* **Design**: all seven part families x 10 Sobol design points (seed 2026) =
  70 parts; one material (AA5754-O - one material per model); tool radius
  4-8 mm, step-down 0.2-1 mm, thickness 0.6-1.5 mm, friction 0.05-0.2;
  z-level contour tool paths; 2 mm grid on the 0.2 m blank. Three commanded
  variants per part (uncompensated, perturbed, one DA step): 210 samples,
  no generation failure. dz over the parts: 0.36-1.34 mm RMS per sample,
  0.64-0.75 mm on average per split.
* **Splits** (grouped by part - a part's three variants stay together): the
  family `saddle` held out of training and calibration entirely; of the other
  60 parts, 12 test, 12 calibration, 36 training ([design.csv](design.csv)).
* **Models**: GBMEnsemble (8 members, 300 iterations), MLPEnsemble (5 members,
  128-128), FieldUNet (40 x 40 grid, 12-24-48 channels, 100 epochs); 2000
  stratified points per training sample; the `time_frac` feature from the
  tool path; target divided by the elastic unloading curvature. Each
  calibrated at 90 % on the same 12 calibration parts - over whole parts,
  with the (k + 1) correction ([docs](../../docs/precomp_ml.md#uncertainty-conformal-intervals)).
* **Metrics**: dz error on the commanded part's nodes, per sample, then
  averaged; relative error = RMS error / RMS dz of the sample. Coverage per
  part (its three variants pooled): the mean over parts (the quantity the
  conformal guarantee is about), the share of parts covered at 90 % or
  better and the worst part, for the single calibration split and as
  expected over 200 random calibration/test partitions (12 + 12) of the
  24 calibration + test parts (`partition_coverage`); the interval width
  against a constant-width interval calibrated on the same parts.
* **Compensation**: GBM surrogate displacement adjustment of the 12
  uncompensated test targets - refused for a target outside the envelope -
  at most 8 predictions, stopping at < 2 % improvement; the compensated shape
  and the target each formed by the proxy; improvement factor = RMS vertical
  deviation over the part, target formed as is / compensated shape formed.
  The `stop_*` columns repeat it with `interval_stop` (stop once the
  predicted error is within the interval half-width).
* **Envelope probes** (GBM): test parts, the held-out family, and the test
  targets with a material (DC04), a thickness (2.5 mm), a release
  ("clamped_only") and a tool path (spiral) the model was not trained on.
  The spiral probe takes only the targets a spiral can form: the first
  freeform run stopped on one it cannot (assess needs the tool path for the
  features) - the one change after a first run, and it does not touch the
  main run's targets, which a spiral can all form.
* **Leave each family out**: for each of the seven families, a
  GBMEnsemble (4 members, 150 iterations, 1000 points per sample) and its
  envelope trained on the other families without a fixed grouped 20 % of all
  parts; flagged share and relative error on the left-out family and on
  those in-distribution parts (errors on at most 1500 evenly spread part
  nodes per sample, the envelope's node set).

**Why revision 2.** A review of the first version found that the conformal
calibration counted nodes, not parts, in its finite-sample correction and
grouped calibration by sample. Calibrated over whole parts, 90 % needs at
least 9 calibration parts, and revision 1 had 8; this revision draws 10
parts per family and calibrates on 12. The design seeds changed too
(families are now seeded by name), so the parts differ from revision 1,
whose tables are in the git history (commit fcfe8b0); the numbers are not
comparable one to one.

The `holdout_freeform` run repeats the GBM part with `freeform` held out
instead - the one family the envelope flags when it is left out (below).

## Results (proxy - not physics)

Held-out prediction error of dz and interval coverage at 90 %
([summary.csv](summary.csv), [calibration.csv](calibration.csv)):

| Model | Split | Parts | RMS dz [mm] | RMS error [mm] | rel. error | coverage per part, one split | parts covered >= 90 % | expected coverage | width / constant width (one split; expected) |
|-------|-------|------:|------:|------:|------:|------:|------:|------:|------:|
| GBMEnsemble | test | 12 | 0.734 | 0.049 | 6.7 % | 99.3 % | 12 / 12 | 97.0 +- 3.1 % | 0.72; 0.80 |
| GBMEnsemble | held-out saddle | 10 | 0.641 | 0.059 | 9.3 % | 97.3 % | 9 / 10 | - | 0.61 |
| MLPEnsemble | test | 12 | 0.734 | 0.023 | 3.1 % | 99.8 % | 12 / 12 | 97.1 +- 2.5 % | 0.67; 0.74 |
| MLPEnsemble | held-out saddle | 10 | 0.641 | 0.024 | 3.9 % | 98.9 % | 10 / 10 | - | 0.62 |
| FieldUNet | test | 12 | 0.734 | 0.071 | 9.5 % | 99.8 % | 12 / 12 | 96.7 +- 3.1 % | 0.81; 0.87 |
| FieldUNet | held-out saddle | 10 | 0.641 | 0.064 | 10.2 % | 99.7 % | 10 / 10 | - | 0.81 |
| GBMEnsemble (freeform run) | test | 12 | 0.729 | 0.047 | 6.8 % | 99.8 % | 12 / 12 | 96.6 +- 3.7 % | 1.20; 1.06 |
| GBMEnsemble (freeform run) | held-out freeform | 10 | 0.755 | 0.128 | 17.6 % | 97.0 % | 9 / 10 | - | 2.05 |

With 12 calibration parts the guarantee - an expected coverage of a new part
of at least 90 % - is conservative by construction: up to 90 % x 13 / 12 =
97.5 % when parts behave alike, and the expected coverages are 96.7-97.1 %.
The worst test part is covered at 96 % (GBM) or better; one calibration split
lands 2-3 % above the expectation. 95 % and 99 % are not supported by 12
parts (19 and 99 needed) and are reported as such. In the main run the
model's std makes the 90 % intervals 13-26 % narrower than a constant width
calibrated on the same parts (expected width ratio 0.74-0.87); in the
freeform run it does not (1.06), nor on the smaller data of the unit tests
(1.2-1.4) - there it only moves width to the parts with the largest errors.
On the held-out freeform family the GBM's intervals are twice as wide as the
constant width (ratio 2.05): its std grows where it extrapolates, and 9 of
the 10 freeform parts are still covered at 90 %. Per region ([per_region.csv](per_region.csv)) the base carries the
largest error for the GBM and MLP (GBM 0.070 mm RMS against 0.035-0.056 mm
on wall and rim).

Envelope - share of samples flagged out of distribution by the GBM surrogate
([ood_probes.csv](ood_probes.csv)):

| Case | main run (saddle held out) | freeform run | reasons given |
|------|------:|------:|------|
| test parts (in distribution) | 3 / 36 (one part) | 4 / 36 | area, radius of gyration, tool radius |
| held-out family | saddle 0 / 30 | freeform 27 / 30 (90 %) | freeform: mean wall angle, volume, radius of gyration |
| test targets with material DC04 (trained on AA5754-O only) | 12 / 12 | 12 / 12 | Young's modulus, r-values |
| test targets with thickness 2.5 mm (trained 0.6-1.5 mm) | 12 / 12 | 12 / 12 | thickness, elastic curvature |
| test targets released "clamped_only" (trained 3-2-1) | 12 / 12 | 12 / 12 | setup.release |
| test targets with a spiral tool path (trained contour; the targets a spiral can form) | 12 / 12 | 10 / 10 | setup.toolpath_style |

Leave each family out - GBMEnsemble (4 x 150) and envelope per left-out
family ([leave_family_out.csv](leave_family_out.csv)):

| Left-out family | flagged: family | flagged: in-distribution test | rel. error: family | rel. error: in-distribution test |
|-----------------|------:|------:|------:|------:|
| dome | 0.00 | 0.03 | 4.7 % | 7.7 % |
| elliptic_cone | 0.00 | 0.09 | 12.4 % | 8.7 % |
| freeform | 0.87 | 0.08 | 17.6 % | 7.2 % |
| pyramid | 0.00 | 0.00 | 10.4 % | 8.9 % |
| saddle | 0.00 | 0.03 | 9.4 % | 8.4 % |
| truncated_cone | 0.20 | 0.00 | 4.4 % | 9.0 % |
| two_level | 0.00 | 0.08 | 6.5 % | 7.7 % |

The envelope detects new **descriptors**, not new family names: of the seven
families left out in turn only freeform - whose dents no other family has -
is flagged (87 %), and it is also the one predicted worst (17.6 % against
7.2 %). Dome, saddle and two-level parts lie inside the envelope of the other
families and are predicted about as well as they are; the elliptic cone is
predicted worse (12.4 % against 8.7 %) without being flagged - a failure the
envelope cannot see. 0-9 % of in-distribution test samples are flagged.

Surrogate compensation of the 12 uncompensated test targets, verified by the
proxy ([compensation.csv](compensation.csv)):

| | main run | freeform run |
|---|---:|---:|
| targets compensated (refused: outside the envelope) | 11 (1) | 10 (2) |
| RMS deviation over the part, target formed as is [mm] | 0.69 (mean) | 0.71 |
| RMS deviation, compensated shape formed [mm] | 0.045 (mean), 0.033 (median) | 0.040, 0.036 |
| improvement factor, median (min - max) | 18.2 (7.1 - 59.2) | 17.7 (11.6 - 45.3) |
| largest deviation of a compensated part [mm] | 0.30 | 0.26 |
| predicted residual RMS [mm] / mean interval half-width [mm] | < 0.0013 / 0.15 | < 0.0001 / 0.19 |
| stopped by | stagnation 10, iteration limit 1 | stagnation 10 |
| with `interval_stop`: predictions per target; RMS deviation [mm]; factor median | 2; 0.049 (mean); 17.9 | 2; 0.040 (mean); 21.2 |

The refused targets are in-distribution parts the envelope flags (a
truncated cone: area, radius of gyration, tool radius; in the freeform run
also a dome: wall angles, volume, step-down): with 99 % thresholds from 36
parts, one or two in twelve in-distribution parts are refused until
overridden.
The predicted residual after compensation is essentially zero - the
surrogate's map is inverted - while the achieved deviation is the model's
error at the compensated shape: the prediction is not the result, which is
why `verify_with_fea` exists, and why the predicted residual is reported
with its interval half-width (0.15-0.19 mm here). Stopping within that
half-width saved two thirds of the predictions; the verified deviation was
9 % larger in the main run, about the same in the freeform run (0.040 mm
either way; median factor 21.2 against 17.7) and 4 % smaller on the
unit-test data, so it is an option, not the default.

## Files

| File | Content |
|------|---------|
| `summary.csv` | per model and split: parts, mean RMS dz and error [mm], relative error, coverage per part at 90 % (single split and expected), share of parts covered, worst part, widths and width ratio, share flagged |
| `per_part.csv` | per model and sample: error metrics [m], dz RMS, coverage and width at 90 %, mean std, envelope verdict, score, reasons |
| `per_region.csv` | per model, split and region (rim, wall, base, all): pooled error metrics [m], share within 0.2 mm, coverage |
| `calibration.csv` | per model, split and level (50-99 %): supported or not, coverage per part and pooled, worst / 10th-percentile / median part, share of parts at the level, width, constant-width baseline and ratio |
| `ood_probes.csv` | envelope verdicts of the probe cases |
| `leave_family_out.csv` | per left-out family: flagged share and relative error, left-out family and in-distribution test |
| `compensation.csv` | per test target: iterations, stop reason, predicted / achieved / uncompensated RMS and max deviation [m], improvement factor, interval half-width; `stop_*` with `interval_stop` |
| `design.csv` | every sample: split, process parameters, depth, dz RMS [m] |
| `run.json` | seed, protocol revision, sizes, timings, command |
| `holdout_freeform/` | the same tables for the freeform hold-out (GBM only, no leave-family-out) |

Every CSV has a `data_source` column reading `proxy - not physics`.
