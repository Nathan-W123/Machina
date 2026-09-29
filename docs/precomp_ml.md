# precomp.ml: learned springback prediction and surrogate pre-compensation

`precomp.ml` predicts, before the first forming attempt, the surface a
commanded tool path will form, and inverts that prediction by displacement
adjustment so that the first part lands within tolerance - and says when the
prediction should not be trusted. It is built on the `precomp` package
([docs/precomp.md](precomp.md)): its height maps, part families, tool paths,
the SparLab bridge, metrology and compensation are used as they are.

**Where the numbers come from.** Every metric the package reports carries its
data source: `SparLab simulation`, `scan`, or `proxy - not physics`. The
C++ forming solver `sparlab_form` does not exist on this branch yet, so
everything measured so far - the tests and the demonstration in
[Results on proxy data](#results-on-proxy-data) - comes from the
`ProxySimulator`, an analytic stand-in that is **not physics**. Those numbers
show that the pipeline works and what its outputs look like; they say nothing
about how well a model predicts SparLab or a real part. The SparLab path is
exercised end to end with the test-double solver (formed = 0.9 x commanded)
and will be benchmarked when the solver lands.

## Contents

* [The learning problem](#the-learning-problem)
* [Assumptions](#assumptions)
* [Modules](#modules)
* [Features](#features)
* [Data](#data)
* [Models](#models)
* [Uncertainty: conformal intervals](#uncertainty-conformal-intervals)
* [The training envelope](#the-training-envelope)
* [Model bundles](#model-bundles)
* [Evaluation protocol](#evaluation-protocol)
* [Compensation on the surrogate](#compensation-on-the-surrogate)
* [Choosing the next runs](#choosing-the-next-runs)
* [Command line and API](#command-line-and-api)
* [Results on proxy data](#results-on-proxy-data)
* [Limitations](#limitations)

## The learning problem

The target is the vertical deviation of the tool-side surface after
springback (unclamped),

    dz(x, y) = z_formed(x, y) - z_commanded(x, y)    [m],

on the grid of the commanded surface, as a function of the **commanded**
surface, the process and the material. It is not a function of the target:
during compensation the model is queried at compensated shapes, which are
deeper and steeper than the target. The data sets therefore contain three
commanded variants per design point - the target itself, a randomly
perturbed target and a once-compensated target - so the model has seen the
neighbourhood it will be asked about. dz is what displacement adjustment
inverts: formed(c) = c + dz_hat(c). The signed normal deviation is derived
for reporting (`precomp.metrology.signed_deviation`).

## Assumptions

* **Height fields.** Parts, commanded surfaces and formed surfaces are height
  fields z(x, y) <= 0 of the tool side (no overhangs), as everywhere in
  `precomp`.
* **One material per model.** Material descriptors are features, but a
  design over a handful of library alloys gives the model a handful of points
  in material space; a new material is outside the envelope and is flagged.
  Train one model per material (and per machine setup) - the demonstration
  does.
* **One fixture.** The blank size, clamp frame and 3-2-1 release enter
  through placement features and the global descriptors; a model trained on
  one fixture has not seen another.
* **Grid.** The grid spacing must not exceed the ring width
  (4 tool radii / 4 rings = one tool radius by default); refused otherwise.
* **Parts centred on the blank**, as `precomp` builds them.
* **Tool path.** The `time_frac` feature uses the setup's own tool path
  (`precomp.fea.make_toolpath`); a spiral path needs one loop per level, so
  multi-pocket parts need `toolpath_style="contour"` (the deck builder would
  refuse them too).

## Modules

| Module | Public API |
|--------|------------|
| `features` | `FEATURE_SCHEMA_VERSION`, `FeatureSpec(name, unit, description, kind)`, `FeatureConfig(scales, rings, ring_max, time_source, part_eps, rim_width)` with `.specs()`, `.names`, `.schema_hash()`; `point_features(commanded, setup, toolpath=None, points=None, *, config) -> (X, names)`; `feature_maps(...) -> FeatureMaps`; `global_features(commanded, setup, config) -> (d, names)`; `region_labels(...)` |
| `dataset` | `Sample(sample_id, commanded, formed, setup, source, fidelity, kind, part, part_id, target, toolpath, provenance)`; `Dataset.create(root, created_at=...)` / `Dataset(root)`: `append`, `load`, `ids`, `index`, `filter`, `record_failure`, `failures`, `to_table(features, target="dz", points_per_sample, mask, rng, ids, stratify, prior, n_jobs) -> Table` (unpacks as `X, y, groups`); `build_table(samples, ...)`; `grouped_split`, `grouped_kfold`, `family_split`, `check_disjoint`; `SOURCE_LABELS`, `source_label` |
| `generate` | `DesignSpace`, `design_points(space, n_per_family, seed) -> [DesignPoint]`, `PerturbationSpec`, `perturbed_commanded`, `compensated_commanded`, `generate(dataset, points, simulator, *, created_at, kinds, seed, compensator) -> GenerationReport`, `simulate_samples`; `SparlabSimulator(work_dir, max_workers)`, `ProxySimulator(params)`, `ProxyParams`, the `Simulator` protocol |
| `models` | `GBMEnsemble`, `MLPEnsemble`, `FieldUNet`, `ResidualModel(prior, learner)`, `TransferModel(base, correction)`, `FEAPrior(work_dir, overrides)`, `model_from_state`, `thread_limit`; protocols `DeviationModel` (point: `fit(X, y, groups)`, `predict(X) -> (mean, std)`, `feature_names`, `to_state`/`from_state`) and `FieldDeviationModel` (`fit_samples`, `predict_field(commanded, setup)`) |
| `uncertainty` | `ConformalCalibrator(eps, eps_fraction, mondrian)`: `fit(mu, std, y, groups, region)`, `quantile(level)`, `intervals(mu, std, level)`; `intervals(...)`, `coverage_report(...)`, `partition_coverage(...)` |
| `ood` | `OODEnvelope(quantile, k, ...)`: `fit(point_X, point_groups, part_D, part_groups, ...)`, `assess_features(point_X, part_d) -> {in_envelope, part_score, fraction_points_out, reasons, ...}` |
| `surrogate` | `DeviationSurrogate(model, features, data_source, training, calibrator, ood, domain)`: `predict_deviation(commanded, setup) -> (mean, std)`, `predict_interval(commanded, setup, level)`, `assess(commanded, setup)`; `train_surrogate(model, train, calibration=..., features, points_per_sample, mask, seed)`, `transfer_surrogate(base, real, calibration)`, `calibrate`, `fit_envelope` |
| `registry` | `save_model(model, directory, manifest)`, `load_model(directory)`, `read_manifest`, `check_schema` |
| `evaluate` | `evaluate_surrogate(surrogate, samples, tolerance, level) -> EvaluationReport` (`per_part`, `per_region`, `calibration`, `ood` frames, `summary()`, `to_csv()`); `cross_validate(samples, make_model, n_splits)`; `family_holdout(samples, family, make_model)` |
| `compensate` | `SurrogatePredictor(model, setup)`, `surrogate_compensate(target, setup, surrogate, iterations, tolerance, stagnation, ...) -> SurrogateCompensation`, `verify_with_fea(result, setup, work_dir)`, `verify_with_simulator(result, setup, simulator)` |
| `active` | `Candidate(commanded, setup, label)`, `rank_candidates(surrogate, candidates, n_select, weight_std, weight_ood, diversity) -> DataFrame` |
| `cli` | `main(argv)`: `precomp dataset generate|info`, `precomp train`, `precomp evaluate`, `precomp active` |

`DeviationSurrogate` implements the `FieldModel` protocol of
`precomp.compensation` (plus `predict_interval` and `assess`), so
`precomp.api.predict` / `compensate`, `precomp compensate --method surrogate`
and `CompositePredictor` use it without knowing about machine learning; the
API and the command line also accept a bundle directory.

## Features

Per grid node, `FEATURE_SCHEMA_VERSION = "1"`, 57 features in five kinds
(`FeatureConfig().specs()` gives each with its unit and description):

| Kind | Features | Units |
|------|----------|-------|
| local (invariant when the part moves in x-y) | `depth`, `depth_frac`; `wall_angle_s1..3` (Gaussian-smoothed at 1, 2.5, 6 tool radii); `mean_curv_s1..3` = H sigma and `gauss_curv_s1..3` = K sigma^2 (the package's moment-normalised derivatives, dimensionless); `rim_dist` (signed distance to the rim, + inside); `radial_norm` (distance to the part centroid / equivalent radius sqrt(A/pi)); `polar_sin`, `polar_cos` (regularised by the tool radius R: dy / sqrt(r^2 + R^2) - the angle is undefined at the centroid); `time_frac` (pseudo-time of the tool path where the tool passes nearest: the node offset by R along its normal, nearest in-contact path point; or the depth fraction with `time_source="depth"`); `ring{1..4}_{mean,min,max}` (height of the neighbours in four annuli out to 4 R, relative to the node) | m, -, rad |
| placement | `clamp_dist` (L-infinity distance to the clamped frame), `blank_radius` (distance to the blank centre / free half-width) | m, - |
| global (per part) | `max_depth`, `area`, `volume`, `mean_wall_angle`, `max_wall_angle`, `perimeter_over_area` (rim contour length), `aspect` (principal axes of the plan area), `depth_gyration_major/minor` (radii of gyration of the depth field) | m, m^2, m^3, rad, 1/m, - |
| process | `tool_radius`, `step_down`, `thickness`, `friction`, `elastic_curvature` = 3 sigma_0.2 / (E t), the unloading curvature of a fully plastic bend | m, -, 1/m |
| material | `youngs_modulus`, `poisson_ratio`, `yield_stress`, `flow_stress_20`, `flow_stress_50` (isotropic flow stress at plastic strain 0.2, 0.5), `hardening_modulus`, `voce_Q`, `voce_rate`, `r0`, `r45`, `r90` (1 when isotropic), `backstress_sat` (monotonic back stress at strain 0.5: Prager + Armstrong-Frederick, tending to sum C / gamma), `yield_over_E` (sigma_0.2 / E) | Pa, - |

`FeatureConfig.schema_hash()` fingerprints the version, names, units and
configuration; a bundle stores it and refuses to load under another. Every
feature is checked finite. Measured: a 200 x 200 grid featurises in 0.6 s
(0.9 s including the tool path) - the specification asked for 2 s; moving a
pyramid by (7, -4) grid spacings leaves every non-placement feature equal to
1e-9 of its range at the corresponding nodes (`test_ml_features`).

Region codes (`region_labels`): 0 flange, 1 rim (part nodes within
`rim_width` = 2 tool radii of the boundary), 2 wall (wall angle > 5 deg,
`metrology.wall_mask`), 3 base.

## Data

**Samples and data sets.** A `Sample` holds the commanded and formed
surfaces on one grid, the target it was derived from, the setup and part as
JSON, the kind of commanded variant, a `part_id` shared by the variants of one
design point, the source and fidelity, and provenance: `created_at` (a string
supplied by the caller - nothing in the library reads the clock), and for
simulations the SparLab version and the deck hash (the content hash of the run
cache), the runtime and whether the run came from the cache. On disk:

    dataset.json            format, created_at, description
    index.jsonl             one line per sample (append-only; the last line of an id wins)
    samples/<id>.npz        commanded / formed / target heights and masks, grid, tool path
    samples/<id>.json       part, setup, provenance, data source label
    failures.jsonl          jobs that produced no sample, with their reason

**Tables.** `to_table` featurises each sample, keeps the nodes of the mask
(`"part"` by default, or `"all"`), and draws `points_per_sample` of them
shared equally among the regions present - rim, wall, base (and flange) - so a
large flat floor does not dominate; a region with fewer nodes gives all it has
and the rest goes to the others. Rows carry the sample id (`groups`), part
id, family, region and source.

**Splits** are grouped: by part (the default - the three variants of a design
point are never on both sides), by sample, or by family (leave-family-out);
`grouped_kfold` tests every group once; `check_disjoint` verifies a split.

**Design of experiments.** Per family, a scrambled Sobol sequence
(`scipy.stats.qmc`, seeded per family) over the family's parameter bounds, the
process ranges (tool radius 4-8 mm, step-down 0.2-1 mm, thickness 0.6-1.5 mm,
friction 0.05-0.2 by default) and the material index; draws the family
refuses, or that do not fit the unclamped window with the tool, are skipped
and the sequence continues, so a seed gives the same design every time.
Freeform parts are drawn by `Freeform.sample` from a seed taken from the Sobol
point.

**Variants.** `uncompensated` (commanded = target); `perturbed` (target + a
Gaussian random field of RMS 0.2-1 mm and correlation length 8-25 mm, faded
in over 2 tool radii inside the rim and out where the wall comes within 10 deg
of the 65 deg forming limit, held at zero on the flange, projected onto the
formable set); `compensated` (one displacement-adjustment step from the
simulated uncompensated part, or from any compensator, e.g. a surrogate of an
earlier round).

**Generation** (`generate`) is resumable - samples present are skipped, and a
SparLab run comes from the content-addressed cache - and records every failure
in `failures.jsonl` with its reason (a compensated variant whose uncompensated
run failed is recorded as "dependency failed"); nothing is dropped silently.
Generating the proxy data set also exposed a stall in
`precomp.compensation.limit_wall_angle`, fixed on this branch: the cone slope
is now bisected towards tan(limit) / sqrt(2) when reducing it by the measured
excess does not move the measure.

**The ProxySimulator** (`generate.ProxySimulator`) - NOT physics. With d the
commanded depth, theta the wall angle, s the signed rim distance and
sigma_0.2 / E, t, R, step-down and friction from the setup:

* rim under-forming (dz > 0, too shallow): A sin(theta at 2.5 R) exp(-s / 2.5 R)
  inside the rim and exp(-(s / R)^2) outside, A = 0.8 mm
  (sigma_0.2 / E / 3.8e-3)^0.75 (1 mm / t)^0.5 (1 + (mu - 0.1)) (1 - e^(-D / 10 mm));
* pillow on flat bases (theta < 12 deg): P (1 - exp(-(e / 12 mm)^2)), e the
  distance from the flat region's edge, P = 0.5 mm (sigma_0.2 / E / 3.8e-3)^0.5
  sqrt(R / 5 mm) sqrt(step-down / 0.5 mm) (1 mm / t)^0.25 min(D / 10 mm, 1);
* both smoothed at R / 2, plus global bending after the 3-2-1 release,
  (kx (a^2 - x^2) + ky (a^2 - y^2)) / 2, zero at the supports x = y = a,
  kx + ky = 2 kappa, kappa = 0.1 /m (sigma_0.2 / E / 3.8e-3)(D / 25 mm)(1 mm / t),
  shared between x and y by the second moments of the depth field (a quarter
  of it, inside the window, with `release="clamped_only"`).

It gives dz of 0.15-2.4 mm on the default five-material design, depends on the commanded
shape (so DA on it is not trivial), and has local and global structure. It is
labelled `proxy - not physics` in every sample, table, manifest and report.

## Models

All models are deterministic for a seed (bootstraps, initialisations,
shuffles and dropout masks derive from it) and cap their threads at
`n_threads` or `$PRECOMP_ML_THREADS` (scikit-learn/OpenMP via threadpoolctl,
and torch). scikit-learn is imported with the module, not lazily:
threadpoolctl caps only runtimes already loaded, and a lazy import inside the
limited block ran unlimited - twenty times slower on a loaded machine.

**Target scale.** The point models and the U-Net learn dz divided by the
`elastic_curvature` feature 3 sigma_0.2 / (E t) (option `target_scale`, stored
in the state; `None` learns dz as is). Springback scales with it to first
order; dividing it out spares the model building the product of material,
thickness and depth from a few dozen parts. On a proxy design with one
material and five process parameters the held-out error fell from 0.14 to
0.09 mm; across five materials from 0.17 to 0.09 mm.

* **`GBMEnsemble`** - M (8) `HistGradientBoostingRegressor` members (300
  iterations, 31 leaves, 40 samples per leaf, 80 % of the features per split),
  each on a bootstrap of whole samples; mean and standard deviation over the
  members; optional quantile members for a 5/95 % band.
* **`MLPEnsemble`** (torch) - M (5) MLPs (128, 128, SiLU) with a
  (mean, log-variance) head trained on the Gaussian NLL after a few epochs of
  squared error, Adam, early stopping on a grouped validation split; the
  mixture's mean and variance; standardisation in the state.
* **`FieldUNet`** (torch) - a three-level U-Net on a fixed training grid over
  the window every training grid covers: channels commanded height / 30 mm,
  depth fraction, wall angle at R, clamp distance, part indicator, and the
  part descriptors broadcast; output dz and its log-variance, Gaussian NLL
  weighted 1 on the part and 0.2 elsewhere; Monte-Carlo dropout (16 seeded
  passes). It sees the whole part, the path to the non-local unclamping
  springback a point model reaches only through its global descriptors.
  Predictions are resampled back to the commanded grid; nodes outside the
  window are reported by `assess` (`fraction_outside_window`).
* **`ResidualModel(prior, learner)`** - the learner fits y - prior, with the
  prior's value as an extra feature column `prior_dz`; the prior is anything
  with `prior_deviation(commanded, setup)`: a coarse SparLab run
  (`FEAPrior(work_dir, {"element_size": 5e-3, "layers": 1})`, needing the
  executable whenever evaluated) or a closed-form estimate.
* **`TransferModel(base, correction)`** - the simulation-trained base frozen,
  plus a correction fitted on a few real samples: Bayesian ridge of the
  residual on the standardised [features, base mean, base std] (columns
  constant over the real data dropped), its predictive std combined with the
  base std in quadrature; or a small, heavily regularised GBM. It reports
  `n_real`; with zero real samples it is the base, bit for bit.
  `transfer_surrogate` keeps the base's envelope and gives intervals only when
  calibrated on held-out real parts - never with a calibration from another
  data source.

## Uncertainty: conformal intervals

The ensemble spread or predicted variance is not calibrated by itself.
`ConformalCalibrator` is fitted on held-out **whole parts** - never on parts
the model was trained on (`train_surrogate` refuses an overlap):

    score r = |y - mu| / (std + eps),   eps = 0.1 median(std) of the calibration set,
    q = the ceil((n + 1) level)-th smallest of the n scores,
    interval = mu -+ q (std + eps).

A level the calibration set cannot support (ceil((n + 1) level) > n) is an
error, not an infinite interval. `mondrian=True` takes one quantile per
region (rim / wall / base), equalising coverage across regions.

**What the guarantee is, and is not.** Split conformal prediction guarantees
coverage >= level on average over calibration draws for a new point
exchangeable with the calibration points. Nodes of one part are strongly
correlated - most of a part's error is one smooth field - so the effective
sample size is the number of parts, and a single calibration set of a dozen
parts lands several per cent from the expected coverage (6-8 % standard
deviation on proxy data). `partition_coverage` measures the expected coverage
honestly: it re-splits the held-out parts many times into calibration and test
parts and reports the mean and spread. The test suite asserts the mean within
3 % of nominal at 80 % and 90 %; reports show both the single-split coverage
and, where computed, the expected one, and the per-part minimum.

## The training envelope

`OODEnvelope` standardises features and measures two distances at two levels:
Mahalanobis under a Ledoit-Wolf shrunk covariance, and the mean distance to
the 5 nearest training vectors; for every node's feature vector and for the
part's descriptor vector (global, process, material). Thresholds are the 99 %
quantiles of distances of **held-out training parts** - nodes scored fold by
fold (parts split by part id into 5 folds), parts leave-one-part-out - so an
in-distribution part scores like a new part from the training distribution. A
part is inside when its part score (the larger distance over its threshold)
is <= 1 and at most `max_fraction_out` of its nodes exceed a node threshold;
that share is itself the 99 % quantile over the held-out training parts (at
least 5 %). An 80/20 hold-out for the node thresholds, tried first, flagged a
third of in-distribution parts; the folds brought that to 7 %, with an unseen
family flagged throughout. `assess` returns `in_envelope`, `part_score`,
`fraction_points_out`, the thresholds, and `reasons`: the features outside
their training range (largest excess first), then the largest contributions
to the Mahalanobis distance.

## Model bundles

`save_model(surrogate, dir, manifest)` writes `manifest.json`,
`model.joblib` (the state; scikit-learn estimators pickled) and
`torch/<name>.pt` (state dicts, tensors only). The manifest holds the model
class, `feature_schema_version`, `feature_schema_hash`, the feature names and
their SHA-256, the feature configuration, the data source of the training
data, the training record (sample ids up to 5000 and always their SHA-256 and
count, parts, families, sources, fidelities, SparLab versions, seed, sizes),
metrics (each block must state its `data_source`), the conformal quantiles
at 50/80/90/95/99 %, the OOD thresholds, package versions, the git revision
(`-dirty` when the tree has changes), `created_at` given by the caller, notes,
and the SHA-256 of every file. `load_model` refuses another schema version,
another schema hash for the same configuration, feature names that do not
match their hash, and any file that does not match its hash - with a message
saying which and what to do. **`model.joblib` is a pickle: load bundles you
trust only.**

## Evaluation protocol

A model is scored only on whole parts it has not seen (`evaluate_surrogate`
refuses samples recorded in its training or calibration): a held-out test
split, grouped K-fold cross-validation (`cross_validate`, each fold trained
and calibrated without it), or a held-out family (`family_holdout`: the family
removed from training and calibration, and a split of in-distribution parts
held out as well). Errors are those of dz on the commanded part's nodes:

* `per_part` - rms, mae, max_abs, p95_abs, bias and the share within a
  tolerance (0.2 mm by default) of the dz error [m]; the rms of dz and the
  ratio; interval coverage and mean width at the report level; the mean
  predicted std; the envelope verdict, score and reasons;
* `per_region` - the same pooled over rim, wall and base nodes;
* `calibration` - coverage against nominal at 50-99 %, pooled, per-part
  minimum and median, mean width, and the data source the calibrator used;
* `ood` - the share of parts flagged, per split and family: a held-out family
  should be flagged more often than in-distribution parts.

Every frame has a `data_source` column.

## Compensation on the surrogate

`surrogate_compensate(target, setup, surrogate)` runs
`precomp.compensation.displacement_adjustment` - with its flange hold, z <= 0
and 65 deg formability projection - one iteration at a time on
formed(c) = c + dz_hat(c), stopping at `tolerance` (predicted RMS over the
part) or when an iteration improves the predicted error by less than 2 %
(`stagnation`), and keeps the best predicted iterate. It returns the
compensated surface, the predicted formed surface, the predicted residual
(formed - target) with its conformal interval, the envelope report of the
compensated shape and of the target, and the history - all predictions,
labelled with the model's data source.

**A prediction is not a result.** `verify_with_fea(result, setup, work_dir)`
forms the compensated surface - and the target, for the improvement factor -
with sparlab_form through the run cache and reports the achieved vertical and
normal deviation per region (all, part, wall, flange), the solver version and
deck hash, the predicted and achieved RMS, and the improvement factor
(uncompensated / compensated RMS vertical deviation over the part).
`verify_with_simulator` does the same with any simulator; with the proxy it is
labelled `proxy - not physics`. With the test-double solver the check reports
an achieved error far above the predicted one, as it must when the model was
trained on something else.

## Choosing the next runs

`rank_candidates(surrogate, candidates, n_select)` scores candidates by the
mean predictive std over the part and the envelope's part score (each
normalised over the candidates, weights 1 and 1), then picks greedily
maximising acquisition + diversity x d / median(d), d the distance in
standardised part-descriptor space to the nearest training part or earlier
pick. `precomp active` draws candidates from the design of experiments.

## Command line and API

```bash
# data (SparLab; --simulator proxy for smoke runs, announced as not physics)
precomp dataset generate --out data --simulator sparlab --work-dir runs \
    --materials AA5754-O --n-per-family 8 --seed 1 --grid-spacing 2e-3 --setup setup.json
precomp dataset info --data data
# train, calibrate, save (grouped test and calibration splits in split.json)
precomp train --data data --model gbm --out models/gbm --param n_members=8 max_iter=300
precomp train --data data --model residual --prior fea --prior-work-dir runs \
    --prior-set element_size=5e-3 layers=1 --out models/hybrid
precomp train --data scans --model transfer --base models/gbm --out models/transfer
precomp train --data data --model gbm --holdout-family saddle --out models/no_saddle
# evaluate (CSV tables with their data source), rank the next runs
precomp evaluate --model models/no_saddle --data data --family saddle --out eval
precomp active --model models/gbm --n-candidates 64 --select 8 --out next.csv
# compensate on the surrogate and verify with one simulation
precomp compensate --target cone.npz --setup setup.json --method surrogate \
    --model models/gbm --iterations 8 --verify-fea --work-dir runs --out comp
```

```python
from precomp import api
from precomp.ml import load_model, surrogate_compensate, verify_with_fea

model = load_model("models/gbm")
pred = api.predict(target, setup, model, method="surrogate")        # or the directory
res = surrogate_compensate(target, setup, model, iterations=8, tolerance=1e-4)
res.ood["in_envelope"], res.residual_interval                       # trust, and how much
check = verify_with_fea(res, setup, "runs")                         # the honest check
```

## Results on proxy data

**All numbers in this section are measured on ProxySimulator data - proxy,
NOT physics.** Script: `python/scripts/ml_proxy_demo.py`; tables:
[benchmarks/ml_proxy](../benchmarks/ml_proxy/README.md).

Design: 56 parts (8 per family, Sobol, seed 2026), AA5754-O, the default
process ranges, three variants each (168 samples); `saddle` held out; of the
other parts 30 train, 8 calibrate, 10 test, grouped by part. Relative error =
RMS dz error over the part / RMS dz of the part, averaged over parts.

| Model | test: rel. error (RMS error) | held-out saddle | coverage 90 %, expected over calibration draws | one split |
|-------|------:|------:|------:|------:|
| GBMEnsemble | 8.7 % (0.071 mm) | 7.4 % | 89.5 +- 6.1 % | 93.2 % |
| MLPEnsemble | 3.4 % (0.026 mm) | 6.0 % | 89.0 +- 5.6 % | 97.6 % |
| FieldUNet (40 x 40) | 11.3 % (0.086 mm) | 12.7 % | 88.4 +- 7.0 % | 97.5 % |

* **Envelope.** No in-distribution test part flagged; the held-out saddle not
  flagged either - it lies inside the training envelope and is predicted as
  well as the in-distribution parts. With `freeform` held out instead (a
  second run, added after that result), 20 of 24 freeform samples are
  flagged (reasons: mean wall angle, depth, radius of gyration) and predicted
  worse (13.3 % against 7.7 %). Test targets with another material (DC04) or
  a thickness outside the trained range (2.5 mm) are flagged 10 of 10, with
  the material or thickness named as the reason.
* **Compensation.** GBM surrogate DA on the 10 uncompensated test targets,
  each compensated shape formed by the proxy: RMS deviation over the part
  from 0.78 mm (target formed as is) to 0.073 mm; improvement factor median
  10.8, range 5.3-41 (freeform run: 14.0, 5.0-36). The predicted residual is
  below 0.2 um - the achieved deviation is the model's error, which only the
  verification run shows.
* **Coverage caveat.** One calibration split of 8 parts over-covered by 4-9 %
  here; the worst test part was covered at 3 % (its error is a part-wide
  offset). The expected coverage is the number to quote.
* **Cost** (two threads, shared machine): training 1.3-2.5 min per point model
  on 180 000 rows, 1-1.5 min for the U-Net; featurising a 101 x 101 grid with
  its contour tool path 0.07-0.15 s; one surrogate DA iteration (featurise, predict,
  condition) about 1 s.

## Limitations

* **No physics yet.** Every number measured so far is on proxy data or the
  test double. Accuracy, coverage and the compensation factor on SparLab data
  are unknown until the solver lands and a design is simulated; the proxy is
  smooth and noise-free and probably easier than FE data in some respects
  (no discretisation noise, no path-dependent local effects) and harder in
  others (strong multiplicative process dependence by construction).
* **Small-data regime.** A few dozen parts; the part, not the node, is the
  unit of evidence. Errors are dominated by a part-wide offset, the ensemble
  spread tracks them only loosely (correlation ~0.3 between a part's mean std
  and its error on proxy data), and a single conformal calibration of a dozen
  parts can land 6-8 % from its expected coverage; per-part coverage ranges
  from 0 to 1. Report the expected coverage (`partition_coverage`) and the
  per-part minimum, not only the pooled number.
* **Point models see global effects only through descriptors.** The
  unclamping springback of the whole sheet is non-local; the GBM and MLP learn
  it from the part descriptors and the placement features, the U-Net from the
  whole map but with a coarse training grid (40 x 40 in the demonstration).
* **Envelope thresholds** are 99 % quantiles over a few dozen held-out parts,
  i.e. close to the largest distance seen; 0-8 % of in-distribution parts
  were flagged in the proxy runs. A new family is flagged only if its
  descriptors differ (freeform: 83 %; saddle, close to the pyramids: 0 %). The reasons are indicative (range excess, Mahalanobis
  contributions), not causal.
* **Extrapolation off the training mask.** Point models are trained on part
  nodes; their prediction on the flange is an extrapolation (displacement
  adjustment holds the flange and measures on the part, so it does not use
  it).
* **Transfer learning** has been tested only on a shifted proxy; the Bayesian
  ridge correction is linear in the features and the base prediction.
* **Hybrid prior.** `FEAPrior` needs sparlab_form at prediction time and its
  cost per call (one coarse run per DA iteration).
* **Pickles.** Bundles use joblib; load only trusted bundles.
* **Tool path.** The `time_frac` feature rebuilds the setup's tool path;
  if the tool-path code changes, features of old samples change with it
  (the schema version does not track it). Record paths with
  `generate(..., store_toolpath=True)` when that matters.
