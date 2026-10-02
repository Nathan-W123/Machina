# precomp.ml: learned springback prediction and surrogate pre-compensation

`precomp.ml` predicts, before the first forming attempt, the surface a
commanded tool path will form, and inverts that prediction by displacement
adjustment so that the first part lands within tolerance - and says when the
prediction should not be trusted. It is built on the `precomp` package
([docs/precomp.md](precomp.md)): its height maps, part families, tool paths,
the SparLab bridge, metrology and compensation are used as they are.

**Where the numbers come from.** Every metric the package reports carries its
data source: `SparLab simulation`, `scan`, or `proxy - not physics`. The
benchmarks on the C++ forming solver `sparlab_form` are summarised in
[Results on SparLab simulations](#results-on-sparlab-simulations): they are
simulations with nominal material data, not experiments, and no model here
has seen a scanned part. The demonstration in
[Results on proxy data](#results-on-proxy-data) and the unit tests use the
`ProxySimulator`, an analytic stand-in that is **not physics**; those numbers
show that the pipeline works and what its outputs look like, nothing more.

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
* [Results on SparLab simulations](#results-on-sparlab-simulations)
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
* **One fixture and machine setup.** The blank size and clamp frame enter
  through the placement features; the release (3-2-1 or clamped), the tool
  path style, spacing and direction, the mesh, the element, the increment
  and the solver settings enter through **no** feature. A model records the
  values its training setups had (`training["setup"]`, `setup_envelope`) and
  its envelope flags a query with any other value (`setup.<field>` in the
  reasons): with the proxy, the release alone changes dz by a factor of
  2.7, which no feature would show.
* **Grid.** The grid spacing must not exceed the ring width
  (4 tool radii / 4 rings = one tool radius by default); refused otherwise.
* **Parts centred on the blank**, as `precomp` builds them.
* **Tool path.** The `time_frac` feature uses the setup's own tool path
  (`precomp.fea.make_toolpath`); a spiral path needs one loop per level, so
  multi-pocket parts need `toolpath_style="contour"`. The deck builder refuses
  a commanded surface without a path, and so does the `ProxySimulator`: the
  job fails at generation with that reason. With the default spiral, 36 of
  the 147 samples of the test design (freeform dents, perturbed pyramids)
  fail; with contour paths none do.

## Modules

| Module | Public API |
|--------|------------|
| `features` | `FEATURE_SCHEMA_VERSION`, `FeatureSpec(name, unit, description, kind)`, `FeatureConfig(scales, rings, ring_max, time_source, part_eps, rim_width)` with `.specs()`, `.names`, `.schema_hash()`; `point_features(commanded, setup, toolpath=None, points=None, *, config) -> (X, names)`; `feature_maps(...) -> FeatureMaps`; `global_features(commanded, setup, config) -> (d, names)`; `region_labels(...)` |
| `dataset` | `Sample(sample_id, commanded, formed, setup, source, fidelity, kind, part, part_id, target, toolpath, provenance)`; `Dataset.create(root, created_at=...)` / `Dataset(root)`: `append`, `load`, `ids`, `index`, `filter`, `record_failure`, `failures`, `to_table(features, target="dz", points_per_sample, mask, rng, ids, stratify, prior, n_jobs) -> Table` (unpacks as `X, y, groups`); `build_table(samples, ...)`; `grouped_split`, `grouped_kfold`, `family_split`, `check_disjoint`; `SOURCE_LABELS`, `source_label` |
| `generate` | `DesignSpace`, `design_points(space, n_per_family, seed) -> [DesignPoint]`, `PerturbationSpec`, `perturbed_commanded`, `compensated_commanded`, `generate(dataset, points, simulator, *, created_at, kinds, seed, compensator) -> GenerationReport`, `simulate_samples`; `SparlabSimulator(work_dir, max_workers)`, `ProxySimulator(params, check_toolpath=True)`, `ProxyParams`, the `Simulator` protocol |
| `models` | `GBMEnsemble`, `MLPEnsemble`, `FieldUNet`, `ResidualModel(prior, learner)`, `TransferModel(base, correction="linear", prior_gain, prior_offset, min_parts)`, `FEAPrior(work_dir, overrides)`, `model_from_state`, `thread_limit`; protocols `DeviationModel` (point: `fit(X, y, groups)`, `predict(X) -> (mean, std)`, `feature_names`, `to_state`/`from_state`) and `FieldDeviationModel` (`fit_samples`, `predict_field(commanded, setup)`) |
| `uncertainty` | `ConformalCalibrator(eps, eps_fraction, mondrian)`: `fit(mu, std, y, part_ids, region)`, `quantile(level)`, `supports(level)`, `intervals(mu, std, level)`, `baseline_halfwidth(level)`; `min_parts(level)`, `intervals(...)`, `coverage_report(...)`, `partition_coverage(...)` |
| `ood` | `OODEnvelope(quantile, k, ...)`: `fit(point_X, point_groups, part_D, part_groups, ...)`, `assess_features(point_X, part_d) -> {in_envelope, part_score, fraction_points_out, reasons, ...}` |
| `surrogate` | `DeviationSurrogate(model, features, data_source, training, calibrator, ood, domain)`: `target` ("dz"), `predict_deviation(commanded, setup) -> (mean, std)`, `predict_interval(commanded, setup, level)`, `assess(commanded, setup)`; `train_surrogate(model, train, calibration=..., features, points_per_sample, mask, seed, group_by)`, `transfer_surrogate(base, real, calibration, correction)`, `calibrate`, `fit_envelope`, `setup_envelope`, `setup_mismatch`, `seen_ids` |
| `registry` | `save_model(model, directory, manifest)`, `load_model(directory)`, `read_manifest`, `check_schema` |
| `evaluate` | `evaluate_surrogate(surrogate, samples, tolerance, level) -> EvaluationReport` (`per_part`, `per_region`, `calibration`, `ood` frames, `summary()`, `to_csv()`); `cross_validate(samples, make_model, n_splits)`; `family_holdout(samples, family, make_model)` |
| `compensate` | `SurrogatePredictor(model, setup)`, `surrogate_compensate(target, setup, surrogate, iterations, tolerance, stagnation, level, allow_out_of_envelope, interval_stop, ...) -> SurrogateCompensation`, `verify_with_fea(result, setup, work_dir)`, `verify_with_simulator(result, setup, simulator)` |
| `active` | `Candidate(commanded, setup, label)`, `rank_candidates(surrogate, candidates, n_select, weight_std, weight_ood, diversity) -> DataFrame` |
| `cli` | `main(argv)`: `precomp dataset generate|info`, `precomp train`, `precomp evaluate`, `precomp active` |

`DeviationSurrogate` implements the `FieldModel` protocol of
`precomp.compensation` (plus `predict_interval` and `assess`), so
`precomp.api.predict` / `compensate` and `precomp compensate --method
surrogate` use it without knowing about machine learning; the API and the
command line also accept a bundle directory. It predicts the **total**
deviation dz (`target = "dz"`) - also a `ResidualModel`, which evaluates its
prior itself - so `method="hybrid"`, which adds the model to a simulation,
refuses it: the springback would be counted twice.

## Features

Per grid node, `FEATURE_SCHEMA_VERSION = "2"`, 57 features in five kinds
(`FeatureConfig().specs()` gives each with its unit and description):

| Kind | Features | Units |
|------|----------|-------|
| local (invariant when the part moves in x-y) | `depth`, `depth_frac`; `wall_angle_s1..3` (Gaussian-smoothed at 1, 2.5, 6 tool radii); `mean_curv_s1..3` = H sigma and `gauss_curv_s1..3` = K sigma^2 (the package's moment-normalised derivatives, dimensionless); `rim_dist` (signed distance to the rim, + inside); `radial_norm` (distance to the part centroid / equivalent radius sqrt(A/pi)); `polar_sin`, `polar_cos` (regularised by the tool radius R: dy / sqrt(r^2 + R^2) - the angle is undefined at the centroid); `time_frac` (pseudo-time of the tool path where the tool passes nearest: the node offset by R along its normal, nearest in-contact path point; or the depth fraction with `time_source="depth"`); `ring{1..4}_{mean,min,max}` (height of the neighbours in four annuli out to 4 R, relative to the node) | m, -, rad |
| placement | `clamp_dist` (L-infinity distance to the clamped frame), `blank_radius` (distance to the blank centre / free half-width) | m, - |
| global (per part) | `max_depth`, `area`, `volume`, `mean_wall_angle`, `max_wall_angle`, `perimeter_over_area` (rim contour length), `aspect` (principal axes of the plan area), `depth_gyration_major/minor` (radii of gyration of the depth field) | m, m^2, m^3, rad, 1/m, - |
| process | `tool_radius`, `step_down`, `thickness`, `friction`, `elastic_curvature` = 3 sigma_f(0.2) / (E t), the unloading curvature of a fully plastic bend | m, -, 1/m |
| material | `youngs_modulus`, `poisson_ratio`, `yield_stress`, `flow_stress_20`, `flow_stress_50` (isotropic flow stress at plastic strain 0.2, 0.5), `hardening_modulus`, `voce_Q`, `voce_rate`, `r0`, `r45`, `r90` (1 when isotropic), `backstress_sat` (monotonic back stress at strain 0.5: Prager + Armstrong-Frederick, tending to sum C / gamma), `flow_stress_20_over_E` (sigma_f(0.2) / E) | Pa, - |

sigma_f(0.2) is the flow stress at 20 % plastic strain - a representative
stress of an incrementally formed sheet - not the 0.2 % proof stress that
"sigma_0.2" usually denotes. Schema version 2 renamed `yield_over_E`, which
was never the yield stress over E, to `flow_stress_20_over_E`; version-1
bundles are refused.

`FeatureConfig.schema_hash()` fingerprints the version, names, units and
configuration; a bundle stores it and refuses to load under another. Every
feature is checked finite. The ring minima and maxima are exact footprint
filters computed as 1-D running extremes over the runs of each footprint
row, so their cost grows with the ring radius in cells, not its area (equal
to `scipy.ndimage.minimum_filter` bit for bit). Measured: a 200 x 200 grid
at 1 mm featurises in 0.25-0.6 s (0.9 s including the tool path) - the
specification asked for 2 s; at 0.5 mm with an 8 mm tool (32 cells to the
last ring) 0.7 s, where the footprint filters took 16.8 s. Moving a pyramid
by (7, -4) grid spacings leaves every non-placement feature equal to 1e-9 of
its range at the corresponding nodes with `time_source="depth"`; with the
tool path as time source a spiral gives the same `time_frac`, a contour path
differs at under 1 % of the nodes by up to 0.03 (where its loops start;
`test_ml_features`).

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
(`scipy.stats.qmc`, seeded by the seed and a hash of the family's **name**)
over the family's parameter bounds, the process ranges (tool radius 4-8 mm,
step-down 0.2-1 mm, thickness 0.6-1.5 mm, friction 0.05-0.2 by default) and
the material index; draws the family refuses, or that do not fit the
unclamped window with the tool, are skipped and the sequence continues, so a
seed gives the same design every time, and a point id (`dome-s1-0003`) names
the same part whatever other families are drawn and however many points per
family. Freeform parts are drawn by `Freeform.sample` from a seed taken from
the Sobol point.

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
Before skipping, it checks that a stored sample was made from the same design
point (part, setup physics, grid); a data set that holds another part under
the same id - a design drawn from another space under the same seed - is
refused before anything runs. The perturbation of a point is seeded by its
point id, not by its position in the list.
Generating the proxy data set also exposed a stall in
`precomp.compensation.limit_wall_angle`, fixed on this branch: the cone slope
is now bisected towards tan(limit) / sqrt(2) when reducing it by the measured
excess does not move the measure.

**The ProxySimulator** (`generate.ProxySimulator`) - NOT physics. With d the
commanded depth, theta the wall angle, s the signed rim distance and
sigma_f(0.2) / E, t, R, step-down and friction from the setup:

* rim under-forming (dz > 0, too shallow): A sin(theta at 2.5 R) exp(-s / 2.5 R)
  inside the rim and exp(-(s / R)^2) outside, A = 0.8 mm
  (sigma_f(0.2) / E / 3.8e-3)^0.75 (1 mm / t)^0.5 (1 + (mu - 0.1)) (1 - e^(-D / 10 mm));
* pillow on flat bases (theta < 12 deg): P (1 - exp(-(e / 12 mm)^2)), e the
  distance from the flat region's edge, P = 0.5 mm (sigma_f(0.2) / E / 3.8e-3)^0.5
  sqrt(R / 5 mm) sqrt(step-down / 0.5 mm) (1 mm / t)^0.25 min(D / 10 mm, 1);
* both smoothed at R / 2, plus global bending after the 3-2-1 release,
  (kx (a^2 - x^2) + ky (a^2 - y^2)) / 2, zero at the supports x = y = a,
  kx + ky = 2 kappa, kappa = 0.1 /m (sigma_f(0.2) / E / 3.8e-3)(D / 25 mm)(1 mm / t),
  shared between x and y by the second moments of the depth field (a quarter
  of it, inside the window, with `release="clamped_only"`).

It gives dz of 0.15-2.4 mm on the default five-material design, depends on the commanded
shape (so DA on it is not trivial), and has local and global structure. It is
labelled `proxy - not physics` in every sample, table, manifest and report.
It knows single-point forming only: as a simulator or a model it refuses a
setup with support from below (`FormingSetup.support`; a failed job for
`run`) rather than give single-point numbers under the supported setup's
name; as a `ResidualModel` prior it gives single-point dz for any setup.

## Models

All models are deterministic for a seed (bootstraps, initialisations,
shuffles and dropout masks derive from it) and cap their threads at
`n_threads` or `$PRECOMP_ML_THREADS` (scikit-learn/OpenMP via threadpoolctl,
and torch). scikit-learn is imported with the module, not lazily:
threadpoolctl caps only runtimes already loaded, and a lazy import inside the
limited block ran unlimited - twenty times slower on a loaded machine.

**Target scale.** The point models and the U-Net learn dz divided by the
`elastic_curvature` feature 3 sigma_f(0.2) / (E t) (option `target_scale`, stored
in the state; `None` learns dz as is). Springback scales with it to first
order; dividing it out spares the model building the product of material,
thickness and depth from a few dozen parts. On a proxy design with one
material and five process parameters the held-out error fell from 0.14 to
0.09 mm; across five materials from 0.17 to 0.09 mm.

* **`GBMEnsemble`** - M (8) `HistGradientBoostingRegressor` members (300
  iterations, 31 leaves, 40 samples per leaf, 80 % of the features per split),
  each on a bootstrap of whole samples; mean and standard deviation over the
  members; optional quantile members for a 5/95 % band. The three variants
  of a part are near duplicates, so a member sees about 95 % of the parts
  and the spread is far below the error on a new part (0.020 against
  0.055 mm on proxy data): it is a relative scale for the conformal
  calibration, not an error estimate. `train_surrogate(group_by="part")`
  bootstraps parts instead; on proxy data that doubled the spread and raised
  its correlation with a part's error (0.57 to 0.70), but raised the held-out
  error by 18 % and, after calibration over parts, widened the 90 % intervals
  from 0.35 to 0.65 mm at the same expected coverage - so samples stay the
  default.
* **`MLPEnsemble`** (torch) - M (5) MLPs (128, 128, SiLU) with a
  (mean, log-variance) head trained on the Gaussian NLL after a few epochs of
  squared error, Adam, early stopping on a validation split of whole samples
  (with `group_by="part"`, of whole parts: the validation loss is then
  honest, but the held-out error was worse in a review experiment, 0.034
  against 0.023 mm, because the validation parts are lost to training); the
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
  executable whenever evaluated) or a closed-form estimate. A prior that
  takes `target=` (`FEAPrior` does) gets the part a support's fixture is
  made for - each training sample's target, and at prediction the target
  given to `predict_deviation(commanded, setup, target)` (as
  `SurrogatePredictor`, `precomp.api.predict` and `compensate` pass it) - so
  it simulates a compensated command on the same backing plate or rim pass
  band as the label runs; without it, on the command's own outline. The prior - its
  class, settings and data source, e.g. `proxy - not physics` for
  `--prior proxy` - is recorded in `training["prior"]`, the manifest and
  `describe()`: a model trained on SparLab data around a proxy prior says so.
* **`TransferModel(base, correction)`** - the simulation-trained base frozen,
  plus a low-dimensional correction fitted on a few real parts:
  `"linear"` (default) y = a + (1 + b) mu, a gain and an offset on the base
  mean, or `"gain"` alone - Bayesian regressions with a fixed prior
  (b ~ N(0, 0.5^2), a ~ N(0, (0.5 s)^2), s the RMS base prediction) in which
  every real part weighs one observation, the parameter uncertainty combined
  with the base std in quadrature; or `"gbm"`, a small, heavily regularised
  GBM on [features, base mean, base std]. A correction is kept only if its
  leave-one-real-part-out error beats the base's; otherwise - and always
  with fewer than two real parts - the model predicts as the base
  (`accepted`, `validation`). It reports `n_real` (parts); with none it is
  the base, bit for bit. The former default, a Bayesian ridge on all ~60
  standardised features, made the base 2-5 times worse from 1-3 real parts
  on proxy data, also when the real data equalled the simulation; in the
  same experiment the validated linear correction was never accepted without
  a shift (error = the base's) and cut a 30-40 % shift from 0.18 to
  0.09-0.15 mm from two parts on. `transfer_surrogate` keeps the base's
  envelope (widened to the real setups), records the real and calibration
  parts and families with the base's, and gives intervals only when
  calibrated on held-out real parts - never with a calibration from another
  data source.

## Uncertainty: conformal intervals

The ensemble spread or predicted variance is not calibrated by itself.
`ConformalCalibrator` is fitted on held-out **whole parts** - never on parts
the model was trained on (`train_surrogate` refuses an overlap) - and counts
**parts**, not nodes: the errors of one part are one smooth field, and the
three variants of a design point are one part (`calibrate` groups by part id).

    score r = |y - mu| / (std + eps),   eps = 0.1 median(std) of the calibration set,
    q = the smallest t with  sum_j F_j(t) >= level (k + 1),
    interval = mu -+ q (std + eps),

F_j the empirical distribution of the scores of calibration part j (each
part weighs one, whatever its node count) and k the number of calibration
parts.

**What the guarantee is.** For a new part exchangeable with the calibration
parts, the **expected share of its nodes inside the interval is at least the
level**. (Had the new part's F been included with the same weight, the
threshold would be symmetric in all k + 1 parts and the new part's expected
coverage would equal the average, which is >= level by construction; leaving
it out - counting its F as 0 - can only raise q.) With one node per part this
is ordinary split conformal prediction. A level needs k >= level / (1 - level)
parts - 9 at 90 %, 19 at 95 %, 99 at 99 % (`min_parts`); asking more is an
error, and tables report such a level as not supported. The price of
validity with few parts is conservatism: when parts behave alike the
expected coverage approaches level (k + 1) / k (97.5 % at a nominal 90 % with
12 parts).

**What it is not.** A guarantee for any one part: a part whose error is a
large offset can be covered far less. The evaluation reports, per level, the
mean of the per-part coverages (the guaranteed quantity), the worst, 10th
percentile and median part, and the share of parts covered at the level.
Nor does validity show that the std is useful - any std, even a random one,
gives valid intervals. The calibrator therefore also calibrates a
**constant-width** interval mu -+ q0 on the same parts, and the tables report
the width ratio: below 1 the std narrows the intervals, at or above 1 it
carries no information about where the errors are. On proxy data it
depends on the data: 0.74-0.87 for the three models of the demonstration
(expected over calibration draws; the std narrows the 90 % intervals by
13-26 %), but 1.06 for the GBM of its freeform run and 1.2-1.4 for the
smaller GBM of the tests, where the std does not narrow the intervals and
only moves width to the parts with the largest errors (the median worst-part coverage at 90 % rises from 0.28 to 0.84); a
std permuted at random over the nodes gives wider intervals still.
`partition_coverage` re-splits held-out parts many times into calibration
and test parts and reports the expected coverage, the share of parts at the
level, the worst part and the width ratio. `mondrian=True` takes one
quantile per region (rim / wall / base), each over the parts that have it.

(The first version counted nodes in the finite-sample correction - about
70 000 of them, from 8 parts - so the correction was void, and it grouped
the calibration by sample, so a part's three variants counted three times.
Its expected coverage over calibration draws was about nominal only because
exchangeability makes that hold for any std; one draw of 9 parts
under-covered, below 85 %, a third of the time.)

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
third of in-distribution parts; the folds brought that to 0-10 %. `assess`
returns `in_envelope`, `part_score`, `fraction_points_out`, the thresholds,
and `reasons`: the features outside their training range (largest excess
first), then the largest contributions to the Mahalanobis distance.

**Setup fields.** Every `FormingSetup` field that no feature describes -
blank and clamp, mesh, element, tool path style, spacing and direction,
contact, increment, release, kinematics, solver, the rim support and its
settings (`support`, `support_settings`; `docs/precomp.md`, "Rim support") -
is recorded from the training setups (a range for numbers, the set of
values otherwise); a query outside puts the part outside the envelope, with
`setup.<field>` first among the reasons and the details in `setup_mismatch`.
A model trained on single-point forming therefore refuses a backing plate
or DSIF, and one trained on a support refuses another; a support setting
that spells out its default counts as that default. A field the model
does not record at all was added after it was trained: its runs had the
field's default, so any other value is outside (the mismatch says so). (A
setup whose tool path cannot form the part at all - a spiral on several
pockets - raises instead: the features need the path.)

**What it detects.** New **descriptors**, not new family names. Leaving each
family out in turn (see [Results](#results-on-proxy-data)), only freeform -
whose dents no other family has - is flagged; a dome, a saddle or a two-level
part is inside the envelope of the other families and, on proxy data,
predicted about as well as they are. A family that is flagged is predicted
worse; the envelope says nothing about a part that is hard to predict for
reasons its descriptors do not show.

## Model bundles

`save_model(surrogate, dir, manifest)` writes `manifest.json`,
`model.joblib` (the state; scikit-learn estimators pickled) and
`torch/<name>.pt` (state dicts, tensors only). The manifest holds the model
class, what it predicts (`target`: "dz"), `feature_schema_version`,
`feature_schema_hash`, the feature names and their SHA-256, the feature
configuration, the data source of the training data, the training record
(sample and part ids up to 5000 and always a SHA-256 and count, calibration
parts, families, sources, fidelities, SparLab versions, the setup fields
seen, grid spacings, the prior of a residual model, the real parts and
validation of a transfer model, seed, sizes), metrics (each block must state
its `data_source`; a metric that could not be computed is null), the
conformal quantiles at 50/80/90/95/99 % (null where the calibration parts
cannot support the level) with the constant-width baseline, the OOD
thresholds, package versions, the git revision (`-dirty` when the tree has
changes), `created_at` given by the caller, notes, and the SHA-256 of every
file. Everything that can fail on the content is checked before a file is
written, and the manifest is written last, so an interrupted save leaves no
loadable bundle. `load_model` refuses another schema version, another schema
hash for the same configuration, feature names that do not match their hash,
a manifest without the hash of `model.joblib` or of a torch file the state
refers to, any file that does not match its hash, and a state whose model
class, kind, domain or data source differ from the manifest's - with a
message saying which and what to do; a neural bundle without torch says that
torch is missing. **`model.joblib` is a pickle: load bundles you trust
only.**

## Evaluation protocol

A model is scored only on whole parts it has not seen (`evaluate_surrogate`
refuses a sample whose sample id **or part id** was used to train, correct -
the real parts of a transfer model - or calibrate it, the base's included):
a held-out test
split, grouped K-fold cross-validation (`cross_validate`, each fold trained
and calibrated without it), or a held-out family (`family_holdout`: the family
removed from training and calibration, and a split of in-distribution parts
held out as well). Errors are those of dz on the commanded part's nodes:

* `per_part` - rms, mae, max_abs, p95_abs, bias and the share within a
  tolerance (0.2 mm by default) of the dz error [m]; the rms of dz and the
  ratio; interval coverage and mean width at the report level; the mean
  predicted std; the envelope verdict, score and reasons;
* `per_region` - the same pooled over rim, wall and base nodes;
* `calibration` - per level (50-99 %): whether the calibration parts support
  it; the mean per-part coverage (the guaranteed quantity) and the pooled
  node coverage; the worst, 10th-percentile and median part and the share of
  parts covered at the level (parts are part ids, the variants of a design
  point pooled); the mean width, the constant-width baseline's width and
  coverage and the width ratio; the number of test and calibration parts and
  the data source the calibrator used;
* `ood` - the share of parts flagged, per split and family: in-distribution
  parts should rarely be; a held-out family is flagged when its descriptors
  are new to the model.

Every frame has a `data_source` column.

## Compensation on the surrogate

`surrogate_compensate(target, setup, surrogate)` runs
`precomp.compensation.displacement_adjustment` - with its flange hold, z <= 0
and 65 deg formability projection - one iteration at a time on
formed(c) = c + dz_hat(c), and keeps the best predicted iterate. It stops at
`tolerance` (predicted RMS over the part) or when an iteration improves the
predicted error by less than 2 % (`stagnation`) - by then the model's own
map is inverted and the predicted residual is near zero, which says nothing
about the achieved one. `interval_stop=True` stops instead once the
predicted RMS error is within the RMS half-width of the model's interval at
`level` - the model cannot resolve less. That saves two thirds of the
predictions, but on proxy data the verified deviation was 4 % smaller on the
23 held-out targets of the tests, 9 % larger on the 11 of the
demonstration and about the same in its freeform run, so it is not the
default. It returns the compensated
surface, the predicted formed surface, the predicted residual (formed -
target) with its conformal interval and mean half-width, the envelope report
of the compensated shape and of the target, and the history - all
predictions, labelled with the model's data source.

**The envelope is enforced.** A target outside the training envelope is
refused (PrecompError naming the reasons), and every iterate is assessed
before the model is evaluated on it: the first one outside ends the loop
(`stopped = "left_envelope"`), keeping the best iterate inside.
`allow_out_of_envelope=True` overrides both and records the verdicts.
`precomp.api.compensate` and `precomp compensate` refuse a target or a
compensated shape outside the envelope unless `allow_out_of_envelope` /
`--allow-out-of-envelope`; the command prints the verdict, puts it in the
report, and reports the **verified** formed surface (the simulation of the
compensated shape) when `--verify-fea` ran, the prediction - labelled as such
- only otherwise. The API and the command line run `--iterations` DA steps;
the stopping rules above are those of `surrogate_compensate`.

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
pick. `precomp active` draws candidates from the design of experiments in the
model's **training setup** by default - its materials, the process ranges it
saw (a fixed process stays fixed), every other setup field as trained, its
grid spacing - so the candidates differ in geometry; `--materials`,
`--process KEY=LO,HI`, `--setup` and `--grid-spacing` override.

## Command line and API

```bash
# data (SparLab; --simulator proxy for smoke runs, announced as not physics)
precomp dataset generate --out data --simulator sparlab --work-dir runs \
    --materials AA5754-O --n-per-family 8 --seed 1 --grid-spacing 2e-3 --setup setup.json
precomp dataset info --data data
# train, calibrate, save (grouped test and calibration splits in split.json)
precomp train --data data --model gbm --out models/gbm --param n_members=8 max_iter=300
precomp train --data data --model residual --prior fea --prior-work-dir runs \
    --prior-set element_size=5e-3 layers=1 --out models/residual   # use with --method surrogate
precomp train --data scans --model transfer --base models/gbm --out models/transfer \
    --param correction=linear             # kept only if it beats the base on held-out real parts
precomp train --data data --model gbm --holdout-family saddle --out models/no_saddle
# evaluate (CSV tables with their data source), rank the next runs
precomp evaluate --model models/no_saddle --data data --family saddle --out eval
precomp active --model models/gbm --n-candidates 64 --select 8 --out next.csv
# compensate on the surrogate and verify with one simulation (refused outside
# the training envelope unless --allow-out-of-envelope)
precomp compensate --target cone.npz --setup setup.json --method surrogate \
    --model models/gbm --iterations 8 --verify-fea --work-dir runs --out comp
```

```python
from precomp import api
from precomp.ml import load_model, surrogate_compensate, verify_with_fea

model = load_model("models/gbm")
pred = api.predict(target, setup, model, method="surrogate")        # or the directory
res = surrogate_compensate(target, setup, model, iterations=8)     # refused outside the envelope
res.stopped, res.residual_interval                                  # why it stopped, how sure
check = verify_with_fea(res, setup, "runs")                         # the honest check
```

## Results on SparLab simulations

**SparLab simulations (`sparlab_form`, AA5754-O with nominal handbook data),
not experiments.** Two benchmarks compare the ML first shot (surrogate DA on
the target, then one simulation of the compensated shape) with forming the
target as it is and with FE displacement adjustment (FE-DA-k: k steps with
the FE model in the loop, k + 1 runs per part), on the same 8 held-out test
parts (two per family: truncated cone, pyramid, dome, elliptic cone), stage
n18 (48 training / 16 calibration parts). Vertical RMS deviation over the
part, mean over the 8 parts [mm]:

| benchmark | uncompensated | FE-DA-1 | FE-DA-2 | ML-GBM | ML-MLP | dz error GBM / MLP |
|---|--:|--:|--:|--:|--:|--:|
| [2 mm mesh, single-point](../benchmarks/springback/README.md) | 0.731 | 0.623 | 0.602 | 0.634 | 0.626 | 0.085 / 0.106 |
| [fine mesh, backing plate](../benchmarks/springback_fine_ml/README.md) (data at penalty 3, test at 10) | 0.300 | 0.263 | 0.272 | 0.253 | 0.250 | 0.049 / 0.052 |

On both, the ML first shot is about as good as one FE-DA step with one FE
run fewer per part; on the plate it is 0.010-0.012 mm better on average
(5 / 6 of 8 parts), which 8 parts do not resolve (paired sd 0.017-0.018 mm).
The remaining deviation is dominated by the rim sag, which no method
removes. The training set is 144 simulated samples on the 2 mm mesh and
128 (about 43 CPU-h) on the fine mesh. Coverage of the 90 % conformal intervals
on the fine-mesh test samples: 95 % (GBM), 88 % (MLP); one of the 8 test
targets was outside the training envelope and compensated with the
override. Details, per-part tables and limitations are in the two READMEs.

## Results on proxy data

**All numbers in this section are measured on ProxySimulator data - proxy,
NOT physics.** Script: `python/scripts/ml_proxy_demo.py`; tables:
[benchmarks/ml_proxy](../benchmarks/ml_proxy/README.md) (protocol revision
2: revision 1 had 8 calibration parts, too few for a 90 % interval once
calibration counts parts).

Design: 70 parts (10 per family, Sobol, seed 2026), AA5754-O, the default
process ranges, contour tool paths, three variants each (210 samples);
`saddle` held out; of the other parts 36 train, 12 calibrate, 12 test,
grouped by part. Relative error = RMS dz error over the part / RMS dz of the
part, averaged over samples.

| Model | test: rel. error (RMS error) | held-out saddle | coverage per part at 90 %: expected over calibration draws / one split | width / constant width (expected) |
|-------|------:|------:|------:|------:|
| GBMEnsemble | 6.7 % (0.049 mm) | 9.3 % | 97.0 +- 3.1 % / 99.3 % | 0.80 |
| MLPEnsemble | 3.1 % (0.023 mm) | 3.9 % | 97.1 +- 2.5 % / 99.8 % | 0.74 |
| FieldUNet (40 x 40) | 9.5 % (0.071 mm) | 10.2 % | 96.7 +- 3.1 % / 99.8 % | 0.87 |

* **Coverage.** Valid and conservative, as it must be with 12 calibration
  parts (the bound is 97.5 %); every test part is covered at 96 % or more.
  95 % and 99 % are not supported by 12 parts. The std narrows the 90 %
  intervals by 13-26 % against a constant width on this data (not on the
  unit-test data).
* **Envelope.** 3 of 36 in-distribution test samples (one part) flagged; the
  held-out saddle not flagged, and predicted about as well. Leaving each
  family out in turn, only freeform is flagged (87 %, predicted at 17.6 %
  against 7.2 %); dome, elliptic cone, pyramid, saddle and two-level are not
  (0 %), truncated cones 20 %, and the elliptic cone is predicted worse
  (12.4 % against 8.7 %) without being flagged. With `freeform` held out in a
  second run, 27 of its 30 samples are flagged and it is predicted at
  17.6 % (6.8 % in distribution). Test targets with another material (DC04), a thickness outside
  the trained range (2.5 mm), another release or a spiral tool path are
  flagged 12 of 12, with the material, thickness or setup field named first.
* **Compensation.** GBM surrogate DA on the 12 uncompensated test targets:
  one refused (outside the envelope); the other 11, each compensated shape
  formed by the proxy: RMS deviation over the part from 0.69 mm (target
  formed as is) to 0.045 mm; improvement factor median 18.2, range 7.1-59
  (freeform run: 17.7, 11.6-45; two targets refused). The predicted residual is below 1.3 um with an
  interval half-width of 0.15 mm - the achieved deviation is the model's
  error, which only the verification run shows.
* **Cost** (two threads, shared machine): training 1.6-3.3 min per model on
  216 000 rows (U-Net: 108 samples), about 26 min for the whole run including
  the seven leave-family-out models; featurising a 101 x 101 grid with its
  contour tool path about 0.15 s; one surrogate DA iteration (featurise,
  predict, condition) about 1 s.

## Limitations

* **Simulation, not experiment.** The SparLab results are simulations with
  nominal material data, on 8 test parts from four families; no model has
  been trained on or checked against scanned parts. The proxy results say
  nothing about SparLab: the proxy is smooth and noise-free and probably
  easier than FE data in some respects (no discretisation noise, no
  path-dependent local effects) and harder in others (strong multiplicative
  process dependence by construction).
* **Small-data regime.** A few dozen parts; the part, not the node, is the
  unit of evidence. Errors are dominated by a part-wide offset, and the
  ensemble spread tracks them only loosely: whether it narrows the intervals
  at all depends on the data (width ratio against a constant width 0.7-0.9
  in the demonstration, 1.06 in its freeform run, 1.2-1.4 on the test
  data). The conformal guarantee
  needs at least 9 calibration parts at 90 % and is conservative with a
  dozen (expected coverage up to 97.5 %); it is a guarantee on average over
  parts, and a single part can still be covered far less. Report the
  per-part distribution and the share of parts at the level, not only the
  pooled number.
* **Point models see global effects only through descriptors.** The
  unclamping springback of the whole sheet is non-local; the GBM and MLP learn
  it from the part descriptors and the placement features, the U-Net from the
  whole map but with a coarse training grid (40 x 40 in the demonstration).
* **Envelope thresholds** are 99 % quantiles over a few dozen held-out parts,
  i.e. close to the largest distance seen; 0-9 % of in-distribution samples
  were flagged in the proxy runs. A new family is flagged only if its
  descriptors differ: of seven families left out in turn, only freeform
  (87 %); an unflagged family can still be predicted worse (the elliptic
  cone). The reasons are indicative (range excess, Mahalanobis
  contributions), not causal.
* **Extrapolation off the training mask.** Point models are trained on part
  nodes; their prediction on the flange is an extrapolation (displacement
  adjustment holds the flange and measures on the part, so it does not use
  it).
* **Transfer learning** has been tested only on proxy data. The correction is
  a gain and an offset of the base prediction: it cannot fix an error whose
  shape the base does not predict, it needs two real parts to be validated
  at all, and a leave-one-part-out check on two or three parts is itself
  noisy (the unvalidated gain alone was accepted wrongly in 2 of 12
  no-shift cases).
* **Refusals.** With 99 % thresholds from a few dozen parts, 0-10 % of
  in-distribution targets fall outside the envelope and are refused by the
  compensation until overridden (`--allow-out-of-envelope`) after reading the
  reasons.
* **Hybrid prior.** `FEAPrior` needs sparlab_form at prediction time and its
  cost per call (one coarse run per DA iteration).
* **Pickles.** Bundles use joblib; load only trusted bundles.
* **Tool path.** The `time_frac` feature rebuilds the setup's tool path;
  if the tool-path code changes, features of old samples change with it
  (the schema version does not track it). Record paths with
  `generate(..., store_toolpath=True)` when that matters.
