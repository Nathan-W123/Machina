#!/usr/bin/env python3
"""End-to-end demonstration of precomp.ml on the ProxySimulator - NOT physics.

usage:
    python3 python/scripts/ml_proxy_demo.py --out benchmarks/ml_proxy \
        [--work scratch_dir] [--n-per-family 10] [--seed 2026] [--threads 2]

Every number this script writes comes from the ProxySimulator, an analytic
stand-in for springback (rim under-forming, pillow, global bending after a
3-2-1 release) that is NOT physics-validated. It exercises the pipeline -
design, data, three model classes, conformal calibration, the envelope,
surrogate compensation and its verification - and shows what the numbers
look like; it says nothing about how well the models predict SparLab or a
real part. The real benchmark uses `precomp dataset generate --simulator
sparlab` once sparlab_form exists.

Protocol (revision 2, fixed before its first run; revision 1 had 8 parts per
family and 8 calibration parts, too few for a 90 % interval once calibration
counts parts - see benchmarks/ml_proxy/README.md):
  * design: all seven part families, `--n-per-family` Sobol design points each,
    one material (AA5754-O: one material per model), the default process
    ranges (tool radius 4-8 mm, step-down 0.2-1 mm, thickness 0.6-1.5 mm,
    friction 0.05-0.2), z-level contour tool paths, a 2 mm grid on the 0.2 m
    blank; three commanded variants per point (uncompensated, perturbed, one
    DA step); the time_frac feature from the tool path;
  * the family "saddle" is held out of training and calibration entirely;
  * of the other parts, grouped by part: 20 % test, then 25 % of the rest for
    conformal calibration, the remainder for training;
  * GBMEnsemble (8 members), MLPEnsemble (5 members) and FieldUNet (small),
    each calibrated on the same calibration parts at 90 %;
  * evaluation on the test parts and on the held-out family: per-part dz
    error on the part nodes; interval coverage per part (the single
    calibration split, and the expected coverage over 200 random
    calibration/test partitions of the calibration + test parts, half each),
    the share of parts covered at the level, and the width against a
    constant-width interval calibrated on the same parts; the share of parts
    flagged out of the envelope;
  * surrogate displacement adjustment (GBM) of every uncompensated test
    target inside the envelope (targets outside are refused, and counted),
    at most 8 predictions, stopping at < 2 % improvement (the default) -
    and, for comparison, also within the interval half-width
    (interval_stop, columns stop_*) - verified by forming the compensated
    shape with the proxy; improvement factor = RMS vertical deviation over
    the part of the target formed as is / of the compensated shape formed;
  * envelope probes (GBM): the share of samples flagged for the test parts,
    the held-out family, and the test targets with a material (DC04), a
    thickness (2.5 mm), a release ("clamped_only") and a tool path (spiral;
    only the targets a spiral can form) the model was not trained on;
  * leave each family out (all seven): a GBMEnsemble (4 members, 150
    iterations, 1000 points per sample) and its envelope trained on the other
    families without a fixed grouped 20 % of the parts, then the share
    flagged and the relative error on the left-out family and on those
    in-distribution parts (errors on at most 1500 evenly spread part nodes
    per sample, the envelope's own node set).

`--holdout-family freeform --skip mlp unet lofo` repeats the GBM part with
the freeform family held out instead (benchmarks/ml_proxy/holdout_freeform).
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from collections import Counter
from pathlib import Path

import _bootstrap  # noqa: F401  (puts python/ on sys.path)

import numpy as np
import pandas as pd

LABEL = "proxy - not physics"


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def held_out_arrays(sur, samples):
    """mu, std, y and part id of every part node of `samples`."""
    mus, sds, ys, gs = [], [], [], []
    for s in samples:
        mu, sd = sur.predict_deviation(s.commanded, s.setup)
        sel = sur.feature_maps(s.commanded, s.setup).part & s.valid
        mus.append(mu[sel])
        sds.append(sd[sel])
        ys.append(s.dz[sel])
        gs.append(np.full(sel.sum(), s.part_id, dtype=object))
    return [np.concatenate(v) for v in (mus, sds, ys, gs)]


def leave_each_family_out(samples, cfg, seed, points_per_sample):
    """The LOFO table (see the module docstring): one GBM and envelope per family."""
    from precomp.ml import GBMEnsemble, grouped_split
    from precomp.ml.dataset import Table, index_frame, sample_table
    from precomp.ml.features import feature_maps, global_features
    from precomp.ml.surrogate import fit_envelope

    byid = {s.sample_id: s for s in samples}
    _, test_in = grouped_split(index_frame(samples), test_fraction=0.2, seed=seed + 7)
    test_in = set(test_in)
    # featurise every sample once: training rows (seeded per sample) and the
    # evaluation / envelope nodes
    t0 = time.perf_counter()
    rows, probe = {}, {}
    for i, s in enumerate(samples):
        rows[s.sample_id] = sample_table(s, cfg, points_per_sample=points_per_sample,
                                         rng=np.random.default_rng([seed, i]))
        fm = feature_maps(s.commanded, s.setup, None, cfg)
        nodes = np.flatnonzero(fm.part.ravel())
        nodes = nodes[np.linspace(0, nodes.size - 1, min(1500, nodes.size)).round()
                      .astype(int)]
        d, _ = global_features(s.commanded, s.setup, cfg)
        probe[s.sample_id] = (fm.gather(nodes), s.dz.ravel()[nodes], d)
    log(f"LOFO: featurised {len(samples)} samples in {time.perf_counter() - t0:.0f} s")
    out = []
    for fam in sorted({s.family for s in samples}):
        t0 = time.perf_counter()
        train = [i for i, s in byid.items() if s.family != fam and i not in test_in]
        inside = [i for i in sorted(test_in) if byid[i].family != fam]
        left = [i for i, s in byid.items() if s.family == fam]
        table = Table.concatenate([rows[i] for i in train])
        model = GBMEnsemble(4, max_iter=150, seed=seed).fit(
            table.X, table.y, table.groups, feature_names=table.feature_names)
        rng = np.random.default_rng(seed + 101)          # as train_surrogate thins it
        keep = np.sort(np.concatenate([rng.permutation(np.flatnonzero(table.groups == g))[:400]
                                       for g in np.unique(table.groups)]))
        env = fit_envelope([byid[i] for i in train], cfg, seed=seed,
                           points=(table.X[keep], table.part_id[keep]))
        for split, ids in (("left-out family", left), ("in-distribution test", inside)):
            flags, rel = [], []
            for i in ids:
                X, y, d = probe[i]
                mu, _ = model.predict(X)
                rel.append(np.sqrt(np.mean((mu - y) ** 2)) / np.sqrt(np.mean(y * y)))
                flags.append(not env.assess_features(X, d)["in_envelope"])
            out.append({"left_out_family": fam, "split": split, "n_samples": len(ids),
                        "n_parts": len({byid[i].part_id for i in ids}),
                        "flagged_share": float(np.mean(flags)),
                        "rel_rms_mean": float(np.mean(rel)),
                        "train_parts": len({byid[i].part_id for i in train}),
                        "model": "GBMEnsemble(4 x 150)", "data_source": LABEL})
        log(f"LOFO {fam}: flagged {out[-2]['flagged_share']:.2f} (in-dist "
            f"{out[-1]['flagged_share']:.2f}), rel. error {out[-2]['rel_rms_mean']:.3f} "
            f"(in-dist {out[-1]['rel_rms_mean']:.3f}) [{time.perf_counter() - t0:.0f} s]")
    return pd.DataFrame(out)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--work", help="scratch directory for the data set and bundles")
    ap.add_argument("--n-per-family", type=int, default=10)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--holdout-family", default="saddle")
    ap.add_argument("--skip", nargs="*", default=[], choices=["mlp", "unet", "lofo"])
    ap.add_argument("--quick", action="store_true",
                    help="tiny models, for a smoke run of the script (not for the tables)")
    args = ap.parse_args(argv)
    os.environ["PRECOMP_ML_THREADS"] = str(args.threads)
    try:
        import torch
        torch.set_num_threads(args.threads)
    except ImportError:
        args.skip = sorted(set(args.skip) | {"mlp", "unet"})

    from precomp._util import PrecompError
    from precomp.materials import get_material
    from precomp.ml import (Dataset, DesignSpace, FeatureConfig, FieldUNet, GBMEnsemble,
                            MLPEnsemble, ProxySimulator, design_points, evaluate_surrogate,
                            family_split, generate, grouped_split, save_model,
                            surrogate_compensate, train_surrogate, verify_with_simulator)
    from precomp.ml.uncertainty import partition_coverage

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    work = Path(args.work) if args.work else out / "_work"
    work.mkdir(parents=True, exist_ok=True)
    created = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    t_start = time.perf_counter()
    timings = {}

    # -- data ----------------------------------------------------------------------
    t0 = time.perf_counter()
    # z-level contour paths: freeform parts with several dents and perturbed
    # shapes have several loops per level, which a spiral cannot follow
    space = DesignSpace(grid_spacing=2e-3, materials=("AA5754-O",),
                        base_setup={"toolpath_style": "contour"})
    points = design_points(space, args.n_per_family, args.seed)
    ds = Dataset.create(work / "data", created_at=created, exist_ok=True,
                        description="proxy demonstration - NOT physics")
    rep = generate(ds, points, ProxySimulator(), created_at=created, seed=args.seed)
    if rep.failed:
        log(f"{len(rep.failed)} generation failures: {rep.summary()['failures']}")
    idx = ds.index()
    timings["generate_s"] = time.perf_counter() - t0
    log(f"data: {len(idx)} samples, {idx['part_id'].nunique()} parts, "
        f"{len(rep.failed)} failures")

    rest, fam_ids = family_split(idx, args.holdout_family)
    rest_idx = idx[idx["sample_id"].isin(rest)]
    trcal, test_ids = grouped_split(rest_idx, test_fraction=0.2, seed=args.seed)
    train_ids, cal_ids = grouped_split(rest_idx[rest_idx["sample_id"].isin(trcal)],
                                       test_fraction=0.25, seed=args.seed + 1)
    split_of = {**{i: "train" for i in train_ids}, **{i: "calibration" for i in cal_ids},
                **{i: "test" for i in test_ids}, **{i: "held_out_family" for i in fam_ids}}
    samples = {s.sample_id: s for s in ds}
    design_rows = []
    for sid, s in samples.items():
        st = s.forming_setup()
        sel = s.commanded.z < -1e-6
        design_rows.append({"sample_id": sid, "part_id": s.part_id, "family": s.family,
                            "kind": s.kind, "split": split_of[sid],
                            "tool_radius_m": st.tool_radius, "step_down_m": st.step_down,
                            "thickness_m": st.thickness, "friction": st.friction,
                            "material": st.material.name,
                            "commanded_depth_m": s.commanded.depth,
                            "dz_rms_part_m": float(np.sqrt(np.mean(s.dz[sel] ** 2))),
                            "data_source": LABEL})
    pd.DataFrame(design_rows).sort_values("sample_id").to_csv(out / "design.csv", index=False,
                                                              float_format="%.6g")
    pick = lambda ids: [samples[i] for i in ids]  # noqa: E731
    train, cal, test, fam = pick(train_ids), pick(cal_ids), pick(test_ids), pick(fam_ids)
    nparts = lambda xs: len({s.part_id for s in xs})  # noqa: E731
    log(f"split: train {nparts(train)} / calibration {nparts(cal)} / test {nparts(test)} "
        f"parts, held-out family {args.holdout_family}: {nparts(fam)}")

    # -- models --------------------------------------------------------------------
    cfg = FeatureConfig()                                   # tool-path pseudo-time
    q = args.quick
    makers = {"GBMEnsemble": lambda: GBMEnsemble(2 if q else 8, max_iter=40 if q else 300,
                                                 seed=args.seed)}
    if "mlp" not in args.skip:
        makers["MLPEnsemble"] = lambda: MLPEnsemble(2 if q else 5, hidden=(128, 128),
                                                    epochs=3 if q else 40, seed=args.seed)
    if "unet" not in args.skip:
        makers["FieldUNet"] = lambda: FieldUNet(resolution=40, channels=(12, 24, 48),
                                                epochs=4 if q else 100, features=cfg,
                                                seed=args.seed)
    summary, per_part, per_region, calib, surrogates = [], [], [], [], {}
    for name, make in makers.items():
        t0 = time.perf_counter()
        sur = train_surrogate(make(), train, calibration=cal, features=cfg,
                              points_per_sample=300 if q else 2000, seed=args.seed, n_jobs=2)
        timings[f"train_{name}_s"] = time.perf_counter() - t0
        log(f"{name}: trained and calibrated on {sur.calibrator.n_parts} parts in "
            f"{timings[f'train_{name}_s']:.0f} s")
        surrogates[name] = sur
        t0 = time.perf_counter()
        ev_test = evaluate_surrogate(sur, test, split="test")
        ev_fam = evaluate_surrogate(sur, fam, split="held_out_family")
        timings[f"evaluate_{name}_s"] = time.perf_counter() - t0
        # expected coverage over random calibration/test partitions of the
        # calibration + test parts (the model never saw them), half each
        try:
            pc = partition_coverage(*held_out_arrays(sur, cal + test), level=0.9,
                                    n_partitions=200, seed=args.seed)
        except PrecompError as exc:                  # too few parts for 90 % (--quick)
            log(f"{name}: no expected coverage: {exc}")
            pc = {k: np.nan for k in ("mean", "std", "share_parts_at_level", "width_ratio")}
        for ev in (ev_test, ev_fam):
            for frame, sink in ((ev.per_part, per_part), (ev.per_region, per_region),
                                (ev.calibration, calib)):
                f = frame.copy()
                f.insert(0, "model", name)
                sink.append(f)
            split = ev.per_part["split"].iloc[0]
            s = ev.summary()[split]
            c90 = ev.calibration.set_index("level").loc[0.9]
            row = {"model": name, "split": split, "n_parts": s["n_parts"],
                   "n_samples": s["n_samples"], "dz_rms_mean_mm": 1e3 * s["dz_rms_mean_m"],
                   "err_rms_mean_mm": 1e3 * s["err_rms_mean_m"],
                   "err_rms_median_mm": 1e3 * s["err_rms_median_m"],
                   "err_max_abs_mean_mm": 1e3 * s["err_max_abs_mean_m"],
                   "rel_rms_mean": s["rel_rms_mean"],
                   "calibration_parts": sur.calibrator.n_parts,
                   "coverage90_part_mean": float(c90["coverage_part_mean"]),
                   "share_parts_covered90": float(c90["share_parts_at_level"]),
                   "part_coverage90_min": float(c90["part_coverage_min"]),
                   "width90_mm": 1e3 * float(c90["mean_width_m"]),
                   "baseline_width90_mm": 1e3 * float(c90["baseline_width_m"]),
                   "width_ratio90": float(c90["width_ratio"]),
                   "expected_coverage90": np.nan, "expected_coverage90_sd": np.nan,
                   "expected_share_parts_covered90": np.nan,
                   "expected_width_ratio90": np.nan,
                   "ood_flagged_share": s["ood_flagged_share"],
                   "train_parts": sur.training["n_parts"], "data_source": LABEL}
            if ev is ev_test:
                row.update(expected_coverage90=pc["mean"], expected_coverage90_sd=pc["std"],
                           expected_share_parts_covered90=pc["share_parts_at_level"],
                           expected_width_ratio90=pc["width_ratio"])
            summary.append(row)
        log(f"{name}: test rel. error {summary[-2]['rel_rms_mean']:.3f}, coverage90 per part "
            f"{summary[-2]['coverage90_part_mean']:.3f} (expected {pc['mean']:.3f} "
            f"+- {pc['std']:.3f}, width ratio {pc['width_ratio']:.2f}), flagged "
            f"{summary[-2]['ood_flagged_share']:.2f} / {summary[-1]['ood_flagged_share']:.2f} "
            "(held-out family)")
        save_model(sur, work / "models" / name,
                   {"created_at": created, "notes": "proxy demonstration - NOT physics",
                    "metrics": {"test": {**ev_test.summary()["test"], "data_source": LABEL}}},
                   overwrite=True)
    pd.DataFrame(summary).to_csv(out / "summary.csv", index=False, float_format="%.4g")
    pd.concat(per_part).to_csv(out / "per_part.csv", index=False, float_format="%.6g")
    pd.concat(per_region).to_csv(out / "per_region.csv", index=False, float_format="%.6g")
    pd.concat(calib).to_csv(out / "calibration.csv", index=False, float_format="%.6g")

    # -- compensation ----------------------------------------------------------------
    t0 = time.perf_counter()
    sur = surrogates["GBMEnsemble"]
    proxy = ProxySimulator()
    comp_rows = []
    for s in [s for s in test if s.kind == "uncompensated"]:
        row = {"sample_id": s.sample_id, "family": s.family, "model": "GBMEnsemble",
               "verified_by": LABEL, "data_source": LABEL}
        for tag, stop in (("", False), ("stop_", True)):
            try:
                res = surrogate_compensate(s.target, s.setup, sur, iterations=2 if q else 8,
                                           interval_stop=stop)
            except PrecompError as exc:
                row.update({f"{tag}stopped": "refused: target outside the envelope",
                            f"{tag}iterations": 0})
                row["target_reasons"] = str(exc).split("reasons: ")[-1].split(";")[0]
                continue
            v = verify_with_simulator(res, s.setup, proxy)
            vd = v["compensated"]["vertical_deviation"]["part"]
            row.update({
                f"{tag}iterations": len(res.history), f"{tag}stopped": res.stopped,
                f"{tag}predicted_rms_m": v["predicted_part_rms_m"],
                f"{tag}achieved_rms_m": v["achieved_part_rms_m"],
                f"{tag}achieved_max_abs_m": vd["max_abs"],
                f"{tag}improvement_factor": v["improvement_factor"],
                f"{tag}in_envelope": res.in_envelope})
            if not tag:
                pm = res.predicted_metrics()
                row.update(interval_halfwidth_mean_m=pm.get("interval_halfwidth_mean_m"),
                           uncompensated_rms_m=v["uncompensated"]["vertical_deviation"][
                               "part"]["rms"],
                           uncompensated_max_abs_m=v["uncompensated"]["vertical_deviation"][
                               "part"]["max_abs"],
                           achieved_normal_rms_m=v["compensated"]["normal_deviation"]["part"][
                               "rms"])
        comp_rows.append(row)
    comp = pd.DataFrame(comp_rows)
    comp.to_csv(out / "compensation.csv", index=False, float_format="%.6g")
    timings["compensation_s"] = time.perf_counter() - t0
    done = comp[comp["iterations"] > 0]
    log(f"compensation: {len(done)} of {len(comp)} targets compensated, improvement factor "
        f"median {done['improvement_factor'].median():.1f} (min "
        f"{done['improvement_factor'].min():.1f}); with the interval stop "
        f"{done['stop_improvement_factor'].median():.1f}")

    # -- envelope probes (GBM surrogate) -------------------------------------------------
    from precomp.fea.deck import make_toolpath

    def spiral(s):
        """The target's setup with a spiral path, or None when a spiral cannot
        form it (several pockets on one level): such a target is not probed."""
        st = s.forming_setup().replace(toolpath_style="spiral")
        try:
            make_toolpath(st, s.commanded)
        except (ValueError, PrecompError):
            return None
        return st

    base = [s for s in test if s.kind == "uncompensated"]
    cases = {
        "test parts (in distribution)": [(s.commanded, s.setup) for s in test],
        f"held-out family {args.holdout_family}": [(s.commanded, s.setup) for s in fam],
        "test targets, material DC04 (not trained)":
            [(s.commanded, s.forming_setup().replace(material=get_material("DC04")))
             for s in base],
        "test targets, thickness 2.5 mm (trained 0.6-1.5 mm)":
            [(s.commanded, s.forming_setup().replace(thickness=2.5e-3)) for s in base],
        "test targets, release clamped_only (trained 3-2-1)":
            [(s.commanded, s.forming_setup().replace(release="clamped_only")) for s in base],
        "test targets, spiral tool path (trained contour)":
            [(s.commanded, st) for s, st in ((s, spiral(s)) for s in base) if st is not None]}
    probe_rows = []
    for case, items in cases.items():
        reps = [sur.assess(c, st) for c, st in items]
        reasons = Counter(r for a in reps for r in a["reasons"][:3])
        probe_rows.append({"case": case, "n_samples": len(reps),
                           "flagged_share": float(np.mean([not a["in_envelope"] for a in reps])),
                           "part_score_median": float(np.median([a["part_score"] for a in reps])),
                           "fraction_points_out_median":
                               float(np.median([a["fraction_points_out"] for a in reps])),
                           "top_reasons": ";".join(k for k, _ in reasons.most_common(3)),
                           "model": "GBMEnsemble", "data_source": LABEL})
    probes = pd.DataFrame(probe_rows)
    probes.to_csv(out / "ood_probes.csv", index=False, float_format="%.4g")
    log("envelope: " + "; ".join(f"{r['case']}: {r['flagged_share']:.2f} flagged"
                                 for r in probe_rows))

    # -- leave each family out ------------------------------------------------------------
    if "lofo" not in args.skip:
        t0 = time.perf_counter()
        lofo = leave_each_family_out(list(samples.values()), cfg, args.seed,
                                     300 if q else 1000)
        lofo.to_csv(out / "leave_family_out.csv", index=False, float_format="%.4g")
        timings["leave_family_out_s"] = time.perf_counter() - t0

    timings["total_s"] = time.perf_counter() - t_start
    cmd = [Path(sys.argv[0]).name] + [("<scratch>" if a == args.work else a)
                                      for a in sys.argv[1:]]
    meta = {"created_at": created, "data_source": LABEL, "seed": args.seed,
            "protocol_revision": 2,
            "n_per_family": args.n_per_family, "holdout_family": args.holdout_family,
            "threads": args.threads, "timings_s": timings,
            "n_samples": int(len(idx)), "n_parts": int(idx["part_id"].nunique()),
            "split_parts": {"train": nparts(train), "calibration": nparts(cal),
                            "test": nparts(test), "held_out_family": nparts(fam)},
            "generation_failures": rep.summary()["failures"], "command": " ".join(cmd)}
    (out / "run.json").write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")
    log(f"done in {timings['total_s']:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
