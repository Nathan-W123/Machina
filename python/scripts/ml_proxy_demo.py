#!/usr/bin/env python3
"""End-to-end demonstration of precomp.ml on the ProxySimulator - NOT physics.

usage:
    python3 python/scripts/ml_proxy_demo.py --out benchmarks/ml_proxy \
        [--work scratch_dir] [--n-per-family 8] [--seed 2026] [--threads 2]

Every number this script writes comes from the ProxySimulator, an analytic
stand-in for springback (rim under-forming, pillow, global bending after a
3-2-1 release) that is NOT physics-validated. It exercises the pipeline -
design, data, three model classes, conformal calibration, the envelope,
surrogate compensation and its verification - and shows what the numbers
look like; it says nothing about how well the models predict SparLab or a
real part. The real benchmark uses `precomp dataset generate --simulator
sparlab` once sparlab_form exists.

Protocol (fixed before the run):
  * design: all seven part families, `--n-per-family` Sobol design points each,
    one material (AA5754-O: one material per model), the default process
    ranges (tool radius 4-8 mm, step-down 0.2-1 mm, thickness 0.6-1.5 mm,
    friction 0.05-0.2), z-level contour tool paths, a 2 mm grid on the 0.2 m
    blank; three commanded variants per point (uncompensated, perturbed, one
    DA step); the time_frac feature from the tool path;
  * the family "saddle" is held out of training and calibration entirely;
  * of the other parts, grouped by part: 20 % test, then 20 % of the rest for
    conformal calibration, the remainder for training;
  * GBMEnsemble (8 members), MLPEnsemble (5 members) and FieldUNet (small),
    each calibrated on the same calibration parts at 90 %;
  * evaluation on the test parts and on the held-out family: per-part dz
    error on the part nodes, interval coverage (the single calibration split,
    and the expected coverage over 200 random calibration/test partitions of
    the calibration + test parts), the share of parts flagged out of the
    envelope;
  * surrogate displacement adjustment (GBM) of every uncompensated test
    target, verified by forming the compensated shape with the proxy;
    improvement factor = RMS vertical deviation over the part of the target
    formed as is / of the compensated shape formed;
  * envelope probes (GBM): the share of samples flagged for the test parts,
    the held-out family, and the test targets with a material (DC04) and a
    thickness (2.5 mm) the model was not trained on.

`--holdout-family freeform --skip mlp unet` repeats the GBM part with the
freeform family held out instead (benchmarks/ml_proxy/holdout_freeform).
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import _bootstrap  # noqa: F401  (puts python/ on sys.path)

import numpy as np
import pandas as pd

LABEL = "proxy - not physics"


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--work", help="scratch directory for the data set and bundles")
    ap.add_argument("--n-per-family", type=int, default=8)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--holdout-family", default="saddle")
    ap.add_argument("--skip", nargs="*", default=[], choices=["mlp", "unet"])
    ap.add_argument("--quick", action="store_true",
                    help="tiny models, for a smoke run of the script (not for the tables)")
    args = ap.parse_args(argv)
    os.environ["PRECOMP_ML_THREADS"] = str(args.threads)
    try:
        import torch
        torch.set_num_threads(args.threads)
    except ImportError:
        args.skip = sorted(set(args.skip) | {"mlp", "unet"})

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
                                       test_fraction=0.2, seed=args.seed + 1)
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
    log(f"split: train {len(train)} / calibration {len(cal)} / test {len(test)} samples, "
        f"held-out family {args.holdout_family}: {len(fam)}")

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
        log(f"{name}: trained and calibrated in {timings[f'train_{name}_s']:.0f} s")
        surrogates[name] = sur
        t0 = time.perf_counter()
        ev_test = evaluate_surrogate(sur, test, split="test")
        ev_fam = evaluate_surrogate(sur, fam, split="held_out_family")
        timings[f"evaluate_{name}_s"] = time.perf_counter() - t0
        # expected coverage over random calibration/test partitions of the
        # calibration + test parts (the model never saw them)
        mus, sds, ys, gs = [], [], [], []
        for s in cal + test:
            mu, sd = sur.predict_deviation(s.commanded, s.setup)
            sel = sur.feature_maps(s.commanded, s.setup).part & s.valid
            mus.append(mu[sel])
            sds.append(sd[sel])
            ys.append(s.dz[sel])
            gs.append(np.full(sel.sum(), s.part_id, dtype=object))
        pc = partition_coverage(*(np.concatenate(v) for v in (mus, sds, ys, gs)), level=0.9,
                                n_partitions=200, seed=args.seed)
        for ev in (ev_test, ev_fam):
            for frame, sink in ((ev.per_part, per_part), (ev.per_region, per_region),
                                (ev.calibration, calib)):
                f = frame.copy()
                f.insert(0, "model", name)
                sink.append(f)
            s = ev.summary()[ev.per_part["split"].iloc[0]]
            summary.append({"model": name, "split": ev.per_part["split"].iloc[0],
                            "n_parts": s["n_parts"], "n_samples": s["n_samples"],
                            "dz_rms_mean_mm": 1e3 * s["dz_rms_mean_m"],
                            "err_rms_mean_mm": 1e3 * s["err_rms_mean_m"],
                            "err_rms_median_mm": 1e3 * s["err_rms_median_m"],
                            "err_max_abs_mean_mm": 1e3 * s["err_max_abs_mean_m"],
                            "rel_rms_mean": s["rel_rms_mean"],
                            "coverage90": float(ev.calibration.set_index("level")
                                                .loc[0.9, "coverage"]),
                            "part_coverage90_min": float(ev.calibration.set_index("level")
                                                         .loc[0.9, "part_coverage_min"]),
                            "expected_coverage90": pc["mean"] if ev is ev_test else np.nan,
                            "expected_coverage90_sd": pc["std"] if ev is ev_test else np.nan,
                            "ood_flagged_share": s["ood_flagged_share"],
                            "train_parts": sur.training["n_parts"],
                            "data_source": LABEL})
        log(f"{name}: test rel. error {summary[-2]['rel_rms_mean']:.3f}, "
            f"coverage90 {summary[-2]['coverage90']:.3f} (expected {pc['mean']:.3f} "
            f"+- {pc['std']:.3f}), flagged {summary[-2]['ood_flagged_share']:.2f} / "
            f"{summary[-1]['ood_flagged_share']:.2f} (held-out family)")
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
        res = surrogate_compensate(s.target, s.setup, sur, iterations=2 if q else 8)
        v = verify_with_simulator(res, s.setup, proxy)
        comp_rows.append({
            "sample_id": s.sample_id, "family": s.family, "iterations": len(res.history),
            "stopped": res.stopped,
            "predicted_rms_m": v["predicted_part_rms_m"],
            "achieved_rms_m": v["achieved_part_rms_m"],
            "uncompensated_rms_m": v["uncompensated"]["vertical_deviation"]["part"]["rms"],
            "achieved_max_abs_m": v["compensated"]["vertical_deviation"]["part"]["max_abs"],
            "uncompensated_max_abs_m":
                v["uncompensated"]["vertical_deviation"]["part"]["max_abs"],
            "achieved_normal_rms_m": v["compensated"]["normal_deviation"]["part"]["rms"],
            "improvement_factor": v["improvement_factor"],
            "in_envelope": res.ood["in_envelope"], "model": "GBMEnsemble",
            "verified_by": v["data_source"], "data_source": LABEL})
    comp = pd.DataFrame(comp_rows)
    comp.to_csv(out / "compensation.csv", index=False, float_format="%.6g")
    timings["compensation_s"] = time.perf_counter() - t0
    log(f"compensation: {len(comp)} parts, improvement factor median "
        f"{comp['improvement_factor'].median():.1f} (min {comp['improvement_factor'].min():.1f})")

    # -- envelope probes (GBM surrogate) -------------------------------------------------
    from collections import Counter

    from precomp.materials import get_material

    base = [s for s in test if s.kind == "uncompensated"]
    cases = {
        "test parts (in distribution)": [(s.commanded, s.setup) for s in test],
        f"held-out family {args.holdout_family}": [(s.commanded, s.setup) for s in fam],
        "test targets, material DC04 (not trained)":
            [(s.commanded, s.forming_setup().replace(material=get_material("DC04")))
             for s in base],
        "test targets, thickness 2.5 mm (trained 0.6-1.5 mm)":
            [(s.commanded, s.forming_setup().replace(thickness=2.5e-3)) for s in base]}
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
    timings["total_s"] = time.perf_counter() - t_start
    cmd = [Path(sys.argv[0]).name] + [("<scratch>" if a == args.work else a)
                                      for a in sys.argv[1:]]
    meta = {"created_at": created, "data_source": LABEL, "seed": args.seed,
            "n_per_family": args.n_per_family, "holdout_family": args.holdout_family,
            "threads": args.threads, "timings_s": timings,
            "n_samples": int(len(idx)), "n_parts": int(idx["part_id"].nunique()),
            "generation_failures": rep.summary()["failures"], "command": " ".join(cmd)}
    (out / "run.json").write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")
    log(f"done in {timings['total_s']:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
