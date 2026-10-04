#!/usr/bin/env python3
"""One-part smoke test (SparLab simulation, 1 part: not statistically
meaningful) of two ideas against the upper-band error the ML first shot
leaves on the fine mesh:

  (A) optimiser: the stage_n26 MLP's surrogate DA of the target (as
      `springback_fine_ml.ml_shot`), then `surrogate_optimize` from it; the
      optimiser's command verified by ONE run on exactly ml_shot's path
      (backing plate, `springback_fine` preset, penalty 10);
  (B) edge fix: the UNCOMPENSATED target with the DSIF support and its rim
      pass (support_cone's dsif_rim settings, otherwise
      benchmarks/springback_fine's dsif strategy at penalty 10).

usage (from the worktree root, OPENBLAS_NUM_THREADS=1 PYTHONPATH=python):
    python3 python/scripts/smoke_optimizer_edge.py optimize   # (A) prediction only
    python3 python/scripts/smoke_optimizer_edge.py verify     # (A) the SparLab run
    python3 python/scripts/smoke_optimizer_edge.py edge       # (B) the SparLab run
    python3 python/scripts/smoke_optimizer_edge.py report     # results.csv, README.md

Every run goes through the content-addressed cache in
benchmarks/springback_fine_ml/work/runs (resumable).
"""

from __future__ import annotations

import argparse
import copy
import json
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import _bootstrap  # noqa: F401

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import springback_benchmark as sb  # noqa: E402
import springback_fine_ml as fm  # noqa: E402
import support_validation as sv  # noqa: E402

MM = 1e-3
LABEL = "SparLab simulation"
REPO = HERE.parents[1]
OUT_ML = REPO / "benchmarks" / "springback_fine_ml"
WORK = OUT_ML / "work"
OUT = OUT_ML / "smoke_n26"
PART = "pyramid-s2026-0001"
STAGE, N = "stage_n26", 26
OPT_DIR = WORK / "smoke_n26"


def log(msg: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def cfgs():
    cfg = json.loads((OUT_ML / "config.json").read_text())
    old = json.loads((fm.OLD / "config.json").read_text())
    return cfg, old


def exe() -> str:
    return str((WORK / "bin" / "sparlab_form").resolve())


def the_point(cfg, old):
    pts = fm.design(cfg, old, N, exe())
    return next(p for p in pts if p.point_id == PART), pts


def load_mlp(cfg, old, pts):
    """The stage_n26 MLP exactly as the driver loads it (sb.train_or_load
    with the stage's train / calibration samples: a load, never a retrain)."""
    from precomp.ml import Dataset

    ds = Dataset.create(WORK / "data", created_at="", exist_ok=True)
    spl = fm.cfg_split(cfg)
    pids = {p.point_id for p in pts}
    by: Dict[str, List[Any]] = {"train": [], "calibration": []}
    for s in ds:
        if s.part_id in pids:
            by.setdefault(sb.split_of(s.part_id, spl), []).append(s)
    bundle = WORK / "models" / STAGE / "mlp"
    sur, dt = sb.train_or_load("mlp", by["train"], by["calibration"], cfg, bundle, 1)
    if dt is not None:
        raise SystemExit("the MLP was retrained: the bundle did not match the stage")
    return sur


# ---------------------------------------------------------------------------
def cmd_optimize(args) -> None:
    from precomp.geometry.heightmap import HeightMap
    from precomp.ml import surrogate_compensate
    from precomp.ml.compensate import surrogate_optimize

    cfg, old = cfgs()
    p, pts = the_point(cfg, old)
    sur = load_mlp(cfg, old, pts)
    target = p.target()
    o = fm.da_options(cfg)(p)
    sc = cfg["surrogate"]
    kw = dict(max_wall_angle_deg=o["max_wall_angle_deg"], upper_bound=o["upper_bound"],
              hold_mask=o["hold_mask"], adjust_mask=o["adjust_mask"],
              tool_radius=o.get("tool_radius"))
    t0 = time.perf_counter()
    da = surrogate_compensate(target, p.setup, sur, iterations=sc["iterations"],
                              stagnation=sc["stagnation"], alpha=o["alpha"],
                              direction=o["direction"], smoothing=o["smoothing"], **kw)
    t_da = time.perf_counter() - t0
    # the DA command must be ml_shot's (the cached ML-MLP shot)
    shot = HeightMap.load(WORK / "ml" / STAGE / "mlp" / f"{PART}.npz")
    da_diff = float(np.abs(shot.z - da.compensated.z).max())
    log(f"DA: predicted RMS {da.predicted_metrics()['rms'] / MM:.4f} mm "
        f"({da.stopped}, best {da.best_iteration}); max |DA - ml_shot command| "
        f"{da_diff / MM:.2e} mm")
    t0 = time.perf_counter()
    opt = surrogate_optimize(target, p.setup, sur, start=da.compensated,
                             max_time_s=args.max_time, **kw)
    t_opt = time.perf_counter() - t0
    OPT_DIR.mkdir(parents=True, exist_ok=True)
    opt.compensated.save(OPT_DIR / f"{PART}-optimizer.npz")
    opt.save(OPT_DIR / f"{PART}-optimizer")

    def pred(res):
        d = sv.deviation(res.predicted, target, cfg["upper_band_depth"])
        return {k: d[k] for k in ("vert_rms_mm", "upper_band_rms_mm", "upper_band_bias_mm",
                                  "deep_rms_mm", "deep_bias_mm")}

    info = {"point_id": PART, "data_source": "MLP surrogate prediction (trained on "
            "SparLab simulation, penalty 3)",
            "da": {**pred(da), "stopped": da.stopped, "time_s": t_da,
                   "max_abs_diff_to_ml_shot_mm": da_diff / MM},
            "optimizer": {**pred(opt), "stopped": opt.stopped, "time_s": t_opt,
                          "evaluations": len(opt.history),
                          "in_envelope": opt.in_envelope,
                          "std_mean_mm": float(np.mean(opt.std[target.z < -1e-6])) / MM
                          if opt.std is not None else None,
                          "commanded_depth_mm": opt.compensated.depth / MM,
                          "max_abs_change_from_da_mm":
                              float(np.abs(opt.compensated.z - da.compensated.z).max()) / MM},
            "da_std_mean_mm": float(np.mean(da.std[target.z < -1e-6])) / MM
            if da.std is not None else None}
    sb.jdump(OPT_DIR / f"{PART}-optimizer.json", info)
    log(json.dumps(info, indent=1, default=str))


def cmd_verify(args) -> None:
    """ONE run of the optimiser's command on ml_shot's path."""
    from precomp.geometry.heightmap import HeightMap
    from precomp.ml import SparlabSimulator
    from precomp.ml.generate import sim_job

    cfg, old = cfgs()
    p, _ = the_point(cfg, old)
    comp = HeightMap.load(OPT_DIR / f"{PART}-optimizer.npz")
    test_setup = fm.make_setup(cfg, cfg["preset_test"], exe())
    sim = SparlabSimulator(WORK / "runs", max_workers=1, executor="process")
    t0 = time.perf_counter()
    oc = sim.run([sim_job(test_setup, comp, p.target())])[0]
    rec = {"case": "optimizer", "ok": oc.ok, "deck_hash": oc.provenance.get("deck_hash"),
           "runtime_s": oc.provenance.get("runtime_s"),
           "cache_hit": oc.provenance.get("cache_hit"), "wall_s": time.perf_counter() - t0}
    if oc.ok:
        rec.update(sv.deviation(oc.formed, p.target(), cfg["upper_band_depth"]))
    else:
        rec["error"] = (oc.error or "")[:2000]
    sb.jdump(OPT_DIR / f"{PART}-verify.json", rec)
    log(json.dumps(rec, indent=1, default=str))


def edge_setup(rim: bool):
    """benchmarks/springback_fine's dsif strategy (its make_setup), with
    support_cone's dsif_rim settings when `rim`."""
    from precomp.fea import FormingSetup

    fcfg = json.loads((fm.FINE / "config.json").read_text())
    scfg = json.loads((REPO / "benchmarks" / "support_cone" / "config.json").read_text())
    over = copy.deepcopy(fcfg.get("setup_overrides", {}))
    over.update(copy.deepcopy(fcfg["strategies"]["dsif"]))
    if rim:
        over.update(copy.deepcopy(scfg["strategies"]["dsif_rim"]))
    return FormingSetup.preset(fcfg["preset"], executable=exe(), timeout=fcfg["timeout_s"],
                               **over)


def deck_key(setup, commanded, target) -> str:
    import shutil
    import tempfile

    from precomp.fea.deck import build_deck, deck_hash
    from precomp.fea.runner import sparlab_version

    d = Path(tempfile.mkdtemp(dir=OPT_DIR))
    try:
        build_deck(setup, commanded, d, target=target)
        return deck_hash(d, sparlab_version(setup.resolved_executable()))
    finally:
        shutil.rmtree(d, ignore_errors=True)


def cmd_edge(args) -> None:
    """ONE run: the uncompensated target, dsif + rim pass, penalty 10."""
    from precomp.fea.runner import FormingError, simulate

    cfg, old = cfgs()
    p, _ = the_point(cfg, old)
    target = p.target()
    OPT_DIR.mkdir(parents=True, exist_ok=True)
    # consistency: the same setup without the rim pass is springback_fine's dsif run
    fine = pd.read_csv(fm.FINE / "cases.csv")
    ref = fine[(fine.point_id == PART) & (fine.strategy == "dsif")
               & (fine.method == "uncompensated")]["deck_hash"].iloc[0]
    k0 = deck_key(edge_setup(False), target, target)
    log(f"dsif (no rim pass) deck {k0[:12]}; springback_fine's {ref[:12]}: "
        f"{'SAME' if k0 == ref else 'DIFFERENT'}")
    setup = edge_setup(True)
    log(f"dsif + rim pass: support {setup.support} {setup.resolved_support()}")
    t0 = time.perf_counter()
    rec: Dict[str, Any] = {"case": "dsif_rim_uncompensated", "dsif_deck_matches_fine": k0 == ref}
    try:
        res = simulate(setup, target, WORK / "runs", target=target)
        formed = res.formed_surface(-1, grid=target.grid)
        rec.update(ok=True, deck_hash=res.provenance.get("key"),
                   cache_hit=res.provenance.get("cache_hit"))
        rec.update(sv.physics(res))
        rec.update(sv.deviation(formed, target, cfg["upper_band_depth"]))
    except FormingError as exc:
        rec.update(ok=False, error=str(exc)[:2000],
                   record={k: v for k, v in (exc.record or {}).items()
                           if isinstance(v, (str, int, float, bool))})
    rec["wall_s"] = time.perf_counter() - t0
    sb.jdump(OPT_DIR / f"{PART}-edge.json", rec)
    log(json.dumps(rec, indent=1, default=str))


# ---------------------------------------------------------------------------
COLS = ["vert_rms_mm", "upper_band_rms_mm", "upper_band_bias_mm", "deep_rms_mm",
        "deep_bias_mm"]


def cmd_report(args) -> None:
    rows: List[Dict[str, Any]] = []
    fine = pd.read_csv(fm.FINE / "cases.csv")
    fine = fine[fine.point_id == PART]
    for strat, meth, label in (("backing_plate", "uncompensated", "plate uncompensated"),
                               ("backing_plate", "FE-DA-1", "plate FE-DA-1"),
                               ("dsif", "uncompensated", "dsif uncompensated")):
        r = fine[(fine.strategy == strat) & (fine.method == meth)].iloc[0]
        rows.append({"case": label, "kind": "simulated", "fe_runs": int(r.fe_run),
                     **{c: r[c] for c in COLS}, "deck_hash": r.deck_hash,
                     "source": "benchmarks/springback_fine"})
    hd = pd.read_csv(OUT_ML / STAGE / "headline.csv")
    r = hd[(hd.point_id == PART) & (hd.method == "ML-MLP")].iloc[0]
    pinfo = {}
    pf = OPT_DIR / f"{PART}-optimizer.json"
    if pf.exists():
        pinfo = json.loads(pf.read_text())
    if pinfo:
        rows.append({"case": "ML-MLP n26 (surrogate DA)", "kind": "predicted (MLP)",
                     "fe_runs": 0, **{c: pinfo["da"][c] for c in COLS},
                     "source": "this smoke test (prediction)"})
    rows.append({"case": "ML-MLP n26 (surrogate DA)", "kind": "simulated", "fe_runs": 1,
                 **{c: r[c] for c in COLS}, "deck_hash": r.deck_hash,
                 "source": "benchmarks/springback_fine_ml/stage_n26"})
    if pinfo:
        rows.append({"case": "optimizer (from ML-MLP DA)", "kind": "predicted (MLP)",
                     "fe_runs": 0, **{c: pinfo["optimizer"][c] for c in COLS},
                     "source": "this smoke test (prediction)"})
    for name, case, src in (("verify", "optimizer (from ML-MLP DA)", "this smoke test"),
                            ("edge", "dsif + rim pass uncompensated", "this smoke test")):
        f = OPT_DIR / f"{PART}-{name}.json"
        if not f.exists():
            rows.append({"case": case, "kind": "simulated", "fe_runs": 1,
                         "note": "not run yet", "source": src})
            continue
        x = json.loads(f.read_text())
        rows.append({"case": case, "kind": "simulated", "fe_runs": 1,
                     **{c: x.get(c) for c in COLS}, "deck_hash": x.get("deck_hash"),
                     "note": "" if x.get("ok") else "FAILED: " + str(x.get("error"))[:300],
                     "source": src})
    df = pd.DataFrame(rows)
    df.insert(0, "point_id", PART)
    df["data_source"] = LABEL
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "results.csv", index=False, float_format="%.4f")

    def f3(v):
        return "-" if v is None or (isinstance(v, float) and not np.isfinite(v)) \
            else f"{v:+.3f}" if isinstance(v, float) and v < 0 else f"{v:.3f}"

    L = [f"# Smoke test on {PART}: optimiser and DSIF rim pass", "",
         f"**{LABEL}** (sparlab_form; predictions are the stage_n26 MLP surrogate's). "
         "**One part: not statistically meaningful.** The part was chosen in advance "
         "because the ML first shot barely helped it.", "",
         "All simulated rows at penalty 10 (`springback_fine` preset). Deviation formed - "
         "target, vertical, in mm; upper band = the top 1 mm of depth (next to the clamp), "
         "deep = the rest of the part.", "",
         "| case | kind | FE runs | vertical RMS | upper-band RMS | upper-band bias | "
         "deep RMS | deep bias | note |", "|---|---|--:|--:|--:|--:|--:|--:|---|"]
    for _, x in df.iterrows():
        note = x.get("note")
        L.append(f"| {x['case']} | {x['kind']} | {x['fe_runs']} | "
                 + " | ".join(f3(x.get(c)) for c in COLS)
                 + f" | {'' if pd.isna(note) else note} |")
    if pinfo:
        o = pinfo["optimizer"]
        L += ["", f"Optimiser: {o['evaluations']} evaluation batches, {o['time_s']:.0f} s, "
              f"stopped '{o['stopped']}', in envelope {o['in_envelope']}, largest change "
              f"of the command from the DA command {o['max_abs_change_from_da_mm']:.3f} mm."]
    L += ["", "Made by `python/scripts/smoke_optimizer_edge.py` (optimize, verify, edge, "
          "report). The optimiser's command was verified on the ML first shot's path "
          "(backing plate, 1 mm clearance); the rim-pass run is "
          "benchmarks/springback_fine's dsif setup with `rim_pass: true` "
          "(benchmarks/support_cone's dsif_rim).", ""]
    (OUT / "README.md").write_text("\n".join(L))
    print(df.to_string(index=False))


def main() -> int:
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("PRECOMP_ML_THREADS", "1")
    try:
        import torch
        torch.set_num_threads(1)
    except ImportError:
        pass
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("what", choices=["optimize", "verify", "edge", "report"])
    ap.add_argument("--max-time", type=float, default=240.0,
                    help="optimiser wall-time cap [s]")
    a = ap.parse_args()
    {"optimize": cmd_optimize, "verify": cmd_verify, "edge": cmd_edge,
     "report": cmd_report}[a.what](a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
