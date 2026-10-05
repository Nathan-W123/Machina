#!/usr/bin/env python3
"""ML first shot + an ML-based second round, on the 8 test parts.

Parts above TOLERANCE after the ML first shot (combined v1: DSIF + rim pass,
transfer MLP v1, DA + optimiser, std weight 0.25; already simulated in
work/combined_rim) get a second round that keeps the ML in the loop:

1. the simulated first shot tells where the model was wrong for THIS part:
   bias = formed_sim - (command + predicted dz at the command), smoothed over
   the part at the tool scale (BIAS_SMOOTHING) so ring noise the tool cannot
   form is not fed back;
2. the surrogate is corrected for this part (predicted dz + bias), which is
   exact at the first-shot command by construction;
3. the same compensation (DA + optimiser, same masks and limits) is re-run on
   the corrected surrogate, and the new command is simulated once.

The round is kept only if it is better than the first shot. SparLab
simulation, penalty 10. Resumable (per-part JSON + the run cache).

    python3 python/scripts/combined_ml2.py all | report
"""
from __future__ import annotations

import argparse
import json
import sys
import threading

import _bootstrap  # noqa: F401

import numpy as np
import pandas as pd

import combined_rim2 as c2
import combined_rim_ml as cr

TOLERANCE_MM = 0.15
BIAS_SMOOTHING = 1.0e-3   # [m], the optimiser's smoothing scale
CWM = c2.WORK / "combined_ml2"
OUTM = c2.OUT_ML / "combined_ml2"
CW1 = c2.WORK / "combined_rim"
FB = c2.WORK / "combined_fallback"       # the stopped engine-round run, for the table


class BiasCorrected:
    """A surrogate whose predicted dz is shifted by a fixed per-part bias map."""

    def __init__(self, base, bias: np.ndarray):
        self.base, self.bias = base, bias

    def predict_deviation(self, commanded, setup, target=None):
        mu, sd = self.base.predict_deviation(commanded, setup, target)
        return mu + self.bias, sd

    def __getattr__(self, name):
        return getattr(self.base, name)


def corrected_model(ctx, base, pid: str):
    from precomp.compensation import smooth_update
    from precomp.geometry.heightmap import HeightMap

    target = ctx.points[pid].target()
    cmd = HeightMap.load(CW1 / f"{pid}-combined.npz")
    formed = HeightMap.load(CW1 / f"{pid}-verify-formed.npz")
    mu, _ = base.predict_deviation(cmd, ctx.setup, target)
    raw = formed.z - (cmd.z + mu)
    part = target.mask & (target.z < -1e-6)
    bias = smooth_update(np.where(part, raw, 0.0), part, BIAS_SMOOTHING, target.grid.h)
    return BiasCorrected(base, bias), float(np.sqrt(np.mean(raw[part] ** 2)))


def one(ctx, base, pid: str, max_time: float) -> None:
    f = CWM / f"{pid}.json"
    if f.exists():
        return
    rec1 = json.loads((CW1 / f"{pid}-verify.json").read_text())
    out = {"point_id": pid, "tolerance_mm": TOLERANCE_MM, "data_source": "SparLab simulation",
           "first_shot": {k: rec1.get(k) for k in c2.COLS}}
    if rec1["vert_rms_mm"] <= TOLERANCE_MM:
        out.update(decision="within tolerance after the ML first shot",
                   final=out["first_shot"], simulations=1)
        c2.log(f"ml2 {pid}: first shot {rec1['vert_rms_mm']:.3f} mm - within tolerance")
    else:
        sur, bias_rms = corrected_model(ctx, base, pid)
        c2.log(f"ml2 {pid}: first shot {rec1['vert_rms_mm']:.3f} mm; model error at the shot "
               f"{bias_rms * 1e3:.3f} mm RMS; re-optimising on the corrected model")
        info, cmd2 = c2.command(ctx, sur, pid, CWM, max_time,
                                "ML second round: transfer v1 + per-part bias from the first shot")
        rec2 = c2.sim(ctx, pid, cmd2, CWM / f"{pid}-ml2-verify.json", as_sample=False)
        r2 = {k: rec2.get(k) for k in c2.COLS} if rec2.get("ok") else None
        better = r2 is not None and r2["vert_rms_mm"] < rec1["vert_rms_mm"]
        out.update(model_error_at_shot_mm=bias_rms * 1e3,
                   predicted=info["optimizer"].get("vert_rms_mm"),
                   ml2_round=r2 if r2 else {"error": rec2.get("termination")},
                   decision="ML round " + ("kept (better)" if better else
                                           "rejected (not better: keep the first shot)"),
                   final=r2 if better else out["first_shot"], simulations=2)
        c2.log(f"ml2 {pid}: LANDED - first shot {rec1['vert_rms_mm']:.3f} -> ML round "
               f"{(r2 or {}).get('vert_rms_mm', float('nan')):.3f} mm "
               f"(predicted {out['predicted']:.3f}; {'kept' if better else 'rejected'})")
    CWM.mkdir(parents=True, exist_ok=True)
    f.write_text(json.dumps(out, indent=1, default=str) + "\n")


def cmd_all(args) -> None:
    from precomp.ml import load_model

    c2.STD_WEIGHT = c2.STD_WEIGHT_V1           # the first shot's optimiser setting (0.25)
    ctx = cr.Ctx()
    base = load_model(CW1 / "model" / "transfer_mlp")
    todo = list(c2.TEST)
    lock = threading.Lock()
    c2.log(f"combined_ml2: tolerance {TOLERANCE_MM} mm, bias smoothing "
           f"{BIAS_SMOOTHING * 1e3:.1f} mm, std weight {c2.STD_WEIGHT}, "
           f"{c2.SLOTS} runs x {c2.SIM_THREADS} threads")

    def worker() -> None:
        while True:
            with lock:
                if not todo:
                    return
                pid = todo.pop(0)
            try:
                one(ctx, base, pid, args.max_time)
            except Exception as exc:  # recorded, the others go on
                c2.log(f"ml2 {pid}: ERROR {type(exc).__name__}: {exc}")

    ths = [threading.Thread(target=worker, daemon=True) for _ in range(c2.SLOTS)]
    for t in ths:
        t.start()
    for t in ths:
        t.join()
    cmd_report(args)
    (CWM / "done").write_text(json.dumps({"finished": c2.now_iso()}) + "\n")
    c2.log("combined_ml2: done")


def cmd_report(args) -> None:
    rows = []
    for p in c2.TEST:
        d = json.loads((CWM / f"{p}.json").read_text()) if (CWM / f"{p}.json").exists() else {}
        fb = FB / f"{p}.json"
        fbd = json.loads(fb.read_text()) if fb.exists() else {}
        rows.append({"point_id": p,
                     "first_shot_mm": (d.get("first_shot") or {}).get("vert_rms_mm", np.nan),
                     "engine_round_mm": (fbd.get("engine_round") or {}).get("vert_rms_mm", np.nan),
                     "ml_round_mm": (d.get("ml2_round") or {}).get("vert_rms_mm", np.nan),
                     "ml_round_predicted_mm": d.get("predicted", np.nan),
                     "final_mm": (d.get("final") or {}).get("vert_rms_mm", np.nan),
                     "simulations": d.get("simulations", np.nan),
                     "decision": d.get("decision", "not run"), "data_source": "SparLab simulation"})
    df = pd.DataFrame(rows)
    for c in df.columns:
        if c.endswith("_mm"):
            df[c] = pd.to_numeric(df[c], errors="coerce")
    OUTM.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTM / "results.csv", index=False)
    f3 = lambda x: "-" if pd.isna(x) else f"{x:.3f}"  # noqa: E731
    L = ["# ML first shot + ML second round (8 test parts)", "",
         f"**SparLab simulation**, DSIF + rim pass, penalty 10. Tolerance {TOLERANCE_MM} mm "
         "(vertical RMS). A part above it after the ML first shot gets a second round: the "
         "surrogate is corrected by the first shot's simulated model error (smoothed at "
         f"{BIAS_SMOOTHING * 1e3:.0f} mm) and the compensation is re-run on it; kept only if "
         "better. 'Engine round' = the stopped full-step displacement-adjustment round "
         "(combined_fallback), for comparison. 8 parts: a small sample.", "",
         "| part | ML first shot | engine round (full DA step) | ML round predicted | ML round "
         "| final | simulations | decision |", "|---|--:|--:|--:|--:|--:|--:|---|"]
    for r in rows:
        L.append(f"| {r['point_id']} | {f3(r['first_shot_mm'])} | {f3(r['engine_round_mm'])} | "
                 f"{f3(r['ml_round_predicted_mm'])} | {f3(r['ml_round_mm'])} | "
                 f"{f3(r['final_mm'])} | {r['simulations']} | {r['decision']} |")
    m = df.mean(numeric_only=True)
    L.append(f"| **mean** | **{f3(m['first_shot_mm'])}** | | | | **{f3(m['final_mm'])}** | "
             f"**{m['simulations']:.2f}** | |")
    L += ["", "Made by `python/scripts/combined_ml2.py`."]
    (OUTM / "README.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("what", choices=["all", "report"])
    ap.add_argument("--max-time", type=float, default=240.0)
    a = ap.parse_args()
    {"all": cmd_all, "report": cmd_report}[a.what](a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
