#!/usr/bin/env python3
"""Combined v3: the v2 transfer correction (fix 1) with the v1 optimiser std
weight 0.25 (fix 2 undone), on the 8 test parts - the ablation that splits
combined_rim2's v1 -> v2 change between the two fixes.

SparLab simulation, DSIF + rim pass, penalty 10; one verifying run per test
part. Resumable (per-part JSON + the run cache). Reuses combined_rim2's
`command` and `sim` unchanged, with its module std weight set to 0.25.

    python3 python/scripts/combined_rim3.py all      # commands + runs + report
    python3 python/scripts/combined_rim3.py report
"""
from __future__ import annotations

import argparse
import json
import sys
import threading
import time
from pathlib import Path

import _bootstrap  # noqa: F401

import numpy as np
import pandas as pd

import combined_rim2 as c2
import combined_rim_ml as cr

STD_WEIGHT_V3 = 0.25
CW3 = c2.WORK / "combined_rim3"
OUT3 = c2.OUT_ML / "combined_rim3"
MODEL = "v3: transfer v2 (refitted on 8 uncompensated + 8 compensated fit runs), std weight 0.25"


def cmd_all(args) -> None:
    from precomp.ml import load_model

    c2.STD_WEIGHT = STD_WEIGHT_V3      # `command` reads the module value
    ctx = cr.Ctx()
    sur = load_model(c2.CW2 / "model" / "transfer_mlp_v2")
    d = CW3 / "test"
    todo = list(c2.TEST)
    lock = threading.Lock()
    c2.log(f"combined_rim3: std weight {STD_WEIGHT_V3}, model transfer_mlp_v2, "
           f"{c2.SLOTS} runs x {c2.SIM_THREADS} threads")

    def worker() -> None:
        while True:
            with lock:
                if not todo:
                    return
                pid = todo.pop(0)
            try:
                _, comp = c2.command(ctx, sur, pid, d, args.max_time, MODEL)
                c2.sim(ctx, pid, comp, d / f"{pid}-verify.json", as_sample=False)
            except Exception as exc:  # recorded, the others go on
                c2.log(f"v3 {pid}: ERROR {type(exc).__name__}: {exc}")

    ths = [threading.Thread(target=worker, daemon=True) for _ in range(c2.SLOTS)]
    for t in ths:
        t.start()
    for t in ths:
        t.join()
    cmd_report(args)
    n_ok = sum(json.loads(f.read_text()).get("ok", False) for f in d.glob("*-verify.json"))
    (CW3 / "done").write_text(json.dumps({"finished": c2.now_iso(), "test_runs_ok": n_ok}) + "\n")
    c2.log(f"combined_rim3: done ({n_ok} of {len(c2.TEST)} test runs ok)")


def cmd_report(args) -> None:
    prev = pd.read_csv(c2.OUT_ML / "combined_rim2" / "results.csv").set_index("point_id")
    rows = []
    for p in c2.TEST:
        f, fp = CW3 / "test" / f"{p}-verify.json", CW3 / "test" / f"{p}-combined.json"
        v = json.loads(f.read_text()) if f.exists() else {}
        pr = json.loads(fp.read_text())["optimizer"] if fp.exists() else {}
        r = {"point_id": p}
        for k in ("plate_unc", "mlp_plate", "v1", "v2"):
            for q in ("vert_rms_mm", "upper_band_rms_mm", "deep_rms_mm"):
                r[f"{k}_{q}"] = prev.loc[p, f"{k}_{q}"]
        r["v1_pred_vert_rms_mm"] = prev.loc[p, "v1_pred_vert_rms_mm"]
        r["v2_pred_vert_rms_mm"] = prev.loc[p, "v2_pred_vert_rms_mm"]
        r["v3_pred_vert_rms_mm"] = pr.get("vert_rms_mm", np.nan)
        for q in ("vert_rms_mm", "upper_band_rms_mm", "deep_rms_mm"):
            r[f"v3_{q}"] = v.get(q) if v.get("ok") else np.nan
        r["v3_status"] = "ok" if v.get("ok") else ("failed: " + str(v.get("termination"))
                                                   if v else "not run")
        r["data_source"] = "SparLab simulation"
        rows.append(r)
    df = pd.DataFrame(rows)
    for c in df.columns:
        if c.endswith("_mm"):
            df[c] = pd.to_numeric(df[c], errors="coerce")
    OUT3.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT3 / "results.csv", index=False)
    f3 = lambda x: "-" if pd.isna(x) else f"{x:.3f}"  # noqa: E731
    L = ["# Combined v3: retrained transfer correction (fix 1) + the v1 std weight 0.25",
         "", "**SparLab simulation**, penalty 10, DSIF + rim pass. 8 test parts, one verifying "
         "run each: a small sample. Predictions are the surrogate's, not results.", "",
         "The ablation of combined_rim2: v1 = transfer v1 (linear, 8 uncompensated fit runs), "
         "std weight 0.25; v2 = transfer v2 (gbm, 8 uncompensated + 8 compensated fit runs), "
         "std weight 1.0; **v3 = transfer v2, std weight 0.25**. v1 -> v3 is fix 1 alone; "
         "v3 -> v2 is fix 2 alone (on the v2 model).", "",
         "| part | plate uncompensated | ML-MLP plate | v1 | v2 | v3 | v3 predicted | v3 status |",
         "|---|--:|--:|--:|--:|--:|--:|---|"]
    for r in rows:
        L.append(f"| {r['point_id']} | {f3(r['plate_unc_vert_rms_mm'])} | "
                 f"{f3(r['mlp_plate_vert_rms_mm'])} | {f3(r['v1_vert_rms_mm'])} | "
                 f"{f3(r['v2_vert_rms_mm'])} | {f3(r['v3_vert_rms_mm'])} | "
                 f"{f3(r['v3_pred_vert_rms_mm'])} | "
                 f"{r['v3_status']} |")
    m = df.mean(numeric_only=True)
    L.append(f"| **mean** | **{f3(m['plate_unc_vert_rms_mm'])}** | **{f3(m['mlp_plate_vert_rms_mm'])}**"
             f" | **{f3(m['v1_vert_rms_mm'])}** | **{f3(m['v2_vert_rms_mm'])}** | "
             f"**{f3(m['v3_vert_rms_mm'])}** | **{f3(m['v3_pred_vert_rms_mm'])}** | |")
    L += ["", f"v3 simulated on {int(df['v3_vert_rms_mm'].notna().sum())} of {len(df)} parts.", "",
          "| mean [mm] | v1 | v2 | v3 |", "|---|--:|--:|--:|"]
    for q, lab in (("upper_band_rms_mm", "upper-band RMS"), ("deep_rms_mm", "deep RMS")):
        L.append(f"| {lab} | {f3(m['v1_' + q])} | {f3(m['v2_' + q])} | {f3(m['v3_' + q])} |")
    L.append(f"| predicted vertical RMS | {f3(m['v1_pred_vert_rms_mm'])} | "
             f"{f3(m['v2_pred_vert_rms_mm'])} | {f3(m['v3_pred_vert_rms_mm'])} |")
    L += ["", "Made by `python/scripts/combined_rim3.py` (reuses combined_rim2's `command` and "
          "`sim`)."]
    (OUT3 / "README.md").write_text("\n".join(L) + "\n")
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
