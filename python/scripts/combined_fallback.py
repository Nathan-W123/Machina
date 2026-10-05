#!/usr/bin/env python3
"""ML first shot + automatic engine fallback, on the 8 test parts.

The process: the ML first shot (combined v1: DSIF + rim pass, transfer MLP,
DA + optimiser; already simulated in work/combined_rim) is checked; if its
vertical RMS error is above TOLERANCE, ONE engine correction follows - a
displacement-adjustment step from the simulated error (deep region only, the
rim band held as the ML shot does, the same wall-angle / upper-bound / tool
reach conditioning) - and that command is simulated. Parts already within
TOLERANCE stop after the first shot.

SparLab simulation, penalty 10, DSIF + rim pass. Resumable (per-part JSON +
the run cache). Logs each part as it lands.

    python3 python/scripts/combined_fallback.py all
    python3 python/scripts/combined_fallback.py report
"""
from __future__ import annotations

import argparse
import json
import sys
import threading
from pathlib import Path

import _bootstrap  # noqa: F401

import numpy as np
import pandas as pd

import combined_rim2 as c2
import combined_rim_ml as cr

TOLERANCE_MM = 0.15
ALPHA = 1.0
CWF = c2.WORK / "combined_fallback"
OUTF = c2.OUT_ML / "combined_fallback"
CW1 = c2.WORK / "combined_rim"          # the ML first shot (combined v1)


def first_shot(pid: str):
    from precomp.geometry.heightmap import HeightMap
    rec = json.loads((CW1 / f"{pid}-verify.json").read_text())
    return rec, HeightMap.load(CW1 / f"{pid}-combined.npz"), \
        HeightMap.load(CW1 / f"{pid}-verify-formed.npz")


def engine_command(ctx, pid: str, cmd, formed):
    """One DA step from the simulated first shot (the deep region adjusted)."""
    from precomp.compensation import _apply_update, _condition, error_field
    from precomp.fea.support import command_upper_bound

    target = ctx.points[pid].target()
    hold, adjust, _, _ = cr.masks(ctx, target)
    da = ctx.cfg["da"]
    err = error_field(formed, target, da["direction"])
    z_new = _apply_update(cmd, target, err, adjust, ALPHA, da["direction"], da["smoothing"])
    return _condition(cmd, z_new, hold, target, da["max_wall_angle_deg"],
                      upper_bound=command_upper_bound(ctx.setup, target),
                      tool_radius=ctx.setup.tool_radius)


def one(ctx, pid: str) -> None:
    f = CWF / f"{pid}.json"
    if f.exists():
        return
    rec1, cmd, formed = first_shot(pid)
    out = {"point_id": pid, "tolerance_mm": TOLERANCE_MM,
           "first_shot": {k: rec1.get(k) for k in c2.COLS}, "data_source": "SparLab simulation"}
    if not rec1.get("ok", True):
        out.update(decision="first shot failed", final=None, simulations=1)
    elif rec1["vert_rms_mm"] <= TOLERANCE_MM:
        out.update(decision="within tolerance after the ML first shot", final=out["first_shot"],
                   simulations=1)
        c2.log(f"fallback {pid}: first shot {rec1['vert_rms_mm']:.3f} mm <= "
               f"{TOLERANCE_MM} mm - no engine round")
    else:
        c2.log(f"fallback {pid}: first shot {rec1['vert_rms_mm']:.3f} mm > {TOLERANCE_MM} mm "
               "- one engine correction round")
        CWF.mkdir(parents=True, exist_ok=True)
        cmd2 = engine_command(ctx, pid, cmd, formed)
        cmd2.save(CWF / f"{pid}-engine.npz")
        rec2 = c2.sim(ctx, pid, cmd2, CWF / f"{pid}-engine-verify.json", as_sample=False)
        r2 = {k: rec2.get(k) for k in c2.COLS} if rec2.get("ok") else None
        better = r2 is not None and r2["vert_rms_mm"] < rec1["vert_rms_mm"]
        out.update(decision="engine round " + ("kept (better)" if better else
                                                "rejected (not better: keep the first shot)"),
                   engine_round=r2 if r2 else {"error": rec2.get("termination")},
                   final=r2 if better else out["first_shot"], simulations=2)
        c2.log(f"fallback {pid}: LANDED - first shot {rec1['vert_rms_mm']:.3f} -> engine round "
               f"{(r2 or {}).get('vert_rms_mm', float('nan')):.3f} mm "
               f"({'kept' if better else 'rejected'})")
    CWF.mkdir(parents=True, exist_ok=True)
    f.write_text(json.dumps(out, indent=1, default=str) + "\n")


def cmd_all(args) -> None:
    ctx = cr.Ctx()
    todo = list(c2.TEST)
    lock = threading.Lock()
    c2.log(f"combined_fallback: tolerance {TOLERANCE_MM} mm, {c2.SLOTS} runs x "
           f"{c2.SIM_THREADS} threads")

    def worker() -> None:
        while True:
            with lock:
                if not todo:
                    return
                pid = todo.pop(0)
            try:
                one(ctx, pid)
            except Exception as exc:  # recorded, the others go on
                c2.log(f"fallback {pid}: ERROR {type(exc).__name__}: {exc}")

    ths = [threading.Thread(target=worker, daemon=True) for _ in range(c2.SLOTS)]
    for t in ths:
        t.start()
    for t in ths:
        t.join()
    cmd_report(args)
    (CWF / "done").write_text(json.dumps({"finished": c2.now_iso()}) + "\n")
    c2.log("combined_fallback: done")


def cmd_report(args) -> None:
    rows = []
    for p in c2.TEST:
        f = CWF / f"{p}.json"
        d = json.loads(f.read_text()) if f.exists() else {}
        rows.append({"point_id": p,
                     "first_shot_mm": (d.get("first_shot") or {}).get("vert_rms_mm", np.nan),
                     "engine_round_mm": (d.get("engine_round") or {}).get("vert_rms_mm", np.nan),
                     "final_mm": (d.get("final") or {}).get("vert_rms_mm", np.nan),
                     "simulations": d.get("simulations", np.nan),
                     "decision": d.get("decision", "not run"),
                     "data_source": "SparLab simulation"})
    df = pd.DataFrame(rows)
    OUTF.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTF / "results.csv", index=False)
    f3 = lambda x: "-" if pd.isna(x) else f"{x:.3f}"  # noqa: E731
    L = ["# ML first shot + automatic engine fallback (8 test parts)", "",
         f"**SparLab simulation**, DSIF + rim pass, penalty 10. Tolerance {TOLERANCE_MM} mm "
         "(vertical RMS): a part above it after the ML first shot (combined v1) gets ONE "
         "displacement-adjustment step from its simulated error, simulated again; the result "
         "is kept only if better. 8 parts: a small sample.", "",
         "| part | ML first shot | engine round | final | simulations | decision |",
         "|---|--:|--:|--:|--:|---|"]
    for r in rows:
        L.append(f"| {r['point_id']} | {f3(r['first_shot_mm'])} | {f3(r['engine_round_mm'])} | "
                 f"{f3(r['final_mm'])} | {r['simulations']} | {r['decision']} |")
    m = df.mean(numeric_only=True)
    L.append(f"| **mean** | **{f3(m['first_shot_mm'])}** | | **{f3(m['final_mm'])}** | "
             f"**{m['simulations']:.2f}** | |")
    L += ["", "Made by `python/scripts/combined_fallback.py`."]
    (OUTF / "README.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("what", choices=["all", "report"])
    a = ap.parse_args()
    {"all": cmd_all, "report": cmd_report}[a.what](a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
