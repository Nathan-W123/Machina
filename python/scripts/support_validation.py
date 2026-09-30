#!/usr/bin/env python3
"""Rim support on one small cone: single-point forming vs a backing plate vs
double-sided incremental forming (DSIF), with and without the DSIF rim pass,
uncompensated and after one FE displacement-adjustment step.

usage:
    python3 python/scripts/support_validation.py --out benchmarks/support_cone \
        --work benchmarks/support_cone/work [--workers 3] [--stage uncompensated|da|all]

Every formed shape is a `sparlab_form` simulation (implicit, quasi-static,
finite_logarithmic kinematics, Hill48 AA5754-O with nominal handbook data):
**simulation, not experiment**. The part, blank, mesh, tool and path are
those of the springback benchmark (`benchmarks/springback`): its smallest
held-out cone, truncated_cone-s2026-0000. Resumable and idempotent: every
run goes through precomp's content-addressed cache in `<work>/runs` (an
identical deck is never run twice), so a rerun after an interruption picks up
where it stopped and rewrites the tables from the cache.

Per case the script records the deviation of the released part from the
target (vertical RMS / max / bias over the part, normal RMS; upper band -
the part less than 1 mm deep, where the rim sags - deeper part and
flange), the run's cost (runtime on one thread, increments, Newton
iterations, cut-backs) and its physics checks: every step completed, the
largest penetration and peak force of every tool, the plate's reaction
(the force the sheet puts on it), the reactions the clamp holds in the
unload steps and what the 3-2-1 release leaves (must be ~0), the warnings.

Outputs in --out: config.json (the protocol), cases.csv (one row per case),
regions.csv, tools.csv, profiles.csv (azimuthal means of target and formed
surfaces against the radius), tables.md, run.json (versions, timings, load).
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import copy
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

MM = 1e-3
PROTOCOL: Dict[str, Any] = {
    "protocol_revision": 1,
    "part": {"family": "truncated_cone", "point_id": "truncated_cone-s2026-0000",
             "top_radius": 0.007322111932560801, "wall_angle_deg": 38.61650181468576,
             "depth": 0.0030903661344200374, "top_fillet": 0.0010835243463516236,
             "bottom_fillet": 0.002076483284123242},
    "grid_spacing": 0.00025,
    "setup": {"blank_size": 0.04, "clamp_margin": 0.005, "element": "hex8",
              "element_size": 0.002, "friction": 0.1, "kinematics": "finite_logarithmic",
              "layers": 2, "material": "AA5754-O", "max_tool_travel": 0.001,
              "release": "321", "step_down": 0.001, "thickness": 0.001, "threads": 1,
              "tool_radius": 0.004, "toolpath_spacing": 0.001, "toolpath_style": "spiral"},
    "strategies": {
        "spif": {"support": "none"},
        "plate": {"support": "backing_plate", "support_settings": {"clearance": 0.001}},
        "dsif": {"support": "dsif", "support_settings": {}},
        "dsif_rim": {"support": "dsif", "support_settings": {"rim_pass": True}},
        # sensitivity: a gap of the initial thickness t instead of the sine
        # law's t cos(theta) - the model's wall comes out thicker than the
        # sine law, so the sine-law gap squeezes it
        "dsif_t0": {"support": "dsif", "support_settings": {"thickness_law": "initial"}},
    },
    "da": {"alpha": 1.0, "direction": "vertical", "smoothing": None,
           "max_wall_angle_deg": 65.0},
    "upper_band_depth": 0.001,
}


def log(msg: str) -> None:
    print(time.strftime("%H:%M:%S"), msg, flush=True)


def make_setup(cfg: Dict[str, Any], name: str, executable: Optional[str]):
    from precomp.fea import FormingSetup
    from precomp.materials import get_material

    doc = dict(cfg["setup"])
    doc["material"] = get_material(doc["material"])
    doc.update(copy.deepcopy(cfg["strategies"][name]))
    doc["executable"] = executable
    doc["timeout"] = 4 * 3600.0
    return FormingSetup(**doc)


def make_target(cfg: Dict[str, Any], setup):
    from precomp.geometry import Grid
    from precomp.geometry.parts import part_from_dict

    p = {k: v for k, v in cfg["part"].items() if k != "point_id"}
    part = part_from_dict(p)
    return part.heightmap(Grid.centered(setup.meshed_blank_size, cfg["grid_spacing"]))


def deviation(formed, target, band: float) -> Dict[str, float]:
    from precomp.metrology import (flange_mask, metrics, part_mask, signed_deviation,
                                   vertical_deviation)

    part = part_mask(target)
    v = vertical_deviation(formed, target)
    n = signed_deviation(formed, target)
    mv, mn = metrics(v, part), metrics(n, part)
    upper = part & (target.z > -band)
    out = {"vert_rms_mm": mv["rms"] / MM, "vert_max_mm": mv["max_abs"] / MM,
           "vert_bias_mm": mv["bias"] / MM, "normal_rms_mm": mn["rms"] / MM}
    tot = float(np.sum(v.z[part & v.mask] ** 2))
    for reg, m in {"upper_band": upper, "deep": part & ~upper,
                   "flange": flange_mask(target)}.items():
        sel = m & v.mask
        out[f"{reg}_rms_mm"] = float(np.sqrt(np.mean(v.z[sel] ** 2))) / MM
        out[f"{reg}_bias_mm"] = float(np.mean(v.z[sel])) / MM
        out[f"{reg}_max_mm"] = float(np.abs(v.z[sel]).max()) / MM
    out["upper_band_share"] = float(np.sum(v.z[upper & v.mask] ** 2)) / tot
    return out


def physics(res) -> Dict[str, Any]:
    """The run's cost and its physics checks, from summary.json and
    tool_forces.csv."""
    s = res.summary
    t = s.get("timing", {})
    steps = {st["name"]: st for st in s["steps"]}
    out: Dict[str, Any] = {
        "completed": bool(s["completed"]), "termination": s["termination"],
        "runtime_s": float(s["runtime_s"]), "increments": t.get("increments"),
        "iterations": t.get("iterations"), "cuts": t.get("cuts"),
        "max_plastic_strain": max(st.get("max_plastic_strain", 0.0) for st in s["steps"]),
        "warnings": len(s.get("warnings", [])),
        "warning_text": " | ".join(w[:160] for w in s.get("warnings", [])),
        "steps": ",".join(st["name"] for st in s["steps"]),
        "unload_clamp_reaction_N": steps.get("unload", {}).get("reaction_norm_N"),
        "release_reaction_N": steps.get("release", {}).get("reaction_norm_N"),
        "release_reference_force_N": steps.get("release", {}).get("reference_force_N"),
    }
    if "rim_unload" in steps:
        out["rim_unload_clamp_reaction_N"] = steps["rim_unload"].get("reaction_norm_N")
    return out


def tool_rows(case: str, res) -> List[Dict[str, Any]]:
    f = res.forming_forces()
    rows = []
    for tool in res.summary.get("tools", []):
        name = tool["name"]
        g = f[f["tool"] == name]
        row = {"case": case, "tool": name, "shape": tool.get("shape"),
               "max_active_nodes": tool.get("max_active_nodes"),
               "max_penetration_um": 1e6 * float(tool.get("max_penetration_m", 0.0)),
               "peak_force_N": float(np.linalg.norm(tool.get("peak_force_N", [0, 0, 0]))),
               "rows": int(len(g))}
        if len(g):
            touching = g[g["active_nodes"] > 0]
            row["increments_in_contact"] = int(len(touching))
            row["fz_max_N"] = float(g["fz"].max())
            row["fz_min_N"] = float(g["fz"].min())
            row["mean_abs_fz_in_contact_N"] = float(touching["fz"].abs().mean()) \
                if len(touching) else 0.0
            row["steps"] = ",".join(str(int(k)) for k in sorted(g["step"].unique()))
        rows.append(row)
    return rows


def profile(hm, bins: np.ndarray) -> np.ndarray:
    X, Y = hm.grid.mesh()
    r = np.hypot(X, Y)
    out = np.full(len(bins), np.nan)
    for i, b in enumerate(bins):
        sel = hm.mask & (np.abs(r - b) <= 0.125e-3)
        if sel.any():
            out[i] = float(hm.z[sel].mean())
    return out


def run_case(cfg: Dict[str, Any], name: str, work: Path, executable: Optional[str],
             da_steps: int) -> Dict[str, Any]:
    """Simulate the target (and, with `da_steps`, the FE-DA iterates) with
    strategy `name` through the cache; return the per-iterate records."""
    from precomp.compensation import FEAPredictor, displacement_adjustment
    from precomp.fea.support import command_upper_bound, compensation_masks

    setup = make_setup(cfg, name, executable)
    target = make_target(cfg, setup)
    pred = FEAPredictor(setup, work / "runs", target=target)
    upper = command_upper_bound(setup, target) if setup.support != "none" else None
    hold, adjust = compensation_masks(setup, target)
    records: List[Dict[str, Any]] = []

    def cb(k, commanded, formed, error):
        res = pred.results[-1]
        rec = {"case": name if k == 0 else f"{name}+DA{k}", "strategy": name,
               "fe_run": k + 1, "support": setup.support,
               "commanded_max_mm": float(commanded.z.max()) / MM,
               "commanded_depth_mm": commanded.depth / MM,
               "key": res.provenance.get("key"), "cache_hit": res.provenance.get("cache_hit"),
               "wall_runtime_s": res.provenance.get("runtime_s")}
        rec.update(deviation(formed, target, cfg["upper_band_depth"]))
        rec.update(physics(res))
        # the wall's thickness after forming against the sine law t cos(theta)
        # (the DSIF support's gap), and what the support reports
        theta = target.wall_angle()
        wall = (target.z < -3e-4) & (np.degrees(theta) > 30.0)
        th = res.thickness_map("form", grid=target.grid)
        sel = wall & th.mask
        rec["wall_thickness_mm"] = float(np.median(th.z[sel])) / MM
        rec["sine_law_thickness_mm"] = float(np.median(setup.thickness * np.cos(theta[sel]))) / MM
        prov = json.loads((res.directory.parent / "precomp_deck.json").read_text())
        plate = prov.get("support", {}).get("plate")
        if plate:
            rec["plate_realised_clearance_mm"] = plate["realised_clearance_m"] / MM
        rec["_tools"] = tool_rows(rec["case"], res)
        rec["_formed"] = formed
        rec["_commanded"] = commanded
        records.append(rec)
        log(f"{rec['case']}: RMS {rec['vert_rms_mm']:.3f} mm, upper band bias "
            f"{rec['upper_band_bias_mm']:+.3f} mm, {rec['runtime_s']:.0f} s"
            + (" (cached)" if rec["cache_hit"] else ""))

    da = cfg["da"]
    try:
        displacement_adjustment(target, pred, iterations=da_steps + 1, alpha=da["alpha"],
                                direction=da["direction"], smoothing=da["smoothing"],
                                max_wall_angle_deg=da["max_wall_angle_deg"],
                                upper_bound=upper, hold_mask=hold, adjust_mask=adjust,
                                callback=cb)
    except Exception as exc:  # recorded, reported
        records.append({"case": f"{name}+DA{len(records)}" if records else name,
                        "strategy": name, "fe_run": len(records) + 1, "completed": False,
                        "error": f"{type(exc).__name__}: {exc}".splitlines()[0][:400]})
        log(f"{name}: FAILED - {records[-1]['error']}")
    return {"name": name, "records": records, "target": target, "setup": setup}


def versions(executable: str) -> Dict[str, Any]:
    import scipy

    from precomp import __version__
    from precomp.fea import sparlab_version

    v = {"precomp": __version__, "python": platform.python_version(), "numpy": np.__version__,
         "scipy": scipy.__version__, "pandas": pd.__version__, "platform": platform.platform(),
         "sparlab_form": sparlab_version(executable)}
    try:
        v["git_revision"] = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                                           capture_output=True, text=True,
                                           cwd=Path(__file__).parent).stdout.strip()
    except OSError:
        pass
    return v


def fmt(x: Any, nd: int = 3) -> str:
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "-"
    return f"{x:.{nd}f}" if isinstance(x, (float, np.floating)) else str(x)


def write_tables(out: Path, cases: pd.DataFrame, tools: pd.DataFrame) -> None:
    lines = ["# Rim support on one small cone (SparLab simulations)", "",
             "Deviation of the released part from the target [mm]; rim sag = the upper "
             "band's bias (part < 1 mm deep; negative = too deep).", "",
             "| case | FE runs | vertical RMS | max | bias | normal RMS | upper band RMS | "
             "rim sag (upper band bias) | deep RMS | flange RMS | runtime [s] | increments | "
             "iterations | cuts | warnings |",
             "|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|"]
    for _, r in cases.iterrows():
        if not r.get("completed", False):
            lines.append(f"| {r['case']} | {r.get('fe_run')} | failed: {r.get('error', '')} "
                         "|||||||||||||")
            continue
        lines.append(
            f"| {r['case']} | {int(r['fe_run'])} | {fmt(r['vert_rms_mm'])} | "
            f"{fmt(r['vert_max_mm'])} | {fmt(r['vert_bias_mm'])} | {fmt(r['normal_rms_mm'])} | "
            f"{fmt(r['upper_band_rms_mm'])} | {fmt(r['upper_band_bias_mm'])} | "
            f"{fmt(r['deep_rms_mm'])} | {fmt(r['flange_rms_mm'])} | {fmt(r['runtime_s'], 0)} | "
            f"{int(r['increments'])} | {int(r['iterations'])} | {int(r['cuts'])} | "
            f"{int(r['warnings'])} |")
    lines += ["", "Tools (force of the sheet on the tool; fz < 0: pushed down):", "",
              "| case | tool | max nodes in contact | max penetration [um] | peak force [N] | "
              "fz min [N] | fz max [N] | steps |", "|---|---|--:|--:|--:|--:|--:|---|"]
    for _, r in tools.iterrows():
        lines.append(f"| {r['case']} | {r['tool']} | {r['max_active_nodes']} | "
                     f"{fmt(r['max_penetration_um'], 2)} | {fmt(r['peak_force_N'], 0)} | "
                     f"{fmt(r.get('fz_min_N'), 0)} | {fmt(r.get('fz_max_N'), 0)} | "
                     f"{r.get('steps', '')} |")
    (out / "tables.md").write_text("\n".join(lines) + "\n")


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--work", required=True, help="run cache (large, git-ignored)")
    ap.add_argument("--workers", type=int, default=3, help="concurrent simulations")
    ap.add_argument("--da-steps", type=int, default=1, help="FE-DA steps per strategy")
    ap.add_argument("--strategies", nargs="*", help="subset of the protocol's strategies")
    ap.add_argument("--executable", help="sparlab_form (default: build/bin/sparlab_form)")
    args = ap.parse_args(argv)
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    out, work = Path(args.out), Path(args.work)
    out.mkdir(parents=True, exist_ok=True)
    work.mkdir(parents=True, exist_ok=True)
    cfg_path = out / "config.json"
    if cfg_path.exists():
        cfg = json.loads(cfg_path.read_text())
    else:
        cfg = copy.deepcopy(PROTOCOL)
        cfg_path.write_text(json.dumps(cfg, indent=2, sort_keys=True) + "\n")
    names = args.strategies or list(cfg["strategies"])
    exe = args.executable or str(Path(__file__).resolve().parents[2] / "build" / "bin"
                                 / "sparlab_form")
    t0 = time.time()
    load0 = os.getloadavg()
    log(f"strategies {names}, {args.da_steps} DA step(s), {args.workers} at a time")
    # one DA chain per strategy; the chains run side by side, each run on one thread
    with cf.ThreadPoolExecutor(max_workers=args.workers) as pool:
        results = list(pool.map(lambda n: run_case(cfg, n, work, exe, args.da_steps), names))
    rows, trows, prof = [], [], {}
    bins = np.arange(0.0, 16.01, 0.5) * MM
    for r in results:
        if "target" in r and "target" not in prof:
            prof["r_mm"] = bins / MM
            prof["target_mm"] = profile(r["target"], bins) / MM
        for rec in r["records"]:
            trows += rec.pop("_tools", [])
            f = rec.pop("_formed", None)
            c = rec.pop("_commanded", None)
            if f is not None:
                prof[f"{rec['case']}_formed_mm"] = profile(f, bins) / MM
                prof[f"{rec['case']}_commanded_mm"] = profile(c, bins) / MM
            rows.append(rec)
    cases = pd.DataFrame(rows)
    tools = pd.DataFrame(trows)
    cases.to_csv(out / "cases.csv", index=False, float_format="%.6g")
    tools.to_csv(out / "tools.csv", index=False, float_format="%.6g")
    pd.DataFrame(prof).to_csv(out / "profiles.csv", index=False, float_format="%.5g")
    write_tables(out, cases, tools)
    meta = {"command": " ".join([Path(sys.argv[0]).name] + list(argv or sys.argv[1:])),
            "versions": versions(exe), "started_utc": time.strftime(
                "%Y-%m-%dT%H:%M:%SZ", time.gmtime(t0)),
            "elapsed_s": time.time() - t0, "load_average_start": load0,
            "load_average_end": os.getloadavg(), "workers": args.workers,
            "cpu_count": os.cpu_count(), "data_source": "SparLab simulation"}
    (out / "run.json").write_text(json.dumps(meta, indent=2, sort_keys=True, default=str) + "\n")
    log(f"done in {time.time() - t0:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
