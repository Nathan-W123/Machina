#!/usr/bin/env python3
"""Springback compensation by displacement adjustment with SparLab's forming
FEA - the first physics number of precomp.

usage:
    python3 python/scripts/fea_da_demo.py --out benchmarks/fea_da_cone/<run> \
        --work <scratch_dir> [--clamp-margin 0.009] [--threads 1]

It runs the `precomp` command exactly as a user would, and records what it
ran and what it measured:

    precomp part --family truncated_cone --param top_radius=0.009 \
        wall_angle_deg=45 depth=0.003 top_fillet=0.002 bottom_fillet=0.002 \
        --size 0.04 --spacing 0.00025 --out target.npz --json-out part.json
    precomp setup --material AA5754-O --set blank_size=0.04 \
        clamp_margin=<clamp margin> element_size=0.002 layers=2 thickness=0.001 \
        tool_radius=0.004 step_down=0.001 threads=<threads> --out setup.json
    precomp compensate --target target.npz --setup setup.json --method fea \
        --iterations 2 --work-dir <work> --out <work>/compensation

(every other setup field at its default: finite_logarithmic kinematics, a
spiral tool path with 1 mm spacing, 1 mm tool travel per increment,
friction 0.1, release onto 3-2-1 supports). Two sparlab_form simulations:
the target formed as it is (iteration 0) and the once-compensated shape
(iteration 1). It writes into --out: part.json, setup.json,
compensation.json (as `precomp compensate` wrote it), history.csv (per
iteration: the vertical and the normal deviation over the part, and the
cost of its simulation), profiles.csv (azimuthal means of the target, and
per iteration of the commanded shape and of the formed one after each step,
every 0.5 mm of radius) and run.json (commands, versions, runtimes).
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

import numpy as np  # noqa: E402

from precomp import __version__  # noqa: E402
from precomp.cli import main as precomp  # noqa: E402
from precomp.fea import load_result, sparlab_version  # noqa: E402
from precomp.fea.setup import FormingSetup  # noqa: E402
from precomp.geometry.heightmap import HeightMap  # noqa: E402
from precomp.metrology import metrics, part_mask, signed_deviation, vertical_deviation  # noqa: E402

PART = ["top_radius=0.009", "wall_angle_deg=45", "depth=0.003", "top_fillet=0.002",
        "bottom_fillet=0.002"]
GRID = ["--size", "0.04", "--spacing", "0.00025"]
SETUP = ["blank_size=0.04", "element_size=0.002", "layers=2", "thickness=0.001",
         "tool_radius=0.004", "step_down=0.001"]
ITERATIONS = 2


def mm(metric: dict) -> dict:
    return {k: (v * 1e3 if k != "n" else v) for k, v in metric.items()}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True, help="where the record goes")
    ap.add_argument("--work", required=True, help="scratch directory (run cache, maps)")
    ap.add_argument("--clamp-margin", type=float, default=0.009,
                    help="width of the clamped frame [m]: 0.009 leaves a 22 mm square "
                         "window around the 19.7 mm part, like a backing plate; 0.005 a "
                         "30 mm one")
    ap.add_argument("--threads", type=int, default=1, help="OpenMP threads per simulation")
    args = ap.parse_args()
    out = Path(args.out)
    work = Path(args.work)
    out.mkdir(parents=True, exist_ok=True)
    work.mkdir(parents=True, exist_ok=True)
    commands = [
        ["part", "--family", "truncated_cone", "--param", *PART, *GRID,
         "--out", str(work / "target.npz"), "--json-out", str(work / "part.json")],
        ["setup", "--material", "AA5754-O", "--set", *SETUP,
         f"clamp_margin={args.clamp_margin:g}", f"threads={args.threads}",
         "--out", str(work / "setup.json")],
        ["compensate", "--target", str(work / "target.npz"), "--setup",
         str(work / "setup.json"), "--method", "fea", "--iterations", str(ITERATIONS),
         "--work-dir", str(work), "--out", str(work / "compensation")],
    ]
    timings = {}
    started = time.time()
    for argv in commands:
        start = time.perf_counter()
        code = precomp(argv)
        timings[argv[0] + "_s"] = time.perf_counter() - start
        if code:
            print(f"precomp {' '.join(argv)} exited with {code}", file=sys.stderr)
            return code

    setup = FormingSetup.from_dict(json.loads((work / "setup.json").read_text()))
    target = HeightMap.load(work / "target.npz")
    part = part_mask(target)
    comp = json.loads((work / "compensation" / "compensation.json").read_text())
    rows = []
    # The runs of the cache, in the order DA made them: iteration k commanded c_k.
    entries = sorted((work / "runs").glob("*/*/COMPLETE"), key=lambda p: p.stat().st_mtime)
    if len(entries) != ITERATIONS:
        print(f"expected {ITERATIONS} simulations in the cache, found {len(entries)} "
              "(reuse of an old cache?)", file=sys.stderr)
    X, Y = target.grid.mesh()
    radius = np.hypot(X, Y)
    rings = np.arange(0.0, 0.5 * setup.meshed_blank_size - 1e-9, 5e-4)

    def ring_means(z):
        return [float(z[np.abs(radius - r) <= 1.25e-4].mean()) for r in rings]

    profiles = {"r_mm": [r * 1e3 for r in rings], "target_mm": ring_means(target.z * 1e3)}
    for k, (h, done) in enumerate(zip(comp["history"], entries)):
        entry = done.parent
        res = load_result(entry / "output")
        run = json.loads((entry / "run.json").read_text())
        s = res.summary
        formed = res.formed_surface(-1, grid=target.grid)
        commanded = HeightMap.load(entry / "commanded.npz")
        profiles[f"commanded_{k}_mm"] = ring_means(commanded.z * 1e3)
        for name in res.step_names:
            profiles[f"{name}_{k}_mm"] = ring_means(res.formed_surface(name, grid=target.grid).z
                                                    * 1e3)
        vert = metrics(vertical_deviation(formed, target), part)
        normal = metrics(signed_deviation(formed, target), part)
        assert abs(vert["rms"] - h["error"]["rms"]) < 1e-12, "cache order mismatch"
        forces = res.forming_forces()
        rows.append({
            "iteration": k, "key": entry.name[:16],
            "commanded_depth_mm": h["commanded_depth_m"] * 1e3,
            "formed_depth_mm": formed.depth * 1e3,
            **{f"vertical_{m}_mm": v for m, v in mm(vert).items() if m != "n"},
            **{f"normal_{m}_mm": v for m, v in mm(normal).items() if m != "n"},
            "part_nodes": vert["n"],
            "cache_hit": done.stat().st_mtime < started,
            "runtime_s": run["runtime_s"], "increments": s["timing"]["increments"],
            "newton_iterations": s["timing"]["iterations"], "cuts": s["timing"]["cuts"],
            "max_plastic_strain": max(st["max_plastic_strain"] for st in s["steps"]),
            "peak_tool_force_z_N": float(forces["fz"].max()),
            "release_reaction_N": s["steps"][-1]["reaction_norm_N"],
            "warnings": len(s["warnings"]),
        })
    with open(out / "profiles.csv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(list(profiles))
        writer.writerows(zip(*[[f"{v:.6g}" for v in col] for col in profiles.values()]))
    with open(out / "history.csv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    shutil.copy(work / "part.json", out / "part.json")
    shutil.copy(work / "setup.json", out / "setup.json")
    shutil.copy(work / "compensation" / "compensation.json", out / "compensation.json")
    record = {
        "created_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "data_source": "SparLab simulation (sparlab_form)",
        "commands": ["precomp " + " ".join(a).replace(str(work), "<work>") for a in commands],
        "script": "python3 python/scripts/fea_da_demo.py --out <out> --work <work> "
                  f"--clamp-margin {args.clamp_margin:g} --threads {args.threads}",
        "precomp_version": __version__,
        "sparlab_version": sparlab_version(setup.resolved_executable()),
        "clamp_margin_m": args.clamp_margin, "threads": args.threads,
        "cpu_count": os.cpu_count(),
        "load_average": list(os.getloadavg()),
        "timings_s": timings,
        "deviation_before_mm": {"vertical_rms": rows[0]["vertical_rms_mm"],
                                "vertical_max_abs": rows[0]["vertical_max_abs_mm"],
                                "normal_rms": rows[0]["normal_rms_mm"],
                                "normal_max_abs": rows[0]["normal_max_abs_mm"]},
        "deviation_after_mm": {"vertical_rms": rows[-1]["vertical_rms_mm"],
                               "vertical_max_abs": rows[-1]["vertical_max_abs_mm"],
                               "normal_rms": rows[-1]["normal_rms_mm"],
                               "normal_max_abs": rows[-1]["normal_max_abs_mm"]},
    }
    (out / "run.json").write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    for r in rows:
        print(f"iteration {r['iteration']}: vertical RMS {r['vertical_rms_mm']:.4f} mm, "
              f"max {r['vertical_max_abs_mm']:.4f} mm; normal RMS {r['normal_rms_mm']:.4f} mm, "
              f"max {r['normal_max_abs_mm']:.4f} mm; {r['increments']} increments, "
              f"{r['newton_iterations']} Newton iterations, {r['runtime_s']:.0f} s")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
