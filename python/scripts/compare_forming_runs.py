"""Compare two sparlab_form runs of the same deck: an explicit forming run
(`form_explicit` steps) against the implicit reference (`form` steps).

Both runs must share the mesh and the step names. Reported (in mm and N):

  * the shape of the tool-side surface (the nodes on the top face, z = z_max)
    at the end of the last step - after the release onto 3-2-1 supports -
    as the largest and the RMS distance between the deformed positions of
    the two runs, next to the springback of the reference (its largest nodal
    displacement over the release step), since the springback is what the
    pre-compensation must predict;
  * the springback itself (displacement over the release step) of both runs
    and the largest / RMS difference between them;
  * the shape at the end of the forming step (before release);
  * the tool force history: explicit forces are noisy (the penalty contact
    rattles at the scale of the time step), so both histories are resampled
    on the pseudo-time of the step, averaged over a moving window (default
    0.1 s of pseudo-time, a quarter of a contour here) and compared as the
    RMS and largest difference of the moving averages over the contact
    phase, relative to the reference's mean force, and as the difference of
    the mean forces.

usage: python3 compare_forming_runs.py <reference_dir> <explicit_dir>
       [--form-step form] [--release-step unclamp] [--window 0.1] [--json out.json]
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np


def read_nodes(path: Path) -> np.ndarray:
    """Rows node, X, Y, Z, ux, uy, uz."""
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    return data[np.argsort(data[:, 0])]


def step_file(run: Path, name: str, suffix: str) -> Path:
    hits = sorted(run.glob(f"step_*_{name}{suffix}"))
    hits = [h for h in hits if "_inc_" not in h.name]
    if not hits:
        raise SystemExit(f"{run}: no step file for step '{name}' ({suffix})")
    return hits[-1]


def tool_forces(run: Path, step_index: int):
    """(t, fz, active nodes) of the rows of step `step_index` (1-based)."""
    t, fz, fn = [], [], []
    with open(run / "tool_forces.csv", newline="") as f:
        for row in csv.DictReader(f):
            if int(row["step"]) != step_index:
                continue
            t.append(float(row["t"]))
            fz.append(float(row["fz"]))
            fn.append(int(row["active_nodes"]))
    return np.array(t), np.array(fz), np.array(fn)


def moving_average(t: np.ndarray, f: np.ndarray, grid: np.ndarray, window: float) -> np.ndarray:
    """The average of the piecewise-linear f(t) over [g - w/2, g + w/2]."""
    fine = np.linspace(grid[0] - window, grid[-1] + window, 20 * len(grid) + 1)
    values = np.interp(fine, t, f)
    cumulative = np.concatenate([[0.0], np.cumsum(0.5 * (values[1:] + values[:-1]) *
                                                  np.diff(fine))])
    hi = np.interp(grid + 0.5 * window, fine, cumulative)
    lo = np.interp(grid - 0.5 * window, fine, cumulative)
    return (hi - lo) / window


def step_index(run: Path, name: str) -> int:
    summary = json.loads((run / "summary.json").read_text())
    for k, s in enumerate(summary["steps"]):
        if s["name"] == name:
            return k + 1
    raise SystemExit(f"{run}: no step '{name}'")


def compare(reference: Path, explicit: Path, form: str, release: str, window: float) -> dict:
    out: dict = {"reference": str(reference), "explicit": str(explicit)}
    ref_end = read_nodes(step_file(reference, release, "_nodes.csv"))
    exp_end = read_nodes(step_file(explicit, release, "_nodes.csv"))
    ref_form = read_nodes(step_file(reference, form, "_nodes.csv"))
    exp_form = read_nodes(step_file(explicit, form, "_nodes.csv"))
    # The state the release starts from: the step before it.
    summary = json.loads((reference / "summary.json").read_text())
    names = [s["name"] for s in summary["steps"]]
    before = names[names.index(release) - 1]
    ref_before = read_nodes(step_file(reference, before, "_nodes.csv"))
    exp_before = read_nodes(step_file(explicit, before, "_nodes.csv"))
    if ref_end.shape != exp_end.shape or np.any(ref_end[:, :4] != exp_end[:, :4]):
        raise SystemExit("the two runs do not share the mesh")
    top = np.isclose(ref_end[:, 3], ref_end[:, 3].max())
    mm = 1.0e3

    def distance(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        return np.linalg.norm(a[top, 4:7] - b[top, 4:7], axis=1)

    def stats(d: np.ndarray) -> dict:
        return {"max_mm": float(d.max() * mm), "rms_mm": float(math.sqrt(np.mean(d ** 2)) * mm)}

    ref_spring = ref_end[:, 4:7] - ref_before[:, 4:7]
    exp_spring = exp_end[:, 4:7] - exp_before[:, 4:7]
    spring_mag = np.linalg.norm(ref_spring[top], axis=1)
    out["top_nodes"] = int(top.sum())
    out["reference_springback"] = {"max_mm": float(spring_mag.max() * mm),
                                   "rms_mm": float(math.sqrt(np.mean(spring_mag ** 2)) * mm)}
    exp_mag = np.linalg.norm(exp_spring[top], axis=1)
    out["explicit_springback"] = {"max_mm": float(exp_mag.max() * mm),
                                  "rms_mm": float(math.sqrt(np.mean(exp_mag ** 2)) * mm)}
    out["final_shape_difference"] = stats(distance(ref_end, exp_end))
    out["springback_difference"] = stats(np.linalg.norm(ref_spring[top] - exp_spring[top], axis=1))
    out["formed_shape_difference"] = stats(distance(ref_form, exp_form))
    depth = -ref_form[top, 6].min()
    out["reference_depth_mm"] = float(depth * mm)
    out["explicit_depth_mm"] = float(-exp_form[top, 6].min() * mm)

    # Tool forces of the forming step.
    k_ref = step_index(reference, form)
    k_exp = step_index(explicit, form)
    t_r, f_r, n_r = tool_forces(reference, k_ref)
    t_e, f_e, n_e = tool_forces(explicit, k_exp)
    lo = max(t_r[n_r > 0].min(), t_e[n_e > 0].min()) + window
    hi = min(t_r.max(), t_e.max()) - window
    grid = np.linspace(lo, hi, 400)
    a_r = moving_average(t_r, f_r, grid, window)
    a_e = moving_average(t_e, f_e, grid, window)
    mean_r = float(np.mean(a_r))
    mean_e = float(np.mean(a_e))
    out["tool_force"] = {
        "window_s": window,
        "pseudo_time_s": [float(lo), float(hi)],
        "reference_mean_N": mean_r,
        "explicit_mean_N": mean_e,
        "mean_difference": (mean_e - mean_r) / mean_r,
        "moving_average_rms_difference": float(math.sqrt(np.mean((a_e - a_r) ** 2)) / mean_r),
        "moving_average_max_difference": float(np.max(np.abs(a_e - a_r)) / mean_r),
        "explicit_raw_rms_scatter": float(
            math.sqrt(np.mean((np.interp(grid, t_e, f_e) - a_e) ** 2)) / mean_r),
    }
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("reference", type=Path)
    p.add_argument("explicit", type=Path)
    p.add_argument("--form-step", default="form")
    p.add_argument("--release-step", default="unclamp")
    p.add_argument("--window", type=float, default=0.1, help="moving average [s pseudo-time]")
    p.add_argument("--json", type=Path)
    a = p.parse_args(argv)
    result = compare(a.reference, a.explicit, a.form_step, a.release_step, a.window)
    text = json.dumps(result, indent=2)
    print(text)
    if a.json:
        a.json.write_text(text + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
