#!/usr/bin/env python3
"""Figures of the convergence study (from profiles.csv and diffs.csv)."""
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def style(ax):
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2)


def profiles(cases, labels, out, field="release_dev_mm", title=None):
    p = pd.read_csv(HERE / "profiles.csv")
    fig, ax = plt.subplots(figsize=(8, 4.2), dpi=130)
    for i, (c, lab) in enumerate(zip(cases, labels)):
        q = p[p.case == c]
        if q.empty:
            continue
        ax.plot(q.r_mm, q[field], color=SERIES[i], lw=2, label=lab)
    t = p[p.case == cases[0]]
    ax2 = ax.twinx() if False else None
    ax.axvline(15.0, color=INK2, lw=0.8, ls=":")
    ax.text(14.9, ax.get_ylim()[0], "clamp", color=INK2, ha="right", va="bottom", fontsize=8)
    ax.set_xlabel("radius r [mm]", color=INK)
    ax.set_ylabel({"release_dev_mm": "released - target, vertical [mm]",
                   "sb_tot_mm": "released - formed (tool on) [mm]"}.get(field, field), color=INK)
    style(ax)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK)
    if title:
        ax.set_title(title, color=INK, fontsize=10, loc="left")
    # target outline for orientation, scaled into the axes as a thin grey line
    if field == "release_dev_mm":
        ax.plot(t.r_mm, t.target_z_mm * 0.1, color="#b5b3ad", lw=1)
        ax.text(0.2, float(t.target_z_mm.min()) * 0.1, "target z / 10", color=INK2, fontsize=8,
                va="bottom")
    fig.tight_layout()
    fig.savefig(HERE / out)
    print("wrote", out)


def error_profiles(cases, labels, ref, out, title=None):
    """Azimuthal mean of (released surface of case - released surface of ref)."""
    p = pd.read_csv(HERE / "profiles.csv")
    r = p[p.case == ref].set_index("r_mm")
    fig, ax = plt.subplots(figsize=(8, 4.2), dpi=130)
    for i, (c, lab) in enumerate(zip(cases, labels)):
        q = p[p.case == c].set_index("r_mm")
        if q.empty:
            continue
        ax.plot(q.index, q.release_dev_mm - r.release_dev_mm, color=SERIES[i], lw=2, label=lab)
    ax.axhline(0.0, color=INK2, lw=0.8)
    for y in (-0.05, 0.05):
        ax.axhline(y, color=INK2, lw=0.8, ls="--")
    ax.text(0.2, 0.052, "+-0.05 mm", color=INK2, fontsize=8, va="bottom")
    ax.axvline(15.0, color=INK2, lw=0.8, ls=":")
    ax.set_xlabel("radius r [mm]", color=INK)
    ax.set_ylabel(f"released shape minus {ref} [mm]", color=INK)
    style(ax)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK)
    if title:
        ax.set_title(title, color=INK, fontsize=10, loc="left")
    fig.tight_layout()
    fig.savefig(HERE / out)
    print("wrote", out)


if __name__ == "__main__":
    ref = sys.argv[1] if len(sys.argv) > 1 else "h0.625_L1_im_tp7"
    profiles(["h2_L2_std", "h1_L2_std", "h2_L1_im_tp7", "h1_L1_im_tp7", ref],
             ["2 mm, 2 std layers (benchmark)", "1 mm, 2 std layers", "2 mm, 1 IM layer 2x2x7",
              "1 mm, 1 IM layer 2x2x7", "0.625 mm, 1 IM layer 2x2x7 (reference)"],
             "fig_profile_release.png",
             title="Truncated cone t-0000: azimuthal mean of the released part's deviation")
    error_profiles(["h2_L2_std", "h2_L1_im_tp7", "h1_L2_std", "h1_L1_im_tp7", "h2_L2_std_x05"],
                   ["2 mm, 2 std layers (benchmark)", "2 mm, 1 IM layer 2x2x7",
                    "1 mm, 2 std layers", "1 mm, 1 IM layer 2x2x7",
                    "2 mm, 2 std layers, explicit 0.5 m/s"], ref, "fig_profile_error.png",
                   title="Discretisation error of the released shape (azimuthal mean)")
