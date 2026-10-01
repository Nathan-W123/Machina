#!/usr/bin/env python3
"""Richardson estimates of the in-plane discretisation error (convergence study).

Scalars: Q(h) = Q* + C h^p over three meshes h1 > h2 > h3 (any ratios):
p solves (Q1 - Q2) / (Q2 - Q3) = (h1^p - h2^p) / (h2^p - h3^p); valid only
when the differences have one sign (monotone convergence) and p comes out in a
plausible range (0.5 - 4). Fields: the same with the norms of the differences
of the released surfaces on the target grid, ||z1 - z2|| and ||z2 - z3||
(over the part), assuming e(h) = C h^p with one error shape.

Writes richardson.csv.
"""
from pathlib import Path
import sys

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


def solve_p(d12, d23, h1, h2, h3):
    """p with (h1^p - h2^p) / (h2^p - h3^p) = d12 / d23 (d's > 0)."""
    if not (d12 > 0 and d23 > 0):
        return np.nan
    target = d12 / d23
    f = lambda p: (h1 ** p - h2 ** p) / (h2 ** p - h3 ** p) - target  # noqa: E731
    lo, hi = 0.05, 8.0
    if f(lo) * f(hi) > 0:
        return np.nan
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def scalar(q, h):
    (q1, q2, q3), (h1, h2, h3) = q, h
    d12, d23 = q1 - q2, q2 - q3
    if d12 * d23 <= 0:
        return {"p": np.nan, "q_ext": np.nan, "err_h3": np.nan, "err_h2": np.nan,
                "err_h1": np.nan, "monotone": False}
    p = solve_p(abs(d12), abs(d23), h1, h2, h3)
    C = d23 / (h2 ** p - h3 ** p)
    qx = q3 - C * h3 ** p
    return {"p": p, "q_ext": qx, "err_h1": q1 - qx, "err_h2": q2 - qx, "err_h3": q3 - qx,
            "monotone": True}


def main():
    import os
    part = os.environ.get("CONV_PART", "t0")
    sfx = "" if part == "t0" else "_" + part
    res = pd.read_csv(HERE / f"results{sfx}.csv").set_index("case")
    rows = []
    fams = {
        "IM 1 layer 2x2x7": (["h2_L1_im_tp7", "h1_L1_im_tp7", "h0.625_L1_im_tp7"], [2, 1, 0.625]),
        "IM 1 layer 2x2x7 (h 2, 1, 0.5)": (["h2_L1_im_tp7", "h1_L1_im_tp7", "h0.5_L1_im_tp7"],
                                            [2, 1, 0.5]),
        "IM 1 layer 2x2x5 (h 2, 1, 0.833)": (["h2_L1_im_tp5", "h1_L1_im_tp5", "h0.833_L1_im_tp5"],
                                              [2, 1, 40 / 48]),
    }
    if part != "t0":
        fams = {
            "c1: IM 1 layer 2x2x5": ([f"{part}_h2_L1_im_tp5", f"{part}_h1_L1_im_tp5",
                                      f"{part}_h0.625_L1_im_tp5"], [2, 1, 0.625]),
            "c1: IM 1 layer 2x2x5 (h 2, 1, 0.833)": ([f"{part}_h2_L1_im_tp5", f"{part}_h1_L1_im_tp5",
                                                     f"{part}_h0.833_L1_im_tp5"], [2, 1, 40 / 48]),
            "c1: IM 1 layer 2x2x5 (h 1, 0.833, 0.625)": ([f"{part}_h1_L1_im_tp5",
                                                         f"{part}_h0.833_L1_im_tp5",
                                                         f"{part}_h0.625_L1_im_tp5"],
                                                        [1, 40 / 48, 0.625]),
        }
    qs = ["vert_rms_mm", "vert_max_mm", "rim_sag_mm", "deep_bias_mm", "free_flange_bias_mm",
          "depth_release_mm", "sb_tot_rms_mm", "sb_rel_rms_mm", "fz_mean_N"]
    for fam, (cs, h) in fams.items():
        if not all(c in res.index for c in cs):
            continue
        for q in qs:
            vals = [float(res.loc[c, q]) for c in cs]
            r = scalar(vals, h)
            rows.append({"family": fam, "quantity": q, "h": "/".join(map(str, h)),
                         "values": " / ".join(f"{v:.4f}" for v in vals), **r})
        # fields
        S = {c: np.load(HERE / "surfaces" / f"{c}.npz") for c in cs}
        from analyze import PART  # noqa: E402
        for key in ("release", "sb_tot"):
            def fld(c):
                s = S[c]
                return s["release"] if key == "release" else s["release"] - s["form"]
            d12 = np.sqrt(np.nanmean(((fld(cs[0]) - fld(cs[1]))[PART]) ** 2)) * 1e3
            d23 = np.sqrt(np.nanmean(((fld(cs[1]) - fld(cs[2]))[PART]) ** 2)) * 1e3
            p = solve_p(d12, d23, *h)
            C = d23 / (h[1] ** p - h[2] ** p) if np.isfinite(p) else np.nan
            rows.append({"family": fam, "quantity": f"{key} field RMS diff (part)",
                         "h": "/".join(map(str, h)),
                         "values": f"d(h1,h2) {d12:.4f}, d(h2,h3) {d23:.4f}", "p": p,
                         "err_h1": C * h[0] ** p, "err_h2": C * h[1] ** p,
                         "err_h3": C * h[2] ** p, "monotone": np.isfinite(p)})
    df = pd.DataFrame(rows)
    df.to_csv(HERE / f"richardson{sfx}.csv", index=False, float_format="%.5g")
    with pd.option_context("display.width", 250, "display.max_columns", 20):
        print(df.to_string(index=False, float_format=lambda v: f"{v:.4f}"))


if __name__ == "__main__":
    main()
