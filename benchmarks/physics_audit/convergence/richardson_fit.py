#!/usr/bin/env python3
"""Robustness of the Richardson limit (added on resume).

The observed order of the IM (penalty 10) family depends on which three meshes
are used (p = 1.6 ... 4), i.e. the convergence is not cleanly asymptotic
(the tool path's 1 mm spacing is not commensurate with 0.833 mm elements, etc.).
This script (a) evaluates every monotone triple of the IM 1-layer meshes
(2 x 2 x 5 at 2 / 1 / 0.833 mm, 2 x 2 x 7 at 0.625 mm - the 5 -> 7 point change is
<= 0.003 mm, Table 2) and (b) least-squares fits Q = Q* + C h^2 over h <= 1 mm
and Q = Q* + C h^p over all four meshes. The spread of Q* over these estimates is
the uncertainty of the limit, and so of every "error" quoted against it.
Writes richardson_fit.csv. CONV_PART=c1 does the same for the second part.
"""
import itertools
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from richardson import scalar

HERE = Path(__file__).resolve().parent
PART = os.environ.get("CONV_PART", "t0")
SFX = "" if PART == "t0" else "_" + PART
R = pd.read_csv(HERE / f"results{SFX}.csv").set_index("case")
if PART == "t0":
    FAM = {2.0: "h2_L1_im_tp5", 1.0: "h1_L1_im_tp5", 40 / 48: "h0.833_L1_im_tp5",
           40 / 56: "h0.714_L1_im_tp5", 0.625: "h0.625_L1_im_tp7"}
    EVAL = ["h2_L2_std", "h1_L2_std", "h1_L1_im_tp5", "h0.833_L1_im_tp5",
            "h0.833_L1_im_tp5_p3", "h0.714_L1_im_tp5_p3", "h0.714_L1_im_tp5", "h0.625_L1_im_tp7"]
else:
    FAM = {2.0: f"{PART}_h2_L1_im_tp5", 1.0: f"{PART}_h1_L1_im_tp5",
           40 / 48: f"{PART}_h0.833_L1_im_tp5", 0.625: f"{PART}_h0.625_L1_im_tp5"}
    EVAL = ["bench_cache", f"{PART}_h2_L2_std", f"{PART}_h1_L1_im_tp5",
            f"{PART}_h0.833_L1_im_tp5", f"{PART}_h0.833_L1_im_tp5_p3", f"{PART}_h0.625_L1_im_tp5"]
FAM = {h: c for h, c in FAM.items() if c in R.index}
QS = ["vert_rms_mm", "vert_max_mm", "rim_sag_mm", "deep_bias_mm", "free_flange_bias_mm"]

rows = []
for q in QS:
    lims = {}
    hs = sorted(FAM, reverse=True)
    for tri in itertools.combinations(hs, 3):
        r = scalar([float(R.loc[FAM[h], q]) for h in tri], list(tri))
        if r["monotone"] and np.isfinite(r["p"]):
            lims["triple " + "/".join(f"{h:.3g}" for h in tri) + f" (p={r['p']:.2f})"] = r["q_ext"]
    fine = [h for h in hs if h <= 1.0]
    if len(fine) >= 2:
        A = np.column_stack([np.ones(len(fine)), np.array(fine) ** 2])
        y = np.array([float(R.loc[FAM[h], q]) for h in fine])
        c, *_ = np.linalg.lstsq(A, y, rcond=None)
        lims["LSQ p=2, h<=1"] = c[0]
    if len(hs) >= 4:
        y = np.array([float(R.loc[FAM[h], q]) for h in hs])
        h = np.array(hs)
        sol = least_squares(lambda x: x[0] + x[1] * h ** x[2] - y, [y[-1], (y[0] - y[-1]) / 4, 2.0])
        lims[f"LSQ free p (p={sol.x[2]:.2f}), all"] = sol.x[0]
    v = np.array(list(lims.values()))
    row = {"quantity": q, "limit_min": v.min(), "limit_max": v.max(),
           "limit_median": float(np.median(v)), "n_estimates": len(v),
           "estimates": "; ".join(f"{k}: {x:.4f}" for k, x in lims.items())}
    for c in EVAL:
        if c in R.index:
            val = float(R.loc[c, q])
            e = sorted([val - v.min(), val - v.max()], key=abs)
            row[f"err[{c}]"] = f"{e[0]:+.3f}..{e[1]:+.3f}"
    rows.append(row)
df = pd.DataFrame(rows)
df.to_csv(HERE / f"richardson_fit{SFX}.csv", index=False, float_format="%.4f")
with pd.option_context("display.width", 250, "display.max_colwidth", 400):
    for _, r in df.iterrows():
        print(r["quantity"], f"limit {r.limit_min:.3f}..{r.limit_max:.3f} (median {r.limit_median:.3f})")
        print("   ", r["estimates"])
        print("   ", {k: r[k] for k in df.columns if k.startswith("err[") and isinstance(r[k], str)})
