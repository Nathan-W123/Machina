#!/usr/bin/env python3
"""Markdown tables of the convergence study (from results.csv, diffs.csv,
richardson.csv) -> tables.md. The README quotes them."""
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
R = pd.read_csv(HERE / "results.csv").set_index("case")
D = pd.read_csv(HERE / "diffs.csv")
REF = "h0.625_L1_im_tp7"

LABEL = {
    "h2_L2_std": "2 mm, 2 std layers (**benchmark**)",
    "h2_L4_std": "2 mm, 4 std layers",
    "h2_L2_im": "2 mm, 2 IM layers 2x2x2",
    "h2_L1_im_tp5": "2 mm, 1 IM layer 2x2x5",
    "h2_L1_im_tp7": "2 mm, 1 IM layer 2x2x7",
    "h2_L2_std_tp5": "2 mm, 2 std layers 2x2x5",
    "h2_L2_std_t0.5": "2 mm, 2 std, travel 0.5 mm",
    "h2_L2_std_p30": "2 mm, 2 std, penalty 30",
    "h2_L2_std_p3": "2 mm, 2 std, penalty 3",
    "h2_L2_std_tol1e-8": "2 mm, 2 std, Newton tol 1e-8",
    "h2_L2_std_rel40": "2 mm, 2 std, 40 release increments",
    "h2_L2_std_x05": "2 mm, 2 std, explicit 0.5 m/s",
    "h1_L2_std": "1 mm, 2 std layers",
    "h1_L2_std_x05": "1 mm, 2 std, explicit 0.5 m/s",
    "h1_L1_im_tp5": "1 mm, 1 IM layer 2x2x5",
    "h1_L1_im_tp7": "1 mm, 1 IM layer 2x2x7",
    "h1_L2_im_tp3": "1 mm, 2 IM layers 2x2x3",
    "h0.833_L1_im_tp5": "0.833 mm, 1 IM layer 2x2x5",
    "h0.833_L1_im_tp5_p3": "0.833 mm, 1 IM layer 2x2x5, penalty 3",
    "h0.714_L1_im_tp5": "0.714 mm, 1 IM layer 2x2x5",
    "h0.714_L1_im_tp5_p3": "0.714 mm, 1 IM layer 2x2x5, penalty 3",
    "h0.625_L1_im_tp7": "0.625 mm, 1 IM layer 2x2x7 (**reference**)",
    "h0.625_L2_std_x05": "0.625 mm, 2 std, explicit 0.5 m/s",
}
ORDER = list(LABEL)


def f(x, nd=3):
    return "-" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{nd}f}"


def diff(a, b, col):
    q = D[(D.case == a) & (D.reference == b)]
    return float(q[col].iloc[0]) if len(q) and col in q else np.nan


out = []
out.append("### Table 1 - every case (truncated_cone-s2026-0000, released part vs target) [mm]\n")
out.append("| Case | DOFs | vert. RMS | vert. max | rim sag (upper-band bias) | deep-part bias | "
           "springback RMS (release - form) | springback max | formed depth | "
           "released vs reference RMS / max | wall [s] (load) |")
out.append("|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|")
import json  # noqa: E402
for c in ORDER:
    if c not in R.index:
        continue
    r = R.loc[c]
    rr = HERE / "cases" / c / "run.json"
    load = ""
    if rr.exists():
        j = json.loads(rr.read_text())
        load = f"{j['wall_s']:.0f} ({j['loadavg_start'][0]:.1f}-{j['loadavg_end'][0]:.1f})"
    dv = "0 (reference)" if c == REF else \
        f"{f(diff(c, REF, 'release_rms_part_mm'))} / {f(diff(c, REF, 'release_max_part_mm'))}"
    out.append(f"| {LABEL[c]} | {int(r.dofs)} | {f(r.vert_rms_mm)} | {f(r.vert_max_mm)} | "
               f"{f(r.rim_sag_mm)} | {f(r.deep_bias_mm)} | {f(r.sb_tot_rms_mm)} | "
               f"{f(r.sb_tot_max_mm)} | {f(r.depth_form_mm)} | {dv} | {load} |")

out.append("\n### Table 2 - one refinement at a time: change of each quantity [mm]\n")
out.append("| Refinement (from -> to) | d vert. RMS | d vert. max | d rim sag | d springback RMS | "
           "released-shape change RMS / max | springback-field change RMS / max | "
           "formed-shape change RMS |")
out.append("|---|--:|--:|--:|--:|--:|--:|--:|")
STEPS = [
    ("in-plane 2 -> 1 mm (2 std layers)", "h2_L2_std", "h1_L2_std"),
    ("in-plane 2 -> 1 mm (1 IM layer 2x2x7)", "h2_L1_im_tp7", "h1_L1_im_tp7"),
    ("in-plane 1 -> 0.833 mm (1 IM layer 2x2x5)", "h1_L1_im_tp5", "h0.833_L1_im_tp5"),
    ("in-plane 1 -> 0.625 mm (1 IM layer 2x2x7)", "h1_L1_im_tp7", "h0.625_L1_im_tp7"),
    ("layers 2 -> 4 (std, 2 mm)", "h2_L2_std", "h2_L4_std"),
    ("std 2 layers -> IM 2 layers (2 mm)", "h2_L2_std", "h2_L2_im"),
    ("std 2 layers -> IM 1 layer 2x2x5 (2 mm)", "h2_L2_std", "h2_L1_im_tp5"),
    ("std 2 layers -> IM 1 layer 2x2x7 (1 mm)", "h1_L2_std", "h1_L1_im_tp7"),
    ("IM 1 layer -> 2 layers 2x2x3 (1 mm)", "h1_L1_im_tp7", "h1_L2_im_tp3"),
    ("std 2x2x2 -> 2x2x5 per layer (2 mm)", "h2_L2_std", "h2_L2_std_tp5"),
    ("thickness points 5 -> 7 (IM, 2 mm)", "h2_L1_im_tp5", "h2_L1_im_tp7"),
    ("thickness points 5 -> 7 (IM, 1 mm)", "h1_L1_im_tp5", "h1_L1_im_tp7"),
    ("tool travel 1 -> 0.5 mm per increment", "h2_L2_std", "h2_L2_std_t0.5"),
    ("penalty scale 10 -> 30", "h2_L2_std", "h2_L2_std_p30"),
    ("penalty scale 3 -> 10", "h2_L2_std_p3", "h2_L2_std"),
    ("penalty scale 3 -> 10 (0.833 mm IM)", "h0.833_L1_im_tp5_p3", "h0.833_L1_im_tp5"),
    ("Newton tolerances 1e-6 -> 1e-8", "h2_L2_std", "h2_L2_std_tol1e-8"),
    ("release increments 10 -> 40", "h2_L2_std", "h2_L2_std_rel40"),
    ("implicit -> explicit 0.5 m/s (2 mm)", "h2_L2_std", "h2_L2_std_x05"),
    ("implicit -> explicit 0.5 m/s (1 mm)", "h1_L2_std", "h1_L2_std_x05"),
    ("benchmark -> reference", "h2_L2_std", REF),
    ("1 mm IM 2x2x5 -> reference", "h1_L1_im_tp5", REF),
    ("0.833 mm IM 2x2x5 -> reference", "h0.833_L1_im_tp5", REF),
    ("0.833 mm IM 2x2x5 penalty 3 -> reference", "h0.833_L1_im_tp5_p3", REF),
    ("in-plane 0.833 -> 0.714 mm (IM 2x2x5, penalty 3)", "h0.833_L1_im_tp5_p3", "h0.714_L1_im_tp5_p3"),
    ("penalty scale 3 -> 10 (0.714 mm IM)", "h0.714_L1_im_tp5_p3", "h0.714_L1_im_tp5"),
    ("0.714 mm IM 2x2x5 penalty 3 -> reference", "h0.714_L1_im_tp5_p3", REF),
    ("implicit -> explicit 0.5 m/s (0.625 mm; std 2 layers vs IM reference)", REF, "h0.625_L2_std_x05"),
]
for lab, a, b in STEPS:
    if a not in R.index or b not in R.index:
        continue
    ra, rb = R.loc[a], R.loc[b]
    q = D[((D.case == a) & (D.reference == b))]
    sgn = 1
    if q.empty:
        q = D[((D.case == b) & (D.reference == a))]
    if q.empty:
        cells = "- | - | -"
    else:
        q = q.iloc[0]
        cells = (f"{f(q.release_rms_part_mm)} / {f(q.release_max_part_mm)} | "
                 f"{f(q.sb_tot_rms_part_mm)} / {f(q.sb_tot_max_part_mm)} | "
                 f"{f(q.form_rms_part_mm)}")
    out.append(f"| {lab} | {rb.vert_rms_mm - ra.vert_rms_mm:+.3f} | "
               f"{rb.vert_max_mm - ra.vert_max_mm:+.3f} | {rb.rim_sag_mm - ra.rim_sag_mm:+.3f} | "
               f"{rb.sb_tot_rms_mm - ra.sb_tot_rms_mm:+.3f} | {cells} |")

rp = HERE / "richardson.csv"
if rp.exists():
    rc = pd.read_csv(rp)
    rc = rc[rc.family == "IM 1 layer 2x2x7"]
    out.append("\n### Table 3 - Richardson extrapolation, 1 IM layer 2x2x7, h = 2 / 1 / 0.625 mm\n")
    out.append("| Quantity | values at 2 / 1 / 0.625 mm | observed order p | extrapolated (h -> 0) | "
               "error at 2 mm | at 1 mm | at 0.625 mm |")
    out.append("|---|---|--:|--:|--:|--:|--:|")
    for _, r in rc.iterrows():
        out.append(f"| {r.quantity} | {r['values']} | {f(r.p, 2)} | {f(r.q_ext)} | "
                   f"{f(r.err_h1)} | {f(r.err_h2)} | {f(r.err_h3)} |")

if rp.exists():
    rc = pd.read_csv(rp)
    rc = rc[rc.family == "IM 1 layer 2x2x7"].set_index("quantity")
    lim = {q: float(rc.loc[q, "q_ext"]) for q in ("vert_rms_mm", "vert_max_mm", "rim_sag_mm")}
    fe_ref = float(rc.loc["release field RMS diff (part)", "err_h3"])
    out.append("\n### Table 4 - estimated discretisation error of each setup against the "
               "Richardson limit [mm]\n")
    out.append(f"Limits (h -> 0): vertical RMS {lim['vert_rms_mm']:.3f}, vertical max "
               f"{lim['vert_max_mm']:.3f}, rim sag {lim['rim_sag_mm']:.3f}. Released-shape "
               f"field error = RMS difference to the reference plus the reference's own "
               f"Richardson error ({fe_ref:.3f}), an upper bound (the two add when the errors "
               "have one shape, as in the IM family).\n")
    out.append("| Setup | error of vert. RMS | error of vert. max | error of rim sag | "
               "released-shape field error RMS (bound) | wall [s] (load) | cost / benchmark run |")
    out.append("|---|--:|--:|--:|--:|--:|--:|")
    base = json.loads((HERE / "cases" / "h2_L2_std" / "run.json").read_text())["wall_s"]
    for c in ORDER:
        if c not in R.index or c in ("h2_L2_std_t0.5", "h2_L2_std_p30", "h2_L2_std_p3",
                                      "h2_L2_std_tol1e-8", "h2_L2_std_rel40"):
            continue
        r = R.loc[c]
        dref = 0.0 if c == REF else diff(c, REF, "release_rms_part_mm")
        j = json.loads((HERE / "cases" / c / "run.json").read_text())
        out.append(f"| {LABEL[c]} | {r.vert_rms_mm - lim['vert_rms_mm']:+.3f} | "
                   f"{r.vert_max_mm - lim['vert_max_mm']:+.3f} | "
                   f"{r.rim_sag_mm - lim['rim_sag_mm']:+.3f} | {dref + fe_ref:.3f} | "
                   f"{j['wall_s']:.0f} ({j['loadavg_start'][0]:.1f}-{j['loadavg_end'][0]:.1f}) | "
                   f"{j['wall_s'] / base:.1f} |")

    out.append("\nThe explicit rows look better than their implicit twins only because the "
               "explicit runs' smaller springback (finding 5) happens to offset part of the "
               "mesh error on this part; that is not a reason to prefer them. Wall times "
               "were measured at different load averages (in brackets) on a shared machine "
               "and are indicative to about +-50 %.")

(HERE / "tables.md").write_text("\n".join(out) + "\n")
print("\n".join(out))
