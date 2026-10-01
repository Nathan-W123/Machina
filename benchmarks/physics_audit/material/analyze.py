#!/usr/bin/env python3
"""Tables of the material lens from the finished runs (runs/<part>/<variant>).

For every variant of a part, against that part's reference ("base", or
"base_im5" for variants ending in _im5):
  d_release : released shape of the variant minus released shape of the
              reference, over the part (target below the sheet plane):
              RMS, max |.|, mean; and RMS in the rim band (target < 1 mm deep)
              and the deeper part  [mm]
  d_loaded  : the same for the loaded shape (end of "form", tool in place)
  d_springback : change of the free springback (released - loaded), RMS [mm]
  dev_rms   : the variant's RMS vertical deviation from the target [mm]
              (the benchmark's metric) and its change
Writes sensitivity.csv and sensitivity.md.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs"


def surfaces(part, variant):
    f = RUNS / part / variant / "surfaces.npz"
    return np.load(f) if f.exists() else None


def stats(v):
    return float(np.sqrt(np.mean(v ** 2))), float(np.abs(v).max()), float(np.mean(v))


def main():
    recs = {}
    for line in (HERE / "results.jsonl").read_text().splitlines():
        r = json.loads(line)
        if r.get("ok"):
            recs[(r["part"], r["variant"])] = r          # last record wins
    rows = []
    for (part, var), r in sorted(recs.items()):
        ref = "base_im5" if var.endswith("_im5") else "base"
        S, B = surfaces(part, var), surfaces(part, ref)
        if S is None or B is None:
            continue
        t = S["target"]
        partm = t < -1e-6
        upper = partm & (t > -1e-3)
        deep = partm & ~upper
        mm = 1e3
        dr = (S["release"] - B["release"]) * mm
        dl = (S["form"] - B["form"]) * mm
        dsb = ((S["release"] - S["form"]) - (B["release"] - B["form"])) * mm
        rr, rmax, rmean = stats(dr[partm])
        rows.append({
            "part": part, "variant": var, "reference": ref,
            "d_release_rms_mm": rr, "d_release_max_mm": rmax, "d_release_mean_mm": rmean,
            "d_release_rim_rms_mm": stats(dr[upper])[0], "d_release_deep_rms_mm": stats(dr[deep])[0],
            "d_release_flange_rms_mm": stats(dr[~partm])[0],
            "d_loaded_rms_mm": stats(dl[partm])[0], "d_loaded_mean_mm": stats(dl[partm])[2],
            "d_springback_rms_mm": stats(dsb[partm])[0], "d_springback_mean_mm": stats(dsb[partm])[2],
            "springback_rms_mm": r["springback_part"]["rms"],
            "springback_mean_mm": r["springback_part"]["bias"],
            "springback_change_pct": 100.0 * (r["springback_part"]["rms"]
                                              / recs[(part, ref)]["springback_part"]["rms"] - 1.0)
            if (part, ref) in recs else np.nan,
            "dev_rms_mm": r["dev_part"]["rms"], "dev_bias_mm": r["dev_part"]["bias"],
            "dev_rim_rms_mm": r["dev_upper"]["rms"], "dev_deep_rms_mm": r["dev_deep"]["rms"],
            "d_dev_rms_mm": r["dev_part"]["rms"] - recs[(part, ref)]["dev_part"]["rms"]
            if (part, ref) in recs else np.nan,
            "released_depth_mm": r["released_depth_mm"], "loaded_depth_mm": r["loaded_depth_mm"],
            "peak_fz_N": r["peak_fz_N"], "max_eq_plastic_strain": r["max_plastic_strain"],
            "min_thickness_mm": r["min_thickness_mm"], "runtime_s": r["runtime_s"],
            "increments": r["increments"], "cuts": r["cuts"],
        })
    df = pd.DataFrame(rows)
    df.to_csv(HERE / "sensitivity.csv", index=False, float_format="%.4f")
    lines = []
    for part, g in df.groupby("part"):
        g = g.sort_values("d_release_rms_mm", ascending=False)
        lines += [f"### {part}", "",
                  "| variant | ref | d released RMS | max | mean | rim RMS | deep RMS | "
                  "d loaded RMS | d springback RMS | springback RMS (change) | dev RMS (d) | "
                  "peak Fz [N] | max eps_p | runtime [s] |",
                  "|---|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|"]
        for _, r in g.iterrows():
            lines.append(
                f"| {r.variant} | {r.reference} | {r.d_release_rms_mm:.3f} | {r.d_release_max_mm:.3f} | "
                f"{r.d_release_mean_mm:+.3f} | {r.d_release_rim_rms_mm:.3f} | "
                f"{r.d_release_deep_rms_mm:.3f} | {r.d_loaded_rms_mm:.3f} | "
                f"{r.d_springback_rms_mm:.3f} | {r.springback_rms_mm:.3f} ({r.springback_change_pct:+.0f} %) | "
                f"{r.dev_rms_mm:.3f} ({r.d_dev_rms_mm:+.3f}) | {r.peak_fz_N:.0f} | "
                f"{r.max_eq_plastic_strain:.3f} | {r.runtime_s:.0f} |")
        lines.append("")
    (HERE / "sensitivity.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
