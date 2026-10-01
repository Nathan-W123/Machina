#!/usr/bin/env python3
"""Material share of the released-shape error, from sensitivity.csv.

Groups the one-at-a-time variants by the material input they probe, takes
the largest released-shape change (RMS over the part) per group and part,
and sums the groups in quadrature (a conservative bound: the groups are not
independent - the locus group contains von Mises and the stress-ratio bound;
signs are ignored). The literature combination (`combo_lit*`) is the
realistic single estimate. Writes error_budget.json.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
GROUPS = {
    "yield locus (vm, r sets, stress-ratio bound)": ["vm", "hill_iadicola", "hill_coer", "hill_sb110"],
    "Bauschinger": ["kin_shutov", "kin_mild"],
    "hardening beyond tests": ["hard_voce_coer", "hard_pow_iadicola"],
    "E (constant -10 %)": ["E63"],
    "E(eps_p)": ["Edeg"],
    "initial yield +-15 MPa": ["ys85", "ys115"],
    "friction 0.05-0.2": ["mu005", "mu02"],
    "sheet thickness -5 %": ["t095"],
}


def main():
    d = pd.read_csv(HERE / "sensitivity.csv")
    out = {}
    for part, g in d.groupby("part"):
        for elem, suffix in (("standard", ""), ("im5", "_im5")):
            rows = {}
            for name, vs in GROUPS.items():
                sel = g[g.variant.isin([v + suffix for v in vs])]
                if elem == "im5" and sel.empty:          # fall back to the standard element
                    sel = g[g.variant.isin(vs)]
                if len(sel):
                    rows[name] = float(sel.d_release_rms_mm.max())
            if elem == "im5" and not any(v.endswith("_im5") for v in g.variant):
                continue
            q = float(np.sqrt(sum(v ** 2 for v in rows.values())))
            combo = g[g.variant == "combo_lit" + suffix]
            out[f"{part} [{elem}]"] = {
                "largest_per_group_mm": rows, "quadrature_sum_mm": q,
                "combo_lit_mm": float(combo.d_release_rms_mm.iloc[0]) if len(combo) else None,
                "base_deviation_rms_mm": float(g[g.variant == "base" + suffix].dev_rms_mm.iloc[0])
                if len(g[g.variant == "base" + suffix]) else None,
            }
    (HERE / "error_budget.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
