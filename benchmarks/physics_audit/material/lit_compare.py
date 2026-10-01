#!/usr/bin/env python3
"""Library AA5754-O vs published AA5754-O data: flow curves and Hill48 yield ratios.

Writes lit_curves.csv (flow stress [MPa] vs plastic strain for the library fit
and each published law) and lit_yield_ratios.csv (directional, equibiaxial and
plane-strain yield stress ratios implied by each r-value set under Hill48, and
the measured ones). Sources: README.md, section 'Literature'.
"""
from pathlib import Path

import numpy as np
import pandas as pd

import matlens

HERE = Path(__file__).resolve().parent
EPS = np.array([0.0, 0.002, 0.01, 0.02, 0.05, 0.10, 0.15, 0.19, 0.25, 0.30, 0.40, 0.50, 0.60, 1.0])


def main():
    m = matlens.library_material()
    curves = {
        "library (Hollomon 420/0.30 from 100 MPa, lin+Voce fit)": m.flow_stress(EPS) / 1e6,
        "library Hollomon itself (Swift form)": 420 * (m.hardening_source["eps0"] + EPS) ** 0.30,
        "Iadicola 2008 RD power law K=474 n=0.317 (fit to eps_u 0.19)":
            474 * np.maximum(EPS, 1e-12) ** 0.317,
        "Iadicola 2008 RD Voce S=289 A=0.686 B=-12.3": 289 * (1 - 0.686 * np.exp(-12.3 * EPS)),
        "Coer 2018 RD Voce Y0=102.75 Ysat=292.14 Cy=13.5":
            102.75 + (292.14 - 102.75) * (1 - np.exp(-13.5 * EPS)),
        "Coer 2018 RD Hockett-Sherby Y0=91.74 Ysat=308.63 Cy=7.98 n=0.831":
            308.63 - (308.63 - 91.74) * np.exp(-7.98 * EPS ** 0.831),
    }
    df = pd.DataFrame({"eps_p": EPS, **curves})
    lib = df.iloc[:, 1]
    for c in df.columns[3:]:
        df[f"ratio {c.split(' ')[0]} {c.split(' ')[2] if len(c.split(' ')) > 2 else ''}".strip()] = \
            df[c] / lib
    df.to_csv(HERE / "lit_curves.csv", index=False, float_format="%.4f")
    print(df.iloc[:, :7].to_string(index=False, float_format=lambda x: f"{x:7.1f}"))

    def hill(r0, r45, r90):
        F = r0 / (r90 * (1 + r0)); G = 1 / (1 + r0); H = r0 / (1 + r0)
        N = (r0 + r90) * (1 + 2 * r45) / (2 * r90 * (1 + r0))
        return F, G, H, N

    def ratios(F, G, H, N):
        def s(th):
            c, s_ = np.cos(th), np.sin(th)
            return (F * s_ ** 4 + G * c ** 4 + H * np.cos(2 * th) ** 2 + 2 * N * s_ ** 2 * c ** 2) ** -0.5
        # plane strain along RD (d eps_TD = 0): s22 = H s11 / (F + H)
        k = H / (F + H)
        ps_rd = (F * k ** 2 + G + H * (1 - k) ** 2) ** -0.5
        k2 = H / (G + H)                    # along TD, d eps_RD = 0: s11 = H s22 / (G + H)
        ps_td = (F + G * k2 ** 2 + H * (1 - k2) ** 2) ** -0.5
        return {"s45/s0": s(np.pi / 4), "s90/s0": s(np.pi / 2), "sb/s0": (F + G) ** -0.5,
                "ps_RD/s0 (major stress)": ps_rd, "ps_TD/s0 (major stress)": ps_td}

    sets = {"von Mises": (0.5, 0.5, 0.5, 1.5),
            "library r 0.75/0.70/0.80": hill(0.75, 0.70, 0.80),
            "Iadicola 2008 r 0.69/0.73/0.87": hill(0.69, 0.73, 0.87),
            "Coer 2018 r 0.663/0.860/0.717": hill(0.663, 0.860, 0.717)}
    F, G, H = 0.4349, 0.3915, 0.6085                    # stress-ratio Hill48, sb = 1.10
    sets["stress-ratio Hill48 (s45 .966, s90 .979, sb 1.10)"] = (
        F, G, H, 0.5 * (4 / 0.9660 ** 2 - F - G))
    rows = [{"model": k, **ratios(*v)} for k, v in sets.items()]
    rows.append({"model": "measured, Iadicola 2008 (0.2 % yield; BB at 1 % / 10-15 % work)",
                 "s45/s0": 90.9 / 94.1, "s90/s0": 92.1 / 94.1, "sb/s0": np.nan,
                 "ps_RD/s0 (major stress)": np.nan, "ps_TD/s0 (major stress)": np.nan})
    yr = pd.DataFrame(rows)
    yr.to_csv(HERE / "lit_yield_ratios.csv", index=False, float_format="%.4f")
    print(yr.to_string(index=False, float_format=lambda x: f"{x:6.3f}"))


if __name__ == "__main__":
    main()
