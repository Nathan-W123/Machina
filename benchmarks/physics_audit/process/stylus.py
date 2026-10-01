"""Forming-stylus bending (cantilever, steel E = 210 GPa): delta = F L^3 / (3 E I),
I = pi d^4 / 64, under the simulated lateral tool force (8 test parts, cached
runs: mean 136 N, peak 321 N; measure_states.csv) and the 2-20 kN peak forces
Machina's patent quotes (lateral share taken as 1/4 of the peak, the ratio in
our runs: 321/924 N ~ 0.35 at peak). Writes stylus.json."""
import json
from math import pi
from pathlib import Path

E = 210e9
out = []
for d_mm in (10, 16, 25):
    for L_mm in (60, 100, 150):
        I = pi * (d_mm * 1e-3) ** 4 / 64
        c = (L_mm * 1e-3) ** 3 / (3 * E * I)      # m/N
        out.append({"d_mm": d_mm, "L_mm": L_mm, "compliance_um_per_N": c * 1e6,
                    "delta_mm_at_136N": c * 136 * 1e3, "delta_mm_at_321N": c * 321 * 1e3,
                    "delta_mm_at_1kN": c * 1000 * 1e3})
        print(out[-1])
Path(__file__).with_name("stylus.json").write_text(json.dumps(out, indent=1))
