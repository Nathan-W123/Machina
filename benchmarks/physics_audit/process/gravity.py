"""Self-weight sag of a flat sheet (upper bound for flat flanges; a formed shell
is stiffer): CalculiX 2.21 (/usr/bin/ccx) S4 shells, linear static, gravity
normal to the sheet, on the benchmark's 3-2-1 support pattern (three corner
points of the former clamp: xyz, yz, z) and clamped on all edges (the frame),
plus the Timoshenko closed form for the clamped square plate
(w = 0.00126 q a^4 / D, nu = 0.3) as a check. Cases: the benchmark blank (40 mm,
1 mm, window 30 mm) and Machina-scale panels (0.6-1.5 m, 1-3.2 mm; Machina
quotes 0.080" and 0.125" sheet). Writes gravity.json."""
import json
import subprocess
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent / "gravity"
E, NU, RHO, G = 70e9, 0.33, 2670.0, 9.81


def run(name, L, t, support, n=40, a_frac=None):
    h = L / n
    nodes = {}
    k = 1
    for j in range(n + 1):
        for i in range(n + 1):
            nodes[(i, j)] = k
            k += 1
    lines = ["*NODE, NSET=NALL"]
    for (i, j), idx in nodes.items():
        lines.append(f"{idx}, {i*h - L/2:.9g}, {j*h - L/2:.9g}, 0")
    lines.append("*ELEMENT, TYPE=S4, ELSET=EALL")
    e = 1
    for j in range(n):
        for i in range(n):
            lines.append(f"{e}, {nodes[(i,j)]}, {nodes[(i+1,j)]}, {nodes[(i+1,j+1)]}, {nodes[(i,j+1)]}")
            e += 1
    lines += ["*MATERIAL, NAME=AL", "*ELASTIC", f"{E}, {NU}", "*DENSITY", f"{RHO}",
              "*SHELL SECTION, ELSET=EALL, MATERIAL=AL", f"{t}"]
    lines.append("*BOUNDARY")
    if support == "321":
        m = a_frac if a_frac is not None else 0
        i0, i1 = m, n - m
        lines += [f"{nodes[(i0,i0)]}, 1, 3", f"{nodes[(i1,i0)]}, 2, 3", f"{nodes[(i0,i1)]}, 3, 3"]
    else:  # clamped edges
        for (i, j), idx in nodes.items():
            if i in (0, n) or j in (0, n):
                lines.append(f"{idx}, 1, 6")
    lines += ["*STEP", "*STATIC", "*DLOAD", f"EALL, GRAV, {G}, 0., 0., -1.",
              "*NODE FILE", "U", "*NODE PRINT, NSET=NALL", "U", "*END STEP"]
    d = HERE / name
    d.mkdir(parents=True, exist_ok=True)
    (d / "job.inp").write_text("\n".join(lines) + "\n")
    subprocess.run(["ccx", "-i", "job"], cwd=d, check=True, capture_output=True,
                   env={"OMP_NUM_THREADS": "1", "PATH": "/usr/bin:/bin"})
    txt = (d / "job.dat").read_text().split("\n")
    w = []
    for l in txt:
        p = l.split()
        if len(p) == 4:
            try:
                w.append(float(p[3]))
            except ValueError:
                pass
    w = np.array(w)
    return float(np.abs(w).max())


out = []
cases = [("bench_40mm_t1", 0.040, 0.001), ("panel_0.6m_t1", 0.6, 0.001),
         ("panel_1.2m_t2", 1.2, 0.002), ("panel_1.2m_t3.2", 1.2, 0.0032),
         ("panel_1.5m_t3.2", 1.5, 0.0032)]
for name, L, t in cases:
    D = E * t ** 3 / (12 * (1 - NU ** 2))
    q = RHO * G * t
    win = L - 2 * (0.005 if L < 0.1 else 0.05)        # frame window: 5 mm / 50 mm clamp
    rec = {"case": name, "L_m": L, "t_m": t, "q_Pa": q,
           "w_321_corners_mm": run(name + "_321", L, t, "321", n=40, a_frac=2 if L < 0.1 else 2) * 1e3,
           "w_clamped_window_fe_mm": run(name + "_clamped", win, t, "clamped") * 1e3,
           "w_clamped_window_timoshenko_mm": 0.00126 * q * win ** 4 / (E * t ** 3 / (12 * 0.91)) * 1e3}
    out.append(rec)
    print(rec)
(HERE / "gravity.json").write_text(json.dumps(out, indent=1))
