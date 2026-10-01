"""Mesh / element study of the 3D plate (plunge 3 mm + drag 3 mm) in SparLab,
checked against CalculiX on the 2 mm mesh.

Every case is sampled at the top-surface nodes of the 2 mm grid (which every
finer grid contains). For each case and each reference: rms / max of the
released-shape difference overall and per region (core: max(|x|,|y|) < 8 mm,
wall: 8..15, rim outside the clamp window: >= 15) and the springback
(release - form) rms.
usage: python3 plate_mesh.py
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import case_punch as cp
import plate_compare as pc
from common import ROOT

PD = ROOT / "punch"
CASES = ["plate_h2_std_mu0.1_drag3", "plate_h2_std_mu0.1_drag3_tp5", "plate_h2_std_mu0.1_drag3_trav1",
         "plate_h2_im_mu0.1_drag3", "plate_h2_im_mu0.1_drag3_trav1", "plate_h2_nobbar_mu0.1_drag3",
         "plate_h1_std_mu0.1_drag3", "plate_h1_im_mu0.1_drag3", "plate_h1_std_nz4_mu0.1_drag3"]


def sample(case):
    p = json.load(open(PD / case / "params.json"))
    shp, F, pe, ss = pc.sparlab_data(PD / case, p)
    out = {}
    for st, a in shp.items():
        key = {(round(x, 6), round(y, 6)): i for i, (x, y) in enumerate(a[:, :2])}
        out[st] = a
        out[st + "_key"] = key
    return p, out, pc.force_stats(F, p["drag"]), pe


def on_h2(data, st, xy):
    a, key = data[st], data[st + "_key"]
    return np.array([a[key[(round(x, 6), round(y, 6))], 4] for x, y in xy])


def main():
    avail = [c for c in CASES if (PD / c / "sparlab" / "out" / "summary.json").exists()]
    D = {c: sample(c) for c in avail}
    g = cp.grid_mm({"kind": "plate", "h": 2.0, "nz": 2})
    xy = g.xyz[g.xyz[:, 2] > -1e-6][:, :2]
    r = np.max(np.abs(xy), axis=1)
    regions = {"all": r >= 0, "core": r < 8, "wall": (r >= 8) & (r < 15), "rim": r >= 15}
    # CalculiX shapes on the 2 mm grid
    refs = {}
    for name in ("plate_h2_im_mu0.1_drag3", "plate_h2_nobbar_mu0.1_drag3"):
        p = json.load(open(PD / name / "params.json"))
        shp, F, pe = pc.ccx_data(PD / name, p)
        refs["CCX_" + p["elem"]] = {st: shp[st][:, 4] for st in shp}
    fine = [c for c in ("plate_h1_std_nz4_mu0.1_drag3", "plate_h1_im_mu0.1_drag3") if c in D]
    for c in fine:
        refs["SL_" + c] = {st: on_h2(D[c][1], st, xy) for st in ("form", "unload", "release")}
    rows = []
    for c in avail:
        p, data, fs, pe = D[c]
        z = {st: on_h2(data, st, xy) for st in ("form", "unload", "release")}
        row = {"case": c, "h": p["h"], "nz": p["nz"], "elem": p["elem"], "travel": p["travel"], "tp": p.get("tp", 0),
               **fs, "peeq_max": pe, "release_uz_min": float(z["release"].min()),
               "springback_rms": float(np.sqrt(np.mean((z["release"] - z["form"]) ** 2)))}
        for rn, ref in refs.items():
            d = z["release"] - ref["release"]
            for gn, m in regions.items():
                row[f"rel_vs_{rn}_{gn}_rms"] = float(np.sqrt(np.mean(d[m] ** 2)))
            row[f"rel_vs_{rn}_all_max"] = float(np.max(np.abs(d)))
            sb = (z["release"] - z["form"]) - (ref["release"] - ref["form"])
            row[f"sb_vs_{rn}_rms"] = float(np.sqrt(np.mean(sb ** 2)))
        rows.append(row)
    keys = []
    for rr in rows:
        keys += [k for k in rr if k not in keys]
    with open(PD / "plate_mesh_study.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
    for rr in rows:
        print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in rr.items()}))


if __name__ == "__main__":
    main()
