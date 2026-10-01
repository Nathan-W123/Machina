"""Compare SparLab plate runs (any element / tool travel) against one
CalculiX plate run (default: C3D8I, plate_h2_im_mu0.1_drag).

usage: python3 plate_compare.py [REF_CASE]
Writes punch/plate_compare_<REF_CASE>.csv: per SparLab case and step, the
rms / max |dz| of the top surface against CalculiX, the deepest point, and
the punch forces (end of plunge, mean over the second half of the drag).
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import case_punch as cp
from common import ROOT, last_per_step, parse_dat


def ccx_data(case, p):
    blocks = parse_dat(case / "ccx" / "job.dat")
    g = cp.grid_mm(p)
    top = np.nonzero(g.xyz[:, 2] > -1e-6)[0]
    U = last_per_step(blocks, "displacements")
    nseg = 2 if p["drag"] > 0 else 1
    steps = {"form": nseg, "unload": nseg + 1, "release": nseg + 2}
    shp = {}
    for st, k in steps.items():
        cu = {int(r[0]): r[1:4] for r in U[k][3]}
        shp[st] = np.array([[g.xyz[i, 0], g.xyz[i, 1], *cu[i + 1]] for i in top])
    _, pts = cp.tool_path_mm(p)
    F = []
    for b in blocks:
        if b[0].startswith("total force") and b[4] <= nseg:
            a, bb = pts[b[4] - 1], pts[b[4]]
            c = a + (b[2] - (b[4] - 1)) * (bb - a)  # printed time is the total time
            rf = np.array(b[3][0])
            F.append((b[4], cp.R_TOOL - c[2], c[0], -rf[2], -rf[0]))
    peeq = last_per_step(blocks, "equivalent plastic")
    pe = float(np.max(np.array(peeq[nseg][3])[:, 2])) if nseg in peeq else float("nan")
    return shp, np.array(F), pe


def sparlab_data(case, p):
    g = cp.grid_mm(p)
    top = np.nonzero(g.xyz[:, 2] > -1e-6)[0]
    shp = {}
    for k, st in ((1, "form"), (2, "unload"), (3, "release")):
        sn = np.loadtxt(case / "sparlab" / "out" / f"step_{k}_{st}_nodes.csv", delimiter=",", skiprows=1)
        key = {tuple(np.round(q * 1e3, 6)): j for j, q in enumerate(sn[:, 1:4])}
        shp[st] = np.array([[g.xyz[i, 0], g.xyz[i, 1], *(sn[key[tuple(np.round(g.xyz[i], 6))], 4:7] * 1e3)]
                            for i in top])
    tf = [r for r in csv.DictReader(open(case / "sparlab" / "out" / "tool_forces.csv")) if r["step"] == "1"]
    F = np.array([(0, cp.R_TOOL - float(r["cz"]) * 1e3, float(r["cx"]) * 1e3, float(r["fz"]), float(r["fx"]))
                  for r in tf])
    el = np.loadtxt(case / "sparlab" / "out" / "step_1_form_elements.csv", delimiter=",", skiprows=1)
    ss = json.load(open(case / "sparlab" / "out" / "summary.json"))
    return shp, F, float(el[:, 1].max()), ss


def force_stats(F, drag):
    plunge = F[F[:, 2] < 1e-6]
    out = {"F_plunge_end_N": float(plunge[np.argmax(plunge[:, 1]), 3])}
    d = F[F[:, 2] > drag / 2]
    out["drag_fz_mean_N"] = float(d[:, 3].mean())
    out["drag_fx_mean_N"] = float(d[:, 4].mean())
    return out


def main():
    ref = sys.argv[1] if len(sys.argv) > 1 else "plate_h2_im_mu0.1_drag"
    rcase = ROOT / "punch" / ref
    p = json.load(open(rcase / "params.json"))
    cshp, cF, cpe = ccx_data(rcase, p)
    rows = []
    cf = force_stats(cF, p["drag"])
    rows.append({"case": f"CCX:{ref}", **{k: v for k, v in cf.items()}, "peeq_max": cpe,
                 **{f"{st}_uz_min_mm": float(cshp[st][:, 4].min()) for st in cshp}})
    for case in sorted((ROOT / "punch").glob("plate_*")):
        if not (case / "sparlab" / "out" / "summary.json").exists():
            continue
        q = json.load(open(case / "params.json"))
        if (q["kind"], q["h"], q["nz"], q["depth"], q["drag"]) != (p["kind"], p["h"], p["nz"], p["depth"], p["drag"]):
            continue
        sshp, sF, spe, ss = sparlab_data(case, q)
        r = {"case": case.name, "elem": q["elem"], "travel_mm": q["travel"], "tp": q.get("tp", 0)}
        r.update(force_stats(sF, q["drag"]))
        r["peeq_max"] = spe
        for st in ("form", "unload", "release"):
            dz = sshp[st][:, 4] - cshp[st][:, 4]
            r[f"{st}_uz_min_mm"] = float(sshp[st][:, 4].min())
            r[f"{st}_dz_rms_mm"] = float(np.sqrt(np.mean(dz ** 2)))
            r[f"{st}_dz_max_mm"] = float(np.max(np.abs(dz)))
        # springback = released - formed (z), compared as a field
        sb_s = sshp["release"][:, 4] - sshp["form"][:, 4]
        sb_c = cshp["release"][:, 4] - cshp["form"][:, 4]
        r["springback_rms_sparlab_mm"] = float(np.sqrt(np.mean(sb_s ** 2)))
        r["springback_rms_ccx_mm"] = float(np.sqrt(np.mean(sb_c ** 2)))
        r["springback_diff_rms_mm"] = float(np.sqrt(np.mean((sb_s - sb_c) ** 2)))
        r["sparlab_runtime_s"] = ss["runtime_s"]
        rows.append(r)
        np.savetxt(case / f"shape_release_vs_{ref}.csv",
                   np.column_stack([sshp["release"][:, :2], sshp["release"][:, 4], cshp["release"][:, 4]]),
                   delimiter=",", header="X_mm,Y_mm,uz_sparlab,uz_ccx", comments="")
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    out = ROOT / "punch" / f"plate_compare_{ref}.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader(); w.writerows(rows)
    for r in rows:
        print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}))


if __name__ == "__main__":
    main()
