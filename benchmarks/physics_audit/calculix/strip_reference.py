"""An element-independent reference for the plane-strain strip of
case_punch.py: CalculiX C3D20R (quadratic, reduced integration) on fine
meshes, against SparLab at the benchmark's and finer resolutions.

usage: python3 strip_reference.py ccx NAME h=0.5 nz=4 [etype=C3D20R]
       python3 strip_reference.py compare
The strip, punch, material, steps and supports are those of case_punch.py
(kind=strip, mu=0.1, depth=4). Compared: punch force per width along the
plunge, and the top-surface z displacement at the end of form, unload and
release, as a function of x (interpolated onto SparLab's nodes).
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import case_punch as cp
from common import E_PA, ROOT, ccx_material_cards, last_per_step, parse_dat, run_ccx
from punch_geometry import cylinder_ring

OUT = ROOT / "strip_reference"


def quad_grid(lx, w, lz, nx, nz):
    """C3D20 grid on [-lx/2, lx/2] x [0, w] x [-lz, 0], one element in y."""
    xs = np.linspace(-lx / 2, lx / 2, 2 * nx + 1)
    ys = np.linspace(0, w, 3)
    zs = np.linspace(-lz, 0, 2 * nz + 1)
    ids = -np.ones((len(xs), 3, len(zs)), int)
    nodes = []

    def nid(i, j, k):
        if ids[i, j, k] < 0:
            ids[i, j, k] = len(nodes)
            nodes.append([xs[i], ys[j], zs[k]])
        return ids[i, j, k]

    conn = []
    for k in range(nz):
        for i in range(nx):
            I, K = 2 * i, 2 * k
            c = [(I, 0, K), (I + 2, 0, K), (I + 2, 2, K), (I, 2, K),
                 (I, 0, K + 2), (I + 2, 0, K + 2), (I + 2, 2, K + 2), (I, 2, K + 2)]
            m = [(I + 1, 0, K), (I + 2, 1, K), (I + 1, 2, K), (I, 1, K),
                 (I + 1, 0, K + 2), (I + 2, 1, K + 2), (I + 1, 2, K + 2), (I, 1, K + 2),
                 (I, 0, K + 1), (I + 2, 0, K + 1), (I + 2, 2, K + 1), (I, 2, K + 1)]
            conn.append([nid(*q) for q in c + m])
    return np.array(nodes), np.array(conn)


def ccx_ref_deck(h, nz, etype, ctype="NODE TO SURFACE", inc=0.1, unload="remove"):
    nx = int(round(40.0 / h))
    W = h
    x, conn = quad_grid(40.0, W, 1.0, nx, nz)
    nodes_p, conn_p, faces_p = cylinder_ring((0.0, cp.R_TOOL), cp.R_TOOL, -1.0, W + 1.0)
    n0 = len(x)
    L = ["*HEADING", "strip reference", "*NODE, NSET=NSHEET"]
    L += [f"{i+1}, {q[0]:.12g}, {q[1]:.12g}, {q[2]:.12g}" for i, q in enumerate(x)]
    L.append("*NODE, NSET=PUNCH")
    L += [f"{n0+i+1}, {q[0]:.12g}, {q[1]:.12g}, {q[2]:.12g}" for i, q in enumerate(nodes_p)]
    L.append(f"*ELEMENT, TYPE={etype}, ELSET=SHEET")
    for e, c in enumerate(conn):
        ids = [str(n + 1) for n in c]
        L.append(f"{e+1}, " + ", ".join(ids[:15]) + ",")
        L.append(", ".join(ids[15:]))
    ne = len(conn)
    L.append("*ELEMENT, TYPE=C3D8, ELSET=TOOL")
    for e, c in enumerate(conn_p):
        L.append(f"{ne+e+1}, " + ", ".join(str(n0 + n + 1) for n in c))

    def nset(name, mask):
        ids = [str(i + 1) for i in np.nonzero(mask)[0]]
        L.append(f"*NSET, NSET={name}")
        for k in range(0, len(ids), 12):
            L.append(", ".join(ids[k:k + 12]))
    tol = 1e-6
    nset("CLAMP", np.abs(x[:, 0]) >= cp.CLAMP - tol)
    nset("SUPA", (np.abs(x[:, 0] + cp.SUP) < tol) & (x[:, 2] > -tol))
    nset("SUPB", (np.abs(x[:, 0] - cp.SUP) < tol) & (x[:, 2] > -tol))
    nset("NTOP", (x[:, 2] > -tol) & (x[:, 1] < tol))
    L.append("*ELSET, ELSET=ETOPL")
    top = [e + 1 for e, c in enumerate(conn) if np.all(x[c[4:8], 2] > -tol)]
    for k in range(0, len(top), 12):
        L.append(", ".join(str(v) for v in top[k:k + 12]))
    L += ["*SURFACE, NAME=SSLAVE, TYPE=ELEMENT", "ETOPL, S2", "*SURFACE, NAME=SMASTER, TYPE=ELEMENT"]
    L += [f"{ne+e+1}, S{f}" for e, f in faces_p]
    L.append(ccx_material_cards("AL"))
    L.append("*MATERIAL, NAME=STEEL\n*ELASTIC\n210000., 0.3")
    L.append("*SOLID SECTION, ELSET=SHEET, MATERIAL=AL")
    L.append("*SOLID SECTION, ELSET=TOOL, MATERIAL=STEEL")
    kappa = cp.PENALTY_S * E_PA / 1e6 / h
    L += ["*SURFACE INTERACTION, NAME=SI", "*SURFACE BEHAVIOR, PRESSURE-OVERCLOSURE=LINEAR",
          f"{kappa:.6g}, 0.01", "*FRICTION", f"0.1, {kappa:.6g}",
          f"*CONTACT PAIR, INTERACTION=SI, TYPE={ctype}", "SSLAVE, SMASTER"]
    outs = ["*NODE PRINT, NSET=PUNCH, TOTALS=ONLY", "RF", "*NODE PRINT, NSET=NTOP", "U"]
    clamp = ["*BOUNDARY", "NSHEET, 2, 2, 0.0", "CLAMP, 1, 1, 0.0", "CLAMP, 3, 3, 0.0"]
    D = 4.0
    f = inc / D
    L += ["*STEP, NLGEOM, INC=100000", "*STATIC", f"{f:.6g}, 1.0, {f*1e-4:.6g}, {f:.6g}"] + clamp + \
         ["*BOUNDARY", "PUNCH, 1, 2, 0.0", f"PUNCH, 3, 3, {-D}"] + outs + ["*END STEP"]
    # unload: the punch's contact pair is removed (as SparLab's 'release'
    # unload step removes the tool); retracting the punch diverged with
    # C3D20R slaves at the first loss of contact (increment < minimum).
    if unload == "remove":
        L += ["*STEP, NLGEOM, INC=100000", "*STATIC", "0.05, 1.0, 1e-7, 0.05",
              "*MODEL CHANGE, TYPE=CONTACT PAIR, REMOVE", "SSLAVE, SMASTER"] + clamp + \
             ["*BOUNDARY", "PUNCH, 1, 2, 0.0", f"PUNCH, 3, 3, {-D}"] + outs + ["*END STEP"]
    else:  # retract the punch 1 mm above the sheet; contact is removed in the release step
        f = inc / (D + 1)
        L += ["*STEP, NLGEOM, INC=100000", "*STATIC", f"{f:.6g}, 1.0, {f*1e-4:.6g}, {f:.6g}"] + clamp + \
             ["*BOUNDARY", "PUNCH, 1, 2, 0.0", "PUNCH, 3, 3, 1.0"] + outs + ["*END STEP"]
    L += ["*STEP, NLGEOM, INC=100000", "*STATIC", "0.05, 1.0, 1e-7, 0.05", "*BOUNDARY, OP=NEW",
          "NSHEET, 2, 2, 0.0", "PUNCH, 1, 2, 0.0", f"PUNCH, 3, 3, {-D if unload == 'remove' else 1.0}", "SUPA, 1, 1, 0.0", "SUPA, 3, 3, 0.0",
          "SUPB, 3, 3, 0.0"] + outs + ["*END STEP"]
    return "\n".join(L) + "\n", W


def ccx_profiles(cdir, W):
    blocks = parse_dat(cdir / "job.dat")
    inp = (cdir / "job.inp").read_text().splitlines()
    xyz = {}
    i0 = inp.index("*NODE, NSET=NSHEET")
    for line in inp[i0 + 1:]:
        if line.startswith("*"):
            break
        a = line.split(",")
        xyz[int(a[0])] = [float(v) for v in a[1:]]
    U = last_per_step(blocks, "displacements")
    prof = {}
    for st, k in (("form", 1), ("unload", 2), ("release", 3)):
        rows = sorted((xyz[int(r[0])][0], r[3]) for r in U[k][3])
        prof[st] = np.array(rows)
    force = []
    for b in blocks:
        if b[0].startswith("total force") and b[4] == 1:
            force.append((4.0 * b[2], -b[3][0][2] / W))
    return prof, np.array(force)


def sparlab_profiles(case):
    prof = {}
    for st in ("form", "unload", "release"):
        rows = np.loadtxt(case / f"shape_{st}.csv", delimiter=",", skiprows=1)
        rows = rows[rows[:, 1] < 1e-9]
        prof[st] = np.column_stack([rows[:, 0], rows[:, 4]])
    f = np.array([[float(r["depth_mm"]), float(r["fz"])] for r in csv.DictReader(open(case / "force.csv"))
                  if r["code"] == "sparlab" and float(r["x_mm"]) < 1e-6])
    return prof, f


def compare(ref_name):
    rdir = OUT / ref_name
    W = json.load(open(rdir / "params.json"))["h"]
    rprof, rforce = ccx_profiles(rdir / "ccx", W)
    rows = []
    for case in sorted((ROOT / "punch").glob("strip_*")):
        if not (case / "result.json").exists():
            continue
        prof, f = sparlab_profiles(case)
        p = json.load(open(case / "params.json"))
        r = {"case": case.name, "h": p["h"], "nz": p["nz"], "elem": p["elem"], "mu": p["mu"],
             "travel": p["travel"], "tp": p.get("tp", 0)}
        for st in ("form", "unload", "release"):
            ref = np.interp(prof[st][:, 0], rprof[st][:, 0], rprof[st][:, 1])
            d = prof[st][:, 1] - ref
            m = np.abs(prof[st][:, 0]) <= 19.0
            r[f"{st}_dz_rms_mm"] = float(np.sqrt(np.mean(d[m] ** 2)))
            r[f"{st}_dz_max_mm"] = float(np.max(np.abs(d[m])))
            r[f"{st}_uz_min_mm"] = float(prof[st][:, 1].min())
        r["ref_release_uz_min_mm"] = float(rprof["release"][:, 1].min())
        fd = np.interp(4.0, f[:, 0], f[:, 1])
        fr = np.interp(4.0, rforce[:, 0], rforce[:, 1])
        r["F4_N_per_mm"] = float(fd)
        r["F4_ref"] = float(fr)
        r["F4_rel"] = float((fd - fr) / fr)
        rows.append(r)
    with open(OUT / f"compare_{ref_name}.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    for r in rows:
        print(f"{r['case']:34s} F4 {r['F4_N_per_mm']:7.2f} ({r['F4_rel']:+.3f})  form {r['form_dz_rms_mm']:.4f}/{r['form_dz_max_mm']:.4f}"
              f"  unload {r['unload_dz_rms_mm']:.4f}/{r['unload_dz_max_mm']:.4f}  release {r['release_dz_rms_mm']:.4f}/{r['release_dz_max_mm']:.4f}"
              f"  rel-min {r['release_uz_min_mm']:.4f} (ref {r['ref_release_uz_min_mm']:.4f})")


def main():
    if sys.argv[1] == "ccx":
        name = sys.argv[2]
        kv = dict(a.split("=") for a in sys.argv[3:])
        h, nz, et = float(kv.get("h", 0.5)), int(kv.get("nz", 4)), kv.get("etype", "C3D20R")
        d = OUT / name
        ct = {"n2s": "NODE TO SURFACE", "s2s": "SURFACE TO SURFACE"}[kv.get("contact", "n2s")]
        inp, W = ccx_ref_deck(h, nz, et, ct, float(kv.get("inc", 0.1)), kv.get("unload", "remove"))
        d.mkdir(parents=True, exist_ok=True)
        (d / "params.json").write_text(json.dumps({"h": h, "nz": nz, "etype": et, "contact": ct, "argv": sys.argv[3:]}))
        rc, wall = run_ccx(inp, d / "ccx")
        (d / "ccx" / "wall_s.txt").write_text(f"{wall:.1f}\n")
        print(name, rc, f"{wall:.1f}s")
    else:
        compare(sys.argv[2])


if __name__ == "__main__":
    main()
