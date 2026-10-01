"""Shared helpers of the CalculiX cross-check lens.

SparLab decks are written in SI (m, Pa); CalculiX decks in mm, N, MPa.
Both codes get the same nodes (structured hex grids built here, identical
coordinates), the same material law and, where possible, the same contact law.
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[2]
SPARLAB_FORM = REPO / "build" / "bin" / "sparlab_form"
CCX = "/usr/bin/ccx"

# AA5754-O as the benchmark deck carries it (precomp library fit of
# Hollomon K = 420 MPa, n = 0.30 by linear + Voce isotropic hardening).
E_PA = 70.0e9
NU = 0.33
SY0 = 100.0e6
H_LIN = 237610605.32948908
Q_VOCE = 123864572.88204344
B_VOCE = 14.377555530345461


def flow_stress(a):
    """Yield stress [Pa] at equivalent plastic strain a."""
    a = np.asarray(a, dtype=float)
    return SY0 + H_LIN * a + Q_VOCE * (1.0 - np.exp(-B_VOCE * a))


def sparlab_material(hill=False):
    plast = {
        "yield_stress": SY0,
        "hardening_modulus": H_LIN,
        "saturation_stress": Q_VOCE,
        "saturation_rate": B_VOCE,
        "kinematic_hardening_modulus": 0.0,
    }
    if hill:
        plast["yield_criterion"] = "hill48"
        plast["anisotropy"] = {"r0": 0.75, "r45": 0.7, "r90": 0.8,
                               "rolling_direction": [1.0, 0.0, 0.0],
                               "sheet_normal": [0.0, 0.0, 1.0]}
    return {"name": "AA5754-O", "youngs_modulus": E_PA, "poisson_ratio": NU,
            "density": 2670.0, "plasticity": plast}


def ccx_material_cards(name="AL", eps_max=3.0):
    """*MATERIAL with the same hardening as an *PLASTIC table (MPa, -).

    CalculiX interpolates the table linearly and holds the last value
    beyond it, so the table is dense where the Voce term curves and runs
    well past any strain the tests reach."""
    a = np.unique(np.concatenate([np.linspace(0, 0.3, 301), np.linspace(0.3, eps_max, 271)]))
    lines = [f"*MATERIAL, NAME={name}", "*ELASTIC", f"{E_PA/1e6:.6g}, {NU}", "*PLASTIC"]
    for ai, si in zip(a, flow_stress(a) / 1e6):
        lines.append(f"{si:.9g}, {ai:.9g}")
    return "\n".join(lines) + "\n"


def table_interp_error():
    """Largest relative error of the linear interpolation of the table."""
    a = np.unique(np.concatenate([np.linspace(0, 0.3, 301), np.linspace(0.3, 3.0, 271)]))
    x = np.linspace(0, 3.0, 300001)
    return float(np.max(np.abs(np.interp(x, a, flow_stress(a)) - flow_stress(x)) / flow_stress(x)))


# ----------------------------------------------------------------------------
# structured hex grid, node numbering shared by both decks (1-based for ccx)

class Grid:
    """Structured nx x ny x nz hex grid on [x0, x0+lx] x ... (any unit)."""

    def __init__(self, x0, y0, z0, lx, ly, lz, nx, ny, nz):
        self.nx, self.ny, self.nz = nx, ny, nz
        xs = np.linspace(x0, x0 + lx, nx + 1)
        ys = np.linspace(y0, y0 + ly, ny + 1)
        zs = np.linspace(z0, z0 + lz, nz + 1)
        Z, Y, X = np.meshgrid(zs, ys, xs, indexing="ij")
        self.xyz = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
        self.idx = lambda i, j, k: k * (ny + 1) * (nx + 1) + j * (nx + 1) + i
        conn = []
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    n = self.idx
                    conn.append([n(i, j, k), n(i + 1, j, k), n(i + 1, j + 1, k), n(i, j + 1, k),
                                 n(i, j, k + 1), n(i + 1, j, k + 1), n(i + 1, j + 1, k + 1),
                                 n(i, j + 1, k + 1)])
        self.conn = np.array(conn)


# ----------------------------------------------------------------------------
# running

def run(cmd, cwd, log, timeout=4 * 3600):
    env = dict(os.environ)
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["OMP_NUM_THREADS"] = "1"
    t0 = time.time()
    with open(log, "w") as fh:
        p = subprocess.run(cmd, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT, env=env,
                           timeout=timeout)
    return p.returncode, time.time() - t0


def run_sparlab(deck: dict, case_dir: Path):
    case_dir.mkdir(parents=True, exist_ok=True)
    (case_dir / "deck.json").write_text(json.dumps(deck, indent=1))
    rc, wall = run([str(SPARLAB_FORM), "--config", "deck.json", "--output", "out", "--threads", "1",
                    "--strict-config", "--no-vtk"], case_dir, case_dir / "sparlab.log")
    return rc, wall


def run_ccx(inp_text: str, case_dir: Path, job="job"):
    case_dir.mkdir(parents=True, exist_ok=True)
    (case_dir / f"{job}.inp").write_text(inp_text)
    rc, wall = run([CCX, "-i", job], case_dir, case_dir / "ccx.log")
    return rc, wall


def read_sparlab_nodes(path):
    a = np.loadtxt(path, delimiter=",", skiprows=1)
    return a  # node, X, Y, Z, ux, uy, uz


def read_sparlab_elements(path):
    return np.loadtxt(path, delimiter=",", skiprows=1)


def ccx_sta_complete(case_dir: Path, job="job"):
    log = (case_dir / "ccx.log").read_text(errors="replace")
    return ("Job finished" in log) and ("*ERROR" not in log)


def parse_dat(path):
    """Blocks of a CalculiX .dat file: list of (title, set, time, rows, step)."""
    blocks = []
    cur = None
    step = 0
    for line in Path(path).read_text(errors="replace").splitlines():
        ms = re.match(r"\s*S T E P\s+(\d+)", line)
        if ms:
            step = int(ms.group(1))
            continue
        m = re.match(r"\s*(.+?)\s*for set\s+(\S+)\s+and time\s+(\S+)", line)
        if m:
            cur = (m.group(1).strip(), m.group(2), float(m.group(3)), [], step)
            blocks.append(cur)
            continue
        if cur is None:
            continue
        parts = line.split()
        if not parts:
            continue
        try:
            cur[3].append([float(p) for p in parts])
        except ValueError:
            pass
    return blocks


def last_per_step(blocks, prefix):
    """The last block (end of step) whose title starts with prefix, by step."""
    out = {}
    for b in blocks:
        if b[0].startswith(prefix):
            out[b[4]] = b
    return out
