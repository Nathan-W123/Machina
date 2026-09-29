"""A test double of sparlab_form that writes a contract-conforming result.

It is *not* a solver. It reads the deck (mesh, forming steps, trajectory) as
the real executable would, and - a test-only shortcut - the provenance file
`commanded.npz` that `precomp.fea.build_deck` writes beside the deck. Every
node of a through-thickness column moves vertically by the commanded height
at its (X, Y), times a springback factor S after the "form" step:

    form step:           u_z = z_cmd(X, Y)
    unload / release:    u_z = S z_cmd(X, Y),   S = $FAKE_SPARLAB_SPRINGBACK (0.9)

so the formed tool-side surface is exactly S times the commanded one at the
nodes, a linear "springback" whose displacement-adjustment inverse is known.
Environment switches for the tests: FAKE_SPARLAB_COUNTER (append one line per
run to this file), FAKE_SPARLAB_FAIL (exit 3 with a message).
"""

import argparse
import json
import os
import sys

import numpy as np
from scipy.interpolate import RegularGridInterpolator


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if "--version" in argv:
        print("sparlab_form 0.0.0-fake (test double)")
        return 0
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args(argv)
    counter = os.environ.get("FAKE_SPARLAB_COUNTER")
    if counter:
        with open(counter, "a", encoding="utf-8") as handle:
            handle.write(args.config + "\n")
    if os.environ.get("FAKE_SPARLAB_FAIL"):
        print("fake sparlab_form: failure requested by FAKE_SPARLAB_FAIL")
        return 3
    with open(args.config, encoding="utf-8") as handle:
        deck = json.load(handle)
    deck_dir = os.path.dirname(os.path.abspath(args.config))
    out = args.output
    os.makedirs(out, exist_ok=True)
    m = deck["mesh"]
    nx, ny, nz = m["nx"], m["ny"], m["nz"]
    xs = m["x0"] + m["lx"] * np.arange(nx + 1) / nx
    ys = m["y0"] + m["ly"] * np.arange(ny + 1) / ny
    zs = m["z0"] + m["lz"] * np.arange(nz + 1) / nz
    Z, Y, X = np.meshgrid(zs, ys, xs, indexing="ij")          # node(i, j, k): x fastest
    X, Y, Z = X.ravel(), Y.ravel(), Z.ravel()
    with np.load(os.path.join(deck_dir, "commanded.npz")) as data:
        grid = json.loads(str(data["grid"]))
        zc = data["z"]
    gx = grid["x0"] + grid["h"] * np.arange(grid["nx"])
    gy = grid["y0"] + grid["h"] * np.arange(grid["ny"])
    interp = RegularGridInterpolator((gy, gx), zc, bounds_error=False, fill_value=0.0)
    zcmd = interp(np.column_stack([Y, X]))
    s = float(os.environ.get("FAKE_SPARLAB_SPRINGBACK", "0.9"))
    steps = deck["forming"]["steps"]
    summary_steps = []
    for k, step in enumerate(steps):
        factor = 1.0 if step["type"] == "form" else s
        uz = factor * zcmd
        table = np.column_stack([np.arange(len(X)), X, Y, Z, np.zeros_like(X),
                                 np.zeros_like(X), uz])
        name = step["name"]
        np.savetxt(os.path.join(out, f"step_{k}_{name}_nodes.csv"), table, delimiter=",",
                   header="node,X,Y,Z,ux,uy,uz", comments="",
                   fmt=["%d"] + ["%.17g"] * 6)
        ne = nx * ny * nz
        et = np.column_stack([np.arange(ne), np.full(ne, 0.1 * factor), np.full(ne, 1.0e8)])
        np.savetxt(os.path.join(out, f"step_{k}_{name}_elements.csv"), et, delimiter=",",
                   header="element,eq_plastic_strain,von_mises", comments="",
                   fmt=["%d", "%.17g", "%.17g"])
        summary_steps.append({"name": name, "type": step["type"], "completed": True,
                              "increments": 1, "iterations": 1})
    traj = np.loadtxt(os.path.join(deck_dir, deck["forming"]["tools"][0]["trajectory"]["file"]),
                      delimiter=",", skiprows=1)
    rows = []
    for i, (t, cx, cy, cz) in enumerate(traj[::5]):
        fz = -200.0 - 1.0e5 * max(0.0, -cz)
        rows.append(f"0,{i},{t:.17g},tool,{cx:.17g},{cy:.17g},{cz:.17g},0,0,{fz:.17g},3")
    with open(os.path.join(out, "tool_forces.csv"), "w", encoding="utf-8") as handle:
        handle.write("step,increment,t,tool,cx,cy,cz,fx,fy,fz,active_nodes\n")
        handle.write("\n".join(rows) + "\n")
    with open(os.path.join(out, "summary.json"), "w", encoding="utf-8") as handle:
        json.dump({"completed": True, "steps": summary_steps, "runtime_s": 0.0,
                   "warnings": [], "solver": "fake"}, handle)
    return 0


if __name__ == "__main__":
    sys.exit(main())
