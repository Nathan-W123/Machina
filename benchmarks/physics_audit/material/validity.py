#!/usr/bin/env python3
"""Validity limits of the material model in the benchmark runs (lens part c).

For a finished run (runs/<part>/<variant>) this measures, from the solver's
own output:
  * the equivalent plastic strain per element (max, percentiles, share of the
    part beyond the uniaxial test range eps_p > 0.19-0.22, where every
    hardening law is an extrapolation), in the rim band and deeper;
  * the membrane strain state of the mid-surface (principal log strains of
    each cell of the mid-surface node grid; strain ratio beta = eps2/eps1:
    0 plane strain, 1 equibiaxial) - which part of the yield locus the wall
    samples;
  * the through-thickness (transverse) shear: the angle between the deformed
    through-thickness node column (bottom -> top) and the normal of the
    deformed mid-surface, per column and per layer;
  * thickness strain ln(t/t0) from the top-bottom distance.
Writes validity_<part>_<variant>.json and prints a summary.
Usage: validity.py PART VARIANT [STEP]   (STEP default "release")
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def load(run: Path, step: str):
    out = run / "output"
    stem = {"form": "step_1_form", "unload": "step_2_unload", "release": "step_3_release"}[step]
    nodes = np.loadtxt(out / f"{stem}_nodes.csv", delimiter=",", skiprows=1)
    el = np.loadtxt(out / f"{stem}_elements.csv", delimiter=",", skiprows=1)
    mesh = json.loads((out / "mesh.json").read_text())
    return nodes, el, mesh


def grid_nodes(nodes):
    X = nodes[:, 1:4]
    u = nodes[:, 4:7]
    xs = np.unique(np.round(X[:, 0], 9))
    ys = np.unique(np.round(X[:, 1], 9))
    zs = np.unique(np.round(X[:, 2], 9))
    ix = np.searchsorted(xs, np.round(X[:, 0], 9))
    iy = np.searchsorted(ys, np.round(X[:, 1], 9))
    iz = np.searchsorted(zs, np.round(X[:, 2], 9))
    G = np.zeros((len(zs), len(ys), len(xs), 3))
    G0 = np.zeros_like(G)
    G[iz, iy, ix] = X + u
    G0[iz, iy, ix] = X
    return G0, G, xs, ys, zs


def main():
    part, variant = sys.argv[1], sys.argv[2]
    step = sys.argv[3] if len(sys.argv) > 3 else "release"
    run = HERE / "runs" / part / variant
    nodes, el, mesh = load(run, step)
    G0, G, xs, ys, zs = grid_nodes(nodes)
    nz = len(zs)
    t0 = zs[-1] - zs[0]
    # commanded (target) depth at columns, for the rim band / deeper split
    tgt = np.load(run / "surfaces.npz")
    import matlens
    _, target = matlens.part_inputs(part)
    zt = target.interpolate(G0[-1, :, :, 0].ravel(), G0[-1, :, :, 1].ravel(),
                            masked=False).reshape(G0.shape[1:3])
    in_part = zt < -1e-6
    upper = in_part & (zt > -1e-3)
    deep = in_part & ~upper

    # --- element plastic strain, located by element centroid column
    conn = np.array(mesh["elements"], dtype=int)
    P = np.array(mesh["nodes_m"])
    cen = P[conn].mean(axis=1)
    ep = el[np.argsort(el[:, 0]), 1]
    ztc = target.interpolate(cen[:, 0], cen[:, 1], masked=False)
    e_part = ztc < -1e-6
    e_upper = e_part & (ztc > -1e-3)
    layer_top = cen[:, 2] > 0.5 * (zs[0] + zs[-1])
    res = {"part": part, "variant": variant, "step": step}

    def pct(v):
        v = np.asarray(v)
        if not len(v):
            return None
        return {"max": float(v.max()), "p95": float(np.percentile(v, 95)),
                "p50": float(np.percentile(v, 50)), "mean": float(v.mean())}

    res["eq_plastic_strain"] = {
        "all": pct(ep), "part": pct(ep[e_part]), "rim_band": pct(ep[e_upper]),
        "deeper": pct(ep[e_part & ~e_upper]), "flange": pct(ep[~e_part]),
        "top_layer_part": pct(ep[e_part & layer_top]),
        "bottom_layer_part": pct(ep[e_part & ~layer_top]),
        "share_part_beyond_0.19": float(np.mean(ep[e_part] > 0.19)),
        "share_part_beyond_0.1": float(np.mean(ep[e_part] > 0.10)),
    }

    # --- mid-surface membrane strains per cell (average of the layer surfaces)
    mid = G[nz // 2]
    mid0 = G0[nz // 2]
    h = xs[1] - xs[0]
    # cell tangents (forward differences averaged over the cell)
    a1 = 0.5 * ((mid[:-1, 1:] - mid[:-1, :-1]) + (mid[1:, 1:] - mid[1:, :-1])) / h
    a2 = 0.5 * ((mid[1:, :-1] - mid[:-1, :-1]) + (mid[1:, 1:] - mid[:-1, 1:])) / h
    C11 = np.einsum("...i,...i", a1, a1)
    C22 = np.einsum("...i,...i", a2, a2)
    C12 = np.einsum("...i,...i", a1, a2)
    tr = C11 + C22
    disc = np.sqrt(np.maximum(0.25 * (C11 - C22) ** 2 + C12 ** 2, 0))
    l1 = 0.5 * tr + disc
    l2 = 0.5 * tr - disc
    e1 = 0.5 * np.log(l1)
    e2 = 0.5 * np.log(l2)
    beta = np.where(np.abs(e1) > 1e-6, e2 / e1, np.nan)
    ccen = 0.5 * (mid0[:-1, :-1] + mid0[1:, 1:])
    zc = target.interpolate(ccen[..., 0].ravel(), ccen[..., 1].ravel(),
                            masked=False).reshape(ccen.shape[:2])
    cpart = zc < -1e-6
    wall = cpart & (e1 > 0.05)
    # net (path-independent) equivalent strain of the membrane state per cell,
    # against the accumulated plastic strain of the elements of that column
    e3 = -(e1 + e2)
    eq_net = np.sqrt(2.0 / 3.0 * (e1 ** 2 + e2 ** 2 + e3 ** 2))
    col_ep = np.zeros_like(eq_net)
    cnt = np.zeros_like(eq_net)
    ix = np.clip(np.searchsorted(xs, cen[:, 0]) - 1, 0, len(xs) - 2)
    iy = np.clip(np.searchsorted(ys, cen[:, 1]) - 1, 0, len(ys) - 2)
    np.add.at(col_ep, (iy, ix), ep)
    np.add.at(cnt, (iy, ix), 1)
    col_ep /= np.maximum(cnt, 1)
    strained = cpart & (eq_net > 0.02)
    res["accumulated_vs_net"] = {
        "net_membrane_eq_strain_part": pct(eq_net[cpart]),
        "column_mean_eps_p_part": pct(col_ep[cpart]),
        "ratio_accumulated_over_net_where_net>0.02": pct(col_ep[strained] / eq_net[strained])
        if strained.any() else None,
        "ratio_of_maxima": float(col_ep[cpart].max() / eq_net[cpart].max()),
    }
    res["membrane"] = {
        "major_log_strain_part": pct(e1[cpart]),
        "minor_log_strain_part": pct(e2[cpart]),
        "beta_where_eps1>0.05": pct(beta[wall]) if wall.any() else None,
        "share_cells_eps1>0.05_with_|beta|<0.25": float(np.mean(np.abs(beta[wall]) < 0.25))
        if wall.any() else None,
    }

    # --- through-thickness shear: column direction vs mid-surface normal
    n = np.cross(np.gradient(mid, h, axis=1), np.gradient(mid, h, axis=0))
    n /= np.linalg.norm(n, axis=-1, keepdims=True)

    def angle(f):
        f = f / np.linalg.norm(f, axis=-1, keepdims=True)
        c = np.clip(np.einsum("...i,...i", f, n), -1, 1)
        return np.arccos(c)

    g_total = angle(G[-1] - G[0])
    g_layers = [angle(G[k + 1] - G[k]) for k in range(nz - 1)]
    thick = np.linalg.norm(G[-1] - G[0], axis=-1)
    eps_t = np.log(thick / t0)
    res["transverse_shear_rad"] = {
        "column_part": pct(g_total[in_part]), "column_rim_band": pct(g_total[upper]),
        "column_deeper": pct(g_total[deep]),
        "per_layer_part": [pct(g[in_part]) for g in g_layers],
        "layer_difference_part": pct(np.abs(g_layers[-1] - g_layers[0])[in_part]),
    }
    res["thickness_log_strain"] = {"part": pct(-eps_t[in_part]), "min_t_mm": float(thick.min() * 1e3)}
    (HERE / f"validity_{part}_{variant}_{step}.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
