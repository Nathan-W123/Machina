"""Measurement-state study on the 8 benchmark test parts (cached uncompensated
runs, no new simulation): the deviation from the target of
  * form    : loaded, tool in place (what the sheet is at the end of forming)
  * unload  : tool removed, still clamped in the frame - what an in-cell scan
              (Machina: laser profiler on the robot before unclamping) sees
  * release : after unclamping onto 3-2-1 supports (the benchmark metric)
and the release surface after a best fit: 3-DOF (z shift + two tilts, least
squares over the part) and 6-DOF robust point-to-plane ICP over the whole
scanned surface (precomp.metrology.align, mode rigid) - the "aligned with the
nominal CAD model" QA practice. Writes measure_states.csv."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import proclens as P                       # noqa: E402
from precomp.metrology import align, part_mask  # noqa: E402

rows = []
for part in P.TEST_PARTS:
    d = P.RUNS / part / "base"
    m = P.measure(part, d)
    _, target = P.part_inputs(part)
    S = np.load(d / "surfaces.npz")
    X, Y = target.grid.mesh()
    pm = part_mask(target)
    # 6-DOF ICP of the whole released tool-side surface onto the target
    pts = np.column_stack([X.ravel(), Y.ravel(), S["release"].ravel()])
    al = align(pts, target, mode="rigid", initial="identity")
    moved = al.apply(pts)
    from precomp.geometry.heightmap import HeightMap
    hm = HeightMap.from_points(moved, target.grid)
    dev6 = (hm.z - target.z) / 1e-3
    ok = pm & hm.mask
    # ICP on the part only (exclude the flange)
    pts_p = pts[pm.ravel()]
    alp = align(pts_p, target, mode="rigid", initial="identity")
    hmp = HeightMap.from_points(alp.apply(pts), target.grid)
    devp = (hmp.z - target.z) / 1e-3
    okp = pm & hmp.mask
    row = {"part": part,
           "form_rms": m["form_part"]["rms"], "unload_rms": m["unload_part"]["rms"],
           "release_rms": m["release_part"]["rms"],
           "unload_upper_bias": m["unload_upper"]["bias"], "release_upper_bias": m["release_upper"]["bias"],
           "form_upper_bias": m["form_upper"]["bias"],
           "release_bestfit3_rms": m["release_part_bestfit3"]["rms"],
           "release_icp6_all_rms": float(np.sqrt(np.mean(dev6[ok] ** 2))),
           "icp6_all_rot_deg": al.rotation_deg, "icp6_all_tz_mm": al.translation[2] / 1e-3,
           "release_icp6_part_rms": float(np.sqrt(np.mean(devp[okp] ** 2))),
           "icp6_part_rot_deg": alp.rotation_deg, "icp6_part_tz_mm": alp.translation[2] / 1e-3,
           "springback_free_rms": float(np.sqrt(np.mean(((S["release"] - S["form"])[pm] / 1e-3) ** 2))),
           "springback_clamped_rms": float(np.sqrt(np.mean(((S["unload"] - S["form"])[pm] / 1e-3) ** 2))),
           "fz_mean_N": m["force_tool"]["fz_mean"], "fz_max_N": m["force_tool"]["fz_max"],
           "fxy_mean_N": m["force_tool"]["fxy_mean"], "fxy_max_N": m["force_tool"]["fxy_max"]}
    rows.append(row)
    print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
df = pd.DataFrame(rows)
df.loc["mean"] = df.mean(numeric_only=True)
df.to_csv(HERE / "measure_states.csv", float_format="%.4f")
print(df.round(3).to_string())
