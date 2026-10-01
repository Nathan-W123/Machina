"""Check the fixed-point residual of the process lens' robot-compliance iteration:
forces fed into rob_iso1_it5 vs forces it produced, along the path (by path fraction)."""
import sys, json
from pathlib import Path
import numpy as np
P_DIR = Path(__file__).resolve().parent.parent / "process"
sys.path.insert(0, str(P_DIR))
import proclens as P
part = "truncated_cone-s2026-0000"
setup, commanded = P.part_inputs(part)
path = P.make_toolpath(setup, commanded)
out = {}
for it in ["rob_iso1_it4", "rob_iso1_it5"]:
    fin = np.load(P.RUNS / part / it / "robot_forces_in.npy")
    fout = P.input_forces(part, it, path)
    contact = path.level != P.AIR
    s = path.s / path.s[contact].max()
    res = {}
    for lo, hi in [(0, .5), (.5, .8), (.8, 1.01)]:
        m = contact & (s >= lo) & (s < hi)
        res[f"{lo}-{hi}"] = {"fz_in_mean": float(fin[m, 2].mean()), "fz_out_mean": float(fout[m, 2].mean()),
                             "dz_in_um_mean": float(fin[m, 2].mean()),  # 1 um/N
                             "rel_resid_rms": float(np.sqrt(np.mean(np.sum((fout[m]-fin[m])**2, 1))) / np.sqrt(np.mean(np.sum(fin[m]**2, 1))))}
    out[it] = res
print(json.dumps(out, indent=1))
json.dump(out, open(Path(__file__).with_name("robot_fixedpoint.json"), "w"), indent=1)
