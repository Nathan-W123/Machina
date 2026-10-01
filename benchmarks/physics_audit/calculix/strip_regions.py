"""Release-shape error of each SparLab strip case against a CalculiX
reference, split into the formed region (|x| < 8 mm), the wall/clamp region
(8..15) and the flange (15..19).  usage: python3 strip_regions.py REF"""
import sys, csv
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import strip_reference as sr
from common import ROOT
ref = sys.argv[1]
rdir = sr.OUT / ref
import json
W = json.load(open(rdir / "params.json"))["h"]
rprof, rforce = sr.ccx_profiles(rdir / "ccx", W)
rows = []
for case in sorted((ROOT / "punch").glob("strip_*")):
    if not (case / "result.json").exists():
        continue
    prof, f = sr.sparlab_profiles(case)
    r = {"case": case.name}
    for st in ("form", "unload", "release"):
        x = prof[st][:, 0]
        d = prof[st][:, 1] - np.interp(x, rprof[st][:, 0], rprof[st][:, 1])
        for lo, hi in ((0, 8), (8, 15), (15, 19.01)):
            m = (np.abs(x) >= lo) & (np.abs(x) < hi)
            r[f"{st}_rms_{lo}_{int(hi)}"] = float(np.sqrt(np.mean(d[m] ** 2)))
            r[f"{st}_max_{lo}_{int(hi)}"] = float(np.max(np.abs(d[m])))
    rows.append(r)
with open(sr.OUT / f"regions_{ref}.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
for r in rows:
    print(f"{r['case']:28s} release rms/max  core {r['release_rms_0_8']:.4f}/{r['release_max_0_8']:.4f}  wall {r['release_rms_8_15']:.4f}/{r['release_max_8_15']:.4f}  flange {r['release_rms_15_19']:.4f}/{r['release_max_15_19']:.4f}")
