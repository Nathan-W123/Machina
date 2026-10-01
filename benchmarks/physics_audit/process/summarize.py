"""Print / write the variant results (last record per part and variant)."""
import json
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
recs = {}
for l in open(HERE / "results.jsonl"):
    r = json.loads(l)
    recs[(r["part"], r["variant"])] = r
rows = []
for (p, v), r in recs.items():
    if not r.get("ok"):
        rows.append({"part": p, "variant": v, "ok": False})
        continue
    vb = r.get("vs_base_release_part", {})
    row = {"part": p, "variant": v, "ok": True,
           "rel_rms": r["release_part"]["rms"], "rel_max": r["release_part"]["max"],
           "upper_bias": r["release_upper"]["bias"], "deep_bias": r["release_deep"]["bias"],
           "depth": r["release_depth_mm"], "form_rms": r["form_part"]["rms"],
           "form_upper_bias": r["form_upper"]["bias"],
           "unload_rms": r["unload_part"]["rms"], "bestfit3_rms": r["release_part_bestfit3"]["rms"],
           "d_rms": vb.get("rms"), "d_bias": vb.get("bias"), "d_max": vb.get("max"),
           "d_upper_rms": r.get("vs_base_release_upper", {}).get("rms"),
           "d_deep_rms": r.get("vs_base_release_deep", {}).get("rms"),
           "fz_mean": r["force_tool"]["fz_mean"], "fz_max": r["force_tool"]["fz_max"],
           "sup_fz_mean": r.get("force_support", {}).get("fz_mean"),
           "sup_fz_max": r.get("force_support", {}).get("fz_max"),
           "runtime_s": r["runtime_s"], "cuts": r["cuts"]}
    rows.append(row)
df = pd.DataFrame(rows).sort_values(["part", "variant"])
df.to_csv(HERE / "results.csv", index=False, float_format="%.4f")
pd.set_option("display.width", 250)
print(df.round(3).to_string(index=False))
