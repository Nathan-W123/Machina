"""Tool-centre deflection C f along the benchmark paths (no simulation):
forces from the cached uncompensated runs of the 8 test parts (and of all 207
cached runs for the force distribution), compliance models:
  iso1       : 1 um/N isotropic (typical articulated robot, "< 1 N/um")
  iso2       : 2 um/N isotropic (the effective compliance implied by the
               measured deviations of Bharti et al. 2024: ~1 mm at <= ~500 N)
  irb_vert   : IRB 7600 VJM (Bharti Table 3), Machina-style vertical sheet,
               reach 1.8 m, TCP 1.3 m high, 0.30 m tool
  irb_horiz  : IRB 7600 VJM, horizontal sheet below the robot, reach 1.8 m,
               0.3 m high, 0.30 m tool
Writes robot_deflection.csv (per part and model: mean / p95 / max of |delta|,
delta_z, |delta_xy| over the in-contact path) and prints the summary."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import proclens as P            # noqa: E402
import variants as V            # noqa: E402
from precomp.fea.deck import make_toolpath    # noqa: E402
from precomp.robot import CartesianCompliance, forces_on_path  # noqa: E402
from precomp.toolpath import AIR   # noqa: E402

MODELS = {"iso1": np.eye(3) * 1e-6, "iso2": np.eye(3) * 2e-6,
          "irb_vert": np.array(V.C_VERT), "irb_horiz": np.array(V.C_HORIZ)}
rows = []
for part in P.TEST_PARTS:
    s, c = P.part_inputs(part)
    path = make_toolpath(s, c)
    f = forces_on_path(path, P.forces_table(part, "cache"), step=1)
    con = path.level != AIR
    for name, C in MODELS.items():
        d = CartesianCompliance(C).deflection(f)[con] / 1e-3
        n = np.linalg.norm(d, axis=1)
        dxy = np.linalg.norm(d[:, :2], axis=1)
        rows.append({"part": part, "model": name,
                     "norm_mean": n.mean(), "norm_p95": np.percentile(n, 95), "norm_max": n.max(),
                     "dz_mean": d[:, 2].mean(), "dz_max": d[:, 2].max(),
                     "dxy_mean": dxy.mean(), "dxy_max": dxy.max()})
df = pd.DataFrame(rows)
df.to_csv(HERE / "robot_deflection.csv", index=False, float_format="%.4f")
print(df.groupby("model").mean(numeric_only=True).round(3).to_string())
