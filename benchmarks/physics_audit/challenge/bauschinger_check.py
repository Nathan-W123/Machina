"""Uniaxial tension -> compression response of the material lens' kin_shutov card
(two AF backstresses + linear+Voce isotropic, values from its deck.json) vs the
library isotropic card: reverse 0.2%-offset yield after prestrain, ratio to the
forward flow stress. Return-mapping in 1D (rate-independent, small strain)."""
import json
import numpy as np
from pathlib import Path
A = Path(__file__).resolve().parent.parent
def card(v):
    m = json.loads((A / f"material/runs/truncated_cone-s2026-0000/{v}/deck.json").read_text())["material"]
    p = m["plasticity"]
    return dict(E=m["youngs_modulus"] / 1e6, s0=p["yield_stress"] / 1e6, H=p["hardening_modulus"] / 1e6,
                Q=p["saturation_stress"] / 1e6, b=p["saturation_rate"],
                bs=[(x["modulus"] / 1e6, x["recovery"]) for x in p.get("backstresses", [])])
def R(c, p): return c["s0"] + c["H"] * p + c["Q"] * (1 - np.exp(-c["b"] * p))
def run(c, pre, back=0.03, n=20000):
    eps = np.concatenate([np.linspace(0, pre, n), np.linspace(pre, pre - back, n)[1:]])
    ep, p, X = 0.0, 0.0, np.zeros(len(c["bs"]))
    out = []
    for e in eps:
        for _ in range(50):  # simple fixed-point return mapping on dp
            s_tr = c["E"] * (e - ep)
            xi = s_tr - X.sum()
            f = abs(xi) - R(c, p)
            if f <= 1e-9: break
            sg = np.sign(xi)
            Hk = sum(C - g * sg * Xi for (C, g), Xi in zip(c["bs"], X))
            Hi = c["H"] + c["Q"] * c["b"] * np.exp(-c["b"] * p)
            dp = f / (c["E"] + Hk + Hi)
            ep += sg * dp; p += dp
            X = np.array([Xi + (C * sg - g * Xi) * dp for (C, g), Xi in zip(c["bs"], X)])
        out.append((e, c["E"] * (e - ep)))
    return np.array(out)
res = {}
for v in ["base", "kin_shutov"]:
    c = card(v)
    for pre in [0.02, 0.05, 0.1, 0.2]:
        o = run(c, pre)
        k = np.argmax(o[:, 0] >= pre - 1e-12)
        sf = o[k, 1]
        rev = o[k:]
        # 0.2% offset reverse yield: |s| on unloading line minus E*0.002
        e_unl = rev[:, 0]
        lin = sf + c["E"] * (e_unl - pre)
        idx = np.argmax(rev[:, 1] - (lin + c["E"] * 0.002) > 0) if np.any(rev[:, 1] > lin + c["E"] * 0.002) else None
        # reverse yield = first point where stress deviates from elastic line by 0.2% offset
        dev = rev[:, 1] - (sf + c["E"] * (e_unl - pre))
        j = np.argmax(dev > c["E"] * 0.002)
        sr = rev[j, 1]
        # proportional limit: first deviation of 1 MPa
        j1 = np.argmax(dev > 1.0)
        res[f"{v}_pre{pre}"] = {"forward_MPa": round(float(sf), 1), "reverse_0.2pct_MPa": round(float(sr), 1),
                               "ratio_0.2pct": round(float(-sr / sf), 3),
                               "reverse_prop_limit_MPa": round(float(rev[j1, 1]), 1),
                               "elastic_range_MPa": round(float(sf - rev[j1, 1]), 1)}
print(json.dumps(res, indent=1))
json.dump(res, open(Path(__file__).with_name("bauschinger_check.json"), "w"), indent=1)
