"""Fit sigma = s0 + H e + Q (1 - exp(-b e)) to Swift 343.3 (0.0015 + e)^0.184 MPa (AA7075-O, Neto et al. 2016 Table 1)
over e in [0, 1.2] (SPIF reaches ~1). Writes swift_fit.json and the fit error."""
import json, numpy as np
from scipy.optimize import least_squares
K, n, e0 = 343.3e6, 0.184, 0.0015
e = np.concatenate([np.linspace(0, 0.05, 60), np.linspace(0.05, 1.2, 120)])
sw = K * (e0 + e) ** n
s0 = K * e0 ** n
f = lambda p: (s0 + p[0] * 1e6 + p[1] * 1e6 * (1 - np.exp(-p[2] * e)) - 0) if False else (s0 + p[0]*1e6*e + p[1]*1e6*(1-np.exp(-p[2]*e)) - sw) / sw
r = least_squares(f, [50, 150, 10], bounds=([0, 0, 0.1], [2000, 1000, 500]))
H, Q, b = r.x[0]*1e6, r.x[1]*1e6, r.x[2]
err = f(r.x)
out = {"sigma0": s0, "H": H, "Q": Q, "b": b, "max_rel_err": float(np.max(np.abs(err))),
       "rms_rel_err": float(np.sqrt(np.mean(err**2))), "range": [0, 1.2],
       "check_MPa": {str(x): [float(K*(e0+x)**n/1e6), float((s0+H*x+Q*(1-np.exp(-b*x)))/1e6)] for x in [0, 0.002, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0]}}
json.dump(out, open('swift_fit.json', 'w'), indent=1); print(json.dumps(out, indent=1))
