"""Regional normal-distance errors (floor |x|<15, corner 15-25, wall 25-58, rim/flange 58-75 mm) of a case's
profiles against the measurement (markers) and the published FE. usage: python3 regions_bm3.py cases/<a> ... -> regions.csv (appended)"""
import sys, os, csv
import numpy as np
import analyze_bm3 as A
REG = {'floor': (0, 15), 'corner': (15, 25), 'wall': (25, 58), 'rim': (58, 75)}
def nd(pts, q):
    out = []
    for pt in pts:
        a, b = q[:-1], q[1:]; ab = b - a; L2 = np.maximum((ab ** 2).sum(1), 1e-12)
        t = np.clip(((pt - a) * ab).sum(1) / L2, 0, 1); pr = a + t[:, None] * ab
        d = np.sqrt(((pr - pt) ** 2).sum(1)); i = np.argmin(d); out.append(d[i] * (1 if pt[1] >= pr[i, 1] else -1))
    return np.array(out)
rows = []
def evalp(name, p, e, tag, step):
    for refname, ref in (('measured', e['exp_marker']), ('paperFE', None)):
        if ref is None:
            q = e['sim']; xs = np.linspace(-75, 75, 601); ref = np.stack([xs, np.interp(xs, q[:, 0], q[:, 1])], 1)
        for r, (lo, hi) in REG.items():
            m = (np.abs(ref[:, 0]) >= lo) & (np.abs(ref[:, 0]) < hi)
            d = nd(ref[m], p)   # + : reference above the model
            rows.append({'case': name, 'step': step, 'section': tag, 'reference': refname, 'region': r,
                         'normal_rms_mm': round(float(np.sqrt(np.mean(d ** 2))), 3), 'normal_mean_mm': round(float(np.mean(d)), 3), 'n': int(m.sum())})
for tag in ('RD', 'TD'):
    e = A.load_exp(tag)
    evalp('paperFE', e['sim'], e, tag, '-')
for c in sys.argv[1:]:
    for step in ('unload', 'release'):
        for tag in ('RD', 'TD'):
            f = os.path.join(c, f'profile_{step}_{tag}_top.csv')
            if os.path.exists(f): evalp(os.path.basename(c.rstrip('/')), np.loadtxt(f, delimiter=',', skiprows=1), A.load_exp(tag), tag, step)
new = not os.path.exists('regions.csv')
with open('regions.csv', 'a', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    if new: w.writeheader()
    w.writerows(rows)
for r in rows:
    if r['section'] == 'RD': print(r)
