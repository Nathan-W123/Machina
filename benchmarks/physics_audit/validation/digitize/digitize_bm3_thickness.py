"""Digitise Fig. 15 (thickness vs distance from centreline; left half RD, right half TD) of Neto et al. 2016
from the vector paths of PDF page 9 (p9.json from extract_paths.py). Axis: ticks x -90..90 at 81.68..283.59 pt,
y 0.8..1.8 mm at 586.58..732.41 pt."""
import json, csv
d = json.load(open('p9.json'))
gx = lambda X: (X - 81.68) / (283.59 - 81.68) * 180 - 90
gy = lambda Y: (Y - 586.58) / (732.41 - 586.58) * 1.0 + 0.8
r = (60, 620, 300, 740)
inr = lambda pts: all(r[0] <= x <= r[2] and r[1] <= y <= r[3] for x, y in pts)
with open('bm3_fig15_thickness.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['series', 'x_mm', 't_mm'])
    for p in d['paths']:
        pts = [q for s in p['subs'] for q in s]
        if not inr(pts): continue
        if p['op'] == 'S' and p['stroke'] == [1.0, 0.0, 0.0]: name = 'sim'
        elif p['op'] == 'S' and p['stroke'] == [0.0, 0.0, 1.0]: name = 'sine_law'
        elif p['op'] == 'f' and p['fill'] == [0.0, 0.0, 0.0]:
            for s in p['subs']:
                xs = [q[0] for q in s]; ys = [q[1] for q in s]
                w.writerow(['exp_marker', round(gx(sum(xs) / len(xs)), 3), round(gy(sum(ys) / len(ys)), 4)])
            continue
        else: continue
        for s in p['subs']:
            for x, y in s: w.writerow([name, round(gx(x), 3), round(gy(y), 4)])
print('ok')
