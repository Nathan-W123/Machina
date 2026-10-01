"""Digitise Fig. 8 (tool force vs time) and Fig. 10 (profiles after unloading) of
Neto et al., Int J Adv Manuf Technol 85 (2016) 521-534, doi:10.1007/s00170-015-7954-9
(NUMISHEET 2014 BM3 data), from the vector paths of PDF page 7 (exact coordinates, no image digitising).
Axis calibration from the tick marks (see README)."""
import json, csv
d = json.load(open('p7.json'))
def inr(pts, r): return all(r[0] <= x <= r[2] and r[1] <= y <= r[3] for x, y in pts)
def sel(r, op, stroke=None, fill=None, lw=None):
    out = []
    for p in d['paths']:
        pts = [q for s in p['subs'] for q in s]
        if not inr(pts, r) or p['op'] != op: continue
        if stroke is not None and p['stroke'] != stroke: continue
        if fill is not None and p['fill'] != fill: continue
        if lw is not None and abs(p['lw'] - lw) > 0.01: continue
        out.append(p)
    return out
# Fig 8
R8 = (205, 620, 540, 735)
fx = lambda X: (X - 210.98) / (535.84 - 210.98) * 1350.0
fy = lambda Y: (Y - 627.32) / (731.47 - 627.32) * 2000.0
series = {'exp_Fz_line': ([0,0,0], 0.59), 'sim_Fz': ([1,0,0], 0.78), 'sim_Fr': ([0,0.8,0], 0.78), 'sim_Ft': ([0,0,1], 0.78)}
with open('bm3_fig8_forces.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['series', 't_s', 'F_N'])
    for name, (col, lw) in series.items():
        for p in sel(R8, 'S', col, None, lw):
            for s in p['subs']:
                if len(s) < 20: continue
                for x, y in s: w.writerow([name, round(fx(x), 3), round(fy(y), 2)])
    # filled black glyphs: experimental markers? record centroids
    for p in sel(R8, 'f', None, [0,0,0]):
        for s in p['subs']:
            xs = [q[0] for q in s]; ys = [q[1] for q in s]
            w.writerow(['exp_Fz_marker', round(fx(sum(xs)/len(xs)), 3), round(fy(sum(ys)/len(ys)), 2)])
# Fig 10
for tag, r, y0, y1 in [('RD', (330, 240, 545, 330), 247.05, 322.67), ('TD', (330, 125, 545, 215), 131.48, 207.10)]:
    gx = lambda X: (X - 337.7) / (538.81 - 337.7) * 180.0 - 90.0
    gy = lambda Y, y0=y0, y1=y1: (Y - y0) / (y1 - y0) * 60.0 - 50.0
    with open(f'bm3_fig10_profile_{tag}.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['series', 'x_mm', 'z_mm'])
        for name, col, lw in [('sim', [1,0,0], 0.96), ('cad', [0,0,1], 0.96), ('exp_line', [0,0,0], 0.72)]:
            for p in sel(r, 'S', col, None, lw):
                for s in p['subs']:
                    for x, y in s: w.writerow([name, round(gx(x), 3), round(gy(y), 3)])
        for p in sel(r, 'f', None, [0,0,0]):
            for s in p['subs']:
                xs = [q[0] for q in s]; ys = [q[1] for q in s]
                w.writerow(['exp_marker', round(gx(sum(xs)/len(xs)), 3), round(gy(sum(ys)/len(ys)), 3)])
print('ok')
