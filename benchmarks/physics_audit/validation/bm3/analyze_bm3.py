"""Compare a BM3 replica run with the digitised experiment (Neto et al. 2016, Figs. 8 and 10).
usage: python3 analyze_bm3.py cases/<case> [...]; appends a row per case to results.csv and writes
<case>/profile_<step>.csv and <case>/force_by_depth.csv."""
import sys, os, csv, json, math, glob
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); DIG = os.path.join(HERE, '..', 'digitize')
T = 1.63; R = 6.33
def load_exp(tag):
    rows = list(csv.DictReader(open(os.path.join(DIG, f'bm3_fig10_profile_{tag}_clean.csv'))))
    ser = {}
    for r in rows: ser.setdefault(r['series'], []).append((float(r['x_mm']), float(r['z_mm'])))
    return {k: np.array(sorted(v)) for k, v in ser.items()}
def model_profile(case, step, angle_deg, surface):
    f = glob.glob(os.path.join(case, 'output', f'step_*_{step}_nodes.csv'))
    if not f: return None
    a = np.genfromtxt(f[0], delimiter=',', names=True)
    X, Y, Z = a['X'] * 1e3, a['Y'] * 1e3, a['Z'] * 1e3
    x, y, z = X + a['ux'] * 1e3, Y + a['uy'] * 1e3, Z + a['uz'] * 1e3
    zt = Z.max() if surface == 'top' else Z.min()
    sel = np.abs(Z - zt) < 1e-6
    out = []
    for sgn in (1, -1):
        th = math.radians(angle_deg) + (0 if sgn == 1 else math.pi)
        c, s = math.cos(th), math.sin(th)
        # nodes whose reference position lies on the ray (O-grid: exact rays at 0/90/180/270 deg)
        rr = X * c + Y * s; perp = -X * s + Y * c
        m = sel & (np.abs(perp) < 1e-3) & (rr >= -1e-6)
        # deformed coordinate along the section line; a node may twist off it slightly (tangential motion)
        xs = (x[m] * c + y[m] * s) * sgn; zs = z[m] - (0 if surface == 'top' else -T)
        out += list(zip(xs, zs))
    p = np.array(sorted(out))
    return p

def model_thickness(case, step, angle_deg):
    f = glob.glob(os.path.join(case, 'output', f'step_*_{step}_nodes.csv'))
    if not f: return None
    a = np.genfromtxt(f[0], delimiter=',', names=True)
    X, Y, Z = a['X'] * 1e3, a['Y'] * 1e3, a['Z'] * 1e3
    x = np.stack([X + a['ux'] * 1e3, Y + a['uy'] * 1e3, Z + a['uz'] * 1e3], 1)
    top = np.abs(Z - Z.max()) < 1e-6; bot = np.abs(Z - Z.min()) < 1e-6
    key = lambda i: (round(X[i], 5), round(Y[i], 5))
    bi = {key(i): i for i in np.where(bot)[0]}
    out = []
    for sgn in (1, -1):
        th = math.radians(angle_deg) + (0 if sgn == 1 else math.pi); c, s = math.cos(th), math.sin(th)
        for i in np.where(top)[0]:
            rr = X[i] * c + Y[i] * s; perp = -X[i] * s + Y[i] * c
            if abs(perp) < 1e-3 and rr >= -1e-6 and key(i) in bi:
                j = bi[key(i)]
                out.append(((x[i, 0] * c + x[i, 1] * s) * sgn, float(np.linalg.norm(x[i] - x[j]))))
    return np.array(sorted(out))
def compare(p, exp, xlim=75.0):
    xe, ze = exp[:, 0], exp[:, 1]
    m = np.abs(xe) <= xlim
    zm = np.interp(xe[m], p[:, 0], p[:, 1])
    d = zm - ze[m]
    # horizontal error of the wall at mid-depth: x where z crosses -22 on each side
    def xcross(q, zc=-22.0, side=1):
        q = q[q[:, 0] * side > 0]; q = q[np.argsort(np.abs(q[:, 0]))]
        for i in range(len(q) - 1):
            if (q[i, 1] - zc) * (q[i + 1, 1] - zc) <= 0:
                t = (zc - q[i, 1]) / (q[i + 1, 1] - q[i, 1]); return abs(q[i, 0] + t * (q[i + 1, 0] - q[i, 0]))
        return float('nan')
    floor = lambda q: float(np.mean(q[np.abs(q[:, 0]) < 10, 1])) if (np.abs(q[:, 0]) < 10).any() else float('nan')
    # shortest distance from each measured point to the model polyline (normal error; insensitive to the wall slope)
    def pdist(pt, q):
        a, b = q[:-1], q[1:]; ab = b - a; L2 = np.maximum((ab ** 2).sum(1), 1e-12)
        t = np.clip(((pt - a) * ab).sum(1) / L2, 0, 1); proj = a + t[:, None] * ab
        dd = np.sqrt(((proj - pt) ** 2).sum(1)); i = np.argmin(dd)
        # signed: + when the measured point lies above the model
        return dd[i] * (1 if pt[1] >= proj[i, 1] else -1)
    nd = np.array([pdist(pt, p) for pt in exp[m]])
    return {'rms_mm': float(np.sqrt(np.mean(d ** 2))), 'max_mm': float(np.max(np.abs(d))), 'bias_mm': float(np.mean(d)),
            'floor_model_mm': floor(p), 'floor_exp_mm': floor(exp),
            'wall_r_mid_model_mm': 0.5 * (xcross(p, side=1) + xcross(p, side=-1)),
            'wall_r_mid_exp_mm': 0.5 * (xcross(exp, side=1) + xcross(exp, side=-1)),
            'wall_x_mid_model_pos_neg': [xcross(p, side=1), xcross(p, side=-1)], 'wall_x_mid_exp_pos_neg': [xcross(exp, side=1), xcross(exp, side=-1)],
            'normal_rms_mm': float(np.sqrt(np.mean(nd ** 2))), 'normal_max_mm': float(np.max(np.abs(nd))), 'normal_bias_mm': float(np.mean(nd)),
            'n': int(m.sum())}
def forces(case, dz):
    f = os.path.join(case, 'output', 'tool_forces.csv')
    if not os.path.exists(f): return None
    a = np.genfromtxt(f, delimiter=',', names=True, dtype=None, encoding=None)
    a = a[a['tool'] == 'tool']
    tp = np.genfromtxt(os.path.join(case, 'toolpath.csv'), delimiter=',', names=True)
    # contour index from the tool-centre z (tip depth = R - z_c)
    depth = np.ceil((R - a['cz'] * 1e3) / dz - 1e-6) * dz
    rows = []
    for d in sorted(set(depth)):
        if d <= 0 or d > 44.01: continue
        m = (depth == d)
        rows.append((d, float(np.mean(a['fz'][m])), float(np.max(a['fz'][m])), float(np.mean(np.hypot(a['fx'][m], a['fy'][m]))), int(m.sum())))
    return rows
def exp_force_by_depth():
    rows = list(csv.DictReader(open(os.path.join(DIG, 'bm3_fig8_forces.csv'))))
    out = {}
    for ser in ('exp_Fz_line', 'sim_Fz'):
        a = np.array(sorted((float(r['t_s']), float(r['F_N'])) for r in rows if r['series'] == ser))
        # time -> tool depth, assuming constant feed along the 88 contours (0.5 mm step, r_c = 62 - 0.5 (k-1))
        L = np.array([2 * math.pi * (62 - 0.5 * k) for k in range(88)]); cum = np.concatenate([[0], np.cumsum(L)])
        tmax = 1330.0
        s = a[:, 0] / tmax * cum[-1]
        k = np.clip(np.searchsorted(cum, s, side='right'), 1, 88)
        d = 0.5 * k
        tu = np.linspace(0, tmax, 20000); Fu = np.interp(tu, a[:, 0], a[:, 1]); ku = np.clip(np.searchsorted(cum, tu / tmax * cum[-1], side='right'), 1, 88)
        out[ser] = {float(0.5 * kk): float(Fu[ku == kk].mean()) for kk in range(1, 89) if (ku == kk).any()}
    return out
def scallop(dz, alpha=45.0):
    p = dz / math.sin(math.radians(alpha)); return R - math.sqrt(R * R - (p / 2) ** 2)
def aerens(dz, Rm=198.0, t=1.63, dt=12.66, alpha=45.0):
    return 0.0716 * Rm * t ** 1.57 * dt ** 0.41 * scallop(dz) ** 0.09 * alpha * math.cos(math.radians(alpha))
if __name__ == '__main__':
    expF = exp_force_by_depth()
    json.dump(expF, open(os.path.join(HERE, 'exp_force_by_depth.json'), 'w'), indent=0)
    res = []
    for case in sys.argv[1:]:
        deck = json.load(open(os.path.join(case, 'deck.json')))
        dz = float(deck['description'].split('dz=')[1].split()[0])
        row = {'case': os.path.basename(case.rstrip('/')), 'dz_mm': dz}
        summ = os.path.join(case, 'output', 'summary.json')
        if os.path.exists(summ):
            s = json.load(open(summ)); row.update(runtime_s=s.get('runtime_s'), completed=s.get('completed'),
                increments=s['timing'].get('increments'), iterations=s['timing'].get('iterations'), cuts=s['timing'].get('cuts'))
        for step in ('form', 'unload', 'release'):
            for tag, ang in (('RD', 0.0), ('TD', 90.0)):
                exp = load_exp(tag)
                for surf in ('top', 'bottom'):
                    p = model_profile(case, step, ang, surf)
                    if p is None: continue
                    np.savetxt(os.path.join(case, f'profile_{step}_{tag}_{surf}.csv'), p, delimiter=',', header='x_mm,z_mm', comments='')
                    c = compare(p, exp['exp_marker'])
                    for k, v in c.items(): row[f'{step}_{tag}_{surf}_{k}'] = v
                    # code-to-code: against the published FE curve, sampled densely (its polyline)
                    q = exp['sim']; xs = np.linspace(-75, 75, 301); ref = np.stack([xs, np.interp(xs, q[:, 0], q[:, 1])], 1)
                    c2 = compare(p, ref)
                    for k in ('rms_mm', 'max_mm', 'normal_rms_mm', 'normal_max_mm', 'normal_bias_mm'): row[f'{step}_{tag}_{surf}_vsPaperFE_{k}'] = c2[k]
        thk = list(csv.DictReader(open(os.path.join(DIG, 'bm3_fig15_thickness.csv'))))
        et = np.array(sorted((float(r['x_mm']), float(r['t_mm'])) for r in thk if r['series'] == 'exp_marker'))
        for step in ('release',):
            for tag, ang, side in (('RD', 0.0, -1), ('TD', 90.0, 1)):
                q = model_thickness(case, step, ang)
                if q is None or len(q) == 0: continue
                np.savetxt(os.path.join(case, f'thickness_{step}_{tag}.csv'), q, delimiter=',', header='x_mm,t_mm', comments='')
                e = et[et[:, 0] * side > 0]
                e = e[np.abs(e[:, 0]) <= 75]
                tm = np.interp(e[:, 0], q[:, 0], q[:, 1]); dd = tm - e[:, 1]
                row[f'thk_{tag}_rms_mm'] = float(np.sqrt(np.mean(dd ** 2))); row[f'thk_{tag}_min_model_mm'] = float(q[np.abs(q[:, 0]) <= 75, 1].min())
                row[f'thk_{tag}_min_exp_mm'] = float(e[:, 1].min())
        F = forces(case, dz)
        if F:
            with open(os.path.join(case, 'force_by_depth.csv'), 'w') as f:
                f.write('depth_mm,fz_mean_N,fz_max_N,fh_mean_N,n,exp_fz_mean_N_at_depth,paper_sim_fz_N,aerens_factor_to_dz05\n')
                fac = (scallop(0.5) / scallop(dz)) ** 0.09
                for d, fm, fx, fh, n in F:
                    ke = min(expF['exp_Fz_line'], key=lambda q: abs(q - d))
                    f.write(f'{d},{fm:.1f},{fx:.1f},{fh:.1f},{n},{expF["exp_Fz_line"][ke]:.1f},{expF["sim_Fz"][ke]:.1f},{fac:.4f}\n')
            ss = [r for r in F if r[0] >= 24]
            if ss:
                row['fz_steady_model_N'] = float(np.mean([r[1] for r in ss]))
                row['fz_steady_model_scaled_to_dz05_N'] = row['fz_steady_model_N'] * (scallop(0.5) / scallop(dz)) ** 0.09
                row['fz_peak_model_N'] = float(max(r[2] for r in F))
        row['aerens_fz_dz05_N'] = aerens(0.5); row['aerens_fz_dz_N'] = aerens(dz)
        res.append(row)
        print(json.dumps(row, indent=1))
    out = os.path.join(HERE, 'results.jsonl')
    with open(out, 'a') as f:
        for r in res: f.write(json.dumps(r) + '\n')
