"""Unclamping springback (release minus unload, vertical, tool-side surface, |x|<=75 mm) per case -> springback.csv.
usage: python3 springback_table.py cases/<a> ..."""
import sys, numpy as np, csv
w = csv.writer(open('springback.csv', 'w')); w.writerow(['case', 'section', 'springback_z_rms_mm', 'springback_z_max_mm'])
for c in sys.argv[1:]:
    for t in ('RD', 'TD'):
        u = np.genfromtxt(f'{c}/profile_unload_{t}_top.csv', delimiter=',', skip_header=1)
        r = np.genfromtxt(f'{c}/profile_release_{t}_top.csv', delimiter=',', skip_header=1)
        xs = np.linspace(-75, 75, 301); d = np.interp(xs, r[:, 0], r[:, 1]) - np.interp(xs, u[:, 0], u[:, 1])
        w.writerow([c.split('/')[-1], t, round(float(np.sqrt((d ** 2).mean())), 3), round(float(abs(d).max()), 3)])
