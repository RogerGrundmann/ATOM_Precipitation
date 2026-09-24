#!/usr/bin/env python3
# Per-level rms horizontal speed on the fixed full-depth column set (from the first run), every
# level i = 0..12 plus 20/30/39/40. i = 0 is the 200 m TRUNCATION boundary (v/w cubic-extrapolated),
# NOT a prognostic level -- read the profile on i >= 1.
import numpy as np, sys
sys.argv = ['x']; exec(open('ocprofile.py').read().split("def analyse")[0])
runs = [('ctl','output_bn_ctl'),('b1','output_bn_1'),('b01','output_bn_01'),('sh1','output_bn_sh1')]
res, col = {}, None
for tag, d in runs:
    try: it, a = load(f'{d}/hyd_restart_0Ma_1600.bin')
    except Exception as e: print(f'{tag}: {e}'); continue
    if col is None: col = (a[4] > 1e-6).all(axis=0)
    sp = np.sqrt(a[2]**2 + a[3]**2) * U0 * 100
    res[tag] = [np.sqrt((sp[i][col]**2).mean()) for i in range(IM)]
    ru = np.sqrt((a[1][:, col]**2).mean()) * U0
    print(f'{tag}: iter {it}  rms|u_radial| {ru:.3e} m/s  max|u| {abs(a[1]).max()*U0:.3e} m/s')
print(' i   depth ' + ''.join(f'{t:>8}' for t in res) + '   (rms horizontal speed cm/s; i=0 = extrapolated BC)')
for i in list(range(13)) + [20, 30, 39, 40]:
    print(f'{i:2d} {depth(i):7.1f} ' + ''.join(f'{res[t][i]:8.3f}' for t in res))
for t in res: print(f'{t}: i=1 / i=12 = {res[t][1]/res[t][12]:.3f}')
