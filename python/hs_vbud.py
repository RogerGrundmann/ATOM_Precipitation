#!/usr/bin/env python3
"""Geostrophic balance off the zonal-mean v budget: 20-70 deg, above 3 km.
residual = |pgf + coriolis| / |coriolis|  (1 = no balance, 0 = geostrophic)."""
import sys, csv, numpy as np
for f in sys.argv[1:]:
    P=[];C=[]
    for r in csv.DictReader(open(f)):
        la, z = abs(float(r['lat_deg'])), float(r['height_m'])
        c = float(r['coriolis'])
        if 20 <= la <= 70 and z > 3000 and c != 0: P.append(float(r['pgf'])); C.append(c)
    P=np.array(P); C=np.array(C); res=np.abs(P+C)/np.abs(C)
    print(f"{f}: n={len(C)}  median |pgf|/|cor|={np.median(np.abs(P)/np.abs(C)):.3e}  "
          f"residual median={np.median(res):.3f}  p05={np.percentile(res,5):.3f}  "
          f"frac(sign opposes cor)={np.mean(np.sign(P)!=np.sign(C)):.2f}")
