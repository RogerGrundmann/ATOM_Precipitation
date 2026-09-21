#!/usr/bin/env python3
# ocprofile.py -- the ocean velocity PROFILE, which is the score that matters for the
# horizontal-metric question, and `max` is the one that misleads: the retracted 2026-09-04
# "ACC jet appears at 1.12 m/s" was ONE noisy cell while p50/p90/p99 were identical between
# the arms. Scores a FIXED population of full-depth columns -- set by the FIRST argument and
# reused for all the others -- so every level and every arm scores the same cells and a
# difference cannot be the land mask.
#   python3 ocprofile.py ctl:output_oc_ctl/hyd_restart_0Ma_1000.bin met:output_oc_met/...
# Reads the restart binary directly: header {magic 0x4F434D31, im, jm, km, total_iter_count}
# then 17 arrays [im][jm][km] doubles in the order FileIO_Hyd.cpp::restart_arrays() writes.
# i = 0 is the SEAFLOOR and i = im-1 the surface; u is RADIAL, v meridional, w zonal.
import numpy as np, sys, math
IM,JM,KM = 41,181,361
NARR = 17
NAMES = "t u v w c tn un vn wn cn p_dyn p_hydro tke dis tken disn nue".split()
U0 = 200.0/833.3          # m/s  (L_hyd / (L_hyd/u_0)), from the model's own [SCALES]
LH, BETA = 200.0, 2.0
def depth(i):             # i=40 surface, i=0 seafloor; inverse of FileIO_Hyd.cpp:411
    return -LH*math.sinh(BETA*(1.0-i/(IM-1)))/math.sinh(BETA)

def load(path):
    hdr = np.fromfile(path, dtype=np.int32, count=5)
    assert hdr[0] == 0x4F434D31, (path, hex(hdr[0]))
    assert (hdr[1],hdr[2],hdr[3]) == (IM,JM,KM), hdr
    a = np.memmap(path, dtype=np.float64, mode='r', offset=20,
                  shape=(NARR,IM,JM,KM))
    return hdr[4], a

def analyse(tag, path, colmask=None):
    it, a = load(path)
    u,v,w,c = a[1],a[2],a[3],a[4]
    water = c > 1e-6                              # salinity marks fluid cells
    if colmask is None:
        colmask = water.all(axis=0)               # FULL-DEPTH columns only
    speed = np.sqrt(v*v + w*w)*U0*100.0           # cm/s
    print(f"### {tag}   iter {it}   full-depth columns {colmask.sum()}")
    print(f"  {'depth':>8} {'median':>9} {'rms':>9} {'p90':>9}   (horizontal speed, cm/s)")
    for i in (40,39,35,30,20,10,0):
        s = speed[i][colmask]
        print(f"  {depth(i):8.1f} {np.median(s):9.3f} {np.sqrt((s**2).mean()):9.3f} {np.percentile(s,90):9.3f}")
    ur = np.abs(u)[:,colmask]; hz = speed[:,colmask]/(U0*100.0)
    print(f"  rms |u| (radial, nondim) {np.sqrt((ur**2).mean()):.6e}"
          f"   |u|/|horiz| rms {np.sqrt((ur**2).mean())/np.sqrt((hz**2).mean()):.6e}")
    sw = w[40]                                     # grid-scale noise on surface zonal velocity
    lap = (sw[2:,1:-1]+sw[:-2,1:-1]+sw[1:-1,2:]+sw[1:-1,:-2]-4*sw[1:-1,1:-1])
    m = colmask[1:-1,1:-1]
    print(f"  grid-scale noise index {np.sqrt((lap[m]**2).mean())/np.sqrt((sw[1:-1,1:-1][m]**2).mean()):.4f}")
    print()
    return colmask

if __name__ == "__main__":
    cm = None
    for arg in sys.argv[1:]:
        tag, path = arg.split(":",1)
        cm = analyse(tag, path, cm)     # ONE fixed population, set by the first arm
