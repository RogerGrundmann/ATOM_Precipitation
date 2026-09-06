#!/usr/bin/env python3
"""Supercooling census on a hydrosphere restart binary.

Reads hyd_restart_<Ma>Ma_<iter>.bin directly (5 x int32 header -- magic "OCM1", im, jm,
km, total_iter_count -- then the 17 prognostic arrays of save_state/restart_arrays, in
that order) and applies the model's OWN criterion for a supercooled cell: salinity
S = c * c_35 with c_35 = 34.6 psu, cells with S >= 5 psu only (the fresh-cell fallback
in ValueLimitationHyd), and the same UNESCO/Millero freezing point as
seawater_freezing_point_C in UtilsHyd.h.

VALIDATED against the recorded four-arm HYD_T_FREEZE / HYD_T_FREEZE_SFC table before it
was used on anything new: it reproduces f0, f1, s1 and s2 exactly -- every cell count and
every mean T.

ONE CORRECTION IT NEEDED, AND IT IS NOT A DETAIL. A strict T < T_f OVER-counts by 700
cells in f1 and 683 in s2, because a cell the interior floor WROTE sits at a deficit of
~1e-12 C -- round-off through the non-dimensionalisation. A cell held AT the floor is not
supercooled, so the test carries a 1e-9 C tolerance.

Land cells are identified as c == 0 (ValueLimitationHyd forces salinity to zero there).

    python3 census_freeze.py output_f0/hyd_restart_0Ma_500.bin [...]
"""
import numpy as np, sys, os
ARRS = ["t","u","v","w","c","tn","un","vn","wn","cn",
        "p_dyn","p_hydro","tke","dis","tken","disn","nue"]
T0, C35 = 273.15, 34.6

def load(path):
    with open(path,"rb") as f:
        hdr = np.fromfile(f, dtype=np.int32, count=5)
        magic, im, jm, km, it = hdr
        assert magic == 0x4F434D31, hex(magic)
        n = im*jm*km
        d = {}
        for a in ARRS:
            d[a] = np.fromfile(f, dtype=np.float64, count=n).reshape(im,jm,km)
    d["_dims"] = (int(im),int(jm),int(km),int(it))
    return d

def Tf(S):
    S = np.maximum(S, 0.0)
    return -0.0575*S + 1.710523e-3*S**1.5 - 2.154996e-4*S*S

def census(d):
    t, c = d["t"], d["c"]
    im = d["_dims"][0]
    S = c*C35
    fluid = c > 0.0                    # land cells are forced to c = 0
    salty = fluid & (S >= 5.0)
    tC = (t - 1.0)*T0                  # non-dim -> Celsius
    # A cell HELD AT the floor is not supercooled: the interior floor writes t exactly
    # to (t_0 + T_f)/t_0, and the round trip through the non-dimensionalisation leaves a
    # deficit of ~1e-12 C. Without this tolerance those cells (700 in f1, 683 in s2) are
    # miscounted as still supercooled, which is the whole difference between this script
    # and the recorded four-arm table.
    sup = salty & ((tC - Tf(S)) < -1e-9)
    sfc = np.zeros_like(sup); sfc[im-1] = True
    return dict(
        salty        = int(salty.sum()),
        sup_total    = int(sup.sum()),
        sup_surface  = int((sup &  sfc).sum()),
        sup_interior = int((sup & ~sfc).sum()),
        minT         = float(tC[salty].min()),
        worst_defC   = float((tC - Tf(S))[salty].min()),
        meanT_fluid  = float(tC[fluid].mean()),
        meanT_salty  = float(tC[salty].mean()),
        rms_u        = float(np.sqrt((d["u"][fluid]**2).mean())),
        rms_v        = float(np.sqrt((d["v"][fluid]**2).mean())),
        rms_w        = float(np.sqrt((d["w"][fluid]**2).mean())),
    )

if __name__ == "__main__":
    for p in sys.argv[1:]:
        d = census(load(p))
        tag = os.path.basename(os.path.dirname(p))
        print(f"{tag:10s} salty={d['salty']:7d}  sup tot={d['sup_total']:7d} "
              f"(sfc {d['sup_surface']:5d} / int {d['sup_interior']:7d})  "
              f"minT={d['minT']:8.4f}  worst={d['worst_defC']:7.4f}  "
              f"meanT_fluid={d['meanT_fluid']:9.5f}  meanT_salty={d['meanT_salty']:9.5f}  "
              f"rms u/v/w={d['rms_u']:.6e} {d['rms_v']:.6e} {d['rms_w']:.6e}")
