#!/usr/bin/env python3
"""Read an ATM_BUOY_TREF arm pair: the u trajectory, the rhs_u balance, the cells, the climate."""
import re, sys, glob, numpy as np

def logvals(f, pat, cast=float):
    try: txt=open(f, errors='replace').read()
    except OSError: return []
    return [cast(x) for x in re.findall(pat, txt)]

def maxu(f):   return logvals(f, r"max u-component \.+ =\s+(-?[\d.]+)")
def minu(f):   return logvals(f, r"min u-component \.+ =\s+(-?[\d.]+)")
def bal(f):
    t=open(f, errors='replace').read()
    return dict(
      buoy = [float(x) for x in re.findall(r"rms buoy = ([\d.eE+-]+)", t)],
      pgf  = [float(x) for x in re.findall(r"rms pgf = ([\d.eE+-]+)", t)],
      corr = [float(x) for x in re.findall(r"corr\(pgf,buoy\) = (-?[\d.]+)", t)],
      canc = [float(x) for x in re.findall(r"cancellation = (-?[\d.]+)", t)],
      net  = [float(x) for x in re.findall(r"rms NET rhs_u \(all six terms\) = ([\d.eE+-]+)", t)],
    )
def clim(f):
    t=open(f, errors='replace').read()
    def g(p, d=float('nan')):
        m=re.findall(p,t)
        return float(m[-1]) if m else d
    return dict(precip=g(r"model\s+([\d.]+)\s+NASA"),
                r=g(r"pattern r = ([+-][\d.]+)"),
                rms=g(r"centred RMS\s+([\d.]+)"),
                sig=g(r"sigma model/NASA = ([\d.]+)"),
                b35=g(r"35-65\s+([\d.]+) /"),
                b65=g(r"65-90\s+([\d.]+) /"),
                wu=g(r"max w_u \.+ =\s+([\d.]+)"))
def psi_cells(d, it):
    fs=sorted(glob.glob(f"{d}/meridional_streamfunction_{it}.csv"))
    if not fs: return None
    import csv
    rows=list(csv.DictReader(open(fs[0])))
    out={}
    for lat in (75,45,15,-15,-45,-75):
        col=[r for r in rows if abs(float(r['lat_deg'])-lat)<0.5]
        if not col: continue
        p=np.array([float(r['psi_fixdiv_kg_per_s']) for r in col])
        gnd=p[np.argmin([float(r['height_m']) for r in col])]
        out[lat]=(np.abs(p-gnd).max()/1e9, abs(gnd)/max(np.abs(p-gnd).max(),1e-30))
    return out

def show(tag_ctl, tag_on, label):
    print(f"\n{'='*78}\n{label}\n{'='*78}")
    for nm,tag in (("control",tag_ctl),("BUOY_TREF=1",tag_on)):
        f=f"{tag}.log"; mu=maxu(f); mn=minu(f)
        if not mu: print(f"  {nm}: no data yet"); continue
        b=bal(f); c=clim(f)
        print(f"  {nm:12s} n_iter={len(mu):4d}  max|u| first {max(abs(mu[0]),abs(mn[0])):.6f}"
              f"  last {max(abs(mu[-1]),abs(mn[-1])):.6f}"
              f"  peak {max(max(map(abs,mu)),max(map(abs,mn))):.6f}")
        if b['buoy']:
            print(f"               rms buoy {b['buoy'][-1]:.4e}   rms pgf {b['pgf'][-1]:.4e}"
                  f"   corr {b['corr'][-1]:+.4f}   canc {b['canc'][-1]:.4f}   NET {b['net'][-1]:.4e}")
        print(f"               Precip {c['precip']:.1f}  r {c['r']:+.3f}  cRMS {c['rms']:.1f}"
              f"  sigma {c['sig']:.2f}  35-65 {c['b35']:.1f}  65-90 {c['b65']:.1f}  max w_u {c['wu']:.6f}")
    # u trajectory side by side
    a,bb=maxu(f"{tag_ctl}.log"), maxu(f"{tag_on}.log")
    an,bn=minu(f"{tag_ctl}.log"), minu(f"{tag_on}.log")
    if a and bb:
        n=min(len(a),len(bb)); step=max(1,n//12)
        print(f"\n  {'iter':>6} {'ctl max|u|':>12} {'tref max|u|':>12} {'ratio':>8}")
        for k in range(0,n,step):
            A=max(abs(a[k]),abs(an[k])); B=max(abs(bb[k]),abs(bn[k]))
            print(f"  {k:>6} {A:>12.6f} {B:>12.6f} {B/A if A else float('nan'):>8.4f}")
for pair,label in ((("bt_ctl","bt_on"),"RESTART PAIR 600->700, buoyancy_ramp = 1.0 (FULL STRENGTH)"),
                   (("bs_ctl","bs_on"),"FROM-SCRATCH PAIR 0->100, ramp 0->0.33 (avg ~0.17)")):
    show(*pair, label)
