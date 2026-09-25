#!/usr/bin/env python3
# Score run_sgzb600.sh: sgzb_ctl / _e / _sg / _sge. Climate from the logs, parcel from the zonal-87 slice.
import re, sys, os, numpy as np
SLICE=int(os.environ.get('SLICE','520'))   # zonal slices exist at 20, 120, ..., 520 (VTK stride 5)
arms = sys.argv[1:] or ['sgzb_ctl','sgzb_e','sgzb_sg','sgzb_sge']
def logscan(a):
    P=[];R=[];S=[];B=[];LO=[];gp=None;cap=None;mx={};nan=0
    for l in open(a+'.log', errors='replace'):
        m=re.search(r'Precip mean = ([\d.e+-]+)',l);  P.append(float(m.group(1))) if m else None
        m=re.search(r'pattern r = ([+-][\d.]+).*sigma model/NASA = ([\d.]+)',l)
        if m: R.append(float(m.group(1))); S.append(float(m.group(2)))
        if 'by |latitude|' in l: B.append('/'.join(x.split('.')[0] for x in re.findall(r'(\d+\.\d) /',l)))
        m=re.search(r'land\s+([\d.]+) /\s+[\d.]+\s+ocean\s+([\d.]+)',l); LO.append((m.group(1),m.group(2))) if m else None
        m=re.search(r'g_p generated ([\d.e+-]+).*P_conv\(ground\) ([\d.e+-]+)',l)
        if m: gp=(round(float(m.group(1)),1),round(float(m.group(2)),3))
        m=re.search(r'MC_t cap=\S+\s+truncated \d+ \(([\d.]+) %\)',l); cap=float(m.group(1)) if m else cap
        m=re.search(r'max ([uvw])-component\s+=\s+([-\d.]+)',l)
        if m: mx[m.group(1)]=float(m.group(2))
        if 'NaN/Inf DETECTED' in l: nan+=1
    return P,R,S,B,LO,gp,cap,mx,nan
print('== CLIMATE (high/low parity alternate; iteration index = print index)')
print(f'{"arm":9s} ' + ' '.join(f'P{i:>4d}' for i in (100,200,300,400,500,600)) + '  slope400-600  r600   sig  bands@600            land/ocean   g_p/P_conv      MC_t trunc%  max u/v/w   NaN')
for a in arms:
    try: P,R,S,B,LO,gp,cap,mx,nan = logscan(a)
    except FileNotFoundError: print(a,'no log'); continue
    n=len(P)-1; g=lambda i: f'{P[i]:5.0f}' if i<=n else '    -'
    sl=(P[min(n,600)]-P[400])/(min(n,600)-400) if n>400 else float('nan')
    i=min(n,600)
    print(f'{a:9s} '+' '.join(g(k) for k in (100,200,300,400,500,600))+f'  {sl:11.2f}  {R[i]:+.3f} {S[i]:5.2f}  {B[i]:20s} {"/".join(LO[-1]) if LO else "":12s} {str(gp):15s} {cap}  {mx.get("u",0):.3f}/{mx.get("v",0):.2f}/{mx.get("w",0):.2f}  {nan}')
print()
print(f'== PARCEL, zonal 87E slice at iteration {SLICE}, cells where the updraft recurrence runs (|M_u| > 0.1, file x1e3)')
N=7421; t0=273.15
def slice_(f):
    L=open(f).read().split('\n')
    pts=np.array([list(map(float,L[6+n].split())) for n in range(N)])
    def get(name):
        for n,l in enumerate(L):
            if l.startswith('SCALARS '+name+' '): return np.array(L[n+2:n+2+N],float)
    return pts,{k:get(k) for k in ['height','s','s_u','M_u','D_u','Temperature','g_p']}
for a in arms:
    f=f'output_{a}/0Ma_smooth_Atm_zonal_87_{SLICE}.vtk'
    try: pts,d=slice_(f)
    except FileNotFoundError: print(a,'no slice'); continue
    xs=np.unique(np.round(pts[:,0],6)); lev=np.searchsorted(xs,np.round(pts[:,0],6)); col=np.round(pts[:,1],6)
    act=np.abs(d['M_u'])>100.0
    dT=(d['s_u']-d['s'])*t0
    # layer thickness per level from the height field of one column
    h=np.array([np.median(d['height'][lev==i]) for i in range(41)])*1000; dz=np.append(np.diff(h),np.nan)
    mixf=np.where(act, np.abs(d['D_u'])*dz[lev]/np.maximum(np.abs(d['M_u']),1e-12), np.nan)
    print(f'{a}: active cells {act.sum()} in {len(np.unique(col[act]))} columns')
    for lo,hi in ((0,2),(2,5),(5,8),(8,12),(12,17)):
        m=act&(d['height']>=lo)&(d['height']<hi)
        if m.sum()<3: print(f'   {lo:2d}-{hi:2d} km   n={m.sum()}'); continue
        print(f'   {lo:2d}-{hi:2d} km   n={m.sum():4d}   s_u-s median {np.median(dT[m]):+7.2f} K  (p10 {np.percentile(dT[m],10):+7.2f}, p90 {np.percentile(dT[m],90):+7.2f})   mixing D_u*dz/M_u median {np.nanmedian(mixf[m]):.3f}')
