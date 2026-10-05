import sys,numpy as np
# tropical-ocean rain structure from radial VTK slices at 220: Pacific band width / peak by sector, rain binned by initial SST,
# and ocean surface RH / evaporation by latitude.  usage: pacband.py <run> [<run> ...]   (reads output_<run>/)
W=('Temperature','WaterVapour','PressureStatic','Evaporation','Topography','Precipitation','PrecipitationConv','Precipitation_NASA')
def rd(f):
    L=open(f,errors='replace').read().split('\n'); o={}
    for n,l in enumerate(L):
        if l.startswith('SCALARS '):
            nm=l.split()[1]
            if nm in W and nm not in o: o[nm]=np.array(L[n+2:n+2+65341],float).reshape(181,361)
    return o
lat=90-np.arange(181); lon=np.where(np.arange(361)<=180,np.arange(361),np.arange(361)-360); w=np.cos(np.radians(lat))[:,None]*np.ones((1,361))
for n in sys.argv[1:]:
    o=rd(f'output_{n}/0Ma_smooth_Atm_radial_0_220.vtk'); T0=rd(f'output_{n}/0Ma_smooth_Atm_radial_0_0.vtk')['Temperature']
    P,C,N=o['Precipitation'],o['PrecipitationConv'],o['Precipitation_NASA']; oc=o['Topography']==0
    print(f'== {n}')
    for nm,(a,b) in {'E Pacific 150W-100W':(-150,-100),'C Pacific 180-150W':(-180,-150),'W Pacific 150E-180':(150,180)}.items():
        k=(lon>=a)&(lon<=b); pm=np.array([P[j,k&oc[j]].mean() for j in range(60,121)]); nn=np.array([N[j,k&oc[j]].mean() for j in range(60,121)])
        eq=lambda A:A[25:36].mean()   # 5N..5S
        print(f'  {nm}: peak {pm.max():4.1f} at {90-60-pm.argmax():+3d} (NASA {nn.max():4.1f} at {90-60-nn.argmax():+3d}); width above half peak {(pm>pm.max()/2).sum():2d} deg (NASA {(nn>nn.max()/2).sum():2d}); 5S-5N {eq(pm):4.1f} (NASA {eq(nn):4.1f}); 30S-30N mean {pm.mean():4.2f} (NASA {nn.mean():4.2f})')
    m=oc&(abs(lat)<=30)[:,None]
    print('  |lat|<=30 ocean by initial SST: area % | NASA | model | conv   mm/d')
    for lo,hi in [(10,24),(24,25),(25,26),(26,27),(27,28),(28,29),(29,31)]:
        mm=m&(T0>=lo)&(T0<hi); g=lambda A:(A*w)[mm].sum()/w[mm].sum()
        print(f'    {lo:2d}-{hi:2d} C: {100*w[mm].sum()/w[m].sum():5.1f} | {g(N):5.2f} | {g(P):5.2f} | {g(C):5.2f}')
    print(f'  |lat|<=30 ocean: mean {(P*w)[m].sum()/w[m].sum():.2f} (NASA {(N*w)[m].sum()/w[m].sum():.2f}), r(model, NASA) {np.corrcoef(P[m],N[m])[0,1]:.3f}, r(model, SST) {np.corrcoef(P[m],T0[m])[0,1]:.3f}, rainless (< 0.3 mm/d) {100*w[m&(P<0.3)].sum()/w[m].sum():.0f} % (NASA {100*w[m&(N<0.3)].sum()/w[m].sum():.0f} %), p50/p90/p99 {np.percentile(P[m],50):.1f}/{np.percentile(P[m],90):.1f}/{np.percentile(P[m],99):.1f} (NASA {np.percentile(N[m],50):.1f}/{np.percentile(N[m],90):.1f}/{np.percentile(N[m],99):.1f})')
    T=o['Temperature'];es=6.1078*np.exp(17.2694*T/(T+273.15-35.86)); rh=o['WaterVapour']/(0.622*es/(o['PressureStatic']-es)*1000)
    print('  ocean by latitude  SST C | surface RH mean (min..max) | E mm/d:')
    for c in (30,24,18,12,6,0,-6,-12,-18,-24,-30):
        j=90-c; mm=oc[j]; print(f'    {c:4d}: {T[j,mm].mean():5.1f} | {rh[j,mm].mean():.3f} ({rh[j,mm].min():.3f}..{rh[j,mm].max():.3f}) | {o["Evaporation"][j,mm].mean():4.2f}')
