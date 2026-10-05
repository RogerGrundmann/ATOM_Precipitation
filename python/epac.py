import sys,numpy as np
from scipy.ndimage import gaussian_filter
# RC-RESID, READ-ONLY (2026-10-05): why is the E Pacific rain band 28 deg wide against NASA's 9?   usage: epac.py <run> [iteration]
run=sys.argv[1]; it=sys.argv[2] if len(sys.argv)>2 else '520'
W=('Temperature','Topography','Precipitation','PrecipitationConv','Precipitation_NASA','Evaporation','HumidityRel','WaterVapour','PressureStatic','w-Component','v-Component')
def rd(f):
    L=open(f,errors='replace').read().split('\n'); o={}
    for n,l in enumerate(L):
        if l.startswith('SCALARS '):
            nm=l.split()[1]
            if nm in W and nm not in o: o[nm]=np.array(L[n+2:n+2+65341],float).reshape(181,361)
    return o
o=rd(f'output_{run}/0Ma_smooth_Atm_radial_0_{it}.vtk'); o0=rd(f'output_{run}/0Ma_smooth_Atm_radial_0_0.vtk')
lat=90-np.arange(181); lon=np.where(np.arange(361)<=180,np.arange(361),np.arange(361)-360); w=np.cos(np.radians(lat))[:,None]*np.ones((1,361))
P,C,N=o['Precipitation'],o['PrecipitationConv'],o['Precipitation_NASA']; oc=o['Topography']==0; T0=o0['Temperature']
T=o['Temperature']; es=6.1078*np.exp(17.2694*T/(T+273.15-35.86)); rh0=o['WaterVapour']/(0.622*es/(o['PressureStatic']-es)*1000)   # surface RH at the iteration read (the level-0 value of the iteration-0 slice is the skin before the profile acts)
print(f'== {run} iteration {it}.  E Pacific 150W-100W, ocean, by latitude: SST C | surface RH | model rain (conv + strat) | NASA | E mm/d')
k=(lon>=-150)&(lon<=-100)
for c in range(20,-22,-2):
    j=90-c; m=k&oc[j]
    print(f'  {c:4d}: {T0[j,m].mean():5.1f} | {rh0[j,m].mean():.3f} | {P[j,m].mean():4.1f} ({C[j,m].mean():4.1f} + {(P-C)[j,m].mean():4.1f}) | {N[j,m].mean():4.1f} | {o["Evaporation"][j,m].mean():4.2f}')
# the same SST class rains differently in different basins: NASA and model by SST class in four sectors (|lat|<=25 ocean)
sec={'E Pacific S of eq (150W-90W, 0-25S)':((-25,0),(-150,-90)),'E Pacific N of eq (150W-90W, 0-25N)':((0,25),(-150,-90)),'W Pacific (130E-180, 25S-25N)':((-25,25),(130,180)),'Atlantic S of eq (35W-10E, 0-25S)':((-25,0),(-35,10)),'Indian (50-100E, 25S-25N)':((-25,25),(50,100))}
print('\n  rain by SST class per sector: NASA / model mm/d (area %)')
print('  sector                                  ' + ' '.join(f'{a}-{b} C        ' for a,b in [(22,25),(25,26),(26,27),(27,28),(28,31)]))
for nm,((a,b),(c,d)) in sec.items():
    m=oc&((lat>=a)&(lat<=b))[:,None]&((lon>=c)&(lon<=d))[None,:]; row=[]
    for lo,hi in [(22,25),(25,26),(26,27),(27,28),(28,31)]:
        mm=m&(T0>=lo)&(T0<hi)
        row.append(f'{(N*w)[mm].sum()/w[mm].sum():4.1f}/{(P*w)[mm].sum()/w[mm].sum():4.1f} ({100*w[mm].sum()/w[m].sum():2.0f})' if mm.sum()>20 else '      -      ')
    print(f'  {nm:38s}  '+'  '.join(row))
# what predicts NASA rain over the tropical ocean besides SST itself?
m=oc&(abs(lat)<=25)[:,None]
def sm(A,s):                                             # ocean-only Gaussian mean, sigma in degrees
    a=gaussian_filter(np.where(oc,A,0.),s,mode=('nearest','wrap')); b=gaussian_filter(oc.astype(float),s,mode=('nearest','wrap')); return a/np.maximum(b,1e-3)
Tz=np.array([T0[j][oc[j]].mean() if oc[j].any() else np.nan for j in range(181)])[:,None]*np.ones((1,361))
Ttrop=(T0*w)[m].sum()/w[m].sum()
pred={'SST':T0,'SST - zonal ocean mean':T0-Tz,'SST - 10 deg neighbourhood mean':T0-sm(T0,10),'SST - 20 deg neighbourhood mean':T0-sm(T0,20),
      'SST - warmest SST within 10 deg':T0-np.maximum.reduce([np.roll(np.roll(np.where(oc,T0,-99),dj,0),dk,1) for dj in range(-10,11,2) for dk in range(-10,11,2)]),
      '-(Laplacian of 5-deg-smoothed SST)':-(np.gradient(np.gradient(sm(T0,5),axis=0),axis=0)+np.gradient(np.gradient(sm(T0,5),axis=1),axis=1))}
print(f'\n  |lat|<=25 ocean: correlation with NASA rain | with the model rain | NASA rain explained by SST + this (multiple r)')
x0=T0[m]; y=N[m]; ym=P[m]
for nm,A in pred.items():
    x=A[m]; X=np.c_[np.ones(len(x0)),x0,x]; b=np.linalg.lstsq(X,y,rcond=None)[0]; R=np.corrcoef(X@b,y)[0,1]
    print(f'    {nm:36s} {np.corrcoef(x,y)[0,1]:+.3f} | {np.corrcoef(x,ym)[0,1]:+.3f} | {R:.3f}')
print(f'    model rain itself vs NASA: {np.corrcoef(ym,y)[0,1]:+.3f}')
# E Pacific south of the equator: how far below its neighbourhood is the SST where NASA is dry?
k2=oc&((lat>=-20)&(lat<=2))[:,None]&((lon>=-150)&(lon<=-85))[None,:]
d10=(T0-sm(T0,10)); dz=T0-Tz
print(f'\n  E Pacific 20S-2N, 150W-85W (NASA {N[k2].mean():.2f}, model {P[k2].mean():.2f} mm/d): SST {T0[k2].mean():.1f} C, minus zonal mean {dz[k2].mean():+.2f} K, minus 10-deg neighbourhood {d10[k2].mean():+.2f} K, share of area colder than 25.5 C {100*np.mean(T0[k2]<25.5):.0f} %, 25.5-26.5 C {100*np.mean((T0[k2]>=25.5)&(T0[k2]<26.5)):.0f} %, warmer {100*np.mean(T0[k2]>=26.5):.0f} %')
for lo,hi in [(10,24),(24,25.5),(25.5,26.5),(26.5,27.5),(27.5,31)]:
    mm=k2&(T0>=lo)&(T0<hi)
    if mm.sum()>10: print(f'     SST {lo}-{hi} C: area {100*mm.sum()/k2.sum():3.0f} % | NASA {N[mm].mean():4.2f} | model {P[mm].mean():4.2f} (conv {C[mm].mean():4.2f}, strat {(P-C)[mm].mean():4.2f}) | surface RH {rh0[mm].mean():.3f}')
