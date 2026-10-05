import sys,numpy as np
def rd(f):
    L=open(f,errors='replace').read().split('\n'); o={}
    for n,l in enumerate(L):
        if l.startswith('SCALARS '):
            nm=l.split()[1]
            if nm in('Precipitation','PrecipitationConv','Precipitation_NASA','Topography','Landscape') and nm not in o: o[nm]=np.array(L[n+2:n+2+65341],float).reshape(181,361)
    return o
lat=90-np.arange(181); lon=np.where(np.arange(361)<=180,np.arange(361),np.arange(361)-360); w=np.cos(np.radians(lat))[:,None]*np.ones((1,361))
R={n:rd(f'output_{n}/0Ma_smooth_Atm_radial_0_220.vtk') for n in sys.argv[1:]}
T=R[sys.argv[1]]['Topography']; land=T>0; N=R[sys.argv[1]]['Precipitation_NASA']
print('land rain by latitude, mm/d: total (conv) per run | NASA')
for c in range(-22,23,4):
    m=land&((lat>=c-2)&(lat<c+2))[:,None]; g=lambda A:(A*w)[m].sum()/w[m].sum()
    print(f'{c:4d} '+' '.join(f"{g(R[n]['Precipitation']):5.2f} ({g(R[n]['PrecipitationConv']):5.2f})" for n in R)+f' | {g(N):5.2f}')
reg={'Amazon':(-10,2,-70,-52),'Congo':(-5,5,15,28),'Maritime':(-8,6,95,150),'C America':(8,18,-95,-78),'W Africa':(4,12,-12,12),'E Africa':(-10,8,30,42),'N S-America':(2,10,-75,-55),'SE Asia':(10,20,95,108),'N Austral':(-18,-11,125,145),'India':(10,28,72,88)}
H=R[sys.argv[1]]['Landscape']   # ground elevation, m (Topography is only a land flag)
print('regions mm/d: total (conv) per run | NASA | land elevation mean m')
for nm,(a,b,c,d) in reg.items():
    m=land&((lat>=a)&(lat<=b))[:,None]&((lon>=c)&(lon<=d))[None,:]; g=lambda A:(A*w)[m].sum()/w[m].sum()
    print(f'{nm:12s} '+' '.join(f"{g(R[n]['Precipitation']):5.2f} ({g(R[n]['PrecipitationConv']):5.2f})" for n in R)+f' | {g(N):5.2f} | {H[m].mean():5.0f}')
# tropical land (< 15 deg) by ground elevation: stratiform | convective | NASA, mm/d
m=land&(abs(lat)<15)[:,None]
for n in R:
    S=R[n]['Precipitation']-R[n]['PrecipitationConv']; C=R[n]['PrecipitationConv']
    for lo,hi in[(-1,150),(150,400),(400,800),(800,1300),(1300,2000),(2000,9000)]:
        mm=m&(H>lo)&(H<=hi); g=lambda A:(A*w)[mm].sum()/w[mm].sum()
        if mm.any(): print(f'{n} elev {lo:5d}-{hi:5d} m: area {w[mm].sum()/w[m].sum():.2f}  strat {g(S):5.2f}  conv {g(C):5.2f}  NASA {g(N):5.2f}')
