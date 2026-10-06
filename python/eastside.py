import sys,heapq,glob,re,numpy as np
# SUBTROP, READ-ONLY (2026-10-06): what separates the east sides that rain too little (S China, SE US, S Brazil) from the one that rains too much
# (E Australia) once the cap ATM_RH_LAND_EAST_MAX is on?  Rebuilds the initial land surface RH (as eastfetch.py) and splits every region into the
# cells AT the cap and the cells below it.   usage: eastside.py <run> [strength] [cap]
run=sys.argv[1]; S=float(sys.argv[2]) if len(sys.argv)>2 else 1.6; CAP=float(sys.argv[3]) if len(sys.argv)>3 else 0.010
RH_OCEAN,RH_DRY,L_C=0.82,0.35,1.5e6
def rd(f,want):
    L=open(f,errors='replace').read().split('\n'); o={}
    for n,l in enumerate(L):
        if l.startswith('SCALARS '):
            nm=l.split()[1]
            if nm in want and nm not in o: o[nm]=np.array(L[n+2:n+2+65341],float).reshape(181,361)
    return o
f=sorted(glob.glob(f'output_{run}/0Ma_smooth_Atm_radial_0_*.vtk'),key=lambda s:int(re.findall(r'_(\d+)\.vtk',s)[0]))[-1]
o=rd(f,('Topography','Landscape','Precipitation','PrecipitationConv','Precipitation_NASA','Temperature','Temperature_NASA','HumidityRel','PrecipitableWater'))
land=o['Topography']>0; jm,km=181,361; a=6.371e6; dlon=2*np.pi/(km-1); dlat=np.pi/(jm-1); Le=2.0e6
lat=90-np.arange(181); lon=np.where(np.arange(361)<=180,np.arange(361),np.arange(361)-360); w=np.cos(np.radians(lat))[:,None]*np.ones((1,361))
me=np.zeros((jm,km))
for j in range(1,jm-1):
    dx=a*dlon*np.cos(np.radians(lat[j])); n=min(km-2,int(3*Le/dx)); q=np.arange(1,n+1); wq=np.exp(-q*dx/Le)
    for k in range(km):
        kk=(k%(km-1)+q)%(km-1); me[j,k]=(wq*(~land[j,kk])).sum()/wq.sum()
d=np.full((jm,km),1e30); pq=[]
for j in range(jm):
    for k in range(km):
        if not land[j,k]: d[j,k]=0.; pq.append((0.,j,k))
heapq.heapify(pq)
while pq:
    dd,j,k=heapq.heappop(pq)
    if dd>d[j,k]: continue
    for dj in(-1,0,1):
        for dk in(-1,0,1):
            if not dj and not dk: continue
            jj=j+dj
            if jj<0 or jj>=jm: continue
            kk=((k+dk)%(km-1)+(km-1))%(km-1); cl=np.cos(np.radians(90-0.5*(j+jj))); st=a*np.hypot(dj*dlat,dk*dlon*cl)
            if dd+st<d[jj,kk]: d[jj,kk]=dd+st; heapq.heappush(pq,(dd+st,jj,kk))
d[:,km-1]=d[:,0]
al=np.abs(lat)[:,None]*np.ones((1,km)); sdesc=np.exp(-((al-25)/10)**2); c=1-np.exp(-d/L_C)
wet=0.75+np.clip((45-al)/15,0,1)*(RH_OCEAN-0.75)
r=np.minimum(wet-(wet-RH_DRY)*sdesc*(0.5+0.5*c)*(1-S*me),wet+CAP); atcap=r>=wet+CAP-1e-9
P,C,N,T,TN,H,RHn,PW=(o[x] for x in('Precipitation','PrecipitationConv','Precipitation_NASA','Temperature','Temperature_NASA','Landscape','HumidityRel','PrecipitableWater'))
reg={'E Australia':(-35,-18,145,153),'Madagascar':(-25,-13,43,50),'SE Africa':(-28,-15,28,38),'S Brazil':(-30,-18,-55,-42),'S China':(20,32,102,120),'SE US':(27,36,-95,-78),'India':(10,28,72,88)}
print(f'{run} ({f.split("_")[-1][:-4]}), strength {S}, cap +{CAP}.  Per region: |lat| | ground m | coast km | m_east | T sfc C | init RH (wet end) | RH now % | PW mm | cells at cap % | rain model (conv+strat) / NASA mm/d')
for nm,(a_,b_,c_,d_) in reg.items():
    m=land&((lat>=a_)&(lat<=b_))[:,None]&((lon>=c_)&(lon<=d_))[None,:]; g=lambda A,mm=None:(A*w)[m if mm is None else mm].sum()/w[m if mm is None else mm].sum()
    print(f'  {nm:12s} {g(al):4.1f} | {g(H):5.0f} | {g(d)/1e3:5.0f} | {g(me):.2f} | {g(T):5.1f} | {g(r):.3f} ({g(wet):.3f}) | {g(RHn):4.1f} | {g(PW):4.1f} | {100*w[m&atcap].sum()/w[m].sum():3.0f} | {g(P):5.2f} ({g(C):4.2f}+{g(P-C):4.2f}) / {g(N):4.2f}')
    for lab,mm in(('at cap',m&atcap),('below',m&~atcap)):
        if mm.sum()>3: print(f'      {lab:6s} n={mm.sum():3d}: |lat| {g(al,mm):4.1f}  ground {g(H,mm):5.0f} m  coast {g(d,mm)/1e3:4.0f} km  T {g(T,mm):5.1f}  init RH {g(r,mm):.3f}  RH now {g(RHn,mm):4.1f}  rain {g(P,mm):5.2f} ({g(C,mm):4.2f}+{g(P-C,mm):4.2f}) / NASA {g(N,mm):4.2f}')
m=land&(al>=15)&(al<35)
print('\nall land 15-35 deg, cells AT the cap, by surface temperature:   T class | n | ground m | rain (conv+strat) | NASA')
for lo,hi in[(-20,10),(10,15),(15,20),(20,24),(24,27),(27,40)]:
    mm=m&atcap&(T>=lo)&(T<hi)
    if mm.sum()>3: print(f'   {lo:3d}..{hi:2d} C | {mm.sum():4d} | {(H*w)[mm].sum()/w[mm].sum():5.0f} | {(P*w)[mm].sum()/w[mm].sum():5.2f} ({(C*w)[mm].sum()/w[mm].sum():4.2f}+{((P-C)*w)[mm].sum()/w[mm].sum():4.2f}) | {(N*w)[mm].sum()/w[mm].sum():4.2f}')
print('all land 15-35 deg, cells AT the cap, by |latitude| (the wet end falls from 0.82 at 30 deg to 0.75 at 45):')
for lo,hi in[(15,20),(20,25),(25,30),(30,35)]:
    mm=m&atcap&(al>=lo)&(al<hi)
    if mm.sum()>3: print(f'   {lo}-{hi} | {mm.sum():4d} | ground {(H*w)[mm].sum()/w[mm].sum():5.0f} | T {(T*w)[mm].sum()/w[mm].sum():5.1f} | init RH {(r*w)[mm].sum()/w[mm].sum():.3f} | {(P*w)[mm].sum()/w[mm].sum():5.2f} ({(C*w)[mm].sum()/w[mm].sum():4.2f}+{((P-C)*w)[mm].sum()/w[mm].sum():4.2f}) | NASA {(N*w)[mm].sum()/w[mm].sum():4.2f}')
print(f'land 15-35: at cap {100*w[m&atcap].sum()/w[m].sum():.0f} % of the area, rain {(P*w)[m&atcap].sum()/w[m&atcap].sum():.2f} / NASA {(N*w)[m&atcap].sum()/w[m&atcap].sum():.2f};  below cap rain {(P*w)[m&~atcap].sum()/w[m&~atcap].sum():.2f} / NASA {(N*w)[m&~atcap].sum()/w[m&~atcap].sum():.2f}')
print('below the cap, by initial RH:   RH class | area % of land 15-35 | rain | NASA')
for lo,hi in[(0.3,0.5),(0.5,0.6),(0.6,0.7),(0.7,0.75),(0.75,0.80),(0.80,0.83)]:
    mm=m&~atcap&(r>=lo)&(r<hi)
    if mm.sum()>3: print(f'   {lo:.2f}-{hi:.2f} | {100*w[mm].sum()/w[m].sum():4.1f} | {(P*w)[mm].sum()/w[mm].sum():5.2f} | {(N*w)[mm].sum()/w[mm].sum():4.2f}')
