import sys,heapq,numpy as np
# SUBTROP, READ-ONLY (2026-10-05): the east-fetch factor of ATM_RH_LAND_EAST and the initial land surface RH it gives, rebuilt from
# InitValues_Atm.cpp (m_east, rh_land_of).  usage: eastfetch.py <run> [strength]     Why does E Australia rain 3x NASA?
run=sys.argv[1]; S=float(sys.argv[2]) if len(sys.argv)>2 else 1.35
RH_OCEAN,RH_DRY,L_C=0.82,0.35,1.5e6
def rd(f,want):
    L=open(f,errors='replace').read().split('\n'); o={}
    for n,l in enumerate(L):
        if l.startswith('SCALARS '):
            nm=l.split()[1]
            if nm in want and nm not in o: o[nm]=np.array(L[n+2:n+2+65341],float).reshape(181,361)
    return o
import glob,re
f=sorted(glob.glob(f'output_{run}/0Ma_smooth_Atm_radial_0_*.vtk'),key=lambda s:int(re.findall(r'_(\d+)\.vtk',s)[0]))[-1]
o=rd(f,('Topography','Landscape','Precipitation','PrecipitationConv','Precipitation_NASA'))
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
def rhs(strength,clamp=False):
    east=1-strength*me
    if clamp: east=np.maximum(east,0.)
    return wet-(wet-RH_DRY)*sdesc*(0.5+0.5*c)*east, east
r,east=rhs(S); rc,_=rhs(S,True); r0,_=rhs(0.)
P,C,N=o['Precipitation'],o['PrecipitationConv'],o['Precipitation_NASA']
reg={'E Australia':(-35,-18,145,153),'Madagascar':(-25,-13,43,50),'E Mexico coast':(18,26,-100,-96),'S Brazil':(-30,-18,-55,-42),'SE Africa':(-28,-15,28,38),'S China':(20,32,102,120),'SE US':(27,36,-95,-78),'India':(10,28,72,88),'Arabia':(16,28,42,56),'Sahara':(18,30,-10,25),'W Australia':(-30,-20,115,128)}
print(f'{run} ({f.split("_")[-1][:-4]}), ATM_RH_LAND_EAST = {S}:  ocean fetch to the east m_east | factor 1 - s*m | cells with factor < 0 | initial surface RH: no fetch -> with fetch (clamped at the wet end {RH_OCEAN}) | max | rain model / NASA mm/d')
for nm,(a_,b_,c_,d_) in reg.items():
    m=land&((lat>=a_)&(lat<=b_))[:,None]&((lon>=c_)&(lon<=d_))[None,:]; g=lambda A:(A*w)[m].sum()/w[m].sum()
    print(f'  {nm:15s} {g(me):.2f} | {g(east):+.2f} | {100*np.mean(east[m]<0):3.0f} % | {g(r0):.3f} -> {g(r):.3f} ({g(rc):.3f}) | {r[m].max():.3f} | {g(P):5.2f} / {g(N):4.2f}')
m=land&((al>=15)&(al<35))
print(f'\nland 15-35 deg: cells with a NEGATIVE factor (RH above the wet end) {100*w[m&(east<0)].sum()/w[m].sum():.0f} % of the area; they carry {100*(P*w)[m&(east<0)].sum()/(P*w)[m].sum():.0f} % of the band\x27s rain ({(P*w)[m&(east<0)].sum()/w[m&(east<0)].sum():.2f} mm/d there, NASA {(N*w)[m&(east<0)].sum()/w[m&(east<0)].sum():.2f}); the rest {(P*w)[m&(east>=0)].sum()/w[m&(east>=0)].sum():.2f} mm/d, NASA {(N*w)[m&(east>=0)].sum()/w[m&(east>=0)].sum():.2f}')
print('rain against the initial surface RH, land 15-35 deg, ground below 500 m:   RH class | area % | model (conv + strat) | NASA')
H=o['Landscape']
for lo,hi in[(0.3,0.6),(0.6,0.7),(0.7,0.75),(0.75,0.80),(0.80,0.84),(0.84,0.88),(0.88,1.0)]:
    mm=m&(H<500)&(r>=lo)&(r<hi)
    if mm.sum()>5: print(f'   {lo:.2f}-{hi:.2f} | {100*w[mm].sum()/w[m&(H<500)].sum():4.1f} | {(P*w)[mm].sum()/w[mm].sum():5.2f} ({(C*w)[mm].sum()/w[mm].sum():4.2f} + {((P-C)*w)[mm].sum()/w[mm].sum():5.2f}) | {(N*w)[mm].sum()/w[mm].sum():4.2f}')
# ---- rough estimate: what a clamp at the wet end would do, using the measured rain-vs-RH response of this run (lowland classes above)
xs=np.array([0.60,0.72,0.775,0.82,0.86,0.90,0.95]); ys=np.array([0.0,0.02,0.45,2.46,6.68,10.48,12.0])
low=m&(H<500)
resp=lambda RH: np.interp(RH,xs,ys)
act=(P*w)[low].sum()/w[low].sum(); est=(resp(r)*w)[low].sum()/w[low].sum()
print(f'\nESTIMATE (response curve applied to the initial RH; lowland < 500 m at 15-35 deg, {100*w[low].sum()/w[m].sum():.0f} % of the band): actual {act:.2f} mm/d, curve gives {est:.2f} (NASA {(N*w)[low].sum()/w[low].sum():.2f})')
print('   strength, clamp | lowland mean mm/d | E Australia | S China | SE US | S Brazil | Madagascar | Arabia | area raining > 1 mm/d %')
def regm(A,nm):
    a_,b_,c_,d_=reg[nm]; mm=land&(H<500)&((lat>=a_)&(lat<=b_))[:,None]&((lon>=c_)&(lon<=d_))[None,:]; return (A*w)[mm].sum()/w[mm].sum() if mm.any() else np.nan
for s_,cl in[(1.35,False),(1.35,True),(1.6,True),(2.0,True),(3.0,True)]:
    rr,_=rhs(s_,cl); R=resp(rr)
    print(f'   {s_:4.2f} {"clamp" if cl else "     "} | {(R*w)[low].sum()/w[low].sum():5.2f} | '+' | '.join(f'{regm(R,n):5.2f}' for n in('E Australia','S China','SE US','S Brazil','Madagascar','Arabia'))+f' | {100*w[low&(R>1)].sum()/w[low].sum():4.1f}')
print('   NASA            | '+f'{(N*w)[low].sum()/w[low].sum():5.2f} | '+' | '.join(f'{regm(N,n):5.2f}' for n in('E Australia','S China','SE US','S Brazil','Madagascar','Arabia'))+f' | {100*w[low&(N>1)].sum()/w[low].sum():4.1f}')
