import sys,numpy as np
# STORM-SHAPE, READ-ONLY (2026-10-06): the northern ocean rows at 54-62 deg rain 1590 / 1531 mm/a against NASA 1100 / 1048 while the southern rows at
# the same latitude are on NASA. Which basin, which water, and what separates them?   usage: nrows.py <run> [iter]
run=sys.argv[1]; it=sys.argv[2] if len(sys.argv)>2 else '60'
def rd(f,want):
    L=open(f,errors='replace').read().split('\n'); o={}
    for n,l in enumerate(L):
        if l.startswith('SCALARS '):
            nm=l.split()[1]
            if nm in want and nm not in o: o[nm]=np.array(L[n+2:n+2+65341],float).reshape(181,361)
    return o
s=rd(f'output_{run}/0Ma_smooth_Atm_radial_0_{it}.vtk',('Topography','Precipitation','Precipitation_NASA','Temperature','HumidityRel','PrecipitableWater','PressureStatic'))
lat=90-np.arange(181); lon=np.where(np.arange(361)<=180,np.arange(361),np.arange(361)-360); LA=lat[:,None]*np.ones((1,361)); LO=lon[None,:]*np.ones((181,1))
w=np.cos(np.radians(LA)); w[:,360]=0; oc=s['Topography']==0; land=~oc
P,N,T,RH,PW=(s[x] for x in('Precipitation','Precipitation_NASA','Temperature','HumidityRel','PrecipitableWater'))
g=lambda A,m:(A*w)[m].sum()/w[m].sum() if w[m].sum()>0 else np.nan
# distance to land along the latitude circle to the WEST (upwind in the westerlies), in degrees of longitude
west=np.zeros((181,361))
for j in range(181):
    d=999
    for k in list(range(360))*2:
        d=0 if land[j,k] else d+1
        west[j,k]=d
west[:,360]=west[:,0]
B={'N Atlantic 60W-10E':(-60,10),'N Pacific 140E-125W':None,'S Indian 20E-120E':(20,120),'S Pacific 150E-75W':None,'S Atlantic 60W-20E':(-60,20)}
def bm(nm,sg):
    if nm.startswith('N Pacific'): return (LO>=140)|(LO<=-125)
    if nm.startswith('S Pacific'): return (LO>=150)|(LO<=-75)
    a,b=B[nm]; return (LO>=a)&(LO<=b)
print(f'{run} {it}: ocean by basin and row.  row | basin | area share of the row % | T C | PW mm | RH % | model mm/a | NASA | model / NASA')
for a in(46,50,54,58,62):
    for nm in B:
        sg=1 if nm[0]=='N' else -1; m=oc&(sg*LA>=a)&(sg*LA<a+4)&bm(nm,sg); row=oc&(sg*LA>=a)&(sg*LA<a+4)
        if w[m].sum()/w[row].sum()>0.05: print(f'  {a}-{a+4}{nm[0]} | {nm:22s} | {100*w[m].sum()/w[row].sum():3.0f} | {g(T,m):5.1f} | {g(PW,m):4.1f} | {g(RH,m):4.1f} | {365*g(P,m):5.0f} | {365*g(N,m):5.0f} | {g(P,m)/g(N,m):4.2f}')
print('\nocean 50-66 deg, both hemispheres, by surface temperature:   T class | S: area %, model, NASA, ratio | N: area %, model, NASA, ratio')
for lo,hi in[(-6,-2),(-2,1),(1,3),(3,5),(5,7),(7,9),(9,14)]:
    out=[]
    for sg in(-1,1):
        tot=oc&(sg*LA>=50)&(sg*LA<66); m=tot&(T>=lo)&(T<hi)
        out.append(f'{100*w[m].sum()/w[tot].sum():4.1f}, {365*g(P,m):5.0f}, {365*g(N,m):5.0f}, {g(P,m)/g(N,m):4.2f}' if m.sum()>10 else '  --')
    print(f'  {lo:3d}..{hi:2d} C | '+' | '.join(out))
print('\nnorthern ocean 50-66N by the distance to land upwind (to the west along the row):   deg of longitude | area % | T C | model | NASA | ratio')
tot=oc&(LA>=50)&(LA<66)
for lo,hi in[(1,5),(5,10),(10,20),(20,40),(40,400)]:
    m=tot&(west>=lo)&(west<hi)
    if m.sum()>10: print(f'  {lo:3d}-{hi:3d} | {100*w[m].sum()/w[tot].sum():4.1f} | {g(T,m):5.1f} | {365*g(P,m):5.0f} | {365*g(N,m):5.0f} | {g(P,m)/g(N,m):4.2f}')
tot=oc&(LA<=-50)&(LA>-66)
print('the same for the southern ocean 50-66S:')
for lo,hi in[(1,5),(5,10),(10,20),(20,40),(40,400)]:
    m=tot&(west>=lo)&(west<hi)
    if m.sum()>10: print(f'  {lo:3d}-{hi:3d} | {100*w[m].sum()/w[tot].sum():4.1f} | {g(T,m):5.1f} | {365*g(P,m):5.0f} | {365*g(N,m):5.0f} | {g(P,m)/g(N,m):4.2f}')
