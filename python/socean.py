import sys,numpy as np
# RC-RESID, READ-ONLY (2026-10-06): why does the Southern Ocean poleward of 65 deg rain 184 mm/a against NASA 490 when the Arctic Ocean is at 328 / 379?
# Ocean rows of both hemispheres side by side.   usage: socean.py <run> [iter]
run=sys.argv[1]; it=sys.argv[2] if len(sys.argv)>2 else '520'
def rd(f,want):
    L=open(f,errors='replace').read().split('\n'); o={}
    for n,l in enumerate(L):
        if l.startswith('SCALARS '):
            nm=l.split()[1]
            if nm in want and nm not in o: o[nm]=np.array(L[n+2:n+2+65341],float).reshape(181,361)
    return o
s=rd(f'output_{run}/0Ma_smooth_Atm_radial_0_{it}.vtk',('Topography','Precipitation','PrecipitationSnow','Precipitation_NASA','Temperature','Temperature_NASA','HumidityRel','PressureStatic','PrecipitableWater','Evaporation'))
lat=90-np.arange(181); LA=lat[:,None]*np.ones((1,361)); w=np.cos(np.radians(LA)); w[:,360]=0; oc=s['Topography']==0
P,N,T,TN,RH,ps,PW,E,SN=(s[x] for x in('Precipitation','Precipitation_NASA','Temperature','Temperature_NASA','HumidityRel','PressureStatic','PrecipitableWater','Evaporation','PrecipitationSnow'))
g=lambda A,m:(A*w)[m].sum()/w[m].sum() if w[m].sum()>0 else np.nan
print(f'{run} {it}: OCEAN rows.  |lat| || SOUTH: ocean % of the row | T sfc C | PW mm | RH % | p hPa | rain mm/a | NASA | rain / PW per day || NORTH: the same')
for a in range(34,86,4):
    out=[]
    for sgn in(-1,1):
        row=(sgn*LA>=a)&(sgn*LA<a+4); m=row&oc
        if w[m].sum()/w[row].sum()<0.03: out.append('   --   (no ocean)'+' '*40); continue
        out.append(f'{100*w[m].sum()/w[row].sum():3.0f} | {g(T,m):6.1f} | {g(PW,m):4.1f} | {g(RH,m):4.1f} | {g(ps,m):4.0f} | {365*g(P,m):4.0f} | {365*g(N,m):4.0f} | {g(P,m)/g(PW,m):5.3f} (NASA {g(N,m)/g(PW,m):5.3f})')
    print(f'  {a:2d}-{a+4:2d} || '+' || '.join(out))
print('\nrain against the column water, ocean poleward of 55 deg (the same PW class in both hemispheres):  PW mm | S: area %, T C, RH %, rain, NASA | N: area %, T C, RH %, rain, NASA   (mm/a)')
for lo,hi in[(0,1.5),(1.5,2.5),(2.5,3.5),(3.5,5),(5,7),(7,10),(10,14)]:
    out=[]
    for sgn in(-1,1):
        m=oc&(sgn*LA>=55)&(PW>=lo)&(PW<hi); tot=oc&(sgn*LA>=55)
        out.append(f'{100*w[m].sum()/w[tot].sum():4.1f}, {g(T,m):6.1f}, {g(RH,m):4.1f}, {365*g(P,m):5.0f}, {365*g(N,m):5.0f}' if m.sum()>5 else '   --')
    print(f'  {lo:4.1f}-{hi:4.1f} | '+' | '.join(out))
m=oc&(LA<=-65); print(f'\nS ocean 65-90: surface T model {g(T,m):.1f} C, the NASA surface temperature field {g(TN,m):.1f} C; snow share of the rain {100*g(SN,m)/max(g(P,m),1e-9):.0f} %')
m=oc&(LA>=65);  print(f'N ocean 65-90: surface T model {g(T,m):.1f} C, the NASA surface temperature field {g(TN,m):.1f} C; snow share of the rain {100*g(SN,m)/max(g(P,m),1e-9):.0f} %')
