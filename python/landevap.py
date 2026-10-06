import sys,glob,re,numpy as np
# RC-RESID / P-E, READ-ONLY (2026-10-06): land evaporation with ATM_LAND_BUCKET (wb53c: land E 728 mm/a > land P 653; Earth ~500, E/P ~0.6).
# Rebuilds beta and the potential (open-water Meyer) rate from the surface slice and tries a Budyko closure on the model's own rain.
# usage: landevap.py <run with the bucket on> [iter]
run=sys.argv[1]
def rd(f,want):
    L=open(f,errors='replace').read().split('\n'); o={}
    for n,l in enumerate(L):
        if l.startswith('SCALARS '):
            nm=l.split()[1]
            if nm in want and nm not in o: o[nm]=np.array(L[n+2:n+2+65341],float).reshape(181,361)
    return o
fs=sorted(glob.glob(f'output_{run}/0Ma_smooth_Atm_radial_0_*.vtk'),key=lambda s:int(re.findall(r'_(\d+)\.vtk',s)[0]))
f=[x for x in fs if len(sys.argv)>2 and x.endswith(f'_{sys.argv[2]}.vtk')] or [fs[-1]]; f=f[0]
o=rd(f,('Topography','Precipitation','Precipitation_NASA','Temperature','HumidityRel','Evaporation','Landscape'))
land=o['Topography']>0; lat=90-np.arange(181); al=np.abs(lat)[:,None]*np.ones((1,361)); w=np.cos(np.radians(lat))[:,None]*np.ones((1,361)); w[:,360]=0
E,P,N,T,RH=o['Evaporation'],np.maximum(o['Precipitation'],0),o['Precipitation_NASA'],o['Temperature'],o['HumidityRel']
CAP=150.; beta=np.minimum(1.,(CAP*np.minimum(1.,N*365/1000.))/(0.75*CAP)); Ep=np.where(beta>1e-6,E/np.maximum(beta,1e-6),np.nan)
g=lambda A,m:np.nansum((A*w)[m])/w[m&np.isfinite(A)].sum()
aw=lambda m:w[m].sum()/w.sum()
def budyko(Ep,P):
    phi=np.where(P>1e-6,Ep/np.maximum(P,1e-6),np.inf)
    return np.where(P>1e-6,P*np.sqrt(phi*np.tanh(1/np.maximum(phi,1e-9))*(1-np.exp(-phi))),0.)
Eb=budyko(np.nan_to_num(Ep,nan=0.),P)
oc=~land
print(f'{run} {f.split("_")[-1][:-4]}: global P {365*g(P,w>0):.0f}, E {365*g(E,w>0):.0f};  ocean E {365*g(E,oc):.0f};  land ({100*aw(land):.0f} % of the area) P {365*g(P,land):.0f}, E {365*g(E,land):.0f}  (Earth land: P ~750-800, E ~480-550, E/P ~0.62)')
print('land by |latitude|:  area % of land | P model | P NASA | E bucket | E/P | beta | potential (open water) | cells with E > P % | T C | RH % || Budyko on the model rain: E | E/P')
for lo,hi in[(0,15),(15,35),(35,65),(65,90)]:
    m=land&(al>=lo)&(al<hi)
    print(f'  {lo:2d}-{hi:2d} | {100*w[m].sum()/w[land].sum():4.1f} | {365*g(P,m):5.0f} | {365*g(N,m):5.0f} | {365*g(E,m):5.0f} | {g(E,m)/max(g(P,m),1e-9):4.2f} | {g(beta,m):4.2f} | {365*g(Ep,m):5.0f} | {100*w[m&(E>P)].sum()/w[m].sum():3.0f} | {g(T,m):5.1f} | {g(RH,m):4.1f} || {365*g(Eb,m):5.0f} | {g(Eb,m)/max(g(P,m),1e-9):4.2f}')
R={'Amazon':(-10,2,-70,-52),'Congo':(-6,5,14,28),'Sahara':(18,30,-10,25),'Sahel':(10,16,-10,30),'W Europe':(43,55,-5,15),'Siberia':(55,68,70,120),'E US':(32,42,-92,-78),'India':(10,28,72,88),'Australia':(-30,-18,118,145)}
lon=np.where(np.arange(361)<=180,np.arange(361),np.arange(361)-360)
print('regions:  P model | E bucket | potential | Budyko E | observed E (rough, mm/a)')
obs={'Amazon':1250,'Congo':1150,'Sahara':40,'Sahel':450,'W Europe':500,'Siberia':300,'E US':750,'India':700,'Australia':350}
for nm,(a,b,c,d) in R.items():
    m=land&((lat>=a)&(lat<=b))[:,None]&((lon>=c)&(lon<=d))[None,:]
    print(f'  {nm:10s} {365*g(P,m):5.0f} | {365*g(E,m):5.0f} | {365*g(Ep,m):5.0f} | {365*g(Eb,m):5.0f} | ~{obs[nm]}')
EoB=g(E,oc)*aw(oc)+g(Eb,land)*aw(land)
print(f'\nBudyko closure, land total: E {365*g(Eb,land):.0f} mm/a (E/P {g(Eb,land)/g(P,land):.2f});  global E {365*EoB:.0f}, P/E {g(P,w>0)/EoB:.2f}')
for s_ in(0.5,0.7,1.0):
    Ebs=budyko(s_*np.nan_to_num(Ep,nan=0.),P); Eg=g(E,oc)*aw(oc)+g(Ebs,land)*aw(land)
    print(f'   potential x {s_:.1f}: land E {365*g(Ebs,land):.0f} (E/P {g(Ebs,land)/g(P,land):.2f}), global E {365*Eg:.0f}, P/E {g(P,w>0)/Eg:.2f}')
