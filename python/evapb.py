import sys,glob,re,numpy as np
# RC-RESID / P-E, READ-ONLY (2026-10-06): where is the evaporation missing?  Rebuilds the Meyer bulk formula of ThermoAtm::waterVapourEvaporation
# (closure branch) from the surface slice and splits E by band, by what limits it (wind term, humidity deficit), and by land / ocean.
# usage: evapb.py <run> [iteration]
run=sys.argv[1]
def rd(f,want):
    L=open(f,errors='replace').read().split('\n'); o={}; seen={}
    for n,l in enumerate(L):
        if l.startswith('SCALARS '):
            nm=l.split()[1]; seen[nm]=seen.get(nm,0)+1; key=nm if seen[nm]==1 else f'{nm}#{seen[nm]}'
            if nm in want: o[key]=np.array(L[n+2:n+2+65341],float).reshape(181,361)
    return o
fs=sorted(glob.glob(f'output_{run}/0Ma_smooth_Atm_radial_0_*.vtk'),key=lambda s:int(re.findall(r'_(\d+)\.vtk',s)[0]))
f=[x for x in fs if len(sys.argv)>2 and x.endswith(f'_{sys.argv[2]}.vtk')] or [fs[-1]]; f=f[0]
o=rd(f,('Topography','Precipitation','Precipitation_NASA','Temperature','HumidityRel','Evaporation','u-Component','v-Component','w-Component','v-velocity_NASA','w-velocity_NASA','WaterVapour','PressureStatic'))
print(f'{run} {f.split("_")[-1][:-4]}; fields:',' '.join(f'{k}[{v.min():.3g},{v.max():.3g}]' for k,v in o.items()))
land=o['Topography']>0; oc=~land; lat=90-np.arange(181); al=np.abs(lat)[:,None]*np.ones((1,361)); w=np.cos(np.radians(lat))[:,None]*np.ones((1,361)); w[:,360]=0
E=o['Evaporation']; P=o['Precipitation']; T=o['Temperature']; RH=o['HumidityRel']; q=o['WaterVapour']*1e-3; p=o['PressureStatic']
V=np.sqrt(o['u-Component']**2+o['v-Component']**2+o['w-Component']**2); Vc=V/np.sqrt(3.)          # the code's speed: sqrt((u2+v2+w2)/3)
Vn=np.hypot(o['v-velocity_NASA'],o['w-velocity_NASA'])
Es=np.where(T>=0,6.1078*np.exp(17.2694*T/(T+273.15-35.86)),6.1078*np.exp(21.8746*T/(T+273.15-7.66)))
cM=lambda v: 11.0*0.750062*(1+3.6*v/16.)/30.
sd=np.where(cM(Vc)>0,E/cM(Vc),0)                    # deficit the model used (level-1 humidity), hPa
g=lambda A,m:(A*w)[m].sum()/w[m].sum()
aw=lambda m: w[m].sum()/w.sum()
print(f'global: P {365*g(P,w>0):.0f}  E {365*g(E,w>0):.0f} mm/a;  ocean ({100*aw(oc):.0f} % of the area) E {365*g(E,oc):.0f}, P {365*g(P,oc):.0f};  land E {365*g(E,land):.0f}, P {365*g(P,land):.0f}  (Earth: ocean E ~1250-1400, land E ~480-550, global ~1000)')
print('ocean by |latitude|:  E mm/a | P mm/a | SST C | surface RH % | deficit used hPa (of E_sat) | wind in formula m/s | |V| level 0 m/s | Meyer wind factor 1+W/16')
for lo,hi in[(0,5),(5,15),(15,25),(25,35),(35,50),(50,65),(65,90)]:
    m=oc&(al>=lo)&(al<hi)
    print(f'  {lo:2d}-{hi:2d} | {365*g(E,m):5.0f} | {365*g(P,m):5.0f} | {g(T,m):5.1f} | {g(RH,m):4.1f} | {g(sd,m):5.2f} ({100*g(sd,m)/g(Es,m):3.0f} %) | {g(Vc,m):4.2f} | {g(V,m):4.2f} | {g(1+3.6*Vc/16,m):4.2f}')
print('\nWHAT-IF, same humidity deficit, ocean E in mm/a (global mean in brackets, land E = 0):')
for nm,vv in(('as run: sqrt((u2+v2+w2)/3) of level 0',Vc),('|V| of level 0 (no /3)',V),('|V| + 3 m/s gustiness in quadrature',np.hypot(V,3.0)),('|V| + 5 m/s gustiness in quadrature',np.hypot(V,5.0)),('7 m/s everywhere (Earth mean scalar wind)',np.full_like(V,7.0))):
    En=cM(vv)*sd
    print(f'  {nm:42s} {365*g(En,oc):5.0f} ({365*g(En,oc)*aw(oc):4.0f})   bands 0-15 / 15-35 / 35-65 / 65-90: '+' / '.join(f'{365*g(En,oc&(al>=a)&(al<b)):5.0f}' for a,b in((0,15),(15,35),(35,65),(65,90))))
print('  Earth (OAFlux-like, rough)                  ~1300 ( ~920)   bands ~1450 / ~1650 / ~950 / ~250')
print(f'\nland: P {365*g(P,land):.0f} mm/a on {100*aw(land):.0f} % of the area = {365*g(P,land)*aw(land):.0f} mm/a global; E_land = 0 (ATM_LAND_BUCKET off). If land returned 60 % of its rain (Earth ~0.6-0.65): +{0.6*365*g(P,land)*aw(land):.0f} mm/a global')
