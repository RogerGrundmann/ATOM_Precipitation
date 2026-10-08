import sys,numpy as np
# RC-RESID, READ-ONLY (2026-10-08): polar land rain by ground height and by surface temperature (N land 65-90 243 / 418, Antarctica 282 / 204 on wb68:
# lowland starved, high ground ~2x too wet), and the evaporation rebuilt for other gustiness values (true deficit under ATM_EVAP_WIND=1 + ATM_EVAP_GUST=<G>).
# usage: [G=<gust the run used, default 3>] polarland.py <run> [iter]
import os
run=sys.argv[1]; it=sys.argv[2] if len(sys.argv)>2 else '520'
def rd(f,want):
    L=open(f,errors='replace').read().split('\n'); o={}
    for n,l in enumerate(L):
        if l.startswith('SCALARS '):
            nm=l.split()[1]
            if nm in want and nm not in o: o[nm]=np.array(L[n+2:n+2+65341],float).reshape(181,361)
    return o
s=rd(f'output_{run}/0Ma_smooth_Atm_radial_0_{it}.vtk',('Topography','Landscape','Precipitation','PrecipitationSnow','PrecipitationRain','Precipitation_NASA','Temperature','HumidityRel','PressureStatic','Evaporation','PrecipitableWater','v-Component','w-Component','u-Component','WaterVapour'))
lat=90-np.arange(181); LA=lat[:,None]*np.ones((1,361)); al=np.abs(LA); w=np.cos(np.radians(LA)); w[:,360]=0; land=s['Topography']>0; oc=~land
P,N,T,RH,ps,H,E=(s[x] for x in('Precipitation','Precipitation_NASA','Temperature','HumidityRel','PressureStatic','Landscape','Evaporation'))
g=lambda A,m:(A*w)[m].sum()/max(w[m].sum(),1e-30)
aw=lambda m:w[m].sum()/w.sum()
print('fields',{k:(round(float(v.min()),3),round(float(v.max()),3)) for k,v in s.items()})
print('\n== POLAR LAND (|lat|>=60) by ground height: area % of hemisphere polar land | model | NASA mm/a | T | p_sfc | RH')
for hem,mh in(('N',LA>=60),('S',LA<=-60)):
    ml=mh&land
    for a,b in[(0,200),(200,500),(500,1000),(1000,2000),(2000,5000)]:
        m=ml&(H>=a)&(H<b)
        if w[m].sum()<=0: continue
        print(f'  {hem} {a:4d}-{b:4d} m | {100*w[m].sum()/w[ml].sum():4.1f} | {365*g(P,m):5.0f} | {365*g(N,m):5.0f} | {g(T,m):6.1f} | {g(ps,m):5.0f} | {g(RH,m):4.1f}')
print('\n== POLAR LAND by surface temperature (|lat|>=60): hem | T class | area % | model | NASA | ground m')
for hem,mh in(('N',LA>=60),('S',LA<=-60)):
    ml=mh&land
    for a,b in[(-60,-35),(-35,-25),(-25,-18),(-18,-12),(-12,-6),(-6,5)]:
        m=ml&(T>=a)&(T<b)
        if w[m].sum()/w[ml].sum()<0.01: continue
        print(f'  {hem} {a:4d}..{b:4d} C | {100*w[m].sum()/w[ml].sum():4.1f} | {365*g(P,m):5.0f} | {365*g(N,m):5.0f} | {g(H,m):5.0f}')
print('\n== N land rows 50-80N by ground height < 300 m / >= 300 m: model | NASA')
for a,b in[(50,55),(55,60),(60,65),(65,70),(70,75),(75,80)]:
    r=''
    for lab,mm in(('<300',H<300),('>=300',H>=300)):
        m=(LA>=a)&(LA<b)&land&mm
        r+=f'   {lab}: area {100*w[m].sum()/max(w[(LA>=a)&(LA<b)&land].sum(),1e-9):3.0f} %  {365*g(P,m):5.0f} / {365*g(N,m):5.0f}'
    print(f'  {a}-{b}N'+r)
# regions
lon=np.where(np.arange(361)<=180,np.arange(361),np.arange(361)-360)
R={'W Siberia 60-72N 60-90E':(60,72,60,90),'E Siberia 60-72N 100-160E':(60,72,100,160),'N Europe 60-70N 5-40E':(60,70,5,40),'Canada 60-72N 130-70W':(60,72,-130,-70),'Alaska 60-70N 165-141W':(60,70,-165,-141),'Greenland 62-82N 60-20W':(62,82,-60,-20),'E Antarctica 70-85S 0-150E':(-85,-70,0,150),'W Antarctica 75-85S 150-60W':(-85,-75,-150,-60)}
print('\n== regions (land): model | NASA mm/a | ground m | T C')
for nm,(a,b,c,d) in R.items():
    m=land&((lat>=a)&(lat<=b))[:,None]&((lon>=c)&(lon<=d))[None,:]
    print(f'  {nm:30s} {365*g(P,m):5.0f} | {365*g(N,m):5.0f} | {g(H,m):5.0f} | {g(T,m):6.1f}')
# ---------- evaporation: true deficit under WIND=1 GUST=3
V=np.sqrt(s['u-Component']**2+s['v-Component']**2+s['w-Component']**2)
cM=lambda v: 11.0*0.750062*(1+3.6*v/16.)/30.
G0=float(os.environ.get('G','3')); W3=np.hypot(V,G0); sd=np.where(oc,E/cM(W3),0)
Es=np.where(T>=0,6.1078*np.exp(17.2694*T/(T+273.15-35.86)),6.1078*np.exp(21.8746*T/(T+273.15-7.66)))
print(f'\n== EVAPORATION as run (WIND=1, GUST={G0:g}): global P {365*g(P,w>0):.0f} E {365*g(E,w>0):.0f}; ocean E {365*g(E,oc):.0f} P {365*g(P,oc):.0f}; land E {365*g(E,land):.0f} P {365*g(P,land):.0f} (E/P {g(E,land)/g(P,land):.2f})')
print('ocean band | E | deficit hPa (% of E_sat) | |V| | W=hypot(V,G)')
for lo,hi in[(0,15),(15,25),(25,35),(35,50),(50,65),(65,90)]:
    m=oc&(al>=lo)&(al<hi)
    print(f'  {lo:2d}-{hi:2d} | {365*g(E,m):5.0f} | {g(sd,m):5.2f} ({100*g(sd,m)/g(Es,m):3.0f} %) | {g(V,m):4.2f} | {g(W3,m):4.2f}')
# land potential: invert Budyko numerically per cell to get pot from E,P
def bud(phi): return np.sqrt(phi*np.tanh(1/np.maximum(phi,1e-9))*(1-np.exp(-phi)))
Pm=np.maximum(P,0); ratio=np.where(Pm>1e-6,np.clip(E/np.maximum(Pm,1e-6),0,0.9999),0)
phi=np.full(P.shape,1.0); lo_=np.full(P.shape,1e-4); hi_=np.full(P.shape,200.)
for _ in range(60):
    mid=0.5*(lo_+hi_); f=bud(mid); up=f<ratio; lo_=np.where(up,mid,lo_); hi_=np.where(up,hi_,mid)
phi=0.5*(lo_+hi_); pot=phi*Pm
print(f'land potential rate (inverted): {365*g(np.where(land,pot,0),land):.0f} mm/a, mean phi-weighted; land P {365*g(Pm,land):.0f}')
print('\nWHAT-IF gustiness g (same deficit, same rain): ocean E | land E (Budyko, pot scaled by wind factor) | global E | P/E as is | P/E if P = 980')
Pg=365*g(P,w>0)
for gg in(3,4,5,6,7):
    Wg=np.hypot(V,gg); fac=cM(Wg)/cM(W3)
    Eo=E*fac; potg=pot*fac; phig=np.where(Pm>1e-6,potg/np.maximum(Pm,1e-6),0); El=np.where(Pm>1e-6,Pm*bud(np.maximum(phig,1e-9)),0)
    Et=np.where(oc,Eo,El); ge=365*g(Et,w>0)
    print(f'  g={gg}: {365*g(Eo,oc):5.0f} | {365*g(El,land):4.0f} (E/P {g(El,land)/g(Pm,land):.2f}) | {ge:5.0f} | {Pg/ge:.3f} | {980/ge:.3f}   ocean bands 0-15/15-35/35-65/65-90: '+' / '.join(f'{365*g(Eo,oc&(al>=a)&(al<b)):5.0f}' for a,b in((0,15),(15,35),(35,65),(65,90))))
print('  Earth rough: ocean ~1250-1400, land ~480-550 (E/P ~0.62), bands ~1450 / ~1650 / ~950 / ~250')
print('\nland E/P by band (as run): '+' | '.join(f'{a}-{b}: P {365*g(Pm,land&(al>=a)&(al<b)):.0f} NASA {365*g(N,land&(al>=a)&(al<b)):.0f} E {365*g(E,land&(al>=a)&(al<b)):.0f} pot {365*g(pot,land&(al>=a)&(al<b)):.0f}' for a,b in((0,15),(15,35),(35,65),(65,90))))
