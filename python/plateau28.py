import sys,numpy as np
# Why do plateau columns not convect?  Rebuilds the scheme's own trigger (MoistConvection::mlParcel + findCloudBaseLFS, ML-LCL branch) column by
# column on the 28N section (longal_62) of a run.  usage: plateau28.py <run> [iteration]     READ-ONLY, 2026-10-05 (SUBTROP)
run=sys.argv[1]; it=sys.argv[2] if len(sys.argv)>2 else '520'
f=f'output_{run}/0Ma_smooth_Atm_longal_62_{it}.vtk'
W=('Temperature','WaterVapour','PressureStatic','height','Topography','HumidityRel','PrecipitationConv','Precipitation','M_u','Cloud_Base','CloudWater')
L=open(f,errors='replace').read().split('\n'); o={}
for n,l in enumerate(L):
    if l.startswith('SCALARS '):
        nm=l.split()[1]
        if nm in W and nm not in o: o[nm]=np.array(L[n+2:n+2+14801],float).reshape(41,361)
T=o['Temperature']+273.15; q=o['WaterVapour']*1e-3; p=o['PressureStatic']; z=o['height']*1e3; rock=o['Topography']>0
lon=np.where(np.arange(361)<=180,np.arange(361),np.arange(361)-360)
cp,Rd,g,Lv,ep=1005.,287.,9.81,2.5e6,0.622
qsat=lambda T,p: ep*(6.1078*np.exp(17.2694*(T-273.15)/(T-35.86)))/(p-6.1078*np.exp(17.2694*(T-273.15)/(T-35.86)))
thE=lambda T,p,q: T*(1000./p)**(Rd/cp)*np.exp(Lv*q/(cp*T))
T_ADD,Q_ADD,T_LAND=0.15,4.5e-4,2.3
def col(k):
    air=np.where(~rock[:,k])[0]
    i0=air[0] if rock[0,k] else 0            # first air level = i_topography (0 over ocean)
    land=i0>0; zg=z[i0,k]
    ml=[i for i in range(i0+1,40) if z[i,k]-zg<=500.] or [i0+1]
    w=np.gradient(z[:,k])[ml]*p[ml,k]/T[ml,k]
    Tm=(w*T[ml,k]).sum()/w.sum(); qm=(w*q[ml,k]).sum()/w.sum(); pm=(w*p[ml,k]).sum()/w.sum(); zm=(w*z[ml,k]).sum()/w.sum()
    Tp=Tm+T_ADD+(T_LAND if land else 0); qp=qm+Q_ADD; the=thE(Tp,pm,qp)
    base=None
    for i in range(i0+1,39):
        if p[i,k]<500: break
        if qp>=qsat(Tp-g*(z[i,k]-zm)/cp,p[i,k]): base=i; break
    lev=range((base if base else i0)+1,40)
    thes=np.array([thE(T[i,k],p[i,k],qsat(T[i,k],p[i,k])) for i in lev]); zz=np.array([z[i,k] for i in lev])
    sel=zz<12000; imin=thes[sel].argmin()
    cape=0.; b=False
    for n,i in enumerate(lev):
        if the-thes[n]>0: b=True; cape+=g*(the-thes[n])/thes[n]*(z[i,k]-z[i-1,k])
        elif b: break
    # parcel humidity that would just reach buoyancy at the scheme's own parcel temperature
    qn=cp*Tp/Lv*np.log(thes[sel][imin]/(Tp*(1000./pm)**(Rd/cp)))
    return dict(lon=lon[k],land=land,zg=zg,Tg=T[i0,k]-273.15,Tp=Tp-273.15,qp=qp*1e3,rh=qm/qsat(Tm,pm),pm=pm,the=the,thes=thes[sel][imin],zmin=zz[sel][imin],
                marg=the-thes[sel][imin],cape=cape,base=(z[base,k]-zg if base else np.nan),qn=qn*1e3,rhn=(qn-Q_ADD)/qsat(Tm,pm),
                pc=o['PrecipitationConv'][i0,k] if 'PrecipitationConv' in o else np.nan,T5=np.interp(5000,z[:,k],T[:,k])-273.15,theta=Tp*(1000./pm)**(Rd/cp))
C=[col(k) for k in range(360)]
reg={'E Pacific ocean (130-118W)':(-130,-118),'Mexican plateau (108-101W)':(-108,-101),'Gulf coast / S Texas (99-97W)':(-99,-97),'Gulf of Mexico (95-85W)':(-95,-85),'Florida (82-81W)':(-82,-81),
     'Sahara (0-25E)':(0,25),'Arabia / Iran (45-60E)':(45,60),'N India plain (75-83E)':(75,83),'Himalaya / Tibet edge (84-98E)':(84,98),'S China (105-119E)':(105,119),'W Pacific ocean (125-150E)':(125,150)}
print(f'{run} iteration {it}, 28N.  ground m | parcel T C, theta K, q g/kg, ML RH | parcel theta_e K | env theta_es min K at m | margin K | CAPE J/kg | LCL above ground m | q (ML RH) needed for margin 0 | T(5 km) C | P_conv mm/d')
for nm,(a,b) in reg.items():
    S=[c for c in C if a<=c['lon']<=b]
    if not S: continue
    m=lambda key: np.nanmean([c[key] for c in S])
    print(f'{nm:32s} {m("zg"):5.0f} | {m("Tp"):5.1f} {m("theta"):6.1f} {m("qp"):5.2f} {m("rh"):.2f} | {m("the"):6.1f} | {m("thes"):6.1f} at {m("zmin"):5.0f} | {m("marg"):+6.1f} | {m("cape"):6.1f} | {m("base"):5.0f} | {m("qn"):5.2f} ({m("rhn"):.2f}) | {m("T5"):5.1f} | {m("pc"):5.2f}')
print('\nland columns at 28N by ground elevation:  n | parcel theta K | q g/kg | ML RH | theta_e | env theta_es min | margin K | buoyant % | ML RH needed')
Ld=[c for c in C if c['land']]
for lo,hi in[(0,300),(300,800),(800,1300),(1300,2000),(2000,3000),(3000,9000)]:
    S=[c for c in Ld if lo<c['zg']<=hi]
    if S: m=lambda key: np.mean([c[key] for c in S]); print(f'  {lo:5d}-{hi:5d} m: {len(S):3d} | {m("theta"):6.1f} | {m("qp"):5.2f} | {m("rh"):.2f} | {m("the"):6.1f} | {m("thes"):6.1f} | {m("marg"):+6.1f} | {100*np.mean([c["marg"]>0 for c in S]):3.0f} | {m("rhn"):.2f}')
