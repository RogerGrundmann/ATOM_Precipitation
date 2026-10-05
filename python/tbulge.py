import sys,numpy as np
from scipy.ndimage import gaussian_filter1d
# SUBTROP option 2, READ-ONLY ESTIMATE (2026-10-05): how much warmer is the free air over elevated land than over the lowlands around it, and
# what would the convection trigger see if the column aloft followed the lowland reference instead?
# The model builds every column from its own sea-level temperature T_sl (modern: the NASA ground temperature projected down with the COSMO
# profile T(z) = T_sl*sqrt(1 - 2*beta*g*z/(R*T_sl^2)), beta = 44 K), so a warm plateau carries a warm column.
# usage: tbulge.py <run>        (reads output_<run>/..radial_0_0.vtk and ..longal_62_<last>.vtk)
run=sys.argv[1]; beta,R,g,cp,Lv,ep=44.,287.,9.81,1005.,2.5e6,0.622
def rd(f,want,n):
    L=open(f,errors='replace').read().split('\n'); o={}
    for i,l in enumerate(L):
        if l.startswith('SCALARS '):
            nm=l.split()[1]
            if nm in want and nm not in o: o[nm]=np.array(L[i+2:i+2+n[0]*n[1]],float).reshape(n)
    return o
o=rd(f'output_{run}/0Ma_smooth_Atm_radial_0_0.vtk',('Temperature','Topography','Landscape','Precipitation_NASA'),(181,361))
Tg=o['Temperature']+273.15; h=o['Landscape']; land=o['Topography']>0; h=np.where(land,h,0.)
lat=90-np.arange(181); lon=np.where(np.arange(361)<=180,np.arange(361),np.arange(361)-360); w=np.cos(np.radians(lat))[:,None]*np.ones((1,361))
c=2*beta*g/R
Tsl=np.sqrt(Tg**2+c*h)                                   # the column's sea-level temperature
prof=lambda T0,z: T0*np.sqrt(np.maximum(0.,1-c*z/T0**2))
# reference: sea-level temperature of the LOWLAND (ocean, or land below 300 m), smoothed with a Gaussian of sigma km (normalised convolution)
def ref(sigma_km):
    m=((~land)|(h<300.)).astype(float); A=Tsl*m; B=m.copy()
    s_lat=sigma_km/111.2
    for X in (A,B):
        for j in range(181):
            s=s_lat/max(np.cos(np.radians(lat[j])),0.05); X[j,:360]=gaussian_filter1d(X[j,:360],min(s,120.),mode='wrap'); X[j,360]=X[j,0]
        X[:]=gaussian_filter1d(X,s_lat,axis=0,mode='nearest')
    return A/np.maximum(B,1e-6)
for sig in (500.,1000.,2000.):
    Tr=ref(sig); d=Tsl-Tr
    m=land&(h>=300.)
    print(f'\n== reference = lowland sea-level temperature, Gaussian sigma {sig:.0f} km.  dT = column minus reference (K), land above 300 m')
    print('   |lat| band:  area % of land | mean dT | p10 | p90 | share with dT > +1 K')
    for a,b in[(0,15),(15,35),(35,52),(52,90)]:
        mm=m&((abs(lat)>=a)&(abs(lat)<b))[:,None]; x=d[mm]
        print(f'   {a:2d}-{b:2d}: {100*w[mm].sum()/w[land&((abs(lat)>=a)&(abs(lat)<b))[:,None]].sum():5.1f} | {np.average(x,weights=w[mm]):+5.1f} | {np.percentile(x,10):+5.1f} | {np.percentile(x,90):+5.1f} | {100*np.mean(x>1):3.0f} %')
    if sig==1000.:
        print('   by ground elevation, |lat| < 35:  mean dT | p10 | p90')
        for lo,hi in[(300,800),(800,1300),(1300,2000),(2000,3000),(3000,9000)]:
            mm=land&(h>lo)&(h<=hi)&(abs(lat)<35)[:,None]; x=d[mm]; print(f'   {lo:5d}-{hi:5d} m: {np.average(x,weights=w[mm]):+5.1f} | {np.percentile(x,10):+5.1f} | {np.percentile(x,90):+5.1f}')
        reg={'Mexican plateau':(18,30,-108,-96),'Highveld':(-32,-22,24,32),'E Africa':(-10,8,30,42),'Ethiopia':(5,14,35,42),'Congo':(-5,5,15,28),'S Brazil plateau':(-30,-18,-55,-42),'Altiplano / Andes':(-25,-10,-72,-64),
             'Iran':(28,37,48,62),'Tibet':(28,38,78,100),'Deccan':(12,24,74,82),'Sahara':(18,30,-10,25),'Arabia':(16,28,42,56),'W US / Rockies':(35,48,-120,-104),'Mongolia':(42,50,90,115),'Greenland':(65,80,-50,-30),'E Antarctica':(-85,-70,0,140)}
        print('   regions (land above 300 m):  mean ground m | ground T C | column T_sl C | reference T_sl C | dT K | NASA rain mm/d')
        for nm,(a,b,cc,dd) in reg.items():
            mm=land&(h>=300)&((lat>=a)&(lat<=b))[:,None]&((lon>=cc)&(lon<=dd))[None,:]; g_=lambda A:(A*w)[mm].sum()/w[mm].sum()
            print(f'   {nm:18s} {g_(h):5.0f} | {g_(Tg)-273.15:5.1f} | {g_(Tsl)-273.15:5.1f} | {g_(Tr)-273.15:5.1f} | {g_(d):+5.1f} | {g_(o["Precipitation_NASA"]):4.2f}')
        Tr1000=Tr
# ---- what the trigger would see at 28N (true model profiles from the longal_62 section)
import glob,re
f=sorted(glob.glob(f'output_{run}/0Ma_smooth_Atm_longal_62_*.vtk'),key=lambda s:int(re.findall(r'_(\d+)\.vtk',s)[0]))[-1]
q_=rd(f,('Temperature','WaterVapour','PressureStatic','height','Topography','PrecipitationConv'),(41,361))
T=q_['Temperature']+273.15; q=q_['WaterVapour']*1e-3; p=q_['PressureStatic']; z=q_['height']*1e3; rock=q_['Topography']>0
qsat=lambda T,p: ep*(6.1078*np.exp(17.2694*(T-273.15)/(T-35.86)))/(p-6.1078*np.exp(17.2694*(T-273.15)/(T-35.86)))
thE=lambda T,p,q: T*(1000./p)**(R/cp)*np.exp(Lv*q/(cp*T))
j28=90-28
def trig(k,D0,D1):
    """margin (K), CAPE, largest lapse K/km in the transition; environment replaced by the reference column above zg+D1, blended from zg+D0"""
    i0=np.where(~rock[:,k])[0][0] if rock[0,k] else 0; zg=z[i0,k]
    ml=[i for i in range(i0+1,40) if z[i,k]-zg<=500.] or [i0+1]
    ww=np.gradient(z[:,k])[ml]*p[ml,k]/T[ml,k]
    Tm=(ww*T[ml,k]).sum()/ww.sum(); qm=(ww*q[ml,k]).sum()/ww.sum(); pm=(ww*p[ml,k]).sum()/ww.sum(); zm=(ww*z[ml,k]).sum()/ww.sum()
    Tp=Tm+0.15+(2.3 if i0>0 else 0); qp=qm+4.5e-4; the=thE(Tp,pm,qp)
    wt=np.clip((z[:,k]-zg-D0)/max(D1-D0,1.),0,1) if D1>0 else np.zeros(41)
    Te=T[:,k]-wt*(prof(Tsl[j28,k],z[:,k])-prof(Tr1000[j28,k],z[:,k]))*(1 if i0>0 else 0)
    base=None
    for i in range(i0+1,39):
        if p[i,k]<500: break
        if qp>=qsat(Tp-g*(z[i,k]-zm)/cp,p[i,k]): base=i; break
    lev=[i for i in range((base if base else i0)+1,40) if z[i,k]<12000]
    thes=np.array([thE(Te[i],p[i,k],qsat(Te[i],p[i,k])) for i in lev])
    cape=0.; b=False
    for n,i in enumerate(lev):
        if the-thes[n]>0: b=True; cape+=g*(the-thes[n])/thes[n]*(z[i,k]-z[i-1,k])
        elif b: break
    up=[i for i in range(i0,39) if z[i,k]-zg<=max(D1,500)+600]
    lapse=max(-(Te[i+1]-Te[i])/(z[i+1,k]-z[i,k])*1e3 for i in up)
    return the-thes.min(),cape,lapse,i0>0,zg
print(f'\n== trigger at 28N ({f.split("/")[-1]}), environment aloft replaced by the sigma-1000-km lowland reference; blend from D0 to D1 above ground')
reg={'Mexican plateau (108-101W)':(-108,-101),'Gulf coast (99-97W)':(-99,-97),'Sahara (0-25E)':(0,25),'Arabia / Iran (45-60E)':(45,60),'N India (75-83E)':(75,83),'Himalaya / Tibet edge (84-98E)':(84,98),'S China (105-119E)':(105,119)}
print('   region: dT_sl K | margin now | margin with (D0,D1) = (500,1500) m | (0,1000) m | buoyant % now -> (500,1500) | CAPE J/kg | steepest lapse K/km now -> (500,1500)')
for nm,(a,b) in reg.items():
    ks=[k for k in range(360) if a<=lon[k]<=b]
    N=[trig(k,0,0) for k in ks]; A=[trig(k,500,1500) for k in ks]; B=[trig(k,0,1000) for k in ks]
    mean=lambda X,n: np.mean([x[n] for x in X])
    print(f'   {nm:30s} {np.mean([Tsl[j28,k]-Tr1000[j28,k] for k in ks if land[j28,k]] or [0]):+5.1f} | {mean(N,0):+6.1f} | {mean(A,0):+6.1f} | {mean(B,0):+6.1f} | {100*np.mean([x[0]>0 for x in N]):3.0f} -> {100*np.mean([x[0]>0 for x in A]):3.0f} | {mean(A,1):6.1f} | {mean(N,2):4.1f} -> {mean(A,2):4.1f}')
Ld=[k for k in range(360) if land[j28,k]]
print('   all land columns at 28N by ground elevation: n | margin now -> (500,1500) | buoyant % now -> then')
for lo,hi in[(0,300),(300,800),(800,1300),(1300,2000),(2000,3000),(3000,9000)]:
    ks=[k for k in Ld if lo<h[j28,k]<=hi]
    if ks:
        N=[trig(k,0,0) for k in ks]; A=[trig(k,500,1500) for k in ks]
        print(f'   {lo:5d}-{hi:5d} m: {len(ks):3d} | {np.mean([x[0] for x in N]):+6.1f} -> {np.mean([x[0] for x in A]):+6.1f} | {100*np.mean([x[0]>0 for x in N]):3.0f} -> {100*np.mean([x[0]>0 for x in A]):3.0f}')
