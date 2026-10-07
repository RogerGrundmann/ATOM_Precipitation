import sys,numpy as np
# READ-ONLY (2026-10-07): the surface wind keeps its initial state (wb68: pattern r 0.989 / 0.994 after 520 iterations) and the initial meridional
# wind is zonally uniform, so the model has no convergence lines. Would a surface wind DERIVED FROM THE SEA TEMPERATURE (Lindzen-Nigam boundary layer:
# p_s' = -G*T_s', eps*u - f*v = -p_x/rho, eps*v + f*u = -p_y/rho) carry the part of the observed rain that the sea temperature itself does not?
# usage: lnwind.py <run> [iter]
run=sys.argv[1]; it=sys.argv[2] if len(sys.argv)>2 else '520'
L=open(f'output_{run}/0Ma_smooth_Atm_radial_0_{it}.vtk',errors='replace').read().split('\n'); o={}
for i,l in enumerate(L):
    if l.startswith('SCALARS '):
        nm=l.split()[1]
        if nm not in o and nm in('Topography','Temperature','Precipitation','PrecipitationConv','Precipitation_NASA','v-Component','w-Component','WaterVapour'): o[nm]=np.array(L[i+2:i+2+65341],float).reshape(181,361)[:,:360]
lat=90.0-np.arange(181); phi=np.radians(lat)[:,None]; LA=lat[:,None]*np.ones((1,360)); LO=np.arange(360)[None,:]*np.ones((181,1))
w8=np.cos(phi)*np.ones((1,360)); oc=o['Topography']==0; T=o['Temperature']; N=o['Precipitation_NASA']; P=o['Precipitation']; Pc=o['PrecipitationConv']
a=6.371e6; d=np.radians(1.0); cs=np.maximum(np.cos(phi),0.05); om=7.292e-5; f=2*om*np.sin(phi)
def fill(A,m,n=400):                      # Laplace fill of the land cells from the ocean values (periodic in longitude)
    B=np.where(m,A,np.nan); z=np.nanmean(B,axis=1,keepdims=True); z=np.where(np.isnan(z),np.nanmean(B),z); B=np.where(m,A,z*np.ones_like(A))
    for _ in range(n):
        S=(np.roll(B,1,1)+np.roll(B,-1,1)+np.vstack([B[:1],B[:-1]])+np.vstack([B[1:],B[-1:]]))/4; B=np.where(m,A,S)
    return B
def sm(A,n):
    if n==0: return A
    B=np.zeros_like(A); c=0
    for dj in range(-n,n+1):
        Aj=A[np.clip(np.arange(181)+dj,0,180)]
        for dk in range(-n,n+1): B+=np.roll(Aj,dk,1); c+=1
    return B/c
ddx=lambda A:(np.roll(A,-1,1)-np.roll(A,1,1))/(2*d*a*cs)            # eastward
ddy=lambda A:-np.gradient(A,axis=0)/(d*a)                            # northward (j runs southward)
div=lambda u,v:ddx(u)+ddy(v*cs)/cs
def ln(Ts,days,G=50.0,rho=1.2):
    eps=1.0/(days*86400.0); p=-G*Ts; px,py=ddx(p),ddy(p); D=rho*(eps*eps+f*f)
    u=-(eps*px+f*py)/D; v=-(eps*py-f*px)/D; return u,v
def wr(A,B,m):
    ww=w8[m]; x=A[m]-(A[m]*ww).sum()/ww.sum(); y=B[m]-(B[m]*ww).sum()/ww.sum(); return (ww*x*y).sum()/np.sqrt((ww*x*x).sum()*(ww*y*y).sum())
Tf=fill(T,oc)
# keep away from the coast (the fill is not a sea temperature): ocean cells with ocean within +-3 cells
far=oc.copy()
for dj in range(-3,4):
    for dk in range(-3,4): far&=np.roll(oc[np.clip(np.arange(181)+dj,0,180)],dk,1)
m=far&(abs(LA)<=28)
res=np.zeros_like(N); cls=np.zeros_like(N); bins=np.floor(T*2)/2
for b in np.unique(bins[m]):
    k=m&(bins==b); cls[k]=N[k].mean(); res[k]=N[k]-cls[k]
resM=np.zeros_like(N)
for b in np.unique(bins[m]):
    k=m&(bins==b); resM[k]=P[k]-P[k].mean()
print(f'{run} {it}: tropical ocean |lat| <= 28, >= 3 cells from land: {m.sum()} cells. r(NASA, SST) {wr(N,T,m):.3f}; r(NASA, its SST-class mean) {wr(N,cls,m):.3f}; r(model rain, NASA) {wr(P,N,m):.3f}')
print(f'  NASA rain not explained by the SST class: rms {np.sqrt((w8[m]*res[m]**2).sum()/w8[m].sum()):.2f} of {np.sqrt((w8[m]*(N[m]-(N[m]*w8[m]).sum()/w8[m].sum())**2).sum()/w8[m].sum()):.2f} mm/d; the model rain has {np.sqrt((w8[m]*resM[m]**2).sum()/w8[m].sum()):.2f} mm/d inside the classes; r(model residual, NASA residual) {wr(resM,res,m):.3f}')
vm,wm=o['v-Component'],o['w-Component']
Cm=sm(-div(wm,-vm),2)      # model wind: v is southward-positive (checked: gives convergence at the equator)
print(f'  model wind at level 0:   convergence vs NASA rain r {wr(Cm,N,m):+.3f}, vs the NASA residual r {wr(Cm,res,m):+.3f}')
print('  wind derived from the sea temperature (smoothing of T +-n deg, damping time):  r(conv, NASA rain) | r(conv, NASA residual) | r(conv, SST) | rms conv 1e-6/s | rms u, v m/s')
best=None
for n in(1,2,3):
    Ts=sm(Tf,n)
    for days in(1.0,2.5,5.0):
        u,v=ln(Ts,days); C=sm(-div(u,v),1)
        rN,rR,rT=wr(C,N,m),wr(C,res,m),wr(C,T,m)
        print(f'    n {n}, {days:3.1f} d: {rN:+.3f} | {rR:+.3f} | {rT:+.3f} | {1e6*np.sqrt((C[m]**2).mean()):5.2f} | {np.sqrt((u[m]**2).mean()):.2f}, {np.sqrt((v[m]**2).mean()):.2f}')
        if best is None or rR>best[0]: best=(rR,n,days,C,u,v)
    Lp=sm(-(ddx(ddx(Ts))+ddy(ddy(Ts)*cs)/cs),1)
    print(f'    n {n}, -Laplacian(SST): {wr(Lp,N,m):+.3f} | {wr(Lp,res,m):+.3f} | {wr(Lp,T,m):+.3f}')
rR,n,days,C,u,v=best
print(f'  best on the residual: n {n}, {days} d (r {rR:+.3f}).  NASA = class mean + b*conv: r {wr(cls+np.polyfit(C[m],res[m],1)[0]*C,N,m):.3f} (class mean alone {wr(N,cls,m):.3f})')
b=np.polyfit(C[m],res[m],1)[0]; P2=np.maximum(P+b*C,0)
print(f'  model rain + b*conv (same b, floored at 0): r with NASA {wr(P2,N,m):.3f} (model alone {wr(P,N,m):.3f}); mean {(P2*w8)[m].sum()/w8[m].sum():.2f} vs {(P*w8)[m].sum()/w8[m].sum():.2f} mm/d')
print('  regions (ocean, >= 3 cells from land):  NASA | model | NASA residual after SST | derived convergence 1e-6/s | model + b*conv')
for nm,la0,la1,lo0,lo1 in(('E Pacific 4-16S 130-90W',-16,-4,230,270),('E Pacific 4-12N 130-90W',4,12,230,270),('C Pacific 5S-5N 180-150W',-5,5,180,210),('W Pacific 0-10N 140-170E',0,10,140,170),('SPCZ 5-15S 160E-170W',-15,-5,160,190),('S Atlantic 4-16S 30W-0',-16,-4,330,360),('Arabian Sea 8-20N 55-70E',8,20,55,70),('Indian 0-10S 60-90E',-10,0,60,90)):
    k=far&(LA>=la0)&(LA<=la1)&(LO>=lo0)&(LO<lo1)
    if k.sum()>20:
        g=lambda A:(A*w8)[k].sum()/w8[k].sum(); print(f'    {nm:26s} {g(N):5.2f} | {g(P):5.2f} | {g(res):+5.2f} | {1e6*g(C):+6.2f} | {g(P2):5.2f}')
# --- extension: stronger damping (the limit eps -> infinity is -Laplacian(SST)) and wider smoothing
print('\n  wider smoothing / stronger damping:  r(pred, NASA rain) | r(pred, NASA residual) | NASA = class mean + b*pred: r')
keep={}
for n in(3,4,5,6,8):
    Ts=sm(Tf,n); Lp=-(ddx(ddx(Ts))+ddy(ddy(Ts)*cs)/cs); keep[n]=Lp
    row=f'    n {n}: -Laplacian {wr(Lp,N,m):+.3f} | {wr(Lp,res,m):+.3f} | {wr(cls+np.polyfit(Lp[m],res[m],1)[0]*Lp,N,m):.3f}'
    for days in(0.25,0.5):
        u,v=ln(Ts,days); C=-div(u,v); row+=f'   LN {days} d {wr(C,N,m):+.3f} | {wr(C,res,m):+.3f}'
    print(row)
Lp=keep[5]; b=np.polyfit(Lp[m],res[m],1)[0]; P2=np.maximum(P+b*Lp,0)
print(f'  -Laplacian, n 5: model rain + b*pred (floored at 0): r with NASA {wr(P2,N,m):.3f} (model alone {wr(P,N,m):.3f}); inside the SST classes the model would then have {np.sqrt((w8[m]*(b*Lp[m])**2).sum()/w8[m].sum()):.2f} mm/d rms (NASA {np.sqrt((w8[m]*res[m]**2).sum()/w8[m].sum()):.2f})')
print('  regions:  NASA | model | NASA residual after SST | b*pred mm/d | model + b*pred')
for nm,la0,la1,lo0,lo1 in(('E Pacific 4-16S 130-90W',-16,-4,230,270),('E Pacific 4-12N 130-90W',4,12,230,270),('C Pacific 5S-5N 180-150W',-5,5,180,210),('W Pacific 0-10N 140-170E',0,10,140,170),('SPCZ 5-15S 160E-170W',-15,-5,160,190),('S Atlantic 4-16S 30W-0',-16,-4,330,360),('Arabian Sea 8-20N 55-70E',8,20,55,70),('Indian 0-10S 60-90E',-10,0,60,90)):
    k=far&(LA>=la0)&(LA<=la1)&(LO>=lo0)&(LO<lo1)
    if k.sum()>20:
        g=lambda A:(A*w8)[k].sum()/w8[k].sum(); print(f'    {nm:26s} {g(N):5.2f} | {g(P):5.2f} | {g(res):+5.2f} | {g(b*Lp):+5.2f} | {g(P2):5.2f}')
# is the residual zonal or meridional? the zonal-mean part of the NASA residual and what is left
zr=np.where(m,res,np.nan); zmean=np.nanmean(zr,axis=1,keepdims=True)*np.ones_like(res); zmean=np.where(np.isnan(zmean),0,zmean)
tot=(w8[m]*res[m]**2).sum(); print(f'  of the NASA residual variance, the zonal-mean (latitude-only) part is {100*(w8[m]*zmean[m]**2).sum()/tot:.0f} %; r(-Laplacian n 5, zonal-mean part) {wr(Lp,zmean,m):+.3f}, r(-Laplacian, the zonally varying part) {wr(Lp,res-zmean,m):+.3f}')
