import sys,glob,re,numpy as np
# RC-RESID, READ-ONLY (2026-10-06): tropical-ocean rain at 28-29 C is 6.5 mm/d against NASA 4.8. Where, and what is it made of?   usage: sst2829.py <run> [iter]
run=sys.argv[1]; it=sys.argv[2] if len(sys.argv)>2 else '520'
def rd(f,want,lev0=None):
    L=open(f,errors='replace').read().split('\n'); o={}
    for n,l in enumerate(L):
        if l.startswith('SCALARS '):
            nm=l.split()[1]
            if nm in want and nm not in o: o[nm]=np.array(L[n+2:n+2+65341],float).reshape(181,361)
    return o
o=rd(f'output_{run}/0Ma_smooth_Atm_radial_0_{it}.vtk',('Topography','Precipitation','PrecipitationConv','Precipitation_NASA','Temperature','HumidityRel','Evaporation','PrecipitableWater','w-Component','v-Component'))
o0=rd(f'output_{run}/0Ma_smooth_Atm_radial_0_0.vtk',('Temperature','HumidityRel'))
lat=90-np.arange(181); lon=np.where(np.arange(361)<=180,np.arange(361),np.arange(361)-360); LA=lat[:,None]*np.ones((1,361)); LO=lon[None,:]*np.ones((181,1))
w=np.cos(np.radians(LA)); w[:,360]=0; oc=(o['Topography']==0)&(np.abs(LA)<=30)
P,C,N,T0,RH0,PW=o['Precipitation'],o['PrecipitationConv'],o['Precipitation_NASA'],o0['Temperature'],o0['HumidityRel'],o['PrecipitableWater']
g=lambda A,m:(A*w)[m].sum()/w[m].sum()
def wq(A,m,qs):
    a=A[m]; ww=w[m]; i=np.argsort(a); c=np.cumsum(ww[i])/ww.sum(); return [a[i][np.searchsorted(c,q)] for q in qs]
print(f'{run} {it}: |lat| <= 30 ocean by initial SST:  area % | NASA | model (conv + strat) | initial surface RH % | model p10/p50/p90 | NASA p10/p50/p90 | r(model, NASA) inside the class | model sd / NASA sd')
for lo,hi in[(26,27),(27,28),(28,28.5),(28.5,29),(29,29.5),(29.5,31)]:
    m=oc&(T0>=lo)&(T0<hi)
    if m.sum()<20: continue
    r=np.corrcoef(P[m],N[m])[0,1]
    print(f'  {lo:4.1f}-{hi:4.1f} | {100*w[m].sum()/w[oc].sum():4.1f} | {g(N,m):4.2f} | {g(P,m):4.2f} ({g(C,m):4.2f}+{g(P-C,m):4.2f}) | {g(RH0,m):4.1f} | '+'/'.join(f'{x:.1f}' for x in wq(P,m,(.1,.5,.9)))+' | '+'/'.join(f'{x:.1f}' for x in wq(N,m,(.1,.5,.9)))+f' | {r:+.2f} | {P[m].std():.2f} / {N[m].std():.2f}')
m28=oc&(T0>=28)&(T0<29)
print('\n28-29 C by region:   area % of the class | SST | NASA | model | model - NASA | NASA cells < 3 mm/d %')
B={'W Pacific warm pool 120E-180':(120,180),'C Pacific 180-140W':(-180,-140),'E Pacific 140W-78W':(-140,-78),'Atlantic 70W-15E':(-70,15),'W Indian 35-75E':(35,75),'E Indian 75-120E':(75,120)}
for nm,(a,b) in B.items():
    for lab,mm in(('N',LA>3),('eq',np.abs(LA)<=3),('S',LA<-3)):
        m=m28&(LO>=a)&(LO<b)&mm
        if w[m].sum()/w[m28].sum()>0.02: print(f'  {nm:28s} {lab:2s} | {100*w[m].sum()/w[m28].sum():4.1f} | {g(T0,m):5.2f} | {g(N,m):4.2f} | {g(P,m):4.2f} | {g(P-N,m):+5.2f} | {100*w[m&(N<3)].sum()/w[m].sum():3.0f}')
print('\n28-29 C by |latitude|:   area % | NASA | model | conv | initial RH | PW mm')
AL=np.abs(LA)
for lo,hi in[(0,3),(3,6),(6,10),(10,15),(15,20),(20,30)]:
    m=m28&(AL>=lo)&(AL<hi)
    if m.sum()>10: print(f'  {lo:2d}-{hi:2d} | {100*w[m].sum()/w[m28].sum():4.1f} | {g(N,m):4.2f} | {g(P,m):4.2f} | {g(C,m):4.2f} | {g(RH0,m):4.1f} | {g(PW,m):4.1f}')
print('\n28-29 C by NASA rain class (does the model know the dry cells?):   NASA class | area % | NASA | model | initial RH | |lat|')
for lo,hi in[(0,2),(2,4),(4,6),(6,8),(8,20)]:
    m=m28&(N>=lo)&(N<hi)
    if m.sum()>10: print(f'  {lo:2d}-{hi:2d} | {100*w[m].sum()/w[m28].sum():4.1f} | {g(N,m):4.2f} | {g(P,m):4.2f} | {g(RH0,m):4.1f} | {g(AL,m):4.1f}')
print(f'excess budget of the class: model - NASA = {g(P-N,m28):+.2f} mm/d; from cells where NASA < 4 mm/d: {(( P-N)*w)[m28&(N<4)].sum()/w[m28].sum():+.2f}, from NASA >= 4: {((P-N)*w)[m28&(N>=4)].sum()/w[m28].sum():+.2f}')
