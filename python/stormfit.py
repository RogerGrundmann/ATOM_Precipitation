import sys,itertools,numpy as np
# STORM-SHAPE, READ-ONLY (2026-10-08): OFFLINE estimate of the ocean rain poleward of 48 deg for other settings of the high-latitude humidity knobs
# (ATM_RH_STORM_SST slope / _MAX / _REF / _LAT, ATM_RH_STORM_POLAR), from a run's surface slice: each ocean cell's initial RH factor is rebuilt for the
# run's settings and for the trial, and the rain is scaled by exp(S * (f_new / f_old - 1)), S = 10 (wb62b/c: 0.01 in RH ~ 11 % of the rain near the
# storm peak; wb56: Arctic +52 % for +5..7 %). An ESTIMATE to choose screen arms, not a result.
# usage: [S=10] [PARS='slope,max,ref,lat,polar;...'] stormfit.py <run> [iter]   (without PARS: a scan, the 14 best by the row misfit)
run=sys.argv[1]; it=sys.argv[2] if len(sys.argv)>2 else '520'; S=float(__import__('os').environ.get('S','10'))
def rd(f,want):
    L=open(f,errors='replace').read().split('\n'); o={}
    for n,l in enumerate(L):
        if l.startswith('SCALARS '):
            nm=l.split()[1]
            if nm in want and nm not in o: o[nm]=np.array(L[n+2:n+2+65341],float).reshape(181,361)
    return o
s=rd(f'output_{run}/0Ma_smooth_Atm_radial_0_{it}.vtk',('Topography','Precipitation','Precipitation_NASA','Temperature'))
lat=90-np.arange(181); LA=lat[:,None]*np.ones((1,361)); al=np.abs(LA); w=np.cos(np.radians(LA)); w[:,360]=0; oc=s['Topography']<=0
P,N,T=s['Precipitation'],s['Precipitation_NASA'],s['Temperature']
g=lambda A,m:(A*w)[m].sum()/max(w[m].sum(),1e-30)
ss=lambda x:np.where(x<=0,0.,np.where(x>=1,1.,x*x*(3-2*x)))
def fac(k,mx,ref,l0,pol,st=1.13,wid=17.,c=55.):
    G=1+(st-1)*np.exp(-((al-c)/wid)**2); f=np.where(al>c,np.maximum(G,pol),G)
    return f*(1-ss((al-l0)/10.)*np.clip(k*(T-ref),-mx,mx))
BASE=(0.005,0.04,4.,48.,1.08); f0=fac(*BASE)
def trial(par):
    Pn=np.where(oc,P*np.exp(S*(fac(*par)/f0-1)),P); return Pn
def r_(Pn):
    m=w>0; a=Pn-g(Pn,m); b=N-g(N,m); return ((a*b*w)[m].sum())/np.sqrt(((a*a*w)[m].sum())*((b*b*w)[m].sum()))
rows=[(50,54),(54,58),(58,62),(62,66),(66,70),(70,74),(74,82)]
def line(par,Pn):
    o=f'{par[0]:.4f} {par[1]:.2f} {par[2]:4.1f} {par[3]:2.0f} {par[4]:.2f} | {365*g(Pn,w>0):6.1f} {r_(Pn):.4f} |'
    for hem in(-1,1):
        o+=' '+' '.join(f'{g(Pn,oc&(LA*hem>=a)&(LA*hem<b))/max(g(N,oc&(LA*hem>=a)&(LA*hem<b)),1e-9):4.2f}' for a,b in rows)+' |'
    o+=f' {365*g(Pn,oc&(LA<=-65)):4.0f} {365*g(Pn,oc&(LA>=65)):4.0f} |'
    m5=oc&(al>=50)&(al<66)
    o+=' '+' '.join(f'{g(Pn,m5&(T>=a)&(T<b))/g(N,m5&(T>=a)&(T<b)):4.2f}' for a,b in((-6,-2),(-2,1),(1,3),(3,5),(5,7),(7,9),(9,14)))
    # misfit: rms log ratio over rows of both hemispheres, area weighted
    e=0;ws=0
    for hem in(-1,1):
        for a in range(48,86,2):
            m=oc&(LA*hem>=a)&(LA*hem<a+2)
            if w[m].sum()>0: e+=w[m].sum()*np.log(g(Pn,m)/g(N,m))**2; ws+=w[m].sum()
    return o+f' | {np.sqrt(e/ws):.3f}',np.sqrt(e/ws)
print(f'{run} {it}, S={S:g}.  slope MAX REF LAT POLAR | global mm/a  r | S rows 50-54 54-58 58-62 62-66 66-70 70-74 74-82 (model/NASA) | N rows the same | ocean 65-90 S N (NASA 490 379) | 50-66 deg by T class -6..-2 -2..1 1..3 3..5 5..7 7..9 9..14 | rms log misfit of 2-deg rows')
print(line(BASE,P)[0],' <- as run')
import os
if os.environ.get('PARS'):
    for q in os.environ['PARS'].split(';'): par=tuple(float(x) for x in q.split(',')); print(line(par,trial(par))[0])
    sys.exit()
res=[]
for k,mx,ref,l0,pol in itertools.product((0.005,0.007,0.009,0.011),(0.04,0.06,0.08,0.10),(2.,3.,4.),(48.,),(1.0,1.04,1.06,1.08)):
    par=(k,mx,ref,l0,pol); Pn=trial(par); l,e=line(par,Pn); res.append((e,l))
res.sort()
for e,l in res[:14]: print(l)
