import sys,numpy as np
# rain-rate distribution over the tropical ocean (|lat| <= 30), model vs NASA, from a radial VTK slice (mm/d)
for f in sys.argv[1:]:
    L=open(f,errors='replace').read().split('\n'); r={}
    for i,l in enumerate(L):
        if l.startswith('SCALARS '):
            n=l.split()[1]
            if n in('Precipitation','Precipitation_NASA','Topography') and n not in r: r[n]=np.array(L[i+2:i+2+65341],float).reshape(181,361)
    lat=90-np.arange(181); oc=(r['Topography']<=0)&(abs(lat)<=30)[:,None]; W=(np.cos(np.radians(lat))[:,None]*np.ones((1,361)))[oc]
    for name,A in(('model',r['Precipitation'][oc]),('NASA',r['Precipitation_NASA'][oc])):
        o=np.argsort(A); c=np.cumsum(W[o])/W.sum(); q=lambda x:A[o][min(np.searchsorted(c,x),len(A)-1)]
        print(f'  ocean<=30 {name:5s}: mean {(A*W).sum()/W.sum():.2f} mm/d  p10/p25/p50/p75/p90/p99 '+' / '.join('%.1f'%q(x) for x in(.1,.25,.5,.75,.9,.99))+f'  zero-rain (<0.1) area {100*W[A<0.1].sum()/W.sum():.0f} %  area > 8 mm/d {100*W[A>8].sum()/W.sum():.0f} %')
