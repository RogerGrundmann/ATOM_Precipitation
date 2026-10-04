import sys,numpy as np
# land precipitation by latitude (mm/d, model / NASA), N and S, plus the wettest land cell poleward of 50 deg -- from a radial VTK slice
for f in sys.argv[1:]:
    L=open(f,errors='replace').read().split('\n'); r={}
    for i,l in enumerate(L):
        if l.startswith('SCALARS '):
            n=l.split()[1]
            if n in('Precipitation','Precipitation_NASA','Topography') and n not in r: r[n]=np.array(L[i+2:i+2+65341],float).reshape(181,361)
    P,N,land=r['Precipitation'],r['Precipitation_NASA'],r['Topography']>0; lat=90-np.arange(181)
    s='  land N (model/NASA mm/d):'
    for la in(40,44,48,52,56,60,64,68,72,76):
        m=land[90-la]; s+=f' {la}: {P[90-la][m].mean():.2f}/{N[90-la][m].mean():.2f}'
    print(s); s='  land S:'
    for la in(40,48,68,72,76,80):
        m=land[90+la]; s+=f' {la}: {P[90+la][m].mean():.2f}/{N[90+la][m].mean():.2f}' if m.sum()>3 else ''
    print(s)
    hi=land&(abs(lat)>=50)[:,None]; A=np.where(hi,P,-1); j,k=divmod(A.argmax(),361)
    print(f'  wettest land cell poleward of 50: {A.max():.1f} mm/d at {lat[j]} {k if k<=180 else k-360}; land cells there > 10 mm/d: {(A>10).sum()}')
