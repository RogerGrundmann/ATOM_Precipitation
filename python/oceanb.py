import sys,numpy as np
# ocean / land precipitation by |latitude| band from a radial VTK slice (mm/a), model vs NASA; companion of landb.py
for f in sys.argv[1:]:
    L=open(f,errors='replace').read().split('\n'); o={}
    for n,l in enumerate(L):
        if l.startswith('SCALARS '):
            nm=l.split()[1]
            if nm in('Precipitation','PrecipitationConv','Precipitation_NASA','Topography') and nm not in o: o[nm]=np.array(L[n+2:n+2+65341],float).reshape(181,361)
    P,C,N,T=o['Precipitation'],o['PrecipitationConv'],o['Precipitation_NASA'],o['Topography']
    lat=90-np.arange(181); w=np.cos(np.radians(lat))[:,None]*np.ones((1,361)); land=T>0
    for a,b in[(0,15),(15,35),(35,65),(65,90)]:
        B=((abs(lat)>=a)&(abs(lat)<b))[:,None]; g=lambda A,m:(A*w)[m].sum()/w[m].sum()*365
        print(f'  {a:2d}-{b}: ocean {g(P,B&~land):6.0f} (conv {g(C,B&~land):6.0f}) NASA {g(N,B&~land):6.0f} | land {g(P,B&land):6.0f} (conv {g(C,B&land):6.0f}) NASA {g(N,B&land):6.0f}')
