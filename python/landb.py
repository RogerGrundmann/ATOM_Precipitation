import sys,numpy as np
for f in sys.argv[1:]:
    L=open(f,errors='replace').read().split('\n'); o={}
    for n,l in enumerate(L):
        if l.startswith('SCALARS '):
            nm=l.split()[1]
            if nm in('Precipitation','PrecipitationConv','Precipitation_NASA','Topography') and nm not in o: o[nm]=np.array(L[n+2:n+2+65341],float).reshape(181,361)
    P,C,N,T=o['Precipitation'],o['PrecipitationConv'],o['Precipitation_NASA'],o['Topography']
    lat=90-np.arange(181); lon=np.where(np.arange(361)<=180,np.arange(361),np.arange(361)-360); w=np.cos(np.radians(lat))[:,None]*np.ones((1,361)); land=T>0
    print(f)
    for a,b in[(0,15),(15,35),(35,65),(65,90)]:
        m=land&((abs(lat)>=a)&(abs(lat)<b))[:,None]; g=lambda A:(A*w)[m].sum()/w[m].sum()*365
        print(f'  land {a:2d}-{b}: model {g(P):6.0f}  NASA {g(N):6.0f}')
    reg={'India':(10,28,72,88),'S China':(20,32,102,120),'SE US':(27,36,-95,-78),'E Austral':(-35,-18,145,153),'SE Africa':(-28,-15,28,38),'S Brazil':(-30,-18,-55,-42),
         'Sahara':(18,30,-10,25),'Arabia':(16,28,42,56),'Oman':(17,24,52,59),'W Austral':(-30,-20,115,128),'Namib':(-28,-18,12,18),'Amazon':(-10,2,-70,-52),'Congo':(-5,5,15,28),'Mexico plt':(18,30,-108,-96),'Highveld':(-32,-22,24,32)}
    for nm,(a,b,c,d) in reg.items():
        m=land&((lat>=a)&(lat<=b))[:,None]&((lon>=c)&(lon<=d))[None,:]; g=lambda A:(A*w)[m].sum()/w[m].sum()
        print(f'  {nm:10s} model {g(P):6.2f} mm/d   NASA {g(N):5.2f}   max cell {P[m].max():6.1f}')
    j,k=divmod(P.argmax(),361); print('  global max',P.max().round(1),'mm/d at',lat[j],lon[k],'land' if T[j,k]>0 else 'ocean')
