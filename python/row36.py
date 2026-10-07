import sys,numpy as np
# STORM-SHAPE, READ-ONLY (2026-10-07): ocean rows 34-38 deg rain ~565 mm/a against NASA 1000-1200 and sit on a cliff of the storm-track RH factor
# (1.030 -> 424, 1.037 -> 565, 1.079 -> 1887 mm/a). Is the rain not generated, or generated and evaporated below the cloud?
# usage: row36.py <run> [iter] [run2 ...]   (1) ocean rows from the level-0 slice; (2) the 87E section (S Indian ocean): rain flux profile per column
runs=[a for a in sys.argv[1:] if not a.isdigit()]; it=[a for a in sys.argv[1:] if a.isdigit()]; it=it[0] if it else '520'
def rd(f,n,shape,want=None):
    L=open(f,errors='replace').read().split('\n'); o={}
    for i,l in enumerate(L):
        if l.startswith('SCALARS '):
            nm=l.split()[1]
            if nm not in o and (want is None or nm in want): o[nm]=np.array(L[i+2:i+2+n],float).reshape(shape)
    return o
lat=90-np.arange(181)
def fac(a,F=1.13,W=17.0,C=55.0): return 1+(F-1)*np.exp(-((a-C)/W)**2)
def mlw(a):
    x=np.clip((a-30)/10,0,1); return 0.4*x*x*(3-2*x)
for run in runs:
    s=rd(f'output_{run}/0Ma_smooth_Atm_radial_0_{it}.vtk',65341,(181,361),('Topography','Precipitation','PrecipitationConv','Precipitation_NASA','Temperature','HumidityRel','PrecipitableWater','Evaporation'))
    LA=lat[:,None]*np.ones((1,361)); w=np.cos(np.radians(LA)); w[:,360]=0; oc=s['Topography']==0
    g=lambda A,m:(A*w)[m].sum()/w[m].sum()
    print(f'\n== {run} {it}: OCEAN rows (S / N).  row | T C | PW mm | RH0 % | model mm/a | of it convective | NASA | model/NASA   [storm factor wb66, ML weight at row centre]')
    for a in range(22,54,4):
        o=[]
        for sg in(-1,1):
            m=oc&(sg*LA>=a)&(sg*LA<a+4)
            o.append(f"{g(s['Temperature'],m):5.1f} {g(s['PrecipitableWater'],m):5.1f} {g(s['HumidityRel'],m):5.1f} {365*g(s['Precipitation'],m):5.0f} {365*g(s['PrecipitationConv'],m):5.0f} {365*g(s['Precipitation_NASA'],m):5.0f} {g(s['Precipitation'],m)/g(s['Precipitation_NASA'],m):4.2f}")
        print(f'  {a}-{a+4} | '+' | '.join(o)+f'   [{fac(a+2):.3f}, {mlw(a+2):.2f}]')
    # inside the row 34-38: how is the rain distributed over the cells?
    for sg,nm in((-1,'S'),(1,'N')):
        m=oc&(sg*LA>=34)&(sg*LA<38); p=365*s['Precipitation'][m]; n=365*s['Precipitation_NASA'][m]; c=365*s['PrecipitationConv'][m]; T=s['Temperature'][m]
        print(f'  row 34-38{nm}: cells {m.sum()}, model p10/p50/p90 {np.percentile(p,10):.0f}/{np.percentile(p,50):.0f}/{np.percentile(p,90):.0f}, NASA {np.percentile(n,10):.0f}/{np.percentile(n,50):.0f}/{np.percentile(n,90):.0f}, cells < 100 mm/a {100*(p<100).mean():.0f} %, cells with conv > 100 {100*(c>100).mean():.0f} %')
        for lo,hi in((10,16),(16,18),(18,20),(20,22),(22,26)):
            q=(T>=lo)&(T<hi)
            if q.sum()>5: print(f'      T {lo}-{hi} C: {100*q.mean():3.0f} % of cells, model {p[q].mean():5.0f} (conv {c[q].mean():4.0f}), NASA {n[q].mean():5.0f}')
    z=rd(f'output_{run}/0Ma_smooth_Atm_zonal_87_{it}.vtk',7421,(41,181),('Topography','height','HumidityRel','CloudWater','CloudIce','PrecipitationRain','PrecipitationSnow','PrecipitationConv','Temperature','S_r','S_c_c','WaterVapour'))
    h=1000*z['height'][:,0]; F=z['PrecipitationRain']+z['PrecipitationSnow']
    print(f'  87E section, ocean columns: lat | rain flux max [units of the slice] at height m | flux at the sea | reaching % | conv at sea | cloud water max g/kg at m | lowest cloud level m (qc>0.005) | RH % at 0 / 440 / 940 / 1370 / 2160 / 3320 m')
    ih=[int(np.argmin(abs(h-x))) for x in(0,440,940,1370,2160,3320)]
    for la in(-28,-30,-32,-34,-35,-36,-37,-38,-40,-42,-46,-50):
        j=90-la
        if z['Topography'][0,j]!=0: print(f'   {la}: land'); continue
        f=F[:,j]; im=int(np.argmax(f)); qc=z['CloudWater'][:,j]+z['CloudIce'][:,j]; ic=int(np.argmax(qc)); lowc=np.where(qc>0.005)[0]
        print(f"   {la} | {f[im]:8.4g} at {h[im]:5.0f} | {f[0]:8.4g} | {100*f[0]/f[im] if f[im]>0 else 0:5.1f} | {z['PrecipitationConv'][0,j]:7.3g} | {qc[ic]:.4f} at {h[ic]:5.0f} | {h[lowc[0]] if len(lowc) else -1:5.0f} | "+' / '.join(f"{z['HumidityRel'][i,j]:.0f}" for i in ih))
