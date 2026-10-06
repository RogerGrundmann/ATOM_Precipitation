import sys,numpy as np
# RC-RESID, READ-ONLY (2026-10-06): the 65-90 deg band rains 220 mm/a against NASA 364 (ocean 190 / 428).  Where, and what does the column look like?
# usage: polar.py <run> [iter]
run=sys.argv[1]; it=sys.argv[2] if len(sys.argv)>2 else '520'
def rd(f,want,n0,shape):
    L=open(f,errors='replace').read().split('\n'); o={}
    for n,l in enumerate(L):
        if l.startswith('SCALARS '):
            nm=l.split()[1]
            if nm in want and nm not in o: o[nm]=np.array(L[n+2:n+2+n0],float).reshape(shape)
    return o
s=rd(f'output_{run}/0Ma_smooth_Atm_radial_0_{it}.vtk',('Topography','Landscape','Precipitation','PrecipitationSnow','PrecipitationRain','Precipitation_NASA','Temperature','HumidityRel','PressureStatic','Evaporation','PrecipitableWater'),65341,(181,361))
lat=90-np.arange(181); LA=lat[:,None]*np.ones((1,361)); w=np.cos(np.radians(LA)); w[:,360]=0; land=s['Topography']>0
P,N,T,RH,ps,PW,H,E=(s[x] for x in('Precipitation','Precipitation_NASA','Temperature','HumidityRel','PressureStatic','PrecipitableWater','Landscape','Evaporation'))
g=lambda A,m:(A*w)[m].sum()/w[m].sum()
print(f'{run} {it}: rows of latitude, both hemispheres.   lat | surface | area % of 65-90 | model mm/a | NASA mm/a | deficit share % | T sfc C | RH sfc % | p sfc hPa | ground m | PW mm | E mm/a')
m65=np.abs(LA)>=65; tot=((N-P)*w)[m65].sum()
for a,b in[(65,70),(70,75),(75,80),(80,90),(-70,-65),(-75,-70),(-80,-75),(-90,-80)]:
    for lab,mm in(('ocean',~land),('land ',land)):
        m=(LA>=a)&(LA<b)&mm if a>0 else (LA>a)&(LA<=b)&mm
        if w[m].sum()/w[m65].sum()<0.01: continue
        print(f'  {a:3d}..{b:3d} | {lab} | {100*w[m].sum()/w[m65].sum():4.1f} | {365*g(P,m):5.0f} | {365*g(N,m):5.0f} | {100*((N-P)*w)[m].sum()/tot:5.1f} | {g(T,m):6.1f} | {g(RH,m):4.1f} | {g(ps,m):5.0f} | {g(H,m):5.0f} | {g(PW,m):4.1f} | {365*g(E,m):4.0f}')
for lab,m in(('N ocean',m65&~land&(LA>0)),('N land',m65&land&(LA>0)),('S ocean',m65&~land&(LA<0)),('S land',m65&land&(LA<0))):
    print(f'  {lab}: area {100*w[m].sum()/w[m65].sum():4.1f} %  model {365*g(P,m):4.0f}  NASA {365*g(N,m):4.0f}  share of the deficit {100*((N-P)*w)[m].sum()/tot:4.0f} %')
# ---- column profiles on the 87E section
z=rd(f'output_{run}/0Ma_smooth_Atm_zonal_87_{it}.vtk',('Temperature','WaterVapour','CloudWater','CloudIce','PressureStatic','height','Topography','Precipitation','PrecipitationRain','PrecipitationSnow','S_v','S_r','S_s','HumidityRel'),7421,(41,181))
Tz,q,p,h,topo=z['Temperature'],z['WaterVapour'],z['PressureStatic'],z['height']*1000,z['Topography']
es=lambda Tc: 6.1078*np.exp(17.2694*Tc/(Tc+237.29)); esi=lambda Tc: 6.1078*np.exp(21.8746*Tc/(Tc+265.49))
def hcrit(p): x=np.maximum(p/1000.,0.55); return np.clip(1-(0.7/(0.55*0.45))*x*(1-x),0,1)
def i0f(j): c=np.where(topo[:,j]>0)[0]; return (c.max()+1) if len(c) else 0
G={'ARCTIC OCEAN 78-88N':[j for j in range(181) if 78<=lat[j]<=88],'SIBERIA 66-74N':[j for j in range(181) if 66<=lat[j]<=74],'S OCEAN 60-66S':[j for j in range(181) if -66<=lat[j]<=-60],'ANTARCTICA 70-80S':[j for j in range(181) if -80<=lat[j]<=-70]}
for nm,js in G.items():
    js=[j for j in js]; i0s=[i0f(j) for j in js]
    print(f'\n{nm} (87E): ground {np.mean([h[i0,j] for i0,j in zip(i0s,js)]):.0f} m, surface precip {365*np.mean([z["Precipitation"][0,j] for j in js]):.0f} mm/a')
    print('  z above gnd |   p  |   T   | RH water  RH ice  H_crit  RH-Hc | q_v    q_c    q_i  g/kg | S_v   S_r   S_s (x1e-9) | rain  snow flux mm/d')
    for zt in[50,250,600,1000,1500,2200,3000,4000,5000,6500]:
        v=[]
        for i0,j in zip(i0s,js):
            i=min(range(i0,41),key=lambda i:abs(h[i,j]-h[i0,j]-zt)); e=q[i,j]*1e-3*p[i,j]/0.622
            v.append((p[i,j],Tz[i,j],e/es(Tz[i,j]),e/esi(Tz[i,j]),hcrit(p[i,j]),q[i,j],z['CloudWater'][i,j],z['CloudIce'][i,j],z['S_v'][i,j]*1e9,z['S_r'][i,j]*1e9,z['S_s'][i,j]*1e9,z['PrecipitationRain'][i,j],z['PrecipitationSnow'][i,j]))
        v=np.mean(v,axis=0)
        print(f'  {zt:6d} m   | {v[0]:4.0f} | {v[1]:5.1f} | {v[2]:.2f}  {v[3]:.2f}  {v[4]:.2f}  {v[2]-v[4]:+.2f} | {v[5]:.3f} {v[6]:.4f} {v[7]:.4f} | {v[8]:7.1f} {v[9]:7.1f} {v[10]:7.1f} | {v[11]:.3g} {v[12]:.3g}')
