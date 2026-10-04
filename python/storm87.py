import sys,numpy as np
f=sys.argv[1]; L=open(f,errors='replace').read().split('\n'); o={}
want=['Temperature','WaterVapour','CloudWater','CloudIce','PressureStatic','height','Topography','Precipitation','PrecipitationRain','PrecipitationSnow','S_c_c','S_v','S_c','S_i','S_r','S_s','r_humid','Precipitation_NASA']
for n,l in enumerate(L):
    if l.startswith('SCALARS '):
        nm=l.split()[1]
        if nm in want and nm not in o: o[nm]=np.array(L[n+2:n+2+7421],float).reshape(41,181)
T=o['Temperature'];q=o['WaterVapour'];p=o['PressureStatic'];h=o['height']*1000;topo=o['Topography']
lat=90-np.arange(181)
def es(Tc): return 6.112*np.exp(17.67*Tc/(Tc+243.5))
def esi(Tc): return 6.112*np.exp(22.46*Tc/(Tc+272.62))
def hcrit(p):
    x=np.maximum(p/1000.,0.55); return np.clip(1-(0.7/(0.55*0.45))*x*(1-x),0,1)
print('units: P_rain max',o['PrecipitationRain'].max(),'S_r max',abs(o['S_r']).max(),'cloud max',o['CloudWater'].max())
def col(j):
    c=np.where(topo[:,j]>0)[0]; i0=(c.max()+1) if len(c) else 0; return i0
groups={'OCEAN 40-56S (87E)':[j for j in range(181) if -56<=lat[j]<=-40],'LAND 40-52N (87E)':[j for j in range(181) if 40<=lat[j]<=52],'OCEAN 36S':[90+36],'LAND 56-64N':[j for j in range(181) if 56<=lat[j]<=64]}
for nm,js in groups.items():
    i0s=[col(j) for j in js]; print(f'\n{nm}: ground levels {min(i0s)}..{max(i0s)}, ground height {np.mean([h[i0,j] for i0,j in zip(i0s,js)]):.0f} m, surface P {np.mean([o["Precipitation"][0,j] for j in js]):.2f} mm/d  rain {np.mean([o["PrecipitationRain"][i0,j] for i0,j in zip(i0s,js)]):.3g} snow {np.mean([o["PrecipitationSnow"][i0,j] for i0,j in zip(i0s,js)]):.3g}')
    print('  z above gnd |  p   |   T   |  RH   H_crit  RH-Hc | q_c    q_i  g/kg | S_c_c    S_r     S_s     S_v (x1e-9) | P_rain  P_snow (flux)')
    for zt in [50,250,600,1000,1500,2200,3000,4000,5000,6500,8000]:
        v=[]
        for i0,j in zip(i0s,js):
            i=min(range(i0,41),key=lambda i:abs(h[i,j]-h[i0,j]-zt))
            qs=0.622*es(T[i,j])/(p[i,j]-es(T[i,j]))*1e3
            v.append((p[i,j],T[i,j],q[i,j]/qs,hcrit(p[i,j]),o['CloudWater'][i,j],o['CloudIce'][i,j],o['S_c_c'][i,j]*1e9,o['S_r'][i,j]*1e9,o['S_s'][i,j]*1e9,o['S_v'][i,j]*1e9,o['PrecipitationRain'][i,j],o['PrecipitationSnow'][i,j]))
        v=np.mean(v,axis=0)
        print(f'  {zt:6d} m   | {v[0]:4.0f} | {v[1]:5.1f} | {v[2]:.2f}  {v[3]:.2f}  {v[2]-v[3]:+.2f} | {v[4]:.4f} {v[5]:.4f} | {v[6]:7.2f} {v[7]:7.2f} {v[8]:7.2f} {v[9]:7.2f} | {v[10]:.3g} {v[11]:.3g}')
