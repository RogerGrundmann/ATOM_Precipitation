#!/usr/bin/env python3
# Score the bisection arms (run_bisect.sh) against each other and def600. One Precip print per iteration.
import re, sys
def parse(f):
    P=[];R=[];S=[];PW=[];B=[];gp=None;ev=None
    for l in open(f, errors='replace'):
        m=re.search(r'Precip mean = ([\d.e+-]+)',l);                    P.append(float(m.group(1))) if m else None
        m=re.search(r'pattern r = ([+-][\d.]+).*sigma model/NASA = ([\d.]+)',l)
        if m: R.append(float(m.group(1))); S.append(float(m.group(2)))
        m=re.search(r'precipitable water average =\s+([\d.]+).*NASA',l); PW.append(float(m.group(1))) if m else None
        if 'by |latitude|' in l: B.append('/'.join(x.split('.')[0] for x in re.findall(r'(\d+\.\d) /',l)))
        m=re.search(r'g_p generated ([\d.e+-]+).*P_conv\(ground\) ([\d.e+-]+)',l)
        if m: gp=(float(m.group(1)),float(m.group(2)))
        m=re.search(r'\[CWB\] evaporation\s+([-\d.e+]+)',l)
        if m: ev=float(m.group(1))
    return P,R,S,PW,B,gp,ev
arms=sys.argv[1:] or ['def600','bis_ctl','bis_wc','bis_sf','bis_sp','bis_ms','bis_hs','bis_rt','bis_old']
print(f'{"arm":8s} {"P20":>6s} {"P40":>6s} {"P60":>6s} {"P80":>6s} {"P100":>6s} {"slope40-100":>11s}  {"r":>6s} {"sig":>5s} {"PW":>5s}  bands@last          g_p/P_conv      CWB evap')
for a in arms:
    try: P,R,S,PW,B,gp,ev=parse(a+'.log')
    except FileNotFoundError: continue
    n=len(P)-1
    g=lambda i: f'{P[i]:6.0f}' if i<len(P) else '     -'
    last=min(n,100)
    sl=(P[last]-P[40])/(last-40) if last>40 else float('nan')
    print(f'{a:8s} {g(20)} {g(40)} {g(60)} {g(80)} {g(100)} {sl:11.2f}  {R[last] if last<len(R) else 0:+.3f} {S[last] if last<len(S) else 0:5.2f} {PW[min(last,len(PW)-1)] if PW else 0:5.1f}  {B[last] if last<len(B) else "":18s} {str(gp):15s} {ev}')
