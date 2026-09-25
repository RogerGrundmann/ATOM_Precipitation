#!/usr/bin/env python3
# ATM_CWB_BANDS consistency: for the LAST [CWB] window in a log, the four band-mean rates weighted by the
# bands' cos-lat weights must reproduce the stage's global total. Weights replicate ColumnWaterBudget.
import sys, re, math
f = sys.argv[1]; jm, km = 181, 361
lat_w = lambda j: math.cos((90 - j) * math.pi / 180) if j <= 90 else math.cos((j - 90) * math.pi / 180)
def band(j):
    a = abs(90.0 - j * 180.0 / (jm - 1)); return 0 if a < 15 else 1 if a < 35 else 2 if a < 65 else 3
wb = [0.0]*4
for j in range(jm): wb[band(j)] += lat_w(j) * km
W = sum(wb)
L = open(f, errors='replace').read().split('\n')
ib = max(i for i, l in enumerate(L) if '[CWB-BANDS] stage' in l)
ig = max(i for i, l in enumerate(L[:ib]) if '[CWB] stage' in l)
glob = {}
for l in L[ig+1:ib]:
    m = re.match(r'\s*AGCM: \[CWB\] (\S.*?)\s{2,}([-+.\de]+)\s+([-+.\de]+)\s+([-+.\de]+)\s*$', l)
    if m: glob[m.group(1).strip()] = float(m.group(2))
worst = 0.0
print(f'{"stage":22s} {"global":>12s} {"from bands":>12s}   0-15 / 15-35 / 35-65 / 65-90')
for l in L[ib+1:]:
    m = re.match(r'\s*AGCM: \[CWB-BANDS\] (\S.*?)\s{2,}([-+.\de]+)\s+([-+.\de]+)\s+([-+.\de]+)\s+([-+.\de]+)\s*$', l)
    if not m: break
    name = m.group(1).strip(); v = [float(m.group(k)) for k in range(2, 6)]
    rec = sum(wb[b] * v[b] for b in range(4)) / W
    g = glob.get(name)
    if g is not None:
        scale = max(abs(g), 1e-6 * max(abs(x) for x in v + [1.0]))
        worst = max(worst, abs(rec - g) / max(abs(g), 1.0))
    print(f'{name:22s} {g if g is not None else float("nan"):12.4e} {rec:12.4e}   ' + ' / '.join(f'{x:.3e}' for x in v))
print(f'worst relative mismatch (vs |global|, floor 1 mm/a): {worst:.2e}   -> {"PASS" if worst < 1e-3 else "FAIL"}  (print precision 4 digits)')
