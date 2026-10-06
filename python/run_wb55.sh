#!/bin/bash
# 2026-10-06: wb55a / b / c = working branch (wb49 stack) with the ocean RH-vs-SST relation moved, for RC-RESID "rain at 28-29 C 6.51 vs 4.79 mm/d"
# (user: "pursue the next RC-RESID residual, rain at 28-29 C"). Existing knobs only:
#   a ATM_RH_OCEAN_SST_REF 29.5 -> 30.5;   b REF 31;   c REF 30.5 + ATM_RH_OCEAN_SST_MAX 0.014 -> 0.018.
# SCREENING AT nm 60 from scratch, cli/atm_evw (-O2, md5 b88a0004; the evaporation knobs OFF = wb49's physics), 3 x 5 threads (wb54 is running on 8).
# WHY (sst2829.py on wb49): by SST class the model is +18 % at 27-28 C, +40 % at 28-28.5, +33 % at 28.5-29, +16 % at 29-29.5; inside a class the model
# is almost uniform (28-29 C: p10/p90 6.1/7.2 against NASA 2.5/7.2, sd 0.7 vs 1.7) -- its p90 is right and two thirds of the excess sit in cells where
# NASA has < 4 mm/d (Arabian Sea 5.7 vs 2.4, C Pacific at and south of the equator 6.5-6.8 vs 3.4-3.8): subsidence, which no local quantity knows.
# So a uniform reduction can fix the class MEAN only, at the cost of the wettest cells. RH now: 0.82 at >= 29.5 C, -0.007 /K below, at most -0.014.
# With REF 30.5: -0.014 at <= 28.5 C, -0.012 at 28.75, -0.009 at 29.25 (extra -0.005 / -0.007 / -0.007 at 28.25 / 28.75 / 29.25; 27-28 C unchanged).
# Control = wb49 at 60 (= wb47b): 971.8 mm/a, r .609, sigma 1.24, 0-15 1761, ocean 0-15 1788 (NASA 1440); classes 27-28 / 28-29 / 29-31: 4.36 / 6.51 / 7.64
# (NASA 3.71 / 4.79 / 6.55); tropical ocean p50/p90/p99 2.7/6.9/8.3 (NASA 2.5/6.4/8.4), r(model, NASA) over the tropical ocean .708.
# PRE-REGISTERED (0.02 in RH = x1.65 on warm water): a -- 28-29 C 5.5-5.9, 29-31 6.6-7.0, 27-28 within 3 %, global -2 to -3 %, ocean 0-15 1650-1720,
# p99 7.3-7.8, r >= .605;  b -- 28-29 C 5.2-5.6, 29-31 6.2-6.6, global -3 to -4 %, p99 7.0-7.5;  c -- as a at >= 28.5 C, 27-28 C 3.6-4.0, 26-27 C may
# COLLAPSE (cool-water cliff; rainless area above 5 %) -- then c is out.
# RESULT (10:40-10:46, NaN 0; classes 27-28 / 28-29 / 29-31 C mm/d | global | r | sigma | ocean 0-15 | tropical-ocean p50/p90/p99):
#   wb49:            4.36 / 6.51 / 7.64 | -0.7 % | .609 | 1.24 | 1788 | 2.7/6.9/8.3      NASA 3.71 / 4.79 / 6.55, 1440, 2.5/6.4/8.4
#   a REF 30.5:      4.18 / 5.95 / 6.93 | -3.5 % | .599 | 1.18 | 1668 | 2.7/6.3/7.4
#   b REF 31:        4.18 / 5.87 / 6.60 | -4.0 % | .595 | 1.16 | 1643 | 2.7/6.1/7.1
#   c 30.5, MAX .018: 3.64 / 5.85 / 6.92 | -6.4 % | .598 | 1.15 | 1573 | 2.3/6.3/7.4   (26-27 C 2.65 -> 2.23, no collapse; ocean 15-35 613 -> 562)
# The class mean falls less than pre-registered (5.95 vs 5.5-5.9: these cells sit at the mass-flux ceiling, where the humidity lever is weak) and still
# stands 24 % above NASA, while the wettest cells fall below NASA (p99 7.4 vs 8.4), r loses .010 and the global mean 2.8 points. The spread inside a class
# does not change (28-28.5 C p10/p90 5.3/6.0 against NASA 2.1/7.0). No arm is better than wb49.
set -u; cd "$(dirname "$0")"; rm -f WB55_DONE
for t in wb55a wb55b wb55c; do
  mkdir output_$t || { touch WB55_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB55_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>60<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>60<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_evw $(md5sum < ../cli/atm_evw | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_VTK_STRIDE=2"
env OMP_NUM_THREADS=5 $K ATM_RH_OCEAN_SST_REF=30.5 ../cli/atm_evw config_wb55a.xml > wb55a.log 2>&1 &
env OMP_NUM_THREADS=5 $K ATM_RH_OCEAN_SST_REF=31 ../cli/atm_evw config_wb55b.xml > wb55b.log 2>&1 &
env OMP_NUM_THREADS=5 $K ATM_RH_OCEAN_SST_REF=30.5 ATM_RH_OCEAN_SST_MAX=0.018 ../cli/atm_evw config_wb55c.xml > wb55c.log 2>&1 &
wait
for t in wb55a wb55b wb55c; do
  V=output_$t/0Ma_smooth_Atm_radial_0_60.vtk
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)  slices: $(ls output_$t | grep -c radial_0)"
  echo "banner diff vs wb49 (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb49.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "nm=\|output" | head -8
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90
  python3 oceanb.py $V
  ITER=60 python3 pacband.py $t | grep -v "^ \+-\?[0-9]\+:\|ocean by latitude"
  python3 sst2829.py $t 60 | sed -n "2,7p" | cut -c1-150
done
touch WB55_DONE
