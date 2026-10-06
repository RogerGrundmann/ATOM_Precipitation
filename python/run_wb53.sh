#!/bin/bash
# 2026-10-06: wb53a / b / c = working branch (wb49 stack) + the evaporation repairs for RC-RESID P/E 1.80 (user: "pursue RC-RESID, start with P/E"):
#   a ATM_EVAP_WIND=1 (wind speed |V| instead of |V|/sqrt(3));  b + ATM_EVAP_GUST=3;  c + ATM_LAND_BUCKET=150 (land evaporates; Manabe bucket).
# SCREENING AT nm 60 from scratch, cli/atm_evw (-O2, HEAD 1170bc8 + the two knobs), 3 x 6 threads, ATM_VTK_STRIDE=2; -O0 byte check (run_vevw.sh) alongside.
# WHY (evapb.py on wb49): E 539 mm/a against P 970. Ocean E 752 (Earth ~1300) because the Meyer formula is given sqrt((u2+v2+w2)/3) = 3.0-3.6 m/s in the
# tropics where the surface wind is 5.2-6.2, and 0.9 m/s at 25-35 deg (annual-mean vector wind; Earth's evaporation maximum is there); land E is 0.
# Control = wb49 at 60 (= wb47b): 971.8 mm/a, r .609, sigma 1.24, land 15-35 217, land / ocean 651 / 1099, E 539, P/E 1.80.
# PRE-REGISTERED (at 60): E a 680-710, b 745-785, c 830-900 (P/E 1.40 / 1.27 / 1.12); rain a and b within 0.5 % of wb49 in every band and region, r .609
# (E is 0.002 mm in 120 s against a 32 mm column); c: land rain may rise (land evaporation moistens the layer the convection reads) -- land 0-15 within
# 5 %, deserts < 0.3 mm/d, no cell > 14 mm/d, else the bucket is not usable on this branch.
# RESULT (10:16-10:22, NaN 0; E global / ocean / land mm/a, P/E; rain):
#   wb49: 539 / 752 / 0, 1.80
#   a WIND:                 706 /  986 /   0, P/E 1.38;  rain 971.9, r .609, every band and region as wb49 (to 0.1 %)
#   b WIND + GUST 3:        776 / 1083 /   0, P/E 1.25;  rain 972.0, identical
#   c + LAND_BUCKET 150:    982 / 1083 / 728, P/E 0.99;  rain 972.0, land 650.7 (wb49 650.5), SE US 1.19 -> 1.21, nothing else moves
# As pre-registered for a and b (E 706 / 776 vs 680-710 / 745-785). c: land rain did NOT rise, but land E 728 exceeds land P 653 (Earth land E ~500,
# E/P ~0.6): the bucket starts full wherever NASA rain > 750 mm/a and then evaporates at the open-water rate, so P/E 0.99 is partly a compensation
# (ocean 1083 against Earth ~1300, land 728 against ~500). Ocean E by band with b: 1570 / 1150 / 780 / 184 (Earth ~1450 / ~1650 / ~950 / ~250).
set -u; cd "$(dirname "$0")"; rm -f WB53_DONE
for t in wb53a wb53b wb53c; do
  mkdir output_$t || { touch WB53_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB53_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>60<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>60<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_evw $(md5sum < ../cli/atm_evw | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_VTK_STRIDE=2 ATM_EVAP_WIND=1"
env OMP_NUM_THREADS=6 $K ../cli/atm_evw config_wb53a.xml > wb53a.log 2>&1 &
env OMP_NUM_THREADS=6 $K ATM_EVAP_GUST=3 ../cli/atm_evw config_wb53b.xml > wb53b.log 2>&1 &
env OMP_NUM_THREADS=6 $K ATM_EVAP_GUST=3 ATM_LAND_BUCKET=150 ../cli/atm_evw config_wb53c.xml > wb53c.log 2>&1 &
wait
for t in wb53a wb53b wb53c; do
  V=output_$t/0Ma_smooth_Atm_radial_0_60.vtk
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)  slices: $(ls output_$t | grep -c radial_0)"
  grep -a "water budget closure" $t.log | tail -1 | cut -c1-110; grep -a "LAND-BUCKET\|land bucket" $t.log | tail -1 | cut -c1-200
  python3 evapb.py $t | sed -n "2,10p" | cut -c1-120
  echo "banner diff vs wb49 (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb49.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "nm=\|output" | head -8
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90
  python3 oceanb.py $V
  python3 landb.py $V | grep -v "Oman\|Namib\|land \|output"
done
touch WB53_DONE
