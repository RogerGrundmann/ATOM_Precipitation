#!/bin/bash
# 2026-10-07: wb67a / b / c = working branch (wb66 stack, 54 knobs) with the ATM_RH_OCEAN_ML latitude taper moved EQUATORWARD (new knobs
# ATM_RH_OCEAN_ML_LAT0 / _LAT1, shipped 30 / 40; strength 0.4 and depth 1500 m unchanged):   a 26 / 38;   b 22 / 38;   c 18 / 36.
# SCREENING AT nm 60 from scratch, cli/atm_mlt (-O2, HEAD 0cb20ae + the two knobs, md5 44f6333f), 3 x 6 threads.
# WHY (row36.py on wb66): ocean rows 22-38 deg rain at half of NASA in both hemispheres (S: 424 / 443 / 462 / 564 against 728 / 841 / 950 / 1001), all
# stratiform. At 87E the column makes 4.4-6.7 mm/d at ~2.2 km at 28-34S and 19-26 % reaches the sea; 43 % at 36S, 66 % at 40S -- along the taper.
# Mixed-layer weight at 24 / 28 / 32 / 36 / 40 deg: now 0 / 0 / .04 / .26 / .40;  a 0 / .03 / .20 / .37 / .40;  b .02 / .13 / .27 / .38 / .40;  c .10 / .23 / .35 / .40 / .40.
# Control = wb66 (600; wb62b/wb65a at 60 agree to 0.3 %): 970.1 mm/a (-0.8 %), r .631, sigma 1.16, ocean 15-35 693 (NASA 809), ocean 35-65 1070 (1107),
# ocean rows S / N 22-26 424 / 697 (conv 110 / 477), 26-30 443 / 527, 30-34 462 / 478, 34-38 564 / 574, 38-42 1074 / 1075; E Austral 3.20 mm/d.
# PRE-REGISTERED (a cliff, so wide; southern rows, mm/a):  a -- 26-30 450-520, 30-34 560-750, 34-38 700-950;   b -- 22-26 430-480, 26-30 520-700,
# 30-34 750-1050, 34-38 750-1000;   c -- 22-26 480-650, 26-30 700-1000, 30-34 950-1350, 34-38 800-1050.  Rows poleward of 40 deg within 2 % of wb66.
# Global mean UP by 2-3 % (a), 3-6 % (b), 5-9 % (c): the rows hold ~18 % of the globe and nothing else comes down.
# USABLE if r >= .631, no ocean row 22-38 above 1.25x NASA, convective rain in rows 22-30 up by < 20 %, 22-26N <= 1.2x NASA (765), E Austral <= 3.3.
set -u; cd "$(dirname "$0")"; rm -f WB67_DONE
for t in wb67a wb67b wb67c; do
  mkdir output_$t || { touch WB67_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB67_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>60<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>60<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_mlt $(md5sum < ../cli/atm_mlt | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_VTK_STRIDE=2"
env OMP_NUM_THREADS=6 $K ATM_RH_OCEAN_ML_LAT0=26 ATM_RH_OCEAN_ML_LAT1=38 ../cli/atm_mlt config_wb67a.xml > wb67a.log 2>&1 &
env OMP_NUM_THREADS=6 $K ATM_RH_OCEAN_ML_LAT0=22 ATM_RH_OCEAN_ML_LAT1=38 ../cli/atm_mlt config_wb67b.xml > wb67b.log 2>&1 &
env OMP_NUM_THREADS=6 $K ATM_RH_OCEAN_ML_LAT0=18 ATM_RH_OCEAN_ML_LAT1=36 ../cli/atm_mlt config_wb67c.xml > wb67c.log 2>&1 &
wait
for t in wb67a wb67b wb67c; do
  V=output_$t/0Ma_smooth_Atm_radial_0_60.vtk
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)  slices: $(ls output_$t | grep -c radial_0)"
  echo "banner diff vs wb65a (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb65a.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "output" | head -8
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90
  python3 oceanb.py $V
  python3 row36.py $t 60 | grep -v "T [0-9]" | cut -c1-230
  python3 landb.py $V | grep "E Austral\|S China\|SE US\|S Brazil\|Arabia\|Sahara\|Madagascar\|global max"
done
touch WB67_DONE
