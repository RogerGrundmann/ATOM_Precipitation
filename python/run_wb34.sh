#!/bin/bash
# 2026-10-04: wb34 = wb33b run to 600 from scratch -- the VERIFICATION against the closing criteria of RAIN-CONV and STORM.
# Stack = working branch (wb27) + ATM_RH_OCEAN=0.82 + ATM_MC_MB_SAT_OCEAN = ATM_MC_MB_SAT_LAND = 0.045 + ATM_RH_OCEAN_ML=1500 + ATM_RH_OCEAN_ML_STRENGTH=0.4.
# cli/atm_oms (-O2, d8193cd + the strength knob = this commit), 8 threads, alone on the machine.
# wb33b at 220: 971.2 mm/a (-0.7 %), r .582, sigma 1.36, bands 2036/515/878/37, land/ocean 561/1134, ocean 0-15 1998, ocean 35-65 1049, land 35-65 420,
# max cell 23.9 mm/d, P/E 1.83.
# RAIN-CONV criteria (at 600): global within 10 % of NASA, ocean 0-15 < 2200, sigma < 2.0, r >= .55, no cell > 50 mm/d.
# STORM criteria (proposed): ocean 35-65 > 500, land 35-65 <= 770, r >= .55.
# PRE-REGISTERED: all eight met at 600. Tropics as wb27 behaved (+0.7 % over 100-600): global 960-1010, r .57-.59, sigma 1.3-1.45, ocean 0-15 1990-2060.
# RISK, the reason for this run: the storm-track rain is drawn from a prescribed initial layer -- at RH_STORM 1.25 the band rained out
# (1847 at iteration 20 -> 859 at 600). 35-65 may DECAY: band 878 -> 700-880 is a pass, below 600 means the layer is being consumed.
# Exit 0, zero NaN through 155 / 357 / 483, max |v| < 4 m/s, converged 1.
set -u; cd "$(dirname "$0")"; rm -f WB34_DONE
t=wb34
mkdir output_$t || { touch WB34_DONE; exit 1; }
[ -e config_$t.xml ] && { touch WB34_DONE; exit 1; }
sed -e "s#output_wb7/#output_$t/#" config_wb7.xml > config_$t.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_oms $(md5sum < ../cli/atm_oms | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_RH_OCEAN=0.82 ATM_MC_MB_SAT_OCEAN=0.045 ATM_MC_MB_SAT_LAND=0.045 ATM_RH_OCEAN_ML=1500 ATM_RH_OCEAN_ML_STRENGTH=0.4"
env OMP_NUM_THREADS=8 $K ../cli/atm_oms config_$t.xml > $t.log 2>&1
echo "== $t  exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
grep -a "RUN CONFIG" $t.log | grep -o "MC_MB_SAT_OCEAN=[^ ]*\|MC_MB_SAT_LAND=[^ ]*\|RH_OCEAN=[^ ]*\|RH_OCEAN_ML=[^ ]*\|RH_OCEAN_ML_STRENGTH=[^ ]*" | sort -u | tr '\n' ' '; echo
echo "-- trajectory (every 100: bands; score)"
grep -a "by |latitude|" $t.log | grep -v MFC | awk 'NR%100==0' | cut -c1-150
grep -a "model .*NASA .*bias" $t.log | awk 'NR%100==0' | cut -c1-150
echo "-- final"
grep -a "by |latitude|" $t.log | grep -v MFC | tail -2; grep -a "model .*NASA .*bias" $t.log | tail -2 | cut -c1-150
grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1
grep -a "P_conv mean\|P_rain mean" $t.log | tail -2
grep -a "max v-component\|max w-component\|max u-component" $t.log | tail -3 | cut -c1-170
grep -a "\[MC-LO\] ocean<30\|\[MC-LO\] land<30" $t.log | tail -2 | cut -c1-330
grep -a "\[MC-RG\]" $t.log | tail -6 | cut -c1-260
tail -1 output_$t/convergence.csv
python3 oceanb.py output_$t/0Ma_smooth_Atm_radial_0_520.vtk
python3 landb.py output_$t/0Ma_smooth_Atm_radial_0_520.vtk | grep -E "Amazon|Congo|global max"
python3 tropdist.py output_$t/0Ma_smooth_Atm_radial_0_520.vtk | head -1
touch WB34_DONE
