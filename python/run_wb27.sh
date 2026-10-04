#!/bin/bash
# 2026-10-04: wb27 = wb25b run to 600 from scratch (user) -- the VERIFICATION of the screening stack. wb25b = working branch + ML parcel stack +
# ATM_RH_LAND=2 + ATM_RH_LAND_QCAP=1 + ATM_MC_GATE_BLEND=3 + ATM_MC_GP_AREA=1 + ATM_RH_OCEAN=0.80 + ATM_MC_T_ADD_LAND=2.3.
# cli/atm_g3 (-O2, = cli/atm at 759f496), 8 threads (same count as the nm 220 screen), alone on the machine.
# wb25b at 220: 912 mm/a (-6.8 %), r .562, sigma 2.85, bands 3001/290/204/4.5, land/ocean 443/1098, land 0-15 1641, Amazon 8.25, Congo 3.12 mm/d,
# P/E 1.56, wettest cell 38 mm/d (ocean 4S 145E), max |v| 2.5 m/s. wb14 (the last 600 verification): 849 at 200 -> 864 at 600 (+1.8 %).
# PRE-REGISTERED: within ~3 % of the 220 values -- global 900-960, r .55-.575, sigma 2.8-3.0, ocean 1085-1135, land 430-475, land 0-15 1600-1750
# (read at 520, the last radial slice), Amazon 8.0-8.8, Congo 3.0-3.4, bands 0-15 2950-3150, 15-35 280-310, 35-65 200-215, 65-90 4.5;
# no cell > 50 mm/d; exit 0, zero NaN through 155 / 357 / 483; max |v| < 4 m/s (no seam mode); converged = 1.
# RISK: a slow drift past 220 (the mean still moving at 600), as the pre-EVAP_LIMIT branch had; GATE_BLEND=3 has never been run past 220.
set -u; cd "$(dirname "$0")"; rm -f WB27_DONE
t=wb27
mkdir output_$t || { touch WB27_DONE; exit 1; }
[ -e config_$t.xml ] && { touch WB27_DONE; exit 1; }
sed -e "s#output_wb7/#output_$t/#" config_wb7.xml > config_$t.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_g3 $(md5sum < ../cli/atm_g3 | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_MC_ML_PARCEL=1 ATM_MC_T_ADD=0.15 ATM_MC_Q_ADD=4.5e-4 ATM_MC_ML_LCL=1 ATM_RH_LAND=2 ATM_MC_GATE_BLEND=3 ATM_MC_GP_AREA=1 ATM_RH_OCEAN=0.80 ATM_MC_T_ADD_LAND=2.3 ATM_RH_LAND_QCAP=1"
env OMP_NUM_THREADS=8 $K ../cli/atm_g3 config_$t.xml > $t.log 2>&1
echo "== $t  exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
grep -a "RUN CONFIG" $t.log | grep -o "MC_GP_AREA=[^ ]*\|MC_GATE_BLEND=[^ ]*\|RH_OCEAN=[^ ]*\|MC_T_ADD_LAND=[^ ]*" | sort -u | tr '\n' ' '; echo
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
python3 landb.py output_$t/0Ma_smooth_Atm_radial_0_520.vtk
touch WB27_DONE
