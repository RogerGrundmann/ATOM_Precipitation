#!/bin/bash
# 2026-10-03: wb14 = wb10a stack + ATM_MC_GATE_BLEND=1 + ATM_RH_LAND_QCAP=1 at ATM_RH_OCEAN 0.805 (user). 600 from scratch, cli/atm_gq (-O2,
# = 9aba340 source), 8 threads. Two changes against wb13 (cap, RH_OCEAN between the wb13 arms).
# wb13a (0.82): 1642 mm/a, r .445, sigma 5.78, bands 5408/657/238/4.6, land/ocean 635/2041.  wb13b (0.79): 347, r .289, 2.30, 1108/9.1/198/4.6, 265/380.
# PRE-REGISTERED: (1) the Horn of Africa is dry: P_conv 0 in 1-5N 40-44E, global max P_conv an ocean cell < 60 mm/d (wb13a 129 mm/d at 3N 36E);
# (2) global 700-1000 mm/a (geometric mean of the wb13 arms is 755; CAPE is convex in RH so rather below the midpoint 995), ocean 800-1300,
# 0-15 band 2400-3400; (3) land BELOW the wb13 interpolation (~450) because the cap dries hot land: 300-450; r >= .40.
set -u; cd "$(dirname "$0")"; rm -f WB14_DONE
mkdir output_wb14 || { touch WB14_DONE; exit 1; }
[ -e config_wb14.xml ] && { touch WB14_DONE; exit 1; }
sed "s#output_wb7/#output_wb14/#" config_wb7.xml > config_wb14.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_gq $(md5sum < ../cli/atm_gq | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_RH_OCEAN=0.805 ATM_MC_ML_PARCEL=1 ATM_MC_T_ADD=0.15 ATM_MC_Q_ADD=4.5e-4 ATM_MC_ML_LCL=1 ATM_RH_LAND=2 \
    ATM_MC_GATE_BLEND=1 ATM_RH_LAND_QCAP=1 ../cli/atm_gq config_wb14.xml > wb14.log 2>&1
t=wb14
echo "== $t  exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
grep -a "RH-LAND-QCAP" $t.log | head -1
grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1
grep -a "P_conv mean" $t.log | tail -1
grep -a "\[MC-GATE\]" $t.log | tail -2 | cut -c1-260
grep -a "\[MC-LO\] ocean<30\|\[MC-LO\] land<30" $t.log | tail -2 | cut -c1-330
grep -a "\[MC-RG\]" $t.log | tail -6 | cut -c1-260
tail -1 output_$t/convergence.csv
touch WB14_DONE
