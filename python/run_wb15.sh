#!/bin/bash
# 2026-10-03: wb15 = wb14 + ATM_RH_SIGMA_SFC=1 (user). 600 from scratch, cli/atm_ss (-O2), 8 threads; starts when wb14 has finished.
# wb14 (interim, iteration ~220): land 0-15 / 15-35 / 35-65 / 65-90 = 354 / 1 / 385 / 7 mm/a (NASA 1653 / 643 / 643 / 295), land<30 CAPE
# 2.3 J/kg, P_conv 0.1 mm/a; parcel theta_e - env min: Amazon -2.5 K (RH 0.72), Congo -3.8, SE Asia -2.9, ocean<15 +2.3 (RH 0.78).
# PRE-REGISTERED: forest parcel RH 0.72 -> 0.77-0.79, theta_e + ~3 K: Amazon and SE Asia margins >= 0 and they convect, Congo still
# slightly negative (weak rain through the gate blend); land 0-15 354 -> 800-1500 mm/a; global land 195 -> 350-500; ocean unchanged
# (1107 +- 3 %); 15-35 land still ~0 (RH_LAND dries it); Horn of Africa still dry (cap); deserts still dry; r >= wb14's.
set -u; cd "$(dirname "$0")"; rm -f WB15_DONE
until [ -e WB14_DONE ]; do sleep 20; done
mkdir output_wb15 || { touch WB15_DONE; exit 1; }
[ -e config_wb15.xml ] && { touch WB15_DONE; exit 1; }
sed "s#output_wb7/#output_wb15/#" config_wb7.xml > config_wb15.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_ss $(md5sum < ../cli/atm_ss | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_RH_OCEAN=0.805 ATM_MC_ML_PARCEL=1 ATM_MC_T_ADD=0.15 ATM_MC_Q_ADD=4.5e-4 ATM_MC_ML_LCL=1 ATM_RH_LAND=2 \
    ATM_MC_GATE_BLEND=1 ATM_RH_LAND_QCAP=1 ATM_RH_SIGMA_SFC=1 ../cli/atm_ss config_wb15.xml > wb15.log 2>&1
t=wb15
echo "== $t  exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
grep -a "RH-LAND-QCAP" $t.log | head -1
grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1
grep -a "P_conv mean" $t.log | tail -1
grep -a "\[MC-GATE\]" $t.log | tail -2 | cut -c1-260
grep -a "\[MC-LO\] ocean<30\|\[MC-LO\] land<30" $t.log | tail -2 | cut -c1-330
grep -a "\[MC-RG\]" $t.log | tail -6 | cut -c1-260
tail -1 output_$t/convergence.csv
touch WB15_DONE
