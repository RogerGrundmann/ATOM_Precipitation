#!/bin/bash
# 2026-10-02: wb11 = wb10a with ATM_RH_OCEAN 0.79 (was 0.82; user). 600 from scratch, cli/atm_rl2 (-O2, = 310455e source), 8 threads.
# wb10a: 1148 mm/a, land/ocean 486/1410, r .407, sigma 4.84, bands 3883/301/238/4.6; ocean<15 theta_e - bar +3.3 K (91 % buoyant),
# Amazon -1.6 K (2.1 mm/d), Congo -3.0 K (3.9 mm/d).
# PRE-REGISTERED: ocean BL q -0.75 g/kg, theta_e -~1.9 K -> ocean<15 margin ~+1.4 K, ocean rain falls 30-50 %, 0-15 band toward NASA;
# but RH_LAND=2 ties wet land to the same value, so Amazon/Congo margins fall ~2 K and their convective rain falls toward 0.
set -u; cd "$(dirname "$0")"; rm -f WB11_DONE
mkdir output_wb11 || { touch WB11_DONE; exit 1; }
[ -e config_wb11.xml ] && { touch WB11_DONE; exit 1; }
sed "s#output_wb7/#output_wb11/#" config_wb7.xml > config_wb11.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_rl2 $(md5sum < ../cli/atm_rl2 | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_RH_OCEAN=0.79 ATM_MC_ML_PARCEL=1 ATM_MC_T_ADD=0.15 ATM_MC_Q_ADD=4.5e-4 ATM_MC_ML_LCL=1 ATM_RH_LAND=2 \
    ../cli/atm_rl2 config_wb11.xml > wb11.log 2>&1
t=wb11
echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
grep -a "RH-LAND" $t.log | head -1
grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1
grep -a "P_conv mean" $t.log | tail -1
grep -a "\[MC-RG\]" $t.log | tail -6 | cut -c1-260
tail -1 output_$t/convergence.csv
touch WB11_DONE
