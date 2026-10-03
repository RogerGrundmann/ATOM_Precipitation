#!/bin/bash
# 2026-10-03: wb18 = wb14 + ATM_TEQ_WTG=1 (user). SCREENING arm: nm 220 from scratch, cli/atm_ea (-O2), 8 threads. One variable against wb14.
# wb14 at 200: 849 mm/a, r .458, sigma 3.74, land/ocean 194/1108, bands 2942/131/210/4.5; land 0-15 / 15-35 / 35-65 / 65-90 = 354 / 1 / 385 / 7;
# parcel theta_e - env theta_es min: Amazon -2.5 K, Congo -3.8, SE Asia -2.9, Sahara -14.4, Arabia -12.1, ocean<15 +2.3; env min 348-355 K over land, 345.3 ocean.
# ATM_TEQ_WTG was refuted on 10-02 (lo3: land P 1000 -> 3187; wb10: Arabia 35 mm/d) BEFORE ATM_RH_LAND_QCAP and the gate blend.
# PRE-REGISTERED: env theta_es min over tropical land falls to ~345-346 K; Amazon and Congo margins ~0 (+-1 K), SE Asia ~-2: forests convect
# weakly-to-moderately through the blend -- Amazon 1-5 mm/d, Congo 2-6; land 0-15 354 -> 800-1800; Sahara stays 0 (margin ~-8 K);
# Arabia margin ~-2 K: 0-1 mm/d (the RISK: > 2 mm/d = the old refutation again); land 15-35 < 150; ocean within 5 % of wb14; no land cell > 50 mm/d.
set -u; cd "$(dirname "$0")"; rm -f WB18_DONE
mkdir output_wb18 || { touch WB18_DONE; exit 1; }
[ -e config_wb18.xml ] && { touch WB18_DONE; exit 1; }
sed -e "s#output_wb7/#output_wb18/#" -e "s#<nm>600<#<nm>220<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>220<#" config_wb7.xml > config_wb18.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_ea $(md5sum < ../cli/atm_ea | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_RH_OCEAN=0.805 ATM_MC_ML_PARCEL=1 ATM_MC_T_ADD=0.15 ATM_MC_Q_ADD=4.5e-4 ATM_MC_ML_LCL=1 ATM_RH_LAND=2 \
    ATM_MC_GATE_BLEND=1 ATM_RH_LAND_QCAP=1 ATM_TEQ_WTG=1 ../cli/atm_ea config_wb18.xml > wb18.log 2>&1
t=wb18
echo "== $t  exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
grep -a "RH-LAND" $t.log | head -3
grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1
grep -a "P_conv mean" $t.log | tail -1
grep -a "\[MC-LO\] ocean<30\|\[MC-LO\] land<30" $t.log | tail -2 | cut -c1-330
grep -a "\[MC-RG\]" $t.log | tail -6 | cut -c1-260
tail -1 output_$t/convergence.csv
python3 landb.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk
touch WB18_DONE
