#!/bin/bash
# 2026-10-02: wb10a / wb10, 600 from scratch, cli/atm_rl2 (-O2), 2 x 8 threads. = wb9 (RH_OCEAN 0.82, ML_PARCEL 1 T .15 Q 4.5e-4, ML_LCL 1)
#   with ATM_RH_LAND=2 (wet tropical land starts like the tropical ocean);  wb10 adds ATM_TEQ_WTG=1 (no land-ocean T contrast aloft).
# wb9: 1029 mm/a, land/ocean 122/1389, r .436, sigma 4.39, bands 3440/296/220/4.6; Amazon/Congo P_conv 0 (rg1: theta_e - min -5.4 / -6.7 K).
# PRE-REGISTERED: wb10a -- Amazon/Congo parcels gain ~3.5 K but stay mostly below the bar (-2..-3 K): little forest rain. wb10 --
# forests near +1 / 0 K and convecting (P_conv > 0, land well above wb9's 122), deserts still dry, ocean roughly unchanged.
set -u; cd "$(dirname "$0")"; rm -f WB10_DONE
for t in wb10a wb10; do
  mkdir output_$t || { touch WB10_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB10_DONE; exit 1; }
  sed "s#output_wb7/#output_$t/#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_rl2 $(md5sum < ../cli/atm_rl2 | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_RH_OCEAN=0.82 ATM_MC_ML_PARCEL=1 ATM_MC_T_ADD=0.15 ATM_MC_Q_ADD=4.5e-4 ATM_MC_ML_LCL=1 ATM_RH_LAND=2"
env OMP_NUM_THREADS=8 $K                ../cli/atm_rl2 config_wb10a.xml > wb10a.log 2>&1 &
env OMP_NUM_THREADS=8 $K ATM_TEQ_WTG=1  ../cli/atm_rl2 config_wb10.xml  > wb10.log  2>&1 &
wait
for t in wb10a wb10; do
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
  grep -a "RH-LAND\|TEQ-WTG" $t.log | head -2
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1
  grep -a "P_conv mean" $t.log | tail -1
  grep -a "\[MC-RG\]" $t.log | tail -6 | cut -c1-260
  tail -1 output_$t/convergence.csv
done
touch WB10_DONE
