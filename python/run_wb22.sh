#!/bin/bash
# 2026-10-03: wb22a / wb22b = wb20b (wb14 + ATM_MC_T_ADD_LAND 2.6 K, no taper) + ATM_MC_MB_SAT_LAND 0.06 / 0.04 kg/(m2 s) (user: cure the inland
# hotspots at the closest tropical fit). SCREENING: nm 220 from scratch, cli/atm_mb (-O2), 2 x 8 threads; starts when wb21 has finished.
# wb20b at 200: 927 mm/a, r .472, sigma 4.12, land/ocean 459/1113, land 0-15 1714 (NASA 1653), Amazon 8.78 mm/d (max cell 74), Congo 2.01,
# global max 104.6 mm/d (2S 139E, land); tropical land: 14.2 % of cells convect, 2.0 % rain > 50 mm/d and carry 43 % of the convective rain.
# Ocean: mean M_b 0.037 -> a near-uniform 31-32 mm/d. The gate blend weights the parcel by smoothstep(|M| 0.01..0.1): 0.26 at 0.04, 0.60 at 0.06.
# PRE-REGISTERED: no land cell above 60 mm/d (a) / 40 mm/d (b) (wb20b: 105); tropical-land cells > 50 mm/d: 2.0 % -> < 0.3 % (a), 0 (b);
# the convecting fraction stays ~14 % (the saturation changes how hard a column rains, not whether); so the land mean FALLS: land 0-15
# 1714 -> 1000-1400 (a) / 600-1000 (b), Amazon 8.8 -> 5-7 / 3-5 mm/d; ocean unchanged (1113 +- 0.5 %), deserts 0; sigma 4.12 -> 3.8-4.0.
set -u; cd "$(dirname "$0")"; rm -f WB22_DONE
until [ -e WB21_DONE ]; do sleep 20; done
for t in wb22a wb22b; do
  mkdir output_$t || { touch WB22_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB22_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>220<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>220<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_mb $(md5sum < ../cli/atm_mb | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_RH_OCEAN=0.805 ATM_MC_ML_PARCEL=1 ATM_MC_T_ADD=0.15 ATM_MC_Q_ADD=4.5e-4 ATM_MC_ML_LCL=1 ATM_RH_LAND=2 ATM_MC_GATE_BLEND=1 ATM_RH_LAND_QCAP=1"
env OMP_NUM_THREADS=8 $K ATM_MC_T_ADD_LAND=2.6 ATM_MC_MB_SAT_LAND=0.06 ../cli/atm_mb config_wb22a.xml > wb22a.log 2>&1 &
env OMP_NUM_THREADS=8 $K ATM_MC_T_ADD_LAND=2.6 ATM_MC_MB_SAT_LAND=0.04 ../cli/atm_mb config_wb22b.xml > wb22b.log 2>&1 &
wait
for t in wb22a wb22b; do
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1
  grep -a "P_conv mean" $t.log | tail -1
  grep -a "\[MC-LO\] ocean<30\|\[MC-LO\] land<30" $t.log | tail -2 | cut -c1-330
  grep -a "\[MC-RG\]" $t.log | tail -6 | cut -c1-260
  tail -1 output_$t/convergence.csv
  python3 landb.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk
done
touch WB22_DONE
