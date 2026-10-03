#!/bin/bash
# 2026-10-03: wb20a / wb20b = wb14 + ATM_MC_T_ADD_LAND 2.3 / 2.6 K (user). SCREENING: nm 220 from scratch, cli/atm_tl (-O2), 2 x 8 threads.
# wb20a (2.0 K) at 200: 872 mm/a, r .464, sigma 3.83, land/ocean 268/1111, land 0-15 732, Amazon 1.56 mm/d, Congo 1.98, max cell 96 mm/d (7S 114E).
# wb20b (3.0 K):        995 mm/a, r .473, sigma 4.52, land/ocean 695/1114, land 0-15 2930, Amazon 18.15, Congo 2.16, max cell 113 mm/d (4S 144E).
# NASA: 978, land/ocean 782/1056, land 0-15 1653, Amazon 6.80, Congo 4.75. The Amazon margin moved 1.0 K per K (-0.5 at 2.0, +0.5 at 3.0).
# PRE-REGISTERED (the response is steep, convex): 2.3 K -- Amazon 3-7 mm/d, land 0-15 1000-1600, land 330-450, global 885-925;
# 2.6 K -- Amazon 7-13 mm/d, land 0-15 1600-2300, land 450-580, global 925-965. Both: ocean 1108-1116, deserts 0, Congo 2.0-2.2,
# land 15-35 = 1, r .465-.475, sigma between 3.9 and 4.4, single land cells still > 80 mm/d.
set -u; cd "$(dirname "$0")"; rm -f WB20_DONE
for t in wb20a wb20b; do
  mkdir output_$t || { touch WB20_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB20_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>220<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>220<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_tl $(md5sum < ../cli/atm_tl | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_RH_OCEAN=0.805 ATM_MC_ML_PARCEL=1 ATM_MC_T_ADD=0.15 ATM_MC_Q_ADD=4.5e-4 ATM_MC_ML_LCL=1 ATM_RH_LAND=2 ATM_MC_GATE_BLEND=1 ATM_RH_LAND_QCAP=1"
env OMP_NUM_THREADS=8 $K ATM_MC_T_ADD_LAND=2.3 ../cli/atm_tl config_wb20a.xml > wb20a.log 2>&1 &
env OMP_NUM_THREADS=8 $K ATM_MC_T_ADD_LAND=2.6 ../cli/atm_tl config_wb20b.xml > wb20b.log 2>&1 &
wait
for t in wb20a wb20b; do
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1
  grep -a "P_conv mean" $t.log | tail -1
  grep -a "\[MC-LO\] ocean<30\|\[MC-LO\] land<30" $t.log | tail -2 | cut -c1-330
  grep -a "\[MC-RG\]" $t.log | tail -6 | cut -c1-260
  tail -1 output_$t/convergence.csv
  python3 landb.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk
done
touch WB20_DONE
