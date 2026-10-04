#!/bin/bash
# 2026-10-04: wb25a / wb25b = the wb24 stack (wb20b + ATM_MC_GATE_BLEND=3 + ATM_MC_GP_AREA=1) at ATM_RH_OCEAN=0.80 with
# ATM_MC_T_ADD_LAND 2.6 (a) / 2.3 (b) (user). SCREENING: nm 220 from scratch, cli/atm_g3 (-O2, = cli/atm at 759f496), 2 x 8 threads.
# wb24a (0.805, 2.6): 1128 mm/a, r .572, sigma 3.38, bands 3716/400/210/4.5, land/ocean 577/1346, land 0-15 2318, Amazon 11.48, Congo 4.80 mm/d.
# wb24b (0.79, 2.6):   615 mm/a, r .554, sigma 2.03, bands 2010/146/194/4.5, land/ocean 364/714,  land 0-15 1249, Amazon 6.53,  Congo 2.13 mm/d.
# NASA: 978, bands 1487/761/981/364, land/ocean 782/1056, land 0-15 1653, Amazon 6.80, Congo 4.75.
# PRE-REGISTERED (0.80 is 2/3 of the way from 0.79 to 0.805; the response is convex, so below the linear value):
# a -- global 860-960 (linear 957), ocean 1000-1140, land 440-510, land 0-15 1700-2000, Amazon 8.3-9.9, Congo 3.2-4.0 mm/d, 0-15 2900-3200,
# 15-35 250-320, r .56-.575, sigma 2.7-3.0, P/E 1.5-1.7; b -- ocean equal to a within 1 %, land 380-450, land 0-15 1100-1500,
# Amazon 4.5-6.5, Congo 2.5-3.6, global 20-45 mm/a below a, r within .005 of a. Both: no cell > 50 mm/d, deserts 0, land 15-35 < 10,
# 35-65 200-210, 65-90 4.5, exit 0, zero NaN. wb20 series (blend 1) had Amazon 4.2 / 8.8 at 2.3 / 2.6 K -- steeper than expected here.
set -u; cd "$(dirname "$0")"; rm -f WB25_DONE
for t in wb25a wb25b; do
  mkdir output_$t || { touch WB25_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB25_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>220<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>220<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_g3 $(md5sum < ../cli/atm_g3 | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_MC_ML_PARCEL=1 ATM_MC_T_ADD=0.15 ATM_MC_Q_ADD=4.5e-4 ATM_MC_ML_LCL=1 ATM_RH_LAND=2 ATM_MC_GATE_BLEND=3 ATM_MC_GP_AREA=1 ATM_RH_OCEAN=0.80 ATM_RH_LAND_QCAP=1"
env OMP_NUM_THREADS=8 $K ATM_MC_T_ADD_LAND=2.6 ../cli/atm_g3 config_wb25a.xml > wb25a.log 2>&1 &
env OMP_NUM_THREADS=8 $K ATM_MC_T_ADD_LAND=2.3 ../cli/atm_g3 config_wb25b.xml > wb25b.log 2>&1 &
wait
for t in wb25a wb25b; do
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
  grep -a "RUN CONFIG" $t.log | grep -o "MC_GP_AREA=[^ ]*\|MC_GATE_BLEND=[^ ]*\|RH_OCEAN=[^ ]*\|MC_T_ADD_LAND=[^ ]*" | sort -u | tr '\n' ' '; echo
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1
  grep -a "P_conv mean\|P_rain mean" $t.log | tail -2
  grep -a "\[MC-GATE\]" $t.log | tail -2 | cut -c1-260
  grep -a "\[MC-LO\] ocean<30\|\[MC-LO\] land<30" $t.log | tail -2 | cut -c1-330
  grep -a "\[MC-RG\]" $t.log | tail -6 | cut -c1-260
  tail -1 output_$t/convergence.csv
  python3 landb.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk
done
touch WB25_DONE
