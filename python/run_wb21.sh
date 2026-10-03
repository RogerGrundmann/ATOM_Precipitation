#!/bin/bash
# 2026-10-03: wb21a / wb21b = wb14 + ATM_MC_T_ADD_LAND 2.6 / 3.0 K with the coastal taper ATM_MC_T_ADD_LAND_L=300 km (user). SCREENING: nm 220 from
# scratch, cli/atm_ct (-O2), 2 x 8 threads; starts when wb20 has finished. One variable (the taper) against wb20b (2.6 K) and wb21b (3.0 K).
# wb21b (3.0 K, no taper) at 200: 995 mm/a, r .473, sigma 4.52, land/ocean 695/1114, land 0-15 2930, Amazon 18.15 mm/d, Congo 2.16; Maritime
# Continent (10S-10N 95-155E) land mean 23.9 mm/d (NASA 8.1), 40 of 173 land cells > 50 mm/d, coastal 39 m cells 96-113 mm/d.
# Taper factor 1 - exp(-d/300 km): 0.31 at one cell from the coast (111 km), 0.52 at two, 0.81 at 500 km, 0.96 at 1000 km.
# PRE-REGISTERED: coastal 39 m cells fall to the neighbouring ocean's level (30-45 mm/d); Maritime-Continent land cells > 50 mm/d: 40 -> < 10 (b);
# its land mean 23.9 -> 5-12 mm/d; the INLAND hot lowland cells (Amazon 2-4S 75-76W, South Sudan 6-7N 32-33E; d >= 500 km) keep ~100 mm/d,
# so the global land maximum stays ~100; Amazon mean falls by less than a third; land (b) 695 -> 450-600; ocean and deserts unchanged; r >= wb21b's.
set -u; cd "$(dirname "$0")"; rm -f WB21_DONE
until [ -e WB20_DONE ]; do sleep 20; done
for t in wb21a wb21b; do
  mkdir output_$t || { touch WB21_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB21_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>220<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>220<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_ct $(md5sum < ../cli/atm_ct | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_RH_OCEAN=0.805 ATM_MC_ML_PARCEL=1 ATM_MC_T_ADD=0.15 ATM_MC_Q_ADD=4.5e-4 ATM_MC_ML_LCL=1 ATM_RH_LAND=2 ATM_MC_GATE_BLEND=1 ATM_RH_LAND_QCAP=1"
env OMP_NUM_THREADS=8 $K ATM_MC_T_ADD_LAND=2.6 ATM_MC_T_ADD_LAND_L=300 ../cli/atm_ct config_wb21a.xml > wb21a.log 2>&1 &
env OMP_NUM_THREADS=8 $K ATM_MC_T_ADD_LAND=3.0 ATM_MC_T_ADD_LAND_L=300 ../cli/atm_ct config_wb21b.xml > wb21b.log 2>&1 &
wait
for t in wb21a wb21b; do
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1
  grep -a "P_conv mean" $t.log | tail -1
  grep -a "\[MC-LO\] ocean<30\|\[MC-LO\] land<30" $t.log | tail -2 | cut -c1-330
  grep -a "\[MC-RG\]" $t.log | tail -6 | cut -c1-260
  tail -1 output_$t/convergence.csv
  python3 landb.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk
done
touch WB21_DONE
