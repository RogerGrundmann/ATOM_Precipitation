#!/bin/bash
# 2026-10-04: wb24a / wb24b = wb23a / wb23b + ATM_MC_GP_AREA=1 (user), i.e. the wb20b stack with ATM_MC_GATE_BLEND=3 (undiluted parcel,
# gate threshold 0.01) and g_p = a_u*K_p*q_c_u (a_u = 0.03) instead of K_p*q_c_u. ATM_RH_OCEAN 0.805 (a) and 0.79 (b, the cliff test).
# ONE VARIABLE against wb23a / wb23b. SCREENING: nm 220 from scratch, cli/atm_g3 (-O2, = cli/atm at 759f496), 2 x 8 threads.
# wb23a at 200: 4335 mm/a, r .397, sigma 11.97, bands 13616/2703/210/4.5, P_conv ~4250; every convecting column 80-130 mm/d
# (Amazon 24374, SE Asia 21075, ocean<15 11113 mm/a). wb23b: 3434 mm/a, P_conv 3349 (b/a 0.79).
# wb20b (GATE_BLEND=1, no GP_AREA) at 200: 927 mm/a, r .472, sigma 4.12, land/ocean 459/1113, bands 3253/131/210/4.5, P_conv 807,
# 2.0 % of tropical land cells > 50 mm/d carrying 43 % of the rain, wettest land cell 105 mm/d.
# PRE-REGISTERED: the conversion rate falls 33x, but q_c_u rises toward its 10 g/kg ceiling and the unrained condensate leaves by
# detrainment (e_l up), so P_conv falls 6-15x, not 33x: a -- P_conv 280-700, global 380-800 (BELOW NASA), 0-15 1000-2400,
# ocean<15 / Amazon / SE Asia P_conv 700-2500 mm/a each (2-7 mm/d), wettest cell < 40 mm/d, no cell > 50 mm/d, sigma 1.5-3.5, r .40-.47,
# land 15-35 stays ~1, 35-65 / 65-90 equal to wb23a (210 / 4.5), deserts 0; b -- smooth, wb24b/wb24a 0.7-0.95. Exit 0, zero NaN.
# RISK: the retained condensate re-evaporates aloft (MC_q moistening) and stratiform rain picks it up -> P_rain rises, total falls less.
set -u; cd "$(dirname "$0")"; rm -f WB24_DONE
for t in wb24a wb24b; do
  mkdir output_$t || { touch WB24_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB24_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>220<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>220<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_g3 $(md5sum < ../cli/atm_g3 | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_MC_ML_PARCEL=1 ATM_MC_T_ADD=0.15 ATM_MC_Q_ADD=4.5e-4 ATM_MC_ML_LCL=1 ATM_RH_LAND=2 ATM_MC_GATE_BLEND=3 ATM_MC_GP_AREA=1 ATM_MC_T_ADD_LAND=2.6 ATM_RH_LAND_QCAP=1"
env OMP_NUM_THREADS=8 $K ATM_RH_OCEAN=0.805 ../cli/atm_g3 config_wb24a.xml > wb24a.log 2>&1 &
env OMP_NUM_THREADS=8 $K ATM_RH_OCEAN=0.79 ../cli/atm_g3 config_wb24b.xml > wb24b.log 2>&1 &
wait
for t in wb24a wb24b; do
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
  grep -a "RUN CONFIG" $t.log | grep -o "MC_GP_AREA=[^ ]*\|MC_GATE_BLEND=[^ ]*" | sort -u | tr '\n' ' '; echo
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1
  grep -a "P_conv mean\|P_rain mean" $t.log | tail -2
  grep -a "\[MC-GATE\]" $t.log | tail -2 | cut -c1-260
  grep -a "\[MC-LO\] ocean<30\|\[MC-LO\] land<30" $t.log | tail -2 | cut -c1-330
  grep -a "\[MC-RG\]" $t.log | tail -6 | cut -c1-260
  tail -1 output_$t/convergence.csv
  python3 landb.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk
done
touch WB24_DONE
