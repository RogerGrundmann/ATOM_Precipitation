#!/bin/bash
# 2026-09-28: the ocean seam blow-up. osf_40 / osf_26 / ohs_sh1 (shipped metric, restart 1000 from oc_ctl, -O2) all
# grow max|u| 0.13 -> ~6.5 -> ~80 m/s at 55N, k = 0 (the phi seam) by iter 1200. oc_ctl (made 09-21) was clean.
# Suspect: HYD_BC_SECOND_ORDER, default 1 since 0571bc1 (09-23) -- the atmosphere's seam story.
# Two arms, cli/hyd_rp (-O2, c1aae23), 1000 -> 1200 from oc_ctl, config_osf_40 derivative (floor 0.4) with
# checkpoint 20 for the trajectory, 4 threads each; starts when run_rp0928.sh's three ocean runs have finished.
#   sm_bc1  HYD_BC_SECOND_ORDER=1 (default)  -- same-binary control; must reproduce the blow-up
#   sm_bc0  HYD_BC_SECOND_ORDER=0 (first-order, the pre-09-23 branch)
# PRE-REGISTERED: sm_bc1 reaches max|u| > 10 m/s at 55N 0E by 1200; sm_bc0 stays O(0.1) m/s -> the second-order BC
# reactivated a seam mode in the ocean and needs reverting/seam repair as in the atmosphere. If sm_bc0 blows up too,
# the suspect is refuted and the cause lies elsewhere in 09-21..09-28.
set -u; cd "$(dirname "$0")"; rm -f SM0928_DONE
until [ "$(grep -cE '^rp_(oc|od240|od2400) exit' run_rp0928.out)" = 3 ]; do sleep 60; done
for t in bc1 bc0; do mkdir output_sm_$t || { touch SM0928_DONE; exit 1; }
  sed -e "s#output_osf_40/#output_sm_$t/#" -e "s#<checkpoint>100</checkpoint>#<checkpoint>20</checkpoint>#" config_osf_40.xml > config_sm_$t.xml
  cp output_oc_ctl/hyd_restart_0Ma_1000.bin output_sm_$t/; cp output_twctl/0Ma_smooth_Transfer_Atm_600.vwtp output_sm_$t/; done
echo "start $(date +%H:%M)  hyd $(md5sum < ../cli/hyd_rp | cut -c1-8)"
run(){ ( env OMP_NUM_THREADS=${NT:-4} $2 ../cli/hyd_rp config_sm_$1.xml > sm_$1.log 2>&1
         echo "sm_$1 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' sm_$1.log)  $(date +%H:%M)" ) & }
run bc1 "HYD_BC_SECOND_ORDER=1"
run bc0 "HYD_BC_SECOND_ORDER=0"
wait
for t in bc1 bc0; do echo "== sm_$t $(grep -o 'BC_SECOND_ORDER=[^ ]*' sm_$t.log | head -1)"; grep 'max u-component' sm_$t.log | cut -c1-85; done
touch SM0928_DONE
