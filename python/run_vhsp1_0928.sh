#!/bin/bash
# 2026-09-28: HYD_SEAM_PERIODIC flipped ON (d5d0272) -- byte check both directions, then re-run every
# shipped-metric ocean arm the seam blow-up spoiled, on the new default.
# 1) -O0, 1 thread, restart 1000 -> 1020 from oc_ctl: old = cli/hsp_new_hyd (ead0be3, knob default 0), new = d5d0272
#      A  new clean == old HYD_SEAM_PERIODIC=1;  B  new HYD_SEAM_PERIODIC=0 == old clean;  C  new clean != old clean
#    (RUN_CONFIG may differ only by the SEAM_PERIODIC token)
# 2) -O2 build of d5d0272 in a throw-away worktree -> cli/hyd_ro; 5 arms x 4 threads, shipped metric, from oc_ctl:
#      ro_osf40  1000->1200  floor 0.4 (control)        ro_osf26  1000->1200  HYD_METRIC_SIN_FLOOR=0.26
#      ro_shctl  1000->1600  control                    ro_ohssh1 1000->1600  HYD_HYDRO_SPLIT=1.0
#      ro_bnsh1  1000->1600  HYD_BUOY_CONSISTENT=1.0
#    Starts only after run_sq20928.sh is done (cores).
# PRE-REGISTERED: no arm blows up at the seam (max|u| stays O(0.1) m/s unless the knob under test itself runs away).
#   osf: equatorward of 66 deg identical to thread noise, 75-90 moves a few %. ohs_sh1: stable (split keeps the
#   buoyancy out of rhs_u). bn_sh1: the 07-14 buoyancy blow-up EITHER reappears (then it was real) OR not (then
#   the recorded 77.7 m/s was the seam) -- no expectation registered.
set -u; cd "$(dirname "$0")"; rm -f VHSP1_0928_DONE
. ./verify_lib.sh
[ -e ../cli/hsp1_new_hyd ] || ./build_o0.sh d5d0272 hsp1_new || { touch VHSP1_0928_DONE; exit 1; }
for t in a b c d; do mkdir output_vhsp1_$t || { touch VHSP1_0928_DONE; exit 1; }
  sed -e "s#output_sm_bc1/#output_vhsp1_$t/#" -e "s#<nm>200</nm>#<nm>20</nm>#" \
      -e "s#<checkpoint_save_iter>1200</checkpoint_save_iter>#<checkpoint_save_iter>1020</checkpoint_save_iter>#" config_sm_bc1.xml > config_vhsp1_$t.xml
  cp output_oc_ctl/hyd_restart_0Ma_1000.bin output_vhsp1_$t/; cp output_twctl/0Ma_smooth_Transfer_Atm_600.vwtp output_vhsp1_$t/; done
arm vhsp1_a ../cli/hsp1_new_hyd config_vhsp1_a.xml
arm vhsp1_b ../cli/hsp_new_hyd  config_vhsp1_b.xml HYD_SEAM_PERIODIC=1
arm vhsp1_c ../cli/hsp1_new_hyd config_vhsp1_c.xml HYD_SEAM_PERIODIC=0
arm vhsp1_d ../cli/hsp_new_hyd  config_vhsp1_d.xml
wait_arms
for t in a b c d; do [ -f output_vhsp1_$t/RUN_CONFIG.txt ] && sed -i 's#output_vhsp1_[a-d]/#OUT/#' output_vhsp1_$t/RUN_CONFIG.txt; done
cmp_dirs vhsp1_a vhsp1_b "A  FLIP == old+knob"
cmp_dirs vhsp1_c vhsp1_d "B  new+knob=0 == old clean"
want_differ vhsp1_a vhsp1_d "C  CONTROL new clean vs old clean -- MUST differ"
BAD=$(grep 'DIFFERS:' run_vhsp1_0928.out | grep -vc 'RUN_CONFIG.txt')
if [ "$BAD" != 0 ] || ! grep -q "C  CONTROL.*PASS" run_vhsp1_0928.out; then
    echo "byte check did not pass -- arms NOT started"; touch VHSP1_0928_DONE; exit 1; fi
W=$(mktemp -d /tmp/atom_o2_XXXXXX)
(cd .. && git worktree add --detach $W d5d0272 >/dev/null && cd $W && make -j6 hyd > build.log 2>&1) \
  && cp -n $W/cli/hyd ../cli/hyd_ro; (cd .. && git worktree remove --force $W; git worktree prune)
[ -e ../cli/hyd_ro ] || { echo "O2 build failed"; touch VHSP1_0928_DONE; exit 1; }
until [ -e SQ20928_DONE ]; do sleep 60; done
for t in osf40 osf26 shctl ohssh1 bnsh1; do mkdir output_ro_$t || { touch VHSP1_0928_DONE; exit 1; }
  case $t in osf*) src=config_osf_40.xml; from=output_osf_40/ ;; *) src=config_ohs_sh1.xml; from=output_ohs_sh1/ ;; esac
  sed "s#$from#output_ro_$t/#" $src > config_ro_$t.xml
  cp output_oc_ctl/hyd_restart_0Ma_1000.bin output_ro_$t/; cp output_twctl/0Ma_smooth_Transfer_Atm_600.vwtp output_ro_$t/; done
echo "arms start $(date +%H:%M)  hyd_ro $(md5sum < ../cli/hyd_ro | cut -c1-8)"
run(){ ( env OMP_NUM_THREADS=${NT:-4} $2 ../cli/hyd_ro config_ro_$1.xml > ro_$1.log 2>&1
         echo "ro_$1 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' ro_$1.log)  $(date +%H:%M)" ) & }
run osf40  ""
run osf26  "HYD_METRIC_SIN_FLOOR=0.26"
run shctl  ""
run ohssh1 "HYD_HYDRO_SPLIT=1.0"
run bnsh1  "HYD_BUOY_CONSISTENT=1.0"
wait
for t in osf40 osf26 shctl ohssh1 bnsh1; do echo "== ro_$t $(grep -o 'SEAM_PERIODIC=[^ ]*' ro_$t.log | head -1)"; grep 'max u-component' ro_$t.log | tail -3 | cut -c1-85; done
touch VHSP1_0928_DONE
