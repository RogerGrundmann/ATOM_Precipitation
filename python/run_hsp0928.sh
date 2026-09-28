#!/bin/bash
# 2026-09-28: HYD_SEAM_PERIODIC (branch hyd-seam-periodic, ead0be3), the fix for the ocean phi-seam blow-up
# measured by run_sm0928.sh (HYD_BC_SECOND_ORDER=1: max|u| 0.128 -> 79 m/s at 55N 0E; =0: 0.128 flat).
# 1) -O0 byte check, 1 thread, restart 1000 -> 1020 from oc_ctl (shipped metric):
#      old = cli/conv_new_hyd (c1aae23 -O0), new = cli/hsp_new_hyd (hyd-seam-periodic -O0)
#      A  new clean == old clean (RUN_CONFIG may differ only by the new SEAM_PERIODIC=0* token)
#      B  CONTROL new SEAM_PERIODIC=1 != new clean -- the knob fires within 20 iterations
# 2) -O2 build of the branch (cli/hyd_hsp2), sp_ship: shipped metric, 1000 -> 1200, checkpoint 20, 4 threads,
#    HYD_SEAM_PERIODIC=1 (HYD_BC_SECOND_ORDER default 1) -- against sm_bc1 / sm_bc0.
# PRE-REGISTERED: sp_ship stays at max|u| ~0.128 m/s (as sm_bc0) and its fields match sm_bc0 far away from the
# seam more closely than sm_bc1; exit 0, zero NaN. If it still grows, the PressureSolverHyd p_dyn seam
# extrapolation (:374) is the next site.
set -u; cd "$(dirname "$0")"; rm -f HSP0928_DONE
. ./verify_lib.sh
W=/tmp/claude-1000/-home-roger-SynologyDrive-Cloudstation-Notebook-ATOM-Precipitation/c7b1af68-f74f-4048-957c-3452bc3a136a/scratchpad/hsp_wt
[ -e ../cli/hsp_new_hyd ] || ./build_o0.sh hyd-seam-periodic hsp_new || { touch HSP0928_DONE; exit 1; }
for t in old new on; do mkdir output_vhsp_$t || { touch HSP0928_DONE; exit 1; }
  sed -e "s#output_sm_bc1/#output_vhsp_$t/#" -e "s#<nm>200</nm>#<nm>20</nm>#" \
      -e "s#<checkpoint_save_iter>1200</checkpoint_save_iter>#<checkpoint_save_iter>1020</checkpoint_save_iter>#" config_sm_bc1.xml > config_vhsp_$t.xml
  cp output_oc_ctl/hyd_restart_0Ma_1000.bin output_vhsp_$t/; cp output_twctl/0Ma_smooth_Transfer_Atm_600.vwtp output_vhsp_$t/; done
arm vhsp_old ../cli/conv_new_hyd config_vhsp_old.xml
arm vhsp_new ../cli/hsp_new_hyd  config_vhsp_new.xml
arm vhsp_on  ../cli/hsp_new_hyd  config_vhsp_on.xml HYD_SEAM_PERIODIC=1
wait_arms
for t in old new on; do [ -f output_vhsp_$t/RUN_CONFIG.txt ] && sed -i 's#output_vhsp_[a-z]*/#OUT/#' output_vhsp_$t/RUN_CONFIG.txt; done
cmp_dirs vhsp_new vhsp_old "A  OFF BRANCH at -O0 (branch clean vs HEAD)"
want_differ vhsp_on vhsp_new "B  CONTROL SEAM_PERIODIC=1 -- MUST differ"
BAD=$(sed -n '1,/A  OFF BRANCH/p' run_hsp0928.out | grep 'DIFFERS:' | grep -vc 'RUN_CONFIG.txt')
if [ "$BAD" != 0 ] || ! grep -q "B  CONTROL.*PASS" run_hsp0928.out; then
    echo "byte check did not pass -- science arm NOT started"; touch HSP0928_DONE; exit 1; fi
(cd $W && make -j4 hyd > build_o2.log 2>&1) && cp -n $W/cli/hyd ../cli/hyd_hsp2 || { echo "O2 build failed"; touch HSP0928_DONE; exit 1; }
readelf --debug-dump=info ../cli/hyd_hsp2 | grep -m1 -o ' -O[0-3]'
mkdir output_sp_ship || { touch HSP0928_DONE; exit 1; }
sed "s#output_sm_bc1/#output_sp_ship/#" config_sm_bc1.xml > config_sp_ship.xml
cp output_oc_ctl/hyd_restart_0Ma_1000.bin output_sp_ship/; cp output_twctl/0Ma_smooth_Transfer_Atm_600.vwtp output_sp_ship/
echo "sp_ship start $(date +%H:%M)"
env OMP_NUM_THREADS=${NT:-4} HYD_SEAM_PERIODIC=1 ../cli/hyd_hsp2 config_sp_ship.xml > sp_ship.log 2>&1
echo "sp_ship exit $?  NaN $(grep -c 'NaN/Inf DETECTED' sp_ship.log)  $(date +%H:%M)"
echo "== sp_ship $(grep -o 'SEAM_PERIODIC=[^ ]*' sp_ship.log | head -1)"; grep 'max u-component' sp_ship.log | cut -c1-85
touch HSP0928_DONE
