#!/bin/bash
# 2026-09-28: OCN-METRIC -- should the ocean adopt the correct horizontal metric? The 09-21 ladder (oc_*) is now
# stale: it ran FIRST-order BCs, floor 0.4 and no seam fix. Re-run the decisive pair from scratch on TODAY's
# defaults (second-order BCs + HYD_SEAM_PERIODIC=1 + floor 0.26 + all 09-28 flips), 1000 iterations, -O2.
#   om_ship  shipped metric (current default)
#   om_met   HYD_METRIC_RADIUS=6370 HYD_RUN_NEUMANN=1 HYD_A_H_BIHARM=3.0e18  (oc_vis's trio)
# Binary: -O2 build of 9438b0c in a throw-away worktree (cli/hyd_om). 2 x 6 threads. Starts only after the
# two pending -O0 byte checks (run_vsqf0928.sh, run_vosf0928.sh) have PASSED, so it runs on verified defaults.
# SCORED WITH ocprofile.py on one fixed full-depth column set + convergence_hyd.csv.
# PRE-REGISTERED (from the 09-21 ladder): om_met has rms radial u ~1e-4 of om_ship's (2.4e-2 m/s there), grid
# noise at or below om_ship, NO outlier (max surface speed ~ om_ship's), a bottom-intensified rms profile
# (+50 % from -65 m to -200 m, the missing deep outlet), KE drift ~1.6 % vs ~1.1 %, converged 0 in both.
# DECISION RULE, written before the run: recommend the flip if the radial-velocity and noise results hold and
# the profile/KE-drift costs are no worse than on 09-21; recommend NOT if om_met is unstable or noisier.
set -u; cd "$(dirname "$0")"; rm -f OM0928_DONE
until [ -e VSQF0928_DONE ] && [ -e VOSF0928_DONE ]; do sleep 60; done
for f in run_vsqf0928.out run_vosf0928.out; do
  if [ "$(grep 'DIFFERS:' $f | grep -vc 'RUN_CONFIG.txt')" != 0 ] || ! grep -q 'C  CONTROL.*PASS' $f; then
     echo "byte check $f did not pass -- NOT started"; touch OM0928_DONE; exit 1; fi; done
W=$(mktemp -d /tmp/atom_o2_XXXXXX)
(cd .. && git worktree add --detach $W 9438b0c >/dev/null && cd $W && make -j6 hyd > build.log 2>&1) \
  && cp -n $W/cli/hyd ../cli/hyd_om; (cd .. && git worktree remove --force $W; git worktree prune)
[ -e ../cli/hyd_om ] || { echo "O2 build failed"; touch OM0928_DONE; exit 1; }
for t in ship met; do mkdir output_om_$t || { touch OM0928_DONE; exit 1; }
  sed "s#output_oc_ctl/#output_om_$t/#" config_oc_ctl.xml > config_om_$t.xml
  cp output_twctl/0Ma_smooth_Transfer_Atm_600.vwtp output_om_$t/; done
echo "start $(date +%H:%M)  hyd_om $(md5sum < ../cli/hyd_om | cut -c1-8)"
run(){ ( env OMP_NUM_THREADS=${NT:-6} $2 ../cli/hyd_om config_om_$1.xml > om_$1.log 2>&1
         echo "om_$1 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' om_$1.log)  $(date +%H:%M)" ) & }
run ship ""
run met  "HYD_METRIC_RADIUS=6370 HYD_RUN_NEUMANN=1 HYD_A_H_BIHARM=3.0e18"
wait
python3 ocprofile.py ship:output_om_ship/hyd_restart_0Ma_1000.bin met:output_om_met/hyd_restart_0Ma_1000.bin > om0928_profile.txt 2>&1
cat om0928_profile.txt
for t in ship met; do echo "== om_$t conv: $(tail -1 output_om_$t/convergence_hyd.csv)"; grep 'max u-component' om_$t.log | tail -2 | cut -c1-85; done
touch OM0928_DONE
