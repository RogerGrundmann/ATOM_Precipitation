#!/bin/bash
# 2026-09-29: OCN-METRIC with the GRID-SCALED biharmonic. om_met (metric + floor 0.26 + B = 3e18) NaN'd at 89N at
# iteration 27: the explicit nabla^4 limit at floor 0.26 is ~9e17 (zonal dx 28.9 km). omy_b0 / omy_b1e18 (100 from
# scratch, B = 0 / 1e18 at floor 0.26) both ran clean, which confirms the diagnosis.
#   om_scl  HYD_METRIC_RADIUS=6370 HYD_RUN_NEUMANN=1 HYD_A_H_BIHARM=3.0e18 HYD_A_H_BIHARM_SCALED=1, floor 0.26 (default)
# compared with om_ship (shipped metric) and om_met40 (metric, floor 0.4, unscaled B = 3e18), both 1000 from scratch
# on 09-28, -O2, 6 threads. Same thread count here. Starts only if run_vbsc.sh PASSED.
# PRE-REGISTERED: stable through 1000; rms radial u ~1e-4 of om_ship; grid noise and max surface speed at or below
# om_met40's (the scaling weakens the damping poleward of ~50 deg, so polar noise is the thing to watch);
# bottom-intensified profile as in om_met40 (the scaling cannot touch it).
set -u; cd "$(dirname "$0")"; rm -f OMSC_DONE
until [ -e VBSC_DONE ]; do sleep 30; done
if [ "$(grep 'DIFFERS:' run_vbsc.out | grep -vc 'RUN_CONFIG.txt')" != 0 ] || ! grep -q 'C  CONTROL.*PASS' run_vbsc.out \
   || grep -q 'exit [^0]' run_vbsc.out; then echo "byte check did not pass -- NOT started"; touch OMSC_DONE; exit 1; fi
[ -e ../cli/hyd_omsc ] && { echo "cli/hyd_omsc exists"; touch OMSC_DONE; exit 1; }
W=$(mktemp -d /tmp/atom_o2_XXXXXX)
(cd .. && git worktree add --detach $W HEAD >/dev/null && git diff HEAD -- hydrosphere > $W/.wt.diff \
   && cd $W && { [ -s .wt.diff ] && git apply .wt.diff; true; } && make -j6 hyd > build.log 2>&1) \
  && cp -n $W/cli/hyd ../cli/hyd_omsc; (cd .. && git worktree remove --force $W; git worktree prune)
[ -e ../cli/hyd_omsc ] || { echo "O2 build failed"; touch OMSC_DONE; exit 1; }
mkdir output_om_scl || { touch OMSC_DONE; exit 1; }
sed "s#output_om_ship/#output_om_scl/#" config_om_ship.xml > config_om_scl.xml
cp output_twctl/0Ma_smooth_Transfer_Atm_600.vwtp output_om_scl/
echo "start $(date +%H:%M)  hyd_omsc $(md5sum < ../cli/hyd_omsc | cut -c1-8)"
env OMP_NUM_THREADS=6 HYD_METRIC_RADIUS=6370 HYD_RUN_NEUMANN=1 HYD_A_H_BIHARM=3.0e18 HYD_A_H_BIHARM_SCALED=1 \
    ../cli/hyd_omsc config_om_scl.xml > om_scl.log 2>&1
echo "om_scl exit $?  NaN $(grep -c 'NaN/Inf DETECTED' om_scl.log)  $(date +%H:%M)"
python3 ocprofile.py ship:output_om_ship/hyd_restart_0Ma_1000.bin met40:output_om_met40/hyd_restart_0Ma_1000.bin \
    scl:output_om_scl/hyd_restart_0Ma_1000.bin > omsc_profile.txt 2>&1
cat omsc_profile.txt
for t in ship met40 scl; do echo "== om_$t conv: $(tail -1 output_om_$t/convergence_hyd.csv)"; grep 'max u-component' om_$t.log | tail -1 | cut -c1-85; done
touch OMSC_DONE
