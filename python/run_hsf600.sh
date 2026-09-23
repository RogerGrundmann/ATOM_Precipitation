#!/bin/bash
# ATM_HYDRO_SPLIT=1.0 -- the 600-iteration FROM-SCRATCH pair a default flip owes.
# Binary cli/atm_hsf (md5 514e9c17 = cli/atm_bc2): physical-height split, AND second-order
# Neumann BCs as the new default (2026-09-23). ATM_RAD_TOPO=0 on both arms, because the working
# tree it was built from also carries the uncommitted RAD_TOPO flip.
#   hsf_ctl  new default, split OFF  -- ALSO the 600-iteration arm of the BC_SECOND_ORDER re-flip:
#            compare with output_el2 (same config, first-order BCs) for the k=1 seam mode.
#   hsf_on   new default + ATM_HYDRO_SPLIT=1.0
# WAITS FOR: the BC byte check to finish AND PASS on all six lines (banner-only RUN_CONFIG diffs
# accepted), and rt600 to release its 20 threads. Then 11 threads each, concurrent.
# PRE-REGISTERED (written before the run):
#   1. p_dyn (non-hydrostatic) PLATEAUS below the 3.0 ceiling; zero "CEILING BINDING" lines. The
#      restart grew 0.014 -> 0.116 nd over 80 iterations, decelerating. Clipping = FAIL.
#   2. max|u| stays within noise of hsf_ctl at every checkpoint (no radial runaway; the
#      =BUOY_CONSISTENT branch grew +2.8e-3 m/s per iteration from ~iteration 60).
#   3. max|v| and its LOCATION: hsf_ctl is expected to show the k=1 seam mode (2026-09-21:
#      2.28 -> 26 m/s, 11N 1E, 236 m). hsf_on must not make it worse or add a second locus.
#   4. v-budget 20-70 deg above 3 km: p05 residual ~0.06 once buoyancy_ramp reaches 1 (iter 300).
#   5. Precip / r / sigma / bands null to ~1 %; mean KE and `converged` comparable.
#   6. exit 0, zero NaN through 155 / 357 / 483.
set -u; cd "$(dirname "$0")"; rm -f HSF600_DONE
until [ -f BC2_VERIFY_DONE ] && [ -f RT600_DONE ]; do sleep 60; done
n_pass=$(grep -c "PASS" run_verify_bc2.out)
n_bad=$(grep "DIFFERS" run_verify_bc2.out | grep -vc "RUN_CONFIG.txt")
if [ "$n_pass" -lt 2 ] || [ "$n_bad" -ne 0 ]; then
    echo "BC byte check not clean (PASS=$n_pass, non-banner diffs=$n_bad) -- pair NOT started"
    touch HSF600_DONE; exit 1; fi
mkdir -p output_hsf_ctl output_hsf_on
echo "start $(date +%H:%M)"
( env OMP_NUM_THREADS=11 ATM_RAD_TOPO=0 ATM_UBUD_BALANCE=1 \
      ../cli/atm_hsf config_hsf_ctl.xml > hsf_ctl.log 2>&1; echo "hsf_ctl exit $?" ) &
( env OMP_NUM_THREADS=11 ATM_RAD_TOPO=0 ATM_UBUD_BALANCE=1 ATM_HYDRO_SPLIT=1.0 \
      ../cli/atm_hsf config_hsf_on.xml  > hsf_on.log  2>&1; echo "hsf_on exit $?" ) &
wait
echo "NaN ctl $(grep -c 'NaN/Inf DETECTED' hsf_ctl.log)  on $(grep -c 'NaN/Inf DETECTED' hsf_on.log)  $(date +%H:%M)"
touch HSF600_DONE
