#!/bin/bash
# ATM_HYDRO_SPLIT -- the hydrostatic / non-hydrostatic pressure split (AtmHydroSplit.h).
# Binary cli/atm_hs (md5 6f8c752b). ⚠ It was built from a working tree that also carries the
# UNCOMMITTED ATM_RAD_TOPO default flip, so EVERY arm sets ATM_RAD_TOPO=0 to stay on the committed
# default. cli/atm_el (md5 d5116dd7) predates the 2026-09-22 flip, so its arm sets
# MC_EVAP_LIMIT=1 MC_T_NDIM=1 explicitly, as output_el2 did.
#
# PART 1, BYTE CHECK (1 thread, nm = 20 from scratch):
#   hsv_old  cli/atm_el, MC flags              = shipped branch
#   hsv_off  cli/atm_hs, split unset           -- MUST equal hsv_old (RUN_CONFIG.txt: banner only)
#   hsv_on   cli/atm_hs, ATM_HYDRO_SPLIT=1     -- CONTROL, MUST differ
# PART 2, RESTART PAIR, 600 -> 700 from output_el2/atm_restart_0Ma_600.bin (buoyancy ramp = 1),
# 4 threads each, concurrent, ATM_UBUD_BALANCE=1 on both:
#   hs_ctl   split unset
#   hs_1     ATM_HYDRO_SPLIT=1.0
# PRE-REGISTERED (written before the run):
#   1. rhs_u: ubud_buoy rms ~3.9 (the consistent b), corr(pgf,buoy) ~ -1 and cancellation ~ 1,
#      because dp_hb/dr = b exactly; rms NET rhs_u NOT larger than the control's ~0.9. The
#      2026-09-11 =1 branch read net 0.796 of the buoyancy; anything near that is a FAIL.
#   2. max|u| does not grow: decays or stays level like hs_ctl. Growth like the =1 branch
#      (+2.8e-03 m/s per iteration) is a FAIL.
#   3. v-budget 20-70 deg above 3 km: median |pgf|/|cor| rises from ~1e-4 to O(0.3-3); band p05
#      ageostrophic residual falls well below 0.997 (ATM_HYDRO_PGF gave 0.079,
#      BUOY_CONSISTENT 0.437 at 100 iterations).
#   4. Velocity, jet and precipitation: null to ~1 % (1/f = 66 500 iterations; this is 100).
#   5. exit 0, zero NaN.
set -u; cd "$(dirname "$0")"; rm -f HS_DONE
. ./verify_lib.sh
rm -rf output_hsv_old output_hsv_off output_hsv_on
mkdir -p output_hsv_old output_hsv_off output_hsv_on output_hs_ctl output_hs_1
for d in output_hs_ctl output_hs_1; do cp -n output_el2/atm_restart_0Ma_600.bin $d/; done
arm hsv_old ../cli/atm_el config_hsv_old.xml ATM_MC_EVAP_LIMIT=1 ATM_MC_T_NDIM=1
arm hsv_off ../cli/atm_hs config_hsv_off.xml ATM_RAD_TOPO=0
arm hsv_on  ../cli/atm_hs config_hsv_on.xml  ATM_RAD_TOPO=0 ATM_HYDRO_SPLIT=1
( env OMP_NUM_THREADS=4 ATM_RAD_TOPO=0 ATM_UBUD_BALANCE=1 ATM_MC_CAP_DIAG=1 \
      ../cli/atm_hs config_hs_ctl.xml > hs_ctl.log 2>&1; echo "hs_ctl exit $?" > .hs_ctl.exit ) &
( env OMP_NUM_THREADS=4 ATM_RAD_TOPO=0 ATM_UBUD_BALANCE=1 ATM_MC_CAP_DIAG=1 ATM_HYDRO_SPLIT=1.0 \
      ../cli/atm_hs config_hs_1.xml > hs_1.log 2>&1; echo "hs_1 exit $?" > .hs_1.exit ) &
wait_arms
cmp_dirs    hsv_off hsv_old "A  OFF branch (atm_hs unset vs atm_el)"
want_differ hsv_on  hsv_off "C  CONTROL split on vs off -- MUST differ"
wait
cat .hs_ctl.exit .hs_1.exit
echo "NaN hs_ctl $(grep -c 'NaN/Inf DETECTED' hs_ctl.log)  hs_1 $(grep -c 'NaN/Inf DETECTED' hs_1.log)"
touch HS_DONE
