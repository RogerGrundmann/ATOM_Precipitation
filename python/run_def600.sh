#!/bin/bash
# THE CURRENT DEFAULT, 600 FROM SCRATCH, 24 threads (2026-09-24, at the user's instruction).
# Binary cli/atm_def600 (= atm_msflip: HEAD efcbdaf + ATM_MC_S_NDIM default 1). No knob set, so every
# default of the day is on: ATM_RAD_TOPO, ATM_HYDRO_SPLIT, ATM_SEAM_PERIODIC, ATM_WATER_CLOSURE (+ fade 2,
# RK sync, mass filter), the evaporation stride fix, ATM_RH_MIN_PTOP 482, ATM_MC_S_NDIM. Print-only:
# ATM_CWB_DIAG, ATM_MC_CAP_DIAG. Waits for the MC_S_NDIM byte check and the last ocean run.
# PRE-REGISTERED (the three questions this run was owed):
#   1. SEAM: max|v| returns to ~2.3 m/s in the upper troposphere, NOT ~26 m/s at 11N 1E 236 m (hsf_ctl).
#   2. PRECIP on the closure default: global mean, bands, land/ocean, r, sigma at 600 vs NASA; whether it
#      drifts (PW and precip by checkpoint). Expect BELOW 978 (closure dried the column 10 % in 4 s; PTOP
#      482 was fitted without it) -- the size of the miss sets the PTOP re-fit.
#   3. MC_t cap: 0 % truncated throughout; P_conv revives or not (0.17 mm/a on the old branch).
#   Also: jet (max|w|) keeps growing or levels off; p_dyn non-hydrostatic stays well below 3.0;
#   CWB leapfrog_reset 0, NET small; exit 0, zero NaN through 155/357/483.
set -u; cd "$(dirname "$0")"; rm -f DEF600_DONE
until [ -f MSF_VERIFY_DONE ] && [ -f BN600_DONE ]; do sleep 30; done
if [ "$(grep -c 'PASS' run_verify_msflip.out)" -lt 1 ]; then echo "byte check not complete -- NOT started"; touch DEF600_DONE; exit 1; fi
echo "start $(date +%H:%M)"
env OMP_NUM_THREADS=24 ATM_CWB_DIAG=1 ATM_MC_CAP_DIAG=1 ../cli/atm_def600 config_def600.xml > def600.log 2>&1
echo "def600 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' def600.log)  $(date +%H:%M)"
touch DEF600_DONE
