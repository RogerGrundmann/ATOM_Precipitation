#!/bin/bash
# ATM_RAD_TOPO=1 -- the 600-iteration from-scratch arm the flip is judged on. ONE VARIABLE
# against output_el2 (the current default, PTOP 475, same binary cli/atm_el md5 d5116dd7, same
# env). Waits for the PTOP sweep to release the cores.
# WHY IT WAS REVERTED 2026-08-28: 622 cells at layer eps > 0.9, 270 > 0.99, land radiation
# 327.1 -> 306.4 W/m2 -- attributed to cwp_cap_col sharing one column total among fewer, thinner
# layers. Since then ATM_CLOUD_TAU_MAX (layer bound, 0.885-0.892 max eps on both branches) is
# on and the cap is disabled. MLR is NOT purely diagnostic under mode 5: apply_co2_perturbation
# feeds MLR(CO2) - MLR(280) into t_eq, so the climate can move over land.
# PRE-REGISTERED:
#   1. max epsilon stays <= ~0.89 (the tau ceiling); ZERO cells > 0.99. If not, the blocker stands.
#   2. [co2-perturb DIAG] surface shift changes over land only; ocean identical in pattern.
#   3. clear-sky / all-sky OLR and cloud LW forcing move < ~1 W/m2 globally.
#   4. precipitation, r, sigma, bands: null to ~1 % (MLR reaches t only via the CO2 increment).
#   5. stability: exit 0, zero NaN through 155/357/483; max u/v/w and mean KE within noise of el2.
set -u; cd "$(dirname "$0")"; rm -f RT600_DONE
until [ -f PE600_DONE ]; do sleep 60; done
mkdir -p output_rt600
env OMP_NUM_THREADS=20 ATM_MC_CAP_DIAG=1 ATM_MC_EVAP_LIMIT=1 ATM_MC_T_NDIM=1 ATM_RAD_TOPO=1 \
    ../cli/atm_el config_rt600.xml > rt600.log 2>&1
echo "rt600 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' rt600.log)  $(date +%H:%M)"
touch RT600_DONE
