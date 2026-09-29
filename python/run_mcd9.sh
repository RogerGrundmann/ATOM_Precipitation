#!/bin/bash
# 2026-09-29: CONV-DEAD -- why does MoistConvection rain ~0.15 mm/a? ATM_MC_DIAG census + generation/evaporation
# budget, 600 -> 640 restarts (atm nm = TOTAL, so nm 640), cli/atm_sdr, 2 x 8 threads.
#   mcd9_stm  from stm_2 (energy-conserving working branch: closure sync 2, moisture filter + vertical T filter off,
#             UPWIND + SNOW_WINDOW=2)
#   mcd9_ra   from ra_10 (shipped branch, same UPWIND + SNOW_WINDOW=2)
set -u; cd "$(dirname "$0")"; rm -f MCD9_DONE
for t in stm ra; do mkdir output_mcd9_$t || { touch MCD9_DONE; exit 1; }
  [ -e config_mcd9_$t.xml ] && { echo "config exists"; touch MCD9_DONE; exit 1; }
  sed -e "s#output_o2val/#output_mcd9_$t/#" -e "s#<restart_from_iter>-1<#<restart_from_iter>600<#" -e "s#<nm>600<#<nm>640<#" \
      config_o2val.xml > config_mcd9_$t.xml; done
cp output_stm_2/atm_restart_0Ma_600.bin output_mcd9_stm/; cp output_ra_10/atm_restart_0Ma_600.bin output_mcd9_ra/
BR="ATM_MC_DIAG=1 ATM_MC_CAP_DIAG=1 ATM_CWB_DIAG=1 ATM_PRECIP_UPWIND=1 ATM_SNOW_WINDOW=2"
echo "start $(date +%H:%M)"
( env OMP_NUM_THREADS=8 $BR ATM_WATER_CLOSURE=1 ATM_DAMP_Q_VERT=0 ATM_DAMP_Q_HORIZ=0 ATM_DAMP_T_VERT=0 \
    ../cli/atm_sdr config_mcd9_stm.xml > mcd9_stm.log 2>&1; echo "mcd9_stm exit $? $(date +%H:%M)" ) &
( env OMP_NUM_THREADS=8 $BR ATM_RAIN_AREA=0.10 ../cli/atm_sdr config_mcd9_ra.xml > mcd9_ra.log 2>&1; echo "mcd9_ra exit $? $(date +%H:%M)" ) &
wait; touch MCD9_DONE
