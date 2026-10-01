#!/bin/bash
# 2026-10-01: RAIN-CONV probe -- do the ice scheme's RATE ARRAYS carry the water its own ground FLUX does?
# Print-only: ATM_SNOW_DIAG rows 8-13 (S_s_dep deposition/sublimation, S_i_cri, S_r_cri, SUM S_r+S_s, flux removed,
# identity residual) and the CWB 'rest = MC_q + transport' split. Restart 600 -> 640 from output_rhst_115 (working
# branch + ATM_RH_STORM=1.15), cli/atm_sprb (-O2), 8 threads.
# PRE-REGISTERED: sublimation ~ -5.8e3 mm/a global, mostly 0-35 deg; identity residual ~0; MC_q column != -P_conv (~ -296).
set -u; cd "$(dirname "$0")"; rm -f SPRB_DONE
mkdir output_sprb || { touch SPRB_DONE; exit 1; }
[ -e config_sprb.xml ] && { touch SPRB_DONE; exit 1; }
ln -s ../output_rhst_115/atm_restart_0Ma_600.bin output_sprb/atm_restart_0Ma_600.bin
sed -e "s#output_rhst_115/#output_sprb/#" \
    -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>-1<#" \
    -e "s#<restart_from_iter>-1<#<restart_from_iter>600<#" \
    -e "s#<nm>600<#<nm>640<#" config_rhst_115.xml > config_sprb.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_sprb $(md5sum < ../cli/atm_sprb | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_SNOW_DIAG=1 ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_CWB_BANDS=1 ATM_SR_DIAG=1 \
    ../cli/atm_sprb config_sprb.xml > sprb.log 2>&1
echo "sprb exit $?  NaN $(grep -c 'NaN/Inf DETECTED' sprb.log)  $(date +%H:%M)"
touch SPRB_DONE
