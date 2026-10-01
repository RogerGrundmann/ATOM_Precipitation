#!/bin/bash
# 2026-10-01: ATM_Q_DIFF_FLUX probe -- flux-form vertical moisture diffusion. Restart 600 -> 640 from output_rhst_115,
# current working branch + ATM_Q_DIFF_FLUX=1, cli/atm_qdf (-O2), 8 threads. Compare td3 (knob off, same window).
# PRE-REGISTERED: radial diffusion row ~0 (exact telescoping to the zero boundary fluxes), interior rows ~0,
# the transport/diffusion remainder falls from +206 to about the advection (-14); split residual still ~1.
set -u; cd "$(dirname "$0")"; rm -f TD4_DONE
mkdir output_td4 || { touch TD4_DONE; exit 1; }
[ -e config_td4.xml ] && { touch TD4_DONE; exit 1; }
ln -s ../output_rhst_115/atm_restart_0Ma_600.bin output_td4/atm_restart_0Ma_600.bin
sed -e "s#output_rhst_115/#output_td4/#" \
    -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>-1<#" \
    -e "s#<restart_from_iter>-1<#<restart_from_iter>600<#" \
    -e "s#<nm>600<#<nm>640<#" config_rhst_115.xml > config_td4.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_td4 $(md5sum < ../cli/atm_qdf | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_CWB_DIAG=1 ATM_Q_DIFF_FLUX=1 \
    ../cli/atm_qdf config_td4.xml > td4.log 2>&1
echo "td4 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' td4.log)  $(date +%H:%M)"
touch TD4_DONE
