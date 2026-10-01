#!/bin/bash
# 2026-10-01: ATM_MC_Q_NDIM probe -- RK4 applies MC_q at L/u_0. Restart 600 -> 640 from output_rhst_115 on the CURRENT
# working branch (now incl. SNOW_DEP_FLUX + MC_COND_DEBIT), cli/atm_qnd (-O2), 8 threads.
# PRE-REGISTERED: [MC-Q] applied == raw (ratio 1.0000) ~ -P_conv; CWB "rest = MC_q" moves from ~0.713 to 1.0 x raw.
set -u; cd "$(dirname "$0")"; rm -f QND1_DONE
mkdir output_qnd1 || { touch QND1_DONE; exit 1; }
[ -e config_qnd1.xml ] && { touch QND1_DONE; exit 1; }
ln -s ../output_rhst_115/atm_restart_0Ma_600.bin output_qnd1/atm_restart_0Ma_600.bin
sed -e "s#output_rhst_115/#output_qnd1/#" \
    -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>-1<#" \
    -e "s#<restart_from_iter>-1<#<restart_from_iter>600<#" \
    -e "s#<nm>600<#<nm>640<#" config_rhst_115.xml > config_qnd1.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_qnd1 $(md5sum < ../cli/atm_qnd | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 ATM_MC_CAP_DIAG=1 ATM_CWB_DIAG=1 ATM_MC_Q_NDIM=1 \
    ../cli/atm_qnd config_qnd1.xml > qnd1.log 2>&1
echo "qnd1 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' qnd1.log)  $(date +%H:%M)"
touch QND1_DONE
