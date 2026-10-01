#!/bin/bash
# 2026-10-01: MC-Q-LEAK knob probe -- ATM_MC_COND_DEBIT=1 (c_u = g_p + e_l). Restart 600 -> 640 from output_rhst_115,
# working branch, cli/atm_cdb (-O2), 8 threads. Compare with mcq.log / cuz.log (same window, knob off).
# PRE-REGISTERED: [MC-Q] c_u - e_l - g_p = 0; MC_q column ~ -P_conv (~ -300) unless MCq_max binds (printed);
# latent heating ~ +46 W/m2 in convecting columns, so MCt_max truncation rises; P_conv itself unchanged at 40 iters.
set -u; cd "$(dirname "$0")"; rm -f CDB1_DONE
mkdir output_cdb1 || { touch CDB1_DONE; exit 1; }
[ -e config_cdb1.xml ] && { touch CDB1_DONE; exit 1; }
ln -s ../output_rhst_115/atm_restart_0Ma_600.bin output_cdb1/atm_restart_0Ma_600.bin
sed -e "s#output_rhst_115/#output_cdb1/#" \
    -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>-1<#" \
    -e "s#<restart_from_iter>-1<#<restart_from_iter>600<#" \
    -e "s#<nm>600<#<nm>640<#" config_rhst_115.xml > config_cdb1.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_cdb1 $(md5sum < ../cli/atm_cdb | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 ATM_MC_CAP_DIAG=1 ATM_CWB_DIAG=1 ATM_MC_COND_DEBIT=1 \
    ../cli/atm_cdb config_cdb1.xml > cdb1.log 2>&1
echo "cdb1 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' cdb1.log)  $(date +%H:%M)"
touch CDB1_DONE
