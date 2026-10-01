#!/bin/bash
# 2026-10-01: ATM_MC_QC_DETRAIN probe pair -- the updraft condensate detrains its own q_c_u instead of the environment's
# cloud. Restart 600 -> 640 from output_rhst_115, CURRENT working branch (incl. the three MC-TV knobs), one binary
# cli/atm_qc (-O2), 2 x 8 threads concurrent: qc0 knob off, qc1 knob on. [MC-QC] prints the condensate growth above
# cloud base and the recurrence's M_u[i-1]/M_u[i].
# PRE-REGISTERED: qc0 growth p50 >> 1 (tens) with geometric-mean ratio > 1; qc1 growth ~<= 1 (no source while c_u = 0),
# g_p and P_conv fall sharply; MC_q column still = -P_conv (closure intact).
set -u; cd "$(dirname "$0")"; rm -f QC_DONE
for t in qc0 qc1; do
  mkdir output_$t || { touch QC_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch QC_DONE; exit 1; }
  ln -s ../output_rhst_115/atm_restart_0Ma_600.bin output_$t/atm_restart_0Ma_600.bin
  sed -e "s#output_rhst_115/#output_$t/#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>-1<#" \
      -e "s#<restart_from_iter>-1<#<restart_from_iter>600<#" -e "s#<nm>600<#<nm>640<#" config_rhst_115.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_qc $(md5sum < ../cli/atm_qc | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 ATM_MC_QC_DETRAIN=0 ../cli/atm_qc config_qc0.xml > qc0.log 2>&1 &
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 ATM_MC_QC_DETRAIN=1 ../cli/atm_qc config_qc1.xml > qc1.log 2>&1 &
wait
for t in qc0 qc1; do
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(grep -ao 'MC_QC_DETRAIN=[^ ]*' output_$t/RUN_CONFIG.txt | head -1)"
  grep -a "\[MC-QC\]" $t.log | tail -1; grep -a "\[MC DIAG\] column budget" $t.log | tail -1 | cut -c1-230
  grep -a "\[MC-Q\] closure" $t.log | tail -1 | cut -c1-140; grep -a "Precip mean\|P_conv mean" $t.log | tail -2 | cut -c1-70
done
touch QC_DONE
