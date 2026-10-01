#!/bin/bash
# 2026-10-01: let the convective parcel condense -- ATM_MC_SGZ=1 (s = cp*T + g*z) + ATM_MC_ENTR=1e-4 (Tiedtke deep
# entrainment, shipped 2e-3) on the CURRENT working branch (incl. ATM_MC_QC_DETRAIN). Restart 600 -> 640 from
# output_rhst_115, cli/atm_qc (-O2), 8 threads. Baseline: qc1 (same binary, same window, = the current working branch).
# PRE-REGISTERED: the [MC-Q] census finds supersaturated parcels where the condensation step runs (was 0) and c_u > 0;
# [MC-QC] condensate grows above cloud base (growth p50 > 1, was 0); g_p and P_conv rise from qc1's 179 / 101;
# MC_q = -P_conv still; no NaN.
set -u; cd "$(dirname "$0")"; rm -f PCD1_DONE
mkdir output_pcd1 || { touch PCD1_DONE; exit 1; }
[ -e config_pcd1.xml ] && { touch PCD1_DONE; exit 1; }
ln -s ../output_rhst_115/atm_restart_0Ma_600.bin output_pcd1/atm_restart_0Ma_600.bin
sed -e "s#output_rhst_115/#output_pcd1/#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>-1<#" \
    -e "s#<restart_from_iter>-1<#<restart_from_iter>600<#" -e "s#<nm>600<#<nm>640<#" config_rhst_115.xml > config_pcd1.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_qc $(md5sum < ../cli/atm_qc | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 ATM_MC_SGZ=1 ATM_MC_ENTR=1e-4 ../cli/atm_qc config_pcd1.xml > pcd1.log 2>&1
echo "pcd1 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' pcd1.log)  $(date +%H:%M)  $(grep -ao 'MC_SGZ=[^ ]*\|MC_ENTR=[^ ]*' output_pcd1/RUN_CONFIG.txt | tr '\n' ' ')"
for t in qc1 pcd1; do
  echo "== $t"; grep -a "\[MC-Q\] updraftRecurrence" $t.log | tail -1 | cut -c1-330
  grep -a "\[MC-QC\]" $t.log | tail -1; grep -a "\[MC DIAG\] column budget\|\[MC DIAG\] cloud base\|\[MC DIAG\] as shares" $t.log | tail -3 | cut -c1-260
  grep -a "\[MC-Q\] closure" $t.log | tail -1 | cut -c1-140; grep -a "Precip mean\|P_conv mean" $t.log | tail -2 | cut -c1-70
  grep -a "max temperature\|max w-component\|max u-component" $t.log | tail -3 | cut -c1-80
done
touch PCD1_DONE
