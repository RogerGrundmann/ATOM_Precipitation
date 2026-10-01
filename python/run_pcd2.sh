#!/bin/bash
# 2026-10-01: ATM_MC_ED_ABOVE_BASE probe -- the pcd1 setting (working branch + ATM_MC_SGZ=1 + ATM_MC_ENTR=1e-4) with e_d
# kept strictly above cloud base. Restart 600 -> 640 from output_rhst_115, cli/atm_eab (-O2), 8 threads. Baseline pcd1.
# PRE-REGISTERED: surviving convective rain 8.5 % -> ~25-30 % of generation, P_conv 109 -> ~320-380 mm/a; e_d share
# falls by its below-base part (16-24 %); g_p roughly unchanged at first; MC_q = -P_conv; no NaN.
set -u; cd "$(dirname "$0")"; rm -f PCD2_DONE
mkdir output_pcd2 || { touch PCD2_DONE; exit 1; }
[ -e config_pcd2.xml ] && { touch PCD2_DONE; exit 1; }
ln -s ../output_rhst_115/atm_restart_0Ma_600.bin output_pcd2/atm_restart_0Ma_600.bin
sed -e "s#output_rhst_115/#output_pcd2/#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>-1<#" \
    -e "s#<restart_from_iter>-1<#<restart_from_iter>600<#" -e "s#<nm>600<#<nm>640<#" config_rhst_115.xml > config_pcd2.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_eab $(md5sum < ../cli/atm_eab | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 ATM_MC_SGZ=1 ATM_MC_ENTR=1e-4 ATM_MC_ED_ABOVE_BASE=1 ../cli/atm_eab config_pcd2.xml > pcd2.log 2>&1
echo "pcd2 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' pcd2.log)  $(date +%H:%M)  $(grep -ao 'MC_SGZ=[^ ]*\|MC_ENTR=[^ ]*\|MC_ED_ABOVE_BASE=[^ ]*' output_pcd2/RUN_CONFIG.txt | tr '\n' ' ')"
for t in pcd1 pcd2; do
  echo "== $t"; grep -a "\[MC-Q\] updraftRecurrence" $t.log | tail -1 | cut -c1-330
  grep -a "\[MC-QC\]" $t.log | tail -1; grep -a "\[MC DIAG\] column budget\|\[MC DIAG\] cloud base\|\[MC DIAG\] as shares" $t.log | tail -3 | cut -c1-260
  grep -a "\[MC-Q\] closure" $t.log | tail -1 | cut -c1-140; grep -a "Precip mean\|P_conv mean" $t.log | tail -2 | cut -c1-70
  grep -a "max temperature\|max w-component\|max u-component" $t.log | tail -3 | cut -c1-80
done
touch PCD2_DONE
