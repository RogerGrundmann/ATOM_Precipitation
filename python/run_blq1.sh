#!/bin/bash
# 2026-10-02: BL-Q probe -- what sets the marine BL humidity (ocean |lat|<30); restart wb7 600 -> 640, working branch, cli/atm_blq.
# PRE-REGISTERED: moisture diffusivity K_q is >100x below the eddy viscosity at 500 m (moisture mixed molecularly, heat turbulently);
# E (~1000 mm/a) is deposited in levels 1..3 and leaves layer A mainly by resolved w, not by diffusion; RH falls sharply above ~130 m.
set -u; cd "$(dirname "$0")"; rm -f BLQ1_DONE
mkdir output_blq1 || { touch BLQ1_DONE; exit 1; }
[ -e config_blq1.xml ] && { touch BLQ1_DONE; exit 1; }
ln -s ../output_wb7/atm_restart_0Ma_600.bin output_blq1/atm_restart_0Ma_600.bin
sed -e "s#output_wb7/#output_blq1/#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>-1<#" \
    -e "s#<restart_from_iter>-1<#<restart_from_iter>600<#" -e "s#<nm>600<#<nm>640<#" config_wb7.xml > config_blq1.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_blq $(md5sum < ../cli/atm_blq | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 ../cli/atm_blq config_blq1.xml > blq1.log 2>&1
echo "blq1 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' blq1.log)  $(date +%H:%M)"
grep -a "\[BL-Q\]" blq1.log | tail -8
grep -a -i "land .*ocean" blq1.log | tail -1 | cut -c1-100
touch BLQ1_DONE
