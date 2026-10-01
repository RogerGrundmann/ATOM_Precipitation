#!/bin/bash
# 2026-10-01: DIFF-LEAK probe -- [CWB-TD] splits the moisture diffusion (+220 mm/a, +261 above i0+1 in td1) by term
# (radial 2nd / radial metric 2/r / merid 2nd / merid metric cot / zonal) and the interior part by height and |lat|.
# Restart 600 -> 640 from output_rhst_115, current working branch, cli/atm_td2 (-O2), 8 threads.
# PRE-REGISTERED (from the code, not measured): a metric term built on the shell coordinate rm (1-2) rather than the
# Earth radius, or the one-sided lid stencil -- i.e. the top-3-levels row or a metric row carries most of it.
set -u; cd "$(dirname "$0")"; rm -f TD2_DONE
mkdir output_td2 || { touch TD2_DONE; exit 1; }
[ -e config_td2.xml ] && { touch TD2_DONE; exit 1; }
ln -s ../output_rhst_115/atm_restart_0Ma_600.bin output_td2/atm_restart_0Ma_600.bin
sed -e "s#output_rhst_115/#output_td2/#" \
    -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>-1<#" \
    -e "s#<restart_from_iter>-1<#<restart_from_iter>600<#" \
    -e "s#<nm>600<#<nm>640<#" config_rhst_115.xml > config_td2.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_td2 $(md5sum < ../cli/atm_td2 | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_CWB_DIAG=1 \
    ../cli/atm_td2 config_td2.xml > td2.log 2>&1
echo "td2 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' td2.log)  $(date +%H:%M)"
touch TD2_DONE
