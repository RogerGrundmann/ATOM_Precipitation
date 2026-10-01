#!/bin/bash
# 2026-10-01: TRANSPORT/DIFFUSION leak probe -- [CWB-TD] splits the +208 mm/a remainder into advection as applied
# (minmod) / centred / -q div(u) / diffusion by level. Restart 600 -> 640 from output_rhst_115, current working branch
# (incl. SNOW_DEP_FLUX + MC_COND_DEBIT + MC_Q_NDIM), cli/atm_td (-O2), 8 threads.
# PRE-REGISTERED: -q div(u) carries most of the +208; limiter part tens of mm/a; ocean i0+1 diffusion ~0;
# residual vs the remainder ~0 (the self-check).
set -u; cd "$(dirname "$0")"; rm -f TD1_DONE
mkdir output_td1 || { touch TD1_DONE; exit 1; }
[ -e config_td1.xml ] && { touch TD1_DONE; exit 1; }
ln -s ../output_rhst_115/atm_restart_0Ma_600.bin output_td1/atm_restart_0Ma_600.bin
sed -e "s#output_rhst_115/#output_td1/#" \
    -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>-1<#" \
    -e "s#<restart_from_iter>-1<#<restart_from_iter>600<#" \
    -e "s#<nm>600<#<nm>640<#" config_rhst_115.xml > config_td1.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_td1 $(md5sum < ../cli/atm_td | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_CWB_DIAG=1 \
    ../cli/atm_td config_td1.xml > td1.log 2>&1
echo "td1 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' td1.log)  $(date +%H:%M)"
touch TD1_DONE
