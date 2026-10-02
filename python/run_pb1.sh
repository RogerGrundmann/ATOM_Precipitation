#!/bin/bash
# 2026-10-02: MC-PB probe -- mixed-layer parcel vs the scheme base parcel; restart wb7 600 -> 640, working branch, cli/atm_pb.
# PRE-REGISTERED: the scheme CAPE is reproduced by CAPE_base; ocean<30 ML parcel theta_e is high (SST 27 C, RH ~0.8) so CAPE_ML >> CAPE_base over ocean -- the base parcel discards the boundary layer.
set -u; cd "$(dirname "$0")"; rm -f PB1_DONE
mkdir output_pb1 || { touch PB1_DONE; exit 1; }
[ -e config_pb1.xml ] && { touch PB1_DONE; exit 1; }
ln -s ../output_wb7/atm_restart_0Ma_600.bin output_pb1/atm_restart_0Ma_600.bin
sed -e "s#output_wb7/#output_pb1/#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>-1<#" \
    -e "s#<restart_from_iter>-1<#<restart_from_iter>600<#" -e "s#<nm>600<#<nm>640<#" config_wb7.xml > config_pb1.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_pb $(md5sum < ../cli/atm_pb | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 ../cli/atm_pb config_pb1.xml > pb1.log 2>&1
echo "pb1 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' pb1.log)  $(date +%H:%M)"
grep -a "\[MC-PB\]" pb1.log | tail -4
grep -a "\[MC-LO\]" pb1.log | tail -6
grep -a -i "land .*ocean" pb1.log | tail -1 | cut -c1-100
touch PB1_DONE
