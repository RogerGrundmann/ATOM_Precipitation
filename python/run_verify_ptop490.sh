#!/bin/bash
# BOTH-DIRECTIONS BYTE CHECK, 2026-09-23: ATM_RH_MIN_PTOP default 475 -> 490. 1 thread, nm = 20
# FROM SCRATCH (the knob is read only in the initial condition; a restart cannot see it).
#   ptv_old  cli/atm_bc2 (514e9c17, default 475), clean
#   ptv_off  cli/atm_pt490, ATM_RH_MIN_PTOP=475  -- MUST equal ptv_old
#   ptv_new  cli/atm_pt490, clean                -- the new default
#   ptv_on   cli/atm_bc2, ATM_RH_MIN_PTOP=490    -- MUST equal ptv_new
#   CONTROL  ptv_new vs ptv_off MUST differ. RUN_CONFIG.txt differs by path and '*' only.
set -u; cd "$(dirname "$0")"; rm -f PT490_VERIFY_DONE
. ./verify_lib.sh
for t in old off new on; do rm -rf output_ptv_$t; mkdir -p output_ptv_$t; done
arm ptv_old ../cli/atm_bc2   config_ptv_old.xml
arm ptv_off ../cli/atm_pt490 config_ptv_off.xml ATM_RH_MIN_PTOP=475
arm ptv_new ../cli/atm_pt490 config_ptv_new.xml
arm ptv_on  ../cli/atm_bc2   config_ptv_on.xml  ATM_RH_MIN_PTOP=490
wait_arms
cmp_dirs    ptv_off ptv_old "A  RESTORE"
cmp_dirs    ptv_on  ptv_new "B  FLIP"
want_differ ptv_new ptv_off "C  CONTROL"
touch PT490_VERIFY_DONE
