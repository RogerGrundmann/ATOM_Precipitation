#!/bin/bash
# OFF-BRANCH BYTE CHECK for ATM_SNOW_WINDOW (default 0) and ATM_SNOW_DIAG (print-only), 2026-09-27.
# 1 thread, nm = 20 from scratch. old = cli/atm_rkp (38694a9 source), new = cli/atm_snw. NOT QUEUED -- launch by hand.
#   A   new clean == old
#   D   new + ATM_SNOW_DIAG=1 == new clean  (print-only: every written file identical, the tables go to the log)
#   C1  SNOW_WINDOW=1 MUST differ from new clean;   C2  SNOW_WINDOW=2 MUST differ
set -u; cd "$(dirname "$0")"; rm -f VSNW_VERIFY_DONE
. ./verify_lib.sh
for t in old new w1 w2 dg; do mkdir output_vsnw_$t || { touch VSNW_VERIFY_DONE; exit 1; }; done
arm vsnw_old ../cli/atm_rkp config_vsnw_old.xml
arm vsnw_new ../cli/atm_snw config_vsnw_new.xml
arm vsnw_dg  ../cli/atm_snw config_vsnw_dg.xml ATM_SNOW_DIAG=1
arm vsnw_w1  ../cli/atm_snw config_vsnw_w1.xml ATM_SNOW_WINDOW=1
arm vsnw_w2  ../cli/atm_snw config_vsnw_w2.xml ATM_SNOW_WINDOW=2
wait_arms
cmp_dirs    vsnw_new vsnw_old "A  OFF BRANCH (new clean vs old)"
cmp_dirs    vsnw_dg  vsnw_new "D  SNOW_DIAG print-only (diag vs new clean)"
want_differ vsnw_w1  vsnw_new "C1 CONTROL SNOW_WINDOW=1 -- MUST differ"
want_differ vsnw_w2  vsnw_new "C2 CONTROL SNOW_WINDOW=2 -- MUST differ"
grep -h "SNOW DIAG" vsnw_dg.log | tail -9
touch VSNW_VERIFY_DONE
