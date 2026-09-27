#!/bin/bash
# RAIN PASS DIAGNOSTIC (2026-09-27; NOT QUEUED -- launch by hand). Both builds -O2 (the tree is -O2 since dec8cec).
# 1) byte check, 1 thread, nm = 20 from scratch: new = cli/atm_ocfix (4aec4dd + ATM_RAIN_PASS_DIAG + OneCat S_au fix),
#    old = cli/atm_O2h (4aec4dd at -O2).  A  new clean == old (default branch: TwoCat, OneCat fix cannot reach it)
#    D  new + ATM_RAIN_PASS_DIAG=1 == new clean (print-only)
# 2) only if both pass: 600 -> 620 restarts from output_sgzb_ctl's checkpoint, 2 x ${NT:-4} threads, the diagnostic on:
#    rpd_w0 shipped, rpd_w1 ATM_SNOW_WINDOW=1 (the arm whose rain collapsed in snw40).
# QUESTION: do the three iter_prec passes agree? On the probe column (j=90, k=180) the ground rain cycles
# 1.97 -> 0.0 -> 14.7 mm/d every call. If |pass3 - pass2| is a large fraction of the band's rain, the reported
# precipitation is set by the pass count, and the snow-window warm repair cannot be judged until that is repaired.
set -u; cd "$(dirname "$0")"; rm -f RPD20_DONE
. ./verify_lib.sh
for t in old new dg; do mkdir output_vrpd_$t || { touch RPD20_DONE; exit 1; }; done
arm vrpd_old ../cli/atm_O2h   config_vrpd_old.xml
arm vrpd_new ../cli/atm_ocfix config_vrpd_new.xml
arm vrpd_dg  ../cli/atm_ocfix config_vrpd_dg.xml ATM_RAIN_PASS_DIAG=1
wait_arms
cmp_dirs vrpd_new vrpd_old "A  OFF BRANCH (new clean vs old, both -O2)"
cmp_dirs vrpd_dg  vrpd_new "D  RAIN_PASS_DIAG print-only"
BAD=$(sed -n '1,/D  RAIN_PASS_DIAG/p' run_rpd20.out 2>/dev/null | grep 'DIFFERS:' | grep -vc 'RUN_CONFIG.txt')
if [ "$BAD" != 0 ]; then echo "byte check did not pass -- restarts NOT started"; touch RPD20_DONE; exit 1; fi
grep -h "RAIN PASS DIAG" vrpd_dg.log | tail -9
for t in w0 w1; do mkdir output_rpd_$t || { touch RPD20_DONE; exit 1; }; cp output_sgzb_ctl/atm_restart_0Ma_600.bin output_rpd_$t/; done
export OMP_NUM_THREADS=${NT:-4}
( env ATM_RAIN_PASS_DIAG=1 ATM_SNOW_WINDOW=0 ../cli/atm_ocfix config_rpd_w0.xml > rpd_w0.log 2>&1; echo "rpd_w0 exit $?" ) &
( env ATM_RAIN_PASS_DIAG=1 ATM_SNOW_WINDOW=1 ../cli/atm_ocfix config_rpd_w1.xml > rpd_w1.log 2>&1; echo "rpd_w1 exit $?" ) &
wait
for t in w0 w1; do echo "== rpd_$t"; grep "RAIN PASS DIAG" rpd_$t.log | tail -9; done
touch RPD20_DONE
