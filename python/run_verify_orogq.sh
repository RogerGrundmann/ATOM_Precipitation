#!/bin/bash
# ATM_OROG_Q_MASS (B+4, new 2026-09-24, default 0). SHORT CHECKS ONLY.
# (1) off-branch byte check, 1 thread, nm = 20: boq_new (atm_oqm unset) MUST equal boq_old (atm_msn),
#     RUN_CONFIG.txt by the banner token; boq_on (=1) MUST differ.
# (2) 600 -> 620 restart from output_hsf_ctl, ATM_CWB_DIAG=1 ATM_OROG_Q_MASS=1 on today's default
#     (oqd); control is bcd (same default, knob off): the orographic_shapiro row must fall to ~0.
set -u; cd "$(dirname "$0")"; rm -f OQ_VERIFY_DONE
. ./verify_lib.sh
for t in old new on; do rm -rf output_boq_$t; mkdir -p output_boq_$t; done
arm boq_old ../cli/atm_msn config_boq_old.xml
arm boq_new ../cli/atm_oqm config_boq_new.xml
arm boq_on  ../cli/atm_oqm config_boq_on.xml ATM_OROG_Q_MASS=1
( env OMP_NUM_THREADS=3 ATM_CWB_DIAG=1 ATM_OROG_Q_MASS=1 ../cli/atm_oqm config_oqd.xml > oqd.log 2>&1 ) &
wait_arms; wait
cmp_dirs    boq_new boq_old "ATM A  OFF-BRANCH (atm_oqm unset vs atm_msn)"
diff <(sed 's#output_boq_new#X#g' output_boq_new/RUN_CONFIG.txt) <(sed 's#output_boq_old#X#g' output_boq_old/RUN_CONFIG.txt) | grep -o "OROG_Q_MASS=[^ ]*"
want_differ boq_on  boq_new "ATM C  CONTROL =1 vs unset -- MUST differ"
for f in bcd oqd; do echo "== $f  NaN=$(grep -c 'NaN/Inf DETECTED' $f.log)"
  grep "\[CWB\]" $f.log | grep -E "orographic|NET|unattributed|BC:phi" | tail -4 | cut -c1-110
  grep -E "model .* NASA" $f.log | tail -1 | cut -c1-120; done
touch OQ_VERIFY_DONE
