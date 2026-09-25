#!/bin/bash
# PRINT-ONLY CHECK for ATM_CWB_BANDS (2026-09-25). 1 thread, nm = 20 from scratch.
# old = cli/atm_qvd (6b49292), new = cli/atm_cwbb. All three arms must write IDENTICAL files (apart from
# RUN_CONFIG's path): the instrument prints and writes nothing. The on-arm must print the band table, and
# its bands (weighted by band weight) must reproduce the global rows -- checked by check_cwbb.py.
set -u; cd "$(dirname "$0")"; rm -f CWBB_VERIFY_DONE
. ./verify_lib.sh
for t in old new on; do mkdir output_cwbb_$t || exit 1; done
arm cwbb_old ../cli/atm_qvd  config_cwbb_old.xml ATM_CWB_DIAG=1
arm cwbb_new ../cli/atm_cwbb config_cwbb_new.xml ATM_CWB_DIAG=1
arm cwbb_on  ../cli/atm_cwbb config_cwbb_on.xml  ATM_CWB_DIAG=1 ATM_CWB_BANDS=1
wait_arms
cmp_dirs cwbb_new cwbb_old "A  OFF (new vs old, CWB on, bands off)"
cmp_dirs cwbb_on  cwbb_new "B  ON  (bands on vs off) -- print-only, MUST MATCH"
diff <(grep -v 'CWB-BANDS' cwbb_on.log | grep '\[CWB\]') <(grep '\[CWB\]' cwbb_new.log) > /dev/null && echo "C  [CWB] rows identical with bands on: PASS" || echo "C  [CWB] rows differ: FAIL"
echo "band rows printed: $(grep -c 'CWB-BANDS' cwbb_on.log)"
touch CWBB_VERIFY_DONE
