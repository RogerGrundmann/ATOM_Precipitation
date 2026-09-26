#!/bin/bash
# STEP B (2026-09-26): ATM_MC_ALF1 on the qvd2 stack (SGZ + ENTR 1e-4 + BASE_SAT + QVD=2). 600 from scratch,
# cli/atm_alf, 24 threads, SERIAL, queued behind run_qh100.sh (step A) so A is not disturbed, gated on
# run_verify_dalf.sh. Compare qvd2 (alf_1 = 0.05) and sgzb_ctl (the corrected default).
#   alf_tie  ATM_MC_ALF1=5.44e-4   Tiedtke (1989), 92x below shipped
#   alf_mid  ATM_MC_ALF1=5.0e-3    halfway in log, 10x below shipped
# PRE-REGISTERED (qvd2 at 600: g_p 10943, e_d 18 %, e_p 82 %, P_conv 2.68 mm/a): e_p's share falls and P_conv
# rises by an order of magnitude or more -- the first convective rain carrying a real fraction of the total.
# Whether total precip, r, sigma and the 0-15 band improve is OPEN (every repaired-parcel arm so far was a
# climate null to slightly worse, sigma 2.42-2.44, drift +0.06/iter). Watch the drift 400-600 against
# sgzb_ctl's -0.01. Stability: exit 0, zero NaN through 155/357/483; MC_t / MC_q truncation.
set -u; cd "$(dirname "$0")"; rm -f ALF600_DONE
until [ -f QH100_DONE ] && [ -f DALF_VERIFY_DONE ]; do sleep 30; done
BAD=$(sed -n '1,/A  OFF BRANCH/p' run_verify_dalf.out | grep 'DIFFERS:' | grep -vc 'RUN_CONFIG.txt')
if [ "$BAD" != 0 ] || ! grep -q "A  OFF BRANCH" run_verify_dalf.out || ! grep -q "CONTROL.*PASS" run_verify_dalf.out; then
    echo "byte check did not pass -- NOT started"; cat run_verify_dalf.out; touch ALF600_DONE; exit 1; fi
for t in alf_tie alf_mid; do
  mkdir output_$t || { touch ALF600_DONE; exit 1; }
  echo "$t start $(date +%H:%M)"
  env OMP_NUM_THREADS=24 ATM_CWB_DIAG=1 ATM_MC_CAP_DIAG=1 ATM_MC_SGZ=1 ATM_MC_ENTR=1.0e-4 ATM_MC_BASE_SAT=1 ATM_MC_QVD=2 \
      ATM_MC_ALF1=$( [ $t = alf_tie ] && echo 5.44e-4 || echo 5.0e-3 ) \
      ../cli/atm_alf config_$t.xml > $t.log 2>&1
  echo "$t exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
done
touch ALF600_DONE
