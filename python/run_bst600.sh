#!/bin/bash
# B.10c: ATM_MC_BASE_SAT on top of the B.10b pair (2026-09-25). 600 from scratch, cli/atm_bst, 6 threads,
# queued behind qvt_on, gated on run_verify_bst.sh. Pairs with sgzb_sge (same, BASE_SAT off) and sgzb_ctl.
#   bst_sgeb  ATM_MC_SGZ=1 ATM_MC_ENTR=1e-4 ATM_MC_BASE_SAT=1
# PRE-REGISTERED (sgzb_sge @120, 87E: parcel -3.9/-4.1/-5.7 K at 2-5/5-8/8-12 km; unsaturated for ~0.9 km above
# base): the dry segment vanishes -- parcel saturated from base (q_v_u >= q_sat), s_u - s at 2-5 km within
# ~+-1 K of zero or positive, the aloft deficit shrinks by ~4 K; g_p rises further (sgzb_sge 1682 mm/a early);
# THE open question is whether P_conv leaves 0.0. Stability: exit 0, zero NaN through 155/357/483; MC_t trunc.
set -u; cd "$(dirname "$0")"; rm -f BST600_DONE
until [ -f QVT100_DONE ] && [ -f BST_VERIFY_DONE ]; do sleep 30; done
BAD=$(sed -n '1,/A  OFF BRANCH/p' run_verify_bst.out | grep 'DIFFERS:' | grep -vc 'RUN_CONFIG.txt')
if [ "$BAD" != 0 ] || ! grep -q "A  OFF BRANCH" run_verify_bst.out || ! grep -q "CONTROL.*PASS" run_verify_bst.out; then
    echo "byte check did not pass -- NOT started"; cat run_verify_bst.out; touch BST600_DONE; exit 1; fi
mkdir output_bst_sgeb || { touch BST600_DONE; exit 1; }
echo "start $(date +%H:%M)"
env OMP_NUM_THREADS=6 ATM_CWB_DIAG=1 ATM_MC_CAP_DIAG=1 ATM_MC_SGZ=1 ATM_MC_ENTR=1.0e-4 ATM_MC_BASE_SAT=1 \
    ../cli/atm_bst config_bst_sgeb.xml > bst_sgeb.log 2>&1
echo "bst_sgeb exit $?  NaN $(grep -c 'NaN/Inf DETECTED' bst_sgeb.log)  $(date +%H:%M)"
touch BST600_DONE
