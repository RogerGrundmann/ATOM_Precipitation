#!/bin/bash
# ATM_MC_QVD (downdraft humidity) on top of SGZ + ENTR + BASE_SAT, 2026-09-25. 600 from scratch, cli/atm_qvd,
# 2 arms x 6 threads, queued behind the sgzb arms, gated on run_verify_qvd.sh. Control = bst_sgeb (QVD=0).
#   qvd1  ATM_MC_QVD=1   downdraft = 0.5*(c + 0.98 q_sat)   (the equal-parts mix, presumably intended)
#   qvd2  ATM_MC_QVD=2   downdraft saturated 0.98 q_sat     (Tiedtke)
# PRE-REGISTERED (bst_sgeb ~iter 120: g_p 1.05e+04 mm/a, e_d takes 100 %, P_conv 0.00): e_d's share falls
# (qvd2 most: deficit 0.02 q_sat instead of ~0.5), P_conv(ground) leaves zero -- the first working
# convective rain in this tree; total precip may RISE (less re-evaporation feeding the stratiform path is
# not guaranteed to balance) -- score r, sigma, bands, drift vs bst_sgeb; 0-15 band is where it must show.
# Stability: exit 0, zero NaN through 155/357/483; MC_t / MC_q truncation.
set -u; cd "$(dirname "$0")"; rm -f QVD600_DONE
until [ -f SGZB600_DONE ] && [ -f QVD_VERIFY_DONE ]; do sleep 30; done
BAD=$(sed -n '1,/A  OFF BRANCH/p' run_verify_qvd.out | grep 'DIFFERS:' | grep -vc 'RUN_CONFIG.txt')
if [ "$BAD" != 0 ] || ! grep -q "A  OFF BRANCH" run_verify_qvd.out || ! grep -q "CONTROL.*PASS" run_verify_qvd.out; then
    echo "byte check did not pass -- NOT started"; cat run_verify_qvd.out; touch QVD600_DONE; exit 1; fi
for t in qvd1 qvd2; do mkdir output_$t || { touch QVD600_DONE; exit 1; }; done
echo "start $(date +%H:%M)"
run(){ local tag=$1; shift
       ( env OMP_NUM_THREADS=6 ATM_CWB_DIAG=1 ATM_MC_CAP_DIAG=1 ATM_MC_SGZ=1 ATM_MC_ENTR=1.0e-4 ATM_MC_BASE_SAT=1 "$@" \
             ../cli/atm_qvd config_$tag.xml > $tag.log 2>&1
         echo "$tag exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $tag.log)  $(date +%H:%M)" ) & }
run qvd1 ATM_MC_QVD=1
run qvd2 ATM_MC_QVD=2
wait
touch QVD600_DONE
