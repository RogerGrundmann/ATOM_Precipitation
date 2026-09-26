#!/bin/bash
# (b) + (c) on the alf_tie stack (2026-09-26). 600 from scratch, cli/atm_gpa, 24 threads, gated on
# run_verify_dgpa.sh (itself behind qth_on). Reference: alf_tie (same stack with BASE_SAT=1, unweighted g_p).
#   gpa_bc  ATM_MC_SGZ=1 ATM_MC_ENTR=1.0e-4 ATM_MC_QVD=2 ATM_MC_ALF1=5.44e-4 ATM_MC_BASE_SAT=2 ATM_MC_GP_AREA=1
# PRE-REGISTERED (alf_tie at 600: 2956 mm/a, land 7620, r 0.209, sigma 9.81, g_p 10007, P_conv 1947; sgzb_ctl 987,
# r 0.471, sigma 2.37): (1) in-updraft q_c_u at 28N land falls from 6.9 g/kg toward the unseeded 0.13-1 g/kg (b);
# (2) g_p falls by well over an order of magnitude (c: a_u = 0.03 alone is 33x), so P_conv drops from 1947 to
# O(10-100) mm/a and total precip returns near sgzb_ctl's ~990; land back toward ~850 (NASA 782). (3) the open
# question: does a SMALL but non-zero convective rain (e_p now Tiedtke-sized) improve r / sigma / 0-15 vs sgzb_ctl,
# or is it another null? (4) exit 0, zero NaN through 155/357/483; MC_t truncation; max T location.
set -u; cd "$(dirname "$0")"; rm -f GPA600_DONE
until [ -f DGPA_VERIFY_DONE ]; do sleep 30; done
BAD=$(sed -n '1,/A  OFF BRANCH/p' run_verify_dgpa.out | grep 'DIFFERS:' | grep -vc 'RUN_CONFIG.txt')
if [ "$BAD" != 0 ] || ! grep -q "A  OFF BRANCH" run_verify_dgpa.out || [ "$(grep -c 'CONTROL.*PASS' run_verify_dgpa.out)" != 2 ]; then
    echo "byte check did not pass -- NOT started"; cat run_verify_dgpa.out; touch GPA600_DONE; exit 1; fi
mkdir output_gpa_bc || { touch GPA600_DONE; exit 1; }
echo "gpa_bc start $(date +%H:%M)"
env OMP_NUM_THREADS=24 ATM_CWB_DIAG=1 ATM_MC_CAP_DIAG=1 ATM_MC_SGZ=1 ATM_MC_ENTR=1.0e-4 ATM_MC_QVD=2 ATM_MC_ALF1=5.44e-4 \
    ATM_MC_BASE_SAT=2 ATM_MC_GP_AREA=1 ../cli/atm_gpa config_gpa_bc.xml > gpa_bc.log 2>&1
echo "gpa_bc exit $?  NaN $(grep -c 'NaN/Inf DETECTED' gpa_bc.log)  $(date +%H:%M)"
touch GPA600_DONE
