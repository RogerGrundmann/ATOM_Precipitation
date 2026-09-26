#!/bin/bash
# The (b)/(c) split (2026-09-26), queued behind run_gpa600.sh (GPA600_DONE). 600 from scratch each, cli/atm_gpa,
# 24 threads, serial. Same alf_tie stack (SGZ + ENTR 1e-4 + QVD=2 + ALF1 5.44e-4); references alf_tie and gpa_bc.
#   gpa_b  ATM_MC_BASE_SAT=2                      -- (b) alone: fraction-weighted base seed, g_p unweighted
#   gpa_c  ATM_MC_BASE_SAT=1 ATM_MC_GP_AREA=1     -- (c) alone: full saturation seed, g_p area-weighted
# PRE-REGISTERED: (b) alone lowers q_c_u (28N land 6.9 g/kg -> toward ~1) but g_p stays unweighted, so convective
# rain remains well above gpa_bc; (c) alone cuts g_p by a_u = 0.03 (33x) from alf_tie's 1e4 mm/a but keeps the
# 5.3 g/kg seeded water, which then detrains / is carried as condensate instead of raining. Which one removes the
# land excess (alf_tie 7620 mm/a, NASA 782) is the question. Stability as for gpa_bc.
set -u; cd "$(dirname "$0")"; rm -f GPASPLIT_DONE
until [ -f GPA600_DONE ]; do sleep 60; done
if ! grep -q "gpa_bc exit" run_gpa600.out; then echo "gpa_bc did not run -- split NOT started"; cat run_gpa600.out; touch GPASPLIT_DONE; exit 1; fi
STK="ATM_MC_SGZ=1 ATM_MC_ENTR=1.0e-4 ATM_MC_QVD=2 ATM_MC_ALF1=5.44e-4"
for t in gpa_b gpa_c; do
  mkdir output_$t || { touch GPASPLIT_DONE; exit 1; }
  if [ $t = gpa_b ]; then X="ATM_MC_BASE_SAT=2"; else X="ATM_MC_BASE_SAT=1 ATM_MC_GP_AREA=1"; fi
  echo "$t start $(date +%H:%M)  $X"
  env OMP_NUM_THREADS=24 ATM_CWB_DIAG=1 ATM_MC_CAP_DIAG=1 $STK $X ../cli/atm_gpa config_$t.xml > $t.log 2>&1
  echo "$t exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
done
touch GPASPLIT_DONE
