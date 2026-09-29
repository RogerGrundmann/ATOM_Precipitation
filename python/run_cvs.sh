#!/bin/bash
# 2026-09-29: CONV-DEAD -- the convection stack on the energy-conserving working branch, 600 from scratch.
# Working branch = closure (sync 2) + ATM_DAMP_Q_VERT=0 ATM_DAMP_Q_HORIZ=0 ATM_DAMP_T_VERT=0 + UPWIND + SNOW_WINDOW=2.
# Control = stm_2 (run_stm.sh, cli/atm_sdr, same setup, 6 threads). Binary cli/atm_ced (fresh -O2, = 55489e9 + ED_AREA;
# SURF_DRAG default 1.0 here vs 0.0 in stm_2 -- measured climate-null). 2 x 8 threads.
#   cvs_all  ATM_MC_QVD=2 ATM_MC_ED_AREA=1 ATM_MC_ALF1=5.44e-4 ATM_MC_BASE_SAT=2
#   cvs_ne   same without ATM_MC_ED_AREA (is the new knob needed once QVD=2 and ALF1 are in?)
# PRE-REGISTERED (from the 40-iter restarts ced_q2all): P_conv at the ground O(100) mm/a (0.15 in stm_2), total
# precipitation well above stm_2's 55 mm/a; stable (exit 0, zero NaN, max|v| ~2.5); watch land over-rain (09-26 qvd2
# stack generated ~46 000 mm/a per convecting land column with BASE_SAT=1 -- here BASE_SAT=2) and MC_t truncation.
set -u; cd "$(dirname "$0")"; rm -f CVS_DONE
for t in all ne; do mkdir output_cvs_$t || { touch CVS_DONE; exit 1; }
  [ -e config_cvs_$t.xml ] && { touch CVS_DONE; exit 1; }
  sed "s#output_o2val/#output_cvs_$t/#" config_o2val.xml > config_cvs_$t.xml; done
BR="ATM_MC_DIAG=1 ATM_MC_CAP_DIAG=1 ATM_CWB_DIAG=1 ATM_PRECIP_UPWIND=1 ATM_SNOW_WINDOW=2 ATM_WATER_CLOSURE=1 ATM_DAMP_Q_VERT=0 ATM_DAMP_Q_HORIZ=0 ATM_DAMP_T_VERT=0 ATM_MC_QVD=2 ATM_MC_ALF1=5.44e-4 ATM_MC_BASE_SAT=2"
echo "start $(date +%H:%M)  atm_ced $(md5sum < ../cli/atm_ced | cut -c1-8)"
run(){ ( env OMP_NUM_THREADS=8 $BR $2 ../cli/atm_ced config_cvs_$1.xml > cvs_$1.log 2>&1
         echo "cvs_$1 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' cvs_$1.log)  $(date +%H:%M)" ) & }
run all "ATM_MC_ED_AREA=1"
run ne  ""
wait
for t in all ne; do echo "== cvs_$t"; grep -a "by |latitude|" cvs_$t.log | tail -1; grep -a "model .*NASA .*bias" cvs_$t.log | tail -1 | cut -c1-150
  grep -a -i "land .*ocean" cvs_$t.log | tail -1 | cut -c1-100; grep -a 'P_conv mean' cvs_$t.log | tail -1
  grep -a '\[MC DIAG\] as shares' cvs_$t.log | tail -1 | cut -c1-200; grep -a 'water budget closure' cvs_$t.log | tail -1
  grep -a 'max v-component' cvs_$t.log | tail -1 | cut -c1-80; tail -1 output_cvs_$t/convergence.csv; done
touch CVS_DONE
