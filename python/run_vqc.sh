#!/bin/bash
# 2026-10-03: ATM_RH_LAND_QCAP (default 0) -- -O0 byte check, 1 thread, nm 20 from scratch. Waits for the vgb check to finish.
# old = cli/gbo_atm (4960020 + ATM_MC_GATE_BLEND), new = cli/qco_atm (+ ATM_RH_LAND_QCAP).
#   A  new clean == old clean;  B  new knob=0 == old clean;  C  CONTROL new knob=1 != old clean
#   e  (information) the wb10a initial state with the cap: census line + does the Horn still convect at iteration 20
set -u; cd "$(dirname "$0")"; rm -f VQC_DONE
until [ -e VGB_DONE ]; do sleep 20; done
. ./verify_lib.sh
for t in a b c d e; do mkdir output_vqc_$t || { touch VQC_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vqc_$t/#" config_vconv_new.xml > config_vqc_$t.xml; done
arm vqc_a ../cli/qco_atm config_vqc_a.xml
arm vqc_b ../cli/qco_atm config_vqc_b.xml ATM_RH_LAND_QCAP=0
arm vqc_c ../cli/gbo_atm config_vqc_c.xml
arm vqc_d ../cli/qco_atm config_vqc_d.xml ATM_RH_LAND_QCAP=1
( . ./working_branch.env; arm vqc_e ../cli/qco_atm config_vqc_e.xml ATM_RH_OCEAN=0.82 ATM_MC_ML_PARCEL=1 ATM_MC_T_ADD=0.15 ATM_MC_Q_ADD=4.5e-4 ATM_MC_ML_LCL=1 ATM_RH_LAND=2 ATM_RH_LAND_QCAP=1; wait )
wait_arms
for t in a b c d; do sed -i 's#output_vqc_[a-e]/#OUT/#' output_vqc_$t/RUN_CONFIG.txt; done
cmp_dirs vqc_a vqc_c "A  new clean == old clean"
cmp_dirs vqc_b vqc_c "B  new knob=0 == old clean"
for t in a b; do diff <(sed -E "s/  RH_LAND_QCAP=[^ ]*//g" output_vqc_$t/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vqc_c/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$t: RUN_CONFIG differs only by the banner token" || echo "$t: RUN_CONFIG differs BEYOND the banner token"; done
want_differ vqc_d vqc_c "C  CONTROL new knob=1 vs old clean -- MUST differ"
grep -a "RH-LAND-QCAP" vqc_d.log vqc_e.log
touch VQC_DONE
