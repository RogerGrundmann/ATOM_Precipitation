#!/bin/bash
# 2026-10-02: MC-GATE probe -- restarts of wb12a (RH_OCEAN 0.82) and wb12b (0.79), 600 -> 620, wb12 knobs, cli/atm_gate, 2 x 8 threads.
# PRE-REGISTERED: the coeff_recurr gate (|M_u| > 0.1, i.e. CAPE > ~460 J/kg at the base) is the cliff: ocean<30 columns with an open
# level ~20-40 % at 0.82 and a few % at 0.79, and the open columns carry ~all of the column g_p in both.
set -u; cd "$(dirname "$0")"; rm -f GATE_DONE
. ./working_branch.env
K="ATM_MC_DIAG=1 ATM_MC_ML_PARCEL=1 ATM_MC_T_ADD=0.15 ATM_MC_Q_ADD=4.5e-4 ATM_MC_ML_LCL=1 ATM_RH_LAND=2 ATM_MC_DEPTH_RAMP=1"
for t in a b; do
  mkdir output_gate_$t || { touch GATE_DONE; exit 1; }
  ln -s ../output_wb12$t/atm_restart_0Ma_600.bin output_gate_$t/atm_restart_0Ma_600.bin
  sed -e "s#output_wb12$t/#output_gate_$t/#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>-1<#" \
      -e "s#<restart_from_iter>-1<#<restart_from_iter>600<#" -e "s#<nm>600<#<nm>620<#" config_wb12$t.xml > config_gate_$t.xml
done
env OMP_NUM_THREADS=8 $K ATM_RH_OCEAN=0.82 ../cli/atm_gate config_gate_a.xml > gate_a.log 2>&1 &
env OMP_NUM_THREADS=8 $K ATM_RH_OCEAN=0.79 ../cli/atm_gate config_gate_b.xml > gate_b.log 2>&1 &
wait
for t in a b; do echo "== gate_$t  NaN $(grep -c 'NaN/Inf DETECTED' gate_$t.log)"; grep -a "\[MC-GATE\]" gate_$t.log | tail -2; done
touch GATE_DONE
