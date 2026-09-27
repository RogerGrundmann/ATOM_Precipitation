#!/bin/bash
# -O2 VALIDATION (2026-09-27, added at the user's instruction: run it only if run_scaling.sh shows -O2 >= 2x faster).
# 600 from scratch, cli/atm_O2 (eb2181c source, -O2 per its DW_AT_producer), DEFAULT configuration, the exact
# config of sgzb_ctl (the reference run for the current default: cli/atm_wcr -O0, 6 threads) with only the output
# path changed, and the same print-only diagnostics (ATM_CWB_DIAG=1 ATM_MC_CAP_DIAG=1).
# Between atm_wcr and eb2181c the default branch is byte-identical by the chain of byte checks (closure revert is
# already in sgzb_ctl's banner; later commits are default-off knobs and a console print).
# PRE-REGISTERED against sgzb_ctl at iteration 600 (987 mm/a, r 0.471, sigma 2.37, bands 3388/211/205/24,
# land/ocean 854/1039, drift -0.01/iter, converged 1): every scored metric agrees to the fixed-thread-count
# non-determinism (~0.1-0.5 %); exit 0, zero NaN through 155/357/483; max|u|/|v|/|w| at the same cells.
# -O2 -ffast-math reorders floating point, so byte identity is NOT expected and not tested. A difference beyond
# the parity noise in any band or in r/sigma means -O2 is NOT validated and the Makefile stays at -O0.
set -u; cd "$(dirname "$0")"; rm -f O2VAL600_DONE
mkdir output_o2val || { touch O2VAL600_DONE; exit 1; }
echo "o2val start $(date +%H:%M)"
env OMP_NUM_THREADS=${NT:-24} ATM_CWB_DIAG=1 ATM_MC_CAP_DIAG=1 ../cli/atm_O2 config_o2val.xml > o2val.log 2>&1
echo "o2val exit $?  NaN $(grep -c 'NaN/Inf DETECTED' o2val.log)  $(date +%H:%M)"
touch O2VAL600_DONE
