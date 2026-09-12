#!/bin/bash
# DOES A RELEASED p_dyn CEILING SURVIVE 600 ITERATIONS FROM SCRATCH?
#
# The 2026-09-12 2x2 showed the ceiling binds from iteration ~400 on 4.8-7.6 % of all cells
# whenever ATM_BUOY_CONSISTENT=1, and that releasing it is what lets the line solve's amplitude
# gain appear (cancellation 0.204 -> 0.291). Every one of those arms was 80 iterations = 16 s of
# physical time. The open question is stability: p_dyn_ceiling was lowered 10.0 -> 3.0 for the
# iteration-483 runaway (p_dyn pegged -10 over the Tian Shan/Pamir, the -10 -> 0 vertical gradient
# became a -16.8 pgr that drove the velocity), and NO arm so far goes anywhere near 483.
#
# ONE VARIABLE against output_tr600d: ATM_PDYN_CEILING 0 (shipped phase-dependent) -> 50.0.
# Same config (nm=600 from scratch, moist gate 0, checkpoint_save_iter 600), same knob set
#   ATM_CELLS_FROM_PSI=1 ATM_RADIAL_SHAPIRO_STRENGTH_VW=0.25 ATM_V_MASSBAL_STRIDE=1
#   ATM_CELL_ROT_DIAG=1 ATM_TROPO_INDEX_FIX=1
# tr600d IS the shipped-ceiling arm, so the pair costs 62 min and not two hours.
#
# ⚠ TWO HONEST CAVEATS, recorded before the run.
# (1) BINARY. tr600d used cli/atm_rot600, which has no ATM_UBUD_BALANCE, so this arm uses
#     cli/atm with ATM_PRESS_LINE_SOLVE UNSET. Off-branch identity between the two was
#     established yesterday at nm=4 / 1 thread, 6 of 7 written files byte-identical (the 7th is
#     RUN_CONFIG.txt, the output path and the new banner entry). That is a real check but a SHORT
#     one; ATM_UBUD_BALANCE is print-only, so the exposure is the line-solve code being inert
#     when unset, which is what the check tested.
# (2) THE OVERRIDE HITS BOTH PHASES. The shipped ceiling is 10.0 through the dry spin-up and 3.0
#     once total_iter_count > 300, and ATM_PDYN_CEILING=50 replaces both. That is a NULL during
#     spin-up and tr600d's own trajectory says so: max|p_dyn| is 1.077e-01 at iteration 320,
#     100x below the 10.0 it would have been clipped at. The override only acts after ~400.
#
# EXPECTED, recorded before the run so it cannot be rationalised afterwards:
#   - if it exits 0 with zero NaN through 155, 357 AND 483, a released ceiling is stable on this
#     configuration and the 80-iteration amplitude result can be read as real;
#   - if it NaNs, the 3.0 value is load-bearing and the amplitude repair needs the clamp RESIZED
#     with something else (a smooth saturation was already tried and FALSIFIED in 2026-06), not
#     simply raised -- and that is the answer, not a failure of the arm;
#   - max|p_dyn| is expected to exceed 3.0 from ~400 and keep growing; the number to watch is
#     whether it SATURATES or runs to the new 50.
# WATCH: max|u| (the 51x runaway), max w_u against the +-100 clamp, iterations 155 / 357 / 483,
# and whether cancellation keeps climbing past the 0.2908 reached at iteration 680 in pdc_ls.
set -u
cd "$(dirname "$0")"
export OMP_NUM_THREADS=24
echo "=== [pdc600] PDYN_CEILING=50, nm=600 from scratch  $(date +%H:%M:%S) ==="
env ATM_CELLS_FROM_PSI=1 ATM_RADIAL_SHAPIRO_STRENGTH_VW=0.25 ATM_V_MASSBAL_STRIDE=1 \
    ATM_CELL_ROT_DIAG=1 ATM_TROPO_INDEX_FIX=1 ATM_UBUD_BALANCE=1 ATM_PDYN_CEILING=50.0 \
    ../cli/atm config_pdc600.xml > pdc600.log 2>&1
echo "  exit $?  nan $(grep -ci nan pdc600.log)  $(date +%H:%M:%S)"
