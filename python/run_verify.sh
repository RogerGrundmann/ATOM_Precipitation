#!/bin/bash
# BOTH-DIRECTIONS VERIFICATION OF THE FOUR DEFAULT FLIPS (2026-09-12).
#   ATM_CELLS_FROM_PSI            0   -> 1
#   ATM_TROPO_INDEX_FIX           0   -> 1
#   ATM_RADIAL_SHAPIRO_STRENGTH_VW (=STRENGTH) -> 0.0
#   ATM_V_MASSBAL_STRIDE          0   -> 1
# 1 THREAD, nm=20 FROM SCRATCH -- from scratch is mandatory: load_state() overwrites everything
# CELLS_FROM_PSI and TROPO_INDEX_FIX produce, so a restart cannot exercise either of them.
# 1 thread because this model is not bit-reproducible under OpenMP and a 24-thread pair could not
# tell a null flip from the threads.
#   A: new binary CLEAN     == old binary with the four SET      (the flip is what was measured)
#   B: new binary with four SET BACK == old binary CLEAN         (the revert restores the branch)
set -u
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
NEW=../cli/atm; OLD=../cli/atm_preflip
ON="ATM_CELLS_FROM_PSI=1 ATM_TROPO_INDEX_FIX=1 ATM_RADIAL_SHAPIRO_STRENGTH_VW=0 ATM_V_MASSBAL_STRIDE=1"
OFF="ATM_CELLS_FROM_PSI=0 ATM_TROPO_INDEX_FIX=0 ATM_RADIAL_SHAPIRO_STRENGTH_VW=1.0 ATM_V_MASSBAL_STRIDE=0 ATM_VTK_STRIDE=1"
echo "=== A: old+SET  $(date +%H:%M:%S)"; env $ON  $OLD config_va_old.xml > va_old.log 2>&1; echo "  exit $?"
echo "=== A: new+clean $(date +%H:%M:%S)"; env      $NEW config_va_new.xml > va_new.log 2>&1; echo "  exit $?"
echo "=== B: old+clean $(date +%H:%M:%S)"; env      $OLD config_vb_old.xml > vb_old.log 2>&1; echo "  exit $?"
echo "=== B: new+SETBACK $(date +%H:%M:%S)"; env $OFF $NEW config_vb_new.xml > vb_new.log 2>&1; echo "  exit $?"
for pair in "va_old va_new A(flip)" "vb_old vb_new B(revert)"; do
  set -- $pair; ok=0; bad=0
  for f in output_$1/*; do b="output_$2/$(basename "$f")"
    if cmp -s "$f" "$b"; then ok=$((ok+1)); else echo "  DIFFERS: $(basename "$f")"; bad=$((bad+1)); fi; done
  echo "$3: identical $ok / differing $bad"
done
echo "=== verify done $(date +%H:%M:%S)"
