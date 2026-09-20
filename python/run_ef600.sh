#!/bin/bash
# B+3 -- ATM_EVAP_FLUX, 600 from scratch. CONCURRENT, 12 THREADS EACH.
# (2026-09-20). See the header note at the bottom on why this is the slower arrangement; it is
# so sequential-24 is 160 min against 100 min for concurrent-12. DO NOT RE-SERIALISE.
#
# THE DEFECT. waterVapourEvaporation computes Dalton/Meyer/Rohwer into m.Evaporation.y and never
# uses it -- print, VTK, and nothing that drives the model. What moistens the air is c_eq, the
# humidity at which evaporation would BALANCE precipitation, as a relaxation at i=0 and an
# ABSOLUTE ADDITION at i=1..3. A humidity PRESCRIPTION wearing a flux's name.
# ColumnWaterBudget charges the stage +5.4459e+06 mm/a against the model's own E of 500.4 --
# ~10 900x -- and 30.7 % of ocean cells sit ON the c_sat_i cap.
# Moisture twin of omega_teq and worse: w_norm = 0.6439/iteration is an e-folding of 0.97
# iterations = 0.194 s, against a bulk flux needing ~1e5 iterations. omega_teq is 3290x.
#
# LEVEL 0 IS NOT TOUCHED: the 2026-09-06 decision -- level 0 is a PRESCRIBED SKIN, a surface flux
# belongs at the first AIR level, and level 0 is not RK4-integrated so removing its relaxation
# would leave it unconstrained. Same choice as ATM_SFC_FLUX and HYD_SFC_FLUX.
#
# ALREADY CONFIRMED on the discarded 10-thread partial (output_*_partial10t), at iteration 44:
# injection 1.282e-01 -> 1.55e-06 mm/call (83 000x) and the c_sat cap 31.06 % -> 0.00 % of ocean
# cells. So the knob engages exactly as designed; what these arms add is the climate response.
# Early matched-iteration precipitation was -1 to -5 %, NOT the collapse pre-registered below.
#
# PRE-REGISTERED (kept verbatim from the 10-thread launch so the record is honest):
#  1. A LARGE DRYING -- shipped increment ~2.4e-03 kg/kg per iteration against the flux's ~1.6e-08.
#  2. PRECIPITATION FALLS a lot; how much says how much of this model's rain rests on the
#     injection. ⚠ The partial already suggests this is WRONG in magnitude.
#  3. THE c_sat CAP STOPS BINDING -- CONFIRMED, 31.06 % -> 0.00 %.
#  4. Watch 155, 357, 483, and watch Precip: exit 0 with zero NaN is not a stability criterion.
#
# ON THE THREAD ARRANGEMENT, measured today and recorded so it is not re-litigated:
#   24 thr 6.2 s/iter | 12 thr 10.0 | 10 thr 11.4 | 8 thr 14.0
# Doubling 12 -> 24 buys 1.6x, not 2x, so a PAIR is faster run concurrently at 12 each (~100 min)
# than sequentially at 24 each (~124 min). The arms are independent processes and the marginal
# thread is worth less than the marginal process -- the same finding that parallelised the byte
# checks on 2026-09-14. This script does it sequentially at 24 because that is what was asked.
set -u
cd "$(dirname "$0")"
export OMP_NUM_THREADS=12
rm -f EF600_DONE
echo "=== ef600 CONCURRENT, 12 threads each $(date +%H:%M:%S)"
( env ATM_CWB_DIAG=1                 ../cli/atm_ef config_ef_ctl.xml > ef_ctl.log 2>&1
  echo "  ef_ctl exit $?  nan $(grep -ci nan ef_ctl.log)  $(date +%H:%M:%S)" ) &
( env ATM_CWB_DIAG=1 ATM_EVAP_FLUX=1 ../cli/atm_ef config_ef_on.xml  > ef_on.log  2>&1
  echo "  ef_on  exit $?  nan $(grep -ci nan ef_on.log)   $(date +%H:%M:%S)" ) &
wait
echo "=== ef600 done $(date +%H:%M:%S)"
touch EF600_DONE
