#!/bin/bash
# B+1 -- ATM_MICRO_NDIM, the 2783x microphysics coefficient, measured FROM SCRATCH TO 600.
#
# WHY THIS ARM EXISTS. The only measurement of this knob is 20 iterations (= 4.0 s of physical
# time) from `output_twctl`'s checkpoint, and that checkpoint was made 2026-09-03 -- BEFORE the
# four defaults flipped on 2026-09-12 (CELLS_FROM_PSI, TROPO_INDEX_FIX, _VW=0, V_MASSBAL_STRIDE)
# and before the 2026-09-14 ATM_BUOY_CONSISTENT revert. So it is a short arm on a configuration
# that is no longer the default, and its own write-up says "its sign over a long run is not what
# its sign over four seconds is". This pair is 600 from scratch on the CURRENT defaults.
#
# THE DEFECT, RE-DERIVED RATHER THAN QUOTED (2026-09-20, no run):
#   * The rates ARE per-second. `max_cloud_loss = cloud / dt_snow_dim` with
#     dt_snow_dim = step[i]/0.96 -- metres / (m/s) -- so the limiter that scales every cloud-side
#     rate divides a mixing ratio by a TIME; and c_c_au = 1.0e-3 is 1/s. Independent of the
#     "kg/(kg*s)" comment on the assignment, which in this tree is not evidence.
#   * The TIME UNIT IS metricShellLength()/u_0, pinned by the HORIZONTAL advection and not by
#     the comment: transport_c carries v_invrm * dcdthe, and inv_rm gives the 1/r factors the
#     planetary radius in units of L_atm (rmet ~ 397.5), so the term is
#     (v/u_0)/397.5 * dc/dthe = (a/(397.5*u_0)) * [v dc/(a dthe)], and a/397.5 = 16028 m
#     = metricShellLength(). So the scalar equation's tendency unit is L/u_0 with L = 16023 m,
#     which is what the knob's endpoint uses. (The RADIAL advection carries the separate
#     ATM_METRIC_EXACT 1/(rm+1)-vs-L/J defect and does not enter this.)
#   * AND THE TWO ENDS OF THE MICROPHYSICS DISAGREE, WHICH IS THE PART WORTH STATING.
#     The PRECIPITATION path is dimensionally RIGHT: dP_rain = S_x * mass_layer with
#     mass_layer = r_humid*step[i] [kg/m2], giving kg/(m2 s). The COLUMN path is not:
#     rhs_c gets S_x * r_humid [kg/m3] where it needs S_x * (L/u_0) [s]. At s = 1 the column
#     debit becomes EXACTLY the precipitation-flux integrand -- same quantity, same units -- so
#     the two can close against each other; on the shipped branch they differ by ~2800 and
#     cannot. That is the same shape as the max(0,.) floor that manufactured 8129 mm/a before
#     2026-09-01: water leaves the domain as rain while the reservoir it came from is not
#     debited for it.
#
# ⚠ THIS IS A MICROPHYSICS ARM, SO [[project-dynamics-cannot-move-precip]] DOES NOT APPLY.
# That rule -- 600 iterations is 120 s, a parcel moves 36 m, no dynamical change can move a
# precipitation band -- is about ADVECTIVE displacement. This knob acts on a LOCAL rate law in
# the column it is already in, and it already moved precipitation +0.5 % in twenty iterations.
# Precipitation is the right score here; for a dynamics knob it would not be.
#
# PRE-REGISTERED EXPECTATIONS, recorded because this file records them when they are wrong
# (run_btref.sh got its sign backwards five days ago):
#   1. PRECIPITATION RISES, and by MORE than the 0.5 % measured at 20 iterations. The net q-side
#      term is POSITIVE -- the schemes net-ADD to vapour+cloud+ice because S_ev dominates -- so
#      correcting the coefficient adds water 2783x harder. Anyone expecting a 2800x SINK to dry
#      the model out has the sign of the net term backwards.
#   2. PRECIPITABLE WATER BARELY MOVES. 1920 mm/a on a 30.7 mm reservoir is a 5.8-day e-folding
#      = 2.5e+06 iterations; 600 iterations is 2.4e-04 of one. So if precipitation moves, it is
#      a CONVERSION-EFFICIENCY change and not a moistening -- the same signature the post-600
#      precipitation drift has. If PW moves materially, expectation 2 is wrong and the e-folding
#      arithmetic needs re-deriving.
#   3. STABILITY IS THE REAL RISK AND IT IS NOT COVERED BY THE 20-ITERATION ARM. That arm was a
#      RESTART; this one runs the amplified source against initCloudIce's field from iteration 0,
#      through the three documented failure points 155, 357 and 483. Watch for NaN -- and watch
#      `Precip mean` too, because this tree's own record says exit 0 with zero NaN is NOT a
#      stability criterion (the dt ladder reached 8733 mm/a without ever NaN-ing).
#
# ATM_CWB_DIAG=1 on BOTH arms: it is the instrument that measured the defect, it is print-only
# and byte-identical off, and at checkpoint = 20 it reports 30 windows -- the trajectory of the
# water budget, not two endpoints. It also prints the applied microphysics row directly, so the
# knob's connectivity is re-confirmed in the arm itself rather than assumed from a 1-thread check.
set -u
cd "$(dirname "$0")"
export OMP_NUM_THREADS=12
rm -f MN600_DONE
echo "=== mn600 0->600 from scratch, both arms parallel, 12 threads each $(date +%H:%M:%S)"
( env ATM_CWB_DIAG=1                      ../cli/atm_mn config_mn_ctl.xml > mn_ctl.log 2>&1
  echo "  mn_ctl exit $?  nan $(grep -ci 'nan' mn_ctl.log)  $(date +%H:%M:%S)" ) &
( env ATM_CWB_DIAG=1 ATM_MICRO_NDIM=1.0   ../cli/atm_mn config_mn_on.xml  > mn_on.log  2>&1
  echo "  mn_on  exit $?  nan $(grep -ci 'nan' mn_on.log)   $(date +%H:%M:%S)" ) &
wait
echo "=== mn600 done $(date +%H:%M:%S)"
touch MN600_DONE
