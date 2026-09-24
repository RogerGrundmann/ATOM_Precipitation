#!/bin/bash
# HYD_BUOY_CONSISTENT -- the B.6 arm. 1000 -> 1600 RESTARTS from the 2026-09-21 from-scratch ladder
# seeds (correct salinity), binary cli/hyd_buoy, 4 arms x 6 threads concurrent (~15 s/iter, ~2.5 h).
#   bn_ctl  metric branch (METRIC_RADIUS=6370 RUN_NEUMANN=1 A_H_BIHARM=3e18), knob 0   seed oc_vis
#   bn_1    metric branch + HYD_BUOY_CONSISTENT=1.0                                   seed oc_vis
#   bn_01   metric branch + HYD_BUOY_CONSISTENT=0.1  (strength fallback)              seed oc_vis
#   bn_sh1  SHIPPED metric + HYD_BUOY_CONSISTENT=1.0 -- re-test of the 2026-07-14 blow-up
#           (0.04 -> 19.5 m/s radial in 30 iters) on today's model                    seed oc_ctl
# NOT combined with HYD_PHYDRO_SALT/HYD_BAROCLINIC_PGF: with a consistent buoyancy the projection
# pressure p_dyn carries the hydrostatic anomaly and its horizontal gradient itself, so adding the
# p_hydro PGF would count the same force twice.
# QUESTION (B.6): does the metric branch's lower-column rise go away? i = 0 is NOT a seafloor: in
# shallow mode it is the 200 m TRUNCATION, and v/w there are CUBIC EXTRAPOLATION, so the i = 0 row
# is a boundary value, not a prognostic one. Scored on i = 1..12 ONLY. Ladder seeds at 1000:
#   vis 1.333 cm/s at -105 m (i=12) rising to 1.558 at -190 m (i=1), +17 %; ctl DECAYS 1.166 -> 0.951.
#   (The recorded "+52 %" used i = 0, 2.047 -- the extrapolated value.)
# PRE-REGISTERED:
#   1. stability: bn_1 max|u| stays O(1e-4..1e-3) m/s (metric branch sits at 1.4e-4); runaway
#      (>1e-2 and growing) = FAIL at 1.0, then read bn_01. bn_sh1 may blow up (07-14); the NaN
#      detector stops it, which is itself the answer.
#   2. profile (fixed full-depth column set, 1600, i = 1 / i = 12): the +17 % rise in bn_ctl shrinks
#      or reverses in bn_1. No change = the buoyancy is not what the lower column lacks.
#   3. radial velocity physical: rms|u| within ~10x of bn_ctl (68 m/yr), not the shipped 48 758 m/yr.
#   4. mean T and KE: no runaway; KE drift reported but not judged (1/f = 324 000 iterations).
set -u; cd "$(dirname "$0")"; rm -f BN600_DONE
until [ -f PR600_DONE ] && [ -f BUOYND_VERIFY_DONE ]; do sleep 60; done
if ! grep -q "OFF-BRANCH.*PASS" run_verify_buoynd.out; then
    echo "byte check did not pass -- arm NOT started"; cat run_verify_buoynd.out; touch BN600_DONE; exit 1; fi
VIS="HYD_METRIC_RADIUS=6370 HYD_RUN_NEUMANN=1 HYD_A_H_BIHARM=3.0e18"
echo "start $(date +%H:%M)"
run(){ ( env OMP_NUM_THREADS=6 $2 ../cli/hyd_buoy config_$1.xml > $1.log 2>&1
         echo "$1 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $1.log)  $(date +%H:%M)" ) & }
run bn_ctl "$VIS"
run bn_1   "$VIS HYD_BUOY_CONSISTENT=1.0"
run bn_01  "$VIS HYD_BUOY_CONSISTENT=0.1"
run bn_sh1 "HYD_BUOY_CONSISTENT=1.0"
wait
python3 ocprofile_levels.py > bn600_profile.txt 2>&1
touch BN600_DONE
