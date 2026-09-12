#!/bin/bash
# DOES max|w_u| = 40.6 m/s STAY PUT PAST ITERATION 600?
#
# EXACT CONTINUATION of output_vw0pdc: restart from its own atm_restart_0Ma_600.bin
# (md5 3d875ea1), nm=1200 -- in the ATMOSPHERE nm is the TOTAL iteration count, so this runs
# iterations 601-1200. 62 min instead of the 124 a from-scratch 1200 would cost, and a restart is
# the right instrument here because the question is about the trajectory BEYOND 600, not about
# the initial condition.
#
# THE FOUR KNOBS ARE SET EXPLICITLY even though all four are compiled-in defaults as of today.
# Deliberate: it makes this run a pure continuation of vw0pdc regardless of whether the default
# flips are correct, and the default flips are being certified separately by run_verify.sh. The
# banner will print them WITHOUT the `*`.
#   (A restart could not exercise CELLS_FROM_PSI or TROPO_INDEX_FIX anyway -- load_state()
#    overwrites everything the init path produces -- but the checkpoint was BUILT with both.)
#
# NEW HERE: panorama .vts ON (flag was false in every arm to date) and both the VTK slices and
# the panorama on a 100-ITERATION cadence -- ATM_VTK_STRIDE=5 against checkpoint=20, and
# panorama_print=100. vw0pdc wrote 93 slices / 6.1 GB at stride 1; this writes ~6 sets plus 6
# panoramas. Every CSV diagnostic stays on checkpoint=20 and is unaffected.
#
# WHAT vw0pdc ESTABLISHED, for reference at 600: max|w_u| -29.8, -36.2, -35.4, -41.1, -40.8,
# -40.6 -- rose to ~41 by iteration ~400 then flat-to-falling, at 41 % of the +-100 clamp.
# mean KE 36.434 flat to 0.13 % over 575 iterations, converged=1 from ~225.
#
# EXPECTED, recorded before the run:
#   - if max|w_u| stays in the high 30s / low 40s, the unfiltered zonal wind has a real
#     equilibrium and the convergence at 600 was not a plateau on the way up;
#   - if it resumes climbing, the filter was holding back a slow runaway and `converged`=1 was
#     premature -- which is the outcome that would matter, and it is exactly the shape the
#     precipitation showed when it "was stable to 600" and then climbed 18 % to 1200;
#   - watch mean KE and `converged` with it: a rising w_u with flat KE would be a redistribution.
# ALSO WATCH: max|u| (the radial runaway, 1.2965 at 600 and set by BUOY_CONSISTENT, not by _VW).
set -u
cd "$(dirname "$0")"
export OMP_NUM_THREADS=24
echo "=== [vw1200] _VW=0, 600->1200, panorama on, VTK+VTS every 100  $(date +%H:%M:%S) ==="
env ATM_CELLS_FROM_PSI=1 ATM_TROPO_INDEX_FIX=1 ATM_RADIAL_SHAPIRO_STRENGTH_VW=0 \
    ATM_V_MASSBAL_STRIDE=1 ATM_VTK_STRIDE=5 ATM_PDYN_CEILING=50.0 ATM_UBUD_BALANCE=1 \
    ../cli/atm config_vw1200.xml > vw1200.log 2>&1
echo "  exit $?  nan $(grep -ci nan vw1200.log)  $(date +%H:%M:%S)"
