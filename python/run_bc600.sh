#!/bin/bash
# THE OWED 600-ITERATION FROM-SCRATCH ARM FOR ATM_BC_SECOND_ORDER -- and for ATM_NUE_GRAD with it.
#
# Both were flipped to default ON on 2026-09-21 on CORRECTNESS, with the both-directions byte
# check passing and NO climate measurement behind either. That is the ATM_METRIC_SIN_FLOOR
# precedent, and that knob got a 600-iteration clearance; these have not had one.
#
# ⚠ THE TWO ARE NOT IN THE SAME EVIDENTIAL STATE AND THIS ARM EXISTS FOR THE SECOND.
# `ATM_NUE_GRAD` was measured CONNECTED and stable at full strength on 2026-09-07 (velocity
# 5e-06..1.3e-05 rms, scalars 2.9e-04 / 8.1e-04 after 20 iterations).
# `ATM_BC_SECOND_ORDER` has NEVER been run past 20 iterations in this model. Its own header says
# "UNMEASURED HERE" and warns that the OCEAN's 300-iteration null must not be assumed to carry
# over: this model has 41 levels over a 16 km shell, a moist boundary layer, and
# `IceSchemeCommon` uses c43/c13 on the PRECIPITATION-FLUX edges -- a field with sharp gradients.
# 29 call sites move at once (BC_Atm.h 13, PressureSolverAtm.h 9, IceSchemeCommon.h 6).
#
# A TRIO RATHER THAN A PAIR, because it costs the same wall clock -- 3 arms x 8 threads fills the
# machine either way -- and it attributes the two flips SEPARATELY instead of as a bundle:
#   bc_ctl   NUE_GRAD=0.0 BC_SECOND_ORDER=0   the pre-flip branch
#   bc_bc    NUE_GRAD=0.0                     BC_SECOND_ORDER alone   (bc_ctl -> bc_bc)
#   bc_new   clean                            = the CURRENT DEFAULT    (bc_bc -> bc_new = NUE_GRAD)
# One pinned binary `cli/atm_bc` (md5 14a1e717). 600 from scratch, moist physics from iteration 0.
#
# AND bc_ctl IS ALSO A CHAIN CHECK. It should reproduce `output_mn_on` -- same configuration
# (MICRO_NDIM on, both correctness knobs off), different binary generation (cli/atm_mn, which
# predates SATADJ_FADE, FREEZE_LATENT, EVAP_FLUX and both flips) and without ATM_CWB_DIAG. If it
# does not reproduce it to ~three figures, something in four binary generations touched the
# default branch and the byte checks missed it.
#
# PRE-REGISTERED EXPECTATIONS:
#  1. BOTH ARE SMALL. Both forms of the Neumann condition are VALID -- this is an accuracy
#     change, not a wrong-physics one -- and the cross term completes a Laplacian rather than
#     adding a new force. Expect sub-percent on the precipitation and the bands.
#  2. THE RISK IS CONCENTRATED, NOT DIFFUSE. If ATM_BC_SECOND_ORDER bites, it bites through
#     IceSchemeCommon's precipitation-flux edges, i.e. at COASTS and over TERRAIN -- so read the
#     LAND/OCEAN split and the land value, not just the global mean. A big move there with a
#     small global mean is the signature to look for, and it would mean the boundary treatment is
#     load-bearing in the microphysics.
#  3. STABILITY through 155, 357 and 483, this tree's three documented failure points. And
#     ⚠ SCORE IT ON max|w|, max|u| AND mean KE -- NOT on `max w_u`, which is an UNCONSUMED
#     diagnostic measuring the updraft recurrence's 1/M_u amplification (2026-09-12), and which
#     this file leaned on for years of acceptances.
#  4. AND CHECK `Precip mean` AGAINST THE NASA LINE BEFORE CALLING ANYTHING AN IMPROVEMENT --
#     this tree's own standing rule, recorded after a cloud arm was written up as a success with
#     the precipitation never looked at, in a model named ATOM_Precipitation.
#  5. WHAT WOULD ARGUE FOR REVERTING: a material precipitation move, a band move, or any
#     instability. A NULL is the expected and acceptable outcome -- it is what clears a
#     correctness flip, exactly as the sin-floor's null did.
set -u; cd "$(dirname "$0")"
export OMP_NUM_THREADS=8
rm -f BC600_DONE
echo "=== bc600 0->600 from scratch, 3 arms concurrent, 8 threads each $(date +%H:%M:%S)"
( env ATM_NUE_GRAD=0.0 ATM_BC_SECOND_ORDER=0 ../cli/atm_bc config_bc_ctl.xml > bc_ctl.log 2>&1
  echo "  bc_ctl exit $?  nan $(grep -ci 'nan' bc_ctl.log)  $(date +%H:%M:%S)" ) &
( env ATM_NUE_GRAD=0.0                       ../cli/atm_bc config_bc_bc.xml  > bc_bc.log  2>&1
  echo "  bc_bc  exit $?  nan $(grep -ci 'nan' bc_bc.log)   $(date +%H:%M:%S)" ) &
( env                                        ../cli/atm_bc config_bc_new.xml > bc_new.log 2>&1
  echo "  bc_new exit $?  nan $(grep -ci 'nan' bc_new.log)  $(date +%H:%M:%S)" ) &
wait
echo "=== bc600 done $(date +%H:%M:%S)"
touch BC600_DONE
