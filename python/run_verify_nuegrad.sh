#!/bin/bash
# BOTH-DIRECTIONS BYTE CHECK for the 2026-09-21 correctness flip of TWO defaults together:
#   ATM_NUE_GRAD        0.0 -> 1.0   (the grad(nu).grad(phi) half of div(nu grad phi))
#   ATM_BC_SECOND_ORDER   0 -> 1     (the 2nd-order one-sided Neumann the comments describe)
#
# ⚠ THIS IS A CORRECTNESS FLIP, NOT A MEASURED ONE, AND THE TWO KNOBS ARE NOT IN THE SAME
# EVIDENTIAL STATE. `ATM_NUE_GRAD` was measured connected (2026-09-07: every field written
# after iteration 0 differs, velocity 5e-06..1.3e-05 rms, scalars 2.9e-04 / 8.1e-04) and stable
# at full strength. `ATM_BC_SECOND_ORDER` has NEVER been run in this model -- its own header
# says "UNMEASURED HERE" and warns that the ocean's null must not be assumed to carry over,
# because this model has a moist boundary layer and IceSchemeCommon uses c43/c13 on the
# PRECIPITATION FLUX edges, a field with sharp gradients. 29 call sites move at once.
# So this check clears BYTE-IDENTITY OFF and CONNECTIVITY ON; it does NOT clear the climate,
# and the 600-iteration from-scratch arm through 155/357/483 is the thing that does.
#
# FOUR ARMS, so BOTH directions are checked rather than one (both knobs already existed in the
# pre-flip binary at default 0, which the 3-arm form used for brand-new knobs cannot exploit):
#   ng_old  cli/atm_preng (pre-flip HEAD, md5 1b775b8e), clean      = the SHIPPED branch
#   ng_off  new binary, both explicitly 0      -- MUST equal ng_old   (restore direction)
#   ng_new  new binary, clean                  = the NEW default
#   ng_on   cli/atm_preng, both explicitly on  -- MUST equal ng_new   (flip direction)
#   CONTROL ng_new vs ng_off MUST DIFFER, or both passes above are vacuous.
#
# Safe to run while the PTOP sweep is executing: that sweep runs `cli/atm_mn`, a different file,
# so rebuilding `cli/atm` cannot pull the rug from under it. (run_verify_latent.sh had to wait
# for exactly that reason -- its arms were executing the binary being rebuilt.)
set -u; cd "$(dirname "$0")"
rm -f NUEGRAD_DONE
. ./verify_lib.sh
rm -rf output_ng_new output_ng_old output_ng_off output_ng_on
mkdir -p output_ng_new output_ng_old output_ng_off output_ng_on
arm ng_old ../cli/atm_preng config_ng_old.xml
arm ng_off ../cli/atm       config_ng_off.xml ATM_NUE_GRAD=0.0 ATM_BC_SECOND_ORDER=0
arm ng_new ../cli/atm       config_ng_new.xml
arm ng_on  ../cli/atm_preng config_ng_on.xml  ATM_NUE_GRAD=1.0 ATM_BC_SECOND_ORDER=1
wait_arms
cmp_dirs    ng_off ng_old "A  RESTORE direction (new binary, both=0  vs  pre-flip clean)"
cmp_dirs    ng_on  ng_new "B  FLIP direction    (pre-flip, both set  vs  new binary clean)"
want_differ ng_new ng_off "C  CONTROL new default vs off -- MUST differ"
touch NUEGRAD_DONE
