#!/bin/bash
# OFF-BRANCH BYTE CHECK, 2026-09-26, for three default-off atmosphere knobs written together:
#   ATM_ONECAT_CLOUD_LIMIT (item 1), ATM_TURB_SIN_FLOOR (item 2), ATM_SEAM_Q_CONSERVE (item 6).
# 1 thread, nm = 20 from scratch. old = cli/atm_gpa (c3b78a7), new = cli/atm_fx4. LAUNCH BY HAND (postponed).
# (ColumnWaterBudget now counts the seam once -- print-only, ATM_CWB_DIAG is off here, so it cannot show.)
set -u; cd "$(dirname "$0")"; rm -f VFX4_VERIFY_DONE
. ./verify_lib.sh
for t in old new tsf sqc oc1 oc2; do mkdir output_vfx4_$t || { touch VFX4_VERIFY_DONE; exit 1; }; done
arm vfx4_old ../cli/atm_gpa config_vfx4_old.xml
arm vfx4_new ../cli/atm_fx4 config_vfx4_new.xml
arm vfx4_tsf ../cli/atm_fx4 config_vfx4_tsf.xml ATM_TURB_SIN_FLOOR=1
arm vfx4_sqc ../cli/atm_fx4 config_vfx4_sqc.xml ATM_SEAM_Q_CONSERVE=1
arm vfx4_oc1 ../cli/atm_fx4 config_vfx4_oc1.xml
arm vfx4_oc2 ../cli/atm_fx4 config_vfx4_oc2.xml ATM_ONECAT_CLOUD_LIMIT=1
wait_arms
cmp_dirs    vfx4_new vfx4_old "A  OFF BRANCH (new clean vs old)"
want_differ vfx4_tsf vfx4_new "C1 CONTROL TURB_SIN_FLOOR=1 -- MUST differ"
want_differ vfx4_sqc vfx4_new "C2 CONTROL SEAM_Q_CONSERVE=1 -- MUST differ"
want_differ vfx4_oc2 vfx4_oc1 "C3 CONTROL OneCat + ONECAT_CLOUD_LIMIT=1 -- MUST differ"
touch VFX4_VERIFY_DONE
