#!/bin/bash
# BOTH-DIRECTIONS BYTE CHECK for the 2026-09-21 flip of ATM_MICRO_NDIM 0.0 -> 1.0.
#
# Flipped on a MEASUREMENT, not an argument: 600 from scratch moved the two starved bands
# 35-65 deg +14.6 % and 65-90 deg +101 %, both saturating by iteration ~450, with the
# circulation identical to six figures -- and the PTOP re-sweep the same day showed the fitted
# constant it was coupled to does not have to move (35-65 varies 1.0 % across 425..500 and
# 65-90 is 20.3 in all four arms).
#
# ⚠ THE KNOB IS READ IN THREE PLACES AND ALL THREE DEFAULTS MOVED TOGETHER:
# RHS_Atm_Turb.cpp:1448 (the physics) and ColumnWaterBudget.h:281 and :371 (ATM_CWB_DIAG's
# "coefficient RK4 ACTUALLY uses" row and its printed strength). Flipping only the physics would
# leave the instrument reporting `r_humid` while RK4 applies L/u_0 -- an instrument lying about
# the one quantity it exists to measure, which is this tree's most frequent defect class.
# This check does NOT exercise that: ATM_CWB_DIAG is off in these arms. It is verified by
# reading the source, and by the 600-iteration arm whose applied row read 1712.05 against its
# own independent L/u_0 reference of 1712.0.
#
#   mnf_old  cli/atm_premn (pre-flip, md5 4af88963), clean    = the shipped branch
#   mnf_off  new binary, ATM_MICRO_NDIM=0.0  -- MUST equal mnf_old   (restore direction)
#   mnf_new  new binary, clean               = the new default
#   mnf_on   cli/atm_premn, ATM_MICRO_NDIM=1.0 -- MUST equal mnf_new (flip direction)
#   CONTROL  mnf_new vs mnf_off MUST DIFFER.
#
# Expect RUN_CONFIG.txt to differ in every comparison: the banner prints MICRO_NDIM=1.0* for the
# new compiled-in default, 0.0* for the old, and unstarred where the environment sets it. That
# is the documented exception; the physics files are what must match.
#
# Safe alongside the ocean ladder, which runs cli/hyd_oc.
set -u; cd "$(dirname "$0")"
rm -f MICRONDIM_DONE
. ./verify_lib.sh
rm -rf output_mnf_new output_mnf_old output_mnf_off output_mnf_on
mkdir -p output_mnf_new output_mnf_old output_mnf_off output_mnf_on
arm mnf_old ../cli/atm_premn config_mnf_old.xml
arm mnf_off ../cli/atm       config_mnf_off.xml ATM_MICRO_NDIM=0.0
arm mnf_new ../cli/atm       config_mnf_new.xml
arm mnf_on  ../cli/atm_premn config_mnf_on.xml  ATM_MICRO_NDIM=1.0
wait_arms
cmp_dirs    mnf_off mnf_old "A  RESTORE direction (new binary =0.0  vs  pre-flip clean)"
cmp_dirs    mnf_on  mnf_new "B  FLIP direction    (pre-flip =1.0    vs  new binary clean)"
want_differ mnf_new mnf_off "C  CONTROL new default vs off -- MUST differ"
touch MICRONDIM_DONE
