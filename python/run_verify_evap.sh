#!/bin/bash
# ATM_EVAP_FLUX byte check. The off-branch arm is compared against cli/atm_sf, which predates
# BOTH ATM_SATADJ_FREEZE_LATENT and ATM_EVAP_FLUX -- so this is a CHAIN check: passing it clears
# both off-branches at once, across two binary generations.
set -u; cd "$(dirname "$0")"; . ./verify_lib.sh
rm -rf output_ve_new output_ve_old output_ve_on; mkdir -p output_ve_new output_ve_old output_ve_on
arm ve_new ../cli/atm_ef config_ve_new.xml
arm ve_old ../cli/atm_sf config_ve_old.xml
arm ve_on  ../cli/atm_ef config_ve_on.xml  ATM_EVAP_FLUX=1
wait_arms
cmp_dirs    ve_new ve_old "A  off-branch CHAIN (new unset vs cli/atm_sf, 2 generations back)"
want_differ ve_on  ve_new "C  CONTROL EVAP_FLUX=1 vs off -- MUST differ"
touch VEVAP_DONE
