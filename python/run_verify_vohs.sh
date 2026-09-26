#!/bin/bash
# OFF-BRANCH BYTE CHECK for HYD_HYDRO_SPLIT (default 0.0), 2026-09-26. 1 thread, nm = 20 from scratch, atm forcing
# seeded. old = cli/hyd_prehsp (HEAD 4556dc2), new = cli/hyd_hsp. Queued behind nothing -- launch by hand (postponed 2026-09-26 to the next day).
set -u; cd "$(dirname "$0")"; rm -f VOHS_VERIFY_DONE
# (no gate: launched by hand)
. ./verify_lib.sh
for t in old new ctl; do mkdir output_vohs_$t || { touch VOHS_VERIFY_DONE; exit 1; }
  cp output_twctl/0Ma_smooth_Transfer_Atm_600.vwtp output_vohs_$t/; done
arm vohs_old ../cli/hyd_prehsp config_vohs_old.xml
arm vohs_new ../cli/hyd_hsp    config_vohs_new.xml
arm vohs_ctl ../cli/hyd_hsp    config_vohs_ctl.xml HYD_HYDRO_SPLIT=1.0
wait_arms
cmp_dirs    vohs_new vohs_old "HYD A  OFF-BRANCH (new unset vs HEAD)"
want_differ vohs_ctl vohs_new "HYD C  CONTROL HYDRO_SPLIT=1 -- MUST differ"
grep -h "HYDRO_SPLIT\] \|WARNING\|NaN" vohs_ctl.log | head -5
touch VOHS_VERIFY_DONE
