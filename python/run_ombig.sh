#!/bin/bash
# 2026-09-30: OCN-METRIC -- om_scl with a LARGER equatorial biharmonic. om_scl (floor 0.26, B = 3e18, sin^4-scaled) was
# clean with om_met40's profile / radial u / KE drift, but grid noise 0.797 against om_met40's 0.512 and om_ship's 0.549:
# the sin^4 scaling weakens the damping poleward. Explicit nabla^4 limit with scaling ~2e20, so 1e19 and 3e19 are safe.
#   om_b1e19 / om_b3e19  HYD_METRIC_RADIUS=6370 HYD_RUN_NEUMANN=1 HYD_A_H_BIHARM=<B> HYD_A_H_BIHARM_SCALED=1, floor 0.26
# Same binary as om_scl (cli/hyd_omsc), 1000 from scratch, 6 threads each, both arms concurrent (atm rhst_115 also running).
# PRE-REGISTERED: noise falls with B toward <= 0.55; profile, radial u and KE drift stay at om_scl's; stable through 1000.
set -u; cd "$(dirname "$0")"; rm -f OMBIG_DONE
for b in 1e19 3e19; do
  mkdir output_om_b$b || { touch OMBIG_DONE; exit 1; }
  [ -e config_om_b$b.xml ] && { touch OMBIG_DONE; exit 1; }
  sed "s#output_om_ship/#output_om_b$b/#" config_om_ship.xml > config_om_b$b.xml
  cp output_twctl/0Ma_smooth_Transfer_Atm_600.vwtp output_om_b$b/
done
echo "start $(date +%H:%M)  hyd_omsc $(md5sum < ../cli/hyd_omsc | cut -c1-8)"
for b in 1e19 3e19; do
  ( env OMP_NUM_THREADS=6 HYD_METRIC_RADIUS=6370 HYD_RUN_NEUMANN=1 HYD_A_H_BIHARM=${b/e/.0e} HYD_A_H_BIHARM_SCALED=1 \
      ../cli/hyd_omsc config_om_b$b.xml > om_b$b.log 2>&1
    echo "om_b$b exit $?  NaN $(grep -c 'NaN/Inf DETECTED' om_b$b.log)  $(date +%H:%M)" ) &
done
wait
python3 ocprofile.py ship:output_om_ship/hyd_restart_0Ma_1000.bin met40:output_om_met40/hyd_restart_0Ma_1000.bin \
    scl:output_om_scl/hyd_restart_0Ma_1000.bin b1e19:output_om_b1e19/hyd_restart_0Ma_1000.bin \
    b3e19:output_om_b3e19/hyd_restart_0Ma_1000.bin > ombig_profile.txt 2>&1
cat ombig_profile.txt
for t in ship met40 scl b1e19 b3e19; do echo "== om_$t conv: $(tail -1 output_om_$t/convergence_hyd.csv)"; grep 'max u-component' om_$t.log | tail -1 | cut -c1-85; done
touch OMBIG_DONE
