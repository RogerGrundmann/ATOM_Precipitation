#!/bin/bash
# 2026-09-28: the runs the repaired-but-unrun bug-list items are waiting for. ALL -O2, one pinned binary pair from
# c1aae23 (cli/atm_rp, cli/hyd_rp); every arm has its own control on the SAME binary (the old -O0 controls are not
# comparable). 8 runs x 3 threads concurrent.
# ATMOSPHERE, 600 from scratch (config_o2val derivative = current defaults):
#   rp_ctl   default                          -- control for rp_oqm
#   rp_oqm   ATM_OROG_Q_MASS=1 (B+4)
#   rp_uw    ATM_PRECIP_UPWIND=1              -- control for the snow window on the CONVERGED rain iteration
#   rp_uw2   UPWIND + ATM_SNOW_WINDOW=2 (cold side only)
#   rp_uw3   UPWIND + ATM_SNOW_WINDOW=3 (cold side + melt)
# OCEAN, 1000 -> 1600 restart from oc_vis (metric branch, where B.6 lives), as run_ohs600.sh:
#   rp_oc    metric branch control
#   rp_od240   HYD_DEEP_DRAG tau = 2.315e-4 d (20 s  = e-folding ~240 iterations; numerical, not physical)
#   rp_od2400  HYD_DEEP_DRAG tau = 2.315e-3 d (200 s = ~2400 iterations)
# PRE-REGISTERED:
#  1. rp_oqm: CWB orographic_shapiro water bucket -> ~0 (conservation); climate change small (default branch, filter
#     leak is ~1e4 mm/a only with the closure on). All scores quoted against rp_ctl at 600.
#  2. rp_uw vs 09-27 upw restarts: converged rain stays ~10 % of NASA through 600 from scratch (not a restart transient).
#  3. rp_uw2/rp_uw3: with a converged iteration the snow window ADDS rain (09-27 restart: 92 -> 177 mm/a with =3);
#     bit 1 (melt) no longer "collapses" rain. Bands judged, not the mean.
#  4. rp_od*: the B.6 lower-column rise (i=1 / i=12 rms speed; 1.268 in bn_ctl at -O0) falls monotonically with
#     stronger drag; exit 0, no runaway. If it does not move, the missing deep outlet is not the cause of B.6.
#  5. every run exit 0, zero NaN.
set -u; cd "$(dirname "$0")"; rm -f RP0928_DONE
A=(ctl oqm uw uw2 uw3); O=(oc od240 od2400)
for t in "${A[@]}" "${O[@]}"; do mkdir output_rp_$t || { echo "output_rp_$t exists"; touch RP0928_DONE; exit 1; }; done
for t in "${A[@]}"; do sed "s#output_o2val/#output_rp_$t/#" config_o2val.xml > config_rp_$t.xml; done
for t in "${O[@]}"; do sed "s#output_ohs_1/#output_rp_$t/#" config_ohs_1.xml > config_rp_$t.xml
  cp output_oc_vis/hyd_restart_0Ma_1000.bin output_rp_$t/; cp output_twctl/0Ma_smooth_Transfer_Atm_600.vwtp output_rp_$t/; done
echo "start $(date +%H:%M)  atm $(md5sum < ../cli/atm_rp | cut -c1-8)  hyd $(md5sum < ../cli/hyd_rp | cut -c1-8)"
run(){ ( env OMP_NUM_THREADS=${NT:-3} $3 ../cli/$2 config_rp_$1.xml > rp_$1.log 2>&1
         echo "rp_$1 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' rp_$1.log)  $(date +%H:%M)" ) & }
D="ATM_CWB_DIAG=1 ATM_CWB_BANDS=1 ATM_SNOW_DIAG=1 ATM_SR_DIAG=1 ATM_RAIN_PASS_DIAG=1"
VIS="HYD_METRIC_RADIUS=6370 HYD_RUN_NEUMANN=1 HYD_A_H_BIHARM=3.0e18"
run ctl  atm_rp "$D"
run oqm  atm_rp "$D ATM_OROG_Q_MASS=1"
run uw   atm_rp "$D ATM_PRECIP_UPWIND=1"
run uw2  atm_rp "$D ATM_PRECIP_UPWIND=1 ATM_SNOW_WINDOW=2"
run uw3  atm_rp "$D ATM_PRECIP_UPWIND=1 ATM_SNOW_WINDOW=3"
run oc     hyd_rp "$VIS"
run od240  hyd_rp "$VIS HYD_DEEP_DRAG=2.315e-4"
run od2400 hyd_rp "$VIS HYD_DEEP_DRAG=2.315e-3"
wait
sed -e "s#runs = \[.*#runs = [('oc','output_rp_oc'),('od240','output_rp_od240'),('od2400','output_rp_od2400'),('bn_ctl_O0','output_bn_ctl')]#" ocprofile_ohs.py > ocprofile_rp.py
python3 ocprofile_rp.py > rp0928_profile.txt 2>&1
for t in "${A[@]}"; do echo "== rp_$t"; grep "by |latitude|" rp_$t.log | tail -1; grep "model .*NASA .*bias" rp_$t.log | tail -1 | cut -c1-150
  grep "orographic" rp_$t.log | tail -1; grep "SR DIAG\] as shares" rp_$t.log | tail -1; done
touch RP0928_DONE
