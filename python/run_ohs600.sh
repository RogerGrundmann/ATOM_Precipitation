#!/bin/bash
# HYD_HYDRO_SPLIT -- the ocean analogue of ATM_HYDRO_SPLIT (2026-09-26). 1000 -> 1600 RESTARTS from the ladder
# seeds, exactly as run_bn600.sh, binary cli/hyd_hsp, 2 arms x 12 threads concurrent, gated on run_verify_vohs.sh.
#   ohs_1    metric branch (METRIC_RADIUS=6370 RUN_NEUMANN=1 A_H_BIHARM=3e18) + HYD_HYDRO_SPLIT=1.0   seed oc_vis
#           -- against bn_ctl (knob 0) and bn_1 (HYD_BUOY_CONSISTENT=1.0, the radial runaway)
#   ohs_sh1  SHIPPED metric + HYD_HYDRO_SPLIT=1.0                                                seed oc_ctl
#           -- against bn_sh1 (shipped + BUOY_CONSISTENT=1.0: max|u| 77.7 m/s, the 07-14 blow-up)
# PRE-REGISTERED:
#   1. the radial runaway is GONE: hs_1 rms|u_radial| within ~10x of bn_ctl (bn_1 was 7.4e-3 m/s, 3000x);
#      hs_sh1 does not run away (bn_sh1 77.7 m/s). The buoyancy is not in rhs_u by construction.
#   2. the force reaches the horizontal flow: [HYDRO_SPLIT] max|p_hb| printed; v/w differ from bn_ctl. The
#      geostrophic ADJUSTMENT itself is out of reach (1/f = 324 000 iterations) -- do not score balance.
#   3. B.6 profile (i = 1 / i = 12 on the fixed column set, ocprofile_hs.py): bn_ctl 1.268. Open whether
#      a baroclinic pressure gradient changes the lower-column rise; no expectation registered.
#   4. exit 0, zero NaN; mean T and KE, no runaway.
set -u; cd "$(dirname "$0")"; rm -f OHS600_DONE
until [ -f VOHS_VERIFY_DONE ]; do sleep 60; done
if ! grep -q "OFF-BRANCH.*PASS" run_verify_vohs.out || ! grep -q "CONTROL.*PASS" run_verify_vohs.out; then
    echo "byte check did not pass -- NOT started"; cat run_verify_vohs.out; touch OHS600_DONE; exit 1; fi
mkdir output_ohs_1 output_ohs_sh1 || { touch OHS600_DONE; exit 1; }
for d in ohs_1:oc_vis ohs_sh1:oc_ctl; do t=${d%%:*}; s=${d##*:}
  cp output_$s/hyd_restart_0Ma_1000.bin output_$t/; cp output_twctl/0Ma_smooth_Transfer_Atm_600.vwtp output_$t/; done
VIS="HYD_METRIC_RADIUS=6370 HYD_RUN_NEUMANN=1 HYD_A_H_BIHARM=3.0e18"
echo "start $(date +%H:%M)"
run(){ ( env OMP_NUM_THREADS=12 $2 ../cli/hyd_hsp config_$1.xml > $1.log 2>&1
         echo "$1 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $1.log)  $(date +%H:%M)" ) & }
run ohs_1   "$VIS HYD_HYDRO_SPLIT=1.0"
run ohs_sh1 "HYD_HYDRO_SPLIT=1.0"
wait
python3 ocprofile_ohs.py > ohs600_profile.txt 2>&1
touch OHS600_DONE
