#!/bin/bash
# 2026-10-04: wb29a / wb29b = working branch (wb27 stack) + ATM_MC_MB_SAT_OCEAN 0.045 / 0.03 kg/(m2 s) (user; RAIN-CONV: tropical-ocean rain too peaked).
# SCREENING: nm 220 from scratch, cli/atm_mbs (-O2, HEAD 755f191 + the knob), 2 x 8 threads. GATED on the -O0 byte check (run_vmso.sh) and on wb28.
# Control = wb25b / wb27 (same knob set, saturation off; cli/atm_g3 -- the knob's off branch is byte-identical at -O0).
# wb25b at 220: 912 mm/a, r .562, sigma 2.85, bands 3001/290/204/4.5, land/ocean 443/1098, P/E 1.56; ocean 0-15 3385 (NASA 1440), 15-35 406 (809);
# |lat|<=30 ocean: mean 5.40 mm/d (NASA 3.05), p50/p90/p99 2.1/18.3/28.6 (2.6/6.5/8.4); mean M_b 0.032 at 5.6 mm/d -> M_s 0.045 / 0.03 ~ 7.7 / 5.2 mm/d.
# PRE-REGISTERED (assumes rain ~linear in M_b; post-hoc tanh on the wb27 field): a -- ocean 0-15 1600-1950, ocean 15-35 330-400, ocean total 640-740,
# global 580-660, 0-15 band 1600-1900, sigma 1.5-1.9; b -- ocean 0-15 1200-1500, ocean 15-35 290-370, ocean total 520-620, global 500-580, 0-15 band
# 1300-1600, sigma 1.2-1.6. Both: land 443 +-2 %, Amazon / Congo unchanged, 35-65 and 65-90 unchanged, r .56-.62 (the flattened peak helps), P/E 0.9-1.15,
# wettest ocean cell < 12 mm/d (a) / < 9 (b), no NaN. The mean then sits 30-45 % below NASA: the next step is a higher ATM_RH_OCEAN.
# RISK: rain is not linear in M_b (detrainment, gate at 0.01) -- the reduction is smaller than the tanh estimate.
set -u; cd "$(dirname "$0")"; rm -f WB29_DONE
until [ -f VMSO_DONE ] && [ -f WB28_DONE ]; do sleep 20; done
# A, B, D: at most RUN_CONFIG.txt differs (the new banner token), and only by that token; C (control) must differ.
[ "$(grep -c 'differing [01]  ->' run_vmso.out)" = 3 ] && [ "$(grep -c 'differs only by the banner token' run_vmso.out)" = 3 ] && grep -q 'CONTROL.*-> PASS' run_vmso.out \
  || { echo "byte check not clean -- wb29 NOT started"; touch WB29_DONE; exit 1; }
for t in wb29a wb29b; do
  mkdir output_$t || { touch WB29_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB29_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>220<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>220<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_mbs $(md5sum < ../cli/atm_mbs | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1"
env OMP_NUM_THREADS=8 $K ATM_MC_MB_SAT_OCEAN=0.045 ../cli/atm_mbs config_wb29a.xml > wb29a.log 2>&1 &
env OMP_NUM_THREADS=8 $K ATM_MC_MB_SAT_OCEAN=0.03 ../cli/atm_mbs config_wb29b.xml > wb29b.log 2>&1 &
wait
for t in wb29a wb29b; do
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
  echo "banner diff vs wb25b (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb25b.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | head -8
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1
  grep -a "P_conv mean\|P_rain mean" $t.log | tail -2
  grep -a "\[MC-LO\] ocean<30\|\[MC-LO\] land<30" $t.log | tail -2 | cut -c1-330
  grep -a "\[MC-RG\] ocean<15" $t.log | tail -1 | cut -c1-260
  tail -1 output_$t/convergence.csv
  python3 oceanb.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk
  python3 landb.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk | grep -E "Amazon|Congo|global max"
done
touch WB29_DONE
