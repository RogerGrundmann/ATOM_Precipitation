#!/bin/bash
# RUN 2026-10-08 (user: "try to finish RC-RESID" -> "run step 1"). cli/atm_cmb (-O2, HEAD 5035541, md5 bdfbc89f) built today.
# ARMS CHANGED from the 10-07 proposal (0.85 / 0.75 / 0.65 alone): the evaporation gustiness is screened in the same arms, because E does not move
# the rain on a 120 s run:   a ATM_MC_CMB_OCEAN=0.85 + ATM_EVAP_GUST=5;   b 0.80 + 5;   c 0.85 + 6.   (working branch: 1 / 3)
# WHY THE GUST (read-only on wb68, same deficit, same rain): gust 3 leaves the formula wind at 3.5 m/s at 25-35 deg; 5 -> ocean E 1185, land E 376
# (E/P 0.58), global E 955, P/E 1.07 (1.03 at P = 980); 6 -> 1258 / 396 / 1014, P/E 1.01 (0.97). PRE-REGISTERED P/E: a 1.01-1.05, b 1.00-1.04, c 0.95-0.99.
# Rain ranges: a as the 0.85 row below; b between the 0.85 and 0.75 rows (global -2 to 0 %); c = a.
# 2026-10-07: wb69a / b / c = working branch (wb68 stack, 56 knobs) + ATM_MC_CMB_OCEAN (new knob: factor on the cloud-base mass-flux coefficient over
# OCEAN, applied before the ceiling ATM_MC_MB_SAT_OCEAN=0.05, which is unchanged):   a 0.85;   b 0.75;   c 0.65.
# SCREENING AT nm 60 from scratch, cli/atm_cmb (-O2, HEAD + the knob), 3 x 6 threads.
# WHY (wb68, |lat| <= 30 ocean): convective rain is a fixed function of the sea temperature (r 0.886; NASA's rain 0.598; iteration 20 against 520
# pattern r 0.9999); at 87E the flux is 93-100 % of the ceiling from 9S to 15N. Convective medians by SST class 26-27 / 27-28 / 28-29 / 29-31 C:
# 2.1 / 3.2 / 5.5 / 6.3 mm/d against NASA's total 1.7 / 2.6 / 4.3 / 6.0. The model's own surface moisture convergence explains none of what the SST does
# not (r -0.04 with the NASA residual), so a convergence closure is not the route; this lowers the medians below the warm pool and keeps the peak.
# Control = wb68 (600; the nm 60 screen wb67b agrees to 0.1 %): 1016.9 mm/a (+4.0 %), r .648, sigma 1.14, ocean 0-15 1796 (NASA 1440), tropical ocean
# mean 3.59 mm/d (3.05), rainless 0 % (7 %), p50/p90/p99 2.8/6.9/8.3 (2.5/6.4/8.4), ocean 15-35 810 (809; conv 406), N row 22-26 728 (conv 498; NASA 637).
# PRE-REGISTERED from the 87E flux-rain relation applied to the wb68 map (global | r | sigma | ocean 0-15 | medians 26-27 / 27-28 / 28-29 / 29-31 | p99):
#   a 0.85: +0.3 % (-1 to +1.5) | .643-.649 | 1.10-1.12 | 1650-1720 | 1.8 / 2.8 / 5.0 / 6.2 | 8.0-8.3
#   b 0.75: -2.4 % (-4 to -1)   | .640-.647 | 1.08-1.10 | 1560-1640 | 1.5 / 2.4 / 4.6 / 6.1 | 7.9-8.2
#   c 0.65: -5.5 % (-7.5 to -4) | .632-.642 | 1.05-1.08 | 1450-1540 | 1.3 / 2.1 / 4.2 / 6.0 | 7.8-8.1
# Ocean 15-35 falls with it (conv 406 is half of 810): a ~780, b ~755, c ~725. Land, the extratropics and the 65-90 band unchanged.
# USABLE if r >= .645, rainless tropical ocean <= 7 % (NASA), p99 >= 7.9, ocean 15-35 >= 750, global within 3 % of NASA.
# RESULT (2026-10-08 08:23, nm 60, 3 x 6 threads; wb69.out), every pre-registered range met:
#   a 0.85 / gust 5: 981.6 mm/a (+0.3 %), r .646, sigma 1.10, bands 1666/608/1015/320, ocean 0-15 / 15-35 1666 / 754, tropical-ocean mean 3.31 mm/d (3.05),
#                    p50/p90/p99 2.4/6.7/8.2, E 964.8 (ocean 1205, land 359), P/E 1.02  -> USABLE
#   b 0.80 / gust 5: 967.4 (-1.1 %), r .645, sigma 1.09, ocean 0-15 / 15-35 1616 / 734 (fails >= 750), P/E 1.00
#   c 0.85 / gust 6: rain = a; E 1022 (ocean 1280, land 374), P/E 0.96
# Land, the extratropics and the 65-90 band identical in all three. NOT adopted (the user decides); a 600 is owed before adoption.
set -u; cd "$(dirname "$0")"; rm -f WB69_DONE
for t in wb69a wb69b wb69c; do
  mkdir output_$t || { touch WB69_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB69_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>60<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>60<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_cmb $(md5sum < ../cli/atm_cmb | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_VTK_STRIDE=2"
env OMP_NUM_THREADS=6 $K ATM_MC_CMB_OCEAN=0.85 ATM_EVAP_GUST=5 ../cli/atm_cmb config_wb69a.xml > wb69a.log 2>&1 &
env OMP_NUM_THREADS=6 $K ATM_MC_CMB_OCEAN=0.80 ATM_EVAP_GUST=5 ../cli/atm_cmb config_wb69b.xml > wb69b.log 2>&1 &
env OMP_NUM_THREADS=6 $K ATM_MC_CMB_OCEAN=0.85 ATM_EVAP_GUST=6 ../cli/atm_cmb config_wb69c.xml > wb69c.log 2>&1 &
wait
for t in wb69a wb69b wb69c; do
  V=output_$t/0Ma_smooth_Atm_radial_0_60.vtk
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)  slices: $(ls output_$t | grep -c radial_0)"
  echo "banner diff vs wb68 (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb68.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "output\|nm=\|VTK_STRIDE" | head -8
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90
  grep -a 'water budget closure' $t.log | tail -2 | cut -c1-120
  python3 evapb.py $t 60 | sed -n '2p' | cut -c1-200
  python3 oceanb.py $V
  ITER=60 python3 pacband.py $t | grep -v "^ \+-\?[0-9]\+:\|ocean by latitude"
  python3 row36.py $t 60 | sed -n "3,6p" | cut -c1-170
  python3 landb.py $V | grep "E Austral\|Amazon\|Congo\|Madagascar\|India\|global max"
done
touch WB69_DONE
