#!/bin/bash
# *** PROPOSAL, NOT RUN YET (2026-10-08). ***
# wb73a / b / c = working branch (wb72 stack, 60 knobs) with the high-latitude ocean humidity re-set -- EXISTING knobs only (working branch:
# ATM_RH_STORM_SST=0.005, _MAX 0.04, _REF 4, ATM_RH_STORM_POLAR=1.08):
#   a  ATM_RH_STORM_SST=0.009, _MAX 0.10, _REF 3, ATM_RH_STORM_POLAR=1.0 (the polar floor off)
#   b  0.007, 0.08, 3, floor off (the milder dose)
#   c  0.009, 0.10, REF 2, floor off (drier warm water)
# SCREENING AT nm 60 from scratch, cli/atm_plm (-O2, = HEAD, md5 fe6987f6), 3 x 6 threads.
# WHY (STORM-SHAPE; socean.py, nrows.py on wb72): poleward of 50 deg the ocean rain still follows the sea temperature and NASA's does not. Ocean 50-66 deg
# by sea temperature, model / NASA, S and N: -6..-2 C 0.65 / 0.47, -2..1 C 0.85 / 0.95, 1..3 C 0.97 / 1.12, 3..5 C 1.03 / 1.36, 5..7 C 1.29 / 1.30,
# 7..9 C 1.40 / 1.21, 9..14 C 1.57 / 1.12. At the SAME latitude NASA is the same in both hemispheres and the model is not: 70-74 deg S 251 / 453 at -12.2 C,
# N 515 / 444 at -5.9 C. The polar floor lifts by latitude alone and puts the Arctic at 74-86N 1.2-1.5x over; the cap of 0.04 stops the cold-water lift at -4 C.
# OFFLINE ESTIMATE (stormfit.py on wb72, rain x exp(10 x relative RH change); an extrapolation -- the largest RH change tested so far is 0.05):
#   arm | global r | S rows 54-58 58-62 62-66 66-70 70-74 | N rows 54-58 58-62 62-66 66-70 70-74 74-82 | ocean 65-90 S / N (NASA 490 / 379) | rms log row misfit
#   now | .6496 | 1.06 1.12 0.95 0.58 0.55 | 1.19 1.30 1.12 0.90 1.16 1.24 | 316 / 428 | 0.234
#   a   | .651  | 0.95 1.09 1.09 0.87 0.75 | 0.98 1.12 1.03 0.86 1.02 1.09 | 435 / 385 | 0.124
#   b   | .652  | 0.98 1.08 1.00 0.73 0.62 | 1.05 1.17 1.05 0.84 0.95 0.95 | 368 / 354 | 0.153
#   c   | .651  | 0.89 0.99 1.00 0.83 0.74 | 0.91 1.02 0.94 0.80 0.96 1.06 | 418 / 367 | 0.129
# USABLE if the row misfit falls below 0.17, ocean 65-90 S >= 360 and N within 330-430, no ocean row between 46 and 74 deg above 1.25x NASA,
# global r >= .650, global mean within 2 % of NASA; land and everything equatorward of 48 deg unchanged (to 0.5 %).
set -u; cd "$(dirname "$0")"; rm -f WB73_DONE
for t in wb73a wb73b wb73c; do
  mkdir output_$t || { touch WB73_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB73_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>60<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>60<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_plm $(md5sum < ../cli/atm_plm | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_VTK_STRIDE=2 ATM_RH_STORM_POLAR=1.0"
env OMP_NUM_THREADS=6 $K ATM_RH_STORM_SST=0.009 ATM_RH_STORM_SST_MAX=0.10 ATM_RH_STORM_SST_REF=3 ../cli/atm_plm config_wb73a.xml > wb73a.log 2>&1 &
env OMP_NUM_THREADS=6 $K ATM_RH_STORM_SST=0.007 ATM_RH_STORM_SST_MAX=0.08 ATM_RH_STORM_SST_REF=3 ../cli/atm_plm config_wb73b.xml > wb73b.log 2>&1 &
env OMP_NUM_THREADS=6 $K ATM_RH_STORM_SST=0.009 ATM_RH_STORM_SST_MAX=0.10 ATM_RH_STORM_SST_REF=2 ../cli/atm_plm config_wb73c.xml > wb73c.log 2>&1 &
wait
for t in wb73a wb73b wb73c; do
  V=output_$t/0Ma_smooth_Atm_radial_0_60.vtk
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)  slices: $(ls output_$t | grep -c radial_0)"
  echo "banner diff vs wb72 (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb72.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "output\|nm=\|VTK_STRIDE" | head -10
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'water budget closure' $t.log | tail -1 | cut -c1-120
  python3 oceanb.py $V | sed -n '1,4p'
  PARS='0.005,0.04,4,48,1.08' python3 stormfit.py $t 60 | tail -1 | cut -c26-250
  python3 socean.py $t 60 | sed -n '2,14p' | cut -c1-190
  python3 nrows.py $t 60 | sed -n '/by surface temperature/,/^$/p' | cut -c1-150
  python3 landb.py $V | grep "global max"
done
touch WB73_DONE
