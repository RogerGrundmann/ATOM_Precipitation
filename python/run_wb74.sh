#!/bin/bash
# *** PROPOSAL, NOT RUN YET (2026-10-08). ***
# wb74 = wb73b at 600 from scratch = working branch (wb72 stack, 60 knobs) with ATM_RH_STORM_SST 0.005 -> 0.007, _MAX 0.04 -> 0.08, _REF 4 -> 3 and
# ATM_RH_STORM_POLAR 1.08 -> 1.0 (the polar floor off). cli/atm_plm (-O2, = HEAD, md5 fe6987f6), 8 threads.
# PRE-REGISTERED from the nm 60 screen wb73b (977.4 mm/a, r .652, sigma 1.08, bands 1666/608/1002/322, ocean 35-65 / 65-90 1137 / 374, P/E 1.01):
#   global -0.5 to +0.5 % | r >= .650 | sigma 1.06-1.10 | ocean 35-65 1120-1155 | ocean 65-90 S 400-440, N 320-355 | no ocean row 46-74 deg above 1.25x NASA
#   row misfit (stormfit.py) <= 0.15 | land 35-65 / 65-90 625-655 / 250-280 | P/E 1.00-1.05 | no cell > 12 mm/d
#   HOLDS if also: zero NaN, no drift (global mean within 0.5 % over 100-600), max |v| not growing.
set -u; cd "$(dirname "$0")"; rm -f WB74_DONE; t=wb74
mkdir output_$t || { touch WB74_DONE; exit 1; }
[ -e config_$t.xml ] && { touch WB74_DONE; exit 1; }
sed -e "s#output_wb7/#output_$t/#" config_wb7.xml > config_$t.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_plm $(md5sum < ../cli/atm_plm | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_RH_STORM_SST=0.007 ATM_RH_STORM_SST_MAX=0.08 ATM_RH_STORM_SST_REF=3 ATM_RH_STORM_POLAR=1.0"
env OMP_NUM_THREADS=8 $K ../cli/atm_plm config_$t.xml > $t.log 2>&1
echo "== $t  exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
echo "banner diff vs wb72 (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb72.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "output"
echo "-- trajectory (every 100: bands; score; land / ocean)"
grep -a "by |latitude|" $t.log | grep -v MFC | awk 'NR%100==0' | cut -c1-150
grep -a "model .*NASA .*bias" $t.log | awk 'NR%100==0' | cut -c1-150
grep -a 'land .*ocean .*(model / NASA)' $t.log | awk 'NR%100==0' | cut -c1-90
echo "-- final"
grep -a "by |latitude|" $t.log | grep -v MFC | tail -2; grep -a "model .*NASA .*bias" $t.log | tail -2 | cut -c1-150
grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1 | cut -c1-110; grep -a "LAND-EVAP" $t.log | tail -1
grep -a "P_conv mean\|P_rain mean" $t.log | tail -2
grep -a "max v-component\|max w-component\|max u-component" $t.log | awk 'NR%300<3' | cut -c1-170
tail -1 output_$t/convergence.csv
for it in 120 520; do V=output_$t/0Ma_smooth_Atm_radial_0_$it.vtk; echo "-- slice $it"; python3 oceanb.py $V; python3 landb.py $V | grep -v "output"; done
python3 row36.py $t 520 | grep -v "T [0-9]" | cut -c1-230
ITER=520 python3 pacband.py $t | grep -v "^ \+-\?[0-9]\+:\|ocean by latitude"
G=5 python3 polarland.py $t 520 | sed -n '3,$p' | cut -c1-200
python3 polar.py $t 520 | sed -n "2,19p" | cut -c1-110
PARS='0.005,0.04,4,48,1.08' python3 stormfit.py $t 520 | tail -1 | cut -c26-250
python3 socean.py $t 520 | sed -n '2,14p' | cut -c1-190
python3 nrows.py $t 520 | cut -c1-150
touch WB74_DONE
