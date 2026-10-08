#!/bin/bash
# RUN 2026-10-08 (user: "run wb71, then the 600, adopt and push when it holds"): wb72 = wb71a at 600 from scratch = working branch (wb68 stack, 56 knobs)
# + ATM_MC_CMB_OCEAN=0.85 + ATM_EVAP_GUST=5 (wb69a) + ATM_HCRIT_SFC_POLAR=1 + ATM_RH_LAND_ML=1500 + ATM_RH_LAND_ML_STRENGTH=0.45 (wb71a).
# cli/atm_plm (-O2, HEAD + the polar-land knobs, md5 fe6987f6), 8 threads.
# PRE-REGISTERED from the nm 60 screen wb71a (984.4 mm/a, r .650, sigma 1.10, bands 1666/608/1021/326, land 35-65 / 65-90 641 / 264, P/E 1.02):
#   global +0.2 to +1.0 % | r >= .648 | sigma 1.08-1.12 | ocean 0-15 / 15-35 1640-1690 / 740-770 | land 35-65 625-655 | land 65-90 250-280 | 65-90 band 315-340
#   P/E 1.00-1.04 | N polar land 0-200 m 420-470, 1000-2000 m 440-500 | E Antarctica <= 160 | N Europe <= 1150 | no cell > 12 mm/d
#   HOLDS if also: zero NaN, no drift (global mean within 0.5 % over 100-600), max |v| not growing (the new code touches no dynamics).
set -u; cd "$(dirname "$0")"; rm -f WB72_DONE; t=wb72
mkdir output_$t || { touch WB72_DONE; exit 1; }
[ -e config_$t.xml ] && { touch WB72_DONE; exit 1; }
sed -e "s#output_wb7/#output_$t/#" config_wb7.xml > config_$t.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_plm $(md5sum < ../cli/atm_plm | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_MC_CMB_OCEAN=0.85 ATM_EVAP_GUST=5 ATM_HCRIT_SFC_POLAR=1 ATM_RH_LAND_ML=1500 ATM_RH_LAND_ML_STRENGTH=0.45"
env OMP_NUM_THREADS=8 $K ../cli/atm_plm config_$t.xml > $t.log 2>&1
echo "== $t  exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
echo "banner diff vs wb68 (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb68.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "output"
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
touch WB72_DONE
