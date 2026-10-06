#!/bin/bash
# 2026-10-06: wb58 = working branch (the wb57 stack) + ATM_LAND_EVAP=1 (new knob; user: "pursue the land evaporation").
# SCREENING AT nm 60 from scratch, cli/atm_lev (-O2, HEAD 2aeae17 + the knob), 8 threads, ATM_VTK_STRIDE=2; -O0 byte check (run_vlev.sh) alongside.
# WHY (landevap.py on wb53c): land E is 0; ATM_LAND_BUCKET=150 gives land E 728 mm/a against land P 653 (15-35 deg: 911 against 218) because the bucket
# is filled from the OBSERVED rain and evaporates at the open-water rate. The knob takes E from the Budyko curve on the model's own rain (E <= P, E <= E_pot).
# Control = wb57 at 60 (= wb56b rain, wb53b evaporation): 975.8 mm/a, r .607, sigma 1.24, land / ocean 651 / 1104, E 776 (ocean 1083, land 0), P/E 1.25.
# PRE-REGISTERED (offline estimate on wb53c): land E 300-345 mm/a (E/P 0.46-0.53), global E 860-875, P/E 1.11-1.14; rain within 0.3 % of wb57 in the
# global mean, every band, land / ocean and every listed region (wb53c: 728 mm/a of land E moved land rain by 0.2 mm/a); r .607; NaN 0.
# RESULT (11:46-11:50, exit 0, NaN 0): land P 651, potential 1124, land E 325 mm/a (E/P 0.50); global E 867 (ocean 1082), P/E 1.13. Rain 976.1 mm/a
# (-0.2 %), r .607, sigma 1.24, bands 1761/500/1045/265, land / ocean 651 / 1105, every listed region as wb57 -- all inside the pre-registered ranges.
set -u; cd "$(dirname "$0")"; rm -f WB58_DONE
t=wb58
mkdir output_$t || { touch WB58_DONE; exit 1; }
[ -e config_$t.xml ] && { touch WB58_DONE; exit 1; }
sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>60<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>60<#" config_wb7.xml > config_$t.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_lev $(md5sum < ../cli/atm_lev | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_VTK_STRIDE=2 ATM_LAND_EVAP=1 ../cli/atm_lev config_$t.xml > $t.log 2>&1
V=output_$t/0Ma_smooth_Atm_radial_0_60.vtk
echo "== $t  exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)  slices: $(ls output_$t | grep -c radial_0)"
echo "banner diff vs wb57 (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb57.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "nm=\|output" | head -8
grep -a "LAND-EVAP" $t.log | tail -1; grep -a "water budget closure" $t.log | tail -1 | cut -c1-110
grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90
python3 oceanb.py $V
python3 landb.py $V | grep -v "Oman\|Namib\|land \|output"
touch WB58_DONE
