#!/bin/bash
# 2026-10-05: wb46a / wb46b = working branch (the wb45 stack) + ATM_RH_SIGMA_SFC=1 (user: "try option 1"; SUBTROP), ATM_HCRIT_SFC_LAT 25 (a) / 35 (b).
# SCREENING AT nm 60 from scratch, cli/atm_ssc (-O2, md5 e69847fd), 2 x 8 threads, ATM_VTK_STRIDE=2. No code change.
# WHY (plateau28.py on wb45, 28N): plateau columns find a cloud base but miss buoyancy by 8-14 K theta_e. (a) the air aloft over elevated ground is
# 2-6 K warmer than over the neighbours (each column's profile is built from its own ground temperature); (b) the plateau boundary layer is dry --
# Mexican plateau ML RH 0.65 against 0.84 needed, N India 0.44 -- partly because Manabe-Wetherald's sigma = p/p_0 dries elevated ground by ~0.86.
# ATM_RH_SIGMA_SFC=1 uses the column's own surface pressure on every land column. It flooded the highlands on wb15 through cloud AT the ground,
# which ATM_HCRIT_SFC has since removed -- but only equatorward of 25 deg (taper to 35): arm b extends that to 35 (taper to 45).
# Control = wb45 (600; = wb44a at 60): 975.4 mm/a (-0.3 %), r .605, sigma 1.26, land 0-15 / 15-35 / 35-65 / 65-90 1627 / 319 / 643 / 253,
# Mexico plt 1.13 (NASA 1.98), Highveld 1.23 (1.72), S Brazil 1.69 (4.55), E Africa 2.9 (2.58), Congo 4.79, E Austral 6.07, Arabia 0.17, max cell 18.2.
# PRE-REGISTERED (at 60): a -- Mexican plateau ML RH 0.65 -> 0.74-0.78 (still below 0.84: convective rain there 0-0.5 mm/d), tropical highlands
# (E Africa, Congo) convective +0.3-1.5 mm/d, land 0-15 1650-1800, land 15-35 330-430, land 35-65 643 -> 700-950 with cells > 20 mm/d on elevated land
# at 30-45 deg (Tibet, Iran, Rockies: no terrain-following threshold there), global 985-1040, r .59-.61; b -- land 35-65 600-700, land 15-35 300-360,
# no land cell > 20 mm/d poleward of 25 deg, global 975-1000.
# RISK: the highland flood of wb15 returns between 25 and 45 deg (a); arm b removes the step-1 land rain at 35-45 deg.
# RESULT (13:55-13:59, NaN 0): NOT USABLE. a: 1354 mm/a (+38.4 %), r .173, sigma 2.75, land 35-65 643 -> 3000+ (land at 40N 23.9 mm/d vs 1.5), cell 50 mm/d at 44N 77E;
# b: 1184 (+21.1 %), r .349, land at 44-52N still 7-15 mm/d (the knob also moistens the cold low-pressure columns at 35-52 deg, outside any HCRIT limit).
# Mexican plateau: ML RH 0.65 -> 0.74, margin -8.0 -> -4.2 K, still no convection (as predicted). Tropical highlands DO convect: land < 15 deg at
# 800-1300 m conv 0.74 -> 2.82 mm/d (NASA total 3.44), 1300-2000 m 0.05 -> 0.82 -- but a stratiform deck returns above the ground (0.03 -> 2.80; E Africa 8.1 vs 2.58).
set -u; cd "$(dirname "$0")"; rm -f WB46_DONE
for t in wb46a wb46b; do
  mkdir output_$t || { touch WB46_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB46_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>60<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>60<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_ssc $(md5sum < ../cli/atm_ssc | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_VTK_STRIDE=2 ATM_RH_SIGMA_SFC=1"
env OMP_NUM_THREADS=8 $K ../cli/atm_ssc config_wb46a.xml > wb46a.log 2>&1 &
env OMP_NUM_THREADS=8 $K ATM_HCRIT_SFC_LAT=35 ../cli/atm_ssc config_wb46b.xml > wb46b.log 2>&1 &
wait
for t in wb46a wb46b; do
  V=output_$t/0Ma_smooth_Atm_radial_0_60.vtk
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)  slices: $(ls output_$t | grep -c radial_0)"
  echo "banner diff vs wb45 (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb45.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "nm=\|output" | head -8
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90
  python3 oceanb.py $V
  python3 landb.py $V | grep -v "Oman\|Namib\|W Austral\|land \|output"
  python3 landlat.py $V
  python3 plateau28.py $t 60 | sed -n '2,12p;13,20p' | cut -c1-200
done
ITER=60 python3 tropland.py wb44a wb46a wb46b
touch WB46_DONE
