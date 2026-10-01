#!/bin/bash
# 2026-10-01: SNOW-SUBL -- ATM_SNOW_DEP_FLUX 0 vs 1 on the working branch (incl. ATM_RH_STORM=1.15), 600 from scratch,
# one -O2 binary cli/atm_sdx, 2 x 8 threads concurrent. Waits for the -O0 byte check (VSDX_DONE) first.
# PRE-REGISTERED: sublimation in the arrays unchanged in sign but the rate arrays then sum to the ground flux
# (SUM S_r+S_s ~ +250 instead of -5515); the CWB microphysics row turns from +5.5e3 to about -(stratiform P);
# column-water gain falls from ~7200 mm/a toward E - P; 0-15 and 35-65 snow at ground falls; precipitation most
# likely FALLS (less vapour made aloft), r/sigma/bands to be read, not predicted.
set -u; cd "$(dirname "$0")"; rm -f SDX_DONE
while [ ! -e VSDX_DONE ]; do sleep 30; done
for t in ctl on; do
  mkdir output_sdx_$t || { touch SDX_DONE; exit 1; }
  [ -e config_sdx_$t.xml ] && { touch SDX_DONE; exit 1; }
  sed "s#output_rhst_115/#output_sdx_$t/#" config_rhst_115.xml > config_sdx_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_sdx $(md5sum < ../cli/atm_sdx | cut -c1-8)"
D="ATM_SNOW_DIAG=1 ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_CWB_BANDS=1 ATM_SR_DIAG=1"
env OMP_NUM_THREADS=8 $D ATM_SNOW_DEP_FLUX=0 ../cli/atm_sdx config_sdx_ctl.xml > sdx_ctl.log 2>&1 &
env OMP_NUM_THREADS=8 $D ATM_SNOW_DEP_FLUX=1 ../cli/atm_sdx config_sdx_on.xml  > sdx_on.log  2>&1 &
wait
for t in ctl on; do
  echo "== sdx_$t  exit-log NaN $(grep -c 'NaN/Inf DETECTED' sdx_$t.log)  $(date +%H:%M)"
  grep -a "by |latitude|" sdx_$t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" sdx_$t.log | tail -1 | cut -c1-150
  grep -a -i "land .*ocean" sdx_$t.log | tail -1 | cut -c1-100; grep -a 'water budget closure' sdx_$t.log | tail -1
  grep -a "P_rain mean\|P_snow mean\|P_conv mean" sdx_$t.log | tail -3
  grep -a "SNOW DIAG\] GLOBAL" sdx_$t.log | tail -1
  grep -a "as RK4 applies them\|RungeKutta bucket\|\[CWB\] NET\|rest = MC_q" sdx_$t.log | tail -4
  grep -a 'precipitable water average' sdx_$t.log | tail -1 | cut -c1-60
  grep -a 'max v-component' sdx_$t.log | tail -1 | cut -c1-80; tail -1 output_sdx_$t/convergence.csv
done
touch SDX_DONE
