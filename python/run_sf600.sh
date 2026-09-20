#!/bin/bash
# B+2's REMAINDER -- clampAndFade's cold fade, 600 from scratch, three arms, one binary.
#
# With ATM_SATADJ_PHASE default ON, this fade is the ONLY non-conservation left in
# SaturationAdjustment: ATM_SATADJ_DIAG on the current default reads neg 0 / supersat 0 / cap 0 /
# fade -4.3075e-04 mm/call with unattributed 0.0000, and ColumnWaterBudget's SaturationAdjust row
# (-3.3e+04 mm/a) agrees independently. That is 36x the model's own precipitation, deleted.
#
# AND IT IS NOT A COLD-CLOUD DEVICE. Measured on output_mn_ctl's zonal slice at iteration 520:
# 92.5 % of the removal happens ABOVE t_00 -- 38.4 % between -27 and -17 C, 26.8 % between -17 and
# 0 C, and 2.6 % ABOVE 0 C, because the sigmoid never reaches 1. Below t_00 it has nothing left to
# do: adjustSaturation's hard freeze already moved that cloud to ice under ATM_ICE_COLD.
#
# THREE ARMS BECAUSE TWO WOULD CONFOUND TWO THINGS:
#   sf_ctl  shipped            -- the liquid is deleted
#   sf_frz  ATM_SATADJ_FADE=1  -- same cells, same alpha, liquid -> ice + latent heat of fusion.
#                                 Conserves MASS only. Isolates "where the water goes".
#   sf_off  ATM_SATADJ_FADE=2  -- no fade. Conserves mass AND phase. The physically right one.
# frz minus ctl is the mass; off minus frz is the phase. Reading only off-vs-ctl would blame the
# whole difference on conservation when most of it may be the phase the fade was mis-assigning.
#
# PRE-REGISTERED, and recorded because this tree records them when they are wrong:
#   1. BOTH REPAIRS RAISE THE CLOUD WATER, and sf_off by more than sf_frz, since sf_frz moves the
#      liquid to ice where sf_off leaves it liquid.
#   2. sf_off ADDS CIRRUS ONLY WEAKLY and sf_frz ADDS IT STRONGLY -- sf_frz is manufacturing ice at
#      -10 to -25 C, which is exactly the band ATM_RH_CRIT_ICE was calibrated in. If IWP moves a
#      lot in sf_frz, that is an argument AGAINST mode 1, not for it.
#   3. PRECIPITATION RISES IN BOTH, because the fade is a sink on the condensate that feeds
#      autoconversion. Whether the starved bands move is the open question -- B+1 has just shown
#      that a water-conservation repair CAN move them where no dynamical change has.
#   4. STABILITY: the fade's stated original purpose was containing the NZ Alps cold-zone
#      supersaturation runaway. Its ROOT FIX (the always-on supersaturation removal) is in and
#      reads EXACTLY 0 on the current default, so removing the fade should be safe -- but watch
#      155, 357, 483 and watch Precip, because exit 0 with zero NaN is not a stability criterion
#      in this tree.
set -u
cd "$(dirname "$0")"
export OMP_NUM_THREADS=8
rm -f SF600_DONE
echo "=== sf600 0->600 from scratch, THREE arms parallel, 8 threads each $(date +%H:%M:%S)"
( env ATM_CWB_DIAG=1 ATM_SATADJ_DIAG=1                    ../cli/atm_sf config_sf_ctl.xml > sf_ctl.log 2>&1
  echo "  sf_ctl exit $?  nan $(grep -ci 'nan' sf_ctl.log)  $(date +%H:%M:%S)" ) &
( env ATM_CWB_DIAG=1 ATM_SATADJ_DIAG=1 ATM_SATADJ_FADE=1  ../cli/atm_sf config_sf_frz.xml > sf_frz.log 2>&1
  echo "  sf_frz exit $?  nan $(grep -ci 'nan' sf_frz.log)  $(date +%H:%M:%S)" ) &
( env ATM_CWB_DIAG=1 ATM_SATADJ_DIAG=1 ATM_SATADJ_FADE=2  ../cli/atm_sf config_sf_off.xml > sf_off.log 2>&1
  echo "  sf_off exit $?  nan $(grep -ci 'nan' sf_off.log)  $(date +%H:%M:%S)" ) &
wait
echo "=== sf600 done $(date +%H:%M:%S)"
touch SF600_DONE
