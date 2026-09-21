#!/bin/bash
# THE ATM_RH_MIN_PTOP RE-SWEEP, WITH ATM_MICRO_NDIM=1.0 ON.
#
# WHY IT EXISTS. B+1 (2026-09-20, `9e1aadf`) is the only repair that has ever moved this tree's
# starved precipitation bands -- 35-65 deg +14.6 %, 65-90 deg +101 %, both SATURATING by
# iteration ~450, with `max u-component` identical to six figures. What holds its flip back is
# one coupling and nothing else: `ATM_RH_MIN_PTOP` = 475 hPa was FITTED to land on NASA on the
# branch WITHOUT this correction, and B+1 moves the global mean -3.1 % -> +0.9 %. Same pair
# logic as ATM_ICE_LIMIT_ARRIVING + ATM_RAIN_AREA, and as ATM_SATADJ_PHASE a day earlier.
#
# ⚠ THE SWEEP'S DIRECTION IS REVERSED FROM 2026-09-09, AND THAT IS THE WHOLE POINT.
# That sweep ran 475 / 500 / 525 -- UPWARD -- because SATADJ_PHASE had taken the mean to -4.7 %
# and the question was whether raising PTOP bought it back. It does, and it hands the shape gain
# straight away (sigma 2.25 / 2.48 / 2.66). With MICRO_NDIM on the mean is ALREADY +0.9 % of
# NASA, so there is nothing to buy back. The open question is the other one: does LOWERING PTOP
# -- which takes the RH floor off the 425-675 hPa liquid deck -- cut the tropical spike
# (3267.5 against NASA's 1487.0) and with it sigma and the centred RMS, WITHOUT giving back the
# bands? Arms 425 and 450 ask that; 500 is the direction control.
#
# ⚠ JUDGE THIS SWEEP ON SIGMA AND THE FOUR BANDS, NOT ON THE MEAN. CLAUDE.md, twice:
# "the global mean was never the problem and is not the target -- quote r, sigma and the four
# bands, or quote nothing." Two configurations have matched the mean to 2 % with pattern
# correlations of 0.21 and 0.46. Expect the mean to FALL below NASA in the downward arms; that
# is the knob working, not a regression.
# And DISCOUNT `r` where it disagrees with sigma: a pattern correlation is insensitive to
# amplitude, so it can improve while the field's variance runs away (measured 2026-09-09, where
# PTOP = 500 gave the best r of the sweep at its second-worst sigma).
#
# THE FREE ARM. `output_mn_on` IS the PTOP = 475 arm of this sweep -- same binary `cli/atm_mn`
# (md5 59ffc67b), same config but for output_path, same env but for PTOP, 600 from scratch,
# 12 threads. It is reused rather than re-run, which is why the binary is PINNED to atm_mn and
# not to HEAD's cli/atm: the three knobs added since (SATADJ_FADE, SATADJ_FREEZE_LATENT,
# EVAP_FLUX) are all default 0 and verified byte-identical off, so atm_mn's default branch IS
# HEAD's -- but a sweep whose arms share one binary needs no such argument.
# ATM_CWB_DIAG=1 is carried on the new arms for the same reason: mn_on had it.
#
# PRE-REGISTERED EXPECTATIONS, recorded because this tree records them when they are wrong
# (three of four were wrong yesterday, and run_btref.sh had its sign backwards on 2026-09-14):
#   1. THE MEAN FALLS, MONOTONE, ~90-100 mm/a PER 25 hPa. From the 475->500 step of +96.6
#      measured 2026-09-09. So 450 ~ 890 (-9 %) and 425 ~ 800 (-18 %). If the response is NOT
#      monotone, or is much flatter than this, the knob interacts with the 2783x coefficient and
#      the old sweep does not transfer -- which arm 500 is there to catch.
#   2. SIGMA AND THE CENTRED RMS IMPROVE DOWNWARD. The floor's whole effect below ~500 hPa is
#      LIQUID (measured 2026-08-31: IWP saturates by 500 hPa and never moves again, LWP never
#      saturates), tropical liquid is what the tropical spike is made of, and sigma is that spike
#      measured a second way. This is the arm's hypothesis; if sigma does NOT fall, the spike is
#      not the floor's doing and this whole line is closed.
#   3. ★ THE BANDS HOLD. 35-65 = 183.2 and 65-90 = 20.3 are B+1's gain and the reason for the
#      flip. Across 475->525 in 2026-09-09 the 35-65 band moved 160.4 -> 162.1, i.e. ~1 %, so
#      PTOP is not a storm-track lever in either direction. IF LOWERING PTOP COSTS THE BANDS,
#      THE PAIR IS NOT SEPARABLE AND THAT IS THE FINDING OF THIS SWEEP -- it would mean B+1's
#      gain cannot be kept at a defensible mean, and the flip decision changes shape.
#   4. THE RADIATION IMPROVES DOWNWARD, AND THIS IS A SCORE NOBODY ASKED FOR. The unconditional
#      [cloud-rad DIAG] line reads clear 263.92 / all-sky 234.74 / forcing 29.18 W/m2 in mn_on
#      against Earth's ~265 / ~240 / ~25. Removing liquid should take the forcing DOWN toward 25
#      and the all-sky OLR UP toward 240 -- so expectation 2's shape gain and the radiation
#      should move together. If they move OPPOSITE ways there is a second effect in the floor.
#   5. AND THE FLOOR BELOW WHICH IT COSTS: 2026-08-31 measured PTOP = 400 starving the cirrus to
#      IWP 14.85 against an observed 20-30. 425 is the lowest arm here for that reason. If 425
#      is the best point on shape, the NEXT arm is not 400 -- it is a check on the IWP.
#   6. STABILITY: from scratch through 155, 357 and 483, this tree's three documented failure
#      points, on a configuration whose microphysics source is 2783x the shipped one. Exit 0 and
#      zero NaN is NOT a stability criterion here (the dt ladder reached 8733 mm/a without ever
#      NaN-ing), so read `Precip mean`'s trajectory too.
#
# THREADS: 3 arms x 8 threads = 24, concurrent. Measured yesterday: 24 thr 8.0 s/iter,
# 12 thr 10.0, 8 thr 14.0 -- so three concurrent at 8 is ~140 min for the whole sweep against
# ~240 min run as pairs. Do not re-litigate; and MEASURE the rate from checkpoint mtimes rather
# than quoting this comment.
set -u
cd "$(dirname "$0")"
export OMP_NUM_THREADS=8
rm -f PTOP600_DONE
echo "=== ptop600 sweep 0->600 from scratch, 3 arms concurrent, 8 threads each $(date +%H:%M:%S)"
for V in 425 450 500; do
  ( env ATM_CWB_DIAG=1 ATM_MICRO_NDIM=1.0 ATM_RH_MIN_PTOP=$V \
        ../cli/atm_mn config_pt${V}.xml > pt${V}.log 2>&1
    echo "  pt${V} exit $?  nan $(grep -ci 'nan' pt${V}.log)  $(date +%H:%M:%S)" ) &
done
wait
echo "=== ptop600 done $(date +%H:%M:%S)"
touch PTOP600_DONE
