#pragma once

#include "Utils.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>

// ============================================================================
// CloudFraction — the uniform-PDF (Smith / Sundqvist) sub-grid cloud closure.
//
// Three modules now need the SAME f and the SAME H_crit: initCloudIce creates
// the condensate with it, SaturationAdjustment must not evaporate what the
// closure created, and MultiLayerRadiation must not treat a cell of fraction f
// as overcast. It was written twice already (InitValues_Atm.cpp and
// SaturationAdjustment.h) and a third copy in the radiation would be the point
// at which the three drift apart, so it lives here once.
//
// Total water is uniform over [q_t - D, q_t + D] with D = (1 - H_crit)*q_sat,
// and the cell is cloudy where the local total exceeds saturation:
//
//     s   = q_t + D - q_sat = q_sat*(RH - H_crit)
//     f   = clamp(s / 2D, 0, 1)
//     q_c = f^2 * D          (f < 1),      q_t - q_sat   (f = 1)
//
// f = 0 exactly at RH = H_crit and rises smoothly; q_c is the GRID-MEAN
// condensate, so the IN-CLOUD value q_c/f stays finite as f -> 0 (q_c ~ f^2).
// One parameter, H_crit, which the shipped scheme already had.
//
// EVERYTHING HERE REDUCES TO THE SHIPPED GRID-MEAN LIMIT AT H_crit = 1: D = 0,
// f = 1, q_c = max(0, q_t - q_sat). That is why the off-branch is unchanged by
// CONSTRUCTION rather than by a guard, and it is the property to preserve in
// anything added below.
// ============================================================================

namespace CloudFraction {

    // ATM_CLOUD_FRAC=1 — use the sub-grid closure instead of the grid-mean
    // supersaturation. Default 0 = shipped. See CLAUDE.md: this knob is one of
    // a set (ATM_RH_PROFILE, ATM_RH_CRIT, ATM_CWP_CAP) that must move together.
    inline bool enabled(){
        static const bool v = [](){
            const char* e = getenv("ATM_CLOUD_FRAC"); return e && atoi(e) != 0; }();
        return v;
    }

    // ATM_RH_CRIT — the critical-humidity midpoint. Default 0.8 = shipped.
    inline double critMid(){
        static const double v = [](){
            const char* e = getenv("ATM_RH_CRIT");
            const double x = e ? atof(e) : 0.8;
            return (x > 0.0 && x < 1.0) ? x : 0.8; }();
        return v;
    }

    // ATM_RH_CRIT_ICE -- a SEPARATE critical humidity for the ice-cloud branch aloft.
    // Default 0 = disabled: hCrit() is exactly the single-parameter curve below.
    //
    // WHY IT IS NEEDED, MEASURED. With one H_crit the closure sets the liquid deck and the
    // cirrus with the same number, so ice can only be bought with liquid at a fixed rate: the
    // ATM_RH_MIN sweep gives LWP 74.7 -> 98.9 -> 263.6 as IWP goes 5.0 -> 13.5 -> 45.3, and the
    // three targets never agree -- the radiation is best at RH_MIN 0.40 (all-sky OLR 242.7
    // against ~240), the precipitation wants ~0.45-0.47, and the IWP wants higher still.
    // Splitting the threshold decouples them: q_c = f^2*D with f = (RH - H_crit)/(2(1-H_crit))
    // and D = (1-H_crit)*q_sat, so lowering H_crit ALOFT raises the ice at fixed humidity while
    // leaving the liquid deck below p_mid untouched.
    //
    // Physically the discriminator is the PHASE, but hCrit is a function of pressure in every
    // call site in the tree, so the split is made in pressure: the parabola is unchanged for
    // p >= p_mid, then ramps linearly from critMid() at p_mid to critIce() at p_ice = 300 hPa
    // and stays flat above. Continuous at p_mid by construction, because the parabola's minimum
    // IS critMid() there -- a jump in H_crit would be a jump in cloud fraction.
    inline double critIce(){
        static const double v = [](){
            const char* e = getenv("ATM_RH_CRIT_ICE");
            const double x = e ? atof(e) : 0.0;
            return (x > 0.0 && x < 1.0) ? x : 0.0; }();
        return v;
    }

    // Parabola H_crit(p): roots at p = 0 and p = p_crit, minimum critMid() at
    // p = p_mid. The same curve initCloudIce builds from Hu_cr_max / Hu_curv.
    //
    // FLAT ABOVE THE MINIMUM, BECAUSE THE PARABOLA'S UPPER ROOT IS AN ARTEFACT OF THE FIT.
    // The curve was built to reproduce the shipped Hu_cr profile over the TROPOSPHERE, and a
    // parabola with roots at p = 0 and p = 1000 necessarily returns to H_crit = 1 at the top:
    // measured on the flip arm it runs 0.298 at 6.1 km, 0.362 at 8.2 km, 0.499 at 10.9 km and
    // 0.683 at 14.6 km -- it climbs back toward saturation exactly where the humidity is
    // falling (RH 0.32 / 0.26 / 0.19 / 0.17 at those levels). The two curves cross at ~4.5 km
    // and diverge above it, so NO cell aloft can be cloudy and the ice water path is 1.15 g/m2
    // against an observed 20-30.
    //
    // A threshold going to 1 at p -> 0 also says the sub-grid variability of total water
    // VANISHES at the tropopause, which is backwards, and it contradicts the closure's own
    // purpose: H_crit = 1 is the grid-mean limit this scheme exists to replace. Every
    // operational form (Sundqvist, ECMWF) has RH_crit decreasing from ~1 at the ground to a
    // constant in the free troposphere. Holding x at x_mid below p_mid does exactly that,
    // leaves the calibrated tropospheric branch untouched, and still degenerates to H_crit = 1
    // everywhere as critMid -> 1, so the grid-mean limit is preserved.
    inline double hCrit(double p_hPa){
        constexpr double p_crit = 1000.0, x_mid = 0.55, p_ice = 300.0;
        const double curv = (1.0 - critMid()) / (x_mid * (1.0 - x_mid));
        double x = p_hPa / p_crit;
        if (x < x_mid) x = x_mid;                     // no upturn toward the lid
        double h = 1.0 - curv * x * (1.0 - x);        // == critMid() for p <= p_mid
        const double h_ice = critIce();
        constexpr double p_mid = x_mid * p_crit;
        if (h_ice > 0.0 && p_hPa < p_mid) {           // ATM_RH_CRIT_ICE: ramp to the ice branch
            double u = (p_mid - p_hPa) / (p_mid - p_ice);
            if (u < 0.0) u = 0.0; else if (u > 1.0) u = 1.0;
            h = critMid() + u * (h_ice - critMid());
        }
        return (h > 1.0) ? 1.0 : ((h < 0.0) ? 0.0 : h);
    }

    // Cloudy area fraction of the cell, from TOTAL water q_t = q_v + q_c + q_i.
    inline double fraction(double q_t, double q_s, double p_hPa){
        const double D = (1.0 - hCrit(p_hPa)) * q_s;
        if (!(D > 0.0)) return (q_t > q_s) ? 1.0 : 0.0;   // H_crit = 1 -> overcast or clear
        const double f = (q_t + D - q_s) / (2.0 * D);
        return (f <= 0.0) ? 0.0 : ((f >= 1.0) ? 1.0 : f);
    }

    // Grid-mean condensate the closure supports for total water q_t.
    inline double qcEquilibrium(double q_t, double q_s, double p_hPa){
        const double D = (1.0 - hCrit(p_hPa)) * q_s;
        if (!(D > 0.0)) return std::max(0.0, q_t - q_s);  // H_crit = 1 -> grid-mean limit
        const double f = fraction(q_t, q_s, p_hPa);
        if (f <= 0.0) return 0.0;
        if (f >= 1.0) return std::max(0.0, q_t - q_s);
        return f * f * D;
    }

    // THE FRACTION A PROCESS SHOULD SEE WHEN CONDENSATE IS PRESENT.
    //
    // fraction() returns 0 wherever the grid mean sits at or below H_crit, and a cell can hold
    // condensate there anyway -- advected in, or one the adjustment has since dried. Dividing a
    // real condensate by a zero fraction is not physical, so the floor is the fraction that
    // would hold this much water under the same closure: q_c = f^2*D, hence f = sqrt(q_c/D).
    //
    // Returns 1 when the closure is off, when there is no condensate, or when D degenerates, so
    // every caller reduces EXACTLY to its grid-mean form on the shipped branch. That is the
    // property that keeps the off-branch bit-identical without a guard at each call site.
    inline double effectiveFraction(double q_t, double q_s, double p_hPa, double q_cond){
        if (!enabled() || !(q_cond > 0.0)) return 1.0;
        double f = fraction(q_t, q_s, p_hPa);
        if (f > 0.0) return (f > 1.0) ? 1.0 : f;
        const double D = (1.0 - hCrit(p_hPa)) * q_s;
        if (!(D > 0.0)) return 1.0;
        f = std::sqrt(q_cond / D);
        if (!(f > 0.0)) return 1.0;
        return (f > 1.0) ? 1.0 : f;
    }

    // Saturation specific humidity, branched water/ice exactly as initCloudIce
    // branches it — the closure has to see the same q_sat that made the field.
    inline double qSat(double T_K, double p_hPa, double t_0, double hp, double ep){
        const double E_sat = (T_K >= t_0)
            ? hp * AtomUtils::exp_func(T_K, 17.2694, 35.86)
            : hp * AtomUtils::exp_func(T_K, 21.8747, 7.66);
        return (p_hPa > E_sat) ? ep * E_sat / (p_hPa - E_sat) : ep * 1e-5;
    }
}

// ============================================================================
// ColdCloud -- ATM_ICE_COLD: let ICE exist below t_00.
//
// t_00 = 236.15 K (-37 C) is the homogeneous-freezing temperature: the point below which
// supercooled LIQUID cannot exist. This tree applies it to the ice as well, in five places,
// and one of them deletes the water vapour too:
//
//   InitValues_Atm.cpp:1265   if (t_u <= t_00) { cloud = ice = gr = 0; }
//   SaturationAdjustment:251  if (T   <  t_00) { cloud = ice = 0; }
//   SaturationAdjustment:340  a logistic fade to 0 below t_00 applied to c, cloud AND ice
//   SaturationAdjustment:147  the same fade on the adjustment's entry weight
//   UtilsAtm.h:254            if (t_u <= t_00) { c = cloud = ice = gr = 0; }
//
// SO THE MODEL HAS NO CIRRUS BY CONSTRUCTION: cirrus lives at -40 to -70 C, and the whole
// water cycle is switched off below -37 C. It is latent today only because this tree's upper
// troposphere is ~20 K TOO WARM -- the lid reaches only -37.0 C, so the cutoff sits at the top
// of the model instead of at 8-9 km where a correct profile would put it. Repair the radiation
// (ATM_RAD_EQUIL + ATM_SW_INSOL, both still default off) and this cutoff removes EVERY cloud
// above ~8 km. The warm bias and the cutoff are a cancelling pair, and fixing either alone
// makes the cloud field worse -- the same shape as ATM_RAD_EQUIL / ATM_SW_INSOL themselves.
//
// Below t_00 the physics is: liquid cannot exist, so freeze it into ice rather than deleting
// it; deposition continues (SaturationAdjustment's own CND/DEP split already sends 100 % to
// the ice branch there); and vapour is not a condensate and must never be zeroed at all.
// Default 0 = shipped, every site unchanged.
// ============================================================================
namespace ColdCloud {

    inline bool enabled(){
        static const bool v = [](){
            const char* e = getenv("ATM_ICE_COLD"); return e && atoi(e) != 0; }();
        return v;
    }
}
