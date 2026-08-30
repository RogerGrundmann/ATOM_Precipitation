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

    // Parabola H_crit(p): roots at p = 0 and p = p_crit, minimum critMid() at
    // p = p_mid. The same curve initCloudIce builds from Hu_cr_max / Hu_curv.
    inline double hCrit(double p_hPa){
        constexpr double p_crit = 1000.0, x_mid = 0.55;
        const double curv = (1.0 - critMid()) / (x_mid * (1.0 - x_mid));
        const double x    = p_hPa / p_crit;
        const double h    = 1.0 - curv * x * (1.0 - x);
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
