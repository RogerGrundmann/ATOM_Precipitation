#pragma once

#include <vector>
#include <cstdlib>

// ======================================================================================
// HYD_BUOY_CONSISTENT=<strength>, default 0.0 = shipped and byte-identical.
//
// The ocean's only route from temperature or salinity into momentum is `buoy_nd` in rhs_u,
// and it is `buoyancy * g * dt / u_0 * ((t-1) - alpha_S*(c-1))`. rhs_u is a TENDENCY that
// RungeKutta_Hyd_Turb already multiplies by dt, so the `dt` in that coefficient is a second,
// spurious one -- the pre-5f71571 convention the comment beside it names -- and the term is
// ~1e-7 too small: buoyancy is effectively OFF and the ocean is dynamically barotropic.
// Same defect as the atmosphere's ATM_BUOY_CONSISTENT (5.0e5) and surf_drag (4.0e5).
//
// This knob replaces the coefficient with the consistent Boussinesq one and the anomaly
// with a DENSITY anomaly from the model's own equation of state:
//
//     b_nd = - s * (g * L_hyd / u_0^2) * (rho(T,S) - rho_ref[i]) / r_0_water
//
//   * rho(T,S) is the Gill surface form ThermoHyd uses for r_salt_water (C_p, beta_p,
//     alfa_t_p(T), gamma_t_p(T)), WITHOUT the pressure terms. Those depend on depth alone,
//     so they are removed by the level reference below and would only add round-off.
//   * rho_ref[i] is the sin(theta)-area-weighted mean over WATER cells of level i. The
//     level-mean part of the buoyancy is a function of depth alone -- the hydrostatic
//     reference -- and belongs in the pressure, not in rhs_u. The shipped term measures
//     against a GLOBAL constant (t = 1, c = 1) instead, so on a stratified ocean its
//     horizontally uniform part is the whole stratification. Subtracting the level mean is
//     what the atmosphere's t_ref_level does for the same reason. What is left is the
//     HORIZONTAL density anomaly at each depth -- exactly the part that drives baroclinic
//     overturning and cannot be absorbed by a pressure that depends on depth alone.
//
// The reference is built on the TIME-LEVEL-n fields tn/cn, before the sweep, for the reason
// HydHorizViscosity.h gives: this RK4 is pointwise, so running t/c at another cell are at a
// sweep-order-dependent stage. Each level is summed by one thread, so the reference is
// deterministic. The local rho uses the running t/c of the cell itself.
//
// WHY A STRENGTH AND NOT A FLAG. The full repair was TESTED on 2026-07-14 and destabilised:
// max radial velocity 0.04 -> 19.5 m/s in 30 iterations. That was BEFORE the Dirichlet u BC
// at both radial walls (575032f), HYD_RUN_NEUMANN, HYD_BAROCLINIC_PGF and the level
// reference here, so it is a different model -- but the endpoint is ~5e5 x the shipped
// coefficient and a strength lets the arm find where it stops being safe.
//
// Scratch lives in inline variables, not as members of cHydrosphereModel: adding members
// moves sizeof(cHydrosphereModel), the stack-canary hazard (CLAUDE.md, "The build hazard").
// ======================================================================================
namespace HydBuoy {

inline double strength(){
    static const double s = [](){
        const char* e = std::getenv("HYD_BUOY_CONSISTENT"); return e ? std::atof(e) : 0.0; }();
    return s;
}

inline std::vector<double> rho_ref;          // per level, kg/m^3, water cells only

// Gill (1982) surface form, as ThermoHyd builds r_salt_water, pressure terms omitted.
inline double rho_eos(double t_Celsius, double S_psu){
    const double C_p       = 999.83;
    const double beta_p    = 0.808;
    const double alfa_t_p  = 0.0708 * (1.0 + 0.068 * t_Celsius);
    const double gamma_t_p = 0.003  * (1.0 - 0.012 * t_Celsius);
    return C_p + beta_p * S_psu - alfa_t_p * t_Celsius
         - gamma_t_p * (35.0 - S_psu) * t_Celsius;
}

} // namespace HydBuoy

// ======================================================================================
// HYD_HYDRO_SPLIT=<strength>, default 0.0 = OFF and byte-identical: nothing below is read.
// The ocean analogue of ATM_HYDRO_SPLIT (atmosphere/AtmHydroSplit.h), 2026-09-26.
//
// WHY. HYD_BUOY_CONSISTENT gives the buoyancy its correct size and lets the projection pressure
// answer it. It does not: bn_* (2026-09-25) measured an unphysical RADIAL velocity exactly linear
// in the strength, rms 7.4e-3 m/s at 1.0 (3000x the control) -- the atmosphere's "the pressure
// cannot answer a body force" in the ocean. Hydrostatic balance needs no elliptic solve; it is a
// COLUMN INTEGRAL. So split p = p_hb + p_nh with dp_hb/dz = -g rho'/r_0 exactly, and then
//     rhs_u :  b - dp_hb/dz  = 0     -- not formed; nothing is left to drive a radial runaway
//     rhs_v :  - (1/R) dp_hb/dthe    -- the baroclinic pressure gradient (thermal wind)
//     rhs_w :  - (1/(R sin)) dp_hb/dphi
// with p_dyn from the Poisson solve playing p_nh. The horizontal force stays inside aux_v/aux_w
// (as HYD_BAROCLINIC_PGF's does), so the projection removes only its divergent part.
//
// THE DENSITY IS HYD_BUOY_CONSISTENT's: rho_eos(T,S) against rho_ref[i], the per-level water
// mean, times the strength -- so the radial balance is exact for the same anomaly the radial
// term would have used. HYD_BAROCLINIC_PGF differences p_hydro, a DIFFERENT field (full Gill with
// pressure terms, integrated in bar); setting both double-counts the thermal wind, and the
// banner warns.
//
// THE INTEGRAL RUNS DOWN FROM THE SURFACE, p_hb(im-1) = 0 -- the rigid-lid ocean convention: the
// surface pressure is the barotropic mode, which project_barotropic / the mode split impose
// separately. (The atmosphere integrates UP from the ground for the mirror-image reason.) Physical
// depth: dz = (rad.z[i+1] - rad.z[i]) * L_hyd. Non-dimensional p_hb = p'/(r_0 u_0^2), so
//     dp_hb = (g / u_0^2) * (rho' / r_0) * dz[m]      (denser water -> pressure grows downward)
//
// THE HORIZONTAL METRIC IS THE EARTH'S, as in HYD_BAROCLINIC_PGF and for the same reason: the
// force must balance Coriolis, and the ocean's inv_rm is 2e4 off (HYD_METRIC_RADIUS default 0).
// Non-dimensional tendency = (L_hyd / R_Earth[m]) * dp_hb/dthe. Land neighbours: one-sided
// against the cell itself, zero if both are land (never difference against rock).
//
// TIME LEVEL: tn/cn, once per iteration before the pointwise RK4 sweep (deterministic), exactly
// like rho_ref. Storage is namespace-scope, not a member (sizeof / stack-canary hazard).
// ======================================================================================
namespace HydSplit {

inline double strength(){
    static const double s = [](){
        const char* e = std::getenv("HYD_HYDRO_SPLIT"); return e ? std::atof(e) : 0.0; }();
    return s;
}
inline bool enabled(){ return strength() != 0.0; }

inline std::vector<double> p_hb;             // im*jm*km, (i*jm + j)*km + k, non-dim p'/(r_0 u_0^2)
inline int n_j = 0, n_k = 0;
inline double at(int i, int j, int k){ return p_hb[((size_t)i * n_j + j) * n_k + k]; }

} // namespace HydSplit
