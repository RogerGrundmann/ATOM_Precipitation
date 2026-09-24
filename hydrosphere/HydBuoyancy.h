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
