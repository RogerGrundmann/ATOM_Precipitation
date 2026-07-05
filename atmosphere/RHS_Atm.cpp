/*
 * Atmosphere General Circulation Modell (AGCM) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to combine the right hand sides of the differential equations for the Runge-Kutta scheme
*/

#include <iostream>
#include <cmath>
#include "cAtmosphereModel.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;

// minmod slope limiter on two undivided differences (see RHS_Atm_Turb.cpp).
// Returns 0 at a local extremum so the limited advective gradient is
// positivity-preserving for the moisture scalars.
static inline double minmod(double a, double b){
    return (a * b <= 0.0) ? 0.0 : (std::fabs(a) < std::fabs(b) ? a : b);
}


void cAtmosphereModel::RHS_Atmosphere(int i, int j, int k, const CellGeometry& geo){

    // All geometric quantities come from the precomputed struct —
    // NO sin(), cos(), division, or reciprocal computation here.
    const double exp_rm        = geo.exp_rm;
    const double exp_2_rm      = geo.exp_2_rm;
    const double sinthe        = geo.sinthe;
    const double sinthe2       = geo.sinthe2;
    const double costhe        = geo.costhe;
    const double inv_rm        = geo.inv_rm;
    const double inv_rm2       = geo.inv_rm2;
    const double inv_rmsinthe  = geo.inv_rmsinthe;
    const double inv_rm2sinthe = geo.inv_rm2sinthe;
    const double inv_rm2sinthe2 = geo.inv_rm2sinthe2;
    const double costhe_inv_rm2sinthe = geo.costhe_inv_rm2sinthe;

    const double inv_2dr   = geo.inv_2dr;
    const double inv_2dthe = geo.inv_2dthe;
    const double inv_2dphi = geo.inv_2dphi;
    const double inv_dr2   = geo.inv_dr2;
    const double inv_dthe2 = geo.inv_dthe2;
    const double inv_dphi2 = geo.inv_dphi2;

    // Cache local cell velocities
    double u_ijk = u.x[i][j][k];
    double v_ijk = v.x[i][j][k];
    double w_ijk = w.x[i][j][k];

    bool land_ijk = is_land(h, i, j, k);

    bool r_flag = false, the_flag = false, phi_flag = false;

    // We compute derivatives for 10 first-order arrays and 9 second-order arrays.
    // Instead of looping over Array* with virtual indirection, compute directly.

    // Helper lambdas to compute stencil values for ALL variables at once.
    // Each reads from the 3D arrays directly — no pointer-array indirection.

    // ---- First-order derivatives storage ----
    double dudr, dvdr, dwdr, dtdr, dpdr, dcdr, dclouddr, dicedr, dgdr, dcodr;
    double dudthe, dvdthe, dwdthe, dtdthe, dpdthe, dcdthe, dclouddthe, dicedthe, dgdthe, dcodthe;
    double dudphi, dvdphi, dwdphi, dtdphi, dpdphi, dcdphi, dclouddphi, dicedphi, dgdphi, dcodphi;

    // ---- Second-order derivatives storage ----
    double d2udr2, d2vdr2, d2wdr2, d2tdr2, d2cdr2, d2clouddr2, d2icedr2, d2gdr2, d2codr2;
    double d2udthe2, d2vdthe2, d2wdthe2, d2tdthe2, d2cdthe2, d2clouddthe2, d2icedthe2, d2gdthe2, d2codthe2;
    double d2udphi2, d2vdphi2, d2wdphi2, d2tdphi2, d2cdphi2, d2clouddphi2, d2icedphi2, d2gdphi2, d2codphi2;


    // ===== R-direction derivatives =====
    if(i < im-2 && land_ijk && is_air(h, i+1, j, k)){
        // Forward-biased stencil near land-air boundary
        #define COMPUTE_DR_B(FIELD, d1, d2) \
            d1 = (-3.0*FIELD.x[i][j][k] + 4.0*FIELD.x[i+1][j][k] - FIELD.x[i+2][j][k]) * inv_2dr; \
            d2 = (FIELD.x[i][j][k] - 2.0*FIELD.x[i+1][j][k] + FIELD.x[i+2][j][k]) * inv_dr2;

        COMPUTE_DR_B(u, dudr, d2udr2)
        COMPUTE_DR_B(v, dvdr, d2vdr2)
        COMPUTE_DR_B(w, dwdr, d2wdr2)
        COMPUTE_DR_B(t, dtdr, d2tdr2)
        COMPUTE_DR_B(c, dcdr, d2cdr2)
        COMPUTE_DR_B(cloud, dclouddr, d2clouddr2)
        COMPUTE_DR_B(ice, dicedr, d2icedr2)
        COMPUTE_DR_B(gr, dgdr, d2gdr2)
        COMPUTE_DR_B(co2, dcodr, d2codr2)
        dpdr = (-3.0*p_dyn.x[i][j][k] + 4.0*p_dyn.x[i+1][j][k] - p_dyn.x[i+2][j][k]) * inv_2dr;
        r_flag = true;
        #undef COMPUTE_DR_B
    }
    else if(i == 0){
        #define COMPUTE_DR_F(FIELD, d1, d2) \
            d1 = (-3.0*FIELD.x[0][j][k] + 4.0*FIELD.x[1][j][k] - FIELD.x[2][j][k]) * inv_2dr; \
            d2 = (FIELD.x[0][j][k] - 2.0*FIELD.x[1][j][k] + FIELD.x[2][j][k]) * inv_dr2;

        COMPUTE_DR_F(u, dudr, d2udr2)
        COMPUTE_DR_F(v, dvdr, d2vdr2)
        COMPUTE_DR_F(w, dwdr, d2wdr2)
        COMPUTE_DR_F(t, dtdr, d2tdr2)
        COMPUTE_DR_F(c, dcdr, d2cdr2)
        COMPUTE_DR_F(cloud, dclouddr, d2clouddr2)
        COMPUTE_DR_F(ice, dicedr, d2icedr2)
        COMPUTE_DR_F(gr, dgdr, d2gdr2)
        COMPUTE_DR_F(co2, dcodr, d2codr2)
        dpdr = (-3.0*p_dyn.x[0][j][k] + 4.0*p_dyn.x[1][j][k] - p_dyn.x[2][j][k]) * inv_2dr;
        r_flag = true;
        #undef COMPUTE_DR_F
    }

    if(!r_flag){
        // Central differences (default path — most common)
        #define COMPUTE_DR_C(FIELD, d1, d2) \
            d1 = (FIELD.x[i+1][j][k] - FIELD.x[i-1][j][k]) * inv_2dr; \
          d2 = (FIELD.x[i+1][j][k] - 2.0*FIELD.x[i][j][k] + FIELD.x[i-1][j][k]) * inv_dr2;

        COMPUTE_DR_C(u, dudr, d2udr2)
        COMPUTE_DR_C(v, dvdr, d2vdr2)
        COMPUTE_DR_C(w, dwdr, d2wdr2)
        COMPUTE_DR_C(t, dtdr, d2tdr2)
        COMPUTE_DR_C(c, dcdr, d2cdr2)
        COMPUTE_DR_C(cloud, dclouddr, d2clouddr2)
        COMPUTE_DR_C(ice, dicedr, d2icedr2)
        COMPUTE_DR_C(gr, dgdr, d2gdr2)
        COMPUTE_DR_C(co2, dcodr, d2codr2)
        dpdr = (p_dyn.x[i+1][j][k] - p_dyn.x[i-1][j][k]) * inv_2dr;
        #undef COMPUTE_DR_C
    }


    // ===== THETA-direction derivatives =====
    if(j > 2 && j < jm-3){
        if(land_ijk && is_air(h, i, j+1, k) && is_air(h, i, j+2, k)){
            #define COMPUTE_DTHE_B(FIELD, d1, d2) \
                d1 = (-3.0*FIELD.x[i][j][k] + 4.0*FIELD.x[i][j+1][k] - FIELD.x[i][j+2][k]) * inv_2dthe; \
                d2 = (FIELD.x[i][j][k] - 2.0*FIELD.x[i][j+1][k] + FIELD.x[i][j+2][k]) * inv_dthe2;

            COMPUTE_DTHE_B(u, dudthe, d2udthe2) COMPUTE_DTHE_B(v, dvdthe, d2vdthe2)
            COMPUTE_DTHE_B(w, dwdthe, d2wdthe2) COMPUTE_DTHE_B(t, dtdthe, d2tdthe2)
            COMPUTE_DTHE_B(c, dcdthe, d2cdthe2) COMPUTE_DTHE_B(cloud, dclouddthe, d2clouddthe2)
            COMPUTE_DTHE_B(ice, dicedthe, d2icedthe2) COMPUTE_DTHE_B(gr, dgdthe, d2gdthe2)
            COMPUTE_DTHE_B(co2, dcodthe, d2codthe2)
            dpdthe = (-3.0*p_dyn.x[i][j][k] + 4.0*p_dyn.x[i][j+1][k] - p_dyn.x[i][j+2][k]) * inv_2dthe;
            the_flag = true;
            #undef COMPUTE_DTHE_B
        }
        if(!the_flag && land_ijk && is_air(h, i, j-1, k) && is_air(h, i, j-2, k)){
            #define COMPUTE_DTHE_C(FIELD, d1, d2) \
                d1 = -(-3.0*FIELD.x[i][j][k] + 4.0*FIELD.x[i][j-1][k] - FIELD.x[i][j-2][k]) * inv_2dthe; \
                d2 = -(FIELD.x[i][j][k] - 2.0*FIELD.x[i][j-1][k] + FIELD.x[i][j-2][k]) * inv_dthe2;

            COMPUTE_DTHE_C(u, dudthe, d2udthe2) COMPUTE_DTHE_C(v, dvdthe, d2vdthe2)
            COMPUTE_DTHE_C(w, dwdthe, d2wdthe2) COMPUTE_DTHE_C(t, dtdthe, d2tdthe2)
            COMPUTE_DTHE_C(c, dcdthe, d2cdthe2) COMPUTE_DTHE_C(cloud, dclouddthe, d2clouddthe2)
            COMPUTE_DTHE_C(ice, dicedthe, d2icedthe2) COMPUTE_DTHE_C(gr, dgdthe, d2gdthe2)
            COMPUTE_DTHE_C(co2, dcodthe, d2codthe2)
            dpdthe = -(-3.0*p_dyn.x[i][j][k] + 4.0*p_dyn.x[i][j-1][k] - p_dyn.x[i][j-2][k]) * inv_2dthe;
            the_flag = true;
            #undef COMPUTE_DTHE_C
        }
    }
    else if(j == 0){
        #define COMPUTE_DTHE_NP(FIELD, d1, d2) \
            d1 = (-3.0*FIELD.x[i][0][k] + 4.0*FIELD.x[i][1][k] - FIELD.x[i][2][k]) * inv_2dthe; \
            d2 = (FIELD.x[i][0][k] - 2.0*FIELD.x[i][1][k] + FIELD.x[i][2][k]) * inv_dthe2;

        COMPUTE_DTHE_NP(u, dudthe, d2udthe2) COMPUTE_DTHE_NP(v, dvdthe, d2vdthe2)
        COMPUTE_DTHE_NP(w, dwdthe, d2wdthe2) COMPUTE_DTHE_NP(t, dtdthe, d2tdthe2)
        COMPUTE_DTHE_NP(c, dcdthe, d2cdthe2) COMPUTE_DTHE_NP(cloud, dclouddthe, d2clouddthe2)
        COMPUTE_DTHE_NP(ice, dicedthe, d2icedthe2) COMPUTE_DTHE_NP(gr, dgdthe, d2gdthe2)
        COMPUTE_DTHE_NP(co2, dcodthe, d2codthe2)
        dpdthe = (-3.0*p_dyn.x[i][0][k] + 4.0*p_dyn.x[i][1][k] - p_dyn.x[i][2][k]) * inv_2dthe;
        the_flag = true;
        #undef COMPUTE_DTHE_NP
    }
    else if(j == jm-1){
        #define COMPUTE_DTHE_SP(FIELD, d1, d2) \
            d1 = -(-3.0*FIELD.x[i][jm-1][k] + 4.0*FIELD.x[i][jm-2][k] - FIELD.x[i][jm-3][k]) * inv_2dthe; \
            d2 = -(FIELD.x[i][jm-1][k] - 2.0*FIELD.x[i][jm-2][k] + FIELD.x[i][jm-3][k]) * inv_dthe2;

        COMPUTE_DTHE_SP(u, dudthe, d2udthe2) COMPUTE_DTHE_SP(v, dvdthe, d2vdthe2)
        COMPUTE_DTHE_SP(w, dwdthe, d2wdthe2) COMPUTE_DTHE_SP(t, dtdthe, d2tdthe2)
        COMPUTE_DTHE_SP(c, dcdthe, d2cdthe2) COMPUTE_DTHE_SP(cloud, dclouddthe, d2clouddthe2)
        COMPUTE_DTHE_SP(ice, dicedthe, d2icedthe2) COMPUTE_DTHE_SP(gr, dgdthe, d2gdthe2)
        COMPUTE_DTHE_SP(co2, dcodthe, d2codthe2)
        dpdthe = -(-3.0*p_dyn.x[i][jm-1][k] + 4.0*p_dyn.x[i][jm-2][k] - p_dyn.x[i][jm-3][k]) * inv_2dthe;
        the_flag = true;
        #undef COMPUTE_DTHE_SP
    }

    if(!the_flag){
        #define COMPUTE_DTHE_A(FIELD, d1, d2) \
            d1 = (FIELD.x[i][j+1][k] - FIELD.x[i][j-1][k]) * inv_2dthe; \
            d2 = (FIELD.x[i][j+1][k] - 2.0*FIELD.x[i][j][k] + FIELD.x[i][j-1][k]) * inv_dthe2;

        COMPUTE_DTHE_A(u, dudthe, d2udthe2) COMPUTE_DTHE_A(v, dvdthe, d2vdthe2)
        COMPUTE_DTHE_A(w, dwdthe, d2wdthe2) COMPUTE_DTHE_A(t, dtdthe, d2tdthe2)
        COMPUTE_DTHE_A(c, dcdthe, d2cdthe2) COMPUTE_DTHE_A(cloud, dclouddthe, d2clouddthe2)
        COMPUTE_DTHE_A(ice, dicedthe, d2icedthe2) COMPUTE_DTHE_A(gr, dgdthe, d2gdthe2)
        COMPUTE_DTHE_A(co2, dcodthe, d2codthe2)
        dpdthe = (p_dyn.x[i][j+1][k] - p_dyn.x[i][j-1][k]) * inv_2dthe;
        #undef COMPUTE_DTHE_A
    }


    // ===== PHI-direction derivatives =====
    if(k > 2 && k < km-3){
        if(land_ijk && is_air(h, i, j, k+1) && is_air(h, i, j, k+2)){
          #define COMPUTE_DPHI_B(FIELD, d1, d2) \
                d1 = (-3.0*FIELD.x[i][j][k] + 4.0*FIELD.x[i][j][k+1] - FIELD.x[i][j][k+2]) * inv_2dphi; \
                d2 = (FIELD.x[i][j][k] - 2.0*FIELD.x[i][j][k+1] + FIELD.x[i][j][k+2]) * inv_dphi2;

            COMPUTE_DPHI_B(u, dudphi, d2udphi2) COMPUTE_DPHI_B(v, dvdphi, d2vdphi2)
            COMPUTE_DPHI_B(w, dwdphi, d2wdphi2) COMPUTE_DPHI_B(t, dtdphi, d2tdphi2)
            COMPUTE_DPHI_B(c, dcdphi, d2cdphi2) COMPUTE_DPHI_B(cloud, dclouddphi, d2clouddphi2)
            COMPUTE_DPHI_B(ice, dicedphi, d2icedphi2) COMPUTE_DPHI_B(gr, dgdphi, d2gdphi2)
            COMPUTE_DPHI_B(co2, dcodphi, d2codphi2)
            dpdphi = (-3.0*p_dyn.x[i][j][k] + 4.0*p_dyn.x[i][j][k+1] - p_dyn.x[i][j][k+2]) * inv_2dphi;
            phi_flag = true;
            #undef COMPUTE_DPHI_B
        }
        if(!phi_flag && land_ijk && is_air(h, i, j, k-1) && is_air(h, i, j, k-2)){
            #define COMPUTE_DPHI_C(FIELD, d1, d2) \
                d1 = -(-3.0*FIELD.x[i][j][k] + 4.0*FIELD.x[i][j][k-1] - FIELD.x[i][j][k-2]) * inv_2dphi; \
                d2 = -(FIELD.x[i][j][k] - 2.0*FIELD.x[i][j][k-1] + FIELD.x[i][j][k-2]) * inv_dphi2;

            COMPUTE_DPHI_C(u, dudphi, d2udphi2) COMPUTE_DPHI_C(v, dvdphi, d2vdphi2)
            COMPUTE_DPHI_C(w, dwdphi, d2wdphi2) COMPUTE_DPHI_C(t, dtdphi, d2tdphi2)
            COMPUTE_DPHI_C(c, dcdphi, d2cdphi2) COMPUTE_DPHI_C(cloud, dclouddphi, d2clouddphi2)
            COMPUTE_DPHI_C(ice, dicedphi, d2icedphi2) COMPUTE_DPHI_C(gr, dgdphi, d2gdphi2)
            COMPUTE_DPHI_C(co2, dcodphi, d2codphi2)
            dpdphi = -(-3.0*p_dyn.x[i][j][k] + 4.0*p_dyn.x[i][j][k-1] - p_dyn.x[i][j][k-2]) * inv_2dphi;
            phi_flag = true;
            #undef COMPUTE_DPHI_C
        }
    }
    else if(k == 0){
        #define COMPUTE_DPHI_W(FIELD, d1, d2) \
            d1 = (-3.0*FIELD.x[i][j][0] + 4.0*FIELD.x[i][j][1] - FIELD.x[i][j][2]) * inv_2dphi; \
            d2 = (FIELD.x[i][j][km-1] - 2.0*FIELD.x[i][j][km+2] + FIELD.x[i][j][km+3]) * inv_dphi2;

        COMPUTE_DPHI_W(u, dudphi, d2udphi2) COMPUTE_DPHI_W(v, dvdphi, d2vdphi2)
        COMPUTE_DPHI_W(w, dwdphi, d2wdphi2) COMPUTE_DPHI_W(t, dtdphi, d2tdphi2)
        COMPUTE_DPHI_W(c, dcdphi, d2cdphi2) COMPUTE_DPHI_W(cloud, dclouddphi, d2clouddphi2)
        COMPUTE_DPHI_W(ice, dicedphi, d2icedphi2) COMPUTE_DPHI_W(gr, dgdphi, d2gdphi2)
        COMPUTE_DPHI_W(co2, dcodphi, d2codphi2)
        dpdphi = (-3.0*p_dyn.x[i][j][0] + 4.0*p_dyn.x[i][j][1] - p_dyn.x[i][j][2]) * inv_2dphi;
        phi_flag = true;
        #undef COMPUTE_DPHI_W
    }
    else if(k == km-1){
        #define COMPUTE_DPHI_E(FIELD, d1, d2) \
            d1 = -(-3.0*FIELD.x[i][j][km-1] + 4.0*FIELD.x[i][j][km-2] - FIELD.x[i][j][km-3]) * inv_2dphi; \
            d2 = -(FIELD.x[i][j][km-1] - 2.0*FIELD.x[i][j][km-2] + FIELD.x[i][j][km-3]) * inv_dphi2;

        COMPUTE_DPHI_E(u, dudphi, d2udphi2) COMPUTE_DPHI_E(v, dvdphi, d2vdphi2)
        COMPUTE_DPHI_E(w, dwdphi, d2wdphi2) COMPUTE_DPHI_E(t, dtdphi, d2tdphi2)
        COMPUTE_DPHI_E(c, dcdphi, d2cdphi2) COMPUTE_DPHI_E(cloud, dclouddphi, d2clouddphi2)
        COMPUTE_DPHI_E(ice, dicedphi, d2icedphi2) COMPUTE_DPHI_E(gr, dgdphi, d2gdphi2)
        COMPUTE_DPHI_E(co2, dcodphi, d2codphi2)
        dpdphi = -(-3.0*p_dyn.x[i][j][km-1] + 4.0*p_dyn.x[i][j][km-2] - p_dyn.x[i][j][km-3]) * inv_2dphi;
        phi_flag = true;
        #undef COMPUTE_DPHI_E
    }

    if(!phi_flag){
        #define COMPUTE_DPHI_A(FIELD, d1, d2) \
            d1 = (FIELD.x[i][j][k+1] - FIELD.x[i][j][k-1]) * inv_2dphi; \
            d2 = (FIELD.x[i][j][k+1] - 2.0*FIELD.x[i][j][k] + FIELD.x[i][j][k-1]) * inv_dphi2;

        COMPUTE_DPHI_A(u, dudphi, d2udphi2) COMPUTE_DPHI_A(v, dvdphi, d2vdphi2)
        COMPUTE_DPHI_A(w, dwdphi, d2wdphi2) COMPUTE_DPHI_A(t, dtdphi, d2tdphi2)
        COMPUTE_DPHI_A(c, dcdphi, d2cdphi2) COMPUTE_DPHI_A(cloud, dclouddphi, d2clouddphi2)
        COMPUTE_DPHI_A(ice, dicedphi, d2icedphi2) COMPUTE_DPHI_A(gr, dgdphi, d2gdphi2)
        COMPUTE_DPHI_A(co2, dcodphi, d2codphi2)
        dpdphi = (p_dyn.x[i][j][k+1] - p_dyn.x[i][j][k-1]) * inv_2dphi;
        #undef COMPUTE_DPHI_A
    }


    // ===== Physics: Coriolis, centrifugal, coefficients =====
    // Coriolis acceleration -2 Ω × v in (r, θ, φ) with θ = colatitude (the0 = 0
    // at N pole), v south-positive (= +ê_θ), w east-positive. TRADITIONAL
    // APPROXIMATION: keep only the locally-vertical rotation Ω_r = Ω cosθ; drop
    // the non-traditional ("Eötvös") terms from the horizontal Ω_θ = -Ω sinθ:
    //   F_r = 0                          (was +2Ω sinθ w — dropped)
    //   F_θ = +2Ω cos(θ) w
    //   F_φ = -2Ω cos(θ) v               (the -2Ω sinθ u part dropped)
    // The kept sign on F_φ's cosθ v term still matters (it gave the N/S symmetry
    // and removed the iter-234 blow-up). The non-traditional terms are tiny in
    // reality but, in this non-hydrostatic explicit solver with a periodically-
    // fired pressure projection, +2Ω sinθ w drove a persistent unbalanced upward
    // force on the eastward jets → the jet pattern slowly advected up the column.
    double two_omega_Latm_dt_over_u0 = 2.0 * omega * L_atm / u_0 * dt;
    double Coriolis_rad =  0.0;
    double Coriolis_the =  two_omega_Latm_dt_over_u0 * costhe * w_ijk;
    double Coriolis_phi = -two_omega_Latm_dt_over_u0 * costhe * v_ijk;

    double rad_dist = (double)i * L_atm * exp_rm;
    double rad_Earth_m = rad_dist + r_Earth * 1e3;
    double dt_over_u0_sq = (dt / u_0) * (dt / u_0);
    double omega2 = omega * omega;
    double centrifugal_rad = omega2 * rad_Earth_m * dt_over_u0_sq;
    double centrifugal_the = omega2 * rad_Earth_m * fabs(sinthe) * dt_over_u0_sq;

    double coeff_energy   = L_atm / (u_0 * cp_l * t_0);
    double coeff_energy_p = u_0 * u_0 / (cp_l * t_0);
    double coeff_trans    = L_atm / (c_0 * u_0);

    double coeff_MC_t     = L_atm / (u_0 * t_0);                        // 0.18304   s/K
    double coeff_MC_q     = L_atm / (u_0 * c_0);                        // 1428.57   s/(kg/kg)
    double coeff_MC_vel   = L_atm / (u_0 * u_0);                        // 6.25      s/(m/s)

    // Inviscid spin-up: diffusion_ramp = 0 during the Euler phase, ramped to 1 afterwards.
    // Multiplying every diff_*_re coefficient by it preserves the RK4 structure — only the
    // Laplacian contributions in the rhs assembly are scaled.
    double diff_t_re   = diffusion_ramp * 1.0 / (re * pr);
    double diff_vel_re = diffusion_ramp * 1.0 / re;
    double diff_prec_re_inv = diffusion_ramp * 1.0 / (sc_WaterVapour * re);
    double diff_co2_re_inv  = diffusion_ramp * 1.0 / (sc_CO2 * re);

    // ===== Transport terms (advection) =====
    // Precompute velocity * metric factors
    double u_exp   = u_ijk * exp_rm;
    double v_invrm = v_ijk * inv_rm;
    double w_invrs = w_ijk * inv_rmsinthe;

    double pressure_t = coeff_energy_p * (u_exp * dpdr + v_invrm * dpdthe + w_invrs * dpdphi);

    // Positivity-preserving (minmod-limited) advective gradients for the moisture
    // scalars — see RHS_Atm_Turb.cpp for the rationale. Limits only the interior
    // central stencils; boundary/pole/seam one-sided stencils keep the centered value.
    const double inv_dr   = 2.0 * inv_2dr;
    const double inv_dthe = 2.0 * inv_2dthe;
    const double inv_dphi = 2.0 * inv_2dphi;

    double dcdr_adv = dcdr, dclouddr_adv = dclouddr, dicedr_adv = dicedr;
    double dcdthe_adv = dcdthe, dclouddthe_adv = dclouddthe, dicedthe_adv = dicedthe;
    double dcdphi_adv = dcdphi, dclouddphi_adv = dclouddphi, dicedphi_adv = dicedphi;
    // Temperature gets the same minmod-limited advective gradient as the moisture
    // scalars — see RHS_Atm_Turb.cpp: a coastal rhs_t budget showed the growing cold
    // bubble is driven by centered -transport_t as a 2Δt checkerboard. diffusion_t keeps
    // its centered stencil.
    double dtdr_adv = dtdr, dtdthe_adv = dtdthe, dtdphi_adv = dtdphi;
    // Velocity gets the same minmod-limited advective gradient (transport_{u,v,w} only) —
    // see RHS_Atm_Turb.cpp: a momentum budget at the Greenwich-seam blow-up cell showed a
    // 2Δ grid-scale velocity checkerboard driven by centered zonal self-advection -w·∂w/∂φ.
    double dudr_adv = dudr, dudthe_adv = dudthe, dudphi_adv = dudphi;
    double dvdr_adv = dvdr, dvdthe_adv = dvdthe, dvdphi_adv = dvdphi;
    double dwdr_adv = dwdr, dwdthe_adv = dwdthe, dwdphi_adv = dwdphi;

    if(!r_flag){
        dtdr_adv     = minmod(t.x[i][j][k]     - t.x[i-1][j][k],     t.x[i+1][j][k]     - t.x[i][j][k])     * inv_dr;
        dudr_adv     = minmod(u.x[i][j][k]     - u.x[i-1][j][k],     u.x[i+1][j][k]     - u.x[i][j][k])     * inv_dr;
        dvdr_adv     = minmod(v.x[i][j][k]     - v.x[i-1][j][k],     v.x[i+1][j][k]     - v.x[i][j][k])     * inv_dr;
        dwdr_adv     = minmod(w.x[i][j][k]     - w.x[i-1][j][k],     w.x[i+1][j][k]     - w.x[i][j][k])     * inv_dr;
        dcdr_adv     = minmod(c.x[i][j][k]     - c.x[i-1][j][k],     c.x[i+1][j][k]     - c.x[i][j][k])     * inv_dr;
        dclouddr_adv = minmod(cloud.x[i][j][k] - cloud.x[i-1][j][k], cloud.x[i+1][j][k] - cloud.x[i][j][k]) * inv_dr;
        dicedr_adv   = minmod(ice.x[i][j][k]   - ice.x[i-1][j][k],   ice.x[i+1][j][k]   - ice.x[i][j][k])   * inv_dr;
    }
    if(!the_flag){
        dtdthe_adv     = minmod(t.x[i][j][k]     - t.x[i][j-1][k],     t.x[i][j+1][k]     - t.x[i][j][k])     * inv_dthe;
        dudthe_adv     = minmod(u.x[i][j][k]     - u.x[i][j-1][k],     u.x[i][j+1][k]     - u.x[i][j][k])     * inv_dthe;
        dvdthe_adv     = minmod(v.x[i][j][k]     - v.x[i][j-1][k],     v.x[i][j+1][k]     - v.x[i][j][k])     * inv_dthe;
        dwdthe_adv     = minmod(w.x[i][j][k]     - w.x[i][j-1][k],     w.x[i][j+1][k]     - w.x[i][j][k])     * inv_dthe;
        dcdthe_adv     = minmod(c.x[i][j][k]     - c.x[i][j-1][k],     c.x[i][j+1][k]     - c.x[i][j][k])     * inv_dthe;
        dclouddthe_adv = minmod(cloud.x[i][j][k] - cloud.x[i][j-1][k], cloud.x[i][j+1][k] - cloud.x[i][j][k]) * inv_dthe;
        dicedthe_adv   = minmod(ice.x[i][j][k]   - ice.x[i][j-1][k],   ice.x[i][j+1][k]   - ice.x[i][j][k])   * inv_dthe;
    }
    if(!phi_flag){
        dtdphi_adv     = minmod(t.x[i][j][k]     - t.x[i][j][k-1],     t.x[i][j][k+1]     - t.x[i][j][k])     * inv_dphi;
        dudphi_adv     = minmod(u.x[i][j][k]     - u.x[i][j][k-1],     u.x[i][j][k+1]     - u.x[i][j][k])     * inv_dphi;
        dvdphi_adv     = minmod(v.x[i][j][k]     - v.x[i][j][k-1],     v.x[i][j][k+1]     - v.x[i][j][k])     * inv_dphi;
        dwdphi_adv     = minmod(w.x[i][j][k]     - w.x[i][j][k-1],     w.x[i][j][k+1]     - w.x[i][j][k])     * inv_dphi;
        dcdphi_adv     = minmod(c.x[i][j][k]     - c.x[i][j][k-1],     c.x[i][j][k+1]     - c.x[i][j][k])     * inv_dphi;
        dclouddphi_adv = minmod(cloud.x[i][j][k] - cloud.x[i][j][k-1], cloud.x[i][j][k+1] - cloud.x[i][j][k]) * inv_dphi;
        dicedphi_adv   = minmod(ice.x[i][j][k]   - ice.x[i][j][k-1],   ice.x[i][j][k+1]   - ice.x[i][j][k])   * inv_dphi;
    }

    double transport_u     = u_exp * dudr_adv     + v_invrm * dudthe_adv     + w_invrs * dudphi_adv;
    double transport_v     = u_exp * dvdr_adv     + v_invrm * dvdthe_adv     + w_invrs * dvdphi_adv;
    double transport_w     = u_exp * dwdr_adv     + v_invrm * dwdthe_adv     + w_invrs * dwdphi_adv;
    double transport_t     = u_exp * dtdr_adv     + v_invrm * dtdthe_adv     + w_invrs * dtdphi_adv;
    double transport_c     = u_exp * dcdr_adv     + v_invrm * dcdthe_adv     + w_invrs * dcdphi_adv;
    double transport_cloud = u_exp * dclouddr_adv + v_invrm * dclouddthe_adv + w_invrs * dclouddphi_adv;
    double transport_ice   = u_exp * dicedr_adv   + v_invrm * dicedthe_adv   + w_invrs * dicedphi_adv;
    double transport_g     = u_exp * dgdr     + v_invrm * dgdthe     + w_invrs * dgdphi;
    double transport_co2   = u_exp * dcodr    + v_invrm * dcodthe    + w_invrs * dcodphi;

    // ===== Diffusion terms =====
    // Precompute shared sub-expressions
    double two_over_rm_exp = 2.0 * inv_rm * exp_rm;
    double cos_rm2sin = costhe_inv_rm2sinthe;
    // Curvature term on v_θ / v_φ in the spherical vector Laplacian: −v/(r²sin²θ),
    // i.e. v_metric = (1/sin²θ)/r² = (1 + cot²θ)/r² = (1 + cos²θ/sin²θ)/r².  The
    // cosine MUST be squared — with a bare costhe the bracket goes negative where
    // cosθ < −sin²θ (poleward of ~38° in the SOUTHERN hemisphere, where costhe<0),
    // turning −v_metric·v from damping into anti-diffusion → exponential growth of
    // v,w. That was the persistent southern-hemisphere (≈68°S) near-surface seed;
    // the sign flips only in the south (costhe<0), explaining the N/S asymmetry.
    double v_metric = (1.0 + costhe * costhe / sinthe2) * inv_rm2;

    double diffusion_t = (d2tdr2 * exp_2_rm + dtdr * two_over_rm_exp
        + d2tdthe2 * inv_rm2 + dtdthe * cos_rm2sin
        + d2tdphi2 * inv_rm2sinthe2) * diff_t_re;

    double diffusion_u = (d2udr2 * exp_2_rm + 2.0 * u_ijk * inv_rm2
        + d2udthe2 * inv_rm2 + 4.0 * dudr * inv_rm * exp_rm
        + dudthe * cos_rm2sin + d2udphi2 * inv_rm2sinthe2) * diff_vel_re;

    double diffusion_v = (d2vdr2 * exp_2_rm + dvdr * two_over_rm_exp
        + d2vdthe2 * inv_rm2 + dvdthe * cos_rm2sin
        - v_metric * v_ijk + d2vdphi2 * inv_rm2sinthe2
        + 2.0 * dudthe * inv_rm2
        - dwdphi * 2.0 * costhe * inv_rm2sinthe2) * diff_vel_re;

    double diffusion_w = (d2wdr2 * exp_2_rm + dwdr * two_over_rm_exp
        + d2wdthe2 * inv_rm2 + dwdthe * cos_rm2sin
        - v_metric * w_ijk + d2wdphi2 * inv_rm2sinthe2
        + 2.0 * dudphi * inv_rm2sinthe
        + dvdphi * 2.0 * costhe * inv_rm2sinthe2) * diff_vel_re;

    double diffusion_c = (d2cdr2 * exp_2_rm + dcdr * two_over_rm_exp
        + d2cdthe2 * inv_rm2 + dcdthe * cos_rm2sin
        + d2cdphi2 * inv_rm2sinthe2) * diff_prec_re_inv;

    double diffusion_cloud = (d2clouddr2 * exp_2_rm + dclouddr * two_over_rm_exp
        + d2clouddthe2 * inv_rm2 + dclouddthe * cos_rm2sin
        + d2clouddphi2 * inv_rm2sinthe2) * diff_prec_re_inv;

    double diffusion_ice = (d2icedr2 * exp_2_rm + dicedr * two_over_rm_exp
        + d2icedthe2 * inv_rm2 + dicedthe * cos_rm2sin
        + d2icedphi2 * inv_rm2sinthe2) * diff_prec_re_inv;

    double diffusion_g = (d2gdr2 * exp_2_rm + dgdr * two_over_rm_exp  // NOTE: original had dicedr here — likely a bug
        + d2gdthe2 * inv_rm2 + dgdthe * cos_rm2sin
        + d2gdphi2 * inv_rm2sinthe2) * diff_prec_re_inv;

    double diffusion_co2 = (d2codr2 * exp_2_rm + dcodr * two_over_rm_exp
        + d2codthe2 * inv_rm2 + dcodthe * cos_rm2sin
        + d2codphi2 * inv_rm2sinthe2) * diff_co2_re_inv;


    // ===== RHS assembly =====
    double dpdr_exp = dpdr * exp_rm;
    double dpdthe_invrm = dpdthe * inv_rm;
    double dpdphi_invrs = dpdphi * inv_rmsinthe;

    rhs_t.x[i][j][k] = pressure_t - transport_t + diffusion_t
        + coeff_energy * lv * (S_c.x[i][j][k] + S_r.x[i][j][k])
        + coeff_energy * ls * (S_i.x[i][j][k] + S_s.x[i][j][k] + S_g.x[i][j][k])
        + coeff_MC_t * MC_t.x[i][j][k];

    // ----- Held-Suarez Newtonian thermal relaxation (PROTOTYPE) -----
    // The general circulation has no maintained energy source in ATOM: rhs_t had no
    // radiative forcing, so the imposed APE is consumed and the flow spins down (a
    // forced-equilibrium IC like the Grotjahn mean cells cannot survive). Add a Newtonian
    // relaxation that nudges T back toward the radiative-equilibrium target t_eq (the
    // snapshot of the initial Scotese+topography+land-sea field), continuously regenerating
    // APE. Held-Suarez (1994) rate: k_T = k_a + (k_s-k_a)*max(0,(sigma-sigma_b)/(1-sigma_b))
    // *cos^4(lat), with sigma = p/p_surf, cos(lat)=sin(colatitude)=sinthe. Faster near the
    // surface (k_s, ~4 day) than in the free atmosphere (k_a, ~40 day). Non-dimensionalised
    // exactly like surf_drag/Coriolis (k * L_atm/u_0 * dt).
    constexpr bool   HELD_SUAREZ_RELAX = true;        // prototype toggle
    if (HELD_SUAREZ_RELAX) {
        constexpr double k_a     = 1.0 / ( 4.0 * 86400.0);   // free-atmosphere relax rate [1/s] (10x H-S: matches ATOM's 10x friction scaling)
        constexpr double k_s     = 1.0 / ( 0.4 * 86400.0);   // near-surface relax rate   [1/s] (10x H-S)
        constexpr double sigma_b = 0.7;                       // boundary-layer top in sigma
        const double p_surf = p_stat.x[i_topography[j][k]][j][k];
        const double sigma  = (p_surf > 0.0) ? p_stat.x[i][j][k] / p_surf : 0.0;
        double bl = (sigma - sigma_b) / (1.0 - sigma_b);     // 0 at sigma_b -> 1 at surface
        if (bl < 0.0) bl = 0.0;
        const double cos4lat = sinthe2 * sinthe2;             // cos^4(lat) = sin^4(colatitude)
        const double k_T = k_a + (k_s - k_a) * bl * cos4lat;  // [1/s]
        // HS-dt TEST (2026-07-04): removed the extra *dt. Advective-time nondim of a Newtonian
        // relaxation is k_T*L_atm/u_0 (dimensionless), integrated by RK4's own *dt like every
        // other rhs_t term (advection/diffusion/latent, all dt^1). The extra *dt made HS enter
        // as dt^2 = ~1e-4 too weak -> inert (t_eq shift never reached T; CO2 forcing invisible).
        rhs_t.x[i][j][k] -= (k_T * L_atm / u_0) * (t.x[i][j][k] - t_eq.x[i][j][k]);
    }

    // ----- Near-surface Rayleigh (boundary-layer) drag on the horizontal wind -----
    // The free-slip wall (bcSolidGround) and the init-only bcVelSurfSur leave the
    // near-surface tangential wind with NO momentum sink, and molecular diffusion is
    // ~1e-6·v at re≈1000. Without surface friction the coastal-jet v,w accelerate
    // unopposed → the i=1 CFL blow-up that recurs at steep coasts (Antarctic
    // Peninsula / Baffin / Norway / Sea of Okhotsk). Add Held-Suarez-style linear drag
    // on v,w (NOT radial u): full strength at the first air cell above the LOCAL
    // surface (i_topography[j][k]), ramping linearly to zero over the boundary layer.
    // Smooth volumetric sink — no coastal discontinuity (unlike bcVelSurfSur's
    // 0.01/0.9 sawtooth). Non-dimensionalised exactly like the Coriolis term.
    constexpr double rayleigh_kf   = 1.0 / 86400.0;   // [EXPERIMENT 2026-06-16: cut 10x to ~1/day H-S to relieve momentum over-damping (T->wind decoupling); orig 1/8640]  surface drag rate ≈ Ekman strength [1/s], 10× the old ≈1/day Held-Suarez baseline. Tuned 2026-06-13: baseline 1/day gave ~34 m/s eastward w off W-coast S-America; 10× cut the surface to ~28 m/s and 20× gained nothing (the jet max sits at ~1 km, above the drag layer), so 10× is the settled strength.
    constexpr double drag_n_layers = 5.0;             // boundary-layer depth [air cells]. Tried 10 (2026-06-13) to reach the ~1 km coastal-jet max — no effect (that ~27 m/s max at 27°S/71°W is an Andes orographic/pressure-gradient feature, immune to friction), so kept at the physical 5.
    double drag_profile = 1.0 - (double)(i - i_topography[j][k]) / drag_n_layers;
    if(drag_profile < 0.0) drag_profile = 0.0;
    if(drag_profile > 1.0) drag_profile = 1.0;
    double surf_drag = rayleigh_kf * L_atm / u_0 * dt * drag_profile;

    // Boussinesq buoyancy in perturbation-pressure form: the radial body force is the
    // density-anomaly contribution, +g·(t − t̄(i)), referenced to the per-level horizontal
    // mean temperature t_ref_level[i] (NOT a global t_0=273.15 K). The hydrostatic part is
    // absorbed into p_stat; referencing to the level mean makes the body force zero in the
    // horizontal mean at every height, so it no longer fights the lapse-rate p_stat/p_dyn
    // split — only horizontal thermal contrasts drive vertical motion. Earlier (t − 1) left
    // a large standing column force (≈ −0.08 at 5.5 km) that loaded the pressure solver.
    rhs_u.x[i][j][k] = -dpdr_exp - transport_u + diffusion_u
        + buoyancy_ramp * buoyancy * g * dt / u_0 * (t.x[i][j][k] - t_ref_level[i])
        + Coriolis * Coriolis_rad - centrifugal * centrifugal_rad;

    rhs_v.x[i][j][k] = -dpdthe_invrm - transport_v + diffusion_v
        + Coriolis * Coriolis_the - centrifugal * centrifugal_the
        + coeff_MC_vel * MC_v.x[i][j][k]
        - surf_drag * v_ijk;

    rhs_w.x[i][j][k] = -dpdphi_invrs - transport_w + diffusion_w
        + Coriolis * Coriolis_phi
        + coeff_MC_vel * MC_w.x[i][j][k]
        - surf_drag * w_ijk;

    rhs_c.x[i][j][k] = -transport_c + diffusion_c
        + coeff_trans * S_v.x[i][j][k] 
        + coeff_MC_q * MC_q.x[i][j][k];

    rhs_cloud.x[i][j][k] = -transport_cloud + diffusion_cloud
        + coeff_trans * S_c.x[i][j][k];

    rhs_ice.x[i][j][k] = -transport_ice + diffusion_ice
        + coeff_trans * S_i.x[i][j][k];

    rhs_g.x[i][j][k] = -transport_g + diffusion_g
        + coeff_trans * S_g.x[i][j][k];

    rhs_co2.x[i][j][k] = -transport_co2 + diffusion_co2;

    aux_u.x[i][j][k] = rhs_u.x[i][j][k] + dpdr_exp;
    aux_v.x[i][j][k] = rhs_v.x[i][j][k] + dpdthe_invrm;
    aux_w.x[i][j][k] = rhs_w.x[i][j][k] + dpdphi_invrs;
 }
