/*
 * Ocean General Circulation Modell(OGCM) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to combine the right hand sides of the differential equations for the Runge-Kutta scheme
*/

#include "cHydrosphereModel.h"
#include <iostream>
#include <cmath>
#include <vector>
#include "Utils.h"

using namespace std;
using namespace AtomUtils;




void cHydrosphereModel::RHS_Hydrosphere(int i, int j, int k, const CellGeometry& geo){

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

    const double inv_2dthe = geo.inv_2dthe;
    const double inv_2dphi = geo.inv_2dphi;
    const double inv_dthe2 = geo.inv_dthe2;
    const double inv_dphi2 = geo.inv_dphi2;
    // radial derivatives use the per-i non-uniform stencil coeffs (rc*/rf*), not inv_2dr/inv_dr2

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
    double dudr, dvdr, dwdr, dtdr, dpdr, dcdr;
    double dudthe, dvdthe, dwdthe, dtdthe, dpdthe, dcdthe;
    double dudphi, dvdphi, dwdphi, dtdphi, dpdphi, dcdphi;

    // ---- Second-order derivatives storage ----
    double d2udr2, d2vdr2, d2wdr2, d2tdr2, d2cdr2;
    double d2udthe2, d2vdthe2, d2wdthe2, d2tdthe2, d2cdthe2;
    double d2udphi2, d2vdphi2, d2wdphi2, d2tdphi2, d2cdphi2;


    // ===== R-direction derivatives =====
    if(i < im-2 && land_ijk && is_air(h, i+1, j, k)){
        // Forward-biased stencil near land-air boundary
        // Forward-biased stencil (points i, i+1, i+2) on the stretched radial grid
        #define COMPUTE_DR_B(FIELD, d1, d2) \
            d1 = rf10[i]*FIELD.x[i][j][k] + rf11[i]*FIELD.x[i+1][j][k] + rf12[i]*FIELD.x[i+2][j][k]; \
            d2 = rf20[i]*FIELD.x[i][j][k] + rf21[i]*FIELD.x[i+1][j][k] + rf22[i]*FIELD.x[i+2][j][k];

        COMPUTE_DR_B(u, dudr, d2udr2)
        COMPUTE_DR_B(v, dvdr, d2vdr2)
        COMPUTE_DR_B(w, dwdr, d2wdr2)
        COMPUTE_DR_B(t, dtdr, d2tdr2)
        COMPUTE_DR_B(c, dcdr, d2cdr2)
        dpdr = rf10[i]*p_dyn.x[i][j][k] + rf11[i]*p_dyn.x[i+1][j][k] + rf12[i]*p_dyn.x[i+2][j][k];
        r_flag = true;
        #undef COMPUTE_DR_B
    }
    else if(i == 0){
        // Forward stencil at the bottom (points 0, 1, 2) on the stretched radial grid
        #define COMPUTE_DR_F(FIELD, d1, d2) \
            d1 = rf10[0]*FIELD.x[0][j][k] + rf11[0]*FIELD.x[1][j][k] + rf12[0]*FIELD.x[2][j][k]; \
            d2 = rf20[0]*FIELD.x[0][j][k] + rf21[0]*FIELD.x[1][j][k] + rf22[0]*FIELD.x[2][j][k];

        COMPUTE_DR_F(u, dudr, d2udr2)
        COMPUTE_DR_F(v, dvdr, d2vdr2)
        COMPUTE_DR_F(w, dwdr, d2wdr2)
        COMPUTE_DR_F(t, dtdr, d2tdr2)
        COMPUTE_DR_F(c, dcdr, d2cdr2)
        dpdr = rf10[0]*p_dyn.x[0][j][k] + rf11[0]*p_dyn.x[1][j][k] + rf12[0]*p_dyn.x[2][j][k];
        r_flag = true;
        #undef COMPUTE_DR_F
    }

    if(!r_flag){
        // Central differences (default path — most common)
        // Central stencil (points i-1, i, i+1) on the stretched radial grid
        #define COMPUTE_DR_C(FIELD, d1, d2) \
            d1 = rc1m[i]*FIELD.x[i-1][j][k] + rc10[i]*FIELD.x[i][j][k] + rc1p[i]*FIELD.x[i+1][j][k]; \
            d2 = rc2m[i]*FIELD.x[i-1][j][k] + rc20[i]*FIELD.x[i][j][k] + rc2p[i]*FIELD.x[i+1][j][k];

        COMPUTE_DR_C(u, dudr, d2udr2)
        COMPUTE_DR_C(v, dvdr, d2vdr2)
        COMPUTE_DR_C(w, dwdr, d2wdr2)
        COMPUTE_DR_C(t, dtdr, d2tdr2)
        COMPUTE_DR_C(c, dcdr, d2cdr2)
        dpdr = rc1m[i]*p_dyn.x[i-1][j][k] + rc10[i]*p_dyn.x[i][j][k] + rc1p[i]*p_dyn.x[i+1][j][k];
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
            COMPUTE_DTHE_C(c, dcdthe, d2cdthe2) 
            
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
        COMPUTE_DTHE_NP(c, dcdthe, d2cdthe2) 
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
        COMPUTE_DTHE_SP(c, dcdthe, d2cdthe2) 
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
        COMPUTE_DTHE_A(c, dcdthe, d2cdthe2) 
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
            COMPUTE_DPHI_B(c, dcdphi, d2cdphi2) 
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
            COMPUTE_DPHI_C(c, dcdphi, d2cdphi2) 
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
        COMPUTE_DPHI_W(c, dcdphi, d2cdphi2) 
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
        COMPUTE_DPHI_E(c, dcdphi, d2cdphi2) 
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
        COMPUTE_DPHI_A(c, dcdphi, d2cdphi2) 
        dpdphi = (p_dyn.x[i][j][k+1] - p_dyn.x[i][j][k-1]) * inv_2dphi;
        #undef COMPUTE_DPHI_A
    }


    // ===== Physics: Coriolis, centrifugal, coefficients =====
    // Coriolis acceleration -2 Ω × v in (r, θ, φ) with θ = colatitude
    // (the0 = 0 at N pole), v south-positive (= +ê_θ), w east-positive:
    //   F_r = +2Ω sin(θ) w               (Eötvös: east → up at equator)
    //   F_θ = +2Ω cos(θ) w
    //   F_φ = -2Ω cos(θ) v - 2Ω sin(θ) u
    // EXPERIMENT 2026-06-26 (uncommitted): the inverse Rossby number in the
    // advective-time nondim is 2*Omega*L/u_0 with NO dt (the RK4 supplies the
    // timestep). The previous coefficient carried an extra *dt, leaving Coriolis
    // ~1/dt (~2000x) too weak vs advection/diffusion (w-budget confirmed) -> no
    // geostrophic gyre balance -> SH spin-down. Test: drop the *dt.
    // See project_hydro_coriolis_dt_scaling.
    double two_omega_Lhyd_over_u0 = 2.0 * omega * L_hyd / u_0;
    double Coriolis_rad =  two_omega_Lhyd_over_u0 * sinthe * w_ijk;
    double Coriolis_the =  two_omega_Lhyd_over_u0 * costhe * w_ijk;
    double Coriolis_phi = -two_omega_Lhyd_over_u0 * (costhe * v_ijk + sinthe * u_ijk);

    double rad_dist = (double)i * L_hyd * exp_rm;
    double rad_Earth_m = rad_dist + r_Earth * 1e3;
    double dt_over_u0_sq = (dt / u_0) * (dt / u_0);
    double omega2 = omega * omega;
    double centrifugal_rad = omega2 * rad_Earth_m * dt_over_u0_sq;
    double centrifugal_the = omega2 * rad_Earth_m * fabs(sinthe) * dt_over_u0_sq;

//    double coeff_energy   = L_hyd / (u_0 * cp_l * t_0);
    double coeff_energy_p = u_0 * u_0 / (cp_w * t_0);

    // Inviscid spin-up: diffusion_ramp = 0 during the Euler phase, ramped to 1 afterwards.
    double diff_t_re   = diffusion_ramp * 1.0 / (re * pr);
    double diff_vel_re = diffusion_ramp * 1.0 / re;
    double diff_prec_re_inv = diffusion_ramp * 1.0 / (sc * re);

    // ===== Transport terms (advection) =====
    // Precompute velocity * metric factors
    double u_exp   = u_ijk * exp_rm;
    double v_invrm = v_ijk * inv_rm;
    double w_invrs = w_ijk * inv_rmsinthe;

    double pressure_t = coeff_energy_p * (u_exp * dpdr + v_invrm * dpdthe + w_invrs * dpdphi);

    double transport_t     = u_exp * dtdr     + v_invrm * dtdthe     + w_invrs * dtdphi;
    double transport_u     = u_exp * dudr     + v_invrm * dudthe     + w_invrs * dudphi;
    double transport_v     = u_exp * dvdr     + v_invrm * dvdthe     + w_invrs * dvdphi;
    double transport_w     = u_exp * dwdr     + v_invrm * dwdthe     + w_invrs * dwdphi;
    double transport_c     = u_exp * dcdr     + v_invrm * dcdthe     + w_invrs * dcdphi;

    // ===== Diffusion terms =====
    // Precompute shared sub-expressions
    double two_over_rm_exp = 2.0 * inv_rm * exp_rm;
    double cos_rm2sin = costhe_inv_rm2sinthe;
    // Curvature term on v_θ / v_φ in the spherical vector Laplacian: −v/(r²sin²θ),
    // i.e. (1 + cos²θ/sin²θ)/r². The cosine MUST be squared — a bare costhe makes
    // the bracket negative poleward of ~38°S (costhe<0), turning −v_metric·v into
    // anti-diffusion. Same bug/fix as atmosphere RHS_Atm.cpp; mirror convention.
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


    // ===== RHS assembly =====
    double dpdr_exp = dpdr * exp_rm;
    double dpdthe_invrm = dpdthe * inv_rm;
    double dpdphi_invrs = dpdphi * inv_rmsinthe;

    rhs_t.x[i][j][k] = pressure_t - transport_t + diffusion_t;

    // Boussinesq buoyancy with thermal + haline anomaly.  Only the anomaly
    // relative to the reference state (t = 1, c = 1) drives radial motion;
    // the hydrostatic reference is absorbed into the pressure split.
    //
    // Linearised EOS:  ρ/ρ_0 = 1 - α_T·T_0·(t-1) + β_S·S_0·(c-1)
    // Buoyancy:        -g·(ρ-ρ_0)/ρ_0 = +g·α_T·T_0·(t-1) - g·β_S·S_0·(c-1)
    //
    // Warm  (t > 1) → less dense → rises  (positive contribution).
    // Salty (c > 1) → denser     → sinks (negative contribution).
    //
    // alpha_S ≈ (β_S·S_0)/(α_T·T_0) ≈ (7.6e-4 × 35)/(2e-4 × 283) ≈ 0.5
    // sets the relative weighting of haline vs. thermal forcing.  Adjust if
    // the model uses different α_T, β_S, or reference T_0/S_0.
    constexpr double alpha_S = 0.5;
    rhs_u.x[i][j][k] = -dpdr_exp - transport_u + diffusion_u
        + buoyancy * g * dt / u_0
            * ((t.x[i][j][k] - 1.0) - alpha_S * (c.x[i][j][k] - 1.0))
        + Coriolis * Coriolis_rad - centrifugal * centrifugal_rad;

    rhs_v.x[i][j][k] = -dpdthe_invrm - transport_v + diffusion_v
        + Coriolis * Coriolis_the - centrifugal * centrifugal_the;

    rhs_w.x[i][j][k] = -dpdphi_invrs - transport_w + diffusion_w
        + Coriolis * Coriolis_phi;

    // w-momentum-budget capture (checkpoint iters only): store the four rhs_w
    // contributions so write_w_momentum_budget can attribute the zonal-mean
    // zonal-velocity tendency by latitude/depth. The last RK4 sub-stage written
    // wins, which is representative for the diagnostic.
    if (wbudget_capture) {
        wbud_pgf.x[i][j][k]  = -dpdphi_invrs;
        wbud_cor.x[i][j][k]  =  Coriolis * Coriolis_phi;
        wbud_adv.x[i][j][k]  = -transport_w;
        wbud_diff.x[i][j][k] =  diffusion_w;
    }

    // Surface wind-stress forcing. The ocean momentum equation otherwise has NO
    // sustained wind input — the atmosphere wind only INITIALISES the currents
    // (EkmanSpiral, once), so every wind-driven current (West Wind Drift / ACC,
    // subtropical gyres) spins down (Coriolis+diffusion erode it, nothing drives
    // it; see project_hydro_no_wind_stress_forcing). Apply the bulk wind stress
    // tau = r_air*C_D*|U_wind|*U_wind as a body force in the topmost integrated
    // cell (i=im-2, whose thickness ~ the mixed layer; the surface row im-1 is a
    // BC extrapolated from here). Non-dim by L_hyd/u_0^2 (advective time, matching
    // PGF/advection/diffusion). v_wind/w_wind [m/s] = atmosphere surface wind,
    // set once in EkmanSpiral and persistent (resetArrays runs before it).
    constexpr bool WIND_STRESS_FORCING = true;
    if (WIND_STRESS_FORCING && i == im - 2 && is_water(h, i, j, k)) {
        constexpr double C_D = 2.6e-3;                                    // bulk drag coefficient
        const double dz_top = (rad.z[im-1] - rad.z[im-2]) * L_hyd;        // [m] top-cell (mixed-layer) thickness
        const double U_wind = sqrt(v_wind.y[j][k] * v_wind.y[j][k]
                                 + w_wind.y[j][k] * w_wind.y[j][k]);      // [m/s]
        const double accel_nd = r_air * C_D * U_wind / (r_0_water * dz_top)
                              * L_hyd / (u_0 * u_0);                      // nondim accel per (m/s) of wind component
        rhs_w.x[i][j][k] += accel_nd * w_wind.y[j][k];                    // zonal stress      (East+)
        rhs_v.x[i][j][k] += accel_nd * v_wind.y[j][k];                    // meridional stress (South+)
    }

    rhs_c.x[i][j][k] = -transport_c + diffusion_c;


    aux_u.x[i][j][k] = rhs_u.x[i][j][k] + dpdr_exp;
    aux_v.x[i][j][k] = rhs_v.x[i][j][k] + dpdthe_invrm;
    aux_w.x[i][j][k] = rhs_w.x[i][j][k] + dpdphi_invrs;
 }
