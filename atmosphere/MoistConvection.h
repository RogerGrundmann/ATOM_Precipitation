#pragma once

#include "cAtmosphereModel.h"

#include <vector>
#include <cmath>
#include <chrono>
#include <iostream>
#include <cstdio>
#include <algorithm>

using namespace AtomUtils;



namespace AtomMoistConvection {
    constexpr double a_ev = 1.0e-3;
    constexpr double b_ev = 5.9;
    constexpr double t_00 = 236.15;
    constexpr double t_00_minus2 = t_00 - 2.0;
    constexpr double bet = 42.0;
    constexpr double b_u = 0.33;
    constexpr double alf_1 = 0.05;                                      // evaporation rate coefficient [1/s]
    constexpr double alf_2 = 0.05;                                      // evaporation rate coefficient [1/s]
    constexpr double bet_p = 2.0e-3;                                    // in [1/s]
    constexpr double R_cloud = 100.0;                                   // cloud radius in [m] by ECMWF
    double eps_u = 0.2 / R_cloud;                                       // in [1/m] by ECMWF
    double del_u = eps_u;                                               // in [1/m] by ECMWF
//    double eps_u = 1.0e-4;                                              // in [1/m] by COSMO for penetrative (deep) and midlevel clouds
//    double del_u = eps_u;                                               // in [1/m] by COSMO for penetrative (deep) and midlevel clouds
//    constexpr double eps_u = 1.0e-4;                                    // in [1/m] by COSMO for shallow clouds
//    constexpr double del_u = eps_u;                                     // in [1/m] by COSMO for shallow clouds
    constexpr double eps_d = 1.0e-4;                                    // turbulent mixing by COSMO only applied in lower levels for deep clouds
    constexpr double del_d = eps_d;                                     // turbulent mixing by COSMO only applied in lower levels for deep clouds
    constexpr double gam_d = -0.3;                                      // by COSMO
    constexpr double a_u = 3.0e-2;                                      // [·/·]
    constexpr double a_d = 3.0e-2;                                      // [·/·]
    constexpr double vel = 1.1;                                         // [·/·]
    constexpr double p_stat_beg = 1000.0;                               // in [hPa]
    constexpr double p_stat_end = 970.0;                                // in [hPa]
    constexpr double p_stat_Cloud_Base = 900.0;                         // in [hPa]
    constexpr double p_stat_base = 800.0;                               // in [hPa]
    constexpr double p_stat_diff = 200.0;                               // in [hPa]
    constexpr double p_stat_midlevel = 700.0;                           // cloud-base pressure threshold for midlevel convection [hPa] (Bechtold 2001)
    constexpr double t_add_u = 0.2;                                     // in [K]
    constexpr double q_v_u_add = 1.0e-4;                                // in [kg/kg]
    constexpr double coeff_recurr = 0.1; 
    constexpr double scale = 0.98;                                      // this value creates deep cloud initial values q_c_u because of t_add_u and t_add_u
}
/*
*
*/
// ===================== CALCULATION OF DEEP MOIST CONVECTION =======================================
class MoistConvection {
public:
    explicit MoistConvection(cAtmosphereModel& model)
        : m(model)
        , cloud(model.cloud)
        , R_W_R_A(model.R_WaterVapour / model.R_Air)
        , eps_thermo(1.0 / R_W_R_A)
        , alf((1.0 - eps_thermo) / eps_thermo)
    {}

    int iter;

    void run(int iter) {
        using namespace std;
        using namespace AtomMoistConvection;

        cout << endl << endl << "      MoistConvection  class starts" << endl;
        auto begin = std::chrono::high_resolution_clock::now();

        iter = m.iter_n;

        precompute();
        findCloudBaseLFS();
        surfaceFluxPerturbations();
        initUpdraft();
        CloudBaseUpdate();


        for(int iter_prec = 1; iter_prec <= m.iter_prec; iter_prec++){
            updraftEntrainment(iter_prec);
            CloudBaseUpdate();
            updraftRecurrence();

            downdraftEntrainment(iter_prec);
            downdraftRecurrence();
        }

        computeCAPE();
        rhsForcing();
        subTerrainFill();



        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for MoistConvection  class ended\n", elapsed.count() * 1e-9);
    }

private:
    cAtmosphereModel& m;
    const Array& cloud;
    const double R_W_R_A;
    const double eps_thermo;
    const double alf;


    // Precomputed data
    std::vector<double> height_table;
    std::vector<double> step;
    std::vector<int8_t> land_surf;
    std::vector<int8_t> air_surf;

    // Column indices (local to this run)
    std::vector<std::vector<int>> i_Base_local;
    std::vector<std::vector<int>> i_LFS_local;
    std::vector<std::vector<int>> i_deep_beg_local;
    std::vector<std::vector<int>> i_deep_end_local;
    std::vector<std::vector<double>> M_d_LFS_local;
    std::vector<std::vector<double>> delta_T_sfp;    // surface-flux T perturbation per column [K]      (Bechtold 2008)
    std::vector<std::vector<double>> delta_q_sfp;    // surface-flux q perturbation per column [kg/kg]  (Bechtold 2008)


// ==================== NaN-SAFE HELPERS ====================
    // Saturation specific humidity: q_sat = ep * E_sat / (p - E_sat).
    // Guards against p <= E_sat (fully saturated / super-saturated column).
    static double safe_q_sat(double ep, double E_sat, double p_u) noexcept {
        return ep * E_sat / std::max(p_u - E_sat, 1e-10);
    }

    // Moist-air density floor: prevents 1/r_humid blow-up when r_humid is
    // clamped to 0 by the limiter in UtilsAtm.
    static double safe_r_humid(double r_h) noexcept {
        return std::max(r_h, 1e-6);
    }

// ==================== ENTRAINMENT GRADIENTS (shared helper) ====================
    // Returns the moisture gradients in r, theta and phi directions plus c_ijk = c.x[i][j][k] (floored at 1e-4).
    struct EntrainmentGradients { double dcdr, dcdthe, dcdphi, c_ijk; };  // c_ijk = c.x[i][j][k], floored at 1e-4

    // Neumann (zero-normal-gradient) extrapolation of c into a terrain cell,
    // using the same second-order formula as bcScalarSurfSur:
    //   c_land = (4/3)*c_near_air - (1/3)*c_far_air
    // Called when one neighbor of the centered-difference stencil is terrain.
    // 'near' is the current (air) cell value; 'far' is the cell on the opposite
    // side of the stencil (which may itself be terrain — fall back to 'near' then).
    static inline double neumannC(double c_near, double c_far) noexcept {
        return (4.0/3.0) * c_near - (1.0/3.0) * c_far;
    }

    inline EntrainmentGradients entrainmentGradients(int i, int j, int k) {
        const double costhe       = cos(m.the.z[90] - m.the.z[j]);
        const double height       = m.get_layer_height(i);              // [m]
        const double r_the_height = (m.r_Earth * 1.0e3 + height) * costhe;
        const double U_Earth      = r_the_height * 2.0 * M_PI;          // [m]
        const double step_dr      = step[i];
        const double inv_step_dr  = 0.5 / step_dr;
        const double inv_step_the = 0.5 / (U_Earth / (m.jm-1));         // [1/m]
        const double inv_step_phi = 1.0 / (U_Earth / (m.km-1));
        const int    i_p          = std::min(i+1, m.im-1);
        const int    i_m          = std::max(i-1, 0);

        // Surface-gradient BC (bcScalarSurfSur approach): when a stencil neighbor
        // is terrain, replace its c value with a second-order Neumann extrapolation
        // from the current (air) cell and the cell on the opposite side, so the
        // normal gradient at the land face evaluates to the air-side gradient
        // rather than a spurious c_ocean − c_land=0 difference.
        // If both neighbors are terrain the component is set to zero.
        const double c0 = m.c.x[i][j][k];

        // ---- radial (i) ----
        const bool land_ip = is_land(m.h, i_p, j, k);
        const bool land_im = is_land(m.h, i_m, j, k);
        const double c_ip = land_ip
            ? (land_im ? c0 : neumannC(c0, m.c.x[i_m][j][k]))
            : m.c.x[i_p][j][k];
        const double c_im = land_im
            ? (land_ip ? c0 : neumannC(c0, m.c.x[i_p][j][k]))
            : m.c.x[i_m][j][k];
        const double dcdr = (c_ip - c_im) * inv_step_dr;

        // ---- meridional (j / theta) ----
        const bool land_jp = is_land(m.h, i, j+1, k);
        const bool land_jm = is_land(m.h, i, j-1, k);
        const double c_jp = land_jp
            ? (land_jm ? c0 : neumannC(c0, m.c.x[i][j-1][k]))
            : m.c.x[i][j+1][k];
        const double c_jm = land_jm
            ? (land_jp ? c0 : neumannC(c0, m.c.x[i][j+1][k]))
            : m.c.x[i][j-1][k];
        const double dcdthe = (c_jp - c_jm) * inv_step_the;

        // ---- zonal (k / phi) ----
        const bool land_kp = is_land(m.h, i, j, k+1);
        const bool land_km = is_land(m.h, i, j, k-1);
        const double c_kp = land_kp
            ? (land_km ? c0 : neumannC(c0, m.c.x[i][j][k-1]))
            : m.c.x[i][j][k+1];
        const double c_km = land_km
            ? (land_kp ? c0 : neumannC(c0, m.c.x[i][j][k+1]))
            : m.c.x[i][j][k-1];
        const double dcdphi = (c_kp - c_km) * inv_step_phi;

        return {
            dcdr,
            dcdthe,
            dcdphi,
            std::max(c0, 1e-4)   // floor: was exact ==0.0 check, misses small positives
        };
    }

// ==================== PRECOMPUTE ====================
    void precompute() {
        height_table.resize(m.im);
        step.assign(m.im, 0.0);
        for(int i = 0; i < m.im; i++){
            height_table[i] = m.get_layer_height(i);
        }
        for(int i = 0; i < m.im - 1; i++){
            step[i] = height_table[i+1] - height_table[i];
        }
        if(m.im >= 2) step[m.im-1] = step[m.im-2];

        land_surf.resize(m.jm * m.km);
        air_surf.resize(m.jm * m.km);

        #pragma omp parallel for collapse(2)
        for(int j = 0; j < m.jm; j++){
            for(int k = 0; k < m.km; k++){
                int idx = j * m.km + k;
                land_surf[idx] = is_land(m.h, 0, j, k) ? 1 : 0;
                air_surf[idx]  = is_air(m.h, 0, j, k)  ? 1 : 0;
            }
        }

        i_Base_local.assign(m.jm, std::vector<int>(m.km, 0));
        i_LFS_local.assign(m.jm, std::vector<int>(m.km, 0));
        i_deep_beg_local.assign(m.jm, std::vector<int>(m.km, -1));
        i_deep_end_local.assign(m.jm, std::vector<int>(m.km, -1));
        M_d_LFS_local.assign(m.jm, std::vector<double>(m.km, 0));
        delta_T_sfp.assign(m.jm, std::vector<double>(m.km, 0.0));
        delta_q_sfp.assign(m.jm, std::vector<double>(m.km, 0.0));

        m.K_u = std::vector<double>(m.im, 0.0);
        m.K_d = std::vector<double>(m.im, 0.0);
        m.CAPE = std::vector<double>(m.im, 0.0);

        // Reset 3-D convection fields so stale values from a previous call
        // cannot corrupt M_u/M_d when the saturation
        // condition is no longer met in cells that were active before.
        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                for (int i = 0; i < m.im; i++) {
                    // Scalar fields driven by s = cp_l*T/s_0.
                    // s must be computed from current t before any other function reads it;
                    // if left stale, downdraftRecurrence propagates the error downward.
                    m.s.x[i][j][k]      = m.cp_l * m.t.x[i][j][k] * m.t_0 / m.s_0;

                    // Mass-flux fields — reset so cells outside the active [i_base..i_lfs]
                    // range cannot contribute to rhsForcing flux divergences.
                    m.M_u.x[i][j][k]    = 0.0;
                    m.M_d.x[i][j][k]    = 0.0;

                    // Updraft scalars — stale values corrupt MC_t / MC_q (causes #5, #2, #3).
                    m.s_u.x[i][j][k]    = 0.0;
                    m.q_v_u.x[i][j][k]  = 0.0;
                    m.q_c_u.x[i][j][k]  = 0.0;

                    // Source terms used in rhsForcing conv_src:
                    // (L/cp)*c_u*t_0 ≈ 68 K/s per 1e-4 kg/kg/s; stale c_u drives 1e7 K spikes.
                    m.c_u.x[i][j][k]    = 0.0;
                    m.e_l.x[i][j][k]    = 0.0;
                    m.e_d.x[i][j][k]    = 0.0;
                    m.e_p.x[i][j][k]    = 0.0;
                    m.g_p.x[i][j][k]    = 0.0;

                    // Downdraft scalars — stale s_d / q_v_d corrupt the downdraft recurrence
                    // and the flux divergence in rhsForcing.
                    m.s_d.x[i][j][k]    = 0.0;
                    m.q_v_d.x[i][j][k]  = 0.0;

                    // Convective precipitation — stale P_conv propagates through sigma_p
                    // and e_p into MC_q/PrecipitationConv.
                    m.P_conv.x[i][j][k] = 0.0;
                }
            }
        }

        m.rad.Coordinates(m.im, m.r0, m.dr);
    }
/*
*
*/
// ==================== SURFACE-FLUX PERTURBATIONS (BECHTOLD 2008) ====================
    // Computes per-column T and q perturbations for shallow convection based on surface fluxes.
    // Populates delta_T_sfp / delta_q_sfp for every column, controlled by
    // m.convection_perturbation (set in param.py / config XML):
    //
    //   0 – fixed offsets: δT = t_add_u, δq = q_v_u_add  for all columns
    //   1 – Bechtold (2008) for shallow columns (cloud depth < p_stat_diff):
    //           δT = α · H_s / (ρ · c_p · w*)
    //           δq = α · E  / (ρ · w*)
    //       where w* = (g/T · z_BL · H_s/(ρ·c_p))^(1/3)  [Deardorff 1970]
    //       and z_BL ≈ cloud-base height as proxy for BL top.
    //       Deep columns and shallow columns with H_s ≤ 0 receive zero perturbation,
    //       suppressing deep convection triggering.
    //
    // Reference: Bechtold et al. (2008) J. Atmos. Sci. 65:3077, Sect. 2b.
    void surfaceFluxPerturbations() {
        using namespace AtomMoistConvection;
        constexpr double alpha_sfp = 0.2;   // perturbation amplitude (Bechtold 2008, ECMWF)

        #pragma omp parallel for collapse(2) schedule(static)
        for(int j = 1; j < m.jm-1; j++){
            for(int k = 1; k < m.km-1; k++){

                if(m.convection_perturbation == 0){
                    // switch 0: fixed offsets for all columns (original behaviour)
                    delta_T_sfp[j][k] = t_add_u;
                    delta_q_sfp[j][k] = q_v_u_add;
                    continue;
                }

                // switch 1: zero by default → deep columns get no perturbation,
                // suppressing deep-convection triggering
                delta_T_sfp[j][k] = 0.0;
                delta_q_sfp[j][k] = 0.0;

                int i_base = i_Base_local[j][k];
                int i_lfs  = i_LFS_local[j][k];
                if(i_base < 0 || i_lfs <= i_base) continue;

                // Midlevel convection (Bechtold 2001): triggered by free tropospheric
                // instability, not surface fluxes → always use fixed offsets,
                // independent of convection_mode.  convection_mode controls whether
                // midlevel clouds are non-precipitating (updraftEntrainment), not
                // how they are triggered.
                if(m.p_stat.x[i_base][j][k] < p_stat_midlevel){
                    delta_T_sfp[j][k] = t_add_u;
                    delta_q_sfp[j][k] = q_v_u_add;
                    continue;
                }

                // Deep convection: leave at zero (Bechtold 2008: α = 0 for deep)
                if((m.p_stat.x[i_base][j][k] - m.p_stat.x[i_lfs][j][k]) >= p_stat_diff)
                    continue;

                // Shallow convection branch: zero if surface flux is non-positive
                double H_s = m.Q_sensible_2D.y[j][k];   // [W/m²]
                if(H_s <= 0.0) continue;

                // Moisture flux derived from latent heat flux
                double E = m.Q_latent_2D.y[j][k] / m.lv;  // [kg/(m²s)]

                // Surface state (i = 0 is the lowest model level)
                double T_sfc   = m.t.x[0][j][k] * m.t_0;             // [K]
                double rho_sfc = safe_r_humid(m.r_humid.x[0][j][k]); // [kg/m³]

                // BL depth: cloud-base height as proxy for BL top
                double z_BL = height_table[i_base];     // [m]
                if(z_BL <= 0.0) continue;

                // Deardorff (1970) convective velocity scale: w* = (g/T · z_BL · H_s/(ρ·c_p))^(1/3)
                double B_sfc  = (m.g / T_sfc) * H_s / (rho_sfc * m.cp_l);  // [m²/s³]
                double w_star = std::cbrt(B_sfc * z_BL);                    // [m/s]
                if(w_star <= 0.0) continue;

                double inv_rho_w  = 1.0 / (rho_sfc * w_star);
                // Cap at fixed offsets: Bechtold gives smaller values for strong BL
                // (w* ~ 2 m/s → δq ~ 5e-6 kg/kg); the cap prevents blow-up when
                // w* is near zero but positive (weak-flux edge case).
                delta_T_sfp[j][k] = std::min(alpha_sfp * H_s * inv_rho_w / m.cp_l, t_add_u);
                delta_q_sfp[j][k] = std::min(alpha_sfp * E  * inv_rho_w, q_v_u_add);
            }
        }
    }
/*
*
*/
// ==================== UPDRAFT INITIAL VALUES ====================
    void initUpdraft() {
        using namespace AtomMoistConvection;

        #pragma omp parallel for collapse(2)
        for(int j = 1; j < m.jm-1; j++){
            for(int k = 1; k < m.km-1; k++){

                int local_i_Base = i_Base_local[j][k];
                int local_i_beg = i_deep_beg_local[j][k];
                int local_i_end = i_deep_end_local[j][k];
                if(local_i_beg == -1 || local_i_end == -1) continue;

                int local_i_LFS = i_LFS_local[j][k];

                // Perturbations set by surfaceFluxPerturbations() for all columns
                const double t_pert = delta_T_sfp[j][k];
                const double q_pert = delta_q_sfp[j][k];

                m.M_u.x[local_i_beg][j][k] = 
                    m.r_humid.x[local_i_beg][j][k] * m.u.x[local_i_beg][j][k];

                for(int i = local_i_beg+1; i <= local_i_end; i++){

                    double t_u     = m.t.x[i][j][k] * m.t_0;
                    double t_u_add = t_u + t_pert;
                    double r_h_i   = m.r_humid.x[i][j][k];
                    double p_u     = m.p_stat.x[i][j][k];

                    double E_sat_add = m.hp * AtomUtils::exp_func(t_u_add, 17.2694, 35.86);
                    double q_sat_add = safe_q_sat(m.ep, E_sat_add, p_u);


                    m.q_v_u.x[i][j][k]     = m.c.x[i][j][k] + q_pert;
                    if(m.q_v_u.x[i][j][k] >= scale * q_sat_add && t_u >= t_00)
                        m.q_c_u.x[i][j][k] = m.q_v_u.x[i][j][k] - scale * q_sat_add;
                    if(m.q_v_u.x[i][j][k] <= 0.0) m.q_v_u.x[i][j][k] = 0.0;
                    if(m.q_c_u.x[i][j][k] <= 0.0) m.q_c_u.x[i][j][k] = 0.0;


                    // Entrainment
                    auto [dcdr, dcdthe, dcdphi, c_ijk] = entrainmentGradients(i, j, k);

                    if((m.p_stat.x[local_i_Base][j][k]
                        - m.p_stat.x[local_i_LFS][j][k]) >= p_stat_diff) {
                        m.E_u.x[i][j][k] =                              // local moisture convergence by entrainment
                            - r_h_i / c_ijk
                            * (m.u.x[i][j][k] * dcdr
                            +  m.v.x[i][j][k] * dcdthe
                            +  m.w.x[i][j][k] * dcdphi);
                    } else  m.E_u.x[i][j][k] = 0.0;

                    if(t_u <= t_00_minus2) m.E_u.x[i][j][k] = 0.0;


                    // Detrainment
                    m.D_u.x[i][j][k] = 0.0;                             // no detrainment


                    double d_Mu      = m.E_u.x[i][j][k] - m.D_u.x[i][j][k];// in [kg/(m³s)]   updraft mass flux based on the local moisture convergence
                    m.M_u.x[i][j][k] = m.M_u.x[i-1][j][k] + d_Mu * step[i];// in [kg/(m³s)]
                    if(is_land(m.h, i, j, k)) m.M_u.x[i][j][k] = 0.0;

                    m.s.x[i][j][k]   = m.cp_l * t_u / m.s_0;
                    m.s_u.x[i][j][k] = m.cp_l * t_u_add / m.s_0;

                    m.u_u.x[i][j][k] = vel * m.u.x[i][j][k];
                    m.v_u.x[i][j][k] = vel * m.v.x[i][j][k];
                    m.w_u.x[i][j][k] = vel * m.w.x[i][j][k];

                }

                M_d_LFS_local[j][k] = gam_d * m.M_u.x[local_i_beg][j][k];// downdraft mass flux based on the mass flux at cloud base
                m.M_d.x[local_i_LFS][j][k] = M_d_LFS_local[j][k];


                for(int i = local_i_end+1; i <= local_i_LFS-1; i++){

                    double inv_a_u = 1.0 / (a_u * m.u_0);

                    double t_u     = m.t.x[i][j][k] * m.t_0;
                    double r_h_i   = m.r_humid.x[i][j][k];
                    double t_u_add = t_u + t_pert;
                    double p_u     = m.p_stat.x[i][j][k];

                    double E_sat_add = m.hp * AtomUtils::exp_func(t_u_add, 17.2694, 35.86);
                    double q_sat_add = safe_q_sat(m.ep, E_sat_add, p_u);

                    if(i == local_i_end+1)  m.M_u.x[i-1][j][k] = m.M_u.x[local_i_beg][j][k];

                    m.s.x[i][j][k]   = m.cp_l * t_u / m.s_0;
                    m.s_u.x[i][j][k] = m.cp_l * t_u / m.s_0;

                    m.q_v_u.x[i][j][k]     = m.c.x[i][j][k] + q_pert;
                    if(m.q_v_u.x[i][j][k] >= scale * q_sat_add && t_u_add >= t_00)
                        m.q_c_u.x[i][j][k] = m.q_v_u.x[i][j][k] - scale * q_sat_add;
//                    if(m.q_v_u.x[i][j][k] <= 0.0) m.q_v_u.x[i][j][k] = 0.0;
//                    if(m.q_c_u.x[i][j][k] <= 0.0) m.q_c_u.x[i][j][k] = 0.0;


                    // Entrainment
                    auto [dcdr, dcdthe, dcdphi, c_ijk] = entrainmentGradients(i, j, k);

                    if((m.p_stat.x[local_i_Base][j][k]
                        - m.p_stat.x[local_i_LFS][j][k]) >= p_stat_diff) {
                        m.E_u.x[i][j][k] =                              // local moisture convergence by entrainment
                            - r_h_i / c_ijk
                            * (m.u.x[i][j][k] * dcdr
                            +  m.v.x[i][j][k] * dcdthe
                            +  m.w.x[i][j][k] * dcdphi);
                    } else  m.E_u.x[i][j][k] = 0.0;

                    if(t_u <= t_00_minus2) m.E_u.x[i][j][k] = 0.0;


                    // Detrainment
                    m.D_u.x[i][j][k] = 0.0;                             // no detrainment


                    double d_Mu = m.E_u.x[i][j][k] - m.D_u.x[i][j][k];  // in [kg/(m³s)]   updraft mass flux based on the local moisture convergence
                    m.M_u.x[i][j][k] = m.M_u.x[i-1][j][k] + d_Mu * step[i];  // in [kg/(m²s)]
                    if(is_land(m.h, i, j, k)) m.M_u.x[i][j][k] = 0.0;

                    m.u_u.x[i][j][k] = m.u.x[i][j][k] + m.M_u.x[i][j][k] * inv_a_u / safe_r_humid(r_h_i);
                    m.v_u.x[i][j][k] = m.v.x[i][j][k];
                    m.w_u.x[i][j][k] = m.w.x[i][j][k];
                }
            }
        }
    }
/*
*
*/
// ==================== FIND CLOUD BASE & LFS, DEEP BEG/END ====================
void findCloudBaseLFS() {                                                                                                                                                      
    using namespace AtomMoistConvection;                                                                                                                                                      
                                                                                                                                                                                              
    // ---- Phase 1: compute raw indices (parallel) ----                                                                                                                                    
    #pragma omp parallel
    {
        std::vector<double> q_sat_col(m.im);

        #pragma omp for collapse(2) schedule(dynamic, 4)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {

                for (int i = 0; i < m.im; i++) {
                    if (m.p_stat.x[i][j][k] <= p_stat_beg) {
                        i_deep_beg_local[j][k] = i;
                        m.Deep_beg.x[i][j][k]  = height_table[i];
                        break;
                    }
                }

                for (int i = m.im - 1; i >= 0; i--) {
                    if (m.p_stat.x[i][j][k] >= p_stat_end) {
                        i_deep_end_local[j][k] = i;
                        m.Deep_end.x[i][j][k]  = height_table[i];
                        break;
                    }
                }



                for (int i = i_deep_beg_local[j][k]; i < m.im; i++) {
                    const double t_u     = m.t.x[i][j][k] * m.t_0;
                    const double t_u_add = t_u + t_add_u;
                    const double p_u     = m.p_stat.x[i][j][k];

                    double E_sat_add = m.hp * AtomUtils::exp_func(t_u_add, 17.2694, 35.86);

                    q_sat_col[i] = safe_q_sat(m.ep, E_sat_add, p_u);
                }

                for (int i = i_deep_beg_local[j][k]; i < m.im; i++) {
                    const double t_u            = m.t.x[i][j][k] * m.t_0;
                    const double cloud_t_weight = 1.0 / (1.0 + std::exp(-(t_u - t_00) / 2.0));

                    m.q_v_u.x[i][j][k] = m.c.x[i][j][k] + q_v_u_add;
                    if(m.q_v_u.x[i][j][k] >= scale * q_sat_col[i] && t_u >= t_00)
                        m.q_c_u.x[i][j][k] = m.q_v_u.x[i][j][k] - scale * q_sat_col[i];
                    if(m.q_v_u.x[i][j][k] <= 0.0) m.q_v_u.x[i][j][k] = 0.0;
                    if(m.q_c_u.x[i][j][k] <= 0.0) m.q_c_u.x[i][j][k] = 0.0;

                    if (m.q_v_u.x[i][j][k] >= scale * q_sat_col[i] && cloud_t_weight > 0.01
                            && m.p_stat.x[i][j][k] <= p_stat_Cloud_Base) {  // cloud base must be above 900 hPa
                        i_Base_local[j][k]     = i;
                        m.CloudBase.x[i][j][k] = height_table[i] * cloud_t_weight;
                        break;
                    }
                }

                for (int i = m.im - 1; i >= 0; i--) {
                    const double q_buoy_lfs =
                        0.5 * (cloud.x[i][j][k] + scale * q_sat_col[i]) - m.c.x[i][j][k];

                    if (q_buoy_lfs <= 0.0) {
                        i_LFS_local[j][k] = std::min(i, m.im-2);        // guard: i+1 must stay in bounds
                        m.LevelFreeSinking.x[i][j][k] = height_table[i];
                        break;
                    }
                }
            }
        }
    }

    for(int j = 0; j < m.jm; j++){
        for(int k = 0; k < m.km; k++){
            for(int i = 0; i < m.im; i++){
                if(is_land(m.h, i, j, k)){
                    m.Deep_beg.x[i][j][k]         = 0.0;
                    m.Deep_end.x[i][j][k]         = 0.0;
                    m.CloudBase.x[i][j][k]        = 0.0;
                    m.LevelFreeSinking.x[i][j][k] = 0.0;
                }
            }
        }
    }
}
/*
*
*/
// ==================== UPDRAFT ENTRAINMENT ====================
    void updraftEntrainment(int iter_prec) {
        using namespace AtomMoistConvection;

        #pragma omp parallel for collapse(2) schedule(static)
        for(int j = 1; j < m.jm-1; j++){
            for(int k = 1; k < m.km-1; k++){

                int i_base = i_Base_local[j][k];
                int i_lfs  = i_LFS_local[j][k];

                int idx_surf = j * m.km + k;
                bool is_land_surf_jk = land_surf[idx_surf];

                if(i_base < 0 || i_lfs <= i_base) continue;

                // Shallow and midlevel convection are non-precipitating:
                // set delta_i_c beyond any reachable cloud depth so K_p stays zero.
                // Midlevel: cloud base above p_stat_midlevel (700 hPa), Bechtold (2001).
                bool is_shallow   = m.convection_mode >= 1
                                    && (m.p_stat.x[i_base][j][k]
                                        - m.p_stat.x[i_lfs][j][k]) < p_stat_diff;
                bool is_midlevel  = m.convection_mode >= 2
                                    && m.p_stat.x[i_base][j][k] < p_stat_midlevel;
                double delta_i_c  = (is_shallow || is_midlevel) ? 1.0e9
//                                  : (is_land_surf_jk ? 3000.0 : 1500.0); // Tiedtke (1989)
                                  : (is_land_surf_jk ? 1000.0 : 500.0); // Bechtold (2001)
                double height_base = height_table[i_base];

                // Cause #1: initUpdraft left non-zero M_u in [0..i_base-1].  Those cells
                // are never touched by this function, so their values create a flux
                // discontinuity at i_base in rhsForcing.  Zero them here on every
                // iter_prec call so the chain below i_base is clean.
                for (int ii = 0; ii < i_base; ii++)
                    m.M_u.x[ii][j][k] = 0.0;

                m.M_u.x[i_base][j][k] =
                    m.r_humid.x[i_base][j][k] * m.u.x[i_base][j][k];    // [kg/(m²s)]
                if(is_land(m.h, i_base, j, k)) m.M_u.x[i_base][j][k] = 0.0;

                // Causes #2 and #3: updraftRecurrence starts at i_base+1 and reads
                // s_u[i_base] and q_v_u[i_base] as i-1 seed values, but no function
                // ever assigns them at exactly i_base.  Set them here.
                {
                    const double t_base = m.t.x[i_base][j][k] * m.t_0;
                    const double t_pert = delta_T_sfp[j][k];
                    const double q_pert = delta_q_sfp[j][k];
                    m.s_u.x[i_base][j][k]   = m.cp_l * (t_base + t_pert) / m.s_0;
                    m.q_v_u.x[i_base][j][k] = std::max(m.c.x[i_base][j][k] + q_pert, 0.0);
                }

                M_d_LFS_local[j][k] = gam_d * m.M_u.x[i_base][j][k];

                for(int i = i_base+1; i <= i_lfs; i++){

                    double height_i = height_table[i];
                        double t_u      = m.t.x[i][j][k] * m.t_0;
                        double r_h_i    = m.r_humid.x[i][j][k];

                        // Entrainment
                        auto [dcdr, dcdthe, dcdphi, c_ijk] = entrainmentGradients(i, j, k);

                        if((m.p_stat.x[i_base][j][k]
                            - m.p_stat.x[i_lfs][j][k]) >= p_stat_diff) {
                            m.E_u.x[i][j][k] =                          // local moisture convergence by entrainment
                                - r_h_i / c_ijk
                                * (m.u.x[i][j][k] * dcdr
                                +  m.v.x[i][j][k] * dcdthe
                                +  m.w.x[i][j][k] * dcdphi)
                                + eps_u * m.M_u.x[i][j][k];
                        } else  m.E_u.x[i][j][k] = 0.0;

                        if(t_u <= t_00_minus2) m.E_u.x[i][j][k] = 0.0;



                        // Detrainment
                        // Turbulent detrainment: applies on all iterations below LFS
                        m.D_u.x[i][j][k] = (i <= i_lfs)
                            ? eps_u * m.M_u.x[i][j][k] : 0.0;           // turbulent detrainment [kg/(m³s)]
                            

                        if(iter_prec > 1 && i == i_lfs)
                            m.D_u.x[i][j][k] = del_u * m.M_u.x[i][j][k] // turbulent detrainment
                                + b_u * m.M_u.x[i][j][k]/step[i];       // dynamic detrainment

                        if(iter_prec > 1 && i == i_lfs-1)
                            m.D_u.x[i][j][k] = del_u * m.M_u.x[i][j][k] // turbulent detrainment
                                + (1.0 - b_u) * m.M_u.x[i][j][k]/step[i];// dynamic detrainment

                        if(t_u <= t_00_minus2) m.D_u.x[i][j][k] = 0.0;


                        double d_Mu = m.E_u.x[i][j][k] - m.D_u.x[i][j][k];// in [(kg/kg)/s]
                        m.M_u.x[i][j][k] = m.M_u.x[i-1][j][k] + d_Mu * step[i];// in [(kg/kg)/s]
                        if(is_land(m.h, i, j, k)) m.M_u.x[i][j][k] = 0.0;


                        // Condensation
                        if(iter_prec == 1)
                            m.c_u.x[i][j][k] = 0.0;
                        else{
                            // Cause #4: dividing by r_humid at high altitude (where
                            // r_humid → 0) amplifies c_u by up to 1/safe_r_humid = 1e6.
                            // Floor at 0.01 kg/m³ (≈ air density at 50 hPa) so the
                            // result stays physical for all tropospheric levels.
                            const double r_h_floor = std::max(r_h_i, 0.01);
                            // Units: q_c_u [kg/kg] * M_u [kg/(m²s)] / r_h [kg/m³] / step [m]
                            //      = (kg/kg)/s  ✓
                            m.c_u.x[i][j][k] = m.q_c_u.x[i][j][k]
                                * m.M_u.x[i][j][k] / (r_h_floor * step[i]);                      // in [(kg/kg)/s]
                            if(m.c_u.x[i][j][k] <= 0.0 || m.q_c_u.x[i][j][k] <= 0.0)
                                m.c_u.x[i][j][k] = 0.0;
                        }

                        // Evaporation of detrained cloud water in updraft
                        if(iter_prec == 1)
                            m.e_l.x[i][j][k] = 0.0;
                        else{
                            const double r_h_floor = std::max(r_h_i, 0.01);
                            m.e_l.x[i][j][k] = m.D_u.x[i][j][k] * m.q_c_u.x[i][j][k] / r_h_floor;// in [(kg/kg)/s]
                            if(m.e_l.x[i][j][k] <= 0.0) m.e_l.x[i][j][k] = 0.0;
                        }

                        // Precipitation formation
                        double K_p = (height_i > height_base + delta_i_c) ? bet_p : 0.0;

                        if(t_u >= m.t_0){
                            m.g_p.x[i][j][k] = (iter_prec == 1) ? 0.0 : K_p * m.q_c_u.x[i][j][k];// in [(kg/kg)/s]
                        }else{
                            m.g_p.x[i][j][k] = 0.0;
                        }
                } // end i

                // Zero M_u at the LFS so rhsForcing sees no updraft flux jump at the
                // cloud top: the updraft fully detrains at this level.
                if (i_base >= 0 && i_lfs > i_base)
                    m.M_u.x[i_lfs][j][k] = 0.0;

           } // end k
        } // end j
    }
/*
*
*/
// ==================== DOWNDRAFT ENTRAINMENT ====================
    void downdraftEntrainment(int iter_prec) {
        using namespace AtomMoistConvection;

        #pragma omp parallel for collapse(2) schedule(static)
        for(int j = 0; j < m.jm; j++){
            for(int k = 0; k < m.km; k++){

                int i_base = i_Base_local[j][k];
                int i_lfs  = i_LFS_local[j][k];
                // Cause #7: without this guard, columns with no detected cloud
                // (i_lfs still at its initialised value 0) still set M_d[0],
                // which rhsForcing then applies to every row/column.
                if(i_lfs <= 0 || i_lfs <= i_base) continue;

                double t_u = m.t.x[i_lfs][j][k] * m.t_0;
                double p_u = m.p_stat.x[i_lfs][j][k];

                m.M_d.x[i_lfs][j][k] = M_d_LFS_local[j][k];

                m.E_d.x[i_lfs][j][k] = eps_d * fabs(m.M_d.x[i_lfs][j][k]);
                m.D_d.x[i_lfs][j][k] = del_d * fabs(m.M_d.x[i_lfs][j][k]);

                double E_sat = m.hp * AtomUtils::exp_func(t_u, 17.2694, 35.86);
                double q_sat = safe_q_sat(m.ep, E_sat, p_u);

                m.q_v_d.x[i_lfs][j][k] = 0.5 * (cloud.x[i_lfs][j][k] + scale * q_sat);

                double r_humid_lfs = 1e2 * p_u
                    / (m.R_Air * (1.0 + (R_W_R_A - 1.0) * m.c.x[i_lfs][j][k]) * t_u);
                double inv_a_d = 1.0 / (a_d * m.u_0);

                m.u_d.x[i_lfs][j][k] = m.u.x[i_lfs][j][k]
                    + m.M_d.x[i_lfs][j][k] * inv_a_d / safe_r_humid(r_humid_lfs);

                m.v_d.x[i_lfs][j][k] = m.v.x[i_lfs][j][k];
                m.w_d.x[i_lfs][j][k] = m.w.x[i_lfs][j][k];
                m.s_d.x[i_lfs][j][k] = m.s.x[i_lfs][j][k];

                m.P_conv.x[i_lfs][j][k] = 0.0;

                // Nordeng (1994): convective precipitation area fraction σ = M_u / (ρ·w_u)
                // evaluated at cloud base; w_u is in model units → multiply by u_0 for [m/s]
                double sigma_p = 0.0;
                {
                    double w_u_cb = fabs(m.w_u.x[i_base][j][k]) * m.u_0;
                    if(w_u_cb > 0.0)
                        sigma_p = fabs(m.M_u.x[i_base][j][k])
                                  / (safe_r_humid(m.r_humid.x[i_base][j][k]) * w_u_cb);
                    sigma_p = std::min(sigma_p, 1.0);
                }

                for(int i = i_lfs-1; i >= 0; i--){

                    p_u = m.p_stat.x[i][j][k];
                    t_u = m.t.x[i][j][k] * m.t_0;

                    double r_h_i = m.r_humid.x[i][j][k];

                    m.s_d.x[i][j][k] = m.cp_l * t_u / m.s_0;

                    double E_sat = m.hp * AtomUtils::exp_func(t_u, 17.2694, 35.86);
                    double q_sat = safe_q_sat(m.ep, E_sat, p_u);

                    m.q_v_d.x[i][j][k] = 0.5 * (cloud.x[i][j][k] + scale * q_sat);

                    m.E_d.x[i][j][k] = eps_d * fabs(m.M_d.x[i][j][k]);
                    m.D_d.x[i][j][k] = del_d * fabs(m.M_d.x[i][j][k]);

                    double d_Md      = m.E_d.x[i][j][k] - m.D_d.x[i][j][k];  // in [kg/(m³s)]
                    m.M_d.x[i][j][k] = m.M_d.x[i+1][j][k] - d_Md * step[i];  // in [kg/(m²s)]
                    if(is_land(m.h, i, j, k)) m.M_d.x[i][j][k] = 0.0;

                    double inv_a_d   = 1.0 / (a_d * m.u_0);
                    m.u_d.x[i][j][k] = m.u.x[i][j][k]
                        + m.M_d.x[i][j][k] * inv_a_d / safe_r_humid(r_h_i);
                    m.v_d.x[i][j][k] = m.v.x[i][j][k];
                    m.w_d.x[i][j][k] = m.w.x[i][j][k];

                    // Evaporation below cloud base
                    if(i <= i_base){
                        if(iter_prec == 1)
                          m.e_p.x[i][j][k] = 0.0;
                        else{
                            double precip_i = m.P_conv.x[i][j][k];
                            if(precip_i > 0.0 && sigma_p > 0.0)
                                m.e_p.x[i][j][k] = sigma_p * alf_1      // in [(kg/kg)/s]
                                    * (q_sat - m.q_v_d.x[i][j][k])
                                    * sqrt(1e1 / alf_2 * precip_i / sigma_p);
                            if(m.e_p.x[i][j][k] <= 0.0) m.e_p.x[i][j][k] = 0.0;
                        }
                    }

                    // Evaporation within downdraft
                    if(iter_prec == 1)
                        m.e_d.x[i][j][k] = 0.0;
                    else{
                        double precip_i = m.P_conv.x[i][j][k];
                        if(precip_i > 0.0){
                          double P_16 = pow(precip_i, 1.0/6.0);
                          double P_49 = pow(precip_i, 4.0/9.0);
                          m.e_d.x[i][j][k] = a_ev * (1.0 + b_ev * P_16) // in [(kg/kg)/s]
                              * (q_sat - m.q_v_d.x[i][j][k]) * P_49;
                        }else{
                          m.e_d.x[i][j][k] = 0.0;
                        }
                        if(m.e_d.x[i][j][k] <= 0.0) m.e_d.x[i][j][k] = 0.0;
                    }

                } // end i downdraft
            }
        }
    }
/*
*
*/
// ==================== UPDRAFT MAIN RECURRENCE ====================
    void updraftRecurrence() {
        using namespace AtomMoistConvection;

        #pragma omp parallel for collapse(2) schedule(static)
        for(int j = 1; j < m.jm-1; j++){
            for(int k = 1; k < m.km-1; k++){

                int i_base = i_Base_local[j][k];
                int i_LFS  = i_LFS_local[j][k];
                if(i_base < 0 || i_LFS <= i_base) continue;

                for(int i = i_base+1; i <= i_LFS-1; i++){

                    double t_u       = m.t.x[i-1][j][k] * m.t_0;
                    double r_h_prev  = m.r_humid.x[i-1][j][k];
                    double step_prev = step[i-1];
                    double r_h_i     = m.r_humid.x[i][j][k];

                    double d_Mu = m.E_u.x[i-1][j][k] - m.D_u.x[i-1][j][k];

                    m.M_u.x[i][j][k] = m.M_u.x[i-1][j][k] + d_Mu * step_prev;
                    if(is_land(m.h, i, j, k)) m.M_u.x[i][j][k] = 0.0;
                    if(t_u <= t_00_minus2) m.M_u.x[i][j][k] = 0.0;

                    double M_u_prev = m.M_u.x[i-1][j][k];

                    double dummy_q_v_u = M_u_prev * m.q_v_u.x[i-1][j][k]
                        + step_prev * (m.E_u.x[i-1][j][k] * m.c.x[i-1][j][k]
                        - m.D_u.x[i-1][j][k] * m.q_v_u.x[i-1][j][k]
                        - r_h_prev * m.c_u.x[i-1][j][k]);

                    // c_u is already included in the q_c_u seed via updraftEntrainment;
                    // adding r_h * c_u here again would double-count condensation and
                    // create a runaway positive feedback (amplification ~step[m] per level).
                    double dummy_q_c_u = M_u_prev * m.q_c_u.x[i-1][j][k]
                        + step_prev * (-m.D_u.x[i-1][j][k] * cloud.x[i-1][j][k]
                        - r_h_prev * m.g_p.x[i-1][j][k]);

                    double dummy_vel_v_u = M_u_prev * m.v_u.x[i-1][j][k]
                        + step_prev * m.E_u.x[i-1][j][k] * m.v.x[i-1][j][k];
                    double dummy_vel_w_u = M_u_prev * m.w_u.x[i-1][j][k]
                        + step_prev * m.E_u.x[i-1][j][k] * m.w.x[i-1][j][k];

                    double L_latent = (t_u >= m.t_0) ? m.lv : m.ls;

                    double dummy_s_u = (M_u_prev * m.s_u.x[i-1][j][k]
                        + step_prev * (m.E_u.x[i-1][j][k] * m.s.x[i-1][j][k]
                        - m.D_u.x[i-1][j][k] * m.s_u.x[i-1][j][k]
                        + L_latent * r_h_prev * m.c_u.x[i-1][j][k])) / m.s_0;

                    double M_u_i = m.M_u.x[i][j][k];
                    if(fabs(M_u_i) > coeff_recurr){
                        double inv_M_u = 1.0 / M_u_i;
                        m.q_v_u.x[i][j][k] = dummy_q_v_u * inv_M_u;
                        m.q_c_u.x[i][j][k] = dummy_q_c_u * inv_M_u;
                        m.v_u.x[i][j][k]   = dummy_vel_v_u * inv_M_u;
                        m.w_u.x[i][j][k]   = dummy_vel_w_u * inv_M_u;
                        m.s_u.x[i][j][k]   = dummy_s_u * inv_M_u;
                    } else {
                        // When M_u ≤ coeff_recurr the updraft is negligible; set
                        // updraft scalars to their environment counterparts so every
                        // flux term M_u*(scalar_u - env) evaluates to ~0 in rhsForcing.
                        // Setting s_u=0 with s≈0.8 gave flux=M_u*(0-0.8)→13M °C blow-up.
                        m.s_u.x[i][j][k]   = m.s.x[i][j][k];
                        m.q_v_u.x[i][j][k] = m.c.x[i][j][k];
                        m.q_c_u.x[i][j][k] = 0.0;
                        m.v_u.x[i][j][k]   = m.v.x[i][j][k];
                        m.w_u.x[i][j][k]   = m.w.x[i][j][k];
                    }

                    double inv_a_u = 1.0 / (a_u * m.u_0);
                    m.u_u.x[i][j][k] = m.u.x[i][j][k] + m.M_u.x[i][j][k] * inv_a_u / safe_r_humid(r_h_i);

                    if(m.q_v_u.x[i][j][k] <= 0.0) m.q_v_u.x[i][j][k] = 0.0;
                    if(m.q_c_u.x[i][j][k] <= 0.0) m.q_c_u.x[i][j][k] = 0.0;
                    if(m.s_u.x[i][j][k] <= 0.0)   m.s_u.x[i][j][k]   = 0.0;

                    if(is_land(m.h, i, j, k) || t_u <= t_00){
                        m.E_u.x[i][j][k]   = 0.0; m.D_u.x[i][j][k] = 0.0; m.M_u.x[i][j][k] = 0.0;
                        m.q_v_u.x[i][j][k] = 0.0; m.q_c_u.x[i][j][k] = 0.0;
                        m.u_u.x[i][j][k]   = 0.0; m.v_u.x[i][j][k] = 0.0; m.w_u.x[i][j][k] = 0.0;
                        m.s_u.x[i][j][k]   = 1.0;
                        m.g_p.x[i][j][k]   = 0.0; m.c_u.x[i][j][k] = 0.0; m.e_l.x[i][j][k] = 0.0;
                    }
                } // end i
            }
        }
    }
/*
*
*/
// ==================== DOWNDRAFT MAIN RECURRENCE ====================
    void downdraftRecurrence() {
        using namespace AtomMoistConvection;

        #pragma omp parallel for collapse(2) schedule(static)
        for(int j = 0; j < m.jm; j++){
            for(int k = 0; k < m.km; k++){

                int i_lfs_loc = i_LFS_local[j][k];

                for(int i = i_lfs_loc-1; i >= 0; i--){

                    double t_u      = m.t.x[i][j][k] * m.t_0;
                    double step_ip1 = step[i+1];
                    double r_h_ip1  = m.r_humid.x[i+1][j][k];
                    double r_h_i    = m.r_humid.x[i][j][k];

                    double d_Md      = m.E_d.x[i][j][k] - m.D_d.x[i][j][k];
                    m.M_d.x[i][j][k] = m.M_d.x[i+1][j][k] - d_Md * step_ip1;
                    if(t_u <= t_00_minus2) m.M_d.x[i][j][k] = 0.0;

                    double L_latent = (t_u >= m.t_0) ? m.lv : m.ls;
                    double M_d_ip1  = m.M_d.x[i+1][j][k];

                    double dummy_q_v_d = M_d_ip1 * m.q_v_d.x[i+1][j][k]
                        - step_ip1 * (m.E_d.x[i+1][j][k] * m.c.x[i+1][j][k]
                        - m.D_d.x[i+1][j][k] * m.q_v_d.x[i+1][j][k]
                        + r_h_ip1 * m.e_d.x[i+1][j][k]);

                    double dummy_v_d = M_d_ip1 * m.v_d.x[i+1][j][k]
                        - step_ip1 * (m.E_d.x[i+1][j][k] * m.v.x[i+1][j][k]
                        - m.D_d.x[i+1][j][k] * m.v_d.x[i+1][j][k]);

                    double dummy_w_d = M_d_ip1 * m.w_d.x[i+1][j][k]
                        - step_ip1 * (m.E_d.x[i+1][j][k] * m.w.x[i+1][j][k]
                        - m.D_d.x[i+1][j][k] * m.w_d.x[i+1][j][k]);

                    double dummy_s_d = (M_d_ip1 * m.s_d.x[i+1][j][k]
                        - step_ip1 * (m.E_d.x[i+1][j][k] * m.s.x[i+1][j][k]
                        - m.D_d.x[i+1][j][k] * m.s_d.x[i+1][j][k]
                        - L_latent * r_h_ip1 * m.e_d.x[i+1][j][k])) / m.s_0;

                    double M_d_i = m.M_d.x[i][j][k];
                    if(fabs(M_d_i) > coeff_recurr){
                        double inv_M_d = 1.0 / M_d_i;
                        m.q_v_d.x[i][j][k] = dummy_q_v_d * inv_M_d;
                        m.v_d.x[i][j][k]   = dummy_v_d   * inv_M_d;
                        m.w_d.x[i][j][k]   = dummy_w_d   * inv_M_d;
                        m.s_d.x[i][j][k]   = dummy_s_d   * inv_M_d;
                    } else {
                        // Mirrors the updraft else-branch above. When |M_d| ≤ coeff_recurr
                        // the downdraft is negligible, but skipping the assignment leaves
                        // stale q_v_d / s_d from the previous MoistConvection call. The
                        // next iteration then feeds that stale value through the D_d*q_v_d
                        // term in dummy_q_v_d and through the M_d*(q_v_d-c) flux in
                        // rhsForcing — amplifying it by the M_d ratio across marginal cells
                        // (observed blow-up at 1°N 119°E, 0–38 m once spin-up was off).
                        m.q_v_d.x[i][j][k] = m.c.x[i][j][k];
                        m.v_d.x[i][j][k]   = m.v.x[i][j][k];
                        m.w_d.x[i][j][k]   = m.w.x[i][j][k];
                        m.s_d.x[i][j][k]   = m.s.x[i][j][k];
                    }

                    double inv_a_d = 1.0 / (a_d * m.u_0);
                    m.u_d.x[i][j][k] = m.u.x[i][j][k] + m.M_d.x[i][j][k] * inv_a_d / safe_r_humid(r_h_i);

                    if(m.q_v_d.x[i][j][k] <= 0.0) m.q_v_d.x[i][j][k] = 0.0;
                    if(m.s_d.x[i][j][k] <= 0.0)   m.s_d.x[i][j][k]   = 0.0;

                    if(is_land(m.h, i, j, k) || t_u <= t_00){
                        m.E_d.x[i][j][k]   = 0.0; m.D_d.x[i][j][k] = 0.0; m.M_d.x[i][j][k] = 0.0;
                        m.q_v_d.x[i][j][k] = 0.0;
                        m.u_d.x[i][j][k]   = 0.0; m.v_d.x[i][j][k] = 0.0; m.w_d.x[i][j][k] = 0.0;
                        m.s_d.x[i][j][k]   = 1.0; m.e_d.x[i][j][k] = 0.0; m.e_p.x[i][j][k] = 0.0;
                    }

                    // Precipitation from convection (skip i=0: flux arrives unchanged)
                    if(i > 0)
                        m.P_conv.x[i][j][k] = std::max(0.0, m.P_conv.x[i+1][j][k]
                            + step[i] * m.r_humid.x[i][j][k]
                            * (m.g_p.x[i][j][k] - m.e_d.x[i][j][k] - m.e_p.x[i][j][k]));

                } // end i

                m.P_conv.x[0][j][k] = m.P_conv.x[1][j][k];
            }
        }
    }
/*
*
*/
// ==================== CAPE ====================
    void computeCAPE() {
        using namespace AtomMoistConvection;

        m.CAPE.assign(m.im, 0.0);

        for(int k = 0; k < m.km; k++){
            for(int j = 0; j < m.jm; j++){

                int i_base = i_Base_local[j][k];
                int i_lfs  = i_LFS_local[j][k];

                // cloud base contribution (mirrors findCloudBase() in MoistConvShall)
                {
                    double t_u       = m.t.x[i_base][j][k] * m.t_0;
                    double t_u_add   = t_u + t_add_u;
                    double t_vir     = t_u_add * (1.0 + alf * m.q_v_u.x[i_base][j][k]
                                       - m.q_c_u.x[i_base][j][k]);
                    double t_vir_env = std::max(t_u * (1.0 + alf * m.c.x[i_base][j][k]
                                       - cloud.x[i_base][j][k]), 1.0);
                    double rm        = m.rad.z[i_base];
                    double exp_rm    = 1.0 / (rm + 1.0);
                    m.CAPE[i_base]   = m.g * step[i_base]
                                       / exp_rm * (t_vir - t_vir_env) / t_vir_env;
                }

                // accumulate CAPE upward through the cloud layer (mirrors convectiveUpdraft() in MoistConvShall)
                for(int i = i_base; i < i_lfs && i < m.im-1; i++){
                    double t_u       = m.t.x[i][j][k] * m.t_0;
                    double t_u_add   = t_u + t_add_u;
                    double t_vir     = t_u_add * (1.0 + alf * m.q_v_u.x[i][j][k]
                                       - m.q_c_u.x[i][j][k]);
                    double t_vir_env = std::max(t_u * (1.0 + alf * m.c.x[i][j][k]
                                       - cloud.x[i][j][k]), 1.0);
                    double rm        = m.rad.z[i];
                    double exp_rm    = 1.0 / (rm + 1.0);
                    m.CAPE[i+1]      = m.CAPE[i] + m.g * step[i]
                                       / exp_rm * (t_vir - t_vir_env) / t_vir_env;
                }
            }
        }
    }
/*
*
*/
// ==================== RHS FORCING ====================
    void rhsForcing() {
        #pragma omp parallel for collapse(2)
        for(int j = 0; j < m.jm; j++){
            for(int k = 0; k < m.km; k++){
                for(int i = 0; i < m.im - 1; i++){

                    double L_latent = (m.t.x[i][j][k] * m.t_0 >= m.t_0) ? m.lv : m.ls;                  // J/kg
                    // Use the same physical floor as c_u/e_l (0.01 kg/m³ ≈ air at 50 hPa).
                    // safe_r_humid's 1e-6 floor allowed 1/r → 1e6 amplification of flux
                    // divergence at near-vacuum altitudes if M_u/M_d leaked upward.
                    double inv_step_rh = 1.0 / (step[i] * std::max(m.r_humid.x[i][j][k], 0.01));        // m²/kg

                    double flux_s_ip1 = m.M_u.x[i+1][j][k] * (m.s_u.x[i+1][j][k] - m.s.x[i+1][j][k])    // kg/(m²s)
                                    + m.M_d.x[i+1][j][k] * (m.s_d.x[i+1][j][k] - m.s.x[i+1][j][k]);
                    double flux_s_i   = m.M_u.x[i][j][k] * (m.s_u.x[i][j][k] - m.s.x[i][j][k])          // kg/(m²s)    s non-dimensional
                                    + m.M_d.x[i][j][k] * (m.s_d.x[i][j][k] - m.s.x[i][j][k]);

                    double conv_src = m.c_u.x[i][j][k] - m.e_d.x[i][j][k]                               // (kg/kg)/s
                                  - m.e_l.x[i][j][k] - m.e_p.x[i][j][k];

                    m.MC_t.x[i][j][k] = -(flux_s_ip1 - flux_s_i) * inv_step_rh * m.t_0                  // K/s
                        + (L_latent / m.cp_l) * conv_src* m.t_0;                                        // K/s



                    double flux_q_ip1 = m.M_u.x[i+1][j][k] * (m.q_v_u.x[i+1][j][k] - m.c.x[i+1][j][k])  // kg/(m²s)(kg/kg)
                                    + m.M_d.x[i+1][j][k] * (m.q_v_d.x[i+1][j][k] - m.c.x[i+1][j][k]);
                    double flux_q_i   = m.M_u.x[i][j][k] * (m.q_v_u.x[i][j][k] - m.c.x[i][j][k])
                                    + m.M_d.x[i][j][k] * (m.q_v_d.x[i][j][k] - m.c.x[i][j][k]);

                    m.MC_q.x[i][j][k] = - (flux_q_ip1 - flux_q_i) * inv_step_rh                         // (kg/kg)/s
                                        - conv_src;                                                     // (kg/kg)/s




                    double flux_v_ip1 = m.M_u.x[i+1][j][k] * (m.v_u.x[i+1][j][k] - m.v.x[i+1][j][k])    // kg/(m²s)
                                    + m.M_d.x[i+1][j][k] * (m.v_d.x[i+1][j][k] - m.v.x[i+1][j][k]);
                    double flux_v_i   = m.M_u.x[i][j][k] * (m.v_u.x[i][j][k] - m.v.x[i][j][k])
                                    + m.M_d.x[i][j][k] * (m.v_d.x[i][j][k] - m.v.x[i][j][k]);

                    m.MC_v.x[i][j][k] = -(flux_v_ip1 - flux_v_i) * inv_step_rh * m.u_0;                 // (m/s)(1/s)



                    double flux_w_ip1 = m.M_u.x[i+1][j][k] * (m.w_u.x[i+1][j][k] - m.w.x[i+1][j][k])
                                    + m.M_d.x[i+1][j][k] * (m.w_d.x[i+1][j][k] - m.w.x[i+1][j][k]);
                    double flux_w_i   = m.M_u.x[i][j][k] * (m.w_u.x[i][j][k] - m.w.x[i][j][k])
                                    + m.M_d.x[i][j][k] * (m.w_d.x[i][j][k] - m.w.x[i][j][k]);

                    m.MC_w.x[i][j][k] = -(flux_w_ip1 - flux_w_i) * inv_step_rh * m.u_0;                 // (m/s)(1/s)
                }
            }
        }
    }
/*
*
*/
// ==================== SIMPSON INTEGRATION OF -c * E_u OVER i ====================
    // Accumulates ∫ -c(i,j,k)·E_u(i,j,k) dz from i_beg upward in a single (j,k) column.
    // Uses composite non-uniform Simpson's 1/3 rule on the variable vertical steps step[i].
    // Falls back to the trapezoidal rule for one remaining interval when the count is odd.
    // Returns the first i-index where the cumulative sum becomes positive, or -1 if never.
    int simpsonCxEu(int j, int k, int i_beg, int i_end) const {
        if(i_end <= i_beg) return -1;
        double cumsum = 0.0;
        int i = i_beg;
        // --- composite non-uniform Simpson pairs ---
        for(; i + 2 <= i_end; i += 2){
            double h1 = step[i];
            double h2 = step[i+1];
            double f0 = -(m.c.x[i  ][j][k] * m.E_u.x[i  ][j][k]);
            double f1 = -(m.c.x[i+1][j][k] * m.E_u.x[i+1][j][k]);
            double f2 = -(m.c.x[i+2][j][k] * m.E_u.x[i+2][j][k]);
            // non-uniform 3-point quadratic rule (reduces to h/3*(f0+4f1+f2) when h1==h2)
            double pair_contrib = (h1 + h2) / 6.0
                * ((2.0 - h2 / h1) * f0
                +  (h1 + h2) * (h1 + h2) / (h1 * h2) * f1
                +  (2.0 - h1 / h2) * f2);
            if(cumsum + pair_contrib > 0.0){
                // Refine: check trapezoidal partial [i, i+1] to find which sub-interval crosses
                return (cumsum + 0.5 * h1 * (f0 + f1) > 0.0) ? i + 1 : i + 2;
            }
            cumsum += pair_contrib;
        }
        // --- one remaining interval: trapezoidal ---
        if(i + 1 <= i_end){
            double f0 = -(m.c.x[i  ][j][k] * m.E_u.x[i  ][j][k]);
            double f1 = -(m.c.x[i+1][j][k] * m.E_u.x[i+1][j][k]);
            if(cumsum + 0.5 * step[i] * (f0 + f1) > 0.0) return i + 1;
        }
        return -1;
    }
/*
*
*/
    // Compute i_Base_local[j][k] = first i where ∫ -c·E_u dz becomes positive, per column.
    // Integrates from i_deep_beg_local[j][k] to i_LFS_local[j][k].
    // Rejects high-altitude results: cloud base must be below p_stat_base (800 hPa).
    void CloudBaseUpdate() {
        using namespace AtomMoistConvection;
        #pragma omp parallel for collapse(2) schedule(static)
        for(int j = 1; j < m.jm-1; j++){
            for(int k = 1; k < m.km-1; k++){
                const int i_beg = i_deep_beg_local[j][k];
                const int i_end = i_LFS_local[j][k];
                if(i_beg < 0 || i_end <= i_beg) continue;
                const int i_simpson = simpsonCxEu(j, k, i_beg, i_end);
                if(i_simpson >= 0 && i_simpson < m.im
                        && m.p_stat.x[i_simpson][j][k] <= p_stat_Cloud_Base   // reject high-altitude outliers
                        // Cause #8: old condition (i_simpson > i_Base_local) only raised
                        // i_base, widening the stale-initUpdraft region below cloud base
                        // on every call.  Take the lowest (surface-nearest) valid base instead.
                        && (i_Base_local[j][k] <= 0 || i_simpson <= i_Base_local[j][k])){
                    i_Base_local[j][k] = i_simpson;
                    m.CloudBase.x[i_simpson][j][k] = height_table[i_simpson];
                }
            }
        }
    }
/*
*
*/
// ==================== SUB-TERRAIN FILL ====================
    void subTerrainFill() {
        #pragma omp parallel for collapse(2)
        for(int j = 0; j < m.jm; j++){
            for(int k = 0; k < m.km; k++){
              if(!land_surf[j * m.km + k]) continue;

                int i_mount = m.i_topography[j][k];

                // Cache mount-level values once
                double pc_m = m.P_conv.x[i_mount][j][k];
                double eu_m = m.E_u.x[i_mount][j][k], ed_m = m.E_d.x[i_mount][j][k];
                double du_m = m.D_u.x[i_mount][j][k], dd_m = m.D_d.x[i_mount][j][k];
                double mu_m = m.M_u.x[i_mount][j][k], md_m = m.M_d.x[i_mount][j][k];
                double mct_m = m.MC_t.x[i_mount][j][k], mcq_m = m.MC_q.x[i_mount][j][k];
                double mcv_m = m.MC_v.x[i_mount][j][k], mcw_m = m.MC_w.x[i_mount][j][k];
                double qvu_m = m.q_v_u.x[i_mount][j][k], qvd_m = m.q_v_d.x[i_mount][j][k];
                double qcu_m = m.q_c_u.x[i_mount][j][k];
                double cu_m = m.c_u.x[i_mount][j][k], gp_m = m.g_p.x[i_mount][j][k];
                double el_m = m.e_l.x[i_mount][j][k], ep_m = m.e_p.x[i_mount][j][k];
                double edm = m.e_d.x[i_mount][j][k];
                double uu_m = m.u_u.x[i_mount][j][k], ud_m = m.u_d.x[i_mount][j][k];
                double vu_m = m.v_u.x[i_mount][j][k], vd_m = m.v_d.x[i_mount][j][k];
                double wu_m = m.w_u.x[i_mount][j][k], wd_m = m.w_d.x[i_mount][j][k];
                double s_m = m.s.x[i_mount][j][k], su_m = m.s_u.x[i_mount][j][k];
                double sd_m = m.s_d.x[i_mount][j][k];

                for(int i = i_mount - 1; i >= 0; i--){
                    m.P_conv.x[i][j][k] = pc_m;
                    m.E_u.x[i][j][k] = eu_m; m.E_d.x[i][j][k] = ed_m;
                    m.D_u.x[i][j][k] = du_m; m.D_d.x[i][j][k] = dd_m;
                    m.M_u.x[i][j][k] = mu_m; m.M_d.x[i][j][k] = md_m;
                    m.MC_t.x[i][j][k] = mct_m; m.MC_q.x[i][j][k] = mcq_m;
                    m.MC_v.x[i][j][k] = mcv_m; m.MC_w.x[i][j][k] = mcw_m;
                    m.q_v_u.x[i][j][k] = qvu_m; m.q_v_d.x[i][j][k] = qvd_m;
                    m.q_c_u.x[i][j][k] = qcu_m;
                    m.c_u.x[i][j][k] = cu_m; m.g_p.x[i][j][k] = gp_m;
                    m.e_l.x[i][j][k] = el_m; m.e_p.x[i][j][k] = ep_m;
                    m.e_d.x[i][j][k] = edm;
                    m.u_u.x[i][j][k] = uu_m; m.u_d.x[i][j][k] = ud_m;
                    m.v_u.x[i][j][k] = vu_m; m.v_d.x[i][j][k] = vd_m;
                    m.w_u.x[i][j][k] = wu_m; m.w_d.x[i][j][k] = wd_m;
                    m.s.x[i][j][k] = s_m; m.s_u.x[i][j][k] = su_m;
                    m.s_d.x[i][j][k] = sd_m;
                }
            }
        }
    }
};
/*
*
*/
