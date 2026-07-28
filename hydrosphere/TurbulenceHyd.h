/*
 * Ocean General Circulation Model (OGCM) applied to turbulent flow
 * Standalone turbulence class: k-epsilon (Chien 1982),
 *                              k-omega   (Wilcox 1988),
 *                              k-omega SST (Menter 1994)
 * Turbulence arrays (tke, dis, nue, prod, tke_source, dis_source) are owned
 * by cHydrosphereModel and accessed here through the model reference m.
 *
 * Coordinate convention:
 *   i=0 is the deepest level (ocean floor / land cells).
 *   i=im-1 is the ocean surface.
 *   i_bathymetry[j][k] is the topmost land cell; i_bath+1 is the first water cell.
 *   "wall" = seafloor at i_bathymetry[j][k]; BC is asserted at i_bath+1.
 *
 * Non-dimensionalisation:
 *   k*  = k_phys / u_0²
 *   ω*  = ω_phys · L_hyd / u_0        (k-ω, k-ω SST)
 *   ε*  = ε_phys · L_hyd / u_0³       (k-ε)
 *   ν*  = ν_T / (u_0 · L_hyd)
*/

#pragma once

#include "cHydrosphereModel.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>

using namespace AtomUtils;

class TurbulenceHyd {
public:
    enum Model { k_epsilon, k_omega, k_omega_SST };

    explicit TurbulenceHyd(cHydrosphereModel& model,
                           double             vel_star_ref = 0.01,      // reference friction velocity [m/s]
                           double             z_0          = 1.0e-3,    // roughness length [m]
                           double             nue_water    = 1.0e-6)    // kin. viscosity seawater [m²/s]
        : m(model),
          turb_model(parse_model(model.turb_model)),
          vel_star_ref(vel_star_ref),
          z_0(z_0),
          nue_water(nue_water),
          re_turb(vel_star_ref * z_0 / nue_water)
    {}

    // -----------------------------------------------------------------------
    void init() {
        using namespace std;
        cout << endl << "      OGCM: TurbulenceHyd::init" << endl;
        print_model_name();

        auto begin = chrono::high_resolution_clock::now();

        m.re_turb = re_turb;
        compute_vel_star();
        init_fields();
        apply_wall_bc();
        compute_sources();
        zero_land_cells();
        clamp_nue();

        auto end     = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for TurbulenceHyd::init\n",
               elapsed.count() * 1e-9);
        cout << "      OGCM: TurbulenceHyd::init ended" << endl;
    }

    // -----------------------------------------------------------------------
    void run() {
        using namespace std;
        cout << endl << "      OGCM: TurbulenceHyd::run" << endl;
        print_model_name();

        auto begin = chrono::high_resolution_clock::now();

        compute_vel_star();
        compute_sources();
        zero_land_cells();
        clamp_nue();

        auto end     = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for TurbulenceHyd::run\n",
               elapsed.count() * 1e-9);
        cout << "      OGCM: TurbulenceHyd::run ended" << endl;
    }

    // -----------------------------------------------------------------------
    // Compute the friction velocity u_τ per (j,k) column from the horizontal
    // current speed at the first water cell above the seafloor.
    //   U_horiz = u_0 · sqrt(v² + w²)  at i = i_bath + 1
    //   u_τ     = max(vel_star_min, sqrt(C_D) · U_horiz)
    // C_D = 0.002 (neutral drag coefficient); vel_star_min = 0.001 m/s.
    void compute_vel_star() {
        constexpr double C_D          = 0.002;
        constexpr double vel_star_min = 0.001;   // [m/s]

        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                const int i_bath  = m.i_bathymetry[j][k];
                const int i_surf  = std::min(i_bath + 1, m.im - 1);
                const double v_s  = m.v.x[i_surf][j][k];
                const double w_s  = m.w.x[i_surf][j][k];
                const double U_horiz = m.u_0 * std::sqrt(v_s * v_s + w_s * w_s);
                m.vel_star.y[j][k] = std::max(vel_star_min, std::sqrt(C_D) * U_horiz);
            }
        }
    }

    // -----------------------------------------------------------------------
    // Wall boundary condition for the turbulence scalars at the seafloor.
    // "Wall" = first water cell above the floor at i = i_bathymetry[j][k] + 1.
    // Public so run_3D_loop can reassert the wall values after bcRadius
    // overwrites tke[i_w1] / dis[i_w1] with cubic Neumann extrapolation.
    void apply_wall_bc() {
        const double nd_omega = m.L_hyd / m.u_0;   // ω_phys [s⁻¹] → ω*

        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                const int i_bath = m.i_bathymetry[j][k];
                // Skip columns that have no water cells above the seafloor.
                // When i_bath >= im-1 the floor is at or above the surface, so
                // i_w1 = min(i_bath+1, im-1) = im-1 = i_bath → y1 = 0 → floor
                // at 1e-6 m → enormous ω_wall written into a land cell.
                if (i_bath >= m.im - 1) continue;
                const int    i_w1 = i_bath + 1;   // always < im-1 after the guard above
                const double vs   = m.vel_star.y[j][k];

                // Height of first water cell above seafloor [m].
                // rad.z[i] is the dimensionless stretched coordinate; multiply by
                // L_hyd to get a physical length scale.  Floor at 1e-6 m to prevent
                // division by zero in the wall-omega formula.
                const double y1 = std::max(
                    (m.rad.z[i_w1] - m.rad.z[i_bath]) * m.L_hyd, 1.0e-6);

                // tke wall BC — log-layer equilibrium k = u_τ²/√Cμ
                m.tke.x[i_w1][j][k] = vs * vs / (std::sqrt(C_nue) * m.u_0 * m.u_0);

                // nue wall BC — turbulent viscosity vanishes at the no-slip seafloor
                m.nue.x[i_w1][j][k] = 0.0;
                // Also zero the land ghost cell itself
                if (i_bath >= 0 && i_bath < m.im)
                    m.nue.x[i_bath][j][k] = 0.0;

                if (turb_model == k_epsilon) {
                    // Chien (1982): ε_wall = 0 at the no-slip wall
                    m.dis.x[i_w1][j][k] = 0.0;
                } else if (turb_model == k_omega) {
                    // Wilcox 1988: ω_wall = 6ν / (β₁ y₁²),  β₁ = 0.075
                    const double om_phys = 6.0 * nue_water / (0.075 * y1 * y1);
                    m.dis.x[i_w1][j][k] = om_phys * nd_omega;
                } else if (turb_model == k_omega_SST) {
                    // Menter 1994: ω_wall = 60ν / (β₁ y₁²),  β₁ = 0.0333
                    const double om_phys = 60.0 * nue_water / (0.0333 * y1 * y1);
                    m.dis.x[i_w1][j][k] = om_phys * nd_omega;
                }
            }
        }
    }

private:
    cHydrosphereModel& m;
    Model  turb_model;
    double vel_star_ref;
    double z_0;
    double nue_water;
    double re_turb;

    // ABL (benthic boundary layer) model constants
    static constexpr double C_nue   = 0.028;
    static constexpr double Karman  = 0.42;
    static constexpr double dis_min = 1.0e-10;

    // Map string parameter to enum
    static Model parse_model(const std::string& s) {
        if (s == "k_epsilon")   return k_epsilon;
        if (s == "k_omega")     return k_omega;
        if (s == "k_omega_SST") return k_omega_SST;
        return k_omega_SST;
    }

    static double blend(double inner, double outer, double F1) {
        return F1 * inner + (1.0 - F1) * outer;
    }

    void print_model_name() const {
        using namespace std;
        switch (turb_model) {
            case k_epsilon:   cout << "      k-epsilon turbulence model (Chien 1982)" << endl; break;
            case k_omega:     cout << "      k-omega turbulence model (Wilcox 1988)" << endl;  break;
            case k_omega_SST: cout << "      k-omega SST turbulence model (Menter 1994)" << endl; break;
        }
    }

    // -----------------------------------------------------------------------
    // Initialise tke and dis throughout the ocean domain using a log-layer
    // profile from the seafloor upward.
    // Physical height above seafloor: z_above = (rad.z[i] - rad.z[i_bath]) * L_hyd [m].
    void init_fields() {
        const double inv_u02  = 1.0 / (m.u_0 * m.u_0);
        const double nd_omega = m.L_hyd / m.u_0;
        const double nd_eps   = m.L_hyd / (m.u_0 * m.u_0 * m.u_0);

        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 1; j < m.jm - 1; j++) {
            for (int k = 1; k < m.km - 1; k++) {
                const int    i_bath = m.i_bathymetry[j][k];
                const double vs     = m.vel_star.y[j][k];

                // Background floor: 1% of peak near-floor TKE
                const double k_bg     = 0.01 * vs * vs / std::sqrt(C_nue);  // [m²/s²]
                const double om_bg    = std::pow(C_nue, -0.25) * std::sqrt(k_bg) / (Karman * m.L_hyd);
                const double eps_bg   = std::pow(C_nue, 0.75) * std::pow(k_bg, 1.5) / (Karman * m.L_hyd);
                const double z_floor  = m.rad.z[std::max(i_bath, 0)] * m.L_hyd;

                for (int i = i_bath + 1; i < m.im; i++) {
                    const double z_i     = m.rad.z[i] * m.L_hyd;
                    const double z_above = std::max(z_i - z_floor, 1.0e-6);   // height above floor [m]

                    // Parabolic profile (height above floor / total water column)
                    const double z_col   = std::max(m.rad.z[m.im - 1] * m.L_hyd - z_floor, 1.0e-6);
                    const double z_frac  = std::min(z_above, z_col) / z_col;

                    const double k_bbl   = vs * vs / std::sqrt(C_nue) * std::pow(1.0 - z_frac, 2);
                    const double k_phys  = std::max(k_bbl, k_bg);

                    double dis_nd = 0.0;
                    if (turb_model == k_epsilon) {
                        const double eps_phys = std::pow(C_nue, 0.75)
                            * std::pow(k_phys, 1.5) / (Karman * z_above);
                        dis_nd = std::max(eps_phys, eps_bg) * nd_eps;
                    } else {
                        const double om_phys = std::pow(C_nue, -0.25)
                            * std::sqrt(k_phys) / (Karman * z_above);
                        dis_nd = std::max(om_phys, om_bg) * nd_omega;
                    }

                    m.tke.x[i][j][k]  = std::max(0.0, k_phys * inv_u02);
                    m.dis.x[i][j][k]  = std::max(dis_min, dis_nd);
                    m.tken.x[i][j][k] = m.tke.x[i][j][k];
                    m.disn.x[i][j][k] = m.dis.x[i][j][k];
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    void clamp_nue() {
        const double nue_max = 1000.0 / (m.u_0 * m.L_hyd);
        #pragma omp parallel for collapse(3) schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                for (int k = 0; k < m.km; k++) {
                    m.nue.x[i][j][k] = std::clamp(m.nue.x[i][j][k], 0.0, nue_max);
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    void zero_land_cells() {
        #pragma omp parallel for collapse(3) schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                for (int k = 0; k < m.km; k++) {
                    if (is_land(m.h, i, j, k)) {
                        m.tke.x[i][j][k]        = 0.0;
                        m.tken.x[i][j][k]        = 0.0;
                        m.dis.x[i][j][k]         = 0.0;
                        m.disn.x[i][j][k]        = 0.0;
                        m.nue.x[i][j][k]         = 0.0;
                        m.prod.x[i][j][k]        = 0.0;
                        m.tke_source.x[i][j][k]  = 0.0;
                        m.dis_source.x[i][j][k]  = 0.0;
                    }
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    void compute_sources() {
        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 1; j < m.jm - 1; j++) {
            for (int k = 1; k < m.km - 1; k++) {
                const int    i_bath    = m.i_bathymetry[j][k];
                const double z_floor   = m.rad.z[std::max(i_bath, 0)] * m.L_hyd;

                const double vs        = m.vel_star.y[j][k];
                // Background tke floor [dimensionless]
                const double k_bg_nd   = 0.01 * vs * vs / (std::sqrt(C_nue) * m.u_0 * m.u_0);

                for (int i = i_bath + 1; i < m.im - 1; i++) {
                    if (is_land(m.h, i, j, k)) continue;

                    // ---- geometry ----
                    const double rm       = m.rad.z[i];
                    const double exp_rm   = 1.0 / (rm + 1.0);
                    double sinthe         = sin(m.the.z[j]);
                    if (sinthe == 0.0) sinthe = 1.0e-5;
                    const double rmsinthe = rm * sinthe;
                    const double inv_2dthe    = 1.0 / (2.0 * m.dthe);
                    const double inv_2dphi    = 1.0 / (2.0 * m.dphi);
                    const double inv_rm2dthe  = inv_2dthe / rm;
                    const double inv_rmsinthe2dphi = inv_2dphi / rmsinthe;

                    m.tke.x[i][j][k] = std::max(0.0,     m.tke.x[i][j][k]);
                    m.dis.x[i][j][k] = std::max(dis_min, m.dis.x[i][j][k]);

                    // ---- velocity gradients ----
                    const double dudr   = (m.rc1m[i]*m.u.x[i-1][j][k] + m.rc10[i]*m.u.x[i][j][k] + m.rc1p[i]*m.u.x[i+1][j][k]) * exp_rm;
                    const double dvdr   = (m.rc1m[i]*m.v.x[i-1][j][k] + m.rc10[i]*m.v.x[i][j][k] + m.rc1p[i]*m.v.x[i+1][j][k]) * exp_rm;
                    const double dwdr   = (m.rc1m[i]*m.w.x[i-1][j][k] + m.rc10[i]*m.w.x[i][j][k] + m.rc1p[i]*m.w.x[i+1][j][k]) * exp_rm;
                    const double dudthe = (m.u.x[i][j+1][k] - m.u.x[i][j-1][k]) * inv_rm2dthe;
                    const double dvdthe = (m.v.x[i][j+1][k] - m.v.x[i][j-1][k]) * inv_rm2dthe;
                    const double dwdthe = (m.w.x[i][j+1][k] - m.w.x[i][j-1][k]) * inv_rm2dthe;
                    const double dudphi = (m.u.x[i][j][k+1] - m.u.x[i][j][k-1]) * inv_rmsinthe2dphi;
                    const double dvdphi = (m.v.x[i][j][k+1] - m.v.x[i][j][k-1]) * inv_rmsinthe2dphi;
                    const double dwdphi = (m.w.x[i][j][k+1] - m.w.x[i][j][k-1]) * inv_rmsinthe2dphi;

                    // ---- Neumann BC for tke/dis gradients at land faces ----
                    const bool neumann_bot = (i > 0      && is_land(m.h, i-1, j, k));
                    const bool neumann_jm1 = (j > 0      && is_land(m.h, i, j-1, k));
                    const bool neumann_jp1 = (j < m.jm-1 && is_land(m.h, i, j+1, k));
                    const bool neumann_km1 = (k > 0      && is_land(m.h, i, j, k-1));
                    const bool neumann_kp1 = (k < m.km-1 && is_land(m.h, i, j, k+1));

                    const double tke_im1 = neumann_bot ? m.tke.x[i][j][k] : m.tke.x[i-1][j][k];
                    const double dis_im1 = neumann_bot ? m.dis.x[i][j][k] : m.dis.x[i-1][j][k];
                    const double tke_jm1 = neumann_jm1 ? m.tke.x[i][j][k] : m.tke.x[i][j-1][k];
                    const double tke_jp1 = neumann_jp1 ? m.tke.x[i][j][k] : m.tke.x[i][j+1][k];
                    const double dis_jm1 = neumann_jm1 ? m.dis.x[i][j][k] : m.dis.x[i][j-1][k];
                    const double dis_jp1 = neumann_jp1 ? m.dis.x[i][j][k] : m.dis.x[i][j+1][k];
                    const double tke_km1 = neumann_km1 ? m.tke.x[i][j][k] : m.tke.x[i][j][k-1];
                    const double tke_kp1 = neumann_kp1 ? m.tke.x[i][j][k] : m.tke.x[i][j][k+1];
                    const double dis_km1 = neumann_km1 ? m.dis.x[i][j][k] : m.dis.x[i][j][k-1];
                    const double dis_kp1 = neumann_kp1 ? m.dis.x[i][j][k] : m.dis.x[i][j][k+1];

                    const double dtkedr   = (m.rc1m[i]*tke_im1 + m.rc10[i]*m.tke.x[i][j][k] + m.rc1p[i]*m.tke.x[i+1][j][k]) * exp_rm;
                    const double ddisdr   = (m.rc1m[i]*dis_im1 + m.rc10[i]*m.dis.x[i][j][k] + m.rc1p[i]*m.dis.x[i+1][j][k]) * exp_rm;
                    const double dtkedthe = (tke_jp1 - tke_jm1) * inv_rm2dthe;
                    const double ddisdthe = (dis_jp1 - dis_jm1) * inv_rm2dthe;
                    const double dtkedphi = (tke_kp1 - tke_km1) * inv_rmsinthe2dphi;
                    const double ddisdphi = (dis_kp1 - dis_km1) * inv_rmsinthe2dphi;

                    // ---- production P_k ----
                    double cnue;
                    {
                        const double dis_here = std::max(m.dis.x[i][j][k], dis_min);
                        const double tke_here = std::max(m.tke.x[i][j][k], 0.0);
                        const double nue_max  = 1000.0 / (m.u_0 * m.L_hyd);
                        cnue = (turb_model == k_epsilon)
                             ? C_nue * tke_here * tke_here / dis_here
                             : tke_here / dis_here;
                        cnue = std::min(cnue, nue_max);
                    }
                    const double der = 0.66667 * (dudr + dvdthe + dwdphi);
                    m.prod.x[i][j][k] = std::max(0.0,
                          (cnue * (2.0 * dudr   - der) - 0.66667 * m.tke.x[i][j][k]) * dudr
                        + (cnue * (dudthe + dvdr))                                   * dudthe
                        + (cnue * (dudphi + dwdr))                                   * dudphi
                        + (cnue * (2.0 * dvdthe - der) - 0.66667 * m.tke.x[i][j][k]) * dvdthe
                        + (cnue * (dvdr + dudthe))                                   * dvdr
                        + (cnue * (dvdphi + dwdthe))                                 * dvdphi
                        + (cnue * (2.0 * dwdphi - der) - 0.66667 * m.tke.x[i][j][k]) * dwdphi
                        + (cnue * (dwdr + dudphi))                                   * dwdr
                        + (cnue * (dwdthe + dvdphi))                                 * dwdthe);

                    // ---- vorticity magnitude ----
                    const double W12 = dudthe - dvdr;
                    const double W13 = dudphi - dwdr;
                    const double W23 = dvdphi - dwdthe;
                    const double Omega = std::sqrt(W12*W12 + W13*W13 + W23*W23);

                    // ================================================================
                    const double z_above = std::max(m.rad.z[i] * m.L_hyd - z_floor, 1.0e-6);

                    if (turb_model == k_epsilon) {
                        compute_k_epsilon(i, j, k, z_floor, dtkedr, ddisdr,
                            dtkedthe, ddisdthe, dtkedphi, ddisdphi, Omega);
                    } else if (turb_model == k_omega) {
                        compute_k_omega(i, j, k,
                            dtkedr, ddisdr, dtkedthe, ddisdthe, dtkedphi, ddisdphi,
                            Omega, Omega);
                    } else if (turb_model == k_omega_SST) {
                        compute_k_omega_SST(i, j, k, z_floor, rm, sinthe,
                            dtkedr, ddisdr, dtkedthe, ddisdthe, dtkedphi, ddisdphi,
                            Omega, Omega);
                    }

                    const double nue_max    = 1000.0 / (m.u_0 * m.L_hyd);
                    m.tke.x[i][j][k] = std::max(k_bg_nd, m.tke.x[i][j][k]);
                    m.dis.x[i][j][k] = std::max(dis_min, m.dis.x[i][j][k]);
                    m.nue.x[i][j][k] = std::clamp(m.nue.x[i][j][k], 0.0, nue_max);
                    const double dis_src_max = 20.0 * m.dis.x[i][j][k] * m.dis.x[i][j][k];
                    m.dis_source.x[i][j][k] = std::clamp(m.dis_source.x[i][j][k],
                                                          -dis_src_max, dis_src_max);
                    // Clamp tke_source symmetrically — 20×tke×dis covers the full
                    // range of P_k-Y_k without saturating legitimate physics.
                    const double tke_src_max = 20.0 * m.tke.x[i][j][k] * m.dis.x[i][j][k];
                    m.tke_source.x[i][j][k] = std::clamp(m.tke_source.x[i][j][k],
                                                          -tke_src_max, tke_src_max);
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // z_floor: physical height of the seafloor [m] = rad.z[i_bath] * L_hyd
    void compute_k_epsilon(int i, int j, int k, double z_floor,
        double dtkedr, double ddisdr,
        double dtkedthe, double ddisdthe,
        double dtkedphi, double ddisdphi,
        double /*Omega*/)
    {
        const double C_eps_1 = 1.35;
        const double C_eps_2 = 1.80;

        // Physical wall distance [m]; floor at L_hyd so 1/y_star² stays bounded
        const double y_phys = std::max(m.rad.z[i] * m.L_hyd - z_floor, m.L_hyd);
        const double y_star = y_phys / m.L_hyd;   // dimensionless, always ≥ 1

        m.nue.x[i][j][k] = C_nue * m.tke.x[i][j][k] * m.tke.x[i][j][k]
                          / std::max(m.dis.x[i][j][k], dis_min);

        const double d_plus = y_phys * m.vel_star.y[j][k] / nue_water;
        const double f_nue  = 1.0 - std::exp(-0.0115 * d_plus);

        const double Re_T = m.tke.x[i][j][k] * m.tke.x[i][j][k]
                            * m.u_0 * m.L_hyd
                            / std::max(m.dis.x[i][j][k] * nue_water, dis_min * nue_water);
        const double f_2  = 1.0 - 0.4 / 1.8 * std::exp(-Re_T * Re_T / 36.0);

        m.nue.x[i][j][k] *= f_nue;

        const double y_star2 = y_star * y_star;
        const double L_k = -2.0 * m.tke.x[i][j][k] / y_star2;
        const double L_w = -2.0 * m.dis.x[i][j][k] / y_star2 * std::exp(-d_plus / 2.0);

        const double P_k      = m.prod.x[i][j][k];
        const double Y_k      = m.dis.x[i][j][k];
        const double tke_safe = std::max(m.tke.x[i][j][k], dis_min);
        const double P_w      = C_eps_1 * m.dis.x[i][j][k] / tke_safe * P_k;
        const double Y_w      = C_eps_2 * f_2 * m.dis.x[i][j][k] * m.dis.x[i][j][k]
                                / tke_safe;

        m.tke_source.x[i][j][k] = P_k - Y_k + L_k / re_turb;
        m.dis_source.x[i][j][k] = P_w - Y_w + L_w / re_turb;
    }

    // -----------------------------------------------------------------------
    void compute_k_omega(int i, int j, int k,
        double dtkedr, double ddisdr,
        double dtkedthe, double ddisdthe,
        double dtkedphi, double ddisdphi,
        double Omega_mag, double /*W_unused*/)
    {
        const double bet_star = 0.09;
        const double gam      = 0.52;
        const double C_lim    = 0.875;
        const double bet_0    = 0.0708;

        const double rm       = m.rad.z[i];
        double sinthe         = sin(m.the.z[j]);
        if (sinthe == 0.0) sinthe = 1.0e-5;
        const double rmsinthe = rm * sinthe;
        const double exp_rm   = 1.0 / (rm + 1.0);
        const double inv_rm   = 1.0 / rm;

        const bool neumann_bot = (i > 0      && is_land(m.h, i-1, j, k));
        const bool neumann_jm1 = (j > 0      && is_land(m.h, i, j-1, k));
        const bool neumann_jp1 = (j < m.jm-1 && is_land(m.h, i, j+1, k));
        const bool neumann_km1 = (k > 0      && is_land(m.h, i, j, k-1));
        const bool neumann_kp1 = (k < m.km-1 && is_land(m.h, i, j, k+1));

        const double tke_im1 = neumann_bot ? m.tke.x[i][j][k] : m.tke.x[i-1][j][k];
        const double dis_im1 = neumann_bot ? m.dis.x[i][j][k] : m.dis.x[i-1][j][k];
        const double tke_jm1 = neumann_jm1 ? m.tke.x[i][j][k] : m.tke.x[i][j-1][k];
        const double tke_jp1 = neumann_jp1 ? m.tke.x[i][j][k] : m.tke.x[i][j+1][k];
        const double dis_jm1 = neumann_jm1 ? m.dis.x[i][j][k] : m.dis.x[i][j-1][k];
        const double dis_jp1 = neumann_jp1 ? m.dis.x[i][j][k] : m.dis.x[i][j+1][k];
        const double tke_km1 = neumann_km1 ? m.tke.x[i][j][k] : m.tke.x[i][j][k-1];
        const double tke_kp1 = neumann_kp1 ? m.tke.x[i][j][k] : m.tke.x[i][j][k+1];
        const double dis_km1 = neumann_km1 ? m.dis.x[i][j][k] : m.dis.x[i][j][k-1];
        const double dis_kp1 = neumann_kp1 ? m.dis.x[i][j][k] : m.dis.x[i][j][k+1];

        const double dtkedr_neu  = neumann_bot
            ? (m.tke.x[i+1][j][k] - m.tke.x[i][j][k]) / (m.rad.z[i+1] - m.rad.z[i]) * exp_rm
            : dtkedr;
        const double ddisdr_neu  = neumann_bot
            ? (m.dis.x[i+1][j][k] - m.dis.x[i][j][k]) / (m.rad.z[i+1] - m.rad.z[i]) * exp_rm
            : ddisdr;
        const double dtkedthe_neu = (tke_jp1 - tke_jm1) / (rm  * 2.0 * m.dthe);
        const double ddisdthe_neu = (dis_jp1 - dis_jm1) / (rm  * 2.0 * m.dthe);
        const double dtkedphi_neu = (tke_kp1 - tke_km1) / (rmsinthe * 2.0 * m.dphi);
        const double ddisdphi_neu = (dis_kp1 - dis_km1) / (rmsinthe * 2.0 * m.dphi);

        const double dudr   = (m.rc1m[i]*m.u.x[i-1][j][k] + m.rc10[i]*m.u.x[i][j][k] + m.rc1p[i]*m.u.x[i+1][j][k]) * exp_rm;
        const double dvdthe = (m.v.x[i][j+1][k] - m.v.x[i][j-1][k]) / (rm * 2.0 * m.dthe);
        const double dwdphi = (m.w.x[i][j][k+1] - m.w.x[i][j][k-1]) / (rmsinthe * 2.0 * m.dphi);

        const double dudthe_u = (m.u.x[i][j+1][k] - m.u.x[i][j-1][k]) / (rm       * 2.0 * m.dthe);
        const double dudphi_u = (m.u.x[i][j][k+1] - m.u.x[i][j][k-1]) / (rmsinthe * 2.0 * m.dphi);
        const double dvdr_u   = (m.rc1m[i]*m.v.x[i-1][j][k] + m.rc10[i]*m.v.x[i][j][k] + m.rc1p[i]*m.v.x[i+1][j][k]) * exp_rm;
        const double dwdr_u   = (m.rc1m[i]*m.w.x[i-1][j][k] + m.rc10[i]*m.w.x[i][j][k] + m.rc1p[i]*m.w.x[i+1][j][k]) * exp_rm;
        const double dvdphi_u = (m.v.x[i][j][k+1] - m.v.x[i][j][k-1]) / (rmsinthe * 2.0 * m.dphi);
        const double dwdthe_u = (m.w.x[i][j+1][k] - m.w.x[i][j-1][k]) / (rm       * 2.0 * m.dthe);

        const double S11 = dudr, S22 = dvdthe, S33 = dwdphi;
        const double S12 = 0.5 * (dudthe_u + dvdr_u);
        const double S13 = 0.5 * (dudphi_u + dwdr_u);
        const double S23 = 0.5 * (dvdphi_u + dwdthe_u);
        const double S_mag = std::sqrt(2.0 * (S11*S11 + S22*S22 + S33*S33
                                            + 2.0*(S12*S12 + S13*S13 + S23*S23)));

        const double w_hat = std::max(m.dis.x[i][j][k], C_lim * S_mag / std::sqrt(bet_star));

        double sig_d = dtkedr_neu * ddisdr_neu + dtkedthe_neu * ddisdthe_neu + dtkedphi_neu * ddisdphi_neu;
        sig_d = (sig_d <= 0.0) ? 0.0 : 0.125;

        const double chi_w = std::fabs(Omega_mag * Omega_mag * S_mag
            / std::pow(bet_star * std::max(m.dis.x[i][j][k], 1.0e-20), 3));
        const double f_bet  = (1.0 + 85.0 * chi_w) / (1.0 + 100.0 * chi_w);
        const double bet_wc = bet_0 * f_bet;

        m.nue.x[i][j][k] = m.tke.x[i][j][k] / std::max(w_hat, 1.0e-20);

        const double P_k = std::min(m.prod.x[i][j][k],
            20.0 * bet_star * m.tke.x[i][j][k] * m.dis.x[i][j][k]);
        const double Y_k = bet_star * m.tke.x[i][j][k] * m.dis.x[i][j][k];
        const double P_w = gam * m.dis.x[i][j][k]
                           / std::max(m.tke.x[i][j][k], 1.0e-20) * P_k;
        const double Y_w = bet_wc * m.dis.x[i][j][k] * m.dis.x[i][j][k];
        const double D_w = sig_d / std::max(m.dis.x[i][j][k], 1.0e-20)
            * (dtkedr_neu * ddisdr_neu + dtkedthe_neu * ddisdthe_neu + dtkedphi_neu * ddisdphi_neu);

        m.tke_source.x[i][j][k] = P_k - Y_k;
        m.dis_source.x[i][j][k] = P_w - Y_w + D_w;
    }

    // -----------------------------------------------------------------------
    // z_floor: physical height of the seafloor [m]
    void compute_k_omega_SST(int i, int j, int k,
        double z_floor, double rm, double sinthe,
        double dtkedr, double ddisdr,
        double dtkedthe, double ddisdthe,
        double dtkedphi, double ddisdphi,
        double Omega, double /*W*/)
    {
        const double a1     = 0.31;
        const double bet_star = 0.09;
        const double bet1   = 0.0333;
        const double bet2   = 0.0368;
        const double gam1   = 0.413;
        const double gam2   = 0.2;
        const double sig_k1 = 1.176;
        const double sig_k2 = 1.0;
        const double sig_w1 = 2.0;
        const double sig_w2 = 1.168;

        const double nue_water_nd = nue_water / (m.u_0 * m.L_hyd);

        const double y_phys = std::max(rm * m.L_hyd - z_floor, 1.0e-6);
        const double y_star = y_phys / m.L_hyd;

        const double rmsinthe = rm * sinthe;
        const double exp_rm   = 1.0 / (rm + 1.0);

        const bool neumann_bot = (i > 0      && is_land(m.h, i-1, j, k));
        const bool neumann_jm1 = (j > 0      && is_land(m.h, i, j-1, k));
        const bool neumann_jp1 = (j < m.jm-1 && is_land(m.h, i, j+1, k));
        const bool neumann_km1 = (k > 0      && is_land(m.h, i, j, k-1));
        const bool neumann_kp1 = (k < m.km-1 && is_land(m.h, i, j, k+1));

        const double tke_im1 = neumann_bot ? m.tke.x[i][j][k] : m.tke.x[i-1][j][k];
        const double dis_im1 = neumann_bot ? m.dis.x[i][j][k] : m.dis.x[i-1][j][k];
        const double tke_jm1 = neumann_jm1 ? m.tke.x[i][j][k] : m.tke.x[i][j-1][k];
        const double tke_jp1 = neumann_jp1 ? m.tke.x[i][j][k] : m.tke.x[i][j+1][k];
        const double dis_jm1 = neumann_jm1 ? m.dis.x[i][j][k] : m.dis.x[i][j-1][k];
        const double dis_jp1 = neumann_jp1 ? m.dis.x[i][j][k] : m.dis.x[i][j+1][k];
        const double tke_km1 = neumann_km1 ? m.tke.x[i][j][k] : m.tke.x[i][j][k-1];
        const double tke_kp1 = neumann_kp1 ? m.tke.x[i][j][k] : m.tke.x[i][j][k+1];
        const double dis_km1 = neumann_km1 ? m.dis.x[i][j][k] : m.dis.x[i][j][k-1];
        const double dis_kp1 = neumann_kp1 ? m.dis.x[i][j][k] : m.dis.x[i][j][k+1];

        const double dtkedr_neu  = neumann_bot
            ? (m.tke.x[i+1][j][k] - m.tke.x[i][j][k]) / (m.rad.z[i+1] - m.rad.z[i]) * exp_rm
            : dtkedr;
        const double ddisdr_neu  = neumann_bot
            ? (m.dis.x[i+1][j][k] - m.dis.x[i][j][k]) / (m.rad.z[i+1] - m.rad.z[i]) * exp_rm
            : ddisdr;
        const double dtkedthe_neu = (tke_jp1 - tke_jm1) / (rm       * 2.0 * m.dthe);
        const double ddisdthe_neu = (dis_jp1 - dis_jm1) / (rm       * 2.0 * m.dthe);
        const double dtkedphi_neu = (tke_kp1 - tke_km1) / (rmsinthe * 2.0 * m.dphi);
        const double ddisdphi_neu = (dis_kp1 - dis_km1) / (rmsinthe * 2.0 * m.dphi);

        const double CD_kw = std::max(
            2.0 * sig_w2 / std::max(m.dis.x[i][j][k], 1.0e-20)
            * (dtkedr_neu * ddisdr_neu + dtkedthe_neu * ddisdthe_neu + dtkedphi_neu * ddisdphi_neu),
            1.0e-20);

        const double dis_safe = std::max(m.dis.x[i][j][k], 1.0e-20);
        const double arg1 = std::min(
            std::max(std::sqrt(m.tke.x[i][j][k]) / (bet_star * dis_safe * y_star),
                     500.0 * nue_water_nd / (y_star * y_star * dis_safe)),
            4.0 * sig_w2 * m.tke.x[i][j][k] / (CD_kw * y_star * y_star));
        const double arg2 = std::max(
            2.0 * std::sqrt(m.tke.x[i][j][k]) / (bet_star * dis_safe * y_star),
            500.0 * nue_water_nd / (y_star * y_star * dis_safe));

        const double F1 = tanh(std::pow(arg1, 4));
        const double F2 = tanh(std::pow(arg2, 2));

        m.nue.x[i][j][k] = a1 * m.tke.x[i][j][k]
            / std::max(a1 * dis_safe, Omega * F2);

        const double tke_safe = std::max(m.tke.x[i][j][k], 1.0e-20);
        const double P_k = std::min(m.prod.x[i][j][k],
            20.0 * bet_star * tke_safe * dis_safe);
        const double Y_k = bet_star * tke_safe * dis_safe;
        const double P_w = blend(gam1, gam2, F1) * P_k * dis_safe / tke_safe;
        const double Y_w = blend(bet1, bet2, F1) * dis_safe * dis_safe;
        const double D_w = 2.0 * (1.0 - F1) * sig_w2 / dis_safe
            * (dtkedr_neu * ddisdr_neu + dtkedthe_neu * ddisdthe_neu + dtkedphi_neu * ddisdphi_neu);

        m.tke_source.x[i][j][k] = P_k - Y_k;
        m.dis_source.x[i][j][k] = P_w - Y_w + D_w;
    }
};
