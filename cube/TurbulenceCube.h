/*
 * Atmosphere General Circulation Model (AGCM) applied to turbulent flow
 * Standalone turbulence class: k-epsilon (Chien 1982),
 *                              k-omega   (Wilcox 1988),
 *                              k-omega SST (Menter 1994)
 * Turbulence arrays (tke, dis, nue, prod, tke_source, dis_source) are owned
 * by cCubeModel and accessed here through the model reference m.
*/

#pragma once

#include "cCubeModel.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>

using namespace AtomUtils;

class TurbulenceCube {
public:
    enum Model { k_epsilon, k_omega, k_omega_SST };

    explicit TurbulenceCube(cCubeModel& model,
                            double           vel_star_ref = 0.4,        // reference friction velocity for re_turb [m/s]
                            double           z_0          = 0.1,        // roughness length    [m]
                            double           nue_air      = 1.5e-5)     // kin. viscosity air  [m²/s]
        : m(model),
          turb_model(parse_model(model.turb_model)),
          vel_star_ref(vel_star_ref),
          z_0(z_0),
          nue_air(nue_air),
          re_turb(vel_star_ref * z_0 / nue_air)
    {}

    // -----------------------------------------------------------------------
    void init() {
        using namespace std;
/*
        cout << endl << "      AGCM: TurbulenceAtm::init" << endl;
        print_model_name();

        auto begin = chrono::high_resolution_clock::now();
*/
        m.re_turb = re_turb;
        compute_vel_star();
        init_fields();
        apply_wall_bc();
        compute_sources();   // prime nue, tke_source, dis_source from the ABL profile
        zero_land_cells();
        clamp_nue();
/*
        auto end     = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for TurbulenceAtm::init\n",
               elapsed.count() * 1e-9);
        cout << "      AGCM: TurbulenceAtm::init ended" << endl;
*/
    }

    // -----------------------------------------------------------------------
    void run() {
        using namespace std;
/*
        cout << endl << "      AGCM: TurbulenceAtm::run" << endl;
        print_model_name();

        auto begin = chrono::high_resolution_clock::now();
*/
        compute_vel_star();
        compute_sources();
        zero_land_cells();
        clamp_nue();
/*
        auto end     = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for TurbulenceAtm::run\n",
               elapsed.count() * 1e-9);
        cout << "      AGCM: TurbulenceAtm::run ended" << endl;
*/
    }

    // -----------------------------------------------------------------------
    // Compute the friction velocity u_τ per (j,k) column from the horizontal
    // wind speed at the first air cell above the topographic surface.
    //   U_horiz = u_0 · sqrt(v² + w²)  at i = i_mount + 1
    //   u_τ     = max(vel_star_min, sqrt(C_D) · U_horiz)
    // C_D = 0.002 (neutral drag coefficient); vel_star_min = 0.05 m/s so
    // that calm-wind columns still receive a small but non-zero ABL seed.
    void compute_vel_star() {
        constexpr double C_D          = 0.002;
        constexpr double vel_star_min = 0.05;   // [m/s]

        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                const int i_surf = std::min(m.i_topography[j][k] + 1, m.im - 1);
                const double v_s = m.v.x[i_surf][j][k];
                const double w_s = m.w.x[i_surf][j][k];
                const double U_horiz = m.u_0 * std::sqrt(v_s * v_s + w_s * w_s);
                m.vel_star.y[j][k] = std::max(vel_star_min, std::sqrt(C_D) * U_horiz);
            }
        }
    }

    // -----------------------------------------------------------------------
    // Wall boundary condition for ω at i = 0 (ground surface).
    // Physical formulas produce dimensional ω [s⁻¹]; normalise to ω* = ω · L_atm/u_0.
    // Public so run_3D_loop can reassert the wall values after bcRadius
    // overwrites tke[0] / dis[0] with cubic Neumann extrapolation.
    void apply_wall_bc() {
        const double y1      = m.get_layer_height(1);                   // first-cell height above surface [m]
        const double y1_safe = std::max(y1, 1.0e-6);
        const double nd_omega = m.L_atm / m.u_0;                        // ω_phys [s⁻¹] → ω*

        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                const double vs = m.vel_star.y[j][k];

                // tke wall BC — log-layer equilibrium k = u_τ²/√C_wall.
                // k-ε uses C_μ = C_nue; k-ω and k-ω SST use β* = 0.09
                // (the destruction coefficient that sets the wall TKE level for ω models).
                const double c_wall = (turb_model == k_epsilon) ? C_nue : 0.09;
                m.tke.x[0][j][k] = vs * vs / (std::sqrt(c_wall) * m.u_0 * m.u_0);

                // nue wall BC — turbulent viscosity vanishes at the no-slip wall
                // (viscous sublayer: ν_T → 0).  bcRadius cubic extrapolation can
                // push nue[0] negative; clamp it to 0 here unconditionally.
                m.nue.x[0][j][k] = 0.0;

                if (turb_model == k_epsilon) {
                    // Chien (1982) low-Re k-ε: ε_wall = 0 (Dirichlet).
                    // bcRadius cubic extrapolation can overshoot to negative ε;
                    // the no-slip condition requires ε→0 as y→0.
                    m.dis.x[0][j][k] = 0.0;
                } else if (turb_model == k_omega) {
                    // Wilcox 1988: ω_wall = 6ν / (β₁ y₁²),  β₁ = 0.075
                    const double om_phys = 6.0 * nue_air / (0.075 * y1_safe * y1_safe);
                    m.dis.x[0][j][k] = om_phys * nd_omega;
                } else if (turb_model == k_omega_SST) {
                    // Menter 1994: ω_wall = 60ν / (β₁ y₁²),  β₁ = 0.0333
                    const double om_phys = 60.0 * nue_air / (0.0333 * y1_safe * y1_safe);
                    m.dis.x[0][j][k] = om_phys * nd_omega;
                }
            }
        }
    }

private:
    cCubeModel& m;
    Model  turb_model;
    double vel_star_ref;   // reference value for re_turb; per-column u_τ lives in m.vel_star
    double z_0;
    double nue_air;
    double re_turb;    // = vel_star · z_0 / nue_air  (computed, not a free parameter)

    // ABL model constants
    static constexpr double C_nue   = 0.028;
    static constexpr double Karman  = 0.42;
    // Safety floors / caps (dimensionless)
    static constexpr double dis_min = 1.0e-10;                          // minimum dis* to avoid nue → ∞

    // Map the string parameter m.turb_model to the internal enum
    static Model parse_model(const std::string& s) {
        if (s == "k_epsilon")   return k_epsilon;
        if (s == "k_omega")     return k_omega;
        if (s == "k_omega_SST") return k_omega_SST;
        return k_omega_SST;                                             // default / "none" fallback
    }

    // SST blending: inner (zone 1, near wall) and outer (zone 2, free stream)
    static double blend(double inner, double outer, double F1) {
        return F1 * inner + (1.0 - F1) * outer;
    }

    // -----------------------------------------------------------------------
    void print_model_name() const {
        using namespace std;
        switch(turb_model) {
            case k_epsilon:   cout << "      k-epsilon turbulence model (Chien 1982)" << endl; break;
            case k_omega:     cout << "      k-omega turbulence model (Wilcox 1988)" << endl;  break;
            case k_omega_SST: cout << "      k-omega SST turbulence model (Menter 1994)" << endl; break;
        }
    }

    // -----------------------------------------------------------------------
    // Initialise tke and dis throughout the domain (ABL profile, Narjisse).
    // All source quantities (vel_star, get_layer_height) are dimensional;
    // results are normalised to the model's dimensionless convention:
    //   k*  = k_phys / u_0²
    //   ω*  = ω_phys · L_atm / u_0        (k-ω, k-ω SST)
    //   ε*  = ε_phys · L_atm / u_0³       (k-ε)
    void init_fields() {
        const double inv_u02  = 1.0 / (m.u_0 * m.u_0);
        const double nd_omega = m.L_atm / m.u_0;                        // ω_phys → ω*
        const double nd_eps   = m.L_atm / (m.u_0 * m.u_0 * m.u_0);      // ε_phys → ε*

        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 1; j < m.jm - 1; j++) {
            for (int k = 1; k < m.km - 1; k++) {
                // Per-column friction velocity and derived background floors.
                // k_bg: 1% of peak surface TKE — prevents k→0 and ω/ε→0 at altitude.
                const double vs    = m.vel_star.y[j][k];
                const double k_bg  = 0.01 * vs * vs / sqrt(C_nue);      // [m²/s²]
                const double om_bg = pow(C_nue, -0.25) * sqrt(k_bg) / (Karman * m.L_atm); // [s⁻¹]
                const double eps_bg = pow(C_nue, 0.75) * pow(k_bg, 1.5) / (Karman * m.L_atm); // [m²/s³]

                for (int i = 1; i < m.im - 1; i++) {
                    const double z_i    = m.get_layer_height(i);        // [m]
                    const double z_safe = std::max(z_i, 1.0e-6);

                    // Clamp z_frac to [0,1] so the parabolic profile never
                    // exceeds 1 above the ABL height L_atm (prevents huge TKE aloft).
                    const double z_frac = std::min(z_i, double(m.L_atm)) / m.L_atm;

                    // Physical TKE [m²/s²]: ABL parabolic profile with background floor
                    const double k_abl  = vs * vs / sqrt(C_nue)
                        * pow(1.0 - z_frac, 2);
                    const double k_phys = std::max(k_abl, k_bg);

                    double dis_nd = 0.0;
                    if (turb_model == k_epsilon) {
                        // Physical ε [m²/s³]: ε = C_μ^0.75 · k^1.5 / (κ · z)
                        const double eps_phys = pow(C_nue, 0.75)
                            * pow(k_phys, 1.5) / (Karman * z_safe);
                        dis_nd = std::max(eps_phys, eps_bg) * nd_eps;
                    } else {
                        // Physical ω [s⁻¹]: ω = C_μ^-0.25 · k^0.5 / (κ · z)
                        // Floor at om_bg to prevent ω→0 aloft driving nue=k/ω → ∞.
                        const double om_phys = pow(C_nue, -0.25)
                            * sqrt(k_phys) / (Karman * z_safe);
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
    // Global nue clamp: enforce 0 ≤ nue* ≤ nue_max everywhere, including pole
    // rows (j=0, j=jm-1) that are not touched by compute_sources and may receive
    // negative values from cubic BC extrapolation in bcTheta.
    void clamp_nue() {
        const double nue_max = 1000.0 / (m.u_0 * m.L_atm);
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
    // Zero all turbulence quantities inside land cells
    void zero_land_cells() {
        #pragma omp parallel for collapse(3) schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                for (int k = 0; k < m.km; k++) {
                    if (is_land(m.h, i, j, k)) {
                        m.tke.x[i][j][k]        = 0.0;
                        m.tken.x[i][j][k]       = 0.0;
                        m.dis.x[i][j][k]        = 0.0;
                        m.disn.x[i][j][k]       = 0.0;
                        m.nue.x[i][j][k]        = 0.0;
                        m.prod.x[i][j][k]       = 0.0;
                        m.tke_source.x[i][j][k] = 0.0;
                        m.dis_source.x[i][j][k] = 0.0;
                    }
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // Main turbulence source computation
    void compute_sources() {
        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 1; j < m.jm - 1; j++) {
            for (int k = 1; k < m.km - 1; k++) {
                const int    i_mount    = m.i_topography[j][k];
                const double y_mount    = m.get_layer_height(i_mount);  // surface physical height [m]

                // Background tke floor for this column [dimensionless k* = k/u_0²].
                // Prevents tke from decaying below the ABL seed level in low-shear
                // regions where P_k is insufficient to offset destruction.
                const double vs      = m.vel_star.y[j][k];
                const double k_bg_nd = 0.01 * vs * vs / (std::sqrt(C_nue) * m.u_0 * m.u_0);

                for (int i = 1; i < m.im - 1; i++) {
                    if (is_land(m.h, i, j, k)) continue;

                    // ---- geometry ----
                    const double inv_2dz = 1.0 / (2.0 * m.dr);
                    const double inv_2dy = 1.0 / (2.0 * m.dy);
                    const double inv_2dx = 1.0 / (2.0 * m.dx);

                    // Enforce floors on this cell before any reads — the RK4 can push
                    // tke or dis below their safe limits between compute_sources calls.
                    m.tke.x[i][j][k] = std::max(0.0,     m.tke.x[i][j][k]);
                    m.dis.x[i][j][k] = std::max(dis_min, m.dis.x[i][j][k]);

                    // ---- velocity gradients ----
                    const double dudz   = (m.u.x[i+1][j][k] - m.u.x[i-1][j][k]) * inv_2dz;
                    const double dvdz   = (m.v.x[i+1][j][k] - m.v.x[i-1][j][k]) * inv_2dz;
                    const double dwdz   = (m.w.x[i+1][j][k] - m.w.x[i-1][j][k]) * inv_2dz;

                    const double dudy   = (m.u.x[i][j+1][k] - m.u.x[i][j-1][k]) * inv_2dy;
                    const double dvdy   = (m.v.x[i][j+1][k] - m.v.x[i][j-1][k]) * inv_2dy;
                    const double dwdy   = (m.w.x[i][j+1][k] - m.w.x[i][j-1][k]) * inv_2dy;

                    const double dudx   = (m.u.x[i][j][k+1] - m.u.x[i][j][k-1]) * inv_2dx;
                    const double dvdx   = (m.v.x[i][j][k+1] - m.v.x[i][j][k-1]) * inv_2dx;
                    const double dwdx   = (m.w.x[i][j][k+1] - m.w.x[i][j][k-1]) * inv_2dx;

                    // ---- tke / dis gradients with Neumann BC at land faces ----
                    // Replace each land neighbour's value with the current cell's value
                    // (zero-gradient) so stencils do not read the zeroed land-cell values.
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

                    const double dtkedz   = (m.tke.x[i+1][j][k] - tke_im1) * inv_2dz;
                    const double ddisdz   = (m.dis.x[i+1][j][k] - dis_im1) * inv_2dz;
                    const double dtkedy   = (tke_jp1 - tke_jm1) * inv_2dy;
                    const double ddisdy   = (dis_jp1 - dis_jm1) * inv_2dy;
                    const double dtkedx   = (tke_kp1 - tke_km1) * inv_2dx;
                    const double ddisdx   = (dis_kp1 - dis_km1) * inv_2dx;

                    // ---- production tensor P_k ----
                    // Compute νt fresh from current tke/dis rather than reading m.nue,
                    // which lags one call behind (nue is updated inside compute_k_*
                    // AFTER prod is needed here).  On the first call from init() m.nue
                    // is still 0, giving prod≈0 and breaking the turbulence bootstrap.
                    // k-ε: νt = Cμ·k²/ε   k-ω / k-ω SST: νt ≈ k/ω  (F2 limiter is
                    // a minor correction; the full SST value is recomputed below).
                    double cnue;
                    {
                        // Use dis_min (not 1e-20) so that RK4-driven dips below dis_min
                        // do not produce cnue=k/1e-20 → inf in prod.
                        const double dis_here = std::max(m.dis.x[i][j][k], dis_min);
                        const double tke_here = std::max(m.tke.x[i][j][k], 0.0);
                        const double nue_max  = 1000.0 / (m.u_0 * m.L_atm);
                        cnue = (turb_model == k_epsilon)
                             ? C_nue * tke_here * tke_here / dis_here
                             : tke_here / dis_here;
                        cnue = std::min(cnue, nue_max);                 // guard against large tke spikes
                    }
                    const double der  = 0.66667 * (dudz + dvdy + dwdx);

                    m.prod.x[i][j][k] = std::max(0.0,
                          (cnue * (2.0 * dudz - der) - 0.66667 * m.tke.x[i][j][k]) * dudz
                        + (cnue * (dudy + dvdz))                                   * dudy
                        + (cnue * (dudx + dwdz))                                   * dudx
                        + (cnue * (2.0 * dvdy - der) - 0.66667 * m.tke.x[i][j][k]) * dvdy
                        + (cnue * (dvdz + dudy))                                   * dvdz
                        + (cnue * (dvdx + dwdy))                                   * dvdx
                        + (cnue * (2.0 * dwdx - der) - 0.66667 * m.tke.x[i][j][k]) * dwdx
                        + (cnue * (dwdz + dudx))                                   * dwdz
                        + (cnue * (dwdy + dvdx))                                   * dwdy);

                    // ---- vorticity magnitude Ω = |∇×u| ----
                    // Each cross-derivative pair must be squared and summed — not
                    // added with its negated partner (which would cancel to zero).
                    const double W12 = dudy - dvdz;
                    const double W13 = dudx - dwdz;
                    const double W23 = dvdx - dwdy;
                    const double Omega = sqrt(W12 * W12 + W13 * W13 + W23 * W23);

                    // ================================================================
                    if (turb_model == k_epsilon) {
                        compute_k_epsilon(i, j, k, y_mount,
                            dudz, dvdz, dwdz, dudy, dvdy, dwdy, dudx, dvdx, dwdx,
                            Omega);
                    }
                    // ================================================================
                    else if (turb_model == k_omega) {
                        compute_k_omega(i, j, k,
                            dtkedz, ddisdz, dtkedy, ddisdy, dtkedx, ddisdx,
                            Omega, Omega);
                    }
                    // ================================================================
                    else if (turb_model == k_omega_SST) {
                        compute_k_omega_SST(i, j, k, y_mount,
                            dtkedz, ddisdz, dtkedy, ddisdy, dtkedx, ddisdx,
                            Omega, Omega);
                    }

                    // nue_max: 1000 m²/s physical cap converted to dimensionless ν* = ν/(u_0·L_atm)
                    const double nue_max = 1000.0 / (m.u_0 * m.L_atm);
                    m.tke.x[i][j][k] = std::max(k_bg_nd, m.tke.x[i][j][k]);
                    m.dis.x[i][j][k] = std::max(dis_min, m.dis.x[i][j][k]);
                    m.nue.x[i][j][k] = std::clamp(m.nue.x[i][j][k], 0.0, nue_max);
                    // dis_source_max: 20·ω*² cap — prevents polar blow-up from the 1/sinθ²
                    // amplification in the cross-diffusion term D_w near j→0/jm-1.
                    // Factor 20 gives ~20×β_max ≈ 1.8× max destruction, ample headroom.
                    const double dis_src_max = 20.0 * m.dis.x[i][j][k] * m.dis.x[i][j][k];
                    m.dis_source.x[i][j][k] = std::clamp(m.dis_source.x[i][j][k],
                                                          -dis_src_max, dis_src_max);
                    // tke_source_max: 20·k*·ω* cap — mirrors dis_source logic;
                    // bounds the P_k-Y_k balance without suppressing legitimate physics.
                    const double tke_src_max = 20.0 * m.tke.x[i][j][k] * m.dis.x[i][j][k];
                    m.tke_source.x[i][j][k] = std::clamp(m.tke_source.x[i][j][k],
                                                          -tke_src_max, tke_src_max);
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // y_mount: physical height of the topographic surface [m] (get_layer_height(i_mount))
    void compute_k_epsilon(int i, int j, int k, double y_mount,
        double dudz, double dvdz, double dwdz,
        double dudy, double dvdy, double dwdy,
        double dudx, double dvdx, double dwdx,
        double /*Omega*/)
    {
        const double C_eps_1 = 1.35;
        const double C_eps_2 = 1.80;

        // Physical wall distance [m] and its dimensionless form y* = y/L_atm.
        // Floor at L_atm (one grid cell) so that y_star ≥ 1 and the 1/y_star² in
        // the Chien D term stays bounded for cells straddling mountain tops.
        const double y_phys = std::max(m.get_layer_height(i) - y_mount, double(m.L_atm));  // [m]
        const double y_star = y_phys / m.L_atm;  // dimensionless, always ≥ 1

        // ν_T* = C_μ · k*² / ε*   (already dimensionless — see header)
        m.nue.x[i][j][k] = C_nue * m.tke.x[i][j][k] * m.tke.x[i][j][k]
                          / std::max(m.dis.x[i][j][k], dis_min);

        // y⁺ = y_phys · u_τ / ν_air  (viscous wall units, dimensionless)
        const double d_plus = y_phys * m.vel_star.y[j][k] / nue_air;
        const double f_nue  = 1.0 - exp(-0.0115 * d_plus);

        // Turbulent Reynolds number Re_T = k*² · u_0 · L_atm / (ε* · ν_air)
        // (uses molecular viscosity, not eddy viscosity)
        const double Re_T = m.tke.x[i][j][k] * m.tke.x[i][j][k]
                            * m.u_0 * m.L_atm
                            / std::max(m.dis.x[i][j][k] * nue_air, dis_min * nue_air);
        const double f_2  = 1.0 - 0.4 / 1.8 * exp(-Re_T * Re_T / 36.0);

        m.nue.x[i][j][k] *= f_nue;

        // Chien (1982) viscous near-wall damping terms D and E.
        // y_star ≥ 1 (enforced above), so y_star² ≥ 1 and L_k is bounded.
        // For ABL cells (d_plus ≈ 10⁷), exp(-d_plus/2) ≈ 0 so L_w ≈ 0.
        const double y_star2 = y_star * y_star;
        const double L_k = -2.0 * m.tke.x[i][j][k] / y_star2;
        const double L_w = -2.0 * m.dis.x[i][j][k] / y_star2 * exp(-d_plus / 2.0);

        const double P_k = m.prod.x[i][j][k];
        const double Y_k = m.dis.x[i][j][k];
        // Standard k-ε: P_ε = C_ε1 · (ε/k) · P_k  — NOT k·ε·P_k (wrong by k²).
        const double tke_safe = std::max(m.tke.x[i][j][k], dis_min);
        const double P_w = C_eps_1 * m.dis.x[i][j][k] / tke_safe * P_k;
        const double Y_w = C_eps_2 * f_2 * m.dis.x[i][j][k] * m.dis.x[i][j][k]
                           / tke_safe;

        // Chien viscous correction: D_k* = −2 ν* k*/y*², factor = ν* = ν/(u₀ L_atm)
        const double nu_star = nue_air / (m.u_0 * m.L_atm);
        m.tke_source.x[i][j][k] = P_k - Y_k + L_k * nu_star;
        m.dis_source.x[i][j][k] = P_w - Y_w + L_w * nu_star;
    }

    // -----------------------------------------------------------------------
    void compute_k_omega(int i, int j, int k,
        double dtkedz, double ddisdz,
        double dtkedy, double ddisdy,
        double dtkedx, double ddisdx,
        double Omega_mag, double /*W_unused*/)
    {
        const double bet_star = 0.09;
        const double gam      = 0.52;
        const double C_lim    = 0.875;
        const double bet_0    = 0.0708;

        // Neumann BC at all six faces
        const bool neumann_bot = (i > 0      && is_land(m.h, i-1, j, k));
        const bool neumann_jm1 = (j > 0      && is_land(m.h, i, j-1, k));
        const bool neumann_jp1 = (j < m.jm-1 && is_land(m.h, i, j+1, k));
        const bool neumann_km1 = (k > 0      && is_land(m.h, i, j, k-1));
        const bool neumann_kp1 = (k < m.km-1 && is_land(m.h, i, j, k+1));

        const double tke_jm1 = neumann_jm1 ? m.tke.x[i][j][k] : m.tke.x[i][j-1][k];
        const double tke_jp1 = neumann_jp1 ? m.tke.x[i][j][k] : m.tke.x[i][j+1][k];
        const double dis_jm1 = neumann_jm1 ? m.dis.x[i][j][k] : m.dis.x[i][j-1][k];
        const double dis_jp1 = neumann_jp1 ? m.dis.x[i][j][k] : m.dis.x[i][j+1][k];
        const double tke_km1 = neumann_km1 ? m.tke.x[i][j][k] : m.tke.x[i][j][k-1];
        const double tke_kp1 = neumann_kp1 ? m.tke.x[i][j][k] : m.tke.x[i][j][k+1];
        const double dis_km1 = neumann_km1 ? m.dis.x[i][j][k] : m.dis.x[i][j][k-1];
        const double dis_kp1 = neumann_kp1 ? m.dis.x[i][j][k] : m.dis.x[i][j][k+1];

        const double dtkedz_neu = neumann_bot
            ? (m.tke.x[i+1][j][k] - m.tke.x[i][j][k]) / (2.0 * m.dr)
            : dtkedz;
        const double ddisdz_neu = neumann_bot
            ? (m.dis.x[i+1][j][k] - m.dis.x[i][j][k]) / (2.0 * m.dr)
            : ddisdz;
        const double dtkedy_neu = (tke_jp1 - tke_jm1) / (2.0 * m.dy);
        const double ddisdy_neu = (dis_jp1 - dis_jm1) / (2.0 * m.dy);
        const double dtkedx_neu = (tke_kp1 - tke_km1) / (2.0 * m.dx);
        const double ddisdx_neu = (dis_kp1 - dis_km1) / (2.0 * m.dx);

        const double dudz   = (m.u.x[i+1][j][k] - m.u.x[i-1][j][k]) / (2.0 * m.dr);
        const double dvdy   = (m.v.x[i][j+1][k] - m.v.x[i][j-1][k]) / (2.0 * m.dy);
        const double dwdx   = (m.w.x[i][j][k+1] - m.w.x[i][j][k-1]) / (2.0 * m.dx);

        // Strain-rate magnitude S = sqrt(2 S_ij S_ij)
        const double dudy_u  = (m.u.x[i][j+1][k] - m.u.x[i][j-1][k]) / (2.0 * m.dy);
        const double dudx_u  = (m.u.x[i][j][k+1] - m.u.x[i][j][k-1]) / (2.0 * m.dx);
        const double dvdz_u  = (m.v.x[i+1][j][k] - m.v.x[i-1][j][k]) / (2.0 * m.dr);
        const double dwdz_u  = (m.w.x[i+1][j][k] - m.w.x[i-1][j][k]) / (2.0 * m.dr);
        const double dvdx_u  = (m.v.x[i][j][k+1] - m.v.x[i][j][k-1]) / (2.0 * m.dx);
        const double dwdy_u  = (m.w.x[i][j+1][k] - m.w.x[i][j-1][k]) / (2.0 * m.dy);

        const double S11 = dudz,  S22 = dvdy, S33 = dwdx;
        const double S12 = 0.5 * (dudy_u + dvdz_u);
        const double S13 = 0.5 * (dudx_u + dwdz_u);
        const double S23 = 0.5 * (dvdx_u + dwdy_u);
        const double S_mag = sqrt(2.0 * (S11*S11 + S22*S22 + S33*S33
                                       + 2.0*(S12*S12 + S13*S13 + S23*S23)));

        const double w_hat = std::max(m.dis.x[i][j][k],
            C_lim * S_mag / sqrt(bet_star));

        double sig_d = dtkedz_neu * ddisdz_neu + dtkedy_neu * ddisdy_neu + dtkedx_neu * ddisdx_neu;
        sig_d = (sig_d <= 0.0) ? 0.0 : 0.125;

        const double chi_w  = fabs(Omega_mag * Omega_mag * S_mag
            / pow(bet_star * std::max(m.dis.x[i][j][k], 1.0e-20), 3));
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
            * (dtkedz_neu * ddisdz_neu + dtkedy_neu * ddisdy_neu + dtkedx_neu * ddisdx_neu);

        m.tke_source.x[i][j][k] = P_k - Y_k;
        m.dis_source.x[i][j][k] = P_w - Y_w + D_w;
    }

    // -----------------------------------------------------------------------
    // y_mount: physical height of the topographic surface [m] (get_layer_height(i_mount))
    void compute_k_omega_SST(int i, int j, int k,
        double y_mount,
        double dtkedz, double ddisdz,
        double dtkedy, double ddisdy,
        double dtkedx, double ddisdx,
        double Omega, double /*W*/)
    {
        // Menter 1994 ABL constants
        const double a1       = 0.31;
        const double bet_star = 0.09;
        const double bet1     = 0.0333;
        const double bet2     = 0.0368;
        const double gam1     = 0.413;
        const double gam2     = 0.2;
        const double sig_w2   = 1.168;

        const double nue_air_nd = nue_air / (m.u_0 * m.L_atm);

        const double y_phys = std::max(m.get_layer_height(i) - y_mount, 1.0e-6);
        const double y_star = y_phys / m.L_atm;

        // Neumann BC at all six faces
        const bool neumann_bot = (i > 0      && is_land(m.h, i-1, j, k));
        const bool neumann_jm1 = (j > 0      && is_land(m.h, i, j-1, k));
        const bool neumann_jp1 = (j < m.jm-1 && is_land(m.h, i, j+1, k));
        const bool neumann_km1 = (k > 0      && is_land(m.h, i, j, k-1));
        const bool neumann_kp1 = (k < m.km-1 && is_land(m.h, i, j, k+1));

        const double tke_jm1 = neumann_jm1 ? m.tke.x[i][j][k] : m.tke.x[i][j-1][k];
        const double tke_jp1 = neumann_jp1 ? m.tke.x[i][j][k] : m.tke.x[i][j+1][k];
        const double dis_jm1 = neumann_jm1 ? m.dis.x[i][j][k] : m.dis.x[i][j-1][k];
        const double dis_jp1 = neumann_jp1 ? m.dis.x[i][j][k] : m.dis.x[i][j+1][k];
        const double tke_km1 = neumann_km1 ? m.tke.x[i][j][k] : m.tke.x[i][j][k-1];
        const double tke_kp1 = neumann_kp1 ? m.tke.x[i][j][k] : m.tke.x[i][j][k+1];
        const double dis_km1 = neumann_km1 ? m.dis.x[i][j][k] : m.dis.x[i][j][k-1];
        const double dis_kp1 = neumann_kp1 ? m.dis.x[i][j][k] : m.dis.x[i][j][k+1];

        const double dtkedz_neu = neumann_bot
            ? (m.tke.x[i+1][j][k] - m.tke.x[i][j][k]) / (2.0 * m.dr)
            : dtkedz;
        const double ddisdz_neu = neumann_bot
            ? (m.dis.x[i+1][j][k] - m.dis.x[i][j][k]) / (2.0 * m.dr)
            : ddisdz;
        const double dtkedy_neu = (tke_jp1 - tke_jm1) / (2.0 * m.dy);
        const double ddisdy_neu = (dis_jp1 - dis_jm1) / (2.0 * m.dy);
        const double dtkedx_neu = (tke_kp1 - tke_km1) / (2.0 * m.dx);
        const double ddisdx_neu = (dis_kp1 - dis_km1) / (2.0 * m.dx);

        const double CD_kw = std::max(
            2.0 * sig_w2 / std::max(m.dis.x[i][j][k], 1.0e-20)
            * (dtkedz_neu * ddisdz_neu + dtkedy_neu * ddisdy_neu + dtkedx_neu * ddisdx_neu),
            1.0e-20);

        const double dis_safe = std::max(m.dis.x[i][j][k], 1.0e-20);
        const double arg1 = std::min(
            std::max(sqrt(m.tke.x[i][j][k]) / (bet_star * dis_safe * y_star),
                     500.0 * nue_air_nd / (y_star * y_star * dis_safe)),
            4.0 * sig_w2 * m.tke.x[i][j][k] / (CD_kw * y_star * y_star));
        const double arg2 = std::max(
            2.0 * sqrt(m.tke.x[i][j][k]) / (bet_star * dis_safe * y_star),
            500.0 * nue_air_nd / (y_star * y_star * dis_safe));

        const double F1 = tanh(pow(arg1, 4));
        const double F2 = tanh(pow(arg2, 2));

        m.nue.x[i][j][k] = a1 * m.tke.x[i][j][k]
            / std::max(a1 * dis_safe, Omega * F2);

        const double tke_safe = std::max(m.tke.x[i][j][k], 1.0e-20);
        const double P_k = std::min(m.prod.x[i][j][k],
            20.0 * bet_star * tke_safe * dis_safe);
        const double Y_k = bet_star * tke_safe * dis_safe;
        const double P_w = blend(gam1, gam2, F1) * P_k * dis_safe / tke_safe;
        const double Y_w = blend(bet1, bet2, F1) * dis_safe * dis_safe;
        const double D_w = 2.0 * (1.0 - F1) * sig_w2 / dis_safe
            * (dtkedz_neu * ddisdz_neu + dtkedy_neu * ddisdy_neu + dtkedx_neu * ddisdx_neu);

        m.tke_source.x[i][j][k] = P_k - Y_k;
        m.dis_source.x[i][j][k] = P_w - Y_w + D_w;
    }
};
