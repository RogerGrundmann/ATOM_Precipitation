#pragma once

#include "cHydrosphereModel.h"
#include "Utils.h"

#include <vector>
#include <cstdint>
#include <cmath>
#include <chrono>
#include <cstdio>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace AtomUtils;

// ============================================================================
// PressureSolverHyd — friend class of cHydrosphereModel
//
// Solves the Poisson equation for dynamic pressure on the hydrosphere grid
// using a Gauss-Seidel-style sweep with one-sided (land-contact) and
// central-difference divergence terms.
//
// Usage:
//   PressureSolverHyd(*this).run();
// ============================================================================
class PressureSolverHyd {
public:
    explicit PressureSolverHyd(cHydrosphereModel& model)
        : m(model)
    {}

    void run(int n_sweeps = 1, bool verbose = true)
    {
        using namespace std;
        if (verbose) cout << endl << endl << endl
            << "      OGCM: PressureSolverHyd (" << n_sweeps << " Jacobi sweeps)" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        // ====================================================================
        // Precompute sin(the) table — avoids redundant sin() calls inside loops
        // ====================================================================
        std::vector<double> sinthe_table(m.jm);
        for (int j = 0; j < m.jm; j++) {
            sinthe_table[j] = sin(m.the.z[j]);
            if (sinthe_table[j] < 0.4)  sinthe_table[j] = 0.4;
        }

        // ====================================================================
        // Flat land mask — eliminates repeated is_land() call overhead
        // ====================================================================
        std::vector<int8_t> land(m.im * m.jm * m.km);
        #pragma omp parallel for collapse(3)
        for (int i = 0; i < m.im; i++)
            for (int j = 0; j < m.jm; j++)
                for (int k = 0; k < m.km; k++)
                    land[i * m.jm * m.km + j * m.km + k] =
                        is_land(m.h, i, j, k) ? 1 : 0;

        #define LAND(i,j,k) land[(i) * m.jm * m.km + (j) * m.km + (k)]

        // ====================================================================
        // Radial boundary conditions for aux velocity fields
        // ====================================================================
        #pragma omp parallel for collapse(2)
        for (int j = 1; j < m.jm-1; j++) {
            for (int k = 1; k < m.km-1; k++) {
                m.aux_u.x[0][j][k]      = m.c43 * m.aux_u.x[1][j][k]
                                          - m.c13 * m.aux_u.x[2][j][k];
                m.aux_u.x[m.im-1][j][k] = m.c43 * m.aux_u.x[m.im-2][j][k]
                                          - m.c13 * m.aux_u.x[m.im-3][j][k];
                m.aux_v.x[0][j][k]      = m.c43 * m.aux_v.x[1][j][k]
                                          - m.c13 * m.aux_v.x[2][j][k];
                m.aux_v.x[m.im-1][j][k] = m.c43 * m.aux_v.x[m.im-2][j][k]
                                          - m.c13 * m.aux_v.x[m.im-3][j][k];
                m.aux_w.x[0][j][k]      = m.c43 * m.aux_w.x[1][j][k]
                                          - m.c13 * m.aux_w.x[2][j][k];
                m.aux_w.x[m.im-1][j][k] = m.c43 * m.aux_w.x[m.im-2][j][k]
                                          - m.c13 * m.aux_w.x[m.im-3][j][k];
            }
        }

        // ====================================================================
        // Grid-spacing reciprocals — constant for the entire grid
        // ====================================================================
        const double inv_2dr   = 1.0 / (2.0 * m.dr);
        const double inv_2dthe = 1.0 / (2.0 * m.dthe);
        const double inv_2dphi = 1.0 / (2.0 * m.dphi);
        const double inv_dr2   = 1.0 / (m.dr   * m.dr);
        const double inv_dthe2 = 1.0 / (m.dthe * m.dthe);
        const double inv_dphi2 = 1.0 / (m.dphi * m.dphi);
        const double inv_dthe  = 1.0 / m.dthe;
        const double inv_dphi  = 1.0 / m.dphi;

        // ====================================================================
        // Main pressure sweep — iterated n_sweeps times so the pressure-Poisson
        // equation actually converges each call. A single Jacobi sweep barely
        // dents the residual; the hydro solver previously did just ONE sweep per
        // call (and only every 2nd iter), so continuity (div u -> 0) was never
        // enforced and the wind-driven convergence in the tropical Indian Ocean
        // drove a radial-velocity runaway to NaN (~iter 635, see
        // project_hydro_polar_blowup). aux (the divergence source) and the land
        // mask are fixed across sweeps, so only the p_dyn update + its BCs repeat.
        // ====================================================================
        for (int sweep = 0; sweep < n_sweeps; sweep++) {
        #pragma omp parallel for collapse(2) schedule(dynamic, 4)
        for (int i = 1; i < m.im-1; i++) {
            for (int j = 1; j < m.jm-1; j++) {

                // Build geometry struct once per (i,j) — reused for every k
                cHydrosphereModel::CellGeometry geo;

                geo.rm           = m.rad.z[i];
                geo.rm2          = geo.rm * geo.rm;
                geo.exp_rm       = 1.0 / (geo.rm + 1.0);
                geo.exp_2_rm     = geo.exp_rm * geo.exp_rm;
                geo.sinthe       = sinthe_table[j];
                geo.sinthe2      = geo.sinthe * geo.sinthe;
                geo.costhe       = cos(m.the.z[j]);
                geo.inv_rm       = 1.0 / geo.rm;
                geo.inv_rm2      = 1.0 / geo.rm2;
                geo.inv_rmsinthe         = 1.0 / (geo.rm * geo.sinthe);
                geo.inv_rm2sinthe        = geo.inv_rm2 / geo.sinthe;
                geo.inv_rm2sinthe2       = geo.inv_rm2 / geo.sinthe2;
                geo.costhe_inv_rm2sinthe = geo.costhe * geo.inv_rm2sinthe;
                geo.inv_2dr   = inv_2dr;
                geo.inv_2dthe = inv_2dthe;
                geo.inv_2dphi = inv_2dphi;
                geo.inv_dr2   = inv_dr2;
                geo.inv_dthe2 = inv_dthe2;
                geo.inv_dphi2 = inv_dphi2;

                // Radial terms carry the exp_rm metric (exp_2_rm = exp_rm² on the
                // Laplacian, exp_rm on the divergence source below) so this Poisson
                // operator is consistent with the RHS_Hyd momentum pressure gradient
                // (dpdr * exp_rm) and with the working atmosphere solver. Dropping
                // exp_rm here left the solved p_dyn unable to cancel the radial
                // divergence the gradient correction applied -> radial-velocity
                // runaway (see project_hydro_polar_blowup metric-mismatch audit).
                // Radial Laplacian uses the per-i non-uniform 2nd-derivative
                // coefficients (stretched grid); theta/phi grids are uniform so keep
                // the symmetric inv_d*2 form. rc20 = -(rc2m+rc2p) is the radial
                // diagonal contribution.
                const double num1m = geo.exp_2_rm * m.rc2m[i];   // coeff of p[i-1]
                const double num1p = geo.exp_2_rm * m.rc2p[i];   // coeff of p[i+1]
                const double rad_diag = geo.exp_2_rm * (m.rc2m[i] + m.rc2p[i]);
                const double denom = rad_diag
                                   + 2.0 * geo.inv_rm       * inv_dthe2
                                   + 2.0 * geo.inv_rmsinthe * inv_dphi2;
                const double inv_denom = 1.0 / denom;
                const double num2 = geo.inv_rm       * inv_dthe2;
                const double num3 = geo.inv_rmsinthe * inv_dphi2;

                const bool i_in_range = (i < m.im-2);
                const bool j_inner    = (j > 2) && (j < m.jm-2);

                // Sliding window for k-direction land status
                int8_t lnd_k0 = LAND(i,j,0), lnd_k1 = LAND(i,j,1);

                for (int k = 1; k < m.km-1; k++) {

                    const int8_t lnd_ijk = lnd_k1;
                    const int8_t lnd_kp1 = LAND(i,j,k+1);
                    const int8_t lnd_km1 = lnd_k0;
                    const int8_t lnd_kp2 = (k < m.km-2) ? LAND(i,j,k+2) : 0;
                    const int8_t lnd_km2 = (k > 2)      ? LAND(i,j,k-2) : 0;

                    lnd_k0 = lnd_k1;
                    lnd_k1 = lnd_kp1;

                    double du_dr, dv_dthe, dw_dphi;
                    bool r_flag   = false;
                    bool the_flag = false;
                    bool phi_flag = false;

                    // ---- r direction -----------------------------------
                    if (i_in_range && lnd_ijk && !LAND(i+1,j,k)) {
                        du_dr = m.rf10[i] * m.aux_u.x[i][j][k]
                              + m.rf11[i] * m.aux_u.x[i+1][j][k]
                              + m.rf12[i] * m.aux_u.x[i+2][j][k];
                        r_flag = true;
                    }

                    // ---- theta direction --------------------------------
                    if (j_inner) {
                        const int8_t air_jp1 = !LAND(i,j+1,k);
                        const int8_t air_jm1 = !LAND(i,j-1,k);

                        if (lnd_ijk && air_jp1 && !LAND(i,j+2,k)) {
                            dv_dthe = (-3.0 * m.aux_v.x[i][j][k]
                                       + 4.0 * m.aux_v.x[i][j+1][k]
                                       - m.aux_v.x[i][j+2][k]) * inv_2dthe;
                            the_flag = true;
                        }
                        if (lnd_ijk && air_jm1 && !LAND(i,j-2,k)) {
                            dv_dthe = -(-3.0 * m.aux_v.x[i][j][k]
                                        + 4.0 * m.aux_v.x[i][j-1][k]
                                        - m.aux_v.x[i][j-2][k]) * inv_2dthe;
                            the_flag = true;
                        }
                        if ((lnd_ijk && air_jp1 && LAND(i,j+2,k))
                            || (j == m.jm-2 && !lnd_ijk && LAND(i,j+1,k))) {
                            dv_dthe = (m.aux_v.x[i][j+1][k]
                                      - m.aux_v.x[i][j][k]) * inv_dthe;
                            the_flag = true;
                        }
                        if ((lnd_ijk && air_jm1 && LAND(i,j-2,k))
                            || (j == 1 && lnd_ijk && air_jm1)) {
                            dv_dthe = (m.aux_v.x[i][j-1][k]
                                      - m.aux_v.x[i][j][k]) * inv_dthe;
                            the_flag = true;
                        }
                    }

                    // ---- phi direction ----------------------------------
                    const bool k_inner = (k > 2) && (k < m.km-2);
                    if (k_inner) {
                        const int8_t air_kp1 = !lnd_kp1;
                        const int8_t air_km1 = !lnd_km1;

                        if (lnd_ijk && air_kp1 && !lnd_kp2) {
                            dw_dphi = (-3.0 * m.aux_w.x[i][j][k]
                                       + 4.0 * m.aux_w.x[i][j][k+1]
                                       - m.aux_w.x[i][j][k+2]) * inv_2dphi;
                            phi_flag = true;
                        }
                        if (lnd_ijk && air_km1 && !lnd_km2) {
                            dw_dphi = -(-3.0 * m.aux_w.x[i][j][k]
                                        + 4.0 * m.aux_w.x[i][j][k-1]
                                        - m.aux_w.x[i][j][k-2]) * inv_2dphi;
                            phi_flag = true;
                        }
                        if ((lnd_ijk && air_kp1 && lnd_kp2)
                            || (k == m.km-2 && !lnd_ijk && lnd_kp1)) {
                            dw_dphi = (m.aux_w.x[i][j][k+1]
                                      - m.aux_w.x[i][j][k]) * inv_dphi;
                            phi_flag = true;
                        }
                        if ((lnd_ijk && air_km1 && lnd_km2)
                            || (k == 1 && lnd_ijk && air_km1)) {
                            dw_dphi = (m.aux_w.x[i][j][k-1]
                                      - m.aux_w.x[i][j][k]) * inv_dphi;
                            phi_flag = true;
                        }
                    }

                    // ---- Central-difference fallbacks -------------------
                    if (!r_flag)
                        du_dr   = m.rc1m[i] * m.aux_u.x[i-1][j][k]
                                + m.rc10[i] * m.aux_u.x[i][j][k]
                                + m.rc1p[i] * m.aux_u.x[i+1][j][k];
                    if (!the_flag)
                        dv_dthe = (m.aux_v.x[i][j+1][k]
                                  - m.aux_v.x[i][j-1][k]) * inv_2dthe;
                    if (!phi_flag)
                        dw_dphi = (m.aux_w.x[i][j][k+1]
                                  - m.aux_w.x[i][j][k-1]) * inv_2dphi;

                    // ---- Pressure update --------------------------------
                    m.p_dyn.x[i][j][k] =
                        ( num1m * m.p_dyn.x[i-1][j][k] + num1p * m.p_dyn.x[i+1][j][k]
                       + (m.p_dyn.x[i][j+1][k] + m.p_dyn.x[i][j-1][k]) * num2
                       + (m.p_dyn.x[i][j][k+1] + m.p_dyn.x[i][j][k-1]) * num3
                       - du_dr   * geo.exp_rm
                       - dv_dthe * geo.inv_rm
                       - dw_dphi * geo.inv_rmsinthe) * inv_denom;

                }  // k
            }  // j
        }  // i

        // ====================================================================
        // Boundary conditions on p_dyn (re-applied after every sweep)
        // ====================================================================

        // Radial: cubic extrapolation at both ends
        #pragma omp parallel for collapse(2)
        for (int k = 0; k < m.km; k++) {
            for (int j = 0; j < m.jm; j++) {
                m.p_dyn.x[0][j][k]      = m.p_dyn.x[3][j][k]
                    - 3.0 * m.p_dyn.x[2][j][k]
                    + 3.0 * m.p_dyn.x[1][j][k];
                m.p_dyn.x[m.im-1][j][k] = m.p_dyn.x[m.im-4][j][k]
                    - 3.0 * m.p_dyn.x[m.im-3][j][k]
                    + 3.0 * m.p_dyn.x[m.im-2][j][k];
            }
        }

        // Theta: von Neumann at both poles
        #pragma omp parallel for collapse(2)
        for (int k = 0; k < m.km; k++) {
            for (int i = 0; i < m.im; i++) {
                m.p_dyn.x[i][0][k]      = m.c43 * m.p_dyn.x[i][1][k]
                                          - m.c13 * m.p_dyn.x[i][2][k];
                m.p_dyn.x[i][m.jm-1][k] = m.c43 * m.p_dyn.x[i][m.jm-2][k]
                                          - m.c13 * m.p_dyn.x[i][m.jm-3][k];
            }
        }

        // Phi: von Neumann + periodic average at Greenwich meridian
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                m.p_dyn.x[i][j][0]      = m.c43 * m.p_dyn.x[i][j][1]
                                          - m.c13 * m.p_dyn.x[i][j][2];
                m.p_dyn.x[i][j][m.km-1] = m.c43 * m.p_dyn.x[i][j][m.km-2]
                                          - m.c13 * m.p_dyn.x[i][j][m.km-3];
                m.p_dyn.x[i][j][0] = m.p_dyn.x[i][j][m.km-1]
                    = 0.5 * (m.p_dyn.x[i][j][0] + m.p_dyn.x[i][j][m.km-1]);
            }
        }

        }  // sweep

        #undef LAND

        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        if (verbose) {
            printf(" time measured: %.3f seconds for PressureSolverHyd\n",
                   elapsed.count() * 1e-9);
            cout << "      OGCM: PressureSolverHyd ended" << endl;
        }
    }

    // ====================================================================
    // Project the initial velocity field divergence-free (mirror of the
    // atmosphere's PressureSolverAtm::project_initial_velocity). The ocean IC
    // (EkmanSpiral + thermohaline currents) is NOT divergence-free and was never
    // projected, so the single-sweep in-loop solver could never catch up and the
    // residual divergence drove the radial-velocity runaway. Feed the velocity
    // itself as the Poisson source (aux <- v), iterate to approximate the
    // projection pressure, apply v <- v - grad p in the same metric form as the
    // RHS pressure-gradient term, then clear p_dyn/aux so the time loop starts
    // fresh. See project_hydro_polar_blowup.
    // ====================================================================
    void project_initial_velocity(int n_sweeps = 200)
    {
        using namespace std;
        cout << endl << endl << "      OGCM: project_initial_velocity ("
             << n_sweeps << " Jacobi sweeps)" << endl;
        auto t0 = std::chrono::high_resolution_clock::now();

        // Step 1 — seed Poisson source from current velocity, zero p_dyn.
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < m.im; i++)
            for (int j = 0; j < m.jm; j++)
                for (int k = 0; k < m.km; k++) {
                    m.aux_u.x[i][j][k] = m.u.x[i][j][k];
                    m.aux_v.x[i][j][k] = m.v.x[i][j][k];
                    m.aux_w.x[i][j][k] = m.w.x[i][j][k];
                    m.p_dyn.x[i][j][k] = 0.0;
                }

        // Step 2 — Jacobi iteration until p_dyn approximates the projection pressure.
        run(n_sweeps, false);

        // Step 3 — gradient correction v <- v - grad p_dyn in the interior.
        // Metric factors match rhs_u/v/w (exp_rm = 1/(rm+1), inv_rm, inv_rmsinthe).
        const double inv_2dthe = 1.0 / (2.0 * m.dthe);
        const double inv_2dphi = 1.0 / (2.0 * m.dphi);
        // radial dp/dr uses the per-i non-uniform central coeffs (rc1*), not inv_2dr

        std::vector<double> sinthe_tab(m.jm);
        for (int j = 0; j < m.jm; j++) {
            sinthe_tab[j] = sin(m.the.z[j]);
            if (sinthe_tab[j] < 0.4) sinthe_tab[j] = 0.4;   // metric floor — match run()
        }

        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 1; i < m.im-1; i++) {
            for (int j = 1; j < m.jm-1; j++) {
                const double rm           = m.rad.z[i];
                const double exp_rm       = 1.0 / (rm + 1.0);
                const double inv_rm       = 1.0 / rm;
                const double inv_rmsinthe = 1.0 / (rm * sinthe_tab[j]);

                for (int k = 1; k < m.km-1; k++) {
                    if (is_land(m.h, i, j, k)) continue;
                    const double dpdr   = m.rc1m[i]*m.p_dyn.x[i-1][j][k]
                                        + m.rc10[i]*m.p_dyn.x[i][j][k]
                                        + m.rc1p[i]*m.p_dyn.x[i+1][j][k];
                    const double dpdthe = (m.p_dyn.x[i][j+1][k] - m.p_dyn.x[i][j-1][k]) * inv_2dthe;
                    const double dpdphi = (m.p_dyn.x[i][j][k+1] - m.p_dyn.x[i][j][k-1]) * inv_2dphi;

                    m.u.x[i][j][k] -= dpdr   * exp_rm;
                    m.v.x[i][j][k] -= dpdthe * inv_rm;
                    m.w.x[i][j][k] -= dpdphi * inv_rmsinthe;
                }
            }
        }

        // Step 4 — clear p_dyn and aux so the time loop starts fresh.
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < m.im; i++)
            for (int j = 0; j < m.jm; j++)
                for (int k = 0; k < m.km; k++) {
                    m.p_dyn.x[i][j][k] = 0.0;
                    m.aux_u.x[i][j][k] = 0.0;
                    m.aux_v.x[i][j][k] = 0.0;
                    m.aux_w.x[i][j][k] = 0.0;
                }

        auto t1 = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0);
        printf("      OGCM: project_initial_velocity ended (%.3fs)\n", elapsed.count() * 1e-9);
    }

private:
    cHydrosphereModel& m;
};
