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

    // ====================================================================
    // Face-consistent fractional-step projection (the real Rhie-Chow fix).
    //
    // The first cut (advect with faces, but project with run()'s rc1 cell
    // divergence + central cell correction) BLEW UP: the cell-level checkerboard
    // feedback survived (source read the cell odd-even mode -> phi -> correction
    // re-injected it). Cure = make the SOURCE, LAPLACIAN and GRADIENT a single
    // consistent triple, all face-based, so div_face(corrected face flux) = 0 by
    // construction and the source is BLIND to the cell checkerboard:
    //
    //   source  D(u*) = exp_rm*(2/hs)*(Uf~[i]-Uf~[i-1]) + inv_rm/dthe*(Vf~-Vf~)
    //                   + inv_rmsinthe/dphi*(Wf~-Wf~),   Uf~=0.5(u[i]+u[i+1])
    //   Laplacian L(p) = exp_2_rm*rc2 (radial, factors as the same face diff)
    //                   + inv_rm2/dthe2 * d2_theta + inv_rm2sinthe2/dphi2 * d2_phi
    //   face flux Uf   = Uf~ - exp_rm_face*(p[i+1]-p[i])/hp   (etc.)
    //
    // Crucially the theta/phi Laplacian uses inv_rm2 / inv_rm2sinthe2 (= the
    // metric of div o grad), NOT run()'s single inv_rm / inv_rmsinthe — that
    // single-power form is the latent inconsistency that made the collocated
    // projection non-idempotent. All operators are LAND-MASKED as homogeneous
    // Neumann (drop land neighbours; no flux through land faces), which also
    // fixes the coastal blow-up. The div-free faces uf/vf/wf are the advecting
    // velocity (RHS_Hyd); the cell velocity gets the pointwise gradient too.
    // Uses aux_u as the per-cell source scratch. See
    // project_hydro_continuity_checkerboard.
    // ====================================================================
    void project_velocity(int n_sweeps = 60)
    {
        using namespace std;
        cout << endl << "      OGCM: project_velocity (face-consistent, "
             << n_sweeps << " sweeps)" << endl;

        const double inv_dthe  = 1.0 / m.dthe;
        const double inv_dphi  = 1.0 / m.dphi;
        const double inv_dthe2 = inv_dthe * inv_dthe;
        const double inv_dphi2 = inv_dphi * inv_dphi;

        std::vector<double> sinthe_tab(m.jm);
        for (int j = 0; j < m.jm; j++) {
            sinthe_tab[j] = sin(m.the.z[j]);
            if (sinthe_tab[j] < 0.4) sinthe_tab[j] = 0.4;   // metric floor — match run()
        }

        auto water = [&](int i, int j, int k){ return is_water(m.h, i, j, k); };

        // ---- Step 1: source = face-average divergence of the current velocity,
        // with no-flux masking through land faces. Stored in aux_u; p_dyn zeroed.
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < m.im; i++)
            for (int j = 0; j < m.jm; j++)
                for (int k = 0; k < m.km; k++)
                    m.p_dyn.x[i][j][k] = 0.0;

        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 1; i < m.im-1; i++) {
            for (int j = 1; j < m.jm-1; j++) {
                const double rm           = m.rad.z[i];
                const double exp_rm       = 1.0 / (rm + 1.0);
                const double inv_rm       = 1.0 / rm;
                const double inv_rmsinthe = 1.0 / (rm * sinthe_tab[j]);
                const double hs_r         = m.rad.z[i+1] - m.rad.z[i-1];
                const double two_over_hs  = 2.0 / hs_r;
                for (int k = 1; k < m.km-1; k++) {
                    if (!water(i, j, k)) { m.aux_u.x[i][j][k] = 0.0; continue; }
                    // intermediate face fluxes (0 through a land face)
                    const double ufRp = water(i+1,j,k) ? 0.5*(m.u.x[i][j][k]+m.u.x[i+1][j][k]) : 0.0;
                    const double ufRm = water(i-1,j,k) ? 0.5*(m.u.x[i-1][j][k]+m.u.x[i][j][k]) : 0.0;
                    const double ufTp = water(i,j+1,k) ? 0.5*(m.v.x[i][j][k]+m.v.x[i][j+1][k]) : 0.0;
                    const double ufTm = water(i,j-1,k) ? 0.5*(m.v.x[i][j-1][k]+m.v.x[i][j][k]) : 0.0;
                    const double ufPp = water(i,j,k+1) ? 0.5*(m.w.x[i][j][k]+m.w.x[i][j][k+1]) : 0.0;
                    const double ufPm = water(i,j,k-1) ? 0.5*(m.w.x[i][j][k-1]+m.w.x[i][j][k]) : 0.0;
                    m.aux_u.x[i][j][k] = exp_rm       * two_over_hs * (ufRp - ufRm)
                                       + inv_rm       * inv_dthe    * (ufTp - ufTm)
                                       + inv_rmsinthe * inv_dphi    * (ufPp - ufPm);
                }
            }
        }

        // ---- Step 2: masked-Neumann Jacobi solve  L(p) = source.
        for (int sweep = 0; sweep < n_sweeps; sweep++) {
            #pragma omp parallel for collapse(2) schedule(static)
            for (int i = 1; i < m.im-1; i++) {
                for (int j = 1; j < m.jm-1; j++) {
                    const double rm      = m.rad.z[i];
                    const double exp_2_rm = 1.0 / ((rm + 1.0) * (rm + 1.0));
                    const double inv_rm2  = 1.0 / (rm * rm);
                    const double inv_rm2sinthe2 = inv_rm2 / (sinthe_tab[j] * sinthe_tab[j]);
                    const double cT = inv_rm2 * inv_dthe2;
                    const double cP = inv_rm2sinthe2 * inv_dphi2;
                    for (int k = 1; k < m.km-1; k++) {
                        if (!water(i, j, k)) continue;
                        const double aRm = water(i-1,j,k) ? exp_2_rm * m.rc2m[i] : 0.0;
                        const double aRp = water(i+1,j,k) ? exp_2_rm * m.rc2p[i] : 0.0;
                        const double aTm = water(i,j-1,k) ? cT : 0.0;
                        const double aTp = water(i,j+1,k) ? cT : 0.0;
                        const double aPm = water(i,j,k-1) ? cP : 0.0;
                        const double aPp = water(i,j,k+1) ? cP : 0.0;
                        const double diag = aRm + aRp + aTm + aTp + aPm + aPp;
                        if (diag <= 0.0) continue;
                        m.p_dyn.x[i][j][k] =
                            ( aRm*m.p_dyn.x[i-1][j][k] + aRp*m.p_dyn.x[i+1][j][k]
                            + aTm*m.p_dyn.x[i][j-1][k] + aTp*m.p_dyn.x[i][j+1][k]
                            + aPm*m.p_dyn.x[i][j][k-1] + aPp*m.p_dyn.x[i][j][k+1]
                            - m.aux_u.x[i][j][k] ) / diag;
                    }
                }
            }
            // zero-gradient (Neumann) pressure BCs on the domain faces
            #pragma omp parallel for collapse(2)
            for (int j = 0; j < m.jm; j++)
                for (int k = 0; k < m.km; k++) {
                    m.p_dyn.x[0][j][k]      = m.p_dyn.x[1][j][k];
                    m.p_dyn.x[m.im-1][j][k] = m.p_dyn.x[m.im-2][j][k];
                }
            #pragma omp parallel for collapse(2)
            for (int i = 0; i < m.im; i++)
                for (int k = 0; k < m.km; k++) {
                    m.p_dyn.x[i][0][k]      = m.p_dyn.x[i][1][k];
                    m.p_dyn.x[i][m.jm-1][k] = m.p_dyn.x[i][m.jm-2][k];
                }
            #pragma omp parallel for collapse(2)
            for (int i = 0; i < m.im; i++)
                for (int j = 0; j < m.jm; j++) {
                    m.p_dyn.x[i][j][0]      = m.p_dyn.x[i][j][1];
                    m.p_dyn.x[i][j][m.km-1] = m.p_dyn.x[i][j][m.km-2];
                }
        }

        // ---- Step 3: divergence-free face fluxes  Uf = Uf~ - grad_face(p).
        // Same metric on the face gradient as in the source (exp_rm/inv_rm/
        // inv_rmsinthe), so div_face(Uf) = source - L(p) -> 0. Zero on land faces.
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                const double rm           = m.rad.z[i];
                const double inv_rm       = 1.0 / rm;
                const double inv_rmsinthe = 1.0 / (rm * sinthe_tab[j]);
                for (int k = 0; k < m.km; k++) {
                    if (i < m.im-1 && water(i,j,k) && water(i+1,j,k)) {
                        const double rmf = 0.5*(m.rad.z[i]+m.rad.z[i+1]);
                        const double exp_rm_face = 1.0/(rmf+1.0);
                        const double hp = m.rad.z[i+1]-m.rad.z[i];
                        m.uf.x[i][j][k] = 0.5*(m.u.x[i][j][k]+m.u.x[i+1][j][k])
                            - exp_rm_face*(m.p_dyn.x[i+1][j][k]-m.p_dyn.x[i][j][k])/hp;
                    } else m.uf.x[i][j][k] = 0.0;

                    if (j < m.jm-1 && water(i,j,k) && water(i,j+1,k))
                        m.vf.x[i][j][k] = 0.5*(m.v.x[i][j][k]+m.v.x[i][j+1][k])
                            - inv_rm*(m.p_dyn.x[i][j+1][k]-m.p_dyn.x[i][j][k])*inv_dthe;
                    else m.vf.x[i][j][k] = 0.0;

                    if (k < m.km-1 && water(i,j,k) && water(i,j,k+1))
                        m.wf.x[i][j][k] = 0.5*(m.w.x[i][j][k]+m.w.x[i][j][k+1])
                            - inv_rmsinthe*(m.p_dyn.x[i][j][k+1]-m.p_dyn.x[i][j][k])*inv_dphi;
                    else m.wf.x[i][j][k] = 0.0;
                }
            }
        }

        // ---- Step 4: reconstruct the cell velocity as the average of its two
        // bounding divergence-free faces (uf/vf/wf from Step 3), instead of a
        // SEPARATE central-gradient pressure correction.
        //
        // The old central-gradient form (u -= dpdr*exp_rm) is a DIFFERENT discrete
        // operator from the face flux, so the collocated cell velocity did NOT
        // inherit the faces' divergence-free property — it accumulated spurious
        // divergence that Step 1 re-read as the source next call and the projection
        // re-amplified: a self-feeding u->source->u loop that ran the deepest, most-
        // convergent cell away (first the metric-amplified North Pole ~iter600, then
        // — once the pole was filtered — a deep Southern-Ocean cell at 52 S ~iter1150;
        // the instability is latitude-independent, see project_hydro_polar_blowup).
        //
        // Face-averaging is the standard staggered->collocated reconstruction: it
        // makes the prognostic cell velocity consistent with the div-free faces, so
        // its divergence (hence the next source) stays small and the loop is broken.
        // uf[i-1]/vf[j-1]/wf[k-1] are the low-side faces of cell (i,j,k), uf[i]/
        // vf[j]/wf[k] the high-side faces; faces are 0 on land (no-normal-flow), so a
        // cell against a boundary correctly gets a halved (reduced) velocity. The
        // average also carries an implicit 1-2-1 smoothing of the pre-projection
        // velocity (mild, stabilising — inherent to Rhie-Chow collocation).
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 1; i < m.im-1; i++) {
            for (int j = 1; j < m.jm-1; j++) {
                for (int k = 1; k < m.km-1; k++) {
                    if (!water(i, j, k)) continue;
                    m.u.x[i][j][k] = 0.5 * (m.uf.x[i-1][j][k] + m.uf.x[i][j][k]);
                    m.v.x[i][j][k] = 0.5 * (m.vf.x[i][j-1][k] + m.vf.x[i][j][k]);
                    m.w.x[i][j][k] = 0.5 * (m.wf.x[i][j][k-1] + m.wf.x[i][j][k]);
                }
            }
        }
    }

    // ------------------------------------------------------------------------
    // project_barotropic — rigid-lid barotropic (Stommel) streamfunction solve.
    // Supplies the wind-driven barotropic transport the Chorin/face projection
    // cannot build (project_hydro_eastern_boundary_current: interior is Coriolis-
    // vs-diffusion, PGF ~0). Solves, on the ocean surface sheet (psi=0 on land &
    // |lat|>78, periodic in longitude):
    //     eps*Lap(psi) + (2*omega/a^2)*dpsi/dphi = curl(tau_wind)/(rho*H)
    // then stores the barotropic (depth-mean) VELOCITY u_bt = curl^-1... i.e.
    // (v_bt,w_bt) = (-dpsi/dx, dpsi/dy)/H, non-dim by u_0, in m.v_bt/w_bt.
    // apply_barotropic_mode_split() imposes it each iter as the column depth-mean.
    // NB: injecting it as a geostrophic PGF *force* (the earlier E.3 port) is
    // self-defeating — f x u_bt is curl-free, so the divergence projection deletes
    // it (verified: 1.5% flow change, ACC stays deep-westward). Setting the
    // velocity mode directly cannot be projected out.
    // Called once after the surface wind is read (wind is fixed during a hyd run).
    // Validated offline: scratchpad/{stommel_test,baro_e2}.py  (E.1/E.2).
    // ------------------------------------------------------------------------
    void project_barotropic(double eps = 3.0e-6, double gain = 1.0, int n_sweeps = 12000)
    {
        using namespace std;
        cout << endl << "      OGCM: project_barotropic (Stommel eps=" << eps
             << ", gain=" << gain << ", " << n_sweeps << " SOR sweeps)" << endl;
        const double a     = 6371000.0;                  // earth radius [m]
        const double omega = m.omega;                    // [rad/s]
        const double rho   = m.r_0_water;                // [kg/m^3]
        const double C_D   = 2.6e-3, r_air = m.r_air;
        const double dthe  = m.dthe, dphi = m.dphi;
        const int jm = m.jm, km = m.km, Kp = km - 1;     // longitude period = km-1
        auto kper = [&](int k){ k %= Kp; if (k < 0) k += Kp; return k; };
        auto sfl  = [&](double s){ return s < 1.0e-3 ? 1.0e-3 : s; };
        auto latd = [&](int j){ return 90.0 - (double)j; };
        auto ocean= [&](int j,int k){ return is_water(m.h, m.im-2, j, k)
                                        && j>0 && j<jm-1 && std::fabs(latd(j))<=78.0; };
        std::vector<double> sinth(jm), costh(jm);
        for (int j=0;j<jm;j++){ sinth[j]=std::sin(m.the.z[j]); costh[j]=std::cos(m.the.z[j]); }
        std::vector<std::vector<double>> H(jm,std::vector<double>(km,0.0)),
                                         W(jm,std::vector<double>(km,0.0)),
                                         Psi(jm,std::vector<double>(km,0.0));
        // Barotropic transport is a FULL-DEPTH quantity, so convert it to a
        // depth-mean velocity u_bt = transport/H with the TRUE seafloor depth
        // (m.Bathymetry, metres — raw topo, populated for both modes), NOT the
        // model column depth c*dz. In shallow mode c*dz caps at L_hyd = 200 m,
        // which over-speeds u_bt ~20x (transport/200 vs transport/~4000) and
        // drove the CFL runaway that forced the deep-only gate. Using the true
        // depth makes u_bt mode-independent (deep mode is unchanged: c*dz there
        // already ≈ Bathymetry over the resolved column). Floor 200 m keeps
        // shelves/land from blowing up, matching the old floor.
        for (int j=0;j<jm;j++) for(int k=0;k<km;k++)
            H[j][k] = std::max(m.Bathymetry.y[j][k], 200.0);
        auto tauP=[&](int j,int k){ double v=m.v_wind.y[j][k], w=m.w_wind.y[j][k];
                                    return r_air*C_D*std::sqrt(v*v+w*w)*w; };   // east
        auto tauT=[&](int j,int k){ double v=m.v_wind.y[j][k], w=m.w_wind.y[j][k];
                                    return r_air*C_D*std::sqrt(v*v+w*w)*v; };   // colatitude (south+)
        for (int j=1;j<jm-1;j++){ double st=sfl(sinth[j]);
            for(int k=0;k<km;k++){ int kp=kper(k+1), kmn=kper(k-1);
                double dsinTauP=(sinth[j+1]*tauP(j+1,k)-sinth[j-1]*tauP(j-1,k))/(2*dthe);
                double dTauT=(tauT(j,kp)-tauT(j,kmn))/(2*dphi);
                // transport-streamfunction source curl(tau)/rho (H enters only in u_bt=transport/H)
                W[j][k]=((dsinTauP-dTauT)/(a*st))/rho;
            }
        }
        const double relax=1.5;
        for (int sweep=0; sweep<n_sweeps; sweep++){
            for (int j=1;j<jm-1;j++){
                double st=sfl(sinth[j]);
                double sTp=sfl(std::sin(m.the.z[j]+dthe/2)), sTm=sfl(std::sin(m.the.z[j]-dthe/2));
                double cTp=eps*sTp/(st*a*a*dthe*dthe), cTm=eps*sTm/(st*a*a*dthe*dthe);
                double cP =eps/(st*st*a*a*dphi*dphi);
                double bP =2.0*omega/(a*a*2.0*dphi);
                double diag=cTp+cTm+2.0*cP;
                if (diag<=0.0) continue;
                for (int k=0;k<km;k++){
                    if(!ocean(j,k)){ Psi[j][k]=0.0; continue; }
                    int kp=kper(k+1), kmn=kper(k-1);
                    double nb = cTp*Psi[j+1][k]+cTm*Psi[j-1][k]
                              + (cP+bP)*Psi[j][kp] + (cP-bP)*Psi[j][kmn];
                    Psi[j][k] += relax*((nb - W[j][k])/diag - Psi[j][k]);
                }
            }
        }
        for (int j=0;j<jm;j++) for(int k=0;k<km;k++){
            m.psi_bt.y[j][k]=Psi[j][k]; m.v_bt.y[j][k]=0.0; m.w_bt.y[j][k]=0.0; }
        double mp=0.0, mg=0.0;
        for (int j=1;j<jm-1;j++){ double st=sfl(sinth[j]);
            for(int k=0;k<km;k++){
                if(!ocean(j,k)) continue;
                int kp=kper(k+1), kmn=kper(k-1);
                double U_east =((Psi[j+1][k]-Psi[j-1][k])/(2*dthe))/a;            // -dPsi/dy [m2/s]
                double V_north=((Psi[j][kp]-Psi[j][kmn])/(2*dphi))/(a*st);        // dPsi/dx  [m2/s]
                double wbt =  U_east / H[j][k];           // barotropic zonal  (E+) [m/s]
                double vbt = -V_north / H[j][k];          // barotropic merid  (S+) [m/s]
                // Cap the barotropic speed. The discrete wind-stress curl over
                // coastal/shelf topography produces spurious transport gradients
                // that, divided by the shallow floor depth, give unphysical u_bt
                // (~80 cm/s) — the CFL runaway (~iter 65, SH shelves) that forced
                // the deep-only gate. A real barotropic depth-mean current stays
                // O(10-40 cm/s); deep mode's validated max was 37.6 cm/s, below
                // this cap, so deep behaviour is unchanged and only the shelf
                // spikes are clipped — letting the split run in shallow mode too.
                const double u_bt_max = 0.40;             // m/s
                double sp = std::sqrt(wbt*wbt + vbt*vbt);
                if (sp > u_bt_max) { double fac = u_bt_max/sp; wbt *= fac; vbt *= fac; }
                // Equatorial taper: the Sverdrup/Stommel barotropic balance is
                // ill-posed as f->0 (geostrophy fails, the source ~1/f), and the
                // wind-driven gyres + ACC we want are all extratropical. Fade u_bt
                // to zero within |lat|<10 deg (full by 20 deg) — this also removes
                // the tropical island-strait forcing (Indonesia/Timor) that drove
                // the surface CFL runaway ~iter 230.
                double al = std::fabs(latd(j));
                double eqf = al <= 10.0 ? 0.0 : (al >= 20.0 ? 1.0 : (al-10.0)/10.0);
                m.w_bt.y[j][k] = eqf*gain*wbt / m.u_0;    // non-dim; imposed as column depth-mean
                m.v_bt.y[j][k] = eqf*gain*vbt / m.u_0;    //          (velocity mode split)
                mp=std::max(mp,std::fabs(Psi[j][k]));
            }
        }
        // Smooth u_bt to reduce grid-scale divergence. A barotropic mode is (by
        // construction) nearly non-divergent, but the discrete d(Psi) + the cap +
        // the tapers leave 2*dx noise in u_bt; since the imposed depth-mean sets
        // the surface flow and EkmanPumping = div(surface velocity), that noise
        // showed up as spurious pumping at the ACC/topography choke points
        // (Kerguelen 50 S 66 E, Drake 64 S 57 W ~ +-600 cm/d). A few land-aware
        // 1-2-1 passes over ocean cells (periodic in lon; land neighbours use the
        // centre value so nothing leaks onto land) damp the 2*dx noise while
        // preserving the large-scale gyre/ACC pattern.
        {
            const int n_smooth = 6;
            std::vector<std::vector<double>> t(jm, std::vector<double>(km, 0.0));
            auto smooth = [&](Array_2D& F){
                for(int j=0;j<jm;j++) for(int k=0;k<km;k++) t[j][k]=F.y[j][k];
                for(int j=1;j<jm-1;j++) for(int k=0;k<km;k++){
                    if(!ocean(j,k)) continue;
                    const double c = t[j][k];
                    const double jp = ocean(j+1,k)      ? t[j+1][k]      : c;
                    const double jn = ocean(j-1,k)      ? t[j-1][k]      : c;
                    const double wp = ocean(j,kper(k+1))? t[j][kper(k+1)]: c;
                    const double wn = ocean(j,kper(k-1))? t[j][kper(k-1)]: c;
                    F.y[j][k] = 0.5*c + 0.125*(jp+jn+wp+wn);
                }
            };
            for(int s=0;s<n_smooth;s++){ smooth(m.v_bt); smooth(m.w_bt); }
        }
        for (int j=1;j<jm-1;j++) for(int k=0;k<km;k++) if(ocean(j,k))
            mg=std::max(mg,std::max(std::fabs(m.w_bt.y[j][k]),std::fabs(m.v_bt.y[j][k])));
        cout << "      OGCM: project_barotropic done  max|Psi|=" << mp
             << "  max|u_bt(nd)|=" << mg << "  (=" << mg*m.u_0*100.0 << " cm/s)" << endl;
    }

    // ------------------------------------------------------------------------
    // apply_barotropic_mode_split — impose the wind-driven barotropic (depth-mean)
    // velocity u_bt on the horizontal flow, preserving the baroclinic deviation:
    //     v(i) += v_bt - <v>_column ,   w(i) += w_bt - <w>_column
    // Called every iter AFTER project_velocity. This is the piece the divergence
    // projection cannot supply — a barotropic mode is (nearly) non-divergent, so
    // the projection has nothing to correct and never builds it. Setting it as a
    // velocity (not a force) is immune to being projected out. The ACC — the most
    // barotropic current — is the pass/fail signal: it needs u_bt to fill the
    // column at 55-60 S (see project_hydro_eastern_boundary_current).
    // ------------------------------------------------------------------------
    void apply_barotropic_mode_split()
    {
        // Two stabilisers, needed once the split runs from a COLD START (the deep-
        // mode validation only ever ran from a spun-up restart). Imposing the full
        // depth-mean u_bt on the near-zero from-scratch field, or onto a coastal
        // column, shocks the fine shallow grid -> CFL runaway ~iter 65-69, with the
        // epicentre at land-adjacent cells (e.g. the Australian E coast, w->600 m/s):
        //   (1) cold-start RAMP: bring the amplitude in linearly over N_ramp iters;
        //   (2) coastal TAPER: skip land-adjacent columns — the Stommel western
        //       boundary current is a thin jet the coarse grid cannot resolve at the
        //       coast. The open-ocean interior, incl. the circumpolar ACC band (the
        //       pass/fail signal), keeps the full barotropic mode.
        const int    N_ramp = 120;
        const double ramp   = std::min(1.0, (double)m.iter_n / (double)N_ramp);
        if (ramp <= 0.0) return;                       // seed call (iter_n==0) is a no-op
        const int km = m.km, Kp = km - 1, isurf = m.im - 2;
        auto kper = [&](int k){ k %= Kp; if (k < 0) k += Kp; return k; };
        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 1; j < m.jm-1; j++) {
            for (int k = 0; k < km; k++) {
                int c = 0; double mv = 0.0, mw = 0.0;
                for (int i = 0; i < m.im; i++)
                    if (is_water(m.h, i, j, k)) { mv += m.v.x[i][j][k]; mw += m.w.x[i][j][k]; ++c; }
                if (c == 0) continue;
                // coastal taper (graded, smoothstep to zero over N_taper cells):
                // a thin/strait channel carrying full-strength u_bt shocks the coarse
                // grid (CFL runaway, Indonesian straits). A HARD cutoff instead makes
                // a surface-velocity step -> divergence spike -> spurious Ekman pumping
                // (EkmanPumping = div of the surface velocity, ThermoHyd). So fade
                // u_bt smoothly to zero: weight 0 at a land-adjacent column, ->1 at
                // >=N_taper cells from land (Chebyshev distance).
                const int N_taper = 3;
                int dmin = N_taper + 1;
                for (int dj=-N_taper; dj<=N_taper; ++dj)
                    for (int dk=-N_taper; dk<=N_taper; ++dk){
                        int jj=j+dj;
                        bool land = (jj<0||jj>=m.jm) || !is_water(m.h, isurf, jj, kper(k+dk));
                        if (land){ int d=std::max(std::abs(dj),std::abs(dk)); if(d<dmin)dmin=d; }
                    }
                double cw = std::min(1.0, std::max(0.0, (double)(dmin-1)/N_taper));
                cw = cw*cw*(3.0-2.0*cw);                 // smoothstep
                if (cw <= 0.0) continue;
                const double dv = ramp * cw * (m.v_bt.y[j][k] - mv / c);
                const double dw = ramp * cw * (m.w_bt.y[j][k] - mw / c);
                for (int i = 0; i < m.im; i++)
                    if (is_water(m.h, i, j, k)) { m.v.x[i][j][k] += dv; m.w.x[i][j][k] += dw; }
            }
        }
    }

private:
    cHydrosphereModel& m;
};
