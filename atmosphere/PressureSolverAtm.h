#pragma once

#include "cAtmosphereModel.h"
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

class PressureSolverAtm {
public:
    explicit PressureSolverAtm(cAtmosphereModel& model)
        : m(model)
    {}

    void run(bool verbose = true)
    {
        using namespace std;
        if (verbose) cout << endl << endl << endl << "      ATOM: PressureSolverAtm" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        // precompute sin(the) table — only depends on j, avoids redundant sin() calls
        std::vector<double> sinthe_table(m.jm);
        for (int j = 0; j < m.jm; j++) {
            sinthe_table[j] = sin(m.the.z[j]);
            if (sinthe_table[j] < 0.4) sinthe_table[j] = 0.4;
        }

        // Precompute land mask — eliminates repeated function call overhead
        // Allocate flat mask: 1 = land, 0 = air
        std::vector<int8_t> land(m.im * m.jm * m.km);
        #pragma omp parallel for collapse(3)
        for (int i = 0; i < m.im; i++)
            for (int j = 0; j < m.jm; j++)
                for (int k = 0; k < m.km; k++)
                    land[i*m.jm*m.km + j*m.km + k] = is_land(m.h, i, j, k) ? 1 : 0;

        #define LAND(i,j,k) land[(i)*m.jm*m.km + (j)*m.km + (k)]

        // Fuse the three boundary loops into one pass
        #pragma omp parallel for collapse(2)
        for (int j = 1; j < m.jm-1; j++) {
            for (int k = 1; k < m.km-1; k++) {
                m.aux_u.x[0][j][k]      = m.c43 * m.aux_u.x[1][j][k]      - m.c13 * m.aux_u.x[2][j][k];
                m.aux_u.x[m.im-1][j][k] = m.c43 * m.aux_u.x[m.im-2][j][k] - m.c13 * m.aux_u.x[m.im-3][j][k];
                m.aux_v.x[0][j][k]      = m.c43 * m.aux_v.x[1][j][k]      - m.c13 * m.aux_v.x[2][j][k];
                m.aux_v.x[m.im-1][j][k] = m.c43 * m.aux_v.x[m.im-2][j][k] - m.c13 * m.aux_v.x[m.im-3][j][k];
                m.aux_w.x[0][j][k]      = m.c43 * m.aux_w.x[1][j][k]      - m.c13 * m.aux_w.x[2][j][k];
                m.aux_w.x[m.im-1][j][k] = m.c43 * m.aux_w.x[m.im-2][j][k] - m.c13 * m.aux_w.x[m.im-3][j][k];
            }
        }

        // Grid-spacing reciprocals — constant for the entire grid
        const double inv_2dr   = 1.0 / (2.0 * m.dr);
        const double inv_2dthe = 1.0 / (2.0 * m.dthe);
        const double inv_2dphi = 1.0 / (2.0 * m.dphi);
        const double inv_dr2   = 1.0 / (m.dr   * m.dr);
        const double inv_dthe2 = 1.0 / (m.dthe * m.dthe);
        const double inv_dphi2 = 1.0 / (m.dphi * m.dphi);
        const double inv_dthe  = 1.0 / m.dthe;
        const double inv_dphi  = 1.0 / m.dphi;

        // Main compute loop — land mask lookups + hoisted j-invariants + k sliding window
        #pragma omp parallel for collapse(2) schedule(dynamic, 4)
        for (int i = 1; i < m.im-1; i++) {
            for (int j = 1; j < m.jm-1; j++) {

                // Build geometry struct ONCE per (i,j)
                cAtmosphereModel::CellGeometry geo;

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

                const double denom = 2.0 * geo.exp_2_rm    * inv_dr2
                                   + 2.0 * geo.inv_rm       * inv_dthe2
                                   + 2.0 * geo.inv_rmsinthe * inv_dphi2;
                const double inv_denom = 1.0 / denom;
                const double num1 = geo.exp_2_rm    * inv_dr2;
                const double num2 = geo.inv_rm      * inv_dthe2;
                const double num3 = geo.inv_rmsinthe * inv_dphi2;

                const bool i_in_range = (i < m.im-2);
                const bool j_inner    = (j > 2) && (j < m.jm-2);

                // sliding window for k-direction land/air status
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

                    // r direction
                    if (i_in_range && lnd_ijk && !LAND(i+1,j,k)) {
                        du_dr = (-3.0 * m.aux_u.x[i][j][k] + 4.0 * m.aux_u.x[i+1][j][k]
                                 - m.aux_u.x[i+2][j][k]) * inv_2dr;
                        r_flag = true;
                    }

                    // theta direction
                    if (j_inner) {
                        const int8_t air_jp1 = !LAND(i,j+1,k);
                        const int8_t air_jm1 = !LAND(i,j-1,k);

                        if (lnd_ijk && air_jp1 && !LAND(i,j+2,k)) {
                            dv_dthe = (-3.0 * m.aux_v.x[i][j][k] + 4.0 * m.aux_v.x[i][j+1][k]
                                       - m.aux_v.x[i][j+2][k]) * inv_2dthe;
                            the_flag = true;
                        }
                        if (lnd_ijk && air_jm1 && !LAND(i,j-2,k)) {
                            dv_dthe = -(-3.0 * m.aux_v.x[i][j][k] + 4.0 * m.aux_v.x[i][j-1][k]
                                        - m.aux_v.x[i][j-2][k]) * inv_2dthe;
                            the_flag = true;
                        }
                        if ((lnd_ijk && air_jp1 && LAND(i,j+2,k))
                            || (j == m.jm-2 && !lnd_ijk && LAND(i,j+1,k))) {
                            dv_dthe = (m.aux_v.x[i][j+1][k] - m.aux_v.x[i][j][k]) * inv_dthe;
                            the_flag = true;
                        }
                        if ((lnd_ijk && air_jm1 && LAND(i,j-2,k))
                            || (j == 1 && lnd_ijk && air_jm1)) {
                            dv_dthe = (m.aux_v.x[i][j-1][k] - m.aux_v.x[i][j][k]) * inv_dthe;
                            the_flag = true;
                        }
                    }

                    // phi direction
                    const bool k_inner = (k > 2) && (k < m.km-2);
                    if (k_inner) {
                        const int8_t air_kp1 = !lnd_kp1;
                        const int8_t air_km1 = !lnd_km1;

                        if (lnd_ijk && air_kp1 && !lnd_kp2) {
                            dw_dphi = (-3.0*m.aux_w.x[i][j][k] + 4.0*m.aux_w.x[i][j][k+1]
                                       - m.aux_w.x[i][j][k+2]) * inv_2dphi;
                            phi_flag = true;
                        }
                        if (lnd_ijk && air_km1 && !lnd_km2) {
                            dw_dphi = -(-3.0*m.aux_w.x[i][j][k] + 4.0*m.aux_w.x[i][j][k-1]
                                        - m.aux_w.x[i][j][k-2]) * inv_2dphi;
                            phi_flag = true;
                        }
                        if ((lnd_ijk && air_kp1 && lnd_kp2)
                            || (k == m.km-2 && !lnd_ijk && lnd_kp1)) {
                            dw_dphi = (m.aux_w.x[i][j][k+1] - m.aux_w.x[i][j][k]) * inv_dphi;
                            phi_flag = true;
                        }
                        if ((lnd_ijk && air_km1 && lnd_km2)
                            || (k == 1 && lnd_ijk && air_km1)) {
                            dw_dphi = (m.aux_w.x[i][j][k-1] - m.aux_w.x[i][j][k]) * inv_dphi;
                            phi_flag = true;
                        }
                    }

                    // central-difference fallbacks
                    if (!r_flag)
                        du_dr   = (m.aux_u.x[i+1][j][k] - m.aux_u.x[i-1][j][k]) * inv_2dr;
                    if (!the_flag)
                        dv_dthe = (m.aux_v.x[i][j+1][k] - m.aux_v.x[i][j-1][k]) * inv_2dthe;
                    if (!phi_flag)
                        dw_dphi = (m.aux_w.x[i][j][k+1] - m.aux_w.x[i][j][k-1]) * inv_2dphi;

                    // pressure update
                    m.p_dyn.x[i][j][k] =
                        ((m.p_dyn.x[i+1][j][k] + m.p_dyn.x[i-1][j][k]) * num1
                       + (m.p_dyn.x[i][j+1][k] + m.p_dyn.x[i][j-1][k]) * num2
                       + (m.p_dyn.x[i][j][k+1] + m.p_dyn.x[i][j][k-1]) * num3
                       - du_dr   * geo.exp_rm
                       - dv_dthe * geo.inv_rm
                       - dw_dphi * geo.inv_rmsinthe) * inv_denom;
                } // k
            } // j
        } // i

        #undef LAND

        // Radial boundary extrapolation.
        //
        // At i=0 the old cubic p[0] = p[3] - 3*p[2] + 3*p[1] has condition number ~7
        // and amplifies any cliff-cell pressure into the reference layer by up to 3×.
        // Over many pressure-solve calls this compounds into huge ± spikes at steep
        // topography (e.g. Himalaya). Replaced by:
        //   - i_topography[j][k] >= 1 (i=0 inside the mountain): hold p_dyn = 0, matching
        //     bcSolidGround's treatment of fully buried interior cells; acts as a Dirichlet
        //     pin for the otherwise all-Neumann pressure Poisson.
        //   - i_topography[j][k] == 0 (i=0 is a real ocean surface): von Neumann zero-gradient
        //     ∂p/∂n = 0 — physically correct at a free-slip wall.
        // Tried fully zero-gradient (no Dirichlet anywhere) — Poisson becomes singular and
        // diverges to NaN within 100 iters at the pole.
        // At i=im-1 the column is always air, so the cubic is fine and is kept.
        #pragma omp parallel for collapse(2)
        for (int k = 0; k < m.km; k++) {
            for (int j = 0; j < m.jm; j++) {
                if (m.i_topography[j][k] >= 1) {
                    m.p_dyn.x[0][j][k] = 0.0;
                } else {
                    m.p_dyn.x[0][j][k] = m.c43 * m.p_dyn.x[1][j][k]
                                       - m.c13 * m.p_dyn.x[2][j][k];
                }
                m.p_dyn.x[m.im-1][j][k] = m.p_dyn.x[m.im-4][j][k]
                    - 3.0 * m.p_dyn.x[m.im-3][j][k]
                    + 3.0 * m.p_dyn.x[m.im-2][j][k];
            }
        }

        // Theta-pole BC: plain copy (zero-gradient, no extrapolation).
        // The 2nd-order Neumann form p[0] = (4/3)·p[1] − (1/3)·p[2] amplifies polar
        // grid noise by 4/3 per call. The Poisson iteration is rerun every moist_stride,
        // and the metric 1/sin²θ already blows up near the pole — compounding the 4/3
        // factor drove p_dyn to NaN at 90°N within ~150 iters, then dpdr/dpdthe/dpdphi
        // propagated NaN into u,v,w in the next RHS call. Plain copy = amplification 1
        // and matches the axisymmetric-pole assumption used for the same fields in bcTheta.
        #pragma omp parallel for collapse(2)
        for (int k = 0; k < m.km; k++) {
            for (int i = 0; i < m.im; i++) {
                m.p_dyn.x[i][0][k]      = m.p_dyn.x[i][1][k];
                m.p_dyn.x[i][m.jm-1][k] = m.p_dyn.x[i][m.jm-2][k];
            }
        }

        // Phi boundary average
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                m.p_dyn.x[i][j][0]      = m.c43 * m.p_dyn.x[i][j][1]      - m.c13 * m.p_dyn.x[i][j][2];
                m.p_dyn.x[i][j][m.km-1] = m.c43 * m.p_dyn.x[i][j][m.km-2] - m.c13 * m.p_dyn.x[i][j][m.km-3];
                m.p_dyn.x[i][j][0] = m.p_dyn.x[i][j][m.km-1]
                    = (m.p_dyn.x[i][j][0] + m.p_dyn.x[i][j][m.km-1]) / 2.0;
            }
        }

        auto end = std::chrono::high_resolution_clock::now();
        if (verbose) {
            auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
            printf(" time measured: %.3f seconds for PressureSolverAtm\n", elapsed.count() * 1e-9);
            cout << "      ATOM: PressureSolverAtm ended" << endl;
        }
    }

    // One-shot Helmholtz projection of the prescribed initial velocity onto its
    // divergence-free subspace.  Without this, the incremental projection that runs
    // every moist_stride iters spends the first ~100 iters silently destroying the
    // dilatational component of the Hadley/Ferrel profile, taking the prescribed
    // global circulation with it (max u decays from 30 m/s to <0.5 m/s by iter 150).
    // Performing the projection once at t=0 lets the simulation start from a clean
    // divergence-free state; the surviving solenoidal part is preserved and only the
    // unphysical dilatational artefact of the analytical profile is removed.
    //
    // The standard solver loop assembles ∇²p = ∇·v* where v* is the post-RHS
    // intermediate velocity; here we feed v itself as the source by copying v into
    // aux_u/aux_v/aux_w. After n_sweeps Jacobi passes, p_dyn approximates the
    // projection pressure; we then apply v ← v − ∇p in the same metric form used
    // by the time-stepping RHS, and reset p_dyn to 0 so the next RK4 call does not
    // double-correct via its own −∂p/∂r term.
    void project_initial_velocity(int n_sweeps = 200)
    {
        using namespace std;
        cout << endl << endl << "      ATOM: project_initial_velocity ("
             << n_sweeps << " Jacobi sweeps)" << endl;
        auto t0 = std::chrono::high_resolution_clock::now();

        // Step 1 — seed Poisson source from current velocity, zero p_dyn.
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                for (int k = 0; k < m.km; k++) {
                    m.aux_u.x[i][j][k] = m.u.x[i][j][k];
                    m.aux_v.x[i][j][k] = m.v.x[i][j][k];
                    m.aux_w.x[i][j][k] = m.w.x[i][j][k];
                    m.p_dyn.x[i][j][k] = 0.0;
                }
            }
        }

        // Step 2 — quiet Jacobi iteration until p_dyn approximates the projection
        // pressure.  Each call to run() also re-applies the i, theta, and phi BCs on
        // p_dyn, so polar/topographic anchors stay consistent with the time loop.
        for (int s = 0; s < n_sweeps; s++) {
            run(false);
        }

        // Step 3 — gradient correction v ← v − ∇p_dyn in the interior.
        // Metric factors match the rhs_u/v/w pressure-gradient term so the magnitudes
        // are consistent with the rest of the code.  Boundary cells (i, theta, phi
        // outer faces) are left untouched; bcRadius / bcTheta / bcPhi will re-impose
        // their patterns at the next call.
        const double inv_2dr   = 1.0 / (2.0 * m.dr);
        const double inv_2dthe = 1.0 / (2.0 * m.dthe);
        const double inv_2dphi = 1.0 / (2.0 * m.dphi);

        std::vector<double> sinthe_tab(m.jm);
        for (int j = 0; j < m.jm; j++) {
            sinthe_tab[j] = sin(m.the.z[j]);
            if (sinthe_tab[j] < 0.4) sinthe_tab[j] = 0.4;
        }

        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 1; i < m.im-1; i++) {
            for (int j = 1; j < m.jm-1; j++) {
                const double rm           = m.rad.z[i];
                const double exp_rm       = 1.0 / (rm + 1.0);
                const double inv_rm       = 1.0 / rm;
                const double inv_rmsinthe = 1.0 / (rm * sinthe_tab[j]);

                for (int k = 1; k < m.km-1; k++) {
                    const double dpdr   = (m.p_dyn.x[i+1][j][k] - m.p_dyn.x[i-1][j][k]) * inv_2dr;
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
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                for (int k = 0; k < m.km; k++) {
                    m.p_dyn.x[i][j][k] = 0.0;
                    m.aux_u.x[i][j][k] = 0.0;
                    m.aux_v.x[i][j][k] = 0.0;
                    m.aux_w.x[i][j][k] = 0.0;
                }
            }
        }

        auto t1 = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0);
        printf("      ATOM: project_initial_velocity ended (%.3fs)\n", elapsed.count() * 1e-9);
    }

private:
    cAtmosphereModel& m;
};
