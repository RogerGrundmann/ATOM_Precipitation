#pragma once

#include "cHydrosphereModel.h"
#include "Utils.h"

#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace AtomUtils;

// ============================================================================
// BC_Hyd — friend class of cHydrosphereModel
//
// Provides all boundary condition logic for the hydrosphere grid:
//   bcRadius()          — deep-ocean (i=0) and sea-surface (i=im-1) BCs
//   bcTheta()           — north/south pole BCs
//   bcPhi()             — Greenwich meridian BCs (periodic wrap)
//   bcSolidGround()     — zero land cells; extrapolate p_dyn and salinity at
//                         land-ocean interfaces
//   boundaryCondition() — damp velocity components adjacent to solid walls
// ============================================================================
class BC_Hyd {
public:
    explicit BC_Hyd(cHydrosphereModel& model)
        : m(model)
    {}

    // ------------------------------------------------------------------
    // bcRadius — radial boundary conditions
    //
    // Extrapolation key:
    //   cubic:      x[a] = x[a+3d] - 3·x[a+2d] + 3·x[a+d]
    //   von Neumann: x[a] = c43·x[a+d] - c13·x[a+2d]
    // ------------------------------------------------------------------
    void bcRadius()
    {
        // Pattern A: cubic at i=0 AND i=im-1
        Array* both_cubic[] = {
            &m.u, &m.v, &m.w, &m.c,
            &m.BuoyancyForce, &m.CoriolisForce,
            &m.CentrifugalForce, &m.PresGradForce
        };
        constexpr int n_both = sizeof(both_cubic) / sizeof(both_cubic[0]);

        // Pattern B: von Neumann at i=0, cubic at i=im-1
        Array* vn_bot_cubic_top[] = {
            &m.Salt_Finger, &m.Salt_Diffusion, &m.Salt_Balance
        };
        constexpr int n_vn_cubic = sizeof(vn_bot_cubic_top) / sizeof(vn_bot_cubic_top[0]);

        const int iml = m.im - 1;

        #pragma omp parallel for schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {

                // t: cubic at i=0 only (i=im-1 is the sea surface — set from data)
                m.t.x[0][j][k] = m.t.x[3][j][k]
                    - 3.0 * m.t.x[2][j][k] + 3.0 * m.t.x[1][j][k];

                // Pattern A
                for (int f = 0; f < n_both; f++) {
                    double*** xf = both_cubic[f]->x;
                    xf[0][j][k]   = xf[3][j][k]       - 3.0 * xf[2][j][k]       + 3.0 * xf[1][j][k];
                    xf[iml][j][k] = xf[iml-3][j][k]   - 3.0 * xf[iml-2][j][k]   + 3.0 * xf[iml-1][j][k];
                }

                // Pattern B
                for (int f = 0; f < n_vn_cubic; f++) {
                    double*** xf = vn_bot_cubic_top[f]->x;
                    xf[0][j][k]   = m.c43 * xf[1][j][k]     - m.c13 * xf[2][j][k];
                    xf[iml][j][k] = xf[iml-3][j][k] - 3.0 * xf[iml-2][j][k] + 3.0 * xf[iml-1][j][k];
                }
            }
        }
    }
/*
* 
*/
    // ------------------------------------------------------------------
    // bcTheta — polar boundary conditions (j = 0, j = jm-1)
    // ------------------------------------------------------------------
    void bcTheta()
    {
        // Pole BC: zero-gradient (copy from first interior cell). At the spherical
        // singularity sin θ → 0, any extrapolation amplifies grid noise:
        //   - (4/3, -1/3) form has amplification factor 4/3
        //   - cubic p[0] = p[3] − 3 p[2] + 3 p[1] has condition number ~7
        // With strong bulk flow these are washed out, but with a weak / spinning-up
        // velocity field the amplified noise dominates and reaches NaN at the pole
        // (atmosphere experience: NaN at 90°N within ~150 iter). Plain copy gives
        // factor 1.0 — same as an axisymmetric-pole assumption to first order.
        // Turbulence scalars (tke, dis, nue, prod, tke_source, dis_source) included
        // so the polar metric singularity does not amplify them either.
        Array* fields_vn[] = {
            &m.t, &m.u, &m.v, &m.w, &m.c,
            &m.Salt_Finger, &m.Salt_Diffusion, &m.Salt_Balance,
            &m.BuoyancyForce, &m.CoriolisForce,
            &m.CentrifugalForce, &m.PresGradForce,
            &m.tke, &m.dis, &m.nue, &m.prod, &m.tke_source, &m.dis_source
        };
        constexpr int n_vn = sizeof(fields_vn) / sizeof(fields_vn[0]);

        const int jml = m.jm - 1;

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int k = 0; k < m.km; k++) {
                for (int f = 0; f < n_vn; f++) {
                    double*** xf = fields_vn[f]->x;
                    xf[i][0][k]   = xf[i][1][k];
                    xf[i][jml][k] = xf[i][jml-1][k];
                }
            }
        }
    }
/*
* 
*/
    // ------------------------------------------------------------------
    // bcPhi — Greenwich meridian boundary conditions (k = 0, k = km-1)
    // All fields: von Neumann extrapolation + periodic averaging
    // ------------------------------------------------------------------
    void bcPhi()
    {
        Array* fields_avg[] = {
            &m.t, &m.u, &m.v, &m.w, &m.c,
            &m.Salt_Finger, &m.Salt_Diffusion, &m.Salt_Balance,
            &m.BuoyancyForce, &m.CoriolisForce,
            &m.CentrifugalForce, &m.PresGradForce
        };
        constexpr int n_avg = sizeof(fields_avg) / sizeof(fields_avg[0]);

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                for (int f = 0; f < n_avg; f++) {
                    double** xij = fields_avg[f]->x[i];
                    double v0   = m.c43 * xij[j][1]      - m.c13 * xij[j][2];
                    double vend = m.c43 * xij[j][m.km-2] - m.c13 * xij[j][m.km-3];
                    xij[j][0] = xij[j][m.km-1] = 0.5 * (v0 + vend);
                }
            }
        }
    }
/*
* 
*/
    // ------------------------------------------------------------------ 
    // bcSolidGround — enforce no-flow on land; extrapolate scalars at
    // land-ocean interfaces in all three coordinate directions
    // ------------------------------------------------------------------
    void bcSolidGround()
    {
        using namespace std;
        cout << endl << "      OGCM: BC_SolidGround" << endl;

        // 6-neighbour water averaging.
        //
        // Replaces the previous (4/3)·x_near − (1/3)·x_far one-sided Neumann
        // extrapolation, which amplified any boundary noise by 4/3 and was
        // direction-dependent (only the first matching side was used).  At steep
        // bathymetric features (mid-ocean ridges, continental shelves) that
        // single-direction stencil could overshoot, producing spurious p_dyn
        // anomalies that propagated into the interior pressure field.  This
        // mirrors the atmosphere fix for bcSolidGround Pass 2 (Himalaya peak).
        //
        // For each land ghost cell, look at the 6 face neighbours (i±1, j±1, k±1
        // with k periodic).  Average whichever of those are water cells, with
        // equal weights.  Amplification factor 1.0, direction-independent, and
        // naturally degrades to plain copy when only one neighbour is water.
        // Fully buried cells (no water neighbour) are left at the Step-1 default
        // state and continue to anchor the Poisson Dirichlet for p_dyn.
        //
        // Non-negative fields (tke, dis, nue, c, t in non-dim form) are wrapped
        // in std::max(0.0, ...) defensively — averaging non-negative values can
        // only produce non-negative results, but the clamp protects against any
        // pre-existing negative value in a contributing neighbour.
        auto average_from_water_neighbors = [&](int i, int j, int k) {
            int kp = (k + 1)        % m.km;
            int kn = (k - 1 + m.km) % m.km;
            bool w_ip = (i + 1 < m.im) && is_water(m.h, i+1, j, k);
            bool w_im = (i - 1 >= 0)    && is_water(m.h, i-1, j, k);
            bool w_jp = (j + 1 < m.jm) && is_water(m.h, i, j+1, k);
            bool w_jm = (j - 1 >= 0)    && is_water(m.h, i, j-1, k);
            bool w_kp = is_water(m.h, i, j, kp);
            bool w_km = is_water(m.h, i, j, kn);

            int count = (int)w_ip + (int)w_im + (int)w_jp + (int)w_jm + (int)w_kp + (int)w_km;
            if (count == 0) return;
            const double inv = 1.0 / static_cast<double>(count);

            auto avg = [&](const Array& f) -> double {
                double sum = 0.0;
                if (w_ip) sum += f.x[i+1][j][k];
                if (w_im) sum += f.x[i-1][j][k];
                if (w_jp) sum += f.x[i][j+1][k];
                if (w_jm) sum += f.x[i][j-1][k];
                if (w_kp) sum += f.x[i][j][kp];
                if (w_km) sum += f.x[i][j][kn];
                return sum * inv;
            };

            // Turbulence scalars (non-negative).
            m.tke.x[i][j][k] = std::max(0.0, avg(m.tke));
            m.dis.x[i][j][k] = std::max(0.0, avg(m.dis));
            m.nue.x[i][j][k] = std::max(0.0, avg(m.nue));

            // Thermodynamic / dynamic scalars.
            m.t.x[i][j][k]            = avg(m.t);
            m.c.x[i][j][k]            = std::max(0.0, avg(m.c));
            m.p_dyn.x[i][j][k]        = avg(m.p_dyn);
            m.p_hydro.x[i][j][k]      = avg(m.p_hydro);
            m.r_water.x[i][j][k]      = avg(m.r_water);
            m.r_salt_water.x[i][j][k] = avg(m.r_salt_water);

            // Body-force diagnostics — averaged so cross-sections show smooth
            // fields across the seafloor instead of a sharp zero step.
            m.BuoyancyForce.x[i][j][k]    = avg(m.BuoyancyForce);
            m.CoriolisForce.x[i][j][k]    = avg(m.CoriolisForce);
            m.CentrifugalForce.x[i][j][k] = avg(m.CentrifugalForce);
            m.PresGradForce.x[i][j][k]    = avg(m.PresGradForce);
        };

        // ---- Step 1: Zero all land cells --------------------------------
        #pragma omp parallel for collapse(3)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                for (int k = 0; k < m.km; k++) {
                    if (is_land(m.h, i, j, k)) {
                        m.u.x[i][j][k]                = 0.0;
                        m.v.x[i][j][k]                = 0.0;
                        m.w.x[i][j][k]                = 0.0;
                        m.un.x[i][j][k]               = 0.0;
                        m.vn.x[i][j][k]               = 0.0;
                        m.wn.x[i][j][k]               = 0.0;
                        m.p_dyn.x[i][j][k]            = 0.0;
                        m.p_hydro.x[i][j][k]          = 0.0;
                        m.r_water.x[i][j][k]          = m.r_0_water;
                        m.r_salt_water.x[i][j][k]     = m.r_0_saltwater;
//                        m.c.x[i][j][k]                = 0.0;
                        m.BuoyancyForce.x[i][j][k]    = m.r_0_saltwater * m.g;
                        m.CoriolisForce.x[i][j][k]    = 0.0;
                        m.CentrifugalForce.x[i][j][k] = 0.0;
                        m.PresGradForce.x[i][j][k]    = 0.0;

                        // Turbulence: zero all scalars at solid surfaces.
                        m.tke.x[i][j][k]              = 0.0;
                        m.tken.x[i][j][k]             = 0.0;
                        m.dis.x[i][j][k]              = 0.0;
                        m.disn.x[i][j][k]             = 0.0;
                        m.nue.x[i][j][k]              = 0.0;
                        m.prod.x[i][j][k]             = 0.0;
                        m.tke_source.x[i][j][k]       = 0.0;
                        m.dis_source.x[i][j][k]       = 0.0;
                    }
                }
            }
        }

        // ---- Step 1b (inviscid spin-up only): free-slip override on surface land cells ---
        // Step 1 above zeroed velocity in every land cell. Here we overwrite the velocity
        // of land cells that border water with a reflecting condition: wall-normal component
        // zero, tangential components mirrored from the adjacent water cell. The choice of
        // "normal" follows the first water neighbour found (priority i+1, then j-faces, then
        // k-faces with Greenwich wrap). One face per cell is enough for spin-up.
        if (m.inviscid_phase) {
            #pragma omp parallel for collapse(3)
            for (int i = 0; i < m.im; i++) {
                for (int j = 0; j < m.jm; j++) {
                    for (int k = 0; k < m.km; k++) {
                        if (!is_land(m.h, i, j, k)) continue;

                        bool water_above = (i + 1 <= m.im - 1) && is_water(m.h, i+1, j, k);
                        if (water_above) {
                            m.u.x[i][j][k] = 0.0;
                            m.v.x[i][j][k] = m.v.x[i+1][j][k];
                            m.w.x[i][j][k] = m.w.x[i+1][j][k];
                        } else if (j > 0 && is_water(m.h, i, j-1, k)) {
                            m.u.x[i][j][k] = m.u.x[i][j-1][k];
                            m.v.x[i][j][k] = 0.0;
                            m.w.x[i][j][k] = m.w.x[i][j-1][k];
                        } else if (j < m.jm-1 && is_water(m.h, i, j+1, k)) {
                            m.u.x[i][j][k] = m.u.x[i][j+1][k];
                            m.v.x[i][j][k] = 0.0;
                            m.w.x[i][j][k] = m.w.x[i][j+1][k];
                        } else {
                            int k_prev = (k > 0)        ? k - 1 : m.km - 1;
                            int k_next = (k < m.km - 1) ? k + 1 : 0;
                            if (is_water(m.h, i, j, k_prev)) {
                                m.u.x[i][j][k] = m.u.x[i][j][k_prev];
                                m.v.x[i][j][k] = m.v.x[i][j][k_prev];
                                m.w.x[i][j][k] = 0.0;
                            } else if (is_water(m.h, i, j, k_next)) {
                                m.u.x[i][j][k] = m.u.x[i][j][k_next];
                                m.v.x[i][j][k] = m.v.x[i][j][k_next];
                                m.w.x[i][j][k] = 0.0;
                            }
                            // else: fully buried — leave at the zero set by Step 1.
                        }
                        m.un.x[i][j][k] = m.u.x[i][j][k];
                        m.vn.x[i][j][k] = m.v.x[i][j][k];
                        m.wn.x[i][j][k] = m.w.x[i][j][k];
                    }
                }
            }
        }

        // ---- Step 2: 6-neighbour water averaging at every land ghost cell ---
        // Replaces direction-biased (4/3, -1/3) Neumann extrapolation.  Buried
        // interior cells (no water neighbour) are left at the Step-1 zero/default
        // state and continue to serve as the Poisson Dirichlet anchor for p_dyn.
        #pragma omp parallel for schedule(static)
        for (int i = 1; i < m.im-1; i++) {
            for (int j = 0; j < m.jm; j++) {
                for (int k = 0; k < m.km; k++) {
                    if (is_land(m.h, i, j, k)) {
                        average_from_water_neighbors(i, j, k);
                    }
                }
            }
        }

        // ---- Step 3: defensive bit-level sanitization of turbulence fields ---
        // Under -ffast-math, std::max / std::clamp are unreliable on NaN inputs
        // (the compiler assumes operands are finite).  A NaN intermediate in
        // `prod` (computed in TurbulenceHyd / RHS_Hyd_Turb near the polar metric
        // singularity or a steep bathymetric corner) could leak through max(0,...)
        // and contaminate the next RK4 cycle.  Use the same bit-level non-finite
        // check as Paraview_Hyd::safe_val to reset non-finite values to physically
        // sane defaults, then clip to physical ceilings.  Last operation here so
        // all fields are bounded by the time the next solver step reads them.
        const double tke_max_phys = 1000.0;
        const double tke_max_nd   = tke_max_phys / (m.u_0 * m.u_0);
        const double nue_max      = 1000.0 / (m.u_0 * m.L_hyd);
        const double prod_max     = 1.0e4;            // [m²/s³] — generous
        const double dis_max      = 1.0e6;            // [m²/s³] — generous for ω*

        auto safe_clamp = [](double v, double lo, double hi) -> double {
            std::uint64_t bits;
            std::memcpy(&bits, &v, sizeof(bits));
            if ((bits & 0x7FF0000000000000ULL) == 0x7FF0000000000000ULL) return lo;
            return (v < lo) ? lo : (v > hi) ? hi : v;
        };

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                for (int k = 0; k < m.km; k++) {
                    m.tke.x[i][j][k]        = safe_clamp(m.tke.x[i][j][k],        0.0,        tke_max_nd);
                    m.tken.x[i][j][k]       = safe_clamp(m.tken.x[i][j][k],       0.0,        tke_max_nd);
                    m.dis.x[i][j][k]        = safe_clamp(m.dis.x[i][j][k],        1.0e-10,    dis_max);
                    m.disn.x[i][j][k]       = safe_clamp(m.disn.x[i][j][k],       1.0e-10,    dis_max);
                    m.nue.x[i][j][k]        = safe_clamp(m.nue.x[i][j][k],        0.0,        nue_max);
                    m.prod.x[i][j][k]       = safe_clamp(m.prod.x[i][j][k],       0.0,        prod_max);
                    m.tke_source.x[i][j][k] = safe_clamp(m.tke_source.x[i][j][k], -prod_max,  prod_max);
                    m.dis_source.x[i][j][k] = safe_clamp(m.dis_source.x[i][j][k], -prod_max,  prod_max);
                }
            }
        }

        cout << "      OGCM: BC_SolidGround ended" << endl;
    }
/*
* 
*/
    // ------------------------------------------------------------------
    // boundaryCondition — damp velocity components adjacent to solid walls
    // ------------------------------------------------------------------
    void boundaryCondition()
    {
        using namespace std;
        cout << endl << "      OGCM: BoundaryCondition" << endl;

        constexpr double coeff   = 0.1;
        constexpr double coeff_5 = 5.0 * coeff;

        #pragma omp parallel for collapse(3)
        for (int i = 1; i < m.im-1; i++) {
            for (int j = 1; j < m.jm-1; j++) {
                for (int k = 1; k < m.km-1; k++) {

                    // i-direction: damp u near vertical land walls
                    if (i < m.im-2 && is_land(m.h, i, j, k) && is_water(m.h, i+1, j, k)) {
                        #pragma omp atomic
                        m.u.x[i+1][j][k] *= coeff;
                        #pragma omp atomic
                        m.u.x[i+2][j][k] *= coeff_5;
                    }

                    // j-direction: damp v near meridional land walls
                    if (j == 0 && j < m.jm-1) {
                        if (is_land(m.h, i, j, k) && is_water(m.h, i, j+1, k) && is_water(m.h, i, j+2, k)) {
                            #pragma omp atomic
                            m.v.x[i][j+1][k] *= coeff;
                            #pragma omp atomic
                            m.v.x[i][j+2][k] *= coeff_5;
                        }
                        if (is_land(m.h, i, j, k) && is_water(m.h, i, j-1, k) && is_water(m.h, i, j-2, k)) {
                            #pragma omp atomic
                            m.v.x[i][j-1][k] *= coeff;
                            #pragma omp atomic
                            m.v.x[i][j-2][k] *= coeff_5;
                        }
                    }

                    // k-direction: damp w near zonal land walls
                    if (k >= 0 && k < m.km-1) {
                        if (is_land(m.h, i, j, k) && is_water(m.h, i, j, k+1) && is_water(m.h, i, j, k+2)) {
                            #pragma omp atomic
                            m.w.x[i][j][k+1] *= coeff;
                            #pragma omp atomic
                            m.w.x[i][j][k+2] *= coeff_5;
                        }
                        if (is_land(m.h, i, j, k) && is_water(m.h, i, j, k-1) && is_water(m.h, i, j, k-2)) {
                            #pragma omp atomic
                            m.w.x[i][j][k-1] *= coeff;
                            #pragma omp atomic
                            m.w.x[i][j][k-2] *= coeff_5;
                        }
                    }
                }
            }
        }
        cout << "      OGCM: BoundaryCondition ended" << endl;
    }

private:
    cHydrosphereModel& m;
};
