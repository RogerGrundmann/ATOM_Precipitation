#pragma once

#include "cHydrosphereModel.h"
#include "Utils.h"

#include <iostream>
#include <algorithm>

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
    // All fields: von Neumann at both poles
    // ------------------------------------------------------------------
    void bcTheta()
    {
        Array* fields_vn[] = {
            &m.t, &m.u, &m.v, &m.w, &m.c,
            &m.Salt_Finger, &m.Salt_Diffusion, &m.Salt_Balance,
            &m.BuoyancyForce, &m.CoriolisForce,
            &m.CentrifugalForce, &m.PresGradForce
        };
        constexpr int n_vn = sizeof(fields_vn) / sizeof(fields_vn[0]);

        const int jml = m.jm - 1;

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int k = 0; k < m.km; k++) {
                for (int f = 0; f < n_vn; f++) {
                    double*** xf = fields_vn[f]->x;
                    xf[i][0][k]   = m.c43 * xf[i][1][k]       - m.c13 * xf[i][2][k];
                    xf[i][jml][k] = m.c43 * xf[i][jml-1][k]   - m.c13 * xf[i][jml-2][k];
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


        // Local c43/c13 as double — the class members are int (truncated),
        // so shadow them explicitly for correct extrapolation here.
        constexpr double c43 = 4.0 / 3.0;
        constexpr double c13 = 1.0 / 3.0;

        auto extrap2 = [&](int i, int j, int k,
                           int i1, int j1, int k1,
                           int i2, int j2, int k2) {
            m.p_dyn.x[i][j][k]   = c43 * m.p_dyn.x[i1][j1][k1]   - c13 * m.p_dyn.x[i2][j2][k2];
            m.p_hydro.x[i][j][k] = c43 * m.p_hydro.x[i1][j1][k1] - c13 * m.p_hydro.x[i2][j2][k2];
//            m.c.x[i][j][k]       = c43 * m.c.x[i1][j1][k1]       - c13 * m.c.x[i2][j2][k2];
        };

        auto copy2 = [&](int i, int j, int k,
                         int i1, int j1, int k1) {
            m.p_dyn.x[i][j][k]   = m.p_dyn.x[i1][j1][k1];
            m.p_hydro.x[i][j][k] = m.p_hydro.x[i1][j1][k1];
//            m.c.x[i][j][k]       = m.c.x[i1][j1][k1];
        };

        // ---- Step 1: Zero all land cells --------------------------------
        #pragma omp parallel for collapse(3)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                for (int k = 0; k < m.km; k++) {
                    if (is_land(m.h, i, j, k)) {
                        m.u.x[i][j][k]                    = 0.0;
                        m.v.x[i][j][k]                    = 0.0;
                        m.w.x[i][j][k]                    = 0.0;
                        m.un.x[i][j][k]                   = 0.0;
                        m.vn.x[i][j][k]                   = 0.0;
                        m.wn.x[i][j][k]                   = 0.0;
                        m.p_dyn.x[i][j][k]                = 0.0;
                        m.p_hydro.x[i][j][k]              = 0.0;
                        m.r_water.x[i][j][k]              = m.r_0_water;
                        m.r_salt_water.x[i][j][k]         = m.r_0_saltwater;
//                        m.c.x[i][j][k]                    = 0.0;
                        m.BuoyancyForce.x[i][j][k]        = m.r_0_saltwater * m.g;
                        m.CoriolisForce.x[i][j][k]        = 0.0;
                        m.CentrifugalForce.x[i][j][k]     = 0.0;
                        m.PresGradForce.x[i][j][k]        = 0.0;
                        // Turbulence: zero all scalars at solid surfaces.
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

        // ---- Step 2: Extrapolate p_dyn, p_hydro and c at land-ocean interfaces ---
        #pragma omp parallel for collapse(3)
        for (int i = 1; i < m.im-1; i++) {
            for (int j = 1; j < m.jm-1; j++) {
                for (int k = 1; k < m.km-1; k++) {

                    const bool land_ijk = is_land(m.h, i, j, k);

                    // i-direction
                    if (i < m.im-3) {
                        if (land_ijk && is_water(m.h, i+1, j, k))
                            extrap2(i,j,k, i+1,j,k, i+2,j,k);
                    } else if (i == m.im-2) {
                        if (land_ijk && is_water(m.h, i+1, j, k))
                            copy2(i,j,k, i+1,j,k);
                    }

                    // j-direction — interior
                    if (j > 2 && j < m.jm-3) {
                        if (land_ijk && is_water(m.h, i, j+1, k) && is_water(m.h, i, j+2, k))
                            extrap2(i,j,k, i,j+1,k, i,j+2,k);
                        if (land_ijk && is_water(m.h, i, j-1, k) && is_water(m.h, i, j-2, k))
                            extrap2(i,j,k, i,j-1,k, i,j-2,k);
                    }

                    // j-direction — poles
                    if (j == 0)        copy2(i,j,k, i,j+1,k);
                    if (j == m.jm-1)   copy2(i,j,k, i,j-1,k);

                    // k-direction — interior
                    if (k > 2 && k < m.km-3) {
                        if (land_ijk && is_water(m.h, i, j, k+1) && is_water(m.h, i, j, k+2))
                            extrap2(i,j,k, i,j,k+1, i,j,k+2);
                        if (land_ijk && is_water(m.h, i, j, k-1) && is_water(m.h, i, j, k-2))
                            extrap2(i,j,k, i,j,k-1, i,j,k-2);
                    }

                    // k-direction — Greenwich
                    if (k == 0)        copy2(i,j,k, i,j,k+1);
                    if (k == m.km-1)   copy2(i,j,k, i,j,k-1);
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
