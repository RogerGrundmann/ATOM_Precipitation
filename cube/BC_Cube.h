#pragma once

#include "cCubeModel.h"
#include "Utils.h"

#include <iostream>
#include <algorithm>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace AtomUtils;

class BC_Cube {
public:
    explicit BC_Cube(cCubeModel& model)
        : m(model)
    {}

    // ------------------------------------------------------------------
    void bcSolidGround()
    {
        using namespace std;
//        cout << endl << endl << endl << "      AGCM: BC_SolidGround" << endl;


        auto extrapolation = [&](int i, int j, int k,
                           int i1, int j1, int k1,
                           int i2, int j2, int k2) {
            m.tke.x[i][j][k] = std::max(0.0, m.c43 * m.tke.x[i1][j1][k1] - m.c13 * m.tke.x[i2][j2][k2]);
            m.dis.x[i][j][k] = std::max(0.0, m.c43 * m.dis.x[i1][j1][k1] - m.c13 * m.dis.x[i2][j2][k2]);
            m.nue.x[i][j][k] = std::max(0.0, m.c43 * m.nue.x[i1][j1][k1] - m.c13 * m.nue.x[i2][j2][k2]);

            m.p_dyn.x[i][j][k]   = m.c43 * m.p_dyn.x[i1][j1][k1]   - m.c13 * m.p_dyn.x[i2][j2][k2];
        };

        // Zero-gradient (Neumann) copy of the same six fields from one cell to another.
        // Used at grid boundaries (poles, Greenwich seam, top layer) where a second
        // air neighbour for extrapolation is not available.
        auto copy = [&](int i, int j, int k,
                         int i1, int j1, int k1) {
            m.tke.x[i][j][k] = m.tke.x[i1][j1][k1];
            m.dis.x[i][j][k] = m.dis.x[i1][j1][k1];
            m.nue.x[i][j][k] = m.nue.x[i1][j1][k1];

            m.p_dyn.x[i][j][k] = m.p_dyn.x[i1][j1][k1];
        };

        // Pass 1 — apply solid-wall conditions to every land cell from i=0 up to
        // and including the topographic surface layer i_mount.
        // Two distinct cases are handled in separate branches:
        //   interior : i < i_mount AND no air neighbour in j or k
        //              (fully buried — cannot affect the air column in any direction)
        //   surface  : i == i_mount OR adjacent to air in j or k
        //              (topographic top or cliff face — boundary to the air column)
        #pragma omp parallel for schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                int i_mount = m.i_topography[j][k];

                for (int i = 0; i <= i_mount; i++) {
                    if (!is_land(m.h, i, j, k)) continue;

                    // Check for air neighbours in j (latitude) and k (longitude).
                    // k is periodic: the Greenwich seam (k=0) wraps to k=km-1.
                    bool adj_air_j =
                        (j > 0       && is_air(m.h, i, j-1, k)) ||
                        (j < m.jm-1  && is_air(m.h, i, j+1, k));
                    int k_prev = (k > 0) ? k-1 : m.km-1;
                    int k_next = (k < m.km-1) ? k+1 : 0;
                    bool adj_air_k = is_air(m.h, i, j, k_prev) || is_air(m.h, i, j, k_next);

                    bool on_surface = (i == i_mount) || adj_air_j || adj_air_k;

                    if (!on_surface) {
                        // Fully buried interior cell — zero every prognostic field so
                        // the orography cannot contaminate the overlying air column.

                        // No-slip velocity (current and next time level).
                        m.u.x[i][j][k]  = 0.0;
                        m.v.x[i][j][k]  = 0.0;
                        m.w.x[i][j][k]  = 0.0;
                        m.un.x[i][j][k] = 0.0;
                        m.vn.x[i][j][k] = 0.0;
                        m.wn.x[i][j][k] = 0.0;

                        m.p_dyn.x[i][j][k] = 0.0;

                        // Turbulence: zero all scalars inside solid cells.
                        m.tke.x[i][j][k]        = 0.0;
                        m.tken.x[i][j][k]       = 0.0;
                        m.dis.x[i][j][k]        = 0.0;
                        m.disn.x[i][j][k]       = 0.0;
                        m.nue.x[i][j][k]        = 0.0;
                        m.prod.x[i][j][k]       = 0.0;
                        m.tke_source.x[i][j][k] = 0.0;
                        m.dis_source.x[i][j][k] = 0.0;

                    } else {
                        // Surface (wall) cell — topographic top or cliff face adjacent to air.
                        // No-slip: zero velocity on the wall.
                        // Turbulence scalars and p_dyn are NOT zeroed here; they receive
                        // values from the surrounding flow via the turbulence model and BCs.
                        m.u.x[i][j][k]  = 0.0;
                        m.v.x[i][j][k]  = 0.0;
                        m.w.x[i][j][k]  = 0.0;
                        m.un.x[i][j][k] = 0.0;
                        m.vn.x[i][j][k] = 0.0;
                        m.wn.x[i][j][k] = 0.0;
                    }
                }
            }
        }

        // Pass 2 — extrapolate p_dyn, microphysics, and CO2 from air into land ghost
        // cells at every land/air interface in the j- and k-directions, and at the
        // uppermost land level in the i-direction.  This gives the pressure-gradient
        // and tracer-advection stencils smooth values across the topography boundary.
        #pragma omp parallel for schedule(static)
        for (int i = 1; i < m.im-1; i++) {
            for (int j = 1; j < m.jm-1; j++) {
                for (int k = 1; k < m.km-1; k++) {

                    bool land_ijk = is_land(m.h, i, j, k);

                    // i-direction (radial): extrapolate outward from the first two air
                    // cells above the surface.  At the second-to-last level use copy
                    // because there is only one air cell available above.
                    if (i < m.im-3) {
                        if (land_ijk && is_air(m.h, i+1, j, k))
                            extrapolation(i,j,k, i+1,j,k, i+2,j,k);
                    } else if (i == m.im-2) {
                        if (land_ijk && is_air(m.h, i+1, j, k))
                            copy(i,j,k, i+1,j,k);
                    }

                    // j-direction interior (latitude): extrapolate from whichever
                    // side has two consecutive air cells.
                    if (j > 2 && j < m.jm-3) {
                        if (land_ijk && is_air(m.h, i, j+1, k) && is_air(m.h, i, j+2, k))
                            extrapolation(i,j,k, i,j+1,k, i,j+2,k);
                        if (land_ijk && is_air(m.h, i, j-1, k) && is_air(m.h, i, j-2, k))
                            extrapolation(i,j,k, i,j-1,k, i,j-2,k);
                    }

                    // j-direction poles: single neighbour only, use copy.
                    if (j == 0)        copy(i,j,k, i,j+1,k);
                    if (j == m.jm-1)   copy(i,j,k, i,j-1,k);

                    // k-direction interior (longitude): same logic as j-direction.
                    if (k > 2 && k < m.km-3) {
                        if (land_ijk && is_air(m.h, i, j, k+1) && is_air(m.h, i, j, k+2))
                            extrapolation(i,j,k, i,j,k+1, i,j,k+2);
                        if (land_ijk && is_air(m.h, i, j, k-1) && is_air(m.h, i, j, k-2))
                            extrapolation(i,j,k, i,j,k-1, i,j,k-2);
                    }

                    // k-direction Greenwich seam (k=0 wraps to k=km-1): copy from inner neighbour.
                    if (k == 0)        copy(i,j,k, i,j,k+1);
                    if (k == m.km-1)   copy(i,j,k, i,j,k-1);
                }
            }
        }

        // Pass 3 — copy surface-level (i_mount) values of p_dyn and velocity down
        // to the i=0 reference layer for land/cube columns only.  Air columns
        // (i_topography == 0) have no solid surface at i>0 and must not be
        // touched here; their i=0 boundary is set by bcRadius.
        // Turbulence scalars are only copied when i=0 is a surface cell (borders
        // air in j or k).  Interior land cells at i=0 (e.g. cube body) must stay
        // zero — copying from i_mount would overwrite the zeros set by Pass 1.
        #pragma omp parallel for schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                int i_mount = m.i_topography[j][k];
                if (i_mount == 0) continue;                             // pure air column — skip

                m.p_dyn.x[0][j][k]      = m.p_dyn.x[i_mount][j][k];

                m.u.x[0][j][k]          = m.u.x[i_mount][j][k];
                m.v.x[0][j][k]          = m.v.x[i_mount][j][k];
                m.w.x[0][j][k]          = m.w.x[i_mount][j][k];

                bool adj_air_j_0 =
                    (j > 0      && is_air(m.h, 0, j-1, k)) ||
                    (j < m.jm-1 && is_air(m.h, 0, j+1, k));
                int k_prev = (k > 0) ? k-1 : m.km-1;
                int k_next = (k < m.km-1) ? k+1 : 0;
                bool adj_air_k_0 = is_air(m.h, 0, j, k_prev) || is_air(m.h, 0, j, k_next);

                if (adj_air_j_0 || adj_air_k_0) {
                    m.tke.x[0][j][k]        = m.tke.x[i_mount][j][k];
                    m.dis.x[0][j][k]        = m.dis.x[i_mount][j][k];
                    m.tke_source.x[0][j][k] = m.tke_source.x[i_mount][j][k];
                    m.dis_source.x[0][j][k] = m.dis_source.x[i_mount][j][k];
                    m.nue.x[0][j][k]        = m.nue.x[i_mount][j][k];
                    m.prod.x[0][j][k]       = m.prod.x[i_mount][j][k];
                } else {
                    m.tke.x[0][j][k]        = 0.0;
                    m.tken.x[0][j][k]       = 0.0;
                    m.dis.x[0][j][k]        = 0.0;
                    m.disn.x[0][j][k]       = 0.0;
                    m.nue.x[0][j][k]        = 0.0;
                    m.prod.x[0][j][k]       = 0.0;
                    m.tke_source.x[0][j][k] = 0.0;
                    m.dis_source.x[0][j][k] = 0.0;
                }
            }
        }

        // Pass 4 — clip any negative microphysics values to zero across the full
        // domain.  Negative values can arise from the advection scheme near sharp
        // gradients; clamping here prevents unphysical condensate.
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                for (int k = 0; k < m.km; k++) {
                    if (m.tke.x[i][j][k]   < 0.0)  m.tke.x[i][j][k]   = 0.0;
                    if (m.dis.x[i][j][k]   < 0.0)  m.dis.x[i][j][k]   = 0.0;
/*
                    if (m.tke_source.x[i][j][k]   < 0.0)  m.tke_source.x[i][j][k] = 0.0;
                    if (m.dis_source.x[i][j][k]   < 0.0)  m.dis_source.x[i][j][k] = 0.0;
                    if (m.nue.x[i][j][k]   < 0.0)  m.nue.x[i][j][k]   = 0.0;
                    if (m.prod.x[i][j][k]   < 0.0)  m.prod.x[i][j][k] = 0.0;
*/
                }
            }
        }
//        cout << "      AGCM: BC_SolidGround ended" << endl;
    }
/*
* 
*/
    // ------------------------------------------------------------------
    void bcRadius()
    {
        #pragma omp parallel for schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {

                m.u.x[0][j][k] = 0.0;
                m.v.x[0][j][k] = 0.0;
                m.w.x[0][j][k] = 0.0;

                m.u.x[m.im-1][j][k] = 0.0;
                m.v.x[m.im-1][j][k] = 0.0;
                m.w.x[m.im-1][j][k] = m.w.x[m.im-4][j][k]
                    - 3.0 * m.w.x[m.im-3][j][k] + 3.0 * m.w.x[m.im-2][j][k];

                if (m.use_k_epsilon_turbulence_model) {
                    m.nue.x[0][j][k]        = m.nue.x[3][j][k]        - 3.0 * m.nue.x[2][j][k]        + 3.0 * m.nue.x[1][j][k];
                    m.prod.x[0][j][k]       = m.prod.x[3][j][k]       - 3.0 * m.prod.x[2][j][k]       + 3.0 * m.prod.x[1][j][k];
                    m.tke.x[0][j][k]        = m.tke.x[3][j][k]        - 3.0 * m.tke.x[2][j][k]        + 3.0 * m.tke.x[1][j][k];
                    m.dis.x[0][j][k]        = m.dis.x[3][j][k]        - 3.0 * m.dis.x[2][j][k]        + 3.0 * m.dis.x[1][j][k];
                    m.tke_source.x[0][j][k] = m.tke_source.x[3][j][k] - 3.0 * m.tke_source.x[2][j][k] + 3.0 * m.tke_source.x[1][j][k];
                    m.dis_source.x[0][j][k] = m.dis_source.x[3][j][k] - 3.0 * m.dis_source.x[2][j][k] + 3.0 * m.dis_source.x[1][j][k];

                    m.nue.x[m.im-1][j][k]        = m.c43 * m.nue.x[m.im-2][j][k]        - m.c13 * m.nue.x[m.im-3][j][k];
                    m.prod.x[m.im-1][j][k]       = m.c43 * m.prod.x[m.im-2][j][k]       - m.c13 * m.prod.x[m.im-3][j][k];
                    m.tke.x[m.im-1][j][k]        = m.c43 * m.tke.x[m.im-2][j][k]        - m.c13 * m.tke.x[m.im-3][j][k];
                    m.dis.x[m.im-1][j][k]        = m.c43 * m.dis.x[m.im-2][j][k]        - m.c13 * m.dis.x[m.im-3][j][k];
                    m.tke_source.x[m.im-1][j][k] = m.c43 * m.tke_source.x[m.im-2][j][k] - m.c13 * m.tke_source.x[m.im-3][j][k];
                    m.dis_source.x[m.im-1][j][k] = m.c43 * m.dis_source.x[m.im-2][j][k] - m.c13 * m.dis_source.x[m.im-3][j][k];
                }

                if (m.use_k_omega_turbulence_model) {
                    m.nue.x[0][j][k]        = 0.0;
                    m.prod.x[0][j][k]       = m.prod.x[3][j][k]       - 3.0 * m.prod.x[2][j][k]       + 3.0 * m.prod.x[1][j][k];
                    m.tke.x[0][j][k]        = 0.0;
                    m.dis.x[0][j][k]        = m.dis.x[3][j][k]        - 3.0 * m.dis.x[2][j][k]        + 3.0 * m.dis.x[1][j][k];
                    m.tke_source.x[0][j][k] = m.tke_source.x[3][j][k] - 3.0 * m.tke_source.x[2][j][k] + 3.0 * m.tke_source.x[1][j][k];
                    m.dis_source.x[0][j][k] = m.dis_source.x[3][j][k] - 3.0 * m.dis_source.x[2][j][k] + 3.0 * m.dis_source.x[1][j][k];

                    m.nue.x[m.im-1][j][k]        = m.c43 * m.nue.x[m.im-2][j][k]        - m.c13 * m.nue.x[m.im-3][j][k];
                    m.prod.x[m.im-1][j][k]       = m.c43 * m.prod.x[m.im-2][j][k]       - m.c13 * m.prod.x[m.im-3][j][k];
                    m.tke.x[m.im-1][j][k]        = m.c43 * m.tke.x[m.im-2][j][k]        - m.c13 * m.tke.x[m.im-3][j][k];
                    m.dis.x[m.im-1][j][k]        = m.c43 * m.dis.x[m.im-2][j][k]        - m.c13 * m.dis.x[m.im-3][j][k];
                    m.tke_source.x[m.im-1][j][k] = m.c43 * m.tke_source.x[m.im-2][j][k] - m.c13 * m.tke_source.x[m.im-3][j][k];
                    m.dis_source.x[m.im-1][j][k] = m.c43 * m.dis_source.x[m.im-2][j][k] - m.c13 * m.dis_source.x[m.im-3][j][k];
                }

                if (m.use_k_omega_SST_turbulence_model) {
                    m.nue.x[0][j][k]        = m.nue.x[3][j][k]        - 3.0 * m.nue.x[2][j][k]        + 3.0 * m.nue.x[1][j][k];
                    m.prod.x[0][j][k]       = m.prod.x[3][j][k]       - 3.0 * m.prod.x[2][j][k]       + 3.0 * m.prod.x[1][j][k];
                    m.tke.x[0][j][k]        = m.tke.x[3][j][k]        - 3.0 * m.tke.x[2][j][k]        + 3.0 * m.tke.x[1][j][k];
                    m.dis.x[0][j][k]        = m.dis.x[3][j][k]        - 3.0 * m.dis.x[2][j][k]        + 3.0 * m.dis.x[1][j][k];
                    m.tke_source.x[0][j][k] = m.tke_source.x[3][j][k] - 3.0 * m.tke_source.x[2][j][k] + 3.0 * m.tke_source.x[1][j][k];
                    m.dis_source.x[0][j][k] = m.dis_source.x[3][j][k] - 3.0 * m.dis_source.x[2][j][k] + 3.0 * m.dis_source.x[1][j][k];

                    m.nue.x[m.im-1][j][k]        = m.c43 * m.nue.x[m.im-2][j][k]        - m.c13 * m.nue.x[m.im-3][j][k];
                    m.prod.x[m.im-1][j][k]       = m.c43 * m.prod.x[m.im-2][j][k]       - m.c13 * m.prod.x[m.im-3][j][k];
                    m.tke.x[m.im-1][j][k]        = m.c43 * m.tke.x[m.im-2][j][k]        - m.c13 * m.tke.x[m.im-3][j][k];
                    m.dis.x[m.im-1][j][k]        = m.c43 * m.dis.x[m.im-2][j][k]        - m.c13 * m.dis.x[m.im-3][j][k];
                    m.tke_source.x[m.im-1][j][k] = m.c43 * m.tke_source.x[m.im-2][j][k] - m.c13 * m.tke_source.x[m.im-3][j][k];
                    m.dis_source.x[m.im-1][j][k] = m.c43 * m.dis_source.x[m.im-2][j][k] - m.c13 * m.dis_source.x[m.im-3][j][k];
                }
            }
        }
    }
/*
* 
*/
    // ------------------------------------------------------------------
    void bcTheta()
    {
        #pragma omp parallel for schedule(static)
        for (int k = 0; k < m.km; k++) {
            for (int i = 0; i < m.im; i++) {

                m.u.x[i][0][k]      = m.c43 * m.u.x[i][1][k]      - m.c13 * m.u.x[i][2][k];
                m.u.x[i][m.jm-1][k] = m.c43 * m.u.x[i][m.jm-2][k] - m.c13 * m.u.x[i][m.jm-3][k];

                m.v.x[i][0][k]      = m.c43 * m.v.x[i][1][k]      - m.c13 * m.v.x[i][2][k];
                m.v.x[i][m.jm-1][k] = m.c43 * m.v.x[i][m.jm-2][k] - m.c13 * m.v.x[i][m.jm-3][k];

                m.w.x[i][0][k]      = m.c43 * m.w.x[i][1][k]      - m.c13 * m.w.x[i][2][k];
                m.w.x[i][m.jm-1][k] = m.c43 * m.w.x[i][m.jm-2][k] - m.c13 * m.w.x[i][m.jm-3][k];

                m.tke.x[i][0][k]        = m.c43 * m.tke.x[i][1][k]        - m.c13 * m.tke.x[i][2][k];
                m.dis.x[i][0][k]        = m.c43 * m.dis.x[i][1][k]        - m.c13 * m.dis.x[i][2][k];
                m.nue.x[i][0][k]        = m.c43 * m.nue.x[i][1][k]        - m.c13 * m.nue.x[i][2][k];
                m.prod.x[i][0][k]       = m.c43 * m.prod.x[i][1][k]       - m.c13 * m.prod.x[i][2][k];
                m.tke_source.x[i][0][k] = m.c43 * m.tke_source.x[i][1][k] - m.c13 * m.tke_source.x[i][2][k];
                m.dis_source.x[i][0][k] = m.c43 * m.dis_source.x[i][1][k] - m.c13 * m.dis_source.x[i][2][k];

                m.tke.x[i][m.jm-1][k]        = m.c43 * m.tke.x[i][m.jm-2][k]        - m.c13 * m.tke.x[i][m.jm-3][k];
                m.dis.x[i][m.jm-1][k]        = m.c43 * m.dis.x[i][m.jm-2][k]        - m.c13 * m.dis.x[i][m.jm-3][k];
                m.nue.x[i][m.jm-1][k]        = m.c43 * m.nue.x[i][m.jm-2][k]        - m.c13 * m.nue.x[i][m.jm-3][k];
                m.prod.x[i][m.jm-1][k]       = m.c43 * m.prod.x[i][m.jm-2][k]       - m.c13 * m.prod.x[i][m.jm-3][k];
                m.tke_source.x[i][m.jm-1][k] = m.c43 * m.tke_source.x[i][m.jm-2][k] - m.c13 * m.tke_source.x[i][m.jm-3][k];
                m.dis_source.x[i][m.jm-1][k] = m.c43 * m.dis_source.x[i][m.jm-2][k] - m.c13 * m.dis_source.x[i][m.jm-3][k];
            }
        }
    }
/*
* 
*/
    // ------------------------------------------------------------------
    void bcPhi()
    {
        // false — quadratic Neumann: (4/3)·f[1] − (1/3)·f[2]  (zero-gradient, 2nd-order)
        // true  — cubic Lagrange:   4·f[1] − 6·f[2] + 4·f[3] − f[4]  (zero-curvature, 4th-order)
        constexpr bool cubic_extrap = false;

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                const int K = m.km - 1;

                double* rows[] = {
                    m.u.x[i][j],          m.v.x[i][j],          m.w.x[i][j],
                    m.tke.x[i][j],        m.dis.x[i][j],        m.nue.x[i][j],
                    m.prod.x[i][j],       m.tke_source.x[i][j], m.dis_source.x[i][j]
                };

                for (double* r : rows) {
                    if (cubic_extrap) {
                        r[0] = 4.0*r[1] - 6.0*r[2] + 4.0*r[3] - r[4];
                        r[K] = 4.0*r[K-1] - 6.0*r[K-2] + 4.0*r[K-3] - r[K-4];
                    } else {
                        r[0] = m.c43 * r[1] - m.c13 * r[2];
                        r[K] = m.c43 * r[K-1] - m.c13 * r[K-2];
                    }
                }
            }
        }
    }
/*
* 
*/
    // ------------------------------------------------------------------
    void bcVelSurfSur()
    {
        using namespace std;
//        cout << endl << "      AGCM: BC_vel_surf_sur" << endl;


        constexpr double coeff  = 0.01;                                 // represents the flow better
        constexpr double coeff5 = 0.9;                                  // represents the flow better

        #pragma omp parallel for schedule(static)
        for (int i = 1; i < m.im-1; i++) {
            for (int j = 1; j < m.jm-1; j++) {
                for (int k = 1; k < m.km-1; k++) {

                    bool land_ijk = is_land(m.h, i, j, k);

                    // i-direction (radial: land→air outward)
                    if (i < m.im-2 && land_ijk && is_air(m.h, i+1, j, k)) {
                        m.u.x[i][j][k]   = 0.0;
                        m.u.x[i+1][j][k] = 0.0;

                        m.v.x[i][j][k]    = 0.0;
                        m.v.x[i+1][j][k] *= coeff;
                        m.v.x[i+2][j][k] *= coeff5;

                        m.w.x[i][j][k]    = 0.0;
                    }

                    // j-direction (meridional: land→air south/north)
                    if (land_ijk) {
                        if (j < m.jm-2 && is_air(m.h, i, j+1, k) && is_air(m.h, i, j+2, k)) {  // south
                            m.u.x[i][j][k]    = 0.0;
                            m.u.x[i][j+1][k] *= coeff;
                            m.u.x[i][j+2][k] *= coeff5;

                            m.v.x[i][j][k]    = 0.0;

                            m.w.x[i][j][k]    = 0.0;
                        }

                        if (j >= 2 && is_air(m.h, i, j-1, k) && is_air(m.h, i, j-2, k)) {  // north
                            m.u.x[i][j][k]    = 0.0;
                            m.u.x[i][j-1][k] *= coeff;
                            m.u.x[i][j-2][k] *= coeff5;

                            m.v.x[i][j][k]    = 0.0;

                            m.w.x[i][j][k]    = 0.0;
                        }
                    }

                    // k-direction (zonal: land→air east/west)
                    if (land_ijk) {
                        if (k < m.km-2 && is_air(m.h, i, j, k+1) && is_air(m.h, i, j, k+2)) {  // east
                            m.u.x[i][j][k]    = 0.0;
                            m.u.x[i][j][k+1] *= coeff;
                            m.u.x[i][j][k+2] *= coeff5;

                            m.v.x[i][j][k]    = 0.0;
                            m.v.x[i][j][k+1] *= coeff;
                            m.v.x[i][j][k+2] *= coeff5;

                            m.w.x[i][j][k]    = 0.0;
                            m.w.x[i][j][k+1]  = 0.0;
                        }

                        if (k >= 2 && is_air(m.h, i, j, k-1) && is_air(m.h, i, j, k-2)) {  // west
                            m.u.x[i][j][k]    = 0.0;
                            m.u.x[i][j][k-1] *= coeff;
                            m.u.x[i][j][k-2] *= coeff5;

                            m.v.x[i][j][k]    = 0.0;
                            m.v.x[i][j][k-1] *= coeff;
                            m.v.x[i][j][k-2] *= coeff5;

                            m.w.x[i][j][k]    = 0.0;
                        }
                    }
                }
            }
        }
//        cout << "      AGCM: BC_vel_surf_sur ended" << endl;
    }
/*
*
*/
    // ------------------------------------------------------------------
    // Neumann (zero normal gradient) boundary condition for all scalar
    // 3-D fields at land-air interfaces.  For each interior land cell the
    // value is extrapolated from the two nearest air cells in the normal
    // direction: x[land] = c43*x[n1] - c13*x[n2].  Velocity components
    // from MoistConvection (u_u, v_u, w_u, u_d, v_d, w_d, MC_v, MC_w)
    // are excluded; the main-model velocities are handled by bcVelSurfSur.
    // At corners where multiple directions border air, the vertical (i)
    // direction is preferred.
    void bcScalarSurfSur()
    {
        using namespace std;
//        cout << endl << "      AGCM: BC_scalar_surf_sur" << endl;

        // S-functions (S_v, S_c, S_i, S_r, S_s, S_g, S_c_c), precipitation
        // (P_rain, P_snow, Precipitation), and p_dyn are excluded: bcSolidGround
        // zeros them at land cells and projects x[i_mount] to i=0, so Neumann
        // extrapolation here would only produce spuriously high land-cell values.
        Array* scalars[] = {
            &m.p_stat
        };
        constexpr int ns = sizeof(scalars) / sizeof(scalars[0]);

        #pragma omp parallel for schedule(static)
        for (int i = 1; i < m.im-1; i++) {
            for (int j = 1; j < m.jm-1; j++) {
                for (int k = 1; k < m.km-1; k++) {

                    if (!is_land(m.h, i, j, k)) continue;

                    // ---- i-direction (vertical, preferred at corners) ----
                    if (is_air(m.h, i+1, j, k)) {
                        if (i + 2 < m.im) {
                            for (int f = 0; f < ns; f++) {
                                double*** x = scalars[f]->x;
                                x[i][j][k] = m.c43 * x[i+1][j][k] - m.c13 * x[i+2][j][k];
                            }
                        } else {
                            for (int f = 0; f < ns; f++)
                                scalars[f]->x[i][j][k] = scalars[f]->x[i+1][j][k];
                        }
                        continue;                                    // i preferred: skip j and k
                    }

                    // ---- j-direction ----
                    if (j < m.jm-2 && is_air(m.h, i, j+1, k) && is_air(m.h, i, j+2, k)) {
                        for (int f = 0; f < ns; f++) {
                            double*** x = scalars[f]->x;
                            x[i][j][k] = m.c43 * x[i][j+1][k] - m.c13 * x[i][j+2][k];
                        }
                    } else if (j >= 2 && is_air(m.h, i, j-1, k) && is_air(m.h, i, j-2, k)) {
                        for (int f = 0; f < ns; f++) {
                            double*** x = scalars[f]->x;
                            x[i][j][k] = m.c43 * x[i][j-1][k] - m.c13 * x[i][j-2][k];
                        }
                    }

                    // ---- k-direction ----
                    if (k < m.km-2 && is_air(m.h, i, j, k+1) && is_air(m.h, i, j, k+2)) {
                        for (int f = 0; f < ns; f++) {
                            double*** x = scalars[f]->x;
                            x[i][j][k] = m.c43 * x[i][j][k+1] - m.c13 * x[i][j][k+2];
                        }
                    } else if (k >= 2 && is_air(m.h, i, j, k-1) && is_air(m.h, i, j, k-2)) {
                        for (int f = 0; f < ns; f++) {
                            double*** x = scalars[f]->x;
                            x[i][j][k] = m.c43 * x[i][j][k-1] - m.c13 * x[i][j][k-2];
                        }
                    }
                }
            }
        }

        #pragma omp parallel for schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                int i_mount = m.i_topography[j][k];

                m.p_dyn.x[0][j][k] = m.p_dyn.x[i_mount][j][k];
            }
        }

//        cout << "      AGCM: BC_scalar_surf_sur ended" << endl;
    }

    // ------------------------------------------------------------------
    // Sommerfeld radiation BC at the downstream outflow boundary (k=km-1).
    // Applied every time step after RK4. Advects the boundary value out of the
    // domain at the local flow speed w, removing reflections that Neumann mirrors:
    //   f_new[km-1] = f_old[km-1] - CFL*(f_old[km-1] - f_old[km-2])
    // CFL = w*dt/dx, clamped to [0,1] for stability. Only positive w (true outflow)
    // drives the update; negative w (backflow) leaves the ghost cell unchanged.
    void bcRadiationOutflow()
    {
        const double inv_dx = 1.0 / m.dx;
        const int    kml    = m.km - 1;

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                const double c   = std::max(0.0, m.w.x[i][j][kml-1]);
                const double cfl = std::min(1.0, c * m.dt * inv_dx);

                m.u.x[i][j][kml] -= cfl * (m.u.x[i][j][kml] - m.u.x[i][j][kml-1]);
                m.v.x[i][j][kml] -= cfl * (m.v.x[i][j][kml] - m.v.x[i][j][kml-1]);
                m.w.x[i][j][kml] -= cfl * (m.w.x[i][j][kml] - m.w.x[i][j][kml-1]);

                m.un.x[i][j][kml] = m.u.x[i][j][kml];
                m.vn.x[i][j][kml] = m.v.x[i][j][kml];
                m.wn.x[i][j][kml] = m.w.x[i][j][kml];
            }
        }
    }

    // ------------------------------------------------------------------
    // Rayleigh sponge layer at the upstream inflow (k=0 … k_sponge-1).
    // Relaxes u, v, w toward the background parabolic profile each time step
    // with quadratic strength σ(k) = σ_max · ((k_sponge-k)/k_sponge)² so that
    // σ = σ_max at k=0 and σ = 0 at k = k_sponge.
    void bcSpongeInflow()
    {
        constexpr int    k_sponge = 20;     // cells upstream of cube_k_beg covered by sponge
        constexpr double sigma_max = 0.15;  // max relaxation fraction per time step

        const double d_i_max = double(m.im - 1);

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < m.im; i++) {
            const double x    = double(i) / d_i_max;
            const double w_bg = 2.0*x - x*x;

            for (int j = 0; j < m.jm; j++) {
                for (int k = 0; k < k_sponge; k++) {
                    const double xi    = double(k_sponge - k) / double(k_sponge);
                    const double sigma = sigma_max * xi * xi;

                    m.u.x[i][j][k]  = (1.0 - sigma) * m.u.x[i][j][k];
                    m.v.x[i][j][k]  = (1.0 - sigma) * m.v.x[i][j][k];
                    m.w.x[i][j][k]  = (1.0 - sigma) * m.w.x[i][j][k] + sigma * w_bg;
                    m.un.x[i][j][k] = (1.0 - sigma) * m.un.x[i][j][k];
                    m.vn.x[i][j][k] = (1.0 - sigma) * m.vn.x[i][j][k];
                    m.wn.x[i][j][k] = (1.0 - sigma) * m.wn.x[i][j][k] + sigma * w_bg;
                }
            }
        }
    }

private:
    cCubeModel& m;
};
