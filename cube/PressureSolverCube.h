#pragma once

#include "cCubeModel.h"
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

class PressureSolverCube {
public:
    explicit PressureSolverCube(cCubeModel& model)
        : m(model)
    {}

    void run()
    {
        using namespace std;
/*
        cout << endl << endl << endl << "      ATOM: PressureSolverAtm" << endl;

        auto begin = std::chrono::high_resolution_clock::now();
*/
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
        const double inv_2dz = 1.0 / (2.0 * m.dr);
        const double inv_2dy = 1.0 / (2.0 * m.dy);
        const double inv_2dx = 1.0 / (2.0 * m.dx);
        const double inv_dz2 = 1.0 / (m.dr * m.dr);
        const double inv_dy2 = 1.0 / (m.dy * m.dy);
        const double inv_dx2 = 1.0 / (m.dx * m.dx);
        const double inv_dy  = 1.0 / m.dy;
        const double inv_dx  = 1.0 / m.dx;

        // Main compute loop — land mask lookups + k sliding window
        #pragma omp parallel for collapse(2) schedule(dynamic, 4)
        for (int i = 1; i < m.im-1; i++) {
            for (int j = 1; j < m.jm-1; j++) {

                // Poisson coefficients — uniform Cartesian.
                const double denom = 2.0 * inv_dz2
                                   + 2.0 * inv_dy2
                                   + 2.0 * inv_dx2;
                const double inv_denom = 1.0 / denom;
                const double num1 = inv_dz2;
                const double num2 = inv_dy2;
                const double num3 = inv_dx2;

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

                    double du_dz, dv_dy, dw_dx;
                    bool z_flag   = false;
                    bool y_flag = false;
                    bool x_flag = false;

                    // z direction
                    if (i_in_range && lnd_ijk && !LAND(i+1,j,k)) {
                        du_dz = (-3.0 * m.aux_u.x[i][j][k] + 4.0 * m.aux_u.x[i+1][j][k]
                                 - m.aux_u.x[i+2][j][k]) * inv_2dz;
                        z_flag = true;
                    }

                    // y direction
                    if (j_inner) {
                        const int8_t air_jp1 = !LAND(i,j+1,k);
                        const int8_t air_jm1 = !LAND(i,j-1,k);

                        if (lnd_ijk && air_jp1 && !LAND(i,j+2,k)) {
                            dv_dy = (-3.0 * m.aux_v.x[i][j][k] + 4.0 * m.aux_v.x[i][j+1][k]
                                       - m.aux_v.x[i][j+2][k]) * inv_2dy;
                            y_flag = true;
                        }
                        if (lnd_ijk && air_jm1 && !LAND(i,j-2,k)) {
                            dv_dy = -(-3.0 * m.aux_v.x[i][j][k] + 4.0 * m.aux_v.x[i][j-1][k]
                                        - m.aux_v.x[i][j-2][k]) * inv_2dy;
                            y_flag = true;
                        }
                        if ((lnd_ijk && air_jp1 && LAND(i,j+2,k))
                            || (j == m.jm-2 && !lnd_ijk && LAND(i,j+1,k))) {
                            dv_dy = (m.aux_v.x[i][j+1][k] - m.aux_v.x[i][j][k]) * inv_dy;
                            y_flag = true;
                        }
                        if ((lnd_ijk && air_jm1 && LAND(i,j-2,k))
                            || (j == 1 && lnd_ijk && air_jm1)) {
                            dv_dy = (m.aux_v.x[i][j-1][k] - m.aux_v.x[i][j][k]) * inv_dy;
                            y_flag = true;
                        }
                    }

                    // x direction
                    const bool k_inner = (k > 2) && (k < m.km-2);
                    if (k_inner) {
                        const int8_t air_kp1 = !lnd_kp1;
                        const int8_t air_km1 = !lnd_km1;

                        if (lnd_ijk && air_kp1 && !lnd_kp2) {
                            dw_dx = (-3.0*m.aux_w.x[i][j][k] + 4.0*m.aux_w.x[i][j][k+1]
                                       - m.aux_w.x[i][j][k+2]) * inv_2dx;
                            x_flag = true;
                        }
                        if (lnd_ijk && air_km1 && !lnd_km2) {
                            dw_dx = -(-3.0*m.aux_w.x[i][j][k] + 4.0*m.aux_w.x[i][j][k-1]
                                        - m.aux_w.x[i][j][k-2]) * inv_2dx;
                            x_flag = true;
                        }
                        if ((lnd_ijk && air_kp1 && lnd_kp2)
                            || (k == m.km-2 && !lnd_ijk && lnd_kp1)) {
                            dw_dx = (m.aux_w.x[i][j][k+1] - m.aux_w.x[i][j][k]) * inv_dx;
                            x_flag = true;
                        }
                        if ((lnd_ijk && air_km1 && lnd_km2)
                            || (k == 1 && lnd_ijk && air_km1)) {
                            dw_dx = (m.aux_w.x[i][j][k-1] - m.aux_w.x[i][j][k]) * inv_dx;
                            x_flag = true;
                        }
                    }

                    // central-difference fallbacks
                    if (!z_flag)
                        du_dz   = (m.aux_u.x[i+1][j][k] - m.aux_u.x[i-1][j][k]) * inv_2dz;
                    if (!y_flag)
                        dv_dy = (m.aux_v.x[i][j+1][k] - m.aux_v.x[i][j-1][k]) * inv_2dy;
                    if (!x_flag)
                        dw_dx = (m.aux_w.x[i][j][k+1] - m.aux_w.x[i][j][k-1]) * inv_2dx;

                    // Interior solid cell — no air neighbour in any direction:
                    // p_dyn stays zero inside the cube body.
                    if (lnd_ijk && !z_flag && !y_flag && !x_flag) {
                        m.p_dyn.x[i][j][k] = 0.0;
                        continue;
                    }

                    // pressure update
                    m.p_dyn.x[i][j][k] =
                        ((m.p_dyn.x[i+1][j][k] + m.p_dyn.x[i-1][j][k]) * num1
                       + (m.p_dyn.x[i][j+1][k] + m.p_dyn.x[i][j-1][k]) * num2
                       + (m.p_dyn.x[i][j][k+1] + m.p_dyn.x[i][j][k-1]) * num3
                       - du_dz
                       - dv_dy
                       - dw_dx) * inv_denom;
                } // k
            } // j
        } // i

        #undef LAND

        // z boundary extrapolation
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

        // y boundary extrapolation
        #pragma omp parallel for collapse(2)
        for (int k = 0; k < m.km; k++) {
            for (int i = 0; i < m.im; i++) {
                m.p_dyn.x[i][0][k]      = m.c43 * m.p_dyn.x[i][1][k]      - m.c13 * m.p_dyn.x[i][2][k];
                m.p_dyn.x[i][m.jm-1][k] = m.c43 * m.p_dyn.x[i][m.jm-2][k] - m.c13 * m.p_dyn.x[i][m.jm-3][k];
            }
        }

        // x boundary extrapolation
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                m.p_dyn.x[i][j][0]      = m.c43 * m.p_dyn.x[i][j][1]      - m.c13 * m.p_dyn.x[i][j][2];
                m.p_dyn.x[i][j][m.km-1] = m.c43 * m.p_dyn.x[i][j][m.km-2] - m.c13 * m.p_dyn.x[i][j][m.km-3];
            }
        }
/*
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for PressureSolverCube\n", elapsed.count() * 1e-9);
        cout << "      ATOM: PressureSolverCube ended" << endl;
*/
    }

private:
    cCubeModel& m;
};
