#pragma once

#include "cCubeModel.h"
#include "Utils.h"

#include <iostream>
#include <cmath>
#include <chrono>
#include <cstdio>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace AtomUtils;

class UtilsCube {
public:
    explicit UtilsCube(cCubeModel& model)
        : m(model)
    {}

    // ------------------------------------------------------------------
    void storeIntermediateData3D(float coeff = 1.0f)
    {
        const double cf = coeff;
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                #pragma omp simd
                for (int k = 0; k < m.km; k++) {
                    m.un.x[i][j][k]     = cf * m.u.x[i][j][k];
                    m.vn.x[i][j][k]     = cf * m.v.x[i][j][k];
                    m.wn.x[i][j][k]     = cf * m.w.x[i][j][k];
                    m.tken.x[i][j][k]   = cf * m.tke.x[i][j][k];
                    m.disn.x[i][j][k]   = cf * m.dis.x[i][j][k];
                }
            }
        }
    }

    // ------------------------------------------------------------------
    void findResiduumAtm()
    {
        using namespace std;
/*
        cout << endl << "      AGCM: find_residuum_atm" << endl;

        auto begin = std::chrono::high_resolution_clock::now();
*/
        struct ErrorState { double val; int i, j, k; };
        ErrorState global_max = {0.0, 0, 0, 0};
        double global_sum = 0.0;
        long   global_cnt = 0;

        #pragma omp parallel
        {
            ErrorState local_max = {0.0, 0, 0, 0};
            double local_sum = 0.0;
            long   local_cnt = 0;

            #pragma omp for collapse(2) schedule(static)
            for (int j = 1; j < m.jm - 1; j++) {
                for (int k = 1; k < m.km - 1; k++) {

                    for (int i = 10; i <= 30; i++) {
                        double dudz = (m.u.x[i+1][j][k] - m.u.x[i-1][j][k])
                                      / (2.0 * m.dr);
                        double dvdy = (m.v.x[i][j+1][k] - m.v.x[i][j-1][k])
                                      / (2.0 * m.dy);
                        double dwdx = (m.w.x[i][j][k+1] - m.w.x[i][j][k-1])
                                      / (2.0 * m.dx);

                        double res = sqrt((dudz * dudz
                                        + dvdy * dvdy
                                        + dwdx * dwdx) / 3.0);

                        local_sum += res;
                        local_cnt++;

                        if(res > local_max.val) {
                            local_max = {res, i, j, k};
                            m.residuum_old = res;
                        }
                   }
                }
            }

            #pragma omp critical
            {
                if (local_max.val > global_max.val)
                    global_max = local_max;
                global_sum += local_sum;
                global_cnt += local_cnt;
            }
        }

        const double avg_abs = (global_cnt > 0) ? global_sum / global_cnt : 0.0;
        const double avg_rel = (global_max.val > 0.0) ? avg_abs / global_max.val : 0.0;

        cout.precision(8);
        const bool declining = (m.residuum_old - global_max.val) > 0.0;
        cout << endl
             << (declining
                 ? "      AGCM: find_residuum_atm, absolute error declining .......................\n"
                 : "      AGCM: find_residuum_atm, absolute error is too high .....................\n")
             << "      residuum_atm = " << global_max.val
             << "      residuum_old = " << m.residuum_old
             << "      eps_residuum = " << m.eps_residuum << endl << endl
             << "      i_error = " << global_max.i
             << "   j_error = "    << global_max.j
             << "   k_error = "    << global_max.k << endl << endl
             << "      absolute error = " << fabs(m.residuum_old - global_max.val) << endl
             << "      relative error = " << fabs(global_max.val / m.residuum_old - 1.0) << endl << endl
             << "      error location: y = " << m.the.z[global_max.j]
             << "   x = " << m.phi.z[global_max.k]
             << "   height = " << global_max.i * 400 << " m" << endl
             << "      i_topography at error = " << m.i_topography[global_max.j][global_max.k]
             << "   surface height = " << m.get_layer_height(m.i_topography[global_max.j][global_max.k]) << " m" << endl
             << "      h.x[i][j][k] = " << m.h.x[global_max.i][global_max.j][global_max.k]
             << "   (0=air, 1=land)" << endl << endl
             << "      grid average absolute error = " << avg_abs
             << "   (" << global_cnt << " cells)" << endl
             << "      grid average relative error = " << avg_rel
             << "   (avg/max)" << endl << endl;
/*
        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for find_residuum_atm\n", elapsed.count() * 1e-9);
        cout << "      AGCM: find_residuum_atm ended" << endl;
*/
    }

    // ------------------------------------------------------------------
    void valueLimitationAtm()
    {
        using namespace std;
//        cout << endl << endl << endl << "      ValueLimitationAtm" << endl;

        for (int k = 0; k < m.km; k++) {
            for (int j = 0; j < m.jm; j++) {
                for (int i = 0; i < m.im; i++) {
/*
                    if (m.u.x[i][j][k] >=  0.5)   m.u.x[i][j][k] =  0.5;
                    if (m.u.x[i][j][k] <= -0.5)   m.u.x[i][j][k] = -0.5;
                    if (m.v.x[i][j][k] >=  0.5)   m.v.x[i][j][k] =  0.5;
                    if (m.v.x[i][j][k] <= -0.5)   m.v.x[i][j][k] = -0.5;
                    if (m.w.x[i][j][k] >=  6.25)  m.w.x[i][j][k] =  6.25;
                    if (m.w.x[i][j][k] <= -2.0)   m.w.x[i][j][k] = -2.0;
*/
                    if (m.p_stat.x[i][j][k] <= 0.0)  m.p_stat.x[i][j][k] = 0.0;

                }
            }
        }

//        cout << "      ValueLimitationAtm ended" << endl;
    };

    // ------------------------------------------------------------------
    void writeFile(std::string& bathymetry_name, std::string& output_path, bool is_final_result)
    {
        using namespace std;
//        cout << endl << endl << endl << "      AGCM: write_file" << endl;

        // Skip VTK slices on the final-result call when iter_n==0 (Ma==0 case):
        // both writeFile calls would write to the same filename; the second
        // truncates the file before writing, so a failure there leaves 0 bytes.
        // For Ma>0 the 3D loop advances iter_n, so the filenames differ and
        // neither call is skipped.

        if (!is_final_result || m.iter_n != 0) {
            int i_radial = 3;                                           // 3 above sea level
            m.paraview_vtk_radial(bathymetry_name, i_radial, m.iter_n);

            int j_longal = m.jm / 2+3;                                    // mid-plane (j)
            m.paraview_vtk_longal(bathymetry_name, j_longal, m.iter_n);

            int k_zonal = 85;                                           // middle of the cube
            m.paraview_vtk_zonal(bathymetry_name, k_zonal, m.iter_n);
        }

        if (m.paraview_panorama_vts_flag && m.panorama_cnt == m.panorama_print) {
            m.paraview_panorama_vts(bathymetry_name, m.iter_n);
        }

//        cout << endl << "      AGCM: write_file ended " << endl;
    }

    // ------------------------------------------------------------------
    void resetArrays()
    {
        using namespace std;
//        cout << endl << endl << endl << "      AGCM: reset_arrays" << endl;

        // 1D arrays
        m.rad.initArray_1D(m.im, 1.0);
        m.the.initArray_1D(m.jm, 2.0);
        m.phi.initArray_1D(m.km, 3.0);
        m.aux_grad_v.initArray_1D(m.im, 4.0);

        // 2D arrays — all initialized to 0.0
        Array_2D* arrays_2d[] = {
            &m.Topography, 
        };
        // vel_star is initialised separately — must not be zero when TurbulenceAtm::init()
        // reads it before the first compute_vel_star() call.
        m.vel_star.initArray_2D(m.jm, m.km, 0.4);

        #pragma omp parallel for schedule(static)
        for (size_t n = 0; n < sizeof(arrays_2d) / sizeof(arrays_2d[0]); n++)
            arrays_2d[n]->initArray_2D(m.jm, m.km, 0.0);

        // 3D arrays initialized to 0.0
        Array* arrays_3d_zero[] = {
            &m.h, &m.u, &m.v, &m.w, 
            &m.un, &m.vn, &m.wn, 
            &m.p_dyn, &m.p_stat,
            &m.rhs_u, &m.rhs_v, &m.rhs_w, 
            &m.aux_u, &m.aux_v, &m.aux_w, &m.aux_t,
            // turbulence arrays
            &m.tke, &m.tken, &m.dis, &m.disn,
            &m.nue, &m.prod, &m.tke_source, &m.dis_source,
            &m.rhs_tke, &m.rhs_dis,
        };

        #pragma omp parallel for schedule(static)
        for (size_t n = 0; n < sizeof(arrays_3d_zero) / sizeof(arrays_3d_zero[0]); n++)
            arrays_3d_zero[n]->initArray(m.im, m.jm, m.km, 0.0);


        // 3D arrays initialized to 1.0 (extend as needed)
        Array* arrays_3d_one[] = {};
        (void)arrays_3d_one;   // suppress unused-variable warning when empty

        // integer 2D arrays
        for (auto& i : m.i_topography)
            std::fill(i.begin(), i.end(), 0);
        for (auto& i : m.i_tropopause)
            std::fill(i.begin(), i.end(), 0);

//        cout << "      AGCM: reset_arrays ended" << endl;
    }

    // ------------------------------------------------------------------
    void initCubeContour()
    {
        using namespace std;
//        cout << endl << "      ATOM_cube: init_cube_contour" << endl;

        #pragma omp parallel for collapse(3) schedule(static)
        for (int i = 0; i < m.im; i++)
            for (int j = 0; j < m.jm; j++)
                for (int k = 0; k < m.km; k++)
                    m.h.x[i][j][k] = 0.0;

        #pragma omp parallel for collapse(3) schedule(static)
        for (int k = m.cube_k_beg; k <= m.cube_k_end; k++)
            for (int j = m.cube_j_beg; j <= m.cube_j_end; j++)
                for (int i = m.cube_i_beg; i <= m.cube_i_end; i++)
                    m.h.x[i][j][k] = 1.0;

//        cout << "      ATOM_cube: init_cube_contour ended" << endl;
    }

    // ------------------------------------------------------------------
    void initVelocities()
    {
        using namespace std;
//        cout << endl << "      AGCM: init_velocities" << endl;
        cout << endl << " .. initVelocities:  extention of the flow field .. "
             << m.im << " ..... " << m.jm << " ..... " << m.km << endl << endl;

        const double d_i_max = double(m.im - 1);

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < m.im; i++) {
            const double x = double(i) / d_i_max;
            for (int j = 0; j < m.jm; j++) {
                #pragma omp simd
                for (int k = 0; k < m.km; k++) {
                    m.u.x[i][j][k] = 0.0;
                    m.v.x[i][j][k] = 0.0;
                    m.w.x[i][j][k] = 2.0 * x - x * x;
                }
            }
        }
//        cout << "      AGCM: init_velocities ended" << endl;
    }

private:
    cCubeModel& m;
};
