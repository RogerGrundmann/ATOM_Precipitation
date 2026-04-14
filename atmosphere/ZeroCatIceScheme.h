#pragma once

#include "cAtmosphereModel.h"

#include <algorithm>
#include <cmath>
#include <chrono>
#include <iostream>
#include <cstdio>
#include <vector>

using namespace AtomUtils;


namespace ZeroCatIce {
    constexpr double tau_r = 3.3e3;                                     // s
    constexpr double tau_r_inv = 1.0 / tau_r;
    constexpr double b_ev = 8.05;                                       // m2*s/kg
    constexpr double c_ac = 0.24;                                       // m2/kg
    constexpr int iter_prec_end = 2;                                    // COSMO iterations
    constexpr double conv_mmd = 8.64e4;                                 // s/d conversion to mm/d
}


// Zero-Category-Ice-Scheme — rain-only warm-phase precipitation scheme.
// No ice, snow or graupel; uses simple autoconversion + accretion + evaporation.
class ZeroCatIceScheme {
public:
    explicit ZeroCatIceScheme(cAtmosphereModel& model)
        : m(model)
    {}

    void run() {
        using namespace std;
        using namespace ZeroCatIce;

        cout << endl << endl << endl << "      ZeroCategoryIceScheme" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        initArrays();
        computeColumns();
        applyBoundaryConditions();
        applyTopography();

        printReport();

        cout << "      ZeroCategoryIceScheme ended" << endl;

        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for ZeroCategoryIceScheme\n", elapsed.count() * 1e-9);
    }

private:
    cAtmosphereModel& m;

    // diagnostic output state (set during computeColumns)
    bool rain = false;
    double global_max_rain = 0.0;
    int global_i_rain = 0, global_j_rain = 0, global_k_rain = 0;

    // layer geometry tables (filled in initArrays)
    std::vector<double> step_table;
    std::vector<double> height_table;


    // ==================== INIT ====================
    void initArrays() {
        // Precompute layer thickness and height tables
        step_table.assign(m.im, 0.0);
        height_table.resize(m.im);
        for(int i = 0; i < m.im; i++){
            height_table[i] = m.get_layer_height(i);
            if(i < m.im - 1)
                step_table[i] = m.get_layer_height(i + 1) - m.get_layer_height(i);
        }

        // Clamping pass
        #pragma omp parallel for collapse(2)
        for(int j = 0; j < m.jm; j++){
            for(int k = 0; k < m.km; k++){
                #pragma omp simd
                for(int i = 0; i < m.im; i++){
                    m.c.x[i][j][k]     = std::max(0.0, m.c.x[i][j][k]);
                    m.cloud.x[i][j][k] = std::max(0.0, m.cloud.x[i][j][k]);
                    m.P_rain.x[i][j][k] = std::max(0.0, m.P_rainn.x[i][j][k]);
                    m.P_snow.x[i][j][k] = std::max(0.0, m.P_snown.x[i][j][k]);
                }
            }
        }
    }


    // ==================== COLUMN COMPUTATION ====================
    void computeColumns() {
        using namespace ZeroCatIce;

        double local_max_rain = 0.0;

        #pragma omp parallel for collapse(2) schedule(static) \
            reduction(max:local_max_rain)
        for(int j = 1; j < m.jm - 1; j++){
            for(int k = 1; k < m.km - 1; k++){

                double local_max = 0.0;

                for(int iter_prec = 1; iter_prec <= iter_prec_end; iter_prec++){

                    double Rain_check = m.P_rain.x[23][j][k];

                    m.P_rain.x[m.im-1][j][k] = 0.0;
                    m.P_snow.x[m.im-1][j][k] = 0.0;

                    for(int i = m.im - 2; i >= 0; i--){

                        double Rain = m.P_rain.x[i][j][k];
                        double t_u  = m.t.x[i][j][k] * m.t_0;
                        double cl_i = m.cloud.x[i][j][k];

                        double E_sat = m.hp * AtomUtils::exp_func(t_u, 17.2694, 35.86);
                        double q_sat = m.ep * E_sat / (m.p_stat.x[i][j][k] - E_sat);

                        // Accretion: cloud water collected by raindrops
                        double S_ac = (t_u >= m.t_0 && Rain > 0.0)
                            ? c_ac * cl_i * pow(Rain, 7.0/9.0)
                            : 0.0;

                        // Autoconversion: cloud water -> rain
                        double S_au = (cl_i > 0.0)
                            ? cl_i * tau_r_inv
                            : 0.0;

                        // Evaporation of rain
                        double S_ev = 0.0;
                        if(t_u >= m.t_0 && Rain > 0.0){
                            double a_ev = 2.76e-3 * exp(0.055 * (m.t_0 - t_u));
                            double Rain_1_6 = pow(Rain, 1.0/6.0);
                            double Rain_4_9 = pow(Rain, 4.0/9.0);
                            S_ev = a_ev * (1.0 + b_ev * Rain_1_6)
                                 * (q_sat - m.c.x[i][j][k]) * Rain_4_9;
                        }

                        // Sinks and sources
                        double S_c_c_ijk = m.S_c_c.x[i][j][k];
                        m.S_v.x[i][j][k] = -S_c_c_ijk + S_ev;
                        m.S_c.x[i][j][k] =  S_c_c_ijk - S_au - S_ac;
                        m.S_r.x[i][j][k] =  S_au + S_ac - S_ev;

                        // Column integration (top-down)
                        if(t_u >= m.t_0)
                            m.P_rain.x[i][j][k] = std::max(0.0,
                                m.P_rain.x[i+1][j][k]
                                + m.r_humid.x[i+1][j][k] * m.S_r.x[i+1][j][k] * step_table[i]);
                        else
                            m.P_rain.x[i][j][k] = 0.0;

                        if(m.P_rain.x[i][j][k] > local_max)
                            local_max = m.P_rain.x[i][j][k];

                        m.Precipitation.x[i][j][k] = m.P_rain.x[i][j][k];

                    } // end i

                    double P_rain_diff = fabs(m.P_rain.x[23][j][k] - Rain_check) * conv_mmd;
                    if(P_rain_diff <= 1.0e-3) break;

                } // end iter_prec

                if(local_max > local_max_rain)
                    local_max_rain = local_max;

            } // end k
        } // end j

        global_max_rain = local_max_rain;

        // Find exact location of global max (single pass, outside hot loop)
        if(global_max_rain > 0.0){
            rain = true;
            #pragma omp parallel for collapse(2)
            for(int j = 1; j < m.jm - 1; j++){
                for(int k = 1; k < m.km - 1; k++){
                    for(int i = m.im - 2; i >= 0; i--){
                        if(m.P_rain.x[i][j][k] == global_max_rain){
                            #pragma omp critical
                            {
                                global_i_rain = i;
                                global_j_rain = j;
                                global_k_rain = k;
                            }
                        }
                    }
                }
            }
        }
    }


    // ==================== BOUNDARY CONDITIONS ====================
    void applyBoundaryConditions() {
        // Top boundary (extrapolation)
        #pragma omp parallel for collapse(2)
        for(int j = 0; j < m.jm; j++){
            for(int k = 0; k < m.km; k++){
                m.P_rain.x[m.im-1][j][k] = m.c43 * m.P_rain.x[m.im-2][j][k]
                    - m.c13 * m.P_rain.x[m.im-3][j][k];
                m.P_snow.x[m.im-1][j][k] = m.c43 * m.P_snow.x[m.im-2][j][k]
                    - m.c13 * m.P_snow.x[m.im-3][j][k];
            }
        }

        // Latitude (theta) boundaries
        #pragma omp parallel for collapse(2)
        for(int k = 0; k < m.km; k++){
            for(int i = 0; i < m.im; i++){
                m.P_rain.x[i][0][k]       = m.c43 * m.P_rain.x[i][1][k]       - m.c13 * m.P_rain.x[i][2][k];
                m.P_rain.x[i][m.jm-1][k]  = m.c43 * m.P_rain.x[i][m.jm-2][k] - m.c13 * m.P_rain.x[i][m.jm-3][k];
                m.P_snow.x[i][0][k]       = m.c43 * m.P_snow.x[i][1][k]       - m.c13 * m.P_snow.x[i][2][k];
                m.P_snow.x[i][m.jm-1][k]  = m.c43 * m.P_snow.x[i][m.jm-2][k] - m.c13 * m.P_snow.x[i][m.jm-3][k];
            }
        }

        // Longitude (phi) boundaries – von Neumann + periodicity
        #pragma omp parallel for collapse(2)
        for(int i = 0; i < m.im; i++){
            for(int j = 0; j < m.jm; j++){
                m.P_rain.x[i][j][0]      = m.c43 * m.P_rain.x[i][j][1]      - m.c13 * m.P_rain.x[i][j][2];
                m.P_rain.x[i][j][m.km-1] = m.c43 * m.P_rain.x[i][j][m.km-2] - m.c13 * m.P_rain.x[i][j][m.km-3];
                m.P_rain.x[i][j][0] = m.P_rain.x[i][j][m.km-1]
                    = (m.P_rain.x[i][j][0] + m.P_rain.x[i][j][m.km-1]) * 0.5;
                m.P_snow.x[i][j][0]      = m.c43 * m.P_snow.x[i][j][1]      - m.c13 * m.P_snow.x[i][j][2];
                m.P_snow.x[i][j][m.km-1] = m.c43 * m.P_snow.x[i][j][m.km-2] - m.c13 * m.P_snow.x[i][j][m.km-3];
            }
        }
    }


    // ==================== TOPOGRAPHY FILL ====================
    void applyTopography() {
        #pragma omp parallel for collapse(2)
        for(int j = 0; j < m.jm; j++){
            for(int k = 0; k < m.km; k++){
                int i_mount = m.i_topography[j][k];
                double pr_m  = m.P_rain.x[i_mount][j][k];
                double ps_m  = m.P_snow.x[i_mount][j][k];
                double sr_m  = m.S_r.x[i_mount][j][k];
                double ss_m  = m.S_s.x[i_mount][j][k];
                double sv_m  = m.S_v.x[i_mount][j][k];
                double sc_m  = m.S_c.x[i_mount][j][k];
                double scc_m = m.S_c_c.x[i_mount][j][k];

                for(int i = i_mount - 1; i >= 0; i--){
                    if(is_land(m.h, i, j, k)){
                        m.P_rain.x[i][j][k] = pr_m;
                        m.P_snow.x[i][j][k] = ps_m;
                        m.S_r.x[i][j][k]    = sr_m;
                        m.S_s.x[i][j][k]    = ss_m;
                        m.S_v.x[i][j][k]    = sv_m;
                        m.S_c.x[i][j][k]    = sc_m;
                        m.S_c_c.x[i][j][k]  = scc_m;
                    }
                }
            }
        }
    }


    // ==================== DIAGNOSTIC REPORT ====================
    void printReport() const {
        using namespace std;
        using namespace ZeroCatIce;

        if(!rain)
            cout << "      no rain fall in ZeroCategoryIceScheme found" << endl;
        else
            cout << "      rain fall in ZeroCategoryIceScheme found" << endl
                 << "      i_rain = " << global_i_rain
                 << "      j_rain = " << global_j_rain
                 << "   k_rain = " << global_k_rain
                 << "   height_rain[m] = " << height_table[global_i_rain]
                 << "   P_sat[mm/d] = " << global_max_rain * conv_mmd << endl;
    }
};
