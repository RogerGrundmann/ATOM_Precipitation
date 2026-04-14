#pragma once

#include "cAtmosphereModel.h"

#include <algorithm>
#include <cmath>
#include <chrono>
#include <iostream>
#include <cstdio>
#include <vector>

using namespace AtomUtils;


namespace OneCatIce {
    // collection/diffusion coefficients
    constexpr double c_ac = 0.24;                                       // m2/kg
    constexpr double c_rim = 0.69;                                      // m2/kg
    constexpr double b_ev = 5.98;                                       // m2*s/kg
    constexpr double a_melt = 3.90e-6;                                  // K/(kg/kg)
    constexpr double b_melt = 10.50;                                    // m2*s/kg
    constexpr double b_dep = 10.50;                                     // m2*s/kg
    constexpr double a_if = 1.92e-6;
    constexpr double a_cf = 3.97e-5;
    constexpr double E_cf = 5.0e-3;
    constexpr double N_cf_0_surf = 2.0e5;                               // 1/m3
    constexpr double N_cf_0_top = 1.0e4;                                // 1/m3
    constexpr double tau_r = 3.3e3;                                     // s, adjusted for NASA avg 2.68 mm/d
    constexpr double tau_s = 1.0e3;                                     // s
    constexpr double a_mc = 0.08;                                       // kg/m2
    constexpr double a_mv = 0.02;                                       // kg/m2
    constexpr int iter_prec_end = 2;                                    // COSMO iterations
}


// One-Category-Ice-Scheme, COSMO-module from the German Weather Forecast,
// resulting the precipitation distribution formed by rain and snow.
class OneCatIceScheme {
public:
    explicit OneCatIceScheme(cAtmosphereModel& model)
        : m(model)
        , t_m1(0.5 * (model.t_0 + model.t_000))
        , t_m2(0.5 * (model.t_0 + model.t_00))
    {}

    void run() {
        using namespace std;
        using namespace OneCatIce;

        cout << endl << endl << endl << "      OneCategoryIceScheme" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        initArrays();
        computeColumns();
        applyBoundaryConditions();
        applyTopography();

        printReport();

        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for OneCategoryIceScheme\n", elapsed.count() * 1e-9);

        cout << "      OneCategoryIceScheme ended" << endl << endl << endl;
    }

private:
    cAtmosphereModel& m;
    const double t_m1;
    const double t_m2;

    // diagnostic output state (set during computeColumns)
    bool rain = false;
    bool snow = false;
    double P_sat = 0.0;
    double P_Snow = 0.0;
    int i_rain = 0, j_radain = 0, k_rain = 0;
    int i_snow = 0, j_snow = 0,   k_snow = 0;
    double height_rain = 0.0;
    double height_snow = 0.0;


    // ==================== INIT ====================
    void initArrays() {
        for(int k = 0; k < m.km; k++){
            for(int j = 0; j < m.jm; j++){
                for(int i = 0; i < m.im; i++){
                    if(m.c.x[i][j][k] < 0.0)       m.c.x[i][j][k] = 0.0;
                    if(m.cloud.x[i][j][k] < 0.0)   m.cloud.x[i][j][k] = 0.0;
                    if(m.ice.x[i][j][k] < 0.0)     m.ice.x[i][j][k] = 0.0;
                    if(m.P_rain.x[i][j][k] < 0.0)  m.P_rain.x[i][j][k] = 0.0;
                    if(m.P_snow.x[i][j][k] < 0.0)  m.P_snow.x[i][j][k] = 0.0;

                    m.P_rain.x[i][j][k] = m.P_rainn.x[i][j][k];
                    m.P_snow.x[i][j][k] = m.P_snown.x[i][j][k];
                }
            }
        }
    }


    // ==================== COLUMN COMPUTATION ====================
    void computeColumns() {
        using namespace std;
        using namespace OneCatIce;

        double maxValue_rain = 0.0;
        double maxValue_snow = 0.0;

        double N_cf_0 = 0.0;
        double eps_t = 0.0;
        double a_m = 0.0;
        double a_ev_local = 0.0;
        double a_dep = 0.0;
        double N_cf = 0.0;

        double q_sat = 0.0, E_sat = 0.0;
        double q_Ice = 0.0,  E_Ice = 0.0;

        double S_nuc, S_frz, S_cf_frz, S_if_frz, S_dep, S_au, S_ac, S_rim, S_shed, S_ev, S_melt;

        std::vector<double> step(m.im, 0.0);

        double Rain_check = 0.0;
        double P_rain_diff = 0.0;

        for(int k = 1; k < m.km-1; k++){
            for(int j = 1; j < m.jm-1; j++){

                for(int iter_prec = 1; iter_prec <= iter_prec_end; iter_prec++){

                    Rain_check = m.P_rain.x[23][j][k];

                    m.P_rain.x[m.im-1][j][k] = 0.0;
                    m.P_snow.x[m.im-1][j][k] = 0.0;

                    for(int i = m.im-2; i >= 0; i--){

                        double Rain = m.P_rain.x[i][j][k];
                        double Snow = m.P_snow.x[i][j][k];

                        double t_u = m.t.x[i][j][k] * m.t_0;

                        step[i] = m.get_layer_height(i+1) - m.get_layer_height(i); // local atmospheric shell thickness in m

                        E_sat = m.hp * AtomUtils::exp_func(t_u, 17.2694, 35.86); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                        q_sat = m.ep * E_sat/(m.p_stat.x[i][j][k] - E_sat);   // relativ water vapour contents on ocean surface reduced by factor in kg/kg

                        E_Ice = m.hp * AtomUtils::exp_func(t_u, 21.8746, 7.66);
                        q_Ice = m.ep * E_Ice/(m.p_stat.x[i][j][k] - E_Ice);      // relativ water vapour contents on ocean surface reduced by factor in kg/kg


                        // mass size relation of circular plates
                        if((t_u < m.t_0)&&(t_u >= t_m1))
                            a_m = a_mc - a_mv * (1.0 + cos(2.0 * M_PI  // kg/m²
                                * (t_u - t_m1)/(m.t_0 - m.t_000)));
                        else  a_m = a_mc;


                        // autoconversion and nucleation process, epsilon(T)
                        if((t_u > m.t_00)&&(t_u <= m.t_0))
                            eps_t = 0.5 * (1.0 + sin(M_PI*(t_m2 - t_u)/(m.t_0 - m.t_00)));
                        else
                            eps_t = 0.0;

                        if(m.cloud.x[i][j][k] > 0.0){
                            S_au  = (1.0 - eps_t)/tau_r * max(0.0, m.cloud.x[i][j][k]);
                            S_nuc = eps_t/tau_s * max(0.0, m.cloud.x[i][j][k]);
                        }


                        // collection mechanisms: accretion, riming, shedding
                        // accretion of cloud water by raindrops
                        S_ac = (1.0 - eps_t) * c_ac * m.cloud.x[i][j][k] // c_ac = 0.24, in m²/kg
                            * pow(Rain,(7.0/9.0));                      // in kg/s

                        if(t_u < m.t_0)
                            S_rim = c_rim/a_m * m.cloud.x[i][j][k] * Snow; // c_rim = 18.6, m²/kg
                        else  S_rim = 0.0;                              // riming rate of snow mass due to collection of supercooled cloud droplets, < VIII >

                        if(t_u >= m.t_0)
                            S_shed = c_rim/a_m * m.cloud.x[i][j][k] * Snow;
                        else  S_shed = 0.0;                             // rate of water shed by melting wet snow particles, < IX >


                        // diffusional growth of rain and snow
                        // evaporation of rain water
                        a_ev_local = 2.76e-3 * exp(0.055 * (m.t_0 - t_u));
                        S_ev = a_ev_local * (1.0 + b_ev * pow(Rain,(1.0/6.0))) // evaporation of rain due to water vapour diffusion, < XIII >
                            * (q_sat - m.c.x[i][j][k])
                            * pow(Rain, (4.0/9.0));

                        // deposition growth and sublimation of cloud ice
                        a_dep = 1.13e-3 * exp(0.073 * (m.t_0 - t_u));
                        S_dep = a_dep/pow(a_m, - 0.5) * (1.0 + b_dep   // melting rate of snow to form rain, < XVI >
                            * pow(a_m, - 0.25) * pow(Snow, (0.9/4.3))) // c_s_melt = 8.43e-5, (m²*s)/(K*kg)
                            * ((m.c.x[i][j][k] - q_Ice) * pow(Snow, (5.0/8.6)));

                        // melting of snow to form cloud water
                        S_melt = a_melt/pow(a_mc, - 0.5) * (1.0 + b_melt // melting rate of snow to form rain, < XVI >
                            * pow(a_mc, - 0.25) * pow(Snow, (0.9/4.3))) // c_s_melt = 8.43e-5, (m²*s)/(K*kg)
                            * ((t_u - m.t_0) * pow(Snow, (5.0/8.6)));


                        // freezing of rain to form snow
                        N_cf_0 = N_cf_0_surf + (N_cf_0_surf - N_cf_0_top)
                            /(m.p_stat.x[0][j][k] - 500.0)
                            * (m.p_stat.x[i][j][k] - m.p_stat.x[0][j][k]);

                        if(m.p_stat.x[i][j][k] <= 500.0)  N_cf_0 = N_cf_0_top;

                        if(t_u < 270.16){
                            N_cf = N_cf_0 * pow((270.16 - t_u), 1.3);
                        }
                        else  N_cf = 0.0;

                        if((t_u <= m.t_0)&&(t_u = m.t_00)){
                            S_if_frz = a_if * (exp(a_if * (t_u - m.t_0)) - 1.0) // immersion freezing
                                * pow(Rain, (14.0/9.0));
                            S_cf_frz = a_cf * E_cf * N_cf                  // contact freezing nucleation
                                * pow(Rain, (13.0/9.0));
                            S_frz = S_if_frz + S_cf_frz;                   // mixing ratio of snow
                        }
                        else  S_frz = 0.0;


                        // sinks and sources
                        m.S_v.x[i][j][k] = - m.S_c_c.x[i][j][k] + S_ev - S_dep; // in kg/(kg*s)
                        m.S_c.x[i][j][k] =   m.S_c_c.x[i][j][k] - S_au - S_ac
                                           - S_nuc - S_rim - S_shed;
                        m.S_r.x[i][j][k] =   S_au + S_ac - S_ev + S_shed
                                           - S_frz + S_melt;
                        m.S_s.x[i][j][k] =   S_nuc + S_rim + S_dep
                                           + S_frz - S_melt;


                        // rain integration
                        if(t_u >= m.t_0)
                              m.P_rain.x[i][j][k] = m.P_rain.x[i+1][j][k]
                                  + m.r_humid.x[i+1][j][k] * m.S_r.x[i+1][j][k]
                                  * step[i];                            // in kg/(m² * s) == mm/s
                        else  m.P_rain.x[i][j][k] = 0.0;

                        if(m.P_rain.x[i][j][k] < 0.0)  m.P_rain.x[i][j][k] = 0.0;

                        if(m.P_rain.x[i][j][k] > 0.0){
                            rain = true;

                            if(m.P_rain.x[i][j][k] > maxValue_rain){
                                maxValue_rain = m.P_rain.x[i][j][k];
                                P_sat = maxValue_rain * 8.64e4;        // in mm/d
                                i_rain = i;
                                j_radain = j;
                                k_rain = k;
                                height_rain = m.get_layer_height(i);
                            }
                        }


                        // snow integration
                        if((t_u < m.t_0)&&(t_u >= m.t_000))
                              m.P_snow.x[i][j][k] = m.P_snow.x[i+1][j][k]
                                  + m.r_humid.x[i+1][j][k] * m.S_s.x[i+1][j][k]
                                  * step[i];                            // in kg/(m² * s) == mm/s
                        else  m.P_snow.x[i][j][k] = 0.0;

                        if(m.P_snow.x[i][j][k] < 0.0)  m.P_snow.x[i][j][k] = 0.0;

                        if(m.P_snow.x[i][j][k] > 0.0){
                            snow = true;
                            if(m.P_snow.x[i][j][k] > maxValue_snow){
                                maxValue_snow = m.P_snow.x[i][j][k];
                                P_Snow = maxValue_snow * 8.64e4;        // in mm/d
                                i_snow = i;
                                j_snow = j;
                                k_snow = k;
                                height_snow = m.get_layer_height(i);
                            }
                        }

                        m.Precipitation.x[i][j][k] =
                            m.P_rain.x[i][j][k] + m.P_snow.x[i][j][k]; // in mm/s

                    }  // end i


                    P_rain_diff = fabs(m.P_rain.x[23][j][k] - Rain_check) * 8.64e4;

                    if(P_rain_diff <= 1.0e-3){
                        std::cout.precision(10);
                        std::cout.setf(std::ios::fixed);
                        if((j == 90)&&(k == 180))  std::cout << std::endl
                            << " .... iteration ended " << std::endl
                            << " .... iter_prec = " << iter_prec
                            << " .... iter_prec_end = " << iter_prec_end << std::endl
                            << " .... Rain_check = " << Rain_check * 8.64e4
                            << " .... P_rain_diff = " << P_rain_diff
                            << " .... P_rain[23] = " << m.P_rain.x[23][j][k] * 8.64e4 << std::endl << std::endl << std::endl;
                        break;
                    }else{
                        std::cout.precision(10);
                        std::cout.setf(std::ios::fixed);
                        if((j == 90)&&(k == 180))  std::cout << std::endl
                            << " .... iter_prec = " << iter_prec
                            << " .... iter_prec_end = " << iter_prec_end << std::endl
                            << " .... Rain_check = " << Rain_check * 8.64e4
                            << " .... P_rain_diff = " << P_rain_diff
                            << " .... P_rain[23] = " << m.P_rain.x[23][j][k] * 8.64e4 << std::endl << std::endl;
                    }
                }  // end iter_prec
            }  // end j
        }  // end k
    }


    // ==================== BOUNDARY CONDITIONS ====================
    void applyBoundaryConditions() {
        // Top boundary (extrapolation)
        #pragma omp parallel for
        for(int j = 0; j < m.jm; j++){
            for(int k = 0; k < m.km; k++){
                m.P_rain.x[m.im-1][j][k] = m.c43 * m.P_rain.x[m.im-2][j][k]
                    - m.c13 * m.P_rain.x[m.im-3][j][k];
                m.P_snow.x[m.im-1][j][k] = m.c43 * m.P_snow.x[m.im-2][j][k]
                    - m.c13 * m.P_snow.x[m.im-3][j][k];
            }
        }

        // Latitude (theta) boundaries
        #pragma omp parallel for
        for(int k = 0; k < m.km; k++){
            for(int i = 0; i < m.im; i++){
                m.P_rain.x[i][0][k] = m.c43 * m.P_rain.x[i][1][k]
                    - m.c13 * m.P_rain.x[i][2][k];
                m.P_rain.x[i][m.jm-1][k] = m.c43 * m.P_rain.x[i][m.jm-2][k]
                    - m.c13 * m.P_rain.x[i][m.jm-3][k];
                m.P_snow.x[i][0][k] = m.c43 * m.P_snow.x[i][1][k]
                    - m.c13 * m.P_snow.x[i][2][k];
                m.P_snow.x[i][m.jm-1][k] = m.c43 * m.P_snow.x[i][m.jm-2][k]
                    - m.c13 * m.P_snow.x[i][m.jm-3][k];
            }
        }

        // Longitude (phi) boundaries – von Neumann + periodicity
        #pragma omp parallel for
        for(int i = 0; i < m.im; i++){
            for(int j = 0; j < m.jm; j++){
                m.P_rain.x[i][j][0] = m.c43 * m.P_rain.x[i][j][1]
                    - m.c13 * m.P_rain.x[i][j][2];                     // von Neumann boundary condition dt/dphi = 0.0
                m.P_rain.x[i][j][m.km-1] = m.c43 * m.P_rain.x[i][j][m.km-2]
                    - m.c13 * m.P_rain.x[i][j][m.km-3];                // von Neumann boundary condition dt/dphi = 0.0
                m.P_rain.x[i][j][0] = m.P_rain.x[i][j][m.km-1]
                    = (m.P_rain.x[i][j][0] + m.P_rain.x[i][j][m.km-1])/2.0;

                m.P_snow.x[i][j][0] = m.c43 * m.P_snow.x[i][j][1]
                    - m.c13 * m.P_snow.x[i][j][2];                     // von Neumann boundary condition dt/dphi = 0.0
                m.P_snow.x[i][j][m.km-1] = m.c43 * m.P_snow.x[i][j][m.km-2]
                    - m.c13 * m.P_snow.x[i][j][m.km-3];                // von Neumann boundary condition dt/dphi = 0.0
                m.P_snow.x[i][j][0] = m.P_snow.x[i][j][m.km-1]
                    = (m.P_snow.x[i][j][0] + m.P_snow.x[i][j][m.km-1])/2.0;
            }
        }
    }


    // ==================== TOPOGRAPHY FILL ====================
    void applyTopography() {
        for(int j = 0; j < m.jm; j++){
            for(int k = 0; k < m.km; k++){
                int i_mount = m.i_topography[j][k];
                for(int i = i_mount; i >= 0; i--){
                    if((is_land(m.h, i, j, k))&&(i < i_mount)){
                        m.P_rain.x[i][j][k] = m.P_rain.x[i_mount][j][k];
                        m.P_snow.x[i][j][k] = m.P_snow.x[i_mount][j][k];

                        m.S_r.x[i][j][k]   = m.S_r.x[i_mount][j][k];
                        m.S_s.x[i][j][k]   = m.S_s.x[i_mount][j][k];
                        m.S_v.x[i][j][k]   = m.S_v.x[i_mount][j][k];
                        m.S_c.x[i][j][k]   = m.S_c.x[i_mount][j][k];
                        m.S_c_c.x[i][j][k] = m.S_c_c.x[i_mount][j][k];
                    }
                }
            }
        }
    }


    // ==================== DIAGNOSTIC REPORT ====================
    void printReport() const {
        using namespace std;

        if(!rain)
            cout << "      no rain fall in OneCategoryIceScheme found" << endl;
        else
            cout << "      rain fall in OneCategoryIceScheme found"
            << endl
            << "      i_rain = " << i_rain
            << "      j_radain = " << j_radain
            << "   k_rain = " << k_rain
            << "   height_rain[m] = " << height_rain
            << "   P_sat[mm/d] = " << P_sat << endl;

        if(!snow)
            cout << "      no snow fall in OneCategoryIceScheme found" << endl;
        else
            cout << "      snow fall in OneCategoryIceScheme found"
            << endl
            << "      i_snow = " << i_snow
            << "      j_snow = " << j_snow
            << "   k_snow = " << k_snow
            << "   height_snow[m] = " << height_snow
            << "   P_Snow[mm/d] = " << P_Snow << endl;
    }
};
