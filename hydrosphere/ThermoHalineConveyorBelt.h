/*
 * Ocean General Circulation Modell (OGCM) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 *
 * class to prepare the boundary and initial conditions for diverse variables
*/
#pragma once

#include "cHydrosphereModel.h"

#include <chrono>
#include <cstdio>
#include <iostream>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

// ============================================================================
// ThermoHalineConveyorBelt — friend class of cHydrosphereModel
//
// Provides initial conditions for thermohaline circulation:
//   IC_v_w_WestEastCoast()    — v/w velocity near east/west coasts
//   IC_Atlantischer_Ozean()   — Atlantic ocean currents
//   IC_Indischer_Ozean()      — Indian ocean currents
//   IC_Pazifischer_Ozean()    — Pacific ocean currents
//   IC_Nord_Polar_Meer()      — Arctic ocean currents
//   IC_South_Polar_Sea()      — Southern ocean currents
//   IC_EquatorialCurrents()   — equatorial current systems
//   IC_DeepWater()            — deep-water thermohaline conveyor belt
// ============================================================================

class ThermoHalineConveyorBelt {
public:
    explicit ThermoHalineConveyorBelt(cHydrosphereModel& model)
        : i_max   (model.im - 1)            // sea surface (top layer)
        , i_beg   (model.im / 2)            // == 100 m depth, Ekman-layer base
        , i_middle(model.im / 2)            // lower bound of coast downwelling loop
        , i_half  (model.im / 2 + 10)       // location of u-velocity maximum
        , i_bottom(0)                        // ocean floor
        , i_deep  (model.im / 4)            // deep-flow reference level
        , d_i_half(static_cast<double>(model.im / 2 + 10))
        , m(model)
    {}

    // configuration — derived from model in constructor, may be overridden by caller
    int    i_max;
    int    i_beg;
    int    i_middle;
    int    i_half;
    int    i_bottom;
    int    i_deep;
    double IC_water = 1.0;
    double ca       = 0.0;
    double ca_max   = 1.0;
    double d_i_half;

    // working state — written and read within each method call
    int    k_grad = 0, k_a = 0, k_b = 0, k_water = 0, k_sequel = 0, flip = 0;
    int    j_beg = 0, j_end = 0, k_beg = 0, k_end = 0;
    int    j_run = 0, k_run = 0, j_step = 0, k_step = 0;
    int    k_exp = 0, k_w = 0, j_z = 0, j_n = 0, k_z = 0, k_n = 0, l = 0;
    double v_grad = 0.0, d_i = 0.0;
    int    mi = 0;   // formerly 'm' — renamed to avoid clash with model reference
    int    i_EIC_u = 0, i_EIC_o = 0, i_SCC_u = 0, i_SCC_o = 0;
    int    i_ECC_u = 0, i_ECC_o = 0;
    int    k1 = 0, k2 = 0, k3 = 0, kd = 0, kn = 0;
    int    j1 = 0, j2 = 0, j3 = 0, jd = 0, jn = 0;
/*
* 
*/
    // ------------------------------------------------------------------
    void run() {
        using namespace std;

        cout << endl << endl << "      ThermoHalineConveyorBelt  class starts" << endl;
        auto begin = std::chrono::high_resolution_clock::now();

        IC_v_w_WestEastCoast    (m.h, m.u, m.v, m.w, m.c);
        IC_Atlantischer_Ozean   (m.h, m.u, m.v, m.w, m.c);
        IC_Indischer_Ozean      (m.h, m.u, m.v, m.w);
        IC_Pazifischer_Ozean    (m.h, m.u, m.v, m.w);
        IC_Nord_Polar_Meer      (m.h, m.u, m.v, m.w);
        IC_South_Polar_Sea      (m.h, m.u, m.v, m.w, m.c);
        IC_EquatorialCurrents   (m.h, m.u, m.v, m.w);
        IC_DeepWater            (m.h, m.u, m.v, m.w, m.c);

        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for ThermoHalineConveyorBelt  class ended\n",
               elapsed.count() * 1e-9);
    }
/*
* 
*/
    void IC_v_w_WestEastCoast (Array &h, Array &u, Array &v, Array &w, Array &c) {
    // initial conditions for v and w velocity components at the sea surface close to east or west coasts
    // reversal of v velocity component between north and south equatorial current ommitted at respectively 10°
    // w component unchanged
    // search for east coasts and associated velocity components to close the circulations
    // transition between coast flows and open sea flows included
    // northern hemisphere: east coast

    //  k_grad = 10;                                                      // extension of velocity change
        k_grad = 8;                                                     // extension of velocity change
        k_a = k_grad;                                                   // left distance
        k_b = 0;                                                        // right distance

        k_water = 0;                                                    // on water closest to coast
        k_sequel = 1;                                                   // on solid ground

        for (int j = 0; j < 91; j++) {                                  // outer loop: latitude
            for (int k = 0; k < m.km - k_grad; k++) {                   // inner loop: longitude
                if (h.x[i_max][j][k] == 1.) k_sequel = 0;               // if solid ground: k_sequel = 0

                if ((h.x[i_max][j][k] == 0.) && (k_sequel == 0)) k_water = 0; // if water and and k_sequel = 0 then is water closest to coast
                else k_water = 1;                                       // somewhere on water

                if ((h.x[i_max][j][k] == 0.) && (k_water == 0)) {       // if water is closest to coast, change of velocity components begins
                    for (int l = 0; l < k_grad; l++) {                  // extension of change, sign change in v-velocity and distribution of u-velocity with depth
                        v.x[i_max][j][std::max(0, std::min(k + l, m.km - 1))] = - v.x[i_max][j][std::max(0, std::min(k + l, m.km - 1))];  // existing velocity changes sign

                        for (int i = i_middle; i <= i_half; i++) {      // loop in radial direction, extension for u -velocity component, downwelling here
                            mi = i_half - i;
                            d_i = static_cast<double>(i);
                            c.x[i][j][k] = ca_max;
                            u.x[i][j][std::max(0, std::min(k + l, m.km - 1))] = - 10. * d_i / d_i_half * IC_water / (static_cast<double>(l + 1));  // increase with depth, decrease with distance from coast
                            u.x[mi][j][std::max(0, std::min(k + l, m.km - 1))] = - 10. * d_i / d_i_half * IC_water / (static_cast<double>(l + 1));  // decrease with depth, decrease with distance from coast
                        }
                    }

/*
                    for (int l = (k + k_grad - k_a); l < (k + k_grad + k_b + 1) && l >= 0 && l < m.km; l++) {  // starting at local longitude + max extension - begin of smoothing k_a  until ending at  + k_b
                        v.x[i_max][j][l] = (v.x[i_max][j][std::max(0, std::min(k + k_grad + k_b, m.km - 1))] - v.x[i_max][j][std::max(0, std::min(k + k_grad - k_a, m.km - 1))]) / static_cast<double>((k + k_grad + k_b) -  (k + k_grad - k_a)) * static_cast<double>(l -  (k + k_grad - k_a)) + v.x[i_max][j][std::max(0, std::min(k + k_grad - k_a, m.km - 1))]; // extension of v-velocity, smoothing algorithm by a linear equation
                    }

                    for (int l = k; l < (k + k_grad + k_b + 1) && l >= 0 && l < m.km; l++) {  // smoothing algorithm by a linear equation, starting at local longitude until ending at max extension + k_b
                        w.x[i_max][j][l] = w.x[i_max][j][std::max(0, std::min(k + k_grad + k_b, m.km - 1))]  / static_cast<double>((k + k_grad + k_b) -  k) * static_cast<double>(l - k); // extension of v-velocity
                    }
*/
                    k_sequel = 1;                                       // looking for another east coast
                }
            }                                                           // end of longitudinal loop
            k_water = 0;                                                // starting at another latitude
        }                                                               // end of latitudinal loop


    // southern hemisphere: east coast
        k_water = 0;
        k_sequel = 1;

        for (int j = 91; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                if (h.x[i_max][j][k] == 1.) k_sequel = 0;

                if ((h.x[i_max][j][k] == 0.) && (k_sequel == 0)) k_water = 0;
                else k_water = 1;

                if ((h.x[i_max][j][k] == 0.) && (k_water == 0)) {
                    for (int l = 0; l < k_grad; l++) {
                        v.x[i_max][j][std::max(0, std::min(k + l, m.km - 1))] = - v.x[i_max][j][std::max(0, std::min(k + l, m.km - 1))];

                        for (int i = i_middle; i <= i_half; i++) {
                            mi = i_half - i;
                            d_i = static_cast<double>(i);
                            c.x[i][j][k] = ca_max;
                            u.x[i][j][std::max(0, std::min(k + l, m.km - 1))] = - 10. * d_i / d_i_half * IC_water / (static_cast<double>(l + 1)); // increase with depth, decrease with distance from coast
                            u.x[mi][j][std::max(0, std::min(k + l, m.km - 1))] = - 10. * d_i / d_i_half * IC_water / (static_cast<double>(l + 1)); // decrease with depth, decrease with distance from coast
                        }
                    }
/*
                    for (int l = (k + k_grad - k_a); l < (k + k_grad + k_b + 1) && l >= 0 && l < m.km; l++) {
                        v.x[i_max][j][l] = (v.x[i_max][j][std::max(0, std::min(k + k_grad + k_b, m.km - 1))] - v.x[i_max][j][std::max(0, std::min(k + k_grad - k_a, m.km - 1))]) / static_cast<double>((k + k_grad + k_b) -  (k + k_grad - k_a)) * static_cast<double>(l -  (k + k_grad - k_a)) + v.x[i_max][j][std::max(0, std::min(k + k_grad - k_a, m.km - 1))];
                    }

                    for (int l = k; l < (k + k_grad + k_b + 1) && l >= 0 && l < m.km; l++) {
                        w.x[i_max][j][l] = w.x[i_max][j][std::max(0, std::min(k + k_grad + k_b, m.km - 1))]  / static_cast<double>((k + k_grad + k_b) -  k) * static_cast<double>(l - k);
                    }
*/
                    k_sequel = 1;
                }
            }
            k_water = 0;
        }


    // search for west coasts and associated velocity components to close the circulations
    // transition between coast flows and open sea flows included
    // northern hemisphere: west coast
    //  k_grad = 10;                                                      // extension of velocity change
        k_grad = 8;                                                     // extension of velocity change
        k_a = 0;                                                        // left distance

        k_water = 0;                                                    // somewhere on water
        flip = 0;                                                       // somewhere on water

            #pragma omp parallel for private(k_water, flip, mi, d_i)
        for (int j = 0; j < 91; j++) {  // outer loop: latitude
            for (int k = k_grad; k < m.km; k++) {  // inner loop: longitude
                if (h.x[i_max][j][k] == 0.) {  // if somewhere on water
                    k_water = 0;                                        // somewhere on water: k_water = 0
                    flip = 0;                                           // somewhere on water: flip = 0
                }
                else k_water = 1;                                       // first time on land

                if ((flip == 0) && (k_water == 1)) {  // on water closest to land
                    for (int l = k; l > (k - k_grad - 1) && l >= 0 && l < m.km; l--) {        // backward extention of velocity change: nothing changes
                        w.x[i_max][j][l] = - w.x[i_max][j][l];

                        for (int i = i_middle; i <= i_half; i++) {      // loop in radial direction, extension for u -velocity component, upwelling here
                            mi = i_half - i;
                            d_i = static_cast<double>(i);
                            c.x[i][j][k] = ca_max;
                            u.x[i][j][l] = + 10. * d_i / d_i_half * IC_water / (static_cast<double>(k - l + 1)); // increase with depth, decrease with distance from coast
                            u.x[mi][j][l] = + 10. * d_i / d_i_half * IC_water / (static_cast<double>(k - l + 1)); // decrease with depth, decrease with distance from coast
                        }
                    }
/*
                    for (int l = k; l > (k - k_grad - k_a + 1) && l >= 0 && l < m.km; l--) {  // smoothing algorithm by a linear equation, starting at local longitude until ending at max extension + k_b
                        v.x[i_max][j][l] = v.x[i_max][j][std::max(0, std::min(k - k_grad - k_a, m.km - 1))] / static_cast<double>((k - k_grad - k_a) - k) * static_cast<double>(l - k); // extension of v-velocity
                    }

                    for (int l = (k - k_grad - 3); l < (k - k_grad + 3) && l >= 0 && l < m.km; l++) {  // smoothing algorithm by a linear equation, starting at local longitude until ending at max extension + k_b
                        w.x[i_max][j][l] = (- w.x[i_max][j][std::max(0, std::min(k - k_grad - 3, m.km - 1))] + w.x[i_max][j][std::max(0, std::min(k - k_grad + 3, m.km - 1))]) * static_cast<double>(l - (k - k_grad - 3)) / static_cast<double>((k - k_grad + 3) - (k - k_grad - 3)) - w.x[i_max][j][std::max(0, std::min(k - k_grad + 3, m.km - 1))];
                    }
*/
                    flip = 1;
                }
            }
            flip = 0;
        }


    // southern hemisphere: west coast
        k_water = 0;
        flip = 0;

        #pragma omp parallel for private(k_water, flip, mi, d_i)
        for (int j = 91; j < m.jm; j++) {
            for (int k = k_grad; k < m.km; k++) {
                if (h.x[i_max][j][k] == 0.) {
                    k_water = 0;
                    flip = 0;
                }
                else k_water = 1;

                if ((flip == 0) && (k_water == 1)) {
                    for (int l = k; l > (k - k_grad + 1) && l >= 0 && l < m.km; l--) {
                        w.x[i_max][j][l] = - w.x[i_max][j][l];

                        for (int i = i_middle; i <= i_half; i++) {
                            mi = i_half - i;
                            d_i = static_cast<double>(i);
                            c.x[i][j][k] = ca_max;
                            u.x[i][j][l] = + 10. * d_i / d_i_half * IC_water / (static_cast<double>(k - l + 1));
                            u.x[mi][j][l] = + 10. * d_i / d_i_half * IC_water / (static_cast<double>(k - l + 1));
                        }
                    }
/*
                    for (int l = k; l > (k - k_grad - k_a - 1) && l >= 0 && l < m.km; l--) {
                        v.x[i_max][j][l] = v.x[i_max][j][std::max(0, std::min(k - k_grad - k_a, m.km - 1))] / static_cast<double>((k - k_grad - k_a) - k) * static_cast<double>(l - k);
                    }

                    for (int l = (k - k_grad - 3); l < (k - k_grad + 3) && l >= 0 && l < m.km; l++) {  // smoothing algorithm by a linear equation, starting at local longitude until ending at max extension + k_b
                        w.x[i_max][j][l] = (- w.x[i_max][j][std::max(0, std::min(k - k_grad - 3, m.km - 1))] + w.x[i_max][j][std::max(0, std::min(k - k_grad + 3, m.km - 1))]) * static_cast<double>(l - (k - k_grad - 3)) / static_cast<double>((k - k_grad + 3) - (k - k_grad - 3)) - w.x[i_max][j][std::max(0, std::min(k - k_grad + 3, m.km - 1))];
                    }
*/
                    flip = 1;
                }
            }
            flip = 0;
        }
    }
/*
* 
*/
    void IC_Atlantischer_Ozean (Array &h, Array &u, Array &v, Array &w, Array &c) {
    // Currents along the coasts
    // Closing the polar, subpolar and subtropical atmospheric circulation systems
    // Thermohaline Conveyor Belt
    // Atlantic
    // South Atlantic      diagonal from Kap Agulhas in South Africa until Kap St. Roque in South America
    // South Equatorial Current as surface flow (from j=96 until j=112 compares to 6°S until 22°S,
    //                                               from k=325 until k=m.km compares to 35°W until 0°)
        j_beg = 96;
        j_end = 113;
        k_beg = 325;
        k_end = m.km;
        j_run = 0;
        k_run = 0;
        j_step = 2;
        k_step = 10;

        while ((j_beg + j_run) <= j_end && (k_beg + k_run + k_step) <= (k_end + k_step)) {
            for (int j = j_beg + j_run; j < (j_beg + j_step + j_run) && j < m.jm; j++) {
                for (int k = k_beg + k_run; k < (k_beg + k_step + k_run) && k < m.km && k >= 0; k++) {
                    for (int i = i_beg; i < m.im; i++) {
                        if (h.x[i][j][k] == 0.) {
                            v.x[i][j][k] = - IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                            w.x[i][j][k] = - IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        }
                    }
                }
            k_run++;
            }
        j_run++;
        }

    // west of the constant velocity strip diagonal from Kap Agulhas in South Africa until Kap St. Roque in South America
    // South Equatorial Current as surface flow (from j=96 until j=112 compares to 6°S until 22°S,
    //                                               from k=325 until k=m.km compares to 35°W until 0°)
        j_beg = 96;
        j_end = 113;
        k_beg = 325;
        k_end = m.km;
        j_run = 0;
        k_run = 0;
        j_step = 2;
        k_step = 10;
        k_exp = 10;

        while ((j_beg + j_run) <= j_end && (k_beg + k_run - k_exp) <= (k_end - k_exp)) {
            for (int j = j_beg + j_run; j < (j_beg + j_step + j_run) && j < m.jm; j++) {
                for (int k = (k_beg + k_run - k_exp); k < (k_beg + k_run) && k < m.km && k >= 0; k++) {
                    k_z = k - (k_beg + k_run - k_exp);
                    k_n = (k_beg + k_run) - (k_beg + k_run - k_exp);

                    for (int i = i_beg; i < m.im; i++) {
                        if (h.x[i][j][k] == 0.) {
                            v.x[i][j][k] = (v.x[i][j][std::min(k_beg + k_run, m.km - 1)] - v.x[i][j][std::max(k_beg + k_run - k_exp, 0)]) * static_cast<double>(k_z) / static_cast<double>(k_n) + v.x[i][j][std::max(k_beg + k_run - k_exp, 0)];
                            w.x[i][j][k] = (w.x[i][j][std::min(k_beg + k_run, m.km - 1)] - w.x[i][j][std::max(k_beg + k_run - k_exp, 0)]) * static_cast<double>(k_z) / static_cast<double>(k_n) + w.x[i][j][std::max(k_beg + k_run - k_exp, 0)];
                        }
                    }
                }
                k_run++;
            }
        j_run++;
        }

    // east of the constant velocity strip diagonal from Kap Agulhas in South Africa until Kap St. Roque in South America
    // South Equatorial Current as surface flow (from j=96 until j=112 compares to 6°S until 22°S,
    //                                                                              from k=325 until k=m.km compares to 35°W until 0°)
        j_beg = 96;
        j_end = 113;
        k_beg = 325;
        k_end = m.km;
        j_run = 0;
        k_run = 0;
        j_step = 2;
        k_step = 10;
        k_exp = 3;

        while ((j_beg + j_run) <= j_end && (k_beg + k_run + k_step + k_exp) <= (k_end + k_exp)) {
            for (int j = j_beg + j_run; j < (j_beg + j_step + j_run) && j < m.jm; j++) {
                for (int k = (k_beg + k_run + k_step); k < (k_beg + k_run + k_step + k_exp) && k < m.km && k >= 0; k++) {
                    if (k >= k_end) break;
                    k_z = k - (k_beg + k_run + k_step);
                    k_n = (k_beg + k_run + k_step + k_exp) - (k_beg + k_run + k_step) ;

                    for (int i = i_beg; i < m.im; i++) {
                        if (h.x[i][j][k] == 0.) {
                            v.x[i][j][k] = (v.x[i][j][std::min(k_beg + k_run + k_step + k_exp, m.km - 1)] - v.x[i][j][std::min(k_beg + k_run + k_step, m.km - 1)]) * static_cast<double>(k_z) / static_cast<double>(k_n) + v.x[i][j][std::min(k_beg + k_run + k_step, m.km - 1)];
                            w.x[i][j][k] = (w.x[i][j][std::min(k_beg + k_run + k_step + k_exp, m.km - 1)] - w.x[i][j][std::min(k_beg + k_run + k_step, m.km - 1)]) * static_cast<double>(k_z) / static_cast<double>(k_n) + w.x[i][j][std::min(k_beg + k_run + k_step, m.km - 1)];
                        }
                    }
                }
                k_run++;
            }
        j_run++;
        }


    // Thermohalin Conveyor Belt
    // Atlantic
    // South Atlantic      diagonal from Kap Agulhas in South Africa until Kap St. Roque in South America
    // South Equatorial Current as surface flow (from j=117 until j=128 compares to 27°S until 38°S,
    //                                                                               from k=0 until k=20 compares to 0° until 20°O)
        j_beg = 117;
        j_end = 129;
        k_beg = 0;
        k_end = 21;
        j_run = 0;
        k_run = 0;
        j_step = 2;
        k_step = 10;

        while ((j_beg + j_run) <= j_end && (k_beg + k_run) <= k_end) {
            for (int j = j_beg + j_run; j < (j_beg + j_step + j_run) && j < m.jm; j++) {
                for (int k = k_beg + k_run; k < (k_beg + k_step + k_run) && k < m.km && k >= 0; k++) {
                    if (k >= k_end) break;
                    for (int i = i_beg; i < m.im; i++) {
                        if (h.x[i][j][k] == 0.) {
                            v.x[i][j][k] = - IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                            w.x[i][j][k] = - IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        }
                    }
                }
            k_run++;
            }
        j_run++;
        }

    // west of the constant velocity strip diagonal from Kap Agulhas in South Africa until Kap St. Roque in South America
    // South Equatorial Current as surface flow (from j=117 until j=134 compares to 27°S until 44°S,
    //                                                                               from k=0 until k=20 compares to 0° until 20°O)
        j_beg = 117;
        j_end = 135;
        k_beg = 0;
        k_end = 21;
        j_run = 0;
        k_run = 0;
        j_step = 2;
        k_step = 10;
        k_exp = 10;

        while ((j_beg + j_run) <= j_end && (k_beg + k_run - k_exp) <= (k_end - k_exp)) {
            for (int j = j_beg + j_run; j < (j_beg + j_step + j_run) && j < m.jm; j++) {
                for (int k = (k_beg + k_run - k_exp); k < (k_beg + k_run) && k < m.km && k >= 0; k++) {
                    if (k >= k_end) break;
                    k_z = k - (k_beg + k_run - k_exp);
                    k_n = (k_beg + k_run) - (k_beg + k_run - k_exp);

                    for (int i = i_beg; i < m.im; i++) {
                        if (h.x[i][j][k] == 0.) {
                            v.x[i][j][k] = (v.x[i][j][std::min(k_beg + k_run, m.km - 1)] - v.x[i][j][std::max(k_beg + k_run - k_exp, 0)]) * static_cast<double>(k_z) / static_cast<double>(k_n) + v.x[i][j][std::max(k_beg + k_run - k_exp, 0)];
                            w.x[i][j][k] = (w.x[i][j][std::min(k_beg + k_run, m.km - 1)] - w.x[i][j][std::max(k_beg + k_run - k_exp, 0)]) * static_cast<double>(k_z) / static_cast<double>(k_n) + w.x[i][j][std::max(k_beg + k_run - k_exp, 0)];
                        }
                    }
                }
                k_run++;
            }
        j_run++;
        }

    // east of the constant velocity strip diagonal from Kap Agulhas in South Africa until Kap St. Roque in South America
    // South Equatorial Current as surface flow (from j=117 until j=134 compares to 27°S until 44°S,
    //                                                                               from k=0 until k=20 compares to 0° until 20°O)
        j_beg = 117;
        j_end = 113;
        k_beg = 0;
        k_end = 21;
        j_run = 0;
        k_run = - 1;
        j_step = 2;
        k_step = 0;
        k_exp = 6;

        while ((j_beg + j_run) <= j_end && (k_beg + k_run + k_step + k_exp) <= (k_end + k_exp)) {
            for (int j = j_beg + j_run; j < (j_beg + j_step + j_run) && j < m.jm; j++) {
                for (int k = (k_beg + k_run + k_step); k < (k_beg + k_run + k_step + k_exp) && k < m.km && k >= 0; k++) {
                    if (k >= k_end) break;
                    k_z = k - (k_beg + k_run + k_step);
                    k_n = (k_beg + k_run + k_step + k_exp) - (k_beg + k_run + k_step) ;

                    for (int i = i_beg; i < m.im; i++) {
                        if (h.x[i][j][k] == 0.) {
                            v.x[i][j][k] = (v.x[i][j][std::min(k_beg + k_run + k_step + k_exp, m.km - 1)] - v.x[i][j][std::min(k_beg + k_run + k_step, m.km - 1)]) * static_cast<double>(k_z) / static_cast<double>(k_n) + v.x[i][j][std::min(k_beg + k_run + k_step, m.km - 1)];
                            w.x[i][j][k] = (w.x[i][j][std::min(k_beg + k_run + k_step + k_exp, m.km - 1)] - w.x[i][j][std::min(k_beg + k_run + k_step, m.km - 1)]) * static_cast<double>(k_z) / static_cast<double>(k_n) + w.x[i][j][std::min(k_beg + k_run + k_step, m.km - 1)];
                        }
                    }
                }
                k_run++;
            }
        j_run++;
        }


    // Atlantic
    // Coastal currents
    // South America       North coasts
    // Guyana Current in the north-east of South America (from j=71 until j=96 compares to 6°S until 19°N,
    //                                                                                  from k=280 until k=325 compares to 80°W until 35°W)
        j_beg = 71;
        j_end = 97;
        k_beg = 280;
        k_end = 326;

        k_step = 30;

        v_grad = + 0.001;

        k_a = k_b = 0;

        flip = 0;

        for (int k = std::min(k_end, m.km - 1); k > k_beg; k--) {
            if (h.x[i_max][j_beg][k] == 1.) {
                k_a = k;
                flip = 1;
            }
            if (flip == 1) break;
        }
        flip = 0;

        for (int j = j_beg+1; j < j_end; j++) {
            for (int k = std::min(k_end, m.km - 1); k > k_beg; k--) {
                if (h.x[i_max][j][k] == 1.) {
                    k_b = k;
                    k_grad = k_a - k_b;
                    if (k_grad >= 2) k_grad = 1;
                    if (k_grad <= 0) k_grad = 1;
                    flip = 1;
                }
            if (flip == 1) break;
            }


            for (int k = k_b; k < (k_b + k_step) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg; i < m.im; i++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = - v_grad * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = - v_grad * static_cast<double>(k_grad) * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }
            k_a = k_b;
            flip = 0;

            for (int k = k_b + 8; k < (k_b + k_step) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg; i < m.im; i++) {
                    v.x[i][j][k] = (v.x[i][j][std::max(0, std::min(k_b + k_step +1, m.km - 1))] - v.x[i][j][std::max(0, std::min(k_b + 8, m.km - 1))]) * static_cast<double>(k - (k_b + 8)) / static_cast<double>((k_b + k_step -1) - (k_b + 8)) + v.x[i][j][std::max(0, std::min(k_b + 8, m.km - 1))];
                    w.x[i][j][k] = (w.x[i][j][std::max(0, std::min(k_b + k_step +1, m.km - 1))] - w.x[i][j][std::max(0, std::min(k_b + 8, m.km - 1))]) * static_cast<double>(k - (k_b + 8)) / static_cast<double>((k_b + k_step -1) - (k_b + 8)) + w.x[i][j][std::max(0, std::min(k_b + 8, m.km - 1))];
                }
            }
        }


    // Atlantic Ocean
    // Gulf of Mexico                                     (from j=53 until j=58 compares to 23°N until 27°N,
    // parallel to the equator eastward                 from k=234 until k=252 compares to 77°W until 96°W)
        #pragma omp parallel for collapse(3)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = 53; j < 59; j++) {
                for (int k = 234; k < 253; k++) {
                    if (h.x[i][j][k] == 0.) {
                        {
                            v.x[i][j][k] = + 0.0;
                            w.x[i][j][k] = + IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        }
                    }
                }
            }
        }


    // Atlantic Ocean
    // Gulf of Mexico                                          (from j=63 until j=67 compares to 23°N until 27°N,
    // parallel to the equator westward                    from k=264 until k=283 compares to 77°W until 96°W)
        #pragma omp parallel for collapse(3)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = 63; j < 68; j++) {
                for (int k = 264; k < 284; k++) {
                    if (h.x[i][j][k] == 0.) {
                        {
                            v.x[i][j][k] = + 0.0;
                            w.x[i][j][k] = - IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        }
                    }
                }
            }
        }


    // Atlantic Ocean
    // Gulf of Mexico                                          (from j=63 until j=67 compares to 23°N until 27°N,
    // parallel to the equator northward                    from k=264 until k=283 compares to 77°W until 96°W)
        #pragma omp parallel for collapse(3)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = 63; j < 68; j++) {
                for (int k = 264; k < 284; k++) {
                    if (h.x[i][j][k] == 0.) {
                        {
                            v.x[i][j][k] = - IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                            w.x[i][j][k] = - 0.00;
                        }
                    }
                }
            }
        }


    // Atlantic
    // Connection east coast Gulf of Mexico
    // Gulf of Mexico                                     (from j=63 until j=67 compares to 23°N until 27°N,
    // parallel to the equator east coast                 from k=264 until k=283 compares to 77°W until 96°W)
        j_beg = 63;
        j_end = 68;
        k_beg = 264;
        k_end = 284;

        #pragma omp parallel for collapse(3) private(k_z, k_n)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    k_z = k - k_beg;
                    k_n = k_end - k_beg;

                    if (h.x[i][j][k] == 0.) {
//                        u.x[i][j][k] = (u.x[i][j][k_end] - u.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + u.x[i][j][k_beg];
                        v.x[i][j][k] = (v.x[i][j][k_end] - v.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + v.x[i][j][k_beg];
                        w.x[i][j][k] = (w.x[i][j][k_end] - w.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + w.x[i][j][k_beg];
                    }
                }
            }
        }


    // Atlantic
    // Connection between Guyana and Gulf of Mexico current
    // North America      East coasts
    // Florida Strait at the southern tip of Florida (from j=63 until j=67 compares to 23°N until 27°N,
    // parallel to the equator                         from k=264 until k=283 compares to 77°W until 96°W)
        j_beg = 63;
        j_end = 68;
        k_beg = 264;
        k_end = 284;

        #pragma omp parallel for collapse(3) private(j_z, j_n)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    j_z = j - j_beg;
                    j_n = j_end - j_beg;

                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = (v.x[i][j_end][k] - v.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + v.x[i][j_beg][k];
                        w.x[i][j][k] = (w.x[i][j_end][k] - w.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + w.x[i][j_beg][k];
                    }
                }
            }
        }


    // Atlantic
    // North America      East coasts
    // Gulf Stream in the north-east of North America (from j=40 until j=70 compares to 20°N until 50°N,
    // south of Newfoundland                    from k=278 until k=308 compares to 82°W until 52°W)
        j_beg = 40;
        j_end = 71;
        k_beg = 278;
        k_end = 309;
        k_step = 20;

        v_grad = + 0.001;

        k_a = k_b = 0;

        flip = 0;

        for (int k = std::min(k_end, m.km - 1); k > k_beg; k--) {
            if (h.x[i_max][j_beg][k] == 1.) {
                k_a = k;
                flip = 1;
            }
            if (flip == 1) break;
        }

        flip = 0;

        for (int j = j_beg+1; j < j_end; j++) {
            for (int k = std::min(k_end, m.km - 1); k > k_beg; k--) {
                if (h.x[i_max][j][k] == 1.) {
                    k_b = k;
                    k_grad = k_a - k_b;
                    if (k_grad >= 2) k_grad = 1;
                    if (k_grad <= 0) k_grad = 1;
                    flip = 1;
                }
            if (flip == 1) break;
            }

            for (int k = k_b; k < (k_b + k_step) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg; i < m.im; i++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = - v_grad * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = + v_grad * static_cast<double>(k_grad) * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }
            k_a = k_b;
            flip = 0;

            for (int k = k_b + 4; k < (k_b + k_step) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg; i < m.im; i++) {
                    v.x[i][j][k] = (v.x[i][j][std::max(0, std::min(k_b + k_step +1, m.km - 1))] - v.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))]) * static_cast<double>(k - (k_b + 4)) / static_cast<double>((k_b + k_step -1) - (k_b + 4)) + v.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))];
                    w.x[i][j][k] = (w.x[i][j][std::max(0, std::min(k_b + k_step +1, m.km - 1))] - w.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))]) * static_cast<double>(k - (k_b + 4)) / static_cast<double>((k_b + k_step -1) - (k_b + 4)) + w.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))];
                }
            }

        }


    // Atlantic
    // Connection between Guyana and Gulf Stream
    // North America      East coasts
    // Gulf Stream in the north-east of North America (from j=44 until j=70 compares to 20°N until 46°N,
    // south of Newfoundland                    from k=280 until k=308 compares to 80°W until 52°W)
        j_beg = 44;
        j_end = 71;
        k_beg = 280;
        k_end = 309;

        #pragma omp parallel for collapse(3) private(j_z, j_n)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    j_z = j - j_beg;
                    j_n = j_end - j_beg;

                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = (v.x[i][j_end][k] - v.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + v.x[i][j_beg][k];
                        w.x[i][j][k] = (w.x[i][j_end][k] - w.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + w.x[i][j_beg][k];
                    }
                }
            }
        }


    // Atlantic
    // Gulf of St. Lawrence        velocity reduction       west coast
    // Gulf of St. Lawrence (from j=38 until j=45 compares to 45°N until 52°N,
    //                                 from k=295 until k=305 compares to 65°W until 55°W)
        #pragma omp parallel for collapse(3)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = 38; j < 46; j++) {
                for (int k = 295; k < 306; k++) {
                    if (h.x[i][j][k] == 0.) {
    //                  u.x[i][j][k] = - IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        v.x[i][j][k] = + IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = - IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }
        }


    // Atlantic
    // Thermohalin Conveyor Belt
    // below Newfoundland and the St. Lawrence Current as connecting piece
    // Gulf Stream on the Atlantic from south-west in north-east direction (from j=43 until j=51 compares to 39°N until 47°N,
    //                                                                  from k=295 until k=320 compares to 40°W until 65°W)
        j_beg = 43;
        j_end = 52;
        k_beg = 295;
        k_end = 321;
        j_run = 0;
        k_run = 0;
        j_step = 1;
        k_step = 50;

        while ((j_end - j_run) >= j_beg && (k_beg + k_run) <= k_end) {
            for (int j = std::min(j_end - j_run, m.jm - 1); j > (j_end - j_step - j_run); j--) {
                for (int k = (k_beg + k_run); k < (k_beg + k_step + k_run) && k < m.km && k >= 0; k++) {
                    for (int i = i_beg; i < m.im; i++) {
                        if (h.x[i][j][k] == 0.) {
                            v.x[i][j][k] = - IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                            w.x[i][j][k] = + IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        }
                    }
                }
            k_run++;
            }
        j_run++;
        }


    // Atlantic
    // Thermohalin Conveyor Belt
    // Connection below Newfoundland and the St. Lawrence Current and the North Atlantic
    // Gulf Stream on the Atlantic from south-west in north-east direction (from j=48 until j=54 compares to 36°N until 42°N,
    //                                                                                                           from k=295 until k=340 compares to 20°W until 65°W
        j_beg = 48;
        j_end = 55;
        k_beg = 295;
        k_end = 341;

        #pragma omp parallel for collapse(3) private(j_z, j_n)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    j_z = j - j_beg;
                    j_n = j_end - j_beg;

                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = (v.x[i][j_end][k] - v.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + v.x[i][j_beg][k];
                        w.x[i][j][k] = (w.x[i][j_end][k] - w.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + w.x[i][j_beg][k];
                    }
                }
            }
        }


    // Atlantic
    // Thermohalin Conveyor Belt
    // diagonal from the Labrador tip until northern Norway
    // Gulf Stream on the Atlantic from south-west in north-east direction (from j=37 until j=44 compares to 63°N until 46°N,
    // extension                                                                                     from k=304 until k=m.km compares to 56°W until 0°W)
        j_beg = 37;
        j_end = 45;
        k_beg = 304;
        k_end = m.km;
        j_run = 0;
        k_run = 0;
        j_step = 1;
        k_step = 50;

        while ((j_end - j_run) >= j_beg && (k_beg + k_run) <= k_end) {
            for (int j = std::min(j_end - j_run, m.jm - 1); j > (j_end - j_step - j_run); j--) {
                for (int k = (k_beg + k_run); k < (k_beg + k_step + k_run) && k < m.km && k >= 0; k++) {
                    for (int i = i_beg; i < m.im; i++) {
                        if (h.x[i][j][k] == 0.) {
                            v.x[i][j][k] = - IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                            w.x[i][j][k] = + IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        }
                    }
                }
            k_run++;
            }
        j_run++;
        }


    // Atlantic
    // Transition from Gulf Stream to North Atlantic
    // Greenland Current in the east                                      (from j=15 until j=27 compares to 75°N until 63°N,
    // parallel to the equator                                                from k=350 until k=m.km compares to 10°W until 0°)
        j_beg = 15;
        j_end = 28;
        k_beg = 350;
        k_end = m.km;

            #pragma omp parallel for collapse(3) private(j_z, j_n)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    j_z = j - j_beg;
                    j_n = j_end - j_beg;

                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = (v.x[i][j_end][k] - v.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + v.x[i][j_beg][k];
                        w.x[i][j][k] = (w.x[i][j_end][k] - w.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + w.x[i][j_beg][k];
                    }
                }
            }
        }


    // Atlantic
    // Connection between Atlantic and Norwegian Sea, Gulf Stream toward Norway
    // North Atlantic
    // diagonal from the Labrador tip until northern Norway
    // Gulf Stream on the Atlantic from south-west in north-east direction (from j=27 until j=44 compares to 63°N until 46°N,
    // extension                                                                                    from k=304 until k=m.km compares to 56°W until 0°W)
        j_beg = 27;
        j_end = 45;
        k_beg = 304;
        k_end = m.km;

        #pragma omp parallel for collapse(3) private(j_z, j_n)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    j_z = j - j_beg;
                    j_n = j_end - j_beg;

                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = (v.x[i][j_end][k] - v.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + v.x[i][j_beg][k];
                        w.x[i][j][k] = (w.x[i][j_end][k] - w.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + w.x[i][j_beg][k];
                    }
                }
            }
        }


    // Atlantic
    // Thermohalin Conveyor Belt
    // diagonal as extension beyond 0° until northern Norwa
    // Gulf Stream on the Atlantic from south-west in north-east direction (from j=18 until j=30 compares to 72°N until 60°N,
    //                                                                                                           from k=0 until k=20 compares to 0°W until 20°O)
        j_beg = 18;
        j_end = 31;
        k_beg = 0;
        k_end = 21;
        k_step = 10;

        k_exp = 10;
        k_w = 4;

        v_grad = 0.0008;

        k_a = k_b = 0;

        flip = 0;

        for (int k = k_beg; k < k_end; k++) {
            if (h.x[i_max][j_beg][k] == 1.) {
                k_a = k;
                flip = 1;
            }
            if (flip == 1) break;
        }

        flip = 0;

        for (int j = j_beg+1; j < j_end; j++) {
            for (int k = k_beg; k < k_end; k++) {
                if (h.x[i_max][j][k] == 1.) {
                    k_b = k;
                    k_grad = k_a - k_b;
                    if (k_grad >= 2) k_grad = 1;
                    if (k_grad <= 0) k_grad = 1;
                    flip = 1;
                }
            if (flip == 1) break;
            }

            for (int k = k_b; k > (k_b - k_step) && k >= 0; k--) {
                for (int i = i_beg; i < m.im; i++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = - v_grad * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = + v_grad * static_cast<double>(k_grad) * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }
            k_a = k_b;
            flip = 0;

            for (int k = k_b - k_step - k_exp; k < (k_b - k_w) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg; i < m.im; i++) {
                    v.x[i][j][k] = (v.x[i][j][std::max(k_b - k_w, 0)] - v.x[i][j][std::max(0, std::min(k_b - k_step - k_exp, m.km - 1))]) * static_cast<double>(k - (k_b - k_step - k_exp)) / static_cast<double>((k_b - k_w) - (k_b - k_step - k_exp)) + v.x[i][j][std::max(0, std::min(k_b - k_step - k_exp, m.km - 1))];
                    w.x[i][j][k] = (w.x[i][j][std::max(k_b - k_w, 0)] - w.x[i][j][std::max(0, std::min(k_b - k_step - k_exp, m.km - 1))]) * static_cast<double>(k - (k_b - k_step - k_exp)) / static_cast<double>((k_b - k_w) - (k_b - k_step - k_exp)) + w.x[i][j][std::max(0, std::min(k_b - k_step - k_exp, m.km - 1))];
                }
            }
        }


    // Atlantic
    // Greenland     West coasts
    // Greenland Current in the west until into Baffin Bay (from j=15 until j=31 compares to 75°N until 59°N,
    //                                                   from k=300 until k=315 compares to 60°W until 45°W)
        j_beg  = 15;
        j_end  = 32;
        k_beg  = 300;
        k_end   = 316;
        k_step = 4;
        k_w    = 2;
        k_exp  = 4;

        v_grad = + 0.0010;

        k_a = k_b = 0;

        flip = 0;

        for (int k = k_beg; k < k_end; k++) {
            if (h.x[i_max][j_beg][k] == 1.) {
                k_a = k;
                flip = 1;
            }
            if (flip == 1) break;
        }

        flip = 0;

        for (int j = j_beg+1; j < j_end; j++) {
            for (int k = k_beg; k < k_end; k++) {
                if (h.x[i_max][j][k] == 1.) {
                    k_b = k;
                    k_grad = k_a - k_b;
                    if (k_grad >= 2) k_grad = 1;
                    if (k_grad <= 0) k_grad = 1;
                    flip = 1;
                }
            if (flip == 1) break;
            }

            for (int k = k_b; k > (k_b - k_step) && k >= 0; k--) {
                for (int i = i_beg; i < m.im; i++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = - v_grad * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = - v_grad * static_cast<double>(k_grad) * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }
            k_a = k_b;
            flip = 0;

            for (int k = k_b - k_step - k_exp; k < (k_b - k_w) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg; i < m.im; i++) {
                    v.x[i][j][k] = (v.x[i][j][std::max(k_b - k_w, 0)] - v.x[i][j][std::max(0, std::min(k_b - k_step - k_exp, m.km - 1))]) * static_cast<double>(k - (k_b - k_step - k_exp)) / static_cast<double>((k_b - k_w) - (k_b - k_step - k_exp)) + v.x[i][j][std::max(0, std::min(k_b - k_step - k_exp, m.km - 1))];
                    w.x[i][j][k] = (w.x[i][j][std::max(k_b - k_w, 0)] - w.x[i][j][std::max(0, std::min(k_b - k_step - k_exp, m.km - 1))]) * static_cast<double>(k - (k_b - k_step - k_exp)) / static_cast<double>((k_b - k_w) - (k_b - k_step - k_exp)) + w.x[i][j][std::max(0, std::min(k_b - k_step - k_exp, m.km - 1))];
                }
            }
        }


    // Atlantic
    // Southern tip of Greenland      West coasts
    // Greenland Current in the west until into Baffin Bay (from j=15 until j=31 compares to 75°N until 59°N,
    // parallel to the equator                                                from k=310 until k=315 compares to 50°W until 45°W)

        #pragma omp parallel for collapse(3)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = 15; j < 32; j++) {
                for (int k = 310; k < 316; k++) {
                    if (h.x[i][j][k] == 0.) {
                        {
                            v.x[i][j][k] = + 0.0;
                            w.x[i][j][k] = - IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        }
                    }
                }
            }
        }


    // Atlantic
    // Transition at the southern tip of Greenland to the Atlantic     West coasts
    // Greenland Current in the west until into Baffin Bay (from j=15 until j=31 compares to 75°N until 59°N,
    // parallel to the equator                                                from k=300 until k=315 compares to 60°W until 45°W)
        j_beg = 15;
        j_end = 32;
        k_beg = 300;
        k_end = 316;

        #pragma omp parallel for collapse(3) private(j_z, j_n)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    j_z = j - j_beg;
                    j_n = j_end - j_beg;

                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = (v.x[i][j_end][k] - v.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + v.x[i][j_beg][k];
                        w.x[i][j][k] = (w.x[i][j_end][k] - w.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + w.x[i][j_beg][k];
                    }
                }
            }
        }


    // Atlantic
    // Greenland     East coasts
    // Thermohalin Conveyor Belt
    // East Greenland Current in the north-east of Greenland (from j=10 until j=31 compares to 80°N until 59°N,
    //                                                                                    from k=315 until k=342 compares to 45°W until 18°W)
        j_beg = 10;
        j_end = 32;
        k_beg = 315;
        k_end = 343;
        k_step = 8;
        k_w = 4;
        k_exp = 4;

        v_grad = + 0.0008;

        k_a = k_b = 0;

        flip = 0;

        for (int k = std::min(k_end, m.km - 1); k > k_beg; k--) {
            if (h.x[i_max][j_beg][k] == 1.) {
                k_a = k;
                flip = 1;
            }
            if (flip == 1) break;
        }

        flip = 0;

        for (int j = j_beg+1; j < j_end; j++) {
            for (int k = std::min(k_end, m.km - 1); k > k_beg; k--) {
                if (h.x[i_max][j][k] == 1.) {
                    k_b = k;
                    k_grad = k_a - k_b;
                    if (k_grad >= 2) k_grad = 1;
                    if (k_grad <= 0) k_grad = 1;
                    flip = 1;
                }
            if (flip == 1) break;
            }

            for (int k = k_b; k < (k_b + k_step) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg-3; i < m.im; i++) {
                    if (h.x[i][j][k] == 0.) {
//                        u.x[i][j][k] = - IC_water;
                        v.x[i][j][k] = + v_grad * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = - v_grad * static_cast<double>(k_grad) * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }
            k_a = k_b;
            flip = 0;

            for (int k = k_b + k_w; k < (k_b + k_step) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg; i < m.im; i++) {
                    v.x[i][j][k] = (v.x[i][j][std::max(0, std::min(k_b + k_step, m.km - 1))] - v.x[i][j][std::max(0, std::min(k_b + k_w, m.km - 1))]) * static_cast<double>(k - (k_b + k_w)) / static_cast<double>((k_b + k_step) - (k_b + k_w)) + v.x[i][j][std::max(0, std::min(k_b + k_w, m.km - 1))];
                    w.x[i][j][k] = (w.x[i][j][std::max(0, std::min(k_b + k_step, m.km - 1))] - w.x[i][j][std::max(0, std::min(k_b + k_w, m.km - 1))]) * static_cast<double>(k - (k_b + k_w)) / static_cast<double>((k_b + k_step) - (k_b + k_w)) + w.x[i][j][std::max(0, std::min(k_b + k_w, m.km - 1))];
                }
            }

        }


    // Atlantic
    // Iceland       South coast
    // Thermohalin Conveyor Belt
    // South Iceland Current (from j=26 until j=29 compares to 64°N until 61°N,
    //                                 from k=326 until k=347 compares to 34°W until 13°W)
        #pragma omp parallel for collapse(3)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = 26; j < 30; j++) {
                for (int k = 326; k < 348; k++) {
                    if (h.x[i][j][k] == 0.) {
                        u.x[i][j][k] = - IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        v.x[i][j][k] = + IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = - IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }
        }


    // Atlantic
    // southern transition at Iceland        South coast
    // Thermohalin Conveyor Belt
    // South Iceland Current (from j=24 until j=29 compares to 64°N until 61°N,
    //                                 from k=326 until k=347 compares to 34°W until 13°W)
        j_beg = 24;
        j_end = 30;
        k_beg = 326;
        k_end = 348;

        #pragma omp parallel for collapse(3) private(j_z, j_n)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    j_z = j - j_beg;
                    j_n = j_end - j_beg;

                    if (h.x[i][j][k] == 0.) {
                        u.x[i][j][k] = (u.x[i][j_end][k] - u.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + u.x[i][j_beg][k];
                        v.x[i][j][k] = (v.x[i][j_end][k] - v.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + v.x[i][j_beg][k];
                        w.x[i][j][k] = (w.x[i][j_end][k] - w.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + w.x[i][j_beg][k];
                    }
                }
            }
        }


    // Atlantic
    // western transition at Iceland       South coast
    // Thermohalin Conveyor Belt
    // South Iceland Current (from j=24 until j=29 compares to 64°N until 61°N,
    //                                 from k=326 until k=347 compares to 34°W until 13°W)
        j_beg = 24;
        j_end = 30;
        k_beg = 326;
        k_end = 348;

        #pragma omp parallel for collapse(3) private(k_z, k_n)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    k_z = k - k_beg;
                    k_n = k_end - k_beg;

                    if (h.x[i][j][k] == 0.) {
                        u.x[i][j][k] = (u.x[i][j][k_end] - u.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + u.x[i][j][k_beg];
                        v.x[i][j][k] = (v.x[i][j][k_end] - v.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + v.x[i][j][k_beg];
                        w.x[i][j][k] = (w.x[i][j][k_end] - w.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + w.x[i][j][k_beg];
                    }
                }
            }
        }


    // Atlantic
    // eastern transition at Iceland        South coast
    // Thermohalin Conveyor Belt
    // South Iceland Current (from j=24 until j=29 compares to 64°N until 61°N,
    //                                from k=326 until k=347 compares to 34°W until 13°W)
        j_beg = 24;
        j_end = 30;
        k_beg = 326;
        k_end = 348;

        #pragma omp parallel for collapse(3) private(k_z, k_n)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    k_z = k - k_beg;
                    k_n = k_end - k_beg;

                    if (h.x[i][j][k] == 0.) {
                        u.x[i][j][k] = (u.x[i][j][k_end] - u.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + u.x[i][j][k_beg];
                        v.x[i][j][k] = (v.x[i][j][k_end] - v.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + v.x[i][j][k_beg];
                        w.x[i][j][k] = (w.x[i][j][k_end] - w.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + w.x[i][j][k_beg];
                    }
                }
            }
        }


    // Atlantic
    // Thermohalin Conveyor Belt
    // Preparation for the outflow of the Gulf Stream south of Iceland
    // Gulf Stream on the Atlantic toward the south (from j=32 until j=40 compares to 58°N until 50°N,
    //                                                                     from k=317 until k=327 compares to 43°W until 33°W)
        #pragma omp parallel for collapse(3)
        for (int i = i_bottom; i < m.im-1; i++) {
            for (int j = 32; j < 41; j++) {
                for (int k = 317; k < 328; k++) {
                    if (h.x[i][j][k] == 0.) {
                        c.x[i][j][k] = c.x[m.im-1][j][k];
                        u.x[i][j][k] = - IC_water;
//                        v.x[i][j][k] = - IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
//                        w.x[i][j][k] = - IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        v.x[i][j][k] = + IC_water;
                        w.x[i][j][k] = + IC_water;
                    }
                }
            }
        }


    // Atlantic
    // South America       East coasts
    // Brazil Current in the south-east of South America (from j=98 until j=127 compares to 8°S until 37°S,
    // surface flow                                       from k=305 until k=325 compares to 55°W until 35°W)
        j_beg = 98;
        j_end = 128;
        k_beg = 305;
        k_end = 326;
        k_step = 20;

        v_grad = + 0.001;

        k_a = k_b = 0;

        flip = 0;

        for (int k = std::min(k_end, m.km - 1); k > k_beg; k--) {
            if (h.x[i_max][j_beg][k] == 1.) {
                k_a = k;
                flip = 1;
            }
            if (flip == 1) break;
        }

        flip = 0;

        for (int j = j_beg+1; j < j_end; j++) {
            for (int k = std::min(k_end, m.km - 1); k > k_beg; k--) {
                if (h.x[i_max][j][k] == 1.) {
                    k_b = k;
                    k_grad = k_a - k_b;
                    if (k_grad >= 2) k_grad = 1;
                    if (k_grad <= 0) k_grad = 1;
                    flip = 1;
                }
            if (flip == 1) break;
            }

            for (int k = k_b; k < (k_b + k_step) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg; i < m.im; i++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = + v_grad * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = - v_grad * static_cast<double>(k_grad) * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }

            for (int k = k_b + 4; k < (k_b + k_step) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg; i < m.im; i++) {
                    v.x[i][j][k] = (v.x[i][j][std::max(0, std::min(k_b + k_step +1, m.km - 1))] - v.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))]) * static_cast<double>(k - (k_b + 4)) / static_cast<double>((k_b + k_step -1) - (k_b + 4)) + v.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))];
                    w.x[i][j][k] = (w.x[i][j][std::max(0, std::min(k_b + k_step +1, m.km - 1))] - w.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))]) * static_cast<double>(k - (k_b + 4)) / static_cast<double>((k_b + k_step -1) - (k_b + 4)) + w.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))];
                }
            }
            k_a = k_b;
            flip = 0;
        }


    // Atlantic
    // South America       East coasts
    // Falkland Current in the south-east of South America (from j=127 until j=146 compares to 37°S until 56°S,
    //                                                                               from k=290 until k=330 compares to 70°W until 30°W)

        j_beg = 127;
        j_end = 147;
        k_beg = 290;
        k_end = 331;
        k_step = 20;

        v_grad = + 0.001;

        k_a = k_b = 0;

        flip = 0;

        for (int k = std::min(k_end, m.km - 1); k > k_beg; k--) {
            if (h.x[i_max][j_beg][k] == 1.) {
                k_a = k;
                flip = 1;
            }
            if (flip == 1) break;
        }

        flip = 0;

        for (int j = j_beg+1; j < j_end; j++) {
            for (int k = std::min(k_end, m.km - 1); k > k_beg; k--) {
                if (h.x[i_max][j][k] == 1.) {
                    k_b = k;
                    k_grad = k_a - k_b;
                    if (k_grad >= 2) k_grad = 1;
                    if (k_grad <= 0) k_grad = 1;
                    flip = 1;
                }
            if (flip == 1) break;
            }

            for (int k = k_b; k < (k_b + k_step) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg; i < m.im; i++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = - v_grad * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = + v_grad * static_cast<double>(k_grad) * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }
            k_a = k_b;
            flip = 0;

            for (int k = k_b + 4; k < (k_b + k_step) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg; i < m.im; i++) {
                    v.x[i][j][k] = (v.x[i][j][std::max(0, std::min(k_b + k_step +1, m.km - 1))] - v.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))]) * static_cast<double>(k - (k_b + 4)) / static_cast<double>((k_b + k_step -1) - (k_b + 4)) + v.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))];
                    w.x[i][j][k] = (w.x[i][j][std::max(0, std::min(k_b + k_step +1, m.km - 1))] - w.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))]) * static_cast<double>(k - (k_b + 4)) / static_cast<double>((k_b + k_step -1) - (k_b + 4)) + w.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))];
                }
            }
        }


    // Connection between Brazil and Falkland current
    // south of the Brazil Current in the south-east of South America (from j=98 until j=127 compares to 8°S until 37°S,
    // surface flow with Falkland current                            from k=300 until k=325 compares to 60°W until 35°W)
        j_beg = 198;
        j_end = 128;
        k_beg = 300;
        k_end = 326;

        #pragma omp parallel for collapse(3) private(j_z, j_n)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = j_beg; j < j_end + 1; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    j_z = j - j_beg;
                    j_n = j_end - j_beg;

                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = (v.x[i][j_end][k] - v.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + v.x[i][j_beg][k];
                        w.x[i][j][k] = (w.x[i][j_end][k] - w.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + w.x[i][j_beg][k];
                    }
                }
            }
        }


    // Atlantic
    // North Africa       West coasts
    // Canary Current in the west of North Africa (from j=50 until j=77 compares to 13°N until 40°N,
    //                                                                        from k=340 until k=352 compares to 20°W until 8°W)
        j_beg = 50;
        j_end = 78;
        k_beg = 340;
        k_end = 353;

        v_grad = + 0.001;
        k_step = 8;
        k_w = 4;
        k_exp = 4;

        k_a = k_b = 0;

        flip = 0;

        for (int k = k_beg; k < k_end; k++) {
            if (h.x[i_max][j_beg][k] == 1.) {
                k_a = k;
                flip = 1;
            }
            if (flip == 1) break;
        }

        flip = 0;

        for (int j = j_beg+1; j < j_end; j++) {
            for (int k = k_beg; k < k_end; k++) {
                if (h.x[i_max][j][k] == 1.) {
                    k_b = k;
                    k_grad = k_a - k_b;
                    if (k_grad >= 2) k_grad = 1;
                    if (k_grad <= 0) k_grad = 1;
                    flip = 1;
                }
            if (flip == 1) break;
            }

            for (int k = k_b; k > (k_b - k_step) && k >= 0; k--) {
                for (int i = i_beg; i < m.im; i++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = + v_grad * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = - v_grad * static_cast<double>(k_grad) * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }
            k_a = k_b;
            flip = 0;

            for (int k = k_b - k_step - k_exp; k < (k_b - k_w) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg; i < m.im; i++) {
                    v.x[i][j][k] = (v.x[i][j][std::max(k_b - k_w, 0)] - v.x[i][j][std::max(0, std::min(k_b - k_step - k_exp, m.km - 1))]) * static_cast<double>(k - (k_b - k_step - k_exp)) / static_cast<double>((k_b - k_w) - (k_b - k_step - k_exp)) + v.x[i][j][std::max(0, std::min(k_b - k_step - k_exp, m.km - 1))];
                    w.x[i][j][k] = (w.x[i][j][std::max(k_b - k_w, 0)] - w.x[i][j][std::max(0, std::min(k_b - k_step - k_exp, m.km - 1))]) * static_cast<double>(k - (k_b - k_step - k_exp)) / static_cast<double>((k_b - k_w) - (k_b - k_step - k_exp)) + w.x[i][j][std::max(0, std::min(k_b - k_step - k_exp, m.km - 1))];
                }
            }
        }


    // Atlantic
    // North Africa       West coasts
    // Guinea Current in the west of Africa (from j=77 until j=89 compares to 13°N until 1°N,
    //                                                              from k=340 until k=352 compares to 20°W until 8°W)
        j_beg = 77;
        j_end = 90;
        k_beg = 340;
        k_end = 353;

        v_grad = + 0.001;
        k_step = 8;
        k_w = 4;
        k_exp = 4;

        k_a = k_b = 0;

        flip = 0;

        for (int k = k_beg; k < k_end; k++) {
            if (h.x[i_max][j_beg][k] == 1.) {
                k_a = k;
                flip = 1;
            }
            if (flip == 1) break;
        }

        flip = 0;

        for (int j = j_beg+1; j < j_end; j++) {
            for (int k = k_beg; k < k_end; k++) {
                if (h.x[i_max][j][k] == 1.) {
                    k_b = k;
                    k_grad = k_a - k_b;
                    if (k_grad >= 2) k_grad = 1;
                    if (k_grad <= 0) k_grad = 1;
                    flip = 1;
                }
            if (flip == 1) break;
            }

            for (int k = k_b; k > (k_b - k_step) && k >= 0; k--) {
                for (int i = i_beg; i < m.im; i++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = + v_grad * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = + v_grad * static_cast<double>(k_grad) * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }
            k_a = k_b;
            flip = 0;

            for (int k = k_b - k_step - k_exp; k < (k_b - k_w) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg; i < m.im; i++) {
                    v.x[i][j][k] = (v.x[i][j][std::max(k_b - k_w, 0)] - v.x[i][j][std::max(0, std::min(k_b - k_step - k_exp, m.km - 1))]) * static_cast<double>(k - (k_b - k_step - k_exp)) / static_cast<double>((k_b - k_w) - (k_b - k_step - k_exp)) + v.x[i][j][std::max(0, std::min(k_b - k_step - k_exp, m.km - 1))];
                    w.x[i][j][k] = (w.x[i][j][std::max(k_b - k_w, 0)] - w.x[i][j][std::max(0, std::min(k_b - k_step - k_exp, m.km - 1))]) * static_cast<double>(k - (k_b - k_step - k_exp)) / static_cast<double>((k_b - k_w) - (k_b - k_step - k_exp)) + w.x[i][j][std::max(0, std::min(k_b - k_step - k_exp, m.km - 1))];
                }
            }
        }


    // Atlantic
    // North Africa       West coasts
    // Connection between Canary and Guinea current
    // Canary Current in the west of North Africa (from j=50 until j=77 compares to 13°N until 40°N,
    //                                                                        from k=340 until k=352 compares to 20°W until 8°W)
        j_beg = 50;
        j_end = 78;
        k_beg = 340;
        k_end = 353;

        #pragma omp parallel for collapse(3) private(j_z, j_n)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = j_beg; j < j_end + 1; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    j_z = j - j_beg;
                    j_n = j_end - j_beg;

                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = (v.x[i][j_end][k] - v.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + v.x[i][j_beg][k];
                        w.x[i][j][k] = (w.x[i][j_end][k] - w.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + w.x[i][j_beg][k];
                    }
                }
            }
        }


    // Atlantic
    // North Africa       West coasts
    // Connection between Canary current and North Atlantic
    // Canary Current in the west of North Africa (from j=50 until j=77 compares to 13°N until 40°N,
    //                                                                        from k=340 until k=352 compares to 20°W until 8°W)
        j_beg = 50;
        j_end = 78;
        k_beg = 340;
        k_end = 353;

        #pragma omp parallel for collapse(3) private(j_z, j_n)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = j_beg; j < j_end + 1; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    j_z = j - j_beg;
                    j_n = j_end - j_beg;

                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = (v.x[i][j_end][k] - v.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + v.x[i][j_beg][k];
                        w.x[i][j][k] = (w.x[i][j_end][k] - w.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + w.x[i][j_beg][k];
                    }
                }
            }
        }


    // Atlantic
    // North Africa       West coasts
    // Connection between Guinea current and South Atlantic
    // Guinea Current in the west of Africa (from j=77 until j=89 compares to 13°N until 1°N,
    //                                                              from k=340 until k=352 compares to 20°W until 8°W)
        j_beg = 77;
        j_end = 90;
        k_beg = 340;
        k_end = 353;

        #pragma omp parallel for collapse(3) private(j_z, j_n)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = j_beg; j < j_end + 1; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    j_z = j - j_beg;
                    j_n = j_end - j_beg;

                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = (v.x[i][j_end][k] - v.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + v.x[i][j_beg][k];
                        w.x[i][j][k] = (w.x[i][j_end][k] - w.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + w.x[i][j_beg][k];
                    }
                }
            }
        }
    }
/*
* 
*/
    void IC_Indischer_Ozean (Array &h, Array &u, Array &v, Array &w) {
    // Currents along the coasts
    // Closing the polar, subpolar and subtropical atmospheric circulation systems
    // Thermohalin Conveyor Belt
    // Indian Ocean
    // diagonal from Indonesia/Australia until Kap Agulhas in South Africa
    // Westward flowing current as surface flow (from j=104 until j=135 compares to 14°S until 45°S,
    //                                                                  from k=10 until k=120 compares to 10°O until 120°O)

        j_beg = 104;
        j_end = 136;
        k_beg = 10;
        k_end = 121;
        j_run = 0;
        k_run = 0;
        j_step = 1;
        k_step = 20;


        while ((j_end - j_run) >= j_beg && (k_beg + k_run) <= k_end) {
            for (int j = std::min(j_end - j_run, m.jm - 1); j > (j_end - j_step - j_run); j--) {
                for (int k = (k_beg + k_run); k < (k_beg + k_step + k_run) && k < m.km && k >= 0; k++) {
                    for (int i = i_beg; i < m.im; i++) {
                        if (h.x[i][j][k] == 0.) {
                            v.x[i][j][k] = IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                            w.x[i][j][k] = - IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        }
                    }
                }
            k_run = k_run + 3;
            }
        j_run++;
        }


    // west of the constant velocity strip diagonal from Indonesia/Australia until Kap Agulhas in South Africa
    // Westward flowing current as surface flow (from j=104 until j=135 compares to 14°S until 45°S,
    //                                                                                      from k=10 until k=120 compares to 10°O until 120°O)
        j_beg = 104;
        j_end = 136;
        k_beg = 10;
        k_end = 121;
        j_run = 0;
        k_run = 0;
        j_step = 1;
        k_step = 20;
        k_exp = 30;

        while ((j_end - j_run) >= j_beg && (k_beg + k_run + k_step - k_exp) <= (k_end - k_exp)) {
            for (int j = std::min(j_end - j_run, m.jm - 1); j > (j_end - j_step - j_run); j--) {
                for (int k = (k_beg + k_run + k_step - k_exp); k < (k_beg + k_run + k_step) && k < m.km && k >= 0; k++) {
                    k_z = k - (k_beg + k_run + k_step - k_exp);
                    k_n = (k_beg + k_run + k_step) - (k_beg + k_run + k_step - k_exp) ;

                    for (int i = i_beg; i < m.im; i++) {
                        if (h.x[i][j][k] == 0.) {
                            v.x[i][j][k] = (v.x[i][j][std::max(0, std::min(k_beg + k_run + k_step - 1, m.km - 1))] - v.x[i][j][std::max(0, std::min(k_beg + k_run + k_step - k_exp, m.km - 1))]) * static_cast<double>(k_z) / static_cast<double>(k_n) + v.x[i][j][std::max(0, std::min(k_beg + k_run + k_step - k_exp, m.km - 1))];
                            w.x[i][j][k] = (w.x[i][j][std::max(0, std::min(k_beg + k_run + k_step - 1, m.km - 1))] - w.x[i][j][std::max(0, std::min(k_beg + k_run + k_step - k_exp, m.km - 1))]) * static_cast<double>(k_z) / static_cast<double>(k_n) + w.x[i][j][std::max(0, std::min(k_beg + k_run + k_step - k_exp, m.km - 1))];
                        }
                    }
                }
                k_run = k_run + 3;
            }
        j_run++;
        }


    // east of the constant velocity strip diagonal from Indonesia/Australia until Kap Agulhas in South Africa
    // Westward flowing current as surface flow (from j=104 until j=135 compares to 14°S until 45°S,
    //                                                                                      from k=10 until k=120 compares to 10°O until 120°O)
        j_beg = 104;
        j_end = 136;
        k_beg = 10;
        k_end = 121;
        j_run = 0;
        k_run = 0;
        j_step = 1;
        k_step = 20;
        k_exp = 30;

        while ((j_end - j_run) >= j_beg && (k_beg + k_run + k_step + k_exp) <= (k_end + k_exp)) {
            for (int j = std::min(j_end - j_run, m.jm - 1); j > (j_end - j_step - j_run); j--) {
                for (int k = (k_beg + k_run + k_step); k < (k_beg + k_run + k_step + k_exp) && k < m.km && k >= 0; k++) {
                    k_z = k - (k_beg + k_run + k_step);
                    k_n = (k_beg + k_run + k_step + k_exp) - (k_beg + k_run + k_step) ;

                    for (int i = i_beg; i < m.im; i++) {
                        if (h.x[i][j][k] == 0.) {
                            v.x[i][j][k] = (v.x[i][j][std::min(k_beg + k_run + k_step + k_exp, m.km - 1)] - v.x[i][j][std::max(0, std::min(k_beg + k_run + k_step  - 1, m.km - 1))]) * static_cast<double>(k_z) / static_cast<double>(k_n) + v.x[i][j][std::max(0, std::min(k_beg + k_run + k_step - 1, m.km - 1))];
                            w.x[i][j][k] = (w.x[i][j][std::min(k_beg + k_run + k_step + k_exp, m.km - 1)] - w.x[i][j][std::max(0, std::min(k_beg + k_run + k_step  - 1, m.km - 1))]) * static_cast<double>(k_z) / static_cast<double>(k_n) + w.x[i][j][std::max(0, std::min(k_beg + k_run + k_step - 1, m.km - 1))];
                        }
                    }
                }
                k_run = k_run + 3;
            }
        j_run++;
        }


    // Coastal currents
    // Indian Ocean
    // South Africa/Madagascar     East coasts
    // Agulhas Current in the south-east of South Africa (from j=100 until j=126 compares to 10°S until 36°S,
    //                                                                          from k=20 until k=42 compares to 20°O until 42°O)
        j_beg = 100;
        j_end = 127;
        k_beg = 20;
        k_end = 43;
        k_step = 14;

        v_grad = + 0.001;

        k_a = k_b = 0;

        flip = 0;

        for (int k = std::min(k_end, m.km - 1); k > k_beg; k--) {
            if (h.x[i_max][j_beg][k] == 1.) {
                k_a = k;
                flip = 1;
            }
            if (flip == 1) break;
        }

        flip = 0;

        for (int j = j_beg+1; j < j_end; j++) {
            for (int k = std::min(k_end, m.km - 1); k > k_beg; k--) {
                if (h.x[i_max][j][k] == 1.) {
                    k_b = k;
                    k_grad = k_a - k_b;
                    if (k_grad >= 2) k_grad = 1;
                    if (k_grad <= 0) k_grad = 1;
                    flip = 1;
                }
            if (flip == 1) break;
            }

            for (int k = k_b; k < (k_b + k_step) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg; i < m.im; i++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = + v_grad * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = - v_grad * static_cast<double>(k_grad) * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }
            k_a = k_b;
            flip = 0;

            for (int k = k_b + 4; k < (k_b + k_step) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg; i < m.im; i++) {
                    v.x[i][j][k] = (v.x[i][j][std::max(0, std::min(k_b + k_step +1, m.km - 1))] - v.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))]) * static_cast<double>(k - (k_b + 4)) / static_cast<double>((k_b + k_step -1) - (k_b + 4)) + v.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))];
                    w.x[i][j][k] = (w.x[i][j][std::max(0, std::min(k_b + k_step +1, m.km - 1))] - w.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))]) * static_cast<double>(k - (k_b + 4)) / static_cast<double>((k_b + k_step -1) - (k_b + 4)) + w.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))];
                }
            }
        }


    // Coastal currents
    // Indian Ocean
    // South Africa/Madagascar     East coasts
    // Madagascar Current in the east (from j=103 until j=125 compares to 13°S until 35°S,
    //                                                  from k=45 until k=54 compares to 45°O until 54°O)
        j_beg = 103;
        j_end = 126;
        k_beg = 45;
        k_end = 55;
        k_step = 20;

        v_grad = + 0.001;

        k_a = k_b = 0;

        flip = 0;

        for (int k = std::min(k_end, m.km - 1); k > k_beg; k--) {
            if (h.x[i_max][j_beg][k] == 1.) {
                k_a = k;
                flip = 1;
            }
            if (flip == 1) break;
        }

        flip = 0;

        for (int j = j_beg+1; j < j_end; j++) {
            for (int k = std::min(k_end, m.km - 1); k > k_beg; k--) {
                if (h.x[i_max][j][k] == 1.) {
                    k_b = k;
                    k_grad = k_a - k_b;
                    if (k_grad >= 2) k_grad = 1;
                    if (k_grad <= 0) k_grad = 1;
                    flip = 1;
                }
            if (flip == 1) break;
            }

            for (int k = k_b; k < (k_b + k_step) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg; i < m.im; i++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = + v_grad * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = - v_grad * static_cast<double>(k_grad) * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }
            k_a = k_b;
            flip = 0;

            for (int k = k_b + 4; k < (k_b + k_step) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg; i < m.im; i++) {
                    v.x[i][j][k] = (v.x[i][j][std::max(0, std::min(k_b + k_step +1, m.km - 1))] - v.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))]) * static_cast<double>(k - (k_b + 4)) / static_cast<double>((k_b + k_step -1) - (k_b + 4)) + v.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))];
                    w.x[i][j][k] = (w.x[i][j][std::max(0, std::min(k_b + k_step +1, m.km - 1))] - w.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))]) * static_cast<double>(k - (k_b + 4)) / static_cast<double>((k_b + k_step -1) - (k_b + 4)) + w.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))];
                }
            }
        }


    // Indian Ocean
    // Agulhas Current at the southern tip of South Africa (from j=125 until j=128 compares to 35°S until 38°S,
    // parallel to the equator                                               from k=18 until k=27 compares to 18°O until 27°O)
        #pragma omp parallel for collapse(3)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = 125; j < 129; j++) {
                for (int k = 18; k < 28; k++) {
                    if (h.x[i][j][k] == 0.) {
                        {
                            v.x[i][j][k] = + 0.0;
                            w.x[i][j][k] = - IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        }
                    }
                }
            }
        }


    // Indian Ocean
    // Agulhas Current at the southern tip of South Africa (from j=125 until j=128 compares to 35°S until 38°S,
    // western transition parallel to the equator             from k=18 until k=27 compares to 18°O until 27°O)
        j_beg = 125;
        j_end = 129;
        k_beg = 18;
        k_end = 28;

        #pragma omp parallel for collapse(3) private(k_z, k_n)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    k_z = k - k_beg;
                    k_n = k_end - k_beg;

                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = (v.x[i][j][k_end] - v.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + v.x[i][j][k_beg];
                        w.x[i][j][k] = (w.x[i][j][k_end] - w.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + w.x[i][j][k_beg];
                    }
                }
            }
        }


    // Indian Ocean
    // Agulhas Current at the southern tip of South Africa (from j=125 until j=128 compares to 35°S until 38°S,
    // eastern transition parallel to the equator                from k=18 until k=27 compares to 18°O until 27°O)
        j_beg = 125;
        j_end = 129;
        k_beg = 18;
        k_end = 28;

        #pragma omp parallel for collapse(3) private(k_z, k_n)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    k_z = k - k_beg;
                    k_n = k_end - k_beg;

                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = (v.x[i][j][k_end] - v.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + v.x[i][j][k_beg];
                        w.x[i][j][k] = (w.x[i][j][k_end] - w.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + w.x[i][j][k_beg];
                    }
                }
            }
        }


    // Indian Ocean
    // Agulhas Current at the southern tip of South Africa (from j=125 until j=128 compares to 35°S until 38°S,
    // southern transition parallel to the equator               from k=18 until k=27 compares to 18°O until 27°O)
        j_beg = 125;
        j_end = 129;
        k_beg = 18;
        k_end = 28;

        #pragma omp parallel for collapse(3) private(j_z, j_n)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    j_z = j - j_beg;
                    j_n = j_end - j_beg;

                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = (v.x[i][j_end][k] - v.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + v.x[i][j_beg][k];
                        w.x[i][j][k] = (w.x[i][j_end][k] - w.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + w.x[i][j_beg][k];
                    }
                }
            }
        }


    // Indian Ocean
    // North Africa       East coasts
    // Somali Current in the north-east of North Africa (from j=80 until j=100 compares to 10°N until 10°S,
    //                                                                           from k=38 until k=58 compares to 38°O until 58°O)
        j_beg = 80;
        j_end = 101;
        k_beg = 38;
        k_end = 59;
        k_step = 20;

        v_grad = + 0.001;

        k_a = k_b = 0;

        flip = 0;

        for (int k = std::min(k_end, m.km - 1); k > k_beg; k--) {
            if (h.x[i_max][j_beg][k] == 1.) {
                k_a = k;
                flip = 1;
            }
            if (flip == 1) break;
        }

        flip = 0;

        for (int j = j_beg+1; j < j_end; j++) {
            for (int k = std::min(k_end, m.km - 1); k > k_beg; k--) {
                if (h.x[i_max][j][k] == 1.) {
                    k_b = k;
                    k_grad = k_a - k_b;
                    if (k_grad >= 2) k_grad = 1;
                    if (k_grad <= 0) k_grad = 1;
                    flip = 1;
                }
            if (flip == 1) break;
            }

            for (int k = k_b; k < (k_b + k_step) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg; i < m.im; i++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = - v_grad * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = + v_grad * static_cast<double>(k_grad) * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }
            k_a = k_b;
            flip = 0;

            for (int k = k_b + 4; k < (k_b + k_step) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg; i < m.im; i++) {
                    v.x[i][j][k] = (v.x[i][j][std::max(0, std::min(k_b + k_step +1, m.km - 1))] - v.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))]) * static_cast<double>(k - (k_b + 4)) / static_cast<double>((k_b + k_step -1) - (k_b + 4)) + v.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))];
                    w.x[i][j][k] = (w.x[i][j][std::max(0, std::min(k_b + k_step +1, m.km - 1))] - w.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))]) * static_cast<double>(k - (k_b + 4)) / static_cast<double>((k_b + k_step -1) - (k_b + 4)) + w.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))];
                }
            }
        }


    // Indian Ocean
    // North Africa       East coasts
    // Connection between Somali and Agulhas Current
    // Somali Current in the north-east of North Africa (from j=80 until j=100 compares to 10°N until 10°S,
    //                                                                           from k=38 until k=58 compares to 38°O until 58°O)
        j_beg = 81;
        j_end = 101;
        k_beg = 38;
        k_end = 59;

        #pragma omp parallel for collapse(3) private(j_z, j_n)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    j_z = j - j_beg;
                    j_n = j_end - j_beg;

                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = (v.x[i][j_end][k] - v.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + v.x[i][j_beg][k];
                        w.x[i][j][k] = (w.x[i][j_end][k] - w.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + w.x[i][j_beg][k];
                    }
                }
            }
        }


    // Indian Ocean
    // Indonesia       West coast
    // Sumatra/Java Current in the west of Indonesia (from j=85 until j=99 compares to 5°N until 9°S,
    //                                                                                 from k=95 until k=105 compares to 95°O until 105°O)
        j_beg = 85;
        j_end = 100;
        k_beg = 95;
        k_end = 106;

        k_exp = 5;
        k_w = 4;

        v_grad = + 0.0008;
        k_step = 8;

        k_a = k_b = 0;

        flip = 0;

        for (int k = k_beg; k < k_end; k++) {
            if (h.x[i_max][j_beg][k] == 1.) {
                k_a = k;
                flip = 1;
            }
            if (flip == 1) break;
        }

        flip = 0;

        for (int j = j_beg+1; j < j_end; j++) {
            for (int k = k_beg; k < k_end; k++) {
                if (h.x[i_max][j][k] == 1.) {
                    k_b = k;
                    k_grad = k_a - k_b;
                    if (k_grad >= 2) k_grad = 1;
                    if (k_grad <= 0) k_grad = 1;
                    flip = 1;
                }
            if (flip == 1) break;
            }

            for (int k = k_b; k > (k_b - k_step) && k >= 0; k--) {
                for (int i = i_beg; i < m.im; i++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = + v_grad * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = + v_grad * static_cast<double>(k_grad) * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }
            k_a = k_b;
            flip = 0;

            for (int k = k_b - k_step - k_exp; k < (k_b - k_w) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg; i < m.im; i++) {
                    v.x[i][j][k] = (v.x[i][j][std::max(k_b - k_w, 0)] - v.x[i][j][std::max(0, std::min(k_b - k_step - k_exp, m.km - 1))]) * static_cast<double>(k - (k_b - k_step - k_exp)) / static_cast<double>((k_b - k_w) - (k_b - k_step - k_exp)) + v.x[i][j][std::max(0, std::min(k_b - k_step - k_exp, m.km - 1))];
                    w.x[i][j][k] = (w.x[i][j][std::max(k_b - k_w, 0)] - w.x[i][j][std::max(0, std::min(k_b - k_step - k_exp, m.km - 1))]) * static_cast<double>(k - (k_b - k_step - k_exp)) / static_cast<double>((k_b - k_w) - (k_b - k_step - k_exp)) + w.x[i][j][std::max(0, std::min(k_b - k_step - k_exp, m.km - 1))];
                }
            }
        }


    // Indian Ocean
    // Australia       West coast
    // West Australian Current in the west of Australia (from j=112 until j=125 compares to 22°S until 35°S,
    //                                                                                 from k=103 until k=115 compares to 103°O until 115°O)
        j_beg = 112;
        j_end = 126;
        k_beg = 103;
        k_end = 116;

        k_exp = 5;
        k_w = 4;

        v_grad = 0.0008;
        k_step = 8;

        k_a = k_b = 0;

        flip = 0;

        for (int k = k_beg; k < k_end; k++) {
            if (h.x[i_max][j_beg][k] == 1.) {
                k_a = k;
                flip = 1;
            }
            if (flip == 1) break;
        }

        flip = 0;

        for (int j = j_beg+1; j < j_end; j++) {
            for (int k = k_beg; k < k_end; k++) {
                if (h.x[i_max][j][k] == 1.) {
                    k_b = k;
                    k_grad = k_a - k_b;
                    if (k_grad >= 2) k_grad = 1;
                    if (k_grad <= 0) k_grad = 1;
                    flip = 1;
                }
            if (flip == 1) break;
            }

            for (int k = k_b; k > (k_b - k_step) && k >= 0; k--) {
                for (int i = i_beg; i < m.im; i++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = + v_grad * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = + 0.1 * v_grad * static_cast<double>(k_grad) * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }
            k_a = k_b;
            flip = 0;

            for (int k = k_b - k_step - k_exp; k < (k_b - k_w) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg; i < m.im; i++) {
                    v.x[i][j][k] = (v.x[i][j][std::max(k_b - k_w, 0)] - v.x[i][j][std::max(0, std::min(k_b - k_step - k_exp, m.km - 1))]) * static_cast<double>(k - (k_b - k_step - k_exp)) / static_cast<double>((k_b - k_w) - (k_b - k_step - k_exp)) + v.x[i][j][std::max(0, std::min(k_b - k_step - k_exp, m.km - 1))];
                    w.x[i][j][k] = (w.x[i][j][std::max(k_b - k_w, 0)] - w.x[i][j][std::max(0, std::min(k_b - k_step - k_exp, m.km - 1))]) * static_cast<double>(k - (k_b - k_step - k_exp)) / static_cast<double>((k_b - k_w) - (k_b - k_step - k_exp)) + w.x[i][j][std::max(0, std::min(k_b - k_step - k_exp, m.km - 1))];
                }
            }
        }
    }
/*
* 
*/
    void IC_Pazifischer_Ozean (Array &h, Array &u, Array &v, Array &w) {
    // Thermohalin Conveyor Belt and currents along the coasts
    // Closing the hydrospheric circulation systems
    // Pacific Ocean
    // Thermohalin Conveyor Belt
    // Surface flow west of California (from j=52 until j=80 compares to 10°N until 38°N
    //                                                                  from k=220 until k=250 compares to 110°W until 140°W)
        j_beg = 52;
        j_end = 81;
        k_beg = 220;
        k_end = 251;

        v_grad = 0.0010;
        k_step = 8;

        k_a = k_b = 0;

        flip = 0;

        for (int k = k_beg; k < k_end; k++) {
            if (h.x[i_max][j_beg][k] == 1.) {
                k_a = k;
                flip = 1;
            }
            if (flip == 1) break;
        }

        flip = 0;

        for (int j = j_beg+1; j < j_end; j++) {
            for (int k = k_beg; k < k_end; k++) {
                if (h.x[i_max][j][k] == 1.) {
                    k_b = k;
                    k_grad = k_b - k_a;
                    if (k_grad <= - 3) k_grad = - 2;
                    flip = 1;
                }
            if (flip == 1) break;
            }

            for (int k = k_b; k > (k_b - k_step) && k >= 0; k--) {
                for (int i = i_beg; i < m.im; i++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = + v_grad * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = + v_grad * static_cast<double>(k_grad) * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }
            k_a = k_b;
            flip = 0;
        }


    // Pacific Ocean
    // Thermohalin Conveyor Belt
    // Surface flow south of the Aleutians (from j=50 until j=56 compares to 34°N until 40°N
    //                                           from k=169 until k=230 compares to 169°O until 130°W)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = 50; j < 57; j++) {
                for (int k = 169; k < 231; k++) {
                    v.x[i][j][k] = + 0.;
                    w.x[i][j][k] = + IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                }
            }
        }


    // Pacific
    // Current from Japan until California in North America
    // T-junction
    // Eastward flowing current as surface flow (from j=50 until j=55 compares to 35°N until 40°N,
    //                                                                                           from k=160 until k=220 compares to 160°O until 140°W)
        #pragma omp parallel for collapse(3)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = 50; j < 56; j++) {
                for (int k = 160; k < 221; k++) {
                    v.x[i][j][k] = + IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    w.x[i][j][k] = + IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                }
            }
        }


    // Pacific
    // Current from Japan until California in North America
    // T-junction
    // Eastward flowing current as surface flow (from j=50 until j=55 compares to 35°N until 40°N,
    // western transition parallel to the equator                       from k=160 until k=220 compares to 160°O until 140°W)

        j_beg = 50;
        j_end = 56;
        k_beg = 160;
        k_end = 221;

        #pragma omp parallel for collapse(3) private(k_z, k_n)
        for (int i = 5; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    k_z = k - k_beg;
                    k_n = k_end - k_beg;

                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = (v.x[i][j][k_end] - v.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + v.x[i][j][k_beg];
                        w.x[i][j][k] = (w.x[i][j][k_end] - w.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + w.x[i][j][k_beg];
                    }
                }
            }
        }


    // Pacific
    // Current from Japan until California in North America
    // T-junction
    // Eastward flowing current as surface flow (from j=50 until j=55 compares to 35°N until 40°N,
    // eastern transition                                                            from k=160 until k=220 compares to 160°O until 140°W)
        j_beg = 50;
        j_end = 56;
        k_beg = 160;
        k_end = 221;

        #pragma omp parallel for collapse(3) private(k_z, k_n)
        for (int i = 5; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    k_z = k - k_beg;
                    k_n = k_end - k_beg;

                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = (v.x[i][j][k_end] - v.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + v.x[i][j][k_beg];
                        w.x[i][j][k] = (w.x[i][j][k_end] - w.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + w.x[i][j][k_beg];
                    }
                }
            }
        }


    // Pacific
    // Current from Japan until California in North America
    // T-junction
    // Eastward flowing current as surface flow (from j=50 until j=55 compares to 35°N until 40°N,
    // northern transition                                                         from k=160 until k=220 compares to 160°O until 140°W)
        j_beg = 50;
        j_end = 56;
        k_beg = 160;
        k_end = 221;

        #pragma omp parallel for collapse(3) private(j_z, j_n)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    j_z = j - j_beg;
                    j_n = j_end - j_beg;

                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = (v.x[i][j_end][k] - v.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + v.x[i][j_beg][k];
                        w.x[i][j][k] = (w.x[i][j_end][k] - w.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + w.x[i][j_beg][k];
                    }
                }
            }
        }


    // Current from Japan until California in North America
    // T-junction
    // Eastward flowing current as surface flow (from j=50 until j=55 compares to 35°N until 40°N,
    // southern transition                                                           from k=160 until k=220 compares to 160°O until 140°W)
        j_beg = 50;
        j_end = 56;
        k_beg = 160;
        k_end = 221;

        #pragma omp parallel for collapse(3) private(j_z, j_n)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    j_z = j - j_beg;
                    j_n = j_end - j_beg;

                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = (v.x[i][j_end][k] - v.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + v.x[i][j_beg][k];
                        w.x[i][j][k] = (w.x[i][j_end][k] - w.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + w.x[i][j_beg][k];
                    }
                }
            }
        }


    // Pacific
    // Current from the Aleutians until the eastward-flowing current
    // T-junction
    // Eastward flowing current as surface flow (from j=27 until j=50 compares to 40°N until 63°N,
    // north-south, east branch                                                         from k=195 until k=200 compares to 165°W until 160°W)

        #pragma omp parallel for collapse(3)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = 27; j < 51; j++) {
                for (int k = 195; k < 201; k++) {
                    v.x[i][j][k] = + IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    w.x[i][j][k] = + IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                }
            }
        }


    // Pacific
    // Current from the eastward-flowing current until the Aleutians
    // T-junction
    // Eastward flowing current as surface flow (from j=27 until j=50 compares to 40°N until 63°N,
    // south-north, west branch                                                      from k=195 until k=200 compares to 170°W until 165°W)
   
        #pragma omp parallel for collapse(3)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = 27; j < 51; j++) {
                for (int k = 195; k < 201; k++) {
                    v.x[i][j][k] = - IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    w.x[i][j][k] = + IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                }
            }
        }


    // Pacific
    // Current from the eastward-flowing current until the Aleutians
    // Current from the Aleutians until the eastward-flowing current
    // T-junction
    // Eastward flowing current as surface flow (from j=40 until j=50 compares to 40°N until 50°N,
    // transition                                                                          from k=194 until k=198 compares to 166°W until 162°W)
        j_beg = 40;
        j_end = 51;
        k_beg = 194;
        k_end = 199;

        #pragma omp parallel for collapse(3) private(k_z, k_n)
        for (int i = 5; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    k_z = k - k_beg;
                    k_n = k_end - k_beg;

                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = (v.x[i][j][k_end] - v.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + v.x[i][j][k_beg];
                        w.x[i][j][k] = (w.x[i][j][k_end] - w.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + w.x[i][j][k_beg];
                    }
                }
            }
        }


    // Pacific
    // Current from Japan until California in North America
    // north-south
    // T-junction
    // Eastward flowing current as surface flow (from j=27 until j=50 compares to 40°N until 63°N,
    // north-south, east branch eastern transition                         from k=198 until k=204 compares to 162°W until 156°W)
        j_beg = 27;
        j_end = 51;
        k_beg = 198;
        k_end = 205;

        #pragma omp parallel for collapse(3) private(k_z, k_n)
        for (int i = 5; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    k_z = k - k_beg;
                    k_n = k_end - k_beg;

                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = (v.x[i][j][k_end] - v.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + v.x[i][j][k_beg];
                        w.x[i][j][k] = (w.x[i][j][k_end] - w.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + w.x[i][j][k_beg];
                    }
                }
            }
        }


    // Pazifik
    // Current from Japan until California in North America
    // south-north
    // T-junction
    // Eastward flowing current as surface flow (from j=27 until j=50 compares to 40°N until 63°N,
    // south-north, west branch  western transition                   from k=183 until k=192 compares to 177°W until 168°W)
        j_beg = 27;
        j_end = 51;
        k_beg = 183;
        k_end = 193;

        #pragma omp parallel for collapse(3) private(k_z, k_n)
        for (int i = 5; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    k_z = k - k_beg;
                    k_n = k_end - k_beg;

                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = (v.x[i][j][k_end] - v.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + v.x[i][j][k_beg];
                        w.x[i][j][k] = (w.x[i][j][k_end] - w.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + w.x[i][j][k_beg];
                    }
                }
            }
        }


    // Pazifischer Ozean
    // Thermohalin Conveyor Belt
    // southerly, parallel to the equator from Indonesia/Australia until California in North America
    // Westward flowing current as surface flow (from j=102 until j=107 compares to 12°S until 17°S,
    //                                                                                     from k=170 until k=260 compares to 170°O until 90°W)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = 102; j < 108; j++) {
                for (int k = 701; k < 261; k++) {
                    v.x[i][j][k] = + IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    w.x[i][j][k] = - IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                }
            }
        }


    // Coastal currents
    // Pacific
    // South America       West coasts
    // Peru or Humboldt Current in the west of South America (from j=95 until j=130 compares to 5°S until 40°S,
    // extension as surface flow                                    from k=270 until k=290 compares to 70°W until 90°W)
        j_beg = 95;
        j_end = 131;
        k_beg = 270;
        k_end = 291;

        v_grad = + 0.0010;
        k_step = 8;

        k_a = k_b = 0;

        flip = 0;

        for (int k = k_beg; k < k_end; k++) {
            if (h.x[i_max][j_beg][k] == 1.) {
                k_a = k;
                flip = 1;
            }
            if (flip == 1) break;
        }

        flip = 0;

        for (int j = j_beg+1; j < j_end; j++) {
            for (int k = k_beg; k < k_end; k++) {
                if (h.x[i_max][j][k] == 1.) {
                    k_b = k;
                    k_grad = k_b - k_a;
                    if (k_grad <= - 3) k_grad = - 2;
                    flip = 1;
                }
            if (flip == 1) break;
            }

            for (int k = k_b; k > (k_b - k_step) && k >= 0; k--) {
                for (int i = i_beg; i < m.im; i++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = + v_grad * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = + v_grad * static_cast<double>(k_grad) * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }
            k_a = k_b;
            flip = 0;
        }


    // Pazifik
    // Central America        West coasts
    // "Central America" current in the west of Central America (West Indies) (from j=60 until j=85 compares to 5°N until 30°N,
    //                                                                                                              from k=245 until k=275 compares to 85°W until 115°W)
        j_beg = 60;
        j_end = 86;
        k_beg = 245;
        k_end = 276;

        v_grad = 0.0010;
        k_step = 6;

        k_a = k_b = 0;

        flip = 0;

        for (int k = k_beg; k < k_end; k++) {
            if (h.x[i_max][j_beg][k] == 1.) {
                k_a = k;
                flip = 1;
            }
            if (flip == 1) break;
        }

        flip = 0;

        for (int j = j_beg+1; j < j_end; j++) {
            for (int k = k_beg; k < k_end; k++) {
                if (h.x[i_max][j][k] == 1.) {
                    k_b = k;
                    k_grad = k_b - k_a;
                    if (k_grad <= - 3) k_grad = - 2;
                    flip = 1;
                }
            if (flip == 1) break;
            }

            for (int k = k_b; k > (k_b - k_step) && k >= 0; k--) {
                for (int i = i_beg; i < m.im; i++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = + v_grad * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = + v_grad * static_cast<double>(k_grad) * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }
            k_a = k_b;
            flip = 0;
        }


    // Pazifik
    // North America      West coasts
    // Current in the Gulf of Alaska (from j=29 until j=50 compares to 40°N until 61°N,
    //                                              from k=215 until k=240 compares to 145°W until 120°W)
        j_beg = 29;
        j_end = 51;
        k_beg = 215;
        k_end = 241;

        k_exp = 10;
        k_w = 4;

        v_grad = 0.0008;
        k_step = 20;

        k_a = k_b = 0;

        flip = 0;

        for (int k = k_beg; k < k_end; k++) {
            if (h.x[i_max][j_beg][k] == 1.) {
                k_a = k;
                flip = 1;
            }
            if (flip == 1) break;
        }

        flip = 0;

        for (int j = j_beg+1; j < j_end; j++) {
            for (int k = k_beg; k < k_end; k++) {
                if (h.x[i_max][j][k] == 1.) {
                    k_b = k;
                    k_grad = k_a - k_b;
                    if (k_grad >= 2) k_grad = 1;
                    if (k_grad <= 0) k_grad = 1;
                    flip = 1;
                }
            if (flip == 1) break;
            }

            for (int k = k_b; k > (k_b - k_step) && k >= 0; k--) {
                for (int i = i_beg; i < m.im; i++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = - v_grad * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = - v_grad * static_cast<double>(k_grad) * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }
            k_a = k_b;
            flip = 0;

            for (int k = k_b - k_step - k_exp; k < (k_b - k_w) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg; i < m.im; i++) {
                    v.x[i][j][k] = (v.x[i][j][std::max(k_b - k_w, 0)] - v.x[i][j][std::max(0, std::min(k_b - k_step - k_exp, m.km - 1))]) * static_cast<double>(k - (k_b - k_step - k_exp)) / static_cast<double>((k_b - k_w) - (k_b - k_step - k_exp)) + v.x[i][j][std::max(0, std::min(k_b - k_step - k_exp, m.km - 1))];
                    w.x[i][j][k] = (w.x[i][j][std::max(k_b - k_w, 0)] - w.x[i][j][std::max(0, std::min(k_b - k_step - k_exp, m.km - 1))]) * static_cast<double>(k - (k_b - k_step - k_exp)) / static_cast<double>((k_b - k_w) - (k_b - k_step - k_exp)) + w.x[i][j][std::max(0, std::min(k_b - k_step - k_exp, m.km - 1))];
                }
            }
        }


    // Pazifik
    // Connection between Gulf of Alaska current and NEC current
    // North-West America        West coasts
    // Current in the Gulf of Alaska (from j=29 until j=50 compares to 40°N until 61°N,
    //                                              from k=215 until k=240 compares to 145°W until 120°W)
        j_beg = 29;
        j_end = 51;
        k_beg = 215;
        k_end = 241;

            #pragma omp parallel for collapse(3) private(j_z, j_n)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    j_z = j - j_beg;
                    j_n = j_end - j_beg;

                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = (v.x[i][j_end][k] - v.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + v.x[i][j_beg][k];
                        w.x[i][j][k] = (w.x[i][j_end][k] - w.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + w.x[i][j_beg][k];
                    }
                }
            }
        }


    // Pazifik
    // North-East Asia       East coasts
    // Kamchatka Current in the north-east of Japan (from j=25 until j=40 compares to 65°N until 50°N,
    //                                                                            from k=159 until k=180 compares to 159°O until 180°O)
        j_beg = 25;
        j_end = 41;
        k_beg =159;
        k_end = 181;
        k_step = 15;

        v_grad = + 0.001;

        k_a = k_b = 0;

        flip = 0;

        for (int k = std::min(k_end, m.km - 1); k > k_beg; k--) {
            if (h.x[i_max][j_beg][k] == 1.) {
                k_a = k;
                flip = 1;
            }
            if (flip == 1) break;
        }

        flip = 0;

        for (int j = j_beg+1; j < j_end; j++) {
            for (int k = std::min(k_end, m.km - 1); k > k_beg; k--) {
                if (h.x[i_max][j][k] == 1.) {
                    k_b = k;
                    k_grad = k_a - k_b;
                    if (k_grad >= 2) k_grad = 1;
                    if (k_grad <= 0) k_grad = 1;
                    flip = 1;
                }
            if (flip == 1) break;
            }

            for (int k = k_b; k < (k_b + k_step) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg; i < m.im; i++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = + 0.5 * v_grad * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = - v_grad * static_cast<double>(k_grad) * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }

            for (int k = k_b + 4; k < (k_b + k_step) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg; i < m.im; i++) {
                    v.x[i][j][k] = (v.x[i][j][std::max(0, std::min(k_b + k_step +1, m.km - 1))] - v.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))]) * static_cast<double>(k - (k_b + 4)) / static_cast<double>((k_b + k_step -1) - (k_b + 4)) + v.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))];
                    w.x[i][j][k] = (w.x[i][j][std::max(0, std::min(k_b + k_step +1, m.km - 1))] - w.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))]) * static_cast<double>(k - (k_b + 4)) / static_cast<double>((k_b + k_step -1) - (k_b + 4)) + w.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))];
                }
            }
            k_a = k_b;
            flip = 0;
        }


    // Pazifik
    // North-East Asia       East coasts
    // Oyashio Current in the north-east of Japan (from j=39 until j=48 compares to 51°N until 42°N,
    //                                                                    from k=142 until k=161 compares to 142°O until 161°O)

        j_beg = 39;
        j_end = 49;
        k_beg =142;
        k_end = 162;
        k_step = 30;

        v_grad = + 0.001;

        k_a = k_b = 0;

        flip = 0;

        for (int k = std::min(k_end, m.km - 1); k > k_beg; k--) {
            if (h.x[i_max][j_beg][k] == 1.) {
                k_a = k;
                flip = 1;
            }
            if (flip == 1) break;
        }

        flip = 0;

        for (int j = j_beg+1; j < j_end; j++) {
            for (int k = std::min(k_end, m.km - 1); k > k_beg; k--) {
                if (h.x[i_max][j][k] == 1.) {
                    k_b = k;
                    k_grad = k_a - k_b;
                    if (k_grad >= 2) k_grad = 1;
                    if (k_grad <= 0) k_grad = 1;
                    flip = 1;
                }
            if (flip == 1) break;
            }

            for (int k = k_b; k < (k_b + k_step) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg; i < m.im; i++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = + v_grad * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = - v_grad * static_cast<double>(k_grad) * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }

            for (int k = k_b + 4; k < (k_b + k_step) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg; i < m.im; i++) {
                    v.x[i][j][k] = (v.x[i][j][std::max(0, std::min(k_b + k_step +1, m.km - 1))] - v.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))]) * static_cast<double>(k - (k_b + 4)) / static_cast<double>((k_b + k_step -1) - (k_b + 4)) + v.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))];
                    w.x[i][j][k] = (w.x[i][j][std::max(0, std::min(k_b + k_step +1, m.km - 1))] - w.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))]) * static_cast<double>(k - (k_b + 4)) / static_cast<double>((k_b + k_step -1) - (k_b + 4)) + w.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))];
                }
            }
            k_a = k_b;
            flip = 0;
        }


    // Pazifik
    // North-East Asia       East coasts
    // Kuroshio Current in the south-east of Japan (from j=48 until j=68 compares to 22°N until 42°N,
    //                                                                    from k=120 until k=142 compares to 120°O until 142°O)
        j_beg = 48;
        j_end = 69;
        k_beg = 120;
        k_end = 143;
        k_step = 30;

        v_grad = + 0.001;

        k_a = k_b = 0;

        flip = 0;

        for (int k = std::min(k_end, m.km - 1); k > k_beg; k--) {
            if (h.x[i_max][j_beg][k] == 1.) {
                k_a = k;
                flip = 1;
            }
            if (flip == 1) break;
        }

        flip = 0;

        for (int j = j_beg+1; j < j_end; j++) {
            for (int k = std::min(k_end, m.km - 1); k > k_beg; k--) {
                if (h.x[i_max][j][k] == 1.) {
                    k_b = k;
                    k_grad = k_a - k_b;
                    if (k_grad >= 2) k_grad = 1;
                    if (k_grad <= 0) k_grad = 1;
                    flip = 1;
                }
            if (flip == 1) break;
            }


            for (int k = k_b; k < (k_b + k_step) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg; i < m.im; i++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = - v_grad * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = + v_grad * static_cast<double>(k_grad) * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }
            k_a = k_b;
            flip = 0;

            for (int k = k_b + 4; k < (k_b + k_step) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg; i < m.im; i++) {
                    v.x[i][j][k] = (v.x[i][j][std::max(0, std::min(k_b + k_step +1, m.km - 1))] - v.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))]) * static_cast<double>(k - (k_b + 4)) / static_cast<double>((k_b + k_step -1) - (k_b + 4)) + v.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))];
                    w.x[i][j][k] = (w.x[i][j][std::max(0, std::min(k_b + k_step +1, m.km - 1))] - w.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))]) * static_cast<double>(k - (k_b + 4)) / static_cast<double>((k_b + k_step -1) - (k_b + 4)) + w.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))];
                }
            }
        }


    // Pazifik
    // Connection between Kamchatka and Oyashio current
    // North-East Asia       East coasts
    // Kamchatka and Oyashio current in the north of Japan (from j=45 until j=75 compares to 15°N until 45°N,
    //                                                                                              from k=120 until k=160 compares to 120°O until 160°O)
        j_beg = 45;
        j_end = 76;
        k_beg = 120;
        k_end = 161;

        #pragma omp parallel for collapse(3) private(j_z, j_n)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    j_z = j - j_beg;
                    j_n = j_end - j_beg;

                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = (v.x[i][j_end][k] - v.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + v.x[i][j_beg][k];
                        w.x[i][j][k] = (w.x[i][j_end][k] - w.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + w.x[i][j_beg][k];
                    }
                }
            }
        }


    // Pazifik
    // Connection between Oyashio and Kuroshio current
    // North-East Asia       East coasts
    // Oyashio and Kuroshio current in the east of Japan (from j=45 until j=75 compares to 15°N until 45°N,
    //                                                                                     from k=120 until k=160 compares to 120°O until 160°O)
        j_beg = 45;
        j_end = 76;
        k_beg = 120;
        k_end = 161;

       #pragma omp parallel for collapse(3) private(j_z, j_n)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    j_z = j - j_beg;
                    j_n = j_end - j_beg;

                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = (v.x[i][j_end][k] - v.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + v.x[i][j_beg][k];
                        w.x[i][j][k] = (w.x[i][j_end][k] - w.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + w.x[i][j_beg][k];
                    }
                }
            }
        }


    // Pazifik
    // Connection between Kuroshio and East China Sea
    // Kuroshio and East China Sea       East coasts
    // Kuroshio current and East China Sea in the east of Japan (from j=45 until j=75 compares to 15°N until 45°N,
    //                                                                                                          from k=120 until k=160 compares to 120°O until 160°O)
        j_beg = 45;
        j_end = 76;
        k_beg = 120;
        k_end = 161;

        #pragma omp parallel for collapse(3) private(k_z, k_n)
        for (int i = 5; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    k_z = k - k_beg;
                    k_n = k_end - k_beg;

                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = (v.x[i][j][k_end] - v.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + v.x[i][j][k_beg];
                        w.x[i][j][k] = (w.x[i][j][k_end] - w.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + w.x[i][j][k_beg];
                    }
                }
            }
        }


    // Connection between Kuroshio and NEC current
    // North-East Asia       East coasts
    // Oyashio and Kuroshio current in the east of Japan (from j=45 until j=75 compares to 15°N until 45°N,
    //                                                                                     from k=120 until k=160 compares to 120°O until 160°O)
        j_beg = 45;
        j_end = 76;
        k_beg = 120;
        k_end = 161;

        #pragma omp parallel for collapse(3) private(j_z, j_n)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = j_beg; j < j_end + 1; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    j_z = j - j_beg;
                    j_n = j_end - j_beg;

                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = (v.x[i][j_end][k] - v.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + v.x[i][j_beg][k];
                        w.x[i][j][k] = (w.x[i][j_end][k] - w.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + w.x[i][j_beg][k];
                    }
                }
            }
        }


    // Pacific
    // North-East Asia       East coasts
    // South China Sea in the east of Vietnam (from j=60 until j=80 compares to 23°N until 0°,
    //                                             from k=92 until k=107 compares to 104°O until 120°O)
        #pragma omp parallel for collapse(3)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = 60; j < 81; j++) {
                for (int k = 92; k < 108; k++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = + IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = - IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }
        }


    // Pacific
    // North-East Asia       East coasts
    // Tsushima Current in the west of Japan until Indonesia (from j=40 until j=53 compares to 45°N until 30°N,
    //                                                     from k=116 until k= 129 compares to 130°O until 145°O)
        l = 1;

        for (int i = i_beg; i < m.im; i++) {
            for (int j = 40; j < 54; j++) {
                for (int k = 116; k < 130; k++) {
                    if (h.x[i][j][k] == 0.) {
                        if (l <= 3) {
                            v.x[i][j][k] = - IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                            w.x[i][j][k] = - IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                            l++;
                        }
                    }
                if (h.x[i][j][k] == 1.) l = 1;
                }
            }
        }


    // Pacific
    // Philippines (Indonesia)         East coasts
    // Current in the east of the Philippines (from j=60 until j=71 compares to 10°N until 22°N,
    //                                  from k=111 until k=116 compares to 125°O until 130°O)
        j_beg = 60;
        j_end = 72;
        k_beg = 111;
        k_end = 117;
        k_step = 20;

        v_grad = 0.001;

        k_a = k_b = 0;

        flip = 0;

        for (int k = std::min(k_end, m.km - 1); k > k_beg; k--) {
            if (h.x[i_max][j_beg][k] == 1.) {
                k_a = k;
                flip = 1;
            }
            if (flip == 1) break;
        }

        flip = 0;

        for (int j = j_beg+1; j < j_end; j++) {
            for (int k = std::min(k_end, m.km - 1); k > k_beg; k--) {
                if (h.x[i_max][j][k] == 1.) {
                    k_b = k;
                    k_grad = k_a - k_b;
                    if (k_grad >= 2) k_grad = 1;
                    if (k_grad <= 0) k_grad = 1;
                    flip = 1;
                }
            if (flip == 1) break;
            }

            for (int k = k_b; k < (k_b + k_step) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg; i < m.im; i++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = - v_grad * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = - v_grad * static_cast<double>(k_grad) * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }
            k_a = k_b;
            flip = 0;
        }


    // Pacific
    // Philippines (Indonesia) New Guinea       East coasts
    // Throughflow between Indonesia and New Guinea (from j=80 until j=97 compares to 10°N until 7°S
    // continuation toward the Pacific                                   from k=125 until k=135 compares to 125°O until 135°O)

        #pragma omp parallel for collapse(3)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = 80; j < 98; j++) {
                for (int k = 125; k < 136; k++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = + IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = - IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }
        }


    // Pacific
    // Philippines (Indonesia) New Guinea       East coasts
    // Throughflow between Indonesia and New Guinea (from j=80 until j=97 compares to 10°N until 7°S
    // western transition parallel to the equator             from k=125 until k=135 compares to 125°O until 135°O)
        j_beg = 80;
        j_end = 98;
        k_beg = 125;
        k_end = 136;

        #pragma omp parallel for collapse(3) private(k_z, k_n)
        for (int i = 5; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    k_z = k - k_beg;
                    k_n = k_end - k_beg;

                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = (v.x[i][j][k_end] - v.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + v.x[i][j][k_beg];
                        w.x[i][j][k] = (w.x[i][j][k_end] - w.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + w.x[i][j][k_beg];
                    }
                }
            }
        }


    // Pacific
    // Philippines (Indonesia) New Guinea       East coasts
    // Throughflow between Indonesia and New Guinea (from j=80 until j=97 compares to 10°N until 7°S
    // eastern transition parallel to the equator                from k=125 until k=135 compares to 125°O until 135°O)

        j_beg = 80;
        j_end = 98;
        k_beg = 125;
        k_end = 136;

        #pragma omp parallel for collapse(3) private(k_z, k_n)
        for (int i = 5; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    k_z = k - k_beg;
                    k_n = k_end - k_beg;

                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = (v.x[i][j][k_end] - v.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + v.x[i][j][k_beg];
                        w.x[i][j][k] = (w.x[i][j][k_end] - w.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + w.x[i][j][k_beg];
                    }
                }
            }
        }


    // Pacific
    // Philippines (Indonesia) New Guinea       East coasts
    // Throughflow between Indonesia and New Guinea (from j=80 until j=97 compares to 10°N until 7°S
    // northern transition parallel to the equator              from k=120 until k=135 compares to 125°O until 135°O)

        j_beg = 80;
        j_end = 98;
        k_beg = 120;
        k_end = 136;

        #pragma omp parallel for collapse(3) private(j_z, j_n)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    j_z = j - j_beg;
                    j_n = j_end - j_beg;

                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = (v.x[i][j_end][k] - v.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + v.x[i][j_beg][k];
                        w.x[i][j][k] = (w.x[i][j_end][k] - w.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + w.x[i][j_beg][k];
                    }
                }
            }
        }


    // Pacific
    // Philippines (Indonesia) New Guinea       East coasts
    // Throughflow between Indonesia and New Guinea (from j=80 until j=97 compares to 10°N until 7°S
    // southern transition parallel to the equator               from k=125 until k=135 compares to 125°O until 135°O)
        j_beg = 80;
        j_end = 98;
        k_beg = 125;
        k_end = 136;

        #pragma omp parallel for collapse(3) private(j_z, j_n)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    j_z = j - j_beg;
                    j_n = j_end - j_beg;

                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = (v.x[i][j_end][k] - v.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + v.x[i][j_beg][k];
                        w.x[i][j][k] = (w.x[i][j_end][k] - w.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + w.x[i][j_beg][k];
                    }
                }
            }
        }


    // Pacific
    // Papua New Guinea (Indonesia)     North coasts of New Guinea
    // Current in the north-east of New Guinea (from j=90 until j=100 compares to 0°S until 10°S,
    //                                                                from k=135 until k=150 compares to 135°O until 150°O)
        j_beg = 90;
        j_end = 100;
        k_beg = 135;
        k_end = 151;
        k_step = 20;

        v_grad = + 0.001;

        k_a = k_b = 0;

        flip = 0;

        for (int k = std::min(k_end, m.km - 1); k > k_beg; k--) {
            if (h.x[i_max][j_beg][k] == 1.) {
                k_a = k;
                flip = 1;
            }
            if (flip == 1) break;
        }

        flip = 0;

        for (int j = j_beg+1; j < j_end; j++) {
            for (int k = std::min(k_end, m.km - 1); k > k_beg; k--) {
                if (h.x[i_max][j][k] == 1.) {
                    k_b = k;
                    k_grad = k_a - k_b;
                    if (k_grad >= 2) k_grad = 1;
                    if (k_grad <= 0) k_grad = 1;
                    flip = 1;
                }
            if (flip == 1) break;
            }

            for (int k = k_b; k < (k_b + k_step) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg; i < m.im; i++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = - v_grad * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = - v_grad * static_cast<double>(k_grad) * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }
            k_a = k_b;
            flip = 0;

            for (int k = k_b + 4; k < (k_b + k_step) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg; i < m.im; i++) {
                    v.x[i][j][k] = (v.x[i][j][std::max(0, std::min(k_b + k_step +1, m.km - 1))] - v.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))]) * static_cast<double>(k - (k_b + 4)) / static_cast<double>((k_b + k_step -1) - (k_b + 4)) + v.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))];
                    w.x[i][j][k] = (w.x[i][j][std::max(0, std::min(k_b + k_step +1, m.km - 1))] - w.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))]) * static_cast<double>(k - (k_b + 4)) / static_cast<double>((k_b + k_step -1) - (k_b + 4)) + w.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))];
                }
            }
        }


    // Pacific
    // Australia       East coasts in the Coral and Tasman Sea
    // Current in the east of Australia (from j=101 until j=140 compares to 11°S until 50°S,
    //                                               from k=142 until k=158 compares to 142°O until 158°O)

        j_beg = 101;
        j_end = 141;
        k_beg = 142;
        k_end = 159;
        k_step = 20;

        v_grad = + 0.001;

        k_a = k_b = 0;

        flip = 0;

        for (int k = std::min(k_end, m.km - 1); k > k_beg; k--) {
            if (h.x[i_max][j_beg][k] == 1.) {
                k_a = k;
                flip = 1;
            }
            if (flip == 1) break;
        }

        flip = 0;

        for (int j = j_beg+1; j < j_end; j++) {
            for (int k = std::min(k_end, m.km - 1); k > k_beg; k--) {
                if (h.x[i_max][j][k] == 1.) {
                    k_b = k;
                    k_grad = k_a - k_b;
                    if (k_grad >= 2) k_grad = 1;
                    if (k_grad <= 0) k_grad = 1;
                    flip = 1;
                }
            if (flip == 1) break;
            }

            for (int k = k_b; k < (k_b + k_step) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg; i < m.im; i++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = + v_grad * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = - 0.4 * v_grad * static_cast<double>(k_grad) * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }
            k_a = k_b;
            flip = 0;

            for (int k = k_b + 4; k < (k_b + k_step) && k >= 0 && k < m.km; k++) {
                for (int i = i_beg; i < m.im; i++) {
                    v.x[i][j][k] = (v.x[i][j][std::max(0, std::min(k_b + k_step +1, m.km - 1))] - v.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))]) * static_cast<double>(k - (k_b + 4)) / static_cast<double>((k_b + k_step -1) - (k_b + 4)) + v.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))];
                    w.x[i][j][k] = (w.x[i][j][std::max(0, std::min(k_b + k_step +1, m.km - 1))] - w.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))]) * static_cast<double>(k - (k_b + 4)) / static_cast<double>((k_b + k_step -1) - (k_b + 4)) + w.x[i][j][std::max(0, std::min(k_b + 4, m.km - 1))];
                }
            }
        }


    // Pacific
    // Connection in the north of Australia and New Zealand
    // Australia       East coasts in the Coral and Tasman Sea
    // Current in the east of Australia (from j=110 until j=128 compares to 20°S until 38°S,
    //                                                     from k=142 until k=158 compares to 142°O until 158°O)

        j_beg = 110;
        j_end = 129;
        k_beg = 142;
        k_end = 159;

        #pragma omp parallel for collapse(3) private(j_z, j_n)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    j_z = j - j_beg;
                    j_n = j_end - j_beg;

                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = (v.x[i][j_end][k] - v.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + v.x[i][j_beg][k];
                        w.x[i][j][k] = (w.x[i][j_end][k] - w.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + w.x[i][j_beg][k];
                    }
                }
            }
        }


    // Pacific
    // Australia       South coasts in the South Australian Basin
    // Current in the south of Australia (from j=123 until j=132 compares to 33°S until 42°S,
    //                                                      from k=116 until k=148 compares to 116°O until 148°O)
        #pragma omp parallel for collapse(3)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = 123; j < 133; j++) {
                for (int k = 116; k < 149; k++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = + IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = + IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }
        }


    // Connection between Australia and circumpolar current
    // Australia       South coast

        j_beg = 128;
        j_end = 133;
        k_beg = 116;
        k_end = 149;

        #pragma omp parallel for collapse(3) private(j_z, j_n)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    j_z = j - j_beg;
                    j_n = j_end - j_beg;

                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = (v.x[i][j_end][k] - v.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + v.x[i][j_beg][k];
                        w.x[i][j][k] = (w.x[i][j_end][k] - w.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + w.x[i][j_beg][k];
                    }
                }
            }
        }
    }
/*
* 
*/
    void IC_Nord_Polar_Meer (Array &h, Array &u, Array &v, Array &w) {
    // Currents along the coasts
    // Closing the polar, subpolar and subtropical atmospheric circulation systems
    // Arctic currents
    // Extension to the North Pole and widening (depth -45 m)
    // Throughflow of the Bering Strait (from j=0 until j=30 compares to 60°N until 90°N,
    //                                                           from k=180 until k=194 compares to 180°W until 166°W)
        #pragma omp parallel for collapse(3)
        for (int i = m.im-4; i < m.im; i++) {
            for (int j = 0; j < 31; j++) {
                for (int k = 180; k < 195; k++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = - IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = - IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }
        }


    // Arctic currents
    // Connection between Bering Sea and North Pacific
    // Extension until the North Pole and widening (depth -45 m)
    // Throughflow of the Bering Strait (from j=0 until j=30 compares to 60°N until 90°N,
    //                                                           from k=180 until k=194 compares to 180°W until 166°W)
        j_beg = 31;
        j_end = 36;
        k_beg = 180;
        k_end = 195;

        #pragma omp parallel for collapse(3) private(j_z, j_n)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    j_z = j - j_beg;
                    j_n = j_end - j_beg;

                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = (v.x[i][j_end][k] - v.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + v.x[i][j_beg][k];
                        w.x[i][j][k] = (w.x[i][j_end][k] - w.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + w.x[i][j_beg][k];
                    }
                }
            }
        }


    // Arctic currents
    // western connection between Bering Sea and North Pacific
    // Extension until the North Pole and widening (depth -45 m)
    // Throughflow of the Bering Strait (from j=0 until j=30 compares to 60°N until 90°N,
    //                                                           from k=180 until k=194 compares to 180°W until 166°W)
        j_beg = 0;
        j_end = 36;
        k_beg = 180;
        k_end = 195;

        #pragma omp parallel for collapse(3) private(k_z, k_n)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    k_z = k - k_beg;
                    k_n = k_end - k_beg;

                    if (h.x[i][j][k] == 0.) {
//                        u.x[i][j][k] = (u.x[i][j][k_end] - u.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + u.x[i][j][k_beg];
                        v.x[i][j][k] = (v.x[i][j][k_end] - v.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + v.x[i][j][k_beg];
                        w.x[i][j][k] = (w.x[i][j][k_end] - w.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + w.x[i][j][k_beg];
                    }
                }
            }
        }


    // Arctic currents
    // eastern connection between Bering Sea and North Pacific
    // Extension until the North Pole and widening (depth -45 m)
    // Throughflow of the Bering Strait (from j=0 until j=30 compares to 60°N until 90°N,
    //                                                           from k=180 until k=194 compares to 180°W until 166°W)
        j_beg = 0;
        j_end = 31;
        k_beg = 180;
        k_end = 195;

        #pragma omp parallel for collapse(3) private(k_z, k_n)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    k_z = k - k_beg;
                    k_n = k_end - k_beg;

                    if (h.x[i][j][k] == 0.) {
//                        u.x[i][j][k] = (u.x[i][j][k_end] - u.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + u.x[i][j][k_beg];
                        v.x[i][j][k] = (v.x[i][j][k_end] - v.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + v.x[i][j][k_beg];
                        w.x[i][j][k] = (w.x[i][j][k_end] - w.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + w.x[i][j][k_beg];
                    }
                }
            }
        }


    // Continuation of the Bering Strait current beyond the North Pole (depth -45 m)
    // Extension of the Bering Strait until the North Pole and widening  (from j=0 until j=20 compares to 90°N until 70°N,
    //                                                                                                                  from k=320 until k=m.km compares to 20°W until 0°)

        #pragma omp parallel for collapse(3)
        for (int i = m.im-4; i < m.im; i++) {
            for (int j = 0; j < 21; j++) {
                for (int k = 320; k < m.km; k++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = + IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = - IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }
        }
    }
/*
* 
*/
    void IC_South_Polar_Sea (Array &h, Array &u, Array &v, Array &w, Array &c) {
    // flow along coasts
    // closing the polar, subpolar and subtropic circulation systems
    // south polar sea
    // setting the South Pole with zero velocity
    // antarctic circumpolar current (-5000m deep) (from j=147 until j=152 compares to 57°S until 62°S,
    //                                                                            from k=0 until k=m.km compares to 0° until 360°)

        #pragma omp parallel for collapse(3)
        for (int i = 10; i < m.im-2; i++) {
            for (int j = 147; j < 153; j++) {
                for (int k = 0; k < m.km; k++) {
                    if (h.x[i][j][k] == 0.) {
                        c.x[i][j][k] = ca_max;
//                  v.x[i][j][k] = - IC_water * static_cast<double>(i - 5) / static_cast<double>(m.im-2 - 5);
//                  w.x[i][j][k] = IC_water * static_cast<double>(i - 5) / static_cast<double>(m.im-2 - 5);
//                  v.x[i][j][k] = - IC_water;
                        w.x[i][j][k] = IC_water;
                    }
                }
            }
        }


    // antarctic circumpolar current (-5000m deep)
    // north of the Antarctic Circumpolar Current
    // antarctic circumpolar current (from j=147 until j=152 compares to 57°S until 62°S,
    //                                                        from k=0 until k=m.km compares to 0° until 360°)
        j_beg = 144;
        j_end = 150;
        k_beg = 0;
        k_end = m.km;

        #pragma omp parallel for collapse(3) private(j_z, j_n)
        for (int i = 10; i < m.im - 2; i++) {
            for (int j = j_beg; j < j_end + 1; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    j_z = j - j_beg;
                    j_n = j_end - j_beg;

                    if (h.x[i][j][k] == 0.) {
                        c.x[i][j][k] = (c.x[i][j_end][k] - c.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + c.x[i][j_beg][k];
                        v.x[i][j][k] =  .5 * IC_water * (v.x[i][j_end][k] - v.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + v.x[i][j_beg][k];
//                        w.x[i][j][k] = .5 * IC_water * (w.x[i][j_end][k] - w.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + w.x[i][j_beg][k];
                    }
                }
            }
        }


    // antarctic circumpolar current (-5000m deep)
    // south of antarctic circumpolar current
    // antarctic circumpolar current (from j=147 until j=152 compares to 57°S until 62°S,
    //                                                 from k=0 until k=m.km compares to 0° until 360°)
        j_beg = 150;
        j_end = 155;
        k_beg = 0;
        k_end = m.km;

        #pragma omp parallel for collapse(3) private(j_z, j_n)
        for (int i = 10; i < m.im - 2; i++) {
            for (int j = j_beg; j < j_end + 1; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    j_z = j - j_beg;
                    j_n = j_end - j_beg;

                    if (h.x[i][j][k] == 0.) {
                        c.x[i][j][k] = (c.x[i][j_end][k] - c.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + c.x[i][j_beg][k];
                        v.x[i][j][k] = .5 * IC_water * (v.x[i][j_end][k] - v.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + v.x[i][j_beg][k];
                        w.x[i][j][k] = .5 * IC_water * (w.x[i][j_end][k] - w.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + w.x[i][j_beg][k];
                    }
                }
            }
        }


    // south polar sea
    // subantarctic  mode water (-1000m deep)
    // subantarctic front parallel to the circumpolar current (from j=144 until j=146 compares to 54°S until 56°S,
    //                                                                                        from k=0 until k=m.km compares to 0° until 360°)
        #pragma omp parallel for collapse(3)
        for (int i = i_beg+10; i < m.im-2; i++) {
            for (int j = 144; j < 147; j++) {
                for (int k = 0; k < m.km; k++) {
                    if (h.x[i][j][k] == 0.) {
                        c.x[i][j][k] = ca_max;
    //                  v.x[i][j][k] = - IC_water * static_cast<double>(i - i_beg + 5) / static_cast<double>(i_max - i_beg + 5);
                        w.x[i][j][k] = .5 * IC_water * static_cast<double>(i - i_beg + 5) / static_cast<double>(i_max - i_beg + 5);
                    }
                }
            }
        }


    // north of the constant velocity belt of the subantarctic mode water (-1000m deep)
    // subantarctic front parallel to the circumpolar current (from j=144 until j=146 compares to 54°S until 56°S,
    //                                                                                        from k=0 until k=m.km compares to 0° until 360°)
        j_beg = 144;
        j_end = 147;
        k_beg = 0;
        k_end = m.km;

        #pragma omp parallel for collapse(3) private(j_z, j_n)
        for (int i = i_beg+10; i < m.im-2; i++) {
            for (int j = j_beg; j < j_end + 1; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    j_z = j - j_beg;
                    j_n = j_end - j_beg;

                    if (h.x[i][j][k] == 0.) {
                        c.x[i][j][k] = ca_max;
                        v.x[i][j][k] = .5 * IC_water * (v.x[i][j_end][k] - v.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + v.x[i][j_beg][k];
                        w.x[i][j][k] = .5 * IC_water * (w.x[i][j_end][k] - w.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + w.x[i][j_beg][k];
                    }
                }
            }
        }


    // north of the constant velocity belt of the subantarctic mode water (-1000m deep)
    // subantarctic front parallel to the circumpolar current (from j=144 until j=146 compares to 54°S until 56°S,
    //                                                                                        from k=0 until k=m.km compares to 0° until 360°)
        j_beg = 144;
        j_end = 147;
        k_beg = 0;
        k_end = m.km;

        #pragma omp parallel for collapse(3) private(j_z, j_n)
        for (int i = i_beg+10; i < m.im-2; i++) {
            for (int j = j_beg; j < j_end + 1; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    j_z = j - j_beg;
                    j_n = j_end - j_beg;

                    if (h.x[i][j][k] == 0.) {
                        c.x[i][j][k] = ca_max;
                        v.x[i][j][k] = .5 * IC_water * (v.x[i][j_end][k] - v.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + v.x[i][j_beg][k];
                        w.x[i][j][k] = .5 * IC_water * (w.x[i][j_end][k] - w.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + w.x[i][j_beg][k];
                    }
                }
            }
        }


    // south polar sea
    // Weddell sea
    // diagonal from Lasarew sea until antarktic peninsula
    // south-west directed flow at the surface (from j=155 until j=165 compares to 65°S until 75°S,
    //                                                                  from k=299 until k=m.km compares to 61°W until 0°)
        j_beg = 153;
        j_end = 166;
        k_beg = 299;
        k_end = m.km;
        j_run = 0;
        k_run = 0;
        j_step = 1;
        k_step = 50;


        while ((j_end - j_run) >= j_beg && (k_beg + k_run) <=  k_end) {
            for (int j = std::min(j_end - j_run, m.jm - 1); j > (j_end - j_step - j_run); j--) {
                for (int k = (k_beg + k_run); k < (k_beg + k_step + k_run) && k < m.km && k >= 0; k++) {
                    for (int i = 10; i < m.im-2; i++) {
                        if (h.x[i][j][k] == 0.) {
                            c.x[i][j][k] = ca_max;
//                            v.x[i][j][k] = + IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                            w.x[i][j][k] = - IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        }
                    }
                }
            k_run = k_run+1;
            }
        j_run++;
        }


    // Weddell sea
    // north of Lasarew sea until antarktic peninsula
    // south-west directed flow at the surface (from j=155 until j=165 compares to 65°S until 75°S,
    //                                                                  from k=299 until k=m.km compares to 61°W until 0°)
        j_beg = 153;
        j_end = 166;
        k_beg = 299;
        k_end = m.km;

        #pragma omp parallel for collapse(3) private(j_z, j_n)
        for (int i = 10; i < m.im-2; i++) {
            for (int j = j_beg; j < j_end + 1; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    j_z = j - j_beg;
                    j_n = j_end - j_beg;

                    if (h.x[i][j][k] == 0.) {
                        c.x[i][j][k] = ca_max;
                        v.x[i][j][k] = .5 * IC_water * (v.x[i][j_end][k] - v.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + v.x[i][j_beg][k];
                        w.x[i][j][k] = .5 * IC_water * (w.x[i][j_end][k] - w.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + w.x[i][j_beg][k];
                    }
                }
            }
        }


    // Weddell sea
    // east of of Lasarew sea until antarktic peninsula
    // south-west directed flow at the surface (from j=155 until j=165 compares to 65°S until 75°S,
    //                                                                  from k=352 until k=m.km compares to 8°W until 0°)
        j_beg = 153;
        j_end = 166;
        k_beg = m.km - 8;
        k_end = m.km;

        #pragma omp parallel for collapse(3)
        for (int i = 10; i < m.im-2; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    if (h.x[i][j][k] == 0.) {
                        c.x[i][j][k] = ca_max;
                        v.x[i][j][k] = + .5 * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = + .5 * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }
        }

        j_beg = 153;
        j_end = 166;
        k_beg = m.km - 8;
        k_end = m.km;

        #pragma omp parallel for collapse(3)
        for (int i = 10; i < m.im-2; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    if (h.x[i][j][k] == 0.) {
                        c.x[i][j][k] = ca_max;
                        v.x[i][j][k] = + .5 * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = - .5 * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }
        }

        j_beg = 153;
        j_end = 166;
        k_beg = m.km - 43;
        k_end = m.km - 1;

        #pragma omp parallel for collapse(3) private(k_z, k_n)
        for (int i = 10; i < m.im-2; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    k_z = k - k_beg;
                    k_n = k_end - k_beg;

                    if (h.x[i][j][k] == 0.) {
                        c.x[i][j][k] = ca_max;
                        v.x[i][j][k] = .5 * IC_water * (v.x[i][j][k_end] - v.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + v.x[i][j][k_beg];
                        w.x[i][j][k] = .5 * IC_water * (w.x[i][j][k_end] - w.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + w.x[i][j][k_beg];
                    }
                }
            }
        }


    // south polar sea
    // Ross sea
    // diagonally from Martin peninsula sea until Kap Adare
    // south-west directed flow at the surface (from j=155 until j=168 compares to 65°S until 78°S,
    //                                                                  from k=180 until k=280 compares to 180°W until 80°W)
        j_beg = 153;
        j_end = 169;
        k_beg = 180;
        k_end = 281;
        j_run = 0;
        k_run = 0;
        j_step = 1;
        k_step = 90;

        while ((j_end - j_run) >= j_beg && (k_beg + k_run) <= k_end) {
            for (int j = std::min(j_end - j_run, m.jm - 1); j > (j_end - j_step - j_run); j--) {
                for (int k = (k_beg + k_run); k < (k_beg + k_step + k_run) && k < m.km && k >= 0; k++) {
                    for (int i = 10; i < m.im-2; i++) {
                        if (h.x[i][j][k] == 0.) {
                            c.x[i][j][k] = ca_max;
    //                      v.x[i][j][k] = + IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                            w.x[i][j][k] = - IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        }
                    }
                }
            k_run = k_run+1;
            }
        j_run++;
        }


    // Ross sea
    // north of Martin peninsula sea until Kap Adare
    // south-west directed flow at the surface (from j=155 until j=168 compares to 65°S until 78°S,
    //                                                                  from k=180 until k=245 compares to 180°W until 115°W)
        j_beg = 153;
        j_end = 169;
        k_beg = 180;
        k_end = 246;

        #pragma omp parallel for collapse(3) private(j_z, j_n)
        for (int i = 10; i < m.im-2; i++) {
            for (int j = j_beg; j < j_end + 1; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    j_z = j - j_beg;
                    j_n = j_end - j_beg;

                    if (h.x[i][j][k] == 0.) {
                        c.x[i][j][k] = ca_max;
                        v.x[i][j][k] = .5 * IC_water * (v.x[i][j_end][k] - v.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + v.x[i][j_beg][k];
                        w.x[i][j][k] = .5 * IC_water * (w.x[i][j_end][k] - w.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + w.x[i][j_beg][k];
                    }
                }
            }
        }


    // Ross sea
    // east of Martin peninsula sea until Kap Adare
    // south-west directed flow at the surface (from j=155 until j=168 compares to 65°S until 78°S,
    //                                                                  from k=180 until k=245 compares to 180°W until 115°W)
        j_beg = 153;
        j_end = 169;
        k_beg = 254;
        k_end = 262;

        #pragma omp parallel for collapse(3)
        for (int i = 10; i < m.im-2; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    if (h.x[i][j][k] == 0.) {
                        c.x[i][j][k] = ca_max;
                        v.x[i][j][k] = + .5 * IC_water * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = + .5 * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }
        }

        j_beg = 153;
        j_end = 169;
        k_beg = 254;
        k_end = 262;

        #pragma omp parallel for collapse(3)
        for (int i = 10; i < m.im-2; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    if (h.x[i][j][k] == 0.) {
                        c.x[i][j][k] = ca_max;
                        v.x[i][j][k] = + .5 * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        w.x[i][j][k] = - .5 * IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                    }
                }
            }
        }

        j_beg = 153;
        j_end = 169;
        k_beg = 190;
        k_end = 250;

        #pragma omp parallel for collapse(3) private(k_z, k_n)
        for (int i = 10; i < m.im-2; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    k_z = k - k_beg;
                    k_n = k_end - k_beg;

                    if (h.x[i][j][k] == 0.) {
                        c.x[i][j][k] = ca_max;
                        v.x[i][j][k] = .5 * IC_water * (v.x[i][j][k_end] - v.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + v.x[i][j][k_beg];
                        w.x[i][j][k] = .5 * IC_water * (w.x[i][j][k_end] - w.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + w.x[i][j][k_beg];
                    }
                }
            }
        }

        j_beg = 153;
        j_end = 169;
        k_beg = 246;
        k_end = 255;

        #pragma omp parallel for collapse(3) private(k_z, k_n)
        for (int i = 10; i < m.im-2; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    k_z = k - k_beg;
                    k_n = k_end - k_beg;

                    if (h.x[i][j][k] == 0.) {
                        c.x[i][j][k] = ca_max;
                        v.x[i][j][k] = .5 * IC_water * (v.x[i][j][k_end] - v.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + v.x[i][j][k_beg];
                        w.x[i][j][k] = .5 * IC_water * (w.x[i][j][k_end] - w.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + w.x[i][j][k_beg];
                    }
                }
            }
        }


    // south polar sea
    // antarctic Indic basin
    // diagonally from Adelieland until Queen Mary Land
    // south-west directed flow at the surface (from j=152 until j=158 compares to 62°S until 68°S,
    //                                                                  from k=0 until k=160 compares to 0°O until 160°O)
        j_beg = 152;
        j_end = 159;
        k_beg = 0;
        k_end = 161;
        j_run = 0;
        k_run = 0;
        j_step = 1;
        k_step = 110;

        while ((j_end - j_run) >= j_beg && (k_beg + k_run) <= k_end) {
            for (int j = std::min(j_end - j_run, m.jm - 1); j > (j_end - j_step - j_run); j--) {
                for (int k = (k_beg + k_run); k < (k_beg + k_step + k_run) && k < m.km && k >= 0; k++) {
                    for (int i = 10; i < m.im-2; i++) {
                        if (h.x[i][j][k] == 0.) {
                            c.x[i][j][k] = ca_max;
//                            v.x[i][j][k] = + IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                            w.x[i][j][k] = - IC_water * static_cast<double>(i - i_beg) / static_cast<double>(i_max - i_beg);
                        }
                    }
                }
            k_run = k_run+1;
            }
        j_run++;
        }


    // antarctic Indic basin
    // north of Adelieland until Queen Mary Land
    // south-west directed flow at the surface (from j=152 until j=158 compares to 62°S until 68°S,
    //                                                                  from k=0 until k=160 compares to 0°O until 160°O)
        j_beg = 152;
        j_end = 159;
        k_beg = 0;
        k_end = 161;

        #pragma omp parallel for collapse(3) private(j_z, j_n)
        for (int i = 10; i < m.im-2; i++) {
            for (int j = j_beg; j < j_end + 1; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    j_z = j - j_beg;
                    j_n = j_end - j_beg;

                    if (h.x[i][j][k] == 0.) {
                        c.x[i][j][k] = ca_max;
                        v.x[i][j][k] = .5 * IC_water * (v.x[i][j_end][k] - v.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + v.x[i][j_beg][k];
                        w.x[i][j][k] = .5 * IC_water * (w.x[i][j_end][k] - w.x[i][j_beg][k]) * static_cast<double>(j_z) / static_cast<double>(j_n) + w.x[i][j_beg][k];
                    }
                }
            }
        }

        j_beg = 152;
        j_end = 159;
        k_beg = 132;
        k_end = 138;

        #pragma omp parallel for collapse(3) private(k_z, k_n)
        for (int i = 10; i < m.im-2; i++) {
            for (int j = j_beg; j < j_end; j++) {
                for (int k = k_beg; k < k_end; k++) {
                    k_z = k - k_beg;
                    k_n = k_end - k_beg;

                    if (h.x[i][j][k] == 0.) {
                        c.x[i][j][k] = ca_max;
                        v.x[i][j][k] = .5 * IC_water * (v.x[i][j][k_end] - v.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + v.x[i][j][k_beg];
                        w.x[i][j][k] = .5 * IC_water * (w.x[i][j][k_end] - w.x[i][j][k_beg]) * static_cast<double>(k_z) / static_cast<double>(k_n) + w.x[i][j][k_beg];
                    }
                }
            }
        }

    }
/*
* 
*/
    void IC_EquatorialCurrents (Array &h, Array &u, Array &v, Array &w) {
    // currents along the equator
    // equatorial undercurrent - Cromwell flow, EUC
    // upwelling at the end of the equatorial undercurrent - Cromwell flow, EUC
    // equatorial intermediate current, EIC
    // nothern and southern equatorial subsurface counter-currents, NSCC und SSCC
    // nothern and southern equatorial counter-currents, NECC und SECC
    // domes at the end of nothern and southern equatorial counter-currents, NECC und SECC

    // for the depth of the sea compares i=1 to a depth of 150 m    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        i_EIC_u = m.im - 14;
        i_EIC_o = m.im - 6;
        i_SCC_u = m.im - 10;
        i_SCC_o = m.im - 6;
        i_ECC_u = m.im - 6;

        i_ECC_o = m.im;



    // equatorial currents and counter-currents

    //    §§§§§§§§§§§§§§§§§§§§§§§§§     Pacific ocean       §§§§§§§§§§§§§§§§§§§§§§§§§§

    // Pacific ocean
    // equatorial undercurrent - Cromwell current (EUC, i=m.im-2 until i=m.im-1 compares to -100 until -200m depth)
    // equatorial undercurrent - Cromwell current (from j=87 until j=93 compares to 3°N until 3°S,
    //                                                                         from k=145 until k=270 compares to 145°O until 90°W)
            #pragma omp parallel for collapse(3)
        for (int i = m.im-4; i < m.im-2; i++) {
            for (int j = 87; j < 94; j++) {
                for (int k = 145; k < 271; k++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = 0.;
                        w.x[i][j][k] = + .05 * IC_water;
                    }
                }
            }
        }


    // Pacific ocean
    // downward flow at the end of the equatorial undercurrent - Cromwell current (EUC, -200m, 3°N - 3°S, 85°W - 90°W)

        #pragma omp parallel for collapse(3)
        for (int i = m.im-14; i < m.im-2; i++) {
            for (int j = 87; j < 94; j++) {
                for (int k = 265; k < 271; k++) {
                    if (h.x[i][j][k] == 0.) {
                        {
                            u.x[i][j][k] = + .01 * IC_water;
                            w.x[i][j][k] = - .05 * IC_water;
                        }
                    }
                }
            }
        }


    // Pacific ocean
    // upward flow at the end of the equatorial undercurrent - Cromwell current (EUC, -200m, 3°N - 3°S, 140°O - 145°O)

        #pragma omp parallel for collapse(3)
        for (int i = m.im-14; i < m.im-2; i++) {
            for (int j = 87; j < 94; j++) {
                for (int k = 140; k < 146; k++) {
                    if (h.x[i][j][k] == 0.) {
                        {
                            u.x[i][j][k] = - .01 *IC_water;
                            w.x[i][j][k] = + .05 * IC_water;
                        }
                    }
                }
            }
        }


    // Pacific ocean
    // equatorial intermediate current (EIC, i=m.im-4 until i=m.im-2 compares to -300 until -1000m depth)
    // equatorial intermediate current (from j=88 until j=92 compares to 2°N until 2°S,
    //                                                      from k=145 until k=270 compares to 145°O until 90°W)
        #pragma omp parallel for collapse(3)
        for (int i = i_EIC_u; i < i_EIC_o; i++) {
            for (int j = 88; j < 93; j++) {
                for (int k = 145; k < 271; k++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = 0.;
                        w.x[i][j][k] = - .05 * IC_water * static_cast<double>(i - i_EIC_u) / static_cast<double>(i_EIC_o - i_EIC_u);
                    }
                }
            }
        }


    // Pacific ocean
    // equatorial northern and southern subsurface counter-current
    // equatorial northern subsurface counter-current (NSCC, i=m.im-3 until i=m.im-2 compares to -300 until -800m depth)
    // equatorial northern subsurface counter-current (from j=86 until j=87 compares to 3°N until 4°N,
    //                                                                               from k=135 until k=270 compares to 135°O until 90°W)
        #pragma omp parallel for collapse(3)
        for (int i = i_SCC_u; i < i_SCC_o; i++) {
            for (int j = 86; j < 88; j++) {
                for (int k = 135; k < 271; k++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = 0.;
                        w.x[i][j][k] = .05 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_SCC_o - i_SCC_u);
                    }
                }
            }
        }


    // Pacific ocean
    // equatorial southern subsurface counter-current (SSCC, i=m.im-3 until i=m.im-2 compares to -300 until -800m depth)
    // equatorial southern subsurface counter-current (from j=93 until j=94 compares to 3°S until 4°S,
    //                                                                               from k=165 until k=270 compares to 165°O until 90°W)
        #pragma omp parallel for collapse(3)
        for (int i = i_SCC_u; i < i_SCC_o; i++) {
            for (int j = 93; j < 95; j++) {
                for (int k = 165; k < 271; k++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = 0.;
                        w.x[i][j][k] = .05 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_SCC_o - i_SCC_u);
                    }
                }
            }
        }


    // Pacific ocean
    // equatorial current at the surface
    // equatorial northern counter-current (NECC, i=m.im-1 until i=m.im compares to 0 until -200m depth)
    // equatorial northern counter-current (from j=83 until j=87 compares to 3°N until 7°N,
    //                                                            from k=135 until k=270 compares to 135°O until 90°W)
        #pragma omp parallel for collapse(3)
        for (int i = i_ECC_u; i < i_ECC_o; i++) {
            for (int j =83; j < 88; j++) {
                for (int k = 135; k < 271; k++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = 0.;
                        w.x[i][j][k] = + .05 * IC_water * static_cast<double>(i - i_ECC_u) / static_cast<double>(i_ECC_o - i_ECC_u);
                    }
                }
            }
        }


    // Pacific ocean
    // equatorial northern and southern counter-currents at the surface
    // equatorial currents between NECC and SECC
    // equatorial northern counter-current (from j=87 until j=93 compares to 3°N until 3°S,
    // equatorial northern counter-current   from k=135 until k=270 compares to 135°O until 90°W)
        #pragma omp parallel for collapse(3)
        for (int i = m.im-2; i < m.im; i++) {
            for (int j = 87; j < 94; j++) {
                for (int k = 135; k < 271; k++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = 0.;
                        w.x[i][j][k] = - .05 * IC_water;
                    }
                }
            }
        }


    // Pacific ocean
    // equatorial northern and southern counter-currents at the surface
    // equatorial currents between NECC and SECC
    // equatorial southern counter-current (SECC, i=m.im-1 until i=m.im compares to 0 until -200m depth)
    // equatorial southern counter-current (from j=93 until j=96 compares to 3°S until 6°S,
    //                                                             from k=155 until k=270 compares to 155°O until 90°W)
        #pragma omp parallel for collapse(3)
        for (int i = i_ECC_u; i < i_ECC_o; i++) {
            for (int j = 93; j < 97; j++) {
                for (int k = 155; k < 271; k++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = 0.;
                        w.x[i][j][k] = + .05 * IC_water * static_cast<double>(i - i_ECC_u) / static_cast<double>(i_ECC_o - i_ECC_u);
                    }
                }
            }
        }


    // Pacific ocean
    // Guatemala dome at the end of the northern equatorial counter-current at the surface, Guinea current = NECC
    // Guatemala-dome (belonging to NECC, i=m.im-1 until i=m.im compares to 0 until -200m depth)
    // Guatemala-dome (from j=83 until j=87 compares to 3°N until 7°N,
    //                                from k=265 until k=275 compares to 85°W until 95°W)
        #pragma omp parallel for collapse(3)
        for (int i = i_SCC_u; i < i_ECC_o - 1; i++) {
            for (int j = 83; j < 88; j++) {
                for (int k = 265; k < 276; k++) {
                    if (h.x[i][j][k] == 0.) {
                        {
                            u.x[i][j][k] = + .01 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_ECC_o - 1 - i_SCC_u);
                            v.x[i][j][k] = - .05 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_ECC_o - 1 - i_SCC_u);
//                            w.x[i][j][k] = - .05 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_ECC_o - 1 - i_SCC_u);
                        }
                    }
                }
            }
        }


    // Pacific ocean
    // Equador dome at the end of the southern equatorial counter-current at the surface
    // Equador-dome (belonging to SECC, i=m.im-1 until i=m.im compares to 0 until -200m depth)
    // Equador-dome (from j=93 until j=96 compares to 3°S until 6°S,
    //                           from k=265 until k=275 compares to 85°W until 95°W)
        #pragma omp parallel for collapse(3)
        for (int i = i_SCC_u; i < i_ECC_o - 1; i++) {
            for (int j = 93; j < 97; j++) {
                for (int k = 265; k < 276; k++) {
                    if (h.x[i][j][k] == 0.) {
                        {
                            u.x[i][j][k] = + .01 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_ECC_o - 1 - i_SCC_u);
                            v.x[i][j][k] = + .05 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_ECC_o - 1 - i_SCC_u);
//                            w.x[i][j][k] = - .05 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_ECC_o - 1 - i_SCC_u);
                        }
                    }
                }
            }
        }


    // Pacific ocean
    // downward flow opposit the Guatemala dome at the end of the northern equatorial counter-current at the surface, Guinea current = NECC
    // Guatemala-dome (from j=83 until j=87 compares to 3°N until 7°N,
    //                                from k=135 until k=140 compares to 135°W until 140°W)
        #pragma omp parallel for collapse(3)
        for (int i = i_SCC_u; i < i_ECC_o - 1; i++) {
            for (int j = 83; j < 88; j++) {
                for (int k = 135; k < 141; k++) {
                    if (h.x[i][j][k] == 0.) {
                        {
                            u.x[i][j][k] = - .01 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_ECC_o - 1 - i_SCC_u);
                            v.x[i][j][k] = + .05 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_ECC_o - 1 - i_SCC_u);
//                            w.x[i][j][k] = - .05 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_ECC_o - 1 - i_SCC_u);
                        }
                    }
                }
            }
        }


    // Pacific ocean
    // downward flow opposit the Equador-dome dome at the end of the northern equatorial counter-current at the surface, Guinea current = NECC
    // Equador-dome (from j=93 until j=96 compares to 3°S until 6°S,
    //                            from k=155 until k=160 compares to 155°O until 160°O)
        #pragma omp parallel for collapse(3)
        for (int i = i_SCC_u; i < i_ECC_o - 1; i++) {
            for (int j = 93; j < 97; j++) {
                for (int k = 155; k < 161; k++) {
                    if (h.x[i][j][k] == 0.) {
                        {
                            u.x[i][j][k] = - .01 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_ECC_o - 1 - i_SCC_u);
                            v.x[i][j][k] = - .05 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_ECC_o - 1 - i_SCC_u);
//                            w.x[i][j][k] = - .05 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_ECC_o - 1 - i_SCC_u);
                        }
                    }
                }
            }
        }



    //    §§§§§§§§§§§§§§§§§§§§§§§§§     Indic ocean     §§§§§§§§§§§§§§§§§§§§§§§§§§


    // equator currents and counter-currents
    // Indic ocean
    // equatorial under current - Cromwell current (EUC, i=m.im-2 until i=m.im-1 compares to -100 until -200m depth)
    // equatorial under current - Cromwell current (from j=89 until j=91 compares to 1°N until 1°S,
    //                                                                             from k=55 until k=90 compares to 55°O until 90°O)
        #pragma omp parallel for collapse(3)
        for (int i = m.im-2; i < m.im-1; i++) {
            for (int j = 89; j < 92; j++) {
                for (int k = 55; k < 91; k++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = 0.;
                        w.x[i][j][k] = .05 * IC_water;
                    }
                }
            }
        }


    // Indic ocean
    // downward flow at the end of the equatorial under-current - Cromwell-current (EUC, -200m, 3°N - 3°S, 85°W - 90°W)
    // equatorial under current - Cromwell current (from j=89 until j=91 compares to 1°N until 1°S,
    //                                                                              from k=85 until k=90 compares to 85°O until 90°O)
        #pragma omp parallel for collapse(3)
        for (int i = m.im-14; i < m.im-2; i++) {
            for (int j = 89; j < 92; j++) {
                for (int k = 84; k < 91; k++) {
                    if (h.x[i][j][k] == 0.) {
                        {
                            u.x[i][j][k] = + .01 * IC_water;
                            w.x[i][j][k] = - .05 * IC_water;
                        }
                    }
                }
            }
        }


    // Indic ocean
    // upward flow at the end of the equatorial under-current - Cromwell-current (EUC, -200m, 3°N - 3°S, 140°W - 145°W)
    // equatorial under current - Cromwell current (from j=89 until j=91 compares to 1°N until 1°S,
    //                                                                             from k=50 until k=55 compares to 50°O until 55°O)
        #pragma omp parallel for collapse(3)
        for (int i = m.im-1; i < m.im-2; i++) {
            for (int j = 89; j < 92; j++) {
                for (int k = 50; k < 56; k++) {
                    if (h.x[i][j][k] == 0.) {
                        {
                            u.x[i][j][k] = - .01 * IC_water;
                            w.x[i][j][k] = + .05 * IC_water;
                        }
                    }
                }
            }
        }


    // Indic ocean
    // equatorial intermediate current (EIC, i=m.im-4 until i=m.im-2 compares to -300 until -1000m depth)
    // equatorial intermediate current (from j=89 until j=91 compares to 1°N until 1°S,
    //                                                      from k=55 until k=90 compares to 55°O until 90°O)
        #pragma omp parallel for collapse(3)
        for (int i = i_EIC_u; i < i_EIC_o; i++) {
            for (int j = 89; j < 92; j++) {
                for (int k = 55; k < 91; k++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = 0.;
                        w.x[i][j][k] = - .05 * IC_water * static_cast<double>(i - i_EIC_u) / static_cast<double>(i_EIC_o - i_EIC_u);
                    }
                }
            }
        }

    // Indic ocean
    // northern and southern equatorial subsurface counter-currents
    // equatorial northern subsurface counter-current (NSCC, i=m.im-3 until i=m.im-2 compares to -300 until -800m depth)
    // equatorial northern subsurface counter-current (from j=85 until j=88 compares to 2°N until 5°N,
    //                                                                         from k=55 until k=90 compares to 55°O until 90°O)
        #pragma omp parallel for collapse(3)
        for (int i = i_SCC_u; i < i_SCC_o; i++) {
            for (int j = 85; j < 89; j++) {
                for (int k = 55; k < 91; k++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = 0.;
                        w.x[i][j][k] = .05 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_SCC_o - i_SCC_u);
                    }
                }
            }
        }


    // Indic ocean
    // equatorial northern subsurface counter-current (SSCC, i=m.im-3 until i=m.im-2 compares to -300 until -800m depth)
    // equatorial northern subsurface counter-current (from j=92 until j=95 compares to 2°S until 5°S,
    //                                                                        from k=55 until k=90 compares to 55°O until 90°O)
        #pragma omp parallel for collapse(3)
        for (int i = i_SCC_u; i < i_SCC_o; i++) {
            for (int j = 92; j < 96; j++) {
                for (int k = 55; k < 91; k++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = 0.;
                        w.x[i][j][k] = .05 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_SCC_o - i_SCC_u);
                    }
                }
            }
        }


    // Indic ocean
    // equatorial northern counter-current at the surface, SW Monsun-current = NECC
    // equatorial northern counter-current (NECC, i=m.im-1 until i=m.im compares to 0 until -200m depth)
    // equatorial northern counter-current (from j=85 until j=88 compares to 2°N until 5°N,
    //                                                                 from k=55 until k=90 compares to 55°O until 90°O)
        #pragma omp parallel for collapse(3)
        for (int i = i_ECC_u; i < i_ECC_o; i++) {
            for (int j = 85; j < 99; j++) {
                for (int k = 55; k < 91; k++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = 0.;
                        w.x[i][j][k] = + .05 * IC_water * static_cast<double>(i - i_ECC_u) / static_cast<double>(i_ECC_o - i_ECC_u);
                    }
                }
            }
        }


    // Indic ocean
    // Northern and southern equatorial counter-current at the surface
    // Equatorial flow between NECC and SECC
    // Equatorial counter-currents (from j=89 until j=91 compares to 1°N until 1°S,
    //                                                     from k=55 until k=90 compares to 55°O until 90°O)
        #pragma omp parallel for collapse(3)
        for (int i = m.im-2; i < m.im; i++) {
            for (int j = 89; j < 92; j++) {
                for (int k = 55; k < 91; k++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = 0.;
                        w.x[i][j][k] = - .05 * IC_water;
                    }
                }
            }
        }


    // Indic ocean
    // Southern equatorial counter-current at the surface
    // equatorial southern counter-current (SECC, i=m.im-1 until i=m.im compares to 0 until -200m depth)
    // equatorial southern counter-current (from j=92 until j=95 compares to 2°S until 5°S,
    //                                                                from k=55 until k=90 compares to 55°O until 90°O)
        #pragma omp parallel for collapse(3)
        for (int i = i_ECC_u; i < i_ECC_o; i++) {
            for (int j = 92; j < 95; j++) {
                for (int k = 55; k < 91; k++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = 0.;
                        w.x[i][j][k] = + .05 * IC_water * static_cast<double>(i - i_ECC_u) / static_cast<double>(i_ECC_o - i_ECC_u);
                    }
                }
            }
        }


    // Indic ocean
    // North Sumatra dome at the end of the northern equatorial counter-current at the surface, Guinea current = NECC
    // North Sumatra dome (belonging to NECC, i=m.im-1 until i=m.im compares to 0 until -200m depth)
    // North Sumatra dome (from j=85 until j=88 compares to 2°N until 5°N,
    //                                  from k=88 until k=95 compares to 88°O until 95°O)
        #pragma omp parallel for collapse(3)
        for (int i = i_SCC_u; i < i_ECC_o - 1; i++) {
            for (int j = 85; j < 89; j++) {
                for (int k = 88; k < 96; k++) {
                    if (h.x[i][j][k] == 0.) {
                        {
                            u.x[i][j][k] = + .01 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_ECC_o - 1 - i_SCC_u);
                            v.x[i][j][k] = - .05 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_ECC_o - 1 - i_SCC_u);
//                            w.x[i][j][k] = .05 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_ECC_o - 1 - i_SCC_u);
                        }
                    }
                }
            }
        }


    // Indic ocean
    // South Sumatra dome at the end of the southern equatorial counter-current at the surface
    // South Sumatra dome (belonging to SECC, i=m.im-1 until i=m.im compares to 0 until -200m depth)
    // South Sumatra dome (from j=92 until j=95 compares to 2°S until 5°S,
    //                                 from k=90 until k=98 compares to 90°O until 98°O)
        #pragma omp parallel for collapse(3)
        for (int i = i_SCC_u; i < i_ECC_o - 1; i++) {
            for (int j = 92; j < 96; j++) {
                for (int k = 90; k < 99; k++) {
                    if (h.x[i][j][k] == 0.) {
                        {
                            u.x[i][j][k] = + .01 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_ECC_o - 1 - i_SCC_u);
                            v.x[i][j][k] = + .05 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_ECC_o - 1 - i_SCC_u);
//                            w.x[i][j][k] = + .05 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_ECC_o - 1 - i_SCC_u);
                        }
                    }
                }
            }
        }


    // Indic ocea
    // Outflow opposite the North Sumatra dome at the end of the northern equatorial counter-current at the surface, Guinea current = NECC
    // Outflow opposite the North Sumatra dome (belonging to NECC, i=m.im-1 until i=m.im compares to 0 until -200m depth)
    // Outflow opposite the North Sumatra dome (from j=84 until j=87 compares to 3°N until 6°N,
    //                                                                                 from k=55 until k=59 compares to 55°O until 59°O)
        #pragma omp parallel for collapse(3)
        for (int i = i_SCC_u; i < i_ECC_o - 1; i++) {
            for (int j = 84; j < 88; j++) {
                for (int k = 55; k < 60; k++) {
                    if (h.x[i][j][k] == 0.) {
                        {
                            u.x[i][j][k] = - .01 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_ECC_o - 1 - i_SCC_u);
                            v.x[i][j][k] = + .05 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_ECC_o - 1 - i_SCC_u);
//                            w.x[i][j][k] = .05 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_ECC_o - 1 - i_SCC_u);
                        }
                    }
                }
            }
        }


    // Indic ocean
    // Outflow opposite the South Sumatra dome at the end of the southern equatorial counter-current at the surface
    // Outflow opposite the South Sumatra dome (belonging to SECC, i=m.im-1 until i=m.im compares to 0 until -200m depth)
    // Outflow opposite the South Sumatra dome (from j=93 until j=96 compares to 3°S until 6°S,
    //                                                                                from k=55 until k=59 compares to 55°O until 59°O)
        #pragma omp parallel for collapse(3)
        for (int i = i_SCC_u; i < i_ECC_o - 1; i++) {
            for (int j = 93; j < 97; j++) {
                for (int k = 55; k < 60; k++) {
                    if (h.x[i][j][k] == 0.) {
                        {
                            u.x[i][j][k] = - .01 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_ECC_o - 1 - i_SCC_u);
                            v.x[i][j][k] = - .05 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_ECC_o - 1 - i_SCC_u);
//                            w.x[i][j][k] = + .05 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_ECC_o - 1 - i_SCC_u);
                        }
                    }
                }
            }
        }



    //    §§§§§§§§§§§§§§§§§§§§§§§§§     Atlantic ocean      §§§§§§§§§§§§§§§§§§§§§§§§§§

    // Atlantic ocean
    // equator currents and counter-currents
    // equatorial under current - Cromwell current (EUC, i=m.im-2 until i=m.im-1 compares to -100 until -200m depth)
    // equatorial under current - Cromwell current (from j=89 until j=91 compares to 1°N until 1°S,
    //                                                                             from k=326 until k=355 compares to 34°W until 5°W)
        #pragma omp parallel for collapse(3)
        for (int i = m.im-4; i < m.im-2; i++) {
            for (int j = 89; j < 92; j++) {
                for (int k = 326; k < 356; k++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = 0.;
                        w.x[i][j][k] = .05 * IC_water;
                    }
                }
            }
        }


    // Atlantic ocean
    // upward flow at the end of the equatorial under-current - Cromwell-current (EUC, -200m, 3°N - 3°S, 30°W - 35°W)
    // equatorial under current - Cromwell current (from j=89 until j=91 compares to 1°N until 1°S,
    //                                                                             from k=326 until k=330 compares to 34°W until 30°W)
        #pragma  omp parallel for collapse(3)
        for (int i = m.im-14; i < m.im-2; i++) {
            for (int j = 89; j < 92; j++) {
                for (int k = 326; k < 330; k++) {
                    if (h.x[i][j][k] == 0.) {
                        {
                            u.x[i][j][k] = - .01 * IC_water;
                            w.x[i][j][k] = + .05 * IC_water;
                        }
                    }
                }
            }
        }


    // Atlantic ocean
    // downward flow at the end of the equatorial under-current - Cromwell-current (EUC, -200m, 2°N - 0°, 3°W - 7°W)
    // equatorial under current - Cromwell current (from j=89 until j=91 compares to 1°N until 1°S,
    //                                                                             from k=350 until k=355 compares to 10°W until 5°W)
        #pragma  omp parallel for collapse(3)
        for (int i = m.im-14; i < m.im-2; i++) {
            for (int j = 89; j < 92; j++) {
                for (int k = 350; k < 356; k++) {
                    if (h.x[i][j][k] == 0.) {
                        {
                            u.x[i][j][k] = + .01 * IC_water;
                            w.x[i][j][k] = - .05 * IC_water;
                        }
                    }
                }
            }
        }


    // Atlantic ocean
    // equatorial intermediate current (EIC, i=m.im-4 until i=m.im-2 compares to -300 until -1000m depth,)
    // equatorial intermediate current (from j=88 until j=92 compares to 2°N until 2°S,
    //                                                      from k=330 until k=355 compares to 30°W until 5°W)
        #pragma  omp parallel for collapse(3)
        for (int i = i_EIC_u; i < i_EIC_o; i++) {
            for (int j = 88; j < 93; j++) {
                for (int k = 330; k < 356; k++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = 0.;
                        w.x[i][j][k] = - .05 * IC_water * static_cast<double>(i - i_EIC_u) / static_cast<double>(i_EIC_o - i_EIC_u);
                    }
                }
            }
        }


    // Atlantic ocean
    // northern and southern equatorial subsurface counter-currents
    // equatorial northern subsurface counter-current (NSCC, i=m.im-3 until i=m.im-2 compares to -300 until -800m depth)
    // equatorial northern subsurface counter-current (from j=86 until j=87 compares to 3°N until 4°N,
    //                                                                         from k=320 until k=340 compares to 40°W until 20°W)
        #pragma  omp parallel for collapse(3)
        for (int i = i_SCC_u; i < i_SCC_o; i++) {
            for (int j = 86; j < 88; j++) {
                for (int k = 320; k < 341; k++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = 0.;
                        w.x[i][j][k] = .05 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_SCC_o - i_SCC_u);
                    }
                }
            }
        }


    // Atlantic ocean
    // equatorial northern subsurface counter-current (SSCC, i=m.im-3 until i=m.im-2 compares to -300 until -800m depth)
    // equatorial northern subsurface counter-current (from j=93 until j=94 compares to 3°S until 4°S,
    //                                                                        from k=329 until k=355 compares to 31°W until 5°W)
        #pragma  omp parallel for collapse(3)
        for (int i = i_SCC_u; i < i_SCC_o; i++) {
            for (int j = 93; j < 95; j++) {
                for (int k = 329; k < 355; k++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = 0.;
                        w.x[i][j][k] = .05 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_SCC_o - i_SCC_u);
                    }
                }
            }
        }


    // Atlantic ocean
    // equatorial currents at the surface
    // equatorial northern counter-current (NECC, i=m.im-1 until i=m.im compares to 0 until -200m depth)
    // equatorial northern counter-current (from j=84 until j=87 compares to 3°N until 6°N,
    //                                                                 from k=330 until k=340 compares to 30°W until 20°W)
        #pragma  omp parallel for collapse(3)
        for (int i = i_ECC_u; i < i_ECC_o; i++) {
            for (int j = 84; j < 88; j++) {
                for (int k = 330; k < 340; k++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = 0.;
                        w.x[i][j][k] = + .05 * IC_water * static_cast<double>(i - i_ECC_u) / static_cast<double>(i_ECC_o - i_ECC_u);
                    }
                }
            }
        }


    // Atlantic ocean
    // northern and southern equatorial currents at the surface
    // Equatorial flow between NECC and SECC
    // equatorial northern counter-current (from j=87 until j=93 compares to 3°N until 3°S,
    //                                                                 from k=330 until k=355 compares to 30°W until 5°W)
        #pragma  omp parallel for collapse(3)
        for (int i = m.im-2; i < m.im; i++) {
            for (int j = 87; j < 94; j++) {
                for (int k = 330; k < 356; k++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = 0.;
                        w.x[i][j][k] = - .05 * IC_water;
                    }
                }
            }
        }


    // Atlantic ocean
    // northern and southern equatorial currents at the surface
    // equatorial southern counter-current (SECC, i=m.im-1 until i=m.im compares to 0 until -200m depth)
    // equatorial southern counter-current (from j=83 until j=85 compares to 3°S until 6°S,
    //                                                               from k=293 until k=316 compares to 30°W until 5°W)
        #pragma  omp parallel for collapse(3)
        for (int i = i_ECC_u; i < i_ECC_o; i++) {
            for (int j = 83; j < 86; j++) {
                for (int k = 293; k < 317; k++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = 0.;
                        w.x[i][j][k] = + .05 * IC_water * static_cast<double>(i - i_ECC_u) / static_cast<double>(i_ECC_o - i_ECC_u);
                    }
                }
            }
        }


    // Atlantic ocean
    // Guinea-dome at the end of the northern equatorial counter-current at the surface, Guinea-current = NECC
    // Guinea-dome (zum NECC gehörend, i=m.im-1 until i=m.im compares to 0 until -200m depth)
    // Guinea-dome (from j=84 until j=87 compares to 3°N until 6°N,
    //                         from k=340 until k=345 compares to 15°W until 20°W)
        #pragma  omp parallel for collapse(3)
        for (int i = i_SCC_u; i < i_ECC_o - 1; i++) {
            for (int j = 84; j < 88; j++) {
                for (int k = 340; k < 345; k++) {
                    if (h.x[i][j][k] == 0.) {
                        {
                            u.x[i][j][k] = + .01 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_ECC_o - 1 - i_SCC_u);
                            v.x[i][j][k] = - .05 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_ECC_o - 1 - i_SCC_u);
//                            w.x[i][j][k] = + .05 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_ECC_o - 1 - i_SCC_u);
                        }
                    }
                }
            }
        }


    // Atlantic ocean
    // Angola-dome at the end of the southern equatorial counter-current at the surface
    // Angola-dome (zum SECC gehörend, i=m.im-1 until i=m.im compares to 0 until -200m depth)
    // Angola-dome (from j=93 until j=96 compares to 3°S until 6°S,
    //                        from k=330 until k=355 compares to 30°W until 5°W)
        #pragma  omp parallel for collapse(3)
        for (int i = i_SCC_u; i < i_ECC_o - 1; i++) {
            for (int j = 93; j < 97; j++) {
                for (int k = 330; k < 355; k++) {
                    if (h.x[i][j][k] == 0.) {
                        {
                            u.x[i][j][k] = + .01 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_ECC_o - 1 - i_SCC_u);
                            v.x[i][j][k] = + .05 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_ECC_o - 1 - i_SCC_u);
//                            w.x[i][j][k] = + .05 * IC_water * static_cast<double>(i - i_SCC_u) / static_cast<double>(i_ECC_o - 1 - i_SCC_u);
                        }
                    }
                }
            }
        }

    }
/*
* 
*/
    void IC_DeepWater (Array &h, Array &u, Array &v, Array &w, Array &c) {
    // initial conditions for v- and w-velocity components in deep flows
    // Atlantic ocean
    // Thermohaline Conveyor Belt
    // downwards flow to deep flow                (from j=32 until j=36 compares to 58°N until 54°N,
    // Greenland until middel atlantic ridge       from k=317 until k=327 compares to 43°W until 33°W)
        k1 = 322;
        k3 = 334;
        k2 = (k3 - k1) / 2 + k1;
        kd = (k1 - k2) * (k1 - k2);

        #pragma  omp parallel for collapse(3) private(kn)
        for (int i = i_bottom; i < m.im-1; i++) {
            for (int j = 32; j < 37; j++) {
                for (int k = k1; k < k3; k++) {
                    if (h.x[i][j][k] == 0.) {
                        kn = k * k - 2 * k2 * k - k1 * k1 + 2 * k1 * k2;

                        c.x[i][j][k] = (ca_max - ca) * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom) * static_cast<double>(- kn) / static_cast<double>(kd) + ca;
                        u.x[i][j][k] = - .1 * IC_water * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom) * static_cast<double>(- kn) / static_cast<double>(kd);
                    }
                }
            }
        }


    // Atlantic ocean
    // Thermohaline Conveyor Belt
    // southerly directed current as deep flow (from j=32 until j=50 compares to 58°N until 40°N,
    // Greenland until middel atlantic ridge       from k=317 until k=327 compares to 43°W until 33°W)
        k1 = 322;
        k3 = 332;
        k2 = (k3 - k1) / 2 + k1;
        kd = (k1 - k2) * (k1 - k2);

        #pragma  omp parallel for collapse(3) private(kn)
        for (int i = i_bottom; i < i_deep + 1; i++) {
            for (int j = 32; j < 51; j++) {
                for (int k = k1; k <= k3; k++) {
                    if (h.x[i][j][k] == 0.) {
                        kn = k * k - 2 * k2 * k - k1 * k1 + 2 * k1 * k2;

                        c.x[i][j][k] = (ca_max - ca) * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom) * static_cast<double>(- kn) / static_cast<double>(kd) + ca;
                        v.x[i][j][k] = + .5 * IC_water * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom) * static_cast<double>(- kn) / static_cast<double>(kd);
                    }
                }
            }
        }


    // Atlantic ocean
    // Thermohaline Conveyor Belt
    // south of Newfoundland along the middle atlantic ridge
    // south-west directed flow as deep flow (from j=45 until j=60 compares to 30°N until 45°N,
    //                                                                   from k=306 until k=320 compares to 54°W until 40°W)
        j_beg = 45;
        j_end = 61;
        k_beg = 308;
        k_end = 323;
        j_run = 0;
        k_run = 0;
        j_step = 1;
        k_step = 11;

        while ((j_end - j_run) >= j_beg && (k_beg + k_run) <= k_end) {
            for (int j = std::min(j_end - j_run, m.jm - 1); j > (j_end - j_step - j_run); j--) {
                k1 = k_beg + k_run;
                k3 = k_beg + k_step + k_run;
                k2 = (k3 - k1) / 2 + k1;
                kd = (k1 - k2) * (k1 - k2);

                for (int k = (k_beg + k_run); k < (k_beg + k_step + k_run) && k < m.km && k >= 0; k++) {
                    for (int i = i_bottom; i < i_deep + 1; i++) {
                        if (h.x[i][j][k] == 0.) {
                            kn = k * k - 2 * k2 * k - k1 * k1 + 2 * k1 * k2;

                            c.x[i][j][k] = (ca_max - ca) * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom) * static_cast<double>(- kn) / static_cast<double>(kd) + ca;
                            v.x[i][j][k] = + .5 * IC_water * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom) * static_cast<double>(- kn) / static_cast<double>(kd);
                            w.x[i][j][k] = - .5 * IC_water * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom) * static_cast<double>(- kn) / static_cast<double>(kd);
                        }
                    }
                }
            k_run++;
            }
        j_run++;
        }


    // Atlantic
    // parallel to the north atlantic ridge
    // south of Newfoundland Guayana-fault (from j=60 until j=88 compares to 2°N until 30°N,
    //                                                                 from k=307 until k=316 compares to 53°W until 44°W)
        k1 = 309;
        k3 = 319;
        k2 = (k3 - k1) / 2 + k1;
        kd = (k1 - k2) * (k1 - k2);

        #pragma  omp parallel for collapse(3) private(kn)
        for (int i = i_bottom; i < i_deep + 1; i++) {
            for (int j = 60; j < 89; j++) {
                for (int k = k1; k <= k3; k++) {
                    if (h.x[i][j][k] == 0.) {
                        kn = k * k - 2 * k2 * k - k1 * k1 + 2 * k1 * k2;

                        c.x[i][j][k] = (ca_max - ca) * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom) * static_cast<double>(- kn) / static_cast<double>(kd) + ca;
                        v.x[i][j][k] = + IC_water * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom) * static_cast<double>(- kn) / static_cast<double>(kd);
                    }
                }
            }
        }


    // Atlantic
    // Southamerica     north coasts
    // below the Guayana-current as deep flow (from j=80 until j=100 compares to 10°N until 10°S,
    //                                                                     from k=307 until k=327 compares to 53°W until 33°W)
        j_beg = 80;
        j_end = 100;
        k_beg = 310;
        k_end = 331;
        j_run = 0;
        k_run = 0;
        j_step = 1;
        k_step = 13;

        while ((j_beg + j_run) <= j_end && (k_beg + k_run) <= k_end) {
            for (int j = j_beg + j_run; j < (j_beg + j_step + j_run) && j < m.jm; j++) {
                k1 = k_beg + k_run;
                k3 = k_beg + k_step + k_run;
                k2 = (k3 - k1) / 2 + k1;
                kd = (k1 - k2) * (k1 - k2);

                for (int k = k_beg + k_run; k < (k_beg + k_step + k_run) && k < m.km && k >= 0; k++) {
                    for (int i = i_bottom; i < i_deep + 1; i++) {
                        if (h.x[i][j][k] == 0.) {
                            kn = k * k - 2 * k2 * k - k1 * k1 + 2 * k1 * k2;

                            c.x[i][j][k] = (ca_max - ca) * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom) * static_cast<double>(- kn) / static_cast<double>(kd) + ca;
                            v.x[i][j][k] = + .5 * IC_water * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom) * static_cast<double>(- kn) / static_cast<double>(kd);
                            w.x[i][j][k] = + .5 * IC_water * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom) * static_cast<double>(- kn) / static_cast<double>(kd);
                        }
                    }
                }
            k_run++;
            }
        j_run++;
        }


    // Atlantic
    // Southamerica     east coasts
    // Brasil-current bending east (from j=107 until j=111 compares to 30°S until 35°S,
    //                                               from k=273 until k=298 compares to 53°W until 25°W)
        #pragma  omp parallel for collapse(3)
        for (int i = 0; i < i_beg; i++) {
            for (int j = 107; j < 112; j++) {
                for (int k = 273; k < 299; k++) {
                    if (h.x[i][j][k] == 0.) {
                        v.x[i][j][k] = IC_water * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom);
                        w.x[i][j][k] = IC_water * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom);
                    }
                }
            }
        }


    // Atlantic
    // Southamerica     east coasts
    // Parallel to the Brasil-current as deep flow (from j=99 until j=149 compares to 9°S until 59°S,
    //                                                                      from k=334 until k=339 compares to 26°W until 21°W)
        k1 = 330;
        k3 = 340;
        k2 = (k3 - k1) / 2 + k1;
        kd = (k1 - k2) * (k1 - k2);

        #pragma  omp parallel for collapse(3) private(kn)
        for (int i = i_bottom; i < i_deep + 1; i++) {
            for (int j = 99; j < 153; j++) {
                for (int k = k1; k <= k3; k++) {
                    if (h.x[i][j][k] == 0.) {
                        kn = k * k - 2 * k2 * k - k1 * k1 + 2 * k1 * k2;

                        c.x[i][j][k] = (ca_max - ca) * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom) * static_cast<double>(- kn) / static_cast<double>(kd) + ca;
                        v.x[i][j][k] = + IC_water * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom) * static_cast<double>(- kn) / static_cast<double>(kd);
                    }
                }
            }
        }


    //      %%%%%%%%%%%%%%%%%%%%%%%%%%%%     Pacific ocean     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    // Pacific ocean
    // Peru-current (Humboldt-current) coming from circumpolar current (from j=120 until j=148 compares to 30°S until 58°S,
    //                                                                                                              from k=277 until k=290 compares to 70°W until 83°W)
        k1 = 277;
        k3 = 291;
        k2 = (k3 - k1) / 2 + k1;
        kd = (k1 - k2) * (k1 - k2);

        #pragma  omp parallel for collapse(3) private(kn)
        for (int i = i_bottom; i < i_deep + 1; i++) {
            for (int j = 118; j < 153; j++) {
                for (int k = k1; k <= k3; k++) {
                    if (h.x[i][j][k] == 0.) {
                        kn = k * k - 2 * k2 * k - k1 * k1 + 2 * k1 * k2;

                        c.x[i][j][k] = (ca_max - ca) * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom) * static_cast<double>(- kn) / static_cast<double>(kd) + ca;
                        v.x[i][j][k] = - IC_water * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom) * static_cast<double>(- kn) / static_cast<double>(kd);
                    }
                }
            }
        }


    // Pacific ocean
    // deep flow to west out of Peru-current (Humboldt-current) (from j=120 until j=125 compares to 30°S until 35°S,
    // parallel to the equator                                                             from k=190 until k=285 compares to 75°W until 170°W)
        j1 = 118;
        j3 = 126;
        j2 = (j3 - j1) / 2 + j1;
        jd = (j1 - j2) * (j1 - j2);

        #pragma  omp parallel for collapse(3) private(jn)
        for (int i = i_bottom; i < i_deep + 1; i++) {
            for (int j = j1; j <= j3; j++) {
                for (int k = 193; k < 286; k++) {
                    if (h.x[i][j][k] == 0.) {
                        jn = j * j - 2 * j2 * j - j1 * j1 + 2 * j1 * j2;

                        c.x[i][j][k] = (ca_max - ca) * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom) * static_cast<double>(- jn) / static_cast<double>(jd) + ca;
                        w.x[i][j][k] = - IC_water * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom) * static_cast<double>(- jn) / static_cast<double>(jd);
                    }
                }
            }
        }


    // Pacific ocean
    // Newseeland in the east coming from circumpolar current (from j=92 until j=152 compares to 2°S until 62°S,
    // in direction to Japan                                                             from k=188 until k=197 compares to 172°W until 163°W)
        k1 = 188;
        k3 = 198;
        k2 = (k3 - k1) / 2 + k1;
        kd = (k1 - k2) * (k1 - k2);

        #pragma  omp parallel for collapse(3) private(kn)
        for (int i = i_bottom; i < i_deep + 1; i++) {
            for (int j = 92; j < 153; j++) {
                for (int k = k1; k <= k3; k++) {
                    if (h.x[i][j][k] == 0.) {
                        kn = k * k - 2 * k2 * k - k1 * k1 + 2 * k1 * k2;

                        c.x[i][j][k] = (ca_max - ca) * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom) * static_cast<double>(- kn) / static_cast<double>(kd) + ca;
                        v.x[i][j][k] = - IC_water * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom) * static_cast<double>(- kn) / static_cast<double>(kd);
                    }
                }
            }
        }


    // Pacific ocean
    // equator until Japan as deep flow (from j=68 until j=93 compares to 3°S until 22°N,
    //                                                        from k=162 until k=193 compares to 167°W until 162°O)

        j_beg = 68;
        j_end = 94;
        k_beg = 162;
        k_end = 194;
        j_run = 0;
        k_run = 0;
        j_step = 1;
        k_step = 10;

        while ((j_beg + j_run) <= j_end && (k_beg + k_run) <= k_end) {
            for (int j = j_beg + j_run; j < (j_beg + j_step + j_run) && j < m.jm; j++) {
                k1 = k_beg + k_run;
                k3 = k_beg + k_step + k_run;
                k2 = (k3 - k1) / 2 + k1;
                kd = (k1 - k2) * (k1 - k2);

                for (int k = k_beg + k_run; k < (k_beg + k_step + k_run) && k < m.km && k >= 0; k++) {
                    for (int i = i_bottom; i < i_deep + 1; i++) {
                        if (h.x[i][j][k] == 0.) {
                            kn = k * k - 2 * k2 * k - k1 * k1 + 2 * k1 * k2;

                            c.x[i][j][k] = (ca_max - ca) * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom) * static_cast<double>(- kn) / static_cast<double>(kd) + ca;
                            v.x[i][j][k] = - .5 * IC_water * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom) * static_cast<double>(- kn) / static_cast<double>(kd);
                            w.x[i][j][k] = - .5 * IC_water * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom) * static_cast<double>(- kn) / static_cast<double>(kd);
                        }
                    }
                }
            k_run++;
            }
        j_run++;
        }


    // Pacific ocean
    // east of Japan
    // upward flow east of Japan (from j=54 until j=67 compares to 23°N until 36°N
    //                                              from k=162 until k=171 compares to 162°O until 171°O)
        k1 = 162;
        k3 = 172;
        k2 = (k3 - k1) / 2 + k1;
        kd = (k1 - k2) * (k1 - k2);

        #pragma  omp parallel for collapse(3) private(kn)
        for (int i = i_bottom; i < i_deep + 1; i++) {
            for (int j = 54; j < 68; j++) {
                for (int k = k1; k <= k3; k++) {
                    if (h.x[i][j][k] == 0.) {
                        kn = k * k - 2 * k2 * k - k1 * k1 + 2 * k1 * k2;

                        c.x[i][j][k] = (ca_max - ca) * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom) * static_cast<double>(- kn) / static_cast<double>(kd) + ca;
                        v.x[i][j][k] = - .5 * IC_water * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom) * static_cast<double>(- kn) / static_cast<double>(kd);
                    }
                }
            }
        }


    // Pacific ocean
    // east of Japan
    // upward flow east of Japan (from j=50 until j=58 compares to 32°N until 40°N
    //                                              from k=162 until k=210 compares to 162°O until 150°W)
        j1 = 54;
        j3 = 64;
        j2 = (j3 - j1) / 2 + j1;
        jd = (j1 - j2) * (j1 - j2);

        #pragma  omp parallel for collapse(3) private(jn)
        for (int i = i_deep-1; i < m.im - 1; i++) {
            for (int j = j1; j <= j3; j++) {
                for (int k = 171; k < 210; k++) {
                    if (h.x[i][j][k] == 0.) {
                        jn = j * j - 2 * j2 * j - j1 * j1 + 2 * j1 * j2;

                        c.x[i][j][k] = (ca_max - ca) * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom) * static_cast<double>(- jn) / static_cast<double>(jd) + ca;
                        u.x[i][j][k] = + .1 * IC_water * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom) * static_cast<double>(- kn) / static_cast<double>(kd);
                        w.x[i][j][k] = + .5 *  IC_water * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom) * static_cast<double>(- jn) / static_cast<double>(jd);
                    }
                }
            }
        }


    //  &&&&&&&&&&&&&&&&&&&&&&&&&&&     Indic ocean     &&&&&&&&&&&&&&&&&&&&&&&&&&

    // Indic ocean
    // south east of Southafrica        east coasts
    // south east of Southafrica coming from  the circumpolar current (from j=83 until j=103 compares to 3°S until 26°S,
    // prolongation of the circumpolar current                                         from k=55 until k=58 compares to 62°O until 65°O)

        #pragma  omp parallel for collapse(3)
        for (int i = i_bottom; i < i_deep + 1; i++) {
            for (int j = 83; j < 104; j++) {
                for (int k = 55; k < 59; k++) {
                if (h.x[i][j][k] == 0.) {
                        c.x[i][j][k] = ca_max;
//                        u.x[i][j][k] = + IC_water * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom);
                        v.x[i][j][k] = - IC_water * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom);
                        w.x[i][j][k] = + IC_water * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom);
                    }
                }
            }
        }


    // Indic ocean
    // south east of Southafrica        east coasts
    // south east of Southafrica coming from  the circumpolar current (from j=93 until j=152 compares to 3°S until 62°S,
    //                                                                                                       from k=55 until k=60 compares to 55°O until 60°O)
        k1 = 55;
        k3 = 61;
        k2 = (k3 - k1) / 2 + k1;
        kd = (k1 - k2) * (k1 - k2);

        #pragma  omp parallel for collapse(3) private(kn)
        for (int i = i_bottom; i < i_deep + 1; i++) {
            for (int j = 93; j < 153; j++) {
                for (int k = k1; k <= k3; k++) {
                if (h.x[i][j][k] == 0.) {
                        kn = k * k - 2 * k2 * k - k1 * k1 + 2 * k1 * k2;

                        c.x[i][j][k] = (ca_max - ca) * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom) * static_cast<double>(- kn) / static_cast<double>(kd) + ca;
                        v.x[i][j][k] = - IC_water * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom) * static_cast<double>(- kn) / static_cast<double>(kd);
                    }
                }
            }
        }


    // Indic ocean
    // couth of India along the equator
    // deep flow south of India (from j=93 until j=97 compares to 3°S until 7°S,
    // upwards flow                    from k=64 until k=95 compares to 64°O until 95°O)
        j1 = 93;
        j3 = 99;
        j2 = (j3 - j1) / 2 + j1;
        jd = (j1 - j2) * (j1 - j2);

        #pragma  omp parallel for collapse(3) private(jn)
        for (int i = i_deep-1; i < m.im - 2; i++) {
            for (int j = j1; j <= j3; j++) {
                for (int k = 60; k < 96; k++) {
                    if (h.x[i][j][k] == 0.) {
                        jn = j * j - 2 * j2 * j - j1 * j1 + 2 * j1 * j2;

                        c.x[i][j][k] = (ca_max - ca) * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom) * static_cast<double>(- jn) / static_cast<double>(jd) + ca;
                        u.x[i][j][k] = + .1 * IC_water * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom) * static_cast<double>(- kn) / static_cast<double>(kd);
                        w.x[i][j][k] = + .5 *  IC_water * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom) * static_cast<double>(- jn) / static_cast<double>(jd);
                    }
                }
            }
        }


    // Indic ocean
    // west of Corea along the equator
    // upwards flow south of India (from j=84 until j=89 compares to 5°S until 10°S,
    // 90°W surface flow                  from k=78 until k=82 compares to 88°O until 92°O)
        #pragma  omp parallel for collapse(3)
        for (int i = 0; i < m.im; i++) {
            for (int j = 84; j < 90; j++) {
                for (int k = 78; k < 83; k++) {
                    { {
                            u.x[i][j][k] = + IC_water * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom);
                            v.x[i][j][k] = + IC_water * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom);
                            w.x[i][j][k] = + IC_water * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom);
                        }
                    }
                }
            }
        }


    // Indic ocean
    // south of Corea
    // 90°W surface flow south of India (from j=84 until j=96 compares to 5°S until 18°S,
    //                                                        from k=82 until k=85 compares to 92°O until 96°O)
        #pragma  omp parallel for collapse(3)
        for (int i = i_beg; i < m.im; i++) {
            for (int j = 84; j < 97; j++) {
                for (int k = 82; k < 86; k++) {
                    {
                        v.x[i][j][k] = + IC_water * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom);
                        w.x[i][j][k] = - IC_water * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom);
                    }
                }
            }
        }


    // Indic ocean
    // south-east of Southafric     east coasts
    // 90°W deep flow from circumpolar current (from j=89 until j=132 compares to 10°S until 59°S,
    //                                                                      from k=74 until k=78 compares to 83°O until 88°O)
        #pragma  omp parallel for collapse(3)
        for (int i = 0; i < i_beg; i++) {
            for (int j = 89; j < 133; j++) {
                for (int k = 74; k < 79; k++) {
                    { 
                        v.x[i][j][k] = - IC_water * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom);
                        w.x[i][j][k] = - IC_water * static_cast<double>(i - i_bottom) / static_cast<double>(i_deep - i_bottom);
                    }
                }
            }
        }
    }

private:
    cHydrosphereModel& m;
};
