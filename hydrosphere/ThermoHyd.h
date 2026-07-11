#pragma once

#include "cHydrosphereModel.h"
#include "Utils.h"

#include <vector>
#include <cmath>
#include <chrono>
#include <cstdio>
#include <iostream>
#include <iomanip>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace AtomUtils;

class ThermoHyd {
public:
    explicit ThermoHyd(cHydrosphereModel& model)
        : m(model)
    {}

    // ------------------------------------------------------------------
    void SaltWaterDens()
    {
        using namespace std;
        cout << endl << "      PresStat_SaltWaterDens" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        // hydrostatic pressure, equations of state for water and salt water density
        // as functions of salinity, temperature and hydrostatic pressure
        constexpr double E_water   = 2.15e9;                            // [N/m²]
        constexpr double beta_water = 8.8e-5;                           // [m³/(m³·°C)]

        const double step_m   = m.L_hyd / (double)(m.im-1);
        const double step_km  = step_m / 1000.0;

        // --- hydrostatic pressure, water and salt water density at the surface ---
        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                // Surface density initialised first; surface depth is 0,
                // so hydrostatic pressure equals atmospheric pressure only.
                m.r_water.x[m.im-1][j][k] = m.r_0_water;                // [kg/m³]

                m.p_hydro.x[m.im-1][j][k] = m.p_0 / 1000.0;             // [bar] atmospheric pressure
                if(is_land(m.h, m.im-1, j, k)) m.p_hydro.x[m.im-1][j][k] = 0.0;


                double t_Celsius   = m.t.x[m.im-1][j][k] * m.t_0 - m.t_0;
                double C_p         = 999.83;
                double beta_p      = 0.808;
                double alfa_t_p    = 0.0708 * (1.0 + 0.068 * t_Celsius);
                double gamma_t_p   = 0.003  * (1.0 - 0.012 * t_Celsius);

                m.r_salt_water.x[m.im-1][j][k] =                        // [kg/m³], Gill approximation
                    C_p + beta_p * m.c.x[m.im-1][j][k] * m.c_35
                    - alfa_t_p * t_Celsius
                    - gamma_t_p * (35.0 - m.c.x[m.im-1][j][k] * m.c_35) * t_Celsius;
            }
        }

        // --- hydrostatic pressure, water and salt water density in the flow field ---
        // outer (j,k) parallel; inner i-loop is sequential (downward dependency chain)
        #pragma omp parallel for collapse(2) schedule(dynamic, 4)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                for (int i = m.im-2; i >= 0; i--) {

                    double d_i         = (double)(m.im-1 - i);
                    double t_Celsius_1 = m.t.x[i][j][k]   * m.t_0 - m.t_0;
                    double t_Celsius_0 = m.t.x[i+1][j][k] * m.t_0 - m.t_0;

                    // Incremental hydrostatic step: p[i] = p[i+1] + rho * g * dz
                    // (previous formula used rho[i+1]*g*d_i*step_m, which sets the
                    // entire column density to a single layer value and — combined
                    // with the wrong surface pressure — yielded negative ΔP at
                    // i=im-2, inverting the compressibility term and producing
                    // r_water values ~1 kg/m³ too low throughout the column.)
                    m.p_hydro.x[i][j][k] =                              // hydrostatic pressure in [bar]
                        m.p_hydro.x[i+1][j][k]
                        + m.r_water.x[i+1][j][k] * m.g * step_m / 100000.0;
                    if(is_land(m.h, i, j, k)) m.p_hydro.x[i][j][k] = 0.0;

                    m.r_water.x[i][j][k] =                              // [kg/m³]
                        m.r_water.x[i+1][j][k]
                        / (1.0 + beta_water * (t_Celsius_1 - t_Celsius_0))
                        / (1.0 - (m.p_hydro.x[i][j][k]
                                - m.p_hydro.x[i+1][j][k]) / E_water * 1e5);

                    double p_km      = d_i * step_km;                   // depth in [km]
                    double C_p       = 999.83 + 5.053 * p_km - 0.048 * p_km * p_km;
                    double beta_p    = 0.808 - 0.0085 * p_km;
                    double alfa_t_p  = 0.0708 * (1.0 + 0.351 * p_km
                                     + 0.068 * (1.0 - 0.0683 * p_km) * t_Celsius_1);
                    double gamma_t_p = 0.003 * (1.0 - 0.059  * p_km
                                     - 0.012 * (1.0 - 0.064  * p_km) * t_Celsius_1);

                    m.r_salt_water.x[i][j][k] =                         // [kg/m³], Gill approximation
                        C_p + beta_p * m.c.x[i][j][k] * m.c_35
                        - alfa_t_p * t_Celsius_1
                        - gamma_t_p * (35.0 - m.c.x[i][j][k] * m.c_35) * t_Celsius_1;

                    if (is_land(m.h, i, j, k))
                        m.r_salt_water.x[i][j][k] = m.r_0_saltwater;
                }  // i
            }  // k
        }  // j

        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for PresStat_SaltWaterDens\n", elapsed.count() * 1e-9);
        cout << "      PresStat_SaltWaterDens ended" << endl;
    }
/*
* 
*/
    // ------------------------------------------------------------------
    void SalinityEvaporation()
    {
        using namespace std;
        cout << endl << "      SalinityEvaporation" << endl;

        // Virtual salt flux BC (Huang, Ocean Circulation, p. 165 — 1-D model in r-direction):
        //   Dc = (rho_sw/rho_w) * c_surf * (E-P)/(E-P)_ref * (1 - 2*c_s)
        // where c = S/c_35 [dimensionless], c_s = (c_35/1000)*c [salt mass fraction ~0.035].
        // Units: [kg/m3]/[kg/m3] * [1] * [mm/d]/[mm/d] * [1] = [1]  (dimensionless)
        auto begin = std::chrono::high_resolution_clock::now();

        // Reference E-P MAGNITUDE for the non-dimensional flux. The old code used
        // the SIGNED global mean E-P, which is ~0 in a balanced climate and
        // EXACTLY 0 with no atmosphere coupling (Evaporation/Precipitation_2D are
        // reset to 0 when the transfer is skipped). That made 1/(rho*E_P_ref) blow
        // up to Inf and drove salinity_evaporation to NaN, and a negative mean
        // flipped the flux sign. Use the mean of |E-P| over ocean as a robust
        // positive scale, floored so the division is always finite; with no E-P
        // data every local (E-P) is 0 so the flux is 0.
        double ep_abs_sum = 0.0;
        long   ep_n       = 0;
        for (int j = 0; j < m.jm; j++)
            for (int k = 0; k < m.km; k++)
                if (is_water(m.h, m.im-1, j, k)) {
                    ep_abs_sum += fabs(m.Evaporation.y[j][k] - m.Precipitation_2D.y[j][k]);
                    ++ep_n;
                }
        const double E_P_ref   = std::max(ep_n ? ep_abs_sum / (double)ep_n : 0.0,
                                          1.0e-3);                      // [mm/d], floored positive
        const double coeff_c   = m.c_35 / 1000.0;                       // salt mass fraction at c=1 [1]
        const double inv_r0_E0 = 1.0 / (m.r_0_water * E_P_ref);        // [m3/(kg*mm/d)]

        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {

                if (is_land(m.h, m.im-1, j, k)) {
                    m.c.x[m.im-1][j][k] = 0.0;
                    continue;
                }

                // NOTE: c_fix is the FIXED IC surface-salinity baseline captured
                // once at init (cHydrosphereModel), NOT the running value. Do not
                // overwrite it here — re-reading the current surface salinity each
                // call turns the flux into an unbounded accumulation that pins the
                // surface to the 45-psu clamp (E>P) or 0 psu (P>E).

                double evap_precip = m.Evaporation.y[j][k]
                                   - m.Precipitation_2D.y[j][k];        // [mm/d]
                double c_surf      = m.c.x[m.im-1][j][k];               // [1]

                // Dc = (rho_sw/rho_w) * c_surf * (E-P)/(E-P)_ref * (1 - 2*c_s)
                m.salinity_evaporation.y[j][k] =
                    m.r_salt_water.x[m.im-1][j][k] * inv_r0_E0
                  * coeff_c * c_surf * (1.0 - 2.0 * coeff_c * c_surf)
                  * evap_precip;

                // clamp to ±1.0 psu
                if (m.salinity_evaporation.y[j][k] * m.c_35 >=  1.0)
                    m.salinity_evaporation.y[j][k] =  1.0 / m.c_35;
                if (m.salinity_evaporation.y[j][k] * m.c_35 <= -1.0)
                    m.salinity_evaporation.y[j][k] = -1.0 / m.c_35;

                // Apply surface flux BC: only the surface cell is updated.
                // Dynamics (advection/diffusion) propagate the signal downward.
                double c_eq = m.c_fix.y[j][k] + m.salinity_evaporation.y[j][k];
                m.c.x[m.im-1][j][k] = std::min(45.0 / m.c_35,
                    std::max(0.0, c_eq));
            }  // k
        }  // j



        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for SalinityEvaporation\n", elapsed.count() * 1e-9);
        cout << "      SalinityEvaporation ended" << endl;
    }
/*
* 
*/
    // ------------------------------------------------------------------
    void forces()
    {
        using namespace std;
        cout << endl << endl << endl << "      ATOM: Forces" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        const int    i_max   = m.im - 1;
        const double dr_phys = m.L_hyd / i_max;                         // [m] uniform depth step
        const double two_omega = 2.0 * m.omega;                         // [rad/s]
        const double omega2    = m.omega * m.omega;
        const double pres_coeff = 1e2 * m.p_0;  // [Pa] hPa→Pa; p_dyn is dimensionless, dpdr is [1/m]
        const double r_Earth_m = m.r_Earth * 1e3;      // [m]

        // precompute trig tables to avoid redundant sin/cos in the inner loop
        std::vector<double> sinthe_table(m.jm);
        std::vector<double> costhe_table(m.jm);
        for (int j = 0; j < m.jm; j++) {
            sinthe_table[j] = sin(m.the.z[j]);
            costhe_table[j] = cos(m.the.z[j]);
        }

        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 1; i < m.im-1; i++) {
            for (int j = 1; j < m.jm-1; j++) {

                double rm           = m.rad.z[i] * m.L_hyd / (double)(m.im - 1);// [m]
                double inv_rm       = 1.0 / rm;
                double sinthe       = sinthe_table[j];
                double costhe       = costhe_table[j];
                double inv_rmsinthe = 1.0 / (rm * sinthe);              // [m]
                double abs_sinthe   = fabs(sinthe);
                double depth        = i * dr_phys;                      // [m] depth at level i
                double r_surface_eq = m.r_Earth * 1.0e3 + depth;// [m] radius at equator and at any latitude, level i
                double U_Earth      = r_surface_eq * costhe * 2.0 * M_PI;// [m]
                double inv_2dr      = 0.5 * (double)(m.im - 1) / m.L_hyd;// [1/m]
                double inv_2dthe    = 0.5 * (m.jm-1) / U_Earth;         // [1/m]
                double inv_2dphi    = 0.5 * (m.km-1) / U_Earth;         // [1/m]

                double rad_dist  = (double)i * m.L_hyd / (double)(m.im - 1); // [m]
                double rad_Earth = rad_dist + r_Earth_m;                // [m]

                for (int k = 1; k < m.km-1; k++) {

                    // Coriolis force components
                    double w_ijk   = m.w.x[i][j][k] * m.u_0;            // [m/s]
                    double Cor_r   = - two_omega * sinthe * w_ijk;      // [rad/s·m/s]
                    double Cor_the = + two_omega * costhe * w_ijk;
                    double Cor_phi = + two_omega * (-costhe * m.v.x[i][j][k]
                                                   + sinthe * m.u.x[i][j][k]) * m.u_0;

                    m.CoriolisForce.x[i][j][k] = m.Coriolis * m.r_0_saltwater
                        * sqrt(Cor_r * Cor_r + Cor_the * Cor_the + Cor_phi * Cor_phi); // [N/m³]

                    m.CentrifugalForce.x[i][j][k] = m.centrifugal
                        * m.r_water.x[i][j][k] * omega2 * rad_Earth
                        * (1.0 + abs_sinthe);                           // [N/m³]

                    m.BuoyancyForce.x[i][j][k] =
                        - 1.0e-3 * m.buoyancy * m.r_salt_water.x[i][j][k] * m.g;  // [kN/m³]

                    // pressure gradient (central differences)
                    double dpdr   = (m.p_dyn.x[i+1][j][k] - m.p_dyn.x[i-1][j][k]) * inv_2dr; // [1/m]
                    double dpdthe = (m.p_dyn.x[i][j+1][k] - m.p_dyn.x[i][j-1][k]) * inv_2dthe * inv_rm;
                    double dpdphi = (m.p_dyn.x[i][j][k+1] - m.p_dyn.x[i][j][k-1]) * inv_2dphi * inv_rmsinthe;

                    m.PresGradForce.x[i][j][k] =
                        sqrt((dpdr * dpdr + dpdthe * dpdthe + dpdphi * dpdphi) / 3.0)
                        * pres_coeff * 1.0e-3;                          // [kN/m³]
                    if(j == 90)  m.PresGradForce.x[i][j][k] = 
                          0.5 * (m.PresGradForce.x[i][j+1][k] 
                               + m.PresGradForce.x[i][j-1][k]);

                    if (is_land(m.h, i, j, k)) {
                        m.BuoyancyForce.x[i][j][k]         = 0.0;
                        m.PresGradForce.x[i][j][k] = 0.0;
                        m.CoriolisForce.x[i][j][k]         = 0.0;
                    }
                }  // k
            }  // j
        }  // i


        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for Forces\n", elapsed.count() * 1e-9);
        cout << "      ATOM: Forces ended" << endl;
    }
/*
* 
*/
    // ------------------------------------------------------------------
    // Ekman pumping, upwelling/downwelling, salt column diagnostics,
    // and force boundary conditions at the ocean surface (i = im-1).
    //
    // Notes:
    //   • The shared Array_1D scratch buffers (aux_grad_v/w) are replaced
    //     with thread-local std::vector<double> so the outer (j,k) loop
    //     can be safely parallelised with OpenMP.
    //   • Local c43/c13 are double to avoid the truncation bug in
    //     cHydrosphereModel's const int members.
    // ------------------------------------------------------------------
    void runDataHyd()
    {
        using namespace std;
        cout << endl << "      RunDataHyd" << endl;
        auto begin = chrono::high_resolution_clock::now();

        const int    i_max  = m.im-1;
        const double c43    = 4.0 / 3.0;
        const double c13    = 1.0 / 3.0;

        // ----------------------------------------------------------------
        // 1. Initialise 2D surface diagnostic arrays
        // ----------------------------------------------------------------
        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                m.Upwelling.y[j][k]       = 0.0;
                m.Downwelling.y[j][k]     = 0.0;
                m.EkmanPumping.y[j][k]    = 0.0;
            }
        }

        // ----------------------------------------------------------------
        // 2. Depth-integrated Ekman transport (trapezoidal rule)
        //    Thread-local buffers replace the shared Array_1D members.
        // ----------------------------------------------------------------
        const double dr_rm = m.dr * m.L_hyd / (double)i_max;

        #pragma omp parallel
        {
            vector<double> grad_v(m.im), grad_w(m.im);

            #pragma omp for collapse(2) schedule(dynamic, 4)
            for (int j = 0; j < m.jm; j++) {
                for (int k = 0; k < m.km; k++) {
                    for (int i = 0; i <= i_max; i++) {
                        if (is_land(m.h, i, j, k)) {
                            grad_v[i] = 0.;
                            grad_w[i] = 0.;
                        } else {
                            grad_v[i] = m.r_salt_water.x[i][j][k] * m.v.x[i][j][k] * m.u_0;
                            grad_w[i] = m.r_salt_water.x[i][j][k] * m.w.x[i][j][k] * m.u_0;
                        }
                    }
                    double sv = 0., sw = 0.;
                    for (int i = 0; i < i_max; i++) {
                        sv += grad_v[i] + grad_v[i + 1];
                        sw += grad_w[i] + grad_w[i + 1];
                    }
                    m.aux_v.x[i_max][j][k] = 0.5 * sv * dr_rm;
                    m.aux_w.x[i_max][j][k] = 0.5 * sw * dr_rm;

                    if (is_land(m.h, i_max, j, k)) {
                        m.aux_v.x[i_max][j][k] = 0.;
                        m.aux_w.x[i_max][j][k] = 0.;
                    }
                }
            }
        } // end parallel

        // ----------------------------------------------------------------
        // 3. Boundary conditions for aux_v / aux_w at the ocean surface
        // ----------------------------------------------------------------
        for (int k = 0; k < m.km; k++) {
            m.aux_v.x[i_max][0][k]        = c43 * m.aux_v.x[i_max][1][k]        - c13 * m.aux_v.x[i_max][2][k];
            m.aux_v.x[i_max][m.jm-1][k]   = c43 * m.aux_v.x[i_max][m.jm-2][k]   - c13 * m.aux_v.x[i_max][m.jm-3][k];
            m.aux_w.x[i_max][0][k]        = c43 * m.aux_w.x[i_max][1][k]        - c13 * m.aux_w.x[i_max][2][k];
            m.aux_w.x[i_max][m.jm-1][k]   = c43 * m.aux_w.x[i_max][m.jm-2][k]   - c13 * m.aux_w.x[i_max][m.jm-3][k];
        }
        for (int j = 0; j < m.jm; j++) {
            double av0 = c43 * m.aux_v.x[i_max][j][1]      - c13 * m.aux_v.x[i_max][j][2];
            double av1 = c43 * m.aux_v.x[i_max][j][m.km-2] - c13 * m.aux_v.x[i_max][j][m.km-3];
            m.aux_v.x[i_max][j][0] = m.aux_v.x[i_max][j][m.km-1] = (av0 + av1) * 0.5;

            double aw0 = c43 * m.aux_w.x[i_max][j][1]      - c13 * m.aux_w.x[i_max][j][2];
            double aw1 = c43 * m.aux_w.x[i_max][j][m.km-2] - c13 * m.aux_w.x[i_max][j][m.km-3];
            m.aux_w.x[i_max][j][0] = m.aux_w.x[i_max][j][m.km-1] = (aw0 + aw1) * 0.5;
        }

        // ----------------------------------------------------------------
        // 4. Ekman pumping = divergence of depth-integrated transport
        //    Wall-adjacent cells use one-sided 3-point stencils;
        //    all other cells use central differences.
        // ----------------------------------------------------------------
        const double coeff_pumping = 864.0;  // converts m/s to cm/d

        #pragma omp parallel for collapse(2) schedule(dynamic, 4)
        for (int j = 1; j < m.jm-1; j++) {
            for (int k = 1; k < m.km-1; k++) {
                const double sinthe   = sin(m.the.z[j]);
                const double rm       = m.L_hyd;
                const double rmsinthe = rm * sinthe;
                const double h_d      = is_water(m.h, i_max, j, k) ? 1.0 : 0.0;

                // Helper: surface-layer transport shorthands
                auto Vj = [&](int dj) -> double { return m.aux_v.x[i_max][j + dj][k]; };
                auto Wk = [&](int dk) -> double { return m.aux_w.x[i_max][j][k + dk]; };

                // d(aux_v)/dthe — prefer backward stencil when available at walls
                double dvdthe;
                bool the_flag = false;
                if ((j >= 2) && (j < m.jm-3)) {
                    if (is_land(m.h, 0, j, k) && is_air(m.h, 0, j+1, k) && is_air(m.h, 0, j+2, k)) {
                        dvdthe = h_d * (-3.0 * Vj(0) + 4.0 * Vj(+1) - Vj(+2)) / (2.0 * m.dthe);
                        the_flag = true;
                    }
                    if (is_land(m.h, 0, j, k) && is_air(m.h, 0, j-1, k) && is_air(m.h, 0, j-2, k)) {
                        dvdthe = -h_d * (-3.0 * Vj(0) + 4.0 * Vj(-1) - Vj(-2)) / (2.0 * m.dthe);
                        the_flag = true;
                    }
                }
                if (!the_flag)
                    dvdthe = h_d * (Vj(+1) - Vj(-1)) / (2.0 * m.dthe);

                // d(aux_w)/dphi — same pattern
                double dwdphi;
                bool phi_flag = false;
                if ((k >= 2) && (k < m.km-3)) {
                    if (is_land(m.h, 0, j, k) && is_air(m.h, 0, j, k+1) && is_air(m.h, 0, j, k+2)) {
                        dwdphi = h_d * (-3.0 * Wk(0) + 4.0 * Wk(+1) - Wk(+2)) / (2.0 * m.dphi);
                        phi_flag = true;
                    }
                    if (is_land(m.h, 0, j, k) && is_air(m.h, 0, j, k-1) && is_air(m.h, 0, j, k-2)) {
                        dwdphi = -h_d * (-3.0 * Wk(0) + 4.0 * Wk(-1) - Wk(-2)) / (2.0 * m.dphi);
                        phi_flag = true;
                    }
                }
                if (!phi_flag)
                    dwdphi = h_d * (Wk(+1) - Wk(-1)) / (2.0 * m.dphi);

                m.EkmanPumping.y[j][k] = -coeff_pumping
                    / m.r_salt_water.x[i_max][j][k]
                    * (dvdthe / rm + dwdphi / rmsinthe);

                if (is_land(m.h, i_max, j, k))
                    m.EkmanPumping.y[j][k] = 0.0;

                m.Upwelling.y[j][k]   = (m.EkmanPumping.y[j][k] > 0.0) ? m.EkmanPumping.y[j][k] : 0.0;
                m.Downwelling.y[j][k] = (m.EkmanPumping.y[j][k] < 0.0) ? m.EkmanPumping.y[j][k] : 0.0;
            }
        }

        // ----------------------------------------------------------------
        // 5. Boundary conditions for EkmanPumping / Upwelling / Downwelling
        // ----------------------------------------------------------------
        for (int k = 0; k < m.km; k++) {
            m.EkmanPumping.y[0][k]       = c43 * m.EkmanPumping.y[1][k]       - c13 * m.EkmanPumping.y[2][k];
            m.EkmanPumping.y[m.jm-1][k]  = c43 * m.EkmanPumping.y[m.jm-2][k]  - c13 * m.EkmanPumping.y[m.jm-3][k];
            m.Upwelling.y[0][k]          = c43 * m.Upwelling.y[1][k]          - c13 * m.Upwelling.y[2][k];
            m.Upwelling.y[m.jm-1][k]     = c43 * m.Upwelling.y[m.jm-2][k]     - c13 * m.Upwelling.y[m.jm-3][k];
            m.Downwelling.y[0][k]        = c43 * m.Downwelling.y[1][k]        - c13 * m.Downwelling.y[2][k];
            m.Downwelling.y[m.jm-1][k]   = c43 * m.Downwelling.y[m.jm-2][k]   - c13 * m.Downwelling.y[m.jm-3][k];
        }
        for (int j = 0; j < m.jm; j++) {
            auto periodic_bc_2D = [&](Array_2D& a) {
                double v0 = c43 * a.y[j][1]      - c13 * a.y[j][2];
                double v1 = c43 * a.y[j][m.km-2] - c13 * a.y[j][m.km-3];
                a.y[j][0] = a.y[j][m.km-1] = (v0 + v1) * 0.5;
            };
            periodic_bc_2D(m.EkmanPumping);
            periodic_bc_2D(m.Upwelling);
            periodic_bc_2D(m.Downwelling);
        }

        // ----------------------------------------------------------------
        // 7a. Salt / force BCs at i = 0 (linear) and i = im-1 (cubic),
        //     plus surface forces at i_max — all (j,k)-independent.
        // ----------------------------------------------------------------
        const double coeff_u_p = m.r_0_water * m.u_0 * m.u_0 / m.L_hyd;
        const double coeff_Cor = m.r_0_water * m.u_0 * m.omega;

        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {

                // i = 0: cubic extrapolation of the surface forces
                m.BuoyancyForce.x[0][j][k] =
                      m.BuoyancyForce.x[3][j][k]
                    - 3.0 * m.BuoyancyForce.x[2][j][k]
                    + 3.0 * m.BuoyancyForce.x[1][j][k];
                m.CoriolisForce.x[0][j][k] =
                      m.CoriolisForce.x[3][j][k]
                    - 3.0 * m.CoriolisForce.x[2][j][k]
                    + 3.0 * m.CoriolisForce.x[1][j][k];
                m.PresGradForce.x[0][j][k] =
                      m.PresGradForce.x[3][j][k]
                    - 3.0 * m.PresGradForce.x[2][j][k]
                    + 3.0 * m.PresGradForce.x[1][j][k];

                // Surface Coriolis force
                const double sinthe = sin(m.the.z[j]);
                double costhe = cos(m.the.z[j]);
                if (j > 90) costhe = -costhe;

                const double cor_r   =  2.0 *  costhe * m.w.x[i_max][j][k];
                const double cor_the = -2.0 *  sinthe * m.w.x[i_max][j][k];
                const double cor_phi =  2.0 * (sinthe * m.v.x[i_max][j][k]
                                             - costhe * m.u.x[i_max][j][k]);

                double dpdr = m.p_hydro.x[i_max][j][k] = m.p_hydro.x[m.im-4][j][k] 
                    - 3.0 * m.p_hydro.x[m.im-3][j][k] + 3.0 * m.p_hydro.x[m.im-2][j][k];  // extrapolation

                m.CoriolisForce.x[i_max][j][k] = coeff_Cor
                    * sqrt((cor_r*cor_r + cor_the*cor_the + cor_phi*cor_phi) / 3.0);

                m.BuoyancyForce.x[i_max][j][k] =
                    m.r_0_water * (m.t.x[i_max][j][k] - 1.0) * m.g;

                m.PresGradForce.x[i_max][j][k] = -coeff_u_p * dpdr;

                if (is_land(m.h, i_max, j, k)) {
                    m.BuoyancyForce.x[i_max][j][k]         = 0.0;
                    m.PresGradForce.x[i_max][j][k]  = 0.0;
                    m.CoriolisForce.x[i_max][j][k]          = 0.0;
                }
            }
        }

        // ----------------------------------------------------------------
        // 7b. j-direction (pole) BCs for all 3D force arrays
        // ----------------------------------------------------------------
        #pragma omp parallel for schedule(static)
        for (int k = 0; k < m.km; k++) {
            for (int i = 0; i < m.im; i++) {
                m.BuoyancyForce.x[i][0][k]             = c43 * m.BuoyancyForce.x[i][1][k]       - c13 * m.BuoyancyForce.x[i][2][k];
                m.BuoyancyForce.x[i][m.jm-1][k]        = c43 * m.BuoyancyForce.x[i][m.jm-2][k]  - c13 * m.BuoyancyForce.x[i][m.jm-3][k];
                m.CoriolisForce.x[i][0][k]             = c43 * m.CoriolisForce.x[i][1][k]       - c13 * m.CoriolisForce.x[i][2][k];
                m.CoriolisForce.x[i][m.jm-1][k]        = c43 * m.CoriolisForce.x[i][m.jm-2][k]  - c13 * m.CoriolisForce.x[i][m.jm-3][k];
                m.PresGradForce.x[i][0][k]     = c43 * m.PresGradForce.x[i][1][k]      - c13 * m.PresGradForce.x[i][2][k];
                m.PresGradForce.x[i][m.jm-1][k] = c43 * m.PresGradForce.x[i][m.jm-2][k] - c13 * m.PresGradForce.x[i][m.jm-3][k];
            }
        }

        // ----------------------------------------------------------------
        // 7c. k-direction (periodic longitude) BCs for all 3D arrays
        // ----------------------------------------------------------------
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                auto periodic_bc_3D = [&](Array& A) {
                    double v0 = c43 * A.x[i][j][1]      - c13 * A.x[i][j][2];
                    double v1 = c43 * A.x[i][j][m.km-2] - c13 * A.x[i][j][m.km-3];
                    A.x[i][j][0] = A.x[i][j][m.km-1] = (v0 + v1) * 0.5;
                };
                periodic_bc_3D(m.BuoyancyForce);
                periodic_bc_3D(m.CoriolisForce);
                periodic_bc_3D(m.PresGradForce);
            }
        }

        auto end     = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for RunDataHyd\n", elapsed.count() * 1e-9);
        cout << "      RunDataHyd ended" << endl;
    }

private:
    cHydrosphereModel& m;
};
