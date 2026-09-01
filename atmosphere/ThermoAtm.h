#pragma once

#include "cAtmosphereModel.h"
#include "Utils.h"

#include <vector>
#include <cmath>
#include <chrono>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace AtomUtils;

// Physical constants for precipitable water column integration
namespace PrecipWaterConstants {
    constexpr double HPA_TO_PA        = 100.0;
    constexpr double MIN_SAFE_TEMP    = 100.0;   // [K]
    constexpr double MIN_SAFE_PRESSURE = 1e-6;   // [hPa]
}

class ThermoAtm {
public:
    explicit ThermoAtm(cAtmosphereModel& model)
        : m(model)
    {}

    // ------------------------------------------------------------------
    void latentSensibleHeat()
    {
        using namespace std;
        cout << endl << endl << endl << "      LatentSensibleHeat" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        // Q_SENSIBLE WAS DELETED HERE (2026-08-30), AND THIS FUNCTION IS WHAT IS LEFT OF IT.
        //
        // The array was written twice with two different formulas and read by nothing:
        // RHS_Atm_Turb.cpp had `coeff_S*lap(T)`, a conductive flux DIVERGENCE, and this block
        // had `cp*T_0/sqrt(3)*|grad T|`, an unsigned gradient MAGNITUDE. Neither is a sensible
        // heat flux; the surface flux is `c_H*(T_s - T_air1)` in MultiLayerRadiation, which is
        // debited from the surface balance and credited to nothing. UtilsAtm then clamped the
        // array to >= 0, discarding half of whichever quantity had last been written, and
        // Results_Atm printed its min/max in kJ/kg as though it meant something. Wiring it to
        // the temperature equation would have double-counted `diffusion_t_re`.
        //
        // What remains here is the Q_Latent land-zeroing this loop also did. It is kept rather
        // than deleted because bcSolidGround's identical zeroing runs only AFTER bcRadius /
        // bcTheta / bcPhi, whose extrapolation lists write Q_Latent at the boundaries -- so the
        // two are not trivially redundant and proving it is a separate job.

        #pragma omp parallel for collapse(2)
        for (int j = 0; j < m.jm; j++)
            for (int k = 0; k < m.km; k++)
                m.Q_Latent.x[0][j][k] = 0.0;

        #pragma omp parallel for collapse(2) schedule(dynamic, 4)
        for (int j = 1; j < m.jm-1; j++)
            for (int k = 1; k < m.km-1; k++)
                for (int i = 1; i < m.im-1; i++)
                    if (is_land(m.h, i, j, k))
                        m.Q_Latent.x[i][j][k] = 0.0;

        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for LatentSensibleHeat\n", elapsed.count() * 1e-9);
        cout << "      LatentSensibleHeat ended" << endl;
    }

    // ------------------------------------------------------------------
    void waterVapourEvaporation()
    {
        using namespace std;
        cout << "\n\n\n      WaterVapourEvaporation  (model: " << m.evap_model << ")" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        const double conv_factor  = 8.64e4;                             // [mm/s] -> [mm/d]
        const double hPa_to_mmHg  = 0.750062;                           // [mmHg/hPa]
        const double ms_to_kmh    = 3.6;                                // [(km/h)/(m/s)]
        // Meyer (1915): E[mm/month] = C·(1 + W/16)·(e_s − e_a), with the deficit in
        // mmHg and C the open-water coefficient ≈ 11 (deep/large water, e.g. ocean) or
        // 15 (shallow ponds).  C is unit-identical for E in mm/month & deficit in mmHg
        // because inch→mm and inHg→mmHg share the 25.4 factor.  The previous value 0.36
        // was ~30× too small, giving ~0.18 mm/d (≈66 mm/yr) instead of the realistic
        // ~5 mm/d (≈1800 mm/yr) and far below the Dalton/Rohwer diagnostics.
        const double K_Meyer      = 11.0;                               // Meyer (1915) deep open-water coefficient [mm/month/mmHg]

        // Vertical spread: distribute c_eq over n_spread+1 levels with exp(-i) weights,
        // normalised so that sum_{i=0}^{n_spread} exp(-i) = 1 (moisture conserved).
//        const int    n_spread = 5;
        const int    n_spread = 3;
//        const int    n_spread = 2;
        const double r        = std::exp(-1.0);
        const double w_norm   = (1.0 - r) / (1.0 - std::pow(r, n_spread + 1)); // 1/sum

        // Determine which formula drives the surface-humidity (c) update.
        // All three are always computed for diagnostic output.
        enum class EvapModel { Dalton, Meyer, Rohwer } active;

        if      (m.evap_model == "Meyer")  active = EvapModel::Meyer;
        else if (m.evap_model == "Rohwer") active = EvapModel::Rohwer;
        else                               active = EvapModel::Dalton;

        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {

                if (is_land(m.h, 0, j, k)) {
                    m.Evaporation_Dalton.y[j][k] = 0.0;
                    m.Evaporation_Meyer.y[j][k]  = 0.0;
                    m.Evaporation_Rohwer.y[j][k] = 0.0;
                    m.Evaporation.y[j][k]        = 0.0;
                    continue;
                }

                double c_Dalton    = AtomUtils::C_Dalton(0, j, k, m.coeff_Dalton, m.u_0, m.u, m.v, m.w);
                                                                        // [mm/(h*hPa)]
                double p_stat_0jk  = m.p_stat.x[0][j][k];               // [hPa]
                // Floor the surface temperature at 180 K (below Earth's coldest-ever ~184 K, so
                // it never touches real physics) to backstop the marginal Gulf-of-Alaska crash:
                // an undamped i=0 coastal-surface T oscillation drove t.x[0] to ~160 K, making
                // E_sat/c_eq in the evaporation formula singular → c NaN. See
                // [[project_upper_velocity_secular_growth]].
                double t_u_base    = std::max(180.0, m.t.x[0][j][k] * m.t_0);  // [K], floored
                double precip_term = conv_factor * m.Precipitation.x[0][j][k]; // [mm/d]

                double vel_ms = sqrt((m.u.x[0][j][k] * m.u.x[0][j][k]
                                    + m.v.x[0][j][k] * m.v.x[0][j][k]
                                    + m.w.x[0][j][k] * m.w.x[0][j][k]) / 3.0) * m.u_0; // [m/s]
                double u_kmh  = vel_ms * ms_to_kmh;                     // [km/h]
                double p_mmHg = p_stat_0jk * hPa_to_mmHg;               // [mmHg]

                double E_sat  = (t_u_base >= m.t_0)                     // [hPa]
                    ? m.hp * AtomUtils::exp_func(t_u_base, 17.2694, 35.86)
                    : m.hp * AtomUtils::exp_func(t_u_base, 21.8746,  7.66);

                // Per-formula coefficients [mm/(d*hPa)]: E = coeff * sat_deficit_hPa
                //   Dalton: from wind-dependent mass-transfer coefficient
                //   Meyer (1915): K_M*(1+u/16)/30 [mm/month/mmHg -> mm/d/hPa]
                //   Rohwer (1931): 0.771*(1.465-0.000732*p)*(0.44+0.0733*u) [mm/d/mmHg -> mm/d/hPa]
                double coeff_D = c_Dalton * 24.0;
                double coeff_M = K_Meyer * hPa_to_mmHg * (1.0 + u_kmh / 16.0) / 30.0;
                double coeff_R = 0.771 * (1.465 - 0.000732 * p_mmHg)
                               * (0.44  + 0.0733   * u_kmh) * hPa_to_mmHg;

                // Calm conditions (zero wind): skip c update, Dalton = 0.
                // Meyer and Rohwer retain their still-air (u=0) terms.
                if (c_Dalton <= 0.0) {
                    double e  = m.c.x[0][j][k] * p_stat_0jk / m.ep;     // [hPa]  vapour pressure at c_Dalton = 0
                    double sd = std::max(0.0, E_sat - e);               // [hPa]  saturation deficit at c_Dalton = 0
                    m.Evaporation_Dalton.y[j][k] = 0.0;
                    m.Evaporation_Meyer.y[j][k]  = coeff_M * sd;        // [mm/d]
                    m.Evaporation_Rohwer.y[j][k] = coeff_R * sd;        // [mm/d]
                    m.Evaporation.y[j][k] = (active == EvapModel::Meyer)  ? m.Evaporation_Meyer.y[j][k]
                                          : (active == EvapModel::Rohwer) ? m.Evaporation_Rohwer.y[j][k]
                                          :                                 m.Evaporation_Dalton.y[j][k];
                    // Calm ocean cell: the wind-driven evaporation FLUX is zero (above), but the
                    // surface-layer vapour CONCENTRATION still equilibrates toward saturation (a
                    // windy no-precip cell converges to ~c_sat too). Relax c.x[0] toward c_sat with
                    // the same w_norm as the windy branch below, instead of leaving the
                    // RK4-unintegrated i=0 layer at a stale near-zero value — which showed as
                    // near-zero water-vapour patches over calm ocean (with normal "spots" at windy
                    // cells) in zonal cross-sections. Concentration only; the flux stays wind-limited.
                    double denom_calm = p_stat_0jk - (1.0 - m.ep) * E_sat;      // [hPa]
                    double c_sat_calm = (denom_calm > 0.0) ? m.ep * E_sat / denom_calm
                                                           : m.ep * E_sat / p_stat_0jk; // [kg/kg]
                    m.c_fix.y[j][k] = m.c.x[0][j][k];
                    m.c.x[0][j][k]  = m.c_fix.y[j][k] + (c_sat_calm - m.c_fix.y[j][k]) * w_norm;
                    continue;
                }

                // Saturation specific humidity (exact, not dilute approximation == Verdünnungsnäherung)
                double denom_guard = p_stat_0jk - (1.0 - m.ep) * E_sat; // [hPa]  guards that static pressure is greater than reduced current saturation pressure
                double c_sat       = (denom_guard > 0.0)
                    ? m.ep * E_sat / denom_guard                        // [kg/kg]
                    : m.ep * E_sat / p_stat_0jk;                        // dilute fallback

                // Active formula solves for c_eq balancing precipitation:
                //   E_active = coeff_active * (E_sat - c_eq*p/ep) = precip_term
                //   => c_eq = ep/p * (E_sat - precip_term / coeff_active)
                // If c_eq < 0, heavy rain saturates air -> clamp to c_sat.
                double coeff_active = (active == EvapModel::Meyer)  ?   coeff_M
                                    : (active == EvapModel::Rohwer) ?   coeff_R
                                    :                                   coeff_D;

                double c_eq = m.ep * (E_sat - precip_term / coeff_active) / p_stat_0jk;  // equilibrium humidity (dilute-approx saturation)
                c_eq = (c_eq < 0.0) ? c_sat : std::min(c_eq, c_sat);

                m.c_fix.y[j][k] = m.c.x[0][j][k];

                // Evaporations from current (pre-update) saturation deficit
                double e_cur  = m.c_fix.y[j][k] * p_stat_0jk / m.ep;    // [hPa]  current vapour pressure
                double sd_cur = std::max(0.0, E_sat - e_cur);           // [hPa]  current saturation deficit
 
                m.Evaporation_Dalton.y[j][k] = coeff_D * sd_cur;        // [mm/d]
                m.Evaporation_Meyer.y[j][k]  = coeff_M * sd_cur;        // [mm/d]
                m.Evaporation_Rohwer.y[j][k] = coeff_R * sd_cur;        // [mm/d]
                m.Evaporation.y[j][k] = (active == EvapModel::Meyer)  ? m.Evaporation_Meyer.y[j][k]
                                      : (active == EvapModel::Rohwer) ? m.Evaporation_Rohwer.y[j][k]
                                      :                                 m.Evaporation_Dalton.y[j][k];

                // c update: distribute c_eq across surface + n_spread levels above,
                // conserving total moisture (weights sum to 1).
                m.c.x[0][j][k] = m.c_fix.y[j][k] + (c_eq - m.c_fix.y[j][k]) * w_norm;  // i=0: weight = exp(0)*w_norm

                if (c_eq > 0.0) {
                    for (int i = 1; i <= n_spread; i++) {
                        double weight  = std::pow(r, i) * w_norm;
                        double t_i     = m.t.x[i][j][k] * m.t_0;
                        double p_i     = m.p_stat.x[i][j][k];
                        double E_i     = (t_i >= m.t_0)
                            ? m.hp * AtomUtils::exp_func(t_i, 17.2694, 35.86)
                            : m.hp * AtomUtils::exp_func(t_i, 21.8746,  7.66);
                        double denom_i = p_i - (1.0 - m.ep) * E_i;
                        double c_sat_i = (denom_i > 0.0) ? m.ep * E_i / denom_i : m.ep * E_i / p_i;
                        m.c.x[i][j][k] = std::min(m.c.x[i][j][k] + c_eq * weight, c_sat_i);
                    }
                }


/*
                cout.precision(8);
                cout.setf(ios::fixed);
                if ((j == 90) && (k == 180)) cout << endl
                    << "  WaterVapourEvaporation" << endl
                    << "  j = " << j << "  k = " << k << endl
                    << "  p_stat_0jk = "     << p_stat_0jk
                    << "  p_mmHg = "         << p_mmHg << endl
                    << "  c_fix = "          << m.c_fix.y[j][k] * 1e3
                    << "  c_eq = "           << c_eq * 1e3
                    << "  c_sat = "          << c_sat * 1e3
                    << "  c = "              << m.c.x[0][j][k] * 1e3 << endl
                    << "  E_sat = "          << E_sat
                    << "  e_cur = "          << e_cur
                    << "  sd_cur = "         << sd_cur << endl
                    << "  c_Dalton = "       << c_Dalton
                    << "  coeff_D = "        << coeff_D
                    << "  coeff_M = "        << coeff_M
                    << "  coeff_R = "        << coeff_R << endl
                    << "  Evap_Dalton = "    << m.Evaporation_Dalton.y[j][k]
                    << "  Evap_Meyer = "     << m.Evaporation_Meyer.y[j][k]
                    << "  Evap_Rohwer = "    << m.Evaporation_Rohwer.y[j][k]
                    << "  Evap = "           << m.Evaporation.y[j][k]
                    << "  Prec = "           << precip_term << endl;
*/
            }  // k
        }  // j

        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for WaterVapourEvaporation\n", elapsed.count() * 1e-9);
        cout << "      WaterVapourEvaporation ended" << endl;
    }

    // ------------------------------------------------------------------
    void standAtm_DewPoint_HumidRel()
    {
        using namespace std;
        cout << endl << endl << endl << "      StandAtm_DewPoint_HumidRel" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        std::vector<double> height_table(m.im);
        for (int i = 0; i < m.im; i++)
            height_table[i] = m.get_layer_height(i);

        const double lapse_rate = 6.5e-3;
        const double inv_ep     = 1.0 / m.ep;

        std::vector<int> i_trop_table(m.jm);
        for (int j = 0; j < m.jm; j++) {
            i_trop_table[j] = (j == 90)
                ? (int)m.tropopause_layers[j-1]
                : (int)m.tropopause_layers[j];
        }

        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {

                int i_trop = i_trop_table[j];

                double TempStand_surface   = m.t.x[0][j][k] * m.t_0 - m.t_0;
                m.TempStand.x[0][j][k]    = TempStand_surface;

                for (int i = 0; i <= i_trop; i++) {

                    double t_u = m.t.x[i][j][k] * m.t_0;

                    double E = (t_u > m.t_0)
                        ? m.hp * AtomUtils::exp_func(t_u, 17.2694, 35.86)
                        : m.hp * AtomUtils::exp_func(t_u, 21.8746,  7.66);

                    double e           = m.c.x[i][j][k] * m.p_stat.x[i][j][k] * inv_ep;
                    bool   zero_vapour = (e == 0.0);
                    if (zero_vapour) e = 1e-3;

                    m.TempStand.x[i][j][k]    = TempStand_surface - lapse_rate * height_table[i];

                    double L    = std::log(e / m.hp);
                    double a_dp = (t_u > m.t_0) ? 21.8746 : 17.2694;
                    double b_dp = (t_u > m.t_0) ?  7.66   : 35.86;
                    m.TempDewPoint.x[i][j][k] = L * (m.t_0 - b_dp) / (a_dp - L);

                    m.HumidityRel.x[i][j][k]  = zero_vapour
                        ? 0.0
                        : std::min(e / E * 100.0, 100.0);

                }  // i

                // Above tropopause: copy tropopause values
                double ts_trop = m.TempStand.x[i_trop][j][k];
                double td_trop = m.TempDewPoint.x[i_trop][j][k];
                double hr_trop = m.HumidityRel.x[i_trop][j][k];

                for (int i = i_trop + 1; i < m.im; i++) {
                    m.TempStand.x[i][j][k]    = ts_trop;
                    m.TempDewPoint.x[i][j][k] = td_trop;
                    m.HumidityRel.x[i][j][k]  = hr_trop;
                }

                // Below-terrain fill: copy from topography surface downward
                int    i_mount   = m.i_topography[j][k];
                double ts_mount  = m.TempStand.x[i_mount][j][k];
                double td_mount  = m.TempDewPoint.x[i_mount][j][k];
                double hr_mount  = m.HumidityRel.x[i_mount][j][k];

                for (int i = i_mount - 1; i >= 0; i--) {
                    if (is_land(m.h, i, j, k)) {
                        m.TempStand.x[i][j][k]    = ts_mount;
                        m.TempDewPoint.x[i][j][k] = td_mount;
                        m.HumidityRel.x[i][j][k]  = hr_mount;
                    }
                }
            }  // k
        }  // j

        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for StandAtm_DewPoint_HumidRel\n", elapsed.count() * 1e-9);
        cout << "      StandAtm_DewPoint_HumidRel ended" << endl;
    }

    // ------------------------------------------------------------------
    void forces()
    {
        using namespace std;
        cout << endl << endl << endl << "      ATOM: Forces" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        const double inv_2dr    = 1.0 / (2.0 * m.dr);
        const double inv_2dthe  = 1.0 / (2.0 * m.dthe);
        const double inv_2dphi  = 1.0 / (2.0 * m.dphi);
        const double two_omega  = 2.0 * m.omega;
        const double omega2     = m.omega * m.omega;
        const double pres_coeff = 1e2 * m.p_0 / m.L_atm;
        const double r_Earth_m  = m.r_Earth * 1e3;

        std::vector<double> sinthe_table(m.jm);
        std::vector<double> costhe_table(m.jm);
        for (int j = 0; j < m.jm; j++) {
            sinthe_table[j] = sin(m.the.z[j]);
            costhe_table[j] = cos(m.the.z[j]);
        }

        #pragma omp parallel for collapse(2)
        for (int i = 1; i < m.im-1; i++) {
            for (int j = 1; j < m.jm-1; j++) {

                double rm           = m.rad.z[i];
                double exp_rm       = m.metricExpRm(rm);
                double inv_rm       = 1.0 / rm;
                double sinthe       = sinthe_table[j];
                double costhe       = costhe_table[j];
                double inv_rmsinthe = 1.0 / (rm * sinthe);
                double abs_sinthe   = fabs(sinthe);

                double rad_dist  = (double)i * m.L_atm * exp_rm;
                double rad_Earth = rad_dist + r_Earth_m;

                for (int k = 1; k < m.km-1; k++) {

                    double w_ijk   = m.w.x[i][j][k] * m.u_0;
                    double Cor_r   = -two_omega * sinthe * w_ijk;
                    double Cor_the =  two_omega * costhe * w_ijk;
                    double Cor_phi =  two_omega * (-costhe * m.v.x[i][j][k]
                                                   + sinthe * m.u.x[i][j][k]) * m.u_0;

                    m.CoriolisForce.x[i][j][k] = m.Coriolis * m.r_air
                        * sqrt(Cor_r*Cor_r + Cor_the*Cor_the + Cor_phi*Cor_phi);

                    m.CentrifugalForce.x[i][j][k] =
                        m.centrifugal * m.r_air * omega2 * rad_Earth * (1.0 + abs_sinthe);

                    // Diagnostic buoyancy force for ParaView/Results — must match the
                    // perturbation-form body force applied in RHS_Atm.cpp (rhs_u), i.e.
                    // +g·ρ·(t − t̄(i)) referenced to the per-level mean t_ref_level[i].
                    const double t_ref_b = ((int)m.t_ref_level.size() == m.im)
                                           ? m.t_ref_level[i] : 1.0;
                    m.BuoyancyForce.x[i][j][k] = 1.0e-3 * m.buoyancy * m.r_humid.x[i][j][k] * m.g
                                                 * (m.t.x[i][j][k] - t_ref_b);

                    double dpdr   = (m.p_dyn.x[i+1][j][k] - m.p_dyn.x[i-1][j][k])
                                    * inv_2dr * exp_rm;
                    double dpdthe = (m.p_dyn.x[i][j+1][k] - m.p_dyn.x[i][j-1][k])
                                    * inv_2dthe * inv_rm;
                    double dpdphi = (m.p_dyn.x[i][j][k+1] - m.p_dyn.x[i][j][k-1])
                                    * inv_2dphi * inv_rmsinthe;

                    m.PresGradForce.x[i][j][k] =
                        - 1.0e-3 * sqrt(dpdr*dpdr + dpdthe*dpdthe + dpdphi*dpdphi) * pres_coeff;
                }  // k
            }  // j
        }  // i

        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for Forces\n", elapsed.count() * 1e-9);
        cout << "      ATOM: Forces ended" << endl;
    }

    // ------------------------------------------------------------------
    void adjustTemperatureIC(double** t, int jm, int km)
    {
        const double inv_t0 = 1.0 / m.t_0;
        const int    k_half = (km - 1) / 2;

        #pragma omp parallel for
        for (int j = 0; j < jm; j++) {
            for (int k = 0; k < km; k++)
                t[j][k] = (t[j][k] + m.t_0) * inv_t0;
            t[j][k_half] = (t[j][k_half + 1] + t[j][k_half - 1]) * 0.5;
            m.temperature_NASA.y[j][k_half] =
                (m.temperature_NASA.y[j][k_half + 1] +
                 m.temperature_NASA.y[j][k_half - 1]) * 0.5;
        }
    }

    // ------------------------------------------------------------------
    void precipitableWater()
    {
        using namespace std;
        using namespace PrecipWaterConstants;

        cout << "\n\n\n      PrecipitableWater" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        #pragma omp parallel for collapse(2)
        for (int j = 0; j < m.jm; j++)
            for (int k = 0; k < m.km; k++)
                m.precipitable_water.y[j][k] = 0.0;

        #pragma omp parallel for collapse(2)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {

                double column_sum = 0.0;

                for (int i = 0; i < m.im - 1; i++) {

                    const double t_actual  = m.t.x[i][j][k] * m.t_0;
                    const double p_actual  = m.p_stat.x[i][j][k];
                    const double q_mixing  = m.c.x[i][j][k];

                    if (t_actual < MIN_SAFE_TEMP || p_actual < MIN_SAFE_PRESSURE) {
                        m.PrecipitableWaterLocal.x[i][j][k] = 0.0;
                        continue;
                    }

                    const double e          = HPA_TO_PA * q_mixing * p_actual / m.ep;
                    const double a          = e / (m.R_WaterVapour * t_actual);
                    const double step       = m.get_layer_height(i + 1) - m.get_layer_height(i);
                    const double local_mass = a * step;

                    m.PrecipitableWaterLocal.x[i][j][k] = local_mass;
                    column_sum += local_mass / m.r_0_water;
                }

                m.precipitable_water.y[j][k] = column_sum * 1e3;   // m -> mm
            }  // k
        }  // j

        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" Time measured: %.3f seconds for PrecipitableWater\n", elapsed.count() * 1e-9);
        cout << "      PrecipitableWater ended" << endl;
    }

    // ------------------------------------------------------------------
    void vegetationLand()
    {
        using namespace std;
        cout << "\n\n\n      VegetationLand" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        constexpr double max_vegetation_height_m = 4400.0;
        constexpr double min_vegetation_temp_c   = -40.0;

        std::vector<double> height_table(m.im);
        for (int i = 0; i < m.im; i++)
            height_table[i] = m.get_layer_height(i);

        double max_Precipitation = 0.0;

        #pragma omp parallel for collapse(2) reduction(max:max_Precipitation)
        for (int j = 0; j < m.jm; j++)
            for (int k = 0; k < m.km; k++)
                max_Precipitation = std::max(max_Precipitation, m.Precipitation.x[0][j][k]);

        double inv_max_Precipitation = 1.0 / std::max(max_Precipitation, 1e-5);

        #pragma omp parallel for collapse(2)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                int i_mount = m.i_topography[j][k];
                if (is_land(m.h, 0, j, k)
                    && height_table[i_mount] < max_vegetation_height_m
                    && m.t.x[i_mount][j][k] * m.t_0 - m.t_0 > min_vegetation_temp_c) {
                    m.Vegetation.y[j][k] = m.Precipitation.x[0][j][k] * inv_max_Precipitation;
                } else {
                    m.Vegetation.y[j][k] = 0.0;
                }
            }
        }

        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        cout << " time measured: " << fixed << setprecision(3)
             << elapsed.count() * 1e-9 << " seconds for VegetationLand" << endl;
        cout << "      VegetationLand ended" << endl;
    }

    // ------------------------------------------------------------------
    void printDataAtm()
    {
        using namespace std;
        cout << endl << endl << endl << "      print_data_atm" << endl;
        cout << endl << " properties of the atmosphere at the surface: " << endl;

        double temperature_NASA_average, temperature_average, temperature_expected_average,
               precipitablewater_average, precipitation_average, precipitation_NASA_average,
               co2_average, Evaporation_average, Evaporation_Dalton_average,
               Evaporation_Meyer_average, Evaporation_Rohwer_average;

        #pragma omp parallel sections
        {
            #pragma omp section
                { temperature_NASA_average    = AtomUtils::GetMean_2D(m.jm, m.km, m.temperature_NASA); }
            #pragma omp section
                { temperature_average         = (AtomUtils::GetMean_3D(m.jm, m.km, m.t) - 1.0) * m.t_0; }
            #pragma omp section
                { precipitablewater_average   = AtomUtils::GetMean_2D(m.jm, m.km, m.precipitable_water); }
            #pragma omp section
                { precipitation_average       = 365.0 * AtomUtils::GetMean_3D(m.jm, m.km, m.Precipitation); }
            #pragma omp section
                { precipitation_NASA_average  = 365.0 * AtomUtils::GetMean_2D(m.jm, m.km, m.precipitation_NASA); }
            #pragma omp section
                { co2_average                 = AtomUtils::GetMean_2D(m.jm, m.km, m.co2_total); }
            #pragma omp section
                { Evaporation_Dalton_average  = 365.0 * AtomUtils::GetMean_2D(m.jm, m.km, m.Evaporation_Dalton); }
            #pragma omp section
                { Evaporation_Meyer_average   = 365.0 * AtomUtils::GetMean_2D(m.jm, m.km, m.Evaporation_Meyer); }
            #pragma omp section
                { Evaporation_Rohwer_average  = 365.0 * AtomUtils::GetMean_2D(m.jm, m.km, m.Evaporation_Rohwer); }
            #pragma omp section
                { Evaporation_average         = 365.0 * AtomUtils::GetMean_2D(m.jm, m.km, m.Evaporation); }
        }
        temperature_expected_average =
            m.get_temperatures_from_curve(*m.get_current_time(), m.m_global_temperature_curve);

        cout.precision(2);
        cout << endl << endl;

        auto row = [](const char* n1, double v1, const char* u1,
                      const char* n2, double v2, const char* u2,
                      const char* n3, double v3, const char* u3) {
            cout << setiosflags(ios::left) << setw(40) << setfill('.')
                 << n1 << " = " << resetiosflags(ios::left) << setw(7)
                 << fixed << setfill(' ') << v1 << setw(6) << u1 << "   "
                 << setiosflags(ios::left) << setw(40) << setfill('.')
                 << n2 << " = " << resetiosflags(ios::left) << setw(7)
                 << fixed << setfill(' ') << v2 << setw(6) << u2 << "   "
                 << setiosflags(ios::left) << setw(40) << setfill('.')
                 << n3 << " = " << resetiosflags(ios::left) << setw(7)
                 << fixed << setfill(' ') << v3 << setw(6) << u3 << endl;
        };

        double precip_mm_a = precipitation_average * 8.64e4;

        row(" precipitable water average", precipitablewater_average, " mm",
            " precipitation average per year", precip_mm_a, " mm/a",
            " precipitation average per day", precip_mm_a / 365.0, " mm/d");

        // Per-component breakdown: which of P_rain / P_snow / P_graupel / P_conv
        // is climbing? Same units as the total: mean × 365·86400 gives mm/a.
        // Max scan reports the worst single cell with its (i,j,k) so a runaway
        // can be localised geographically.
        {
            const double s_per_year = 365.0 * 8.64e4;
            auto component = [&](const char* name, const Array& F){
                double mean_v = AtomUtils::GetMean_3D(m.jm, m.km, const_cast<Array&>(F));
                double mx = 0.0; int mi = 0, mj = 0, mk = 0;
                for(int i = 0; i < m.im; i++)
                    for(int j = 0; j < m.jm; j++)
                        for(int k = 0; k < m.km; k++){
                            double v = F.x[i][j][k];
                            if(v > mx){ mx = v; mi = i; mj = j; mk = k; }
                        }
                cout << "    " << setw(8) << setfill(' ') << name
                     << " mean = " << scientific << setprecision(3) << mean_v * s_per_year
                     << " mm/a   max = " << mx * s_per_year
                     << " mm/a  @(i=" << mi << ",j=" << mj
                     << ",k=" << mk << ")" << fixed << endl;
            };
            const double mean_precip_mm_a =
                AtomUtils::GetMean_3D(m.jm, m.km, m.Precipitation) * s_per_year;
            cout << " precipitation by component (mm/a, surface-equivalent):" << endl;
            component("Precip",    m.Precipitation);   // sanity check vs model total above
            component("P_rain",    m.P_rain);
            component("P_snow",    m.P_snow);
            component("P_graupel", m.P_graupel);
            component("P_conv",    m.P_conv);

            // WHAT THE FLOOR MANUFACTURED, printed beside the total it belongs to.
            //
            // `P_x[i] = max(0, P_x[i+1] + dP)` is not a guard, it is a SOURCE: wherever the
            // sinks at a level exceed the flux present it creates the difference out of nothing,
            // and the flux continues down the column as if the water had been there. On the
            // accepted configuration it manufactures 8116 mm/a -- eight times NASA's entire
            // precipitation -- and the total above is what survived of it. Enforcing mass
            // conservation (ATM_ICE_LIMIT_ARRIVING, TwoCat) takes the same configuration from
            // 992 mm/a to 369.
            //
            // This is printed unconditionally because the reported total was quoted as agreeing
            // with NASA to 7 % for a week before anyone asked what it was made of.
            const auto& fa = IceSchemeCommon::floorAudit();
            if(fa.valid){
                const double tot = mean_precip_mm_a;
                cout << "    " << setw(8) << setfill(' ') << "floor"
                     << " injected " << scientific << setprecision(3) << fa.injected_mm_a
                     << " mm/a of water no source produced";
                if(tot > 0.0)
                    cout << "  = " << fixed << setprecision(1)
                         << 1e2 * fa.injected_mm_a / tot << " % of the reported total";
                cout << fixed << endl;
            }

            // ================= SCORED AGAINST THE NASA FIELD, NOT ITS MEAN =================
            //
            // `precipitation_NASA` has been read into memory, written to VTK and printed as a
            // GLOBAL MEAN beside the model's since the beginning, and the two fields have never
            // been differenced. A global mean is one scalar against a model with several free
            // constants per module: this tree has matched it twice -- 1046 shipped and 992
            // accepted -- and BOTH matches turned out to be a floor injecting eight times NASA's
            // precipitation. A matched mean is not agreement, and nothing printed until now
            // could tell the difference.
            //
            // Cos-lat weighted throughout, on the same weights GetMean_3D uses, so the model
            // mean below is the `Precip` mean above. NASA is mm/d in the file; x365 -> mm/a.
            {
                const double yr = 365.0;
                double w_t = 0.0, m_s = 0.0, n_s = 0.0;
                double mm = 0.0, nn = 0.0, mn = 0.0, se = 0.0;
                // 0-15, 15-35, 35-65, 65-90 degrees of |latitude|, then land / ocean
                double bw[4] = {0,0,0,0}, bm[4] = {0,0,0,0}, bn[4] = {0,0,0,0};
                double lw = 0.0, lm = 0.0, ln_ = 0.0, ow = 0.0, om = 0.0, on = 0.0;

                for(int j = 0; j < m.jm; j++){
                    const double w = cos((j / (double)(m.jm - 1) - 0.5) * M_PI);
                    const double alat = fabs(90.0 - j * 180.0 / (double)(m.jm - 1));
                    const int b = (alat < 15.0) ? 0 : (alat < 35.0) ? 1 : (alat < 65.0) ? 2 : 3;
                    for(int k = 0; k < m.km; k++){
                        const double mv = m.Precipitation.x[0][j][k] * s_per_year;
                        const double nv = m.precipitation_NASA.y[j][k] * yr;
                        w_t += w;  m_s += w * mv;  n_s += w * nv;
                        bw[b] += w;  bm[b] += w * mv;  bn[b] += w * nv;
                        if(is_land(m.h, 0, j, k)){ lw += w; lm += w * mv; ln_ += w * nv; }
                        else                     { ow += w; om += w * mv; on  += w * nv; }
                    }
                }
                const double mbar = m_s / w_t, nbar = n_s / w_t;
                for(int j = 0; j < m.jm; j++){
                    const double w = cos((j / (double)(m.jm - 1) - 0.5) * M_PI);
                    for(int k = 0; k < m.km; k++){
                        const double dm = m.Precipitation.x[0][j][k] * s_per_year - mbar;
                        const double dn = m.precipitation_NASA.y[j][k] * yr - nbar;
                        mm += w * dm * dm;  nn += w * dn * dn;  mn += w * dm * dn;
                        const double d = dm - dn;
                        se += w * d * d;
                    }
                }
                const double sm = sqrt(mm / w_t), sn = sqrt(nn / w_t);
                const double r  = (sm > 0.0 && sn > 0.0) ? (mn / w_t) / (sm * sn) : 0.0;
                const double rmse = sqrt(se / w_t + (mbar - nbar) * (mbar - nbar));

                cout << " precipitation scored against SurfacePrecipitation_NASA.xyz"
                        " (cos-lat weighted, surface, mm/a):" << endl;
                printf("      model %8.1f   NASA %8.1f   bias %+8.1f (%+.1f %%)"
                       "   pattern r = %+.3f   centred RMS %7.1f   sigma model/NASA = %.2f\n",
                       mbar, nbar, mbar - nbar,
                       nbar > 0.0 ? 1e2 * (mbar - nbar) / nbar : 0.0,
                       r, rmse, sn > 0.0 ? sm / sn : 0.0);
                printf("      by |latitude|   0-15 %7.1f /%7.1f   15-35 %7.1f /%7.1f"
                       "   35-65 %7.1f /%7.1f   65-90 %7.1f /%7.1f   (model / NASA)\n",
                       bw[0] > 0 ? bm[0]/bw[0] : 0.0, bw[0] > 0 ? bn[0]/bw[0] : 0.0,
                       bw[1] > 0 ? bm[1]/bw[1] : 0.0, bw[1] > 0 ? bn[1]/bw[1] : 0.0,
                       bw[2] > 0 ? bm[2]/bw[2] : 0.0, bw[2] > 0 ? bn[2]/bw[2] : 0.0,
                       bw[3] > 0 ? bm[3]/bw[3] : 0.0, bw[3] > 0 ? bn[3]/bw[3] : 0.0);
                printf("      land %8.1f /%8.1f      ocean %8.1f /%8.1f   (model / NASA)\n",
                       lw > 0 ? lm/lw : 0.0, lw > 0 ? ln_/lw : 0.0,
                       ow > 0 ? om/ow : 0.0, ow > 0 ? on/ow : 0.0);
            }
        }

        row(" precipitable water average", precipitablewater_average, " mm",
            " precipitation NASA average per year", precipitation_NASA_average, " mm/a",
            " precipitation NASA average per day", precipitation_NASA_average / 365.0, " mm/d");

        row(" precipitable water average", precipitablewater_average, " mm",
            " Evaporation_Dalton_average per year", Evaporation_Dalton_average, " mm/a",
            " Evaporation_Dalton_average per day", Evaporation_Dalton_average / 365.0, " mm/d");

        row(" co2_average", co2_average, " ppm",
            " Evaporation_average per year", Evaporation_Dalton_average, " mm/a",
            " Evaporation_average per day", Evaporation_Dalton_average / 365.0, " mm/d");

        row(" co2_average", co2_average, " ppm",
            " Evaporation_Meyer_average per year", Evaporation_Meyer_average, " mm/a",
            " Evaporation_Meyer_average per day", Evaporation_Meyer_average / 365.0, " mm/d");

        row(" co2_average", co2_average, " ppm",
            " Evaporation_Rohwer_average per year", Evaporation_Rohwer_average, " mm/a",
            " Evaporation_Rohwer_average per day", Evaporation_Rohwer_average / 365.0, " mm/d");

        cout << endl;

        cout << setiosflags(ios::left) << setw(40) << setfill('.')
             << " temperature_NASA_average" << " = " << resetiosflags(ios::left)
             << setw(7) << fixed << setfill(' ') << temperature_NASA_average
             << setw(6) << " deg" << "   "
             << setiosflags(ios::left) << setw(40) << setfill('.')
             << " temperature_average" << " = " << resetiosflags(ios::left)
             << setw(7) << fixed << setfill(' ') << temperature_average
             << setw(6) << " deg" << "   "
             << setiosflags(ios::left) << setw(40) << setfill('.')
             << " temperature_expected_average" << " = " << resetiosflags(ios::left)
             << setw(7) << fixed << setfill(' ') << temperature_expected_average
             << setw(6) << " deg" << endl << endl << endl;

        cout << "      print_data_atm ended" << endl;
    }

    // ------------------------------------------------------------------
    void co2Atmosphere()
    {
        using namespace std;
        cout << endl << endl << endl << "      AGCM: co2_atmosphere" << endl;

        double t_paleo_add = 0.0;

        // The paleo CO2 increment is defined relative to the preceding slice. On a
        // single-Ma run there is no foregoing slice to difference against (and
        // get_previous_time() would throw), so leave t_paleo_add = 0 in that case.
        if ((!m.use_NASA_temperature || *m.get_current_time() > 0) && !m.is_first_time_slice())
            t_paleo_add =
                m.get_temperatures_from_curve(*m.get_current_time(),  m.m_global_temperature_curve)
              - m.get_temperatures_from_curve(*m.get_previous_time(), m.m_global_temperature_curve);

        m.co2_paleo = t_paleo_add * (3.2886 * t_paleo_add
            + 6.5772 * m.t_equat_modern - 32.8859);

        double co2_average = 3.2886 * m.t_equat_modern * m.t_equat_modern
            - 32.8859 * m.t_equat_modern + 102.2148;

        cout.precision(3);

        const char* temperature_comment = "      temperature increase at paleo times: ";
        const char* temperature_gain    = " t increase";
        const char* temperature_unit    = "°C ";
        const char* co_comment          = "      co2 increase at paleo times: ";
        const char* co_gain             = " co2 increase";
        const char* co_modern           = "      mean co2 at modern times: ";
        const char* co_paleo_str        = "      mean co2 at paleo times: ";
        const char* co_average_str      = " co2 modern";
        const char* co_average_pal      = " co2 paleo";
        const char* co_unit             = "ppm ";

        cout << setiosflags(ios::left) << setw(55) << setfill('.')
            << temperature_comment << resetiosflags(ios::left) << setw(13)
            << temperature_gain << " = " << setw(7) << setfill(' ')
            << t_paleo_add << setw(5) << temperature_unit << endl
            << setiosflags(ios::left) << setw(55) << setfill('.')
            << co_comment << resetiosflags(ios::left) << setw(12) << co_gain << " = "
            << setw(7) << setfill(' ') << m.co2_paleo << setw(5) << co_unit
            << endl << setw(55) << setfill('.') << setiosflags(ios::left) << co_modern
            << resetiosflags(ios::left) << setw(13) << co_average_str << " = "
            << setw(7) << setfill(' ') << co2_average << setw(5) << co_unit
            << endl << setw(55) << setfill('.') << setiosflags(ios::left)
            << co_paleo_str << resetiosflags(ios::left) << setw(13) << co_average_pal
            << " = " << setw(7) << setfill(' ') << co2_average + m.co2_paleo
            << setw(5) << co_unit << endl;

        double co2_max   = 397.0;
        double inv_h_top = 1.0 / m.get_layer_height(m.im - 1);

        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                double co2_surf      = m.t.x[0][j][k] * co2_max + m.co2_paleo;
                m.co2.x[0][j][k]    = co2_surf;

                int    i_mount       = m.i_topography[j][k];
                double co2_at_mount  = 0.0;

                if (i_mount > 0) {
                    double x_mount  = m.get_layer_height(i_mount) * inv_h_top;
                    co2_at_mount    = parabola_interp(m.co2_tropopause, co2_surf, x_mount);
                }

                for (int i = 1; i < m.im; i++) {
                    if (i <= i_mount && is_land(m.h, i, j, k)) {
                        m.co2.x[i][j][k] = co2_at_mount;
                    } else {
                        double x = m.get_layer_height(i) * inv_h_top;
                        m.co2.x[i][j][k] = parabola_interp(m.co2_tropopause, co2_surf, x);
                    }
                }

                if (i_mount > 0 && is_land(m.h, 0, j, k))
                    m.co2.x[0][j][k] = co2_at_mount;
            }
        }

        // CO2 sensitivity multiplier: scale the whole (ppm) field uniformly. co2_scale=1 leaves
        // the field as built; co2_scale=2 is a clean CO2 doubling (modern base + paleo increment),
        // so the mode-3 radiation sees 2x the ppm and the perturbation MLR(2x field) - MLR(280 ppm)
        // carries the doubling forcing. Leaves the paleo-CO2 construction (co2_max, co2_paleo)
        // untouched — this is a dedicated experiment knob (default 1.0).
        if (m.co2_scale != 1.0) {
            #pragma omp parallel for collapse(2) schedule(static)
            for (int j = 0; j < m.jm; j++)
                for (int k = 0; k < m.km; k++)
                    for (int i = 0; i < m.im; i++)
                        m.co2.x[i][j][k] *= m.co2_scale;
        }
        cout << "      AGCM: co2_atmosphere  co2_scale = " << m.co2_scale << endl;
        cout << "      AGCM: co2_atmosphere ended" << endl;
    }

    // ------------------------------------------------------------------
    void densities()
    {
        using namespace std;
        cout << "\n\n\n      PressureDensity" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        const double R_W_R_A_m1       = m.R_WaterVapour / m.R_Air - 1.0;
        const double inv_R_Air        = 1.0 / m.R_Air;
        const double scale            = 1e2 * inv_R_Air;
        const double p_sl_factor      = 1e-2 * m.r_air * m.R_Air;
        const double beta             = 42.0;                           // K, COSMO
        const double inv_beta         = 1.0 / beta;
        const double two_beta_g_inv_R = 2.0 * beta * m.g * inv_R_Air;

        std::vector<double> height_table(m.im);
        for (int i = 0; i < m.im; i++)
            height_table[i] = m.get_layer_height(i);

        #pragma omp parallel for collapse(2) schedule(dynamic, 4)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {

                // Floor surface T at 180 K (below Earth's coldest-ever ~184 K): coeff ∝ 1/T²,
                // so an unphysical cold surface (the coastal/Pamir T oscillation) blows coeff up
                // and drives the barometric sqrt negative. Backstops the sqrt guard at the
                // source. See [[project_upper_velocity_secular_growth]].
                const double t_u_0        = std::max(180.0, m.t.x[0][j][k] * m.t_0);
                const double p_sl         = p_sl_factor * t_u_0;
                const double t_0_inv_beta = t_u_0 * inv_beta;
                const double coeff        = two_beta_g_inv_R / (t_u_0 * t_u_0);

                for (int i = 0; i < m.im; i++) {
                    const double h_i = height_table[i];
                    const double t_u = m.t.x[i][j][k] * m.t_0;
                    // Guard the COSMO barometric sqrt: coeff = 2βg/(R_Air·t_u_0²) ∝ 1/T², so an
                    // anomalously cold surface column drives (1 - coeff·h_i) negative aloft →
                    // sqrt(NaN) → p_i/r_humid NaN → whole-field blow-up. The init-time copy of
                    // this formula (InitValues_Atm.cpp:644) already clamps with max(0,…); this
                    // in-loop version had dropped it. See [[project_upper_velocity_secular_growth]].
                    const double p_i = p_sl * exp(-t_0_inv_beta
                        * (1.0 - sqrt(std::max(0.0, 1.0 - coeff * h_i))));  // COSMO barometric formula

                    m.p_stat.x[i][j][k]  = p_i;
                    m.r_dry.x[i][j][k]   = scale * p_i / t_u;
                    m.r_humid.x[i][j][k] = scale * p_i
                        / ((1.0 + R_W_R_A_m1 * m.c.x[i][j][k]
                            - m.cloud.x[i][j][k] - m.ice.x[i][j][k]) * t_u);
                }

                const int    i_m = m.i_topography[j][k];
                m.p_stat_landscape.y[j][k] = m.p_stat.x[i_m][j][k];

                if (i_m >= 0 && i_m < m.im) {
                    m.p_stat.x[0][j][k]  = m.p_stat.x[i_m][j][k];
                    m.r_dry.x[0][j][k]   = m.r_dry.x[i_m][j][k];
                    m.r_humid.x[0][j][k] = m.r_humid.x[i_m][j][k];
                }
            }
        }

        // ---- Anelastic base state (ported from ATHAD) -----------------------------------
        //
        // rho_bar(i), the cos-latitude weighted horizontal mean of r_humid, and
        // d ln(rho_bar)/d(rad.z). Needed by ATM_ANELASTIC in PressureSolverAtm.h; costs one
        // sweep and is computed unconditionally so the knob has no ordering dependence.
        //
        // Land is skipped: a sub-surface cell carries a barometric r_humid that is not air, and
        // averaging it in would bias the base state toward sea level exactly where the terrain
        // is highest.
        m.m_rho_base.assign(m.im, 0.0);
        m.m_dlnrho_dr.assign(m.im, 0.0);
        {
            std::vector<double> lnrho(m.im, 0.0);
            for (int i = 0; i < m.im; i++) {
                double num = 0.0, den = 0.0;
                for (int j = 0; j < m.jm; j++) {
                    const double coslat = cos((j / (double)(m.jm - 1) - 0.5) * M_PI);
                    for (int k = 0; k < m.km; k++) {
                        if (i < m.i_topography[j][k]) continue;          // inside terrain
                        const double rho = m.r_humid.x[i][j][k];
                        if (!(rho > 0.0) || !std::isfinite(rho)) continue;
                        num += coslat * rho;
                        den += coslat;
                    }
                }
                m.m_rho_base[i] = (den > 0.0) ? num / den : 0.0;
                lnrho[i] = (m.m_rho_base[i] > 0.0) ? log(m.m_rho_base[i]) : 0.0;
            }
            // rad.z is uniform, so the centred difference is the stencil the solver uses on
            // every other field; the ends get the one-sided form so the first and last interior
            // cells -- exactly where the projection's wall BC acts -- are not fed a half-step.
            const double inv_2dr_b = 1.0 / (2.0 * m.dr);
            const double inv_dr_b  = 1.0 / m.dr;
            for (int i = 1; i < m.im - 1; i++)
                m.m_dlnrho_dr[i] = (lnrho[i+1] - lnrho[i-1]) * inv_2dr_b;
            if (m.im > 1) {
                m.m_dlnrho_dr[0]        = (lnrho[1] - lnrho[0]) * inv_dr_b;
                m.m_dlnrho_dr[m.im - 1] = (lnrho[m.im-1] - lnrho[m.im-2]) * inv_dr_b;
            }
        }

        // ---- Brunt-Vaisala frequency squared -------------------------------------------
        //
        // N^2 = (g/theta) d(theta)/dz, potential temperature theta = T*(p0/p)^(R/cp).
        // Ported from ATHAD, where it was written to test a claim of neutral stratification
        // that had no instrument. Here it has a use that tree could not make of it: Earth's
        // troposphere is near 1e-4 s^-2 and the stratosphere near 4e-4, so this is an ABSOLUTE
        // check, not a comparison between two arms of the same model.
        //
        // kappa = R_Air/cp_l = 0.2855, the dry-air value. ATHAD builds it from the local
        // mixture because its cp varies ~2x across 300-1500 K; here the variation is small
        // enough that a constant is honest, and saying so is better than importing machinery
        // this tree does not have.
        //
        // dz FROM get_layer_height(), NOT from the core's exp_rm -- the point of the field is
        // that it is physical and can be checked against the numbers above, which it could not
        // be if it inherited the 23.2x metric error checkRadialMetric() reports. Centred
        // difference inside, one-sided at the ends.
        {
            const double kap = (m.cp_l > 0.0) ? m.R_Air / m.cp_l : 0.2855;
            #pragma omp parallel for collapse(2) schedule(static)
            for (int j = 0; j < m.jm; j++) {
                for (int k = 0; k < m.km; k++) {
                    // THE GROUND IS THE BOTTOM OF THIS FIELD, NOT LEVEL 0. Over topography the
                    // cells below i_topography are rock, and they are not empty: the barometric
                    // loop above writes a real p_stat into every one of them, then overwrites
                    // p_stat.x[0] with the GROUND value -- so p_stat is non-monotonic across the
                    // sub-surface column (level 0 lower than level 1). theta inherits that jump,
                    // and a centred difference at the first air level straddles the terrain and
                    // differences an air theta against a rock one. That is the source of the
                    // +-0.03 to 0.047 s^-2 boundary-layer extrema the port commit recorded as
                    // appearing "every time": they are the terrain, not the boundary layer.
                    //
                    // So: build theta from the ground up, make the ground level ONE-SIDED
                    // exactly as the lid already is, and fill the rock below with the ground
                    // value (IceSchemeCommon::fillTopography's convention, applied inline here
                    // because this loop already owns the column). Over ocean i_m = 0 and every
                    // value is bit-identical to before.
                    const int i_m = m.i_topography[j][k];
                    std::vector<double> theta(m.im, 0.0);
                    for (int i = i_m; i < m.im; i++) {
                        const double T_i = m.t.x[i][j][k] * m.t_0;
                        const double p_i = m.p_stat.x[i][j][k];
                        theta[i] = (T_i > 0.0 && p_i > 0.0) ? T_i * pow(m.p_0 / p_i, kap) : 0.0;
                    }
                    for (int i = i_m; i < m.im; i++) {
                        const int il = (i > i_m)       ? i - 1 : i;
                        const int iu = (i < m.im - 1)  ? i + 1 : i;
                        const double dz = m.get_layer_height(iu) - m.get_layer_height(il);
                        const double th = theta[i];
                        m.brunt_N2.x[i][j][k] = (dz > 0.0 && th > 0.0)
                                        ? (m.g / th) * (theta[iu] - theta[il]) / dz : 0.0;
                    }
                    for (int i = i_m - 1; i >= 0; i--)
                        m.brunt_N2.x[i][j][k] = m.brunt_N2.x[i_m][j][k];
                }
            }
        }

        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for PressureDensity\n", elapsed.count() * 1e-9);
        cout << "      PressureDensity ended" << endl;
    }

private:
    cAtmosphereModel& m;
};
