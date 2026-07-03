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

// ============================================================================
// MultiLayerRadiation — friend class of cAtmosphereModel
//
// Multi-layer long/short-wave radiation model for the surface-temperature
// distribution (computation of the local temperature from short- and long-wave
// radiation). For each (j,k) column it builds:
//   - the grey-body emission of every layer (Stefan-Boltzmann, sigma T^4),
//   - the layer emissivity epsilon from the water-vapour law (Bignami 1995,
//     valid -40..45 C) plus an optional CO2 term (Atwater & Ball) — the CO2
//     contribution is presently disabled (eps_co2 = 0; the CO2 coupling is
//     rough and to be updated),
// then solves the layer radiative-balance as a tridiagonal system with the
// Thomas algorithm and inverts the resulting radiation back to temperature.
// The surface layer receives the albedo-reduced incoming short-wave flux.
//
// Refactored from the former free member function
// cAtmosphereModel::RadiationMultiLayer() into the project's friend-class idiom
// (mirrors ThermoAtm / TurbulenceAtm / PressureSolverAtm). The (j,k) columns are
// independent, so the outer j loop is OpenMP-parallel with ALL Thomas-solver
// scratch (step/alfa/beta/AA/CA/CC and the per-column radiation_original)
// thread-local — the original reused a single shared set of scratch vectors,
// which would race under OpenMP.
//
// Reads : t, c, co2, p_stat, h, layer heights, albedo/short-wave/emissivity
//         constants (albedo_pole/equator, rad_pole/equator_short, sigma, ep,
//         co2_0, t_0, p_0).
// Writes: t (updated), radiation, epsilon, epsilon_2D, albedo, short_wave_radiation.
//
// Usage:  MultiLayerRadiation(*this).run();
// ============================================================================
class MultiLayerRadiation {
public:
    explicit MultiLayerRadiation(cAtmosphereModel& model)
        : m(model)
    {}

    void run()
    {
        using namespace std;
        cout << endl << "      RadiationMultiLayer" << endl;
        auto begin = std::chrono::high_resolution_clock::now();

        const int j_max  = m.jm - 1;
        const int j_half = j_max / 2;

        // ---- latitude profiles (computed once; read-only in the parallel loop) ----

        // Effective surface albedo: pole -> equator parabola over any surface cell.
        const double albedo_eff = m.albedo_pole - m.albedo_equator;
        #pragma omp parallel for schedule(static)
        for (int j = 0; j < m.jm; j++)
            for (int k = 0; k < m.km; k++)
                for (int i = 0; i < m.im - 1; i++)
                    if (is_ocean_surface(m.h, i, j, k) || is_land_surface(m.h, i, j, k)) {
                        m.albedo.y[j][k] =
                            albedo_eff * parabola((double)j / (double)j_half) + m.albedo_pole;
                        break;   // albedo depends only on latitude — first surface cell suffices
                    }

        // Incoming short-wave radiation: pole -> equator parabola, hemispherically symmetric.
        m.short_wave_radiation = std::vector<double>(m.jm, m.rad_pole_short);
        const double rad_short_eff = m.rad_pole_short - m.rad_equator_short;
        for (int j = j_half; j >= 0; j--)
            m.short_wave_radiation[j] =
                rad_short_eff * parabola((double)j / (double)j_half) + m.rad_pole_short;
        for (int j = j_max; j > j_half; j--)
            m.short_wave_radiation[j] = m.short_wave_radiation[j_max - j];

        // ---- per-column radiative balance (columns independent -> OpenMP over j) ----
        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < m.jm; j++) {

            // Thread-local Thomas-solver scratch; reused across k within this j.
            // Every column fully overwrites the entries it later reads (i_mount = 0,
            // i_trop = im-1), so reuse is race-free.
            std::vector<double> dp_col(m.im, 0.0), vpath_col(m.im, 0.0);
            std::vector<double> alfa(m.im, 0.0), beta(m.im, 0.0);
            std::vector<double> AA(m.im, 0.0), CA(m.im, 0.0);
            std::vector<double> radiation_original(m.im, 0.0);
            std::vector<std::vector<double> > CC(m.im, std::vector<double>(m.im, 0.0));

            const int i_trop  = m.im - 1;   // top layer (tropopause proxy)
            const int i_mount = 0;          // surface / bottom layer

            for (int k = 0; k < m.km; k++) {

                // Grey-body emission of each layer and its "original" reference.
                for (int i = i_mount; i <= i_trop; i++) {
                    m.radiation.x[i][j][k] = m.sigma * pow(m.t.x[i][j][k] * m.t_0, 4.0);
                    radiation_original[i]  = m.radiation.x[i][j][k];
                }

                // Layer emissivity — DE-SATURATED per-layer grey optical depth.
                //
                // The original code applied the Bignami (1995) clear-sky COLUMN emissivity
                // eps = 0.684 + 0.0056*e per layer, so every one of the ~40 layers was
                // ~0.7-opaque -> the troposphere was a near-perfect blackbody, OLR collapsed
                // to ~125 W/m2 (Earth ~240) and the emission level pinned to the cold model
                // top (see project_multilayer_radiation). Bignami is a whole-column quantity,
                // not a per-layer one.
                //
                // Fix: preserve the Bignami column value but SPLIT it and distribute each
                // part physically as an optical depth, so thin layers are optically thin:
                //   - dry baseline (0.684, well-mixed CO2/continuum) -> by layer mass  dp
                //   - water-vapour part (0.0056*e_surf)              -> by vapour path c*dp
                //   eps_i = 1 - exp(-tau_i),  tau_i = tau_dry*dp_i/Sum(dp) + tau_wv*vp_i/Sum(vp)
                // Gives OLR ~263 W/m2 with the water-vapour greenhouse feedback retained.
                // The CO2 anomaly forcing is handled separately (5.35*ln(C/C0)); no CO2
                // emissivity term here.
                const double eps_dry = 0.684;                     // Bignami dry-air baseline
                double sum_dp = 0.0, sum_vp = 0.0;
                for (int i = i_mount; i <= i_trop; i++) {
                    double dp = (i < i_trop) ? (m.p_stat.x[i][j][k] - m.p_stat.x[i+1][j][k])
                                             : m.p_stat.x[i][j][k];   // top layer: all mass above
                    if (dp < 0.0) dp = 0.0;
                    double cw    = (m.c.x[i][j][k] > 0.0) ? m.c.x[i][j][k] : 0.0;
                    dp_col[i]    = dp;
                    vpath_col[i] = cw * dp;                          // ~ layer precipitable water
                    sum_dp      += dp;
                    sum_vp      += vpath_col[i];
                }
                const double e_surf  = m.c.x[i_mount][j][k] * m.p_stat.x[i_mount][j][k] / m.ep; // [hPa]
                double eps_col       = eps_dry + 0.0056 * e_surf;    // Bignami column emissivity
                if (eps_col > 0.999) eps_col = 0.999;
                const double tau_dry = -log(1.0 - eps_dry);          // dry-baseline column optical depth
                const double tau_col = -log(1.0 - eps_col);
                const double tau_wv  = (tau_col > tau_dry) ? (tau_col - tau_dry) : 0.0;
                const double inv_dp  = (sum_dp > 0.0) ? 1.0 / sum_dp : 0.0;
                const double inv_vp  = (sum_vp > 0.0) ? 1.0 / sum_vp : 0.0;
                for (int i = i_mount; i <= i_trop; i++) {
                    // CO2 band contribution (the former MLR CO2 integration, restored and
                    // un-zeroed). Well-mixed CO2 partial pressure P_c -> layer absorber path
                    // u_c [atm*cm] -> Atwater & Ball band emissivity (curve-of-growth 0.185,
                    // 0.3919; Byun & Chen 2013). Expressed as an OPTICAL DEPTH so it composes
                    // additively with the de-saturated dry/water-vapour optical depths:
                    // eps_i = 1 - exp(-(tau_dry_i + tau_wv_i + tau_co2_i)).
                    //
                    // co2_band_scale replaces the original Atwater 0.5 fit factor: 0.5 trapped
                    // ~65 W/m2 (too strong vs real CO2), so it is TUNED to 0.17, giving ~30 W/m2
                    // of CO2 greenhouse on a US-standard column (OLR ~263 -> ~233), matching
                    // Earth's clear-sky CO2 contribution. This grey band saturates against the
                    // 16 km model top so it does NOT give 3.7 W/m2/doubling — the calibrated
                    // per-doubling forcing is the separate 5.35*ln(C/C0) t_eq shift
                    // (cAtmosphereModel.cpp). If MLR is ever wired AND that t_eq forcing is on,
                    // CO2 acts twice; reconcile then (see project_multilayer_radiation).
                    const double co2_band_scale = 0.17;
                    const double dz  = (i < i_trop) ? (m.get_layer_height(i+1) - m.get_layer_height(i))
                                                    : (m.get_layer_height(i) - m.get_layer_height(i-1));
                    const double P_c = 1e-6 * m.p_stat.x[i][j][k] * m.co2.x[i][j][k] * m.co2_0 / m.p_0; // [atm]
                    double u_c = P_c * dz * 100.0;                                                       // [atm*cm]
                    if (u_c < 0.0) u_c = 0.0;
                    double eps_co2 = co2_band_scale * 0.185 * (1.0 - exp(-0.3919 * pow(u_c, 0.4)));       // Atwater & Ball
                    if (eps_co2 > 0.999) eps_co2 = 0.999;
                    const double tau_co2 = -log(1.0 - eps_co2);

                    double tau = tau_dry * dp_col[i] * inv_dp + tau_wv * vpath_col[i] * inv_vp + tau_co2;
                    m.epsilon.x[i][j][k] = 1.0 - exp(-tau);
                }
                m.epsilon_2D.y[j][k] = m.epsilon.x[i_mount][j][k];

                // Transmitted (AA) / absorbed (CC diagonal) radiation, and the sum CA
                // of all radiations transmitted through each layer.
                AA[i_mount]          = m.radiation.x[i_mount][j][k];              // surface radiation
                CC[i_mount][i_mount] = m.epsilon.x[i_mount][j][k] * m.radiation.x[i_mount][j][k];
                for (int i = i_mount + 1; i <= i_trop; i++) {
                    AA[i]    = AA[i - 1] * (1.0 - m.epsilon.x[i][j][k]);          // transmitted from each layer
                    CC[i][i] = m.epsilon.x[i][j][k] * m.radiation.x[i][j][k];     // absorbed in each layer
                }
                for (int i = i_mount + 2; i <= i_trop; i++) {
                    CA[i] = 0.0;
                    for (int l = 1; l <= i - 1; l++) {
                        CC[l][i] = CC[l][i - 1] * (1.0 - m.epsilon.x[i][j][k]);   // transmitted past layer i
                        CA[i] += CC[l][i];                                        // sum over all l
                    }
                }

                // Thomas algorithm — forward elimination (alfa/beta recurrence).
                // The sweep now INCLUDES the top row i = i_trop. The old code stopped at
                // i_trop-1 and closed the system with an ad-hoc top formula that divided by
                // (CA[i_trop]-CA[i_trop-1]) — a difference of two near-equal cumulative sums
                // that collapses to ~0 and flips sign in the quasi-isothermal upper
                // atmosphere, seeding a huge negative top radiation that back-substitution
                // smeared into NaN temperatures aloft (see project_multilayer_radiation).
                double aa, bb, cc, dd;
                for (int i = i_mount; i <= i_trop; i++) {
                    if (i == i_mount) {
                        aa = 0.0;
                        bb = -2.0 * m.radiation.x[i][j][k];
                        cc = m.epsilon.x[i + 1][j][k] * m.radiation.x[i + 1][j][k];
                        dd = -(1.0 - m.albedo.y[j][k]) * m.short_wave_radiation[j];
                        alfa[i] = -cc / bb;
                        beta[i] = +dd / bb;
                    }
                    if (i == i_mount + 1) {
                        aa = m.radiation.x[i - 1][j][k];
                        bb = -2.0 * m.epsilon.x[i][j][k] * m.radiation.x[i][j][k];
                        cc = m.epsilon.x[i + 1][j][k] * m.radiation.x[i + 1][j][k];
                        dd = AA[i];
                        alfa[i] = -cc / (bb + aa * alfa[i - 1]);
                        beta[i] = +(dd - aa * beta[i - 1]) / (bb + aa * alfa[i - 1]);
                    }
                    if (i == i_mount + 2) {
                        aa = m.epsilon.x[i - 1][j][k] * m.radiation.x[i - 1][j][k];
                        bb = -2.0 * m.epsilon.x[i][j][k] * m.radiation.x[i][j][k];
                        cc = m.epsilon.x[i + 1][j][k] * m.radiation.x[i + 1][j][k];
                        dd = -AA[i - 1] + AA[i] + CC[i - 1][i];
                        alfa[i] = -cc / (bb + aa * alfa[i - 1]);
                        beta[i] = +(dd - aa * beta[i - 1]) / (bb + aa * alfa[i - 1]);
                    }
                    if (i > i_mount + 2) {
                        aa = m.epsilon.x[i - 1][j][k] * m.radiation.x[i - 1][j][k];
                        bb = -2.0 * m.epsilon.x[i][j][k] * m.radiation.x[i][j][k];
                        // Top row (i == i_trop) has no layer above -> cc = 0. Guarding this
                        // also avoids the out-of-bounds read of radiation[i_trop+1] that the
                        // extended sweep would otherwise make.
                        cc = (i < i_trop) ? m.epsilon.x[i + 1][j][k] * m.radiation.x[i + 1][j][k]
                                          : 0.0;
                        dd = -AA[i - 1] + AA[i] - CA[i - 1] + CA[i];
                        alfa[i] = -cc / (bb + aa * alfa[i - 1]);
                        beta[i] = +(dd - aa * beta[i - 1]) / (bb + aa * alfa[i - 1]);   // FIX: aa (sub-diagonal), was alfa[i]
                    }
                }

                // Back-substitution (Thomas). The top unknown is beta[i_trop] because the
                // top row's alfa[i_trop] = 0 (cc = 0 — no layer above).
                m.radiation.x[i_trop][j][k] = beta[i_trop];
                for (int i = i_trop - 1; i >= 0; i--)
                    m.radiation.x[i][j][k] = alfa[i] * m.radiation.x[i + 1][j][k] + beta[i];

                // Radiation -> temperature (add back the reference emission, invert sigma T^4).
                for (int i = 0; i <= i_trop; i++) {
                    m.radiation.x[i][j][k] = radiation_original[i] + m.radiation.x[i][j][k];
                    m.t.x[i][j][k] = pow(m.radiation.x[i][j][k] / m.sigma, 0.25) / m.t_0;
                }

                // ---- Surface energy balance (radiative-CONVECTIVE) ----
                // The tridiagonal inversion left the surface T ~insensitive to the longwave
                // opacity, so CO2 did not warm the surface. Replace the surface value with an
                // explicit balance:  (1-albedo)*SW + L_down = sigma*T_s^4 + c_H*(T_s - T_air1).
                //  - L_down = downwelling longwave (back-radiation) = sum of atmospheric-layer
                //    emissions transmitted down to the surface; it RISES with CO2/H2O emissivity,
                //    so more greenhouse -> warmer surface (the physically-correct response).
                //  - c_H*(T_s - T_air1): bulk turbulent (sensible+latent) flux to the lowest air
                //    layer, which keeps the surface off the pure-radiative overheating —
                //    radiative-convective, not pure radiative. sigma*T_s^4 is linearised about
                //    the current surface T_s0 (one Newton step). See project_multilayer_radiation.
                double L_down = 0.0, trans = 1.0;
                for (int i = i_mount + 1; i <= i_trop; i++) {
                    L_down += m.epsilon.x[i][j][k] * m.sigma
                            * pow(m.t.x[i][j][k] * m.t_0, 4.0) * trans;     // layer i emission reaching surface
                    trans  *= (1.0 - m.epsilon.x[i][j][k]);                // attenuation through layer i
                }
                const double SW_abs = (1.0 - m.albedo.y[j][k]) * m.short_wave_radiation[j];
                const double T_air1 = m.t.x[i_mount + 1][j][k] * m.t_0;     // lowest air-layer T [K]
                const double T_s0   = m.t.x[i_mount][j][k] * m.t_0;         // linearisation point [K]
                const double c_H    = 15.0;                                // bulk turbulent transfer [W/m2/K]
                const double dsigT4 = 4.0 * m.sigma * T_s0 * T_s0 * T_s0;
                const double T_s    = (SW_abs + L_down - m.sigma * pow(T_s0, 4.0)
                                       + dsigT4 * T_s0 + c_H * T_air1) / (dsigT4 + c_H);
                m.radiation.x[i_mount][j][k] = m.sigma * pow(T_s, 4.0);
                m.t.x[i_mount][j][k]         = T_s / m.t_0;
            }  // k
        }  // j

        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for MultiLayerRadiation\n", elapsed.count() * 1e-9);
        cout << "      RadiationMultiLayer ended" << endl;
    }

private:
    cAtmosphereModel& m;
};
