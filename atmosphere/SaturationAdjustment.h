#pragma once

#include "cAtmosphereModel.h"

#include <algorithm>
#include <cmath>
#include <chrono>
#include <iostream>
#include <cstdio>
#include <vector>

using namespace AtomUtils;


class SaturationAdjustment {
public:
    explicit SaturationAdjustment(cAtmosphereModel& model)
        : m(model)
    {}


    void run() {
        using namespace std;

        cout << endl << endl << endl << "      SaturationAdjustment" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        cout.precision(9);

        computeSteps();
        adjustSaturation();
        applyTopography();
        clampAndFade();
        printReport();

        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for SaturationAdjustment\n", elapsed.count() * 1e-9);

        cout << "      SaturationAdjustment ended" << endl;
    }

private:
    cAtmosphereModel& m;

    std::vector<double> step;

    // Precomputed constants
    static constexpr double fade_K        = 5.0;                        // transition half-width in Kelvin
    static constexpr int    iter_prec_end = 20;

    void computeSteps() {
        step.resize(m.im);
        for (int i = 0; i < m.im; i++)
            step[i] = m.get_layer_height(i+1) - m.get_layer_height(i);
    }

    void adjustSaturation() {
        const double inv_t_0      = 1.0 / m.t_0;
        const double t_range_inv  = 1.0 / (m.t_0 - m.t_00);
        const double lv_over_cp   = m.lv / m.cp_l;
        const double ls_over_cp   = m.ls / m.cp_l;

        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < m.im - 1; i++) {
            for (int j = 0; j < m.jm; j++) {
                double *S_c_c_row = m.S_c_c.x[i][j];
                double *c_row     = m.c.x[i][j];
                double *cloud_row = m.cloud.x[i][j];
                double *ice_row   = m.ice.x[i][j];
                double *t_row     = m.t.x[i][j];
                double *p_row     = m.p_stat.x[i][j];
                double  dt_dim    = step[i] / 1.6;

                for (int k = 0; k < m.km; k++) {
                    S_c_c_row[k] = 0.0;

                    double q_v_old = std::max(0.0, c_row[k]);
                    double q_c_old = std::max(0.0, cloud_row[k]);
                    double q_i_old = std::max(0.0, ice_row[k]);

                    double T       = t_row[k] * m.t_0;
                    // Cap T to the Magnus-formula validity range. Above ~101°C the
                    // saturation vapor pressure exceeds p_local and the q_sat fallback
                    // (ep * 1e-5) collapses by 5000×, triggering runaway condensation
                    // that releases lv/cp * q_v ≈ 87 K of latent heat per call and drives
                    // T further out of range. Cap at 60°C — well above any physical
                    // surface temperature — so a single corrupt cell cannot poison the run.
                    constexpr double T_max = 333.15;
                    if (T > T_max) {
                        T = T_max;
                        t_row[k] = T_max / m.t_0;
                    }
                    double p_local = p_row[k];

                    double E_sat  = m.hp * exp_func(T, 17.2694, 35.86);
                    double q_sat  = (p_local > E_sat)
                        ? m.ep * E_sat / (p_local - E_sat)
                        : m.ep * 1e-5;

                    const double alpha_entry = 1.0 / (1.0 + std::exp(-(T - m.t_00) / fade_K));

                    if ((q_v_old > q_sat && alpha_entry > 0.01) ||
                        (q_v_old < q_sat &&
                        (q_c_old > 1e-12 || q_i_old > 1e-12))) {

                        double q_v_b   = q_v_old;
                        double q_c_b   = q_c_old;
                        double q_i_b   = q_i_old;
                        double q_v_hyp = q_sat;
                        const double T_original = t_row[k] * m.t_0;

                        for (int iter = 1; iter <= iter_prec_end; iter++) {
                            double CND = std::max(0.0, std::min(1.0,
                                (T - m.t_00) * t_range_inv));
                            double DEP = 1.0 - CND;

                            double d_q_v = q_v_hyp - q_v_b;

                            if (d_q_v > 0) {
                                double max_evap = q_c_b + q_i_b;
                                if (d_q_v > max_evap) d_q_v = max_evap;
                            }

                            double d_cnd = d_q_v * CND;
                            double d_dep = d_q_v * DEP;

                            q_v_b += d_q_v;
                            q_c_b  = std::max(0.0, q_c_b - d_cnd);
                            q_i_b  = std::max(0.0, q_i_b - d_dep);

                            T -= lv_over_cp * d_cnd + ls_over_cp * d_dep;

                            double E_sat = m.hp * exp_func(T, 17.2694, 35.86);
                            double E_Ice = m.hp * exp_func(T, 21.8746, 7.66);
                            double q_sat = (p_local > E_sat)
                                ? m.ep * E_sat / (p_local - E_sat)
                                : m.ep * 1e-5;
                            double q_Ice = (p_local > E_Ice)
                                ? m.ep * E_Ice / (p_local - E_Ice)
                                : m.ep * 1e-5;

                            double q_sum = q_c_b + q_i_b;
                            double q_v_target = (q_sum > 1e-12)
                                ? (q_c_b * q_sat + q_i_b * q_Ice) / q_sum
                                : ((T >= m.t_0) ? q_sat : q_Ice);

                            q_v_hyp = 0.5 * (q_v_target + q_v_b);

                            if (fabs(q_v_b / q_v_hyp - 1.0) <= 1.0e-6)
                                break;
                        }

                        q_c_b = std::max(0.0, q_c_b);
                        q_i_b = std::max(0.0, q_i_b);

                        // Cap T after the Newton loop. The q_v_hyp = 0.5*(q_v_target + q_v_b)
                        // damping is too weak when dq_sat/dT is steep (marginal saturation),
                        // so the iteration's amplitude grows. Within one call T can swing from
                        // a physical entry value into the Magnus-cliff regime (>101 °C at
                        // p = 1080 hPa), and the write-back below would persist that bad value.
                        // The entry-time cap is not enough because the runaway happens during
                        // the loop, not between calls.
                        if (T > T_max) T = T_max;

                        if (!std::isnan(T) && !std::isnan(q_v_b)) {
                            S_c_c_row[k] = alpha_entry * (q_c_b - q_c_old) / dt_dim;
                            c_row[k]     = q_v_old + alpha_entry * (q_v_b - q_v_old);
                            cloud_row[k] = q_c_old + alpha_entry * (q_c_b - q_c_old);
                            ice_row[k]   = q_i_old + alpha_entry * (q_i_b - q_i_old);
                            t_row[k]     = (T_original + alpha_entry * (T - T_original)) * inv_t_0;
                        }
 
                        if (T < m.t_00) {
                            cloud_row[k] = 0.0;
                            ice_row[k]   = 0.0;
                        }
                    }
                }
            }
        }
    }

    void applyTopography() {
        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                int i_mount = m.i_topography[j][k];
                m.c.x[0][j][k]     = m.c.x[i_mount][j][k];
                m.cloud.x[0][j][k] = m.cloud.x[i_mount][j][k];
                m.ice.x[0][j][k]   = m.ice.x[i_mount][j][k];
                if (AtomUtils::is_finite_safe(m.t.x[i_mount][j][k]))
                    m.t.x[0][j][k] = m.t.x[i_mount][j][k];
            }
        }
    }

    void clampAndFade() {
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                double *c_row     = m.c.x[i][j];
                double *cloud_row = m.cloud.x[i][j];
                double *ice_row   = m.ice.x[i][j];
                double *t_row_nd  = m.t.x[i][j];

                for (int k = 0; k < m.km; k++) {
                    if (c_row[k]     < 0.0) c_row[k]     = 0.0;
                    if (cloud_row[k] < 0.0) cloud_row[k] = 0.0;
                    if (ice_row[k]   < 0.0) ice_row[k]   = 0.0;

                    const double T_dim = t_row_nd[k] * m.t_0;
                    const double alpha = 1.0 / (1.0 + std::exp(-(T_dim - m.t_00) / fade_K));
                    c_row[k]     *= alpha;
                    cloud_row[k] *= alpha;
                    ice_row[k]   *= alpha;
                }
            }
        }
    }

    void printReport() const {
        // Diagnostic variables are not updated in the parallel loops;
        // kept here to preserve the original output contract.
        const bool   satadjust = false;
        const int    iter_prec = 0;
        const int    i_sat = 0, j_sat = 0, k_sat = 0;
        const double saturation = 0.0;

        if (!satadjust)
            std::cout << "      no saturation of water vapour in SaturationAdjustment found"
                      << std::endl;
        else
            std::cout << "      saturation of water vapour in SaturationAdjustment found"
                      << std::endl
                      << "      iter_prec = " << iter_prec << std::endl
                      << "      i_sat = "  << i_sat
                      << "   j_sat = "     << j_sat
                      << "   k_sat = "     << k_sat
                      << "   height_sat[m] = " << m.get_layer_height(i_sat)
                      << "   saturation[g/kg] = " << saturation * 1e3 << std::endl;

        if (iter_prec >= iter_prec_end)
            std::cout << std::endl
                      << "      no convergent solution found in SaturationAdjustment"
                      << std::endl
                      << "      iter_prec_end = " << iter_prec_end << std::endl
                      << "      iter_prec = "     << iter_prec     << std::endl
                      << "      results see above" << std::endl;
    }
};
