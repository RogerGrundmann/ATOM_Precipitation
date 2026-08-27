/*
 * Dry convective adjustment — Manabe & Strickler (1964).
 *
 * Ported from the family's shared planet/ConvectiveAdjustment.h (ATURAN/ATJUP/ATSAT), whose
 * header comment is the reference for why a model needs this at all and is worth reading. The
 * two properties that matter: it removes the instability COMPLETELY rather than damping it, and
 * it does not create or destroy energy.
 *
 * WHY THIS TREE NEEDS IT, MEASURED 2026-08-27. ATOM_Precipitation had NO dry convective
 * adjustment at all -- no such file, no call site, only MoistConvection.h -- so nothing in the
 * model could remove a superadiabatic layer. Measured over 100 iterations at 24 threads: at
 * iteration 0 the surface layer is STABLE (-3.93 K/km) and brunt_N2 < 0 in 0.000 of columns; by
 * iteration 20 the lowest 39 m runs at -18.05 K/km, TWICE the dry adiabat, with brunt_N2 < 0 in
 * 59 % of columns at i = 0 and i = 1 and 0 % at i = 2 -- one unstable layer, exactly the
 * surface-to-first-level step. It is manufactured by the time loop: the surface warms +0.76 K
 * over the run while the air above warms +0.15 K, and nothing mixes them.
 *
 * NOTE WHAT IT IS NOT. MultiLayerRadiation is called seven times, ALL during setup (lines
 * 713-889) while the iteration loop starts at line 985, so MLR cannot be maintaining this and a
 * surface-balance repair does not touch it (ATM_SFC_COUPLED was written on that wrong premise,
 * measured a null, and removed). The separate finding that MLR's surface balance debits
 * c_H*(T_s - T_air1) which `Q_Sensible` never credits to the air -- it is written in
 * RHS_Atm_Turb.cpp:485 and read by NOTHING -- is real, but it is an INITIALISATION defect and
 * not the cause of the standing instability.
 *
 * WHAT WAS ADAPTED FROM ATHAD'S COPY. Two of its assumptions do not hold here:
 *
 *   1. ATHAD HAS NO TOPOGRAPHY (its invariant 1), so its fluid column is the whole column and it
 *      hard-codes i0 = 0 -- with a comment saying it was kept as i0/i1 precisely so a sibling
 *      with topography could pick it up. Here i0 = i_topography[j][k], so the adjustment starts
 *      at the ground and never mixes rock into air.
 *
 *   2. ATHAD TAKES cp LOCALLY from AtmMixture::cp_of(), because its cp varies ~2x across
 *      300-1500 K and follows the composition. This tree has no MixtureAtm and its cp variation
 *      is small, so cp_l is used and saying so is better than importing machinery that is not
 *      here. That is the same judgement brunt_N2's port made about kappa.
 *
 * WHAT WAS KEPT, and it is the reason this is a port rather than a copy of the family's shared
 * planet/ConvectiveAdjustment.h: the per-layer critical drop from get_layer_height() and the
 * cumulative-sum segment adiabat. The shared file computes dz_m = L_atm*1e3/(im-1), a layer
 * thickness only on a UNIFORM grid; this tree's radial coordinate is exponentially stretched, so
 * a single critical drop would be far too strict at the bottom and too lax at the top.
 *
 * WHAT IT DELIBERATELY DOES NOT DO, carried over unchanged: it mixes temperature only, not
 * composition, and it uses the DRY adiabat even where cloud is present. A moist adjustment would
 * use the saturated lapse — ATHAD_COND has SaturationH2O::moistLapse for exactly that — but that
 * is a modelling decision with consequences for the microphysics already running in the
 * saturation adjustment, and it is not made here.
 */

#pragma once

#include "cAtmosphereModel.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

class ConvectiveAdjustment {
public:
    explicit ConvectiveAdjustment(cAtmosphereModel& model) : m(model) {}

    // ATM_CONV_ADJ -- default OFF. This tree has never had a dry adjustment, so switching one on
    // changes the thermal structure of every unstable column; it is measured before any flip.
    static bool enabled(){
        static const bool v = [](){
            const char* e = getenv("ATM_CONV_ADJ"); return e && atoi(e) != 0; }();
        return v;
    }

    void run()
    {
        using namespace std;
        cout << endl << "      ATOM: ConvectiveAdjustment (dry, Manabe-Strickler)" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        // ATM_CONV_ADJ_LAPSE scales the critical lapse rate: 0 gives an isothermal criterion,
        // values above 1 make the scheme stricter than the dry adiabat, which is one crude way
        // to stand in for a moist adiabat in a condensing region.
        static const double lapse_fac = [](){ const char* e = getenv("ATM_CONV_ADJ_LAPSE");
                                              return e ? atof(e) : 1.0; }();
        // Sweeps repeat until the column is stable; the cap only exists so a pathological
        // column cannot spin here.
        static const int max_pass = [](){ const char* e = getenv("ATM_CONV_ADJ_PASSES");
                                          const int v = e ? atoi(e) : 64;
                                          return v > 0 ? v : 64; }();

        // A column is left alone unless it is superadiabatic by more than this, so round-off
        // does not make the scheme fire on a column that is already neutral.
        constexpr double tol_nd = 1.0e-12;

        // Layer thicknesses of the stretched grid, once.
        std::vector<double> dz(m.im, 0.0);
        for (int i = 0; i + 1 < m.im; i++)
            dz[i] = m.get_layer_height(i+1) - m.get_layer_height(i);

        long long n_columns_adjusted = 0, n_layers_mixed = 0;
        int    max_passes_used = 0;
        double max_dT_K = 0.0, max_rel_drift = 0.0;

        std::vector<double> w(m.im), t_col(m.im), dT_ad(m.im);

        #pragma omp parallel for collapse(2) schedule(static) firstprivate(w, t_col, dT_ad) \
            reduction(+:n_columns_adjusted, n_layers_mixed) \
            reduction(max:max_passes_used, max_dT_K, max_rel_drift)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {

                // The fluid column starts at the ground, not at level 0 -- this tree has
                // topography. Below i_topography the cells are rock and carry a barometric
                // p_stat that would otherwise be mixed in as if it were air.
                const int i0 = m.i_topography[j][k], i1 = m.im - 1;
                if (i1 - i0 < 1) continue;

                // Mass weights and the per-layer critical drop.
                //
                // The enthalpy per unit area of a layer is (cp/g)*dp, so the conserved quantity
                // is the dp-weighted temperature. dp comes from the hydrostatic p_stat with
                // faces at the midpoints; the end layers get their half-cell. Using dp rather
                // than a density avoids a trap: rho = p/(R*T) makes rho*T identically p/R, so a
                // density weight would conserve nothing at all.
                bool ok = true;
                for (int i = i0; i <= i1; i++) {
                    const double p_lo = (i > i0) ? 0.5 * (m.p_stat.x[i-1][j][k] + m.p_stat.x[i][j][k])
                                                 : m.p_stat.x[i0][j][k];
                    const double p_hi = (i < i1) ? 0.5 * (m.p_stat.x[i][j][k] + m.p_stat.x[i+1][j][k])
                                                 : m.p_stat.x[i1][j][k];
                    w[i]     = p_lo - p_hi;
                    t_col[i] = m.t.x[i][j][k];
                    if (!(w[i] > 0.0) || !std::isfinite(w[i]))     ok = false;
                    if (!std::isfinite(t_col[i]) || t_col[i] <= 0.0) ok = false;

                    // Critical drop across layer i -> i+1, non-dimensional, on the LOCAL cp.
                    if (i < i1) {
                        const double cp_loc = (m.cp_l > 0.0) ? m.cp_l : 1005.0;   // constant here; see header
                        dT_ad[i] = lapse_fac * (m.g / cp_loc) * dz[i] / m.t_0;
                    }
                }
                // A column with a non-monotonic p_stat or a bad temperature is left untouched:
                // this routine is not the place to repair either, and mixing across a NaN would
                // spread it through the whole column.
                if (!ok) continue;

                double sum_before = 0.0;
                for (int i = i0; i <= i1; i++) sum_before += w[i] * t_col[i];

                // Whole unstable SEGMENTS are mixed at once, not adjacent pairs: pairwise mixing
                // moves heat one layer per sweep, so a deep unstable block needs as many sweeps
                // as it has layers. A segment [a..b] is put on the adiabat T_q = C - D_q, where
                // D_q is the CUMULATIVE adiabatic drop from a to q (the shared file's linear
                // dT_ad*(q-a), generalised to a stretched grid), with C fixed by conserving the
                // dp-weighted temperature:  C = [sum w_q T_q + sum w_q D_q] / sum w_q.
                int passes = 0;
                long long layers_mixed = 0;
                bool changed = true;
                std::vector<double> D;
                while (changed && passes < max_pass) {
                    changed = false;
                    passes++;
                    int i = i0;
                    while (i < i1) {
                        if (t_col[i] - t_col[i+1] <= dT_ad[i] + tol_nd) { i++; continue; }

                        int a = i, b = i + 1;
                        for (;;) {
                            D.assign(b - a + 1, 0.0);
                            for (int q = a + 1; q <= b; q++)
                                D[q - a] = D[q - a - 1] + dT_ad[q - 1];

                            double sw = 0.0, swt = 0.0, swd = 0.0;
                            for (int q = a; q <= b; q++) {
                                sw  += w[q];
                                swt += w[q] * t_col[q];
                                swd += w[q] * D[q - a];
                            }
                            const double C = (swt + swd) / sw;
                            for (int q = a; q <= b; q++) {
                                const double t_new = C - D[q - a];
                                const double dK = std::fabs(t_new - t_col[q]) * m.t_0;
                                if (dK > max_dT_K) max_dT_K = dK;
                                t_col[q] = t_new;
                            }

                            // Mixing can destabilise the joint with the layer below or above, so
                            // the segment grows in whichever direction is still too steep and is
                            // re-mixed. It can only grow, and only within the column, so this
                            // terminates.
                            bool extended = false;
                            if (a > i0 && t_col[a-1] - t_col[a] > dT_ad[a-1] + tol_nd) { a--; extended = true; }
                            if (b < i1 && t_col[b] - t_col[b+1] > dT_ad[b]   + tol_nd) { b++; extended = true; }
                            if (!extended) break;
                        }

                        layers_mixed += (b - a + 1);
                        changed = true;
                        i = b;                          // carry on above the block just mixed
                    }
                }

                if (layers_mixed == 0) continue;

                double sum_after = 0.0;
                for (int i = i0; i <= i1; i++) sum_after += w[i] * t_col[i];
                const double drift = (sum_before != 0.0)
                                   ? std::fabs(sum_after - sum_before) / std::fabs(sum_before) : 0.0;
                if (drift > max_rel_drift) max_rel_drift = drift;

                for (int i = i0; i <= i1; i++) m.t.x[i][j][k] = t_col[i];

                n_columns_adjusted++;
                n_layers_mixed += layers_mixed;
                if (passes > max_passes_used) max_passes_used = passes;
            }
        }

        const double frac = 100.0 * double(n_columns_adjusted) / double(m.jm * m.km);
        printf("      ATOM: convective adjustment — %lld of %d columns (%.2f %%), %lld layers,"
               " worst column %d sweeps of %d, max dT %.3f K, enthalpy drift %.2e\n",
               n_columns_adjusted, m.jm * m.km, frac, n_layers_mixed,
               max_passes_used, max_pass, max_dT_K, max_rel_drift);
        if (max_passes_used >= max_pass)
            cout << "      ATOM: WARNING - the sweep cap was reached, a column may still be "
                    "superadiabatic (raise ATM_CONV_ADJ_PASSES)" << endl;

        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for ConvectiveAdjustment\n", elapsed.count() * 1e-9);
        cout << "      ATOM: ConvectiveAdjustment ended" << endl;
    }

private:
    cAtmosphereModel& m;
};
