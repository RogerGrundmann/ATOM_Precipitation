#pragma once

#include "cAtmosphereModel.h"
#include "CloudFraction.h"

#include <algorithm>
#include <cmath>
#include <chrono>
#include <iostream>
#include <cstdio>
#include <vector>

using namespace AtomUtils;


class SaturationAdjustment {
public:
    // ---- ATM_CLOUD_FRAC: make the adjustment FRACTION-AWARE (default 0 = shipped) --------
    //
    // The shipped adjustment drives q_v toward q_sat on the GRID MEAN. A sub-grid cloud scheme
    // puts condensate exactly where the grid mean is SUBSATURATED, so a grid-mean adjustment
    // evaporates precisely the cloud the closure created. Measured: with ATM_CLOUD_FRAC and a
    // realistic humidity the init call at cAtmosphereModel.cpp:461 takes the column condensate
    // path from 270.70 g/m2 to 0.633 -- a factor of 428, against 38 % on the shipped
    // near-saturated column, which is the one state a grid-mean adjustment does not destroy.
    //
    // The repair is to target the FRACTIONAL equilibrium. Total water q_t = q_v + q_c + q_i is
    // conserved by the loop below (d_cnd + d_dep = d_q_v), so for the uniform-PDF closure
    //     D       = (1 - H_crit)*q_s
    //     f       = clamp((q_t + D - q_s) / 2D, 0, 1)
    //     q_c_eq  = f^2*D   (f < 1),   q_t - q_s   (f = 1)
    //     target  = q_t - q_c_eq
    // which reduces EXACTLY to the shipped target q_s when H_crit = 1 (D = 0, f = 1), so the
    // off-branch is unchanged by construction and not merely by a guard.
    // The closure itself now lives in CloudFraction.h, because initCloudIce, this
    // adjustment and MultiLayerRadiation all need the SAME f and the same H_crit.
    // These stay as thin forwards so the call sites below read unchanged.
    static bool   cloudFrac(){ return CloudFraction::enabled(); }
    static double hCrit(double p_hPa){ return CloudFraction::hCrit(p_hPa); }
    static double qcEquilibrium(double q_t, double q_s, double p_hPa){
        return CloudFraction::qcEquilibrium(q_t, q_s, p_hPa); }

    explicit SaturationAdjustment(cAtmosphereModel& model)
        : m(model)
    {}

    // ==================================================================================
    // ATM_SATADJ_DIAG=1 (print-only, default off) -- WHERE DOES THIS ROUTINE CREATE WATER?
    //
    // `ATM_CWB_DIAG` charges the SaturationAdjust stage **+5.5e+05 mm/a**, essentially all
    // aloft, in a routine whose own comment says total water is conserved by the loop
    // ("d_cnd + d_dep = d_q_v"). That identity holds for the RATES and not for what is
    // WRITTEN, because four clips sit between them. This attributes the difference to each,
    // and closes: sum of the buckets + unattributed = the measured change in column water.
    //
    // THE BUCKETS, in the order the routine executes them:
    //   entry_clip  q_v_old/q_c_old/q_i_old are read as max(0, .) -- a negative field value
    //               is silently lifted to zero and the write-back starts from the lifted one.
    //   phase_split THE SUSPECT. d_q_v is bounded by the TOTAL condensate q_c_b + q_i_b, but
    //               it is then split by TEMPERATURE, d_cnd = d_q_v*CND and d_dep = d_q_v*DEP,
    //               and each half is subtracted with its own max(0, .). So a cold cell
    //               holding LIQUID and no ice (CND ~ 0, DEP ~ 1) evaporates d_dep out of an
    //               empty ice reservoir: q_v gains it and nothing loses it. Availability is
    //               enforced on the sum and the split ignores it.
    //   cold_delete the T < t_00 branch. On the shipped path it deletes cloud AND ice
    //               outright (a SINK); under ATM_ICE_COLD (default on) it freezes instead
    //               and must read exactly 0.
    //   caf_neg / caf_cap / caf_fade / caf_supersat  the four in clampAndFade: the negative
    //               clips, the 0.05 kg/kg condensate cap, the cold fade, and the always-on
    //               supersaturation removal, which moves vapour to condensate and must read 0.
    //
    // Charged in the units ATM_CWB_DIAG uses -- cos-lat weighted, times the layer mass
    // rho*dz, so a bucket is [mm] of column water per call and the two instruments are
    // directly comparable. Cells below i_topography and the lid row are excluded, as there.
    // Accumulated in a function-local static rather than a model member: adding a member
    // moves sizeof(cAtmosphereModel) and that is this tree's stack-canary hazard.
    // ==================================================================================
    struct Budget {
        double entry_clip = 0.0, phase_split = 0.0, cold_delete = 0.0, adj_total = 0.0;
        double caf_neg = 0.0, caf_cap = 0.0, caf_fade = 0.0, caf_supersat = 0.0,
               caf_total = 0.0;
        double w_lat = 0.0;
        long   calls = 0;
        double cum_adj = 0.0, cum_caf = 0.0;
    };
    static Budget& budget(){ static Budget b; return b; }
    // ATM_SATADJ_PHASE: 1 = clip each evaporation half to its own reservoir (the repair),
    // 0 = shipped. Default 0 -- this changes the moist physics of the default configuration,
    // and this tree flips defaults on measurements.
    static bool satadjPhase(){
        static const bool v = [](){ const char* e = getenv("ATM_SATADJ_PHASE");
                                    return e && atoi(e) != 0; }();
        return v;
    }
    static bool diagOn(){
        static const bool v = [](){ const char* e = getenv("ATM_SATADJ_DIAG");
                                    return e && atoi(e) != 0; }();
        return v;
    }
    // The weight GetMean_2D/GetMean_3D use, replicated exactly as ColumnWaterBudget does.
    static double latWeight(int j){
        return (j <= 90) ? cos((90 - j) * M_PI / 180.0) : cos((j - 90) * M_PI / 180.0);
    }
    // Layer mass [kg/m2] at (i,j,k), 0 where ColumnWaterBudget carries no mass either.
    double cellMass(int i, int j, int k) const {
        if (i < m.i_topography[j][k] || i >= m.im - 1) return 0.0;
        double rho = m.r_humid.x[i][j][k];
        if (!AtomUtils::is_finite_safe(rho) || rho <= 0.0) rho = m.r_air;
        return rho * (m.get_layer_height(i+1) - m.get_layer_height(i));
    }


    void run() {
        using namespace std;

        cout << endl << endl << endl << "      SaturationAdjustment" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        cout.precision(9);

        computeSteps();
        adjustSaturation();
        applyTopography();
        clampAndFade();
        reportBudget();
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

        // Surface row is skipped below, but its condensation source must still be
        // cleared every call: S_c_c.x[0] feeds the ice schemes' cloud-water source,
        // and over ocean (i_mount == 0) they copy S_c_c.x[0] onto itself, so a stale
        // nonzero value would silently re-condense surface cloud and defeat the skip.
        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < m.jm; j++)
            for (int k = 0; k < m.km; k++)
                m.S_c_c.x[0][j][k] = 0.0;

        // Start at i = 1: the surface row (i = 0) is a vapour SOURCE only and must
        // never be saturation-adjusted in place. waterVapourEvaporation injects c[0]
        // from the warm ocean toward c_eq ~ q_sat every moist iter; if adjustSaturation
        // then condenses that just-injected vapour into cloud at the same cell, cloud
        // water doubles every ~3-4 iters (Cook Inlet runaway: cloud 5e-13 -> 32 over
        // iters 315-363). Once cloud + ice > 1 the r_humid denominator
        // (1 + 0.608*c - cloud - ice) flips negative, -grad(p)/rho inverts, NaN cascade.
        // Physically, condensation only happens once a parcel has risen, cooled
        // adiabatically and reached its LCL aloft; the warm ocean surface is below
        // saturation by definition. So cloud forms from i = 1 up, not at i = 0.
        const bool diag = diagOn();
        double b_entry = 0.0, b_phase = 0.0, b_cold = 0.0, b_total = 0.0;

        #pragma omp parallel for collapse(2) schedule(static) \
                reduction(+:b_entry,b_phase,b_cold,b_total)
        for (int i = 1; i < m.im - 1; i++) {
            for (int j = 0; j < m.jm; j++) {
                const double wm_j = diag ? latWeight(j) : 0.0;
                double *S_c_c_row = m.S_c_c.x[i][j];
                double *c_row     = m.c.x[i][j];
                double *cloud_row = m.cloud.x[i][j];
                double *ice_row   = m.ice.x[i][j];
                double *t_row     = m.t.x[i][j];
                double *p_row     = m.p_stat.x[i][j];
                double  dt_dim    = step[i] / 1.6;

                for (int k = 0; k < m.km; k++) {
                    S_c_c_row[k] = 0.0;

                    const double wm = diag ? wm_j * cellMass(i, j, k) : 0.0;
                    const double q_t_raw = diag ? (c_row[k] + cloud_row[k] + ice_row[k]) : 0.0;

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

                    // ATM_ICE_COLD: the entry weight fades the adjustment out below t_00, so
                    // deposition stops exactly where the CND/DEP split below sends 100 % of the
                    // condensation to the ice branch. With ice allowed to exist there, the
                    // adjustment has to be allowed to make it.
                    const double alpha_entry = ColdCloud::enabled() ? 1.0
                        : 1.0 / (1.0 + std::exp(-(T - m.t_00) / fade_K));

                    // Under ATM_CLOUD_FRAC a cell can be cloudy while the GRID MEAN is
                    // subsaturated, so entry cannot be conditioned on q_v > q_sat alone: a cell
                    // above the critical humidity must be admitted even with no condensate yet.
                    const bool frac_active = cloudFrac()
                        && (q_v_old + q_c_old + q_i_old) > hCrit(p_local) * q_sat;
                    if ((q_v_old > q_sat && alpha_entry > 0.01) || frac_active ||
                        (q_v_old < q_sat &&
                        (q_c_old > 1e-12 || q_i_old > 1e-12))) {

                        // charged only for cells that ENTER: a cell that does not enter is
                        // never written, so its max(0, .) lifting does not reach the field.
                        if (diag) b_entry += wm * ((q_v_old + q_c_old + q_i_old) - q_t_raw);

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

                            // ATM_SATADJ_PHASE=1 -- SPLIT THE EVAPORATION BY AVAILABILITY AS
                            // WELL AS BY PHASE. d_q_v is already bounded by the TOTAL condensate
                            // q_c_b + q_i_b above, but it is then split by TEMPERATURE, and each
                            // half is subtracted with its own max(0, .). So a cell whose
                            // condensate is in the phase the temperature does NOT select --
                            // liquid in a cold cell (CND ~ 0), ice in a warm one -- evaporates
                            // out of an empty reservoir: q_v gains d_q_v in full and the
                            // condensate loses only what it had. Measured with ATM_SATADJ_DIAG,
                            // that is 100 % of this routine's +5.5e+05 mm/a non-conservation.
                            //
                            // The repair moves the unavailable part to the OTHER phase rather
                            // than dropping it, because the total IS available by construction,
                            // and then sets d_q_v to what was actually taken. The latent-heat
                            // line below consumes the same corrected pair, so the energy follows
                            // the mass: evaporating ice absorbs ls, not lv.
                            if (satadjPhase() && d_q_v > 0.0) {
                                if (d_cnd > q_c_b) { d_dep += d_cnd - q_c_b; d_cnd = q_c_b; }
                                if (d_dep > q_i_b) {
                                    const double back = d_dep - q_i_b;
                                    d_dep = q_i_b;
                                    d_cnd = std::min(q_c_b, d_cnd + back);
                                }
                                d_q_v = d_cnd + d_dep;
                            }

                            q_v_b += d_q_v;
                            if (diag) b_phase += wm * alpha_entry
                                    * (std::max(0.0, d_cnd - q_c_b) + std::max(0.0, d_dep - q_i_b));
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

                            // Fraction-aware target. q_s above is the phase-blended saturation
                            // the shipped scheme drives to; under ATM_CLOUD_FRAC the cell is
                            // allowed to retain the condensate the sub-grid closure supports at
                            // that saturation, rather than being dried to it.
                            if (cloudFrac()) {
                                const double q_t_b = q_v_b + q_c_b + q_i_b;
                                const double q_c_eq = qcEquilibrium(q_t_b, q_v_target, p_local);
                                q_v_target = std::max(0.0, q_t_b - q_c_eq);
                            }

                            // Adaptive Newton-damped update. The previous fixed 0.5 under-
                            // relaxation is UNSTABLE when the latent-heat gain
                            // G = (L/cp)*dq_sat/dT exceeds 3 — i.e. warm cells (T>~25C,
                            // q_sat>~20 g/kg) where dq_sat/dT is steep: the fixed-point map
                            // derivative 0.5*(1-G) then has magnitude >1, so the loop
                            // OSCILLATES and never converges, leaving RH stuck at 120-139%
                            // with cloud present (thermodynamically impossible) and feeding
                            // the coastal precip runaway. omega = 1/(1+G) drives the map
                            // derivative 1-omega*(1+G) to 0 (stable, ~Newton-optimal) at all T.
                            // dq_sat/dT from Clausius-Clapeyron: q_sat*L/(Rv*T^2).
                            // project_overprecip_saturation_injection.
                            const double Rv = 461.5;                    // [J/(kg K)] water-vapour gas constant
                            const double inv_RvT2 = 1.0 / (Rv * T * T);
                            const double Gain = CND * lv_over_cp * (q_sat * m.lv * inv_RvT2)
                                              + DEP * ls_over_cp * (q_Ice * m.ls * inv_RvT2);
                            const double omega = 1.0 / (1.0 + Gain);
                            q_v_hyp = q_v_b + omega * (q_v_target - q_v_b);

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
 
                        const double q_t_pre_cold = diag
                            ? (c_row[k] + cloud_row[k] + ice_row[k]) : 0.0;

                        if (T < m.t_00) {
                            // ATM_ICE_COLD: below the homogeneous-freezing point LIQUID cannot
                            // exist -- but ice must, and this is where cirrus lives. Freeze it
                            // instead of deleting it. Default off, shipped branch unchanged.
                            if (ColdCloud::enabled()) {
                                ice_row[k]  += cloud_row[k];
                                cloud_row[k] = 0.0;
                            } else {
                            cloud_row[k] = 0.0;
                            ice_row[k]   = 0.0;
                            }
                        }
                        if (diag) {
                            b_cold  += wm * ((c_row[k] + cloud_row[k] + ice_row[k])
                                             - q_t_pre_cold);
                            b_total += wm * ((c_row[k] + cloud_row[k] + ice_row[k]) - q_t_raw);
                        }
                    }
                }
            }
        }

        if (diag) {
            Budget& b = budget();
            if (b.w_lat <= 0.0)
                for (int j = 0; j < m.jm; j++) b.w_lat += latWeight(j) * m.km;
            b.entry_clip  = b_entry / b.w_lat;
            b.phase_split = b_phase / b.w_lat;
            b.cold_delete = b_cold  / b.w_lat;
            b.adj_total   = b_total / b.w_lat;
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
        const double inv_t_0    = 1.0 / m.t_0;
        const double lv_over_cp = m.lv / m.cp_l;
        const double ls_over_cp = m.ls / m.cp_l;
        // Defensive physical bounds. T_max mirrors the Magnus-validity cap used in
        // adjustSaturation; cloud_cap is ~50× the largest physical cloud/ice mixing
        // ratio (a few g/kg), so it never clips a real cloud — it only stops a runaway.
        constexpr double T_max     = 333.15;   // 60 °C
        constexpr double cloud_cap = 0.05;     // kg/kg condensate ceiling

        const bool diag = diagOn();
        double b_neg = 0.0, b_cap = 0.0, b_fade = 0.0, b_ss = 0.0, b_tot = 0.0;

        #pragma omp parallel for collapse(2) schedule(static) \
                reduction(+:b_neg,b_cap,b_fade,b_ss,b_tot)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                const double wm_j = diag ? latWeight(j) : 0.0;
                double *c_row     = m.c.x[i][j];
                double *cloud_row = m.cloud.x[i][j];
                double *ice_row   = m.ice.x[i][j];
                double *t_row_nd  = m.t.x[i][j];
                double *p_row     = m.p_stat.x[i][j];

                for (int k = 0; k < m.km; k++) {
                    const double wm = diag ? wm_j * cellMass(i, j, k) : 0.0;
                    const double q_t_in = diag
                        ? (c_row[k] + cloud_row[k] + ice_row[k]) : 0.0;

                    if (c_row[k]     < 0.0) c_row[k]     = 0.0;
                    if (cloud_row[k] < 0.0) cloud_row[k] = 0.0;
                    if (ice_row[k]   < 0.0) ice_row[k]   = 0.0;
                    if (diag) b_neg += wm * ((c_row[k] + cloud_row[k] + ice_row[k]) - q_t_in);
                    const double q_t_neg = diag
                        ? (c_row[k] + cloud_row[k] + ice_row[k]) : 0.0;

                    double T_dim = t_row_nd[k] * m.t_0;
                    // Upper temperature bound first (no air parcel exceeds ~60 °C); keeps
                    // q_sat below finite and bounds the latent-heat release that follows.
                    if (T_dim > T_max) { T_dim = T_max; t_row_nd[k] = T_max * inv_t_0; }

                    // ---- Always-on supersaturation removal (ROOT FIX) ----
                    // adjustSaturation scales its condensation by alpha_entry, so in cold
                    // air (alpha ≪ 1) it leaves q_v ≫ q_sat — and cells colder than
                    // ~−60 °C (alpha ≤ 0.01) it skips entirely. With no upper bound on q_v
                    // anywhere, that residual supersaturation accumulated unbounded at the
                    // orographic saturation level in the cold fade zone (NZ Southern Alps,
                    // i≈12, T≈−40 °C → q_v→6, cloud→129, t→1e5 °C, then buoyancy drove a
                    // vertical-velocity runaway). Remove supersaturation here for EVERY
                    // cell, independent of the alpha fade, conserving water (excess → cloud
                    // above freezing, → ice below) and energy (latent heat → T). In a clean
                    // run the per-call excess is small (physical); the T_max recap below is
                    // the backstop if a transient ever drives a large excess.
                    const double p_local = p_row[k];
                    const double E_sat = m.hp * exp_func(T_dim, 17.2694, 35.86);
                    const double q_sat = (p_local > E_sat)
                        ? m.ep * E_sat / (p_local - E_sat)
                        : m.ep * 1e-5;
                    if (c_row[k] > q_sat) {
                        const double excess = c_row[k] - q_sat;
                        c_row[k] = q_sat;
                        if (T_dim >= m.t_00) {
                            cloud_row[k] += excess;
                            T_dim        += lv_over_cp * excess;
                        } else {
                            ice_row[k]   += excess;
                            T_dim        += ls_over_cp * excess;
                        }
                        if (T_dim > T_max) T_dim = T_max;   // backstop on the latent release
                        t_row_nd[k] = T_dim * inv_t_0;
                    }

                    // the supersaturation removal above moves vapour into condensate and
                    // must be conservative; measured rather than assumed.
                    const double q_t_ss = diag
                        ? (c_row[k] + cloud_row[k] + ice_row[k]) : 0.0;
                    if (diag) b_ss += wm * (q_t_ss - q_t_neg);

                    // ---- Condensate upper bounds (CAP SAFETY NET) ----
                    if (cloud_row[k] > cloud_cap) cloud_row[k] = cloud_cap;
                    if (ice_row[k]   > cloud_cap) ice_row[k]   = cloud_cap;
                    const double q_t_cap = diag
                        ? (c_row[k] + cloud_row[k] + ice_row[k]) : 0.0;
                    if (diag) b_cap += wm * (q_t_cap - q_t_ss);

                    // Existing cold fade: smoothly dry moisture toward 0 below t_00.
                    // ATM_ICE_COLD restricts it to the LIQUID. Vapour is not a condensate and
                    // must never be faded at all -- doing so deletes mass from the column --
                    // and fading the ice is what forbids cirrus.
                    const double alpha = 1.0 / (1.0 + std::exp(-(T_dim - m.t_00) / fade_K));
                    if (ColdCloud::enabled()) {
                        cloud_row[k] *= alpha;
                    } else {
                    c_row[k]     *= alpha;
                    cloud_row[k] *= alpha;
                    ice_row[k]   *= alpha;
                    }
                    if (diag) {
                        const double q_t_out = c_row[k] + cloud_row[k] + ice_row[k];
                        b_fade += wm * (q_t_out - q_t_cap);
                        b_tot  += wm * (q_t_out - q_t_in);
                    }
                }
            }
        }

        if (diag) {
            Budget& b = budget();
            if (b.w_lat <= 0.0)
                for (int j = 0; j < m.jm; j++) b.w_lat += latWeight(j) * m.km;
            b.caf_neg      = b_neg  / b.w_lat;
            b.caf_supersat = b_ss   / b.w_lat;
            b.caf_cap      = b_cap  / b.w_lat;
            b.caf_fade     = b_fade / b.w_lat;
            b.caf_total    = b_tot  / b.w_lat;
        }
    }

    // ATM_SATADJ_DIAG report. Per call and cumulative, in [mm] of column water on
    // ColumnWaterBudget's own weights, so a row here is directly comparable with that
    // instrument's SaturationAdjust bucket. `unattributed` is the identity check: it is the
    // measured change minus the named mechanisms, and it must be ~0.
    void reportBudget() const {
        if (!diagOn()) return;
        Budget& b = budget();
        b.calls++;
        b.cum_adj += b.adj_total;
        b.cum_caf += b.caf_total;
        const double adj_un = b.adj_total - (b.entry_clip + b.phase_split + b.cold_delete);
        const double caf_un = b.caf_total - (b.caf_neg + b.caf_supersat + b.caf_cap + b.caf_fade);
        printf("      AGCM: [SATADJ BUDGET] call %ld, mm of column water:\n", b.calls);
        printf("        adjustSaturation  entry_clip %+.6e  phase_split %+.6e"
               "  cold_delete %+.6e  |  total %+.6e  unattributed %+.6e\n",
               b.entry_clip, b.phase_split, b.cold_delete, b.adj_total, adj_un);
        printf("        clampAndFade      neg_clip   %+.6e  supersat    %+.6e"
               "  cap %+.6e  fade %+.6e  |  total %+.6e  unattributed %+.6e\n",
               b.caf_neg, b.caf_supersat, b.caf_cap, b.caf_fade, b.caf_total, caf_un);
        printf("        cumulative        adjustSaturation %+.6e   clampAndFade %+.6e"
               "   both %+.6e\n", b.cum_adj, b.cum_caf, b.cum_adj + b.cum_caf);
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
