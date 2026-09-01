#pragma once

#include "cAtmosphereModel.h"
#include "IceSchemeCommon.h"
#include "CloudFraction.h"

#include <algorithm>
#include <cmath>
#include <chrono>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>

using namespace AtomUtils;


namespace ThreeCatIce {
    // ice particle parameters
    constexpr double N_i_0 = 1.0e2;                                    // 1/m3
    constexpr double N_g_0 = 4.0e6;                                    // 1/m3
    constexpr double m_i_0 = 1.0e-12;                                  // kg
    constexpr double m_i_max = 1.0e-9;                                  // kg
    constexpr double m_s_0 = 3.0e-9;                                    // kg
    constexpr double N_r_0 = 8.0e+6;                                   // 1/m3
    constexpr double N_s_0 = 8.0e+5;                                   // 1/m3

    // fall velocity prefactors
    constexpr double v_s_0 = 130.0;
    constexpr double v_r_0 = 4.9;
    constexpr double v_g_0 = 442.0;

    // microphysical rate coefficients
    // Cloud->rain autoconversion rate. Raised from the COSMO 4.0e-4 to 1.0e-3 as the rain-side
    // partner of the reduced snow riming (c_rim_snow) -- PORTED FROM TwoCatIceScheme (`c4d1bd1`,
    // 2026-07-07), which ThreeCat never received: Zero/One/ThreeCat were frozen at the
    // 2026-04-14 initial commit and only TwoCat carries the tuning. Reducing riming frees cloud
    // water that must reach the ground as RAIN rather than accumulate as cloud, so the
    // cloud->rain drain is sped up to match. project_snow_overproduction.
    constexpr double c_c_au = 1.0e-3;                                  // 1/s (COSMO original 4.0e-4)
    constexpr double c_i_au = 1.0e-3;                                  // 1/s COSMO
    constexpr double c_ac = 0.24;                                       // m2/kg
    constexpr double c_rim = 18.6;                                      // m2/kg
    // Reduced SNOW-side riming coefficient, ported from TwoCatIceScheme (`c4d1bd1`). Un-reduced
    // riming (S_s_rim = c_rim*cloud*Snow) dominated snow production ~10x every other source and
    // monopolised the proportional cloud-water limiter, routing ~90 % of cloud water into SNOW
    // and starving the autoconversion->rain pathway. In TwoCat it took the snow fraction from
    // ~45 % to ~16 % at fixed total precipitation. Measured here before the port, from the
    // accepted configuration's iteration-600 checkpoint: P_snow 2470 mm/a against TwoCat's 0.29
    // and a total of 3749 against NASA's 978.
    // WARM-SIDE SHEDDING (S_s_shed) KEEPS THE FULL c_rim, exactly as in TwoCat, and so does the
    // GRAUPEL riming S_g_rim -- TwoCat has no graupel, so the port has nothing to say about it.
    constexpr double c_rim_snow = 18.6 / 5.0;                           // m2/kg (reduced snow riming)
    // Reduced GRAUPEL-side riming, the same treatment applied to S_g_rim. The snow port left
    // graupel at 919 mm/a -- 94 % of NASA's entire precipitation from one species -- because
    // TwoCat has no graupel and `c4d1bd1` therefore says nothing about it. The mechanism is the
    // same one: an un-reduced riming term monopolises the proportional cloud-water limiter and
    // starves the autoconversion->rain pathway, so the same factor of 5 applies.
    //
    // AND THE COLD SIDE WAS USING THE SNOW CONSTANT, WHICH IS FIXED HERE TOO. `c_g_rim` = 4.43
    // exists, is named for graupel riming, and was used ONLY by the warm-side shedding
    // S_g_shed; the cold-side S_g_rim reached for the SNOW constant c_rim = 18.6, 4.2x larger.
    // The two corrections were measured SEPARATELY rather than together, so neither is
    // confounded with the other: c_rim/5 = 3.72 first (graupel 919 -> 578), then c_g_rim/5 =
    // 0.886 (see CLAUDE.md for the second number).
    constexpr double c_rim_graupel = 4.43 / 5.0;                        // m2/kg = c_g_rim/5
    constexpr double c_agg = 10.3;                                      // m2/kg
    constexpr double c_i_cri = 0.24;                                    // m2
    constexpr double c_r_cri = 3.2e-5;                                  // m2
    constexpr double b_ev = 5.9;                                        // m2*s/kg
    constexpr double c_r_frz = 3.75e-2;                                 // m2/(K*kg)
    constexpr double c_i_dep = 1.3e-5;                                  // m3/(s*kg)

    // graupel-specific coefficients
    constexpr double z_csg = 0.5;
    constexpr double c_g_rim = 4.43;
    constexpr double c_g_agg = 2.46;

    // mass-size relation prefactors
    constexpr double a_s_m = 0.038;
    constexpr double a_g_m = 169.6;

    // temperature thresholds
    constexpr double t_nuc = 267.15;                                    // K  (-6 C)
    constexpr double t_d = 248.15;                                      // K  (-25 C)
    constexpr double t_hn = 236.15;                                     // K  (-37 C)
    constexpr double t_r_frz = 271.15;                                  // K  (-2 C)

    // iteration / conversion
    constexpr double conv_mmd = 8.64e4;                                 // s/d conversion to mm/d
    constexpr int iter_prec_end = 2;                                    // COSMO iterations
}


// Three-Category-Ice-Scheme, COSMO-module from the German Weather Forecast,
// yielding the precipitation distribution from rain, snow, and graupel.
class ThreeCatIceScheme {
public:
    explicit ThreeCatIceScheme(cAtmosphereModel& model)
        : m(model)
    {}

    void run() {
        using namespace std;
        using namespace ThreeCatIce;

        cout << endl << endl << endl << "      ThreeCategoryIceScheme" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        initArrays();
        computeColumns();
        applyBoundaryConditions();
        applyTopography();

        printReport();

        cout << "      ThreeCategoryIceScheme ended" << endl;

        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for ThreeCategoryIceScheme\n", elapsed.count() * 1e-9);
    }

private:
    cAtmosphereModel& m;

    // diagnostic output state (set during computeColumns)
    bool rain = false;
    bool snow = false;
    bool graupel = false;
    double global_max_rain = 0.0;
    double global_max_snow = 0.0;
    double global_max_graupel = 0.0;
    int i_rain = 0, j_rain = 0, k_rain = 0;
    int i_snow = 0, j_snow = 0, k_snow = 0;
    int i_graupel = 0, j_graupel = 0, k_graupel = 0;

    // layer geometry tables (filled in initArrays, read in computeColumns/printReport)
    std::vector<double> step_table;
    std::vector<double> height_table;


    // ==================== INIT ====================
    void initArrays() {
        // Precompute layer thickness table
        step_table.assign(m.im, 0.0);
        height_table.resize(m.im);
        for(int i = 0; i < m.im; i++){
            height_table[i] = m.get_layer_height(i);
            if(i < m.im - 1)
                step_table[i] = m.get_layer_height(i + 1) - m.get_layer_height(i);
        }

        // Clamping + copy pass
        #pragma omp parallel for collapse(2)
        for(int j = 0; j < m.jm; j++){
            for(int k = 0; k < m.km; k++){
                #pragma omp simd
                for(int i = 0; i < m.im; i++){
                    m.c.x[i][j][k]     = std::max(0.0, m.c.x[i][j][k]);
                    m.cloud.x[i][j][k] = std::max(0.0, m.cloud.x[i][j][k]);
                    m.ice.x[i][j][k]   = std::max(0.0, m.ice.x[i][j][k]);
                    m.gr.x[i][j][k]    = std::max(0.0, m.gr.x[i][j][k]);
                    m.P_rain.x[i][j][k] = m.P_rainn.x[i][j][k];
                    m.P_snow.x[i][j][k] = m.P_snown.x[i][j][k];
                }
            }
        }
    }


    // ==================== ATM_SS_DIAG: the snow / graupel / rain flux budgets ====================
    // ATM_SS_DIAG=1 -- print the COLUMN BUDGET of all three precipitation fluxes, decomposed
    // into every source and sink term that feeds them. Print-only, default off, no field is
    // written and no physics changes.
    //
    // The scheme integrates three fluxes down the column,
    //
    //     P_x[i] = min(P_max_flux, max(0, P_x[i+1] + rho[i+1]*S_x[i+1]*dz_i))       (x = s, g, r)
    //
    // and only the SUM S_x is written to any output file, so the surface value cannot say which
    // of its eight or nine terms decided it. That is the open question left by the riming ports:
    // snow is 1046 mm/a, 39 % of a total precipitation 2.74x NASA, against Earth's 5-10 %, and
    // the riming fix that took TwoCat's snow fraction to 16 % has ALREADY been applied here. So
    // the residual snow is coming from terms TwoCat does not have -- `S_s_agg`, `S_s_dep`, the
    // graupel/snow exchange `S_csg` -- and guessing which one costs a run per coefficient.
    //
    // As with ATM_MC_DIAG and ATM_SR_DIAG, a raw sum of the terms is NOT a budget: the
    // max(0, ...) floor, the P_max_flux cap and the temperature window each truncate the flux at
    // every level. The truncation is therefore charged to its own bucket, signed, and then
    //
    //     sources - sinks + clamp = P_x(ground)
    //
    // is an IDENTITY that closes against the model's own field to every printed digit. The
    // residual is printed too, so a broken decomposition announces itself.
    //
    // Two buckets have no counterpart in the TwoCat print. `stale top` is the very first
    // increment, P_x[im-2] <- rho*S_x[im-1]*dz, which reads S_x at the LID -- a value this call
    // never wrote, left behind by the previous one. And the budget is accumulated only down to
    // `i_topography`, so the ground value is the physical surface rather than sea level.
    // ATM_ICE_RAW_FLUX=1 -- feed the precipitation fluxes to the microphysics in SI units
    // instead of normalised by the column's own surface flux. Default 0 = shipped.
    //
    // ThreeCat divides every flux by the value at the ground before using it,
    //
    //     Rain = P_rain[i] / max(P_rain[0], 1e-6),
    //
    // and then hands that RATIO to relations that are dimensional in kg/(m2 s). The mass
    // content of falling rain, `r_q_r = A_r*B_rad^(-8/9)*Rain^(8/9)`, is the clearest case: with
    // the raw flux, a typical 1e-5 kg/(m2 s) gives r_q_r = 7.8e-5 kg/m3, an ordinary 0.08 g/m3
    // of rain water. With the ratio, which the 1e-6 floor lets reach P_max_flux/1e-6 = 3000, the
    // same expression returns ~2.5e3 kg/m3 -- denser than liquid water. Every collection term is
    // a power of these quantities, so all of them inflate together.
    //
    // TWOCAT USES THE RAW FLUX (`TwoCatIceScheme.h:241`, `double Rain = m.P_rain.x[i][j][k];`)
    // and so does the COSMO documentation the constants come from; the normalisation is
    // ThreeCat's alone and has been there since the initial commit. It is also the reason the
    // scheme once produced NaN: `P_snow_0` in the denominator with `S_s_rim` proportional to
    // `Snow` is a geometric amplifier down the column, which is what the `P_norm_floor` and the
    // `P_max_flux` cap were added to contain. Both treat the symptom.
    //
    // Measured with ATM_SS_DIAG on the shipped branch, six iterations from the accepted
    // configuration's iteration-600 checkpoint: the snow budget's gross demand is 7.35e10 mm/a
    // and the clamp removes 7.35e10, leaving 1166. The ground flux is a residual of a demand
    // seven orders of magnitude larger, so it is set by `P_max_flux` and the max(0, ...) floor
    // rather than by any microphysical rate -- which is why cutting a riming coefficient by 5
    // and then by a further 4.2 bought 37 % and then 11 %.
    static bool rawFlux(){
        static const bool v = [](){
            const char* e = getenv("ATM_ICE_RAW_FLUX"); return e && atoi(e) != 0; }();
        return v;
    }

    // ATM_ICE_LIMITERS=1 -- charge every sink against the water that is actually there.
    // Default 0 = shipped.
    //
    // TwoCat carries five availability limiters; ThreeCat has NONE, and its own comments cite a
    // "proportional cloud-water limiter" that exists only in the other file. Nothing stops a
    // ThreeCat sink from removing more cloud water, cloud ice or falling flux than the cell
    // holds: the only thing that keeps the scheme finite is the max(0, ...) floor and the
    // `P_max_flux` cap on the INTEGRATED flux, applied one level later. Measured with
    // ATM_SS_DIAG on the dimensionally-repaired branch, rain evaporation `S_ev` demands
    // 5074 mm/a against 1664 mm/a of rain sources -- three times the water present -- and the
    // floor then discards the difference.
    //
    // Ported term for term from `TwoCatIceScheme.h:463`, with the two categories TwoCat does not
    // have (graupel, and the snow<->graupel conversion) folded into the same totals, and the
    // per-species clips written against the RAW flux so they mean the same thing on both
    // branches of ATM_ICE_RAW_FLUX.
    //
    // The non-negativity of melting is part of the same repair: `melt_driver` carries a
    // `(c - q_sat(t_0))` term that goes negative in subsaturated air, and a negative `S_s_melt`
    // is an unphysical rain->snow conversion, not a melting rate.
    static bool iceLimiters(){
        static const bool v = [](){
            const char* e = getenv("ATM_ICE_LIMITERS"); return e && atoi(e) != 0; }();
        return v;
    }

    static bool ssDiag(){
        static const bool v = [](){
            const char* e = getenv("ATM_SS_DIAG"); return e && atoi(e) != 0; }();
        return v;
    }

    // Bucket layout. The first 26 are TERM values carrying the sign with which they enter their
    // own S_x, so that summing a species' range reproduces S_x exactly.
    enum {
        SS_i_au = 0, SS_d_au, SS_agg, SS_rim, SS_dep, SS_i_cri, SS_r_cri, SS_melt, SS_csg,
        SG_agg, SG_rim, SG_dep, SG_i_cri, SG_r_cri, SG_r_frz, SG_melt, SG_csg,
        SR_c_au, SR_ac, SR_ev, SR_s_shed, SR_g_shed, SR_r_cri, SR_r_frz, SR_s_melt, SR_g_melt,
        SS_END = SG_agg, SG_END = SR_c_au, SR_END = 26,
        // The clamp is three different mechanisms and they are not interchangeable: the FLOOR
        // injects mass where the net tendency would drive the flux negative, the CAP truncates
        // at P_max_flux, and the WINDOW deletes the whole flux at a phase boundary. Their sum is
        // what a single clamp bucket used to report.
        SX_floor_s = 26, SX_floor_g, SX_floor_r,
        SX_cap_s, SX_cap_g, SX_cap_r,
        SX_win_s, SX_win_g, SX_win_r,
        SX_top_s, SX_top_g, SX_top_r,
        SX_gnd_s, SX_gnd_g, SX_gnd_r,
        NSSD
    };

    // ==================== COLUMN COMPUTATION ====================
    void computeColumns() {
        using namespace ThreeCatIce;

        // Precompute gamma and composite constants
        const double gamma_4_5 = tgamma(4.5);
        const double gamma_4_0 = tgamma(4.0);
        const double gamma_3_25 = tgamma(3.25);

        const double A_r = m.r_0_water * M_PI * N_r_0;
        const double A_s = 2.0 * a_s_m * N_s_0;
        const double A_g = 2.0 * a_g_m * N_g_0;
        const double B_rad = m.r_0_water * M_PI * N_r_0 * v_r_0 * gamma_4_5 / gamma_4_0;
        const double B_s = N_s_0 * v_s_0 * a_s_m * gamma_3_25;
        const double B_g = N_g_0 * v_g_0 * a_g_m * gamma_3_25;

        // Precompute constant pow() bases — called im×jm×km times otherwise
        const double B_rad_neg89 = pow(B_rad, -8.0/9.0);
        const double B_s_neg1213 = pow(B_s, -12.0/13.0);
        const double B_g_neg1213 = pow(B_g, -12.0/13.0);

        // Precompute E_sat at t_0 (constant, used in every cell)
        const double E_sat_t_0_Pa = 1e2 * m.hp * AtomUtils::exp_func(m.t_0, 17.2694, 35.86);

        // Kessler autoconversion threshold, in-cloud [kg/kg]; ATM_QC_CRIT, default 0.05 g/kg.
        const double q_c_crit = IceSchemeCommon::qcCrit();

        // Thread-local accumulators for the OMP reduction
        double local_max_rain = 0.0, local_max_snow = 0.0, local_max_graupel = 0.0;

        // Per-column budget slots; each (j,k) is owned by one thread, so no synchronisation.
        const bool raw_flux = rawFlux();
        const bool limiters = iceLimiters();
        const bool ssd_on = ssDiag();
        std::vector<double> ssd;
        if(ssd_on) ssd.assign((size_t)NSSD * m.jm * m.km, 0.0);

        #pragma omp parallel for collapse(2) schedule(static) \
            reduction(max:local_max_rain, local_max_snow, local_max_graupel)
        for(int j = 1; j < m.jm - 1; j++){
            for(int k = 1; k < m.km - 1; k++){

                m.P_rain.x[23][j][k] = 0.0;

                // Budget state, carried DOWN the column: the flux increment at level i is
                // built from the terms at level i+1, so each level saves its own for the next.
                double pv[SR_END];
                double pv_rho = 0.0;
                bool   pv_valid = false;
                const int i_gnd = ssd_on
                    ? std::min(std::max(m.i_topography[j][k], 0), m.im - 1) : 0;

                for(int iter_prec = 1; iter_prec <= iter_prec_end; iter_prec++){

                    double Rain_check = m.P_rain.x[23][j][k];

                    if(ssd_on){
                        for(int q = 0; q < NSSD; q++)
                            ssd[(size_t)q * m.jm * m.km + (size_t)j * m.km + k] = 0.0;
                        for(int q = 0; q < SR_END; q++) pv[q] = 0.0;
                        pv_valid = false;
                    }

                    m.P_rain.x[m.im-1][j][k] = 0.0;
                    m.P_snow.x[m.im-1][j][k] = 0.0;
                    m.P_graupel.x[m.im-1][j][k] = 0.0;

                    for(int i = m.im - 2; i >= 0; i--){

                        // Normalize precipitation by the surface flux. FLOORED denominator: a
                        // tiny-but-nonzero surface flux made Snow = P_snow[i]/P_snow_0 explode,
                        // and S_s_rim ∝ Snow then amplified P_snow geometrically down the column
                        // to overflow (the ThreeCat NaN blow-up). Flooring at P_norm_floor bounds
                        // the normalized ratios; the flux cap below guarantees finiteness.
                        constexpr double P_norm_floor = 1.0e-6;      // kg/(m2*s) ~0.09 mm/d
                        double P_rain_0    = std::max(m.P_rain.x[0][j][k],    P_norm_floor);
                        double P_snow_0    = std::max(m.P_snow.x[0][j][k],    P_norm_floor);
                        double P_graupel_0 = std::max(m.P_graupel.x[0][j][k], P_norm_floor);

                        double Rain    = raw_flux ? m.P_rain.x[i][j][k]
                                                  : m.P_rain.x[i][j][k]    / P_rain_0;
                        double Snow    = raw_flux ? m.P_snow.x[i][j][k]
                                                  : m.P_snow.x[i][j][k]    / P_snow_0;
                        double Graupel = raw_flux ? m.P_graupel.x[i][j][k]
                                                  : m.P_graupel.x[i][j][k] / P_graupel_0;

                        // pow() calls guarded — skip entirely when base is zero
                        double Rain_79  = (Rain > 0.0) ? pow(Rain, 7.0/9.0)  : 0.0;
                        double Rain_89  = (Rain > 0.0) ? pow(Rain, 8.0/9.0)  : 0.0;
                        double Rain_139 = (Rain > 0.0) ? pow(Rain, 13.0/9.0) : 0.0;
                        double Rain_16  = (Rain > 0.0) ? pow(Rain, 1.0/6.0)  : 0.0;
                        double Rain_49  = (Rain > 0.0) ? pow(Rain, 4.0/9.0)  : 0.0;

                        double r_q_r = A_r * B_rad_neg89 * Rain_89;
                        double r_q_s = A_s * B_s_neg1213 * ((Snow > 0.0)    ? pow(Snow,    12.0/13.0) : 0.0);
                        double r_q_g = A_g * B_g_neg1213 * ((Graupel > 0.0) ? pow(Graupel, 12.0/13.0) : 0.0);

                        double t_u   = m.t.x[i][j][k] * m.t_0;
                        double p_u   = m.p_stat.x[i][j][k];
                        double p_u_0 = 1e2 * p_u;
                        double r_h_i = m.r_humid.x[i][j][k];
                        double c_ijk = m.c.x[i][j][k];
                        double cl_i  = m.cloud.x[i][j][k];
                        double ice_i = m.ice.x[i][j][k];

                        double step_i = step_table[i];

                        double E_sat = m.hp * AtomUtils::exp_func(t_u, 17.2694, 35.86);
                        double q_sat = m.ep * E_sat / (p_u - E_sat);
                        double E_Ice  = m.hp * AtomUtils::exp_func(t_u, 21.8746, 7.66);
                        double q_Ice  = m.ep * E_Ice / (p_u - E_Ice);

                        double dt_rain_dim = step_i / 1.6;
                        double dt_snow_dim = step_i / 0.96;

                        // Precompute shared pow() terms for snow/graupel deposition/melting
                        double rh_rqs_08 = (r_q_s > 0.0) ? pow(r_h_i * r_q_s, 0.8)     : 0.0;
                        double rh_rqg_06 = (r_q_g > 0.0) ? pow(r_h_i * r_q_g, 0.6)     : 0.0;
                        double rh_rqg_95 = (r_q_g > 0.0) ? pow(r_h_i * r_q_g, 0.94878) : 0.0;

                        // --- Ice particle properties + vapour->ice->snow throttle (shared) ---
                        IceSchemeCommon::IceSnowRates thr =
                            IceSchemeCommon::depositionThrottle(m, i, j, k, t_u, q_Ice, dt_snow_dim);
                        double m_i = thr.m_i;
                        double N_i = thr.N_i;

                        // --- Nucleation ---
                        double S_nuc = 0.0;
                        if(ice_i == 0.0){
                            if((t_u < t_d && c_ijk >= q_Ice)
                                || (t_d <= t_u && t_u <= t_nuc && c_ijk >= q_sat))
                                S_nuc = m_i_0 / (r_h_i * dt_snow_dim) * N_i;
                        }

                        double S_c_frz = (t_u < t_hn && cl_i > 0.0)
                            ? cl_i / dt_rain_dim : 0.0;

                        // --- Deposition growth of cloud ice (shared throttle) ---
                        double S_i_dep = thr.S_i_dep;

                        // --- Autoconversion, IN-CLOUD (ported from `3ea78e3`) ---
                        // The shipped test was `c_c_au*(q_c - 0.0002)` on the GRID MEAN, a
                        // hardcoded 0.2 g/kg threshold. TwoCat, OneCat and ZeroCat were all made
                        // fraction-aware in `3ea78e3` and ThreeCat, then not the default scheme,
                        // was missed. Autoconversion is a LOCAL process, so the grid-mean
                        // tendency is `f*R(q_c/f)`; only the thresholded, nonlinear terms move
                        // under that transform. Under the now-default `ATM_CLOUD_FRAC` the
                        // grid-mean cloud water peaks near 0.06 g/kg, so `q_c - 0.0002` is
                        // negative almost everywhere and this scheme's autoconversion was
                        // effectively OFF -- its rain came from accretion and melting alone.
                        const double q_t_frac = c_ijk + std::max(0.0, cl_i)
                                              + std::max(0.0, ice_i);
                        const double f_cld = CloudFraction::effectiveFraction(
                                q_t_frac, q_sat, p_u,
                                std::max(0.0, cl_i) + std::max(0.0, ice_i));
                        const double cloud_in = cl_i / f_cld;            // in-cloud water
                        double S_c_au = (t_u >= m.t_0 && cloud_in > q_c_crit)
                            ? f_cld * c_c_au * (cloud_in - q_c_crit) : 0.0;

                        double S_i_au = thr.S_i_au;   // shared throttle (ice aggregation -> snow)
                        double S_d_au = thr.S_d_au;   // shared throttle (ice deposition -> snow)

                        // --- Collection ---
                        // Accretion, in-cloud and above the same reservoir threshold, exactly
                        // as TwoCat carries it: only cloud water in excess of `q_c_crit` is
                        // collected by falling rain, so the reservoir that autoconversion cannot
                        // drain is not drained by accretion round the back either.
                        double S_ac = (t_u >= m.t_0 && cloud_in > q_c_crit)
                            ? f_cld * c_ac * (cloud_in - q_c_crit) * Rain_79 : 0.0;

                        double S_s_rim, S_g_rim, S_s_shed, S_g_shed;
                        if(t_u < m.t_0){
                            S_s_rim  = c_rim_snow * cl_i * Snow;   // reduced riming (see c_rim_snow)
                            S_g_rim  = c_rim_graupel * cl_i * rh_rqg_95;  // reduced riming
                            S_s_shed = 0.0;
                            S_g_shed = 0.0;
                        }else{
                            S_s_rim  = 0.0;
                            S_g_rim  = 0.0;
                            S_s_shed = c_rim * cl_i * Snow;
                            S_g_shed = c_g_rim * cl_i * rh_rqg_95;
                        }

                        double S_s_agg, S_g_agg, S_i_cri, S_r_cri;
                        if(t_u <= m.t_0){
                            S_s_agg = c_agg * ice_i * Snow;
                            S_g_agg = c_g_agg * ice_i * rh_rqg_95;
                            S_i_cri = c_i_cri * ice_i * Rain_79;
                            S_r_cri = c_r_cri * ice_i / m_i * Rain_139;
                        }else{
                            S_s_agg = 0.0;
                            S_g_agg = 0.0;
                            S_i_cri = 0.0;
                            S_r_cri = 0.0;
                        }

                        // --- Evaporation ---
                        double S_ev = 0.0;
                        if(t_u >= m.t_0 && Rain > 0.0){
                            double a_ev = 2.76e-3 * exp(0.055 * (m.t_0 - t_u));
                            S_ev = a_ev * (1.0 + b_ev * Rain_16)
                                 * (q_sat - c_ijk) * Rain_49;
                        }

                        // --- Deposition/sublimation of snow and graupel ---
                        double S_s_dep = 0.0, S_g_dep = 0.0;
                        double inv_p_u_0 = 1.0 / p_u_0;
                        double c_minus_qIce = c_ijk - q_Ice;

                        // Compute d_v, l_h, t_crit for melting branch
                        double d_v = (t_u >= m.t_0)
                            ? 101325.0 * inv_p_u_0 * (2.22e-5 + 1.46e-7 * (t_u - m.t_0))
                            : 101325.0 * inv_p_u_0 * (2.22e-5 + 1.25e-7 * (t_u - m.t_0));
                        double l_h = 0.024 + 8.0e-5 * (t_u - m.t_0);
                        double q_sat_t_0 = m.ep * E_sat_t_0_Pa / (p_u_0 - E_sat_t_0_Pa);

                        double t_crit = (c_ijk >= q_sat_t_0)
                            ? m.t_0 - (1.0 / l_h) * m.ls * d_v * r_h_i * (c_ijk - q_sat_t_0)
                            : m.t_0;

                        if(t_u < m.t_0){
                            S_s_dep = (2.91955 - 0.0109928*t_u + 15871.3*inv_p_u_0 + 1.74744e-6*p_u_0)
                                    * c_minus_qIce * rh_rqs_08;
                            S_g_dep = 0.0;  // as in original
                        }else if(t_u < t_crit){
                            S_s_dep = (0.28003 - 0.146293e-6*p_u_0)
                                    * c_minus_qIce * rh_rqs_08;
                            S_g_dep = (0.0418521 - 4.7524e-8*p_u_0)
                                    * c_minus_qIce * rh_rqg_06;
                        }else{
                            S_s_dep = (2.41897 + 31282.3*inv_p_u_0)
                                    * c_minus_qIce * rh_rqs_08;
                            S_g_dep = (0.153907 - 7.86703e-7*p_u_0)
                                    * c_minus_qIce * rh_rqg_06;
                        }

                        // --- Melting ---
                        double S_s_melt = 0.0, S_g_melt = 0.0;
                        if(t_u >= m.t_0){
                            constexpr double a_melt = 2.95e+3;
                            double melt_driver = (t_u - m.t_0) + a_melt * (c_ijk - q_sat_t_0);
                            S_s_melt = melt_driver * (0.612654e-3 + 79.6863 * inv_p_u_0) * rh_rqs_08;
                            S_g_melt = melt_driver * (7.39441e-5  + 12.31698 * inv_p_u_0) * rh_rqg_06;
                        }

                        double S_i_melt = (t_u > m.t_0 && ice_i > 0.0)
                            ? ice_i / dt_snow_dim : 0.0;

                        // --- Freezing ---
                        double S_r_frz;
                        if(t_u > t_hn){
                            double t_frz = std::max(t_r_frz - t_u, 0.0);
                            S_r_frz = c_r_frz * pow(t_frz, 1.5)
                                    * pow(r_h_i * r_q_r, 27.0/16.0);
                        }else{
                            // Below the homogeneous-nucleation temperature every raindrop
                            // freezes within the fall time. The quantity that has to be divided
                            // by that time is a MIXING RATIO, kg/kg, so on the raw-flux branch
                            // it is the rain mass content over the air density; the shipped
                            // branch's `Rain` is a dimensionless ratio and only looks like one
                            // by accident.
                            S_r_frz = raw_flux ? r_q_r / (r_h_i * dt_rain_dim)
                                               : Rain / dt_rain_dim;
                        }

                        // --- Snow-to-graupel conversion ---
                        double S_csg = (t_u <= m.t_0 && cl_i > 0.0002)
                            ? z_csg * cl_i * ((r_q_s > 0.0) ? pow(r_h_i * r_q_s, 0.75) : 0.0)
                            : 0.0;

                        // --- Availability limiters (ATM_ICE_LIMITERS) ---
                        if(limiters){
                            // Melting is a one-way process.
                            S_s_melt = std::max(0.0, S_s_melt);
                            S_g_melt = std::max(0.0, S_g_melt);

                            // Cloud water: seven sinks share what the cell holds.
                            double tot_c = S_c_au + S_ac + S_c_frz
                                         + S_s_rim + S_g_rim + S_s_shed + S_g_shed;
                            double max_c = cl_i / dt_snow_dim;
                            if(tot_c > max_c && tot_c > 0.0){
                                const double f = max_c / tot_c;
                                S_c_au *= f; S_ac *= f; S_c_frz *= f;
                                S_s_rim *= f; S_g_rim *= f; S_s_shed *= f; S_g_shed *= f;
                            }

                            // Cloud ice: six sinks share what the cell holds.
                            double tot_i = S_i_au + S_d_au + S_s_agg + S_g_agg
                                         + S_i_cri + S_i_melt;
                            double max_i = ice_i / dt_snow_dim;
                            if(tot_i > max_i && tot_i > 0.0){
                                const double f = max_i / tot_i;
                                S_i_au *= f; S_d_au *= f; S_s_agg *= f; S_g_agg *= f;
                                S_i_cri *= f; S_i_melt *= f;
                            }

                            // Falling species: a sink cannot take more than the flux carries.
                            //
                            // CHARGED AGAINST THE ARRIVING FLUX, NOT THE STALE ONE. The first
                            // version of this clip read `m.P_rain.x[i][j][k]`, which at this
                            // point in the loop is the PREVIOUS pass's value at this level --
                            // and the sink computed here is not charged against that at all. It
                            // is charged one level LOWER, as
                            //
                            //     P_x[i-1] = clamp(P_x[i] + rho[i]*S_x[i]*dz[i-1]),
                            //
                            // so the flux it can take from is `P_x[i]`, which this iteration is
                            // about to compute from `P_x[i+1]` and `S_x[i+1]` -- both already
                            // known -- and the mass it is spread over is `rho[i]*dz[i-1]`, NOT
                            // `rho[i]*dz[i]`. Measured with the stale form: the floor still
                            // injected 5930 mm/a of rain against 1699 of sources WITH the
                            // limiter on, i.e. the clip did not bind at all.
                            //
                            // The temperature window is applied here too, because a species
                            // whose window is shut at this level arrives with nothing.
                            const double dz_cons  = step_table[(i > 0) ? (i - 1) : 0];
                            const double mass_arr = r_h_i * dz_cons;
                            const double rho_up   = m.r_humid.x[i+1][j][k];
                            auto arriving = [&](double P_up, double S_up, bool window_open){
                                const double raw = P_up + rho_up * S_up * step_i;
                                return window_open
                                    ? std::min(IceSchemeCommon::P_max_flux, std::max(0.0, raw))
                                    : 0.0;
                            };
                            const double P_r = arriving(m.P_rain.x[i+1][j][k],
                                                        m.S_r.x[i+1][j][k],
                                                        (t_u >= m.t_0)) / mass_arr;
                            const double P_s = arriving(m.P_snow.x[i+1][j][k],
                                                        m.S_s.x[i+1][j][k],
                                                        (t_u < m.t_0 && t_u >= m.t_000)) / mass_arr;
                            const double P_g = arriving(m.P_graupel.x[i+1][j][k],
                                                        m.S_g.x[i+1][j][k],
                                                        (t_u < m.t_0 && t_u >= m.t_00)) / mass_arr;

                            double tot_r = S_ev + S_r_cri + S_r_frz;
                            if(tot_r > P_r && tot_r > 0.0){
                                const double f = P_r / tot_r;
                                S_ev *= f; S_r_cri *= f; S_r_frz *= f;
                            }
                            double tot_s = S_s_melt + S_csg + std::max(0.0, -S_s_dep);
                            if(tot_s > P_s && tot_s > 0.0){
                                const double f = P_s / tot_s;
                                S_s_melt *= f; S_csg *= f;
                                if(S_s_dep < 0.0) S_s_dep *= f;
                            }
                            double tot_g = S_g_melt + std::max(0.0, -S_g_dep);
                            if(tot_g > P_g && tot_g > 0.0){
                                const double f = P_g / tot_g;
                                S_g_melt *= f;
                                if(S_g_dep < 0.0) S_g_dep *= f;
                            }
                        }

                        // --- Sinks and sources ---
                        double S_c_c_ijk = m.S_c_c.x[i][j][k];

                        m.S_v.x[i][j][k] = -S_c_c_ijk + S_ev - S_i_dep - S_s_dep - S_g_dep - S_nuc;
                        m.S_c.x[i][j][k] =  S_c_c_ijk - S_c_au - S_ac - S_c_frz + S_i_melt
                                           - S_s_rim - S_g_rim - S_s_shed - S_g_shed;
                        m.S_i.x[i][j][k] =  S_nuc + S_c_frz + S_i_dep - S_i_melt
                                           - S_i_au - S_d_au - S_s_agg - S_g_agg - S_i_cri;
                        m.S_r.x[i][j][k] =  S_c_au + S_ac - S_ev + S_s_shed + S_g_shed
                                           - S_r_cri - S_r_frz + S_s_melt + S_g_melt;
                        m.S_s.x[i][j][k] =  S_i_au + S_d_au + S_s_agg + S_s_rim + S_s_dep
                                           + S_i_cri + S_r_cri - S_s_melt - S_csg;
                        m.S_g.x[i][j][k] =  S_g_agg + S_g_rim + S_g_dep + S_i_cri + S_r_cri
                                           + S_r_frz - S_g_melt + S_csg;

                        // --- Column integration (top-down) ---
                        // Each flux is clamped to [0, P_max_flux]. The upper cap is a hard backstop
                        // against the P_snow riming runaway (see the floored normalization above):
                        // even if a source spikes, the accumulated flux can no longer overflow to
                        // inf, so the scheme stays finite (over-precipitates at worst, like the
                        // other untuned schemes — usable, not NaN). ~260 mm/d, well above any
                        // physical precip.
                        constexpr double P_max_flux = 3.0e-3;        // kg/(m2*s) ~260 mm/d hard cap
                        m.P_rain.x[i][j][k] = (t_u >= m.t_0)
                            ? std::min(P_max_flux, std::max(0.0, m.P_rain.x[i+1][j][k]
                                + m.r_humid.x[i+1][j][k] * m.S_r.x[i+1][j][k] * step_i))
                            : 0.0;

                        m.P_snow.x[i][j][k] = (t_u < m.t_0 && t_u >= m.t_000)
                            ? std::min(P_max_flux, std::max(0.0, m.P_snow.x[i+1][j][k]
                                + m.r_humid.x[i+1][j][k] * m.S_s.x[i+1][j][k] * step_i))
                            : 0.0;

                        m.P_graupel.x[i][j][k] = (t_u < m.t_0 && t_u >= m.t_00)
                            ? std::min(P_max_flux, std::max(0.0, m.P_graupel.x[i+1][j][k]
                                + m.r_humid.x[i+1][j][k] * m.S_g.x[i+1][j][k] * step_i))
                            : 0.0;

                        // Track column maxima (thread-local via reduction)
                        local_max_rain    = std::max(local_max_rain,    m.P_rain.x[i][j][k]);
                        local_max_snow    = std::max(local_max_snow,    m.P_snow.x[i][j][k]);
                        local_max_graupel = std::max(local_max_graupel, m.P_graupel.x[i][j][k]);

                        m.Precipitation.x[i][j][k] = m.P_rain.x[i][j][k]
                            + m.P_snow.x[i][j][k] + m.P_graupel.x[i][j][k];

                        // ---- ATM_SS_DIAG: charge this level's increment to its terms ----
                        // The increment that produced P_x[i] was built from the terms at i+1,
                        // which the previous pass of this loop saved in pv[]. The raw value is
                        // re-formed from the model's own S_x[i+1] rather than from the sum of
                        // pv[], so the clamp bucket absorbs any floating-point reassociation and
                        // the identity closes exactly.
                        if(ssd_on && i >= i_gnd){
                            const size_t b = (size_t)j * m.km + k, N = (size_t)m.jm * m.km;
                            const double rho_up = m.r_humid.x[i+1][j][k];
                            const double raw_s = m.P_snow.x[i+1][j][k]
                                               + rho_up * m.S_s.x[i+1][j][k] * step_i;
                            const double raw_g = m.P_graupel.x[i+1][j][k]
                                               + rho_up * m.S_g.x[i+1][j][k] * step_i;
                            const double raw_r = m.P_rain.x[i+1][j][k]
                                               + rho_up * m.S_r.x[i+1][j][k] * step_i;
                            if(pv_valid){
                                for(int q = 0;      q < SS_END; q++) ssd[q*N+b] += pv_rho * pv[q] * step_i;
                                for(int q = SG_agg; q < SG_END; q++) ssd[q*N+b] += pv_rho * pv[q] * step_i;
                                for(int q = SR_c_au;q < SR_END; q++) ssd[q*N+b] += pv_rho * pv[q] * step_i;
                            }else{
                                // i = im-2: S_x at the lid was never written by this call.
                                ssd[SX_top_s*N+b] += rho_up * m.S_s.x[i+1][j][k] * step_i;
                                ssd[SX_top_g*N+b] += rho_up * m.S_g.x[i+1][j][k] * step_i;
                                ssd[SX_top_r*N+b] += rho_up * m.S_r.x[i+1][j][k] * step_i;
                            }
                            // Walked exactly as the integration applies them, in order:
                            // floor, then cap, then the temperature window. Their sum is
                            // kept - raw, which is what closes the budget.
                            auto split = [&](double raw, bool window_open,
                                             int q_floor, int q_cap, int q_win){
                                const double floored = std::max(0.0, raw);
                                const double capped  = std::min(P_max_flux, floored);
                                const double kept    = window_open ? capped : 0.0;
                                ssd[q_floor*N+b] += floored - raw;
                                ssd[q_cap  *N+b] += capped  - floored;
                                ssd[q_win  *N+b] += kept    - capped;
                            };
                            split(raw_s, (t_u < m.t_0 && t_u >= m.t_000),
                                  SX_floor_s, SX_cap_s, SX_win_s);
                            split(raw_g, (t_u < m.t_0 && t_u >= m.t_00),
                                  SX_floor_g, SX_cap_g, SX_win_g);
                            split(raw_r, (t_u >= m.t_0),
                                  SX_floor_r, SX_cap_r, SX_win_r);
                        }
                        if(ssd_on){
                            pv[SS_i_au]  = S_i_au;   pv[SS_d_au]  = S_d_au;
                            pv[SS_agg]   = S_s_agg;  pv[SS_rim]   = S_s_rim;
                            pv[SS_dep]   = S_s_dep;  pv[SS_i_cri] = S_i_cri;
                            pv[SS_r_cri] = S_r_cri;  pv[SS_melt]  = -S_s_melt;
                            pv[SS_csg]   = -S_csg;
                            pv[SG_agg]   = S_g_agg;  pv[SG_rim]   = S_g_rim;
                            pv[SG_dep]   = S_g_dep;  pv[SG_i_cri] = S_i_cri;
                            pv[SG_r_cri] = S_r_cri;  pv[SG_r_frz] = S_r_frz;
                            pv[SG_melt]  = -S_g_melt; pv[SG_csg]  = S_csg;
                            pv[SR_c_au]  = S_c_au;   pv[SR_ac]    = S_ac;
                            pv[SR_ev]    = -S_ev;    pv[SR_s_shed]= S_s_shed;
                            pv[SR_g_shed]= S_g_shed; pv[SR_r_cri] = -S_r_cri;
                            pv[SR_r_frz] = -S_r_frz; pv[SR_s_melt]= S_s_melt;
                            pv[SR_g_melt]= S_g_melt;
                            pv_rho   = m.r_humid.x[i][j][k];
                            pv_valid = true;
                        }

                    } // end i

                    if(ssd_on){
                        const size_t b = (size_t)j * m.km + k, N = (size_t)m.jm * m.km;
                        ssd[SX_gnd_s*N+b] = m.P_snow.x[i_gnd][j][k];
                        ssd[SX_gnd_g*N+b] = m.P_graupel.x[i_gnd][j][k];
                        ssd[SX_gnd_r*N+b] = m.P_rain.x[i_gnd][j][k];
                    }

                    double P_rain_diff = fabs(m.P_rain.x[23][j][k] - Rain_check) * conv_mmd;
                    if(P_rain_diff <= 1.0e-3) break;

                } // end iter_prec
            } // end k
        } // end j

        global_max_rain    = local_max_rain;
        global_max_snow    = local_max_snow;
        global_max_graupel = local_max_graupel;

        if(ssd_on) printBudget(ssd);

        // Locate global max positions (single cheap pass, outside hot loop)
        if(global_max_rain > 0.0 || global_max_snow > 0.0 || global_max_graupel > 0.0){
            #pragma omp parallel for collapse(2)
            for(int j = 1; j < m.jm - 1; j++){
                for(int k = 1; k < m.km - 1; k++){
                    for(int i = 0; i < m.im; i++){
                        if(m.P_rain.x[i][j][k] == global_max_rain && global_max_rain > 0.0){
                            #pragma omp critical(rain_loc)
                            { i_rain = i; j_rain = j; k_rain = k; rain = true; }
                        }
                        if(m.P_snow.x[i][j][k] == global_max_snow && global_max_snow > 0.0){
                            #pragma omp critical(snow_loc)
                            { i_snow = i; j_snow = j; k_snow = k; snow = true; }
                        }
                        if(m.P_graupel.x[i][j][k] == global_max_graupel && global_max_graupel > 0.0){
                            #pragma omp critical(graupel_loc)
                            { i_graupel = i; j_graupel = j; k_graupel = k; graupel = true; }
                        }
                    }
                }
            }
        }
    }


    // ==================== ATM_SS_DIAG REPORT ====================
    void printBudget(const std::vector<double>& ssd) const {
        using namespace ThreeCatIce;
        const double yr = 365.0 * 8.64e4;              // [kg/(m2 s)] = [mm/s]  ->  [mm/a]
        const size_t N = (size_t)m.jm * m.km;
        std::vector<double> a((size_t)NSSD, 0.0);
        double w_tot = 0.0;
        for(int j = 0; j < m.jm; j++){
            const double w = cos((j / (double)(m.jm - 1) - 0.5) * M_PI);
            for(int k = 0; k < m.km; k++){
                w_tot += w;
                for(int q = 0; q < NSSD; q++) a[q] += w * ssd[(size_t)q * N + (size_t)j * m.km + k];
            }
        }
        for(int q = 0; q < NSSD; q++) a[q] *= yr / w_tot;

        auto sum = [&](int lo, int hi){ double t = 0.0; for(int q = lo; q < hi; q++) t += a[q]; return t; };
        const double net_s = sum(0, SS_END), net_g = sum(SG_agg, SG_END), net_r = sum(SR_c_au, SR_END);

        printf("\n      AGCM: [SS DIAG] iter %d.  ThreeCat flux budgets, cos-lat GLOBAL mean, mm/a."
               "  Every number is that term's contribution to the ground flux.\n", m.iter_n);

        printf("      AGCM: [SS DIAG] SNOW      S_i_au %9.1f  S_d_au %9.1f  S_s_agg %9.1f"
               "  S_s_rim %9.1f  S_s_dep %9.1f\n"
               "      AGCM: [SS DIAG]           S_i_cri %8.1f  S_r_cri %8.1f  S_s_melt %8.1f"
               "  S_csg %11.1f  ->  net %9.1f\n"
               "      AGCM: [SS DIAG]           stale top %6.1f  floor +%.1f  cap %.1f"
               "  window %.1f  =  ground %9.1f   (residual %.3e)\n",
               a[SS_i_au], a[SS_d_au], a[SS_agg], a[SS_rim], a[SS_dep],
               a[SS_i_cri], a[SS_r_cri], a[SS_melt], a[SS_csg], net_s,
               a[SX_top_s], a[SX_floor_s], a[SX_cap_s], a[SX_win_s], a[SX_gnd_s],
               net_s + a[SX_top_s] + a[SX_floor_s] + a[SX_cap_s] + a[SX_win_s] - a[SX_gnd_s]);

        printf("      AGCM: [SS DIAG] GRAUPEL   S_g_agg %8.1f  S_g_rim %8.1f  S_g_dep %8.1f"
               "  S_i_cri %8.1f  S_r_cri %8.1f\n"
               "      AGCM: [SS DIAG]           S_r_frz %8.1f  S_g_melt %7.1f  S_csg %11.1f"
               "  ->  net %9.1f\n"
               "      AGCM: [SS DIAG]           stale top %6.1f  floor +%.1f  cap %.1f"
               "  window %.1f  =  ground %9.1f   (residual %.3e)\n",
               a[SG_agg], a[SG_rim], a[SG_dep], a[SG_i_cri], a[SG_r_cri],
               a[SG_r_frz], a[SG_melt], a[SG_csg], net_g,
               a[SX_top_g], a[SX_floor_g], a[SX_cap_g], a[SX_win_g], a[SX_gnd_g],
               net_g + a[SX_top_g] + a[SX_floor_g] + a[SX_cap_g] + a[SX_win_g] - a[SX_gnd_g]);

        printf("      AGCM: [SS DIAG] RAIN      S_c_au %9.1f  S_ac %11.1f  S_ev %11.1f"
               "  S_s_shed %7.1f  S_g_shed %7.1f\n"
               "      AGCM: [SS DIAG]           S_r_cri %8.1f  S_r_frz %8.1f  S_s_melt %8.1f"
               "  S_g_melt %7.1f  ->  net %9.1f\n"
               "      AGCM: [SS DIAG]           stale top %6.1f  floor +%.1f  cap %.1f"
               "  window %.1f  =  ground %9.1f   (residual %.3e)\n",
               a[SR_c_au], a[SR_ac], a[SR_ev], a[SR_s_shed], a[SR_g_shed],
               a[SR_r_cri], a[SR_r_frz], a[SR_s_melt], a[SR_g_melt], net_r,
               a[SX_top_r], a[SX_floor_r], a[SX_cap_r], a[SX_win_r], a[SX_gnd_r],
               net_r + a[SX_top_r] + a[SX_floor_r] + a[SX_cap_r] + a[SX_win_r] - a[SX_gnd_r]);

        // Shares of the GROSS sources, which is the number a coefficient change acts on: a term
        // worth 2 % of the sources cannot be the lever however wrong its constant is.
        double gross_s = 0.0;
        for(int q = 0; q < SS_END; q++) if(a[q] > 0.0) gross_s += a[q];
        if(gross_s > 0.0)
            printf("      AGCM: [SS DIAG] snow sources by share:  S_i_au %.1f %%  S_d_au %.1f %%"
                   "  S_s_agg %.1f %%  S_s_rim %.1f %%  S_s_dep %.1f %%  S_i_cri %.1f %%"
                   "  S_r_cri %.1f %%\n",
                   1e2*std::max(0.0,a[SS_i_au])/gross_s, 1e2*std::max(0.0,a[SS_d_au])/gross_s,
                   1e2*std::max(0.0,a[SS_agg])/gross_s, 1e2*std::max(0.0,a[SS_rim])/gross_s,
                   1e2*std::max(0.0,a[SS_dep])/gross_s, 1e2*std::max(0.0,a[SS_i_cri])/gross_s,
                   1e2*std::max(0.0,a[SS_r_cri])/gross_s);
    }


    // ==================== BOUNDARY CONDITIONS ====================
    void applyBoundaryConditions() {
        IceSchemeCommon::extrapolateBC(m, m.P_rain,    true);   // rain:    phi-seam periodicity averaged
        IceSchemeCommon::extrapolateBC(m, m.P_snow,    true);   // snow:    phi-seam periodicity averaged
        IceSchemeCommon::extrapolateBC(m, m.P_graupel, true);   // graupel: phi-seam periodicity averaged
    }


    // ==================== TOPOGRAPHY FILL ====================
    void applyTopography() {
        IceSchemeCommon::fillTopography(m, m.P_rain);
        IceSchemeCommon::fillTopography(m, m.P_snow);
        IceSchemeCommon::fillTopography(m, m.P_graupel);
    }


    // ==================== DIAGNOSTIC REPORT ====================
    void printReport() const {
        using namespace std;
        using namespace ThreeCatIce;

        if(!rain)
            cout << endl << "      no rain fall in ThreeCategoryIceScheme found" << endl << endl;
        else
            cout << endl << "      rain fall found  i=" << i_rain << " j=" << j_rain << " k=" << k_rain
                 << " h=" << height_table[i_rain] << "m  P=" << global_max_rain * conv_mmd << "mm/d" << endl;

        if(!snow)
            cout << "      no snow fall found" << endl;
        else
            cout << "      snow fall found  i=" << i_snow << " j=" << j_snow << " k=" << k_snow
                 << " h=" << height_table[i_snow] << "m  P=" << global_max_snow * conv_mmd << "mm/d" << endl;

        if(!graupel)
            cout << "      no graupel fall found" << endl;
        else
            cout << "      graupel fall found  i=" << i_graupel << " j=" << j_graupel << " k=" << k_graupel
                 << " h=" << height_table[i_graupel] << "m  P=" << global_max_graupel * conv_mmd << "mm/d" << endl;
    }
};
