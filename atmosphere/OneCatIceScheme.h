#pragma once

#include "cAtmosphereModel.h"
#include "CloudFraction.h"
#include "IceSchemeCommon.h"

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
    constexpr double a_if = 1.92e-6;
    constexpr double a_cf = 3.97e-5;
    constexpr double E_cf = 5.0e-3;
    constexpr double N_cf_0_surf = 2.0e5;                               // 1/m3
    constexpr double N_cf_0_top = 1.0e4;                                // 1/m3
    constexpr double tau_r = 3.3e3;                                     // s, adjusted for NASA avg 2.68 mm/d
    const double q_c_crit = IceSchemeCommon::qcCrit();                                 // [kg/kg] Kessler autoconversion threshold (~0.5 g/kg). Ported from TwoCat — without it the tau_r tuning above is defeated (ALL cloud autoconverts) -> ~25x over-precip. Applied to the rain autoconversion only, which is now the scheme's ONLY direct cloud-water conversion (see the S_s note below).
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
        double N_cf = 0.0;

        double q_sat = 0.0, E_sat = 0.0;
        double q_Ice = 0.0,  E_Ice = 0.0;

        double S_frz, S_cf_frz, S_if_frz, S_au, S_ac, S_rim, S_shed, S_ev, S_melt;

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

                        // Sub-terrain guard, ported from TwoCatIceScheme.h:272 on 2026-09-07.
                        // Cells with i < i_topography are INSIDE the mountain: their t, p,
                        // r_humid, cloud and ice are sub-terrain copies, not air, and feeding
                        // them to the rate laws produces an unphysical dP_rain that the
                        // downward flux integration then carries up into the real column. In
                        // TwoCat that is what drove the iteration-323 P_rain runaway at
                        // (i=6, j=30, k=209), the Gulf of Alaska / Cook Inlet inside-mountain
                        // cell. THIS SCHEME HAD NO SUCH GUARD AT ALL -- `i_topography` did not
                        // appear in the file -- so it computed rates at every level including
                        // rock, and its flux integration accumulated from the bottom of the
                        // grid rather than from the ground.
                        //
                        // Zeroing S_i as well matters here specifically: until today S_i was
                        // never assigned in this scheme, so a stale sub-terrain value could not
                        // arise. It is assigned now (the ice->snow debit), so it needs clearing
                        // inside the terrain like every other rate.
                        if (i < m.i_topography[j][k]) {
                            m.P_rain.x[i][j][k]        = 0.0;
                            m.P_snow.x[i][j][k]        = 0.0;
                            m.Precipitation.x[i][j][k] = 0.0;
                            m.S_v.x[i][j][k] = 0.0;
                            m.S_c.x[i][j][k] = 0.0;
                            m.S_i.x[i][j][k] = 0.0;
                            m.S_r.x[i][j][k] = 0.0;
                            m.S_s.x[i][j][k] = 0.0;
                            continue;
                        }

                        double Rain = m.P_rain.x[i][j][k];
                        double Snow = m.P_snow.x[i][j][k];

                        double t_u = m.t.x[i][j][k] * m.t_0;

                        step[i] = m.get_layer_height(i+1) - m.get_layer_height(i); // local atmospheric shell thickness in m

                        E_sat = m.hp * AtomUtils::exp_func(t_u, 17.2694, 35.86); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                        q_sat = m.ep * E_sat/(m.p_stat.x[i][j][k] - E_sat);   // relativ water vapour contents on ocean surface reduced by factor in kg/kg

                        E_Ice = m.hp * AtomUtils::exp_func(t_u, 21.8746, 7.66);
                        q_Ice = m.ep * E_Ice/(m.p_stat.x[i][j][k] - E_Ice);      // relativ water vapour contents on ocean surface reduced by factor in kg/kg

                        // Snow is grown through the cloud-ICE reservoir via the shared, bounded
                        // vapour->ice->snow throttle (TwoCat's), replacing OneCat's unbounded
                        // direct vapour->snow deposition. thr.S_i_au + thr.S_d_au feed S_s below.
                        const double dt_snow_dim = step[i] / 0.96;               // snow fall time step (v_s = 0.96 m/s)
                        IceSchemeCommon::IceSnowRates thr =
                            IceSchemeCommon::depositionThrottle(m, i, j, k, t_u, q_Ice, dt_snow_dim);


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

                        // ATM_CLOUD_FRAC: every q_c_crit test below is an IN-CLOUD one, so it
                        // is applied to cloud/f and the rate scaled by f. See TwoCatIceScheme for
                        // the argument and the measurement. f = 1 off-branch, so all three are
                        // the shipped expressions exactly. S_ac is LINEAR in cloud water and
                        // therefore invariant under the transform -- untouched.
                        const double f_cld = CloudFraction::effectiveFraction(
                                max(0.0, m.c.x[i][j][k]) + max(0.0, m.cloud.x[i][j][k])
                                    + max(0.0, m.ice.x[i][j][k]),
                                q_sat, m.p_stat.x[i][j][k],
                                max(0.0, m.cloud.x[i][j][k]) + max(0.0, m.ice.x[i][j][k]));
                        const double cloud_in = m.cloud.x[i][j][k] / f_cld;   // in-cloud water

                        if(m.cloud.x[i][j][k] > 0.0)
                            S_au  = f_cld * (1.0 - eps_t)/tau_r * max(0.0, cloud_in - q_c_crit); // Kessler threshold: only cloud excess rains


                        // collection mechanisms: accretion, riming, shedding
                        // accretion of cloud water by raindrops
                        S_ac = (1.0 - eps_t) * c_ac * m.cloud.x[i][j][k] // c_ac = 0.24, in m²/kg
                            * pow(Rain,(7.0/9.0));                      // in kg/s

                        // Riming S_rim = c_rim*cloud*Snow is ∝ Snow, a positive feedback that
                        // ran P_snow away to 1e7+ once the accumulation band was widened (same
                        // instability TwoCat cured / ThreeCat had). Two guards, matching TwoCat:
                        // (1) q_c_crit threshold — in the cold snow-forming layers cloud water is
                        //     ~0.09 g/kg < q_c_crit, so riming (the runaway term) switches OFF there
                        //     while the ice->snow throttle still makes snow; (2) c_rim reduced 5x.
                        //     Snow now forms (from thr.S_i_au + thr.S_d_au) and stays bounded.
                        // Fraction-aware for the same reason as S_au. Snow is left as the GRID
                        // MEAN, exactly as Rain is in TwoCat's accretion: scaling it too needs an
                        // assumption about precipitation fraction versus cloud fraction, which
                        // this closure does not imply. Note these two are thresholded and so
                        // NONLINEAR -- an unthresholded c_rim*cloud*Snow would need no change.
                        if(t_u < m.t_0)
                            S_rim = (cloud_in > q_c_crit)
                                ? f_cld * (c_rim/5.0)/a_m * (cloud_in - q_c_crit) * Snow : 0.0;
                        else  S_rim = 0.0;                              // riming rate of snow mass due to collection of supercooled cloud droplets, < VIII >

                        if(t_u >= m.t_0)
                            S_shed = (cloud_in > q_c_crit)
                                ? f_cld * (c_rim/5.0)/a_m * (cloud_in - q_c_crit) * Snow : 0.0;
                        else  S_shed = 0.0;                             // rate of water shed by melting wet snow particles, < IX >


                        // diffusional growth of rain and snow
                        // evaporation of rain water
                        a_ev_local = 2.76e-3 * exp(0.055 * (m.t_0 - t_u));
                        S_ev = a_ev_local * (1.0 + b_ev * pow(Rain,(1.0/6.0))) // evaporation of rain due to water vapour diffusion, < XIII >
                            * (q_sat - m.c.x[i][j][k])
                            * pow(Rain, (4.0/9.0));

                        // (Old direct vapour->snow deposition S_dep removed: it was ∝ Snow^0.58,
                        //  an unbounded feedback the model's excess vapour ran to the flux cap.
                        //  Snow deposition now goes through the shared ice reservoir — thr above.)

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

                        if((t_u <= m.t_0)&&(t_u >= m.t_00)){        // FIX: was (t_u = m.t_00) — an assignment that corrupted t_u to -37C for sub-freezing cells
                            S_if_frz = a_if * (exp(a_if * (t_u - m.t_0)) - 1.0) // immersion freezing
                                * pow(Rain, (14.0/9.0));
                            S_cf_frz = a_cf * E_cf * N_cf                  // contact freezing nucleation
                                * pow(Rain, (13.0/9.0));
                            S_frz = S_if_frz + S_cf_frz;                   // mixing ratio of snow
                        }
                        else  S_frz = 0.0;


                        // ICE CONSUMPTION LIMITER, ported from TwoCatIceScheme.h:558. It is part
                        // of the S_i repair below rather than an extra: writing the ice sink
                        // without it would turn "snow from nothing" into "negative ice", and the
                        // RK4 stages clamp `ice` with std::max(0.0, ...), which is itself a water
                        // SOURCE -- the shape ATM_CWB_DIAG charges to the RungeKutta bucket. The
                        // two rates need it for different reasons: S_i_au = c_i_au*ice is
                        // proportional and safe on its own, but S_d_au is sized by the DEPOSITION
                        // rate (S_i_dep/1.5*((m_s_0/m_i)^(2/3)-1)), which is supersaturation-
                        // limited and knows nothing about how much ice is present.
                        {
                            const double S_ice_total = thr.S_i_au + thr.S_d_au;
                            const double max_ice_loss = m.ice.x[i][j][k] / dt_snow_dim;
                            if(S_ice_total > max_ice_loss && S_ice_total > 0.0){
                                const double factor = max_ice_loss / S_ice_total;
                                thr.S_i_au *= factor;
                                thr.S_d_au *= factor;
                            }
                        }

                        // sinks and sources. Snow now grows from the cloud-ICE reservoir
                        // (thr.S_i_au aggregation + thr.S_d_au depositional autoconversion),
                        // not directly from vapour/cloud — the shared, bounded throttle.
                        // Vapour->ice is owned by SaturationAdjustment upstream, so S_i_dep is
                        // used only to size S_d_au and is not re-applied to the vapour budget.
                        //
                        // S_i IS THE DEBIT FOR THE TWO ICE->SNOW TERMS, AND UNTIL 2026-09-07 IT
                        // WAS NEVER ASSIGNED AT ALL -- so `S_s` drew `thr.S_i_au + thr.S_d_au`
                        // from a reservoir nothing depleted and the scheme MANUFACTURED SNOW.
                        // Every other term here is a conserving pair (S_c_c between S_v and S_c;
                        // S_ev between S_v and S_r; S_au, S_ac, S_shed between S_c and S_r;
                        // S_rim between S_c and S_s; S_frz and S_melt between S_r and S_s), so
                        // these two were the whole of the non-closure. With the line below,
                        // S_v + S_c + S_i + S_r + S_s = 0 identically, term by term.
                        //
                        // WHY THE DEBIT IS `ice` AND NOT VAPOUR, WHICH IS THE ONE JUDGEMENT IN
                        // IT. S_d_au is sized by S_i_dep, a vapour->ice deposition rate, so it
                        // looks like a vapour term. It is not: this scheme deliberately does NOT
                        // apply S_i_dep to the budget, because SaturationAdjustment owns
                        // vapour->ice and writes `ice` DIRECTLY (SaturationAdjustment.h:252/:261)
                        // rather than through S_i. The deposited mass therefore does land in
                        // `ice` first, and S_d_au moves it on to snow. Debiting `ice` is the
                        // reading consistent with the comment above; debiting `c` would take the
                        // same water twice, once here and once in the adjustment.
                        //
                        // NOTE the sub-terrain guard TwoCatIceScheme.h:273 has and this scheme
                        // does not: OneCat computes its rates at every i including cells inside
                        // mountains. That is pre-existing and applies equally to S_v/S_c/S_r/S_s,
                        // so it is left alone here rather than half-fixed under an S_i change.
                        m.S_v.x[i][j][k] = - m.S_c_c.x[i][j][k] + S_ev;         // in kg/(kg*s)
                        m.S_c.x[i][j][k] =   m.S_c_c.x[i][j][k] - S_au - S_ac
                                           - S_rim - S_shed;
                        m.S_r.x[i][j][k] =   S_au + S_ac - S_ev + S_shed
                                           - S_frz + S_melt;
                        m.S_i.x[i][j][k] = -(thr.S_i_au + thr.S_d_au);       // the debit for the two ice->snow terms
                        m.S_s.x[i][j][k] =   thr.S_i_au + thr.S_d_au + S_rim
                                           + S_frz - S_melt;
                        // THERE IS DELIBERATELY NO DIRECT cloud-water -> snow TERM HERE, AND IT
                        // MUST NOT BE ADDED BACK. Until 633e9c6 this scheme carried
                        // S_nuc = eps_t/tau_s*cloud as a conserving pair -- credited to S_s,
                        // debited from S_c -- and that commit removed BOTH sides together when it
                        // routed snow through the bounded ice reservoir instead. The physical
                        // chain S_nuc stood for is still here, in two steps: SaturationAdjustment
                        // freezes cloud water into `ice` (:252 partitions by the ice fraction,
                        // :261 freezes it outright below the homogeneous-freezing floor), and the
                        // throttle turns that ice into snow via S_i_au + S_d_au. Re-adding S_nuc
                        // would run the SAME conversion a second time in parallel, and re-adding
                        // it to S_s alone -- the naive reading of a compiler "set but not used" --
                        // would create snow from nothing, because its matching -S_nuc in S_c is
                        // gone. The leftover assignment and the tuned tau_s were deleted
                        // 2026-09-06; the constant had been fitted when S_nuc WAS the seed and
                        // governed a term that no longer ran.
                        // (The note that used to stand here said cloud ice is DIAGNOSTIC --
                        // "the scheme's S_i tendency is not integrated into it, so ice cannot be
                        // depleted and m_i stays saturated". That was an accurate description of
                        // the defect repaired above, not of a design choice: S_i now carries the
                        // ice->snow debit and reaches `ice` through rhs_ice like every other
                        // species. Whether OneCat then reaches a realistic snow fraction is a
                        // separate question and is still open -- see project_ice_scheme_states.)


                        // rain integration
                        constexpr double P_max_flux = 3.0e-3;         // kg/(m2*s) ~260 mm/d hard cap (finiteness backstop)
                        if(t_u >= m.t_0)
                              m.P_rain.x[i][j][k] = m.P_rain.x[i+1][j][k]
                                  + m.r_humid.x[i+1][j][k] * m.S_r.x[i+1][j][k]
                                  * step[i];                            // in kg/(m² * s) == mm/s
                        else  m.P_rain.x[i][j][k] = 0.0;

                        if(m.P_rain.x[i][j][k] < 0.0)  m.P_rain.x[i][j][k] = 0.0;
                        if(m.P_rain.x[i][j][k] > P_max_flux) m.P_rain.x[i][j][k] = P_max_flux;

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


                        // snow integration. Accumulate over the FULL sub-freezing column
                        // [t_00, t_0) = [-37,0)C, not just [t_000, t_0) = [-20,0)C: the old lower
                        // bound t_000 discarded the -37..-20C layers where the ice fraction eps_t
                        // (and hence the ice reservoir the snow grows from) is STRONGEST, so the
                        // snow that formed was zeroed before it could fall — OneCat produced zero
                        // snow everywhere.
                        if((t_u < m.t_0)&&(t_u >= m.t_00))
                              m.P_snow.x[i][j][k] = m.P_snow.x[i+1][j][k]
                                  + m.r_humid.x[i+1][j][k] * m.S_s.x[i+1][j][k]
                                  * step[i];                            // in kg/(m² * s) == mm/s
                        else  m.P_snow.x[i][j][k] = 0.0;

                        if(m.P_snow.x[i][j][k] < 0.0)  m.P_snow.x[i][j][k] = 0.0;
                        if(m.P_snow.x[i][j][k] > P_max_flux) m.P_snow.x[i][j][k] = P_max_flux;

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
        IceSchemeCommon::extrapolateBC(m, m.P_rain, true);   // rain: phi-seam periodicity averaged
        IceSchemeCommon::extrapolateBC(m, m.P_snow, true);   // snow: phi-seam periodicity averaged
    }


    // ==================== TOPOGRAPHY FILL ====================
    void applyTopography() {
        IceSchemeCommon::fillTopography(m, m.P_rain);
        IceSchemeCommon::fillTopography(m, m.P_snow);
        IceSchemeCommon::fillTopography(m, m.S_r);
        IceSchemeCommon::fillTopography(m, m.S_s);
        IceSchemeCommon::fillTopography(m, m.S_v);
        IceSchemeCommon::fillTopography(m, m.S_c);
        IceSchemeCommon::fillTopography(m, m.S_c_c);
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
