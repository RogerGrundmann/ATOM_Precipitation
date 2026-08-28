#pragma once

#include "cAtmosphereModel.h"
#include "Utils.h"

#include <vector>
#include <cmath>
#include <chrono>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstdlib>
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

        // Surface albedo with a temperature-dependent ICE/SNOW FEEDBACK (replaces the former
        // fixed pole->equator parabola, which froze the feedback out). Base value by surface
        // type (ocean vs land); where the surface is cold, blend toward a bright snow/sea-ice
        // albedo over a smooth ramp around freezing. Recomputed every MLR call on the CURRENT
        // surface T, so it is a LIVE feedback: warming -> less ice -> lower albedo -> more
        // absorbed shortwave -> extra (polar-amplified) warming that the dynamics cannot mix
        // away (unlike a forcing/nudge). Constants are tunable. NOTE the cloud SW bump below
        // overwrites this per column where cloud is present, so its net reach is cloud-limited.
        constexpr double alb_ocean  = 0.08;    // open water
        constexpr double alb_land   = 0.20;    // generic snow-free land
        constexpr double alb_ice    = 0.60;    // snow / sea-ice
        constexpr double T_ice_none = 275.15;  // surface T >= this -> ice-free   [+2 C]
        constexpr double T_ice_full = 265.15;  // surface T <= this -> full ice   [-8 C]
        #pragma omp parallel for schedule(static)
        for (int j = 0; j < m.jm; j++)
            for (int k = 0; k < m.km; k++)
                for (int i = 0; i < m.im - 1; i++) {
                    const bool ocean = is_ocean_surface(m.h, i, j, k);
                    if (ocean || is_land_surface(m.h, i, j, k)) {
                        const double a_base = ocean ? alb_ocean : alb_land;
                        const double T_s    = m.t.x[i][j][k] * m.t_0;          // surface T [K]
                        double f_ice = (T_ice_none - T_s) / (T_ice_none - T_ice_full);
                        f_ice = std::max(0.0, std::min(1.0, f_ice));           // 0 (warm) -> 1 (cold)
                        m.albedo.y[j][k] = a_base + (alb_ice - a_base) * f_ice;
                        break;   // first surface cell per column
                    }
                }

        // Incoming short-wave radiation: pole -> equator parabola, hemispherically symmetric.
        m.short_wave_radiation = std::vector<double>(m.jm, m.rad_pole_short);
        const double rad_short_eff = m.rad_pole_short - m.rad_equator_short;
        for (int j = j_half; j >= 0; j--)
            m.short_wave_radiation[j] =
                rad_short_eff * parabola((double)j / (double)j_half) + m.rad_pole_short;
        for (int j = j_max; j > j_half; j--)
            m.short_wave_radiation[j] = m.short_wave_radiation[j_max - j];

        // ATM_RAD_TOPO -- see the i_mount comment inside the column loop below.
        //
        // FLIPPED ON 2026-08-28 AND REVERTED THE SAME DAY, BY ITS OWN VERIFICATION ARM. The
        // off-branch case for flipping is real and stands: at 100 iterations max tau_layer is
        // 7.885 over the Himalaya against 0.022 over ocean, and the original 4-iteration
        // comparison (0.088 -> 0.119) was taken before that defect even appears -- it emerges as
        // a STEP between iterations 20 and 40. But the on-branch arm showed the cure is worse:
        // 622 cells reach layer emissivity > 0.9 and 270 exceed 0.99 (max 0.99973 at 28N 88E),
        // against ZERO such cells off-branch, and land-mean radiation falls 327.1 -> 306.4 W/m2.
        // That is the saturation pathology the de-saturation split above exists to prevent.
        // See the tau_dry mass scaling below for why, and what has to be true before it flips.
        static const bool topo_rad = [](){
            const char* e = getenv("ATM_RAD_TOPO"); return e && atoi(e) != 0; }();
        // ATM_SFC_COUPLED -- see the surface/column consistency block in the column loop.
        static const bool sfc_coupled = [](){
            const char* e = getenv("ATM_SFC_COUPLED"); return e && atoi(e) != 0; }();

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

            for (int k = 0; k < m.km; k++) {

                // THE RADIATION COLUMN STARTS AT THE GROUND, NOT AT LEVEL 0.
                //
                // i_mount was the literal 0 with the comment "surface / bottom layer", which is
                // true only over ocean. Over topography levels 0 .. i_topography-1 are rock, and
                // ThermoAtm's barometric loop fills every one of them with a real p_stat, so they
                // carry MASS: dp = p_stat[i] - p_stat[i+1] is positive through the rock and enters
                // sum_dp. Since each layer's optical depth is tau_dry * dp_i / sum_dp, every AIR
                // layer over every mountain is diluted by the rock beneath it -- so tau_layer and
                // tau_above are wrong ABOVE the ground too, not merely inside it. The surface
                // energy balance, e_surf, epsilon_2D and T_air1 are all placed at sea level by the
                // same literal.
                //
                // This is the family's Earth-constant pattern with terrain in place of a number:
                // the code assumes the surface is level 0, which holds at home only over water.
                //
                // ATM_RAD_TOPO=1 puts the column on the real ground. DEFAULT OFF because it moves
                // the radiation over all land and this tree measures before it flips. Off-branch
                // is bit-identical: i_mount is 0 and every loop below reduces to what it was.
                const int i_mount = topo_rad
                    ? std::min(std::max(m.i_topography[j][k], 0), i_trop - 3)
                    : 0;

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
                // tau_dry IS A SEA-LEVEL WHOLE-COLUMN VALUE, AND A MOUNTAIN COLUMN IS NOT ONE.
                // eps_dry = 0.684 is Bignami's clear-sky emissivity for a full column; the
                // optical depth it implies, 1.152, was handed unchanged to every column and then
                // split by dp_i/Sum(dp). With i_mount on the ground Sum(dp) telescopes to
                // p_stat[i_mount], the GROUND pressure, so a 5 km column carries ~55 % of the air
                // mass and was still given 100 % of the dry optical depth. Scaling by the actual
                // mass, Sum(dp)/p_0, makes the Sum(dp) cancel -- each layer's dry opacity becomes
                // tau_dry_full*dp_i/p_0, proportional to its OWN mass and independent of where the
                // column starts. That is the right form and it is kept for that reason.
                //
                // IT IS NOT WHAT CAUSED THE SATURATION, AND THIS COMMENT ONCE SAID IT WAS.
                // Measured: with the scaling in, cells at layer emissivity > 0.9 go 622 -> 639 and
                // > 0.99 goes 270 -> 271, land radiation 306.4 -> 305.9 W/m2. No improvement.
                // The driver is tau_cloud: the saturated cells carry mean CloudWater 0.4016
                // against 0.0328 over land as a whole, TWELVE times, while their water vapour is
                // slightly LOWER than average. cwp_cap_col caps the COLUMN condensate path, and
                // with the rock excluded that same 250 g/m2 is shared among fewer, thinner air
                // layers -- same column total, larger per-layer LWP_i, and eps_i = 1-exp(-tau_i)
                // saturates. Capping a column does not bound a layer. That is what has to be
                // fixed before ATM_RAD_TOPO can be the default.
                //
                // Gated on topo_rad so the sea-level branch stays bit-identical: there Sum(dp) is
                // the surface pressure and the ratio would be ~1.01 rather than exactly 1, which
                // would silently move every recorded off-branch number.
                const double tau_dry_full = -log(1.0 - eps_dry);     // Bignami, sea-level column
                const double mass_frac    = (m.p_0 > 0.0) ? (sum_dp / m.p_0) : 1.0;
                const double tau_dry = topo_rad ? tau_dry_full * mass_frac : tau_dry_full;
                const double tau_col = -log(1.0 - eps_col);
                const double tau_wv  = (tau_col > tau_dry) ? (tau_col - tau_dry) : 0.0;
                const double inv_dp  = (sum_dp > 0.0) ? 1.0 / sum_dp : 0.0;
                const double inv_vp  = (sum_vp > 0.0) ? 1.0 / sum_vp : 0.0;

                // Physical cap on the column cloud condensate the RADIATION sees. The model
                // over-condenses (column LWP ~1500 g/m2 vs observed ~100), which saturates BOTH
                // the LW cloud greenhouse (each layer's k_liq*LWP_i ~5 -> ~opaque) AND the SW
                // albedo bump (pinned at the cloud value everywhere). Compute the raw column path
                // and a single uniform scale so the radiation treats the column as a physically
                // thick cloud (<= cwp_cap_col), preserving the vertical cloud DISTRIBUTION.
                // Applied to LWP_i/IWP_i below, it feeds BOTH the LW tau_cloud and the SW path, so
                // the two stay BALANCED — fixing only one (e.g. SW albedo alone) removes the
                // excess cooling but leaves the excess greenhouse and tips the climate hot.
                constexpr double cwp_cap_col = 250.0;                // g/m2 physical thick-cloud column condensate
                double cwp_raw = 0.0;
                for (int i = i_mount; i <= i_trop; i++) {
                    const double dz_i   = (i < i_trop) ? (m.get_layer_height(i+1) - m.get_layer_height(i))
                                                       : (m.get_layer_height(i) - m.get_layer_height(i-1));
                    const double T_ii   = m.t.x[i][j][k] * m.t_0;
                    const double rho_ii = (T_ii > 0.0) ? (m.p_stat.x[i][j][k] * 100.0) / (287.0 * T_ii) : 0.0;
                    const double cwl    = (m.cloud.x[i][j][k] > 0.0) ? m.cloud.x[i][j][k] : 0.0;
                    const double cwi    = (m.ice.x[i][j][k]   > 0.0) ? m.ice.x[i][j][k]   : 0.0;
                    cwp_raw += (cwl + cwi) * rho_ii * dz_i * 1000.0;
                }
                const double cloud_scale = (cwp_raw > cwp_cap_col) ? (cwp_cap_col / cwp_raw) : 1.0;

                double lwp_col = 0.0, iwp_col = 0.0;                  // accumulated (scaled) condensate paths [g/m2] (for the SW albedo bump)
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
                    // co2.x is stored in ppm (ThermoAtm init, BCs, transport all use ppm): a
                    // mixing ratio co2*1e-6 times the local pressure fraction p_stat/p_0 gives the
                    // CO2 partial pressure in atm. (Earlier this multiplied by co2_0, which assumed
                    // co2.x was a ~1 ratio — wrong for the ppm field; it over-saturated the band
                    // regardless of co2_0. See project_multilayer_radiation CO2-units fix.)
                    const double P_c = 1e-6 * m.p_stat.x[i][j][k] * m.co2.x[i][j][k] / m.p_0; // [atm], co2 in ppm
                    double u_c = P_c * dz * 100.0;                                                       // [atm*cm]
                    if (u_c < 0.0) u_c = 0.0;
                    double eps_co2 = co2_band_scale * 0.185 * (1.0 - exp(-0.3919 * pow(u_c, 0.4)));       // Atwater & Ball
                    if (eps_co2 > 0.999) eps_co2 = 0.999;
                    const double tau_co2 = -log(1.0 - eps_co2);

                    // Cloud liquid + ice LONGWAVE greenhouse (suspended condensate only; the
                    // precipitation fluxes P_rain/P_snow are NOT radiatively active here). Layer
                    // water/ice paths [g/m2] = mixing ratio [kg/kg] * air density [kg/m3] * layer
                    // thickness dz [m] * 1000; times a mass-absorption coefficient [m2/g] gives a
                    // dimensionless optical depth that composes additively with the dry/vapour/CO2
                    // depths. k_liq/k_ice after Stephens (1978): liquid ~0.10-0.15, ice ~0.05-0.06
                    // (ice less absorbing). Density computed locally from the ideal-gas law (p_stat
                    // in hPa -> Pa) so this needs no populated r_humid, matching the CO2 band's use
                    // of p_stat. LW ONLY: clouds here can only ADD greenhouse (raise L_down / warm
                    // the surface). The compensating SHORTWAVE (cloud-albedo cooling) is applied
                    // below the column loop as an albedo bump on the accumulated condensate path
                    // (lwp_col/iwp_col), so low thick cloud can NET-cool while thin cirrus stays
                    // net-warming. See project_multilayer_radiation cloud/ice plan.
                    constexpr double k_liq = 0.12, k_ice = 0.055;                 // LW mass absorption [m2/g]
                    const double T_i   = m.t.x[i][j][k] * m.t_0;                  // [K]
                    const double rho_i = (T_i > 0.0) ? (m.p_stat.x[i][j][k] * 100.0) / (287.0 * T_i) : 0.0; // [kg/m3]
                    const double cw_l  = (m.cloud.x[i][j][k] > 0.0) ? m.cloud.x[i][j][k] : 0.0; // [kg/kg]
                    const double cw_i  = (m.ice.x[i][j][k]   > 0.0) ? m.ice.x[i][j][k]   : 0.0; // [kg/kg]
                    double LWP_i = cloud_scale * cw_l * rho_i * dz * 1000.0;  // liquid water path [g/m2], capped
                    double IWP_i = cloud_scale * cw_i * rho_i * dz * 1000.0;  // ice   water path [g/m2], capped
                    double tau_cloud = k_liq * LWP_i + k_ice * IWP_i;

                    // PER-LAYER CLOUD OPTICAL-DEPTH CEILING. cwp_cap_col above bounds the COLUMN
                    // condensate path; it does not bound a LAYER, and the two are not the same
                    // constraint. With ATM_RAD_TOPO on, the same capped 250 g/m2 is shared among
                    // fewer, thinner air layers once the rock is excluded, so a single layer goes
                    // optically black: 622 cells at eps > 0.9 and 270 above 0.99 (max 0.99973),
                    // against ZERO off-branch, and land radiation down 20.7 W/m2. The saturated
                    // cells carry twelve times the mean land CloudWater while their water vapour
                    // is slightly BELOW average, so it is the condensate and not the vapour.
                    // A near-blackbody layer is exactly what the de-saturation split above exists
                    // to prevent -- it pins the emission level and collapses the OLR.
                    //
                    // Scale LWP_i AND IWP_i by the same factor rather than clipping tau_cloud, so
                    // the LW (tau_cloud) and the SW (lwp_col/iwp_col, the albedo bump below) see
                    // the same condensate. Capping one and not the other removes the excess
                    // cooling but leaves the excess greenhouse -- the imbalance the cwp_cap_col
                    // note above already warns about.
                    //
                    // GATED, AND THE REASON IS A PREDICTION OF MINE THAT WAS WRONG. This was first
                    // written ungated, argued "inert on the shipped branch" from the off-branch max
                    // layer emissivity of 0.67683. That number is the maximum AT LEVEL 0 -- the
                    // radial slice -- and it is not the column maximum. Measured over all 41 levels
                    // on one latitude slice, the SHIPPED branch reaches eps = 0.99962 at level 28
                    // with 1773 cells above 0.9. The saturation this ceiling exists to stop is
                    // ALREADY PRESENT ALOFT with ATM_RAD_TOPO off; it is not something the topo
                    // branch introduces, only something the topo branch brings down to the ground
                    // where a surface diagnostic finally showed it.
                    //
                    // So the ceiling changes shipped behaviour -- land radiation 327.08 -> 326.82,
                    // ocean 362.26 -> 362.19, and 1773 -> 0 saturated cells per slice -- and goes
                    // DEFAULT 2.0 (ON) since 2026-08-28; ATM_CLOUD_TAU_MAX=0 disables it exactly.
                    // Flipped on the top-of-atmosphere evidence, which is where a radiation change
                    // has to be judged: clear-sky OLR is BIT-IDENTICAL across the flip
                    // (180.33882 W/m2 both ways) and cloudy OLR moves +0.053 (177.607 -> 177.660),
                    // while 1773 near-blackbody layers per latitude slice become zero. Removing a
                    // pathology the design already forbids, at 0.05 W/m2, is worth taking.
                    static const double tau_cloud_max = [](){
                        const char* e = getenv("ATM_CLOUD_TAU_MAX");
                        const double v = e ? atof(e) : 2.0;
                        return v > 0.0 ? v : 0.0; }();
                    if (tau_cloud_max > 0.0 && tau_cloud > tau_cloud_max) {
                        const double f = tau_cloud_max / tau_cloud;
                        LWP_i *= f;  IWP_i *= f;  tau_cloud = tau_cloud_max;
                    }
                    lwp_col += LWP_i;  iwp_col += IWP_i;                          // column paths for the SW albedo bump

                    double tau = tau_dry * dp_col[i] * inv_dp + tau_wv * vpath_col[i] * inv_vp
                               + tau_co2 + tau_cloud;
                    m.epsilon.x[i][j][k] = 1.0 - exp(-tau);

                    // Stash the LAYER optical depth; the downward pass below turns tau_above
                    // into the cumulative-from-the-lid value. Not recoverable from epsilon
                    // afterwards, which saturates at tau ~ 37.
                    m.tau_above.x[i][j][k] = tau;
                    m.tau_layer.x[i][j][k] = tau;
                }

                // tau_above: walk DOWNWARD from the lid, accumulating, so tau_above[i] is the
                // optical depth of everything ABOVE level i -- 0 at the lid, increasing
                // downward, and tau_above = 1 is the photosphere. One thread owns this whole
                // column (the parallel for is over j, k inner), so this is race-free.
                {
                    double acc = 0.0;
                    for (int i = i_trop; i >= i_mount; i--) {
                        const double layer = m.tau_above.x[i][j][k];
                        m.tau_above.x[i][j][k] = acc;
                        acc += layer;
                    }
                }
                m.epsilon_2D.y[j][k] = m.epsilon.x[i_mount][j][k];

                // Cloud/ice SHORTWAVE albedo bump (stage 2): reflective clouds raise the column
                // albedo toward a cloud value, cutting absorbed SW (cooling) to compete with the LW
                // greenhouse above. alpha_eff = alpha_surf + (alpha_cloud - alpha_surf)*
                // (1 - exp(-k_sw*CWP_sw)), CWP_sw = LWP + f_ice*IWP (ice is optically thinner / less
                // reflective per unit path, so weighted down -> thin cirrus stays net-warming while
                // thick low liquid cloud net-cools). Overwrites the clear-sky latitude albedo for
                // this column; it is read below by the SW terms (tridiagonal source dd + surface
                // energy balance) and re-derived fresh each MLR call. See project_multilayer_radiation.
                {
                    // Two-stream-like cloud SW reflectivity on a CAPPED condensate path. The old
                    // exp(-k_sw*CWP) with k_sw=0.030 saturated at CWP~100 g/m2, so — and especially
                    // with the model's excessive cloud water (LWP ~1500 g/m2 vs observed ~100) —
                    // EVERY cloudy column was pinned at the asymptotic 0.60, wiping out the latitude
                    // gradient and MASKING the surface ice-albedo feedback (poles read the same 0.60
                    // as tropical ocean). Instead cap the path the radiation sees at a physical
                    // thick-cloud value and let reflectivity rise GENTLY as refl = tau/(tau+2)
                    // (0.5 at tau=2), composited over the surface (ice-feedback) albedo. Cloudy
                    // tropics now land ~0.3, thin/clear cells relax toward the surface value, and
                    // the polar ice albedo shows through — a physical gradient. (The same cloud-water
                    // excess still inflates the LW tau_cloud above; capping that is the next step.)
                    constexpr double alpha_cloud = 0.50;   // thick cloud-top SW albedo
                    constexpr double f_ice_sw    = 0.50;   // ice SW reflectivity weight vs liquid
                    constexpr double cwp_tau     = 100.0;  // g/m2 per unit effective optical thickness
                    const double cwp_sw = lwp_col + f_ice_sw * iwp_col;          // already capped via cloud_scale above
                    const double tau    = cwp_sw / cwp_tau;
                    const double refl   = tau / (tau + 2.0);                     // gentle saturation (0.5 at tau=2)
                    const double a0     = m.albedo.y[j][k];                      // surface (ice-feedback) albedo
                    if (alpha_cloud > a0)
                        m.albedo.y[j][k] = a0 + (alpha_cloud - a0) * refl;
                }

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
                    for (int l = i_mount + 1; l <= i - 1; l++) {
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
                for (int i = i_trop - 1; i >= i_mount; i--)
                    m.radiation.x[i][j][k] = alfa[i] * m.radiation.x[i + 1][j][k] + beta[i];

                // Radiation -> temperature (add back the reference emission, invert sigma T^4).
                for (int i = i_mount; i <= i_trop; i++) {
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

                // NOTE (2026-08-27): ATM_SFC_COUPLED lived here and is REMOVED as a measured
                // null written on a wrong premise. It solved the surface balance together with a
                // neutral-lapse constraint, to cure brunt_N2 < 0 in the bottom layer. It never
                // fired: at setup the column is STABLE (-3.93 K/km) and the two arms were
                // identical at iteration 0, because the instability is manufactured by the TIME
                // LOOP -- and MLR's seven calls are all at lines 713-889 while that loop starts
                // at 985, so MLR cannot maintain it. The real gap was that this tree had no dry
                // convective adjustment at all; see ConvectiveAdjustment.h.
                //
                // What survives, and is NOT fixed here: the balance above debits the surface
                // c_H*(T_s - T_air1) and nothing credits it to the air. `Q_Sensible` is written
                // in RHS_Atm_Turb.cpp:485 and read by NOTHING. A real conservation defect in the
                // initialisation, recorded rather than silently repaired.

                // De-kink the surface radiative step in the DIAGNOSTIC radiation profile only.
                // The 1-point surface energy balance (sigma T_s^4 at i_mount) and the column
                // Thomas solve above it are computed separately, leaving a sharp discontinuity
                // at the surface. One light 1-2-1 pass over the lowest layers softens it for
                // plotting. This touches ONLY radiation.x — t.x / the surface-balance T_s (and
                // hence the CO2 surface sensitivity) are left exactly as computed above.
                {
                    const int i_top_sm = std::min(i_mount + 4, i_trop);
                    double r_orig[5];                            // originals (i_mount .. i_top_sm)
                    for (int i = i_mount; i <= i_top_sm; i++)
                        r_orig[i - i_mount] = m.radiation.x[i][j][k];
                    for (int i = i_mount + 1; i < i_top_sm; i++)  // interior 1-2-1
                        m.radiation.x[i][j][k] = 0.25 * r_orig[i - 1 - i_mount]
                                               + 0.5  * r_orig[i - i_mount]
                                               + 0.25 * r_orig[i + 1 - i_mount];
                    if (i_top_sm > i_mount)                       // surface: one-sided blend toward air
                        m.radiation.x[i_mount][j][k] = 0.5 * r_orig[0] + 0.5 * r_orig[1];
                }

                // Sub-surface land cells. With i_mount on the ground the loops above no longer
                // write i < i_mount at all, so those cells would keep whatever the PREVIOUS call
                // left -- and t is read by the rest of the model, not just plotted. Copy the
                // ground value down, which is IceSchemeCommon::fillTopography's convention for
                // every other field. Empty loop when i_mount == 0, so the off-branch is untouched.
                for (int i = i_mount - 1; i >= 0; i--) {
                    m.t.x[i][j][k]         = m.t.x[i_mount][j][k];
                    m.radiation.x[i][j][k] = m.radiation.x[i_mount][j][k];
                    m.epsilon.x[i][j][k]   = m.epsilon.x[i_mount][j][k];
                }

                // tau_above / tau_layer are fed to NOTHING -- print_min_max_atm and the four
                // ParaView writers, and that is all; epsilon carries the physics. So their
                // ground boundary condition is applied on BOTH branches: with ATM_RAD_TOPO off
                // the loops above still walk the rock and leave an optical depth inside the
                // mountain, which is what makes the plotted field wrong there. Filling from the
                // ground costs nothing and cannot move a result.
                {
                    const int i_g = std::min(std::max(m.i_topography[j][k], 0), i_trop);
                    for (int i = i_g - 1; i >= 0; i--) {
                        m.tau_above.x[i][j][k] = m.tau_above.x[i_g][j][k];
                        m.tau_layer.x[i][j][k] = m.tau_layer.x[i_g][j][k];
                    }
                }
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
