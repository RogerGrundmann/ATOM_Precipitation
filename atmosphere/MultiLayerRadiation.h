#pragma once

#include "cAtmosphereModel.h"
#include "Utils.h"
#include "CloudFraction.h"

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
    // ATM_RAD_EQUIL -- solve the layer radiative equilibrium EXACTLY instead of with the
    // hand-rolled tridiagonal. Default 0 = the shipped path, bit-identical.
    //
    // WHY. The shipped Thomas solve does not satisfy the balance it is written to solve.
    // Measured on a US-standard column (test/rad_iterate), its output violates
    // absorbed = emitted by up to -80 W/m2 (-17 K) at 1.4 km, and applying the scheme to its
    // OWN output diverges -- surface 288.15 -> 280.76 -> 268.31 -> 243.81 -> 193.67 -> NaN in
    // five passes, with the lapse inverting by pass 3. So it has no fixed point and the single
    // pass everything was tuned against is one step of a divergent map.
    //
    // THE REPLACEMENT IS EXACT AND CHEAPER. In radiative equilibrium the net upward flux is
    // constant with height. Care with the two flux arrays: U[i] crosses the interface BELOW
    // layer i and D[i] the one ABOVE it, so U[i] - D[i] is not a net flux at one level -- the
    // net at layer i's upper interface is U[i+1] - D[i], and that is what is held at F. With
    // B_i = (U[i] + D[i])/2 (absorbed = emitted; eps cancels) and the transfer step
    // U[i+1] = U[i](1 - eps_i) + eps_i*B_i, eliminating B and D gives
    //     U[i+1] = U[i] - eps_i*F/(2 - eps_i)
    //     B_i    = U[i] - F/2 - eps_i*F/(2*(2 - eps_i))
    // integrated downward from U above the lid = F = absorbed shortwave, with B at the ground
    // equal to U[i_mount+1] because a blackbody surface does not attenuate its own emission.
    // The eps/(2 - eps) rather than eps/2 is the layer's opacity to its own emission; assuming
    // eps/2 costs 2.23 K. Verified against a 200 000-sweep Jacobi solve of the same balance:
    // max |dT| = 0.000000 K, with TOA closure OLR = SW_abs coming out to +0.000 W/m2 although
    // it is enforced nowhere.
    //
    // O(im) and non-iterative, against the shipped O(im^2) CC construction -- so this is
    // cheaper than what it replaces. On this branch the separate surface energy balance and the
    // de-kink smoother are SKIPPED: the surface is solved by the same physics as every other
    // level (more optical depth -> larger U[i_mount+1] -> warmer ground, so CO2 warms the
    // surface by construction rather than through a bolted-on Newton step), and the smoother
    // exists to hide the discontinuity that override created.
    // ATM_SW_INSOL=<solar constant W/m2> -- use the real annual-mean insolation instead of the
    // fitted pole/equator parabola. Default 0 = off, shipped path bit-identical. 1361 = Earth.
    //
    // WHY. `rad_equator_short` = 163.3 and `rad_pole_short` = 100.0 give a cos-lat GLOBAL MEAN
    // of 151.3 W/m2 against Earth's 340.3 -- the shortwave forcing is 2.25x too weak. It went
    // unnoticed because the shipped tridiagonal is biased warm by a compensating amount and was
    // never energy-closed: it emits 255 W/m2 while absorbing 150, a 105 W/m2 TOA imbalance that
    // no diagnostic in the tree reports. Two errors cancelling, the family's mue_ch4 pattern.
    //
    // AND THE PARABOLA IS THE WRONG SHAPE, NOT MERELY THE WRONG SCALE. `parabola(j/j_half)` is a
    // quadratic in LATITUDE; the annual-mean insolation is exactly a quadratic in SIN(latitude),
    //     S(phi) = (S0/4) * (1 - s2*P2(sin phi)),   P2(x) = (3x^2 - 1)/2,   s2 = 0.477,
    // which is solar geometry (the 0.477 is the obliquity's second Legendre coefficient), not a
    // fit. Putting the true endpoints 421/178 into the parabola still overshoots the global mean
    // by 8.7 % because the parabola is too flat in mid-latitudes; a least-squares parabola gets
    // the integral right only by driving the pole to 63.6 W/m2 against a true 178. So the two
    // constants cannot both be right in this form -- and in the correct form there are no free
    // constants at all, only S0.
    //
    // Endpoints come out at 421.4 (equator) and 178.0 (pole) with a global mean of 340.3 by
    // construction. `rad_equator_short` / `rad_pole_short` are ignored on this branch.
    static double swInsol(){
        static const double v = [](){
            const char* e = getenv("ATM_SW_INSOL"); return e ? atof(e) : 0.0; }();
        return v;
    }

    // ATM_CWP_CENSUS=1 -- print the DISTRIBUTION of the raw column condensate path. Print-only.
    //
    // `cwp_cap_col` scales every column's condensate to the same 20 g/m2, and its own note says
    // what that is: "a fit, not a mechanism ... the cap stands in for a cloud fraction the scheme
    // does not have, so every column is treated as fully overcast with a thin cloud rather than
    // 65 % of them with a thick one." The number that decides the replacement is the DISTRIBUTION
    // of cwp_raw across columns -- how many are cloudy at all, and how thick the cloudy ones are
    // -- and nothing in the tree prints it. A single mean cannot separate "every column at 1500"
    // from "a fifth of them at 7500", and those imply completely different cloud fractions.
    static bool cwpCensus(){
        static const bool v = [](){
            const char* e = getenv("ATM_CWP_CENSUS"); return e && atoi(e) != 0; }();
        return v;
    }

    // ATM_CLOUD_RAD_FRAC=1 -- weight the radiation by the sub-grid CLOUD FRACTION.
    // Default 0 = shipped: every column overcast.
    //
    // ATM_CLOUD_FRAC gives every cell an area fraction f and a GRID-MEAN condensate
    // q_c = f^2*D, and the radiation then reads that grid mean as if it filled the
    // cell. That is wrong in BOTH directions at once and the errors do not cancel:
    // the cover is overstated (a cell with f = 0.2 shades the whole box) while the
    // in-cloud water is understated by 1/f (the cloud that is there is five times
    // thicker than the mean). `cwp_cap_col`'s own note names this -- "the cap stands
    // in for a cloud fraction the scheme does not have, so every column is treated as
    // fully overcast with a thin cloud rather than 65 % of them with a thick one" --
    // and this is the fraction it was standing in for.
    //
    // LONGWAVE: the layer emissivity becomes the area-weighted mean of a clear and a
    // cloudy sub-column, with the cloudy one carrying the IN-CLOUD path LWP_i/f:
    //
    //     eps_i = (1-f)*(1 - exp(-tau_gas)) + f*(1 - exp(-tau_gas - tau_cloud_in))
    //
    // SHORTWAVE: MAXIMUM overlap -- the column cover is max_i f_i and the reflecting
    // sub-column carries the summed in-cloud path -- so the albedo bump becomes
    // alpha_surf + (alpha_cloud - alpha_surf)*refl*cover. Maximum overlap is the
    // conservative choice of the three: random overlap over 41 layers at f ~ 0.2 each
    // drives the cover to ~1 and would put the overcast treatment back by another
    // route. Maximum-random is the better scheme and is NOT written here.
    //
    // At f = 1 every line above is the shipped expression exactly (LWP_i/1, cover 1,
    // the eps mixture collapsing to its second term), so the off-branch is unchanged
    // by CONSTRUCTION. It is also a null unless ATM_CLOUD_FRAC is on, because the
    // shipped grid-mean closure does not produce an f -- hence the warning below.
    static bool cloudRadFrac(){
        static const bool v = [](){
            // DEFAULT ON since 2026-08-31 (the accepted configuration). Set the variable to 0 to restore the old branch.
            const char* e = getenv("ATM_CLOUD_RAD_FRAC"); return e ? atoi(e) != 0 : true; }();
        return v;
    }

    // ATM_RAD_COLDIAG=1 -- print the COLUMN DECOMPOSITION at the cell of maximum layer
    // emissivity. Print-only, default off.
    //
    // Results_Atm.cpp:40 has been printing `max epsilon` and its height in every run log for
    // months, and that print says WHERE the maximum is and nothing about WHAT it is made of.
    // Every attempt to explain it -- cloud, water vapour, the rock-humidity copy -- has been an
    // argument from the source rather than a measurement, which is the mistake this tree has now
    // recorded three times. epsilon is a sum of four terms and the print collapses them into one
    // number; this prints the four, plus the cell's water, so the attribution is read rather
    // than reasoned.
    static bool radColDiag(){
        static const bool v = [](){
            const char* e = getenv("ATM_RAD_COLDIAG"); return e && atoi(e) != 0; }();
        return v;
    }

    static bool radEquil(){
        static const bool v = [](){
            const char* e = getenv("ATM_RAD_EQUIL"); return e && atoi(e) != 0; }();
        return v;
    }

    explicit MultiLayerRadiation(cAtmosphereModel& model)
        : m(model)
    {}

    // The two repairs are ONLY correct together, and each alone is worse than shipping neither.
    // Measured on a US-standard column, iterating the scheme on its own output (test/rad_iterate):
    //
    //   arm                    pass 1    fixed point        lapse 0-9.9 km   TOA closure
    //   shipped                280.76    NaN by pass 6          3.965        never closes
    //   ATM_SW_INSOL only      297.08    140.92, still falling  5.490        --
    //   ATM_RAD_EQUIL only     269.85    219.12 snowball        6.331        +0.000
    //   BOTH                   342.02    342.02 stationary      8.024        +0.000
    //
    // A solver biased warm and a shortwave 2.25x too weak have been cancelling. Warn loudly if
    // exactly one is set, because that arm is the worst of the four and looks like a repair.
    static void warnIfHalfRepaired() {
        static bool done = false;
        if (done) return;
        done = true;
        if ((swInsol() > 0.0) != radEquil())
            std::cout << "      AGCM: [RADIATION WARNING] ATM_RAD_EQUIL and ATM_SW_INSOL are a PAIR"
                      << " and only one is set. A warm-biased solver and a 2.25x-too-weak shortwave"
                      << " have been cancelling; either alone is worse than neither." << std::endl;

        // The radiation's cloud fraction is the SAME f the condensate was made with. Read
        // without ATM_CLOUD_FRAC it would be diagnosed from a humidity field the grid-mean
        // closure never consulted, so it would rescale a condensate that does not correspond
        // to it -- a fraction invented after the fact rather than the one in use.
        if (cloudRadFrac() && !CloudFraction::enabled())
            std::cout << "      AGCM: [RADIATION WARNING] ATM_CLOUD_RAD_FRAC is set without"
                      << " ATM_CLOUD_FRAC. The radiation is weighting by a fraction that did not"
                      << " make the condensate it is weighting; set both or neither." << std::endl;
    }

    void run()
    {
        warnIfHalfRepaired();
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
        if (swInsol() > 0.0) {
            // Real annual-mean insolation; see ATM_SW_INSOL above. Latitude uses the same
            // j -> phi convention as every cos-lat weight in the tree.
            m.short_wave_radiation = std::vector<double>(m.jm, 0.0);
            for (int j = 0; j < m.jm; j++) {
                const double sp = sin(((double)j / (double)(m.jm - 1) - 0.5) * M_PI);
                const double P2 = 0.5 * (3.0 * sp * sp - 1.0);
                m.short_wave_radiation[j] = 0.25 * swInsol() * (1.0 - 0.477 * P2);
            }
        } else {
        m.short_wave_radiation = std::vector<double>(m.jm, m.rad_pole_short);
        const double rad_short_eff = m.rad_pole_short - m.rad_equator_short;
        for (int j = j_half; j >= 0; j--)
            m.short_wave_radiation[j] =
                rad_short_eff * parabola((double)j / (double)j_half) + m.rad_pole_short;
        for (int j = j_max; j > j_half; j--)
            m.short_wave_radiation[j] = m.short_wave_radiation[j_max - j];
        }

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

        std::vector<double> cwp_census;                 // ATM_CWP_CENSUS: cwp_raw per column
        if (cwpCensus()) cwp_census.assign((size_t)m.jm * m.km, 0.0);
        // Per-LEVEL condensate, so the column total can be split into "too much per level" and
        // "spread over too many levels" -- different defects with different fixes. Indexed
        // [j][i] so each OpenMP thread owns its own j and no reduction is needed.
        std::vector<double> cwp_prof;
        if (cwpCensus()) cwp_prof.assign((size_t)m.jm * m.im, 0.0);

        // ATM_RAD_COLDIAG: the maximum-emissivity cell and what its optical depth is made of.
        struct EpsMax {
            double eps = -1.0;
            int    i = 0, j = 0, k = 0, i_mount = 0;
            double t_dry = 0.0, t_wv = 0.0, t_co2 = 0.0, t_cld = 0.0;
            double cf = 1.0, cloud_scale = 1.0, LWP = 0.0, IWP = 0.0;
            double q_v = 0.0, q_c = 0.0, q_i = 0.0, p = 0.0, T = 0.0, cwp_raw = 0.0;
        } eps_max;

        // ---- per-column radiative balance (columns independent -> OpenMP over j) ----
        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < m.jm; j++) {

            // Thread-local Thomas-solver scratch; reused across k within this j.
            // Every column fully overwrites the entries it later reads (i_mount = 0,
            // i_trop = im-1), so reuse is race-free.
            std::vector<double> dp_col(m.im, 0.0), vpath_col(m.im, 0.0), wdp_col(m.im, 0.0);
            std::vector<double> alfa(m.im, 0.0), beta(m.im, 0.0);
            std::vector<double> AA(m.im, 0.0), CA(m.im, 0.0);
            std::vector<double> Uc(m.im + 1, 0.0);          // ATM_RAD_EQUIL upward-flux sweep
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
                // ATM_EPS_DRY / ATM_CO2_BAND -- the two clear-sky constants, made sweepable so
                // the offline harness (test/rad_selftest, ~0.07 s per case) can scan the scheme's
                // radiative equilibrium against a US-standard column without running the model.
                // Defaults are the shipped values, so unset is bit-identical.
                static const double eps_dry = [](){
                    const char* e = getenv("ATM_EPS_DRY");
                    const double v = e ? atof(e) : 0.684;
                    return (v > 0.0 && v < 1.0) ? v : 0.684; }();     // Bignami dry-air baseline
                double sum_dp = 0.0, sum_vp = 0.0, sum_wdp = 0.0, sum_wvp = 0.0;
                // ATM_TAU_PBROAD -- PRESSURE BROADENING of the optical-depth distribution.
                // Exponent n on (p_i/p_0); 0.0 (default) is the shipped mass-only weighting and
                // is bit-identical, 1.0 is full broadening. Real longwave absorption goes roughly
                // as p*dp, not dp, because line widths scale with pressure. Mass-only weighting
                // spreads optical depth evenly through the column, which flattens the radiative
                // equilibrium; the shipped scheme sits at 3.97 K/km against a US-standard 6.49,
                // and no setting of eps_dry or co2_band_scale reaches the latter without
                // collapsing the OLR. The column TOTAL is renormalised, so this redistributes
                // optical depth downward without adding any.
                static const double pbroad = [](){
                    const char* e = getenv("ATM_TAU_PBROAD");
                    const double v = e ? atof(e) : 0.0;
                    return (v >= 0.0 && v <= 4.0) ? v : 0.0; }();
                for (int i = i_mount; i <= i_trop; i++) {
                    double dp = (i < i_trop) ? (m.p_stat.x[i][j][k] - m.p_stat.x[i+1][j][k])
                                             : m.p_stat.x[i][j][k];   // top layer: all mass above
                    if (dp < 0.0) dp = 0.0;
                    double cw    = (m.c.x[i][j][k] > 0.0) ? m.c.x[i][j][k] : 0.0;
                    dp_col[i]    = dp;
                    vpath_col[i] = cw * dp;                          // ~ layer precipitable water
                    const double pw = (pbroad > 0.0 && m.p_0 > 0.0)
                                    ? pow(m.p_stat.x[i][j][k] / m.p_0, pbroad) : 1.0;
                    wdp_col[i]   = dp * pw;
                    sum_dp      += dp;
                    sum_vp      += vpath_col[i];
                    sum_wdp     += wdp_col[i];
                    sum_wvp     += vpath_col[i] * pw;
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
                const double inv_dp  = (sum_wdp > 0.0) ? 1.0 / sum_wdp : 0.0;
                const double inv_vp  = (sum_wvp > 0.0) ? 1.0 / sum_wvp : 0.0;

                // Physical cap on the column cloud condensate the RADIATION sees. The model
                // over-condenses (column LWP ~1500 g/m2 vs observed ~100), which saturates BOTH
                // the LW cloud greenhouse (each layer's k_liq*LWP_i ~5 -> ~opaque) AND the SW
                // albedo bump (pinned at the cloud value everywhere). Compute the raw column path
                // and a single uniform scale so the radiation treats the column as a physically
                // thick cloud (<= cwp_cap_col), preserving the vertical cloud DISTRIBUTION.
                // Applied to LWP_i/IWP_i below, it feeds BOTH the LW tau_cloud and the SW path, so
                // the two stay BALANCED — fixing only one (e.g. SW albedo alone) removes the
                // excess cooling but leaves the excess greenhouse and tips the climate hot.
                // ATM_CWP_CAP -- the column condensate cap in g/m2. **20.0 since 2026-08-28**,
                // was 250.0; ATM_CWP_CAP=250 restores the old branch exactly.
                //
                // 250 was never tunable against anything: the "OLR" printed next to it was the
                // LID TEMPERATURE (sigma*236.15^4 = 176.3 W/m2) and moved 2.7 W/m2 when cloud was
                // removed entirely. With column_olr() integrating the real upward flux, the cap
                // became measurable and 250 was three times too opaque. Sweep at nm = 21,
                // cos-lat-mean, cloud LW forcing against Earth's ~25:
                //
                //     cap    15      17      20      25      50     100     250
                //     forc 20.38   22.56   25.58   30.05   45.39   61.55   80.59
                //
                // 20 is the value, and it is CONFIRMED AT nm = 100 rather than fitted at 21 and
                // hoped for: the forcing runs 25.58239 / 25.59069 / 25.59188 / 25.59284 /
                // 25.59398 at iterations 20/40/60/80/100 -- a drift of +0.05 % over 80 iterations.
                // The worry it was run to answer -- that cloud_scale = cap/cwp_raw would over-clip
                // as the condensate field grew -- does not materialise, because cwp_raw (~1500)
                // exceeds the cap in every column, so the cap sets the column path outright and
                // only its vertical distribution can move.
                //
                // WHAT IT DOES NOT FIX, and must not be tuned to hide: clear-sky OLR is 273.9
                // against Earth's ~265, so all-sky lands at 248.3 rather than ~240. That ~9 W/m2
                // lives in the CLEAR column -- eps_dry, the Bignami 0.0056, co2_band_scale -- and
                // absorbing it into a cloud constant would make this number mean nothing.
                // AND IT IS A FIT, NOT A MECHANISM: the cap stands in for a cloud fraction the
                // scheme does not have, so every column is treated as fully overcast with a thin
                // cloud rather than 65 % of them with a thick one. The global mean is matched;
                // the regional distribution is not.
                static const double cwp_cap_col = [](){
                    // DEFAULT DISABLED since 2026-08-31: the cap was compensating for a
                    // condensate 20x too large and it INVERTS the geography (see CLAUDE.md).
                    // ATM_CWP_CAP=20 restores the shipped cap.
                    const char* e = getenv("ATM_CWP_CAP");
                    const double v = e ? atof(e) : 1.0e9;
                    return v > 0.0 ? v : 1.0e9; }();
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
                if (!cwp_census.empty()) cwp_census[(size_t)j * m.km + k] = cwp_raw;
                if (!cwp_prof.empty()) {
                    const double wj = cos((j / (double)(m.jm - 1) - 0.5) * M_PI);
                    for (int i = i_mount; i <= i_trop; i++) {
                        const double dz_i = (i < i_trop) ? (m.get_layer_height(i+1) - m.get_layer_height(i))
                                                         : (m.get_layer_height(i) - m.get_layer_height(i-1));
                        const double T_ii = m.t.x[i][j][k] * m.t_0;
                        const double rho_ii = (T_ii > 0.0) ? (m.p_stat.x[i][j][k] * 100.0) / (287.0 * T_ii) : 0.0;
                        const double cwl = (m.cloud.x[i][j][k] > 0.0) ? m.cloud.x[i][j][k] : 0.0;
                        const double cwi = (m.ice.x[i][j][k]   > 0.0) ? m.ice.x[i][j][k]   : 0.0;
                        cwp_prof[(size_t)j * m.im + i] += wj * (cwl + cwi) * rho_ii * dz_i * 1000.0;
                    }
                }
                const double cloud_scale = (cwp_raw > cwp_cap_col) ? (cwp_cap_col / cwp_raw) : 1.0;

                double lwp_col = 0.0, iwp_col = 0.0;                  // accumulated (scaled) condensate paths [g/m2] (for the SW albedo bump)
                double lwp_in_col = 0.0, iwp_in_col = 0.0;            // ATM_CLOUD_RAD_FRAC: the same paths IN-CLOUD (grid mean / f)
                double cf_col = 0.0;                                  // ATM_CLOUD_RAD_FRAC: column cover, MAXIMUM overlap
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
                    static const double co2_band_scale = [](){
                        const char* e = getenv("ATM_CO2_BAND");
                        const double v = e ? atof(e) : 0.17;
                        return (v >= 0.0) ? v : 0.17; }();
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
                    const double LWP_gm = cloud_scale * cw_l * rho_i * dz * 1000.0;  // GRID-MEAN liquid water path [g/m2], capped
                    const double IWP_gm = cloud_scale * cw_i * rho_i * dz * 1000.0;  // GRID-MEAN ice   water path [g/m2], capped

                    // ATM_CLOUD_RAD_FRAC: cloudy area fraction of this cell, from the SAME
                    // uniform-PDF closure that made the condensate (CloudFraction.h). The
                    // radiation sees the IN-CLOUD path LWP_gm/f over that fraction, not the
                    // grid mean over the whole cell. f = 1 off-branch, so LWP_i is unchanged
                    // bit for bit and every expression below collapses to the shipped one.
                    double cf = 1.0;
                    if (cloudRadFrac()) {
                        const double p_hPa = m.p_stat.x[i][j][k];
                        const double q_s   = CloudFraction::qSat(T_i, p_hPa, m.t_0, m.hp, m.ep);
                        const double q_t   = std::max(0.0, m.c.x[i][j][k]) + cw_l + cw_i;
                        cf = CloudFraction::effectiveFraction(q_t, q_s, p_hPa, cw_l + cw_i);
                    }
                    double LWP_i = LWP_gm / cf;                               // IN-CLOUD paths [g/m2]
                    double IWP_i = IWP_gm / cf;
                    double tau_cloud = k_liq * LWP_i + k_ice * IWP_i;         // IN-CLOUD optical depth

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
                    lwp_col += cf * LWP_i;  iwp_col += cf * IWP_i;                // grid-mean column paths
                    lwp_in_col += LWP_i;    iwp_in_col += IWP_i;                  // in-cloud column paths for the SW albedo bump
                    if (cf > cf_col) cf_col = cf;                                 // maximum overlap

                    const double pw_i = (dp_col[i] > 0.0) ? wdp_col[i] / dp_col[i] : 1.0;
                    const double tau_gas = tau_dry * wdp_col[i] * inv_dp
                                         + tau_wv * vpath_col[i] * pw_i * inv_vp
                                         + tau_co2;
                    double tau, eps;
                    if (cf >= 1.0) {                       // shipped path, bit for bit
                        tau = tau_gas + tau_cloud;
                        eps = 1.0 - exp(-tau);
                    } else {
                        // Area-weighted mean of a clear and a cloudy sub-column. Mixing the
                        // EMISSIVITIES and not the optical depths is the point: a partly
                        // cloudy layer is two columns seen side by side, and averaging tau
                        // instead would give one column of intermediate opacity, which emits
                        // differently because 1 - exp(-tau) is not linear in tau.
                        eps = (1.0 - cf) * (1.0 - exp(-tau_gas))
                            + cf * (1.0 - exp(-tau_gas - tau_cloud));
                        // Effective layer optical depth, for the tau_above/tau_layer
                        // instruments only. Equals tau_gas + tau_cloud when cf = 1.
                        tau = (eps < 1.0) ? -log(1.0 - eps) : (tau_gas + tau_cloud);
                    }
                    m.epsilon.x[i][j][k] = eps;

                    if (radColDiag() && eps > eps_max.eps) {
                        #pragma omp critical (rad_coldiag)
                        if (eps > eps_max.eps) {
                            eps_max.eps = eps;
                            eps_max.i = i; eps_max.j = j; eps_max.k = k; eps_max.i_mount = i_mount;
                            eps_max.t_dry = tau_dry * wdp_col[i] * inv_dp;
                            eps_max.t_wv  = tau_wv * vpath_col[i] * pw_i * inv_vp;
                            eps_max.t_co2 = tau_co2;
                            eps_max.t_cld = tau_cloud;
                            eps_max.cf = cf; eps_max.cloud_scale = cloud_scale;
                            eps_max.LWP = LWP_i; eps_max.IWP = IWP_i;
                            eps_max.q_v = m.c.x[i][j][k] * 1000.0;
                            eps_max.q_c = cw_l * 1000.0;
                            eps_max.q_i = cw_i * 1000.0;
                            eps_max.p = m.p_stat.x[i][j][k];
                            eps_max.T = T_i;
                            eps_max.cwp_raw = cwp_raw;
                        }
                    }

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
                    // ATM_CLOUD_RAD_FRAC: the reflecting sub-column carries the IN-CLOUD path
                    // and covers cf_col of the column (maximum overlap), so the bump is
                    // refl*cover rather than refl. Both reduce to the shipped values at f = 1.
                    const double cwp_sw = cloudRadFrac() ? (lwp_in_col + f_ice_sw * iwp_in_col)
                                                         : (lwp_col + f_ice_sw * iwp_col);
                    const double cover  = cloudRadFrac() ? cf_col : 1.0;
                    const double tau    = cwp_sw / cwp_tau;
                    const double refl   = tau / (tau + 2.0);                     // gentle saturation (0.5 at tau=2)
                    const double a0     = m.albedo.y[j][k];                      // surface (ice-feedback) albedo
                    if (alpha_cloud > a0)
                        m.albedo.y[j][k] = a0 + (alpha_cloud - a0) * refl * cover;
                }

                if (radEquil()) {
                    // ---- exact layer radiative equilibrium (see ATM_RAD_EQUIL above) ----
                    const double F = (1.0 - m.albedo.y[j][k]) * m.short_wave_radiation[j];
                    Uc[i_trop + 1] = F;
                    for (int i = i_trop; i >= i_mount + 1; i--) {
                        const double e = m.epsilon.x[i][j][k];
                        Uc[i] = Uc[i + 1] + e * F / (2.0 - e);
                    }
                    for (int i = i_mount + 1; i <= i_trop; i++) {
                        const double e = m.epsilon.x[i][j][k];
                        m.radiation.x[i][j][k] = Uc[i] - 0.5 * F - e * F / (2.0 * (2.0 - e));
                    }
                    m.radiation.x[i_mount][j][k] = Uc[i_mount + 1];
                    for (int i = i_mount; i <= i_trop; i++)
                        m.t.x[i][j][k] = pow(std::max(1.0, m.radiation.x[i][j][k]) / m.sigma, 0.25) / m.t_0;
                } else {

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
                //
                // THE FLOOR IS NOT COSMETIC: the Thomas back-substitution can return a NEGATIVE
                // emission on a cold, optically thin column, and pow(negative, 0.25) is NaN.
                // The identical expression 75 lines above already carries this guard; this copy
                // did not. It is latent on the shipped branch only because the initial
                // temperature is clamped at t_00 = 236.15 K, so no column is ever cold enough --
                // lower that floor to a physical tropopause value (ATM_T_FLOOR=216.65) and the
                // model NaNs at initialisation, in apply_co2_perturbation, before the time loop.
                // Another face of the register's item 3: this solver is not energy-closed, so
                // nothing prevents a negative emission. ATM_RAD_EQUIL=1 replaces it outright.
                for (int i = i_mount; i <= i_trop; i++) {
                    m.radiation.x[i][j][k] = radiation_original[i] + m.radiation.x[i][j][k];
                    m.t.x[i][j][k] = pow(std::max(1.0, m.radiation.x[i][j][k]) / m.sigma, 0.25)
                                     / m.t_0;
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
                // c_H*(T_s - T_air1) and nothing credits it to the air. A real conservation
                // defect, recorded rather than silently repaired. (`Q_Sensible`, long named here
                // as the array that should have carried it, was never that quantity -- two
                // writers, two formulas, neither a flux -- and was deleted on 2026-08-30.)

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

                }  // end of the shipped tridiagonal branch (ATM_RAD_EQUIL = 0)

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

                // Ground boundary condition for the DIAGNOSTIC radiation fields, applied on
                // BOTH branches. With ATM_RAD_TOPO off the loops above start at level 0 and walk
                // the rock, so they leave an optical depth AND an emissivity inside the mountain.
                //
                // EPSILON WAS THE ONE FIELD LEFT OUT, and it is the one that is plotted against
                // terrain. The note that used to stand here said tau is "fed to NOTHING ...
                // epsilon carries the physics", and used that to justify filling tau and not
                // epsilon. But the fill CANNOT reach the physics either way: MLR recomputes
                // epsilon from i_mount at the top of every call, before anything reads it, and
                // the only reader outside this file is column_olr, whose loop starts strictly
                // ABOVE i_topography (cAtmosphereModel.cpp) precisely so it never integrates
                // through rock. So the sub-surface values are write-only, they are what ParaView
                // draws along a mountain surface, and leaving them at the emissivity of a
                // sub-terrain "air" column is a boundary-condition gap, not a physics choice.
                // Measured on the shipped branch, 87E section: mean epsilon in rock 0.031,
                // max 0.24, against 0.023 in the ocean column at the same level.
                {
                    const int i_g = std::min(std::max(m.i_topography[j][k], 0), i_trop);
                    for (int i = i_g - 1; i >= 0; i--) {
                        m.tau_above.x[i][j][k] = m.tau_above.x[i_g][j][k];
                        m.tau_layer.x[i][j][k] = m.tau_layer.x[i_g][j][k];
                        m.epsilon.x[i][j][k]   = m.epsilon.x[i_g][j][k];
                    }
                }
            }  // k
        }  // j

        if (radColDiag() && eps_max.eps >= 0.0) {
            const double lat = 90.0 - eps_max.j;
            const double lon = (eps_max.k <= 180) ? eps_max.k : eps_max.k - 360;
            std::cout.precision(6);
            std::cout << "      [RAD COLDIAG] max epsilon = " << std::fixed << eps_max.eps
                      << "  at i=" << eps_max.i << " (" << (int)m.get_layer_height(eps_max.i) << " m)"
                      << "  lat " << lat << "  lon " << lon
                      << "  i_mount=" << eps_max.i_mount << std::endl;
            const double t_tot = eps_max.t_dry + eps_max.t_wv + eps_max.t_co2 + eps_max.t_cld;
            const double pc = (t_tot > 0.0) ? 100.0 / t_tot : 0.0;
            std::cout << "      [RAD COLDIAG] tau = " << t_tot
                      << "   dry " << eps_max.t_dry << " (" << eps_max.t_dry * pc << " %)"
                      << "   wv " << eps_max.t_wv << " (" << eps_max.t_wv * pc << " %)"
                      << "   co2 " << eps_max.t_co2 << " (" << eps_max.t_co2 * pc << " %)"
                      << "   cloud " << eps_max.t_cld << " (" << eps_max.t_cld * pc << " %)" << std::endl;
            std::cout << "      [RAD COLDIAG] q_v " << eps_max.q_v << " g/kg   q_c " << eps_max.q_c
                      << "   q_i " << eps_max.q_i
                      << "   p " << eps_max.p << " hPa   T " << eps_max.T << " K" << std::endl;
            std::cout << "      [RAD COLDIAG] cloud fraction f = " << eps_max.cf
                      << "   in-cloud LWP " << eps_max.LWP << " IWP " << eps_max.IWP << " g/m2"
                      << "   column cwp_raw " << eps_max.cwp_raw
                      << "   cwp cloud_scale " << eps_max.cloud_scale << std::endl;
        }

        // ---- ATM_CWP_CENSUS: the distribution cwp_cap_col is standing in for --------------
        if (!cwp_census.empty()) {
            // Cos-lat weighted, because a column at 80N covers far less area than one at the
            // equator and an unweighted percentile would over-count the poles.
            std::vector<std::pair<double,double> > v;                 // (cwp, weight)
            v.reserve(cwp_census.size());
            double w_tot = 0.0, w_cloudy = 0.0, mean = 0.0;
            for (int j = 0; j < m.jm; j++) {
                const double w = cos((j / (double)(m.jm - 1) - 0.5) * M_PI);
                for (int k = 0; k < m.km; k++) {
                    const double c = cwp_census[(size_t)j * m.km + k];
                    v.push_back(std::make_pair(c, w));
                    w_tot += w; mean += w * c;
                    if (c > 5.0) w_cloudy += w;                       // 5 g/m2 ~ a just-visible cloud
                }
            }
            std::sort(v.begin(), v.end());
            const double pct[] = {0.05, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99};
            double q[7] = {0,0,0,0,0,0,0};
            double acc = 0.0; int pi = 0;
            for (size_t n = 0; n < v.size() && pi < 7; n++) {
                acc += v[n].second;
                while (pi < 7 && acc >= pct[pi] * w_tot) q[pi++] = v[n].first;
            }
            const double cap = [](){ const char* e = getenv("ATM_CWP_CAP");
                                     const double x = e ? atof(e) : 1.0e9; return x > 0.0 ? x : 1.0e9; }();
            std::cout << "      AGCM: [CWP CENSUS] column condensate path g/m2, cos-lat weighted."
                      << "  mean=" << mean / w_tot
                      << "  cloudy(>5)=" << 100.0 * w_cloudy / w_tot << " %" << std::endl;
            std::cout << "      AGCM: [CWP CENSUS] p05=" << q[0] << " p25=" << q[1] << " p50=" << q[2]
                      << " p75=" << q[3] << " p90=" << q[4] << " p95=" << q[5] << " p99=" << q[6]
                      << "  max=" << v.back().first << std::endl;
            // What the cap is doing, stated as the cloud fraction it implicitly asserts: it
            // rescales every column to `cap`, so a column whose real path is cwp carries a cloud
            // of optical thickness cap/cwp of the real one. Reading that as "fraction of the sky
            // covered by a REALISTIC 100 g/m2 cloud" is the comparison the replacement needs.
            // WHERE the condensate sits. A realistic cloud occupies a LAYER; if the profile
            // is flat the column total is an integration over depth that should not exist.
            if (!cwp_prof.empty()) {
                std::vector<double> prof(m.im, 0.0);
                for (int j = 0; j < m.jm; j++)
                    for (int i = 0; i < m.im; i++) prof[i] += cwp_prof[(size_t)j * m.im + i];
                double tot = 0.0; int nlev = 0;
                for (int i = 0; i < m.im; i++) { prof[i] /= w_tot; tot += prof[i]; if (prof[i] > 0.01) nlev++; }
                std::cout << "      AGCM: [CWP CENSUS] per-level g/m2 (cos-lat mean), "
                          << nlev << " of " << m.im << " levels carrying > 0.01:" << std::endl;
                std::cout << "      AGCM: [CWP CENSUS]";
                for (int i = 0; i < m.im; i += 4)
                    std::cout << "  " << (int)m.get_layer_height(i) << "m:" << prof[i];
                std::cout << std::endl;
                std::cout << "      AGCM: [CWP CENSUS] column total from the profile = " << tot
                          << " g/m2;  mean per carrying level = " << (nlev ? tot / nlev : 0.0)
                          << " g/m2" << std::endl;
            }
            std::cout << "      AGCM: [CWP CENSUS] cap=" << cap
                      << " g/m2 applied to every column; implied cover if the cloudy part carried"
                      << " 100 g/m2: mean f = " << std::min(1.0, cap / 100.0) * 100.0
                      << " % everywhere, against a physical f = "
                      << 100.0 * std::min(1.0, (mean / w_tot) / 100.0) << " % from the mean path"
                      << std::endl;
        }

        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for MultiLayerRadiation\n", elapsed.count() * 1e-9);
        cout << "      RadiationMultiLayer ended" << endl;
    }

private:
    cAtmosphereModel& m;
};
