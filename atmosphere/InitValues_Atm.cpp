#include <queue>
#include "Knobs.h"
#include "cAtmosphereModel.h"
#include "CloudFraction.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;


// ============================================================================
// Water Vapor and Cloud Initialization Module - Final Improved Version
// ============================================================================

#include <chrono>
#include <iostream>
#include <algorithm>
#include <cmath>

// ============================================================================
// Physical and Numerical Constants
// ============================================================================
namespace VaporCloudConstants {
    // Surface evaporation coefficients
    constexpr double COEFF_LAND = 0.74;                                 // Land surface evaporation coefficient
    constexpr double COEFF_OCEAN = 0.98;                                // Ocean surface evaporation coefficient
    
    // Moisture limits
    constexpr double Q_LIMIT = 3.0e-3;                                  // Maximum specific humidity [kg/kg]
    constexpr double Q_SCALING = 0.84;                                  // Empirical moisture scaling factor
    constexpr double FALLBACK_Q_FACTOR = 1e-5;                          // Fallback saturation factor

    // Relative humidity and cloud thresholds
    constexpr double RH_THRESHOLD = 85.0;                               // Cloud formation RH threshold [%]
    constexpr double CLOUD_SCALING = 0.0005;                            // Cloud density scaling factor
    
    // Stability parameters
    constexpr double LAPSE_RATE_REF = -0.0065;                          // Reference lapse rate [K/m]
    constexpr double STABILITY_SCALING = 500.0;                         // Stability weight scaling
    constexpr double STABILITY_IMPACT = 0.3;                            // Max stability influence
    
    // Dewpoint spread threshold
    constexpr double SPREAD_THRESHOLD = 2.0;                            // Temperature-dewpoint spread [K]

    // Magnus-Tetens formula coefficients
    constexpr double MAGNUS_A_WATER = 17.2694;
    constexpr double MAGNUS_B_WATER = 35.86;
    constexpr double MAGNUS_A_ICE = 21.8747;
    constexpr double MAGNUS_B_ICE = 7.66;

    // Dewpoint calculation (inverse Magnus)
    constexpr double MAG_A = 17.27;
    constexpr double MAG_B = 237.3;
    constexpr double E0_HPA = 6.1078;                                   // Reference vapor pressure [hPa]
    
    // Safety limits
    constexpr double MIN_PRESSURE = 1e-10;
    constexpr double MIN_Q_SATUR = 1e-12;
}

using namespace VaporCloudConstants;

// ============================================================================
// Main Water Vapor and Cloud Initialization
// ============================================================================
void cAtmosphereModel::init_vapour_cloud() {                            // calculates initial water vapour and cloud distribution, no ice clouds are prepared
    std::cout << "\n\n\n      AGCM: init_vapour_cloud" << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    // Thread-local variables
    double p_u = 0.0;
    double t_u = 0.0;

    // Precompute inverse for efficiency
    const double inv_rh_range = 1.0 / (100.0 - RH_THRESHOLD);

    // ========================================================================
    // Main Computation Loop: Calculate Vapor and Cloud Fields
    // ========================================================================
    #pragma omp parallel for collapse(2) private(p_u, t_u)
    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {

            // Temperature thresholds (local constants for each thread)
            const double T_ice_end = t_000;
            const double T_freeze = t_0;

            const int i_mount = i_topography[j][k];

            // Process vertical column
            for (int i = 0; i < im-1; i++) {

                // Local variables for this grid point
                double E_sat_loc = 0.0;
                double q_sat_loc = 0.0;

                // Get local temperature and pressure
                t_u = t.x[i][j][k] * t_0;                               // [K]
                p_u = p_stat.x[i][j][k];                                // [hPa]

                // ------------------------------------------------------------
                // Calculate Saturation Vapor Pressure (Magnus-Tetens)
                // ------------------------------------------------------------
                const double E_wat = hp * AtomUtils::exp_func(t_u, MAGNUS_A_WATER, MAGNUS_B_WATER);
                const double E_ice = hp * AtomUtils::exp_func(t_u, MAGNUS_A_ICE, MAGNUS_B_ICE);
                double E_sat;
                                                                                                                                                                                                        
                // Temperature-dependent phase transition
                if (t_u >= T_freeze) {
                    E_sat = E_wat;                                      // Pure liquid water
                } else if (t_u <= T_ice_end) {
                    E_sat = E_ice;                                      // Pure ice
                } else {
                    // Linear interpolation in mixed phase region
                    const double w = (t_u - T_ice_end) / (T_freeze - T_ice_end);
                    E_sat = w * E_wat + (1.0 - w) * E_ice;
                }                                                                                                                                                                                                     

                // Saturation specific humidity with safety fallback
                const double q_Satur = (p_u > E_sat) ? 
                                       (ep * E_sat / (p_u - E_sat)) : 
                                       (ep * FALLBACK_Q_FACTOR);

                // ------------------------------------------------------------
                // Surface Evaporation (only above topography)
                // ------------------------------------------------------------
                if (i >= i_mount) {
                    const double current_coeff = is_land(h, i, j, k) ? COEFF_LAND : COEFF_OCEAN;
                    E_sat_loc = E_sat * current_coeff;
                } else {
                    E_sat_loc = 0.0;
                }

                // Calculate specific humidity from evaporation
                if (p_u > E_sat_loc && E_sat_loc > 0.0) {
                    q_sat_loc = ep * E_sat_loc / (p_u - E_sat_loc);
                } else if (E_sat_loc <= 0.0) {
                    q_sat_loc = 0.0;
                } else {
                    q_sat_loc = ep * FALLBACK_Q_FACTOR;
                }

                // ------------------------------------------------------------
                // Dewpoint Temperature and Spread
                // ------------------------------------------------------------
                // (Disabled with the dewpoint/stability block below, 2026-09-07: e_actual fed
                // only log_val, which fed only commented-out code. Restore all four together.)
//                const double e_actual   = (p_u * q_sat_loc) / (ep + q_sat_loc);
//                const double log_val    = std::log(std::max(e_actual, MIN_PRESSURE) / E0_HPA);
//                const double t_dewpoint = (MAG_B * log_val) / (MAG_A - log_val) + 273.15;
//                const double spread     = t_u - t_dewpoint;

                // ------------------------------------------------------------
                // Atmospheric Stability Assessment
                // ------------------------------------------------------------
                // Bounds check: ensure i+1 is valid
//                const double dT = (i < im-2) ? (t.x[i+1][j][k] - t.x[i][j][k]) : 0.0;
//                const double diff = dT - LAPSE_RATE_REF;

                // Stability weight: reduces moisture in stable conditions
//                const double stability_weight = 1.0 - STABILITY_IMPACT * std::tanh(diff * STABILITY_SCALING);

                // ------------------------------------------------------------
                // Cloud Formation Signal (based on relative humidity)
                // ------------------------------------------------------------
                const double rel_hum = (q_Satur > MIN_Q_SATUR) ? 
                                       (q_sat_loc / q_Satur * 100.0) : 0.0;
                double cloud_signal = (rel_hum - RH_THRESHOLD) * inv_rh_range;
                cloud_signal = std::clamp(cloud_signal, 0.0, 1.0);

                // ------------------------------------------------------------
                // Final Moisture Content (with stability and spread constraints)
                // ------------------------------------------------------------
//                const double q_final = Q_SCALING * q_sat_loc;

                // Apply moisture only when conditions are favorable:
                // - Small dewpoint spread (near saturation)
                // - Above topography
/*
                if (spread < SPREAD_THRESHOLD * stability_weight && i >= i_mount) {
                    const double q_weighted = q_final * stability_weight * cloud_signal;
                    c.x[i][j][k] = std::min(Q_LIMIT, q_weighted);
                } else {
                    c.x[i][j][k] = 0.0;
                }
*/

                const double RH_init = is_land(h, i, j, k) ? 0.60 : 0.75;
                c.x[i][j][k] = (i >= i_mount) ? RH_init * q_Satur : 0.0;                                                                                                                                         

                // Cloud density field
                cloud.x[i][j][k] = cloud_signal * CLOUD_SCALING;
                if (t_u < t_00)  cloud.x[i][j][k] = 0.0;
            }  // end i loop

            // Boundary conditions at top of domain
            c.x[im-1][j][k] = 0.0;
            cloud.x[im-1][j][k] = 0.0;

        }  // end k loop
    }  // end j loop

    // ========================================================================
    // Apply Surface Boundary Condition (copy topography values to surface)
    // ========================================================================
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {
            const int i_mount = i_topography[j][k];
            
            // Bounds check before accessing array
            if (i_mount >= 0 && i_mount < im && is_land(h, i_mount, j, k)) {
                c.x[0][j][k] = c.x[i_mount][j][k];
            }
        }
    }
/*
    // ========================================================================                                                                                                             
    // Post-processing: smooth c and cloud horizontally
    // ========================================================================                                                                                                             
    {           
        // Precompute Gaussian weights for the stencil (constant for all levels/cells)
        const int R = 2;
        const int D = 2 * R + 1;

        std::vector<double> gauss_w(D * D);

        for (int dj = -R; dj <= R; dj++)
            for (int dk = -R; dk <= R; dk++)
                gauss_w[(dj + R) * D + (dk + R)] =
                    std::exp(-0.5 * (dj*dj + dk*dk) / double(R*R));


        #pragma omp parallel
        {
            std::vector<double> tmp(jm * km);       // thread-local scratch buffer

            auto smoothLevel = [&](auto& field, int i) {
                for (int j = 0; j < jm; j++)
                    for (int k = 0; k < km; k++)
                        tmp[j * km + k] = field.x[i][j][k];

                for (int j = 0; j < jm; j++) {
                    for (int k = 0; k < km; k++) {
                        if (i < i_topography[j][k]) continue;

                        double sum = 0.0;
                        double cnt = 0.0;

                        for (int dj = -R; dj <= R; dj++) {
                            const int jj = std::clamp(j + dj, 0, jm - 1);

                            for (int dk = -R; dk <= R; dk++) {
                                const int kk = (k + dk + km) % km;
//                                const double w = std::exp(-0.5 * (dj*dj + dk*dk) / double(R*R));   // if not use Gaussian weights
                                const double w = gauss_w[(dj + R) * D + (dk + R)];

                                if (i >= i_topography[jj][kk]) {
                                    sum += w * tmp[jj * km + kk];
                                    cnt += w;
                                }
                            }
                        }
                        field.x[i][j][k] = cnt > 0.0 ? sum / cnt : 0.0;
                    }
                }
            };

            // Each i-level is independent: threads take different levels
            #pragma omp for schedule(dynamic, 4)
            for (int i = 0; i < im - 1; i++) {     // im-1 stays 0 (top BC)
                smoothLevel(c,     i);
                smoothLevel(cloud, i);
            }

        } // end omp parallel

        // Re-apply surface BC
        #pragma omp parallel for collapse(2)
        for (int j = 0; j < jm; j++) {
            for (int k = 0; k < km; k++) {
                const int i_mount = i_topography[j][k];

                if (i_mount >= 0 && i_mount < im && is_land(h, i_mount, j, k))
                    c.x[0][j][k] = c.x[i_mount][j][k];
            }
        }
    }
*/

    // ========================================================================
    // Performance Timing and Output
    // ========================================================================
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" Time measured: %.3f seconds for init_vapour_cloud\n", elapsed.count() * 1e-9);

    std::cout << "      AGCM: init_vapour_cloud ended" << std::endl;
}
/*
*
*/
// ============================================================================
// Optional: Debug Output Function (compile with -DDEBUG_VAPOR_CLOUD)
// ============================================================================
#ifdef DEBUG_VAPOR_CLOUD

void cAtmosphereModel::debug_vapor_output(int i, int j, int k, 
                                          double t_u, double p_u, 
                                          double E_sat, double q_sat) {
    if ((j == 60) && (k == 87)) {
        std::cout.precision(5);
        std::cout.setf(std::ios::fixed);
        std::cout << "\n  WaterVapour Debug Output"
                  << "\n  Position: i=" << i << " j=" << j << " k=" << k
                  << "\n  p_static = " << p_u
                  << "\n  T = " << t_u << " K  (" << t_u - t_0 << " °C)"
                  << "\n  E_sat = " << E_sat
                  << "\n  0.84 * q_sat = " << q_sat
                  << "\n  Moisture c = " << c.x[i][j][k]
                  << "\n  Diff_c = " << c.x[i][j][k] - q_sat
                  << "\n" << std::endl;
    }
}
#endif
/*
*
*/
// ==================================================================
// THE TROPOPAUSE HEIGHT IS CONVERTED TO A GRID INDEX BY DIVIDING BY L_atm
// ATM_TROPO_INDEX_FIX=1, DEFAULT 0 = SHIPPED and bit-identical unset.
//
// THE DEFECT. `round(h / L_atm)` is the conversion for a UNIFORM grid of 400 m layers, and
// cAtmosphereModel.h:695 says in as many words that L_atm = 400 m is "the AMPLITUDE OF THE
// EXPONENTIAL STRETCH, not a grid step". The grid is z_i = (exp(zeta*(r-r0)) - 1)*L_atm with
// zeta = 3.715, spread 23.21x -- 38.9 m at the bottom, over a kilometre aloft -- so the index
// this returns lands nowhere near the height asked for, and it gets worse poleward because the
// map from index to height is convex:
//
//     lat   intended    shipped index -> height        fixed index -> height
//       0     15000 m        38 -> 13239 m  (88 %)        39 -> 14567 m  ( 97 %)
//      15     14644 m        37 -> 12030 m  (82 %)        39 -> 14567 m  ( 99 %)
//      30     13671 m        34 ->  9007 m  (66 %)        38 -> 13239 m  ( 97 %)
//      45     12308 m        31 ->  6719 m  (55 %)        37 -> 12030 m  ( 98 %)
//      60     10800 m        27 ->  4510 m  (42 %)        36 -> 10927 m  (101 %)
//      75      9330 m        23 ->  2987 m  (32 %)        34 ->  9007 m  ( 97 %)
//      90      8000 m        20 ->  2163 m  (27 %)        33 ->  8173 m  (102 %)
//
// WHERE IT DOES MEASURED DAMAGE, AND IT IS THE SOUTHERN POLAR CELL.
// install_cells_from_streamfunction() confines each cell to i_topography..tropopause and ZEROES
// any column with no room between them. At 75 deg the false tropopause is 2987 m and the
// Antarctic plateau stands at 2500-4000 m, i.e. AT OR ABOVE IT, so those columns get no cell at
// all: the run log reads "installed over 61 451 columns" of 65 341, and the 3 890 missing are
// concentrated there. Measured consequences at iteration 200 (ATM_V_MASSBAL_STRIDE=1):
//
//   - the southern polar cell is 2.90e9 kg/s against 75N's 10.58e9, a factor of 3.6, because it
//     is built from only the low-lying part of its latitude circle;
//   - with a 3.6x smaller cell against a comparable absolute Psi(ground) offset, 75S is the
//     worst-closing band at EVERY mass-balance stride -- 0.0612 where the other five reach
//     0.0016 or better;
//   - and the streamfunction's own two divisor columns confirm the terrain independently:
//     max|psi_old|/max|psi_fixdiv| is 1.000 at 45N/15N/15S/45S/85N and 1.975 at 75S, 2.005 at
//     85S -- about half of those latitude circles is land.
//
// With the corrected index (34, 9007 m) every one of those Antarctic columns has 9-15 levels of
// room and carries a cell.
//
// ⚠ WHAT THE FIX COSTS, AND IT IS WHY THIS IS NOT FLIPPED ON. In the tropics the corrected
// index is 39 of 40 -- ONE LEVEL BELOW THE LID. That is not the conversion's fault: this tree's
// shell is 16 023 m and the intended equatorial tropopause is 15 000 m, so there is genuinely
// only ~1 km of stratosphere above it, and the top layer is 1456 m thick, so 0 deg and 15 deg
// both land on 39 and the tropopause becomes FLAT across the tropics in index space. Any
// consumer that treats "tropopause" as "with stratosphere above it" -- ThermoAtm's lapse
// construction, balance_thermal_wind's shear integration, install_cells_from_streamfunction's
// cell span -- meets a very different field. The index is clamped to im-2 so at least one level
// always remains above it, and no further judgement is applied.
//
// THE CONSUMERS, CORRECTED 2026-09-11 BY GREP AGAINST HEAD -- the list first written here was
// wrong in both directions, and the difference matters because it moves the fix's reach from the
// thermodynamics onto the DEFAULT initial velocity field.
//
//   - balance_thermal_wind DOES NOT READ THE TROPOPAUSE AT ALL. It was named for three sites it
//     does not have: `grep -n get_tropopause_layer VelocityInitializer.h` returns 396, 859, 879
//     and 900, and balance_thermal_wind spans 702-811.
//   - ThermoAtm.h:440 is DIAGNOSTIC ONLY. Its i_trop bounds the fill of TempStand, TempDewPoint
//     and HumidityRel, and those three reach Results_Atm's print and the VTK writers and nothing
//     dynamical. Tropopause.y is likewise VTK-only, and it carries the HEIGHT, which this does
//     not change.
//   - THE THREE REAL SITES ARE init_u, init_v_or_w AND init_v_or_w_above_tropopause -- the
//     analytic velocity IC, which is LIVE ON EVERY RUN behind no knob at all. init_u builds the
//     RADIAL profile; init_v_or_w is called for the MERIDIONAL v (the Hadley/Ferrel/polar ramps)
//     AND for the ZONAL w (the jet), so all three components move with the index.
//   - install_cells_from_streamfunction (site 396) is gated on ATM_CELLS_FROM_PSI and overwrites
//     v only, so under that knob the v-ramp change is masked and the u and w changes are not.
//
// MEASURED 2026-09-11, nm = 20 from scratch, 24 threads, ATM_CELLS_FROM_PSI=1 + _VW=0.25, one
// pinned binary, all arms exit 0 with zero NaN. AT INITIALISATION the predicted repair lands and
// one effect was unpredicted: columns given a cell 61 451 -> 64 619 of 65 341, and the initial
// max|v| 13.583 -> 3.026 m/s. That factor of 4.5 is the defect stated as a velocity -- a 26e9
// kg/s polar cell forced into a 2987 m column needs 13.6 m/s of meridional wind to carry it, and
// with 8173 m of room it needs 3.0.
//
// AT ITERATION 20, WITH THE STANDING LEAK REMOVED (ATM_V_MASSBAL_STRIDE=1, so this compares the
// CELLS and not the offset sitting on them). closure = |Psi(ground)|/cell:
//
//     closure        75N     45N     15N     15S     45S     75S
//     shipped     0.0019  0.0021  0.0012  0.0001  0.0011  0.0484
//     fixed       0.0007  0.0004  0.0041  0.0051  0.0004  0.0052
//
//     detrended cell, 1e9 kg/s
//     shipped      11.78   36.95  118.87  117.05   39.59    4.25
//     fixed        12.58   39.00  109.32  108.05   39.96   11.03
//
// 75S WAS THE OUTLIER AT EVERY STRIDE AND IS NOT ANY MORE: closure 0.0484 -> 0.0052 and its cell
// +160 %, so the N/S polar asymmetry goes 2.77x -> 1.14x. That is this knob's whole purpose,
// measured. Two costs come with it: the TROPICAL cell loses 8 % (109 against 119), because the
// same prescribed Psi is now spread over a span reaching 14 567 m; and tropical closure is ~3x
// worse (0.0041 against 0.0012) -- small in absolute terms, and three orders below the 0.163 the
// same comparison shows WITHOUT the stride, which is the leak and not this knob. The two knobs
// are coupled and must be read together.
// max w_u falls 37.608 -> 26.438 m/s on BOTH stride branches, so that is this fix and not the
// stride. The jet weakens 27.19 -> 23.56 m/s and its core rises 9007 -> 9923 m, which is init_u
// and init_v_or_w placing their tropopause values higher.
// Precipitation is IDENTICAL to the printed digit in all four arms -- 541.3 mm/a, r +0.349,
// centred RMS 1280.4, every band -- which is what 20 iterations = 4 s of physical time requires
// of any dynamical arm in this tree.
//
// ORDERING HAZARD, FIXED WITH IT: init_tropopause_layers() used to run in an
// `omp parallel sections` block ALONGSIDE init_layer_heights(). The shipped conversion reads
// only L_atm so the two were independent; this one reads m_layer_heights, which the other
// section builds. The two calls are now sequential (cAtmosphereModel.cpp) -- two O(im)/O(jm)
// loops, so the parallelism bought nothing and the race would have been real.
double cAtmosphereModel::tropopause_index(double h_m){
    // ATM_TROPO_INDEX_FIX (on since 2026-09-12) -- retired 2026-09-30 (KNOB-INV plan C); the switch and its old branch are in git history.
    if(m_layer_heights.size() < (std::size_t)im) return round(h_m / L_atm);   // not built yet
    int best = 1; double bd = 1.0e30;
    for(int i = 1; i <= im - 2; i++){                   // clamp: always leave one level above
        const double d = std::fabs((double)m_layer_heights[i] - h_m);
        if(d < bd){ bd = d; best = i; }
    }
    return (double)best;
}
/*
*
*/
void cAtmosphereModel::init_tropopause_layers(){                                                                                                                                                         
    cout << endl << endl << endl << "      AGCM: init_tropopause_layers" << endl;                                                                                                                        

    int j_max = jm - 1;                                                                                                                                                                                  
    int j_half = j_max / 2;                                                                                                                                                                              
                  
    // Derive x_max so that Agnesi(tropopause_equator, x_max) == tropopause_pole exactly.                                                                                                                
    // Agnesi: a^3/(a^2+x^2) = b  =>  x = a * sqrt(a/b - 1)
    // Requires tropopause_equator > tropopause_pole (always true physically).                                                                                                                           
    double x_max = tropopause_equator                                                                                                                                                                    
                   * std::sqrt(tropopause_equator / tropopause_pole - 1.0);                                                                                                                              
                  
  // The pole index is reported from tropopause_index(), i.e. from the conversion actually in
  // force, and with the HEIGHT that index lands on beside it. It used to print
  // round(tropopause_pole/L_atm) unconditionally, so under ATM_TROPO_INDEX_FIX=1 it reported 20
  // where the model was using 33 -- an instrument that lies on the branch under test.
  cout << "tropopause_pole=" << tropopause_pole
       << " x_max=" << x_max
       << " pole_index=" << (int)tropopause_index(tropopause_pole)
       << " -> " << get_layer_height((int)tropopause_index(tropopause_pole)) << " m"
       << "   equator_index=" << (int)tropopause_index(tropopause_equator)
       << " -> " << get_layer_height((int)tropopause_index(tropopause_equator)) << " m"
       << endl;



                                                                                                                                                                                         
    // Build symmetric cache of heights [m] and grid indices in one pass.                                                                                                                                
    std::vector<double> tropo_height_cache(jm);
    tropopause_layers = std::vector<double>(jm);                                                                                                                                                         
  
    for(int j = 0; j <= j_half; j++){                                                                                                                                                                    
        double x = x_max * (double)(j_half - j) / (double)j_half;
        double h = AtomUtils::Agnesi(tropopause_equator, x);                                                                                                                                             
        tropo_height_cache[j]       = h;
        tropo_height_cache[j_max-j] = h;                                                                                                                                                               
        tropopause_layers[j]        = tropopause_index(h);
        tropopause_layers[j_max-j]  = tropopause_layers[j];                                                                                                                                             
    }                                                                                                                                                                                                    
                                                                                                                                                                                                           
    #pragma omp parallel for schedule(static)                                                                                                                                                            
   for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){                                                                                                                                                                     
            Tropopause.y[j][k] = tropo_height_cache[j];
        }                                                                                                                                                                                                
    }           
                                                                                                                                                                                                           
    cout << "      AGCM: init_tropopause_layers ended" << endl;                                                                                                                                          
}
/*
*
*/
// ============================================================================
// Atmosphere Model Initialization - Improved Version
// ============================================================================

#include <chrono>
#include <iostream>
#include <algorithm>
#include <cmath>

// ============================================================================
// Physical Constants (should ideally be in a separate constants header)
// ============================================================================
namespace AtmosphereConstants {
    // Temperature constants
    constexpr double BETA_COSMO = 44.0;              // [K] COSMO parameter
//    constexpr double BETA_COSMO = 42.0;              // [K] COSMO parameter
//    constexpr double BETA_COSMO = 38.0;              // [K] COSMO parameter
//    constexpr double BETA_COSMO = 35.0;              // [K] COSMO parameter
//    constexpr double BETA_COSMO = 30.0;              // [K] COSMO parameter
    constexpr double MIN_SAFE_TEMP = 100.0;          // [K] Minimum safe temperature
    
    // Humidity thresholds
    constexpr double MIN_HUMIDITY = 0.0;             // [%]
    constexpr double MAX_HUMIDITY = 100.0;           // [%]
    
    // Water/air ratio factors
    constexpr double MIN_WATER_FACTOR = 0.5;         // Minimum total water factor
}

using namespace AtmosphereConstants;

// ============================================================================
// Helper Function: Project temperature to sea level (algebraic COSMO inversion)
// ============================================================================
/**
 * Inverts the COSMO temperature profile T(h) = sqrt(T0^2 - 2*beta*g*h/R)
 * exactly, so that feeding T0 back into the vertical reconstruction reproduces
 * t_u_init at h_mt without any round-trip error.
 *
 * T0 = sqrt(t_u_init^2 + 2*beta*g*h_mt / R_Air)
 *
 * @param t_u_init  Surface temperature at mountain height h_mt [K]
 * @param h_mt      Mountain height [m]
 * @return pair<T0_sea_level [K], p0_sea_level [hPa]>
 */
std::pair<double, double> project_to_sea_level(
    double t_u_init,
    double h_mt,
    double beta,
    double R_Air,
    double r_air,
    double g)
{
    double T0 = sqrt(t_u_init * t_u_init + (2.0 * beta * g * h_mt) / R_Air);
    double p0 = 1e-2 * (r_air * R_Air * T0);
    return {T0, p0};
}

// ============================================================================
// Main Initialization Function
// ============================================================================
void cAtmosphereModel::initTemperatureData(int Ma) {
    std::cout << "\n\n\n      AGCM: initTemperatureData" << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    // ========================================================================
    // Temperature Variables Declaration
    // ========================================================================
    double t_equat = 0.0; 
    double t_pole = 0.0; 
    double t_equat_add = 0.0; 
    double t_pole_add = 0.0; 
    double t_global_mean_exp = 0.0;
    double t_equat_curr = 0.0; 
    double t_pole_curr = 0.0; 
    double t_equat_prev = 0.0; 
    double t_pole_prev = 0.0;

    // ========================================================================
    // Step 1: Fix NASA Temperature Data Artifact at 180°E
    // ========================================================================
    if (is_first_time_slice()) {
        int k_half = (km - 1) / 2;

        // Interpolate bad data at dateline
        #pragma omp parallel for
        for (int j = 0; j < jm; j++) {
            temperature_NASA.y[j][k_half] = 
                0.5 * (temperature_NASA.y[j][k_half + 1] + 
                       temperature_NASA.y[j][k_half - 1]);              // [°C]
        }
    }

    // ========================================================================
    // Step 2: Initialize Temperature Field
    // ========================================================================
    const double inv_t0 = 1.0 / t_0;

    if (is_first_time_slice()) {
        // First time slice: Use NASA data directly
        #pragma omp parallel for collapse(2)
        for (int k = 0; k < km; k++) {
            for (int j = 0; j < jm; j++) {
                t.x[0][j][k] = (temperature_NASA.y[j][k] + t_0) * inv_t0;// Non-dimensional
            }
        }
    }

    // ========================================================================
    // Step 3: Apply EarthByte Reconstruction (if enabled)
    // ========================================================================
    if (!is_first_time_slice() && use_earthbyte_reconstruction) {
        #pragma omp parallel for collapse(2)
        for (int k = 0; k < km; k++) {
            for (int j = 0; j < jm; j++) {
                double val_nd = (t.x[0][j][k] + t_0) * inv_t0;
                t.x[0][j][k] = val_nd;                                  // Non-dimensional
                temp_reconst.y[j][k] = val_nd;                          // Store for later use
            }
        }
    }

    // ========================================================================
    // Step 4: Extract Temperature Values from Curves
    // ========================================================================
    t_equat_modern = get_temperatures_from_curve(0, m_equat_temperature_curve);
    t_pole_modern = get_temperatures_from_curve(0, m_pole_temperature_curve);
    t_global_mean_exp = get_temperatures_from_curve(*get_current_time(), 
                                                     m_global_temperature_curve);

    if (is_first_time_slice()) {
        t_global_mean_exp = get_temperatures_from_curve(0, m_global_temperature_curve);
        t_global_mean = GetMean_2D(jm, km, temperature_NASA);
    }

    t_global_mean = get_temperatures_from_curve(*get_current_time(), 
                                                 m_global_temperature_curve);

    // ========================================================================
    // Step 5: Calculate Temperature Increments (for non-first time slices)
    // ========================================================================
    // Current-Ma equatorial/polar temperatures from the Scotese curves. These
    // ALONE define the parabolic paleo profile applied in Step 7 — there is no
    // dependence on any foregoing Ma. Computed for every paleo slice, whether or
    // not a preceding slice exists, so a single-Ma run (time_start == time_end)
    // works without prior-slice reconstruction.
    if (*get_current_time() > 0) {
        t_equat = get_temperatures_from_curve(*get_current_time(), m_equat_temperature_curve);
        t_pole  = get_temperatures_from_curve(*get_current_time(), m_pole_temperature_curve);
    }

    // The inter-slice increments below feed ONLY the optional EarthByte
    // reconstruction correction and need a preceding slice — get_previous_time()
    // throws on the first slice — so keep them guarded by !is_first_time_slice().
    if (!is_first_time_slice()) {
        // Temperature changes between time steps
        t_equat_add = get_temperatures_from_curve(*get_current_time(), m_equat_temperature_curve)
                    - get_temperatures_from_curve(*get_previous_time(), m_equat_temperature_curve);

        t_pole_add = get_temperatures_from_curve(*get_current_time(), m_pole_temperature_curve)
                   - get_temperatures_from_curve(*get_previous_time(), m_pole_temperature_curve);

        // Current and previous values
        t_equat_curr = get_temperatures_from_curve(*get_current_time(), m_equat_temperature_curve);
        t_equat_prev = get_temperatures_from_curve(*get_previous_time(), m_equat_temperature_curve);
        t_pole_curr = get_temperatures_from_curve(*get_current_time(), m_pole_temperature_curve);
        t_pole_prev = get_temperatures_from_curve(*get_previous_time(), m_pole_temperature_curve);
    }

    // ========================================================================
    // Step 6: Print Diagnostics
    // ========================================================================
    std::cout.precision(3);
    std::cout << "\n       Time slice of Paleo-AGCM: ...................... Ma = " << Ma << " million years\n";
    std::cout << "\n       Equatorial temperature increase: ................ t_equat_add      = " << t_equat_add << " °C";
    std::cout << "\n       Polar temperature increase: ..................... t_pole_add       = " << t_pole_add << " °C";
    std::cout << "\n       Equatorial temperature at paleo times: .......... t_equat_paleo    = " << t_equat << " °C";
    std::cout << "\n       Polar temperature at paleo times: ............... t_pole_paleo     = " << t_pole << " °C";
    std::cout << "\n       Mean temperature at paleo times: ................ t_global_mean    = " << t_global_mean << " °C";
    std::cout << "\n       Expected mean temperature at paleo times: ....... t_global_mean_exp= " << t_global_mean_exp << " °C";
    std::cout << "\n       Equatorial temperature at modern times: ......... t_modern_equat   = " << t_equat_modern << " °C";
    std::cout << "\n       Polar temperature at modern times: .............. t_modern_pole    = " << t_pole_modern << " °C\n\n";

    // ========================================================================
    // Step 7: Apply Latitudinal Temperature Distribution
    // ========================================================================
    const double d_j_half = 0.5 * (jm - 1);
    const double t_0_inv = 1.0 / t_0;

    // Convert to non-dimensional
    t_equat = (t_equat + t_0) * t_0_inv;
    t_pole = (t_pole + t_0) * t_0_inv;
    t_equat_add = (t_equat_add + t_0) * t_0_inv;
    t_pole_add = (t_pole_add + t_0) * t_0_inv;

    // Effective temperature gradients
    const double delta_equat_nd = (t_equat_curr - t_equat_prev) / t_0;
    const double delta_pole_nd = (t_pole_curr - t_pole_prev) / t_0;
    const double delta_t_eff = delta_pole_nd - delta_equat_nd;
    const double t_eff = t_pole - t_equat;

    // Modern slice (Ma == 0) retains the observed NASA field assigned in Step 2;
    // paleo slices (Ma > 0) get the Scotese pole→equator parabola.
    const bool modern = (*get_current_time() == 0);

    #pragma omp parallel for collapse(2)
    for (int k = 0; k < km; k++) {
        for (int j = 0; j < jm; j++) {
            double ratio = (double)j / d_j_half;

            if (!use_earthbyte_reconstruction) {
                // Standard parabolic pole-to-pole distribution, built solely from
                // the current Ma's Scotese equator/pole temperatures — no foregoing
                // Ma required. Ma == 0 keeps the NASA field set in Step 2.
                if (!modern) {
                    t.x[0][j][k] = t_eff * AtomUtils::parabola(ratio) + t_pole;
                }
            } else {
                // EarthByte reconstruction active
                if (*get_current_time() == 0) {
                    // Initial state from NASA data
                    t.x[0][j][k] = (temperature_NASA.y[j][k] + t_0) * t_0_inv;
                } else {
                    // Apply correction to maintain latitudinal gradient evolution
                    double correction_nd = delta_t_eff * AtomUtils::parabola(ratio) + delta_pole_nd;
                    t.x[0][j][k] = temp_reconst.y[j][k] + correction_nd;
                }
            }
        }
    }

    // ========================================================================
    // Step 8: Vertical Temperature Profile & Potential Temperature
    // ========================================================================
    const double beta = BETA_COSMO;
    const double R_W_R_A = R_WaterVapour / R_Air;

    #pragma omp parallel for collapse(2)
    for (int k = 0; k < km; k++) {
        for (int j = 0; j < jm; j++) {
            int i_mount = i_topography[j][k];
            double t_u_init = t.x[0][j][k] * t_0;                       // [K]

            // OPTION B: broad land-sea thermal contrast (stacks on top of Option A).
            // The high-heat-capacity ocean stays near the zonal parabola; land equilibrates
            // closer to radiative equilibrium -> WARMER than the zonal reference in low
            // latitudes (subtropical/tropical continents -> thermal lows / monsoons) and
            // COLDER at high latitudes (continental interiors). cos(2*lat) gives +amp at the
            // equator, 0 near 45 deg, -amp at the poles. Applied to the SEA-LEVEL reference,
            // so Option A's elevation cooling still stacks on top. Paleo land only; the modern
            // NASA field and all ocean cells are untouched.
            constexpr double LANDSEA_AMP = 8.0;                        // [K] land-sea contrast amplitude (prototype, tunable)
            if (*get_current_time() != 0 && is_land(h, 0, j, k)) {
                const double lat_rad = (90.0 - (double)j * 180.0 / (double)(jm - 1)) * M_PI / 180.0;
                t_u_init += LANDSEA_AMP * cos(2.0 * lat_rad);
            }

            // Surface-temperature anchor.
            // Modern (Ma==0): the NASA field is observed AT the terrain top, so project it
            // dry-adiabatically to sea level; the COSMO column build below then reproduces it
            // exactly at i_mount (legacy behaviour, unchanged).
            // Paleo (Ma>0, OPTION A): the Scotese parabola is the SEA-LEVEL latitudinal
            // reference, NOT the mountain-top value. Anchor it at i=0 WITHOUT projecting, so
            // the COSMO vertical profile cools each column by its own DEM elevation:
            //     t.x[i_mount] = sqrt(T_sl^2 - 2*beta*g*h_mount/R)      (~5.5 K/km)
            // -> cold plateaus / warm lowlands -> zonal thermal structure that restores the
            // topographically-anchored baroclinic eddies the zonally-constant profile killed.
            // Ocean columns (i_mount=0) carry no elevation, so they stay at the parabola.
            if (*get_current_time() == 0) {                             // modern: project mountain-top NASA value up to sea level
                auto [t_pot, p_stat_0] = project_to_sea_level(
                    t_u_init,
                    get_layer_height(i_mount),
                    beta, R_Air, r_air, g
                );
                t.x[0][j][k] = t_pot;
                p_stat.x[0][j][k] = p_stat_0;
            } else {                                                    // paleo OPTION A: parabola IS the sea-level reference (no projection)
                p_stat.x[0][j][k] = 1e-2 * (r_air * R_Air * t_u_init);
                t.x[0][j][k]      = t_u_init;
            }

            temp_pot.y[j][k]     = t.x[0][j][k];
            temp_reconst.y[j][k] = t.x[0][j][k];

            // ================================================================
            // Surface Humidity Calculation (Magnus Formula)
            // ================================================================
            const double t_u_0 = t_u_init;
            const double p_u_0 = p_stat.x[0][j][k];
            
            // Magnus coefficients depend on phase (water vs ice)
            const double a_loc = (t_u_0 >= t_0) ? MAGNUS_A_WATER : MAGNUS_A_ICE;
            const double b_loc = (t_u_0 >= t_0) ? MAGNUS_B_WATER : MAGNUS_B_ICE;

            const double E_Satur = hp * exp(a_loc * (t_u_0 - t_0) / (t_u_0 - b_loc));
            const double e_curr  = c.x[0][j][k] * p_u_0 / (c.x[0][j][k] + ep);

            relative_humidity.y[j][k] = std::clamp(
                (e_curr / E_Satur) * 100.0,
                MIN_HUMIDITY,
                MAX_HUMIDITY
            );

            // ================================================================
            // Vertical Profile: Temperature, Pressure, Density
            // ================================================================
            const double t_safe = std::max(MIN_SAFE_TEMP, t.x[0][j][k]);
            const double tu_be     = t_safe / beta;
            const double c_inv_tu2 = (2.0 * beta * g) / (R_Air * t_safe * t_safe);
            const double inv_R_Air = 1.0 / R_Air;
            const double p_basis   = p_stat.x[0][j][k];

            // Removed inner #pragma omp simd to avoid nested parallelization issues
            for (int i = 0; i < im; i++) {
                const double height = get_layer_height(i);
                const double s_i    = sqrt(std::max(0.0, 1.0 - height * c_inv_tu2));

                const double t_curr = t_safe * s_i;  // [K]
                const double p_val  = p_basis * exp(-tu_be * (1.0 - s_i));

                // ATM_T_FLOOR=<K> -- the floor on the initial barometric temperature profile.
                // Default 236.15 = t_00 = shipped, bit-identical.
                //
                // THE FLOOR WAS THE HOMOGENEOUS-FREEZING POINT, AND THAT IS THE SAME CONSTANT
                // THE MOIST PHYSICS USES TO DELETE ALL CONDENSATE. t_00 is a PHASE-TRANSITION
                // temperature -- the point below which supercooled liquid cannot exist -- and
                // using it as a bound on the ATMOSPHERE guarantees the model is never colder
                // than the cutoff at which it throws its cloud away. Measured: 22.2 % of the
                // initial air column sits at EXACTLY 236.15 K, every level from 2987 m to the
                // lid in some columns and 100 % of the lid; `t_top_init` then snapshots the
                // clamped lid and BC_Atm re-imposes it every iteration, so at iteration 100 the
                // coldest cell in the whole atmosphere is still exactly -37.0000 C.
                //
                // CONSEQUENCE: NO CIRRUS IS POSSIBLE. Cirrus lives at -40 to -70 C. With the
                // floor at -37 C, q_sat aloft is 3-14x too large (measured against US-standard:
                // +21.9 K at 10.9 km, e_ice ratio 13.8x), so RH_ice reads 0.19 where the SAME
                // vapour at a correct temperature would be at 2.5 -- the model's upper-level
                // water is not too little, its saturation ceiling is too high.
                //
                // The unclamped profile is sound over the troposphere (-49.1 C at 10.9 km
                // against a US-standard -56.0) and overshoots only in the stratosphere, where it
                // has no tropopause and keeps lapsing to -86 C at the 16 km lid. So the floor is
                // doing real work and the fault is its VALUE: a tropopause temperature (~216.65,
                // the US-standard stratosphere) is the physical bound, not a freezing point.
                static const double t_floor_env = [](){
                    // DEFAULT 216.65 K since 2026-08-31 (the US-standard stratosphere).
                    // ATM_T_FLOOR=236.15 restores the shipped t_00 clamp.
                    const double v = knob::real(knob::ATM_T_FLOOR);
                    return (v > 0.0) ? v : 216.65; }();
                const double t_floor = t_floor_env;

                // The soft bound that used to stand here was DEAD: it was overwritten by the
                // hard clamp on the very next line, so the comment promising an asymptotic
                // approach described code that never ran. Removed, not revived -- removing a
                // dead store cannot change a result, and the identity test confirms it.
                t.x[i][j][k] = std::max(t_floor, t_curr);
                p_stat.x[i][j][k] = p_val;

                // Density calculations
                const double rho_base   = (p_val * 100.0 * inv_R_Air);
                const double inv_t_curr = 1.0 / std::max(MIN_SAFE_TEMP, t_curr);

                const double virtual_mult = 1.0 + (R_W_R_A - 1.0) * c.x[i][j][k];
                const double total_water_factor = 
                    std::max(MIN_WATER_FACTOR, 1.0 - (cloud.x[i][j][k] + ice.x[i][j][k]));

                const double mask = is_land(h, i, j, k) ? r_air : 1.0;

                r_dry.x[i][j][k]   = (rho_base * inv_t_curr) * mask;
                r_humid.x[i][j][k] = (rho_base / (t_curr * virtual_mult * total_water_factor)) * mask;
            }

            temp_landscape.y[j][k]   = t.x[i_mount][j][k] - t_0;        // [°C]
            p_stat_landscape.y[j][k] = p_stat.x[i_mount][j][k];
        }
    }

    // ========================================================================
    // Step 9: Convert to Non-Dimensional and Apply Boundary Conditions
    // ========================================================================
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {
            temp_reconst.y[j][k] = t.x[0][j][k] - t_0;  // [°C]

            // Convert all vertical levels to non-dimensional
            for (int i = 0; i < im; i++) {
                t.x[i][j][k] *= inv_t0;
            }

            // Apply topography boundary condition
            int i_mount = i_topography[j][k];
            if (i_mount >= 0 && i_mount < im) {
                t.x[0][j][k] = t.x[i_mount][j][k];
                p_stat.x[0][j][k] = p_stat.x[i_mount][j][k];
                r_dry.x[0][j][k] = r_dry.x[i_mount][j][k];
                r_humid.x[0][j][k] = r_humid.x[i_mount][j][k];
            }
        }
    }

    // ========================================================================
    // Step 9b (OPTION A): smooth the topography-draped surface temperature.
    // The DEM enters the lapse through the INTEGER level i_topography, so the draped
    // terrain-surface temperature has grid-scale steps (~one layer ~400 m ~2 K) at
    // cliffs/coasts -- exactly where this model's 2dx coastal/pressure modes historically
    // ignite. Apply a few light 1-2-1 passes (phi periodic, poles fixed) to the 2D surface
    // field and write it back into the PROGNOSTIC surface cell t.x[i_mount]; bcSolidGround
    // copies i_mount->0 every step, so smoothing i=0 alone would be undone on iteration 1.
    // Paleo only -- the modern NASA field and ocean cells (i_mount=0) are left untouched.
    if (*get_current_time() != 0) {
        const int n_passes = 4;
        std::vector<double> Tsurf(jm * km), Ttmp(jm * km);
        for (int j = 0; j < jm; j++)
            for (int k = 0; k < km; k++)
                Tsurf[j * km + k] = t.x[i_topography[j][k]][j][k];

        for (int p = 0; p < n_passes; p++) {
            for (int j = 0; j < jm; j++)                                // phi (k) pass, periodic
                for (int k = 0; k < km; k++) {
                    int km1 = (k - 1 + km) % km, kp1 = (k + 1) % km;
                    Ttmp[j * km + k] = 0.25 * Tsurf[j * km + km1]
                                     + 0.50 * Tsurf[j * km + k]
                                     + 0.25 * Tsurf[j * km + kp1];
                }
            for (int k = 0; k < km; k++) {                              // theta (j) pass, poles held
                Tsurf[k]              = Ttmp[k];
                Tsurf[(jm - 1) * km + k] = Ttmp[(jm - 1) * km + k];
                for (int j = 1; j < jm - 1; j++)
                    Tsurf[j * km + k] = 0.25 * Ttmp[(j - 1) * km + k]
                                      + 0.50 * Ttmp[j * km + k]
                                      + 0.25 * Ttmp[(j + 1) * km + k];
            }
        }

        for (int j = 0; j < jm; j++)
            for (int k = 0; k < km; k++) {
                int i_mount = i_topography[j][k];
                t.x[i_mount][j][k]     = Tsurf[j * km + k];             // prognostic surface cell
                t.x[0][j][k]           = Tsurf[j * km + k];             // sea-level reference layer
                temp_landscape.y[j][k] = Tsurf[j * km + k] * t_0 - t_0; // keep the diagnostic in sync [degC]
            }
    }

/*
    // ========================================================================
    // Post-processing: smooth t, p_stat, r_dry, r_humid horizontally                                        
    // ========================================================================
    {
        // Precompute Gaussian weights for the stencil (constant for all levels/cells)
        const int R = 2;
        const int D = 2 * R + 1;

        std::vector<double> gauss_w(D * D);

        for (int dj = -R; dj <= R; dj++)
            for (int dk = -R; dk <= R; dk++)

              gauss_w[(dj + R) * D + (dk + R)] =
                    std::exp(-0.5 * (dj*dj + dk*dk) / double(R*R));


        #pragma omp parallel
        {
            std::vector<double> tmp(jm * km);                           // thread-local scratch: one slice per thread

          // Lambda captures thread-local tmp — safe to call from any thread
            auto smoothLevel = [&](auto& field, int i, auto pred) {
                // Snapshot level i before writing (read-then-write on same array)
                for (int j = 0; j < jm; j++)
                    for (int k = 0; k < km; k++)
                        tmp[j * km + k] = field.x[i][j][k];

                for (int j = 0; j < jm; j++) {
                    for (int k = 0; k < km; k++) {
                        if (!pred(j, k)) continue;

                        double sum = 0.0;
                        double cnt = 0;

                        for (int dj = -R; dj <= R; dj++) {
                            const int jj = std::clamp(j + dj, 0, jm - 1);

                            for (int dk = -R; dk <= R; dk++) {
                                const int kk = (k + dk + km) % km;

                                if (i >= i_topography[jj][kk]) {
//                                    const double w = std::exp(-0.5 * (dj*dj + dk*dk) / double(R*R));   // if not use Gaussian weights
                                    const double w = gauss_w[(dj + R) * D + (dk + R)];
                                    sum += w * tmp[jj * km + kk];
                                    cnt += w;
                                }
                            }
                        }
                        field.x[i][j][k] = cnt > 0.0 ? sum / cnt : 0.0;
                    }
                }
            };

            // --- t and p_stat: above topography only ---
            // Each i-level is independent: no data dependency between levels
            #pragma omp for schedule(dynamic, 4)
            for (int i = 1; i < im - 1; i++) {
                auto above_topo = [&](int j, int k) {
                    return i >= i_topography[j][k];
                };

                smoothLevel(t,      i, above_topo);
                smoothLevel(p_stat, i, above_topo);
            }
            // implicit barrier: all t/p_stat levels done before r_dry/r_humid start

            // --- r_dry and r_humid: ocean cells only (land stays 0) ---
            #pragma omp for schedule(dynamic, 4)
            for (int i = 1; i < im - 1; i++) {
                auto ocean_above_topo = [&](int j, int k) {
                    return i >= i_topography[j][k] && !is_land(h, i, j, k);
                };

                smoothLevel(r_dry,   i, ocean_above_topo);
                smoothLevel(r_humid, i, ocean_above_topo);
            }
        } // end omp parallel — implicit barrier before BC loop

        // Re-apply i=0 boundary conditions from smoothed i_mount values
        #pragma omp parallel for collapse(2)
        for (int j = 0; j < jm; j++) {
            for (int k = 0; k < km; k++) {
                const int im0 = i_topography[j][k];

                if (im0 >= 0 && im0 < im) {
                    t.x[0][j][k]       = t.x[im0][j][k];
                    p_stat.x[0][j][k]  = p_stat.x[im0][j][k];
                    r_dry.x[0][j][k]   = r_dry.x[im0][j][k];
                    r_humid.x[0][j][k] = r_humid.x[im0][j][k];
                }

                temp_reconst.y[j][k]     = t.x[0][j][k] * t_0 - t_0;
                temp_landscape.y[j][k]   = t.x[im0][j][k] * t_0 - t_0;
                p_stat_landscape.y[j][k] = p_stat.x[im0][j][k];
            }
        }
    }
*/


    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "      Initialization completed in " << duration.count() << " ms\n";
}
/*
*
*/
void cAtmosphereModel::load_global_temperature_curve(){
    load_map_from_file(temperature_global_file, m_global_temperature_curve);
/*
    cout << "   m_global_temperature_curve" << endl;
    for(const auto &printout : m_global_temperature_curve){
        cout << printout.first << " ..... " << printout.second << '\n';
    }
*/
}
/*
*
*/
void cAtmosphereModel::load_equat_temperature_curve(){
    load_map_from_file(temperature_equat_file, m_equat_temperature_curve);
/*
    cout << "   m_equat_temperature_curve" << endl;
    for(const auto &printout : m_equat_temperature_curve){
        cout << printout.first << " ..... " << printout.second << '\n';
    }
*/
}
/*
*
*/
void cAtmosphereModel::load_pole_temperature_curve(){
    load_map_from_file(temperature_pole_file, m_pole_temperature_curve);
/*
    cout << "   m_pole_temperature_curve" << endl;
    for(const auto &printout : m_pole_temperature_curve){
        cout << printout.first << " ..... " << printout.second << '\n';
    }
*/
}
/*
*
*/
float cAtmosphereModel::get_temperatures_from_curve(float time, 
    std::map<float, float>& m) const{
    // THE SIZE TEST MUST COME FIRST. It used to sit BELOW the range test, which
    // dereferences m.begin() and decrements m.end() -- both undefined behaviour on an
    // empty map, and (--m.end()) is UB whether or not the map is empty when begin()==end().
    // The guard that was written to catch a too-small map could not run until after the
    // code it was guarding. Found in ATHAD, where the curve machinery was deleted outright
    // (one epoch, no time slices); here the curves are real, so the fix is the order.
    if(m.size() < 2){
        std::cout << "No enough data in map m" << std::endl;
        return NAN;
    }
    if(time < m.begin()->first 
        || time > (--m.end())->first){
        std::cout << "Input time out of range: " << time << std::endl;    
        return NAN;
    }
    map<float, float>::const_iterator upper = m.begin(), 
        bottom = ++m.begin(); 
    for(map<float, float>::const_iterator it = m.begin();
            it != m.end(); ++it){
        if(time < it->first){
            bottom = it;
            break;
        }else{
            upper = it;
        }
    }
/*
    std::cout << "   get_temperatures_from_curve" << std::endl;
    std::cout << "   Ma ->   " << upper->first << " " 
        << bottom->first << std::endl;
    std::cout << "   temp-range ->   "<< upper->second 
        << " " << bottom->second << std::endl;
    std::cout << "   temp-interpolation ->   " 
        << upper->second + (time - upper->first) 
       /(bottom->first - upper->first) 
        * (bottom->second - upper->second) << std::endl << std::endl;
*/
    return upper->second + (time - upper->first) 
       /(bottom->first - upper->first) 
        * (bottom->second - upper->second);
}
/*
*
*/
void cAtmosphereModel::LandOceanFraction(){
// calculation of the ratio ocean to land, also addition and subtraction of CO2 of land, ocean and vegetation

    cout << endl << endl << endl << "      AGCM: LandOceanFraction" << endl;

    int h_point_max = (jm-1) * (km-1);
    int h_land = 0;

    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            if(is_land(h, 0, j, k))  h_land = h_land + h.x[0][j][k];
        }
    }

    int h_ocean = h_point_max - h_land;
    double ocean_land = (double)h_ocean/(double)h_land;

    cout.precision(3);
    cout << endl;
    cout << setiosflags(ios::left) << setw(50) << setfill('.') 
        << "      total number of points at constant height " << " = " 
        << resetiosflags(ios::left) << setw(7) << fixed << setfill(' ') 
        << h_point_max << endl << setiosflags(ios::left) << setw(50) 
        << setfill('.') << "      number of points on the ocean surface " 
        << " = " << resetiosflags(ios::left) << setw(7) << fixed 
        << setfill(' ') << h_ocean << endl << setiosflags(ios::left) 
        << setw(50) << setfill('.') << "      number of points on the land surface " 
        << " = " << resetiosflags(ios::left) << setw(7) << fixed 
        << setfill(' ') << h_land << endl << setiosflags(ios::left) 
        << setw(50) << setfill('.') << "      ocean/land ratio " 
        << " = " << resetiosflags(ios::left) << setw(7) << fixed 
        << setfill(' ') << ocean_land 
        << endl << endl;
    cout << setiosflags(ios::left) << setw(50) << setfill('.') 
        << "      addition of CO2 by ocean surface " << " = " 
        << resetiosflags(ios::left) << setw(7) << fixed << setfill(' ') 
        << co2_ocean << endl << setiosflags(ios::left) << setw(50) 
        << setfill('.') << "      addition of CO2 by land surface " 
        << " = " << resetiosflags(ios::left) << setw(7) << fixed 
        << setfill(' ') << co2_land << endl << setiosflags(ios::left) 
        << setw(50) << setfill('.') << "      subtraction of CO2 by vegetation " 
        << " = " << resetiosflags(ios::left) << setw(7) << fixed 
        << setfill(' ') << co2_vegetation << endl << setiosflags(ios::left) 
        << setw(50) << "      valid for one single point on the surface"<< endl << endl;
    cout << endl;


    cout << "      AGCM: LandOceanFraction ended" << endl;
}
/*
*
*/
void cAtmosphereModel::initWaterWapour() {
    std::cout << "\n\n\n      AGCM: initWaterWapour" << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();


    // ========================================================================
    // Main Computation Loop: water vapour field
    // ========================================================================
    const double rh_ocean = knob::real(knob::ATM_RH_OCEAN);   // read once, outside the parallel region
    const double rh_ocean_ml = knob::real(knob::ATM_RH_OCEAN_ML);   // mixed-layer depth [m], 0 = off; read once, outside the parallel region
    const double rh_ocean_ml_s = knob::real(knob::ATM_RH_OCEAN_ML_STRENGTH);   // 0..1 partial mixing, default 1
    const double rh_land_east_ml   = knob::real(knob::ATM_RH_LAND_EAST_ML);            // mixed-layer depth above ground [m], 0 = off
    const double rh_land_east_ml_s = knob::real(knob::ATM_RH_LAND_EAST_ML_STRENGTH);   // 0..1 partial mixing, default 1
    const double rh_land_east_ml_T = knob::real(knob::ATM_RH_LAND_EAST_ML_T);          // deg C, no mixed layer on ground at least this warm
    const double rh_ocean_sst     = knob::real(knob::ATM_RH_OCEAN_SST);       // 1/K, 0 = off
    const double rh_ocean_sst_ref = knob::real(knob::ATM_RH_OCEAN_SST_REF);   // deg C
    const double rh_ocean_sst_max = knob::real(knob::ATM_RH_OCEAN_SST_MAX);
    const double rh_ocean_sst_cold = knob::real(knob::ATM_RH_OCEAN_SST_COLD); // deg C, -99 = no cut-off
    const double rh_sigma_lat = knob::real(knob::ATM_RH_SIGMA_LAT);             // deg, 0 = off; read once, outside the parallel region
    const bool   rh_sigma_lat_ocean = knob::on(knob::ATM_RH_SIGMA_LAT_OCEAN);
    // ATM_RH_LAND=<0|1> (2026-10-02), default 0 = shipped (every land column starts at 0.75 -- see the LATENT DEFECT
    // note below). [MC-LO]/radial-slice census (lcl_1, rho2_82ml, wb7): desert boundary layers start as humid as rain
    // forests (Arabia / Sahara q 17.5 g/kg at the ground, Amazon 16.2) and nothing changes that on an affordable run, so
    // the hot deserts out-rain the Amazon and Congo. =1 dries land from quantities the model knows at ANY time slice --
    // no vegetation or desert map: the subtropical DESCENT of the Hadley cell (a property of rotation, prescribed by the
    // model at every Ma) and the distance to the nearest ocean of the loaded land mask (continentality):
    //   RH = RH_wet - (RH_wet - RH_dry) * s(lat) * (0.5 + 0.5 * c(d)),
    //   s = exp(-((|lat| - 25)/10)^2),  c = 1 - exp(-d / L),  RH_wet 0.75 (today's value), RH_dry = ATM_RH_LAND_DRY,
    //   L = ATM_RH_LAND_L km.
    // Known miss: dry regions outside 15-35 deg made by monsoon dynamics (the Horn of Africa) stay moist.
    const int    rh_land_mode = knob::integer(knob::ATM_RH_LAND);
    const double rh_land_dry  = knob::real(knob::ATM_RH_LAND_DRY);
    const double rh_land_L    = knob::real(knob::ATM_RH_LAND_L) * 1.0e3;     // [m]
    std::vector<double> d_ocean;                                              // distance to the nearest ocean [m]
    if (rh_land_mode) {
        // multi-source Dijkstra on the lat-lon grid (8 neighbours, periodic in longitude), ocean columns = sources
        // (AtomLand::distanceToOcean in LandDistance.h is a copy for MoistConvection; calling it from here instead flipped
        // last bits at -O0 in the byte check run_vct.sh, so this block stays in place)
        const double a = 6.371e6, dlat = M_PI / (jm - 1), dlon = 2.0 * M_PI / (km - 1);
        d_ocean.assign((size_t)jm * km, 1e30);
        using QE = std::pair<double, int>;
        std::priority_queue<QE, std::vector<QE>, std::greater<QE>> pq;
        for (int j = 0; j < jm; j++) for (int k = 0; k < km; k++)
            if (i_topography[j][k] == 0) { d_ocean[(size_t)j * km + k] = 0.0; pq.push({0.0, j * km + k}); }
        while (!pq.empty()) {
            auto [dd, id] = pq.top(); pq.pop();
            if (dd > d_ocean[id]) continue;
            const int j = id / km, k = id % km;
            for (int dj = -1; dj <= 1; dj++) for (int dk = -1; dk <= 1; dk++) {
                if (!dj && !dk) continue;
                const int jj = j + dj; if (jj < 0 || jj >= jm) continue;
                const int kk = ((k + dk) % (km - 1) + (km - 1)) % (km - 1);
                const double cl = cos((90.0 - 0.5 * (j + jj) * 180.0 / (jm - 1)) * M_PI / 180.0);
                const double step = a * std::sqrt((dj * dlat) * (dj * dlat) + (dk * dlon * cl) * (dk * dlon * cl));
                const size_t nid = (size_t)jj * km + kk;
                if (dd + step < d_ocean[nid]) { d_ocean[nid] = dd + step; pq.push({d_ocean[nid], (int)nid}); }
            }
        }
        for (int j = 0; j < jm; j++) d_ocean[(size_t)j * km + km - 1] = d_ocean[(size_t)j * km];   // seam copy
    }
    // ATM_RH_LAND_EAST=<strength 0..1> (2026-10-03), default 0 = off (factor exactly 1). The descent drying above is
    // symmetric in longitude, so it also dries the subtropical EAST sides of the continents (India, south China, the
    // south-eastern US, eastern Australia, south-east Africa / Brazil): wb14 land 15-35 deg rains 1 mm/a against NASA 643.
    // Equatorward of ~30 deg the surface wind is the trade easterly at ANY time slice -- a property of the rotation, like
    // the descent itself -- so land with open ocean to its EAST is fed marine air and land with a continent to its east
    // is not (the west-coast deserts). m_east = exponentially weighted ocean fraction of the fetch to the east along the
    // latitude circle (e-folding 2000 km, summed to 6000 km), so a narrow sea (Red Sea, Persian Gulf) counts for little.
    // The drying is multiplied by (1 - strength * m_east). Needs only the land mask. Known miss: east coasts kept dry by
    // a coast-parallel monsoon jet (Oman / Somalia).
    const double rh_land_east = knob::real(knob::ATM_RH_LAND_EAST);
    // ATM_RH_LAND_EAST_MAX=<dRH> (2026-10-05), default -1 = no clamp (the value returned untouched). The factor (1 - strength*m_east)
    // goes NEGATIVE once strength > 1 where the fetch is nearly complete, so the land starts MOISTER than the wet end it is fed from:
    // at 1.35 (working branch, wb45) E Australia has m_east 0.85, factor -0.15, surface RH 0.846 (max 0.908) against the wet end 0.82
    // and rains 6.07 mm/d against NASA 1.86 (5.04 of it stratiform); Madagascar 7.3 / 4.3. Rain on subtropical lowland is a cliff in
    // this RH (0.80-0.84: 2.5 mm/d, 0.84-0.88: 6.7, above: 10.5; NASA 3.1-3.7), and the 16 % of the 15-35 deg land area with a
    // negative factor carries 87 % of the band's rain (python/eastfetch.py). With the knob the RH is at most wet + <dRH>.
    const double rh_land_east_max = knob::real(knob::ATM_RH_LAND_EAST_MAX);
    std::vector<double> m_east;
    if (rh_land_mode && rh_land_east > 0.0) {
        const double a = 6.371e6, dlon = 2.0 * M_PI / (km - 1), Le = 2.0e6;
        m_east.assign((size_t)jm * km, 0.0);
        for (int j = 1; j < jm - 1; j++) {
            const double dx = a * dlon * cos((90.0 - j * 180.0 / (jm - 1)) * M_PI / 180.0);
            const int n = std::min(km - 2, (int)(3.0 * Le / dx));
            for (int k = 0; k < km; k++) {
                double so = 0.0, sw = 0.0;
                for (int q = 1; q <= n; q++) {
                    const int kk = (k % (km - 1) + q) % (km - 1);
                    const double wq = std::exp(-q * dx / Le);
                    sw += wq; if (i_topography[j][kk] == 0) so += wq;
                }
                m_east[(size_t)j * km + k] = sw > 0.0 ? so / sw : 0.0;
            }
        }
    }
    auto rh_land_of = [&](int j, int k) {
        const double alat = fabs(90.0 - j * 180.0 / (double)(jm - 1));
        const double sdesc = std::exp(-((alat - 25.0) / 10.0) * ((alat - 25.0) / 10.0));
        const double c = 1.0 - std::exp(-d_ocean[(size_t)j * km + k] / rh_land_L);
        // =2: the humid end follows ATM_RH_OCEAN with the same tropical taper (full |lat| <= 30, none from 45), so wet
        // tropical land (Amazon, Congo, SE Asia) starts like the tropical ocean. [MC-RG] (rg1, mode 1): the Amazon parcel was
        // 5.7 K theta_e poorer than the tropical ocean's (q 15.0 vs 18.3 g/kg) because only the ocean got the boost.
        double wet = 0.75;
        if (rh_land_mode == 2) wet += ((alat <= 30.0) ? 1.0 : (alat >= 45.0) ? 0.0 : (45.0 - alat) / 15.0) * (rh_ocean - 0.75);
        const double east = m_east.empty() ? 1.0 : 1.0 - rh_land_east * m_east[(size_t)j * km + k];
        const double v = wet - (wet - rh_land_dry) * sdesc * (0.5 + 0.5 * c) * east;
        return (rh_land_east_max >= 0.0) ? std::min(v, wet + rh_land_east_max) : v; };
    // ATM_RH_LAND_EAST_ML=<depth m> (2026-10-06, user: "option A, write the knob and screen it"), default 0 = off (weight exactly 0,
    // block skipped, byte-identical). The land counterpart of ATM_RH_OCEAN_ML for the subtropical EAST sides that are too cool to convect.
    // WHY (SUBTROP, python/eastside.py on wb49): at the cap ATM_RH_LAND_EAST_MAX every east-side cell starts with the same surface RH
    // (0.83) and what it rains is set by the ground temperature -- at or above 24 C it convects and is on NASA (3.88 / 3.54 mm/d), below
    // 24 C the convective rain is ZERO and 1.2-1.3 mm/d of stratiform rain stand against 3.0-3.3 (S China 19.3 C 1.43 / 4.47, SE US
    // 18.6 C 1.66 / 3.93, S Brazil 22.3 C 1.09 / 4.71); the S China parcel at 28N misses buoyancy by 4.6 K (ML RH 0.71, needs 0.85).
    // The real rain there is summer convection and fronts, which an annual-mean column has not; the one route this model has to rain
    // in a cool column is the stratiform one (STORM: ATM_RH_OCEAN_ML made the storm-track rain reach the sea).
    // Below <depth> above the GROUND the vapour moves toward the ground value carried upward (constant q, RH capped 0.98) by the weight
    //   strength * m_east * lat(|lat|: 0 at <= 10, 1 from 20 to 30, 0 from 40) * cool(T_ground: 1 at <= T-2, 0 at >= T),
    // m_east the ocean fetch to the east of ATM_RH_LAND_EAST (needs that knob > 0). ATM_RH_LAND_QCAP still caps the result at the
    // zonal-mean ocean vapour at the same height above the surface. A scaffold, like ATM_RH_OCEAN_ML.
    // KNOWN MISS, by construction: nothing local separates E Australia (NASA 1.86 mm/d) from S Brazil (4.55) -- same latitude,
    // temperature and fetch -- so E Australia rises with the others.
    auto rh_land_east_ml_w = [&](int j, int k) {
        if (rh_land_east_ml <= 0.0 || m_east.empty() || i_topography[j][k] <= 0 || i_topography[j][k] >= im) return 0.0;
        const double alat = fabs(90.0 - j * 180.0 / (double)(jm - 1));
        auto sstep = [](double x) { return (x <= 0.0) ? 0.0 : (x >= 1.0) ? 1.0 : x * x * (3.0 - 2.0 * x); };
        const double T_g = t.x[i_topography[j][k]][j][k] * t_0 - 273.15;
        return rh_land_east_ml_s * m_east[(size_t)j * km + k] * sstep((alat - 10.0) / 10.0) * (1.0 - sstep((alat - 30.0) / 10.0))
               * sstep((rh_land_east_ml_T - T_g) / 2.0); };
    if (rh_land_east_ml > 0.0 && !m_east.empty()) {
        double sw = 0, ww = 0; long n = 0;
        for (int j = 0; j < jm; j++) { const double al = fabs(90.0 - j * 180.0 / (jm - 1)); if (al < 15.0 || al >= 35.0) continue;
            const double wt = cos((90.0 - j * 180.0 / (jm - 1)) * M_PI / 180.0);
            for (int k = 0; k < km - 1; k++) if (i_topography[j][k] > 0) { const double v = rh_land_east_ml_w(j, k); sw += wt * v; ww += wt; if (v > 0.05) n++; } }
        printf("      AGCM: [RH-LAND-EAST-ML] depth %.0f m, strength %.2f, T %.1f C: mean weight over 15-35 deg land %.3f, above 0.05 in %ld cells\n",
               rh_land_east_ml, rh_land_east_ml_s, rh_land_east_ml_T, ww > 0 ? sw / ww : 0.0, n);
    }
    if (rh_land_mode) {
        double s = 0, w = 0, dmax = 0; long n_dry = 0, n_land = 0;
        for (int j = 0; j < jm; j++) for (int k = 0; k < km - 1; k++) if (i_topography[j][k] > 0) {
            const double wt = cos((90.0 - j * 180.0 / (jm - 1)) * M_PI / 180.0), r = rh_land_of(j, k);
            s += wt * r; w += wt; n_land++; if (r < 0.5) n_dry++; dmax = std::max(dmax, d_ocean[(size_t)j * km + k]);
        }
        printf("      AGCM: [RH-LAND] land initial surface RH: mean %.3f, below 0.5 in %.1f %% of land cells, max distance to ocean %.0f km"
               "  (RH_dry %.2f, L %.0f km)\n", w > 0 ? s / w : 0.0, n_land ? 1e2 * n_dry / n_land : 0.0, dmax / 1e3,
               rh_land_dry, rh_land_L / 1e3);
        if (!m_east.empty()) {
            double se = 0, we = 0;
            for (int j = 0; j < jm; j++) { const double al = fabs(90.0 - j * 180.0 / (jm - 1)); if (al < 15.0 || al >= 35.0) continue;
                for (int k = 0; k < km - 1; k++) if (i_topography[j][k] > 0) { se += m_east[(size_t)j * km + k]; we += 1.0; } }
            printf("      AGCM: [RH-LAND-EAST] strength %.2f: mean ocean fetch to the east over 15-35 deg land %.3f\n",
                   rh_land_east, we > 0 ? se / we : 0.0);
        }
    }
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {
            int i_mount = i_topography[j][k];

            // ATM_RH_PROFILE -- give the initial humidity a VERTICAL PROFILE. Default 0 =
            // the shipped constant-RH column, bit-identical.
            //
            // The shipped column sets RH = RH_init at EVERY level and then multiplies by 1.25,
            // so the ocean column stands at RH = 0.9375 from the surface to the lid -- within
            // 6 % of saturation through the entire troposphere. Earth's mean RH is ~80 % in the
            // boundary layer and falls to 40-60 % in the mid-troposphere; this profile does not
            // fall, it does not vary at all.
            //
            // MEASURED CONSEQUENCE (ATM_CLOUD_INIT_DIAG + ATM_CWP_CENSUS): mean RH reaches 0.937
            // at 5 km with RH > 0.8 in 100 % of cells, while H_crit FALLS with height
            // (0.983 -> 0.801), so the two cross at ~1.4 km and above that essentially every
            // cell condenses -- cloud in 92 % of cells over 38 of 41 levels and a column
            // condensate path of 1584 g/m2 against an observed ~50-100. The per-cell values are
            // ordinary (peak 0.72 g/kg, ~0.25 g/kg at the profile peak); the excess is entirely
            // that cloud exists EVERYWHERE. So the defect is the HUMIDITY, and cwp_cap_col = 20
            // -- a factor of 79 -- has been compensating for it three modules downstream.
            //
            // The replacement is Manabe-Wetherald: RH(sigma) = RH_s*(sigma - 0.02)/0.98 with
            // sigma = p/p_0, the standard idealised profile, RH_s at the ground falling to zero
            // at the top, introducing no constant beyond the surface value already here. The
            // 1.25 multiplier goes with it: its own comment says it exists to give "a nice cloud
            // around 1 km height", i.e. a cloud deck manufactured by a fudge factor.
            // ATM_RH_PROFILE (Manabe-Wetherald profile, on since 2026-08-31) -- retired 2026-09-30 (KNOB-INV plan C); the switch and its old branch are in git history.
            // ⚠ LATENT DEFECT, NOT REPAIRED HERE (found 2026-10-02): i_mount = i_topography is the FIRST AIR level
            // (FileIO_Atm.cpp: the air cell above the land transition), so is_land(h, i_mount, j, k) is false in EVERY
            // column and the 0.60 land branch has never fired -- land starts at 0.75 like the ocean. Fixing it is a
            // separate, larger change (land 0.75 -> 0.60); the expression is left verbatim.
            double RH_init = is_land(h, i_mount, j, k) ? 0.60 : 0.75;
            // ATM_RH_OCEAN=<RH> (2026-10-02), default 0.75 = shipped: the initial surface RH over TROPICAL OCEAN
            // (i_topography == 0, full for |lat| <= 30, linear to the shipped 0.75 at 45). [BL-Q] (blq1): the tropical
            // marine BL (0-500 m, RH 0.71-0.74) is still this initial value at iteration 640 -- evaporation (~1400 mm/a)
            // needs ~2.3 days (~1e6 iterations) to change it -- and it is ~4 K theta_e short of a real tropical sounding,
            // so no surface parcel is buoyant ([MC-PB]). A scaffold in the sense ATM_RH_STORM is one. Applied as an
            // offset on the shipped value, so the default adds exactly 0.
            if (i_topography[j][k] == 0) {
                const double alat = fabs(90.0 - j * 180.0 / (double)(jm - 1));
                const double f = (alat <= 30.0) ? 1.0 : (alat >= 45.0) ? 0.0 : (45.0 - alat) / 15.0;
                RH_init += f * (rh_ocean - 0.75);
                // ATM_RH_OCEAN_SST=<1/K> (2026-10-05, user: "humidity must follow the surface temperature"), default 0 = off (block
                // skipped, byte-identical). ATM_RH_OCEAN puts ONE relative humidity on every tropical ocean column: wb35b at 220 has
                // surface RH 0.81 +- 0.01 from 18S to 18N, evaporation flat at 2.7-3.0 mm/d with no trade-wind maximum, and rain bands
                // 30-43 deg wide with no centre line against NASA's 9-30 (the equatorial cold tongue rains 2.5-4 mm/d against 1).
                // Here the surface RH is ATM_RH_OCEAN where the sea is at least ATM_RH_OCEAN_SST_REF warm and falls by <slope> per K
                // below it, by at most ATM_RH_OCEAN_SST_MAX; same tropical taper f, so the extratropical ocean (STORM) is untouched.
                // FIT (wb35b, |lat| <= 30, by SST): model / NASA rain 1.77 / 1.53 / 1.35 / 1.07 at 26-27 / 27-28 / 28-29 / 29-30 C and
                // ~0.5 below 25 C; 0.02 in RH is a factor ~1.65 in tropical-ocean rain (wb29a vs wb31a) -> slope ~0.007 /K from
                // 29.5 C, floor 0.03. A scaffold like ATM_RH_OCEAN: the real pattern also needs moisture convergence, which no
                // affordable run has (rain follows SST with r 0.88 here, 0.60 in NASA).
                // ATM_RH_OCEAN_SST_COLD=<deg C> (2026-10-05), default -99 = no cut-off (factor exactly skipped): the reduction is
                // ramped linearly from 0 at <cold> to its full value at <cold> + 1 K. WHY (wb41, slope 0.007, max 0.02): the
                // response is a cliff on cool water -- rain at 25-26 C 2.31 -> 0.16 mm/d (NASA 1.77), the all-stratiform rain below
                // 25 C 0.9 -> 0.05-0.09 (NASA 1.5-1.7), rainless tropical ocean 1 % -> 38 % (NASA 7 %), ocean 15-35 741 -> 305 --
                // where 0.02 in RH is a factor 14, not the 1.65 of the warm ocean. The excess to remove sits at 26-29 C only.
                // A FIT to this model's rain, though real surface RH is also lowest in the trade-wind belts and higher over both
                // the warm pool and the cool stratus regions.
                if (rh_ocean_sst > 0.0) {
                    const double T_s = t.x[0][j][k] * t_0 - 273.15;
                    const double dT = rh_ocean_sst_ref - T_s;
                    if (dT > 0.0) {
                        double red = std::min(rh_ocean_sst_max, rh_ocean_sst * dT);
                        if (rh_ocean_sst_cold > -90.0) {
                            const double x = T_s - rh_ocean_sst_cold;
                            red *= (x <= 0.0) ? 0.0 : (x >= 1.0) ? 1.0 : x;
                        }
                        RH_init -= f * red;
                    }
                }
            }
            if (rh_land_mode && i_topography[j][k] > 0) RH_init = rh_land_of(j, k);
            // ATM_RH_SIGMA_SFC=<0|1> (2026-10-03), default 0 = shipped. The Manabe-Wetherald profile below is
            // RH_s*(sigma - 0.02)/0.98 with sigma = p/p_0, p_0 the SEA-LEVEL reference -- so RH_s is reached only by ground
            // at sea level, and elevated ground starts drier by its own pressure: 5 % at 400 m, 7 % at 640 m. [MC-RG] (wb14):
            // Amazon / Congo / SE Asia are given the tropical ocean's surface RH (0.805) and their mixed-layer parcels hold
            // 0.72 against the ocean's 0.78; they fall 2.5 / 3.8 / 2.9 K theta_e short of the environment and do not convect
            // (land<30 CAPE 2.3 J/kg, P_conv 0.1 mm/a). MW's sigma is p over the SURFACE pressure. =1 uses the column's own
            // ground pressure p_stat(i_topography) on land; ocean columns are unchanged.
            static const bool rh_sigma_sfc = [](){ return knob::on(knob::ATM_RH_SIGMA_SFC); }();
            const double p_sfc = (rh_sigma_sfc && i_mount > 0 && i_mount < im) ? p_stat.x[i_mount][j][k] : 0.0;
            // ATM_RH_SIGMA_LAT=<deg> (2026-10-04), default 0 = off (block skipped, byte-identical): poleward of <deg> (smoothstep
            // over 10 deg) the Manabe-Wetherald sigma is p over the COLUMN'S OWN surface pressure instead of p/p_0 -- land columns,
            // and ocean columns too with ATM_RH_SIGMA_LAT_OCEAN=1. WHY (LAND-POLAR, wb34): this model's surface pressure follows the
            // surface temperature (p_sl = rho_0 R T_s), so cold columns stand at 903-911 hPa AT SEA LEVEL and p/p_0 reads them as
            // elevated: the initial surface RH is 0.76 at 60N land, 0.71 at 68N land, 0.68 at 76N ocean -- at or below H_crit(p)
            // (0.76 / 0.74 / 0.74), so the lowest 1-1.5 km is cloud-free; land 56-64N rains 0.08-0.7 mm/d against NASA 1.6-1.9,
            // land 65-90 8 mm/a against 295, ocean 65-90 62 against 428. With the column's own pressure: 0.84 / 0.80 / 0.76.
            // NOT the continentality term: ATM_RH_LAND's drying is multiplied by the descent Gaussian (25 +- 10 deg) and vanishes here.
            // The latitude limit keeps it off elevated land equatorward of it (ATM_RH_SIGMA_SFC on every land column flooded
            // the Tibetan and tropical highlands, wb15).
            double sig_lat_w = 0.0, p_col_sfc = 0.0;
            if (rh_sigma_lat > 0.0 && i_mount < im && (i_mount > 0 || rh_sigma_lat_ocean)) {
                const double x = (fabs(90.0 - j * 180.0 / (double)(jm - 1)) - rh_sigma_lat) / 10.0;
                sig_lat_w = (x <= 0.0) ? 0.0 : (x >= 1.0) ? 1.0 : x * x * (3.0 - 2.0 * x);
                p_col_sfc = p_stat.x[i_mount][j][k];
            }
            // ATM_RH_OCEAN_ML=<depth m> (2026-10-04), default 0 = off (block skipped, byte-identical). A well-mixed marine boundary
            // layer over EXTRATROPICAL ocean: below <depth> the vapour is the surface value carried upward (constant q), so RH
            // RISES with height to the 0.98 cap, instead of Manabe-Wetherald's fall with pressure. Weight 0 for |lat| <= 30,
            // smoothstep to 1 at 40 deg, ocean columns only -- the tropics, where the convection stack is tuned, are untouched.
            // WHY (STORM, wb27 at 87E): ocean 35-65 rains 138 mm/a against NASA 1107 while evaporating 579. The stratiform deck sits
            // at 1.5-4 km (RH > H_crit there); below ~1.2 km the prescribed RH falls 0.76 -> 0.67 while H_crit(p) rises to 0.91, so the
            // layer is cloud-free and evaporates the rain: 30 % reaches the sea at 40-56S, 6 % at 36S. Land at 1200 m, where the deck
            // touches the ground, keeps 96 %. A global cut of the evaporation (ATM_RAIN_AREA 0.05 / 0.02, wb28) floods the tropics.
            // A scaffold for the missing boundary-layer mixing, like ATM_RH_OCEAN: evaporation cannot build this layer in 120 s.
            double q_ml0 = -1.0, ml_w = 0.0;
            if (rh_ocean_ml > 0.0 && i_mount == 0) {
                const double alat = fabs(90.0 - j * 180.0 / (double)(jm - 1));
                const double x = (alat - 30.0) / 10.0;
                ml_w = (x <= 0.0) ? 0.0 : (x >= 1.0) ? 1.0 : x * x * (3.0 - 2.0 * x);
                // ATM_RH_OCEAN_ML_STRENGTH=<0..1>, default 1 (x1.0, the ON branch unchanged bit for bit): PARTIAL mixing -- RH moves
                // this fraction of the way from the Manabe-Wetherald value to the well-mixed one. Fully mixed (wb32a/b, 1000 / 1500 m)
                // the layer sits 0.13-0.22 above H_crit, holds 0.14-0.35 g/kg of cloud and drizzles 10 / 31 mm/d at 40-56S (ocean 35-65
                // 2984 / 9732 mm/a against NASA 1107); depth is too coarse a dose (one model level ~ +3 mm/d).
                ml_w *= rh_ocean_ml_s;
            }
            const double lml_w = rh_land_east_ml_w(j, k);            // ATM_RH_LAND_EAST_ML, see above rh_land_of's print
            double q_lml0 = -1.0;
            for (int i = 0; i < im; i++) {
                double t_u = t.x[i][j][k] * t_0;
                double p_u = p_stat.x[i][j][k];

                const double E_sat = (t_u >= t_0)
                    ? hp * AtomUtils::exp_func(t_u, MAGNUS_A_WATER, MAGNUS_B_WATER)
                    : hp * AtomUtils::exp_func(t_u, MAGNUS_A_ICE,   MAGNUS_B_ICE);

                const double q_sat = (p_u > E_sat) ? ep * E_sat / (p_u - E_sat)
                                                    : ep * FALLBACK_Q_FACTOR;

                // ATM_RH_MIN=<RH> -- floor under the Manabe-Wetherald profile aloft.
                // Default 0 = unfloored, bit-identical to the ATM_RH_PROFILE branch as shipped.
                //
                // MW is linear in sigma and therefore drives RH to ZERO at p -> 0: at sigma =
                // 0.23 (10.9 km) it prescribes 0.75*0.21 = 0.16, against an observed
                // upper-tropospheric RH over ice of 0.4-0.7. That is what forbids cirrus, and
                // it is NOT reachable by fixing the temperature: the initial humidity is set as
                // c = RH*q_sat, so cooling the column lowers q_sat and the vapour together and
                // leaves RH where it was -- measured, ATM_T_FLOOR=216.65 took 10.9 km from
                // -34.1 to -42.4 C and moved RH there from 26.6 % to 27.1 %.
                //
                // MW's own paper applies the profile to a troposphere and pins a stratospheric
                // minimum separately; taking the linear branch all the way to the lid is an
                // extrapolation past where it was meant to hold. The floor restores that.
                static const double rh_min = [](){
                    // DEFAULT 0.65 since 2026-08-31 (the TROPICAL value; see ATM_RH_MIN_LAT).
                    // ATM_RH_MIN=0 restores the unfloored Manabe-Wetherald profile.
                    const double v = knob::real(knob::ATM_RH_MIN);
                    return (v > 0.0 && v < 1.0) ? v : 0.0; }();

                // ATM_RH_MIN_LAT=1 -- give the floor the observed LATITUDE structure instead of
                // applying one number everywhere. Default 0 = the flat floor, bit-identical.
                //
                // WHY A FLAT FLOOR CANNOT WORK, MEASURED. It lifts every column to the same RH
                // aloft, so every column sits the same distance above H_crit and every column
                // makes the same cirrus: high-cloud cover comes out at 96.8 / 97.9 / 100.0 % over
                // the ATM_RH_CRIT_ICE sweep, against Earth's ~20-25 %. The model then reaches the
                // right TOTAL ice (IWP 21.4, observed 20-30) by spreading a thin cirrus over the
                // whole globe, and a global thin cirrus forces about twice what a patchy one does
                // -- 51.2 W/m2 against ~25. **This is the shipped model's "cloud in 99 % of
                // columns" pathology reproduced one layer up**: there it came from a constant-RH
                // column, here from a constant-RH FLOOR. A uniform field cannot make patchy cloud.
                //
                // The real structure is the overturning circulation: upper-tropospheric humidity
                // is high where air RISES (the ITCZ, and the mid-latitude storm track) and low
                // where it DESCENDS (the subtropical highs, the Hadley descending branch). The
                // annual-mean 300 hPa RH over ice runs ~0.60 at the equator, ~0.27 near 25 deg,
                // ~0.42 near 55 deg. ATM_RH_MIN sets the TROPICAL value and the other two follow
                // it at the observed RATIOS 1 : 0.45 : 0.70, so the knob keeps one meaning and no
                // second amplitude is introduced.
                //
                // This is a SCAFFOLD, not the physics. The humidity aloft should be carried there
                // by the model's own circulation; prescribing its latitude structure puts the
                // answer in by hand, and the tell is that these Gaussians know nothing about
                // where THIS model's ITCZ actually is. Read any cloud-cover agreement it buys as
                // assumed, not predicted.
                // ATM_RH_MIN_LAT (on since 2026-08-31) -- retired 2026-09-30 (KNOB-INV plan C); the switch and its old branch are in git history.
                double rh_floor = rh_min;
                if (rh_min > 0.0) {
                    const double phi_deg = std::fabs((j / (double)(jm - 1) - 0.5) * 180.0);
                    constexpr double phi_itcz = 12.0, w_itcz = 14.0;   // ITCZ centre / width [deg]
                    constexpr double phi_storm = 55.0, w_storm = 15.0; // storm track centre/width
                    constexpr double r_sub = 0.45, r_mid = 0.70;       // observed ratios to tropics
                    const double a = (phi_deg - 0.0) / w_itcz;
                    const double b = (phi_deg - phi_storm) / w_storm;
                    const double trop  = std::exp(-a * a);
                    const double storm = std::exp(-b * b);
                    // subtropical floor, plus a tropical bump, plus a storm-track bump
                    rh_floor = rh_min * (r_sub
                                         + (1.0   - r_sub) * trop
                                         + (r_mid - r_sub) * storm);
                    (void)phi_itcz;                                    // centre is the equator
                }
                // ATM_RH_MIN_PTOP=<hPa> -- confine the floor to the UPPER troposphere, ramped
                // in over the 200 hPa below it. Default 0 = unconfined, bit-identical.
                //
                // ATM_RH_MIN is meant to be an upper-tropospheric humidity, and unconfined it is
                // not one: Manabe-Wetherald gives RH = 0.665 at 900 hPa and 0.589 at 800, so a
                // floor of 0.65 BINDS FROM 800 hPa UPWARD -- through the entire free troposphere,
                // including the liquid deck it was never meant to touch. Measured, with the
                // latitude floor at a tropical 0.65: the cover structure comes out right
                // (subtropics 18.2 % against Earth's ~15 %, mid-latitudes 45.7 % against ~35 %,
                // global 46.4 % from 97.9 %) and the LW forcing falls 51.2 -> 30.9, but LWP goes
                // to 208 against a 50-80 band and precipitation to 1558 against NASA's 978.
                // Confining the floor separates the two: the cirrus levels keep it, the liquid
                // deck goes back to Manabe-Wetherald.
                static const double rh_min_ptop = [](){
                    // DEFAULT 490 hPa SINCE 2026-09-23 (was 475 from 2026-08-31), at the user's
                    // instruction, on the re-sweep after ATM_MC_EVAP_LIMIT stopped the drift
                    // (output_pe450/490/500/525 + output_el2 as 475, 600 from scratch). On that
                    // branch PTOP is a PURE SCALE knob: bands flat (35-65 183.5-186.0, 65-90 20.5
                    // in all five), sigma/mean flat to 0.8 %, r 0.460-0.470 with no optimum, the
                    // land/ocean RATIO invariant 0.689-0.696 (NASA 0.741). Only the mean moves:
                    // parity-mean 907.1 at 475 (-7.3 % of NASA 978.3), 965.2 at 490 (-1.3 %).
                    // A FITTED CONSTANT, like 475 was. ATM_RH_MIN_PTOP=475 restores the old
                    // default; =0 the unconfined floor.
                    // 490 -> 482 ON 2026-09-24, at the user's instruction, after ATM_RAD_TOPO went
                    // on (+5.5 % precip): output_pr465/475/482/490, 600 from scratch, parity-mean
                    // 927.4 / 953.5 / 979.4 / 1010.8 (NASA 978.3); still a pure scale knob (35-65
                    // 205.0-205.8, 65-90 24.4, land/ocean ratio 0.766-0.773). ⚠ Fitted on the branch
                    // WITHOUT ATM_WATER_CLOSURE (default since 5c7b001) and without the evaporation
                    // stride fix; re-fit on the current default is owed. =490 restores.
                    const double v = knob::real(knob::ATM_RH_MIN_PTOP);
                    return (v >= 0.0) ? v : 482.0; }();
                if (rh_min_ptop > 0.0) {
                    const double p_lo = rh_min_ptop + 200.0;          // no floor below this
                    double u = (p_lo - p_u) / (p_lo - rh_min_ptop);
                    if (u < 0.0) u = 0.0; else if (u > 1.0) u = 1.0;
                    rh_floor *= u;
                }
                double rh_i = RH_init;
                {
                    double sig = (p_0 > 0.0) ? p_u / p_0 : 1.0;
                    if (p_sfc > 0.0) sig = std::min(1.0, p_u / p_sfc);   // ATM_RH_SIGMA_SFC, land columns only
                    if (sig_lat_w > 0.0 && p_col_sfc > 0.0) sig += sig_lat_w * (std::min(1.0, p_u / p_col_sfc) - sig);   // ATM_RH_SIGMA_LAT
                    rh_i = RH_init * std::max(0.0, (sig - 0.02) / 0.98);
                    if (rh_i < rh_floor) rh_i = rh_floor;
                }
                // ATM_RH_STORM=<factor>, default 1.0 = OFF (code skipped, byte-identical).
                // Multiplies the INITIAL relative humidity by 1 + (factor-1)*G, G the storm-track
                // Gaussian of ATM_RH_MIN_LAT (55 deg, width 15 deg), capped at 0.98. A test of the
                // 2026-09-29 fill hypothesis: on the energy-conserving working branch the 35-65 deg
                // column gains +2471 mm/a of water (below saturation) while raining 69 mm/a, and
                // filling it takes ~1 day of physical time (~4.3e5 iterations). If starting it
                // moister makes the band rain from the start, the deficit is spin-up.
                static const double rh_storm = [](){
                    return knob::real(knob::ATM_RH_STORM); }();
                // ATM_RH_STORM_LAT / ATM_RH_STORM_WIDTH (2026-10-06, STORM-SHAPE), defaults 55 / 15 deg = shipped (the literal expression
                // below is kept for the default, so it is byte-identical). WHY (python/socean.py on wb57): ocean rain by latitude is
                // 424 mm/a at 34-38 deg (NASA 1000-1200), 1390 at 46-50 (1030-1230), 1540-1980 at 54-58 (1100-1160), 1280-1810 at 58-62
                // (1050-1090) -- the observed storm-track rain is flat from 34 to 62 deg, the factor peaks at 55. Lowering the peak
                // (wb60) fixes the shape (r .607 -> .617) and loses the band mean; the centre is the lever.
                static const double rh_storm_lat = [](){ return knob::real(knob::ATM_RH_STORM_LAT); }();
                static const double rh_storm_wid = [](){ return knob::real(knob::ATM_RH_STORM_WIDTH); }();
                static const bool   rh_storm_shipped_shape = (rh_storm_lat == 55.0 && rh_storm_wid == 15.0);
                if (rh_storm != 1.0) {
                    const double phi_deg = std::fabs((j / (double)(jm - 1) - 0.5) * 180.0);
                    const double b = rh_storm_shipped_shape ? (phi_deg - 55.0) / 15.0 : (phi_deg - rh_storm_lat) / rh_storm_wid;
                    rh_i *= 1.0 + (rh_storm - 1.0) * std::exp(-b * b);
                    if (rh_i > 0.98) rh_i = 0.98;
                }
                // ATM_RH_STORM_POLAR=<factor> (2026-10-06, RC-RESID 65-90 band), default 1.0 = OFF (block skipped, byte-identical).
                // Over OCEAN columns poleward of 55 deg the initial RH is multiplied by max(<factor>, the ATM_RH_STORM Gaussian) / the
                // Gaussian, i.e. the storm-track factor does not fall below <factor> toward the pole (it is 1.15 at 55 deg and 1.01 at 80).
                // WHY (python/polar.py on wb49): the 65-90 band rains 220 mm/a against NASA 364 and 85 % of the deficit is polar OCEAN
                // (S 154 / 490, N 215 / 379; N land 235 / 418; Antarctica 267 / 204 is over). 87E, Arctic Ocean 78-88N: snow forms at
                // 600-1000 m (0.28 mm/d) and 76 % of it sublimates in the cloud-free 600 m below, where the surface RH is 0.68-0.72
                // against H_crit 0.74 (S Ocean 60-66S: snow 0.79 -> 0.29 mm/d). The same loss STORM had at 35-65 deg, where the storm
                // factor and ATM_RH_OCEAN_ML closed it; the factor fades before it reaches the polar ocean. Capped at 0.98. A scaffold.
                static const double rh_storm_polar = [](){ return knob::real(knob::ATM_RH_STORM_POLAR); }();
                if (rh_storm_polar > 1.0 && i_mount == 0) {
                    const double phi_deg = std::fabs((j / (double)(jm - 1) - 0.5) * 180.0);
                    if (phi_deg > (rh_storm_shipped_shape ? 55.0 : rh_storm_lat)) {
                        const double b = rh_storm_shipped_shape ? (phi_deg - 55.0) / 15.0 : (phi_deg - rh_storm_lat) / rh_storm_wid;
                        const double gs = 1.0 + (rh_storm - 1.0) * std::exp(-b * b);
                        if (rh_storm_polar > gs) { rh_i *= rh_storm_polar / gs; if (rh_i > 0.98) rh_i = 0.98; }
                    }
                }
                if (ml_w > 0.0) {                                    // ATM_RH_OCEAN_ML, see above the loop
                    if (i == 0) q_ml0 = rh_i * q_sat;
                    else if (q_ml0 >= 0.0 && get_layer_height(i) <= rh_ocean_ml) {
                        const double rh_ml = std::min(0.98, q_ml0 / q_sat);
                        if (rh_ml > rh_i) rh_i += ml_w * (rh_ml - rh_i);
                    }
                }
                if (lml_w > 0.0 && i >= i_mount) {                   // ATM_RH_LAND_EAST_ML
                    if (i == i_mount) q_lml0 = rh_i * q_sat;
                    else if (q_lml0 >= 0.0 && get_layer_height(i) - get_layer_height(i_mount) <= rh_land_east_ml) {
                        const double rh_ml = std::min(0.98, q_lml0 / q_sat);
                        if (rh_ml > rh_i) rh_i += lml_w * (rh_ml - rh_i);
                    }
                }
                c.x[i][j][k]     = (i >= i_mount) ? rh_i * q_sat : 0.0;
//                c.x[i][j][k]     = 1.5 * c.x[i][j][k];                  // a very big cloud at 2 km and tends to reach the ground in higher latitudes
                cloud.x[i][j][k] = 0.0;
            }
        }
    }

    // ATM_RH_LAND_QCAP=<0|1> (2026-10-03), default 0 = shipped (code skipped). Land air gets its vapour from the ocean, so a
    // land boundary layer cannot hold more vapour than the marine one it is fed from; prescribing RH on land that is
    // HOTTER than the sea does exactly that. Found on the Horn of Africa (output_rg1 = wb9 state): 3-4N 42E, NASA surface
    // T 32.5 C x wet-end RH -> 22.9-23.4 g/kg against 20.1 over the adjacent ocean (28.7 C); these two cells alone
    // convect, at 96-98 mm/d (global maximum, neighbours 0); with ATM_RH_LAND=2 the whole Horn rains 75-130 mm/d
    // (wb10a/wb11/wb12a). =1 caps the initial land vapour, at each height ABOVE THE LOCAL GROUND, at the zonal-mean
    // ocean value of the same latitude at that height above the sea: q_land(i) <= <q_ocean>(i - i_topography, j).
    // RH then falls where land is hotter than the sea; land cooler than the sea, and elevated (cold) land, is untouched.
    // Paleo-safe: needs only the land mask and the model's own initial ocean column. Rows without ocean are not capped.
    if (knob::on(knob::ATM_RH_LAND_QCAP)) {
        std::vector<double> q_ocn((size_t)jm * im, -1.0);
        for (int j = 0; j < jm; j++) {
            int n = 0;
            for (int k = 0; k < km - 1; k++) if (i_topography[j][k] == 0) n++;
            if (!n) continue;
            for (int i = 0; i < im; i++) {
                double sum = 0.0;
                for (int k = 0; k < km - 1; k++) if (i_topography[j][k] == 0) sum += c.x[i][j][k];
                q_ocn[(size_t)j * im + i] = sum / n;
            }
        }
        long n_land = 0, n_cap = 0; double dq_sum = 0.0, w_sum = 0.0, dq_max = 0.0; int j_max = -1, k_max = -1;
        for (int j = 0; j < jm; j++) {
            if (q_ocn[(size_t)j * im] < 0.0) continue;
            const double wt = cos((90.0 - j * 180.0 / (jm - 1)) * M_PI / 180.0);
            for (int k = 0; k < km; k++) {
                const int i_mount = i_topography[j][k];
                if (i_mount <= 0 || i_mount >= im) continue;
                bool capped = false;
                for (int i = i_mount; i < im; i++) {
                    const double cap = q_ocn[(size_t)j * im + (i - i_mount)];
                    if (c.x[i][j][k] > cap) {
                        if (i == i_mount) { const double dq = c.x[i][j][k] - cap; capped = true;
                            if (k < km - 1) { dq_sum += wt * dq; if (dq > dq_max) { dq_max = dq; j_max = j; k_max = k; } } }
                        c.x[i][j][k] = cap;
                    }
                }
                if (k < km - 1) { n_land++; w_sum += wt; if (capped) n_cap++; }
            }
        }
        printf("      AGCM: [RH-LAND-QCAP] land surface vapour capped at the zonal-mean ocean value in %.1f %% of land columns;"
               " mean reduction %.2f g/kg (all land), max %.2f g/kg at %dN %dE\n", n_land ? 1e2 * n_cap / n_land : 0.0,
               w_sum > 0 ? 1e3 * dq_sum / w_sum : 0.0, 1e3 * dq_max, 90 - j_max, k_max <= 180 ? k_max : k_max - 360);
    }

    // ========================================================================
    // Surface Boundary Condition (copy topography values to surface)
    // ========================================================================
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {
            const int i_mount = i_topography[j][k];
            if (i_mount >= 0 && i_mount < im && is_land(h, 0, j, k)) {
                c.x[0][j][k] = c.x[i_mount][j][k];
                cloud.x[0][j][k] = cloud.x[i_mount][j][k];
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" Time measured: %.3f seconds for initWaterWapour\n", elapsed.count() * 1e-9);

    std::cout << "      AGCM: initWaterWapour ended" << std::endl;
}
/*
*
*/
void cAtmosphereModel::initCloudIce() {
    std::cout << "\n\n\n      AGCM: initCloudIce" << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    const double alfa_s    = 1.5;
    // H_crit is CloudFraction::hCrit() (ATM_RH_CRIT); the inline parabola went with ATM_CLOUD_FRAC, 2026-09-30.
//    const double det_T_0   = t_0 - 3.0;
    const double det_T_0   = t_0;

    // ========================================================================
    // Pass 1: cloud_max[i] — parallel over i, sum thread-local
    // ========================================================================
    cloud_max = std::vector<double>(im, 0.0);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < im; i++) {
        double sum = 0.0;
        for (int j = 0; j < jm; j++) {
            for (int k = 0; k < km; k++) {
                const double t_u = t.x[i][j][k] * t_0;
                const double p_u = p_stat.x[i][j][k];

                const double E_sat = (t_u >= t_0)
                    ? hp * AtomUtils::exp_func(t_u, MAGNUS_A_WATER, MAGNUS_B_WATER)
                    : hp * AtomUtils::exp_func(t_u, MAGNUS_A_ICE,   MAGNUS_B_ICE);
                const double q_sat = ep * E_sat / (p_u - E_sat);

//                sum += std::max(0.0, c.x[i][j][k] - q_sat); }
//                sum += std::max(0.0, c.x[i][j][k] - 0.84 * q_sat); }
                sum += std::max(0.0, c.x[i][j][k] - 0.74 * q_sat); }
        }
        cloud_max[i] = sum / ((jm-1) * (km-1));
        if (cloud_max[i] <= 1e-4)  cloud_max[i] = 1e-4;
    }

    // ========================================================================
    // Precompute per-level quantities (sequential: im is small)
    // ========================================================================
    std::vector<double> step(im, 0.0);
    std::vector<double> dt_dim(im, 0.0);
    std::vector<double> alfa_over_cmax(im, 0.0);
    std::vector<double> two_step(im, 0.0);

    for (int i = 0; i < im - 1; i++) {
        step[i]           = get_layer_height(i+1) - get_layer_height(i);
        dt_dim[i]         = step[i] / 1.6;
        alfa_over_cmax[i] = alfa_s / cloud_max[i];
        two_step[i]       = 2.0 * step[i];
    }
    step[im-1]           = step[im-2];
    dt_dim[im-1]         = dt_dim[im-2];
    alfa_over_cmax[im-1] = alfa_over_cmax[im-2];
    two_step[im-1]       = two_step[im-2];

    // ========================================================================
    // Pass 2: cloud, ice, graupel fields — collapse(j, k), i inner
    // ========================================================================
    #pragma omp parallel for collapse(2) schedule(static)
    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {
            for (int i = 0; i < im; i++) {
                const double t_u = t.x[i][j][k] * t_0;
                const double p_u = p_stat.x[i][j][k];

                const double E_sat = (t_u >= t_0)
                    ? hp * AtomUtils::exp_func(t_u, MAGNUS_A_WATER, MAGNUS_B_WATER)
                    : hp * AtomUtils::exp_func(t_u, MAGNUS_A_ICE,   MAGNUS_B_ICE);
                const double q_sat = ep * E_sat / (p_u - E_sat);

                // ONE H_crit CURVE, NOT TWO. The parabola is duplicated here and in
                // CloudFraction.h, and the header's own note says a third copy is where the
                // three modules would drift apart -- they already had, because the flattening
                // above p_mid that lets cloud exist aloft went into the header only. On the
                // fractional branch take the header's curve; the shipped branch keeps its own,
                // bit for bit.
                double H_crit = CloudFraction::hCrit(CloudFraction::pEff(*this, p_u, j, k));   // pEff: ATM_HCRIT_SFC, == p_u when off
                if (H_crit > 1.0)  H_crit = 1.0;

                // ---- ATM_CLOUD_FRAC: sub-grid cloud fraction (default 0 = shipped) ----
                //
                // The shipped closure diagnoses condensate from the GRID-MEAN supersaturation
                // and caps it at cloud_max[i], itself the horizontal MEAN of
                // max(0, c - 0.74*q_sat). In a realistic atmosphere that mean is essentially
                // never positive, so cloud_max floors at 1e-4 and almost nothing condenses --
                // measured: with a Manabe-Wetherald humidity the column path collapses to
                // 0.0006 g/m2 at ANY threshold. The shipped model only makes cloud because its
                // humidity is held at RH = 0.9375 everywhere. Real cloud forms from SUB-GRID
                // variability: parts of a cell saturated while the mean is not.
                //
                // Uniform-PDF (Smith / Sundqvist) closure. Total water is uniform over
                // [qbar - D, qbar + D] with D = (1 - H_crit)*q_sat, cloud where q_t > q_sat:
                //     s   = qbar + D - q_sat = q_sat*(RH - H_crit)
                //     f   = clamp(s / 2D, 0, 1)
                //     q_c = f^2 * D           (f < 1),      qbar - q_sat   (f = 1)
                // f = 0 exactly at RH = H_crit and rises smoothly; q_c is the GRID-MEAN
                // condensate, so the in-cloud value is q_c/f and stays physical as f -> 0.
                // One parameter, H_crit, which already exists. No new field: f is recomputed
                // where needed rather than stored, so sizeof(cAtmosphereModel) is untouched and
                // the stack-canary hazard in the README does not apply.
                // ATM_CLOUD_FRAC retired 2026-09-30 (KNOB-INV plan C); the switch and its old branch are in git history.
                double cloud_ls;
                {
                    const double D = (1.0 - H_crit) * q_sat;
                    if (D > 0.0) {
                        const double sfrac = c.x[i][j][k] + D - q_sat;
                        double f = sfrac / (2.0 * D);
                        if (f < 0.0) f = 0.0; else if (f > 1.0) f = 1.0;
                        cloud_ls = (f >= 1.0) ? std::max(0.0, c.x[i][j][k] - q_sat) : f * f * D;
                    } else cloud_ls = std::max(0.0, c.x[i][j][k] - q_sat);
                }

                double cloud_conv = 0.0;
                if (P_rain.x[i][j][k] > 0.0 && i < im - 1) {
                    const double del_q_conv = std::max(0.0,
                        (P_rain.x[i+1][j][k] - P_rain.x[i][j][k])
                        / (two_step[i] * r_humid.x[i][j][k]) * dt_dim[i]);
                    cloud_conv = cloud_max[i] * (1.0 - exp(-alfa_over_cmax[i] * del_q_conv));
                }

                double cloud_val = cloud_ls + cloud_conv;
                if (is_land(h, i, j, k))  cloud_val = 0.0;
                cloud.x[i][j][k] = cloud_val;
 
                double h_T = 0.0;
                if (t_u < t_0) {
                    const double ratio = (t_u - t_0) / det_T_0;
                    h_T = 1.0 - exp(-0.5 * ratio * ratio);
                }

                ice.x[i][j][k] = cloud_val * h_T;
                gr.x[i][j][k]  = 0.1 * cloud_val * h_T;
 
                if (t_u <= t_00) {
                    // ATM_ICE_COLD: liquid cannot exist below the homogeneous-freezing point;
                    // ice can, and must -- this is the cirrus range.
                    ice.x[i][j][k]  += cloud.x[i][j][k];
                    gr.x[i][j][k]    = 0.1 * ice.x[i][j][k];
                    cloud.x[i][j][k] = 0.0;
                }
            }  // i
        }  // k
    }  // j

    // ========================================================================
    // Surface Boundary Condition
    // ========================================================================
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {
            const int i_mount = i_topography[j][k];
//            if (i_mount >= 0 && i_mount < im && is_land(h, 0, j, k)) {
            if (i_mount >= 0 && i_mount < im && is_land(h, i_mount, j, k)) {
                cloud.x[0][j][k] = cloud.x[i_mount][j][k];
                ice.x[0][j][k]   = ice.x[i_mount][j][k];
                gr.x[0][j][k]    = gr.x[i_mount][j][k];
            }
            for (int i = i_mount-1; i >= 0; i--) {
                if (is_land(h, i, j, k)) {
                    cloud.x[i][j][k] = 0.0;
                    ice.x[i][j][k]   = 0.0;
                    gr.x[i][j][k]    = 0.0;
                }
            }
        }
    }

    // ---- ATM_CLOUD_INIT_DIAG: which term sets the condensate magnitude? -------------------
    //
    // The column condensate path is 1584 g/m2 in 99 % of columns against an observed ~50-100,
    // spread over 38 of 41 levels, while the PER-CELL values are physically ordinary (peak
    // 0.72 g/kg, mean per carrying level ~0.09 g/kg). So the excess is not "too much in a
    // cloud", it is "cloud everywhere". Two terms could set that and they need different fixes:
    //
    //   cloud = cloud_max[i] * (1 - exp(-alfa_s * del_q / cloud_max[i]))
    //
    // For del_q << cloud_max/alfa_s this is LINEAR, cloud ~ alfa_s*del_q, and the supersaturation
    // sets the amount. For del_q >> cloud_max/alfa_s it SATURATES at cloud_max[i], and the
    // ceiling sets it -- a ceiling that is itself the horizontal MEAN of max(0, c - 0.74*q_sat)
    // over the level, with a threshold that DISAGREES with the 0.8-1.0 H_crit used two passes
    // below. Printing the ratio cloud/cloud_max per level says which regime the field is in.
    // Print-only, default off.
    if (knob::on(knob::ATM_CLOUD_INIT_DIAG)) {
        std::cout << "      AGCM: [CLOUD-INIT] lvl      z[m]   cloud_max[g/kg]   mean del_q   "
                  << "mean cloud   mean RH   H_crit   RH>0.8   cells with cloud" << std::endl;
        for (int i = 0; i < im; i += 2) {
            double sum_dq = 0.0, sum_cl = 0.0, w_tot = 0.0, sum_rh = 0.0, sum_hc = 0.0;
            long n_cl = 0, n_tot = 0, n_rh80 = 0;
            for (int j = 0; j < jm; j++) {
                const double w = cos((j / (double)(jm - 1) - 0.5) * M_PI);
                for (int k = 0; k < km; k++) {
                    const double t_u = t.x[i][j][k] * t_0, p_u = p_stat.x[i][j][k];
                    const double E_sat = (t_u >= t_0)
                        ? hp * AtomUtils::exp_func(t_u, MAGNUS_A_WATER, MAGNUS_B_WATER)
                        : hp * AtomUtils::exp_func(t_u, MAGNUS_A_ICE,   MAGNUS_B_ICE);
                    const double q_sat  = ep * E_sat / (p_u - E_sat);
                    double H_crit = CloudFraction::hCrit(CloudFraction::pEff(*this, p_u, j, k));        // the curve the field was built with
                    if (H_crit > 1.0) H_crit = 1.0;
                    if (q_sat > 0.0) { sum_rh += w * (c.x[i][j][k] / q_sat);
                                       if (c.x[i][j][k] > 0.8 * q_sat) n_rh80++; }
                    sum_hc += w * H_crit;
                    sum_dq += w * std::max(0.0, c.x[i][j][k] - H_crit * q_sat);
                    sum_cl += w * cloud.x[i][j][k];
                    w_tot  += w;
                    if (cloud.x[i][j][k] > 1e-8) n_cl++;
                    n_tot++;
                }
            }
            const double mdq = (w_tot > 0.0) ? sum_dq / w_tot : 0.0;
            const double mcl = (w_tot > 0.0) ? sum_cl / w_tot : 0.0;
            printf("      AGCM: [CLOUD-INIT] %3d %9.0f   %13.6f  %11.6f  %11.6f  %7.3f  %7.3f  %5.1f %%  %6.2f %%\n",
                   i, get_layer_height(i), cloud_max[i] * 1000.0, mdq * 1000.0, mcl * 1000.0,
                   (w_tot > 0.0) ? sum_rh / w_tot : 0.0, (w_tot > 0.0) ? sum_hc / w_tot : 0.0,
                   100.0 * (double)n_rh80 / (double)std::max(1L, n_tot),
                   100.0 * (double)n_cl / (double)std::max(1L, n_tot));
        }
        // The SAME quantity ATM_CWP_CENSUS reports in MLR, computed here at init, so the two
        // are directly comparable and any loss between init and the first radiation call is
        // localised rather than inferred.
        {
            double cw_sum = 0.0, w_tot2 = 0.0;
            for (int j = 0; j < jm; j++) {
                const double w = cos((j / (double)(jm - 1) - 0.5) * M_PI);
                for (int k2 = 0; k2 < km; k2++) {
                    double col = 0.0;
                    for (int i2 = 0; i2 < im; i2++) {
                        const double dz_i = (i2 < im-1) ? (get_layer_height(i2+1) - get_layer_height(i2))
                                                        : (get_layer_height(i2) - get_layer_height(i2-1));
                        const double T_ii = t.x[i2][j][k2] * t_0;
                        const double rho_ii = (T_ii > 0.0) ? (p_stat.x[i2][j][k2] * 100.0) / (287.0 * T_ii) : 0.0;
                        const double cwl = (cloud.x[i2][j][k2] > 0.0) ? cloud.x[i2][j][k2] : 0.0;
                        const double cwi = (ice.x[i2][j][k2]   > 0.0) ? ice.x[i2][j][k2]   : 0.0;
                        col += (cwl + cwi) * rho_ii * dz_i * 1000.0;
                    }
                    cw_sum += w * col; w_tot2 += w;
                }
            }
            printf("      AGCM: [CLOUD-INIT] COLUMN PATH AT INIT = %.6f g/m2  (compare ATM_CWP_CENSUS in MLR)\n",
                   (w_tot2 > 0.0) ? cw_sum / w_tot2 : 0.0);
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" Time measured: %.3f seconds for initCloudIce\n", elapsed.count() * 1e-9);

    std::cout << "      AGCM: initCloudIce ended" << std::endl;
}
/*
*
*/
