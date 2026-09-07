/*
 * THE COLUMN WATER BUDGET -- which stage supplies the water this model precipitates.
 *
 * WHY IT EXISTS. `printDataAtm` prints, in every run, `P = 1014.6 mm/a  E = 512.5 mm/a`. In any
 * closed atmospheric water budget the long-run global means of precipitation and surface
 * evaporation are EQUAL: this atmosphere holds ~30 mm of precipitable water against ~1000 mm/a
 * of rainfall, so it turns its entire reservoir over in ~11 days and cannot run a persistent
 * imbalance. Measured on a spun-up field (2026-09-06, restart from iteration 600, both parities
 * of the 2dt sawtooth) the gap is 390-500 mm/a and it is FLAT for ninety iterations, so it is
 * not a transient. About 40-50 % of everything this model precipitates has no surface source.
 *
 * The one named candidate has been REFUTED. `ATM_EVAP_SPREAD`'s absolute moisture addition to
 * levels 1..n_spread is a genuine defect read off the source -- a moisture source no flux
 * produced, bounded only by the `c_sat_i` clamp -- but re-running all three of its modes from
 * the same spun-up checkpoint moves the settled evaporation by 0.2 %. The FORM of the moistening
 * is not what supplies the missing ~500 mm/a. So the source is unidentified, and CLAUDE.md's own
 * conclusion was that the next instrument is a column water budget rather than another knob:
 * attribute the change in the column's water to the STAGE that made it, and read the table
 * instead of the control flow.
 *
 * That is the same instrument as `ATM_T0_ATTRIB`, which found that `apply_teq_relaxation` is the
 * only thing that writes ocean level 0 -- after this file had twice recorded a mechanism taken
 * from where a routine is DEFINED instead of where it is CALLED. The pattern earns its place.
 *
 * ------------------------------------------------------------------------------------------
 * WHAT IS MEASURED, AND WHY IT IS NOT SIMPLY THE WATER PATH
 *
 * The quantity is the column TOTAL water -- vapour + cloud + ice + graupel -- because every
 * stage but the microphysics moves water BETWEEN those reservoirs, and a budget on vapour alone
 * would report condensation as a sink. `SaturationAdjustment` conserving total water
 * (`d_cnd + d_dep = d_q_v`) then becomes a CONTROL ROW that must read zero, not a term.
 *
 * The layer mass is FROZEN for the whole reporting window, and that is the difference between
 * an instrument and noise. The obvious formulation, W = SUM(rho*q*dz) differenced between marks,
 * fails: `densities()` rewrites `r_humid` inside every moist iteration, so W moves when nothing
 * about the WATER moved. The window here is ~20 s of physical time and a rate is quoted in mm/a,
 * a factor of 1.6e6 -- so a relative density drift of 1e-6 on a 30 mm reservoir would print as
 * ~470 mm/a, which is the size of the signal being chased. Freezing `rho*dz` at the start of the
 * window makes every stage's contribution SUM(M*dq), the sums telescope exactly, and the
 * identity closes on arithmetic rather than on the density staying put. The true water path is
 * printed beside it, and the difference between the two IS the density drift -- named, and not
 * charged to any stage.
 *
 * The column starts at `i_topography`, not at level 0. Below it the cells are rock; they carry a
 * barometric `p_stat` and receive copied-down surface values from `BC_Atm` Pass 3, and counting
 * them would charge the budget for water that is not in the atmosphere. Over ocean
 * `i_topography == 0` and nothing changes.
 *
 * Each stage is split into the SURFACE BAND (i_top .. i_top+3, exactly the levels
 * `waterVapourEvaporation` writes) and everything above it, because "the manufactured water
 * enters at the surface" and "it enters aloft" are different defects and the total cannot
 * separate them.
 *
 * ------------------------------------------------------------------------------------------
 * THE THREE REFERENCE ROWS, WHICH ARE NOT PART OF THE IDENTITY
 *
 * P and E are the model's own printed means, accumulated over the same window with the same
 * `GetMean_3D`/`GetMean_2D` calls `printDataAtm` uses, so the numbers in this table are the
 * numbers in that one. They are reference rows because precipitation LEAVES the domain at the
 * ground: it is not a term in a column budget, it is what the budget has to be consistent with.
 *
 * The third row exists because of where the microphysics actually acts. The ice schemes do NOT
 * write `c`/`cloud`/`ice` -- they write the rate arrays `S_v`/`S_c`/`S_i`/`S_g`, which reach the
 * water only through `rhs_c`/`rhs_cloud`/`rhs_ice`/`rhs_g` in the RK4 (`RHS_Atm_Turb.cpp:1296`,
 * `coeff_trans * S_x * r_humid`). So the ice scheme's own bucket is near zero and its sink
 * appears inside the `RungeKutta` bucket, mixed with transport. The reference row is the exact
 * term RK4 adds -- `(S_v + S_c + S_i + S_g) * r_humid * dt`, integrated with the same frozen
 * mass -- and it is exact rather than approximate because the S-terms are constant across the
 * four stages, so the RK4 weights (k1 + 2k2 + 2k3 + k4)/6 return S*dt identically.
 *
 *     RungeKutta bucket - microphysics reference = transport non-conservation
 *                                                  + the max(0,...) floor in every RK4 stage
 *
 * and that floor (`RungeKutta_Atm_Turb.cpp:188` and its siblings, one per stage per species) is
 * a water SOURCE of exactly the shape that manufactured 8129 mm/a in the microphysics before
 * 2026-09-01. This instrument cannot split those two, and says so rather than implying it can.
 *
 * ------------------------------------------------------------------------------------------
 * ATM_CWB_DIAG=1 enables. Default off. It never writes a field -- every hook reads -- and the
 * two 3-D scratch arrays (~21 MB each on a 41x181x361 grid) are allocated only when it is on.
 * The state is a function-local static rather than a member of `cAtmosphereModel`, deliberately:
 * adding a member moves `sizeof(cAtmosphereModel)` and that is this tree's stack-canary hazard.
 * Sums are accumulated per latitude row and added serially, so the table does not depend on the
 * thread count -- this tree's own caveat about OpenMP reduction order applies to a diagnostic
 * whose whole claim is that a residual is zero.
 */

#pragma once

#include "cAtmosphereModel.h"
#include "Utils.h"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

class ColumnWaterBudget {
public:
    static bool enabled(){
        static const bool v = [](){ const char* e = getenv("ATM_CWB_DIAG");
                                    return e && atoi(e) != 0; }();
        return v;
    }

    // Difference the column water against the previous mark and charge it to `stage`.
    static void mark(cAtmosphereModel& m, const char* stage){
        if(enabled()) state().mark(m, stage);
    }
    // Once per iteration, before the RK4 that consumes the S-terms: accumulate the window's
    // elapsed physical time and the reference rows. `iter` is passed rather than read off the
    // model because `cAtmosphereModel::iter_n` is private.
    static void tick(cAtmosphereModel& m, int iter, double sec_per_iter){
        if(enabled()) state().tick(m, iter, sec_per_iter);
    }
    // Print the window and open the next one.
    static void report(cAtmosphereModel& m, int iter){
        if(enabled()) state().report(m, iter);
    }

private:
    struct Impl {
        std::vector<double>      mass;      // rho*dz, frozen for the window          [kg/m2]
        std::vector<double>      q_prev;    // total water at the previous mark       [kg/kg]
        std::vector<std::string> names;     // stage names, in first-seen order
        std::vector<double>      sums;      // 2 per stage: surface band, aloft       [mm]
        double w_lat   = 0.0;               // sum of the cos-lat weights over (j,k)
        double elapsed = 0.0;               // physical seconds in this window
        double P_mm = 0.0, E_mm = 0.0;
        double Sq_mm = 0.0, Sq_nd = 0.0;   // S_v+S_c+S_i+S_g, as shipped / at L/u_0
        double Sp_mm = 0.0, Sp_nd = 0.0;   // S_r+S_s,         as shipped / at L/u_0
        double W_start = 0.0;               // true water path at the window start    [mm]
        int    it_start = 0;
        bool   it_valid = false;
        bool   open = false;

        // The stages, in the order the time loop executes them, so the table reads in execution
        // order however the moist stride falls. A stage that did not run this window prints
        // zero, which is information: `IceScheme` reading 0.0 is how one sees that the ice
        // schemes write the S-rate arrays and not the water.
        static const std::vector<std::string>& loop_order(){
            static const std::vector<std::string> v = {
                "pressure+project", "SaturationAdjust", "damp_wiggles(q)", "IceScheme",
                "MoistConvection", "cap_S+clamp", "ConvectiveAdjust", "ThermoAtm(pre)",
                "evaporation", "ThermoAtm(rest)", "BC_Atm", "RungeKutta",
                "orographic_shapiro", "radiation/teq", "unattributed" };
            return v;
        }

        // The weight `GetMean_2D`/`GetMean_3D` use, replicated rather than called so that the
        // means in this table are on exactly the weights `Precip mean` and the NASA score are
        // on. The hard-coded 90 is theirs (lib/Utils.cpp:37), not a latitude computed from jm.
        static double lat_weight(int j){
            return (j <= 90) ? cos((90 - j) * M_PI / 180.0)
                             : cos((j - 90) * M_PI / 180.0);
        }

        static double q_total(cAtmosphereModel& m, int i, int j, int k){
            return m.c.x[i][j][k] + m.cloud.x[i][j][k] + m.ice.x[i][j][k] + m.gr.x[i][j][k];
        }

        // The true column water path, cos-lat mean [mm]. Not the budget quantity -- it moves
        // with the density -- but it is the reservoir the rates have to be read against.
        double water_path(cAtmosphereModel& m) const {   // non-const: get_layer_height() is
            std::vector<double> row(m.jm, 0.0);
            #pragma omp parallel for schedule(static)
            for(int j = 0; j < m.jm; j++){
                double s = 0.0;
                const double wj = lat_weight(j);
                for(int k = 0; k < m.km; k++){
                    const int i0 = m.i_topography[j][k];
                    for(int i = i0; i < m.im - 1; i++){
                        double rho = m.r_humid.x[i][j][k];
                        if(!AtomUtils::is_finite_safe(rho) || rho <= 0.0) rho = m.r_air;
                        s += wj * rho * (m.get_layer_height(i+1) - m.get_layer_height(i))
                                * q_total(m, i, j, k);
                    }
                }
                row[j] = s;
            }
            double tot = 0.0;
            for(int j = 0; j < m.jm; j++) tot += row[j];
            return (w_lat > 0.0) ? tot / w_lat : 0.0;
        }

        size_t idx(cAtmosphereModel& m, int i, int j, int k) const {
            return ((size_t)i * m.jm + j) * m.km + k;
        }

        void begin_window(cAtmosphereModel& m, int iter){
            const size_t n = (size_t)m.im * m.jm * m.km;
            if(mass.size()   != n) mass.assign(n, 0.0);
            if(q_prev.size() != n) q_prev.assign(n, 0.0);
            if(w_lat <= 0.0){
                for(int j = 0; j < m.jm; j++) w_lat += lat_weight(j) * m.km;
            }
            #pragma omp parallel for collapse(2) schedule(static)
            for(int j = 0; j < m.jm; j++){
                for(int k = 0; k < m.km; k++){
                    const int i0 = m.i_topography[j][k];
                    for(int i = 0; i < m.im; i++){
                        const size_t p = idx(m, i, j, k);
                        q_prev[p] = q_total(m, i, j, k);
                        if(i < i0 || i >= m.im - 1){ mass[p] = 0.0; continue; }
                        double rho = m.r_humid.x[i][j][k];
                        if(!AtomUtils::is_finite_safe(rho) || rho <= 0.0) rho = m.r_air;
                        mass[p] = rho * (m.get_layer_height(i+1) - m.get_layer_height(i));
                    }
                }
            }
            names = loop_order();
            sums.assign(2 * names.size(), 0.0);
            elapsed = 0.0;
            P_mm = E_mm = 0.0;
            Sq_mm = Sq_nd = Sp_mm = Sp_nd = 0.0;
            W_start  = water_path(m);
            it_start = iter;
            open     = true;
        }

        void mark(cAtmosphereModel& m, const char* stage){
            if(!open){ begin_window(m, 0); it_valid = false; return; }

            // Per-latitude partials, summed serially, so the table is thread-count independent.
            std::vector<double> row_s(m.jm, 0.0), row_a(m.jm, 0.0);
            #pragma omp parallel for schedule(static)
            for(int j = 0; j < m.jm; j++){
                double ss = 0.0, sa = 0.0;
                const double wj = lat_weight(j);
                for(int k = 0; k < m.km; k++){
                    const int i0 = m.i_topography[j][k];
                    for(int i = i0; i < m.im - 1; i++){
                        const size_t p = idx(m, i, j, k);
                        const double q = q_total(m, i, j, k);
                        const double d = wj * mass[p] * (q - q_prev[p]);
                        if(i <= i0 + 3) ss += d; else sa += d;
                        q_prev[p] = q;
                    }
                }
                row_s[j] = ss;
                row_a[j] = sa;
            }
            double ds = 0.0, da = 0.0;
            for(int j = 0; j < m.jm; j++){ ds += row_s[j]; da += row_a[j]; }

            size_t s = 0;
            while(s < names.size() && names[s] != stage) s++;
            if(s == names.size()){ names.emplace_back(stage); sums.resize(2 * (s + 1), 0.0); }
            sums[2*s]     += ds;
            sums[2*s + 1] += da;
        }

        void tick(cAtmosphereModel& m, int iter, double sec_per_iter){
            if(!open) return;
            if(!it_valid){ it_start = iter - 1; it_valid = true; }
            elapsed += sec_per_iter;

            // P and E through the model's own means, so these rows ARE the printed ones.
            P_mm += AtomUtils::GetMean_3D(m.jm, m.km, m.Precipitation) * sec_per_iter; // [mm/s]
            E_mm += AtomUtils::GetMean_2D(m.jm, m.km, m.Evaporation) / 8.64e4 * sec_per_iter;

            // The microphysics term exactly as RK4 adds it: coeff_trans = 1.0 at
            // RHS_Atm_Turb.cpp:372, applied over the non-dimensional step m.dt.
            //
            // AND THE SAME TERM UNDER THE NEIGHBOURING NON-DIMENSIONALISATION. `S_v` is a rate
            // in kg/(kg*s) -- TwoCatIceScheme.h:571 says so on the assignment -- so reaching a
            // tendency in non-dimensional time it needs the factor L/u_0, which is what
            // `coeff_MC_q = ndimLength()/(u_0*c_0)` gives the convective moisture source two
            // lines below it (RHS_Atm_Turb.cpp:374, 1298). The S-terms instead carry `r_humid`,
            // a DENSITY. Both numbers are printed and no claim is made here about which is
            // right: the ratio is (L/u_0)/rho, ~1670 at the surface and larger aloft, and a
            // number that large belongs in a measurement rather than in a comment.
            // L/u_0 without reaching for the private metricShellLength(): the caller
            // already passes dt*L/u_0 as the seconds per iteration.
            const double L_over_u0 = (m.dt != 0.0) ? sec_per_iter / m.dt : 0.0;
            // The coefficient RK4 ACTUALLY uses, which is not `r_humid` once ATM_MICRO_NDIM is
            // set -- reading the shipped one would make this row constant across a sweep on the
            // very coefficient it is meant to measure. Mirrors RHS_Atm_Turb.cpp exactly.
            static const double micro_s = [](){ const char* e = getenv("ATM_MICRO_NDIM");
                                                return e ? atof(e) : 0.0; }();
            //
            // AND THE RAIN/SNOW RATES BESIDE THEM, AS A CONSERVATION CHECK ON THE SCHEME'S OWN
            // ARRAYS. A scheme that conserves has (S_v + S_c + S_i + S_g) + (S_r + S_s) = 0 at
            // every cell, so the two rows sum to zero or they do not, and that is printed
            // rather than assumed.
            //
            // IT IS NOT THE GROUND FLUX AND MUST NOT BE READ AS ONE. The flux is integrated
            // from a DIFFERENT SUBSET of the same terms: `dP_rain` omits `S_r_cri`
            // (TwoCatIceScheme.h:613, whose own comment says "in S_r, NOT in dP_rain") and
            // `dP_snow` omits `S_s_dep`, `S_i_cri` and `S_r_cri`, each for a stated
            // double-counting reason. So SUM(rho*(S_r+S_s)*dz) and P(ground) are different
            // quantities and can differ in SIGN, which they do here. What the two rows are for
            // is the pair of magnitudes: the microphysical exchange as RK4 APPLIES it, against
            // the same rates read as the kg/(kg*s) their assignment says they are.
            std::vector<double> row(m.jm, 0.0), row_nd(m.jm, 0.0),
                                rowp(m.jm, 0.0), rowp_nd(m.jm, 0.0);
            #pragma omp parallel for schedule(static)
            for(int j = 0; j < m.jm; j++){
                double s = 0.0, s_nd = 0.0, sp = 0.0, sp_nd = 0.0;
                const double wj = lat_weight(j);
                for(int k = 0; k < m.km; k++){
                    const int i0 = m.i_topography[j][k];
                    for(int i = i0; i < m.im - 1; i++){
                        const double S = m.S_v.x[i][j][k] + m.S_c.x[i][j][k]
                                       + m.S_i.x[i][j][k] + m.S_g.x[i][j][k];
                        const double R = m.S_r.x[i][j][k] + m.S_s.x[i][j][k];
                        const double M = wj * mass[idx(m, i, j, k)];
                        double cm = m.r_humid.x[i][j][k];
                        if(micro_s != 0.0) cm += micro_s * (L_over_u0 - cm);
                        s     += M * S * cm         * m.dt;   // as RK4 applies it
                        s_nd  += M * S * L_over_u0  * m.dt;   // the same rates at L/u_0
                        sp    += M * R * cm         * m.dt;
                        sp_nd += M * R * L_over_u0  * m.dt;
                    }
                }
                row[j]     = s;
                row_nd[j]  = s_nd;
                rowp[j]    = sp;
                rowp_nd[j] = sp_nd;
            }
            double tot = 0.0, tot_nd = 0.0, totp = 0.0, totp_nd = 0.0;
            for(int j = 0; j < m.jm; j++){
                tot  += row[j];  tot_nd  += row_nd[j];
                totp += rowp[j]; totp_nd += rowp_nd[j];
            }
            const double iw = (w_lat > 0.0) ? 1.0 / w_lat : 0.0;
            Sq_mm += tot * iw;  Sq_nd += tot_nd * iw;
            Sp_mm += totp * iw; Sp_nd += totp_nd * iw;
        }

        void report(cAtmosphereModel& m, int iter){
            using namespace std;
            if(!open || names.empty() || elapsed <= 0.0) return;

            const double per_year = 365.0 * 8.64e4 / elapsed;   // [mm] over the window -> [mm/a]
            const double inv_w    = (w_lat > 0.0) ? 1.0 / w_lat : 0.0;

            cout << "      AGCM: [CWB] column water budget, iterations " << it_start + 1
                 << ".." << iter << "  (" << fixed << setprecision(2) << elapsed
                 << " s of physical time; layer mass frozen at iteration " << it_start
                 << "; total water = vapour + cloud + ice + graupel, from i_topography up)" << endl;
            cout << "      AGCM: [CWB] " << left << setw(22) << "stage" << right
                 << setw(15) << "total mm/a" << setw(15) << "sfc i0..i0+3"
                 << setw(15) << "aloft" << endl;

            double tot_s = 0.0, tot_a = 0.0;
            for(size_t s = 0; s < names.size(); s++){
                const double a = sums[2*s] * inv_w * per_year;
                const double b = sums[2*s + 1] * inv_w * per_year;
                cout << "      AGCM: [CWB] " << left << setw(22) << names[s] << right
                     << setw(15) << scientific << setprecision(4) << (a + b)
                     << setw(15) << a << setw(15) << b << endl;
                tot_s += sums[2*s];
                tot_a += sums[2*s + 1];
            }
            const double net_s = tot_s * inv_w * per_year;
            const double net_a = tot_a * inv_w * per_year;
            cout << "      AGCM: [CWB] " << left << setw(22) << "NET" << right
                 << setw(15) << scientific << setprecision(4) << (net_s + net_a)
                 << setw(15) << net_s << setw(15) << net_a << endl;

            // Reference rows. P leaves the domain at the ground, so it is not a budget term; the
            // microphysics row is the part of the RungeKutta bucket that the ice scheme asked for.
            cout << "      AGCM: [CWB] reference (NOT part of the identity):"
                 << "  P = " << fixed << setprecision(1) << P_mm * per_year
                 << "   E = " << E_mm * per_year
                 << "   P - E = " << (P_mm - E_mm) * per_year
                 << " mm/a" << endl;
            static const double micro_ndim = [](){ const char* e = getenv("ATM_MICRO_NDIM");
                                                   return e ? atof(e) : 0.0; }();
            cout << "      AGCM: [CWB] microphysics rate arrays over the column [mm/a], NOT the"
                 << " ground flux (different subset -- see the header);"
                 << "  ATM_MICRO_NDIM = " << scientific << setprecision(3)
                 << micro_ndim << fixed << endl;
            cout << "      AGCM: [CWB]     as RK4 applies them (S*coeff):    vap+cld+ice+grp "
                 << setprecision(2) << Sq_mm * per_year << "   rain+snow " << Sp_mm * per_year
                 << "   sum " << (Sq_mm + Sp_mm) * per_year
                 << "  <- conservation check" << endl;
            cout << "      AGCM: [CWB]     the same rates at L/u_0:          vap+cld+ice+grp "
                 << setprecision(1) << Sq_nd * per_year << "   rain+snow " << Sp_nd * per_year
                 << "   sum " << (Sq_nd + Sp_nd) * per_year << endl;

            // The reservoir, and the part of its change that is density rather than water.
            const double W_now = water_path(m);
            cout << "      AGCM: [CWB] reservoir: total water path " << setprecision(4) << W_now
                 << " mm now, " << W_start << " mm at the window start;"
                 << "  budget says " << scientific << setprecision(3)
                 << (net_s + net_a) / (365.0 * 8.64e4) * elapsed
                 << " mm, the rest (" << ((W_now - W_start)
                        - (net_s + net_a) / (365.0 * 8.64e4) * elapsed)
                 << " mm) is the density moving, not water" << endl;

            cout << defaultfloat;
            begin_window(m, iter);
        }
    };

    static Impl& state(){ static Impl s; return s; }
};
