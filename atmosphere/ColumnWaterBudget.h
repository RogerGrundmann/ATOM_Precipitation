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
 * 2026-09-01.
 *
 * ------------------------------------------------------------------------------------------
 * SPLITTING THE RUNGEKUTTA BUCKET (2026-09-22) -- the paragraph above used to end "this
 * instrument cannot split those two, and says so rather than implying it can." It can now.
 *
 * The bucket is the largest row in the table (-7.8e+06 mm/a against a NET of +1.6e+05) and was
 * the only unattributed one, which matters because the whole column water is the ~2 % residue
 * of it and the +5.0e+06 evaporation term: a 0.01 % asymmetry between them is 500 mm/a, the
 * size of the P - E gap this file was written to chase. The split is
 *
 *     bucket = INT M*(cn - c)   [leapfrog_reset, measured at RK4 entry]
 *            + INT M*D          [the RK4 tendency: microphysics + transport + MC_q + diffusion]
 *            + INT M*floor      [the FINAL-stage max(0,...) clips, and only those]
 *
 * because RK4 writes `x = max(0, xn + D)` where `xn` is the time-level-n array.
 *
 * TWO THINGS IN IT ARE NOT OBVIOUS FROM THE ARITHMETIC.
 *
 * (1) `leapfrog_reset` is not a rounding term. `storeIntermediateData3D` writes `cn := c` at
 * `cAtmosphereModel.cpp:1940`, which is AFTER the RK4 call at :1762 -- so during RK4 at
 * iteration N, `cn` holds the state saved at the END of iteration N-1, while `c` also carries
 * everything `SaturationAdjustment`, the ice scheme, `MoistConvection`, `ConvectiveAdjustment`,
 * `waterVapourEvaporation` and `BC_Atm` wrote into it earlier in iteration N. RK4 integrates
 * from `cn` and overwrites `c`, so those direct state writes are discarded and survive only
 * through their effect on `rhs_*` and on the other prognostic fields. Whether that is a defect
 * or the intended leapfrog structure is NOT decided here -- the term is named, measured and
 * printed, and the arithmetic says how large it is. The pre-registered expectation, from the
 * published table alone, is that it is LARGE and negative: evaporation +5.0e+06, damp_wiggles(q)
 * +2.5e+06 and SaturationAdjust +5.5e+05 sum to ~+8.0e+06 against a bucket of -7.8e+06, and
 * that near-cancellation is exactly what a reset of the pre-RK4 increments would produce. If it
 * comes out SMALL instead, the reading in this paragraph is wrong and the cancellation is
 * something else.
 *
 * (2) ONLY THE FINAL-STAGE CLIPS ARE IN THE IDENTITY. The k1/k2/k3 clips fire on the
 * intermediate stage values, which the NEXT stage's RHS reads (`c.x` is the live array all four
 * stages evaluate from), so they perturb the integrator rather than adding to the answer: their
 * effect on the final field is not their own magnitude, and no sum of them belongs in the
 * budget. They are printed as a reference row with their own clip counts, because "the floor is
 * a positivity guard that fires rarely at sharp coastal gradients" is a claim
 * `RungeKutta_Atm_Turb.cpp:185` makes and nothing has ever measured -- and in this tree that
 * claim has been orders out three times (`P_max_flux`, the 8129 mm/a microphysics floor,
 * `SaturationAdjustment`'s phase-split clip).
 *
 * `rest` is a REMAINDER and is printed with a control: transport is a divergence and integrates
 * to zero over the sphere, so |rest| should be a small fraction of the largest named term. If
 * it is not, the split has not accounted for the bucket and the table says so instead of
 * leaving a large residual looking like a result.
 *
 * No second environment knob. The floor and the reset are part of THIS budget and are
 * meaningless outside it, and a half-on state is the failure mode this repository keeps
 * rediscovering (`warnIfHalfRepaired`). `ATM_CWB_DIAG=1` turns on all of it; unset, every hook
 * is a branch on a `static const bool` and nothing is allocated.
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

#include <algorithm>
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
    static bool bands(){
        static const bool v = [](){ const char* e = getenv("ATM_CWB_BANDS"); return e && atoi(e) != 0; }();
        return v;
    }
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

    // ---------------------------------------------------------------------------------------
    // THE RUNGEKUTTA BUCKET, SPLIT. The two hooks below exist because that bucket is the
    // largest row in the table and was the only one this instrument could not attribute -- the
    // header above says so in as many words. See "SPLITTING THE RUNGEKUTTA BUCKET" there.

    // Called from RungeKutta_Atm_Turb.cpp at each of the 16 `std::max(0.0, ...)` clip sites,
    // ONLY when the pre-clip value is negative. `sp` 0..3 = vapour/cloud/ice/graupel,
    // `st` 0..3 = k1/k2/k3/final, `inj` = the (positive) mixing ratio the clip adds.
    static void floor_hit(cAtmosphereModel& m, int sp, int st, int i, int j, int k, double inj){
        state().floor_hit(m, sp, st, i, j, k, inj);
    }
    // Called immediately before solveRungeKutta_Atmosphere_Turb(). Measures INT M*(cn - c),
    // the pre-RK4 stage increments the integrator is about to discard by starting from `cn`.
    static void mark_leapfrog(cAtmosphereModel& m){
        if(enabled()) state().mark_leapfrog(m);
    }

private:
    struct Impl {
        std::vector<double>      mass;      // rho*dz, frozen for the window          [kg/m2]
        std::vector<double>      q_prev;    // total water at the previous mark       [kg/kg]
        std::vector<std::string> names;     // stage names, in first-seen order
        std::vector<double>      sums;      // 2 per stage: surface band, aloft       [mm]
        // ATM_CWB_BANDS=1 (2026-09-25, print-only): the same stage totals split by |latitude| into
        // 0-15 / 15-35 / 35-65 / 65-90, the bands the precipitation score prints (ThermoAtm.h,
        // alat = |90 - j*180/(jm-1)|, edges < 15 / < 35 / < 65). Asked for by the closure's
        // polar engine: with ATM_WATER_CLOSURE on and both filter vertical passes off, 65-90
        // still climbs 28 -> 669 mm/a over 100 iterations (qvt_on) and the global table cannot
        // say which stage feeds it. Each band is normalised by its OWN cos-lat weight, so a row
        // reads as a band-mean rate in mm/a, comparable with that band's printed precipitation.
        std::vector<double>      sums_b;    // 4 per stage: bands 0-15, 15-35, 35-65, 65-90   [mm*wj]
        double w_band[4] = {0.0, 0.0, 0.0, 0.0};
        double w_lat   = 0.0;               // sum of the cos-lat weights over (j,k)
        double elapsed = 0.0;               // physical seconds in this window
        double P_mm = 0.0, E_mm = 0.0;
        double Sq_mm = 0.0, Sq_nd = 0.0;   // S_v+S_c+S_i+S_g, as shipped / at L/u_0
        double Sp_mm = 0.0, Sp_nd = 0.0;   // S_r+S_s,         as shipped / at L/u_0
        double W_start = 0.0;               // true water path at the window start    [mm]
        int    it_start = 0;
        bool   it_valid = false;
        bool   open = false;

        // The RungeKutta split. `flr_buf` is indexed ((st*4 + sp)*im + i)*jm + j so that each
        // (i,j) belongs to exactly one thread under the RK4 loop's `collapse(2)` over (i,j) --
        // the same reason `mark()` accumulates per latitude row: a diagnostic whose claim is
        // that a residual is zero must not depend on the thread count. Reduced serially at
        // report time. Counts are kept separately and are exactly reproducible; the masses are
        // sums of doubles in a fixed order and so are too.
        std::vector<double>    flr_buf;    // injected water, per species/stage    [kg/m2 * wj]
        std::vector<long long> flr_cnt;    // clip events,    per species/stage
        double lf_s = 0.0, lf_a = 0.0;     // INT M*(cn - c) at RK4 entry, sfc / aloft  [mm*wj]
        int    rk_calls = 0;               // RK4 sweeps in this window

        // The stages, in the order the time loop executes them, so the table reads in execution
        // order however the moist stride falls. A stage that did not run this window prints
        // zero, which is information: `IceScheme` reading 0.0 is how one sees that the ice
        // schemes write the S-rate arrays and not the water.
        static const std::vector<std::string>& loop_order(){
            static const std::vector<std::string> v = {
                "pressure+project", "SaturationAdjust", "damp_wiggles(q)", "IceScheme",
                "MoistConvection", "cap_S+clamp", "ConvectiveAdjust", "ThermoAtm(pre)",
                "evaporation", "ThermoAtm(rest)", "BC:radius", "BC:theta", "BC:phi(seam)",
                "BC:scalarSurfSur", "BC:solidGround", "RungeKutta",
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

        static int band_of(cAtmosphereModel& m, int j){
            const double alat = fabs(90.0 - j * 180.0 / (double)(m.jm - 1));
            return (alat < 15.0) ? 0 : (alat < 35.0) ? 1 : (alat < 65.0) ? 2 : 3;
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
                for(int k = 0; k < m.km - 1; k++){   // k = km-1 IS k = 0 (the seam): count it once
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
                for(int j = 0; j < m.jm; j++) w_lat += lat_weight(j) * (m.km - 1);   // seam counted once (2026-09-26)
            }
            #pragma omp parallel for collapse(2) schedule(static)
            for(int j = 0; j < m.jm; j++){
                for(int k = 0; k < m.km - 1; k++){   // k = km-1 IS k = 0 (the seam): count it once
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
            if(bands()){
                sums_b.assign(4 * names.size(), 0.0);
                if(w_band[0] + w_band[1] + w_band[2] + w_band[3] <= 0.0)
                    for(int j = 0; j < m.jm; j++) w_band[band_of(m, j)] += lat_weight(j) * (m.km - 1);   // seam once
            }
            elapsed = 0.0;
            P_mm = E_mm = 0.0;
            Sq_mm = Sq_nd = Sp_mm = Sp_nd = 0.0;
            W_start  = water_path(m);
            it_start = iter;
            open     = true;

            const size_t nf = 16 * (size_t)m.im * m.jm;
            if(flr_buf.size() != nf){ flr_buf.assign(nf, 0.0); flr_cnt.assign(16, 0); }
            else { std::fill(flr_buf.begin(), flr_buf.end(), 0.0);
                   std::fill(flr_cnt.begin(), flr_cnt.end(), 0LL); }
            lf_s = lf_a = 0.0;
            rk_calls = 0;
        }

        // ---- the RungeKutta split -------------------------------------------------------
        size_t fidx(cAtmosphereModel& m, int sp, int st, int i, int j) const {
            return (((size_t)st * 4 + sp) * m.im + i) * m.jm + j;
        }

        void floor_hit(cAtmosphereModel& m, int sp, int st, int i, int j, int k, double inj){
            if(!open || mass.empty()) return;
            const double M = lat_weight(j) * mass[idx(m, i, j, k)];   // 0 in rock, by construction
            if(M == 0.0) return;
            flr_buf[fidx(m, sp, st, i, j)] += M * inj;
            #pragma omp atomic
            flr_cnt[st * 4 + sp] += 1;
        }

        // INT M*(cn - c) at RK4 entry. RK4 integrates from `cn`, which `storeIntermediateData3D`
        // last wrote at the END of the previous iteration (cAtmosphereModel.cpp:1940, i.e. AFTER
        // this call site), so everything the pre-RK4 stages wrote into `c` this iteration is
        // about to be overwritten. Whether that is a defect or the intended leapfrog structure
        // is NOT decided here: the term is measured, named and printed, and the table says which
        // of the two the arithmetic supports. Read from `cn`/`c` directly rather than from
        // `q_prev` so the number is exact whatever happened since the last mark.
        void mark_leapfrog(cAtmosphereModel& m){
            if(!open || mass.empty()) return;
            rk_calls++;
            std::vector<double> row_s(m.jm, 0.0), row_a(m.jm, 0.0);
            #pragma omp parallel for schedule(static)
            for(int j = 0; j < m.jm; j++){
                double ss = 0.0, sa = 0.0;
                const double wj = lat_weight(j);
                for(int k = 0; k < m.km - 1; k++){   // k = km-1 IS k = 0 (the seam): count it once
                    const int i0 = m.i_topography[j][k];
                    for(int i = i0; i < m.im - 1; i++){
                        const double qn = m.cn.x[i][j][k] + m.cloudn.x[i][j][k]
                                        + m.icen.x[i][j][k] + m.grn.x[i][j][k];
                        const double d  = wj * mass[idx(m, i, j, k)] * (qn - q_total(m, i, j, k));
                        if(i <= i0 + 3) ss += d; else sa += d;
                    }
                }
                row_s[j] = ss;
                row_a[j] = sa;
            }
            for(int j = 0; j < m.jm; j++){ lf_s += row_s[j]; lf_a += row_a[j]; }
        }

        void mark(cAtmosphereModel& m, const char* stage){
            if(!open){ begin_window(m, 0); it_valid = false; return; }

            // Per-latitude partials, summed serially, so the table is thread-count independent.
            std::vector<double> row_s(m.jm, 0.0), row_a(m.jm, 0.0);
            #pragma omp parallel for schedule(static)
            for(int j = 0; j < m.jm; j++){
                double ss = 0.0, sa = 0.0;
                const double wj = lat_weight(j);
                for(int k = 0; k < m.km - 1; k++){   // k = km-1 IS k = 0 (the seam): count it once
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
            if(s == names.size()){ names.emplace_back(stage); sums.resize(2 * (s + 1), 0.0);
                                   if(bands()) sums_b.resize(4 * (s + 1), 0.0); }
            sums[2*s]     += ds;
            sums[2*s + 1] += da;
            if(bands())                                   // serial over j: thread-count independent
                for(int j = 0; j < m.jm; j++) sums_b[4*s + band_of(m, j)] += row_s[j] + row_a[j];
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
            // Default 1.0 since 2026-09-21 -- MUST track RHS_Atm_Turb.cpp's default exactly,
            // or this row reports the shipped coefficient while RK4 applies the corrected one.
            static const double micro_s = [](){ const char* e = getenv("ATM_MICRO_NDIM");
                                                return e ? atof(e) : 1.0; }();
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
                for(int k = 0; k < m.km - 1; k++){   // k = km-1 IS k = 0 (the seam): count it once
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
            if(bands() && sums_b.size() >= 4 * names.size()){
                cout << "      AGCM: [CWB-BANDS] band-mean rates [mm/a], each band on its own cos-lat weight"
                     << " (compare the band's printed precipitation)" << endl;
                cout << "      AGCM: [CWB-BANDS] " << left << setw(22) << "stage" << right
                     << setw(13) << "0-15" << setw(13) << "15-35" << setw(13) << "35-65"
                     << setw(13) << "65-90" << endl;
                double nb[4] = {0.0, 0.0, 0.0, 0.0};
                for(size_t s = 0; s < names.size(); s++){
                    cout << "      AGCM: [CWB-BANDS] " << left << setw(22) << names[s] << right;
                    for(int b = 0; b < 4; b++){
                        const double v = (w_band[b] > 0.0) ? sums_b[4*s + b] / w_band[b] * per_year : 0.0;
                        nb[b] += sums_b[4*s + b];
                        cout << setw(13) << scientific << setprecision(3) << v;
                    }
                    cout << endl;
                }
                cout << "      AGCM: [CWB-BANDS] " << left << setw(22) << "NET" << right;
                for(int b = 0; b < 4; b++)
                    cout << setw(13) << scientific << setprecision(3)
                         << ((w_band[b] > 0.0) ? nb[b] / w_band[b] * per_year : 0.0);
                cout << endl << fixed;
            }

            // Reference rows. P leaves the domain at the ground, so it is not a budget term; the
            // microphysics row is the part of the RungeKutta bucket that the ice scheme asked for.
            cout << "      AGCM: [CWB] reference (NOT part of the identity):"
                 << "  P = " << fixed << setprecision(1) << P_mm * per_year
                 << "   E = " << E_mm * per_year
                 << "   P - E = " << (P_mm - E_mm) * per_year
                 << " mm/a" << endl;
            static const double micro_ndim = [](){ const char* e = getenv("ATM_MICRO_NDIM");
                                                   return e ? atof(e) : 1.0; }();
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

            // ---- THE RUNGEKUTTA BUCKET, SPLIT -------------------------------------------
            // bucket = INT M*[ max(0, cn + D) - c_at_previous_mark ]
            //        = INT M*(cn - c)            <- leapfrog_reset, measured at RK4 entry
            //        + INT M*D                   <- the RK4 tendency: microphysics + the rest
            //        + INT M*(floor, FINAL stage) <- the only clip that writes the field mark() sees
            //
            // The k1/k2/k3 clips are NOT in that identity and are printed as a reference row.
            // They fire on the INTERMEDIATE stage values, which the next stage's RHS reads, so
            // they perturb the integrator rather than adding to the answer: their effect on the
            // final field is not their own magnitude and this instrument cannot price it.
            if(rk_calls > 0 && !flr_buf.empty()){
                const char* spn[4] = {"vapour", "cloud", "ice", "graupel"};
                const char* stn[4] = {"k1", "k2", "k3", "final"};
                double fs[4][4];
                for(int st = 0; st < 4; st++)
                    for(int sp = 0; sp < 4; sp++){
                        double acc = 0.0;
                        for(int i = 0; i < m.im; i++)
                            for(int j = 0; j < m.jm; j++)
                                acc += flr_buf[fidx(m, sp, st, i, j)];
                        fs[st][sp] = acc * inv_w * per_year;
                    }

                double f_fin = 0.0, f_int = 0.0;
                long long n_fin = 0, n_int = 0;
                for(int sp = 0; sp < 4; sp++){
                    f_fin += fs[3][sp];              n_fin += flr_cnt[3*4 + sp];
                    for(int st = 0; st < 3; st++){ f_int += fs[st][sp]; n_int += flr_cnt[st*4 + sp]; }
                }

                size_t rk = 0;
                while(rk < names.size() && names[rk] != "RungeKutta") rk++;
                const double bucket = (rk < names.size())
                                    ? (sums[2*rk] + sums[2*rk + 1]) * inv_w * per_year : 0.0;
                const double lf    = (lf_s + lf_a) * inv_w * per_year;
                const double micro = Sq_mm * per_year;
                const double rest  = bucket - lf - micro - f_fin;

                cout << "      AGCM: [CWB] RK4 positivity floor -- water INJECTED by"
                     << " max(0,...) [mm/a], " << rk_calls << " RK4 sweeps:" << endl;
                cout << "      AGCM: [CWB]     " << left << setw(10) << "stage" << right;
                for(int sp = 0; sp < 4; sp++) cout << setw(14) << spn[sp];
                cout << setw(14) << "TOTAL" << setw(14) << "clips/call" << endl;
                for(int st = 0; st < 4; st++){
                    double tot = 0.0; long long cn_ = 0;
                    for(int sp = 0; sp < 4; sp++){ tot += fs[st][sp]; cn_ += flr_cnt[st*4 + sp]; }
                    cout << "      AGCM: [CWB]     " << left << setw(10) << stn[st] << right
                         << scientific << setprecision(4);
                    for(int sp = 0; sp < 4; sp++) cout << setw(14) << fs[st][sp];
                    cout << setw(14) << tot << fixed << setprecision(1)
                         << setw(14) << (double)cn_ / rk_calls << endl;
                }
                cout << "      AGCM: [CWB]     k1+k2+k3 (NOT part of the identity -- they"
                     << " perturb the integrator, not the answer): " << scientific
                     << setprecision(4) << f_int << " mm/a, " << fixed << setprecision(1)
                     << (double)n_int / rk_calls << " clips/call" << endl;

                cout << "      AGCM: [CWB] RungeKutta bucket " << scientific << setprecision(4)
                     << bucket << " = leapfrog_reset " << lf << " + microphysics " << micro
                     << " + floor(final) " << f_fin << " + rest " << rest << " mm/a" << endl;
                cout << "      AGCM: [CWB]     leapfrog_reset = INT M*(cn - c) at RK4 entry: the"
                     << " pre-RK4 stage increments the integrator starts over the top of."
                     << "  rest = transport + MC_q + diffusion, a divergence-dominated term"
                     << " whose global mean should be SMALL against the others." << endl;
                // `rest` is a REMAINDER, not a measurement, so it needs a control. Transport
                // is a divergence and integrates to zero over the sphere; MC_q and diffusion
                // are small. So a `rest` that is COMPARABLE to the named terms means either the
                // transport is not conservative or one of the three named terms is wrong, and a
                // `rest` that is a rounding fraction of them means the split has accounted for
                // the bucket. Printed as that ratio so the table says which, rather than
                // leaving a large residual looking like a result.
                double big = std::abs(lf);
                if(std::abs(micro) > big) big = std::abs(micro);
                if(std::abs(f_fin) > big) big = std::abs(f_fin);
                cout << "      AGCM: [CWB]     |rest| / largest named term = " << fixed
                     << setprecision(4) << ((big > 0.0) ? std::abs(rest) / big : 0.0)
                     << "   <- small means the split accounts for the bucket" << endl;
            }

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
