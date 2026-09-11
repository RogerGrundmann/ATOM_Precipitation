#pragma once

#include "cAtmosphereModel.h"
#include "Utils.h"

#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <algorithm>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace AtomUtils;

class VelocityInitializer {
public:
    explicit VelocityInitializer(cAtmosphereModel& model)
        : m(model)
    {}

    void compute()
    {
        using namespace std;
        cout << endl << "      AGCM: init_velocities" << endl;

        // u-component up to tropopause and back on half distance
        init_u(m.u,  0);
        init_u(m.u, 30);
        init_u(m.u, 60);
        init_u(m.u, 90);
        init_u(m.u,120);
        init_u(m.u,150);
        init_u(m.u,180);

        // initialise v: tropopause and surface values per latitude
        // equator
        init_v_or_w(m.v,  90,  0.0,  0.0);                              // lat:   0   j=90
        // ---- THE POLAR CELLS ARE THE ONLY ROWS WHOSE TWO ENDPOINTS HAVE THE SAME SIGN ----
        //
        // `init_v_or_w(v, j, coeff_trop, coeff_sl)` ramps linearly from the SURFACE value to the
        // TROPOPAUSE value, so what makes a row an overturning CELL is the DIFFERENCE between the
        // two -- the shear -- and every other row straddles zero:
        //
        //     Hadley  j=75   (-3.0, +4.0)   shear 7.0, sign reversal
        //     Ferrel  j=45   (+4.0, -1.5)   shear 5.5, sign reversal
        //     polar   j=15   (+0.5, +0.6)   shear 0.1, NO sign reversal
        //
        // 0.1 m/s against 7.0 is not a weak polar cell, it is an almost uniform 0.55 m/s drift
        // through the whole troposphere. What appears in `Psi` as a polar cell is `ATM_V_MASSBAL`
        // (default on since 2026-08-28) subtracting the density-weighted column mean afterwards,
        // which turns ANY same-signed column into a sign-reversing one of amplitude ~= the shear.
        // So the polar cell's strength is set by an accidental 0.1, and it measures 1.9 % of the
        // Hadley cell's `Psi` where the real atmosphere has 10-20 % (2026-09-09; the sense is
        // right -- equatorward at the surface, poleward aloft -- only the amplitude is not).
        //
        // ATM_POLAR_CELL_SHEAR=<m/s> sets that shear about the SAME column mean of 0.55, so the
        // default 0.1 reproduces (0.5, 0.6) EXACTLY and is bit-identical. ~0.55 puts `Psi` at
        // roughly a tenth of the Hadley cell, which is where observations put it.
        // The unset branch returns the shipped LITERALS rather than recomputing them, because
        // 0.55 + 0.5*0.1 is 0.6000000000000001 in double and not 0.6 -- so the arithmetic form
        // would NOT be bit-identical off, which is this tree's standing requirement for a knob.
        static const bool   polar_set   = (getenv("ATM_POLAR_CELL_SHEAR") != nullptr);
        static const double polar_shear = [](){
            const char* e = getenv("ATM_POLAR_CELL_SHEAR"); return e ? atof(e) : 0.1; }();
        const double pc_mean = 0.55;
        const double pc_trop = polar_set ? (pc_mean - 0.5 * polar_shear) : 0.5;
        const double pc_sl   = polar_set ? (pc_mean + 0.5 * polar_shear) : 0.6;
        // northern polar cell
        init_v_or_w(m.v,   0,  0.5,  0.0);                              // lat:  90   j=0
        init_v_or_w(m.v,  15,  pc_trop, pc_sl);                         // lat:  75   j=15
        // southern polar cell
        init_v_or_w(m.v, 180,  0.5,  0.0);                              // lat: -90   j=180
        init_v_or_w(m.v, 165,  pc_trop, pc_sl);                         // lat: -75   j=165
        // northern Ferrel cell
        init_v_or_w(m.v,  30, -0.2,  0.0);                              // lat:  60   j=30
        init_v_or_w(m.v,  45,  4.0, -1.5);                              // lat:  45   j=45
        // southern Ferrel cell
        init_v_or_w(m.v, 150, -0.2,  0.0);                              // lat: -60   j=150
        init_v_or_w(m.v, 135,  4.0, -1.5);                              // lat: -45   j=135
        // ---- THE TWO HADLEY CELLS ARE ASYMMETRIC, AND IT IS A HALF-APPLIED EDIT ----
        //
        // `git log --follow` on this file, every version: the pair was SYMMETRIC at 3.0/3.0 in
        // `f03ff0b` (initial commit) and `ee821ab`, and the asymmetry arrives in `24ff23a`
        // ("revive Hadley/Ferrel cells", 2026-06-19), which raised the NORTHERN surface branch
        // 3.0 -> 4.0 and left the southern one at 3.0 -- with the intended southern replacement
        // WRITTEN AND COMMENTED OUT directly above it. Both were meant to go to 4.0; the north
        // got the edit and the south got the comment. Fifth instance of this tree's "a comment
        // describing code that does not run" class.
        //
        // Measured consequence: `Psi` at iteration 20 is -142.0 north against +123.4 south
        // (1.15x), and precipitation at iteration 20 is NH/SH **1.42** where NASA is 1.04 --
        // 2.9x in the subtropics and 18x over the polar caps. Most of that is a spin-up
        // transient (by iteration 600 the ratio is 1.07), but the initial condition is lopsided
        // by construction and the transient is what a plot at low iteration counts shows.
        //
        // ATM_HADLEY_SL=<m/s> gives BOTH hemispheres the same surface branch. Unset returns the
        // shipped LITERALS (4.0 north, 3.0 south) and is bit-identical; 4.0 completes what
        // `24ff23a` intended, 3.0 restores the symmetric pair this model ran with for its first
        // two months. Do not read 3.5 as "the compromise" -- it splits the difference between a
        // value and a mistake.
        static const bool   hadley_set = (getenv("ATM_HADLEY_SL") != nullptr);
        static const double hadley_sl  = [](){
            const char* e = getenv("ATM_HADLEY_SL"); return e ? atof(e) : 4.0; }();
        // northern Hadley cell
        init_v_or_w(m.v,  60,  0.0,  0.5);                              // lat:  30   j=60
        init_v_or_w(m.v,  75, -3.0, hadley_set ? hadley_sl : 4.0);      // lat:  15   j=75
        // southern Hadley cell
        init_v_or_w(m.v, 120,  0.0,  0.5);                              // lat: -30   j=120
        init_v_or_w(m.v, 105, -3.0, hadley_set ? hadley_sl : 3.0);      // lat: -15   j=105

        // initialise w: tropopause and surface values per latitude.
        // w is the ZONAL jet (East+). The SURFACE value (2nd coeff) is what the
        // atm->ocean transfer hands to the ocean, so it must reproduce the observed
        // surface wind BANDS: trade EASTERLIES (w<0) through the tropics/subtropics
        // (0-30deg, peak ~15deg), mid-latitude WESTERLIES (w>0) peaking ~45deg, then
        // weakening poleward. The wind-stress CURL between the trade easterlies and
        // the mid-latitude westerlies is what drives the subtropical (anticyclonic)
        // gyre; the curl between the westerly max and the pole drives the subpolar
        // (cyclonic) gyre. The tropopause value (1st coeff) keeps the upper-level
        // westerly jets (subtropical jet strongest at 30deg). See wind-IC diagnosis
        // in project_hydro_ekman_sh_gyre.
        // equator
        init_v_or_w(m.w,  90, -3.0, -5.0);                             // lat:   0   j=90   easterly (equatorial)
        // northern polar cell
        init_v_or_w(m.w,   0,  0.0,  0.0);                              // lat:  90   j=0
        // southern polar cell
        init_v_or_w(m.w, 180,  0.0,  0.0);                              // lat: -90   j=180
        // northern Ferrel cell (mid-latitude westerlies, weakening to the pole)
        init_v_or_w(m.w,  30, 10.0,  6.0);                             // lat:  60   j=30   westerly
        // southern Ferrel cell
        init_v_or_w(m.w, 150, 10.0,  6.0);                             // lat: -60   j=150  westerly
        // northern subtropics — horse latitudes (trade/westerly transition, calm)
        init_v_or_w(m.w,  60, 30.0, -1.0);                             // lat:  30   j=60   weak easterly
        // southern subtropics
        init_v_or_w(m.w, 120, 30.0, -1.0);                             // lat: -30   j=120  weak easterly
        // northern westerly max at j=45
        init_v_or_w(m.w,  45, 15.0, 10.0);                             // lat:  45   j=45   westerly max
        // southern westerly max at j=135
        init_v_or_w(m.w, 135, 15.0, 10.0);                             // lat: -45   j=135  westerly max
        // northern trade-easterly max at j=75 (15N)
        init_v_or_w(m.w,  75,  5.0, -7.0);                             // lat:  15   j=75   easterly (trade max)
        // southern trade-easterly max at j=105 (15S)
        init_v_or_w(m.w, 105,  5.0, -7.0);                             // lat: -15   j=105  easterly (trade max)

        // forming diagonals — northern hemisphere
        form_diagonals(m.u,  0,  30);
        form_diagonals(m.w,  0,  30);
        form_diagonals(m.w, 30,  45);
        form_diagonals(m.v,  0,  15);
        form_diagonals(m.v, 15,  30);

        form_diagonals(m.u, 30,  60);
        form_diagonals(m.w, 45,  60);
        form_diagonals(m.v, 30,  45);
        form_diagonals(m.v, 45,  60);

        form_diagonals(m.u, 60,  90);
        form_diagonals(m.w, 60,  75);                                   // 30N->15N (trade node at j=75)
        form_diagonals(m.w, 75,  90);                                   // 15N->0
        form_diagonals(m.v, 60,  75);
        form_diagonals(m.v, 75,  90);

        // forming diagonals — southern hemisphere
        form_diagonals(m.u,  90, 120);
        form_diagonals(m.w,  90, 105);                                  // 0->15S (trade node at j=105)
        form_diagonals(m.w, 105, 120);                                 // 15S->30S
        form_diagonals(m.w, 120, 135);
        form_diagonals(m.v,  90, 105);
        form_diagonals(m.v, 105, 120);

        form_diagonals(m.u, 120, 150);
        form_diagonals(m.w, 135, 150);
        form_diagonals(m.v, 120, 135);
        form_diagonals(m.v, 135, 150);

        form_diagonals(m.u, 150, 180);
        form_diagonals(m.w, 150, 180);
        form_diagonals(m.v, 150, 165);
        form_diagonals(m.v, 165, 180);

        // Zero land cells; non-dimensionalise air cells — single fused pass
        const double inv_u_0 = 1.0 / m.u_0;

        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int k = 0; k < m.km; k++) {
                for (int j = 0; j < m.jm; j++) {
                    if (is_land(m.h, i, j, k)) {
                        m.u.x[i][j][k] = 0.0;
                        m.v.x[i][j][k] = 0.0;
                        m.w.x[i][j][k] = 0.0;
                    } else {
                        m.u.x[i][j][k] *= inv_u_0;
                        m.w.x[i][j][k] *= inv_u_0;
                        if (!m.use_NASA_velocity && j > 90) {
                            m.v.x[i][j][k] = -m.v.x[i][j][k] * inv_u_0;
                        } else {
                            m.v.x[i][j][k] *= inv_u_0;
                        }
                    }
                }
            }
        }
/*
        // Surface taper on v and w (the two HORIZONTAL components in this model's
        // (r,θ,φ) convention: v = meridional, w = zonal): linearly damp from the
        // local value at i=5 down to zero at i=0, so the lowest five layers carry no
        // horizontal wind at the ground reference and grow smoothly into the
        // prescribed profile above. u is left alone because u is the RADIAL/VERTICAL
        // velocity here (NOT the zonal jet — that is w); init_u already gives it a
        // small profile that ramps to zero at the surface.
        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                for (int i = 0; i <= 5; i++) {
                    const double factor = static_cast<double>(i) / 5.0;
                    m.v.x[i][j][k] *= factor;
                    m.w.x[i][j][k] *= factor;
                }
            }
        }
*/
    cout << "      AGCM: init_velocities ended" << endl;
    }

    // Linear blend of u/v/w across j in [lat-3, lat+3].
    // Currently not called by compute() (dead code in the original),
    // but kept here as it logically belongs with velocity initialisation.
    void smooth_transition(int lat)
    {
        const int    start     = lat - 3;
        const int    end       = lat + 3;
        const double inv_range = 1.0 / (double)(end - start);

        #pragma omp parallel for collapse(2) schedule(static)
        for (int k = 0; k < m.km; k++) {
            for (int i = 0; i < m.im; i++) {
                const double u_start = m.u.x[i][start][k];
                const double v_start = m.v.x[i][start][k];
                const double w_start = m.w.x[i][start][k];
                const double u_slope = (m.u.x[i][end][k] - u_start) * inv_range;
                const double v_slope = (m.v.x[i][end][k] - v_start) * inv_range;
                const double w_slope = (m.w.x[i][end][k] - w_start) * inv_range;
                for (int j = start; j <= end; j++) {
                    const double t     = (double)(j - start);
                    m.u.x[i][j][k] = u_slope * t + u_start;
                    m.v.x[i][j][k] = v_slope * t + v_start;
                    m.w.x[i][j][k] = w_slope * t + w_start;
                }
            }
        }
    }

    // ==================================================================
    // COLUMN MASS-FLUX BALANCE (ATM_V_MASSBAL, DEFAULT ON since 2026-08-28; =0 restores)
    //
    // THE PRESCRIBED CELL IS NOT MASS-BALANCED, AND NOTHING DOWNSTREAM FIXES IT.
    // init_v_or_w() builds v as a LINEAR RAMP IN HEIGHT from coeff_sl at the surface to
    // coeff_trop at the tropopause, then a linear decay to zero at the lid. Nothing in that
    // construction constrains INT(rho*v*dz) = 0, or even INT(v*dz) = 0: for a linear ramp the
    // volume integral is H*(v_s + v_t)/2, which vanishes only if the two hand-set endpoints are
    // exact opposites, and the MASS integral needs something different again because rho decays
    // roughly exponentially while the ramp is linear in z.
    //
    // Measured consequence (ATOM_Precipitation, 2026-08-27): Psi(lid) = 0.0000e+00 exactly and
    // Psi(ground) rms = 1.55e+11 kg/s -- about 40 % of the Hadley cell's own strength -- with
    // u pinned to 0 at BOTH walls (BC_Atm.h:679, :697, verified in the field). With no flux
    // through either wall, div(rho v) = 0 would force Psi(ground) = 0. At 15N the column runs
    // +3.67 m/s at the ground and -0.54 m/s at 7.4 km: the poleward branch sits in the dense
    // lower 6 km and the return in thin air above, so the two cannot cancel. INT(v*dz) is
    // +9151 m^2/s there, so it fails to close in VOLUME as well -- two defects stacked, and the
    // density weighting would still be wrong if the endpoints were opposites.
    //
    // This is the family's rho-blindness, one file upstream of where it was already caught:
    // MinMax_Atm.cpp's `const double rho = r_air` was the DIAGNOSTIC version (dabbc94 in ATHAD,
    // ede4810 here, 2.28x overweight at the cell core). The instrument was corrected; the thing
    // it measures never was.
    //
    // Neither existing lever touches it, both measured on this quantity rather than inferred:
    // ATM_PROJ_SWEEPS -0.085 % at 100x (ATHAD gets -52.5 % at 10), and ATM_RHIE_CHOW +0.005 %,
    // which is structural -- D4 annihilates smooth fields by construction, and Psi(ground) is a
    // domain-scale quantity, so that knob CANNOT act on it. The checkerboard and this
    // non-closure are not the same defect.
    //
    // The repair is the initial-condition analogue of ATHAD's initBalancedState: impose the
    // constraint rather than hope the projection removes it. Per fluid column, subtract the
    // density-weighted column mean,
    //
    //     v <- v - INT(rho*v*dz) / INT(rho*dz)
    //
    // which makes INT(rho*v*dz) = 0 exactly while changing the profile by a single constant, so
    // the SHEAR that defines the cell -- and hence the overturning -- is untouched. It is the
    // minimum-norm correction that satisfies the constraint.
    //
    // Runs after densities(), because it needs r_humid; r_air is the documented fallback, the
    // same one write_meridional_streamfunction uses for pre-densities() cells. v only: the
    // constraint is on the meridional overturning, and w is zonal.
    // ==================================================================
    // ==================================================================
    // SIX CLOSED CELLS, BUILT FROM A STREAMFUNCTION
    // ATM_CELLS_FROM_PSI=1, DEFAULT 0 = OFF and byte-identical unset.
    //
    // WHY. The prescribed v is a LINEAR RAMP in height between two hand-set endpoints at nine
    // latitudes (init_v_or_w above). Nothing in that construction constrains INT(rho*v*dz), so
    // Psi(ground) starts at 1.55e+11 kg/s where it must be ZERO, and balance_column_mass_flux()
    // then removes 94.8 % of it AFTER the fact. A ramp also ignores rho, so the return branch
    // aloft is not the mass-compensating counterpart of the surface branch -- the cell cannot
    // close even in principle, only be corrected. And the polar rows (0.5, 0.6) have a shear of
    // 0.1 m/s against Hadley's 7.0 and no sign reversal at all: unchanged in EVERY commit since
    // `f03ff0b`, so this model has never had a closed polar cell.
    //
    // THE FIX IS TO BUILD v FROM Psi INSTEAD OF BUILDING Psi FROM v. Prescribe
    //
    //     Psi(phi,z) = A(phi) * sin( pi * (z - z_g) / (H - z_g) )
    //
    // per column, where z_g is that column's own ground (get_layer_height(i_topography)) and H
    // the model lid. Psi vanishes at BOTH ends by construction, so INT(rho*v*dz) = 0 EXACTLY in
    // every column and the cells close at every latitude at iteration 0. Differentiating the
    // diagnostic's own definition, Psi(z) = 2*pi*a*cos(phi) * INT_z^H rho*v dz (MinMax_Atm.cpp),
    //
    //     v(phi,z) = -(1 / (2*pi*a*cos(phi)*rho(z))) * dPsi/dz
    //
    // which carries 1/rho, so the return branch is automatically STRONGER in thin air aloft --
    // the property a linear ramp cannot express and the reason the shipped cell needs a
    // correction to close.
    //
    // A(phi) gives three cells per hemisphere with edges at 0/30/60/90 degrees, Psi = 0 at every
    // edge so the cells are genuinely separate, and alternating sign so Hadley is thermally
    // direct, Ferrel indirect and polar direct:
    //
    //     A(phi) = -sign(phi) * (-1)^n * A_n * sin( pi*(|phi| - phi_lo) / 30 )
    //
    // with n = 0,1,2 for Hadley/Ferrel/polar. The leading -sign(phi) is what makes the northern
    // surface branch equatorward; note the shipped path gets its hemispheric antisymmetry from a
    // separate `if (j > 90) v = -v` in the non-dimensionalisation pass, which this function does
    // NOT go through -- it writes non-dimensional v directly, after that pass.
    //
    // AMPLITUDES ARE OBSERVED VALUES, not the model's current ones: 120 / 40 / 12 * 1e9 kg/s
    // against a shipped Hadley of ~200e9 and a Ferrel of ~150e9 at iteration 20. So this makes
    // the two strong cells WEAKER and the polar cell exist for the first time. Each is a knob.
    //
    // THE BUILT-IN CHECK: balance_column_mass_flux() runs immediately after this and prints the
    // correction it applies. If the cells really close, that correction collapses to round-off.
    // A non-zero value there means this function is wrong, and it is measured rather than argued.
    void install_cells_from_streamfunction()
    {
        static const bool on = [](){ const char* e = getenv("ATM_CELLS_FROM_PSI");
                                     return e && atoi(e) != 0; }();
        if (!on) return;
        auto amp = [](const char* n, double d){ const char* e = getenv(n);
                                                return (e ? atof(e) : d) * 1.0e9; };
        const double A_had = amp("ATM_PSI_HADLEY", 120.0);
        const double A_fer = amp("ATM_PSI_FERREL",  40.0);
        // 26 rather than an observed 12: the taper and the ground-to-tropopause span cut the
        // installed cell to ~46 % of the prescribed value in the north and ~22 % in the south,
        // so the knob is not the cell. Measured at 26: 75N reaches 12.10e9 = 10.1 % of Hadley
        // and 75S 5.85e9 = 4.9 %, both inside the observed 5-15 %. The asymmetry is Antarctica
        // shortening the span and must NOT be tuned out -- forcing the south to 10 % needs ~53,
        // which puts the north at 20 %.
        const double A_pol = amp("ATM_PSI_POLAR",   26.0);
        const double a_E   = m.r_Earth * 1000.0;                 // r_Earth is in km
        const double inv_u0 = 1.0 / m.u_0;
        const double cos60 = cos(60.0 * M_PI / 180.0);
        double vmax = 0.0, worst_res = 0.0; long ncol = 0;

        #pragma omp parallel for collapse(2) schedule(static) \
                reduction(max:vmax) reduction(max:worst_res) reduction(+:ncol)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                const double lat  = 90.0 - (double)j * 180.0 / (double)(m.jm - 1);
                const double cphi = cos(lat * M_PI / 180.0);
                if (fabs(cphi) < 1.0e-6) {                       // the pole rows carry no cell
                    for (int i = 0; i < m.im; i++) m.v.x[i][j][k] = 0.0;   // v -> 0 at the pole
                    continue;
                }
                const double alat = fabs(lat);
                const int    n    = (alat < 30.0) ? 0 : (alat < 60.0) ? 1 : 2;
                const double lo   = (n == 0) ? 0.0 : (n == 1) ? 30.0 : 60.0;
                const double A_n  = (n == 0) ? A_had : (n == 1) ? A_fer : A_pol;
                const double sgn  = (lat >= 0.0 ? -1.0 : 1.0) * ((n % 2 == 0) ? 1.0 : -1.0);
                // POLAR TAPER. v carries 1/cos(phi), so a cell of fixed Psi amplitude implies a
                // wind that runs away toward the pole. Psi must vanish at the pole in any case,
                // so taper the amplitude with cos(phi), normalised to 1 at 60 deg: equatorward
                // of 60 this is inert, and it bounds v in the polar band without a cap (a cap
                // would break the closure this whole construction exists to provide).
                const double taper = (alat <= 60.0) ? 1.0 : (cphi / cos60);
                const double A     = sgn * A_n * taper * sin(M_PI * (alat - lo) / 30.0);

                const int i0 = m.i_topography[j][k];
                // CONFINE THE CELL TO THE TROPOSPHERE. The overturning is a tropospheric
                // circulation, and carrying it to the 16 km lid also divides by a stratospheric
                // rho ~0.1, which is where the unphysical 28.5 m/s in the first version came
                // from. Psi = 0 from the tropopause up, so the cell closes below it.
                int it = m.get_tropopause_layer(j);
                if (it >= m.im) it = m.im - 1;
                // WHERE THE GROUND REACHES THE TROPOPAUSE THERE IS NO ROOM FOR A CELL, so set
                // the column to zero and move on. It must not be SKIPPED -- a skipped column
                // keeps the old ramp, whose integral does not vanish, and balance_column_mass_flux
                // then has something to remove after all (the first version installed 61 451 of
                // 65 341 columns and its 4.3e-02 correction was entirely those 3 890). And it
                // must not fall back to the LID either, which is what the second version did:
                // over the Antarctic plateau that gave those columns a full-depth cell, and in
                // the zonal mean it appeared as a second circulation ABOVE the tropopause
                // reaching 2.04e9 kg/s at 10.9 km at 75S against only 0.24e9 at 75N -- the
                // hemispheric asymmetry being the tell, since the Arctic has no such terrain.
                // Zero has a vanishing column integral, so the mass-balance check still holds.
                if (it <= i0) {
                    for (int i = 0; i < m.im; i++) m.v.x[i][j][k] = 0.0;
                    continue;
                }
                const double z_g  = m.get_layer_height(i0);
                const double span = m.get_layer_height(it) - z_g;
                if (!(span > 0.0)) continue;
                ncol++;

                // Pass 1: the analytic profile, rv = rho*v at each level.
                //
                // ATM_PSI_SHAPE selects the VERTICAL shape of the cell, and it matters more than
                // it looks:
                //
                //   0  Psi = A*sin(pi*zeta),  v ~ cos(pi*zeta)/rho
                //      The obvious choice and it is WRONG AT THE TOP. cos is maximal at BOTH
                //      ends and 1/rho is also largest at the top, so the return branch is
                //      amplified twice at the same boundary: measured at 15N it grows
                //      monotonically to -2.28 m/s at 12 030 m -- the fastest flow in the column
                //      -- and then drops DISCONTINUOUSLY to zero above the tropopause. In
                //      streamlines that is a thin fast sheet pinned to the lid, which is not
                //      what a Hadley cell looks like. (Found by the user looking at the plot;
                //      no closure or mass check could see it, they all passed.)
                //
                //   1  Psi = A*(27/4)*zeta*(1-zeta)^2,  v ~ (1-zeta)*(1-3*zeta)/rho   [DEFAULT]
                //      v is maximal at the surface, reverses at zeta = 1/3 so the return flow
                //      occupies the upper TWO THIRDS as a real cell does, and goes to zero
                //      SMOOTHLY at the tropopause -- no discontinuity, no sheet. Psi peaks at
                //      zeta = 1/3 with value (4/27)*A, so the 27/4 normalises the knob to mean
                //      the cell's actual maximum. Closure is untouched: Psi = 0 at both ends.
                static const int shape = [](){ const char* e = getenv("ATM_PSI_SHAPE");
                                               return e ? atoi(e) : 1; }();
                std::vector<double> rv(m.im, 0.0), wq(m.im, 0.0);
                for (int i = i0; i <= it; i++) {
                    double rho = m.r_humid.x[i][j][k];
                    if (!AtomUtils::is_finite_safe(rho) || rho <= 0.0) rho = m.r_air;
                    const double zeta = (m.get_layer_height(i) - z_g) / span;
                    const double dPsi = (shape == 0)
                        ? A * (M_PI / span) * cos(M_PI * zeta)
                        : A * (27.0 / 4.0) / span * (1.0 - zeta) * (1.0 - 3.0 * zeta);
                    rv[i] = -dPsi / (2.0 * M_PI * a_E * cphi);                     // = rho*v
                }
                // Trapezoid weights, exactly those balance_column_mass_flux and the Psi
                // diagnostic use, so "the discrete integral vanishes" means the same thing here
                // as it does there.
                for (int i = i0; i < it; i++) {
                    const double dz = m.get_layer_height(i+1) - m.get_layer_height(i);
                    if (!(dz > 0.0)) continue;
                    wq[i] += 0.5 * dz;  wq[i+1] += 0.5 * dz;
                }
                // THE TOP HALF-INTERVAL, which the first version omitted and which was the whole
                // of the 3.6e-02 the balance still had to remove. Psi = A*sin(pi*zeta) puts v at
                // its EXTREMUM at the tropopause (cos(pi) = -1) and this function sets v = 0
                // above it, so balance_column_mass_flux -- which integrates i0..im-2 -- picks up
                // a trapezoid term 0.5*rv[it]*dz over [it, it+1] that was not in these weights.
                // Include it, and the integral this function zeroes becomes exactly the integral
                // that routine computes.
                if (it + 1 < m.im) {
                    const double dz_top = m.get_layer_height(it+1) - m.get_layer_height(it);
                    if (dz_top > 0.0) wq[it] += 0.5 * dz_top;
                }
                // Pass 2: MAKE THE DISCRETE INTEGRAL VANISH BY BALANCING THE TWO BRANCHES.
                // The analytic derivative integrates to zero in the CONTINUUM; on a 41-level
                // grid stretched 23x it does not, and the first version left 4.1e-02 (non-dim v,
                // ~10 % of the cell) for balance_column_mass_flux to remove as a CONSTANT SHIFT
                // -- which moves the zero-crossing and distorts the cell. Scaling the poleward
                // and equatorward branches instead makes them carry equal and opposite mass
                // flux, which is what a closed cell IS, and leaves the sign structure and the
                // crossing height untouched.
                double P = 0.0, N = 0.0;
                for (int i = i0; i <= it; i++) {
                    const double c = rv[i] * wq[i];
                    if (c > 0.0) P += c; else N -= c;
                }
                // A column needs BOTH branches present for this to work. With shape 0 that was
                // automatic -- cos(pi*zeta) is +1 at zeta=0 and -1 at zeta=1, so even a 2-level
                // column has one of each. With shape 1, v VANISHES at zeta=1, so a column with
                // very few levels above its ground can sample only one sign, P or N comes out
                // zero, no balancing is possible and the column's integral does not vanish:
                // measured, that took the worst relative residual from 3.6e-16 straight to
                // 1.00e+00 and put 7.0e-02 back into balance_column_mass_flux. Zero those
                // columns -- a zero column has a vanishing integral, so closure survives, and
                // they are shallow columns that cannot hold a cell in any case.
                if (!(P > 0.0 && N > 0.0)) {
                    for (int i = 0; i < m.im; i++) m.v.x[i][j][k] = 0.0;
                    continue;
                }
                {
                    const double half = 0.5 * (P + N);
                    const double fp = half / P, fn = half / N;
                    for (int i = i0; i <= it; i++) rv[i] *= (rv[i] > 0.0 ? fp : fn);
                }
                // residual of the discrete integral, reported rather than assumed
                double res = 0.0, scale = 0.0;
                for (int i = i0; i <= it; i++) { res += rv[i]*wq[i]; scale += fabs(rv[i])*wq[i]; }
                if (scale > 0.0 && fabs(res)/scale > worst_res) worst_res = fabs(res)/scale;

                for (int i = 0; i < m.im; i++) {
                    if (i < i0 || i > it) { m.v.x[i][j][k] = 0.0; continue; }
                    double rho = m.r_humid.x[i][j][k];
                    if (!AtomUtils::is_finite_safe(rho) || rho <= 0.0) rho = m.r_air;
                    const double v_phys = rv[i] / rho;
                    m.v.x[i][j][k] = v_phys * inv_u0;
                    if (fabs(v_phys) > vmax) vmax = fabs(v_phys);
                }
            }
        }
        std::cout << "      AGCM: [CELLS FROM PSI] installed over " << ncol
                  << " columns;  Psi Hadley/Ferrel/polar = "
                  << A_had/1e9 << " / " << A_fer/1e9 << " / " << A_pol/1e9 << " e9 kg/s;  max |v| = "
                  << vmax << " m/s;  worst relative residual of the discrete column integral = "
                  << std::scientific << std::setprecision(2) << worst_res << std::defaultfloat
                  << "  (balance_column_mass_flux should now find ~nothing to remove)"
                  << std::endl;
    }


    // ================================================================================
    // THE CELLS HAD NO VERTICAL VELOCITY OF THEIR OWN: ATM_CELLS_U_FROM_PSI
    // Default 0 = off and byte-identical. Requires ATM_CELLS_FROM_PSI=1; inert without it.
    //
    // THE DEFECT, found by the user looking at the glyphs (2026-09-11): the arrows point UPWARD
    // over most of the polar band, and at the pole itself they should descend. They do descend
    // -- but only poleward of 77.7 deg, and the cells are scored at 75.
    //
    // install_cells_from_streamfunction() writes v AND NOTHING ELSE. u keeps whatever init_u
    // left, and init_u prescribes only four latitudes per hemisphere -- 0, 30, 60, 90 -- with
    //     ua_00, ua_30, ua_60, ua_90 = 0.02894, 0.02315, 0.01736, 0.011574
    // which is 0.02894 * (1, 0.8, 0.6, 0.4): a smooth ramp with an alternating sign, not a
    // balance of anything. form_diagonals(m.u, 0, 30) then blends LINEARLY from -ua_90 at the
    // pole to +ua_60 at 60 deg, so the vertical velocity changes sign at
    //     j = 30 * 0.011574 / (0.011574 + 0.01736) = 12.0  ->  78.0 deg
    // against a MEASURED 77.7 deg (87E, i=14, post-projection initial state). The descending
    // branch of the polar cell therefore gets 12 degrees of latitude and the ascending branch
    // 18, and cos-weighted the descent is only ~16 % of the cell's area.
    //
    // So in the ATM_CELLS_FROM_PSI configuration the glyphs draw a streamfunction-derived v
    // against an unrelated analytic u, and project_initial_velocity(200) does not repair it --
    // the 77.7 deg crossing is measured AFTER the projection.
    //
    // THE REPAIR IS CONTINUITY, NOT THE ANALYTIC Psi. Deriving u from A*dPsi/dphi would be
    // wrong here: Pass 2 above RESCALES the two branches by fp/fn to make the discrete column
    // integral vanish, so after that rv is no longer exactly -dPsi/dz, and each column carries
    // its own i0, span and taper. Integrating the model's OWN discrete continuity against the v
    // that was actually written makes u consistent with the field the model holds rather than
    // with the field the formula describes.
    //
    //   v here is SOUTHWARD-positive, so with V = -v northward and f = rho*v*cos(phi):
    //     (1/(a cos phi)) d(rho V cos phi)/dphi + d(rho u)/dz = 0
    //   and phi(j) decreases with j (j=0 is 90N, dphi/dj = -PI/(jm-1) = -Delta), so
    //     d(rho u)/dz = -(1/(a cos phi)) * (1/Delta) * df/dj
    //   integrated upward from rho*u = 0 at the ground.
    //
    // THE BUILT-IN CHECK. Integrating that over the whole column gives
    // d/dj of INT(rho v cos phi dz), and install_cells has just made INT(rho v dz) vanish in
    // every column -- so rho*u must return to ZERO at the top of a column whose neighbours also
    // closed. The residual at the top is printed rather than assumed, as the top half-interval
    // omission was caught that way before.
    //
    // NOT TOUCHED: the pole rows (cos phi -> 0 divides here, and they carry no cell), land
    // cells, and columns the loop above zeroed. u above the tropopause is left at whatever the
    // integration returns, which is the residual and should be ~0.
    void install_u_from_cells()
    {
        static const bool on = [](){ const char* e = getenv("ATM_CELLS_U_FROM_PSI");
                                     return e && atoi(e) != 0; }();
        if (!on) return;
        static const bool cells_on = [](){ const char* e = getenv("ATM_CELLS_FROM_PSI");
                                           return e && atoi(e) != 0; }();
        if (!cells_on) {
            std::cout << "      AGCM: [CELLS U] ATM_CELLS_U_FROM_PSI set but ATM_CELLS_FROM_PSI"
                      << " is off -- nothing to derive u from, ignored." << std::endl;
            return;
        }
        const double a_E    = m.r_Earth * 1000.0;
        const double inv_u0 = 1.0 / m.u_0;
        const double Delta  = M_PI / (double)(m.jm - 1);          // radians of latitude per j
        const double dlam   = 2.0 * M_PI / (double)(m.km - 1);    // radians of longitude per k
        double umax = 0.0, worst_top = 0.0, scale_top = 0.0;
        std::vector<double> resid;   // per-column |rho*u(top)| / peak |rho*u|, for the quantiles
        std::vector<double> imbal;   // per-column |P-N|/(P+N) BEFORE balancing -- how much was needed
        long nzeroed = 0;

        auto rho_at = [&](int i, int j, int k){
            double r = m.r_humid.x[i][j][k];
            if (!AtomUtils::is_finite_safe(r) || r <= 0.0) r = m.r_air;
            return r;
        };
        // f = rho * v_phys * cos(phi); land and out-of-cell cells carry v = 0, so f = 0 there
        auto f_at = [&](int i, int j, int k){
            const double lat = 90.0 - (double)j * 180.0 / (double)(m.jm - 1);
            return rho_at(i,j,k) * m.v.x[i][j][k] * m.u_0 * cos(lat * M_PI / 180.0);
        };
        // g = rho * w_phys, the ZONAL mass flux. Omitting this was wrong and it is not a small
        // term: w is the jet, 10-30 m/s against v's ~1, and it carries strong zonal structure
        // over topography. Continuity in 3-D is
        //     d(rho u)/dz = -(1/(a cos phi)) * [ d(rho V cos phi)/dphi + d(rho w)/dlambda ]
        // and a meridional-only version is a ZONAL-MEAN statement being applied per column.
        auto g_at = [&](int i, int j, int k){
            return rho_at(i,j,k) * m.w.x[i][j][k] * m.u_0;
        };

        #pragma omp parallel for collapse(2) schedule(static) \
                reduction(max:umax) reduction(max:worst_top) reduction(max:scale_top)
        for (int j = 1; j < m.jm - 1; j++) {                       // pole rows excluded
            for (int k = 0; k < m.km; k++) {
                const double lat  = 90.0 - (double)j * 180.0 / (double)(m.jm - 1);
                const double cphi = cos(lat * M_PI / 180.0);
                if (fabs(cphi) < 1.0e-6) continue;
                const int i0 = m.i_topography[j][k];

                // d(rho u)/dz at each level: centred in j for the meridional part and in k
                // (periodic) for the zonal part.
                const int km1 = (k == 0) ? m.km - 2 : k - 1;          // phi is periodic
                const int kp1 = (k == m.km - 1) ? 1 : k + 1;
                std::vector<double> src(m.im, 0.0);
                for (int i = i0; i < m.im; i++) {
                    const double dfdj = 0.5 * (f_at(i, j+1, k) - f_at(i, j-1, k));
                    const double dgdk = 0.5 * (g_at(i, j, kp1) - g_at(i, j, km1));
                    src[i] = -(dfdj / Delta + dgdk / dlam) / (a_E * cphi);
                }
                // MAKE THE COLUMN CLOSE AT BOTH RADIAL WALLS, WITH A BOUNDED CORRECTION.
                // Each column's INT(rho v dz) vanishes over ITS OWN range, but i0 and it differ
                // between neighbouring j, so on the landward side of every terrain step the
                // j-derivative integrates over a range the neighbour does not share and the
                // column does not close: measured without any correction the median column
                // closes exactly (0.0000) and the p95 leaves 100 % of its peak at the lid. A
                // non-zero rho*u at a RIGID LID is a worse error than a redistribution, and
                // both radial walls are no-flow, so force both endpoints.
                //
                // NOT by scaling the two signs the way Pass 2 does it for v. That was tried and
                // it is catastrophic here: with P or N near zero the factor half/P is unbounded,
                // and a worst-column imbalance of 0.9998 produced max|u| = 84 m/s. Pass 2 gets
                // away with it because v's two branches are set by an analytic shape that always
                // has both; this source is a j-DERIVATIVE and can legitimately be nearly
                // one-signed next to terrain.
                //
                // Subtracting a correction proportional to |src| instead is bounded by
                // construction: the removed amount at each level is R*|src[i]|/INT|src|, which
                // is at most |src[i]|, so no level can more than double. It also cannot act
                // where src is zero.
                {
                    double R = 0.0, S = 0.0;
                    for (int i = i0; i < m.im - 1; i++) {
                        const double dz = m.get_layer_height(i+1) - m.get_layer_height(i);
                        if (!(dz > 0.0)) continue;
                        R += 0.5 * (src[i] + src[i+1]) * dz;
                        S += 0.5 * (fabs(src[i]) + fabs(src[i+1])) * dz;
                    }
                    if (S > 0.0) {
                        const double imb = fabs(R) / S;
                        #pragma omp critical
                        { imbal.push_back(imb); }
                        // A column that is essentially one-signed cannot hold a closed vertical
                        // circulation; zero it rather than flatten it to noise.
                        if (imb > 0.9) {
                            for (int i = 0; i < m.im; i++) m.u.x[i][j][k] = 0.0;
                            #pragma omp critical
                            { nzeroed++; }
                            continue;
                        }
                        const double corr = R / S;
                        for (int i = i0; i < m.im; i++) src[i] -= corr * fabs(src[i]);
                    }
                }
                // integrate upward, trapezoid, from rho*u = 0 at the ground
                double ru = 0.0, ru_peak = 0.0;
                for (int i = 0; i <= i0; i++) m.u.x[i][j][k] = 0.0;
                for (int i = i0; i < m.im - 1; i++) {
                    const double dz = m.get_layer_height(i+1) - m.get_layer_height(i);
                    if (dz > 0.0) ru += 0.5 * (src[i] + src[i+1]) * dz;
                    const double up = ru / rho_at(i+1, j, k);
                    m.u.x[i+1][j][k] = up * inv_u0;
                    if (fabs(up) > umax) umax = fabs(up);
                    if (fabs(ru) > ru_peak) ru_peak = fabs(ru);     // the cell's own mass flux
                }
                m.u.x[m.im-1][j][k] = 0.0;                          // rigid lid, as BC_Atm sets
                // THE CLOSURE RESIDUAL IS |rho*u| LEFT AT THE TOP, against that column's OWN
                // peak. Taking the peak instead of the endpoint -- which the first version of
                // this print did -- measures the cell rather than the residual and always reads
                // ~1. Seventh instrument-shaped defect in this tree, and mine.
                if (ru_peak > 0.0) {
                    const double rel = fabs(ru) / ru_peak;
                    if (rel > worst_top) worst_top = rel;
                    #pragma omp critical
                    { resid.push_back(rel); }
                }
                if (ru_peak > scale_top) scale_top = ru_peak;
            }
        }
        std::sort(resid.begin(), resid.end());
        auto q = [&](double f){ return resid.empty() ? 0.0
                     : resid[std::min(resid.size()-1, (std::size_t)(f*resid.size()))]; };
        std::cout << "      AGCM: [CELLS U] u rebuilt from discrete continuity against the"
                  << " installed v over " << resid.size() << " columns;  max |u| = "
                  << umax << " m/s;  peak column |rho*u| = " << std::scientific
                  << std::setprecision(2) << scale_top << std::defaultfloat
                  << ";  CLOSURE |rho*u(top)|/peak: median " << std::fixed << std::setprecision(4)
                  << q(0.5) << "  p95 " << q(0.95) << "  worst " << worst_top
                  << std::defaultfloat
                  << std::endl;
        std::sort(imbal.begin(), imbal.end());
        auto qi = [&](double f){ return imbal.empty() ? 0.0
                     : imbal[std::min(imbal.size()-1, (std::size_t)(f*imbal.size()))]; };
        std::cout << "      AGCM: [CELLS U] ascent/descent imbalance removed |P-N|/(P+N):"
                  << std::fixed << std::setprecision(4)
                  << "  median " << qi(0.5) << "  p95 " << qi(0.95) << "  worst " << qi(1.0)
                  << std::defaultfloat << ";  columns zeroed for carrying one sign only = "
                  << nzeroed << std::endl;
    }

    void balance_column_mass_flux()
    {
        if (!massBalance()) return;
        long n_cols = 0;
        const double worst = apply_column_mass_flux_balance(&n_cols);
        std::cout << "      ATOM: column mass-flux balance applied to " << n_cols
                  << " columns, largest correction " << std::scientific << std::setprecision(3)
                  << worst << " (non-dim v)" << std::endl;
    }

    // ==================================================================
    // THE SAME CONSTRAINT, RE-IMPOSED INSIDE THE TIME LOOP
    // ATM_V_MASSBAL_STRIDE=<N>, DEFAULT 0 = OFF and byte-identical unset.
    //
    // WHY. balance_column_mass_flux() above is called ONCE, before the loop, and it removes
    // 94.8 % of Psi(ground) at initialisation. It does not stay removed. Measured on
    // `output_cellpsi` -- the 600-iteration closed-cell run, whose Psi(ground) starts at
    // 3.5e-17 and does not -- the eroding tendency `dv_dyn` in the model's own v-momentum
    // budget is dominated by its COLUMN MEAN at every latitude that carries a cell:
    //
    //     lat        75N    45N    15N    15S    45S    75S
    //     mean/rms   2.57   3.70   1.30   1.36   8.23   3.18      (iteration 600)
    //
    // A column-mean tendency is not a cell, it is an OFFSET, and it is exactly the mode this
    // routine removes. The accounting closes: <dv_dyn> x 600 iterations is 0.055 m/s at 75N
    // and 0.126 at 45N, which through Psi = 2*pi*a*cos(phi)*INT(rho dz)*dv predicts 5.9e9 and
    // 3.7e10 kg/s against a MEASURED Psi(ground) of 4.97e9 and 3.23e10. Both within ~20 %.
    // So the closure loss is this drift accumulating, and the polar cells fail first because
    // they are the smallest -- not because anything polar is acting on them. (Two candidates
    // that ARE polar were checked and refuted: `dv_polar`, the polar zonal filter, is 1e-16 in
    // the zonal mean at every latitude including 75N/75S -- it conserves the zonal mean even
    // through its solid-neighbour substitution; and the radial filter is not preferentially
    // eating the polar cell, whose profile is SMOOTHER in index space than the Hadley one.)
    //
    // WHY THE CONSTRAINT IS LEGITIMATE EVERY ITERATION AND NOT ONLY AT SETUP. In a real
    // atmosphere INT(rho*v*dz) per column need not vanish instantaneously -- its divergence is
    // d(p_s)/dt. THIS MODEL HAS NO PROGNOSTIC SURFACE PRESSURE: `p_stat` is diagnosed
    // barometrically from the temperature by densities() and does not respond to the flow's
    // mass convergence at all (CLAUDE.md, the `dt` subsection). So a column-integrated mass
    // flux here has nothing to raise and nowhere to go, and letting one accumulate is
    // unphysical by construction rather than merely inconvenient.
    //
    // WHAT IT IS NOT. It is not a projection: it removes one number per column, the
    // density-weighted column mean of v, and leaves the SHEAR -- which is what defines the
    // cell -- untouched. It writes `v` only and lets storeIntermediateData3D sync `vn`, which
    // is the convention every filter in the loop already follows.
    //
    // CAVEAT ON THE INSTRUMENT: the call site sits AFTER write_v_momentum_budget, so the four
    // captured stages (dyn/polar/orog/radial) and their `dv_net` are the PRE-balance net. The
    // correction this routine applies is reported on its own line instead.
    static int massBalanceStride(){
        static const int v = [](){
            const char* e = getenv("ATM_V_MASSBAL_STRIDE"); return e ? atoi(e) : 0; }();
        return v;
    }

    void balance_column_mass_flux_in_loop(int iter_n)
    {
        const int stride = massBalanceStride();
        if (stride <= 0 || iter_n % stride != 0) return;
        const double worst = apply_column_mass_flux_balance(nullptr);
        // One line per VTK checkpoint rather than one per application: at stride 1 over 600
        // iterations the per-call print is 600 lines of log for one number.
        static int calls = 0;
        static double worst_max = 0.0;
        calls++;
        if (worst > worst_max) worst_max = worst;
        if (calls == 1 || iter_n % 100 == 0) {
            std::cout << "      ATOM: [V_MASSBAL_STRIDE] iter " << iter_n << ", call " << calls
                      << ", correction " << std::scientific << std::setprecision(3) << worst
                      << " (max so far " << worst_max << ", non-dim v)"
                      << std::defaultfloat << std::endl;
        }
    }

    // The shared core. Subtracts the density-weighted column mean of v from every fluid column
    // and returns the largest correction applied, in non-dimensional v.
    double apply_column_mass_flux_balance(long* n_cols_out)
    {
        long n_cols = 0; double worst = 0.0;
        #pragma omp parallel for collapse(2) schedule(static) reduction(+:n_cols) reduction(max:worst)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                const int i0 = m.i_topography[j][k];
                if (i0 >= m.im - 1) continue;
                double I1 = 0.0, I0 = 0.0;
                for (int i = i0; i < m.im - 1; i++) {
                    const double dz = m.get_layer_height(i+1) - m.get_layer_height(i);
                    if (!(dz > 0.0)) continue;
                    double r1 = m.r_humid.x[i][j][k], r2 = m.r_humid.x[i+1][j][k];
                    if (!AtomUtils::is_finite_safe(r1) || r1 <= 0.0) r1 = m.r_air;
                    if (!AtomUtils::is_finite_safe(r2) || r2 <= 0.0) r2 = m.r_air;
                    const double v1 = m.v.x[i][j][k], v2 = m.v.x[i+1][j][k];
                    if (!AtomUtils::is_finite_safe(v1) || !AtomUtils::is_finite_safe(v2)) continue;
                    I1 += 0.5 * (r1*v1 + r2*v2) * dz;
                    I0 += 0.5 * (r1    + r2   ) * dz;
                }
                if (!(I0 > 0.0)) continue;
                const double dv = I1 / I0;
                if (std::fabs(dv) > worst) worst = std::fabs(dv);
                for (int i = i0; i < m.im; i++) {
                    m.v.x[i][j][k]  -= dv;
                    m.vn.x[i][j][k]  = m.v.x[i][j][k];
                }
                n_cols++;
            }
        }
        // The walls are re-imposed by bcRadius afterwards; this only shifts interior v.
        if (n_cols_out) *n_cols_out = n_cols;
        return worst;
    }

    // DEFAULT ON since 2026-08-28. ATM_V_MASSBAL=0 restores the unbalanced prescribed profile
    // exactly, which is the branch every measurement recorded before that date was made on.
    static bool massBalance(){
        static const bool v = [](){
            const char* e = getenv("ATM_V_MASSBAL"); return e ? (atoi(e) != 0) : true; }();
        return v;
    }

    // ==================================================================
    // THERMAL-WIND BALANCED INITIAL STATE
    // ATM_TW_BALANCE=<strength>, DEFAULT 0.0 = OFF and byte-identical unset.
    //
    // WHY THIS IS AN INITIAL CONDITION AND NOT A FORCE. CLAUDE.md records that nothing in
    // rhs_u/rhs_v/rhs_w carries the temperature field: p_stat appears in no momentum equation,
    // p_dyn is a projection pressure, and the Boussinesq buoyancy -- the one surrogate route --
    // measures 0.03 % of the radial PGF on the shipped branch. Two repairs have been tried
    // against that from the FORCE side. ATM_HYDRO_PGF injects the hydrostatic gradient and takes
    // the jet-core ageostrophic residual 0.9999 -> 0.2498; ATM_BUOY_CONSISTENT restores the
    // buoyancy's non-dimensionalisation and takes the band p05 residual 0.996 -> 0.437 through
    // the model's OWN elliptic pressure. NEITHER MOVES THE VELOCITY, and the reason is not the
    // size of either term: geostrophic adjustment takes 1/f, which at 31 deg is 13 315 s =
    // 66 500 iterations at dt = 0.2 s. The longest run in this tree is 1600. The `dt` route to
    // reach it was tried and failed -- the usable ceiling is ~4x on the model's own
    // precipitation, where 1/f still needs 16 600 iterations.
    //
    // SO THE ADJUSTMENT MUST NOT BE WAITED FOR. IT MUST BE SUPPLIED. This is the same move
    // balance_column_mass_flux() makes one function above -- impose the constraint on the
    // initial state rather than hope the dynamics find it -- and that one removed 94.8 % of
    // Psi(ground) at initialisation where two solver knobs had moved it 0.085 % and 0.005 %.
    //
    // WHAT IT COMPUTES. Thermal wind for the ZONAL wind, which in this model's (r,theta,phi)
    // convention is w (East+), NOT u (u is radial):
    //
    //     dw/dz = -(g/(f*T)) * dT/dy,     y = northward distance
    //
    // the.z[j] is COLATITUDE (MinMax_Atm.cpp:235: "cos(latitude) = sin(colatitude)"), so
    // j increases southward, dy = -a*dtheta and f = 2*omega*cos(theta). Hence
    //
    //     dw/dz = +(g/(a*f*T)) * dT/dtheta
    //
    // Checked against the model's own field before it was written: at 31S on the 87E section
    // dT/dy = 0.683 K/deg over 1-10 km, which this relation turns into +28.1 m/s of westerly
    // shear where the model carries +4.56. The shear is present in the temperature and absent
    // from the wind, and that difference is exactly what this function installs.
    //
    // THE SURFACE VALUE IS THE ANCHOR AND IS NOT TOUCHED. Integration starts at
    // i_topography[j][k] from the prescribed w, so the surface wind BANDS survive intact --
    // they are what the atm->ocean transfer hands to the hydrosphere, and the trade/westerly
    // curl between them is what drives the gyres (project_hydro_ekman_sh_gyre). This function
    // adds the SHEAR its own temperature implies; it does not rewrite the wind.
    //
    // THE TROPICS ARE TAPERED OUT, because f -> 0 there and thermal wind is not the balance
    // that holds. weight = sin^2(lat)/(sin^2(lat) + sin^2(lat_min)) is smooth, 0 at the equator
    // and 0.5 at lat_min (ATM_TW_LATMIN, default 15 deg), so the prescribed Hadley/trade
    // profile is retained where it is the right structure and the geostrophic shear takes over
    // polewards. There is no 1/f blow-up anywhere: the weight vanishes faster than f does.
    //
    // A LEVEL WHOSE MERIDIONAL NEIGHBOURS ARE NOT BOTH FLUID CONTRIBUTES ZERO SHEAR. A centred
    // dT/dtheta straddling a terrain step is the defect that gave brunt_N2 its +-0.03 s^-2
    // "boundary-layer" extrema, which turned out to be the Andes and the Himalaya. Not repeated
    // here.
    //
    // ORDERING: this runs BEFORE balance_column_mass_flux(), so if the v-component is enabled
    // the mass constraint is applied to the final v rather than to a profile this then shifts.
    // ==================================================================
    void balance_thermal_wind()
    {
        const double s = twStrength();
        if (s == 0.0) return;

        const double a        = m.r_Earth * 1000.0;              // [m] (r_Earth is in km)
        const double lat_min  = twLatMin() * M_PI / 180.0;
        const double s2_min   = sin(lat_min) * sin(lat_min);
        const double w_max    = twWmax();                        // [m/s] sanity cap
        const double inv_u_0  = 1.0 / m.u_0;
        const bool   do_v     = twDoV();

        long n_cols = 0, n_capped = 0;
        double worst_dw = 0.0;                                   // largest shear added [m/s]

        #pragma omp parallel for collapse(2) schedule(static) \
                reduction(+:n_cols,n_capped) reduction(max:worst_dw)
        for (int j = 1; j < m.jm - 1; j++) {
            for (int k = 0; k < m.km; k++) {
                const double theta = m.the.z[j];
                const double slat  = cos(theta);                 // sin(latitude)
                const double f     = 2.0 * m.omega * slat;
                if (f == 0.0) continue;
                const double wgt = (slat*slat) / (slat*slat + s2_min);
                if (wgt <= 0.0) continue;

                const int i0 = m.i_topography[j][k];
                if (i0 >= m.im - 1) continue;

                const int kp = (k + 1) % m.km;
                const int km1 = (k + m.km - 1) % m.km;
                const double sth = sin(theta);

                // shear at a level, [ (m/s) / m ], zero where the stencil is not all fluid
                auto shear_w = [&](int i)->double {
                    if (is_land(m.h, i, j+1, k) || is_land(m.h, i, j-1, k)) return 0.0;
                    const double T = m.t.x[i][j][k] * m.t_0;
                    if (!(T > 0.0) || !AtomUtils::is_finite_safe(T)) return 0.0;
                    const double dTdthe = (m.t.x[i][j+1][k] - m.t.x[i][j-1][k])
                                        * m.t_0 / (2.0 * m.dthe);
                    if (!AtomUtils::is_finite_safe(dTdthe)) return 0.0;
                    return m.g * dTdthe / (a * f * T);
                };
                auto shear_v = [&](int i)->double {
                    if (!do_v) return 0.0;
                    if (is_land(m.h, i, j, kp) || is_land(m.h, i, j, km1)) return 0.0;
                    if (!(sth > 1.0e-3)) return 0.0;
                    const double T = m.t.x[i][j][k] * m.t_0;
                    if (!(T > 0.0) || !AtomUtils::is_finite_safe(T)) return 0.0;
                    const double dTdphi = (m.t.x[i][j][kp] - m.t.x[i][j][km1])
                                        * m.t_0 / (2.0 * m.dphi);
                    if (!AtomUtils::is_finite_safe(dTdphi)) return 0.0;
                    // v is meridional SOUTH-positive here, so it is minus the northward
                    // geostrophic component dv_n/dz = +(g/(f*T)) * (1/(a sin(theta))) dT/dphi
                    return -m.g * dTdphi / (a * sth * f * T);
                };

                const double w0 = m.w.x[i0][j][k];               // anchor: prescribed surface
                const double v0 = m.v.x[i0][j][k];
                double acc_w = 0.0, acc_v = 0.0;                 // accumulated shear [m/s]
                double prev_w = shear_w(i0), prev_v = shear_v(i0);

                for (int i = i0 + 1; i < m.im; i++) {
                    const double dz = m.get_layer_height(i) - m.get_layer_height(i-1);
                    if (!(dz > 0.0)) continue;
                    const double cur_w = shear_w(i), cur_v = shear_v(i);
                    acc_w += 0.5 * (cur_w + prev_w) * dz;
                    acc_v += 0.5 * (cur_v + prev_v) * dz;
                    prev_w = cur_w; prev_v = cur_v;

                    double dw = s * wgt * acc_w;                 // [m/s]
                    double dv = s * wgt * acc_v;
                    if (!AtomUtils::is_finite_safe(dw)) dw = 0.0;
                    if (!AtomUtils::is_finite_safe(dv)) dv = 0.0;

                    double wnew = w0 + dw * inv_u_0;             // stored non-dimensional
                    double vnew = v0 + dv * inv_u_0;
                    const double w_lim = w_max * inv_u_0;
                    if (fabs(wnew) > w_lim) { wnew = (wnew > 0.0 ? w_lim : -w_lim); n_capped++; }
                    if (do_v && fabs(vnew) > w_lim) { vnew = (vnew > 0.0 ? w_lim : -w_lim); n_capped++; }

                    if (fabs(dw) > worst_dw) worst_dw = fabs(dw);
                    m.w.x[i][j][k] = wnew;  m.wn.x[i][j][k] = wnew;
                    if (do_v) { m.v.x[i][j][k] = vnew;  m.vn.x[i][j][k] = vnew; }
                }
                n_cols++;
            }
        }

        double wmax = 0.0;
        for (int i = 0; i < m.im; i++)
            for (int j = 0; j < m.jm; j++)
                for (int k = 0; k < m.km; k++)
                    if (!is_land(m.h, i, j, k) && fabs(m.w.x[i][j][k]) > wmax)
                        wmax = fabs(m.w.x[i][j][k]);

        std::cout << "      ATOM: thermal-wind balance strength " << std::fixed
                  << std::setprecision(3) << s
                  << "  lat_min " << std::setprecision(1) << twLatMin() << " deg"
                  << "  v-component " << (do_v ? "ON" : "off")
                  << "  applied to " << n_cols << " columns"
                  << "  largest shear added " << std::fixed << std::setprecision(2)
                  << worst_dw << " m/s"
                  << "  max|w| now " << (wmax * m.u_0) << " m/s"
                  << "  capped cells " << n_capped << std::endl;
    }

    // Strength on the thermal-wind shear. 0.0 = off (default, byte-identical); 1.0 = the full
    // shear the model's own temperature implies. A strength rather than a flag, for the same
    // reason ATM_HYDRO_PGF is one: it makes a partial arm possible if the full one is unstable.
    static double twStrength(){
        static const double v = [](){
            const char* e = getenv("ATM_TW_BALANCE"); return e ? atof(e) : 0.0; }();
        return v;
    }
    static double twLatMin(){
        static const double v = [](){
            const char* e = getenv("ATM_TW_LATMIN"); return e ? atof(e) : 15.0; }();
        return v;
    }
    static double twWmax(){
        static const double v = [](){
            const char* e = getenv("ATM_TW_WMAX"); return e ? atof(e) : 80.0; }();
        return v;
    }
    // The MERIDIONAL component, from the ZONAL temperature gradient. Default OFF even when
    // ATM_TW_BALANCE is on: v is the component balance_column_mass_flux() constrains and the
    // one Psi(ground) is built from, so enabling it changes two tracked quantities at once.
    // The jet -- and the 35-65 deg storm track that needs it -- is entirely in w.
    static bool twDoV(){
        static const bool v = [](){
            const char* e = getenv("ATM_TW_BALANCE_V"); return e && atoi(e) != 0; }();
        return v;
    }

private:
    cAtmosphereModel& m;

    void init_u(Array& u, int j)
    {
        const double ua_00  = 0.02894;
        const double ua_30  = 0.02315;
        const double ua_60  = 0.01736;
        const double ua_90  = 0.011574;

        double coeff;
        switch (j) {
            case  90: coeff =  ua_00; break;
            case  60: coeff = -ua_30; break;
            case 120: coeff = -ua_30; break;
            case  30: coeff =  ua_60; break;
            case 150: coeff =  ua_60; break;
            case   0: coeff = -ua_90; break;
            case 180: coeff = -ua_90; break;
            default:  return;
        }

        const int    tl           = m.get_tropopause_layer(j);
        const double tropo_h      = m.get_layer_height(tl);
        const double half_tropo_h = tropo_h / 3.0;
        const double inv_ascent   = 3.0 / half_tropo_h;
        const double inv_descent  = 1.0 / half_tropo_h;

        #pragma omp parallel for schedule(static)
        for (int k = 0; k < m.km; k++) {
            for (int i = 0; i < tl; i++) {
                const double h     = m.get_layer_height(i);
                const double ratio = (h < half_tropo_h)
                    ? h * inv_ascent
                    : (tropo_h - h) * inv_descent;
                u.x[i][j][k] = coeff * ratio;
            }
        }
    }

    void init_v_or_w(Array& v_or_w, int j, double coeff_trop, double coeff_sl)
    {
        const int    tl          = m.get_tropopause_layer(j);
        const double inv_tropo_h = 1.0 / m.get_layer_height(tl);

        #pragma omp parallel for schedule(static)
        for (int k = 0; k < m.km; k++) {
            double sl = coeff_sl;
            if (m.use_NASA_velocity && is_ocean_surface(m.h, 0, j, k)) {
                sl = v_or_w.x[0][j][k];
            }
            const double slope = (coeff_trop - sl) * inv_tropo_h;

            for (int i = 0; i < tl; i++) {
                v_or_w.x[i][j][k] = slope * m.get_layer_height(i) + sl;
            }
        }

        init_v_or_w_above_tropopause(v_or_w, j, coeff_trop);
    }

    void init_v_or_w_above_tropopause(Array& v_or_w, int j, double coeff)
    {
        const int tl = m.get_tropopause_layer(j);
        if (tl >= m.im - 1) return;

        const double h_top     = m.get_layer_height(m.im - 1);
        const double inv_range = 1.0 / (h_top - m.get_layer_height(tl));

        #pragma omp parallel for schedule(static)
        for (int k = 0; k < m.km; k++) {
            for (int i = tl; i < m.im; i++) {
                v_or_w.x[i][j][k] = coeff * (h_top - m.get_layer_height(i)) * inv_range;
            }
        }
    }

    void form_diagonals(Array& a, int start, int end)
    {
        const double inv_range = 1.0 / (double)(end - start);

        #pragma omp parallel for collapse(2) schedule(static)
        for (int k = 0; k < m.km; k++) {
            for (int i = 0; i < m.im; i++) {
                const double a_start = a.x[i][start][k];
                const double slope   = (a.x[i][end][k] - a_start) * inv_range;
                for (int j = start; j < end; j++) {
                    a.x[i][j][k] = slope * (double)(j - start) + a_start;
                }
            }
        }
    }
};
