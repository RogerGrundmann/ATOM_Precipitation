#pragma once

#include "cAtmosphereModel.h"
#include "Utils.h"

#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstdlib>

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
        // northern polar cell
        init_v_or_w(m.v,   0,  0.5,  0.0);                              // lat:  90   j=0
        init_v_or_w(m.v,  15,  0.5,  0.6);                              // lat:  75   j=15
        // southern polar cell
        init_v_or_w(m.v, 180,  0.5,  0.0);                              // lat: -90   j=180
        init_v_or_w(m.v, 165,  0.5,  0.6);                              // lat: -75   j=165
        // northern Ferrel cell
        init_v_or_w(m.v,  30, -0.2,  0.0);                              // lat:  60   j=30
        init_v_or_w(m.v,  45,  4.0, -1.5);                              // lat:  45   j=45
        // southern Ferrel cell
        init_v_or_w(m.v, 150, -0.2,  0.0);                              // lat: -60   j=150
        init_v_or_w(m.v, 135,  4.0, -1.5);                              // lat: -45   j=135
        // northern Hadley cell
        init_v_or_w(m.v,  60,  0.0,  0.5);                              // lat:  30   j=60
//        init_v_or_w(m.v,  75, -3.0,  3.0);                              // lat:  15   j=75
        init_v_or_w(m.v,  75, -3.0,  4.0);                              // lat:  15   j=75
        // southern Hadley cell
        init_v_or_w(m.v, 120,  0.0,  0.5);                              // lat: -30   j=120
//        init_v_or_w(m.v, 105, -3.0,  4.0);                              // lat: -15   j=105
        init_v_or_w(m.v, 105, -3.0,  3.0);                              // lat: -15   j=105

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
    void balance_column_mass_flux()
    {
        if (!massBalance()) return;
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
        std::cout << "      ATOM: column mass-flux balance applied to " << n_cols
                  << " columns, largest correction " << std::scientific << std::setprecision(3)
                  << worst << " (non-dim v)" << std::endl;
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
