#include "cHydrosphereModel.h"
#include "Utils.h"
#include "cAtmosphereModel.h"

using namespace std;
using namespace AtomUtils;

// ============================================================================
// Hydrosphere Initialization Module
// ============================================================================

#include <chrono>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <sys/stat.h>

// ============================================================================
// Physical Constants — Ekman Spiral
// ============================================================================
namespace EkmanConstants {
    constexpr double DRAG_COEFF      = 2.6e-3;                          // Drag coefficient [.]
    constexpr double EKMAN_ANGLE_DEG = 45.0;                            // Ekman surface turning angle [deg]
    constexpr double V0_COEFF        = 0.0127;                          // Surface velocity coefficient [.]
    constexpr double WIND_EPS        = 1.0e-5;                          // Min |w_wind| to avoid division by zero
    constexpr double SINTHE_EPS      = 1.0e-5;                          // Min |sinthe| to avoid division by zero
}

// ============================================================================
// Physical Constants — Salinity
// ============================================================================
namespace SalinityConstants {
    constexpr double C_P             = 999.83;                          // Pure-water density term [kg/m³]
    constexpr double BETA_P          = 0.808;                           // Haline density coefficient [.]
    constexpr double ALPHA_T         = 0.0708;                          // Thermal expansion coefficient [.]
    constexpr double ALPHA_T2        = 0.068;                           // Thermal expansion sub-coefficient [.]
    constexpr double GAMMA_T         = 0.003;                           // Haline contraction coefficient [.]
    constexpr double GAMMA_T2        = 0.012;                           // Haline contraction sub-coefficient [.]
    constexpr double SAL_MAX_PSU     = 45.0;                            // Maximum salinity clamp [psu]
    constexpr double DEPTH_TEMP_DIFF = 2.0;                             // Sea-surface to deep-ocean ΔT [°C]
}

using namespace EkmanConstants;
using namespace SalinityConstants;


// ============================================================================
// Ekman Spiral
// ============================================================================
void cHydrosphereModel::EkmanSpiral() {

    std::cout << "\n\n\n      OGCM: EkmanSpiral" << std::endl;

    // Ekman spiral: 45° turning at sea surface, rotating to 90° at base of shear layer.
    // Ocean surface speed ≈ 3 % of 10-m wind speed (V0_COEFF encodes the exact scaling).

    auto begin = std::chrono::high_resolution_clock::now();

    // ========================================================================
    // Precompute radius levels (non-dimensional, surface = 1)
    // ========================================================================
    std::vector<double> radius(im, 0.0);
    for (int i = 0; i < im; i++)
        radius[i] = rad.z[im - 1] - rad.z[i];                          // depth below surface [./.], correct for stretched and uniform grids

    const double Ekman_angle = EKMAN_ANGLE_DEG / pi180;   // [rad]
    const double CD          = DRAG_COEFF;

    // ========================================================================
    // Step 1: Surface wind field (dimensional, m/s)
    // ========================================================================
    #pragma omp parallel for collapse(2)
    for (int k = 0; k < km; k++) {
        for (int j = 0; j < jm; j++) {

            v_wind.y[j][k] = v.x[im-1][j][k] * u_0_wind;                // [m/s]
            w_wind.y[j][k] = w.x[im-1][j][k] * u_0_wind;

            if (is_land(h, im - 1, j, k)) {
                v_wind.y[j][k] = 0.0;
                w_wind.y[j][k] = 0.0;
            }
        }
    }

    // ========================================================================
    // Step 2: Surface Ekman current + vertical spiral (vector form).
    //
    // Ekman (1905): the wind-driven surface current is deflected Ekman_angle
    // (45 deg) cum sole - to the RIGHT of the wind in the northern hemisphere,
    // to the LEFT in the southern - then rotates further the same way and
    // decays exponentially with depth: magnitude V_0*exp(-a*z), direction
    // turning by a*z over the Ekman scale a = sqrt(f/2Az).
    //
    // This replaces an earlier per-latitude quadrant lookup that (a) left the
    // surface-current `angle` UNINITIALISED whenever the (w,v) wind sign pair
    // fell outside the three hard-coded cases per hemisphere - the missing case
    // being the easterly-trade quadrant, i.e. exactly the subtropics - and (b)
    // treated the two hemispheres inconsistently. The result was a spurious,
    // hemisphere-asymmetric surface current and no Southern-Hemisphere
    // subtropical gyre. The vector form below is symmetric by construction and
    // has no uninitialised path. Model frame: w = east(+), v = south(+), so the
    // northward wind/current component is -v.
    // ========================================================================
    const int j_eq = (jm - 1) / 2;

    #pragma omp parallel for
    for (int k = 1; k < km - 1; k++) {
        for (int j = 1; j < jm - 1; j++) {

            // Equator: Coriolis -> 0 makes the Ekman balance singular. Keep the
            // historical treatment: no meridional flow, zonal current copied
            // from the adjacent (already-computed) cell toward the equator.
            if (j == j_eq) {
                for (int i = 0; i < im; i++) {
                    v.x[i][j][k] = 0.0;
                    w.x[i][j][k] = w.x[i][j-1][k];
                }
                continue;
            }

            double sinthe = sin(the.z[j]);
            if (sinthe == 0.0) sinthe = SINTHE_EPS;

            const double We   = w_wind.y[j][k];                            // east  [m/s]
            const double Wn   = -v_wind.y[j][k];                           // north [m/s]
            const double U_10 = sqrt(We * We + Wn * Wn);                   // [m/s]

            // Physical Ekman scales (hemisphere-symmetric via |sin(theta)|).
            const double V_0    = V0_COEFF * U_10 / sqrt(fabs(sinthe));    // [m/s] Stewart eq. 9.16
            const double T_yz   = r_air * CD * U_10 * U_10;                // [kg/(m s^2)]
            const double f      = 2.0 * omega * fabs(sinthe);             // [1/s]
            const double Az     = pow(T_yz / (r_0_water * (V_0 > 0.0 ? V_0 : 1.0)), 2) / f; // [m^2/s]
            const double a      = sqrt(f / (2.0 * Az));                    // [1/m]
            const double DE     = sqrt(2.0 * M_PI * M_PI * Az / f) * U_10; // [m] Ekman depth
            const double DE_180 = M_PI / a;                               // [m] depth of reversed current

            // Hemisphere handedness: +1 north (deflect/rotate right = clockwise
            // = decreasing math-angle), -1 south.
            const double hsign = (j < j_eq) ? 1.0 : -1.0;

            // Surface current direction = wind direction deflected by Ekman_angle.
            const double phi0 = atan2(Wn, We) - hsign * Ekman_angle;

            const bool kill = (U_10 < WIND_EPS) || (DE < DE_180);

            for (int i = 0; i < im; i++) {
                if (kill || is_land(h, i, j, k)) {
                    u.x[i][j][k] = 0.0;
                    v.x[i][j][k] = 0.0;
                    w.x[i][j][k] = 0.0;
                    continue;
                }
                const double depth = radius[i] * L_hyd;                   // [m] >= 0 below surface
                const double mag   = V_0 * exp(-a * depth);               // [m/s]
                const double ang   = phi0 - hsign * a * depth;            // spiral turning with depth
                w.x[i][j][k] =  mag * cos(ang);                           // east
                v.x[i][j][k] = -mag * sin(ang);                           // south = -north
            }
        }  // j
    }  // k

    // ========================================================================
    // Step 3: Polar boundary conditions (extrapolation, j = 0 and j = jm-1)
    // ========================================================================
    #pragma omp parallel for collapse(2)
    for (int k = 0; k < km; k++) {
        for (int i = 0; i < im; i++) {
            v.x[i][0][k]     = c43 * v.x[i][1][k]     - c13 * v.x[i][2][k];
            w.x[i][0][k]     = c43 * w.x[i][1][k]     - c13 * w.x[i][2][k];
            v.x[i][jm-1][k]  = c43 * v.x[i][jm-2][k]  - c13 * v.x[i][jm-3][k];
            w.x[i][jm-1][k]  = c43 * w.x[i][jm-2][k]  - c13 * w.x[i][jm-3][k];
/*
            v.x[i][0][k] = v.x[i][3][k]
                - 3.0 * v.x[i][2][k] + 3.0 * v.x[i][1][k];              // extrapolation
            w.x[i][jm-1][k] = w.x[i][3][k]
                - 3.0 * w.x[i][2][k] + 3.0 * w.x[i][1][k];              // extrapolation

            v.x[i][jm-1][k] = v.x[i][jm-4][k]
                - 3.0 * v.x[i][jm-3][k] + 3.0 * v.x[i][jm-2][k];        // extrapolation
            w.x[i][jm-1][k] = w.x[i][jm-4][k]
                - 3.0 * w.x[i][jm-3][k] + 3.0 * w.x[i][jm-2][k];        // extrapolation
*/
        }
    }

    // ========================================================================
    // Step 4: Longitudinal boundary conditions — periodic wrap at Greenwich
    // ========================================================================
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < im; i++) {
        for (int j = 0; j < jm; j++) {
            v.x[i][j][0]    = c43 * v.x[i][j][1]    - c13 * v.x[i][j][2];
            w.x[i][j][0]    = c43 * w.x[i][j][1]    - c13 * w.x[i][j][2];
            v.x[i][j][km-1] = c43 * v.x[i][j][km-2] - c13 * v.x[i][j][km-3];
            w.x[i][j][km-1] = c43 * w.x[i][j][km-2] - c13 * w.x[i][j][km-3];

            // Average to enforce periodicity
            v.x[i][j][0] = v.x[i][j][km-1]
                = 0.5 * (v.x[i][j][0] + v.x[i][j][km-1]);
            w.x[i][j][0] = w.x[i][j][km-1]
                = 0.5 * (w.x[i][j][0] + w.x[i][j][km-1]);
        }
    }

    // ========================================================================
    // Step 5: Vertical u-velocity from continuity (Ekman pumping)
    // ========================================================================
    {
        const int    i_max   = im - 1;
        const double dr_phys = L_hyd / i_max;                           // [m] uniform depth step

        // Surface boundary condition: no flux through the ocean surface
        #pragma omp parallel for collapse(2)
        for (int j = 0; j < jm; j++)
            for (int k = 0; k < km; k++)
                u.x[0][j][k] = 0.0;

        // Integrate continuity equation downward, one level at a time.
        // i-loop is serial: u.x[i] depends on u.x[i-1].
        // j-loop is OpenMP-parallel: each column (j,k) is independent at fixed i.
        for (int i = 1; i < im; i++) {

            const double depth        = i * dr_phys;                    // [m] depth at level i
            const double r_surface_eq = r_Earth * 1.0e3 + depth;        // [m] radius at equator and at any latitude, level i

            #pragma omp parallel for schedule(dynamic)
            for (int j = 1; j < jm-1; j++) {

                const double costhe       = cos(the.z[90] - the.z[j]);  // sin(latitude)
                const double U_Earth      = r_surface_eq * costhe * 2.0 * M_PI; // [m]
                const double inv_step_the = 0.5 * (jm-1) / U_Earth;     // [1/m]
                const double inv_step_phi = 0.5 * (km-1) / U_Earth;     // [1/m]

                for (int k = 1; k < km; k++) {

                    // Land point: zero and skip
                    if (is_land(h, i, j, k)) {
                        u.x[i][j][k] = 0.0;
                        continue;
                    }

                    // ── dvdthe: surface divergence, same for all i ────────────
                    double dvdthe;
                    bool the_flag = false;
                    if (j > 2 && j < jm - 3) {
                        if (is_land(h, i_max, j, k)
                            && is_water(h, i_max, j+1, k)
                            && is_water(h, i_max, j+2, k)) {            // pointing south
                            dvdthe = (- 3.0 * v.x[i_max][j  ][k]
                                      + 4.0 * v.x[i_max][j+1][k]
                                      -       v.x[i_max][j+2][k]) * inv_step_the; // [1/s]
                            the_flag = true;
                        } else if (is_land(h, i_max, j, k)
                                   && is_water(h, i_max, j-1, k)
                                   && is_water(h, i_max, j-2, k)) {     // pointing north
                            dvdthe = -(- 3.0 * v.x[i_max][j  ][k]
                                       + 4.0 * v.x[i_max][j-1][k]
                                       -       v.x[i_max][j-2][k]) * inv_step_the; // [1/s]
                            the_flag = true;
                        }
                    }
                    if (!the_flag)
                        dvdthe = (v.x[i_max][j+1][k] - v.x[i_max][j-1][k]) * inv_step_the; // [1/s]

                    // ── dwdphi: surface divergence, same for all i ────────────
                    double dwdphi;
                    bool phi_flag = false;
                    if (k > 2 && k < km - 3) {
                        if (is_land(h, i_max, j, k)
                            && is_water(h, i_max, j, k+1)
                            && is_water(h, i_max, j, k+2)) {            // pointing east
                            dwdphi = (- 3.0 * w.x[i_max][j][k  ]
                                      + 4.0 * w.x[i_max][j][k+1]
                                      -       w.x[i_max][j][k+2]) * inv_step_phi; // [1/s]
                            phi_flag = true;
                        } else if (is_land(h, i_max, j, k)
                                   && is_water(h, i_max, j, k-1)
                                   && is_water(h, i_max, j, k-2)) {     // pointing west
                            dwdphi = -(- 3.0 * w.x[i_max][j][k  ]
                                       + 4.0 * w.x[i_max][j][k-1]
                                       -       w.x[i_max][j][k-2]) * inv_step_phi; // [1/s]
                            phi_flag = true;
                        }
                    } else if (k == km - 1) {                           // Greenwich east boundary
                        dwdphi = -(- 3.0 * w.x[i_max][j][km-1]
                                   + 4.0 * w.x[i_max][j][km-2]
                                   -       w.x[i_max][j][km-3]) * inv_step_phi;  // [1/s]
                        phi_flag = true;
                    }
                    if (!phi_flag)
                        dwdphi = (w.x[i_max][j][k+1] - w.x[i_max][j][k-1]) * inv_step_phi; // [1/s]

                    // ── integrate u downward ──────────────────────────────────
                    u.x[i][j][k] = u.x[i-1][j][k]
                        - dr_phys * (dvdthe + dwdphi);                  // [m/s] continuity equation

                }  // end k
            }  // end j (OMP)
        }  // end i

        // Equator: zero u
        const int j_eq = (jm-1) / 2;
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < im; i++)
            for (int k = 0; k < km; k++)
                u.x[i][j_eq][k] = 0.0;
    }


    // ========================================================================
    // Step 6: Non-dimensionalise velocities; zero land points
    // ========================================================================
    #pragma omp parallel for collapse(3)
    for (int k = 0; k < km; k++) {
        for (int j = 0; j < jm; j++) {
            for (int i = 0; i < im; i++) {
                u.x[i][j][k] /= u_0;
                v.x[i][j][k] /= u_0;
                w.x[i][j][k] /= u_0;

                if (is_land(h, i, j, k)) {
                    u.x[i][j][k] = 0.0;
                    v.x[i][j][k] = 0.0;
                    w.x[i][j][k] = 0.0;
                }
            }
        }
    }

    auto end     = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" Time measured: %.3f seconds for EkmanSpiral\n", elapsed.count() * 1e-9);

    std::cout << "      OGCM: EkmanSpiral ended" << std::endl;
}


// ============================================================================
// Sverdrup / Stommel wind-driven gyre initial condition
// ============================================================================
//
// A more physically complete ocean IC than the local Ekman spiral. The observed
// surface currents are dominated by the wind-stress-CURL-driven geostrophic
// gyres (Sverdrup interior + western-boundary intensification), NOT by the thin
// ~2 %-of-wind ageostrophic Ekman drift. This routine seeds those gyres so the
// wind-stress body force (RHS_Hyd_Turb.cpp) spins the ocean up toward its own
// attractor from a state that already resembles the target circulation, instead
// of first unwinding an analytic spiral that disagrees with the discrete
// dynamics.
//
// Method (barotropic transport streamfunction Psi, per latitude row):
//   tau        = r_air * C_D * |U_wind| * U_wind            wind stress [N/m^2]
//   curl_z tau = (1/(R sin th)) [ d(tau_e sin th)/d th - d tau_s/d phi ]
//   V_north    = curl_z tau / (r_0_water * beta)            Sverdrup transport [m^2/s]
//   Psi(x)     = integral_{x_E}^{x} V_north dx'             from each basin's eastern wall
//   U_east     = -dPsi/dy = (1/R) dPsi/d th                 zonal transport [m^2/s]
//   (u,v,w)    = transports / H_ref                         barotropic velocities [m/s]
//
// The eastern-boundary integration keys off the land mask h alone, so it adapts
// to whatever basins a given Ma bathymetry provides (multiple basins, arbitrary
// coastlines). Circumpolar latitudes (no meridional walls, e.g. a Drake-open
// Southern Ocean) have no Sverdrup eastern boundary; they are left to the wind-
// stress body force to spin up and are zeroed here. The vertical velocity u is
// diagnostic and is set to zero by the caller, so it is not computed here.
namespace SverdrupConstants {
    constexpr double DRAG_COEFF    = 2.6e-3;    // bulk drag coefficient [.]
    constexpr double H_REF         = 4000.0;    // reference barotropic ocean depth [m]
    constexpr double V_MAX         = 1.5;       // cap on |velocity| component [m/s]
    constexpr double SINTHE_FLOOR  = 0.26;      // floor on sin(colat) (~poleward of 75 deg lat) — guards curl 1/sin th and beta->0
    constexpr double LAT_TAPER_HI  = 80.0;      // taper velocities to 0 poleward of this latitude [deg]
    constexpr double LAT_TAPER_LO  = 3.0;       // taper to 0 within this latitude of the equator [deg]
}

void cHydrosphereModel::SverdrupGyre() {

    // Explicit qualification (not `using namespace`): DRAG_COEFF also lives in
    // EkmanConstants, which is `using`-imported at file scope.
    using SverdrupConstants::H_REF;
    using SverdrupConstants::V_MAX;
    using SverdrupConstants::SINTHE_FLOOR;
    using SverdrupConstants::LAT_TAPER_HI;
    using SverdrupConstants::LAT_TAPER_LO;

    std::cout << "\n\n\n      OGCM: SverdrupGyre" << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    const double R      = r_Earth * 1.0e3;                 // Earth radius [m]
    const double CD     = SverdrupConstants::DRAG_COEFF;
    const int    j_eq   = (jm - 1) / 2;

    // ========================================================================
    // Step 1: Surface wind field (dimensional, m/s) — identical to EkmanSpiral.
    // These 2D fields are also consumed by the sustained wind-stress body force
    // in RHS_Hyd_Turb.cpp, so they must be populated here.
    // ========================================================================
    #pragma omp parallel for collapse(2)
    for (int k = 0; k < km; k++) {
        for (int j = 0; j < jm; j++) {
            v_wind.y[j][k] = v.x[im-1][j][k] * u_0_wind;               // [m/s]
            w_wind.y[j][k] = w.x[im-1][j][k] * u_0_wind;
            if (is_land(h, im - 1, j, k)) {
                v_wind.y[j][k] = 0.0;
                w_wind.y[j][k] = 0.0;
            }
        }
    }

    // ========================================================================
    // Step 2: Wind stress on the surface (j,k) grid.
    //   model frame: w = east(+), v = south(+)
    //   tau_east  along +phi (east),  tau_south along +theta (increasing colat)
    // Zero over land.
    // ========================================================================
    std::vector<std::vector<double>> tau_e(jm, std::vector<double>(km, 0.0));
    std::vector<std::vector<double>> tau_s(jm, std::vector<double>(km, 0.0));

    #pragma omp parallel for collapse(2)
    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {
            if (is_land(h, im - 1, j, k)) continue;
            const double We = w_wind.y[j][k];                          // east  [m/s]
            const double Ws = v_wind.y[j][k];                          // south [m/s]
            const double U  = sqrt(We * We + Ws * Ws);                 // [m/s]
            tau_e[j][k] = r_air * CD * U * We;                         // [N/m^2]
            tau_s[j][k] = r_air * CD * U * Ws;                         // [N/m^2]
        }
    }

    // ========================================================================
    // Step 3: Vertical component of the wind-stress curl and the Sverdrup
    // meridional (northward) transport V_north = curl_z(tau)/(rho*beta) [m^2/s].
    //   curl_z = (1/(R sin th)) [ d(tau_e sin th)/d th  -  d tau_s/d phi ]
    //   beta   = df/dy = 2 Omega sin(th)/R   (th = colatitude)
    // ========================================================================
    std::vector<std::vector<double>> Vn(jm, std::vector<double>(km, 0.0));

    #pragma omp parallel for
    for (int j = 1; j < jm - 1; j++) {

        double sinthe = sin(the.z[j]);
        double sinthe_g = (fabs(sinthe) < SINTHE_FLOOR)
                        ? (sinthe < 0.0 ? -SINTHE_FLOOR : SINTHE_FLOOR) : sinthe;

        const double sin_jp = sin(the.z[j+1]);
        const double sin_jm = sin(the.z[j-1]);
        const double beta   = 2.0 * omega * fabs(sinthe_g) / R;        // [1/(m s)]

        for (int k = 0; k < km; k++) {
            if (is_land(h, im - 1, j, k)) continue;

            const int kp = (k + 1) % km;
            const int km_ = (k - 1 + km) % km;

            const double dtau_e_dthe = (tau_e[j+1][k] * sin_jp
                                      - tau_e[j-1][k] * sin_jm) / (2.0 * dthe);
            const double dtau_s_dphi = (tau_s[j][kp] - tau_s[j][km_]) / (2.0 * dphi);

            const double curl_z = (dtau_e_dthe - dtau_s_dphi) / (R * sinthe_g); // [N/m^3]
            Vn[j][k] = curl_z / (r_0_water * beta);                    // [m^2/s] northward
        }
    }

    // ========================================================================
    // Step 4: Barotropic transport streamfunction Psi [m^3/s] by zonal
    // integration of dPsi/dx = V_north, westward from each basin's eastern wall
    // (Psi = 0 there). Marching starts at a land cell so basins are traversed
    // whole and the Greenwich seam is handled naturally. All-water (circumpolar)
    // rows have no eastern wall — left at Psi = 0.
    // ========================================================================
    std::vector<std::vector<double>> Psi(jm, std::vector<double>(km, 0.0));

    #pragma omp parallel for
    for (int j = 1; j < jm - 1; j++) {

        const double dx = R * std::max(fabs(sin(the.z[j])), SINTHE_FLOOR) * dphi; // [m] zonal step

        // locate a land cell to anchor the westward march
        int k_land = -1;
        for (int k = 0; k < km; k++)
            if (is_land(h, im - 1, j, k)) { k_land = k; break; }
        if (k_land < 0) continue;                                      // circumpolar row: no eastern wall

        bool prev_land = true;
        int  k_east    = -1;
        for (int s = 0; s < km; s++) {
            const int k = (k_land - s % km + km) % km;                 // march west from land
            if (is_land(h, im - 1, j, k)) {
                Psi[j][k] = 0.0;
                prev_land = true;
                continue;
            }
            if (prev_land) {
                Psi[j][k] = 0.0;                                       // eastern boundary of a basin
            } else {
                Psi[j][k] = Psi[j][k_east] - 0.5 * (Vn[j][k] + Vn[j][k_east]) * dx;
            }
            k_east    = k;
            prev_land = false;
        }
    }

    // ========================================================================
    // Step 5: Transports -> barotropic velocities, assigned over the whole water
    // column. w (east) = U_east/H_ref, v (south) = -V_north/H_ref.
    //   U_east = -dPsi/dy = (1/R) dPsi/d th  (central in colatitude)
    // Latitude taper zeros the equatorial (geostrophically singular) band and the
    // high latitudes. Velocity components are capped at V_MAX.
    // ========================================================================
    // start from zero so land and untouched cells are clean
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < im; i++)
        for (int j = 0; j < jm; j++)
            for (int k = 0; k < km; k++) {
                v.x[i][j][k] = 0.0;
                w.x[i][j][k] = 0.0;
            }

    #pragma omp parallel for
    for (int j = 2; j < jm - 2; j++) {

        const double lat_deg = (the.z[j_eq] - the.z[j]) * pi180;      // + north, - south
        const double alat    = fabs(lat_deg);

        // smooth taper: 0 in |lat|<LO, ramp to 1, then ramp back to 0 beyond HI
        double w_lat = 1.0;
        if (alat < LAT_TAPER_LO)      w_lat = 0.0;
        else if (alat < LAT_TAPER_LO + 2.0)
            w_lat = (alat - LAT_TAPER_LO) / 2.0;
        if (alat > LAT_TAPER_HI)      w_lat = 0.0;
        else if (alat > LAT_TAPER_HI - 5.0)
            w_lat = std::min(w_lat, (LAT_TAPER_HI - alat) / 5.0);
        if (w_lat <= 0.0) continue;

        for (int k = 0; k < km; k++) {
            if (is_land(h, im - 1, j, k)) continue;

            // zonal transport from the meridional streamfunction gradient;
            // only where both meridional neighbours are ocean (interior).
            double U_east = 0.0;
            if (is_water(h, im - 1, j+1, k) && is_water(h, im - 1, j-1, k))
                U_east = (Psi[j+1][k] - Psi[j-1][k]) / (2.0 * dthe * R); // [m^2/s]

            double v_south =  -Vn[j][k] / H_REF * w_lat;               // [m/s]
            double w_east  =   U_east   / H_REF * w_lat;               // [m/s]

            v_south = std::max(-V_MAX, std::min(V_MAX, v_south));
            w_east  = std::max(-V_MAX, std::min(V_MAX, w_east));

            for (int i = 0; i < im; i++) {
                if (is_land(h, i, j, k)) continue;
                v.x[i][j][k] = v_south;                                // barotropic
                w.x[i][j][k] = w_east;
            }
        }
    }

    // ========================================================================
    // Step 6: Polar boundary conditions (extrapolation, j = 0 and j = jm-1)
    // ========================================================================
    #pragma omp parallel for collapse(2)
    for (int k = 0; k < km; k++) {
        for (int i = 0; i < im; i++) {
            v.x[i][0][k]    = c43 * v.x[i][1][k]    - c13 * v.x[i][2][k];
            w.x[i][0][k]    = c43 * w.x[i][1][k]    - c13 * w.x[i][2][k];
            v.x[i][jm-1][k] = c43 * v.x[i][jm-2][k] - c13 * v.x[i][jm-3][k];
            w.x[i][jm-1][k] = c43 * w.x[i][jm-2][k] - c13 * w.x[i][jm-3][k];
        }
    }

    // ========================================================================
    // Step 7: Longitudinal periodic wrap at Greenwich
    // ========================================================================
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < im; i++) {
        for (int j = 0; j < jm; j++) {
            v.x[i][j][0] = v.x[i][j][km-1] = 0.5 * (v.x[i][j][0] + v.x[i][j][km-1]);
            w.x[i][j][0] = w.x[i][j][km-1] = 0.5 * (w.x[i][j][0] + w.x[i][j][km-1]);
        }
    }

    // ========================================================================
    // Step 8: Non-dimensionalise horizontal velocities; zero land and vertical.
    // ========================================================================
    double vmax_diag = 0.0, wmax_diag = 0.0, psi_min = 0.0, psi_max = 0.0;
    #pragma omp parallel for collapse(3)
    for (int k = 0; k < km; k++) {
        for (int j = 0; j < jm; j++) {
            for (int i = 0; i < im; i++) {
                u.x[i][j][k] = 0.0;                                    // diagnostic; caller re-zeros too
                v.x[i][j][k] /= u_0;
                w.x[i][j][k] /= u_0;
                if (is_land(h, i, j, k)) {
                    u.x[i][j][k] = 0.0;
                    v.x[i][j][k] = 0.0;
                    w.x[i][j][k] = 0.0;
                }
            }
        }
    }
    for (int j = 0; j < jm; j++)
        for (int k = 0; k < km; k++) {
            vmax_diag = std::max(vmax_diag, fabs(v.x[im-1][j][k]) * u_0);
            wmax_diag = std::max(wmax_diag, fabs(w.x[im-1][j][k]) * u_0);
            psi_min   = std::min(psi_min, Psi[j][k]);
            psi_max   = std::max(psi_max, Psi[j][k]);
        }

    auto end     = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf("      OGCM: SverdrupGyre  Psi[Sv] min/max = %.1f / %.1f   surface |v|,|w| max = %.3f / %.3f m/s\n",
           psi_min * 1e-6, psi_max * 1e-6, vmax_diag, wmax_diag);
    printf(" Time measured: %.3f seconds for SverdrupGyre\n", elapsed.count() * 1e-9);

    std::cout << "      OGCM: SverdrupGyre ended" << std::endl;
}
/*
            if((j == 45) &&(k == 180)) cout << endl << "Ekman-Layer north" << endl
                << "   j = " << j << "   k = " << k << endl
                << "   sinthe = " << sinthe << "   sqrt(sinthe) = " << sqrt(sinthe) << endl
                << "   v_wind = " << v_wind.y[j][k]
                << "   w_wind = " << w_wind.y[j][k] << endl
                << "   v = " << v.x[im-1][j][k] * u_0_wind
                << "   w = " << w.x[im-1][j][k] * u_0_wind << endl
                << "   Ekman_angle = " << Ekman_angle
                << "   Ekman_angle_deg = " << Ekman_angle * pi180 << endl
                << "   alfa = " << alfa
                << "   alfa_deg = " << alfa * pi180 << endl
                << "   U_10 = " << U_10
                << "   V_0 = " << V_0 << endl
                << "   T_yz = " << T_yz << endl
                << "   f = " << f
                << "   Az = " << Az
                << "   a = " << a << endl
                << "   DE = " << DE
                << "   DE_180 = " << DE_180
                << endl << endl;

            if((j == 135) &&(k == 180)) cout << endl << "Ekman-Layer south" << endl
                << "   j = " << j << "   k = " << k << endl
                << "   sinthe = " << sinthe << "   sqrt(sinthe) = " << sqrt(sinthe) << endl
                << "   v_wind = " << v_wind.y[j][k]
                << "   w_wind = " << w_wind.y[j][k] << endl
                << "   v = " << v.x[im-1][j][k] * u_0_wind
                << "   w = " << w.x[im-1][j][k] * u_0_wind << endl
                << "   Ekman_angle = " << Ekman_angle
                << "   Ekman_angle_deg = " << Ekman_angle * pi180 << endl
                << "   alfa = " << alfa
                << "   alfa_deg = " << alfa * pi180 << endl
                << "   U_10 = " << U_10
                << "   V_0 = " << V_0 << endl
                << "   T_yz = " << T_yz << endl
                << "   f = " << f
                << "   Az = " << Az
                << "   a = " << a << endl
                << "   DE = " << DE
                << "   DE_180 = " << DE_180
                << endl << endl;
*/
/*
                    if((j == 75) &&(k == 180)) cout << "north" << endl
                        << "   i = " << i << "   j = " << j << "   k = " << k  << endl
                        << "   rad = " << radius[i] << endl
                        << "   sinthe = " << sinthe << "   sqrt(sinthe) = " << sqrt(sinthe) << endl
                        << "   a = " << a << "   a_z = " << a_z << "   exp_a_z = " << exp_a_z
                        << "   sin_a_z = " << sin_a_z << "   cos_a_z = " << cos_a_z << endl
                        << "   alfa = " << alfa
                        << "   angle = " << angle << endl
                        << "   alfa_deg = " << alfa * pi180
                        << "   angle_deg = " << angle * pi180 << endl
                        << "   T_yz = " << T_yz
                        << "   f = " << f
                        << "   Az = " << Az
                        << "   a = " << a << endl
                        << "   U_10 = " << U_10
                        << "   V_0 = " << V_0 << endl
                        << "   v_wind = " << v_wind.y[j][k]
                        << "   w_wind = " << w_wind.y[j][k] << endl
                        << "   u = " << u.x[i][j][k]
                        << "   v = " << v.x[i][j][k]
                        << "   w = " << w.x[i][j][k] << endl << endl;
*/
/*
                    if((j == 105) &&(k == 180)) cout << "south" << endl
                        << "   i = " << i << "   j = " << j << "   k = " << k  << endl
                        << "   rad = " << radius[i] << endl
                        << "   sinthe = " << sinthe << "   sqrt(sinthe) = " << sqrt(sinthe) << endl
                        << "   a = " << a << "   a_z = " << a_z << "   exp_a_z = " << exp_a_z
                        << "   sin_a_z = " << sin_a_z << "   cos_a_z = " << cos_a_z << endl
                        << "   alfa = " << alfa
                        << "   angle = " << angle << endl
                        << "   alfa_deg = " << alfa * pi180
                        << "   angle_deg = " << angle * pi180 << endl
                        << "   T_yz = " << T_yz
                        << "   f = " << f
                        << "   Az = " << Az
                        << "   a = " << a << endl
                        << "   U_10 = " << U_10
                        << "   V_0 = " << V_0 << endl
                        << "   v_wind = " << v_wind.y[j][k]
                        << "   w_wind = " << w_wind.y[j][k] << endl
                        << "   v = " << v.x[i][j][k]
                        << "   w = " << w.x[i][j][k] << endl << endl;
*/

// ============================================================================
// Ocean Temperature Initialization
// ============================================================================
void cHydrosphereModel::initTemperature(int Ma) {

    std::cout << "\n\n\n      OGCM: initTemperature" << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    // ========================================================================
    // Temperature Variables
    // ========================================================================
    double t_paleo_add      = 0.0;
    double t_pole_add       = 0.0;
    double t_global_mean_exp = 0.0;
    double t_average        = 0.0;

    // ========================================================================
    // Step 1: Fix NASA temperature artifact at 180°E (first time slice only)
    // ========================================================================
    if (is_first_time_slice()) {
        const int k_half = (km - 1) / 2;

        #pragma omp parallel for
        for (int j = 0; j < jm; j++) {
            t.x[im-1][j][k_half] =
                0.5 * (t.x[im-1][j][k_half + 1] + t.x[im-1][j][k_half - 1]);
            temperature_NASA.y[j][k_half] =
                0.5 * (temperature_NASA.y[j][k_half + 1] +
                       temperature_NASA.y[j][k_half - 1]);
        }
    }

    // ========================================================================
    // Step 2: Apply EarthByte reconstruction temperature (non-first slices)
    // ========================================================================
    if (Ma != 0 && use_earthbyte_reconstruction) {
        #pragma omp parallel for collapse(2)
        for (int k = 0; k < km; k++) {
            for (int j = 0; j < jm; j++) {
                t.x[im-1][j][k] = (t.x[im-1][j][k] + t_0) / t_0;        // non-dimensional
                temp_reconst.y[j][k] = t.x[im-1][j][k];
            }
        }
    }

    if (is_first_time_slice() && use_earthbyte_reconstruction && !use_NASA_temperature) {
        #pragma omp parallel for collapse(2)
        for (int k = 0; k < km; k++) {
            for (int j = 0; j < jm; j++) {
                t.x[im-1][j][k] = (temperature_NASA.y[j][k] + t_0) / t_0;
            }
        }
    }

    // ========================================================================
    // Step 3: Extract temperatures from paleo curves
    // ========================================================================
    t_average        = get_temperatures_from_curve(0, m_equat_temperature_curve);
    t_pole_modern    = get_temperatures_from_curve(0, m_pole_temperature_curve);
    t_global_mean_exp = get_temperatures_from_curve(*get_current_time(),
                            m_global_temperature_curve);

    if (is_first_time_slice()) {
        t_global_mean_exp = get_temperatures_from_curve(0, m_global_temperature_curve);
        t_global_mean     = GetMean_2D(jm, km, temperature_NASA);
    }

    // ========================================================================
    // Step 4: Temperature increments between adjacent time steps
    // ========================================================================
    if (!is_first_time_slice()) {
        t_paleo_add =                                                   // equatorial increment
            get_temperatures_from_curve(*get_current_time(),  m_equat_temperature_curve)
          - get_temperatures_from_curve(*get_previous_time(), m_equat_temperature_curve);

        t_pole_add =                                                    // polar increment
            get_temperatures_from_curve(*get_current_time(),  m_pole_temperature_curve)
          - get_temperatures_from_curve(*get_previous_time(), m_pole_temperature_curve);

        t_paleo_add /= t_0;                                             // non-dimensional
        t_pole_add  /= t_0;
    }

    t_paleo_total += t_paleo_add;
    t_pole_total  += t_pole_add;

    // ========================================================================
    // Step 5: Print diagnostics
    // ========================================================================
    std::cout.precision(3);
    std::cout << "\n       Time slice of Paleo-OGCM: ...................... Ma = " << Ma << " million years\n"
              << "\n       Equatorial temperature increase: ................ t equat increase   = "
                  << t_paleo_add * t_0 << " °C"
              << "\n       Polar temperature increase: ..................... t pole increase    = "
                  << t_pole_add * t_0 << " °C"
              << "\n       Equatorial temperature at paleo times: .......... t paleo            = "
                  << ((t_average + t_0) + t_paleo_total * t_0) - t_0 << " °C"
              << "\n       Polar temperature at paleo times: ............... t pole             = "
                  << ((t_pole_modern + t_0) + t_pole_total * t_0) - t_0 << " °C"
              << "\n       Mean temperature at paleo times: ................ t global mean      = "
                  << t_global_mean << " °C"
              << "\n       Expected mean temperature at paleo times: ....... t global mean exp  = "
                  << t_global_mean_exp << " °C"
              << "\n       Equatorial temperature at modern times: ......... t modern equat     = "
                  << t_average << " °C"
              << "\n       Polar temperature at modern times: .............. t modern pole      = "
                  << t_pole_modern << " °C\n\n";

    // ========================================================================
    // Step 6: Latitudinal surface temperature distribution
    // ========================================================================
    const double d_j_half = 0.5 * (jm - 1);

    const double t_equator = (get_temperatures_from_curve(*get_current_time(),
                                  m_equat_temperature_curve) + t_0) / t_0;
    const double t_pole    = (get_temperatures_from_curve(*get_current_time(),
                                  m_pole_temperature_curve)  + t_0) / t_0;
    const double t_eff     = t_pole - t_equator;

    #pragma omp parallel for collapse(2)
    for (int k = 0; k < km; k++) {
        for (int j = 0; j < jm; j++) {
            const double ratio = (double)j / d_j_half;

            if (!use_NASA_temperature && !use_earthbyte_reconstruction) {
                // Parabolic pole-to-pole distribution
                t.x[im-1][j][k] = t_eff * parabola(ratio) + t_pole;
            }

            if (use_NASA_temperature && use_earthbyte_reconstruction) {
                if (*get_current_time() == 0) {
                    t.x[im-1][j][k] = (temperature_NASA.y[j][k] + t_0) / t_0;
                } else {
                    t.x[im-1][j][k] = temp_reconst.y[j][k];
                }
                if (*get_current_time() >= Ma_switch) {
                    t.x[im-1][j][k] = t_eff * parabola(ratio) + t_pole;
                }
            }

            temp_landscape.y[j][k] = t.x[im-1][j][k] * t_0 - t_0;       // [°C]

            if (t.x[im-1][j][k] <= t_pole_salt)
                t.x[im-1][j][k] = t_pole_salt;
        }
    }

    // ========================================================================
    // Step 7: Vertical temperature profile
    // ========================================================================
    // Shallow mode (200 m):
    //   Single linear zone — surface to (surface − 2 °C) at the ocean floor.
    //
    // Deep mode (6000 m) — three-zone piecewise profile with C¹ continuity
    //   at every zone boundary (no kinks):
    //
    //   Zone 1 [  0 – 200 m]  Linear: T_surface → T_surface − 2 °C.
    //                          Slope at 200 m = −(DEPTH_TEMP_DIFF/t_0) / z_s.
    //
    //   Zone 2 [200 – 1000 m] Cubic Hermite: T_surface−2 °C → 4 °C.
    //                          Matches zone 1 slope at 200 m; arrives at 4 °C
    //                          with zero slope at 1000 m (smooth hand-off to zone 3).
    //
    //   Zone 3 [1000 m – floor] Cubic Hermite: 4 °C → 2 °C (Antarctic Bottom Water).
    //                          Zero slope at both ends; at high latitudes the deep
    //                          water is warmer than the surface, which is physically
    //                          correct for thermohaline overturning.

    const int    i_max   = im - 1;
    const double d_i_max = (double)i_max;

    // Zone boundary depths [m]
    const double z_s = 200.0;                                          // shallow / thermocline
    const double z_t = 1000.0;                                         // thermocline / deep water

    // Anchor temperatures [non-dimensional]
    const double t_4C = (4.0 + t_0) / t_0;                            // thermocline base
    const double t_2C = (2.0 + t_0) / t_0;                            // Antarctic Bottom Water

    #pragma omp parallel for collapse(2)
    for (int k = 0; k < km; k++) {
        for (int j = 0; j < jm; j++) {

            const double t_surf = t.x[im-1][j][k];

            if (ocean_depth_mode != "deep") {

                // ── Shallow: single linear zone ───────────────────────────
                const double t_floor = t_surf - DEPTH_TEMP_DIFF / t_0;
                for (int i = 0; i < im; i++) {
                    t.x[i][j][k] = (t_surf - t_floor) * (double)i / d_i_max + t_floor;
                }

            } else {

                // ── Deep: three-zone piecewise Hermite ────────────────────
                // Temperature and slope at the zone 1/2 boundary (200 m)
                const double t_s = t_surf - DEPTH_TEMP_DIFF / t_0;    // [non-dim]
                const double m_s = -(DEPTH_TEMP_DIFF / t_0) / z_s;    // [non-dim / m]

                for (int i = 0; i < im; i++) {
                    const double depth_m = (rad.z[im-1] - rad.z[i]) * L_hyd; // [m]
                    double t_val;

                    if (depth_m <= z_s) {
                        // Zone 1 — linear (0–200 m)
                        t_val = t_surf - (DEPTH_TEMP_DIFF / t_0) * depth_m / z_s;

                    } else if (depth_m <= z_t) {
                        // Zone 2 — cubic Hermite (200–1000 m)
                        // BC: T=t_s, dT/dz=m_s at z_s;  T=t_4C, dT/dz=0 at z_t
                        const double dz = z_t - z_s;
                        const double s  = (depth_m - z_s) / dz;
                        const double s2 = s * s, s3 = s2 * s;
                        t_val = ( 2.0*s3 - 3.0*s2 + 1.0) * t_s
                              + (     s3 - 2.0*s2 + s  ) * dz * m_s
                              + (-2.0*s3 + 3.0*s2      ) * t_4C;
                                                                        // h11 term omitted (m1 = 0)
                    } else {
                        // Zone 3 — cubic Hermite (1000 m → ocean floor)
                        // BC: T=t_4C, dT/dz=0 at z_t;  T=t_2C, dT/dz=0 at L_hyd
                        const double dz = L_hyd - z_t;
                        const double s  = (depth_m - z_t) / dz;
                        const double s2 = s * s, s3 = s2 * s;
                        t_val = ( 2.0*s3 - 3.0*s2 + 1.0) * t_4C
                              + (-2.0*s3 + 3.0*s2      ) * t_2C;
                                                                        // h10, h11 omitted (m0 = m1 = 0)
                    }
                    t.x[i][j][k] = t_val;
                }
            }
        }
    }

    auto end     = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" Time measured: %.3f seconds for initTemperature\n", elapsed.count() * 1e-9);

    std::cout << "      OGCM: initTemperature ended" << std::endl;
}
/*
*
*/
// ============================================================================
// Ocean Salinity Initialization
// ============================================================================
void cHydrosphereModel::initSalinity() {

    std::cout << "\n\n\n      OGCM: initSalinity" << std::endl;
    auto begin = std::chrono::high_resolution_clock::now();

    // ========================================================================
    // Step 1: Paleo salinity correction
    // ========================================================================
    // Empirical salinity–temperature relation (Ocean Circulation, Open University):
    //   S [psu] = (T_global_mean [°C] + 346) / 10
    // t_paleo_add [°C] is the change in global mean temperature between time steps;
    // it is read directly from the curve and is NOT non-dimensionalised here.

    // Paleo time slice. NOTE: initSalinity takes no Ma argument, so the bare
    // member `Ma` is NOT the running slice (RunTimeSlice sets m_current_time, not
    // the member) — it reads 0 at every slice, which silently disabled the paleo
    // offset and the Ma>0 SSS branch below. Use get_current_time() as
    // initTemperature does.
    const int Ma_slice = (int)*get_current_time();

    const double t_global_mean_exp = get_temperatures_from_curve(*get_current_time(),
                                         m_global_temperature_curve);      // [°C]

    double t_paleo_add = 0.0;
    if (Ma_slice > 0) {
        // Sequential NASA run: increment relative to the previous slice. But a
        // single-Ma paleo run (the supported workflow) has no previous slice —
        // get_previous_time() then throws — so fall back to the absolute paleo-
        // minus-modern difference, as the non-NASA branch already does. Guard
        // with is_first_time_slice() (mirrors initTemperature).
        if (use_NASA_temperature && !is_first_time_slice())
            t_paleo_add =
                get_temperatures_from_curve(*get_current_time(),  m_global_temperature_curve)
              - get_temperatures_from_curve(*get_previous_time(), m_global_temperature_curve);
        else
            t_paleo_add =
                get_temperatures_from_curve(*get_current_time(), m_global_temperature_curve)
              - t_global_mean_exp;
    }
    // t_paleo_add is in °C — do NOT multiply by t_0

    const double c_average  = (t_global_mean_exp              + 346.0) / 10.0; // [psu] modern mean
    const double c_paleo    = (Ma_slice == 0) ? 0.0
                            : (t_global_mean_exp + t_paleo_add + 346.0) / 10.0
                              - c_average;                                      // [psu] paleo offset
    const double c_paleo_nd = c_paleo / c_35;                                   // [non-dim]

    // ========================================================================
    // Step 2: Print diagnostics
    // ========================================================================
    std::cout.precision(3);
    std::cout << "\n       Temperature increase at paleo times: ............ t increase         = "
                  << t_paleo_add << " °C"
              << "\n       Salinity increase at paleo times: ............... salinity increase   = "
                  << c_paleo << " psu"
              << "\n       Mean salinity at modern times: ................... salinity modern     = "
                  << c_average << " psu"
              << "\n       Mean salinity at paleo times: .................... salinity paleo      = "
                  << c_average + c_paleo << " psu\n\n";

    // ========================================================================
    // Step 3: Sea-surface salinity
    //   use_NASA_salinity: observed NASA sea-surface salinity climatology.
    //     Ma == 0 -> the full 2-D field (modern geography matches).
    //     Ma  > 0 -> its zonal-mean latitude profile. The E-P-driven latitudinal
    //                band structure (salty subtropics, fresh equator/poles) is
    //                set by the latitude-organised atmospheric circulation and
    //                transfers across paleogeographies; the modern LONGITUDINAL
    //                pattern is tied to today's coastlines and does not.
    //   Otherwise (toggle off, or the file is missing) -> invert the Gill (1982)
    //   density equation for surface salinity (the historical synthetic IC).
    //   All branches add the paleo mean-salinity offset c_paleo_nd.
    // ========================================================================
    const int i_max = im - 1;

    bool nasa_sss = false;
    if (use_NASA_salinity) {
        struct stat info;
        if (stat(salinity_file.c_str(), &info) == 0) {
            Array_2D sss(jm, km, 0.0);
            read_IC(salinity_file, sss.y, jm, km);                          // [psu]

            if (Ma_slice == 0) {
                #pragma omp parallel for collapse(2)
                for (int k = 0; k < km; k++)
                    for (int j = 0; j < jm; j++)
                        c.x[i_max][j][k] = sss.y[j][k] / c_35 + c_paleo_nd;
            } else {
                // zonal-mean latitude profile over ocean cells
                std::vector<double> zmean(jm, 0.0);
                for (int j = 0; j < jm; j++) {
                    double sum = 0.0; int n = 0;
                    for (int k = 0; k < km; k++)
                        if (is_water(h, i_max, j, k)) { sum += sss.y[j][k]; ++n; }
                    zmean[j] = (n > 0) ? sum / (double)n : 0.0;
                }
                // fill all-land latitude rows from the nearest valid neighbour
                for (int j = 0; j < jm; j++) if (zmean[j] == 0.0)
                    for (int d = 1; d < jm; d++) {
                        if (j-d >= 0   && zmean[j-d] > 0.0) { zmean[j] = zmean[j-d]; break; }
                        if (j+d < jm   && zmean[j+d] > 0.0) { zmean[j] = zmean[j+d]; break; }
                    }
                #pragma omp parallel for collapse(2)
                for (int k = 0; k < km; k++)
                    for (int j = 0; j < jm; j++)
                        c.x[i_max][j][k] = zmean[j] / c_35 + c_paleo_nd;
            }
            nasa_sss = true;
            std::cout << "       sea-surface salinity IC: NASA "
                      << (Ma_slice == 0 ? "2-D field" : "zonal-mean latitude profile")
                      << " (+ " << c_paleo << " psu paleo offset)\n";
        } else {
            std::cout << "       use_NASA_salinity set but " << salinity_file
                      << " not found — falling back to the Gill density inversion\n";
        }
    }

    if (!nasa_sss) {
        #pragma omp parallel for collapse(2)
        for (int k = 0; k < km; k++) {
            for (int j = 0; j < jm; j++) {
                const double t_Celsius = t.x[i_max][j][k] * t_0 - t_0;      // [°C]
                const double alfa_t_p  = ALPHA_T * (1.0 + ALPHA_T2 * t_Celsius);
                const double gamma_t_p = GAMMA_T * (1.0 - GAMMA_T2 * t_Celsius);
                const double rho       = (iter_n == 1) ? r_0_saltwater
                                                       : r_salt_water.x[i_max][j][k];

                c.x[i_max][j][k] =
                    ((rho - C_P + (alfa_t_p + 35.0 * gamma_t_p) * t_Celsius)
                     / (BETA_P + gamma_t_p * t_Celsius)) / c_35
                    + c_paleo_nd;                                            // apply paleo offset
            }
        }
    }

    // ========================================================================
    // Step 4: Linear salinity profile from surface to ocean floor
    // ========================================================================
    // Deep-ocean salinity is anchored at 34.7 psu (Antarctic Bottom Water /
    // North Atlantic Deep Water average), independent of the local surface value.
    // Where the surface is saltier (tropics), salinity decreases with depth.
    // Where the surface is fresher (poles), it increases — both physically correct.
    // The paleo correction shifts the entire water column uniformly.

    const double c_floor_nd = 34.7 / c_35 + c_paleo_nd;                       // [non-dim] floor anchor

    #pragma omp parallel for collapse(2)
    for (int k = 0; k < km; k++) {
        for (int j = 0; j < jm; j++) {
            const double c_surf = c.x[i_max][j][k];

            for (int i = 0; i < im; i++) {
                const double s = (double)i / (double)i_max;                   // 0 = floor, 1 = surface
                c.x[i][j][k]  = c_floor_nd * (1.0 - s) + c_surf * s;

                if (c.x[i][j][k] >= SAL_MAX_PSU / c_35)  c.x[i][j][k] = SAL_MAX_PSU / c_35;
                if (c.x[i][j][k] < 0.0)                 c.x[i][j][k] = 0.0;
                if (is_land(h, i, j, k))                 c.x[i][j][k] = 0.0;
            }
        }
    }

    auto end     = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" Time measured: %.3f seconds for initSalinity\n", elapsed.count() * 1e-9);

    std::cout << "      OGCM: initSalinity ended" << std::endl;
}
/*
*
*/
// ============================================================================
// Ocean / Land Fraction and CO2 Budget
// ============================================================================
void cHydrosphereModel::LandOceanFraction() {

    std::cout << "\n\n\n      OGCM: LandOceanFraction" << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    // ========================================================================
    // Step 1: Count ocean and land grid points at the sea surface
    // ========================================================================
    const int h_point_max = (jm - 1) * (km - 1);
    int       h_land      = 0;

    for (int j = 0; j < jm; j++)
        for (int k = 0; k < km; k++)
            if (is_land(h, 0, j, k))  h_land += h.x[0][j][k];

    const int    h_ocean    = h_point_max - h_land;
    const double ocean_land = (double)h_ocean / (double)h_land;

    // ========================================================================
    // Step 2: Print diagnostics
    // ========================================================================
    std::cout.precision(3);
    std::cout << "\n       Total number of points at constant height: ....... = " << h_point_max
              << "\n       Number of points on the ocean surface: ........... = " << h_ocean
              << "\n       Number of points on the land surface: ............. = " << h_land
              << "\n       Ocean / land ratio: ............................... = " << ocean_land
              << "\n"
              << "\n       Addition of CO2 by ocean surface: ................. = " << co2_ocean
              << "\n       Addition of CO2 by land surface: .................. = " << co2_land
              << "\n       Subtraction of CO2 by vegetation: ................. = " << co2_vegetation
              << "\n       (valid for one single point on the surface)\n\n";

    auto end     = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" Time measured: %.3f seconds for LandOceanFraction\n", elapsed.count() * 1e-9);

    std::cout << "      OGCM: LandOceanFraction ended" << std::endl;
}
