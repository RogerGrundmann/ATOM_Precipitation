/*
 * Atmosphere General Circulation Modell(AGCM) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
*/

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cstring>
#include <cfloat>
#include <cmath>
#include "cAtmosphereModel.h"
#include "Utils.h"     // AtomUtils::is_finite_safe — bit-level NaN/Inf test that survives -ffast-math

using namespace std;

namespace{
    string heading_1 = " printout of maximum and minimum values of properties at their locations: latitude, longitude, level";
    string heading_2 = " results based on three dimensional considerations of the problem";
    string level = "m";
    struct HemisphereCoords{
        double lat, lon;
        string east_or_west, north_or_south;
    };
    HemisphereCoords convert_coords(double lon, double lat){
        HemisphereCoords ret;
        if(lat > 90){
            ret.lat = lat - 90;
            ret.north_or_south = "°S";
        }else{
            ret.lat = 90 - lat;
            ret.north_or_south = "°N";
        }
        if(lon > 180){
            ret.lon = 360 - lon;
            ret.east_or_west = "°W";
        }else{
            ret.lon = lon;
            ret.east_or_west = "°E";
        }
        return ret;
    }
}
/*
*
*/
void cAtmosphereModel::searchMinMax_3D(const string &name_maxValue, const string &name_minValue,
    const string &name_unitValue, Array &value_D, double coeff,
    std::function< double(double) > lambda, bool print_heading){
    double maxValue = -DBL_MAX;
    double minValue =  DBL_MAX;
    int imax = 0, jmax = 0, kmax = 0;
    int imin = 0, jmin = 0, kmin = 0;
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            double *row = value_D.x[i][j];
            for(int k = 0; k < km; k++){
                double val = row[k];
                if(!std::isfinite(val)) continue;
                if(val > maxValue){
                    maxValue = val;
                    imax = i; jmax = j; kmax = k;
                }
                if(val < minValue){
                    minValue = val;
                    imin = i; jmin = j; kmin = k;
                }
            }
        }
    }
    if(maxValue == -DBL_MAX) maxValue = 0.0;
    if(minValue ==  DBL_MAX) minValue = 0.0;
 
    int imax_level = get_layer_height(imax);
    int imin_level = get_layer_height(imin);
    //  maximum latitude and longitude units recalculated
    HemisphereCoords coords = convert_coords(kmax, jmax);
    int jmax_deg = coords.lat;
    string deg_lat_max = coords.north_or_south;
    int kmax_deg = coords.lon;
    string deg_lon_max = coords.east_or_west;
    //  minimum latitude and longitude units recalculated
    coords = convert_coords(kmin, jmin);
    int jmin_deg = coords.lat;
    string deg_lat_min= coords.north_or_south;
    int kmin_deg = coords.lon;
    string deg_lon_min = coords.east_or_west;
    cout.precision(6);
    if(print_heading){
        cout << endl << heading_1 << endl << heading_2 << endl << endl;
    }
    maxValue = lambda(maxValue * coeff);
    minValue = lambda(minValue * coeff);
    cout << setiosflags(ios::left) << setw(26) << setfill('.') << name_maxValue << " = " << 
        resetiosflags(ios::left) << setw(12) << fixed << setfill(' ') << maxValue << setw(6) << 
        name_unitValue << setw(5) << jmax_deg << setw(3) << deg_lat_max << setw(4) << kmax_deg << 
        setw(3) << deg_lon_max << setw(6) << imax_level << setw(2) << level << "          " << 
        setiosflags(ios::left) << setw(26) << setfill('.') << name_minValue << " = "<< 
        resetiosflags(ios::left) << setw(12) << fixed << setfill(' ') << minValue << setw(6) << 
        name_unitValue << setw(5)  << jmin_deg << setw(3) << deg_lat_min << setw(4) << kmin_deg << 
        setw(3) << deg_lon_min  << setw(6) << imin_level << setw(2) << level << endl;
    return;
}
/*
*
*/
void cAtmosphereModel::searchMinMax_2D(const string &name_maxValue, const string &name_minValue,
    const string &name_unitValue, Array_2D &value, double coeff){
    double minValue = value.y[0][0];
    double maxValue = value.y[0][0];
    int jmax = 0, kmax = 0;
    int jmin = 0, kmin = 0;
    for(int j = 1; j < jm-1; j++){
        double *row = value.y[j];
        for(int k = 1; k < km-1; k++){
            double val = row[k];
            if(val > maxValue){
                maxValue = val;
                jmax = j; kmax = k;
            }else if(val < minValue){
                minValue = val;
                jmin = j; kmin = k;
            }
        }
    }

    int imax_level = 0;
    int imin_level = 0;
    //  maximum latitude and longitude units recalculated
    HemisphereCoords coords = convert_coords(kmax, jmax);
    int jmax_deg = coords.lat;
    string deg_lat_max = coords.north_or_south;
    int kmax_deg = coords.lon;
    string deg_lon_max = coords.east_or_west;
    //  minimum latitude and longitude units recalculated
    coords = convert_coords(kmin, jmin);
    int jmin_deg = coords.lat;
    string deg_lat_min= coords.north_or_south;
    int kmin_deg = coords.lon;
    string deg_lon_min = coords.east_or_west;
    cout.precision(6);
    maxValue = maxValue * coeff;
    minValue = minValue * coeff;
    cout << setiosflags(ios::left) << setw(26) << setfill('.') << name_maxValue << " = " <<
        resetiosflags(ios::left) << setw(12) << fixed << setfill(' ') << maxValue << setw(6) <<
        name_unitValue << setw(5) << jmax_deg << setw(3) << deg_lat_max << setw(4) << kmax_deg <<
        setw(3) << deg_lon_max << setw(6) << imax_level << setw(2) << level << "          " <<
        setiosflags(ios::left) << setw(26) << setfill('.') << name_minValue << " = "<<
        resetiosflags(ios::left) << setw(12) << fixed << setfill(' ') << minValue << setw(6) <<
        name_unitValue << setw(5)  << jmin_deg << setw(3) << deg_lat_min << setw(4) << kmin_deg <<
        setw(3) << deg_lon_min  << setw(6) << imin_level << setw(2) << level << endl;
    return;
}

/*
* Zonal-mean meridional wind [v] and meridional mass streamfunction Ψ — the standard
* diagnostic for the Hadley/Ferrel/Polar overturning cells. ATOM convention: v is the
* meridional component, POSITIVE SOUTHWARD; v.x is non-dimensional, so ×u_0 → m/s. The
* zonal mean is taken over FLUID cells only (longitudes where the cell sits above the
* local terrain, i >= i_topography[j][k]).
*
*   Ψ(i,j) = 2π·a·cosφ·∫_z^{z_top} [ρ·v]_phys dz'       [kg/s]
*
* integrated DOWNWARD from the model lid (Ψ = 0 at i = im-1). With v > 0 southward, the
* cell STRENGTH is max|Ψ|; a tropical Ψ extremum is the Hadley cell. Writes a long-format
* CSV (lat × height) per checkpoint and logs the strongest cells so collapse can be read
* straight from the run log without ParaView. See [[project_vw_drag_cut_hadley]].
*
* THE DENSITY IS INSIDE THE INTEGRAL AND INSIDE THE ZONAL MEAN, since 2026-08-26. Ported
* from ATHAD dabbc94. What was here was
*
*     const double rho = r_air;               // Boussinesq reference density
*     coeff = two_pi * a * cosphi * rho;      // ... applied at EVERY level
*
* i.e. one constant, the SURFACE density, used from the ground to the lid, on an integral
* labelled kg/s. That is a volume flux, and the meridional MASS streamfunction is the
* mass-weighted one.
*
* ON THIS SHELL THE ERROR DISTORTS RATHER THAN INVERTS, and the honest statement is that it
* is a factor of ~6 (rho ~1.2 at the surface against ~0.2 at 16 km), about 2x at the Hadley
* cell's own core near 500 hPa and worse above it. It moves the reported Psi_max and can move
* the height it is reported at; it will not manufacture a cell. Measured on landing:
* Psi_max 851.68 -> 373.84 (1e9 kg/s), 2.28x, with the location unmoved.
*
* Repaired anyway, because this diagnostic is what 24ff23a ("revive Hadley/Ferrel cells") and
* 1e59daa (jet spin-down) judge cell strength WITH, and an instrument used to decide whether
* a circulation is spinning down should not be 2x overweight at the branch it is measuring.
*
* AND IN ATHAD (dabbc94) THE SAME LINE HID A CLOSED CELL OUTRIGHT, because rho spans four
* orders of magnitude over its 300 km shell: corrected there, the peak moved from
* +45 deg / 44.9 km to +15 deg / surface and a return branch appeared at ~60 km, while what
* the old weighting had reported as the Hadley cell was a one-signed drift with NO RETURN
* FLOW. That is the difference between an instrument that is inaccurate and one that is
* lying, and which of the two you have depends on how far rho spans -- worth knowing before
* trusting this diagnostic in a new tree.
*
* EVERY Psi NUMBER RECORDED IN THIS TREE BEFORE 2026-08-26 IS THE OLD, VOLUME-FLUX-LIKE
* QUANTITY. They are not comparable with what this function prints now.
*
* Averaging rho*v zonally rather than multiplying the two zonal means also keeps the
* correlation term <rho'v'>, which is the part a warm rising branch carries.
*/

// ============================================================================================
// ATM_CELL_ROT_DIAG=1 -- WHICH WAY DO THE CELLS ACTUALLY TURN. Print-only, default off.
//
// WHY THIS EXISTS, AND IT IS A GAP THIS TREE HAS HAD ALL ALONG. The user reported (2026-09-11)
// that the polar cells turn the wrong way in a glyph plot. Nothing in this model could check
// it:
//
//   - Psi is built from the ZONAL-MEAN MERIDIONAL WIND ALONE (write_meridional_streamfunction
//     integrates rho*v and nothing else), so it is BLIND to the radial velocity by
//     construction. A cell can have a perfectly correct Psi and a vertical velocity that turns
//     the plotted arrows the other way, and every cell table in CLAUDE.md would still read
//     "correct".
//   - the only (v,u) field written is the zonal VTK slice at ONE longitude, k_zonal = 87, which
//     cuts Tibet and the Antarctic plateau. At 75S that slice has 25 of its 41 levels in rock,
//     so the southern polar cell -- the one this tree scores worst everywhere -- cannot be read
//     off it at all.
//
// So the rotation of the cells has never been measured. This measures it.
//
// WHAT IT COMPUTES. The zonal means of the RADIAL and MERIDIONAL velocity first, on the same
// masking and the same rock-contributes-zero divisor the corrected Psi uses, and then the
// rotation of that mean field in the meridional plane:
//
//     omega = d(ubar)/dy - d(vbar_north)/dz          [1/s], y NORTHWARD, z UP
//
// v is stored SOUTHWARD-positive here, so vbar_north = -vbar; and j increases SOUTHWARD, so
// d/dy = -(1/(a*Delta)) d/dj. This is the plane-Cartesian curl of the field a glyph plot draws,
// which is the quantity the eye judges -- NOT the full spherical vorticity, whose curvature
// terms are irrelevant to which way an arrow points. Said plainly so nobody quotes it as
// vorticity.
//
// HOW TO READ IT. A thermally DIRECT cell (Hadley, polar) and a thermally INDIRECT one (Ferrel)
// must have OPPOSITE omega, and the two hemispheres mirror, so the test is always against the
// HADLEY CELL OF THE SAME HEMISPHERE. Earth: 75 and 15 agree, 45 opposes. The print also gives
// the sign Psi implies, and flags DISAGREE wherever the two differ -- which is the whole point,
// because that is the case Psi cannot see.
//
// It also reports, per hemisphere, the latitude at which the zonal-mean RADIAL velocity changes
// sign between 60 and 90 degrees -- the ascent/descent boundary of the polar cell. Continuity
// against a cell whose amplitude peaks at 75 puts it near 70; init_u's linear blend between
// -ua_90 at the pole and +ua_60 at 60 puts it at 78.0, and the projection leaves it near 77.7.
void cAtmosphereModel::report_cell_rotation(int iter){
    static const bool on = [](){ const char* e = getenv("ATM_CELL_ROT_DIAG");
                                 return e && atoi(e) != 0; }();
    if(!on) return;
    using namespace std;
    const double a     = r_Earth * 1000.0;
    const double Delta = M_PI / (double)(jm - 1);         // radians of latitude per j index
    const double dy    = a * Delta;                       // metres per j index

    // Zonal means over fluid cells, rock contributing ZERO rather than being excluded -- the
    // same divisor the corrected streamfunction uses, and for the same reason: no mass crosses
    // rock, and dividing by the fluid count inflates levels that clip terrain by km/n.
    vector<vector<double> > ubar(im, vector<double>(jm, 0.0));
    vector<vector<double> > vbar(im, vector<double>(jm, 0.0));   // NORTHWARD positive
    vector<vector<double> > rvbar(im, vector<double>(jm, 0.0));  // rho*v northward, for Psi
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            double us = 0.0, vs = 0.0, rs = 0.0;
            for(int k = 0; k < km; k++){
                if(i < i_topography[j][k]) continue;
                const double uu = u.x[i][j][k], vv = v.x[i][j][k];
                double rho = r_humid.x[i][j][k];
                if(!AtomUtils::is_finite_safe(uu) || !AtomUtils::is_finite_safe(vv)) continue;
                if(!AtomUtils::is_finite_safe(rho) || rho <= 0.0) rho = r_air;
                us += uu; vs += -vv; rs += rho * (-vv);      // -v = northward
            }
            ubar[i][j]  = us / (double)km * u_0;
            vbar[i][j]  = vs / (double)km * u_0;
            rvbar[i][j] = rs / (double)km * u_0;
        }
    }
    // Psi from the same means, integrated downward from the lid, so the two senses are formed
    // from ONE field and a disagreement cannot be a difference of convention.
    vector<vector<double> > psi(im, vector<double>(jm, 0.0));
    for(int j = 0; j < jm; j++){
        const double coeff = 2.0 * M_PI * a * sin(the.z[j]);       // cos(latitude)
        for(int i = im - 2; i >= 0; i--)
            psi[i][j] = psi[i+1][j] + coeff * 0.5 * (rvbar[i][j] + rvbar[i+1][j])
                        * (get_layer_height(i+1) - get_layer_height(i));
    }

    // Function-local, NOT members: adding a member moves sizeof(cAtmosphereModel) and that is
    // this tree's stack-canary hazard (see the build-hazard note in CLAUDE.md). Same reason the
    // floor-injection counter in IceSchemeCommon is a function-local static.
    vector<double> rot_om(jm, 0.0), rot_psi(jm, 0.0), rot_a(jm, 0.0), rot_b(jm, 0.0);
    vector<int>    rot_n(jm, 0);

    // THE AVERAGING WINDOW IS FIXED IN HEIGHT, AND THAT IS NOT A DETAIL.
    //
    // The first version of this took the window from tropopause_layers[j] and from the HIGHEST
    // ground on the latitude circle. Both were wrong, and the print said so on its face by
    // reporting `levels = 0` at 75N and 75S -- the two bands it exists to measure:
    //
    //   - the highest ground anywhere on a latitude circle is the Antarctic plateau at 75S and a
    //     single mountain elsewhere, so one column deleted the whole band. It also made the
    //     level count wildly asymmetric between mirror latitudes (15N 20, 15S 9), which is a
    //     terrain artefact and not a hemispheric difference.
    //   - and the SHIPPED tropopause index at 75 deg is 23 = 2987 m, BELOW that plateau, so the
    //     window ran from 26 down to 22 and was empty. Worse, tropopause_layers is exactly what
    //     ATM_TROPO_INDEX_FIX changes: a window taken from it would differ between the two
    //     branches, and the diagnostic could not compare them at all.
    //
    // So the window is ground-to-~10 km in METRES, identical on every branch and at every
    // latitude. Terrain needs no masking here because ubar/vbar are zonal means in which rock
    // already contributes zero.
    int itop_fixed = im - 2;
    for(int i = 0; i < im - 1; i++) if(get_layer_height(i) >= 10000.0){ itop_fixed = i; break; }
    if(itop_fixed > im - 2) itop_fixed = im - 2;

    // ---- PLOT-NORMALISED, AND THE PHYSICAL CURL IS THE WRONG MEASURE ----
    //
    // The first version averaged omega = du/dy - dv/dz in METRES and it is not a discriminator
    // of cell sense at all. A tropospheric cell is ~1000:1 FLAT, so the two terms differ by
    // three orders of magnitude, their difference is whichever one is larger plus noise, and
    // every band comes out looking alike. Measured, it read -2.7e-06 at 45N -- indistinguishable
    // from zero -- and fired a DISAGREE against a Psi that was perfectly healthy (-3.95e+10).
    //
    // What the eye judges, and what uv_plot in Paraview_Atm.cpp draws, is the rotation of the
    // field AFTER the two axes are brought to comparable length. So normalise the coordinates:
    // Z' = z/H over the shell depth and Y' = y/L pole to pole, giving
    //
    //     omega' = (L/H) du/dy  -  (H/L) dv/dz
    //
    // with L/H = pi*a / h_top ~ 1249 here. This is the SHAPE-PRESERVING normalisation and it is
    // deliberately NOT tied to Paraview's dx = 0.1 / dy = 0.05 plot constants, which are
    // arbitrary and would couple this print to a writer's cosmetics; the picture's own ratio is
    // 555, a factor 2.25 from this one, which matters only where the two terms nearly cancel.
    // BOTH TERMS ARE THEREFORE PRINTED, and a band whose terms are within 3x AND of opposite
    // sign is flagged `marginal` rather than being given a confident sense.
    const double H_shell = get_layer_height(im - 1);
    const double L_merid = M_PI * a;
    const double aspect  = (H_shell > 0.0) ? (L_merid / H_shell) : 1.0;

    auto band = [&](int j)->void{
        const int itop = itop_fixed;
        double accA = 0.0, accB = 0.0; int n = 0;
        for(int i = 1; i <= itop - 1; i++){
            if(j < 1 || j > jm - 2) continue;
            const double dudy = -(ubar[i][j+1] - ubar[i][j-1]) / (2.0 * dy);   // j is southward
            const double dz   = get_layer_height(i+1) - get_layer_height(i-1);
            if(!(dz > 0.0)) continue;
            const double dvdz = (vbar[i+1][j] - vbar[i-1][j]) / dz;
            accA += aspect * dudy;          // the du/dy half, stretched to plot proportions
            accB += dvdz / aspect;          // the dv/dz half, compressed to match
            n++;
        }
        const double A = (n > 0) ? accA / n : 0.0;
        const double B = (n > 0) ? accB / n : 0.0;
        // Psi's own sense: the sign of the column extremum measured from the GROUND value, so
        // the standing Psi(ground) offset cannot flip it (the 2026-09-10 detrending lesson).
        const double pg = psi[0][j];
        double pk = 0.0;
        for(int i = 0; i < im; i++) if(fabs(psi[i][j] - pg) > fabs(pk)) pk = psi[i][j] - pg;
        rot_om[j] = A - B; rot_a[j] = A; rot_b[j] = B; rot_psi[j] = pk; rot_n[j] = n;
    };

    const int JS[6] = {15, 45, 75, 105, 135, 165};       // 75N 45N 15N 15S 45S 75S
    for(int q = 0; q < 6; q++) band(JS[q]);

    cout << "      AGCM: [CELL ROT] iter " << iter
         << " -- rotation of the PLOTTED (v,u) zonal-mean field in ASPECT-NORMALISED"
         << " coordinates, omega\' = (L/H)du/dy - (H/L)dv/dz, L/H = " << fixed
         << setprecision(0) << aspect << defaultfloat
         << ".  Sense is vs the HADLEY cell of the SAME hemisphere." << endl;
    cout << "        lat   (L/H)du/dy    (H/L)dv/dz       omega\'    sense(omega)  sense(Psi)"
         << "   Psi-Psi(gnd) e9" << endl;
    for(int q = 0; q < 6; q++){
        const int j = JS[q];
        const double lat = 90.0 - (double)j * 180.0 / (double)(jm - 1);
        const int jref = (lat > 0.0) ? 75 : 105;                  // 15N / 15S
        const bool isref = (j == jref);
        const bool so_dir = (rot_om[j]  * rot_om[jref]  > 0.0);
        const bool sp_dir = (rot_psi[j] * rot_psi[jref] > 0.0);
        // A band whose two terms nearly cancel has a sense that depends on the normalisation
        // constant rather than on the flow; say so instead of asserting one.
        const double mA = fabs(rot_a[j]), mB = fabs(rot_b[j]);
        const bool marginal = (rot_a[j] * rot_b[j] > 0.0) && mA > 0.0 && mB > 0.0
                              && (mA < 3.0 * mB) && (mB < 3.0 * mA);
        const char* so = isref ? "reference" : (marginal ? "marginal" : (so_dir ? "direct  " : "INDIRECT"));
        const char* sp = isref ? "reference" : (sp_dir ? "direct  " : "INDIRECT");
        const bool disagree = !isref && !marginal && (so_dir != sp_dir);
        cout << "      " << setw(5) << (int)lat << "  " << scientific << setprecision(3)
             << setw(12) << rot_a[j] << " " << setw(13) << rot_b[j] << " "
             << setw(12) << rot_om[j] << "   " << so << "     " << sp
             << "   " << setw(11) << rot_psi[j] / 1.0e9 << defaultfloat
             << (disagree ? "   <== DISAGREE: Psi cannot see this" : "")
             << (marginal ? "   (terms within 3x: sense depends on the normalisation)" : "")
             << endl;
    }
    // The polar ascent/descent boundary, from the zonal-mean radial velocity.
    for(int hemi = 0; hemi < 2; hemi++){
        const int ja = (hemi == 0) ? 0  : jm - 1;        // pole
        const int jb = (hemi == 0) ? 30 : jm - 31;       // 60 deg
        const int st = (hemi == 0) ? 1 : -1;
        // A FIXED reference HEIGHT, ~5 km, for the same reason the window above is fixed: taken
        // from the terrain it gave 7412 m in the north against 9007 m in the south and the two
        // hemispheres were not comparable.
        int iref = im - 3;
        for(int i = 0; i < im - 2; i++) if(get_layer_height(i) >= 5000.0){ iref = i; break; }
        double prev = ubar[iref][ja]; int jprev = ja; bool found = false;
        for(int j = ja + st; j != jb + st; j += st){
            const double cur = ubar[iref][j];
            if(prev != 0.0 && cur * prev < 0.0){
                const double f = fabs(prev) / (fabs(prev) + fabs(cur));
                const double jz = (double)jprev + st * f;
                cout << "      [CELL ROT] " << (hemi == 0 ? "N" : "S")
                     << " polar ascent/descent boundary (zonal-mean u sign change at "
                     << (int)get_layer_height(iref) << " m): "
                     << fixed << setprecision(1)
                     << fabs(90.0 - jz * 180.0 / (double)(jm - 1)) << " deg"
                     << defaultfloat
                     << "   (continuity against a cell peaking at 75 wants ~70; init_u's blend"
                     << " gives 78.0)" << endl;
                found = true; break;
            }
            prev = cur; jprev = j;
        }
        if(!found) cout << "      [CELL ROT] " << (hemi == 0 ? "N" : "S")
                        << " polar band: zonal-mean u does not change sign between 60 and 90 deg"
                        << endl;
    }
}

void cAtmosphereModel::write_meridional_streamfunction(int iter){
    const double a      = r_Earth * 1000.0;   // Earth radius [m] (r_Earth is in km)
    const double two_pi = 2.0 * M_PI;

    // zonal-mean meridional wind [m/s] and zonal-mean meridional MASS flux [kg/(m²·s)],
    // both over fluid cells. r_humid is dimensional [kg/m³], so rho*v*u_0 is already a
    // mass flux.
    // ---- THE DIVISOR, AND WHY IT DECIDED WHETHER CELLS CLOSE ----
    //
    // The zonal mean below used to divide by `n`, the number of FLUID cells at that level. `n`
    // changes with height wherever terrain varies with longitude, so each level was averaged
    // over a DIFFERENT population and the vertical integral of `rvbar` stopped being the zonal
    // mean of the column-integrated mass flux. balance_column_mass_flux() makes every COLUMN's
    // integral vanish -- with ATM_CELLS_FROM_PSI it vanishes to 3.5e-17 -- yet Psi(ground) did
    // not, and the gap was this.
    //
    // Physically the fixed divisor is also the correct one: no mass crosses rock, so a
    // sub-terrain cell contributes ZERO to the zonal mean rather than being excluded from it.
    // Dividing by `n` inflates levels that clip terrain by km/n.
    //
    // MEASURED, and it is the sharpest evidence: with the streamfunction initial condition five
    // of six latitude bands closed and only the SOUTHERN polar one did not (ratio 1.75 against
    // 0.137 for its northern twin). That band is Antarctica -- high terrain varying strongly
    // with longitude, so `n` varies most with height. The Arctic at the same latitudes is ocean,
    // `n` is constant, and it closed.
    //
    // BOTH are written: `psi_kg_per_s` keeps the historical definition so every Psi number
    // recorded in this tree stays comparable, and `psi_fixdiv_kg_per_s` is the corrected one.
    // This tree has already had one silent Psi redefinition (the constant-density fix of
    // 2026-08-26, which made every earlier number a different quantity) and does not need a
    // second one.
    vector<vector<double> > vbar(im, vector<double>(jm, 0.0));
    vector<vector<double> > rvbar(im, vector<double>(jm, 0.0));
    vector<vector<double> > rvbar_fix(im, vector<double>(jm, 0.0));
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            double sum = 0.0, rsum = 0.0; int n = 0;
            for(int k = 0; k < km; k++){
                if(i < i_topography[j][k]) continue;            // inside terrain
                double vv  = v.x[i][j][k];
                double rho = r_humid.x[i][j][k];
                if(!AtomUtils::is_finite_safe(vv)) continue;
                if(!AtomUtils::is_finite_safe(rho) || rho <= 0.0) rho = r_air;  // pre-densities()
                sum += vv; rsum += rho * vv; n++;
            }
            vbar[i][j]      = (n > 0) ? (sum  / n) * u_0 : 0.0;
            rvbar[i][j]     = (n > 0) ? (rsum / n) * u_0 : 0.0;
            rvbar_fix[i][j] = (rsum / (double)km) * u_0;        // rock contributes zero flux
        }
    }

    // meridional mass streamfunction, integrated downward from the lid (Ψ_top = 0)
    vector<vector<double> > psi(im, vector<double>(jm, 0.0));
    vector<vector<double> > psi_fix(im, vector<double>(jm, 0.0));
    for(int j = 0; j < jm; j++){
        const double cosphi = sin(the.z[j]);                    // cos(latitude) = sin(colatitude)
        const double coeff  = two_pi * a * cosphi;
        for(int i = im - 2; i >= 0; i--){
            const double dz   = get_layer_height(i + 1) - get_layer_height(i);   // [m] > 0
            const double rvm  = 0.5 * (rvbar[i][j] + rvbar[i + 1][j]);           // [kg/(m²·s)]
            psi[i][j] = psi[i + 1][j] + coeff * rvm * dz;                        // [kg/s]
            (void)rvbar_fix;                                  // superseded, see below
        }
    }

    // THE EXACT ZONAL MEAN: INTEGRATE EACH COLUMN, THEN AVERAGE -- not the reverse.
    //
    // Averaging the profiles first and applying uniform trapezoid weights afterwards is still
    // not the zonal mean of the column integrals, because each column's integral carries a
    // BOTTOM HALF-WEIGHT 0.5*dz at its OWN ground level i_topography, and that level differs
    // with longitude. Fixing the divisor alone took the southern polar band from 1.750 to 0.496;
    // this takes it the rest of the way, because it is the definition rather than an
    // approximation to it:
    //
    //     Psi(phi,z) = 2*pi*a*cos(phi) * (1/km) * SUM_k INT_z^top rho*v dz
    //
    // With balance_column_mass_flux() making every column's integral vanish, Psi(ground) is then
    // zero to round-off BY CONSTRUCTION -- which is the check, and it is why this is worth the
    // extra loop rather than being left as a 0.5 ratio nobody can interpret.
    for(int j = 0; j < jm; j++){
        const double coeff = two_pi * a * sin(the.z[j]);
        vector<double> acc(im, 0.0);
        for(int k = 0; k < km; k++){
            double running = 0.0;
            for(int i = im - 2; i >= 0; i--){
                if(i < i_topography[j][k]){ acc[i] += running; continue; }   // rock: flux unchanged
                const double dz = get_layer_height(i + 1) - get_layer_height(i);
                double r1 = r_humid.x[i][j][k], r2 = r_humid.x[i+1][j][k];
                if(!AtomUtils::is_finite_safe(r1) || r1 <= 0.0) r1 = r_air;
                if(!AtomUtils::is_finite_safe(r2) || r2 <= 0.0) r2 = r_air;
                const double v1 = v.x[i][j][k], v2 = v.x[i+1][j][k];
                if(AtomUtils::is_finite_safe(v1) && AtomUtils::is_finite_safe(v2))
                    running += 0.5 * (r1 * v1 + r2 * v2) * dz * u_0;
                acc[i] += running;
            }
        }
        for(int i = 0; i < im; i++) psi_fix[i][j] = coeff * acc[i] / (double)km;
    }

    auto lat_of = [&](int j){ return 90.0 - (double)j * 180.0 / (double)(jm - 1); };  // °N positive

    // long-format CSV: lat × height
    ostringstream fname;
    fname << output_path << "meridional_streamfunction_" << iter << ".csv";
    ofstream f(fname.str().c_str());
    if(f.is_open()){
        f << "lat_deg,height_m,vbar_mps,psi_kg_per_s,psi_fixdiv_kg_per_s\n";
        for(int j = 0; j < jm; j++){
            const double lat = lat_of(j);
            for(int i = 0; i < im; i++){
                f << lat << "," << get_layer_height(i) << ","
                  << vbar[i][j] << "," << psi[i][j] << "," << psi_fix[i][j] << "\n";
            }
        }
        f.close();
    }

    // Publish Psi as a field. It is a zonal mean, so it is replicated across k; the point is
    // that the cells become visible in ParaView and appear in the Results min/max, instead of
    // the whole cell-structure argument resting on the scalar Psi_max below.
    for(int i = 0; i < im; i++)
        for(int j = 0; j < jm; j++)
            for(int k = 0; k < km; k++)
                Psi.x[i][j][k] = psi[i][j];

    // log the strongest overturning cells so the structure is readable from the run log
    double psimax = -DBL_MAX, psimin = DBL_MAX;
    int jmx = 0, imx = 0, jmn = 0, imn = 0;
    for(int j = 0; j < jm; j++){
        for(int i = 0; i < im; i++){
            if(psi[i][j] > psimax){ psimax = psi[i][j]; jmx = j; imx = i; }
            if(psi[i][j] < psimin){ psimin = psi[i][j]; jmn = j; imn = i; }
        }
    }
    cout << "      [streamfn] iter=" << iter << fixed << setprecision(2)
         << "  Psi_max=" << psimax / 1.0e9 << " (1e9 kg/s) @ lat=" << lat_of(jmx)
         << " z=" << setprecision(0) << get_layer_height(imx) << "m"
         << setprecision(2) << "   Psi_min=" << psimin / 1.0e9 << " @ lat=" << lat_of(jmn)
         << " z=" << setprecision(0) << get_layer_height(imn) << "m"
         << "   -> " << fname.str() << endl;
}


// Zonal-mean meridional wind vbar[i][j] in m/s, averaged over fluid cells only
// (i >= i_topography, finite) — identical masking/units to write_meridional_streamfunction
// so the two diagnostics are directly comparable.
void cAtmosphereModel::zonal_mean_v(std::vector<std::vector<double> >& vbar){
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            double sum = 0.0; int n = 0;
            for(int k = 0; k < km; k++){
                if(i < i_topography[j][k]) continue;            // inside terrain
                double vv = v.x[i][j][k];
                if(!AtomUtils::is_finite_safe(vv)) continue;
                sum += vv; n++;
            }
            vbar[i][j] = (n > 0) ? (sum / n) * u_0 : 0.0;
        }
    }
}


// Zonal-mean meridional-wind momentum budget: attribute the per-iteration change of
// vbar to each algorithmic step. The four contributions are differenced from snapshots
// of zonal-mean v taken in run_3D_loop and sum (to rounding) to the net Δvbar that iter:
//   dv_dyn    — RK4 net of all rhs_v physics (PGF + Coriolis + advection + diffusion + drag + MC)
//   dv_polar  — polar zonal (φ) filter
//   dv_orog   — orographic Shapiro filter
//   dv_radial — radial (vertical) Shapiro filter  [prime spin-down suspect]
// Units: m/s per iteration. A term that is NEGATIVE where vbar<0 (or positive where
// vbar>0) is ERODING the meridional flow / overturning at that (lat,z).
void cAtmosphereModel::write_v_momentum_budget(int iter,
    const std::vector<std::vector<double> >& dv_dyn,
    const std::vector<std::vector<double> >& dv_polar,
    const std::vector<std::vector<double> >& dv_orog,
    const std::vector<std::vector<double> >& dv_radial){

    std::vector<std::vector<double> > vbar(im, std::vector<double>(jm, 0.0));
    zonal_mean_v(vbar);   // current (post-iteration) vbar

    auto lat_of = [&](int j){ return 90.0 - (double)j * 180.0 / (double)(jm - 1); };

    // RK4 term-split: zonal-mean of each rhs_v contribution captured in vbud_*, scaled
    // from nondim tendency to m/s per iteration (× dt × u_0, the Euler-step contribution)
    // so it is directly comparable to dv_dyn from the differencing budget. Same fluid-cell
    // masking as zonal_mean_v. Their sum (pgf+cor+advv+advh+diff+other) ≈ dv_dyn.
    const double term_scale = dt * u_0;
    auto zmean = [&](const Array& A, std::vector<std::vector<double> >& out){
        for(int i = 0; i < im; i++)
            for(int j = 0; j < jm; j++){
                double sum = 0.0; int n = 0;
                for(int k = 0; k < km; k++){
                    if(i < i_topography[j][k]) continue;
                    double a = A.x[i][j][k];
                    if(!AtomUtils::is_finite_safe(a)) continue;
                    sum += a; n++;
                }
                out[i][j] = (n > 0) ? (sum / n) * term_scale : 0.0;
            }
    };
    std::vector<std::vector<double> > t_pgf(im, std::vector<double>(jm,0.0)), t_cor=t_pgf,
        t_advv=t_pgf, t_advh=t_pgf, t_diff=t_pgf, t_other=t_pgf;
    zmean(vbud_pgf, t_pgf);   zmean(vbud_cor, t_cor);     zmean(vbud_advv, t_advv);
    zmean(vbud_advh, t_advh); zmean(vbud_diff, t_diff);   zmean(vbud_other, t_other);

    // long-format CSV: lat × height × per-step contributions + RK4 term-split (m/s per iter)
    ostringstream fname;
    fname << output_path << "v_momentum_budget_" << iter << ".csv";
    ofstream f(fname.str().c_str());
    if(f.is_open()){
        f << "lat_deg,height_m,vbar_mps,dv_dyn,dv_polar,dv_orog,dv_radial,dv_net,"
          << "pgf,coriolis,adv_vert,adv_horiz,diffusion,drag_mc,dyn_sum\n";
        for(int j = 0; j < jm; j++){
            const double lat = lat_of(j);
            for(int i = 0; i < im; i++){
                const double net = dv_dyn[i][j] + dv_polar[i][j] + dv_orog[i][j] + dv_radial[i][j];
                const double dyn_sum = t_pgf[i][j] + t_cor[i][j] + t_advv[i][j]
                                     + t_advh[i][j] + t_diff[i][j] + t_other[i][j];
                f << lat << "," << get_layer_height(i) << "," << vbar[i][j] << ","
                  << dv_dyn[i][j] << "," << dv_polar[i][j] << "," << dv_orog[i][j] << ","
                  << dv_radial[i][j] << "," << net << ","
                  << t_pgf[i][j] << "," << t_cor[i][j] << "," << t_advv[i][j] << ","
                  << t_advh[i][j] << "," << t_diff[i][j] << "," << t_other[i][j] << ","
                  << dyn_sum << "\n";
            }
        }
        f.close();
    }

    // Log the budget at the mid-latitude lower-troposphere return branch. Track a FIXED
    // probe (lat ~44°, z ~500 m) where the negative-vbar surface layer decays, so the
    // per-term time series is readable. Report the RK4 dynamical split there.
    int jp = (int)round((90.0 - 44.0) * (jm - 1) / 180.0);   // ~lat 44°N
    int ip = get_layer_index(500.0);                         // ~500 m AGL (sea level)
    if(ip < 1) ip = 1;
    cout << "      [vbudget] iter=" << iter << fixed << setprecision(5)
         << "  @lat=" << setprecision(0) << lat_of(jp) << " z=" << get_layer_height(ip) << "m"
         << setprecision(5) << "  vbar=" << vbar[ip][jp]
         << "  | dv_dyn=" << dv_dyn[ip][jp] << " (radial=" << dv_radial[ip][jp] << ")"
         << "  || split: pgf=" << t_pgf[ip][jp] << " cor=" << t_cor[ip][jp]
         << " advV=" << t_advv[ip][jp] << " advH=" << t_advh[ip][jp]
         << " diff=" << t_diff[ip][jp] << " drag/mc=" << t_other[ip][jp]
         << "   -> " << fname.str() << endl;
}


void cAtmosphereModel::zonal_mean_w(std::vector<std::vector<double> >& wbar){
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            double sum = 0.0; int n = 0;
            for(int k = 0; k < km; k++){
                if(i < i_topography[j][k]) continue;            // inside terrain
                double ww = w.x[i][j][k];
                if(!AtomUtils::is_finite_safe(ww)) continue;
                sum += ww; n++;
            }
            wbar[i][j] = (n > 0) ? (sum / n) * u_0 : 0.0;
        }
    }
}


// Zonal-mean ZONAL-wind (w) momentum budget — the trade / Walker component the v budget
// cannot see. Mirror of write_v_momentum_budget: the four step contributions difference to
// the net Δwbar/iter, and the RK4 split (pgf+cor+advv+advh+diff+other) ≈ dw_dyn.
//   dw_dyn    — RK4 net of all rhs_w physics (PGF + Coriolis + advection + diffusion + drag + MC)
//   dw_polar  — polar zonal (φ) filter   (starts at 30° — reaches the poleward trade edge)
//   dw_orog   — orographic Shapiro filter
//   dw_radial — radial (vertical) Shapiro filter  [prime spin-down suspect]
// Units: m/s per iteration. Trades are EASTERLY (wbar<0 in the tropics): a term POSITIVE
// where wbar<0 ERODES the trades; NEGATIVE there SUSTAINS them. wbud_cor is the trade
// source (Coriolis on the meridional Hadley inflow).
void cAtmosphereModel::write_w_momentum_budget(int iter,
    const std::vector<std::vector<double> >& dw_dyn,
    const std::vector<std::vector<double> >& dw_polar,
    const std::vector<std::vector<double> >& dw_orog,
    const std::vector<std::vector<double> >& dw_radial){

    std::vector<std::vector<double> > wbar(im, std::vector<double>(jm, 0.0));
    zonal_mean_w(wbar);   // current (post-iteration) wbar

    auto lat_of = [&](int j){ return 90.0 - (double)j * 180.0 / (double)(jm - 1); };

    const double term_scale = dt * u_0;
    auto zmean = [&](const Array& A, std::vector<std::vector<double> >& out){
        for(int i = 0; i < im; i++)
            for(int j = 0; j < jm; j++){
                double sum = 0.0; int n = 0;
                for(int k = 0; k < km; k++){
                    if(i < i_topography[j][k]) continue;
                    double a = A.x[i][j][k];
                    if(!AtomUtils::is_finite_safe(a)) continue;
                    sum += a; n++;
                }
                out[i][j] = (n > 0) ? (sum / n) * term_scale : 0.0;
            }
    };
    std::vector<std::vector<double> > t_pgf(im, std::vector<double>(jm,0.0)), t_cor=t_pgf,
        t_advv=t_pgf, t_advh=t_pgf, t_diff=t_pgf, t_other=t_pgf;
    zmean(wbud_pgf, t_pgf);   zmean(wbud_cor, t_cor);     zmean(wbud_advv, t_advv);
    zmean(wbud_advh, t_advh); zmean(wbud_diff, t_diff);   zmean(wbud_other, t_other);

    ostringstream fname;
    fname << output_path << "w_momentum_budget_" << iter << ".csv";
    ofstream f(fname.str().c_str());
    if(f.is_open()){
        f << "lat_deg,height_m,wbar_mps,dw_dyn,dw_polar,dw_orog,dw_radial,dw_net,"
          << "pgf,coriolis,adv_vert,adv_horiz,diffusion,drag_mc,dyn_sum\n";
        for(int j = 0; j < jm; j++){
            const double lat = lat_of(j);
            for(int i = 0; i < im; i++){
                const double net = dw_dyn[i][j] + dw_polar[i][j] + dw_orog[i][j] + dw_radial[i][j];
                const double dyn_sum = t_pgf[i][j] + t_cor[i][j] + t_advv[i][j]
                                     + t_advh[i][j] + t_diff[i][j] + t_other[i][j];
                f << lat << "," << get_layer_height(i) << "," << wbar[i][j] << ","
                  << dw_dyn[i][j] << "," << dw_polar[i][j] << "," << dw_orog[i][j] << ","
                  << dw_radial[i][j] << "," << net << ","
                  << t_pgf[i][j] << "," << t_cor[i][j] << "," << t_advv[i][j] << ","
                  << t_advh[i][j] << "," << t_diff[i][j] << "," << t_other[i][j] << ","
                  << dyn_sum << "\n";
            }
        }
        f.close();
    }

    // Log the budget in the tropical lower-troposphere trade layer (lat ~18°, z ~500 m),
    // where the easterly trades live and decay. Report the RK4 dynamical split there.
    int jp = (int)round((90.0 - 18.0) * (jm - 1) / 180.0);   // ~lat 18°N
    int ip = get_layer_index(500.0);                         // ~500 m AGL (sea level)
    if(ip < 1) ip = 1;
    cout << "      [wbudget] iter=" << iter << fixed << setprecision(5)
         << "  @lat=" << setprecision(0) << lat_of(jp) << " z=" << get_layer_height(ip) << "m"
         << setprecision(5) << "  wbar=" << wbar[ip][jp]
         << "  | dw_dyn=" << dw_dyn[ip][jp] << " (radial=" << dw_radial[ip][jp]
         << " polar=" << dw_polar[ip][jp] << ")"
         << "  || split: pgf=" << t_pgf[ip][jp] << " cor=" << t_cor[ip][jp]
         << " advV=" << t_advv[ip][jp] << " advH=" << t_advh[ip][jp]
         << " diff=" << t_diff[ip][jp] << " drag/mc=" << t_other[ip][jp]
         << "   -> " << fname.str() << endl;
}

