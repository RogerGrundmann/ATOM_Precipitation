/*
 * Ocean General Circulation Modell(OGCM) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
*/

#include "cHydrosphereModel.h"

using namespace std;
using namespace AtomUtils;

// NaN/Inf test must survive -ffast-math. The project compiles with -ffast-math
// (=> -ffinite-math-only), under which the compiler assumes no non-finite values
// exist and constant-folds std::isfinite()/std::isnan() to a constant, so a scanner
// built on them silently never fires. Use the project's canonical bit-level test
// AtomUtils::is_finite_safe (lib/Utils.h) instead — it inspects the IEEE-754
// exponent field with no FP comparison the optimiser can defeat.
static inline bool nonfinite_bits(double v){ return !AtomUtils::is_finite_safe(v); }

void cHydrosphereModel::print_min_max_hyd(){
    cout << endl << " flow properties: " << endl << endl;
    searchMinMax_3D(" max temperature ", " min temperature ", 
        " deg", t, 273.15, [](double i)->double{return i - 273.15;}, true);
    searchMinMax_3D(" max u-component ", " min u-component ", "m/s", u, u_0);
    searchMinMax_3D(" max v-component ", " min v-component ", "m/s", v, u_0);
    searchMinMax_3D(" max w-component ", " min w-component ", "m/s", w, u_0);
    searchMinMax_3D(" max pressure static ", " min pressure static ", 
        "bar", p_hydro, 1.0);
    searchMinMax_3D(" max pressure dynamic ", " min pressure dynamic ", 
        "hPa", p_dyn, p_0);
    searchMinMax_3D(" max water density ", " min water density ", "kg/m3", r_water, 1.0);
    searchMinMax_3D(" max salt water density ", " min salt water density ", "kg/m3", r_salt_water, 1.0);

    cout << endl << " salinity based results in the three dimensional space: " << endl << endl;
    searchMinMax_3D(" max saltinity ", " min saltinity ", "psu", c, c_35);
    searchMinMax_3D(" max salt balance ", " min salt balance ", "kg/m3", Salt_Balance, 1.0);
    searchMinMax_3D(" max salt finger ", " min salt finger ", "kg/m3", Salt_Finger, 1.0);
    searchMinMax_3D(" max salt diffusion ", " min salt diffusion ", "kg/m3", Salt_Diffusion, 1.0);

    cout << endl << " forces per unit volume: " << endl << endl;
    searchMinMax_3D(" max pressure force ", " min pressure force ", "kN/m3", PresGradForce, 1.0);
    searchMinMax_3D(" max buoyancy force ", " min buoyancy force ", "kN/m3", BuoyancyForce, 1.0);
    searchMinMax_3D(" max Coriolis force ", " min Coriolis force ", "N/m3", CoriolisForce, 1.0);
    searchMinMax_3D(" max centrifugal force ", " min centrifugal force ", "N/m3", CentrifugalForce, 1.0);

    cout << endl << " salinity increase due to evaporation: " << endl << endl;
    searchMinMax_2D(" max precipitation ", " min precipitation ", "mm/d", Precipitation_2D, 1.0);
    searchMinMax_2D(" max evaporation ", " min evaporation ", "mm/d", Evaporation, 1.0);
    searchMinMax_2D(" max sal_evap ", " min sal_evap ", "psu", salinity_evaporation, c_35);

    cout << endl << " mass flows like Ekman pumping, up- and downwelling near coasts: " << endl << endl;
    searchMinMax_2D(" max EkmanPumping ", " min EkmanPumping ", "kg/m²s", EkmanPumping, 1.0);
    searchMinMax_2D(" max upwelling ", " min upwelling ", "kg/m²s", Upwelling, 1.0);
    searchMinMax_2D(" max downwelling ", " min downwelling ", "kg/m²s", Downwelling, 1.0);

    cout << endl << " solid soil contours like sea mounts: " << endl << endl;
    searchMinMax_2D(" max bathymetry ", " min bathymetry ", "m", Bathymetry, 1.0);

    if (turb_model != "none") {
        cout << endl << " turbulence: " << endl;
        searchMinMax_3D(" max TKE ",           " min TKE ",           "m²/s²", tke,        1.0);
        searchMinMax_3D(" max dissipation ",   " min dissipation ",   "m²/s³", dis,        1.0);
        searchMinMax_3D(" max eddy viscosity "," min eddy viscosity ","m²/s",  nue,        1.0);
        searchMinMax_3D(" max turb production "," min turb production ","m²/s³", prod,     1.0);
        searchMinMax_3D(" max TKE source ",    " min TKE source ",    "m²/s³", tke_source, 1.0);
        searchMinMax_3D(" max dis source ",    " min dis source ",    "m²/s⁴", dis_source, 1.0);
    }
}
/*
 *
*/



// First-NaN scanner: sweep the prognostic fields (u, v, w, t, c, p_dyn) for the
// first non-finite cell and report field + (i,j,k) + lat/lon + iteration, then
// abort. Mirrors the atmosphere's scanForNaN. The hydro 0Ma spin-up dies ~iter700
// in a whack-a-mole where each field-fix relocates the NaN; MinMax only exposes
// u/v/w indirectly (via EkmanPumping at 90N), so we cannot see *which* prognostic
// field and *which* cell fails first. Called right after the RK4 solve each iter
// (before BC/filters can move or mask the blow-up). Grid: i=0 deepest .. im-1
// surface; lat = 90 - j (deg); lon = k (deg E). Returns true if a NaN was found.
bool cHydrosphereModel::scanForNaN_hyd(int iter, const char *stage){
    bool found = false;

    // 3D fields: prognostic (u,v,w,t,c,p_dyn) + derived ones in the EkmanPumping
    // chain (aux_v, aux_w, r_salt_water, r_water) so an end-of-iter scan after the
    // heavy block can attribute a diagnostic NaN to its real source.
    struct Field3D { const char *name; Array *a; };
    Field3D fields3d[] = {
        {"u",            &u},
        {"v",            &v},
        {"w",            &w},
        {"t",            &t},
        {"c",            &c},
        {"p_dyn",        &p_dyn},
        {"aux_v",        &aux_v},
        {"aux_w",        &aux_w},
        {"r_salt_water", &r_salt_water},
        {"r_water",      &r_water},
    };
    for(auto &f : fields3d){
        for(int i = 0; i < im && !found; i++)
            for(int j = 0; j < jm && !found; j++)
                for(int k = 0; k < km && !found; k++){
                    double val = f.a->x[i][j][k];
                    if(nonfinite_bits(val)){
                        cout << endl
                             << "      OGCM: *** NaN/Inf DETECTED *** stage=" << stage
                             << "  iter=" << iter
                             << "  field=" << f.name
                             << "  value=" << val << endl
                             << "      OGCM:   cell (i,j,k) = (" << i << "," << j << "," << k << ")"
                             << "  lat=" << (90 - j) << " deg"
                             << "  lon=" << k << " deg E"
                             << "  level i=" << i << " of " << (im-1) << " (im-1=surface)"
                             << "  water=" << (is_water(h, i, j, k) ? "yes" : "land") << endl;
                        found = true;
                    }
                }
    }

    // 2D surface diagnostics (EkmanPumping etc.) — these have a polar 1/sinthe
    // factor and are output-only, but report them so we can confirm whether the
    // first NaN is genuinely born here or inherited from a 3D field above.
    struct Field2D { const char *name; Array_2D *a; };
    Field2D fields2d[] = {
        {"EkmanPumping", &EkmanPumping},
        {"Upwelling",    &Upwelling},
    };
    for(auto &f : fields2d){
        for(int j = 0; j < jm && !found; j++)
            for(int k = 0; k < km && !found; k++){
                double val = f.a->y[j][k];
                if(nonfinite_bits(val)){
                    cout << endl
                         << "      OGCM: *** NaN/Inf DETECTED (2D) *** stage=" << stage
                         << "  iter=" << iter
                         << "  field=" << f.name
                         << "  value=" << val << endl
                         << "      OGCM:   cell (j,k) = (" << j << "," << k << ")"
                         << "  lat=" << (90 - j) << " deg"
                         << "  lon=" << k << " deg E" << endl;
                    found = true;
                }
            }
    }

    return found;
}
/*
 *
*/

