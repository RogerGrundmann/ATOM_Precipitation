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

