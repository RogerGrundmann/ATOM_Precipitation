/*
 * Atmosphere General Circulation Modell(AGCM)applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
*/

#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>

#include "Utils.h"
#include "cCubeModel.h"

using namespace std;
using namespace AtomUtils;

void cCubeModel::print_min_max_atm(){
/*
    cout << endl << endl << endl << "      print_min_max_atm" << endl;

    auto begin = std::chrono::high_resolution_clock::now();
*/
    cout << endl << endl << endl << " flow properties: " << endl << endl;


    cout << endl << " velocity components: " << endl;
    searchMinMax_3D(" max u-component ", " min u-component ", "m/s", u, u_0);
    searchMinMax_3D(" max v-component ", " min v-component ", "m/s", v, u_0);
    searchMinMax_3D(" max w-component ", " min w-component ", "m/s", w, u_0);



    cout << endl << " pressures: " << endl;
    searchMinMax_3D(" max pressure static ", " min pressure static ", 
        "hPa", p_stat, 1.0);
    searchMinMax_3D(" max pressure dynamic ", " min pressure dynamic ", 
        "hPa", p_dyn, p_0);



    if(turb_model != "none"){
        cout << endl << " turbulence: " << endl;
        searchMinMax_3D(" max TKE ",          " min TKE ",          "m²/s²",     tke,        1.0);
        searchMinMax_3D(" max dissipation ",  " min dissipation ",  "m²/s³",     dis,        1.0);
        searchMinMax_3D(" max eddy viscosity "," min eddy viscosity ","m²/s",    nue,        1.0);
        searchMinMax_3D(" max turb production "," min turb production ","m²/s³", prod,       1.0);
        searchMinMax_3D(" max TKE source ",   " min TKE source ",   "m²/s³",     tke_source, 1.0);
        searchMinMax_3D(" max dis source ",   " min dis source ",   "m²/s⁴",     dis_source, 1.0);
    }

    cout << endl << endl;
/*
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for print_min_max_atm\n", elapsed.count() * 1e-9);

    cout << "      print_min_max_atm ended" << endl;
*/
}
