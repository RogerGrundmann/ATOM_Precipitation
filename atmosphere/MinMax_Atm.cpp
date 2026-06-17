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
*   Ψ(i,j) = 2π·a·cosφ·ρ0·∫_z^{z_top} [v]_phys dz'      [kg/s]
*
* integrated DOWNWARD from the model lid (Ψ = 0 at i = im-1). With v > 0 southward, the
* cell STRENGTH is max|Ψ|; a tropical Ψ extremum is the Hadley cell. Writes a long-format
* CSV (lat × height) per checkpoint and logs the strongest cells so collapse can be read
* straight from the run log without ParaView. See [[project_vw_drag_cut_hadley]].
*/
void cAtmosphereModel::write_meridional_streamfunction(int iter){
    const double a      = r_Earth * 1000.0;   // Earth radius [m] (r_Earth is in km)
    const double rho    = r_air;              // Boussinesq reference density [kg/m³]
    const double two_pi = 2.0 * M_PI;

    // zonal-mean meridional wind [m/s] over fluid cells
    vector<vector<double> > vbar(im, vector<double>(jm, 0.0));
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

    // meridional mass streamfunction, integrated downward from the lid (Ψ_top = 0)
    vector<vector<double> > psi(im, vector<double>(jm, 0.0));
    for(int j = 0; j < jm; j++){
        const double cosphi = sin(the.z[j]);                    // cos(latitude) = sin(colatitude)
        const double coeff  = two_pi * a * cosphi * rho;
        for(int i = im - 2; i >= 0; i--){
            const double dz   = get_layer_height(i + 1) - get_layer_height(i);   // [m] > 0
            const double vmid = 0.5 * (vbar[i][j] + vbar[i + 1][j]);             // [m/s]
            psi[i][j] = psi[i + 1][j] + coeff * vmid * dz;                       // [kg/s]
        }
    }

    auto lat_of = [&](int j){ return 90.0 - (double)j * 180.0 / (double)(jm - 1); };  // °N positive

    // long-format CSV: lat × height
    ostringstream fname;
    fname << output_path << "meridional_streamfunction_" << iter << ".csv";
    ofstream f(fname.str().c_str());
    if(f.is_open()){
        f << "lat_deg,height_m,vbar_mps,psi_kg_per_s\n";
        for(int j = 0; j < jm; j++){
            const double lat = lat_of(j);
            for(int i = 0; i < im; i++){
                f << lat << "," << get_layer_height(i) << ","
                  << vbar[i][j] << "," << psi[i][j] << "\n";
            }
        }
        f.close();
    }

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

