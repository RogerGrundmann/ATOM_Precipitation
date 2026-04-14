/*
 * Atmosphere General Circulation Modell(AGCM) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
*/

#include <iostream>
#include <iomanip>
#include <cstring>
#include <cfloat>
#include <cmath>
#include "cAtmosphereModel.h"

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

