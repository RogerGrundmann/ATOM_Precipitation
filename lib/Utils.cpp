#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstring>

#include <Utils.h>
#include "cAtmosphereModel.h"
#include "cHydrosphereModel.h"

using namespace AtomUtils;
using namespace std;

std::vector<std::vector<double> > m_node_weights;

std::ofstream& AtomUtils::get_logger(){
    static std::ofstream o("atom_log.txt", std::ofstream::out);
    return o;
}
/*
*
*/
HemisphereCoords AtomUtils::convert_coords(double lon, double lat){
    HemisphereCoords ret;
    if ( lat > 90 ){
        ret.lat = lat - 90;
        ret.north_or_south = "°S";
    }else{
        ret.lat = 90 - lat;
        ret.north_or_south = "°N";
    }
    if ( lon > 180 ){
        ret.lon = 360 - lon;
        ret.east_or_west = "°W";
    }else{
        ret.lon = lon;
        ret.east_or_west = "°E";
    }
    return ret;
}
/*
* [-180, 180] to [0, 360]
*/
int AtomUtils::lon_2_index(double lon){
    if(lon<0)
        return int(round(360+lon));
    else
        return int(round(lon));
}
/*
* [-90, 90] to [0, 180]
*/
int AtomUtils::lat_2_index(double lat){
    return int(round(90-lat));
}
/*
*
*/
//change data coordinate system from -180° _ 0° _ +180° to 0°- 360°
void AtomUtils::move_data(double* data, int len){
    for(int i=0; i<len/2; i++){
        std::iter_swap(data+i, data+len/2+i);
    }
    data[len-1] = data[0];
    return;
}
/*
*
*/
void AtomUtils::move_data(std::vector<int>& data, int len){   
    for(int i=0; i<len/2; i++){
        std::iter_swap(data.begin()+i, data.begin()+len/2+i);
    }
    data[len-1] = data[0];
    return;
}
/*
*
*/
//the size of arrays must be the same with the size of new_arrays
//the values in the new_arrays will be assigned to coeff * old_values
void AtomUtils::move_data_to_new_arrays(int im, int jm, int km, double coeff, 
    std::vector<Array*>& arrays, std::vector<Array*>& new_arrays){
    assert(arrays.size() == new_arrays.size());
    for ( int i = 0; i < im; i++ ){
        move_data_to_new_arrays(jm, km, coeff, arrays, new_arrays, i);
    }
    return;
}
/*
*
*/
void AtomUtils::move_data_to_new_arrays(int jm, int km, double coeff, 
    std::vector<Array*>& arrays, std::vector<Array*>& new_arrays, int i){
    assert(arrays.size() == new_arrays.size());
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++ ){
            for(std::size_t n=0; n<arrays.size(); n++){
                new_arrays[n]->x[ i ][ j ][ k ] = coeff * arrays[n]->x[ i ][ j ][ k ];
            }
        }
    }
    return;
}
/*
*
*/
std::tuple<double, int, int, int>
AtomUtils::max_diff(int im, int jm, int km, const Array &a1, const Array &a2){
    std::tuple<double, int, int, int> ret;
    for ( int i = 0; i < im; i++ ){
        for ( int j = 0; j < jm; j++ ){
            for ( int k = 0; k < km; k++ ){
                double tmp = fabs ( a1.x[ i ][ j ][ k ] - a2.x[ i ][ j ][ k ] );
                if(tmp > std::get<0>(ret)){
                    std::get<0>(ret) = tmp;
                    std::get<1>(ret) = i;
                    std::get<2>(ret) = j;
                    std::get<3>(ret) = k;
                }
            }
        }
    }
    return ret;
}
/*
*
*/
double AtomUtils::C_Dalton(int i, int j, int k, double coeff_Dalton, 
    double u_0, Array &u, Array &v, Array &w){
    // variation of the evaporation coefficient in Dalton's evaporation law, parabola
    // air velocity measured 2m above sea level
    // for v_max = 10 m/s, but C_Dalton is function of v, should be included
    // Geiger ( 1961 ) by > Zmarsly, Kuttler, Pethe in mm/( h * hPa ), p. 133
    // for the modern world the global precipitation is 10% higher than evaporation

    double v_max = 10.0;                                                // in [m/s] ..... Geiger ( 1961 ) by Zmarsly, Kuttler, Pethe in m/s, p. 133
    double c_Dalton_max = 0.053;                                        // in [mm/(h*hPa)]
    double vel_magnitude = sqrt((u.x[i][j][k] * u.x[i][j][k] 
                               + v.x[i][j][k] * v.x[i][j][k] 
                               + w.x[i][j][k] * w.x[i][j][k]) / 3.0) * u_0;// in [m/s]
    return coeff_Dalton 
        * sqrt(c_Dalton_max * c_Dalton_max * vel_magnitude/v_max);      // result in [mm/(h*hPa)]
}
/*
*
*/
void AtomUtils::read_IC(const string& fn, double** a, int jm, int km){

    cout << endl << "      read_IC" << endl;

    cout << "      " << fn << endl;

    ifstream ifs(fn);

    if(!ifs.is_open()){
        cerr << "ERROR: unable to open " << fn << "\n";
        abort();
    }

    double lat, lon, d;

    for(int k=0; k < km && !ifs.eof(); k++){
        for(int j=0; j < jm; j++){
            ifs >> lat >> lon >> d;
            a[j][k] = d;
        }
    }

    cout << "      read_IC ended" << endl;
    return;
}
/*
*
*/
void AtomUtils::smooth_steps(int k, int im, int jm, Array &value){
    double coeff = 0.5;
    std::vector<std::vector<double> > inter(im, std::vector<double>(jm, 0));
    for(int i = 0; i < im ; i++){
        for(int j = 0; j < jm; j++){
            inter[i][j] = value.x[i][j][k];
        }
    }
    for(int i = 1; i < im-1 ; i++){
        for(int j = 1; j < jm-1; j++){
            value.x[i][j][k] = (coeff * inter[i][j] 
                + (1.0 - coeff) * (inter[i+1][j] 
                + inter[i-1][j] + inter[i][j+1] 
                + inter[i][j-1])/4.0);
        } 
    }
    return;
}
/*
*
*/
void AtomUtils::smooth_tropopause(int jm, std::vector<double> &value){
    double coeff = 0.5;
    std::vector<double>inter(jm, 0);
    for(int j = 0; j < jm; j++){
        inter[j] = value[j];
    }
    for(int j = 2; j < jm-2; j++){
        value[j] = (coeff * inter[j] + (1.0 - coeff) 
            * (inter[j+1] + inter[j-1] + inter[j+2] + inter[j-2])/4.0);
    } 
    return;
}
/*
*
*/
double AtomUtils::simpson(int n1, int n2, double dstep, Array_1D &value){
        double sum_even=0, sum_odd=0;
        if(n2 % 2 == 0){
            for(int i = n1+1; i < n2; i+=2){sum_odd += 4 * value.z[i];}
            for(int i = n1+2; i < n2; i+=2){sum_even += 2 * value.z[i];}
        }else cout << "       n2    must be an even number to use the \
            Simpson integration method" << endl;
    return dstep/3 * (value.z[n1] + sum_odd + sum_even + value.z[n2]); // Simpson Rule integration
}
/*
*
*/
double AtomUtils::trapezoidal(int n1, int n2, double dstep, Array_1D &value){
        double sum = 0;
        for(int i = n1+1; i < n2; i++){sum += 2 * value.z[i];}
        return dstep/2 * (value.z[n1] + sum + value.z[n2]);  // Trapezoidal Rule integration
}
/*
*
*/
double AtomUtils::rectangular(int n1, int n2, double dstep, Array_1D &value){
        double sum = 0;
        for(int i = n1; i <= n2; i++){sum += value.z[i];}
        return dstep * sum;  // Rectangular Rule integration
}
/*
*
*/
void AtomUtils::load_map_from_file(const std::string& fn, 
                                   std::map<float, float>& m) {
    std::ifstream f(fn);
    if (!f.is_open()) {
        std::cerr << "error occurred while opening file: " << fn << std::endl;
        return;
    }

    std::string line;
    while (std::getline(f, line)) {
//        if (line.empty() || line[0] == '#') continue;
        size_t first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos || line[first] == '#') continue;
        float key = 0, val = 0;
        std::stringstream ss(line);
        
        if (ss >> key >> val) {
            m[key] = val;
        }
    }
}
/*
*
*/
double AtomUtils::GetMean_3D(int jm, int km, Array &val_3D){
    if(m_node_weights.size() != (unsigned)jm){
        CalculateNodeWeights(jm, km);
    }
    double ret=0.0, weight=0.0;
    for(int j=0; j<jm; j++){
        for(int k=0; k<km; k++){
            //std::cout << (val_3D.x[0][j][k]-1)*t_0 << "  " << m_node_weights[j][k] << std::endl;
            ret += val_3D.x[0][j][k] * m_node_weights[j][k];
            weight += m_node_weights[j][k];
        }
    }
    return ret/weight;
}
/*
*
*/
double AtomUtils::GetMean_2D(int jm, int km, Array_2D &val_2D){
    if(m_node_weights.size() != (unsigned)jm){
        CalculateNodeWeights(jm, km);
    }
    double ret=0.0, weight=0.0;
    for(int j=0; j<jm; j++){
        for(int k=0; k<km; k++){
            //std::cout << (val_2D.y[j][k]-1)*t_0 << "  " << m_node_weights[j][k] << std::endl;
            ret += val_2D.y[j][k] * m_node_weights[j][k];
            weight += m_node_weights[j][k];
        }
    }
    return ret/weight;
}
/*
*
*/
void AtomUtils::CalculateNodeWeights(int jm, int km){
    //use cosine of latitude as weights for now
    //longitudes: 0-360(km) latitudes: 90-(-90)(jm)
    double weight = 0.0;
    m_node_weights.clear();
    for(int i=0; i<jm; i++){
        if(i<=90){
            weight = cos((90-i) * M_PI/180.0);
        }else{
            weight = cos((i-90) * M_PI/180.0);
        }
        m_node_weights.push_back(std::vector<double>());
        m_node_weights[i].resize(km, weight);
    }
    return;
}
/*
*
*/
int AtomUtils::RunStart(string comment){
    string at = "AGCM";
    string hy = "OGCM";
    string comment_1 = "";
    string comment_2 = "";
    if(comment.compare(at) == 0){
        comment_1 = " ... AGCM: time and date at run time begin:   ";
        comment_2 = "    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   \
             Atmosphere General Circulation Model \
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n";
    }
    if(comment.compare(hy) == 0){
        comment_1 = " ... OGCM: time and date at run time begin:   ";
        comment_2 = "    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   \
            Ocean General Circulation Model \
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n";
    }
    std::time_t Run_start;
    struct tm * timeinfo_begin;
    std::time(&Run_start);
    timeinfo_begin = std::localtime(&Run_start);
    std::cout << std::endl << std::endl;
    std::cout << comment_2 
        << std::endl << std::endl;
    std::cout << comment_1 
        << std::asctime(timeinfo_begin);
    return Run_start;
}
/*
*
*/
int AtomUtils::RunEnd(string comment, int Ma, int Run_start){
    string at = "AGCM";
    string hy = "OGCM";
    string comment_1 = "";
    string comment_2 = "";
        if(comment.compare(at) == 0){
        comment_1 = " ... AGCM: time and date at run time begin:   ";
        comment_2 = "    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   \
             Atmosphere General Circulation Model \
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n";
    }
    if(comment.compare(hy) == 0){
        comment_1 = " ... OGCM: time and date at run time begin:   ";
        comment_2 = "    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   \
            Ocean General Circulation Model \
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n";
    }
    std::time_t Run_end;
    std::time(&Run_end);
    struct tm * timeinfo_end;
    timeinfo_end = std::localtime(&Run_end);
    std::cout << std::endl << std::endl;
    std::cout << comment_2 
        << std::endl << std::endl;
    std::cout << comment_1 
        << std::asctime(timeinfo_end);
    int Run_total = Run_end - Run_start;
    int Run_total_minutes = Run_total/60;
    int Run_total_hours = Run_total_minutes/60;
    std::cout << std::endl << " ... computed time slice:"
        << "  Ma = " << Ma << std::endl << std::endl
        << " ... computer time needed:" << std::endl << std::endl
        << setw(20) << setfill(' ') << Run_total 
        << " seconds" << std::endl
        << " ... compares to:" << std::endl << std::endl
        << setw(20) << setfill(' ') << Run_total_hours 
        << " hours" << std::endl
        << setw(20) << setfill(' ') << Run_total_minutes 
        << " minutes" << std::endl
        << setw(20) << setfill(' ') << Run_total%60 
        << " seconds" << std::endl << std::endl
        << std::endl;
    return 0;
}
/*
*
*/
void AtomUtils::use_presets(bool preset_1, bool preset_2, bool preset_3){
    if(preset_1)
        cout << endl << "      use_earthbyte_reconstruction is in use" << endl;
    else
        cout << endl << "      use_earthbyte_reconstruction is NOT in use" << endl;

    if(preset_2)
        cout << "      use_NASA_temperature is in use" << endl;
    else
        cout << "      use_NASA_temperature is NOT in use" << endl;

    if(preset_3)
        cout << "      use_NASA_velocity is in use" << endl;
    else
        cout << "      use_NASA_velocity is NOT in use" << endl;
    return;
}

// ============================================================================
// Shapiro (1-2-1) wiggle-damping filter
// ============================================================================
void AtomUtils::damp_wiggles(Array& field,
                              const std::vector<std::vector<int>>* i_surface,
                              bool along_i, bool along_j, bool along_k,
                              double strength,
                              int passes)
{
    const int im = field.im;
    const int jm = field.jm;
    const int km = field.km;

    // Coefficient for the Shapiro update: f += coeff*(f_{-1} - 2f + f_{+1})
    const double coeff = strength * 0.25;

    // Helper: true when cell (i,j,k) is in the fluid domain.
    // Works for both i_topography (atmosphere) and i_bathymetry (hydrosphere):
    // in both cases i_surface[j][k] is the first fluid cell; solid cells have i < i_surface[j][k].
    auto in_fluid = [&](int i, int j, int k) -> bool {
        if (!i_surface) return true;
        // std::max guards against i_surface[j][k] < 0 (bathymetry deeper than L_hyd)
        return i >= std::max((*i_surface)[j][k], 0);
    };

    // Scratch buffer — snapshot before each axis pass to avoid in-place corruption.
    std::vector<double> tmp(im * jm * km);

    auto idx = [&](int i, int j, int k) { return i * jm * km + j * km + k; };

    for (int pass = 0; pass < passes; ++pass) {

        // --- along k (longitude, periodic) ---
        if (along_k) {
            for (int i = 0; i < im; ++i)
                for (int j = 0; j < jm; ++j)
                    for (int k = 0; k < km; ++k)
                        tmp[idx(i,j,k)] = field.x[i][j][k];

            #pragma omp parallel for collapse(2) schedule(static)
            for (int i = 0; i < im; ++i) {
                for (int j = 0; j < jm; ++j) {
                    for (int k = 0; k < km; ++k) {
                        if (!in_fluid(i, j, k)) continue;
                        const int km1 = (k - 1 + km) % km;
                        const int kp1 = (k + 1)      % km;
                        // no-flux at solid neighbours: substitute current cell value
                        const double v_km1 = in_fluid(i, j, km1) ? tmp[idx(i,j,km1)] : tmp[idx(i,j,k)];
                        const double v_kp1 = in_fluid(i, j, kp1) ? tmp[idx(i,j,kp1)] : tmp[idx(i,j,k)];
                        field.x[i][j][k] = tmp[idx(i,j,k)]
                            + coeff * (v_km1 - 2.0 * tmp[idx(i,j,k)] + v_kp1);
                    }
                }
            }
        }

        // --- along j (latitude, clamped at poles) ---
        if (along_j) {
            for (int i = 0; i < im; ++i)
                for (int j = 0; j < jm; ++j)
                    for (int k = 0; k < km; ++k)
                        tmp[idx(i,j,k)] = field.x[i][j][k];

            #pragma omp parallel for collapse(2) schedule(static)
            for (int i = 0; i < im; ++i) {
                for (int j = 0; j < jm; ++j) {
                    const int jm1 = std::max(j - 1, 0);
                    const int jp1 = std::min(j + 1, jm - 1);
                    for (int k = 0; k < km; ++k) {
                        if (!in_fluid(i, j, k)) continue;
                        // no-flux at solid neighbours: substitute current cell value
                        const double v_jm1 = in_fluid(i, jm1, k) ? tmp[idx(i,jm1,k)] : tmp[idx(i,j,k)];
                        const double v_jp1 = in_fluid(i, jp1, k) ? tmp[idx(i,jp1,k)] : tmp[idx(i,j,k)];
                        field.x[i][j][k] = tmp[idx(i,j,k)]
                            + coeff * (v_jm1 - 2.0 * tmp[idx(i,j,k)] + v_jp1);
                    }
                }
            }
        }

        // --- along i (vertical, clamped at fluid floor and top) ---
        if (along_i) {
            for (int i = 0; i < im; ++i)
                for (int j = 0; j < jm; ++j)
                    for (int k = 0; k < km; ++k)
                        tmp[idx(i,j,k)] = field.x[i][j][k];

            #pragma omp parallel for collapse(2) schedule(static)
            for (int j = 0; j < jm; ++j) {
                for (int k = 0; k < km; ++k) {
                    // clamp to [0, im-1]: guards against i_bathymetry < 0 when depth > L_hyd
                    const int i_surf = i_surface
                        ? std::max((*i_surface)[j][k], 0) : 0;
                    for (int i = i_surf; i < im; ++i) {
                        const int im1 = std::max(i - 1, i_surf);
                        const int ip1 = std::min(i + 1, im - 1);
                        field.x[i][j][k] = tmp[idx(i,j,k)]
                            + coeff * (tmp[idx(im1,j,k)]
                                     - 2.0 * tmp[idx(i,j,k)]
                                     + tmp[idx(ip1,j,k)]);
                    }
                }
            }
        }

    } // passes
}



// ============================================================================
// Extreme-peak removal filter
// ============================================================================
// Along each enabled axis the two direct neighbours define a local background
// (mean) and variability scale (half-spread = |f[n+1]-f[n-1]|/2).
// Cells that deviate from the background by more than threshold*half_spread
// are replaced by the background.  Smooth gradients are untouched.
void AtomUtils::remove_peaks(Array& field,
                              const std::vector<std::vector<int>>* i_surface,
                              bool along_i, bool along_j, bool along_k,
                              double threshold,
                              int passes)
{
    const int im = field.im;
    const int jm = field.jm;
    const int km = field.km;

    auto in_fluid = [&](int i, int j, int k) -> bool {
        if (!i_surface) return true;
        return i >= std::max((*i_surface)[j][k], 0);
    };

    std::vector<double> tmp(im * jm * km);
    auto idx = [&](int i, int j, int k) { return i * jm * km + j * km + k; };

    for (int pass = 0; pass < passes; ++pass) {

        // --- along k (longitude, periodic) ---
        if (along_k) {
            for (int i = 0; i < im; ++i)
                for (int j = 0; j < jm; ++j)
                    for (int k = 0; k < km; ++k)
                        tmp[idx(i,j,k)] = field.x[i][j][k];

            #pragma omp parallel for collapse(2) schedule(static)
            for (int i = 0; i < im; ++i) {
                for (int j = 0; j < jm; ++j) {
                    for (int k = 0; k < km; ++k) {
                        if (!in_fluid(i, j, k)) continue;
                        const int km1 = (k - 1 + km) % km;
                        const int kp1 = (k + 1)      % km;
                        const double v   = tmp[idx(i,j,k)];
                        const double v_m = in_fluid(i, j, km1) ? tmp[idx(i,j,km1)] : v;
                        const double v_p = in_fluid(i, j, kp1) ? tmp[idx(i,j,kp1)] : v;
                        const double mean        = 0.5 * (v_m + v_p);
                        const double half_spread = 0.5 * std::abs(v_p - v_m);
                        if (std::abs(v - mean) > threshold * half_spread)
                            field.x[i][j][k] = mean;
                    }
                }
            }
        }

        // --- along j (latitude, clamped at poles) ---
        if (along_j) {
            for (int i = 0; i < im; ++i)
                for (int j = 0; j < jm; ++j)
                    for (int k = 0; k < km; ++k)
                        tmp[idx(i,j,k)] = field.x[i][j][k];

            #pragma omp parallel for collapse(2) schedule(static)
            for (int i = 0; i < im; ++i) {
                for (int j = 0; j < jm; ++j) {
                    const int jm1 = std::max(j - 1, 0);
                    const int jp1 = std::min(j + 1, jm - 1);
                    for (int k = 0; k < km; ++k) {
                        if (!in_fluid(i, j, k)) continue;
                        const double v   = tmp[idx(i,j,k)];
                        const double v_m = in_fluid(i, jm1, k) ? tmp[idx(i,jm1,k)] : v;
                        const double v_p = in_fluid(i, jp1, k) ? tmp[idx(i,jp1,k)] : v;
                        const double mean        = 0.5 * (v_m + v_p);
                        const double half_spread = 0.5 * std::abs(v_p - v_m);
                        if (std::abs(v - mean) > threshold * half_spread)
                            field.x[i][j][k] = mean;
                    }
                }
            }
        }

        // --- along i (vertical, clamped at fluid floor and top) ---
        if (along_i) {
            for (int i = 0; i < im; ++i)
                for (int j = 0; j < jm; ++j)
                    for (int k = 0; k < km; ++k)
                        tmp[idx(i,j,k)] = field.x[i][j][k];

            #pragma omp parallel for collapse(2) schedule(static)
            for (int j = 0; j < jm; ++j) {
                for (int k = 0; k < km; ++k) {
                    const int i_surf = i_surface
                        ? std::max((*i_surface)[j][k], 0) : 0;
                    for (int i = i_surf; i < im; ++i) {
                        const int im1 = std::max(i - 1, i_surf);
                        const int ip1 = std::min(i + 1, im - 1);
                        const double v   = tmp[idx(i,j,k)];
                        const double v_m = tmp[idx(im1,j,k)];
                        const double v_p = tmp[idx(ip1,j,k)];
                        const double mean        = 0.5 * (v_m + v_p);
                        const double half_spread = 0.5 * std::abs(v_p - v_m);
                        if (std::abs(v - mean) > threshold * half_spread)
                            field.x[i][j][k] = mean;
                    }
                }
            }
        }

    } // passes
}

// ============================================================================
// Shapiro (1-2-1) wiggle-damping filter — Array_2D overload
// ============================================================================
void AtomUtils::damp_wiggles(Array_2D& field,
                              bool along_j, bool along_k,
                              double strength,
                              int passes)
{
    const int jm = field.jm;
    const int km = field.km;
    const double coeff = strength * 0.25;

    // Flat scratch buffer: tmp[j*km + k]
    std::vector<double> tmp(jm * km);

    for (int pass = 0; pass < passes; ++pass) {

        // --- along k (longitude, periodic) ---
        if (along_k) {
            #pragma omp parallel for schedule(static)
            for (int j = 0; j < jm; ++j)
                for (int k = 0; k < km; ++k)
                    tmp[j * km + k] = field.y[j][k];

            #pragma omp parallel for schedule(static)
            for (int j = 0; j < jm; ++j) {
                for (int k = 0; k < km; ++k) {
                    const int km1 = (k - 1 + km) % km;
                    const int kp1 = (k + 1)      % km;
                    field.y[j][k] = tmp[j * km + k]
                        + coeff * (tmp[j * km + km1]
                                 - 2.0 * tmp[j * km + k]
                                 + tmp[j * km + kp1]);
                }
            }
        }

        // --- along j (latitude, clamped at poles) ---
        if (along_j) {
            #pragma omp parallel for schedule(static)
            for (int j = 0; j < jm; ++j)
                for (int k = 0; k < km; ++k)
                    tmp[j * km + k] = field.y[j][k];

            #pragma omp parallel for schedule(static)
            for (int j = 0; j < jm; ++j) {
                const int jm1 = std::max(j - 1, 0);
                const int jp1 = std::min(j + 1, jm - 1);
                for (int k = 0; k < km; ++k) {
                    field.y[j][k] = tmp[j * km + k]
                        + coeff * (tmp[jm1 * km + k]
                                 - 2.0 * tmp[j   * km + k]
                                 + tmp[jp1 * km + k]);
                }
            }
        }

    } // passes
}

// ============================================================================
// Extreme-peak removal filter — Array_2D overload
// ============================================================================
void AtomUtils::remove_peaks(Array_2D& field,
                              bool along_j, bool along_k,
                              double threshold,
                              int passes)
{
    const int jm = field.jm;
    const int km = field.km;

    std::vector<double> tmp(jm * km);

    for (int pass = 0; pass < passes; ++pass) {

        // --- along k (longitude, periodic) ---
        if (along_k) {
            #pragma omp parallel for schedule(static)
            for (int j = 0; j < jm; ++j)
                for (int k = 0; k < km; ++k)
                    tmp[j * km + k] = field.y[j][k];

            #pragma omp parallel for schedule(static)
            for (int j = 0; j < jm; ++j) {
                for (int k = 0; k < km; ++k) {
                    const int km1 = (k - 1 + km) % km;
                    const int kp1 = (k + 1)      % km;
                    const double v   = tmp[j * km + k];
                    const double v_m = tmp[j * km + km1];
                    const double v_p = tmp[j * km + kp1];
                    const double mean        = 0.5 * (v_m + v_p);
                    const double half_spread = 0.5 * std::abs(v_p - v_m);
                    if (std::abs(v - mean) > threshold * half_spread)
                        field.y[j][k] = mean;
                }
            }
        }

        // --- along j (latitude, clamped at poles) ---
        if (along_j) {
            #pragma omp parallel for schedule(static)
            for (int j = 0; j < jm; ++j)
                for (int k = 0; k < km; ++k)
                    tmp[j * km + k] = field.y[j][k];

            #pragma omp parallel for schedule(static)
            for (int j = 0; j < jm; ++j) {
                const int jm1 = std::max(j - 1, 0);
                const int jp1 = std::min(j + 1, jm - 1);
                for (int k = 0; k < km; ++k) {
                    const double v   = tmp[j   * km + k];
                    const double v_m = tmp[jm1 * km + k];
                    const double v_p = tmp[jp1 * km + k];
                    const double mean        = 0.5 * (v_m + v_p);
                    const double half_spread = 0.5 * std::abs(v_p - v_m);
                    if (std::abs(v - mean) > threshold * half_spread)
                        field.y[j][k] = mean;
                }
            }
        }

    } // passes
}

