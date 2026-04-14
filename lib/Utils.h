#ifndef _UTILS_
#define _UTILS_

#include <string>
#include <set>
#include <map>
#include <vector>

#include <iostream>
#include <fstream>
#include <limits>

#include "Array.h"
#include "Array_1D.h"
#include "Array_2D.h"

#define logger() \
if (false) ; \
else get_logger()

namespace AtomUtils{
    using namespace std;

    struct HemisphereCoords{
        double lat, lon;
        string east_or_west, north_or_south;
    }; 

    HemisphereCoords convert_coords(double lon, double lat);

    std::ofstream& get_logger();

    std::tuple<double, int, int, int>
        max_diff(int i, int j, int k, const Array &a1, const Array &a2);

    inline bool is_land(const Array& h, int i, int j, int k){
        return fabs(h.x[i][j][k] - 1) < std::numeric_limits<double>::epsilon();
    }
    inline bool is_air(const Array& h, int i, int j, int k){
        return !is_land(h,i,j,k);
    }
    inline bool is_water(const Array& h, int i, int j, int k){
        return !is_land(h,i,j,k);
    }
    inline bool is_ocean_surface(const Array& h, int i, int j, int k){
        return i==0 && !is_land(h, i, j, k);
    }
    inline bool is_land_surface(const Array& h, int i, int j, int k){
        return is_land(h, i, j, k) && !is_land(h, i+1, j, k);
    }
    inline bool is_east_coast(const Array& h, int j, int k){
        if(k == h.get_km()-1 ) return false;//on grid boundary
        return is_land(h, 0, j, k) && !is_land(h, 0, j, k+1);
    }
    inline bool is_west_coast(const Array& h, int j, int k){
        if( k == 0 ) return false;//on grid boundary
        return is_land(h, 0, j, k) && !is_land(h, 0, j, k-1);
    }
    inline bool is_north_coast(const Array& h, int j, int k){
        if( j == 0 ) return false;//on grid boundary
        return is_land(h, 0, j, k) && !is_land(h, 0, j-1, k);
    }
    inline bool is_south_coast(const Array& h, int j, int k){
        if(j == h.get_jm()-1 ) return false;//on grid boundary
        return is_land(h, 0, j, k) && !is_land(h, 0, j+1, k);
    }
    inline double parabola(double x){
        return x * x - 2.0 * x;
    }
    // x^2-2x mapped to [lower_bound, upper_bound]:
    // x=1 → lower_bound, x=0 or 2 → upper_bound
    inline double parabola_interp(double lower_bound, double upper_bound, double x){
        return upper_bound + (upper_bound - lower_bound) * (x * x - 2.0 * x);
    }
    inline double cosine_profile(double ratio){
    // ratio 0 -> Pol (Süd), ratio 1 -> Äquator, ratio 2 -> Pol (Nord)
    // Wir mappen ratio [0, 2] auf Cosinus-Argument [-PI/2, PI/2]
    // Ein einfacherer Weg für das gleiche Profil:
    // f(1) = 1 (Äquator), f(0) = 0 (Pol), f(2) = 0 (Pol)
    
//        return sin(M_PI_2 * ratio) * (2.0 - ratio); // Alternative Formel
    // ODER klassisch:
        return cos((ratio - 1.0) * M_PI_2) - 1.0; 
    }    
    inline double Agnesi(double a, double x){
        return pow(a, 3.0) // Versiera di Agnesi
            /(pow(a, 2.0) + pow(x, 2.0));
    }
    inline double exp_func(double T_K, const double co_1, const double co_2){
        return exp(co_1 * (T_K - 273.15) / (T_K - co_2));  // temperature in °K
    }


    int lon_2_index(double lon);
    int lat_2_index(double lat);
    int RunStart(string comment);
    int RunEnd(string comment, int Ma, int Run_start);

    double C_Dalton(int i, int j, int k, double coeff_Dalton, 
        double u_0, Array &u, Array &v, Array &w);
    double simpson(int n1, int n2, double dstep, Array_1D &value);
    double trapezoidal(int n1, int n2, double dstep, Array_1D &value);
    double rectangular(int n1, int n2, double dstep, Array_1D &value);
    double GetMean_3D(int jm, int km, Array &val_3D);
    double GetMean_2D(int jm, int km, Array_2D &val_2D);

    void read_IC(const string& fn, double** a, int jm, int km);
    void smooth_steps(int k, int im, int jm, Array &value);
    void smooth_tropopause(int jm, std::vector<double> &value);
    void load_map_from_file(const std::string& fn, std::map<float, float>& m);
    void CalculateNodeWeights(int jm, int km);
    void calculate_node_weights();
    //do the fft in all 3 directions
    void fft_gaussian_filter_3d(Array& data, int sigma);
    // do the fft in one direction
    // the valid directions are 'i', 'j' and 'k'
    // sigma is the Standard deviation for Gaussian kernel and controls how blurry the result will be
    //the result will replace the input array
    void fft_gaussian_filter_3d(Array& data, int sigma, char direction);
    void mirror_padding(double* data, size_t i_len, size_t p_len);
    void fft_gaussian_filter(double* _data, double* kernel, size_t len);
    //change data coordinate system from -180° _ 0° _ +180° to 0°- 360°
    void move_data(double* data, int len);
    void move_data(std::vector<int>& data, int len);
    void move_data_to_new_arrays(int im, int jm, int km, double coeff, 
        std::vector<Array*>& arrays, std::vector<Array*>& new_arrays);
    void move_data_to_new_arrays( int jm, int km, double coeff, 
        std::vector<Array*>& arrays, 
        std::vector<Array*>& new_arrays, int i = 0);
    void use_presets(bool preset_1, bool preset_2, bool preset_3);





    // Shapiro (1-2-1) wiggle-damping filter for an Array field.
    //
    // Each enabled axis gets one filter pass per call to damp_wiggles():
    //   f_new[n] = f[n] + strength * 0.25 * (f[n-1] - 2*f[n] + f[n+1])
    //
    // Axis boundary conventions:
    //   k (longitude) — periodic wrap
    //   j (latitude)  — clamped (no-flux) at poles
    //   i (vertical)  — clamped at top/bottom
    //
    // i_surface: pointer to the 2D fluid-domain floor index array.
    //   Accepts either i_topography[j][k] (atmosphere: first cell above terrain)
    //   or i_bathymetry[j][k] (hydrosphere: first cell above the seafloor).
    //   In both cases cells with i < i_surface[j][k] are skipped (solid body).
    //   Pass nullptr to skip the mask and filter all cells.
    //
    // strength ∈ [0, 1]: 0 = no change, 1 = full Shapiro step.
    // passes: number of times the filter sequence is repeated.
    void damp_wiggles(Array& field,
                      const std::vector<std::vector<int>>* i_surface,
                      bool along_i, bool along_j, bool along_k,
                      double strength = 1.0,
                      int passes = 1);

    // Extreme-peak removal filter for an Array field.
    //
    // Along each enabled axis the two direct neighbours of every cell define a
    // local background (their mean) and a variability scale (half-spread,
    // |f[n+1] - f[n-1]| / 2).  If the cell's deviation from the background
    // exceeds threshold * half_spread the cell is replaced by the background.
    // Smooth gradients are left intact; only isolated spikes are removed.
    //
    // Axis boundary conventions and i_surface masking are identical to
    // damp_wiggles(): k periodic, j/i clamped, solid cells skipped.
    //
    // threshold: sensitivity — lower values remove more; typical range 1–4.
    // passes:    number of times the scan is repeated per call.
    void remove_peaks(Array& field,
                      const std::vector<std::vector<int>>* i_surface,
                      bool along_i, bool along_j, bool along_k,
                      double threshold = 2.0,
                      int passes = 1);

    // 2-D overloads (Array_2D, axes j and k only).
    // Boundary conventions: k periodic, j clamped at poles.
    void damp_wiggles(Array_2D& field,
                      bool along_j, bool along_k,
                      double strength = 1.0,
                      int passes = 1);
    void remove_peaks(Array_2D& field,
                      bool along_j, bool along_k,
                      double threshold = 2.0,
                      int passes = 1);

    template<class T>
    void set_values(T* a, T value, int len){
        for(int i = 0 ; i < len; i++){
            a[i] = value;
        }
    }
    template<class T>
    void chessboard_grid(T** a, int j, int k, int jm, int km){
        T value_1 = 0.9, value_2 = 1.1;
        T* line_1 = new T[km];
        T* line_2 = new T[km];
        int i = 0;
        bool flip = true;
        while(i+k < km){
            if(flip){
                set_values(line_1+i, value_1, k);
                set_values(line_2+i, value_2, k);
            }else{
                set_values(line_1+i, value_2, k);
                set_values(line_2+i, value_1, k);
            }
            i+=k;
            flip = !flip;
        }
        if(i<km-1){
            if(flip){
                set_values(line_1+i, value_1, km-1-i);
                set_values(line_2+i, value_2, km-1-i);
            }else{
                set_values(line_1+i, value_2, km-1-i);
                set_values(line_2+i, value_1, km-1-i);
            }
        }
        int cnt=0;
        flip = true;
        for(i = 0; i<jm; i++){
            if(flip)
                memcpy(a[i], line_1, km*sizeof(T));
            else
                memcpy(a[i], line_2, km*sizeof(T));
            cnt++;
            if(cnt==j){
                flip = !flip;
                cnt = 0;
            }
        }
        delete[] line_1;
        delete[] line_2;
        return;
    }
}

#endif
