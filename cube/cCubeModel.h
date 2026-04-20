#ifndef CCUBEMODEL_H
#define CCUBEMODEL_H

#include <fenv.h>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <ctime>    
#include <cmath>
#include <map>
#include <set>
#include <csignal>
#include <cstring>
#include <limits>
#include <functional>
#include <stdexcept>
#include <sys/stat.h>
#include <sys/types.h>
#include "Array.h"
#include "Array_1D.h"
#include "Array_2D.h"
#include "Utils.h"

#ifdef _OPENMP
#include <omp.h>
#endif


using namespace std;

class PressureSolverCube;
class UtilsCube;
class BC_Cube;
class TurbulenceCube;

class cCubeModel{

    friend class PressureSolverCube;
    friend class UtilsCube;
    friend class BC_Cube;
    friend class TurbulenceCube;

public:

    double re_turb;    // set by TurbulenceAtm::init() as vel_star * z_0 / nue_air

    cCubeModel();
    ~cCubeModel();

    cCubeModel(const cCubeModel&) = delete;
    cCubeModel& operator=(const cCubeModel&) = delete;

    static cCubeModel* get_model(){
        if(!m_model){
            m_model = new cCubeModel();
        }
        return m_model;
    }

    static const double dy, dx, dr;     // uniform grid spacing in y, x; stretched in z (dr)
    double dt = 0.0;
    static const double y0, x0, r0, residuum_ref_atm;

    const int c43 = 4.0/3.0, c13 = 1.0/3.0;

/*
    static const int im = 41, jm = 21, km = 81, nm = 2000;

    static const int cube_i_beg = 0, cube_j_beg = 8, cube_k_beg = 15;
    static const int cube_i_end = 8, cube_j_end = 12, cube_k_end = 25;
*/
    static const int im = 41, jm = 101, km = 241;

//    static const int cube_i_beg = 0, cube_j_beg = 16, cube_k_beg = 40;
//    static const int cube_i_end = 20, cube_j_end = 24, cube_k_end = 60;

//    static const int cube_i_beg = 0, cube_j_beg = 10, cube_k_beg = 40;
//    static const int cube_i_end = 20, cube_j_end = 30, cube_k_end = 50;

//    static const int cube_i_beg = 0, cube_j_beg = 15, cube_k_beg = 40;
//    static const int cube_i_end = 10, cube_j_end = 25, cube_k_end = 50;

    static const int cube_i_beg = 0, cube_j_beg = 45, cube_k_beg = 80;
    static const int cube_i_end = 10, cube_j_end = 55, cube_k_end = 90;

    double residuum_old = 1.0e-5;

    int Ma;
    int n, n_print, n_paraview, panorama, panorama_step;


    struct CellGeometry {
        // Grid-spacing reciprocals — uniform in all three directions.
        double inv_2dz, inv_2dy, inv_2dx;
        double inv_dz2, inv_dy2, inv_dx2;
    };


    bool is_final_result = false;
    bool is_print_mode = false;

    bool is_first_time_slice() const{
        if(m_time_list.empty()){
            throw("The time list is empty. It is likely the model has not started yet.");
        }
        return (m_current_time == m_time_list.begin());
    }

    std::vector<std::vector<int> > i_topography;
    std::vector<std::vector<int> > i_tropopause;
    std::vector<std::vector<int> > i_landscape;

    std::map<float,float> temperature_curve;
    float get_temperatures_from_curve(float time, std::map<float, float>& curve) const;

    std::set<float>::const_iterator get_current_time() const{
        if(m_time_list.empty()){
            throw("The time list is empty. It is likely the model has not started yet.");
        }else{
            return m_current_time;
        }
    }

    std::set<float>::const_iterator get_previous_time() const{
        if(m_time_list.empty()){
            throw("The time list is empty. It is likely the model has not started yet.");
        }
        if(m_current_time != m_time_list.begin()){
            std::set<float>::const_iterator ret = m_current_time;
            ret--;
            return ret;
        }
        else{
            throw("The current time is the only time slice for now. There is no previous time yet.");
        }
    }

    /*
     * Given a latitude, return the layer index of tropopause
    */
    int get_tropopause_layer(int j){
        assert(j>=0);
        assert(j<jm);
        //refer to  BC_Thermo::TropopauseLocation and BC_Thermo::GetTropopauseHightAdd
        //tropopause height is proportional to the mean tropospheric temperature.
        //higher near the equator - warm troposphere
        //lower at the poles - cold troposphere
        return tropopause_layers[j];
    }
    /*
     *
    */
    int get_surface_layer(int j, int k){
        return i_topography[j][k];
    }
    /*
     * This function must be called after init_layer_heights()
     * Given a layer index i, return the height of this layer
    */
    float get_layer_height(int i){
        if(0>i || i>im-1){
            return -1;
        }
        return m_layer_heights[i];
    }
    std::vector<float> get_layer_heights(){
        return m_layer_heights;
    }   
    /*
    * Given a altitude, return the layer index
    */
    int get_layer_index(float height){
        std::size_t i = 0;
        for(; i<m_layer_heights.size(); i++){
            if(height<m_layer_heights[i])
                return i-1;
        }
        return i;
    }


    void LoadConfig(const char *filename);
    void Run();
    void RunTimeSlice(int time_slice);

    // Generated configuration parameters (from param.py → CubeParameters.h)
    #include "CubeParameters.h"

private:

    void SetDefaultConfig();

    static cCubeModel* m_model;

    string at = "AGCM";
    string bathymetry_name;

    bool has_printed_welcome_msg;

    bool is_global_temperature_curve_loaded(){
        return !m_global_temperature_curve.empty();
    }
    bool is_equat_temperature_curve_loaded(){
        return !m_equat_temperature_curve.empty();
    }
    bool is_pole_temperature_curve_loaded(){
        return !m_pole_temperature_curve.empty();
    }

    int panorama_cnt, iter_n;


    double t_paleo_total = 0.0;
    double t_pole_total = 0.0;
    double t_global_mean = 0.0;

    // Turbulence model selection flags
    bool use_turbulence_model          = false;
    bool use_k_epsilon_turbulence_model    = false;
    bool use_k_omega_turbulence_model      = false;
    bool use_k_omega_SST_turbulence_model  = false;
    bool use_stretched_coordinate_system   = false;

    // SST blending: inner (zone 1, near wall) and outer (zone 2, free stream)
    static double blend(double inner, double outer, double F1) {
        return F1 * inner + (1.0 - F1) * outer;
    }

    std::set<float> m_time_list;
    std::set<float>::const_iterator m_current_time;

    std::map<float,float> m_global_temperature_curve;
    std::map<float,float> m_equat_temperature_curve;
    std::map<float,float> m_pole_temperature_curve;

    std::vector<double> alfa;
    std::vector<double> beta;

    std::vector<std::vector<double> > M_u_Base;
    std::vector<std::vector<double> > M_d_LFS;

    std::vector<std::vector<int> > i_Base;
    std::vector<std::vector<int> > i_LFS;

    std::vector<std::vector<int> > i_deep_beg;
    std::vector<std::vector<int> > i_deep_end;

    std::vector<double> short_wave_radiation;                           // lateral short wave radiation

    std::vector<double> tropopause_layers;                              // keep the tropopause layer index

    std::vector<float> m_layer_heights;

    void print_min_max_atm();
    void run_3D_loop();
    void calculate_node_weights();
    void init_steps();
    void init_tropopause_layers();
    void RHS_Cube_Turb(int i, int j, int k, const CellGeometry& geo);
    void solveRungeKutta_Cube();
    void solveRungeKutta_Cube_Turb();
    void fft(Array &);

    void init_layer_heights(){
        for(int i=0; i<im; i++)
            m_layer_heights.push_back(float(i) * float(dr) * float(L_atm));
    }

    void init_topography(const string &topo_filename);
    void save_data();
    void save_array(const string& fn, const Array& a);
    void print_welcome_msg();
    void print_final_remarks();
    void print_loop_3D_headings();
    void paraview_panorama_vts(string &Name_Bathymetry_File, int n);
    void paraview_sphere_vts(string &Name_Bathymetry_File, int n);
    void paraview_vtk_radial(string &Name_Bathymetry_File, int i_radial, int n);
    void paraview_vtk_zonal(string &Name_Bathymetry_File, int k_zonal, int n);
    void paraview_vtk_longal(string &Name_Bathymetry_File, int j_longal, int n); 
    void AtmospherePlotData(const string &Name_Bathymetry_File);
    void AtmosphereDataTransfer(const string &Name_Bathymetry_File);
    void read_Atmosphere_Surface_Data(int Ma);
    void searchMinMax_2D(const string &, const string &,
        const string &, Array_2D &, double coeff = 1.0);

    void searchMinMax_3D(const string &, const string &,
        const string &, Array &, double coeff = 1.0,
        std::function< double(double) > lambda = [](double i) -> double{return i;},
        bool print_heading = false);

public:
    Array_1D rad;                                                       // radial coordinate direction
    Array_1D the;                                                       // lateral coordinate direction
    Array_1D phi;                                                       // longitudinal coordinate direction
    Array_1D aux_grad_v;                                                // auxilliar array

    Array_2D Topography;                                                // topography
    Array_2D value_top;                                                 // auxiliar topography
    Array_2D vel_star;                                                  // friction velocity u_τ [m/s], per (j,k) column

    Array h;                                                            // bathymetry, depth from sea level
    Array u;                                                            // u-component velocity component in r-direction
    Array v;                                                            // v-component velocity component in theta-direction
    Array w;                                                            // w-component velocity component in phi-direction
    Array un;                                                           // u-velocity component in r-direction new
    Array vn;                                                           // v-velocity component in theta-direction new
    Array wn;                                                           // w-velocity component in phi-direction new
    Array p_stat;                                                       // static pressure
    Array p_dyn;                                                        // dynamic pressure
    Array rhs_u;                                                        // auxilliar field RHS u-velocity component
    Array rhs_v;                                                        // auxilliar field RHS v-velocity component
    Array rhs_w;                                                        // auxilliar field RHS w-velocity component
    Array aux_u;                                                        // auxilliar field u-velocity component
    Array aux_v;                                                        // auxilliar field v-velocity component
    Array aux_w;                                                        // auxilliar field w-velocity component
    Array aux_t;                                                        // auxilliar field t

    // Turbulence arrays
    Array tke;                                                          // turbulent kinetic energy
    Array tken;                                                         // turbulent kinetic energy (new time level)
    Array dis;                                                          // dissipation rate (ε or specific dissipation ω)
    Array disn;                                                         // dissipation rate (new time level)
    Array nue;                                                          // turbulent eddy viscosity
    Array prod;                                                         // turbulence production P_k
    Array tke_source;                                                   // RHS source term for k equation
    Array dis_source;                                                   // RHS source term for ε/ω equation
    Array rhs_tke;                                                      // auxiliary RHS field for tke
    Array rhs_dis;                                                      // auxiliary RHS field for dis
};
#endif
