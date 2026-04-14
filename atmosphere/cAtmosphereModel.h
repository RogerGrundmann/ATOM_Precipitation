#ifndef CATMOSPHEREMODEL_H
#define CATMOSPHEREMODEL_H

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
#include "tinyxml2.h"
#include "Utils.h"
#include "Config.h"

#ifdef _OPENMP
#include <omp.h>
#endif


using namespace std;

class MoistConvection;
class ZeroCatIceScheme;
class OneCatIceScheme;
class TwoCatIceScheme;
class ThreeCatIceScheme;
class SaturationAdjustment;
class VelocityInitializer;
class PressureSolverAtm;
class ThermoAtm;
class UtilsAtm;
class BC_Atm;
class TurbulenceAtm;

class cAtmosphereModel{

    friend class MoistConvection;
    friend class ZeroCatIceScheme;
    friend class OneCatIceScheme;
    friend class TwoCatIceScheme;
    friend class ThreeCatIceScheme;
    friend class SaturationAdjustment;
    friend class VelocityInitializer;
    friend class PressureSolverAtm;
    friend class ThermoAtm;
    friend class UtilsAtm;
    friend class BC_Atm;
    friend class TurbulenceAtm;

public:

    #include "AtmosphereParams.h.inc"

    double re_turb;    // set by TurbulenceAtm::init() as vel_star * z_0 / nue_air

    cAtmosphereModel();
    ~cAtmosphereModel();

    cAtmosphereModel(const cAtmosphereModel&) = delete;
    cAtmosphereModel& operator=(const cAtmosphereModel&) = delete;

    static cAtmosphereModel* get_model(){
        if(!m_model){
            m_model = new cAtmosphereModel();
        }
        return m_model;
    }

    static const double pi180, the_degree, phi_degree, dthe, dphi, dr;
    double dt = 0.0;
    static const double the0, phi0, r0, residuum_ref_atm;

    const int c43 = 4.0/3.0, c13 = 1.0/3.0;

    static const int im = 41, jm = 181, km = 361;

    double residuum_old = 1.0e-5;

    int Ma;
    int n, n_print, n_paraview, panorama, panorama_step;


    struct CellGeometry {                                                                                                                                                                                            
        double rm, rm2, exp_rm, exp_2_rm;                                                                                                                                                                              
        double sinthe, sinthe2, costhe;                                                                                                                                                                                
        double inv_rm, inv_rm2;
        double inv_rmsinthe, inv_rm2sinthe, inv_rm2sinthe2;
        double costhe_inv_rm2sinthe;
        double inv_2dr, inv_2dthe, inv_2dphi;
        double inv_dr2, inv_dthe2, inv_dphi2;
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

    std::map<float,float> m;
    float get_temperatures_from_curve(float time, std::map<float, float>& m) const;

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

private:

    static cAtmosphereModel* m_model;

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

    // Turbulence model parameters
    double zeta     = 3.715;  // coordinate-stretching factor

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

    std::vector<double> radiation_original;                             // original radiation
    std::vector<double> tropopause_layers;                              // keep the tropopause layer index

    std::vector<double> cloud_loc;                                      // lateral cloudwater distribution
    std::vector<double> cloud_max;                                      // lateral cloud_max distribution
    std::vector<double> ice_loc;                                        // lateral ice distribution
    std::vector<double> ice_max;                                        // lateral ice_max distribution

    std::vector<double> c_land_red;
    std::vector<double> c_ocean_red;

    std::vector<double> CAPE;
    std::vector<double> K_u;
    std::vector<double> K_d;

    std::vector<double> lapse_rate;

    std::vector<float> m_layer_heights;

    void SetDefaultConfig();
    void print_min_max_atm();
    void run_3D_loop(int Ma);
    void load_global_temperature_curve();
    void load_equat_temperature_curve();
    void load_pole_temperature_curve();
    void calculate_node_weights();
    void init_steps();
    void init_tropopause_layers();
    void RHS_Atmosphere(int i, int j, int k, const CellGeometry& geo);
    void RHS_Atmosphere_Turb(int i, int j, int k, const CellGeometry& geo);
    void solveRungeKutta_Atmosphere();
    void solveRungeKutta_Atmosphere_Turb();
    void fft(Array &);
    void LandOceanFraction();
    void initTemperatureData(int Ma);

    void initWaterWapour();
    void initCloudIce();
    void init_vapour_cloud();
    void cloudiness_backup();
    void init_Maxwell();

    void init_layer_heights(){
        const float zeta = 3.715;
        float h = L_atm;
        for(int i=0; i<im; i++){
            m_layer_heights.push_back((exp(zeta 
            * (rad.z[i] - 1.0)) - 1) * h);                              // in m      local atmospheric shell thickness
//            std::cout << m_layer_heights.back() << std::endl;
        } 
        return;
    }

    void init_topography(const string &topo_filename);
    void save_data();
    void save_array(const string& fn, const Array& a);
    void BC_pole();
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
    Array_2D Vegetation;                                                // vegetation via precipitation
    Array_2D precipitable_water;                                        // areas of precipitable water at the surface
    Array_2D precipitation_NASA;                                        // surface precipitation by NASA
    Array_2D temperature_NASA;                                          // surface temperature by NASA
    Array_2D velocity_v_NASA;                                           // surface v-velocity by NASA
    Array_2D velocity_w_NASA;                                           // surface w-velocity by NASA
    Array_2D temp_pot;                                                  // surface temperature by barotropic projection
    Array_2D temp_reconst;                                              // surface temperature by reconstuction tool
    Array_2D temp_landscape;                                            // landscape temperature
    Array_2D p_stat_landscape;                                          // landscape static pressure
    Array_2D relative_humidity;                                         // relative_humidity
    Array_2D Q_radiation_2D;                                            // heat from the radiation balance in [W/m2]
    Array_2D Q_latent_2D;                                               // latent heat from bottom values by the energy transport equation
    Array_2D Q_sensible_2D;                                             // sensible heat from bottom values by the energy transport equation
    Array_2D Q_bottom_2D;                                               // difference by Q_radiation - Q_latent - Q_sensible
    Array_2D vapour_evaporation;                                        // water vapour by evaporation in [mm/d]
    Array_2D Evaporation_Dalton;                                        // evaporation by Dalton in [mm/d]
    Array_2D Evaporation_Meyer;                                         // evaporation by Meyer (1915) in [mm/d]
    Array_2D Evaporation_Rohwer;                                        // evaporation by Rohwer (1931) in [mm/d]
    Array_2D Evaporation;                                               // evaporation by active model in [mm/d]
    Array_2D co2_total;                                                 // areas of higher co2 concentration
    Array_2D dew_point_temperature;                                     // dew point temperature
    Array_2D condensation_level;                                        // local condensation level
    Array_2D c_fix;                                                     // local surface water vapour fixed for iterations
    Array_2D Landscape;                                                 // local landscape
    Array_2D Tropopause;                                                // local tropopause
    Array_2D vel_star;                                                  // friction velocity u_tau per (j,k) column [m/s]

    Array h;                                                            // bathymetry, depth from sea level
    Array t;                                                            // temperature
    Array u;                                                            // u-component velocity component in r-direction
    Array v;                                                            // v-component velocity component in theta-direction
    Array w;                                                            // w-component velocity component in phi-direction
    Array c;                                                            // water vapour
    Array cloud;                                                        // cloud water
    Array ice;                                                          // cloud ice
    Array gr;                                                           // cloud graupel
    Array co2;                                                          // CO2
    Array tn;                                                           // temperature new
    Array un;                                                           // u-velocity component in r-direction new
    Array vn;                                                           // v-velocity component in theta-direction new
    Array wn;                                                           // w-velocity component in phi-direction new
    Array cn;                                                           // water vapour new
    Array cloudn;                                                       // cloud water new
    Array icen;                                                         // cloud ice new
    Array grn;                                                          // cloud ice new
    Array co2n;                                                         // CO2 new
    Array p_dyn;                                                        // dynamic pressure
    Array p_stat;                                                       // static pressure
    Array TempStand;                                                    // US Standard Atmosphere Temperature
    Array TempDewPoint;                                                 // Dew Point Temperature
    Array HumidityRel;                                                  // relative humidity
    Array rhs_t;                                                        // auxilliar field RHS temperature
    Array rhs_u;                                                        // auxilliar field RHS u-velocity component
    Array rhs_v;                                                        // auxilliar field RHS v-velocity component
    Array rhs_w;                                                        // auxilliar field RHS w-velocity component
    Array rhs_c;                                                        // auxilliar field RHS water vapour
    Array rhs_cloud;                                                    // auxilliar field RHS cloud water
    Array rhs_ice;                                                      // auxilliar field RHS cloud ice
    Array rhs_g;                                                        // auxilliar field RHS cloud graupel
    Array rhs_co2;                                                      // auxilliar field RHS CO2
    Array aux_u;                                                        // auxilliar field u-velocity component
    Array aux_v;                                                        // auxilliar field v-velocity component
    Array aux_w;                                                        // auxilliar field w-velocity component
    Array aux_t;                                                        // auxilliar field t
    Array Q_Latent;                                                     // latent heat
    Array Q_Sensible;                                                   // sensible heat
    Array BuoyancyForce;                                                // buoyancy force, Boussinesque approximation
    Array CoriolisForce;                                                // Coriolis force terms
    Array CentrifugalForce;                                             // centrifugal force terms
    Array PresGradForce;                                                // Force caused by normal pressure gradient
    Array epsilon;                                                      // emissivity/ absorptivity
    Array radiation;                                                    // radiation
    Array P_rain;                                                       // rain precipitation mass rate
    Array P_snow;                                                       // snow precipitation mass rate
    Array P_rainn;                                                      // rain precipitation mass rate   new
    Array P_snown;                                                      // snow precipitation mass rate   new
    Array P_graupel;                                                    // graupel precipitation mass rate
    Array P_conv;                                                       // rain formation by deep-level cloud convection
    Array Precipitation;                                                // areas of higher precipitation
    Array PrecipitableWaterLocal;                                       // precipitable water at each level
    Array S_v;                                                          // water vapour mass rate
    Array S_c;                                                          // cloud water mass rate
    Array S_i;                                                          // cloud ice mass rate
    Array S_r;                                                          // rain mass rate
    Array S_s;                                                          // snow mass rate
    Array S_g;                                                          // graupel mass rate
    Array S_c_c;                                                        // cloud water mass rate due to condensation and evaporation in the saturation adjustment technique
    Array M_u;                                                          // moist convection within the updraft
    Array M_d;                                                          // moist convection within the downdraft
    Array MC_t;                                                         // moist convection acting on dry static energy
    Array MC_q;                                                         // moist convection acting on water vapour development
    Array MC_v;                                                         // moist convection acting on v-velocity component
    Array MC_w;                                                         // moist convection acting on w-velocity component
    Array r_dry;                                                        // density of dry air
    Array r_humid;                                                      // density of humid air
    Array g_p;                                                          // conversion cloud droplets to raindrops
    Array c_u;                                                          // condensation in the updraft
    Array e_d;                                                          // evaporation of precipitation in the downdraft
    Array e_l;                                                          // evaporation of cloud water in the environment
    Array e_p;                                                          // evaporation of cloud water in the environment
    Array s;                                                            // dry static energy
    Array s_u;                                                          // dry static energy in the updraft
    Array s_d;                                                          // dry static energy in the downdraft
    Array u_u;                                                          // u-velocity component in the updraft
    Array u_d;                                                          // u-velocity component in the downdraft
    Array v_u;                                                          // u-velocity component in the updraft
    Array v_d;                                                          // u-velocity component in the downdraft
    Array w_u;                                                          // u-velocity component in the updraft
    Array w_d;                                                          // u-velocity component in the downdraft
    Array q_v_u;                                                        // water vapour in the updraft
    Array q_v_d;                                                        // water vapour in the downdraft
    Array q_c_u;                                                        // cloud water in the updraft
    Array E_u;                                                          // moist entrainment in the updraft
    Array D_u;                                                          // moist detrainment in the updraft
    Array E_d;                                                          // moist entrainment in the downdraft
    Array D_d;                                                          // moist detrainment in the downdraft

    Array CloudBase;                                                    // cloud base
    Array LevelFreeSinking;                                             // level of free sinking

    Array Deep_beg;                                                     // cloud base
    Array Deep_end;                                                     // level of free sinking

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
    Array p_hydro;                                                      // hydrostatic pressure
    Array PressureGradientForce;                                        // pressure gradient force magnitude
};
#endif
