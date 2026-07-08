#ifndef CHYDROSPHEREMODEL_H
#define CHYDROSPHEREMODEL_H

#include <fenv.h>
#include <algorithm>
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
#include <cstdio>
#include <cstdlib>
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

class BoundaryConditionHyd;        // forward declaration
class PressureSolverHyd;           // forward declaration
class ThermoHyd;                   // forward declaration
class ThermoHalineConveyorBelt;    // forward declaration
class UtilsHyd;                    // forward declaration
class TurbulenceHyd;               // forward declaration

class cHydrosphereModel{

public:

    #include "HydrosphereParams.h.inc"

    friend class BoundaryConditionHyd;
    friend class PressureSolverHyd;
    friend class ThermoHyd;
    friend class ThermoHalineConveyorBelt;
    friend class UtilsHyd;
    friend class TurbulenceHyd;

    cHydrosphereModel();
    ~cHydrosphereModel();
 
    cHydrosphereModel(const cHydrosphereModel&) = delete;
    cHydrosphereModel& operator=(const cHydrosphereModel&) = delete;

    static cHydrosphereModel* get_model(){
        if(!m_model){
            m_model = new cHydrosphereModel();
        }
        return m_model;
    }

    static const double pi180, the_degree, phi_degree, dthe, dphi, dr;
    double dt = 0.0;
    static const double the0, phi0, r0, residuum_ref_hyd;
    static const double dr_stretch;                                         // sinh stretching parameter for radial coordinate (cosh(dr_stretch) = step ratio surface/bottom)

    const int c43 = 4.0/3.0, c13 = 1.0/3.0;

    static const int im = 41, jm = 181, km = 361;

    double residuum_old = 1.0e-5;

    int Ma;



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

    // Turbulence model flags — set by TurbulenceHyd::init() based on turb_model string
    bool use_turbulence_model              = false;
    bool use_k_epsilon_turbulence_model    = false;
    bool use_k_omega_turbulence_model      = false;
    bool use_k_omega_SST_turbulence_model  = false;
    double re_turb = 1.0;    // set by TurbulenceHyd::init() as vel_star * z_0 / nue_water
    double pr_turb = 0.9;    // turbulent Prandtl number for heat/salinity transport

    bool is_first_time_slice() const{
        if(m_time_list.empty()){
            throw("The time list is empty. It is likely the model has not started yet.");
        }
        return (m_current_time == m_time_list.begin());
    }

    std::vector<std::vector<int> > i_bathymetry;
    std::vector<std::vector<int> > wind;

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

    void LoadConfig(const char *filename);
    void Run();
    void RunTimeSlice(int time_slice);

private:

    static cHydrosphereModel* m_model;


    string bathymetry_name;
    string hy = "OGCM";

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

    std::set<float> m_time_list;
    std::set<float>::const_iterator m_current_time;

    std::map<float,float> m_global_temperature_curve;
    std::map<float,float> m_equat_temperature_curve;
    std::map<float,float> m_pole_temperature_curve;

    void SetDefaultConfig();
    void load_temperature_curve();
    void save_data();
    void IC_t_WestEastCoast();
    void IC_u_WestEastCoast();
    void IC_Equatorial_Currents();
    void CircumPolarCurrent();
    void solveRungeKutta_Hydrosphere_Turb();   // single dynamical core (laminar RHS_Hydrosphere dropped 2026-07-08)
    void RHS_Hydrosphere_Turb(int i, int j, int k, const CellGeometry& geo);
    void print_welcome_msg();
    void print_final_remarks();
    void print_loop_3D_headings();
    void print_min_max_hyd();
    bool scanForNaN_hyd(int iter, const char *stage);
    void LandOceanFraction();
    void init_bathymetry(const string &bathymetry_file);
    void BC_Surface_Salinity_NASA(const string &Name_SurfaceSalinity_File);
    void BC_Surface_Temperature_NASA
       (const string &Name_SurfaceTemperature_File);
    void EkmanSpiral();
    void initTemperature(int Ma);
    void initSalinity();
    void Pressure_Limitation_Hyd();
    void paraview_vts(const string &Name_Bathymetry_File, int n);
    void paraview_panorama_vts(const string &Name_Bathymetry_File, int n);
    void paraview_sphere_vts(const string &Name_Bathymetry_File, int n);
    void paraview_vtk_zonal(const string &Name_Bathymetry_File, int k_zonal, int n);
    void paraview_vtk_radial(const string &Name_Bathymetry_File, int i_radial, int n);
    void paraview_vtk_longal(const string &Name_Bathymetry_File, int j_longal, int n);
    void HydrospherePlotData(const string &Name_Bathymetry_File);
    void HydrosphereDataTransfer(const string &Name_Bathymetry_File);
    void run_3D_loop();
    std::vector<Array*> restart_arrays();
    void save_state(int iter);
    bool load_state(int iter);
    void zonal_mean_w(std::vector<std::vector<double> >& wbar);
    void write_w_momentum_budget(int iter);
    void write_deep_momentum_budget(int iter);
    void locate_blowup_cell();
    void record_stage(int s);
    void write_stage_budget(int iter);
    void searchMinMax_2D(string, string, 
        string, Array_2D &, double coeff=1.);
    void searchMinMax_3D(string, string, 
        string, Array &, double coeff=1.,
        std::function< double(double) > lambda = [](double i)->double{return i;},
        bool print_heading=false);
    void load_global_temperature_curve();
    void load_equat_temperature_curve();
    void load_pole_temperature_curve();
    void read_Hydrosphere_Surface_Data(int Ma);

public:

    // Inviscid spin-up state (mirror of cAtmosphereModel — see RHS_Hyd.cpp and BC_Hyd.h).
    int total_iter_count = 0;
    double diffusion_ramp = 1.0;
    bool inviscid_phase = false;

    // Zonal-mean zonal-velocity (w) momentum-budget diagnostic — attributes the
    // SH surface westward spin-down. When wbudget_capture is set (checkpoint
    // iters), RHS_Hyd stores each rhs_w contribution into wbud_*; the run loop
    // snapshots wbar_before and write_w_momentum_budget emits the per-latitude
    // CSV. See project_hydro_ekman_sh_gyre.
    bool wbudget_capture = false;
    std::vector<std::vector<double> > wbar_before;

    // Stage-resolved velocity-change decomposition at the deep blow-up cell.
    // The radial-momentum budget showed the momentum RHS is net-DAMPING at the
    // near-pole runaway cell, so a step AFTER the RK4 must inject the growth.
    // On checkpoint iters (stage_capture), locate_blowup_cell() picks the
    // max-|u| interior ocean cell and record_stage() snapshots u,v,w at it after
    // each iteration stage; write_stage_budget attributes the per-iteration net
    // Δ(u,v,w) to RK4 / projection / BC+clamp / polar-filter. See
    // project_hydro_polar_blowup.
    bool stage_capture = false;
    int  scell_i = -1, scell_j = -1, scell_k = -1;
    double su[5] = {0}, sv[5] = {0}, sw[5] = {0};   // velocity [nondim] after stages 0..4

    Array_1D rad;                                                       // radial coordinate direction
    Array_1D the;                                                       // lateral coordinate direction
    Array_1D phi;                                                       // longitudinal coordinate direction

    // Per-i non-uniform 3-point radial finite-difference coefficients, derived
    // from the actual (stretched) rad.z spacing. The radial grid is stretched
    // (sinh, fine at the surface i=im-1, coarse at the bottom i=0), but the FD
    // stencils were written for a uniform grid (constant dr). These coefficients
    // make every radial derivative correct on the stretched grid. Filled once by
    // setupRadialStencilCoeffs() after rad is built; const for the whole run.
    //   central (points i-1,i,i+1), valid i in [1,im-2]
    std::vector<double> rc1m, rc10, rc1p;   // 1st derivative
    std::vector<double> rc2m, rc20, rc2p;   // 2nd derivative
    //   forward (points i,i+1,i+2), valid i in [0,im-3]
    std::vector<double> rf10, rf11, rf12;   // 1st derivative
    std::vector<double> rf20, rf21, rf22;   // 2nd derivative
    void setupRadialStencilCoeffs();
    Array_1D aux_grad_v;                                                // auxilliar array
    Array_1D aux_grad_w;                                                // auxilliar array
    
    Array_2D Bathymetry;                                                // Bathymetry in m
    Array_2D value_top;                                                 // auxiliar field for bathymetzry
    Array_2D Upwelling;                                                 // upwelling
    Array_2D Downwelling;                                               // downwelling
    Array_2D EkmanPumping;                                              // 2D bottom water summed up in a vertical column
    Array_2D BuoyancyForce_2D;                                          // radiation balance at the surface
    Array_2D salinity_evaporation;                                      // additional salinity by evaporation
    Array_2D Evaporation_Dalton;                                        // evaporation by Dalton in [mm/d]
    Array_2D Evaporation_Penman;                                        // evaporation by Penman in [mm/d]
    Array_2D Evaporation;                                               // evaporation transferred from atmosphere in [mm/d]
    Array_2D Precipitation_2D;                                          // areas of higher precipitation
    Array_2D precipitation_NASA;                                        // surface precipitation from NASA
    Array_2D temperature_NASA;                                          // surface temperature from NASA
    Array_2D temp_landscape;                                            // landscape temperature
    Array_2D temp_reconst;                                              // surface temperature from reconstuction tool
    Array_2D c_fix;                                                     // local surface salinity fixed for iterations
    Array_2D t_surf_fix;                                                // prescribed surface SST, fixed for iterations (sustained surface heat-flux forcing)
    Array_2D v_wind;                                                    // v-component of surface wind
    Array_2D w_wind;                                                    // w-component of surface wind
    Array_2D velocity_v_NASA;                                           // surface v-velocity from NASA
    Array_2D velocity_w_NASA;                                           // surface w-velocity from NASA

    Array h;                                                            // bathymetry, depth from sea level
    Array t;                                                            // temperature
    Array u;                                                            // u-component velocity component in r-direction
    Array v;                                                            // v-component velocity component in theta-direction
    Array w;                                                            // w-component velocity component in phi-direction
    Array c;                                                            // water vapour

    Array tn;                                                           // temperature new
    Array un;                                                           // u-velocity component in r-direction new
    Array vn;                                                           // v-velocity component in theta-direction new
    Array wn;                                                           // w-velocity component in phi-direction new
    Array cn;                                                           // water vapour new

    Array p_dyn;                                                        // dynamic pressure
    Array p_hydro;                                                      // static pressure

    Array rhs_t;                                                        // auxilliar field RHS temperature
    Array rhs_u;                                                        // auxilliar field RHS u-velocity component
    Array rhs_v;                                                        // auxilliar field RHS v-velocity component
    Array rhs_w;                                                        // auxilliar field RHS w-velocity component
    Array rhs_c;                                                        // auxilliar field RHS water vapour

    Array aux_u;                                                        // auxilliar field u-velocity component
    Array aux_v;                                                        // auxilliar field v-velocity component
    Array aux_w;                                                        // auxilliar field w-velocity component

    // Rhie-Chow face mass fluxes (divergence-free transporting velocity).
    // uf.x[i][j][k] = radial flux on face i+1/2, vf on face j+1/2, wf on face
    // k+1/2 (last face index in each direction unused). Built by
    // PressureSolverHyd::project_velocity so continuity is enforced on faces
    // (the collocated cell-centre projection cannot — div/grad central vs the
    // compact Poisson Laplacian decouple odd/even -> checkerboard). Step 2 uses
    // these as the advecting velocity in RHS_Hyd. See project_hydro_continuity_checkerboard.
    Array uf;                                                           // radial face flux (face i+1/2)
    Array vf;                                                           // meridional face flux (face j+1/2)
    Array wf;                                                           // zonal face flux (face k+1/2)

    Array Salt_Finger;                                                  // salt bulge of higher density
    Array Salt_Diffusion;                                               // salt bulge of lowerer density and temperature
    Array Salt_Balance;                                                 // +/- salt balance

    Array r_water;                                                      // water density as function of pressure
    Array r_salt_water;                                                 // salt water density as function of pressure and temperature
    Array BuoyancyForce;                                                // 3D buoyancy force
    Array CoriolisForce;                                                // Coriolis force terms
    Array CentrifugalForce;                                             // Coriolis force terms
    Array PresGradForce;                                        // Force caused by normal pressure gradient

    // Turbulence fields (allocated with the main arrays; zero unless turbulence is active)
    Array tke;                                                          // turbulent kinetic energy k*
    Array tken;                                                         // tke at time level n
    Array dis;                                                          // dissipation rate ε* or specific dissipation ω*
    Array disn;                                                         // dis at time level n
    Array nue;                                                          // turbulent eddy viscosity ν_T*
    Array prod;                                                         // turbulence production P_k
    Array tke_source;                                                   // RHS source term for k equation
    Array dis_source;                                                   // RHS source term for ε/ω equation
    Array rhs_tke;                                                      // auxiliary RHS field for tke
    Array rhs_dis;                                                      // auxiliary RHS field for dis

    Array wbud_pgf;                                                     // w-budget: -dp/dphi /(rm*sinthe)  zonal pressure gradient (~0 in zonal mean)
    Array wbud_cor;                                                     // w-budget: Coriolis  -2*Omega*(costhe*v + sinthe*u)
    Array wbud_adv;                                                     // w-budget: advection -transport_w
    Array wbud_diff;                                                    // w-budget: diffusion (incl. curvature/metric)

    // Radial (u) momentum-budget diagnostic — attributes the DEEP polar
    // velocity blow-up that makes the spin-up not long-run stable (clean to
    // ~300-700, blows up by ~1500: N Pole -> 84N -> 51S). Captured into ubud_*
    // on the same wbudget_capture checkpoint iters; write_deep_momentum_budget
    // self-locates the max-|u| interior ocean cell and dumps the full split so
    // the driving term (metric-singular advection/diffusion, buoyancy, ...) is
    // isolated wherever the blow-up relocates. See project_hydro_polar_blowup.
    Array ubud_pgf;                                                     // u-budget: -dp/dr * exp_rm  radial pressure gradient
    Array ubud_adv;                                                     // u-budget: advection -transport_u
    Array ubud_diff;                                                    // u-budget: diffusion (eddy visc + 1/sin^2 metric — polar singular)
    Array ubud_buoy;                                                    // u-budget: Boussinesq buoyancy (thermal + haline anomaly)
    Array ubud_cor;                                                     // u-budget: Coriolis(radial) + centrifugal
    Array_2D vel_star;                                                  // friction velocity u_τ per (j,k) column [m/s]
};
#endif
