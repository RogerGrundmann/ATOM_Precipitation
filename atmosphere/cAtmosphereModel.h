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
    friend class MultiLayerRadiation;
    friend class RadiationSelfTest;   // offline standalone driver (test/ dir) — no run-loop call site

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
        double curv;            // J'/J, the curvature term in d2f/dz2; 0 on the legacy metric                                                                                                                                                                              
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

    // THE definition of the grid coordinates. rad/the/phi are filled here and nowhere else;
    // every other site that used to write them (UtilsAtm::resetArrays, MoistConvection::precompute,
    // test/rad_selftest) calls this instead, so there is one place to change and no way for two
    // definitions to drift apart. See checkMetricConsistency() for why that matters.
    void initGridCoordinates();

    // Startup check that the radius the CORE metric uses and the radius the PHYSICS uses are the
    // same length. See the definition for the numbers and the reasoning.
    void checkMetricConsistency() const;

    // Startup check that exp_rm is the Jacobian of the radial stretch -- it is not. Print-only
    // and self-silencing once the spread drops below 1.05. Call AFTER init_layer_heights.
    void checkRadialMetric() const;

    // ATM_METRIC_RADIUS: reads the knob once and fills m_metric_r0. Call after the coordinates.
    void initMetricRadius();

    // The radius the HORIZONTAL metric should use, in rad.z units, given a grid coordinate rm.
    // Identity when the knob is off, so every call site is bit-identical by default.
    double metricRadius(double rm) const {
        return (m_metric_r0 > 0.0) ? (m_metric_r0 + (rm - rad.z[0])) : rm;
    }

    // Planetary radius in rad.z units when ATM_METRIC_RADIUS is on; 0.0 means off.
    double m_metric_r0 = 0.0;

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

    // panorama_cnt IS DEAD AND WAS UNINITIALISED UNTIL 2026-08-30. UtilsAtm.h records that the
    // counter it named was replaced by an `iter_n % panorama_print` test, and the member was left
    // behind: no assignment exists anywhere in the tree, yet PrintMsg prints it on every
    // iteration line. `cAtmosphereModel model;` is a STACK LOCAL in main, so the printed value
    // was stack garbage and reading it was undefined behaviour -- it was found because deleting
    // Q_Sensible shifted the object layout and the printed value changed from 0 to 4 with no
    // other difference in the run. Initialised here rather than removed, because dropping it
    // from the print changes the format of every iteration line.
    int panorama_cnt = 0, iter_n;

    // Inviscid spin-up state.
    // total_iter_count accumulates across time slices so the inviscid window is global, not per-slice.
    // diffusion_ramp ∈ [0,1] multiplies every diff_*_re coefficient in RHS_Atm and the no-slip flag.
    int total_iter_count = 0;
    double diffusion_ramp = 1.0;
    bool inviscid_phase = false;
    bool ubudget_capture = false;   // when true, rhs_u stores its per-term split into ubud_* (set on checkpoint iters)
    bool vbudget_capture = false;   // when true, rhs_v stores its per-term split into vbud_* (set on checkpoint iters)
    bool wbudget_capture = false;   // when true, rhs_w stores its per-term split into wbud_* (set on checkpoint iters)

    // Buoyancy ramp ∈ [0,1] — linearly increases the Boussinesq body force in rhs_u
    // from 0 at iter 0 to 1 at iter buoyancy_ramp_iters.  Set ramp_iters = 0 to
    // disable (buoyancy_ramp stays 1.0).  The earlier "500-iter ramp starved the
    // system" note applied to the OLD ~336× too-weak buoyancy; with the 2026-06-19
    // g/(omega*L_atm) scaling fix the body force is now ~336× stronger, so switching
    // it on cold CFL-blows up. Ramp it in over a few hundred iters so the pressure
    // field and circulation adjust gradually to the corrected forcing.
    double buoyancy_ramp = 1.0;
    static constexpr int buoyancy_ramp_iters = 300;

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

    // Anelastic base state, ported from ATHAD 2026-08-27. m_rho_base is the cos-latitude
    // weighted horizontal mean of r_humid at each level; m_dlnrho_dr is d ln(rho_bar)/d(rad.z),
    // i.e. differentiated in the GRID radial coordinate, so it composes with exp_rm the same way
    // every other radial derivative in the solver does. One-dimensional: the horizontal density
    // contrast is small against the vertical span, so a 1-D base state loses almost nothing and
    // keeps the elliptic operator constant in time.
    std::vector<double> m_rho_base;
    std::vector<double> m_dlnrho_dr;

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
    std::vector<double> m_layer_J;   // dz/d(rad.z) per level; only filled for the pressure grid

    // Per-level horizontal-mean (non-dim) temperature, used as the Boussinesq
    // buoyancy base state so the body force has zero mean at every height and
    // does not fight the hydrostatic p_stat/p_dyn split. Refilled once per
    // RK4 step by computeLevelMeanTemperature().
    std::vector<double> t_ref_level;
    // Virtual-temperature counterpart of t_ref_level, for ATM_BUOY_MOIST. It must exist and be
    // built the same way: the buoyancy's design property is that the body force has ZERO MEAN at
    // every height, so moistening the parcel without moistening the reference would add a
    // uniform upward force everywhere instead of a buoyancy.
    std::vector<double> tv_ref_level;
    void computeLevelMeanTemperature();

    // Initial (non-dim) temperature at the model lid (i=im-1), snapshotted once
    // from the IC. bcRadius pins t at the lid to this so the isothermal-floor top
    // stays constant instead of drifting upward through the old cubic top
    // extrapolation. Sized [jm][km]; empty until RunStart populates it.
    std::vector<std::vector<double>> t_top_init;

    void SetDefaultConfig();
    void print_min_max_atm();
    void write_meridional_streamfunction(int iter);   // zonal-mean v + meridional mass streamfunction Ψ (Hadley/Ferrel cell diagnostic)
    // Zonal-mean meridional-wind momentum budget (Hadley/Ferrel spin-down attribution):
    // zonal_mean_v fills vbar[i][j] in m/s; write_v_momentum_budget writes the per-step
    // Δv̄ contributions (RK4 dynamics + each post-RK4 filter) differenced across one iteration.
    void zonal_mean_v(std::vector<std::vector<double> >& vbar);
    void write_v_momentum_budget(int iter,
        const std::vector<std::vector<double> >& dv_dyn,
        const std::vector<std::vector<double> >& dv_polar,
        const std::vector<std::vector<double> >& dv_orog,
        const std::vector<std::vector<double> >& dv_radial);
    // Zonal-mean ZONAL-wind (w) momentum budget — the trade / Walker component that the
    // meridional (v) budget cannot see. zonal_mean_w fills wbar[i][j] in m/s; the writer
    // attributes the per-iteration Δw̄ to RK4 dynamics + each post-RK4 filter, so the
    // trade-easterly spin-down can be pinned to a source (weakening Coriolis) or a sink
    // (a specific filter / diffusion). Same masking + scaling as the v budget.
    void zonal_mean_w(std::vector<std::vector<double> >& wbar);
    void write_w_momentum_budget(int iter,
        const std::vector<std::vector<double> >& dw_dyn,
        const std::vector<std::vector<double> >& dw_polar,
        const std::vector<std::vector<double> >& dw_orog,
        const std::vector<std::vector<double> >& dw_radial);
    void run_3D_loop(int Ma);
    void load_global_temperature_curve();
    void load_equat_temperature_curve();
    void load_pole_temperature_curve();
    void calculate_node_weights();
    void init_steps();
    void init_tropopause_layers();
    void RHS_Atmosphere_Turb(int i, int j, int k, const CellGeometry& geo);   // single dynamical core (laminar RHS_Atmosphere dropped 2026-07-08)
    void solveRungeKutta_Atmosphere_Turb();
    void fft(Array &);
    void LandOceanFraction();
    void initTemperatureData(int Ma);

    void initWaterWapour();
    void initCloudIce();
    void init_vapour_cloud();
    void cloudiness_backup();
    void init_Maxwell();

    // ==================================================================
    // ATM_GRID_PRESSURE -- place the levels uniformly in ln p on a reference hydrostatic
    // column instead of exponentially in height. Ported from ATHAD (README items 82/83),
    // DEFAULT OFF, and the two constants below are FITTED FOR THIS TREE: they are not
    // ATHAD's and must not be copied between trees.
    //
    //   ATM_GRID_PTOP  p_top/p_0, default 0.08538 -- the pressure at THIS tree's own legacy
    //                  lid (16023 m, 86.5 hPa). ATHAD ships 1e-6, which on Earth's column
    //                  sits at 81.9 km: importing that default would not regrid this model,
    //                  it would EXTEND the shell 5.1x and quietly make it a different model.
    //   ATM_GRID_BETA  ln-p stretch, default 3.988 -- fitted so the bottom layer matches the
    //                  legacy 38.9 m. At beta = zeta = 3.715 the bottom is 1.23x coarser.
    //
    // EXPECT NO PAYOFF HERE, and the reason is structural rather than a matter of tuning.
    // The ladder has only Lambda = -ln(p_top/p_0) = 2.46 e-foldings to redistribute across
    // this shell, against ATHAD's 13.8, because shell/H = 2.01 here (ATHAD 5.06, ATHAD_COND
    // 7.74, ATHAD_PERID 15.54, and PERID is the one tree where the branch paid). For an
    // ISOTHERMAL column a ln-p ladder IS a height ladder; the two differ only through the
    // temperature contrast up the column, and Earth's (288 -> 214 K) is the mildest in the
    // family. Measured beta table for this tree, dz_0 relative to legacy:
    //
    //     beta      1.0    2.0    3.0    3.715   4.5
    //     dz_0     7.42x  4.05x  2.06x   1.23x  0.68x       (matched at beta = 3.988)
    //
    // The lid is invariant under beta -- only ATM_GRID_PTOP moves it.
    // ==================================================================
    static bool gridPressure(){
        static const bool v = [](){
            const char* e = getenv("ATM_GRID_PRESSURE"); return e && atoi(e) != 0; }();
        return v;
    }
    static double gridPTop(){
        static const double v = [](){
            const char* e = getenv("ATM_GRID_PTOP");
            const double d = e ? atof(e) : 0.08538;
            return (d > 0.0 && d < 1.0) ? d : 0.08538; }();
        return v;
    }
    static double gridBeta(){
        static const double v = [](){
            const char* e = getenv("ATM_GRID_BETA");
            const double d = e ? atof(e) : 3.988;
            return (d > 0.0) ? d : 3.988; }();
        return v;
    }

    // A dry hydrostatic reference column: ln(p/p_0) against height, on THIS tree's own
    // initialisation law -- the 6.5 K/km reference lapse of InitValues_Atm.cpp, isothermal
    // above the mean tropopause. ATHAD integrates a DRY ADIABAT here because its column is
    // one by construction; using that law on Earth would put the 16 km lid at 117 K and the
    // ladder would be built for an atmosphere this model does not have.
    bool buildReferenceColumn(std::vector<double>& z_ref,
                              std::vector<double>& lnp_ref) const {
        const double T_s    = 0.5 * ((t_equat_modern + t_0) + (t_pole_modern + t_0));
        const double lapse  = 0.0065;                                  // LAPSE_RATE_REF [K/m]
        const double z_trop = 0.5 * (tropopause_equator + tropopause_pole);
        const double T_trop = T_s - lapse * z_trop;
        if(!(T_s > 0.0) || !(R_Air > 0.0) || !(p_0 > 0.0) || !(g > 0.0) || !(T_trop > 0.0))
            return false;

        const double dz      = 25.0;          // m, fine enough that the ladder is smooth
        const double z_limit = 2.0e5;         // m, a guard and nothing more
        const double lnp_target = std::log(gridPTop());

        z_ref.clear(); lnp_ref.clear();
        z_ref.push_back(0.0); lnp_ref.push_back(0.0);        // ln(p/p_0) = 0 at the ground

        double T = T_s, lnp = 0.0, z = 0.0;
        while(z < z_limit && lnp > lnp_target){
            double T_next = T_s - lapse * (z + dz);
            if(T_next < T_trop) T_next = T_trop;             // isothermal stratosphere
            const double T_mid = 0.5 * (T + T_next);
            if(!(T_mid > 0.0)) break;
            lnp -= g * dz / (R_Air * T_mid);                 // d(ln p) = -g dz /(R T)
            z   += dz;
            T    = T_next;
            z_ref.push_back(z); lnp_ref.push_back(lnp);
        }
        return (lnp <= lnp_target) && (z_ref.size() > 2);
    }

    // dz/d(rad.z) per level, central-differenced from the height table. A pressure-placed
    // grid has no closed-form Jacobian; the legacy branch keeps its analytic one, so that
    // branch stays bit-identical.
    void buildMetricTable(){
        m_layer_J.assign(im, 0.0);
        const double dr_loc = (im > 1) ? (rad.z[im-1] - rad.z[0]) / (double)(im - 1) : 1.0;
        for(int i = 0; i < im; i++){
            const int lo = (i == 0) ? 0 : i - 1;
            const int hi = (i == im - 1) ? im - 1 : i + 1;
            const double dz = (double)m_layer_heights[hi] - (double)m_layer_heights[lo];
            const double dn = (double)(hi - lo) * dr_loc;
            m_layer_J[i] = (dn > 0.0) ? (dz / dn) : 1.0;
        }
    }

    void init_layer_heights(){
        m_layer_heights.clear();
        m_layer_J.clear();

        if(gridPressure()){
            std::vector<double> z_ref, lnp_ref;
            if(buildReferenceColumn(z_ref, lnp_ref)){
                const double Lambda = -std::log(gridPTop());       // e-foldings, > 0
                const double beta   = gridBeta();
                const double denom  = std::exp(beta) - 1.0;
                std::size_t k = 0;
                for(int i = 0; i < im; i++){
                    const double t = (double)i / (double)(im - 1);
                    const double x = (denom > 0.0)
                                   ? (std::exp(beta * t) - 1.0) / denom : t;
                    const double want = -x * Lambda;               // target ln(p/p_0), <= 0
                    while(k + 1 < lnp_ref.size() && lnp_ref[k + 1] > want) k++;
                    double zi = z_ref.back();
                    if(k + 1 < lnp_ref.size()){
                        const double d = lnp_ref[k] - lnp_ref[k + 1];
                        const double f = (d > 0.0) ? (lnp_ref[k] - want) / d : 0.0;
                        zi = z_ref[k] + f * (z_ref[k + 1] - z_ref[k]);
                    }
                    m_layer_heights.push_back((float)zi);
                }
                // Strictly increasing, or every radial derivative divides by zero.
                bool ok = true;
                for(int i = 1; i < im; i++)
                    if(!(m_layer_heights[i] > m_layer_heights[i-1])) ok = false;
                if(ok){ buildMetricTable(); return; }
                std::cout << "      AGCM: ATM_GRID_PRESSURE produced a non-monotonic grid"
                          << " - falling back to the legacy stretch" << std::endl;
                m_layer_heights.clear();
            } else {
                std::cout << "      AGCM: ATM_GRID_PRESSURE could not reach p_top = "
                          << gridPTop() << " p_0 - falling back to the legacy stretch"
                          << std::endl;
            }
        }

        float h = L_atm;
        m_layer_heights.clear();
        for(int i=0; i<im; i++){
            // rad.z[0] instead of a hardcoded 1.0: the surface is wherever the radial coordinate
            // starts, which ATM_METRIC_RADIUS may move.
            m_layer_heights.push_back((exp(zeta
            * (rad.z[i] - rad.z[0])) - 1) * h);                         // in m      local atmospheric shell thickness
//            std::cout << m_layer_heights.back() << std::endl;
        }
        return;
    }

    // ==================================================================
    // IS exp_rm THE JACOBIAN OF THE RADIAL STRETCH? It is documented as one in
    // TurbulenceAtm.h and written as one in PressureSolverAtm.h, and it is not. Ported from
    // ATHAD README item 80; checkRadialMetric() has printed the spread at every startup since
    // the instruments went in, and in THIS tree it is 23.21x -- the worst in the family.
    //
    // exp_rm = 1/(rm+1) is the Jacobian of a QUADRATIC stretch z ~ (rm+1)^2/2, while
    // init_layer_heights() builds an EXPONENTIAL one, z = (exp(zeta*(r-r0)) - 1)*L_atm. The
    // true Jacobian is J(r) = dz/d(rad.z) = zeta*L_atm*exp(zeta*(r-r0)) [m per rad.z unit],
    // and the core wants it as a DIMENSIONLESS factor against its own length unit, which is
    // metricShellLength() -- the metres one rad.z unit represents on average.
    //
    // AND THE SECOND DERIVATIVE NEEDS A TERM THE CODE DOES NOT HAVE. With e = U/J,
    //     U^2 * d2f/dz2 = e^2 * ( d2f/dr2 - (J'/J) * df/dr )
    // and J'/J is `zeta` for an exponential stretch. The core computes d2f/dr2 * exp_2_rm and
    // stops, so the curvature term is missing in BOTH metrics -- small under the legacy one
    // (J'/J = 1/(rm+1) <= 0.5) and the same order as the retained term under the true one,
    // since zeta = 3.715 here. metricCurv() returns 0 on the legacy branch so that branch
    // stays bit-identical; the missing legacy term is a separate, smaller defect, recorded
    // rather than silently fixed.
    //
    // ATM_METRIC_EXACT=1 switches to the true Jacobian. DEFAULT OFF, and in ATHAD it is
    // default off for a MEASURED reason rather than caution: making the metric exact there
    // left the OLR untouched but made the pressure projection WORSE (div(rho u)/rho rms
    // 2.739e-02 -> 7.722e-02, Psi_max jumping to the ground). Whatever this tree does with
    // it, it must be measured here before any flip.
    // ==================================================================
    static bool metricExact(){
        static const bool v = [](){
            const char* e = getenv("ATM_METRIC_EXACT"); return e && atoi(e) != 0; }();
        return v;
    }

    // dz/d(rad.z) at rm, in metres per rad.z unit. Analytic: this tree has one grid.
    int metricLevelOf(double rm) const {
        const double dr_loc = (im > 1) ? (rad.z[im-1] - rad.z[0]) / (double)(im - 1) : 1.0;
        int i = (dr_loc > 0.0) ? (int)std::lround((rm - rad.z[0]) / dr_loc) : 0;
        if(i < 0) i = 0;
        if(i > im - 1) i = im - 1;
        return i;
    }
    double metricJ(double rm) const {
        if(gridPressure() && (int)m_layer_J.size() == im) return m_layer_J[metricLevelOf(rm)];
        return zeta * L_atm * exp(zeta * (rm - rad.z[0]));
    }

    // The dimensionless radial Jacobian factor a FIRST derivative is multiplied by.
    double metricExpRm(double rm) const {
        if(!metricExact()) return 1.0 / (rm + 1.0);
        const double J = metricJ(rm);
        return (J > 0.0) ? (metricShellLength() / J) : (1.0 / (rm + 1.0));
    }

    // ATM_METRIC_NOCURV -- exact Jacobian, curvature term OFF. An ATTRIBUTION knob, not a
    // physics option: ATM_METRIC_EXACT changes two things at once (the first-derivative factor
    // and the second-derivative curvature term), and measuring it showed the vertical wind
    // halving with no integrated diagnostic registering it. This isolates which half does that.
    // Inert unless ATM_METRIC_EXACT is also set, since metricCurv is 0 on the legacy branch
    // anyway.
    static bool metricNoCurv(){
        static const bool v = [](){
            const char* e = getenv("ATM_METRIC_NOCURV"); return e && atoi(e) != 0; }();
        return v;
    }

    // J'/J, the coefficient of the curvature term in d2f/dz2 = e^2*(f'' - curv*f').
    // Zero on the legacy branch, so the legacy operator is unchanged to the bit.
    double metricCurv(double rm) const {
        if(!metricExact() || metricNoCurv()) return 0.0;
        if(gridPressure() && (int)m_layer_J.size() == im){
            // J'/J by central difference on the same table, in rad.z units.
            const int i = metricLevelOf(rm);
            const int lo = (i == 0) ? 0 : i - 1, hi = (i == im - 1) ? im - 1 : i + 1;
            const double dr_loc = (im > 1) ? (rad.z[im-1] - rad.z[0]) / (double)(im - 1) : 1.0;
            const double dn = (double)(hi - lo) * dr_loc;
            const double Ji = m_layer_J[i];
            return (dn > 0.0 && Ji > 0.0) ? ((m_layer_J[hi] - m_layer_J[lo]) / dn) / Ji : 0.0;
        }
        return zeta;
    }

    // Physical length that ONE unit of rad.z represents, in metres.
    //
    // This is the number the metric needs and the one that is easy to get wrong. rad.z runs
    // 1.0 .. 2.0 and init_layer_heights above maps that span onto 0 .. 16.02 km, so one rad.z
    // unit is the SHELL THICKNESS, ~16 km — NOT L_atm. L_atm = 400 m is the amplitude of the
    // exponential stretch, not a grid step: the actual layer spacing runs from 39 m at the
    // surface to 1457 m at the top, and the "400 m for 1 radial step" in the dr comment
    // (cAtmosphereModel.cpp) is nominal only.
    //
    // The two readings differ by exactly 1/dr = 40, which is why the planetary radius expressed
    // in rad.z units is 6370/16.02 = 397.5 and not 6.37e6/400 = 15925.
    // ==================================================================
    // ATM_LENGTH_NDIM -- one length for the whole non-dimensionalisation. DEFAULT OFF.
    //
    // L_atm = 400 m is the AMPLITUDE OF THE EXPONENTIAL STRETCH, not a grid step: one rad.z
    // unit is metricShellLength() ~ 16023 m. force_nd was moved onto the latter on 2026-07-28
    // (the "40" in 16000 = 397.5 x 40, see the note in RHS_Atm_Turb.cpp), and that note lists
    // the terms left behind: the Held-Suarez relaxation, the Rayleigh surface drag, coeff_MC_*,
    // coeff_S, coeff_L and nue_max -- plus nue_air_nd, which has the same shape and was not
    // named. This knob moves all of them together, because a nondimensionalisation is only
    // meaningful as a whole.
    //
    // THE DIRECTION IS NOT UNIFORM, AND "40x too weak" IS ONLY HALF THE STORY. Where the length
    // is in the NUMERATOR -- rates, k*L/u_0: Held-Suarez, surf_drag, coeff_MC_* -- the shipped
    // terms are 40x too WEAK and this makes them stronger. Where it is in the DENOMINATOR --
    // diffusivities, nu/(u_0*L): nue_max, nue_air_nd, coeff_S, coeff_L -- the shipped terms are
    // 40x too STRONG and this makes them weaker. Both follow from the same substitution.
    //
    // NOT FLIPPED ON BY DEFAULT, and not because of caution: a 40x change to the surface drag
    // and the eddy-viscosity ceiling is a change to the shipped dynamics, and this tree's rule
    // is that such a thing is measured behind a knob first. Off-branch bit-identical by
    // construction -- ndimLength() returns L_atm unchanged.
    // ==================================================================
    static bool lengthNdimConsistent(){
        static const bool v = [](){
            const char* e = getenv("ATM_LENGTH_NDIM"); return e && atoi(e) != 0; }();
        return v;
    }
    double ndimLength() const {
        return lengthNdimConsistent() ? metricShellLength() : L_atm;
    }

    double metricShellLength() const {
        const double span = rad.z[im-1] - rad.z[0];
        if(!(span > 0.0)) return L_atm;
        return (exp(zeta * span) - 1.0) * L_atm / span;
    }

    void init_topography(const string &topo_filename);
    void save_data();
    void save_array(const string& fn, const Array& a);
    std::vector<Array*> restart_arrays();   // the prognostic 3D fields a checkpoint serializes
    void save_state(int iter, int Ma);      // dump restart_arrays() + total_iter_count to a binary file (name carries Ma+iter)
    bool load_state(int iter, int Ma);      // restore them; returns false (and runs from scratch) if absent/mismatched
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
    void read_Hydrosphere_SST(int Ma);   // reverse coupling: blend hydrosphere SST into t.x[0] (Picard loop)
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
    Array_2D albedo;                                                    // surface albedo (pole->equator parabola) — MultiLayerRadiation
    Array_2D epsilon_2D;                                                // surface emissivity (bottom layer of epsilon) — MultiLayerRadiation
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
    Array t_eq;                                                         // Held-Suarez radiative-equilibrium temperature target (snapshot of the initial field)
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
    Array BuoyancyForce;                                                // buoyancy force, Boussinesque approximation
    Array CoriolisForce;                                                // Coriolis force terms
    Array CentrifugalForce;                                             // centrifugal force terms
    Array PresGradForce;                                                // Force caused by normal pressure gradient
    // Brunt-Vaisala frequency squared, N^2 = (g/theta) d(theta)/dz. Ported from ATHAD, and
    // it is the FIRST diagnostic in this tree with an externally known right answer: Earth's
    // troposphere sits near 1e-4 s^-2 and the stratosphere near 4e-4, so this can be checked
    // against reality rather than against another arm of the same model.
    //
    // THE VERTICAL DERIVATIVE USES TRUE HEIGHTS, get_layer_height(), NOT the core's exp_rm.
    // That is deliberate: it makes the field a PHYSICAL quantity, comparable with the numbers
    // above, and therefore usable as a check ON the metric instead of a victim of it. See
    // checkRadialMetric() for the 23.2x spread that would otherwise contaminate it.
    // Long-wave optical depth, ported from ATHAD (README item 42, where the photosphere had
    // been quoted for months from an offline Python column because the model computed it
    // NOWHERE). Both are already formed inside MultiLayerRadiation and were being thrown away
    // after conversion to epsilon -- and they cannot be recovered from epsilon afterwards,
    // because epsilon saturates at tau ~ 37.
    //
    //   tau_layer  per-layer d(tau). THE RESOLUTION MEASURE: if one layer carries d(tau) of
    //              order 1 or more, the two-stream sweep (first order in d(tau)) is resolving
    //              the emission level with a single cell.
    //   tau_above  cumulative optical depth of everything ABOVE a level: 0 at the lid,
    //              increasing downward. tau_above = 1 is the photosphere, i.e. where the
    //              atmosphere actually radiates to space.
    Array tau_above;
    Array tau_layer;
    Array brunt_N2;
    // Psi -- the meridional mass streamfunction as a FIELD [kg/s]. Ported from ATHAD 2026-08-27.
    // write_meridional_streamfunction() already formed psi[i][j] and threw it away except for a
    // CSV and a printed Psi_max; publishing it puts the cells in ParaView and in the Results
    // min/max, instead of the whole cell-structure argument resting on one scalar. That
    // mattered in ATHAD: Psi_max was reporting the SPURIOUS SURFACE FLUX rather than the
    // circulation (README item 68), which a scalar cannot show and a field can. It is a zonal
    // mean, so it is replicated across k.
    Array Psi;
    // RADIAL momentum-budget term capture (diagnostic): per-cell rhs_u contributions,
    // stored when ubudget_capture is set. Ported from ATHAD (README item 42): the vertical
    // component was the one with NO instrument, so a spurious radial acceleration is
    // invisible until something decomposes rhs_u. Two things the split makes checkable at a
    // glance -- with coriolis_nontraditional() false, ubud_cor must be identically zero, and
    // ubud_buoy is the size of the buoyancy body force, which is the prior question behind
    // ATM_BUOY_TREF / ATM_BUOY_CONSISTENT. Diagnostic only: no physics reads these.
    Array ubud_pgf;                                                     // -∂p/∂r ·exp_rm (radial pressure gradient)
    Array ubud_cor;                                                     // Coriolis (non-traditional; off by default)
    Array ubud_advv;                                                    // vertical advection  -u·∂u/∂r
    Array ubud_advh;                                                    // horizontal advection
    Array ubud_diff;                                                    // diffusion (molecular + turbulent, + metric)
    Array ubud_buoy;                                                    // buoyancy body force
    // Zonal-mean v momentum-budget term capture (diagnostic): per-cell rhs_v
    // contributions, stored when vbudget_capture is set so write_v_momentum_budget
    // can attribute the Hadley/Ferrel spin-down to a specific dynamical term.
    Array vbud_pgf;                                                     // -∂p/∂θ /rm (meridional pressure gradient)
    Array vbud_cor;                                                     // Coriolis term (2·cosθ·w coupling to zonal wind)
    Array vbud_advv;                                                    // vertical advection  -u·∂v/∂r
    Array vbud_advh;                                                    // horizontal advection -(v/rm)∂v/∂θ -(w/rmsinθ)∂v/∂φ
    Array vbud_diff;                                                    // diffusion (molecular + turbulent, + metric terms)
    Array vbud_other;                                                   // surface drag + moist-convection momentum
    // Zonal-mean w (zonal-wind / trade) momentum-budget term capture — mirror of vbud_*,
    // stored when wbudget_capture is set so write_w_momentum_budget can attribute the
    // trade-easterly spin-down to a specific rhs_w term.
    Array wbud_pgf;                                                     // -∂p/∂φ /(rm·sinθ) (zonal pressure gradient)
    Array wbud_cor;                                                     // Coriolis term (coupling to radial/meridional wind)
    Array wbud_advv;                                                    // vertical advection  -u·∂w/∂r
    Array wbud_advh;                                                    // horizontal advection -(v/rm)∂w/∂θ -(w/rmsinθ)∂w/∂φ
    Array wbud_diff;                                                    // diffusion (molecular + turbulent, + metric terms)
    Array wbud_other;                                                   // surface drag + moist-convection momentum
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
