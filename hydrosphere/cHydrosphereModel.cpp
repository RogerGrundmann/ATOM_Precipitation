/*
 * Ocean General Circulation Modell(OGCM) applied to laminar flow
 * program for the computation of geo-atmospherical circulating flows in a spherical shell
 * finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 1 additional transport equations to describe the salinity
 * 4th order Runge-Kutta scheme to solve 2nd order differential equations inside an inner iterational loop
 * Poisson equation for the pressure solution in an outer iterational loop
 * temperature distribution given as a parabolic distribution from pole to pole, zonaly constant
 * code developed by Roger Grundmann, Zum Marktsteig 1, D-01728 Bannewitz(roger.grundmann@web.de)
*/

#include "cHydrosphereModel.h"
#include "BC_Hyd.h"
#include "PressureSolverHyd.h"
#include "ThermoHyd.h"
#include "ThermoHalineConveyorBelt.h"
#include "UtilsHyd.h"
#include "TurbulenceHyd.h"

using namespace std;
using namespace tinyxml2;
using namespace AtomUtils;

extern std::vector<std::vector<double> > m_node_weights;

cHydrosphereModel* cHydrosphereModel::m_model = NULL;

const double cHydrosphereModel::pi180 = 180.0/M_PI;                     // pi180 = 57.3

const double cHydrosphereModel::the_degree = 1.0;                       // compares to 1° step size laterally
const double cHydrosphereModel::phi_degree = 1.0;                       // compares to 1° step size longitudinally

const double cHydrosphereModel::dr = 0.025;                             // 0.025 x 40 = 1.0 compares to 200m : 40 = 5m for 1 radial step

const double cHydrosphereModel::dthe = the_degree/pi180;                // dthe = the_degree/pi180 = 1.0/57.3 = 0.01745, 180 * .01745 = 3.141
const double cHydrosphereModel::dphi = phi_degree/pi180;                // dphi = phi_degree/pi180 = 1.0/57.3 = 0.01745, 360 * .01745 = 6.282
    
const double cHydrosphereModel::the0 = 0.0;                             // North Pole
const double cHydrosphereModel::phi0 = 0.0;                             // zero meridian in Greenwich

//earth's radius is r_earth = 6731 km, here it is assumed to be infinity, circumference of the earth 40074 km 
const double cHydrosphereModel::r0 = 1.0;                               // non-dimensional
const double cHydrosphereModel::dr_stretch = 2.0;                       // sinh stretching for radial grid: cosh(2)≈3.76 step ratio surface/bottom
//const double cHydrosphereModel::r_Earth = 6731.001;                      // in km

//const double cHydrosphereModel::residuum_ref_hyd = 1.0e-4;            // criterium to finish global iterations

cHydrosphereModel::cHydrosphereModel():
    i_bathymetry(std::vector<std::vector<int> >(jm, std::vector<int>(km, 0))),

    has_printed_welcome_msg(false){

    // If Ctrl-C is pressed, quit
    signal(SIGINT, exit);

    // set default configuration
    SetDefaultConfig();

    m_model = this;

    //  Coordinate system in form of a spherical shell
    //  rad for r-direction bottom to the surface of the earth, the for lateral and phi for longitudinal direction
    rad.initArray_1D(im, 0);
    the.initArray_1D(jm, 0);
    phi.initArray_1D(km, 0);
}

cHydrosphereModel::~cHydrosphereModel(){
    m_model = NULL;
    logger().close();
}

#include "cHydrosphereDefaults.cpp.inc"

void cHydrosphereModel::LoadConfig(const char *filename){
    XMLDocument doc;
    XMLError err = doc.LoadFile(filename);
    try{
        if(err){
            doc.PrintError();
            throw std::invalid_argument(std::string("unable to load config file:  ") 
                + filename);
        }
        XMLElement *atom = doc.FirstChildElement("atom"), *elem_common = NULL, 
            *elem_hydrosphere = NULL;
        if(!atom){
            throw std::invalid_argument(std::string
                ("Failed to find the 'atom' element in config file: ") 
                + filename);
        }else{
            elem_common = atom->FirstChildElement("common");
            if(!elem_common){
                throw std::invalid_argument(std::string
                    ("Failed to find the 'common' element in 'atom' element in config file: ") 
                    + filename);
            }
            elem_hydrosphere = atom->FirstChildElement("hydrosphere");
            if(!elem_hydrosphere){
                throw std::invalid_argument(std::string(
                    "Failed to find the 'hydrosphere' element in 'atom' element in config file: ") 
                    + filename);
            }
        }
        #include "HydrosphereLoadConfig.cpp.inc"
    }catch(const std::exception &exc){
        std::cerr << exc.what() << std::endl;
        abort();
    }
}
/*
*
*/
void cHydrosphereModel::RunTimeSlice(int Ma){

    #ifdef _OPENMP
        printf("\n\n   number of processors: %d\n\n", omp_get_num_procs());

    #pragma omp parallel
        {
            printf("   thread %d of %d in \"Desktop-Dell-XPS 8960\"\n", 
                omp_get_thread_num(), omp_get_num_threads());
        }
    #else
        printf("   OpenMP is not supported\n");
    #endif
        printf("   ended\n\n");

    int start = AtomUtils::RunStart(hy);

    m_current_time = m_time_list.insert(int(Ma)).first;

    use_presets(use_earthbyte_reconstruction,
        use_NASA_temperature, use_NASA_velocity);

    // Apply depth mode.  "shallow" targets Ekman/near-surface flows (200 m, 5 m/step).
    // "deep" targets the thermohaline conveyor belt (6000 m, 150 m/step).
    // An empty or unrecognised mode leaves L_hyd at the value loaded from config.
    if (ocean_depth_mode == "shallow") {
        L_hyd = 200.0;
    } else if (ocean_depth_mode == "deep") {
        L_hyd = 6000.0;
    } else if (!ocean_depth_mode.empty()) {
        std::cerr << "WARNING: unknown ocean_depth_mode \"" << ocean_depth_mode
                  << "\" — using L_hyd = " << L_hyd << " m from config.\n";
    }

    UtilsHyd(*this).resetArrays();

//    dt = 0.0005;                                                        //  no dimension      prevents wiggles
    dt = dt_visc;                                                       //  no dimension; toggled to dt_inviscid during spin-up

    iter_n = 0;
    total_iter_count = 0;                                               // reset per time slice so the inviscid spin-up fires at the start of every Ma slice

    rad.StretchedCoordinates(im, r0, dr, dr_stretch);
    the.Coordinates(jm, the0, dthe);
    phi.Coordinates(km, phi0, dphi);

    read_Hydrosphere_Surface_Data(Ma);

//    if(use_NASA_velocity)
//        IC_vwt_WestEastCoast();                                         // horizontal velocity initial condition and temperature adjustments along coasts

    initTemperature(Ma);
    AtomUtils::damp_wiggles(t, &i_bathymetry, true, true, true);

//    if(!use_NASA_temperature)  IC_t_WestEastCoast();

    initSalinity();
    AtomUtils::damp_wiggles(c, &i_bathymetry, true, true, true);

    ThermoHyd(*this).SaltWaterDens();
//    UtilsHyd(*this).valueLimitationHyd();
    ThermoHyd(*this).SalinityEvaporation();
    ThermoHyd(*this).forces();
    LandOceanFraction();

    EkmanSpiral();                                                      // computes local velocity components by Ekmans spiral

    if (ocean_depth_mode == "deep") {
        ThermoHalineConveyorBelt(*this).run();                          // deep-water thermohaline circulation — only meaningful at full ocean depth
    }
    AtomUtils::damp_wiggles(u, &i_bathymetry, true, true, true);
    AtomUtils::damp_wiggles(v, &i_bathymetry, true, true, true);
    AtomUtils::damp_wiggles(w, &i_bathymetry, true, true, true);

    UtilsHyd(*this).valueLimitationHyd();

    BC_Hyd(*this).bcRadius();
    BC_Hyd(*this).bcTheta();
    BC_Hyd(*this).bcPhi();
    BC_Hyd(*this).bcSolidGround();
    BC_Hyd(*this).boundaryCondition();

    if (!turb_model.empty() && turb_model != "none") {
        TurbulenceHyd(*this).init();
        TurbulenceHyd(*this).apply_wall_bc();
    }

    UtilsHyd(*this).storeIntermediateData3D(1.0);


//    goto Printout;
    ThermoHyd(*this).runDataHyd();
    print_min_max_hyd();
    UtilsHyd(*this).writeFile(bathymetry_name, output_path, false);     // printing files for ParaView, AtmosphereDataTransfer and AtmospherePlotData
    if(Ma == 0) run_3D_loop();                                          // iterational 3D loop to solve variables in 4-step Runge-Kutta time scheme


    cout << endl << endl;
    if(Ma > 0) run_3D_loop();                                           // iterational 3D loop to solve variables in 4-step Runge-Kutta time scheme
    cout << endl << endl;

/*
    Printout:
    print_min_max_hyd();
    UtilsHyd(*this).writeFile(bathymetry_name, output_path, true);
*/

    RunEnd(hy, Ma, start);
    print_final_remarks();
    return;
/*
    t.printArray("OGCM", im, jm, km);
    v.printArray("OGCM", im, jm, km);
    w.printArray("OGCM", im, jm, km);
    temperature_NASA.printArray_2D(jm, km);
    temp_reconst.printArray_2D(jm, km);
    temp_landscape.printArray_2D(jm, km);
*/
}
/*
*
*/
void cHydrosphereModel::Run(){

    std::time_t Run_start;

    struct tm * timeinfo_begin;

    std::time(&Run_start);

    timeinfo_begin = std::localtime(&Run_start);

    std::cout << std::endl << std::endl;
    std::cout << " ... ATOM: time and date at run time begin:   " 
        << std::asctime(timeinfo_begin);
    std::cout << std::endl << std::endl;

    mkdir(output_path.c_str(), 0777);

    if(!has_printed_welcome_msg)  print_welcome_msg();

    for(int i = time_start; i <= time_end; i+=time_step){
        RunTimeSlice(i);
    }

    print_final_remarks();

    std::time_t Run_end;

    struct tm * timeinfo_end;

    std::time(&Run_end);

    timeinfo_end = std::localtime(&Run_end);

    std::cout << " ... ATOM: time and date at run time end:   " 
        << std::asctime(timeinfo_end) << std::endl;

    int Run_total = Run_end - Run_start;
    int Run_total_minutes = Run_total/60;
    int Run_total_hours = Run_total_minutes/60;

    std::cout <<  " ... ATOM: computer time needed for Ma = " << time_start 
        << " to Ma = " << time_end << " in " << time_step << " Ma steps:     "
        << Run_total << " seconds" << std::endl
        << " ... compares to:" << std::endl << std::endl
        << setw(30) << setfill(' ') << Run_total_hours 
        << " hours" << std::endl
        << setw(30) << setfill(' ') << Run_total_minutes 
        << " minutes" << std::endl
        << setw(30) << setfill(' ') << Run_total%60 
        << " seconds" << std::endl << std::endl
        << std::endl;

    return;
}
/*
*
*/
void cHydrosphereModel::run_3D_loop(){

cout << endl << endl << endl << "      OGCM: run_3D_loop ..........................." << endl;

    // Set turbulence model flags once, before the iteration loop
    if (!turb_model.empty() && turb_model != "none") {
        use_turbulence_model             = true;
        use_k_epsilon_turbulence_model   = (turb_model == "k_epsilon");
        use_k_omega_turbulence_model     = (turb_model == "k_omega");
        use_k_omega_SST_turbulence_model = (turb_model == "k_omega_SST");
    }

    // Resume from a binary checkpoint if requested (mirror of cAtmosphereModel).
    // load_state restores total_iter_count so the inviscid-phase schedule below
    // continues from where the saved run left off; a missing/mismatched file is a
    // no-op and the run proceeds from the IC.
    if(restart_from_iter >= 0)
        load_state(restart_from_iter);

    for(iter_n = 1; iter_n <= nm; iter_n++){

        print_loop_3D_headings();

        // Inviscid spin-up: zero diffusion for the first `inviscid_spinup_iters`
        // cumulative iterations, then ramp linearly to full over `inviscid_ramp_iters`.
        total_iter_count++;
        if(inviscid_spinup_iters <= 0){
            inviscid_phase = false;
            diffusion_ramp = 1.0;
        } else if(total_iter_count <= inviscid_spinup_iters){
            inviscid_phase = true;
            diffusion_ramp = 0.0;
        } else {
            inviscid_phase = false;
            int over = total_iter_count - inviscid_spinup_iters;
            diffusion_ramp = (inviscid_ramp_iters > 0)
                ? std::min(1.0, static_cast<double>(over) / inviscid_ramp_iters)
                : 1.0;
        }
        dt = inviscid_phase ? dt_inviscid : dt_visc;
        cout << "      OGCM: iter = " << total_iter_count
             << "  inviscid_phase = " << (inviscid_phase ? "true" : "false")
             << "  diffusion_ramp = " << std::fixed << std::setprecision(3) << diffusion_ramp
             << "  dt = " << std::scientific << std::setprecision(4) << dt
             << std::defaultfloat << endl;

        // w-momentum-budget diagnostic: on checkpoint iters, snapshot the
        // zonal-mean zonal velocity before the RK4 and have RHS_Hyd store its
        // per-term split this step (see write_w_momentum_budget).
        wbudget_capture = (iter_n % checkpoint == 0);
        if(wbudget_capture){
            wbar_before.assign(im, std::vector<double>(jm, 0.0));
            zonal_mean_w(wbar_before);
        }

        if (use_turbulence_model && !inviscid_phase) {
            solveRungeKutta_Hydrosphere_Turb();
        } else {
            solveRungeKutta_Hydrosphere();
        }

        wbudget_capture = false;   // wbud_* now hold this iter's term split

        // Heavy block (pressure solver, salinity/thermo, turbulence update) runs every 2
        // iterations in the viscous phase, but only every 10 iterations during the inviscid
        // spin-up — dt_inviscid is 10× smaller so the physical-time cadence stays comparable.
        int heavy_block_stride = inviscid_phase ? 10 : 2;
        if(iter_n % heavy_block_stride == 0){

            PressureSolverHyd(*this).run();
            AtomUtils::damp_wiggles(p_dyn, &i_bathymetry, true, true, true);

            ThermoHyd(*this).SaltWaterDens();
//            UtilsHyd(*this).valueLimitationHyd();
            ThermoHyd(*this).SalinityEvaporation();
            ThermoHyd(*this).runDataHyd();
            ThermoHyd(*this).forces();

            if (use_turbulence_model && !inviscid_phase) {
                TurbulenceHyd(*this).run();
                AtomUtils::damp_wiggles(tke, &i_bathymetry, true, true, true);
                AtomUtils::damp_wiggles(dis, &i_bathymetry, true, true, true);
                AtomUtils::damp_wiggles(nue, &i_bathymetry, true, true, true);
            }
        }  // iter_n % heavy_block_stride == 0

        UtilsHyd(*this).valueLimitationHyd();

        BC_Hyd(*this).bcRadius();
        BC_Hyd(*this).bcTheta();
        BC_Hyd(*this).bcPhi();
        BC_Hyd(*this).bcSolidGround();

        if (use_turbulence_model && !inviscid_phase) {
            TurbulenceHyd(*this).apply_wall_bc();
        }

        // Polar zonal (φ) filter — root-cause stabiliser for the high-latitude blow-up,
        // ported from the atmosphere (cAtmosphereModel.cpp). The ocean grid has the same
        // lat-lon polar singularity (zonal cell width rm·sinθ·dφ → 0 at the pole), but the
        // hydrosphere had no filter, so a long spin-up diverged at the North Pole (90°N
        // velocities → NaN ~iter 1500-1700, see project_hydro_polar_blowup). Damp the short
        // zonal wavelengths in u,v,w with a latitude-dependent pass count, every iteration
        // on the freshly-solved + BC'd field; filtered values flow into un,vn,wn via
        // storeIntermediateData3D below. Same params as the atmosphere (start=30°, gain=3).
        AtomUtils::polar_zonal_filter(u, the.z, &i_bathymetry, 30.0, 12, 3.0);
        AtomUtils::polar_zonal_filter(v, the.z, &i_bathymetry, 30.0, 12, 3.0);
        AtomUtils::polar_zonal_filter(w, the.z, &i_bathymetry, 30.0, 12, 3.0);

//        UtilsHyd(*this).findResiduumHyd();

        UtilsHyd(*this).storeIntermediateData3D(1.0);

//        print_min_max_hyd();

        if(iter_n % checkpoint == 0){
            print_min_max_hyd();

            UtilsHyd(*this).writeFile(bathymetry_name, output_path, true);
            cout << endl << "      OGCM: write_file in run_3D_loop ......................." << endl;

            // Binary restart checkpoint alongside the vtk/vts output, every
            // `checkpoint` iters, so a long spin-up can be resumed.
            save_state(total_iter_count);

            // Zonal-mean w-momentum budget CSV (uses the wbud_* term split
            // captured this iter + the wbar_before snapshot for the net Δwbar).
            write_w_momentum_budget(total_iter_count);
        }
    }  // end iter_n


    cout << endl << "      OGCM: run_3D_loop ended ..........................." << endl;
    return;
}
/*
*
*/

