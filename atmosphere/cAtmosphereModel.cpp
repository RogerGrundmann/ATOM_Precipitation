/*
 * Atmosphere General Circulation Modell(AGCM) applied to laminar flow
 * program for the computation of geo-atmospherical circulating flows in a spherical shell
 * including the magnetic field of the Earth
 * finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 1 additional transport equations to describe the salinity
 * 4th order Runge-Kutta scheme to solve 2nd order differential equations inside an inner iterational loop
 * Poisson equation for the pressure solution in an outer iterational loop
 * temperature distribution given as a parabolic distribution from pole to pole, zonaly constant
 * code developed by Roger Grundmann, Zum Marktsteig 1, D-01728 Bannewitz(roger.grundmann@web.de)
*/


#include "cAtmosphereModel.h"

#include "MoistConvection.h"

#include "ZeroCatIceScheme.h"
#include "OneCatIceScheme.h"
#include "TwoCatIceScheme.h"
#include "ThreeCatIceScheme.h"
#include "SaturationAdjustment.h"
#include "VelocityInitializer.h"
#include "PressureSolverAtm.h"
#include "ThermoAtm.h"
#include "UtilsAtm.h"
#include "BC_Atm.h"
#include "TurbulenceAtm.h"


using namespace std;
using namespace tinyxml2;
using namespace AtomUtils;


extern std::vector<std::vector<double> > m_node_weights;


cAtmosphereModel* cAtmosphereModel::m_model = NULL;

const double cAtmosphereModel::pi180 = 180.0/M_PI;                      // pi180 = 57.3

const double cAtmosphereModel::the_degree = 1.0;                        // compares to 1° step size laterally
const double cAtmosphereModel::phi_degree = 1.0;                        // compares to 1° step size longitudinally

const double cAtmosphereModel::dr = 0.025;                              // 0.025 x 40 = 1.0 compares to 16 km : 40 == 400 m for 1 radial step

const double cAtmosphereModel::dthe = the_degree/pi180; 
const double cAtmosphereModel::dphi = phi_degree/pi180;
    
const double cAtmosphereModel::the0 = 0.0;                              // North Pole
const double cAtmosphereModel::phi0 = 0.0;                              // zero meridian in Greenwich

//earth's radius is r_earth = 6731 km, here it is assumed to be infinity, circumference of the earth 40074 km 
const double cAtmosphereModel::r0 = 1.0;                                // non-dimensional
//const double cAtmosphereModel::r_Earth = 6731.001;                      // in km

const double cAtmosphereModel::residuum_ref_atm = 1.0e-2;               // criterium to finish global iterations




cAtmosphereModel::cAtmosphereModel():
    i_topography(std::vector<std::vector<int> >(jm, std::vector<int>(km, 0))),
    i_tropopause(std::vector<std::vector<int> >(jm, std::vector<int>(km, 0))),
    i_landscape(std::vector<std::vector<int> >(jm, std::vector<int>(km, 0))),

    has_printed_welcome_msg(false){

    // If Ctrl-C is pressed, quit
    signal(SIGINT, exit);

    // set default configuration
    SetDefaultConfig();

    m_model = this;

    //  Coordinate system in form of a spherical shell
    //  rad for r-direction normal to the surface of the earth, the for lateral and phi for longitudinal direction
    rad.initArray_1D(im, 0);                                            // radial coordinate direction
    the.initArray_1D(jm, 0);                                            // lateral coordinate direction
    phi.initArray_1D(km, 0);                                            // longitudinal coordinate direction
}

cAtmosphereModel::~cAtmosphereModel(){
    m_model = NULL;
    logger().close();
}
 
#include "cAtmosphereDefaults.cpp.inc"
/*
*
*/
void cAtmosphereModel::LoadConfig(const char *filename){
    XMLDocument doc;
    XMLError err = doc.LoadFile(filename);
    if(err){
        doc.PrintError();
        throw std::invalid_argument("   couldn't load config file inside cAtmosphereModel");
    }
    XMLElement *atom = doc.FirstChildElement("atom");
    if(!atom){
        return;
    }
    XMLElement* elem_common = doc.FirstChildElement("atom")->FirstChildElement("common");
    if(!elem_common){
        return;
    }
    XMLElement* elem_atmosphere = doc.FirstChildElement("atom")->FirstChildElement("atmosphere");
    if(!elem_atmosphere){
        return;
    }
#include "AtmosphereLoadConfig.cpp.inc"
}
/*
*
*/
void cAtmosphereModel::RunTimeSlice(int Ma){

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




    cout << endl << endl << endl 
        << " ----------------------------- Ma = "  << Ma 
        << " million years back in time ------------------------------ " 
        << endl << endl;

    int start = AtomUtils::RunStart(at);

    m_current_time = m_time_list.insert(int(Ma)).first;

    use_presets(use_earthbyte_reconstruction, 
        use_NASA_temperature, use_NASA_velocity);

    UtilsAtm(*this).resetArrays();

    dt = dt_visc;                                                       //  no dimension; toggled to dt_inviscid during spin-up

    iter_n = 0;
    total_iter_count = 0;                                               // reset per time slice so the inviscid spin-up fires at the start of every Ma slice

    use_turbulence_model             = (turb_model != "none");
    use_k_epsilon_turbulence_model   = (turb_model == "k_epsilon");
    use_k_omega_turbulence_model     = (turb_model == "k_omega");
    use_k_omega_SST_turbulence_model = (turb_model == "k_omega_SST");

    rad.Coordinates(im, r0, dr);
    the.Coordinates(jm, the0, dthe);
    phi.Coordinates(km, phi0, dphi);

    #pragma omp parallel sections
    {
        #pragma omp section
        { init_layer_heights(); }

        #pragma omp section
        { init_tropopause_layers(); }
    }

    read_Atmosphere_Surface_Data(Ma);                                   // reading topography data and NASA measurements

    LandOceanFraction();                                                // ratio of land to sea surface of the various time slices

    if(!use_NASA_velocity){
        VelocityInitializer(*this).compute();                           // construction of zonal initial velocities from measurements
    }

    // Skip the boundary-layer damping and the wiggle filter on the initial velocity
    // field when an inviscid spin-up is requested: those routines suppress exactly the
    // large-scale near-surface gradients that the Euler phase needs in order to organise
    // a global circulation. The viscous phase that follows will re-apply both.
    if(inviscid_spinup_iters <= 0){
        BC_Atm(*this).bcVelSurfSur();                                   // velocities close to surfaces, resembling a boundary layer
        AtomUtils::damp_wiggles(u, &i_topography, true, true, true);
        AtomUtils::damp_wiggles(v, &i_topography, true, true, true);
        AtomUtils::damp_wiggles(w, &i_topography, true, true, true);
    }

//    BC_Atm(*this).coastalCurrents();

    initTemperatureData(Ma);                                            // initialization of temperature, hydrostatic pressure and density of dry air, reconstruction of potential surface values
    AtomUtils::damp_wiggles(t, &i_topography, true, true, true);

    initWaterWapour();                                                  // initWaterWapour() and initCloudIce() belong together, init_vapour_cloud() stands alone
    initCloudIce();
//    init_vapour_cloud();                                                // initialisation of water vapour and cloud/ice formation based on the temperature profile

//    goto Printout;

    SaturationAdjustment(*this).run();                                  // based on the initial distribution, recomputation of the cloud water and cloud ice formation in case of saturated water vapour detected
    AtomUtils::damp_wiggles(ice,   &i_topography, true, true, true);
    AtomUtils::damp_wiggles(c,     &i_topography, true, true, true);
    AtomUtils::damp_wiggles(cloud, &i_topography, true, true, true);

    ThermoAtm(*this).precipitableWater();

    switch(CategoryIceScheme){                                          // rain, snow graupel and precipitation production and reduction
        case -1: cout << endl << endl << endl                           // no CategoryIceScheme used
            << "  no CategoryIceScheme used" << endl;
                break;
        case 0: ZeroCatIceScheme(*this).run();                          // development of rain fall, water vapour and cloud water
                break;
        case 1: OneCatIceScheme(*this).run();                           // development of rain and snow fall, water vapour, cloud water and ice
                break;
        case 2: TwoCatIceScheme(*this).run();                           // development of rain and snow fall, water vapour, cloud water and ice
                break;
        case 3: ThreeCatIceScheme(*this).run();                         // development of rain and snow fall, water vapour, cloud water, ice and graupel
                break;
    }
    AtomUtils::damp_wiggles(P_rain, &i_topography, true, true, true);
    AtomUtils::damp_wiggles(P_snow, &i_topography, true, true, true);

//    goto Printout;

    MoistConvection(*this).run(iter_n);                                 // rainfall from convecting clouds

    AtomUtils::damp_wiggles(P_conv, &i_topography, true, true, true);
    AtomUtils::damp_wiggles(E_u,    &i_topography, true, true, true);
    AtomUtils::damp_wiggles(M_u,    &i_topography, true, true, true);
    AtomUtils::damp_wiggles(q_v_u,  &i_topography, true, true, true);
    AtomUtils::damp_wiggles(q_c_u,  &i_topography, true, true, true);
    AtomUtils::damp_wiggles(E_d,    &i_topography, true, true, true);
    AtomUtils::damp_wiggles(M_d,    &i_topography, true, true, true);
    AtomUtils::damp_wiggles(q_v_d,  &i_topography, true, true, true);

//    goto Printout;

    if(turb_model != "none") {
        TurbulenceAtm(*this).init();                                    // turbulence field initialisation (skipped in laminar mode)
        TurbulenceAtm(*this).run();                                     // call for turbulence models

        AtomUtils::damp_wiggles(tke,        &i_topography, true, true, true);
        AtomUtils::damp_wiggles(dis,        &i_topography, true, true, true);
        AtomUtils::damp_wiggles(nue,        &i_topography, true, true, true);
        AtomUtils::damp_wiggles(prod,       &i_topography, true, true, true);
        AtomUtils::damp_wiggles(tke_source, &i_topography, true, true, true);
        AtomUtils::damp_wiggles(dis_source, &i_topography, true, true, true);

    }

    UtilsAtm(*this) .precipitationSum();
    ThermoAtm(*this).densities();
    ThermoAtm(*this).forces();
    ThermoAtm(*this).standAtm_DewPoint_HumidRel();                      // International Standard Atmosphere temperature profile, dew point temperature, relative humidity profile
    ThermoAtm(*this).waterVapourEvaporation();                          // correction of surface water vapour by evaporation
    ThermoAtm(*this).latentSensibleHeat();                              // latent and sensible heat
    ThermoAtm(*this).vegetationLand();                                  // vegetation on land
    ThermoAtm(*this).co2Atmosphere();                                   // greenhouse gas co2 as function of temperature
//    UtilsAtm(*this) .valueLimitationAtm();                              // value limitation prevents local formation of NANs

    BC_Atm(*this).bcRadius();                                           // extrapolation in i-direction alomg grid boundaries
    BC_Atm(*this).bcTheta();                                            // extrapolation in j-direction alomg grid boundaries
    BC_Atm(*this).bcPhi();                                              // extrapolation in k-direction alomg grid boundaries

    BC_Atm(*this).bcScalarSurfSur();                                    // scalar variable at surfaces extrapolated by von Neumann
    BC_Atm(*this).bcSolidGround();                                      // values inside mountains

    UtilsAtm(*this).storeIntermediateData3D(1.0);
    UtilsAtm(*this).findResiduumAtm();


//    goto Printout;
    print_min_max_atm();
    ThermoAtm(*this).printDataAtm();

//    Printout:
    UtilsAtm(*this).writeFile(bathymetry_name, output_path, false);     // printing files for ParaView, AtmosphereDataTransfer and AtmospherePlotData


    if(Ma == 0) run_3D_loop(Ma);                                        // iterational 3D loop to solve variables in 4-step Runge-Kutta time scheme

    cout << endl << endl;
    if(Ma > 0) run_3D_loop(Ma);                                         // iterational 3D loop to solve variables in 4-step Runge-Kutta time scheme
    cout << endl << endl;
/*
    print_min_max_atm();
    ThermoAtm(*this).printDataAtm();
    UtilsAtm(*this).writeFile(bathymetry_name, output_path, true);      // printing files for ParaView, AtmosphereDataTransfer and AtmospherePlotData
*/

    RunEnd(at, Ma, start);

    print_final_remarks();
}


/*
  Usage examples

  using namespace AtomUtils;

  // Horizontal wiggle damping (most common use case)
  AtomUtils::damp_wiggles(t, &i_topography, false, true, true);

  // Full 3D, two passes, reduced strength
  AtomUtils::damp_wiggles(cloud, &i_topography, true, true, true, 0.5, 2);

  // Vertical only, no topography mask
  AtomUtils::damp_wiggles(p_stat, nullptr, true, false, false);


    AtomUtils::damp_wiggles(u, &i_topography, true, true, true, 0.5, 2);
    AtomUtils::damp_wiggles(v, &i_topography, true, true, true, 0.5, 2);
    AtomUtils::damp_wiggles(w, &i_topography, true, true, true, 0.5, 2);

    AtomUtils::remove_peaks(BuoyancyForce, &i_bathymetry, true, true, true, 1.0, 2);


    temperature_NASA.printArray_2D(jm, km);
    temp_reconst.printArray_2D(jm, km);
    temp_landscape.printArray_2D(jm, km);

    t.printArray("AGCM", im, jm, km);
    c.printArray("AGCM", im, jm, km);
    cloud.printArray("AGCM", im, jm, km);
    P_rain.printArray("AGCM", im, jm, km);

    fft_gaussian_filter_3d(t,2);                                        // filter for smoothing temperature profiles, three grades of smoothing

    fft_gaussian_filter_3d(c,1);
    fft_gaussian_filter_3d(cloud,1);
    fft_gaussian_filter_3d(ice,1);

    fft_gaussian_filter_3d(P_rain,1);
    fft_gaussian_filter_3d(P_snow,1);
    fft_gaussian_filter_3d(P_graupel,1);
*/

/*
*
*/
void cAtmosphereModel::run_3D_loop(int Ma){ 

cout << endl << endl << endl << "      AGCM: run_3D_loop atm ..........................." << endl;

    for(iter_n = 1; iter_n <= nm; iter_n++){
        print_loop_3D_headings();

        // Inviscid spin-up: zero diffusion for the first `inviscid_spinup_iters`
        // cumulative iterations, then ramp linearly to full over `inviscid_ramp_iters`.
        // total_iter_count is per-RHS-call and persists across time slices.
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
        cout << "      AGCM: iter = " << total_iter_count
             << "  inviscid_phase = " << (inviscid_phase ? "true" : "false")
             << "  diffusion_ramp = " << std::fixed << std::setprecision(3) << diffusion_ramp
             << "  dt = " << std::scientific << std::setprecision(4) << dt
             << std::defaultfloat << endl;

        // Heavy block (pressure solver, saturation, ice scheme, moist convection, BCs, output)
        // runs every 2 iterations in the viscous phase, but only every 10 iterations during
        // the inviscid spin-up — dt_inviscid is 10× smaller so the physical-time cadence
        // stays comparable, and the spin-up doesn't spend most of its cycles in the heavy block.
        int heavy_block_stride = inviscid_phase ? 10 : 2;
        if(iter_n % heavy_block_stride == 0){
            PressureSolverAtm(*this).run();

            SaturationAdjustment(*this).run();                          // based on the initial distribution, recomputation of the cloud water and cloud ice formation in case of saturated water vapour detected
            AtomUtils::damp_wiggles(ice,   &i_topography, true, true, true);
            AtomUtils::damp_wiggles(c,     &i_topography, true, true, true);
            AtomUtils::damp_wiggles(cloud, &i_topography, true, true, true);

//            AtomUtils::remove_peaks(cloud, &i_topography, true, true, true);
//            AtomUtils::remove_peaks(ice,   &i_topography, true, true, true);

            switch(CategoryIceScheme){                                  // rain, snow graupel and precipitation production and reduction
                case -1: cout << endl << endl << endl                   // no CategoryIceScheme used
                    << "  no CategoryIceScheme used" << endl;
                    break;
                case 0: ZeroCatIceScheme(*this).run();                  // development of rain and snow fall, water vapour and cloud water
                    break;
                case 1: OneCatIceScheme(*this).run();                   // development of rain and snow fall, water vapour, cloud water and ice
                    break;
                case 2: TwoCatIceScheme(*this).run();                   // development of rain and snow fall, water vapour, cloud water and ice
                    break;
                case 3: ThreeCatIceScheme(*this).run();                 // development of rain and snow fall, water vapour, cloud water, ice and graupel
                    break;
            }

            AtomUtils::damp_wiggles(P_rain, &i_topography, true, true, true);
            AtomUtils::damp_wiggles(P_snow, &i_topography, true, true, true);

            MoistConvection(*this).run(iter_n);                         // rainfall from convecting clouds

            AtomUtils::damp_wiggles(P_conv, &i_topography, true, true, true);
            AtomUtils::damp_wiggles(E_u,    &i_topography, true, true, true);
            AtomUtils::damp_wiggles(M_u,    &i_topography, true, true, true);
            AtomUtils::damp_wiggles(q_v_u,  &i_topography, true, true, true);
            AtomUtils::damp_wiggles(q_c_u,  &i_topography, true, true, true);
            AtomUtils::damp_wiggles(E_d,    &i_topography, true, true, true);
            AtomUtils::damp_wiggles(M_d,    &i_topography, true, true, true);
            AtomUtils::damp_wiggles(q_v_d,  &i_topography, true, true, true);



            UtilsAtm(*this).precipitationSum();
            ThermoAtm(*this).densities();
            ThermoAtm(*this).forces();
            ThermoAtm(*this).standAtm_DewPoint_HumidRel();              // International Standard Atmosphere temperature profile, dew point temperature, relative humidity profile
            ThermoAtm(*this).waterVapourEvaporation();                  // correction of surface water vapour by evaporation
            ThermoAtm(*this).latentSensibleHeat();                      // latent and sensible heat
            ThermoAtm(*this).vegetationLand();                          // vegetation on land
            ThermoAtm(*this).co2Atmosphere();                           // greenhouse gas co2 as function of temperature

            if(turb_model != "none") {
                TurbulenceAtm(*this).run();        // update turbulence sources (skipped in laminar mode)

            AtomUtils::damp_wiggles(tke,        &i_topography, true, true, true);
            AtomUtils::damp_wiggles(dis,        &i_topography, true, true, true);
            AtomUtils::damp_wiggles(nue,        &i_topography, true, true, true);
            AtomUtils::damp_wiggles(prod,       &i_topography, true, true, true);
            AtomUtils::damp_wiggles(tke_source, &i_topography, true, true, true);
            AtomUtils::damp_wiggles(dis_source, &i_topography, true, true, true);
            }

//            UtilsAtm(*this).valueLimitationAtm();                       // value limitation prevents local formation of NANs

            BC_Atm(*this).bcRadius();                                   // extrapolation in i-direction alomg grid boundaries
            if(turb_model != "none") TurbulenceAtm(*this).apply_wall_bc();  // reassert ω_wall at i=0 after bcRadius cubic extrapolation
            BC_Atm(*this).bcTheta();                                    // extrapolation in j-direction alomg grid boundaries
            BC_Atm(*this).bcPhi();                                      // extrapolation in k-direction alomg grid boundaries

            BC_Atm(*this).bcScalarSurfSur();                            // scalar variable at surfaces extrapolated by von Neumann
            BC_Atm(*this).bcSolidGround();                              // values inside mountains

            print_min_max_atm();
            UtilsAtm(*this).writeFile(bathymetry_name, output_path, false);
            cout << endl << "      AGCM: write_file in run_3D_loop atm ......................." << endl;

        }  // iter_n % heavy_block_stride == 0

        if(turb_model == "none" || inviscid_phase) solveRungeKutta_Atmosphere();   // laminar: standard RK4 on transport equations (also during inviscid spin-up)
        else                                       solveRungeKutta_Atmosphere_Turb();  // turbulent: RK4 extended with k and ω equations

        UtilsAtm(*this).findResiduumAtm();

        UtilsAtm(*this).storeIntermediateData3D(1.0);

        ThermoAtm(*this).printDataAtm();

    }  // end iter_n


    cout << endl << "      AGCM: run_3D_loop atm ended ..........................." << endl;
}
/*
*
*/
void cAtmosphereModel::Run(){

    std::time_t Run_start;

    struct tm * timeinfo_begin;

    std::time(&Run_start);

    timeinfo_begin = std::localtime(&Run_start);

    std::cout << std::endl << std::endl;
    std::cout << " ... AGCM: time and date at run time begin:   " 
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

    std::cout << " ... AGCM: time and date at run time end:   " 
        << std::asctime(timeinfo_end) << std::endl;

    int Run_total = Run_end - Run_start;
    int Run_total_minutes = Run_total/60;
    int Run_total_hours = Run_total_minutes/60;

    std::cout <<  " ... AGCM: computer time needed for Ma = " << time_start 
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
}
/*
*
*/
