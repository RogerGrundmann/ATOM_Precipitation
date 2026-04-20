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


#include "cCubeModel.h"

#include "PressureSolverCube.h"
#include "UtilsCube.h"
#include "BC_Cube.h"
#include "TurbulenceCube.h"


using namespace std;
using namespace AtomUtils;


extern std::vector<std::vector<double> > m_node_weights;


cCubeModel* cCubeModel::m_model = NULL;

const double cCubeModel::dr = 0.025;                                    // uniform z step (vertical)
const double cCubeModel::dy = dr;                                       // uniform y step (lateral) — equal to dr so cube face is square
const double cCubeModel::dx = 2.0 * dr;                                 // uniform x step (streamwise) — 0.05: doubles downstream extent to 30H, reduces outlet reflection

const double cCubeModel::y0 = -((jm-1) / 2.0) * dy;                     // domain centred at y=0
const double cCubeModel::x0 = 0.0;

const double cCubeModel::r0 = 1.0;                                      // non-dimensional

const double cCubeModel::residuum_ref_atm = 1.0e-2;                     // criterium to finish global iterations




cCubeModel::cCubeModel():
    i_topography(std::vector<std::vector<int> >(jm, std::vector<int>(km, 0))),

    has_printed_welcome_msg(false){


    // set default configuration
    SetDefaultConfig();

    m_model = this;

    rad.initArray_1D(im, 0);                                            // radial coordinate direction
    the.initArray_1D(jm, 0);                                            // lateral coordinate direction
    phi.initArray_1D(km, 0);                                            // longitudinal coordinate direction
}

cCubeModel::~cCubeModel(){
    m_model = NULL;
}

void cCubeModel::LoadConfig(const char* /*filename*/)
{
    // All parameters are set by SetDefaultConfig() in the constructor.
    // XML loading is not yet implemented for the cube model.
}
 

void cCubeModel::RunTimeSlice(int Ma){

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


    cout << endl << endl << endl;

    int start = AtomUtils::RunStart(at);

    UtilsCube(*this).resetArrays();

    dt = 0.01;                                                          //  no dimension      prevents wiggles in laminar flows

    iter_n = 0;

    use_turbulence_model             = (turb_model != "none");
    use_k_epsilon_turbulence_model   = (turb_model == "k_epsilon");
    use_k_omega_turbulence_model     = (turb_model == "k_omega");
    use_k_omega_SST_turbulence_model = (turb_model == "k_omega_SST");

    rad.Coordinates(im, r0, dr);                                         // z: uniform
    the.Coordinates(jm, y0, dy);                                        // y: uniform
    phi.Coordinates(km, x0, dx);                                      
      // x: uniform

    #pragma omp parallel sections
    {
        #pragma omp section
        { init_layer_heights(); }
    }

    UtilsCube(*this).initCubeContour();                                 // initialise cube contour (h=1 inside cube, 0 elsewhere)
    UtilsCube(*this).initVelocities();                                  // initialise cube contour (h=1 inside cube, 0 elsewhere)

    BC_Cube(*this).bcVelSurfSur();                                      // velocities close to surfaces, resembling a boundary layer

    AtomUtils::damp_wiggles(u, &i_topography, true, true, true);
    AtomUtils::damp_wiggles(v, &i_topography, true, true, true);
    AtomUtils::damp_wiggles(w, &i_topography, true, true, true);

//    goto Printout;

    if(turb_model != "none") {
        TurbulenceCube(*this).init();                                   // turbulence field initialisation (skipped in laminar mode)
        TurbulenceCube(*this).run();                                    // call for turbulence models

        AtomUtils::damp_wiggles(tke,        &i_topography, true, true, true);
        AtomUtils::damp_wiggles(dis,        &i_topography, true, true, true);
        AtomUtils::damp_wiggles(nue,        &i_topography, true, true, true);
        AtomUtils::damp_wiggles(prod,       &i_topography, true, true, true);
        AtomUtils::damp_wiggles(tke_source, &i_topography, true, true, true);
        AtomUtils::damp_wiggles(dis_source, &i_topography, true, true, true);

    }

//    UtilsAtm(*this) .valueLimitationAtm();                              // value limitation prevents local formation of NANs

    BC_Cube(*this).bcRadius();                                          // extrapolation in i-direction alomg grid boundaries
    BC_Cube(*this).bcTheta();                                           // extrapolation in j-direction alomg grid boundaries
    BC_Cube(*this).bcPhi();                                             // extrapolation in k-direction alomg grid boundaries

    BC_Cube(*this).bcScalarSurfSur();                                   // scalar variable at surfaces extrapolated by von Neumann
    BC_Cube(*this).bcSolidGround();                                     // values inside mountains

    UtilsCube(*this).storeIntermediateData3D(1.0);
//    UtilsAtm(*this).findResiduumAtm();


//    goto Printout;
    print_min_max_atm();
//    UtilsCube(*this).printDataAtm();

//    Printout:
    UtilsCube(*this).writeFile(bathymetry_name, output_path, false);    // printing files for ParaView, AtmosphereDataTransfer and AtmospherePlotData


    cout << endl << endl;
    if(iter_n == 0) run_3D_loop();                                      // iterational 3D loop to solve variables in 4-step Runge-Kutta time scheme
    cout << endl << endl;


    RunEnd(at, start);
    print_final_remarks();
    return;
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
void cCubeModel::run_3D_loop(){ 

//cout << endl << endl << endl << "      CUBE: run_3D_loop cube ..........................." << endl;
 
    panorama_cnt = 0;

    for(iter_n = 1; iter_n <= nm; iter_n++){
//        print_loop_3D_headings();

        if(iter_n % 2 == 0){
            PressureSolverCube(*this).run();

            if(turb_model != "none") {
                TurbulenceCube(*this).run();        // update turbulence sources (skipped in laminar mode)
/*
            AtomUtils::damp_wiggles(u,          &i_topography, true, true, true);
            AtomUtils::damp_wiggles(v,          &i_topography, true, true, true);
            AtomUtils::damp_wiggles(w,          &i_topography, true, true, true);
*/
            AtomUtils::damp_wiggles(tke,        &i_topography, true, true, true);
            AtomUtils::damp_wiggles(dis,        &i_topography, true, true, true);
            AtomUtils::damp_wiggles(nue,        &i_topography, true, true, true);
            AtomUtils::damp_wiggles(prod,       &i_topography, true, true, true);
            AtomUtils::damp_wiggles(tke_source, &i_topography, true, true, true);
            AtomUtils::damp_wiggles(dis_source, &i_topography, true, true, true);

            }

//            UtilsCube(*this).valueLimitationAtm();                      // value limitation prevents local formation of NANs

            BC_Cube(*this).bcRadius();                                  // extrapolation in i-direction alomg grid boundaries
            if(turb_model != "none") TurbulenceCube(*this).apply_wall_bc();  // reassert ω_wall at i=0 after bcRadius cubic extrapolation
            BC_Cube(*this).bcTheta();                                   // extrapolation in j-direction alomg grid boundaries
            BC_Cube(*this).bcPhi();                                     // extrapolation in k-direction alomg grid boundaries

            BC_Cube(*this).bcScalarSurfSur();                           // scalar variable at surfaces extrapolated by von Neumann
            BC_Cube(*this).bcSolidGround();                             // values inside mountains

        }  // iter_n % 2 == 0


        solveRungeKutta_Cube_Turb();                                    // turbulent: RK4 extended with k and ω equations

        AtomUtils::damp_wiggles(u,          &i_topography, true, true, true);
        AtomUtils::damp_wiggles(v,          &i_topography, true, true, true);
        AtomUtils::damp_wiggles(w,          &i_topography, true, true, true);

//        UtilsCube(*this).findResiduumAtm();

        BC_Cube(*this).bcRadiationOutflow();                            // Sommerfeld radiation: advect k=km-1 out at local flow speed (u,v,w)
        BC_Cube(*this).bcSpongeInflow();                                // Rayleigh sponge: relax u,v,w to background near k=0

        UtilsCube(*this).storeIntermediateData3D(1.0);

        panorama_cnt++;


        if(iter_n % checkpoint == 0){
            print_min_max_atm();

            UtilsCube(*this).writeFile(bathymetry_name, output_path, true);
//            cout << endl << "      Cube: write_file in run_3D_loop ......................." << endl;
        }

        if((iter_n % 2 == 0)&&(panorama_cnt == panorama_print)){
            panorama_cnt = 0;
        }


/*
        if(panorama_cnt == panorama_print){
            print_min_max_atm();
            UtilsCube(*this).writeFile(bathymetry_name, output_path, false);
            cout << endl << "      CUBE: write_file in run_3D_loop cube ......................." << endl;
            panorama_cnt = 0;
        }
*/
    }  // end iter_n


//    cout << endl << "      CUBE: run_3D_loop cube ended ..........................." << endl;
}
/*
*
*/
void cCubeModel::Run(){

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

    RunTimeSlice(0);

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

    std::cout <<  " ... AGCM: computer time needed:     "
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
