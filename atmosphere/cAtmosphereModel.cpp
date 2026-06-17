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
//    if(inviscid_spinup_iters <= 0){
        BC_Atm(*this).bcVelSurfSur();                                   // velocities close to surfaces, resembling a boundary layer

        AtomUtils::damp_wiggles(u, &i_topography, true, true, true);
        AtomUtils::damp_wiggles(v, &i_topography, true, true, true);
        AtomUtils::damp_wiggles(w, &i_topography, true, true, true);
//    }

//    BC_Atm(*this).coastalCurrents();

    initTemperatureData(Ma);                                            // initialization of temperature, hydrostatic pressure and density of dry air, reconstruction of potential surface values
    AtomUtils::damp_wiggles(t, &i_topography, true, true, true);

    // Snapshot the lid temperature (i=im-1) from the IC. bcRadius pins t at the lid
    // to this fixed reference so the isothermal stratospheric top stays constant
    // instead of drifting up via the former cubic top extrapolation.
    t_top_init.assign(jm, std::vector<double>(km, 0.0));
    for(int j = 0; j < jm; j++)
        for(int k = 0; k < km; k++)
            t_top_init[j][k] = t.x[im-1][j][k];

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

    // One-shot Helmholtz projection of the prescribed initial velocity.  The
    // analytical Hadley/Ferrel profile in VelocityInitializer is NOT divergence-
    // free; without this step the incremental projection in PressureSolverAtm
    // strips the dilatational component over the first ~100 iters and the
    // prescribed global circulation decays to noise (max u → <0.5 m/s by iter 150).
    // Doing the projection once here lets the time loop start from a clean
    // ∇·v = 0 state with the solenoidal part of the prescribed flow intact.
    PressureSolverAtm(*this).project_initial_velocity(200);

    UtilsAtm(*this).storeIntermediateData3D(1.0);
    UtilsAtm(*this).findResiduumAtm();


//    goto Printout;


    print_min_max_atm();
    ThermoAtm(*this).printDataAtm();

    UtilsAtm(*this).writeFile(bathymetry_name, output_path, false);     // printing files for ParaView, AtmosphereDataTransfer and AtmospherePlotData


    if(Ma == 0) run_3D_loop(Ma);                                        // iterational 3D loop to solve variables in 4-step Runge-Kutta time scheme

    cout << endl << endl;
    if(Ma > 0) run_3D_loop(Ma);                                         // iterational 3D loop to solve variables in 4-step Runge-Kutta time scheme
    cout << endl << endl;

/*
    Printout:
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

    // Restart shortcut: overwrite the just-built initial conditions with a saved 3D
    // state and resume from its total_iter_count (skips re-running the dry spin-up).
    const bool restarted = (restart_from_iter >= 0) && load_state(restart_from_iter);

    // On a successful restart, start the loop counter at the restart point so the
    // printed iter_n matches where the run actually resumed (not back at 1).
    const int iter_start = restarted ? restart_from_iter : 1;
    for(iter_n = iter_start; iter_n <= nm; iter_n++){
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

        // Buoyancy ramp: linearly from 0 at iter 0 to 1 at buoyancy_ramp_iters.
        // Spreads the build-up of the thermal driving force across many iters so the
        // multi-sweep Jacobi can converge p_dyn against the growing forcing instead
        // of falling behind and ringing.
        buoyancy_ramp = (buoyancy_ramp_iters > 0)
            ? std::min(1.0, static_cast<double>(total_iter_count) / buoyancy_ramp_iters)
            : 1.0;

        cout << "      AGCM: iter = " << total_iter_count
             << "  inviscid_phase = " << (inviscid_phase ? "true" : "false")
             << "  diffusion_ramp = " << std::fixed << std::setprecision(3) << diffusion_ramp
             << "  buoyancy_ramp = " << std::fixed << std::setprecision(3) << buoyancy_ramp
             << "  dt = " << std::scientific << std::setprecision(4) << dt
             << std::defaultfloat << endl;

        // Two cadences during inviscid spin-up:
        //   moist_stride    — PressureSolverAtm, SaturationAdjustment, ice scheme,
        //                     MoistConvection, ThermoAtm, turbulence sources. Throttled
        //                     because (a) the moist/thermo sources feed back unbounded
        //                     in the undamped (diffusion_ramp=0) flow, and (b) the
        //                     pressure solver is a single Jacobi sweep whose target
        //                     balances the full body-force field — calling it every step
        //                     races p_dyn toward that (very large) balance value.
        //   momentum_stride — all BCs (radial/theta/phi extrapolation, scalar surface,
        //                     solid ground). Run every step during spin-up so the
        //                     topographic walls are re-enforced before each RK4 cycle —
        //                     this is what lets the flow "see" the mountains and form
        //                     lee/windward circulations.
        // In the viscous phase moist_stride collapses back to 2 (original behaviour).
        // momentum_stride is kept at 1 in BOTH phases: with stride 2 the solid-ground
        // wall BC was re-enforced only every other RK4 cycle, leaving a 2Δt computational
        // mode pinned at steep high-latitude coasts (Antarctic, ×2.5 via 1/sinθ at 68°S)
        // that grew to a near-surface NaN seed by ~iter 142. Re-enforcing the walls every
        // step (the existing, tested BCs) suppresses it; cost is BCs running 2× as often.
        int moist_stride    = inviscid_phase ? 10  : 2;
        int momentum_stride = 1;

        if(iter_n % moist_stride == 0){
            // Multi-sweep Jacobi: the original single-sweep cadence (one sweep every
            // moist_stride iters) couldn't keep up with the divergence accumulated by
            // a fully-spun-up global circulation, letting p_dyn drift to ±750 hPa by
            // iter 120 and the simulation crash into the velocity caps by iter 300.
            // 20 sweeps per call (~10× the previous rate) actually converges the
            // Poisson before the next time step, so dpdr/dpdthe/dpdphi in the RHS
            // reflect the real pressure each iter.  Verbose only on the first sweep.
            {
                constexpr int n_pressure_sweeps = 20;
                PressureSolverAtm solver(*this);
                for(int s = 0; s < n_pressure_sweeps; s++){
                    solver.run(s == 0);
                }
            }

            // Moist-physics gate: the implicit microphysics (SaturationAdjustment,
            // ice scheme, MoistConvection) is too stiff during the initial velocity
            // transient — a single tropical maritime column repeatedly drove q_c, q_i
            // and S_s into runaway. Delay activation until the circulation has formed
            // on a dry field; the moist physics then sees a stable advection pattern
            // it was designed for. 0 disables the gate (always on).
            const bool moist_phys_active = (moist_phys_start_iter <= 0)
                                        || (total_iter_count > moist_phys_start_iter);

            if(moist_phys_active){
                SaturationAdjustment(*this).run();                      // based on the initial distribution, recomputation of the cloud water and cloud ice formation in case of saturated water vapour detected
                AtomUtils::damp_wiggles(ice,   &i_topography, true, true, true);
                AtomUtils::damp_wiggles(c,     &i_topography, true, true, true);
                AtomUtils::damp_wiggles(cloud, &i_topography, true, true, true);

                switch(CategoryIceScheme){                              // rain, snow graupel and precipitation production and reduction
                    case -1: cout << endl << endl << endl               // no CategoryIceScheme used
                        << "  no CategoryIceScheme used" << endl;
                        break;
                    case 0: ZeroCatIceScheme(*this).run();              // development of rain and snow fall, water vapour and cloud water
                        break;
                    case 1: OneCatIceScheme(*this).run();               // development of rain and snow fall, water vapour, cloud water and ice
                        break;
                    case 2: TwoCatIceScheme(*this).run();               // development of rain and snow fall, water vapour, cloud water and ice
                        break;
                    case 3: ThreeCatIceScheme(*this).run();             // development of rain and snow fall, water vapour, cloud water, ice and graupel
                        break;
                }

                AtomUtils::damp_wiggles(P_rain, &i_topography, true, true, true);
                AtomUtils::damp_wiggles(P_snow, &i_topography, true, true, true);

                MoistConvection(*this).run(iter_n);                     // rainfall from convecting clouds

                AtomUtils::damp_wiggles(P_conv, &i_topography, true, true, true);
                AtomUtils::damp_wiggles(E_u,    &i_topography, true, true, true);
                AtomUtils::damp_wiggles(M_u,    &i_topography, true, true, true);
                AtomUtils::damp_wiggles(q_v_u,  &i_topography, true, true, true);
                AtomUtils::damp_wiggles(q_c_u,  &i_topography, true, true, true);
                AtomUtils::damp_wiggles(E_d,    &i_topography, true, true, true);
                AtomUtils::damp_wiggles(M_d,    &i_topography, true, true, true);
                AtomUtils::damp_wiggles(q_v_d,  &i_topography, true, true, true);
                // The convective tendencies fed straight into rhs_t/rhs_v/rhs_w
                // (coeff_MC_*·MC in RHS_Atm) enter spatially rough over steep orography;
                // smooth them like the other MC fields so the momentum/heat source does
                // not inject a jagged divergence the pressure solver amplifies (Andes seed).
                AtomUtils::damp_wiggles(MC_t,   &i_topography, true, true, true);
                AtomUtils::damp_wiggles(MC_v,   &i_topography, true, true, true);
                AtomUtils::damp_wiggles(MC_w,   &i_topography, true, true, true);

                UtilsAtm(*this).precipitationSum();

                // Physical caps on the microphysics source terms.
                // The ice-scheme S-terms feed rhs_t (latent heat: S_c,S_r,S_i,S_s,S_g)
                // and rhs_c (moisture: S_v) UNCAPPED. They depend on c/cloud/ice, so
                // c,cloud,ice ↔ S_* ↔ rhs_t/rhs_c is a closed loop. Over the steep Andes
                // (≈50°S, 1/sin²θ-amplified orographic ascent) the condensation↔latent-
                // heat↔ascent feedback drives them to ~1e260 by iter ~309 — the leg of
                // the moist runaway the MoistConvection caps (rhsForcing) don't cover.
                // Capping S_* bounds rhs_{c,cloud,ice,t}, hence c/cloud/ice, hence S_*
                // next call. Healthy |S| ≲ 1e-4 (kg/kg)/s (e.g. S_c_c, S_s diagnostics),
                // so 1e-3 is ~10× headroom and clips only the runaway. NaN/Inf→0 via the
                // bit-level is_finite_safe (std::min/max are unreliable under -ffast-math).
                {
                    constexpr double S_max = 1.0e-3;                     // [(kg/kg)/s]
                    auto cap_S = [&](Array& S){
                        #pragma omp parallel for collapse(2) schedule(static)
                        for(int i = 0; i < im; i++){
                            for(int j = 0; j < jm; j++){
                                for(int k = 0; k < km; k++){
                                    double v = S.x[i][j][k];
                                    if(!AtomUtils::is_finite_safe(v)) S.x[i][j][k] = 0.0;
                                    else if(v >  S_max)               S.x[i][j][k] =  S_max;
                                    else if(v < -S_max)               S.x[i][j][k] = -S_max;
                                }
                            }
                        }
                    };
                    cap_S(S_v); cap_S(S_c); cap_S(S_i); cap_S(S_r);
                    cap_S(S_s); cap_S(S_g); cap_S(S_c_c);
                }
            }  // moist_phys_active

            ThermoAtm(*this).densities();
            ThermoAtm(*this).forces();
            ThermoAtm(*this).standAtm_DewPoint_HumidRel();              // International Standard Atmosphere temperature profile, dew point temperature, relative humidity profile
            ThermoAtm(*this).waterVapourEvaporation();                  // correction of surface water vapour by evaporation
            ThermoAtm(*this).latentSensibleHeat();                      // latent and sensible heat
            ThermoAtm(*this).vegetationLand();                          // vegetation on land
            ThermoAtm(*this).co2Atmosphere();                           // greenhouse gas co2 as function of temperature

            if(turb_model != "none") {
                TurbulenceAtm(*this).run();                             // update turbulence sources (skipped in laminar mode)

                AtomUtils::damp_wiggles(tke,        &i_topography, true, true, true);
                AtomUtils::damp_wiggles(dis,        &i_topography, true, true, true);
                AtomUtils::damp_wiggles(nue,        &i_topography, true, true, true);
                AtomUtils::damp_wiggles(prod,       &i_topography, true, true, true);
                AtomUtils::damp_wiggles(tke_source, &i_topography, true, true, true);
                AtomUtils::damp_wiggles(dis_source, &i_topography, true, true, true);
            }

//            UtilsAtm(*this).valueLimitationAtm();                       // value limitation prevents local formation of NANs
        }  // iter_n % moist_stride == 0

        if(iter_n % momentum_stride == 0){
            BC_Atm(*this).bcRadius();                                   // extrapolation in i-direction alomg grid boundaries
            if(turb_model != "none") TurbulenceAtm(*this).apply_wall_bc();  // reassert ω_wall at i=0 after bcRadius cubic extrapolation
            BC_Atm(*this).bcTheta();                                    // extrapolation in j-direction alomg grid boundaries
            BC_Atm(*this).bcPhi();                                      // extrapolation in k-direction alomg grid boundaries

            BC_Atm(*this).bcScalarSurfSur();                            // scalar variable at surfaces extrapolated by von Neumann
            BC_Atm(*this).bcSolidGround();                              // values inside mountains

            if(iter_n % checkpoint == 0){
                print_min_max_atm();
                write_meridional_streamfunction(iter_n);   // Hadley/Ferrel cell strength (zonal-mean v + Ψ) per vtk checkpoint
                UtilsAtm(*this).writeFile(bathymetry_name, output_path, false);
                cout << endl << "      AGCM: write_file in run_3D_loop atm ......................." << endl;
            }
        }  // iter_n % momentum_stride == 0

        if(turb_model == "none" || inviscid_phase) solveRungeKutta_Atmosphere();   // laminar: standard RK4 on transport equations (also during inviscid spin-up)
        else                                       solveRungeKutta_Atmosphere_Turb();  // turbulent: RK4 extended with k and ω equations

        // Polar zonal (φ) filter — root-cause stabiliser for the high-latitude blow-up.
        // Near the poles the zonal cell width rm·sinθ·dφ → 0, so the explicit zonal CFL
        // diverges and short φ-waves in u,v,w explode within a single RK4 step (the clamp
        // in bcSolidGround only bounds RK4 inputs, never the in-step growth). This damps
        // those short zonal wavelengths with a latitude-dependent pass count; runs every
        // iteration on the freshly-solved field, and the filtered values flow into
        // un,vn,wn via storeIntermediateData3D below.
        // start=30°, gain=3 strengthens the high-mid-latitude band (≈45–67°) where the
        // dry near-surface w-mode grows (seed 52°N, explosion 60°N/150°W Gulf of Alaska);
        // gives ~3 passes at 52°, ~6 at 60°, while ≤35° (subtropical jet) stays at 0.
        AtomUtils::polar_zonal_filter(u, the.z, &i_topography, 30.0, 12, 3.0);
        AtomUtils::polar_zonal_filter(v, the.z, &i_topography, 30.0, 12, 3.0);
        AtomUtils::polar_zonal_filter(w, the.z, &i_topography, 30.0, 12, 3.0);

        // Orographic Shapiro filter — same purpose as the polar filter, but for the
        // steep-orography CFL hot-spots (Andes/Atlas/NZ Alps) the latitude-keyed polar
        // filter doesn't reach. Damps the near-surface 2Δ velocity oscillation at cliffs
        // that the bounded pressure projection alone could not contain. passes=3 (was 1)
        // restores near-surface damping after n_layers_above was trimmed 12→3 to keep
        // Tibet aloft; heavier configs (n_layers_above=5, passes=6) tested 2026-05-31
        // SPREAD the Gulf-of-Alaska coastal blow-up across more cells (each pass averages
        // the ±100 cap into neighbours) — wrong tool, the runaway is forcing-driven, not
        // a 2Δ checkerboard. Use this gentle setting and attack the MC forcing instead.
        AtomUtils::orographic_shapiro_filter(u, i_topography, 3, 3, 3);
        AtomUtils::orographic_shapiro_filter(v, i_topography, 3, 3, 3);
        AtomUtils::orographic_shapiro_filter(w, i_topography, 3, 3, 3);

        // Radial (i) de-checkerboarding. The 2Δ grid-scale checkerboard that seeds the
        // near-surface CFL blow-up has its purest component on the radial axis, and no other
        // filter (polar=k, orographic=j+k steep-only) reaches it. The vertical velocity u gets
        // a plain 1-2-1 pass — it carries little resolved shear and this annihilates the i=1
        // high-latitude 2Δ spike. For v and w a 4th-order Shapiro is used instead: it still
        // removes the 2Δ mode completely (CFL guard) but preserves the resolved vertical shear.
        // A 1-2-1 on v,w run EVERY iteration measurably erodes the jet's thermal-wind shear and
        // the Hadley/Ferrel cells' two-branch v structure, spinning the circulation down despite
        // a steady baroclinic (equator-pole) temperature gradient; the 4th-order filter does not.
        AtomUtils::radial_shapiro_filter   (u, i_topography, /*passes=*/2);
        AtomUtils::radial_shapiro_filter_ho(v, i_topography, /*passes=*/2);
        AtomUtils::radial_shapiro_filter_ho(w, i_topography, /*passes=*/2);

        // (Upper-troposphere Rayleigh sponge tested 2026-06-01: did not help. The iter-359
        // "first NaN at i=35, k=0, 10 km" reported by scanForNaN was the overflow location,
        // not the cause — by iter 341 the Gulf-of-Alaska surface cell (60°N/151°W/298m)
        // already had t=-14254°C, p_dyn=-10132 hPa, density=-117 kg/m³, with the coastal
        // sponge pegging u/v/w at ±30 m/s every step. High-i cells overflow first because
        // they're produced by RK4 from the corrupted surface dynamics. The fix belongs at
        // the coastal cell, not aloft — see [[project-coastal-sponge-session]] continuation.)
        // AtomUtils::upper_rayleigh_sponge(u, 25, 0.5);
        // AtomUtils::upper_rayleigh_sponge(v, 25, 0.5);
        // AtomUtils::upper_rayleigh_sponge(w, 25, 0.5);

        // Near-surface coastal velocity HARD CAP — last-resort ceiling for the dry-seeded
        // velocity runaway at steep coasts (Gulf-of-Alaska 60°N/151°W). The slope cap
        // skips ocean↔land cliffs and bcVelSurfSur is init-only; without an in-loop RHS
        // friction term, coastal pressure-gradient force drives a linear velocity growth
        // that saturates at the ±100 m/s BC cap. At the bottom `layers` air cells of a
        // coastal/sloped column, if |u_nondim| > 3.0 (=30 m/s), clamp to ±3.0. Rate=1.0
        // selects the hard-clamp branch (Rayleigh form rounds 6-7 was too weak — per-step
        // RHS forcing >> %-damping). layers=8: the runaway is at i≈5 (298 m) above the
        // OCEAN side of the Cook Inlet coast, not at the land cliff itself — layers=3
        // (round 8) only reached i=0..2 over ocean columns and missed the actual blow-up.
        // 8 covers i_topo..i_topo+7 = 0..7 (≈440 m) above ocean / 298 m..1.5 km over land.
        AtomUtils::coastal_velocity_sponge(u, i_topography, 3.0, 1.0, 8);
        AtomUtils::coastal_velocity_sponge(v, i_topography, 3.0, 1.0, 8);
        AtomUtils::coastal_velocity_sponge(w, i_topography, 3.0, 1.0, 8);

        scanForNaN();                                                       // debug: report field/location/iter of the FIRST non-finite cell to localise the blow-up

        UtilsAtm(*this).findResiduumAtm();

        UtilsAtm(*this).storeIntermediateData3D(1.0);

        ThermoAtm(*this).printDataAtm();

        // Checkpoint the full 3D state once total_iter_count reaches the requested iter.
        if(checkpoint_save_iter >= 0 && total_iter_count == checkpoint_save_iter)
            save_state(checkpoint_save_iter);

        // Periodic restart checkpoints every 100 iters, but ONLY once the run is clean
        // (no non-finite cell in any serialized prognostic field) so we never write a
        // NaN state that would later poison restart_from_iter. Files land at
        // output_path/atm_restart_<total_iter_count>.bin (per-iter names, no clobber).
        {
            constexpr int restart_save_stride = 100;
            if(restart_save_stride > 0 && total_iter_count > 0
               && total_iter_count % restart_save_stride == 0){
                bool clean = true;
                for(Array* a : restart_arrays()){
                    for(int i = 0; i < im && clean; i++)
                        for(int j = 0; j < jm && clean; j++)
                            for(int k = 0; k < km && clean; k++)
                                if(!AtomUtils::is_finite_safe(a->x[i][j][k])) clean = false;
                }
                if(clean)
                    save_state(total_iter_count);
                else
                    cout << "      AGCM: restart checkpoint SKIPPED at iter "
                         << total_iter_count << " — non-finite cell present (not clean)" << endl;
            }
        }

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
