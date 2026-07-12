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
// Build the per-i non-uniform 3-point radial finite-difference coefficients from
// the actual (stretched) rad.z spacing. On a uniform grid these reduce exactly to
// the textbook (f[i+1]-f[i-1])/(2dr) and (f[i+1]-2f[i]+f[i-1])/dr2 stencils; on the
// stretched hydro grid they give the correct derivative w.r.t. the rad coordinate,
// which the constant-dr stencils did not. exp_rm and all metric factors are left
// untouched — they multiply d/drad and are unaffected by how d/drad is discretised.
void cHydrosphereModel::setupRadialStencilCoeffs(){
    rc1m.assign(im, 0.0); rc10.assign(im, 0.0); rc1p.assign(im, 0.0);
    rc2m.assign(im, 0.0); rc20.assign(im, 0.0); rc2p.assign(im, 0.0);
    rf10.assign(im, 0.0); rf11.assign(im, 0.0); rf12.assign(im, 0.0);
    rf20.assign(im, 0.0); rf21.assign(im, 0.0); rf22.assign(im, 0.0);

    // Central stencil (points i-1, i, i+1)
    for(int i = 1; i < im-1; i++){
        const double hm = rad.z[i]   - rad.z[i-1];     // spacing below
        const double hp = rad.z[i+1] - rad.z[i];       // spacing above
        const double hs = hm + hp;
        rc1m[i] = -hp / (hm * hs);
        rc10[i] =  (hp - hm) / (hm * hp);
        rc1p[i] =  hm / (hp * hs);
        rc2m[i] =  2.0 / (hm * hs);
        rc20[i] = -2.0 / (hm * hp);
        rc2p[i] =  2.0 / (hp * hs);
    }

    // Forward stencil (points i, i+1, i+2)
    for(int i = 0; i < im-2; i++){
        const double h1 = rad.z[i+1] - rad.z[i];
        const double h2 = rad.z[i+2] - rad.z[i+1];
        const double hs = h1 + h2;
        rf10[i] = -(2.0 * h1 + h2) / (h1 * hs);
        rf11[i] =  hs / (h1 * h2);
        rf12[i] = -h1 / (h2 * hs);
        rf20[i] =  2.0 / (h1 * hs);
        rf21[i] = -2.0 / (h1 * h2);
        rf22[i] =  2.0 / (h2 * hs);
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

    dt = dt_visc;                                                       //  no dimension; toggled to dt_inviscid during spin-up

    iter_n = 0;
    total_iter_count = 0;                                               // reset per time slice so the inviscid spin-up fires at the start of every Ma slice

    rad.StretchedCoordinates(im, r0, dr, dr_stretch);
    the.Coordinates(jm, the0, dthe);
    phi.Coordinates(km, phi0, dphi);

    setupRadialStencilCoeffs();   // non-uniform radial FD coefficients for the stretched rad grid

    read_Hydrosphere_Surface_Data(Ma);

    initTemperature(Ma);

    // Snapshot the prescribed surface SST as a sustained surface heat-flux forcing.
    // The atmosphere only INITIALISES the ocean surface temperature; rhs_t has no
    // surface heat-flux source term and bcRadius never re-pins the surface (it sets
    // only the insulating bottom t[0]=t[1]), so the ocean thermal state is unforced
    // and drifts — a hot grid-mode runaway without t-smoothing, or a cold collapse
    // with it (the warm surface bleeds into the vast cold deep ocean, nothing
    // re-warms it). Re-pinning t.x[im-1] to this fixed field every iteration in
    // run_3D_loop is a standard ocean-spin-up surface data BC; it is the thermal
    // analogue of the surface wind-stress forcing. Captured here, before the IC
    // damp_wiggles, so it holds the raw prescribed SST.
    for (int j = 0; j < jm; j++)
        for (int k = 0; k < km; k++)
            t_surf_fix.y[j][k] = t.x[im-1][j][k];

    AtomUtils::damp_wiggles(t, &i_bathymetry, true, true, true);

    initSalinity();

    // Fixed surface-salinity baseline for the E-P virtual salt flux, captured
    // ONCE from the prescribed IC (NASA SSS / paleo profile). SalinityEvaporation
    // applies the (bounded ±1 psu) evaporation offset ON TOP OF this fixed field
    // every iteration — it must NOT re-read the running surface salinity, or the
    // flux accumulates unboundedly (E>P pins to the 45-psu clamp, P>E to 0). This
    // is the salinity analogue of the t_surf_fix surface-temperature forcing above.
    for (int j = 0; j < jm; j++)
        for (int k = 0; k < km; k++)
            c_fix.y[j][k] = c.x[im-1][j][k];

    AtomUtils::damp_wiggles(c, &i_bathymetry, true, true, true);

    ThermoHyd(*this).SaltWaterDens();
//    UtilsHyd(*this).valueLimitationHyd();
    ThermoHyd(*this).SalinityEvaporation();
    ThermoHyd(*this).forces();
    LandOceanFraction();

    // Wind-driven surface-current IC. SverdrupGyre seeds the wind-stress-curl
    // geostrophic gyres (Sverdrup interior + western intensification), which is
    // what actually dominates the observed circulation and what the sustained
    // wind-stress body force spins up toward — a better-matched start than the
    // local Ekman spiral. EkmanSpiral() is the alternative (thin ageostrophic
    // drift only); swap the call to A/B them. Both are general over Ma bathymetry.
    SverdrupGyre();                                                     // wind-stress-curl gyres (Sverdrup/Stommel)
//  EkmanSpiral();                                                      // local Ekman drift (alternative IC)

    // THC seeds the deep thermohaline circulation (salinity c at depth; its
    // prescribed radial velocity u is discarded below), but it ALSO overwrites v/w
    // with hard-coded rectangular "surface flow" currents (spanning i_beg..surface)
    // that clobber the wind-driven SverdrupGyre gyres — producing 90-deg-cornered
    // blocks and a southward-shifted NH/SH separation in the iter-0 surface field.
    // Snapshot the clean gyre v/w, run THC for its deep c seeding, then RESTORE v/w:
    // this strips only THC's surface-current overwrite while keeping its deep
    // seeding. See project_hydro_gyre_vertical_diffusion.
    if (ocean_depth_mode == "deep") {
        Array v_gyre(v), w_gyre(w);                                    // clean SverdrupGyre currents
        ThermoHalineConveyorBelt(*this).run();                         // deep thermohaline seeding (c)
        v = v_gyre;                                                    // restore wind-driven gyres
        w = w_gyre;                                                    // (Array::operator= deep-copies)
    }

    // Vertical (radial) velocity is DIAGNOSTIC, not prescribed. SverdrupGyre sets
    // only the horizontal wind-driven currents; ThermoHalineConveyorBelt was
    // painting the radial velocity u directly with O(1 m/s) strips (IC_water=1,
    // factors up to -10) — ~1e4 too large for an ocean (real w ~ 1e-4 m/s) — which
    // seeded the near-surface radial-velocity runaway. The horizontal (v,w) currents
    // are the wind-driven SverdrupGyre gyres (THC's surface overwrite was stripped
    // above); discard THC's prescribed u; the divergence-free projection
    // (project_initial_velocity) then derives u from div(v,w) via continuity, the
    // physically correct source of ocean vertical velocity.
    for(int i = 0; i < im; i++)
        for(int j = 0; j < jm; j++)
            for(int k = 0; k < km; k++)
                u.x[i][j][k] = 0.0;

    AtomUtils::damp_wiggles(u, &i_bathymetry, true, true, true);
    AtomUtils::damp_wiggles(v, &i_bathymetry, true, true, true);
    AtomUtils::damp_wiggles(w, &i_bathymetry, true, true, true);

    UtilsHyd(*this).valueLimitationHyd();

    BC_Hyd(*this).bcRadius();
    BC_Hyd(*this).bcTheta();
    BC_Hyd(*this).bcPhi();
    BC_Hyd(*this).bcSolidGround();
    BC_Hyd(*this).boundaryCondition();

    // "laminar" is the descriptive name for the no-turbulence path; accept legacy "none"
    // and empty as synonyms (runs before the turbulence init and run_3D_loop flag setup).
    if (turb_model.empty() || turb_model == "none") turb_model = "laminar";

    if (turb_model != "laminar") {
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
    if (!turb_model.empty() && turb_model != "laminar") {
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

    // Project the initial velocity field divergence-free before time-stepping
    // (mirror of the atmosphere, cAtmosphereModel.cpp). The ocean IC (EkmanSpiral
    // + thermohaline currents) is NOT divergence-free, and the in-loop projection
    // can never fully catch up — leaving residual divergence that drove the
    // tropical-Indian-Ocean radial-velocity runaway (see project_hydro_polar_blowup).
    // Only on a fresh start; a restart resumes an already-evolved (projected) field.
    if(restart_from_iter < 0)
        PressureSolverHyd(*this).project_initial_velocity(200);

    // Seed the divergence-free face fluxes (uf/vf/wf) once before the loop so the
    // first advection step (before the first heavy block) transports with them
    // rather than a zero field. Also projects the (restored, on restart) field.
    PressureSolverHyd(*this).project_velocity(60);

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

        // Stage-resolved velocity-change decomposition: on checkpoint iters,
        // locate the deep blow-up cell and snapshot u,v,w after each iteration
        // stage (record_stage 0..4) to attribute the net Δ to RK4 / projection /
        // BC / polar-filter. See write_stage_budget / project_hydro_polar_blowup.
        stage_capture = (iter_n % checkpoint == 0);
        if(stage_capture){ locate_blowup_cell(); record_stage(0); }

        // Single hyd dynamical core (2026-07-08): the turbulent RHS reduces to LAMINAR when
        // use_turbulence_model is false (molecular diffusion, zeroed turbulence) and to EULER
        // when diffusion_ramp=0 (inviscid spin-up). The former separate laminar
        // solveRungeKutta_Hydrosphere / RHS_Hyd.cpp path was dropped — inviscid is now an
        // independent switch (diffusion_ramp), decoupled from the turbulence selection.
        solveRungeKutta_Hydrosphere_Turb();

        wbudget_capture = false;   // wbud_* now hold this iter's term split
        record_stage(1);           // stage 1: after RK4

        // First-NaN scanner — pinpoint the prognostic field + cell + iter of the
        // ~iter700 polar-velocity blow-up (project RESUME 2026-06-28). Runs right
        // after the RK4 solve, before BC/filters can relocate or mask the NaN, so
        // the report names the field that fails first (suspected v/w near N Pole).
        if(scanForNaN_hyd(total_iter_count, "post-RK4")){
            print_min_max_hyd();
            cout << endl << "      OGCM: aborting at first NaN (post-RK4 scan)." << endl;
            std::exit(1);
        }

        // Heavy block (pressure solver, salinity/thermo, turbulence update) runs every 2
        // iterations in the viscous phase, but only every 10 iterations during the inviscid
        // spin-up — dt_inviscid is 10× smaller so the physical-time cadence stays comparable.
        int heavy_block_stride = inviscid_phase ? 10 : 2;
        if(iter_n % heavy_block_stride == 0){

            // Continuity via the face-consistent fractional-step projection: solve
            // the masked-Neumann Poisson (source/Laplacian/gradient a single face-
            // consistent triple), correct the cell velocity, and fill the divergence-
            // free face fluxes uf/vf/wf that RHS_Hyd(_Turb) advect with. Replaces the
            // pressure-as-force run() (which left div(u)~2, radial u~0.5). Earlier
            // half-Rhie-Chow (face advection + run()'s rc1 cell projection) blew up;
            // the consistent source is blind to the cell checkerboard so it no longer
            // feeds back. See project_hydro_continuity_checkerboard.
            PressureSolverHyd(*this).project_velocity(60);
            record_stage(2);   // stage 2: after the pressure projection
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

        // Post-heavy-block NaN scan — the post-RK4 scan above showed the prognostic
        // fields stay finite while the EkmanPumping diagnostic goes NaN ~iter700, so
        // the NaN is born in the heavy block (PressureSolver / SaltWaterDens /
        // runDataHyd). This scan covers the derived 3D fields (aux_v/aux_w/
        // r_salt_water/r_water) + the 2D EkmanPumping to attribute the first NaN to
        // its real source, before valueLimitation/BC/filters can move or mask it.
        if(scanForNaN_hyd(total_iter_count, "post-heavy")){
            print_min_max_hyd();
            cout << endl << "      OGCM: aborting at first NaN (post-heavy scan)." << endl;
            std::exit(1);
        }

        UtilsHyd(*this).valueLimitationHyd();

        BC_Hyd(*this).bcRadius();
        BC_Hyd(*this).bcTheta();
        BC_Hyd(*this).bcPhi();
        BC_Hyd(*this).bcSolidGround();

        if (use_turbulence_model && !inviscid_phase) {
            TurbulenceHyd(*this).apply_wall_bc();
        }
        record_stage(3);   // stage 3: after valueLimitation + BC + wall BC

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

        record_stage(4);   // stage 4: after the polar zonal filter (final u,v,w)

        // Temperature wiggle damping — cure the deep high-latitude T instability.
        // The momentum/pressure/turbulence fields are de-checkerboarded every iter
        // (damp_wiggles on p_dyn/tke/dis/nue, polar filter on u/v/w) but T had NO
        // smoothing at all, so with the stronger wind-driven circulation a growing
        // 2Δ (checkerboard) temperature mode at the deepest cells ran away to
        // hundreds of degrees -> NaN (the T-clamp only masked it). A 1-2-1 Shapiro
        // step is curvature-selective: it damps the 2Δ mode strongly while leaving
        // the smooth thermocline essentially unchanged. Applied every iter, all
        // axes (the deep mode is not purely vertical), like the p_dyn treatment.
        AtomUtils::damp_wiggles(t, &i_bathymetry, true, true, true);

        // Sustained surface heat-flux forcing: re-pin the ocean surface temperature
        // to the prescribed atmospheric SST (snapshot in t_surf_fix at init) every
        // iteration — a Dirichlet surface data BC. rhs_t carries no surface heat
        // source and the surface BC never re-pins t.x[im-1], so without this the
        // ocean thermal state is unforced and drifts (cold collapse / hot runaway,
        // the corvalid/windtest failure). Applied after damp_wiggles so the warm
        // surface is the final word each iter; the smoother then only mixes the
        // interior (deep-checkerboard protection). Ocean cells only — land surface
        // temperature is handled by bcSolidGround. Thermal analogue of the wind-stress
        // forcing (see project_hydro_no_surface_heat_flux).
        for (int j = 0; j < jm; j++)
            for (int k = 0; k < km; k++)
                if (is_water(h, im-1, j, k))
                    t.x[im-1][j][k] = t_surf_fix.y[j][k];

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

            // Local per-longitude, depth-integrated MERIDIONAL (v) momentum
            // budget at gyre latitudes (uses the vbud_* split) — diagnoses why
            // the gyres close on the EASTERN boundary instead of forming a
            // western-boundary current. See project_hydro_eastern_boundary_current.
            write_v_momentum_budget(total_iter_count);

            // Radial (u) momentum budget at the self-located deep blow-up cell
            // (uses the ubud_* split captured this iter) — appends one row to
            // deep_momentum_budget.csv + a console [ubudget] line. Diagnoses the
            // ~1500-iter polar radial-velocity runaway. See project_hydro_polar_blowup.
            write_deep_momentum_budget(total_iter_count);

            // Stage-resolved net-Δ(u,v,w) decomposition at the same cell —
            // attributes the growth to RK4 / projection / BC / polar-filter.
            write_stage_budget(total_iter_count);
        }
    }  // end iter_n


    cout << endl << "      OGCM: run_3D_loop ended ..........................." << endl;
    return;
}
/*
*
*/

