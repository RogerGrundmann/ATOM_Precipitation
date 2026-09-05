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

#include <cstdlib>
#include <string>

#include "cHydrosphereModel.h"
#include "BC_Hyd.h"
#include "PressureSolverHyd.h"
#include "ThermoHyd.h"
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
    initMetricRadius();                      // HYD_METRIC_RADIUS; must follow the coordinates
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

    // ==================================================================================
    // HYD_T_FREEZE_SFC -- give the PRESCRIBED SURFACE FORCING a salinity-dependent
    // freezing floor, in place of the constant t_pole_salt. Default 0 = shipped and
    // bit-identical.
    //
    // THE SURFACE ALREADY HAS A FREEZING FLOOR AND IT IS A CONSTANT. InitValues_Hyd's
    // initTemperature clamps the prescribed SST at t_pole_salt = 0.993 nd = -1.912 C,
    // which is T_f for about 34.8 psu, and t_surf_fix then re-pins the surface row to
    // that field every iteration -- AFTER valueLimitationHyd has run. So the interior
    // floor (HYD_T_FREEZE) cannot reach the surface at all: measured 2026-09-05 over
    // iterations 300-500, mean dT at i = im-1 is +0.00000 C with 101 cells changed,
    // against +0.274 C at the seafloor, and 2 859 of the 2 961 cells sitting at exactly
    // -1.912 C are the surface pin.
    //
    // AND A CONSTANT FLOOR IS WRONG WHEREVER THE WATER IS NOT 34.8 psu. Fresher water
    // freezes WARMER: at 31.6 psu T_f is -1.731 C, so the pin holds such a cell 0.18 C
    // BELOW its own freezing point and manufactures supercooled surface water. That is
    // where the residual supercooling of the HYD_T_FREEZE arm lives -- 2 731 of its
    // 15 282 remaining supercooled cells are AT the surface. This raises the pin to
    // max(t_surf_fix, T_f(S)), so it can only ever warm, never cool.
    //
    // Applied to the PRESCRIBED fields and once, here rather than in the loop, because
    // t_surf_fix is a fixed forcing snapshot by design: it is floored against c_fix, the
    // prescribed surface salinity captured on the line above, not against the running
    // c.x[im-1] that SalinityEvaporation moves. (The floor is idempotent, so re-reading
    // the running field would not accumulate the way c_fix's offset would -- it would
    // simply make the forcing time-dependent, which is the one property this field has
    // to keep.) Note the capture of t_surf_fix itself is 20 lines above and runs BEFORE
    // initSalinity, so this could not have been folded into it: there is no surface
    // salinity to floor against at that point.
    //
    // Cells below 5 psu are left alone, the same guard the interior floor carries: ~12 %
    // of this ocean is spuriously fresh (the salinity-IC defect) and T_f(0) = 0 C would
    // warm the Arctic surface by four degrees.
    // ==================================================================================
    {
        const char* e_sfc = getenv("HYD_T_FREEZE_SFC");
        if (e_sfc && atoi(e_sfc) != 0){
            int    n_seen = 0, n_fresh = 0, n_raised = 0;
            double raise_max = 0.0;
            double t_sfc_min = 1.0e30;            // coldest prescribed SST, [C]
            double margin_max = -1.0e30;          // largest T_f - t_surf_fix, [C], even if <= 0
            for (int j = 0; j < jm; j++)
                for (int k = 0; k < km; k++){
                    if (!is_water(h, im-1, j, k))  continue;
                    const double t_c = (t_surf_fix.y[j][k] - 1.0) * t_0;
                    if (t_c < t_sfc_min)  t_sfc_min = t_c;
                    const double S = c_fix.y[j][k] * c_35;               // psu
                    if (S < 5.0){ n_fresh++; continue; }                 // fresh: leave alone
                    n_seen++;
                    const double Tf   = seawater_freezing_point_C(S);
                    const double cand = (t_0 + Tf) / t_0;
                    const double d    = (cand - t_surf_fix.y[j][k]) * t_0;
                    if (d > margin_max)  margin_max = d;
                    if (cand > t_surf_fix.y[j][k]){
                        if (d > raise_max)  raise_max = d;
                        t_surf_fix.y[j][k] = cand;
                        n_raised++;
                    }
                }
            // Report the whole census, not just the count. A null here is ambiguous
            // otherwise -- "the floor found nothing to do" and "the floor is not connected"
            // print the same line -- and the FIRST arm of this knob was a null for a third
            // reason again: the atmosphere transfer file was absent, so the prescribed SST
            // came from the NASA field whose minimum is -0.01 C, 1.9 C above any freezing
            // point. margin_max says which: it is the largest T_f - t_surf_fix over the
            // scored cells, so a value near 0 means the floor is only just inactive and a
            // large negative one means the prescribed SST is nowhere near freezing.
            cout << "      OGCM: HYD_T_FREEZE_SFC=1  surface SST forcing floored at T_f(S): "
                 << n_raised << " of " << n_seen << " salty ocean cells raised (" << n_fresh
                 << " skipped as < 5 psu), largest raise " << raise_max
                 << " C, coldest prescribed SST " << t_sfc_min
                 << " C, max (T_f - SST) " << margin_max << " C" << endl;
        }
    }

    AtomUtils::damp_wiggles(c, &i_bathymetry, true, true, true);

    ThermoHyd(*this).SaltWaterDens();
//    UtilsHyd(*this).valueLimitationHyd();
    ThermoHyd(*this).SalinityEvaporation();
    ThermoHyd(*this).Forces();
    LandOceanFraction();

    // Wind-driven surface-current IC. SverdrupGyre seeds the wind-stress-curl
    // geostrophic gyres (Sverdrup interior + western intensification), which is
    // what actually dominates the observed circulation and what the sustained
    // wind-stress body force spins up toward — a better-matched start than the
    // local Ekman spiral. EkmanSpiral() is the alternative (thin ageostrophic
    // drift only); swap the call to A/B them. Both are general over Ma bathymetry.
    SverdrupGyre();                                                     // wind-stress-curl gyres (Sverdrup/Stommel)
//  EkmanSpiral();                                                      // local Ekman drift (alternative IC)

    // ThermoHalineConveyorBelt REMOVED 2026-07-12. After its surface v/w overwrite
    // was stripped (204f50d) and its prescribed radial u discarded (below), its ONLY
    // remaining effect was a hard-coded deep-salinity seeding painted on top of the
    // observed/paleo initSalinity field — an externally-forced "conveyor" that the
    // thermohaline circulation should develop on its own. Dropped entirely; the
    // salinity is initSalinity's (NASA SSS / paleo + Gill interior).

    // Vertical (radial) velocity is DIAGNOSTIC, not prescribed: SverdrupGyre sets
    // only the horizontal wind-driven currents, so zero any u here; the
    // divergence-free projection (project_initial_velocity) then derives u from
    // div(v,w) via continuity, the physically correct source of ocean vertical velocity.
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
    // Barotropic (rigid-lid Stommel) solve — the wind-driven barotropic transport
    // the Chorin projection cannot build. Wind is fixed during a hyd run, so solve
    // once here for u_bt (m.v_bt/w_bt); apply_barotropic_mode_split() then imposes it
    // as the column depth-mean each iter (velocity mode split, not a force).
    // Runs in BOTH depth modes: u_bt = transport/H now uses the true seafloor
    // depth (Bathymetry), so shallow mode gets the same physical u_bt as deep
    // (~38 cm/s, not the ~80 cm/s CFL runaway that came from dividing by the
    // truncated 200 m column). See project_barotropic().
    PressureSolverHyd(*this).project_barotropic();
    if(Ma == 0) run_3D_loop(Ma);                                        // iterational 3D loop to solve variables in 4-step Runge-Kutta time scheme


    cout << endl << endl;
    if(Ma > 0) run_3D_loop(Ma);                                         // iterational 3D loop to solve variables in 4-step Runge-Kutta time scheme
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
void cHydrosphereModel::initMetricRadius(){
    // One rad.z unit in metres. The ocean grid maps a radial step to a depth as
    // (rad.z[i+1]-rad.z[i])*L_hyd (see the dz lines in this file), and rad.z spans exactly 1.0
    // from r0 = 1.0, so one unit is L_hyd metres EXACTLY and uniformly -- simpler than the
    // atmosphere, whose exponential stretch makes its metricShellLength() an average.
    static const double r_km = [this](){
        const char* e = getenv("HYD_METRIC_RADIUS");
        return e ? atof(e) : 0.0; }();          // DEFAULT 0 = off, unlike the atmosphere's

    if(r_km <= 0.0){
        m_metric_r0 = 0.0;
        cout << "      OGCM: HYD_METRIC_RADIUS off - horizontal metric left on the grid"
             << " coordinate (rad.z ~ 1..2, i.e. a " << L_hyd << " m planet)" << endl;
        return;
    }
    m_metric_r0 = r_km * 1.0e3 / L_hyd;
    cout << "      OGCM: HYD_METRIC_RADIUS = " << r_km << " km, one rad.z unit = " << L_hyd
         << " m  ->  metric r0 = " << m_metric_r0
         << " (grid r0 stays " << rad.z[0] << ")" << endl;
}
/*
*
*/
void cHydrosphereModel::run_3D_loop(int Ma){

cout << endl << endl << endl << "      OGCM: run_3D_loop ..........................." << endl;

    // [RUN CONFIG] for the hydrosphere's environment knobs, unconditional, so no run has to be
    // reconstructed from a shell history. The atmosphere gained this banner on 2026-09-02 after
    // two "flipped on" claims in CLAUDE.md turned out to be wrong about a default; the ocean had
    // no knobs to print until now. `*` marks a compiled-in default not set in the environment.
    {
        auto ev = [](const char* k, const char* dflt){
            const char* e = getenv(k);
            return e ? std::string(e) : (std::string(dflt) + "*");
        };
        cout << "      OGCM: [RUN CONFIG] knobs:  PHYDRO_SALT=" << ev("HYD_PHYDRO_SALT", "0")
             << "  BAROCLINIC_PGF="                            << ev("HYD_BAROCLINIC_PGF", "0.0")
             << "  METRIC_RADIUS="                             << ev("HYD_METRIC_RADIUS", "0")
             << "  RUN_NEUMANN="                               << ev("HYD_RUN_NEUMANN", "0")
             << "  BC_SECOND_ORDER="                           << ev("HYD_BC_SECOND_ORDER", "0")
             << "  LINE_SOLVE="                                << ev("HYD_LINE_SOLVE", "0")
             << "  T_FREEZE="                                  << ev("HYD_T_FREEZE", "0")
             << "  T_FREEZE_SFC="                              << ev("HYD_T_FREEZE_SFC", "0")
             << "  A_H="                                       << ev("HYD_A_H", "0")
             << "  A_H_BIHARM="                                << ev("HYD_A_H_BIHARM", "0")
             << "  SFC_FLUX="                                  << ev("HYD_SFC_FLUX", "0")
             << "   (* = compiled-in default, not set in the environment)" << endl;
        cout << "      OGCM: [SCALES] L_hyd = " << L_hyd << " m   u_0 = " << u_0
             << " m/s   L_hyd/u_0 = " << L_hyd / u_0 << " s   one iteration = "
             << dt * L_hyd / u_0 << " s" << endl;

        // Horizontal-mixing scales, printed unconditionally because every argument about
        // HYD_METRIC_RADIUS and HYD_A_H has had to re-derive them. dx is the equatorial
        // 1-degree cell; the 2*dx e-folding is 1/(A_h*k^2) with k = pi/dx.
        {
            const double dx    = 2.0 * M_PI * r_Earth * 1.0e3 / 360.0;
            const double k2    = (M_PI / dx) * (M_PI / dx);
            const double s_per = dt * L_hyd / u_0;
            // what the GRID metric implies today: rm runs 1..2, so the horizontal metric
            // radius is rm*L_hyd unless HYD_METRIC_RADIUS moved it.
            const double rh_sfc = metricRadius(2.0) * L_hyd;                  // [m]
            cout << "      OGCM: [SCALES] dx(eq) = " << dx / 1.0e3
                 << " km   horizontal metric radius at the surface = " << rh_sfc
                 << " m (Earth: " << r_Earth * 1.0e3 << ")" << endl;
            // domain-scale wavenumber, for the selectivity line below
            const double kL  = 2.0 * M_PI / 2.0e7;
            const double a_h = [](){ const char* e = getenv("HYD_A_H"); return e ? atof(e) : 0.0; }();
            if(a_h != 0.0)
                cout << "      OGCM: [SCALES] A_H = " << a_h << " m2/s -> 2dx e-folding "
                     << 1.0 / (a_h * k2) << " s = " << 1.0 / (a_h * k2) / s_per
                     << " iterations,  domain-scale e-folding "
                     << 1.0 / (a_h * kL * kL) / 86400.0 << " d" << endl;
            // air-sea flux timescale: tau = rho*cp*dz/c_H over the TOP PROGNOSTIC layer,
            // which is where HYD_SFC_FLUX puts it (im-1 is the prescribed skin).
            const double dz_top = (rad.z[im-1] - rad.z[im-2]) * L_hyd;
            const double c_H = [](){ const char* e = getenv("HYD_SFC_FLUX"); return e ? atof(e) : 0.0; }();
            cout << "      OGCM: [SCALES] top prognostic layer (i = im-2) = " << dz_top
                 << " m" << endl;
            if(c_H != 0.0){
                const double tau = r_0_water * cp_w * dz_top / c_H;
                cout << "      OGCM: [SCALES] SFC_FLUX c_H = " << c_H
                     << " W/m2/K -> tau = " << tau << " s = " << tau / 86400.0
                     << " d = " << tau / s_per << " iterations" << endl;
            }
            const double b_h = [](){ const char* e = getenv("HYD_A_H_BIHARM"); return e ? atof(e) : 0.0; }();
            if(b_h != 0.0)
                cout << "      OGCM: [SCALES] A_H_BIHARM = " << b_h << " m4/s -> 2dx e-folding "
                     << 1.0 / (b_h * k2 * k2) << " s = " << 1.0 / (b_h * k2 * k2) / s_per
                     << " iterations,  domain-scale e-folding "
                     << 1.0 / (b_h * kL * kL * kL * kL) / 86400.0 << " d" << endl;
        }
    }

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
        load_state(restart_from_iter, Ma);

    // ---- Dump-only mode: regenerate a checkpoint's ParaView output, don't advance ----
    // nm == 0 with a restart re-emits the VTK/panorama for the RESTORED state and returns.
    // Placed immediately after load_state and BEFORE the pre-loop projections below, so
    // what is written is exactly the saved field, untouched. The stamp is total_iter_count
    // (restored by load_state), so the files carry their true cumulative iteration.
    // iter_n is set to match only so writeFile's `iter_n % panorama_print` cadence test
    // fires — writeFile itself stamps from total_iter_count (see UtilsHyd.h::writeFile).
    // Exists because the old iter_n-based naming let a restart overwrite the original
    // run's VTK, and the surviving restart .bin files are the only intact record; this
    // rebuilds the lost/mislabelled output from them. nm == 0 was previously a no-op
    // (`for(iter_n = 1; iter_n <= 0; ...)` never executes), so this steals nothing.
    if(nm == 0 && restart_from_iter >= 0){
        iter_n = total_iter_count;
        cout << endl << "      OGCM: DUMP-ONLY — re-emitting ParaView output for restored"
             << " total_iter_count " << total_iter_count << " (no time stepping)" << endl;
        ThermoHyd(*this).runDataHyd();
        print_min_max_hyd();
        UtilsHyd(*this).writeFile(bathymetry_name, output_path, true);
        return;
    }

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
    PressureSolverHyd(*this).apply_barotropic_mode_split();   // seed the barotropic mode into the IC (both depth modes)

    // ---- Convergence monitor (report-only), mirror of the atmosphere's ------------
    // Every conv_stride iters, sample two slow volume-integral metrics — mean ocean
    // temperature (thermodynamics) and mean kinetic energy (circulation) — and report the
    // drift of the last conv_window samples' mean vs the previous window's. Window
    // averaging cancels per-iter wobble so the drift is the SECULAR trend. Never stops the
    // run (nm is the hard ceiling); also logged to output_path/convergence.csv.
    //
    // Two deliberate DIFFERENCES from the atmosphere's monitor:
    //  1. OCEAN MASK. Land/unused cells are excluded (is_water). The atmosphere averages
    //     its whole domain; here land cells hold 0 and would drag both means.
    //  2. TEMPERATURE DRIFT IS REPORTED IN ABSOLUTE K/window, and the converged flag uses
    //     it. A PERCENT drift relative to Kelvin is meaningless for an ocean: the mean sits
    //     near 287 K, so a real cooling of ~0.07 K/100 iters is only ~0.024 % — under the
    //     atmosphere's 1 % tolerance, which would report CONVERGED while the ocean visibly
    //     cools ~1 K over 900 iters. KE keeps the percent test (its zero IS physical, so a
    //     relative drift is meaningful there). The percent T drift is still logged for
    //     comparability with the atmosphere's CSV.
    // NOTE the pinned surface row (t.x[im-1] = t_surf_fix every iter) is INCLUDED: it is
    // ~1/41 of the column and constant, so it only damps the reported T drift by ~5%.
    constexpr int    conv_stride   = 25;     // sample cadence [iters]
    constexpr int    conv_window   = 4;      // samples per window (-> 100-iter trailing window)
    constexpr double conv_tol_pct  = 1.0;    // KE convergence threshold [% drift per window]
    constexpr double conv_tol_T_K  = 0.01;   // T convergence threshold [K drift per window]
    std::vector<double> conv_histT, conv_histKE;
    auto ocean_mean_T = [&]() -> double {
        double sum = 0.0, wsum = 0.0;
        #pragma omp parallel for reduction(+:sum,wsum) schedule(static)
        for(int j = 0; j < jm; j++){
            const double coslat = sin(the.z[j]);                 // the.z = colatitude
            for(int i = 0; i < im; i++){
                const double dz = (i < im-1) ? (rad.z[i+1] - rad.z[i]) * L_hyd : 0.0;
                const double wt = coslat * dz;
                for(int k = 0; k < km; k++)
                    if(is_water(h, i, j, k)){ sum += wt * t.x[i][j][k] * t_0; wsum += wt; }
            }
        }
        return (wsum > 0.0) ? sum / wsum : 0.0;                  // [K]
    };
    auto ocean_mean_KE = [&]() -> double {
        double sum = 0.0, wsum = 0.0;
        #pragma omp parallel for reduction(+:sum,wsum) schedule(static)
        for(int j = 0; j < jm; j++){
            const double coslat = sin(the.z[j]);
            for(int i = 0; i < im; i++){
                const double dz = (i < im-1) ? (rad.z[i+1] - rad.z[i]) * L_hyd : 0.0;
                const double wt = coslat * dz;
                for(int k = 0; k < km; k++){
                    if(!is_water(h, i, j, k)) continue;
                    const double uu = u.x[i][j][k], vv = v.x[i][j][k], ww = w.x[i][j][k];
                    sum += wt * 0.5 * (uu*uu + vv*vv + ww*ww); wsum += wt;
                }
            }
        }
        return (wsum > 0.0) ? (sum / wsum) * (u_0 * u_0) : 0.0;  // [m2/s2] (velocities are /u_0)
    };
    // Returns the trailing-window drift of `metric`, or -1 while still warming up
    // (fewer than 2 full windows). `abs_out` receives the drift in the metric's own units.
    auto conv_drift = [&](std::vector<double>& hist, double metric, double& abs_out) -> double {
        hist.push_back(metric);
        const int n = (int)hist.size();
        if(n < 2*conv_window){ abs_out = -1.0; return -1.0; }
        double cur = 0.0, prev = 0.0;
        for(int i = n-conv_window;   i < n;             i++) cur  += hist[i];
        for(int i = n-2*conv_window; i < n-conv_window; i++) prev += hist[i];
        cur /= conv_window; prev /= conv_window;
        abs_out = std::fabs(cur - prev);
        return 100.0 * abs_out / std::max(1e-12, std::fabs(cur));
    };
    if(restart_from_iter < 0){                                   // fresh run: (re)create with header
        std::ofstream cf(output_path + "/convergence_hyd.csv");
        if(cf) cf << "iter,mean_T_K,drift_T_pct,drift_T_K,mean_KE_m2s2,drift_KE_pct,converged\n";
    }

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
            PressureSolverHyd(*this).apply_barotropic_mode_split();   // impose wind-driven barotropic mode (both depth modes)
            record_stage(2);   // stage 2: after the pressure projection + barotropic mode split
            AtomUtils::damp_wiggles(p_dyn, &i_bathymetry, true, true, true);

            ThermoHyd(*this).SaltWaterDens();
//            UtilsHyd(*this).valueLimitationHyd();
            ThermoHyd(*this).SalinityEvaporation();
            ThermoHyd(*this).runDataHyd();
            ThermoHyd(*this).Forces();

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

        // Sustained surface heat-flux forcing: relax the ocean surface temperature
        // toward the prescribed atmospheric SST (snapshot in t_surf_fix at init) every
        // iteration. rhs_t carries no surface heat source and the surface BC never
        // re-pins t.x[im-1], so without this the ocean thermal state is unforced and
        // drifts (cold collapse / hot runaway, the corvalid/windtest failure). Applied
        // after damp_wiggles so the surface forcing is the final word each iter; the
        // smoother then only mixes the interior (deep-checkerboard protection). Ocean
        // cells only — land surface temperature is handled by bcSolidGround. Thermal
        // analogue of the wind-stress forcing (see project_hydro_no_surface_heat_flux).
        //
        // sst_relax_alpha selects the BC strength:
        //   1.0 (default) = hard Dirichlet re-pin, i.e. exactly t.x = t_surf_fix — the
        //                   long-standing behaviour, reproduced BIT-IDENTICALLY.
        //   <1            = Haney-type flux BC. The surface keeps a fraction of its own
        //                   advective/diffusive tendency, so ocean dynamics can build an
        //                   SST anomaly instead of having it erased every iteration. The
        //                   steady anomaly scales as ~dt*tendency/alpha.
        //   0             = unforced; this is the cold-collapse regime, do not ship it.
        // NOTE this is a NUMERICAL knob, not a physical restoring timescale. A calibrated
        // Haney tau (~59 d at 40 W/m^2/K over a 50 m mixed layer) is tau_nd ~6.1e3 given
        // L_hyd/u_0 = 833 s, hence alpha = dt/tau_nd ~1.6e-8 with dt = 1e-4 — numerically
        // identical to unforced. Physical restoring is unreachable in this model's
        // iteration budget, so alpha is tuned for stability like the atmosphere's
        // omega_teq (which exists for exactly the same reason: Held-Suarez at physical
        // k_s is inert here).
        for (int j = 0; j < jm; j++)
            for (int k = 0; k < km; k++)
                if (is_water(h, im-1, j, k)){
                    if (sst_relax_alpha >= 1.0)
                        t.x[im-1][j][k] = t_surf_fix.y[j][k];       // exact legacy path
                    else
                        t.x[im-1][j][k] += sst_relax_alpha
                                         * (t_surf_fix.y[j][k] - t.x[im-1][j][k]);
                }

//        UtilsHyd(*this).findResiduumHyd();

        UtilsHyd(*this).storeIntermediateData3D(1.0);

        // Convergence monitor: sample the slow integral metrics and report secular drift.
        // Report-only — never stops the run. See the setup block before the loop for why
        // T is gated on an ABSOLUTE K drift while KE keeps the percent test.
        if(conv_stride > 0 && total_iter_count > 0 && total_iter_count % conv_stride == 0){
            double dT_K = 0.0, dKE_abs = 0.0;
            const double mT   = ocean_mean_T();
            const double mKE  = ocean_mean_KE();
            const double dT   = conv_drift(conv_histT,  mT,  dT_K);
            const double dKE  = conv_drift(conv_histKE, mKE, dKE_abs);
            const bool warming = (dT < 0.0) || (dKE < 0.0);      // < 2 windows of samples yet
            const bool converged = !warming && dT_K < conv_tol_T_K && dKE < conv_tol_pct;
            cout << "      [conv] iter " << total_iter_count
                 << "  T=" << mT << " K  KE=" << mKE << " m2/s2";
            if(warming) cout << "  (warming up)";
            else cout << "  drift: T=" << dT_K << " K (" << dT << "%)  KE=" << dKE
                      << "%/window" << (converged ? "  -> CONVERGED (both < tol)" : "");
            cout << std::endl;
            std::ofstream cf(output_path + "/convergence_hyd.csv", std::ios::app);
            if(cf) cf << total_iter_count << "," << mT << "," << (warming?0.0:dT) << ","
                      << (warming?0.0:dT_K) << "," << mKE << "," << (warming?0.0:dKE) << ","
                      << (converged?1:0) << "\n";
        }

//        print_min_max_hyd();

        if(iter_n % checkpoint == 0){
            print_min_max_hyd();

            UtilsHyd(*this).writeFile(bathymetry_name, output_path, true);
            cout << endl << "      OGCM: write_file in run_3D_loop ......................." << endl;

            // Binary restart checkpoint. Grouped with the panorama VTS cadence
            // (every panorama_print iters) rather than the VTK cadence (`checkpoint`),
            // so the large .bin + .vts share the coarser interval while the VTK
            // slices write more often. Falls back to every checkpoint if disabled.
            if (panorama_print <= 0 || iter_n % panorama_print == 0)
                save_state(total_iter_count, Ma);

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

