#pragma once

#include "cHydrosphereModel.h"
#include "Utils.h"

#include <iostream>
#include <cmath>
#include <chrono>
#include <cstdio>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace AtomUtils;

class UtilsHyd {
public:
    explicit UtilsHyd(cHydrosphereModel& model)
        : m(model)
    {}

    // ------------------------------------------------------------------
    void resetArrays()
    {
        using namespace std;
        cout << endl << "      OGCM: reset_arrays" << endl;

        m.rad.initArray_1D(m.im, 1.0);
        m.the.initArray_1D(m.jm, 2.0);
        m.phi.initArray_1D(m.km, 3.0);
        m.aux_grad_v.initArray_1D(m.im, 4.0);
        m.aux_grad_w.initArray_1D(m.im, 5.0);

        // 2D arrays — zero-initialised, launched as parallel tasks
        Array_2D* arrays_2d[] = {
            &m.Bathymetry, &m.Upwelling, &m.Downwelling, &m.EkmanPumping,
            &m.salinity_evaporation, &m.Evaporation_Dalton, &m.Evaporation_Penman, &m.Evaporation,
            &m.Precipitation_2D, &m.precipitation_NASA, &m.temperature_NASA,
            &m.temp_reconst, &m.c_fix, &m.t_surf_fix, &m.v_wind, &m.w_wind,
            &m.psi_bt, &m.v_bt, &m.w_bt,
            &m.velocity_v_NASA, &m.velocity_w_NASA, &m.temp_landscape,
            &m.BuoyancyForce_2D,
            &m.vel_star     // friction velocity per (j,k) column — turbulence only
        };
        constexpr int n_2d = sizeof(arrays_2d) / sizeof(arrays_2d[0]);

        // 3D arrays — zero-initialised
        Array* arrays_3d_zero[] = {
            &m.h, &m.u, &m.v, &m.w,
            &m.un, &m.vn, &m.wn,
            &m.p_dyn,
            &m.rhs_t, &m.rhs_u, &m.rhs_v, &m.rhs_w, &m.rhs_c,
            &m.aux_u, &m.aux_v, &m.aux_w,
            &m.uf, &m.vf, &m.wf,   // Rhie-Chow face mass fluxes
            &m.BuoyancyForce, &m.CoriolisForce,
            &m.CentrifugalForce, &m.PresGradForce,
            // Turbulence fields — zero-initialised (inactive until TurbulenceHyd::init())
            &m.tke, &m.tken, &m.dis, &m.disn,
            &m.nue, &m.prod, &m.tke_source, &m.dis_source,
            &m.rhs_tke, &m.rhs_dis,
            // w-momentum-budget diagnostic fields (zero-initialised)
            &m.wbud_pgf, &m.wbud_cor, &m.wbud_adv, &m.wbud_diff,
            // v (meridional) momentum-budget diagnostic fields (zero-initialised)
            &m.vbud_pgf, &m.vbud_cor, &m.vbud_adv, &m.vbud_diff, &m.vbud_wind,
            // u (radial) momentum-budget diagnostic fields (zero-initialised)
            &m.ubud_pgf, &m.ubud_adv, &m.ubud_diff, &m.ubud_buoy, &m.ubud_cor
        };
        constexpr int n_3d_zero = sizeof(arrays_3d_zero) / sizeof(arrays_3d_zero[0]);

        // Serial initialisation. Array allocation is a one-time cost; the previous
        // OpenMP task-parallel version carried a data race (confirmed by
        // ThreadSanitizer, which concentrated its reports here) that was benign
        // only because the since-removed Salt_Finger/Diffusion/Balance arrays
        // happened to absorb the errant access — removing them let it corrupt a
        // live dynamics field, tipping a marginal deep-Southern-Ocean cell into a
        // CFL blow-up at iter 2 (multi-thread only; single-thread was always
        // stable). Serial init is deterministic and negligibly slower.
        for (int n = 0; n < n_2d; n++)
            arrays_2d[n]->initArray_2D(m.jm, m.km, 0.0);
        for (int n = 0; n < n_3d_zero; n++)
            arrays_3d_zero[n]->initArray(m.im, m.jm, m.km, 0.0);

        m.t.initArray(m.im, m.jm, m.km, 1.0);
        m.c.initArray(m.im, m.jm, m.km, 1.0);
        m.tn.initArray(m.im, m.jm, m.km, 1.0);
        m.cn.initArray(m.im, m.jm, m.km, 1.0);
        m.p_hydro.initArray(m.im, m.jm, m.km, 1.0);
        m.r_water.initArray(m.im, m.jm, m.km, m.r_0_water);
        m.r_salt_water.initArray(m.im, m.jm, m.km, m.r_0_water);

        for (auto& row : m.i_bathymetry)
            std::fill(row.begin(), row.end(), 0);

        cout << "      OGCM: reset_arrays ended" << endl;
    }

    // ------------------------------------------------------------------
    void storeIntermediateData3D(double coeff = 1.0)
    {
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                #pragma omp simd
                for (int k = 0; k < m.km; k++) {
                    m.un.x[i][j][k] = coeff * m.u.x[i][j][k];
                    m.vn.x[i][j][k] = coeff * m.v.x[i][j][k];
                    m.wn.x[i][j][k] = coeff * m.w.x[i][j][k];
                    m.tn.x[i][j][k] = coeff * m.t.x[i][j][k];
                    m.cn.x[i][j][k] = coeff * m.c.x[i][j][k];
                    if (m.use_turbulence_model) {
                        m.tken.x[i][j][k] = coeff * m.tke.x[i][j][k];
                        m.disn.x[i][j][k] = coeff * m.dis.x[i][j][k];
                    }
                }
            }
        }
    }

    // ------------------------------------------------------------------
    void writeFile(std::string& bathymetry_name, std::string& output_path,
                   bool is_final_result)
    {
        using namespace std;
        cout << endl << endl << endl << "      OGCM: write_file" << endl;

        // Radial slice at the top INTEGRATED cell (im-2), NOT the surface row
        // i=im-1: the latter is a boundary row that bcRadius overwrites by cubic
        // extrapolation of the three cells below (v,w) / zero-gradient copy (u) —
        // it carries no RHS dynamics, so a "surface" plot at im-1 is just an
        // extrapolated echo of im-2 (appears one level low). im-2 is the real
        // wind-forced dynamical top. Also emit im-1 for the surface scalar fields.
        // Files are STAMPED WITH total_iter_count, NOT iter_n. The hydrosphere loops
        // `for(iter_n = 1; iter_n <= nm; ...)` (cHydrosphereModel.cpp:421) — iter_n RESTARTS
        // AT 1 on a resume — so stamping with iter_n made a restart write ..._100/200/...
        // straight ON TOP of the original run's files (resuming from 400 silently destroyed
        // the 0->400 VTK; hit for real at Ma=100 and Ma=200 on 2026-07-17). total_iter_count
        // is the true cumulative iteration, and is what save_state already stamps restarts
        // with, so it never collides. On a FRESH run total_iter_count == iter_n, so filenames
        // are unchanged — a no-op except on restarts, where it is the fix.
        // (The ATMOSPHERE needs no such change: it loops `for(iter_n = iter_start; ...)` with
        // iter_start = restart_from_iter, so ITS iter_n is already absolute. Same root
        // asymmetry as the nm-semantics split: nm is ABSOLUTE for the atm, ADDITIONAL here.)
        m.paraview_vtk_radial(bathymetry_name, m.im - 2, m.total_iter_count);
        m.paraview_vtk_radial(bathymetry_name, m.im - 1, m.total_iter_count);

        int j_longal = 75;
        m.paraview_vtk_longal(bathymetry_name, j_longal, m.total_iter_count);

        int k_zonal = 185;
        m.paraview_vtk_zonal(bathymetry_name, k_zonal, m.total_iter_count);

        // Panorama VTS fires whenever iter_n is a multiple of panorama_print (see comment
        // in UtilsAtm.h::writeFile for why the panorama_cnt counter was replaced, and finally
        // removed on 2026-08-30 -- iter_n is now the only counter). The FIRING
        // stays on iter_n — that is the run-local cadence the caller asked for; only the
        // STAMP becomes total_iter_count, per the note above.
        if (m.paraview_panorama_vts_flag && m.panorama_print > 0
            && m.iter_n > 0 && m.iter_n % m.panorama_print == 0)
            m.paraview_panorama_vts(bathymetry_name, m.total_iter_count);

        m.HydrospherePlotData(bathymetry_name);
        m.HydrosphereSSTTransfer(bathymetry_name);   // reverse coupling channel for the atm<->hyd Picard loop

        cout << endl << "      OGCM: write_file ended" << endl;
    }

    // ------------------------------------------------------------------
    void valueLimitationHyd()
    {
        using namespace std;
        cout << endl << "      OGCM: ValueLimitationHyd" << endl;

        const double c_max = 45.0 / m.c_35;
        // Physical temperature bounds (non-dim t = T_K/t_0): -4 C .. 40 C. The
        // limiter previously bounded p_dyn and c but NOT t, so a deep-cell
        // temperature instability could run to hundreds of degrees -> NaN (seen
        // with the stronger-Coriolis experiment: deep high-lat T -> 275 C ->
        // NaN). Clamp t to a generous ocean range so a local instability stays
        // finite instead of poisoning the whole field. See
        // project_hydro_coriolis_dt_scaling.
        const double t_min = (m.t_0 -  4.0) / m.t_0;   //  -4 C
        const double t_max = (m.t_0 + 40.0) / m.t_0;   //  40 C

        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                #pragma omp simd
                for (int k = 0; k < m.km; k++) {

                    // p_dyn clamp loosened 0.1 -> 10 (2026-06-29). The old ±0.1
                    // cap was the proximate ENABLER of the iter-635 tropical-Indian-
                    // Ocean radial-velocity runaway: the [usplit] term-split showed
                    // p_dyn pegged at ±0.1 across the wind-driven-convergence zone,
                    // flattening the pressure field so the radial pressure gradient
                    // (pgf=-dpdr) was EXACTLY ZERO — no restoring force, letting the
                    // nonlinear self-advection u·∂u/∂r blow up. Opposing adv (~-300
                    // non-dim) needs p_dyn ~7-8 (dr~0.025), 75× over the old cap.
                    // ±10 lets the projection build a real radial PGF while still
                    // catching a genuine unbounded blow-up. See project_hydro_polar_blowup.
                    constexpr double p_dyn_cap = 10.0;
                    if (m.p_dyn.x[i][j][k] >=  p_dyn_cap)  m.p_dyn.x[i][j][k] =  p_dyn_cap;
                    if (m.p_dyn.x[i][j][k] < - p_dyn_cap)  m.p_dyn.x[i][j][k] = - p_dyn_cap;

                    if (m.c.x[i][j][k] >= c_max)  m.c.x[i][j][k] = c_max;
                    if (m.c.x[i][j][k] <  0.0)    m.c.x[i][j][k] = 0.0;

                    if (m.t.x[i][j][k] > t_max)  m.t.x[i][j][k] = t_max;
                    if (m.t.x[i][j][k] < t_min)  m.t.x[i][j][k] = t_min;
 
                    if (is_land(m.h, i, j, k)) {
                        m.c.x[i][j][k]     = 0.0;
                    }

                }
            }
        }


        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {

                int i_mount = m.i_bathymetry[j][k];

                if (is_land(m.h, i_mount, j, k)) {

                    m.t.x[0][j][k]     = m.t.x[i_mount][j][k];
                }
            }
        }

        cout << "      OGCM: ValueLimitationHyd ended" << endl;
    }

    // ------------------------------------------------------------------
    void findResiduumHyd()
    {
        using namespace std;
        cout << endl << "      OGCM: find_residuum_hyd" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        // The PREVIOUS call's residuum, captured before anything overwrites it -- see the
        // race note at the local_max update below. Same defect as the atmosphere's
        // findResiduumAtm; ATHAD has no ocean, so this site is the port's own find.
        const double residuum_prev = m.residuum_old;

        struct ErrorState { double val; int i, j, k; };
        ErrorState global_max = {0.0, 0, 0, 0};

        #pragma omp parallel
        {
            ErrorState local_max = {0.0, 0, 0, 0};

            #pragma omp for collapse(2) schedule(static) nowait
//            for (int j = 75; j <= 75; j++) {
//                for (int k = 180; k <= 180; k++) {
            for (int j = 1; j < m.jm; j++) {
                for (int k = 1; k < m.km; k++) {

                    double sinthe = sin(m.the.z[j]);
                    if (sinthe == 0.0) sinthe = 1.0e-5;

                    for (int i = 1; i < m.im-3; i++) {
                        if (!(is_water(m.h, i,   j, k) && is_water(m.h, i-1, j, k)
                           && is_water(m.h, i, j-1, k) && is_water(m.h, i, j+1, k)
                           && is_water(m.h, i, j, k-1) && is_water(m.h, i, j, k+1)))
                            continue;

                        double rm       = m.rad.z[i];
                        double rmsinthe = rm * sinthe;

                        double dudr   = (m.u.x[i+1][j][k] - m.u.x[i-1][j][k]) / m.dr;
                        double dvdthe = (m.v.x[i][j+1][k] - m.v.x[i][j-1][k])
                                        / (2.0 * rm * m.dthe);
                        double dwdphi = (m.w.x[i][j][k+1] - m.w.x[i][j][k-1])
                                        / (2.0 * rmsinthe * m.dphi);

                        double res = sqrt((dudr * dudr 
                                         + dvdthe * dvdthe 
                                         + dwdphi * dwdphi) / 3.0);

                        // `m.residuum_old = res;` used to sit here: a SHARED model member
                        // written from inside every thread's loop with no synchronisation.
                        // Whichever thread wrote last survived, so the value reported as
                        // "old" was not the previous call's residuum at all -- and it drives
                        // the declining/too-high line a human reads to judge convergence.
                        // Set once, after the reduction, below.
                        if (res > local_max.val) local_max = {res, i, j, k};
                    }
                }
            }

            #pragma omp critical
            {
                // Deterministic tie-break: a plain `>` lets the first thread into the
                // critical section win an equal maximum, so the reported error LOCATION
                // depended on thread arrival order. Prefer the lexicographically smallest
                // (i,j,k) so the choice is a property of the field, not of the schedule.
                const bool better = (local_max.val > global_max.val)
                                 || (local_max.val == global_max.val
                                     && (local_max.i < global_max.i
                                         || (local_max.i == global_max.i
                                             && (local_max.j < global_max.j
                                                 || (local_max.j == global_max.j
                                                     && local_max.k < global_max.k)))));
                if (better)
                    global_max = local_max;
            }
        }

        cout.precision(8);
        const bool declining = (residuum_prev - global_max.val) > 0.0;
        cout << endl
             << (declining
                 ? "      OGCM: find_residuum_hyd, absolute error declining .......................\n"
                 : "      OGCM: find_residuum_hyd, absolute error is too high .....................\n")
             << "      residuum_hyd = " << global_max.val
             << "      residuum_old = " << residuum_prev
             << "      eps_residuum = " << m.eps_residuum << endl
             << "      i_error = " << global_max.i
             << "   j_error = "    << global_max.j
             << "   k_error = "    << global_max.k << endl
             << "      relative error = " << fabs(global_max.val / residuum_prev - 1.0) << endl
             << "      absolute error = " << fabs(residuum_prev - global_max.val) << endl;

        // Hand this call's residuum to the next one. Written ONCE, outside every parallel
        // region, which is the whole of the fix.
        m.residuum_old = global_max.val;

        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for find_residuum_hyd\n", elapsed.count() * 1e-9);
        cout << "      OGCM: find_residuum_hyd ended" << endl;
    }

private:
    cHydrosphereModel& m;
};
