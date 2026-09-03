#pragma once

#include "cAtmosphereModel.h"
#include "Utils.h"

#include <vector>
#include <cstdint>
#include <cstdlib>   // getenv/atof for the ATM_POISSON_METRIC_FIX A/B knob
#include <cmath>
#include <chrono>
#include <cstdio>
#include <iostream>
#include <iomanip>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace AtomUtils;

class PressureSolverAtm {
public:
    explicit PressureSolverAtm(cAtmosphereModel& model)
        : m(model)
    {}

    // sweeps_override > 0 replaces the ATM_PRESS_SWEEPS knob for THIS call only. It exists so
    // the initial projection can relax harder than the time loop without the two counts being
    // entangled -- see project_initial_velocity and ATHAD README item 68, where setting the
    // time-loop knob silently changed the startup state as well and an attribution was lost.
    void run(bool verbose = true, int sweeps_override = -1)
    {
        using namespace std;
        if (verbose) cout << endl << endl << endl << "      ATOM: PressureSolverAtm" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        // precompute sin(the) table — only depends on j, avoids redundant sin() calls
        std::vector<double> sinthe_table(m.jm);
        for (int j = 0; j < m.jm; j++) {
            sinthe_table[j] = sin(m.the.z[j]);
            if (sinthe_table[j] < 0.55) sinthe_table[j] = 0.55;   // metric floor ~57° — keep in sync with the RHS geometry floor (RungeKutta_Atm*.cpp)
        }

        // Precompute land mask — eliminates repeated function call overhead
        // Allocate flat mask: 1 = land, 0 = air. Also enforce no-penetration on the
        // provisional velocity: bcSolidGround masks u/v/w at land but NOT aux_u/v/w, and
        // RHS_Atm writes aux unconditionally — so land cells carry stale aux (advection of
        // the adjacent air jet + the blown-up p_dyn gradient). The fluid-cell divergence
        // source uses central differences at cliff edges, differencing against that stale
        // land aux. Zero it here so every sweep sees a clean solid wall.
        std::vector<int8_t> land(m.im * m.jm * m.km);
        #pragma omp parallel for collapse(3)
        for (int i = 0; i < m.im; i++)
            for (int j = 0; j < m.jm; j++)
                for (int k = 0; k < m.km; k++) {
                    const bool isl = is_land(m.h, i, j, k);
                    land[i*m.jm*m.km + j*m.km + k] = isl ? 1 : 0;
                    if (isl) {
                        m.aux_u.x[i][j][k] = 0.0;
                        m.aux_v.x[i][j][k] = 0.0;
                        m.aux_w.x[i][j][k] = 0.0;
                    }
                }

        #define LAND(i,j,k) land[(i)*m.jm*m.km + (j)*m.km + (k)]

        // Fuse the three boundary loops into one pass
        #pragma omp parallel for collapse(2)
        for (int j = 1; j < m.jm-1; j++) {
            for (int k = 1; k < m.km-1; k++) {
                m.aux_u.x[0][j][k]      = m.c43 * m.aux_u.x[1][j][k]      - m.c13 * m.aux_u.x[2][j][k];
                m.aux_u.x[m.im-1][j][k] = m.c43 * m.aux_u.x[m.im-2][j][k] - m.c13 * m.aux_u.x[m.im-3][j][k];
                m.aux_v.x[0][j][k]      = m.c43 * m.aux_v.x[1][j][k]      - m.c13 * m.aux_v.x[2][j][k];
                m.aux_v.x[m.im-1][j][k] = m.c43 * m.aux_v.x[m.im-2][j][k] - m.c13 * m.aux_v.x[m.im-3][j][k];
                m.aux_w.x[0][j][k]      = m.c43 * m.aux_w.x[1][j][k]      - m.c13 * m.aux_w.x[2][j][k];
                m.aux_w.x[m.im-1][j][k] = m.c43 * m.aux_w.x[m.im-2][j][k] - m.c13 * m.aux_w.x[m.im-3][j][k];
            }
        }

        // Grid-spacing reciprocals — constant for the entire grid
        const double inv_2dr   = 1.0 / (2.0 * m.dr);
        const double inv_2dthe = 1.0 / (2.0 * m.dthe);
        const double inv_2dphi = 1.0 / (2.0 * m.dphi);
        const double inv_dr2   = 1.0 / (m.dr   * m.dr);
        const double inv_dthe2 = 1.0 / (m.dthe * m.dthe);
        const double inv_dphi2 = 1.0 / (m.dphi * m.dphi);
        const double inv_dthe  = 1.0 / m.dthe;
        const double inv_dphi  = 1.0 / m.dphi;

        // Cap on the divergence source of the pressure Poisson update. Over steep
        // orography (cliff faces: Patagonian Andes 49°S/74°W, Atlas 36°N/2°E) the
        // central-difference divergence of the provisional velocity can spike, inject a
        // huge Poisson source → runaway p_dyn → pressure-gradient force → velocity that
        // pegs the ±100 m/s clamp → larger divergence (a closed dry instability,
        // independent of moist physics). Bounding the per-cell source to denom·p_dyn_cap
        // caps its contribution to p_dyn at p_dyn_cap (discrete max principle). p_dyn is
        // non-dimensionalised by p_0=1013.25 hPa, so healthy dynamic pressure is ≪1
        // non-dim (tens of hPa); the runaway reached 18–43 non-dim (≈18000–44000 hPa).
        // p_dyn_cap is a GENEROUS backstop only. A tight value (1.0) under-removed the
        // velocity divergence (the source IS the divergence the projection must cancel),
        // breaking incompressibility and triggering a worse dry CFL blow-up at the NZ Alps
        // (46°S/170°E, iter 288). 2.0 leaves the projection intact while still catching the
        // gross Andes injection; the orographic Shapiro filter (post-solve, on u/v/w) is now
        // the primary stabiliser for the steep-orography CFL mode, not this clamp.
        // ⚠️ ATM_PDYN_CAP (2026-09-03, default 0 = the shipped 2.0, bit-identical unset). The
        // companion of ATM_PDYN_CEILING below; read that note for why 2.0 is 87 Pa and not
        // 2000 hPa, and why both are a latent barrier rather than a present cause.
        const double p_dyn_cap = [](){
            const char* e = getenv("ATM_PDYN_CAP");
            const double v = e ? atof(e) : 0.0;
            return (v > 0.0) ? v : 2.0;
        }();

        // ⚠️ A/B KNOB 2026-07-21 (ATM_POISSON_METRIC_FIX, default 0 = bit-identical). The θ/φ
        // Poisson Laplacian coefficients (num2/num3/denom) use SINGLE-power inv_rm /
        // inv_rmsinthe, which is inconsistent with this solver's own divergence source and
        // gradient correction (both single-power). The consistent discrete div·grad needs
        // inv_rm² (θ) and inv_rm²/sin²θ (φ) — exactly the pattern the RADIAL term already
        // uses (exp_2_rm in the Laplacian vs single exp_rm in the divergence). The ocean
        // diagnosed this single-power form as the non-idempotent collocated projection and
        // moved its per-iter solve to project_velocity with inv_rm2 / inv_rm2sinthe2
        // (PressureSolverHyd.h:428). This knob brings the atm θ/φ Laplacian to the same
        // consistent metric so it can be A/B tested against the tuned steep-orography
        // stabilisers (p_dyn_cap, p_dyn_ceiling, topo Dirichlet pins) that were calibrated on
        // the old operator. NB: this repairs the metric POWER only; the collocated checkerboard
        // (Rhie-Chow face reconstruction) is a separate, larger port not done here.
        const bool poisson_metric_fix = [](){ const char* e = getenv("ATM_POISSON_METRIC_FIX");
                                              return e ? (atof(e) != 0.0) : false; }();

        // ATOM_METRIC_DIVERGENCE — hoisted out of the cell loop; see lib/Utils.h.
        const bool metric_div = AtomUtils::metric_divergence();

        // ---- ATM_ANELASTIC (ported from ATHAD 2026-08-27, default OFF here) ---------------
        //
        // WHY THIS TREE NEEDS IT. The projection had NO DENSITY IN IT AT ALL -- no dlnrho, no
        // rho_bar -- so div_src was pure div(u) and the solver enforced VOLUME continuity. With
        // u pinned to zero at both walls that forces INT(v dz) = 0. But `Psi` is a MASS
        // streamfunction, INT(rho*v*dz), and rho varies ~4x over the 16 km column, so zero
        // volume flux does not give zero mass flux: Psi(ground) had no reason to close, and did
        // not. Measured over 100 iterations the two integrals move in OPPOSITE directions --
        // with ATM_V_MASSBAL on, volume -21.8 % while mass +51.2 % -- which is why imposing the
        // mass constraint at t = 0 decays.
        //
        // The anelastic form, from ATHAD's copy. The pieces are not separable:
        //   source     D = (1/rho_bar) div(rho_bar u*) = div(u*) + u*_r dln(rho_bar)/dr
        //   operator   (1/rho_bar) div(rho_bar grad p) = lap(p) + dln(rho_bar)/dr dp/dr
        // The operator's extra first-derivative term folds into the EXISTING num_a
        // off-diagonal, which here already carries the metric curvature -- so the coefficient
        // becomes (dlnrho - curv) and no new stencil is introduced.
        //
        // Off-branch bit-identical: dlnrho_i is 0 when the knob is off, so num_a and div_src
        // reduce exactly to what they were.
        static const bool anelastic_knob = [](){
            const char* e = getenv("ATM_ANELASTIC"); return e && atoi(e) != 0; }();
        const bool anelastic = anelastic_knob && ((int)m.m_dlnrho_dr.size() == m.im);
        const double* const dlnrho = anelastic ? m.m_dlnrho_dr.data() : nullptr;
        // ATM_RHIE_CHOW -- fourth-difference pressure smoothing, ported from ATHAD 2026-08-27.
        // 0 (default) restores the old branch exactly. See the term at the update below.
        static const double rc_alpha = [](){
            const char* e = getenv("ATM_RHIE_CHOW"); return e ? atof(e) : 0.0; }();

        // Main compute loop — land mask lookups + hoisted j-invariants + k sliding window
        // ATM_PRESS_SWEEPS -- relaxation sweeps per call, default 1, which is what this solver
        // has always done. Ported from ATHAD, which took it from ATURAN's shared
        // PressureSolver.h (<TAG>_PRESS_SWEEPS). Independent of ATM_PROJ_SWEEPS below, so
        // varying one does not move the other.
        static const int n_sweeps_knob = [](){ const char* e = getenv("ATM_PRESS_SWEEPS");
                                               const int v = e ? atoi(e) : 1;
                                               return v > 0 ? v : 1; }();
        const int n_press_sweeps = (sweeps_override > 0) ? sweeps_override : n_sweeps_knob;

        for (int sweep = 0; sweep < n_press_sweeps; sweep++) {

        // ==================================================================
        // THE RELAXATION IS RED-BLACK, AND WAS NOT ALWAYS. Ported from ATHAD.
        //
        // Each solve is two passes over a checkerboard colouring of (i+j+k): every cell of
        // one colour has all six of its stencil neighbours in the other, so within a pass
        // nothing is read while it is being written.
        //
        // This replaced `#pragma omp parallel for collapse(2) schedule(dynamic, 4)` over
        // (i,j) writing p_dyn IN PLACE while reading p_dyn[i+-1][j+-1] -- over the very two
        // indices the stencil reads across. Cell (i,j,k) was read by the thread owning
        // (i+1,j) or (i,j+1) while its owner was writing it, and schedule(dynamic) made it
        // worse than a thread-count dependence: which thread got which chunk varied with
        // timing, so the SAME binary at the SAME thread count gave different answers run to
        // run. Measured in ATHAD at 20 iterations before the fix:
        //
        //     1 thread                  residuum_atm = 0.75451284
        //     4 threads, run A                       = 0.75455454
        //     4 threads, run B                       = 0.75458447
        //
        // Red-black rather than Jacobi because of what the old loop was reaching for: k runs
        // serially inside a thread, so k-1 is current and k+1 one sweep old, i.e.
        // lexicographic Gauss-Seidel -- correct in serial, broken only by the (i,j)
        // parallelism. Jacobi would have been the easier fix and would have cost the
        // convergence rate. The colour is selected with a `continue` rather than by striding
        // k, because the k loop carries a sliding window over the land mask that assumes
        // consecutive k.
        //
        // This is the family's third encounter with the defect (ATURAN `ffd0e0e` cured the
        // same thing in the shared PressureSolver.h by serialising, which was right there --
        // its computePressure is 0.003 s of a 5.3 s step; here the projection is far too
        // expensive to serialise).
        //
        // NOTE this does NOT reproduce the old 1-thread answer: red-black is a different
        // sweep order from lexicographic Gauss-Seidel, so it converges to the same solution
        // by a different path. Every measurement in this README taken before this commit
        // moves in its last digits.
        // ==================================================================
        for (int colour = 0; colour < 2; colour++) {

        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 1; i < m.im-1; i++) {
            for (int j = 1; j < m.jm-1; j++) {

                // Build geometry struct ONCE per (i,j)
                cAtmosphereModel::CellGeometry geo;

                geo.rm           = m.rad.z[i];
                geo.rm2          = geo.rm * geo.rm;
                geo.exp_rm       = m.metricExpRm(geo.rm);
                geo.curv         = m.metricCurv(geo.rm);
                geo.exp_2_rm     = geo.exp_rm * geo.exp_rm;
                geo.sinthe       = sinthe_table[j];
                geo.sinthe2      = geo.sinthe * geo.sinthe;
                geo.costhe       = cos(m.the.z[j]);
                // ATM_METRIC_RADIUS — must match RungeKutta_Atm_Turb exactly; identity when off.
                const double rmet  = m.metricRadius(geo.rm);
                const double rmet2 = rmet * rmet;
                geo.inv_rm       = 1.0 / rmet;
                geo.inv_rm2      = 1.0 / rmet2;
                geo.inv_rmsinthe         = 1.0 / (rmet * geo.sinthe);
                geo.inv_rm2sinthe        = geo.inv_rm2 / geo.sinthe;
                geo.inv_rm2sinthe2       = geo.inv_rm2 / geo.sinthe2;
                geo.costhe_inv_rm2sinthe = geo.costhe * geo.inv_rm2sinthe;
                geo.inv_2dr   = inv_2dr;
                geo.inv_2dthe = inv_2dthe;
                geo.inv_2dphi = inv_2dphi;
                geo.inv_dr2   = inv_dr2;
                geo.inv_dthe2 = inv_dthe2;
                geo.inv_dphi2 = inv_dphi2;

                // θ/φ Laplacian metric: single-power (legacy) or the consistent double-power
                // (inv_rm2, inv_rm2sinthe2) when ATM_POISSON_METRIC_FIX is set. See knob note above.
                const double m_the = poisson_metric_fix ? geo.inv_rm2        : geo.inv_rm;
                const double m_phi = poisson_metric_fix ? geo.inv_rm2sinthe2 : geo.inv_rmsinthe;
                const double denom = 2.0 * geo.exp_2_rm * inv_dr2
                                   + 2.0 * m_the * inv_dthe2
                                   + 2.0 * m_phi * inv_dphi2;
                const double inv_denom = 1.0 / denom;
                const double num1 = geo.exp_2_rm * inv_dr2;
                // README item 80's second half (ported from ATHAD): the radial Laplacian on a
                // stretched grid is exp_2_rm*(p'' - curv*p'), so the operator needs a
                // FIRST-derivative coefficient of -curv that this stencil never had. curv is 0
                // on the legacy metric, so the added term is exactly +0.0 there and the branch
                // is bit-identical; under ATM_METRIC_EXACT it is the same order as the term it
                // sits beside. Diagonal dominance is unaffected -- the ratio to num1 is
                // curv*dr/2, which is 0.046 at zeta = 3.715.
                const double dlnrho_i = anelastic ? dlnrho[i] : 0.0;
                const double num_a = geo.exp_2_rm * (dlnrho_i - geo.curv) * geo.inv_2dr;
                const double num2 = m_the * inv_dthe2;
                const double num3 = m_phi * inv_dphi2;

                const bool i_in_range = (i < m.im-2);
                const bool j_inner    = (j > 2) && (j < m.jm-2);

                // sliding window for k-direction land/air status
                int8_t lnd_k0 = LAND(i,j,0), lnd_k1 = LAND(i,j,1);

                for (int k = 1; k < m.km-1; k++) {
                    const int8_t lnd_ijk = lnd_k1;
                    const int8_t lnd_kp1 = LAND(i,j,k+1);
                    const int8_t lnd_km1 = lnd_k0;
                    const int8_t lnd_kp2 = (k < m.km-2) ? LAND(i,j,k+2) : 0;
                    const int8_t lnd_km2 = (k > 2)      ? LAND(i,j,k-2) : 0;

                    lnd_k0 = lnd_k1;
                    lnd_k1 = lnd_kp1;

                    // Cells of the other colour are skipped AFTER the window bookkeeping
                    // above -- lnd_k0/lnd_k1 slide with k and assume every k is visited.
                    if (((i + j + k) & 1) != colour) continue;

                    double du_dr, dv_dthe, dw_dphi;
                    bool r_flag   = false;
                    bool the_flag = false;
                    bool phi_flag = false;

                    // EXPERIMENT 2 (2026-06-01): fluid-cell one-sided differences AWAY from
                    // a land neighbour. Without this, the central-difference fallback below
                    // samples aux_u/v/w at the zeroed land cell (line 54), producing a
                    // ~3× overestimate of div_src at coastal fluid cells → spuriously large
                    // negative p_dyn → −∂p/∂r drives u upward → divergence grows → loop.
                    // Use one-sided 2nd-order forward/backward stencils that don't touch
                    // the wall value at all.
                    //
                    // r direction (fluid above land or below overhang)
                    if (!lnd_ijk && i >= 1 && i+2 < m.im) {
                        if (LAND(i-1,j,k)) {
                            du_dr = (-3.0 * m.aux_u.x[i][j][k] + 4.0 * m.aux_u.x[i+1][j][k]
                                     - m.aux_u.x[i+2][j][k]) * inv_2dr;
                            r_flag = true;
                        } else if (i-2 >= 0 && LAND(i+1,j,k)) {
                            du_dr = (3.0 * m.aux_u.x[i][j][k] - 4.0 * m.aux_u.x[i-1][j][k]
                                     + m.aux_u.x[i-2][j][k]) * inv_2dr;
                            r_flag = true;
                        }
                    }

                    // r direction (land cell with air above — existing land-stencil case;
                    // value gets discarded by the Neumann mean update at line 241)
                    if (!r_flag && i_in_range && lnd_ijk && !LAND(i+1,j,k)) {
                        du_dr = (-3.0 * m.aux_u.x[i][j][k] + 4.0 * m.aux_u.x[i+1][j][k]
                                 - m.aux_u.x[i+2][j][k]) * inv_2dr;
                        r_flag = true;
                    }

                    // theta direction
                    if (j_inner) {
                        const int8_t air_jp1 = !LAND(i,j+1,k);
                        const int8_t air_jm1 = !LAND(i,j-1,k);

                        // Fluid cell with land neighbour in theta
                        if (!the_flag && !lnd_ijk) {
                            if (!air_jm1 && air_jp1 && !LAND(i,j+2,k)) {
                                dv_dthe = (-3.0 * m.aux_v.x[i][j][k] + 4.0 * m.aux_v.x[i][j+1][k]
                                           - m.aux_v.x[i][j+2][k]) * inv_2dthe;
                                the_flag = true;
                            } else if (!air_jp1 && air_jm1 && !LAND(i,j-2,k)) {
                                dv_dthe = (3.0 * m.aux_v.x[i][j][k] - 4.0 * m.aux_v.x[i][j-1][k]
                                           + m.aux_v.x[i][j-2][k]) * inv_2dthe;
                                the_flag = true;
                            }
                        }

                        if (lnd_ijk && air_jp1 && !LAND(i,j+2,k)) {
                            dv_dthe = (-3.0 * m.aux_v.x[i][j][k] + 4.0 * m.aux_v.x[i][j+1][k]
                                       - m.aux_v.x[i][j+2][k]) * inv_2dthe;
                            the_flag = true;
                        }
                        if (lnd_ijk && air_jm1 && !LAND(i,j-2,k)) {
                            dv_dthe = -(-3.0 * m.aux_v.x[i][j][k] + 4.0 * m.aux_v.x[i][j-1][k]
                                        - m.aux_v.x[i][j-2][k]) * inv_2dthe;
                            the_flag = true;
                        }
                        if ((lnd_ijk && air_jp1 && LAND(i,j+2,k))
                            || (j == m.jm-2 && !lnd_ijk && LAND(i,j+1,k))) {
                            dv_dthe = (m.aux_v.x[i][j+1][k] - m.aux_v.x[i][j][k]) * inv_dthe;
                            the_flag = true;
                        }
                        if ((lnd_ijk && air_jm1 && LAND(i,j-2,k))
                            || (j == 1 && lnd_ijk && air_jm1)) {
                            dv_dthe = (m.aux_v.x[i][j-1][k] - m.aux_v.x[i][j][k]) * inv_dthe;
                            the_flag = true;
                        }
                    }

                    // phi direction
                    const bool k_inner = (k > 2) && (k < m.km-2);
                    if (k_inner) {
                        const int8_t air_kp1 = !lnd_kp1;
                        const int8_t air_km1 = !lnd_km1;

                        // Fluid cell with land neighbour in phi
                        if (!phi_flag && !lnd_ijk) {
                            if (!air_km1 && air_kp1 && !lnd_kp2) {
                                dw_dphi = (-3.0 * m.aux_w.x[i][j][k] + 4.0 * m.aux_w.x[i][j][k+1]
                                           - m.aux_w.x[i][j][k+2]) * inv_2dphi;
                                phi_flag = true;
                            } else if (!air_kp1 && air_km1 && !lnd_km2) {
                                dw_dphi = (3.0 * m.aux_w.x[i][j][k] - 4.0 * m.aux_w.x[i][j][k-1]
                                           + m.aux_w.x[i][j][k-2]) * inv_2dphi;
                                phi_flag = true;
                            }
                        }

                        if (lnd_ijk && air_kp1 && !lnd_kp2) {
                            dw_dphi = (-3.0*m.aux_w.x[i][j][k] + 4.0*m.aux_w.x[i][j][k+1]
                                       - m.aux_w.x[i][j][k+2]) * inv_2dphi;
                            phi_flag = true;
                        }
                        if (lnd_ijk && air_km1 && !lnd_km2) {
                            dw_dphi = -(-3.0*m.aux_w.x[i][j][k] + 4.0*m.aux_w.x[i][j][k-1]
                                        - m.aux_w.x[i][j][k-2]) * inv_2dphi;
                            phi_flag = true;
                        }
                        if ((lnd_ijk && air_kp1 && lnd_kp2)
                            || (k == m.km-2 && !lnd_ijk && lnd_kp1)) {
                            dw_dphi = (m.aux_w.x[i][j][k+1] - m.aux_w.x[i][j][k]) * inv_dphi;
                            phi_flag = true;
                        }
                        if ((lnd_ijk && air_km1 && lnd_km2)
                            || (k == 1 && lnd_ijk && air_km1)) {
                            dw_dphi = (m.aux_w.x[i][j][k-1] - m.aux_w.x[i][j][k]) * inv_dphi;
                            phi_flag = true;
                        }
                    }

                    // central-difference fallbacks (only for fluid cells fully surrounded by fluid)
                    if (!r_flag)
                        du_dr   = (m.aux_u.x[i+1][j][k] - m.aux_u.x[i-1][j][k]) * inv_2dr;
                    if (!the_flag)
                        dv_dthe = (m.aux_v.x[i][j+1][k] - m.aux_v.x[i][j-1][k]) * inv_2dthe;
                    if (!phi_flag)
                        dw_dphi = (m.aux_w.x[i][j][k+1] - m.aux_w.x[i][j][k-1]) * inv_2dphi;

                    // Solid-wall treatment. A land cell must NOT act as a Poisson pressure
                    // source: the one-sided divergence stencils above sample the adjacent
                    // ocean velocity (e.g. the Drake jet), which injected a spurious high
                    // +p_dyn onto coastal land cells, bled into the ocean pressure and
                    // blocked the through-flow. The first attempt (zero the source ⇒ harmonic
                    // mean of ALL neighbours) over-corrected: it pulled the gap pressure toward
                    // the land average and flattened the venturi gradient, so the Drake flow
                    // stagnated (−u/−w piled up upstream). Instead impose a TRUE Neumann
                    // ∂p/∂n=0 wall — set the land cell to the mean of its FLUID neighbours
                    // only. The fluid then develops its pressure (and the through-gap gradient)
                    // freely; the wall mirrors it, with no spurious source and no cross-coast
                    // force, so the passage flow is preserved. Buried cell (no fluid neighbour) → 0.
                    if (lnd_ijk) {
                        double psum = 0.0; int pn = 0;
                        if (!LAND(i+1,j,k)) { psum += m.p_dyn.x[i+1][j][k]; pn++; }
                        if (!LAND(i-1,j,k)) { psum += m.p_dyn.x[i-1][j][k]; pn++; }
                        if (!LAND(i,j+1,k)) { psum += m.p_dyn.x[i][j+1][k]; pn++; }
                        if (!LAND(i,j-1,k)) { psum += m.p_dyn.x[i][j-1][k]; pn++; }
                        if (!lnd_kp1)       { psum += m.p_dyn.x[i][j][k+1]; pn++; }
                        if (!lnd_km1)       { psum += m.p_dyn.x[i][j][k-1]; pn++; }
                        m.p_dyn.x[i][j][k] = (pn > 0) ? psum / pn : 0.0;
                    } else {
                        // pressure update (fluid interior) — clamp the divergence source
                        // so a steep-orography velocity spike cannot drive p_dyn unbounded.
                        double div_src = du_dr   * geo.exp_rm
                                       + dv_dthe * geo.inv_rm
                                       + dw_dphi * geo.inv_rmsinthe;

                        // ATOM_METRIC_DIVERGENCE (lib/Utils.h) — the spherical divergence also
                        // carries +2u/r and +v*cot(theta)/r. Formed from aux_*, because the source
                        // is the divergence of the PROVISIONAL velocity that the projection has to
                        // remove. Read the warning in lib/Utils.h before switching this on.
                        if (metric_div) {
                            div_src += (2.0 * m.aux_u.x[i][j][k]
                                        + m.aux_v.x[i][j][k] * geo.costhe / geo.sinthe)
                                       * geo.inv_rm;
                        }
                        // Anelastic source term. rho_bar depends on r only, so
                        // (1/rho_bar) div(rho_bar u*) - div(u*) is this one extra piece. Formed
                        // from aux_u, because the source is the divergence of the PROVISIONAL
                        // velocity the projection has to remove.
                        if (anelastic)
                            div_src += m.aux_u.x[i][j][k] * dlnrho[i] * geo.exp_rm;

                        const double src_max = denom * p_dyn_cap;

                        if (!is_finite_safe(div_src))   div_src = 0.0;
                        else if (div_src >  src_max)    div_src =  src_max;
                        else if (div_src < -src_max)    div_src = -src_max;

                        // ---- RHIE-CHOW PRESSURE SMOOTHING (ATM_RHIE_CHOW) --------------
                        //
                        // The operator on the left is the COMPACT 7-point Laplacian at dr, while
                        // div_src and the gradient correction in project_initial_velocity Step 3
                        // are 2*dr central differences. A 2*dr difference annihilates the Nyquist
                        // mode exactly, so a checkerboard in p_dyn is invisible to the correction
                        // and unconstrained by the source. Interpolating velocities to faces does
                        // NOT help -- plain averaging reproduces the 2*dr stencil identically,
                        // which is the whole reason Rhie-Chow exists.
                        //
                        // Add the fourth difference the face reconstruction implies,
                        // D4 = (L_compact - L_wide) p. It annihilates smooth fields, so the
                        // smooth solution is untouched; on the Nyquist mode L_wide = 0. Applied
                        // AFTER the source clamp, so the clamp still bounds the velocity
                        // divergence and never sees this term. Interior only.
                        //
                        // MEASURED IN ATHAD: zonal Nyquist rms 2.42e-3 -> 9.49e-4 at alpha = 1,
                        // then a hard instability edge before 1.5, because a fourth difference
                        // reaches i+-2 -- the SAME COLOUR in a red-black sweep, hence necessarily
                        // lagged. NOT yet measured in THIS tree, where the mode is near-isotropic
                        // (share 0.443/0.331/0.660, absolute Nyquist ~1.6e-06 on all three axes)
                        // rather than confined to k, and where c_phi/c_r is 0.0322 against 0.59.
                        // Do not assume alpha = 1 transfers.
                        double rc = 0.0;
                        if (rc_alpha != 0.0) {
                            if (i >= 2 && i < m.im-2)
                                rc += num1 * ((m.p_dyn.x[i+1][j][k] - 2.0*m.p_dyn.x[i][j][k]
                                             + m.p_dyn.x[i-1][j][k])
                                            - 0.25*(m.p_dyn.x[i+2][j][k] - 2.0*m.p_dyn.x[i][j][k]
                                                  + m.p_dyn.x[i-2][j][k]));
                            if (j >= 2 && j < m.jm-2)
                                rc += num2 * ((m.p_dyn.x[i][j+1][k] - 2.0*m.p_dyn.x[i][j][k]
                                             + m.p_dyn.x[i][j-1][k])
                                            - 0.25*(m.p_dyn.x[i][j+2][k] - 2.0*m.p_dyn.x[i][j][k]
                                                  + m.p_dyn.x[i][j-2][k]));
                            if (k >= 2 && k < m.km-2)
                                rc += num3 * ((m.p_dyn.x[i][j][k+1] - 2.0*m.p_dyn.x[i][j][k]
                                             + m.p_dyn.x[i][j][k-1])
                                            - 0.25*(m.p_dyn.x[i][j][k+2] - 2.0*m.p_dyn.x[i][j][k]
                                                  + m.p_dyn.x[i][j][k-2]));
                            rc *= rc_alpha;
                        }

                        m.p_dyn.x[i][j][k] =
                            ((m.p_dyn.x[i+1][j][k] + m.p_dyn.x[i-1][j][k]) * num1
                           + (m.p_dyn.x[i][j+1][k] + m.p_dyn.x[i][j-1][k]) * num2
                           + (m.p_dyn.x[i][j][k+1] + m.p_dyn.x[i][j][k-1]) * num3
                           + (m.p_dyn.x[i+1][j][k] - m.p_dyn.x[i-1][j][k]) * num_a
                           - div_src - rc) * inv_denom;
                    }
                } // k
            } // j
        } // i

        } // colour

        } // sweep

        #undef LAND

        // Radial boundary extrapolation.
        //
        // At i=0 the old cubic p[0] = p[3] - 3*p[2] + 3*p[1] has condition number ~7
        // and amplifies any cliff-cell pressure into the reference layer by up to 3×.
        // Over many pressure-solve calls this compounds into huge ± spikes at steep
        // topography (e.g. Himalaya). Replaced by:
        //   - i_topography[j][k] >= 1 (i=0 inside the mountain): hold p_dyn = 0, matching
        //     bcSolidGround's treatment of fully buried interior cells; acts as a Dirichlet
        //     pin for the otherwise all-Neumann pressure Poisson.
        //   - i_topography[j][k] == 0 (i=0 is a real ocean surface): von Neumann zero-gradient
        //     ∂p/∂n = 0 — physically correct at a free-slip wall.
        // Tried fully zero-gradient (no Dirichlet anywhere) — Poisson becomes singular and
        // diverges to NaN within 100 iters at the pole.
        // At i=im-1 the column is always air, so the cubic is fine and is kept.
        #pragma omp parallel for collapse(2)
        for (int k = 0; k < m.km; k++) {
            for (int j = 0; j < m.jm; j++) {
                if (m.i_topography[j][k] >= 1) {
                    m.p_dyn.x[0][j][k] = 0.0;
                } else {
                    m.p_dyn.x[0][j][k] = m.c43 * m.p_dyn.x[1][j][k]
                                       - m.c13 * m.p_dyn.x[2][j][k];
                }
                // Top (rigid lid): zero-gradient plain copy. The old cubic
                // p[im-1] = p[im-4] − 3·p[im-3] + 3·p[im-2] is the SAME high-condition (~7,
                // up to 3× amplification) stencil removed at i=0 above, and it has the same
                // failure mode here: over high orography (user-observed Himalaya k=87 / j=62)
                // the strong near-top vertical p_dyn gradient is amplified 3× into a spurious
                // pressure MAXIMUM pegged at the |p_dyn| ceiling (+3 non-dim ≈ 3000 hPa) right
                // at the lid, which then imprints on the orography-oriented velocities. The old
                // "column is air so the cubic is fine" rationale was wrong — the amplification
                // is independent of air/land. At a rigid lid ∂p_dyn/∂r ≈ 0 (no through-lid
                // flow), so a non-amplifying zero-gradient copy is both physical and stable
                // (amplification 1, matching the θ-pole plain-copy BC below). 2026-06-10.
                m.p_dyn.x[m.im-1][j][k] = m.p_dyn.x[m.im-2][j][k];
            }
        }

        // Theta-pole BC: plain copy (zero-gradient, no extrapolation).
        // The 2nd-order Neumann form p[0] = (4/3)·p[1] − (1/3)·p[2] amplifies polar
        // grid noise by 4/3 per call. The Poisson iteration is rerun every moist_stride,
        // and the metric 1/sin²θ already blows up near the pole — compounding the 4/3
        // factor drove p_dyn to NaN at 90°N within ~150 iters, then dpdr/dpdthe/dpdphi
        // propagated NaN into u,v,w in the next RHS call. Plain copy = amplification 1
        // and matches the axisymmetric-pole assumption used for the same fields in bcTheta.
        #pragma omp parallel for collapse(2)
        for (int k = 0; k < m.km; k++) {
            for (int i = 0; i < m.im; i++) {
                m.p_dyn.x[i][0][k]      = m.p_dyn.x[i][1][k];
                m.p_dyn.x[i][m.jm-1][k] = m.p_dyn.x[i][m.jm-2][k];
            }
        }

        // Phi boundary average
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                m.p_dyn.x[i][j][0]      = m.c43 * m.p_dyn.x[i][j][1]      - m.c13 * m.p_dyn.x[i][j][2];
                m.p_dyn.x[i][j][m.km-1] = m.c43 * m.p_dyn.x[i][j][m.km-2] - m.c13 * m.p_dyn.x[i][j][m.km-3];
                m.p_dyn.x[i][j][0] = m.p_dyn.x[i][j][m.km-1]
                    = (m.p_dyn.x[i][j][0] + m.p_dyn.x[i][j][m.km-1]) / 2.0;
            }
        }

        // Hard ceiling on |p_dyn| (applied after every sweep, incl. boundaries). The
        // gradient of p_dyn is the dominant RHS force (presgrad ≫ buoyancy/Coriolis), so
        // bounding p_dyn directly bounds the velocity forcing and breaks the dry steep-
        // orography loop p_dyn → ∇p → w pegs ±100 m/s → divergence → p_dyn. The source
        // clamp above limits injection, but p_dyn still accumulates across the cliff via
        // the Laplacian terms (reached ~7.7 non-dim ≈ 7800 hPa with the source clamp
        // alone); this caps the accumulated result. p_dyn is non-dim'd by p_0=1013.25 hPa
        // p_dyn_ceiling is a NaN/extreme backstop only — NOT the primary control. A tight
        // value (1.0) clipped p_dyn below the level the projection legitimately needs during
        // the violent dry spin-up, leaving residual divergence that blew up at the NZ Alps
        // (iter 288). 10.0 (≈10000 hPa) sits above the accumulated steep-orography value
        // (~7.7 with p_dyn_cap=2) so it never clips normal operation, only the runaway and
        // NaN/Inf (→0).
        //
        // PHASE-DEPENDENT TEST (2026-06-09): the iter-483 runaway is in the VISCOUS/moist
        // phase, where p_dyn pegs −10 over the Tian Shan/Pamir and the −10→0 vertical gradient
        // becomes a −16.8 pgr that drives the velocity (see [[project_upper_velocity_secular_growth]]).
        // Bound the FORCE there by lowering the ceiling to 3.0 once past moist onset (iter 300),
        // while KEEPING 10.0 through the dry spin-up so we don't re-trigger the NZ-Alps failure
        // the tight uniform value caused. Risk: the divergence may legitimately need a larger
        // p_dyn (converging the Poisson gave an even bigger gradient), so a too-tight cap can
        // leave residual velocity divergence that blows up elsewhere — this run tests that.
        // NOTE (2026-06-11): a SMOOTH tanh saturation p=C·tanh(p/C) was tried here in
        // place of this hard clamp (hypothesis: the flat ±ceiling block's discontinuous
        // edge feeds a spurious vertical pgr). FALSIFIED — from atm_restart_500.bin it
        // still NaN'd at the SAME iter 532, merely RELOCATING the seed from i=10/44°N/128°E
        // (near-surface) to i=26/42°N/68°E (Pamir, ~4 km, the documented secular-growth
        // band). Symptom-only / whack-a-mole; the crash time is p_dyn-clip-independent.
        // Reverted to the hard clamp. See [[project_upper_velocity_secular_growth]].
        // ⚠️ ATM_PDYN_CEILING (2026-09-03, default 0 = the shipped phase-dependent value, so an
        // unset environment is BIT-IDENTICAL). Both clamps are LATENT BARRIERS to any repair that
        // gives this model a real thermal pressure field, and the size of the barrier is not the
        // size the comments above claim. p_dyn is non-dimensionalised by rho*u_0^2 ~ 43.7 Pa and
        // NOT by p_0 = 1013.25 hPa -- there is no Euler number in the pressure-gradient term of
        // rhs_v, so the momentum equation is only dimensionally consistent under rho*u_0^2 (see
        // CLAUDE.md, "p_dyn is not in the units this tree has always said it is", verified three
        // independent ways). So this ceiling is 131 Pa, not 3000 hPa, and p_dyn_cap is 87 Pa.
        // A geostrophically balanced mid-latitude pressure field in these units is p_dyn ~ 70 --
        // 23x ABOVE the ceiling. Nothing binds today (max|p_dyn| ~ 0.017, 176x below it), which
        // is exactly why this needs a COUNTER and not an argument: the first repair that works
        // will hit the clamp before it shows a circulation, and would otherwise look like a null.
        const double p_dyn_ceiling = [&](){
            const char* e = getenv("ATM_PDYN_CEILING");
            const double v = e ? atof(e) : 0.0;
            return (v > 0.0) ? v : ((m.total_iter_count > 300) ? 3.0 : 10.0);
        }();
        long   n_clip   = 0;     // cells the ceiling actually truncated
        long   n_nan    = 0;     // cells zeroed as non-finite
        double pmax_pre = 0.0;   // max |p_dyn| BEFORE the ceiling -- the number the ceiling hides
        #pragma omp parallel for collapse(3) reduction(+:n_clip,n_nan) reduction(max:pmax_pre)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                for (int k = 0; k < m.km; k++) {
                    double p = m.p_dyn.x[i][j][k];
                    if (!is_finite_safe(p))    { m.p_dyn.x[i][j][k] = 0.0; n_nan++; continue; }
                    const double a = fabs(p);
                    if (a > pmax_pre) pmax_pre = a;
                    if (p >  p_dyn_ceiling)      { m.p_dyn.x[i][j][k] =  p_dyn_ceiling; n_clip++; }
                    else if (p < -p_dyn_ceiling) { m.p_dyn.x[i][j][k] = -p_dyn_ceiling; n_clip++; }
                }
            }
        }
        auto end = std::chrono::high_resolution_clock::now();
        if (verbose) {
            reportDivergence("at pressure solve");
            // Print-only. No new class member -- deliberately, because adding one changes
            // sizeof() and this tree's build hazard is an ODR size mismatch that AddressSanitizer
            // does not see (see CLAUDE.md, "The build hazard").
            {
                const long ncell = (long)m.im * m.jm * m.km;
                cout << "      ATOM: p_dyn clamp   ceiling = " << fixed << setprecision(3)
                     << p_dyn_ceiling
                     << "   max|p_dyn| pre-clip = " << scientific << setprecision(3) << pmax_pre
                     << "   clipped = " << n_clip << " cells (" << fixed << setprecision(4)
                     << (ncell > 0 ? 100.0 * n_clip / ncell : 0.0) << " %)"
                     << "   non-finite zeroed = " << n_nan
                     << (n_clip > 0 ? "   *** CEILING BINDING ***" : "") << endl;
            }
            auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
            printf(" time measured: %.3f seconds for PressureSolverAtm\n", elapsed.count() * 1e-9);
            cout << "      ATOM: PressureSolverAtm ended" << endl;
        }
    }

    // One-shot Helmholtz projection of the prescribed initial velocity onto its
    // divergence-free subspace.  Without this, the incremental projection that runs
    // every moist_stride iters spends the first ~100 iters silently destroying the
    // dilatational component of the Hadley/Ferrel profile, taking the prescribed
    // global circulation with it (max u decays from 30 m/s to <0.5 m/s by iter 150).
    // Performing the projection once at t=0 lets the simulation start from a clean
    // divergence-free state; the surviving solenoidal part is preserved and only the
    // unphysical dilatational artefact of the analytical profile is removed.
    //
    // The standard solver loop assembles ∇²p = ∇·v* where v* is the post-RHS
    // intermediate velocity; here we feed v itself as the source by copying v into
    // aux_u/aux_v/aux_w. After n_sweeps Jacobi passes, p_dyn approximates the
    // projection pressure; we then apply v ← v − ∇p in the same metric form used
    // by the time-stepping RHS, and reset p_dyn to 0 so the next RK4 call does not
    // double-correct via its own −∂p/∂r term.
    // ATM_PROJ_SWEEPS -- relaxation sweeps PER PASS of the initial projection.
    //
    // DEFAULT 1 HERE, WHICH IS UNMEASURED IN THIS TREE. ATHAD ships 10 (README item 68),
    // where one sweep left the meridional streamfunction not closing at the ground -- Psi
    // there must be zero, and at one sweep it was 2.09x the interior circulation, falling
    // 52.5 % at 10 sweeps and 55.9 % at 100, i.e. a knee rather than a converged value. That
    // measurement has not been made here, so the default stays at the historical 1 and the
    // knob exists to make it. The cost is a one-time startup expense: the projection is 200
    // passes and a relaxation sweep is ~0.1 % of a time step.
    //
    // It is deliberately separate from ATM_PRESS_SWEEPS. This routine calls run() 200 times,
    // so a single knob would have made "10 sweeps in the time loop" mean 2000 relaxations at
    // startup as well, and the arms of any comparison would have differed in their INITIAL
    // STATE as well as in the quantity under test. ATHAD lost an attribution that way.

    // ==================================================================
    // div(u), PRINTED AT THE POINTS WHERE IT MEANS SOMETHING. Ported from ATHAD
    // (PressureSolverAtm.h), which prints it every iteration and judges the projection by it;
    // this tree had no divergence diagnostic at all, so the closure question had to be
    // answered offline from the streamfunction CSV -- which is circular, because Psi does not
    // close at the ground either.
    //
    // READ THE CALL SITES, NOT JUST THE NUMBER. In the TIME LOOP `run()` computes p_dyn and
    // NOTHING ELSE: the velocity is never explicitly projected there, it feels the pressure
    // through the -dp/dr term of the next RK4 stage, and `pressure_stride = 4` means even that
    // happens on one iteration in four. The only place a velocity correction v <- v - grad(p)
    // is actually applied is project_initial_velocity, Step 3. So:
    //
    //   "after initial projection"  tests whether the projection ACHIEVES div(u) = 0;
    //   "at pressure solve"         is the working level of divergence during the run, after
    //                               RK4 and the polar / orographic / radial filters have had
    //                               it -- a different question, and not the solver's fault.
    //
    // The divergence is formed in the model's OWN metric, term for term as the Poisson source
    // is, so a small number here means the solver hit its own target rather than a target
    // rewritten by the diagnostic. Print-only.
    // ==================================================================
    void reportDivergence(const char* tag) const {
        using namespace std;
        const double inv_2dr   = 1.0 / (2.0 * m.dr);
        const double inv_2dthe = 1.0 / (2.0 * m.dthe);
        const double inv_2dphi = 1.0 / (2.0 * m.dphi);
        const bool metric_div  = AtomUtils::metric_divergence();

        // div(rho u)/rho ALONGSIDE div(u). ATHAD prints this every iteration and judges the
        // projection by it; this tree had only the volume form, and CLAUDE.md's own note says
        // the metric question "is NOT DECIDABLE BY MORE A/B RUNS -- the prerequisite is the
        // div(rho u)/rho print, not another arm". The two differ because rho varies ~4x over
        // the 16 km column, so a field can be volume-divergence-free and carry a net MASS flux,
        // which is exactly the Psi(ground) non-closure. Term for term as the anelastic Poisson
        // source forms it (line 444, `div_src += aux_u * dlnrho[i] * exp_rm`), so the number is
        // the model's own quantity and not one rewritten by the diagnostic. The base state is
        // built in densities() and exists whether or not ATM_ANELASTIC is set -- so this prints
        // on BOTH branches, which is the point: it measures what the shipped solver is leaving
        // behind. Zero extra sweep; folded into the existing one. Print-only.
        const bool have_rho = ((int)m.m_dlnrho_dr.size() == m.im);
        const double* const dlnrho_d = have_rho ? m.m_dlnrho_dr.data() : nullptr;

        double d2 = 0.0, dmax = 0.0, rad2 = 0.0, dm2 = 0.0, dmmax = 0.0;
        long   n  = 0;

        #pragma omp parallel for collapse(2) schedule(static) \
                reduction(+:d2,rad2,dm2,n) reduction(max:dmax,dmmax)
        for (int i = 1; i < m.im-1; i++) {
            for (int j = 1; j < m.jm-1; j++) {
                const double rm       = m.rad.z[i];
                const double exp_rm   = m.metricExpRm(rm);
                const double rmet     = m.metricRadius(rm);
                double sinthe = sin(m.the.z[j]);
                if (sinthe < 0.55) sinthe = 0.55;          // the model's own metric floor
                const double inv_rm   = 1.0 / rmet;
                const double inv_rms  = 1.0 / (rmet * sinthe);
                const double cotanthe = cos(m.the.z[j]) / sinthe;

                for (int k = 1; k < m.km-1; k++) {
                    // fluid cells with fluid neighbours only -- a one-sided difference across
                    // a coast is a different operator and would be counted as divergence
                    if (m.h.x[i][j][k] == 1.0 || m.h.x[i+1][j][k] == 1.0 || m.h.x[i-1][j][k] == 1.0
                     || m.h.x[i][j+1][k] == 1.0 || m.h.x[i][j-1][k] == 1.0
                     || m.h.x[i][j][k+1] == 1.0 || m.h.x[i][j][k-1] == 1.0) continue;
                    const double d_r = (m.u.x[i+1][j][k] - m.u.x[i-1][j][k]) * inv_2dr * exp_rm;
                    double d = d_r
                             + (m.v.x[i][j+1][k] - m.v.x[i][j-1][k]) * inv_2dthe * inv_rm
                             + (m.w.x[i][j][k+1] - m.w.x[i][j][k-1]) * inv_2dphi * inv_rms;
                    if (metric_div)
                        d += (2.0 * m.u.x[i][j][k] + m.v.x[i][j][k] * cotanthe) * inv_rm;
                    if (!is_finite_safe(d)) continue;
                    d2 += d * d; rad2 += d_r * d_r;
                    if (fabs(d) > dmax) dmax = fabs(d);
                    const double dm = have_rho
                        ? d + m.u.x[i][j][k] * dlnrho_d[i] * exp_rm : d;
                    if (is_finite_safe(dm)) {
                        dm2 += dm * dm;
                        if (fabs(dm) > dmmax) dmmax = fabs(dm);
                    }
                    n++;
                }
            }
        }
        // ---- CHECKERBOARD INDEX ON p_dyn -------------------------------------------
        //
        // Roger reports diagonal stripes in p_dyn. Diagonal is the signature to take
        // seriously: the red-black colouring below is on (i+j+k)&1, so on ANY single slice
        // the two colours lie along diagonals -- a field that differs between colours draws
        // exactly diagonal stripes and nothing else.
        //
        // Two mechanisms can put a difference there and they need separating, so this
        // measures rather than assumes:
        //
        //  (1) ITERATION. Red-black Gauss-Seidel updates each colour from the other; stop
        //      mid-cycle, or stop far from convergence, and the colours sit at different
        //      iteration levels. Cured by sweeping more.
        //  (2) ODD-EVEN DECOUPLING, i.e. the collocated-grid checkerboard. The Poisson
        //      SOURCE and the gradient correction in project_initial_velocity Step 3 are
        //      2*dr central differences, while the operator inverted here is the COMPACT
        //      7-point Laplacian at dr. A 2*dr central difference annihilates the Nyquist
        //      mode exactly, so a checkerboard in p_dyn is invisible to the velocity
        //      correction and unconstrained by the source -- the classical reason Rhie-Chow
        //      exists, and this file already records div and grad as non-adjoint. NOT cured
        //      by sweeping more.
        //
        // The index: chk = p - (mean of the 6 neighbours). Smooth field -> O(dr^2 * lap p),
        // small next to rms(p). Pure checkerboard -> chk = 2p, so the ratio tends to 2.
        // Printed with the red/black mean split, which is the same thing seen globally.
        double c2 = 0.0, pp2 = 0.0, sred = 0.0, sblk = 0.0;
        long   nc = 0, nred = 0, nblk = 0;
        #pragma omp parallel for collapse(2) schedule(static) \
                reduction(+:c2,pp2,nc,sred,sblk,nred,nblk)
        for (int i = 1; i < m.im-1; i++) {
            for (int j = 1; j < m.jm-1; j++) {
                for (int k = 1; k < m.km-1; k++) {
                    if (m.h.x[i][j][k] == 1.0 || m.h.x[i+1][j][k] == 1.0 || m.h.x[i-1][j][k] == 1.0
                     || m.h.x[i][j+1][k] == 1.0 || m.h.x[i][j-1][k] == 1.0
                     || m.h.x[i][j][k+1] == 1.0 || m.h.x[i][j][k-1] == 1.0) continue;
                    const double pc = m.p_dyn.x[i][j][k];
                    if (!is_finite_safe(pc)) continue;
                    const double nb = (m.p_dyn.x[i+1][j][k] + m.p_dyn.x[i-1][j][k]
                                     + m.p_dyn.x[i][j+1][k] + m.p_dyn.x[i][j-1][k]
                                     + m.p_dyn.x[i][j][k+1] + m.p_dyn.x[i][j][k-1]) / 6.0;
                    const double chk = pc - nb;
                    c2 += chk * chk; pp2 += pc * pc; nc++;
                    if (((i + j + k) & 1) == 0) { sred += pc; nred++; }
                    else                        { sblk += pc; nblk++; }
                }
            }
        }
        // ---- DIRECTIONAL NYQUIST SHARE (ported from ATHAD 2026-08-27) ------------------
        //
        // THE GLOBAL INDEX ABOVE IS A FALSE NULL, and it is kept beside this one rather than
        // replaced so the two can be read together. It divides by rms(p_dyn). In ATHAD that
        // denominator is the large smooth radial/latitudinal structure, so a checkerboard at
        // 1e-6 of it reads 0.0003 and looks clean -- while the SAME field, measured per
        // direction, is 96 % Nyquist in k. This tree reads 0.63 not because its checkerboard is
        // worse but because its rms(p_dyn) has collapsed to ~5e-06; in ABSOLUTE terms ATHAD's
        // grid-scale component is ~15 000x larger. A ratio whose denominator moves between
        // trees cannot be compared between them.
        //
        // So: the 2-delta component measured against the field's own variation ALONG EACH AXIS,
        // and the ABSOLUTES beside the ratio -- because a ratio alone cannot tell a cure from a
        // uniform shrink of everything in that direction. In ATHAD that distinction mattered:
        // the share moved only 0.961 -> 0.885 while the absolute Nyquist fell 2.55x.
        //
        // Land is skipped the same way the divergence above skips it: fluid cells whose
        // neighbours along the axis being differenced are also fluid.
        double share[3] = {0.0, 0.0, 0.0};
        double anom_rms[3] = {0.0, 0.0, 0.0}, nyq_rms[3] = {0.0, 0.0, 0.0};
        for (int dir = 0; dir < 3; dir++) {
            const int ni = (dir == 0) ? m.im : ((dir == 1) ? m.jm : m.km);
            const int n1 = (dir == 0) ? m.jm : m.im;
            const int n2 = (dir == 2) ? m.jm : m.km;
            if (ni < 5) continue;
            double an2 = 0.0, ny2 = 0.0; long ncd = 0;
            #pragma omp parallel for collapse(2) schedule(static) reduction(+:an2,ny2,ncd)
            for (int p1 = 1; p1 < n1 - 1; p1++) {
                for (int p2 = 1; p2 < n2 - 1; p2++) {
                    double sum = 0.0; long cnt = 0;
                    for (int q = 1; q < ni - 1; q++) {
                        const int i = (dir == 0) ? q : p1;
                        const int j = (dir == 1) ? q : ((dir == 0) ? p1 : p2);
                        const int k = (dir == 2) ? q : p2;
                        if (m.h.x[i][j][k] == 1.0) continue;
                        const double v = m.p_dyn.x[i][j][k];
                        if (is_finite_safe(v)) { sum += v; cnt++; }
                    }
                    if (cnt < 3) continue;
                    const double mean = sum / (double)cnt;
                    for (int q = 1; q < ni - 1; q++) {
                        const int i = (dir == 0) ? q : p1;
                        const int j = (dir == 1) ? q : ((dir == 0) ? p1 : p2);
                        const int k = (dir == 2) ? q : p2;
                        const int ip = (dir == 0) ? i+1 : i, im1 = (dir == 0) ? i-1 : i;
                        const int jp = (dir == 1) ? j+1 : j, jm1 = (dir == 1) ? j-1 : j;
                        const int kp = (dir == 2) ? k+1 : k, km1 = (dir == 2) ? k-1 : k;
                        if (m.h.x[i][j][k] == 1.0 || m.h.x[ip][jp][kp] == 1.0
                                                  || m.h.x[im1][jm1][km1] == 1.0) continue;
                        const double v  = m.p_dyn.x[i][j][k];
                        const double vp = m.p_dyn.x[ip][jp][kp];
                        const double vm = m.p_dyn.x[im1][jm1][km1];
                        if (!is_finite_safe(v) || !is_finite_safe(vp) || !is_finite_safe(vm)) continue;
                        const double anom = v - mean;
                        const double nyq  = 0.5 * (v - 0.5 * (vp + vm));
                        an2 += anom * anom; ny2 += nyq * nyq; ncd++;
                    }
                }
            }
            share[dir]    = (an2 > 0.0) ? sqrt(ny2 / an2) : 0.0;
            anom_rms[dir] = (ncd > 0)   ? sqrt(an2 / (double)ncd) : 0.0;
            nyq_rms[dir]  = (ncd > 0)   ? sqrt(ny2 / (double)ncd) : 0.0;
        }

        const double rms_chk = (nc > 0) ? sqrt(c2 / nc)  : 0.0;
        const double rms_p   = (nc > 0) ? sqrt(pp2 / nc) : 0.0;
        const double mred    = (nred > 0) ? sred / nred : 0.0;
        const double mblk    = (nblk > 0) ? sblk / nblk : 0.0;

        const double rms  = (n > 0) ? sqrt(d2 / n)   : 0.0;
        const double rrad = (n > 0) ? sqrt(rad2 / n) : 0.0;
        const ios::fmtflags f = cout.flags();
        const streamsize    p = cout.precision();
        cout << "      ATOM: div(u) " << setw(24) << left << tag << right
             << " rms = " << scientific << setprecision(3) << rms
             << "   max = " << dmax
             << "   radial term = " << rrad
             << fixed << setprecision(3)
             << "   ratio rms/radial = " << ((rrad > 0.0) ? rms / rrad : 0.0)
             << "   (" << n << " cells)" << endl;
        {
            const double rmsm = (n > 0) ? sqrt(dm2 / n) : 0.0;
            cout << "      ATOM: div(rho u)/rho " << setw(17) << left << tag << right
                 << " rms = " << scientific << setprecision(3) << rmsm
                 << "   max = " << dmmax
                 << fixed << setprecision(3)
                 << "   mass/volume = " << ((rms > 0.0) ? rmsm / rms : 0.0)
                 << (have_rho ? "" : "   [NO BASE STATE -- equals div(u)]") << endl;
        }
        cout << "      ATOM: p_dyn checkerboard " << setw(16) << left << tag << right
             << " index = " << fixed << setprecision(4)
             << ((rms_p > 0.0) ? rms_chk / rms_p : 0.0)
             << "   (0 = smooth, 2 = pure)   rms p_dyn = "
             << scientific << setprecision(3) << rms_p
             << "   red-black mean split = " << (mred - mblk) << endl;
        // Operator weights at a mid-column cell. Printed because the direction the Nyquist mode
        // lives in should be the direction the operator constrains LEAST, and that is checkable
        // rather than assumable. In ATHAD it was NOT the explanation: c_phi/c_r = 0.59 there,
        // comparable, so the mode lives in k because k carries no physical signal, not because
        // k is weakly damped.
        {
            const int i = m.im/2, j = m.jm/2;
            const double rm = m.rad.z[i], e = m.metricExpRm(rm);
            double sth = sin(m.the.z[j]); if (sth < 0.55) sth = 0.55;
            const double rmet = m.metricRadius(rm);
            const double c1 = e*e/(m.dr*m.dr);
            const double c2d = (1.0/rmet)/(m.dthe*m.dthe);
            const double c3 = (1.0/(rmet*sth))/(m.dphi*m.dphi);
            cout << "      ATOM: Poisson weights at i=" << i << ", j=" << j
                 << "   c_r = " << scientific << setprecision(3) << c1
                 << "   c_the = " << c2d << "   c_phi = " << c3
                 << fixed << setprecision(4) << "   c_phi/c_r = " << (c1>0.0 ? c3/c1 : 0.0) << endl;
        }
        static const char* dname[3] = {"radial(i)", "merid(j)", "zonal(k)"};
        cout << "      ATOM: p_dyn Nyquist share of the anomaly along each axis: ";
        for (int d = 0; d < 3; d++)
            cout << "  " << dname[d] << " = " << fixed << setprecision(3) << share[d];
        cout << endl;
        cout << "      ATOM: p_dyn absolute rms  ";
        for (int d = 0; d < 3; d++)
            cout << "  " << dname[d] << ": anom = " << scientific << setprecision(3)
                 << anom_rms[d] << " nyq = " << nyq_rms[d];
        cout << endl;
        cout.flags(f); cout.precision(p);
    }

    // ==================================================================
    // ATM_PROJECT_IN_LOOP=<sweeps> -- a real velocity projection inside the time loop.
    // Default 0 = off.
    //
    // WHY THIS IS NOT JUST "CALL run() AND SUBTRACT THE GRADIENT". In the time loop `aux_*` does
    // NOT hold a velocity: RHS_Atm_Turb.cpp:1140 sets `aux_u = rhs_u + dpdr_exp`, a momentum
    // TENDENCY. So `run()` there computes the pressure that makes the ACCELERATION
    // divergence-free, which is a legitimate fractional-step variant -- if div(u) starts at zero
    // and every tendency is divergence-free, div(u) stays zero -- but it means the p_dyn sitting
    // in the array is not the pressure that projects the VELOCITY, and subtracting its gradient
    // from u would be wrong. That is the trap this routine exists to avoid.
    //
    // What it does instead is the same three steps project_initial_velocity uses, on the actual
    // velocity: seed aux from u/v/w, relax the Poisson, then apply v <- v - grad(p) with the
    // metric factors the RHS uses. p_dyn is SAVED and RESTORED around the call, because the
    // time loop's own p_dyn is read by the next RK4 stage through -dp/dr and Step 4 of the
    // initial projection would otherwise wipe it.
    //
    // The question it exists to answer: Psi(ground) does not close, and ATM_ANELASTIC measured a
    // null because changing which continuity the PRESSURE solves for cannot change the VELOCITY
    // when the velocity is never projected. This is the arm that can actually test that.
    // ==================================================================
    static int projectInLoopSweeps(){
        static const int v = [](){
            const char* e = getenv("ATM_PROJECT_IN_LOOP");
            const int n = e ? atoi(e) : 0;
            return n > 0 ? n : 0; }();
        return v;
    }

    void project_velocity_in_loop(int n_sweeps)
    {
        if (n_sweeps <= 0) return;

        // Save p_dyn: the time loop's pressure is live and the projection overwrites it.
        std::vector<double> p_save((size_t)m.im * m.jm * m.km);
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < m.im; i++)
            for (int j = 0; j < m.jm; j++)
                for (int k = 0; k < m.km; k++)
                    p_save[((size_t)i*m.jm + j)*m.km + k] = m.p_dyn.x[i][j][k];

        // Steps 1-3, on the velocity.
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < m.im; i++)
            for (int j = 0; j < m.jm; j++)
                for (int k = 0; k < m.km; k++) {
                    m.aux_u.x[i][j][k] = m.u.x[i][j][k];
                    m.aux_v.x[i][j][k] = m.v.x[i][j][k];
                    m.aux_w.x[i][j][k] = m.w.x[i][j][k];
                    m.p_dyn.x[i][j][k] = 0.0;
                }

        for (int s = 0; s < n_sweeps; s++) run(false, 1);

        // div(u) EITHER SIDE OF THE CORRECTION, which is the one instrument this arm was
        // missing. CLAUDE.md: "`div(u)` is only ever reported AT THE PRESSURE SOLVE, never
        // immediately after `project_velocity_in_loop`. So 'the projection zeroes `div(u)` and
        // one time step puts all of it back' and 'the projection never reduces `div(u)` at all'
        // have the IDENTICAL signature in this table, and they are different defects."
        // Reporting before AND after separates them in one call: if the pair is (large, small)
        // the projection works and the time step undoes it; if it is (large, large) the
        // projection never acts. Only fires when the knob is on, since n_sweeps <= 0 returned
        // above. Print-only.
        reportDivergence("in-loop projection, pre");
        applyGradientCorrection();
        reportDivergence("in-loop projection, post");

        // Restore the time loop's pressure.
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < m.im; i++)
            for (int j = 0; j < m.jm; j++)
                for (int k = 0; k < m.km; k++)
                    m.p_dyn.x[i][j][k] = p_save[((size_t)i*m.jm + j)*m.km + k];
    }

    // v <- v - grad(p_dyn) in the interior, with the metric factors the rhs_u/v/w
    // pressure-gradient terms use. Extracted from project_initial_velocity Step 3 so the
    // in-loop projection cannot drift from it.
    void applyGradientCorrection()
    {
        const double inv_2dr   = 1.0 / (2.0 * m.dr);
        const double inv_2dthe = 1.0 / (2.0 * m.dthe);
        const double inv_2dphi = 1.0 / (2.0 * m.dphi);

        std::vector<double> sinthe_tab(m.jm);
        for (int j = 0; j < m.jm; j++) {
            sinthe_tab[j] = sin(m.the.z[j]);
            if (sinthe_tab[j] < 0.55) sinthe_tab[j] = 0.55;
        }

        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 1; i < m.im-1; i++) {
            for (int j = 1; j < m.jm-1; j++) {
                const double rm           = m.rad.z[i];
                const double exp_rm       = m.metricExpRm(rm);
                const double rmet         = m.metricRadius(rm);
                const double inv_rm       = 1.0 / rmet;
                const double inv_rmsinthe = 1.0 / (rmet * sinthe_tab[j]);
                for (int k = 1; k < m.km-1; k++) {
                    if (m.h.x[i][j][k] == 1.0) continue;          // solid cell
                    const double dpdr   = (m.p_dyn.x[i+1][j][k] - m.p_dyn.x[i-1][j][k]) * inv_2dr;
                    const double dpdthe = (m.p_dyn.x[i][j+1][k] - m.p_dyn.x[i][j-1][k]) * inv_2dthe;
                    const double dpdphi = (m.p_dyn.x[i][j][k+1] - m.p_dyn.x[i][j][k-1]) * inv_2dphi;
                    m.u.x[i][j][k] -= dpdr   * exp_rm;
                    m.v.x[i][j][k] -= dpdthe * inv_rm;
                    m.w.x[i][j][k] -= dpdphi * inv_rmsinthe;
                }
            }
        }
    }

    void project_initial_velocity(int n_sweeps = 200)
    {
        static const int proj_sweeps = [](){
            const char* e = getenv("ATM_PROJ_SWEEPS");
            const int v = e ? atoi(e) : 1;
            return v > 0 ? v : 1; }();
        using namespace std;
        cout << endl << endl << "      ATOM: project_initial_velocity ("
             << n_sweeps << " Jacobi sweeps)" << endl;
        auto t0 = std::chrono::high_resolution_clock::now();

        // Step 1 — seed Poisson source from current velocity, zero p_dyn.
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                for (int k = 0; k < m.km; k++) {
                    m.aux_u.x[i][j][k] = m.u.x[i][j][k];
                    m.aux_v.x[i][j][k] = m.v.x[i][j][k];
                    m.aux_w.x[i][j][k] = m.w.x[i][j][k];
                    m.p_dyn.x[i][j][k] = 0.0;
                }
            }
        }

        // Step 2 — quiet Jacobi iteration until p_dyn approximates the projection
        // pressure.  Each call to run() also re-applies the i, theta, and phi BCs on
        // p_dyn, so polar/topographic anchors stay consistent with the time loop.
        for (int s = 0; s < n_sweeps; s++) {
            run(false, proj_sweeps);
        }

        // Step 3 — gradient correction v ← v − ∇p_dyn in the interior.
        // Metric factors match the rhs_u/v/w pressure-gradient term so the magnitudes
        // are consistent with the rest of the code.  Boundary cells (i, theta, phi
        // outer faces) are left untouched; bcRadius / bcTheta / bcPhi will re-impose
        // their patterns at the next call.
        const double inv_2dr   = 1.0 / (2.0 * m.dr);
        const double inv_2dthe = 1.0 / (2.0 * m.dthe);
        const double inv_2dphi = 1.0 / (2.0 * m.dphi);

        std::vector<double> sinthe_tab(m.jm);
        for (int j = 0; j < m.jm; j++) {
            sinthe_tab[j] = sin(m.the.z[j]);
            if (sinthe_tab[j] < 0.55) sinthe_tab[j] = 0.55;   // metric floor ~57° — keep in sync with the RHS geometry floor (RungeKutta_Atm*.cpp)
        }

        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 1; i < m.im-1; i++) {
            for (int j = 1; j < m.jm-1; j++) {
                const double rm           = m.rad.z[i];
                const double exp_rm       = m.metricExpRm(rm);   // grid coordinate, not the radius
                // The gradient correction must use the SAME metric the source and the RHS use,
                // or the projection stops being a projection. ATM_METRIC_RADIUS; identity when off.
                const double rmet         = m.metricRadius(rm);
                const double inv_rm       = 1.0 / rmet;
                const double inv_rmsinthe = 1.0 / (rmet * sinthe_tab[j]);

                for (int k = 1; k < m.km-1; k++) {
                    const double dpdr   = (m.p_dyn.x[i+1][j][k] - m.p_dyn.x[i-1][j][k]) * inv_2dr;
                    const double dpdthe = (m.p_dyn.x[i][j+1][k] - m.p_dyn.x[i][j-1][k]) * inv_2dthe;
                    const double dpdphi = (m.p_dyn.x[i][j][k+1] - m.p_dyn.x[i][j][k-1]) * inv_2dphi;

                    m.u.x[i][j][k] -= dpdr   * exp_rm;
                    m.v.x[i][j][k] -= dpdthe * inv_rm;
                    m.w.x[i][j][k] -= dpdphi * inv_rmsinthe;
                }
            }
        }

        // THE test of the projection: this is the one place a velocity correction is applied.
        reportDivergence("after initial projection");

        // Step 4 — clear p_dyn and aux so the time loop starts fresh.
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                for (int k = 0; k < m.km; k++) {
                    m.p_dyn.x[i][j][k] = 0.0;
                    m.aux_u.x[i][j][k] = 0.0;
                    m.aux_v.x[i][j][k] = 0.0;
                    m.aux_w.x[i][j][k] = 0.0;
                }
            }
        }

        auto t1 = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0);
        printf("      ATOM: project_initial_velocity ended (%.3fs)\n", elapsed.count() * 1e-9);
    }

private:
    cAtmosphereModel& m;
};
