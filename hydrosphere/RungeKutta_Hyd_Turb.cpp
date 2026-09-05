/*
 * Ocean General Circulation Model (OGCM) applied to turbulent flow
 * 4th-order Runge-Kutta for 7 prognostic fields: t, u, v, w, c, tke, dis
*/

#include "cHydrosphereModel.h"
#include "HydHorizViscosity.h"

#include <cstdint>
#include <cstring>

using namespace std;


void cHydrosphereModel::solveRungeKutta_Hydrosphere_Turb(){
    cout << endl << "      solveRungeKutta_Hydrosphere_Turb begin" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    const double half_dt  = 0.5 * dt;
    const double dt_sixth = dt / 6.0;

    // TKE cap.  Mirrors the atmosphere RK4 — after removing the conservative-form
    // divergence correction the dominant tke growth path is still the
    // −0.667·tke·∇·v term hidden inside `prod`, which is positive and linear in tke
    // in convergent regions.  Until the production assembly is re-derived in
    // strain-rate-only form, hard-clamp |k| to a physically generous bound at every
    // RK4 stage so the simulation cannot diverge.  1000 m²/s² is well above any
    // realistic ocean turbulence (open-ocean ε·L ≪ 1 m²/s² in the mixed layer).
    const double tke_max_phys = 1000.0;                                 // [m²/s²]
    const double tke_max_nd   = tke_max_phys / (u_0 * u_0);

    // NaN-safe clamp.  std::clamp uses unguarded `<` comparisons; under -ffast-math
    // (-ffinite-math-only) the compiler assumes operands are finite and the result
    // on NaN input is undefined.  Detect non-finite values by IEEE-754 exponent bits
    // and replace with lo (no ffast-math hazard).  Same trick used in
    // Paraview_Hyd.cpp::safe_val and lib/Utils.h::is_finite_safe.
    auto safe_clamp = [](double v, double lo, double hi) -> double {
        std::uint64_t bits;
        std::memcpy(&bits, &v, sizeof(bits));
        if ((bits & 0x7FF0000000000000ULL) == 0x7FF0000000000000ULL) return lo;
        return (v < lo) ? lo : (v > hi) ? hi : v;
    };

    const double inv_2dr   = 1.0 / (2.0 * dr);
    const double inv_2dthe = 1.0 / (2.0 * dthe);
    const double inv_2dphi = 1.0 / (2.0 * dphi);
    const double inv_dr2   = 1.0 / (dr * dr);
    const double inv_dthe2 = 1.0 / (dthe * dthe);
    const double inv_dphi2 = 1.0 / (dphi * dphi);

    std::vector<double> sinthe_tbl(jm), costhe_tbl(jm);
    for (int j = 0; j < jm; j++) {
        sinthe_tbl[j] = sin(the.z[j]);
        if (sinthe_tbl[j] < 0.4) sinthe_tbl[j] = 0.4;
        costhe_tbl[j] = cos(the.z[j]);
    }

    // ==================================================================================
    // HYD_A_H_BIHARM -- INNER PASS of the horizontal biharmonic viscosity.
    //
    // Fills H(vn,wn), the angular part of the horizontal vector Laplacian, so the RHS can
    // apply the SAME operator a second time. See HydHorizViscosity.h for why the buffer is
    // there and why it is built on the time-level-n fields.
    //
    // The operator is the 2-D HORIZONTAL vector Laplacian:
    //   H_v = d2v/dthe2 + cot*dv/dthe - (1+cot^2)*v + d2v/dphi2/sin^2 - 2*cos/sin^2*dw/dphi
    //   H_w = d2w/dthe2 + cot*dw/dthe - (1+cot^2)*w + d2w/dphi2/sin^2 + 2*cos/sin^2*dv/dphi
    // It deliberately OMITS the u-coupling terms (2*dudthe, 2*dudphi/sin) that the shipped
    // diffusion_v/diffusion_w carry, because a biharmonic has to be the same operator twice
    // and the outer application has no lap_u to couple to. So this is the horizontal
    // vector Laplacian an ocean model's biharmonic uses, not the full 3-D one -- the one
    // place in this pair of knobs where HYD_A_H (which matches the shipped term exactly)
    // and HYD_A_H_BIHARM differ in shape.
    //
    // Computed only where the full 5-point horizontal stencil is water, and `ok` records
    // that, so the OUTER pass can require its own five stencil points to be marked. The
    // biharmonic therefore switches itself off within two cells of any coast rather than
    // differencing against a buffer that is zero because it was never computed -- which
    // would manufacture exactly the coastal gradient the term exists to remove.
    // ==================================================================================
    if (HydHorizVisc::biharm_strength() != 0.0) {
        HydHorizVisc::s_jm = jm;
        HydHorizVisc::s_km = km;
        const std::size_t n_cells = static_cast<std::size_t>(im) * jm * km;
        HydHorizVisc::lap_v.assign(n_cells, 0.0);
        HydHorizVisc::lap_w.assign(n_cells, 0.0);
        HydHorizVisc::ok.assign(n_cells, 0);

        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 1; i < im-1; i++) {
            for (int j = 1; j < jm-1; j++) {
                const double sthe   = sinthe_tbl[j];
                const double cthe   = costhe_tbl[j];
                const double cot    = cthe / sthe;
                const double inv_s2 = 1.0 / (sthe * sthe);

                for (int k = 1; k < km-1; k++) {
                    if (AtomUtils::is_land(h, i, j,   k) || AtomUtils::is_land(h, i, j+1, k)
                     || AtomUtils::is_land(h, i, j-1, k) || AtomUtils::is_land(h, i, j,   k+1)
                     || AtomUtils::is_land(h, i, j,   k-1))  continue;

                    const double vC = vn.x[i][j][k],   wC = wn.x[i][j][k];
                    const double vT = vn.x[i][j+1][k], vB = vn.x[i][j-1][k];
                    const double wT = wn.x[i][j+1][k], wB = wn.x[i][j-1][k];
                    const double vP = vn.x[i][j][k+1], vM = vn.x[i][j][k-1];
                    const double wP = wn.x[i][j][k+1], wM = wn.x[i][j][k-1];

                    const double d2v_the = (vT - 2.0*vC + vB) * inv_dthe2;
                    const double d2w_the = (wT - 2.0*wC + wB) * inv_dthe2;
                    const double d2v_phi = (vP - 2.0*vC + vM) * inv_dphi2;
                    const double d2w_phi = (wP - 2.0*wC + wM) * inv_dphi2;
                    const double dv_the  = (vT - vB) * inv_2dthe;
                    const double dw_the  = (wT - wB) * inv_2dthe;
                    const double dv_phi  = (vP - vM) * inv_2dphi;
                    const double dw_phi  = (wP - wM) * inv_2dphi;

                    const std::size_t p = HydHorizVisc::idx(i, j, k);
                    HydHorizVisc::lap_v[p] = d2v_the + cot * dv_the
                                           - (1.0 + cot*cot) * vC
                                           + d2v_phi * inv_s2
                                           - 2.0 * cthe * inv_s2 * dw_phi;
                    HydHorizVisc::lap_w[p] = d2w_the + cot * dw_the
                                           - (1.0 + cot*cot) * wC
                                           + d2w_phi * inv_s2
                                           + 2.0 * cthe * inv_s2 * dv_phi;
                    HydHorizVisc::ok[p] = 1;
                }
            }
        }
    }

    #pragma omp parallel for collapse(2) schedule(static)
    for (int i = 1; i < im-1; i++) {
        for (int j = 1; j < jm-1; j++) {

            CellGeometry geo;
            geo.rm       = rad.z[i];
            geo.rm2      = geo.rm * geo.rm;
            geo.exp_rm   = 1.0 / (geo.rm + 1.0);
            geo.exp_2_rm = geo.exp_rm * geo.exp_rm;

            geo.sinthe  = sinthe_tbl[j];
            geo.sinthe2 = geo.sinthe * geo.sinthe;
            geo.costhe  = costhe_tbl[j];

            // HORIZONTAL metric radius (HYD_METRIC_RADIUS). Identity when off.
            // rm / rm2 / exp_rm stay on the GRID coordinate: they are the radial
            // stretch, the layer thickness and the wall distance, not a planet radius.
            const double rh  = metricRadius(geo.rm);
            const double rh2 = rh * rh;
            geo.inv_rm              = 1.0 / rh;
            geo.inv_rm2             = 1.0 / rh2;
            geo.inv_rmsinthe        = 1.0 / (rh * geo.sinthe);
            geo.inv_rm2sinthe       = geo.inv_rm2 / geo.sinthe;
            geo.inv_rm2sinthe2      = geo.inv_rm2 / geo.sinthe2;
            geo.costhe_inv_rm2sinthe = geo.costhe * geo.inv_rm2sinthe;

            geo.inv_2dr   = inv_2dr;
            geo.inv_2dthe = inv_2dthe;
            geo.inv_2dphi = inv_2dphi;
            geo.inv_dr2   = inv_dr2;
            geo.inv_dthe2 = inv_dthe2;
            geo.inv_dphi2 = inv_dphi2;

            for (int k = 1; k < km-1; k++) {

                // Cache time-level-n values
                double tn_ijk   = tn.x[i][j][k];
                double un_ijk   = un.x[i][j][k];
                double vn_ijk   = vn.x[i][j][k];
                double wn_ijk   = wn.x[i][j][k];
                double cn_ijk   = cn.x[i][j][k];
                double tken_ijk = tken.x[i][j][k];
                double disn_ijk = disn.x[i][j][k];

                // --- Stage 1 ---
                cHydrosphereModel::RHS_Hydrosphere_Turb(i, j, k, geo);

                double kt1   = rhs_t.x[i][j][k];
                double ku1   = rhs_u.x[i][j][k];
                double kv1   = rhs_v.x[i][j][k];
                double kw1   = rhs_w.x[i][j][k];
                double kc1   = rhs_c.x[i][j][k];
                double ktke1 = rhs_tke.x[i][j][k];
                double kdis1 = rhs_dis.x[i][j][k];

                t.x[i][j][k]   = tn_ijk   + kt1   * half_dt;
                u.x[i][j][k]   = un_ijk   + ku1   * half_dt;
                v.x[i][j][k]   = vn_ijk   + kv1   * half_dt;
                w.x[i][j][k]   = wn_ijk   + kw1   * half_dt;
                c.x[i][j][k]   = cn_ijk   + kc1   * half_dt;
                tke.x[i][j][k] = safe_clamp(tken_ijk + ktke1 * half_dt, 0.0, tke_max_nd);
                dis.x[i][j][k] = std::max(1.0e-10,  disn_ijk + kdis1 * half_dt);

                // --- Stage 2 ---
                cHydrosphereModel::RHS_Hydrosphere_Turb(i, j, k, geo);

                double kt2   = rhs_t.x[i][j][k];
                double ku2   = rhs_u.x[i][j][k];
                double kv2   = rhs_v.x[i][j][k];
                double kw2   = rhs_w.x[i][j][k];
                double kc2   = rhs_c.x[i][j][k];
                double ktke2 = rhs_tke.x[i][j][k];
                double kdis2 = rhs_dis.x[i][j][k];

                t.x[i][j][k]   = tn_ijk   + kt2   * half_dt;
                u.x[i][j][k]   = un_ijk   + ku2   * half_dt;
                v.x[i][j][k]   = vn_ijk   + kv2   * half_dt;
                w.x[i][j][k]   = wn_ijk   + kw2   * half_dt;
                c.x[i][j][k]   = cn_ijk   + kc2   * half_dt;
                tke.x[i][j][k] = safe_clamp(tken_ijk + ktke2 * half_dt, 0.0, tke_max_nd);
                dis.x[i][j][k] = std::max(1.0e-10,  disn_ijk + kdis2 * half_dt);

                // --- Stage 3 ---
                cHydrosphereModel::RHS_Hydrosphere_Turb(i, j, k, geo);

                double kt3   = rhs_t.x[i][j][k];
                double ku3   = rhs_u.x[i][j][k];
                double kv3   = rhs_v.x[i][j][k];
                double kw3   = rhs_w.x[i][j][k];
                double kc3   = rhs_c.x[i][j][k];
                double ktke3 = rhs_tke.x[i][j][k];
                double kdis3 = rhs_dis.x[i][j][k];

                t.x[i][j][k]   = tn_ijk   + kt3   * dt;
                u.x[i][j][k]   = un_ijk   + ku3   * dt;
                v.x[i][j][k]   = vn_ijk   + kv3   * dt;
                w.x[i][j][k]   = wn_ijk   + kw3   * dt;
                c.x[i][j][k]   = cn_ijk   + kc3   * dt;
                tke.x[i][j][k] = safe_clamp(tken_ijk + ktke3 * dt, 0.0, tke_max_nd);
                dis.x[i][j][k] = std::max(1.0e-10,  disn_ijk + kdis3 * dt);

                // --- Stage 4 + Final combination ---
                cHydrosphereModel::RHS_Hydrosphere_Turb(i, j, k, geo);

                double kt4   = rhs_t.x[i][j][k];
                double ku4   = rhs_u.x[i][j][k];
                double kv4   = rhs_v.x[i][j][k];
                double kw4   = rhs_w.x[i][j][k];
                double kc4   = rhs_c.x[i][j][k];
                double ktke4 = rhs_tke.x[i][j][k];
                double kdis4 = rhs_dis.x[i][j][k];

                t.x[i][j][k]   = tn_ijk   + (kt1   + 2.0*kt2   + 2.0*kt3   + kt4)   * dt_sixth;
                u.x[i][j][k]   = un_ijk   + (ku1   + 2.0*ku2   + 2.0*ku3   + ku4)   * dt_sixth;
                v.x[i][j][k]   = vn_ijk   + (kv1   + 2.0*kv2   + 2.0*kv3   + kv4)   * dt_sixth;
                w.x[i][j][k]   = wn_ijk   + (kw1   + 2.0*kw2   + 2.0*kw3   + kw4)   * dt_sixth;
                c.x[i][j][k]   = cn_ijk   + (kc1   + 2.0*kc2   + 2.0*kc3   + kc4)   * dt_sixth;
                tke.x[i][j][k] = safe_clamp(
                    tken_ijk + (ktke1 + 2.0*ktke2 + 2.0*ktke3 + ktke4) * dt_sixth,
                    0.0, tke_max_nd);
                dis.x[i][j][k] = std::max(1.0e-10,
                    disn_ijk + (kdis1 + 2.0*kdis2 + 2.0*kdis3 + kdis4) * dt_sixth);
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for solveRungeKutta_Hydrosphere_Turb\n",
           elapsed.count() * 1e-9);

    cout << "      solveRungeKutta_Hydrosphere_Turb end" << endl;
}
