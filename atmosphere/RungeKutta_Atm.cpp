/*
 * Atmosphere General Circulation Modell (AGCM) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * class to produce results by the Runge-Kutta solution scheme
*/

#include "cAtmosphereModel.h"

using namespace std;


void cAtmosphereModel::solveRungeKutta_Atmosphere(){
    cout << endl << " .................... solveRungeKutta_Atmosphere begin" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    const double half_dt = 0.5 * dt;
    const double dt_sixth = dt / 6.0;

    // Velocity caps removed.  They didn't preserve mass: every clipped cell injected
    // a divergence that the (now-multi-sweep) pressure solver absorbed by relaxing
    // p_dyn → 0, which collapsed the whole system into a saturated-velocity /
    // zero-pressure degenerate state by iter 100.  Multi-sweep PressureSolverAtm is
    // the structural fix; caps were a band-aid that did more harm than good.

    // Grid-spacing reciprocals — constant for the entire grid
    const double inv_2dr   = 1.0 / (2.0 * dr);
    const double inv_2dthe = 1.0 / (2.0 * dthe);
    const double inv_2dphi = 1.0 / (2.0 * dphi);
    const double inv_dr2   = 1.0 / (dr * dr);
    const double inv_dthe2 = 1.0 / (dthe * dthe);
    const double inv_dphi2 = 1.0 / (dphi * dphi);

    // Precompute sin/cos tables — only depend on j
    std::vector<double> sinthe_tbl(jm), costhe_tbl(jm);

    for(int j = 0; j < jm; j++){
        sinthe_tbl[j] = sin(the.z[j]);
        if(sinthe_tbl[j] < 0.4) sinthe_tbl[j] = 0.4;
        costhe_tbl[j] = cos(the.z[j]);
    }

    #pragma omp parallel for collapse(2) schedule(static)
    for(int i = 1; i < im-1; i++){
        for(int j = 1; j < jm-1; j++){

            // ---- Build geometry struct ONCE per (i,j) ----
            CellGeometry geo;

            geo.rm     = rad.z[i];
            geo.rm2      = geo.rm * geo.rm;
            geo.exp_rm   = 1.0 / (geo.rm + 1.0);
            geo.exp_2_rm = geo.exp_rm * geo.exp_rm;

            geo.sinthe   = sinthe_tbl[j];
            geo.sinthe2  = geo.sinthe * geo.sinthe;
            geo.costhe   = costhe_tbl[j];

            geo.inv_rm          = 1.0 / geo.rm;
            geo.inv_rm2         = 1.0 / geo.rm2;
            geo.inv_rmsinthe    = 1.0 / (geo.rm * geo.sinthe);
            geo.inv_rm2sinthe   = geo.inv_rm2 / geo.sinthe;
            geo.inv_rm2sinthe2  = geo.inv_rm2 / geo.sinthe2;
            geo.costhe_inv_rm2sinthe = geo.costhe * geo.inv_rm2sinthe;

            geo.inv_2dr   = inv_2dr;
            geo.inv_2dthe = inv_2dthe;
            geo.inv_2dphi = inv_2dphi;
            geo.inv_dr2   = inv_dr2;
            geo.inv_dthe2 = inv_dthe2;
            geo.inv_dphi2 = inv_dphi2;

            for(int k = 1; k < km-1; k++){

                // Cache time-level-n values (read once, used 4×)
                double tn_ijk   = tn.x[i][j][k];
                double un_ijk   = un.x[i][j][k];
                double vn_ijk   = vn.x[i][j][k];
                double wn_ijk   = wn.x[i][j][k];
                double cn_ijk   = cn.x[i][j][k];
                double cldn_ijk = cloudn.x[i][j][k];
                double icen_ijk = icen.x[i][j][k];
                double grn_ijk  = grn.x[i][j][k];
                double co2n_ijk = co2n.x[i][j][k];

                // --- Stage 1 ---
                cAtmosphereModel::RHS_Atmosphere(i, j, k, geo);

                double kt1  = rhs_t.x[i][j][k];
                double ku1  = rhs_u.x[i][j][k];
                double kv1  = rhs_v.x[i][j][k];
                double kw1  = rhs_w.x[i][j][k];
                double kc1  = rhs_c.x[i][j][k];
                double kcl1 = rhs_cloud.x[i][j][k];
                double ki1  = rhs_ice.x[i][j][k];
                double kg1  = rhs_g.x[i][j][k];
                double kco1 = rhs_co2.x[i][j][k];

                t.x[i][j][k]     = tn_ijk   + kt1  * half_dt;
                u.x[i][j][k]     = un_ijk   + ku1  * half_dt;
                v.x[i][j][k]     = vn_ijk   + kv1  * half_dt;
                w.x[i][j][k]     = wn_ijk   + kw1  * half_dt;
                // Positivity clamp on the moisture scalars at every RK4 stage (see
                // RungeKutta_Atm_Turb.cpp): mass mixing ratios are physically >= 0, so the
                // prognostic field must never be driven negative by the spherical-operator
                // update at sharp coastal gradients.
                c.x[i][j][k]     = std::max(0.0, cn_ijk   + kc1  * half_dt);
                cloud.x[i][j][k] = std::max(0.0, cldn_ijk + kcl1 * half_dt);
                ice.x[i][j][k]   = std::max(0.0, icen_ijk + ki1  * half_dt);
                gr.x[i][j][k]    = std::max(0.0, grn_ijk  + kg1  * half_dt);
                co2.x[i][j][k]   = co2n_ijk + kco1 * half_dt;

                // --- Stage 2 ---
                cAtmosphereModel::RHS_Atmosphere(i, j, k, geo);

                double kt2  = rhs_t.x[i][j][k];
                double ku2  = rhs_u.x[i][j][k];
                double kv2  = rhs_v.x[i][j][k];
                double kw2  = rhs_w.x[i][j][k];
                double kc2  = rhs_c.x[i][j][k];
                double kcl2 = rhs_cloud.x[i][j][k];
                double ki2  = rhs_ice.x[i][j][k];
                double kg2  = rhs_g.x[i][j][k];
                double kco2 = rhs_co2.x[i][j][k];

                t.x[i][j][k]     = tn_ijk   + kt2  * half_dt;
                u.x[i][j][k]     = un_ijk   + ku2  * half_dt;
                v.x[i][j][k]     = vn_ijk   + kv2  * half_dt;
                w.x[i][j][k]     = wn_ijk   + kw2  * half_dt;
                c.x[i][j][k]     = std::max(0.0, cn_ijk   + kc2  * half_dt);
                cloud.x[i][j][k] = std::max(0.0, cldn_ijk + kcl2 * half_dt);
                ice.x[i][j][k]   = std::max(0.0, icen_ijk + ki2  * half_dt);
                gr.x[i][j][k]    = std::max(0.0, grn_ijk  + kg2  * half_dt);
                co2.x[i][j][k]   = co2n_ijk + kco2 * half_dt;

                // --- Stage 3 ---
                cAtmosphereModel::RHS_Atmosphere(i, j, k, geo);

                double kt3  = rhs_t.x[i][j][k];
                double ku3  = rhs_u.x[i][j][k];
                double kv3  = rhs_v.x[i][j][k];
                double kw3  = rhs_w.x[i][j][k];
                double kc3  = rhs_c.x[i][j][k];
                double kcl3 = rhs_cloud.x[i][j][k];
                double ki3  = rhs_ice.x[i][j][k];
                double kg3  = rhs_g.x[i][j][k];
                double kco3 = rhs_co2.x[i][j][k];

                t.x[i][j][k]     = tn_ijk   + kt3  * dt;
                u.x[i][j][k]     = un_ijk   + ku3  * dt;
                v.x[i][j][k]     = vn_ijk   + kv3  * dt;
                w.x[i][j][k]     = wn_ijk   + kw3  * dt;
                c.x[i][j][k]     = std::max(0.0, cn_ijk   + kc3  * dt);
                cloud.x[i][j][k] = std::max(0.0, cldn_ijk + kcl3 * dt);
                ice.x[i][j][k]   = std::max(0.0, icen_ijk + ki3  * dt);
                gr.x[i][j][k]    = std::max(0.0, grn_ijk  + kg3  * dt);
                co2.x[i][j][k]   = co2n_ijk + kco3 * dt;

                // --- Stage 4 + Final combination ---
                cAtmosphereModel::RHS_Atmosphere(i, j, k, geo);

                double kt4  = rhs_t.x[i][j][k];
                double ku4  = rhs_u.x[i][j][k];
                double kv4  = rhs_v.x[i][j][k];
                double kw4  = rhs_w.x[i][j][k];
                double kc4  = rhs_c.x[i][j][k];
                double kcl4 = rhs_cloud.x[i][j][k];
                double ki4  = rhs_ice.x[i][j][k];
                double kg4  = rhs_g.x[i][j][k];
                double kco4 = rhs_co2.x[i][j][k];

                t.x[i][j][k]     = tn_ijk   + (kt1  + 2.0*kt2  + 2.0*kt3  + kt4)  * dt_sixth;
                u.x[i][j][k]     = un_ijk   + (ku1  + 2.0*ku2  + 2.0*ku3  + ku4)  * dt_sixth;
                v.x[i][j][k]     = vn_ijk   + (kv1  + 2.0*kv2  + 2.0*kv3  + kv4)  * dt_sixth;
                w.x[i][j][k]     = wn_ijk   + (kw1  + 2.0*kw2  + 2.0*kw3  + kw4)  * dt_sixth;
                c.x[i][j][k]     = std::max(0.0, cn_ijk   + (kc1  + 2.0*kc2  + 2.0*kc3  + kc4)  * dt_sixth);
                cloud.x[i][j][k] = std::max(0.0, cldn_ijk + (kcl1 + 2.0*kcl2 + 2.0*kcl3 + kcl4) * dt_sixth);
                ice.x[i][j][k]   = std::max(0.0, icen_ijk + (ki1  + 2.0*ki2  + 2.0*ki3  + ki4)  * dt_sixth);
                gr.x[i][j][k]    = std::max(0.0, grn_ijk  + (kg1  + 2.0*kg2  + 2.0*kg3  + kg4)  * dt_sixth);
                co2.x[i][j][k]   = co2n_ijk + (kco1 + 2.0*kco2 + 2.0*kco3 + kco4) * dt_sixth;
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for solveRungeKutta_Atmosphere\n", elapsed.count() * 1e-9);

    cout << " .................... solveRungeKutta_Atmosphere end" << endl;
}
