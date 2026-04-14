/*
 * Ocean General Circulation Model (OGCM) applied to turbulent flow
 * 4th-order Runge-Kutta for 7 prognostic fields: t, u, v, w, c, tke, dis
*/

#include "cHydrosphereModel.h"

using namespace std;


void cHydrosphereModel::solveRungeKutta_Hydrosphere_Turb(){
    cout << endl << "      solveRungeKutta_Hydrosphere_Turb begin" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    const double half_dt  = 0.5 * dt;
    const double dt_sixth = dt / 6.0;

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

            geo.inv_rm              = 1.0 / geo.rm;
            geo.inv_rm2             = 1.0 / geo.rm2;
            geo.inv_rmsinthe        = 1.0 / (geo.rm * geo.sinthe);
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
                tke.x[i][j][k] = std::max(0.0,      tken_ijk + ktke1 * half_dt);
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
                tke.x[i][j][k] = std::max(0.0,      tken_ijk + ktke2 * half_dt);
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
                tke.x[i][j][k] = std::max(0.0,      tken_ijk + ktke3 * dt);
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
                tke.x[i][j][k] = std::max(0.0,
                    tken_ijk + (ktke1 + 2.0*ktke2 + 2.0*ktke3 + ktke4) * dt_sixth);
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
