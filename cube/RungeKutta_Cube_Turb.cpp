/*
 * Atmosphere General Circulation Modell (AGCM) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * class to produce results by the Runge-Kutta solution scheme
*/

#include "cCubeModel.h"

using namespace std;



void cCubeModel::solveRungeKutta_Cube_Turb(){
/*
    cout << endl << "      solveRungeKutta_Cube_Turb begin" << endl;

    auto begin = std::chrono::high_resolution_clock::now();
*/
    const double half_dt  = 0.5 * dt;
    const double dt_sixth = dt / 6.0;

    // Grid-spacing reciprocals — constant for the entire grid
    const double inv_2dz = 1.0 / (2.0 * dr);
    const double inv_2dy = 1.0 / (2.0 * dy);
    const double inv_2dx = 1.0 / (2.0 * dx);
    const double inv_dz2 = 1.0 / (dr * dr);
    const double inv_dy2 = 1.0 / (dy * dy);
    const double inv_dx2 = 1.0 / (dx * dx);

    #pragma omp parallel for collapse(2) schedule(static)
    for(int i = 1; i < im-1; i++){
        for(int j = 1; j < jm-1; j++){

            // ---- Build geometry struct ONCE per (i,j) ----
            CellGeometry geo;

            geo.inv_2dz = inv_2dz;
            geo.inv_2dy = inv_2dy;
            geo.inv_2dx = inv_2dx;
            geo.inv_dz2 = inv_dz2;
            geo.inv_dy2 = inv_dy2;
            geo.inv_dx2 = inv_dx2;

            for(int k = 1; k < km-1; k++){

                // Solid cell: enforce no-slip and skip integration.
                if (AtomUtils::is_land(h, i, j, k)) {
                    u.x[i][j][k]  = un.x[i][j][k]  = 0.0;
                    v.x[i][j][k]  = vn.x[i][j][k]  = 0.0;
                    w.x[i][j][k]  = wn.x[i][j][k]  = 0.0;
                    continue;
                }

                // Cache time-level-n values (read once, used 4×)
                double un_ijk   = un.x[i][j][k];
                double vn_ijk   = vn.x[i][j][k];
                double wn_ijk   = wn.x[i][j][k];
                double tken_ijk = tken.x[i][j][k];
                double disn_ijk = disn.x[i][j][k];

                // --- Stage 1 ---
                cCubeModel::RHS_Cube_Turb(i, j, k, geo);

                double ku1   = rhs_u.x[i][j][k];
                double kv1   = rhs_v.x[i][j][k];
                double kw1   = rhs_w.x[i][j][k];
                double ktke1 = rhs_tke.x[i][j][k];
                double kdis1 = rhs_dis.x[i][j][k];

                u.x[i][j][k]   = un_ijk   + ku1   * half_dt;
                v.x[i][j][k]   = vn_ijk   + kv1   * half_dt;
                w.x[i][j][k]   = wn_ijk   + kw1   * half_dt;
                tke.x[i][j][k] = std::max(0.0,           tken_ijk + ktke1 * half_dt);
                dis.x[i][j][k] = std::max(1.0e-10, disn_ijk + kdis1 * half_dt);

                // --- Stage 2 ---
                cCubeModel::RHS_Cube_Turb(i, j, k, geo);

                double ku2   = rhs_u.x[i][j][k];
                double kv2   = rhs_v.x[i][j][k];
                double kw2   = rhs_w.x[i][j][k];
                double ktke2 = rhs_tke.x[i][j][k];
                double kdis2 = rhs_dis.x[i][j][k];

                u.x[i][j][k]   = un_ijk   + ku2   * half_dt;
                v.x[i][j][k]   = vn_ijk   + kv2   * half_dt;
                w.x[i][j][k]   = wn_ijk   + kw2   * half_dt;
                tke.x[i][j][k] = std::max(0.0,      tken_ijk + ktke2 * half_dt);
                dis.x[i][j][k] = std::max(1.0e-10,  disn_ijk + kdis2 * half_dt);

                // --- Stage 3 ---
                cCubeModel::RHS_Cube_Turb(i, j, k, geo);

                double ku3   = rhs_u.x[i][j][k];
                double kv3   = rhs_v.x[i][j][k];
                double kw3   = rhs_w.x[i][j][k];
                double ktke3 = rhs_tke.x[i][j][k];
                double kdis3 = rhs_dis.x[i][j][k];

                u.x[i][j][k]   = un_ijk   + ku3   * dt;
                v.x[i][j][k]   = vn_ijk   + kv3   * dt;
                w.x[i][j][k]   = wn_ijk   + kw3   * dt;
                tke.x[i][j][k] = std::max(0.0,      tken_ijk + ktke3 * dt);
                dis.x[i][j][k] = std::max(1.0e-10,  disn_ijk + kdis3 * dt);

                // --- Stage 4 + Final combination ---
                cCubeModel::RHS_Cube_Turb(i, j, k, geo);

                double ku4   = rhs_u.x[i][j][k];
                double kv4   = rhs_v.x[i][j][k];
                double kw4   = rhs_w.x[i][j][k];
                double ktke4 = rhs_tke.x[i][j][k];
                double kdis4 = rhs_dis.x[i][j][k];

                u.x[i][j][k]   = un_ijk   + (ku1   + 2.0*ku2   + 2.0*ku3   + ku4)   * dt_sixth;
                v.x[i][j][k]   = vn_ijk   + (kv1   + 2.0*kv2   + 2.0*kv3   + kv4)   * dt_sixth;
                w.x[i][j][k]   = wn_ijk   + (kw1   + 2.0*kw2   + 2.0*kw3   + kw4)   * dt_sixth;
                tke.x[i][j][k] = tken_ijk + (ktke1 + 2.0*ktke2 + 2.0*ktke3 + ktke4) * dt_sixth;
                dis.x[i][j][k] = disn_ijk + (kdis1 + 2.0*kdis2 + 2.0*kdis3 + kdis4) * dt_sixth;
            }
        }
    }
/*
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for solveRungeKutta_Cube_Turb\n", elapsed.count() * 1e-9);

    cout << "      solveRungeKutta_Cube_Turb end" << endl;
*/
}
