/*
 * Atmosphere General Circulation Modell (AGCM) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 *
 * class to combine the right hand sides of the differential equations for the Runge-Kutta scheme
*/

#include <iostream>
#include <cmath>
#include "cCubeModel.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;


void cCubeModel::RHS_Cube_Turb(int i, int j, int k, const CellGeometry& geo){

    const double inv_2dz  = geo.inv_2dz;
    const double inv_2dy  = geo.inv_2dy;
    const double inv_2dx  = geo.inv_2dx;
    const double inv_dz2  = geo.inv_dz2;
    const double inv_dy2  = geo.inv_dy2;
    const double inv_dx2  = geo.inv_dx2;

    // Cache local cell velocities
    double u_ijk = u.x[i][j][k];
    double v_ijk = v.x[i][j][k];
    double w_ijk = w.x[i][j][k];

    bool land_ijk = is_land(h, i, j, k);

    bool z_flag = false, y_flag = false, x_flag = false;

    // ---- First-order derivatives storage ----
    double dudz, dvdz, dwdz, dpdz;
    double dudy, dvdy, dwdy, dpdy;
    double dudx, dvdx, dwdx, dpdx;

    // ---- Second-order derivatives storage ----
    double d2udz2, d2vdz2, d2wdz2;
    double d2udy2, d2vdy2, d2wdy2;
    double d2udx2, d2vdx2, d2wdx2;


    // ===== Z-direction derivatives =====
    if(i < im-2 && land_ijk && is_air(h, i+1, j, k)){
        // Forward-biased stencil near land-air boundary
        #define COMPUTE_DZ_B(FIELD, d1, d2) \
            d1 = (-3.0*FIELD.x[i][j][k] + 4.0*FIELD.x[i+1][j][k] - FIELD.x[i+2][j][k]) * inv_2dz; \
            d2 = (FIELD.x[i][j][k] - 2.0*FIELD.x[i+1][j][k] + FIELD.x[i+2][j][k]) * inv_dz2;

        COMPUTE_DZ_B(u, dudz, d2udz2)
        COMPUTE_DZ_B(v, dvdz, d2vdz2)
        COMPUTE_DZ_B(w, dwdz, d2wdz2)
        dpdz = (-3.0*p_dyn.x[i][j][k] + 4.0*p_dyn.x[i+1][j][k] - p_dyn.x[i+2][j][k]) * inv_2dz;
        z_flag = true;
        #undef COMPUTE_DZ_B
    }
    else if(i == 0){
        #define COMPUTE_DZ_F(FIELD, d1, d2) \
            d1 = (-3.0*FIELD.x[0][j][k] + 4.0*FIELD.x[1][j][k] - FIELD.x[2][j][k]) * inv_2dz; \
            d2 = (FIELD.x[0][j][k] - 2.0*FIELD.x[1][j][k] + FIELD.x[2][j][k]) * inv_dz2;

        COMPUTE_DZ_F(u, dudz, d2udz2)
        COMPUTE_DZ_F(v, dvdz, d2vdz2)
        COMPUTE_DZ_F(w, dwdz, d2wdz2)
        dpdz = (-3.0*p_dyn.x[0][j][k] + 4.0*p_dyn.x[1][j][k] - p_dyn.x[2][j][k]) * inv_2dz;
        z_flag = true;
        #undef COMPUTE_DZ_F
    }
    if(!z_flag){
        // Central differences (default path — most common)
        #define COMPUTE_DZ_C(FIELD, d1, d2) \
            d1 = (FIELD.x[i+1][j][k] - FIELD.x[i-1][j][k]) * inv_2dz; \
          d2 = (FIELD.x[i+1][j][k] - 2.0*FIELD.x[i][j][k] + FIELD.x[i-1][j][k]) * inv_dz2;

        COMPUTE_DZ_C(u, dudz, d2udz2)
        COMPUTE_DZ_C(v, dvdz, d2vdz2)
        COMPUTE_DZ_C(w, dwdz, d2wdz2)
        dpdz = (p_dyn.x[i+1][j][k] - p_dyn.x[i-1][j][k]) * inv_2dz;
        #undef COMPUTE_DZ_C
    }


    // ===== Y-direction derivatives =====
    if(j > 2 && j < jm-3){
        if(land_ijk && is_air(h, i, j+1, k) && is_air(h, i, j+2, k)){
            #define COMPUTE_DY_B(FIELD, d1, d2) \
                d1 = (-3.0*FIELD.x[i][j][k] + 4.0*FIELD.x[i][j+1][k] - FIELD.x[i][j+2][k]) * inv_2dy; \
                d2 = (FIELD.x[i][j][k] - 2.0*FIELD.x[i][j+1][k] + FIELD.x[i][j+2][k]) * inv_dy2;

            COMPUTE_DY_B(u, dudy, d2udy2) COMPUTE_DY_B(v, dvdy, d2vdy2)
            COMPUTE_DY_B(w, dwdy, d2wdy2)
            dpdy = (-3.0*p_dyn.x[i][j][k] + 4.0*p_dyn.x[i][j+1][k] - p_dyn.x[i][j+2][k]) * inv_2dy;
            y_flag = true;
            #undef COMPUTE_DY_B
        }
        if(!y_flag && land_ijk && is_air(h, i, j-1, k) && is_air(h, i, j-2, k)){
            #define COMPUTE_DY_C(FIELD, d1, d2) \
                d1 = -(-3.0*FIELD.x[i][j][k] + 4.0*FIELD.x[i][j-1][k] - FIELD.x[i][j-2][k]) * inv_2dy; \
                d2 = -(FIELD.x[i][j][k] - 2.0*FIELD.x[i][j-1][k] + FIELD.x[i][j-2][k]) * inv_dy2;

            COMPUTE_DY_C(u, dudy, d2udy2) COMPUTE_DY_C(v, dvdy, d2vdy2)
            COMPUTE_DY_C(w, dwdy, d2wdy2)
            dpdy = -(-3.0*p_dyn.x[i][j][k] + 4.0*p_dyn.x[i][j-1][k] - p_dyn.x[i][j-2][k]) * inv_2dy;
            y_flag = true;
            #undef COMPUTE_DY_C
        }
    }
    else if(j == 0){
        #define COMPUTE_DY_NP(FIELD, d1, d2) \
            d1 = (-3.0*FIELD.x[i][0][k] + 4.0*FIELD.x[i][1][k] - FIELD.x[i][2][k]) * inv_2dy; \
            d2 = (FIELD.x[i][0][k] - 2.0*FIELD.x[i][1][k] + FIELD.x[i][2][k]) * inv_dy2;

        COMPUTE_DY_NP(u, dudy, d2udy2) COMPUTE_DY_NP(v, dvdy, d2vdy2)
        COMPUTE_DY_NP(w, dwdy, d2wdy2)
        dpdy = (-3.0*p_dyn.x[i][0][k] + 4.0*p_dyn.x[i][1][k] - p_dyn.x[i][2][k]) * inv_2dy;
        y_flag = true;
        #undef COMPUTE_DY_NP
    }
    else if(j == jm-1){
        #define COMPUTE_DY_SP(FIELD, d1, d2) \
            d1 = -(-3.0*FIELD.x[i][jm-1][k] + 4.0*FIELD.x[i][jm-2][k] - FIELD.x[i][jm-3][k]) * inv_2dy; \
            d2 = -(FIELD.x[i][jm-1][k] - 2.0*FIELD.x[i][jm-2][k] + FIELD.x[i][jm-3][k]) * inv_dy2;

        COMPUTE_DY_SP(u, dudy, d2udy2) COMPUTE_DY_SP(v, dvdy, d2vdy2)
        COMPUTE_DY_SP(w, dwdy, d2wdy2)
        dpdy = -(-3.0*p_dyn.x[i][jm-1][k] + 4.0*p_dyn.x[i][jm-2][k] - p_dyn.x[i][jm-3][k]) * inv_2dy;
        y_flag = true;
        #undef COMPUTE_DY_SP
    }

    if(!y_flag){
        #define COMPUTE_DY_A(FIELD, d1, d2) \
            d1 = (FIELD.x[i][j+1][k] - FIELD.x[i][j-1][k]) * inv_2dy; \
            d2 = (FIELD.x[i][j+1][k] - 2.0*FIELD.x[i][j][k] + FIELD.x[i][j-1][k]) * inv_dy2;

        COMPUTE_DY_A(u, dudy, d2udy2) COMPUTE_DY_A(v, dvdy, d2vdy2)
        COMPUTE_DY_A(w, dwdy, d2wdy2)
        dpdy = (p_dyn.x[i][j+1][k] - p_dyn.x[i][j-1][k]) * inv_2dy;
        #undef COMPUTE_DY_A
    }


    // ===== X-direction derivatives =====
    if(k > 2 && k < km-3){
        if(land_ijk && is_air(h, i, j, k+1) && is_air(h, i, j, k+2)){
          #define COMPUTE_DX_B(FIELD, d1, d2) \
                d1 = (-3.0*FIELD.x[i][j][k] + 4.0*FIELD.x[i][j][k+1] - FIELD.x[i][j][k+2]) * inv_2dx; \
                d2 = (FIELD.x[i][j][k] - 2.0*FIELD.x[i][j][k+1] + FIELD.x[i][j][k+2]) * inv_dx2;

            COMPUTE_DX_B(u, dudx, d2udx2) COMPUTE_DX_B(v, dvdx, d2vdx2)
            COMPUTE_DX_B(w, dwdx, d2wdx2)
            dpdx = (-3.0*p_dyn.x[i][j][k] + 4.0*p_dyn.x[i][j][k+1] - p_dyn.x[i][j][k+2]) * inv_2dx;
            x_flag = true;
            #undef COMPUTE_DX_B
        }
        if(!x_flag && land_ijk && is_air(h, i, j, k-1) && is_air(h, i, j, k-2)){
            #define COMPUTE_DX_C(FIELD, d1, d2) \
                d1 = -(-3.0*FIELD.x[i][j][k] + 4.0*FIELD.x[i][j][k-1] - FIELD.x[i][j][k-2]) * inv_2dx; \
                d2 = -(FIELD.x[i][j][k] - 2.0*FIELD.x[i][j][k-1] + FIELD.x[i][j][k-2]) * inv_dx2;

            COMPUTE_DX_C(u, dudx, d2udx2) COMPUTE_DX_C(v, dvdx, d2vdx2)
            COMPUTE_DX_C(w, dwdx, d2wdx2)
            dpdx = -(-3.0*p_dyn.x[i][j][k] + 4.0*p_dyn.x[i][j][k-1] - p_dyn.x[i][j][k-2]) * inv_2dx;
            x_flag = true;
            #undef COMPUTE_DX_C
        }
    }
    else if(k == 0){
        #define COMPUTE_DX_W(FIELD, d1, d2) \
            d1 = (-3.0*FIELD.x[i][j][0] + 4.0*FIELD.x[i][j][1] - FIELD.x[i][j][2]) * inv_2dx; \
            d2 = (FIELD.x[i][j][0] - 2.0*FIELD.x[i][j][1] + FIELD.x[i][j][2]) * inv_dx2;

        COMPUTE_DX_W(u, dudx, d2udx2) COMPUTE_DX_W(v, dvdx, d2vdx2)
        COMPUTE_DX_W(w, dwdx, d2wdx2)
        dpdx = (-3.0*p_dyn.x[i][j][0] + 4.0*p_dyn.x[i][j][1] - p_dyn.x[i][j][2]) * inv_2dx;
        x_flag = true;
        #undef COMPUTE_DX_W
    }
    else if(k == km-1){
        #define COMPUTE_DX_E(FIELD, d1, d2) \
            d1 = -(-3.0*FIELD.x[i][j][km-1] + 4.0*FIELD.x[i][j][km-2] - FIELD.x[i][j][km-3]) * inv_2dx; \
            d2 = -(FIELD.x[i][j][km-1] - 2.0*FIELD.x[i][j][km-2] + FIELD.x[i][j][km-3]) * inv_dx2;

        COMPUTE_DX_E(u, dudx, d2udx2) COMPUTE_DX_E(v, dvdx, d2vdx2)
        COMPUTE_DX_E(w, dwdx, d2wdx2)
        dpdx = -(-3.0*p_dyn.x[i][j][km-1] + 4.0*p_dyn.x[i][j][km-2] - p_dyn.x[i][j][km-3]) * inv_2dx;
        x_flag = true;
        #undef COMPUTE_DX_E
    }

    if(!x_flag){
        #define COMPUTE_DX_A(FIELD, d1, d2) \
            d1 = (FIELD.x[i][j][k+1] - FIELD.x[i][j][k-1]) * inv_2dx; \
            d2 = (FIELD.x[i][j][k+1] - 2.0*FIELD.x[i][j][k] + FIELD.x[i][j][k-1]) * inv_dx2;

        COMPUTE_DX_A(u, dudx, d2udx2) COMPUTE_DX_A(v, dvdx, d2vdx2)
        COMPUTE_DX_A(w, dwdx, d2wdx2)
        dpdx = (p_dyn.x[i][j][k+1] - p_dyn.x[i][j][k-1]) * inv_2dx;
        #undef COMPUTE_DX_A
    }



    const double diff_vel_re = 1.0 / re_turb;

    // ===== Advection (u·∇) — uniform Cartesian =====
    const double transport_u = u_ijk * dudz + v_ijk * dudy + w_ijk * dudx;
    const double transport_v = u_ijk * dvdz + v_ijk * dvdy + w_ijk * dvdx;
    const double transport_w = u_ijk * dwdz + v_ijk * dwdy + w_ijk * dwdx;

    // ===== Diffusion ∇²u — uniform Cartesian =====
    const double diffusion_u = (d2udz2 + d2udy2 + d2udx2) * diff_vel_re;
    const double diffusion_v = (d2vdz2 + d2vdy2 + d2vdx2) * diff_vel_re;
    const double diffusion_w = (d2wdz2 + d2wdy2 + d2wdx2) * diff_vel_re;

    // ===== RHS assembly =====
    rhs_u.x[i][j][k] = -dpdz - transport_u + diffusion_u;
    rhs_v.x[i][j][k] = -dpdy - transport_v + diffusion_v;
    rhs_w.x[i][j][k] = -dpdx - transport_w + diffusion_w;

    aux_u.x[i][j][k] = rhs_u.x[i][j][k] + dpdz;
    aux_v.x[i][j][k] = rhs_v.x[i][j][k] + dpdy;
    aux_w.x[i][j][k] = rhs_w.x[i][j][k] + dpdx;
 }
