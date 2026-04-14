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
#include "cAtmosphereModel.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;


void cAtmosphereModel::RHS_Atmosphere_Turb(int i, int j, int k, const CellGeometry& geo){

    // All geometric quantities come from the precomputed struct —
    // NO sin(), cos(), division, or reciprocal computation here.
    const double exp_rm             = geo.exp_rm;
    const double exp_2_rm           = geo.exp_2_rm;
    const double sinthe             = geo.sinthe;
    const double sinthe2            = geo.sinthe2;
    const double costhe             = geo.costhe;
    const double inv_rm             = geo.inv_rm;
    const double inv_rm2            = geo.inv_rm2;
    const double inv_rmsinthe       = geo.inv_rmsinthe;
    const double inv_rm2sinthe      = geo.inv_rm2sinthe;
    const double inv_rm2sinthe2     = geo.inv_rm2sinthe2;
    const double costhe_inv_rm2sinthe = geo.costhe_inv_rm2sinthe;

    const double inv_2dr   = geo.inv_2dr;
    const double inv_2dthe = geo.inv_2dthe;
    const double inv_2dphi = geo.inv_2dphi;
    const double inv_dr2   = geo.inv_dr2;
    const double inv_dthe2 = geo.inv_dthe2;
    const double inv_dphi2 = geo.inv_dphi2;

    // Cache local cell velocities
    double u_ijk = u.x[i][j][k];
    double v_ijk = v.x[i][j][k];
    double w_ijk = w.x[i][j][k];

    bool land_ijk = is_land(h, i, j, k);

    bool r_flag = false, the_flag = false, phi_flag = false;

    // ---- First-order derivatives storage ----
    double dudr, dvdr, dwdr, dtdr, dpdr, dcdr, dclouddr, dicedr, dgdr, dcodr, dtkedr, ddisdr, dnuedr;
    double dudthe, dvdthe, dwdthe, dtdthe, dpdthe, dcdthe, dclouddthe, dicedthe, dgdthe, dcodthe, dtkedthe, ddisdthe, dnuedthe;
    double dudphi, dvdphi, dwdphi, dtdphi, dpdphi, dcdphi, dclouddphi, dicedphi, dgdphi, dcodphi, dtkedphi, ddisdphi, dnuedphi;

    // ---- Second-order derivatives storage ----
    double d2udr2, d2vdr2, d2wdr2, d2tdr2, d2cdr2, d2clouddr2, d2icedr2, d2gdr2, d2codr2, d2tkedr2, d2disdr2;
    double d2udthe2, d2vdthe2, d2wdthe2, d2tdthe2, d2cdthe2, d2clouddthe2, d2icedthe2, d2gdthe2, d2codthe2, d2tkedthe2, d2disdthe2;
    double d2udphi2, d2vdphi2, d2wdphi2, d2tdphi2, d2cdphi2, d2clouddphi2, d2icedphi2, d2gdphi2, d2codphi2, d2tkedphi2, d2disdphi2;


    // ===== R-direction derivatives =====
    if(i < im-2 && land_ijk && is_air(h, i+1, j, k)){
        // Forward-biased stencil near land-air boundary
        #define COMPUTE_DR_B(FIELD, d1, d2) \
            d1 = (-3.0*FIELD.x[i][j][k] + 4.0*FIELD.x[i+1][j][k] - FIELD.x[i+2][j][k]) * inv_2dr; \
            d2 = (FIELD.x[i][j][k] - 2.0*FIELD.x[i+1][j][k] + FIELD.x[i+2][j][k]) * inv_dr2;

        COMPUTE_DR_B(u, dudr, d2udr2)
        COMPUTE_DR_B(v, dvdr, d2vdr2)
        COMPUTE_DR_B(w, dwdr, d2wdr2)
        COMPUTE_DR_B(t, dtdr, d2tdr2)
        COMPUTE_DR_B(c, dcdr, d2cdr2)
        COMPUTE_DR_B(cloud, dclouddr, d2clouddr2)
        COMPUTE_DR_B(ice, dicedr, d2icedr2)
        COMPUTE_DR_B(gr, dgdr, d2gdr2)
        COMPUTE_DR_B(co2, dcodr, d2codr2)
        COMPUTE_DR_B(tke, dtkedr, d2tkedr2)
        COMPUTE_DR_B(dis, ddisdr, d2disdr2)
        dpdr   = (-3.0*p_dyn.x[i][j][k] + 4.0*p_dyn.x[i+1][j][k] - p_dyn.x[i+2][j][k]) * inv_2dr;
        dnuedr = (-3.0*nue.x[i][j][k]   + 4.0*nue.x[i+1][j][k]   - nue.x[i+2][j][k])   * inv_2dr;
        r_flag = true;
        #undef COMPUTE_DR_B
    }
    else if(i == 0){
        #define COMPUTE_DR_F(FIELD, d1, d2) \
            d1 = (-3.0*FIELD.x[0][j][k] + 4.0*FIELD.x[1][j][k] - FIELD.x[2][j][k]) * inv_2dr; \
            d2 = (FIELD.x[0][j][k] - 2.0*FIELD.x[1][j][k] + FIELD.x[2][j][k]) * inv_dr2;

        COMPUTE_DR_F(u, dudr, d2udr2)
        COMPUTE_DR_F(v, dvdr, d2vdr2)
        COMPUTE_DR_F(w, dwdr, d2wdr2)
        COMPUTE_DR_F(t, dtdr, d2tdr2)
        COMPUTE_DR_F(c, dcdr, d2cdr2)
        COMPUTE_DR_F(cloud, dclouddr, d2clouddr2)
        COMPUTE_DR_F(ice, dicedr, d2icedr2)
        COMPUTE_DR_F(gr, dgdr, d2gdr2)
        COMPUTE_DR_F(co2, dcodr, d2codr2)
        COMPUTE_DR_F(tke, dtkedr, d2tkedr2)
        COMPUTE_DR_F(dis, ddisdr, d2disdr2)
        dpdr   = (-3.0*p_dyn.x[0][j][k] + 4.0*p_dyn.x[1][j][k] - p_dyn.x[2][j][k]) * inv_2dr;
        dnuedr = (-3.0*nue.x[0][j][k]   + 4.0*nue.x[1][j][k]   - nue.x[2][j][k])   * inv_2dr;
        r_flag = true;
        #undef COMPUTE_DR_F
    }

    if(!r_flag){
        // Central differences (default path — most common)
        #define COMPUTE_DR_C(FIELD, d1, d2) \
            d1 = (FIELD.x[i+1][j][k] - FIELD.x[i-1][j][k]) * inv_2dr; \
            d2 = (FIELD.x[i+1][j][k] - 2.0*FIELD.x[i][j][k] + FIELD.x[i-1][j][k]) * inv_dr2;

        COMPUTE_DR_C(u, dudr, d2udr2)
        COMPUTE_DR_C(v, dvdr, d2vdr2)
        COMPUTE_DR_C(w, dwdr, d2wdr2)
        COMPUTE_DR_C(t, dtdr, d2tdr2)
        COMPUTE_DR_C(c, dcdr, d2cdr2)
        COMPUTE_DR_C(cloud, dclouddr, d2clouddr2)
        COMPUTE_DR_C(ice, dicedr, d2icedr2)
        COMPUTE_DR_C(gr, dgdr, d2gdr2)
        COMPUTE_DR_C(co2, dcodr, d2codr2)
        COMPUTE_DR_C(tke, dtkedr, d2tkedr2)
        COMPUTE_DR_C(dis, ddisdr, d2disdr2)
        dpdr   = (p_dyn.x[i+1][j][k] - p_dyn.x[i-1][j][k]) * inv_2dr;
        dnuedr = (nue.x[i+1][j][k]   - nue.x[i-1][j][k])   * inv_2dr;
        #undef COMPUTE_DR_C

        // At the first fluid cell directly above terrain (i-1 is a land cell), the
        // centered tke/dis stencil reads tke[i-1]=0 (land), producing a large
        // negative Laplacian that pulls tke → 0 and hence nue → 0 at surface heights.
        // Apply a Neumann (zero-gradient) BC at the terrain surface: replace
        // tke[i-1] with tke[i] so the diffusion flux is one-sided into the atmosphere.
        if (!land_ijk && i > 0 && is_land(h, i-1, j, k)) {
            dtkedr   = (tke.x[i+1][j][k] - tke.x[i][j][k]) * inv_2dr;
            d2tkedr2 = (tke.x[i+1][j][k] - tke.x[i][j][k]) * inv_dr2;
            ddisdr   = (dis.x[i+1][j][k] - dis.x[i][j][k]) * inv_2dr;
            d2disdr2 = (dis.x[i+1][j][k] - dis.x[i][j][k]) * inv_dr2;
        }
    }


    // ===== THETA-direction derivatives =====
    if(j > 2 && j < jm-3){
        if(land_ijk && is_air(h, i, j+1, k) && is_air(h, i, j+2, k)){
            #define COMPUTE_DTHE_B(FIELD, d1, d2) \
                d1 = (-3.0*FIELD.x[i][j][k] + 4.0*FIELD.x[i][j+1][k] - FIELD.x[i][j+2][k]) * inv_2dthe; \
                d2 = (FIELD.x[i][j][k] - 2.0*FIELD.x[i][j+1][k] + FIELD.x[i][j+2][k]) * inv_dthe2;

            COMPUTE_DTHE_B(u, dudthe, d2udthe2) COMPUTE_DTHE_B(v, dvdthe, d2vdthe2)
            COMPUTE_DTHE_B(w, dwdthe, d2wdthe2) COMPUTE_DTHE_B(t, dtdthe, d2tdthe2)
            COMPUTE_DTHE_B(c, dcdthe, d2cdthe2) COMPUTE_DTHE_B(cloud, dclouddthe, d2clouddthe2)
            COMPUTE_DTHE_B(ice, dicedthe, d2icedthe2) COMPUTE_DTHE_B(gr, dgdthe, d2gdthe2)
            COMPUTE_DTHE_B(co2, dcodthe, d2codthe2)
            COMPUTE_DTHE_B(tke, dtkedthe, d2tkedthe2) COMPUTE_DTHE_B(dis, ddisdthe, d2disdthe2)
            dpdthe   = (-3.0*p_dyn.x[i][j][k] + 4.0*p_dyn.x[i][j+1][k] - p_dyn.x[i][j+2][k]) * inv_2dthe;
            dnuedthe = (-3.0*nue.x[i][j][k]   + 4.0*nue.x[i][j+1][k]   - nue.x[i][j+2][k])   * inv_2dthe;
            the_flag = true;
            #undef COMPUTE_DTHE_B
        }
        if(!the_flag && land_ijk && is_air(h, i, j-1, k) && is_air(h, i, j-2, k)){
            #define COMPUTE_DTHE_C(FIELD, d1, d2) \
                d1 = -(-3.0*FIELD.x[i][j][k] + 4.0*FIELD.x[i][j-1][k] - FIELD.x[i][j-2][k]) * inv_2dthe; \
                d2 = -(FIELD.x[i][j][k] - 2.0*FIELD.x[i][j-1][k] + FIELD.x[i][j-2][k]) * inv_dthe2;

            COMPUTE_DTHE_C(u, dudthe, d2udthe2) COMPUTE_DTHE_C(v, dvdthe, d2vdthe2)
            COMPUTE_DTHE_C(w, dwdthe, d2wdthe2) COMPUTE_DTHE_C(t, dtdthe, d2tdthe2)
            COMPUTE_DTHE_C(c, dcdthe, d2cdthe2) COMPUTE_DTHE_C(cloud, dclouddthe, d2clouddthe2)
            COMPUTE_DTHE_C(ice, dicedthe, d2icedthe2) COMPUTE_DTHE_C(gr, dgdthe, d2gdthe2)
            COMPUTE_DTHE_C(co2, dcodthe, d2codthe2)
            COMPUTE_DTHE_C(tke, dtkedthe, d2tkedthe2) COMPUTE_DTHE_C(dis, ddisdthe, d2disdthe2)
            dpdthe   = -(-3.0*p_dyn.x[i][j][k] + 4.0*p_dyn.x[i][j-1][k] - p_dyn.x[i][j-2][k]) * inv_2dthe;
            dnuedthe = -(-3.0*nue.x[i][j][k]   + 4.0*nue.x[i][j-1][k]   - nue.x[i][j-2][k])   * inv_2dthe;
            the_flag = true;
            #undef COMPUTE_DTHE_C
        }
    }
    else if(j == 0){
        #define COMPUTE_DTHE_NP(FIELD, d1, d2) \
            d1 = (-3.0*FIELD.x[i][0][k] + 4.0*FIELD.x[i][1][k] - FIELD.x[i][2][k]) * inv_2dthe; \
            d2 = (FIELD.x[i][0][k] - 2.0*FIELD.x[i][1][k] + FIELD.x[i][2][k]) * inv_dthe2;

        COMPUTE_DTHE_NP(u, dudthe, d2udthe2) COMPUTE_DTHE_NP(v, dvdthe, d2vdthe2)
        COMPUTE_DTHE_NP(w, dwdthe, d2wdthe2) COMPUTE_DTHE_NP(t, dtdthe, d2tdthe2)
        COMPUTE_DTHE_NP(c, dcdthe, d2cdthe2) COMPUTE_DTHE_NP(cloud, dclouddthe, d2clouddthe2)
        COMPUTE_DTHE_NP(ice, dicedthe, d2icedthe2) COMPUTE_DTHE_NP(gr, dgdthe, d2gdthe2)
        COMPUTE_DTHE_NP(co2, dcodthe, d2codthe2)
        COMPUTE_DTHE_NP(tke, dtkedthe, d2tkedthe2) COMPUTE_DTHE_NP(dis, ddisdthe, d2disdthe2)
        dpdthe   = (-3.0*p_dyn.x[i][0][k] + 4.0*p_dyn.x[i][1][k] - p_dyn.x[i][2][k]) * inv_2dthe;
        dnuedthe = (-3.0*nue.x[i][0][k]   + 4.0*nue.x[i][1][k]   - nue.x[i][2][k])   * inv_2dthe;
        the_flag = true;
        #undef COMPUTE_DTHE_NP
    }
    else if(j == jm-1){
        #define COMPUTE_DTHE_SP(FIELD, d1, d2) \
            d1 = -(-3.0*FIELD.x[i][jm-1][k] + 4.0*FIELD.x[i][jm-2][k] - FIELD.x[i][jm-3][k]) * inv_2dthe; \
            d2 = -(FIELD.x[i][jm-1][k] - 2.0*FIELD.x[i][jm-2][k] + FIELD.x[i][jm-3][k]) * inv_dthe2;

        COMPUTE_DTHE_SP(u, dudthe, d2udthe2) COMPUTE_DTHE_SP(v, dvdthe, d2vdthe2)
        COMPUTE_DTHE_SP(w, dwdthe, d2wdthe2) COMPUTE_DTHE_SP(t, dtdthe, d2tdthe2)
        COMPUTE_DTHE_SP(c, dcdthe, d2cdthe2) COMPUTE_DTHE_SP(cloud, dclouddthe, d2clouddthe2)
        COMPUTE_DTHE_SP(ice, dicedthe, d2icedthe2) COMPUTE_DTHE_SP(gr, dgdthe, d2gdthe2)
        COMPUTE_DTHE_SP(co2, dcodthe, d2codthe2)
        COMPUTE_DTHE_SP(tke, dtkedthe, d2tkedthe2) COMPUTE_DTHE_SP(dis, ddisdthe, d2disdthe2)
        dpdthe   = -(-3.0*p_dyn.x[i][jm-1][k] + 4.0*p_dyn.x[i][jm-2][k] - p_dyn.x[i][jm-3][k]) * inv_2dthe;
        dnuedthe = -(-3.0*nue.x[i][jm-1][k]   + 4.0*nue.x[i][jm-2][k]   - nue.x[i][jm-3][k])   * inv_2dthe;
        the_flag = true;
        #undef COMPUTE_DTHE_SP
    }

    if(!the_flag){
        #define COMPUTE_DTHE_A(FIELD, d1, d2) \
            d1 = (FIELD.x[i][j+1][k] - FIELD.x[i][j-1][k]) * inv_2dthe; \
            d2 = (FIELD.x[i][j+1][k] - 2.0*FIELD.x[i][j][k] + FIELD.x[i][j-1][k]) * inv_dthe2;

        COMPUTE_DTHE_A(u, dudthe, d2udthe2) COMPUTE_DTHE_A(v, dvdthe, d2vdthe2)
        COMPUTE_DTHE_A(w, dwdthe, d2wdthe2) COMPUTE_DTHE_A(t, dtdthe, d2tdthe2)
        COMPUTE_DTHE_A(c, dcdthe, d2cdthe2) COMPUTE_DTHE_A(cloud, dclouddthe, d2clouddthe2)
        COMPUTE_DTHE_A(ice, dicedthe, d2icedthe2) COMPUTE_DTHE_A(gr, dgdthe, d2gdthe2)
        COMPUTE_DTHE_A(co2, dcodthe, d2codthe2)
        COMPUTE_DTHE_A(tke, dtkedthe, d2tkedthe2) COMPUTE_DTHE_A(dis, ddisdthe, d2disdthe2)
        dpdthe   = (p_dyn.x[i][j+1][k] - p_dyn.x[i][j-1][k]) * inv_2dthe;
        dnuedthe = (nue.x[i][j+1][k]   - nue.x[i][j-1][k])   * inv_2dthe;
        #undef COMPUTE_DTHE_A

        // Neumann BC for tke/dis at theta-direction land boundaries.
        // Air cells adjacent to land would read tke/dis=0 from the zeroed land
        // cell, creating a false large gradient that amplifies D_w cross-diffusion.
        // Replace the land neighbour's value with the current cell's (zero-gradient).
        if (!land_ijk) {
            const bool neumann_jm1 = (j > 0    && is_land(h, i, j-1, k));
            const bool neumann_jp1 = (j < jm-1 && is_land(h, i, j+1, k));
            if (neumann_jm1 || neumann_jp1) {
                const double tke_jm1_n = neumann_jm1 ? tke.x[i][j][k] : tke.x[i][j-1][k];
                const double tke_jp1_n = neumann_jp1 ? tke.x[i][j][k] : tke.x[i][j+1][k];
                const double dis_jm1_n = neumann_jm1 ? dis.x[i][j][k] : dis.x[i][j-1][k];
                const double dis_jp1_n = neumann_jp1 ? dis.x[i][j][k] : dis.x[i][j+1][k];
                dtkedthe   = (tke_jp1_n - tke_jm1_n) * inv_2dthe;
                d2tkedthe2 = (tke_jp1_n - 2.0*tke.x[i][j][k] + tke_jm1_n) * inv_dthe2;
                ddisdthe   = (dis_jp1_n - dis_jm1_n) * inv_2dthe;
                d2disdthe2 = (dis_jp1_n - 2.0*dis.x[i][j][k] + dis_jm1_n) * inv_dthe2;
            }
        }
    }


    // ===== PHI-direction derivatives =====
    if(k > 2 && k < km-3){
        if(land_ijk && is_air(h, i, j, k+1) && is_air(h, i, j, k+2)){
            #define COMPUTE_DPHI_B(FIELD, d1, d2) \
                d1 = (-3.0*FIELD.x[i][j][k] + 4.0*FIELD.x[i][j][k+1] - FIELD.x[i][j][k+2]) * inv_2dphi; \
                d2 = (FIELD.x[i][j][k] - 2.0*FIELD.x[i][j][k+1] + FIELD.x[i][j][k+2]) * inv_dphi2;

            COMPUTE_DPHI_B(u, dudphi, d2udphi2) COMPUTE_DPHI_B(v, dvdphi, d2vdphi2)
            COMPUTE_DPHI_B(w, dwdphi, d2wdphi2) COMPUTE_DPHI_B(t, dtdphi, d2tdphi2)
            COMPUTE_DPHI_B(c, dcdphi, d2cdphi2) COMPUTE_DPHI_B(cloud, dclouddphi, d2clouddphi2)
            COMPUTE_DPHI_B(ice, dicedphi, d2icedphi2) COMPUTE_DPHI_B(gr, dgdphi, d2gdphi2)
            COMPUTE_DPHI_B(co2, dcodphi, d2codphi2)
            COMPUTE_DPHI_B(tke, dtkedphi, d2tkedphi2) COMPUTE_DPHI_B(dis, ddisdphi, d2disdphi2)
            dpdphi   = (-3.0*p_dyn.x[i][j][k] + 4.0*p_dyn.x[i][j][k+1] - p_dyn.x[i][j][k+2]) * inv_2dphi;
            dnuedphi = (-3.0*nue.x[i][j][k]   + 4.0*nue.x[i][j][k+1]   - nue.x[i][j][k+2])   * inv_2dphi;
            phi_flag = true;
            #undef COMPUTE_DPHI_B
        }
        if(!phi_flag && land_ijk && is_air(h, i, j, k-1) && is_air(h, i, j, k-2)){
            #define COMPUTE_DPHI_C(FIELD, d1, d2) \
                d1 = -(-3.0*FIELD.x[i][j][k] + 4.0*FIELD.x[i][j][k-1] - FIELD.x[i][j][k-2]) * inv_2dphi; \
                d2 = -(FIELD.x[i][j][k] - 2.0*FIELD.x[i][j][k-1] + FIELD.x[i][j][k-2]) * inv_dphi2;

            COMPUTE_DPHI_C(u, dudphi, d2udphi2) COMPUTE_DPHI_C(v, dvdphi, d2vdphi2)
            COMPUTE_DPHI_C(w, dwdphi, d2wdphi2) COMPUTE_DPHI_C(t, dtdphi, d2tdphi2)
            COMPUTE_DPHI_C(c, dcdphi, d2cdphi2) COMPUTE_DPHI_C(cloud, dclouddphi, d2clouddphi2)
            COMPUTE_DPHI_C(ice, dicedphi, d2icedphi2) COMPUTE_DPHI_C(gr, dgdphi, d2gdphi2)
            COMPUTE_DPHI_C(co2, dcodphi, d2codphi2)
            COMPUTE_DPHI_C(tke, dtkedphi, d2tkedphi2) COMPUTE_DPHI_C(dis, ddisdphi, d2disdphi2)
            dpdphi   = -(-3.0*p_dyn.x[i][j][k] + 4.0*p_dyn.x[i][j][k-1] - p_dyn.x[i][j][k-2]) * inv_2dphi;
            dnuedphi = -(-3.0*nue.x[i][j][k]   + 4.0*nue.x[i][j][k-1]   - nue.x[i][j][k-2])   * inv_2dphi;
            phi_flag = true;
            #undef COMPUTE_DPHI_C
        }
    }
    else if(k == 0){
        #define COMPUTE_DPHI_W(FIELD, d1, d2) \
            d1 = (-3.0*FIELD.x[i][j][0] + 4.0*FIELD.x[i][j][1] - FIELD.x[i][j][2]) * inv_2dphi; \
            d2 = (FIELD.x[i][j][km-1] - 2.0*FIELD.x[i][j][km+2] + FIELD.x[i][j][km+3]) * inv_dphi2;

        COMPUTE_DPHI_W(u, dudphi, d2udphi2) COMPUTE_DPHI_W(v, dvdphi, d2vdphi2)
        COMPUTE_DPHI_W(w, dwdphi, d2wdphi2) COMPUTE_DPHI_W(t, dtdphi, d2tdphi2)
        COMPUTE_DPHI_W(c, dcdphi, d2cdphi2) COMPUTE_DPHI_W(cloud, dclouddphi, d2clouddphi2)
        COMPUTE_DPHI_W(ice, dicedphi, d2icedphi2) COMPUTE_DPHI_W(gr, dgdphi, d2gdphi2)
        COMPUTE_DPHI_W(co2, dcodphi, d2codphi2)
        COMPUTE_DPHI_W(tke, dtkedphi, d2tkedphi2) COMPUTE_DPHI_W(dis, ddisdphi, d2disdphi2)
        dpdphi   = (-3.0*p_dyn.x[i][j][0] + 4.0*p_dyn.x[i][j][1] - p_dyn.x[i][j][2]) * inv_2dphi;
        dnuedphi = (-3.0*nue.x[i][j][0]   + 4.0*nue.x[i][j][1]   - nue.x[i][j][2])   * inv_2dphi;
        phi_flag = true;
        #undef COMPUTE_DPHI_W
    }
    else if(k == km-1){
        #define COMPUTE_DPHI_E(FIELD, d1, d2) \
            d1 = -(-3.0*FIELD.x[i][j][km-1] + 4.0*FIELD.x[i][j][km-2] - FIELD.x[i][j][km-3]) * inv_2dphi; \
            d2 = -(FIELD.x[i][j][km-1] - 2.0*FIELD.x[i][j][km-2] + FIELD.x[i][j][km-3]) * inv_dphi2;

        COMPUTE_DPHI_E(u, dudphi, d2udphi2) COMPUTE_DPHI_E(v, dvdphi, d2vdphi2)
        COMPUTE_DPHI_E(w, dwdphi, d2wdphi2) COMPUTE_DPHI_E(t, dtdphi, d2tdphi2)
        COMPUTE_DPHI_E(c, dcdphi, d2cdphi2) COMPUTE_DPHI_E(cloud, dclouddphi, d2clouddphi2)
        COMPUTE_DPHI_E(ice, dicedphi, d2icedphi2) COMPUTE_DPHI_E(gr, dgdphi, d2gdphi2)
        COMPUTE_DPHI_E(co2, dcodphi, d2codphi2)
        COMPUTE_DPHI_E(tke, dtkedphi, d2tkedphi2) COMPUTE_DPHI_E(dis, ddisdphi, d2disdphi2)
        dpdphi   = -(-3.0*p_dyn.x[i][j][km-1] + 4.0*p_dyn.x[i][j][km-2] - p_dyn.x[i][j][km-3]) * inv_2dphi;
        dnuedphi = -(-3.0*nue.x[i][j][km-1]   + 4.0*nue.x[i][j][km-2]   - nue.x[i][j][km-3])   * inv_2dphi;
        phi_flag = true;
        #undef COMPUTE_DPHI_E
    }

    if(!phi_flag){
        #define COMPUTE_DPHI_A(FIELD, d1, d2) \
            d1 = (FIELD.x[i][j][k+1] - FIELD.x[i][j][k-1]) * inv_2dphi; \
            d2 = (FIELD.x[i][j][k+1] - 2.0*FIELD.x[i][j][k] + FIELD.x[i][j][k-1]) * inv_dphi2;

        COMPUTE_DPHI_A(u, dudphi, d2udphi2) COMPUTE_DPHI_A(v, dvdphi, d2vdphi2)
        COMPUTE_DPHI_A(w, dwdphi, d2wdphi2) COMPUTE_DPHI_A(t, dtdphi, d2tdphi2)
        COMPUTE_DPHI_A(c, dcdphi, d2cdphi2) COMPUTE_DPHI_A(cloud, dclouddphi, d2clouddphi2)
        COMPUTE_DPHI_A(ice, dicedphi, d2icedphi2) COMPUTE_DPHI_A(gr, dgdphi, d2gdphi2)
        COMPUTE_DPHI_A(co2, dcodphi, d2codphi2)
        COMPUTE_DPHI_A(tke, dtkedphi, d2tkedphi2) COMPUTE_DPHI_A(dis, ddisdphi, d2disdphi2)
        dpdphi   = (p_dyn.x[i][j][k+1] - p_dyn.x[i][j][k-1]) * inv_2dphi;
        dnuedphi = (nue.x[i][j][k+1]   - nue.x[i][j][k-1])   * inv_2dphi;
        #undef COMPUTE_DPHI_A

        // Neumann BC for tke/dis at phi-direction land boundaries.
        // Without this, an ocean cell adjacent to land reads tke/dis=0 from the
        // zeroed land cell — the 1/(r sinθ)² amplification then makes D_w
        // enormous near the poles, causing blow-up within a few iterations.
        if (!land_ijk) {
            const bool neumann_km1 = (k > 0    && is_land(h, i, j, k-1));
            const bool neumann_kp1 = (k < km-1 && is_land(h, i, j, k+1));
            if (neumann_km1 || neumann_kp1) {
                const double tke_km1_n = neumann_km1 ? tke.x[i][j][k] : tke.x[i][j][k-1];
                const double tke_kp1_n = neumann_kp1 ? tke.x[i][j][k] : tke.x[i][j][k+1];
                const double dis_km1_n = neumann_km1 ? dis.x[i][j][k] : dis.x[i][j][k-1];
                const double dis_kp1_n = neumann_kp1 ? dis.x[i][j][k] : dis.x[i][j][k+1];
                dtkedphi   = (tke_kp1_n - tke_km1_n) * inv_2dphi;
                d2tkedphi2 = (tke_kp1_n - 2.0*tke.x[i][j][k] + tke_km1_n) * inv_dphi2;
                ddisdphi   = (dis_kp1_n - dis_km1_n) * inv_2dphi;
                d2disdphi2 = (dis_kp1_n - 2.0*dis.x[i][j][k] + dis_km1_n) * inv_dphi2;
            }
        }
    }


    // ===== Physics: coefficients =====
    double buoyancy = 1.0;

    double coeff_buoy     = r_air * g;
    double coeff_energy_p = u_0 * u_0 / (cp_l * t_0);
    double coeff_energy   = 1.0 / (c_0 * cp_l * t_0);
    double coeff_u_p      = 0.5 * r_air * u_0 * u_0 / L_atm;
    double coeff_trans    = 1.0;
    double coeff_MC_vel   = L_atm / (u_0 * u_0);
    double coeff_MC_q     = L_atm / (u_0 * c_0);
    double coeff_MC_t     = L_atm / (u_0 * t_0);

    // Coriolis
    double coriolis = 1.0;
    double coeff_Coriolis = r_air * u_0 * omega;
    double coriolis_rad =  2.0 * costhe * w_ijk;
    double coriolis_the = -2.0 * sinthe * w_ijk;
    double coriolis_phi =  2.0 * (sinthe * v_ijk - costhe * u_ijk);

    // Acting forces
    CoriolisForce.x[i][j][k] = coriolis * coeff_Coriolis
        * sqrt((pow(coriolis_rad, 2)
        + pow(coriolis_the, 2)
        + pow(coriolis_phi, 2)) / 3.0);
    BuoyancyForce.x[i][j][k] = buoyancy * coeff_buoy * (t.x[i][j][k] - 1.0);
    PressureGradientForce.x[i][j][k] = -coeff_u_p
        * sqrt((pow(dpdr, 2)
        + pow(dpdthe * inv_rm, 2)
        + pow(dpdphi * inv_rmsinthe, 2)) / 3.0);

    double t_u    = t.x[i][j][k] * t_0;
    double E_Rain = hp * exp_func(t_u, 17.2694, 35.86);
    double E_Ice  = hp * exp_func(t_u, 21.8746, 7.66);
    double q_Rain = ep * E_Rain / (p_hydro.x[i][j][k] - E_Rain);
    double q_Ice  = ep * E_Ice  / (p_hydro.x[i][j][k] - E_Ice);
    double Q_Latent_Ice = 0.0;

    double coeff_S = lamda * t_0 / L_atm;
    double coeff_L = r_air * c_0 * lv * u_0 / L_atm;

    if(c.x[i][j][k] >= 0.85 * q_Rain){
        Q_Latent.x[i][j][k] = u_ijk * dclouddr
            + v_ijk * inv_rm     * dclouddthe
            + w_ijk * inv_rmsinthe * dclouddphi;
    } else Q_Latent.x[i][j][k] = 0.0;

    if(c.x[i][j][k] >= 0.85 * q_Ice){
        Q_Latent_Ice = u_ijk * dicedr
            + v_ijk * inv_rm     * dicedthe
            + w_ijk * inv_rmsinthe * dicedphi;
    } else Q_Latent_Ice = 0.0;

    Q_Latent.x[i][j][k] = coeff_L * (Q_Latent.x[i][j][k] + Q_Latent_Ice);

    Q_Sensible.x[i][j][k] = coeff_S
        * (d2tdr2 + dtdr * 2.0 * inv_rm + d2tdthe2 * inv_rm2
        + dtdthe * costhe_inv_rm2sinthe + d2tdphi2 * inv_rm2sinthe2);

    if(land_ijk){
        BuoyancyForce.x[i][j][k]        = 0.0;
        PressureGradientForce.x[i][j][k] = 0.0;
        CoriolisForce.x[i][j][k]         = 0.0;
        Q_Latent.x[i][j][k]              = 0.0;
        Q_Sensible.x[i][j][k]            = 0.0;
    }


    // ===== Turbulence model: diffusion coefficients =====
    double diffusion_t_re   = 0.0;
    double diffusion_vel_re = 0.0;
    double diffusion_tke_re = 0.0;
    double diffusion_dis_re = 0.0;
    const double diff_prec_re_inv = 1.0 / (sc_WaterVapour * re_turb);
    const double diff_co2_re_inv  = 1.0 / (sc_CO2 * re_turb);

    double d_i       = 0.0;
    double d_i_mount = 0.0;
    int i_mount = i_topography[j][k];
    if(use_stretched_coordinate_system){
        d_i       = log(zeta * (rad.z[i]       + 1.0));
        d_i_mount = log(zeta * (rad.z[i_mount] + 1.0));
    } else {
        d_i       = (double)i;
        d_i_mount = (double)i_mount;
    }

    // ---- Refresh prod.x at every RK4 sub-stage ----
    // compute_sources() updates prod.x only at even iterations, so the stored
    // value is stale for intermediate RK4 stages.  Recompute it here from the
    // velocity derivatives already available above, scaled to match the
    // spherical-coord convention in compute_sources() (exp_rm / inv_rm /
    // inv_rmsinthe factors).  Constants mirror TurbulenceAtm private members.
    if(use_turbulence_model && !land_ijk){
        constexpr double C_nue_loc   = 0.028;    // mirrors TurbulenceAtm::C_nue
        constexpr double dis_min_loc = 1.0e-10;  // mirrors TurbulenceAtm::dis_min
        const double dudr_s   = dudr   * exp_rm;
        const double dvdr_s   = dvdr   * exp_rm;
        const double dwdr_s   = dwdr   * exp_rm;
        const double dudthe_s = dudthe * inv_rm;
        const double dvdthe_s = dvdthe * inv_rm;
        const double dwdthe_s = dwdthe * inv_rm;
        const double dudphi_s = dudphi * inv_rmsinthe;
        const double dvdphi_s = dvdphi * inv_rmsinthe;
        const double dwdphi_s = dwdphi * inv_rmsinthe;

        const double dis_here = std::max(dis.x[i][j][k], dis_min_loc);
        const double tke_here = std::max(tke.x[i][j][k], 0.0);
        const double nue_max  = 1000.0 / (u_0 * L_atm);
        double cnue = use_k_epsilon_turbulence_model
                    ? C_nue_loc * tke_here * tke_here / dis_here
                    : tke_here / dis_here;
        cnue = std::min(cnue, nue_max);

        const double der = 0.66667 * (dudr_s + dvdthe_s + dwdphi_s);

        prod.x[i][j][k] = std::max(0.0,
              (cnue * (2.0 * dudr_s   - der) - 0.66667 * tke.x[i][j][k]) * dudr_s
            + (cnue * (dudthe_s + dvdr_s))                               * dudthe_s
            + (cnue * (dudphi_s + dwdr_s))                               * dudphi_s
            + (cnue * (2.0 * dvdthe_s - der) - 0.66667 * tke.x[i][j][k]) * dvdthe_s
            + (cnue * (dvdr_s + dudthe_s))                               * dvdr_s
            + (cnue * (dvdphi_s + dwdthe_s))                             * dvdphi_s
            + (cnue * (2.0 * dwdphi_s - der) - 0.66667 * tke.x[i][j][k]) * dwdphi_s
            + (cnue * (dwdr_s + dudphi_s))                               * dwdr_s
            + (cnue * (dwdthe_s + dvdphi_s))                             * dwdthe_s);
    }

    if(use_k_epsilon_turbulence_model){                                 // KE_chien 1982
        const double sig_k   = 1.0;
        const double sig_w   = 1.3;
        const double C_eps_1 = 1.35;
        const double C_eps_2 = 1.80;
        // Diffusion coefficients: molecular background 1/re_turb plus turbulent nue*.
        // nue* = ν_T/(u_0·L_atm) enters directly — NOT divided by re_turb — because
        // the non-dimensionalisation already bakes the (u_0·L_atm) scale into nue*.
        // pr_turb (turbulent Prandtl number, configurable) replaces the laminar pr for
        // both the molecular background and the turbulent contribution to heat transport.
        {
            // Diffusion: molecular background 1/re_turb plus turbulent nue* directly.
            // nue* = ν_T/(u_0·L_atm) is already in the same non-dim units as the
            // diffusion coefficient — it must NOT be divided by re_turb again.
            diffusion_t_re   = (1.0/re_turb + nue.x[i][j][k]) / pr_turb;
            diffusion_vel_re = 1.0/re_turb + nue.x[i][j][k];
            diffusion_tke_re = 1.0/re_turb + nue.x[i][j][k] / sig_k;
            diffusion_dis_re = 1.0/re_turb + nue.x[i][j][k] / sig_w;
        }
        // Source terms — computed fresh from current tke/dis so every RK4 sub-stage
        // is consistent (compute_k_epsilon runs only at even iterations).
        // L_k / L_w wall-damping corrections (Chien 1982) are kept from the stored
        // tke_source / dis_source set by compute_k_epsilon; only the P–Y balance is
        // updated here so the damping is not lost.
        {
            const double tke_s = std::max(tke.x[i][j][k], 1.0e-10);
            const double dis_s = std::max(dis.x[i][j][k], 1.0e-10);
            const double P_k   = prod.x[i][j][k];
            const double Y_k   = dis_s;
            const double P_w   = C_eps_1 * (dis_s / tke_s) * P_k;
            const double Y_w   = C_eps_2 * dis_s * dis_s / tke_s;
            const double dis_src_max = 20.0 * dis_s * dis_s;
            const double tke_src_max = 20.0 * tke_s * dis_s;
            tke_source.x[i][j][k] = std::clamp(P_k - Y_k, -tke_src_max, tke_src_max);
            dis_source.x[i][j][k] = std::clamp(P_w - Y_w, -dis_src_max, dis_src_max);
        }
    }

    if(use_k_omega_turbulence_model){                                   // Wilcox 1988
        const double sig_k    = 0.6;
        const double sig_w    = 0.5;
        const double bet_star = 0.09;
        const double gam      = 0.52;
        const double bet_0    = 0.0708;
        {
            diffusion_t_re   = (1.0/re_turb + nue.x[i][j][k]) / pr_turb;
            diffusion_vel_re = 1.0/re_turb + nue.x[i][j][k];
            diffusion_tke_re = 1.0/re_turb + sig_k * nue.x[i][j][k];
            diffusion_dis_re = 1.0/re_turb + sig_w * nue.x[i][j][k];
        }
        // Source terms — fresh tke/dis, stored prod
        {
            const double tke_s  = std::max(tke.x[i][j][k], 1.0e-10);
            const double dis_s  = std::max(dis.x[i][j][k], 1.0e-10);
            const double P_k    = std::min(prod.x[i][j][k], 20.0 * bet_star * tke_s * dis_s);
            const double Y_k    = bet_star * tke_s * dis_s;
            const double P_w    = gam * (dis_s / tke_s) * P_k;
            const double Y_w    = bet_0 * dis_s * dis_s;
            const double cross  = dtkedr   * ddisdr
                                + dtkedthe * inv_rm       * ddisdthe  * inv_rm
                                + dtkedphi * inv_rmsinthe * ddisdphi  * inv_rmsinthe;
            const double sig_d  = (cross <= 0.0) ? 0.0 : 0.125;
            const double D_w    = sig_d * cross / dis_s;
            const double dis_src_max = 20.0 * dis_s * dis_s;
            const double tke_src_max = 20.0 * tke_s * dis_s;
            tke_source.x[i][j][k] = std::clamp(P_k - Y_k, -tke_src_max, tke_src_max);
            dis_source.x[i][j][k] = std::clamp(P_w - Y_w + D_w, -dis_src_max, dis_src_max);
        }
    }

    if(use_k_omega_SST_turbulence_model){                               // Menter 1994
        double bet_star = 0.09;    // β* destruction coefficient — NOT C_μ=0.028
        double sig_k1   = 1.176;
        double sig_k2   = 1.0;
        double sig_w1   = 2.0;
        double sig_w2   = 1.168;

        int j_beg = 0, k_beg = 0;
        double d_j_mount = 0.0, d_k_mount = 0.0;
        double d_j_south = 0.0, d_j_north = 0.0;
        double d_k_east  = 0.0, d_k_west  = 0.0;
        double d_j = (double)j;
        double d_k = (double)k;

        if(land_ijk && is_water(h, i, j+1, k))   j_beg = j;
        for(int jj = j_beg; jj < jm-1; jj++)
            if(is_land(h, i, jj, k) || jj == jm-1) d_j_mount = jj;
        d_j_south = fabs(d_j - d_j_mount);

        if(land_ijk && is_water(h, i, j-1, k))   j_beg = j;
        for(int jj = jm-1; jj > j_beg; jj--)
            if(is_land(h, i, jj, k) || jj == 0) d_j_mount = jj;
        d_j_north = fabs(d_j - d_j_mount);

        if(land_ijk && is_water(h, i, j, k+1))   k_beg = k;
        for(int kk = k_beg; kk < km-1; kk++)
            if(is_land(h, i, j, kk) || kk == km-1) d_k_mount = kk;
        d_k_east = fabs(d_k - d_k_mount);

        if(land_ijk && is_water(h, i, j, k-1))   k_beg = k;
        for(int kk = km-1; kk > k_beg; kk--)
            if(is_land(h, i, j, kk) || kk == 0) d_k_mount = kk;
        d_k_west = fabs(d_k - d_k_mount);

        double d_close = d_i - d_i_mount;
        if(d_j_south <= d_close) d_close = d_j_south;
        if(d_j_north <= d_close) d_close = d_j_north;
        if(d_k_east  <= d_close) d_close = d_k_east;
        if(d_k_west  <= d_close) d_close = d_k_west;

        double CD_kw = std::max(
            (2.0 * sig_w2) / dis.x[i][j][k]
            * (dtkedr   * ddisdr
            + dtkedthe * inv_rm      * ddisdthe  * inv_rm
            + dtkedphi * inv_rmsinthe * ddisdphi * inv_rmsinthe),
            1.0e-20);

        // Menter (1994): arg1 uses molecular viscosity ν_air, not turbulent ν_T.
        // nue_air_nd = ν_air / (u_0·L_atm) is the non-dimensional molecular viscosity
        // (same constant as TurbulenceAtm uses for nue_air_nd).
        constexpr double nue_air_phys = 1.5e-5;   // [m²/s]
        const double nue_air_nd = nue_air_phys / (u_0 * L_atm);
        double arg1 = std::min(std::max(
            sqrt(tke.x[i][j][k]) / (bet_star * dis.x[i][j][k] * d_close),
            (500.0 * nue_air_nd) / (d_close * d_close * dis.x[i][j][k])),
            (4.0 * sig_w2 * tke.x[i][j][k]) / (CD_kw * d_close * d_close));
        double F1 = tanh(pow(arg1, 4));

        {
            diffusion_t_re   = (1.0/re_turb + nue.x[i][j][k]) / pr_turb;
            diffusion_vel_re =  1.0/re_turb + nue.x[i][j][k];
            diffusion_tke_re =  1.0/re_turb + nue.x[i][j][k] / blend(sig_k1, sig_k2, F1);
            diffusion_dis_re =  1.0/re_turb + nue.x[i][j][k] / blend(sig_w1, sig_w2, F1);
        }

        // Source terms — computed fresh from current tke/dis at every RK4 sub-stage.
        // compute_k_omega_SST runs only at even iterations and would otherwise leave
        // tke_source/dis_source one outer iteration behind.
        {
            const double bet1   = 0.0333;
            const double bet2   = 0.0368;
            const double gam1   = 0.413;
            const double gam2   = 0.2;
            const double tke_s  = std::max(tke.x[i][j][k], 1.0e-10);
            const double dis_s  = std::max(dis.x[i][j][k], 1.0e-10);
            const double P_k    = std::min(prod.x[i][j][k], 20.0 * bet_star * tke_s * dis_s);
            const double Y_k    = bet_star * tke_s * dis_s;
            const double P_w    = blend(gam1, gam2, F1) * P_k * dis_s / tke_s;
            const double Y_w    = blend(bet1,  bet2, F1) * dis_s * dis_s;
            const double D_w    = 2.0 * (1.0 - F1) * sig_w2 / dis_s
                                * (dtkedr   * ddisdr
                                + dtkedthe * inv_rm       * ddisdthe  * inv_rm
                                + dtkedphi * inv_rmsinthe * ddisdphi  * inv_rmsinthe);
            const double dis_src_max = 20.0 * dis_s * dis_s;
            const double tke_src_max = 20.0 * tke_s * dis_s;
            tke_source.x[i][j][k] = std::clamp(P_k - Y_k, -tke_src_max, tke_src_max);
            dis_source.x[i][j][k] = std::clamp(P_w - Y_w + D_w, -dis_src_max, dis_src_max);
        }
    }

    if(!use_turbulence_model){
        diffusion_t_re        = 1.0 / (re_turb * pr_turb);
        diffusion_vel_re      = 1.0 / re_turb;
        diffusion_tke_re      = 0.0;
        diffusion_dis_re      = 0.0;
        tke.x[i][j][k]        = 0.0;
        dis.x[i][j][k]        = 0.0;
        nue.x[i][j][k]        = 0.0;
        prod.x[i][j][k]       = 0.0;
        tke_source.x[i][j][k] = 0.0;
        dis_source.x[i][j][k] = 0.0;
        tken.x[i][j][k]       = 0.0;
        disn.x[i][j][k]       = 0.0;
    }


    // ===== Transport terms (advection) =====
    double u_exp   = u_ijk * exp_rm;
    double v_invrm = v_ijk * inv_rm;
    double w_invrs = w_ijk * inv_rmsinthe;

    double pressure_t = coeff_energy_p * (u_exp * dpdr + v_invrm * dpdthe + w_invrs * dpdphi);

    double transport_t     = u_exp * dtdr     + v_invrm * dtdthe     + w_invrs * dtdphi;
    double transport_u     = u_exp * dudr     + v_invrm * dudthe     + w_invrs * dudphi;
    double transport_v     = u_exp * dvdr     + v_invrm * dvdthe     + w_invrs * dvdphi;
    double transport_w     = u_exp * dwdr     + v_invrm * dwdthe     + w_invrs * dwdphi;
    double transport_c     = u_exp * dcdr     + v_invrm * dcdthe     + w_invrs * dcdphi;
    double transport_cloud = u_exp * dclouddr + v_invrm * dclouddthe + w_invrs * dclouddphi;
    double transport_ice   = u_exp * dicedr   + v_invrm * dicedthe   + w_invrs * dicedphi;
    double transport_g     = u_exp * dgdr     + v_invrm * dgdthe     + w_invrs * dgdphi;
    double transport_co2   = u_exp * dcodr    + v_invrm * dcodthe    + w_invrs * dcodphi;
    double transport_tke   = u_exp * dtkedr   + v_invrm * dtkedthe   + w_invrs * dtkedphi;
    double transport_dis   = u_exp * ddisdr   + v_invrm * ddisdthe   + w_invrs * ddisdphi;

    // Divergence correction for tke/dis transport
    double div_vel = dudr + dvdthe * inv_rm + dwdphi * inv_rmsinthe;
    transport_tke += tke.x[i][j][k] * div_vel;
    transport_dis += dis.x[i][j][k] * div_vel;


    // ===== Diffusion terms =====
    double two_over_rm_exp  = 2.0 * inv_rm * exp_rm;
    double cos_rm2sin       = costhe_inv_rm2sinthe;
    double v_metric         = (1.0 + costhe / sinthe2) * inv_rm2;

    double diffusion_t = (d2tdr2 * exp_2_rm + dtdr * two_over_rm_exp
        + d2tdthe2 * inv_rm2 + dtdthe * cos_rm2sin
        + d2tdphi2 * inv_rm2sinthe2) * diffusion_t_re;

    double diffusion_u = (d2udr2 * exp_2_rm + 2.0 * u_ijk * inv_rm2
        + d2udthe2 * inv_rm2 + 4.0 * dudr * inv_rm * exp_rm
        + dudthe * cos_rm2sin + d2udphi2 * inv_rm2sinthe2) * diffusion_vel_re;

    double diffusion_v = (d2vdr2 * exp_2_rm + dvdr * two_over_rm_exp
        + d2vdthe2 * inv_rm2 + dvdthe * cos_rm2sin
        - v_metric * v_ijk + d2vdphi2 * inv_rm2sinthe2
        + 2.0 * dudthe * inv_rm2
        - dwdphi * 2.0 * costhe * inv_rm2sinthe2) * diffusion_vel_re;

    double diffusion_w = (d2wdr2 * exp_2_rm + dwdr * two_over_rm_exp
        + d2wdthe2 * inv_rm2 + dwdthe * cos_rm2sin
        - v_metric * w_ijk + d2wdphi2 * inv_rm2sinthe2
        + 2.0 * dudphi * inv_rm2sinthe
        + dvdphi * 2.0 * costhe * inv_rm2sinthe2) * diffusion_vel_re;

    double diffusion_c = (d2cdr2 * exp_2_rm + dcdr * two_over_rm_exp
        + d2cdthe2 * inv_rm2 + dcdthe * cos_rm2sin
        + d2cdphi2 * inv_rm2sinthe2) * diff_prec_re_inv;

    double diffusion_cloud = (d2clouddr2 * exp_2_rm + dclouddr * two_over_rm_exp
        + d2clouddthe2 * inv_rm2 + dclouddthe * cos_rm2sin
        + d2clouddphi2 * inv_rm2sinthe2) * diff_prec_re_inv;

    double diffusion_ice = (d2icedr2 * exp_2_rm + dicedr * two_over_rm_exp
        + d2icedthe2 * inv_rm2 + dicedthe * cos_rm2sin
        + d2icedphi2 * inv_rm2sinthe2) * diff_prec_re_inv;

    double diffusion_g = (d2gdr2 * exp_2_rm + dgdr * two_over_rm_exp
        + d2gdthe2 * inv_rm2 + dgdthe * cos_rm2sin
        + d2gdphi2 * inv_rm2sinthe2) * diff_prec_re_inv;

    double diffusion_co2 = (d2codr2 * exp_2_rm + dcodr * two_over_rm_exp
        + d2codthe2 * inv_rm2 + dcodthe * cos_rm2sin
        + d2codphi2 * inv_rm2sinthe2) * diff_co2_re_inv;

    double diffusion_tke = (d2tkedr2 * exp_2_rm + dtkedr * two_over_rm_exp
        + d2tkedthe2 * inv_rm2 + dtkedthe * cos_rm2sin
        + d2tkedphi2 * inv_rm2sinthe2) * diffusion_tke_re;

    double diffusion_dis = (d2disdr2 * exp_2_rm + ddisdr * two_over_rm_exp
        + d2disdthe2 * inv_rm2 + ddisdthe * cos_rm2sin
        + d2disdphi2 * inv_rm2sinthe2) * diffusion_dis_re;


    // ===== RHS assembly =====
    double dpdr_exp      = dpdr   * exp_rm;
    double dpdthe_invrm  = dpdthe * inv_rm;
    double dpdphi_invrs  = dpdphi * inv_rmsinthe;

    rhs_t.x[i][j][k] = pressure_t - transport_t + diffusion_t
        + coeff_MC_t * MC_t.x[i][j][k]
        + coeff_energy * (S_c.x[i][j][k] + S_r.x[i][j][k]) * lv * r_humid.x[i][j][k]
        + coeff_energy * (S_i.x[i][j][k] + S_s.x[i][j][k] + S_g.x[i][j][k]) * ls * r_humid.x[i][j][k]
        + Q_Latent.x[i][j][k] / coeff_L;

    rhs_u.x[i][j][k] = -dpdr_exp - transport_u + diffusion_u
        + buoyancy * (t.x[i][j][k] - 1.0)
        + coriolis * coriolis_rad;

    rhs_v.x[i][j][k] = -dpdthe_invrm - transport_v + diffusion_v
        + coriolis * coriolis_the
        + coeff_MC_vel * MC_v.x[i][j][k];

    rhs_w.x[i][j][k] = -dpdphi_invrs - transport_w + diffusion_w
        + coriolis * coriolis_phi
        + coeff_MC_vel * MC_w.x[i][j][k];

    rhs_c.x[i][j][k] = -transport_c + diffusion_c
        + coeff_trans * S_v.x[i][j][k] * r_humid.x[i][j][k]
        + coeff_MC_q * MC_q.x[i][j][k];

    rhs_cloud.x[i][j][k] = -transport_cloud + diffusion_cloud
        + coeff_trans * S_c.x[i][j][k] * r_humid.x[i][j][k];

    rhs_ice.x[i][j][k] = -transport_ice + diffusion_ice
        + coeff_trans * S_i.x[i][j][k] * r_humid.x[i][j][k];

    rhs_g.x[i][j][k] = -transport_g + diffusion_g
        + coeff_trans * S_g.x[i][j][k] * r_humid.x[i][j][k];

    rhs_co2.x[i][j][k] = -transport_co2 + diffusion_co2;

    rhs_tke.x[i][j][k] = -transport_tke + diffusion_tke + tke_source.x[i][j][k];
    rhs_dis.x[i][j][k] = -transport_dis + diffusion_dis + dis_source.x[i][j][k];

    aux_u.x[i][j][k] = rhs_u.x[i][j][k] + dpdr_exp;
    aux_v.x[i][j][k] = rhs_v.x[i][j][k] + dpdthe_invrm;
    aux_w.x[i][j][k] = rhs_w.x[i][j][k] + dpdphi_invrs;
}
