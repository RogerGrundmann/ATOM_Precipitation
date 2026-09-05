/*
 * Ocean General Circulation Model (OGCM) applied to turbulent flow
 * RHS for 7 prognostic fields: t, u, v, w, c (salinity), tke, dis
 * Called by solveRungeKutta_Hydrosphere_Turb at every RK4 sub-stage.
*/

#include <iostream>
#include <cmath>
#include <cstdlib>
#include "cHydrosphereModel.h"
#include "HydHorizViscosity.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;


void cHydrosphereModel::RHS_Hydrosphere_Turb(int i, int j, int k, const CellGeometry& geo){

    const double exp_rm              = geo.exp_rm;
    const double exp_2_rm            = geo.exp_2_rm;
    const double sinthe              = geo.sinthe;
    const double sinthe2             = geo.sinthe2;
    const double costhe              = geo.costhe;
    const double inv_rm              = geo.inv_rm;
    const double inv_rm2             = geo.inv_rm2;
    const double inv_rmsinthe        = geo.inv_rmsinthe;
    const double inv_rm2sinthe       = geo.inv_rm2sinthe;
    const double inv_rm2sinthe2      = geo.inv_rm2sinthe2;
    const double costhe_inv_rm2sinthe = geo.costhe_inv_rm2sinthe;

    const double inv_2dthe = geo.inv_2dthe;
    const double inv_2dphi = geo.inv_2dphi;
    const double inv_dthe2 = geo.inv_dthe2;
    const double inv_dphi2 = geo.inv_dphi2;
    // radial derivatives use the per-i non-uniform stencil coeffs (rc*/rf*), not inv_2dr/inv_dr2

    double u_ijk = u.x[i][j][k];
    double v_ijk = v.x[i][j][k];
    double w_ijk = w.x[i][j][k];

    bool land_ijk = is_land(h, i, j, k);

    bool r_flag = false, the_flag = false, phi_flag = false;

    // ---- First-order derivatives ----
    double dudr, dvdr, dwdr, dtdr, dpdr, dcdr, dtkedr, ddisdr, dnuedr;
    double dudthe, dvdthe, dwdthe, dtdthe, dpdthe, dcdthe, dtkedthe, ddisdthe, dnuedthe;
    double dudphi, dvdphi, dwdphi, dtdphi, dpdphi, dcdphi, dtkedphi, ddisdphi, dnuedphi;

    // ---- Second-order derivatives ----
    double d2udr2, d2vdr2, d2wdr2, d2tdr2, d2cdr2, d2tkedr2, d2disdr2;
    double d2udthe2, d2vdthe2, d2wdthe2, d2tdthe2, d2cdthe2, d2tkedthe2, d2disdthe2;
    double d2udphi2, d2vdphi2, d2wdphi2, d2tdphi2, d2cdphi2, d2tkedphi2, d2disdphi2;


    // ===== R-direction derivatives =====
    if(i < im-2 && land_ijk && is_air(h, i+1, j, k)){
        // Forward-biased stencil (points i, i+1, i+2) on the stretched radial grid
        #define COMPUTE_DR_B(FIELD, d1, d2) \
            d1 = rf10[i]*FIELD.x[i][j][k] + rf11[i]*FIELD.x[i+1][j][k] + rf12[i]*FIELD.x[i+2][j][k]; \
            d2 = rf20[i]*FIELD.x[i][j][k] + rf21[i]*FIELD.x[i+1][j][k] + rf22[i]*FIELD.x[i+2][j][k];
        COMPUTE_DR_B(u, dudr, d2udr2)  COMPUTE_DR_B(v, dvdr, d2vdr2)
        COMPUTE_DR_B(w, dwdr, d2wdr2)  COMPUTE_DR_B(t, dtdr, d2tdr2)
        COMPUTE_DR_B(c, dcdr, d2cdr2)
        COMPUTE_DR_B(tke, dtkedr, d2tkedr2)  COMPUTE_DR_B(dis, ddisdr, d2disdr2)
        dpdr   = rf10[i]*p_dyn.x[i][j][k] + rf11[i]*p_dyn.x[i+1][j][k] + rf12[i]*p_dyn.x[i+2][j][k];
        dnuedr = rf10[i]*nue.x[i][j][k]   + rf11[i]*nue.x[i+1][j][k]   + rf12[i]*nue.x[i+2][j][k];
        r_flag = true;
        #undef COMPUTE_DR_B
    }
    else if(i == 0){
        // Forward stencil at the bottom (points 0, 1, 2) on the stretched radial grid
        #define COMPUTE_DR_F(FIELD, d1, d2) \
            d1 = rf10[0]*FIELD.x[0][j][k] + rf11[0]*FIELD.x[1][j][k] + rf12[0]*FIELD.x[2][j][k]; \
            d2 = rf20[0]*FIELD.x[0][j][k] + rf21[0]*FIELD.x[1][j][k] + rf22[0]*FIELD.x[2][j][k];
        COMPUTE_DR_F(u, dudr, d2udr2)  COMPUTE_DR_F(v, dvdr, d2vdr2)
        COMPUTE_DR_F(w, dwdr, d2wdr2)  COMPUTE_DR_F(t, dtdr, d2tdr2)
        COMPUTE_DR_F(c, dcdr, d2cdr2)
        COMPUTE_DR_F(tke, dtkedr, d2tkedr2)  COMPUTE_DR_F(dis, ddisdr, d2disdr2)
        dpdr   = rf10[0]*p_dyn.x[0][j][k] + rf11[0]*p_dyn.x[1][j][k] + rf12[0]*p_dyn.x[2][j][k];
        dnuedr = rf10[0]*nue.x[0][j][k]   + rf11[0]*nue.x[1][j][k]   + rf12[0]*nue.x[2][j][k];
        r_flag = true;
        #undef COMPUTE_DR_F
    }

    if(!r_flag){
        // Central stencil (points i-1, i, i+1) on the stretched radial grid
        #define COMPUTE_DR_C(FIELD, d1, d2) \
            d1 = rc1m[i]*FIELD.x[i-1][j][k] + rc10[i]*FIELD.x[i][j][k] + rc1p[i]*FIELD.x[i+1][j][k]; \
            d2 = rc2m[i]*FIELD.x[i-1][j][k] + rc20[i]*FIELD.x[i][j][k] + rc2p[i]*FIELD.x[i+1][j][k];
        COMPUTE_DR_C(u, dudr, d2udr2)  COMPUTE_DR_C(v, dvdr, d2vdr2)
        COMPUTE_DR_C(w, dwdr, d2wdr2)  COMPUTE_DR_C(t, dtdr, d2tdr2)
        COMPUTE_DR_C(c, dcdr, d2cdr2)
        COMPUTE_DR_C(tke, dtkedr, d2tkedr2)  COMPUTE_DR_C(dis, ddisdr, d2disdr2)
        dpdr   = rc1m[i]*p_dyn.x[i-1][j][k] + rc10[i]*p_dyn.x[i][j][k] + rc1p[i]*p_dyn.x[i+1][j][k];
        dnuedr = rc1m[i]*nue.x[i-1][j][k]   + rc10[i]*nue.x[i][j][k]   + rc1p[i]*nue.x[i+1][j][k];
        #undef COMPUTE_DR_C

        // Neumann BC at seafloor: first water cell reads land (zeroed) below it.
        // Replace land neighbour's tke/dis with current cell value (zero-gradient),
        // using the local upward spacing hp = rad.z[i+1]-rad.z[i] (stretched grid).
        if (!land_ijk && i > 0 && is_land(h, i-1, j, k)) {
            const double hp = rad.z[i+1] - rad.z[i];
            const double inv_hp = 1.0 / hp;
            const double inv_hp2 = inv_hp * inv_hp;
            dtkedr   = (tke.x[i+1][j][k] - tke.x[i][j][k]) * inv_hp;
            d2tkedr2 = (tke.x[i+1][j][k] - tke.x[i][j][k]) * inv_hp2;
            ddisdr   = (dis.x[i+1][j][k] - dis.x[i][j][k]) * inv_hp;
            d2disdr2 = (dis.x[i+1][j][k] - dis.x[i][j][k]) * inv_hp2;
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
            COMPUTE_DTHE_B(c, dcdthe, d2cdthe2)
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
            COMPUTE_DTHE_C(c, dcdthe, d2cdthe2)
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
        COMPUTE_DTHE_NP(c, dcdthe, d2cdthe2)
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
        COMPUTE_DTHE_SP(c, dcdthe, d2cdthe2)
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
        COMPUTE_DTHE_A(c, dcdthe, d2cdthe2)
        COMPUTE_DTHE_A(tke, dtkedthe, d2tkedthe2) COMPUTE_DTHE_A(dis, ddisdthe, d2disdthe2)
        dpdthe   = (p_dyn.x[i][j+1][k] - p_dyn.x[i][j-1][k]) * inv_2dthe;
        dnuedthe = (nue.x[i][j+1][k]   - nue.x[i][j-1][k])   * inv_2dthe;
        #undef COMPUTE_DTHE_A

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
            COMPUTE_DPHI_B(c, dcdphi, d2cdphi2)
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
            COMPUTE_DPHI_C(c, dcdphi, d2cdphi2)
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
        COMPUTE_DPHI_W(c, dcdphi, d2cdphi2)
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
        COMPUTE_DPHI_E(c, dcdphi, d2cdphi2)
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
        COMPUTE_DPHI_A(c, dcdphi, d2cdphi2)
        COMPUTE_DPHI_A(tke, dtkedphi, d2tkedphi2) COMPUTE_DPHI_A(dis, ddisdphi, d2disdphi2)
        dpdphi   = (p_dyn.x[i][j][k+1] - p_dyn.x[i][j][k-1]) * inv_2dphi;
        dnuedphi = (nue.x[i][j][k+1]   - nue.x[i][j][k-1])   * inv_2dphi;
        #undef COMPUTE_DPHI_A

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


    // ===== Ocean physics: Coriolis, centrifugal, buoyancy =====
    // Coriolis acceleration -2 Ω × v in (r, θ, φ) with θ = colatitude
    // (the0 = 0 at N pole), v south-positive (= +ê_θ), w east-positive:
    //   F_r = +2Ω sin(θ) w               (Eötvös: east → up at equator)
    //   F_θ = +2Ω cos(θ) w
    //   F_φ = -2Ω cos(θ) v - 2Ω sin(θ) u
    // EXPERIMENT 2026-06-26 (uncommitted): inverse Rossby number = 2*Omega*L/u_0
    // with NO dt (RK4 supplies the timestep); the old coefficient's extra *dt made
    // Coriolis ~1/dt too weak vs advection/diffusion -> SH gyre spin-down. Test:
    // drop the *dt. Active (turbulent) path. See project_hydro_coriolis_dt_scaling.
    // ATOM_CORIOLIS_NONTRAD (AtomUtils::coriolis_nontraditional, lib/Utils.h) — the SAME switch
    // the atmosphere reads, so the two spheres can no longer disagree about which approximation
    // they are making. Default 0 = traditional: F_r drops and F_phi keeps only its cos(theta)
    // part. The two dropped terms are one energetically consistent pair and go together.
    const double nontrad = AtomUtils::coriolis_nontraditional() ? 1.0 : 0.0;
    double two_omega_Lhyd_over_u0 = 2.0 * omega * L_hyd / u_0;
    double Coriolis_rad =  nontrad * two_omega_Lhyd_over_u0 * sinthe * w_ijk;
    double Coriolis_the =  two_omega_Lhyd_over_u0 * costhe * w_ijk;
    double Coriolis_phi = -two_omega_Lhyd_over_u0 * (costhe * v_ijk + nontrad * sinthe * u_ijk);

    double rad_dist     = (double)i * L_hyd * exp_rm;
    double rad_Earth_m  = rad_dist + r_Earth * 1e3;
    double dt_over_u0_sq = (dt / u_0) * (dt / u_0);
    double omega2       = omega * omega;

    // ----- Centrifugal acceleration. Three separate errors, all in this one term. -----
    //
    // The centrifugal acceleration in a frame rotating about the polar axis points AWAY from that
    // axis, and the unit vector away from the axis is sin(theta)*e_r + cos(theta)*e_theta, so
    //     a_r     = +Omega^2 * r * sin^2(theta)
    //     a_theta = +Omega^2 * r * sin(theta) * cos(theta)
    // What stood here was Omega^2*r and Omega^2*r*|sin(theta)|, entering rhs_u/rhs_v with a MINUS:
    //   sign   the force pointed toward the rotation axis instead of away from it;
    //   sin^2  the radial part had no sin(theta) at all, so it was at full strength at the poles,
    //          where the true centrifugal force is exactly zero;
    //   |sin|  the meridional part cannot change sign at the equator, so it pushed the same way in
    //          both hemispheres instead of pointing toward the equator in each. cos(theta) is what
    //          carries that hemispheric sign, and it was missing.
    // Identical to the ATJUP defect fixed in 8649675; the two models share this lineage.
    //
    // sinthe is CLAMPED to >= 0.4 (RungeKutta_Hyd_Turb.cpp:54) for the 1/sin^2 metric divisions.
    // That floor is a metric stabiliser and has no business in a body force — with it the term
    // stays finite at the pole. The true sine is recovered from the (unclamped) cosine instead;
    // theta runs 0..pi, so sin(theta) >= 0 and the positive root is the right one.
    const double sinthe_true = sqrt(std::max(0.0, 1.0 - costhe * costhe));

    // NOT changed here: the magnitude scaling (dt/u_0)^2. Whether that is the right
    // nondimensionalisation is the separate, model-wide body-force unit-system question; this
    // fixes the direction and the angular structure only, so the two can be judged apart.
    //
    // A WARNING that no coefficient can fix, and the reason to watch what this does to the
    // meridional circulation. On the real planet nothing has to balance these: they are absorbed
    // into the geopotential, and the answer is the Earth's oblateness — the surfaces of constant
    // effective gravity ARE that shape. This model has a spherical grid, a spherical bottom
    // boundary and a constant radial g, so the meridional part has nothing to work against and
    // will drive a permanent pole-to-equator acceleration that is entirely an artefact. The
    // physical way to carry it is to fold it into an effective gravity and never write it as a
    // force at all. Until that is done, <centrifugal>0</centrifugal> is the more defensible
    // setting than the present 1 — but that is a modelling decision and is left alone here.
    double centrifugal_rad = omega2 * rad_Earth_m * sinthe_true * sinthe_true * dt_over_u0_sq;
    double centrifugal_the = omega2 * rad_Earth_m * sinthe_true * costhe      * dt_over_u0_sq;

    double coeff_energy_p = u_0 * u_0 / (cp_w * t_0);


    // ===== Turbulence: diffusion coefficients and source refresh =====
    double diffusion_t_re   = 0.0;
    double diffusion_vel_re = 0.0;
    double diffusion_tke_re = 0.0;
    double diffusion_dis_re = 0.0;
    const double diff_prec_re_inv = 1.0 / (sc * re_turb);

    // Wall-distance index for SST blending — use uniform grid index
    // (no stretched-coordinate table in the hydrosphere model).
    const double d_i       = (double)i;
    const int    i_bath    = i_bathymetry[j][k];
    const double d_i_floor = (double)i_bath;

    // Refresh prod.x at every RK4 sub-stage from current velocity derivatives.
    if (use_turbulence_model && !land_ijk) {
        constexpr double C_nue_loc   = 0.028;
        constexpr double dis_min_loc = 1.0e-10;
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
        const double nue_max  = 1000.0 / (u_0 * L_hyd);
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

    if (use_k_epsilon_turbulence_model) {
        const double sig_k   = 1.0;
        const double sig_w   = 1.3;
        const double C_eps_1 = 1.35;
        const double C_eps_2 = 1.80;
        {
            diffusion_t_re   = (1.0/re_turb + nue.x[i][j][k]) / pr_turb;
            diffusion_vel_re = 1.0/re_turb + nue.x[i][j][k];
            diffusion_tke_re = 1.0/re_turb + nue.x[i][j][k] / sig_k;
            diffusion_dis_re = 1.0/re_turb + nue.x[i][j][k] / sig_w;
        }
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

    if (use_k_omega_turbulence_model) {
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

    if (use_k_omega_SST_turbulence_model) {
        const double bet_star = 0.09;
        const double sig_k1   = 1.176;
        const double sig_k2   = 1.0;
        const double sig_w1   = 2.0;
        const double sig_w2   = 1.168;

        // Blending uses sigma values as denominators (divide convention).
        // sig_k1=1.176=1/0.85 etc. are reciprocals of Menter 1994 Prandtl numbers.
        auto blend = [](double inner, double outer, double F1){ return F1*inner + (1.0-F1)*outer; };

        // Molecular viscosity of seawater [m²/s], non-dimensionalised
        constexpr double nue_water_phys = 1.0e-6;
        const double nue_water_nd = nue_water_phys / (u_0 * L_hyd);

        int j_beg = 0, k_beg = 0;
        double d_j_mount = 0.0, d_k_mount = 0.0;
        double d_j_south = 0.0, d_j_north = 0.0;
        double d_k_east  = 0.0, d_k_west  = 0.0;
        double d_j = (double)j;
        double d_k = (double)k;

        if (land_ijk && is_water(h, i, j+1, k))   j_beg = j;
        for (int jj = j_beg; jj < jm-1; jj++)
            if (is_land(h, i, jj, k) || jj == jm-1) d_j_mount = jj;
        d_j_south = fabs(d_j - d_j_mount);

        if (land_ijk && is_water(h, i, j-1, k))   j_beg = j;
        for (int jj = jm-1; jj > j_beg; jj--)
            if (is_land(h, i, jj, k) || jj == 0) d_j_mount = jj;
        d_j_north = fabs(d_j - d_j_mount);

        if (land_ijk && is_water(h, i, j, k+1))   k_beg = k;
        for (int kk = k_beg; kk < km-1; kk++)
            if (is_land(h, i, j, kk) || kk == km-1) d_k_mount = kk;
        d_k_east = fabs(d_k - d_k_mount);

        if (land_ijk && is_water(h, i, j, k-1))   k_beg = k;
        for (int kk = km-1; kk > k_beg; kk--)
            if (is_land(h, i, j, kk) || kk == 0) d_k_mount = kk;
        d_k_west = fabs(d_k - d_k_mount);

        double d_close = d_i - d_i_floor;
        if (d_j_south <= d_close) d_close = d_j_south;
        if (d_j_north <= d_close) d_close = d_j_north;
        if (d_k_east  <= d_close) d_close = d_k_east;
        if (d_k_west  <= d_close) d_close = d_k_west;
        d_close = std::max(d_close, 1.0e-6);

        double CD_kw = std::max(
            (2.0 * sig_w2) / dis.x[i][j][k]
            * (dtkedr   * ddisdr
            + dtkedthe * inv_rm       * ddisdthe  * inv_rm
            + dtkedphi * inv_rmsinthe * ddisdphi  * inv_rmsinthe),
            1.0e-20);

        double arg1 = std::min(std::max(
            sqrt(tke.x[i][j][k]) / (bet_star * dis.x[i][j][k] * d_close),
            (500.0 * nue_water_nd) / (d_close * d_close * dis.x[i][j][k])),
            (4.0 * sig_w2 * tke.x[i][j][k]) / (CD_kw * d_close * d_close));
        double F1 = tanh(pow(arg1, 4));

        {
            diffusion_t_re   = (1.0/re_turb + nue.x[i][j][k]) / pr_turb;
            diffusion_vel_re =  1.0/re_turb + nue.x[i][j][k];
            diffusion_tke_re =  1.0/re_turb + nue.x[i][j][k] / blend(sig_k1, sig_k2, F1);
            diffusion_dis_re =  1.0/re_turb + nue.x[i][j][k] / blend(sig_w1, sig_w2, F1);
        }

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
            const double Y_w    = blend(bet1, bet2, F1) * dis_s * dis_s;
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

    if (!use_turbulence_model) {
        diffusion_t_re   = 1.0 / (re_turb * pr_turb);
        diffusion_vel_re = 1.0 / re_turb;
        diffusion_tke_re = 0.0;
        diffusion_dis_re = 0.0;
        tke.x[i][j][k]        = 0.0;
        dis.x[i][j][k]        = 0.0;
        nue.x[i][j][k]        = 0.0;
        prod.x[i][j][k]       = 0.0;
        tke_source.x[i][j][k] = 0.0;
        dis_source.x[i][j][k] = 0.0;
        tken.x[i][j][k]       = 0.0;
        disn.x[i][j][k]       = 0.0;
    }

    // Inviscid spin-up: ramp the diffusion (molecular background + eddy nue) from 0, mirroring
    // the former laminar RHS_Hyd.cpp. Lets the turbulent path also serve the inviscid phase
    // (diffusion_ramp=0 => Euler) and the laminar case (use_turbulence_model=false => molecular
    // only, above), so it is now the single hyd dynamical core. No-op post-spinup (ramp=1).
    diffusion_t_re   *= diffusion_ramp;
    diffusion_vel_re *= diffusion_ramp;
    diffusion_tke_re *= diffusion_ramp;
    diffusion_dis_re *= diffusion_ramp;


    // ===== Transport (advection) =====
    // Rhie-Chow: advect with the divergence-free FACE fluxes (uf/vf/wf, filled by
    // PressureSolverHyd::project_velocity in the heavy block), reconstructed to the
    // cell centre as the average of the two bracketing faces. This is what removes
    // the collocated checkerboard: after the projection the cell-centre velocity
    // still carries a small odd-even pressure mode, but the transporting velocity
    // does not, so central advection no longer re-reads and amplifies it. Fall back
    // to the pointwise cell velocity at the domain edges (no bracketing face). uf is
    // a bare radial velocity (units of u), so the metric factor exp_rm still applies
    // once below. See project_hydro_continuity_checkerboard.
    // Advect with the divergence-free face fluxes from the face-consistent
    // projection (PressureSolverHyd::project_velocity). See
    // project_hydro_continuity_checkerboard.
    constexpr bool RHIE_CHOW_ADVECTION = true;
    double u_tr = u_ijk, v_tr = v_ijk, w_tr = w_ijk;
    if (RHIE_CHOW_ADVECTION) {
        if (i >= 1) u_tr = 0.5 * (uf.x[i-1][j][k] + uf.x[i][j][k]);
        if (j >= 1) v_tr = 0.5 * (vf.x[i][j-1][k] + vf.x[i][j][k]);
        if (k >= 1) w_tr = 0.5 * (wf.x[i][j][k-1] + wf.x[i][j][k]);
    }
    double u_exp   = u_tr * exp_rm;
    double v_invrm = v_tr * inv_rm;
    double w_invrs = w_tr * inv_rmsinthe;

    double pressure_t = coeff_energy_p * (u_exp * dpdr + v_invrm * dpdthe + w_invrs * dpdphi);

    double transport_t   = u_exp * dtdr   + v_invrm * dtdthe   + w_invrs * dtdphi;
    double transport_u   = u_exp * dudr   + v_invrm * dudthe   + w_invrs * dudphi;
    double transport_v   = u_exp * dvdr   + v_invrm * dvdthe   + w_invrs * dvdphi;
    double transport_w   = u_exp * dwdr   + v_invrm * dwdthe   + w_invrs * dwdphi;

    // ATOM_METRIC_CURVATURE — same terms and the same warning as the atmosphere; see lib/Utils.h.
    // The hydrosphere has no metric-radius knob yet, so here this is a measurement tool, not a fix.
    if(AtomUtils::metric_curvature()){
        const double cotanthe = costhe / sinthe;
        transport_u += -(v_ijk * v_ijk + w_ijk * w_ijk) * inv_rm;
        transport_v +=  (u_ijk * v_ijk - w_ijk * w_ijk * cotanthe) * inv_rm;
        transport_w +=  (w_ijk * u_ijk + v_ijk * w_ijk * cotanthe) * inv_rm;
    }
    double transport_c   = u_exp * dcdr   + v_invrm * dcdthe   + w_invrs * dcdphi;
    double transport_tke = u_exp * dtkedr + v_invrm * dtkedthe + w_invrs * dtkedphi;
    double transport_dis = u_exp * ddisdr + v_invrm * ddisdthe + w_invrs * ddisdphi;

    // Advective form only (v·∇k, v·∇ω) — matches the standard k-ε / k-ω derivations
    // (Pope, Wilcox, Menter).  A previous version added the conservative-form
    // correction `+tke·∇·v` and `+dis·∇·v`; in an incompressible flow ∇·v=0 and the
    // two forms agree, but the sea-water solver here is not enforced divergence-free
    // until the pressure projection completes, so the term produced spurious
    // exponential growth of tke wherever the divergence persisted between projections.
    // The atmosphere lost ~1e98 m²/s² at one cell to the same bug; removing the
    // correction restores the standard derivation.  Y_w = β₀·dis² is quadratic and
    // self-limiting, so dis was less affected.


    // ===== Diffusion =====
    double two_over_rm_exp = 2.0 * inv_rm * exp_rm;
    double cos_rm2sin      = costhe_inv_rm2sinthe;
    // v_metric = csc²θ/r² = (1 + cot²θ)/r² = (1 + cos²θ/sin²θ)/r². 2026-07-21: this hyd
    // turbulent core still carried the pre-fix costhe (not costhe²) form of the
    // project_vmetric_antidiffusion bug — the same error corrected in RHS_Atm_Turb.cpp on
    // 2026-06-19 but never propagated here. With bare costhe the metric under-damps v,w in
    // one hemisphere and anti-diffuses poleward of ~38° in the other. A/B validated from the
    // Ma=100 iter-1000 restart (150 steps): stable, 0 NaN, KE differs ~0.01% (T bit-identical).
    double v_metric        = (1.0 + costhe * costhe / sinthe2) * inv_rm2;

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

    // ==================================================================================
    // HYD_A_H -- AN EXPLICIT HORIZONTAL EDDY VISCOSITY, in m^2/s. Default 0 = OFF and
    // bit-identical.
    //
    // WHY THIS EXISTS. This ocean has no horizontal eddy viscosity of its own. It has ONE
    // isotropic coefficient, `diffusion_vel_re` = 1/re_turb + nue, applied to the whole
    // spherical Laplacian, whose horizontal part carries `inv_rm2` -- and `rm` is the GRID
    // coordinate, running 1.000..2.000 with one unit = L_hyd = 200 m. So the horizontal
    // diffusion has been acting over a metric radius of 200-400 m instead of 6370 km, i.e.
    // it is (r_Earth/(rm*L_hyd))^2 ~ 4e8 too strong. Correcting the metric
    // (HYD_METRIC_RADIUS) removes that factor and leaves NOTHING to damp grid-scale
    // structure: measured 2026-09-04 at iteration 1000, the repaired arm carries 2.3x the
    // control's grid-scale noise (rms 5-point Laplacian of surface w over rms w, 0.72 ->
    // 1.69) and a 112 cm/s outlier where the control's max is 49.6. That is why
    // HYD_METRIC_RADIUS is not usable, and this term is the missing piece.
    //
    // WHAT THE MODEL HAS BEEN GETTING BY ACCIDENT, MEASURED RATHER THAN ARGUED. From
    // output_za's iteration-300 field the k-epsilon `nue` has a median of 9.54e-05 nd =
    // **0.0046 m^2/s** -- a molecular-scale number, as an ocean turbulence closure tuned on
    // the VERTICAL should give. Put that through the broken horizontal metric and the
    // IMPLIED horizontal viscosity is **1.2e6 m^2/s at the surface and 4.6e6 at the
    // seafloor**, against the 1e2..1e4 m^2/s a 1-degree ocean model actually uses. The
    // metric error is not a small correction to the horizontal mixing; it IS the horizontal
    // mixing, at ~100x the largest defensible value.
    //
    // AND THE PHYSICAL VALUE CANNOT ACT HERE, WHICH IS THE FINDING AND IT NEEDED NO RUN.
    // A Laplacian viscosity damps a 2*dx mode with e-folding 1/(A_h*k^2), k = pi/dx and
    // dx = 111.2 km at the equator, so k^2 = 7.99e-10:
    //
    //     A_h = 1e2 m^2/s  ->  145 days   = 1.5e8 iterations
    //     A_h = 1e4 m^2/s  ->  1.45 days  = 1.5e6 iterations
    //     A_h = 1e6 m^2/s  ->  21 min     = 1.5e4 iterations
    //     A_h = 1.5e7      ->  83 s       = 1.0e3 iterations   <- one run of this model
    //
    // One iteration is 0.0833 s and the longest run in this tree is 1000. **So a
    // physically-sized horizontal eddy viscosity is three orders of magnitude too slow to
    // damp anything on any run this tree can afford**, and the accidental 1e6 the metric
    // error supplies is the smallest value that comes close. This is the same wall as
    // ATM_HYDRO_PGF's 1/f and HYD_BAROCLINIC_PGF's: the model's iteration budget, not the
    // physics. Read A_h at the values that ACT as a NUMERICAL knob in the sense
    // sst_relax_alpha and the atmosphere's omega_teq are numerical, and do not call a value
    // chosen to control noise an eddy viscosity.
    //
    // STABILITY IS NEVER THE BINDING CONSTRAINT, which is worth recording because it is the
    // usual reason such a term is kept small. Explicit stability wants
    // A_h/(u_0*L_hyd)*(L_hyd/r_E)^2 <= 0.5*dthe^2/dt = 1.52, i.e. A_h <= 7e10 m^2/s. Every
    // value above is 3 to 9 orders under it.
    //
    // THE OPERATOR is the HORIZONTAL part of the existing `diffusion_v`/`diffusion_w`
    // expressions, term for term including the curvature and the u/v/w cross-couplings, so
    // this is the same vector Laplacian with a different coefficient rather than a second
    // and differently-shaped operator. Horizontal momentum only: `u` is the RADIAL
    // component, and the point of the metric repair is that its spurious value collapses --
    // giving it horizontal diffusion would re-couple exactly what was separated.
    //
    // THE METRIC HERE IS THE EARTH'S, DELIBERATELY, and does NOT follow `inv_rm`. Written
    // with `inv_rm` the coefficient would mean 4e8 different things depending on whether
    // HYD_METRIC_RADIUS is set, and A_h would stop being a viscosity in m^2/s. Same
    // reasoning, and the same exception to "match the neighbouring term's structure", as
    // HYD_BAROCLINIC_PGF above.
    //
    // Carries `diffusion_ramp` like every other diffusion coefficient, so an inviscid
    // spin-up stays inviscid.
    //
    // NEXT, AND NOT WRITTEN: a BIHARMONIC form. At a strength that damps 2*dx in one run
    // length, a Laplacian also damps the domain-scale flow by (k_2dx/k_L)^2 ~ 8.3e3, while
    // a biharmonic damps it by (k_2dx/k_L)^4 ~ 6.9e7 -- four orders more scale-selective,
    // which is exactly why ocean models at this resolution use one. That is the instrument
    // that could remove the noise without removing the circulation with it.
    // ==================================================================================
    static const double a_h_phys = [](){
        const char* e = getenv("HYD_A_H"); return e ? atof(e) : 0.0; }();

    if(a_h_phys != 0.0){
        const double inv_r_true = L_hyd / (r_Earth * 1.0e3);          // L_hyd / R_Earth [m/m]
        const double a_h_nd     = a_h_phys / (u_0 * L_hyd)
                                * inv_r_true * inv_r_true * diffusion_ramp;

        const double cot_the = costhe / sinthe;
        const double inv_s2  = 1.0 / sinthe2;

        diffusion_v += a_h_nd * (d2vdthe2 + dvdthe * cot_the
                               - (1.0 + cot_the * cot_the) * v_ijk
                               + d2vdphi2 * inv_s2
                               + 2.0 * dudthe
                               - 2.0 * dwdphi * costhe * inv_s2);

        diffusion_w += a_h_nd * (d2wdthe2 + dwdthe * cot_the
                               - (1.0 + cot_the * cot_the) * w_ijk
                               + d2wdphi2 * inv_s2
                               + 2.0 * dudphi / sinthe
                               + 2.0 * dvdphi * costhe * inv_s2);
    }

    // ==================================================================================
    // HYD_A_H_BIHARM=<B in m^4/s> -- OUTER PASS of the horizontal biharmonic viscosity.
    // Default 0.0 = OFF and bit-identical.
    //
    // WHY A BIHARMONIC AT ALL, AND THE NUMBER THAT DECIDES IT. HYD_A_H above can damp
    // grid-scale noise only at ~1.5e7 m^2/s, three orders above any physical value, and a
    // Laplacian that strong reaches every scale: at the strength that gives a 2*dx
    // e-folding of 83 s (one 1000-iteration run) its DOMAIN-scale e-folding is 7.8 days.
    // The biharmonic at the equivalent strength, B = 1.89e16 m^4/s, gives a domain-scale
    // e-folding of **62 878 days** -- the same grid-scale control, and **8 090x** less
    // damping of the circulation, because selectivity goes as (k_2dx/k_L)^4 = 6.5e7
    // against the Laplacian's (k_2dx/k_L)^2 = 8.1e3. That ratio is the whole argument for
    // the term, and it is why ocean models at this resolution use one.
    //
    //     B = 1e11 m^4/s (a physical 1-degree value) -> 2*dx e-folding 181 days
    //     B = 1e12                                   -> 18.2 days
    //     B = 1.89e16                                -> 83 s = 1000 iterations
    //
    // So the physical value is as unreachable here as the Laplacian's, for the same reason
    // -- one iteration is 0.0833 s -- and B is a NUMERICAL knob at the values that act.
    // What changes is the COST of using a numerical value, and that is what the biharmonic
    // buys: at equal grid-scale control it leaves the circulation alone by 8 000x.
    //
    // STABILITY, again not the binding constraint: explicit nabla^4 wants
    // dt*B_nd/dthe^4 <= 1/8, i.e. B <= **2.3e20 m^4/s**, four orders above the strongest
    // value in the table.
    //
    // NON-DIMENSIONALISATION. Physical dv/dt = -B*nabla^4_h v, and nabla^4_h = H^2/r^4 with
    // H the angular operator the inner pass stored, so the coefficient is
    // B*L_hyd/(u_0*r_Earth^4) -- the same shape as HYD_A_H's A_h*L_hyd/(u_0*r_Earth^2), one
    // more power of the radius. As there, the radius is the EARTH'S and not inv_rm, so B
    // means m^4/s whether or not HYD_METRIC_RADIUS is set.
    //
    // The MINUS sign is the physics: nabla^4 with a plus sign is anti-diffusive.
    // ==================================================================================
    if(HydHorizVisc::biharm_strength() != 0.0 && !HydHorizVisc::ok.empty()
       && j > 0 && j < jm-1 && k > 0 && k < km-1){

        const std::size_t p  = HydHorizVisc::idx(i, j,   k);
        const std::size_t pT = HydHorizVisc::idx(i, j+1, k);
        const std::size_t pB = HydHorizVisc::idx(i, j-1, k);
        const std::size_t pP = HydHorizVisc::idx(i, j,   k+1);
        const std::size_t pM = HydHorizVisc::idx(i, j,   k-1);

        // every stencil point of the OUTER application must have had its INNER application
        // computed, or the difference is taken against a zero that means "not evaluated".
        if(HydHorizVisc::ok[p] && HydHorizVisc::ok[pT] && HydHorizVisc::ok[pB]
                               && HydHorizVisc::ok[pP] && HydHorizVisc::ok[pM]){

            // B * L_hyd / (u_0 * r_Earth^4), dimensionless; 5.06e-25 per (m^4/s) here.
            const double r_m    = r_Earth * 1.0e3;
            const double r_m2   = r_m * r_m;
            const double b_h_nd = HydHorizVisc::biharm_strength() * L_hyd
                                / (u_0 * r_m2 * r_m2) * diffusion_ramp;

            const double cot_the = costhe / sinthe;
            const double inv_s2  = 1.0 / sinthe2;

            const double LvC = HydHorizVisc::lap_v[p],  LwC = HydHorizVisc::lap_w[p];
            const double LvT = HydHorizVisc::lap_v[pT], LvB = HydHorizVisc::lap_v[pB];
            const double LwT = HydHorizVisc::lap_w[pT], LwB = HydHorizVisc::lap_w[pB];
            const double LvP = HydHorizVisc::lap_v[pP], LvM = HydHorizVisc::lap_v[pM];
            const double LwP = HydHorizVisc::lap_w[pP], LwM = HydHorizVisc::lap_w[pM];

            const double H2v = (LvT - 2.0*LvC + LvB) * inv_dthe2
                             + cot_the * (LvT - LvB) * inv_2dthe
                             - (1.0 + cot_the*cot_the) * LvC
                             + (LvP - 2.0*LvC + LvM) * inv_dphi2 * inv_s2
                             - 2.0 * costhe * inv_s2 * (LwP - LwM) * inv_2dphi;

            const double H2w = (LwT - 2.0*LwC + LwB) * inv_dthe2
                             + cot_the * (LwT - LwB) * inv_2dthe
                             - (1.0 + cot_the*cot_the) * LwC
                             + (LwP - 2.0*LwC + LwM) * inv_dphi2 * inv_s2
                             + 2.0 * costhe * inv_s2 * (LvP - LvM) * inv_2dphi;

            diffusion_v -= b_h_nd * H2v;
            diffusion_w -= b_h_nd * H2w;
        }
    }

    double diffusion_c = (d2cdr2 * exp_2_rm + dcdr * two_over_rm_exp
        + d2cdthe2 * inv_rm2 + dcdthe * cos_rm2sin
        + d2cdphi2 * inv_rm2sinthe2) * diff_prec_re_inv;

    double diffusion_tke = (d2tkedr2 * exp_2_rm + dtkedr * two_over_rm_exp
        + d2tkedthe2 * inv_rm2 + dtkedthe * cos_rm2sin
        + d2tkedphi2 * inv_rm2sinthe2) * diffusion_tke_re;

    double diffusion_dis = (d2disdr2 * exp_2_rm + ddisdr * two_over_rm_exp
        + d2disdthe2 * inv_rm2 + ddisdthe * cos_rm2sin
        + d2disdphi2 * inv_rm2sinthe2) * diffusion_dis_re;


    // ===== RHS assembly =====
    double dpdr_exp     = dpdr   * exp_rm;
    double dpdthe_invrm = dpdthe * inv_rm;
    double dpdphi_invrs = dpdphi * inv_rmsinthe;

    // ==================================================================================
    // HYD_BAROCLINIC_PGF -- THE HORIZONTAL HYDROSTATIC PRESSURE GRADIENT IN rhs_v/rhs_w.
    // A STRENGTH, not a flag. Default 0.0 = OFF and bit-identical.
    //
    // WHY. The comment at the buoyancy term below has been asking for this term by name:
    // "a real baroclinic ocean needs the horizontal grad(p_hydro) in rhs_v/rhs_w, not just
    // this term." The two `dpd*` lines above read `p_dyn`, a PROJECTION pressure whose Poisson
    // source is the divergence of the tendency -- it removes the divergent part of forces
    // already present and cannot manufacture the pressure a mass field implies. `p_hydro`,
    // which ThermoHyd rebuilds every iteration and which carries the whole density structure,
    // appears in NO momentum equation. The only other route from t or c into momentum is
    // `buoy_nd` below, which is RADIAL and which its own comment records as carrying a
    // spurious dt, "~1e-7 too small, so buoyancy is effectively OFF and the ocean is
    // dynamically barotropic (wind + Coriolis only)". This is the atmosphere's missing
    // thermal wind, one model over; see CLAUDE.md.
    //
    // HORIZONTAL ONLY, for the same reason the atmosphere's ATM_HYDRO_PGF is. The radial
    // hydrostatic gradient IS gravity to within a per cent, and this model handles the radial
    // direction the Boussinesq way -- `buoy_nd` is the residual and `rhs_u` carries no -g.
    // Adding dp_hydro/dr without the matching -g would be a ~1e5 spurious sinking. Do not
    // "complete" this with a third line.
    //
    // ---- THE METRIC, AND WHY THIS TERM DELIBERATELY DOES NOT MATCH ITS NEIGHBOURS ----
    //
    // MEASURED 2026-09-04, AND IT IS A DEFECT IN THE OCEAN THAT THE ATMOSPHERE DOES NOT HAVE.
    // The horizontal derivatives above are scaled by `inv_rm` with rm = rad.z[i], and the
    // ocean's rad.z runs 1.000 .. 2.000 (StretchedCoordinates from r0 = 1.0 over a length
    // (im-1)*dr = 1.0). The implied metric radius is therefore rm*L_hyd = 200 .. 400 m. In the
    // atmosphere the same expression carries rmet ~ 397.5, and 397.5 * L_atm = 397.5 * 16024 m
    // = 6.370e6 m -- the Earth's radius, exactly. THE OCEAN IS MISSING THAT NORMALISATION.
    // Consequence, measured from the iteration-300 field of `output_ps0` rather than argued:
    // the code's median |horizontal advection| / |Coriolis| is 18 972x the same ratio formed
    // physically, and R_Earth/(rm*L_hyd) over rm in [1,2] is 31 850 .. 15 925 -- the measured
    // number sits inside that band. So the ocean's horizontal advection, diffusion and p_dyn
    // gradient are all ~2e4 too large relative to Coriolis, which is the one horizontal term
    // that needs no metric (it is 2*omega*L_hyd/u_0 * cos(theta) * w, a proper inverse Rossby
    // number) and is therefore the one that is right.
    //
    // THIS TERM USES THE EARTH'S RADIUS, NOT `inv_rm`, AND THAT IS A JUDGEMENT WORTH STATING.
    // Matching the neighbours would be the usual rule and it is wrong here for a specific
    // reason: the entire purpose of this force is to balance CORIOLIS -- that balance is what
    // geostrophy and thermal wind ARE -- and Coriolis is the metric-free term. Written with
    // `inv_rm` this force would be ~2e4 times the Coriolis force it is meant to oppose, and the
    // balance it exists to create would be impossible by construction. Written with R_Earth it
    // is correctly sized against Coriolis and wrongly sized against advection and diffusion,
    // which are themselves wrong. That is the better of two bad options and it is not a repair
    // of the metric: fixing rad.z would move advection, diffusion, the Poisson operator and the
    // projection together, and is a far larger change than this one.
    //
    // THE EULER NUMBER. The two `dpd*` lines above carry none, which fixes `p_dyn`'s
    // normalisation at rho*u_0^2 (the same argument that corrected the atmosphere's p_dyn label
    // from p_0). `p_hydro` is a REAL pressure in BAR, and a physical acceleration enters this
    // RHS multiplied by L_hyd/u_0^2 -- read off the Coriolis term, which is
    // (2*omega*L_hyd/u_0)*cos(theta)*w_nd = [2*omega*cos(theta)*w_phys] * L_hyd/u_0^2. Hence
    //     dv*/dt* = -(1e5 * L_hyd / (rho_0_water * u_0^2 * R_Earth)) * d(p_hydro)/dthe
    // = 0.0547 per bar per radian at the shipped constants.
    //
    // NO NORMALISATION VARIANT IS NEEDED, WHICH THE ATMOSPHERE COULD NOT SAY. There,
    // differencing p_stat raw put a spurious barotropic gradient at the ground because
    // densities() sets p_sl proportional to the surface temperature, and the default had to
    // renormalise every column at a common height. Here the surface value is
    // p_hydro[im-1] = p_0/1000, the SAME CONSTANT in every column, so the raw field already has
    // zero horizontal gradient at the surface and carries only the baroclinic part beneath it.
    // And the ocean grid is a z-coordinate -- level i is the same geometric depth in every
    // column -- so differencing at constant i is differencing at constant depth and the
    // sigma-coordinate pressure-gradient error cannot arise at all.
    //
    // LAND IS THE ONE HAZARD, AND IT IS REAL: ThermoHyd sets p_hydro = 0 in rock, so a centred
    // difference across a coast or over the seafloor would see 0 against ~12 bar and
    // manufacture an enormous gradient. Fall back to a one-sided difference at a land face,
    // exactly as ThermoHyd::Forces() already does for its PresGradForce diagnostic.
    //
    // IT SURVIVES THE PROJECTION, WHICH THE BAROTROPIC MODE DID NOT. The note further down
    // records that f x u_bt was not injected as a force because it is curl-free and the
    // divergence projection deletes it. A BAROCLINIC gradient is not curl-free -- its curl is
    // the thermal wind, which is the whole point -- so the solenoidal part survives. It is
    // therefore left inside aux_v/aux_w (only the p_dyn gradient is added back below), so the
    // projection removes its divergent part and keeps the rest.
    //
    // WHAT TO EXPECT ON THE FIRST RUN. The v/w budgets change immediately; the VELOCITY will
    // not, on any run this tree can afford. Geostrophic adjustment takes 1/f, which at 30 deg
    // is ~2.7e4 s = 2.7e8 iterations at dt = 1e-4 with L_hyd/u_0 = 833 s. That is the same
    // ceiling ATM_HYDRO_PGF met and the reason its default is also 0.0.
    // ==================================================================================
    static const double bcl_pgf = [](){
        const char* e = getenv("HYD_BAROCLINIC_PGF"); return e ? atof(e) : 0.0; }();

    double bcl_the = 0.0, bcl_phi = 0.0;
    if(bcl_pgf != 0.0 && is_water(h, i, j, k)
       && j > 0 && j < jm-1 && k > 0 && k < km-1){
        const double pC = p_hydro.x[i][j][k];
        // One-sided at a land face: substitute the centre value, which zeroes that half of the
        // difference rather than differencing against rock.
        const double pTp = is_land(h,i,j+1,k) ? pC : p_hydro.x[i][j+1][k];
        const double pTm = is_land(h,i,j-1,k) ? pC : p_hydro.x[i][j-1][k];
        const double pPp = is_land(h,i,j,k+1) ? pC : p_hydro.x[i][j][k+1];
        const double pPm = is_land(h,i,j,k-1) ? pC : p_hydro.x[i][j][k-1];

        // 1e5 Pa/bar, / (rho u_0^2), x L_hyd/u_0^2 folded into one constant, and the metric
        // radius is the EARTH'S, for the reason argued at length above.
        const double eu = bcl_pgf * 1.0e5 * L_hyd
                        / (r_0_water * u_0 * u_0 * r_Earth * 1.0e3);

        // sinthe here is the CLAMPED metric sine (>= 0.4, RungeKutta_Hyd_Turb.cpp) -- this is a
        // metric division, not a body force's angular structure, so it takes the same floor
        // every other 1/sin in this file takes.
        bcl_the = eu * (pTp - pTm) * inv_2dthe;
        bcl_phi = eu * (pPp - pPm) * inv_2dphi / sinthe;

        if(!AtomUtils::is_finite_safe(bcl_the)) bcl_the = 0.0;
        if(!AtomUtils::is_finite_safe(bcl_phi)) bcl_phi = 0.0;
    }

    // ----- Shear (velocity-gradient) dissipation heating -----
    // Dissipation function Φ = ν|S|² built DIRECTLY from the RESOLVED strain-rate tensor
    // (the velocity gradients), not the parameterized ε. This is the KE→internal-energy
    // return that the k-ε closure computed but never gave back to T, so the energy budget
    // was not closed. Converted to a temperature tendency by the same KE→T coefficient the
    // pressure-work term uses (coeff_energy_p = u_0²/(cp_w·t_0)). Strain invariant, same
    // spherical components as the turbulence `prod` (the gradients are already to hand, so
    // this costs ~a dozen multiply-adds):
    //   |S|² = 2(S_rr² + S_θθ² + S_φφ²) + S_rθ² + S_rφ² + S_θφ²      (Φ/ν, ≥0)
    // ⚠ resolved (∇u)² is noisiest exactly at coasts (the known checkerboard/seam zones),
    // hence two guardrails: floor nue at 0, and hard-cap the tendency (dissip_cap).
    //
    // MAGNITUDE (measured 2026-07-17, Ma=200 A/B from the iter-400 restart, single-threaded
    // so the noise floor is exactly zero): at a 1e5× gain the peak ΔT is 7.9e-5 K after 20
    // iters, i.e. **~1e-9 K at the physical strength used here** — correct, and far below
    // anything observable. It is kept because it closes the budget, not because it does
    // work: the term never approaches dissip_cap (peak is ~10× under it). The SST lever is
    // ADVECTION, whose tendency is O(1) nondim against this term's ~1e-7. A Python-settable
    // C_dissip gain used to sit here to force an artificial coastal warming for a Picard
    // sweep; it was removed once measurement showed amplification buys nothing (and that the
    // surface is re-pinned every iter anyway, so surface ΔT is identically 0 by construction
    // — see cHydrosphereModel.cpp and project_dissipation_heating).
    constexpr double dissip_cap = 1.0e-3;             // guardrail cap on the nondim tendency

    // resolved strain-rate components (metric-scaled), matching the `prod` block
    const double Srr = dudr   * exp_rm;
    const double Stt = dvdthe * inv_rm;
    const double Spp = dwdphi * inv_rmsinthe;
    const double Srt = dudthe * inv_rm       + dvdr   * exp_rm;
    const double Srp = dudphi * inv_rmsinthe + dwdr   * exp_rm;
    const double Stp = dvdphi * inv_rmsinthe + dwdthe * inv_rm;
    const double strain2 = 2.0 * (Srr*Srr + Stt*Stt + Spp*Spp)
                         + Srt*Srt + Srp*Srp + Stp*Stp;
    const double nue_here = std::max(nue.x[i][j][k], 0.0);
    const double dissip_heat_t = std::min(coeff_energy_p * nue_here * strain2, dissip_cap);

    rhs_t.x[i][j][k] = pressure_t - transport_t + diffusion_t + dissip_heat_t;

    // Boussinesq buoyancy with thermal + haline anomaly.  Only the anomaly
    // relative to the reference state (t = 1, c = 1) drives radial motion;
    // the hydrostatic reference is absorbed into the pressure split.
    // Warm (t > 1) rises; salty (c > 1) sinks.  alpha_S ≈ (β_S·S_0)/(α_T·T_0)
    // sets the haline / thermal weighting; see RHS_Hyd.cpp for derivation.
    constexpr double alpha_S = 0.5;
    // NB the g*dt/u_0 scaling carries the spurious dt of the pre-5f71571 convention
    // (~1e-7 too small), so buoyancy is effectively OFF and the ocean is dynamically
    // barotropic (wind + Coriolis only). Repairing it to g*alpha_T*t_0*L/u_0^2 was
    // TESTED (2026-07-14) and DESTABILISES: max radial velocity 0.04 -> 19.5 m/s in
    // 30 iters, a polar vertical-velocity runaway (project_hydro_polar_blowup) — the
    // vertical buoyancy has nothing to balance it (no baroclinic/hydrostatic PGF in
    // the momentum eq). Left at the stable (barotropic) scaling; a real baroclinic
    // ocean needs the horizontal grad(p_hydro) in rhs_v/rhs_w, not just this term.
    const double buoy_nd = buoyancy * g * dt / u_0
        * ((t.x[i][j][k] - 1.0) - alpha_S * (c.x[i][j][k] - 1.0));
    rhs_u.x[i][j][k] = -dpdr_exp - transport_u + diffusion_u
        + buoy_nd
        + Coriolis * Coriolis_rad + centrifugal * centrifugal_rad;

    // Radial (u) momentum-budget capture — the DEEP polar blow-up is a
    // radial-velocity runaway, so split rhs_u into its five contributions on
    // checkpoint iters. write_deep_momentum_budget self-locates the max-|u|
    // interior ocean cell and dumps this split, so the driver is isolated
    // wherever the blow-up relocates (pole -> 84N -> 51S). The sum of the five
    // equals rhs_u. See project_hydro_polar_blowup.
    if (wbudget_capture) {
        ubud_pgf.x[i][j][k]  = -dpdr_exp;
        ubud_adv.x[i][j][k]  = -transport_u;
        ubud_diff.x[i][j][k] =  diffusion_u;
        ubud_buoy.x[i][j][k] =  buoy_nd;
        ubud_cor.x[i][j][k]  =  Coriolis * Coriolis_rad
                              + centrifugal * centrifugal_rad;
    }

    rhs_v.x[i][j][k] = -dpdthe_invrm - bcl_the - transport_v + diffusion_v
        + Coriolis * Coriolis_the + centrifugal * centrifugal_the;

    // v-momentum-budget capture (meridional). Attributes why the subtropical
    // gyres form their boundary current on the EASTERN coast instead of a
    // western-boundary current. vbud_wind is added in the wind-stress block
    // below (i=im-2 only). Sum of the five vbud_* equals rhs_v.
    if (wbudget_capture) {
        // With HYD_BAROCLINIC_PGF on, this bucket is the TOTAL pressure gradient --
        // projection plus hydrostatic -- so dyn_sum stays equal to the tendency and the
        // budget remains an identity. Separate the two with a HYD_BAROCLINIC_PGF=0 arm.
        vbud_pgf.x[i][j][k]  = -dpdthe_invrm - bcl_the;
        vbud_cor.x[i][j][k]  =  Coriolis * Coriolis_the
                              + centrifugal * centrifugal_the;
        vbud_adv.x[i][j][k]  = -transport_v;
        vbud_diff.x[i][j][k] =  diffusion_v;
        vbud_wind.x[i][j][k] =  0.0;   // set below at i=im-2
    }

    rhs_w.x[i][j][k] = -dpdphi_invrs - bcl_phi - transport_w + diffusion_w
        + Coriolis * Coriolis_phi;

    // w-momentum-budget capture (turbulent path = active path; mirror of
    // RHS_Hyd.cpp). Stores the four rhs_w contributions on checkpoint iters so
    // write_w_momentum_budget can attribute the zonal-mean zonal-velocity
    // tendency by latitude/depth. diffusion_w here carries the eddy viscosity.
    if (wbudget_capture) {
        wbud_pgf.x[i][j][k]  = -dpdphi_invrs - bcl_phi;   // total: projection + hydrostatic
        wbud_cor.x[i][j][k]  =  Coriolis * Coriolis_phi;
        wbud_adv.x[i][j][k]  = -transport_w;
        wbud_diff.x[i][j][k] =  diffusion_w;
    }

    // Surface wind-stress forcing (active turbulent path; mirror of RHS_Hyd.cpp).
    // The ocean momentum equation otherwise has NO sustained wind input — the
    // atmosphere wind only initialises the currents, so wind-driven currents spin
    // down (see project_hydro_no_wind_stress_forcing). Bulk stress
    // tau = r_air*C_D*|U_wind|*U_wind as a body force in the topmost integrated
    // cell (i=im-2 ~ mixed layer), non-dim by L_hyd/u_0^2 (advective time).
    constexpr bool WIND_STRESS_FORCING = true;
    if (WIND_STRESS_FORCING && i == im - 2 && is_water(h, i, j, k)) {
        constexpr double C_D = 2.6e-3;
        const double dz_top = (rad.z[im-1] - rad.z[im-2]) * L_hyd;        // [m]
        const double U_wind = sqrt(v_wind.y[j][k] * v_wind.y[j][k]
                                 + w_wind.y[j][k] * w_wind.y[j][k]);      // [m/s]
        const double accel_nd = r_air * C_D * U_wind / (r_0_water * dz_top)
                              * L_hyd / (u_0 * u_0);
        rhs_w.x[i][j][k] += accel_nd * w_wind.y[j][k];                    // zonal     (East+)
        rhs_v.x[i][j][k] += accel_nd * v_wind.y[j][k];                    // meridional (South+)
        if (wbudget_capture) vbud_wind.x[i][j][k] = accel_nd * v_wind.y[j][k];
    }

    // NB: the barotropic mode is NOT injected here as a force. The geostrophic PGF
    // f x u_bt is curl-free, so the divergence projection deletes it (verified
    // vestigial: 1.5% flow change, ACC stayed deep-westward). Instead the barotropic
    // (depth-mean) VELOCITY u_bt is imposed after project_velocity by
    // PressureSolverHyd::apply_barotropic_mode_split() — immune to projection.

    rhs_c.x[i][j][k] = -transport_c + diffusion_c;

    rhs_tke.x[i][j][k] = -transport_tke + diffusion_tke + tke_source.x[i][j][k];
    rhs_dis.x[i][j][k] = -transport_dis + diffusion_dis + dis_source.x[i][j][k];

    aux_u.x[i][j][k] = rhs_u.x[i][j][k] + dpdr_exp;
    aux_v.x[i][j][k] = rhs_v.x[i][j][k] + dpdthe_invrm;
    aux_w.x[i][j][k] = rhs_w.x[i][j][k] + dpdphi_invrs;
}
