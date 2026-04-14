#pragma once

#include "cAtmosphereModel.h"
#include "Utils.h"

#include <iostream>
#include <algorithm>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace AtomUtils;

class BC_Atm {
public:
    explicit BC_Atm(cAtmosphereModel& model)
        : m(model)
    {}

    // ------------------------------------------------------------------
    void bcSolidGround()
    {
        using namespace std;
        cout << endl << endl << endl << "      AGCM: BC_SolidGround" << endl;

        // 2nd-order one-sided (von Neumann) extrapolation into a land ghost cell from
        // two adjacent air cells: x_ghost = (4/3)*x_near - (1/3)*x_far.
        // Applied to the six fields that need smooth values at land/air interfaces:
        // dynamic pressure, water vapour, cloud water, cloud ice, graupel, CO2.
        //
        // Non-negative bounded scalars (cloud, ice, gr, c, convective moisture/mass-flux)
        // use a monotone clamp: result is restricted to [0, max(near, far)] to prevent
        // overshoot at steep windward gradients from feeding spurious cloud into the air.
        auto clamp_hydro = [](double extrap, double v1, double v2) -> double {
            return std::clamp(extrap, 0.0, std::max(v1, v2));
        };

        auto extrapolation = [&](int i, int j, int k,
                           int i1, int j1, int k1,
                           int i2, int j2, int k2) {
            m.tke.x[i][j][k] = std::max(0.0, m.c43 * m.tke.x[i1][j1][k1] - m.c13 * m.tke.x[i2][j2][k2]);
            m.dis.x[i][j][k] = std::max(0.0, m.c43 * m.dis.x[i1][j1][k1] - m.c13 * m.dis.x[i2][j2][k2]);
            m.nue.x[i][j][k] = std::max(0.0, m.c43 * m.nue.x[i1][j1][k1] - m.c13 * m.nue.x[i2][j2][k2]);

            m.t.x[i][j][k]       = m.c43 * m.t.x[i1][j1][k1]       - m.c13 * m.t.x[i2][j2][k2];

            m.p_dyn.x[i][j][k]   = m.c43 * m.p_dyn.x[i1][j1][k1]   - m.c13 * m.p_dyn.x[i2][j2][k2];

            m.r_dry.x[i][j][k]   = m.c43 * m.r_dry.x[i1][j1][k1]   - m.c13 * m.r_dry.x[i2][j2][k2];
            m.r_humid.x[i][j][k] = m.c43 * m.r_humid.x[i1][j1][k1] - m.c13 * m.r_humid.x[i2][j2][k2];

            m.c.x[i][j][k]     = clamp_hydro(m.c43 * m.c.x[i1][j1][k1]     - m.c13 * m.c.x[i2][j2][k2],
                                              m.c.x[i1][j1][k1],     m.c.x[i2][j2][k2]);
            m.cloud.x[i][j][k] = clamp_hydro(m.c43 * m.cloud.x[i1][j1][k1] - m.c13 * m.cloud.x[i2][j2][k2],
                                              m.cloud.x[i1][j1][k1], m.cloud.x[i2][j2][k2]);
            m.ice.x[i][j][k]   = clamp_hydro(m.c43 * m.ice.x[i1][j1][k1]   - m.c13 * m.ice.x[i2][j2][k2],
                                              m.ice.x[i1][j1][k1],   m.ice.x[i2][j2][k2]);
            m.gr.x[i][j][k]    = clamp_hydro(m.c43 * m.gr.x[i1][j1][k1]    - m.c13 * m.gr.x[i2][j2][k2],
                                              m.gr.x[i1][j1][k1],    m.gr.x[i2][j2][k2]);

            m.co2.x[i][j][k]   = m.c43 * m.co2.x[i1][j1][k1]   - m.c13 * m.co2.x[i2][j2][k2];

            m.S_c_c.x[i][j][k] = m.c43 * m.S_c_c.x[i1][j1][k1] - m.c13 * m.S_c_c.x[i2][j2][k2];
            m.S_v.x[i][j][k]   = m.c43 * m.S_v.x[i1][j1][k1]   - m.c13 * m.S_v.x[i2][j2][k2];
            m.S_c.x[i][j][k]   = m.c43 * m.S_c.x[i1][j1][k1]   - m.c13 * m.S_c.x[i2][j2][k2];
            m.S_i.x[i][j][k]   = m.c43 * m.S_i.x[i1][j1][k1]   - m.c13 * m.S_i.x[i2][j2][k2];
            m.S_r.x[i][j][k]   = m.c43 * m.S_r.x[i1][j1][k1]   - m.c13 * m.S_r.x[i2][j2][k2];
            m.S_s.x[i][j][k]   = m.c43 * m.S_s.x[i1][j1][k1]   - m.c13 * m.S_s.x[i2][j2][k2];

            m.q_v_u.x[i][j][k] = clamp_hydro(m.c43 * m.q_v_u.x[i1][j1][k1] - m.c13 * m.q_v_u.x[i2][j2][k2],
                                              m.q_v_u.x[i1][j1][k1], m.q_v_u.x[i2][j2][k2]);
            m.q_c_u.x[i][j][k] = clamp_hydro(m.c43 * m.q_c_u.x[i1][j1][k1] - m.c13 * m.q_c_u.x[i2][j2][k2],
                                              m.q_c_u.x[i1][j1][k1], m.q_c_u.x[i2][j2][k2]);
            m.q_v_d.x[i][j][k] = clamp_hydro(m.c43 * m.q_v_d.x[i1][j1][k1] - m.c13 * m.q_v_d.x[i2][j2][k2],
                                              m.q_v_d.x[i1][j1][k1], m.q_v_d.x[i2][j2][k2]);

            m.E_u.x[i][j][k] = clamp_hydro(m.c43 * m.E_u.x[i1][j1][k1] - m.c13 * m.E_u.x[i2][j2][k2],
                                            m.E_u.x[i1][j1][k1], m.E_u.x[i2][j2][k2]);
            m.E_d.x[i][j][k] = clamp_hydro(m.c43 * m.E_d.x[i1][j1][k1] - m.c13 * m.E_d.x[i2][j2][k2],
                                            m.E_d.x[i1][j1][k1], m.E_d.x[i2][j2][k2]);
            m.D_u.x[i][j][k] = clamp_hydro(m.c43 * m.D_u.x[i1][j1][k1] - m.c13 * m.D_u.x[i2][j2][k2],
                                            m.D_u.x[i1][j1][k1], m.D_u.x[i2][j2][k2]);
            m.D_d.x[i][j][k] = clamp_hydro(m.c43 * m.D_d.x[i1][j1][k1] - m.c13 * m.D_d.x[i2][j2][k2],
                                            m.D_d.x[i1][j1][k1], m.D_d.x[i2][j2][k2]);
            m.M_u.x[i][j][k] = clamp_hydro(m.c43 * m.M_u.x[i1][j1][k1] - m.c13 * m.M_u.x[i2][j2][k2],
                                            m.M_u.x[i1][j1][k1], m.M_u.x[i2][j2][k2]);
            m.M_d.x[i][j][k] = clamp_hydro(m.c43 * m.M_d.x[i1][j1][k1] - m.c13 * m.M_d.x[i2][j2][k2],
                                            m.M_d.x[i1][j1][k1], m.M_d.x[i2][j2][k2]);
        };

        // Zero-gradient (Neumann) copy of the same six fields from one cell to another.
        // Used at grid boundaries (poles, Greenwich seam, top layer) where a second
        // air neighbour for extrapolation is not available.
        auto copy = [&](int i, int j, int k,
                         int i1, int j1, int k1) {
            m.tke.x[i][j][k] = m.tke.x[i1][j1][k1];
            m.dis.x[i][j][k] = m.dis.x[i1][j1][k1];
            m.nue.x[i][j][k] = m.nue.x[i1][j1][k1];

            m.t.x[i][j][k]     = m.t.x[i1][j1][k1];

            m.p_dyn.x[i][j][k] = m.p_dyn.x[i1][j1][k1];

            m.r_dry.x[i][j][k] = m.r_dry.x[i1][j1][k1];
            m.r_humid.x[i][j][k] = m.r_humid.x[i1][j1][k1];

            m.c.x[i][j][k]     = m.c.x[i1][j1][k1];
            m.cloud.x[i][j][k] = m.cloud.x[i1][j1][k1];
            m.ice.x[i][j][k]   = m.ice.x[i1][j1][k1];
            m.gr.x[i][j][k]    = m.gr.x[i1][j1][k1];

            m.co2.x[i][j][k]   = m.co2.x[i1][j1][k1];

            m.S_c_c.x[i][j][k] = m.S_c_c.x[i1][j1][k1];
            m.S_v.x[i][j][k]   = m.S_v.x[i1][j1][k1];
            m.S_c.x[i][j][k]   = m.S_c.x[i1][j1][k1];
            m.S_i.x[i][j][k]   = m.S_i.x[i1][j1][k1];
            m.S_r.x[i][j][k]   = m.S_r.x[i1][j1][k1];
            m.S_s.x[i][j][k]   = m.S_s.x[i1][j1][k1];

            m.q_v_u.x[i][j][k] = m.q_v_u.x[i1][j1][k1];
            m.q_c_u.x[i][j][k] = m.q_c_u.x[i1][j1][k1];
            m.q_v_d.x[i][j][k] = m.q_v_d.x[i1][j1][k1];

            m.E_u.x[i][j][k] = m.E_u.x[i1][j1][k1];
            m.E_d.x[i][j][k] = m.E_d.x[i1][j1][k1];
            m.D_u.x[i][j][k] = m.D_u.x[i1][j1][k1];
            m.D_d.x[i][j][k] = m.D_d.x[i1][j1][k1];
            m.M_u.x[i][j][k] = m.M_u.x[i1][j1][k1];
            m.M_d.x[i][j][k] = m.M_d.x[i1][j1][k1];
        };

        // Pass 1 — apply solid-wall conditions to every land cell from i=0 up to
        // and including the topographic surface layer i_mount.
        // Two distinct cases are handled in separate branches:
        //   interior : i < i_mount AND no air neighbour in j or k
        //              (fully buried — cannot affect the air column in any direction)
        //   surface  : i == i_mount OR adjacent to air in j or k
        //              (topographic top or cliff face — boundary to the air column)
        #pragma omp parallel for schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                int i_mount = m.i_topography[j][k];

                for (int i = 0; i <= i_mount; i++) {
                    if (!is_land(m.h, i, j, k)) continue;

                    // Check for air neighbours in j (latitude) and k (longitude).
                    // k is periodic: the Greenwich seam (k=0) wraps to k=km-1.
                    bool adj_air_j =
                        (j > 0       && is_air(m.h, i, j-1, k)) ||
                        (j < m.jm-1  && is_air(m.h, i, j+1, k));
                    int k_prev = (k > 0) ? k-1 : m.km-1;
                    int k_next = (k < m.km-1) ? k+1 : 0;
                    bool adj_air_k = is_air(m.h, i, j, k_prev) || is_air(m.h, i, j, k_next);

                    bool on_surface = (i == i_mount) || adj_air_j || adj_air_k;

                    if (!on_surface) {
                        // Fully buried interior cell — zero every prognostic field so
                        // the orography cannot contaminate the overlying air column.

                        // No-slip velocity (current and next time level).
                        m.u.x[i][j][k]  = 0.0;
                        m.v.x[i][j][k]  = 0.0;
                        m.w.x[i][j][k]  = 0.0;
                        m.un.x[i][j][k] = 0.0;
                        m.vn.x[i][j][k] = 0.0;
                        m.wn.x[i][j][k] = 0.0;

                        m.t.x[i][j][k]     = 1.0;
                        m.p_dyn.x[i][j][k] = 0.0;

                        // Density reset to dry-air reference; CO2 to background.
                        m.r_dry.x[i][j][k]   = m.r_air;
                        m.r_humid.x[i][j][k] = m.r_air;
                        m.co2.x[i][j][k]     = m.co2_0;

                        // Microphysics (current and time-level n+1).
                        m.c.x[i][j][k]      = 0.0;
                        m.cloud.x[i][j][k]  = 0.0;
                        m.ice.x[i][j][k]    = 0.0;
                        m.gr.x[i][j][k]     = 0.0;
                        m.cn.x[i][j][k]     = 0.0;
                        m.cloudn.x[i][j][k] = 0.0;
                        m.icen.x[i][j][k]   = 0.0;
                        m.grn.x[i][j][k]    = 0.0;

                        // Microphysical source/sink terms.
                        m.S_c_c.x[i][j][k] = 0.0;
                        m.S_v.x[i][j][k]   = 0.0;
                        m.S_c.x[i][j][k]   = 0.0;
                        m.S_i.x[i][j][k]   = 0.0;
                        m.S_r.x[i][j][k]   = 0.0;
                        m.S_s.x[i][j][k]   = 0.0;

                        // All body forces zero inside the solid body.
                        m.BuoyancyForce.x[i][j][k]    = 0.0;
                        m.PresGradForce.x[i][j][k]    = 0.0;
                        m.CoriolisForce.x[i][j][k]    = 0.0;
                        m.CentrifugalForce.x[i][j][k] = 0.0;

                        m.Q_Latent.x[i][j][k]   = 0.0;
                        m.Q_Sensible.x[i][j][k] = 0.0;

                        // Turbulence: zero all scalars inside solid cells.
                        m.tke.x[i][j][k]        = 0.0;
                        m.tken.x[i][j][k]       = 0.0;
                        m.dis.x[i][j][k]        = 0.0;
                        m.disn.x[i][j][k]       = 0.0;
                        m.nue.x[i][j][k]        = 0.0;
                        m.prod.x[i][j][k]       = 0.0;
                        m.tke_source.x[i][j][k] = 0.0;
                        m.dis_source.x[i][j][k] = 0.0;

                        // Moist-convection scalars: dry static energy set to reference (1),
                        // all convective velocities and moisture quantities zeroed.
                        m.s.x[i][j][k]   = 1.0;
                        m.s_u.x[i][j][k] = 1.0;
                        m.s_d.x[i][j][k] = 1.0;
                        m.u_u.x[i][j][k] = 0.0;
                        m.u_d.x[i][j][k] = 0.0;
                        m.v_u.x[i][j][k] = 0.0;
                        m.v_d.x[i][j][k] = 0.0;
                        m.w_u.x[i][j][k] = 0.0;
                        m.w_d.x[i][j][k] = 0.0;

                        // Convective condensate and precipitation conversion terms.
                        m.c_u.x[i][j][k] = 0.0;
                        m.e_d.x[i][j][k] = 0.0;
                        m.e_l.x[i][j][k] = 0.0;
                        m.e_p.x[i][j][k] = 0.0;
                        m.g_p.x[i][j][k] = 0.0;

                        m.q_c_u.x[i][j][k] = 0.0;
                        m.q_v_u.x[i][j][k] = 0.0;
                        m.q_v_d.x[i][j][k] = 0.0;

                        // Mass-flux and entrainment/detrainment rates.
                        m.D_d.x[i][j][k] = 0.0;
                        m.D_u.x[i][j][k] = 0.0;
                        m.E_d.x[i][j][k] = 0.0;
                        m.E_u.x[i][j][k] = 0.0;
                        m.M_d.x[i][j][k] = 0.0;
                        m.M_u.x[i][j][k] = 0.0;

                        // Moist-convection tendencies (temperature, humidity, momentum).
                        m.MC_q.x[i][j][k] = 0.0;
                        m.MC_t.x[i][j][k] = 0.0;
                        m.MC_v.x[i][j][k] = 0.0;
                        m.MC_w.x[i][j][k] = 0.0;
                    } else {
                        // Surface cell — either the topographic top (i == i_mount) or a
                        // cliff face (i < i_mount but adjacent to air in j or k).

                        // No-slip velocity (current and next time level).
                        m.u.x[i][j][k]  = 0.0;
                        m.v.x[i][j][k]  = 0.0;
                        m.w.x[i][j][k]  = 0.0;
                        m.un.x[i][j][k] = 0.0;
                        m.vn.x[i][j][k] = 0.0;
                        m.wn.x[i][j][k] = 0.0;

                        m.t.x[i][j][k]     = 1.0;
                        m.p_dyn.x[i][j][k] = 0.0;

                        // Density reset to dry-air reference; CO2 to background.
                        m.r_dry.x[i][j][k]   = m.r_air;
                        m.r_humid.x[i][j][k] = m.r_air;
                        m.co2.x[i][j][k]     = m.co2_0;

                        // Microphysics (current and time-level n+1).
                        m.c.x[i][j][k]      = 0.0;
                        m.cloud.x[i][j][k]  = 0.0;
                        m.ice.x[i][j][k]    = 0.0;
                        m.gr.x[i][j][k]     = 0.0;
                        m.cn.x[i][j][k]     = 0.0;
                        m.cloudn.x[i][j][k] = 0.0;
                        m.icen.x[i][j][k]   = 0.0;
                        m.grn.x[i][j][k]    = 0.0;

                        // Microphysical source/sink terms.
                        m.S_c_c.x[i][j][k] = 0.0;
                        m.S_v.x[i][j][k]   = 0.0;
                        m.S_c.x[i][j][k]   = 0.0;
                        m.S_i.x[i][j][k]   = 0.0;
                        m.S_r.x[i][j][k]   = 0.0;
                        m.S_s.x[i][j][k]   = 0.0;

                        // All body forces zero inside the solid body.
                        m.BuoyancyForce.x[i][j][k]    = 0.0;
                        m.PresGradForce.x[i][j][k]    = 0.0;
                        m.CoriolisForce.x[i][j][k]    = 0.0;
                        m.CentrifugalForce.x[i][j][k] = 0.0;

                        m.Q_Latent.x[i][j][k]   = 0.0;
                        m.Q_Sensible.x[i][j][k] = 0.0;

                        // Turbulence: zero all scalars at solid surfaces.
                        m.tke.x[i][j][k]        = 0.0;
                        m.tken.x[i][j][k]       = 0.0;
                        m.dis.x[i][j][k]        = 0.0;
                        m.disn.x[i][j][k]       = 0.0;
                        m.nue.x[i][j][k]        = 0.0;
                        m.prod.x[i][j][k]       = 0.0;
                        m.tke_source.x[i][j][k] = 0.0;
                        m.dis_source.x[i][j][k] = 0.0;

                        // Moist-convection scalars: dry static energy set to reference (1),
                        // all convective velocities and moisture quantities zeroed.
                        m.s.x[i][j][k]   = 1.0;
                        m.s_u.x[i][j][k] = 1.0;
                        m.s_d.x[i][j][k] = 1.0;
                        m.u_u.x[i][j][k] = 0.0;
                        m.u_d.x[i][j][k] = 0.0;
                        m.v_u.x[i][j][k] = 0.0;
                        m.v_d.x[i][j][k] = 0.0;
                        m.w_u.x[i][j][k] = 0.0;
                        m.w_d.x[i][j][k] = 0.0;

                        // Convective condensate and precipitation conversion terms.
                        m.c_u.x[i][j][k] = 0.0;
                        m.e_d.x[i][j][k] = 0.0;
                        m.e_l.x[i][j][k] = 0.0;
                        m.e_p.x[i][j][k] = 0.0;
                        m.g_p.x[i][j][k] = 0.0;

                        m.q_c_u.x[i][j][k] = 0.0;
                        m.q_v_u.x[i][j][k] = 0.0;
                        m.q_v_d.x[i][j][k] = 0.0;

                        // Mass-flux and entrainment/detrainment rates.
                        m.D_d.x[i][j][k] = 0.0;
                        m.D_u.x[i][j][k] = 0.0;
                        m.E_d.x[i][j][k] = 0.0;
                        m.E_u.x[i][j][k] = 0.0;
                        m.M_d.x[i][j][k] = 0.0;
                        m.M_u.x[i][j][k] = 0.0;

                        // Moist-convection tendencies (temperature, humidity, momentum).
                        m.MC_q.x[i][j][k] = 0.0;
                        m.MC_t.x[i][j][k] = 0.0;
                        m.MC_v.x[i][j][k] = 0.0;
                        m.MC_w.x[i][j][k] = 0.0;
                    }
                }
            }
        }

        // Pass 2 — extrapolate p_dyn, microphysics, and CO2 from air into land ghost
        // cells at every land/air interface in the j- and k-directions, and at the
        // uppermost land level in the i-direction.  This gives the pressure-gradient
        // and tracer-advection stencils smooth values across the topography boundary.
        #pragma omp parallel for schedule(static)
        for (int i = 1; i < m.im-1; i++) {
            for (int j = 1; j < m.jm-1; j++) {
                for (int k = 1; k < m.km-1; k++) {

                    bool land_ijk = is_land(m.h, i, j, k);

                    // i-direction (radial): extrapolate outward from the first two air
                    // cells above the surface.  At the second-to-last level use copy
                    // because there is only one air cell available above.
                    if (i < m.im-3) {
                        if (land_ijk && is_air(m.h, i+1, j, k))
                            extrapolation(i,j,k, i+1,j,k, i+2,j,k);
                    } else if (i == m.im-2) {
                        if (land_ijk && is_air(m.h, i+1, j, k))
                            copy(i,j,k, i+1,j,k);
                    }

                    // j-direction interior (latitude): extrapolate from whichever
                    // side has two consecutive air cells.
                    if (j > 2 && j < m.jm-3) {
                        if (land_ijk && is_air(m.h, i, j+1, k) && is_air(m.h, i, j+2, k))
                            extrapolation(i,j,k, i,j+1,k, i,j+2,k);
                        if (land_ijk && is_air(m.h, i, j-1, k) && is_air(m.h, i, j-2, k))
                            extrapolation(i,j,k, i,j-1,k, i,j-2,k);
                    }

                    // j-direction poles: single neighbour only, use copy.
                    if (j == 0)        copy(i,j,k, i,j+1,k);
                    if (j == m.jm-1)   copy(i,j,k, i,j-1,k);

                    // k-direction interior (longitude): same logic as j-direction.
                    if (k > 2 && k < m.km-3) {
                        if (land_ijk && is_air(m.h, i, j, k+1) && is_air(m.h, i, j, k+2))
                            extrapolation(i,j,k, i,j,k+1, i,j,k+2);
                        if (land_ijk && is_air(m.h, i, j, k-1) && is_air(m.h, i, j, k-2))
                            extrapolation(i,j,k, i,j,k-1, i,j,k-2);
                    }

                    // k-direction Greenwich seam (k=0 wraps to k=km-1): copy from inner neighbour.
                    if (k == 0)        copy(i,j,k, i,j,k+1);
                    if (k == m.km-1)   copy(i,j,k, i,j,k-1);
                }
            }
        }

        // Pass 3 — copy surface-level (i_mount) values of temperature, microphysics,
        // microphysical sources, and precipitation down to the i=0 reference layer.
        // This keeps the sea-level layer consistent with the lowest atmospheric level
        // so that surface fluxes and output diagnostics see the correct surface state.
        #pragma omp parallel for schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                int i_mount = m.i_topography[j][k];

                if (std::isfinite(m.t.x[i_mount][j][k]))
                    m.t.x[0][j][k] = m.t.x[i_mount][j][k];
 
                m.p_dyn.x[0][j][k] = m.p_dyn.x[i_mount][j][k];

                m.r_dry.x[0][j][k] = m.r_dry.x[i_mount][j][k];
                m.r_humid.x[0][j][k] = m.r_humid.x[i_mount][j][k];

                m.c.x[0][j][k]     = m.c.x[i_mount][j][k];
                m.cloud.x[0][j][k] = m.cloud.x[i_mount][j][k];
                m.ice.x[0][j][k]   = m.ice.x[i_mount][j][k];
                m.gr.x[0][j][k]    = m.gr.x[i_mount][j][k];

                m.co2.x[0][j][k]   = m.co2.x[i_mount][j][k];

                m.S_c_c.x[0][j][k] = m.S_c_c.x[i_mount][j][k];
                m.S_v.x[0][j][k]   = m.S_v.x[i_mount][j][k];
                m.S_c.x[0][j][k]   = m.S_c.x[i_mount][j][k];
                m.S_i.x[0][j][k]   = m.S_i.x[i_mount][j][k];
                m.S_r.x[0][j][k]   = m.S_r.x[i_mount][j][k];
                m.S_s.x[0][j][k]   = m.S_s.x[i_mount][j][k];

                m.q_v_u.x[0][j][k] = m.q_v_u.x[i_mount][j][k];
                m.q_c_u.x[0][j][k] = m.q_c_u.x[i_mount][j][k];
                m.q_v_d.x[0][j][k] = m.q_v_d.x[i_mount][j][k];

                m.E_u.x[0][j][k] = m.E_u.x[i_mount][j][k];
                m.E_d.x[0][j][k] = m.E_d.x[i_mount][j][k];
                m.D_u.x[0][j][k] = m.D_u.x[i_mount][j][k];
                m.D_d.x[0][j][k] = m.D_d.x[i_mount][j][k];
                m.M_u.x[0][j][k] = m.M_u.x[i_mount][j][k];
                m.M_d.x[0][j][k] = m.M_d.x[i_mount][j][k];

                m.tke.x[0][j][k] = m.tke.x[i_mount][j][k];
                m.dis.x[0][j][k] = m.dis.x[i_mount][j][k];
            }
        }

        // Pass 4 — clip any negative microphysics values to zero across the full
        // domain.  Negative values can arise from the advection scheme near sharp
        // gradients; clamping here prevents unphysical condensate.
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                for (int k = 0; k < m.km; k++) {
                    if (m.c.x[i][j][k]     < 0.0)  m.c.x[i][j][k]     = 0.0;
                    if (m.cloud.x[i][j][k] < 0.0)  m.cloud.x[i][j][k] = 0.0;
                    if (m.ice.x[i][j][k]   < 0.0)  m.ice.x[i][j][k]   = 0.0;
                    if (m.tke.x[i][j][k]   < 0.0)  m.tke.x[i][j][k]   = 0.0;
                    if (m.dis.x[i][j][k]   < 0.0)  m.dis.x[i][j][k]   = 0.0;
                }
            }
        }
        cout << "      AGCM: BC_SolidGround ended" << endl;
    }
/*
* 
*/
    // ------------------------------------------------------------------
    void bcRadius()
    {
        // cubic extrapolation: x[a] = x[a+3d] - 3*x[a+2d] + 3*x[a+d]
        // von Neumann:         x[a] = c43*x[a+d] - c13*x[a+2d]

        // Pattern A: cubic at both i=0 and i=im-1
        Array* both_cubic[] = {
            &m.p_stat, &m.r_humid, &m.r_dry,
            &m.PrecipitableWaterLocal,
            &m.CoriolisForce, &m.CentrifugalForce, &m.BuoyancyForce, &m.PresGradForce,
            &m.Q_Latent, &m.Q_Sensible,
            &m.S_c_c, &m.S_v, &m.S_c, &m.S_i, &m.S_r, &m.S_s, &m.S_g,
            &m.q_v_u, &m.q_c_u, &m.u_u, &m.v_u, &m.w_u,
            &m.s, &m.s_u, &m.s_d,
            &m.MC_t, &m.MC_q, &m.MC_v, &m.MC_w,
            &m.M_u, &m.M_d
        };
        constexpr int n_both = sizeof(both_cubic) / sizeof(both_cubic[0]);

        // Pattern B: von Neumann at i=0, cubic at i=im-1
        Array* vn_bot_cubic_top[] = { 
            &m.u, &m.v, &m.w, 
            &m.P_conv 
//            &m.tke, &m.dis 
            };
        constexpr int n_vn_cubic = sizeof(vn_bot_cubic_top) / sizeof(vn_bot_cubic_top[0]);

        // Pattern C: cubic at i=0, von Neumann at i=im-1
        // Turbulence scalars use this pattern: at i=im-1 turbulence decays to its
        // background floor, making the profile concave. Cubic extrapolation there
        // amplifies the concavity and produces values much larger than nue[im-2].
        // Von Neumann (zero gradient) is the correct physical BC at the model top.
        // At i=0 cubic is used; apply_wall_bc() reasserts the correct dis value.
//        Array* cubic_bot_vn_top[] = { &m.c, &m.cloud, &m.ice, &m.gr, &m.co2,
        Array* cubic_bot_vn_top[] = {
            &m.cloud, &m.ice, &m.gr, &m.co2,
            &m.c_u, &m.e_d, &m.e_l, &m.e_p, &m.g_p,
            &m.tke, &m.dis, &m.nue, &m.prod, &m.tke_source, &m.dis_source
        };
        constexpr int n_cubic_vn = sizeof(cubic_bot_vn_top) / sizeof(cubic_bot_vn_top[0]);

        const int iml = m.im - 1;

        #pragma omp parallel for schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {

                // t: no extrapolation at i=0, cubic at i=im-1
                m.t.x[iml][j][k] = m.t.x[iml-3][j][k]
                    - 3.0 * m.t.x[iml-2][j][k] + 3.0 * m.t.x[iml-1][j][k];

                // Pattern A
                for (int f = 0; f < n_both; f++) {
                    double*** xf = both_cubic[f]->x;
                    xf[0][j][k]   = xf[3][j][k]     - 3.0 * xf[2][j][k]     + 3.0 * xf[1][j][k];
                    xf[iml][j][k] = xf[iml-3][j][k] - 3.0 * xf[iml-2][j][k] + 3.0 * xf[iml-1][j][k];
                }

                // Pattern B
                for (int f = 0; f < n_vn_cubic; f++) {
                    double*** xf = vn_bot_cubic_top[f]->x;
                    xf[0][j][k]   = m.c43 * xf[1][j][k]     - m.c13 * xf[2][j][k];
                    xf[iml][j][k] = xf[iml-3][j][k] - 3.0 * xf[iml-2][j][k] + 3.0 * xf[iml-1][j][k];
                }

                // Pattern C
                for (int f = 0; f < n_cubic_vn; f++) {
                    double*** xf = cubic_bot_vn_top[f]->x;
                    xf[0][j][k]   = xf[3][j][k]     - 3.0 * xf[2][j][k]     + 3.0 * xf[1][j][k];
                    xf[iml][j][k] = m.c43 * xf[iml-1][j][k] - m.c13 * xf[iml-2][j][k];
                }

                // updraft moisture at the surface must be non-negative
                if (m.q_v_u.x[0][j][k] < 0.0) m.q_v_u.x[0][j][k] = 0.0;
                if (m.q_c_u.x[0][j][k] < 0.0) m.q_c_u.x[0][j][k] = 0.0;
            }
        }
    }
/*
* 
*/
    // ------------------------------------------------------------------
    void bcTheta()
    {
        // von Neumann at both j=0 and j=jm-1
        // Turbulence scalars use this too: the cubic formula amplifies concave
        // profiles near the pole singularity (sinθ→0), producing ghost-cell
        // values larger than the interior and inflating prod and D_w.
        Array* fields_vn[] = {
            &m.t, &m.p_stat, &m.r_humid, &m.r_dry,
            &m.u, &m.v, &m.w,
            &m.c, &m.cloud, &m.ice, &m.gr,
            &m.PrecipitableWaterLocal, &m.co2,
            &m.tke, &m.dis, &m.nue, &m.prod, &m.tke_source, &m.dis_source
        };
        constexpr int n_vn = sizeof(fields_vn) / sizeof(fields_vn[0]);

        // cubic extrapolation at both j=0 and j=jm-1
        Array* fields_cubic[] = {
            &m.S_c_c, &m.S_v, &m.S_c, &m.S_i, &m.S_r, &m.S_s, &m.S_g,
            &m.q_v_u, &m.q_c_u, &m.u_u, &m.v_u, &m.w_u,
            &m.s, &m.s_u, &m.s_d,
            &m.CoriolisForce, &m.CentrifugalForce, &m.BuoyancyForce, &m.PresGradForce,
            &m.Q_Latent, &m.Q_Sensible,
            &m.c_u, &m.e_d, &m.e_l, &m.e_p, &m.g_p
        };
        constexpr int n_cubic = sizeof(fields_cubic) / sizeof(fields_cubic[0]);

        const int jml = m.jm - 1;

        #pragma omp parallel for schedule(static)
        for (int k = 0; k < m.km; k++) {
            for (int i = 0; i < m.im; i++) {

                for (int f = 0; f < n_vn; f++) {
                    double*** xf = fields_vn[f]->x;
                    xf[i][0][k]   = m.c43 * xf[i][1][k]     - m.c13 * xf[i][2][k];
                    xf[i][jml][k] = m.c43 * xf[i][jml-1][k] - m.c13 * xf[i][jml-2][k];
                }

                for (int f = 0; f < n_cubic; f++) {
                    double*** xf = fields_cubic[f]->x;
                    xf[i][0][k]   = xf[i][3][k]     - 3.0 * xf[i][2][k]     + 3.0 * xf[i][1][k];
                    xf[i][jml][k] = xf[i][jml-3][k] - 3.0 * xf[i][jml-2][k] + 3.0 * xf[i][jml-1][k];
                }
            }
        }
    }
/*
* 
*/
    // ------------------------------------------------------------------
    void bcPhi()
    {
        // fields: von Neumann extrapolation + periodic averaging
        Array* fields_avg[] = {
            &m.t, &m.p_stat, &m.r_humid, &m.r_dry,
            &m.u, &m.v, &m.w,
            &m.c, &m.cloud, &m.ice, &m.gr, &m.co2,
            &m.S_c_c, &m.S_v, &m.S_c, &m.S_i, &m.S_r, &m.S_s, &m.S_g,
            &m.q_v_u, &m.q_c_u, &m.u_u, &m.v_u, &m.w_u,
            &m.s, &m.s_u, &m.s_d,
            &m.CoriolisForce, &m.CentrifugalForce, &m.PresGradForce, &m.BuoyancyForce,
            &m.Q_Latent, &m.Q_Sensible,
            &m.c_u, &m.e_d, &m.e_l, &m.e_p, &m.g_p
//            &m.tke, &m.dis, &m.nue, &m.prod, &m.tke_source, &m.dis_source
        };
        constexpr int n_avg = sizeof(fields_avg) / sizeof(fields_avg[0]);

        // fields: von Neumann only (no periodic averaging)
        Array* fields_extrap[] = { 
            &m.PrecipitableWaterLocal, 
            &m.tke, &m.dis, &m.nue, &m.prod, &m.tke_source, &m.dis_source
        };
        constexpr int n_extrap = sizeof(fields_extrap) / sizeof(fields_extrap[0]);

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {

                for (int f = 0; f < n_avg; f++) {
                    double** xij = fields_avg[f]->x[i];
                    double v0   = m.c43 * xij[j][1]      - m.c13 * xij[j][2];
                    double vend = m.c43 * xij[j][m.km-2] - m.c13 * xij[j][m.km-3];
                    xij[j][0] = xij[j][m.km-1] = (v0 + vend) * 0.5;
                }

                for (int f = 0; f < n_extrap; f++) {
                    double** xij = fields_extrap[f]->x[i];
                    xij[j][0]      = m.c43 * xij[j][1]      - m.c13 * xij[j][2];
                    xij[j][m.km-1] = m.c43 * xij[j][m.km-2] - m.c13 * xij[j][m.km-3];
                }
            }
        }
    }
/*
* 
*/
    // ------------------------------------------------------------------
    void bcVelSurfSur()
    {
        using namespace std;
        cout << endl << "      AGCM: BC_vel_surf_sur" << endl;


        constexpr double coeff  = 0.01;                                 // represents the flow better
        constexpr double coeff5 = 0.9;                                  // represents the flow better

        #pragma omp parallel for schedule(static)
        for (int i = 1; i < m.im-1; i++) {
            for (int j = 1; j < m.jm-1; j++) {
                for (int k = 1; k < m.km-1; k++) {

                    bool land_ijk = is_land(m.h, i, j, k);

                    // i-direction (radial: land→air outward)
                    if (i < m.im-2 && land_ijk && is_air(m.h, i+1, j, k)) {
                        m.u.x[i][j][k]   = 0.0;
                        m.u.x[i+1][j][k] = 0.0;

                        m.v.x[i][j][k]    = 0.0;
                        m.v.x[i+1][j][k] *= coeff;
                        m.v.x[i+2][j][k] *= coeff5;

                        m.w.x[i][j][k]    = 0.0;
                        m.w.x[i+1][j][k] *= coeff;
                        m.w.x[i+2][j][k] *= coeff5;
                    }

                    // j-direction (meridional: land→air south/north)
                    if (land_ijk) {
                        if (j < m.jm-2 && is_air(m.h, i, j+1, k) && is_air(m.h, i, j+2, k)) {  // south
                            m.u.x[i][j][k]    = 0.0;
                            m.u.x[i][j+1][k] *= coeff;
                            m.u.x[i][j+2][k] *= coeff5;

                            m.v.x[i][j][k]    = 0.0;
                            m.v.x[i][j+1][k]  = 0.0;

                            m.w.x[i][j][k]    = 0.0;
                            m.w.x[i][j+1][k] *= coeff;
                            m.w.x[i][j+2][k] *= coeff5;
                        }

                        if (j >= 2 && is_air(m.h, i, j-1, k) && is_air(m.h, i, j-2, k)) {  // north
                            m.u.x[i][j][k]    = 0.0;
                            m.u.x[i][j-1][k] *= coeff;
                            m.u.x[i][j-2][k] *= coeff5;

                            m.v.x[i][j][k]    = 0.0;
                            m.v.x[i][j-1][k]  = 0.0;

                            m.w.x[i][j][k]    = 0.0;
                            m.w.x[i][j-1][k] *= coeff;
                            m.w.x[i][j-2][k] *= coeff5;
                        }
                    }

                    // k-direction (zonal: land→air east/west)
                    if (land_ijk) {
                        if (k < m.km-2 && is_air(m.h, i, j, k+1) && is_air(m.h, i, j, k+2)) {  // east
                            m.u.x[i][j][k]    = 0.0;
                            m.u.x[i][j][k+1] *= coeff;
                            m.u.x[i][j][k+2] *= coeff5;

                            m.v.x[i][j][k]    = 0.0;
                            m.v.x[i][j][k+1] *= coeff;
                            m.v.x[i][j][k+2] *= coeff5;

                            m.w.x[i][j][k]    = 0.0;
                            m.w.x[i][j][k+1]  = 0.0;
                        }

                        if (k >= 2 && is_air(m.h, i, j, k-1) && is_air(m.h, i, j, k-2)) {  // west
                            m.u.x[i][j][k]    = 0.0;
                            m.u.x[i][j][k-1] *= coeff;
                            m.u.x[i][j][k-2] *= coeff5;

                            m.v.x[i][j][k]    = 0.0;
                            m.v.x[i][j][k-1] *= coeff;
                            m.v.x[i][j][k-2] *= coeff5;

                            m.w.x[i][j][k]    = 0.0;
                            m.w.x[i][j][k-1]  = 0.0;
                        }
                    }
                }
            }
        }
        cout << "      AGCM: BC_vel_surf_sur ended" << endl;
    }
/*
*
*/
    // ------------------------------------------------------------------
    // Neumann (zero normal gradient) boundary condition for all scalar
    // 3-D fields at land-air interfaces.  For each interior land cell the
    // value is extrapolated from the two nearest air cells in the normal
    // direction: x[land] = c43*x[n1] - c13*x[n2].  Velocity components
    // from MoistConvection (u_u, v_u, w_u, u_d, v_d, w_d, MC_v, MC_w)
    // are excluded; the main-model velocities are handled by bcVelSurfSur.
    // At corners where multiple directions border air, the vertical (i)
    // direction is preferred.
    void bcScalarSurfSur()
    {
        using namespace std;
        cout << endl << "      AGCM: BC_scalar_surf_sur" << endl;

        // S-functions (S_v, S_c, S_i, S_r, S_s, S_g, S_c_c), precipitation
        // (P_rain, P_snow, Precipitation), and p_dyn are excluded: bcSolidGround
        // zeros them at land cells and projects x[i_mount] to i=0, so Neumann
        // extrapolation here would only produce spuriously high land-cell values.
        Array* scalars[] = {
            &m.t,
            &m.c, &m.cloud, &m.ice, &m.gr, &m.co2,
            &m.p_stat,
            &m.r_dry, &m.r_humid,
            &m.Q_Latent, &m.Q_Sensible,
            &m.BuoyancyForce, &m.CoriolisForce, &m.CentrifugalForce, &m.PresGradForce,
            &m.PrecipitableWaterLocal,
            &m.M_u,   &m.M_d,
            &m.MC_t,  &m.MC_q,
            &m.s,     &m.s_u,   &m.s_d,
            &m.q_v_u, &m.q_v_d, &m.q_c_u,
            &m.c_u,   &m.e_d,   &m.e_l,  &m.e_p,  &m.g_p,
            &m.E_u,   &m.D_u,   &m.E_d,  &m.D_d
        };
        constexpr int ns = sizeof(scalars) / sizeof(scalars[0]);

        #pragma omp parallel for schedule(static)
        for (int i = 1; i < m.im-1; i++) {
            for (int j = 1; j < m.jm-1; j++) {
                for (int k = 1; k < m.km-1; k++) {

                    if (!is_land(m.h, i, j, k)) continue;

                    // ---- i-direction (vertical, preferred at corners) ----
                    if (is_air(m.h, i+1, j, k)) {
                        if (i + 2 < m.im) {
                            for (int f = 0; f < ns; f++) {
                                double*** x = scalars[f]->x;
                                x[i][j][k] = m.c43 * x[i+1][j][k] - m.c13 * x[i+2][j][k];
                            }
                        } else {
                            for (int f = 0; f < ns; f++)
                                scalars[f]->x[i][j][k] = scalars[f]->x[i+1][j][k];
                        }
                        continue;                                    // i preferred: skip j and k
                    }

                    // ---- j-direction ----
                    if (j < m.jm-2 && is_air(m.h, i, j+1, k) && is_air(m.h, i, j+2, k)) {
                        for (int f = 0; f < ns; f++) {
                            double*** x = scalars[f]->x;
                            x[i][j][k] = m.c43 * x[i][j+1][k] - m.c13 * x[i][j+2][k];
                        }
                    } else if (j >= 2 && is_air(m.h, i, j-1, k) && is_air(m.h, i, j-2, k)) {
                        for (int f = 0; f < ns; f++) {
                            double*** x = scalars[f]->x;
                            x[i][j][k] = m.c43 * x[i][j-1][k] - m.c13 * x[i][j-2][k];
                        }
                    }

                    // ---- k-direction ----
                    if (k < m.km-2 && is_air(m.h, i, j, k+1) && is_air(m.h, i, j, k+2)) {
                        for (int f = 0; f < ns; f++) {
                            double*** x = scalars[f]->x;
                            x[i][j][k] = m.c43 * x[i][j][k+1] - m.c13 * x[i][j][k+2];
                        }
                    } else if (k >= 2 && is_air(m.h, i, j, k-1) && is_air(m.h, i, j, k-2)) {
                        for (int f = 0; f < ns; f++) {
                            double*** x = scalars[f]->x;
                            x[i][j][k] = m.c43 * x[i][j][k-1] - m.c13 * x[i][j][k-2];
                        }
                    }
                }
            }
        }

        #pragma omp parallel for schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                int i_mount = m.i_topography[j][k];

                if (std::isfinite(m.t.x[i_mount][j][k]))
                    m.t.x[0][j][k] = m.t.x[i_mount][j][k];
                m.c.x[0][j][k]     = m.c.x[i_mount][j][k];
                m.cloud.x[0][j][k] = m.cloud.x[i_mount][j][k];
                m.ice.x[0][j][k]   = m.ice.x[i_mount][j][k];
                m.gr.x[0][j][k]    = m.gr.x[i_mount][j][k];

                m.S_c_c.x[0][j][k] = m.S_c_c.x[i_mount][j][k];
                m.S_v.x[0][j][k]   = m.S_v.x[i_mount][j][k];
                m.S_c.x[0][j][k]   = m.S_c.x[i_mount][j][k];
                m.S_i.x[0][j][k]   = m.S_i.x[i_mount][j][k];
                m.S_r.x[0][j][k]   = m.S_r.x[i_mount][j][k];
                m.S_s.x[0][j][k]   = m.S_s.x[i_mount][j][k];

                m.P_rain.x[0][j][k] = m.P_rain.x[i_mount][j][k];
                m.P_snow.x[0][j][k] = m.P_snow.x[i_mount][j][k];

                m.p_dyn.x[0][j][k] = m.p_dyn.x[i_mount][j][k];
            }
        }

        cout << "      AGCM: BC_scalar_surf_sur ended" << endl;
    }
/*
* 
*/
    // ------------------------------------------------------------------
    void coastalCurrents()
    {
        using namespace std;
        cout << endl << "      AGCM: CoastalCurrents" << endl;

        constexpr int l_end = 20;

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < m.im-1; i++) {
            for (int j = 1; j < m.jm-1; j++) {
                for (int k = 1; k < m.km-1; k++) {

                    if (!is_land(m.h, i, j, k)) continue;

                    if (is_land(m.h, i, j, k-1) && is_air(m.h, i, j, k+1)) {  // east coast
                        int l_max = std::min(l_end, m.km - 1 - k);
                        for (int l = 1; l <= l_max; l++) {
                            m.v.x[i][j][k+l] = -m.v.x[i][j][k+l];
                            m.w.x[i][j][k+l] = -m.w.x[i][j][k+l];
                        }
                    }

                    if (is_land(m.h, i, j, k+1) && is_air(m.h, i, j, k-1)) {  // west coast
                        int l_max = std::min(l_end, k);
                        for (int l = 1; l <= l_max; l++) {
                            m.v.x[i][j][k-l] = -m.v.x[i][j][k-l];
                            m.w.x[i][j][k-l] = -m.w.x[i][j][k-l];
                        }
                    }
                }
            }
        }
        cout << "      AGCM: CoastalCurrents ended" << endl;
    }

private:
    cAtmosphereModel& m;
};
