#pragma once

#include "cAtmosphereModel.h"

// ============================================================================
// Shared machinery for the ice / precipitation microphysics schemes
// (Zero/One/Two/ThreeCatIceScheme). The four schemes share a large amount of
// physics-independent boilerplate — boundary extrapolation, topography fill,
// saturation, the finiteness flux cap — which was previously duplicated (and
// independently mis-fixed) in each file. It is factored out here so a fix lives
// in ONE place. The per-scheme MICROPHYSICS stays in each scheme's
// computeColumns(); only the common scaffolding is shared.
// ============================================================================
namespace IceSchemeCommon {

    // Hard cap on a precipitation flux [kg/(m2*s)] (~260 mm/d, well above any
    // physical precip). Backstop against the ∝Snow riming/deposition runaways.
    constexpr double P_max_flux = 3.0e-3;

    // c43/c13 quadratic edge extrapolation of ONE precipitation-flux field over
    // the top (radial), latitude (theta) and longitude (phi) boundaries. The phi
    // seam optionally gets periodicity averaging: rain uses it in every scheme;
    // snow/graupel use it in the >=2-category schemes but NOT in ZeroCat, so it
    // is a per-field flag to preserve each scheme's exact behaviour.
    inline void extrapolateBC(cAtmosphereModel& m, Array& P, bool periodic_avg) {
        // Top boundary
        #pragma omp parallel for collapse(2)
        for(int j = 0; j < m.jm; j++)
            for(int k = 0; k < m.km; k++)
                P.x[m.im-1][j][k] = m.c43 * P.x[m.im-2][j][k] - m.c13 * P.x[m.im-3][j][k];

        // Latitude (theta) boundaries
        #pragma omp parallel for collapse(2)
        for(int k = 0; k < m.km; k++)
            for(int i = 0; i < m.im; i++){
                P.x[i][0][k]      = m.c43 * P.x[i][1][k]      - m.c13 * P.x[i][2][k];
                P.x[i][m.jm-1][k] = m.c43 * P.x[i][m.jm-2][k] - m.c13 * P.x[i][m.jm-3][k];
            }

        // Longitude (phi) boundaries — von Neumann + optional periodicity average
        #pragma omp parallel for collapse(2)
        for(int i = 0; i < m.im; i++)
            for(int j = 0; j < m.jm; j++){
                P.x[i][j][0]      = m.c43 * P.x[i][j][1]      - m.c13 * P.x[i][j][2];
                P.x[i][j][m.km-1] = m.c43 * P.x[i][j][m.km-2] - m.c13 * P.x[i][j][m.km-3];
                if(periodic_avg)
                    P.x[i][j][0] = P.x[i][j][m.km-1] =
                        (P.x[i][j][0] + P.x[i][j][m.km-1]) * 0.5;
            }
    }

    // Fill the sub-surface land cells (i below i_topography) of ONE field with the
    // surface (i_mount) value. Each field independently copies its own ground value
    // downward, so calling this per field is equivalent to the old snapshot-all form.
    inline void fillTopography(cAtmosphereModel& m, Array& F) {
        #pragma omp parallel for collapse(2)
        for(int j = 0; j < m.jm; j++)
            for(int k = 0; k < m.km; k++){
                int i_mount = m.i_topography[j][k];
                double v = F.x[i_mount][j][k];
                for(int i = i_mount - 1; i >= 0; i--)
                    if(AtomUtils::is_land(m.h, i, j, k)) F.x[i][j][k] = v;
            }
    }

    // Saturation specific humidity over water / ice [kg/kg] (Magnus form used by
    // every scheme). Available for computeColumns to use; kept identical to the
    // inline expressions so a scheme can adopt it without changing results.
    inline double qSatWater(cAtmosphereModel& m, double t_u, int i, int j, int k) {
        const double E = m.hp * AtomUtils::exp_func(t_u, 17.2694, 35.86);
        return m.ep * E / (m.p_stat.x[i][j][k] - E);
    }
    inline double qSatIce(cAtmosphereModel& m, double t_u, int i, int j, int k) {
        const double E = m.hp * AtomUtils::exp_func(t_u, 21.8746, 7.66);
        return m.ep * E / (m.p_stat.x[i][j][k] - E);
    }
}
