// ============================================================================
// rad_selftest — standalone OFFLINE driver for MultiLayerRadiation.
//
// Exercises the REAL cAtmosphereModel::MultiLayerRadiation friend class on a
// single representative (US-standard) atmospheric column, WITHOUT wiring any
// call site into RunTimeSlice / run_3D_loop. Purpose: confirm the scheme runs
// on the model's actual 41-layer stretched grid with the p_stat units fix,
// produces no NaN, and yields physically sane epsilon / radiation / T + OLR.
//
// Build (from repo root, links the existing libatom.a — the friend line in
// cAtmosphereModel.h does not change the class layout so no rebuild is needed):
//   g++ -std=c++17 -march=native -fopenmp -Ilib -Iatmosphere -Ihydrosphere \
//       -Itinyxml2 test/rad_selftest.cpp -L. -latom -o rad_selftest
//   ./rad_selftest
// (compiled WITHOUT -ffast-math so the isnan() output check is reliable.)
// ============================================================================
#include "cAtmosphereModel.h"
#include "MultiLayerRadiation.h"
#include <cmath>
#include <cstdio>
#include <vector>

class RadiationSelfTest {
public:
    void run() {
        cAtmosphereModel m;                         // ctor: SetDefaultConfig + rad/the/phi alloc
        const int im = m.im, jm = m.jm, km = m.km;

        // radial coordinate rad.z[i] = 1 + i*dr (dr = 0.025), then layer heights
        for (int i = 0; i < im; i++) m.rad.z[i] = 1.0 + i * 0.025;
        m.init_layer_heights();                     // private -> reachable via friend

        // allocate the fields the scheme reads / writes (ctor leaves these unallocated)
        m.t.initArray(im, jm, km, 1.0);
        m.c.initArray(im, jm, km, 0.0);
        m.co2.initArray(im, jm, km, 1.0);           // ratio 1.0 -> co2_0 (=380 ppm)
        m.p_stat.initArray(im, jm, km, 0.0);
        m.h.initArray(im, jm, km, 0.0);             // 0 everywhere -> ocean surface at i=0
        m.radiation.initArray(im, jm, km, 0.0);
        m.epsilon.initArray(im, jm, km, 0.0);
        m.albedo.initArray_2D(jm, km, 0.0);
        m.epsilon_2D.initArray_2D(jm, km, 0.0);

        // build a US-standard-ish column (T, p, RH=0.5), replicated over all (j,k)
        for (int i = 0; i < im; i++) {
            double zkm = m.get_layer_height(i) / 1000.0;
            double T   = (zkm < 11.0) ? 288.15 - 6.5 * zkm : 216.65;                 // K
            double p   = (zkm < 11.0) ? 1013.25 * pow(1.0 - 0.0065*zkm*1000.0/288.15, 5.256)
                                      : 226.32 * exp(-(zkm - 11.0)*1000.0/6341.6);    // hPa
            double esat = 6.1078 * exp(17.2694*(T-273.15)/(T-35.86));                 // hPa
            double cc   = m.ep * (0.5*esat) / p;                                      // kg/kg (RH=0.5)
            for (int j = 0; j < jm; j++)
                for (int k = 0; k < km; k++) {
                    m.t.x[i][j][k]      = T / m.t_0;   // nondim temperature
                    m.p_stat.x[i][j][k] = p;
                    m.c.x[i][j][k]      = cc;
                    m.co2.x[i][j][k]    = 1.0;
                }
        }

        const int j0 = jm / 2, k0 = km / 2;           // equatorial-ish probe column
        std::vector<double> Tin(im);
        for (int i = 0; i < im; i++) Tin[i] = m.t.x[i][j0][k0] * m.t_0;

        // ---------------- RUN THE REAL SCHEME ----------------
        MultiLayerRadiation(m).run();

        // ---------------- report probe column ----------------
        printf("\n  i    z[m]   Tin[K]  Tout[K]   eps      rad[W/m2]\n");
        bool nan = false;
        for (int i = 0; i < im; i++) {
            double Tout = m.t.x[i][j0][k0] * m.t_0;
            double eps  = m.epsilon.x[i][j0][k0];
            double rad  = m.radiation.x[i][j0][k0];
            if (std::isnan(Tout) || std::isnan(eps) || std::isnan(rad) ||
                std::isinf(Tout) || std::isinf(eps) || std::isinf(rad)) nan = true;
            printf("%3d %7.0f  %6.1f  %7.2f  %.4f  %9.2f\n",
                   i, m.get_layer_height(i), Tin[i], Tout, eps, rad);
        }
        // OLR: grey Schwarzschild upward integration on the scheme's eps & Tout
        double F = m.sigma * pow(m.t.x[0][j0][k0]*m.t_0, 4);
        for (int i = 1; i < im; i++) {
            double eps = m.epsilon.x[i][j0][k0];
            F = F*(1.0-eps) + eps * m.sigma * pow(m.t.x[i][j0][k0]*m.t_0, 4);
        }
        printf("\n  OLR (grey, scheme eps & Tout)   = %.2f W/m2   [Earth ~240]\n", F);
        printf("  albedo=%.3f  SW_in=%.2f W/m2  eps_2D(surf)=%.4f\n",
               m.albedo.y[j0][k0], m.short_wave_radiation[j0], m.epsilon_2D.y[j0][k0]);
        printf("  surface Tin=%.1f K  ->  Tout=%.1f K   (dT=%.1f K)\n",
               Tin[0], m.t.x[0][j0][k0]*m.t_0, m.t.x[0][j0][k0]*m.t_0 - Tin[0]);
        printf("  NaN/Inf in outputs: %s\n", nan ? "*** YES ***" : "no");
    }
};

int main() { RadiationSelfTest().run(); return 0; }
