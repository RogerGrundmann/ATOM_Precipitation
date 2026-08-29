// ============================================================================
// rad_iterate — is MultiLayerRadiation an EQUILIBRIUM solver, and where does it
// converge?
//
// rad_selftest calls the scheme ONCE on a US-standard column and reports what
// comes back. That answers "does it run", not "what is its equilibrium". The
// tree's open question (CLAUDE.md, the clear-sky OLR item) is that the scheme
// returns 3.97 K/km where Earth's radiative-convective troposphere is 6.5 and a
// grey RADIATIVE equilibrium should be steeper still -- and the named suspect is
// the tridiagonal solve itself rather than its optical-depth inputs.
//
// This driver applies the scheme REPEATEDLY to its own output. Three outcomes,
// all informative:
//   - it converges to a fixed point  -> that fixed point IS the scheme's
//     radiative equilibrium, and can be compared with the analytic grey answer;
//   - it walks without converging    -> it is not an equilibrium solver at all;
//   - it diverges                    -> the balance it solves is unstable.
// Alongside, the analytic grey profile sigma*T^4 = (F/2)(1 + 3*tau/2) is
// evaluated on the scheme's OWN tau_above, so the comparison uses the model's
// optical depth and not a textbook one.
//
// Build: g++ -std=c++17 -march=native -fopenmp -Ilib -Iatmosphere -Ihydrosphere \
//            -Itinyxml2 test/rad_iterate.cpp -L. -latom -o rad_iterate
// ============================================================================
#include "cAtmosphereModel.h"
#include "MultiLayerRadiation.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

// Reuses the RadiationSelfTest friend name declared in cAtmosphereModel.h, so this
// driver reaches the same private members the existing harness does.
class RadiationSelfTest {
public:
  void run(int npass) {
    cAtmosphereModel m;
    const int im = m.im, jm = m.jm, km = m.km;
    m.initGridCoordinates();
    m.init_layer_heights();

    m.t.initArray(im, jm, km, 1.0);
    m.c.initArray(im, jm, km, 0.0);
    m.co2.initArray(im, jm, km, 380.0);
    m.p_stat.initArray(im, jm, km, 0.0);
    m.h.initArray(im, jm, km, 0.0);
    m.radiation.initArray(im, jm, km, 0.0);
    m.epsilon.initArray(im, jm, km, 0.0);
    m.cloud.initArray(im, jm, km, 0.0);
    m.ice.initArray(im, jm, km, 0.0);
    m.albedo.initArray_2D(jm, km, 0.0);
    m.epsilon_2D.initArray_2D(jm, km, 0.0);
    m.tau_above.initArray(im, jm, km, 0.0);
    m.tau_layer.initArray(im, jm, km, 0.0);
    m.i_topography.assign(jm, std::vector<int>(km, 0));
    m.short_wave_radiation.assign(jm, 0.0);

    const int j0 = jm / 2, k0 = km / 2;

    // US-standard column, RH = 0.5, clear sky. Water vapour is held FIXED at the
    // initial profile: the scheme reads c but never writes it, so holding it is
    // what the scheme itself does, and it keeps tau from drifting under us.
    std::vector<double> Tin(im);
    for (int i = 0; i < im; i++) {
        const double zkm = m.get_layer_height(i) / 1000.0;
        const double T   = (zkm < 11.0) ? 288.15 - 6.5 * zkm : 216.65;
        const double p   = (zkm < 11.0) ? 1013.25 * pow(1.0 - 0.0065*zkm*1000.0/288.15, 5.256)
                                        : 226.32 * exp(-(zkm - 11.0)*1000.0/6341.6);
        const double esat = 6.1078 * exp(17.2694*(T-273.15)/(T-35.86));
        const double cc   = m.ep * (0.5*esat) / p;
        Tin[i] = T;
        for (int j = 0; j < jm; j++) for (int k = 0; k < km; k++) {
            m.t.x[i][j][k]      = T / m.t_0;
            m.p_stat.x[i][j][k] = p;
            m.c.x[i][j][k]      = cc;
            m.co2.x[i][j][k]    = 380.0;
        }
    }

    auto lapse = [&](int ia, int ib) {                     // K/km between two levels
        const double dz = (m.get_layer_height(ib) - m.get_layer_height(ia)) / 1000.0;
        return (m.t.x[ia][j0][k0] - m.t.x[ib][j0][k0]) * m.t_0 / dz;
    };
    auto olr = [&]() {
        double F = m.sigma * pow(m.t.x[0][j0][k0]*m.t_0, 4);
        for (int i = 1; i < im; i++) {
            const double e = m.epsilon.x[i][j0][k0];
            F = F*(1.0-e) + e*m.sigma*pow(m.t.x[i][j0][k0]*m.t_0, 4);
        }
        return F;
    };

    printf("\n  pass |  T(0 m)   T(2163)   T(6088)   T(9923)  T(16023) | lapse 0-9.9km |    OLR   | max|dT|\n");
    printf("  -----+-------------------------------------------------+---------------+----------+---------\n");
    printf("   in  | %8.2f %8.2f %8.2f %8.2f %8.2f |    %6.3f     |          |\n",
           Tin[0], Tin[20], Tin[30], Tin[35], Tin[40], (Tin[0]-Tin[35])/9.923);

    std::vector<double> prev(im);
    for (int pass = 1; pass <= npass; pass++) {
        for (int i = 0; i < im; i++) prev[i] = m.t.x[i][j0][k0] * m.t_0;
        MultiLayerRadiation(m).run();
        double dmax = 0.0;
        for (int i = 0; i < im; i++)
            dmax = std::max(dmax, std::fabs(m.t.x[i][j0][k0]*m.t_0 - prev[i]));
        if (pass <= 5 || pass % 5 == 0 || pass == npass)
            printf("  %4d | %8.2f %8.2f %8.2f %8.2f %8.2f |    %6.3f     | %8.2f | %7.4f\n",
                   pass, m.t.x[0][j0][k0]*m.t_0, m.t.x[20][j0][k0]*m.t_0, m.t.x[30][j0][k0]*m.t_0,
                   m.t.x[35][j0][k0]*m.t_0, m.t.x[40][j0][k0]*m.t_0, lapse(0, 35), olr(), dmax);
    }

    // ---- the profile against the ANALYTIC grey answer, on the scheme's OWN tau_above.
    //      NOT normalised or fitted: F is the scheme's own OLR, so the lid value is a
    //      PREDICTION (grey gives (F/2sigma)^0.25) and its agreement is a real test.
    const double F_olr = olr();
    printf("\n  grey radiative equilibrium on the scheme's own tau:  sigma*T^4 = (F/2)(1 + 3*tau/2),  F = %.2f W/m2\n", F_olr);
    printf("  lvl      z[m]   tau_above    T model    T grey     model-grey\n");
    for (int i = im-1; i >= 0; i -= (i > 24 ? 5 : 2)) {
        const double tau  = m.tau_above.x[i][j0][k0];
        const double B    = 0.5 * F_olr * (1.0 + 1.5 * tau);
        const double Tg   = pow(B / m.sigma, 0.25);
        const double Tm   = m.t.x[i][j0][k0] * m.t_0;
        printf("  %3d  %8.0f  %9.4f  %9.2f %9.2f     %+8.2f\n", i, m.get_layer_height(i), tau, Tm, Tg, Tm - Tg);
    }
  }
};

int main(int argc, char** argv) {
    RadiationSelfTest().run((argc > 1) ? atoi(argv[1]) : 40);
    return 0;
}
