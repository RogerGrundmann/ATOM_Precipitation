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
    // Shortwave override, for testing whether the scheme's SW forcing is the reason a CORRECT
    // radiative-equilibrium solve lands in a snowball. Earth's annual-mean TOA insolation is
    // ~416 W/m2 at the equator and ~173 at the pole (global mean 341); the shipped constants are
    // 163.3 and 100.0.
    if (getenv("RAD_SW_EQ")) m.rad_equator_short = atof(getenv("RAD_SW_EQ"));
    if (getenv("RAD_SW_PO")) m.rad_pole_short    = atof(getenv("RAD_SW_PO"));
    printf("\n  shortwave: rad_equator_short = %.1f  rad_pole_short = %.1f W/m2\n",
           m.rad_equator_short, m.rad_pole_short);

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

    // ---- DOES THE OUTPUT SATISFY THE SCHEME'S OWN ENERGY BALANCE? ------------------------
    //
    // The grey comparison below uses an EXTERNAL yardstick, so on its own it cannot separate
    // "the solve is wrong" from "grey is the wrong reference for this scheme". This test uses
    // no reference at all. For a layer of emissivity eps_i between an upward flux U_i entering
    // from below and a downward flux D_i entering from above, radiative equilibrium is
    //     absorbed = emitted   ->   eps_i*(U_i + D_i) = 2*eps_i*B_i   ->   B_i = (U_i + D_i)/2,
    // and eps_i cancels, so the balance is independent of the optical depth. U and D are swept
    // from the scheme's OWN epsilon and OWN output radiation, with the surface as a blackbody
    // below and no layer above the lid. The residual is reported as the temperature the layer
    // WOULD have if it were in balance, minus the temperature the scheme returned.
    {
        std::vector<double> B(im), U(im), D(im), E(im);
        for (int i = 0; i < im; i++) {
            B[i] = m.radiation.x[i][j0][k0];
            E[i] = m.epsilon.x[i][j0][k0];
        }
        // U[i] = upward flux ENTERING layer i from below. Layer 0 is the surface: it emits as a
        // blackbody and does not attenuate its own emission, so U[1] = B[0] exactly.
        U[0] = 0.0;
        U[1] = B[0];
        for (int i = 2; i < im; i++) U[i] = U[i-1]*(1.0 - E[i-1]) + E[i-1]*B[i-1];
        D[im-1] = 0.0;                                      // nothing above the lid
        for (int i = im-2; i >= 0; i--) D[i] = D[i+1]*(1.0 - E[i+1]) + E[i+1]*B[i+1];
        printf("\n  the scheme's OWN balance B_i = (U_i + D_i)/2, on its own epsilon and output:\n");
        printf("  lvl      z[m]      B out    B balance   resid[W/m2]   T out    T balance    dT\n");
        // i = 0 is the SURFACE, which also absorbs shortwave, so absorbed = emitted does not
        // apply to it; the balance is a statement about air layers only.
        for (int i = im-1; i >= 1; i -= (i > 24 ? 5 : 2)) {
            const double Bb = 0.5*(U[i] + D[i]);
            const double To = pow(std::max(1.0, B[i]) / m.sigma, 0.25);
            const double Tb = pow(std::max(1.0, Bb  ) / m.sigma, 0.25);
            printf("  %3d  %8.0f  %9.2f  %9.2f   %+9.2f   %8.2f %8.2f  %+8.2f\n",
                   i, m.get_layer_height(i), B[i], Bb, Bb - B[i], To, Tb, Tb - To);
        }
    }

    // ---- THE REFERENCE SOLUTION: the same balance, solved to convergence -------------------
    //
    // The residual above says the scheme's output does not satisfy absorbed = emitted. It does
    // not say what profile would. So solve that balance directly, on the scheme's OWN epsilon
    // and OWN shortwave, by Jacobi iteration:
    //     air     B_i = (U_i + D_i)/2        (absorbed = emitted, eps cancels)
    //     surface B_0 = SW_abs + D_0         (absorbs shortwave plus the downwelling)
    // No diffusivity factor is assumed and no analytic form is imposed -- this is the fixed
    // point of the scheme's own physics. TOA closure OLR = SW_abs is the validation: it is not
    // enforced anywhere in the iteration, so if it comes out satisfied the solve is right.
    {
        std::vector<double> B(im), U(im), D(im), E(im);
        for (int i = 0; i < im; i++) { B[i] = m.radiation.x[i][j0][k0]; E[i] = m.epsilon.x[i][j0][k0]; }
        const double SW_abs = (1.0 - m.albedo.y[j0][k0]) * m.short_wave_radiation[j0];
        for (int it = 0; it < 200000; it++) {
            U[0] = 0.0; U[1] = B[0];
            for (int i = 2; i < im; i++) U[i] = U[i-1]*(1.0 - E[i-1]) + E[i-1]*B[i-1];
            D[im-1] = 0.0;
            for (int i = im-2; i >= 0; i--) D[i] = D[i+1]*(1.0 - E[i+1]) + E[i+1]*B[i+1];
            for (int i = 1; i < im; i++) B[i] = 0.5*(U[i] + D[i]);
            B[0] = SW_abs + D[0];
        }
        double F = B[0];
        for (int i = 1; i < im; i++) F = F*(1.0 - E[i]) + E[i]*B[i];
        printf("\n  REFERENCE: the same balance solved to convergence on the scheme's own epsilon\n");
        printf("  SW absorbed = %.2f W/m2, converged OLR = %.2f W/m2  (TOA closure error %+.3f, not enforced)\n",
               SW_abs, F, F - SW_abs);
        printf("  lvl      z[m]    T scheme   T reference     diff     lapse ref [K/km]\n");
        for (int i = im-1; i >= 0; i -= (i > 24 ? 5 : 2)) {
            const double Ts = m.t.x[i][j0][k0] * m.t_0;
            const double Tr = pow(std::max(1.0, B[i]) / m.sigma, 0.25);
            double lap = 0.0;
            if (i > 0) {
                const double Trm = pow(std::max(1.0, B[i-1]) / m.sigma, 0.25);
                lap = (Trm - Tr) / ((m.get_layer_height(i) - m.get_layer_height(i-1)) / 1000.0);
            }
            printf("  %3d  %8.0f  %9.2f  %11.2f  %+8.2f   %8.3f\n", i, m.get_layer_height(i), Ts, Tr, Tr - Ts, lap);
        }
        // The Jacobi fixed point above has a CLOSED FORM, and checking it is what turns an
        // expensive reference into a usable solver.
        //
        // Care with the two flux arrays: U[i] crosses the interface BELOW layer i and D[i] the
        // interface ABOVE it, so U[i] - D[i] is NOT a net flux at one level. The net at the
        // upper interface of layer i is U[i+1] - D[i], and THAT is what radiative equilibrium
        // holds constant at F. With B_i = (U[i] + D[i])/2 and the transfer step
        // U[i+1] = U[i](1-eps_i) + eps_i*B_i, eliminating B and D gives
        //     U[i+1] = U[i] - eps_i*F/(2 - eps_i),
        //     B_i    = U[i] - F/2 - eps_i*F/(2*(2 - eps_i)),
        // integrated downward from U above the lid = F, with the surface at B_0 = U[1] because a
        // blackbody surface does not attenuate its own emission. The eps/(2-eps) rather than
        // eps/2 is the top layer's opacity to its own emission; assuming eps/2 costs 2.23 K.
        {
            std::vector<double> Uc(im+1, 0.0), Bc(im, 0.0);
            const double F = SW_abs;
            Uc[im] = F;
            for (int i = im-1; i >= 1; i--) Uc[i] = Uc[i+1] + E[i]*F/(2.0 - E[i]);
            for (int i = 1; i < im; i++)    Bc[i] = Uc[i] - 0.5*F - E[i]*F/(2.0*(2.0 - E[i]));
            Bc[0] = Uc[1];
            double worst = 0.0;
            for (int i = 0; i < im; i++) {
                const double Tj = pow(std::max(1.0, B [i]) / m.sigma, 0.25);
                const double Tc = pow(std::max(1.0, Bc[i]) / m.sigma, 0.25);
                worst = std::max(worst, std::fabs(Tj - Tc));
            }
            printf("  closed form vs the 200 000-sweep Jacobi: max |dT| = %.6f K\n", worst);
        }
        const double T0 = pow(std::max(1.0,B[0])/m.sigma, 0.25);
        const double T20 = pow(std::max(1.0,B[20])/m.sigma, 0.25);
        const double T35 = pow(std::max(1.0,B[35])/m.sigma, 0.25);
        printf("  reference lapse 0 -> 2163 m = %.3f K/km   (dry adiabat 9.8);  0 -> 9923 m = %.3f K/km\n",
               (T0-T20)/2.163, (T0-T35)/9.923);
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
