/*
 * Ocean General Circulation Modell(OGCM) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
*/

#include "cHydrosphereModel.h"

using namespace std;

namespace{
    string heading_1 = " printout of maximum and minimum values of properties at their locations: latitude, longitude, level";
    string heading_2 = " results based on three dimensional considerations of the problem";
    string level = "m";
    struct HemisphereCoords{
        double lat, lon;
        string east_or_west, north_or_south;
    };
    HemisphereCoords convert_coords(double lon, double lat){
        HemisphereCoords ret;
        if(lat > 90){
            ret.lat = lat - 90;
            ret.north_or_south = "°S";
        }else{
            ret.lat = 90 - lat;
            ret.north_or_south = "°N";
        }
        if(lon > 180){
            ret.lon = 360 - lon;
            ret.east_or_west = "°W";
        }else{
            ret.lon = lon;
            ret.east_or_west = "°E";
        }
        return ret;
    }
}
/*
*
*/
// Deterministic tie-break for the parallel min/max reporters below.
//
// Both searchMinMax_2D and searchMinMax_3D reduce over threads with `#pragma omp critical`
// and a plain `>`, so whichever thread reached the section FIRST won an equal extremum: the
// reported LOCATION of a tied max/min depended on thread arrival order while its value did
// not. Measured 2026-08-26 over two 4-iteration runs of the same binary at the same 4
// threads: `max temperature = 0.000000` reported at 0N 80W in one run and 45N 91E in the
// other, and likewise `min water density = 997.000000` (the density floor) and
// `min salinity = 0.000000`. Ties are not rare -- they are what a CAP or a floor produces,
// so the fields most likely to tie are exactly the ones worth watching.
//
// Preferring the lexicographically smallest (i,j,k) -- in the per-thread scan AND in the
// merge -- makes the reported cell a property of the field rather than of the schedule, at
// any thread count. The atmosphere's searchMinMax_* is a serial loop and never had this.
static inline bool lex_less(int i, int j, int k, int i0, int j0, int k0){
    if(i != i0) return i < i0;
    if(j != j0) return j < j0;
    return k < k0;
}

void cHydrosphereModel::searchMinMax_3D(string name_maxValue, string name_minValue, 
      string name_unitValue, Array &value_D, double coeff, 
      std::function< double(double) > lambda, bool print_heading){                                                                                                                   
   
      double minValue = value_D.x[0][0][0];                                                                                                                                          
      double maxValue = value_D.x[0][0][0];
      int imax = 0, jmax = 0, kmax = 0;
      int imin = 0, jmin = 0, kmin = 0;

      #pragma omp parallel
      {
          double localMin = value_D.x[0][0][0];
          double localMax = value_D.x[0][0][0];
          int li_max = 0, lj_max = 0, lk_max = 0;
          int li_min = 0, lj_min = 0, lk_min = 0;

          #pragma omp for collapse(2) nowait
          for(int j = 0; j < jm; j++){
              for(int k = 0; k < km; k++){
                  for(int i = 0; i < im; i++){
                      double val = value_D.x[i][j][k];
                      // Ties are broken toward the lexicographically smallest (i,j,k), here
                      // and in the merge below. See the note above this function.
                      if(val > localMax || (val == localMax
                              && lex_less(i, j, k, li_max, lj_max, lk_max))){
                          localMax = val;
                          li_max = i; lj_max = j; lk_max = k;
                      }
                      if(val < localMin || (val == localMin
                              && lex_less(i, j, k, li_min, lj_min, lk_min))){
                          localMin = val;
                          li_min = i; lj_min = j; lk_min = k;
                      }
                  }
              }
          }

          #pragma omp critical
          {
              if(localMax > maxValue || (localMax == maxValue
                      && lex_less(li_max, lj_max, lk_max, imax, jmax, kmax))){
                  maxValue = localMax;
                  imax = li_max; jmax = lj_max; kmax = lk_max;
              }
              if(localMin < minValue || (localMin == minValue
                      && lex_less(li_min, lj_min, lk_min, imin, jmin, kmin))){
                  minValue = localMin;
                  imin = li_min; jmin = lj_min; kmin = lk_min;
              }
          }
      }

      int imax_level = imax * (int)L_hyd/(im-1) - int(L_hyd);
      int imin_level = imin * (int)L_hyd/(im-1) - int(L_hyd);
      HemisphereCoords coords = convert_coords(kmax, jmax);
      int jmax_deg = coords.lat;
      string deg_lat_max = coords.north_or_south;
      int kmax_deg = coords.lon;
      string deg_lon_max = coords.east_or_west;
      coords = convert_coords(kmin, jmin);
      int jmin_deg = coords.lat;
      string deg_lat_min= coords.north_or_south;
      int kmin_deg = coords.lon;
      string deg_lon_min = coords.east_or_west;
      cout.precision(6);
      if(print_heading){
          cout << endl << heading_1 << endl << heading_2 << endl << endl;
      }
      maxValue = lambda(maxValue * coeff);
      minValue = lambda(minValue * coeff);
      cout << setiosflags(ios::left) << setw(26) << setfill('.') << name_maxValue << " = " <<
          resetiosflags(ios::left) << setw(12) << fixed << setfill(' ') << maxValue << setw(6) <<
          name_unitValue << setw(5) << jmax_deg << setw(3) << deg_lat_max << setw(4) << kmax_deg <<
          setw(3) << deg_lon_max << setw(6) << imax_level << setw(2) << level << "          " <<
          setiosflags(ios::left) << setw(26) << setfill('.') << name_minValue << " = "<<
          resetiosflags(ios::left) << setw(12) << fixed << setfill(' ') << minValue << setw(6) <<
          name_unitValue << setw(5)  << jmin_deg << setw(3) << deg_lat_min << setw(4) << kmin_deg <<
          setw(3) << deg_lon_min  << setw(6) << imin_level << setw(2) << level << endl;
  }
/*
*
*/
void cHydrosphereModel::searchMinMax_2D(string name_maxValue, string name_minValue, 
      string name_unitValue, Array_2D &value, double coeff){
                                                                                                                                                                                     
      double minValue = value.y[0][0];
      double maxValue = value.y[0][0];                                                                                                                                               
      int jmax = 0, kmax = 0;
      int jmin = 0, kmin = 0;

      #pragma omp parallel
      {
          double localMin = value.y[0][0];
          double localMax = value.y[0][0];
          int lj_max = 0, lk_max = 0;
          int lj_min = 0, lk_min = 0;

          #pragma omp for collapse(2) nowait
          for(int j = 1; j < jm-1; j++){
              for(int k = 1; k < km-1; k++){
                  double val = value.y[j][k];
                  if(val > localMax || (val == localMax
                          && lex_less(0, j, k, 0, lj_max, lk_max))){
                      localMax = val;
                      lj_max = j; lk_max = k;
                  }
                  if(val < localMin || (val == localMin
                          && lex_less(0, j, k, 0, lj_min, lk_min))){
                      localMin = val;
                      lj_min = j; lk_min = k;
                  }
              }
          }

          #pragma omp critical
          {
              if(localMax > maxValue || (localMax == maxValue
                      && lex_less(0, lj_max, lk_max, 0, jmax, kmax))){
                  maxValue = localMax;
                  jmax = lj_max; kmax = lk_max;
              }
              if(localMin < minValue || (localMin == minValue
                      && lex_less(0, lj_min, lk_min, 0, jmin, kmin))){
                  minValue = localMin;
                  jmin = lj_min; kmin = lk_min;
              }
          }
      }

      int imax_level = 0;
      int imin_level = 0;
      HemisphereCoords coords = convert_coords(kmax, jmax);
      int jmax_deg = coords.lat;
      string deg_lat_max = coords.north_or_south;
      int kmax_deg = coords.lon;
      string deg_lon_max = coords.east_or_west;
      coords = convert_coords(kmin, jmin);
      int jmin_deg = coords.lat;
      string deg_lat_min= coords.north_or_south;
      int kmin_deg = coords.lon;
      string deg_lon_min = coords.east_or_west;
      cout.precision(6);
      maxValue = maxValue * coeff;
      minValue = minValue * coeff;
      cout << setiosflags(ios::left) << setw(26) << setfill('.') << name_maxValue << " = " <<
          resetiosflags(ios::left) << setw(12) << fixed << setfill(' ') << maxValue << setw(6) <<
          name_unitValue << setw(5) << jmax_deg << setw(3) << deg_lat_max << setw(4) << kmax_deg <<
          setw(3) << deg_lon_max << setw(6) << imax_level << setw(2) << level << "          " <<
          setiosflags(ios::left) << setw(26) << setfill('.') << name_minValue << " = "<<
          resetiosflags(ios::left) << setw(12) << fixed << setfill(' ') << minValue << setw(6) <<
          name_unitValue << setw(5)  << jmin_deg << setw(3) << deg_lat_min << setw(4) << kmin_deg <<
          setw(3) << deg_lon_min  << setw(6) << imin_level << setw(2) << level << endl;
  }
/*
*
*/


// ============================================================================
// Zonal-mean zonal-velocity (w) momentum budget.
//
// Diagnoses why the Southern-Hemisphere surface circulation spins down and
// accelerates westward over a spin-up (the prescribed eastward gyre return /
// West Wind Drift reverses; see project_hydro_ekman_sh_gyre). The zonal-mean
// zonal-momentum equation has (to rounding) no net pressure-gradient term —
// -dp/dphi integrates out over a periodic longitude circle — so a sustained
// zonal-mean acceleration must come from Coriolis (-2*Omega*(cos*v+sin*u)),
// advection, or diffusion. The per-term split (captured in wbud_* during the
// RK4) is compared per latitude/depth, alongside the net per-iteration change
// of wbar from the wbar_before snapshot (which also includes the pressure
// projection and the polar zonal filter applied after the RHS).
//
// Units: cm/s for wbar; cm/s per iteration for the term contributions
// (rhs_w term * dt * u_0, the Euler-step increment, since the RK4 advances
// w_nd by rhs_w * dt). A term NEGATIVE where wbar<0 is ERODING / reversing the
// (eastward) flow toward westward at that (lat, depth).
// ----------------------------------------------------------------------------
void cHydrosphereModel::zonal_mean_w(std::vector<std::vector<double> >& wbar){
    for(int i = 0; i < im; i++)
        for(int j = 0; j < jm; j++){
            double sum = 0.0; int n = 0;
            for(int k = 0; k < km; k++){
                if(!AtomUtils::is_water(h, i, j, k)) continue;
                double a = w.x[i][j][k];
                if(!AtomUtils::is_finite_safe(a)) continue;
                sum += a; n++;
            }
            wbar[i][j] = (n > 0) ? sum / n * u_0 * 100.0 : 0.0;   // cm/s, East+
        }
}

void cHydrosphereModel::write_w_momentum_budget(int iter){
    std::vector<std::vector<double> > wbar(im, std::vector<double>(jm, 0.0));
    zonal_mean_w(wbar);   // post-iteration zonal-mean w [cm/s]

    auto lat_of = [&](int j){ return 90.0 - (double)j * 180.0 / (double)(jm - 1); };
    const double term_scale = dt * u_0 * 100.0;   // nondim rhs_w term -> cm/s per iter

    auto zmean = [&](const Array& A, std::vector<std::vector<double> >& out){
        for(int i = 0; i < im; i++)
            for(int j = 0; j < jm; j++){
                double s = 0.0; int n = 0;
                for(int k = 0; k < km; k++){
                    if(!AtomUtils::is_water(h, i, j, k)) continue;
                    double a = A.x[i][j][k];
                    if(!AtomUtils::is_finite_safe(a)) continue;
                    s += a; n++;
                }
                out[i][j] = (n > 0) ? s / n * term_scale : 0.0;
            }
    };
    std::vector<std::vector<double> > t_pgf(im, std::vector<double>(jm, 0.0)),
        t_cor = t_pgf, t_adv = t_pgf, t_diff = t_pgf;
    zmean(wbud_pgf, t_pgf);  zmean(wbud_cor, t_cor);
    zmean(wbud_adv, t_adv);  zmean(wbud_diff, t_diff);

    const bool have_before = ((int)wbar_before.size() == im);

    std::ostringstream fname;
    fname << output_path << "/w_momentum_budget_" << iter << ".csv";
    std::ofstream f(fname.str().c_str());
    if(f.is_open()){
        f << "lat_deg,level_i,depth_m,wbar_cms,dwbar_net_cms,"
          << "pgf,coriolis,advection,diffusion,dyn_sum_cms\n";
        for(int j = 0; j < jm; j++){
            const double lat = lat_of(j);
            for(int i = 0; i < im; i++){
                const double depth   = (rad.z[im-1] - rad.z[i]) * L_hyd;
                const double dyn_sum = t_pgf[i][j] + t_cor[i][j] + t_adv[i][j] + t_diff[i][j];
                const double dnet    = have_before ? (wbar[i][j] - wbar_before[i][j]) : 0.0;
                f << lat << "," << i << "," << depth << "," << wbar[i][j] << ","
                  << dnet << "," << t_pgf[i][j] << "," << t_cor[i][j] << ","
                  << t_adv[i][j] << "," << t_diff[i][j] << "," << dyn_sum << "\n";
            }
        }
        f.close();
        cout << "      OGCM: wrote " << fname.str() << endl;
    }

    // Console summary at the SH spin-down band (~40 S, near-surface interior).
    // i=im-1 is the sea-surface BC row (no RHS terms), so probe the layer below.
    const int jp = (int)round((90.0 + 40.0) * (jm - 1) / 180.0);   // ~40 S
    const int ip = im - 2;                                          // near-surface interior
    const double dyn_sum = t_pgf[ip][jp] + t_cor[ip][jp] + t_adv[ip][jp] + t_diff[ip][jp];
    cout << "      [wbudget] iter=" << iter << fixed << setprecision(4)
         << "  @lat=" << setprecision(0) << lat_of(jp) << " near-surf" << setprecision(4)
         << "  wbar=" << wbar[ip][jp] << "cm/s"
         << "  pgf=" << t_pgf[ip][jp] << " cor=" << t_cor[ip][jp]
         << " adv=" << t_adv[ip][jp] << " diff=" << t_diff[ip][jp]
         << "  dyn_sum=" << dyn_sum
         << "  dnet=" << (have_before ? wbar[ip][jp] - wbar_before[ip][jp] : 0.0) << endl;
}

// ----------------------------------------------------------------------------
// Meridional (v) momentum budget — LOCAL, per longitude, depth-integrated.
//
// The subtropical gyres close through an EASTERN boundary current instead of a
// western-boundary current (project_hydro_eastern_boundary_current): the
// barotropic streamfunction has essentially all its gradient on the eastern
// coast (N Atlantic 40N: |dPsi| 469 east vs 0.45 west). A ZONAL-MEAN budget
// hides this (a gyre is a longitudinal structure; the meridional PGF integrates
// out zonally). So this dumps the per-longitude, DEPTH-INTEGRATED v-momentum
// term split at a few gyre latitudes: for each ocean column (j,k) it sums each
// rhs_v contribution over the water cells. The dominant term at the eastern
// boundary vs the (empty) western boundary reveals what builds the eastern
// current — geostrophic PGF, Coriolis, friction (diffusion) or advection.
//
// Cols: lat_deg,lon_deg,v_baro_cms,pgf,coriolis,advection,diffusion,wind,dyn_sum
// (terms = depth-summed rhs_v contributions * dt*u_0*100 = cm/s per iter;
//  v_baro = depth-summed v * u_0 * 100, a barotropic-transport proxy in cm/s;
//  lon = -180 + k*360/(km-1) to match psi_diag.py / the 0Ma_smooth.xyz topo).
// ----------------------------------------------------------------------------
void cHydrosphereModel::write_v_momentum_budget(int iter){
    const double term_scale = dt * u_0 * 100.0;     // rhs_v term -> cm/s per iter
    const double v_scale    = u_0 * 100.0;          // v -> cm/s
    auto lat_of = [&](int j){ return 90.0 - (double)j * 180.0 / (double)(jm - 1); };
    auto lon_of = [&](int k){ return -180.0 + (double)k * 360.0 / (double)(km - 1); };

    // gyre latitudes to probe (deg): subtropical gyre cores both hemispheres.
    const double probe_lat[] = {40.0, 30.0, 20.0, -20.0, -30.0, -40.0};

    std::ostringstream fname;
    fname << output_path << "/v_momentum_budget_" << iter << ".csv";
    std::ofstream f(fname.str().c_str());
    if(!f.is_open()) return;
    f << "lat_deg,lon_deg,v_baro_cms,pgf,coriolis,advection,diffusion,wind,dyn_sum\n";

    for(double plat : probe_lat){
        const int j = (int)round((90.0 - plat) * (jm - 1) / 180.0);
        if(j < 0 || j >= jm) continue;
        for(int k = 0; k < km; k++){
            double vb = 0.0, pgf = 0.0, cor = 0.0, adv = 0.0, dif = 0.0, wnd = 0.0;
            int nwater = 0;
            for(int i = 0; i < im; i++){
                if(!AtomUtils::is_water(h, i, j, k)) continue;
                nwater++;
                vb  += v.x[i][j][k];
                pgf += vbud_pgf.x[i][j][k];
                cor += vbud_cor.x[i][j][k];
                adv += vbud_adv.x[i][j][k];
                dif += vbud_diff.x[i][j][k];
                wnd += vbud_wind.x[i][j][k];
            }
            if(nwater == 0) continue;               // land column
            const double sum = (pgf + cor + adv + dif + wnd) * term_scale;
            f << lat_of(j) << "," << lon_of(k) << "," << vb * v_scale << ","
              << pgf * term_scale << "," << cor * term_scale << ","
              << adv * term_scale << "," << dif * term_scale << ","
              << wnd * term_scale << "," << sum << "\n";
        }
    }
    f.close();
    cout << "      OGCM: wrote " << fname.str() << endl;
}

// ----------------------------------------------------------------------------
// Radial (u) momentum budget at the DEEP blow-up cell.
//
// The 0Ma spin-up is clean to ~300-700 iters but blows up by ~1500 as an
// isolated deep-cell radial-velocity runaway that relocates when clamped
// (N Pole -> 84N -> 51S; see project_hydro_polar_blowup). A zonal-mean budget
// smears a single-cell blow-up out, so instead this SELF-LOCATES the worst
// cell — the interior ocean cell with the largest |u| (radial velocity) — and
// dumps its full per-term radial-momentum split, captured in ubud_* during the
// RK4. The five terms sum to rhs_u; comparing them isolates the driver
// (metric-singular advection/diffusion near the pole, buoyancy, Coriolis, ...)
// wherever the blow-up sits this checkpoint.
//
// Units: cm/s per iteration for each term (rhs_u contribution * dt * u_0, the
// Euler-step increment the RK4 applies), u in cm/s. inv_sin2 = 1/sin^2(colat)
// is the polar metric-amplification factor carried by the 1/sin^2 diffusion
// and 1/sin advection terms — it diverges at the poles.
//
// A single trajectory file deep_momentum_budget.csv is appended one row per
// checkpoint so the run-up to the blow-up is visible in one place.
// ----------------------------------------------------------------------------
void cHydrosphereModel::write_deep_momentum_budget(int iter){
    // Locate the worst interior ocean cell by |u| (radial velocity). Skip the
    // seafloor (i=0) / sea-surface (i=im-1) radial BC rows and the pole
    // (j=0 / j=jm-1) latitudinal BC rows — bcTheta sets the pole rows by plain
    // copy of j=1 / j=jm-2 and the RK4 skips them, so ubud_* there is 0 and would
    // mask the real (dynamically-updated) near-pole blow-up cell at j=1.
    int bi = -1, bj = -1, bk = -1;
    double umax = -1.0;
    for(int i = 1; i <= im - 2; i++)
        for(int j = 1; j <= jm - 2; j++)
            for(int k = 0; k < km; k++){
                if(!AtomUtils::is_water(h, i, j, k)) continue;
                double a = u.x[i][j][k];
                if(!AtomUtils::is_finite_safe(a)) a = 1.0e30;   // NaN/Inf ranks worst
                if(fabs(a) > umax){ umax = fabs(a); bi = i; bj = j; bk = k; }
            }
    if(bi < 0) return;   // no interior ocean cell

    const double term_scale = dt * u_0 * 100.0;                 // rhs_u term -> cm/s per iter
    const double lat = 90.0 - (double)bj * 180.0 / (double)(jm - 1);
    const double lon = (double)bk * 360.0 / (double)(km - 1);
    const double depth = (rad.z[im-1] - rad.z[bi]) * L_hyd;     // [m] below surface
    const double colat_rad = M_PI * (double)bj / (double)(jm - 1);
    const double s = sin(colat_rad);
    const double inv_sin2 = (fabs(s) > 1.0e-12) ? 1.0 / (s * s) : 0.0;

    const double pgf  = ubud_pgf.x[bi][bj][bk]  * term_scale;
    const double adv  = ubud_adv.x[bi][bj][bk]  * term_scale;
    const double diff = ubud_diff.x[bi][bj][bk] * term_scale;
    const double buoy = ubud_buoy.x[bi][bj][bk] * term_scale;
    const double cor  = ubud_cor.x[bi][bj][bk]  * term_scale;
    const double rhs_sum = pgf + adv + diff + buoy + cor;       // = rhs_u increment [cm/s/iter]

    const double u_cms = u.x[bi][bj][bk] * u_0 * 100.0;
    const double v_cms = v.x[bi][bj][bk] * u_0 * 100.0;
    const double w_cms = w.x[bi][bj][bk] * u_0 * 100.0;

    std::ostringstream fname;
    fname << output_path << "/deep_momentum_budget.csv";
    std::ofstream f(fname.str().c_str(), std::ios::app);
    if(f.is_open()){
        if(f.tellp() == std::streampos(0))
            f << "iter,i,j,k,lat_deg,lon_deg,depth_m,inv_sin2,"
              << "u_cms,v_cms,w_cms,pgf,advection,diffusion,buoyancy,coriolis,rhs_sum\n";
        f << iter << "," << bi << "," << bj << "," << bk << ","
          << lat << "," << lon << "," << depth << "," << inv_sin2 << ","
          << u_cms << "," << v_cms << "," << w_cms << ","
          << pgf << "," << adv << "," << diff << "," << buoy << "," << cor << ","
          << rhs_sum << "\n";
        f.close();
    }

    cout << "      [ubudget] iter=" << iter << fixed << setprecision(2)
         << "  worst|u| @i=" << bi << " lat=" << setprecision(1) << lat
         << " lon=" << lon << " depth=" << setprecision(0) << depth << "m"
         << setprecision(3)
         << "  u=" << u_cms << "cm/s (v=" << v_cms << " w=" << w_cms << ")"
         << "  inv_sin2=" << setprecision(1) << inv_sin2 << setprecision(4)
         << " | pgf=" << pgf << " adv=" << adv << " diff=" << diff
         << " buoy=" << buoy << " cor=" << cor
         << " => rhs_u=" << rhs_sum << " cm/s/iter" << endl;
}

// ----------------------------------------------------------------------------
// Stage-resolved velocity-change decomposition at the deep blow-up cell.
//
// The radial-momentum budget (write_deep_momentum_budget) showed the momentum
// RHS is net-DAMPING at the near-pole runaway cell, yet u/v/w grow — so a step
// AFTER the RK4 solve injects the growth. These three helpers snapshot u,v,w at
// the worst cell after each iteration stage and attribute the per-iteration net
// Δ to RK4 / projection / BC+clamp / polar-filter, pinning the injector.
//   stage 0: entering the iteration (post previous store)   [loop top]
//   stage 1: after solveRungeKutta_Hydrosphere(_Turb)       (RK4)
//   stage 2: after PressureSolverHyd::project_velocity      (projection)
//   stage 3: after valueLimitationHyd + bc* + apply_wall_bc (BC+clamp)
//   stage 4: after polar_zonal_filter(u,v,w)                (polar filter)
// The cell is located at loop top (the persistent max-|u| interior cell); pole
// BC rows j=0/jm-1 are excluded (bcTheta sets them by copy; the RK4 skips them).
// ----------------------------------------------------------------------------
void cHydrosphereModel::locate_blowup_cell(){
    scell_i = scell_j = scell_k = -1;
    double umax = -1.0;
    for(int i = 1; i <= im - 2; i++)
        for(int j = 1; j <= jm - 2; j++)
            for(int k = 0; k < km; k++){
                if(!AtomUtils::is_water(h, i, j, k)) continue;
                double a = u.x[i][j][k];
                if(!AtomUtils::is_finite_safe(a)) a = 1.0e30;   // NaN/Inf ranks worst
                if(fabs(a) > umax){ umax = fabs(a); scell_i = i; scell_j = j; scell_k = k; }
            }
}

void cHydrosphereModel::record_stage(int s){
    if(!stage_capture || scell_i < 0) return;
    su[s] = u.x[scell_i][scell_j][scell_k];
    sv[s] = v.x[scell_i][scell_j][scell_k];
    sw[s] = w.x[scell_i][scell_j][scell_k];
    // Default the conditional projection stage to "no change" so a checkpoint
    // iter that skips the heavy block still yields proj-Δ = 0 (not stale).
    if(s == 1){ su[2] = su[1]; sv[2] = sv[1]; sw[2] = sw[1]; }
}

void cHydrosphereModel::write_stage_budget(int iter){
    if(scell_i < 0) return;
    const double sc = u_0 * 100.0;   // nondim velocity -> cm/s (per-iteration Δ)

    // Per-stage change in each velocity component [cm/s over this iteration].
    const double rk4_u = (su[1]-su[0])*sc, rk4_v = (sv[1]-sv[0])*sc, rk4_w = (sw[1]-sw[0])*sc;
    const double prj_u = (su[2]-su[1])*sc, prj_v = (sv[2]-sv[1])*sc, prj_w = (sw[2]-sw[1])*sc;
    const double bc_u  = (su[3]-su[2])*sc, bc_v  = (sv[3]-sv[2])*sc, bc_w  = (sw[3]-sw[2])*sc;
    const double flt_u = (su[4]-su[3])*sc, flt_v = (sv[4]-sv[3])*sc, flt_w = (sw[4]-sw[3])*sc;
    const double net_u = (su[4]-su[0])*sc, net_v = (sv[4]-sv[0])*sc, net_w = (sw[4]-sw[0])*sc;

    const double lat = 90.0 - (double)scell_j * 180.0 / (double)(jm - 1);
    const double lon = (double)scell_k * 360.0 / (double)(km - 1);
    const double depth = (rad.z[im-1] - rad.z[scell_i]) * L_hyd;

    std::ostringstream fname;
    fname << output_path << "/stage_budget.csv";
    std::ofstream f(fname.str().c_str(), std::ios::app);
    if(f.is_open()){
        if(f.tellp() == std::streampos(0))
            f << "iter,i,j,k,lat_deg,lon_deg,depth_m,"
              << "rk4_u,proj_u,bc_u,filter_u,net_u,"
              << "rk4_v,proj_v,bc_v,filter_v,net_v,"
              << "rk4_w,proj_w,bc_w,filter_w,net_w\n";
        f << iter << "," << scell_i << "," << scell_j << "," << scell_k << ","
          << lat << "," << lon << "," << depth << ","
          << rk4_u << "," << prj_u << "," << bc_u << "," << flt_u << "," << net_u << ","
          << rk4_v << "," << prj_v << "," << bc_v << "," << flt_v << "," << net_v << ","
          << rk4_w << "," << prj_w << "," << bc_w << "," << flt_w << "," << net_w << "\n";
        f.close();
    }

    cout << "      [stagebud] iter=" << iter << fixed << setprecision(1)
         << "  @i=" << scell_i << " lat=" << lat << " lon=" << lon
         << " depth=" << setprecision(0) << depth << "m" << setprecision(4)
         << "  Du: rk4=" << rk4_u << " proj=" << prj_u << " bc=" << bc_u
         << " filt=" << flt_u << " NET=" << net_u
         << "  | Dw: rk4=" << rk4_w << " proj=" << prj_w << " bc=" << bc_w
         << " filt=" << flt_w << " NET=" << net_w << " cm/s/iter" << endl;
}
