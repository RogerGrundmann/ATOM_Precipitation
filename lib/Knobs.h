#pragma once
// Knobs.h -- THE registry of runtime environment knobs (ATM_*, HYD_*, ATOM_*). KNOB-INV plan B, 2026-09-30.
//
// Every environment switch in atmosphere/, hydrosphere/ and lib/ is declared ONCE, in ATOM_KNOB_LIST below:
// name, default (as text), Result/Diag, the class from docs/knob_inventory.md (F X N R P S I D), a one-line note.
// Code reads a knob only through the accessors -- knob::on / integer / real / text / is_set -- never getenv().
// That removes three defects the inventory found (docs/knob_inventory.md, "Defects found by the inventory"):
//   1. a result-changing knob missing from the [RUN CONFIG] banner: the banner is GENERATED from this list;
//   2. one knob parsed differently at different sites (ATM_POISSON_METRIC_FIX was atof at one site and atoi at
//      two): each accessor parses one way, everywhere;
//   3. a default copied into every site that reads the knob: it lives here once.
// Validation that belongs to one use (clamps, "0 means the shipped value") stays at the site.
//
// TO ADD A KNOB: add one X(...) row, read it with knob::<accessor>(knob::NAME). It is in the banner automatically.
// TO FLIP A DEFAULT: change the text in its row. Nothing else holds a copy.
//
// SOURCES, in order of precedence (plan D, 2026-09-30): the ENVIRONMENT, then the <knobs> section of the XML config
// (<atom><knobs><ATM_RH_MIN_PTOP>482</ATM_RH_MIN_PTOP></knobs></atom>, loaded by each model's LoadConfig), then the
// compiled default below. An unknown name in <knobs> is an error (a typo must not silently run the default), and so
// is loading <knobs> after any knob has already been read (the value would be ignored). The banner marks the source.
//
// Parsing: on() = atoi(value) != 0; integer() = atoi; real() = atof; text() = the raw string. A default that is
// not a number (ATM_METRIC_RADIUS = "r_Earth") must be read with is_set() first -- the site supplies the value --
// and real()/integer() abort if asked to parse it. Values are read from the environment once, at first use, and
// cached; nothing in the tree sets environment variables at runtime.
// Byte identity of the conversion: atof/strtod and the compiler both round a decimal literal correctly, so
// real("0.08538") is the same double as the literal 0.08538 it replaced.

#include <array>
#include <atomic>
#include <stdexcept>
#include <utility>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>   // environ
#include <initializer_list>
#include <sstream>
#include <string>

namespace knob {

enum Kind { Result, Diag };   // Diag = print, dump or output cadence only; cannot change a result

// X(name, default, kind, inventory class, note)
#define ATOM_KNOB_LIST(X) \
    X(ATOM_CORIOLIS_NONTRAD, "0", Result, 'S', "non-traditional Coriolis terms (both models)") \
    X(ATOM_METRIC_CURVATURE, "0", Result, 'S', "spherical curvature terms in the shared metric (both models)") \
    X(ATOM_METRIC_DIVERGENCE, "0", Result, 'S', "metric terms in the shared divergence (both models)") \
    X(ATM_CELL_ROT_DIAG, "0", Diag, 'D', "print/dump only") \
    X(ATM_CLOUD_INIT_DIAG, "0", Diag, 'D', "print/dump only") \
    X(ATM_CO2_BAND, "0.17", Result, 'P', "0.17") \
    X(ATM_CONV_ADJ, "0", Result, 'R', "off; palliative") \
    X(ATM_CONV_ADJ_LAPSE, "1.0", Result, 'P', "sub-parameter of CONV_ADJ") \
    X(ATM_CONV_ADJ_PASSES, "64", Result, 'P', "sub-parameter of CONV_ADJ") \
    X(ATM_CWB_BANDS, "0", Diag, 'D', "print/dump only") \
    X(ATM_CWB_DIAG, "0", Diag, 'D', "print/dump only") \
    X(ATM_CWP_CENSUS, "0", Diag, 'D', "print/dump only") \
    X(ATM_DAMP_Q_HORIZ, "1", Result, 'R', "step A; qh600 scored") \
    X(ATM_DAMP_Q_MASS, "0", Result, 'R', "forced by WATER_CLOSURE") \
    X(ATM_DAMP_Q_VERT, "1", Result, 'R', "step A") \
    X(ATM_DAMP_T_HORIZ, "1", Result, 'R', "qth_on today") \
    X(ATM_DAMP_T_VERT, "1", Result, 'R', "null alone (qvt_on)") \
    X(ATM_EPS_DRY, "0.684", Result, 'P', "0.684") \
    X(ATM_HADLEY_SL, "4.0", Result, 'P', "IC 4.0N/3.0S") \
    X(ATM_EVAP_WIND, "0", Result, 'R', "bulk evaporation wind speed: 0 = shipped sqrt((u2+v2+w2)/3), the rms of the components = |V|/sqrt(3); 1 = the wind speed |V|") \
    X(ATM_EVAP_GUST, "0", Result, 'P', "bulk evaporation: gustiness in m/s added in quadrature to the (annual-mean vector) wind speed, W = sqrt(|V|^2 + g^2); 0 = off") \
    X(ATM_LAND_EVAP, "0", Result, 'R', "land evaporation from the Budyko curve on the model's own rain: E = P*sqrt(phi*tanh(1/phi)*(1-exp(-phi))), phi = E_pot/P (E <= P and E <= E_pot; needs ATM_WATER_CLOSURE, ignored with ATM_LAND_BUCKET)") \
    X(ATM_LAND_BUCKET, "0.0", Result, 'R', "20-iter evidence only") \
    X(ATM_LENGTH_NDIM, "0", Result, 'R', "the 40x L_atm defect (item 4 pending)") \
    X(ATM_LONGAL_J, "62", Diag, 'I', "output slice latitude") \
    X(ATM_MC_ALF1, "0.05", Result, 'R', "sub-cloud rain evaporation rate 1/s (Tiedtke 5.44e-4; working branch)") \
    X(ATM_MC_BASE_SAT, "0", Result, 'R', "=2 in gpa_* today") \
    X(ATM_MC_CAP_DIAG, "0", Diag, 'D', "print/dump only") \
    X(ATM_MC_COND_DEBIT, "0", Result, 'R', "charge the convective condensate to the environment: c_u = g_p + e_l (MC-Q-LEAK)") \
    X(ATM_MC_DIAG, "0", Diag, 'D', "print/dump only") \
    X(ATM_MC_T_COEFF, "0", Result, 'R', "coeff_MC_t = metricShellLength()/(u_0*t_0) (was L_atm: 0.025x; MC-TV)") \
    X(ATM_MC_VEL_COEFF, "0", Result, 'R', "coeff_MC_vel = metricShellLength()/u_0^2 (was L_atm: 0.025x); pair with ATM_MC_UV_DETRAIN (MC-TV)") \
    X(ATM_MC_QC_DETRAIN, "0", Result, 'R', "updraft q_c_u recurrence detrains q_c_u, not the environment's cloud (MC-TV)") \
    X(ATM_MC_UV_DETRAIN, "0", Result, 'R', "updraft v_u/w_u recurrence loses -D_u*v_u, -D_u*w_u like q_v_u and s_u (MC-TV)") \
    X(ATM_MC_Q_NDIM, "0", Result, 'R', "coeff_MC_q = metricShellLength()/u_0, as the microphysics (was L_atm/(u_0*c_0), 0.713x)") \
    X(ATM_MC_DEPTH_RAMP, "0", Result, 'R', "convective precipitation scaled by cloud depth (0 at 100 hPa -> 1 at 300 hPa) instead of the 200 hPa shallow switch") \
    X(ATM_MC_GATE_BLEND, "0", Result, 'R', "smoothstep blend env->parcel over |M| 0.01..0.1 instead of the recurrence gate switch (1 updraft, 2 also downdraft; 3 = no blend, gate threshold 0.01, parcel undiluted)") \
    X(ATM_MC_ED_ABOVE_BASE, "0", Result, 'R', "downdraft evaporation e_d only above cloud base; e_p is the single sub-cloud term") \
    X(ATM_MC_ED_AREA, "0", Result, 'R', "area-weight the convective downdraft evaporation (working branch 1)") \
    X(ATM_MC_MB_SAT_LAND, "0", Result, 'P', "land columns: cloud-base mass flux saturates, M_b -> M_s*tanh(M_b/M_s), M_s in kg/(m2 s) (0 = off)") \
    X(ATM_MC_MB_SAT_OCEAN, "0", Result, 'P', "ocean columns: cloud-base mass flux saturates, M_b -> M_s*tanh(M_b/M_s), M_s in kg/(m2 s) (0 = off)") \
    X(ATM_MC_ML_LCL, "0", Result, 'R', "with ATM_MC_ML_PARCEL: convective base = the mixed-layer parcel's LCL (not the stratiform-cloud base)") \
    X(ATM_MC_ML_PARCEL, "0", Result, 'R', "deep-convection parcel from the 0-500 m mixed layer (1: dq = ATM_MC_Q_ADD, 2: dq = sub-grid sigma_q) instead of the saturated cloud-base env") \
    X(ATM_MC_Q_ADD, "1.0e-4", Result, 'P', "parcel moisture excess kg/kg (shipped constant q_v_u_add)") \
    X(ATM_MC_T_ADD_LAND, "0", Result, 'P', "with ATM_MC_ML_PARCEL: additional parcel temperature excess over land columns, K (surface heating; ocean parcel unchanged)") \
    X(ATM_MC_T_ADD_LAND_L, "0", Result, 'P', "ATM_MC_T_ADD_LAND: coastal taper, excess x (1 - exp(-d/L)), d = distance to the ocean, L in km (0 = no taper)") \
    X(ATM_MC_T_ADD, "0.2", Result, 'P', "parcel temperature excess K (shipped constant t_add_u)") \
    X(ATM_MC_ENTR, "0.0", Result, 'R', "updraft entrainment 1/m; 0 = the shipped 0.2/R_cloud") \
    X(ATM_MC_GP_AREA, "0", Result, 'R', "gpa_* today") \
    X(ATM_MC_QVD, "0", Result, 'R', "B.10c") \
    X(ATM_MC_SGZ, "0", Result, 'R', "B.10b") \
    X(ATM_METRIC_CHECK, "0", Diag, 'D', "print/dump only") \
    X(ATM_METRIC_EXACT, "0", Result, 'E', "kept as a documented experiment (user 2026-09-30); null on integrated quantities; undecidable") \
    X(ATM_METRIC_SIN_FLOOR, "0.26", Result, 'P', "0.26 since 09-10") \
    X(ATM_METRIC_STRICT, "0", Diag, 'I', "abort on metric check") \
    X(ATM_MFC_DIAG, "0", Diag, 'D', "print/dump only") \
    X(ATM_PDYN_CAP, "2.0", Result, 'P', "p_dyn source cap (non-dim)") \
    X(ATM_PDYN_CEILING, "0.0", Result, 'P', "p_dyn clamp; 0 = phase-dependent (10 before iteration 300, 3 after)") \
    X(ATM_POISSON_METRIC_FIX, "0", Result, 'E', "kept as a documented experiment (user 2026-09-30); consistent horizontal Poisson metric; parsed as an integer at all sites since 2026-09-30 (was atof at one)") \
    X(ATM_POLAR_CELL_SHEAR, "0.1", Result, 'P', "IC 0.1") \
    X(ATM_PRECIP_PASSES, "3", Result, 'R', "rain-column passes; default is TwoCatIce::iter_prec_end (3), the site falls back to it") \
    X(ATM_PRECIP_RELAX, "1.0", Result, 'R', "under-relaxation of the rain-column passes, (0,1]") \
    X(ATM_PRECIP_UPWIND, "0", Result, 'R', "converged upwind rain column (working branch 1)") \
    X(ATM_PRESS_LINE_SOLVE, "0", Result, 'E', "kept as a documented experiment (user 2026-09-30); measured through the clamp; no runaway fix") \
    X(ATM_PRESS_SWEEPS, "1", Result, 'P', "1x") \
    X(ATM_PROJECT_IN_LOOP, "0", Result, 'E', "kept as a documented experiment (user 2026-09-30); null at 10 and 200 sweeps") \
    X(ATM_PROJ_CONSISTENCY, "0", Diag, 'D', "print/dump only") \
    X(ATM_PSI_FERREL, "40.0", Result, 'P', "Ferrel cell amplitude for ATM_CELLS_FROM_PSI, 1e9 kg/s") \
    X(ATM_PSI_HADLEY, "120.0", Result, 'P', "Hadley cell amplitude for ATM_CELLS_FROM_PSI, 1e9 kg/s") \
    X(ATM_PSI_POLAR, "26.0", Result, 'P', "polar cell amplitude for ATM_CELLS_FROM_PSI, 1e9 kg/s") \
    X(ATM_PSI_PROJ_DUMP, "0", Diag, 'D', "print/dump only") \
    X(ATM_PSI_SHAPE, "1", Result, 'P', "IC shape 1") \
    X(ATM_Q_DIFF_FLUX, "0", Result, 'R', "conservative (flux-form, rho-weighted) vertical diffusion of c/cloud/ice/gr (DIFF-LEAK)") \
    X(ATM_QC_CRIT, "0.05", Result, 'P', "0.05 g/kg since 08-31") \
    X(ATM_RADIAL_SHAPIRO_STRENGTH, "1.0", Result, 'P', "1.0 (u only)") \
    X(ATM_RADIATION_MODE, "5", Result, 'S', "5 (radiation diagnostic)") \
    X(ATM_RAD_COLDIAG, "0", Diag, 'D', "print/dump only") \
    X(ATM_RAD_EQUIL, "0", Result, 'R', "pair with SW_INSOL; needs a prognostic T") \
    X(ATM_RAIN_AREA, "0.10", Result, 'P', "0.10 since 09-01 (fitted)") \
    X(ATM_RAIN_PASS_DIAG, "0", Diag, 'D', "print the rain-column pass convergence") \
    X(ATM_RESTART_KEEP, "", Diag, 'I', "=all keeps every periodic restart file (default: only the latest)") \
    X(ATM_RESTART_STRIDE, "100", Diag, 'I', "restart cadence") \
    X(ATM_RH_CRIT, "0.30", Result, 'P', "0.30 since 08-31") \
    X(ATM_RH_MIN, "0.65", Result, 'P', "0.65") \
    X(ATM_RH_MIN_PTOP, "482.0", Result, 'P', "482 hPa (fitted)") \
    X(ATM_RH_OCEAN_SST, "0", Result, 'P', "initial ocean surface RH follows the sea-surface temperature: ATM_RH_OCEAN at ATM_RH_OCEAN_SST_REF, lower by this much per K colder (1/K; 0 = off, one RH on every tropical ocean column)") \
    X(ATM_RH_OCEAN_SST_REF, "29.5", Result, 'P', "ATM_RH_OCEAN_SST: sea-surface temperature at and above which the surface RH is ATM_RH_OCEAN, deg C") \
    X(ATM_RH_OCEAN_SST_MAX, "0.03", Result, 'P', "ATM_RH_OCEAN_SST: largest RH reduction (cold-water floor)") \
    X(ATM_RH_OCEAN_SST_COLD, "-99", Result, 'P', "ATM_RH_OCEAN_SST: cool-side cut-off, deg C -- no reduction at and below it, the full one from 1 K above it (-99 = no cut-off)") \
    X(ATM_RH_LAND, "0", Result, 'R', "land initial surface RH from subtropical descent x continentality (paleo-safe; 0 = every land column 0.75; 2 = wet end follows ATM_RH_OCEAN in the tropics)") \
    X(ATM_RH_LAND_QCAP, "0", Result, 'R', "cap initial land vapour at the zonal-mean ocean value at the same height above the surface (land hotter than the sea cannot be moister)") \
    X(ATM_RH_SIGMA_SFC, "0", Result, 'R', "initial RH profile sigma = p / local ground pressure on land (was p/p_0: elevated ground starts drier)") \
    X(ATM_HCRIT_SFC, "0", Result, 'P', "land columns: the stratiform cloud threshold H_crit is evaluated at sigma = p / the column's own surface pressure (H_crit -> 1 at the ground on a plateau as at sea level); not the convection scheme") \
    X(ATM_HCRIT_SFC_LAT, "0", Result, 'P', "ATM_HCRIT_SFC: full weight equatorward of this latitude, none from 10 deg poleward of it (deg; 0 = every latitude)") \
    X(ATM_RH_SIGMA_LAT, "0", Result, 'R', "poleward of this latitude (deg, 10 deg taper) the initial RH profile uses sigma = p / the column's own surface pressure, land columns (0 = off)") \
    X(ATM_RH_SIGMA_LAT_OCEAN, "0", Result, 'R', "ATM_RH_SIGMA_LAT also on ocean columns") \
    X(ATM_RH_LAND_EAST, "0", Result, 'P', "ATM_RH_LAND: strength (0..1) by which open ocean to the EAST (trade-wind fetch, 2000 km) cancels the subtropical descent drying") \
    X(ATM_RH_LAND_EAST_MAX, "-1", Result, 'P', "ATM_RH_LAND_EAST: largest initial land surface RH above the wet end that the east fetch may give (0 = clamp at the wet end; -1 = no clamp)") \
    X(ATM_RH_LAND_EAST_ML, "0", Result, 'R', "ATM_RH_LAND_EAST: subtropical land (|lat| 10..20 taper in, 30..40 out) too cool to convect, weighted by the ocean fetch to the east: well-mixed initial boundary layer of this depth above the ground in m (0 = off)") \
    X(ATM_RH_LAND_EAST_ML_STRENGTH, "1", Result, 'P', "ATM_RH_LAND_EAST_ML: partial mixing, fraction (0..1) of the way from the Manabe-Wetherald RH to the well-mixed one (times the fetch m_east)") \
    X(ATM_RH_LAND_EAST_ML_T, "24", Result, 'P', "ATM_RH_LAND_EAST_ML: no mixed layer where the ground is at least this warm (deg C; such land convects), full weight 2 K below") \
    X(ATM_RH_LAND_DRY, "0.35", Result, 'P', "ATM_RH_LAND: driest land surface RH (descent core, far inland)") \
    X(ATM_RH_LAND_L, "1500", Result, 'P', "ATM_RH_LAND: continentality e-folding distance from the ocean, km") \
    X(ATM_RH_OCEAN, "0.75", Result, 'R', "initial surface RH over tropical ocean (|lat|<=30, taper to 45); the marine BL keeps it on affordable runs") \
    X(ATM_RH_OCEAN_ML, "0", Result, 'R', "extratropical ocean (|lat| 30..40 taper): well-mixed initial boundary layer of this depth in m, constant q from the surface, RH capped 0.98 (0 = off)") \
    X(ATM_RH_OCEAN_ML_STRENGTH, "1", Result, 'P', "ATM_RH_OCEAN_ML: partial mixing, fraction (0..1) of the way from the Manabe-Wetherald RH to the well-mixed one") \
    X(ATM_RH_STORM, "1.0", Result, 'R', "initial storm-track RH factor at 55 deg (working branch 1.15)") \
    X(ATM_RH_STORM_LAT, "55", Result, 'P', "ATM_RH_STORM: latitude of the factor's maximum, deg (also where ATM_RH_STORM_POLAR starts)") \
    X(ATM_RH_STORM_WIDTH, "15", Result, 'P', "ATM_RH_STORM: Gaussian width of the factor, deg") \
    X(ATM_RH_STORM_POLAR, "1.0", Result, 'R', "initial RH factor over OCEAN poleward of 55 deg: the ATM_RH_STORM factor does not fall below this toward the pole (1.0 = off, the Gaussian alone)") \
    X(ATM_RH_STORM_SST, "0", Result, 'P', "extratropical ocean (taper in over 10 deg from ATM_RH_STORM_SST_LAT): the initial RH is lowered by this fraction per K of sea-surface temperature above ATM_RH_STORM_SST_REF and raised per K below it (0 = off)") \
    X(ATM_RH_STORM_SST_LAT, "40", Result, 'P', "ATM_RH_STORM_SST: latitude where it starts, deg; full weight 10 deg poleward") \
    X(ATM_RH_STORM_SST_REF, "4", Result, 'P', "ATM_RH_STORM_SST: reference sea-surface temperature, deg C") \
    X(ATM_RH_STORM_SST_MAX, "0.05", Result, 'P', "ATM_RH_STORM_SST: largest relative RH change, either sign") \
    X(ATM_RK_SCALAR_SYNC, "0", Result, 'R', "forced by WATER_CLOSURE") \
    X(ATM_SATADJ_DIAG, "0", Diag, 'D', "print/dump only") \
    X(ATM_SATADJ_FADE, "0", Result, 'R', "mode 2 if ever flipped; forced by closure") \
    X(ATM_SEAM_Q_CONSERVE, "0", Result, 'R', "conserve water at the phi seam; 2 = also next to land (working branch 2)") \
    X(ATM_SEAM_Q_DIAG, "0", Diag, 'D', "print the seam water-conservation buckets") \
    X(ATM_SNOW_DEP_FLUX, "0", Result, 'R', "snow deposition/sublimation and S_i_cri also move the snow FLUX (SNOW-SUBL)") \
    X(ATM_SNOW_DIAG, "0", Diag, 'D', "print the snow budget") \
    X(ATM_SNOW_WINDOW, "0", Result, 'R', "snow kept on the cold side of -20 C (working branch 2)") \
    X(ATM_SR_DIAG, "0", Diag, 'D', "print/dump only") \
    X(ATM_SURF_DRAG_CONSISTENT, "1.0", Result, 'R', "B.9, default 1.0 since 2026-09-29 (run_sdr.sh: connected, climate null)") \
    X(ATM_SW_INSOL, "0.0", Result, 'R', "pair with RAD_EQUIL") \
    X(ATM_T0_ATTRIB, "0", Diag, 'D', "print/dump only") \
    X(ATM_TEQ_SKIN_ONLY, "0", Result, 'R', "instrument branch, must NOT be flipped") \
    X(ATM_TEQ_WTG, "0", Result, 'R', "blend the initial free troposphere (1.5-3 km above ground up, |lat|<30, taper to 45) to the zonal mean; t_eq inherits it") \
    X(ATM_TW_BALANCE, "0.0", Result, 'R', "off; the only mid-lat jet IC") \
    X(ATM_TW_BALANCE_V, "0", Result, 'P', "sub-option of TW_BALANCE") \
    X(ATM_TW_LATMIN, "15.0", Result, 'P', "sub-parameter of TW_BALANCE") \
    X(ATM_TW_WMAX, "80.0", Result, 'P', "sub-parameter of TW_BALANCE") \
    X(ATM_T_FLOOR, "216.65", Result, 'P', "216.65 K since 08-31") \
    X(ATM_UBUD_BALANCE, "0", Diag, 'D', "print/dump only") \
    X(ATM_VTK_STRIDE, "5", Diag, 'I', "VTK cadence") \
    X(ATM_WATER_CLOSURE, "0", Result, 'R', "reverted 09-25; works with filter off (qh600)") \
    X(HYD_A_H, "0.0", Result, 'R', "Laplacian; biharmonic preferred") \
    X(HYD_A_H_BIHARM, "1.0e19", Result, 'F', "default 1.0e19 since 09-30 (OCN-METRIC, om_b1e19: noise 0.515, radial u 2400x below shipped; 3e19 ran away)") \
    X(HYD_A_H_BIHARM_SCALED, "1", Result, 'F', "default 1 since 09-30 (sin^4 scaling; unscaled 3e18 NaN'd at the pole with floor 0.26)") \
    X(HYD_BAROCLINIC_PGF, "0.0", Result, 'R', "pair with PHYDRO_SALT; HYDRO_SPLIT may supersede") \
    X(HYD_BC_DRAG, "0.0", Result, 'R', "numerical, not physical") \
    X(HYD_BUOY_CONSISTENT, "0.0", Result, 'R', "radial runaway; HYDRO_SPLIT") \
    X(HYD_KE_SPLIT, "0", Diag, 'D', "print/dump only") \
    X(HYD_METRIC_RADIUS, "6370.0", Result, 'F', "default 6370 km since 09-30 (user); profile stays bottom-intensified (structural, OCN-PROF)") \
    X(HYD_METRIC_SIN_FLOOR, "0.26", Result, 'P', "0.26 since 09-28, was 0.4 (ro_osf26/40: 75-90 KE -7.5 %, equatorward null; matches the atmosphere)") \
    X(HYD_NUE_GRAD, "0.0", Result, 'R', "off; would inherit the broken metric") \
    X(HYD_PHYDRO_SALT, "0", Result, 'R', "pair") \
    X(HYD_RESTART_KEEP, "", Diag, 'I', "=all keeps every periodic ocean restart file (default: only the latest)") \
    X(HYD_RUN_NEUMANN, "1", Result, 'F', "default 1 since 09-30, with METRIC_RADIUS") \
    X(HYD_SFC_FLUX, "0.0", Result, 'R', "unmeasured")

enum Id : int {
#define ATOM_KNOB_ID(n, d, k, c, doc) n,
    ATOM_KNOB_LIST(ATOM_KNOB_ID)
#undef ATOM_KNOB_ID
    N_KNOBS
};

struct Spec { const char* name; const char* dflt; Kind kind; char cls; const char* doc; };

inline constexpr Spec specs[N_KNOBS] = {
#define ATOM_KNOB_SPEC(n, d, k, c, doc) {#n, d, k, c, doc},
    ATOM_KNOB_LIST(ATOM_KNOB_SPEC)
#undef ATOM_KNOB_SPEC
};

// The <knobs> section of the XML config (see SOURCES above).
struct ConfigSource {
    std::array<std::string, N_KNOBS> val;
    std::array<bool, N_KNOBS> set{};
    std::atomic<bool> reads_started{false};
    std::string file;
};
inline ConfigSource& config() { static ConfigSource c; return c; }

inline int find(const char* name) {
    for (int i = 0; i < N_KNOBS; i++) if (std::strcmp(specs[i].name, name) == 0) return i;
    return -1;
}

// Called by LoadConfig with the (name, text) pairs of <knobs>. Must run before any knob is read.
inline void load_config(const std::vector<std::pair<std::string, std::string>>& kv, const std::string& file) {
    if (kv.empty()) return;
    ConfigSource& c = config();
    if (c.reads_started.load())
        throw std::logic_error("config " + file + ": <knobs> loaded after a knob was already read -- it would be ignored");
    for (const auto& p : kv) {
        const int i = find(p.first.c_str());
        if (i < 0)
            throw std::invalid_argument("config " + file + ": <knobs> names an unknown knob '" + p.first
                                        + "' (retired, or a typo? the list is lib/Knobs.h)");
        c.val[i] = p.second; c.set[i] = true;
    }
    c.file = file;
}

// The environment, read once for all knobs (thread-safe static init).
inline const char* env(Id id) {
    static const std::array<const char*, N_KNOBS> cache = [](){
        std::array<const char*, N_KNOBS> a{};
        for (int i = 0; i < N_KNOBS; i++) a[i] = std::getenv(specs[i].name);
        return a; }();
    return cache[id];
}

// Where the value in force comes from: 'e' environment, 'x' XML config, 'd' compiled default.
inline char source(Id id) {
    config().reads_started.store(true);
    if (env(id)) return 'e';
    return config().set[id] ? 'x' : 'd';
}

inline bool is_set(Id id) { return source(id) != 'd'; }

// The value in force as text: the environment's, else the XML config's, else the registry default.
inline const char* value(Id id) {
    switch (source(id)) {
        case 'e': return env(id);
        case 'x': return config().val[id].c_str();
        default:  return specs[id].dflt;
    }
}

inline const char* numeric_value(Id id) {
    const char* v = value(id);
    char* end = nullptr;
    std::strtod(v, &end);
    if (end == v) {
        std::fprintf(stderr, "knob::%s: default \"%s\" is not a number -- read it with knob::is_set() first\n",
                     specs[id].name, specs[id].dflt);
        std::abort();
    }
    return v;
}

inline bool        on(Id id)      { return std::atoi(numeric_value(id)) != 0; }
inline int         integer(Id id) { return std::atoi(numeric_value(id)); }
inline double      real(Id id)    { return std::atof(numeric_value(id)); }
inline std::string text(Id id)    { return value(id); }

// [RUN CONFIG] lines for every knob whose name starts with one of the prefixes. Result knobs are listed in full,
// SHORT=value, with '*' when the value is the compiled-in default (not set in the environment); Diag knobs only
// when set. `tag` is the log prefix ("      AGCM: "), `per_line` the entries per line.
inline std::string banner(const char* tag, std::initializer_list<const char*> prefixes, int per_line = 8) {
    auto want = [&](const char* n){
        for (const char* p : prefixes) if (std::strncmp(n, p, std::strlen(p)) == 0) return true;
        return false; };
    auto entry = [](int i){
        const char* n = specs[i].name;
        const char* s = std::strchr(n, '_');
        std::string e = std::string(s ? s + 1 : n) + "=";
        const Id id = static_cast<Id>(i);
        const char* v = value(id);
        e += v[0] ? v : "\"\"";
        const char src = source(id);
        if (src == 'd') e += "*"; else if (src == 'x') e += "+";
        return e; };
    std::ostringstream b;
    int n = 0;
    for (int i = 0; i < N_KNOBS; i++) {
        if (!want(specs[i].name) || specs[i].kind != Result) continue;
        if (n % per_line == 0) b << (n ? "\n" : "") << tag << "[RUN CONFIG] knobs:";
        b << "  " << entry(i);
        n++;
    }
    b << "\n" << tag << "[RUN CONFIG] diagnostic/output knobs set:";
    int d = 0;
    for (int i = 0; i < N_KNOBS; i++)
        if (want(specs[i].name) && specs[i].kind == Diag && is_set(static_cast<Id>(i))) { b << "  " << entry(i); d++; }
    if (!d) b << "  none";
    // An environment variable that LOOKS like a knob but is not one -- retired in plan C, or a typo -- is otherwise
    // ignored in silence, which is how a retired knob in an old run script would quietly run the default.
    {
        std::string unknown;
        for (char** ev = environ; ev && *ev; ev++) {
            const std::string kv(*ev);
            const std::string name = kv.substr(0, kv.find('='));
            if (!want(name.c_str()) || find(name.c_str()) >= 0) continue;
            if (name.rfind("ATM_", 0) == 0 || name.rfind("HYD_", 0) == 0 || name.rfind("ATOM_", 0) == 0)
                unknown += "  " + name;
        }
        if (!unknown.empty())
            b << "\n" << tag << "[RUN CONFIG] *** WARNING: set in the environment but NOT a knob (retired or a typo;"
              << " IGNORED):" << unknown;
    }
    b << "\n" << tag << "[RUN CONFIG] (* = compiled-in default, + = from <knobs> in "
      << (config().file.empty() ? std::string("the XML config (none here)") : config().file)
      << ", no mark = environment; the list is lib/Knobs.h)\n";
    return b.str();
}

}  // namespace knob
