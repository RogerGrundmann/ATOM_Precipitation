# Runtime knob inventory (2026-09-26)

Every `getenv(...)` switch in `atmosphere/`, `hydrosphere/`, `lib/`, extracted from the source (not recited) at
`55d889d`, and classified from `CLAUDE.md` and the bug register. **144 knobs**: 119 `ATM_*`, 22 `HYD_*`, 3 `ATOM_*`.

| class | meaning | count | proposal |
|---|---|---|---|
| F | Decided ON by measurement | 29 | Delete the old branch (git keeps it). The default branch stays byte-identical; a byte check per batch proves it. |
| X | Switched on, then reverted by measurement | 1 | Delete the knob and its code path (it is superseded). |
| N | Off; measured null, refuted or superseded | 22 | Delete, unless it is wanted as a documented experiment -- then move it to the config file. |
| R | Open: repair or experiment, runs owed or decision pending | 41 | Keep as a switch until decided; then it becomes F or N. |
| P | Numerical parameter (a value, not a switch) | 25 | Move into the XML config / param.py, where parameters belong and are recorded with the run. |
| S | Structural option at its shipped value | 4 | Move into the config file. |
| I | Output / run infrastructure | 4 | Move into the config file. |
| D | Diagnostic, print or dump only | 18 | Keep as environment switches -- this is what env knobs are good for. |

**If the proposals are followed**: F + X + N = 52 knobs disappear, P + S + I = 33 move into the config file, R = 41 stay as
switches until decided, D = 18 stay as diagnostic switches.

## Defects found by the inventory itself

1. **32 knobs that can change results are missing from the `[RUN CONFIG]` banner**, so a run log does not record them:
   `ATM_ANELASTIC`, `ATM_BUOY_MOIST`, `ATM_BUOY_TREF`, `ATM_CLOUD_TAU_MAX`, `ATM_CO2_BAND`, `ATM_CONV_ADJ`, `ATM_CONV_ADJ_LAPSE`, `ATM_CONV_ADJ_PASSES`, `ATM_EPS_DRY`, `ATM_GRID_BETA`, `ATM_GRID_PRESSURE`, `ATM_GRID_PTOP`, `ATM_LENGTH_NDIM`, `ATM_METRIC_EXACT`, `ATM_METRIC_NOCURV`, `ATM_METRIC_RADIUS`, `ATM_PRESS_SWEEPS`, `ATM_PROJECT_IN_LOOP`, `ATM_PROJ_SWEEPS`, `ATM_RADIATION_MODE`, `ATM_RHIE_CHOW`, `ATM_SFC_FLUX`, `ATM_TAU_PBROAD`, `ATM_TEQ_SKIN_ONLY`, `ATM_TW_WMAX`, `ATOM_CORIOLIS_NONTRAD`, `ATOM_METRIC_CURVATURE`, `ATOM_METRIC_DIVERGENCE`, `HYD_BC_DRAG`, `HYD_LINE_FOLD`, `HYD_LINE_GAUGE`, `HYD_NUE_GRAD`.
   Among them live physics defaults (`ATM_CLOUD_TAU_MAX` = 2.0 is ON, `ATM_METRIC_RADIUS`, `ATM_RADIATION_MODE` = 5,
   `ATM_EPS_DRY`, `ATM_CO2_BAND`) and the grid options. (Two output-only knobs, `ATM_LONGAL_J` and `ATM_METRIC_STRICT`,
   are also missing and do not affect results.)
2. **`ATM_POISSON_METRIC_FIX` is parsed differently at its three sites**: `PressureSolverAtm.h:315` tests `atof(e) != 0.0`,
   `:137` and `:1324` test `atoi(e) != 0`, so e.g. `=0.5` is ON in one place and OFF in the other two. 0/1 are unaffected.
3. **22 knobs are read at more than one site**, each with its own copy of the default (`ATM_WATER_CLOSURE` 4, `ATM_MICRO_NDIM`
   3, ...). All copies currently agree; a flip that misses one would silently split physics from diagnostics.
4. Knobs are cached in function-local `static` lambdas, so the value is fixed at first use; fine for runs, but it means
   the banner (read separately) and the physics could disagree if a knob were ever set programmatically mid-run.

## Proposed next steps (after tomorrow's arms, so pending byte checks are not compared against moving code)

- **B. One registry header** (`atmosphere/Knobs.h`, `hydrosphere/HydKnobs.h`): each knob defined once -- name, default, type,
  one-line doc -- read through one accessor; the banner generated from it (fixes defects 1-3 by construction).
- **C. Retire F / X / N in batches** of ~10, each with a 1-thread byte check of the default branch.
- **D. Move P / S / I into the XML config** via `param.py`, so a run's configuration lives in one recorded place.

## F -- Decided ON by measurement (29)

| knob | read at | sites | banner default | note |
|---|---|---|---|---|
| `ATM_BC_SECOND_ORDER` | `atmosphere/cAtmosphereModel.h:169` | 1 | `1*` | 1; 2nd-order Neumann, 09-23 (B.5) |
| `ATM_CELLS_FROM_PSI` | `atmosphere/VelocityInitializer.h:352` | 2 | `1*` | 1; closed-cell IC, 09-12 |
| `ATM_CLOUD_FRAC` | `atmosphere/CloudFraction.h:44` | 2 | `1*` | 1; sub-grid cloud, 08-31 |
| `ATM_CLOUD_RAD_FRAC` | `atmosphere/MultiLayerRadiation.h:161` | 1 | `1*` | 1; 08-31 |
| `ATM_CLOUD_TAU_MAX` | `atmosphere/MultiLayerRadiation.h:636` | 1 | **missing** | 2.0; layer tau ceiling, 08-28 -- NOT IN BANNER |
| `ATM_CWP_CAP` | `atmosphere/MultiLayerRadiation.h:498` | 2 | `off*` | disabled (1e9), 08-31 |
| `ATM_EVAP_STRIDE_FIX` | `atmosphere/ThermoAtm.h:219` | 1 | `1*` | 1; 09-24 (acts only under the closure) |
| `ATM_HYDRO_SPLIT` | `atmosphere/AtmHydroSplit.h:72` | 1 | `1.0*` | 1.0; thermal wind, 09-24 (B.2) |
| `ATM_ICE_COLD` | `atmosphere/CloudFraction.h:198` | 1 | `1*` | 1; 08-31 |
| `ATM_ICE_LIMIT_ARRIVING` | `atmosphere/TwoCatIceScheme.h:203` | 1 | `1*` | 1; with RAIN_AREA, 09-01 |
| `ATM_MC_EVAP_LIMIT` | `atmosphere/MoistConvection.h:1398` | 1 | `1*` | 1; drift removed, 09-22 |
| `ATM_MC_S_NDIM` | `atmosphere/MoistConvection.h:1194` | 2 | `1*` | 1; 09-24 |
| `ATM_MC_T_NDIM` | `atmosphere/MoistConvection.h:1604` | 1 | `1.0*` | 1.0; 09-22 |
| `ATM_METRIC_RADIUS` | `atmosphere/cAtmosphereModel.cpp:201` | 1 | **missing** | Earth radius; 07-28 -- NOT IN BANNER |
| `ATM_MICRO_NDIM` | `atmosphere/ColumnWaterBudget.h:441` | 3 | `1.0*` | 1.0; 09-21 (read at 3 sites) |
| `ATM_NUE_GRAD` | `atmosphere/RHS_Atm_Turb.cpp:947` | 1 | `1.0*` | 1.0; 09-21 |
| `ATM_RADIAL_SHAPIRO_STRENGTH_VW` | `atmosphere/cAtmosphereModel.cpp:2105` | 1 | `0.0*` | 0.0; 09-12 |
| `ATM_RAD_TOPO` | `atmosphere/MultiLayerRadiation.h:297` | 1 | `1*` | 1; 09-24 (B.4) |
| `ATM_RH_MIN_LAT` | `atmosphere/InitValues_Atm.cpp:1290` | 1 | `1*` | 1; 08-31 |
| `ATM_RH_PROFILE` | `atmosphere/InitValues_Atm.cpp:1228` | 1 | `1*` | 1; 08-31 |
| `ATM_SATADJ_PHASE` | `atmosphere/SaturationAdjustment.h:101` | 1 | `1*` | 1; 09-09 |
| `ATM_SEAM_PERIODIC` | `atmosphere/BC_Atm.h:815` | 1 | `1*` | 1; 09-24 |
| `ATM_TROPO_INDEX_FIX` | `atmosphere/InitValues_Atm.cpp:435` | 1 | `1*` | 1; 09-12 |
| `ATM_V_MASSBAL` | `atmosphere/VelocityInitializer.h:840` | 1 | `1*` | 1; 08-28 |
| `ATM_V_MASSBAL_STRIDE` | `atmosphere/VelocityInitializer.h:776` | 1 | `1*` | 1; 09-12 |
| `HYD_BC_SECOND_ORDER` | `hydrosphere/cHydrosphereModel.h:106` | 1 | `1` | 1; 09-23 |
| `HYD_SSS_FILL` | `hydrosphere/InitValues_Hyd.cpp:1025` | 1 | `1` | 1; SSS sentinel, 09-08 |
| `HYD_T_FREEZE` | `hydrosphere/UtilsHyd.h:257` | 1 | `1` | 1; 09-05 |
| `HYD_T_FREEZE_SFC` | `hydrosphere/cHydrosphereModel.cpp:274` | 1 | `1` | 1; 09-05 |

## X -- Switched on, then reverted by measurement (1)

| knob | read at | sites | banner default | note |
|---|---|---|---|---|
| `ATM_BUOY_CONSISTENT` | `atmosphere/RHS_Atm_Turb.cpp:1320` | 1 | `0*` | reverted 09-14 (radial u 390x); HYDRO_SPLIT supersedes |

## N -- Off; measured null, refuted or superseded (22)

| knob | read at | sites | banner default | note |
|---|---|---|---|---|
| `ATM_ANELASTIC` | `atmosphere/PressureSolverAtm.h:139` | 2 | **missing** | null on Psi (-0.006 %), structural |
| `ATM_BUOY_MOIST` | `atmosphere/RHS_Atm_Turb.cpp:476` | 1 | **missing** | null unless BUOY_CONSISTENT |
| `ATM_BUOY_TREF` | `atmosphere/RHS_Atm_Turb.cpp:1318` | 1 | **missing** | measured +1.0 % on the term, null on the model |
| `ATM_CELLS_U_FROM_PSI` | `atmosphere/VelocityInitializer.h:576` | 1 | `0*` | off; not scored separately |
| `ATM_EVAP_FLUX` | `atmosphere/ThermoAtm.h:213` | 1 | `0*` | superseded by WATER_CLOSURE |
| `ATM_EVAP_SPREAD` | `atmosphere/ThermoAtm.h:163` | 1 | `0*` | null on a spun-up field (0.2 %) |
| `ATM_GRID_PRESSURE` | `atmosphere/cAtmosphereModel.h:613` | 1 | **missing** | ~1 %, sign against it |
| `ATM_HYDRO_PGF` | `atmosphere/RHS_Atm_Turb.cpp:1089` | 1 | `0*` | superseded by HYDRO_SPLIT |
| `ATM_HYDRO_PGF_RAW` | `atmosphere/RHS_Atm_Turb.cpp:1091` | 1 | `0*` | worse than HYDRO_PGF; superseded |
| `ATM_METRIC_EXACT` | `atmosphere/cAtmosphereModel.h:762` | 1 | **missing** | null on integrated quantities; undecidable |
| `ATM_METRIC_NOCURV` | `atmosphere/cAtmosphereModel.h:794` | 1 | **missing** | attribution arm, done |
| `ATM_POISSON_METRIC_FIX` | `atmosphere/PressureSolverAtm.h:137` | 3 | `0*` | retired 09-02 (+17 % pgf, 2.6x checkerboard); PARSED DIFFERENTLY at its 3 sites |
| `ATM_PRESS_LINE_SOLVE` | `atmosphere/PressureSolverAtm.h:365` | 1 | `0*` | measured through the clamp; no runaway fix |
| `ATM_PROJECT_IN_LOOP` | `atmosphere/PressureSolverAtm.h:1193` | 1 | **missing** | null at 10 and 200 sweeps |
| `ATM_PROJ_SWEEPS` | `atmosphere/PressureSolverAtm.h:1473` | 1 | **missing** | inert (-0.04 % at 10x) -- NOT IN BANNER |
| `ATM_RHIE_CHOW` | `atmosphere/PressureSolverAtm.h:348` | 1 | **missing** | null on Psi (structural) |
| `ATM_RH_CRIT_ICE` | `atmosphere/CloudFraction.h:77` | 1 | `off*` | off; global cirrus -- superseded by RH_MIN_LAT/PTOP |
| `ATM_SATADJ_FREEZE_LATENT` | `atmosphere/SaturationAdjustment.h:195` | 1 | `0*` | null (0.05 % of cloud) |
| `ATM_SFC_FLUX` | `atmosphere/RHS_Atm_Turb.cpp:1247` | 2 | **missing** | null (timescale wall) -- NOT IN BANNER |
| `ATM_TAU_PBROAD` | `atmosphere/MultiLayerRadiation.h:404` | 1 | **missing** | refuted (lapse 3 %) |
| `HYD_LINE_SOLVE` | `hydrosphere/PressureSolverHyd.h:444` | 1 | `0` | not needed (w2 slightly worse) |
| `HYD_VW_BOTTOM_ZG` | `hydrosphere/BC_Hyd.h:125` | 1 | `0` | -17 % of the profile excess; refuted as cause |

## R -- Open: repair or experiment, runs owed or decision pending (41)

| knob | read at | sites | banner default | note |
|---|---|---|---|---|
| `ATM_CONV_ADJ` | `atmosphere/ConvectiveAdjustment.h:70` | 1 | **missing** | off; palliative -- NOT IN BANNER |
| `ATM_DAMP_Q_HORIZ` | `atmosphere/cAtmosphereModel.cpp:1635` | 1 | `1*` | step A; qh600 scored |
| `ATM_DAMP_Q_MASS` | `atmosphere/cAtmosphereModel.cpp:1617` | 1 | `0*` | forced by WATER_CLOSURE |
| `ATM_DAMP_Q_VERT` | `atmosphere/cAtmosphereModel.cpp:1626` | 1 | `1*` | step A |
| `ATM_DAMP_T_HORIZ` | `atmosphere/cAtmosphereModel.cpp:1598` | 1 | `1*` | qth_on today |
| `ATM_DAMP_T_VERT` | `atmosphere/cAtmosphereModel.cpp:1590` | 1 | `1*` | null alone (qvt_on) |
| `ATM_ICE_LIMITERS` | `atmosphere/ThreeCatIceScheme.h:257` | 1 | `0*` | ThreeCat only |
| `ATM_ICE_RAW_FLUX` | `atmosphere/ThreeCatIceScheme.h:231` | 1 | `0*` | ThreeCat only |
| `ATM_LAND_BUCKET` | `atmosphere/ThermoAtm.h:288` | 1 | `0*` | 20-iter evidence only |
| `ATM_LENGTH_NDIM` | `atmosphere/cAtmosphereModel.h:849` | 1 | **missing** | the 40x L_atm defect (item 4 pending) -- NOT IN BANNER |
| `ATM_MC_ALF1` | `atmosphere/MoistConvection.h:32` | 1 | `0.05*` | alf_tie/alf_mid scored |
| `ATM_MC_BASE_SAT` | `atmosphere/MoistConvection.h:908` | 1 | `0*` | =2 in gpa_* today |
| `ATM_MC_ENTR` | `atmosphere/MoistConvection.h:42` | 1 | `2.0e-3*` | B.10b |
| `ATM_MC_GP_AREA` | `atmosphere/MoistConvection.h:1001` | 1 | `0*` | gpa_* today |
| `ATM_MC_QVD` | `atmosphere/MoistConvection.h:177` | 1 | `0*` | B.10c |
| `ATM_MC_SGZ` | `atmosphere/MoistConvection.h:166` | 1 | `0*` | B.10b |
| `ATM_ONECAT_CLOUD_LIMIT` | `atmosphere/OneCatIceScheme.h:301` | 1 | `0*` | new, byte check tomorrow |
| `ATM_OROG_Q_MASS` | `atmosphere/cAtmosphereModel.cpp:2139` | 1 | `0*` | arm owed |
| `ATM_RAD_EQUIL` | `atmosphere/MultiLayerRadiation.h:183` | 1 | `0*` | pair with SW_INSOL; needs a prognostic T |
| `ATM_RK_SCALAR_SYNC` | `atmosphere/cAtmosphereModel.cpp:1982` | 1 | `0*` | forced by WATER_CLOSURE |
| `ATM_SATADJ_FADE` | `atmosphere/SaturationAdjustment.h:158` | 1 | `0*` | mode 2 if ever flipped; forced by closure |
| `ATM_SEAM_Q_CONSERVE` | `atmosphere/BC_Atm.h:831` | 1 | `0*` | new, byte check tomorrow |
| `ATM_SURF_DRAG_CONSISTENT` | `atmosphere/RHS_Atm_Turb.cpp:1407` | 2 | `0.0*` | B.9, arm postponed |
| `ATM_SW_INSOL` | `atmosphere/MultiLayerRadiation.h:110` | 1 | `0*` | pair with RAD_EQUIL |
| `ATM_TEQ_SKIN_ONLY` | `atmosphere/cAtmosphereModel.cpp:822` | 1 | **missing** | instrument branch, must NOT be flipped -- NOT IN BANNER |
| `ATM_TURB_SIN_FLOOR` | `atmosphere/TurbulenceAtm.h:29` | 1 | `0*` | new, byte check tomorrow |
| `ATM_TW_BALANCE` | `atmosphere/VelocityInitializer.h:1013` | 1 | `0.0*` | off; the only mid-lat jet IC |
| `ATM_WATER_CLOSURE` | `atmosphere/SaturationAdjustment.h:160` | 4 | `0*` | reverted 09-25; works with filter off (qh600) |
| `HYD_A_H` | `hydrosphere/RHS_Hyd_Turb.cpp:762` | 2 | `0` | Laplacian; biharmonic preferred |
| `HYD_A_H_BIHARM` | `hydrosphere/HydHorizViscosity.h:47` | 2 | `0` | metric-branch stack |
| `HYD_BAROCLINIC_PGF` | `hydrosphere/RHS_Hyd_Turb.cpp:969` | 2 | `0.0` | pair with PHYDRO_SALT; HYDRO_SPLIT may supersede |
| `HYD_BC_DRAG` | `hydrosphere/PressureSolverHyd.h:1067` | 1 | **missing** | numerical, not physical -- NOT IN BANNER |
| `HYD_BUOY_CONSISTENT` | `hydrosphere/HydBuoyancy.h:51` | 2 | `0.0` | radial runaway; HYDRO_SPLIT |
| `HYD_DEEP_DRAG` | `hydrosphere/RHS_Hyd_Turb.cpp:1248` | 2 | `0` | B.7, arm owed |
| `HYD_HYDRO_SPLIT` | `hydrosphere/HydBuoyancy.h:108` | 2 | `0.0` | new, runs tomorrow |
| `HYD_METRIC_RADIUS` | `hydrosphere/cHydrosphereModel.cpp:475` | 1 | `0` | blocked by the profile item |
| `HYD_METRIC_SIN_FLOOR` | `hydrosphere/cHydrosphereModel.h:357` | 1 | `0.4` | new, sweep tomorrow |
| `HYD_NUE_GRAD` | `hydrosphere/RHS_Hyd_Turb.cpp:654` | 1 | **missing** | off; would inherit the broken metric -- NOT IN BANNER |
| `HYD_PHYDRO_SALT` | `hydrosphere/ThermoHyd.h:101` | 1 | `0` | pair |
| `HYD_RUN_NEUMANN` | `hydrosphere/PressureSolverHyd.h:39` | 1 | `0` | pair with METRIC_RADIUS |
| `HYD_SFC_FLUX` | `hydrosphere/RHS_Hyd_Turb.cpp:1111` | 2 | `0` | unmeasured |

## P -- Numerical parameter (a value, not a switch) (25)

| knob | read at | sites | banner default | note |
|---|---|---|---|---|
| `ATM_CO2_BAND` | `atmosphere/MultiLayerRadiation.h:546` | 1 | **missing** | 0.17 -- NOT IN BANNER |
| `ATM_CONV_ADJ_LAPSE` | `atmosphere/ConvectiveAdjustment.h:84` | 1 | **missing** | sub-parameter of CONV_ADJ |
| `ATM_CONV_ADJ_PASSES` | `atmosphere/ConvectiveAdjustment.h:88` | 1 | **missing** | sub-parameter of CONV_ADJ |
| `ATM_EPS_DRY` | `atmosphere/MultiLayerRadiation.h:390` | 1 | **missing** | 0.684 -- NOT IN BANNER |
| `ATM_GRID_BETA` | `atmosphere/cAtmosphereModel.h:625` | 1 | **missing** | 3.988 (with GRID_PRESSURE) |
| `ATM_GRID_PTOP` | `atmosphere/cAtmosphereModel.h:618` | 1 | **missing** | 0.08538 (with GRID_PRESSURE) |
| `ATM_HADLEY_SL` | `atmosphere/VelocityInitializer.h:105` | 2 | `4.0N/3.0S*` | IC 4.0N/3.0S |
| `ATM_METRIC_SIN_FLOOR` | `atmosphere/cAtmosphereModel.h:394` | 1 | `0.26*` | 0.26 since 09-10 |
| `ATM_PDYN_CAP` | `atmosphere/PressureSolverAtm.h:297` | 2 | `2.0*` | 2.0 |
| `ATM_PDYN_CEILING` | `atmosphere/PressureSolverAtm.h:840` | 1 | `3.0*` | 3.0; release measured harmless |
| `ATM_POLAR_CELL_SHEAR` | `atmosphere/VelocityInitializer.h:66` | 2 | `0.1*` | IC 0.1 |
| `ATM_PRESS_SWEEPS` | `atmosphere/PressureSolverAtm.h:355` | 1 | **missing** | 1x -- NOT IN BANNER |
| `ATM_PSI_SHAPE` | `atmosphere/VelocityInitializer.h:442` | 1 | `1*` | IC shape 1 |
| `ATM_QC_CRIT` | `atmosphere/IceSchemeCommon.h:35` | 1 | `0.05*` | 0.05 g/kg since 08-31 |
| `ATM_RADIAL_SHAPIRO_STRENGTH` | `atmosphere/cAtmosphereModel.cpp:2080` | 1 | `1.0*` | 1.0 (u only) |
| `ATM_RAIN_AREA` | `atmosphere/IceSchemeCommon.h:81` | 1 | `0.10*` | 0.10 since 09-01 (fitted) |
| `ATM_RH_CRIT` | `atmosphere/CloudFraction.h:52` | 2 | `0.30*` | 0.30 since 08-31 |
| `ATM_RH_MIN` | `atmosphere/InitValues_Atm.cpp:1258` | 1 | `0.65*` | 0.65 |
| `ATM_RH_MIN_PTOP` | `atmosphere/InitValues_Atm.cpp:1336` | 1 | `482*` | 482 hPa (fitted) |
| `ATM_TW_BALANCE_V` | `atmosphere/VelocityInitializer.h:1032` | 1 | `0*` | sub-option of TW_BALANCE |
| `ATM_TW_LATMIN` | `atmosphere/VelocityInitializer.h:1018` | 1 | `15*` | sub-parameter of TW_BALANCE |
| `ATM_TW_WMAX` | `atmosphere/VelocityInitializer.h:1023` | 1 | **missing** | sub-parameter of TW_BALANCE -- NOT IN BANNER |
| `ATM_T_FLOOR` | `atmosphere/InitValues_Atm.cpp:848` | 1 | `216.65*` | 216.65 K since 08-31 |
| `HYD_LINE_FOLD` | `hydrosphere/PressureSolverHyd.h:624` | 1 | **missing** | sub-option of LINE_SOLVE |
| `HYD_LINE_GAUGE` | `hydrosphere/PressureSolverHyd.h:639` | 1 | **missing** | sub-option of LINE_SOLVE |

## S -- Structural option at its shipped value (4)

| knob | read at | sites | banner default | note |
|---|---|---|---|---|
| `ATM_RADIATION_MODE` | `atmosphere/cAtmosphereModel.cpp:740` | 1 | **missing** | 5 (radiation diagnostic) -- NOT IN BANNER |
| `ATOM_CORIOLIS_NONTRAD` | `lib/Utils.h:68` | 1 | **missing** | 0 (lib) -- NOT IN BANNER |
| `ATOM_METRIC_CURVATURE` | `lib/Utils.h:106` | 1 | **missing** | 0 (lib) -- NOT IN BANNER |
| `ATOM_METRIC_DIVERGENCE` | `lib/Utils.h:112` | 1 | **missing** | 0 (lib) -- NOT IN BANNER |

## I -- Output / run infrastructure (4)

| knob | read at | sites | banner default | note |
|---|---|---|---|---|
| `ATM_LONGAL_J` | `atmosphere/UtilsAtm.h:360` | 1 | **missing** | output slice latitude |
| `ATM_METRIC_STRICT` | `atmosphere/cAtmosphereModel.cpp:253` | 1 | **missing** | abort on metric check |
| `ATM_RESTART_STRIDE` | `atmosphere/cAtmosphereModel.cpp:2285` | 1 | `100*` | restart cadence |
| `ATM_VTK_STRIDE` | `atmosphere/cAtmosphereModel.cpp:1877` | 1 | `5*` | VTK cadence |

## D -- Diagnostic, print or dump only (18)

| knob | read at | sites | banner default | note |
|---|---|---|---|---|
| `ATM_CELL_ROT_DIAG` | `atmosphere/MinMax_Atm.cpp:250` | 1 | `0*` | print/dump only |
| `ATM_CLOUD_INIT_DIAG` | `atmosphere/InitValues_Atm.cpp:1596` | 2 | **missing** | print/dump only |
| `ATM_CWB_BANDS` | `atmosphere/ColumnWaterBudget.h:163` | 1 | **missing** | print/dump only |
| `ATM_CWB_DIAG` | `atmosphere/ColumnWaterBudget.h:156` | 1 | **missing** | print/dump only |
| `ATM_CWP_CENSUS` | `atmosphere/MultiLayerRadiation.h:125` | 1 | **missing** | print/dump only |
| `ATM_MC_CAP_DIAG` | `atmosphere/MoistConvection.h:1547` | 1 | **missing** | print/dump only |
| `ATM_MC_DIAG` | `atmosphere/MoistConvection.h:2018` | 1 | **missing** | print/dump only |
| `ATM_METRIC_CHECK` | `atmosphere/cAtmosphereModel.cpp:292` | 1 | **missing** | print/dump only |
| `ATM_MFC_DIAG` | `atmosphere/ThermoAtm.h:1194` | 1 | **missing** | print/dump only |
| `ATM_PROJ_CONSISTENCY` | `atmosphere/PressureSolverAtm.h:1313` | 1 | **missing** | print/dump only |
| `ATM_PSI_PROJ_DUMP` | `atmosphere/cAtmosphereModel.cpp:600` | 1 | **missing** | print/dump only |
| `ATM_RAD_COLDIAG` | `atmosphere/MultiLayerRadiation.h:177` | 1 | **missing** | print/dump only |
| `ATM_SATADJ_DIAG` | `atmosphere/SaturationAdjustment.h:200` | 1 | **missing** | print/dump only |
| `ATM_SR_DIAG` | `atmosphere/TwoCatIceScheme.h:211` | 1 | **missing** | print/dump only |
| `ATM_SS_DIAG` | `atmosphere/ThreeCatIceScheme.h:263` | 1 | **missing** | print/dump only |
| `ATM_T0_ATTRIB` | `atmosphere/cAtmosphereModel.cpp:1216` | 1 | **missing** | print/dump only |
| `ATM_UBUD_BALANCE` | `atmosphere/Results_Atm.cpp:140` | 1 | **missing** | print/dump only |
| `HYD_KE_SPLIT` | `hydrosphere/cHydrosphereModel.cpp:727` | 1 | **missing** | print/dump only |
