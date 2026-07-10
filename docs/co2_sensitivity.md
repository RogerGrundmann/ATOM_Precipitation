# ATOM CO₂ Sensitivity — Summary of Findings

*0 Ma (present-day) configuration. Distilled from the 2026-07 CO₂ investigation.
Restart equilibria and run logs live under `python/output_0Ma_co2x{1,2}_*/`.*

## Headline result: a two-mode answer

The model answers "what does doubling CO₂ do?" **differently depending on which
radiative energy sink drives the column**, and no single mode currently gives
both the right temperature response and the right hydrological response:

| Radiation mode | Energy sink | Global ΔT / doubling | Precip response (CC sign) | Snow fraction |
|---|---|---|---|---|
| **Mode 5** (production) | Newtonian relaxation of `t` toward `t_eq = Scotese + CO₂ perturbation` | **+0.91 K** ✓ | **−37 %** ✗ (wrong sign) | ~17 % ✓ |
| **Mode 2** (experimental) | direct σT⁴ radiative heating from MultiLayerRadiation (MLR) | +0.16–0.19 K ✗ (weak) | **+12–24 %** ✓ (correct, CC-positive) | ~11 % ✓ |

- **Mode 5** reproduces a realistic warming (+0.91 K/doubling, in the low end of
  the IPCC range) but the hydrological cycle spins *down* with warming — the
  wrong sign for Clausius–Clapeyron.
- **Mode 2** reproduces the correct positive precipitation response
  (+12–24 %, robust across four independent equilibria A/B/C/D and an all-fix
  run) but its global warming is structurally weak (~0.17 K).

Both results are **settled and robust** — each was stress-tested against every
lever we could find (below). The tension between them is the standing open
problem.

## Why mode 5 gets precip wrong

Mode 5's column energy sink is Newtonian relaxation, `t ← (1−ω)·t + ω·t_eq`
with ω = 0.20/iter. This cooling depends only on the **departure** `t − t_eq`,
**not on absolute temperature**. Real radiative cooling (σT⁴) *rises* with
warming and is exactly what forces the observed ~2–3 %/K energetic increase in
global precipitation. Newtonian relaxation has no such T-dependence: warming the
whole column (both `t` and `t_eq` rise together) does not increase the cooling,
so latent heating — and therefore precipitation — cannot increase. Tuning ω only
changes how tightly `t` tracks `t_eq`; it cannot manufacture the missing σT⁴
response. The wrong-sign precip is baked into using `t_eq` relaxation as the sink.

Corroborating evidence: precip (~6 mm/d) is ~10× global evaporation (~0.57 mm/d),
i.e. precip does **not** close the surface moisture budget — it is slaved to the
imposed energy relaxation. Ascent (radial velocity `u`) and the meridional
streamfunction are *identical* between 1× and 2× CO₂; column RH stays ~saturated;
model surface water vapour rises +10.8 %/K (CC-correct). The moisture supply is
fine — the precip is energy-constrained, and mode 5 constrains it the wrong way.

## Why mode 2 gets ΔT weak — the realization ceiling

Mode 2 uses the σT⁴ MLR radiative equilibrium (RE) as the sink, nudged with rate
ω_rad. The single-column RE itself gives a healthy **~1–2 K/doubling**, but the
full 3-D model realizes only **~13–15 %** of it (~0.17 K). This ceiling was
proven robust by ruling out, one at a time:

- **Cloud masking** of the CO₂ band — moved the forcing to a cloud-transparent
  L_down increment (mode2C). ΔT unchanged (+0.18 vs +0.19 K).
- **ω_rad strength** — tripled 0.05→0.15 (mode2B). ΔT barely moved (+0.15→+0.19 K).
- **Forcing magnitude** — set to the canonical 3.7 W/m²/doubling. No change.
- **Ice-albedo feedback** — replaced the fixed albedo parabola with a live
  T-dependent surface albedo (mode2D, committed `5862904`). Arctic amplifies
  (+0.27 K local, ~1.4×) but global ΔT is flat: a warm low-ice world has almost
  no marginal ice to retreat, and annual-mean polar insolation is low.
- **Missing convective adjustment** — added a Manabe–Wetherald moist-adiabat
  relaxation (all-fix run). It halved a spurious free-troposphere cooling
  (−0.15→−0.056 K) but left the **surface** ΔT unchanged (+0.16 K).

**Diagnosis (vertical ΔT profile, mode2D iter-2200):** the grey-RE CO₂ response
is dominated by **mid-tropospheric cooling** (−6.5 W/m² ≈ −1.2 K at i≈8), not
surface warming. The boundary layer warms (+0.1 K) but the free troposphere
cools; the cooled free troposphere emits less downward LW, which cancels the
surface CO₂ forcing. The realization fraction is ~13 % *everywhere* — so raising
ω_rad is counterproductive (more realization = more mid-trop cooling). The limit
is the dynamical core mixing away the column radiative signal, not any radiation
formulation detail. See `project_convective_moisture_cc_coupling` and
`project_multilayer_radiation` in the memory store for the full derivation.

## The precipitation payoff and the snow prerequisite

Getting a *correct* precip sign out of any mode required first fixing a
microphysics bias, then choosing a T-dependent energy sink:

1. **Dead convective moisture seed** — `MoistConvection.h::surfaceFluxPerturbations`
   was driven by `Q_latent_2D`/`Q_sensible_2D`, which are never assigned (always
   0), so the convective δq/δT perturbation was dead and CC-blind. Fixed to use
   the T-responsive Dalton evaporation and a Clausius–Clapeyron `cc_factor` seed
   (committed `149bd0a`). Convection is now live and CC-coupled (warm-run peak
   convective precip 2.49 vs 0.78 mm/d control).

2. **Snow over-production** — snow was 45–58 % of total precip (Earth ~5–10 %),
   dominated by riming monopolising the cloud-water limiter. Fixed by reducing
   snow-side riming and accelerating autoconversion so freed cloud water rains
   rather than accretes (committed `c4d1bd1`). Snow fell to ~16 % with total
   precip preserved. This was the *prerequisite*: once snow is realistic and the
   precip is rain-dominated, the healthy rain signal (+15–18 % under warming)
   shows through instead of being swamped by a collapsing snow phase.

3. **σT⁴ sink (mode 2)** — with snow fixed, the mode-2 radiative equilibrium
   gives the correct CC-positive precip response (+12–24 %), and snow stays
   controlled (~11 %). This *is* the payoff, and it is secure.

Later cloud/ice-scheme cleanups: shared vapour→ice→snow throttle (`43f7c55`),
balanced cloud-condensate cap the radiation sees (`e0ea8fd`), shared ice-scheme
helpers header + microphysics fixes (`8b20450`).

## What is committed vs. experimental

- **Production source = mode 5** (`radiation_mode = 5`,
  `cAtmosphereModel.cpp:418`). All microphysics fixes, the ice-albedo feedback,
  the convergence monitor, and the diagnostics dumps are committed on `main`.
- **Mode 2 is an experiment, reverted out of source.** Re-enabling it means
  setting `radiation_mode = 2`, `co2_log_k = 5.0` (cloud-transparent L_down,
  ≈3.7 W/m²/doubling), `ω_rad = 0.15`, plus the convective adjustment and
  cloud-cap fixes. All mode-2 equilibria are preserved on disk:
  `output_0Ma_co2x{1,2}_mode2{A,B,C,D}/`, `_allfix/`, `_icefb/` (`atm_restart_*.bin`).

## Open problems / recommendations

- **The realization ceiling is the deep limiter.** It is a property of the
  dynamical core washing out ~85 % of the column RE signal — not a radiation bug.
  Breaking it is structural work (the earlier architecture note points to an
  anelastic Boussinesq ρ₀(z) upgrade, ~1–2 months). Radiation/feedback tweaks are
  exhausted.
- **Mode-2 absolute precip is ~2–4× too high** (~9–12 mm/d vs NASA ~2.7): the
  grey σT⁴ sink over-cools and over-condenses. Sign is right, amount is a
  grey-climatology bias.
- **Paleo Ma=100/200** have a separate, severe pre-existing over-precip
  (~237 mm/d) and non-convergence — a distinct real problem, not addressed here.
