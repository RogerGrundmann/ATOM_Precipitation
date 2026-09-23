#pragma once
/*
 * ATM_HYDRO_SPLIT=<strength> -- THE HYDROSTATIC / NON-HYDROSTATIC PRESSURE SPLIT.
 * Default 0.0 = OFF and bit-identical: nothing below is read, and nothing is allocated.
 *
 * WHY. B.2: this model has no thermal wind. Two routes were tried and each failed on its own:
 *   - ATM_HYDRO_PGF puts grad_h(p_stat) into rhs_v/rhs_w. It balances Coriolis (jet-core
 *     residual 0.9999 -> 0.25), but p_stat is DIAGNOSED from T with p_sl ∝ T_surface and does
 *     not answer the flow, so it is a forcing and not a pressure.
 *   - ATM_BUOY_CONSISTENT gives the buoyancy its correct size and lets the projection pressure
 *     answer it. The projection answered only 24 % (non-adjoint div/grad, a radially dominated
 *     operator, the ceiling clamp), and the other 76 % became a net radial force: `u` ran to
 *     390x what the model's own Psi requires, and the knob was reverted on 2026-09-14.
 * Hydrostatic balance does not need an elliptic solve. It is a COLUMN INTEGRAL. So split
 *
 *     p_dyn_total = p_hb + p_nh,     exp_rm * d(p_hb)/d(rad.z) = b   (exactly, by construction)
 *
 * and then
 *     rhs_u :  b - dp_hb/dr - dp_nh/dr  =  - dp_nh/dr       -- b cancels; nothing left to run away
 *     rhs_v :  - inv_rm       * dp_hb/dthe - dp_nh/dthe     -- grad_h(p_hb) IS the thermal wind
 *     rhs_w :  - inv_rmsinthe * dp_hb/dphi - dp_nh/dphi
 * with `p_dyn` from the Poisson solve now playing p_nh: `aux` = rhs + grad(p_dyn) carries the
 * divergence of the hydrostatic horizontal force and no buoyancy, so the projection only has to
 * remove what the hydrostatic part leaves divergent. This is how hydrostatic GCMs and ocean
 * models obtain their pressure.
 *
 * THE BUOYANCY IS THE CONSISTENT ONE, times the strength: s*ramp*buoyancy*(g*L_atm/u_0^2) *
 * (tn - t_ref)/t_ref -- exactly ATM_BUOY_CONSISTENT's coefficient, because that is the one whose
 * balance was measured (band p05 residual 0.996 -> 0.336). The shipped g*dt/u_0 carries `dt`
 * twice and would give a thermal wind 5e5 too weak. Temperature only (ATM_BUOY_MOIST is not
 * followed here).
 *
 * THE INTEGRAL IS IN PHYSICAL HEIGHT, dp_nd/dz = (g/u_0^2) dT/T_ref per METRE, because the
 * horizontal gradient is taken with inv_rm / inv_rmsinthe, whose length unit is
 * metricShellLength() -- the same as Coriolis's. It is NOT integrated with the RHS's radial
 * operator exp_rm*d/d(rad.z): with the legacy exp_rm that operator is 2-50x off the physical
 * one (the 23x-spread defect), and a version built on it measured 14x too weak. Since rhs_u no
 * longer contains b, the radial operator does not enter the split at all. Integrated UPWARD by the trapezoid rule in rad.z from p_hb = 0 at
 * i = 0 -- NOT from the lid; see the note at the integration loop for why the lid gauge was
 * wrong. The integral is of the ANOMALY against t_ref_level[i], so its level mean is ~0 and a uniform offset per level has no horizontal gradient anyway.
 *
 * TIME LEVEL: built from `tn`, once per iteration before the RK4 sweep. The RK4 here is POINTWISE
 * and writes `t` in place during the sweep, so a neighbour's `t` is at a sweep-order-dependent
 * stage; `tn` is constant through the sweep, which makes this deterministic. The price is that
 * the hydrostatic pressure lags the stage temperature by < 1 dt -- the radial balance is exact
 * for the field it was built from, and rhs_u no longer contains b at all.
 *
 * TERRAIN. Below i_topography the column is rock and carries b = 0, so p_hb is 0 through the rock
 * and a land column is referenced at its OWN ground. The horizontal difference uses only AIR
 * neighbours: centred if both are air, one-sided against the cell itself if one is, zero if
 * neither -- so a rock value is never differenced. Level i is a constant height in this grid, so
 * a difference at constant i is at constant height. What remains at a plateau edge is that the
 * neighbouring ocean column has accumulated b over the plateau's height and the plateau column
 * has not: a spurious gradient of order INT_0^z_ground grad(b) dz. It is bounded by the surface
 * anomaly times the terrain height and is NOT measured; it is the price of having no
 * prognostic barotropic mode to set each column's constant.
 *
 * Storage is namespace-scope, NOT a cAtmosphereModel member: adding a member moves
 * sizeof(cAtmosphereModel), which is this tree's stack-canary hazard.
 */
#include <vector>
#include <cstdlib>

namespace AtmHydroSplit {

inline double strength(){
    static const double s = [](){ const char* e = getenv("ATM_HYDRO_SPLIT"); return e ? atof(e) : 0.0; }();
    return s;
}
inline bool enabled(){ return strength() != 0.0; }

inline std::vector<double> p_hb;            // im*jm*km, (i*jm + j)*km + k
inline int n_i = 0, n_j = 0, n_k = 0;

inline double at(int i, int j, int k){ return p_hb[((size_t)i * n_j + j) * n_k + k]; }

} // namespace AtmHydroSplit
