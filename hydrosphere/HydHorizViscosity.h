#pragma once

#include <vector>
#include <cstdlib>
#include <cstddef>

// ======================================================================================
// Scratch state for the HORIZONTAL biharmonic viscosity, HYD_A_H_BIHARM.
//
// A biharmonic operator is not computable from one cell's own derivatives: nabla^4 = the
// Laplacian applied twice, so the outer application needs its NEIGHBOURS' Laplacians. The
// RHS here is evaluated cell by cell, so the inner application has to be a separate pass
// over the grid, and its result has to live somewhere both translation units can see.
//
// It lives HERE, in inline variables, rather than as Array members of cHydrosphereModel,
// deliberately: adding members moves sizeof(cHydrosphereModel), which is the stack-canary
// hazard this tree has already been bitten by once (see CLAUDE.md, "The build hazard").
// Two doubles per cell is ~43 MB at 41x181x361, allocated only when the knob is on.
//
// The inner pass is evaluated on the TIME-LEVEL-n fields vn/wn, not on the running v/w.
// That is a real choice and it is the conservative one: this RK4 is POINTWISE -- all four
// stages are taken at one cell before the sweep moves on -- so a neighbour's v is at a
// stage that depends on the sweep order, and an operator built from it would be
// thread-order dependent in a model already documented as not bit-reproducible under
// OpenMP. vn/wn are written once per iteration by storeIntermediateData3D and are constant
// through the sweep, so this operator is deterministic and sweep-order independent. The
// price is that the biharmonic is lagged by one iteration, which for a small linear
// dissipation is standard.
// ======================================================================================
namespace HydHorizVisc {

// H(v,w): the ANGULAR part of the horizontal vector Laplacian -- the 1/r^2 metric is NOT
// folded in, so the buffer is a pure angular operator and each application's metric factor
// is supplied once, in the RHS coefficient. Applying it twice gives r^4 in the denominator.
inline std::vector<double>        lap_v;
inline std::vector<double>        lap_w;
inline std::vector<unsigned char> ok;      // 1 where the inner pass actually ran
inline int s_jm = 0, s_km = 0;

inline std::size_t idx(int i, int j, int k){
    return (static_cast<std::size_t>(i) * s_jm + j) * s_km + k;
}

// B in m^4/s. Default 0.0 = OFF and bit-identical.
inline double biharm_strength(){
    static const double b = [](){
        const char* e = getenv("HYD_A_H_BIHARM"); return e ? atof(e) : 0.0; }();
    return b;
}

}  // namespace HydHorizVisc
