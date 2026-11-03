// SPDX-License-Identifier: Apache-2.0
// difftest shim for unit `blend.fillet_corner_p3`, part 1: SmFilletCorner::Create, the corner
// classification of upstream v1.0.1's source/SMLib/src/SmFilletCorner.cpp, on a real SmBrep
// with an SmFilletExecutive whose solvers carry the cross-section surface generators.
//
// Not linked: the SMLib fork the difftest is pinned to classifies corners for its generic
// blend. It has no SmLinearCrossSectionFSG, SmCircularCrossSectionFSG (tangency-labs/smlib
// 3aae1c4) or SmBlendCurveCrossSectionFSG (0bd3f4d), and no 1x1 or 2x2 (543d25a), Nx1, open or
// degenerate (53c47b8) or Nx2 (74e9beb) corner. The entry point keeps its signature
// (crates/solidcore-ffi/src/units/blend_fillet_corner_p3.rs) and returns SC_ERR_UNLINKED; the
// v1.0.1 shim is in this file's history.
#include "sc_shim.h"

SC_EXPORT int32_t sc_blend_fillet_corner_p3_create(const char* path, const double* in, int64_t n_in, double* out,
                                                   int64_t cap, int64_t* out_len)
{
    (void)path, (void)in, (void)n_in, (void)out, (void)cap, (void)out_len;
    return SC_ERR_UNLINKED;
}
