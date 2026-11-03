// SPDX-License-Identifier: Apache-2.0
// Included at the end of blend_fillet_geom_p1.cpp: the difftest shim of unit
// `blend.fillet_corner_p1`, part 2 (the SmFilletCorner methods), each on a corner of one of
// upstream v1.0.1's derived corner classes. Part 1 (the file-static helpers, on a real
// SmBrep) is blend_fillet_corner_p1.cpp.
//
// Not linked: the SMLib fork the difftest is pinned to has none of those corner classes (the
// 1x1 and 2x2 corners, tangency-labs/smlib 543d25a; the Nx1, open and degenerate ones, 53c47b8;
// the Nx2 one, 74e9beb) and no cross-section surface generators (3aae1c4, 0bd3f4d). The entry
// point keeps its signature (crates/solidcore-ffi/src/units/blend_fillet_corner_p1.rs) and
// returns SC_ERR_UNLINKED; the v1.0.1 shim is in this file's history.

SC_EXPORT int32_t sc_blend_fillet_corner_p1_corner(int32_t op, const double* in, int64_t n_in, double* out,
                                                   int64_t cap, int64_t* out_len)
{
    (void)op, (void)in, (void)n_in, (void)out, (void)cap, (void)out_len;
    return SC_ERR_UNLINKED;
}
