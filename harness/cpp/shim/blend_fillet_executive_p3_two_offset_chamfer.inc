// SPDX-License-Identifier: Apache-2.0
// tangency-labs/solidcore#880: an unequal (two-offset) chamfer of the rim a flat cut
// leaves on a placed cylinder, through upstream v1.0.1's SmConstantRadiusFS with
// pdOptSecondRadius and an SmLinearCrossSectionFSG.
//
// Not linked: the SMLib fork the difftest is pinned to builds chamfers with its generic blend
// and has no SmLinearCrossSectionFSG (tangency-labs/smlib 3aae1c4). The entry point keeps its
// signature (crates/solidcore-ffi/src/units/blend_fillet_executive_p3.rs) and returns
// SC_ERR_UNLINKED; the v1.0.1 shim is in this file's history.

SC_EXPORT int32_t sc_blend_fillet_executive_p3_two_offset_chamfer(
    const double* input, int64_t n, double* out, int64_t cap, int64_t* out_len) {
    (void)input, (void)n, (void)out, (void)cap, (void)out_len;
    return SC_ERR_UNLINKED;
}
