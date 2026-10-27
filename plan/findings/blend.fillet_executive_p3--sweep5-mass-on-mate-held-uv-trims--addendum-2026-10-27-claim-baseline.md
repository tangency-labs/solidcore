# Claim baseline: sweep5 known mismatches on fbd9b3f74 (2026-10-27)

Updates: plan/findings/blend.fillet_executive_p3--sweep5-mass-on-mate-held-uv-trims.md

Claim for the three known mismatches that keep the unit partial (draft PR, branch
`fix/fillet-exec-p3-sweep5-known-mismatches`). Baseline on main fbd9b3f74, 2026-10-27:
PARTIAL, cases=592, passed=587, failed=0, max_rel_err=1e0, with_values=587,
no_values=0, vacuous=0, rejected=2 (native crash traps UB-209 / UB-261, pre-existing,
do not block green), known_mismatch=3.

Scope correction from the claim as filed: the lazy surface-of-extrusion knots fix is
already merged upstream (4e86db739, 2026-10-27;
plan/findings/blend.fillet_executive_p3--sweep5-lazy-extrusion-knots--addendum-2026-10-27-resolved.md).
The remaining work is:

- `sweep5-letter_g-e10-r1.0` — mass only: 855.6511940507597 vs 855.6495364737889;
  topology and geometry agree; the fillet face's UV trims are shifted by 5.7e-14 in u
  and the adaptive area integration amplifies it. Localising and fixing the trim
  source is the in-scope work here.
- `sweep5-ext_w_hole-e15-r1.0` — the to-Rust import drops the planes' NURBs; the lazy
  remake is 1 ulp off, so the fillet vertex UV and the rail natural interval differ by
  1 ulp and the stitch keeps the other rail edge
  (plan/findings/blend.fillet_executive_p3--sweep5-ext-w-hole-e15-face-domain--addendum-2026-10-27-plane-nurb-import.md).
  The deviation is owned by the open upstream PR #1318 (fillet runner's own-reader
  input) and the missing lazy-NURB cache in #1300 / #1349.
- `sweep5-ext_w_hole-e14-r1.0` — same import-NURB rounding family; one near-zero
  vertex coordinate (C++ 1.4453956777755315e-15, Rust 3.6012056774890033e-16),
  everything else agrees
  (plan/findings/blend.fillet_executive_p3--sweep5-ext-w-hole-e14-near-zero-vertex.md).

Coordination: #1300 / #1318 / #1349 work the two deviations behind e14/e15. As each
lands, re-run the case and remove its entry from
crates/difftest/known_mismatches/blend.fillet_executive_p3.toml; when all three are
clear, set plan/units/blend.fillet_executive_p3.toml to green. No line overlap is
expected with the open PRs; if #1318 lands while e10 is being localised, rebase onto it.
