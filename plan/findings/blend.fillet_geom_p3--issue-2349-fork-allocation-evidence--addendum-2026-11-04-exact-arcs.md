# blend.fillet_geom_p3: exact-arc options explain the sampled geometry differences (2026-11-04)

Updates: `plan/findings/blend.fillet_geom_p3--issue-2349-fork-allocation-evidence--addendum-2026-11-04-diagnosis.md`.

Units: `blend.fillet_geom_p3`, `blend.fillet_executive_p3`.
PR: https://github.com/tangency-labs/solidcore/pull/2381

## Resolved sampled-geometry discrepancy

The 54 generated surface-coordinate differences and 11 corner-edge sample
differences came from comparing approximate circular sections with exact ones.
They are not evidence of a missing numerical fix in the exact-arc path.

The independent WSL diagnostic reproduced the default Rust values, then enabled
only the existing `fillet_exact_arc_sections` flag. All 18 native surface points
(54 coordinates, two split surfaces) then agreed with fork a18bbe523b9c at the
unchanged relative tolerance 1e-9. Maximum measured relative error in those
coordinates was 3.720817356480505e-14. The all-fixes-except-geometric-G2 profile
agreed too (3.698603521516443e-14). The diagnostic constructed the operand in Rust;
it did not import native geometry or feed native output into the construction.

The persistent regression
`closed_angle_exact_arc_geometry_samples_match_the_fork` repeats both matching
profiles and also compares all 19 corner-edge interval/sample observations.
It reaches actual solver geometry and corner-edge calculation. Native expected
samples, with their fork revision and case IDs, live in
`crates/solidcore/tests/data/issue_2349_fork_surface_samples.txt` and
`issue_2349_fork_corner_samples.txt`. They are comparison operands only.

The first exact-arc point is
(-3.968755872964239, 2.0020378837388937, 9.624514441544004), versus native
(-3.9687558729642345, 2.002037883738891, 9.624514441544006).
The approximate-arc default remains the earlier
(-3.9687368198336785, 2.0044118846065655, 9.631940486795703).
`fillet_executive::p2` selects `!fillet_exact_arc_sections` for the circular
generator's approximate-arcs argument. This isolates the controlling option;
no production change, new flag, fallback or tolerance relaxation is needed.

## Evidence and remaining boundary

Ubuntu WSL, nice 10, CARGO_BUILD_JOBS=2, RUST_TEST_THREADS=2,
CARGO_INCREMENTAL=0:

```text
cargo test --profile difftest -p solidcore --test issue_2349_closed_angle -- --nocapture
4 passed; 0 failed; 2 ignored
```

Fmt, focused Clippy (`--profile difftest -p solidcore --test
issue_2349_closed_angle -- -D warnings`), wasm check and all four Python probe
tests passed. The unchanged native executive runner also completed:

```text
PASS blend.fillet_executive_p3 cases=798 passed=796 failed=0 deferred=0 max_rel_err=0e0 with_values=796 no_values=0 vacuous=0 rejected=2
```

That native unit result uses the configured 1.0.1 libraries; the fork sample
comparison is the distinct evidence above, not a replacement reference pin.

The two ignored tests still state the unfinished certified self-mitre and
fitting-size contracts. The new sample regression does not implement them.
The original ten default-option checkpoint records and their strict diagnostic
remain unchanged, so unmatched-option differences are still reported rather
than hidden. This result does not certify complete-body topology, every
generated field, all other inputs, or the separately observed point-step state.
It does not establish that the construction is ready for flag deletion.

Three fresh native processes per reference reproduced the r3 records with
zero within-reference differences. Its eight cross-reference mass-property
differences remain a separate diagnosis; the differing source/build/reference
configuration must be isolated before assigning a root cause. Full r0.5 body
and mass comparison under matching options also remains open.

The exact/certified curved self-intersection representation and safe fit-size
semantics at the open curvature limit still await the concrete design answer:
https://github.com/tangency-labs/solidcore/issues/2243#issuecomment-6622359026
No competing construction or shared size-margin change is being implemented
while that question is pending.
