# Final affine-candidate validation for PR #2840

Production candidate: `f79b779f5`. Revision `f290ed2de` changes only a test helper
call for Rust 1.85 compatibility; its production code is identical.
Latest merged main: `a418def44`.
Native reference: pinned v1.0.1, key `v1-3fd5b9548d19-726643fb17a52b7d`.

## Passing checks

- Formatting, workspace/all-target clippy with warnings denied, and wasm check pass.
- 3,566 library unit tests pass; four are ignored.
- All 19 focused regression tests pass; one diagnostic export helper is ignored
  in the ordinary run and passes when explicitly invoked.
- The final test-only revision passes formatting, clippy and both private affine
  parallel-identity tests again.
- All four exported small/control body and edge artifacts are byte-identical to
  the corrected historical `09495e082` exports. The before/after body figure
  therefore still depicts the current production geometry.

## Fresh native comparison

This is an actual Rust/C++ comparison, separate from the Rust-only digest run.
All twelve affected units were explicitly selected, including blocked units.

| Unit | Cases | Passed | Failed |
|---|---:|---:|---:|
| curve.curve_p1 | 3102 | 3049 | 53 |
| curve.curve_p4 | 3263 | 3234 | 24 |
| curve.curve_p5 | 3140 | 3130 | 7 |
| curve.curve_p7 | 3067 | 3051 | 16 |
| curve.curve_p9 | 3126 | 3119 | 7 |
| geom.global_solver | 630 | 628 | 2 |
| intersect.topology_solver_p1 | 3592 | 3568 | 24 |
| intersect.topology_solver_p2 | 2575 | 2532 | 43 |
| topology.brep_p1 | 2636 | 2634 | 1 |
| topology.brep_p10 | 1436 | 1417 | 4 |
| blend.hcr | 481 | 463 | 1 |
| primitive.api_trim_ops | 211 | 196 | 5 |
| **Total** | **27259** | **27021** | **187** |

The other 51 cases retain existing classifications: 13 parked, nine known
mismatches and 29 rejected. This change introduces no new parking and does not
loosen the relative-only parity criterion.

All 173 historical failing-case mismatch payloads are unchanged. Eleven added
reproduction fixtures contribute eleven additional failures. Three existing cases
newly fail: `curve.curve_p1` gen1172 and gen2883, and
`intersect.topology_solver_p1` gen1646. No historical failed case becomes a pass.

The four existing global-query output changes (including two that remain passes)
are audited in [the affine caller record](curve.curve_p1--affine-caller-audit-2840.md).
The newly refused ordinary finite zero-tolerance query is documented in
[the API-contract record](intersect.topology_solver_p1--zero-tolerance-affine-query-2840.md).
Neither record classifies every disagreement as a native defect.

## Limits and blockers

The full integration suite is incomplete. The full 280-unit native sweep was run
on historical `09495e082`, not on this affine candidate; see
[the historical full caller sweep](curve.curve_p1--full-caller-sweep-2840.md).
The fresh comparison above covers the twelve affected units only.

All twelve affected units remain blocked and PR #2840 remains a draft. Maintainer
clarification is pending on out-of-bounds query semantics, zero distance tolerance,
and whether silhouette accuracy repair belongs in this PR or a separate follow-up.
The previously recorded curve-pair maxima and trim/caller uncertainties also remain
unresolved. Passing the focused geometry tests does not establish native parity
or make this candidate ready to merge.
