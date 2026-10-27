# intersect.surface_intersector_p5: an empty trace is not reversed, opt-in (2026-10-27)

Units: `intersect.surface_intersector_p4`, `intersect.surface_intersector_p5`,
`intersect.adv_surface_intersector` (all green). Issue #391, record WR-80
(`plan/upstream-issues/WR-80--addendum-2026-10-27-rust-port-opt-in.md`). PR: #1285.

## What was done

`Kernel::fixes.ssi_empty_trace_not_reversed` (off by default) skips the "trace the other way
too" block of `FindInteriorCurves` / `DoIntersection` when `TraceCurve` left no points, in the
three blocks that had no `!curve_points.is_empty()` guard (`surface_intersector/p4.rs`,
`p5.rs`, `adv_surface_intersector.rs`; `DoPointIntersection` already had it). One predicate,
`SurfaceIntersector::trace_may_go_backwards`; with the flag off it is the old `!curve_is_closed`.

## Evidence (origin/main 4d7ef7021 + this branch, ssi-corpus 119 cases, tol 1e-6)

| Rust engine | main | branch, flag off | branch, flag on |
|---|---|---|---|
| ERROR | 2 | 2 | **0** |
| INCOMPLETE | 25 | 25 | 26 |
| WRONG TOPOLOGY | 11 | 11 | 12 |
| PASS | 81 | 81 | 81 |

- Flag off: the report is identical to a clean `origin/main` run (JSON-normalised compare).
- Flag on: only two cases change, `c2n_quartic_contact_parabola` ERROR -> WRONG TOPOLOGY and
  `y23_t2_ex5_two_loops` ERROR -> INCOMPLETE (both return no curves). C++ and the `ssi` /
  `solidcore-ssi` engines are unchanged. The cases still need the tangent-curve fixes of
  #392/#394/#395 to pass.
- Unit difftests of the four touched units match main's CI counts (1542 / 1731 / 1931 / 16688
  cases, 0 deferred, 0 vacuous).
- `crates/solidcore/tests/issue_391_empty_trace.rs`: two hand-built degree (2, 4) Bezier
  patches in 4th-order contact (formulas in the file, no corpus data). Default:
  `SM_ERR_NULL_POINTER` from `reverse_curve_direction`; flag on: `Ok`.

## Not covered

- Only the `p5` boundary start-point loop has a reproducing test. The `p4` and `adv`
  `FindInteriorCurves` blocks carry the same predicate but no case reaches them with an empty
  trace.
- No C++ comparison with the flag on is possible (C++ has the defect); the flag-on result is
  not shown to be geometrically complete.

## Seen on main, not caused by this change

`c5_cyl_quilt_seam` is the one parity mismatch of the corpus run on main
(`values.curve000.s1_at_uv.len: cpp=0 rust=12003`; both engines PASS).

## Harness

`SSI_SC_FIXES=<name>[,<name>]` enables `UpstreamFixes` for the `rust` engine of
`--ssi-corpus`; such a run writes `*_fixes_<names>` reports and raw directories and records a
`fixes` array, so it never overwrites the parity run.
