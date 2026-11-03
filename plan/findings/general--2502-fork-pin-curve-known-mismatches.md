# curve: cases that move with the difftest pin to the SMLib fork (2026-11-03)

Units: curve.composite_curve_p1, curve.composite_curve_p2, curve.crv_on_surf, curve.curve_p1, curve.curve_p3, curve.curve_p4, curve.curve_p5, curve.curve_p6, curve.curve_p7, curve.curve_p8, curve.curve_p9, curve.ellipse_p1, curve.ellipse_p2, curve.line, curve.offset_curve, curve.projected_curve, curve.tangent_field.

#2512 (#2502) pins the difftest's C++ reference to the SMLib fork, tangency-labs/smlib
`v1.1.0-alpha.1` (`8eac63c`): usd-brep v1.0.1 (`3fd5b95`, what solidcore was ported from)
plus the fork's changes. The Rust side here is the v1.0.1 port with the default flags
(`UpstreamFixes` off). The cases below agreed with v1.0.1 and disagree with the fork, so
`crates/difftest/known_mismatches/<unit>.toml` lists each one (the C++ and Rust value of
its first mismatching field are in its entry). They pass again as the fork's changes are
ported (#2502 ask 3, the flags), and each entry then fails its unit until it is removed.

**Inputs:** the listed cases. **C++:** the fork at `8eac63c`. **Rust:** this PR's head.

**Suspected cause:** the fork commits below, which changed the unit's own upstream files
since v1.0.1 (`git log 3fd5b95..8eac63c -- <upstream files>`). Where none did, the unit
moved through code it calls, and the list is its dependencies' commits instead. This
names candidates; it is not established case by case.

## curve.composite_curve_p1

- 68 of 3096 cases; first mismatching field: `values.out.len` 47, `values.out` 21.
- Fork commits:
  - 19310fa SMLib: composite curves, gaps, pseudo boxes, point sets and the geometry utilities pass on statuses
  - 54148db SMLib: draft corners that keep the curves' curvature

## curve.composite_curve_p2

- 10 of 676 cases; first mismatching field: `values.out.len` 8, `ints.status` 1, `values.out` 1.
- Fork commits:
  - 19310fa SMLib: composite curves, gaps, pseudo boxes, point sets and the geometry utilities pass on statuses
  - 54148db SMLib: draft corners that keep the curves' curvature

## curve.crv_on_surf

- 22 of 2407 cases; first mismatching field: `values.out.len` 20, `values.out` 2.
- Fork commits:
  - e7711bb SMLib: curve classification and trimming tools pass on their statuses

## curve.curve_p1

- 2 of 3096 cases; first mismatching field: `values.out.len` 1, `values.out` 1.
- Fork commits:
  - be45e8d curves: a loop's normal is its area vector, so a draft leans the way it is asked
  - 02fc29b SMLib: a curve's tight box holds the knots inside its interval only
  - 51eefa5 SMLib: a tangential contact at a curve's end is found there; the offset's stitch keeps its tolerance
  - 02e0dde SMLib, capi: the fuzzers' sanitizer findings fixed at their causes
  - 4f5ad88 SMLib: a range solution's end is refined against the end of the range
  - 8ea6a20 SMLib: refining a range solution takes an outside parameter as an answer
  - 1612a7d SMLib: the curve and surface SE() calls pass on or discard with a reason
  - 8088f5e SMLib: a disjoint interval is no work, and new objects are freed on failure
  - 8311f8e SMLib: curves pass on their statuses
  - 04b42b4 SMLib: curves coincide only where their separation vanishes
  - c2e3d9d SMLib: two curves that change sides cross; they don't coincide

## curve.curve_p3

- 250 of 3085 cases; first mismatching field: `values.out` 246, `values.out.len` 4.
- Fork commits:
  - be45e8d curves: a loop's normal is its area vector, so a draft leans the way it is asked
  - 02fc29b SMLib: a curve's tight box holds the knots inside its interval only
  - 51eefa5 SMLib: a tangential contact at a curve's end is found there; the offset's stitch keeps its tolerance
  - 02e0dde SMLib, capi: the fuzzers' sanitizer findings fixed at their causes
  - 4f5ad88 SMLib: a range solution's end is refined against the end of the range
  - 8ea6a20 SMLib: refining a range solution takes an outside parameter as an answer
  - 1612a7d SMLib: the curve and surface SE() calls pass on or discard with a reason
  - 8088f5e SMLib: a disjoint interval is no work, and new objects are freed on failure
  - 8311f8e SMLib: curves pass on their statuses
  - 04b42b4 SMLib: curves coincide only where their separation vanishes
  - c2e3d9d SMLib: two curves that change sides cross; they don't coincide

## curve.curve_p4

- 110 of 3260 cases; first mismatching field: `values.out.len` 109, `values.out` 1.
- Fork commits:
  - be45e8d curves: a loop's normal is its area vector, so a draft leans the way it is asked
  - 02fc29b SMLib: a curve's tight box holds the knots inside its interval only
  - 51eefa5 SMLib: a tangential contact at a curve's end is found there; the offset's stitch keeps its tolerance
  - 02e0dde SMLib, capi: the fuzzers' sanitizer findings fixed at their causes
  - 4f5ad88 SMLib: a range solution's end is refined against the end of the range
  - 8ea6a20 SMLib: refining a range solution takes an outside parameter as an answer
  - 1612a7d SMLib: the curve and surface SE() calls pass on or discard with a reason
  - 8088f5e SMLib: a disjoint interval is no work, and new objects are freed on failure
  - 8311f8e SMLib: curves pass on their statuses
  - 04b42b4 SMLib: curves coincide only where their separation vanishes
  - c2e3d9d SMLib: two curves that change sides cross; they don't coincide

## curve.curve_p5

- 6 of 3140 cases; first mismatching field: `values.out` 6.
- Fork commits:
  - be45e8d curves: a loop's normal is its area vector, so a draft leans the way it is asked
  - 02fc29b SMLib: a curve's tight box holds the knots inside its interval only
  - 51eefa5 SMLib: a tangential contact at a curve's end is found there; the offset's stitch keeps its tolerance
  - 02e0dde SMLib, capi: the fuzzers' sanitizer findings fixed at their causes
  - 4f5ad88 SMLib: a range solution's end is refined against the end of the range
  - 8ea6a20 SMLib: refining a range solution takes an outside parameter as an answer
  - 1612a7d SMLib: the curve and surface SE() calls pass on or discard with a reason
  - 8088f5e SMLib: a disjoint interval is no work, and new objects are freed on failure
  - 8311f8e SMLib: curves pass on their statuses
  - 04b42b4 SMLib: curves coincide only where their separation vanishes
  - c2e3d9d SMLib: two curves that change sides cross; they don't coincide
  - ca7d1f5 capi: no shim drops a status SMLib returns

## curve.curve_p6

- 62 of 3138 cases; first mismatching field: `values.out.len` 61, `values.out` 1.
- Fork commits:
  - be45e8d curves: a loop's normal is its area vector, so a draft leans the way it is asked
  - 02fc29b SMLib: a curve's tight box holds the knots inside its interval only
  - 51eefa5 SMLib: a tangential contact at a curve's end is found there; the offset's stitch keeps its tolerance
  - 02e0dde SMLib, capi: the fuzzers' sanitizer findings fixed at their causes
  - 4f5ad88 SMLib: a range solution's end is refined against the end of the range
  - 8ea6a20 SMLib: refining a range solution takes an outside parameter as an answer
  - 1612a7d SMLib: the curve and surface SE() calls pass on or discard with a reason
  - 8088f5e SMLib: a disjoint interval is no work, and new objects are freed on failure
  - 8311f8e SMLib: curves pass on their statuses
  - 04b42b4 SMLib: curves coincide only where their separation vanishes
  - c2e3d9d SMLib: two curves that change sides cross; they don't coincide

## curve.curve_p7

- 73 of 3067 cases; first mismatching field: `values.out.len` 52, `values.out` 21.
- Fork commits:
  - be45e8d curves: a loop's normal is its area vector, so a draft leans the way it is asked
  - 02fc29b SMLib: a curve's tight box holds the knots inside its interval only
  - 51eefa5 SMLib: a tangential contact at a curve's end is found there; the offset's stitch keeps its tolerance
  - 02e0dde SMLib, capi: the fuzzers' sanitizer findings fixed at their causes
  - 4f5ad88 SMLib: a range solution's end is refined against the end of the range
  - 8ea6a20 SMLib: refining a range solution takes an outside parameter as an answer
  - 1612a7d SMLib: the curve and surface SE() calls pass on or discard with a reason
  - 8088f5e SMLib: a disjoint interval is no work, and new objects are freed on failure
  - 8311f8e SMLib: curves pass on their statuses
  - 04b42b4 SMLib: curves coincide only where their separation vanishes
  - c2e3d9d SMLib: two curves that change sides cross; they don't coincide

## curve.curve_p8

- 71 of 3069 cases; first mismatching field: `values.out.len` 53, `values.out` 18.
- Fork commits:
  - be45e8d curves: a loop's normal is its area vector, so a draft leans the way it is asked
  - 02fc29b SMLib: a curve's tight box holds the knots inside its interval only
  - 51eefa5 SMLib: a tangential contact at a curve's end is found there; the offset's stitch keeps its tolerance
  - 02e0dde SMLib, capi: the fuzzers' sanitizer findings fixed at their causes
  - 4f5ad88 SMLib: a range solution's end is refined against the end of the range
  - 8ea6a20 SMLib: refining a range solution takes an outside parameter as an answer
  - 1612a7d SMLib: the curve and surface SE() calls pass on or discard with a reason
  - 8088f5e SMLib: a disjoint interval is no work, and new objects are freed on failure
  - 8311f8e SMLib: curves pass on their statuses
  - 04b42b4 SMLib: curves coincide only where their separation vanishes
  - c2e3d9d SMLib: two curves that change sides cross; they don't coincide

## curve.curve_p9

- 147 of 3126 cases; first mismatching field: `values.out` 114, `values.out.len` 33.
- Fork commits:
  - be45e8d curves: a loop's normal is its area vector, so a draft leans the way it is asked
  - 02fc29b SMLib: a curve's tight box holds the knots inside its interval only
  - 51eefa5 SMLib: a tangential contact at a curve's end is found there; the offset's stitch keeps its tolerance
  - 02e0dde SMLib, capi: the fuzzers' sanitizer findings fixed at their causes
  - 4f5ad88 SMLib: a range solution's end is refined against the end of the range
  - 8ea6a20 SMLib: refining a range solution takes an outside parameter as an answer
  - 1612a7d SMLib: the curve and surface SE() calls pass on or discard with a reason
  - 8088f5e SMLib: a disjoint interval is no work, and new objects are freed on failure
  - 8311f8e SMLib: curves pass on their statuses
  - 04b42b4 SMLib: curves coincide only where their separation vanishes
  - c2e3d9d SMLib: two curves that change sides cross; they don't coincide

## curve.ellipse_p1

- 156 of 2527 cases; first mismatching field: `values.rec` 89, `values.rec.len` 67.
- Fork commits:
  - cbeb14b SMLib: an open elliptic arc ending near its seam keeps its own end angle when its polar conversion is set up
  - 672a437 SMLib: the analytic surfaces and curves pass on their statuses
  - 5c0b5ed SMLib: an arc moved a turn past its interval keeps a parameter length

## curve.ellipse_p2

- 13 of 2107 cases; first mismatching field: `values.rec` 13.
- Fork commits:
  - cbeb14b SMLib: an open elliptic arc ending near its seam keeps its own end angle when its polar conversion is set up
  - 672a437 SMLib: the analytic surfaces and curves pass on their statuses
  - 5c0b5ed SMLib: an arc moved a turn past its interval keeps a parameter length

## curve.line

- 33 of 2105 cases; first mismatching field: `values.rec.len` 28, `values.rec` 5.
- Fork commits:
  - 672a437 SMLib: the analytic surfaces and curves pass on their statuses

## curve.offset_curve

- 53 of 3123 cases; first mismatching field: `values.out` 53.
- Fork commits:
  - 4803d6d SMLib: offsets and shells pass on the statuses they were dropping

## curve.projected_curve

- 4 of 3047 cases; first mismatching field: `values.out.len` 4.
- Fork commits (none on its own files; its dependencies'):
  - be45e8d curves: a loop's normal is its area vector, so a draft leans the way it is asked
  - 02fc29b SMLib: a curve's tight box holds the knots inside its interval only
  - 51eefa5 SMLib: a tangential contact at a curve's end is found there; the offset's stitch keeps its tolerance
  - 02e0dde SMLib, capi: the fuzzers' sanitizer findings fixed at their causes
  - 4f5ad88 SMLib: a range solution's end is refined against the end of the range
  - 8ea6a20 SMLib: refining a range solution takes an outside parameter as an answer
  - 1612a7d SMLib: the curve and surface SE() calls pass on or discard with a reason
  - 8088f5e SMLib: a disjoint interval is no work, and new objects are freed on failure
  - 8311f8e SMLib: curves pass on their statuses
  - 04b42b4 SMLib: curves coincide only where their separation vanishes
  - c2e3d9d SMLib: two curves that change sides cross; they don't coincide
  - ca7d1f5 capi: no shim drops a status SMLib returns
  - and 3 more

## curve.tangent_field

- 301 of 1355 cases; first mismatching field: `values.out.len` 208, `ints.status` 88, `values.out` 5.
- Fork commits (none on its own files; its dependencies'):
  - be45e8d curves: a loop's normal is its area vector, so a draft leans the way it is asked
  - 02fc29b SMLib: a curve's tight box holds the knots inside its interval only
  - 51eefa5 SMLib: a tangential contact at a curve's end is found there; the offset's stitch keeps its tolerance
  - 02e0dde SMLib, capi: the fuzzers' sanitizer findings fixed at their causes
  - 4f5ad88 SMLib: a range solution's end is refined against the end of the range
  - 8ea6a20 SMLib: refining a range solution takes an outside parameter as an answer
  - 1612a7d SMLib: the curve and surface SE() calls pass on or discard with a reason
  - 8088f5e SMLib: a disjoint interval is no work, and new objects are freed on failure
  - 8311f8e SMLib: curves pass on their statuses
  - 04b42b4 SMLib: curves coincide only where their separation vanishes
  - c2e3d9d SMLib: two curves that change sides cross; they don't coincide
  - ca7d1f5 capi: no shim drops a status SMLib returns
  - and 3 more

Measured in #2512 (WSL, the fork's libraries built with `harness/cpp/CMakeLists.txt`'s flags).
