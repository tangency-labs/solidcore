# blend.fillet_corner_p3: reachable closed teardrop angle (#2349, 2026-11-02)

Units: `blend.fillet_corner_p3`, `blend.fillet_executive_p3`.
PR: https://github.com/tangency-labs/solidcore/pull/2381.
Initial base: `fdc55bbb4`; merged upstream main `04d51a64e` before final verification.
This is preparatory evidence, not the
implementation of the certified mitre event. The executive unit is marked blocked
because the newly reachable checkpoints expose the parity differences below.

## Reproducible operand

`crates/solidcore/tests/issue_2349_closed_angle.rs` constructs a single, nonrational,
clamped cubic Bezier with poles `(0,0,0)`, `(20,-10,0)`, `(20,10,0)`, `(0,0,0)`,
knots `[0,1]`, multiplicities `[4,4]`. Extrude it 10 along Z with both caps.
There are no interior knots, profile joins, booleans, imported upstream assets or
synthetic corner handles. The result has 3 faces, 3 edges and 2 vertices and passes
level-2 validation before filleting. Select its upper closed rim only, circular
section, continuity 1, thumbweight 1.

The endpoint derivatives are `(60,-30,0)` and `(-60,-30,0)`, cosine -0.6.
The outgoing directions at the seam enclose 53.130102354 degrees: this is far from
G1 at the factory's ten-degree tolerance. The rim is one edge whose start and end
are the same vertex, not two edges meeting at an ordinary Turn.

The real `FilletSession` factory selects one `NByOneClosed` corner, S=1, T=2.
At radius 0.5 its real corner-vertex stage completes, with `g1_case=false`,
`blending=false`, two fillet vertices, one `FilletXFillet` edge, a
`FilletXFillet` vertex, and `FilletGeom.need_split=true`. The integration test
checks these observations with `fillet_blend_corner` both off and on.

The ASCII input `fixtures/blend.fillet_executive_p3/data/issue_2349_teardrop.smb`
was written from that Rust-constructed operand, before any filleting. The existing
file-source difftest reads the same input independently on Rust and pinned C++.
No C++ state or result is fed into the Rust construction.
To regenerate, set `SOLIDCORE_2349_WRITE_FIXTURE` to the desired output path and
run the `teardrop_reaches_the_real_closed_angle_topology` test (one test thread).
The committed ASCII copy has trailing spaces stripped; no numbers were edited.

## Observed behavior and failing contract

On the base implementation:

| Radius | Rust public executive result | Level-2 result / diagnostic |
| --- | --- | --- |
| 0.5 | success | 5 faces, no validation messages |
| 1 | success | 5 faces, no validation messages (exploratory run) |
| 2 | success | 5 faces; nonmonotonic surface control net, slow curve derivatives and repeated curve control points (exploratory run) |
| 3 | status 1001 | code 3, `NByOneClosed`, `Fillet Corner Error: Vertex Computation Failure`; no fit size |
| 5 | status 1001 | code 4, `Fillet Surface Computation Failure: Fillet Face Creation Error`; no fit size (exploratory run) |

Radius 2's validation messages alone do not prove that the seam event causes its
failure; a fold in the offset/blend could be responsible. No local remedy is
proposed for it. Radius 0.5 is a positive control, not a failing construction.

Two explicitly ignored future contracts remain failing:

- `closed_angle_is_a_blend_event_instead_of_the_legacy_corner`, radius 0.5,
  `fillet_blend_corner=true`: the real factory still leaves `corner.blend=None`.
  This asserts the missing event classification, not geometric certification.
- `closed_angle_round_builds_validly_or_refuses_with_a_fit_size`, radius 3,
  `fillet_blend_corner=true` and `fillet_gn_section_generator=true`: the latter
  flag enables the existing fit-size machinery. The request is refused with
  `BallDoesNotFit`, code 8, and radius `1.8749999999999993` (`FitKind::Radius`).
  That suggested size builds but fails level-2 validation: nonmonotonic surface
  control net, zero/slow UV derivatives and repeated curve control points. The
  test checks a finite positive suggested size below the request and actually
  builds it on a fresh operand; it fails on that result's validation messages.

The fit-size failure is a fold-boundary diagnostic, not evidence that the
closed-angle event alone fixes it. At the cubic's midpoint the curvature radius
is exactly `15^2 / 120 = 1.875`. The legacy profile without the fit-size flag
instead reaches the unsized code-3 corner failure in the table. Neither ignored
test is counted as passing; enable them as the corresponding work is completed.

```sh
export CARGO_BUILD_JOBS=2
nice -n 10 cargo test -p solidcore --test issue_2349_closed_angle
nice -n 10 cargo test -p solidcore --test issue_2349_closed_angle \
  -- --ignored --nocapture
nice -n 10 cargo run --profile difftest -p difftest --features cpp -- \
  --unit blend.fillet_executive_p3 --filter issue-2349
```

## Implementation boundary and dependencies

### Newly exposed native parity blocker

At contribution `f68c1104a`, after merging main `04d51a64e`, the pinned native
comparison on the ten new file-source fixtures gives:

```text
FAIL blend.fillet_executive_p3 cases=10 passed=8 failed=2 deferred=0 max_rel_err=1e0 with_values=8 no_values=0 vacuous=0
FAIL blend.fillet_executive_p3 cases=798 passed=794 failed=2 deferred=0 max_rel_err=1e0 with_values=794 no_values=0 vacuous=0 rejected=2
PASS blend.fillet_corner_p3 cases=1013 passed=994 failed=0 deferred=0 max_rel_err=0e0 with_values=994 no_values=0 vacuous=0 rejected=19
```

| Case | Key | C++ | Rust |
| --- | --- | --- | --- |
| `issue-2349-teardrop-solver_geometry-r0.5` | `values.all_corner_vertices[8]` | 0 | 1e-5 |
| `issue-2349-teardrop-corner_edges-r0.5` | `values.all_corner_vertices[13]` | 0 | 1e-5 |

`prefix.rs` exports five values per processed corner vertex: x, y, z,
tolerance, original parameter. Both differing entries are tolerances. Corner
selection, initial topology and first vertex geometry compare successfully;
the differences become visible after solver construction and corner-edge work.
Both full-operation fixtures (radii 0.5 and 3) compare successfully. A full
operation pass does not erase an intermediate-state mismatch.

Suspected cause: a fresh or split fillet vertex retains the convenience
`Brep::add_vertex` system tolerance rather than the native constructor's zero.
`plan/findings/blend.fillet_solver_p5--pr-208-maintainer-review-2026-10-19.md`
records the same constructor distinction in other rail-split sites. That is a
lead, not a traced root cause for these two vertices. No kernel correction is
included here. All ten cases remain in the runner; no tolerance is relaxed, no key is
removed, and no `cpp_defect` exemption is claimed. The executive unit remains
blocked until this separate parity issue is resolved.
The full unit run has only those two failures. A `--green` run that omits this
blocked unit is not parity evidence for it; its lost coverage must not be waived
in the exact-base artifact ratchet.

Local verification on the merged-main code: fmt check and focused Clippy with
`-D warnings` pass; the integration test has 2 passes and 2 explicitly ignored
future contracts. Running those contracts explicitly produces 2 failures with
the diagnostics above. No full-workspace, wasm, or all-green pass is claimed.

### Exact event work

The live #2349 brief requires an End-like self-mitre, cutting the two ends on
their bisecting surface with exact pcurves, or a certified refusal naming the size
that builds. #2307 still owns the End and Join event/section interfaces; #2342
and #2353 own deletion of the old corner's G1 paths. Neither implementation is
duplicated here. There is no new kernel flag, numerical change, tolerance change,
old-corner deletion, fit-size formula or claim that the old fitted intersection
meets the exact-pcurve rule.

After those interfaces settle, the implementation can use this operand to test
both ends of one blend at one vertex, distinguish its angle from a closed G1
Join, and enable the failing contract. Exact/certified bisector intersections,
exact pcurves on both sides, geometric fit-size derivation, and boundary-size
tests remain work for the event. General fold handling remains independent.
This fixture covers T=2 only; a T=3 seam and a revolved kinked seam are not covered.
