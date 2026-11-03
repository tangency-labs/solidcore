# blend.fillet_corner_p3: reachable closed teardrop angle (#2349, 2026-11-02)

Units: `blend.fillet_corner_p3`, `blend.fillet_geom_p3`, `blend.fillet_executive_p3`.
PR: https://github.com/tangency-labs/solidcore/pull/2381.
Initial base: `fdc55bbb4`; merged upstream main `04d51a64e`, then `23ea17c1c`.
This contribution provides preparatory evidence and a native allocation-parity
correction, not the certified mitre event. The newly reachable checkpoint failures
below are resolved, and the executive unit's green status is restored on native
evidence.

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

### Native parity failure and correction

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

Traced cause: `FilletSolver::CalcFilletSurfaces` splits the closed fillet through
`SmFilletGeom::TopologySplit` (`SmFilletSolver.cpp:2478`). Its new rail vertices
are `new (pPseudoBrep) SmFilletVertex()` (`SmFilletGeom.cpp:12653`). That constructor
calls `SmVertex()` (`SmFilletGeom.cpp:1316`), which initializes `m_sZoneTol3d` to
zero (`SmVertex.h:118`). Rust's `FilletBrep::topology_split` instead used
`Brep::add_vertex`, initializing tolerance to `1e-5`. It now allocates
`Vertex::default()` through `add_unlinked_vertex_record`, preserving the native zero
without altering later tolerance updates. `SmTopology::operator new(size, SmObject*)`
only supplies the context (`SmTopology.cpp:571`); it does not insert the vertex
into the B-rep's list. `AddFilletVertex` retains it in the corner instead.
The first correction to the tolerance alone exposed this second allocation
error: C++ lists one pseudo vertex, while Rust listed three at both checkpoints.
The corrected allocation preserves both constructor state and list ownership.
The focused topology-split test checks the new record's position, tolerance and
absence from the B-rep's vertex list. It selects the first new vertex: later rail
setup can allocate a separate shell placeholder before returning its error.
This is a general allocation correction,
not a tolerance reset for the teardrop.

All ten cases remain in the runner; no tolerance is relaxed, no key is removed,
and no `cpp_defect` exemption is claimed.

At corrected source commit `4fe560ec4c2c49cad485b47deb396e73a520e574`, the focused
native rerun gives:

```text
PASS blend.fillet_executive_p3 cases=10 passed=10 failed=0 deferred=0 max_rel_err=0e0 with_values=10 no_values=0 vacuous=0
PASS blend.fillet_executive_p3 cases=798 passed=796 failed=0 deferred=0 max_rel_err=0e0 with_values=796 no_values=0 vacuous=0 rejected=2
```

`--status-of blend.fillet_executive_p3` prints `green`, and its plan status is
restored on that evidence. The pre-correction full unit run had only those two
failing cases. The restored unit participates in the all-green gate; its coverage
must not be waived in the exact-base artifact ratchet.

Initial verification before the correction: fmt check and focused Clippy with
`-D warnings` pass; the integration test has 2 passes and 2 explicitly ignored
future contracts. Running those contracts explicitly produces 2 failures with
the diagnostics above. That initial run did not verify the full workspace, wasm
or the all-green gate.

The focused test's final vertex selection is corrected in
`2e143ad0eaf9cc43cfdecb5a5bbf1866a3af8271` (production allocation unchanged).
That regression passes. Native verification of the affected unit and callers
at this source gives:

```text
PASS blend.fillet_geom_p3 cases=325 passed=318 failed=0 deferred=0 max_rel_err=0e0 with_values=318 no_values=0 vacuous=0 rejected=7
PASS blend.fillet_solver_p1 cases=351 passed=351 failed=0 deferred=0 max_rel_err=0e0 with_values=351 no_values=0 vacuous=0
PASS blend.fillet_solver_p4 cases=651 passed=651 failed=0 deferred=0 max_rel_err=0e0 with_values=651 no_values=0 vacuous=0
PASS blend.fillet_solver_p5 cases=100 passed=100 failed=0 deferred=0 max_rel_err=0e0 with_values=100 no_values=0 vacuous=0
PASS blend.fillet_corner_p3 cases=1013 passed=994 failed=0 deferred=0 max_rel_err=0e0 with_values=994 no_values=0 vacuous=0 rejected=19
PASS blend.fillet_executive_p3 cases=798 passed=796 failed=0 deferred=0 max_rel_err=0e0 with_values=796 no_values=0 vacuous=0 rejected=2
```

The full native production gate also passes (the subsequent test-only selection
correction does not change production code):

```text
difftest --green --jobs 2 --baseline bin/issue-2349-baseline-23ea17c1c
difftest: 280 unit(s), 280 passed, 0 partial, 0 failed
```

The baseline is the exact main `23ea17c1c21a587e5bd50653e7c73c490f34bbbe`
`difftest-reports` artifact from successful CI run `41441037947`, downloaded with
Windows `gh`. The unfiltered run enforces the value-share floor and coverage
ratchet; no historical committed reports or waivers are used. Final-source
workspace Clippy (`--all-targets -- -D warnings`) and wasm32 check also pass.
Formatting passes. The unoptimized whole-kernel attempt was interrupted while
the existing `issue_478_full_turn_nested_shell_is_a_cavity` test was still running;
it is not counted as a pass. The CI test profile is `difftest`
(`.github/workflows/ci.yml`), with both overflow checks and debug assertions
enabled (`Cargo.toml`). No test or assertion is removed when using that profile.

Explicitly running both ignored future contracts at final source still gives two
failures: no blend event, and the proposed radius `1.8749999999999993` builds with
level-2 geometry warnings. These remain future implementation contracts, not
passing coverage.

### Exact event work

The live #2349 brief requires an End-like self-mitre, cutting the two ends on
their bisecting surface with exact pcurves, or a certified refusal naming the size
that builds. #2353's deletion of the old corner's G1 paths is merged here through
main `23ea17c1c`. #2307 merged as
`379f9dafe0e64a2c9c104a50fc74b606ef8c0c1e` and is included through main
`ebe6722e09169cce29a20c455fc6ec606b57caa0`; its implementation is not duplicated
here. There is no new
kernel flag, old-corner deletion, fit-size formula or claim that the old fitted intersection
meets the exact-pcurve rule.

The current `blend_corner::exact_mitre_curve` requires straight spines and equal
cylinders; it returns `false` for this curved-spine round. The merged
#2307 End/Join dispatch does not supply the exact/certified bisector section
and its two pcurves or prove a fitting radius. This preparatory contribution does
not implement a competing End interface.

A follow-up implementation can use this operand to test
both ends of one blend at one vertex, distinguish its angle from a closed G1
Join, and enable the failing contract. Exact/certified bisector intersections,
exact pcurves on both sides, geometric fit-size derivation, and boundary-size
tests remain work for the event. General fold handling remains independent.
This fixture covers T=2 only; a T=3 seam and a revolved kinked seam are not covered.


Full kernel verification: `nice -n 10 cargo test --profile difftest -p solidcore` with `CARGO_BUILD_JOBS=2` and `RUST_TEST_THREADS=2` exited 0: 4,164 passed, 0 failed, 34 ignored across library, integration and documentation targets. The library accounts for 3,061 passes and 4 ignored tests. Ignored tests are not counted as passes.

### Latest-main refresh after #2307 merged

Merged main `ebe6722e09169cce29a20c455fc6ec606b57caa0` without conflicts in
`d747975c4149984dd69e0b98c26d5d0b2f141b2d`. No production correction, checkpoint
field or tolerance was changed. Post-merge verification exited 0:

- `cargo fmt --all --check`.
- CI-profile focused topology-split tests: 5 passed, 0 failed.
- CI-profile issue integration target: 2 passed, 0 failed, 2 ignored future contracts.
- Native affected units: geom_p3 318 passes/7 rejections; solver_p1 351;
  solver_p4 651; solver_p5 100; corner_p3 994 passes/19 rejections;
  executive_p3 796 passes/2 rejections. All six pass with zero failures,
  zero deferred cases and `max_rel_err=0e0`.

Commands used Ubuntu WSL, `CARGO_BUILD_JOBS=2`, `RUST_TEST_THREADS=2`, nice 10,
and the unchanged shared C++ build. The earlier full green gate, full kernel
suite, Clippy and wasm results above predate this main refresh; they are not
claimed as rerun at this merged revision. GitHub CI still needs maintainer
approval for the fork workflow; the inspected run 41458895332 has
`action_required` and zero jobs. This account has no upstream push/maintain/admin
permission. Required review is also outstanding.

## Rebase and #2499 audit (2026-11-03)

Rebased the owned PR branch onto `fb937dd8d991072cf5b6fa7d16cf46ea110ea238`
with the user's explicit authorization to replace its published history using
`--force-with-lease`. The rebase was clean. Original commit IDs above describe
historical evidence; in the rebased history the failing fixture commit is
`5860f15a6`, before the allocation correction `e24b4866c`.

#2499 fixes the new-blend overflow judgment in `rails_fit.rs`: an untrimmed rail
past the cutting plane at its closed corner is not an overflow off the face.
Its corner remains the legacy closed Nx1 mitre; #2349 still owns its replacement.
Its teardrop is a six-pole, multi-span cubic with a 90-degree tip and least profile
curvature radius about 4.69, checked by completed bodies and independent volume
at radii 0.5 and 1. This PR's four-pole single-span cubic has a 53.13-degree tip
and curvature radius 1.875 at its midpoint. Its serialized operand and ten
intermediate/full native checkpoints complement that end-to-end test rather than
duplicate its geometry or acceptance measurement.

The added post-rebase regression follows real corner creation, smooth-edge
processing, vertex creation and adjustment, solver geometry and corner edges.
It checks the two new rail-split vertices at the solver-geometry checkpoint for
exactly zero tolerance and absence from the B-rep owning vertex list. It runs
both default flags and #2499's profile (every fix except geometric G2 continuity).
The small-radius body control now runs both profiles as well. No existing
fixture, expected numeric value, checkpoint field or tolerance is changed.

### Reference and construction policy

The normative reference for the redesigned blend is the SMLib fork. Read only
permitted fork sources at v1.1.0-alpha,
`a18bbe523b9c7819f1a9da8c40a54c15f9936125`: `SmFilletGeom.cpp`'s `TopologySplit`
(starts at line 13105) still allocates `new (pPseudoBrep) SmFilletVertex()` for
both split rails and attaches them via `AddFilletVertex`. Its `SmFilletVertex`
constructor calls `SmVertex()`, whose `SmVertex.h` constructor starts tolerance
at zero. `SmTopology.cpp`'s placement allocator supplies context without inserting
into an owning topology list. Thus the allocation correction agrees with the
fork contract too. The configured difftest still uses the unchanged 1.0.1 pin;
its results verify the shared allocation path, not fork executable parity for a
redesigned construction or certification of the self-mitre.

The correction is generic to split allocation, with no per-input branch, new
flag, fallback/retry chain, global or thread-local state, or library thread.
It fixes the still-reached parity path before replacement and is deleted with
its owning old construction, not kept as a permanent fallback. The regression
fixture was committed before the fix. No SolidCraft gallery case changes state
in this contribution; no gallery row or render changes are claimed.

### Fresh post-rebase evidence

`cargo fmt --all --check` and workspace/all-target Clippy with `-D warnings`
pass. The issue integration target has 3 passes and 2 ignored future contracts.
Both flag profiles reach the two zero-tolerance, unlinked split vertices and
complete corner edges through the legacy closed corner. The radius-0.5 body is
valid at Level2 in both profiles. The landed #2499 closed-rim test also passes
at both of its radii. No acceptance expectation changes are needed.

Explicitly running the two ignored contracts still fails: the corner has no
certified blend event, and the suggested radius `1.8749999999999993` produces
Level2 geometry errors. They remain known failing follow-up contracts, not
passes. The first focused library-test build failed before tests ran because
rustc reported `Cannot allocate memory (os error 12)` reading solver `p2.rs`;
this resource failure is separate from a test assertion failure.

The fresh wasm check passes. The configured native executive runner reports:

```text
PASS blend.fillet_executive_p3 cases=10 passed=10 failed=0 deferred=0 max_rel_err=0e0 with_values=10 no_values=0 vacuous=0
```

This is the unchanged `--filter issue-2349` checkpoint set, including solver
geometry, corner edges and completed radius-0.5/radius-3 operations. It confirms
that #2499 did not remove the split-allocation coverage or require new expected
numeric values.

The full executive rerun using that freshly built binary also passes:

```text
PASS blend.fillet_executive_p3 cases=798 passed=796 failed=0 deferred=0 max_rel_err=0e0 with_values=796 no_values=0 vacuous=0 rejected=2
```

The redundant `cargo run` invocation was stopped while waiting on the library
retry's build lock; the completed runner was invoked directly under nice 10.
No native case, field, rejection, or comparison tolerance was changed.

The focused library-test retry also failed during compilation, this time writing
the incremental `query-cache.bin`: `Cannot allocate memory (os error 12)`.
The shared machine's load had risen to 174.68. No library test assertion ran;
the focused library check is not claimed green. A third heavy build was not
started under that load. Re-run `cargo test --profile difftest -p solidcore
--lib topology_split` under the documented machine budget when resources allow.
This does not invalidate the fresh integration and native passes above, but
leaves the focused library verification incomplete locally.
After the query-cache error, the still-running compiler owned by this task was
terminated by its verified PID to release resources; Cargo reports SIGTERM
alongside the prior allocation error. No other thread's process was stopped.
