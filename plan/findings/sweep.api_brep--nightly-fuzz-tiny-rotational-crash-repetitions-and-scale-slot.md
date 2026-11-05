# sweep.api_brep: nightly fuzz failures of 2026-10-25/26: a tiny rotational profile at the axis crashes upstream, repetitions past 360 degrees, and the scale reference's stale solution slot (2026-10-27)

Unit `sweep.api_brep`. The nightly `--fuzz 300` runs 40177686438 and 40298730423 failed on
five inputs (artifact `fuzz-failures`, folder `sweep.api_brep/`); PR #1366.

## fuzz00043 (crash): UB-271

`rotational_reps` of an `SmApiCreateCircle` circle of radius 1.566e-6 centred 3.03e-6 from
the axis, 90°, 3 repetitions, `cap: false` (`gen0158` with `/cap = false`; capped, upstream
returns `SM_ERR` before the sweep because `bCapEnds` needs 1 repetition). The C++ dies in
`SmTopologySweep::SweepEdgeHigh` (`SmTopologySweep.cpp:2165`) reading the owner of the
original edge curve that `MakeFaceWithCurves` discarded: UB-271 has the regime (a profile of
extent below the B-rep tolerance whose start vertex is within about 5e-6 of the axis or
sweeps an arc shorter than the tolerance). The Rust side returns `SM_ERR_NULL_POINTER`
(1006) for the fuzz input and its variants, without a panic.

Fix (PR #502's convention): the FFI guard `check_rotational_curves` refuses, where the
sweep runs, a profile curve of extent below 1e-5 with an end point within 2e-5 of the axis
or sweeping an arc below 2e-5 (the observed crash bounds are 5e-6 and 1e-5; the margins keep
the fuzz mutations of the tiny-scale generated cases out). The runner's Rust side applies
`ffi::validate` to the same call, so the case is a harness error on both sides (corpus
`fuzz00043_tiny_circle_at_axis` and the probed variants). Lost coverage: sub-tolerance
profiles at the axis that upstream turns into `SM_ERR`, `SM_ERR_AXIS_INSIDE_FACE` or a
degenerate one-face result (a tiny line or short arc at the axis); tiny profiles away from
the axis, which build bodies, stay covered.

## fuzz00058 (differs): the repetition bound came before upstream's own failure

`review_rot_reps_big` (a unit square 10° about an axis 5 away, 5000 repetitions) with
`/cap = false`: C++ `SM_ERR` (1001), Rust `SM_ERR_INVALID_INPUT` (1007). Upstream's
`SmTopologySweep::DoSweep` checks `SmRotationalSweepGeometry::IsSweepClosed(ii + 1)` at the
top of every repetition and fails with `SM_ERR` at the 37th (370° > 360°), after 36
repetitions ran; the port refused more than 4096 repetitions before the loop. Fix: the
4096 work bound (a documented deviation of `sweep::topology_sweep::engine`) is checked
inside the loop before the 4097th repetition, so every upstream failure inside the first
4096 comes first. Repetitions 1, 2, 3, 10, 36, 37, 100 and 4096 agreed before and after;
4097 agrees now. A linear sweep of more than 4096 repetitions still fails with
`SM_ERR_INVALID_INPUT` where upstream runs them all (the FFI admits up to 2^20).

## fuzz00089, fuzz00250 (differs) and fuzz00258 (differs): the scale reference's point solve without a solution

`from_faces` and `curve_sweep` with a scale reference or scale curve that a mutation moved
to 1e18 or 1e300 (`/sweep/scale_reference/points/0/2 = -1e+300`,
`/sweep/scale_reference/points/0/1 = 2^63-1`, `/sweep/scale_curve/end/2 = -2^63`): C++
`SM_ERR` (1001) or a body, Rust `SM_ERR_NULL_POINTER` (1006). In
`SmPrimitiveCreation::CreateCurveSweep` (`SmPrimitiveCreation.cpp:5026`, `:5037`) the
scale samples are dropped on the reference with `GlobalPointSolve(.., SM_SR_SINGLE, ..)`,
which returns `SM_SUCCESS` with no solution for a far-away point (the curve/point tree
solver finds none), and `sSolutions[0].m_vStart[0]` is read regardless: the borrowed
`SmSolutionArray` is zeroed at construction and `ReSet` only clears its size, so the
parameter is 0 at first and the last found parameter afterwards. A stand-alone replay of the
scaling block against the reference shows: `fuzz00258` finds one solution at the first
sample and none at the other 30 (slot 0 throughout), builds the scaled curve and the body;
`fuzz00089` finds none from the start (`dScale = 1 / inf = 0`, NaN points) and the
interpolation, reparameterisation and `MakeCurvesCompatible` still succeed, the `SM_ERR`
coming later. The port returned `null_pointer` at the empty solution array. Fix:
`primitive_creation::p3::geometry::scaling` keeps the stale slot as upstream does
(parameter 0 at first, then the last found).

## Verification notes

The macOS arm64 Apple clang build of the v1.0.1 pin was the local reference (Linux CI is
the parity gate). Because the unit's case list grew since the nightly runs (fixtures added
on `main`, and the fixture file of this PR), a local `--fuzz 300` picks other base cases
than the nightly did: its reproducers of the same families were `fuzz00090` (`curve_sweep`
with a scale-curve start at `-2^63`) and `fuzz00196` (`from_faces` with a scale-curve start
at `1e300`), both now agreeing (`SM_ERR` on both sides) and kept as fixtures.

The final local `--fuzz 300` (case list with this PR's fixture file, 606 cases) is
`FUZZ-FAIL cases=300 rejected=17 malformed=100 ran=181 differs=0 crashes=1 panics=0
timeouts=1`: the crash is `fuzz00214`, `review_draft_square_a90` with `/curves/1/start/2 =
4294967296` (`SmApiCreateDraftSweep` of a four-line profile with one line start moved to
z = 2^32), the fuzz-sweep trigger recorded under UB-112 ("Not guarded yet"; the corpus keeps
`fuzz00538` and `fuzz01538` of that family as `cpp_skip`); the timeout is `fuzz00252`,
`gen0055` with `/sweep/path/points/3/0 = 1e+300` (`from_edges` along a path with a 1e300
point), the HANG-43 family (corpus `fuzz01218` / `fuzz01458`, `cpp_skip`). Neither is a
sweep.api_brep regime of this PR; both reproducers are in the nightly artifact format under
the PR's scratch output and belong to the owners of those records.

Merging `origin/main` (2026-11-03) brought the HANG-5 / HANG-43 path-extent guard
(`check_path_extent`, nightly fuzz00293) and the repetition check (fuzz00132), which
`main` applies on the runner's Rust side one check at a time; this branch applies the whole
`ffi::validate` there instead (PR #502 / #1358's convention), which includes both, so the
corpus outcomes of `fuzz00293` and `fuzz00132` are unchanged.

After the merge the gate is `PASS sweep.api_brep cases=705 passed=705 failed=0 deferred=0
max_rel_err=0e0` and a local `--fuzz 300` is `FUZZ-PASS cases=300 rejected=36 malformed=105
ran=159 differs=0 crashes=0 panics=0 timeouts=0` (the larger case list picks other base cases).

## Follow-up (not taken here): a draft sweep of an open profile with a gap crashes upstream

The local `--fuzz 300` of 2026-10-27 (case list of `main` at b76307823) also found
`fuzz00292`: `review_draft_square_a90` with `/curves/0/start/0 = 3`, i.e. `SmApiCreateDraftSweep`
(height 1, angle 90°, fillet corners, caps 3) of four lines whose first starts at `(3,0,0)`
instead of `(0,0,0)` (a chain with a gap of 3). The C++ dies with a null dereference in
`SmCompositeCurve::IntersectAndTrimOffsets` (one recursion deep, called from
`CreateOffsetsOfManyCurves` ← `CreateTrimmedOffsets` ← `SmPrimitiveCreation::CreateTaperExtrude`),
the Rust side returns without a panic. It is the offset-corner family of UB-61 / UB-112,
not the sweep families above; it needs its own regime probe and guard.
