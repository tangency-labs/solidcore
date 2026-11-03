# primitive: cases that move with the difftest pin to the SMLib fork (2026-11-03)

Units: primitive.api_curve_ops, primitive.api_primitives, primitive.api_trim_ops, primitive.primitive_creation_p1, primitive.primitive_creation_p2, primitive.primitive_creation_p3.

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

## primitive.api_curve_ops

- 26 of 638 cases; first mismatching field: `values.rec.len` 19, `values.rec` 7.
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
  - and 15 more

## primitive.api_primitives

- 198 of 2245 cases; first mismatching field: `values.rec` 143, `values.rec.len` 55.
- Fork commits:
  - c8f6a31 SMLib: sweeps, lofts and extrusions pass on their statuses

## primitive.api_trim_ops

- 10 of 211 cases; first mismatching field: `values.rec` 10.
- Fork commits:
  - 449c81c SMLib: a face made before its making fails keeps its surface
  - 54148db SMLib: draft corners that keep the curves' curvature
  - 4803d6d SMLib: offsets and shells pass on the statuses they were dropping
  - 5d3939e SMLib: a curve sweep says where it placed the profile

## primitive.primitive_creation_p1

- 232 of 773 cases; first mismatching field: `error` 165, `values.curves` 29, `values.geometry.len` 22.
- Fork commits:
  - 9846aa0 SMLib compiles clean under the USD build's warning flags (-Wall -Wshadow=local -Wunused-parameter, -Werror): 19 warnings fixed at their causes
  - 3ae7bfa Draft: a drafted arc's wall is its exact cone, an SmCone in its own parametrization, its seams' pcurves cut at the cone's spans
  - e4f1377 revolve: a curve touching the axis inside itself is split there
  - 183ea65 draft: a wall that closes partway up ends at its apex, its neighbours meeting above
  - fd59111 extrude: a text outline's quadratic Béziers sweep exactly and pass SMLib's check
  - 67b8a1e boolean2d: a plane grows to cover a face however far away it lies
  - 38b4212 revolve: a whole turn of a holed profile leaves the hole's ring a void
  - fe21e32 SMLib: a sheet sweep's profile curves share their meeting vertices
  - 449c81c SMLib: a face made before its making fails keeps its surface
  - e8f23e0 SMLib: the tangent draft corner continues lines' and arcs' walls as their tangent planes
  - 54148db SMLib: draft corners that keep the curves' curvature
  - 9853e54 SMLib: a revolved curve is its analytic surface, per profile curve
  - and 9 more

## primitive.primitive_creation_p2

- 54 of 374 cases; first mismatching field: `values.record` 54.
- Fork commits:
  - 9846aa0 SMLib compiles clean under the USD build's warning flags (-Wall -Wshadow=local -Wunused-parameter, -Werror): 19 warnings fixed at their causes
  - 3ae7bfa Draft: a drafted arc's wall is its exact cone, an SmCone in its own parametrization, its seams' pcurves cut at the cone's spans
  - e4f1377 revolve: a curve touching the axis inside itself is split there
  - 183ea65 draft: a wall that closes partway up ends at its apex, its neighbours meeting above
  - fd59111 extrude: a text outline's quadratic Béziers sweep exactly and pass SMLib's check
  - 67b8a1e boolean2d: a plane grows to cover a face however far away it lies
  - 38b4212 revolve: a whole turn of a holed profile leaves the hole's ring a void
  - fe21e32 SMLib: a sheet sweep's profile curves share their meeting vertices
  - 449c81c SMLib: a face made before its making fails keeps its surface
  - e8f23e0 SMLib: the tangent draft corner continues lines' and arcs' walls as their tangent planes
  - 54148db SMLib: draft corners that keep the curves' curvature
  - 9853e54 SMLib: a revolved curve is its analytic surface, per profile curve
  - and 9 more

## primitive.primitive_creation_p3

- 16 of 330 cases; first mismatching field: `values.record` 15, `values.record.len` 1.
- Fork commits:
  - 9846aa0 SMLib compiles clean under the USD build's warning flags (-Wall -Wshadow=local -Wunused-parameter, -Werror): 19 warnings fixed at their causes
  - 3ae7bfa Draft: a drafted arc's wall is its exact cone, an SmCone in its own parametrization, its seams' pcurves cut at the cone's spans
  - e4f1377 revolve: a curve touching the axis inside itself is split there
  - 183ea65 draft: a wall that closes partway up ends at its apex, its neighbours meeting above
  - fd59111 extrude: a text outline's quadratic Béziers sweep exactly and pass SMLib's check
  - 67b8a1e boolean2d: a plane grows to cover a face however far away it lies
  - 38b4212 revolve: a whole turn of a holed profile leaves the hole's ring a void
  - fe21e32 SMLib: a sheet sweep's profile curves share their meeting vertices
  - 449c81c SMLib: a face made before its making fails keeps its surface
  - e8f23e0 SMLib: the tangent draft corner continues lines' and arcs' walls as their tangent planes
  - 54148db SMLib: draft corners that keep the curves' curvature
  - 9853e54 SMLib: a revolved curve is its analytic surface, per profile curve
  - and 9 more

Measured in #2512 (WSL, the fork's libraries built with `harness/cpp/CMakeLists.txt`'s flags).
