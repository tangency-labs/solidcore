# intersect: cases that move with the difftest pin to the SMLib fork (2026-11-03)

Units: intersect.adv_surface_intersector, intersect.api_intersectors, intersect.curve_class_p2, intersect.curve_class_p3, intersect.curve_class_p4, intersect.curve_class_p5, intersect.curve_class_p7, intersect.surface_drop_curve, intersect.surface_drop_lib_p1, intersect.surface_drop_lib_p2, intersect.surface_drop_lib_p3, intersect.surface_intersector_p1, intersect.surface_intersector_p2, intersect.surface_intersector_p3, intersect.surface_intersector_p4, intersect.surface_intersector_p5, intersect.topology_solver_p1, intersect.topology_solver_p2.

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

## intersect.adv_surface_intersector

- 30 of 1949 cases; first mismatching field: `values.rec` 28, `values.rec.len` 2.
- Fork commits:
  - d37ca80 SSI: a curve starting at a singular point starts along the inside direction nearest the given one
  - e7ebc8d SMLib: booleans' merge, stitch and topology intersection pass on statuses
  - 1f5c0f7 SSI: trace small closed loops above tolerance instead of dropping them
  - ef4ab25 SSI: trace tangential intersection curves
  - 9e8bf10 SSI: find isolated tangential contact points inside both surfaces

## intersect.api_intersectors

- 20 of 477 cases; first mismatching field: `values.record` 20.
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
  - and 75 more

## intersect.curve_class_p2

- 27 of 1671 cases; first mismatching field: `values.out` 26, `values.out.len` 1.
- Fork commits:
  - 44e5e30 SMLib: an interval running along a face edge between two points on it classifies to that edge
  - 02e0dde SMLib, capi: the fuzzers' sanitizer findings fixed at their causes
  - e7711bb SMLib: curve classification and trimming tools pass on their statuses

## intersect.curve_class_p3

- 3 of 3708 cases; first mismatching field: `values.out.len` 2, `values.out` 1.
- Fork commits:
  - 44e5e30 SMLib: an interval running along a face edge between two points on it classifies to that edge
  - 02e0dde SMLib, capi: the fuzzers' sanitizer findings fixed at their causes
  - e7711bb SMLib: curve classification and trimming tools pass on their statuses

## intersect.curve_class_p4

- 54 of 5410 cases; first mismatching field: `values.out` 53, `values.out.len` 1.
- Fork commits:
  - 44e5e30 SMLib: an interval running along a face edge between two points on it classifies to that edge
  - 02e0dde SMLib, capi: the fuzzers' sanitizer findings fixed at their causes
  - e7711bb SMLib: curve classification and trimming tools pass on their statuses

## intersect.curve_class_p5

- 19 of 1411 cases; first mismatching field: `values.out.len` 19.
- Fork commits:
  - 44e5e30 SMLib: an interval running along a face edge between two points on it classifies to that edge
  - 02e0dde SMLib, capi: the fuzzers' sanitizer findings fixed at their causes
  - e7711bb SMLib: curve classification and trimming tools pass on their statuses

## intersect.curve_class_p7

- 4 of 1673 cases; first mismatching field: `values.out` 4.
- Fork commits:
  - 44e5e30 SMLib: an interval running along a face edge between two points on it classifies to that edge
  - 02e0dde SMLib, capi: the fuzzers' sanitizer findings fixed at their causes
  - e7711bb SMLib: curve classification and trimming tools pass on their statuses

## intersect.surface_drop_curve

- 5 of 547 cases; first mismatching field: `values.rec.len` 5.
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
  - and 30 more

## intersect.surface_drop_lib_p1

- 4 of 2926 cases; first mismatching field: `values.rec` 2, `values.rec.len` 2.
- Fork commits:
  - e149e05 SMLib: point solves, drops and the curve and surface caches pass on statuses

## intersect.surface_drop_lib_p2

- 69 of 1531 cases; first mismatching field: `values.rec` 69.
- Fork commits:
  - e149e05 SMLib: point solves, drops and the curve and surface caches pass on statuses

## intersect.surface_drop_lib_p3

- 569 of 1532 cases; first mismatching field: `values.rec.len` 379, `values.rec` 190.
- Fork commits:
  - e149e05 SMLib: point solves, drops and the curve and surface caches pass on statuses

## intersect.surface_intersector_p1

- 3 of 1510 cases; first mismatching field: `values.rec` 2, `values.rec.len` 1.
- Fork commits:
  - 8d95600 intersections: a curve traced into a surface pole keeps the pole's tangent
  - bcc6f01 SMLib: surface intersection and proto-topology pass on their statuses
  - e9fde8f SSI: trace crossing branches through tangent contact nodes
  - ef4ab25 SSI: trace tangential intersection curves
  - 1ff9264 SSI: split intersection curves at singular crossing points
  - 920da20 SSI: trace each tangent curve once (no duplicate/overlapping traces)
  - 0dbe681 SSI: don't reverse-trace an empty curve (SM_ERR_NULL_POINTER)

## intersect.surface_intersector_p2

- 61 of 1686 cases; first mismatching field: `values.rec.len` 41, `values.rec` 20.
- Fork commits:
  - 8d95600 intersections: a curve traced into a surface pole keeps the pole's tangent
  - bcc6f01 SMLib: surface intersection and proto-topology pass on their statuses
  - e9fde8f SSI: trace crossing branches through tangent contact nodes
  - ef4ab25 SSI: trace tangential intersection curves
  - 1ff9264 SSI: split intersection curves at singular crossing points
  - 920da20 SSI: trace each tangent curve once (no duplicate/overlapping traces)
  - 0dbe681 SSI: don't reverse-trace an empty curve (SM_ERR_NULL_POINTER)

## intersect.surface_intersector_p3

- 42 of 1527 cases; first mismatching field: `values.rec` 38, `values.rec.len` 4.
- Fork commits:
  - 8d95600 intersections: a curve traced into a surface pole keeps the pole's tangent
  - bcc6f01 SMLib: surface intersection and proto-topology pass on their statuses
  - e9fde8f SSI: trace crossing branches through tangent contact nodes
  - ef4ab25 SSI: trace tangential intersection curves
  - 1ff9264 SSI: split intersection curves at singular crossing points
  - 920da20 SSI: trace each tangent curve once (no duplicate/overlapping traces)
  - 0dbe681 SSI: don't reverse-trace an empty curve (SM_ERR_NULL_POINTER)

## intersect.surface_intersector_p4

- 12 of 1542 cases; first mismatching field: `values.rec.len` 12.
- Fork commits:
  - 8d95600 intersections: a curve traced into a surface pole keeps the pole's tangent
  - bcc6f01 SMLib: surface intersection and proto-topology pass on their statuses
  - e9fde8f SSI: trace crossing branches through tangent contact nodes
  - ef4ab25 SSI: trace tangential intersection curves
  - 1ff9264 SSI: split intersection curves at singular crossing points
  - 920da20 SSI: trace each tangent curve once (no duplicate/overlapping traces)
  - 0dbe681 SSI: don't reverse-trace an empty curve (SM_ERR_NULL_POINTER)

## intersect.surface_intersector_p5

- 26 of 1731 cases; first mismatching field: `values.rec.len` 19, `values.rec` 7.
- Fork commits:
  - 8d95600 intersections: a curve traced into a surface pole keeps the pole's tangent
  - bcc6f01 SMLib: surface intersection and proto-topology pass on their statuses
  - e9fde8f SSI: trace crossing branches through tangent contact nodes
  - ef4ab25 SSI: trace tangential intersection curves
  - 1ff9264 SSI: split intersection curves at singular crossing points
  - 920da20 SSI: trace each tangent curve once (no duplicate/overlapping traces)
  - 0dbe681 SSI: don't reverse-trace an empty curve (SM_ERR_NULL_POINTER)

## intersect.topology_solver_p1

- 15 of 3592 cases; first mismatching field: `values.out` 14, `ints.status` 1.
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
  - and 22 more

## intersect.topology_solver_p2

- 1 of 2575 cases; first mismatching field: `values.static` 1.
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
  - and 7 more

Measured in #2512 (WSL, the fork's libraries built with `harness/cpp/CMakeLists.txt`'s flags).
