# query: cases that move with the difftest pin to the SMLib fork (2026-11-03)

Units: query.api_queries, query.brep, query.edge, query.face_p1, query.face_p2, query.face_p3, query.face_p4, query.poly, query.tess.

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

## query.api_queries

- 1 of 1399 cases; first mismatching field: `values.geometry` 1.
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
  - and 37 more

## query.brep

- 22 of 1172 cases; first mismatching field: `values.steps` 21, `values.classify` 1.
- Fork commits:
  - edf7414 capi tests: the perforated plate's 2D boolean is gated on the work it counts, not on its time
  - 32fd80a Fillet: the mitre of two equal cylinders is their exact plane section, and its faces keep their exact trim curves
  - 4ecfbab SmBrep: a vertex where the body touches itself splits into one vertex per fan, and joins again
  - 829b982 brep: a topological vertex between a line and a seam iso curve is joined, not a crash
  - 2e9e6b2 capi tests: the Menger row join counts the work it does instead of timing it
  - be8ac73 booleans: a new edge from a faceuse shell to a wire keeps the faceuse shell
  - 69c5198 booleans: a new face's closure test searches from both of its sides
  - 60b1c6a boolean2d: faces apart are not intersected, and loops apart are not classified
  - 972604a SMLib: a scaled body's model size is its box, vertices, edges and faces
  - 2670650 SMLib: a rotation or reflection keeps every tolerance; a scale scales them
  - a1ee273 offset faces: neighbours are extended to reach a face moved further than the body's size
  - 1072d1b Fillet: a corner's patch meets a fillet on the fillet's own edge (solidcore#1377)
  - and 14 more

## query.edge

- 25 of 1368 cases; first mismatching field: `values.box` 19, `values.normal` 6.
- Fork commits:
  - 8d0fbff Edge: an edge with one face on both sides, its radial pair in two loopuses, is not removable
  - b518bbe merge: a cover crosses a wrapped seam only on the faces' trims, never retaken after the stored one
  - 345bec1 merge: a cover that spans a wrapped surface's period is taken again on the faces' trims
  - 2e8758b fillets: an edge is smooth within half a degree, not two
  - b7ec55b fillet: a knife edge rounds to the ball touching both faces back from it
  - c343033 SMLib: a UV trim's turn at its ends is judged on its surface; the check no longer kinks a circle's trim
  - 02e0dde SMLib, capi: the fuzzers' sanitizer findings fixed at their causes
  - 47507b6 SMLib: edges pass on statuses; reversing one keeps it whole on failure
  - 5533bbb SMLib: the check orders faces round a tangent edge by where they go

## query.face_p1

- 24 of 2266 cases; first mismatching field: `values.box` 22, `values.pc` 2.
- Fork commits:
  - d80ea93 Brep: a face shrunk to a bound within rounding of a knot ends on that knot, with no sliver span past it
  - 6f811b8 Brep: a face shrunk to its trims keeps them where its surface's trim keeps its parameterization
  - 95bf434 measure: a face's tight box takes its points inside it, its boundary from its edges
  - d907cda measure: the face's tight box no longer takes a zone tolerance it does not use
  - 326de73 measure: a face's tight box divides its cells at the surface's knots, so a screw wall's cells settle
  - beb6061 measure: a face's tight box is its trimmed region's, found inside the trim, not the surface's peak
  - 54ca3e0 face, check: a loop far smaller than its face is oriented, and its trim measured, on the surface's scale
  - 55bffe0 face check: a seam that collapses to a point is not classified
  - fcd5ea5 SMLib: no preprocessor directive inside a macro's argument
  - 44e5e30 SMLib: an interval running along a face edge between two points on it classifies to that edge
  - 02e0dde SMLib, capi: the fuzzers' sanitizer findings fixed at their causes
  - e5606fb SmFace: closing a corner's UV gap keeps the trims on their edges' parameters
  - and 7 more

## query.face_p2

- 14 of 1401 cases; first mismatching field: `values.pc` 14.
- Fork commits:
  - d80ea93 Brep: a face shrunk to a bound within rounding of a knot ends on that knot, with no sliver span past it
  - 6f811b8 Brep: a face shrunk to its trims keeps them where its surface's trim keeps its parameterization
  - 95bf434 measure: a face's tight box takes its points inside it, its boundary from its edges
  - d907cda measure: the face's tight box no longer takes a zone tolerance it does not use
  - 326de73 measure: a face's tight box divides its cells at the surface's knots, so a screw wall's cells settle
  - beb6061 measure: a face's tight box is its trimmed region's, found inside the trim, not the surface's peak
  - bec8a0b imprint: an analytic intersection curve continued to its vertex stays a line, circle or ellipse
  - 54ca3e0 face, check: a loop far smaller than its face is oriented, and its trim measured, on the surface's scale
  - 81cd61a boolean: an intersection edge along a face's own edge is that edge in the other body
  - 9eec159 boolean: an intersection curve inside one ProtoVertex is split at the vertices on it
  - 55bffe0 face check: a seam that collapses to a point is not classified
  - f11d0ac boolean: a coincident edge definition across a split point is split with its edge
  - and 23 more

## query.face_p3

- 20 of 1199 cases; first mismatching field: `values.min` 10, `values.max` 7, `values.point` 1.
- Fork commits:
  - d80ea93 Brep: a face shrunk to a bound within rounding of a knot ends on that knot, with no sliver span past it
  - 6f811b8 Brep: a face shrunk to its trims keeps them where its surface's trim keeps its parameterization
  - 95bf434 measure: a face's tight box takes its points inside it, its boundary from its edges
  - d907cda measure: the face's tight box no longer takes a zone tolerance it does not use
  - 326de73 measure: a face's tight box divides its cells at the surface's knots, so a screw wall's cells settle
  - beb6061 measure: a face's tight box is its trimmed region's, found inside the trim, not the surface's peak
  - 54ca3e0 face, check: a loop far smaller than its face is oriented, and its trim measured, on the surface's scale
  - 55bffe0 face check: a seam that collapses to a point is not classified
  - fcd5ea5 SMLib: no preprocessor directive inside a macro's argument
  - 44e5e30 SMLib: an interval running along a face edge between two points on it classifies to that edge
  - 02e0dde SMLib, capi: the fuzzers' sanitizer findings fixed at their causes
  - e5606fb SmFace: closing a corner's UV gap keeps the trims on their edges' parameters
  - and 6 more

## query.face_p4

- 4 of 503 cases; first mismatching field: `error` 4.
- Fork commits:
  - d80ea93 Brep: a face shrunk to a bound within rounding of a knot ends on that knot, with no sliver span past it
  - 6f811b8 Brep: a face shrunk to its trims keeps them where its surface's trim keeps its parameterization
  - 95bf434 measure: a face's tight box takes its points inside it, its boundary from its edges
  - d907cda measure: the face's tight box no longer takes a zone tolerance it does not use
  - 326de73 measure: a face's tight box divides its cells at the surface's knots, so a screw wall's cells settle
  - beb6061 measure: a face's tight box is its trimmed region's, found inside the trim, not the surface's peak
  - 54ca3e0 face, check: a loop far smaller than its face is oriented, and its trim measured, on the surface's scale
  - 55bffe0 face check: a seam that collapses to a point is not classified
  - fcd5ea5 SMLib: no preprocessor directive inside a macro's argument
  - 44e5e30 SMLib: an interval running along a face edge between two points on it classifies to that edge
  - 02e0dde SMLib, capi: the fuzzers' sanitizer findings fixed at their causes
  - e5606fb SmFace: closing a corner's UV gap keeps the trims on their edges' parameters
  - and 6 more

## query.poly

- 41 of 966 cases; first mismatching field: `values.bbox` 41.
- Fork commits:
  - 3aae1c4 fillet: the chord and the circle are the section generator's orders 0 and 1
  - ae8dd10 SMLib: B-rep data arrays construct and destroy their elements
  - 02e0dde SMLib, capi: the fuzzers' sanitizer findings fixed at their causes
  - ca7d1f5 capi: no shim drops a status SMLib returns
  - e309858 SMLib: tessellated triangles stay within the chord tolerance

## query.tess

- 1 of 676 cases; first mismatching field: `values.node_geom.len` 1.
- Fork commits:
  - bb9f45e SMLib: the tessellator meshes a round by its definition on its net (SmBlendNetScope)
  - e309858 SMLib: tessellated triangles stay within the chord tolerance
  - f98bdf5 build: macOS (Apple Silicon) builds and tests with the documented commands

Measured in #2512 (WSL, the fork's libraries built with `harness/cpp/CMakeLists.txt`'s flags).
