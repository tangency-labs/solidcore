# topology: cases that move with the difftest pin to the SMLib fork (2026-11-03)

Units: topology.api_files, topology.api_projection, topology.assert_valid, topology.brep_constructor, topology.brep_cutting, topology.brep_data, topology.brep_local_operation, topology.brep_mirror, topology.brep_p1, topology.brep_p10, topology.brep_p2, topology.brep_p3, topology.brep_p4, topology.brep_p5, topology.brep_p6, topology.brep_p7, topology.brep_p8, topology.brep_planar_section, topology.brep_validation, topology.edge, topology.edgeuse_p1, topology.edgeuse_p2, topology.edgeuse_p3, topology.face_p1, topology.face_p2, topology.face_p3, topology.face_p4, topology.loopuse, topology.proto_topology_p1, topology.proto_topology_p2, topology.proto_topology_p3, topology.trimming_tools, topology.vertex, topology.vertexuse.

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

## topology.api_files

- 21 of 339 cases; first mismatching field: `ints.rec` 19, `values.record` 2.
- Fork commits:
  - 54148db SMLib: draft corners that keep the curves' curvature
  - 4803d6d SMLib: offsets and shells pass on the statuses they were dropping
  - 5d3939e SMLib: a curve sweep says where it placed the profile

## topology.api_projection

- 7 of 48 cases; first mismatching field: `values.record` 7.
- Fork commits:
  - 54148db SMLib: draft corners that keep the curves' curvature
  - e7711bb SMLib: curve classification and trimming tools pass on their statuses
  - 4803d6d SMLib: offsets and shells pass on the statuses they were dropping
  - 5d3939e SMLib: a curve sweep says where it placed the profile

## topology.assert_valid

- 9 of 965 cases; first mismatching field: `values.edge_test13_max_gap` 3, `values.update_report_values.len` 3, `values.report_values.len` 2.
- Fork commits:
  - edf7414 capi tests: the perforated plate's 2D boolean is gated on the work it counts, not on its time
  - 8d0fbff Edge: an edge with one face on both sides, its radial pair in two loopuses, is not removable
  - d80ea93 Brep: a face shrunk to a bound within rounding of a knot ends on that knot, with no sliver span past it
  - 32fd80a Fillet: the mitre of two equal cylinders is their exact plane section, and its faces keep their exact trim curves
  - 6f811b8 Brep: a face shrunk to its trims keeps them where its surface's trim keeps its parameterization
  - 4ecfbab SmBrep: a vertex where the body touches itself splits into one vertex per fan, and joins again
  - b518bbe merge: a cover crosses a wrapped seam only on the faces' trims, never retaken after the stored one
  - 95bf434 measure: a face's tight box takes its points inside it, its boundary from its edges
  - 345bec1 merge: a cover that spans a wrapped surface's period is taken again on the faces' trims
  - d907cda measure: the face's tight box no longer takes a zone tolerance it does not use
  - 326de73 measure: a face's tight box divides its cells at the surface's knots, so a screw wall's cells settle
  - beb6061 measure: a face's tight box is its trimmed region's, found inside the trim, not the surface's peak
  - and 43 more

## topology.brep_constructor

- 1 of 185 cases; first mismatching field: `ints.ops` 1.
- Fork commits:
  - ae8dd10 SMLib: B-rep data arrays construct and destroy their elements
  - 4ea11a1 SMLib: brep data's file and database I/O pass on their statuses

## topology.brep_cutting

- 26 of 108 cases; first mismatching field: `values.geometry` 21, `values.record` 5.
- Fork commits:
  - 54148db SMLib: draft corners that keep the curves' curvature
  - 4803d6d SMLib: offsets and shells pass on the statuses they were dropping
  - 5d3939e SMLib: a curve sweep says where it placed the profile

## topology.brep_data

- 1681 of 2379 cases; first mismatching field: `ints.ints` 1068, `values.dbls.len` 588, `ints.ints.len` 25.
- Fork commits:
  - ae8dd10 SMLib: B-rep data arrays construct and destroy their elements
  - 4ea11a1 SMLib: brep data's file and database I/O pass on their statuses

## topology.brep_local_operation

- 163 of 165 cases; first mismatching field: `values.edge_curve_owners.len` 77, `values.geometry` 75, `values.geometry.len` 11.
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

## topology.brep_mirror

- 145 of 341 cases; first mismatching field: `values.geometry` 137, `values.record` 7, `values.record.len` 1.
- Fork commits:
  - edf7414 capi tests: the perforated plate's 2D boolean is gated on the work it counts, not on its time
  - 9846aa0 SMLib compiles clean under the USD build's warning flags (-Wall -Wshadow=local -Wunused-parameter, -Werror): 19 warnings fixed at their causes
  - 0c00de7 Surfaces: an analytic surface turned round by SwapUV keeps its definition as made, not re-derived from its net
  - 32fd80a Fillet: the mitre of two equal cylinders is their exact plane section, and its faces keep their exact trim curves
  - 5c275e6 Fillet: build-then-subtract overflow is gone
  - 43fd6b6 offset surface: a surface of revolution whose generator is a circle continues it as its own circle
  - 4ecfbab SmBrep: a vertex where the body touches itself splits into one vertex per fan, and joins again
  - db4a655 fillets: a cliff overflow past a closed rim's wall closes on the wall continued as itself
  - b298cc2 offset surface: a surface of revolution's extension closes its turn and wraps past it
  - 93db777 offset surface: a cylinder, cone, sphere or torus continues as its own circles, span by span
  - de0a6b8 offset surface: an analytic face continues as itself past its domain
  - 829b982 brep: a topological vertex between a line and a seam iso curve is joined, not a crash
  - and 30 more

## topology.brep_p1

- 26 of 2636 cases; first mismatching field: `values.out` 20, `values.out.len` 4, `values.intersection.points.len` 2.
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
  - and 15 more

## topology.brep_p10

- 110 of 1416 cases; first mismatching field: `values.dbls.len` 51, `values.dbls` 46, `ints.ints` 7.
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

## topology.brep_p2

- 14 of 1720 cases; first mismatching field: `values.out` 13, `values.out.len` 1.
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

## topology.brep_p3

- 6 of 1358 cases; first mismatching field: `values.out.len` 6.
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

## topology.brep_p4

- 8 of 327 cases; first mismatching field: `values.faces_values` 5, `values.edges_values` 3.
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

## topology.brep_p5

- 7 of 797 cases; first mismatching field: `values.record` 4, `values.geometry` 3.
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

## topology.brep_p6

- 123 of 595 cases; first mismatching field: `ints.state` 99, `values.geom` 24.
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

## topology.brep_p7

- 13 of 462 cases; first mismatching field: `values.out.len` 12, `values.out` 1.
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

## topology.brep_p8

- 4 of 355 cases; first mismatching field: `values.out.len` 3, `values.out` 1.
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

## topology.brep_planar_section

- 24 of 167 cases; first mismatching field: `values.record` 14, `values.record.len` 10.
- Fork commits:
  - edf7414 capi tests: the perforated plate's 2D boolean is gated on the work it counts, not on its time
  - dda881f Surface: a cover steps across the seam a surface wraps in STEP u, where a swapped surface's NURBS closes in v; the drafted peg's cone walls are one face
  - 8d0fbff Edge: an edge with one face on both sides, its radial pair in two loopuses, is not removable
  - f8e6709 Fillet: every Ball's face is the rational patch on the ball's sphere, carrying that sphere as its definition; the SmSphere path is gone
  - 32fd80a Fillet: the mitre of two equal cylinders is their exact plane section, and its faces keep their exact trim curves
  - 014f66f SMLib: SmProceduralBlend, a blend surface stored as its definition, its net the approximation
  - 4ecfbab SmBrep: a vertex where the body touches itself splits into one vertex per fan, and joins again
  - b518bbe merge: a cover crosses a wrapped seam only on the faces' trims, never retaken after the stored one
  - 345bec1 merge: a cover that spans a wrapped surface's period is taken again on the faces' trims
  - eaa998a surface: a dropped point's foot from a guess is one the local solve converged to
  - 737c392 surface: a dropped point's foot is off along the surface by no more than the face's tolerance
  - 2e8758b fillets: an edge is smooth within half a degree, not two
  - and 36 more

## topology.brep_validation

- 2 of 112 cases; first mismatching field: `ints.reports` 2.
- Fork commits:
  - edf7414 capi tests: the perforated plate's 2D boolean is gated on the work it counts, not on its time
  - 8d0fbff Edge: an edge with one face on both sides, its radial pair in two loopuses, is not removable
  - 32fd80a Fillet: the mitre of two equal cylinders is their exact plane section, and its faces keep their exact trim curves
  - 4ecfbab SmBrep: a vertex where the body touches itself splits into one vertex per fan, and joins again
  - b518bbe merge: a cover crosses a wrapped seam only on the faces' trims, never retaken after the stored one
  - 345bec1 merge: a cover that spans a wrapped surface's period is taken again on the faces' trims
  - 2e8758b fillets: an edge is smooth within half a degree, not two
  - 829b982 brep: a topological vertex between a line and a seam iso curve is joined, not a crash
  - 2e9e6b2 capi tests: the Menger row join counts the work it does instead of timing it
  - be8ac73 booleans: a new edge from a faceuse shell to a wire keeps the faceuse shell
  - 69c5198 booleans: a new face's closure test searches from both of its sides
  - b7ec55b fillet: a knife edge rounds to the ball touching both faces back from it
  - and 22 more

## topology.edge

- 49 of 2949 cases; first mismatching field: `values.state` 49.
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

## topology.edgeuse_p1

- 52 of 5586 cases; first mismatching field: `values.state` 52.
- Fork commits:
  - e36feb1 SMLib UV trims: a drop that comes back in pieces is joined where they meet, never stretched onto the edge (solidcore#1383)
  - 02e0dde SMLib, capi: the fuzzers' sanitizer findings fixed at their causes
  - 4038552 SMLib: shells, vertexuses, loopuses, faceuses and the traverser pass on statuses
  - 666242b SMLib: edgeuses pass on statuses; their checks answer only from what they computed
  - ef9041d SMLib: the validity check takes time linear in the body, not quadratic

## topology.edgeuse_p2

- 3 of 4795 cases; first mismatching field: `values.state` 3.
- Fork commits:
  - e36feb1 SMLib UV trims: a drop that comes back in pieces is joined where they meet, never stretched onto the edge (solidcore#1383)
  - 666242b SMLib: edgeuses pass on statuses; their checks answer only from what they computed
  - ef9041d SMLib: the validity check takes time linear in the body, not quadratic

## topology.edgeuse_p3

- 138 of 5611 cases; first mismatching field: `values.state` 122, `values.state.len` 16.
- Fork commits:
  - e36feb1 SMLib UV trims: a drop that comes back in pieces is joined where they meet, never stretched onto the edge (solidcore#1383)
  - 666242b SMLib: edgeuses pass on statuses; their checks answer only from what they computed
  - ef9041d SMLib: the validity check takes time linear in the body, not quadratic

## topology.face_p1

- 349 of 2569 cases; first mismatching field: `values.out` 349.
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

## topology.face_p2

- 3 of 2692 cases; first mismatching field: `values.state` 2, `values.state.len` 1.
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

## topology.face_p3

- 74 of 2407 cases; first mismatching field: `values.out.len` 42, `values.out` 32.
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

## topology.face_p4

- 516 of 5930 cases; first mismatching field: `values.out` 511, `values.out.len` 5.
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

## topology.loopuse

- 2 of 9013 cases; first mismatching field: `values.out` 2.
- Fork commits:
  - 4038552 SMLib: shells, vertexuses, loopuses, faceuses and the traverser pass on statuses
  - 7a35f33 SMLib: loops and vertices pass on statuses; their queries read only what they computed
  - ef9041d SMLib: the validity check takes time linear in the body, not quadratic

## topology.proto_topology_p1

- 123 of 1922 cases; first mismatching field: `values.out` 118, `values.out.len` 5.
- Fork commits:
  - bec8a0b imprint: an analytic intersection curve continued to its vertex stays a line, circle or ellipse
  - 81cd61a boolean: an intersection edge along a face's own edge is that edge in the other body
  - 9eec159 boolean: an intersection curve inside one ProtoVertex is split at the vertices on it
  - f11d0ac boolean: a coincident edge definition across a split point is split with its edge
  - 432e589 SMLib: a crossing re-places a vertex only when its points meet within tolerance
  - c9db18e SMLib: an edge/face vertex is placed by the crossing that fixes it best
  - 03e7bf1 SMLib: an intersection curve inside one vertex is that vertex, not a closed edge on it
  - 02e0dde SMLib, capi: the fuzzers' sanitizer findings fixed at their causes
  - 4ea11a1 SMLib: brep data's file and database I/O pass on their statuses
  - 72ab5f4 SMLib: a definition removed on an error path is freed, one still listed is not
  - 3545008 SMLib: a split ProtoEdge's upper EdgeDefinition starts without an owner
  - bcc6f01 SMLib: surface intersection and proto-topology pass on their statuses
  - and 8 more

## topology.proto_topology_p2

- 187 of 1767 cases; first mismatching field: `values.out.len` 98, `values.out` 89.
- Fork commits:
  - bec8a0b imprint: an analytic intersection curve continued to its vertex stays a line, circle or ellipse
  - 81cd61a boolean: an intersection edge along a face's own edge is that edge in the other body
  - 9eec159 boolean: an intersection curve inside one ProtoVertex is split at the vertices on it
  - f11d0ac boolean: a coincident edge definition across a split point is split with its edge
  - 432e589 SMLib: a crossing re-places a vertex only when its points meet within tolerance
  - c9db18e SMLib: an edge/face vertex is placed by the crossing that fixes it best
  - 03e7bf1 SMLib: an intersection curve inside one vertex is that vertex, not a closed edge on it
  - 72ab5f4 SMLib: a definition removed on an error path is freed, one still listed is not
  - 3545008 SMLib: a split ProtoEdge's upper EdgeDefinition starts without an owner
  - bcc6f01 SMLib: surface intersection and proto-topology pass on their statuses
  - 0b1b945 SMLib: an intersection edge's end goes to its own vertex when its ProtoVertex is several
  - bf5d3ae SMLib: an intersection edge that stops short of its vertex is continued to it
  - and 5 more

## topology.proto_topology_p3

- 148 of 3013 cases; first mismatching field: `values.out.len` 79, `values.out` 69.
- Fork commits:
  - bec8a0b imprint: an analytic intersection curve continued to its vertex stays a line, circle or ellipse
  - 81cd61a boolean: an intersection edge along a face's own edge is that edge in the other body
  - 9eec159 boolean: an intersection curve inside one ProtoVertex is split at the vertices on it
  - f11d0ac boolean: a coincident edge definition across a split point is split with its edge
  - 432e589 SMLib: a crossing re-places a vertex only when its points meet within tolerance
  - c9db18e SMLib: an edge/face vertex is placed by the crossing that fixes it best
  - 03e7bf1 SMLib: an intersection curve inside one vertex is that vertex, not a closed edge on it
  - 72ab5f4 SMLib: a definition removed on an error path is freed, one still listed is not
  - 3545008 SMLib: a split ProtoEdge's upper EdgeDefinition starts without an owner
  - bcc6f01 SMLib: surface intersection and proto-topology pass on their statuses
  - 0b1b945 SMLib: an intersection edge's end goes to its own vertex when its ProtoVertex is several
  - bf5d3ae SMLib: an intersection edge that stops short of its vertex is continued to it
  - and 5 more

## topology.trimming_tools

- 41 of 3305 cases; first mismatching field: `values.out` 33, `values.out.len` 8.
- Fork commits:
  - b0ee3f9 loops: where loops touch at a vertex, turn about one normal for the whole list
  - 9209e76 loops: where loops touch at a vertex, go on round the vertex
  - 1f2b0c4 SMLib: a hole touching its outline at a corner is one loop of the face
  - 7a6b4ad SMLib: a seam edgeuse with no UV trim is not read as the loop jumping the seam
  - e7711bb SMLib: curve classification and trimming tools pass on their statuses
  - f72c79b SMLib: CheckLoop measures each edge end against the vertex the edges meet at

## topology.vertex

- 14 of 5954 cases; first mismatching field: `values.state` 14.
- Fork commits:
  - 02e0dde SMLib, capi: the fuzzers' sanitizer findings fixed at their causes
  - 7a35f33 SMLib: loops and vertices pass on statuses; their queries read only what they computed

## topology.vertexuse

- 6 of 2686 cases; first mismatching field: `values.state` 6.
- Fork commits:
  - e7711bb SMLib: curve classification and trimming tools pass on their statuses
  - 4038552 SMLib: shells, vertexuses, loopuses, faceuses and the traverser pass on statuses

Measured in #2512 (WSL, the fork's libraries built with `harness/cpp/CMakeLists.txt`'s flags).
