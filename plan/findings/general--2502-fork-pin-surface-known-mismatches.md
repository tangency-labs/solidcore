# surface: cases that move with the difftest pin to the SMLib fork (2026-11-03)

Units: surface.cone_p1, surface.cone_p2, surface.cone_p3, surface.offset_surface, surface.plane_p1, surface.plane_p2, surface.sphere, surface.surf_of_extrusion_p1, surface.surf_of_extrusion_p2, surface.surf_of_revolution_p1, surface.surf_of_revolution_p2, surface.surf_of_revolution_p3, surface.surface_approximate, surface.surface_assert_valid, surface.surface_cache_p1, surface.surface_p1, surface.surface_p10, surface.surface_p11, surface.surface_p2, surface.surface_p3, surface.surface_p4, surface.surface_p5, surface.surface_p6, surface.surface_p7, surface.surface_p8, surface.surface_point_solve_p2, surface.surface_point_solve_p3, surface.torus, surface.volume_p2.

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

## surface.cone_p1

- 3232 of 20736 cases; first mismatching field: `values.status` 1483, `values.out` 696, `values.build` 631.
- Fork commits:
  - c2f1788 SSI: a plane sunk a hair into a bare cylinder crosses it in two lines
  - b22bfbc boolean: a plane within tolerance of a cylinder wall touches it along one line
  - 811ab6e SMLib: a cylinder found in an inside-out NURBS sweeps counter-clockwise about its STEP axis
  - 758829c SMLib: an offset that collapses onto its spine is told apart from a missing one
  - dfc7cb6 SMLib: analytic STEPInversion answers apart from its failures
  - 1279539 SMLib: a cone's iso circle through a point is at the point's height
  - 9fd7e7e SMLib: a cone of semi-angle 0 reports its missing apex as an answer, not a failure
  - 5df5b27 SMLib: an SmCylinder rebuilt from its NURBS keeps one radius
  - 62762b8 SMLib: parallel cylinders intersect along their faces, not the whole cylinders
  - 672a437 SMLib: the analytic surfaces and curves pass on their statuses

## surface.cone_p2

- 505 of 1713 cases; first mismatching field: `values.out` 477, `values.build` 28.
- Fork commits:
  - c2f1788 SSI: a plane sunk a hair into a bare cylinder crosses it in two lines
  - b22bfbc boolean: a plane within tolerance of a cylinder wall touches it along one line
  - 758829c SMLib: an offset that collapses onto its spine is told apart from a missing one
  - dfc7cb6 SMLib: analytic STEPInversion answers apart from its failures
  - 1279539 SMLib: a cone's iso circle through a point is at the point's height
  - 9fd7e7e SMLib: a cone of semi-angle 0 reports its missing apex as an answer, not a failure
  - 5df5b27 SMLib: an SmCylinder rebuilt from its NURBS keeps one radius
  - 62762b8 SMLib: parallel cylinders intersect along their faces, not the whole cylinders
  - 672a437 SMLib: the analytic surfaces and curves pass on their statuses

## surface.cone_p3

- 165 of 1490 cases; first mismatching field: `values.geometry.len` 140, `error` 20, `values.geometry` 5.
- Fork commits:
  - c2f1788 SSI: a plane sunk a hair into a bare cylinder crosses it in two lines
  - b22bfbc boolean: a plane within tolerance of a cylinder wall touches it along one line
  - 758829c SMLib: an offset that collapses onto its spine is told apart from a missing one
  - dfc7cb6 SMLib: analytic STEPInversion answers apart from its failures
  - 1279539 SMLib: a cone's iso circle through a point is at the point's height
  - 9fd7e7e SMLib: a cone of semi-angle 0 reports its missing apex as an answer, not a failure
  - 5df5b27 SMLib: an SmCylinder rebuilt from its NURBS keeps one radius
  - 62762b8 SMLib: parallel cylinders intersect along their faces, not the whole cylinders
  - 672a437 SMLib: the analytic surfaces and curves pass on their statuses

## surface.offset_surface

- 225 of 3182 cases; first mismatching field: `values.out` 190, `values.out.len` 35.
- Fork commits:
  - 5c275e6 Fillet: build-then-subtract overflow is gone
  - 43fd6b6 offset surface: a surface of revolution whose generator is a circle continues it as its own circle
  - db4a655 fillets: a cliff overflow past a closed rim's wall closes on the wall continued as itself
  - b298cc2 offset surface: a surface of revolution's extension closes its turn and wraps past it
  - 93db777 offset surface: a cylinder, cone, sphere or torus continues as its own circles, span by span
  - de0a6b8 offset surface: an analytic face continues as itself past its domain
  - 4803d6d SMLib: offsets and shells pass on the statuses they were dropping
  - 9f54d05 SMLib: an offset surface extends the same way whatever was evaluated first

## surface.plane_p1

- 162 of 2706 cases; first mismatching field: `values.rec.len` 149, `values.rec` 13.
- Fork commits:
  - 0c00de7 Surfaces: an analytic surface turned round by SwapUV keeps its definition as made, not re-derived from its net
  - 6f811b8 Brep: a face shrunk to its trims keeps them where its surface's trim keeps its parameterization
  - 0286682 SMLib: a plane recognized from a B-spline places its origin along its own axes
  - 02e0dde SMLib, capi: the fuzzers' sanitizer findings fixed at their causes
  - 672a437 SMLib: the analytic surfaces and curves pass on their statuses
  - 0df62c8 SMLib: a box whose faces cannot be made fails instead of leaking its planes

## surface.plane_p2

- 55 of 2228 cases; first mismatching field: `values.rec.len` 48, `values.rec` 7.
- Fork commits:
  - 0286682 SMLib: a plane recognized from a B-spline places its origin along its own axes
  - 02e0dde SMLib, capi: the fuzzers' sanitizer findings fixed at their causes
  - 672a437 SMLib: the analytic surfaces and curves pass on their statuses
  - 0df62c8 SMLib: a box whose faces cannot be made fails instead of leaking its planes

## surface.sphere

- 21 of 816 cases; first mismatching field: `values.out` 10, `values.status` 10, `values.out.len` 1.
- Fork commits:
  - 7ea6a96 SMLib: a torus or sphere found in a NURBS keeps that NURBS; a half meridian measured to its rounding is a half meridian
  - 758829c SMLib: an offset that collapses onto its spine is told apart from a missing one
  - dfc7cb6 SMLib: analytic STEPInversion answers apart from its failures
  - 672a437 SMLib: the analytic surfaces and curves pass on their statuses
  - 1f5c0f7 SSI: trace small closed loops above tolerance instead of dropping them

## surface.surf_of_extrusion_p1

- 16 of 2071 cases; first mismatching field: `values.build` 16.
- Fork commits:
  - 6f811b8 Brep: a face shrunk to its trims keeps them where its surface's trim keeps its parameterization
  - fd59111 extrude: a text outline's quadratic Béziers sweep exactly and pass SMLib's check
  - 02e0dde SMLib, capi: the fuzzers' sanitizer findings fixed at their causes
  - 672a437 SMLib: the analytic surfaces and curves pass on their statuses

## surface.surf_of_extrusion_p2

- 25 of 1548 cases; first mismatching field: `values.build` 15, `values.out` 5, `values.out.len` 5.
- Fork commits:
  - fd59111 extrude: a text outline's quadratic Béziers sweep exactly and pass SMLib's check
  - 02e0dde SMLib, capi: the fuzzers' sanitizer findings fixed at their causes
  - 672a437 SMLib: the analytic surfaces and curves pass on their statuses

## surface.surf_of_revolution_p1

- 3 of 895 cases; first mismatching field: `values.out` 3.
- Fork commits:
  - 0c00de7 Surfaces: an analytic surface turned round by SwapUV keeps its definition as made, not re-derived from its net
  - 6f811b8 Brep: a face shrunk to its trims keeps them where its surface's trim keeps its parameterization
  - 8b2090f shell: a round of revolution moved in past its radius has no offset
  - dfc7cb6 SMLib: analytic STEPInversion answers apart from its failures
  - cc33bc3 SMLib: a revolved face's offset is trimmed at the axis on either side of the frame's origin
  - 02e0dde SMLib, capi: the fuzzers' sanitizer findings fixed at their causes
  - 672a437 SMLib: the analytic surfaces and curves pass on their statuses
  - 0cd9d08 SMLib: the torus check tests its minor interval

## surface.surf_of_revolution_p2

- 138 of 1167 cases; first mismatching field: `values.status` 106, `values.out.len` 24, `values.out` 8.
- Fork commits:
  - 6f811b8 Brep: a face shrunk to its trims keeps them where its surface's trim keeps its parameterization
  - 8b2090f shell: a round of revolution moved in past its radius has no offset
  - dfc7cb6 SMLib: analytic STEPInversion answers apart from its failures
  - cc33bc3 SMLib: a revolved face's offset is trimmed at the axis on either side of the frame's origin
  - 02e0dde SMLib, capi: the fuzzers' sanitizer findings fixed at their causes
  - 672a437 SMLib: the analytic surfaces and curves pass on their statuses
  - 0cd9d08 SMLib: the torus check tests its minor interval

## surface.surf_of_revolution_p3

- 204 of 26731 cases; first mismatching field: `values.out` 168, `values.out.len` 36.
- Fork commits:
  - 6f811b8 Brep: a face shrunk to its trims keeps them where its surface's trim keeps its parameterization
  - 8b2090f shell: a round of revolution moved in past its radius has no offset
  - dfc7cb6 SMLib: analytic STEPInversion answers apart from its failures
  - cc33bc3 SMLib: a revolved face's offset is trimmed at the axis on either side of the frame's origin
  - 02e0dde SMLib, capi: the fuzzers' sanitizer findings fixed at their causes
  - 672a437 SMLib: the analytic surfaces and curves pass on their statuses
  - 0cd9d08 SMLib: the torus check tests its minor interval

## surface.surface_approximate

- 42 of 437 cases; first mismatching field: `values.out.len` 26, `values.out` 16.
- Fork commits:
  - dda881f Surface: a cover steps across the seam a surface wraps in STEP u, where a swapped surface's NURBS closes in v; the drafted peg's cone walls are one face
  - f8e6709 Fillet: every Ball's face is the rational patch on the ball's sphere, carrying that sphere as its definition; the SmSphere path is gone
  - 014f66f SMLib: SmProceduralBlend, a blend surface stored as its definition, its net the approximation
  - b518bbe merge: a cover crosses a wrapped seam only on the faces' trims, never retaken after the stored one
  - 701da4b fillet: a Coons corner patch's normal along its sides is made unit by its own length
  - eaa998a surface: a dropped point's foot from a guess is one the local solve converged to
  - 737c392 surface: a dropped point's foot is off along the surface by no more than the face's tolerance
  - 16bf1cb merge: a full-turn face takes in a coincident face across its seam, the seam run on through it
  - 67b8a1e boolean2d: a plane grows to cover a face however far away it lies
  - eae96cb SmSurface::ApproximateSurface puts a B-spline starter in Bezier form first
  - cae8618 fillet: a Coons corner patch is tangent to its sides' surfaces
  - 1f2b0c4 SMLib: a hole touching its outline at a corner is one loop of the face
  - and 6 more

## surface.surface_assert_valid

- 65 of 2902 cases; first mismatching field: `values.out.len` 58, `values.out` 7.
- Fork commits:
  - dda881f Surface: a cover steps across the seam a surface wraps in STEP u, where a swapped surface's NURBS closes in v; the drafted peg's cone walls are one face
  - 0c00de7 Surfaces: an analytic surface turned round by SwapUV keeps its definition as made, not re-derived from its net
  - f8e6709 Fillet: every Ball's face is the rational patch on the ball's sphere, carrying that sphere as its definition; the SmSphere path is gone
  - 6f811b8 Brep: a face shrunk to its trims keeps them where its surface's trim keeps its parameterization
  - af7a4d9 Torus: a horn torus rebuilt from its NURB stays the horn torus it was
  - 014f66f SMLib: SmProceduralBlend, a blend surface stored as its definition, its net the approximation
  - b518bbe merge: a cover crosses a wrapped seam only on the faces' trims, never retaken after the stored one
  - eaa998a surface: a dropped point's foot from a guess is one the local solve converged to
  - 737c392 surface: a dropped point's foot is off along the surface by no more than the face's tolerance
  - 16bf1cb merge: a full-turn face takes in a coincident face across its seam, the seam run on through it
  - 8b2090f shell: a round of revolution moved in past its radius has no offset
  - c2f1788 SSI: a plane sunk a hair into a bare cylinder crosses it in two lines
  - and 37 more

## surface.surface_cache_p1

- 16 of 872 cases; first mismatching field: `values.out` 16.
- Fork commits:
  - 729f893 SMLib: the surface cache is built on a blend's net, its node boxes grown by the net's bound
  - e149e05 SMLib: point solves, drops and the curve and surface caches pass on statuses
  - 02e0dde SMLib, capi: the fuzzers' sanitizer findings fixed at their causes

## surface.surface_p1

- 4 of 1000 cases; first mismatching field: `values.out.len` 2, `values.out` 2.
- Fork commits:
  - dda881f Surface: a cover steps across the seam a surface wraps in STEP u, where a swapped surface's NURBS closes in v; the drafted peg's cone walls are one face
  - f8e6709 Fillet: every Ball's face is the rational patch on the ball's sphere, carrying that sphere as its definition; the SmSphere path is gone
  - 6f811b8 Brep: a face shrunk to its trims keeps them where its surface's trim keeps its parameterization
  - 014f66f SMLib: SmProceduralBlend, a blend surface stored as its definition, its net the approximation
  - b518bbe merge: a cover crosses a wrapped seam only on the faces' trims, never retaken after the stored one
  - eaa998a surface: a dropped point's foot from a guess is one the local solve converged to
  - 737c392 surface: a dropped point's foot is off along the surface by no more than the face's tolerance
  - 16bf1cb merge: a full-turn face takes in a coincident face across its seam, the seam run on through it
  - 67b8a1e boolean2d: a plane grows to cover a face however far away it lies
  - eae96cb SmSurface::ApproximateSurface puts a B-spline starter in Bezier form first
  - e36feb1 SMLib UV trims: a drop that comes back in pieces is joined where they meet, never stretched onto the edge (solidcore#1383)
  - 758829c SMLib: an offset that collapses onto its spine is told apart from a missing one
  - and 8 more

## surface.surface_p10

- 26 of 867 cases; first mismatching field: `values.rec.len` 26.
- Fork commits:
  - dda881f Surface: a cover steps across the seam a surface wraps in STEP u, where a swapped surface's NURBS closes in v; the drafted peg's cone walls are one face
  - f8e6709 Fillet: every Ball's face is the rational patch on the ball's sphere, carrying that sphere as its definition; the SmSphere path is gone
  - 014f66f SMLib: SmProceduralBlend, a blend surface stored as its definition, its net the approximation
  - b518bbe merge: a cover crosses a wrapped seam only on the faces' trims, never retaken after the stored one
  - eaa998a surface: a dropped point's foot from a guess is one the local solve converged to
  - 737c392 surface: a dropped point's foot is off along the surface by no more than the face's tolerance
  - 16bf1cb merge: a full-turn face takes in a coincident face across its seam, the seam run on through it
  - 67b8a1e boolean2d: a plane grows to cover a face however far away it lies
  - eae96cb SmSurface::ApproximateSurface puts a B-spline starter in Bezier form first
  - 1f2b0c4 SMLib: a hole touching its outline at a corner is one loop of the face
  - 158cbe7 SMLib: a surface tree node with no surface cache has no original surface
  - 1612a7d SMLib: the curve and surface SE() calls pass on or discard with a reason
  - and 3 more

## surface.surface_p11

- 31 of 814 cases; first mismatching field: `values.deviation` 18, `error` 12, `values.data.len` 1.
- Fork commits:
  - dda881f Surface: a cover steps across the seam a surface wraps in STEP u, where a swapped surface's NURBS closes in v; the drafted peg's cone walls are one face
  - f8e6709 Fillet: every Ball's face is the rational patch on the ball's sphere, carrying that sphere as its definition; the SmSphere path is gone
  - 014f66f SMLib: SmProceduralBlend, a blend surface stored as its definition, its net the approximation
  - b518bbe merge: a cover crosses a wrapped seam only on the faces' trims, never retaken after the stored one
  - eaa998a surface: a dropped point's foot from a guess is one the local solve converged to
  - 737c392 surface: a dropped point's foot is off along the surface by no more than the face's tolerance
  - 16bf1cb merge: a full-turn face takes in a coincident face across its seam, the seam run on through it
  - 67b8a1e boolean2d: a plane grows to cover a face however far away it lies
  - eae96cb SmSurface::ApproximateSurface puts a B-spline starter in Bezier form first
  - 1f2b0c4 SMLib: a hole touching its outline at a corner is one loop of the face
  - 158cbe7 SMLib: a surface tree node with no surface cache has no original surface
  - 1612a7d SMLib: the curve and surface SE() calls pass on or discard with a reason
  - and 3 more

## surface.surface_p2

- 4139 of 17551 cases; first mismatching field: `values.status` 3006, `values.out.len` 410, `values.build` 371.
- Fork commits:
  - dda881f Surface: a cover steps across the seam a surface wraps in STEP u, where a swapped surface's NURBS closes in v; the drafted peg's cone walls are one face
  - f8e6709 Fillet: every Ball's face is the rational patch on the ball's sphere, carrying that sphere as its definition; the SmSphere path is gone
  - 6f811b8 Brep: a face shrunk to its trims keeps them where its surface's trim keeps its parameterization
  - 014f66f SMLib: SmProceduralBlend, a blend surface stored as its definition, its net the approximation
  - b518bbe merge: a cover crosses a wrapped seam only on the faces' trims, never retaken after the stored one
  - eaa998a surface: a dropped point's foot from a guess is one the local solve converged to
  - 737c392 surface: a dropped point's foot is off along the surface by no more than the face's tolerance
  - 16bf1cb merge: a full-turn face takes in a coincident face across its seam, the seam run on through it
  - 67b8a1e boolean2d: a plane grows to cover a face however far away it lies
  - eae96cb SmSurface::ApproximateSurface puts a B-spline starter in Bezier form first
  - e36feb1 SMLib UV trims: a drop that comes back in pieces is joined where they meet, never stretched onto the edge (solidcore#1383)
  - 758829c SMLib: an offset that collapses onto its spine is told apart from a missing one
  - and 8 more

## surface.surface_p3

- 94 of 2507 cases; first mismatching field: `values.out.len` 63, `values.out` 31.
- Fork commits:
  - dda881f Surface: a cover steps across the seam a surface wraps in STEP u, where a swapped surface's NURBS closes in v; the drafted peg's cone walls are one face
  - f8e6709 Fillet: every Ball's face is the rational patch on the ball's sphere, carrying that sphere as its definition; the SmSphere path is gone
  - 6f811b8 Brep: a face shrunk to its trims keeps them where its surface's trim keeps its parameterization
  - 014f66f SMLib: SmProceduralBlend, a blend surface stored as its definition, its net the approximation
  - b518bbe merge: a cover crosses a wrapped seam only on the faces' trims, never retaken after the stored one
  - eaa998a surface: a dropped point's foot from a guess is one the local solve converged to
  - 737c392 surface: a dropped point's foot is off along the surface by no more than the face's tolerance
  - 16bf1cb merge: a full-turn face takes in a coincident face across its seam, the seam run on through it
  - 67b8a1e boolean2d: a plane grows to cover a face however far away it lies
  - eae96cb SmSurface::ApproximateSurface puts a B-spline starter in Bezier form first
  - e36feb1 SMLib UV trims: a drop that comes back in pieces is joined where they meet, never stretched onto the edge (solidcore#1383)
  - 758829c SMLib: an offset that collapses onto its spine is told apart from a missing one
  - and 8 more

## surface.surface_p4

- 90 of 3377 cases; first mismatching field: `values.out` 90.
- Fork commits:
  - dda881f Surface: a cover steps across the seam a surface wraps in STEP u, where a swapped surface's NURBS closes in v; the drafted peg's cone walls are one face
  - f8e6709 Fillet: every Ball's face is the rational patch on the ball's sphere, carrying that sphere as its definition; the SmSphere path is gone
  - 6f811b8 Brep: a face shrunk to its trims keeps them where its surface's trim keeps its parameterization
  - 014f66f SMLib: SmProceduralBlend, a blend surface stored as its definition, its net the approximation
  - b518bbe merge: a cover crosses a wrapped seam only on the faces' trims, never retaken after the stored one
  - eaa998a surface: a dropped point's foot from a guess is one the local solve converged to
  - 737c392 surface: a dropped point's foot is off along the surface by no more than the face's tolerance
  - 16bf1cb merge: a full-turn face takes in a coincident face across its seam, the seam run on through it
  - 67b8a1e boolean2d: a plane grows to cover a face however far away it lies
  - eae96cb SmSurface::ApproximateSurface puts a B-spline starter in Bezier form first
  - e36feb1 SMLib UV trims: a drop that comes back in pieces is joined where they meet, never stretched onto the edge (solidcore#1383)
  - 758829c SMLib: an offset that collapses onto its spine is told apart from a missing one
  - and 8 more

## surface.surface_p5

- 120 of 3118 cases; first mismatching field: `values.out.len` 89, `values.out` 31.
- Fork commits:
  - dda881f Surface: a cover steps across the seam a surface wraps in STEP u, where a swapped surface's NURBS closes in v; the drafted peg's cone walls are one face
  - f8e6709 Fillet: every Ball's face is the rational patch on the ball's sphere, carrying that sphere as its definition; the SmSphere path is gone
  - 6f811b8 Brep: a face shrunk to its trims keeps them where its surface's trim keeps its parameterization
  - 014f66f SMLib: SmProceduralBlend, a blend surface stored as its definition, its net the approximation
  - b518bbe merge: a cover crosses a wrapped seam only on the faces' trims, never retaken after the stored one
  - eaa998a surface: a dropped point's foot from a guess is one the local solve converged to
  - 737c392 surface: a dropped point's foot is off along the surface by no more than the face's tolerance
  - 16bf1cb merge: a full-turn face takes in a coincident face across its seam, the seam run on through it
  - 67b8a1e boolean2d: a plane grows to cover a face however far away it lies
  - eae96cb SmSurface::ApproximateSurface puts a B-spline starter in Bezier form first
  - e36feb1 SMLib UV trims: a drop that comes back in pieces is joined where they meet, never stretched onto the edge (solidcore#1383)
  - 758829c SMLib: an offset that collapses onto its spine is told apart from a missing one
  - and 8 more

## surface.surface_p6

- 46 of 2104 cases; first mismatching field: `values.out` 24, `values.out.len` 22.
- Fork commits:
  - dda881f Surface: a cover steps across the seam a surface wraps in STEP u, where a swapped surface's NURBS closes in v; the drafted peg's cone walls are one face
  - f8e6709 Fillet: every Ball's face is the rational patch on the ball's sphere, carrying that sphere as its definition; the SmSphere path is gone
  - 014f66f SMLib: SmProceduralBlend, a blend surface stored as its definition, its net the approximation
  - b518bbe merge: a cover crosses a wrapped seam only on the faces' trims, never retaken after the stored one
  - eaa998a surface: a dropped point's foot from a guess is one the local solve converged to
  - 737c392 surface: a dropped point's foot is off along the surface by no more than the face's tolerance
  - 16bf1cb merge: a full-turn face takes in a coincident face across its seam, the seam run on through it
  - 67b8a1e boolean2d: a plane grows to cover a face however far away it lies
  - eae96cb SmSurface::ApproximateSurface puts a B-spline starter in Bezier form first
  - 1f2b0c4 SMLib: a hole touching its outline at a corner is one loop of the face
  - 158cbe7 SMLib: a surface tree node with no surface cache has no original surface
  - 1612a7d SMLib: the curve and surface SE() calls pass on or discard with a reason
  - and 3 more

## surface.surface_p7

- 42 of 3741 cases; first mismatching field: `values.out` 30, `values.out.len` 12.
- Fork commits:
  - dda881f Surface: a cover steps across the seam a surface wraps in STEP u, where a swapped surface's NURBS closes in v; the drafted peg's cone walls are one face
  - f8e6709 Fillet: every Ball's face is the rational patch on the ball's sphere, carrying that sphere as its definition; the SmSphere path is gone
  - 014f66f SMLib: SmProceduralBlend, a blend surface stored as its definition, its net the approximation
  - b518bbe merge: a cover crosses a wrapped seam only on the faces' trims, never retaken after the stored one
  - eaa998a surface: a dropped point's foot from a guess is one the local solve converged to
  - 737c392 surface: a dropped point's foot is off along the surface by no more than the face's tolerance
  - 16bf1cb merge: a full-turn face takes in a coincident face across its seam, the seam run on through it
  - 67b8a1e boolean2d: a plane grows to cover a face however far away it lies
  - eae96cb SmSurface::ApproximateSurface puts a B-spline starter in Bezier form first
  - 1f2b0c4 SMLib: a hole touching its outline at a corner is one loop of the face
  - 158cbe7 SMLib: a surface tree node with no surface cache has no original surface
  - 1612a7d SMLib: the curve and surface SE() calls pass on or discard with a reason
  - and 3 more

## surface.surface_p8

- 34 of 2108 cases; first mismatching field: `values.out.len` 28, `values.out` 6.
- Fork commits:
  - dda881f Surface: a cover steps across the seam a surface wraps in STEP u, where a swapped surface's NURBS closes in v; the drafted peg's cone walls are one face
  - f8e6709 Fillet: every Ball's face is the rational patch on the ball's sphere, carrying that sphere as its definition; the SmSphere path is gone
  - 014f66f SMLib: SmProceduralBlend, a blend surface stored as its definition, its net the approximation
  - b518bbe merge: a cover crosses a wrapped seam only on the faces' trims, never retaken after the stored one
  - eaa998a surface: a dropped point's foot from a guess is one the local solve converged to
  - 737c392 surface: a dropped point's foot is off along the surface by no more than the face's tolerance
  - 16bf1cb merge: a full-turn face takes in a coincident face across its seam, the seam run on through it
  - 67b8a1e boolean2d: a plane grows to cover a face however far away it lies
  - eae96cb SmSurface::ApproximateSurface puts a B-spline starter in Bezier form first
  - 1f2b0c4 SMLib: a hole touching its outline at a corner is one loop of the face
  - 158cbe7 SMLib: a surface tree node with no surface cache has no original surface
  - 1612a7d SMLib: the curve and surface SE() calls pass on or discard with a reason
  - and 3 more

## surface.surface_point_solve_p2

- 35 of 1366 cases; first mismatching field: `values.out` 34, `values.out.len` 1.
- Fork commits:
  - e6f0a9e SMLib: a point solve's suspicious control-net guess stays a candidate against the nearest control point
  - e149e05 SMLib: point solves, drops and the curve and surface caches pass on statuses

## surface.surface_point_solve_p3

- 2 of 1365 cases; first mismatching field: `values.out` 2.
- Fork commits:
  - e6f0a9e SMLib: a point solve's suspicious control-net guess stays a candidate against the nearest control point
  - e149e05 SMLib: point solves, drops and the curve and surface caches pass on statuses

## surface.torus

- 139 of 758 cases; first mismatching field: `values.out.len` 78, `values.out` 60, `values.build` 1.
- Fork commits:
  - af7a4d9 Torus: a horn torus rebuilt from its NURB stays the horn torus it was
  - 6787ae7 torus: a minor circle starting on the axis takes its direction from its middle
  - 40344bf SMLib: a torus found in a NURBS that runs once round its tube has a whole turn of it
  - 6136a04 SMLib: a torus found in a NURBS has its major radius from the torus's origin
  - 7ea6a96 SMLib: a torus or sphere found in a NURBS keeps that NURBS; a half meridian measured to its rounding is a half meridian
  - d9f20b0 SMLib: a torus rebuilt from its NURB keeps a generating interval below its equator
  - 758829c SMLib: an offset that collapses onto its spine is told apart from a missing one
  - 9853e54 SMLib: a revolved curve is its analytic surface, per profile curve
  - 672a437 SMLib: the analytic surfaces and curves pass on their statuses
  - daba678 SMLib: the torus check reports a spindle torus running through its apex
  - 0cd9d08 SMLib: the torus check tests its minor interval

## surface.volume_p2

- 795 of 11220 cases; first mismatching field: `values.out.len` 689, `values.out` 101, `ints.status` 5.
- Fork commits (none on its own files; its dependencies'):
  - 922f333 periodic interval: a value just past the period's top is a period down, not the interval's start
  - d5c5cfd SMLib: an N-D solve whose guess can't be evaluated finds no solution
  - d2757e1 SMLib: the local solvers take no value from a failed evaluation
  - 19310fa SMLib: composite curves, gaps, pseudo boxes, point sets and the geometry utilities pass on statuses
  - ca7d1f5 capi: no shim drops a status SMLib returns
  - f8e6709 Fillet: every Ball's face is the rational patch on the ball's sphere, carrying that sphere as its definition; the SmSphere path is gone
  - 014f66f SMLib: SmProceduralBlend, a blend surface stored as its definition, its net the approximation

Measured in #2512 (WSL, the fork's libraries built with `harness/cpp/CMakeLists.txt`'s flags).
