# blend: cases that move with the difftest pin to the SMLib fork (2026-11-03)

Units: blend.fillet_custom_solver, blend.fillet_executive_p1, blend.fillet_executive_p2, blend.fillet_geom_p1, blend.fillet_geom_p3, blend.fillet_geom_p4, blend.fillet_geom_p5, blend.fillet_geom_p6, blend.fillet_intersector, blend.fillet_removal, blend.fillet_solver_p4, blend.fillet_standard_solver_p2, blend.hcr.

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

## blend.fillet_custom_solver

- 4 of 451 cases; first mismatching field: `values.guess` 2, `values.rail2_cpts.len` 1, `values.jac` 1.
- Fork commits:
  - 66b02de SMLib: the fillet solvers, corners and blends pass on their statuses

## blend.fillet_executive_p1

- 5 of 155 cases; first mismatching field: `values.dbls` 3, `values.dbls.len` 1, `ints.ints` 1.
- Fork commits:
  - 9822462 Fillet: a round that goes second meets the rounds that went first as one blend, not with a smooth-vertex End
  - ef57570 Fillet: an open setback corner is the Setback's hole: its G1 face curves the spring parabolas, no fill
  - 41a39ec Fillet: a round ends at a smooth vertex on the plane square to its edge, and a rail past its side edge is refused naming the radius that reaches it
  - d657afa Fillet: a trim that leaks past a fillet's rails is an error, and an edge of more than two faces bounds it
  - 9846aa0 SMLib compiles clean under the USD build's warning flags (-Wall -Wshadow=local -Wunused-parameter, -Werror): 19 warnings fixed at their causes
  - fe980d7 Fillet: the big-radius conversion of a two-blend corner goes; no corner reaches it
  - b8ec238 SMLib: an equal area thumbweight is solved once per order and angle within a blend operation
  - f1d81f4 Fillet: a G1 setback fill is taken only certified convex; the rule's setbacks grow to the least scale that certifies, within each edge's room, and one asked for that is too short is refused naming the setback that builds
  - 61547af Fillet: a corner of three or more rounds of unequal radii, of curvature sections, or set back on request is the blend corner's Setback, filled by the setback split
  - 5e8c9c1 Fillet: a refusal's bound is bracketed strictly inside its search, not at its end
  - 49f393d Fillet: size the remaining refusals from their causes: rails past a face, two-distance chamfers, chord widths, polygon faces
  - a3d6d81 capi: smlib_blend_fit_size, the largest size that builds, from each refusal's own cause
  - and 57 more

## blend.fillet_executive_p2

- 2 of 278 cases; first mismatching field: `values.out` 1, `error` 1.
- Fork commits:
  - 9822462 Fillet: a round that goes second meets the rounds that went first as one blend, not with a smooth-vertex End
  - ef57570 Fillet: an open setback corner is the Setback's hole: its G1 face curves the spring parabolas, no fill
  - 41a39ec Fillet: a round ends at a smooth vertex on the plane square to its edge, and a rail past its side edge is refused naming the radius that reaches it
  - d657afa Fillet: a trim that leaks past a fillet's rails is an error, and an edge of more than two faces bounds it
  - 9846aa0 SMLib compiles clean under the USD build's warning flags (-Wall -Wshadow=local -Wunused-parameter, -Werror): 19 warnings fixed at their causes
  - fe980d7 Fillet: the big-radius conversion of a two-blend corner goes; no corner reaches it
  - b8ec238 SMLib: an equal area thumbweight is solved once per order and angle within a blend operation
  - f1d81f4 Fillet: a G1 setback fill is taken only certified convex; the rule's setbacks grow to the least scale that certifies, within each edge's room, and one asked for that is too short is refused naming the setback that builds
  - 61547af Fillet: a corner of three or more rounds of unequal radii, of curvature sections, or set back on request is the blend corner's Setback, filled by the setback split
  - 5e8c9c1 Fillet: a refusal's bound is bracketed strictly inside its search, not at its end
  - 49f393d Fillet: size the remaining refusals from their causes: rails past a face, two-distance chamfers, chord widths, polygon faces
  - a3d6d81 capi: smlib_blend_fit_size, the largest size that builds, from each refusal's own cause
  - and 53 more

## blend.fillet_geom_p1

- 29 of 2199 cases; first mismatching field: `values.e00_ivl` 10, `values.status` 7, `values.e01_ivl` 2.
- Fork commits:
  - 63c7962 Fillet: an End's rail ends on its side edge continued across a seam's vertex, and a rail passing over the rim is refused on the whole of it
  - 41a39ec Fillet: a round ends at a smooth vertex on the plane square to its edge, and a rail past its side edge is refused naming the radius that reaches it
  - d657afa Fillet: a trim that leaks past a fillet's rails is an error, and an edge of more than two faces bounds it
  - 2ac7a17 Fillet: a cross section no longer looks for the closed-edge corner, which has none
  - 26b18e4 Fillet: SmFilletEdge::AdjustCrossSection goes; nothing calls it
  - 61547af Fillet: a corner of three or more rounds of unequal radii, of curvature sections, or set back on request is the blend corner's Setback, filled by the setback split
  - 53c47b8 Fillet: a blend ending at a vertex is the blend corner's End; the Nx1, open and degenerate corners are gone
  - 543d25a Fillet: every G1 continuation of blends is the Join's shared section; the 1x1 and 2x2 corners are gone
  - 3e26f76 Fillet: a mixed-convexity End's rail ends where the rail meets the end face's surface
  - 794946e Fillet: a round whose rails both run on its neighbours' rails keeps its face
  - 32fd80a Fillet: the mitre of two equal cylinders is their exact plane section, and its faces keep their exact trim curves
  - 74e9beb Fillet: the Mitre turn and the Join are the blend corner's; SmFilletNx2Corner is gone
  - and 27 more

## blend.fillet_geom_p3

- 45 of 325 cases; first mismatching field: `values.flat` 44, `values.flat.len` 1.
- Fork commits:
  - 63c7962 Fillet: an End's rail ends on its side edge continued across a seam's vertex, and a rail passing over the rim is refused on the whole of it
  - 41a39ec Fillet: a round ends at a smooth vertex on the plane square to its edge, and a rail past its side edge is refused naming the radius that reaches it
  - d657afa Fillet: a trim that leaks past a fillet's rails is an error, and an edge of more than two faces bounds it
  - 2ac7a17 Fillet: a cross section no longer looks for the closed-edge corner, which has none
  - 26b18e4 Fillet: SmFilletEdge::AdjustCrossSection goes; nothing calls it
  - 61547af Fillet: a corner of three or more rounds of unequal radii, of curvature sections, or set back on request is the blend corner's Setback, filled by the setback split
  - 53c47b8 Fillet: a blend ending at a vertex is the blend corner's End; the Nx1, open and degenerate corners are gone
  - 543d25a Fillet: every G1 continuation of blends is the Join's shared section; the 1x1 and 2x2 corners are gone
  - 3e26f76 Fillet: a mixed-convexity End's rail ends where the rail meets the end face's surface
  - 794946e Fillet: a round whose rails both run on its neighbours' rails keeps its face
  - 32fd80a Fillet: the mitre of two equal cylinders is their exact plane section, and its faces keep their exact trim curves
  - 74e9beb Fillet: the Mitre turn and the Join are the blend corner's; SmFilletNx2Corner is gone
  - and 20 more

## blend.fillet_geom_p4

- 13 of 400 cases; first mismatching field: `values.status` 12, `values.uv1_cpts` 1.
- Fork commits:
  - 63c7962 Fillet: an End's rail ends on its side edge continued across a seam's vertex, and a rail passing over the rim is refused on the whole of it
  - 41a39ec Fillet: a round ends at a smooth vertex on the plane square to its edge, and a rail past its side edge is refused naming the radius that reaches it
  - d657afa Fillet: a trim that leaks past a fillet's rails is an error, and an edge of more than two faces bounds it
  - 2ac7a17 Fillet: a cross section no longer looks for the closed-edge corner, which has none
  - 26b18e4 Fillet: SmFilletEdge::AdjustCrossSection goes; nothing calls it
  - 61547af Fillet: a corner of three or more rounds of unequal radii, of curvature sections, or set back on request is the blend corner's Setback, filled by the setback split
  - 53c47b8 Fillet: a blend ending at a vertex is the blend corner's End; the Nx1, open and degenerate corners are gone
  - 543d25a Fillet: every G1 continuation of blends is the Join's shared section; the 1x1 and 2x2 corners are gone
  - 3e26f76 Fillet: a mixed-convexity End's rail ends where the rail meets the end face's surface
  - 794946e Fillet: a round whose rails both run on its neighbours' rails keeps its face
  - 32fd80a Fillet: the mitre of two equal cylinders is their exact plane section, and its faces keep their exact trim curves
  - 74e9beb Fillet: the Mitre turn and the Join are the blend corner's; SmFilletNx2Corner is gone
  - and 20 more

## blend.fillet_geom_p5

- 102 of 300 cases; first mismatching field: `values.r0_c3d_cpts` 51, `values.r0_c3d_cpts.len` 26, `values.r0_uv1_cpts` 19.
- Fork commits:
  - 63c7962 Fillet: an End's rail ends on its side edge continued across a seam's vertex, and a rail passing over the rim is refused on the whole of it
  - 41a39ec Fillet: a round ends at a smooth vertex on the plane square to its edge, and a rail past its side edge is refused naming the radius that reaches it
  - d657afa Fillet: a trim that leaks past a fillet's rails is an error, and an edge of more than two faces bounds it
  - 2ac7a17 Fillet: a cross section no longer looks for the closed-edge corner, which has none
  - 26b18e4 Fillet: SmFilletEdge::AdjustCrossSection goes; nothing calls it
  - 61547af Fillet: a corner of three or more rounds of unequal radii, of curvature sections, or set back on request is the blend corner's Setback, filled by the setback split
  - 53c47b8 Fillet: a blend ending at a vertex is the blend corner's End; the Nx1, open and degenerate corners are gone
  - 543d25a Fillet: every G1 continuation of blends is the Join's shared section; the 1x1 and 2x2 corners are gone
  - 3e26f76 Fillet: a mixed-convexity End's rail ends where the rail meets the end face's surface
  - 794946e Fillet: a round whose rails both run on its neighbours' rails keeps its face
  - 32fd80a Fillet: the mitre of two equal cylinders is their exact plane section, and its faces keep their exact trim curves
  - 74e9beb Fillet: the Mitre turn and the Join are the blend corner's; SmFilletNx2Corner is gone
  - and 20 more

## blend.fillet_geom_p6

- 14 of 400 cases; first mismatching field: `values.flat.len` 9, `values.flat` 5.
- Fork commits:
  - 63c7962 Fillet: an End's rail ends on its side edge continued across a seam's vertex, and a rail passing over the rim is refused on the whole of it
  - 41a39ec Fillet: a round ends at a smooth vertex on the plane square to its edge, and a rail past its side edge is refused naming the radius that reaches it
  - d657afa Fillet: a trim that leaks past a fillet's rails is an error, and an edge of more than two faces bounds it
  - 2ac7a17 Fillet: a cross section no longer looks for the closed-edge corner, which has none
  - 26b18e4 Fillet: SmFilletEdge::AdjustCrossSection goes; nothing calls it
  - 61547af Fillet: a corner of three or more rounds of unequal radii, of curvature sections, or set back on request is the blend corner's Setback, filled by the setback split
  - 53c47b8 Fillet: a blend ending at a vertex is the blend corner's End; the Nx1, open and degenerate corners are gone
  - 543d25a Fillet: every G1 continuation of blends is the Join's shared section; the 1x1 and 2x2 corners are gone
  - 3e26f76 Fillet: a mixed-convexity End's rail ends where the rail meets the end face's surface
  - 794946e Fillet: a round whose rails both run on its neighbours' rails keeps its face
  - 32fd80a Fillet: the mitre of two equal cylinders is their exact plane section, and its faces keep their exact trim curves
  - 74e9beb Fillet: the Mitre turn and the Join are the blend corner's; SmFilletNx2Corner is gone
  - and 20 more

## blend.fillet_intersector

- 7 of 1292 cases; first mismatching field: `values.fillet_point_step` 5, `values.rec` 2.
- Fork commits:
  - 4438859 fillets: a variable round's rails follow the fillet between its traced points
  - 66b02de SMLib: the fillet solvers, corners and blends pass on their statuses
  - 6eecf04 SMLib: a variable round turns where its radius law has a corner

## blend.fillet_removal

- 27 of 170 cases; first mismatching field: `values.record.len` 19, `values.record` 6, `values.edge_curve_owners.len` 2.
- Fork commits:
  - 661f79c SMLib: gluing two vertices joins their shells in every region they share
  - 5b0534b SMLib: a four-sided fillet can be removed exactly only
  - 52bfef0 SMLib: a ring fillet is removed to the sharp edge where its faces meet

## blend.fillet_solver_p4

- 154 of 651 cases; first mismatching field: `values.fillet_point_step` 124, `values.point` 14, `ints.status` 14.
- Fork commits:
  - 74c4ae6 Fillet: the through point at a non-blending edge's end is solved where the trace runs, over its faces' extended surfaces
  - e916d9e Fillet: at an end its corner does not blend, a trace that stops short of the edge's end runs on to it, its section there a through point, before it is extended
  - 5e8c9c1 Fillet: a refusal's bound is bracketed strictly inside its search, not at its end
  - 49f393d Fillet: size the remaining refusals from their causes: rails past a face, two-distance chamfers, chord widths, polygon faces
  - 53c47b8 Fillet: a blend ending at a vertex is the blend corner's End; the Nx1, open and degenerate corners are gone
  - 543d25a Fillet: every G1 continuation of blends is the Join's shared section; the 1x1 and 2x2 corners are gone
  - a3d6d81 capi: smlib_blend_fit_size, the largest size that builds, from each refusal's own cause
  - a34b5e8 Fillet: delete the transversality check, unreachable where the section solve finds a section
  - 8c40818 Fillet: a ball's faces are judged with the solver's radius about its centre, unequal offsets each with its own
  - 20ba5a9 Fillet: a constant ball's limit is its least along the edge, not the least over the sections
  - e568e45 Fillet: whether a blend's section meets the part is decided exactly, not by sampling it
  - aac4288 Fillet: a blend that takes the part to its tip is no cut, and unequal offsets are judged for folds
  - and 37 more

## blend.fillet_standard_solver_p2

- 111 of 525 cases; first mismatching field: `values.out` 110, `values.out.len` 1.
- Fork commits:
  - 410dfe1 Fillet: the rolling ball's exact local checks along every blend, refusing with the radius that works
  - 3aae1c4 fillet: the chord and the circle are the section generator's orders 0 and 1
  - 0bd3f4d fillet: the G^n section replaces the blend-curve cross-section, certified span by span
  - 5fdefb9 fillets: the rails of a round meeting at its end end at one point
  - b7ec55b fillet: a knife edge rounds to the ball touching both faces back from it
  - 5fbffc9 fillet: a gentle corner in a variable radius splits the round
  - a9277fd SMLib: a fillet's trace reaches its rails across the edge
  - 66b02de SMLib: the fillet solvers, corners and blends pass on their statuses
  - 6eecf04 SMLib: a variable round turns where its radius law has a corner
  - 2412a80 SMLib: a radius law that can't be evaluated fails the round instead of offsetting by an unset radius
  - fd3b6a7 SMLib: blend-curve fillet sections can sit at thirds of the way to the corner

## blend.hcr

- 10 of 481 cases; first mismatching field: `values.record` 5, `values.record.len` 4, `error` 1.
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
  - and 71 more

Measured in #2512 (WSL, the fork's libraries built with `harness/cpp/CMakeLists.txt`'s flags).
