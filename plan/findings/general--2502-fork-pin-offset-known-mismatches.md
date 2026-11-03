# offset: cases that move with the difftest pin to the SMLib fork (2026-11-03)

Units: offset.offset_executive_p1, offset.offset_executive_p2, offset.offset_geometry_creation.

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

## offset.offset_executive_p1

- 6 of 1060 cases; first mismatching field: `ints.status` 5, `ints.probe` 1.
- Fork commits:
  - 1b016ca shell: an offset face's loops winding the other way under an opening are the cavity's floor
  - 37af643 offset: an inset keeps the regions that are cavity, not every finite region
  - 9d8ace9 offset: an edge whose faces' offsets could not be mated is intersected by the merge
  - 1c2486f shell: the cone-apex offsets take names of their own, not the function's
  - d384e86 thicken: a flat side wall is on its plane, a cylinder's end rings among them
  - ccbf77a offset faces: a step wall the move closes to nothing is removed
  - 3626d4b shell: a cylinder extended round its axis keeps SMLib's circle knots and comes back a cone
  - 34318f4 shell: a wall that would turn a face inside out past a step is SMLIB_ERR_UNSUPPORTED
  - d5dadaa shell: the apex of a cone (a drilled hole's tip) is offset as a sphere
  - b283ea2 offset: a cone or cylinder is extended as itself
  - cd85d34 shell: edges between analytic faces are their exact lines and conics
  - 758829c SMLib: an offset that collapses onto its spine is told apart from a missing one
  - and 12 more

## offset.offset_executive_p2

- 405 of 589 cases; first mismatching field: `values.executive.geometry_parameters` 287, `ints.result.edge_kinds` 76, `values.result.geometry` 28.
- Fork commits:
  - 1b016ca shell: an offset face's loops winding the other way under an opening are the cavity's floor
  - 37af643 offset: an inset keeps the regions that are cavity, not every finite region
  - 9d8ace9 offset: an edge whose faces' offsets could not be mated is intersected by the merge
  - 1c2486f shell: the cone-apex offsets take names of their own, not the function's
  - d384e86 thicken: a flat side wall is on its plane, a cylinder's end rings among them
  - ccbf77a offset faces: a step wall the move closes to nothing is removed
  - 3626d4b shell: a cylinder extended round its axis keeps SMLib's circle knots and comes back a cone
  - 34318f4 shell: a wall that would turn a face inside out past a step is SMLIB_ERR_UNSUPPORTED
  - d5dadaa shell: the apex of a cone (a drilled hole's tip) is offset as a sphere
  - b283ea2 offset: a cone or cylinder is extended as itself
  - cd85d34 shell: edges between analytic faces are their exact lines and conics
  - 758829c SMLib: an offset that collapses onto its spine is told apart from a missing one
  - and 14 more

## offset.offset_geometry_creation

- 212 of 926 cases; first mismatching field: `values.state` 182, `values.state.len` 30.
- Fork commits:
  - b439d73 Offset: TODO where the shell's exact offset of a Ball patch loses its definition
  - 3aae1c4 fillet: the chord and the circle are the section generator's orders 0 and 1
  - 9d8ace9 offset: an edge whose faces' offsets could not be mated is intersected by the merge
  - 2e1bb8c shell: four faces meeting at a vertex offset to two vertices joined by a new edge
  - 90aabff offset faces: a step wall between lines stored with different weights is ruled as the plane they span
  - cd85d34 shell: edges between analytic faces are their exact lines and conics
  - 7ee0876 Offset: an open edge's second trim point on a closed lifted curve is its other vertex, not the first point's seam twin
  - 98ab1b0 SMLib: an offset's analytic copy grows only along its rulings, not round a turn
  - fddc366 SMLib: an offset edge between two analytic faces is their exact line or conic
  - 758829c SMLib: an offset that collapses onto its spine is told apart from a missing one
  - 5d428ae SMLib: a local operation that moves one of two tangent faces along its normal builds the step
  - 4803d6d SMLib: offsets and shells pass on the statuses they were dropping
  - and 3 more

**Known cause behind many of them:** the fork's offset tangency tolerance is 0.5 degrees, v1.0.1's
2 (`SmOffsetGeometryCreation::m_dTangencyTolRadians`, tangency-labs/smlib 2e1bb8c); the executive
record's `geometry_parameters[2]` shows it. Five runner tests of offset.offset_executive_p2 that
build their own cases and assert parity are `#[ignore]`d with that reason until the flag flips.

Measured in #2512 (WSL, the fork's libraries built with `harness/cpp/CMakeLists.txt`'s flags).
