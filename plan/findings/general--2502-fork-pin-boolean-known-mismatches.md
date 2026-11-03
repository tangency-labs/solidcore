# boolean: cases that move with the difftest pin to the SMLib fork (2026-11-03)

Units: boolean.api_brep, boolean.merge_p3, boolean.merge_p4, boolean.topology_intersector_p1, boolean.topology_intersector_p2.

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

## boolean.api_brep

- 91 of 1134 cases; first mismatching field: `values.attributes` 26, `values.areas` 20, `values.measures` 18.
- Fork commits:
  - 54148db SMLib: draft corners that keep the curves' curvature
  - 4803d6d SMLib: offsets and shells pass on the statuses they were dropping
  - 5d3939e SMLib: a curve sweep says where it placed the profile

## boolean.merge_p3

- 11 of 438 cases; first mismatching field: `error` 11.
- Fork commits:
  - 3ed0dee Boolean: every edge, vertex and face it destroys leaves its lists, sets and mate maps as it goes
  - cac11ea fillets: the clean-up removes a slit, an input edge with one face on both sides
  - 16bf1cb merge: a full-turn face takes in a coincident face across its seam, the seam run on through it
  - 1ccd461 fillet: the clean-up keeps the target's edges and faces by a mark, not by address
  - 58232c7 ManifoldBoolean puts each region's enclosing shell first
  - 9d8e2c5 SMLib: a boolean's clean-up removes the edges its merging leaves inside one face
  - 242ca97 SMLib: a slice keeps a face coincident with its sheet only on the side kept
  - 7f23af3 SMLib: a boolean of solids touching along a curve inside a face leaves no slit loop
  - c340962 SMLib: faces that bound nothing within tolerance are removed after a boolean
  - d98540c SMLib: a boolean keeps the second operand's UV trims
  - e7ebc8d SMLib: booleans' merge, stitch and topology intersection pass on statuses
  - b7d7b96 SMLib: the clean-up after a fillet leaves the target's own faces apart
  - and 1 more

## boolean.merge_p4

- 24 of 311 cases; first mismatching field: `values.primary.dbls.len` 10, `ints.status_and_state` 9, `values.component.0.dbls.len` 5.
- Fork commits:
  - 3ed0dee Boolean: every edge, vertex and face it destroys leaves its lists, sets and mate maps as it goes
  - cac11ea fillets: the clean-up removes a slit, an input edge with one face on both sides
  - 16bf1cb merge: a full-turn face takes in a coincident face across its seam, the seam run on through it
  - 1ccd461 fillet: the clean-up keeps the target's edges and faces by a mark, not by address
  - 58232c7 ManifoldBoolean puts each region's enclosing shell first
  - 9d8e2c5 SMLib: a boolean's clean-up removes the edges its merging leaves inside one face
  - 242ca97 SMLib: a slice keeps a face coincident with its sheet only on the side kept
  - 7f23af3 SMLib: a boolean of solids touching along a curve inside a face leaves no slit loop
  - c340962 SMLib: faces that bound nothing within tolerance are removed after a boolean
  - d98540c SMLib: a boolean keeps the second operand's UV trims
  - e7ebc8d SMLib: booleans' merge, stitch and topology intersection pass on statuses
  - b7d7b96 SMLib: the clean-up after a fillet leaves the target's own faces apart
  - and 1 more

## boolean.topology_intersector_p1

- 6 of 845 cases; first mismatching field: `values.state.len` 5, `values.state` 1.
- Fork commits:
  - 3ed0dee Boolean: every edge, vertex and face it destroys leaves its lists, sets and mate maps as it goes
  - 2e9e6b2 capi tests: the Menger row join counts the work it does instead of timing it
  - 84a7998 booleans: each face is paired only with the faces whose box comes near its own
  - 60b1c6a boolean2d: faces apart are not intersected, and loops apart are not classified
  - e881cdb boolean: closing an intersection loop relates two edges only when one lies on the other
  - fcd5ea5 SMLib: no preprocessor directive inside a macro's argument
  - 03e7bf1 SMLib: an intersection curve inside one vertex is that vertex, not a closed edge on it
  - e7ebc8d SMLib: booleans' merge, stitch and topology intersection pass on statuses
  - a2d0c5e SMLib boolean: a gap fixed on one side only leaves no freed edge in the maps
  - 55cc869 SMLib: a failed imprint stops the boolean

## boolean.topology_intersector_p2

- 13 of 913 cases; first mismatching field: `values.state` 13.
- Fork commits:
  - 3ed0dee Boolean: every edge, vertex and face it destroys leaves its lists, sets and mate maps as it goes
  - 2e9e6b2 capi tests: the Menger row join counts the work it does instead of timing it
  - 84a7998 booleans: each face is paired only with the faces whose box comes near its own
  - 60b1c6a boolean2d: faces apart are not intersected, and loops apart are not classified
  - e881cdb boolean: closing an intersection loop relates two edges only when one lies on the other
  - fcd5ea5 SMLib: no preprocessor directive inside a macro's argument
  - 03e7bf1 SMLib: an intersection curve inside one vertex is that vertex, not a closed edge on it
  - e7ebc8d SMLib: booleans' merge, stitch and topology intersection pass on statuses
  - a2d0c5e SMLib boolean: a gap fixed on one side only leaves no freed edge in the maps
  - 55cc869 SMLib: a failed imprint stops the boolean

Measured in #2512 (WSL, the fork's libraries built with `harness/cpp/CMakeLists.txt`'s flags).
