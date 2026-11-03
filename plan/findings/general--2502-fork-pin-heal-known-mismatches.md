# heal: cases that move with the difftest pin to the SMLib fork (2026-11-03)

Units: heal.face_props, heal.gap, heal.heal_data_p1, heal.heal_data_p2, heal.heal_data_p3, heal.stitch_p1, heal.stitch_p2.

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

## heal.face_props

- 96 of 1873 cases; first mismatching field: `values.v` 96.
- Fork commits:
  - 45e7b7f SMLib: face properties pass on their loops' statuses
  - adf2b49 SMLib: SmFace::HasSeamProblem returns a status; the seam types say the problem

## heal.gap

- 9 of 749 cases; first mismatching field: `values.state` 9.
- Fork commits:
  - 19310fa SMLib: composite curves, gaps, pseudo boxes, point sets and the geometry utilities pass on statuses
  - ca7d1f5 capi: no shim drops a status SMLib returns

## heal.heal_data_p1

- 17 of 495 cases; first mismatching field: `values.out` 17.
- Fork commits:
  - f98bdf5 build: macOS (Apple Silicon) builds and tests with the documented commands

## heal.heal_data_p2

- 38 of 741 cases; first mismatching field: `values.out` 38.
- Fork commits:
  - f98bdf5 build: macOS (Apple Silicon) builds and tests with the documented commands

## heal.heal_data_p3

- 8 of 702 cases; first mismatching field: `values.out` 5, `values.out.len` 3.
- Fork commits:
  - f98bdf5 build: macOS (Apple Silicon) builds and tests with the documented commands

## heal.stitch_p1

- 4 of 674 cases; first mismatching field: `values.out` 4.
- Fork commits:
  - 351c57c stitch: a stitched solid is a nested solid, its cavities void
  - 54148db SMLib: draft corners that keep the curves' curvature
  - 4803d6d SMLib: offsets and shells pass on the statuses they were dropping
  - e7ebc8d SMLib: booleans' merge, stitch and topology intersection pass on statuses
  - 5d3939e SMLib: a curve sweep says where it placed the profile

## heal.stitch_p2

- 13 of 841 cases; first mismatching field: `values.out.len` 9, `values.out` 4.
- Fork commits:
  - 351c57c stitch: a stitched solid is a nested solid, its cavities void
  - e7ebc8d SMLib: booleans' merge, stitch and topology intersection pass on statuses

Measured in #2512 (WSL, the fork's libraries built with `harness/cpp/CMakeLists.txt`'s flags).
