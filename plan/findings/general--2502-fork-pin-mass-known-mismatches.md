# mass: cases that move with the difftest pin to the SMLib fork (2026-11-03)

Units: mass.integrator.

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

## mass.integrator

- 86 of 1628 cases; first mismatching field: `values.out` 82, `values.out.len` 4.
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
  - and 38 more

Measured in #2512 (WSL, the fork's libraries built with `harness/cpp/CMakeLists.txt`'s flags).
