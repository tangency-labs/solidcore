# nurbs: cases that move with the difftest pin to the SMLib fork (2026-11-03)

Units: nurbs.b_spline_curve_n_lib, nurbs.b_spline_curve_p3, nurbs.b_spline_surface_n_lib_p2, nurbs.b_spline_surface_p2.

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

## nurbs.b_spline_curve_n_lib

- 4 of 2436 cases; first mismatching field: `values.this_cpts.len` 4.
- Fork commits:
  - 0cc8cd4 curves: CreateInterpolatingCurve interpolates at given parameters, with its end derivatives
  - e0c9abc SMLib: trimming a UV curve keeps its short end spans

## nurbs.b_spline_curve_p3

- 15 of 2460 cases; first mismatching field: `values.flag` 9, `values.report.len` 4, `values.other_cpts` 2.
- Fork commits:
  - ecf1963 check: a curve's zero-speed end at a vertex where its edge's faces touch is no kink
  - 54ca3e0 face, check: a loop far smaller than its face is oriented, and its trim measured, on the surface's scale
  - b4f8c57 SMLib: a UV trim's radius of curvature is measured on its surface
  - c343033 SMLib: a UV trim's turn at its ends is judged on its surface; the check no longer kinks a circle's trim
  - 4f2f7cb SMLib: the raw B-spline surface constructor copies all of its V knots
  - 8088f5e SMLib: a disjoint interval is no work, and new objects are freed on failure
  - 8843967 SMLib: B-spline curves pass on their statuses
  - cc4289d SMLib: the B-spline curve check reads the minimum radius of curvature

## nurbs.b_spline_surface_n_lib_p2

- 13 of 774 cases; first mismatching field: `values.patch0_cpts.len` 13.
- Fork commits:
  - 701da4b fillet: a Coons corner patch's normal along its sides is made unit by its own length
  - cae8618 fillet: a Coons corner patch is tangent to its sides' surfaces
  - 5d3939e SMLib: a curve sweep says where it placed the profile

## nurbs.b_spline_surface_p2

- 21 of 815 cases; first mismatching field: `values.curve_cpts` 12, `values.status` 6, `values.curve_cpts.len` 3.
- Fork commits:
  - 0c00de7 Surfaces: an analytic surface turned round by SwapUV keeps its definition as made, not re-derived from its net
  - a766bb4 B-spline surface check: an internal pole is sought in the row through the knot
  - 4f2f7cb SMLib: the raw B-spline surface constructor copies all of its V knots
  - 8088f5e SMLib: a disjoint interval is no work, and new objects are freed on failure
  - cb7cedb SMLib: B-spline surfaces pass on their statuses
  - 7af1987 SMLib: an analytic copy of a split face's surface keeps its parameterization

Measured in #2512 (WSL, the fork's libraries built with `harness/cpp/CMakeLists.txt`'s flags).
