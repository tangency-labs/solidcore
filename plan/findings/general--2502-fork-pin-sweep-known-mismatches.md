# sweep: cases that move with the difftest pin to the SMLib fork (2026-11-03)

Units: sweep.api_brep, sweep.api_surfaces, sweep.rotational_sweep_geometry, sweep.topology_sweep, sweep.trans_sweep_geometry.

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

## sweep.api_brep

- 26 of 696 cases; first mismatching field: `values.counts` 19, `values.status` 4, `values.linear_edge_samples` 2.
- Fork commits:
  - 54148db SMLib: draft corners that keep the curves' curvature
  - 4803d6d SMLib: offsets and shells pass on the statuses they were dropping
  - 5d3939e SMLib: a curve sweep says where it placed the profile

## sweep.api_surfaces

- 11 of 885 cases; first mismatching field: `values.faces` 10, `values.revolution_points` 1.
- Fork commits:
  - c8f6a31 SMLib: sweeps, lofts and extrusions pass on their statuses

## sweep.rotational_sweep_geometry

- 104 of 823 cases; first mismatching field: `values.out` 81, `values.out.len` 23.
- Fork commits:
  - c4e378d revolve: a circle reaching the axis within the tolerance sweeps a horn torus
  - b6f88e4 revolve: an arc with its centre across the axis sweeps a B-spline surface of revolution
  - c4209ab SMLib: a revolved arc centred within tolerance of the axis is a sphere through its ends, its pole on the axis
  - 9a160ba SMLib: a revolved arc's sphere is built on its own domain, not trimmed from a whole one
  - 9853e54 SMLib: a revolved curve is its analytic surface, per profile curve
  - c8f6a31 SMLib: sweeps, lofts and extrusions pass on their statuses

## sweep.topology_sweep

- 15 of 703 cases; first mismatching field: `values.points` 11, `values.points.len` 4.
- Fork commits:
  - 80d307d rotational sweep: a base curve the face creator let go is not read again
  - c8f6a31 SMLib: sweeps, lofts and extrusions pass on their statuses
  - f98bdf5 build: macOS (Apple Silicon) builds and tests with the documented commands

## sweep.trans_sweep_geometry

- 17 of 788 cases; first mismatching field: `values.geometry` 11, `values.geometry.len` 6.
- Fork commits:
  - c8f6a31 SMLib: sweeps, lofts and extrusions pass on their statuses

Measured in #2512 (WSL, the fork's libraries built with `harness/cpp/CMakeLists.txt`'s flags).
