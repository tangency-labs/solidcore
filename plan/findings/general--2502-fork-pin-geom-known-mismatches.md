# geom: cases that move with the difftest pin to the SMLib fork (2026-11-03)

Units: geom.assert_array, geom.extent1d, geom.extent2d, geom.global_solver, geom.local_solve1d, geom.local_solve_nd, geom.point_set, geom.pseudo_box, geom.transform.

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

## geom.assert_array

- 4 of 3735 cases; first mismatching field: `ints.message.len` 4.
- Fork commits: none on its own or its dependencies' upstream files.

## geom.extent1d

- 5 of 3091 cases; first mismatching field: `values.values` 5.
- Fork commits:
  - 922f333 periodic interval: a value just past the period's top is a period down, not the interval's start

## geom.extent2d

- 1 of 3489 cases; first mismatching field: `values.periodic_point` 1.
- Fork commits:
  - 19310fa SMLib: composite curves, gaps, pseudo boxes, point sets and the geometry utilities pass on statuses

## geom.global_solver

- 2 of 628 cases; first mismatching field: `values.floats.len` 2.
- Fork commits:
  - 709c5be SMLib: the perspective cull projects each sphere's own radius from its centre
  - a6d8293 SMLib: the solvers drop no status

## geom.local_solve1d

- 26 of 3318 cases; first mismatching field: `values.state_f` 26.
- Fork commits:
  - a6d8293 SMLib: the solvers drop no status
  - d2757e1 SMLib: the local solvers take no value from a failed evaluation

## geom.local_solve_nd

- 391 of 3043 cases; first mismatching field: `ints.found` 336, `values.trace.len` 32, `values.trace` 22.
- Fork commits:
  - d5c5cfd SMLib: an N-D solve whose guess can't be evaluated finds no solution
  - d2757e1 SMLib: the local solvers take no value from a failed evaluation

## geom.point_set

- 1098 of 2032 cases; first mismatching field: `values.c00_tight` 251, `ints.c00_contain` 202, `values.c01_tight` 125.
- Fork commits:
  - 19310fa SMLib: composite curves, gaps, pseudo boxes, point sets and the geometry utilities pass on statuses
  - ca7d1f5 capi: no shim drops a status SMLib returns

## geom.pseudo_box

- 33 of 3012 cases; first mismatching field: `values.c016_intervals` 3, `values.c015_intervals` 3, `values.c018_intervals` 3.
- Fork commits:
  - 19310fa SMLib: composite curves, gaps, pseudo boxes, point sets and the geometry utilities pass on statuses
  - ca7d1f5 capi: no shim drops a status SMLib returns

## geom.transform

- 26 of 3042 cases; first mismatching field: `values.r` 26.
- Fork commits:
  - 19310fa SMLib: composite curves, gaps, pseudo boxes, point sets and the geometry utilities pass on statuses
  - ca7d1f5 capi: no shim drops a status SMLib returns

Measured in #2512 (WSL, the fork's libraries built with `harness/cpp/CMakeLists.txt`'s flags).
