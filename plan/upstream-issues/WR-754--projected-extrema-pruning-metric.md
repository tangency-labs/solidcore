#### WR-754 — Projected extrema prune with bounds in different units

- **Upstream:** `source/SMLib/src/SmGlobalSolver.cpp`,
  `SmGlobalSolver::BranchMayContainAnswers`, projected min/max sphere bounds.
- **Objective:** curve-point projected distance is `|N cross V cross N|`.
  With `q=N dot N`, this equals `q` times orthogonal projected distance.
  The sphere-center distances and radii used for pruning are unscaled, so
  they do not enclose that objective when `q != 1`.
- **Identical-input native reproduction:** two leaf boxes, a point at the
  origin and the segment from `(2,0,7)` to `(4,0,7)`, zero tolerance,
  Single request. The second box's circumsphere has center `(3,0,7)` and
  radius `1`. The shim calls the native `BranchMayContainAnswers` method;
  it does not copy or replace that method.

  | Operation | Normal | Incumbent distance | Actual improving distance | C++ may contain | Rust may contain |
  |---|---|---:|---:|---|---|
  | projected min | `(0,0,0.5)` | `0.625` | `0.5` at `(2,0,7)` | false | true |
  | projected max | `(0,0,2)` | `14` | `16` at `(4,0,7)` | false | true |

  Incumbents are supplied squared (`0.390625` and `196`) to the solver.
  In the first case the old lower bound is `2`, instead of `0.5`; in the
  second the old upper bound is `4`, instead of `16`. Both discard a branch
  containing a better answer.
- **Confirmed:** difftest with the unmodified, checksum-verified pinned
  v1.0.1 reference `v1-3fd5b9548d19-726643fb17a52b7d`. Both native booleans
  are `0`; the corrected Rust booleans are `1`. The inputs are retained in
  `fixtures/geom.global_solver/issue_2840_projected_bounds.json`.
- **How found:** admitting previously discarded endpoints for the small-profile
  mitre follow-up (#2840) tightened incumbents and exposed incorrect pruning
  of better interior extrema with nonunit normals. Independent rational-curve
  regression cases `gen0002` and `gen2419` exercise that interaction.

This record concerns the shared global solver's projected extrema bounds, not
the separate topology pre-pruning code or unsupported surface-pair operations.
