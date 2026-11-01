Updates: WR-62

#### WR-62 (addendum 2026-11-01) — Opt-in fix stopping `SmTLSPointVertex::SolveIt` ray-firing fallthrough

- **Severity:** Wrong result
- **Area:** Topology, attributes, solvers
- **Status in the port:** `Kernel::fixes.point_vertex_rayfire_does_not_fall_through`, off by default. Default `LocalSolver::solve_it` (`!k.fixes.point_vertex_rayfire_does_not_fall_through`) preserves pinned v1.0.1 C++ parity so `intersect.topology_solver_p1` and `intersect.topology_solver_p2` difftests remain unchanged.
- **Confirmed:** yes — once `Scene::object_zone_tol_3d` was wired in `crates/difftest/src/units/intersect_topology_solver_p1.rs` (PR #569), `intersect.topology_solver_p1` promoted to `green` (`3532/3532` pass, `0` deferred), confirming upstream's fallthrough against the C++ reference on `point_vertex_rayfire_behind_start` (`edge.json`) and `adv_rayfire_behind_start_single_False`, `adv_rayfire_behind_start_single_True`, and `adv_rayfire_ahead_single_True` (`adv.json`).
- **Location:** `source/SMLib/src/SmTopologySolver.cpp:206-262` (`SmTLSPointVertex::SolveIt`, def. `:116`; inherited by `SmTLSVertexVertex::SolveIt` at `:365-382`) at v1.0.1 (`3fd5b95`); in Rust, `crates/solidcore/src/intersect/topology_solver/p1.rs` (`LocalSolver::solve_it`, `LocalSolver::solve_point_vertex`, `LocalSolver::solve_vertex_vertex`).
- **What changes:**
  1. Upstream `SmTLSPointVertex::SolveIt` (`SmTopologySolver.cpp:206-237`) handles `case SM_SO_RAYFIRE:` by computing `dParam` along `ray`, rejecting `dParam < -SM_EFF_ZERO * dScale`, snapping `dParam < 0.0` to `0.0`, testing `sRayPnt.Distance(pV->GetPoint()) <= dTol`, and assigning `dDistSq = dParam * dParam;` (with `dSign = 1.0`, since `SM_SO_RAYFIRE` has `m_eOperationCategory = SM_OPERATION_OTHER`). Because `case SM_SO_RAYFIRE:` omits `break;`, execution falls through into `case SM_SO_3D_SIGNED_DIRECTED_MINIMIZE:` (`:239-262`), which recomputes `dParam` without the `0.0` snap, negates `dParam` when `m_bSwapOrder` is `TRUE`, overwrites `dDistSq = dParam * dParam`, and sets `dSign = -1.0` when `dParam < 0.0`. Consequently, a vertex in `[-SM_EFF_ZERO * dScale, 0.0)` behind the ray origin gets a small signed distance instead of `0.0`, and swapped-order `Rayfire` calls negate positive ray distances (which in `SmGlobalSolver::AddSortedSolution` also defeats `far_from_best` culling and sorts farther vertices ahead of nearer ones). Neither `SmTLSPointCurve::SolveIt` (`:458-468`) nor `SmTLSPointFace::SolveIt` (`:754-760`) falls through or negates `Rayfire` on `m_bSwapOrder`.
  2. When `Kernel::fixes.point_vertex_rayfire_does_not_fall_through` is `true`, `LocalSolver::solve_point_vertex` (and `LocalSolver::solve_vertex_vertex`, which delegates to it) sets `dist_sq = param * param` from the `Rayfire`-snapped `param >= 0.0`, keeps `sign = 1.0`, and skips the `Signed3dDirectedMinimize` block when `op == SolverOperation::Rayfire`.
  3. All non-`Rayfire` operations (`Signed3dDirectedMinimize`, `Minimize`, `Maximize`, `Intersect`, `AtDistance`, `SignedDirectedMinimize`, etc.) are unchanged under both modes.
- **Regression evidence:** Unit tests in `crates/solidcore/src/intersect/topology_solver/p1.rs`:
  - `wr62_flag_default_off_and_lookup_by_name`
  - `wr62_rayfire_vertex_just_behind_start_snaps_to_zero_when_enabled`
  - `wr62_rayfire_swapped_order_keeps_positive_distance_when_enabled`
  - `wr62_non_rayfire_operations_unchanged`
- **Tracking:** tangency-labs/solidcore#2206.
