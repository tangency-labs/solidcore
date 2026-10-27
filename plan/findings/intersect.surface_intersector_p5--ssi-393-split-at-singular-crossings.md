# intersect.surface_intersector_p5: split SSI curves at singular crossings, opt-in (2026-10-27)

Units: `intersect.surface_intersector_p5` (green; the call sits in `do_intersection_in`).
Issue #393, record WR-77
(`plan/upstream-issues/WR-77--addendum-2026-10-27-rust-port-opt-in-split.md`). PR: #1304.

## What was done

`Kernel::fixes.ssi_split_curves_at_singular_points` (off by default) runs
`split_curves_at_singular_points` (`surface_intersector/split_singular.rs`, one call in
`p5.rs` before "set outputs") on the final curve list. It is not a port: there is no C++
counterpart, and it is written from the description in #393 only. A crossing or near
tangent curve is sampled at 200 parameters for the sine of the angle between the surface
normals (from the UV curves); each local minimum is refined by golden section, and the
3D curve and both UV curves are split there when the sine is below 1e-6, rises above 1e-4
at both neighbouring samples and the point is more than ten approximation tolerances from
either end. Pieces are trimmed with a zero zone tolerance (no degenerate-segment removal)
and checked to cover the same intervals.

## Evidence (origin/main aaed3de5e + this branch, ssi-corpus 119 cases, tol 1e-6)

| Rust engine | main | branch, flag off | branch, flag on |
|---|---|---|---|
| ERROR | 2 | 2 | 2 |
| INCOMPLETE | 25 | 25 | 25 |
| WRONG TOPOLOGY | 11 | 11 | **5** |
| PASS | 81 | 81 | **87** |

- Flag off: the report is identical to a clean `origin/main` run (JSON-normalised, timing
  fields dropped).
- Flag on: exactly the six cases of #393 change, WRONG TOPOLOGY -> PASS: `c4_steinmetz`,
  `c4_viviani_figure8`, `c4n_viviani_nurbs`, `y23_fig1a_two_loops_cross`,
  `y23_t6_ex3_two_crossing_curves`, `y23_t3_ex4_torus_ellipsoid`. Every other verdict is
  unchanged. C++ and the `ssi` / `solidcore-ssi` engines are unchanged.
- `y23_t2_ex3_loop_two_curves` stays INCOMPLETE (8% of a branch covered) but its output
  changes: the crossings are split and the missed branch remains. This is the one extra
  parity mismatch of the flag-on run next to the six flips and the existing one on main
  (`c5_cyl_quilt_seam`).
- Unit difftests (`intersect.surface_intersector_p5` 1731, `intersect.adv_surface_intersector`
  1931, `surface.surface_p2` 16688, `surface.plane_p2` 2228, `surface.cone_p3` 1442) give the
  same counts on main and on the branch, 0 deferred, 0 vacuous; `p5` also at `--seeds 6000`.
- `crates/solidcore/tests/issue_393_split_at_singular_crossing.rs`: two equal perpendicular
  cylinders (formulas in the file, no corpus data). Default: 4 curves, all ending at
  (0, 1, 0), none at (0, -1, 0). Flag on: 6 curves, four ends at each crossing, the UV curves
  on the 3D curve along every piece. A cylinder of radius 0.6 in place of one of them
  (no crossing) gives the same curves with the flag on.

## Not covered

- General path only: the analytic special cases of `GlobalSurfaceIntersect` (plane, cone,
  sphere, torus) never reach `do_intersection_in`. The pass also acts in
  `ParallelProjectCurve`, which runs the same intersector; the effect there is not measured.
- A curve with a zero normal anywhere (a pole or an apex) is left whole; two crossings
  closer than a 200th of a curve are not told apart; the thresholds are absolute (1e-6,
  1e-4) and only tested at the corpus' tolerance 1e-6.
- No C++ comparison with the flag on is possible (C++ has the defect). Correctness rests on
  the corpus verdicts (exact ground truth) and the Steinmetz test.
- The cases that stay WRONG TOPOLOGY (5) and INCOMPLETE (25) are not touched.

## Seen on main, not caused by this change

`c5_cyl_quilt_seam` is the one parity mismatch of the flag-off corpus run on main
(`values.curve000.s1_at_uv.len: cpp=0 rust=12003`; both engines PASS).

## Harness

The `SSI_SC_FIXES=<name>[,<name>]` switch is the one of PR #1285, carried here
verbatim (`ssi.rs`, `ssi/engine.rs`, `solidcraft.rs`) plus one `enable_fix` arm.
