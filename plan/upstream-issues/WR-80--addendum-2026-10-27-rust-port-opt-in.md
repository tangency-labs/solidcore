Updates: WR-80

#### WR-80. Addendum 2026-10-27: where the NULL pointer comes from, and the Rust port's opt-in fix

- **Location (traced):** `SmSurfaceIntersector::TraceCurve` returns `SM_SUCCESS` with an empty
  `m_vCurvePoints` on two paths: a `SM_IP_TANGENT_POINT` / `SM_IP_COINCIDENCE` start point that
  is already on a found curve (`SmSurfaceIntersector.cpp` ~4979-4993; `m_bCurveIsClosed` is not
  written, so the previous trace's FALSE survives), and a bounded start point whose first step
  leaves the domain while it lies on a found curve (~5017-5023, `m_vCurvePoints.Init()`). The
  "trace the other way too" blocks of `DoIntersection` (~6555/6565) and
  `SmAdvSurfaceIntersector::FindInteriorCurves` (~1346/1356) test `!m_bCurveIsClosed` and call
  `ReverseCurveDirection` (~7272), which is `NER` on an empty list: `SM_ERR_NULL_POINTER`.
  `DoPointIntersection` (~7001/7012) is guarded by `GetLastNode() != NULL`.
- **Trigger in the Rust port (usd-brep 1.0.1 pin, solidcore `origin/main` 4d7ef7021):** both
  ssi-corpus cases of this record, `c2n_quartic_contact_parabola` and `y23_t2_ex5_two_loops`,
  fail the same way in C++ and Rust (`SM_ERR_NULL_POINTER`, parity identical). The boundary
  start-point loop of `AdvSurfaceIntersector::do_intersection_with_caches` meets a cluster of
  tangent start points about 1e-3 apart on the boundary (c2n), or two coincidence points at the
  singular crossing of the two arms (y23): one whose trace comes back empty is reversed.
  A hand-built pair of quartic-contact Bezier patches reproduces it without corpus data
  (`crates/solidcore/tests/issue_391_empty_trace.rs`).
- **Rust port:** reached (this record's earlier "not reached yet" is out of date). The port
  keeps the defect by default. `Kernel::fixes.ssi_empty_trace_not_reversed` skips the reverse
  trace of an empty curve in the three unguarded blocks (`surface_intersector/p4.rs`,
  `p5.rs`, `adv_surface_intersector.rs`), as `DoPointIntersection` already does.
- **What the fix does not do:** with it on, the two corpus cases no longer fail but return no
  curves (`y23_t2_ex5_two_loops` INCOMPLETE, `c2n_quartic_contact_parabola` WRONG TOPOLOGY).
  Finding the curves needs the tangent-curve fixes of the other SSI records (WR-74, WR-75,
  WR-76, WR-79).
- **Confirmed:** yes, by the C++ build (both cases) and by the Rust port with a debug
  backtrace at `reverse_curve_direction` (empty curve, `curve_is_closed = false`).
