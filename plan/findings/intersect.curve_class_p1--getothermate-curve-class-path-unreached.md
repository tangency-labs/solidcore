# intersect.curve_class_p1: no production input reaches SmTopologyIntersector::GetOtherMate on the curve-class path (2026-10-30)

Unit: `intersect.curve_class_p1` (stub at `crates/solidcore/src/intersect/curve_class/p1.rs`,
`get_cci_geometry_in`) and `intersect.curve_class_p5` (stub at
`crates/solidcore/src/intersect/curve_class/p5.rs`, `vertex_edge_edgeuses`). Both are
`UNPORTED(boolean.topology_intersector_p1): SmTopologyIntersector::GetOtherMate`, hit when
the classification carries a topology intersector and the `ClassGeometry` adapter answers
no mate (`topology_intersector_mate` returns `None`, the trait default).

## Question

Wiring these stubs means giving the library adapters that set a classification's
intersector (`boolean::merge::p4::ArenaPair` for `MergeCurveOnSurfaces` in booleans,
`InsertionGeometry` for the fillet insertion) access to the intersector's mate maps. That
plumbing is only worth doing, and can only be tested, if pinned C++ actually reaches
`GetOtherMate` from these two curve-class sites on some input.

## Method

Pinned C++ (v1.0.1, `3fd5b95`) debug build of main `b6186da0e`, under gdb with
non-stopping breakpoints on the two upstream sites (each runs only when
`GetTopologyIntersector()` is non-NULL) and on the line after each that runs only when a
mate is found:

| Site | Upstream | Rust stub |
| --- | --- | --- |
| A | `SmCurveClass.cpp:1212` (`sm_GetCCIGeometry`, cpTI set), `:1222` (mate found) | `p1.rs`, `get_cci_geometry_in` |
| B | `SmCurveClass.cpp:15234` (`IntersectVertices` UV-trim sweep, cpTI set), `:15239` (mate found) | `p5.rs`, `vertex_edge_edgeuses` |

Units run (each `difftest --unit <id>`, its whole case set): every `boolean.*`, `blend.*`,
`offset.*`, `heal.*`, `sweep.*` and `intersect.curve_class_p1`..`p9` unit, 55 in all.

## Result

- Site A: 0 hits in all 55 units.
- Site B: 0 hits in every boolean, blend, offset, heal and sweep unit. Hits only in
  `intersect.curve_class_p5`: 18 with the intersector set, 12 with a mate found. Every
  backtrace is the runner's own topology harness, not a library operation:

  ```
  SmCurveClassification::IntersectVertices (SmCurveClass.cpp:15234)
  (anonymous namespace)::run_vertices (intersect_curve_class_p5_topo.cpp:344)
  sc_intersect_curve_class_p5_topo (intersect_curve_class_p5_topo.cpp:375)
  solidcore_ffi::units::intersect_curve_class_p5::topo
  difftest::units::intersect_curve_class_p5::topo::cpp_checked (topo.rs:856)
  difftest::units::intersect_curve_class_p5::topo::valid (topo.rs:840)
  difftest::units::intersect_curve_class_p5::topo::generate (topo.rs:1550)
  ```

  Those cases do not skip or defer on the Rust side: the runner's adapter
  (`crates/difftest/src/units/intersect_curve_class_p5/topo.rs`, `ClassGeometry` impl)
  answers `topology_intersector_mate` from the same edge mates the C++ intersector
  relates (`GetOtherMate` first B-rep to second, else `GetBrepMate`), so Rust takes the
  same branch as C++ and the mates' deleted UV trim curves are compared
  (`intersect.curve_class_p5` is green with 0 deferred on main). The library path behind
  the stub is not covered there, because no library adapter answers the mate.

In upstream a boolean sets the intersector on a classification only in
`MergeCurveOnSurfaces` when neither `pBrepEdge` nor `pOtherEdge` is given (the `IIRSurfs`
surface/surface curves); both sites then also need a point solution on an edge with
both surfaces set (site A, inside `sm_RefineCCIToSurfaces`), or a vertex refinement on an
interpolated surf/surf/surf point (site B). No fixture or generated case in the units
above produces that.

## Decision

Both stubs stay as documented stubs (team-lead, 2026-10-30): no production input reaches
them, so adapter plumbing would be untested. Wire them when a case (a boolean or fillet
whose intersection curve refines onto a Boolean-mated edge) reaches either site in pinned
C++; rerun the gdb script above on that case first.

Script: two gdb breakpoint commands per site with a hit counter (the `othermate.py` of the
investigating PR); the counts above are from that run.
