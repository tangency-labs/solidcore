Updates: WR-697

#### WR-697 (update, 2026-11-03). A blend whose rail ends all lie on End planes is traced from the edge's End sections

A rail end on an End's plane (an extended-edge-use end, `SM_FV_RAIL_X_EXTENDED_EDGEUSE`) gives `CalcFilletGeom` no point, as upstream: it is found on the surface once that is made (tangency-labs/solidcore#2500).

- **Defect:** a blend both of whose ends are such Ends, with no G1 knot between, had no point at all, and the trace failed with `SM_ERR` (fillet code 4). WR-641's `fillet_through_extended_rail_end`, now removed, had made the edge's end parameters through points on the old construction.
- **Change (with `fillet_blend_corner`):** `c2_trace_to_edge_ends` (`crates/solidcore/src/blend/fillet_solver/p4.rs`) starts such a trace from the edge's sections at its extended ends (`c2_edge_end_section`), solved as it solves an end the trace stops short of; an end past a tangent edge is skipped as there.
- **Case:** #2175's chamfer 0.5 of a ball's rim, from the top face to its tangent end on the bottom face: refused again with `SM_FILERR_NO_ROOM` (8) by the smooth-vertex End, its rail rolling over the bottom face's tangent edge (no radius builds it). Test: `crates/solidcore/tests/issue_2175_chamfer_rail_past_tangent_rim_end.rs`.
