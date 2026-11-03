Updates: WR-650

#### WR-650 (update, 2026-11-03). A trace on offsets extended as themselves is bounded both ways

Extending a surface of revolution round its turn, past the closure, lets two offsets meet in a closed curve along any edge, not only along an arc (tangency-labs/solidcore#2500).

- **Defect:** upstream's `CalcFilletGeom` gives the intersector a stop distance (`SetExtensionDistance`, 5 x the offset) only along an arc, and clips only the extension after the end. Under SolidCraft's fixes, #1945's round of a bar's top front edge, run over a concave round's torus up to the boss it rounds into, meets the torus's offset in a closed oval. The extension before the End's section ran round the oval past the trace's own span: the end-point and loop checks (`sm_IsPointOnCurve`) compare UV too, and the wrapped extension does not repeat the UV. The fillet surface covered its span twice, the End's rail ends dropped onto the second sheet, and its cap edge on the boss failed (1001, "CalcFilletIntSideFace").
- **Change (with `offset_extension_as_itself`):** `c3_make_intersector` (`crates/solidcore/src/blend/fillet_solver/p4.rs`) gives every extended trace the stop distance, and sets `SurfaceIntersector::clips_extension_before_start`, which stops `DoPointIntersection`'s trace back from its first point past that distance from it. Without the flag nothing changes.
- **Tests:** `crates/solidcore/tests/issue_1945_rail_uv_trim.rs` (one lump, `check_brep` clean, the removed volume within its closed-form bracket, the round on its rolling-ball definition) and `extension_distance_both_ways_with_offsets_extended_as_themselves` in `p4.rs`.
