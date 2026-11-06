# HCR projected endpoint classification in #2840

Historical candidate: `09495e082`; exact main base:
`3f82278a14c1a6104db8e090a59b4bcf677ace01`.
Reference: pinned v1.0.1, `v1-3fd5b9548d19-726643fb17a52b7d`.
Unit: `blend.hcr`, case `gen0392`, generated seed `10620584185765354711`.
This records the historical run, not validation of subsequent affine-span changes.

The case asks for curve visibility on a degree-one square wire of side `0.0002`,
with knots `[0,0,.25,.5,.75,1,1]`, along view `(0,0,1)`. The queried interval is
`[0.014419226325001276, 0.17186077553524884]`. A second square is shifted by
`(0.0001,0,0.0001)`; selector `sil=[-1]` includes every edge, including the original
wire. The approximation tolerance is `0.001`, but that is not the vertex/edge
intersection tolerance used below.

Only the start classification changes. The serialized native record is
`[100,0,-1,1e20,0]` (Unknown); the candidate is
`[104,1,0,0.014419226325001276,0]` (Vertex, vertex 0). The fields are class code,
object kind, object index, pre-snap parameter, and point-class parameter; the final
zero is not a measured distance. All other classification records are unchanged,
including the split at parameter `0.125`.

`blend/hcr.rs::intersect_curve_with_silhouettes` projects each silhouette vertex
onto the queried curve interval using `ProjectedMinimize`, target distance zero,
and `SolutionRequested::All`. The intersection tolerance is the sum of the vertex
and curve zones (`geom/tol.rs::xsect_tol_3d_from_zones`). The wire builder uses the
default `1e-5` zones, giving `2e-5` here.

On this entire interval the curve is exactly `(0.0008*t,0,0)`. Its constrained
closest point to original vertex 0 is the lower endpoint, at projected distance
`0.0008 * 0.014419226325001276 = 1.153538106000102e-5`. This is within `2e-5`,
and the equal z coordinates pass the visibility-depth check. Thus the new Vertex
association is supported by the caller's existing tolerant classification rule.
It does not assert exact geometric coincidence or justify changing any tolerance.

Sources: `crates/difftest/src/units/blend_hcr.rs` (selector, builder and `put_pc`),
`crates/difftest/src/own_builder.rs`, `crates/solidcore/src/blend/hcr.rs`, and
`crates/solidcore/src/geom/tol.rs`. This finding neither parks the mismatch as a
`cpp_defect` nor establishes correctness of other changed callers.
