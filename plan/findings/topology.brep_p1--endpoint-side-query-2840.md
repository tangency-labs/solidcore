# Endpoint side query recovered by projected extrema (#2840)

Unit: `topology.brep_p1`; existing fixture `find_side_6` in
`fixtures/topology.brep_p1/constructed.json`.
Candidate `09495e082`; main base `3f82278a14c1a6104db8e090a59b4bcf677ace01`.

The cubic Bezier has control points `(0,0,0)`, `(1,2,0)`, `(2,-1,1)`,
`(3,0,0)`, unit weights and domain `[0,1]`. The query point is `(-1,0,0)`,
projection direction `(0,0,1)`, and side tolerance `1e-6`.

The projected curve is
`(3t, 6t-15t^2+9t^3)`. Its squared distance from the point is
`(1+3t)^2 + (6t-15t^2+9t^3)^2`, uniquely minimized at `t=0` on `[0,1]`.
The projected tangent there is `(3,6,0)` and the point offset is `(-1,0,0)`;
their cross product dotted with the projection direction is `6 > 0`.
Thus the point is on the left and the API's `right` boolean is false.

The pinned native reference returns error status `1001` and side `UNSURE` (`2`).
The revised Rust library returns success (`0` in the harness) and false (`0`),
matching the independent side calculation. The same endpoint-discard cause
prevents the small-profile elbow trim from classifying its reference point;
see [WR-753](../upstream-issues/WR-753--projected-curve-endpoints-discarded.md).

This is an intentional native disagreement, not a parity pass. No side-query
implementation or tolerance was changed.
