# Projected extrema at a plane UV kink

Review [4501904652](https://github.com/tangency-labs/solidcore/pull/2968#discussion_r4501904652)
identifies a genuine extension of the angular-G1 problem to curves on surfaces.
A plane carrying a native UV spline can receive the cache's G1 label at a small
corner. Omitting that parameter can remove the true projected extremum.

The public regression uses exact dyadic values `a=1/1024` and `d=1/524288`.
Its UV controls are `(-a, 1-d)`, `(0,1)`, `(a,1-d)` for the maximum and
`(-a, 1+d)`, `(0,1)`, `(a,1+d)` for the minimum. Knots are
`[0,0,1,2,2]`, degree one. The plane maps UV to `(u,v,0)`, the query is the
origin and the projection normal is Z. Both true join distances are exactly 1.
The test verifies that the existing continuity classifier calls the join G1.

On the maximum's left span, half the derivative of squared distance is
`a^2*(t-1) + d*(1-d+d*t)`. It is increasing and already strictly positive at
`t=0`, since `d*(1-d)>a^2`. The right span is its reflection. For the minimum,
that derivative is `a^2*(t-1) - d*(1+d-d*t)`, strictly negative throughout the
left span; reflection gives the opposite sign on the right. Neither open span
has a stationary point. Both Single and All requests must therefore retain the
join, rather than relying on an interior stationary search to recover it.
The public regression fails against the pre-fix library.

![The same defining objectives, with the formerly omitted join and corrected answer](figures/projected-plane-uv-kink.png)

The plot shows the defining squared objective, not a substituted geometry or
an inferred old returned point. The ring marks the join the old endpoint policy
omitted; the filled dot marks the corrected exact answer. The reproducible
plot helper is `diagnostics/projected_plane_uv_kink.py`.

The fix preserves the native B-spline proof and applies it to a native UV spline
on a plain analytic plane. Positive weights and a UV control hull strictly
inside the plane's analytic domain exclude clamping; the map is then affine,
so an exact native UV tangent identity is preserved. The hull and continuity
facts follow the defining basis, as described in
[MIT's B-spline curve reference](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node17.html).
A Plane label carrying a kept procedural definition does not establish this
identity and cannot omit the join. Private controls cover a genuine kink,
a straight join, artificial subdivisions, the clamping boundary and a kept
corner definition.

This change does not claim an exact smoothness proof for every other curve
wrapper. Other preexisting kind-specific behavior remains outside this narrow
plane-UV repair. The public test is grouped with the existing projected-query
regressions, as requested in
[4496188354](https://github.com/tangency-labs/solidcore/pull/2840#discussion_r4496188354).
Validation: the public regression passes for minimum and maximum in Single and
All modes, returning exactly one solution at parameter 1 with scalar 1. The
private defining-join test passes all five controls. The same public regression
fails against the pre-fix library. Formatting and diff checks pass. Strict
Clippy, compatibility tests and golden replay remain integrated-branch gates.


## Linear wrapper follow-up (review 4503434212)

The same defining corner is also lost through `ProjectedCurve` parallel/rotation
maps and zero-distance `OffsetCurve`. The regression changes the internal knot
to `.5` in the domain `[0,2]`; the geometry and exact derivative proof above are
unchanged after each local affine reparameterization. Baseline maximum queries
return `0.9999985694893212` instead of the exact join value 1. Reversing the offset
moves that join to 1.5, but the original knot API incorrectly still reported .5.

The boundary policy now recursively follows the actual fixed affine projection
map and exact-zero offset parameter map, with the shared 32-level nesting bound.
It does not infer smoothness from the angular continuity label. Offset knot
queries map the requested interval into the base parameterization, reverse the
returned knot/multiplicity order, and mirror the knot values back; continuity
records reverse consistently. Nonfinite, nonordered, or non-invertible mirrored
knots explicitly refuse. For example, a source knot .1 in `[0,1]` cannot be
reported as .9 while claiming its evaluator maps back to the exact defining knot:
`1-(1-.1)` differs from .1 in binary64. This guard does not spend tolerance.

The focused public controls cover four wrappers, min/max, Single/All, reversed
knot multiplicities/continuity order and trimmed knot queries. All 37 projected-query integration tests pass including the roundtrip guard
and stored-label follow-up below; final static/integrated checks remain pending.

A separate broad-retention control remains a confirmed limitation: a straight
parallel-projected spline queried with All and retention 1 returns its true two
endpoint maxima plus a strict minimum and three nonstationary SingleValue points.
The defining right-span objective is `1+((t-.5)/1536)^2`, which is strictly
increasing after .5. Tight retention isolates this wrapper boundary repair and
returns only the two endpoints; it does not repair or justify the broad-retention
classification. The independent reproducer is retained for separate issue work.

This checkpoint does not yet resolve the remaining stored-spline analytic labels
or nonlinear composition policy. In particular, ellipse/hyperbola/parabola
labels evaluate their stored spline and require a definition-based proof; a
label is not itself proof of a smooth native join.


## Evaluator-backed spline labels

The public database readers for `Ellipse`, `Hyperbola` and `Parabola` accept a
stored spline. Their actual point/derivative, knot and continuity dispatchers
use that spline. A public reader regression loads the same degree-one corner
above into each label, independently verifies the evaluated join `(0,1,0)`,
and queries the projected maximum. Before the change all three return only
natural endpoints with value `0.9999985694893212`; the corrected All answer is
the single exact join at .5 with value 1. This does not infer geometry from the
analytic metadata or change database parsing.

Those evaluators now share the native defining-knot test. Rational C0 joins need
more than the previous equal-weight chord shortcut. Write adjacent homogeneous
controls as `(A,wa),(B,wb)` on the left and `(C,wc),(D,wd)` on the right. Exact
position continuity is `B*wc=C*wb`. Up to positive degree/knot/weight factors,
the two tangent numerators are `L=B*wa-A*wb` and `R=D*wc-C*wd`. Each component of
`L cross R` expands into eight signed products of four stored binary64 numbers.
The test proves these sums exactly zero and proves a common nonzero component
has the same sign; opposite directions and unresolved zero tangents are retained.

The existing fixed integer accumulator is generalized, not duplicated. Original
two-factor sums retain 66 words; four-factor sums use 132 words. Four finite
binary64 significands require at most 212 bits, and their summed exponent is in
[-4296,3884]. Eight terms need at most 8395 bits including carry, within 8448.
The only accepted generic sizes are compile-time `(2,66)` and `(4,132)`; no
allocation or recursion depends on a caller-supplied size. Independent Python
Fraction arithmetic agreed on 512 deterministic mixed-sign, extreme-exponent,
permuted-cancellation and subnormal residual cases. Private controls also exercise
full eight-term carry accumulation and one-ulp weighted kinks.

All 37 projected-query integration tests pass at this checkpoint, including the
linear wrappers, non-invertible reflection refusal and all three stored labels.
The complete original p1 private suite plus new arithmetic controls passes all
22 tests. Final static/integrated checks remain pending. The nonzero-offset, perspective and
nonlinear surface/volume composition policies remain separately unresolved;
this evidence must not be described as closing their regularity contract.
