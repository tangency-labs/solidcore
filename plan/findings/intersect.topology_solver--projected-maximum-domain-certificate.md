# Projected point/face maximum: whole-domain certificate

The reported operation is a single projected point/face maximum. Its candidate
curves now evaluate `S(uv(t))`; an independent fitted spatial curve is not measured.
This removes the original extra fit error but is not, by itself, a certificate.
The verified rectangular boundary curves are also candidate curves, so a
boundary maximum does not depend on a silhouette reaching it or on the global
tree visiting edges first.

## Construction

For `Q(x) = N cross x cross N`, bound `f(u,v) = |Q(S(u,v)-P)|^2` on every
nonempty knot rectangle of the face domain. Existing outward-rounded NURBS
partials supply intervals for `S`, `S_u`, and `S_v`. Intersect the direct squared
norm enclosure with the mean-value enclosure

`f(c) + 2 Q(S-P).Q(S_u) (U-c_u) + 2 Q(S-P).Q(S_v) (V-c_v)`.

The centre evaluation is also an interval evaluation. The rectangle is convex,
so the mean-value formula encloses every value without assuming a stationary
point was found. Knot rectangles are checked separately. Subdivision is
iterative and deterministic, with cancellation, a 16,384-box budget and a
48-bisection path limit. A failed denominator enclosure is subdivided; a box
that remains unresolved is never accepted.

Every box must prove its distance upper bound at most the reported scalar plus
the caller's tolerance, with that threshold rounded inward. This accounts for
missed branches without certifying or replacing the tracer. The method follows
the mean-value inclusion principle of interval analysis; see equation (27) in
[Neumaier, Complete Search in Continuous Global Optimization and Constraint Satisfaction](https://arnold-neumaier.at/glopt/mss/Neu04.pdf),
section 11, and
[Patrikalakis, Maekawa and Cho, interval arithmetic](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node45.html).

## Definition and feasibility

Plain NURBS definitions and the existing checked analytic exactness classes are
eligible. Procedural blend/sweep definitions and unsupported surface classes
refuse; their displayed net cannot substitute for their definition.

A newly returned face candidate must be a single point within the UV rectangle.
The stored trim boundary must consist of one connected four-edge loop whose
edges are full, clamped, degree-one, two-control-point UV segments or analytic
lines joining the four exact rectangle corners. Analytic lines are checked from
their defining point/vector/scale/parameter interval with the existing exact
binary64 product predicate; neither their stored net nor rounded endpoint samples
prove that identity. The tolerance-based rectangular flag is never
used as proof. The metric at the returned UV is enclosed independently, and the
entire enclosure must agree with the reported scalar within the tolerance.
A previous answer on another object is retained as the caller's already
established feasible bound; this local certificate proves that this face cannot
improve it by more than the tolerance.

Failure restores the solver state from before candidate generation and returns
`Unsupported`, rather than publishing an unproved candidate. The certificate
uses no enlarged tolerance or sampled accuracy floor.

## Scope

This change is deliberately restricted to `PointFace`/`VertexFace`,
`ProjectedMaximize`, `Single`. Curve/face and face/face queries, other operations,
and requests enumerating all extrema retain the parent candidate implementation.
The shared projected minimum/maximum topology-box culling now encloses the same
double-cross metric, including non-unit normals. The former projected sphere
test used an unscaled plane distance and could discard the face after seeing a
vertex value scaled by `|N|^2`; the public traversal test exposed this at normal
length 1000. Outward-rounded box bounds replace that mismatched comparison.
The pruning comparison keeps the caller's tolerance around the current best,
so requests for all answers do not discard an answer within that tolerance.
The regression at half the tolerance fails without this margin; an answer two
tolerances away is still culled for both projected minimum and maximum.
Earlier
child-branch pair-query measurements are investigation results, not evidence
that this certificate covers those operations.

Focused validation covers local and public
B-rep traversal for the analytic fixtures. The sphere result is
`0.860555127546399068` against analytic `0.860555127546398957`; the original
inaccurate scalar is refused by the independent whole-domain bound.
The full-suite and golden results are recorded in the PR description.

![The same sphere query before and after: the old reported maximum lies outside the caller's tolerance, while the new value agrees with the analytic maximum.](figures/projected-maximum-comparison-2967.png)

The upper panels show the same analytic projection and maximizing point; the
lower panels magnify the exact distance function near its maximum and overlay
the measured scalar. The drawing is explanatory evidence, not the certificate.
The original scalar is recorded in issue #2967; the new scalar is from the
regression output. The plotting source is
`diagnostics/projected_maximum_comparison_2967.py`.
