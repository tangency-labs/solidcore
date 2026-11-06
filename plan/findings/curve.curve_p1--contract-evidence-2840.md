# Projected-query contract evidence and remaining design decisions

Follow-up to PR #2840, candidate `2f091cce4`. This separates questions answered
by the implementation from decisions about replacing its public contract.
No production code or mismatch classification changes in this addendum.

## Zero distance tolerance: existing behavior is answerable

The legacy projected solver does not interpret `dist_tol = 0` as a request for
an exactly representable mathematical optimum. In `geom/global_solver.rs`,
`SolveRequest::tolerance_3d` controls pruning and solution retention. The
projected stationary search in `curve/curve/p1.rs` calls
`curve/curve/p5.rs::local_property_analysis` without passing that tolerance.
That function independently configures `CurvePropertyEfo` and the local solver
with `EFF_ZERO_SQ`, and accepts convergence using `EFF_ZERO_SQRT`.

Therefore the historical meaning need not be left as an unanswered question:
zero retention slack and nonzero numerical solve thresholds coexist. The new
affine path additionally uses this argument to bound the optimization gap.
The refusal in `intersect.topology_solver_p1::gen1646` is consequently a contract
change, not proof that the caller requested exact arithmetic.

The remaining design decision is the numerical accuracy contract for the new
certified path. Retention slack and optimization accuracy should be distinct;
the latter's units and acceptance rule must be specified. Reusing a Newton
stationarity threshold as a distance-gap allowance is dimensionally unjustified.
Setting a positive floor only for zero, or deleting the certificate, would not
resolve that design. No such shortcut is implemented.

## Cached search domain: existing behavior is answerable

`global_point_solve_cached` builds the curve cache independently of the requested
interval. `CpGlobalSolver` then visits its existing leaves, rejects disjoint
leaves, and intersects overlapping leaves with the requested interval before
the affine dispatch. The dispatch does not create extrapolated cache spans.
Thus this path searches the intersection of the request with the available
cached spans; it does not establish that every parameter evaluable outside
those spans is globally searched.

This explains the clipped constant range in `curve.curve_p1::gen2883`. The
restriction precedes the affine path and is not a new domain-clipping policy
introduced by that path. Extending the global search beyond its cache would be
a separate behavior change. No additional extrapolation implementation is needed
to explain this case, and it should not be presented as an unknown local fix.

## Finite trimming and silhouette accuracy

The executed finite-segment/arc vote replays already explain the three trim
outputs through the existing whole-face voting and tie rules. They provide no
proof that replacing those rules is required by the small-elbow fix. A new
finite-trim interpretation remains a separate modeling-contract decision.

The silhouette accuracy defect is a concrete technical limitation, not a
question about what the reported error is: the sphere-band result exceeds its
requested tolerance. `intersect/topology_solver/p1.rs::face_silhouettes` builds
approximate silhouette curves with `approx_tol = tolerance_3d * 10.0`, then
the projected query measures those curves. An endpoint fix cannot certify the
original surface through that approximation. Merely reducing that multiplier
would not prove a bound; the construction needs a bounded error budget or a
solve on the surface definition. Agreement is still needed on ownership/scope
of that construction before adding it to this PR.

## Review status

The broad statement that all of these are unanswerable without a human was too
strong. The existing zero-tolerance behavior, cached-domain restriction and
trim outputs can be established from code and the recorded executions.
What remains is approval of a changed numerical API and the scope of the
surface-construction repair, plus the individually documented native differences.
The full validation record still applies; this addendum does not make the
candidate ready to merge.
