# Certify the final affine projected answer across native spans

The final scalar accumulator can replace a locally certified answer with a
different point whose own enclosure does not contain the global optimum.
In addition, a geometrically degenerate defining span can be absent from the
curve cache while containing the largest projected distance. The cache is not
a complete partition of the defining objective.

The public regressions in `projected_curve_queries.rs` reproduce both failures.
For the minimum, a large-coordinate segment crosses the query at parameter
`1/3`, but its rounded report loses to a separate constant segment at distance
`1e-9`. For the maximum with normal `(1, 0.1, 0)`, the first degenerate span's
exact binary64-input objective is approximately `5.57880165459372927e-10`, while
the wrongly selected second span is only `1.00498756211208906e-10`.
`diagnostics/projected_affine_representability.py` computes these values with
exact rational inputs and high-precision decimal square roots.

The construction selects analytic lines and complete native equal-weight
degree-one spline spans before solving. It visits their defining intervals
directly, including degenerate spans, rather than treating cache leaves as
coverage. Incomplete coverage is not selected; an unexecuted selected traversal
cannot return empty success. Other curve definitions retain their construction.

Each returned candidate keeps its own native evaluation enclosure. Separately,
directed arithmetic on defining factors bounds the whole objective: the norm
of an affine vector is convex, so endpoint bounds enclose its maximum, and the
component-wise endpoint hull lower-bounds its norm. The latter tightens, but
never replaces with a larger unproved bound, each local minimum certificate.
Exact parallelism proves a constant projected vector, so its anchor suffices;
an isolated right-valued knot is always enclosed separately. Neither retention
slack nor another parameter's enclosure width supplies optimization accuracy.

Single preserves the canonical answer when certified and otherwise selects a
certified recorded candidate consistent with the external cutoff. All preserves
the existing canonical identity/range-merging policy only if every potentially
retained recorded candidate remains represented. An ambiguous omission refuses.
The regressions include an omitted true-zero minimum beside a retained zero
range, and an external cutoff that would otherwise discard every rounded report.

The same directed span lower bound is applied before the local minimum proof.
This prevents an avoidable refusal of the ordinary segment from `(-100, .5, 100)`
to `(100, .5, 100)` queried from `(0, 10, 100)` along `-Z`: its exact minimum is
`9.5` at the representable midpoint. The public regression grades the returned
native objective, allowing the existing canonical rounded-distance tie ordering.

Validation: all 31 public projected-query tests pass without changed assertions
or tolerances. The independent exact-rational diagnostic also passes. Strict
Clippy for the library and public regression target passes. Private tests,
affected golden replays and integrated gates remain to be
recorded before this revision is review-ready. No golden output is changed here.
