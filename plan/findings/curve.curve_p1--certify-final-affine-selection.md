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


## Ordinary oblique trim queries: centered carrier proof

The parent hosted sweep exposed three `primitive.api_trim_ops` refusals:
`gen0017`, `gen0021`, and `gen0137`. Each uses an ordinary analytic line cutter.
The first two APIs normalize their direction; `trim_project_parallel` in the third
preserves its supplied nonunit direction. All reference queries succeed, but an
original box vertex or surface corner fails the minimum certificate. The required
caller-status propagation exposes this error rather than silently classifying a
face using an uncertified projection.

The actual analytic Line's stored anchor, direction, scale and STEP/natural map
were dumped, not reconstructed from evaluated endpoints. Exact Fraction arithmetic
proves the minimum and each candidate lie inside the candidate's existing own
point enclosure. The old carrier expression's lower endpoint missed that own lower
endpoint by one or two floating-point spacings: these are avoidable refusals, unlike
the narrow-knot unrepresentability controls.

| Query | Old carrier lower | Unchanged candidate own lower | Centered carrier lower |
|---|---:|---:|---:|
| gen0017 original origin vertex | 4.717646518721153 | 4.717646518721154 | 4.717646518721160 |
| gen0021 original origin vertex | 5.359602813902480 | 5.359602813902481 | 5.359602813902484 |
| gen0137 upper-left corner | 86.59806893767195 | 86.59806893767198 | 86.59806893767202 |
| gen0137 upper-right corner | 85.56319489742386 | 85.56319489742388 | 85.56319489742390 |

For the double-cross projection P, set X=P(C(t)-Q) at the existing stationary
fraction estimate and let D be any nonzero projected defining direction. Completing
the square gives the full infinite-carrier minimum squared:

`min_s |X+sD|^2 = |X|^2 - (X.D)^2 / |D|^2`.

This is the one-variable least-squares normal-equation identity (Boyd and
Vandenberghe, [Convex Optimization, section 1.2.1](https://web.stanford.edu/~boyd/cvxbook/bv_cvxbook.pdf)).
Directed arithmetic encloses all defining operations; a positive denominator is
required. The infinite carrier contains the entire requested continuous span, so
its minimum is a valid lower bound even when the stationary point is outside the
requested interval. At an isolated upper knot, the existing separate point bound
is still combined by minimum. This new independent lower bound is combined with
the existing carrier and hull bounds. No unsuccessful construction selects a
fallback, no candidate enclosure is widened, and retention slack is not used as
solve accuracy.

The [independent Fraction checker](diagnostics/projected_trim_affine_certificate.py)
uses [actual defining coefficient data](diagnostics/projected_trim_affine_coefficients.jsonl).
Its acceptance checks use exact squared rational inequalities; decimal square roots
are display only. It verifies that the centered lower bound is below the exact
minimum, that the unchanged own enclosure contains the exact minimum and candidate,
and that all candidate parameters are admissible. The measured candidate excess
ranges from approximately 2e-36 to 5e-29.

Four public regressions cover both normalized directions and both nonunit corner
queries, for Single/All and zero/positive retention. All four fail with 1021 against
the preserved parent before this correction. The corrected construction passes all
36 public projected-query tests, the side-classification matrix (12 combinations),
and all eight caller-keyed operation tests. The exact Fraction checker passes all
12 candidate checks, and all 20 private p1 tests pass. Strict library/tests Clippy
and workspace formatting also pass. Scoped golden validation remains pending at
this checkpoint; no golden output is changed.
