# Projected pair maxima: All contract and remaining construction (#3060)

`SolutionRequested::All` is not a request for every stationary extremum.
`geom/global_solver.rs::GlobalSolve::add_sorted_solution` rejects projected
maximum candidates below its incumbent minus tolerance, sorts descending, and
trims trailing maxima outside the best-value tolerance band. Only `Single`
additionally truncates to one answer. At `facfa1fe2`, the enum is at line 312; candidate rejection is at
2058, Single truncation at 2280, and best-band trimming follows at 2286.
The enum documents All as every solution within tolerance. Consequently a near-zero stationary answer is not valid when
an admissible pair exceeds it by more than the requested tolerance.

The exact original All inputs and standalone public-API regressions are in
`diagnostics/projected_pair_all_cases_3060.json` and
`diagnostics/projected_pair_all_regressions_3060.rs`. They are deliberately outside
automatic test discovery while the All construction remains incomplete. A
standalone test run against the pre-fix implementation produced five failures,
all at the feasible-witness domination check:

| Case | Returned value | Admissible witness | Tolerance |
|---|---:|---:|---:|
| gen0472 | 5.010571561172592 | 5.346613700814646 | 0.00001 |
| gen1839 | 0.000000000004244 | 4.408015599603469 | 0.0001 |
| gen2050 | 0.000000061371359 | 40.001231972795864 | 0.0001 |
| gen2187 | 0.099557597577128 | 6.310837322899526 | 0.000001 |
| gen2585 | 0.000128676797843 | 0.996226875922088 | 0.5 |

For gen0472, the ellipse constructor's domain is STEP degrees, but solve
parameters are NURBS coordinates. The natural domain is [0,1], contained in the
requested [-12.067776972579878,65.71806952381752]. STEP angles 0 and 90 convert to
0.470498612714079556 and 0.937810025781023349 respectively. Thus the original
witness is admissible; comparing 90 directly to the requested interval is wrong.
This is a lower-bound witness, not a claim that this pair is globally maximal.

gen1839 is a rational quadratic with positive weights. On the union of the two
requested intervals its x derivative numerator has strictly negative quadratic
Bernstein coefficients, and y has strictly positive coefficients. The planar
distance therefore attains its maximum at opposite interval endpoints. The
x numerator Bernstein coefficients are approximately (-6.347414793612,
-7.965986112795,-8.948873955546); y coefficients are (8.012257386824,
5.978103274711,3.708957842806). These are coefficients of N'W-NW'
after restricting the quadratic to the union interval, using the stored
homogeneous coordinates directly. The larger cross corner is
(0.43561427937496744,0.03474215986040488), with distance
4.408015599603469; the other is approximately 3.942004415785032 and lies outside
the tolerance band. The exact-sign calculation is reproducible with
`diagnostics/projected_pair_monotone_witness_3060.py` (standard library only;
`Fraction` preserves the exact stored binary64 controls).

gen2050 lies entirely inside a single degree-one knot span for each operand.
gen2187 and gen2585 are bounded affine lines. Squared projected distance is a
convex quadratic on their parameter rectangle, so its largest value occurs at
a corner. The distinct constrained corner maxima are:

| Case | Larger cross corner | Other cross corner | All requirement |
|---|---:|---:|---|
| gen2050 | 40.001231972795864 | 40.00123193297711 | retain both |
| gen2187 | 6.310837322899526 | 6.300451372620962 | retain larger |
| gen2585 | 0.996226875922088 | 0.983409688864426 | retain both |

The diagnostic checks distinct representatives by injective matching using the
sum of projected endpoint displacements, bounded by the caller's tolerance.
The triangle inequality bounds the corresponding objective perturbation. It
also checks both requested and natural parameter domains and recomputes each
reported scalar. It does not impose an unrelated parameter-space epsilon.

## Construction boundary

The Single construction's exact-definition interval jets, domain clipping, and
candidate evaluation can be shared. Best-value certification alone does not
establish completeness of an All answer. The initial separate All helper can
handle one proven affine piece per operand: classify corners using interval
inward derivatives, bound the global maximum by the four outward corner value
bounds, and retain every certified constrained maximum in the tolerance band.
The initial NURBS class additionally requires strictly increasing active knot
spans: repeated knots can expose a different endpoint value from a one-sided
affine limit, invalidating a bound using only default endpoint evaluations.
An external target seeds the incumbent band; rounded subtraction certifies
inclusion/exclusion, and a proved cutoff miss returns an empty result.
As accepted for p1, retention slack is not arithmetic solve accuracy: zero slack
permits the candidates' own outward point-enclosure widths, never an arbitrary
floor or a geometric box width. Existing `GlobalSolve::add_sorted_solution`
applies reported-value retention and ordered duplicate comparison, including its
scaled natural-domain parameter tolerance and preferred duplicate rule. There is
no swapped-parameter self-pair canonicalization. Uncertain corner classification
or external-cutoff membership refuses rather than discarding a possible answer.

A structural selector chooses this path before certification: native lines or
clamped, equal-weight degree-one native B-splines, finite positive clipped
intervals, strictly increasing active knots, and one intersecting knot span.
Selected proof failures propagate. Other geometric classes remain on their
unchanged implementation; there is no catch-and-fallback path. In particular,
the checked-in `circle_projected_max` fixture (offset circles, isolated maximum
4.5) is outside this partial construction, as are the original rational and
ellipse All cases. Curved pairs need a complete boundary/interior stationary
root construction or a separately proved recognized class in addition to a
best-value bound.

`SolutionType::RangeOfValues` stores two ends in joint parameter space;
`Solution::classify_solution_end` and `GlobalSolve::combine_min_max` do not prove
that an arbitrary unresolved parameter rectangle or nonlinear correspondence is
represented by those ends. Continuous maxima must not be replaced by sampled
points or unproved ranges. The maintainer question is how such continua should
be represented, and whether typed refusal is the intended behavior when that
representation is unavailable.

`SolutionRequested::Nodes` exists, but the pair solve/global solver contains no
distinct request dispatch implementing it. There is no `AllNodes` enum variant.
This finding does not silently assign either a new contract. The topology
wrapper's unconditional internal All request is a separate issue: narrow
forwarding of an outer Single projected maximum belongs with the Single fix.

The distinction between an epsilon-optimal value certificate and complete
stationary-root isolation is also explicit in Neumaier's *Complete Search in
Continuous Global Optimization and Constraint Satisfaction*, sections 12 and 15
([primary paper](https://arnold-neumaier.at/glopt/mss/Neu04.pdf)). Interval
existence/uniqueness tests can isolate finite roots; nonisolated sets need a
different covering/representation argument.
