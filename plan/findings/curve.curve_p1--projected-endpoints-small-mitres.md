# Projected endpoint extrema and small-profile mitres

Related: #2840; the available geometry follow-up in #2811, under #2794.
Units: `curve.curve_p1`, `curve.curve_p4`, `geom.global_solver`;
caller `primitive.primitive_creation_p2`.

## Reproduction and geometry

Sweep an XY circle of radius `r` along `(0,0,0) -> (0,0,L) -> (L,0,L)`,
with both ends capped and the profile moved to the path. Before the fix,
`(r,L,tolerance) = (0.001,0.01,1e-7)` and `(0.001,1,1e-7)` return success
but fail Level2 validity: two vertices have incident-surface gaps
`0.00182842712474619`. The `(0.01,0.1,1e-5)` case has fourteen validation
failures; `(0.1,1,1e-7)` is a valid control.

The exact cylinder equations are `x^2+y^2=r^2` and
`y^2+(z-L)^2=r^2`. Their mitre plane is `x+z=L`; the seam is
`(r cos(t), r sin(t), L-r cos(t))`. The volume is `2*pi*r^2*L`.
This is a nondegenerate 45-degree cylinder/plane cut at every tested scale.

![Actual body edges before and after, with the exact seam](figures/small-elbow-before-after.png)

Both panels sample edges from the actual B-reps, with identical axes and view.
The dashed seam is the independent analytic oracle. The before body cannot be
rendered through the normal diagnostics converter: it refuses a sheet touching
itself along an edge. Edge extraction therefore supplies the comparable view.
The corrected body also passes the normal renderer, with no faces refused.

## Cause and correction

Tracing the construction shows the swept pieces are valid before and after the
attempted trims, but the small case retains the extended ends. The first invalid
incidence appears after union. The trim's `find_side_of_curve` projects a reference
point onto a finite cutting curve; the nearest point can be an endpoint.

`CpGlobalSolver::local_solve_branch` already selects such endpoints when the
stationary-point search finds no root. It computes their projected distance,
but leaves `found_solution` false, discarding them. The trim then fails and the
default sweep path continues with untrimmed extensions.

The correction enumerates eligible finite boundaries as `SingleValue` candidates
independently of stationary-point searches. A stationary point can be a minimum
when a maximum was requested, so finding one must not suppress the boundaries.
Each endpoint is admitted only after its distance evaluates successfully. A
failed endpoint evaluation returns its error locally, rather than admitting an
error code as a numeric distance.
The outer legacy status policy is unchanged. Two further corrections are necessary
to preserve better interior candidates:

- The objective is `D(V)=|N cross V cross N|`. Writing `q=N dot N`, this is
  `q` times orthogonal projected distance. A sphere of radius `r` therefore maps
  inside a projected disk of radius `q*r`. Projected min/max pruning previously
  compared unscaled sphere bounds against a scaled objective incumbent. Both
  projected center distance and radii now scale by `q`, with the query tolerance
  added in objective units. Other projected operations retain their behavior.
- `ProjectedPointMaximize` negated its stationary function but not its derivative.
  For `C(t)=(t,0,0)`, `P=0`, `N=(0,0,s)`, that function is `F(t)=-s^2*t`:
  its derivative must be `-s^2`, not `+s^2`. Both are now negated together.

No new tolerance, geometry-specific case, retry, reference-point adjustment, or
compatibility flag is introduced.

The independent audit exposed the need for these corrections before readiness:
the initial endpoint-only candidate displaced an interior minimum in `gen2419`
(`t=0.19538434989712908`, distance `2.018003877141825`) with an endpoint at
distance `2.0600868346831587`, outside the query tolerance `0.001`.
For `gen0002`, its endpoint distance `6.513722074907321` missed the interior
maximum `9.248675280586827` at `t=0.8258440787616159` (tolerance `0.05`).
These are permanent regressions with independent rational-quadratic oracles:
Bernstein coefficient signs establish a unique stationary root, located by
bisection. The certificates use f64 arithmetic with large sign margins, not
formal outward-rounded interval arithmetic.

## Independent checks and references

The four elbow regressions require Level2 validity and the closed-form volume.
Direct solver regressions cover 1,024 combinations of analytic/B-spline lines,
scales, frames, orientation, full/subinterval domains, interior/exterior query
feet, min/max and Single/All requests. Twelve quadratic Bezier cases use the
exact curve `C(t)=(t,t^2,0)` and its differentiated squared-distance polynomial.
They check returned parameters, distances and curve evaluations. Parameters and
distances use `1e-9` relative comparison with a documented `1e-12` absolute floor;
the elbow volume check retains the existing `1e-4` relative oracle.

- Miller and Goldman, *Using Tangent Balls to Find Plane Sections of Natural
  Quadrics* (1992), pp. 72-75: the cylinder/plane section has semiaxes `r` and
  `r/abs(axis dot normal)`, hence `r` and `sqrt(2)*r` here.
  <https://people.eecs.ku.edu/~jrmiller/Papers/TangentBalls.pdf>
- Boyd and Vandenberghe, *Convex Optimization*, section 4.2.3, and the associated
  Stanford lecture: a constrained optimum at the boundary need not be a stationary
  point. This is the reason the existing endpoint selection must produce a candidate.
  <https://see.stanford.edu/materials/lsocoee364a/transcripts/ConvexOptimizationI-Lecture05.html>
- Boyd and Mattingley, *Branch and Bound Methods*, section 1: exclusion requires
  a bound on the same objective as the incumbent. The sphere-image derivation above
  supplies that bound for the existing double-cross metric.
  <https://web.stanford.edu/class/ee364b/lectures/bb_notes.pdf>
- MIT OpenCourseWare, *Newton's Method*: the update divides the function by its
  derivative. Changing the sign of only one reverses that update.
  <https://ocw.mit.edu/courses/18-s997-introduction-to-matlab-programming-fall-2011/pages/root-finding/newtons-method/>

The standalone pinned C++ reproduction confirms that both analytic and B-spline
lines return success with zero projected extrema for an exterior query point.
[WR-753](../upstream-issues/WR-753--projected-curve-endpoints-discarded.md)
records the exact inputs, expected endpoints, native outputs and build identity.
Full difftest results and final validation are recorded in the pull request;
Rust-only tests are not counted as evidence of C++ parity.

## Changed-case audit and limits

The revised candidate's [45 changed global-query cases](curve.curve_p1--projected-extrema-audit-2840.md)
have independently checked polynomial extrema. Every returned candidate is within
the original query tolerance on the cached domain. The out-of-bounds request in
`gen0975` is explicitly qualified: neither cache searches its extrapolated portion.
No claim of complete requested-domain correctness is made for that case.

[WR-754](../upstream-issues/WR-754--projected-extrema-pruning-metric.md) confirms
native pruning counterexamples; [WR-755](../upstream-issues/WR-755--projected-maximization-derivative-sign.md)
confirms the derivative sign. The [local convergence audit](curve.curve_p4--projected-derivative-local-convergence-2840.md)
names all changed evaluator and local-search cases, including two valid native
roots missed from the same starting guesses after the correction. Their global
queries still find the independently calculated maxima. Those local differences
are not classified as C++ defects.

The shared bounds use the objective of supported curve-point and curve-curve
projected searches. Separate topology pre-pruning and unsupported projected
surface-pair paths are not repaired here. The existing wrapper's status-suppression
policy also remains unchanged; propagating endpoint errors locally is not a claim
that every enclosing legacy caller propagates them.
