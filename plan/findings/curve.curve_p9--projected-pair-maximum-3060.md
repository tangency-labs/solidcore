# Projected curve-pair maximum: incomplete stationary search

The reported cases return a successful projected maximum smaller than a feasible
pair. On the current parent stack, `gen1839` returns three distances no larger
than 4.243690092983992e-12; a corner pair has distance 4.40801559960347.

![The reported cylinder pair and a larger admissible pair, same view.](figures/projected-pair-witness-3060.png)

For the four cylinder/half-cylinder curve/edge-face fixtures, the projected
surface is [-0.5,0.5]^2 and the query segment is [-1,1] x {0}. The old pair
has distance 0.5. Squared distance is convex in each coordinate, so opposite
x endpoints with |y|=0.5 attain the exact maximum sqrt(2.5).

## Contract and scope

The returned maximum must have admissible parameters on both definitions and
satisfy `max |N cross (A(s)-B(t)) cross N| <= reported + tolerance`, without
changing the tolerance. Its scalar must also agree with the returned pair.
The double-cross map preserves the caller's non-unit normal scaling.

`SolutionRequested::Single` asks for the best solution; `All` asks for every
solution within tolerance. Of the six reported generated cases only `gen0844`
is Single; the other five are All. A Single construction must not silently
collapse All into one approximate point. The initial scope is Single plus
Single variants of those witnesses, not a claim that all of #3060 is closed.
The enumeration contract and completeness remain separate work.

## Design research

Neumaier, *Complete Search in Continuous Global Optimization and Constraint
Satisfaction*, sections 11 and 12, gives the mean-value enclosure (27) and the
branch-and-bound invariant: discard a parameter box only when its bound cannot
beat an established feasible candidate. Local stationary solves alone do not
establish global optimality. [Paper](https://arnold-neumaier.at/glopt/mss/Neu04.pdf).
The rounded interval framework is also described by
[Patrikalakis, Maekawa and Cho, section 4.7](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node45.html).

The proposed replacement for Single projected maxima partitions the full
parameter rectangle at both curves' knots. It bounds the squared double-cross
distance with existing outward-rounded curve derivatives and the mean-value
form, and subdivides unresolved boxes deterministically. Feasible evaluations
at box corners/centers provide lower bounds; they are candidates, never a
completeness certificate. Work and depth are bounded with cancellation;
an unresolved bound or an unsupported definition returns a typed refusal.

Use plain NURBS' stored homogeneous controls directly, without dehomogenizing
and rounding them back into a new net. Analytic lines reuse the existing affine
definition proof; Hermite curves use their stored power polynomial. Analytic
or procedural curves must not have an unverified cached approximation treated
as their definition. Supported exactness classes and refusals will be explicit.

This belongs entirely to the kernel. #2791's N-D evaluation-state correction
is separate; this work does not change it. The branch includes #2968 and #2840
so their existing fixes and enclosures are reused rather than duplicated.


## Implemented Single path (work in progress)

The Single solver now constructs feasible parameter pairs and encloses every
parameter rectangle with interval value and mean-value derivative bounds. Only a
proved upper bound within the original distance tolerance may prune a rectangle.
The operation keeps the kernel double-cross metric, including non-unit normals.
Incumbent filtering returns no answer only when the entire certified maximum is
below the cutoff; an ambiguous cutoff refuses. Cancellation and proof failures
propagate independently of the legacy error-passing flag.

Definition enclosures cover native positive-weight B-splines (degree at most 18,
so the Taylor helper's factorials remain exact integers), Hermite power
coefficients, analytic lines, and ellipses whose native rational form passes the
existing `curve_on_quadrics` exact-class certificate used by record construction.
An analytic label alone is not accepted: a changed radius with a stale net is
refused. The new helper reuses that existing policy rather than a sampled
recognizer. Other definition kinds are still unsupported.

Six targeted tests pass locally: constrained corners; trimmed projection at
normal scales 0, 0.001, 1, and 1000 with nonzero Z; same-object disjoint intervals
and incumbent cutoffs; cancellation and precision refusal under both legacy
flags; the rational gen1839 **Single variant**; and a valid ellipse versus a stale
analytic net. The original gen1839 All row remains a separate obligation. Full
case and downstream audits remain pending.
