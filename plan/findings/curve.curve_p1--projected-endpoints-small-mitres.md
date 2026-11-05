# Projected endpoint extrema and small-profile mitres

Related: #2840; the available geometry follow-up in #2811, under #2794.
Units: `curve.curve_p1`, caller `primitive.primitive_creation_p2`.

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

The correction marks an existing endpoint candidate as a `SingleValue` solution
after its distance evaluates successfully. A failed endpoint evaluation propagates
its error, rather than admitting an error code as a numeric distance. The existing
global comparison selects the optimum. No new tolerance, geometry-specific case,
retry, reference-point adjustment, or compatibility flag is introduced.

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

Native comparison results and final validation are recorded in the pull request.
This record does not claim C++ parity from Rust-only tests.
