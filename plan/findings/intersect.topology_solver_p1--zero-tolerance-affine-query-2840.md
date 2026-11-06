# New refusal for an existing zero-tolerance affine query

Unit: `intersect.topology_solver_p1`; case `gen1646`, generated seed
`4917428403172592891`; issue/PR: #2840.
Candidate `f79b779f5` introduces this difference relative to `09495e082` and the
pinned native reference v1.0.1 (`v1-3fd5b9548d19-726643fb17a52b7d`).
`f290ed2de` has the same production code.

**Status: unresolved API-contract change.** The new affine certificate refuses
this ordinary finite query. The exact arithmetic below explains why a literal
zero optimization-gap requirement cannot be met, but does not establish that
the historical API's `dist_tol = 0` promised such a requirement. This is not a
confirmed native defect and must not receive `cpp_defect` parking on this evidence.

## Input and actual outputs

The local PointCurve solver makes the same swapped call `[0,-1]` three times.
It minimizes projected distance, requests All solutions, uses initial squared
incumbent `1e20`, and sets distance tolerance to **zero**. The queried curve is a
nonrational degree-one B-spline with these controls:

```text
(0.1569346301508998, -0.5329538213524829, -0.8439542626476992)
(0.4224662792522742, -0.4068072031595056, -0.7199688610446469)
(0.6052455395473925, -0.3885960084138218, -0.5099311278223793)
(0.7831613767269248, -0.01515768701313519, -0.701567737379753)
knots = [0,0,0.2628985868293973,0.8036901450210183,
         1.3482261889778533,1.3482261889778533]
query = [0,1.3482261889778533]
point = (0.5501272357521876,-0.3940877155300753,-0.5732694057512258)
normal = (0.7405015371827215,0.6291959064853571,0.23615669520912666)
```

The requested interval equals the stored domain; there is no extrapolation
qualification. The other curve and vertex in the generated scene do not participate
in these calls (`pv = -1`, explicit point supplied).

Native and `09495e082` both return status `0` and the same `values.out` array:

```text
[3.0, 1.8615654003497463e-33, 1.0,
 1.0, 1.0, 1.0, 0.0, -1.0,
 4.3145861914553826e-17, 0.6406108473689742]
```

This means three completed calls, final squared incumbent
`1.8615654003497463e-33`, and one SingleValue solution on object 0, at parameter
`0.6406108473689742`, reported distance `4.3145861914553826e-17`.

The affine candidate returns status `1021` (`SM_ERR_NOT_WITHIN_TOLERANCE`) and:

```text
[1.0, 1e20, 0.0]
```

The first call therefore fails, leaves the initial incumbent unchanged, and returns
no solutions. This is a newly observable refusal, not merely a changed parameter or
additional candidate.

## Exact-geometry audit

Treat each supplied binary64 coordinate, knot and normal component as an exact
rational. On each span form `g + s v`, with projection matrix
`M = (N dot N) I - N N^T`, and minimize its squared length by the clamped fraction
`s = -(g dot v)/(v dot v)`. The rational calculations are exact; displayed square
roots are approximations computed with 80 decimal digits.

| Span | Minimum distance | Parameter |
|---|---:|---:|
| First | 0.13842476800178444 | 0.2628985868293973 |
| Middle | 1.5495187277368766e-16 | approximately 0.6406108473689742 |
| Last | 0.00022310753859584456 | approximately 0.9043638303370678 |

The middle-span squared distance is strictly convex, with a unique interior
minimizer. Its exact parameter's reduced rational denominator has 252 bits and is
not a power of two. It is therefore not a representable binary64 parameter.

| Representable parameter | Exact-definition distance, approximately | Excess above the mathematical minimum |
|---:|---:|---:|
| 0.6406108473689741 | 1.5873948487168345e-16 | 3.787612097995788e-18 |
| 0.6406108473689742 | 1.5507659273475355e-16 | 1.247199610658913e-19 |
| 0.6406108473689743 | 1.6189580029951511e-16 | 6.943927525827446e-18 |

The old returned parameter is the nearest representable minimizer. Nevertheless,
every representable parameter has a strictly positive optimization gap under this
exact-definition interpretation. The old reported distance is smaller than its
exact-definition distance because evaluating and subtracting nearly coincident
points in binary64 changes the residual. No strict relative agreement with the
exact-definition metric is claimed, and no absolute floor is introduced.

## Contract decision still needed

Historically, `dist_tol` controls solution retention and pruning while local Newton
solves use separate nonzero convergence tolerances. Consequently, accepting zero
as this caller does cannot automatically be interpreted as requiring a provably
zero optimization gap in the exact input geometry.

The new certificate's refusal is consistent with a literal zero-gap contract,
but whether that is the intended public contract requires maintainer guidance.
The previous successful near-zero interior answer is not evidence of a native-only
bug. Conversely, its numerical residual is not a proof of exact minimization.
This record preserves both facts without changing tolerances, error classifications,
known-mismatch entries, or the caller's inputs.
