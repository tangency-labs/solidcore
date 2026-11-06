# Projected curve-pair maximum changes in #2840

Historical candidate: `09495e082`; exact main base:
`3f82278a14c1a6104db8e090a59b4bcf677ace01`.
Reference: pinned v1.0.1, `v1-3fd5b9548d19-726643fb17a52b7d`.
Unit: `curve.curve_p9`. These seven cases passed native parity on the base and
differ on the candidate. Every case calls `Curve::global_curve_solve` with
`ProjectedMaximize`. This is evidence for that historical candidate, not subsequent
affine-span changes.

## Output semantics and independent evidence

Each single solution contains a distance and two curve parameters. The candidate's
changed endpoint reciprocal solves call `local_property_analysis` with
`ProjectedPointMaximize`; a newly converged stationary point is not necessarily a
maximum. The public global curve-pair maximum remains independently checkable by
constructing another admissible pair with a larger distance.

For straight segments, evaluate their affine definitions at the interval corners.
For B-splines, construct the homogeneous span polynomials from the published controls
and knots, then evaluate the selected corner pair. The projection metric is
`|(N dot N) V - N (N dot V)|`; these seven cases all use unit Z as the normal.
Corner witnesses suffice to disprove maximality; they need not be the exact global
maximum. No result from the native solver supplies these witnesses.

## Every changed case

| Case | Original tolerance | Native result | Candidate result | Independent conclusion |
|---|---:|---|---|---|
| circle_projected_max | 1e-5 | Two distances 4.5, parameters near `(0.5,0)` and `(0.5,1)` | Same distances, exact parameters `(0.5,0)` and `(0.5,1)` | Both attain the maximum: center separation 3 plus radii 1 and 0.5 gives the triangle bound 4.5. Parameter-rounding differences at the closed seam remain native mismatches, not incorrect geometry. |
| gen0472 | 1e-5 | Two distances 5.0105715611725925; one parameter near zero | Two distances 5.010571561172592; exact zero replaces the near-zero parameter | Neither is a global maximum. Admissible ellipse angles 0 and 90 degrees give distance 5.346613700814649. |
| gen0844 | 0.5 | Distance 0 at `(2.340415323454198,2.340415323454198)` | Distance 0 at `(1.63857830683894,1.63857830683894)` | Neither is a maximum. Corner `(-2.44896896704292,3.043140518489156)` gives 10.979584714079621. |
| gen1839 | 0.0001 | One coincident pair, distance 0 | Original pair plus two pairs with distances about 4.24e-12 and 2.51e-15 | Neither candidate set contains a maximum. Corner `(0.43561427937496744,0.03474215986040488)` gives 4.40801559960347. |
| gen2050 | 0.0001 | One pair, distance about 6.1371363e-8 | Original pair plus two endpoint pairs of the same small distance | Neither candidate set contains a maximum. Corner `(0.42689846142600174,0.30573599521598394)` gives 40.00123197279587. |
| gen2187 | 1e-6 | One pair, distance 0.09955759757712898 | Original pair plus two endpoint pairs of the same small distance | Neither candidate set contains a maximum. Corner `(3.254294670614625,-3.304387257364806)` gives 6.310837322899527. |
| gen2585 | 0.5 | One pair, distance about 0.000128676797868 | Original pair plus two endpoint pairs of the same small distance | Neither candidate set contains a maximum. Corner `(0.6795950883084139,1.0610357492278748)` gives 0.9962268759220879, exceeding the returned distance by more than the original tolerance. |

For `gen0472`, the ellipse radii are `4.723327346620286` and
`2.505285780586296`. Both witness angles lie in the source arc
`[-90,101.59684508222117]` degrees, and the requested NURBS interval contains the
whole stored `[0,1]` domain. The witness distance is `sqrt(rx^2 + ry^2)`.
Thus the result's failure is independent of out-of-domain extrapolation.

For `gen0844`, the first requested interval extends below the stored line domain.
The displayed witness deliberately uses the stored lower boundary and lies within
both the requested and stored domains. It already disproves maximality without
assuming extrapolated-search semantics. The two line definitions are identical,
so equal returned parameters have exactly zero distance.

The other four witnesses are inside the original requested and stored intervals.
Their improvements exceed the original tolerances by wide margins. These numerical
polynomial/affine evaluations establish counterexamples, not formal interval
certificates or exhaustive global maximization. Very small reported residuals are
not claimed to meet strict relative agreement with a differently rounded independent
evaluation; no absolute floor or relaxed tolerance is introduced.

## Status

The six unresolved maximum cases already return suboptimal results in the pinned
native reference; the candidate retains that problem, sometimes adding more
near-minimum pairs. These differences are not evidence that the global maximum API
was repaired and must not be parked as confirmed native-only defects. The circle
case is geometrically correct in both implementations. This record makes no new
`cpp_defect` classification and changes no known-mismatch entries.
