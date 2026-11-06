# Existing caller changes from the affine projected-extrema construction

Unit: `curve.curve_p1`; issue/PR: #2840. This audit compares affine candidate
`f79b779f5`'s four changed existing generated cases with historical candidate
`09495e082`. The exact main comparison base remains
`3f82278a14c1a6104db8e090a59b4bcf677ace01`; the pinned native reference is
v1.0.1, `v1-3fd5b9548d19-726643fb17a52b7d`.
`f290ed2de` has the same production code and changes only a test helper call.
The new regression fixtures are separate from these four existing callers.

All four request `ProjectedMinimize` with `SolutionRequested::All`. The unchanged
generated inputs, old/new output arrays, and closed-form geometry determine the
classifications below. No tolerance is loosened and no mismatch is parked by this
record.

## Method

On an affine span write the projected displacement as `g + s v`, using the API's
double-cross projection `M = (N dot N) I - N N^T`. Its constrained minimum is at
`s = clamp(-(g dot v)/(v dot v), 0, 1)`. The line cases use their normalized
direction and scale. The piecewise-linear B-spline uses its original controls and
knots, with the absent-weight sentinel interpreted as nonrational coordinates.

For the near-zero B-spline case, exact rational arithmetic on the input binary64
values computes the projection, foot fraction, and squared distances. Square roots
are then evaluated with 70 decimal digits. This distinguishes a rounded evaluator's
zero from a mathematically zero distance. The other spans' endpoint minima are much
larger, so the first span supplies the global minimum on the stored domain.

## Every changed existing case

| Case | Original tolerance | `09495e082` result `(parameter, distance)` | Affine candidate result | Interpretation |
|---|---:|---|---|---|
| gen0073 | 1e-6 | `(7.424378115297849, 2.4569516620450993)` | `(7.424378115297848, 2.4569516620450993)` | Same interior analytic minimum; one adjacent parameter representation changes. |
| gen1172 | 0.05 | `(0.2860207154511576, 6.678445460217803e-17)` | `(0.28602071545115765, 0)` | Near-zero evaluation cancellation at a neighboring parameter; a strict relative distance mismatch remains. |
| gen2450 | 1e-6 | `(2.554070753458807, 0.2700249248626431)` | `(2.5540707534588063, 0.27002492486264307)` | Same interior analytic minimum; parameter and distance rounding change. |
| gen2883 | 0.001 | Empty | Range `[0.9300363540079557, 1]`, both endpoint distances zero | Correct constant-zero range on the requested interval intersected with the stored domain; the requested extrapolated remainder is not covered. |

### gen0073 and gen2450: interior affine minima

The closed-form distances are approximately `2.4569516620450995` and
`0.2700249248626433`, respectively. Their foot parameters are approximately
`7.424378115297849` and `2.5540707534588067`, both strictly inside the original
requested intervals. Both old and new results are within the original query
tolerance of those minima. Old/new values and parameters also agree within `1e-9`
relative without an absolute floor. These changes are rounding-level differences,
not evidence of a previously incorrect native result.

### gen1172: zero is not an exact geometric certificate

The first span has parameter interval `[0, 1/3]` as represented in binary64.
Exact rational operations on the published controls, target, normal and knots give
minimum distance approximately `2.1765560538145463e-17`. At the old returned
parameter the exact-definition distance is approximately `7.137380185907995e-17`;
at the new parameter it is approximately `3.406054134638226e-17`. The other two
spans have minimum distances approximately `0.03560538710098969` and
`0.286485524436183`.

The new parameter is an excellent approximate optimizer under the caller's original
`0.05` tolerance. However, the reported zero does not equal the positive
exact-definition distance. Native/old Rust reports `6.678445460217803e-17`;
comparing that with zero fails the unchanged relative-only criterion. Neither an
absolute floor nor the optimization tolerance may turn that native mismatch into
a parity pass. This is a numerical evaluation/representative change, not a
confirmed native-only geometric defect, and it is not marked `cpp_defect`.

### gen2883: exactly constant geometry, qualified domain

Both degree-one controls are exactly the same binary64 point as the target. Thus
the projected displacement is identically zero for every parameter in the stored
domain `[0,1]`, for the supplied nonunit normal as well as any finite normal.
Returning the entire searched interval as `RangeOfValues` is justified directly;
the previous empty answer omitted this constant solution set.

The original query is `[0.9300363540079557, 1.119641960869513]` with out-of-bounds
evaluation disabled. The returned range stops at the stored endpoint `1`. This
audit certifies the range on the intersection only. It does not claim the entire
requested domain was searched, change extrapolation policy, or certify a global
repair of out-of-domain queries. Native's empty answer and the new range remain
different output structures; this record does not alter mismatch classifications.
