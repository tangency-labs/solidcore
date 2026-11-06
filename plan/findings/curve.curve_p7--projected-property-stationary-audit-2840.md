# Projected property-analysis changes in #2840

Historical candidate: `09495e082`; exact main base:
`3f82278a14c1a6104db8e090a59b4bcf677ace01`.
Reference: pinned v1.0.1, `v1-3fd5b9548d19-726643fb17a52b7d`.
Unit: `curve.curve_p7`. These 16 generated cases passed native parity on the base
and differ on the candidate. All request property 21, `ProjectedPointMaximize`.
This records that historical candidate, not subsequent affine-span changes.

## Meaning and method

`Curve::global_property_analysis` dispatches this property through
`CurvePropertyGs` and `local_property_analysis`. This is a separate path from
`Curve::global_point_solve`. Its local construction searches stationary equations;
the derivative-sign correction changes which stationary roots it finds. A stationary
root need not be a maximum. This audit distinguishes stationary-root correctness
from constrained global optimality; it does not infer a global-maximum guarantee
from the property name alone or classify every newly returned minimum as a defect.

The independent calculation constructs homogeneous polynomial coordinates on each
stored B-spline span, or the cubic Hermite polynomial, and uses the projection matrix
`(N dot N) I - N N^T`. With squared-distance numerator `A` and weight `W`, stationary
candidates are real roots of `A' W - 2 A W'`. Evaluating those roots and span boundaries
gives the stored-domain maximum. The second derivative distinguishes local minima
and maxima. This is numerical polynomial analysis in f64, not an outward-rounded
certificate. No native result supplies the oracle.

Native returns no stationary solutions in every case below. In five cases it does
return a visualization comb; the candidate preserves that comb exactly. All other
comb tails are also unchanged. The candidate adds 18 stationary records across the
16 cases. Sixteen positive reported distances agree with the independent metric
within `1e-9` relative without an absolute floor. The two zero reports have independent
residual distances approximately `3.51e-16` and `2.99e-15`; those cannot be described
as strict relative matches to zero. They remain explicitly qualified as near-zero
numerical evaluations, not tolerance-adjusted passes.

## Every changed case

Distances and parameters below are rounded for display. “Maximum” and “minimum”
classify the stationary point locally. The final column compares the returned
distances to the best distance on the requested interval intersected with the stored
domain, using the original query tolerance without an added allowance.

| Case | Query tolerance | Returned parameter(s) | Returned distance(s), local classification | Stored-domain maximum | Within query tolerance of maximum? |
|---|---:|---|---|---:|---|
| gen0204 | 0.1 | 0.9197458834128428 | 18.671389006901112, maximum | 33.337131619534404 | No |
| gen0247 | 0.1 | 5.320326948805733 | 3.4109831271925235, maximum | 3.4109831271925244 | Yes |
| gen0275 | 0.001 | 0.14871522759066944 | 29.774575607750659, maximum | 29.774575607750652 | Yes |
| gen0409 | 0.001 | 0.9575602192205206 | 4.6250196404095325, maximum | 4.625019640409532 | Yes |
| gen0509 | 0.001 | 0.9117265070172493 | 23.757505575290121, maximum | 23.757505575290114 | Yes, stored domain only |
| gen0735 | 0.1 | 0.6111736501042818 | 5.4676116781299413, minimum | 5.951094520215661 | No |
| gen0834 | 0.1 | 0.9198480289454729 | 5.0668265331752576, maximum | 5.066826533175257 | Yes |
| gen0933 | 1e-6 | 5.667762158756711 | 0.54203280602605997, maximum | 4.969521868984223 | No |
| gen0936 | 0.001 | 3.0095743060512032 | 5.3749025326382576, maximum | 5.3749025326382585 | Yes |
| gen1531 | 1e-6 | 0.6484668221329991 | 0.50654672270561063, maximum | 4.593165426545604 | No |
| gen1986 | 0.1 | 0.7060548169794392; 0.8337434617502735; 0.5623663475572159 | 0.086458069614975125, maximum; 0.012748673309883067, minimum; 0, minimum | 0.9576757302466963 | No |
| gen2049 | 0.001 | 1.6134778029167118 | 14.493501816701414, maximum | 14.493501816701421 | Yes |
| gen2341 | 1e-6 | 0.7483245754106941 | 13.228091125875489, minimum | 14.63657324597678 | No |
| gen2394 | 0.1 | 0.89496412675213 | 0, minimum | 3.1331587642464043 | No |
| gen2665 | 0.001 | 0.8132747916059548 | 6.621699366041069, maximum | 6.621699366041071 | Yes |
| gen2945 | 1e-6 | 0.1020561860909834 | 29.698809299538972, maximum | 29.698809299538972 | Yes |

`gen0509` and `gen2394` request Hermite parameters beyond its stored `[0,1]`
interval. The comparison above does not certify extrema on polynomial continuation
outside the stored interval. The lower-endpoint witness for `gen2394` is inside both
the requested and stored intervals, so its failure to attain the displayed maximum
does not depend on extrapolation.

Nine cases attain the stored-domain maximum within their query tolerance. Seven
do not; this is evidence about the returned candidate set, not a claim that a local
stationary-analysis contract was violated. No claim of complete root enumeration,
no blanket native-defect classification, and no `cpp_defect` parking follows from
this audit. A future change to the property-analysis contract must state whether it
returns stationary points, local extrema, or constrained global extrema and test that
contract directly.
