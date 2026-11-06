# Project-and-trim caller audit in #2840

Historical candidate: `09495e082`; exact main base:
`3f82278a14c1a6104db8e090a59b4bcf677ace01`.
Reference: pinned v1.0.1, `v1-3fd5b9548d19-726643fb17a52b7d`.
Unit: `primitive.api_trim_ops`. These five cases changed from native parity on the
base. This is evidence about that candidate, not subsequent affine-span changes.

| Case | Native result | Candidate result | Assessment |
|---|---|---|---|
| `pat_misses` | Error 1001 | 26 faces, 48 edges, 24 vertices | Keeping the entire filleted box is supported by the side proof below. |
| `pat_plain_misses` | Error 1001 | 6 faces, 12 edges, 8 vertices | Keeping the entire plain box is supported by the same proof. |
| `gen0070` | 1 face, 4 edges, 4 vertices | Empty | Unresolved; no independent proof supports deleting every face. |
| `gen0084` | 3 faces, 9 edges, 7 vertices | 5 faces, 16 edges, 12 vertices | Unresolved, with opposite-side retained vertices documented below. |
| `gen0086` | 2 faces, 7 edges, 6 vertices | 3 faces, 11 edges, 9 vertices | Unresolved; larger face count is not evidence of correctness. |

## Caller and independent proof for the misses

`primitive/api_trim_ops.rs::project_and_trim` normalizes the direction and calls
the projection pipeline with `TrimType::KeepPoint`. The side is that of the
reference point. `topology/brep/p10.rs::faces_to_delete` votes with face vertices,
then edge midpoints when no vertex votes exist, then a ray-based tie fallback.
The output can be an open collection of retained faces; these counts do not by
themselves describe a closed solid.

Both miss cases use a box from `(0,0,0)` of size `(10,5,15)`, direction `(0,0,-1)`,
reference `(5,1,20)`, and curve points
`(100,-2,20),(101,5,20),(102,0,20),(103,2,20)`. `pat_misses` additionally fillets
the box at radius `0.5`. The fixture is `fixtures/primitive.api_trim_ops/cases.json`
(entries 38 and 46 in the report). The runner's `curves_of` calls
`primitive/api_curves.rs::create_curve`: these four points are cubic Bezier
**controls**, not interpolation samples.

Consequently `x(t)=100+3t` and `y(t)=-2+21t-36t^2+19t^3` on `[0,1]`.
The control hull gives `y in [-2,5]`; derivative controls give `y' in [-15,21]`.
For every projected box point `q` with `qx in [0,10]`, `qy in [0,5]`, squared
projected distance satisfies

```
D'(t)/2 = 3*(x(t)-qx) + (y(t)-qy)*y'(t) >= 270 - 147 = 123 > 0.
```

The reference also satisfies these bounds. The constrained closest point is
therefore always the first endpoint, whose tangent is `(3,21,0)`. Its planar
signed-side expression is `3*(qy+2)-21*(qx-100)`, strictly positive throughout
the box and at the reference. Changing to projection direction `-Z` reverses
both signs together. The curve never intersects the box (`x>=100`), and all box
points belong to the reference side; the filleted box is a subset. This proves
the intended keep-side outcome independently of the native error.

## Generated cases that remain unresolved

`gen0070` (seed `12099798707737135304`) projects a finite arc centered at
`(5.54295238374327,-2.6847898088439193,9.007984145688914)`, radius
`2.911016437613864`, angles `5.717811128273231` to `147.61560969045235` degrees,
along `(-.013317243772994564,.02179333640208525,-1)`. The reference is
`(4.793100698546034,-1.2787795102224686,4.007984145688914)`. Finite arc endpoint
effects preclude treating it as a full circular-cylinder separator. Neither the
native single face nor the candidate empty result is proved correct by this audit.

`gen0084` (seed `12099798707737135318`) uses a segment from
`(.22719607435831102,2.8788985533302607,15.722419741707533)` to
`(6.991255193308396,10.760981248163384,15.722419741707533)`, direction
`(.09905487761716145,.04457169745233888,-1)`, and reference
`(3.6092256338333533,.7212599630351749,10.722419741707533)`.
With `a` the start and `n=(end-a) cross direction`, the plane normal is approximately
`(-7.88208269,6.76405912,-.47927314)`. The reference has `n dot (reference-a)`
approximately `-38.8555`. Three retained candidate vertices instead have signed
values `+50.6405,+56.4933,+29.8359`. All native retained vertices are reference-side
or on-plane. For a segment the side sign is independent of the chosen foot along
its supporting line, so endpoint clamping does not explain away these opposite
signs. However, the finite projected segment ends just inside the bottom/left box
boundaries. This is adverse evidence against calling the larger candidate result
correct, not a proof that the native result is a complete intended trim.

`gen0086` (seed `12099798707737135320`) projects a finite arc centered at
`(6.149438126529027,-2.7504818875388906,11.223905147812282)`, radius
`7.311995406273734`, angles `31.236694384057458` to `174.6285077384253` degrees,
along `(.09148384575384799,.08585647501698548,-2.2695534158749506)`. Its reference
is `(6.674698364224811,-.1814642870635238,6.223905147812281)`. As for `gen0070`,
finite endpoints mean a full-circle inside/outside oracle is insufficient.

For all three generated cases, the serialized primitive record supplies reachable
vertices, edge samples, and surface-domain samples, but not complete face-loop
incidence. A surface-domain sample need not lie inside the trimmed face. Resolution
requires per-face vote traces: query point, projected parameter and distance,
tangent, side, tie-fallback outcome, and deletion decision, checked against the
actual finite projection and retained face regions. No tolerance is relaxed and
none of these five cases is parked as a `cpp_defect` by this finding.
