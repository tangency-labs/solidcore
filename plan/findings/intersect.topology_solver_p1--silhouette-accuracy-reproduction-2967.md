# Silhouette accuracy: first reproduction and representation experiment

Issue #2967; claim #2968. No production implementation changes yet.

The historical failing query belongs to #2840's projected-curve candidate, not
unmodified main. This distinction matters: the latest main's stored golden for
`analytic_face_sphere_band_5_15` is success with zero solutions. An empty success
does not certify the maximum or demonstrate that main returns the same bad scalar.

## Geometry and executable diagnostic

The sphere is centered at `(0, 0, -1)`, radius `0.5`, longitude `[0, 360]`
degrees and latitude `[-60, 60]` degrees. The point is `(0.2, 0.3, 2)` and the
projection normal is unit Z. The silhouette is the equator. Its exact projected
maximum is `0.5 + sqrt(0.2^2 + 0.3^2) = 0.860555127546399`.

`crates/solidcore/examples/silhouette_accuracy.rs` builds this surface and traces
the silhouette using the face caller's `1e-4` approximation tolerance and
20-degree angle tolerance. It queries both the fitted 3D curve and an existing
`CrvOnSurf` made from the returned UV trace and original surface. It also emits
4,001 diagnostic samples. These samples are not a proof of a global bound.

Executed against the previously validated #2840 candidate library:

| Representation | Projected maximum | Curve parameter |
|---|---:|---:|
| Fitted 3D curve | 0.8605419808705617 | 2.0651143077201275 |
| UV trace evaluated on the sphere | 0.8605551275463991 | 2.0647615968273123 |

The first reproduces the recorded failing scalar. The second agrees with the
analytic maximum to floating-point precision, without changing query tolerance.
This is a representation experiment, not a complete production repair.

The sampled 3D-curve maximum is `0.8605419689531334`; the sampled on-surface
maximum is `0.8605550981685833`. The largest sampled curve-to-surface distance
is `8.695878052481907e-5`. Sampling the fitted curve more densely would still
measure that representation, rather than the original surface.

![Projection of the sphere equator and query point beside the sampled silhouette-to-surface error divided by query tolerance.](figures/silhouette-accuracy-2967.png)

The left curves overlap at this scale. The right plot shows sampled positional
deviation; it is not the query's optimization error or a certified bound.

## Research and next construction checks

Patrikalakis, Maekawa and Cho explain [distance stationarity on the original
surface](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node136.html)
and [interval enclosures for geometric computation](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node45.html).
These support checking the defining geometry and bounding errors rather than
grading a fitted curve solely by samples. The repository's
`docs/definitions-and-certificates.md` requires the same distinction.

For the existing double-cross projection `M = |N|^2 I - N N^T`, its Euclidean
operator norm is `|N|^2`. Consequently a proved positional error epsilon can
contribute up to `|N|^2 * epsilon` to the projected distance error, by the reverse
triangle inequality. This is a derived budget constraint, not a measured bound
for the current tracer.

`CrvOnSurf` removes the fitted 3D curve's positional deviation by evaluating
`S(uv(t))`. It does not by itself prove that the UV trace is the exact silhouette,
that trim endpoints are correct, or that the optimizer found every required
candidate. Next checks must cover those errors, nonunit normals, other surfaces,
and the surrounding public topology query. No arbitrary tolerance multiplier,
native mismatch parking or golden blessing is proposed by this experiment.
