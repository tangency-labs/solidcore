Updates: WR-449

## Independent boundary witness

The partial body from `sweep4-sphere_diff_3cylinders-all-r0.1` provides a
geometric witness independent of mass integration and closest-point solvers.
This is a related fillet-output observation, not a rerun of the original
unfilleted corpus measurement in WR-449.

After the native mass query, face 7 / loop 0 / edgeuse 2 holds a cubic UV trim
for edge 30. Evaluating that actual trim on its installed cone at parameter
`0.62499997224447945` returns the point
`(20.46403968821491, -20.397755691950664, -8.101529611135945)`.
Its distance to **every point** of the original stored rational edge is
strictly greater than **0.06**, exceeding the edge tolerance
`0.017345204306469773`.

The native operation returns status 1007 and leaves 25 vertices, 39 edges and
8 faces; its subsequent mass query succeeds. This is evidence about the returned
partial body, not a claim of successful filleting. The reference is upstream
`c9c5979`, x86-64 Linux release configuration with contraction disabled.

## Reproduce the bound

Run `python3 scripts/verify_wr449_boundary.py`; no native build is needed.
The captured binary64 inputs and observation sequence are in
`scripts/evidence/wr449-boundary.json`.

The degree-2 edge has five positive-weight homogeneous control points and knot
multiplicities 3/2/3. Thus its entire domain consists of two rational Bezier
spans, sharing the middle control point. Exact homogeneous de Casteljau
subdivision produces 18 terminal boxes, at maximum depth 12. Each rational
segment lies inside the bounding box of its projected controls; every box's
exact squared distance from the recorded point exceeds `(3/50)^2`.
Replacing each segment with both children preserves coverage of the entire edge.
All calculations after loading the recorded doubles use rational arithmetic.

The control points were obtained directly with
`GetControlPoint(SM_CP_HOMOGENEOUS_RATIONAL)`, which returns weighted XYZ and W;
the proof does not reconstruct them by multiplying rounded Euclidean controls.
An independent review checked this getter convention, recomputed all terminal
bounds and verified full-domain coverage.

## Scope

The point is an observed floating-point surface evaluation. The proof bounds
its distance to the mathematical rational curve defined by the stored double
controls; it does not certify ideal surface evaluation or a global Hausdorff
bound for the whole boundary.

#1215 restores native trim and mass parity, including this witness. Its parity
fix remains valid. This addendum changes no production behavior, tolerances,
fixtures, mismatch classifications or gate accounting. A physical correction
still needs its own diagnosis and opt-in validation.
