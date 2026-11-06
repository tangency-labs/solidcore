#### WR-753 — Projected curve extrema discard selected interval endpoints

- **Upstream:** `source/SMLib/src/SmCurve.cpp:11289-11309`,
  `SmCPGlobalSolver::LocalSolve`, in NVIDIA usd-brep v1.0.1,
  commit `3fd5b9548d19aa8c12c29bc9637a96d42a06f8a5`.
- **Trigger:** a finite line `C(t)=(0.001*t,0,0)`, `t` in `[0,1]`,
  queried from `P=(0.02,0.0004,0.009)` with projection normal `(0,0,1)`,
  tolerance `1e-10`, and `SM_SR_SINGLE`, for either
  `SM_SO_PROJECTED_MINIMIZE` or `SM_SO_PROJECTED_MAXIMIZE`.
- **Expected:** projected squared distance is
  `(0.001*t-0.02)^2+0.0004^2`, strictly decreasing throughout the interval.
  Its minimum is at `t=1`, distance `0.019004210059878839`; its maximum is
  at `t=0`, distance `0.02000399960007998`. Neither is a stationary point.
- **Native confirmation:** a standalone program linked to the unmodified,
  checksum-verified prebuilt reference
  `v1-3fd5b9548d19-726643fb17a52b7d` returns `SM_SUCCESS` (1000) and zero
  solutions for both operations, for both a `SmBSplineCurve` line segment
  and an analytic `SmLine`. No patched object or replacement implementation
  is linked. The permitted `SmCurve.cpp` source was hash-checked against the
  pinned revision before compiling the reproduction.
- **Rust cross-check on identical inputs:** the corrected library returns one
  solution for each operation and curve kind: minimum `t=1`, distance
  `0.01900421005987884`; maximum `t=0`, distance `0.02000399960007998`.
  These are the closed-form endpoints and distances, rather than the empty C++
  result. The standalone Rust reproduction links the same library build used by
  the regression tests.
- **Cause:** after the stationary-point search reports no root, the code
  explicitly selects each eligible interval endpoint and computes its projected
  distance, but does not set `bFoundSolution`. The subsequent found-solution
  guard discards the candidate.
- **How found:** the small-profile mitre follow-up from #2811, implemented in
  #2840. Its cutting-line side query needs a finite-curve boundary extremum;
  losing that answer prevents trimming and leaves extended ends in the union.
  See [the construction finding](../findings/curve.curve_p1--projected-endpoints-small-mitres.md).

This is a confirmed wrong C++ result, not evidence that an empty answer is a
valid baseline. The corrected Rust endpoint behavior is checked against analytic
distance functions and valid closed-body geometry.
