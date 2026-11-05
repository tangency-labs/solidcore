#### WR-755 — Projected maximization negates the function but not its derivative

- **Upstream:** `source/SMLib/src/SmCurve.cpp:15532-15582`,
  `SmCurvePropertyEFO::Evaluate`, `SM_CP_PROJECTED_POINT_MAXIMIZE`, pinned
  v1.0.1 commit `3fd5b9548d19aa8c12c29bc9637a96d42a06f8a5`.
- **Trigger:** line `C(t)=(t,0,0)`, query point zero, projection normal
  `N=(0,0,s)`, evaluated at `t=0.25`. The stationary function is
  `F(t)=-s^2*t`, so its derivative is exactly `-s^2`.
- **Compiled C++ confirmation:** the evaluator is file-local. The difftest
  shim contains its source; the entire projected min/max branch was checked
  byte-for-byte against the pinned source and matches. It reports the correct
  negative function value but a positive derivative. This is a compiled source
  reproduction, not a claim to call a private exported symbol from the library.

  | Normal scale `s` | Expected derivative | C++ derivative | Rust derivative |
  |---:|---:|---:|---:|
  | `0.5` | `-0.25` | `0.25` | `-0.25` |
  | `1` | `-1` | `1` | `-1` |
  | `2` | `-4` | `4` | `-4` |

  Inputs: `fixtures/curve.curve_p4/issue_2840_projected_derivative.json`.
  Native geometric evaluations use the unmodified reference
  `v1-3fd5b9548d19-726643fb17a52b7d`.
- **Cause:** the maximization branch changes `rdFOfT` to `-rdFOfT` and
  leaves `rdFPrimeOfT` unchanged, reversing Newton's intended step direction.
- **How found:** independent extrema audit for the small-profile mitre fix
  (#2840). Correcting this derivative is necessary for the rational-curve
  interior-maximum regression.

Correct derivatives change Newton convergence basins. This record does not
classify every changed local-search result as a wrong C++ answer: two native
roots remain valid but are missed from the same starting guesses after the
correction. See [the local convergence audit](../findings/curve.curve_p4--projected-derivative-local-convergence-2840.md).
