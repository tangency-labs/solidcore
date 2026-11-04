# surface.offset_surface: default third-derivative diagnostic attribution (2026-11-04)

- **Unit / PR:** `surface.offset_surface`, #1914. The approved four-term correction is unconditional; this record does not propose another geometry change.
- **Measured revision:** `8ec4825280f7b2f9a41fc4183d609f8a080146e9`, on base `c4276a13c92ff696d40f331f8830fbb8cc5d269a`.
- **Command:** `CARGO_BUILD_JOBS=2 CARGO_TARGET_DIR=bin/review1914-target SOLIDCORE_CPP_DIR=<main-checkout>/bin/cpp-v101 nice -n 10 cargo +1.98.0 run -p difftest --features cpp -- --unit surface.offset_surface` (with `RUSTUP_TOOLCHAIN` and `RUSTC_WRAPPER` unset).
- **Observation:** 3,182 cases, 2,663 passed, 519 failed, zero deferred; maximum relative error `1.9972685729541615`. This is a failed native diagnostic, not a parity pass. A separate full `curve.crv_on_surf` run reached the shell's combined 1,200-second limit without a summary.

## One individually supported C++ defect

The existing `fixtures/surface.offset_surface/edge_cases.json#0`, `third_derivative_outputs`, is named in the confirmed original WR-160 and WR-161 records. Its base is a regular cubic B-spline patch; its seven steps request rectangular derivatives of total order at most three or equal-order `(3,3)` upper-half derivatives.

The ten reported mismatches map through the unchanged runner's output serializer:

| `values.out` indices | Zero-based step | Derivative / component | Cause |
| --- | --- | --- | --- |
| 19 | 0, `eval` | `D(1,2).z` | WR-160 and/or WR-161 |
| 50, 51 | 2, `eval` | `D(0,3).y,z` | WR-161 |
| 76, 77 | 4, `eval` | `D(0,3).y,z` | WR-161 |
| 86 | 4, `eval` | `D(1,2).z` | WR-160 and/or WR-161 |
| 114 | 5, `eval_simple` | `D(1,2).z` | WR-160 and/or WR-161 |
| 126, 127 | 6, `eval_simple` | `D(0,3).y,z` | WR-161 |
| 136 | 6, `eval_simple` | `D(1,2).z` | WR-160 and/or WR-161 |

For example, `out[19]` is C++ `0.8656646974185226` versus Rust `-0.21398042233584189`; `out[50]` is C++ `-2.316474981309934` versus Rust `6.683525018690066`. The four corrected expressions follow differentiation of `Bu × Bv` and `B + d*N`, and a separate default-kernel polynomial closed-form regression fails before the fix. This supports a case-specific `cpp_defect = "WR-161"` marker explaining both WR-160 and WR-161. The marker parks this case only; it does not count it as a pass or establish correctness of every unrelated output.

## Remaining failures stay visible

`upper_half_slot_overwrite` (`edge_cases.json#1`) is **not parked**. Its `out[26]` is step 1's unequal-order `(3,1)` upper-half request: the serialized `D(2,0).z` slot is overwritten by a `Duvv` write. C++ `0.42964767031822376` versus Rust `-0.055925937843527895` therefore is not evidence that the destination `Duu` is correct. The unequal-order slot issue remains outside this PR. Its other four mismatches do not justify parking the entire case.

The other 517 failures, including pole/step-off, API-state and generated cases, have not been individually attributed. No claim that they are all C++ defects or all Rust regressions is made. Together with the unequal-order case, **518 observed failures remain unparked**. No tolerance, input, runner, guard, fixture or comparison is changed to remove them. The native unit remains blocked pending attribution; Rust's own canonical CI is a separate gate.

The existing passing cases and failing inputs remain in the runner. A follow-up must establish the actual operation, derivative slot and valid base dispatch for each proposed marker; changing a value is not sufficient evidence. No complete generated-case attribution or end-to-end curve-on-surface improvement has been demonstrated.

## Reference and evidence bounds

All 492 SMLib/SM_API source blobs in the preserved reference copies were checked against pin `3fd5b9548d19aa8c12c29bc9637a96d42a06f8a5`. Shared `cpp-v101` library hashes matched preserved historical build receipts, but this is not a reproducible-build proof. Libraries and shared sources were not rebuilt or modified.

The native results above belong to `8ec4825`, not later integration heads. Local full report, logs and source-hash receipts were preserved before integration. Focused final-head checks and canonical CI must be reported separately; this record does not turn a timeout or a historical run into a current-head pass.
