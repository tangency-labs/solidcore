Updates: WR-162

#### WR-162 (update, 2026-11-05). Bicubic unequal request checked independently

- **Severity:** Wrong result.
- **Location:** `source/SMLib/src/SmOffsetSurface.cpp:1442`,
  `SmOffsetSurface::EvaluatePlainOrSimple`; unequal requested coordinates allow
  an unrequested mixed write to alias a requested slot.
- **Trigger:** `surface.offset_surface` case `gen0523`, generated seed
  `7714667220713742484`. A nonrational bicubic base, upper request `(2,0)`,
  `EvaluateSimple` at `(0.5326216929082825, 0.16048492851110552)`.
  The complete unchanged input is in
  `harness/oracles/offset_grid_tensor_cases/gen0523.json`.
- **Confirmed:** yes, fresh isolated native C++ v1.0.1 (`3fd5b95`) returns
  `Duu=(0.8251605633693019,-0.5884877400243681,5.989497909150047)` with success.
  The input-only tensor-series calculation at 140 digits gives approximately
  `(7.71081091620462,4.83384618851690,-0.211527114990600)`.
  A separate native `(3,3,true)` request returns the Rust values exactly.
  Corrected Rust and that larger-request C++ control agree with the independent
  result within `8.33e-14` relative. The base's normal length
  is approximately `4.269293140375243`, away from a pole.
- **Independent check:** `harness/oracles/offset_grid_tensor.py`; no engine
  output enters the calculation. Its offset-paraboloid point, mixed derivative
  and third derivative agree with separate closed forms to absolute `1e-125`.
  The case's 80-/140-digit requested values differ by at most `4.04e-78`
  absolute. This is numerical evidence, not an interval certificate.
- **Limits:** first request and three retained mismatching components only.
  No new mismatch marker or tolerance change. All 315 ordinary failures in
  #2765 remain, and neither draft is released. The separate finding
  `surface.offset_surface--requested-derivative-grid-2717--addendum-2026-11-05-tensor-oracle.md`
  records tiny bilinear derivatives for which relative accuracy of either
  native implementation is not established; do not generalize this
  confirmation to those cases or all 79 moved cases.
