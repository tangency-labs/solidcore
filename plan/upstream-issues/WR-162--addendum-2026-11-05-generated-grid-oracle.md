Updates: WR-162

#### WR-162 (update, 2026-11-05). Generated request checked independently

- **Severity:** Wrong result.
- **Location:** `source/SMLib/src/SmOffsetSurface.cpp:1442`,
  `SmOffsetSurface::EvaluatePlainOrSimple`; the in-bounds `Duv` write aliases
  `Duu` when `hu=2`, `hv=0`, `bOnlyUpperHalf=TRUE`.
- **Confirmed:** yes, native `surface.offset_surface` case `gen0138`,
  seed `7714667220713742099`, pinned C++ v1.0.1 (`3fd5b95`). Both request
  shapes return success. At `(1.9765064514150505, 3.1211719029587046)`, the
  original request returns `Duu=(0.00324656755803084, 0.006764743619619604,
  0.005256450040563134)`; a separate same-prefix `(3,3,true)` request returns
  `(-0.007467434158383878, 0.0010213493368906668, 0.004523470578431292)`.
- **Independent check:** `harness/oracles/offset_grid.py` uses only the captured
  input, polynomial B-spline basis functions, and Decimal power-series arithmetic
  for the signed normal offset. At 140 digits its `Duu` is approximately
  `(-0.007467434158384999, 0.001021349336890052, 0.004523470578434926)`.
  The named larger-request C++ values and corrected Rust request agree within
  `8.04e-13` relative. Analytic plane and curved-extrusion checks pass; 80- and
  140-digit calculations agree within `3.45e-76` relative for this result.
- **Limits:** first request of this case only, at an interior regular parameter.
  The oracle does not use C++ outputs or process state to calculate its result;
  native values are compared afterwards. This does not confirm out-of-bounds
  UB-171, certify the other generated cases, park any failure, change a tolerance,
  or release #1914 / #2765 from draft. See the finding
  `surface.offset_surface--requested-derivative-grid-2717.md` for the unchanged
  diagnostic counts and the separate exact-zero case.
