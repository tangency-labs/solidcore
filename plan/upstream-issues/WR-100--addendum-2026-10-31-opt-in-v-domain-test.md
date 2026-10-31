Updates: WR-100

- **Location:** `source/SMLib/src/SmSurfaceDropLib.cpp:4250`, `SmSurface::DropCurve`; the pinned C++ still sets `bDegenSurfV` from `sUVDomain.XLength()`.
- **Rust fix:** `Kernel::fixes.drop_curve_v_domain_degeneracy` (default `false`) selects `s_uv_domain.y_length()` for the V-degeneracy test in `Surface::drop_curve_base`. The U test, natural-domain fallback and existing failure path are unchanged. With the field off, the duplicated U test remains for pinned C++ parity.
- **Regression:** `crates/solidcore/tests/wr100_drop_curve_v_domain.rs` uses a borrowed bilinear B-spline with a zero-length natural V domain and nonzero U domain. Before wiring the field, the enabled case returned success; with the fix it returns an error and records `SM_ERR_DEGENERATE_SURFACE`, `DegenSurfaceFail`, U=false and V=true. A U-only degenerate domain still fails, with the V flag corrected only when enabled. The default V-only case retains success and U=false/V=false.
- **Confirmation:** source verified against pinned C++; regression executed in Rust. No new C++ reproduction is claimed for the degenerate-domain fixture.
