# General: macOS reference built with GNU g++ and libstdc++, sin/cos fusion made irrelevant (2026-10-25)

PR #1027. Root-cause follow-up to #450 (libc++ `istream` extraction), #749 and #829.

- **Before (main d4e1814c3, Apple clang 21 + libc++, macOS 26 arm64):** `--green` gave
  `261 unit(s), 239 passed, 13 partial, 4 failed, 5 skipped(platform)`. All 2066
  `topology.database_io` cases were `skipped(platform)` (failing its value-share gate), as
  were the `topology.brep_data` db cases (and the file cases with an underflowing ASCII
  double) and the `surface.srf_in_volume` read_db cases. With the skip removed, 61 of the
  2066 database_io cases disagreed (9 aborted on the `>> char[256]` stack smash, UB-5) and
  2005 agreed bit for bit.
- **Change:** `harness/cpp/build.sh` builds the reference with Homebrew's GNU `g++-<N>` on
  macOS and writes `CXX`/`CXX_STDLIB`/`LIB_DIRS` to `ffi.env` (macOS only; Linux's
  `ffi.env` is unchanged). `crates/solidcore-ffi/build.rs` compiles the shims with that
  compiler and links that C++ runtime only (the `cc` crate otherwise adds libc++ as well).
- **sin/cos fusion.** The first GCC run had 11 new failing units (curve.conic,
  cone_p2/p3, volume_p2, offset_executive_p1, fillet_standard_solver_p1, api_primitives,
  surf_of_revolution_p3, surface_p2 `ssi_tor_sph_tan`, ...). Cause: GCC fuses a
  `sin(x)`/`cos(x)` pair into one call, and Darwin has no `sincos`, so it calls
  `cexp(0 + x*i)` (32 sites in libSMLib.a); Apple's `cexp` differs from its `sin`/`cos`.
  Routing `cexp` to Apple's `__sincos_stret` (what LLVM, and so the Rust side, calls for
  its fused pairs) fixed 8 of the 11, but GCC and LLVM fuse different pairs and
  `__sincos_stret` itself differs from `sin`/`cos` for some arguments. The shim
  (`harness/cpp/shim/sc_shim.cpp`, GNU g++ on macOS only) therefore defines both `cexp`
  and `__sincos_stret` in terms of separate `sin` and `cos` calls. The definitions take
  precedence over libSystem's in the difftest binary, so every fused pair on both sides
  returns what separate calls return and fusion sites no longer matter. This mirrors
  Linux, where glibc `sincos` equals `sin`/`cos` (apart from the rare 1-ulp cases of #357).
  Upstream SMLib calls neither function itself.
- **Platform hooks.** `UnitDiff::skipped_platform` was consulted only when the C++ was not
  libstdc++, but `curve.curve_p6` and `nurbs.b_spline_curve_n_lib` use it for the macOS
  allocator use-after-free in `DeformToEndPointTargets` (UB-251, #449/#793), which does not
  depend on the standard library. The runner now consults the hook for every case; the
  three istream units check `cpp_stdlib_is_libstdcxx()` themselves (`topology.brep_data`
  keeps its token classifier as `libcxx_skip_reason`, still unit-tested on Linux).
- **After (same machine, GNU g++ 16.2, libstdc++):**
  `261 unit(s), 245 passed, 9 partial, 2 failed, 5 skipped(platform)`. database_io
  (2066/2066), brep_data (2379/2379) and srf_in_volume (2112/2112) compare every case;
  `geom.local_solve1d` (#654) and `curve.conic` `gen0316` pass. The only
  `skipped(platform)` cases left are the 10 UB-251 deform cases; the 5 unit-level skips
  are the arm64 float-to-integer casts (#448). The 2 failures: `blend.fillet_solver_p5`
  (Linux-only rollover worker, #652/#653) and `topology.brep_p10`
  `read_Solids_fillet_ts` (glibc-only known mismatch, PR #1010).
- **Without GNU GCC** the build falls back to Apple clang and libc++ (with a warning), and
  the previous libc++ skips apply as before.
