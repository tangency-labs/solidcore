Updates: WR-42

The pinned v1.0.1 `source/SMLib/src/SmSurfaceCache.cpp:1711` still uses
`iv + 1 < rVConts.GetSize() - 1`. With three unique v knots, the sole
internal knot is skipped and a C0 crease can be reported as CInfinity.

Rust's default-off `Kernel::fixes.surface_cache_last_v_knot` enables the
corrected bound, `iv + 1 < v_conts.len()`, through the existing
`SurfaceCache::surface_continuity` query, which now takes `&Kernel`.
The flag is read on every query, including caches constructed before it is enabled.

Production callers are `Surface::create_silhouette_curves_with_cache`,
`Surface::create_parallel_projection_curves_in` (both context-backed and standalone
caches), and `surface_tracer::do_trace_base`. Their upstream counterparts are
`SmSurface::CreateSilhouetteCurves` (SmSurface.cpp:9621),
`SmSurface::CreateParallelProjectionCurves` (SmSurface.cpp:9488), and
`SmSurfaceTracer::DoTrace` (SmSurfaceTracer.cpp:1431).
These calls only signal a diagnostic and continue upstream; Rust continues too.
SmSurface.cpp:9857 also queries continuity in `CreateSectionCurves`, but its
conditional body is commented out; the Rust section implementation has no query.
SmSurfaceIntersector.cpp:6183/6187 are commented-out calls.

The regression `crates/solidcore/tests/wr42_surface_cache_last_v_knot.rs`
constructs a degree-1 surface with a crease at its sole internal v knot, 0.5,
and computes its cached continuities from the geometry. The query across the
knot returns C0 with the flag enabled, versus CInfinity by default; domains on
either side remain CInfinity. A second regression calls the production
`SurfaceTracer::do_trace` with an oversized cached continuity array: only the
corrected loop visits the missing final knot and returns the existing bounds
error. This sentinel verifies propagation through the tracer without inventing
new failure behavior for C0 geometry. No C++ reference or difftest tolerance changes.
