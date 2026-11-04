# blend.fillet_geom_p3: review evidence corrections (2026-11-04)

Updates: plan/findings/blend.fillet_geom_p3--issue-2349-fork-allocation-evidence.md

PR: https://github.com/tangency-labs/solidcore/pull/2381

The earlier source checkout and archive SHA256 records do not establish that the
archives were built from a18bbe523b9c7819f1a9da8c40a54c15f9936125. No producer
build stamp tying that revision to both existing archives was found beside the
build outputs. Prior fork-attributed allocation and sampled-geometry observations
are diagnostic results against those recorded archive bytes, not certified
revision-specific parity. This qualification also applies to the diagnosis and
exact-arcs addenda and the unchanged numerical sample fixtures.

The probe now requires --build-stamp before compiling or executing. Its JSON must
contain fork_commit and libraries (exactly libsmlib_api.a and libsmlib_core.a mapped
to SHA256). The revision must equal checked clean source HEAD; both hashes must
match the actual archives. The stamp must come from their build producer, generated
with the archives, not a retrospective assertion about unidentified binaries. A
stamp verifies consistency, not publisher authenticity; trust its producer as with
the reference build itself. The verified manifest is retained in new reports.
Saved reports without matching embedded provenance are refused. Historical reports
are preserved and are no longer accepted as certified inputs by the analyzer.

The analyzer also requires all ten fixture names exactly once before comparing any
record. Empty, missing, duplicate and unknown cases fail. No comparison field or
1e-9 relative tolerance changed. Probe functions now have explanatory docstrings.

Verification in Ubuntu WSL: seven Python tests pass, including analyzer integration
and archive revision/hash regressions. The first test commit deliberately failed
before the validation helpers were implemented. Existing native runner against the
configured reference: cases=10 passed=10 failed=0 deferred=0 max_rel_err=0e0.
No Rust production code or numerical sample changed, and no shared C++ library was
rebuilt. Full Rust checks are historical evidence, not rerun for this Python/docs fix.

Remaining dependency: obtain producer-stamped fork archives and rerun the native
probe before restoring revision-certified fork parity claims. Exact curved self-mitre,
exact two-sided pcurves and safe fit-size semantics remain separate unfinished #2349
work; this PR does not close that issue.
