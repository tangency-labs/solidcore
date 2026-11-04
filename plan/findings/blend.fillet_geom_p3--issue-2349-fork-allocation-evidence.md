# blend.fillet_geom_p3: #2349 fork allocation evidence (2026-11-04)

Updates: `plan/findings/blend.fillet_corner_p3--issue-2349-reachable-teardrop-angle.md`.

PR: https://github.com/tangency-labs/solidcore/pull/2381

## Resolved verification gaps

On source `e60349c7a6cbc86961f5ac622a85c6d8f6da9c9f`, Ubuntu WSL had
load 0.02 and approximately 31 GB available memory. With `CARGO_BUILD_JOBS=2`,
`RUST_TEST_THREADS=2`, `CARGO_INCREMENTAL=0`, and nice 10:

```text
cargo test --profile difftest -p solidcore --lib topology_split
test result: ok. 5 passed; 0 failed; 0 ignored; 0 measured; 3120 filtered out
```

This includes the allocation assertion in
`topology_split_stops_at_the_new_rail_edgeuse`, plus both caller-path tests.
Disabling incremental compilation avoids the query-cache write that exhausted
memory on the overloaded machine; it changes no numerical options or assertions.
The earlier failed attempts remain historical evidence, not current blockers.

GitHub CI run 41614989092 on that exact source has now passed fmt, workspace
Clippy, tests, wasm, deny, CI scripts and publish dry-run. Its native sweep was
not requested (skipped); the earlier local native results remain separately
identified. Workflow approval is no longer pending. Required independent review
still needs a maintainer; the PR author cannot supply it.

## Fork executable comparison of the allocation

Reference: the clean existing SMLib checkout at
`a18bbe523b9c7819f1a9da8c40a54c15f9936125`. Reused its existing
`libsmlib_api.a` and `libsmlib_core.a`; no kernel rebuild, shared checkout edit,
or change to `scripts/upstream-pin.env` was made. The probe records the archive
hashes and compiler command in its JSON output.

`harness/cpp/probes/issue_2349_fork_reference.py` compiles the existing real
executive/checkpoint shim against those fork headers and libraries. The fork's
smaller transitive includes require `SmAttribute.h` explicitly. The unrelated
two-offset chamfer export is omitted because it instantiates the old
`SmLinearCrossSectionFSG`, which the fork deleted. No closed-round checkpoint
body, observation, tolerance or input is changed.

The Rust side is run independently with `--digest-dir --digest-debug`, without
calling C++. Each native case runs in a fresh process. At both `solver_geometry`
and `corner_edges`, the complete prefix topology and all corner-vertex states
match Rust exactly. Both newly split vertices are absent from the owning list
(label -2), have tolerance exactly 0, and their complete five-value records
(point, tolerance, original parameter) match Rust bit for bit. This closes the
fork **allocation** executable-evidence gap; it is not a claim of completed
construction parity.

Reproduce from this worktree in Ubuntu WSL (the two fork paths name existing
shared artifacts on this machine):

```sh
export CARGO_BUILD_JOBS=2
nice -n 10 bin/difftest/difftest --unit blend.fillet_executive_p3 \
  --filter issue-2349 --digest-dir reports/difftest/issue-2349-fork-rust --digest-debug
nice -n 10 python3 harness/cpp/probes/issue_2349_fork_reference.py \
  --fork-source /home/karlvonrichter/github/sc-ec/target/smlib-a18bbe523b9c \
  --fork-lib-dir /home/karlvonrichter/github/sc-ec/target/release/build/solidcraft-smlib-b354ed1a891de0df/out \
  --rust-debug reports/difftest/issue-2349-fork-rust/blend.fillet_executive_p3.debug \
  --out-dir reports/difftest/issue-2349-fork-reference
```

## Broader differences remain visible

The probe retains all ten complete native records. Its optional `--pinned-shim`
diagnostic compares those records with the unchanged 1.0.1 executive shim. Six
records are exactly equal; four differ: radius-0.5 solver geometry, corner edges
and completed round, and the radius-3 completed/refused result. For example:

| Case / observation | Configured 1.0.1 reference | Fork a18bbe |
| --- | --- | --- |
| solver geometry, raw value 367 (generated surface sample) | -3.96873681983367 | -3.9687558729642345 |
| full r3, raw value 201 (mass-property record) | 457.03059824978146 | 457.03060521148683 |

These differences exceed the ordinary relative parity threshold; no threshold
is widened and the optional complete-record comparison exits unsuccessfully.
The complete radius-0.5 result also differs in topology ordering. The split
allocation state itself agrees at both affected checkpoints. The suspected
broader cause is the fork's redesigned geometry versus the default in-progress
port, rather than split allocation; isolating each later difference is follow-up
work, not a verified diagnosis here. In particular, these results cannot justify
promoting or removing construction flags, claiming whole-unit fork parity, or
certifying the unfinished self-mitre. The existing ten-case configured difftest
is retained without field removals or tolerance changes.

For that optional diagnostic, link the existing pinned shim objects without
rebuilding the reference, then add `--pinned-shim bin/issue-2349-pinned-shim.so`:

```sh
pinned_out=bin/difftest/build/solidcore-ffi-4109f6c1bb9c3982/out
nice -n 10 g++ -shared "$pinned_out"/*-blend_fillet_executive_p3.o \
  "$pinned_out"/*-sc_shim.o \
  /home/karlvonrichter/github/solidcore/bin/cpp/lib/libSM_API.a \
  /home/karlvonrichter/github/solidcore/bin/cpp/lib/libSMLib.a \
  -lpthread -Wl,--no-undefined -o bin/issue-2349-pinned-shim.so
```

The exact curved-spine self-mitre, exact pcurves and proven fit-size boundary
remain #2349 follow-ups. This PR does not implement or reserve those scopes.
