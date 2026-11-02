# Setback scale search comparison (#2344a)

This contribution stacks on #2369 (`ee6e1b349`); #2274 is merged. The native target for
the scale change is SMLib `20b1c01`, not the pinned v1.0.1 reference. Only the two
Apache-2.0 files `source/SMLib/src/SmBlendSetback.cpp` and
`source/SMLib/inc/SmBlendSetback.h` were fetched from that commit.

The focused driver is `harness/cpp/tests/blend_setback_scale_search.cpp`. It can
be linked with those two files and existing fork SMLib/SM_API libraries (the
local comparison used the shared fork libraries at `5e8c9c1`). Put the two files
in a scratch include/source directory and compile the new source object before
linking the static libraries, so the old setback object is not selected.
Use the harness defines and `-O2 -ffp-contract=off -fno-fast-math`, under
`nice -n 10`. The driver checks statuses, exact scales, attempts, solve counts
and whether a refusal reason is present. It does not consume Rust outputs.

Native results:

| Corner | Scale | Attempts | QPs |
| --- | ---: | ---: | ---: |
| rounds 2,3,4, rule setbacks | 1.056 | 12 | 35 |
| rounds 1,1.5,2, rule setbacks | 1.056 | 12 | 35 |
| rounds 2,3,4, upright held at 4.5 | 1.046 | 12 | 39 |
| rounds 2,3,4, asked setbacks 4.25 | 1.110 | 12 | 34 |
| rounds 5,5,5, setbacks 8 | 1 | 1 | 2 |
| rounds 2,3,4, every edge held at its rule setback | 0 (refused) | 2 | 10 |

Rust focused gate:

```sh
CARGO_BUILD_JOBS=2 nice -n 10 cargo test --profile difftest -p solidcore --lib thousandth_scales_and_attempt_budget -- --nocapture
```

Explicit expensive monotonicity gate (separate from default tests):

```sh
CARGO_BUILD_JOBS=2 nice -n 10 cargo test --profile difftest -p solidcore --lib setback_monotone_grid -- --ignored --nocapture
```

The smallest certifying scale on the grid whose neighbour below fails; monotonicity checked, not proved.

The matrix in this contribution has seventeen descriptions, including axis
permutations. It is not the capi fixture inventory outside the source allow-list; the pocket rim
fixture remains a validation follow-up. The full twelve-case native executive gate now passes.
The pinned C++ difftests validate default behavior only. Do not treat them as
proof of opt-in fork parity. Streamed rows and their row/step counters are outside
this slice.

## Native executive reference

The executive driver is `harness/cpp/tests/blend_setback_executive.cpp`. It
constructs all inputs independently through SM_API and uses that API's context
for the executive too. A different context has different topology mark ownership
and is not an equivalent input.

The existing fork libraries are at `5e8c9c1`. For the executive, merely replacing
`SmBlendSetback.o` is insufficient: #104 changes the layout of blend corners.
Fetch the changed files under **source/SMLib/inc and source/SMLib/src only** from
`20b1c01` into an ignored overlay directory. The changed-source inventory is
the comparison of Git tree blob IDs at `5e8c9c1` and `20b1c01`; do not rely
on GitHub's compare-files endpoint, which truncates after 300 files. In addition
to the initial compare listing, the complete inventory includes
`SmFilletSection.cpp` and `SmProceduralBlend.cpp`.

Recompile the changed sources and all blend/fillet consumers of their headers,
then link the overlay before the existing static libraries:

```sh
bash harness/cpp/tests/build_setback_executive.sh \
  bin/reference-2344 /home/karlvonrichter/github/smlib-fork \
  /home/karlvonrichter/github/fork-cpp/lib bin/reference-2344/executive
nice -n 10 bash harness/cpp/tests/run_setback_executive.sh \
  bin/reference-2344/executive > /tmp/native-executive.log
CARGO_BUILD_JOBS=2 RUST_TEST_THREADS=2 nice -n 10 cargo test \
  --profile difftest -p solidcore --test fork_smlib104_setback_corner \
  dump_blends -- --ignored --nocapture > /tmp/rust-executive.log 2>&1
python3 harness/cpp/tests/compare_setback_executive.py \
  /tmp/rust-executive.log /tmp/native-executive.log
```

The comparator requires thirteen independently constructed executive results:
matching success/refusal, comments, named fits, topology counts, volumes and
every vertex coordinate, at 1e-9 relative with no absolute floor. It documents
only the existing unsupported-corner code-number mapping (native 10, Rust 9).
The Rust dump also requires Level2 `check_brep` to have no lines.

This reuses the shared library build without modifying it or creating another
checkout. Only the allowed SMLib sources are compiled. The native overlay is
not the pinned reference used for default-behavior difftests.

## Dependency refresh (2026-11-03)

Merged authoritative main `ebe6722e09169cce29a20c455fc6ec606b57caa0`
and #2369 `ee6e1b349e9ebae287246cc42a787a6620ab429e` into the
published branch, without rebasing or changing dependency branches. #2274 is
merged. #2369 is open, ready for review and has thirteen successful checks;
its review and merge remain dependencies of this stacked contribution.

The main merge resolved two documentation conflicts: retain `JoinSection = 16`
and `Precomputed = 14` with main's discriminant explanation; retain the new
Ball patch, End/Join and collapsed-offset documentation alongside the
thousandth-grid scale guarantee. The dependency merge kept that guarantee over
its old `to 1e-2` text and removed a blank line in the enum documentation.
Its added `BlendCorner.end` initializer is preserved.

The executive matrix now includes #2369's unequal (2, 2, 3) rounds on a
20 mm box at (20, 0, 20), as a thirteenth independent native/Rust input.
Its native scale-search volume is `7927.724282512234`, replacing the dependency's
old-search `7928.142563966`. Both states of the added Rust regression pass;
the flag-off expected volume was preserved.

### Fresh merged-code results

Workspace clippy (`--workspace --all-targets -- -D warnings`), wasm32 check
and formatting pass under Ubuntu WSL with two Cargo jobs and low priority.

The six-case scale/attempt/QP gate passes, as do all eight G2 regressions and
all seven direct setback tests in `fork_870c970_blend_setback`.
All six pinned C++ caller units pass on the merged code: executive p1–p3 and
corner p1–p3, 3172 cases including 45 guard rejections, zero mismatches or
deferrals (`difftest: 6 unit(s), 6 passed, 0 partial, 0 failed`).
The executive dump builds seven bodies and passes Rust Level2 checks.
Stored-versus-streamed native executive comparison passes all thirteen cases,
seven bodies and 167 vertices.

The **strict Rust/native executive comparison fails**: cases 9, 10 and 12 have
ten unmatched vertices with coordinate differences near zero, despite matching
volumes and topology counts. No absolute floor or tolerance was introduced.
Details: `plan/findings/blend.fillet_corner--2344-main-refresh-relative-zero-parity.md`.

The setback integration suite also has four flag-off volume-snapshot failures
(three pass, two helpers ignored). Their flag-on checks pass. No old expected
volume or tolerance was changed. Details and reproduction:
`plan/findings/blend.fillet_corner--2344-main-refresh-flag-off-volumes.md`.
These are explicit readiness blockers; historical proof below does not waive them.

### Historical evidence before this refresh

At `274e0a4a6`, the twelve-case native executive gate passed: six bodies,
144 vertices, statuses, comments, named fits, topology counts, volumes and
every coordinate at 1e-9 relative with no absolute floor. Workspace clippy,
wasm, fmt, the focused six-case scale/count gate, setback and G2 suites and
all six pinned C++ executive/corner caller units passed (3172 cases including
45 guard rejections; zero mismatches or deferrals). This is historical proof,
not evidence for the merge resolutions above.

Full `cargo test -p solidcore`, runner self-tests and `--green` were not run
for this slice; there is no unit promotion or CI-green claim. Production Rust
streaming remains unimplemented; the native preparation is documented separately.

### CI permission boundary

The last published-head run before this refresh,
[41450960214](https://github.com/tangency-labs/solidcore/actions/runs/41450960214)
at `274e0a4a6`, completed `action_required`. Earlier run
[41431300743](https://github.com/tangency-labs/solidcore/actions/runs/41431300743)
at `9f5dce3` had zero jobs and requires repository-admin approval. The earlier
approval attempt returned HTTP 403, "Must have admin rights to Repository".
The current account still has pull permission only on the authoritative repo.
No code edit or rerun can supply that permission; current-head workflow status
must be assessed separately after pushing.
