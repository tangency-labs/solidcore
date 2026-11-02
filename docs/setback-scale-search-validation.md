# Setback scale search comparison (#2344a)

This contribution stacks on #2369 (`6c87a7323`) and #2274. The native target for
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

The comparator requires twelve independently constructed executive results:
matching success/refusal, comments, named fits, topology counts, volumes and
every vertex coordinate, at 1e-9 relative with no absolute floor. It documents
only the existing unsupported-corner code-number mapping (native 10, Rust 9).
The Rust dump also requires Level2 `check_brep` to have no lines.

This reuses the shared library build without modifying it or creating another
checkout. Only the allowed SMLib sources are compiled. The native overlay is
not the pinned reference used for default-behavior difftests.

## Dependency update and CI approval

The branch merges #2369 at `7b59ad6f8` and main at `23ea17c1c`. Both merged
automatically, with no manual conflict resolutions. The new G2 regression suite
passes all eight cases; the setback suite passes seven (two print helpers ignored).
The corrected G2 report concerns SolidCraft's wrapper; no refusal workaround was
added here.

Run [41431300743](https://github.com/tangency-labs/solidcore/actions/runs/41431300743)
at `9f5dce3` completed `action_required` with zero jobs. A repository admin
must approve the fork workflow. The parent received HTTP 403 on its approval
attempt because the account has pull permission only. This is not a test failure
and this contribution does not rerun CI to work around that approval gate.

Native executive result: PASS, 12 cases, 6 bodies, 144 vertices. Every status,
comment, named fit, count, volume and coordinate agrees at 1e-9 relative with
no absolute floor. Changed volumes are native-supported regression values.
The focused six-case Rust scale/attempt/QP gate also passes after both merges.

Final merged-branch checks: workspace clippy (`--all-targets -- -D warnings`),
wasm and fmt pass. All six pinned C++ caller units pass (executive p1-p3 and
corner p1-p3): 3172 cases including 45 guard rejections, zero mismatches and
zero deferrals. The first workspace clippy run found two owned test lints;
both were corrected and the complete workspace check was rerun successfully.
Full `cargo test -p solidcore`, runner self-tests and `--green` were not run;
there is no unit promotion or CI-green claim.
