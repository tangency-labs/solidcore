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
fixture and full fork executive result comparison remain validation follow-ups.
The pinned C++ difftests validate default behavior only. Do not treat them as
proof of opt-in fork parity. Streamed rows and their row/step counters are outside
this slice.