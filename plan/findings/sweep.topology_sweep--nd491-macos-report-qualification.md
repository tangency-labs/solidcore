# sweep.topology_sweep: ND-491 report qualification on macOS (2026-10-31)

Unit: `sweep.topology_sweep`. PR: #1423 (reporting only).
Code tested: `716d59eb6e1fc24e998660d17c237b81fa4567b4`, based on
`d1c6fc48a8197e1c0c9db5f1a03e82ad96de0544`.
Reference: clean SMLib 1.0.1, `3fd5b9548d19aa8c12c29bc9637a96d42a06f8a5`;
macOS, Rust 1.98.0, optimized `difftest` profile.

## Unfiltered observation

The seeds-16 run against the validated exact-base artifact did **not** pass locally:

```text
FAIL sweep.topology_sweep cases=703 passed=699 failed=2 deferred=0
with_values=699 no_values=0 vacuous=0 rejected=2
```

The two disagreements were:

| Case / whole-sweep input | C++ | Rust |
| --- | --- | --- |
| `whole_stitch_linear_x_flat_box`, `[3,7,2,1,5,0,0,0,1,0,1]` | `values.points.len = 0`, `ints.whole.len = 1` | `216`, `8` |
| `whole_merge_rotation_y_90`, `[10,10,10,2,0,1,0,90,0,1,1]` | `ints.whole[1..5] = [3,3,20,20]` | `[2,2,19,19]` |

Both case IDs already have `cpp_defect = "ND-491"` entries in
`crates/difftest/known_mismatches/sweep.topology_sweep.toml`, currently restricted to
glibc. No entry, platform restriction, case, tolerance or unit status was changed.
The failed report was preserved; this run is not counted as a pass.

## Native-self replay and reporting evidence

Three distinct cases were each run once with the fixed pad sequence
`1,2,3,5,8,13,21,34,55,89`, retaining every result rather than selecting a passing
attempt. The heap diagnostic reruns C++ against its first native result; it does not
feed a native result into Rust.

| Case | Ordinary parity | Native-self observations | Heap-qualified? |
| --- | --- | --- | --- |
| `whole_stitch_rotation_x_30` | failed: C++ region count 2, Rust 1 | pad 1 changes the C++ count from 2 to 1; `observed=1`, `padded=1`, divergence recorded | no |
| `whole_stitch_linear_x_flat_box` | passed, 216 values compared | pad 1 changes C++ `values.points.len` from 216 to 0; `observed=1`, `padded=1`, divergence recorded | no |
| `whole_merge_rotation_y_90` | failed: C++ region count 3, Rust 2 | all ten requested pads observed and actually padded, no divergence in this sample | yes, but parity still failed |

The flat-box report provides a real lucky-pass counterexample. On the **same**
report and real parked entry, the trusted main-branch report generator advertised:

```text
[('sweep.topology_sweep', 'whole_stitch_linear_x_flat_box', 'ND-491')]
```

The PR generator instead returned no adoptable row and listed that tuple as
unqualified. The ten-pad-stable rotation report also produced no adoptable row,
because heap qualification is not a substitute for a passing Rust/C++ comparison.
Complete qualified-pass eligibility is covered separately by the Python regression
tests; no current-head native qualified-pass result is claimed here.

Reproduce a diagnostic, replacing `<case>` with one of the three names above:

```sh
CARGO_BUILD_JOBS=2 SOLIDCORE_CPP_DIR=<pinned-sdk> nice -n 10 cargo +1.98.0 run \
  --profile difftest -p difftest --features cpp -- \
  --unit sweep.topology_sweep --seeds 16 --filter <case> \
  --heap-pads 1,2,3,5,8,13,21,34,55,89
```

Addresses and outcomes can change between processes; the tabulated outcomes are
observations, not assertions that every future process must produce those variants.

## Interpretation and boundary

The C++-versus-C++ changes directly reproduce heap dependence on this macOS host
for the flat box and x-axis rotation. They are consistent with the existing
ND-491 pointer-keyed silhouette-order diagnosis in
`plan/upstream-issues/ND-491--whole-brep-topology-sweep-stitch-and-merge-follow-heap-layout.md`.
The y-axis rotation did not vary under the ten sampled pads; its local disagreement
is retained, not diagnosed away or repaired by this reporting PR.

No claim is made that the whole sweep unit is locally green, that ND-491 is fixed,
or that a finite pad sample proves independence from all heap layouts. Extending
the existing platform qualification or changing the kernel is separate work. The
canonical Linux CI gate remains required; its pre-existing parked-case policy is
not changed by #1423.
