# Metric and domain audit for #2967

This addendum audits the outputs at `b9691212d` (production code `cac31e200`)
against the saved validated parent and lift-only snapshots. It does not claim
that the subsequent main merge at `f0ca192ef` has finished validation.

## Completed validation

The full `cargo test --profile goldens -p solidcore` run at `b9691212d`
completed with **5,182 passed, zero failed, 45 ignored**. Formatting and wasm
passed. Strict workspace/all-target Clippy found a needless owned `Surface`
argument in the test helper. Borrowing it in `e95977483` passed the Clippy
rerun; `5244ea9d7` applies the corresponding rustfmt layout.

The independently built `cac31e200` affine integration suite also passed all
eleven tests, confirming the earlier parent CI result.

## Changed-record reconstruction

The committed `diagnostics/silhouette_measurement_audit_2967.py` reconstructs
the runner's analytic surfaces and evaluates them at the returned UVs. It
now reports natural-domain and line-interval membership without padding, and
the transverse residual for directed queries. The static runner projects along
Z and directs along X, so the latter is `abs(a.y - b.y)`.

| Check, across 45 changed p2 records | Validated parent | Surface-lifted metric candidate |
|---|---:|---:|
| Magnitude discrepancy above query tolerance `1e-5` | 16 | 0 |
| Maximum magnitude discrepancy | `3.866563226329056e-5` | `0` in executed evaluations |
| Returned UV outside the natural rectangular domain | 0 | 0 |
| Returned line parameter outside `[0,1]` | 0 | 0 |
| Directed transverse residual above `1e-5` | 1 | 0 |
| Largest directed transverse residual | `2.9795892157180592e-5` | `2.8993241202173158e-13` |

The old directed violation is `analytic_face_sphere_band_5_11`; the largest
remaining directed residual is `analytic_face_torus_half_x_5_11`. Both sides
retain one returned solution per changed record, with unchanged solution
types, variable counts and object fields. No case was added or removed.

The magnitude check uses the absolute reported scalar; it is not a signed-value
certificate. The domain check covers these fixtures' natural rectangular
domains, not arbitrary trim loops. Neither check proves stationarity,
candidate completeness, or a global optimum. No goldens were blessed.

## Separating tracing changes from main

Compared with the lift-only snapshot, the metric tracing request changes 33
analytic p2 records. Comparing the metric snapshot with `b9691212d` changes
**zero p2 records**. Thus the main merge and cache-reader/affine review changes
did not introduce further p2 differences in this comparison.

The eleven p1 changes between metric and refreshed snapshots are ray queries:
`adv_rayfire_ahead_single_True`, `adv_rayfire_behind_start_single_False`,
`adv_rayfire_behind_start_single_True`, `gen1018`, `gen1037`, `gen1576`,
`gen1636`, `gen2119`, `gen2323`, `point_vertex_rayfire_behind_start`, and
`point_vertex_rayfire_behind_start_swapped`. The merge brought main's WR-62
vertex-ray correction; these are separate from silhouette measurement.

The refreshed golden comparison still reports p1: 3,592 cases, 24 changed;
p2: 2,575 cases, 65 changed; both have zero added, removed or broken cases.
That comparison includes the parent projected-query changes, and the
reconstruction above is not a justification for blessing all those outputs.

## Kernel-rig baseline failure

The disk/blade failure in the parent CI was also present on main. Merged
[#3040](https://github.com/tangency-labs/solidcore/pull/3040) reverts #2859 and
has a passing kernel-rig job. Parent `a1382547b` and child `f0ca192ef` now
include that revert through main `ea972bedf`. Their refreshed checks must still
finish; a passing revert PR is not substituted for their own CI result.
