# PR #2840: broad validation on rebased main

Candidate `10869974bc520909cb2d03bce38d6119900bde56` is based on main
`513e6fcfcd049d41e9487599880479de4add63f2`. This addendum supersedes the
incomplete-run statements in the earlier rebase validation record. It does not
resolve the documented solver contracts or certify the unresolved callers.

## Rust validation

The complete `cargo test -p solidcore --profile difftest --no-fail-fast` run
finished with exit status 0: **5,104 passed, zero failed, 45 ignored** across
442 result summaries. This includes 3,569 library tests and 1,535 integration
tests; doc tests contain no cases. The ignored diagnostic export was separately
run as described below. The full run used the optimized difftest profile and
two test threads; an earlier debug build was superseded, not counted as a pass.

Formatting, workspace/all-target clippy with warnings denied, wasm and the
20 focused regressions passed on the rebased code as recorded in
`curve.curve_p1--rebase-review-validation-2840.md`. Subsequent changes add only
the executed diagnostic source and evidence records, not production code.

## Full native sweep

The sweep completed with exit status 1: **280 units, 591,960 cases; 244 green,
22 partial, 14 failed**. Case outcomes are **589,585 passed, 196 failed,
1,402 known mismatches, 146 already parked and 631 rejected**.

The failures consist of this candidate's 187 already-audited cases and the nine
failures on exact main (`heal.gap`: seven; `primitive.primitive_creation_p2`:
two). The eleven additional cases compared with main are this PR's previously
recorded reproduction fixtures. No additional failing unit appeared.

Case records outside the twelve affected units match the completed Menger
instrumentation sweep after normalizing checkout prefixes. That sweep was
independently compared with main's CI artifact from run `42162728256`, artifact
`13477299760`: every outcome and mismatch payload matched; one passing
`sweep.topology_sweep::nm_test_rotation_y_30_stitch` case had a different maximum
relative error (`3.62089958464889e-14` here versus `5.940235728088039e-16` on
main), both below the unchanged `1e-9` limit. This is not a claim of bit-identical
floating-point metadata. The full base comparison is recorded in
`topology.validation--menger-broad-validation-2667.md` on PR #2925.

## Affected native units

The twelve affected units were explicitly included alongside `--green` in one
deduplicated sweep, using the pinned v1.0.1 native reference and two workers.
The blocked plan status therefore did not exclude these units from execution.

| Unit | Cases | Passed | Failed |
|---|---:|---:|---:|
| `curve.curve_p1` | 3,102 | 3,049 | 53 |
| `curve.curve_p4` | 3,263 | 3,234 | 24 |
| `curve.curve_p5` | 3,140 | 3,124 | 7 |
| `curve.curve_p7` | 3,067 | 3,051 | 16 |
| `curve.curve_p9` | 3,126 | 3,119 | 7 |
| `geom.global_solver` | 630 | 628 | 2 |
| `intersect.topology_solver_p1` | 3,592 | 3,568 | 24 |
| `intersect.topology_solver_p2` | 2,575 | 2,532 | 43 |
| `topology.brep_p1` | 2,636 | 2,634 | 1 |
| `topology.brep_p10` | 1,436 | 1,417 | 4 |
| `blend.hcr` | 481 | 463 | 1 |
| `primitive.api_trim_ops` | 211 | 196 | 5 |
| **Total** | **27,259** | **27,015** | **187** |

The remaining classifications are 15 known mismatches, 13 already parked cases
and 29 rejected cases. No new mismatch is parked by this PR.

All case records in eleven units match the historical `f79b779f5` affected-unit
reports. In `curve.curve_p5`, 32 case records changed; every one matches current
main's report. Six moved from pass to known mismatch through main's curve-length
work (`gen0561`, `gen2395`, `gen2629`, `gen2638`,
`length_raw_gap1.5e-12_acc0.0001`, `length_raw_gap1.9e-12_acc0.0001`); the other
26 changed only their reported maximum relative error and still pass. None of
the 187 failing-case payloads changed from the historical candidate.

## Body-image provenance

The ignored `export_elbow_for_diagnostics` test was explicitly rerun on the
rebased candidate and passed. Both small/control `.smb` bodies and both edge
JSON exports are byte-identical to the `f79b779f5` exports. The existing figure
therefore still depicts this candidate's body from the same view.

## Remaining limits

Zero distance-tolerance semantics and searches outside the cached domain still
need an agreed contract. The silhouette accuracy finding, curve-pair maxima,
local convergence differences and finite-trim behavior retain their individual
records; passing elbow regressions does not establish their correctness.
The finite-segment and arc vote replays explain the three trim outputs, with
their public-API/UV-binding caveat, but do not settle the intended trim contract.
All twelve affected units remain blocked and this PR remains a draft.
