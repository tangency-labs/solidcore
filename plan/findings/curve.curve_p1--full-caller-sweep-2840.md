# Full projected-extrema caller sweep (#2840)

Historical library candidate `09495e082`, compared with exact main
`3f82278a14c1a6104db8e090a59b4bcf677ace01`. The reference is the unmodified pinned
v1.0.1 build, key `v1-3fd5b9548d19-726643fb17a52b7d`. The later main merge changes
only the operation rig's known-failure data, not this kernel or these fixtures.

The complete `--green --jobs 2` diagnostic traversed **280 units / 591,949 cases**:
248 units passed, 15 were partial, and 17 failed. There were 190 failing cases.
Six units / 17 failures already failed on the exact base. The remaining eleven
units / 173 failures are listed below; a changed output is not automatically a
corrected native defect.

| Unit | Failed cases | Evidence record |
|---|---:|---|
| `curve.curve_p1` | 45 | [Global-query audit](curve.curve_p1--projected-extrema-audit-2840.md) |
| `curve.curve_p4` | 21 | [Derivative and local convergence](curve.curve_p4--projected-derivative-local-convergence-2840.md) |
| `curve.curve_p5` | 7 | [Derivative and local convergence](curve.curve_p4--projected-derivative-local-convergence-2840.md) |
| `curve.curve_p7` | 16 | [Property stationary points](curve.curve_p7--projected-property-stationary-audit-2840.md) |
| `curve.curve_p9` | 7 | [Curve-pair maxima](curve.curve_p9--projected-pair-maxima-audit-2840.md) |
| `intersect.topology_solver_p1` | 23 | [Local topology queries](intersect.topology_solver_p1--projected-endpoint-callers-2840.md) |
| `intersect.topology_solver_p2` | 43 | [Analytic and subdivided callers](intersect.topology_solver_p2--projected-extrema-callers-2840.md) |
| `topology.brep_p1` | 1 | [Endpoint side query](topology.brep_p1--endpoint-side-query-2840.md) |
| `topology.brep_p10` | 4 | [Raw affine parameter conditioning](curve.curve_p1--affine-parameter-conditioning-2840.md) |
| `blend.hcr` | 1 | [Tolerant endpoint classification](blend.hcr--projected-endpoint-classification-2840.md) |
| `primitive.api_trim_ops` | 5 | [Trim caller audit](primitive.api_trim_ops--projected-endpoint-callers-2840.md) |

The six preexisting failing units are `blend.fillet_executive_p3` (1),
`geom.hash_table` (1), `heal.gap` (7), `primitive.primitive_creation_p2` (2),
`topology.assert_valid` (2), and `topology.brep_p9` (4). Their case names match
the exact-base report; this change does not repair them.

## Limits and subsequent work

This sweep does **not** include eleven fixtures committed after `09495e082`:
four projected endpoint queries, two local/global maximum counterparts, three
maximization derivatives, and two projected branch bounds. They were separately
run against the same native reference and produced the recorded disagreements.
Their addition must be reflected in subsequent consolidated unit counts.

The 45 global-curve cases meet their original accuracy under the recorded domain
qualification. Broader callers expose unresolved behavior: a sphere-band answer
outside its tolerance, six curve-pair maximum cases with explicit better witnesses,
three trim results without a geometric correctness certificate, local convergence
differences, and incomplete or domain-qualified searches. Those findings prohibit
a blanket native-defect classification or a readiness claim.

The same candidate passed formatting, workspace clippy and wasm checks, all ten
then-existing focused tests (one export helper ignored), and all **3,564 library
unit tests** (four ignored). The full integration-test run was deliberately stopped
after the affine parameter-conditioning defect was discovered; it is **incomplete**,
not a passing full-suite result. An earlier debug-profile run was also cancelled.

The subsequent affine construction changes the library again. This historical
sweep is an impact baseline, not validation of that new candidate. No tolerance
was relaxed, no new mismatch was parked, and generated reports are not committed.
