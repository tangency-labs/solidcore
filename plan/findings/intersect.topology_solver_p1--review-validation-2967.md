# Silhouette review validation (#2967)

Production snapshot: `cac31e200`, including main `d3cc2f0aa`, the shared
UV-only cache reader in `0bf087f39`, and parent #2840's affine correction
`a62556200`. The later `5707c8472` changes only audit-helper docstrings.

## Executed checks

- `cargo fmt --all --check`: passed.
- `cargo test --profile goldens -p solidcore --lib silhouette_accuracy_tests`:
  five passed, zero failed. The tests cover the populated cache, original
  sphere, tilted sphere at three projection scales, trimmed sphere endpoint,
  and torus capsule projection.
- The resulting library test executable, run with two test threads:
  **3,615 passed, zero failed, four ignored**.

The cache test inserts a fitted 3D curve displaced from the surface and checks
that the shared reader instead evaluates the cached UV curve on the source
surface. It does not assume that a second topology solve retains the first
solve's local cache. The geometric tests compare the query against closed
forms and independently compare its scalar with the surface evaluated at the
returned parameters.

These are current production-code checks. The full integration run that began
at `14404b56f`, the queued refreshed workspace checks, and the separate affine
integration regression are not reported as passing here before completion.

## Remaining construction

No general UV-trace, trim, or candidate-completeness bound follows from these
tests. For a rational patch `S = X / W`, with `W` bounded away from zero, the
silhouette has the zero set of the polynomial

```text
N dot ((W X_u - X W_u) cross (W X_v - X W_v)).
```

Bernstein coefficient and derivative bounds offer a route to branch coverage;
the trim endpoints and objective variation still require enclosures in the
query metric. The existing `ssi::trace2d` field-cell construction is private
and coupled to surface/surface intersection. The kernel currently has no
dependency on that crate, and `docs/layers.md` requires consolidation rather
than a competing implementation. The integration boundary was raised on
[#2967](https://github.com/tangency-labs/solidcore/issues/2967#issuecomment-6695067203)
before starting that extraction. Neither smaller tracing requests nor lifted
measurements alone satisfy the full certificate requirement.

## Parent CI

The [golden job at parent `a62556200`](https://github.com/tangency-labs/solidcore/actions/runs/42289593125/job/128144443176)
failed in the same twelve units covered by the historical caller audits, with
zero broken cases. This does not establish that all changes are correct:
those records still identify unresolved accuracy and query-contract questions.
No golden was blessed by this validation work. The parent remains ready for
review, not certified ready to merge.
