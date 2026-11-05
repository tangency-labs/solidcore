<!-- SPDX-License-Identifier: Apache-2.0 -->

# Remaining duplicate groups from #2722

This audit starts after #2721's 124 removals. The comparison uses exact `calls`,
`bodies` import options and decoded operand records, without rounding. A different
case description or expected answer cannot make identical replay inputs exercise
another kernel path.

The application-side checks below use SolidCraft commit
`8f2645748`. Its kernel pin is `8ddce901671f1cebd9214a31d21895d725e886c8`,
the same kernel pin as the original corpus. The original recipes remain in
SolidCraft; removing a redundant kernel entry does not remove an application test.

## Verified redundant kernel entries

| Removed entry | Retained entry | Evidence |
|---|---|---|
| `fillets/1207-fillet-bar-top-edge-and-boss-arc-together-reduce-fit` and `-taper-fit` | `fillets/1207-fillet-bar-top-edge-and-boss-arc-together` | `doc::eval::fillet_fitted` first calls `fillet_removing` with the requested radius for every fit mode. Reduction or taper runs only after `TooLarge`; this radius succeeds. Calls, operand records and measured results are identical. |
| `fillets/1207-fillet-bar-top-edge-and-boss-arc-together-curvature-corners` | `fillets/1207-fillet-bar-top-edge-and-boss-arc-together` | The recipe's own expected description says no vertex has three rounded edges. `kernel::blend::g2_corners` calls ordinary `fillet` when the corner list is empty. |
| `fillets/4782-fillet-full-round-two-end-verticals-reduce-fit` and `-taper-fit` | `fillets/4782-fillet-full-round-two-end-verticals` | The requested radius already fits. The same `fillet_fitted` success branch returns before consulting fit mode. The separate case `4782-fillet-two-end-verticals-past-a-full-round-reduce-fit`, where the radius must change, stays. |
| `fillets/every-box-edge-r3-curvature` | `fillets/cov-fillet-section-curvature-box-edge-r3` | Recipes are identical. Retain the entry with the stricter volume tolerance and exact G3 continuity expectation; do not weaken either expectation. |
| `surfaces/cov-extrude-line-symmetric-total-distance` | `surfaces/cov-extrude-line-symmetric` | The current command defaults to `total_distance`. Both ask for distance 10, producing the same profile displacement and extrusion distance. The retained expectation explicitly requires z = -5..5. |
| `sweeps/771-sweep-two-guide-rails-cone-frustum` | `sweeps/771-sweep-one-guide-rail-cone-frustum` | The second rail is the first rail reflected across the straight path. `kernel::sweep_path::railed_stations` converts both descriptions to the same affine stations; `smlib::Brep::sweep_stations` records those stations in full. Both profiles, all station matrices/translations, end modes and measured results are identical. No rail parameter is consumed by this kernel operation. |
| `booleans/3887-ring-standing-on-perforated-plate` | `booleans/3092-ring-on-perforated-plate` | Explicit XY frames versus named XY planes, and different picks inside the same plate region, produce identical operands and union calls. Both analytic methods describe the same plate and ring. Retain the more precise analytic volume and carry over the removed source reference. |
| `solids/cov-align-side-face-turns-onto-top-face` | `solids/cov-align-bottom-face-onto-top-face` | Both already have `operation_expressible: false`: SolidCraft changes the aligned body's record outside traced operations. Their calls construct only the same plate. The remaining entry still reports the missing aligned body as unsupported; this is not a fix or successful replay of alignment. Both application alignment recipes remain in SolidCraft. |
| `tessellation/1666-loft-18-gon-ruled-sides` | `sweeps/1666-loft-zigzag-to-flat-rectangle` | Both create the same ruled loft. The tessellation entry already marks its watertightness check ungraded: no mesh operation is in the trace. Retain the more precise volume and bounding-box expectation, and carry over source #1670. The application mesh case remains in SolidCraft; no mesh coverage is claimed here. |

No removed entry occurs in `known_failing.json`. Expected numeric values and
tolerances of retained entries are unchanged. The manifest is rebuilt from the
remaining entries; the same reconstruction was first checked against the
unmodified manifest.

## Distinct operations that must stay

- The two mitre/swept pairs must retain both cases. SolidCraft commit `f74c6605c`
  adds the `turn` command option and records it in `blend_dropping`; the original
  corpus predates it. Refreshing those exports must record different `turn`
  values and preserve each case's analytic expectation.
- `ladder/02-booleans/007-near-tangent-gap-1e-3` must retain its tool construction.
  The complete-trace exporter and handled-refusal replay are owned by
  [SolidCraft #12099](https://github.com/tangency-labs/solidcraft/pull/12099) and
  [solidcore #2742](https://github.com/tangency-labs/solidcore/pull/2742).

## Refresh provenance

Only the four `cov-fillet-turn-*` entries above and the near-tangent ladder entry
are refreshed. Their `expected` objects are checked equal to the original exports;
no analytic value or tolerance changes. Each refreshed entry records
`export_provenance`; the manifest's original pins still describe all other entries.

The export uses SolidCraft `8f2645748` with the exporter patch from
`0a79d0874400c8508f5bb6df1406b783f47090ca` (SolidCraft #12099). To keep the
original corpus's fix set, set `fillet_through_extended_rail_end` and
`fillet_round_on_round_as_definition` to `false` in `crates/smlib/src/core.rs`
before building. This is a local export configuration, not an application change.
The generated `fixes.json` must match the original corpus's flags exactly.

Build `cargo build --locked --profile ci -p solidcraft-cli`, then run one process
per case (trace body IDs are process-local):

```sh
solidcraft-cli kernel-cases tests/kernel-cases/<area>/<case> --export=DIR --json
```

Move each inline `bodies.<id>.record` into the sibling `.bodies.json.gz` (gzip
timestamp zero), leaving the import options in the entry. Preserve the ladder's
existing `golden` metadata. Rebuild manifest counts from the resulting entries.
There is no duplicate allowlist or tolerance in the uniqueness guard.
