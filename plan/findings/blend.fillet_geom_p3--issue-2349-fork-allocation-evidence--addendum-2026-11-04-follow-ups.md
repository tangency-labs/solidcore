# blend.fillet_geom_p3: remaining #2349 verification and construction work (2026-11-04)

Updates: `plan/findings/blend.fillet_geom_p3--issue-2349-fork-allocation-evidence.md`.

Units: `blend.fillet_geom_p3`, `blend.fillet_executive_p3`.
PR: https://github.com/tangency-labs/solidcore/pull/2381

## Native differences requiring diagnosis

Input: `fixtures/blend.fillet_executive_p3/issue_2349_closed_angle.json`.
The four differing complete records are:

- `issue-2349-teardrop-solver_geometry-r0.5`
- `issue-2349-teardrop-corner_edges-r0.5`
- `issue-2349-teardrop-full-r0.5`
- `issue-2349-teardrop-full-r3`

These compare the configured usd-brep 1.0.1 executable with fork
`a18bbe523b9c7819f1a9da8c40a54c15f9936125`; they are not a complete
Rust-versus-fork comparison. The parent finding records reproduction commands
and sample values: generated surface sample -3.96873681983367 versus
-3.9687558729642345, and mass-property record 457.03059824978146 versus
457.03060521148683. The r0.5 completed result also differs in topology ordering.
The allocation fields independently compared with Rust agree bit for bit.

Next evidence needed: decode the first differing geometry/topology field in
each record, compare independent Rust output with the fork under matching
construction options, and trace the responsible construction. Check reference
build configuration, including random-seed initialization, before attributing
mass-property differences. Redesigned geometry is a suspected cause, not a
demonstrated diagnosis. Completion requires unchanged parity tolerances and
retained checkpoint fields; neither old-reference agreement nor the passing
allocation subset proves complete fork parity.

## Unfinished construction contracts

The following tests in `crates/solidcore/tests/issue_2349_closed_angle.rs`
remain explicitly ignored diagnostic contracts:

| Contract | Current evidence | Required completion evidence |
| --- | --- | --- |
| `closed_angle_is_a_blend_event_instead_of_the_legacy_corner` | The reachable closed rim still uses the legacy closed-corner mitre | Certified curved-spine self-mitre event reached on this operand |
| `closed_angle_round_builds_validly_or_refuses_with_a_fit_size` | Radius 3 refuses with suggested size 1.8749999999999993; the suggested boundary still constructs invalid geometry | Valid construction within the proved fitting range and correct refusal outside it, including boundary tests |
| Exact pcurves on both sides | Not implemented or certified by the allocation correction | Exact supporting-surface parameter curves and their consistency with the shared mitre edge |

Reproduce the two known failing contracts in Ubuntu WSL, after checking shared
machine load, with Cargo on PATH:

```sh
export CARGO_BUILD_JOBS=2 CARGO_INCREMENTAL=0 RUST_TEST_THREADS=2
nice -n 10 cargo test --profile difftest -p solidcore \
  --test issue_2349_closed_angle -- --ignored --nocapture
```

This command is expected to fail until those contracts are implemented. Do not
remove ignores as a substitute for implementing them or loosen their checks.
#2307 and #2353 are merged; they are no longer external merge blockers. The
remaining dependency is the certified construction and its proof, not the
split-vertex allocation. #2499's successful rounding complements this fixture
but does not certify these contracts. Check current claims on #2349/#2243
before implementation; this record reserves no follow-up scope.
