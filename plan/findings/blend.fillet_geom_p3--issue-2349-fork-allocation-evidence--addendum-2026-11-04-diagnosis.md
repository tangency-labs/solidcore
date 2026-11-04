# blend.fillet_geom_p3: named fork differences and self-mitre design boundary (2026-11-04)

Updates: `plan/findings/blend.fillet_geom_p3--issue-2349-fork-allocation-evidence--addendum-2026-11-04-follow-ups.md`.

Units: `blend.fillet_geom_p3`, `blend.fillet_executive_p3`.
PR: https://github.com/tangency-labs/solidcore/pull/2381

## Complete observation comparison

The probe now decodes every native wire field using the existing executive
runner's layout and compares independent Rust debug output with the saved fork
records. This closes the earlier allocation-only diagnostic limitation. It
does not fix or suppress any mismatch. Integer state is exact; real values use
relative 1e-9 with no absolute floor. The runner's separately collected
`fillet_point_step` is outside these wire records and remains unverified by
this probe.

Inputs and reference provenance are unchanged from the parent finding: the ten
`issue_2349_closed_angle.json` cases, Rust source e60349c7a, and fork a18bbe523b9c.
Rust matches every decoded observation from the configured old reference.
Against the fork, six complete records match, while four differ:

| Checkpoint | Differing fields (number of observations) |
| --- | --- |
| solver_geometry r0.5 | generated_geometry_samples (54) |
| corner_edges r0.5 | generated_geometry_samples (54), corner_edge_geometry_samples (11) |
| full r0.5 | result_topology (28), result_geometry (64), mass_properties (12) |
| full r3 | mass_properties (8) |

The first divergent surface point is at normalized (u,v)=(0.17,0.17):
Rust (-3.96873681983367, 2.004411884606559, 9.631940486795699),
fork (-3.9687558729642345, 2.002037883738891, 9.624514441544006).
All earlier vertex/ownership and parameter records match. The r3 topology,
geometry and operation status agree; its mass differences alone do not prove
a different constructed body.

Construction options are not equivalent in this comparison. The difftest's
`end_to_end::rust` creates `Kernel::default()`, with fixes off. In
`fillet_executive::p2`, a circular generator then has approximate arcs
(`!fillet_exact_arc_sections`). The fork's `SmFilletExecutive.cpp` constructs
`SmGnSectionFSG(SmArcSection(FALSE, dXSectTol))` for that circular type.
This establishes a construction-selection difference, not that it explains
every numerical difference. No construction flag may be promoted or removed
from these unmatched-option results. The eight r3 mass observations still need
their own integration/build-configuration diagnosis.

Reproduce the complete diagnostic without rebuilding the shared C++ kernel:

```sh
python3 -B harness/cpp/probes/issue_2349_fork_reference.py --analyze-report \
  reports/difftest/issue-2349-fork-reference/checkpoints.json \
  reports/difftest/issue-2349-fork-rust/blend.fillet_executive_p3.debug \
  reports/difftest/issue-2349-named-differences.json
```

Expected exit: 1, with the differences above. The output JSON retains every
individual differing field/index/value pair, including optional old-reference
comparisons. Generate the input records with the parent finding's commands.
Four Python regressions cover complete-field decoding, malformed/trailing
records, strict zero comparisons and missing-field/length rejection. They were
committed failing before implementation and now pass in Ubuntu WSL.

## Fit-size and exact-pcurve design dependency

There is a concrete lead for the size contract: `ball_regular::ball_is_regular`'s
`limit_at` bisects for a regular lower endpoint, then can replace it with
`s_lo.fits`, the curvature limit. At the teardrop profile's t=1/2,
|C'|=15 and |C' cross C''|=1800, hence the curvature radius is
15^3/1800=1.875. The previously reported fitting size 1.8749999999999993
is near this open regularity boundary. This is not yet a proof of the cause
of every invalid face or a certified safe replacement radius.

The existing `blend_corner::exact_mitre_curve` constructs a rational conic
with exact UV curves for two equal cylinders on straight spines. It explicitly
does not construct a curved-spine self-mitre. The latest-main smooth End
(`blend_corner/end.rs`, including #2364/#2406) calls the existing SSI machinery;
it does not establish a shared certificate for this self-mitre and both UV maps.
Deleting the old closed corner or relabelling it as an event would not supply
the missing construction.

The exact/certified intersection representation and the fitting-size contract
at an open curvature limit are awaiting a design answer on the blend hub:
https://github.com/tangency-labs/solidcore/issues/2243#issuecomment-6622359026

The question follows AGENTS.md §4b's requirement to ask before building when
the docs and interfaces do not determine the construction. #2307 and #2353 are
merged. #2401's sizing work is closed via #2431; #2404 covers other near-bound
assembly problems, so no competing shared fit-margin implementation is claimed.
This contribution changes diagnostic tooling only, retaining production flags,
all ten checkpoints, tolerances, allocation assertions and the two visibly
unfinished ignored contracts.
