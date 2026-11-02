# blend.fillet_solver: Discord 70 mm top-rim case (2026-11-02)

Units: `blend.fillet_solver`, `step-brep`.
PR: https://github.com/tangency-labs/solidcore/pull/2390.

## Input and selection

The user supplied Claas Kuhnen's `fillet bug.step` from the Discord comparison
dated 30 October 2026 and confirmed the upper model's requested radius as 70 mm.
The lower screenshot's 2 mm request is a different model and is not covered.
The fixture is `crates/step-brep/tests/fixtures/discord-2026-10-30-fillet.step`,
SHA-256 `beb34135184eaed70872c314ba3a76b39c103c730dd79f3aea553c23978f71bb`.
It is AP214 in millimetres, exported by smlibpy/SMLib.

The reader produces one body, 15 faces, 39 edges and 26 vertices without skipped
geometry. The topology-first build reports no face failures. The cylinder radius
is 304.8 mm and its top is z = 458.266019821167 mm. The four convex top outer arcs
are STEP #163, #170, #185 and #198. The regression selects these geometrically.
This is an inferred selection: STEP does not preserve the original command, and
the screenshot does not establish which arcs were selected together.

## Reproduction and existing correction

On solidcore `23ea17c1c21a587e5bd50653e7c73c490f34bbbe`, fresh bodies were used
for each radius and flag profile. The four-arc operation fails with default fixes
at radii 1, 10, 20, 40, 60, 70, 80 and 100 mm. Preparation succeeds; `do_fillet`
returns status 1001 with error code 3 and comment
`Fillet Corner Error: Corner Edge Computation Failure`.
At R70, each of the first three arcs succeeds individually, while the fourth
fails. Thus these trials do **not** demonstrate a radius threshold or establish
the cause of the screenshot's visible distortion.

Delta reduction of the successful investigative flag profile leaves only
`Kernel::fixes.fillet_exact_round_classes` (WR-669). That flag alone makes the
four-arc R70 operation complete. The four resulting round surfaces are surfaces
of revolution. An 11 by 11 parameter grid on each lies within 2.6e-13 mm of the
expected torus:

```text
(sqrt(x*x + y*y) - (304.8 - 70))^2
    + (z - (458.266019821167 - 70))^2 = 70^2
```

The regression checks all 484 points against the radius at 1e-9 mm absolute.
This finite surface check is not a certificate for the trims, junctions or solid.
No new kernel flag or geometry algorithm is introduced. WR-669 classifies from
the supporting faces rather than the old edge-shape-driven try chain; this case
shows that its existing correction changes operation completion. The precise
failing branch inside the old corner construction has not been isolated.

## C++ control

The imported operand was written as native SMLib ASCII and read independently by
a standalone SM_API driver linked against the pinned prebuilt SMLib reference
(`upstream_commit=3fd5b9548d19aa8c12c29bc9637a96d42a06f8a5`, g++ 13.3,
`-ffp-contract=off`, `-fno-fast-math`). At R70 it selects the same four arcs
geometrically and `SmApiFilletEdges` returns 1001. Selecting only the fourth arc
also returns 1001. This confirms the default failure on the **imported operand**;
it does not independently validate STEP import or compare against the later
fork's opt-in implementation. No C++ output or state is supplied to Rust.

## Unresolved validity

Before the fillet, `check_brep(Level2)` reports 24 failed edge/vertex tolerance
reports (approximately 0.0233 mm versus its 0.01 mm limit) and skips the next
geometry check. Its `failure_count()` is 25 including the synthetic skipped-check
line. The same counts remain after the flag-enabled operation.

Calling `assert_valid` directly as a separate diagnostic also returns false:
besides those tolerance reports it reports two `Face has Coincident Edges` and
eight `Bad Vertex Gap` failures on both operand and result. The regression does
not promote operation success to whole-body validity. Whether these arise from
the original model, its export, or import needs independent investigation.
No tolerance was changed and no healing/retry was used to obtain success.

Next work: establish input validity independently, recover the actual selected
arcs, and compare trimmed junctions against the Parasolid result. Do not declare
the Discord visual defect fixed from these status and surface samples alone.
