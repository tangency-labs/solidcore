# general: shims of C++ the SMLib fork removed or changed return `cpp_unlinked` (2026-11-03)

Units: blend.fillet_corner_p1, blend.fillet_corner_p2, blend.fillet_corner_p3,
blend.fillet_executive_p3, blend.fillet_geom_p2, blend.fillet_solver_p1,
blend.fillet_solver_p2, blend.fillet_standard_solver_p1, primitive.primitive_creation_p1.

#2502 moves the difftest's C++ reference from usd-brep v1.0.1 (`3fd5b95`) to the SMLib fork,
tangency-labs/smlib `v1.1.0-alpha` (`a18bbe5`), in #2512. The fork replaced v1.0.1's fillet
path with its generic blend (`docs/generic-blend.md`) and changed two primitive helpers. The
Rust ports of these units follow v1.0.1, so some of their shims call C++ the fork no longer
has. Those entry points and ops now return `SC_ERR_UNLINKED` (status `cpp_unlinked`,
`harness/cpp/shim/sc_shim.h`); everything else in the shims runs against the fork. A case
that reaches one is a harness error on the C++ side, so it fails, and a `known_mismatches`
entry cannot hold it. The units are therefore `blocked` and leave `--green` until the flags
flip (#2502 ask 3) and the units go with their old paths, or their runners drop the cases.

**Inputs:** every case of the ops below. **C++:** `cpp_unlinked`. **Rust:** unchanged (the
v1.0.1 port).

## Not linked, by entry point or op

The commits are tangency-labs/smlib's.

- **blend.fillet_corner_p1**
  - `sc_blend_fillet_corner_p1_helper` ops `SEAM_OUT` and `CLOSURE`: SmFilletCorner.cpp has no
    `sm_CheckSeamOut` or `sm_SetSurfaceClosure` (caa67d0; SmFilletGeom.cpp's copies remain).
  - Same entry point, op `NX1_DEGENERACY`: no `SmFilletNx1Corner` (53c47b8).
  - `sc_blend_fillet_corner_p1_corner` (every op): each builds a corner of one of v1.0.1's
    derived corner classes. The fork has none: 1x1 and 2x2 (543d25a), Nx1, open and
    degenerate (53c47b8), Nx2 (74e9beb).
- **blend.fillet_corner_p2:** its `UpdateSplitFace` and `MakeAnalyticCornerSurface` cases run
  through `sc_blend_fillet_corner_p3_corner` (below).
- **blend.fillet_corner_p3:** `sc_blend_fillet_corner_p3_create` (SmFilletCorner::Create,
  whose solvers carry the cross-section surface generators) and
  `sc_blend_fillet_corner_p3_corner` (the derived corner classes, as above).
- **blend.fillet_executive_p3:** `sc_blend_fillet_executive_p3_two_offset_chamfer`
  (solidcore#880) gives the executive an `SmLinearCrossSectionFSG` (3aae1c4).
- **blend.fillet_geom_p2:** op `CROSS_SECTION`, `SmFilletVertex::CalcVertOnCrossSection`
  (543d25a).
- **blend.fillet_solver_p1**
  - `sc_blend_fillet_solver_p1_eval` op `FSG`: the linear and circular cross-section
    generators (3aae1c4) and the blend-curve one (0bd3f4d).
  - `sc_blend_fillet_solver_p1_static_eval` ops `CHECK_DOMAINS`, `ANALYTIC` and `ANALYTIC_REV`:
    the fork's SmFilletSolver.cpp has no `sm_TestAnalyticSurface` (7ac037e), and
    `sm_CheckEdgeCurveDomains` is a static of its SmFilletExact.cpp.
- **blend.fillet_solver_p2:** op `GEOM`, `CalcCylinderFilletGeom`, `CalcConeFilletGeom` and
  `CalcTorusFilletGeom` (7ac037e).
- **blend.fillet_standard_solver_p1:** ops `CIRC`, `BLEND` and `LINEAR` (the three
  generators, 3aae1c4 and 0bd3f4d), and `VR` and `VR_SETUP` (their SmVariableRadiusFS needs a
  generator; v1.0.1's `OffsetRadiiCanChange` reads its type).
- **primitive.primitive_creation_p1**
  - Op 1: `sm_FindPivotPoint` takes the corner's offset curve in place of the tolerance
    (69d8857). The Rust port has v1.0.1's signature.
  - Op 5: `sm_CreateAnalyticRotationalSweep` is gone (9853e54).

## Adapted, not unlinked (same inputs, v1.0.1's behaviour)

- `sc_blend_fillet_solver_p2_eval` op `RAIL_RAIL`: the fork's `RailRailIntersect` takes
  `bHoldEndRadii`. `TRUE` sets up `SmTemporaryRadiusChange` on both solvers, as v1.0.1 always
  did.
- `blend_fillet_geom_p1.cpp` includes `SmBlendCorner.h`: the fork's SmFilletGeom.cpp calls
  `sm_SideEdgeuseAcrossSeamVertex`, declared there.
- The ND-3 parent check (`check_parent_static_step`) accepts the fork. The fork picks the
  finite-difference step for each point (9e17f64), so it has no `sdStaticStep` to latch.
  Without that, blend.fillet_executive_p3, blend.fillet_intersector and blend.fillet_solver_p4
  stopped with a harness error before their first case.

Measured in #2512 (WSL, the fork's libraries built with the harness flags).
