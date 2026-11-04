# general: shims of C++ the SMLib fork removed or changed return `cpp_unlinked` (2026-11-03)

Units: blend.fillet_corner_p1, blend.fillet_corner_p2, blend.fillet_corner_p3,
blend.fillet_executive_p3, blend.fillet_geom_p2, blend.fillet_solver_p1,
blend.fillet_solver_p2, blend.fillet_standard_solver_p1, primitive.primitive_creation_p1.

#2502 moves the difftest's C++ reference from usd-brep v1.0.1 (`3fd5b95`) to the SMLib fork,
tangency-labs/smlib `v1.1.0-alpha.1` (`8eac63c`), in #2512. The fork replaced v1.0.1's fillet
path with its generic blend (`docs/generic-blend.md`) and changed two primitive helpers. The
Rust ports of these units follow v1.0.1, so some of their shims call C++ the fork no longer
has. Those entry points and ops now return `SC_ERR_UNLINKED` (status `cpp_unlinked`,
`harness/cpp/shim/sc_shim.h`); everything else in the shims runs against the fork. A case
that reaches one is a harness error on the C++ side, so it fails, and a `known_mismatches`
entry cannot hold it. The blend units are therefore `blocked` and leave `--green` until the
flags flip (#2502 ask 3) and the units go with their old paths, or their runners drop the
cases. primitive.primitive_creation_p1 stays in the set (as #2502 decided): its runner compares
`cpp_unlinked` as an error, so its two ops' cases are known mismatches citing the commits below
(`crates/difftest/known_mismatches/primitive.primitive_creation_p1.toml`).

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
  - Both ops' 165 cases are known mismatches; the unit is `partial`.

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

## Crashes and a hang against the fork

With the shims above, every unit builds and runs against the fork. Ten units then stopped in
the fork's C++, none of which do against v1.0.1. The sites, from a SIGSEGV backtrace (frames
resolved with `addr2line` against the libraries built with `-g`), and what #2512 does about
each:

- **surface.cone_p1, surface.cone_p3, surface.surface_p2** (the last through cone_p1's
  cases): a null dereference in `SmSurfOfRevolution::RebuildPolarConverter`, from
  `MakeNurbWithSweepParams(const SmCircle*)` (faulting address 0). It uses the iso curve
  `CreateIsoParametricCurve` returns with only an `SM_ASSERT` on the status. Any bounded cone or
  cylinder whose rotation range runs backwards does it (v1.0.1 builds them): reported as
  tangency-labs/smlib#132. The cone runners reject those inputs until it is fixed.
- **surface.surface_assert_valid:** a null read in `N_SrfEvalPt`, from `SmSurface::IsClosed`,
  on a malformed input without a NURBS net; the shim refuses such an input before the check.
- **heal.face_props:** `SmCurveClassification::GetLoopTouches`, from
  `SmFaceProps::SetProps_Stage3`. In case `gen0241` the fork's `sms_FindAndClassifyLoopBisector`
  fails (status 1001, measured with a library copy that prints it), `SmFace::ClassifyLoops` checks
  that only with `SM_ASSERT`, and the classification it goes on with holds a freed curve. Reported
  as tangency-labs/smlib#134; the runner runs the unit's scripts in a fresh child process and
  rejects a case whose child crashes.
- **topology.proto_topology_p2, topology.proto_topology_p3:** in proto_topology_p1's shim
  while their fixtures are built; p3's is `SmCurveInterval::GetCurve` on a null interval
  (address `0x18`), from `sm_EdgeCrossingSine` in `FindOrCreateProtoEdge`. Shim mismatches:
  the shim now passes `FindOrCreateProtoEdge` its classification and resolves proto edges
  with the fork's degenerate-length rule (03e7bf1).
- **primitive.primitive_creation_p2:** a failed `CreateTaperExtrude` (status 1005, where
  v1.0.1 succeeds) leaves the input curves held by the half-built B-rep's edges while they stay
  in the caller's array; the shim's delete (v1.0.1's contract) then crashes in `SmEdge::~SmEdge`.
  Reported as tangency-labs/smlib#133. In a second case (`review2_low_triangle_end_cap_normal_down`,
  added to the issue) deleting only the curves without an owner crashes as well, so the shim no
  longer deletes any of them after a failure.
- **blend.fillet_solver_p2:** `SmFilletSolver::RailEdgeuseIntersect` (address `0x20`).
- **blend.fillet_geom_p2** does not crash: case generation was still in
  `sm_SurfacePlaneCurve`, from `sm_SrfSrfIntersection`, after 30 minutes. At
  `e9f0627` it reports (1151 cases, one `cpp_unlinked`); the cause of the change is unmeasured.

The last two are in units blocked above, so they are not in `--green`. For the others the
shims and runners above settle each site, and every `--green` unit reports at `e9f0627`.

Measured in #2512 (WSL, the fork's libraries built with the harness flags).
