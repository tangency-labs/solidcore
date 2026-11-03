Updates: WR-697

#### WR-697 (update, 2026-11-03). Unequal rounds at a vertex beside a cylinder, cone or sphere end on a Setback whose fill is certified no more concave than what it is G1 to
- **Severity:** Wrong result
- **Area:** Fillets (corners)
- **Location:** `source/SMLib/src/SmFilletCorner.cpp`, `SmFilletCorner::Create`. Pinned SMLib builds such a vertex as an `NxN` corner (`N_FillNSidedHole`), not the Setback of smlib#104, whose `ClassifySetBack` takes only planes.
- **What happens:** three unequal constant rounds meet at a vertex where a face is a natural quadric, for example the D prism's top corner, where the top, the flat and the convex cylinder meet, with the chord, the upright and the arc rounded 2, 2.5 and 3. No one ball touches the three faces there. The corner is the old N-sided patch, whose fit is not certified against the curved face.
- **Source:** smlib#85 part C, built in solidcore (#2343): solidcore-native, no fork commit. The fork's work in progress `archive/wip/85c-curved-setback` (tangency-labs/smlib 37c9e1d, `SmBlendSetback.cpp`, `SmBlendCornerSetback.cpp`) served as the reference and the starting point, not as a bit-for-bit target.
- **Rust port:** under the same opt-in `Kernel::fixes.fillet_blend_corner`.
  - **Classification** (`blend::blend_corner::setback::classify_setback`): a corner face may be a natural quadric (`corner_quadrics`, the fork's `sm_SbCornerPlanes` with `pbCurved` and `sm_SbQuadricFace`, from `BlendQuadric::from_surface` at the middle of the face's uv domain). Its normal out of the material at the vertex gives its side. Beside a curved face any edges are taken, but only constant circular rounds and a G1 fill.
  - **Setbacks** (`set_backs`, `curved_least`, `along_edge`): the least setbacks are where the rails of the two rounds beside each face cross on the face as it is. That point solves the faces' closed-form distances by Newton. Each is measured along its edge from the plane square to the edge through the ball's centre. The rule is `max(1.5 r, least + 0.25 r)`, held within each edge's room.
  - **Fill** (`blend::blend_setback::curved`): the graph fill over the corner's axis, with three changes:
    - each spring curve is exact on its face, through the face's rational parameter plane (`spring_curve`);
    - the column next to it is G1 to the face by construction (Hahn's connecting functions);
    - the certificate is the owner's ruling (a) on #85: `II + |k| I` PSD per patch, where k is the most negative normal curvature of its faces and of its round's end row. A setback patch beside a face where `k < 0` is certified by the faces' normal field (`no_more_concave`).
    - Degrees are tried least first, 3 to 9. The rounds' end rows come from their geometry: the arc about the ball's centre.
  - **Growth** (`build_setback_fill`, `grow_curved`): a fill refused as uncertified grows the setbacks with the real rounds: twice first, then bisecting while the certifying end is more than 1.01 times the failing one. Each trial places the setback vertices again (`SetbackSolvers::calc_set_back_vert`, the executive's `CalcSetBackVert`). If none up to twice certifies, the corner is refused naming none. An asked setback that certifies only when grown is refused naming the setback that builds (`fit_setback`).
  - **Containment** (`blend_setback::containment`): beside a curved face it is certified by the faces' implicit forms in Bernstein form. Along a column on the face the form is divided by `v^2`. A fill not certified on the material's side is refused.
  - **Deviations from the fork's work in progress:** the curved faces are found apart from `corner_planes`. The fork's `SBSCALE` environment override is not taken. The fork's process-wide `CurvedCounts` are per-fill `FillStats` fields. A failing call within a growth trial is that call's error; the fork counts it as not certified.
- **Triggers:** the fork's capi cases `setback curved_d_prism_corner`, `curved_d_pocket_corner`, `curved_knife_edge_corner` (SolidCraft #851's part) and `curved_sphere_corner`, rebuilt from solidcore calls.
- **Confirmed:** in solidcore, graded on the faces' closed forms, not on a net. Every case passes `check_brep` Level2 with the Setback's face count. Each of its 2 n patches is one Bezier span.

  | Case | Faces | Setbacks | Worst seam (rad) | Spring edges off the face | Least principal curvature (k) | Past the faces |
  |---|---|---|---|---|---|---|
  | D prism (convex cylinder) | 13 | 3, 4.5, 3.75 | 5.7e-12 | 3.6e-15 | 9.2e-3 (0) | 3.6e-15 |
  | D pocket (concave cylinder) | 18 | 3, 4.5, 3.75 | 1.4e-11 | 9.8e-15 | 8.1e-3 (0) | 7.1e-15 |
  | Sphere (convex) | 14 | 2.25, 3, 3.75 | 5.2e-12 | 5.3e-15 | 9.0e-3 (0) | 3.6e-15 |
  | Knife edge (concave wall R 21) | 15 | 3.312, 3.75, 4.5 | 1.6e-11 | 3.2e-14 | -1.67e-2 (-1/21) | 3.2e-14 |
  | Cone frustum (convex cone) | 15 | 1.875, 2.25, 3 | 4.3e-12 | 4.6e-15 | 4.8e-3 (0) | 1.5e-15 |

  Each case certifies at the rule's setbacks, which match the fork's comments for these cases (D prism 3, 3.75, 4.5; knife 3.312, 3.75, 4.5).
  - **Without the flag:** the D prism's corner is not a Setback.
  - **The convex cone** (the fork's `curved_cone_corner`, known failing there at 37c9e1d) passes.
  - **Open:** the countersink (a concave cone) grows to 1.414 x the rule's setbacks and certifies, but a patch's control net fails `check_brep`'s net test #5. This is the fork's open item (b) (`plan/findings/blend.blend_corner--curved-setback-countersink-net-doubles-back.md`).
  - **Also open:** a G2 fill beside curved faces (#2346, #2347), and N >= 4 (#2365).
  - **Growth:** the D prism asked 3.05 (at or past every least setback, 3 the most) is refused, naming 3.240625 (3.05 x 1.0625) as building. That one, asked for in turn, is refused naming 3.26594. The trial and the fresh build place the same vertices bit for bit, and their fill inputs agree to 1e-15, but the certificate's decision at the threshold flips (`plan/findings/blend.blend_corner--curved-setback-growth-named-setback-refused.md`, open).
- **Tests:** `crates/solidcore/tests/issue_2343_curved_setback.rs`: the five cases and the flag-off case, the countersink ignored as known failing; the asked setback's refusal and the named setback building, both ignored as slow, the latter known failing. The spring curves on a cylinder, cone, sphere and plane are checked against the faces' closed forms in `crates/solidcore/src/blend/blend_setback/spring.rs` (unit tests).
