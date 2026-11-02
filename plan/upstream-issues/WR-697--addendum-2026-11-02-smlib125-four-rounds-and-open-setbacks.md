Updates: WR-697, WR-616

#### WR-697 / WR-616 (update, 2026-11-02). Corners of four or more rounds are the Ball or the Setback (`fillet_blend_corner`, smlib#125)

- **Source:** the SolidCraft fork's tangency-labs/smlib#125, merged into `replay/1.0.1` at
  `18e96f1`. This is a port of its SMLib commits, not solidcore-native work. It runs under
  the existing opt-in `Kernel::fixes.fillet_blend_corner`, with no new flag and no new
  record id. The fork's capi commits (`fc29fdb`, `cb2a478`, `96e58a1`) and the capi parts
  of the others are not ported. They were read for the cases' inputs and expected values.
- **What happens (pinned and before):** a vertex of four or more equal rounds has no blend
  corner. The pinned `SmFilletConvexNxNCorner` (or `SmFilletConcaveNxNCorner`) builds it,
  and its patch only approximates the ball the rounds lie on. WR-502 is the square pyramid's
  case of this.
- **Four or more rounds** (`278bf77`):
  - `SmBlendCorner::ClassifyBall` takes `n >= 3` rounds; for `n >= 4` the faces must all be
    planes. The rounds' offset planes must pass through the point where the first three
    meet (`sm_BallPlanesMeet`, by Cramer's rule, within the rounds' tolerance).
  - `SmFilletCorner::Create` asks for it in the NxN case before `ClassifySetBack`.
  - The Ball's faces are the sphere's n-gon fanned from contact 0 into `n - 2` triangles,
    each the Ball's rational (4, 2) triangle patch with its pole at contact 0
    (`BuildBallSurface`, `CreateBallPatch`). `MakeBallTopology` and `BallCentre` take `n`
    rounds.
  - `ClassifySetBack` leaves only three equal rounds to the Ball, so four or more equal
    rounds without a common ball are a Setback.
- **Exact zeros in the convexity certificate** (`3e10318`): along a setback patch's side
  that is G1 to a face plane by construction (its first two rows of control points on the
  plane, `sm_SbSidesOnPlanes`), `L = -S_u . N_u` and `M = -S_v . N_u` are identically 0. So
  `sm_SbCertifiedConvex` sets their Bernstein rows there to exactly 0 instead of reading
  them from rounding.
- **WR-502** (`fillet_corner_sphere_from_four_arcs`) is kept as it is. With
  `fillet_blend_corner` a four-round apex with a common ball is the Ball's patches before
  `MakeAnalyticCornerSurface` is reached, so the WR-502 sphere is not built there. Without
  `fillet_blend_corner`, WR-502 behaves as before.
- **Rust port:**
  - `crates/solidcore/src/blend/blend_corner/ball.rs` (`classify_ball`, `ball_centre`,
    `ball_corner_edge_geom`, `create_ball_patch`, `planes_meet`; `BallCorner` holds the
    patches and their loops).
  - `crates/solidcore/src/blend/fillet_corner/topology/ball.rs` (`n` rounds).
  - `crates/solidcore/src/blend/fillet_executive/session.rs` (the dispatch, `n` solvers).
  - `crates/solidcore/src/blend/blend_corner/setback.rs` (`n == 3` in `classify_setback`).
  - `crates/solidcore/src/blend/blend_setback/certify.rs` (`sides_on_planes`,
    `certified_convex`) and `graph.rs` (the call).
- **Default unchanged:** without `fillet_blend_corner` every corner is the pinned one.
- **Tests:** `crates/solidcore/tests/fork_smlib125_four_round_corners.rs`, from the fork's
  capi cases (`ball regular_apex_4/5/6`, `setback four_face_apex`,
  `setback four_face_pocket_floor`) and WR-502's pyramid. They are checked by closed form and
  geometry:
  - `check_brep` Level2 passes, and the face counts are the sharp body's plus `2 n - 2`
    (Ball) or `3 n` (Setback).
  - Each Ball patch lies on the closed-form sphere (centre on the axis, `r / n_z` below the
    apex) to 1e-12 of the part's size, along its edges too. Its normal is radial, and it is
    G1 to every neighbour along its edges (1e-9 rad).
  - With the flag off, the regular apex of four is the pinned corner.
