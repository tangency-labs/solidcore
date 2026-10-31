#### Updates: WR-77 (2026-10-31). Singular-crossing corpus, and the fork's fix

- **Difftest corpus (issue #393):** `fixtures/surface.surface_p2/ssi_singular_crossings.json`
  (`harness/fixturegen/ssi_singular_crossings.py`). It has 14 `gssi` cases: 7 pairs of
  exact NURBS surfaces at the ssi corpus tolerance 1e-6, each in both operand orders.
  The rational cylinders, cones, spheres and tori are 9-point rational circles swept or
  revolved, each circle starting 20 degrees off the frame axes so that no node lies on a
  seam or a knot. All 14 are at parity with the pinned C++.
- **What upstream returns** (C++ and the default Rust kernel alike;
  `crates/solidcore/tests/issue_393_singular_crossing_corpus.rs` measures it against the
  closed forms). The point set is complete in every case: each exact curve is covered end
  to end, every returned curve lies on the exact ones (within 1e-5, the closed forms'
  sampling) and its UV curves on both surfaces, and all curves are `CROSSING`. At the nodes:

  | pair (`sx_<name>`) | node | normals | upstream, both orders |
  | --- | --- | --- | --- |
  | `steinmetz_perpendicular` | `(0, +-1, 0)` | parallel | two curves pass through each node, none ends there |
  | `steinmetz_oblique` (axes at 60 degrees) | `(0, +-1, 0)` | parallel | same |
  | `torus_bitangent_plane` (Villarceau circles) | `+-(0, 1.5, sqrt(3)/2)` | parallel | same |
  | `viviani_sphere_cylinder` | `(2, 0, 0)` | parallel | same (the figure eight passes its double point twice) |
  | `saddle_xy_plane` (`z = x y` on `z = 0`) | origin | parallel | split (four ends) as given; two curves pass through with the operands swapped |
  | `cone_plane_apex` (double cone, plane through the apex) | apex | cone has none | the two lines end at the apex (four ends) |
  | `cones_shared_apex` (two double cones) | apex | neither has one | the four lines pass through the apex |

  So the defect is confirmed beyond the ssi corpus cases: wherever the normals are
  parallel at a node, the marcher steps over it (the saddle is split in one operand order
  only; why was not traced).
- **The fork's fix** is SolidCraft's tangency-labs/smlib@3c9d1d3 (`replay/1.0.1`:
  1ff92647e; the function is unchanged on that branch since). After the tracing,
  `DoIntersection` calls `sm_SplitCurvesAtSingularPoints` on its final curve list
  (after `AddTangentNodes`, not for a self-intersection). For each `CROSSING` or
  `NEAR_TANGENT` B-spline curve that is not degenerate, it samples
  `s(t) = |n1 x n2|` (`EvaluateNormal`, unit normals) 8 times per unique knot span of the
  3D curve. Each local minimum is refined by golden section on its two neighbouring
  samples (at most 80 steps, down to `SM_EFF_ZERO * (1 + |a|)`). The curve is split at
  the midpoint `t*` when `s(t*) < 1e-6`, `s` rises above 1e-3 both somewhere since the
  previous split (or the curve start) and somewhere up to the curve end, and the point is
  more than `10 * m_dThisApproxTol3d` from both curve ends. The 3D and both UV curves are
  split with `CreateBySplitCurve` from the last split down, and each tail is appended to
  the end of the lists with the curve's type and deviation.
- **Fork-patched C++ on the corpus** (a local build: the pinned `SmSurfaceIntersector.cpp`
  with 1ff92647e applied, swapped into `libSMLib.a`, run through difftest's `sx_*`
  cases): the 9 cases with a parallel-normal node passed through upstream (both
  Steinmetz pairs and the torus in both orders, Viviani in both orders, the swapped
  saddle) change; the 5 others (the saddle as given and the apex cases) are as upstream.
  The fork splits every parallel-normal node into four curve ends. A line-by-line Rust
  port of 3c9d1d3 (a local prototype) matches that C++ bit for bit on all 14 cases
  (`max_rel_err=0`) and on the ssi corpus' `c4*` cases (6 identical, all PASS, where
  upstream has `c4_viviani_figure8` and `c4n_viviani_nurbs` WRONG TOPOLOGY).
- **A side effect of the fork's split:** in `sx_steinmetz_perpendicular` one head piece
  ends up to 2.5e-5 off both surfaces near the node (25 times the approximation
  tolerance; the unsplit curve is within 1e-6). The C++ and the prototype agree on it.
  The cause is not traced. #1304's split, which trims copies with a zero zone
  tolerance, keeps every piece within 1e-6.
- **Rust:** not on `main`. #1304 (open) adds an opt-in
  `ssi_split_curves_at_singular_points` written from the issue's description, not from
  3c9d1d3; see its review for how it differs.
- tangency-labs/solidcore#393.
