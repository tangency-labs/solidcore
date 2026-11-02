#### Updates: WR-79 (2026-11-02). Part 1's circles inside the tolerance band break booleans (#2412)

- **What happens:** with `sphere_small_circle_not_touching` on alone, a 10^3 box and a
  radius-3 sphere crossing its top face by `gap` give a wrong body: box − sphere has
  volume 0 at gap 1e-6, 1e-7, 1e-8 and 1e-9 (open at 1e-7, check failure `SmTrimmingTools`
  test 1 at 1e-9), and box + sphere at 1e-9 is open with volume 1000. Flag off, all are
  closed and valid with the closed-form volumes (1000 ∓ 36π = 886.9026645 and
  1113.0973355). Test: `box_and_sphere_poking_a_hair_through_a_face_is_a_closed_solid` in
  `crates/solidcore/tests/issue_396_small_loops.rs`. Two spheres of radius 3 overlapping
  by 1e-6 to 1e-8 also go from a closed manifold union to a non-manifold one.
- **Where it goes wrong:** the merge's SSI tolerance here is 5e-6 and the edge
  tolerances 1e-5. The circle comes back as a crossing (radius `sqrt(6 gap)`, 7.7e-5 to
  2.4e-3). At the circle `SmMerge`'s common-edge classification calls
  `SmEdge::RadialSectorClassify`. The faces cross there at `θ ≈ sqrt(2 gap / R)`, 0.002
  to 0.05 degrees. The step-off distance grows until it leaves the plane's small disc
  inside the circle. Then the classifier takes any pair under 0.5 degrees as coincident.
  So the box's top face (all of it outside the circle) is paired as coincident with the
  sphere's lower face, and the disc with the cap. The cap's seam is left in the disc as
  a slit, and step 7 gives every region a void flag. At gap 1e-5 only the cap/disc pair
  is coincident, and the result is right.
- **Cause:** the sphere's cap beyond the plane is `gap` high. When `gap` is within the
  tolerance, the cap lies on the plane within tolerance everywhere. The faces part by
  the tolerance only at about `tol / θ = tol sqrt(R / (2 gap))` from the circle. That is
  inside the disc of radius `sqrt(2 R gap)` only when `gap` is above about `tol / 2`. So
  a boolean at that tolerance can tell the circle from a tangent contact only when
  `gap > tol`, which is upstream's touching test (`dDistToPlane > m_dRadius - dTol`).
  WR-79's test (the circle's radius against the tolerance, the fork's
  tangency-labs/smlib@1f5c0f7) differs from upstream's only in the band
  `tol²/(2R) < gap <= tol`. In that band it returns a circle that is exact as geometry,
  but at the B-rep's tolerance it is a coincident patch, not a crossing. Every case in
  #2412 is in that band. The sphere/sphere test is the same in the band.
- **Other flags:** with `degenerate_proto_edge_in_one_vertex` (WR-306) on as well, as in
  SolidCraft's preset, the circle (length 5e-4 to 1.5e-2, under WR-306's 0.02) is
  absorbed into one vertex. The difference then loses the sphere: 6 faces, volume 1000,
  at gap 1e-6 to 1e-8. With every flag on it is the same.
- **Not checked:** the fork's C++ on this case (no fork build here). The fork's later
  c340962 (a sliver thinner than the tolerance removed after a boolean) acts after this
  classification, so it would not undo the coincident pairing. That is from reading the
  commit only.
- **Status:** not fixed. Patching the radial-sector classifier or the merge for this
  band would be a per-case branch for an SSI result the boolean's tolerance cannot
  represent. The choice is part 1's rule: drop it (upstream's gap test), or keep the
  circle for the SSI as a query and type it as a tangent contact in the band. Asked on
  #2412.
