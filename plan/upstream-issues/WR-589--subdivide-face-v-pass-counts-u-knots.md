#### WR-589. SubdivideFaceAtDiscontinuities cuts along as many V knots as there are U knots
- **Severity:** Wrong result
- **Area:** Topology
- **Location:** `source/SMLib/src/SmBrep.cpp:18483`, `SmBrep::SubdivideFaceAtDiscontinuities` (def. `:18223`).
- **What happens:** The V pass of the loop that builds the iso lines takes its count from the U knots: `lArraySize = (ii==0) ? sUKnotsAtSplit.GetSize() : sUKnotsAtSplit.GetSize()`. When a surface has more V knots than U knots, the face is not split at the extra V knots. When it has no U knot, the face is not split in V at all. The face is then left with V discontinuities inside it.
- **Trigger:** `SubdivideFaceAtDiscontinuities` on a face whose surface has fewer U knots than V knots at or below the requested continuity. For example, a degree-1 sheet with one interior U knot and two interior V knots comes out as 4 faces, not 6. No upstream caller reaches this function. The SolidCraft fork calls it after a variable round (tangency-labs/smlib@6eecf04) and fixes the count to `sVKnotsAtSplit.GetSize()`.
- **How found:** review of the loop while porting it for the fork's variable-round split.
- **Confirmed:** by the `topology.brep_p10` cases `subdivide_u1v2_cinf` and `subdivide_u0v1_cinf`, which compare the default port against the C++. The opt-in `Kernel::fixes.subdivide_face_counts_v_knots` counts the V knots (`crates/solidcore/tests/subdivide_face_v_knots.rs`).
