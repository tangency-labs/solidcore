#### WR-491. `SmSurfaceIntersector::DoIntersection` ends no curve where two intersection branches cross
- **Severity:** Wrong result
- **Area:** Intersection
- **Location:** `SmSurfaceIntersector::DoIntersection` / `TraceCurve` (the marching); `SmAdvSurfaceIntersector::ComputeSingularityPoint`, which only runs on start points, never on points reached while marching.
- **What happens:** Where two branches of an intersection cross (the surface normals are parallel there), the tracer passes through the point without ending its curve. The point set is complete, the topology is not: a curve runs through the crossing, and fewer than four curves end at it. Two equal cylinders with perpendicular axes (radius 1, axes x and z) meet in two ellipses crossing at (0, 1, 0) and (0, -1, 0): SMLib returns four curves from the cylinder seams, all ending at (0, 1, 0), and two of them run through (0, -1, 0).
- **Trigger (ssi-corpus cases):** `c4_steinmetz`, `c4_viviani_figure8`, `c4n_viviani_nurbs`, `y23_fig1a_two_loops_cross`, `y23_t6_ex3_two_crossing_curves`, `y23_t3_ex4_torus_ellipsoid` (verdict `WRONG TOPOLOGY`, "singular point not identified").
- **How found:** ssi-corpus, exact ground truth (reported in tangency-labs/solidcore#393).
- **Confirmed:** yes, by the C++ build (the six cases above give `WRONG TOPOLOGY` on the C++ engine).
- **Rust port:** reached (same output as C++); opt-in `Kernel::fixes.ssi_split_curves_at_singular_points` splits the final curve list at the crossings (see `plan/findings/intersect.surface_intersector_p5--ssi-393-split-at-singular-crossings.md`).
