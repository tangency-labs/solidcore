# Topology-query changes from projected extrema (#2840)

Unit: `intersect.topology_solver_p2`. Candidate `09495e082`; base
`3f82278a14c1a6104db8e090a59b4bcf677ace01`. The base has zero failures; the
full sweep finds 43 changed cases. These remain native disagreements.

**Accuracy limitation: `analytic_face_sphere_band_5_15` returns distance 0.8605419808705617 against the exact face maximum 0.860555127546399. The gap 1.314667583729e-5 exceeds the requested tolerance 1e-5. This is not certified as a corrected answer at the requested tolerance.**

These 15 inputs invoke a local face solver directly, with projected maximum (operation 15), Single request, and unit projection normal Z. Entry 5 uses vertex P=(0.2,0.3,2); entries 6/7 use curve/edge C(t)=(-1+2t,0,0), t in [0,1]. The returned parameters are (u,v) for entry 5 and (t,u,v) for entries 6/7. They are surface/NURBS parameters, not necessarily analytic angles. All inputs use query tolerance 1e-5.

## Scope and construction

The runner (`crates/difftest/src/units/intersect_topology_solver_p2.rs`, `static_rust`) builds actual analytic faces on their natural surface domains, then directly invokes VertexFace, CurveFace or EdgeFace. In `intersect/topology_solver/p1.rs::solve_point_face` and `p2.rs::solve_curve_face`, projected operations measure trimmed silhouette curves. The point-face implementation explicitly states that bounding edges have already been measured by the surrounding topology traversal; this direct fixture does not perform that traversal. Thus an unrestricted full-face optimum is an independent geometric reference, not proof that this local entry must return the complete topology-query answer.

`p1.rs::face_silhouettes` constructs silhouettes with approximation tolerance 10 times the query tolerance, here 1e-4. Distances returned by these paths are distances to the approximated curves. This is a diagnostic explanation for the sphere/torus accuracy gaps, not permission to loosen the query tolerance or label an approximate answer exact.

## Exact geometric references

- Cylinder axis X: x in [-0.5,0.5], radius 0.5. Both full and half-cylinder project onto x in [-0.5,0.5], y in [-0.5,0.5]. Point maximum is sqrt(0.7^2+0.8^2)=1.063014581273465. Unrestricted line-face maximum is sqrt(1.5^2+0.5^2)=1.5811388300841898.
- Cone axis X: height h in [0.5,1.5], x=h-1, radius h*tan(30 degrees). The farthest projected point is at h=1.5, y=-1.5*tan(30 degrees), giving sqrt(0.3^2+(0.3+1.5*tan(30 degrees))^2)=1.2039996853283075.
- Sphere radius 0.5: the full-longitude band includes latitude zero, so maximum is sqrt(0.2^2+0.3^2)+0.5=0.860555127546399. The half-longitude patch has projected y>=0; its maximum is at (-0.5,0), distance sqrt(0.7^2+0.3^2)=0.7615773105863908.
- Torus axis X: major radius 1, minor radius 0.25. Full and half-major domains have the same projected outer boundary. For P, the farther circle center is (0,-1), so the maximum is sqrt(0.2^2+1.3^2)+0.25=1.5652946437965907. Unrestricted line-face maximum is sqrt(2)+0.25=1.6642135623730951.

## Every changed analytic fixture

Values below are the fresh complete native/Rust output arrays, not reconstructed from mismatch fields. Parameter tuples are given in full, including unchanged components.

| Case | Native distance; parameters | Rust distance; parameters | Independent classification |
|---|---|---|---|
| `analytic_face_cyl_x_5_15` | 0.80000000000000004; (1.5707963267948966, 0.69999999999999996) | 1.063014581273465; (1.5707963267948966, 0) | Matches exact point-face maximum; native value is suboptimal even along the silhouette. |
| `analytic_face_cyl_x_6_15` | 0.5; (0.5, 0, 0.5) | 0.5; (0.25, 0, 0) | Both distances remain 0.5 at valid pairs; changed representatives. Below unrestricted face maximum 1.5811388300841898; no new scalar regression and no C++ defect established. |
| `analytic_face_cyl_x_7_15` | 0.5; (0.5, 0, 0.5) | 0.5; (0.25, 0, 0) | Both distances remain 0.5 at valid pairs; changed representatives. Below unrestricted face maximum 1.5811388300841898; no new scalar regression and no C++ defect established. |
| `analytic_face_halfcyl_x_5_15` | 0.80000000000000004; (1.5707963267948966, 0.69999999999999996) | 1.063014581273465; (1.5707963267948966, 0) | Matches exact point-face maximum; native value is suboptimal even along the silhouette. |
| `analytic_face_halfcyl_x_6_15` | 0.5; (0.5, 0, 0.5) | 0.5; (0.25, 0, 0) | Both distances remain 0.5 at valid pairs; changed representatives. Below unrestricted face maximum 1.5811388300841898; no new scalar regression and no C++ defect established. |
| `analytic_face_halfcyl_x_7_15` | 0.5; (0.5, 0, 0.5) | 0.5; (0.25, 0, 0) | Both distances remain 0.5 at valid pairs; changed representatives. Below unrestricted face maximum 1.5811388300841898; no new scalar regression and no C++ defect established. |
| `analytic_face_cone_x_5_15` | 0.85980762113533138; (0.49999999999999989, 0.27009618943233427) | 1.2039996853283075; (0.5, 1) | Matches exact point-face maximum; native value is suboptimal. |
| `analytic_face_sphere_band_5_15` | empty | 0.86054198087056166; (0.65499056183776805, 0.50000000000000011) | Reported-distance gap 1.314667583729e-5 exceeds query tolerance; approximate-silhouette output, not a strict corrected-result claim. |
| `analytic_face_sphere_quarter_5_15` | empty | 0.76157731058639078; (1, 0.50000000000000011) | Matches exact point-face maximum; native empty. |
| `analytic_face_torus_x_5_15` | 1.3756803885793853; (0.49997295102453776, 0.75) | 1.5652926944322216; (0.5, 0.97356185800456529) | Reported-distance gap 1.949364369036e-6 is within query tolerance, but exceeds 1e-9 relative; approximate-silhouette output. |
| `analytic_face_torus_x_6_15` | 1.25; (0.50000000000000011, 0, 4.1645158291059437e-15) | 1.25; (0.5, 0, 0) | Both distances remain 1.25 at valid pairs; only rounding-sized parameter movement. Below unrestricted face maximum 1.6642135623730951; no new scalar regression and no C++ defect established. |
| `analytic_face_torus_x_7_15` | 1.25; (0.50000000000000011, 0, 4.1645158291059437e-15) | 1.25; (0.5, 0, 0) | Both distances remain 1.25 at valid pairs; only rounding-sized parameter movement. Below unrestricted face maximum 1.6642135623730951; no new scalar regression and no C++ defect established. |
| `analytic_face_torus_half_x_5_15` | empty | 1.5652946437965907; (1, 0.97379278043265505) | Matches exact point-face maximum; native empty. |
| `analytic_face_torus_half_x_6_15` | 1.25; (0.50000000000000011, 0, 4.1645158291059437e-15) | 1.25; (0.5, 0, 0) | Both distances remain 1.25 at valid pairs; only rounding-sized parameter movement. Below unrestricted face maximum 1.6642135623730951; no new scalar regression and no C++ defect established. |
| `analytic_face_torus_half_x_7_15` | 1.25; (0.50000000000000011, 0, 4.1645158291059437e-15) | 1.25; (0.5, 0, 0) | Both distances remain 1.25 at valid pairs; only rounding-sized parameter movement. Below unrestricted face maximum 1.6642135623730951; no new scalar regression and no C++ defect established. |

The eight entry-6/7 rows cannot be classified as eight repaired native extrema. Their changed representatives do not improve the returned objective. Their unchanged distances also do not equal the unrestricted full-face reference. Establishing the complete API behavior would require the enclosing topology traversal, including its boundary-edge contributions; this audit does not claim that traversal regressed or fails.

The closed-form maxima above use original analytic face definitions independently of kernel evaluators. Exact-match claims refer to the reported objective versus the closed form; this audit does not independently reconstruct every surface-parameter map or certify silhouette approximation bounds.


## Remaining 28 point/line fixtures

Entry 1 queries B-rep vertices against `C(t)=(-1+2t,0,0)`, `t in [0,1]`,
projected along unit Z. For vertex `(x,y,z)`, the independent distance is
`sqrt((x-(-1+2t))^2+y^2)`. Enumerate the two endpoints and the clamped foot
`t=clamp((x+1)/2,0,1)` for every vertex to obtain the exact extremum.

Every returned scalar meets strict `1e-9` relative accuracy without an absolute
floor, and every returned candidate lies within the original query tolerance
without an added audit allowance. Four maximum queries replace a stationary
minimum with the correct endpoint maximum. The 24 subdivided minimum queries
retain every native record and the same correct best value; their extra records
are valid within-tolerance endpoints. They are not 24 wrong native optima.
All these records are SingleValue solutions. This checks returned-result validity,
not complete enumeration through the separate topology pruning code.

| Case | Native count | Rust count | Native best | Rust / oracle best | Added only endpoints |
|---|---:|---:|---:|---:|---|
| `static_operation_matrix_1_15_0_0` | 1 | 1 | 0.10000000000000001 | 1.7 | yes |
| `static_operation_matrix_1_15_0_1` | 1 | 1 | 0.10000000000000001 | 1.7 | yes |
| `static_operation_matrix_1_15_1_0` | 1 | 1 | 0.10000000000000001 | 1.7 | yes |
| `static_operation_matrix_1_15_1_1` | 1 | 1 | 0.10000000000000001 | 1.7 | yes |
| `subdivided_1_14_1_12_11_2_1` | 5 | 12 | 0.5 | 0.5 | yes |
| `subdivided_1_14_1_12_11_2_3` | 5 | 12 | 0.5 | 0.5 | yes |
| `subdivided_1_14_1_12_11_50_1` | 5 | 12 | 0.5 | 0.5 | yes |
| `subdivided_1_14_1_12_11_50_3` | 5 | 12 | 0.5 | 0.5 | yes |
| `subdivided_1_14_1_40_30_2_1` | 16 | 40 | 0.5 | 0.5 | yes |
| `subdivided_1_14_1_40_30_2_3` | 16 | 40 | 0.5 | 0.5 | yes |
| `subdivided_1_14_1_40_30_50_1` | 16 | 40 | 0.5 | 0.5 | yes |
| `subdivided_1_14_1_40_30_50_3` | 16 | 40 | 0.5 | 0.5 | yes |
| `subdivided_1_14_1000_12_11_2_1` | 1 | 3 | 500 | 500 | yes |
| `subdivided_1_14_1000_12_11_2_3` | 1 | 3 | 500 | 500 | yes |
| `subdivided_1_14_1000_12_11_50_1` | 1 | 3 | 500 | 500 | yes |
| `subdivided_1_14_1000_12_11_50_3` | 1 | 3 | 500 | 500 | yes |
| `subdivided_1_14_1000_40_30_2_1` | 4 | 12 | 500 | 500 | yes |
| `subdivided_1_14_1000_40_30_2_3` | 4 | 12 | 500 | 500 | yes |
| `subdivided_1_14_1000_40_30_50_1` | 4 | 12 | 500 | 500 | yes |
| `subdivided_1_14_1000_40_30_50_3` | 4 | 12 | 500 | 500 | yes |
| `subdivided_1_14_1e+06_12_11_2_1` | 1 | 3 | 500000 | 500000 | yes |
| `subdivided_1_14_1e+06_12_11_2_3` | 1 | 3 | 500000 | 500000 | yes |
| `subdivided_1_14_1e+06_12_11_50_1` | 1 | 3 | 500000 | 500000 | yes |
| `subdivided_1_14_1e+06_12_11_50_3` | 1 | 3 | 500000 | 500000 | yes |
| `subdivided_1_14_1e+06_40_30_2_1` | 4 | 12 | 500000 | 500000 | yes |
| `subdivided_1_14_1e+06_40_30_2_3` | 4 | 12 | 500000 | 500000 | yes |
| `subdivided_1_14_1e+06_40_30_50_1` | 4 | 12 | 500000 | 500000 | yes |
| `subdivided_1_14_1e+06_40_30_50_3` | 4 | 12 | 500000 | 500000 | yes |
