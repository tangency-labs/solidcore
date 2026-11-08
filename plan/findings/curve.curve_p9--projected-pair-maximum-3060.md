# Projected curve-pair maxima on the main baseline

Issue #3060 reports successful projected maxima below admissible witnesses.
The four cylinder/half-cylinder curve/edge-face queries return 0.5 instead of
sqrt(2.5). Their projected surface is [-0.5,0.5]^2 and query line is
[-1,1] x {0}. Squared distance is convex separately in each coordinate;
opposite x endpoints with |y|=0.5 attain sqrt(2.5).

![Reported cylinder pair and larger admissible pair, same view.](figures/projected-pair-witness-3060.png)

## Construction and dependencies

This extraction starts at main `7a593124c9`. It carries the existing reviewed
curve-pair construction without the #2840 or #2968 point-query solvers.
The only topology change forwards outer ProjectedMaximize/Single to the
curve-pair Single construction; other requests still ask for All.

`projected_pair` partitions both native domains at knots, evaluates feasible
candidates, and encloses every parameter rectangle. A returned candidate must
contain the global upper bound in its own arithmetic enclosure. Retention slack
never supplies solve accuracy. Work is bounded to 32,768 boxes, depth 128 and
degree 18; cancellation and unresolved proof errors propagate.

The native analytic Line map comes from shared `curve::native_affine` support:
its start/end definition and interval jet preserve the original arithmetic,
including bounded STEP mapping and inside-out direction. No parent point-query
selection or join policy is required. The existing main NURBS interval helpers,
Hermite power coefficients and global accumulator are reused. Ellipse forms
are checked against analytic support by the existing quadric identity machinery.

All is not replaced with Single. One-span affine pairs and identical native
NURBS with certified coordinate monotonicity use finite candidate constructions;
retention and ordered deduplication use the existing accumulator. Other curved
All classes retain their existing path. Uncertainty inside a selected class
refuses rather than falling back.

The enclosure/branch-and-bound construction follows Neumaier, sections 11–12,
[Complete Search in Continuous Global Optimization and Constraint Satisfaction](https://arnold-neumaier.at/glopt/mss/Neu04.pdf).
Interval arithmetic is described in
[Patrikalakis, Maekawa and Cho, section 4.7](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node45.html).

## Golden changes against main

Measured with fresh goldens-profile builds of main `7a593124c` and this branch. Main passes all 283 units. This branch changes 35 `curve.curve_p9` records and 12 `intersect.topology_solver_p2` records. Nothing is added, removed or broken.

**Why each new value is right.** The construction returns a candidate only when the global upper bound over every parameter rectangle fits inside that candidate's own arithmetic enclosure, so each value below is a certified maximum of the native curve pair. Each case was also checked independently, by its own witness parameters, closed form or exhaustive affine-corner enumeration, in the [frozen pair audit](https://github.com/ChrisNourse/solidcore/blob/20398468a467d06b1f5cd92558f1b2f0b818dcb9/plan/findings/curve.curve_p9--frozen-pair-audit-3060.md), which covers all 35 `curve.curve_p9` cases and the analytic-face records. Every record blessed here is identical to the one that audit validated. Where main returned an empty success, a feasible pair at the new value exists, so main's empty answer was wrong.

**Reading the golden file.** A stored array longer than 16 entries is kept as `{"n", "blocks"}`: per 256-entry block, the sum, absolute sum and index-weighted sum of every output entry (counts, types, object ids, values and parameters). The block numbers are therefore not maxima. The values below are decoded from the full output arrays.

### `intersect.topology_solver_p2` analytic face pairs (entries 6/7, tolerance 1e-5)

| Cases | Main | This branch | Closed form |
|---|---|---|---|
| `analytic_face_{cyl_x,halfcyl_x}_{6,7}_15` | 0.5 | 1.5811388300841898 | sqrt(1.5² + 0.5²) = sqrt(2.5), at opposite projected endpoints |
| `analytic_face_cone_x_{6,7}_15` | no answer | 1.7320508075688772 | sqrt(3) |
| `analytic_face_sphere_{band,quarter}_{6,7}_15` | 0.5000394345169934 | 1.5 | 1.5 |
| `analytic_face_torus_half_x_{6,7}_15` | 1.25 | 1.664213562373095 | sqrt(2) + 1/4 (to the last digit) |

**Known failing, not blessed:** `analytic_face_torus_x_{6,7}_15` (main 1.25). For a face, this construction certifies the maximum on the face's traced silhouette curve, not on the surface itself. The full torus's traced silhouette lies slightly outside the torus, so the certified curve maximum 1.6642154866673036 is above the native surface maximum sqrt(2) + 1/4 = 1.6642135623730951, and its returned point measures 1.6642135459403802 on the surface. Both cases are listed in `known_failing.json` under #3109. Native certification against the surface is #3110's work; it will remove the entries.

### `curve.curve_p9`: 35 records, largest reported value

On main, 8 cases return an empty success. 25 report a smaller maximum than the certified one. gen0421 keeps its value with a different answer record, and gen0880 moves by one rounding step.

| Case | Main | This branch |
|---|---|---|
| `gen0011` | empty (0 answers) | 7.674280832269905 |
| `gen0131` | 27.292406300295426 | 38.353130876659016 |
| `gen0144` | 70.88929067495704 | 90.48291006330427 |
| `gen0155` | 37.87672736773793 | 41.69836616313635 |
| `gen0302` | empty (0 answers) | 4.217760138079574 |
| `gen0313` | 3.5885710823184898 | 3.677522280315252 |
| `gen0373` | 6.113465906772089 | 8.284826190867808 |
| `gen0421` | 17.11468532354081 | 17.11468532354081 |
| `gen0445` | 2.1774745915837537 | 10.052741867014833 |
| `gen0743` | 83.90388228273946 | 87.12927044222975 |
| `gen0844` | 0.0 | 10.979584714079621 |
| `gen0880` | 10.50084800781891 | 10.500848007818902 |
| `gen1086` | 17.609097780432535 | 73.99422438322947 |
| `gen1116` | 9.918860532286931 | 10.3994004442039 |
| `gen1458` | 1.1799020931977433 | 1.7704231028846957 |
| `gen1658` | empty (0 answers) | 8.001863852854438 |
| `gen1744` | empty (0 answers) | 9.878574576753216 |
| `gen1839` | 0.0 | 4.408015599603469 |
| `gen1860` | empty (0 answers) | 1.5824476883128602 |
| `gen2016` | 2.2711480262147923 | 3.4856960129819874 |
| `gen2017` | 19.071210555301338 | 21.504888840235413 |
| `gen2050` | 6.137136311338198e-08 | 40.001231972795864 (2 answers) |
| `gen2131` | 65.62255209558627 | 94.33068554634433 |
| `gen2187` | 0.09955759757712898 | 6.310837322899526 |
| `gen2189` | 0.09248840291760957 | 4.5154998991468 |
| `gen2325` | 0.8055899438915088 | 0.9235168507522138 |
| `gen2440` | 2.44769534587924 | 8.615017213875161 |
| `gen2501` | 0.6112768486652884 | 1.2800395290048499 |
| `gen2585` | 0.0001286767978681419 | 0.9962268759220879 (2 answers) |
| `gen2603` | empty (0 answers) | 7.017259757475725 |
| `gen2768` | 6.160254511384329 | 8.199199517002864 |
| `gen2771` | 12.110728525770247 | 18.087676035977736 |
| `gen2962` | 82.52299106097162 | 107.26946892615773 |
| `gen2969` | empty (0 answers) | 55.95810876408308 |
| `gen2976` | 5.583408252286654 | 9.516387163397034 |

## Work deliberately kept separate

The original ellipse All case gen0472 remains open under #3081. The other #3060 All-request cases keep their main values: general All enumeration and continuum answers remain open. This is not a claim that every All maximum is fixed.

Point-face native proposals, residual certificates and vertex-face incumbent ownership belong to #2968. A maximum on a silhouette curve alone does not prove a maximum over the whole native face; that obligation is #3110's.
