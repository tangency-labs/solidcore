# Finite-segment trim vote replay (#2840)

Addendum to [the trim caller audit](primitive.api_trim_ops--projected-endpoint-callers-2840.md).
Case `gen0084`, seed `12099798707737135318`. The diagnostic was linked against
the rebased candidate library built for the focused validation; production
numerics are unchanged from `65f2f34fd` through `7a71618a6`.

## Executed diagnostic

The [standalone source](diagnostics/trim-votes-2840.rs) reconstructs classification
from public `Split` topology, then independently executes the original public
`KeepPoint` operation with a fresh kernel. Both complete successfully. It
preserves actual loop/vertex multiplicities rather than inferring votes from
surface-domain samples. Its expected uniform face tolerances are checked before
reconstructing the classification tolerance.

Replay predicts retaining faces `{0,1,2,4,7}`. The independent `KeepPoint` call
retains exactly those faces: **5 faces, 16 edges, 12 vertices**. Each retained
face's recorded loop, edgeuse, vertex coordinates and edge intervals match its
pre-deletion counterpart exactly.

| Face | Side-0 / side-1 votes | Replay decision |
|---|---:|---|
| 0, bottom | 2 / 2 | Keep after fallback |
| 1 | 2 / 0 | Keep |
| 2, left | 2 / 2 | Keep after fallback |
| 3 | 0 / 2 | Delete |
| 4 | 4 / 0 | Keep |
| 5 | 0 / 4 | Delete |
| 6 | 0 / 2 | Delete |
| 7 | 2 / 0 | Keep |

The reference is side 0. Classification tolerance is `1e-5`; vertex and
fallback-witness side tests use `1e-4`. The underlying projected-minimum query
uses the cutter's natural interval `[0,1]`, `Single`, and objective tolerance
`1e-6`. No midpoint voting occurs. Both tie fallbacks have empty carried-edge
lists. Two cut vertices on each tied face return status `1001` and contribute no
vote, as the existing ignored-error policy specifies.

The successful fallback witnesses are:

- Bottom: `(4.818870462064307, 1.2655635309733126, -1.4894456134355734)`,
  on face 0, reference-line distance `1.7763568394002505e-15`.
- Left: `(1.918210854095832, 0.7138666010145358, 10.722419741707533)`,
  on boundary edge 27, reference-line distance `1.6827793267523505`.

Both classify as side 0 and agree with independent geometric witnesses below.
Bottom vertices 17/18 and left vertices 17/16 also classify successfully as side
0. Their unconstrained cutter parameters, and that of the left witness, are
negative; these classifications therefore exercise the finite starting endpoint.
The diagnostic records classification status, not the internal selected parameter.

## Independent geometry

Let `n` be normalized projection direction, `v=b-a` the cutter direction, and
`P=I-n*n^T`. The unconstrained cutter parameter is
`((p-a) dot (P*v))/(v dot (P*v))`; signed side is `(v cross n) dot (p-a)`.
The reference signed value is `-38.628256326372636`.

The bottom witness is the reference line's intersection with the bottom plane,
inside the original rectangle; its distance to that line is zero in exact
geometry. Its cutter parameter is `0.005441046124132234`.

The left witness minimizes the convex projected-distance quadratic over the
original left rectangle. On its top edge,
`y-reference.y = n_y*n_x*(xmin-reference.x)/(1-n_y^2)`.
Its y derivative vanishes and its z derivative points beyond the upper boundary.
The cutter parameter is `-0.10056845694727339` and signed side is
`-25.42719516404516`. Both witnesses lie within the fallback reference segment.
These reference-line distances must not be confused with distance to the finite
cutter, approximately `3.736976272070221` and `2.670519282803533` respectively.

## Interpretation and limits

The replay explains the two additional faces through the existing whole-face tie
rule. Their opposite-side corners alone do **not** establish a new deletion
regression. The earlier audit's adverse observation is narrowed by this evidence;
it is not a reason to classify the native result as a defect or park the mismatch.

Public `Split` binds UV owners before returning; original `KeepPoint` classifies
before that binding. This is a public-API replay, not an instrumented execution
trace inside the private classification function. The independent `KeepPoint`
execution and exact retained-topology match cross-check its outcome. Neither
establishes the intended finite-trim contract or resolves `gen0070`/`gen0086`.

Relevant existing logic is in `topology/brep/p10.rs::faces_to_delete`: vertex
votes; midpoint votes only when both counts are zero; equal-vote shape/curve
fallback; and deletion when the resulting side differs from the reference.
No production code, tolerance, status or mismatch classification changes here.

To reproduce, build this revision's library with the `difftest` profile, then
compile the diagnostic using `rustc --edition=2021`, that build's dependency
directory via `-L dependency=...`, and its exact `libsolidcore-*.rlib` via
`--extern solidcore=...`. Run with argument `gen0084` under the repository's
normal machine budget. The same source accepts `gen0070` and `gen0086` for the
separate finite-arc replay.
