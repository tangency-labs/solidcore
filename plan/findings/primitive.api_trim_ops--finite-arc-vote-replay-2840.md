# Finite-arc trim vote replays (#2840)

Addendum to [the original trim audit](primitive.api_trim_ops--projected-endpoint-callers-2840.md),
using the same rebased candidate and [standalone diagnostic](diagnostics/trim-votes-2840.rs)
as [the finite-segment replay](primitive.api_trim_ops--finite-segment-vote-replay-2840.md).
The embedded inputs reproduce the generated boxes, finite arc constructors,
angles, original projection directions and reference points exactly; arcs are not
replaced with full circles.

Both diagnostic processes and both independent original public `KeepPoint`
operations complete successfully. Replayed deletion predictions agree with their
actual outputs.

| Case | Split faces | Actual retained face IDs | Reachable faces / edges / vertices |
|---|---:|---|---:|
| `gen0070` | 7 | none | 0 / 0 / 0 |
| `gen0086` | 7 | 0, 2, 4 | 3 / 11 / 9 |

Both references classify as side 1. For `gen0070`, faces 0–3 and 5 have votes
`(side0,side1)=(4,0)`; faces 4 and 6 have `(2,0)`. There are no ties, so every
face is deleted by the existing rule without invoking its fallback.

For `gen0086`, votes by face are:

| Face | Side-0 / side-1 votes | Decision |
|---|---:|---|
| 0 | 2 / 2 | Keep after fallback |
| 1 | 3 / 1 | Delete |
| 2 | 0 / 2 | Keep |
| 3 | 3 / 1 | Delete |
| 4 | 1 / 3 | Keep |
| 5 | 4 / 0 | Delete |
| 6 | 2 / 0 | Delete |

The one fallback has an empty carried-edge list and succeeds on bottom face 0 at
`(6.931142857722039, 0.0592057467305267, -0.13803199423746282)`, with
reference-line distance `2.7217452008024307e-15`. Its witness is side 1.
Each retained face's recorded loop, edgeuse, vertex coordinates and edge intervals
match its pre-deletion counterpart exactly.

Each case has eight ignored vertex-classification occurrences, all status `1001`.
Neither uses midpoint votes. Reference-side tolerance is `1e-5`, vertex/witness
side tolerance `1e-4`, and the projected-minimum objective tolerance `1e-6`.

These results explain the output changes through the current whole-face policy.
They do not certify that this policy is the intended interpretation for an open
finite cutter. In particular, both reference points lie on original top
rectangles: whether such queries should retain a reference-containing region
instead is a contract question, not answered by matching the voting code.

As with the segment replay, public `Split` binds UV owners before returning;
original classification occurs earlier. The independent original-operation
outcomes cross-check the replay but do not turn it into an internal execution
trace. No production code, tolerance, unit status or mismatch classification is
changed. Neither case is declared native-correct, candidate-correct in all
geometric senses, or eligible for parking by this evidence.
