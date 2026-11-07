# Native affine curve/face maximum construction

Issue #3109; draft #3110. This is the face-owned pair boundary left after #3077,
not a replacement for its bounded curve/curve construction or #2968's point/face
certificate. The initial evidence checkpoint was 4f883e2dd; the following native construction
is a work-in-progress checkpoint, pending traversal-domain accounting.

![Current native returned pair and exact reference in the same projected view.](figures/torus-pair-attainment-3109.png)

The two geometric pairs are visually close at this scale, but the reported fitted
scalar lies above the exact global bound. The right panel is an analytic reference,
not an output of an implemented repair.

## Measured defect

Actual public queries at pair checkpoint 69d8abd33 return the full-torus cases
`analytic_face_torus_x_6_15` and `_7_15` with scalar 1.6642154866673036 and native
parameters `(0, .5, .12494065250506982)`. Independent native evaluation gives
1.6642135459403802: an attainment error 1.940726923388425e-6. The half-torus controls
at `(0,0,.12499999813735485)` remeasure as their reported 1.664213562373095.
Original requests and outputs are preserved in
[the audit artifact](diagnostics/torus_pair_native_audit_3109.json); the standalone
[Rust evaluator](diagnostics/torus_pair_native_audit_3109.rs) reconstructs the exact
runner surfaces and checks native-domain membership without padding.

The parent/child checkpoint without #3077 returned 1.25 at `(.5,0,0)` for all four.
That old candidate has the elementary counterexample `(0,0,0)`, whose native
projected distance is sqrt(41)/4. This witness is not the exact maximum.

For the full domain, write projected torus coordinates as
`x=r sin(v), y=(R+r cos(v)) cos(u)`, R=1, r=.25, and line coordinate l in [-1,1].
Then `|x-l| <= 1+r|sin(v)|`, `|y| <= R+r cos(v)`, and the triangle inequality gives
`distance <= ||(1,R)+r(|sin(v)|,cos(v))|| <= sqrt(2)+.25`.
Equality is attained at l=-1, u=0, v=45 degrees. Thus the global bound is
1.6642135623730951; the full-torus reported scalar exceeds it by 1.924294208421884e-6.
The bound applies to both full and half major-angle domains because both contain u=0.

## Root cause and construction

The topology pair path solves against fitted 3D silhouettes, then changes the
solution's parameters and object ownership to the source face while retaining the
fitted scalar. The strict point/face refinement is not applied on this pair path.
A precise maximum of a fitted curve is insufficient evidence about its owner.

[Boyd and Vandenberghe, *Convex Optimization*, sections 3.1 and 3.2.2](https://web.stanford.edu/~boyd/cvxbook/bv_cvxbook.pdf)
establish norm convexity and preservation under affine composition. Applying that
result directly: for fixed surface point S and affine C(t)=(1-t)A+tB, projected
norm `||P(C(t)-S)|| <= (1-t)||P(A-S)||+t||P(B-S)||`, hence no interior t exceeds the
larger endpoint value. Taking the supremum over the face preserves this bound.
The double-cross map P uses the caller's unchanged normal, including its squared
scale. A finite piecewise-affine curve reduces to the union of span endpoints.

Planned implementation reuses the existing native surface residual certificate
for each structurally proved affine span endpoint, with the complete endpoint
**definition enclosure** carried through it. A rounded general endpoint Vec3 is
not silently substituted for that definition. Candidate construction may use
native active-set refinement, but acceptance requires every endpoint/face domain
upper bound to fit the actual returned candidate's own pair enclosure.

Structural selection precedes numerical proof. Selected failures propagate, never
fall back to the fitted scalar. Native requested-domain intersections, knot-side
ownership, discontinuities, non-unit/zero normals, cancellation, work bounds,
external incumbents, face trimming and returned parameter order require tests.
General non-affine curves are outside this initial construction. Existing All
contracts are not silently converted to Single.

## Validation status

Native evaluation and full-query baseline replays are complete. Dedicated full
and half torus kernel regressions ran on the integrated base: the full-torus test
fails with reported 1.6642154866673036 versus native 1.6642135459403802, while the
half-torus control passes (one pass, one expected failure). Production implementation, endpoint-enclosure interface
review, affected-case audit, final before/after rendering and broad gates remain.


## Native local construction checkpoint

The native endpoint residual construction passes six dedicated tests: the original
full and half torus, reversed operand order, trimmed curve interval, normals scaled
by .001 and 1000, and the zero normal. The seven existing residual tests and twelve
projected-bound tests also pass. The full-torus red regression now succeeds without
using caller slack as accuracy or returning a fitted-silhouette scalar.

This is not the complete repair: topology traversal can have skipped or discarded
domains before the local face solve. A previous candidate's attainment enclosure
does not prove coverage of that history. The next checkpoint records original
domain groups before pruning and independently certifies them against the final
retained pair. Unknown mixed domains must refuse. This draft is not ready until
that traversal proof, cross-owner tests and the affected-case audit are complete.


## Traversal-domain ledger (implementation review)

The local construction alone is insufficient: an earlier fitted scalar can prune a
winning native face, and a retained point only proves attainment, not coverage of
previously discarded domains. The ledger records original operand groups before
pruning and checks their full Cartesian product at completion.

- **Capture:** `Ledger::record_trees` reads complete original tree object lists,
  native edge intervals, face domains and supplied curve intervals before either
  global branch or topology object pruning. Direct `LocalSolver::solve_it` records
  its complete local pair. References retain the original definitions; no fitted
  silhouette enters the ledger. Scoped tree solvers move the same ledger back.
- **Membership:** a retained solution must belong to one recorded original pair,
  with both owners matching and every curve parameter in its requested and natural
  domains. Native face-point feasibility also requires the exact rectangle proof.
  A valid native scalar at an out-of-trim point is not an admissible candidate.
- **Coverage:** affine curve/face groups use endpoint convexity and the shared
  homogeneous face residual. Curve/edge groups reuse the complete curve-pair
  construction, reconstructing its returned point enclosure with the identical
  bounded requested intervals. Point pairs use direct enclosures. Unsupported
  mixed domains refuse; no unvisited domain is inferred from an earlier answer.
- **Selection:** a feasible alternative is retained when scalar ordering discards
  the candidate whose own enclosure covers the domains. The certificate always
  belongs to the returned pair. The external incumbent cutoff is outward rounded
  and is checked separately after coverage; caller slack never stops the solve.
- **Work:** at most 4,096 recorded operands and 16,384 domain-pair checks. Original
  groups avoid Cartesian allocation; empty groups are discarded. The existing
  curve and surface proof budgets and cancellation remain in force.
- **Entry ownership:** `solve_trees`, `solve_two_trees`, the shape/B-rep wrappers,
  and the `GlobalSolve::solve_branch` trait entry own full capture/finalization.
  A raw generic free-function traversal or untracked `local_solve_branch` partial
  call cannot enforce finalization and explicitly refuses an activated native
  pair query. A partial branch return is never presented as a global certificate.

The first ledger cycle passed 15 focused native tests, all 60 topology-solver
unit tests and 33 projected-pair tests. These include direct/tree parity,
vertex/edge/face history, an unattained legacy scalar that would prune the winning
face, different face owners, an out-of-trim retained native point, empty groups,
cone sqrt(3), and both sphere controls at 1.5. The subsequent focused cycle
passed all 23 tests, adding trait/raw/partial entry guards, external-incumbent
empty/accepted/ambiguous controls, cancellation, the 128/129-span budget boundary,
and unchanged single-tree/cull behavior for an unrelated minimum query. The
strict library/tests lint check also passes. Fresh `GlobalSolve::solve_it` reuse
needs a separate lifecycle reset before this checkpoint is complete: branch
accumulation must preserve history, while a new request must not inherit old
domains or an active certificate. Dedicated reuse controls are pending.
Affected-unit output/refusal accounting and final broad gates remain pending.
