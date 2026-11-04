# General: the twisted disc's side rounds miss their radius by up to 1.8e-3, and its top round leaves vertex gaps (2026-11-04)

Unit: the generic blend (hub #2243); `blend::fillet_exact`, the fillet trace's surface. Issue
#2652 (solidcraft#11848), PR #2659.

Inputs: the twisted slotted disc of solidcraft#11848 rebuilt from the issue's script (2 slots,
radius 30, height 15, 30 degree twist; #2649's `tests/issue_2648_twisted_disc_top_round.rs`
builder), its 4 twisted side edges rounded 2 mm, then its 12 top outline edges rounded.
SolidCraft's fix list (`tests/kernel_cases/fixes.json`). Rust only (the fixes are not in
the pinned C++).

## 1. Two of the four side rounds are not 2 mm rounds of their walls

For a point x of a side round with unit normal n, the ball's centre x + 2n should be 2 from
both walls the round joins. On side round `FaceId(28)` (and its mirror image `FaceId(30)`)
the error against wall `FaceId(12)` grows to 1.8e-3 at u = 23.6 (of its domain
[8.82, 31.52]), v = 0.9, z = 10.4, and is 5e-4 near the top (u = 27, z = 12.9); against the
other wall it stays below 1e-5. Elsewhere along the round (u = 11..18) both errors are
below 5e-5. Its rail-to-rail section has 3 control points, 106 along u. The other two side
rounds (`FaceId(29)`, `FaceId(31)`) put the top round's ball 2 from their walls to 1.4e-6.

Effect at the top round R2: over each side round's top edge the ball is as wide as the side
round, so its centre should be fixed (the spine crossing z = 13). With #2659 the two
accurate side rounds' top edges are the ball's piece (their end centres agree to 1e-12).
Over `FaceId(28)`'s top edge the centres at the edge's ends are 1.1e-3 apart (solver
tolerance 1.6e-4; the centre of a collapsed offset is ill-conditioned, so the round's 1e-4
error moves it 10 times as far), so it is no ball's piece and is traced through a nearly
fixed centre: the rail folds and the round fails with 1001 as before (#2652's diagnosis).

Suspected cause: the side round's surface between its traced sections (the error is small
at the section rows and large between them, with alternating sign), so its section spacing
or fit tolerance along the twist. Not yet compared with the fork or the C++.

## 2. The top round leaves vertex gaps above the vertex tolerance

The top outline rounded 1 mm builds a manifold solid, but `check_brep` Level2 reports
`SmVertex::Test #1 Bad Vertex Gap` (offset 2.4957061e-05, tolerance 2e-05). At the bands'
fit, 1.101826756514697 (the size named by the BandsOverlap refusal at 1.5 mm), it builds
with the same test failing at offsets up to 5.5762783e-05. With a ball smaller than the side
rounds no offset collapses, so these gaps are not #2652's event; whether they come from
the side rounds of part 1 is not yet known.

So the issue's 2 mm round is not yet refused with that fit: at 2 mm the trace over
`FaceId(28)`'s top edge fails (part 1) before the band check runs, and the fit size does not
build cleanly (part 2).
