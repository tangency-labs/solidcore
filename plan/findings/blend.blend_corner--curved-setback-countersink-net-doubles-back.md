# blend.blend_corner: the countersink's curved Setback is certified but a patch's control net doubles back (2026-11-03)

- **Unit:** `blend.blend_corner` (the Setback beside curved faces, #2343, PR #2366; `Kernel::fixes.fillet_blend_corner`).
- **Inputs:** the fork's capi case `setback curved_countersink_corner` (37c9e1d), rebuilt from solidcore calls in `crates/solidcore/tests/issue_2343_curved_setback.rs`, `countersink_corner_is_a_setback_no_more_concave_than_its_cone`.
  - The part: the block x in [-20, 4], y in [-20, 20], z in [0, 20], less a cone's frustum about (4, 0) widening at 45 degrees from radius 2 at z = 14.
  - The corner is V = (4, 8, 20): the top, the side and the concave cone. The top-side edge, the rim and the generator are rounded 1.25, 1.5 and 1. That is a convex corner beside a concave cone, k = -cos 45 / 5.5.
- **Result:** the corner is a Setback. The rule's setbacks are not certified there, so they grow to 1.414 x the rule's: 2.6514, 2.1211, 3.1816. The fill certifies at that scale (ruling (a), containment), and the fillet completes.
  - `check_brep` level 2 then reports: "SmBSplineSurface::Test #5 Bad NonMonotonic CPts: SmBSplineSurface Control Net doubles back on itself".
- **Expected:** a sound body (17 faces, 2 n = 6 patches).
- **Same as the fork:** its open item (b) at 37c9e1d: "grown to 1.414 x the rule's setbacks the fill certifies, but a patch's control net doubles back (SmBSplineSurface test #5 in sound())".
- **Suspected cause:** the certificate bounds the surface (`II >= k G` and its normal one-signed), not its control net. Test #5 checks the net. A patch can be regular and certified while its net doubles back.
  - It is open whether the fill should keep its nets monotone (a constraint on the quadratic programme) or the net test is too strict for a certified patch. That is a design question.
- **State:** open. The test is ignored as known failing, with this file as its reason. It is also slow, because its setbacks grow.
- **The convex cone** (the fork's `curved_cone_corner`, known failing in the fork at 37c9e1d) passes in solidcore: `cone_frustum_corner_is_a_setback_beside_the_cone`, at the rule's setbacks 1.875, 2.25 and 3.
