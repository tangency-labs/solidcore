# blend.fillet_executive_p3: distant rim mistaken for endpoint band overlap (2026-11-06)

Issue: #2903. PR: #2924. Ordinary SolidCraft operation from its kernel update #12413.

The captured 60×40×30 shelled box has wall 3, outer vertical rounds R6, inner rounds R3,
and an inner floor round R1. The final inner-top-rim round/chamfer size is 0.47.
The actual application operand passes full integrity checks. Its symmetric mirror
seeds should build the same closed body (35 faces, 80 edges, 48 vertices), but one
seed returns `NoRoom/BandsOverlap` while its mirror completes. Both main `513e6fcfc`
and candidate `9856ffee5` reproduce with a fresh Core kernel importing the captured
SMB. Clearing application round history does not change it. An ops reconstruction
of the recipe does not reproduce; it has different geometry representations.

## Root cause

`FilletSession::geometry_bands_clear_other_rails` intends to exclude projections
onto source-edge endpoints: they are beside the band rather than across it.
It checks only whether the returned parameter is strictly within the edge interval.
The constrained closest-point solve can return an endpoint one parameter ulp inside
that interval. The following chord test then treats a distant point as inside an
unbounded strip instead of excluding the endpoint projection.

For the first rejection, edge B has interval `[2.6532560294479715, 2.713865182121104]`
and returned parameter `2.7138651821211037`. Its evaluated closest point is exactly
`E=(3,34,30)`, the edge endpoint. The rail endpoint is `R=(2.53,34,30)`. The tested
point is `P=(2.544036467113312,5.688205100494965,30)`, 28.311794899505 mm away across
the face. Its projection along the width-0.47 chord is 0.0140364671133, between the
existing margins, so the operation is falsely refused. The mirror's rounding of the
endpoint parameter changes the classification.

Recognize projected endpoints in geometric space using the existing maximum solver
approximation tolerance, before the chord-width test. No parameter epsilon, new
flag, retry, reordered chain, altered band margin or relaxed output assertion is
introduced. The guard remains a sampled test (15 points on each rail), not a complete
intersection proof; endpoint neighborhoods are treated within the same tolerance.

## Independent expected volumes

Let `A=1−π/4`, `M=5/6−π/4`, and the inner rounded-rectangle perimeter
`L=152+6π`. Integrating the cross section and its first moment against the parallel
curve Jacobian `1±κu`, with total turning `2π`, gives:

- Base: `60*40*30 −54*34*27 −4*A*6²*30 +4*A*3²*27 +L*A −2π*M`.
- Rim round removal: `L*A*r² +2π*M*r³`.
- Rim chamfer removal: `L*d²/2 +π*d³/3`.

At size 0.47 the expected round volume is 21737.746007562786 and chamfer volume
21726.89743733982 mm³. The regression checks each of the four mirror seeds against
these closed forms at relative 1e-9, topology counts, manifold status and Level2
integrity. Each starts from a fresh input before any measurement/check warms it.

Literature read: Patrikalakis, Maekawa and Cho, *Shape Interrogation for Computer
Aided Design and Manufacturing*, §11.2.1 (offset derivative `1+κd`) and §5.1
(point projection with an orthogonality condition for an interior minimum):
https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node219.html and
https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node69.html.
A constrained boundary minimum is distinct from the interior closest-point
condition. Geometric endpoint recognition preserves that distinction independently
of curve parameterization. These sources do not prescribe this specific fix.

Validation results belong in the PR description. The captured input and regression
are `crates/solidcore/tests/data/issue_2903_shelled_box_rim.smb` and
`crates/solidcore/tests/issue_2903_mirrored_shell_rim.rs`.
