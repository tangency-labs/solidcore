# sweep.api_brep: changed volume observations after the integrator update (2026-11-06)

Updates: plan/findings/mass.integrator--whole-revolved-sphere-volume-fails--addendum-2026-11-01-shared-at-1e-8.md

Unit: `sweep.api_brep`. PR: #2895, following the integrator change #2890.

On main `f826c16133a44c5d821d2204adfdccb979e3c65c`, the unchanged native
comparison runs 696 cases: 693 pass and these three fail. The relative comparison
threshold remains `1e-9`, with no absolute floors. These observations are not
passes, and this record does not establish a C++ defect or authorize parking.

## Revolved spheres

The inputs are `wr471_rot_half_circle_whole_sphere_off_0` and
`wr471_rot_half_circle_whole_sphere_off_3e-7` in
`fixtures/sweep.api_brep/wr471-rotational-whole-sphere-volume.json`: a quadratic
radius-5 semicircle, angles -90 to 90 degrees, centered at `(offset,0,0)`,
revolved 360 degrees about Y through the origin, without caps, with default flags.

An independent Rust probe reconstructed the exact fixture inputs using
`primitive::api_curves::create_arc` and `sweep::api_brep::create_rotational_sweep`.
Both resulting surfaces are actually `SurfaceKind::Sphere` with radius 5.
In particular the small offset input must not be treated as an unsnapped spindle
surface of revolution. The represented spheres have volume `4*pi*125/3`, or
approximately `523.5987755982989`.

| Offset | Accuracy requested | C++ | Rust |
|---|---|---|---|
| 0 | 1e-6 | 523.5988331275332 | 523.5987784379859 |
| 3e-7 | 1e-6 | 523.598956772222 | 523.5987784379857 |
| 0 | 1e-8 | status 1001, no volume | 523.5987756272117, status 0 |
| 3e-7 | 1e-8 | status 1001, no volume | 523.598775627212, status 0 |

Rust's relative errors are about `5.42e-9` at requested `1e-6` and `5.53e-11`
at requested `1e-8`. The C++ coarse values also satisfy their requested `1e-6`
accuracy, and its finer query refuses rather than returning a wrong volume.
The zero in its difftest value field is the runner's placeholder for that refusal.
Consequently a disagreement against the stricter comparison threshold is not by
itself proof of a wrong C++ volume. With the default body-length scale, Rust still
refuses at `1e-10` for both fixture inputs.

The old finding's statement that both implementations refuse at `1e-8` no longer
holds after the new step control. The separate library test
`sweep::api_brep::tests::wr471_whole_sphere_volume` reproduces an obsolete refusal
expectation on unchanged main: the default scale now returns
`523.5987756272146` at `1e-8`.

## Box-face curve sweep: polynomial flux checked independently

`box_face1_curve_path_cap`, in
`fixtures/sweep.api_brep/box-face-curve-sweep.json`, sweeps face 1 of the box at
`(0,0,0)` with size `(2,3,4)` along a fitted curve through
`(0,0,4),(0,0,6),(1,0,8),(3,0,9)`, with sweep options
`[true,false,false,false,false]`.

Only `values.volume[0]` differs: C++ `30.22825961551991`, Rust
`30.2282600710325`, relative difference `1.5069097266320568e-8`. The query asks
for `1e-6` accuracy.

An external probe rebuilt this exact input on pre-update main (kernel source tree
`b3c30dab24371664bef6d01651456e6b455b2a5b`) and current main with the test-only
expectation update (kernel source tree `55266757953e46d492b74596207d2e19e12a41e0`).
The two complete `solidcore_convert::brep::to_record` exports are identical:
6 faces, 12 edges and 8 vertices, with no changed serialized field. Thus the
changed Rust observation is not explained by changed constructed geometry.

The exported boundary has four non-rational degree-(3,1) spline sides, each
trimmed to its full natural rectangle, and two planes trimmed by straight
pcurves. Its volume can therefore be integrated as a polynomial, independently
of the kernel's adaptive integrator. An external Python calculation converts
every binary64 coefficient and knot to its exact `fractions.Fraction`, constructs
the Cox-de Boor basis polynomials on each distinct knot span, and integrates the
divergence flux `P dot (P_u cross P_v) / 3` by integrating monomial coefficients.
It uses the recorded shell orientation and checks the side trims are the complete
rectangles. The planes use their signed pcurve polygon area, not the larger face
domain rectangle.

On a side, with local span parameters, write `P(s,t)=A(s)+t D(s)`. Integrating
out `t` leaves `[det(A,A',D) + det(A,D',D)/2] / 3`; the remaining integral is
polynomial and evaluated exactly as rational arithmetic. The independent
implementation also checks a 2-by-3-by-4 box (volume 24) and the polynomial patch
`P(s,t)=(3s,t,s^3)` (flux -1/2) against their closed forms.

The exact rational flux, rounded once for display, is `30.228260071009444` on
both exports. The six signed face contributions, rounded for display, are
`[-3.8642169281794123, 0, 18.59998244584193, 10.076086688348184,
13.416407864998739, -8]`. Rust differs by about `2.31e-11` absolute
(`7.63e-13` relative); the C++ observation differs by about `4.55e-7` absolute
(`1.51e-8` relative). This checks the flux of the actual represented polynomial
boundary; it does not certify the sweep's fitted geometry against its intended
motion. Both observations satisfy the query's requested relative accuracy, so the
strict comparison discrepancy alone still does not justify a C++-defect marker.
