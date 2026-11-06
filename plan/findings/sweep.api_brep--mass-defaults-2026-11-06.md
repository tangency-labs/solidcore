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

## Box-face curve sweep: not independently checked

`box_face1_curve_path_cap`, in
`fixtures/sweep.api_brep/box-face-curve-sweep.json`, sweeps face 1 of the box at
`(0,0,0)` with size `(2,3,4)` along a fitted curve through
`(0,0,4),(0,0,6),(1,0,8),(3,0,9)`, with sweep options
`[true,false,false,false,false]`.

Only `values.volume[0]` differs: C++ `30.22825961551991`, Rust
`30.2282600710325`, relative difference `1.5069097266320568e-8`. The query asks
for `1e-6` accuracy. No independent value for this constructed body's volume has
been established. A geometry comparison is also needed before attributing the
change solely to the final volume quadrature: the integrator update also changed
curve length, which a constructor may consume. Neither value is certified correct
by its agreement with the other one.
