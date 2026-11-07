# Exact affine joins in projected point extrema

Review [4501494576](https://github.com/tangency-labs/solidcore/pull/2968#discussion_r4501494576)
identified an inherited #2840 defect: the cache's G1 classification permits a
0.5-degree kink, but the affine solver treated that classification as proof that
the join could not be an endpoint extremum.

![Actual kink and exact objective profiles](diagnostics/projected_g1_join_2840.png)

The raw degree-one polyline has controls `(-1, 1/1024, 0)`, `(0, 0, 0)`,
`(1, 1/1024, 0)` and knots `[0, 0, 1, 2, 2]`. Its turn is about 0.112 degrees,
so the cache marks the middle knot G1. Projection has normal `(0, 0, 1)`.
For query `(0, -1, 0)` the unique minimum is 1 at parameter 1; the pre-fix
public API returns an empty successful result. For query `(0, 1024, 0)` the
unique maximum is 1024 at parameter 1; the pre-fix API instead returns
1023.9995117190992 at parameter 0. Both failures exceed retention slack `1e-8`.
An exactly collinear control retains its two outer All maxima and excludes the
artificial middle endpoint even with retention slack 10.

Writing `u = |x|` and `e = 1/1024`, the squared minimum objective is
`1 + 2 e u + (1 + e^2) u^2`, strictly increasing for `0 <= u <= 1`.
The squared maximum objective minus `1024^2` is
`u ((1 + e^2) u - 2)`, strictly negative for `0 < u <= 1`.
Thus both expected extrema follow from the exact defining dyadic data,
independently of solver samples.

## Construction

Shewchuk's [Adaptive Precision Floating-Point Arithmetic and Fast Robust
Geometric Predicates](https://www.cs.cmu.edu/~quake/robust.html) explains why
near-zero determinant decisions need exact or certified arithmetic. Here the
existing bounded dyadic product accumulator is extended from four to at most
eight terms; no external implementation is imported. Eight products need at
most 4199 bits, within the existing 4224-bit storage.

An approximate G1 boundary is omitted only if its defining adjacent affine
spans have exactly equal projective endpoint positions, exactly parallel raw
chord differences, and matching nonzero orientation. The cross product is
expanded into eight binary64 products before any subtraction rounds the
chords. Constant positive homogeneous weights cancel from the direction test;
endpoint equality uses exact cross-multiplied weights. Repeated knots and
closed seams must pass the same position test. A subdivision strictly within
one affine span is already smooth. Both the early minimum omission and the
endpoint candidate list use this same criterion.

This changes no cache continuity thresholds, solver accuracy allowance, or
retention tolerance. Failure to establish exact smoothness keeps the candidate.
Non-affine stationary solving is unchanged.

## Validation

The standalone public-API reproduction is
`diagnostics/projected_g1_joins_2840.rs`; permanent regressions are consolidated
in `tests/projected_curve_queries.rs`. Internal predicate coverage includes
subnormal/large coordinates, reversed chords, different homogeneous weights,
a true positional gap, and a kink erased by rounded chord subtraction.
All 15 `curve::curve::p1::tests` pass, including the new exact predicate and
the existing adjacent-float interior-minimum refusal. All 24 projected-query
integration tests pass on the final binary span-lookup version. Strict
`cargo clippy -p solidcore --lib --tests -- -D warnings` and modified-file
rustfmt checks pass;
repository-wide and golden gates belong to the parent
integration rather than this isolated review branch.

## Follow-up: cubic joins and retention slack

The same angular-cache omission also affected native cubic C0 knots. The new
public regression uses two cubic Bezier spans with controls
`(x, |x|/1024, 0)` for integer `x = -3..3` and a triple middle knot. It is the
same small kink drawn above, scaled by three. Its projected minimum from
`(0,-1,0)` is 1, and its projected maximum from `(0,4096,0)` is 4096, both at
the middle knot. Before the follow-up the maximum is wrong; the minimum passes.
Native knot multiplicity now proves C1 continuity where applicable. At C0 knots
the defining endpoint chords must pass the exact position/direction predicate;
the existing policies for other geometry types are unchanged.

The [maintainer contract](https://github.com/tangency-labs/solidcore/pull/2840#discussion_r4496188362)
separates retention slack from construction accuracy. A new mutation regression
confirms that the raw `1.5e-12` span's unrepresentable minimum wrongly succeeds
when slack rises to `1e-4`. The affine certificate now compares the whole-objective
bound against the selected candidate's own enclosure: global upper <= candidate
upper for a maximum, global lower >= candidate lower for a minimum. It cannot
borrow slack or another point's enclosure width. A scalar outside its independently
computed candidate enclosure is refused. Retention filtering remains unchanged.

The carrier lower bound is computed directly as
`|N|^2 |(A-Q) dot (N cross D)| / |N cross D|`, using the defining anchor and raw
direction factors. This avoids a dependency-expanded bound at a rounded foot.
The mathematical endpoint difference is also enclosed from those factors before
world-space endpoint subtraction can erase a short direction. A whole leaf may
be skipped only when its objective lower bound exceeds an independently enclosed
incumbent's upper bound plus retention slack. This is a proof that none of its
answers can be retained, not permission to publish an inaccurate local minimum.

Carrier arithmetic steps the result of each correctly rounded binary64 operation
one representable value outward. Exact zero/one identities avoid needless
dependency expansion; nonfinite operands/results and a denominator containing zero
give an unbounded interval and cannot certify success. This includes subnormal
and signed-zero boundaries without the older interval helper's fixed absolute
padding. Only the objective proof is tightened: candidate point enclosures are
unchanged. Unit checks cover a nonrepresentable rational quotient, half a minimum
subnormal, signed zero, overflow, division by zero, cancellation, and a 3-4-5 norm.
When the constrained foot is outside the segment, the analytical endpoint bound
is used directly instead of expanding an interpolation at fraction zero or one.

The exact [representability oracle](diagnostics/projected_affine_representability.py)
uses rational input bits and the quadratic
`5(s-(2+y)/5)^2+(1-2y)^2/5`. The two binary64 parameters bracketing its exact foot
include the best representable parameter. The original raw spans remain covered:

| Knot width | Parameter scale | Query y | Least representable excess |
|---|---|---|---|
| `1.5e-12` | 1 | 0.2 / 0.3 | `8.167385e-9` / `2.756492e-8` |
| `1.9e-12` | 1 | 0.2 / 0.3 | `8.143492e-10` / `9.237773e-9` |
| `1.5e-12` | `1e6` | 0.2 / 0.3 | `8.980330e-11` / `1.091110e-8` |

These exceed rounding at unit geometric scale and now require explicit refusal,
including with nonunit normals and unchanged positive retention slack. The
successful width-`1e-8` controls have excesses `5.05e-17` through `1.05e-15`;
their success still requires the production own-enclosure proof. No test tolerance
is increased to accept the old geometric errors. The adjacent-float refusal and
the zero-retention near-zero minimum remain required controls.

The original `gen1172` kink at `t=1/3` is not an extremum: exact rational
evaluation gives positive left and right squared-distance derivative signs,
with values approximately `0.017863361436596496` and `0.01786536879786583`.
The newly added conservative endpoint at distance `0.03560538710098969` must
therefore be removed even though it lies within slack 0.05. A defining-data
one-sided sign check excludes it after proving exact positional continuity;
the true near-collinear kink extrema remain candidates.

Mutation baseline: the library built from `0406d00e9` (same p1 construction as
parent `ed6119602`) was linked with the new public tests using `rustc --test`.
`retention_slack_cannot_make_an_unrepresentable_minimum_accurate` failed at slack
`1e-4`; `cubic_near_collinear_join_retains_projected_extrema` failed on projected
maximum Single (true value 4096). Green commands are
`cargo test -p solidcore --test projected_curve_queries` and
`cargo test -p solidcore --lib curve::curve::p1::tests`; final results are recorded
with the implementation handoff. The source branch merges parent `ed6119602`
at `bee1aa6ab`; no published history is rewritten.

## Dyadic arithmetic bounds and exact direction factors

The masked binary64 exponent is at most 2047, so conversion to `i32` is exact;
finite inputs limit it to 2046. `parts` therefore produces exponents in
`[-1074,971]` and significands below `2^53`. A product fits in 106 bits of `u128`.
`ea + eb + 2148` lies in `[0,4090]`, so its `usize` cast is lossless on wasm32
and native targets. The high product word is shifted by at most 4154. Eight
same-sign products occupy at most 4199 bits, within 66 words (4224 bits).
`shift % 64` has the literal nonzero divisor 64 and gives a shift from 0 through
63; the other shifts are the constant 64 on `u128`. Adding two low words is at
most `2^65-2`; the `u64` casts intentionally retain low words while the high bits
are carried. All accumulator access is checked. These bounds are also next to
the operations in the source. The panic inventory automatically proves literal
nonzero divisors; no unproved panic site is being parked in its table.

The exact factors implement the mathematical `end - start` direction without
first rounding two evaluated endpoints. For a mapped line it is a scalar times
the defining vector; zero scalar factors and exact vector parallelism determine
constant projection. Two public cases require this distinction:

| Regression | Rounded endpoint subtraction loses | Exact-factor behavior |
|---|---|---|
| `mapped_line_parallelism_survives_endpoint_rounding` | True parallelism for the line from `(0.1,0.2,0.3)` to `(1.7,2.9,4.2)` | Returns the constant projected range at zero slack |
| `mapped_line_does_not_erase_a_small_transverse_direction` | The `1e-16` transverse component added to coordinate 1 | Keeps isolated values rather than claiming a constant range |

For native spans the raw difference is expanded directly into exact product
identities. Its constant positive weight and parameter scale do not change
parallelism. The scalar direction is never inferred from a tolerance or samples.
The local next-parameter helpers retain Rust 1.85 support: standard
[`next_up`/`next_down`](https://doc.rust-lang.org/std/primitive.f64.html#method.next_up)
stabilized in Rust 1.86.

Final strict-construction validation (source branch based on `bee1aa6ab`):

- `cargo test -p solidcore --lib curve::curve::p1::tests`: 18 passed, including
  directed arithmetic boundaries and actual right-valued native point enclosures.
- The optimized `curve.curve_p1` golden replay of the strict construction:
  3102 cases passed; changed 0, added 0, removed 0, broken 0. This includes the
  corrected `gen1172` behavior without changing a golden expectation.
- `python plan/findings/diagnostics/projected_affine_representability.py`:
  exact rational representability checks passed for all twelve controls.

The shared internal `affine_projected_point_enclosure` gives callers an enclosure
of the actual native point, not a neighboring span limit. It selects the right
control at an interior linear-spline knot and the actual first/last control at
natural endpoints. Unsupported curved or variable-weight affine-span interiors
return no enclosure. Its regression distinguishes a discontinuous right value
10 from the left limit 1 and also checks a mapped analytic-line point.

### End-to-end endpoint mutation

The checked construction passed all 28 public projected-curve regressions. To
isolate the elbow fix, the same source on base `bee1aa6ab` was mutated only by
changing `if !omit` to `if false && !omit` in the affine endpoint candidate loop
and `if projected_extremum` to `if false && projected_extremum` immediately before
the generic constrained endpoint loop. No metric, tolerance, or elbow fixture
changed. The SHA-256 hashes of `crates/solidcore/src/curve/curve/p1.rs` for this
experiment were:

- baseline/restored: `778772381c17132f939781a402038404916961b796734201d697756efb3bbe2e`;
- endpoint-disabled mutation: `3493b2a98a4cb228288960257c1b153d46bde55c374b3f17d43676d056d0cdbc`.

Both runs used:

```text
cargo test -p solidcore --test wr_660_planar_path_pose small_profile_mitres::small_radius_short_legs -- --exact --nocapture
```

The mutation failed the body-validity assertion for radius `0.001`, leg `0.01`,
and construction tolerance `1e-7`. It reported two vertex-face gaps of
`0.0018284271` against the validity tolerance `2e-5`. Restoring the baseline bytes,
verifying their hash, and rebuilding made the identical test pass (1 passed,
11 filtered). The mutation is not part of the implementation. Subsequent
repository-edition formatting only reordered imports and formatted assertions;
the final source hash is
`1f0abc62cf3b39915c631bcf2fbc6b451dabdfbb2bdcd872fecbcdc51ebb2308`.

`cargo clippy -p solidcore --lib --tests -- -D warnings` passed on the final
implementation and tests. The modified Rust files pass repository-edition
formatting; `git diff --check` passes.
