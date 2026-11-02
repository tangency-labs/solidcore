# blend.fillet_section: finite-smoothness interior preparation (2026-11-02)

Issue: <https://github.com/tangency-labs/solidcore/issues/2356>.
PR: <https://github.com/tangency-labs/solidcore/pull/2384>.

## Scope and reproduction

This is preparatory work on the approved option-A construction, stacked on
#2322 (`1beec5ccc`, sided definition evaluation) and #2274 (`941770b26`,
procedural rounds). It does not implement a corrected production definition.
The reference remains pinned v1.0.1 C++; no numerical tolerance or unit status
is changed. The construction is fork-only, so the new measurements are Rust
regressions rather than a comparison with a corrected C++ implementation.

Run from the worktree in Ubuntu WSL, with `CARGO_BUILD_JOBS=2`:

```sh
nice -n 10 cargo test --profile difftest -p solidcore --test issue_2356_smooth_interior -- --nocapture
nice -n 10 cargo test --profile difftest -p solidcore --lib smooth_interior
```

The minimal support is `z = 0` for `x <= 5` and
`z = 0.008 (x - 5)^3` for `x >= 5`, represented by a cubic B-spline with
one simple interior knot. The other support is `y = 0`. A radius-2 ball
crosses the knot at leader parameter `s = 0.5`. A control replaces the top
by the plane `z = 0`, retaining the same spline representation.

The tests check that the top is C2 and not C3, measure both sides of the
definition at 33 section heights, and check position, tangent and normal
curvature contact at both rails. A separate test builds the generator's
actual net, checks every interior knot line at 33 heights, and checks that
the net reproduces the definition's sided first s derivatives at the
crossing. The small direct-definition fixtures use a placeholder net only
to construct a definition; that placeholder is never used as net evidence.

Known-crease assertions deliberately require a jump greater than `1e-6`.
They are defect reproductions, not loosened acceptance gates: when the
production correction lands, replace them with the issue's `1e-9` gate.
Analytic control, rail and corrected station-jet checks already use `1e-9`.
The normal angle uses `atan2(|N_left cross N_right|, N_left dot N_right)`;
`acos(dot)` cannot resolve the specified near-zero angle reliably.

## Closed-form kernels and a station derivative detail

The Q=2 endpoint Hermite conditions uniquely give

```text
H(u)  = 3/16 + u/2 + 3u^2/8 - u^4/16
H2(u) = 3u/32 + u^2/4 + 3u^3/16 - u^5/32
G     = -eps (1 - |u|)^3 (3 + |u|) / 16,  |u| < 1
G2    = d G / 2
```

`smooth_interior::correction_kernels` evaluates their values and sided first
two d derivatives. Tests independently differentiate H/H2, check the unit
jumps, zero endpoint jets, invalid widths, and shrinking-window scaling.

At `d = 0`, `G2_s = -3 eps/32` on **both** sides. Therefore

```text
S~_s(s_j,t) = (S_s(before) + S_s(after))/2 + 3 eps Delta2/32.
```

Taking only the average, as one sentence of the issue's net staging brief
suggests, omits the second-jump kernel's derivative. This is a consequence
of the approved endpoint conditions, not a different smoothing construction.
It was reported on the issue before integration. The corrected-station test
includes this term and uses jumps of the full equal-area section, including
the thumbweight's derivatives, to cancel all first and second partial jumps.

## What remains before production use

The current certified bound applies to the **uncorrected** definition. It
cannot be attached to a corrected definition or its changed net. Production
wiring must first supply a correction-aware interval error bound. In
particular, the width tends to zero at a rail, so separate bounds on
`Delta` and inverse powers of `eps` lose the cancellation needed there.
The closed-form derivatives and sampled station cancellation do not prove
a uniform bound, regularity, contact compatibility or whole-face G2.

Remaining work, available independently of this preparatory slice:

- Discover/store crossing one-sided series and fix lambda from neighbouring
  stations before certificate refinement. Preserve maps, copies, transforms
  and definition persistence when introducing correction state.
- Verify the normal jump Taylor orders at both rails; give the designed
  refusal for insufficient face smoothness and name the order that builds.
- Enclose correction contributions uniformly in t, including the shrinking
  windows, and update definition evaluation, the net and its bound together.
- Run the twisted-wedge and net-on-definition acceptance cases, rail contact
  and near-rail curvature measures. The minimal fixture here is not the
  complete capi/rig/ladder acceptance matrix.
- Stage 2: quintic Hermite net with triple knots, second derivatives and the
  sixth-derivative remainder certificate; enable its interior G2 gate.
- Inventory every changed refusal/result and re-derived volume for the
  owner before landing a production behavior change.

The separate arc-half seam tracked by #2330/WR-711 is not modified here.
The approved design cites Hartmann's rolling-ball contact construction;
its abstract was read, and this contribution follows #2356's specified
finite-smoothness correction rather than inventing another construction.
