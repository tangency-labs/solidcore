# mass.integrator: the precise-property steps' fifth-order error estimate moves 187 cases of 3 units (2026-11-05)

- **Units:** `mass.integrator` (89 new known mismatches), `boolean.api_brep` (51 new; the 12
  quarter-annulus cases of the rotational sweeps re-worded) and `topology.face_p2` (47, a new
  list: the unit goes from `green` to `partial`).
- **Change:** #1556's follow-on. The precise-property integration (face boundary, iso-u lines,
  edges) did not hold the accuracy it was asked for: on the quarter annulus of
  `rotational_sweep_face0` (exact 6.25 pi) -7.7e-2 relative at 1e-2, -4.9e-3 at 1e-3 and
  -5.1e-4 at 1e-4. Two causes:
  1. `SmODEAreaIFEO::Evaluate` skips a point whose `|du/dt| <= SM_EFF_ZERO` and leaves the
     derivative array as it was. The integrator's array still holds the previous step's start
     value, so a step starting at the inner arc's u turn (its middle, where a first step of half
     the arc lands) began from a stale slope (-35 where the integrand is 0). `OdeArea` now
     writes the integrand everywhere: zero where `du/dt` is exactly zero.
  2. The Cash-Karp embedded estimate is the difference of the fourth- and fifth-order values, a
     rule zero on cubics: it measures `g^(4)` only and can vanish where `g^(4)` changes sign
     within the step. On the arc's first step it read 3.1e-5 where the step was off by 3.5e-3.
     For a quadrature (`OdeOptions::quadrature`) the six samples also give the fifth divided
     difference, and `h (-1/960) g[c0..c5]` is the leading error of the fifth-order value
     (0.83 to 1.26 of the true error on the arc's steps). A step is accepted on the larger of
     the two. The unachievable-accuracy guard keeps reading the embedded estimate alone.
- After: the four quarter annuli hold every accuracy from 1e-2 to 1e-8 (worst -3.7e-6 at 1e-4;
  test `precise_area_holds_its_accuracy`). Primitives' areas hold theirs at 1e-2 to 1e-7. The
  default path now holds #1435's repro (cylinder side 1.1e-9 at 1e-8, was 5.2e-8) and torus
  Ixx (3.7e-10, was 1.7e-9). A NURBS sphere's area held nothing before (-1.09e-3 at 1e-3) and
  now holds every accuracy.

## Why the cases move, and who is right

Every moved value is the same integral by other steps. Against a reference integrated at 1e-11
(mass.integrator's face and edge ops, the largest component error over the largest
magnitude): Rust holds the requested accuracy in every case audited; C++ misses it in some,
e.g. `issue1435_multispan_sub_eps0.001_*_flags1` (C++ 1.06e-3, Rust 5.0e-5) and
`issue1435_multispan_full_eps1e-06_*_flags1` (C++ 9.1e-7, Rust 7.5e-8). Two edge cases,
`gen0014` and `review_edge_exact_quarter_circle_1e-09`, are at 1.2x the accuracy on both
sides alike (the 1e-11 reference's own error at that level).

## What is left

- `issue1435_multispan_full_eps1e-08_orientation{0,1}_flags63`: Rust fails (SM_ERR, step
  underflow) where C++ succeeds. The surface has a double v knot (C0): the iso-u integrand
  jumps there. On a jump the fifth-order estimate is honest (about 0.2 h J), the embedded one
  about 5x lower; at 1e-8 with the moment components' scales the step across the jump would have
  to be below SM_EFF_ZERO. The C++ accepts a step whose real error exceeds its tolerance.
  Splitting at the surface's knots (WR-494, `precise_properties_split_at_surface_knots`) is the
  fix; it is the only path since 2026-11-05
  (plan/findings/mass.integrator--split-at-surface-knots-default-moves-cases-2026-11-05.md).
- `review_edge_constant_*`: a zero-length edge, values of order 1e-17 (noise) on both sides.
- `topology.face_p2` `ode_area_*_vertical`, `gen0000` and some `extra*` cases: the trim's
  `du/dt` is zero at the point asked; upstream returns the caller's array unchanged, Rust the
  integrand's value (zeros, one plus the components long). The other `extra*` cases move within
  the accuracy times the case's error scale.
- Not changed: `integrate_opt_in` (curve length, `SmIntegrator::IntegrateIt`) keeps the
  embedded estimate alone; the same blind spot applies there. Follow-on.
- NURBS-converted cylinder and cone areas converge to -1.95e-7 and -3.69e-7 of the closed form
  at every accuracy, before and after: the converted geometry, not the integration.
