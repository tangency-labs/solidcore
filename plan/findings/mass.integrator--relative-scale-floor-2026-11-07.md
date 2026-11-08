# mass.integrator: an integral below 1e-12 is controlled relative to its own size, not to an absolute floor (2026-11-07)

- **Unit:** `mass.integrator` (#2886).
- **Input:** difftest `gen0994`, the scalar `integrate` of `c0 sin(c1 t + c2)` with
  `c0 = -4.719540182810223e-16`, `c1 = 1.8640703447720268`, `c2 = 0.938454838283926` over
  `[0.4192727339336013, 2.1527686633843524]`, accuracy 1e-9. Closed form
  `-c0/c1 (cos(c1 b + c2) - cos(c1 a + c2))` = 9.756956687989441e-17.
- **Values:** upstream 9.756913811407506e-17 (relative error -4.39e-6); solidcore before
  9.75691082970234e-17 (-4.70e-6); after: within the bound the control implies (test
  `an_integral_below_the_old_floor_is_held_relative_to_its_size`). Each step's error is
  held to `accuracy * s_i`, not their sum to `accuracy` times the final value (0.2 of the
  amplitude): the total is at most `accuracy * ((N + 2) span G + N tiny)` for `N` accepted
  steps and `G` the largest `|g|`, since `s_i <= peak + 2 |h_i| G + tiny` and
  `peak <= span G`. The same integrand scaled by a power of two gives the same steps and
  value bit for bit.
- **Cause:** `OdeIntegrator::step_in` replaced every component's error scale below
  `EFF_ZERO` (1e-12) by `EFF_ZERO`, and `integrate` added `tiny = EFF_ZERO_SQ = 1e-24` to the scale
  (`|y| + |h y'| + tiny`). For an integral and an integrand far below 1e-12 the floor decided:
  a step was accepted when its error was below `1e-12 * accuracy` absolute. Upstream has the
  same floor (`SM_EFF_ZERO`), so both kernels agreed on the miss.
- **Change:** `integrate` no longer floors at an absolute constant. Each step reads the scale
  as the largest of the caller's `|y| + |h y'|` (and lower bound), `|h| max_j |g_j|` over
  the step's six samples, the running maximum of `|y|` over the accepted nodes, and
  `EFF_ZERO * |x2 - x1| * G`, `G` the largest sample of the component evaluated so far
  (retries included, never falling). The last term is the one that does not shrink with
  the step: at an integrand zero at the start (`k t^4`) the estimate and the sample term
  both fall like `h^5`, their ratio is constant, and the steps would shrink to underflow;
  against `G`, fixed over the retries, the estimate falls below `accuracy * s`. `EFF_ZERO`
  (upstream's number) is now the relative factor of that floor, `1e-12` of `span * G`. The sample
  term does not vanish at a zero of the integrand and the peak does not where `y` passes
  through zero (a centred moment). All four scale exactly with the integrand: a power of
  two gives the same steps and `k` times the value bit for bit (test
  `the_control_is_invariant_under_the_integrand_scale`), above `tiny`. Near the subnormal
  range (amplitudes about 1e-300) `EFF_ZERO * span * G` can underflow and the invariance is
  not claimed. `G` takes the samples of every trial step, so the largest trial step sets
  the floor and the step count depends on the initial trial step as well as on the
  integrand: the invariance is under `k`, not under the first step (test
  `the_steps_depend_on_the_first_trial_step_and_not_on_the_scale`). A scale of exactly zero means
  every sample is zero and the estimates are exactly zero. `tiny` (default `EFF_ZERO_SQ`)
  is unchanged: an additive constant that still holds an integral below about 1e-24 to an
  absolute bound, a known limit (documented on `OdeIntegrator::new`); changing the default
  would move every caller of the plain constructor and is not what #2886 reports. The
  scale-invariance tests build `OdeIntegrator::new(10_000, 0.0)`.
  `OdeIntegrator::step` (public) keeps upstream's in-place floor, which its goldens record.
  The final owed-error check reads the same size as the control (the peak, not the final
  value), so an integral that cancels is not refused for ending at zero (test
  `the_owed_error_of_a_cancelling_integral_is_checked_against_its_peak`, which scans the
  samples' rounding bound for a run whose steps owe part of the accuracy times the peak
  and one whose steps owe more). To read the owed errors the steps of `integrate` are in a
  private `run`, which returns them; `integrate` makes the check.
- **With `rate`:** an integration given `OdeOptions::rate` (the curve length) keeps the step
  control it has, its in-place floor and its share of the range; the relative scale reads only
  where `rate` is `None`. The owed-error check reads the integral and the range's share with
  `rate`, the peak without.
- **Limit:** an integrand analytically zero but evaluated as noise above `EFF_ZERO` of its own
  size cannot converge: a typed `SM_ERR` at `OP_STEP_CONTROL`, nothing loosened.
- **Unchanged:** the rounding floor and the debt check (already relative), and the
  unachievable-accuracy guard. A request below the integrand's own rounding without a
  rounding bound ends in step underflow or the step budget: `SM_ERR` at `OP_STEP_CONTROL`.
- **Cases that move:** every integration with a component whose scale was below 1e-12 (the
  scalar `gen*` cases of small `c0`, the zero-length `review_edge_constant_*`, edges whose
  length approximation over ten is below 1e-12), and any whose scale was dominated by the
  `1e-24` of `tiny` or whose steps start at a zero of the integrand (the sample term raises
  the scale there). Each moved case is graded against its closed form or an independent
  reference in the pull request.
- **Steps, `k t^4` on `[0, 1]`, `tiny = 0`** (measured on this change and on the unchanged integrator):

  | k, accuracy | before | after |
  |---|---|---|
  | 1, 1e-6 | ok, 39 steps | ok, 51 steps |
  | 1, 1e-9 | `step_control` error | ok, 186 steps |
  | 1, 1e-12 | `step_control` error | ok, 724 steps |
  | 2^-80, 1e-6 and 1e-12 | ok, 3 steps (the absolute floor swallowed the integral) | ok, 51 and 724 steps, the same as k = 1 |

  So a tight accuracy on an integrand that is zero at the start now converges where it did
  not, and costs the steps the accuracy needs at any size.

  The first trial step also moves the count (`G` is set by it): the table above is for a
  first trial step of a tenth of the interval; see the test for `0.1` and `0.5`.
