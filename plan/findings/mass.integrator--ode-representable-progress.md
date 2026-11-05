# mass.integrator: nonzero intervals require representable progress (2026-11-05)

**Unit and PR:** `mass.integrator`, PR #2764. For scalar integration of `1e12`
over `[0, 1e-11]` at accuracy `1e-8`, the previous Rust result was `0` without
evaluation; the closed-form result is `10`. The reverse integral is `-10`.
The ODE driver treated its initial step `1e-12` as approximately equal to
`EFF_ZERO` and returned success before integrating. The new driver integrates
distinct finite endpoints and refuses a step that cannot advance to a finite
floating-point value, including a reduced retry. Equal finite endpoints still
preserve the initial array without evaluating unused controls.

This changes results for nonzero intervals previously skipped by approximate
endpoint or step comparisons, and changes stalled steps to named step-control
errors. The existing scalar `< EFF_ZERO` no-evaluation policy and WR-94 default
failure masking remain; named refusal is guaranteed only by direct ODE calls or
the scalar propagation-enabled path. Initial values remain otherwise ignored.
No general accuracy guarantee at large offsets is asserted.

Cash–Karp's existing method is retained (Cash and Karp, 1990,
doi:10.1145/79505.79507); the regression oracle is the constant's closed form
and representable advance. The existing `near_equal_large_domain`, `zero_step`,
and `step_epsilon_*` diagnostic fixtures remain unchanged. Their differences
are unverified: no new C++ comparison or confirmation was performed, and no
known-mismatch marker was added.
