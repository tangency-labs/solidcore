#### WR-762. `SmODEIntegrator` floors every component's error scale at the absolute `SM_EFF_ZERO`, so an integral below 1e-12 is not held to the accuracy asked
- **Severity:** Wrong result
- **Area:** Numerical integration
- **Location:** `source/SMLib/src/SmIntegrator.cpp`, `SmODEIntegrator::RungeKutta5thOrderStep` (the in-place floor of each component's error scale at `SM_EFF_ZERO`) and `RungeKuttaODEIntegrate` (the scale is `|y| + |h y'|` plus the integrator's `tiny`, `SM_EFF_ZERO_SQ` by default). Ported as `OdeIntegrator::step_in` and `OdeIntegrator::integrate` in `crates/solidcore/src/mass/integrator.rs`.
- **What happens:** For an integral and an integrand far below `1e-12` the floor decides the scale, so each step is accepted when its error is below `1e-12 * eps` absolute, whatever the integral's size. The requested relative accuracy does not hold and nothing reports it.
- **Trigger:** `SmIntegrator::IntegrateIt` of `c0 sin(c1 t + c2)` with `c0 = -4.719540182810223e-16`, `c1 = 1.8640703447720268`, `c2 = 0.938454838283926` over `[0.4192727339336013, 2.1527686633843524]`, accuracy 1e-9: the closed form is 9.756956687989441e-17, upstream returns 9.756913811407506e-17 (relative error -4.39e-6, 4400 times the accuracy).
- **How found:** generated difftest case
- **Confirmed:** yes, case `gen0994` of `mass.integrator`
- **Rust port:** solidcore's `OdeIntegrator::integrate` no longer floors at an absolute constant (`plan/findings/mass.integrator--relative-scale-floor-2026-11-07.md`); `OdeIntegrator::step` keeps the in-place floor.
