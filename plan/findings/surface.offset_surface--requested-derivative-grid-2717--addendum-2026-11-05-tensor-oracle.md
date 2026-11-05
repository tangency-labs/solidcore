# surface.offset_surface: independent tensor oracle and tiny derivatives (2026-11-05)

Updates: `surface.offset_surface--requested-derivative-grid-2717.md`.

Unit: `surface.offset_surface`. PR: #2765, still a draft on #1914. Library code and difftest fixtures, generators, tolerances, floors, mismatch
markers and unit statuses are unchanged. New captured inputs serve only the
independent mathematical reproducer.

## Independent calculation

`harness/oracles/offset_grid_tensor.py` evaluates the tensor B-spline definition,
then `S + d*(Su cross Sv)/|Su cross Sv|`, using standard-library Decimal Taylor
coefficients in both variables. It consumes only captured inputs. Degrees 1–3,
orders through three, nonrational nets and interior parameters away from knots
are supported. Only each captured script's first request is evaluated. It does
not reproduce later mutations, extension, pole substitution or the full script.

The basis definition is the Cox–de Boor recurrence documented by
[SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.BSpline.html);
SciPy is not a dependency and no implementation is imported. At 140 digits, the
point, mixed `uv` derivative and third `u` derivative of an offset of
`S=(u,v,u²+v²)` agree with separately derived closed forms to absolute `1e-125`.
A second calculation for bilinear patches forms `Su`, `Sv` and `Suv` as exact
fractions of the binary-float inputs. With `n=Su cross Sv`, `A=n dot n_u`,
`B=n_u dot n_u`, `L=|n|` and `n_uu=0`, it uses
`Duu=d*(-(2*n_u*A+n*B)/L³+3*n*A²/L⁵)` directly. It uses neither spline-basis
recursion nor Taylor arithmetic. For the two bilinear cases below, the Taylor
and closed-form values agree within `1.1e-123` relative at 140 digits.
Precision convergence is a diagnostic, not interval arithmetic or a universal
error certificate.

```sh
python3 harness/oracles/offset_grid_tensor.py --case gen0523
python3 harness/oracles/offset_grid_tensor.py --case gen1233
python3 harness/oracles/offset_grid_tensor.py --case gen2816
```

The three complete inputs and source seeds are separate files in
`harness/oracles/offset_grid_tensor_cases/`.

## Well-conditioned bicubic case

`gen0523`, seed `7714667220713742484`, requests upper `(2,0)` through
`EvaluateSimple` at `(0.5326216929082825, 0.16048492851110552)` on a nonrational
bicubic base. Its normal length is approximately `4.269293140375243`.

| Duu component | Input-only oracle, rounded | Fresh Rust | Fresh C++ request |
| --- | ---: | ---: | ---: |
| x | 7.71081091620462 | 7.710810916204617 | 0.8251605633693019 |
| y | 4.83384618851690 | 4.833846188516903 | -0.5884877400243681 |
| z | -0.211527114990600 | -0.21152711499061783 | 5.989497909150047 |

Rust's maximum relative error against the 140-digit values is `8.33e-14`.
The 80-/140-digit results differ by at most `4.04e-78` absolute across the
requested point and derivatives. This independently confirms the C++ wrong
result in this case (WR-162); it does not establish every generated result.
Only these three components differ in the full native report for this case.
A separately rerun native `(3,3,true)` request on the same prefix returns exactly
the three Rust values in the table. The larger request is a consistency control,
not the mathematical oracle.

## Tiny bilinear derivatives remain unresolved for relative parity

The following are nonzero derivatives of the exact binary-float input nets,
not exact zeros. The independent bilinear closed form confirms their size.
The examples use component x of the first requested `Duu`:

| Case / source seed | Input-only oracle, rounded | Fresh Rust | Fresh C++ request |
| --- | ---: | ---: | ---: |
| gen1233 / 7714667220713743194 | -6.96737327e-28 | -2.327429053988993e-27 | 2.2737367544323206e-12 |
| gen2816 / 7714667220713744777 | -1.13679237e-31 | -8.917345340782167e-31 | -1.5685122927778304e-45 |

Neither engine establishes relative `1e-9` agreement with these independent
values. Their relative error is sensitive to cancellation; this observation
is not a proved roundoff bound, a proposed absolute floor, a kernel-fix claim
or proof that these cases can be parked. In particular, agreement with a
larger C++ request does not prove mathematical accuracy. Unlike `gen0956` in
the earlier record, exact rational `Suv` is nonzero in both inputs.

## Isolated validation and unchanged outcomes

Both source states were rebuilt with private Cargo target **and** build
directories: parent `b670d5f9c`, fixed head `4f409e89a`. The native libraries
are pinned v1.0.1, commit `3fd5b9548d19aa8c12c29bc9637a96d42a06f8a5`.
The fresh parent checkout initially lacked the upstream header path; that
setup failure was retained, the existing pinned headers were linked, and
no source or reference-library build was changed to correct the setup.

| Native diagnostic | Cases | Passed | Ordinary failed | Inherited parked | Deferred |
| --- | ---: | ---: | ---: | ---: | ---: |
| Parent plus regressions | 3182 | 2663 | 236 | 283 | 0 |
| Fixed requested grid | 3182 | 2584 | 315 | 283 | 0 |

All 79 transitions are still `pass -> fail`, with unchanged case/source
identities. The fixed run reproduces all 315 retained full inputs from the
historical run; compact passing records do not establish full output identity.
The separately relinked request-shape probes also reproduce 95 unequal
requests: 81 C++ requests change named coordinates when enlarged, while 1143
finite larger-request C++ coordinates agree with Rust. Three nonfinite
components in gen0881 are excluded, not counted as agreement. These controls
are not independent mathematical oracles.

Both native commands exit 1. The fixed head separately passes all 27 offset
tests, the full kernel library (3531 pass, four existing ignores), strict core
all-target Clippy and WASM. The parent has the two expected red regressions.
These are local diagnostics, not a successful exact-base baseline.

All ordinary failures remain. This record does not release either draft,
certify all 79 moved cases, or take the dependency author's formula scope.
