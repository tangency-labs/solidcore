# surface.offset_surface: requested derivative coordinates (#2717, 2026-11-05)

Unit: `surface.offset_surface`. PR: https://github.com/tangency-labs/solidcore/pull/2765.
Stacked dependency: #1914 at `052366249e42cfe4ad9defb109fea81c9060a0a3`.

## Reproducer and independent result

Use the quadratic Bezier surface `B(u,v)=(u,v,u²)` with offset distance
`0.3` at `u=v=0.5`. For an upper-half derivative request `(hu,hv)=(3,0)`
or `(3,1)`, the parent returns approximately `-1.08e-62` for `Duu.z`.
The closed-form offset normal `(-2u,0,1)/sqrt(1+4u²)` gives
`Duu=(0.6363961030678927,0,2.2121320343559643)` and
`Duuu=(-1.909188309203678,0,0.6363961030678926)`.
The equal-order `(3,3)` and lower-half `(3,0)` controls agree with that
independent result. All four #1914 correction terms vanish on this base.

The parent computes mixed outputs using only the maximum requested order
in upper-half mode. A flat slot `i*(hv+1)+j` can therefore alias another
named derivative when `j>hv`, or exceed the allocated grid when `i>hu`.
It can also read normal-derivative locals whose prerequisite output was
never requested. WR-162 records the original second-order manifestation;
#2717 includes third-order aliases.

## Change and native evidence

Require each mixed derivative to lie inside both the requested rectangle
and the upper-half total-order budget. The same conditions control regular
calculations and pole step-off output copies. Derivative arithmetic,
operation order, recursion bounds and order limits are unchanged.

Four new tests cover the independent closed form in both curved directions,
both offset signs, every supported derivative request through order three
(624 regular requests across `Evaluate` and `EvaluateSimple`), and rejection
of unsupported/extreme orders. Comparisons use relative `1e-9`, without an
absolute floor. Dyadic parameters keep the constructed zero mixed terms exact.
The existing plane `(0,2,true)` request is retained and now checks its exact
closed form instead of expecting the old out-of-bounds error.

On the unmodified parent plus the new tests: two pass and two fail for the
named-coordinate defect and valid-request rejection. After the mask change:
all 27 offset-surface module tests pass in the optimized difftest profile.
The nonregular pole/substitution behavior has no new independent oracle in
this change; the parent PR's diagnostic hold remains in force.

## C++ diagnostic baseline and limits

The unchanged full unit run uses default seed count 16 (the runner's minimum
generated pool is 3000), the existing matching prebuilt reference key
`v1-3fd5b9548d19-726643fb17a52b7d`, and the standard separate-process runner.
Before the grid fix: **3182 cases, 2663 passed, 236 failed, 283 inherited
parked, zero deferred**. This is an ordinary numeric FAIL and reproduces
#1914's published counts. No C++ values or state feed the Rust implementation.

After the fix: **3182 cases, 2584 passed, 315 failed, 283 inherited parked,
zero deferred**. All 79 classification changes are `pass -> fail`; every
changed script contains an unequal upper-half derivative request. This is
localization evidence, not an independent proof that each generated C++
result is wrong. For example `gen0094`, `values.out[62]`, changes from a
passing comparison to C++ `-0.18310631981091818` versus Rust
`-4.104584168595843`. Its script includes upper `(2,0)` and `(2,1)` requests.
The independent regular parabolic oracle proves the requested-grid defect;
it does not explain all generated geometry, pole or substitution outcomes.
All 315 diagnostic failures remain ordinary failures. The other 283 parked
entries are inherited unchanged from #1914.

The original and fixed reports have identical case names, source identities,
seed count, tolerance, absolute floors and case count. Compact passing records
do not establish bitwise equality of every C++ value. Clippy with warnings
denied, wasm32 checking and workspace formatting pass. The full optimized kernel library suite passes: 3531 passed, zero failed,
four existing ignored tests.

No unit promotion,
new mismatch marker, generator change, case exclusion, guard broadening or
tolerance change is part of this follow-up. Neither PR's diagnostic hold is
released by the native regression results.

## Native request-shape localization

A separate Rust-only diagnostic replays the original prefix before every
unequal upper evaluation in the 79 newly failing scripts. It runs that
request and an otherwise identical `(3,3,true)` control from the same prefix,
then compares only the named coordinates requested by the original call.
All 95 requests complete successfully; all 1146 compared components agree
within relative `1e-9`, without an absolute floor. Neither C++ output nor
C++ process state is used by this check. Script prefixes preserve prior
operations and extension state rather than evaluating on a fresh replacement
surface. This is request-shape consistency evidence, not an independent
mathematical oracle for those generated surfaces, and it does not justify
parking their C++ mismatches or releasing either draft hold.

## Newly failing case identities

`gen0094`, `gen0125`, `gen0138`, `gen0153`, `gen0242`, `gen0405`, `gen0450`, `gen0511`, `gen0513`, `gen0523`, `gen0548`, `gen0576`, `gen0597`, `gen0605`, `gen0632`, `gen0649`, `gen0685`, `gen0711`, `gen0764`, `gen0776`, `gen0787`, `gen0790`, `gen0853`, `gen0857`, `gen0861`, `gen0881`, `gen0887`, `gen0934`, `gen0956`, `gen0989`, `gen1062`, `gen1085`, `gen1172`, `gen1183`, `gen1205`, `gen1219`, `gen1233`, `gen1346`, `gen1415`, `gen1457`, `gen1461`, `gen1526`, `gen1533`, `gen1561`, `gen1604`, `gen1657`, `gen1687`, `gen1714`, `gen1755`, `gen1810`, `gen1867`, `gen1869`, `gen1954`, `gen2059`, `gen2078`, `gen2111`, `gen2115`, `gen2125`, `gen2220`, `gen2224`, `gen2275`, `gen2358`, `gen2379`, `gen2432`, `gen2452`, `gen2460`, `gen2468`, `gen2603`, `gen2639`, `gen2708`, `gen2730`, `gen2796`, `gen2816`, `gen2820`, `gen2835`, `gen2912`, `gen2940`, `gen2982`, `gen2998`.
