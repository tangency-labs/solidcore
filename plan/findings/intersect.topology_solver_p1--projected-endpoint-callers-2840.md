# Local topology projected-query callers (#2840)

Candidate `09495e082` changes 23 cases in `intersect.topology_solver_p1` against
the pinned, unmodified v1.0.1 C++ reference. The exact-base report has no failures
in this unit. This record concerns that candidate, before the later affine-span
construction; it is not a certificate for subsequent revisions.

## Independent checks

The changed calls are projected point/curve, point/edge and vertex/edge minima,
including swapped and repeated calls accumulating into one solution array.
For each stored curve span, homogeneous polynomial coordinates produce the
squared-distance numerator; roots of its rational derivative and interval
endpoints provide an independent numerical optimum. This is a numerical
polynomial audit, not interval root isolation or proof of complete All-mode
enumeration.

All 24 returned Rust records agree with independent reevaluation within strict
`1e-9` relative error (largest observed error `1.29544e-12`, no absolute floor).
Every retained answer satisfies the original query tolerance and initial
incumbent threshold without an added allowance. Repeated calls add no exact
duplicate object/parameter records. Final squared incumbents equal the minimum
of the input incumbent and the returned squared distances.

Nineteen cases change an empty native result to a valid endpoint candidate:
`gen0268`, `gen0323`, `gen0662`, `gen0814`, `gen0865`, `gen1046`, `gen1104`,
`gen1240`, `gen1295`, `gen1571`, `gen1729`, `gen1958`, `gen2112`, `gen2163`,
`gen2261`, `gen2287`, `gen2330`, `gen2508`, `gen2966`.

Three native results are farther from the independent optimum than their
original tolerance permits:

| Case | Native distance | Rust distance | Query tolerance |
|---|---:|---:|---:|
| `gen2326` | `0.00040060002529845016` | `0.00018873029626767892` | `0.0001` |
| `gen2603` | `0.30420753099174874` | `0.0027330816133098457` | `1e-6` |
| `gen2818` | `68417.76422059476` | `23547.453225929043` | `100` |

`gen0729` preserves the minimum `0.02387226938059591` and adds the endpoint
`0.045217023584732174`, within its tolerance `0.1`. It also satisfies the input
incumbent zero plus that tolerance. This is an additional tolerated candidate,
not a corrected native optimum.

## Domain qualification

`gen0268`, `gen1240`, `gen1295` and `gen2163` include competing conic calls whose
requested intervals lie entirely above 1, while the constructed rational conic
has stored domain `[0,1]`. The existing cache traversal skips these disjoint
spans. The audit therefore certifies the returned answers only under the
existing stored-domain traversal. It does not settle whether those fixtures'
conic parameters should be interpreted as angular parameters or extrapolated,
and does not prove complete enumeration over a different domain contract.

The public domain question remains open on #2840. No tolerance was loosened and
these cases have not been parked as proven native defects.
