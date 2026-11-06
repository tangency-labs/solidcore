# Projected-property local audit

The public LocalPropertyAnalysis contract searches for a zero from one guess and returns a found flag. It does not certify a global maximum or promise convergence to every existing root.

| Case | Native result (t, distance) | Rust result (t, distance) | Independent classification |
|---|---|---|---|
| gen0042 | not found | (2.568630556952937, 53.05559091626183) | New stationary minimum; normalized stationarity residual 1.7e-16. |
| gen0773 | not found | (0.02253138757000629, 33.08376518405103) | New stationary minimum; a separate maximum also exists. |
| gen1144 | (0.2424349354598795, 47.30990639718397) | (0.7897693109735117, 37.0895839403429) | Both valid roots; native maximum, Rust minimum. Changed basin, not a false native result. |
| gen1760 | (6.160833635900556, 1702.746539725838) | not found | Native valid maximum; Rust misses it from this guess. Lost local convergence. |
| gen2335 | not found | (0.2057706791318854, 1.526401246903282) | New stationary minimum. |
| gen2938 | (-1.416959240152698, 842.2781967874612) | (-1.635855706635469, 223.8362336456594) | Both valid roots; native maximum, Rust minimum. Changed basin, not a false native result. |
| gen2940 | (0.08124076285751645, 38.51742584776892) | not found | Native valid maximum; Rust misses it from this guess. Lost local convergence. |

All five Rust returned distances independently agree with the rational/polynomial curve evaluator within 3e-16 relative error. gen1144 has normalized stationarity residual 2.38e-11; the other four are below 1.5e-14. The two native maxima lost by the local solve are independently valid. The corresponding global queries were checked on the revised library with Single requests and tolerance 1e-7. Both retrieve the correct global maxima: gen1760 at t=6.160833635900567, distance 1702.7465397258381; gen2940 at t=0.08124076285751648, distance 38.517425847689154. These agree with independent per-span polynomial roots and endpoints. The pinned native global queries both return success with zero solutions on these identical inputs. This does not remove the two local convergence differences.

The global fixtures are `lpa_gen1760_global_max` and `lpa_gen2940_global_max` in
`fixtures/curve.curve_p1/issue_2840_local_root_global_maxima.json`. The first request
lies inside the natural domain. The second extends below zero; an independent
degree-five stationary polynomial over the entire requested Hermite interval
`[-0.19271549177165576,0.36770704114935876]` has exactly one real root in the
interval, at `0.08124076285751636`. There are no roots on the negative portion;
endpoint distances `14.503802936837797` and `20.978930958500555` are smaller.
Thus these two oracle comparisons cover their full requested intervals.

For p4, all 21 changes are only exact negations of values.out[1], the derivative. Twenty nonsingular cases independently match the analytic derivative. gen1477 is a zero-tangent endpoint: the evaluator substitutes a tiny tangent from the second derivative. Its sign change is confirmed, but this audit does not certify that regularized derivative as the derivative of the ordinary smooth curve equation.

Candidate `09495e082`; base `3f82278a14c1a6104db8e090a59b4bcf677ace01`. The base had zero failures in either unit (five parked p4 cases, three parked p5 cases). All changes below remain native mismatches; the local-root changes are not parked as C++ defects.

The smooth-curve oracle differentiates `F=-unit(P C') dot P(C-point)` using independently constructed homogeneous polynomial coordinates and quotient-rule derivatives, where `P=(N dot N) I-N N^T`. Stationary-root classification uses the sign of the derivative of the squared-distance numerator. This numerical calculation is not interval certification.

| p4 case | Native derivative | Rust derivative | Qualification |
|---|---:|---:|---|
| circle_property_21 | 0.64309503671411716 | -0.64309503671411716 | Matches independent smooth derivative |
| circle_property_21_null | 0.64309503671411716 | -0.64309503671411716 | Matches independent smooth derivative |
| gen0343 | 48.320386438490615 | -48.320386438490615 | Matches independent smooth derivative |
| gen0446 | 2.0158575156101239 | -2.0158575156101239 | Matches independent smooth derivative |
| gen0563 | 11113.792108902107 | -11113.792108902107 | Matches independent smooth derivative |
| gen0617 | 176.17894853692081 | -176.17894853692081 | Matches independent smooth derivative |
| gen0643 | 1328.1240378492116 | -1328.1240378492116 | Matches independent smooth derivative |
| gen0650 | 3.4752650729169097 | -3.4752650729169097 | Matches independent smooth derivative |
| gen1002 | 0.63172316674649853 | -0.63172316674649853 | Matches independent smooth derivative |
| gen1281 | -4.1947955102473973 | 4.1947955102473973 | Matches independent smooth derivative |
| gen1345 | 6.1274307413620921 | -6.1274307413620921 | Matches independent smooth derivative |
| gen1477 | -0.0042756974771357694 | 0.0042756974771357694 | Regularized zero tangent; sign only |
| gen1506 | 20.304791431505773 | -20.304791431505773 | Matches independent smooth derivative |
| gen1529 | 2313.4048875158164 | -2313.4048875158164 | Matches independent smooth derivative |
| gen1555 | 74.387085994451468 | -74.387085994451468 | Matches independent smooth derivative |
| gen1785 | 17493.959292360487 | -17493.959292360487 | Matches independent smooth derivative |
| gen1866 | -11.485183439048484 | 11.485183439048484 | Matches independent smooth derivative |
| gen2407 | 1.999531352594077 | -1.999531352594077 | Matches independent smooth derivative |
| gen2432 | -2090.9202731335936 | 2090.9202731335936 | Matches independent smooth derivative |
| gen2812 | 79.95773918004042 | -79.95773918004042 | Matches independent smooth derivative |
| gen2828 | -267.3591805617487 | 267.3591805617487 | Matches independent smooth derivative |
