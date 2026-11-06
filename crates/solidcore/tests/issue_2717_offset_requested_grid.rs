//! #2717: offset derivative requests preserve their named coordinates.
//!
//! The parabolic normal is differentiated independently; all four formula corrections
//! owned by #1914 vanish for these controls. Covers 624 regular requests.
#![allow(clippy::unwrap_used, clippy::indexing_slicing)] // Test oracles and fixture indexing.

use solidcore::geom::vector2d::Point2;
use solidcore::math::Vec3;
use solidcore::nurbs::b_spline_surface::BSplineSurface;
use solidcore::nurbs::b_spline_surface_n_lib::p2::EvaluateFlags;
use solidcore::nurbs::knots::Side;
use solidcore::surface::offset_surface::{OffsetSurface, MAX_DERIVS};
use solidcore::surface::surface::Surface;
use solidcore::Kernel;

/// Evaluation options matching the public regular and simple entry points.
fn flags(non_zero_tangents: bool, only_upper_half: bool) -> EvaluateFlags {
    EvaluateFlags {
        u_side: Side::Left,
        v_side: Side::Left,
        only_upper_half,
        non_zero_tangents,
    }
}

/// Exact quadratic x linear (or linear x quadratic) Bezier representation of
/// B(u,v)=(u,v,u²), or B(u,v)=(u,v,v²).
fn issue_2717_parabolic_base(transpose: bool) -> Surface {
    use solidcore::nurbs::math::{CPoint, NOW};
    let quadratic = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
    let linear = [0.0, 0.0, 1.0, 1.0];
    let curved = [(0.0, 0.0), (0.5, 0.0), (1.0, 1.0)];
    let b = if transpose {
        let net: Vec<Vec<CPoint>> = [0.0, 1.0]
            .iter()
            .map(|&u| {
                curved
                    .iter()
                    .map(|&(v, z)| CPoint::new(u, v, z, NOW))
                    .collect()
            })
            .collect();
        BSplineSurface::from_data(1, 2, 1, 2, 3, 5, &net, &linear, &quadratic).unwrap()
    } else {
        let net: Vec<Vec<CPoint>> = curved
            .iter()
            .map(|&(u, z)| {
                [0.0, 1.0]
                    .iter()
                    .map(|&v| CPoint::new(u, v, z, NOW))
                    .collect()
            })
            .collect();
        BSplineSurface::from_data(2, 1, 2, 1, 5, 3, &net, &quadratic, &linear).unwrap()
    };
    Surface::from(b)
}

/// Differentiate (-2t,0,1)/sqrt(1+4t²) directly, independently of the cross-product
/// recurrence used by OffsetSurface. Transposition puts its curved component in y.
#[allow(clippy::disallowed_methods)] // Independent closed-form oracle.
fn issue_2717_parabolic_jet(transpose: bool, distance: f64, u: f64, v: f64) -> [[Vec3; 4]; 4] {
    let t = if transpose { v } else { u };
    let s = 1.0 + 4.0 * t * t;
    let axis = if transpose {
        Vec3::new(0.0, 1.0, 0.0)
    } else {
        Vec3::new(1.0, 0.0, 0.0)
    };
    let z = Vec3::new(0.0, 0.0, 1.0);
    let normal = [
        axis * (-2.0 * t / s.sqrt()) + z / s.sqrt(),
        axis * (-2.0 / s.powf(1.5)) + z * (-4.0 * t / s.powf(1.5)),
        axis * (24.0 * t / s.powf(2.5)) + z * ((-4.0 + 32.0 * t * t) / s.powf(2.5)),
        axis * (24.0 * (1.0 - 16.0 * t * t) / s.powf(3.5))
            + z * (48.0 * t * (3.0 - 8.0 * t * t) / s.powf(3.5)),
    ];
    let mut out = [[Vec3::ZERO; 4]; 4];
    out[0][0] = Vec3::new(u, v, t * t) + distance * normal[0];
    for (order, n) in normal.iter().enumerate().skip(1) {
        let base = match order {
            1 => axis + z * (2.0 * t),
            2 => z * 2.0,
            _ => Vec3::ZERO,
        };
        if transpose {
            out[0][order] = base + distance * *n;
        } else {
            out[order][0] = base + distance * *n;
        }
    }
    if transpose {
        out[1][0] = Vec3::new(1.0, 0.0, 0.0);
    } else {
        out[0][1] = Vec3::new(0.0, 1.0, 0.0);
    }
    out
}

/// Compare finite values under the strict relative rule, including exact zeros.
fn issue_2717_assert_value(got: Vec3, want: Vec3, context: &str) {
    for (a, b) in [(got.x, want.x), (got.y, want.y), (got.z, want.z)] {
        assert!(
            a.is_finite() && b.is_finite() && (a - b).abs() <= 1e-9 * a.abs().max(b.abs()),
            "{context}: {got:?} != {want:?}"
        );
    }
}

/// Named pure derivatives agree with their closed forms even on narrow grids.
#[test]
fn issue_2717_named_upper_derivatives_match_the_parabolic_offset() {
    let k = Kernel::default();
    let off = OffsetSurface::new(&k, 0.3, issue_2717_parabolic_base(false), true, 0).unwrap();
    let expected = issue_2717_parabolic_jet(false, 0.3, 0.5, 0.5);
    let control = off
        .evaluate_simple(&k, Point2::new(0.5, 0.5), 3, 3, flags(false, true))
        .unwrap();
    issue_2717_assert_value(control[2][0], expected[2][0], "equal-order Duu control");
    issue_2717_assert_value(control[3][0], expected[3][0], "equal-order Duuu control");
    for (hu, hv) in [(2, 0), (3, 0), (3, 1)] {
        let out = off
            .evaluate_simple(&k, Point2::new(0.5, 0.5), hu, hv, flags(false, true))
            .unwrap();
        issue_2717_assert_value(out[2][0], expected[2][0], &format!("({hu},{hv}) Duu"));
        if hu == 3 {
            issue_2717_assert_value(out[3][0], expected[3][0], &format!("({hu},{hv}) Duuu"));
        }
    }
}

/// The larger-request control is checked independently in both axes and signs.
#[test]
fn issue_2717_equal_order_controls_match_independent_parabolic_jets() {
    let k = Kernel::default();
    for transpose in [false, true] {
        for distance in [-0.3, 0.3] {
            let off =
                OffsetSurface::new(&k, distance, issue_2717_parabolic_base(transpose), true, 0)
                    .unwrap();
            // Avoid roots of a nonzero derivative polynomial while keeping the base exact.
            for t in [0.125, 0.5, 0.875] {
                let (u, v) = if transpose { (0.5, t) } else { (t, 0.5) };
                let expected = issue_2717_parabolic_jet(transpose, distance, u, v);
                let out = off
                    .evaluate_simple(&k, Point2::new(u, v), 3, 3, flags(false, true))
                    .unwrap();
                for i in 0..=3 {
                    for j in 0..=3 - i {
                        issue_2717_assert_value(
                            out[i][j],
                            expected[i][j],
                            "independent analytic control",
                        );
                    }
                }
            }
        }
    }
}

/// Every supported regular request keeps the same named derivatives.
#[test]
fn issue_2717_every_supported_request_keeps_its_derivative_names() {
    let k = Kernel::default();
    for transpose in [false, true] {
        for distance in [-0.3, 0.3] {
            let off =
                OffsetSurface::new(&k, distance, issue_2717_parabolic_base(transpose), true, 0)
                    .unwrap();
            // Dyadic inputs preserve exact zero mixed derivatives; no absolute floor.
            for t in [0.25, 0.5, 0.75] {
                let (u, v) = if transpose { (0.5, t) } else { (t, 0.5) };
                let uv = Point2::new(u, v);
                let control = off
                    .evaluate_simple(&k, uv, 3, 3, flags(false, true))
                    .unwrap();
                for upper in [false, true] {
                    for hu in 0..=3 {
                        for hv in 0..=3 {
                            if !upper && hu + hv > MAX_DERIVS {
                                continue;
                            }
                            for simple in [false, true] {
                                let out = if simple { off.evaluate_simple(&k, uv, hu, hv, flags(false, upper)) }
                        else { off.evaluate(&k, uv, hu, hv, flags(true, upper)) }
                        .unwrap_or_else(|e| panic!("valid request ({hu},{hv},upper={upper},simple={simple}): {e}"));
                                for i in 0..=hu {
                                    for j in 0..=hv {
                                        if upper && i + j > hu.max(hv) {
                                            continue;
                                        }
                                        let context = format!("transpose={transpose}, d={distance}, t={t}, request=({hu},{hv},{upper},{simple}), D({i},{j})");
                                        issue_2717_assert_value(out[i][j], control[i][j], &context);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

/// Unsupported orders retain their existing bounded refusal.
#[test]
fn issue_2717_unsupported_orders_still_fail_before_allocation() {
    let k = Kernel::default();
    let off = OffsetSurface::new(&k, 0.3, issue_2717_parabolic_base(false), true, 0).unwrap();
    for upper in [false, true] {
        for (hu, hv) in [
            (4, 0),
            (0, 4),
            (usize::MAX, 0),
            (0, usize::MAX),
            (usize::MAX, usize::MAX),
        ] {
            assert!(off
                .evaluate_simple(&k, Point2::new(0.5, 0.5), hu, hv, flags(false, upper))
                .is_err());
            assert!(off
                .evaluate(&k, Point2::new(0.5, 0.5), hu, hv, flags(true, upper))
                .is_err());
        }
    }
}
