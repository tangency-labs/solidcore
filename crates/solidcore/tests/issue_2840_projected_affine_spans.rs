//! Projected extrema of affine spans must not depend on their knot widths (#2840).
//!
//! The middle segment runs from (-1, 0) to (1, 1). Its supporting line is
//! x - 2y + 1 = 0, so (0, y) has distance |1 - 2y| / sqrt(5).
//! For y = 0.2 or 0.3 its perpendicular foot lies strictly inside that segment
//! and is closer than either adjoining horizontal segment. Raw knots are kept:
//! rebuilding the polyline with chord-length parameters would erase the defect.
#![allow(clippy::unwrap_used, clippy::indexing_slicing)]

use solidcore::curve::curve::{p1::GlobalPointSolveOptions, Curve};
use solidcore::geom::extent1d::Extent1d;
use solidcore::geom::global_solver::{SolutionRequested, SolverOperation};
use solidcore::geom::solution_array::{ObjectKind, SolutionObject, SolutionType};
use solidcore::math::{Vec3, EFF_ZERO_SQRT};
use solidcore::nurbs::b_spline_curve::p1::BSplineCurve;
use solidcore::nurbs::crv_basic::Curve as NurbCurve;
use solidcore::nurbs::math::CPoint;
use solidcore::Kernel;

#[derive(Clone, Copy, PartialEq)]
struct Object;

impl SolutionObject for Object {
    fn kind(&self) -> ObjectKind {
        ObjectKind::Other
    }
}

fn raw_polyline(gap: f64, parameter_scale: f64) -> (Curve, Vec<f64>) {
    let knots: Vec<_> = [0.0, 0.0, 1.0, 1.0 + gap, 2.0, 2.0]
        .into_iter()
        .map(|t| t * parameter_scale)
        .collect();
    assert!(knots[2] < knots[3]);
    let control = [(-2.0, 0.0), (-1.0, 0.0), (1.0, 1.0), (2.0, 1.0)]
        .into_iter()
        .map(|(x, y)| CPoint {
            x,
            y,
            z: 0.0,
            w: 1.0,
        })
        .collect();
    let nurb = NurbCurve::new(control, 1, knots.clone()).unwrap();
    assert_eq!(nurb.knots(), knots.as_slice());
    let spline = BSplineCurve::from_nurb(nurb, 3, true);
    assert_eq!(spline.gw_nurb().unwrap().knots(), knots.as_slice());
    (Curve::from(spline), knots)
}

fn check(
    gap: f64,
    parameter_scale: f64,
    point: Vec3,
    normal_scale: f64,
    operation: SolverOperation,
    geometric_distance: f64,
) {
    let (curve, knots) = raw_polyline(gap, parameter_scale);
    let kernel = Kernel::default();
    let interval = Extent1d::from_min_max(knots[0], knots[5]);
    let normal = Vec3::Z * normal_scale;
    // The API's double-cross metric scales by |normal| squared.
    let expected = geometric_distance * normal_scale * normal_scale;
    let solutions = curve
        .global_point_solve::<Object>(
            &kernel,
            &interval,
            operation,
            point,
            EFF_ZERO_SQRT,
            SolutionRequested::Single,
            &GlobalPointSolveOptions {
                opt_vectors: Some(&[normal]),
                ..Default::default()
            },
        )
        .unwrap();
    assert_eq!(solutions.len(), 1);
    let answer = &solutions.as_slice()[0];
    assert_eq!(answer.solution_type, SolutionType::SingleValue);
    let t = answer.start.parameters[0];
    assert!(t >= knots[0] && t <= knots[5]);
    let evaluated = curve.evaluate_point(&kernel, t).unwrap();
    let actual = normal.cross(point - evaluated).cross(normal).length();
    let reported = answer.start.solution_value;
    assert!(actual > 0.0 && reported > 0.0);
    assert!(
        (reported - actual).abs() <= 1e-9 * reported.abs().max(actual.abs()),
        "reported metric differs from reevaluation: gap={gap}, scale={parameter_scale}, t={t}, reported={reported}, actual={actual}"
    );
    assert!(
        (actual - expected).abs() <= EFF_ZERO_SQRT,
        "missed affine extremum: gap={gap}, scale={parameter_scale}, t={t}, actual={actual}, expected={expected}"
    );
}

#[test]
fn unrepresentable_affine_minimum_returns_error_with_default_kernel() {
    let (curve, knots) = raw_polyline(1.5e-12, 1.0);
    let tolerance = 1e-12;
    let expected = 0.4 / 5.0_f64.sqrt();
    // Squared distance is a convex quadratic with its unique minimum at s=.46.
    // Check the rounded parameter and its immediate neighbors: these include
    // both floating-point parameters bracketing the exact perpendicular foot.
    let rounded = knots[2] + (knots[3] - knots[2]) * 0.46;
    let best_representable = [rounded.next_down(), rounded, rounded.next_up()]
        .into_iter()
        .map(|t| {
            let s = (t - knots[2]) / (knots[3] - knots[2]);
            ((2.0 * s - 1.0).powi(2) + (s - 0.3).powi(2)).sqrt()
        })
        .fold(f64::INFINITY, f64::min);
    assert!(best_representable - expected > tolerance);
    let result = curve.global_point_solve::<Object>(
        &Kernel::default(),
        &Extent1d::from_min_max(knots[0], knots[5]),
        SolverOperation::ProjectedMinimize,
        Vec3::new(0.0, 0.3, 0.0),
        tolerance,
        SolutionRequested::Single,
        &GlobalPointSolveOptions {
            opt_vectors: Some(&[Vec3::Z]),
            ..Default::default()
        },
    );
    assert!(
        result.is_err(),
        "an unrepresentable answer must return SolidError, not empty or partial success"
    );
}

#[test]
fn narrow_raw_spans_retain_interior_projected_minima() {
    for gap in [1.5e-12, 1.9e-12, 1e-8] {
        for y in [0.3, 0.2] {
            check(
                gap,
                1.0,
                Vec3::new(0.0, y, 0.0),
                1.0,
                SolverOperation::ProjectedMinimize,
                (1.0 - 2.0 * y) / 5.0_f64.sqrt(),
            );
        }
    }
}

#[test]
fn affine_minimum_obeys_nonunit_projection_metric() {
    for normal_scale in [0.5, -2.0] {
        check(
            1.9e-12,
            1.0,
            Vec3::new(0.0, 0.3, 9.0),
            normal_scale,
            SolverOperation::ProjectedMinimize,
            0.4 / 5.0_f64.sqrt(),
        );
    }
}

#[test]
fn affine_minimum_is_invariant_under_parameter_scaling() {
    for (gap, scale) in [(1e-8, 1e-3), (1.5e-12, 1e6)] {
        check(
            gap,
            scale,
            Vec3::new(0.0, 0.2, 0.0),
            1.0,
            SolverOperation::ProjectedMinimize,
            0.6 / 5.0_f64.sqrt(),
        );
    }
}

#[test]
fn affine_boundary_minimum_and_maximum_remain_candidates() {
    check(
        1.5e-12,
        1.0,
        Vec3::new(-3.0, -0.4, 0.0),
        1.0,
        SolverOperation::ProjectedMinimize,
        1.16_f64.sqrt(),
    );
    check(
        1.5e-12,
        1.0,
        Vec3::new(0.0, 0.3, 0.0),
        2.0,
        SolverOperation::ProjectedMaximize,
        4.49_f64.sqrt(),
    );
}
