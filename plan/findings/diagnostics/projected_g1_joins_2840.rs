//! Approximate cache continuity must not erase actual projected corner extrema.
//!
//! The two segments meet at an angle 2 atan(1/1024), below the cache's 0.5-degree
//! G1 threshold. Their shared point is nevertheless a genuine minimum or maximum
//! of the projected distance, depending on the query point.

use solidcore::curve::curve::{p1::GlobalPointSolveOptions, Curve};
use solidcore::geom::extent1d::Extent1d;
use solidcore::geom::geom_utility::ContinuityType;
use solidcore::geom::global_solver::{SolutionRequested, SolverOperation};
use solidcore::geom::solution_array::{ObjectKind, SolutionObject};
use solidcore::math::Vec3;
use solidcore::nurbs::b_spline_curve::BSplineCurve;
use solidcore::nurbs::crv_basic::Curve as NurbCurve;
use solidcore::nurbs::math::{CPoint, NOW};
use solidcore::Kernel;

#[derive(Clone, Copy, Debug, PartialEq)]
struct Handle;
impl SolutionObject for Handle {
    fn kind(&self) -> ObjectKind {
        ObjectKind::Curve
    }
}

fn polyline(height: f64) -> Curve {
    Curve::from(BSplineCurve::from_nurb(
        NurbCurve::new(
            vec![
                CPoint::new(-1.0, height, 0.0, NOW),
                CPoint::new(0.0, 0.0, 0.0, NOW),
                CPoint::new(1.0, height, 0.0, NOW),
            ],
            1,
            vec![0.0, 0.0, 1.0, 2.0, 2.0],
        )
        .unwrap(),
        3,
        false,
    ))
}

fn check_corner(operation: SolverOperation, point: Vec3, expected: f64) {
    let k = Kernel::default();
    let curve = polyline(1.0 / 1024.0);
    let (_, continuity) = curve.calculate_continuities(&k, 0.5).unwrap();
    assert!((continuity[1] as i32) >= (ContinuityType::G1 as i32));
    for requested in [SolutionRequested::Single, SolutionRequested::All] {
        let answers = curve
            .global_point_solve::<Handle>(
                &k,
                &curve.natural_interval(),
                operation,
                point,
                1e-8,
                requested,
                &GlobalPointSolveOptions {
                    opt_vectors: Some(&[Vec3::Z]),
                    ..Default::default()
                },
            )
            .unwrap();
        println!(
            "{operation:?} {requested:?}: {:?}",
            answers
                .as_slice()
                .iter()
                .map(|s| (s.start.parameters[0], s.start.solution_value))
                .collect::<Vec<_>>()
        );
        assert!(
            !answers.is_empty(),
            "the G1-tagged kink is the unique extremum"
        );
        for answer in answers.as_slice() {
            assert!(
                (answer.start.solution_value - expected).abs() <= 1e-8,
                "reported {}, expected {expected}",
                answer.start.solution_value
            );
            assert_eq!(answer.start.parameters[0], 1.0);
        }
    }
}

#[test]
fn near_collinear_join_is_projected_minimum() {
    check_corner(
        SolverOperation::ProjectedMinimize,
        Vec3::new(0.0, -1.0, 0.0),
        1.0,
    );
}

#[test]
fn near_collinear_join_is_projected_maximum() {
    check_corner(
        SolverOperation::ProjectedMaximize,
        Vec3::new(0.0, 1024.0, 0.0),
        1024.0,
    );
}

#[test]
fn exactly_collinear_join_does_not_become_an_extra_all_maximum() {
    let curve = polyline(0.0);
    let answers = curve
        .global_point_solve::<Handle>(
            &Kernel::default(),
            &Extent1d::from_min_max(0.0, 2.0),
            SolverOperation::ProjectedMaximize,
            Vec3::Y,
            10.0,
            SolutionRequested::All,
            &GlobalPointSolveOptions {
                opt_vectors: Some(&[Vec3::Z]),
                ..Default::default()
            },
        )
        .unwrap();
    assert_eq!(answers.len(), 2);
    assert!(answers
        .as_slice()
        .iter()
        .all(|answer| answer.start.parameters[0] == 0.0 || answer.start.parameters[0] == 2.0));
}
