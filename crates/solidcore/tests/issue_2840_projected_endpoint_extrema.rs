//! Constrained projected extrema include finite curve endpoints (#2840).
//!
//! Oracle: the squared projected distance to our line is
//! `scale^2 * ((t - foot)^2 + transverse^2)`. Thus the minimum is
//! the clamped foot and the maximum is at the farther endpoint.
//! Comparisons use 1e-9 relative error, with a 1e-12 absolute floor for
//! parameters/distances near zero (not a solver or topology tolerance).
#![allow(clippy::unwrap_used, clippy::indexing_slicing)]

use solidcore::curve::curve::{p1::GlobalPointSolveOptions, Curve};
use solidcore::curve::line::Line;
use solidcore::geom::extent1d::Extent1d;
use solidcore::geom::global_solver::{SolutionRequested, SolverOperation};
use solidcore::geom::solution_array::{ObjectKind, SolutionObject, SolutionType};
use solidcore::math::Vec3;
use solidcore::nurbs::b_spline_curve::p1::BSplineCurve;
use solidcore::Kernel;

#[derive(Clone, Copy, PartialEq)]
struct Object;
impl SolutionObject for Object {
    fn kind(&self) -> ObjectKind {
        ObjectKind::Other
    }
}

fn close(actual: f64, expected: f64) -> bool {
    (actual - expected).abs() <= 1e-12_f64.max(1e-9 * actual.abs().max(expected.abs()))
}

#[allow(clippy::too_many_arguments)]
fn check_case(
    scale: f64,
    rotated: bool,
    reversed: bool,
    spline: bool,
    bounds: [f64; 2],
    foot: f64,
    op: SolverOperation,
    request: SolutionRequested,
) {
    let k = Kernel::default();
    let (axis, transverse, normal) = if rotated {
        (Vec3::X, Vec3::new(0.0, 0.8, -0.6), Vec3::new(0.0, 0.6, 0.8))
    } else {
        (Vec3::new(0.6, 0.8, 0.0), Vec3::new(-0.8, 0.6, 0.0), Vec3::Z)
    };
    let origin = Vec3::new(3.0, -2.0, 5.0) * scale;
    let (start, end) = if reversed {
        (origin + axis * scale, origin)
    } else {
        (origin, origin + axis * scale)
    };
    let curve = if spline {
        Curve::from(BSplineCurve::create_line_segment(3, start, end).unwrap())
    } else {
        Curve::from(Line::create_line_segment(3, start, end, None).unwrap())
    };
    // The normal component must not contribute to the reported distance.
    let point = start + (end - start) * foot + transverse * (0.4 * scale) + normal * (9.0 * scale);
    let interval = Extent1d::from_min_max(bounds[0], bounds[1]);
    let expected_parameters = if op == SolverOperation::ProjectedMinimize {
        vec![foot.clamp(bounds[0], bounds[1])]
    } else {
        let left = (foot - bounds[0]).abs();
        let right = (foot - bounds[1]).abs();
        if close(left, right) {
            bounds.to_vec()
        } else if left > right {
            vec![bounds[0]]
        } else {
            vec![bounds[1]]
        }
    };
    let expected_distance = scale * ((expected_parameters[0] - foot).powi(2) + 0.16).sqrt();
    let vectors = [normal];
    let answers = curve
        .global_point_solve::<Object>(
            &k,
            &interval,
            op,
            point,
            1e-10 * scale,
            request,
            &GlobalPointSolveOptions {
                opt_vectors: Some(&vectors),
                ..Default::default()
            },
        )
        .unwrap();
    let context = format!("scale={scale}, rotated={rotated}, reversed={reversed}, spline={spline}, bounds={bounds:?}, foot={foot}, op={op:?}, request={request:?}");
    assert!(
        !answers.is_empty(),
        "missing constrained extremum: {context}"
    );
    if request == SolutionRequested::Single {
        assert_eq!(answers.len(), 1, "Single keeps the best: {context}");
    }
    let mut parameters = Vec::new();
    for answer in answers.as_slice() {
        assert_eq!(answer.solution_type, SolutionType::SingleValue, "{context}");
        let parameter = answer.start.parameters[0];
        assert!(
            expected_parameters.iter().any(|&p| close(parameter, p)),
            "nonoptimal parameter {parameter}; expected {expected_parameters:?}: {context}"
        );
        assert!(
            close(answer.start.solution_value, expected_distance),
            "reported distance {} != {expected_distance}: {context}",
            answer.start.solution_value
        );
        let evaluated = curve.evaluate_point(&k, parameter).unwrap();
        let gap = evaluated - point;
        let projected = gap - normal * gap.dot(normal);
        assert!(
            close(projected.length(), expected_distance),
            "evaluated distance: {context}"
        );
        parameters.push(parameter);
    }
    // All means every equally optimal answer, without assuming an ordering or
    // imposing an arbitrary count on the solver's duplicate handling.
    if request == SolutionRequested::All {
        for expected in expected_parameters {
            assert!(
                parameters.iter().any(|&p| close(p, expected)),
                "missing optimal parameter {expected}: {context}"
            );
        }
    }
}

#[test]
fn projected_line_extrema_cover_boundaries_and_stationary_points() {
    for scale in [0.001, 0.01, 0.1, 10.0] {
        for rotated in [false, true] {
            for reversed in [false, true] {
                for spline in [false, true] {
                    for bounds in [[0.0, 1.0], [0.2, 0.8]] {
                        for foot in [-20.0, 20.0, 0.35, 0.5] {
                            for op in [
                                SolverOperation::ProjectedMinimize,
                                SolverOperation::ProjectedMaximize,
                            ] {
                                for request in [SolutionRequested::Single, SolutionRequested::All] {
                                    check_case(
                                        scale, rotated, reversed, spline, bounds, foot, op, request,
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
