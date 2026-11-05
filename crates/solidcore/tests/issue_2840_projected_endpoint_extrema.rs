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

#[test]
fn projected_maximum_equation_has_its_actual_derivative() {
    use solidcore::curve::curve::{p4::CurvePropertyEfo, CurvePropertyType};
    use solidcore::geom::local_solve1d::{EvalFunction, FunctionValue};

    let k = Kernel::default();
    let curve = Curve::from(BSplineCurve::create_line_segment(3, Vec3::ZERO, Vec3::X).unwrap());
    // C(t)=(t,0,0), P=0, N=(0,0,s): the maximum's signed stationarity
    // equation is F(t)=-s^2*t, hence F'(t)=-s^2 exactly.
    for scale in [0.5, 1.0, 2.0] {
        let vectors = [Vec3::ZERO, Vec3::Z * scale];
        let mut evaluator = CurvePropertyEfo::new(
            &k,
            &curve,
            CurvePropertyType::ProjectedPointMaximize,
            None,
            Some(&vectors),
            1e-20,
        );
        for t in [0.2, 0.25, 0.3] {
            let mut value = FunctionValue::default();
            evaluator.evaluate(t, &mut value, false).unwrap();
            assert!(close(value.f, -scale * scale * t));
            assert!(close(value.f_prime, -scale * scale));
        }
    }
}

#[test]
fn projected_parabola_extrema_include_constrained_endpoints() {
    use solidcore::nurbs::crv_basic::Curve as NurbCurve;
    use solidcore::nurbs::math::CPoint;

    // Quadratic Bernstein control points give exactly C(t) = (t, t^2, 0).
    let control_points = vec![
        CPoint::from_point(Vec3::ZERO),
        CPoint::from_point(Vec3::new(0.5, 0.0, 0.0)),
        CPoint::from_point(Vec3::new(1.0, 1.0, 0.0)),
    ];
    let nurb = NurbCurve::new(control_points, 2, vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]).unwrap();
    let curve = Curve::from(BSplineCurve::from_nurb(nurb, 3, false));
    let k = Kernel::default();
    let interval = Extent1d::from_min_max(0.0, 1.0);
    let normal = [Vec3::Z];

    // For P=(a,b,z), projected squared distance is
    // f(t)=(t-a)^2+(t^2-b)^2, f'(t)=2(t-a)+4t(t^2-b).
    // (-2,-2): f'>0 throughout [0,1]; (3,3): f'<0 throughout.
    // (3/4,0): f' vanishes at 1/2 and f''=2+12t^2>0;
    // comparing the endpoint values gives the unique maximum at t=1.
    // z=7 deliberately exercises removal of the projection-normal component.
    for (point, minimum, maximum) in [
        (Vec3::new(-2.0, -2.0, 7.0), 0.0, 1.0),
        (Vec3::new(3.0, 3.0, 7.0), 1.0, 0.0),
        (Vec3::new(0.75, 0.0, 7.0), 0.5, 1.0),
    ] {
        for (op, expected_parameter) in [
            (SolverOperation::ProjectedMinimize, minimum),
            (SolverOperation::ProjectedMaximize, maximum),
        ] {
            let expected_point = Vec3::new(expected_parameter, expected_parameter.powi(2), 0.0);
            let projected_gap =
                Vec3::new(expected_point.x - point.x, expected_point.y - point.y, 0.0);
            let expected_distance = projected_gap.length();
            for request in [SolutionRequested::Single, SolutionRequested::All] {
                let answers = curve
                    .global_point_solve::<Object>(
                        &k,
                        &interval,
                        op,
                        point,
                        1e-10,
                        request,
                        &GlobalPointSolveOptions {
                            opt_vectors: Some(&normal),
                            ..Default::default()
                        },
                    )
                    .unwrap();
                let context = format!("parabola point={point:?}, op={op:?}, request={request:?}");
                assert!(!answers.is_empty(), "missing curved extremum: {context}");
                if request == SolutionRequested::Single {
                    assert_eq!(answers.len(), 1, "{context}");
                }
                for answer in answers.as_slice() {
                    assert_eq!(answer.solution_type, SolutionType::SingleValue, "{context}");
                    let parameter = answer.start.parameters[0];
                    assert!(
                        close(parameter, expected_parameter),
                        "parameter {parameter} != {expected_parameter}: {context}"
                    );
                    assert!(
                        close(answer.start.solution_value, expected_distance),
                        "distance {} != {expected_distance}: {context}",
                        answer.start.solution_value
                    );
                    let evaluated = curve.evaluate_point(&k, parameter).unwrap();
                    assert!(
                        close(evaluated.x, expected_point.x)
                            && close(evaluated.y, expected_point.y)
                            && close(evaluated.z, 0.0),
                        "curve evaluation: {context}"
                    );
                }
            }
        }
    }
}
