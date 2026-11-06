//! Projected curve queries: endpoints, affine parameter maps, and objective bounds.

mod affine_spans {
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

    #[test]
    fn near_constant_affine_projection_certifies_candidates_without_dividing_by_zero() {
        let next = f64::from_bits(1.0_f64.to_bits() + 1);
        for isolated_upper in [false, true] {
            let mut points = vec![
                CPoint::new(1.0, 0.0, 0.0, 1.0),
                CPoint::new(next, 0.0, 1.0, 1.0),
            ];
            let mut knots = vec![0.0, 0.0, 1.0, 1.0];
            if isolated_upper {
                points.extend([
                    CPoint::new(10.0, 0.0, 0.0, 1.0),
                    CPoint::new(10.0, 1.0, 0.0, 1.0),
                ]);
                knots.extend([2.0, 2.0]);
            }
            let curve = Curve::from(BSplineCurve::from_nurb(
                NurbCurve::new(points, 1, knots).unwrap(),
                3,
                false,
            ));
            for operation in [
                SolverOperation::ProjectedMinimize,
                SolverOperation::ProjectedMaximize,
            ] {
                for requested in [SolutionRequested::Single, SolutionRequested::All] {
                    let answers = curve
                        .global_point_solve::<Object>(
                            &Kernel::default(),
                            &Extent1d::from_min_max(0.0, 1.0),
                            operation,
                            Vec3::ZERO,
                            1e-6,
                            requested,
                            &GlobalPointSolveOptions {
                                opt_vectors: Some(&[Vec3::Z]),
                                ..Default::default()
                            },
                        )
                        .unwrap();
                    assert!(!answers.is_empty());
                    let expected =
                        if isolated_upper && operation == SolverOperation::ProjectedMaximize {
                            10.0
                        } else {
                            1.0
                        };
                    let best = answers
                        .as_slice()
                        .iter()
                        .map(|s| s.start.solution_value)
                        .reduce(|a, b| {
                            if operation == SolverOperation::ProjectedMaximize {
                                a.max(b)
                            } else {
                                a.min(b)
                            }
                        })
                        .unwrap();
                    assert!(
                        (best - expected).abs() <= 1e-6,
                        "{operation:?}: {best} != {expected}"
                    );
                    for answer in answers.as_slice() {
                        assert_eq!(answer.solution_type, SolutionType::SingleValue);
                        let t = answer.start.parameters[0];
                        assert!((0.0..=1.0).contains(&t));
                        if isolated_upper
                            && operation == SolverOperation::ProjectedMaximize
                            && requested == SolutionRequested::Single
                        {
                            assert_eq!(t, 1.0);
                        }
                    }
                }
            }
        }
    }

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
        let best_representable = [
            f64::from_bits(rounded.to_bits() - 1),
            rounded,
            f64::from_bits(rounded.to_bits() + 1),
        ]
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
            matches!(&result, Err(error) if error.kind() == solidcore::error::ErrorKind::NotWithinTolerance),
            "an unrepresentable answer must return a tolerance refusal, not another error or partial success"
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

    #[test]
    fn rotated_constant_projection_retains_its_range_at_zero_tolerance() {
        let nurb = NurbCurve::new(
            vec![
                CPoint::new(0.0, 0.0, 0.0, 1.0),
                CPoint::new(1.0, 1.0, 0.0, 1.0),
            ],
            1,
            vec![0.0, 0.0, 1.0, 1.0],
        )
        .unwrap();
        let curve = Curve::from(BSplineCurve::from_nurb(nurb, 3, true));
        for operation in [
            SolverOperation::ProjectedMinimize,
            SolverOperation::ProjectedMaximize,
        ] {
            let answers = curve
                .global_point_solve::<Object>(
                    &Kernel::default(),
                    &Extent1d::from_min_max(0.0, 1.0),
                    operation,
                    Vec3::Z,
                    0.0,
                    SolutionRequested::All,
                    &GlobalPointSolveOptions {
                        opt_vectors: Some(&[Vec3::new(1.0, 1.0, 0.0)]),
                        ..Default::default()
                    },
                )
                .unwrap();
            assert_eq!(answers.len(), 1);
            let answer = &answers.as_slice()[0];
            assert_eq!(answer.solution_type, SolutionType::RangeOfValues);
            assert_eq!(answer.start.parameters[0], 0.0);
            assert_eq!(answer.end.parameters[0], 1.0);
            assert_eq!(answer.start.solution_value, 2.0);
            assert_eq!(answer.end.solution_value, 2.0);
        }
    }

    fn upper_jump_curve(start_parameter: f64, constant: bool, right_weight: f64) -> Curve {
        let knots = vec![start_parameter, start_parameter, 1.0, 1.0, 2.0, 2.0];
        let first = if constant {
            CPoint::new(1.0, 0.0, -1.0, 1.0)
        } else {
            CPoint::new(100.0, 0.0, 0.0, 1.0)
        };
        let nurb = NurbCurve::new(
            vec![
                first,
                CPoint::new(1.0, 0.0, 0.0, 1.0),
                CPoint::new(10.0 * right_weight, 0.0, 0.0, right_weight),
                CPoint::new(10.0, 1.0, 0.0, 1.0),
            ],
            1,
            knots.clone(),
        )
        .unwrap();
        let spline = BSplineCurve::from_nurb(nurb, 3, true);
        assert_eq!(spline.gw_nurb().unwrap().knots(), knots.as_slice());
        Curve::from(spline)
    }

    #[test]
    fn upper_discontinuity_and_exact_query_bound_keep_the_left_infimum() {
        let kernel = Kernel::default();
        for right_weight in [1.0, 2.0] {
            let curve = upper_jump_curve(0.0, false, right_weight);
            assert_eq!(curve.evaluate_point(&kernel, 1.0).unwrap().x, 10.0);
            assert!(
                (curve
                    .evaluate_point(&kernel, f64::from_bits(1.0_f64.to_bits() - 1))
                    .unwrap()
                    .x
                    - 1.0)
                    .abs()
                    < 1e-6
            );
            for upper in [1.0, f64::from_bits(1.0_f64.to_bits() - 1), 2.0] {
                for operation in [
                    SolverOperation::ProjectedMinimize,
                    SolverOperation::ProjectedMaximize,
                ] {
                    let answers = curve
                        .global_point_solve::<Object>(
                            &kernel,
                            &Extent1d::from_min_max(0.0, upper),
                            operation,
                            Vec3::ZERO,
                            EFF_ZERO_SQRT,
                            SolutionRequested::Single,
                            &GlobalPointSolveOptions {
                                opt_vectors: Some(&[Vec3::Z]),
                                ..Default::default()
                            },
                        )
                        .unwrap();
                    assert_eq!(answers.len(), 1);
                    let answer = &answers.as_slice()[0];
                    let t = answer.start.parameters[0];
                    assert!(t >= 0.0 && t <= upper);
                    let evaluated = curve.evaluate_point(&kernel, t).unwrap();
                    let actual = Vec3::Z.cross(evaluated).cross(Vec3::Z).length();
                    let expected = if operation == SolverOperation::ProjectedMinimize {
                        1.0
                    } else {
                        100.0
                    };
                    assert!((actual - expected).abs() <= EFF_ZERO_SQRT);
                    assert!(
                        (actual - answer.start.solution_value).abs()
                            <= 1e-9 * actual.max(answer.start.solution_value)
                    );
                }
            }
        }
    }

    #[test]
    fn constant_left_span_does_not_extend_range_over_right_valued_knot() {
        let curve = upper_jump_curve(0.0, true, 2.0);
        for operation in [
            SolverOperation::ProjectedMinimize,
            SolverOperation::ProjectedMaximize,
        ] {
            let answers = curve
                .global_point_solve::<Object>(
                    &Kernel::default(),
                    &Extent1d::from_min_max(0.0, 1.0),
                    operation,
                    Vec3::ZERO,
                    EFF_ZERO_SQRT,
                    SolutionRequested::All,
                    &GlobalPointSolveOptions {
                        opt_vectors: Some(&[Vec3::Z]),
                        ..Default::default()
                    },
                )
                .unwrap();
            assert_eq!(answers.len(), 1);
            let answer = &answers.as_slice()[0];
            if operation == SolverOperation::ProjectedMinimize {
                assert_eq!(answer.solution_type, SolutionType::RangeOfValues);
                assert_eq!(answer.start.parameters[0], 0.0);
                assert!(answer.end.parameters[0] < 1.0);
                assert_eq!(answer.start.solution_value, 1.0);
                assert_eq!(answer.end.solution_value, 1.0);
            } else {
                assert_eq!(answer.solution_type, SolutionType::SingleValue);
                assert_eq!(answer.start.parameters[0], 1.0);
                assert_eq!(answer.start.solution_value, 10.0);
            }
        }
    }

    #[test]
    fn unattainable_left_infimum_at_upper_jump_is_refused() {
        let lower = 1.0 - 1.5e-12;
        let curve = upper_jump_curve(lower, false, 1.0);
        // Distance decreases to the open-left limit 1. Even the final representable
        // parameter misses that limit by more than the caller's tolerance.
        let nearest = 1.0 + 99.0 * (1.0 - f64::from_bits(1.0_f64.to_bits() - 1)) / (1.0 - lower);
        assert!(nearest - 1.0 > EFF_ZERO_SQRT);
        let result = curve.global_point_solve::<Object>(
            &Kernel::default(),
            &Extent1d::from_min_max(lower, 1.0),
            SolverOperation::ProjectedMinimize,
            Vec3::ZERO,
            EFF_ZERO_SQRT,
            SolutionRequested::Single,
            &GlobalPointSolveOptions {
                opt_vectors: Some(&[Vec3::Z]),
                ..Default::default()
            },
        );
        assert!(
            matches!(&result, Err(error) if error.kind() == solidcore::error::ErrorKind::NotWithinTolerance)
        );
    }

    /// Reversed and infinite analytic lines retain their parameter maps in projected queries.
    #[test]
    fn projected_line_extrema_cover_reversed_and_unbounded_parameter_maps() {
        use solidcore::curve::line::Line;
        let kernel = Kernel::default();
        for bounded in [true, false] {
            for inside_out in [false, true] {
                let line = if bounded {
                    let xs = if inside_out { [7.0, -3.0] } else { [-3.0, 7.0] };
                    let control = xs
                        .into_iter()
                        .map(|x| CPoint {
                            x,
                            y: 2.0,
                            z: 3.0,
                            w: 1.0,
                        })
                        .collect();
                    let nurb = NurbCurve::new(control, 1, vec![10.0, 10.0, 20.0, 20.0]).unwrap();
                    Line::from_parts(
                        BSplineCurve::from_nurb(nurb, 3, true),
                        Vec3::new(1.0, 2.0, 3.0),
                        Vec3::X,
                        Extent1d::from_min_max(-2.0, 3.0),
                        2.0,
                        inside_out,
                    )
                } else {
                    Line::infinite(Vec3::new(1.0, 2.0, 3.0), Vec3::X * 2.0, 3, inside_out)
                };
                let (lo, hi) = if bounded { (12.0, 18.0) } else { (-1.0, 2.0) };
                let x_at = |t: f64| {
                    let step = if bounded {
                        let fraction = if inside_out {
                            (20.0 - t) / 10.0
                        } else {
                            (t - 10.0) / 10.0
                        };
                        -2.0 + 5.0 * fraction
                    } else if inside_out {
                        -t
                    } else {
                        t
                    };
                    1.0 + 2.0 * step
                };
                for t in [lo, (lo + hi) / 2.0, hi] {
                    let actual = line.evaluate_point(t);
                    assert!((actual.x - x_at(t)).abs() <= 1e-9 * actual.x.abs().max(x_at(t).abs()));
                    assert_eq!((actual.y, actual.z), (2.0, 3.0));
                }
                let max_x = x_at(lo).abs().max(x_at(hi).abs());
                let curve = Curve::from(line);
                for (operation, expected) in [
                    (SolverOperation::ProjectedMinimize, 2.0),
                    (
                        SolverOperation::ProjectedMaximize,
                        (max_x * max_x + 4.0).sqrt(),
                    ),
                ] {
                    let answers = curve
                        .global_point_solve::<Object>(
                            &kernel,
                            &Extent1d::from_min_max(lo, hi),
                            operation,
                            Vec3::new(0.0, 4.0, 9.0),
                            EFF_ZERO_SQRT,
                            SolutionRequested::Single,
                            &GlobalPointSolveOptions {
                                opt_vectors: Some(&[Vec3::Z]),
                                ..Default::default()
                            },
                        )
                        .unwrap();
                    assert_eq!(answers.len(), 1);
                    let answer = &answers.as_slice()[0];
                    let t = answer.start.parameters[0];
                    assert!(t >= lo && t <= hi);
                    let expected_at_t = (x_at(t) * x_at(t) + 4.0).sqrt();
                    for value in [answer.start.solution_value, expected_at_t] {
                        assert!(
                            (value - expected).abs() <= 1e-9 * value.abs().max(expected.abs()),
                            "bounded={bounded} inside_out={inside_out}: {value} versus {expected}"
                        );
                    }
                }
            }
        }
    }
}

mod endpoint_extrema {
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
        let point =
            start + (end - start) * foot + transverse * (0.4 * scale) + normal * (9.0 * scale);
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
                                    for request in
                                        [SolutionRequested::Single, SolutionRequested::All]
                                    {
                                        check_case(
                                            scale, rotated, reversed, spline, bounds, foot, op,
                                            request,
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
    fn projected_constant_distance_keeps_its_solution_range() {
        let k = Kernel::default();
        let curve = Curve::from(BSplineCurve::create_line_segment(3, Vec3::ZERO, Vec3::Z).unwrap());
        let interval = Extent1d::from_min_max(0.2, 0.8);
        let normal = [Vec3::Z];
        for op in [
            SolverOperation::ProjectedMinimize,
            SolverOperation::ProjectedMaximize,
        ] {
            let answers = curve
                .global_point_solve::<Object>(
                    &k,
                    &interval,
                    op,
                    Vec3::new(0.4, 0.0, 0.3),
                    1e-10,
                    SolutionRequested::All,
                    &GlobalPointSolveOptions {
                        opt_vectors: Some(&normal),
                        ..Default::default()
                    },
                )
                .unwrap();
            assert!(
                answers.as_slice().iter().any(|s| {
                    s.solution_type == SolutionType::RangeOfValues
                        && close(s.start.parameters[0], 0.2)
                        && close(s.end.parameters[0], 0.8)
                }),
                "the whole interval has projected distance 0.4: {op:?}"
            );
            for answer in answers.as_slice() {
                assert!(close(answer.start.solution_value, 0.4));
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
                    let context =
                        format!("parabola point={point:?}, op={op:?}, request={request:?}");
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
}

mod bounds {
    //! Endpoint candidates must not hide better interior projected extrema (#2840).
    //!
    //! Independent oracle for these rational quadratic fixtures: construct homogeneous
    //! Bernstein polynomials directly, then certify the stationary polynomial's sign
    //! using Bernstein coefficients and bisect its unique root. No kernel evaluator,
    //! optimizer, or C++ output supplies the oracle. Certificates use ordinary f64
    //! arithmetic with large sign margins, not formally outward-rounded intervals.
    #![allow(clippy::unwrap_used, clippy::indexing_slicing)]

    use solidcore::curve::curve::{p1::GlobalPointSolveOptions, Curve};
    use solidcore::geom::extent1d::Extent1d;
    use solidcore::geom::global_solver::{SolutionRequested, SolverOperation};
    use solidcore::geom::solution_array::{ObjectKind, SolutionObject, SolutionType};
    use solidcore::math::Vec3;
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

    fn eval(p: &[f64], t: f64) -> f64 {
        p.iter()
            .rev()
            .fold(0.0, |value, coefficient| value * t + coefficient)
    }

    fn choose(n: usize, k: usize) -> f64 {
        (0..k).fold(1.0, |v, i| v * (n - i) as f64 / (i + 1) as f64)
    }

    fn bernstein(power: &[f64]) -> Vec<f64> {
        let degree = power.len() - 1;
        (0..=degree)
            .map(|i| {
                (0..=i)
                    .map(|j| power[j] * choose(i, j) / choose(degree, j))
                    .sum()
            })
            .collect()
    }

    struct Span {
        bounds: [f64; 2],
        homogeneous: [[f64; 3]; 4],
    }

    impl Span {
        fn new(control: &[[f64; 4]], original: [f64; 2], bounds: [f64; 2]) -> Self {
            assert_eq!(control.len(), 3);
            assert!(control.iter().all(|p| p[3] > 0.0));
            let alpha = (bounds[0] - original[0]) / (original[1] - original[0]);
            let beta = (bounds[1] - bounds[0]) / (original[1] - original[0]);
            let homogeneous = std::array::from_fn(|axis| {
                let a = control[0][axis];
                let b = 2.0 * (control[1][axis] - a);
                let c = a - 2.0 * control[1][axis] + control[2][axis];
                [
                    a + b * alpha + c * alpha * alpha,
                    beta * (b + 2.0 * c * alpha),
                    c * beta * beta,
                ]
            });
            Self {
                bounds,
                homogeneous,
            }
        }

        fn point(&self, s: f64) -> Vec3 {
            let w = eval(&self.homogeneous[3], s);
            Vec3::new(
                eval(&self.homogeneous[0], s) / w,
                eval(&self.homogeneous[1], s) / w,
                eval(&self.homogeneous[2], s) / w,
            )
        }

        fn distance(&self, s: f64, target: Vec3, normal: Vec3) -> f64 {
            normal.cross(target - self.point(s)).cross(normal).length()
        }

        fn stationary(&self, target: Vec3, normal: Vec3) -> [f64; 5] {
            let w = self.homogeneous[3];
            let projected: [Vec3; 3] = std::array::from_fn(|i| {
                let xyz = Vec3::new(
                    self.homogeneous[0][i],
                    self.homogeneous[1][i],
                    self.homogeneous[2][i],
                );
                normal.cross(target * w[i] - xyz).cross(normal)
            });
            let mut numerator = [0.0; 5];
            for i in 0..3 {
                for j in 0..3 {
                    numerator[i + j] += projected[i].dot(projected[j]);
                }
            }
            // Distance squared is N/W²; its derivative has numerator H=N'W-2NW'.
            let mut h = [0.0; 6];
            for i in 1..5 {
                for j in 0..3 {
                    h[i - 1 + j] += i as f64 * numerator[i] * w[j];
                }
            }
            for i in 0..5 {
                for j in 1..3 {
                    h[i + j - 1] -= 2.0 * numerator[i] * j as f64 * w[j];
                }
            }
            let magnitude = h.iter().fold(1.0_f64, |a, b| a.max(b.abs()));
            assert!(
                h[5].abs() <= 1e-12 * magnitude,
                "degree-five terms cancel algebraically"
            );
            [h[0], h[1], h[2], h[3], h[4]]
        }

        fn extremum(
            &self,
            target: Vec3,
            normal: Vec3,
            maximize: bool,
            monotone: bool,
        ) -> (f64, f64) {
            let h = self.stationary(target, normal);
            let s = if monotone {
                // The first span of gen0002 increases all the way to the join.
                assert!(bernstein(&h).iter().all(|&v| v > 1.0));
                assert!(maximize);
                1.0
            } else {
                let derivative = [h[1], 2.0 * h[2], 3.0 * h[3], 4.0 * h[4]];
                let direction = if maximize { -1.0 } else { 1.0 };
                assert!(
                    bernstein(&derivative).iter().all(|&v| direction * v > 0.1),
                    "H must be strictly monotone"
                );
                assert!(direction * eval(&h, 0.0) < -0.01 && direction * eval(&h, 1.0) > 0.01);
                let (mut lo, mut hi) = (0.0, 1.0);
                for _ in 0..64 {
                    let mid = (lo + hi) / 2.0;
                    if direction * eval(&h, mid) < 0.0 {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                (lo + hi) / 2.0
            };
            (
                self.bounds[0] + s * (self.bounds[1] - self.bounds[0]),
                self.distance(s, target, normal),
            )
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn check(
        control: &[[f64; 4]],
        knots: &[f64],
        bounds: [f64; 2],
        target: Vec3,
        normal: Vec3,
        tolerance: f64,
        maximize: bool,
        spans: &[Span],
    ) {
        let cp = control
            .iter()
            .map(|p| CPoint {
                x: p[0],
                y: p[1],
                z: p[2],
                w: p[3],
            })
            .collect();
        let curve = Curve::from(BSplineCurve::from_nurb(
            NurbCurve::new(cp, 2, knots.to_vec()).unwrap(),
            3,
            false,
        ));
        let candidates: Vec<_> = spans
            .iter()
            .enumerate()
            .map(|(i, s)| s.extremum(target, normal, maximize, spans.len() == 2 && i == 0))
            .collect();
        let expected = candidates
            .iter()
            .copied()
            .reduce(|a, b| {
                if (maximize && b.1 > a.1) || (!maximize && b.1 < a.1) {
                    b
                } else {
                    a
                }
            })
            .unwrap();
        assert!(expected.0 > bounds[0] && expected.0 < bounds[1]);
        let k = Kernel::default();
        for request in [SolutionRequested::Single, SolutionRequested::All] {
            let answers = curve
                .global_point_solve::<Object>(
                    &k,
                    &Extent1d::from_min_max(bounds[0], bounds[1]),
                    if maximize {
                        SolverOperation::ProjectedMaximize
                    } else {
                        SolverOperation::ProjectedMinimize
                    },
                    target,
                    tolerance,
                    request,
                    &GlobalPointSolveOptions {
                        opt_vectors: Some(&[normal]),
                        ..Default::default()
                    },
                )
                .unwrap();
            assert!(!answers.is_empty());
            if request == SolutionRequested::Single {
                assert_eq!(answers.len(), 1);
            }
            for answer in answers.as_slice() {
                assert_eq!(answer.solution_type, SolutionType::SingleValue);
                let t = answer.start.parameters[0];
                let span = spans
                    .iter()
                    .find(|s| t >= s.bounds[0] && t <= s.bounds[1])
                    .unwrap();
                let actual = span.distance(
                    (t - span.bounds[0]) / (span.bounds[1] - span.bounds[0]),
                    target,
                    normal,
                );
                assert!(
                    (actual - answer.start.solution_value).abs()
                        <= 1e-12_f64.max(1e-9 * actual.abs()),
                    "reported metric disagrees with rational polynomial"
                );
                // Preserve the caller's original tolerance for global optimality.
                assert!((actual-expected.1).abs()<=tolerance,
                    "endpoint displaced interior optimum: t={t}, distance={actual}, oracle={expected:?}, tol={tolerance}");
            }
        }
    }

    // Original generated inputs, retained exactly so the native comparison and the
    // independent mathematical regressions exercise the same geometry.

    #[test]
    fn gen0002_interior_extremum_survives_endpoint_candidates() {
        let control = [
            [
                5.477596627102355,
                -4.556269464099508,
                1.114995725656613,
                1.0,
            ],
            [
                9.076471354595526,
                -2.132022809946599,
                1.2093799310872764,
                0.710088106679831,
            ],
            [11.18257991594583, 3.841957587213183, 4.293909950576297, 1.0],
            [
                6.804762646024233,
                7.58827958804344,
                4.8887288430295435,
                0.710088106679831,
            ],
            [2.251369955507407, 9.248272788188828, 6.340319531047773, 1.0],
            [
                -3.6074205877400116,
                5.545897440403772,
                4.115642140064191,
                0.710088106679831,
            ],
            [
                -3.604555041788431,
                0.9414150425449053,
                3.1959906041602584,
                1.0,
            ],
            [
                -1.5116827423535617,
                -4.208922190082522,
                0.42322769408518396,
                0.710088106679831,
            ],
            [
                5.227686992962869,
                -4.605290016641698,
                1.096440365002791,
                1.0,
            ],
        ];
        let knots = [
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ];
        let bounds = [0.6256821477105545, 0.9635580665549801];
        let spans = [
            Span::new(&control[4..7], [0.5, 0.75], [bounds[0], 0.75]),
            Span::new(&control[6..9], [0.75, 1.0], [0.75, bounds[1]]),
        ];
        check(
            &control,
            &knots,
            bounds,
            Vec3::new(3.864483291304881, 2.34600166204466, 3.727657628352193),
            Vec3::new(0.6184107265151202, -0.5901631565932342, -0.695866229492025),
            0.05,
            true,
            &spans,
        );
    }

    #[test]
    fn gen2419_interior_extremum_survives_endpoint_candidates() {
        let control = [
            [
                4.68246458913153,
                0.2683108077599766,
                0.7857380281272487,
                1.0,
            ],
            [
                4.405156176361579,
                0.5572137913733559,
                0.660330035392789,
                0.8275822849855569,
            ],
            [5.18329866623959, 1.417855376481916, 0.820266301512733, 1.0],
            [
                4.1740561315767195,
                1.789570193122566,
                0.6973456848123297,
                0.8275822849855569,
            ],
            [
                4.2999559828848914,
                2.308063461512357,
                0.8470050195809443,
                1.0,
            ],
            [
                2.94307866372967,
                2.030654673617574,
                0.7045870141857388,
                0.8275822849855569,
            ],
            [
                3.145828375026399,
                1.8168898669671454,
                0.832251893331091,
                1.0,
            ],
        ];
        let knots = [
            0.0,
            0.0,
            0.0,
            0.3333333333333333,
            0.3333333333333333,
            0.6666666666666666,
            0.6666666666666666,
            1.0,
            1.0,
            1.0,
        ];
        let bounds = [0.18190734190117286, 0.27563504714764536];
        let spans = [Span::new(&control[0..3], [0.0, 1.0 / 3.0], bounds)];
        check(
            &control,
            &knots,
            bounds,
            Vec3::new(7.695548237243017, 1.1786104807110709, 0.9369883409089445),
            Vec3::new(-0.2870092592131892, -0.84405099161452, -0.2736338163951917),
            0.001,
            false,
            &spans,
        );
    }
}
