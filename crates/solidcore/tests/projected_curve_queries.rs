//! Projected curve queries: endpoints, affine parameter maps, and objective bounds.

mod affine_spans {
    //! Projected extrema of affine spans must not depend on their knot widths.
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
    fn ordinary_long_line_reference_minimum_is_representable() {
        let curve = Curve::from(
            BSplineCurve::create_line_segment(
                3,
                Vec3::new(-100.0, 0.5, 100.0),
                Vec3::new(100.0, 0.5, 100.0),
            )
            .unwrap(),
        );
        for requested in [SolutionRequested::Single, SolutionRequested::All] {
            let answers = curve
                .global_point_solve::<Object>(
                    &Kernel::default(),
                    &curve.natural_interval(),
                    SolverOperation::ProjectedMinimize,
                    Vec3::new(0.0, 10.0, 100.0),
                    EFF_ZERO_SQRT,
                    requested,
                    &GlobalPointSolveOptions {
                        opt_vectors: Some(&[-Vec3::Z]),
                        ..Default::default()
                    },
                )
                .unwrap();
            assert_eq!(answers.len(), 1);
            let native = curve
                .evaluate_point(
                    &Kernel::default(),
                    answers.as_slice()[0].start.parameters[0],
                )
                .unwrap();
            // Canonical ordering may retain a neighbouring parameter whose
            // objective rounds to the same minimum. Grade the native metric,
            // rather than requiring a particular representative of that tie.
            assert_eq!(
                (native.x * native.x + (native.y - 10.0) * (native.y - 10.0)).sqrt(),
                9.5
            );
            assert_eq!(answers.as_slice()[0].start.solution_value, 9.5);
        }
    }

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

    #[test]
    fn collapsed_analytic_line_has_constant_projected_extrema() {
        use solidcore::curve::line::Line;
        let mut line =
            Line::create_line_segment(3, Vec3::new(2.0, 3.0, 0.0), Vec3::new(5.0, 3.0, 0.0), None)
                .unwrap();
        line.set_analytic_domain(Extent1d::from_min_max(0.0, 0.0));
        let curve = Curve::from(line);
        for op in [
            SolverOperation::ProjectedMinimize,
            SolverOperation::ProjectedMaximize,
        ] {
            let answers = curve
                .global_point_solve::<Object>(
                    &Kernel::default(),
                    &curve.natural_interval(),
                    op,
                    Vec3::ZERO,
                    0.0,
                    SolutionRequested::All,
                    &GlobalPointSolveOptions {
                        opt_vectors: Some(&[Vec3::Z]),
                        ..Default::default()
                    },
                )
                .unwrap();
            assert_eq!(answers.len(), 1);
            assert_eq!(
                answers.as_slice()[0].solution_type,
                SolutionType::RangeOfValues
            );
            assert!((answers.as_slice()[0].start.solution_value - 13.0_f64.sqrt()).abs() < 1e-14);
        }
    }

    #[test]
    fn mapped_line_parallelism_survives_endpoint_rounding() {
        use solidcore::curve::line::Line;
        let line =
            Line::create_line_segment(3, Vec3::new(0.1, 0.2, 0.3), Vec3::new(1.7, 2.9, 4.2), None)
                .unwrap();
        let normal = line.line_vector();
        let curve = Curve::from(line);
        for op in [
            SolverOperation::ProjectedMinimize,
            SolverOperation::ProjectedMaximize,
        ] {
            let answers = curve
                .global_point_solve::<Object>(
                    &Kernel::default(),
                    &curve.natural_interval(),
                    op,
                    Vec3::ZERO,
                    0.0,
                    SolutionRequested::All,
                    &GlobalPointSolveOptions {
                        opt_vectors: Some(&[normal]),
                        ..Default::default()
                    },
                )
                .unwrap();
            assert_eq!(answers.len(), 1);
            assert_eq!(
                answers.as_slice()[0].solution_type,
                SolutionType::RangeOfValues
            );
        }
    }

    #[test]
    fn mapped_line_does_not_erase_a_small_transverse_direction() {
        use solidcore::curve::line::Line;
        let line = Line::from_parts(
            BSplineCurve::create_line_segment(3, Vec3::X, Vec3::new(1.0, 0.0, 1.0)).unwrap(),
            Vec3::X,
            Vec3::new(1e-16, 0.0, 1.0),
            Extent1d::from_min_max(0.0, 1.0),
            1.0,
            false,
        );
        let curve = Curve::from(line);
        for op in [
            SolverOperation::ProjectedMinimize,
            SolverOperation::ProjectedMaximize,
        ] {
            let answers = curve
                .global_point_solve::<Object>(
                    &Kernel::default(),
                    &curve.natural_interval(),
                    op,
                    Vec3::ZERO,
                    0.0,
                    SolutionRequested::All,
                    &GlobalPointSolveOptions {
                        opt_vectors: Some(&[Vec3::Z]),
                        ..Default::default()
                    },
                )
                .unwrap();
            assert!(!answers.is_empty());
            assert!(answers
                .as_slice()
                .iter()
                .all(|s| s.solution_type == SolutionType::SingleValue));
        }
    }

    #[test]
    fn zero_retention_slack_accepts_affine_minimum_to_rounding() {
        let points = [
            (0.1569346301508998, -0.5329538213524829, -0.8439542626476992),
            (0.4224662792522742, -0.4068072031595056, -0.7199688610446469),
            (0.6052455395473925, -0.3885960084138218, -0.5099311278223793),
            (0.7831613767269248, -0.01515768701313519, -0.701567737379753),
        ]
        .into_iter()
        .map(|(x, y, z)| CPoint::new(x, y, z, solidcore::nurbs::math::NOW))
        .collect();
        let knots = vec![
            0.0,
            0.0,
            0.2628985868293973,
            0.8036901450210183,
            1.3482261889778533,
            1.3482261889778533,
        ];
        let curve = Curve::from(BSplineCurve::from_nurb(
            NurbCurve::new(points, 1, knots).unwrap(),
            3,
            false,
        ));
        let answers = curve
            .global_point_solve::<Object>(
                &Kernel::default(),
                &curve.natural_interval(),
                SolverOperation::ProjectedMinimize,
                Vec3::new(0.5501272357521876, -0.3940877155300753, -0.5732694057512258),
                0.0,
                SolutionRequested::All,
                &GlobalPointSolveOptions {
                    opt_vectors: Some(&[Vec3::new(
                        0.7405015371827215,
                        0.6291959064853571,
                        0.23615669520912666,
                    )]),
                    ..Default::default()
                },
            )
            .unwrap();
        let best = answers
            .as_slice()
            .iter()
            .min_by(|a, b| a.start.solution_value.total_cmp(&b.start.solution_value))
            .unwrap();
        assert!((best.start.parameters[0] - 0.6406108473689742).abs() < 2e-15);
        assert!(best.start.solution_value < 2e-15);
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
        let solutions = curve.global_point_solve::<Object>(
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
        );
        assert!(solutions.is_ok(), "gap={gap}, parameter_scale={parameter_scale}, point={point:?}, normal_scale={normal_scale}: {:?}", solutions.as_ref().err());
        let solutions = solutions.unwrap();
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
        for y in [0.3, 0.2] {
            check(
                1e-8,
                1.0,
                Vec3::new(0.0, y, 0.0),
                1.0,
                SolverOperation::ProjectedMinimize,
                (1.0 - 2.0 * y) / 5.0_f64.sqrt(),
            );
        }
    }

    #[test]
    fn narrower_raw_spans_refuse_unrepresentable_minima_even_with_retention_slack() {
        // Exact rational parameter arithmetic gives nearest-point objective
        // excesses from 8.14e-10 to 2.76e-8 (scale1), and 8.98e-11 to1.09e-8
        // (scale1e6). These are geometric errors, not arithmetic enclosure widths.
        for (gap, scale) in [(1.5e-12, 1.0), (1.9e-12, 1.0), (1.5e-12, 1e6)] {
            let (curve, knots) = raw_polyline(gap, scale);
            for y in [0.2, 0.3] {
                for normal_scale in [1.0, 0.5, -2.0] {
                    let result = curve.global_point_solve::<Object>(
                        &Kernel::default(),
                        &Extent1d::from_min_max(knots[0], knots[5]),
                        SolverOperation::ProjectedMinimize,
                        Vec3::new(0.0, y, 0.0),
                        EFF_ZERO_SQRT,
                        SolutionRequested::Single,
                        &GlobalPointSolveOptions {
                            opt_vectors: Some(&[Vec3::Z * normal_scale]),
                            ..Default::default()
                        },
                    );
                    assert!(
                        matches!(result,Err(error)if error.kind()==solidcore::error::ErrorKind::NotWithinTolerance),
                        "gap={gap},scale={scale},y={y},normal={normal_scale} must refuse the geometric gap"
                    );
                }
            }
        }
    }

    #[test]
    fn retention_slack_cannot_make_an_unrepresentable_minimum_accurate() {
        let (curve, knots) = raw_polyline(1.5e-12, 1.0);
        for slack in [1e-12, 1e-4, 1.0] {
            let result = curve.global_point_solve::<Object>(
                &Kernel::default(),
                &Extent1d::from_min_max(knots[0], knots[5]),
                SolverOperation::ProjectedMinimize,
                Vec3::new(0.0, 0.3, 0.0),
                slack,
                SolutionRequested::Single,
                &GlobalPointSolveOptions {
                    opt_vectors: Some(&[Vec3::Z]),
                    ..Default::default()
                },
            );
            assert!(
                matches!(result, Err(error) if error.kind() == solidcore::error::ErrorKind::NotWithinTolerance),
                "retention slack {slack} must not loosen the solve certificate"
            );
        }
    }

    #[test]
    fn affine_minimum_obeys_nonunit_projection_metric() {
        for normal_scale in [0.5, -2.0] {
            check(
                1e-8,
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
        for (gap, scale) in [(1e-8, 1e-3), (1e-8, 1e6)] {
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
                        .unwrap_or_else(|error| panic!("right_weight={right_weight},upper={upper},operation={operation:?}: {error}"));
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
                        .unwrap_or_else(|error| panic!("bounded={bounded},inside_out={inside_out},operation={operation:?}: {error}"));
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
    //! Constrained projected extrema include finite curve endpoints.
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
        let context = format!(
            "scale={scale}, rotated={rotated}, reversed={reversed}, spline={spline}, bounds={bounds:?}, foot={foot}, op={op:?}, request={request:?}"
        );
        let answers = curve.global_point_solve::<Object>(
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
        );
        assert!(answers.is_ok(), "{context}: {:?}", answers.as_ref().err());
        let answers = answers.unwrap();
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
    //! Endpoint candidates must not hide better interior projected extrema.
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
                assert!(
                    (actual - expected.1).abs() <= tolerance,
                    "endpoint displaced interior optimum: t={t}, distance={actual}, oracle={expected:?}, tol={tolerance}"
                );
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

mod missed_cutter_keeps_geometry {
    //! Missed projected cutters must preserve the input geometry, not only its counts.
    #![allow(clippy::unwrap_used, clippy::expect_used)]
    use solidcore::blend::fillet_executive::api::circular_fillet;
    use solidcore::curve::curve::Curve;
    use solidcore::math::Vec3;
    use solidcore::primitive::{api_curves, api_primitives, api_trim_ops};
    use solidcore::topology::brep::p1;
    use solidcore::topology::brep_data_write::DATABASE_VERSION;
    use solidcore::topology::database_io::{DatabaseWriter, FileType};
    use solidcore::topology::geometry_write::{write_curve, write_surface};
    use solidcore::topology::Brep;
    use solidcore::Kernel;
    use std::collections::BTreeMap;

    type Geometry = BTreeMap<String, Vec<u8>>;

    // Binary geometry records contain complete definitions (including placement and
    // NURBS control points/weights/knots), not sampled points or digest aggregates.
    // Typed arena identities and oriented topology are kept too: a missed cutter
    // has no new boundaries to introduce. Deliberately omit query memo/cache data.
    fn geometry(k: &Kernel, b: &Brep) -> Geometry {
        let mut out = BTreeMap::new();
        let mut put = |name: String, bytes: Vec<u8>| {
            out.insert(name, bytes);
        };
        for id in p1::vertices(b, None).unwrap() {
            let v = b.vertex(id).unwrap();
            let bytes = [v.point.x, v.point.y, v.point.z, v.tolerance]
                .into_iter()
                .flat_map(f64::to_le_bytes)
                .collect();
            put(format!("vertex/{id:?}"), bytes);
        }
        for id in p1::edges(b, None).unwrap() {
            let e = b.edge(id).unwrap();
            put(
                format!("edge/{id:?}"),
                format!(
                    "{:?}/{:?}/{:?}/{:?}",
                    e.vertices, e.edgeuses, e.interval, e.tolerance
                )
                .into_bytes(),
            );
            if let Some(curve) = &e.curve {
                let mut db = DatabaseWriter::open_for_write(FileType::Binary);
                db.write_type(curve.type_code(), Some(curve.dim()));
                write_curve(curve, &mut db, DATABASE_VERSION).unwrap();
                put(format!("edge-curve/{id:?}"), db.into_bytes());
            }
        }
        for id in p1::faces(b, None).unwrap() {
            let f = b.face(id).unwrap();
            put(
                format!("face/{id:?}"),
                format!(
                    "{:?}/{:?}/{:?}/{:?}/{:?}",
                    f.shell, f.loops, f.faceuse, f.uv_domain, f.tolerance
                )
                .into_bytes(),
            );
            if let Some(surface) = &f.surface {
                assert!(
                    surface.procedural_definition().is_none(),
                    "fixture has a retained definition outside this geometry record"
                );
                // Ordinary Clone starts fresh NURB/cache slots; evaluation copies
                // would share them and could prime the operation before it runs.
                let surface = surface.clone();
                surface.kind_with_nurb(k).unwrap();
                let mut db = DatabaseWriter::open_for_write(FileType::Binary);
                db.write_type(surface.type_code(), None);
                write_surface(&surface, &mut db, DATABASE_VERSION).unwrap();
                put(format!("face-surface/{id:?}"), db.into_bytes());
            }
        }
        for (id, eu) in b.edgeuses().iter() {
            put(
                format!("edgeuse/{id:?}"),
                format!(
                    "{:?}/{:?}/{:?}/{:?}/{:?}/{:?}",
                    eu.edge, eu.owner, eu.ccw, eu.cw, eu.vertexuse, eu.orientation
                )
                .into_bytes(),
            );
            if let Some(curve) = &eu.uv_trim_curve {
                let mut db = DatabaseWriter::open_for_write(FileType::Binary);
                db.write_type(curve.type_code(), Some(curve.dim()));
                write_curve(curve, &mut db, DATABASE_VERSION).unwrap();
                put(format!("uv-trim/{id:?}"), db.into_bytes());
            }
        }
        // These records contain only topology and orientation (no geometry cache).
        for (id, item) in b.loops().iter() {
            put(format!("loop/{id:?}"), format!("{item:?}").into_bytes());
        }
        for (id, item) in b.loopuses().iter() {
            put(format!("loopuse/{id:?}"), format!("{item:?}").into_bytes());
        }
        for (id, item) in b.faceuses().iter() {
            put(format!("faceuse/{id:?}"), format!("{item:?}").into_bytes());
        }
        for (id, item) in b.shells().iter() {
            put(format!("shell/{id:?}"), format!("{item:?}").into_bytes());
        }
        for (id, item) in b.regions().iter() {
            put(format!("region/{id:?}"), format!("{item:?}").into_bytes());
        }
        for (id, item) in b.vertexuses().iter() {
            put(
                format!("vertexuse/{id:?}"),
                format!("{item:?}").into_bytes(),
            );
        }
        out
    }

    fn input_box(k: &mut Kernel, origin: Vec3, filleted: bool) -> Brep {
        let mut b = api_primitives::create_box(k, origin, 10.0, 5.0, 15.0).unwrap();
        if filleted {
            circular_fillet(k, &mut b, 0.5).unwrap();
        }
        b
    }

    fn counts(b: &Brep) -> (usize, usize, usize) {
        (
            p1::faces(b, None).unwrap().len(),
            p1::edges(b, None).unwrap().len(),
            p1::vertices(b, None).unwrap().len(),
        )
    }

    fn preserves_missed_box(filleted: bool) {
        let mut k = Kernel::default();
        let mut b = input_box(&mut k, Vec3::ZERO, filleted);
        let before = geometry(&k, &b);
        let cutter = Curve::from(
            api_curves::create_curve(&[
                Vec3::new(100.0, -2.0, 20.0),
                Vec3::new(101.0, 5.0, 20.0),
                Vec3::new(102.0, 0.0, 20.0),
                Vec3::new(103.0, 2.0, 20.0),
            ])
            .unwrap(),
        );
        // For every projected box point, half the distance-squared derivative to
        // this Bezier cutter is >=123>0. Its first endpoint is nearest, and
        // 3*(qy+2)-21*(qx-100)>0 throughout the box and at the reference. Keeping
        // that side must retain every input surface and boundary exactly.
        api_trim_ops::project_and_trim(
            &mut k,
            &mut b,
            &cutter,
            -Vec3::Z,
            Vec3::new(5.0, 1.0, 20.0),
        )
        .unwrap();
        let after = geometry(&k, &b);
        let changed: Vec<_> = before
            .keys()
            .chain(after.keys())
            .filter(|key| before.get(*key) != after.get(*key))
            .collect();
        assert!(changed.is_empty(), "changed geometry records: {changed:?}");
        // A valid translated body has exactly the same topology counts. The
        // definition comparison must reject it, so counts cannot make this pass.
        let translated = input_box(&mut k, Vec3::new(1.0, 2.0, 3.0), filleted);
        assert_eq!(counts(&b), counts(&translated));
        let translated_geometry = geometry(&k, &translated);
        assert_ne!(after, translated_geometry);
        assert!(
            after.iter().any(|(name, bytes)| {
                name.starts_with("face-surface/")
                    && translated_geometry
                        .get(name)
                        .is_some_and(|other| bytes != other)
            }),
            "negative control must change a complete surface definition, not just IDs"
        );
        println!(
            "filleted={filleted}: {} geometry/topology records unchanged; translated control rejected",
            after.len()
        );
    }

    #[test]
    fn pat_plain_misses_retains_complete_input_geometry() {
        preserves_missed_box(false);
    }
    #[test]
    fn pat_misses_retains_complete_filleted_input_geometry() {
        preserves_missed_box(true);
    }
}

mod g1_joins {
    #![allow(clippy::unwrap_used, clippy::indexing_slicing)]
    // Approximate cache continuity must not erase actual projected corner extrema.
    //
    // The two segments meet at an angle 2 atan(1/1024), below the cache's 0.5-degree
    // G1 threshold. Their shared point is nevertheless a genuine minimum or maximum
    // of the projected distance, depending on the query point.

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

    #[test]
    fn projected_single_preserves_own_enclosure_across_affine_leaves() {
        for maximize in [false, true] {
            let points = if maximize {
                [
                    (1e8, 1e7, 0.0),
                    (1e8, 1e7, 0.0),
                    (0.0, 1e-10, 0.0),
                    (0.0, 1e-10, 0.0),
                ]
            } else {
                [
                    (1e8, 1e8, 0.0),
                    (-2e8, -2e8, 0.0),
                    (1e-9, 0.0, 0.0),
                    (1e-9, 0.0, 1.0),
                ]
            };
            let curve = Curve::from(BSplineCurve::from_nurb(
                NurbCurve::new(
                    points
                        .into_iter()
                        .map(|(x, y, z)| CPoint::new(x, y, z, NOW))
                        .collect(),
                    1,
                    vec![0.0, 0.0, 1.0, 1.0, 2.0, 2.0],
                )
                .unwrap(),
                3,
                false,
            ));
            let normal = if maximize {
                Vec3::new(1.0, 0.1, 0.0)
            } else {
                Vec3::Z
            };
            let operation = if maximize {
                SolverOperation::ProjectedMaximize
            } else {
                SolverOperation::ProjectedMinimize
            };
            let options = GlobalPointSolveOptions {
                opt_vectors: Some(&[normal]),
                ..Default::default()
            };
            let interval = Extent1d::from_min_max(0.0, 2.0);
            let answers = curve
                .global_point_solve::<Handle>(
                    &Kernel::default(),
                    &interval,
                    operation,
                    Vec3::ZERO,
                    0.0,
                    SolutionRequested::Single,
                    &options,
                )
                .unwrap();
            assert_eq!(answers.len(), 1);
            // Minimum: the first span crosses the origin at exact t=1/3.
            // Maximum: exact binary64 .1 gives |1e8*.1-1e7|=1e8/(5*2^55),
            // so its true projected value exceeds the second span's 1e-10.
            // Scalar sorting alone selects t=1 in both cases, whose narrow
            // own enclosure excludes the actual global extremum.
            assert!(
                answers.as_slice()[0].start.parameters[0] < 1.0,
                "maximize={maximize}"
            );
            // Preserve All's canonical retention/range merging rather than
            // inventing a new set when reported ordering is uncertifiable.
            assert!(curve
                .global_point_solve::<Handle>(
                    &Kernel::default(),
                    &interval,
                    operation,
                    Vec3::ZERO,
                    0.0,
                    SolutionRequested::All,
                    &options
                )
                .is_err());
            if maximize {
                // Maximize seeds its incumbent with target², so negative
                // targets have the same magnitude as positive targets.
                let target_options = GlobalPointSolveOptions {
                    target_distance: Some(-5e-11),
                    ..options
                };
                assert!(curve
                    .global_point_solve::<Handle>(
                        &Kernel::default(),
                        &interval,
                        operation,
                        Vec3::ZERO,
                        0.0,
                        SolutionRequested::Single,
                        &target_options
                    )
                    .is_err());
                let excluded_options = GlobalPointSolveOptions {
                    target_distance: Some(-1e-3),
                    ..options
                };
                assert!(curve
                    .global_point_solve::<Handle>(
                        &Kernel::default(),
                        &interval,
                        operation,
                        Vec3::ZERO,
                        0.0,
                        SolutionRequested::Single,
                        &excluded_options
                    )
                    .unwrap()
                    .is_empty());
            } else {
                let target_options = GlobalPointSolveOptions {
                    target_distance: Some(5e-10),
                    ..options
                };
                assert!(
                    curve
                        .global_point_solve::<Handle>(
                            &Kernel::default(),
                            &interval,
                            operation,
                            Vec3::ZERO,
                            0.0,
                            SolutionRequested::Single,
                            &target_options
                        )
                        .is_err(),
                    "a rounded report outside the target cannot certify an empty true-zero minimum"
                );
            }
        }
    }

    #[test]
    fn projected_all_refuses_an_omitted_rounding_equivalent_minimum() {
        let curve = Curve::from(BSplineCurve::from_nurb(
            NurbCurve::new(
                vec![
                    CPoint::new(0.0, 0.0, 0.0, NOW),
                    CPoint::new(0.0, 0.0, 1.0, NOW),
                    CPoint::new(1e8, 1e8, 0.0, NOW),
                    CPoint::new(-2e8, -2e8, 0.0, NOW),
                ],
                1,
                vec![0.0, 0.0, 1.0, 1.0, 2.0, 2.0],
            )
            .unwrap(),
            3,
            false,
        ));
        let interval = Extent1d::from_min_max(0.0, 2.0);
        let options = GlobalPointSolveOptions {
            opt_vectors: Some(&[Vec3::Z]),
            ..Default::default()
        };
        let k = Kernel::default();
        // The retained zero range is an accurate global minimum, but the
        // second span also has exact minimum0 at t=4/3. Its rounded report is
        // larger and canonical zero-slack filtering drops it. Checking only
        // the first retained answer would silently claim an incomplete set.
        assert!(curve
            .global_point_solve::<Handle>(
                &k,
                &interval,
                SolverOperation::ProjectedMinimize,
                Vec3::ZERO,
                0.0,
                SolutionRequested::Single,
                &options
            )
            .is_ok());
        assert!(curve
            .global_point_solve::<Handle>(
                &k,
                &interval,
                SolverOperation::ProjectedMinimize,
                Vec3::ZERO,
                0.0,
                SolutionRequested::All,
                &options
            )
            .is_err());
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
    fn cubic_near_collinear_join_retains_projected_extrema() {
        let e = 1.0 / 1024.0;
        let points = (-3_i32..=3)
            .map(|x| CPoint::new(f64::from(x), f64::from(x.abs()) * e, 0.0, NOW))
            .collect();
        let curve = Curve::from(BSplineCurve::from_nurb(
            NurbCurve::new(
                points,
                3,
                vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0],
            )
            .unwrap(),
            3,
            false,
        ));
        let (_, continuity) = curve
            .calculate_continuities(&Kernel::default(), 0.5)
            .unwrap();
        assert!((continuity[1] as i32) >= (ContinuityType::G1 as i32));
        for (op, point, expected) in [
            (
                SolverOperation::ProjectedMinimize,
                Vec3::new(0.0, -1.0, 0.0),
                1.0,
            ),
            (
                SolverOperation::ProjectedMaximize,
                Vec3::new(0.0, 4096.0, 0.0),
                4096.0,
            ),
        ] {
            for requested in [SolutionRequested::Single, SolutionRequested::All] {
                let answers = curve
                    .global_point_solve::<Handle>(
                        &Kernel::default(),
                        &curve.natural_interval(),
                        op,
                        point,
                        1e-8,
                        requested,
                        &GlobalPointSolveOptions {
                            opt_vectors: Some(&[Vec3::Z]),
                            ..Default::default()
                        },
                    )
                    .unwrap();
                assert!(
                    !answers.is_empty(),
                    "cubic {op:?} {requested:?} lost the kink"
                );
                assert!(
                    answers
                        .as_slice()
                        .iter()
                        .all(|answer| (answer.start.solution_value - expected).abs() <= 1e-8),
                    "cubic {op:?} {requested:?} lost the true value {expected}"
                );
            }
        }
    }

    #[test]
    fn gen1172_monotone_kink_is_not_an_extra_retained_minimum() {
        let points = [
            [0.07940467954627639, 3.432185442926146, -1.7857031322721397],
            [1.07064237589226, 3.2501505967236453, -1.3565118660393338],
            [2.061980072238244, 3.0681157505211445, -0.9273205998065279],
            [3.0531177685842277, 2.8860809043186437, -0.498129333573722],
        ]
        .into_iter()
        .map(|[x, y, z]| CPoint::new(x, y, z, NOW))
        .collect();
        let curve = Curve::from(BSplineCurve::from_nurb(
            NurbCurve::new(points, 1, vec![0.0, 0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0, 1.0]).unwrap(),
            3,
            false,
        ));
        let answers = curve
            .global_point_solve::<Handle>(
                &Kernel::default(),
                &curve.natural_interval(),
                SolverOperation::ProjectedMinimize,
                Vec3::new(0.9299482248193833, 3.275988232082504, -1.4174303531722534),
                0.05,
                SolutionRequested::All,
                &GlobalPointSolveOptions {
                    opt_vectors: Some(&[Vec3::new(
                        0.2768740720669518,
                        0.09239858445816851,
                        0.5062961704696287,
                    )]),
                    ..Default::default()
                },
            )
            .unwrap();
        assert_eq!(answers.len(), 1);
        assert!(answers.as_slice()[0].start.solution_value < 1e-14);
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
}

mod surface_wrappers {
    //! True UV corners remain candidates when an angular cache label says G1.
    #![allow(clippy::unwrap_used, clippy::indexing_slicing)]
    use solidcore::curve::{
        crv_on_surf::CrvOnSurf,
        curve::{p1::GlobalPointSolveOptions, Curve},
    };
    use solidcore::geom::{
        extent1d::Extent1d,
        extent2d::Extent2d,
        global_solver::{SolutionRequested, SolverOperation},
        solution_array::{ObjectKind, SolutionObject},
    };
    use solidcore::math::{Vec2, Vec3};
    use solidcore::nurbs::{
        b_spline_curve::p1::BSplineCurve,
        crv_basic::Curve as NurbCurve,
        math::{CPoint, NOW},
    };
    use solidcore::surface::{plane::Plane, surface::p2::Surface};
    use solidcore::Kernel;
    #[derive(Debug, Clone, Copy, PartialEq)]
    struct Handle;
    impl SolutionObject for Handle {
        fn kind(&self) -> ObjectKind {
            ObjectKind::Other
        }
    }
    #[test]
    fn projected_extrema_keep_small_uv_kinks_on_a_plane() {
        for maximize in [false, true] {
            let k = Kernel::default();
            // Exact dyadic data: delta > a^2 makes distance strictly increase
            // to the maximum join (or decrease to the minimum join). Neither
            // open affine span contains a stationary point.
            let a = 1.0 / 1024.0;
            let delta = 1.0 / 524288.0;
            let y = if maximize { 1.0 - delta } else { 1.0 + delta };
            let uv = Curve::from(BSplineCurve::from_nurb(
                NurbCurve::new(
                    vec![
                        CPoint::new(-a, y, 0.0, NOW),
                        CPoint::new(0.0, 1.0, 0.0, NOW),
                        CPoint::new(a, y, 0.0, NOW),
                    ],
                    1,
                    vec![0.0, 0.0, 1.0, 2.0, 2.0],
                )
                .unwrap(),
                2,
                false,
            ));
            let domain = Extent2d::from_points(Vec2::new(-2.0, -2.0), Vec2::new(2.0, 2.0));
            let surface = Surface::from(Plane::new(
                Vec3::ZERO,
                Vec3::X,
                Vec3::Y,
                Vec2::new(1.0, 1.0),
                &domain,
            ));
            let curve = Curve::from(CrvOnSurf::new(&k, uv, surface, None, 2).unwrap());
            let expected = 1.0;
            let (_, continuity) = curve.calculate_continuities(&k, 0.5).unwrap();
            assert!(
                (continuity[1] as i32)
                    >= (solidcore::geom::geom_utility::ContinuityType::G1 as i32)
            );
            let options = GlobalPointSolveOptions {
                opt_vectors: Some(&[Vec3::Z]),
                ..Default::default()
            };
            let op = if maximize {
                SolverOperation::ProjectedMaximize
            } else {
                SolverOperation::ProjectedMinimize
            };
            for requested in [SolutionRequested::Single, SolutionRequested::All] {
                let answers = curve
                    .global_point_solve::<Handle>(
                        &k,
                        &Extent1d::from_min_max(0.0, 2.0),
                        op,
                        Vec3::ZERO,
                        1e-12,
                        requested,
                        &options,
                    )
                    .unwrap();
                assert_eq!(answers.len(), 1);
                assert_eq!(answers.as_slice()[0].start.parameters[0], 1.0);
                assert_eq!(answers.as_slice()[0].start.solution_value, expected);
            }
        }
    }
}

mod linear_wrappers {
    //! Defining joins survive affine projections and reversed zero offsets.
    #![allow(clippy::unwrap_used, clippy::indexing_slicing)]
    use solidcore::{
        curve::{
            curve::{p1::GlobalPointSolveOptions, Curve},
            offset_curve::OffsetCurve,
            projected_curve::{OwnsCurve, ProjectedCurve},
        },
        geom::{
            extent1d::Extent1d,
            global_solver::{SolutionRequested, SolverOperation},
            solution_array::{ObjectKind, SolutionObject},
        },
        math::Vec3,
        nurbs::{
            b_spline_curve::{p1::BSplineCurve, ProjectionType},
            crv_basic::Curve as Net,
            math::{CPoint, NOW},
        },
        Kernel,
    };
    #[derive(Clone, Copy, PartialEq)]
    struct Handle;
    impl SolutionObject for Handle {
        fn kind(&self) -> ObjectKind {
            ObjectKind::Other
        }
    }

    fn base(maximize: bool, straight: bool) -> Curve {
        let a = 1.0 / 1024.0;
        let delta = if straight { 0.0 } else { 1.0 / 524288.0 };
        let y = 1.0 + if maximize { -delta } else { delta };
        Curve::from(BSplineCurve::from_nurb(
            Net::new(
                vec![
                    CPoint::new(-a, y, 0.0, NOW),
                    CPoint::new(0.0, 1.0, 0.0, NOW),
                    CPoint::new(a, y, 0.0, NOW),
                ],
                1,
                vec![0.0, 0.0, 0.5, 2.0, 2.0],
            )
            .unwrap(),
            3,
            false,
        ))
    }
    fn wrap(base: Curve, kind: usize) -> Curve {
        match kind {
            0 => Curve::from(
                ProjectedCurve::new(
                    base,
                    Some(Vec3::ZERO),
                    Some(Vec3::Z),
                    None,
                    ProjectionType::Parallel,
                    OwnsCurve::Yes,
                )
                .unwrap(),
            ),
            1 => Curve::from(
                ProjectedCurve::new(
                    base,
                    Some(Vec3::ZERO),
                    Some(Vec3::Y),
                    Some(Vec3::X),
                    ProjectionType::Rotation,
                    OwnsCurve::Yes,
                )
                .unwrap(),
            ),
            _ => {
                let mut offset = OffsetCurve::new(
                    3,
                    base,
                    Extent1d::from_min_max(0.0, 2.0),
                    Vec3::Z,
                    0.0,
                    true,
                )
                .unwrap();
                if kind == 3 {
                    let _ = offset.reverse_parameterization(&Extent1d::from_min_max(0.0, 2.0));
                }
                Curve::from(offset)
            }
        }
    }
    #[test]
    fn projected_extrema_follow_linear_wrapper_defining_joins() {
        let k = Kernel::default();
        for maximize in [false, true] {
            for kind in 0..4 {
                let curve = wrap(base(maximize, false), kind);
                let expected_t = if kind == 3 { 1.5 } else { 0.5 };
                for requested in [SolutionRequested::Single, SolutionRequested::All] {
                    let answers = curve
                        .global_point_solve::<Handle>(
                            &k,
                            &Extent1d::from_min_max(0.0, 2.0),
                            if maximize {
                                SolverOperation::ProjectedMaximize
                            } else {
                                SolverOperation::ProjectedMinimize
                            },
                            Vec3::ZERO,
                            1e-12,
                            requested,
                            &GlobalPointSolveOptions {
                                opt_vectors: Some(&[Vec3::Z]),
                                ..Default::default()
                            },
                        )
                        .unwrap();
                    assert_eq!(
                        answers.len(),
                        1,
                        "maximize={maximize} kind={kind} requested={requested:?}"
                    );
                    assert_eq!(answers.as_slice()[0].start.parameters[0], expected_t);
                    assert_eq!(answers.as_slice()[0].start.solution_value, 1.0);
                }
            }
        }
    }
    #[test]
    fn linear_wrappers_do_not_emit_artificial_all_boundaries() {
        for kind in 0..4 {
            let curve = wrap(base(true, true), kind);
            let answers = curve
                .global_point_solve::<Handle>(
                    &Kernel::default(),
                    &Extent1d::from_min_max(0.0, 2.0),
                    SolverOperation::ProjectedMaximize,
                    Vec3::ZERO,
                    1e-12,
                    SolutionRequested::All,
                    &GlobalPointSolveOptions {
                        opt_vectors: Some(&[Vec3::Z]),
                        ..Default::default()
                    },
                )
                .unwrap();
            assert_eq!(answers.len(), 2, "kind={kind}");
            assert!(answers
                .as_slice()
                .iter()
                .all(|s| s.start.parameters[0] == 0.0 || s.start.parameters[0] == 2.0));
        }
    }
    #[test]
    fn reversed_offset_maps_knots_multiplicities_and_continuity_order() {
        let k = Kernel::default();
        let base = Curve::from(BSplineCurve::from_nurb(
            Net::new(
                (0..6)
                    .map(|i| CPoint::new(f64::from(i), 0.0, 0.0, NOW))
                    .collect(),
                2,
                vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 2.0, 2.0, 2.0],
            )
            .unwrap(),
            3,
            false,
        ));
        let (_, before) = base.calculate_continuities(&k, 0.5).unwrap();
        let mut offset = OffsetCurve::new(
            3,
            base,
            Extent1d::from_min_max(0.0, 2.0),
            Vec3::Z,
            0.0,
            true,
        )
        .unwrap();
        let _ = offset.reverse_parameterization(&Extent1d::from_min_max(0.0, 2.0));
        assert_eq!(
            offset.knots(&k, None).unwrap(),
            (vec![0.0, 1.0, 1.5, 2.0], vec![3, 1, 2, 3])
        );
        let (_, after) = offset.calculate_continuities(&k, 0.5).unwrap();
        assert_eq!(after, before.into_iter().rev().collect::<Vec<_>>());
        let (knots, _) = offset
            .knots(&k, Some(&Extent1d::from_min_max(1.25, 1.75)))
            .unwrap();
        assert!(knots.contains(&1.5));
        assert!(knots.iter().all(|t| *t >= 1.25 && *t <= 1.75));
    }
    #[test]
    fn reversed_offset_refuses_an_unrepresentable_defining_join() {
        let base = Curve::from(BSplineCurve::from_nurb(
            Net::new(
                vec![
                    CPoint::new(0.0, 0.0, 0.0, NOW),
                    CPoint::new(1.0, 1.0, 0.0, NOW),
                    CPoint::new(2.0, 0.0, 0.0, NOW),
                ],
                1,
                vec![0.0, 0.0, 0.1, 1.0, 1.0],
            )
            .unwrap(),
            3,
            false,
        ));
        let mut offset = OffsetCurve::new(
            3,
            base,
            Extent1d::from_min_max(0.0, 1.0),
            Vec3::Z,
            0.0,
            true,
        )
        .unwrap();
        let _ = offset.reverse_parameterization(&Extent1d::from_min_max(0.0, 1.0));
        // The mirrored binary64 parameter does not evaluate at the source knot.
        assert_ne!(1.0 - (1.0 - 0.1_f64), 0.1);
        assert!(offset.knots(&Kernel::default(), None).is_err());
    }
}

mod stored_spline_labels {
    #![allow(clippy::unwrap_used, clippy::indexing_slicing)]
    use solidcore::{
        curve::{
            curve::{p1::GlobalPointSolveOptions, Curve},
            ellipse::{Ellipse, PolarConversion},
            hyperbola::Hyperbola,
            parabola::Parabola,
        },
        geom::{
            extent1d::Extent1d,
            global_solver::{SolutionRequested, SolverOperation},
            solution_array::{ObjectKind, SolutionObject},
        },
        math::Vec3,
        nurbs::{
            b_spline_curve::p1::BSplineCurve,
            crv_basic::Curve as Net,
            math::{CPoint, NOW},
        },
        topology::database_io::{DatabaseReader, DatabaseWriter, FileType, InputStream},
        Kernel,
    };
    #[derive(Clone, Copy, PartialEq)]
    struct H;
    impl SolutionObject for H {
        fn kind(&self) -> ObjectKind {
            ObjectKind::Other
        }
    }
    #[test]
    fn stored_spline_labels_keep_the_actual_defining_corner() {
        let k = Kernel::default();
        let a = 1.0 / 1024.0;
        let y = 1.0 - 1.0 / 524288.0;
        let base = BSplineCurve::from_nurb(
            Net::new(
                vec![
                    CPoint::new(-a, y, 0.0, NOW),
                    CPoint::new(0.0, 1.0, 0.0, NOW),
                    CPoint::new(a, y, 0.0, NOW),
                ],
                1,
                vec![0.0, 0.0, 0.5, 2.0, 2.0],
            )
            .unwrap(),
            3,
            false,
        );
        for kind in 0..3 {
            let mut w = DatabaseWriter::open_for_write(FileType::Binary);
            w.write_doubles(&[0., 0., 0., 1., 0., 0., 0., 1., 0., 0., 2., 1.]);
            if kind != 2 {
                w.write_double(1.);
            }
            if kind == 0 {
                w.write_boolean(0);
                PolarConversion::new()
                    .write_to_db(&mut w, 33, None, 0)
                    .unwrap();
            }
            base.write_to_db(&mut w, 33).unwrap();
            let mut r =
                DatabaseReader::with_stream(FileType::Binary, InputStream::new(w.into_bytes()));
            let curve = match kind {
                0 => Curve::from(Ellipse::read_from_db(&k, &mut r, 3, 33, false).unwrap()),
                1 => Curve::from(Hyperbola::read_from_db(&mut r, 3, 33).unwrap()),
                _ => Curve::from(Parabola::read_from_db(&mut r, 3, 33).unwrap()),
            };
            assert_eq!(
                curve.evaluate_point(&k, 0.5).unwrap(),
                Vec3::new(0., 1., 0.)
            );
            let ans = curve.global_point_solve::<H>(
                &k,
                &Extent1d::from_min_max(0., 2.),
                SolverOperation::ProjectedMaximize,
                Vec3::ZERO,
                1e-12,
                SolutionRequested::All,
                &GlobalPointSolveOptions {
                    opt_vectors: Some(&[Vec3::Z]),
                    ..Default::default()
                },
            );
            let answers = ans.unwrap();
            assert_eq!(answers.len(), 1, "kind={kind}");
            assert_eq!(answers.as_slice()[0].start.parameters[0], 0.5);
            assert_eq!(answers.as_slice()[0].start.solution_value, 1.0);
        }
    }
}

mod centered_affine_minima {
    //! Ordinary oblique line projections have attainable minima despite carrier-bound dependency.
    #![allow(clippy::unwrap_used, clippy::indexing_slicing)]
    use solidcore::{
        curve::curve::{p1::GlobalPointSolveOptions, Curve},
        geom::{
            global_solver::{SolutionRequested, SolverOperation},
            solution_array::{ObjectKind, SolutionObject},
        },
        math::Vec3,
        primitive::api_curves,
        Kernel,
    };
    #[derive(Clone, Copy, PartialEq)]
    struct Object;
    impl SolutionObject for Object {
        fn kind(&self) -> ObjectKind {
            ObjectKind::Curve
        }
    }
    fn check(start: Vec3, end: Vec3, normal: Vec3, point: Vec3, expected: f64) {
        let k = Kernel::default();
        let curve = Curve::from(api_curves::create_line_segment(start, end).unwrap());
        for tolerance in [0.0, 1e-6] {
            for requested in [SolutionRequested::Single, SolutionRequested::All] {
                let answers = curve
                    .global_point_solve::<Object>(
                        &k,
                        &curve.natural_interval(),
                        SolverOperation::ProjectedMinimize,
                        point,
                        tolerance,
                        requested,
                        &GlobalPointSolveOptions {
                            opt_vectors: Some(&[normal]),
                            ..Default::default()
                        },
                    )
                    .unwrap();
                assert_eq!(answers.len(), 1);
                let answer = &answers.as_slice()[0];
                let t = answer.start.parameters[0];
                assert!(t > 0.0 && t < 1.0);
                let native = curve.evaluate_point(&k, t).unwrap();
                assert_eq!(
                    answer.start.solution_value,
                    normal.cross(native - point).cross(normal).length()
                );
                // The reference is independently rounded from exact rational defining
                // coefficients. This comparison allows arithmetic rounding only;
                // neither the construction nor this check uses retention as accuracy.
                assert!(
                    (answer.start.solution_value - expected).abs() <= 4.0 * f64::EPSILON * expected
                );
            }
        }
    }
    #[test]
    fn ordinary_oblique_projection_0() {
        check(
            Vec3::new(-3.894398052128492, 3.8587511241420906, 17.729563382732575),
            Vec3::new(4.606891033441505, 3.73721355265872, 17.729563382732575),
            Vec3::new(
                0.01743086555439506,
                -0.09109708179454899,
                -0.1493022056749114,
            )
            .unitize()
            .unwrap(),
            Vec3::new(-1.7690757807359927, -0.8831997460150274, 0.9101188222002574),
            4.71764651872117,
        );
    }
    #[test]
    fn ordinary_oblique_projection_1() {
        check(
            Vec3::new(-5.483666909400865, 3.7359418500677646, 6.857357179858783),
            Vec3::new(16.34892834629792, 3.3584721852969905, 6.857357179858783),
            Vec3::new(0.0889704855343518, -0.008523026748961729, -1.0)
                .unitize()
                .unwrap(),
            Vec3::new(
                -0.02551809547616868,
                -1.7733953083705631,
                -0.9360354725942557,
            ),
            5.35960281390249,
        );
    }
    #[test]
    fn ordinary_oblique_projection_2() {
        check(
            Vec3::new(-12.423555212966619, 9.669598386804747, 5.0),
            Vec3::new(37.27066563889986, 10.032245827791767, 5.0),
            Vec3::new(
                0.017044433970330952,
                -0.0410992424988931,
                -2.3887968887759023,
            ),
            Vec3::new(0.0, 24.847110425933238, -0.0),
            86.59806893767214,
        );
    }
    #[test]
    fn ordinary_oblique_projection_3() {
        check(
            Vec3::new(-12.423555212966619, 9.669598386804747, 5.0),
            Vec3::new(37.27066563889986, 10.032245827791767, 5.0),
            Vec3::new(
                0.017044433970330952,
                -0.0410992424988931,
                -2.3887968887759023,
            ),
            Vec3::new(24.847110425933238, 24.847110425933238, 0.0),
            85.56319489742404,
        );
    }
}
