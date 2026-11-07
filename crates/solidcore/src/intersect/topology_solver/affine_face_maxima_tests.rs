//! Face-owned pair maxima must be attained on both native definitions.
#![allow(clippy::unwrap_used, clippy::indexing_slicing)]
use crate::curve::curve::Curve;
use crate::geom::extent1d::Extent1d;
use crate::geom::extent2d::Extent2d;
use crate::geom::global_solver::{SolutionRequested, SolverOperation};
use crate::geom::solution_array::{ObjectKind, SolutionObject};
use crate::geom::tol::USE_DEFAULT;
use crate::intersect::topology_solver::p1::{
    CurveRef, LocalSolver, LocalSolverKind, TopologyObjects, TopologySolver,
};
use crate::math::{Vec2, Vec3};
use crate::nurbs::b_spline_curve::BSplineCurve;
use crate::surface::{
    surface::Surface,
    torus::{Torus, TorusOptions},
};
use crate::topology::{
    brep::{
        p2::ArenaNotifyState,
        p3::new_brep,
        p8::{create_face_from_surface, FaceBuildEnv},
    },
    brep_cache::{arena_identity, BrepObjects, CacheObject},
    topology_traverser::Marks,
    BrepId,
};
use crate::Kernel;

fn check_torus(major_end: f64) {
    check_torus_options(major_end, 1.0, false, Extent1d::from_min_max(0.0, 1.0));
}

fn check_torus_options(major_end: f64, scale: f64, swapped: bool, interval: Extent1d) {
    check_torus_mode(major_end, scale, swapped, interval, 0);
}
fn check_torus_mode(major_end: f64, scale: f64, swapped: bool, interval: Extent1d, mode: u8) {
    let k = Kernel::default();
    let mut torus = Torus::new(
        &k,
        Vec3::new(0.0, 0.0, -2.0),
        Vec3::Y,
        Vec3::Z,
        Extent2d::from_bounds(0.0, 0.0, major_end, 360.0),
        1.0,
        0.25,
        TorusOptions::default(),
    )
    .unwrap();
    torus.make_nurb(&k).unwrap();
    let surface = Surface::from(torus);
    let reach = (-1.0 + 2.0 * interval.min())
        .abs()
        .max((-1.0 + 2.0 * interval.max()).abs());
    check_surface(
        k,
        &surface,
        (reach * reach + 1.0).sqrt() + 0.25,
        scale,
        swapped,
        interval,
        mode,
    );
}

fn check_surface(
    mut k: Kernel,
    surface: &Surface,
    expected: f64,
    scale: f64,
    swapped: bool,
    interval: Extent1d,
    mode: u8,
) {
    let domain = surface.natural_uv_domain_with_nurb(&k).unwrap();
    let mut b = new_brep(&mut k, USE_DEFAULT, USE_DEFAULT).unwrap();
    let mut state = ArenaNotifyState::new(&b);
    let mut env = state.env(None);
    let mut marks = Marks::default();
    let keys = |f| arena_identity(CacheObject::Face(f));
    let mut fenv = FaceBuildEnv {
        env: &mut env,
        marks: &mut marks,
        this_face: &keys,
    };
    let face = create_face_from_surface(&mut k, &mut b, Some(surface.clone()), &domain, &mut fenv)
        .unwrap();
    let extra_face = if mode == 5 {
        let mut other = Torus::new(
            &k,
            Vec3::new(0.0, 0.0, -4.0),
            Vec3::Y,
            Vec3::Z,
            Extent2d::from_bounds(0.0, 0.0, 360.0, 360.0),
            0.5,
            0.25,
            TorusOptions::default(),
        )
        .unwrap();
        other.make_nurb(&k).unwrap();
        Some(
            create_face_from_surface(
                &mut k,
                &mut b,
                Some(Surface::from(other)),
                &domain,
                &mut fenv,
            )
            .unwrap(),
        )
    } else {
        None
    };
    let mut preparation = crate::query::brep::RayPreparation::default();
    {
        let objects = BrepObjects {
            preparation: Some(&preparation),
            brep: &b,
            brep_id: BrepId::from_index(0),
            cache: None,
            keys: &arena_identity,
            owner_of: &|_| None,
        };
        let mut context = objects.surface_solve_context().unwrap();
        assert!(matches!(
            context.build_surface_cache(
                &k,
                b.face(face).unwrap().surface.as_ref().unwrap(),
                Some(CacheObject::Surface(face))
            ),
            Err(crate::SolidError::PreparationRequired { .. })
        ));
    }
    let request = preparation.take_request().unwrap();
    preparation.prepare(&k, &mut b, request).unwrap();
    let object_keys = |object| match object {
        CacheObject::Edge(id) | CacheObject::Curve(id) => {
            b.owner_identity(crate::error::EntityRef::Edge(id))
        }
        _ => arena_identity(object),
    };
    let objects = BrepObjects {
        preparation: Some(&preparation),
        brep: &b,
        brep_id: BrepId::from_index(0),
        cache: None,
        keys: &object_keys,
        owner_of: &|_| None,
    };
    let curve = Curve::from(
        BSplineCurve::create_line_segment(3, Vec3::new(-1.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0))
            .unwrap(),
    );
    // The query is external to the face B-rep; its handle names that supplied
    // definition and is never a lookup of one of the face's seam curves.
    let query_handle = CacheObject::Other(arena_identity(CacheObject::Brep));
    let normal = Vec3::Z * scale;
    let query_box = crate::geom::extent3d::Extent3d::from_points(
        Vec3::new(-1.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
    );
    let face_box = crate::geom::extent3d::Extent3d::from_points(
        Vec3::new(-4.0, -4.0, -5.0),
        Vec3::new(4.0, 4.0, 4.0),
    );
    let query_tree = object_tree(query_handle, query_box);
    let face_tree = object_tree(CacheObject::Surface(face), face_box);
    let vertex_trees: Vec<_> = b
        .vertices()
        .iter()
        .map(|(id, _)| object_tree(CacheObject::Vertex(id), face_box))
        .collect();
    let edge_trees: Vec<_> = b
        .edges()
        .iter()
        .filter(|(_, edge)| edge.curve.is_some())
        .map(|(id, _)| object_tree(CacheObject::Curve(id), face_box))
        .collect();
    if mode == 10 {
        let token = crate::kernel::CancelToken::new();
        token.cancel();
        k.set_cancel_token(token);
    }
    let mut solver = TopologySolver::new(
        &k,
        &objects,
        SolverOperation::ProjectedMaximize,
        SolutionRequested::Single,
        if mode == 8 { 0.0 } else { 1e-5 },
        match mode {
            7 => 10_000.0,
            8 => (expected * scale * scale) * (expected * scale * scale),
            9 => 0.25,
            _ => 0.0,
        },
        &[normal],
    );
    if mode == 14 {
        use crate::geom::global_solver::GlobalSolve;
        solver.base.operation = SolverOperation::Minimize;
        solver.base.num_trees = 1;
        solver.base.trees[0] = Some(&query_tree);
        let single = [query_tree.top_node(), None, None];
        assert!(GlobalSolve::branch_may_contain_answers(&mut solver, &single).unwrap());
        // A disjoint branch is culled before a local solver is needed. The native
        // pair lifecycle must not demand one for this unrelated operation.
        let distant = object_tree(
            CacheObject::Surface(face),
            crate::geom::extent3d::Extent3d::from_points(
                Vec3::new(100.0, 100.0, 100.0),
                Vec3::new(101.0, 101.0, 101.0),
            ),
        );
        solver.base.num_trees = 2;
        solver.base.trees[1] = Some(&distant);
        let branch = [query_tree.top_node(), distant.top_node(), None];
        GlobalSolve::solve_branch(&mut solver, &branch).unwrap();
        assert!(solver.base.solutions.is_empty());
        return;
    }
    let mut local = LocalSolver::new(LocalSolverKind::CurveFace, swapped, None);
    local.set_curve(
        Some(CurveRef {
            curve: &curve,
            handle: query_handle,
        }),
        interval,
    );
    if mode == 4 {
        let unknown = LocalSolver::new(LocalSolverKind::FaceFace, false, None);
        solver
            .projected_pair_domains
            .record_local(
                &solver.base,
                &objects,
                &unknown,
                Some(CacheObject::Surface(face)),
                Some(CacheObject::Surface(face)),
            )
            .unwrap();
        assert!(matches!(
            local.solve_it(&k, &mut solver, None, Some(CacheObject::Surface(face))),
            Err(crate::SolidError::Unsupported { .. })
        ));
        return;
    }
    if mode == 6 {
        solver
            .projected_pair_domains
            .record_local(
                &solver.base,
                &objects,
                &local,
                None,
                Some(CacheObject::Surface(face)),
            )
            .unwrap();
        let mut outside = crate::geom::solution_array::Solution::new();
        outside.num_objects = 2;
        outside.objects[0] = Some(query_handle);
        outside.objects[1] = Some(CacheObject::Face(face));
        outside.num_variables = 3;
        outside.solution_type = crate::geom::solution_array::SolutionType::SingleValue;
        outside.start.parameters[0] = 0.0;
        outside.start.parameters[1] = 0.5;
        outside.start.parameters[2] = 0.125;
        outside.start.solution_value = normal
            .cross(
                curve.evaluate_point(&k, 0.0).unwrap()
                    - surface.evaluate_point(&k, Vec2::new(0.5, 0.125)).unwrap(),
            )
            .cross(normal)
            .length();
        solver.base.best_answer_so_far_sq =
            outside.start.solution_value * outside.start.solution_value;
        solver.base.solutions.add(&outside);
    }
    if matches!(mode, 11..=13 | 15..=19) {
        solver.base.trees[0] = Some(&query_tree);
        solver.base.trees[1] = Some(&face_tree);
        solver.base.num_trees = 2;
        solver.local_solver = Some(local.clone());
        let branch = [query_tree.top_node(), face_tree.top_node(), None];
        if mode == 12 {
            assert!(matches!(
                crate::geom::global_solver::solve_branch(&mut solver, &branch),
                Err(crate::SolidError::Unsupported { .. })
            ));
            return;
        }
        if mode == 13 {
            assert!(matches!(
                solver.local_solve_branch(&k, &branch),
                Err(crate::SolidError::Unsupported { .. })
            ));
            return;
        }
        crate::geom::global_solver::GlobalSolve::solve_branch(&mut solver, &branch).unwrap();
        if matches!(mode, 15..=19) {
            use crate::geom::global_solver::{GlobalSolve, SolveRequest};
            use crate::geom::solution_array::SolutionArray;
            assert!(solver.projected_pair_domains.active);
            solver.local_solver.as_mut().unwrap().interval = Extent1d::from_min_max(0.2, 0.8);
            if mode == 19 {
                // Intentional branch accumulation must preserve the first full
                // interval and its larger attained maximum.
                GlobalSolve::solve_branch(&mut solver, &branch).unwrap();
                let value = solver.base.solutions.get(0).unwrap().start.solution_value;
                assert!((value - expected).abs() <= 8.0 * f64::EPSILON * expected);
                return;
            }
            let vectors = [normal * 2.0];
            let request = SolveRequest {
                operation: if mode == 16 {
                    SolverOperation::Minimize
                } else {
                    SolverOperation::ProjectedMaximize
                },
                solution_requested: SolutionRequested::Single,
                tolerance_3d: 1e-5,
                best_answer_so_far_sq: 0.0,
                opt_vectors: &vectors,
            };
            let mut answers = SolutionArray::new();
            if mode == 17 {
                answers.add(solver.base.solutions.get(0).unwrap());
            }
            if mode == 16 {
                solver.base.trees[1] = None;
            }
            if mode == 18 {
                let token = crate::kernel::CancelToken::new();
                token.cancel();
                solver.base.cancel = Some(token);
            }
            let result = GlobalSolve::solve_it(&mut solver, &request, &mut answers);
            if mode == 16 {
                assert!(result.is_err());
                assert!(!solver.projected_pair_domains.active);
                return;
            }
            if mode == 17 {
                assert!(matches!(result, Err(crate::SolidError::Unsupported { .. })));
                return;
            }
            if mode == 18 {
                assert!(matches!(result, Err(crate::SolidError::Cancelled { .. })));
                // The newly captured trim can be active, but no old full-domain
                // candidate survives the fresh request reset.
                solver.base.cancel = None;
                GlobalSolve::solve_it(&mut solver, &request, &mut answers).unwrap();
            } else {
                result.unwrap();
            }
            let answer = answers.get(0).unwrap();
            assert!((0.2..=0.8).contains(&answer.start.param(0).unwrap()));
            let exact = ((0.6_f64 * 0.6 + 1.0).sqrt() + 0.25) * 4.0;
            assert!((answer.start.solution_value - exact).abs() <= 8.0 * f64::EPSILON * exact);
            return;
        }
    } else if matches!(mode, 0 | 5 | 6 | 7 | 8 | 9 | 10) {
        let named_face = Some(CacheObject::Surface(face));
        let result = local.solve_it(
            &k,
            &mut solver,
            if swapped { named_face } else { None },
            if swapped { None } else { named_face },
        );
        if mode == 10 {
            assert!(matches!(result, Err(crate::SolidError::Cancelled { .. })));
            return;
        }
        if mode == 8 {
            assert!(matches!(result, Err(crate::SolidError::Unsupported { .. })));
            return;
        }
        result.unwrap();
        if mode == 7 {
            assert!(solver.base.solutions.is_empty());
            return;
        }
        if let Some(extra_face) = extra_face {
            local
                .solve_it(
                    &k,
                    &mut solver,
                    None,
                    Some(CacheObject::Surface(extra_face)),
                )
                .unwrap();
        }
    } else {
        if mode >= 2 {
            let mut vertices = LocalSolver::new(LocalSolverKind::CurveVertex, false, None);
            vertices.set_curve(
                Some(CurveRef {
                    curve: &curve,
                    handle: query_handle,
                }),
                interval,
            );
            for tree in &vertex_trees {
                solver
                    .solve_trees(&query_tree, tree, &mut vertices)
                    .unwrap();
            }
            let mut edges = LocalSolver::new(LocalSolverKind::CurveEdge, false, None);
            edges.set_curve(
                Some(CurveRef {
                    curve: &curve,
                    handle: query_handle,
                }),
                interval,
            );
            for tree in &edge_trees {
                solver.solve_trees(&query_tree, tree, &mut edges).unwrap();
            }
            if mode == 3 {
                // Model an inherited fitted scalar that would prune the winning
                // native face. Original vertex/edge domains remain in the ledger.
                solver
                    .base
                    .solutions
                    .get_mut(0)
                    .unwrap()
                    .start
                    .solution_value = 100.0;
                solver.base.best_answer_so_far_sq = 10_000.0;
                assert!(!solver.may_contain_answer(&query_box, &face_box).unwrap());
            }
        }
        solver
            .solve_trees(&query_tree, &face_tree, &mut local)
            .unwrap();
    }
    assert_eq!(solver.base.solutions.len(), 1);
    let answer = &solver.base.solutions.as_slice()[0];
    let t = answer.start.parameters[if swapped { 2 } else { 0 }];
    let uv = if swapped {
        Vec2::new(answer.start.parameters[0], answer.start.parameters[1])
    } else {
        Vec2::new(answer.start.parameters[1], answer.start.parameters[2])
    };
    assert!(interval.min() <= t && t <= interval.max() && domain.contains_point(uv, 0.0));
    let a = curve.evaluate_point(&k, t).unwrap();
    let returned_face = answer
        .objects
        .iter()
        .flatten()
        .copied()
        .find(|owner| owner.kind() == ObjectKind::Face)
        .unwrap();
    let p = objects
        .surface(objects.face_surface(returned_face).unwrap())
        .unwrap()
        .evaluate_point(&k, uv)
        .unwrap();
    let measured = normal.cross(a - p).cross(normal).length();
    let exact = expected * scale * scale;
    let arithmetic_check = 8.0 * f64::EPSILON * exact;
    assert!(
        (answer.start.solution_value - measured).abs() <= arithmetic_check,
        "reported={}, native={measured}, t={t}, uv={uv:?}",
        answer.start.solution_value
    );
    assert!(
        (measured - exact).abs() <= arithmetic_check,
        "native={measured}, bound={exact}"
    );
}

#[test]
fn full_torus_curve_face_maximum_is_native_and_global() {
    check_torus(360.0);
}
#[test]
fn half_torus_curve_face_maximum_keeps_pair_fix() {
    check_torus(180.0);
}

#[test]
fn full_torus_affine_pair_preserves_swapped_parameter_order() {
    check_torus_options(360.0, 1.0, true, Extent1d::from_min_max(0.0, 1.0));
}
#[test]
fn full_torus_affine_pair_respects_curve_subdomain() {
    check_torus_options(360.0, 1.0, false, Extent1d::from_min_max(0.2, 0.8));
}
#[test]
fn full_torus_affine_pair_preserves_normal_squared_scale() {
    for scale in [0.001, 1000.0] {
        check_torus_options(360.0, scale, false, Extent1d::from_min_max(0.0, 1.0));
    }
}
#[test]
fn full_torus_affine_pair_zero_normal_has_zero_metric() {
    check_torus_options(360.0, 0.0, false, Extent1d::from_min_max(0.0, 1.0));
}

fn object_tree(
    object: CacheObject,
    bbox: crate::geom::extent3d::Extent3d,
) -> crate::geom::tree::Tree<CacheObject> {
    let mut tree = crate::geom::tree::Tree::from_point(Vec3::ZERO, 0.0, true);
    let node = tree.top_node().unwrap();
    tree.node_mut(node).unwrap().bbox = bbox;
    let list = tree.node(node).unwrap().list_head().unwrap();
    let item = tree.list_mut(list).unwrap();
    item.object = Some(object);
    item.bbox = bbox;
    tree
}
#[test]
fn native_affine_face_tree_entry_matches_direct_entry() {
    check_torus_mode(360.0, 1.0, false, Extent1d::from_min_max(0.0, 1.0), 1);
}
#[test]
fn native_affine_face_accounts_for_vertex_edge_and_face_phases() {
    check_torus_mode(360.0, 1.0, false, Extent1d::from_min_max(0.0, 1.0), 2);
}
#[test]
fn native_affine_face_cannot_be_pruned_by_an_unattained_legacy_scalar() {
    check_torus_mode(360.0, 1.0, false, Extent1d::from_min_max(0.0, 1.0), 3);
}

#[test]
fn native_affine_face_refuses_unknown_mixed_domain_history() {
    check_torus_mode(360.0, 1.0, false, Extent1d::from_min_max(0.0, 1.0), 4);
}
#[test]
fn native_affine_face_preserves_coverage_across_different_face_owners() {
    check_torus_mode(360.0, 1.0, false, Extent1d::from_min_max(0.0, 1.0), 5);
}
#[test]
fn native_affine_cone_control_keeps_sqrt_three_maximum() {
    use crate::geom::axis2_placement::Axis2Placement;
    use crate::surface::cone::Cone;
    let k = Kernel::default();
    let tangent = (30.0 * std::f64::consts::PI / 180.0).tan();
    let placement = Axis2Placement::new(Vec3::new(-0.5, 0.0, 0.0), Vec3::Y, Vec3::Z);
    let cone = Cone::create_canonical(
        &k,
        &placement,
        0.5 * tangent,
        1.5 * tangent,
        1.0,
        false,
        false,
    )
    .unwrap();
    check_surface(
        k,
        &Surface::from(cone),
        3.0_f64.sqrt(),
        1.0,
        false,
        Extent1d::from_min_max(0.0, 1.0),
        0,
    );
}
#[test]
fn native_affine_sphere_band_and_quarter_controls_keep_one_point_five() {
    use crate::surface::sphere::{Sphere, SphereOptions};
    for (longitude, latitude) in [(360.0, 60.0), (180.0, 45.0)] {
        let k = Kernel::default();
        let mut sphere = Sphere::new(
            &k,
            Vec3::new(0.0, 0.0, -1.0),
            Vec3::X,
            Vec3::Y,
            Extent2d::from_bounds(0.0, -latitude, longitude, latitude),
            0.5,
            SphereOptions::default(),
        )
        .unwrap();
        sphere.make_nurb(&k).unwrap();
        check_surface(
            k,
            &Surface::from(sphere),
            1.5,
            1.0,
            false,
            Extent1d::from_min_max(0.0, 1.0),
            0,
        );
    }
}

#[test]
fn native_affine_face_rejects_retained_native_point_outside_requested_curve_trim() {
    check_torus_mode(360.0, 1.0, false, Extent1d::from_min_max(0.2, 0.8), 6);
}
#[test]
fn native_affine_face_supports_singleton_curve_domain() {
    check_torus_mode(360.0, 1.0, false, Extent1d::from_min_max(0.25, 0.25), 0);
}

#[test]
fn native_affine_face_external_target_proves_empty_or_retains_native_answer() {
    for mode in [7, 9] {
        check_torus_mode(360.0, 1.0, false, Extent1d::from_min_max(0.0, 1.0), mode);
    }
}
#[test]
fn native_affine_face_refuses_uncertain_external_cutoff_without_using_retention_as_accuracy() {
    check_torus_mode(360.0, 1.0, false, Extent1d::from_min_max(0.0, 1.0), 8);
}

#[test]
fn native_affine_face_propagates_cancellation() {
    check_torus_mode(360.0, 1.0, false, Extent1d::from_min_max(0.0, 1.0), 10);
}

#[test]
fn native_affine_face_trait_branch_entry_runs_complete_ledger_lifecycle() {
    check_torus_mode(360.0, 1.0, false, Extent1d::from_min_max(0.0, 1.0), 11);
}
#[test]
fn native_affine_face_raw_free_entry_refuses_missing_finalization_hook() {
    check_torus_mode(360.0, 1.0, false, Extent1d::from_min_max(0.0, 1.0), 12);
}
#[test]
fn native_affine_face_partial_branch_entry_does_not_claim_complete_coverage() {
    check_torus_mode(360.0, 1.0, false, Extent1d::from_min_max(0.0, 1.0), 13);
}

#[test]
fn affine_endpoint_budget_is_checked_before_collecting_large_definition() {
    for spans in [128usize, 129usize] {
        let points: Vec<_> = (0..=spans)
            .map(|i| crate::nurbs::math::CPoint::new(i as f64, 0.0, 0.0, crate::nurbs::math::NOW))
            .collect();
        let mut knots = vec![0.0, 0.0];
        knots.extend((1..spans).map(|i| i as f64));
        knots.extend([spans as f64, spans as f64]);
        let curve =
            Curve::from(BSplineCurve::from_data(3, spans, spans + 2, 1, &points, &knots).unwrap());
        let result =
            crate::curve::curve::p1::affine_curve_endpoints(&curve, curve.natural_interval());
        if spans == 128 {
            assert_eq!(result.unwrap().unwrap().len(), 129);
        } else {
            assert!(matches!(result, Err(crate::SolidError::Unsupported { .. })));
        }
    }
}

#[test]
fn native_pair_ledger_preserves_unrelated_single_tree_and_culled_minimum_entries() {
    check_torus_mode(360.0, 1.0, false, Extent1d::from_min_max(0.0, 1.0), 14);
}

#[test]
fn native_pair_fresh_request_replaces_domains_and_projection_options() {
    check_torus_mode(360.0, 1.0, false, Extent1d::from_min_max(0.0, 1.0), 15);
}
#[test]
fn native_pair_fresh_ineligible_error_clears_old_active_history() {
    check_torus_mode(360.0, 1.0, false, Extent1d::from_min_max(0.0, 1.0), 16);
}
#[test]
fn native_pair_fresh_request_refuses_uncertified_supplied_answer_history() {
    check_torus_mode(360.0, 1.0, false, Extent1d::from_min_max(0.0, 1.0), 17);
}
#[test]
fn native_pair_cancelled_fresh_request_can_be_retried_without_stale_history() {
    check_torus_mode(360.0, 1.0, false, Extent1d::from_min_max(0.0, 1.0), 18);
}
#[test]
fn native_pair_explicit_branch_accumulation_preserves_prior_domain() {
    check_torus_mode(360.0, 1.0, false, Extent1d::from_min_max(0.0, 1.0), 19);
}
