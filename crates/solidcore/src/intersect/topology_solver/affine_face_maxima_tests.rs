//! Face-owned pair maxima must be attained on both native definitions.
#![allow(clippy::unwrap_used, clippy::indexing_slicing)]
use crate::curve::curve::Curve;
use crate::geom::extent1d::Extent1d;
use crate::geom::extent2d::Extent2d;
use crate::geom::global_solver::{SolutionRequested, SolverOperation};
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
    let mut k = Kernel::default();
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
    let objects = BrepObjects {
        preparation: Some(&preparation),
        brep: &b,
        brep_id: BrepId::from_index(0),
        cache: None,
        keys: &arena_identity,
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
    let mut solver = TopologySolver::new(
        &k,
        &objects,
        SolverOperation::ProjectedMaximize,
        SolutionRequested::Single,
        1e-5,
        0.0,
        &[normal],
    );
    let mut local = LocalSolver::new(LocalSolverKind::CurveFace, swapped, None);
    local.set_curve(
        Some(CurveRef {
            curve: &curve,
            handle: query_handle,
        }),
        interval,
    );
    let named_face = Some(CacheObject::Surface(face));
    local
        .solve_it(
            &k,
            &mut solver,
            if swapped { named_face } else { None },
            if swapped { None } else { named_face },
        )
        .unwrap();
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
    let p = surface.evaluate_point(&k, uv).unwrap();
    let measured = normal.cross(a - p).cross(normal).length();
    let reach = (-1.0 + 2.0 * interval.min())
        .abs()
        .max((-1.0 + 2.0 * interval.max()).abs());
    let exact = ((reach * reach + 1.0).sqrt() + 0.25) * scale * scale;
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
