//! Projected face measurements must use the surface at the returned UV point.

use crate::geom::extent2d::Extent2d;
use crate::geom::global_solver::{SolutionRequested, SolverOperation};
use crate::geom::tol::USE_DEFAULT;
use crate::geom::vector2d::Vec2;
use crate::intersect::topology_solver::p1::{
    LocalSolver, LocalSolverKind, TopologyObjects, TopologySolver,
};
use crate::kernel::Kernel;
use crate::math::Vec3;
use crate::surface::sphere::{Sphere, SphereOptions};
use crate::surface::surface::Surface;
use crate::surface::torus::{Torus, TorusOptions};
use crate::topology::brep::p2::ArenaNotifyState;
use crate::topology::brep::p3::new_brep;
use crate::topology::brep::p8::{create_face_from_surface, FaceBuildEnv};
use crate::topology::brep_cache::{arena_identity, BrepObjects, CacheObject};
use crate::topology::topology_traverser::Marks;
use crate::topology::BrepId;

/// A prepopulated cache must evaluate its UV curve on S, not copy the 3D fit.
#[test]
fn cached_silhouette_measurements_use_the_source_surface() {
    use crate::curve::curve::Curve;
    use crate::nurbs::b_spline_curve::BSplineCurve;
    use crate::surface::plane::Plane;
    use crate::surface::surface::p4::NoObject;
    use crate::surface::surface_cache::p1::{
        AltCacheKey, AltSrfCacheType, SurfaceCache, SurfaceCacheOptions,
    };

    let k = Kernel::default();
    let domain = Extent2d::from_bounds(0.0, 0.0, 1.0, 1.0);
    let surface = Surface::from(Plane::new(
        Vec3::ZERO,
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec2::new(1.0, 1.0),
        &domain,
    ));
    let uv = Curve::from(
        BSplineCurve::create_line_segment(2, Vec3::new(0.0, 0.25, 0.0), Vec3::new(1.0, 0.25, 0.0))
            .unwrap(),
    );
    let fit = Curve::from(
        BSplineCurve::create_line_segment(
            3,
            Vec3::new(0.0, 0.25, 10.0),
            Vec3::new(1.0, 0.25, 10.0),
        )
        .unwrap(),
    );
    let mut cache: SurfaceCache<NoObject> =
        SurfaceCache::new(surface.clone(), domain, SurfaceCacheOptions::default());
    let key = AltCacheKey {
        vec: Some(Vec3::new(0.0, 0.0, 1.0)),
        boolean: Some(false),
        ..AltCacheKey::default()
    };
    cache
        .add_curves_to_cache(
            AltSrfCacheType::TrimmedSilhouette,
            1e-5,
            0.25,
            &key,
            Some(vec![fit]),
            Some(vec![uv]),
        )
        .unwrap();
    let (curves, uv_curves) = super::p1::cached_silhouette_measurement_curves(
        &k, &cache, &surface, domain, 1e-5, 0.25, &key,
    )
    .unwrap()
    .unwrap();
    assert_eq!((curves.len(), uv_curves.len()), (1, 1));
    let point = curves[0].evaluate_point(&k, 0.5).unwrap();
    let uv_point = uv_curves[0].evaluate_point(&k, 0.5).unwrap();
    assert_eq!(uv_point, Vec3::new(0.5, 0.25, 0.0));
    assert_eq!(
        point,
        surface
            .evaluate_point(&k, Vec2::new(uv_point.x, uv_point.y))
            .unwrap()
    );
    assert_eq!(point.z, 0.0);
    // No matching tolerance is a cache miss, not an empty successful query.
    assert!(super::p1::cached_silhouette_measurement_curves(
        &k, &cache, &surface, domain, 1e-6, 0.25, &key,
    )
    .unwrap()
    .is_none());
}

#[test]
fn projected_sphere_face_maximum_matches_surface_and_closed_form() {
    check_sphere_projection(0.0, 1.0);
}

#[test]
fn projected_tilted_sphere_face_maximum_uses_projection_metric() {
    check_sphere_projection(45.0, 0.001);
    check_sphere_projection(45.0, 1.0);
    check_sphere_projection(45.0, 1000.0);
}

fn check_sphere_projection(tilt: f64, scale: f64) {
    let k = Kernel::default();
    let mut sphere = Sphere::new(
        &k,
        Vec3::new(0.0, 0.0, -1.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Extent2d::from_bounds(0.0, -60.0, 360.0, 60.0),
        0.5,
        SphereOptions::default(),
    )
    .unwrap();
    sphere.make_nurb(&k).unwrap();
    let angle = tilt.to_radians();
    let normal = Vec3::new(angle.sin(), 0.0, angle.cos()) * scale;
    let query = Vec3::new(0.2, 0.3, 2.0);
    let expected = 0.5 * normal.length_squared()
        + normal
            .cross(query - Vec3::new(0.0, 0.0, -1.0))
            .cross(normal)
            .length();
    check_face_projection(k, &Surface::from(sphere), normal, expected);
}

#[test]
fn projected_trimmed_sphere_maximum_keeps_the_trim_endpoint() {
    let k = Kernel::default();
    let mut sphere = Sphere::new(
        &k,
        Vec3::new(0.0, 0.0, -1.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Extent2d::from_bounds(0.0, -45.0, 180.0, 45.0),
        0.5,
        SphereOptions::default(),
    )
    .unwrap();
    sphere.make_nurb(&k).unwrap();
    // The unconstrained farthest equator point is outside this half-longitude
    // face. Its farthest included point is the trim endpoint (-0.5, 0, -1).
    let expected = (0.7_f64 * 0.7 + 0.3 * 0.3).sqrt();
    check_face_projection(k, &Surface::from(sphere), Vec3::new(0.0, 0.0, 1.0), expected);
}

#[test]
fn projected_torus_maximum_matches_the_capsule_projection() {
    let k = Kernel::default();
    let mut torus = Torus::new(
        &k,
        Vec3::new(0.0, 0.0, -2.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        Extent2d::from_bounds(0.0, 0.0, 360.0, 360.0),
        1.0,
        0.25,
        TorusOptions::default(),
    )
    .unwrap();
    torus.make_nurb(&k).unwrap();
    // The projected torus is the radius-0.25 capsule around the y-axis segment
    // [-1, 1]. The farther cap is centered at (0, -1).
    let expected = (0.2_f64 * 0.2 + 1.3 * 1.3).sqrt() + 0.25;
    check_face_projection(k, &Surface::from(torus), Vec3::new(0.0, 0.0, 1.0), expected);
}

fn check_face_projection(mut k: Kernel, surface: &Surface, normal: Vec3, expected: f64) {
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
    let query = Vec3::new(0.2, 0.3, 2.0);
    let tol = 1e-5;
    let mut solver = TopologySolver::new(
        &k,
        &objects,
        SolverOperation::ProjectedMaximize,
        SolutionRequested::Single,
        tol,
        0.0,
        &[normal],
    );
    let mut local = LocalSolver::new(LocalSolverKind::PointFace, false, None);
    local.set_point(query);
    local
        .solve_it(&k, &mut solver, None, Some(CacheObject::Surface(face)))
        .unwrap();
    assert_eq!(solver.base.solutions.len(), 1);
    let solution = &solver.base.solutions.as_slice()[0];
    let uv = Vec2::new(
        solution.start.param(0).unwrap(),
        solution.start.param(1).unwrap(),
    );
    let point = surface.evaluate_point(&k, uv).unwrap();
    let measured = normal.cross(point - query).cross(normal).length();
    assert!(
        (solution.start.solution_value - expected).abs() <= tol,
        "value={}, expected={expected}",
        solution.start.solution_value
    );
    assert!(
        (solution.start.solution_value - measured).abs() <= 8.0 * f64::EPSILON * measured.max(1.0),
        "value={}, surface measurement={measured}",
        solution.start.solution_value
    );
}
