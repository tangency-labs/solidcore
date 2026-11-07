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
    let copies = cache
        .has_cached_curves(
            AltSrfCacheType::TrimmedSilhouette,
            1e-5,
            0.25,
            &key,
            true,
            true,
        )
        .unwrap()
        .unwrap();
    let uv_curves = copies.curves_uv.unwrap();
    let curves = super::p1::lift_silhouette_measurements(&k, &uv_curves, &surface, domain).unwrap();
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
    assert!(cache
        .has_cached_curves(
            AltSrfCacheType::TrimmedSilhouette,
            1e-6,
            0.25,
            &key,
            true,
            true,
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
    check_face_projection(
        k,
        &Surface::from(sphere),
        Vec3::new(0.0, 0.0, 1.0),
        expected,
    );
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

#[test]
fn projected_plane_maximum_uses_its_boundary_without_a_silhouette() {
    let k = Kernel::default();
    let domain = Extent2d::from_bounds(-0.5, -0.5, 0.5, 0.5);
    let surface = Surface::from(crate::surface::plane::Plane::new(
        Vec3::ZERO,
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec2::new(1.0, 1.0),
        &domain,
    ));
    check_face_projection(
        k,
        &surface,
        Vec3::new(0.0, 0.0, 1.0),
        (0.7_f64 * 0.7 + 0.8 * 0.8).sqrt(),
    );
}

#[test]
fn projected_sphere_band_maximum_uses_a_trim_away_from_the_silhouette() {
    let k = Kernel::default();
    let mut sphere = Sphere::new(
        &k,
        Vec3::new(0.0, 0.0, -1.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Extent2d::from_bounds(0.0, 15.0, 360.0, 60.0),
        0.5,
        SphereOptions::default(),
    )
    .unwrap();
    sphere.make_nurb(&k).unwrap();
    let expected = 0.5 * 15.0_f64.to_radians().cos() + (0.2_f64 * 0.2 + 0.3 * 0.3).sqrt();
    check_face_projection(
        k,
        &Surface::from(sphere),
        Vec3::new(0.0, 0.0, 1.0),
        expected,
    );
}

#[test]
fn projected_brep_maximum_is_certified_against_the_final_candidate_in_either_order() {
    for (specs, query_x, expected) in [
        (vec![(0.0, 0.5), (0.1, 0.6)], 1e-6, 0.7 - 1e-6),
        (vec![(0.1, 0.6), (0.0, 0.5)], 1e-6, 0.7 - 1e-6),
        // All seam vertices project to (0.25, 0): a zero vertex answer
        // must not terminate a maximum query before its edges and face.
        (vec![(0.0, 0.5)], 0.25, 0.75),
    ] {
        let mut k = Kernel::default();
        let mut b = new_brep(&mut k, USE_DEFAULT, USE_DEFAULT).unwrap();
        let mut faces = Vec::new();
        for &(x, radius) in &specs {
            let mut sphere = Sphere::new(
                &k,
                Vec3::new(x, 0.0, -1.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                Extent2d::from_bounds(0.0, -60.0, 360.0, 60.0),
                radius,
                SphereOptions::default(),
            )
            .unwrap();
            sphere.make_nurb(&k).unwrap();
            let surface = Surface::from(sphere);
            let domain = surface.natural_uv_domain_with_nurb(&k).unwrap();
            let mut state = ArenaNotifyState::new(&b);
            let mut env = state.env(None);
            let mut marks = Marks::default();
            let keys = |f| arena_identity(CacheObject::Face(f));
            let mut fenv = FaceBuildEnv {
                env: &mut env,
                marks: &mut marks,
                this_face: &keys,
            };
            faces.push(
                create_face_from_surface(&mut k, &mut b, Some(surface), &domain, &mut fenv)
                    .unwrap(),
            );
        }
        let mut preparation = crate::query::brep::RayPreparation::default();
        for face in faces {
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
        }
        let mut cache =
            crate::topology::brep_cache::BrepCache::new(crate::geom::a_object::BREP_TYPE).unwrap();
        cache
            .build_trees(
                &k,
                crate::topology::brep_cache::CacheSource::Brep(&b),
                |o| o,
                &arena_identity,
            )
            .unwrap();
        let objects = BrepObjects {
            preparation: Some(&preparation),
            brep: &b,
            brep_id: BrepId::from_index(0),
            cache: Some(&cache),
            keys: &arena_identity,
            owner_of: &|_| None,
        };
        let tol = 1e-8;
        let answer = TopologySolver::brep_point_solve(
            &k,
            &objects,
            Some(CacheObject::Brep),
            Vec3::new(query_x, 0.0, 2.0),
            crate::intersect::topology_solver::p2::TopologySolveOptions {
                operation: SolverOperation::ProjectedMaximize,
                requested: SolutionRequested::Single,
                tolerance_3d: tol,
                best_answer_so_far_sq: 0.0,
                opt_vectors: &[Vec3::new(0.0, 0.0, 1.0)],
            },
        )
        .unwrap_or_else(|error| panic!("specs={specs:?}, query_x={query_x}: {error:?}"));
        assert_eq!(answer.len(), 1);
        assert!((answer.as_slice()[0].start.solution_value - expected).abs() <= tol);
    }
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
    assert!(super::projected_bound::is_exact_rectangle(&k, &b, face, &domain).unwrap());
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
    let net = crate::blend::fillet_section::station::Nurbs::of(
        &crate::blend::fillet_section::station::as_b_spline(&k, surface)
            .unwrap()
            .unwrap(),
    )
    .unwrap();
    let enclosed = super::projected_bound::value_at(&net, uv, query, normal).unwrap();
    assert!(enclosed.lo <= measured && measured <= enclosed.hi);
    assert!(
        super::projected_bound::maximum_is_bounded(
            &k,
            &net,
            &domain,
            query,
            normal,
            measured + tol * 0.9,
        )
        .unwrap(),
        "whole-domain bound did not close: normal={normal:?}"
    );
    if normal == Vec3::new(0.0, 0.0, 1.0) && expected > 0.86 && expected < 0.87 {
        assert!(
            !super::projected_bound::maximum_is_bounded(
                &k,
                &net,
                &domain,
                query,
                normal,
                0.8605419808705617 + tol,
            )
            .unwrap(),
            "the inaccurate original result must not be certified"
        );
    }
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
    // A plausible scalar is not enough: it must agree with the returned point.
    for delta in [-2.0 * tol, 2.0 * tol] {
        let mut incorrect = solver.base.clone();
        let mut changed = *solution;
        changed.start.solution_value += delta;
        incorrect.solutions.set_at(0, &changed).unwrap();
        assert!(super::projected_bound::certify(
            &k,
            &objects,
            CacheObject::Face(face),
            query,
            &incorrect,
        )
        .is_err());
    }
    // Exercise the public B-rep traversal as well as the reported local solve.
    let mut cache =
        crate::topology::brep_cache::BrepCache::new(crate::geom::a_object::BREP_TYPE).unwrap();
    cache
        .build_trees(
            &k,
            crate::topology::brep_cache::CacheSource::Brep(&b),
            |object| object,
            &arena_identity,
        )
        .unwrap();
    let global_objects = BrepObjects {
        cache: Some(&cache),
        ..objects
    };
    let global = TopologySolver::brep_point_solve(
        &k,
        &global_objects,
        Some(CacheObject::Brep),
        query,
        crate::intersect::topology_solver::p2::TopologySolveOptions {
            operation: SolverOperation::ProjectedMaximize,
            requested: SolutionRequested::Single,
            tolerance_3d: tol,
            best_answer_so_far_sq: 0.0,
            opt_vectors: &[normal],
        },
    )
    .unwrap();
    assert_eq!(global.len(), 1);
    assert!(
        (global.as_slice()[0].start.solution_value - expected).abs() <= tol,
        "global={:?}, normal={normal:?}, expected={expected}",
        global.as_slice()
    );
}

fn check_maximum_after_near_seam_vertex(operation: SolverOperation) {
    let mut k = Kernel::default();
    let mut torus = Torus::new(
        &k,
        Vec3::ZERO,
        Vec3::X,
        Vec3::Y,
        Extent2d::from_bounds(0.0, 0.0, 360.0, 360.0),
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
    let face = create_face_from_surface(&mut k, &mut b, Some(surface), &domain, &mut fenv).unwrap();
    let vertices = crate::topology::brep::p1::vertices(&b, None).unwrap();
    assert_eq!(vertices.len(), 1);
    assert_eq!(
        b.vertex(vertices[0]).unwrap().point,
        Vec3::new(1.25, 0.0, 0.0)
    );
    // Face construction records notification identities for its seam edges;
    // arena-only keys would silently omit them from the cache tree.
    let proven = crate::query::brep::BrepIdentity::of_brep_with_id(BrepId::from_index(0), &b);
    let identity = proven.identity();
    let mut preparation = crate::query::brep::RayPreparation::default();
    {
        let objects = BrepObjects {
            preparation: Some(&preparation),
            brep: &b,
            brep_id: BrepId::from_index(0),
            cache: None,
            keys: identity.keys,
            owner_of: identity.owner_of,
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
    let mut cache =
        crate::topology::brep_cache::BrepCache::new(crate::geom::a_object::BREP_TYPE).unwrap();
    cache
        .build_trees(
            &k,
            crate::topology::brep_cache::CacheSource::Brep(&b),
            |o| o,
            identity.keys,
        )
        .unwrap();
    assert_eq!(cache.curve_tree().unwrap().object_list().unwrap().len(), 2);
    assert_eq!(
        cache.surface_tree().unwrap().object_list().unwrap().len(),
        1
    );
    let objects = BrepObjects {
        preparation: Some(&preparation),
        brep: &b,
        brep_id: BrepId::from_index(0),
        cache: Some(&cache),
        keys: identity.keys,
        owner_of: identity.owner_of,
    };
    let tolerance = 1e-6;
    let offset = if operation == SolverOperation::DirectedMaximize {
        tolerance / 2.0
    } else {
        0.0
    };
    let query = Vec3::new(1.25 + offset, 0.0, 0.0);
    // Every torus point has norm at most R+r=1.25. The triangle bound
    // |C-query| <= 1.25+|query| is attained at (-1.25,0,0). That point
    // also lies on the directed X line in the Z projection, so the same
    // closed form applies to both operations.
    let expected = 2.5 + offset;
    let vectors = [Vec3::Z, Vec3::X];
    let answer = TopologySolver::brep_point_solve(
        &k,
        &objects,
        Some(CacheObject::Brep),
        query,
        crate::intersect::topology_solver::p2::TopologySolveOptions {
            operation,
            requested: SolutionRequested::Single,
            tolerance_3d: tolerance,
            best_answer_so_far_sq: 0.0,
            opt_vectors: &vectors,
        },
    )
    .unwrap();
    assert_eq!(answer.len(), 1);
    let actual = answer.as_slice()[0].start.solution_value;
    assert!(
        (actual - expected).abs() <= tolerance,
        "{operation:?}: actual={actual}, expected={expected}"
    );
}

#[test]
fn distance_maximum_does_not_stop_at_coincident_seam_vertex() {
    check_maximum_after_near_seam_vertex(SolverOperation::Maximize);
}

#[test]
fn directed_maximum_does_not_stop_at_nearby_seam_vertex() {
    check_maximum_after_near_seam_vertex(SolverOperation::DirectedMaximize);
}

fn check_zero_slack_plane_candidate(normal: Vec3) {
    let mut k = Kernel::default();
    let domain = Extent2d::from_bounds(-0.5, -0.5, 0.5, 0.5);
    let surface = &Surface::from(crate::surface::plane::Plane::new(
        Vec3::ZERO,
        Vec3::X,
        Vec3::Y,
        Vec2::new(1.0, 1.0),
        &domain,
    ));
    let expected = 2.5_f64.sqrt() * normal.length_squared();
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
    let query = Vec3::new(-1.0, 0.0, 2.0);
    let tol = 0.0;
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
        .unwrap_or_else(|error| {
            panic!(
                "{error:?}; candidates={:?}",
                solver.base.solutions.as_slice()
            )
        });
    assert_eq!(solver.base.solutions.len(), 1);
    let solution = &solver.base.solutions.as_slice()[0];
    let uv = Vec2::new(
        solution.start.param(0).unwrap(),
        solution.start.param(1).unwrap(),
    );
    let point = surface.evaluate_point(&k, uv).unwrap();
    let measured = normal.cross(point - query).cross(normal).length();
    // Closed form: the farthest projected corners are (0.5,+/-0.5,0).
    // Candidate construction includes the exact trim endpoints. The explicit
    // caller slack is zero; all remaining room comes from arithmetic bounds.
    assert_eq!(solution.objects[0], Some(CacheObject::Face(face)));
    assert!(domain.contains_point(uv, 0.0));
    assert!((solution.start.solution_value - expected).abs() <= 8.0 * f64::EPSILON * expected);
    assert!((measured - expected).abs() <= 8.0 * f64::EPSILON * expected);
    let net = crate::blend::fillet_section::station::Nurbs::of(
        &crate::blend::fillet_section::station::as_b_spline(&k, surface)
            .unwrap()
            .unwrap(),
    )
    .unwrap();
    let enclosure = super::projected_bound::value_at(&net, uv, query, normal).unwrap();
    let perturbation = 4.0 * (enclosure.hi - enclosure.lo);
    assert!(perturbation > 0.0);
    for delta in [-perturbation, perturbation] {
        let mut incorrect = solver.base.clone();
        let mut changed = *solution;
        changed.start.solution_value += delta;
        incorrect.solutions.set_at(0, &changed).unwrap();
        assert!(super::projected_bound::certify_candidate(
            &k,
            &objects,
            CacheObject::Face(face),
            query,
            &incorrect,
        )
        .is_err());
        incorrect.tolerance_3d = 10.0 * expected;
        assert!(super::projected_bound::certify_candidate(
            &k,
            &objects,
            CacheObject::Face(face),
            query,
            &incorrect,
        )
        .is_err());
    }
    // An edge incumbent must supply its own verified point metric. Merely
    // trusting its reported value could "certify" this face against an invented
    // larger answer on another object.
    let edge = crate::topology::brep::p1::edges(&b, None).unwrap()[0];
    let edge_object = CacheObject::Edge(edge);
    let curve = objects
        .curve(objects.edge_curve(edge_object).unwrap().unwrap())
        .unwrap();
    let t = objects.edge_interval(edge_object).unwrap().min();
    let edge_point = curve.evaluate_point(&k, t).unwrap();
    let mut incumbent = solver.base.clone();
    let mut edge_answer = *solution;
    edge_answer.objects.fill(None);
    edge_answer.objects[0] = Some(edge_object);
    edge_answer.num_objects = 1;
    edge_answer.num_variables = 1;
    edge_answer.start.parameters[0] = t;
    edge_answer.start.solution_value = normal.cross(edge_point - query).cross(normal).length();
    incumbent.solutions.set_at(0, &edge_answer).unwrap();
    super::projected_bound::certify_candidate(
        &k,
        &objects,
        CacheObject::Face(face),
        query,
        &incumbent,
    )
    .unwrap();
    edge_answer.start.solution_value = 10.0 * expected;
    incumbent.tolerance_3d = 100.0 * expected;
    incumbent.solutions.set_at(0, &edge_answer).unwrap();
    assert!(super::projected_bound::certify(
        &k,
        &objects,
        CacheObject::Face(face),
        query,
        &incumbent,
    )
    .is_err());

    // Retention slack must not make a feasible but suboptimal point a maximum.
    // For these normals parallel to Z, the center has projected distance
    // |normal|^2, below sqrt(2.5)|normal|^2.
    let mut suboptimal = solver.base.clone();
    suboptimal.tolerance_3d = 10.0 * expected;
    let mut changed = *solution;
    changed.start.parameters[0] = 0.0;
    changed.start.parameters[1] = 0.0;
    changed.start.solution_value = normal.length_squared();
    suboptimal.solutions.set_at(0, &changed).unwrap();
    super::projected_bound::certify_candidate(
        &k,
        &objects,
        CacheObject::Face(face),
        query,
        &suboptimal,
    )
    .unwrap();
    assert!(super::projected_bound::certify(
        &k,
        &objects,
        CacheObject::Face(face),
        query,
        &suboptimal,
    )
    .is_err());
}

#[test]
fn zero_slack_projected_face_maximum_accepts_plane_corner() {
    check_zero_slack_plane_candidate(Vec3::Z);
    check_zero_slack_plane_candidate(Vec3::Z * 2.0);
}
