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
use crate::topology::brep::p2::ArenaNotifyState;
use crate::topology::brep::p3::new_brep;
use crate::topology::brep::p8::{create_face_from_surface, FaceBuildEnv};
use crate::topology::brep_cache::{arena_identity, BrepObjects, CacheObject};
use crate::topology::topology_traverser::Marks;
use crate::topology::BrepId;

#[test]
fn projected_sphere_face_maximum_matches_surface_and_closed_form() {
    let mut k = Kernel::default();
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
    let surface = Surface::from(sphere);
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
    let normal = Vec3::new(0.0, 0.0, 1.0);
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
    let expected = 0.5 + (query.x * query.x + query.y * query.y).sqrt();
    assert!(
        (solution.start.solution_value - expected).abs() <= tol,
        "value={}, expected={expected}",
        solution.start.solution_value
    );
    assert!(
        (solution.start.solution_value - measured).abs() <= 1e-12,
        "value={}, surface measurement={measured}",
        solution.start.solution_value
    );
}
