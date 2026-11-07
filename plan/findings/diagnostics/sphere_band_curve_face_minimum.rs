//! Public curve/face minimum query plus an independently feasible native witness.
//! Matches static-analytic-faces.json: analytic_face_sphere_band_6_3.
use solidcore::curve::curve::Curve;
use solidcore::geom::{
    extent1d::Extent1d,
    extent2d::Extent2d,
    global_solver::{SolutionRequested, SolverOperation},
    tol::USE_DEFAULT,
};
use solidcore::intersect::topology_solver::p1::{
    CurveRef, LocalSolver, LocalSolverKind, TopologySolver,
};
use solidcore::math::{Vec2, Vec3, BIG_DOUBLE};
use solidcore::nurbs::{
    b_spline_curve::BSplineCurve,
    crv_basic::Curve as Net,
    math::{CPoint, NOW},
};
use solidcore::surface::{
    sphere::{Sphere, SphereOptions},
    surface::Surface,
};
use solidcore::topology::{
    brep::{
        p2::ArenaNotifyState,
        p3::new_brep,
        p8::{create_face_from_surface, FaceBuildEnv},
    },
    brep_cache::{arena_identity, CacheObject},
    topology_traverser::Marks,
};
use solidcore::Kernel;

/// Run the ordinary minimum request and measure a smaller admissible pair.
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut kernel = Kernel::default();
    let mut sphere = Sphere::new(
        &kernel,
        Vec3::new(0., 0., -1.),
        Vec3::X,
        Vec3::Y,
        Extent2d::from_bounds(0., -60., 360., 60.),
        0.5,
        SphereOptions::default(),
    )?;
    sphere.make_nurb(&kernel)?;
    let surface = Surface::from(sphere);
    let domain = surface.natural_uv_domain_with_nurb(&kernel)?;
    let mut brep = new_brep(&mut kernel, USE_DEFAULT, USE_DEFAULT)?;
    let mut state = ArenaNotifyState::new(&brep);
    let mut environment = state.env(None);
    let mut marks = Marks::default();
    let keys = |face| arena_identity(CacheObject::Face(face));
    let mut face_environment = FaceBuildEnv {
        env: &mut environment,
        marks: &mut marks,
        this_face: &keys,
    };
    let face = create_face_from_surface(
        &mut kernel,
        &mut brep,
        Some(surface.clone()),
        &domain,
        &mut face_environment,
    )?;
    let proven = solidcore::query::brep::BrepIdentity::of_brep(&brep);
    let identity = proven.identity();
    let objects = solidcore::query::brep::LiveBrepObjects::new(&mut brep, &identity, None)?;
    let curve = Curve::from(BSplineCurve::from_nurb(
        Net::new(
            vec![CPoint::new(-1., 0., 0., NOW), CPoint::new(1., 0., 0., NOW)],
            1,
            vec![0., 0., 1., 1.],
        )?,
        3,
        false,
    ));
    let normal = [Vec3::Z, Vec3::X];
    let mut solver = TopologySolver::new(
        &kernel,
        &objects,
        SolverOperation::Minimize,
        SolutionRequested::Single,
        1e-5,
        BIG_DOUBLE,
        &normal,
    );
    let mut local = LocalSolver::new(LocalSolverKind::CurveFace, false, None);
    local.set_curve(
        Some(CurveRef {
            curve: &curve,
            handle: CacheObject::Other(arena_identity(CacheObject::Brep)),
        }),
        Extent1d::from_min_max(0., 1.),
    );
    local.solve_it(&kernel, &mut solver, None, Some(CacheObject::Surface(face)))?;
    let uv = Vec2::new(0., 1.);
    let surface_point = surface.evaluate_point(&kernel, uv)?;
    let t = (surface_point.x + 1.) / 2.;
    let curve_point = curve.evaluate_point(&kernel, t)?;
    let distance = (curve_point - surface_point).length();
    assert!((0.0..=1.0).contains(&t) && domain.contains_point(uv, 0.));
    assert_eq!(solver.base.solutions.len(), 1);
    for answer in solver.base.solutions.as_slice() {
        println!(
            "reported={} parameters={:?}",
            answer.start.solution_value,
            &answer.start.parameters[..3]
        );
    }
    println!("witness={distance} t={t} uv={uv:?} curve={curve_point:?} surface={surface_point:?}");
    assert!(solver
        .base
        .solutions
        .as_slice()
        .iter()
        .all(|answer| answer.start.solution_value > distance));
    Ok(())
}
