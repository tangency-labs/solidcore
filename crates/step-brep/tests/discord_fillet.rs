//! Claas Kuhnen's Discord comparison of 30 October 2026, supplied by the user.
//!
//! The requested outer top-rim radius was 70 mm. STEP does not retain the command's
//! edge selection; this regression selects all four convex outer arcs. It tests
//! operation completion and the new surfaces, not validity of the complete body:
//! the imported operand already has failed tolerance and geometry reports.

use solidcore::blend::fillet_executive::session::FilletSession;
use solidcore::blend::fillet_solver::p1::FilletSurfaceGeneratorType;
use solidcore::geom::vector2d::Point2;
use solidcore::surface::surface::SurfaceKind;
use solidcore::topology::brep::validation::{check_brep, CheckLevel};
use solidcore::topology::{edge, face::p2 as face, vertex, Brep, EdgeId};
use solidcore::{Kernel, SolidError};

const TOP: f64 = 458.266019821167;
const CYLINDER_RADIUS: f64 = 304.8;
const RADIUS: f64 = 70.0;

fn operand() -> (Kernel, Brep, Vec<EdgeId>) {
    let model = step_import::import(
        include_bytes!("fixtures/discord-2026-10-30-fillet.step"),
        &step_import::Limits::default(),
    )
    .unwrap();
    assert_eq!(model.bodies.len(), 1);
    assert!(model.skipped.is_empty());
    let built =
        step_brep::build_body(&model.bodies[0], step_brep::BuildOptions::default()).unwrap();
    assert!(built.face_failures.is_empty());
    let b = built.brep;
    let k = Kernel::default();
    let mut selected = Vec::new();
    for e in b.edges().ids() {
        let (Some(a), Some(z)) = edge::start_end_vertices(&b, e).unwrap() else {
            continue;
        };
        let a = vertex::point(&b, a).unwrap();
        let z = vertex::point(&b, z).unwrap();
        if (a.z - TOP).abs() > 1e-5 || (z.z - TOP).abs() > 1e-5 {
            continue;
        }
        let curve = edge::curve(&b, e).unwrap().unwrap();
        let interval = edge::interval(&b, e).unwrap();
        let mid = curve.evaluate_point(&k, interval.evaluate(0.5)).unwrap();
        if ((mid.x * mid.x + mid.y * mid.y).sqrt() - CYLINDER_RADIUS).abs() < 1e-3 {
            selected.push(e);
        }
    }
    assert_eq!(selected.len(), 4);
    (k, b, selected)
}

fn round(k: &mut Kernel, b: &mut Brep, edges: &[EdgeId]) -> (solidcore::Result<()>, i32) {
    let mut session = FilletSession::new(k, b).unwrap();
    session
        .prepare_per_edge_constant_radius(
            k,
            edges,
            &vec![RADIUS; edges.len()],
            &vec![FilletSurfaceGeneratorType::Circular; edges.len()],
            1,
            1.0,
        )
        .unwrap();
    let result = session.do_fillet(k);
    (result, session.error_info().code)
}

#[test]
fn pinned_construction_fails_at_the_four_arc_corner_edges() {
    let (mut k, mut b, edges) = operand();
    let (result, code) = round(&mut k, &mut b, &edges);
    assert!(matches!(
        result,
        Err(SolidError::Upstream { code: 1001, .. })
    ));
    assert_eq!(code, 3, "corner boundary computation failure");
}

#[test]
fn exact_classes_build_four_rounds_on_the_expected_seventy_mm_torus() {
    let (mut k, mut b, edges) = operand();
    k.fixes.fillet_exact_round_classes = true;
    let input_check = check_brep(&mut k, &mut b, CheckLevel::Level2).unwrap();
    assert!(b.faces().ids().all(|f| !matches!(
        face::surface(&b, f).unwrap().unwrap().kind(),
        SurfaceKind::SurfOfRevolution(_) | SurfaceKind::Torus(_)
    )));
    round(&mut k, &mut b, &edges).0.unwrap();
    let mut count = 0;
    for f in b.faces().ids() {
        let s = face::surface(&b, f).unwrap().unwrap();
        if !matches!(
            s.kind(),
            SurfaceKind::SurfOfRevolution(_) | SurfaceKind::Torus(_)
        ) {
            continue;
        }
        let d = s.natural_uv_domain_with_nurb(&k).unwrap();
        for i in 0..=10 {
            for j in 0..=10 {
                let u = d.min().x + (d.max().x - d.min().x) * f64::from(i) / 10.0;
                let v = d.min().y + (d.max().y - d.min().y) * f64::from(j) / 10.0;
                let p = s.evaluate_point(&k, Point2::new(u, v)).unwrap();
                let rho = (p.x * p.x + p.y * p.y).sqrt();
                let distance = ((rho - (CYLINDER_RADIUS - RADIUS)).powi(2)
                    + (p.z - (TOP - RADIUS)).powi(2))
                .sqrt();
                assert!(
                    (distance - RADIUS).abs() < 1e-9,
                    "round {f:?} is off its torus: {distance}"
                );
            }
        }
        count += 1;
    }
    assert_eq!(count, 4);
    let output_check = check_brep(&mut k, &mut b, CheckLevel::Level2).unwrap();
    eprintln!(
        "whole-body validation failures: input={}, output={}; surface samples are not a validity certificate",
        input_check.failure_count(), output_check.failure_count()
    );
}
