//! Two rounds that meet at a knife-edge tip are refused by name (SolidCraft #11360 line 2).
//!
//! The arc part: a 41 1/6 x 21 rectangle and, beside it, what a circle of radius 21 about
//! (62 1/6, -21) leaves against the top line it touches, extruded 10. The arc's top edge and
//! the back top edge both leave the tip (62 1/6, 0, 10) heading -x: the top face between them
//! pinches to no width there. Rounded together at radius 2, the rounds' contact lines on the
//! top face cross 13.0 short of the tip (at x = 62 1/6 - sqrt(168)), and at any radius r they
//! cross before it, so no corner fits at the tip. The session refuses that corner with
//! `NoRoom` and says why, instead of failing later on the corner it never made
//! (`SM_ERR_NULL_POINTER`, with no status noted).
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)] // test helpers

use solidcore::blend::fillet_executive::session::FilletSession;
use solidcore::blend::fillet_geom::p1::FilletErrorType;
use solidcore::blend::fillet_solver::p1::FilletSurfaceGeneratorType;
use solidcore::curve::curve::Curve;
use solidcore::geom::global_solver::SolverOperation;
use solidcore::math::{Point3, Vec3, BIG_DOUBLE};
use solidcore::primitive::api_curves::{create_arc, create_line_segment};
use solidcore::query::api_queries::get_edges;
use solidcore::sweep::api_brep::create_linear_sweep;
use solidcore::topology::{edge, EdgeId};
use solidcore::{Brep, Kernel};

const X0: f64 = 41.0 + 1.0 / 6.0;
/// The arc's axis.
const CX: f64 = X0 + 21.0;

/// The arc part: counterclockwise (0, 0), (0, -21), (X0, -21), the arc about (CX, -21) up to
/// its tip (CX, 0), back along the top line; extruded 10 along z.
fn arc_part(k: &mut Kernel) -> Brep {
    let p = |x: f64, y: f64| Point3::new(x, y, 0.0);
    let mut arc = create_arc(Point3::new(CX, -21.0, 0.0), 21.0, 90.0, 180.0).unwrap();
    let ivl = arc.natural_interval();
    arc.reverse_parameterization(&ivl).unwrap();
    let curves = [
        Curve::from(create_line_segment(p(0.0, 0.0), p(0.0, -21.0)).unwrap()),
        Curve::from(create_line_segment(p(0.0, -21.0), p(X0, -21.0)).unwrap()),
        Curve::from(arc),
        Curve::from(create_line_segment(p(CX, 0.0), p(0.0, 0.0)).unwrap()),
    ];
    let refs: Vec<Option<&Curve>> = curves.iter().map(Some).collect();
    create_linear_sweep(k, &refs, Vec3::new(0.0, 0.0, 1.0), 10.0, true).unwrap()
}

/// The edge of `b` through `at`.
fn pick(k: &Kernel, b: &Brep, at: Point3) -> EdgeId {
    let mut best: Option<(EdgeId, f64)> = None;
    for e in get_edges(b).unwrap() {
        let Some(c) = edge::curve(b, e).unwrap() else {
            continue;
        };
        let ivl = edge::interval(b, e).unwrap();
        let Ok(d) = c.drop_point(
            k,
            &ivl,
            at,
            None,
            BIG_DOUBLE,
            None,
            SolverOperation::Minimize,
        ) else {
            continue;
        };
        if d.success && d.distance < 0.1 && best.is_none_or(|(_, bd)| d.distance < bd) {
            best = Some((e, d.distance));
        }
    }
    best.expect("an edge through the point").0
}

/// The arc's top edge and the back top edge rounded together at `radius` under `k`'s fixes:
/// the session's status and comment when it fails.
fn arc_and_back(k: &mut Kernel, radius: f64) -> Option<(i32, String)> {
    let mut b = arc_part(k);
    let a = std::f64::consts::FRAC_PI_4;
    let arc_top = Point3::new(CX - 21.0 * a.cos(), -21.0 + 21.0 * a.sin(), 10.0);
    let edges = [
        pick(k, &b, arc_top),
        pick(k, &b, Point3::new(30.0, 0.0, 10.0)),
    ];
    let xs = [FilletSurfaceGeneratorType::Circular; 2];
    let mut s = FilletSession::new(k, &mut b).unwrap();
    s.prepare_per_edge_constant_radius(k, &edges, &[radius; 2], &xs, 1, 1.0)
        .unwrap();
    let failed = s.create_corners(k, None).is_err() || s.do_filleting(k, None).is_err();
    failed.then(|| {
        let info = s.error_info();
        (info.code, info.comment.clone())
    })
}

fn assert_refused_at_the_tip(k: &mut Kernel) {
    for radius in [2.0, 0.5, 0.01] {
        let Some((code, comment)) = arc_and_back(k, radius) else {
            panic!("radius {radius}: the rounds meeting at the knife-edge tip were built");
        };
        assert_eq!(
            code,
            FilletErrorType::NoRoom as i32,
            "radius {radius}: {comment}"
        );
        assert!(comment.contains("zero angle"), "radius {radius}: {comment}");
    }
}

#[test]
fn rounds_meeting_at_a_knife_edge_tip_are_refused_by_name() {
    assert_refused_at_the_tip(&mut Kernel::default());
}

#[test]
fn rounds_meeting_at_a_knife_edge_tip_are_refused_by_name_with_the_blend_corner() {
    let mut k = Kernel::default();
    k.fixes.fillet_blend_corner = true;
    assert_refused_at_the_tip(&mut k);
}
