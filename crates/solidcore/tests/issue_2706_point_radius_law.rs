//! #2706: a variable round whose B-spline radius law is a single point.
//!
//! `FilletLaw::bspline` reads the radius as the x coordinate of its law curve. A law whose
//! control points are all equal, here degree 1 from (3, 0, 0) to (3, 0, 0), is a curve of one
//! point: its x is 3 everywhere, so the law is the constant radius 3. The round of the
//! 10 x 20 x 30 box's edge x = 10, z = 30 (along y, 20 long) with it must be the constant round
//! of 3, which removes `(1 - pi/4) 3^2 20` in closed form.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use solidcore::blend::fillet_executive::session::FilletSession;
use solidcore::blend::fillet_standard_solver::p1::FilletLaw;
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::geom::extent1d::Extent1d;
use solidcore::mass::integrator::properties::{brep_volume, PropertyAttributes};
use solidcore::math::{Point3, Vec3};
use solidcore::nurbs::b_spline_curve::p1::{
    BSplineCurve, BSplineCurveForm, CanonicalData, KnotType,
};
use solidcore::primitive::api_primitives::create_box;
use solidcore::topology::brep::validation::{check_brep, CheckLevel};
use solidcore::topology::{edge, EdgeId};
use solidcore::{Brep, Kernel};

fn volume(k: &Kernel, b: &mut Brep) -> f64 {
    let store = AttributeStore::default();
    let attrs = PropertyAttributes {
        store: &store,
        key_of: &|_| ObjectKey(0),
    };
    brep_volume(k, b, 1e-8, &attrs).unwrap()
}

/// The edge from `p` to `q` (either way).
fn edge_between(k: &Kernel, b: &Brep, p: Point3, q: Point3) -> EdgeId {
    b.edges()
        .ids()
        .find(|&e| {
            let ivl = edge::interval(b, e).unwrap();
            let c = edge::curve(b, e).unwrap().unwrap();
            let s = c.evaluate_point(k, ivl.min()).unwrap();
            let t = c.evaluate_point(k, ivl.max()).unwrap();
            (s.distance(p) < 1e-9 && t.distance(q) < 1e-9)
                || (s.distance(q) < 1e-9 && t.distance(p) < 1e-9)
        })
        .expect("edge")
}

/// The degree 1 law from `a` to `b` over the edge's map [0, 1].
fn linear_law(a: Vec3, b: Vec3) -> FilletLaw {
    let curve = BSplineCurve::create_canonical(&CanonicalData {
        dim: 3,
        degree: 1,
        control_points: &[a, b],
        form: BSplineCurveForm::default(),
        knot_multiplicities: &[2, 2],
        knots: &[0.0, 1.0],
        knot_type: KnotType::default(),
        weights: None,
        trim: None,
    })
    .unwrap();
    FilletLaw::bspline(curve, &Extent1d::from_min_max(0.0, 1.0))
}

/// The box's volume after its edge x = 10, z = 30 is rounded with `law`, checked a closed
/// manifold solid with a clean Level 2 check.
fn rounded_box(law: FilletLaw) -> f64 {
    let mut k = Kernel::default();
    let mut b = create_box(&mut k, Vec3::ZERO, 10.0, 20.0, 30.0).unwrap();
    let e = edge_between(
        &k,
        &b,
        Point3::new(10.0, 0.0, 30.0),
        Point3::new(10.0, 20.0, 30.0),
    );
    {
        let mut session = FilletSession::new(&mut k, &mut b).unwrap();
        session
            .prepare_variable_radius_edges(&mut k, &[e], 3.0, law, 1, 1)
            .unwrap();
        session.do_fillet(&mut k).unwrap();
    }
    assert!(solidcore::topology::brep::p2::is_manifold_solid(&b).unwrap());
    let lines = check_brep(&mut k, &mut b, CheckLevel::Level2)
        .unwrap()
        .lines();
    assert!(lines.is_empty(), "{lines:?}");
    volume(&k, &mut b)
}

/// The constant round of 3 removes `(1 - pi/4) 9 20` of the box's 6000.
fn constant_removed() -> f64 {
    (1.0 - core::f64::consts::FRAC_PI_4) * 9.0 * 20.0
}

#[test]
fn a_single_point_law_is_the_constant_round() {
    let law = linear_law(Vec3::new(3.0, 0.0, 0.0), Vec3::new(3.0, 0.0, 0.0));
    let v = rounded_box(law);
    let removed = 6000.0 - v;
    let want = constant_removed();
    eprintln!("point law removes {removed:.12}, want {want:.12}");
    assert!(
        (removed - want).abs() <= 1e-6 * want,
        "removed {removed}, want {want}"
    );
}
