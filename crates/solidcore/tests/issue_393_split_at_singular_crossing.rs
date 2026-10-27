//! tangency-labs/solidcore#393 (WR-491): two equal cylinders of radius 1 with perpendicular
//! axes, `x^2 + y^2 = 1` (axis z, height -2..2) and `y^2 + z^2 = 1` (axis x, length -2..2).
//! They meet in the two ellipses `x = z` and `x = -z`, which cross at the singular points
//! (0, 1, 0) and (0, -1, 0), where the surface normals are parallel. SMLib's tracer marches
//! through such a point without ending a curve there: it traces the ellipses from the
//! cylinder seams, two curves through (0, -1, 0) and none ending at it, so the points are
//! complete but the topology is not. The opt-in
//! `Kernel::fixes.ssi_split_curves_at_singular_points` splits the curves there.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers
use solidcore::curve::curve::Curve;
use solidcore::geom::tol::approx_tol3d;
use solidcore::geom::vector3d::Vec3;
use solidcore::surface::cone::Cone;
use solidcore::surface::surface::p2::SurfaceIntersection;
use solidcore::surface::surface::Surface;
use solidcore::Kernel;

/// The accuracy asked of the intersection (the corpus' `--tol 1e-6`).
const TOL: f64 = 1e-6;
const CROSSINGS: [[f64; 3]; 2] = [[0.0, 1.0, 0.0], [0.0, -1.0, 0.0]];

/// The cylinder of radius 1 and height 4 whose axis runs from `origin` along `x` cross `y`.
fn cylinder(k: &Kernel, origin: [f64; 3], x: [f64; 3], y: [f64; 3]) -> Surface {
    let mut c = Cone::new_cylinder_bounded(
        Vec3::new(origin[0], origin[1], origin[2]),
        Vec3::new(x[0], x[1], x[2]),
        Vec3::new(y[0], y[1], y[2]),
        1.0,
        0.0,
        360.0,
        4.0,
        false,
        false,
        false,
    );
    c.make_nurb(k).unwrap();
    Surface::from(c)
}

fn steinmetz(k: &Kernel) -> SurfaceIntersection {
    let a = cylinder(k, [0.0, 0.0, -2.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
    let b = cylinder(k, [-2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]);
    let (da, db) = (a.natural_uv_domain(), b.natural_uv_domain());
    a.global_surface_intersect(
        k,
        &da,
        &b,
        &db,
        [true, true],
        Some(approx_tol3d(k, TOL)),
        None,
    )
    .unwrap()
}

/// The two end points of a 3D curve.
fn ends(k: &Kernel, c: &Curve) -> [Vec3; 2] {
    let ivl = c.natural_interval();
    [ivl.min(), ivl.max()].map(|t| c.evaluate_point(k, t).unwrap())
}

fn near(p: Vec3, q: [f64; 3], d: f64) -> bool {
    p.distance(Vec3::new(q[0], q[1], q[2])) < d
}

/// How many curve ends of `x` lie within `d` of `q`. Four branches meet at a crossing, so
/// four curves end at one that is identified.
fn curve_ends_at(k: &Kernel, x: &SurfaceIntersection, q: [f64; 3], d: f64) -> usize {
    x.curves_3d
        .iter()
        .flat_map(|c| ends(k, c))
        .filter(|&e| near(e, q, d))
        .count()
}

/// Every sampled point of every curve lies on both cylinders (the curves are right).
fn assert_on_both_cylinders(k: &Kernel, x: &SurfaceIntersection) {
    for c in &x.curves_3d {
        let ivl = c.natural_interval();
        for i in 0..=20 {
            let t = ivl.min() + ivl.length() * f64::from(i) / 20.0;
            let p = c.evaluate_point(k, t).unwrap();
            assert!(
                (p.x * p.x + p.y * p.y - 1.0).abs() < 1e-5,
                "{p:?} off x^2+y^2=1"
            );
            assert!(
                (p.y * p.y + p.z * p.z - 1.0).abs() < 1e-5,
                "{p:?} off y^2+z^2=1"
            );
        }
    }
}

#[test]
fn issue_393_split_is_off_by_default() {
    assert!(!Kernel::default().fixes.ssi_split_curves_at_singular_points);
}

#[test]
fn issue_393_default_ends_no_curve_at_a_crossing() {
    let k = Kernel::default();
    let x = steinmetz(&k);
    assert!(!x.curves_3d.is_empty());
    assert_on_both_cylinders(&k, &x);
    // as upstream: the ellipses are traced from the seams round to (0, 1, 0), where the
    // four branches end, and through (0, -1, 0) without a curve ending there
    assert_eq!(x.curves_3d.len(), 4);
    assert_eq!(curve_ends_at(&k, &x, CROSSINGS[0], 1e-6), 4);
    assert_eq!(curve_ends_at(&k, &x, CROSSINGS[1], 1e-3), 0);
}

#[test]
fn issue_393_split_ends_a_curve_at_each_crossing() {
    let mut k = Kernel::default();
    k.fixes.ssi_split_curves_at_singular_points = true;
    let x = steinmetz(&k);
    assert_on_both_cylinders(&k, &x);
    for q in CROSSINGS {
        assert_eq!(
            curve_ends_at(&k, &x, q, 1e-6),
            4,
            "four branches meet at {q:?}"
        );
    }
    // the two long arcs are split at (0, -1, 0)
    assert_eq!(x.curves_3d.len(), 6);
    // the arrays stay parallel, and the UV curves are split with the 3D curve
    let n = x.curves_3d.len();
    assert_eq!(x.surface1_uv_curves.len(), n);
    assert_eq!(x.surface2_uv_curves.len(), n);
    assert_eq!(x.curve_types.len(), n);
    assert_eq!(x.deviations.len(), n);
    for (i, c) in x.curves_3d.iter().enumerate() {
        let ivl = c.natural_interval();
        for uv in [&x.surface1_uv_curves[i], &x.surface2_uv_curves[i]] {
            let u = uv.as_ref().unwrap().natural_interval();
            assert!((u.min() - ivl.min()).abs() < 1e-12 && (u.max() - ivl.max()).abs() < 1e-12);
        }
    }
}
