//! WR-80 (the fork's tangency-labs/smlib@0dbe681): `GlobalSurfaceIntersect` failed with
//! `SM_ERR_NULL_POINTER` when `TraceCurve` returned an empty curve (its start point
//! already on a found curve) and the curve was traced back the other way: the empty list
//! was reversed. solidcore traces back only a trace that produced points, so the
//! intersection completes; this test fails on the old path.
//!
//! The case is the SSI corpus' `c2n_quartic_contact_parabola` (jyviko/ssi-corpus
//! `b04c8a6`, Apache-2.0): two NURBS sheets, degree (2, 4), in 4th-order contact along a
//! parabola (`z = 0.5w^3 - 0.5w + 0.2u^2 + w^4` and the same without `w^4`), at
//! approximation tolerance `1e-6`. Both sheets are graphs over `(x, y) = (u, w)`, so they
//! meet where `y^4 = 0`: on the parabola `y = 0`, `z = 0.2 x^2`.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers
use solidcore::curve::curve::Curve;
use solidcore::geom::tol::approx_tol3d;
use solidcore::kernel::Kernel;
use solidcore::nurbs::b_spline_surface::BSplineSurface;
use solidcore::nurbs::math::CPoint;
use solidcore::nurbs::srf_basic::Surface as NurbSurface;
use solidcore::surface::surface::p2::SurfaceIntersection;
use solidcore::surface::surface::Surface;
use solidcore::Point3;

/// A non-rational NURBS sheet of degree (2, 4) from the corpus' 3 x 5 control points.
fn sheet(points: &[[f64; 3]]) -> Surface {
    let net: Vec<Vec<CPoint>> = points
        .chunks(5)
        .map(|row| {
            row.iter()
                .map(|p| CPoint::new(p[0], p[1], p[2], 1.0))
                .collect()
        })
        .collect();
    let sur = NurbSurface::new(
        net,
        2,
        4,
        vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    )
    .unwrap();
    Surface::from(BSplineSurface::from_nurb(sur, false))
}

fn intersect(k: &Kernel) -> solidcore::Result<SurfaceIntersection> {
    let with_quartic = sheet(&[
        [-1.0, -1.0, 1.2],
        [-1.0, -0.5, -0.30000000000000004],
        [-1.0, 0.0, 1.2],
        [-1.0, 0.5, -1.2999999999999998],
        [-1.0, 1.0, 1.1999999999999993],
        [0.0, -1.0, 0.7999999999999999],
        [0.0, -0.5, -0.7000000000000001],
        [0.0, 0.0, 0.7999999999999999],
        [0.0, 0.5, -1.6999999999999997],
        [0.0, 1.0, 0.7999999999999993],
        [1.0, -1.0, 1.2],
        [1.0, -0.5, -0.30000000000000004],
        [1.0, 0.0, 1.2],
        [1.0, 0.5, -1.2999999999999996],
        [1.0, 1.0, 1.1999999999999993],
    ]);
    let without = sheet(&[
        [-1.0, -1.0, 0.2],
        [-1.0, -0.5, 0.7],
        [-1.0, 0.0, 0.19999999999999996],
        [-1.0, 0.5, -0.30000000000000004],
        [-1.0, 1.0, 0.20000000000000018],
        [0.0, -1.0, -0.2],
        [0.0, -0.5, 0.29999999999999993],
        [0.0, 0.0, -0.20000000000000007],
        [0.0, 0.5, -0.7000000000000001],
        [0.0, 1.0, -0.19999999999999984],
        [1.0, -1.0, 0.19999999999999996],
        [1.0, -0.5, 0.7],
        [1.0, 0.0, 0.19999999999999996],
        [1.0, 0.5, -0.30000000000000004],
        [1.0, 1.0, 0.20000000000000018],
    ]);
    let d1 = with_quartic.natural_uv_domain_with_nurb(k).unwrap();
    let d2 = without.natural_uv_domain_with_nurb(k).unwrap();
    let tol = approx_tol3d(k, 1e-6);
    with_quartic.global_surface_intersect(k, &d1, &without, &d2, [true, true], Some(tol), None)
}

/// The curve's ends and middle.
fn samples(k: &Kernel, curve: &Curve) -> Vec<Point3> {
    let interval = curve.natural_interval();
    [interval.min(), interval.mid(), interval.max()]
        .into_iter()
        .map(|t| curve.evaluate_point(k, t).unwrap())
        .collect()
}

/// The intersection completes on the default kernel, and what it finds lies on the
/// contact parabola `y = 0`, `z = 0.2 x^2`.
#[test]
fn wr_80_quartic_contact_completes() {
    let k = Kernel::default();
    let r = match intersect(&k) {
        Ok(r) => r,
        Err(e) => panic!("the intersection fails: {e:?}"),
    };
    assert!(!r.curves_3d.is_empty(), "no curve found");
    // The sheets are `y^4` apart, so a point within the tolerance `1e-6` of both lies
    // within `1e-6^(1/4)` (about 0.032) of the contact line `y = 0`.
    let contact_band = 1e-6_f64.powf(0.25);
    for p in r.curves_3d.iter().flat_map(|c| samples(&k, c)) {
        // on the sheet without the quartic term, which is the graph
        // z = 0.5 y^3 - 0.5 y + 0.2 x^2
        let on_sheet = 0.5 * p.y * p.y * p.y - 0.5 * p.y + 0.2 * p.x * p.x;
        assert!((p.z - on_sheet).abs() < 1e-5, "{p:?} is off the sheet");
        assert!(p.y.abs() < contact_band, "{p:?} is off the contact parabola");
    }
}
