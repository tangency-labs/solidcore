//! A curve's length holds the accuracy it is asked for, from 1e-3 to 1e-10, against its
//! closed form. Each span of the curve's cache (within one knot span) is the adaptive
//! Cash–Karp integral of `|C'|`:
//! - each span's end was evaluated in the knot span above it (upstream's
//!   `bFromLeft = TRUE`). Where the speed jumps at a C0 knot (the double knot between
//!   unequal rational arcs) the span's last sample was the next span's value; its
//!   steps shrank toward the knot and underflowed (-0.19 at 1e-10: a 70 degree span
//!   failed, its length taken as zero). The end is now evaluated in its own span, as the
//!   area integrals' iso-u pieces (WR-670).
//! - its steps were accepted on the embedded error estimate alone, a null rule of
//!   degree 3 that vanishes where `g''''` changes sign; they are now accepted on the
//!   larger of it and the fifth divided difference (`OdeOptions::quadrature`, as the
//!   precise-property integrals).
#![allow(clippy::unwrap_used, clippy::indexing_slicing, clippy::panic)]
use solidcore::curve::curve::Curve;
use solidcore::geom::extent1d::Extent1d;
use solidcore::nurbs::b_spline_curve::BSplineCurve;
use solidcore::nurbs::math::{CPoint, NOW};
use solidcore::Kernel;

const RADIUS: f64 = 5.0;

/// A rational quadratic B-spline of `RADIUS` through the arcs `arcs_deg` in turn, one
/// arc a unit knot span, joined at double knots: G2, but where two arcs differ its
/// speed jumps (C0 in its parameter).
fn arcs(arcs_deg: &[f64]) -> Curve {
    let mut start = 0.0_f64;
    let mut cpw = vec![CPoint::new(RADIUS, 0.0, 0.0, 1.0)];
    let mut knots = vec![0.0; 3];
    for (i, deg) in arcs_deg.iter().enumerate() {
        let half = deg.to_radians() / 2.0;
        let mid = start + half;
        let end = start + 2.0 * half;
        let w = half.cos();
        // the middle point at RADIUS / w, weighted by w
        cpw.push(CPoint::new(RADIUS * mid.cos(), RADIUS * mid.sin(), 0.0, w));
        cpw.push(CPoint::new(
            RADIUS * end.cos(),
            RADIUS * end.sin(),
            0.0,
            1.0,
        ));
        let k = (i + 1) as f64;
        knots.extend([k, k]);
        start = end;
    }
    knots.push(arcs_deg.len() as f64);
    Curve::from(
        BSplineCurve::from_data(3, cpw.len() - 1, knots.len() - 1, 2, &cpw, &knots).unwrap(),
    )
}

fn holds(curve: &Curve, exact: f64) {
    let k = Kernel::default();
    let interval = curve.natural_interval();
    for e in 3..=10 {
        let accuracy = 10f64.powi(-e);
        let length = curve.length(&k, &interval, accuracy).unwrap();
        let rel = (length - exact) / exact;
        assert!(
            rel.abs() <= accuracy,
            "at accuracy {accuracy:e}: length {length} is off by {rel:e}"
        );
    }
}

/// One 120 degree arc: a single smooth span.
#[test]
fn a_smooth_arc_holds_every_requested_accuracy() {
    holds(&arcs(&[120.0]), RADIUS * 120f64.to_radians());
}

/// The parabola `(t, t^2)` for `t` in `[0, 10]`, one quadratic span: its speed
/// `sqrt(1 + 4 t^2)` turns from flat to linear near the vertex.
#[test]
fn a_parabola_holds_every_requested_accuracy() {
    let a = 10.0_f64;
    let cpw = [
        CPoint::new(0.0, 0.0, 0.0, 1.0),
        CPoint::new(a / 2.0, 0.0, 0.0, 1.0),
        CPoint::new(a, a * a, 0.0, 1.0),
    ];
    let knots = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
    let curve = Curve::from(BSplineCurve::from_data(3, 2, 5, 2, &cpw, &knots).unwrap());
    let exact = a * (1.0 + 4.0 * a * a).sqrt() / 2.0 + (2.0 * a).asinh() / 4.0;
    holds(&curve, exact);
}

/// A full circle of six unequal arcs, its speed jumping at each of five double knots.
#[test]
fn a_circle_of_unequal_arcs_holds_every_requested_accuracy() {
    holds(
        &arcs(&[30.0, 50.0, 70.0, 40.0, 80.0, 90.0]),
        RADIUS * 2.0 * std::f64::consts::PI,
    );
}

/// The knots of [`cubic`]: two knot spans, `[0, 0.4]` and `[0.4, 1]`.
const CUBIC_KNOTS: [f64; 9] = [0.0, 0.0, 0.0, 0.0, 0.4, 1.0, 1.0, 1.0, 1.0];

/// A cubic of two knot spans through the control points `points`, weighted by `w`
/// (stored homogeneous) or non-rational for `None`.
fn cubic(points: &[[f64; 3]; 5], w: Option<[f64; 5]>) -> Curve {
    let cpw: Vec<CPoint> = match w {
        Some(w) => points
            .iter()
            .zip(w)
            .map(|(p, w)| CPoint::new(p[0] * w, p[1] * w, p[2] * w, w))
            .collect(),
        None => points
            .iter()
            .map(|p| CPoint::new(p[0], p[1], p[2], NOW))
            .collect(),
    };
    Curve::from(BSplineCurve::from_data(3, 4, 8, 3, &cpw, &CUBIC_KNOTS).unwrap())
}

/// Lengths at 1e-3 .. 1e-12, by default and with failures propagated.
fn lengths(curve: &Curve) -> Vec<(f64, bool, f64)> {
    let mut out = Vec::new();
    for propagate in [false, true] {
        let mut k = Kernel::default();
        k.fixes.integrator_propagates_failure = propagate;
        for e in 3..=12 {
            let accuracy = 10f64.powi(-e);
            let length = curve
                .length(&k, &curve.natural_interval(), accuracy)
                .unwrap_or_else(|err| panic!("at {accuracy:e}, propagate {propagate}: {err:?}"));
            out.push((accuracy, propagate, length));
        }
    }
    out
}

/// A curve whose control points are all one point has length 0, exactly, at every
/// accuracy and whether failures propagate or not. Its computed speed is rounding noise
/// around zero, which the steps used to refuse at 1e-9 (step control 1001: 0 by
/// default, an error with failures propagated). Weights that are powers of two keep the
/// rational points exactly equal.
#[test]
fn a_point_curve_has_length_zero() {
    for m in [1e-3, 1.0, 1e6] {
        let p = [0.3 * m, -0.7 * m, 0.11 * m];
        for w in [None, Some([1.0, 0.5, 2.0, 0.25, 1.0])] {
            for (accuracy, propagate, length) in lengths(&cubic(&[p; 5], w)) {
                assert_eq!(
                    length, 0.0,
                    "|P| ~ {m:e}, weights {w:?}, at {accuracy:e} (propagate {propagate})"
                );
            }
        }
    }
}

/// A curve whose control points differ by a few ulps is a point to rounding: its length
/// (at most its control polygon's, a few ulps) is resolved to the rounding of its
/// computed speed, `speed_rounding_bound` over each knot span, at every accuracy.
/// Rational weights that are not powers of two make the same points differ by rounding.
#[test]
fn a_near_point_curve_has_length_within_its_rounding() {
    for m in [1e-3, 1.0, 1e6] {
        let base: [f64; 3] = [0.3 * m, -0.7 * m, 0.11 * m];
        let mut points = [base; 5];
        for (i, q) in points.iter_mut().enumerate() {
            for _ in 0..i {
                q[0] = q[0].next_up();
                q[2] = q[2].next_down();
            }
        }
        for w in [None, Some([1.0, 0.7, 1.3, 0.9, 1.0])] {
            let curve = cubic(&points, w);
            let polygon: f64 = points
                .iter()
                .zip(points.iter().skip(1))
                .map(|(a, b)| {
                    ((b[0] - a[0]).powi(2) + (b[1] - a[1]).powi(2) + (b[2] - a[2]).powi(2)).sqrt()
                })
                .sum();
            let nurb = curve.b_spline_base().unwrap();
            let rounding: f64 = [(0.0, 0.4), (0.4, 1.0)]
                .iter()
                .map(|&(a, b)| {
                    let span = Extent1d::from_min_max(a, b);
                    nurb.speed_rounding_bound(&span).unwrap() * (b - a)
                })
                .sum();
            for (accuracy, propagate, length) in lengths(&curve) {
                assert!(
                    (0.0..=polygon + rounding).contains(&length),
                    "|P| ~ {m:e}, weights {w:?}, at {accuracy:e} (propagate {propagate}): \
                     length {length:e} beyond polygon {polygon:e} + rounding {rounding:e}"
                );
            }
        }
    }
}
