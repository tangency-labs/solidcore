use super::*;
use crate::geom::vector3d::Vec3;
use crate::nurbs::math::NOW;
use std::f64::consts::{FRAC_1_SQRT_2, FRAC_PI_2};

/// A curve from homogeneous points `(w x, w y, w z, w)`, `rational` false for a polynomial one.
fn curve(points: &[[f64; 4]], p: usize, knots: &[f64], rational: bool) -> BSplineCurve {
    let pts = points
        .iter()
        .map(|&[x, y, z, w]| CPoint::new(x, y, z, if rational { w } else { NOW }))
        .collect();
    let net = Net::new(pts, p, knots.to_vec()).unwrap();
    BSplineCurve {
        nurb: Some(net),
        ..BSplineCurve::new(3)
    }
}

/// The quarter circle of radius 2 about the origin in z = 0, from (2, 0) to (0, 2).
fn quarter() -> BSplineCurve {
    let h = FRAC_1_SQRT_2;
    curve(
        &[
            [2.0, 0.0, 0.0, 1.0],
            [2.0 * h, 2.0 * h, 0.0, h],
            [0.0, 2.0, 0.0, 1.0],
        ],
        2,
        &[0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        true,
    )
}

/// A cubic in the plane z = 0 with two spans.
fn cubic() -> BSplineCurve {
    curve(
        &[
            [0.0, 0.0, 0.0, 1.0],
            [1.0, 1.0, 0.0, 1.0],
            [2.0, 0.5, 0.0, 1.0],
            [3.0, 1.5, 0.0, 1.0],
            [4.0, 1.0, 0.0, 1.0],
        ],
        3,
        &[0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0],
        false,
    )
}

fn bounds(c: &BSplineCurve) -> (f64, f64) {
    c.nurb.as_ref().unwrap().param_bounds()
}

fn at(c: &BSplineCurve, t: f64) -> Vec3 {
    c.evaluate_point(t).unwrap()
}

/// Every source point and knot is in the result, bit for bit, at the kept end.
fn kept(src: &BSplineCurve, out: &BSplineCurve, high: bool) {
    assert!(keeps(
        src.nurb.as_ref().unwrap(),
        out.nurb.as_ref().unwrap(),
        high
    ));
}

#[test]
fn a_line_runs_on_by_the_distance_at_either_end_for_every_continuity() {
    let line = curve(
        &[
            [0.0, 0.0, 0.0, 1.0],
            [1.0, 2.0, 2.0, 1.0],
            [2.0, 4.0, 4.0, 1.0],
        ],
        2,
        &[0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        false,
    );
    let dir = Vec3::new(1.0, 2.0, 2.0) * (1.0 / 3.0);
    for how in [
        Continuity::Natural,
        Continuity::Tangent,
        Continuity::Curvature,
    ] {
        let e = line.extended_by_distance(CurveEnd::End, how, 1.5).unwrap();
        let (_, hi) = bounds(&e.curve);
        let end = at(&e.curve, hi);
        let want = Vec3::new(2.0, 4.0, 4.0) + dir * 1.5;
        assert!((end - want).length() <= 1e-12, "{how:?}: {end:?}");
        assert!((e.across.enclosure.mid() - 1.5).abs() <= 1.5e-9);
        kept(&line, &e.curve, true);

        let s = line
            .extended_by_distance(CurveEnd::Start, how, 0.75)
            .unwrap();
        let (lo, _) = bounds(&s.curve);
        let start = at(&s.curve, lo);
        assert!(
            (start - dir * -0.75).length() <= 1e-12,
            "{how:?}: {start:?}"
        );
        kept(&line, &s.curve, false);
    }
}

#[test]
fn a_circle_continued_stays_on_its_circle_by_the_arc_length() {
    let q = quarter();
    for (end, d) in [(CurveEnd::End, 1.0), (CurveEnd::Start, 2.5)] {
        let e = q.extended_by_distance(end, Continuity::Natural, d).unwrap();
        let (lo, hi) = bounds(&e.curve);
        let p = at(&e.curve, if end == CurveEnd::End { hi } else { lo });
        // on the circle of radius 2, at the arc d / 2 past the quarter's end
        assert!((p.length() - 2.0).abs() <= 1e-13, "{p:?}");
        // its arc past the quarter's end, in closed form, lies in the certified enclosure of
        // the length (to the rounding of atan2), and that is within its contract of d:
        // max(1e-9 d, the enclosure's truncation and rounding parts)
        let angle = p.y.atan2(p.x);
        let arc = 2.0
            * if end == CurveEnd::End {
                angle - FRAC_PI_2
            } else {
                -angle
            };
        let b = e.across.enclosure;
        assert!(
            b.lower - 1e-14 <= arc && arc <= b.upper + 1e-14,
            "{arc} not in [{}, {}]",
            b.lower,
            b.upper
        );
        let parts = e.across.truncation + e.across.rounding;
        assert!(
            (arc - d).abs() <= (1e-9 * d).max(parts),
            "{arc} {d} {parts}"
        );
        // and the curve between too: a point inside the extension
        let mid = at(
            &e.curve,
            if end == CurveEnd::End {
                0.5 * (1.0 + hi)
            } else {
                0.5 * lo
            },
        );
        assert!((mid.length() - 2.0).abs() <= 1e-13);
        kept(&q, &e.curve, end == CurveEnd::End);
    }
}

#[test]
fn a_tangent_extension_is_the_end_piece_mirrored() {
    let c = cubic();
    let e = c
        .extended_by_distance(CurveEnd::End, Continuity::Tangent, 0.5)
        .unwrap();
    kept(&c, &e.curve, true);
    let pn = Vec3::new(4.0, 1.0, 0.0);
    let leg = (pn - Vec3::new(3.0, 1.5, 0.0)) * (1.0 / 1.25_f64.sqrt());
    let (_, hi) = bounds(&e.curve);
    // the extension at 1 + s is the source at 1 - s reflected in the plane through the end
    // point normal to the last leg
    for f in [0.25, 0.5, 1.0] {
        let s = f * (hi - 1.0);
        let x = at(&c, 1.0 - s);
        let mirrored = x - leg * (2.0 * (x - pn).dot(leg));
        assert!((at(&e.curve, 1.0 + s) - mirrored).length() <= 1e-13);
    }
    assert!((e.across.enclosure.mid() - 0.5).abs() <= 0.5e-9);
}

#[test]
fn a_curvature_extension_keeps_the_curvature_across_the_end() {
    let c = cubic();
    let e = c
        .extended_by_distance(CurveEnd::End, Continuity::Curvature, 0.4)
        .unwrap();
    kept(&c, &e.curve, true);
    let net = e.curve.nurb.as_ref().unwrap();
    let left = net
        .derivs_at_knot(1.0, crate::nurbs::knots::Side::Left, 2)
        .unwrap();
    let right = net
        .derivs_at_knot(1.0, crate::nurbs::knots::Side::Right, 2)
        .unwrap();
    let v = |p: &CPoint| Vec3::new(p.x, p.y, p.z);
    let curvature = |d: &[CPoint]| {
        let (d1, d2) = (v(&d[1]), v(&d[2]));
        d1.cross(d2).length() / d1.length().powi(3)
    };
    // the first derivative continuous (reversed by the mirror, reversed back by reading it
    // backwards), the curvature equal
    assert!((v(&left[1]) - v(&right[1])).length() <= 1e-12);
    assert!((curvature(&left) - curvature(&right)).abs() <= 1e-12);
}

#[test]
fn a_bad_distance_or_an_unclamped_end_is_refused() {
    let c = cubic();
    for d in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            c.extended_by_distance(CurveEnd::End, Continuity::Natural, d),
            Err(SolidError::InvalidInput { .. })
        ));
    }
    let open = curve(
        &[
            [0.0, 0.0, 0.0, 1.0],
            [1.0, 1.0, 0.0, 1.0],
            [2.0, 0.0, 0.0, 1.0],
        ],
        2,
        &[0.0, 0.0, 0.0, 1.0, 2.0, 3.0],
        false,
    );
    assert!(matches!(
        open.extended_by_distance(CurveEnd::End, Continuity::Natural, 1.0),
        Err(SolidError::InvalidInput { .. })
    ));
    assert!(BSplineCurve::new(3)
        .extended_by_distance(CurveEnd::End, Continuity::Natural, 1.0)
        .is_err());
}

#[test]
fn a_short_segment_runs_on_far_past_its_length_as_its_line() {
    // a cutting line of a planar-path sweep: 0.2828 long, extended by 50 at either end
    let seg = curve(
        &[[0.0, 5.0, 0.0, 1.0], [-0.2, 5.2, 0.0, 1.0]],
        1,
        &[0.0, 0.0, 1.0, 1.0],
        false,
    );
    let dir = Vec3::new(-1.0, 1.0, 0.0) * FRAC_1_SQRT_2;
    for how in [
        Continuity::Natural,
        Continuity::Tangent,
        Continuity::Curvature,
    ] {
        let e = seg.extended_by_distance(CurveEnd::End, how, 50.0).unwrap();
        let (_, hi) = bounds(&e.curve);
        let want = Vec3::new(-0.2, 5.2, 0.0) + dir * 50.0;
        assert!((at(&e.curve, hi) - want).length() <= 1e-12, "{how:?}");
        kept(&seg, &e.curve, true);
        let s = seg
            .extended_by_distance(CurveEnd::Start, how, 50.0)
            .unwrap();
        let (lo, _) = bounds(&s.curve);
        let want = Vec3::new(0.0, 5.0, 0.0) - dir * 50.0;
        assert!((at(&s.curve, lo) - want).length() <= 1e-12, "{how:?}");
    }
}
