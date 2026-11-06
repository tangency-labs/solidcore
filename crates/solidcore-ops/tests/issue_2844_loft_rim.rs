//! #2844: the round and the chamfer of a ruled loft's circular rim (solidcraft#12034, the
//! kernel case `fillets/12034-round-loft-square-to-circle-circular-rim`). The loft, degree 1,
//! runs from the square of side 30 at z = 0 (four lines) to the circle r 10 at z = 30 (four
//! quarter arcs, each centred on a side's direction): six faces, the top a plane disc. Each wall
//! is the ruled patch from a side to its arc; two walls meet along a straight side edge from a
//! square's corner to an arcs' joint, tangent to each other only at that joint.
//!
//! At the middle of a side (y = 0 for the side x = 15) the loft is symmetric, so the round's
//! section there is the ball's circle in that plane, between the disc and the wall's middle
//! ruling from (15, 0, 0) to (10, 0, 30). They meet at 90 degrees + atan(1/6) inside the
//! material, so the ball of radius r stands k r inside the rim, k = cot of half that angle =
//! (sqrt 37 - 1) / 6, its centre at (10 - k r, 0, 30 - r). A chamfer's setbacks are its
//! distance at its edge's middle (its offsets are taken there), so its chord there joins
//! (10 - d, 0, 30) to the point d down the ruling.
//!
//! Its own test binary, on the kernel's default fixes: none is set.
#![allow(clippy::unwrap_used, clippy::indexing_slicing, clippy::panic)] // test helpers

mod common;

use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

use common::{edge_at, rel, volume};
use solidcore_ops::sweep::{loft, Curve};
use solidcore_ops::{blend, Blend, Body, CheckLevel, Corners, PointClass};

/// The sides' and arcs' mixed area: the loft's section at height 30 s is (1 - s) L(t) + s A(t),
/// L a side and A its arc at the same parameter (the record's rows), so its area is
/// (1 - s)^2 900 + 2 s (1 - s) M + s^2 100 pi. M by Gauss-Legendre on the rows, to 1e-12.
const MIXED: f64 = 537.671338762224;

/// The loft's volume, 30 times the mean of that area: 10 (900 + M + 100 pi).
fn loft_volume() -> f64 {
    10.0 * (900.0 + MIXED + 100.0 * PI)
}

/// How much of the loft lies in the band rho >= 10 - s, z >= 30 - s: every wall leans in, so
/// a round of radius s or a chamfer of distance s takes nothing outside it.
fn band(s: f64) -> f64 {
    let s0 = 1.0 - s / 30.0;
    let square = 900.0 * (1.0 - s0).powi(3) / 3.0;
    let mixed = 2.0 * MIXED * (1.0 / 6.0 - s0 * s0 / 2.0 + s0.powi(3) / 3.0);
    let circle = 100.0 * PI * (1.0 - s0.powi(3)) / 3.0;
    30.0 * (square + mixed + circle) - s * PI * (10.0 - s).powi(2)
}

/// The loft: the sides paired with the arcs centred on their directions.
fn loft_body() -> Body {
    let side = |a: [f64; 2], b: [f64; 2]| Curve::Line {
        start: [a[0], a[1], 0.0],
        end: [b[0], b[1], 0.0],
    };
    let square = vec![
        side([15.0, -15.0], [15.0, 15.0]),
        side([15.0, 15.0], [-15.0, 15.0]),
        side([-15.0, 15.0], [-15.0, -15.0]),
        side([-15.0, -15.0], [15.0, -15.0]),
    ];
    let circle = (0..4)
        .map(|q| Curve::Arc {
            center: [0.0, 0.0, 30.0],
            axis: [0.0, 0.0, 1.0],
            reference: [1.0, 0.0, 0.0],
            radius: 10.0,
            start: f64::from(q) * FRAC_PI_2 - FRAC_PI_4,
            end: f64::from(q) * FRAC_PI_2 + FRAC_PI_4,
        })
        .collect();
    let b = loft(&[square, circle], 1, true).unwrap();
    assert!(b.check(CheckLevel::All).unwrap().is_valid());
    assert!(rel(volume(&b), loft_volume()) < 1e-9, "{}", volume(&b));
    b
}

/// `p` turned `q` quarters about the z axis.
fn turned(p: [f64; 3], q: i32) -> [f64; 3] {
    let (s, c) = (f64::from(q) * FRAC_PI_2).sin_cos();
    [c * p[0] - s * p[1], s * p[0] + c * p[1], p[2]]
}

/// The rim's four arcs, each found at its middle.
fn rim(b: &Body) -> Vec<usize> {
    (0..4)
        .map(|q| edge_at(b, turned([10.0, 0.0, 30.0], q)))
        .collect()
}

/// k: the disc and a wall's middle ruling inside the material, cot of half their angle.
fn inset() -> f64 {
    (37.0_f64.sqrt() - 1.0) / 6.0
}

/// `p`, turned by each quarter about the z axis (to each side's middle or each joint), lies
/// where `class` says.
fn at_each(b: &Body, p: [f64; 3], class: PointClass, what: &str) {
    for q in 0..4 {
        let got = b.classify_point(turned(p, q)).unwrap();
        assert_eq!(got, class, "{what}, quarter {q}");
    }
}

/// The rim of a fresh loft blended by `kind` of size `size`: one closed valid solid with less
/// than the loft and more than the loft less the band; it and its volume.
fn blended(kind: Blend, size: f64) -> (Body, f64) {
    let b = loft_body();
    let edges = rim(&b);
    if let Err(e) = blend(&b, kind, &edges, &[size], &[], Corners::Filled) {
        panic!("{kind:?} {size}: {e}");
    }
    assert!(
        b.check(CheckLevel::All).unwrap().is_valid(),
        "{kind:?} {size}"
    );
    assert_eq!(b.solid_count().unwrap(), 1, "{kind:?} {size}");
    let v = volume(&b);
    let loft = loft_volume();
    assert!(v < loft && loft - v < band(size), "{kind:?} {size}: {v}");
    (b, v)
}

/// What every rounding of the rim leaves at the old joints, the arcs' meeting points: their
/// corner taken, the side edges below the round untouched.
fn joints(b: &Body) {
    let (s, c) = FRAC_PI_4.sin_cos();
    // 0.01 inside the rim and below the disc, in the plane of a side edge: material before
    at_each(
        b,
        [9.99 * c, 9.99 * s, 29.99],
        PointClass::Outside,
        "corner at a joint",
    );
    // on the side edge from (15, 15, 0) to (10 cos 45, 10 sin 45, 30), at z = 25
    let x = 15.0 - 25.0 / 30.0 * (15.0 - 10.0 * c);
    at_each(
        b,
        [x, x, 25.0],
        PointClass::OnBoundary,
        "side edge below the round",
    );
}

#[test]
fn the_round_of_the_rim_is_one_solid_with_the_balls_section_at_each_sides_middle() {
    let k = inset();
    let n = (k * k + 1.0).sqrt();
    let s37 = 37.0_f64.sqrt();
    let mut volumes = Vec::new();
    for r in [0.5, 1.0, 2.0] {
        let (b, v) = blended(Blend::Round, r);
        let c = [10.0 - k * r, 0.0, 30.0 - r];
        // along the ball's radius towards the rim's corner (10, 0, 30)
        let along = |t: f64| [c[0] + t * r * k / n, 0.0, c[2] + t * r / n];
        let m = along(1.0);
        at_each(&b, m, PointClass::OnBoundary, "the section's middle");
        at_each(
            &b,
            [c[0], 0.0, 30.0],
            PointClass::OnBoundary,
            "the disc's rail",
        );
        let wall = [c[0] + 6.0 * r / s37, 0.0, c[2] + r / s37];
        at_each(&b, wall, PointClass::OnBoundary, "the wall's rail");
        at_each(&b, along(0.5), PointClass::Inside, "inside the ball");
        let corner = [(m[0] + 10.0) / 2.0, 0.0, (m[2] + 30.0) / 2.0];
        at_each(
            &b,
            corner,
            PointClass::Outside,
            "between the section and the rim",
        );
        joints(&b);
        volumes.push(v);
    }
    // a larger ball takes more
    assert!(
        volumes[0] > volumes[1] && volumes[1] > volumes[2],
        "{volumes:?}"
    );
}

#[test]
fn the_chamfer_of_the_rim_is_one_solid_with_its_chord_at_each_sides_middle() {
    let d = 1.0;
    let (b, _) = blended(Blend::Chamfer, d);
    let s37 = 37.0_f64.sqrt();
    let disc = [10.0 - d, 0.0, 30.0];
    let wall = [10.0 + d / s37, 0.0, 30.0 - 6.0 * d / s37];
    let mid = [(disc[0] + wall[0]) / 2.0, 0.0, (disc[2] + wall[2]) / 2.0];
    at_each(&b, disc, PointClass::OnBoundary, "the disc's setback");
    at_each(&b, wall, PointClass::OnBoundary, "the wall's setback");
    at_each(&b, mid, PointClass::OnBoundary, "the chord's middle");
    // across the chord from the rim's corner: material; towards it: taken
    let inside = [1.5 * mid[0] - 5.0, 0.0, 1.5 * mid[2] - 15.0];
    at_each(&b, inside, PointClass::Inside, "beyond the chord");
    let corner = [(mid[0] + 10.0) / 2.0, 0.0, (mid[2] + 30.0) / 2.0];
    at_each(
        &b,
        corner,
        PointClass::Outside,
        "between the chord and the rim",
    );
    joints(&b);
}
