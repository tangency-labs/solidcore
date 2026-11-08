//! `sweep::stations` (#3139): a sweep's stations from its laws, against closed forms.
//!
//! - a helix's rotation-minimizing frame after one turn is its start frame turned about the
//!   tangent by `-2 pi b / sqrt(a^2 + b^2)` (radius `a`, rise `b` per radian);
//! - a planar quarter arc turns the frame a quarter turn about its axis, and keeps a direction
//!   in the plane's normal fixed (the direction orientation agrees);
//! - a straight path with twist `T` and end scale `k` turns the profile by `T u` and scales it by
//!   `1 + (k - 1) u` at the fraction `u`;
//! - a right-angle corner: the frame after it is the quarter turn about `t_0 x t_1`, the mitre
//!   at it the eighth turn stretched by `sqrt 2` along `t_1 - t_0`;
//! - a straight rail diverging at slope `s` scales the profile by `(d_0 + s x) / d_0`; an arc
//!   rail of radius `R` about the path's start by `sqrt(R^2 - x^2) / R`;
//! - a rail that stops short is refused, naming it.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers

use solidcore::curve::curve::Curve;
use solidcore::math::{Point3, Vec3};
use solidcore::primitive::api_curves::{create_arc, create_line_segment};
use solidcore::sweep::stations::{stations, Laws, Orientation, Path, Piece, RailFit, Station};
use solidcore::Kernel;
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI, SQRT_2, TAU};

fn apply(s: &Station, v: Vec3) -> Vec3 {
    s.m[0] * v.x + s.m[1] * v.y + s.m[2] * v.z
}

fn near(a: Vec3, b: Vec3, tol: f64) -> bool {
    a.distance(b) <= tol
}

/// The largest of `| |c_i|^2 - 1 |` and `|c_i . c_j|` over a station's columns.
fn defect(s: &Station) -> f64 {
    let [a, b, c] = s.m;
    [
        (a.dot(a) - 1.0).abs(),
        (b.dot(b) - 1.0).abs(),
        (c.dot(c) - 1.0).abs(),
        a.dot(b).abs(),
        b.dot(c).abs(),
        a.dot(c).abs(),
    ]
    .into_iter()
    .fold(0.0, f64::max)
}

/// The signed angle about the unit `axis` from `a` to `b` (both square to it).
fn angle_about(axis: Vec3, a: Vec3, b: Vec3) -> f64 {
    a.cross(b).dot(axis).atan2(a.dot(b))
}

#[test]
fn a_helix_frame_turns_by_its_torsion_against_the_screw() {
    let k = Kernel::default();
    let (a, b) = (5.0, 2.0);
    let path = Path {
        pieces: vec![Piece::Helix {
            centre: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::Z,
            start: Point3::new(a, 0.0, 0.0),
            angle: TAU,
            rise: b,
        }],
    };
    let t0 = Vec3::new(0.0, a, b) * (1.0 / a.hypot(b));
    let st = stations(&k, &path, t0, &Laws::default(), &[0.0, 0.25, 1.0]).unwrap();
    for s in &st {
        assert!(defect(s) < 1e-13, "{}", defect(s));
    }
    let end = st[2];
    // after a whole turn the tangent is the start's; the frame has turned about it
    assert!(near(apply(&end, t0), t0, 1e-12), "{:?}", apply(&end, t0));
    let normal = Vec3::new(-1.0, 0.0, 0.0);
    let want = -TAU * b / a.hypot(b);
    let got = angle_about(t0, normal, apply(&end, normal));
    // the angle modulo a turn
    let diff = (got - want).rem_euclid(TAU);
    assert!(diff.min(TAU - diff) < 1e-12, "{got} vs {want}");
    assert!(near(end.at, Point3::new(a, 0.0, TAU * b), 1e-12));
}

#[test]
fn a_planar_quarter_arc_turns_the_frame_a_quarter_turn() {
    let k = Kernel::default();
    let r = 7.0;
    let path = Path {
        pieces: vec![Piece::Arc {
            centre: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::Z,
            start: Point3::new(r, 0.0, 0.0),
            angle: FRAC_PI_2,
        }],
    };
    let at = [0.0, 0.5, 1.0];
    for orientation in [Orientation::Perpendicular, Orientation::Direction(Vec3::Z)] {
        let laws = Laws {
            orientation,
            ..Laws::default()
        };
        let st = stations(&k, &path, Vec3::Y, &laws, &at).unwrap();
        for (s, u) in st.iter().zip(at) {
            let phi = FRAC_PI_2 * u;
            assert!(near(
                s.at,
                Point3::new(r * phi.cos(), r * phi.sin(), 0.0),
                1e-12
            ));
            assert!(near(apply(s, Vec3::Z), Vec3::Z, 1e-12), "{orientation:?}");
            let want = Vec3::new(phi.cos(), phi.sin(), 0.0);
            assert!(near(apply(s, Vec3::X), want, 1e-12), "{orientation:?} {u}");
            assert!(defect(s) < 1e-13);
        }
    }
}

#[test]
fn a_straight_path_twists_and_scales_linearly() {
    let k = Kernel::default();
    let len = 12.0;
    let path = Path {
        pieces: vec![Piece::Line {
            a: Point3::new(0.0, 0.0, 0.0),
            b: Point3::new(0.0, 0.0, len),
        }],
    };
    let (twist, scale) = (1.5 * PI, 3.0);
    let laws = Laws {
        twist,
        end_scale: scale,
        ..Laws::default()
    };
    let at = [0.0, 0.3, 0.5, 1.0];
    let st = stations(&k, &path, Vec3::Z, &laws, &at).unwrap();
    for (s, u) in st.iter().zip(at) {
        let k_u = 1.0 + (scale - 1.0) * u;
        let turn = twist * u;
        let want = Vec3::new(turn.cos(), turn.sin(), 0.0) * k_u;
        assert!(near(apply(s, Vec3::X), want, 1e-12), "{u}");
        assert!(near(apply(s, Vec3::Z), Vec3::Z, 1e-12), "{u}");
        assert!(near(s.at, Point3::new(0.0, 0.0, len * u), 1e-12));
    }
}

#[test]
fn a_right_angle_corner_turns_the_frame_and_mitres_its_station() {
    let k = Kernel::default();
    let path = Path {
        pieces: vec![
            Piece::Line {
                a: Point3::new(0.0, 0.0, 0.0),
                b: Point3::new(0.0, 0.0, 10.0),
            },
            Piece::Line {
                a: Point3::new(0.0, 0.0, 10.0),
                b: Point3::new(10.0, 0.0, 10.0),
            },
        ],
    };
    let joints = path.corner_fractions().unwrap();
    assert_eq!(joints.len(), 1);
    assert!((joints[0] - 0.5).abs() < 1e-15);
    assert_eq!(path.joints().unwrap(), joints);
    let st = stations(&k, &path, Vec3::Z, &Laws::default(), &[joints[0], 1.0]).unwrap();
    let (mitre, end) = (st[0], st[1]);
    assert!(near(apply(&end, Vec3::Z), Vec3::X, 1e-12));
    assert!(near(apply(&end, Vec3::X), Vec3::Z * -1.0, 1e-12));
    // at the corner: the eighth turn about y, and x (along t_1 - t_0 there) stretched by sqrt 2
    let half = Vec3::new(FRAC_PI_4.sin(), 0.0, FRAC_PI_4.cos());
    assert!(
        near(apply(&mitre, Vec3::Z), half, 1e-12),
        "{:?}",
        apply(&mitre, Vec3::Z)
    );
    // x is the eighth turn of x, (x - z) / sqrt 2, along t_1 - t_0, stretched by sqrt 2
    let stretched = (Vec3::X - Vec3::Z) * (SQRT_2 / SQRT_2);
    assert!(
        near(apply(&mitre, Vec3::X), stretched, 1e-12),
        "{:?}",
        apply(&mitre, Vec3::X)
    );
    assert!(near(apply(&mitre, Vec3::Y), Vec3::Y, 1e-12));
}

/// A straight path along x, the profile in the (y, z) plane.
fn along_x(len: f64) -> Path {
    Path {
        pieces: vec![Piece::Line {
            a: Point3::new(0.0, 0.0, 0.0),
            b: Point3::new(len, 0.0, 0.0),
        }],
    }
}

#[test]
fn a_diverging_straight_rail_scales_the_profile() {
    let k = Kernel::default();
    let (len, d0, slope) = (20.0, 4.0, 0.25);
    let rail = Curve::from(
        create_line_segment(
            Point3::new(0.0, d0, 0.0),
            Point3::new(len, d0 + slope * len, 0.0),
        )
        .unwrap(),
    );
    let at = [0.0, 0.5, 1.0];
    for fit in [RailFit::Scale, RailFit::Stretch] {
        let laws = Laws {
            rails: vec![&rail],
            rail_fit: fit,
            ..Laws::default()
        };
        let st = stations(&k, &along_x(len), Vec3::X, &laws, &at).unwrap();
        for (s, u) in st.iter().zip(at) {
            let x = len * u;
            let k_u = (d0 + slope * x) / d0;
            assert!(near(apply(s, Vec3::Y), Vec3::Y * k_u, 1e-12), "{fit:?} {u}");
            let z = if fit == RailFit::Scale { k_u } else { 1.0 };
            assert!(near(apply(s, Vec3::Z), Vec3::Z * z, 1e-12), "{fit:?} {u}");
            assert!(near(s.at, Point3::new(x, 0.0, 0.0), 1e-12));
        }
    }
}

#[test]
fn an_arc_rail_scales_by_its_distance_across() {
    let k = Kernel::default();
    let (len, r): (f64, f64) = (6.0, 10.0);
    let start = (len / r).acos().to_degrees();
    let rail = Curve::from(create_arc(Point3::new(0.0, 0.0, 0.0), r, start, 90.0).unwrap());
    let laws = Laws {
        rails: vec![&rail],
        ..Laws::default()
    };
    let at = [0.0, 0.25, 0.5, 1.0];
    let st = stations(&k, &along_x(len), Vec3::X, &laws, &at).unwrap();
    for (s, u) in st.iter().zip(at) {
        let x = len * u;
        let k_u = (r * r - x * x).sqrt() / r;
        assert!(
            near(apply(s, Vec3::Y), Vec3::Y * k_u, 1e-12),
            "{u}: {:?}",
            apply(s, Vec3::Y)
        );
    }
}

#[test]
fn a_rail_that_stops_short_is_refused() {
    let k = Kernel::default();
    let rail = Curve::from(
        create_line_segment(Point3::new(0.0, 3.0, 0.0), Point3::new(5.0, 3.0, 0.0)).unwrap(),
    );
    let laws = Laws {
        rails: vec![&rail],
        ..Laws::default()
    };
    let err = stations(&k, &along_x(20.0), Vec3::X, &laws, &[0.0, 1.0]).unwrap_err();
    assert!(format!("{err:?}").contains("does not reach"), "{err:?}");
}
