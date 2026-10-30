//! The sweeps: extrusions (straight and drafted), revolutions, sweeps along lines, arcs and
//! a closed path, lofts, lofts to a point, sweeps through stations, planar faces and the
//! argument checks (the screw sweeps' cases come with their wrappers).
//!
//! Source: tangency-labs/smlib `capi/tests/sweeps.cpp` (Apache-2.0, Copyright SolidCraft
//! contributors), replay/1.0.1: `List`, `circle`, `expect` and the cases named in each test.
//! Its cases that read the flat arrays (`smlib_brep_export`), project curves or measure a
//! screw's walls against their true surfaces need later families and are left for those.

#![allow(clippy::unwrap_used, clippy::panic)] // test helpers

use std::f64::consts::PI;

use solidcore_smlib_api::{
    Brep, Smlib, Status, CURVE_ARC, CURVE_LINE, DRAFT_CORNER_CURVATURE, DRAFT_CORNER_TANGENT,
    ENDS_CAPPED, ENDS_JOINED, SWEEP_FOLLOW,
};

const Z: [f64; 3] = [0.0, 0.0, 1.0];
const O: [f64; 3] = [0.0, 0.0, 0.0];

/// A curve list being built (`List`).
#[derive(Clone, Default)]
struct List {
    kinds: Vec<i64>,
    reals: Vec<f64>,
}

impl List {
    fn line(mut self, a: [f64; 3], b: [f64; 3]) -> Self {
        self.kinds.push(CURVE_LINE);
        self.reals.extend(a);
        self.reals.extend(b);
        self
    }

    fn arc(mut self, c: [f64; 3], n: [f64; 3], r: [f64; 3], radius: f64, a0: f64, a1: f64) -> Self {
        self.kinds.push(CURVE_ARC);
        self.reals.extend(c);
        self.reals.extend(n);
        self.reals.extend(r);
        self.reals.extend([radius, a0, a1]);
        self
    }

    /// A closed polygon through `pts`.
    fn polygon(mut self, pts: &[[f64; 3]]) -> Self {
        for (i, &p) in pts.iter().enumerate() {
            self = self.line(p, pts[(i + 1) % pts.len()]);
        }
        self
    }

    fn append(mut self, o: &Self) -> Self {
        self.kinds.extend(&o.kinds);
        self.reals.extend(&o.reals);
        self
    }
}

/// A full circle of radius `r` about `n` through `c`, starting along `x`.
fn circle(c: [f64; 3], n: [f64; 3], x: [f64; 3], r: f64) -> List {
    List::default().arc(c, n, x, r, 0.0, 2.0 * PI)
}

/// `expect`: the call gave one valid manifold solid of volume `v` within `rel`.
fn expect(api: &mut Smlib, name: &str, r: Result<Brep, Status>, v: f64, rel: f64) {
    let b = match r {
        Ok(b) => b,
        Err(s) => panic!("{name}: status {s}"),
    };
    let got = api.brep_volume(&b, 1e-6).unwrap();
    assert!(
        api.brep_is_manifold_solid(&b).unwrap(),
        "{name}: not a manifold solid"
    );
    let (report, bad) = api.brep_check(&b, 2).unwrap();
    let err = (got - v).abs() / v;
    eprintln!("{name:34} volume {got:.6} expected {v:.6} (rel {err:.2e})");
    assert!(err <= rel, "{name}: volume {got}, expected {v}");
    assert_eq!(bad, 0, "{name}: check\n{report}");
}

fn rect() -> List {
    List::default().polygon(&[
        [0.0, 0.0, 0.0],
        [10.0, 0.0, 0.0],
        [10.0, 6.0, 0.0],
        [0.0, 6.0, 0.0],
    ])
}

fn square(h: f64) -> List {
    List::default().polygon(&[
        [-5.0, -5.0, h],
        [5.0, -5.0, h],
        [5.0, 5.0, h],
        [-5.0, 5.0, h],
    ])
}

/// Extrusions: a 10 x 6 rectangle 5 up and down, a disc with a hole, a plate with a tiny
/// hole, a square drafted 10 degrees in and out, up and down (a frustum of a pyramid).
#[test]
fn extrusions() {
    let mut api = Smlib::default();
    let r = rect();
    let res = api.extrude(&r.kinds, &r.reals, Z, 5.0, 0.0, DRAFT_CORNER_TANGENT, 1);
    expect(&mut api, "extrude rectangle", res, 300.0, 1e-9);
    let down = [0.0, 0.0, -1.0];
    let res = api.extrude(&r.kinds, &r.reals, down, 5.0, 0.0, DRAFT_CORNER_TANGENT, 1);
    expect(&mut api, "extrude rectangle down", res, 300.0, 1e-9);

    let ring = circle(O, Z, [1.0, 0.0, 0.0], 10.0).append(&circle(
        O,
        [0.0, 0.0, -1.0],
        [1.0, 0.0, 0.0],
        4.0,
    ));
    let res = api.extrude(
        &ring.kinds,
        &ring.reals,
        Z,
        3.0,
        0.0,
        DRAFT_CORNER_TANGENT,
        1,
    );
    expect(
        &mut api,
        "extrude annulus",
        res,
        PI * (100.0 - 16.0) * 3.0,
        1e-6,
    );

    let plate = List::default()
        .polygon(&[
            [0.0, 0.0, 0.0],
            [3000.0, 0.0, 0.0],
            [3000.0, 1000.0, 0.0],
            [0.0, 1000.0, 0.0],
        ])
        .append(&circle(
            [30.0, 20.0, 0.0],
            [0.0, 0.0, -1.0],
            [1.0, 0.0, 0.0],
            0.002,
        ));
    let res = api.extrude(
        &plate.kinds,
        &plate.reals,
        Z,
        2.0,
        0.0,
        DRAFT_CORNER_TANGENT,
        1,
    );
    let v = 6.0e6 - PI * 0.002 * 0.002 * 2.0;
    expect(&mut api, "extrude plate with a tiny hole", res, v, 1e-12);

    let sq = square(0.0);
    let t = (10.0 * PI / 180.0).tan();
    let a2 = 10.0 - 2.0 * 4.0 * t;
    let a3 = 10.0 + 2.0 * 4.0 * t;
    let frustum = |a: f64| 4.0 / 3.0 * (100.0 + a * a + 10.0 * a);
    for (name, dir, draft, a) in [
        ("extrude square drafted 10 deg", Z, 10.0, a2),
        ("extrude square drafted -10 deg", Z, -10.0, a3),
        ("extrude square drafted down 10 deg", down, 10.0, a2),
        ("extrude square drafted down -10 deg", down, -10.0, a3),
    ] {
        let d = draft * PI / 180.0;
        let res = api.extrude(&sq.kinds, &sq.reals, dir, 4.0, d, DRAFT_CORNER_TANGENT, 1);
        expect(&mut api, name, res, frustum(a), 1e-6);
    }
}

/// A slot 0.2 wide cut into a 10 x 10 square, drafted 10 degrees outward over 2: the slot's
/// sides close on each other (`drafted slot closes`), with the curvature and tangent rules.
#[test]
fn drafted_slot_closes() {
    let mut api = Smlib::default();
    let slot = List::default().polygon(&[
        [0.0, 0.0, 0.0],
        [10.0, 0.0, 0.0],
        [10.0, 10.0, 0.0],
        [5.1, 10.0, 0.0],
        [5.1, 2.0, 0.0],
        [4.9, 2.0, 0.0],
        [4.9, 10.0, 0.0],
        [0.0, 10.0, 0.0],
    ]);
    let tt = (10.0 * PI / 180.0).tan();
    let v = ((10.0 + 4.0 * tt).powi(3) - 1000.0) / (6.0 * tt) - 0.08 / tt;
    assert!((v - 213.984095421783).abs() < 1e-9);
    for rule in [DRAFT_CORNER_CURVATURE, DRAFT_CORNER_TANGENT] {
        let draft = -10.0 * PI / 180.0;
        let r = api.extrude(&slot.kinds, &slot.reals, Z, 2.0, draft, rule, 1);
        let b = r.clone().unwrap();
        assert_eq!(api.brep_counts(&b).unwrap()[0], 10);
        expect(&mut api, "drafted slot closes", r, v, 1e-9);
    }
}

/// Revolutions about z of a profile in the xz plane: a square off the axis (a ring), a
/// quarter turn of it, a half disc on the axis (a sphere), partial turns, a holed square.
#[test]
fn revolutions() {
    let mut api = Smlib::default();
    let off = List::default().polygon(&[
        [6.0, 0.0, 0.0],
        [8.0, 0.0, 0.0],
        [8.0, 0.0, 2.0],
        [6.0, 0.0, 2.0],
    ]);
    let res = api.revolve(&off.kinds, &off.reals, O, Z, 2.0 * PI, 1);
    expect(
        &mut api,
        "revolve square full",
        res,
        2.0 * PI * 7.0 * 4.0,
        1e-5,
    );
    let res = api.revolve(&off.kinds, &off.reals, O, Z, PI / 2.0, 1);
    expect(
        &mut api,
        "revolve square 90",
        res,
        PI / 2.0 * 7.0 * 4.0,
        1e-6,
    );
    let half = List::default().line([0.0, 0.0, -5.0], [0.0, 0.0, 5.0]).arc(
        O,
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
        5.0,
        0.0,
        PI,
    );
    let res = api.revolve(&half.kinds, &half.reals, O, Z, 2.0 * PI, 1);
    expect(
        &mut api,
        "revolve half disc",
        res,
        4.0 / 3.0 * PI * 125.0,
        1e-5,
    );
    for turn in [0.5, 0.75, 0.9] {
        let res = api.revolve(&half.kinds, &half.reals, O, Z, 2.0 * PI * turn, 1);
        expect(
            &mut api,
            "revolve half disc turn",
            res,
            turn * 4.0 / 3.0 * PI * 125.0,
            1e-5,
        );
        let res = api.revolve(&off.kinds, &off.reals, O, Z, 2.0 * PI * turn, 1);
        expect(
            &mut api,
            "revolve square turn",
            res,
            2.0 * PI * turn * 7.0 * 4.0,
            1e-5,
        );
    }
    let holed = List::default()
        .polygon(&[
            [10.0, 0.0, 0.0],
            [20.0, 0.0, 0.0],
            [20.0, 0.0, 5.0],
            [10.0, 0.0, 5.0],
        ])
        .polygon(&[
            [12.0, 0.0, 1.0],
            [12.0, 0.0, 2.0],
            [14.0, 0.0, 2.0],
            [14.0, 0.0, 1.0],
        ]);
    for turn in [0.25, 0.75, 1.0] {
        let res = api.revolve(&holed.kinds, &holed.reals, O, Z, 2.0 * PI * turn, 1);
        let v = turn * 2.0 * PI * (50.0 * 15.0 - 2.0 * 13.0);
        expect(&mut api, "revolve holed square", res, v, 1e-5);
    }
}

/// Sweeps: a disc along a line (a cylinder), along a quarter arc (a quarter torus), off the
/// path's start, and round a closed stadium of lines and arcs (one wall per piece, no caps).
#[test]
fn sweeps_along_paths() {
    let mut api = Smlib::default();
    let disc = circle(O, Z, [1.0, 0.0, 0.0], 2.0);
    let up = List::default().line(O, [0.0, 0.0, 10.0]);
    let res = api.sweep(
        &disc.kinds,
        &disc.reals,
        &up.kinds,
        &up.reals,
        SWEEP_FOLLOW,
        1,
    );
    expect(
        &mut api,
        "sweep disc along line",
        res,
        PI * 4.0 * 10.0,
        1e-4,
    );

    let bend = List::default().arc(
        [-20.0, 0.0, 0.0],
        [0.0, -1.0, 0.0],
        [1.0, 0.0, 0.0],
        20.0,
        0.0,
        PI / 2.0,
    );
    let res = api.sweep(
        &disc.kinds,
        &disc.reals,
        &bend.kinds,
        &bend.reals,
        SWEEP_FOLLOW,
        1,
    );
    expect(
        &mut api,
        "sweep disc along quarter arc",
        res,
        PI * 4.0 * (PI / 2.0 * 20.0),
        1e-3,
    );
    let offset = circle([5.0, 0.0, 0.0], Z, [1.0, 0.0, 0.0], 2.0);
    let res = api.sweep(
        &offset.kinds,
        &offset.reals,
        &bend.kinds,
        &bend.reals,
        SWEEP_FOLLOW,
        1,
    );
    expect(
        &mut api,
        "sweep offset disc along arc",
        res,
        PI * 4.0 * (PI / 2.0 * 25.0),
        1e-3,
    );

    let stadium = List::default()
        .line([0.0, -10.0, 0.0], [20.0, -10.0, 0.0])
        .arc([20.0, 0.0, 0.0], Z, [0.0, -1.0, 0.0], 10.0, 0.0, PI)
        .line([20.0, 10.0, 0.0], [0.0, 10.0, 0.0])
        .arc(O, Z, [0.0, 1.0, 0.0], 10.0, 0.0, PI);
    let section = circle([0.0, -10.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 2.0);
    let res = api.sweep(
        &section.kinds,
        &section.reals,
        &stadium.kinds,
        &stadium.reals,
        SWEEP_FOLLOW,
        1,
    );
    let b = res.clone().unwrap();
    assert_eq!(api.brep_counts(&b).unwrap()[0], 4);
    expect(
        &mut api,
        "sweep disc round a stadium",
        res,
        PI * 4.0 * (40.0 + 20.0 * PI),
        1e-4,
    );
}

/// `smlib_sweep_placed` along a straight path: the motions are the translations to the
/// start, middle and end.
#[test]
fn sweep_placed_along_a_line() {
    let mut api = Smlib::default();
    let disc = circle(O, Z, [1.0, 0.0, 0.0], 2.0);
    let up = List::default().line(O, [0.0, 0.0, 10.0]);
    let (_, placed) = api
        .sweep_placed(
            &disc.kinds,
            &disc.reals,
            &up.kinds,
            &up.reals,
            SWEEP_FOLLOW,
            1,
        )
        .unwrap();
    for (k, z) in [0.0, 5.0, 10.0].into_iter().enumerate() {
        let m = &placed[12 * k..12 * k + 12];
        assert_eq!(
            m,
            &[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, z]
        );
    }
}

/// Lofts: square to square (a box), circle to circle (a frustum, ruled), ring to ring, frame
/// to frame.
#[test]
fn lofts() {
    let mut api = Smlib::default();
    let boxes = List::default()
        .polygon(&[
            [0.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            [10.0, 10.0, 0.0],
            [0.0, 10.0, 0.0],
        ])
        .polygon(&[
            [0.0, 0.0, 8.0],
            [10.0, 0.0, 8.0],
            [10.0, 10.0, 8.0],
            [0.0, 10.0, 8.0],
        ]);
    let res = api.loft(&boxes.kinds, &boxes.reals, &[4, 4], 1, 1);
    expect(&mut api, "loft square to square", res, 800.0, 1e-6);
    let cone = circle(O, Z, [1.0, 0.0, 0.0], 5.0).append(&circle(
        [0.0, 0.0, 6.0],
        Z,
        [1.0, 0.0, 0.0],
        3.0,
    ));
    let res = api.loft(&cone.kinds, &cone.reals, &[1, 1], 1, 1);
    expect(
        &mut api,
        "loft circle to circle",
        res,
        PI * 6.0 / 3.0 * (25.0 + 9.0 + 15.0),
        1e-4,
    );
    let x = [1.0, 0.0, 0.0];
    let minus_z = [0.0, 0.0, -1.0];
    let tube = circle(O, Z, x, 20.0)
        .append(&circle(O, minus_z, x, 12.0))
        .append(&circle([0.0, 0.0, 20.0], Z, x, 20.0))
        .append(&circle([0.0, 0.0, 20.0], minus_z, x, 12.0));
    let res = api.loft(&tube.kinds, &tube.reals, &[2, 2], 1, 1);
    expect(
        &mut api,
        "loft ring to ring",
        res,
        PI * (400.0 - 144.0) * 20.0,
        1e-4,
    );
    let frame = |h: f64| {
        List::default()
            .polygon(&[
                [-10.0, -10.0, h],
                [10.0, -10.0, h],
                [10.0, 10.0, h],
                [-10.0, 10.0, h],
            ])
            .polygon(&[
                [-5.0, -5.0, h],
                [-5.0, 5.0, h],
                [5.0, 5.0, h],
                [5.0, -5.0, h],
            ])
    };
    let frames = frame(0.0).append(&frame(10.0));
    let res = api.loft(&frames.kinds, &frames.reals, &[8, 8], 1, 1);
    expect(&mut api, "loft frame to frame", res, 300.0 * 10.0, 1e-6);
}

/// Lofts to a point: a square base to an apex above its centre, a triangle to an apex below;
/// an apex in the base's plane, or a curved base, is refused.
#[test]
fn lofts_to_a_point() {
    let mut api = Smlib::default();
    let base = square(0.0);
    let big = List::default().polygon(&[
        [-10.0, -10.0, 0.0],
        [10.0, -10.0, 0.0],
        [10.0, 10.0, 0.0],
        [-10.0, 10.0, 0.0],
    ]);
    let res = api.loft_to_point(&big.kinds, &big.reals, [0.0, 0.0, 15.0]);
    expect(
        &mut api,
        "loft square to a point",
        res,
        400.0 * 15.0 / 3.0,
        1e-9,
    );
    let tri = List::default().polygon(&[[0.0, 0.0, 5.0], [12.0, 0.0, 5.0], [0.0, 9.0, 5.0]]);
    let res = api.loft_to_point(&tri.kinds, &tri.reals, [30.0, -7.0, -1.0]);
    expect(
        &mut api,
        "loft triangle to a point below",
        res,
        54.0 * 6.0 / 3.0,
        1e-9,
    );
    let flat = api.loft_to_point(&base.kinds, &base.reals, [30.0, 30.0, 0.0]);
    assert_eq!(flat.err(), Some(Status::ERR_ARG));
    let disc = circle(O, Z, [1.0, 0.0, 0.0], 5.0);
    let curved = api.loft_to_point(&disc.kinds, &disc.reals, [0.0, 0.0, 15.0]);
    assert_eq!(curved.err(), Some(Status::ERR_ARG));
}

/// Sweeps through stations: a square twisted 90 degrees over 40, a circle grown to a cone
/// frustum, a square carried round a circle and joined; hostile maps.
#[test]
fn sweeps_through_stations() {
    let mut api = Smlib::default();
    let n = 65;
    let sq = square(0.0);
    let twist: Vec<[f64; 12]> = (0..n)
        .map(|i| {
            let f = f64::from(i) / f64::from(n - 1);
            let a = 0.5 * PI * f;
            [
                a.cos(),
                a.sin(),
                0.0,
                -a.sin(),
                a.cos(),
                0.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
                40.0 * f,
            ]
        })
        .collect();
    let res = api.sweep_stations(&sq.kinds, &sq.reals, &twist, ENDS_CAPPED);
    expect(&mut api, "stations: twisted square", res, 4000.0, 1e-5);

    let disc = circle(O, Z, [1.0, 0.0, 0.0], 5.0);
    let grow: Vec<[f64; 12]> = (0..n)
        .map(|i| {
            let f = f64::from(i) / f64::from(n - 1);
            let k = 1.0 + f;
            [k, 0.0, 0.0, 0.0, k, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 30.0 * f]
        })
        .collect();
    let res = api.sweep_stations(&disc.kinds, &disc.reals, &grow, ENDS_CAPPED);
    expect(
        &mut api,
        "stations: cone frustum",
        res,
        10.0 * PI * (25.0 + 100.0 + 50.0),
        1e-5,
    );

    let m = 97;
    let ring = List::default().polygon(&[
        [26.0, 0.0, -4.0],
        [34.0, 0.0, -4.0],
        [34.0, 0.0, 4.0],
        [26.0, 0.0, 4.0],
    ]);
    let round: Vec<[f64; 12]> = (0..m)
        .map(|i| {
            let f = f64::from(i) / f64::from(m - 1);
            let (phi, psi) = (2.0 * PI * f, 0.5 * PI * f);
            let (cp, sp, cs, ss) = (phi.cos(), phi.sin(), psi.cos(), psi.sin());
            let d = [30.0 - 30.0 * cs, 0.0, 30.0 * ss];
            [
                cp * cs,
                sp * cs,
                -ss,
                -sp,
                cp,
                0.0,
                cp * ss,
                sp * ss,
                cs,
                cp * d[0] - sp * d[1],
                sp * d[0] + cp * d[1],
                d[2],
            ]
        })
        .collect();
    let res = api.sweep_stations(&ring.kinds, &ring.reals, &round, ENDS_JOINED);
    expect(
        &mut api,
        "stations: twisted ring, joined",
        res,
        64.0 * 2.0 * PI * 30.0,
        1e-5,
    );

    // Hostile maps: a singular matrix, a NaN, one station, unknown ends.
    let mut flat = twist.clone();
    flat[5][0] = 0.0;
    flat[5][1] = 0.0;
    flat[5][2] = 0.0;
    let r = api.sweep_stations(&sq.kinds, &sq.reals, &flat, ENDS_CAPPED);
    assert_eq!(r.err(), Some(Status::ERR_ARG));
    let mut nan_map = twist.clone();
    nan_map[0][7] = f64::NAN;
    let r = api.sweep_stations(&sq.kinds, &sq.reals, &nan_map, ENDS_CAPPED);
    assert_eq!(r.err(), Some(Status::ERR_ARG));
    let r = api.sweep_stations(&sq.kinds, &sq.reals, &twist[..1], ENDS_CAPPED);
    assert_eq!(r.err(), Some(Status::ERR_ARG));
    let r = api.sweep_stations(&sq.kinds, &sq.reals, &twist, 7);
    assert_eq!(r.err(), Some(Status::ERR_ARG));
}

/// Planar faces from loops: a rectangle with a round hole is one face of two loops.
#[test]
fn planar_faces_from_loops() {
    let mut api = Smlib::default();
    let loops = rect().append(&circle([5.0, 3.0, 0.0], Z, [1.0, 0.0, 0.0], 1.0));
    let b = api.planar_faces(&loops.kinds, &loops.reals).unwrap();
    assert_eq!(api.brep_counts(&b).unwrap()[0], 1);
    let area = api.brep_area(&b, 1e-6).unwrap();
    assert!((area - (60.0 - PI)).abs() <= 1e-6 * 60.0, "area {area}");
}

/// Hostile lists: nothing is built, nothing panics.
#[test]
fn hostile_lists() {
    let mut api = Smlib::default();
    let six = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
    let t = DRAFT_CORNER_TANGENT;
    assert_eq!(
        api.extrude(&[9], &six, Z, 1.0, 0.0, t, 1).err(),
        Some(Status::ERR_ARG)
    );
    assert_eq!(
        api.extrude(&[CURVE_LINE], &six[..5], Z, 1.0, 0.0, t, 1)
            .err(),
        Some(Status::ERR_ARG)
    );
    let nan6 = [0.0, f64::NAN, 0.0, 1.0, 0.0, 0.0];
    assert_eq!(
        api.extrude(&[CURVE_LINE], &nan6, Z, 1.0, 0.0, t, 1).err(),
        Some(Status::ERR_ARG)
    );
    assert_eq!(
        api.extrude(&[3, 3, 1_000_000_000], &six, Z, 1.0, 0.0, t, 1)
            .err(),
        Some(Status::ERR_ARG)
    );
    let r = rect();
    assert_eq!(
        api.extrude(&r.kinds, &r.reals, Z, -1.0, 0.0, t, 1).err(),
        Some(Status::ERR_ARG)
    );
    let off = List::default().polygon(&[
        [6.0, 0.0, 0.0],
        [8.0, 0.0, 0.0],
        [8.0, 0.0, 2.0],
        [6.0, 0.0, 2.0],
    ]);
    assert_eq!(
        api.revolve(&off.kinds, &off.reals, O, Z, 7.0, 1).err(),
        Some(Status::ERR_ARG)
    );
    let boxes = square(0.0).append(&square(8.0));
    assert_eq!(
        api.loft(&boxes.kinds, &boxes.reals, &[4], 1, 1).err(),
        Some(Status::ERR_ARG)
    );
    assert_eq!(
        api.loft(&boxes.kinds, &boxes.reals, &[4, 3], 1, 1).err(),
        Some(Status::ERR_ARG)
    );
}
