//! Several concave rounds at one vertex, and mixed corners with every edge blended
//! (tangency-labs/solidcore#2361, from smlib#84), under `fillet_blend_corner`: the Setback,
//! its fill certified patch by patch to the owner's one rule `II >= k G` (k = 0 beside convex
//! rounds and faces, `-1 / r` beside a concave round of radius r), built through the fillet
//! executive. Checked against closed forms and the B-rep check, not against C++.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers

use solidcore::blend::fillet_executive::session::FilletSession;
use solidcore::blend::fillet_solver::p1::FilletSurfaceGeneratorType;
use solidcore::boolean::api_brep::boolean_difference;
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::geom::global_solver::SolverOperation;
use solidcore::kernel::UpstreamFixes;
use solidcore::mass::integrator::properties::{brep_volume, face_area, PropertyAttributes};
use solidcore::math::{Point3, BIG_DOUBLE};
use solidcore::primitive::api_primitives::{create_box, create_pyramid};
use solidcore::query::api_queries::{get_edges, get_faces};
use solidcore::topology::brep::p5::{stitch_and_orient, StitchAndOrientOptions};
use solidcore::topology::brep::validation::{check_brep, CheckLevel};
use solidcore::topology::{edge, EdgeId, FaceId};
use solidcore::{Brep, Kernel};

/// Every fix but `geometric_g2_continuity`, `fillet_blend_corner` included.
fn kernel() -> Kernel {
    let mut k = Kernel::default();
    let mut f = UpstreamFixes::default();
    for &name in UpstreamFixes::NAMES {
        if name != "geometric_g2_continuity" {
            *f.flag_mut(name).unwrap() = true;
        }
    }
    k.fixes = f;
    k
}

/// The edge of `b` through `at`.
fn pick(k: &Kernel, b: &Brep, at: [f64; 3]) -> EdgeId {
    let p = Point3::new(at[0], at[1], at[2]);
    let mut best: Option<(EdgeId, f64)> = None;
    for e in get_edges(b).unwrap() {
        let Some(c) = edge::curve(b, e).unwrap() else {
            continue;
        };
        let ivl = edge::interval(b, e).unwrap();
        let d = c
            .drop_point(
                k,
                &ivl,
                p,
                None,
                BIG_DOUBLE,
                None,
                SolverOperation::Minimize,
            )
            .unwrap();
        if d.success && d.distance < 0.1 && best.is_none_or(|(_, bd)| d.distance < bd) {
            best = Some((e, d.distance));
        }
    }
    best.unwrap().0
}

/// The edges through `at` rounded at `radii` (circular, G1), the corners the executive's
/// own; the body stitched and checked (Level 2, no complaint).
fn round(k: &mut Kernel, mut b: Brep, at: &[[f64; 3]], radii: &[f64]) -> Brep {
    let edges: Vec<EdgeId> = at.iter().map(|&p| pick(k, &b, p)).collect();
    let xs = vec![FilletSurfaceGeneratorType::Circular; edges.len()];
    let (ok, comment) = {
        let mut s = FilletSession::new(k, &mut b).unwrap();
        s.prepare_per_edge_constant_radius(k, &edges, radii, &xs, 1, 1.0)
            .unwrap();
        let mut ok = s.create_corners(k, None).is_ok();
        if ok {
            ok = s.do_filleting(k, None).is_ok();
        }
        (ok, s.error_info().comment.clone())
    };
    assert!(ok, "the rounds were refused: {comment}");
    stitch_and_orient(k, &mut b, &StitchAndOrientOptions::default()).unwrap();
    let lines = check_brep(k, &mut b, CheckLevel::Level2).unwrap().lines();
    assert!(lines.is_empty(), "{lines:?}");
    b
}

fn volume(k: &Kernel, b: &mut Brep) -> f64 {
    let store = AttributeStore::default();
    let attrs = PropertyAttributes {
        store: &store,
        key_of: &|_| ObjectKey(0),
    };
    brep_volume(k, b, 1e-12, &attrs).unwrap()
}

/// The faces of `b` with every vertex within `d` of `p`: at a set-back vertex, the fill's
/// patches (the rounds and the original faces reach the edges' far ends).
fn faces_near(b: &Brep, p: [f64; 3], d: f64) -> Vec<FaceId> {
    let p = Point3::new(p[0], p[1], p[2]);
    get_faces(b)
        .unwrap()
        .into_iter()
        .filter(|&f| {
            let vs = solidcore::topology::face::p4::vertices(b, f, None, false, None).unwrap();
            !vs.is_empty()
                && vs
                    .iter()
                    .all(|&v| b.vertex(v).unwrap().point.distance(p) <= d)
        })
        .collect()
}

/// A 20 x 20 x 10 block (`z` from 0) less a square pyramid standing under it (base 16 at
/// `z = -5`, apex at `z = 7`): a pit in its underside whose apex is a vertex of four concave
/// edges, each from the apex to a corner of the pit's mouth at `(±14/3, ±14/3, 0)`.
fn pit(k: &mut Kernel) -> Brep {
    let block = create_box(k, Point3::new(-10.0, -10.0, 0.0), 20.0, 20.0, 10.0).unwrap();
    let pyramid = create_pyramid(k, Point3::new(0.0, 0.0, -5.0), 16.0, 12.0).unwrap();
    boolean_difference(k, block, pyramid).unwrap()
}

const APEX: [f64; 3] = [0.0, 0.0, 7.0];
/// The pit's four concave edges, at their middles.
const PIT_AT: [[f64; 3]; 4] = [
    [7.0 / 3.0, 7.0 / 3.0, 3.5],
    [-7.0 / 3.0, 7.0 / 3.0, 3.5],
    [-7.0 / 3.0, -7.0 / 3.0, 3.5],
    [7.0 / 3.0, -7.0 / 3.0, 3.5],
];
/// Their length, from the apex to a corner of the mouth.
fn pit_edge_length() -> f64 {
    let c = 14.0 / 3.0;
    (2.0 * c * c + 49.0_f64).sqrt()
}

/// The material a concave round of radius r adds per length along a straight edge whose
/// void has angle `phi` between its faces: `r^2 (cot(phi / 2) - (pi - phi) / 2)`. Over the
/// whole edge it bounds what the round adds there (it stops short at its ends, and near the
/// apex the rounds' material overlaps).
fn concave_round_area(r: f64, phi: f64) -> f64 {
    r * r * (1.0 / (phi / 2.0).tan() - (std::f64::consts::PI - phi) / 2.0)
}

/// The angle of the void at the pit's concave edges: the pyramid's own interior dihedral
/// angle at a lateral edge, `pi - acos(a^2 / (h^2 + a^2))` (its faces' outward normals
/// `(h, 0, a)` and `(0, h, a)`) for half base a = 8 and height h = 12.
fn pit_dihedral() -> f64 {
    let (a, h): (f64, f64) = (8.0, 12.0);
    std::f64::consts::PI - ((a * a) / (h * h + a * a)).acos()
}

/// Four concave rounds of one radius at the pit's apex: no Ball (it needs three faces), so the
/// Setback, eight fill patches near the apex (four setback patches, four interior parts),
/// four-fold symmetric (each kind's areas equal), its volume between the sharp pit's and
/// that plus every round's whole-edge material.
#[test]
fn four_equal_concave_rounds_at_a_pit_apex_are_the_setback() {
    let mut k = kernel();
    let mut sharp = pit(&mut k);
    let v0 = volume(&k, &mut sharp);
    let b = pit(&mut k);
    let mut b = round(&mut k, b, &PIT_AT, &[1.0; 4]);
    let near = faces_near(&b, APEX, 5.0);
    assert_eq!(near.len(), 8, "fill patches near the apex");
    let mut areas: Vec<f64> = near
        .iter()
        .map(|&f| face_area(&k, &mut b, f, 1e-10).unwrap())
        .collect();
    areas.sort_by(f64::total_cmp);
    for g in areas.chunks(4) {
        for a in g {
            assert!((a - g[0]).abs() <= 1e-6 * g[0], "{areas:?}");
        }
    }
    let v = volume(&k, &mut b);
    let most = 4.0 * concave_round_area(1.0, pit_dihedral()) * pit_edge_length();
    assert!(v > v0 && v < v0 + most, "{v0} {v} {most}");
}

/// Four concave rounds of unequal radii (1, 1.25, 1.5, 1.75) at the pit's apex: the Setback,
/// its fill certified with each setback patch's k = -1 / its round's radius (the interior's
/// the tightest round's), eight patches near the apex.
#[test]
fn four_unequal_concave_rounds_at_a_pit_apex_are_the_setback() {
    let mut k = kernel();
    let mut sharp = pit(&mut k);
    let v0 = volume(&k, &mut sharp);
    let radii = [1.0, 1.25, 1.5, 1.75];
    let b = pit(&mut k);
    let mut b = round(&mut k, b, &PIT_AT, &radii);
    assert_eq!(faces_near(&b, APEX, 6.0).len(), 8);
    let v = volume(&k, &mut b);
    let most: f64 = radii
        .iter()
        .map(|&r| concave_round_area(r, pit_dihedral()) * pit_edge_length())
        .sum();
    assert!(v > v0 && v < v0 + most, "{v0} {v} {most}");
}

/// The box `[0, 20]^3` less `[10, 20]^3`.
fn notched(k: &mut Kernel) -> Brep {
    let part = create_box(k, Point3::new(0.0, 0.0, 0.0), 20.0, 20.0, 20.0).unwrap();
    let cut = create_box(k, Point3::new(10.0, 10.0, 10.0), 10.0, 10.0, 10.0).unwrap();
    boolean_difference(k, part, cut).unwrap()
}

/// The notch's vertex (20, 10, 10): its edge along x is concave (the notch's floor against its
/// wall), its edges along y and z convex (each against the outer face x = 20).
const MIXED: [f64; 3] = [20.0, 10.0, 10.0];
const MIXED_AT: [[f64; 3]; 3] = [[15.0, 10.0, 10.0], [20.0, 15.0, 10.0], [20.0, 10.0, 15.0]];

/// Mixed convexity with every edge blended, rounds r 1 (concave), 1.5 and 2 (convex): the
/// Setback, no spring curve between the concave round and either convex one (their trimlines
/// meet where their rails cross), its fill certified patch by patch (k = 0 beside the convex
/// rounds, -1 beside the concave one): six patches near the vertex.
#[test]
fn a_mixed_corner_with_every_edge_rounded_is_the_setback() {
    let mut k = kernel();
    let b = notched(&mut k);
    let b = round(&mut k, b, &MIXED_AT, &[1.0, 1.5, 2.0]);
    assert_eq!(faces_near(&b, MIXED, 6.0).len(), 6);
}

/// The same with rounds of one radius (r 1.5): still mixed, so the Setback (the Ball is for
/// one convexity).
#[test]
fn a_mixed_corner_of_equal_rounds_is_the_setback() {
    let mut k = kernel();
    let b = notched(&mut k);
    let b = round(&mut k, b, &MIXED_AT, &[1.5; 3]);
    assert_eq!(faces_near(&b, MIXED, 6.0).len(), 6);
}
