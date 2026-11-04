//! Fillets and chamfers against closed forms (the cases of SolidCraft's `core/blend_tests.rs`,
//! on edge, face and vertex indices). Their own test binary: each holds the process kernel for
//! seconds.
#![allow(
    clippy::unwrap_used,
    clippy::indexing_slicing,
    clippy::panic,
    clippy::expect_used
)] // test helpers

use std::f64::consts::PI;

use solidcore::geom::attribute::AttributeStore;
use solidcore::kernel::Kernel;
use solidcore::query::api_queries::{edge_closest_point, face_bounding_box, vertex_get_point};
use solidcore::topology::api_files::read_brep_from_file;
use solidcore::topology::brep::p10::TopologyEnv;
use solidcore::topology::loopuse;
use solidcore::topology::topology::NotifyRecords;
use solidcore::{Brep, Vec3};
use solidcore_ops::{
    blend, chamfer_sides, edge_faces_along, fillet_variable, make_box, session, Blend, Body,
    ChamferSecond, CheckLevel, Corners, Error, FitKind, SetBack,
};

fn all_fixes() {
    session::set_fixes(solidcore::kernel::UpstreamFixes::NAMES.iter().copied()).unwrap();
}

/// The 30 x 20 x 10 block at the origin.
fn block() -> Body {
    make_box([0.0; 3], [30.0, 20.0, 10.0]).unwrap()
}

const ABC: f64 = 6000.0;

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs().max(1.0)
}

fn volume(b: &Body) -> f64 {
    b.volume(1e-8).unwrap()
}

/// `b` as the kernel reads it back from its `.smb` bytes (the same entity order), with a
/// kernel to query it.
fn kernel_copy(b: &Body) -> (Kernel, Brep) {
    let mut k = Kernel::default();
    let mut attributes = AttributeStore::default();
    let mut records = NotifyRecords::new();
    let mut env = TopologyEnv {
        all_attributes: &[],
        attributes: &mut attributes,
        key_of: &loopuse::notify_ref,
        records: &mut records,
        next_key: 1 << 50,
    };
    let bytes = b.to_smb(false).unwrap();
    let brep = read_brep_from_file(&mut k, Some(bytes), false, false, &mut env).unwrap();
    (k, brep)
}

/// The index of the edge of `b` through `p`.
fn edge_at(b: &Body, p: [f64; 3]) -> usize {
    let (k, brep) = kernel_copy(b);
    let q = Vec3::new(p[0], p[1], p[2]);
    brep.edge_list()
        .iter()
        .position(|&e| {
            edge_closest_point(&k, &brep, e, q)
                .unwrap()
                .closest
                .distance
                < 1e-9
        })
        .unwrap()
}

/// The index of the vertex of `b` at `p`.
fn vertex_at(b: &Body, p: [f64; 3]) -> usize {
    let (_, brep) = kernel_copy(b);
    brep.vertex_list()
        .iter()
        .position(|&v| {
            let q = vertex_get_point(&brep, v).unwrap();
            (q.x - p[0]).abs() + (q.y - p[1]).abs() + (q.z - p[2]).abs() < 1e-9
        })
        .unwrap()
}

/// The block's top front edge, 30 long, between the top (z = 10) and the front (y = 0).
fn top_front(b: &Body) -> usize {
    edge_at(b, [15.0, 0.0, 10.0])
}

#[test]
fn rounds_and_chamfers_take_their_closed_forms() {
    all_fixes();
    // a round of r 2 along 30: (1 - pi/4) r^2 off per unit length
    let b = block();
    let e = top_front(&b);
    let made = blend(&b, Blend::Round, &[e], &[2.0], &[], Corners::Filled).unwrap();
    assert_eq!(made.dropped, vec![false]);
    assert!(made.holes.is_empty());
    assert!(rel(volume(&b), ABC - (1.0 - PI / 4.0) * 4.0 * 30.0) < 1e-9);
    // an equal chamfer of 2 along 30: d^2 / 2 off per unit length
    let b = block();
    blend(
        &b,
        Blend::Chamfer,
        &[top_front(&b)],
        &[2.0],
        &[],
        Corners::Filled,
    )
    .unwrap();
    assert!(rel(volume(&b), ABC - 2.0 * 30.0) < 1e-9);
    // a chord round of width 2 sqrt 2 between planes at 90 degrees is the round of r 2
    let b = block();
    blend(
        &b,
        Blend::Chord,
        &[top_front(&b)],
        &[2.0 * std::f64::consts::SQRT_2],
        &[],
        Corners::Filled,
    )
    .unwrap();
    assert!(rel(volume(&b), ABC - (1.0 - PI / 4.0) * 4.0 * 30.0) < 1e-8);
    // the curvature round builds and takes off the ball's area
    let b = block();
    blend(
        &b,
        Blend::Curvature,
        &[top_front(&b)],
        &[2.0],
        &[],
        Corners::Filled,
    )
    .unwrap();
    assert!(rel(volume(&b), ABC - (1.0 - PI / 4.0) * 4.0 * 30.0) < 1e-6);
}

#[test]
fn three_rounds_meet_at_a_corner_filled_set_back_or_open() {
    all_fixes();
    let edges = |b: &Body| {
        [
            edge_at(b, [15.0, 20.0, 10.0]),
            edge_at(b, [30.0, 10.0, 10.0]),
            edge_at(b, [30.0, 20.0, 5.0]),
        ]
    };
    // the ball's corner
    let b = block();
    blend(&b, Blend::Round, &edges(&b), &[2.0], &[], Corners::Filled).unwrap();
    assert!(b.check(CheckLevel::All).unwrap().is_valid());
    let ball = volume(&b);
    // set back 3 at the corner: a different corner, a closed solid
    let b = block();
    let corner = vertex_at(&b, [30.0, 20.0, 10.0]);
    let back = [SetBack {
        vertex: corner,
        distance: 3.0,
    }];
    blend(&b, Blend::Round, &edges(&b), &[2.0], &back, Corners::Filled).unwrap();
    assert!(b.check(CheckLevel::All).unwrap().is_valid());
    assert!(rel(volume(&b), ball) < 1e-3);
    // opened: the corner's hole, its edges closing one loop, each on a face outside it
    let b = block();
    let corner = vertex_at(&b, [30.0, 20.0, 10.0]);
    let back = [SetBack {
        vertex: corner,
        distance: 3.0,
    }];
    let made = blend(&b, Blend::Round, &edges(&b), &[2.0], &back, Corners::Open).unwrap();
    assert_eq!(made.holes.len(), 1);
    let [faces, edges_n, _] = b.counts().unwrap();
    for h in &made.holes[0].edges {
        assert!(h.edge < edges_n && h.face < faces);
    }
    assert!(made.holes[0].edges.len() >= 3);
}

#[test]
fn chamfers_by_sides_take_the_face_asked() {
    all_fixes();
    // 1 back along the top, 3 down the front: d1 d2 / 2 off per unit length either way round,
    // and the top's leg is 1 when the top is asked first
    for up in [true, false] {
        let b = block();
        let e = top_front(&b);
        let faces = edge_faces_along(&b, e, [0.0, 0.0, 1.0]).unwrap();
        let top = if faces[0].1 > faces[1].1 {
            faces[0].0
        } else {
            faces[1].0
        };
        let front = if faces[0].1 > faces[1].1 {
            faces[1].0
        } else {
            faces[0].0
        };
        let first = if up { top } else { front };
        chamfer_sides(&b, &[e], &[first], 1.0, ChamferSecond::Distance(3.0)).unwrap();
        assert!(rel(volume(&b), ABC - 1.5 * 30.0) < 1e-9);
        // the top keeps its leg: (30 x (20 - 1)) when asked first, else 20 - 3
        let leg = if up { 19.0 } else { 17.0 };
        let (k, made) = kernel_copy(&b);
        let top_area = made
            .face_list()
            .iter()
            .find_map(|&f| {
                let bx = face_bounding_box(&k, &made, f, true).ok()?;
                ((bx.min().z - 10.0).abs() < 1e-9 && (bx.max().z - 10.0).abs() < 1e-9)
                    .then(|| bx.max().y - bx.min().y)
            })
            .unwrap();
        assert!((top_area - leg).abs() < 1e-9, "{top_area}");
    }
    // a 45 degree chamfer 2 back along the top is the equal one
    let b = block();
    let e = top_front(&b);
    let top = edge_faces_along(&b, e, [0.0, 0.0, 1.0]).unwrap();
    let top = if top[0].1 > top[1].1 {
        top[0].0
    } else {
        top[1].0
    };
    chamfer_sides(&b, &[e], &[top], 2.0, ChamferSecond::Angle(PI / 4.0)).unwrap();
    assert!(rel(volume(&b), ABC - 2.0 * 30.0) < 1e-9);
}

#[test]
fn a_variable_round_runs_from_its_start() {
    all_fixes();
    // r from 1 to 2 along the top front edge: between the constant rounds of 1 and 2
    let b = block();
    let e = top_front(&b);
    let start = vertex_at(&b, [0.0, 0.0, 10.0]);
    fillet_variable(&b, e, start, 1, &[0.0, 0.0, 1.0, 1.0], &[1.0, 2.0]).unwrap();
    let v = volume(&b);
    let off = |r: f64| (1.0 - PI / 4.0) * r * r * 30.0;
    assert!(v < ABC - off(1.0) && v > ABC - off(2.0), "{v}");
    assert!(b.check(CheckLevel::All).unwrap().is_valid());

    // a law with a corner, one zero inside, one that starts off the edge
    let b = block();
    let e = top_front(&b);
    let start = vertex_at(&b, [0.0, 0.0, 10.0]);
    let bent = fillet_variable(
        &b,
        e,
        start,
        1,
        &[0.0, 0.0, 0.5, 1.0, 1.0],
        &[1.0, 2.0, 1.0],
    );
    assert!(matches!(bent, Err(Error::Unsupported(_))), "{bent:?}");
    let pinched = fillet_variable(
        &b,
        e,
        start,
        1,
        &[0.0, 0.0, 0.5, 1.0, 1.0],
        &[1.0, 0.0, 1.0],
    );
    assert!(matches!(pinched, Err(Error::Unsupported(_))), "{pinched:?}");
    let away = vertex_at(&b, [30.0, 20.0, 0.0]);
    let off_edge = fillet_variable(&b, e, away, 1, &[0.0, 0.0, 1.0, 1.0], &[1.0, 2.0]);
    assert_eq!(off_edge, Err(Error::Argument));
    // refused: the body as it came in
    assert!(rel(volume(&b), ABC) < 1e-12);
}

#[test]
fn arguments_are_checked_before_the_kernel() {
    all_fixes();
    let b = block();
    let e = top_front(&b);
    let round = |edges: &[usize], sizes: &[f64]| {
        blend(&b, Blend::Round, edges, sizes, &[], Corners::Filled)
    };
    assert_eq!(round(&[], &[1.0]), Err(Error::Argument));
    assert_eq!(round(&[e, e], &[1.0]), Err(Error::Argument));
    assert_eq!(round(&[999], &[1.0]), Err(Error::Argument));
    assert_eq!(round(&[e], &[0.0]), Err(Error::Argument));
    assert_eq!(round(&[e], &[f64::NAN]), Err(Error::Argument));
    assert_eq!(round(&[e], &[2e6]), Err(Error::Argument));
    assert_eq!(round(&[e], &[1.0, 2.0]), Err(Error::Argument));
    let corner = vertex_at(&b, [30.0, 0.0, 10.0]);
    let back = [SetBack {
        vertex: corner,
        distance: 1.0,
    }];
    assert!(matches!(
        blend(&b, Blend::Chamfer, &[e], &[1.0], &back, Corners::Filled),
        Err(Error::Unsupported(_))
    ));
    // a setback longer than its edges leaves no room
    let long = [SetBack {
        vertex: corner,
        distance: 50.0,
    }];
    assert_eq!(
        blend(&b, Blend::Round, &[e], &[1.0], &long, Corners::Filled),
        Err(Error::NoRoom)
    );
    assert_eq!(
        chamfer_sides(&b, &[e], &[], 1.0, ChamferSecond::Distance(1.0)),
        Err(Error::Argument)
    );
    assert_eq!(
        chamfer_sides(&b, &[e], &[0], 1.0, ChamferSecond::Angle(PI)),
        Err(Error::Argument)
    );
    assert_eq!(
        edge_faces_along(&b, e, [f64::NAN, 0.0, 0.0]),
        Err(Error::Argument)
    );
    assert_eq!(
        edge_faces_along(&b, 999, [0.0, 0.0, 1.0]),
        Err(Error::Argument)
    );
    // nothing changed the body
    assert!(rel(volume(&b), ABC) < 1e-12);
}

#[test]
fn a_round_larger_than_its_faces_is_refused_with_the_size_that_fits() {
    all_fixes();
    // a 2 thick wall, both its top edges rounded 1.5: the two rounds' bands overlap on the top
    let b = make_box([0.0; 3], [30.0, 2.0, 10.0]).unwrap();
    let edges = [
        edge_at(&b, [15.0, 0.0, 10.0]),
        edge_at(&b, [15.0, 2.0, 10.0]),
    ];
    let r = blend(&b, Blend::Round, &edges, &[1.5], &[], Corners::Filled);
    match r {
        Err(Error::TooLarge { size, kind, at, .. }) => {
            assert_eq!(kind, FitKind::Radius);
            assert!(size > 0.0 && size <= 1.0 + 1e-9, "{size}");
            let at = at.unwrap();
            assert!((at[2] - 10.0).abs() < 1e-9, "{at:?}");
        }
        Err(Error::NoRoom) => {}
        other => panic!("{other:?}"),
    }
    // refused: the body as it came in
    assert!(rel(volume(&b), 600.0) < 1e-12);
}
