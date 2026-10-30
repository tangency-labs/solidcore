//! The blends: rounds, chamfers and blends of each kind on a box against analytic volumes,
//! the variable round, the selection errors and hostile arguments, the corner setback
//! (closed and open) and the unequal chamfer of a cylinder's rim.
//!
//! Source: tangency-labs/smlib `capi/tests/ops.c` (the blend half of `main`) and
//! `capi/tests/ops_kernel.cpp` (`corner_setback`, `corner_setback_open`, the rim chamfers
//! of the round-wider-than-its-face cases), Apache-2.0, Copyright SolidCraft contributors,
//! replay/1.0.1 at `f8fd83e2f`. The open corner's flat-array and history checks need later
//! families and are left for those.

#![allow(clippy::unwrap_used, clippy::panic)] // test helpers

use solidcore_smlib_api::{
    Brep, Smlib, Status, BLEND_CHAMFER, BLEND_CHORD, BLEND_CURVATURE, BLEND_ROUND,
};

const PI: f64 = core::f64::consts::PI;
const QUARTER: f64 = 1.0 - PI / 4.0;

fn bx(api: &mut Smlib, origin: [f64; 3], size: [f64; 3]) -> Brep {
    api.make_box(origin, size).unwrap()
}

fn volume(api: &mut Smlib, b: &Brep) -> f64 {
    api.brep_volume(b, 1e-7).unwrap()
}

/// SMLib's own validity check, every test (level 2), finds nothing wrong with `b`.
fn sound(api: &mut Smlib, b: &Brep) -> bool {
    let (report, bad) = api.brep_check(b, 2).unwrap();
    if bad != 0 {
        eprintln!("  {bad} failed validity tests:\n{report}");
    }
    bad == 0
}

fn near(got: f64, want: f64, rel: f64) -> bool {
    let ok = (got - want).abs() <= rel * want.abs();
    if !ok {
        eprintln!("  got {got:.9e}, want {want:.9e}");
    }
    ok
}

const EDGE: [f64; 3] = [10.0, 20.0, 20.0];
const TOPS: [[f64; 3]; 2] = [[10.0, 20.0, 20.0], [10.0, 0.0, 20.0]];

#[test]
fn box_rounds_chamfers_and_blends() {
    let mut api = Smlib::default();
    let a = bx(&mut api, [0.0; 3], [20.0; 3]);
    let r = api.chamfer(&a, &[EDGE], 0.1, 3.0).unwrap();
    assert!(near(volume(&mut api, &r), 7910.0, 1e-6));
    api.fillet(&a, &[EDGE], 0.1, 2.0).unwrap();

    // Two parallel top edges rounded at 2 and 3 lose (1 - pi/4) r^2 L each.
    let r = api.blend(&a, BLEND_ROUND, &TOPS, &[2.0, 3.0], 0.1).unwrap();
    assert!((volume(&mut api, &r) - (8000.0 - QUARTER * 13.0 * 20.0)).abs() < 1e-3);
    let r = api
        .blend(&a, BLEND_CHAMFER, &TOPS, &[2.0, 3.0], 0.1)
        .unwrap();
    assert!((volume(&mut api, &r) - (8000.0 - 0.5 * 13.0 * 20.0)).abs() < 1e-3);
    // A chord of 2 sqrt(2) on a right-angled edge is the round of radius 2.
    let chord = 2.0 * 2.0_f64.sqrt();
    let r = api.blend(&a, BLEND_CHORD, &[EDGE], &[chord], 0.1).unwrap();
    assert!((volume(&mut api, &r) - (8000.0 - QUARTER * 4.0 * 20.0)).abs() < 1e-2);
    // A curvature-continuous blend runs outside the round and takes less, but something.
    let r = api
        .blend(&a, BLEND_CURVATURE, &[EDGE], &[2.0], 0.1)
        .unwrap();
    let g2 = volume(&mut api, &r);
    assert!(
        g2 > 8000.0 - QUARTER * 4.0 * 20.0 && g2 < 8000.0 - 1.0,
        "{g2}"
    );

    assert_eq!(
        api.blend(&a, 9, &[EDGE], &[2.0], 0.1).err(),
        Some(Status::ERR_ARG)
    );
    assert_eq!(
        api.blend(&a, BLEND_ROUND, &[EDGE], &[], 0.1).err(),
        Some(Status::ERR_ARG)
    );
    assert_eq!(
        api.blend(&a, BLEND_ROUND, &TOPS, &[1.0, -1.0], 0.1).err(),
        Some(Status::ERR_ARG)
    );
    assert!(
        near(volume(&mut api, &a), 8000.0, 1e-9),
        "the input is unchanged"
    );
}

#[test]
fn variable_rounds() {
    let mut api = Smlib::default();
    let a = bx(&mut api, [0.0; 3], [20.0; 3]);
    // A radius running linearly from 3 at x = 0 to 1 at x = 20: (1 - pi/4) L (r0^2 + r0 r1
    // + r1^2) / 3 taken away, whichever way the edge runs.
    let knots = [0.0, 0.0, 1.0, 1.0];
    let from0 = [0.0, 20.0, 20.0];
    let from1 = [20.0, 20.0, 20.0];
    let linear = 8000.0 - QUARTER * 20.0 * 13.0 / 3.0;
    let r = api
        .fillet_variable(&a, EDGE, 0.1, from0, 1, &knots, &[3.0, 1.0])
        .unwrap();
    let v0 = volume(&mut api, &r);
    assert!((v0 - linear).abs() < 0.05, "{v0}");
    let r = api
        .fillet_variable(&a, EDGE, 0.1, from1, 1, &knots, &[3.0, 1.0])
        .unwrap();
    assert!((volume(&mut api, &r) - v0).abs() < 1e-6 * 8000.0);
    // A constant law is the round of radius 2; a quadratic 3, 2, 1 is the linear law.
    let round = 8000.0 - QUARTER * 4.0 * 20.0;
    let r = api
        .fillet_variable(&a, EDGE, 0.1, from0, 1, &knots, &[2.0, 2.0])
        .unwrap();
    assert!((volume(&mut api, &r) - round).abs() < 1e-6 * 8000.0);
    let knots2 = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
    let r = api
        .fillet_variable(&a, EDGE, 0.1, from0, 2, &knots2, &[2.0, 2.0, 2.0])
        .unwrap();
    assert!((volume(&mut api, &r) - round).abs() < 1e-6 * 8000.0);
    let r = api
        .fillet_variable(&a, EDGE, 0.1, from0, 2, &knots2, &[3.0, 2.0, 1.0])
        .unwrap();
    assert!((volume(&mut api, &r) - v0).abs() < 1e-6 * 8000.0);

    let bad = |api: &mut Smlib, degree: i32, knots: &[f64]| {
        api.fillet_variable(&a, EDGE, 0.1, from0, degree, knots, &[3.0, 1.0])
            .err()
    };
    assert_eq!(
        bad(&mut api, 1, &[0.0, 0.5, 1.0, 1.0]),
        Some(Status::ERR_ARG)
    );
    assert_eq!(bad(&mut api, 4, &knots), Some(Status::ERR_ARG));
    assert_eq!(bad(&mut api, 1, &[0.0, 0.0, 1.0]), Some(Status::ERR_ARG));
    // A law that pinches to zero inside the edge.
    assert_eq!(
        api.fillet_variable(
            &a,
            EDGE,
            0.1,
            from0,
            1,
            &[0.0, 0.0, 0.5, 1.0, 1.0],
            &[1.0, 0.0, 1.0]
        )
        .err(),
        Some(Status::ERR_UNSUPPORTED)
    );
}

#[test]
fn selection_errors_and_hostile_arguments() {
    let mut api = Smlib::default();
    let a = bx(&mut api, [0.0; 3], [20.0; 3]);
    // Selection errors name the point.
    let two = [EDGE, [100.0, 100.0, 100.0]];
    assert_eq!(
        api.fillet(&a, &two, 0.1, 1.0).err(),
        Some(Status::ERR_PICK_NONE)
    );
    assert_eq!(api.last_failed_point(), Some(1));
    assert_eq!(
        api.fillet(&a, &[[20.0, 20.0, 20.0]], 0.1, 1.0).err(),
        Some(Status::ERR_PICK_AMBIGUOUS)
    );
    assert_eq!(api.last_failed_point(), Some(0));

    let arg = Some(Status::ERR_ARG);
    assert_eq!(api.fillet(&a, &[], 0.1, 1.0).err(), arg);
    assert_eq!(api.fillet(&a, &vec![EDGE; 1 << 20], 0.1, 1.0).err(), arg);
    assert_eq!(api.fillet(&a, &[[f64::NAN, 0.0, 0.0]], 0.1, 1.0).err(), arg);
    assert_eq!(api.fillet(&a, &[EDGE], f64::NAN, 1.0).err(), arg);
    assert_eq!(api.fillet(&a, &[EDGE], 0.1, f64::NAN).err(), arg);
    assert_eq!(api.fillet(&a, &[EDGE], 0.1, -1.0).err(), arg);
    assert_eq!(api.chamfer(&a, &[EDGE], 0.1, f64::INFINITY).err(), arg);

    // Sizes SMLib cannot build: any status, and a success must be a real body.
    if let Ok(r) = api.fillet(&a, &[EDGE], 0.1, 50.0) {
        assert!(!r.face_list().is_empty());
    }
    if let Ok(r) = api.chamfer(&a, &[EDGE], 0.1, 30.0) {
        assert!(!r.face_list().is_empty());
    }
    assert!(near(volume(&mut api, &a), 8000.0, 1e-9));
}

/// The corner where the box's three edges at (60, 40, 20) meet, each rounded at 5.
const SETBACK_EDGES: [[f64; 3]; 3] = [[30.0, 40.0, 20.0], [60.0, 20.0, 20.0], [60.0, 40.0, 10.0]];
const CORNER: [f64; 3] = [60.0, 40.0, 20.0];

#[test]
fn corner_setback() {
    let mut api = Smlib::default();
    let part = bx(&mut api, [0.0; 3], [60.0, 40.0, 20.0]);
    let radii = [5.0; 3];
    let r = api
        .blend_setback(
            &part,
            BLEND_ROUND,
            &SETBACK_EDGES,
            &radii,
            &[CORNER],
            &[8.0],
            0.1,
        )
        .unwrap();
    let in_c_sphere = 512.0 - 225.0 * QUARTER - 125.0 * (1.0 - PI / 6.0);
    let outside = 47377.120025868 - in_c_sphere;
    let hull = 299.0201;
    let prisms = 512.0 - 225.0 * QUARTER - 125.0 * (2.0_f64.sqrt() - 1.0);
    let v = volume(&mut api, &r);
    assert!(v > outside + hull && v < outside + prisms, "{v}");
    // The box's six faces, the three rounds and the patch's six faces.
    assert_eq!(api.brep_counts(&r).unwrap(), [15, 36, 23]);
    let (lo, hi) = api.brep_bounds(&r).unwrap();
    for k in 0..3 {
        assert!(lo[k].abs() <= 1e-9 && (hi[k] - CORNER[k]).abs() <= 1e-9);
    }
    assert!(sound(&mut api, &r));

    // A point that picks no vertex, a vertex no blended edge meets, a setback as long as
    // an edge there, a chamfer, a corner with an edge not rounded.
    assert_eq!(
        api.blend_setback(
            &part,
            BLEND_ROUND,
            &SETBACK_EDGES,
            &radii,
            &[[59.0, 39.0, 19.0]],
            &[8.0],
            0.1
        )
        .err(),
        Some(Status::ERR_PICK_NONE)
    );
    assert_eq!(api.last_failed_corner(), Some(0));
    assert_eq!(
        api.blend_setback(
            &part,
            BLEND_ROUND,
            &SETBACK_EDGES,
            &radii,
            &[[0.0; 3]],
            &[8.0],
            0.1
        )
        .err(),
        Some(Status::ERR_ARG)
    );
    assert_eq!(
        api.blend_setback(
            &part,
            BLEND_ROUND,
            &SETBACK_EDGES,
            &radii,
            &[CORNER],
            &[20.0],
            0.1
        )
        .err(),
        Some(Status::ERR_NO_ROOM)
    );
    assert_eq!(
        api.blend_setback(
            &part,
            BLEND_CHAMFER,
            &SETBACK_EDGES,
            &radii,
            &[CORNER],
            &[8.0],
            0.1
        )
        .err(),
        Some(Status::ERR_UNSUPPORTED)
    );
    assert_eq!(
        api.blend_setback(
            &part,
            BLEND_ROUND,
            &SETBACK_EDGES[..2],
            &radii[..2],
            &[CORNER],
            &[8.0],
            0.1
        )
        .err(),
        Some(Status::ERR_UNSUPPORTED)
    );
    assert_eq!(api.last_failed_corner(), Some(0));
}

#[test]
fn corner_setback_open() {
    let mut api = Smlib::default();
    let part = bx(&mut api, [0.0; 3], [60.0, 40.0, 20.0]);
    let (r, holes) = api
        .blend_setback_open(
            &part,
            BLEND_ROUND,
            &SETBACK_EDGES,
            &[5.0; 3],
            &[CORNER],
            &[8.0],
            0.1,
        )
        .unwrap();
    // One hole of twelve edges: its count, then a triple each.
    assert_eq!(holes.len(), 1 + 3 * 12);
    assert_eq!(holes[0], 12);
    assert_eq!(api.brep_counts(&r).unwrap(), [9, 30, 22]);
    assert!(!api.brep_is_manifold_solid(&r).unwrap());
    // Each edge of the hole starts where the one before it ends: every index is an edge
    // and a face of the result.
    let [faces, edges, _] = api.brep_counts(&r).unwrap();
    for t in holes[1..].chunks(3) {
        assert!((0..edges).contains(&t[0]) && (0..faces).contains(&t[1]));
        assert!(t[2] == 0 || t[2] == 1);
    }
}

#[test]
fn rim_chamfers_of_a_cylinder() {
    let mut api = Smlib::default();
    let can = api.make_cylinder([0.0; 3], 10.0, 20.0).unwrap();
    let whole = PI * 100.0 * 20.0;
    let at = 100.0 * PI / 180.0;
    let top_rim = [[10.0 * at.cos(), 10.0 * at.sin(), 20.0]];
    let up = [0.0, 0.0, 1.0];
    let down = [0.0, 0.0, -1.0];
    // The ring the triangle with corners (r, z) sweeps about the axis (Pappus).
    let ring = |r0: f64, z0: f64, r1: f64, z1: f64, r2: f64, z2: f64| {
        let area = 0.5 * ((r1 - r0) * (z2 - z0) - (r2 - r0) * (z1 - z0)).abs();
        2.0 * PI * (r0 + r1 + r2) / 3.0 * area
    };
    for d in [10.0, 9.99] {
        let equal = ring(10.0, 20.0, 10.0 - d, 20.0, 10.0, 20.0 - d);
        let sides = ring(10.0, 20.0, 10.0 - d, 20.0, 10.0, 17.0);
        let slope = (PI / 6.0).tan();
        let r = api.chamfer(&can, &top_rim, 0.1, d);
        rim_case(&mut api, "chamfer", d, r, whole - equal);
        let r = api.blend(&can, BLEND_CHAMFER, &top_rim, &[d], 0.1);
        rim_case(&mut api, "blend chamfer", d, r, whole - equal);
        let r = api.chamfer_sides(&can, &top_rim, 0.1, d, 3.0, 0, up);
        rim_case(
            &mut api,
            "chamfer across the top, 3 down",
            d,
            r,
            whole - sides,
        );
        let r = api.chamfer_sides(&can, &top_rim, 0.1, 3.0, d, 0, down);
        rim_case(
            &mut api,
            "chamfer 3 down, across the top",
            d,
            r,
            whole - sides,
        );
        let r = api.chamfer_sides(&can, &top_rim, 0.1, d, PI / 6.0, 1, up);
        let at_30 = ring(10.0, 20.0, 10.0 - d, 20.0, 10.0, 20.0 - d * slope);
        rim_case(
            &mut api,
            "chamfer across the top at 30 degrees",
            d,
            r,
            whole - at_30,
        );
    }
}

/// A rim chamfer of width `d`: at 9.99 a sound body of volume `want`; at 10, as wide as the
/// top, refused for want of room.
fn rim_case(api: &mut Smlib, what: &str, d: f64, r: solidcore_smlib_api::Result<Brep>, want: f64) {
    if d < 10.0 {
        let r = r.unwrap_or_else(|s| panic!("{what} {d}: status {s}"));
        assert!(near(volume(api, &r), want, 1e-7), "{what} {d}");
        assert!(sound(api, &r), "{what} {d}");
    } else {
        assert_eq!(r.err(), Some(Status::ERR_NO_ROOM), "{what} {d}");
    }
}
