//! `offset_analytic_mate_branch_from_vertex`: a draft (a local extended offset) of one face of
//! a block whose tapered, ROUNDED-corner opening cut crosses an inner pocket wall.
//!
//! The opening's corner arc sweeps a cone; where the opening's side wall crosses the pocket's
//! wall, that plane touches the cone along the ruling where the cone meets the side wall. The
//! pocket wall x the cone edge therefore starts at a singular point: the marcher in
//! `MateEdges` cannot begin there and finds only a branch on the faces' tangent extensions,
//! and the edge is reached by extrapolating that branch past its end -- missing, by up to
//! ~2e-5, the vertex the neighbouring edges agree on. The draft is then refused by the check
//! (a vertex/edge Bad Gap) at some opening widths and not others. The exact intersection of the
//! analytic surfaces starts at the edge's own vertex; with the fix its part on the true faces
//! is the edge.
//!
//! Invented geometry (z is the opening's axis): a block 50 x 50 x 30 from z = 2; a pocket
//! |x| <= W (W = 15), 21 <= y <= 28, 8 <= z <= 22; the opening a rectangle |x| <= Y, 16 <= y <= 28
//! in the plane z = 0 with R 4 corners, tapered 2 degrees inward over 34; then the side face
//! x = -25 drafted 1.5 degrees about its line at z = 2.5. Swept over Y - W = 0.35 .. 0.75.
//! Control: the same with sharp opening corners (no cone, no singular start).
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers

#[path = "support/offset_fixture.rs"]
mod offset_fixture;

use offset_fixture::{planar_faces_at, Notifications};
use solidcore::boolean::api_brep::boolean_difference;
use solidcore::curve::curve::Curve;
use solidcore::curve::ellipse::Ellipse;
use solidcore::curve::line::Line;
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::geom::axis2_placement::Axis2Placement;
use solidcore::geom::extent1d::Extent1d;
use solidcore::mass::integrator::properties::{brep_volume, PropertyAttributes};
use solidcore::math::{Point3, Vec3};
use solidcore::primitive::api_primitives::create_box;
use solidcore::sweep::api_brep::{create_taper_extrude_with_corner, OffsetCornerType, TaperExtrude};
use solidcore::topology::brep::local_operation::local_operation;
use solidcore::topology::brep::validation::{check_brep, CheckLevel};
use solidcore::{Brep, Kernel};

const WALL: f64 = 15.0;
const R: f64 = 4.0;
const TAPER: f64 = 2.0;
const DRAFT: f64 = 1.5;
const STEP0: f64 = 0.35;
/// The pocket's back wall, clear of the opening's corner cones (they span 16 <= y <= 20).
const BACK: f64 = 21.0;

fn kernel(fix: bool) -> Kernel {
    let mut k = Kernel::default();
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/kernel_cases/fixes.json");
    let names: Vec<String> = serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    for name in &names {
        *k.fixes.flag_mut(name).unwrap() = true;
    }
    k.fixes.offset_analytic_mate_branch_from_vertex = fix;
    k
}

fn line(x0: f64, y0: f64, x1: f64, y1: f64) -> Curve {
    Curve::from(Line::create_line_segment(3, Point3::new(x0, y0, 0.0), Point3::new(x1, y1, 0.0), None).unwrap())
}

/// A quarter circle about (cx, cy) in z = 0 from the direction `from`, counter-clockwise.
fn quarter(k: &Kernel, cx: f64, cy: f64, from: Vec3) -> Curve {
    let frame = Axis2Placement::new(Point3::new(cx, cy, 0.0), from, Vec3::Z.cross(from));
    let ivl = Extent1d::from_min_max(0.0, 90.0);
    Curve::from(Ellipse::create_canonical_circle(k, &frame, R, Some(&ivl), None).unwrap())
}

/// The opening profile |x| <= y_half, 16 <= y <= 28, counter-clockwise; rounded or sharp.
fn profile(k: &Kernel, y_half: f64, rounded: bool) -> Vec<Curve> {
    let (lo, hi) = (16.0, 28.0);
    if !rounded {
        return vec![
            line(-y_half, lo, y_half, lo),
            line(y_half, lo, y_half, hi),
            line(y_half, hi, -y_half, hi),
            line(-y_half, hi, -y_half, lo),
        ];
    }
    let (a, b) = (y_half - R, -y_half + R);
    vec![
        line(b, lo, a, lo),
        quarter(k, a, lo + R, -Vec3::Y),
        line(y_half, lo + R, y_half, hi - R),
        quarter(k, a, hi - R, Vec3::X),
        line(a, hi, b, hi),
        quarter(k, b, hi - R, Vec3::Y),
        line(-y_half, hi - R, -y_half, lo + R),
        quarter(k, b, lo + R, -Vec3::X),
    ]
}

fn volume(k: &Kernel, b: &mut Brep) -> f64 {
    let store = AttributeStore::default();
    let attrs = PropertyAttributes {
        store: &store,
        key_of: &|_| ObjectKey(0),
    };
    brep_volume(k, b, 1e-9, &attrs).unwrap()
}

/// `back`: the pocket's back wall y, the pocket running from there to y = 28.
fn body(k: &mut Kernel, y_half: f64, rounded: bool, back: f64) -> Brep {
    let block = create_box(k, Point3::new(-25.0, 0.0, 2.0), 50.0, 50.0, 30.0).unwrap();
    let pocket = create_box(k, Point3::new(-WALL, back, 8.0), 2.0 * WALL, 28.0 - back, 14.0).unwrap();
    let block = boolean_difference(k, block, pocket).unwrap();
    let cs = profile(k, y_half, rounded);
    let refs: Vec<Option<&Curve>> = cs.iter().map(Some).collect();
    let tool = match create_taper_extrude_with_corner(k, &refs, 34.0, TAPER, 3, OffsetCornerType::CurveExtension).unwrap() {
        TaperExtrude::Solid(b) => b,
        TaperExtrude::Refused(r) => panic!("opening refused: {r:?}"),
    };
    boolean_difference(k, block, tool).unwrap()
}

/// The side face x = -25 tilted 1.5 degrees about its line at z = 2.5: Ok(volume change) when the
/// result passes the check, Err(the check's lines or the refusal) otherwise.
fn draft(k: &mut Kernel, part: &mut Brep) -> Result<f64, String> {
    let before = volume(k, part);
    let side = planar_faces_at(k, part, 0, -25.0);
    assert_eq!(side.len(), 1, "the side face x = -25");
    let s = part.face(side[0]).unwrap().surface.as_ref().unwrap();
    let mut tilted = s.copy(k).unwrap();
    let mut turn = Axis2Placement::IDENTITY;
    turn.rotate_about_axis_at_point(-DRAFT.to_radians(), Point3::new(-25.0, 0.0, 2.5), Vec3::Y);
    tilted.transform(k, &turn, None).unwrap();
    let mut notifications = Notifications::new(part);
    local_operation(k, part, &side, vec![tilted], &mut notifications.context()).map_err(|e| format!("{e}"))?;
    let lines = check_brep(k, part, CheckLevel::Level2).unwrap().lines();
    if !lines.is_empty() {
        return Err(format!("{} check line(s): {}", lines.len(), lines[0]));
    }
    Ok(volume(k, part) - before)
}

fn sweep(fix: bool, rounded: bool, back: f64) -> Vec<Result<f64, String>> {
    (0..=20)
        .map(|i| {
            let mut k = kernel(fix);
            let y_half = WALL + STEP0 + 0.02 * f64::from(i);
            let mut part = body(&mut k, y_half, rounded, back);
            draft(&mut k, &mut part)
        })
        .collect()
}


/// Without the fix some widths are refused (measured: 4 of 21, each a vertex/edge Bad Gap of
/// ~2e-5); with it none is. Where the unfixed draft already built, the fixed body is the same
/// one, bit for bit; a rescued draft changes the volume as the built ones do (the drafted face is
/// far from the opening, so every width's change is the same to within ~3e-6, measured).
#[test]
fn rounded_opening_drafts_build_with_the_branch_from_the_vertex() {
    let off = sweep(false, true, BACK);
    let on = sweep(true, true, BACK);
    let refused_off = off.iter().filter(|r| r.is_err()).count();
    assert!(refused_off > 0, "the defect reproduces without the fix");
    let mut built: Vec<f64> = off.iter().filter_map(|r| r.as_ref().ok().copied()).collect();
    built.sort_by(f64::total_cmp);
    let median = built[built.len() / 2];
    for (i, (r, o)) in on.iter().zip(off.iter()).enumerate() {
        let d = r.as_ref().unwrap_or_else(|e| panic!("width {i} refused with the fix: {e}"));
        match o {
            Ok(b) => assert_eq!(d, b, "width {i}: the fix changed a draft that already built"),
            Err(_) => assert!((d - median).abs() < 5e-6, "width {i}: rescued draft dV {d} vs the built ones' {median}"),
        }
    }
}

/// The control: sharp opening corners are unaffected by the fix, bit for bit, and all build.
#[test]
fn sharp_opening_drafts_are_unchanged_by_the_fix() {
    let off = sweep(false, false, BACK);
    let on = sweep(true, false, BACK);
    for (i, (a, b)) in off.iter().zip(on.iter()).enumerate() {
        assert_eq!(a.as_ref().ok(), b.as_ref().ok(), "sharp width {i}");
        assert!(a.is_ok(), "sharp width {i} refused: {a:?}");
    }
}

/// The LIMIT of the fix, kept as a measured record (ignored: it asserts a known remaining defect).
/// With the pocket's back wall at y = 19 it runs through the corner cones (16 <= y <= 20), and the
/// back wall, the pocket's side wall and the cone meet at a tangent triple point. There the
/// marched piece passes within ~3e-6 of the vertex (it does not miss it), so the fix does not
/// apply, and the check refuses on an edge/vertex gap of ~2e-5 from the trim at that point: a
/// separate defect. Measured: 10 of 21 refused without the fix, 8 with it, none newly refused,
/// and every body that built without it is the same.
#[test]
#[ignore = "documents the tangent triple-point class the fix does not cover"]
fn back_wall_through_the_corner_cones_is_a_separate_defect() {
    let off = sweep(false, true, 19.0);
    let on = sweep(true, true, 19.0);
    assert_eq!(off.iter().filter(|r| r.is_err()).count(), 10);
    assert_eq!(on.iter().filter(|r| r.is_err()).count(), 8);
    for (i, (r, o)) in on.iter().zip(off.iter()).enumerate() {
        if let Ok(b) = o {
            assert_eq!(r.as_ref().ok(), Some(b), "width {i}: the fix changed or refused a draft that built");
        }
    }
}
