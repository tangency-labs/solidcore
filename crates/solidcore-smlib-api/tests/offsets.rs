//! The offsets: shells and their cavity guard, drafts, thickened sheets, offset sheets and
//! merged faces.
//!
//! Source: tangency-labs/smlib `capi/tests/ops_kernel.cpp` (Apache-2.0, Copyright SolidCraft
//! contributors), branch replay/1.0.1: the inward walls refused before SMLib sees them (the
//! cube's lid), the open box with walls 4 thick (the slab), the draft of the cube's side and
//! the rod's wall, `drafted_cylinders_are_cones`, the thicken argument checks and the L
//! sheet (rig shell/5326), and `offset_sheet` of a sphere's one face (as a solid and as a
//! closed sheet of void regions). Its sheets come from `smlib_extrude` (a later family);
//! here they are made with solidcore's linear sweep. The C API's tests have no
//! `merge_faces` case: that one is new.

#![allow(clippy::unwrap_used, clippy::panic)] // test helpers

use solidcore::curve::curve::Curve;
use solidcore::math::{Point3, Vec3};
use solidcore::primitive::api_curves::create_line_segment;
use solidcore::sweep::api_brep::create_linear_sweep;
use solidcore_smlib_api::{Brep, Smlib, Status, BOOLEAN_SUBTRACT, BOOLEAN_UNION};

fn bx(api: &mut Smlib, lo: [f64; 3], hi: [f64; 3]) -> Brep {
    let [x0, y0, z0] = lo;
    let [x1, y1, z1] = hi;
    api.make_box(lo, [x1 - x0, y1 - y0, z1 - z0]).unwrap()
}

fn volume(api: &mut Smlib, b: &Brep) -> f64 {
    api.brep_volume(b, 1e-8).unwrap()
}

/// SMLib's own validity check, every test (level 2), finds nothing wrong with `b`.
fn sound(api: &mut Smlib, b: &Brep) -> bool {
    let (report, bad) = api.brep_check(b, 2).unwrap();
    if bad != 0 {
        eprintln!("  {bad} failed validity tests:\n{report}");
    }
    bad == 0
}

/// An open sheet: the polyline through `points` (in z = 0) swept 20 up, uncapped.
fn sheet(api: &mut Smlib, points: &[[f64; 3]]) -> Brep {
    let curves: Vec<Curve> = points
        .windows(2)
        .map(|w| {
            let [a, b] = [w[0], w[1]];
            Curve::from(
                create_line_segment(Point3::new(a[0], a[1], a[2]), Point3::new(b[0], b[1], b[2]))
                    .unwrap(),
            )
        })
        .collect();
    let refs: Vec<Option<&Curve>> = curves.iter().map(Some).collect();
    create_linear_sweep(
        api.kernel_mut(),
        &refs,
        Vec3::new(0.0, 0.0, 1.0),
        20.0,
        false,
    )
    .unwrap()
}

#[test]
fn shell_refuses_inward_walls_that_leave_no_cavity() {
    let mut api = Smlib::default();
    let cube = bx(&mut api, [0.0; 3], [20.0; 3]);
    let lid = [[10.0, 10.0, 20.0]];
    // walls as thick as the cube's inradius or thicker: no cavity can remain
    assert_eq!(
        api.shell(&cube, &lid, 0.1, -10.0).err(),
        Some(Status::ERR_NO_ROOM)
    );
    assert_eq!(
        api.shell(&cube, &lid, 0.1, -15.0).err(),
        Some(Status::ERR_NO_ROOM)
    );
    assert_eq!(
        api.shell(&cube, &[], 0.1, -10.0).err(),
        Some(Status::ERR_NO_ROOM)
    );
    // (the C API's test takes this one only when it succeeds)
    if let Ok(r) = api.shell(&cube, &lid, 0.1, -9.5) {
        let v = volume(&mut api, &r);
        assert!((v - (8000.0 - 10.5)).abs() <= 1e-6 * 8000.0, "volume {v}");
    }
}

#[test]
fn open_box_with_walls_four_thick() {
    // A box 20 x 10 x 30 open at the top with walls 4 thick: the sides need room 4 from
    // each other (10 > 8), the floor's is far, so the shell goes through even though the
    // box's smallest half extent is 5; walls of 5 leave none.
    let mut api = Smlib::default();
    let slab = bx(&mut api, [0.0; 3], [20.0, 10.0, 30.0]);
    let top = [[10.0, 5.0, 30.0]];
    let r = api.shell(&slab, &top, 0.1, -4.0).unwrap();
    assert!(sound(&mut api, &r));
    let v = volume(&mut api, &r);
    assert!(
        (v - (6000.0 - 12.0 * 2.0 * 26.0)).abs() < 1e-6,
        "volume {v}"
    );
    assert_eq!(
        api.shell(&slab, &top, 0.1, -5.0).err(),
        Some(Status::ERR_NO_ROOM)
    );
}

#[test]
fn shell_arguments_and_picks() {
    let mut api = Smlib::default();
    let cube = bx(&mut api, [0.0; 3], [20.0; 3]);
    assert_eq!(api.shell(&cube, &[], 0.1, 0.0).err(), Some(Status::ERR_ARG));
    assert_eq!(
        api.shell(&cube, &[], 0.0, -1.0).err(),
        Some(Status::ERR_ARG)
    );
    assert_eq!(
        api.shell(&cube, &[], 0.1, f64::NAN).err(),
        Some(Status::ERR_ARG)
    );
    // the second point is far from every face
    let picks = [[10.0, 10.0, 20.0], [100.0, 100.0, 100.0]];
    assert_eq!(
        api.shell(&cube, &picks, 0.1, -2.0).err(),
        Some(Status::ERR_PICK_NONE)
    );
    assert_eq!(api.last_failed_point(), Some(1));
}

#[test]
fn draft_tilts_a_side_about_its_hinge() {
    // The cube's side x = 20 tilted 5 degrees about its line on z = 0, leaning in going up
    // (+z): a wedge of 20 x 20 tan(5 deg) / 2 x 20 comes off; -5 degrees adds it. The top
    // isn't parallel to the pull: not supported. The rod's wall (r 4 from z = -1 to 21) is
    // drafted about its circle on z = 0 into a frustum.
    let mut api = Smlib::default();
    let cube = bx(&mut api, [0.0; 3], [20.0; 3]);
    let rod = api.make_cylinder([10.0, 10.0, -1.0], 4.0, 22.0).unwrap();
    let side = [[20.0, 10.0, 10.0]];
    let origin = [0.0; 3];
    let up = [0.0, 0.0, 1.0];
    let deg = std::f64::consts::PI / 180.0;
    let wedge = 0.5 * 20.0 * 20.0 * (5.0 * deg).tan() * 20.0;
    let r = api
        .draft(&cube, &side, 0.1, origin, up, up, 5.0 * deg)
        .unwrap();
    let v = volume(&mut api, &r);
    assert!((v - (8000.0 - wedge)).abs() < 1e-6, "volume {v}");
    let r = api
        .draft(&cube, &side, 0.1, origin, up, up, -5.0 * deg)
        .unwrap();
    let v = volume(&mut api, &r);
    assert!((v - (8000.0 + wedge)).abs() < 1e-6, "volume {v}");
    let top = [[10.0, 10.0, 20.0]];
    assert_eq!(
        api.draft(&cube, &top, 0.1, origin, up, up, 5.0 * deg).err(),
        Some(Status::ERR_UNSUPPORTED)
    );
    let wall = [[14.0, 10.0, 5.0]];
    let r = api
        .draft(&rod, &wall, 0.1, origin, up, up, 5.0 * deg)
        .unwrap();
    let (r0, r1) = (4.0 + (5.0 * deg).tan(), 4.0 - 21.0 * (5.0 * deg).tan());
    let want = std::f64::consts::PI * 22.0 / 3.0 * (r0 * r0 + r0 * r1 + r1 * r1);
    let v = volume(&mut api, &r);
    assert!((v - want).abs() < 1e-6, "rod volume {v}, want {want}");
    let flat = [1.0, 0.0, 0.0];
    assert_eq!(
        api.draft(&cube, &side, 0.1, origin, flat, up, 5.0 * deg)
            .err(),
        Some(Status::ERR_ARG)
    );
    assert_eq!(
        api.draft(&cube, &side, 0.1, origin, up, up, 0.0).err(),
        Some(Status::ERR_ARG)
    );
    assert_eq!(
        api.draft(&cube, &side, 0.1, origin, up, up, 90.0 * deg)
            .err(),
        Some(Status::ERR_ARG)
    );
}

fn frustum(h: f64, r0: f64, r1: f64) -> f64 {
    std::f64::consts::PI * h / 3.0 * (r0 * r0 + r0 * r1 + r1 * r1)
}

fn near(got: f64, want: f64, rel: f64) -> bool {
    let ok = (got - want).abs() <= rel * want.abs();
    if !ok {
        eprintln!("  got {got:.12e}, want {want:.12e}");
    }
    ok
}

#[test]
fn drafted_cylinders_are_cones() {
    // A pin r 10, 20 tall, its side drafted 5 degrees about its base circle (neutral z = 0,
    // pull +z): the frustum from r 10 to 10 - 20 tan 5. A 30 x 30 x 10 block's through hole
    // r 5 drafted the same way widens to 5 + 10 tan 5 at the top. A cone r 10 at its base to
    // 8 at its top, 20 tall, drafted the same way turns its generators 5 degrees more:
    // r 10 + 20 tan(atan(-0.1) - 5 deg) at the top.
    let mut api = Smlib::default();
    let neutral = [0.0; 3];
    let up = [0.0, 0.0, 1.0];
    let a = 5.0 * std::f64::consts::PI / 180.0;
    let pin = api.make_cylinder([0.0; 3], 10.0, 20.0).unwrap();
    let r = api
        .draft(&pin, &[[0.0, 10.0, 10.0]], 0.1, neutral, up, up, a)
        .unwrap();
    assert!(sound(&mut api, &r));
    let v = volume(&mut api, &r);
    assert!(
        near(v, frustum(20.0, 10.0, 10.0 - 20.0 * a.tan()), 1e-7),
        "pin"
    );
    let block = bx(&mut api, [-15.0, -15.0, 0.0], [15.0, 15.0, 10.0]);
    let bore = api.make_cylinder([0.0, 0.0, -1.0], 5.0, 12.0).unwrap();
    let holed = api.boolean(&block, &bore, BOOLEAN_SUBTRACT).unwrap();
    let r = api
        .draft(&holed, &[[0.0, 5.0, 5.0]], 0.1, neutral, up, up, a)
        .unwrap();
    assert!(sound(&mut api, &r));
    let v = volume(&mut api, &r);
    assert!(
        near(v, 9000.0 - frustum(10.0, 5.0, 5.0 + 10.0 * a.tan()), 1e-7),
        "hole"
    );
    let cone = api.make_cone([0.0; 3], 10.0, 8.0, 20.0).unwrap();
    let r = api
        .draft(&cone, &[[0.0, 9.0, 10.0]], 0.1, neutral, up, up, a)
        .unwrap();
    assert!(sound(&mut api, &r));
    let top = 10.0 + 20.0 * ((-0.1f64).atan() - a).tan();
    let v = volume(&mut api, &r);
    assert!(near(v, frustum(20.0, 10.0, top), 1e-7), "cone");
}

#[test]
fn thicken_refuses_solids_and_bad_thickness() {
    let mut api = Smlib::default();
    let cube = bx(&mut api, [0.0; 3], [20.0; 3]);
    assert_eq!(api.thicken(&cube, 1.0, false).err(), Some(Status::ERR_ARG));
    assert_eq!(
        api.thicken(&cube, f64::NAN, false).err(),
        Some(Status::ERR_ARG)
    );
}

#[test]
fn thicken_an_l_sheet_into_a_prism() {
    // Rig case shell/5326: an L-shaped sheet, walls 30 and 20 long at a right angle and 20
    // deep, thickened 2 into its inside corner and 2 out of it: an L-shaped prism of 8 faces,
    // (30 + 20 -+ 2) * 2 * 20, one way each (which way is in depends on the sheet's normals,
    // here the linear sweep's). A flat sheet of two faces side by side thickened 1 keeps both
    // faces and their offsets apart: 8 faces, 20 * 20 * 1.
    let mut api = Smlib::default();
    let ell = sheet(
        &mut api,
        &[[0.0, 0.0, 0.0], [30.0, 0.0, 0.0], [30.0, 20.0, 0.0]],
    );
    let split = sheet(
        &mut api,
        &[[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [20.0, 0.0, 0.0]],
    );
    let mut ell_volumes = Vec::new();
    for t in [-2.0, 2.0] {
        let r = api.thicken(&ell, t, false).unwrap();
        let [faces, _, _] = api.brep_counts(&r).unwrap();
        let v = api.brep_volume(&r, 1e-10).unwrap();
        assert_eq!(faces, 8, "L thickened {t}");
        ell_volumes.push(v);
        let r = api.thicken(&split, t / 2.0, false).unwrap();
        let [faces, _, _] = api.brep_counts(&r).unwrap();
        let v = api.brep_volume(&r, 1e-10).unwrap();
        assert_eq!(faces, 8, "split thickened {}", t / 2.0);
        assert!((v - 400.0).abs() <= 1e-12 * 400.0, "split: volume {v}");
    }
    ell_volumes.sort_by(f64::total_cmp);
    for (v, want) in ell_volumes
        .into_iter()
        .zip([48.0 * 2.0 * 20.0, 52.0 * 2.0 * 20.0])
    {
        assert!(
            (v - want).abs() <= 1e-12 * want,
            "L thickened: volume {v}, want {want}"
        );
    }
    // symmetric: the same prism, half on each side
    let r = api.thicken(&split, 1.0, true).unwrap();
    let v = api.brep_volume(&r, 1e-10).unwrap();
    assert!((v - 400.0).abs() <= 1e-9 * 400.0, "symmetric: volume {v}");
}

#[test]
fn offset_sheet_moves_a_sheet() {
    let mut api = Smlib::default();
    let flat = sheet(&mut api, &[[0.0, 0.0, 0.0], [20.0, 0.0, 0.0]]);
    let r = api.offset_sheet(&flat, 1.5).unwrap();
    assert_eq!(api.brep_counts(&r).unwrap()[0], 1);
    let (lo, hi) = api.brep_bounds(&r).unwrap();
    // the face y = 0 moved 1.5 along its normal (+y or -y)
    assert!((lo[1].abs() - 1.5).abs() < 1e-9 && (hi[1] - lo[1]).abs() < 1e-9);
    assert_eq!(api.offset_sheet(&flat, 0.0).err(), Some(Status::ERR_ARG));
}

#[test]
fn offset_sheet_of_a_sphere_face() {
    // A sphere's one face (radius 10) moved out by 2 is the sphere of radius 12, a closed
    // solid; moved in to its centre, refused.
    let mut api = Smlib::default();
    let ball = api.make_sphere([0.0; 3], 10.0).unwrap();
    let r = api.offset_sheet(&ball, 2.0).unwrap();
    assert!(api.brep_is_manifold_solid(&r).unwrap());
    let v = volume(&mut api, &r);
    assert!(
        near(v, 4.0 / 3.0 * std::f64::consts::PI * 1728.0, 1e-6),
        "sphere 12"
    );
    assert!(api.offset_sheet(&ball, -10.0).is_err());
}

#[test]
fn offset_sheet_of_a_closed_sheet_of_voids() {
    // SolidCraft surface.offset of a ball's one face: the face taken off as a sheet is a
    // closed shell whose regions are all voids. Offset 2 out it is the sphere of radius 22,
    // area 4 pi 22^2.
    let mut api = Smlib::default();
    let ball = api.make_sphere([0.0; 3], 20.0).unwrap();
    let mut sheet = api.brep_copy(&ball).unwrap();
    for r in solidcore::topology::brep::p1::regions(&sheet, None).unwrap() {
        sheet.region_mut(r).unwrap().set_is_void(1);
    }
    let moved = api.offset_sheet(&sheet, 2.0).unwrap();
    let area = api.brep_area(&moved, 1e-8).unwrap();
    assert!(
        near(area, 4.0 * std::f64::consts::PI * 22.0 * 22.0, 1e-7),
        "area"
    );
}

#[test]
fn merge_faces_joins_coplanar_faces() {
    // Two boxes side by side: their union's top, bottom and two sides may come in two
    // faces each; merged, no more faces, the same solid.
    let mut api = Smlib::default();
    let a = bx(&mut api, [0.0; 3], [10.0; 3]);
    let b = bx(&mut api, [10.0, 0.0, 0.0], [20.0, 10.0, 10.0]);
    let joined = api.boolean(&a, &b, BOOLEAN_UNION).unwrap();
    let before = api.brep_counts(&joined).unwrap()[0];
    let r = api.merge_faces(&joined, &[], 0.1).unwrap();
    assert!(sound(&mut api, &r));
    let faces = api.brep_counts(&r).unwrap()[0];
    assert!(
        faces <= before && faces >= 6,
        "{before} faces merged to {faces}"
    );
    assert!(api.brep_is_manifold_solid(&r).unwrap());
    let v = volume(&mut api, &r);
    assert!((v - 2000.0).abs() < 1e-6, "volume {v}");
    // a point picking no edge
    assert_eq!(
        api.merge_faces(&joined, &[[100.0, 0.0, 0.0]], 0.1).err(),
        Some(Status::ERR_PICK_NONE)
    );
}
