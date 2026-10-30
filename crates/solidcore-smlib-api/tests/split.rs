//! Split Face, Split Body and Stitch on a box: the top face cut by a plane, the box cut
//! in two, and its faces sewn back into a solid.
//!
//! The fork's capi tests for these functions are in its `ops_kernel.cpp` kernel suite;
//! these cases follow `capi/include/smlib.h`'s contracts (Apache-2.0, Copyright SolidCraft
//! contributors) at fork commit 183ea65 (replay/1.0.1).

#![allow(clippy::unwrap_used)] // test helpers

use solidcore_smlib_api::{Smlib, Status};

const X: [f64; 3] = [1.0, 0.0, 0.0];
const Y: [f64; 3] = [0.0, 1.0, 0.0];
const Z: [f64; 3] = [0.0, 0.0, 1.0];

#[test]
fn split_faces_cuts_the_chosen_face_only() {
    let mut api = Smlib::default();
    let bx = api.make_box([0.0; 3], [10.0; 3]).unwrap();
    // the plane x = 4 (spanned by y and z), reaching over the box
    let tool = api
        .surface_plane([4.0, 0.0, 0.0], Y, Z, [-20.0, 20.0], [-20.0, 20.0])
        .unwrap();
    let faces = api.face_points(&bx).unwrap();
    // the top face (its interior point at z = 10)
    let top = faces
        .iter()
        .position(|p| (p[2] - 10.0).abs() < 1e-9)
        .unwrap();
    let (cut, made) = api.split_faces(&bx, &[top], &[&tool]).unwrap();
    assert_eq!(made, vec![1]);
    assert_eq!(api.brep_counts(&cut).unwrap()[0], 7);
    assert!((api.brep_volume(&cut, 1e-7).unwrap() - 1000.0).abs() < 1e-6);
    // a face the tool misses: no cut
    let bottom = faces.iter().position(|p| p[2].abs() < 1e-9).unwrap();
    let far = api
        .surface_plane([40.0, 0.0, 0.0], Y, Z, [-20.0, 20.0], [-20.0, 20.0])
        .unwrap();
    assert_eq!(
        api.split_faces(&bx, &[bottom], &[&far]).unwrap_err(),
        Status::ERR_NO_CUT
    );
    // rejected arguments
    assert_eq!(
        api.split_faces(&bx, &[], &[&tool]).unwrap_err(),
        Status::ERR_ARG
    );
    assert_eq!(
        api.split_faces(&bx, &[top], &[]).unwrap_err(),
        Status::ERR_ARG
    );
    assert_eq!(
        api.split_faces(&bx, &[top, top], &[&tool]).unwrap_err(),
        Status::ERR_ARG
    );
    assert_eq!(
        api.split_faces(&bx, &[6], &[&tool]).unwrap_err(),
        Status::ERR_ARG
    );
}

#[test]
fn split_body_halves_a_box() {
    let mut api = Smlib::default();
    let bx = api.make_box([0.0; 3], [10.0; 3]).unwrap();
    // the plane z = 4, its normal +z
    let tool = api
        .surface_plane([0.0, 0.0, 4.0], X, Y, [-20.0, 20.0], [-20.0, 20.0])
        .unwrap();
    let (back, front) = api.split_body(&bx, &tool).unwrap();
    let (vb, vf) = (
        api.brep_volume(&back, 1e-7).unwrap(),
        api.brep_volume(&front, 1e-7).unwrap(),
    );
    // one side is 10 x 10 x 4, the other 10 x 10 x 6
    let mut v = [vb, vf];
    v.sort_by(f64::total_cmp);
    assert!(
        (v[0] - 400.0).abs() < 1e-6 && (v[1] - 600.0).abs() < 1e-6,
        "{v:?}"
    );
    assert!(
        api.brep_is_manifold_solid(&back).unwrap() && api.brep_is_manifold_solid(&front).unwrap()
    );
    // a plane beside the box misses it: each part is the whole body (smlib.h: "A tool that
    // misses the body leaves each part the whole body")
    let far = api
        .surface_plane([0.0, 0.0, 40.0], X, Y, [-20.0, 20.0], [-20.0, 20.0])
        .unwrap();
    let (back, front) = api.split_body(&bx, &far).unwrap();
    for part in [&back, &front] {
        assert!((api.brep_volume(part, 1e-7).unwrap() - 1000.0).abs() < 1e-6);
        assert_eq!(api.brep_counts(part).unwrap()[0], 6);
    }
}

#[test]
fn stitch_keeps_a_closed_box_closed() {
    let mut api = Smlib::default();
    let bx = api.make_box([0.0; 3], [10.0; 3]).unwrap();
    let (sewn, closed) = api.stitch(&bx, 1e-3).unwrap();
    assert!(closed);
    assert_eq!(api.brep_counts(&sewn).unwrap()[0], 6);
    assert!((api.brep_volume(&sewn, 1e-7).unwrap() - 1000.0).abs() < 1e-6);
    assert_eq!(api.stitch(&bx, 0.0).unwrap_err(), Status::ERR_ARG);
    assert_eq!(api.stitch(&bx, 2e3).unwrap_err(), Status::ERR_ARG);
}
