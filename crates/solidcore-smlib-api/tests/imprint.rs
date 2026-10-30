//! Imprint: a box cut by a plane sheet through its middle, a sheet cut by a solid that
//! crosses it, a tool that misses, rejected arguments, and a point inside each face.
//!
//! Source: tangency-labs/smlib `capi/tests/imprint.cpp` (Apache-2.0, Copyright SolidCraft
//! contributors) at fork commit 30db2d7 (replay/1.0.1). The null-pointer checks have no counterpart.

#![allow(clippy::unwrap_used)] // test helpers

use solidcore::geom::vector3d::Point3;
use solidcore::primitive::api_primitives::create_plane;
use solidcore_smlib_api::{Brep, Smlib, Status};

/// A `size` x `size` plane sheet at height `z`, its corner at `(x0, x0)`
/// (`SmApiCreatePlane` then `SmApiTranslate`).
fn plane_at_z(api: &mut Smlib, z: f64, x0: f64, size: f64) -> Brep {
    create_plane(api.kernel_mut(), Point3::new(x0, x0, z), size, size).unwrap()
}

#[test]
fn imprint() {
    let mut api = Smlib::default();
    let bx = api.make_box([0.0; 3], [10.0; 3]).unwrap();

    // A plane through the middle cuts the four side faces: 6 + 4 faces, the volume unchanged.
    let mid = plane_at_z(&mut api, 5.0, -5.0, 20.0);
    let cut = api.imprint(&bx, &mid).unwrap();
    assert_eq!(api.brep_counts(&cut).unwrap()[0], 10);
    assert!((api.brep_volume(&cut, 1e-7).unwrap() - 1000.0).abs() < 1e-6);
    // The inputs are unchanged.
    assert_eq!(api.brep_counts(&bx).unwrap()[0], 6);
    assert_eq!(api.brep_counts(&mid).unwrap()[0], 1);

    // The sheet cut by the box: its one face splits where the box's walls cross it.
    let sheet_cut = api.imprint(&mid, &bx).unwrap();
    assert!(api.brep_counts(&sheet_cut).unwrap()[0] >= 2);
    assert!((api.brep_area(&sheet_cut, 1e-7).unwrap() - 400.0).abs() < 1e-6);

    // A plane that misses the box cuts nothing.
    let far = plane_at_z(&mut api, 50.0, -5.0, 20.0);
    assert_eq!(api.imprint(&bx, &far).unwrap_err(), Status::ERR_NO_CUT);
    assert_eq!(api.imprint(&bx, &bx).unwrap_err(), Status::ERR_ARG);

    // A point inside each face of the cut box, on its box.
    let points = api.face_points(&cut).unwrap();
    assert_eq!(points.len(), 10);
    for p in points {
        assert!(p.iter().all(|x| (-1e-9..=10.0 + 1e-9).contains(x)));
    }
}
