//! Extend Sheet: a square plane sheet continued past one side in each mode, a bilinear
//! B-spline sheet continued naturally, and the shapes and arguments it refuses.
//!
//! Source: tangency-labs/smlib `capi/tests/extend.cpp` (Apache-2.0, Copyright SolidCraft
//! contributors) at fork commit b3dd0e5 (replay/1.0.1). The null-pointer check has no
//! counterpart. Its parabolic and cylinder-strip cases need Extrude (F3) and come with it.

#![allow(clippy::unwrap_used)] // test helpers

use solidcore::geom::object::ObjectRef;
use solidcore::geom::tol::USE_DEFAULT;
use solidcore::geom::vector3d::Point3;
use solidcore::primitive::api_primitives::create_plane;
use solidcore::primitive::api_surface_ops::create_surface_from_corner_points;
use solidcore::surface::surface::Surface;
use solidcore::topology::brep::p2::ArenaNotifyState;
use solidcore::topology::brep::p3::new_brep;
use solidcore::topology::brep::p8::{create_face_from_surface, FaceBuildEnv};
use solidcore::topology::topology_traverser::Marks;
use solidcore::topology::FaceId;
use solidcore::EntityRef;
use solidcore_smlib_api::{Brep, Smlib, Status, EXTEND_CURVATURE, EXTEND_NATURAL, EXTEND_TANGENT};

fn area_of(api: &mut Smlib, b: &Brep) -> f64 {
    api.brep_area(b, 1e-7).unwrap_or(-1.0)
}

/// A one-face sheet on `surface` over its natural domain (`CreateFaceFromSurface`).
fn sheet_of(api: &mut Smlib, surface: Surface) -> Brep {
    let k = api.kernel_mut();
    let mut b = new_brep(k, USE_DEFAULT, USE_DEFAULT).unwrap();
    let domain = surface.natural_uv_domain_with_nurb(k).unwrap();
    let mut state = ArenaNotifyState::new(&b);
    let mut env = state.env(None);
    let key_of: &dyn Fn(EntityRef) -> ObjectRef = env.key_of;
    let this_face = |f: FaceId| key_of(EntityRef::Face(f));
    let mut marks = Marks::default();
    let mut fenv = FaceBuildEnv {
        env: &mut env,
        marks: &mut marks,
        this_face: &this_face,
    };
    create_face_from_surface(k, &mut b, Some(surface), &domain, &mut fenv).unwrap();
    b
}

#[test]
fn extend_sheet() {
    let mut api = Smlib::default();

    // The square [0, 10] x [0, 10] at z = 0, continued 5 past its side x = 10: 15 x 10 in
    // each mode (a plane's G1 and G2 run-offs are the plane).
    let square = create_plane(api.kernel_mut(), Point3::new(0.0, 0.0, 0.0), 10.0, 10.0).unwrap();
    let side = [10.0, 5.0, 0.0];
    for mode in [EXTEND_NATURAL, EXTEND_TANGENT, EXTEND_CURVATURE] {
        let r = api.extend_sheet(&square, side, 1e-3, 5.0, mode).unwrap();
        assert!((area_of(&mut api, &r) - 150.0).abs() < 1e-6, "mode {mode}");
        let n = api.brep_counts(&r).unwrap();
        assert_eq!((n[0], n[1]), (1, 4), "mode {mode}");
    }
    // The input is unchanged.
    assert!((area_of(&mut api, &square) - 100.0).abs() < 1e-6);

    // A bilinear B-spline sheet (a flat square 10 x 10) continued 4 past its side y = 10.
    let s = create_surface_from_corner_points(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(10.0, 0.0, 0.0),
        Point3::new(0.0, 10.0, 0.0),
        Point3::new(10.0, 10.0, 0.0),
    )
    .unwrap();
    let spline = sheet_of(&mut api, Surface::from(s));
    let grown = api
        .extend_sheet(&spline, [5.0, 10.0, 0.0], 1e-3, 4.0, EXTEND_NATURAL)
        .unwrap();
    assert!((area_of(&mut api, &grown) - 140.0).abs() < 1e-4);

    // Refused: a box (not one face), a pick on no edge, bad arguments.
    let bx = api.make_box([0.0; 3], [10.0; 3]).unwrap();
    assert_eq!(
        api.extend_sheet(&bx, side, 1e-3, 5.0, EXTEND_NATURAL)
            .unwrap_err(),
        Status::ERR_UNSUPPORTED
    );
    assert_eq!(
        api.extend_sheet(&square, [5.0, 5.0, 0.0], 1e-3, 5.0, EXTEND_NATURAL)
            .unwrap_err(),
        Status::ERR_PICK_NONE
    );
    assert_eq!(
        api.extend_sheet(&square, side, 1e-3, -1.0, EXTEND_NATURAL)
            .unwrap_err(),
        Status::ERR_ARG
    );
    assert_eq!(
        api.extend_sheet(&square, side, 1e-3, 5.0, 7).unwrap_err(),
        Status::ERR_ARG
    );
}
