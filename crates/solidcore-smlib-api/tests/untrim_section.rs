//! Untrim Face and Face Section on a drilled box: the top face's hole taken away, the
//! face in a body of its own; and the top face's section by a plane as a NURBS curve.
//!
//! The fork's capi tests for these functions are in its `ops_kernel.cpp` kernel suite;
//! these cases follow `capi/include/smlib.h`'s contracts (Apache-2.0, Copyright SolidCraft
//! contributors) at fork commit 183ea65 (replay/1.0.1).

#![allow(clippy::unwrap_used)] // test helpers

use solidcore::geom::attribute::{Attribute, AttributeValue, Behavior};
use solidcore::geom::vector2d::Point2;
use solidcore::topology::face::p2 as face2;
use solidcore::topology::FaceId;
use solidcore::EntityRef;
use solidcore_smlib_api::{Brep, Smlib, Status};

/// A 10 x 10 x 10 box with an r 2 hole through it along z at (5, 5).
fn drilled(api: &mut Smlib) -> Brep {
    let bx = api.make_box([0.0; 3], [10.0; 3]).unwrap();
    let pin = api.make_cylinder([5.0, 5.0, -1.0], 2.0, 12.0).unwrap();
    api.boolean(&bx, &pin, solidcore_smlib_api::BOOLEAN_SUBTRACT)
        .unwrap()
}

#[test]
fn untrim_takes_the_hole_away() {
    let mut api = Smlib::default();
    let b = drilled(&mut api);
    let top = api
        .untrim_face(&b, [1.0, 1.0, 10.0], 1e-3, None, false, false)
        .unwrap();
    let n = api.brep_counts(&top).unwrap();
    assert_eq!(n[0], 1);
    assert!((api.brep_area(&top, 1e-7).unwrap() - 100.0).abs() < 1e-6);
    // one hole by its point
    let top = api
        .untrim_face(
            &b,
            [1.0, 1.0, 10.0],
            1e-3,
            Some(&[[7.0, 5.0, 10.0]]),
            false,
            false,
        )
        .unwrap();
    assert!((api.brep_area(&top, 1e-7).unwrap() - 100.0).abs() < 1e-6);
    // nothing near the pick, or nothing to do
    assert_eq!(
        api.untrim_face(&b, [50.0, 50.0, 50.0], 1e-3, None, false, false)
            .unwrap_err(),
        Status::ERR_PICK_NONE
    );
    assert_eq!(
        api.untrim_face(&b, [1.0, 1.0, 10.0], 1e-3, Some(&[]), false, false)
            .unwrap_err(),
        Status::ERR_ARG
    );
    assert_eq!(
        api.untrim_face(&b, [1.0, 1.0, 10.0], 0.0, None, false, false)
            .unwrap_err(),
        Status::ERR_ARG
    );
}

/// Untrim into a body of its own (`in_place` false) gives the new face the picked face's
/// attributes (`MakeCopy`): each face of the drilled box tagged with its index, the result's
/// one face carries the top face's tag, and only it.
#[test]
fn untrim_keeps_the_face_attributes() {
    const TAG: u64 = 10_001;
    let mut api = Smlib::default();
    let mut b = drilled(&mut api);
    let faces: Vec<FaceId> = b.faces().ids().collect();
    for (i, &f) in faces.iter().enumerate() {
        let a = Attribute::long(TAG, i64::try_from(i).unwrap(), Behavior::Copy);
        b.add_attribute(EntityRef::Face(f), a).unwrap();
    }
    // The top face: its surface at z = 10 at its domain's corners and middle.
    let k = api.kernel();
    let top: Vec<i64> = (0..faces.len())
        .filter(|&i| {
            let s = face2::surface(&b, faces[i]).unwrap().unwrap();
            let Ok(d) = s.natural_uv_domain_with_nurb(k) else {
                return false;
            };
            let (lo, hi) = (d.min(), d.max());
            let mid = Point2::new(0.5 * (lo.x + hi.x), 0.5 * (lo.y + hi.y));
            [lo, hi, mid].iter().all(|&uv| {
                s.evaluate_point(k, uv)
                    .is_ok_and(|p| (p.z - 10.0).abs() < 1e-9)
            })
        })
        .map(|i| i64::try_from(i).unwrap())
        .collect();
    assert_eq!(top.len(), 1);
    let r = api
        .untrim_face(&b, [1.0, 1.0, 10.0], 1e-3, None, false, false)
        .unwrap();
    let got: Vec<FaceId> = r.faces().ids().collect();
    assert_eq!(got.len(), 1);
    let tags: Vec<i64> = r
        .entity_attributes(EntityRef::Face(got[0]))
        .into_iter()
        .filter(|(_, a)| a.attribute_id() == TAG)
        .filter_map(|(_, a)| match a.value() {
            AttributeValue::Long(v) => Some(*v),
            _ => None,
        })
        .collect();
    assert_eq!(tags, top);
}

#[test]
fn face_section_of_the_top() {
    let mut api = Smlib::default();
    let bx = api.make_box([0.0; 3], [10.0; 3]).unwrap();
    // the plane x = 3 through the top face: one line, from y = 0 to 10 at z = 10
    let c = api
        .face_section(
            &bx,
            [5.0, 5.0, 10.0],
            1e-3,
            [3.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            1e-6,
        )
        .unwrap();
    let (_, _, n) = c.size();
    assert_eq!(n, 1);
    let kinds = c.kinds();
    assert_eq!(kinds[0], 3); // SMLIB_CURVE_NURBS
    let count = usize::try_from(kinds[2]).unwrap();
    for p in c.reals()[..4 * count].chunks(4) {
        assert!(
            (p[0] - 3.0).abs() < 1e-9 && (p[2] - 10.0).abs() < 1e-9,
            "{p:?}"
        );
    }
    // a plane that misses the face, a bad normal
    assert_eq!(
        api.face_section(
            &bx,
            [5.0, 5.0, 10.0],
            1e-3,
            [30.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            1e-6
        )
        .unwrap_err(),
        Status::ERR_EMPTY
    );
    assert_eq!(
        api.face_section(&bx, [5.0, 5.0, 10.0], 1e-3, [3.0, 0.0, 0.0], [0.0; 3], 1e-6)
            .unwrap_err(),
        Status::ERR_ARG
    );
}
