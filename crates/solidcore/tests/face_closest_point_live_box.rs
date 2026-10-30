//! `face_closest_point_in_shape` (`SmApiFaceClosestPoint`) over a bare `LiveBrepObjects`,
//! as a library caller makes it (no B-rep cache, no `face_entity`), on a `create_box`
//! B-rep as made (not renumbered or renamed), its geometry named by the B-rep's own
//! owner identities (`Brep::owner_identity`). At the centre of each of the box's faces
//! the face's closest point is the centre itself, at distance 0 (the SolidCraft fork's
//! C API picks faces by their centres): plan/findings
//! `query.api_queries--face-closest-point-solve-answers-a-vertex.md`.
//!
//! The contract the solve relies on: the B-rep's curves and surfaces are owned by its
//! edges and faces (the solve finds a surface's face, a curve's edge, through the
//! owner). `Brep::clone` resets every owner (the geometry's copy constructor; see
//! `impl Clone for Brep`), so on a clone no face or edge answers, only the vertices; a
//! copy made by `api_brep_copy` (`SmApiBrepCopy`) owns its geometry and answers as the
//! original.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers

use solidcore::geom::object::{ObjectRef, Shape};
use solidcore::math::Point3;
use solidcore::primitive::api_primitives::create_box;
use solidcore::query::api_queries::face_closest_point_in_shape;
use solidcore::query::brep::{BrepIdentity, LiveBrepObjects};
use solidcore::topology::brep::p1;
use solidcore::topology::brep::p2::ArenaNotifyState;
use solidcore::topology::brep::p3::api_brep_copy;
use solidcore::topology::brep_cache::{arena_identity, CacheObject};
use solidcore::topology::{Brep, BrepId, FaceId};
use solidcore::{EntityRef, Kernel};

/// The closest point of face `f` of `b` to `p` and its distance, through a bare live view
/// whose keys are the B-rep's own owner identities.
fn closest(k: &mut Kernel, b: &mut Brep, f: FaceId, p: Point3) -> (Point3, f64) {
    let id = BrepId::from_index(0);
    let shape = Shape::new(k, &(id, &*b), &[], &[(id, f)], &[], &[]).unwrap();
    let names = b.clone();
    let keys = move |o: CacheObject| match o {
        CacheObject::Face(f) => names.owner_identity(EntityRef::Face(f)),
        CacheObject::Edge(e) => names.owner_identity(EntityRef::Edge(e)),
        CacheObject::Vertex(v) => names.owner_identity(EntityRef::Vertex(v)),
        other => arena_identity(other),
    };
    let owner_of = |_: ObjectRef| -> Option<ObjectRef> { None };
    let identity = BrepIdentity {
        brep_id: id,
        keys: &keys,
        owner_of: &owner_of,
    };
    let objects = LiveBrepObjects::new(b, &identity, None).unwrap();
    let c = face_closest_point_in_shape(k, &shape, &objects, (id, f), p).unwrap();
    (c.point, c.distance)
}

#[test]
fn box_face_centres_are_their_own_closest_points() {
    let mut k = Kernel::default();
    let mut b = create_box(&mut k, Point3::new(0.0, 0.0, 0.0), 10.0, 20.0, 30.0).unwrap();
    let centres = [
        Point3::new(5.0, 10.0, 0.0),
        Point3::new(5.0, 10.0, 30.0),
        Point3::new(0.0, 10.0, 15.0),
        Point3::new(10.0, 10.0, 15.0),
        Point3::new(5.0, 0.0, 15.0),
        Point3::new(5.0, 20.0, 15.0),
    ];
    let faces = p1::faces(&b, None).unwrap();
    for c in centres {
        let mut answers = Vec::new();
        for &f in &faces {
            let (q, d) = closest(&mut k, &mut b, f, c);
            answers.push((f, q, d));
        }
        eprintln!("centre {c:?}: {answers:?}");
        let on = answers
            .iter()
            .filter(|(_, q, d)| *d <= 1e-9 && (*q - c).length() <= 1e-9)
            .count();
        assert_eq!(on, 1, "centre {c:?}: exactly its own face at distance 0");
        // every other face's answer is no farther than its nearest corner
        for (_, _, d) in answers {
            assert!(d <= ((5.0f64).powi(2) + 10.0f64.powi(2) + 30.0f64.powi(2)).sqrt());
        }
    }
}

/// The closest point of face `f` of `b` to `p`, through a bare live view whose keys
/// name each face and edge by its own geometry's owner (as a caller identifies objects
/// by pointer), or by its arena identity where its geometry has no owner.
fn closest_by_owner(k: &mut Kernel, b: &mut Brep, f: FaceId, p: Point3) -> (Point3, f64) {
    let id = BrepId::from_index(0);
    let shape = Shape::new(k, &(id, &*b), &[], &[(id, f)], &[], &[]).unwrap();
    // (read before the live view borrows the B-rep)
    let mut owners: Vec<(CacheObject, ObjectRef)> = Vec::new();
    for (f, face) in b.faces().iter() {
        if let Some(o) = face.surface.as_ref().and_then(|s| s.owner()) {
            owners.push((CacheObject::Face(f), o));
        }
    }
    for (e, edge) in b.edges().iter() {
        if let Some(o) = edge.curve.as_ref().and_then(|c| c.owner()) {
            owners.push((CacheObject::Edge(e), o));
        }
    }
    let keys = move |o: CacheObject| {
        owners
            .iter()
            .find(|(x, _)| *x == o)
            .map_or_else(|| arena_identity(o), |&(_, r)| r)
    };
    let owner_of = |_: ObjectRef| -> Option<ObjectRef> { None };
    let identity = BrepIdentity {
        brep_id: id,
        keys: &keys,
        owner_of: &owner_of,
    };
    let objects = LiveBrepObjects::new(b, &identity, None).unwrap();
    let c = face_closest_point_in_shape(k, &shape, &objects, (id, f), p).unwrap();
    (c.point, c.distance)
}

/// How many of the box's six face centres some face of `b` answers at distance 0, at
/// the centre itself.
fn centres_answered(
    k: &mut Kernel,
    b: &mut Brep,
    solve: fn(&mut Kernel, &mut Brep, FaceId, Point3) -> (Point3, f64),
) -> usize {
    let centres = [
        Point3::new(5.0, 10.0, 0.0),
        Point3::new(5.0, 10.0, 30.0),
        Point3::new(0.0, 10.0, 15.0),
        Point3::new(10.0, 10.0, 15.0),
        Point3::new(5.0, 0.0, 15.0),
        Point3::new(5.0, 20.0, 15.0),
    ];
    let faces = p1::faces(b, None).unwrap();
    centres
        .into_iter()
        .filter(|&c| {
            faces.iter().any(|&f| {
                let (q, d) = solve(k, b, f, c);
                d <= 1e-9 && (q - c).length() <= 1e-9
            })
        })
        .count()
}

#[test]
fn a_clone_answers_no_face_and_an_api_copy_answers_every_centre() {
    let mut k = Kernel::default();
    let mut b = create_box(&mut k, Point3::new(0.0, 0.0, 0.0), 10.0, 20.0, 30.0).unwrap();
    assert_eq!(centres_answered(&mut k, &mut b, closest_by_owner), 6);
    // a clone's geometry has no owner: no face (nor edge) answers, only vertices
    let mut clone = b.clone();
    assert_eq!(centres_answered(&mut k, &mut clone, closest_by_owner), 0);
    // SmApiBrepCopy: the copy owns its geometry
    let mut state = ArenaNotifyState::new(&Brep::new());
    let mut env = state.env(None);
    let mut copy = api_brep_copy(&mut k, &mut env, Some(&b), &|e| b.owner_identity(e)).unwrap();
    assert_eq!(centres_answered(&mut k, &mut copy, closest_by_owner), 6);
}
