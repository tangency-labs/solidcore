//! #2503: a B-rep ray classification keeps the face surfaces' trimmed caches, as
//! upstream's `SmBrep::Point3DClassify` does through `smsurf_GetSurfaceCache`
//! (`SmSurface::CacheMakeOrValidate`, SmSurface.cpp:5784). A cache kept on a face's
//! surface is reused while the B-rep is being edited, even when its face was modified;
//! `SmTrimSrfCache::BuildTree` implants no trims while editing is on, so a cache first
//! built then would count every face as its untrimmed rectangle. Before #2503 every ray
//! built its caches afresh, and once editing was on (the fillet session's trim step,
//! #2365) points inside the bar classified as outside.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)] // test helpers

use solidcore::boolean::api_brep::boolean_union;
use solidcore::intersect::point_class::PointClassType;
use solidcore::math::{Point3, Vec3};
use solidcore::primitive::api_primitives::{create_box, create_cylinder};
use solidcore::query::brep::BrepQuery;
use solidcore::surface::surface::Surface;
use solidcore::topology::brep::p2::delete_stale_caches;
use solidcore::topology::proto_topology::p1::api_transform::{api_rotate, TransformObject};
use solidcore::topology::FaceId;
use solidcore::{Brep, EntityRef, Kernel};

/// Points in the bar, in the post, and outside both.
const POINTS: [[f64; 3]; 6] = [
    [5.0, 4.0, 3.15],
    [20.0, 8.0, 5.0],
    [40.0, 8.0, 15.0],
    [44.0, 8.0, 5.0],
    [-1.0, 8.0, 5.0],
    [20.0, 8.0, 11.0],
];

/// The bar running into the post (#2365): a 40 x 16 x 10 bar joined to an R8 x 20
/// cylinder about (40, 8).
fn bar_into_post(k: &mut Kernel) -> Brep {
    let bar = create_box(k, Point3::new(0.0, 0.0, 0.0), 40.0, 16.0, 10.0).unwrap();
    let post = create_cylinder(k, Point3::new(40.0, 8.0, 0.0), 8.0, 20.0).unwrap();
    boolean_union(k, bar, post).unwrap()
}

/// Whether each point classifies to a material region of `b`.
fn in_material(k: &Kernel, b: &mut Brep) -> Vec<bool> {
    in_material_at(k, b, &POINTS)
}

/// Whether each of `points` classifies to a material region of `b`.
fn in_material_at(k: &Kernel, b: &mut Brep, points: &[[f64; 3]]) -> Vec<bool> {
    points
        .iter()
        .map(|&[x, y, z]| {
            let pc = BrepQuery::of_brep(b)
                .point_3d_classify(k, b, Point3::new(x, y, z), 1e-6, false)
                .unwrap();
            pc.point_class() == PointClassType::Region
                && match pc.object() {
                    Some(EntityRef::Region(r)) => b.region(r).unwrap().is_void() == 0,
                    _ => false,
                }
        })
        .collect()
}

const EXPECTED: [bool; 6] = [true, true, true, true, false, false];

#[test]
fn classification_is_unchanged_once_editing_is_on() {
    let mut k = Kernel::default();
    let mut b = bar_into_post(&mut k);
    assert_eq!(in_material(&k, &mut b), EXPECTED, "editing off");
    // the caches the rays kept are reused while the B-rep is edited
    b.set_editing_enabled(true);
    assert_eq!(in_material(&k, &mut b), EXPECTED, "editing on");
    b.set_editing_enabled(false);
    assert_eq!(in_material(&k, &mut b), EXPECTED, "editing off again");
}

/// The face surfaces holding a trimmed cache.
fn faces_with_caches(b: &Brep) -> Vec<FaceId> {
    b.faces()
        .iter()
        .filter(|(_, f)| f.surface.as_ref().is_some_and(Surface::has_trim_cache))
        .map(|(id, _)| id)
        .collect()
}

#[test]
fn the_rays_keep_the_face_surfaces_trimmed_caches() {
    let mut k = Kernel::default();
    let mut b = bar_into_post(&mut k);
    in_material(&k, &mut b);
    let kept = faces_with_caches(&b);
    assert!(!kept.is_empty(), "the rays kept no trimmed cache");
    // the next rays reuse them, while the B-rep is edited too
    in_material(&k, &mut b);
    b.set_editing_enabled(true);
    in_material(&k, &mut b);
    b.set_editing_enabled(false);
    assert_eq!(faces_with_caches(&b), kept);
    // a cache not marked modified is not stale
    delete_stale_caches(&b).unwrap();
    assert_eq!(faces_with_caches(&b), kept);
}

#[test]
fn delete_stale_caches_drops_the_caches_of_modified_faces() {
    let mut k = Kernel::default();
    let mut b = bar_into_post(&mut k);
    in_material(&k, &mut b);
    let with_cache = faces_with_caches(&b);
    let (&modified, unmodified) = with_cache.split_first().expect("a kept cache");
    // `SmFace::Notify(SM_NO_PRE_EDIT)` marks the face surface's cache modified
    b.face(modified)
        .unwrap()
        .surface
        .as_ref()
        .unwrap()
        .mark_trim_cache_face_modified();
    delete_stale_caches(&b).unwrap();
    let has = |id| {
        b.face(id)
            .unwrap()
            .surface
            .as_ref()
            .is_some_and(Surface::has_trim_cache)
    };
    assert!(!has(modified));
    assert!(unmodified.iter().all(|&id| has(id)));
}

/// `SmBSplineSurface::Transform` notifies `SM_NO_PRE_EDIT`, and `SmSurface::Notify`
/// deletes the surface's trimmed cache (`SmCacheMgr::DeleteObjectsCache(SM_OC_TRIMSRF)`),
/// so a body moved after its caches were built classifies where it now is, also while
/// it is edited (the Boolean of `boolean.api_brep`'s `arc_lengths_cylinder_cylinder_tee_union`
/// kept the moved cylinder where it was made).
#[test]
fn a_moved_body_drops_its_trimmed_caches() {
    let mut k = Kernel::default();
    // along +z from the origin; its stitch's infinite-region ray keeps the caches
    let mut b = create_cylinder(&mut k, Point3::new(0.0, 0.0, 0.0), 1.5, 10.0).unwrap();
    assert!(
        !faces_with_caches(&b).is_empty(),
        "the ray keeps the face surfaces' caches"
    );
    // a quarter turn about x: along -y
    api_rotate(
        &mut k,
        Some(TransformObject::Brep(&mut b)),
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        90.0,
    )
    .unwrap();
    assert_eq!(
        faces_with_caches(&b).len(),
        0,
        "a transformed surface keeps a trimmed cache"
    );
    let points = [[0.0, -5.0, 0.0], [0.0, 0.0, 5.0]];
    b.set_editing_enabled(true);
    assert_eq!(
        in_material_at(&k, &mut b, &points),
        [true, false],
        "editing on"
    );
    b.set_editing_enabled(false);
    assert_eq!(
        in_material_at(&k, &mut b, &points),
        [true, false],
        "editing off"
    );
}
