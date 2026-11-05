//! `SmBrep::LocalOperation`: replaces the surfaces of some faces of a B-rep by given
//! ones, extending or trimming the faces around them to meet the new surfaces.
//!
//! Upstream: `source/SMLib/src/SmBrep.cpp` (`SmBrep::LocalOperation`),
//! `source/SMLib/inc/SmBrep.h`.
//!
//! | Rust                | Upstream                  |
//! |---------------------|---------------------------|
//! | [`local_operation`] | `SmBrep::LocalOperation`  |
//! | [`local_operation_simple`] | (not upstream) `SmBrep::LocalOperation` with its own notification identities |
//!
//! Deviations: the offset executive's local-operation geometry creator
//! (`SmLocalOperationGeometryCreation`) is the face and surface lists the executive
//! holds ([`OffsetExecutive::set_local_surfaces`]). The executive's notifications are
//! the caller's ([`OffsetNotifications`]). Upstream removes every face of the B-rep and
//! copies the faces of the offset B-rep in with a merge (`RemoveFaces`, `CopyFaces`), as
//! here; the copy notifies through an arena-only store. (not upstream) The B-rep's own
//! attributes ([`Brep::add_attribute`]) are lent to the offset under the caller's
//! original identities, the offset B-rep takes what the offset leaves on its objects,
//! and the copy back gives them to the B-rep (SMLib's Copy handler), so each face keeps
//! its attributes as upstream's do. `SmTemporaryChangeValue` of
//! `m_bEditingEnabled` is a set and restore of the flag. The debug drawing code is omitted. No
//! recursion is introduced here. Opt-in (WR-534,
//! `Kernel::fixes.local_operation_box_holds_new_surfaces`, the fork's
//! tangency-labs/smlib@a1ee273): the box that sizes the executive's offset distance also
//! holds the replacement surfaces. (not upstream) The B-rep's cached topology fields are
//! refreshed before the operation returns ([`Brep::refresh_cached_topology`]), as every
//! public operation leaves them.

use std::collections::BTreeSet;

use crate::error::offset::OffsetFailure;
use crate::error::{EntityRef, Failure, Result, SolidError};
use crate::geom::attribute::ObjectKey;
use crate::kernel::Kernel;
use crate::math::Vec3;
use crate::nurbs::b_spline_surface::p1::SurfParamType;
use crate::nurbs::knots::Side;
use crate::offset::offset_executive::p1::{OffsetExecutive, OffsetOperation};
use crate::offset::offset_executive::p2::OffsetNotifications;
use crate::query::brep::{calculate_bounding_box, BoundingBoxOptions};
use crate::surface::surface::{Surface, SurfaceKind};
use crate::topology::brep::p1 as bp1;
use crate::topology::brep::p2;
use crate::topology::brep::p5;
use crate::topology::face::p2 as face_ops;
use crate::topology::{Brep, FaceId};

pub use crate::offset::offset_executive::p2::local_operation_simple;

const OP: &str = "topology::brep::local_operation::local_operation";

fn null() -> SolidError {
    SolidError::null_pointer(OP)
}

/// `SmBrep::LocalOperation(rFaces, rNewSurfaces)`: `surfaces[i]` becomes the surface of
/// `faces[i]` (turned about if its normal opposes the old one's), every other face
/// keeping a copy of its surface (as an analytic one when possible), and the faces are
/// rebuilt by an extended offset that extends the surfaces to intersect.
///
/// To find where `surfaces[i]` went, mark `faces[i]` before the call
/// ([`Brep::add_attribute`], [`Behavior::Copy`](crate::geom::attribute::Behavior::Copy))
/// and look for the mark on the result's faces: upstream callers mark the `SmSurface`
/// itself, but the port names a surface by its owning face, and `faces[i]`'s offset face
/// takes `surfaces[i]`, so the mark travels with it as upstream's does with the surface
/// (tangency-labs/solidcore#1875; the one difference, a merge keeping an unpicked
/// survivor, is `plan/findings/topology.brep_local_operation--1875-face-mark-for-surface-mark.md`).
///
/// # Errors
///
/// `SM_ERR` when the extended offset fails (its own status with
/// `Kernel::fixes.local_operation_offset_status`); `InvalidInput` where upstream dereferences
/// null (a surface list shorter than the face list, a face without surface);
/// [`OffsetFailure::SeamOpened`] (not upstream) when a face that meets itself across an
/// edge (a seam) is given a surface that is not closed, before anything changes; those of
/// the callees ([`crate::offset::offset_executive::p2`], the surfaces' normals,
/// reversal and analytic copy, the bounding box).
pub fn local_operation(
    k: &mut Kernel,
    b: &mut Brep,
    faces: &[FaceId],
    surfaces: Vec<Surface>,
    notifications: &mut OffsetNotifications<'_>,
) -> Result<()> {
    if surfaces.len() < faces.len() {
        return Err(null());
    }
    // Make a new list of all surfaces in this Brep: for every face, the new surface if it
    // is in the list to be replaced (turned to the original's side), else a copy of its own.
    let all_faces = bp1::faces(b, None)?;
    let mut given: Vec<Option<Surface>> = surfaces.into_iter().map(Some).collect();
    let mut new_list = Vec::with_capacity(all_faces.len());
    // (WR-534) where each given surface went in the new list
    let mut placed: Vec<Option<usize>> = vec![None; given.len()];
    // (not upstream) the faces whose surfaces move: a face given its own plane again keeps
    // a copy of its own surface, as the faces not given one do, and so keeps its edges
    let mut moved = Vec::new();
    for &face in &all_faces {
        let index = faces.iter().position(|&f| f == face);
        let own = face_ops::surface(b, face)?.ok_or_else(null)?;
        let stays = match index.and_then(|i| given.get(i)).and_then(Option::as_ref) {
            Some(new) => is_own_plane(own, new),
            None => true,
        };
        let surface = if let (Some(index), false) = (index, stays) {
            moved.push(face);
            if let Some(slot) = placed.get_mut(index) {
                *slot = Some(new_list.len());
            }
            let mut new = given
                .get_mut(index)
                .and_then(Option::take)
                .ok_or_else(null)?;
            // (not upstream) a face that meets itself across an edge cannot take an open
            // surface: refused before anything is changed
            if meets_itself(b, face)? && !is_closed_surface(k, &new, b.tolerance())? {
                return Err(Failure::new(OP, OffsetFailure::SeamOpened { face }).into_error());
            }
            let domain = new.natural_uv_domain_if_built();
            // (`EvaluateNormal(uv, TRUE, TRUE, ..)`: from the left)
            let new_normal =
                new.evaluate_normal(k, domain.evaluate(0.5, 0.5), Side::Left, Side::Left)?;
            let old = own;
            let old_normal = old.evaluate_normal(
                k,
                b.face(face)?.uv_domain.evaluate(0.5, 0.5),
                Side::Left,
                Side::Left,
            )?;
            // when surface normals oppose - reverse the input surface's orientation
            if old_normal.dot(new_normal) < 0.0 {
                new.reverse(k, SurfParamType::U)?;
            }
            new
        } else {
            // copy the face's current surface, when possible as an analytic one
            own.copy_and_add_analytics(k)?
        };
        new_list.push(surface);
    }

    // the Brep's bounding sphere gives the executive an offset distance it does not use
    let mut bbox = calculate_bounding_box(k, &mut *b, BoundingBoxOptions::default())?;
    // Not upstream (WR-534), opt-in `local_operation_box_holds_new_surfaces` (the fork's
    // tangency-labs/smlib@a1ee273): the box also holds the replacement surfaces. Their
    // faces' edges are found by extending the neighbouring surfaces to meet them, by a
    // distance sized from this box: sized from the original body alone, a face moved
    // further than about its size away left the extended neighbours short of it.
    if k.fixes.local_operation_box_holds_new_surfaces {
        for (j, slot) in placed.iter().enumerate() {
            let surface = match (given.get(j).and_then(Option::as_ref), slot) {
                (Some(surface), _) => surface,
                (None, Some(at)) => new_list.get(*at).ok_or_else(null)?,
                (None, None) => continue,
            };
            let domain = surface.natural_uv_domain_with_nurb(k)?;
            let surf_box = surface
                .calculate_bounding_box(k, &domain, [true, false, false], true)?
                .normal_box
                .ok_or_else(null)?;
            bbox = bbox.union(&surf_box);
        }
    }
    let radius = bbox.sphere_bound().radius;
    let tolerance = b.tolerance();
    let mut executive = OffsetExecutive::new(
        OffsetOperation::LocalOperation,
        radius / 2.0,
        tolerance / 2.0,
    );
    // (the geometry creator's constructor tolerance, half a degree under WR-666's opt-in)
    executive.use_creator_tangency_tol(k);
    // A face's surface is extended as far as its neighbours' new surfaces meet it, and as
    // far again as they moved beside it (not upstream, whose margin is five half-radii of
    // the sphere bounding the box: `local_extension_reach`). The edges whose faces and
    // vertices' faces all keep their surfaces are kept as they are, not intersected again
    // (`OffsetExecutive::keeps_edge`).
    executive.set_local_surfaces(&all_faces, new_list, &moved);
    // the shell faces get copied without being offset: none
    executive.set_shell_faces(&[]);
    // extend/intersect to fill new convex edge gaps rather than adding a fillet
    executive.set_extend_convex_edges(true);

    // (not upstream) `b`'s own attributes go to the offset under the original
    // identities, as upstream's live on the faces it copies (`SM_NO_COPY` of each
    // source face); the caller's store does not keep them afterwards
    let seeded: BTreeSet<ObjectKey> = b
        .attributed_entities()
        .into_iter()
        .map(|e| (notifications.original_key_of)(e).key)
        .collect();
    p2::seed_brep_attributes(
        notifications.attributes,
        b.attribute_store().clone(),
        notifications.original_key_of,
    )?;

    // Build a new Brep by merging in all the faces with their new surfaces
    let offset = executive
        .do_extended_offset(k, Some(&mut *b), notifications)
        .map_err(|e| {
            if e.is_unported_or_cancelled() || k.fixes.local_operation_offset_status {
                e
            } else {
                // status `SM_ERR`, as upstream; the offset's own error is kept (#2615)
                SolidError::failed_after(OP, "the offset of the moved faces could not be built", e)
            }
        });
    drop(executive);
    let new_brep = offset.and_then(|mut new_brep| {
        // the offset B-rep's objects carry what the offset gave them
        p2::restore_keyed_attributes(
            &mut new_brep,
            notifications.attributes,
            notifications.result_key_of,
        )?;
        Ok(new_brep)
    });
    if !seeded.is_empty() {
        notifications
            .attributes
            .retain_rekeyed(&|key| (!seeded.contains(&key)).then_some(key))?;
    }
    let mut new_brep = new_brep?;

    // remove all faces from this Brep, then move all faces of the new Brep into it
    // (`CopyFaces` with a merge, which also grows this Brep's tolerance to the new one's)
    let editing = b.editing_enabled();
    b.set_editing_enabled(true);
    let moved = move_faces(k, b, &mut new_brep);
    b.set_editing_enabled(editing);
    // (not upstream: the port's cached topology fields, true again for the caller)
    b.refresh_cached_topology();
    moved
}

/// Whether `face` meets itself across one of its edges: more edge uses of the edge on
/// the face than its two sides' (a seam of a closed surface, where the face continues
/// round onto itself). A moved face's edges are found again by intersecting the new
/// surfaces of the faces on their two sides; across such an edge both sides are the
/// face itself, so an open surface (a plane) has no curve there: before this was refused,
/// the edge closed onto one vertex grown to the gap (13 units on fc0016's cylinder face
/// given its tangent plane), the shell stayed open, and the body kept was the wrong side
/// (manifold, volume -877.6).
fn meets_itself(b: &Brep, face: FaceId) -> Result<bool> {
    use crate::topology::edgeuse::p1 as eu1;
    for e in crate::topology::face::p4::edges(b, face, None, false, None)? {
        let mut on_face = 0_usize;
        for u in crate::topology::edge::edgeuses(b, e, None, None)? {
            if eu1::face(b, u)? == Some(face) {
                on_face = on_face.saturating_add(1);
            }
        }
        // (each side of the face has its own edge use: two for an ordinary edge)
        if on_face > 2 {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Whether `surface` is closed in either direction over its natural domain (a plane
/// never is; `IsClosed` at `tolerance` for the others, as `SmEdge::IsSeam` asks it).
fn is_closed_surface(k: &Kernel, surface: &Surface, tolerance: f64) -> Result<bool> {
    let domain = surface.natural_uv_domain_with_nurb(k)?;
    for direction in [SurfParamType::U, SurfParamType::V] {
        if surface
            .is_closed(k, &domain, direction, Some(tolerance), false)?
            .0
        {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Whether `given` is the plane `own` lies on, exactly: both are planes, and the plane
/// `own` spans (its origin and its two axes) holds the given plane's origin, and its two
/// axes span the same plane: each determinant is computed without rounding (error-free
/// products and sums, [`exact`]). Not a tolerance: a plane off by a rounding is another
/// plane, and its face moves. (The given plane's extent is not compared: a face's surface
/// is extended as its neighbours need, whatever its extent.)
fn is_own_plane(own: &Surface, given: &Surface) -> bool {
    let (SurfaceKind::Plane(own), SurfaceKind::Plane(given)) = (own.kind(), given.kind()) else {
        return false;
    };
    let (a, b, o) = (
        own.position().x_axis(),
        own.position().y_axis(),
        own.position().origin(),
    );
    let (x, y, g) = (
        given.position().x_axis(),
        given.position().y_axis(),
        given.position().origin(),
    );
    let spans = |u: Vec3, v: Vec3| exact::cross(u, v).iter().any(|c| !exact::is_zero(c));
    // own's axes and the given ones each span a plane, the given axes lie in own's, and
    // the given origin lies on own's plane
    spans(a, b)
        && spans(x, y)
        && exact::is_zero(&exact::det3(a, b, x))
        && exact::is_zero(&exact::det3(a, b, y))
        && exact::is_zero(&exact::sub(&exact::det3(a, b, g), &exact::det3(a, b, o)))
}

/// Error-free arithmetic on doubles: a value is an expansion, a sum of doubles whose
/// components do not overlap (Shewchuk's grow-expansion), so sums and products of doubles
/// are held without rounding and a value's sign is read off exactly (an expansion is zero
/// exactly when every component is). Overflow and underflow are out of range (the inputs
/// are coordinates and unit vectors).
mod exact {
    use crate::math::Vec3;

    /// `a + b` as `(s, e)`: `s = fl(a + b)` and its rounding error (Knuth's TwoSum).
    fn two_sum(a: f64, b: f64) -> (f64, f64) {
        let s = a + b;
        let bb = s - a;
        (s, (a - (s - bb)) + (b - bb))
    }

    /// `e + b` as an expansion (Grow-Expansion).
    fn grow(e: &[f64], b: f64) -> Vec<f64> {
        let mut q = b;
        let mut out = Vec::new();
        for &c in e {
            let (s, err) = two_sum(q, c);
            out.push(err);
            q = s;
        }
        out.push(q);
        out
    }

    /// `e + f` as an expansion.
    fn add(e: &[f64], f: &[f64]) -> Vec<f64> {
        f.iter().fold(e.to_vec(), |acc, &c| grow(&acc, c))
    }

    /// `e - f` as an expansion.
    pub(super) fn sub(e: &[f64], f: &[f64]) -> Vec<f64> {
        f.iter().fold(e.to_vec(), |acc, &c| grow(&acc, -c))
    }

    /// `a b` as an expansion: the product and its rounding error, `fma(a, b, -p)`.
    fn product(a: f64, b: f64) -> Vec<f64> {
        let p = a * b;
        vec![a.mul_add(b, -p), p]
    }

    /// `e b` as an expansion: each component's product, summed.
    fn scale(e: &[f64], b: f64) -> Vec<f64> {
        e.iter()
            .fold(Vec::new(), |acc, &c| add(&acc, &product(c, b)))
    }

    /// `u × v`, each component an expansion.
    pub(super) fn cross(u: Vec3, v: Vec3) -> [Vec<f64>; 3] {
        [
            sub(&product(u.y, v.z), &product(u.z, v.y)),
            sub(&product(u.z, v.x), &product(u.x, v.z)),
            sub(&product(u.x, v.y), &product(u.y, v.x)),
        ]
    }

    /// `det(a, b, c) = (a × b) · c` as an expansion.
    pub(super) fn det3(a: Vec3, b: Vec3, c: Vec3) -> Vec<f64> {
        let [x, y, z] = cross(a, b);
        add(&add(&scale(&x, c.x), &scale(&y, c.y)), &scale(&z, c.z))
    }

    /// Whether the expansion is zero exactly.
    pub(super) fn is_zero(e: &[f64]) -> bool {
        e.iter().all(|&c| c == 0.0)
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn a_rounded_sum_is_not_zero_but_its_exact_one_is() {
            // 0.1 + 0.2 - 0.3 rounds to 5.6e-17; exactly it is that too (the doubles are
            // not the decimals), while 1e16 + 1 - 1e16 is 1, not the 0 or 2 of rounding
            assert!(!is_zero(&sub(&add(&[0.1], &[0.2]), &[0.3])));
            assert_eq!(sub(&add(&[1e16], &[1.0]), &[1e16]).iter().sum::<f64>(), 1.0);
            assert!(is_zero(&sub(&product(3.0, 0.1), &product(0.1, 3.0))));
        }

        #[test]
        fn a_determinant_is_zero_exactly_only_for_dependent_vectors() {
            let a = Vec3::new(0.1, 0.7, 0.3);
            let b = Vec3::new(0.9, 0.2, 0.4);
            assert!(is_zero(&det3(a, b, a)));
            assert!(!is_zero(&det3(
                a,
                b,
                Vec3::new(0.1, 0.7, 0.3 + 1e-16 * 0.3)
            )));
        }
    }
}

/// `RemoveFaces(GetFaces)` of `b`, then `pNewBrep->CopyFaces(GetFaces, b)`.
fn move_faces(k: &mut Kernel, b: &mut Brep, new_brep: &mut Brep) -> Result<()> {
    let old_faces = bp1::faces(b, None)?;
    p5::remove_faces(k, b, true, &old_faces)?;
    let new_faces = bp1::faces(new_brep, None)?;
    // (not upstream: upstream's surface owner is the face itself) the copied faces'
    // surfaces are owned under the identity `b`'s naming scheme names a face by when it
    // records none: the notification identity of a `Notify` B-rep, else the arena
    // identity (tangency-labs/solidcore#2047)
    let naming = b.naming();
    let face_ref = |face| crate::topology::scheme_owner_identity(naming, EntityRef::Face(face));
    // the arena entry carries no attribute records: the copy notifies through an
    // arena-only store (a default record per object of both B-reps) whose source
    // identities are distinct from `b`'s, holding both B-reps' attributes; `b` gets the
    // attributes its objects end with (the merge branch of `p7::copy_faces`)
    crate::boolean::merge::p3::copy_component(k, new_brep, &new_faces, b, &face_ref)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::{Point3, Vec3};
    use crate::primitive::api_primitives::create_box;
    use crate::surface::plane::Plane;

    use crate::geom::attribute::{AttributeStore, ObjectKey};
    use crate::geom::object::ObjectRef;
    use crate::topology::topology::Topology;

    /// Runs `f` with throwaway notification stores keyed by entity id, as callers with no
    /// identities of their own do.
    fn with_notifications<R>(
        b: &mut Brep,
        f: impl FnOnce(&mut Brep, &mut OffsetNotifications<'_>) -> R,
    ) -> R {
        fn key(entity: EntityRef, namespace: u64) -> ObjectRef {
            let (index, type_code) = match entity {
                EntityRef::Brep(id) => (id.index(), crate::geom::a_object::BREP_TYPE),
                EntityRef::Region(id) => (id.index(), crate::topology::Region::type_code()),
                EntityRef::Shell(id) => (id.index(), crate::topology::Shell::type_code()),
                EntityRef::Face(id) => (id.index(), crate::topology::face::p2::TYPE),
                EntityRef::Faceuse(id) => (id.index(), 16022),
                EntityRef::Loopuse(id) => (id.index(), crate::topology::loopuse::TYPE),
                EntityRef::Edge(id) => (id.index(), crate::topology::edge::TYPE),
                EntityRef::Edgeuse(id) => (id.index(), crate::topology::edgeuse::p1::TYPE),
                EntityRef::Vertex(id) => (id.index(), crate::topology::vertex::TYPE),
                EntityRef::Vertexuse(id) => (id.index(), crate::topology::vertexuse::TYPE),
                EntityRef::Loop(id) => (id.index(), crate::topology::topology::TYPE),
            };
            ObjectRef {
                key: ObjectKey((namespace << 60) | ((type_code as u64) << 32) | u64::from(index)),
                type_code,
            }
        }
        let mut records = crate::topology::topology::NotifyRecords::new();
        let mut attributes = AttributeStore::default();
        let mut add = |entity| {
            records.insert(key(entity, 1).key, Topology::new());
        };
        add(EntityRef::Brep(crate::topology::BrepId::from_index(0)));
        b.regions()
            .iter()
            .for_each(|(id, _)| add(EntityRef::Region(id)));
        b.shells()
            .iter()
            .for_each(|(id, _)| add(EntityRef::Shell(id)));
        b.faces()
            .iter()
            .for_each(|(id, _)| add(EntityRef::Face(id)));
        b.faceuses()
            .iter()
            .for_each(|(id, _)| add(EntityRef::Faceuse(id)));
        b.loopuses()
            .iter()
            .for_each(|(id, _)| add(EntityRef::Loopuse(id)));
        b.edges()
            .iter()
            .for_each(|(id, _)| add(EntityRef::Edge(id)));
        b.edgeuses()
            .iter()
            .for_each(|(id, _)| add(EntityRef::Edgeuse(id)));
        b.vertices()
            .iter()
            .for_each(|(id, _)| add(EntityRef::Vertex(id)));
        b.vertexuses()
            .iter()
            .for_each(|(id, _)| add(EntityRef::Vertexuse(id)));
        let result_brep = key(EntityRef::Brep(crate::topology::BrepId::from_index(0)), 3);
        records.insert(result_brep.key, Topology::new());
        let original = |e| key(e, 1);
        let result = |e| key(e, 3);
        let original_curve = |e: crate::topology::EdgeId, t: i64| ObjectRef {
            key: ObjectKey((5_u64 << 60) | ((t as u64) << 32) | u64::from(e.index())),
            type_code: t,
        };
        let result_curve = |e: crate::topology::EdgeId, t: i64| ObjectRef {
            key: ObjectKey((4_u64 << 60) | ((t as u64) << 32) | u64::from(e.index())),
            type_code: t,
        };
        let lists = |_: bool, region: crate::topology::RegionId| {
            [0_u64, 1].map(|list| ObjectRef {
                key: ObjectKey((6_u64 << 60) | (list << 32) | u64::from(region.index())),
                type_code: 0,
            })
        };
        let mut notifications = OffsetNotifications {
            records: &mut records,
            attributes: &mut attributes,
            original_key_of: &original,
            result_key_of: &result,
            original_curve_ref: &original_curve,
            result_curve_ref: &result_curve,
            boolean_source_lists: &lists,
            result_brep,
            merge_identities: None,
        };
        f(b, &mut notifications)
    }

    #[test]
    fn moves_the_top_of_a_box_and_keeps_it_a_solid() {
        let mut k = Kernel::default();
        let mut b = create_box(&mut k, Point3::new(0.0, 0.0, 0.0), 2.0, 2.0, 2.0).unwrap();
        let top = bp1::faces(&b, None)
            .unwrap()
            .into_iter()
            .find(|&f| {
                let s = face_ops::surface(&b, f).unwrap().unwrap();
                let p = s
                    .evaluate_point(
                        &k,
                        s.natural_uv_domain_with_nurb(&k)
                            .unwrap()
                            .evaluate(0.5, 0.5),
                    )
                    .unwrap();
                (p.z - 2.0).abs() < 1e-9
            })
            .unwrap();
        let raised = Surface::from(Plane::from_normal(Vec3::new(0.0, 0.0, 3.0), Vec3::Z));
        with_notifications(&mut b, |b, notes| {
            local_operation(&mut k, b, &[top], vec![raised], notes)
        })
        .unwrap();
        assert_eq!(bp1::faces(&b, None).unwrap().len(), 6);
        // the cached topology fields agree with the links (the public integrity check)
        crate::topology::integrity::check(&b).unwrap();
        let e = calculate_bounding_box(&k, &mut b, BoundingBoxOptions::default()).unwrap();
        // the box carries the body's tolerance
        let near = |a: f64, b: f64| (a - b).abs() < 1e-4;
        assert!(near(e.max_at(2), 3.0) && near(e.min_at(2), 0.0) && near(e.max_at(0), 2.0));
    }

    /// A cylinder's side meets itself across its seam: given its tangent plane (open), it
    /// is refused before anything changes (fc0016 of the difftest unit kept the wrong side
    /// of an open shell, volume -877.6).
    #[test]
    fn a_seamed_face_given_a_plane_is_refused_unchanged() {
        let mut k = Kernel::default();
        let mut b =
            crate::primitive::api_primitives::create_cylinder(&mut k, Point3::ZERO, 1.0, 2.0)
                .unwrap();
        let side = bp1::faces(&b, None)
            .unwrap()
            .into_iter()
            .find(|&f| {
                !matches!(
                    face_ops::surface(&b, f).unwrap().unwrap().kind(),
                    SurfaceKind::Plane(_)
                )
            })
            .unwrap();
        assert!(meets_itself(&b, side).unwrap());
        let before = (b.faces().len(), b.edges().len(), b.vertices().len());
        let tangent = Surface::from(Plane::from_normal(Vec3::new(1.0, 0.0, 1.0), Vec3::X));
        let e = with_notifications(&mut b, |b, notes| {
            local_operation(&mut k, b, &[side], vec![tangent], notes)
        })
        .unwrap_err();
        assert_eq!(e.kind(), crate::error::ErrorKind::Unsupported, "{e}");
        assert!(matches!(
            e.failure().map(|f| &f.reason),
            Some(crate::error::Reason::Offset(OffsetFailure::SeamOpened { face })) if *face == side
        ));
        assert_eq!(
            (b.faces().len(), b.edges().len(), b.vertices().len()),
            before
        );
    }

    #[test]
    fn a_short_surface_list_is_a_null_pointer() {
        let mut k = Kernel::default();
        let mut b = create_box(&mut k, Point3::new(0.0, 0.0, 0.0), 1.0, 1.0, 1.0).unwrap();
        let face = bp1::faces(&b, None).unwrap()[0];
        let r = with_notifications(&mut b, |b, notes| {
            local_operation(&mut k, b, &[face], Vec::new(), notes)
        });
        assert!(r.is_err());
    }

    /// WR-534 (the fork's tangency-labs/smlib@a1ee273): a 20 mm cube's top moved up by
    /// 60, further than the cube's size, is a 20 x 20 x 80 box with the opt-in.
    #[test]
    fn moves_the_top_of_a_box_further_than_its_size() {
        let mut k = Kernel::default();
        k.fixes.local_operation_box_holds_new_surfaces = true;
        let mut b = create_box(&mut k, Point3::new(0.0, 0.0, 0.0), 20.0, 20.0, 20.0).unwrap();
        let top = bp1::faces(&b, None)
            .unwrap()
            .into_iter()
            .find(|&f| {
                let s = face_ops::surface(&b, f).unwrap().unwrap();
                let p = s
                    .evaluate_point(
                        &k,
                        s.natural_uv_domain_with_nurb(&k)
                            .unwrap()
                            .evaluate(0.5, 0.5),
                    )
                    .unwrap();
                (p.z - 20.0).abs() < 1e-9
            })
            .unwrap();
        let raised = Surface::from(Plane::new(
            Vec3::new(0.0, 0.0, 80.0),
            Vec3::X,
            Vec3::Y,
            crate::math::Vec2::new(20.0, 20.0),
            &crate::geom::extent2d::Extent2d::from_bounds(0.0, 0.0, 1.0, 1.0),
        ));
        with_notifications(&mut b, |b, notes| {
            local_operation(&mut k, b, &[top], vec![raised], notes)
        })
        .unwrap();
        assert_eq!(bp1::faces(&b, None).unwrap().len(), 6);
        let e = calculate_bounding_box(&k, &mut b, BoundingBoxOptions::default()).unwrap();
        let near = |a: f64, b: f64| (a - b).abs() < 1e-3;
        assert!(
            near(e.max_at(2), 80.0) && near(e.min_at(2), 0.0) && near(e.max_at(0), 20.0),
            "{e:?}"
        );
    }
}
