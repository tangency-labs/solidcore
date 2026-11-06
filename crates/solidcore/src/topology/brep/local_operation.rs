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

use std::collections::{BTreeMap, BTreeSet};

use crate::error::offset::OffsetFailure;
use crate::error::{EntityRef, Failure, Result, SolidError};
use crate::geom::attribute::ObjectKey;
use crate::geom::axis2_placement::Axis2Placement;
use crate::geom::extent2d::Extent2d;
use crate::kernel::Kernel;
use crate::math::Vec3;
use crate::nurbs::b_spline_surface::p1::SurfParamType;
use crate::nurbs::knots::Side;
use crate::offset::offset_executive::p1::{OffsetExecutive, OffsetOperation};
use crate::offset::offset_executive::p2::OffsetNotifications;
use crate::query::brep::{calculate_bounding_box, BoundingBoxOptions};
use crate::surface::implicit::{self, AnalyticSupport};
use crate::surface::plane::Plane;
use crate::surface::surface::Surface;
use crate::surface::surface::SurfaceKind;
use crate::topology::brep::p1 as bp1;
use crate::topology::brep::p2;
use crate::topology::brep::p5;
use crate::topology::brep::validation::{self, FailureClass};
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
/// rebuilt by an extended offset that extends the surfaces to intersect. (Not upstream) a
/// face whose B-spline surface has a pole at a vertex of a moved face, and lies on a plane,
/// continues as that plane (`continued_planes`).
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
/// [`OffsetFailure::SeamOpened`] (not upstream) when a face that meets itself across
/// edges (seams) is given a surface closed in fewer directions than it meets itself in,
/// before anything changes; [`OffsetFailure::Unsound`] (not upstream) when the body built
/// fails what `b` did not (`refuse_unsound`: a check failure class `b`'s check has fewer
/// lines of, a geometry the check could not reach, a manifold solid lost, a volume not
/// positive); those of the callees ([`crate::offset::offset_executive::p2`], the surfaces'
/// normals, reversal and analytic copy, the bounding box). On every error `b` is as it
/// came in, the caches its surfaces had included (the operation edits a copy of it), and
/// so are `notifications.records` and `notifications.attributes`.
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
    // (not upstream) a refusal or a failure leaves `b` as it came in, the caches its
    // surfaces had included (the operation runs on its copy), and the caller's
    // notification stores as they were
    let records = notifications.records.clone();
    let attributes = notifications.attributes.clone();
    let done = b.with_rollback_untouched(
        |b| local_operation_on(k, b, faces, surfaces, notifications),
        Result::is_ok,
    );
    if done.is_err() {
        *notifications.records = records;
        *notifications.attributes = attributes;
    }
    done
}

/// [`local_operation`] on the copy of `b` it edits.
fn local_operation_on(
    k: &mut Kernel,
    b: &mut Brep,
    faces: &[FaceId],
    surfaces: Vec<Surface>,
    notifications: &mut OffsetNotifications<'_>,
) -> Result<()> {
    // (not upstream) the input's standing, measured once before anything changes: the
    // final gate (`refuse_unsound`) holds the operation to what it adds
    let input = validation::standing(k, &mut b.clone_standing_for_self()?)?;
    // Make a new list of all surfaces in this Brep: for every face, the new surface if it
    // is in the list to be replaced (turned to the original's side), else a copy of its own.
    let all_faces = bp1::faces(b, None)?;
    let mut given: Vec<Option<Surface>> = surfaces.into_iter().map(Some).collect();
    let mut new_list = Vec::with_capacity(all_faces.len());
    // (WR-534) where each given surface went in the new list
    let mut placed: Vec<Option<usize>> = vec![None; given.len()];
    // (not upstream) the faces whose surfaces move: a face given its own surface again keeps
    // a copy of its own surface, as the faces not given one do, and so keeps its edges
    let mut moved = Vec::new();
    // (not upstream) the faces given another surface, and the planes faces that keep
    // theirs continue as beside them (`continued_planes`)
    let mut moving = Vec::new();
    for &face in &all_faces {
        let index = faces.iter().position(|&f| f == face);
        let own = face_ops::surface(b, face)?.ok_or_else(null)?;
        if let Some(new) = index.and_then(|i| given.get(i)).and_then(Option::as_ref) {
            if !is_own_surface(own, new) {
                moving.push(face);
            }
        }
    }
    let mut continued = continued_planes(k, b, &all_faces, &moving)?;
    for &face in &all_faces {
        let index = faces.iter().position(|&f| f == face);
        let own = face_ops::surface(b, face)?.ok_or_else(null)?;
        let replacement = match index {
            Some(index) if moving.contains(&face) => {
                if let Some(slot) = placed.get_mut(index) {
                    *slot = Some(new_list.len());
                }
                Some(
                    given
                        .get_mut(index)
                        .and_then(Option::take)
                        .ok_or_else(null)?,
                )
            }
            _ => continued.remove(&face),
        };
        let surface = if let Some(mut new) = replacement {
            moved.push(face);
            // (not upstream) a face that meets itself across seams cannot take a surface
            // closed in fewer directions: refused before anything is changed
            if closed_directions(k, &new, b.tolerance())? < seam_directions(k, b, face)? {
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
    // (`CopyFaces` with a merge, which also grows this Brep's tolerance to the new one's);
    // (not upstream) a body that is not sound, or a failed move, is refused
    // (`local_operation` puts `b` back)
    let editing = b.editing_enabled();
    b.set_editing_enabled(true);
    let moved = move_faces(k, b, &mut new_brep);
    b.set_editing_enabled(editing);
    // (not upstream: the port's cached topology fields, true again for the caller)
    b.refresh_cached_topology();
    moved?;
    refuse_unsound(k, b, &input)
}

/// (not upstream) The local operation's final gate: the body it built is refused when it
/// fails what its input did not ([`validation::soundness`] of both, `input` measured
/// before the edit, the body on a private copy):
///
/// - the check at `SM_LEVEL_2` has more failure lines of some class (the failed test and
///   the type of the object it is about, [`FailureClass`]) than the input's check had
///   ([`BrepCheck::failures_beyond`](validation::BrepCheck::failures_beyond));
/// - the body's pointer check failed, so its geometry was not checked (`AssertValid` is
///   skipped then), even where the input's pointer check failed alike;
/// - the input was a manifold solid and the body is not;
/// - the body is a manifold solid whose volume is not positive (the input's volume is
///   not measured).
///
/// The owners the check tests the geometry against are the identities the body's naming
/// scheme names its entities by ([`validation::standing`]), as `move_faces` sets them.
///
/// Failures are compared by class, not by entity: the operation rebuilds every face (the
/// extended offset builds them all, unmoved ones on copies of their surfaces, and they
/// are copied into `b`), so no identity persists to compare by. A failure the input had
/// and the body still has is the input's; one the input had on a face the operation
/// rebuilt and the body no longer has is gone; a class the body has more of is the
/// operation's. (A failure gone from one entity while one of its class appears on
/// another is counted as the input's.)
///
/// Each moved edge is found again where extended surfaces meet; where that cannot close
/// the body (a surface that meets its neighbours nowhere near the old edges, a vertex
/// whose faces no longer meet in one point), the merge closes what is left by growing
/// tolerances, and the body it leaves is no solid: refused, not returned.
///
/// # Errors
///
/// [`OffsetFailure::Unsound`] naming what the body fails; those of the copy and the
/// measurement.
fn refuse_unsound(k: &mut Kernel, b: &Brep, input: &validation::Standing) -> Result<()> {
    use crate::error::offset::Unsoundness;
    let mut copy = b.clone_standing_for_self()?;
    let found = validation::standing(k, &mut copy)?;
    let added = found.check.failures_beyond(&input.check);
    let what = if let Some(&(class, _)) = added.first() {
        Unsoundness::Check {
            failures: added.iter().fold(0, |n, &(_, more)| n.saturating_add(more)),
            first: match class {
                FailureClass::Report {
                    reporting_type,
                    test_index,
                    ..
                } => Some((reporting_type, test_index)),
                _ => None,
            },
        }
    } else if found.check.assert_valid.is_none() {
        // the pointer check failed, so `AssertValid` did not run: the body's geometry is
        // not checked at all, whatever the input's pointer failures were
        Unsoundness::GeometryUnchecked
    } else if input.manifold && !found.manifold {
        Unsoundness::NotManifoldSolid
    } else if found.manifold
        && validation::volume(k, &mut copy)?.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater)
    {
        Unsoundness::VolumeNotPositive
    } else {
        return Ok(());
    };
    Err(Failure::new(OP, OffsetFailure::Unsound { what }).into_error())
}

/// (not upstream) The planes faces that keep their surfaces continue as: a face whose
/// B-spline surface has a pole (a side collapsed to a point) at a vertex of a moved face
/// (one of `moving`), and whose net lies on a plane ([`implicit::plane_of_net`],
/// certified to the rounding of its stored points), takes that plane, and moves with
/// the moved faces (its edges are intersected again, as those of a face given its own
/// plane are). The patch cannot be extended across its pole: a side's extension is a
/// strip of parameters whose width in space goes to nothing at the pole, so where the
/// moved faces split the vertex there (four faces at a pyramid's apex no longer meet in
/// one point, and two of them meet in a new edge), the split's vertices next to the pole
/// lie on no extension of it, and its edges to them stopped where they left the strip:
/// a vertex too many, and edges on one face alone. The plane holds them. A face that
/// meets itself across an edge (a seam), whose surface is closed, is left as it is.
///
/// # Errors
///
/// Those of the topology queries, [`crate::topology::vertexuse::is_on_pole`] and
/// [`implicit::plane_of_net`].
fn continued_planes(
    k: &Kernel,
    b: &mut Brep,
    all_faces: &[FaceId],
    moving: &[FaceId],
) -> Result<BTreeMap<FaceId, Surface>> {
    use crate::topology::{face::p4 as face4, vertex, vertexuse};
    let mut out = BTreeMap::new();
    for &face in all_faces {
        if moving.contains(&face)
            || !matches!(
                face_ops::surface(b, face)?.map(Surface::kind),
                Some(SurfaceKind::BSpline(_))
            )
            || meets_itself(b, face)?
        {
            continue;
        }
        let mut pole = false;
        for v in face4::vertices(b, face, None, false, None)? {
            if pole || !vertex::faces(b, v)?.iter().any(|f| moving.contains(f)) {
                continue;
            }
            for vu in vertex::vertexuses_of_face(b, v, Some(face), None)? {
                pole |= vertexuse::is_on_pole(k, b, vu)?;
            }
        }
        if !pole {
            continue;
        }
        let Some(SurfaceKind::BSpline(form)) = face_ops::surface(b, face)?.map(Surface::kind)
        else {
            continue;
        };
        if let Some(implicit::Implicit::Plane { origin, normal }) = implicit::plane_of_net(form)? {
            out.insert(face, Surface::from(Plane::from_normal(origin, normal)));
        }
    }
    Ok(out)
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

/// The number of directions `surface` is closed in over `domain` (a plane none; `IsClosed`
/// at `tolerance` for the others, as `SmEdge::IsSeam` asks it).
fn closed_in(k: &Kernel, surface: &Surface, domain: &Extent2d, tolerance: f64) -> Result<usize> {
    let mut n = 0_usize;
    for direction in [SurfParamType::U, SurfParamType::V] {
        if surface
            .is_closed(k, domain, direction, Some(tolerance), false)?
            .0
        {
            n = n.saturating_add(1);
        }
    }
    Ok(n)
}

/// The number of directions `surface` is closed in over its natural domain.
fn closed_directions(k: &Kernel, surface: &Surface, tolerance: f64) -> Result<usize> {
    closed_in(
        k,
        surface,
        &surface.natural_uv_domain_with_nurb(k)?,
        tolerance,
    )
}

/// The number of directions `face` meets itself across seams in: none for a face that
/// does not ([`meets_itself`]), else the directions its own surface is closed in over the
/// face's UV domain (a torus face round both ways two, a cylinder's side, or a torus face
/// cut across its tube, one), and at least the one its seam is in. The new surface must be
/// closed in as many: a seam has no edge on a surface open across it.
fn seam_directions(k: &Kernel, b: &Brep, face: FaceId) -> Result<usize> {
    if !meets_itself(b, face)? {
        return Ok(0);
    }
    let own = face_ops::surface(b, face)?.ok_or_else(null)?;
    Ok(closed_in(k, own, &b.face(face)?.uv_domain, b.tolerance())?.max(1))
}

/// Whether `given` is the surface `own` lies on, exactly: both of one analytic kind
/// ([`AnalyticSupport`]) with equal radii and semi-angle, and the same placement where it
/// matters, each comparison made without rounding (error-free products and sums,
/// [`exact`]):
/// - a plane: the plane `own` spans (its origin and its two axes) holds the given plane's
///   origin, and its two axes span the same plane;
/// - a cylinder: the same axis line (either sense);
/// - a cone: the same axis line and sense and the same origin (where its radius is the
///   base radius), or the opposite sense with the semi-angle negated;
/// - a sphere: the same centre;
/// - a torus: the same centre and axis line (either sense).
///
/// Not a tolerance: a surface off by a rounding is another surface, and its face moves.
/// The extent and parametrisation are not compared (a face's surface is extended as its
/// neighbours need, whatever its extent, and its own copy keeps its own parameters).
fn is_own_surface(own: &Surface, given: &Surface) -> bool {
    use AnalyticSupport as A;
    let (Some(own), Some(given)) = (AnalyticSupport::of(own), AnalyticSupport::of(given)) else {
        return false;
    };
    match (own, given) {
        (A::Plane { placement: p, .. }, A::Plane { placement: q, .. }) => same_plane(&p, &q),
        (
            A::Cylinder {
                placement: p,
                radius: r,
            },
            A::Cylinder {
                placement: q,
                radius: s,
            },
        ) => r == s && same_axis(&p, &q),
        (
            A::Cone {
                placement: p,
                base_radius: r,
                semi_angle_deg: a,
            },
            A::Cone {
                placement: q,
                base_radius: s,
                semi_angle_deg: b,
            },
        ) => {
            // (the axes are parallel exactly, so their unit z axes' dot is +-1 to
            // rounding: its sign is the sense)
            let same_sense = p.z_axis().dot(q.z_axis()) > 0.0;
            r == s
                && (if same_sense { a == b } else { a == -b })
                && same_point(p.origin(), q.origin())
                && same_axis(&p, &q)
        }
        (
            A::Sphere {
                placement: p,
                radius: r,
            },
            A::Sphere {
                placement: q,
                radius: s,
            },
        ) => r == s && same_point(p.origin(), q.origin()),
        (
            A::Torus {
                placement: p,
                major_radius: r1,
                minor_radius: r2,
            },
            A::Torus {
                placement: q,
                major_radius: s1,
                minor_radius: s2,
            },
        ) => r1 == s1 && r2 == s2 && same_point(p.origin(), q.origin()) && same_axis(&p, &q),
        _ => false,
    }
}

fn same_point(a: Vec3, b: Vec3) -> bool {
    a.x == b.x && a.y == b.y && a.z == b.z
}

/// Whether `p`'s axes span a plane, exactly.
fn spans(p: &Axis2Placement) -> bool {
    exact::cross(p.x_axis(), p.y_axis())
        .iter()
        .any(|c| !exact::is_zero(c))
}

/// Whether the plane `p` spans (its origin and its two axes) is the one `q` spans: `q`'s
/// axes lie in `p`'s and its origin on it, each determinant exact.
fn same_plane(p: &Axis2Placement, q: &Axis2Placement) -> bool {
    let (a, b, o) = (p.x_axis(), p.y_axis(), p.origin());
    let (x, y, g) = (q.x_axis(), q.y_axis(), q.origin());
    spans(p)
        && spans(q)
        && exact::is_zero(&exact::det3(a, b, x))
        && exact::is_zero(&exact::det3(a, b, y))
        && exact::is_zero(&exact::sub(&exact::det3(a, b, g), &exact::det3(a, b, o)))
}

/// Whether `p`'s z axis line (through its origin along `x × y`) is `q`'s, either sense:
/// `q`'s axis `x' × y'` is normal to `p`'s x and y axes (so parallel to `x × y`), and the
/// origins' difference is too, each product exact.
fn same_axis(p: &Axis2Placement, q: &Axis2Placement) -> bool {
    let (a, b, o) = (p.x_axis(), p.y_axis(), p.origin());
    let (x, y, g) = (q.x_axis(), q.y_axis(), q.origin());
    spans(p)
        && spans(q)
        && exact::is_zero(&exact::det3(x, y, a))
        && exact::is_zero(&exact::det3(x, y, b))
        && exact::is_zero(&exact::sub(&exact::dot(a, g), &exact::dot(a, o)))
        && exact::is_zero(&exact::sub(&exact::dot(b, g), &exact::dot(b, o)))
}

/// Error-free arithmetic on doubles: a value is an expansion, a sum of doubles whose
/// components do not overlap (Shewchuk's grow-expansion), so sums and products of doubles
/// are held without rounding and a value's sign is read off exactly (an expansion is zero
/// exactly when every component is). Overflow and underflow are out of range (the inputs
/// are coordinates and unit vectors).
mod exact {
    use crate::math::Vec3;

    /// `a + b` as `(s, e)`: `s = fl(a + b)` and its rounding error (Knuth's `TwoSum`).
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

    /// `u · v` as an expansion.
    pub(super) fn dot(u: Vec3, v: Vec3) -> Vec<f64> {
        add(
            &add(&product(u.x, v.x), &product(u.y, v.y)),
            &product(u.z, v.z),
        )
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
    crate::boolean::merge::p3::copy_component(k, new_brep, &new_faces, b, &face_ref)?;
    // (not upstream, as for the surfaces) the copied edges' curves and edge uses' UV trim
    // curves are owned under the copy's notification identities; a B-rep not named by
    // them names its edges and edge uses by the arena identities: rebound to those, so
    // the body's geometry names the entities it lies on (`SmBrep` test #61, Bad Edge)
    crate::topology::brep::p3::rebind_geometry_owners(
        b,
        crate::topology::loopuse::notify_ref,
        |e| crate::topology::scheme_owner_identity(naming, e),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::{Point3, Vec3};
    use crate::primitive::api_primitives::create_box;
    use crate::surface::plane::Plane;
    use crate::surface::surface::SurfaceKind;

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

    /// A body named by its arena identities (as an imported one is): the edges copied in
    /// from the offset own their curves under those identities, not the copy's
    /// notification ones, so the body passes its check under its own naming, and the
    /// operation returns it (the 20 mm box, its top raised to z = 25: 10000 mm³).
    #[test]
    fn an_arena_named_body_keeps_its_naming_on_the_edges_copied_in() {
        let mut k = Kernel::default();
        let mut b = create_box(&mut k, Point3::new(0.0, 0.0, 0.0), 20.0, 20.0, 20.0).unwrap();
        b.use_arena_naming().unwrap();
        let top = box_top(&k, &b);
        let raised = Surface::from(Plane::from_normal(Vec3::new(0.0, 0.0, 25.0), Vec3::Z));
        let why = refused_unchanged_or_returned(&mut k, &mut b, &[top], vec![raised]);
        assert_eq!(why, None);
        let lines = validation::check_brep(
            &mut k,
            &mut b.clone_standing_for_self().unwrap(),
            validation::CheckLevel::Level2,
        )
        .unwrap()
        .lines();
        assert!(lines.is_empty(), "{lines:?}");
        let (v, manifold) = volume_and_manifold(&k, &mut b);
        assert!(manifold && (v - 10000.0).abs() <= 1e-9 * 10000.0, "{v}");
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

    /// The top face of a 20 mm cube.
    fn box_top(k: &Kernel, b: &Brep) -> FaceId {
        bp1::faces(b, None)
            .unwrap()
            .into_iter()
            .find(|&f| {
                let s = face_ops::surface(b, f).unwrap().unwrap();
                let p = s
                    .evaluate_point(
                        k,
                        s.natural_uv_domain_with_nurb(k).unwrap().evaluate(0.5, 0.5),
                    )
                    .unwrap();
                (p.z - 20.0).abs() < 1e-9
            })
            .unwrap()
    }

    /// WR-534 (the fork's tangency-labs/smlib@a1ee273): a 20 mm cube's top moved up by
    /// 60, further than the cube's size, is a sound 20 x 20 x 80 box with the opt-in.
    #[test]
    fn moves_the_top_of_a_box_further_than_its_size() {
        let mut k = Kernel::default();
        k.fixes.local_operation_box_holds_new_surfaces = true;
        let mut b = create_box(&mut k, Point3::new(0.0, 0.0, 0.0), 20.0, 20.0, 20.0).unwrap();
        let top = box_top(&k, &b);
        let raised = Surface::from(Plane::from_normal(Vec3::new(0.0, 0.0, 80.0), Vec3::Z));
        let why = refused_unchanged_or_returned(&mut k, &mut b, &[top], vec![raised]);
        assert_eq!(why, None);
        assert_eq!(bp1::faces(&b, None).unwrap().len(), 6);
        let (v, _) = volume_and_manifold(&k, &mut b);
        assert!((v - 32000.0).abs() <= 1e-9 * 32000.0, "{v}");
    }

    /// The same top given a plane parametrised over the unit square, 20 mm a unit: the
    /// body built leaves the corner over (20, 20) unbuilt (the vertical edge there ends at
    /// the old top, z = 20, and the moved top's two edges into that corner end 60 from
    /// their vertex), so it is refused as not sound, the cube as it came in (a defect of
    /// its own, not of the gate: the plane through the same points parametrised by length
    /// gives the sound box above).
    #[test]
    fn a_far_top_on_a_unit_square_plane_is_refused_as_not_sound() {
        let mut k = Kernel::default();
        k.fixes.local_operation_box_holds_new_surfaces = true;
        let mut b = create_box(&mut k, Point3::new(0.0, 0.0, 0.0), 20.0, 20.0, 20.0).unwrap();
        let top = box_top(&k, &b);
        let raised = Surface::from(Plane::new(
            Vec3::new(0.0, 0.0, 80.0),
            Vec3::X,
            Vec3::Y,
            crate::math::Vec2::new(20.0, 20.0),
            &Extent2d::from_bounds(0.0, 0.0, 1.0, 1.0),
        ));
        let why = refused_unchanged_or_returned(&mut k, &mut b, &[top], vec![raised]);
        assert_eq!(why, Some(unsound_check(3, 61)));
    }

    /// The volume and whether the body is a manifold solid, with the cached topology
    /// fields agreeing with the links.
    fn volume_and_manifold(k: &Kernel, b: &mut Brep) -> (f64, bool) {
        crate::topology::integrity::check(b).unwrap();
        let volume = crate::mass::integrator::properties::brep_volume_default(k, b, 1e-6).unwrap();
        (
            volume,
            crate::query::api_queries::brep_is_manifold_solid(b).unwrap(),
        )
    }

    /// The face of `b` whose surface is not a plane (a primitive's side).
    fn curved_face(b: &Brep) -> FaceId {
        bp1::faces(b, None)
            .unwrap()
            .into_iter()
            .find(|&f| {
                !matches!(
                    face_ops::surface(b, f).unwrap().unwrap().kind(),
                    SurfaceKind::Plane(_)
                )
            })
            .unwrap()
    }

    /// `b`'s curved face given `surface`: Ok, a manifold solid with the input's volume to
    /// rounding (the closed form `volume` to the integrator's accuracy, 1e-6), faces, edges,
    /// vertices and tolerance.
    fn given_its_own(mut k: Kernel, mut b: Brep, surface: Surface, volume: f64) {
        let face = curved_face(&b);
        let before = (b.faces().len(), b.edges().len(), b.vertices().len());
        let tolerance = b.tolerance();
        let (input, _) = volume_and_manifold(&k, &mut b);
        assert!(
            (input - volume).abs() <= 1e-6 * volume,
            "{input} vs {volume}"
        );
        with_notifications(&mut b, |b, notes| {
            local_operation(&mut k, b, &[face], vec![surface], notes)
        })
        .unwrap();
        let (v, manifold) = volume_and_manifold(&k, &mut b);
        assert!(manifold);
        assert!((v - input).abs() <= 1e-12 * input, "{v} vs {input}");
        assert_eq!(
            (b.faces().len(), b.edges().len(), b.vertices().len()),
            before
        );
        assert_eq!(b.tolerance(), tolerance);
    }

    /// Review finding 1: a cylinder's side given the same cylinder (another surface of
    /// it, over another extent) keeps its surface and edges: the body is unchanged.
    #[test]
    fn a_cylinder_given_its_own_cylinder_is_unchanged() {
        use crate::primitive::api_primitives as ap;
        let mut k = Kernel::default();
        let b = ap::create_cylinder(&mut k, Point3::ZERO, 1.0, 2.0).unwrap();
        let s = Surface::from(
            ap::create_cylinder_surface(&k, Point3::new(0.0, 0.0, -1.0), 1.0, 4.0).unwrap(),
        );
        given_its_own(k, b, s, 2.0 * std::f64::consts::PI);
    }

    #[test]
    fn a_cone_given_its_own_cone_is_unchanged() {
        use crate::primitive::api_primitives as ap;
        let mut k = Kernel::default();
        let b = ap::create_cone(&mut k, Point3::ZERO, 2.0, 1.0, 3.0).unwrap();
        let s = Surface::from(ap::create_cone_surface(&k, Point3::ZERO, 2.0, 1.0, 3.0).unwrap());
        given_its_own(
            k,
            b,
            s,
            std::f64::consts::PI * 3.0 / 3.0 * (4.0 + 2.0 + 1.0),
        );
    }

    #[test]
    fn a_sphere_given_its_own_sphere_is_unchanged() {
        use crate::primitive::api_primitives as ap;
        let mut k = Kernel::default();
        let b = ap::create_sphere(&mut k, Point3::new(1.0, 2.0, 3.0), 1.5).unwrap();
        let other = ap::create_sphere(&mut k, Point3::new(1.0, 2.0, 3.0), 1.5).unwrap();
        let s = face_ops::surface(&other, curved_face(&other))
            .unwrap()
            .unwrap()
            .clone();
        given_its_own(k, b, s, 4.0 / 3.0 * std::f64::consts::PI * 1.5_f64.powi(3));
    }

    #[test]
    fn a_torus_given_its_own_torus_is_unchanged() {
        use crate::primitive::api_primitives as ap;
        let mut k = Kernel::default();
        let b = ap::create_torus(&mut k, Point3::ZERO, 3.0, 1.0).unwrap();
        let s = Surface::from(ap::create_torus_surface(&k, Point3::ZERO, 3.0, 1.0).unwrap());
        given_its_own(k, b, s, 2.0 * std::f64::consts::PI.powi(2) * 3.0);
    }

    /// A cylinder's side given the coaxial cylinder one rounding wider moves, by a reach
    /// of rounding size: the extension of the neighbours by that tiny margin is built, and
    /// the body is the wider cylinder.
    #[test]
    fn a_cylinder_given_one_a_rounding_wider_moves_by_it() {
        use crate::primitive::api_primitives as ap;
        let mut k = Kernel::default();
        let mut b = ap::create_cylinder(&mut k, Point3::ZERO, 1.0, 2.0).unwrap();
        let r = f64::from_bits(1.0_f64.to_bits() + 1);
        let s = Surface::from(
            ap::create_cylinder_surface(&k, Point3::new(0.0, 0.0, -1.0), r, 4.0).unwrap(),
        );
        let face = curved_face(&b);
        let tolerance = b.tolerance();
        let (input, _) = volume_and_manifold(&k, &mut b);
        with_notifications(&mut b, |b, notes| {
            local_operation(&mut k, b, &[face], vec![s], notes)
        })
        .unwrap();
        let (v, manifold) = volume_and_manifold(&k, &mut b);
        assert!(manifold);
        // no edge is kept (each touches the side), and the body keeps the input's tolerance
        // all the same (review finding 7)
        assert_eq!(b.tolerance(), tolerance);
        // the new caps' edges are re-intersected B-splines within the body's tolerance of
        // the circles (as for any moved cylinder side), so the volume is within the
        // tolerance times the area (2 pi r h + 2 pi r^2) of the input's
        let bound = b.tolerance() * 6.0 * std::f64::consts::PI;
        assert!((v - input).abs() <= bound, "{v} vs {input}");
    }

    /// Review findings 7 and 8: a box's top raised, its bottom's four edges kept (their
    /// faces and their vertices' faces all stay). Each is the original edge in the result:
    /// the same curve over the same interval, its tolerance and its vertices' not below
    /// the input's; the body's tolerance is the input's; the result is a manifold solid.
    #[test]
    fn a_raised_top_keeps_the_bottom_edges_as_they_were() {
        let mut k = Kernel::default();
        let mut b = create_box(&mut k, Point3::ZERO, 2.0, 2.0, 2.0).unwrap();
        let input = create_box(&mut k, Point3::ZERO, 2.0, 2.0, 2.0).unwrap();
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
        let (v, manifold) = volume_and_manifold(&k, &mut b);
        assert!(manifold && (v - 12.0).abs() <= 1e-6 * 12.0, "{v}");
        assert_eq!(b.tolerance(), input.tolerance());
        // the curve's points at its interval's ends and quarters, bit for bit
        let points = |br: &Brep, e: crate::topology::EdgeId| {
            let c = crate::topology::edge::curve(br, e).unwrap().unwrap();
            let i = crate::topology::edge::interval(br, e).unwrap();
            [0.0, 0.25, 0.5, 0.75, 1.0].map(|t| {
                c.evaluate_point(&k, i.min() + t * (i.max() - i.min()))
                    .unwrap()
            })
        };
        let bottom: Vec<_> = input
            .edges()
            .iter()
            .map(|(e, _)| e)
            .filter(|&e| points(&input, e).iter().all(|p| p.z == 0.0))
            .collect();
        assert_eq!(bottom.len(), 4);
        for e in bottom {
            let want = points(&input, e);
            let interval = crate::topology::edge::interval(&input, e).unwrap();
            let kept: Vec<_> = b
                .edges()
                .iter()
                .map(|(f, _)| f)
                .filter(|&f| points(&b, f) == want)
                .collect();
            let &[f] = kept.as_slice() else {
                panic!("edge {e:?}: {} copies", kept.len());
            };
            let got = crate::topology::edge::interval(&b, f).unwrap();
            assert_eq!(
                (got.min().to_bits(), got.max().to_bits()),
                (interval.min().to_bits(), interval.max().to_bits())
            );
            assert!(b.edge(f).unwrap().tolerance >= input.edge(e).unwrap().tolerance);
            let (s0, s1) = crate::topology::edge::start_end_vertices(&input, e).unwrap();
            let (t0, t1) = crate::topology::edge::start_end_vertices(&b, f).unwrap();
            for (s, t) in [(s0, t0), (s1, t1)] {
                let (s, t) = (
                    input.vertex(s.unwrap()).unwrap(),
                    b.vertex(t.unwrap()).unwrap(),
                );
                assert_eq!(s.point, t.point);
                assert!(t.tolerance >= s.tolerance);
            }
        }
    }

    /// A pyramid as the difftest unit's `Solids/pyramid.smb` is made: a rectangular base
    /// and four triangular sides, each face a bicubic B-spline patch over `[0, 1]^2` (a
    /// side's far row collapsed onto the apex), at the file's B-rep tolerance.
    fn bicubic_pyramid(k: &mut Kernel) -> Brep {
        use crate::geom::extent2d::Extent2d;
        use crate::nurbs::b_spline_surface::BSplineSurface;
        use crate::nurbs::math::{CPoint, NOW};
        use crate::nurbs::srf_basic::Surface as Net;
        use crate::topology::brep::p2::ArenaNotifyState;
        use crate::topology::brep::p8::{create_face_from_surface, FaceBuildEnv};
        let (x, y, h) = (
            83.555_557_250_976_56,
            87.111_114_501_953_12,
            101.333_335_876_464_84,
        );
        let base = [
            Point3::new(-x, -y, 0.0),
            Point3::new(x, -y, 0.0),
            Point3::new(x, y, 0.0),
            Point3::new(-x, y, 0.0),
        ];
        let apex = Point3::new(0.0, 0.0, h);
        let third = |a: Point3, b: Point3, i: u32| a + (b - a) * (f64::from(i) / 3.0);
        // rows along the first side, columns towards the second
        let patch = |row: &dyn Fn(u32) -> (Point3, Point3)| {
            let net: Vec<Vec<CPoint>> = (0..4)
                .map(|i| {
                    let (a, b) = row(i);
                    (0..4)
                        .map(|j| {
                            let p = third(a, b, j);
                            CPoint::new(p.x, p.y, p.z, NOW)
                        })
                        .collect()
                })
                .collect();
            let knots = vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
            let n = Net::new(net, 3, 3, knots.clone(), knots).unwrap();
            Surface::from(BSplineSurface::from_gw_nurb(&n))
        };
        let mut surfaces = vec![patch(&|i| {
            (third(base[0], base[1], i), third(base[3], base[2], i))
        })];
        for s in 0..4 {
            let (a, b) = (base[s], base[(s + 1) % 4]);
            surfaces.push(patch(&|i| (third(a, b, i), apex)));
        }
        let mut b = crate::topology::brep::p3::new_brep(
            k,
            crate::geom::tol::USE_DEFAULT,
            crate::geom::tol::USE_DEFAULT,
        )
        .unwrap();
        {
            let mut state = ArenaNotifyState::new(&b);
            let mut env = state.env(None);
            let mut marks = crate::topology::topology_traverser::Marks::default();
            let this_face = |f: FaceId| crate::topology::loopuse::notify_ref(EntityRef::Face(f));
            let mut fenv = FaceBuildEnv {
                env: &mut env,
                marks: &mut marks,
                this_face: &this_face,
            };
            let domain = Extent2d::from_bounds(0.0, 0.0, 1.0, 1.0);
            for s in surfaces {
                create_face_from_surface(k, &mut b, Some(s), &domain, &mut fenv).unwrap();
            }
        }
        p2::set_tolerance(k, &mut b, 0.001_013_333_358_765, false, true).unwrap();
        p5::stitch_and_orient(k, &mut b, &p5::StitchAndOrientOptions::default()).unwrap();
        b.refresh_cached_topology();
        b
    }

    /// gen0044 of the difftest unit: a pyramid's side replaced by the plane through its old
    /// edges (its middle point and normal there, so the same plane to rounding). Nothing
    /// moves, and the neighbours' new edges must still cross the extended surfaces'
    /// boundaries: the body is a sound manifold solid with the input's volume, the same
    /// faces, edges and vertices, and its B-spline surfaces' knot vectors are valid. The
    /// side's two edges end at the apex at two parameters (its collapsed row, a pole): the
    /// plane's face has no edge there, where it had one of no length on that face alone.
    #[test]
    fn a_pyramid_side_given_the_plane_through_its_edges_stays_sound() {
        use crate::topology::brep::validation::{check_brep, CheckLevel};
        let mut k = Kernel::default();
        let mut b = bicubic_pyramid(&mut k);
        let (input, _) = volume_and_manifold(&k, &mut b);
        let before = (b.faces().len(), b.edges().len(), b.vertices().len());
        let side = bp1::faces(&b, None)
            .unwrap()
            .into_iter()
            .find(|&f| {
                let s = face_ops::surface(&b, f).unwrap().unwrap();
                let uv = b.face(f).unwrap().uv_domain.evaluate(0.5, 0.5);
                s.evaluate_normal(&k, uv, Side::Left, Side::Left).unwrap().z > 0.0
            })
            .unwrap();
        let s = face_ops::surface(&b, side).unwrap().unwrap();
        let uv = b.face(side).unwrap().uv_domain.evaluate(0.5, 0.5);
        let p = s.evaluate_point(&k, uv).unwrap();
        let n = s.evaluate_normal(&k, uv, Side::Left, Side::Left).unwrap();
        let plane = Surface::from(Plane::from_normal(Vec3::new(p.x, p.y, p.z), n));
        with_notifications(&mut b, |b, notes| {
            local_operation(&mut k, b, &[side], vec![plane], notes)
        })
        .unwrap();
        let (v, manifold) = volume_and_manifold(&k, &mut b);
        assert!(manifold);
        assert!((v - input).abs() <= 1e-9 * input, "{v} vs {input}");
        assert_eq!(
            (b.faces().len(), b.edges().len(), b.vertices().len()),
            before
        );
        b.use_arena_naming().unwrap();
        let lines = check_brep(&mut k, &mut b, CheckLevel::Level2)
            .unwrap()
            .lines();
        assert!(lines.is_empty(), "{lines:?}");
        for (_, f) in b.faces().iter() {
            let net = match f.surface.as_ref().map(Surface::kind) {
                Some(SurfaceKind::BSpline(s)) => s,
                Some(SurfaceKind::Plane(p)) => p.base(),
                _ => continue,
            };
            let n = net.gw_nurb().unwrap();
            for (knots, degree) in [(n.u_knots(), n.p()), (n.v_knots(), n.q())] {
                let ordered = knots.iter().zip(knots.iter().skip(1)).all(|(a, b)| a <= b);
                assert!(ordered, "{knots:?}");
                assert!(
                    knots
                        .chunk_by(|a, b| a == b)
                        .all(|run| run.len() <= degree + 1),
                    "{knots:?} degree {degree}"
                );
            }
        }
    }

    /// gen0091 of the difftest unit: a pyramid's side given the plane through its edges and
    /// the side across its apex from it moved out. The four sides no longer meet in one point:
    /// the apex splits into two vertices and a new edge between the two sides the moved one
    /// parts (with the four-face split on). The new vertex is next to the apex, a pole of the
    /// fourth side's B-spline patch, and on no extension of it (a side's extension narrows to
    /// nothing at the pole): the side, which lies on a plane, continues as that plane. The
    /// body is the convex solid of the five planes, a sound manifold of nine edges and six
    /// vertices.
    #[test]
    fn a_pyramid_side_moved_past_its_apex_splits_it_on_the_planes_beside() {
        use crate::topology::brep::validation::{check_brep, CheckLevel};
        let mut k = Kernel::default();
        k.fixes.offset_closes_loop_gaps = true;
        k.fixes.offset_splits_corner_vertices = true;
        k.fixes.offset_four_face_vertex_split = true;
        let mut b = bicubic_pyramid(&mut k);
        let all = bp1::faces(&b, None).unwrap();
        // the planes of the faces, each as (point, outward unit normal)
        let centre = Vec3::new(0.0, 0.0, 25.0);
        let planes: Vec<(Vec3, Vec3)> = all
            .iter()
            .map(|&f| {
                let s = face_ops::surface(&b, f).unwrap().unwrap();
                let uv = b.face(f).unwrap().uv_domain.evaluate(0.5, 0.5);
                let p = s.evaluate_point(&k, uv).unwrap();
                let n = s.evaluate_normal(&k, uv, Side::Left, Side::Left).unwrap();
                (p, if n.dot(p - centre) > 0.0 { n } else { n * -1.0 })
            })
            .collect();
        // faces: the base, then the sides along y = -y0, x = x0, y = y0, x = -x0
        let (kept, moved, distance) = (2_usize, 1_usize, 0.407_738_046_172_282_6);
        let mut given = planes.clone();
        given[moved].0 = given[moved].0 + given[moved].1 * distance;
        let surfaces = [kept, moved]
            .map(|i| Surface::from(Plane::from_normal(given[i].0, given[i].1)))
            .to_vec();
        with_notifications(&mut b, |b, notes| {
            local_operation(&mut k, b, &[all[kept], all[moved]], surfaces, notes)
        })
        .unwrap();
        let (v, manifold) = volume_and_manifold(&k, &mut b);
        assert!(manifold);
        assert_eq!(
            (b.faces().len(), b.edges().len(), b.vertices().len()),
            (5, 9, 6)
        );
        b.use_arena_naming().unwrap();
        let lines = check_brep(&mut k, &mut b, CheckLevel::Level2)
            .unwrap()
            .lines();
        assert!(lines.is_empty(), "{lines:?}");
        // the convex solid of the five planes: its corners, and its volume from its faces'
        // fans about the base's centre
        let meet = |[i, j, l]: [usize; 3]| {
            let [(p, a), (q, c), (r, e)] = [given[i], given[j], given[l]];
            let det = a.dot(c.cross(e));
            (c.cross(e) * a.dot(p) + e.cross(a) * c.dot(q) + a.cross(c) * e.dot(r)) * (1.0 / det)
        };
        let corner = |i: usize, j: usize| meet([0, i, j]);
        let (b0, b1, b2, b3) = (corner(4, 1), corner(1, 2), corner(2, 3), corner(3, 4));
        let (apex, top) = (meet([2, 3, 4]), meet([1, 2, 4]));
        let faces: [&[Vec3]; 5] = [
            &[b0, b3, b2, b1],
            &[b0, b1, top],
            &[b1, b2, apex, top],
            &[b2, b3, apex],
            &[b3, b0, top, apex],
        ];
        let o = (b0 + b2) * 0.5;
        let want: f64 = faces
            .iter()
            .map(|f| {
                f.iter()
                    .skip(1)
                    .zip(f.iter().skip(2))
                    .map(|(&p, &q)| (f[0] - o).dot((p - o).cross(q - o)) / 6.0)
                    .sum::<f64>()
            })
            .sum();
        assert!((v - want).abs() <= 1e-9 * want, "{v} vs {want}");
    }

    /// `b`'s `faces` given `surfaces`: the offset refusal (`Some`), `b` then exactly as it
    /// came in (its `Debug` form, every surface's cache included, the same before and
    /// after), or `None` when the body is returned. What a returned body is, each test
    /// asserts itself.
    fn refused_unchanged_or_returned(
        k: &mut Kernel,
        b: &mut Brep,
        faces: &[FaceId],
        surfaces: Vec<Surface>,
    ) -> Option<OffsetFailure> {
        let before = format!("{b:?}");
        let r = with_notifications(b, |b, notes| local_operation(k, b, faces, surfaces, notes));
        let e = r.err()?;
        assert_eq!(e.kind(), crate::error::ErrorKind::Unsupported, "{e}");
        let Some(crate::error::Reason::Offset(why)) = e.failure().map(|f| &f.reason) else {
            panic!("not an offset refusal: {e}");
        };
        assert!(
            format!("{b:?}") == before,
            "refused ({why:?}), but the input changed"
        );
        Some(*why)
    }

    /// A refusal leaves everything as it came in: `b` in its `Debug` form, the caches its
    /// surfaces had built before the call included (none dropped, none added), and the
    /// caller's notification records and attributes. Refused by the gate after the offset
    /// ran (a cylinder given a wider coaxial one) and by the reach before anything is built
    /// (a side moved off its round); `b` carries an attribute of its own, lent to the
    /// offset, and the caller's store one of the caller's.
    #[test]
    fn a_refusal_leaves_the_input_its_caches_and_the_callers_stores_as_they_were() {
        use crate::geom::attribute::{Attribute, Behavior};
        use crate::primitive::api_primitives as ap;
        type Make = fn(&mut Kernel) -> (Brep, FaceId, Surface);
        let cases: [(&str, Make); 2] = [
            ("gate", |k| {
                let b = ap::create_cylinder(k, Point3::ZERO, 1.0, 2.0).unwrap();
                let face = curved_face(&b);
                let s = ap::create_cylinder_surface(k, Point3::new(0.0, 0.0, -1.0), 1.5, 4.0);
                (b, face, Surface::from(s.unwrap()))
            }),
            ("reach", |k| {
                let b = cube_rounded_r5(k);
                let side = face_on_x(k, &b, 20.0);
                let s = Plane::from_normal(Vec3::new(20.1, 0.0, 0.0), Vec3::X);
                (b, side, Surface::from(s))
            }),
        ];
        for (name, make) in cases {
            let mut k = Kernel::default();
            let (mut b, face, s) = make(&mut k);
            b.add_attribute(EntityRef::Face(face), Attribute::long(7, 1, Behavior::Copy))
                .unwrap();
            // caches built on `b`'s surfaces before the call (the copy shares their slots)
            validation::soundness(&mut k, &mut b.clone_standing_for_self().unwrap()).unwrap();
            let before = format!("{b:?}");
            let (r, stores) = with_notifications(&mut b, |b, notes| {
                let h = notes
                    .attributes
                    .insert(Attribute::long(9, 2, Behavior::Copy))
                    .unwrap();
                notes
                    .attributes
                    .add_attribute(ObjectKey(1 << 62), h)
                    .unwrap();
                let stores = (format!("{:?}", notes.records), notes.attributes.clone());
                let r = local_operation(&mut k, b, &[face], vec![s], notes);
                let after = (format!("{:?}", notes.records), notes.attributes.clone());
                (r, (stores, after))
            });
            let e = r.expect_err(name);
            assert!(
                matches!(
                    e.failure().map(|f| &f.reason),
                    Some(crate::error::Reason::Offset(
                        OffsetFailure::Unsound { .. } | OffsetFailure::FacesApart { .. }
                    ))
                ),
                "{name}: {e}"
            );
            assert!(format!("{b:?}") == before, "{name}: the input changed");
            assert!(stores.0 .0 == stores.1 .0, "{name}: the records changed");
            assert!(stores.0 .1 == stores.1 .1, "{name}: the attributes changed");
        }
    }

    /// The gate checks the owners against the body's naming, not against the owners the
    /// geometry records itself: an edge's curve owned under an identity no edge has
    /// fails `standing`'s owner test, where `soundness` (each geometry's recorded owner
    /// taken as right) cannot see it.
    #[test]
    fn the_gate_checks_curve_owners_against_the_bodys_naming() {
        let mut k = Kernel::default();
        let mut b = create_box(&mut k, Point3::ZERO, 2.0, 2.0, 2.0).unwrap();
        let measure = |k: &mut Kernel, b: &Brep| {
            let mut copy = b.clone_standing_for_self().unwrap();
            let standing = validation::standing(k, &mut copy).unwrap();
            let mut copy = b.clone_standing_for_self().unwrap();
            (
                standing.check,
                validation::soundness(k, &mut copy).unwrap().check,
            )
        };
        let (clean, clean_own) = measure(&mut k, &b);
        assert!(
            clean.is_valid() && clean_own.is_valid(),
            "{:?}",
            clean.lines()
        );
        // an edge's identity no edge of the body has
        let mut wrong = b.owner_identity(EntityRef::Edge(b.edges().ids().next().unwrap()));
        wrong.key = ObjectKey(wrong.key.0 ^ (1 << 40));
        let e0 = b.edges().ids().next().unwrap();
        let curve = b.edge_mut(e0).unwrap().curve.as_mut().unwrap();
        let _ = curve.set_owner(Some(wrong), None);
        let (named, own) = measure(&mut k, &b);
        assert!(
            !named.failures_beyond(&clean).is_empty(),
            "{:?}",
            named.lines()
        );
        assert!(
            own.failures_beyond(&clean_own).is_empty(),
            "{:?}",
            own.lines()
        );
    }

    /// A body whose pointer check fails is not checked further (`AssertValid` is skipped):
    /// the gate refuses it even where its input failed the pointer check alike, rather
    /// than finding nothing beyond the input's failures.
    #[test]
    fn the_gate_refuses_a_body_whose_geometry_it_cannot_check() {
        let mut k = Kernel::default();
        let mut b = create_box(&mut k, Point3::ZERO, 2.0, 2.0, 2.0).unwrap();
        let e = b.edges().ids().next().unwrap();
        b.edge_mut(e).unwrap().curve = None;
        let input =
            validation::standing(&mut k, &mut b.clone_standing_for_self().unwrap()).unwrap();
        assert!(
            input.check.assert_valid.is_none(),
            "{:?}",
            input.check.lines()
        );
        let e = refuse_unsound(&mut k, &b, &input).unwrap_err();
        assert!(
            matches!(
                e.failure().map(|f| &f.reason),
                Some(crate::error::Reason::Offset(OffsetFailure::Unsound {
                    what: crate::error::offset::Unsoundness::GeometryUnchecked
                }))
            ),
            "{e}"
        );
    }

    /// Review finding 5: a torus face meets itself across seams both ways round; a cylinder
    /// is closed one way only, so one seam would have no edge: refused before anything
    /// changes.
    #[test]
    fn a_torus_given_a_cylinder_is_refused_as_a_seam_opened() {
        use crate::primitive::api_primitives as ap;
        let mut k = Kernel::default();
        let mut b = ap::create_torus(&mut k, Point3::ZERO, 3.0, 1.0).unwrap();
        let face = curved_face(&b);
        let s = Surface::from(
            ap::create_cylinder_surface(&k, Point3::new(0.0, 0.0, -2.0), 3.5, 4.0).unwrap(),
        );
        let why = refused_unchanged_or_returned(&mut k, &mut b, &[face], vec![s]);
        assert!(
            matches!(why, Some(OffsetFailure::SeamOpened { face: f }) if f == face),
            "{why:?}"
        );
    }

    /// Review finding 5: a cylinder's side given a coaxial cylinder of another radius (its
    /// ends' planes then meet the new side in circles the old edges are not on): the body
    /// built is not sound: refused, the first new failure an edge's curve owner (#61) for the
    /// wider, a vertex an edge is not on (#82) for the narrower.
    #[test]
    fn a_cylinder_given_a_coaxial_cylinder_of_another_radius_is_refused_as_not_sound() {
        use crate::primitive::api_primitives as ap;
        for (radius, want) in [
            (1.5, Some(unsound_check(3, 61))),
            (0.5, Some(unsound_check(3, 82))),
        ] {
            let mut k = Kernel::default();
            let mut b = ap::create_cylinder(&mut k, Point3::ZERO, 1.0, 2.0).unwrap();
            let face = curved_face(&b);
            let s = Surface::from(
                ap::create_cylinder_surface(&k, Point3::new(0.0, 0.0, -1.0), radius, 4.0).unwrap(),
            );
            let why = refused_unchanged_or_returned(&mut k, &mut b, &[face], vec![s]);
            assert_eq!(why, want, "{radius}");
        }
    }

    /// Review finding 5: a 10 mm cube's vertical edge rounded r 2, the round given a
    /// coaxial cylinder of another radius: refused as not sound, or sound.
    #[test]
    fn a_round_given_a_coaxial_cylinder_is_refused_as_not_sound_or_apart() {
        use crate::blend::fillet_executive::api::fillet_edges;
        use crate::primitive::api_primitives as ap;
        use crate::topology::{edge, vertex};
        for radius in [2.3, 2.05, 1.8] {
            let mut k = Kernel::default();
            let mut b = create_box(&mut k, Point3::ZERO, 10.0, 10.0, 10.0).unwrap();
            let e = b
                .edges()
                .ids()
                .find(|&e| {
                    let (s, t) = edge::start_end_vertices(&b, e).unwrap();
                    [s, t].into_iter().flatten().all(|v| {
                        let p = vertex::point(&b, v).unwrap();
                        p.x == 10.0 && p.y == 10.0
                    })
                })
                .unwrap();
            fillet_edges(&mut k, &mut b, &[e], 2.0, 1, 1, 1.0).unwrap();
            let round = curved_face(&b);
            let s = Surface::from(
                ap::create_cylinder_surface(&k, Point3::new(8.0, 8.0, -1.0), radius, 12.0).unwrap(),
            );
            let why = refused_unchanged_or_sound(&mut k, &mut b, &[round], vec![s]);
            assert!(
                matches!(why, None | Some(OffsetFailure::Unsound { .. })),
                "{radius}: {why:?}"
            );
        }
    }

    /// gen0091 of the difftest unit with the default fixes: the pyramid's side moved past
    /// its apex needs the four-face vertex split (opt-in), without which every end merges
    /// into one vertex grown to cover the split: not sound, refused.
    #[test]
    fn a_pyramid_side_moved_past_its_apex_is_refused_without_the_four_face_split() {
        let mut k = Kernel::default();
        let mut b = bicubic_pyramid(&mut k);
        let all = bp1::faces(&b, None).unwrap();
        let plane_of = |b: &Brep, k: &Kernel, f: FaceId, d: f64| {
            let s = face_ops::surface(b, f).unwrap().unwrap();
            let uv = b.face(f).unwrap().uv_domain.evaluate(0.5, 0.5);
            let p = s.evaluate_point(k, uv).unwrap();
            let n = s.evaluate_normal(k, uv, Side::Left, Side::Left).unwrap();
            let n = if n.dot(p - Vec3::new(0.0, 0.0, 25.0)) > 0.0 {
                n
            } else {
                n * -1.0
            };
            Surface::from(Plane::from_normal(p + n * d, n))
        };
        // faces: the base, then the sides along y = -y0, x = x0, y = y0, x = -x0
        let surfaces = vec![
            plane_of(&b, &k, all[2], 0.0),
            plane_of(&b, &k, all[1], 0.407_738_046_172_282_6),
        ];
        let why = refused_unchanged_or_returned(&mut k, &mut b, &[all[2], all[1]], surfaces);
        assert_eq!(why, Some(unsound_check(2, 61)));
    }

    /// A 20 mm cube of bilinear B-spline faces whose bottom (z = 0) carries a defect of its
    /// own: its u knot 0.5 twice, over the degree, the middle control points doubled so the
    /// face is whole, and placed off the middle so that the patch is not a plane
    /// parametrised by length. Its two edges along x, its iso lines, keep the knot: they
    /// fail `SmBSplineCurve` test 2 (an internal knot of multiplicity above the degree).
    fn cube_with_a_defective_bottom(k: &mut Kernel) -> Brep {
        use crate::nurbs::b_spline_surface::BSplineSurface;
        use crate::nurbs::math::{CPoint, NOW};
        use crate::nurbs::srf_basic::Surface as Net;
        use crate::topology::brep::p2::ArenaNotifyState;
        use crate::topology::brep::p8::{create_face_from_surface, FaceBuildEnv};
        let c = |x: f64, y: f64, z: f64| CPoint::new(x, y, z, NOW);
        let bilinear = |a: Point3, b: Point3, d: Point3| {
            // corners a, a + (b - a), a + (d - a), rows along b - a
            let e = b + (d - a);
            let net = vec![
                vec![c(a.x, a.y, a.z), c(d.x, d.y, d.z)],
                vec![c(b.x, b.y, b.z), c(e.x, e.y, e.z)],
            ];
            let knots = vec![0.0, 0.0, 1.0, 1.0];
            let n = Net::new(net, 1, 1, knots.clone(), knots).unwrap();
            Surface::from(BSplineSurface::from_gw_nurb(&n))
        };
        let p = Point3::new;
        let s = 20.0;
        let bottom = {
            let row = |x: f64| vec![c(x, 0.0, 0.0), c(x, s, 0.0)];
            let net = vec![row(0.0), row(5.0), row(5.0), row(s)];
            let n = Net::new(
                net,
                1,
                1,
                vec![0.0, 0.0, 0.5, 0.5, 1.0, 1.0],
                vec![0.0, 0.0, 1.0, 1.0],
            )
            .unwrap();
            Surface::from(BSplineSurface::from_gw_nurb(&n))
        };
        let surfaces = vec![
            bottom,
            bilinear(p(0.0, 0.0, s), p(s, 0.0, s), p(0.0, s, s)),
            bilinear(p(0.0, 0.0, 0.0), p(s, 0.0, 0.0), p(0.0, 0.0, s)),
            bilinear(p(0.0, s, 0.0), p(s, s, 0.0), p(0.0, s, s)),
            bilinear(p(0.0, 0.0, 0.0), p(0.0, s, 0.0), p(0.0, 0.0, s)),
            bilinear(p(s, 0.0, 0.0), p(s, s, 0.0), p(s, 0.0, s)),
        ];
        let mut b = crate::topology::brep::p3::new_brep(
            k,
            crate::geom::tol::USE_DEFAULT,
            crate::geom::tol::USE_DEFAULT,
        )
        .unwrap();
        {
            let mut state = ArenaNotifyState::new(&b);
            let mut env = state.env(None);
            let mut marks = crate::topology::topology_traverser::Marks::default();
            let this_face = |f: FaceId| crate::topology::loopuse::notify_ref(EntityRef::Face(f));
            let mut fenv = FaceBuildEnv {
                env: &mut env,
                marks: &mut marks,
                this_face: &this_face,
            };
            let domain = Extent2d::from_bounds(0.0, 0.0, 1.0, 1.0);
            for s in surfaces {
                create_face_from_surface(k, &mut b, Some(s), &domain, &mut fenv).unwrap();
            }
        }
        p5::stitch_and_orient(k, &mut b, &p5::StitchAndOrientOptions::default()).unwrap();
        b.refresh_cached_topology();
        b
    }

    /// The failure classes of `b`'s check ([`validation::soundness`], on a copy).
    fn failure_classes(k: &mut Kernel, b: &Brep) -> BTreeMap<FailureClass, usize> {
        let mut copy = b.clone_standing_for_self().unwrap();
        validation::soundness(k, &mut copy)
            .unwrap()
            .check
            .failure_classes()
    }

    /// The gate holds the operation to what it adds: a cube whose bottom fails its check
    /// on its own, its top raised by 5, is returned (the bottom untouched, its failure
    /// the input's), a manifold solid of the raised volume whose check has no failure
    /// class beyond the input's.
    #[test]
    fn a_failure_the_input_has_on_an_untouched_face_is_not_the_operations() {
        let mut k = Kernel::default();
        let mut b = cube_with_a_defective_bottom(&mut k);
        let had = failure_classes(&mut k, &b);
        assert_eq!(
            had.get(&FailureClass::Report {
                reporting_type: crate::geom::assert_array::labels::B_SPLINE_CURVE_TYPE,
                list_index: 0,
                test_index: 2,
                owner_type: crate::geom::assert_array::labels::B_SPLINE_CURVE_TYPE,
            }),
            Some(&2),
            "{had:?}"
        );
        let (v0, manifold) = volume_and_manifold(&k, &mut b);
        assert!(manifold && (v0 - 8000.0).abs() <= 1e-9 * 8000.0, "{v0}");
        let top = box_top(&k, &b);
        let raised = Surface::from(Plane::from_normal(Vec3::new(0.0, 0.0, 25.0), Vec3::Z));
        with_notifications(&mut b, |b, notes| {
            local_operation(&mut k, b, &[top], vec![raised], notes)
        })
        .unwrap();
        let (v, manifold) = volume_and_manifold(&k, &mut b);
        assert!(manifold && (v - 10000.0).abs() <= 1e-9 * 10000.0, "{v}");
        let has = failure_classes(&mut k, &b);
        assert!(
            has.iter().all(|(c, n)| had.get(c).is_some_and(|m| n <= m)),
            "{has:?} vs {had:?}"
        );
    }

    /// The same cube, its top given a roof with a ridge: a degree-1 patch whose interior
    /// knot is a crease. The body built on it fails tests the input did not (at the time
    /// of writing four lines, the first `SmBrep` test 61, at the moved top's edges): the
    /// operation's, refused as a check failure beyond the input's, the cube as it came in
    /// with the input's own failures.
    #[test]
    fn a_new_failure_on_a_moved_face_of_a_defective_input_is_refused() {
        use crate::nurbs::b_spline_surface::BSplineSurface;
        use crate::nurbs::math::{CPoint, NOW};
        use crate::nurbs::srf_basic::Surface as Net;
        let mut k = Kernel::default();
        let mut b = cube_with_a_defective_bottom(&mut k);
        let had = failure_classes(&mut k, &b);
        assert!(!had.is_empty());
        let top = box_top(&k, &b);
        let row = |x: f64, z: f64| vec![CPoint::new(x, 0.0, z, NOW), CPoint::new(x, 20.0, z, NOW)];
        let net = vec![row(0.0, 25.0), row(10.0, 27.0), row(20.0, 25.0)];
        let roof = Net::new(
            net,
            1,
            1,
            vec![0.0, 0.0, 0.5, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
        )
        .unwrap();
        let roof = Surface::from(BSplineSurface::from_gw_nurb(&roof));
        let why = refused_unchanged_or_returned(&mut k, &mut b, &[top], vec![roof]);
        assert_eq!(why, Some(unsound_check(4, 61)));
        assert_eq!(failure_classes(&mut k, &b), had);
    }

    /// [`BrepCheck::failures_beyond`](validation::BrepCheck::failures_beyond) counts per
    /// class: a class the input has as many lines of is not listed, one more line of it
    /// is, and a class the input lacks is.
    #[test]
    fn failures_beyond_counts_lines_per_class() {
        use crate::geom::assert_array::labels::SURFACE_TYPE;
        use crate::geom::assert_array::{AssertArray, AssertReport, ObjectKey, ReportSpec};
        // (each report about its own surface: the array keeps one of equal reports)
        let report = |test_index: u64, owner: u64| {
            AssertReport::new(ReportSpec {
                test_index,
                test_class: SURFACE_TYPE,
                owner: Some(ObjectKey(owner)),
                owner_type: SURFACE_TYPE,
                ..ReportSpec::default()
            })
            .unwrap()
        };
        let check = |tests: &[u64]| {
            let mut reports = AssertArray::new();
            for (owner, &t) in (0..).zip(tests) {
                assert!(reports.add_unique_delete_duplicate(report(t, owner)));
            }
            validation::BrepCheck {
                pointers_status: crate::error::SM_SUCCESS,
                pointer_reports: AssertArray::new(),
                assert_valid: Some(tests.is_empty()),
                reports,
            }
        };
        let input = check(&[4, 4, 6]);
        assert!(check(&[4, 6]).failures_beyond(&input).is_empty());
        assert!(check(&[4, 4, 6]).failures_beyond(&input).is_empty());
        let more = check(&[4, 4, 4, 6]).failures_beyond(&input);
        assert_eq!(more.len(), 1);
        assert_eq!(more.first().map(|&(_, n)| n), Some(1));
        let other = check(&[5]).failures_beyond(&input);
        assert!(matches!(
            other.as_slice(),
            [(FailureClass::Report { test_index: 5, .. }, 1)]
        ));
    }
}
