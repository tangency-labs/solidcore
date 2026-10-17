//! B-rep editing, part 5: coincident-vertex targets and dependency-limited edits.
//!
//! Upstream: `source/SMLib/src/SmBrep.cpp`, `source/SMLib/inc/SmBrep.h`.
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`coincident_vertex_target`], [`combine_coincident_vertices`] | `SmBrep::CombineCoincidentVertices` |
//! | [`create_wire_edge_from_curve`] | `SmBrep::CreateWireEdgeFromCurve` |
//! | [`delete_topological_vertex`] | `SmBrep::DeleteTopologicalVertex` |
//! | [`find_and_set_infinite_region`] | `SmBrep::FindAndSetInfiniteRegion` |
//! | [`glue_vertices`] | `SmBrep::GlueVertices` |
//! | [`make_edge_in_face`] | `SmBrep::MakeEdgeInFace` |
//! | [`make_manifold`] | `SmBrep::MakeManifold` |
//! | [`remove_topological_edges_and_vertices`] | `SmBrep::RemoveTopologicalEdgesAndVertices` |
//! | [`remove_faces`] | `SmBrep::RemoveFaces` |
//! | [`stitch_and_orient`] | `SmBrep::StitchAndOrient` |
//! | [`squeeze_face`] | `SmBrep::SqueezeFace` |
//! | [`transform`] | `SmBrep::Transform` |
//!
//! Deviations: **partial implementation**. Each edit stops at the first missing
//! dependency; continuations behind these barriers are not implemented or certified.
//! `SmBrep::Notify` is owned by `topology.brep_p2`, not this unit. Owning-list
//! traversal and editing state belong to `topology.brep_p1`: arena allocation order
//! is not substituted for list order. Callers supply the actual editing flag where
//! it is tested before any missing call. Pointer inputs are optional typed ids;
//! out-parameters are mutable outputs where failure must clear them. The old-tol
//! target calculation is factored out for direct numerical-stage verification;
//! that verification does not certify vertex deletion or gluing. The internal
//! `SmPendingTopologyDelete` RAII helper belongs to `MakeTopologyFromData`'s pending
//! allocations; it needs no stand-alone arena API (Rust drops owned inputs on error).
//! Attribute detachment and topology insertion remain that reader's responsibility.
//! Debug drawing/assertions are omitted. No recursion is introduced.

use crate::curve::curve::Curve;
use crate::error::{Result, SolidError};
use crate::geom::axis2_placement::Axis2Placement;
use crate::geom::context::{MarkLock, ALL_MARKS};
use crate::geom::extent1d::Extent1d;
use crate::geom::extent3d::Extent3d;
use crate::geom::tol::ZONE_TOL_3D;
use crate::kernel::Kernel;
use crate::math::{max2, Point3, Vec3};
use crate::nurbs::b_spline_curve::BSplineCurve;
use crate::nurbs::b_spline_surface_n_lib::p2::OrientType;
use crate::topology::{vertex, Brep, EdgeId, FaceId, LoopId, RegionId, ShellId, VertexId};

const OP: &str = "topology::brep::p5";

fn missing<T>(symbol: &'static str, unit: &'static str) -> Result<T> {
    Err(SolidError::unported(symbol, unit))
}

fn require_editing(enabled: bool) -> Result<()> {
    if enabled {
        Ok(())
    } else {
        Err(SolidError::sm_err(OP))
    }
}

/// Reads the owning B-rep's editing member once its record support lands. Unlike
/// the explicit-state entry points, opaque existing callers cannot assume TRUE.
pub(crate) fn editing_enabled(_b: &Brep) -> Result<bool> {
    // UNPORTED(topology.brep_p1): SmBrep::m_bEditingEnabled
    missing("SmBrep::m_bEditingEnabled", "topology.brep_p1")
}

/// The position, tolerance and temporary-edge length threshold calculated before
/// `CombineCoincidentVertices` modifies topology (`SM_USE_OLDTOL`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CoincidentVertexTarget {
    /// Position of the vertex with more incident edges; ties keep the survivor.
    pub point: Point3,
    /// The tolerance assigned after combining the vertices.
    pub tolerance: f64,
    /// `dMinDist`, used when deciding whether a common edge is short enough to squeeze.
    pub minimum_edge_length: f64,
}

/// Computes the initial old-tolerance stage of `CombineCoincidentVertices`.
/// Does not edit either vertex. Distances and the strict edge-count/tolerance
/// comparisons retain the native operation order.
pub fn coincident_vertex_target(
    b: &Brep,
    survivor: VertexId,
    deleted: VertexId,
) -> Result<CoincidentVertexTarget> {
    let edges = vertex::edges(b, survivor, None)?;
    let edges2 = vertex::edges(b, deleted, None)?;
    let s = b.vertex(survivor)?;
    let d = b.vertex(deleted)?;
    let mut point = s.point;
    let mut min_dist = point.distance(d.point);
    if min_dist > d.tolerance {
        let doubled = 2.0 * min_dist;
        min_dist = doubled;
    } else {
        min_dist += d.tolerance;
    }
    let mut tolerance = max2(s.tolerance, min_dist);
    if edges2.len() > edges.len() {
        min_dist = point.distance(d.point);
        if min_dist > s.tolerance {
            let doubled = 2.0 * min_dist;
            min_dist = doubled;
        } else {
            min_dist += s.tolerance;
        }
        point = d.point;
        tolerance = max2(d.tolerance, min_dist);
    }
    Ok(CoincidentVertexTarget {
        point,
        tolerance,
        minimum_edge_length: min_dist,
    })
}

/// Combines coincident vertices. Clears the deleted-edge output before reading
/// either vertex, including failure paths. The target stage is available, but
/// successful deletion and gluing remain dependency limited.
pub fn combine_coincident_vertices(
    k: &Kernel,
    b: &mut Brep,
    survivor: VertexId,
    deleted: VertexId,
    deleted_edge: &mut Option<EdgeId>,
) -> Result<()> {
    *deleted_edge = None;
    if survivor == deleted {
        return Err(SolidError::invalid_input(OP, "merge aliases one vertex"));
    }
    let target = coincident_vertex_target(b, survivor, deleted)?;
    let edges = vertex::edges(b, survivor, None)?;
    let edges2 = vertex::edges(b, deleted, None)?;
    if edges2.is_empty() || edges.is_empty() {
        // UNPORTED(topology.brep_p2): SmBrep::Notify
        return missing("SmBrep::Notify", "topology.brep_p2");
    }
    let mut common = None;
    for e in edges {
        if crate::topology::edge::other_vertex(b, e, Some(survivor))? != Some(deleted) {
            continue;
        }
        let eus = &b.edge(e)?.edgeuses;
        let lamina = if let [head, _] = eus.as_slice() {
            crate::topology::edgeuse::p1::is_loop_edgeuse(b, *head)?
        } else {
            false
        };
        if crate::topology::edge::is_manifold(b, e)? || lamina {
            common = Some(e);
            break;
        }
    }
    if let Some(e) = common {
        let record = b.edge(e)?;
        let curve = record
            .curve
            .as_ref()
            .ok_or_else(|| SolidError::invalid_input(OP, "common edge without curve"))?;
        let length = crate::error::ignored(curve.length(k, &record.interval, target.tolerance))?
            .unwrap_or(0.0);
        if length > target.minimum_edge_length {
            // Retain a long common edge; GlueVertices will close it.
        } else {
            *deleted_edge = Some(e);
            // UNPORTED(topology.brep_p3): SmBrep::SqueezeEdge
            return missing("SmBrep::SqueezeEdge", "topology.brep_p3");
        }
        // A long common edge is retained as a closed edge by GlueVertices.
    } else {
        let faces = vertex::faces(b, survivor)?;
        let faces2 = vertex::faces(b, deleted)?;
        if let Some(face) = faces.into_iter().find(|f| faces2.contains(f)) {
            for vu in vertex::vertexuses(b, survivor)? {
                let faceuse = crate::topology::vertexuse::faceuse(b, vu)?;
                let owning_face = match faceuse {
                    Some(fu) => b.faceuse(fu)?.face,
                    None => None,
                };
                if owning_face == Some(face)
                    && matches!(
                        b.vertexuse(vu)?.owner,
                        crate::topology::VertexuseOwner::Loopuse(_)
                    )
                {
                    return Err(SolidError::sm_err(OP));
                }
            }
            let start = b.vertex(survivor)?.point;
            let mut end = b.vertex(deleted)?.point;
            if crate::geom::tol::in_tol(
                start.distance(end),
                crate::geom::tol::scaled_zero2(&start, &end),
            ) {
                end = start + Vec3::new(target.tolerance, target.tolerance, target.tolerance);
            }
            // Preserve failures from the temporary line before the missing edit
            // metadata is consulted. The owned temporary drops on deferral.
            let _line = BSplineCurve::create_line_segment(3, start, end)?;
            // MakeEdgeInFace first tests the editing member, which the arena does
            // not carry yet. Do not assume editing is enabled for an opaque owner.
            // UNPORTED(topology.brep_p1): SmBrep::m_bEditingEnabled
            return missing("SmBrep::m_bEditingEnabled", "topology.brep_p1");
        }
    }
    crate::error::ignore_status(glue_vertices(b, survivor, deleted, None))
}

/// Makes a wire edge over the caller's curve interval. Evaluates start before end,
/// then needs the owning B-rep's infinite region; the owned curve drops on error.
pub fn create_wire_edge_from_curve(
    _k: &Kernel,
    _b: &mut Brep,
    curve: Curve,
    interval: Extent1d,
) -> Result<EdgeId> {
    curve.evaluate_point(interval.min())?;
    curve.evaluate_point(interval.max())?;
    drop(curve);
    // UNPORTED(topology.brep_p1): SmBrep::GetInfiniteRegion
    missing("SmBrep::GetInfiniteRegion", "topology.brep_p1")
}

/// Deletes a degree-two topological vertex, optionally selecting the shorter edge.
/// The native editing-disabled refusal occurs before Notify and input checks.
pub fn delete_topological_vertex(
    _k: &Kernel,
    _b: &mut Brep,
    editing_enabled: bool,
    _vertex: Option<VertexId>,
    _edge: Option<EdgeId>,
) -> Result<()> {
    require_editing(editing_enabled)?;
    // UNPORTED(topology.brep_p2): SmBrep::Notify
    missing("SmBrep::Notify", "topology.brep_p2")
}

/// Finds the infinite region, then resets its void state. The ray-selection
/// operation belongs to part 9; its optional shell and bounding box are retained
/// in the interface for that integration. `void_setting_code` is the native
/// `SmSetRegionVoidType` code (0 nested-solid, 1 infinite-only, 2 swap). Its common
/// enum is being ported in part 1's independent PR; this deferred boundary does
/// not duplicate that type.
pub fn find_and_set_infinite_region(
    _b: &mut Brep,
    _shell: Option<ShellId>,
    _void_setting_code: i32,
    _bbox: Option<&mut Extent3d>,
) -> Result<()> {
    // UNPORTED(topology.brep_p9): SmBrep::FindInfiniteRegion
    missing("SmBrep::FindInfiniteRegion", "topology.brep_p9")
}

/// Glues vertices, leaving the survivor in place. `wire_edges` is the optional
/// caller-supplied wire list, rather than an empty-list substitute for no list.
pub fn glue_vertices(
    _b: &mut Brep,
    _survivor: VertexId,
    _deleted: VertexId,
    _wire_edges: Option<&[EdgeId]>,
) -> Result<()> {
    glue_vertices_entry()
}

/// Common entry barrier for arena and opaque-handle callers. The first native
/// action is the B-rep notification, before dereferencing either vertex; no
/// temporary or substitute B-rep is constructed for handle-only adapters.
pub(crate) fn glue_vertices_entry() -> Result<()> {
    // UNPORTED(topology.brep_p2): SmBrep::Notify
    missing("SmBrep::Notify", "topology.brep_p2")
}

/// Outputs of `MakeEdgeInFace`, initialized to null after its editing and face checks.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct EdgeInFace {
    /// The inserted edge.
    pub edge: Option<EdgeId>,
    /// The new loop when insertion splits an existing loop.
    pub loop_: Option<LoopId>,
    /// The new face when insertion splits the face.
    pub face: Option<FaceId>,
}

/// Finds both vertex sectors then inserts an edge into a face. The native output
/// lifecycle is preserved: editing-disabled/null-face failures leave outputs alone;
/// a valid face clears all three before sector classification.
pub fn make_edge_in_face(
    b: &mut Brep,
    editing_enabled: bool,
    face: Option<FaceId>,
    _start: Option<VertexId>,
    _end: Option<VertexId>,
    _curve: Curve,
    _uv_curve: Option<BSplineCurve>,
    _interval: Extent1d,
    _orientation: OrientType,
    _max_gap_3d: f64,
    output: &mut EdgeInFace,
) -> Result<()> {
    require_editing(editing_enabled)?;
    let f = face.ok_or_else(|| SolidError::null_pointer(OP))?;
    b.face(f)?;
    *output = EdgeInFace::default();
    // UNPORTED(topology.face_p4): SmFace::FindVertexSector
    missing("SmFace::FindVertexSector", "topology.face_p4")
}

/// Makes the selected regions manifold. `None` means all non-infinite regions;
/// `Some(&[])` means no regions and returns the native `SM_ERR`.
pub fn make_manifold(
    _b: &mut Brep,
    regions_to_keep: Option<&[RegionId]>,
    _keep_interior_faces: bool,
) -> Result<()> {
    if let Some(regions) = regions_to_keep {
        if regions.is_empty() {
            return Err(SolidError::sm_err(OP));
        }
        // UNPORTED(topology.brep_p1): SmBrep::GetFaces
        missing("SmBrep::GetFaces", "topology.brep_p1")
    } else {
        // UNPORTED(topology.brep_p1): SmBrep::GetRegions
        missing("SmBrep::GetRegions", "topology.brep_p1")
    }
}

/// Removes redundant edges/vertices, preserving optional keep lists. The native
/// owning-list traversal is required even for apparently empty arena contents.
pub fn remove_topological_edges_and_vertices(
    _k: &Kernel,
    _b: &mut Brep,
    _keep_edges: Option<&[EdgeId]>,
    _keep_vertices: Option<&[VertexId]>,
) -> Result<()> {
    // UNPORTED(topology.brep_p1): SmBrep::GetEdges
    missing("SmBrep::GetEdges", "topology.brep_p1")
}

/// Removes the requested faces then their wire edges and shell vertices.
/// Even an empty list advances and releases one unlocked context traversal mark.
/// Editing-disabled refusal occurs before the mark is acquired.
pub fn remove_faces(
    k: &mut Kernel,
    b: &mut Brep,
    editing_enabled: bool,
    faces: &[FaceId],
) -> Result<()> {
    require_editing(editing_enabled)?;
    let mut mark = MarkLock::new(core::slice::from_mut(&mut k.context), 0, ALL_MARKS)?;
    let result = (|| {
        if let Some(&f) = faces.first() {
            b.face(f)?;
            // UNPORTED(topology.face_p4): SmFace::GetEdges
            return missing("SmFace::GetEdges", "topology.face_p4");
        }
        Ok(())
    })();
    mark.release(core::slice::from_mut(&mut k.context))?;
    result
}

/// Stitching and orientation options with the native defaults.
#[derive(Debug, Clone, PartialEq)]
pub struct StitchAndOrientOptions {
    /// The raw stitching tolerance; nonpositive values select twice the B-rep tolerance.
    pub stitch_tolerance: f64,
    /// Allows formation of spine edges.
    pub allow_spine_edges: bool,
    /// Optional faces to swap UV after import.
    pub faces_to_swap: Option<Vec<FaceId>>,
    /// Overrides the lamina-count-based decision when present, including `Some(false)`.
    pub just_unify: Option<bool>,
}

impl Default for StitchAndOrientOptions {
    fn default() -> Self {
        Self {
            stitch_tolerance: -0.1,
            allow_spine_edges: false,
            faces_to_swap: None,
            just_unify: None,
        }
    }
}

/// Computes the first stitching tolerance exactly as `StitchAndOrient` does.
/// Only a nonpositive requested tolerance is capped at `100 * SM_ZONE_TOL_3D`.
pub fn initial_stitch_tolerance(brep_tolerance: f64, requested: f64) -> f64 {
    let local = 2.0 * brep_tolerance;
    if requested <= 0.0 {
        if local > 100.0 * ZONE_TOL_3D {
            100.0 * ZONE_TOL_3D
        } else {
            local
        }
    } else {
        requested
    }
}

/// Stitches twice if needed, then orients surfaces. No success is manufactured
/// from the statuses ignored upstream when a required callee is missing.
pub fn stitch_and_orient(_b: &mut Brep, _options: &StitchAndOrientOptions) -> Result<()> {
    // UNPORTED(topology.brep_p4): SmBrep::StitchFaces
    missing("SmBrep::StitchFaces", "topology.brep_p4")
}

/// Output lists of `SqueezeFace`. All are reset even when editing is disabled.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct SqueezedFace {
    /// Deleted vertices, as stale identities.
    pub deleted_vertices: Vec<VertexId>,
    /// New edges.
    pub new_edges: Vec<EdgeId>,
    /// Deleted edges, as stale identities.
    pub deleted_edges: Vec<EdgeId>,
}

/// Deletes a degenerate face and stitches its former boundary topology.
pub fn squeeze_face(
    b: &mut Brep,
    editing_enabled: bool,
    face: FaceId,
    _stitch_tolerance: Option<f64>,
    output: &mut SqueezedFace,
) -> Result<()> {
    *output = SqueezedFace::default();
    require_editing(editing_enabled)?;
    b.face(face)?;
    // UNPORTED(topology.face_p4): SmFace::GetVertices
    missing("SmFace::GetVertices", "topology.face_p4")
}

/// Scales about the current origin then applies the placement. The identity exit
/// precedes notifications, geometry conversion and tolerance scaling, even with
/// `scale_tolerances = true`. Nonidentity transformation requires part 2's Notify
/// before any point or geometry is changed; those continuations are uncertified.
pub fn transform(
    _k: &Kernel,
    _b: &mut Brep,
    placement: &Axis2Placement,
    scale: Option<Vec3>,
    _scale_tolerances: bool,
) -> Result<()> {
    if placement.is_identity(crate::math::EFF_ZERO)
        && scale.is_none_or(|s| s.approx_eq(Vec3::new(1.0, 1.0, 1.0)))
    {
        return Ok(());
    }
    // UNPORTED(topology.brep_p2): SmBrep::Notify
    missing("SmBrep::Notify", "topology.brep_p2")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::topology::Vertex;

    #[test]
    fn failure_outputs_follow_native_lifecycle() {
        let mut b = Brep::new();
        let v = b.add_vertex_record(Vertex::default()).unwrap();
        let mut edge = Some(EdgeId::from_index(3));
        assert!(combine_coincident_vertices(&Kernel::default(), &mut b, v, v, &mut edge).is_err());
        assert!(edge.is_none());
        let mut out = SqueezedFace {
            deleted_vertices: vec![v],
            new_edges: vec![EdgeId::from_index(1)],
            deleted_edges: vec![EdgeId::from_index(2)],
        };
        assert!(squeeze_face(&mut b, false, FaceId::from_index(0), None, &mut out).is_err());
        assert_eq!(out, SqueezedFace::default());
    }

    #[test]
    fn identity_precedes_tolerance_update_and_notifications() {
        let mut b = Brep::new();
        let mut k = Kernel::default();
        let v = b
            .add_vertex_record(Vertex {
                point: Vec3::new(1., 2., 3.),
                tolerance: 42.,
                ..Vertex::default()
            })
            .unwrap();
        let original = b.vertex(v).unwrap().clone();
        transform(&k, &mut b, &Axis2Placement::IDENTITY, None, true).unwrap();
        assert_eq!(b.vertex(v).unwrap(), &original);
        assert_eq!(b.vertices().len(), 1);
        assert!(remove_faces(&mut k, &mut b, true, &[]).is_ok());
        assert!(remove_faces(&mut k, &mut b, false, &[]).is_err());
    }

    #[test]
    fn empty_face_removal_preserves_an_existing_mark_lock() {
        use crate::geom::context::MarkType;
        let mut k = Kernel::default();
        let mut b = Brep::new();
        let mut outer = MarkLock::new(core::slice::from_mut(&mut k.context), 0, ALL_MARKS).unwrap();
        let held = outer.mark_type();
        let marks = [
            MarkType::Mark,
            MarkType::Mark2,
            MarkType::Mark3,
            MarkType::MarkIo,
            MarkType::MarkAssert,
        ];
        let before = marks.map(|m| k.context.current_mark(m));
        remove_faces(&mut k, &mut b, true, &[]).unwrap();
        assert!(k.context.is_mark_locked(held));
        let advanced: Vec<_> = marks
            .into_iter()
            .zip(before)
            .filter(|(m, n)| k.context.current_mark(*m) != *n)
            .collect();
        assert_eq!(advanced.len(), 1);
        assert_ne!(advanced[0].0, held);
        assert_eq!(k.context.current_mark(advanced[0].0), advanced[0].1 + 2);
        assert!(!k.context.is_mark_locked(advanced[0].0));
        let before_refusal = marks.map(|m| k.context.current_mark(m));
        assert!(remove_faces(&mut k, &mut b, false, &[]).is_err());
        assert_eq!(marks.map(|m| k.context.current_mark(m)), before_refusal);
        outer
            .release(core::slice::from_mut(&mut k.context))
            .unwrap();
        assert!(!k.context.is_mark_locked(held));
    }
}
