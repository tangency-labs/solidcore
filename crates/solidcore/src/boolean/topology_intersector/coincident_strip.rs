//! Not upstream (WR-614, tangency-labs/solidcore#1794; opt-in
//! `Kernel::fixes.boolean_coincident_strip_drops_crossing`): the crossing line of two
//! planes that stay within tolerance of each other across a strip is not an edge.
//!
//! Two planar faces `fa` (B-rep 1) and `fb` (B-rep 2) crossing at a small angle meet in a
//! line fixed only to `tol / sin(angle)`. An edge of one face lying on the other (an
//! edge-on-face edge definition, within tolerance by construction), parallel to the
//! crossing and covering it but farther from it than the tolerance, bounds a strip where
//! the planes are within tolerance of each other (the gap between two planes is linear
//! across it). Upstream still makes the crossing an edge definition. Where the proto edge
//! that holds it also holds that edge-on-face definition, its end proto vertices hold the
//! crossing's and the edge's end vertices, one picked per B-rep
//! (`SmProtoVertex::GetFinalVertex`), which `SetRelationships` never relates, and
//! `SmMerge::MakeFaceuse` fails on the edge's inconsistent vertex mates (`SM_ERR`). The
//! same happens when the strip is bounded on both sides by such edges in proto edges of
//! their own. The fix removes the crossing's edge definition before `Resolve` in either
//! case (and its proto edge when that is left without one), so the strip is bounded by
//! the existing edges.
//!
//! | Rust | Upstream |
//! |------|----------|
//! | `TopologyIntersector::drop_crossings_in_coincident_strips` (crate) | none (own design) |
//!
//! Deviations: not upstream at all; off by default.

use crate::error::Result;
use crate::geom::solution_array::SolutionObject;
use crate::intersect::point_class::{ClassObject, PointClassType};
use crate::kernel::Kernel;
use crate::math::{Point3, Vec3};
use crate::surface::surface::p2::SurfaceKind;
use crate::topology::proto_topology::{EdgeDefId, ProtoEdgeId};

use super::p1::{ignore_soft, TopologyIntersector};
use super::p2::IntersectorBrepP2;

/// An edge definition as the strip test reads it: its proto edge, its source objects and
/// their classification types in B-reps 1 and 2, the ends of its curve over its interval,
/// and its tolerance.
struct Def<O> {
    pe: ProtoEdgeId,
    ed: EdgeDefId,
    src: [(PointClassType, Option<O>); 2],
    ends: (Point3, Point3),
    tol: f64,
}

impl<O: ClassObject + SolutionObject + Ord> TopologyIntersector<O> {
    /// Removes every face/face edge definition of two planes whose crossing lies in a
    /// strip bounded by an edge of one face lying on the other (module docs), and every
    /// proto edge that is left without an edge definition, and returns the curves deleted
    /// with them.
    ///
    /// # Errors
    ///
    /// Those of the curve evaluations and the face/edge connection tests.
    pub(crate) fn drop_crossings_in_coincident_strips(
        &mut self,
        k: &Kernel,
        brep: &impl IntersectorBrepP2<O>,
    ) -> Result<Vec<O>> {
        let Some(mgr) = self.proto_topo_mgr.as_mut() else {
            return Ok(Vec::new());
        };
        let mut defs = Vec::new();
        for &pe in mgr.proto_edges() {
            let Some(edge) = mgr.proto_edge(pe) else {
                continue;
            };
            for &ed in edge.edge_definitions() {
                let Some(def) = mgr.edge_def(ed) else {
                    continue;
                };
                let Some(curve) = def.xsect_curve_3d().and_then(|c| brep.curve(c)) else {
                    continue;
                };
                let ivl = def.xsect_interval();
                defs.push(Def {
                    pe,
                    ed,
                    src: [
                        (def.src_topo_type(1), def.src_topo(1)),
                        (def.src_topo_type(2), def.src_topo(2)),
                    ],
                    ends: (
                        curve.evaluate_point(k, ivl.min())?,
                        curve.evaluate_point(k, ivl.max())?,
                    ),
                    tol: def.tolerance(k),
                });
            }
        }
        let mut doomed: Vec<(ProtoEdgeId, EdgeDefId)> = Vec::new();
        for crossing in &defs {
            if in_coincident_strip(brep, crossing, &defs)? {
                doomed.push((crossing.pe, crossing.ed));
            }
        }
        let mut deleted = Vec::new();
        for (pe, ed) in doomed {
            if let Some(edge) = mgr.proto_edge_mut(pe) {
                ignore_soft(edge.remove_edge_definition(ed))?;
            }
            deleted.extend(mgr.delete_edge_definition(ed));
            if mgr
                .proto_edge(pe)
                .is_some_and(|e| e.num_edge_definitions() == 0)
            {
                ignore_soft(mgr.remove_proto_edge(pe))?;
                mgr.forget_proto_edge(pe);
                deleted.extend(mgr.delete_proto_edge_object(pe));
            }
        }
        Ok(deleted)
    }
}

/// Whether `face` lies on a plane.
fn is_planar<O: ClassObject>(brep: &impl IntersectorBrepP2<O>, face: O) -> bool {
    brep.face_surface(face)
        .flatten()
        .and_then(|s| brep.surface(s))
        .is_some_and(|s| matches!(s.kind(), SurfaceKind::Plane(_)))
}

/// Whether `crossing`, a face/face edge definition of two planes, lies in a strip bounded
/// by an edge of one face lying on the other face: one in the crossing's own proto edge,
/// or one on each side in other proto edges (module docs).
fn in_coincident_strip<O: ClassObject>(
    brep: &impl IntersectorBrepP2<O>,
    crossing: &Def<O>,
    defs: &[Def<O>],
) -> Result<bool> {
    let [(PointClassType::Face, Some(fa)), (PointClassType::Face, Some(fb))] = crossing.src else {
        return Ok(false);
    };
    if !is_planar(brep, fa) || !is_planar(brep, fb) {
        return Ok(false);
    }
    let (p0, p1) = crossing.ends;
    let axis = p1 - p0;
    let len = axis.length();
    if len <= crossing.tol {
        return Ok(false);
    }
    let dir = axis / len;
    // the offsets from the crossing of the edges of `fa` on `fb` and of `fb` on `fa` held
    // by other proto edges
    let mut on_fb = Vec::new();
    let mut on_fa = Vec::new();
    for def in defs.iter().filter(|d| d.ed != crossing.ed) {
        let bounding = match def.src {
            [(PointClassType::Edge, Some(e)), (PointClassType::Face, Some(f))] => {
                f == fb && brep.edge_is_connected_to_face(e, fa)?
            }
            [(PointClassType::Face, Some(f)), (PointClassType::Edge, Some(e))] => {
                f == fa && brep.edge_is_connected_to_face(e, fb)?
            }
            _ => false,
        };
        if !bounding {
            continue;
        }
        let Some(off) = band_offset(p0, dir, len, def, crossing.tol) else {
            continue;
        };
        if def.pe == crossing.pe {
            return Ok(true);
        }
        if matches!(def.src[0].0, PointClassType::Edge) {
            on_fb.push(off);
        } else {
            on_fa.push(off);
        }
    }
    Ok(on_fb.iter().any(|a| on_fa.iter().any(|b| a.dot(*b) < 0.0)))
}

/// The offset from the crossing line (through `p0` along unit `dir`, `len` long) of the
/// edge definition `def`, when it runs parallel to the line within `tol`, covers it and
/// lies farther from it than `tol`; else `None`.
fn band_offset<O>(p0: Point3, dir: Vec3, len: f64, def: &Def<O>, tol: f64) -> Option<Vec3> {
    let (q0, q1) = def.ends;
    let (s0, s1) = ((q0 - p0).dot(dir), (q1 - p0).dot(dir));
    let r0 = (q0 - p0) - dir * s0;
    let r1 = (q1 - p0) - dir * s1;
    let parallel = (r0 - r1).length() <= tol;
    let covers = s0.min(s1) <= tol && s0.max(s1) >= len - tol;
    let off = (r0 + r1) / 2.0;
    (parallel && covers && off.length() > tol).then_some(off)
}
