//! Topology-first construction: a body built from the file's own topology.
//!
//! The faces, loops, edges and vertices the file shares are made once each and joined
//! as the file joins them, with no geometric stitching, through the kernel's
//! `BrepConstructor` (`SmBrepConstructor`). The calls are those of
//! [`step_import::script::from_body`]: per shell its region and shell, per face its
//! surface (over the face's STEP parameter box) and loops, every oriented edge joined
//! to the edge an earlier face made, and edge limits at the vertex points.
//!
//! A face whose calls fail is left out: the body is built again without it (the edges
//! and vertices it would have made first are then made by the next face that uses
//! them), up to [`skip_limit`] faces, the first failing face in file order each time.
//! Each left-out face is reported with the call that failed.

use solidcore::geom::extent2d::Extent2d;
use solidcore::geom::vector3d::Point3;
use solidcore::kernel::Kernel;
use solidcore::nurbs::b_spline_surface_n_lib::p2::OrientType;
use solidcore::topology::brep_constructor::{BrepConstructor, VertexNode};
use solidcore::topology::{Brep, EdgeId, FaceId, VertexId};
use solidcore::SolidError;
use step_import::model::{Body, BodyKind, Bound, Face, Shell};

/// A face the topology-first build left out.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FaceFailure {
    /// `#id` of the body.
    pub body: u64,
    /// `#id` of the face.
    pub face: u64,
    /// The call that failed and its error.
    pub reason: String,
}

/// Why a build stopped: the failing face (its position among the faces built, in file
/// order) or none for a call outside any face, and the failing call.
struct Failure {
    face: Option<usize>,
    reason: String,
}

/// The most faces of a body of `faces` the build may leave out: 16, or 1% of a larger
/// body.
pub(crate) fn skip_limit(faces: usize) -> usize {
    16usize.max(faces / 100)
}

/// What [`construct`] built.
pub(crate) struct Constructed {
    /// The kernel it was built with.
    pub kernel: Kernel,
    /// The B-rep.
    pub brep: Brep,
    /// Each built face with its file face's `#id`.
    pub faces: Vec<(FaceId, u64)>,
    /// The faces left out.
    pub failures: Vec<FaceFailure>,
}

/// Builds `body` from its topology, leaving out faces whose calls fail.
///
/// # Errors
///
/// The reason when a call outside any face fails, more than [`skip_limit`] faces fail
/// or no face is built.
pub(crate) fn construct(body: &Body) -> Result<Constructed, String> {
    let mut trial = body.clone();
    let mut failures = Vec::new();
    let limit = skip_limit(body.faces.len());
    loop {
        let mut kernel = Kernel::default();
        let mut faces = Vec::new();
        match build(&mut kernel, &trial, &mut faces) {
            Ok(_) if faces.is_empty() => return Err("no face could be built".to_string()),
            Ok(brep) => {
                return Ok(Constructed {
                    kernel,
                    brep,
                    faces,
                    failures,
                })
            }
            Err(Failure {
                face: Some(position),
                reason,
            }) => {
                let face = remove_face(&mut trial, position)
                    .ok_or_else(|| format!("failing face {position} not found: {reason}"))?;
                failures.push(FaceFailure {
                    body: body.step_id,
                    face,
                    reason,
                });
                if failures.len() > limit {
                    return Err(format!("more than {limit} faces failed"));
                }
            }
            Err(Failure { face: None, reason }) => return Err(reason),
        }
    }
}

/// The faces [`build`] builds, in order: shell by shell, in each shell's face order,
/// skipping missing faces.
fn faces_in_order(body: &Body) -> impl Iterator<Item = (usize, usize, &Face)> + '_ {
    body.shells.iter().enumerate().flat_map(move |(s, shell)| {
        shell
            .faces
            .iter()
            .enumerate()
            .filter_map(move |(slot, &fi)| body.faces.get(fi).map(|face| (s, slot, face)))
    })
}

/// Removes the `position`th face [`build`] builds from its shell; returns its `#id`.
fn remove_face(body: &mut Body, position: usize) -> Option<u64> {
    let (s, slot, id) = faces_in_order(body)
        .nth(position)
        .map(|(s, slot, face)| (s, slot, face.step_id))?;
    body.shells.get_mut(s)?.faces.remove(slot);
    Some(id)
}

fn point(p: [f64; 3]) -> Point3 {
    Point3::new(p[0], p[1], p[2])
}

fn orient(same: bool) -> OrientType {
    if same {
        OrientType::Same
    } else {
        OrientType::Opposite
    }
}

/// `what` failed with `e`.
fn call(what: &'static str) -> impl Fn(SolidError) -> String {
    move |e| format!("{what}: {e}")
}

/// One build of `body` on a new `BrepConstructor`; `faces` gets each built face with
/// its file face's `#id`.
fn build(k: &mut Kernel, body: &Body, faces: &mut Vec<(FaceId, u64)>) -> Result<Brep, Failure> {
    let outside = |reason: String| Failure { face: None, reason };
    let solid = !matches!(body.kind, BodyKind::ShellModel);
    let mut c = BrepConstructor::new();
    c.start_brep(k, solidcore::geom::tol::USE_DEFAULT)
        .map_err(call("StartBrep"))
        .map_err(outside)?;
    let mut edges: Vec<Option<EdgeId>> = vec![None; body.edges.len()];
    let mut vertices: Vec<Option<VertexId>> = vec![None; body.vertices.len()];
    let mut first_region = None;
    let mut position = 0usize;
    for (si, shell) in body.shells.iter().enumerate() {
        let void = solid && si > 0 && matches!(body.kind, BodyKind::SolidWithVoids);
        if solid {
            let region = c
                .start_region(true, if void { first_region } else { None })
                .map_err(call("StartRegion"))
                .map_err(outside)?;
            first_region = first_region.or(Some(region));
        }
        c.start_shell()
            .map_err(call("StartShell"))
            .map_err(outside)?;
        for &fi in &shell.faces {
            let Some(face) = body.faces.get(fi) else {
                continue;
            };
            let made = face_calls(
                k,
                &mut c,
                body,
                shell,
                face,
                void,
                &mut edges,
                &mut vertices,
            )
            .map_err(|reason| Failure {
                face: Some(position),
                reason,
            })?;
            faces.push((made, face.step_id));
            position += 1;
        }
        c.end_shell(k, false)
            .map_err(call("EndShell"))
            .map_err(outside)?;
        if solid {
            c.end_region().map_err(call("EndRegion")).map_err(outside)?;
        }
    }
    c.end_brep().map_err(call("EndBrep")).map_err(outside)
}

/// A face loop ready for B-rep construction: oriented edge uses `(edge_index, forward)`
/// or a single-vertex loop.
pub(crate) enum NormalizedLoop {
    /// `VERTEX_LOOP`: index into [`Body::vertices`].
    Vertex(usize),
    /// Oriented `EDGE_LOOP` uses `(edge_index, forward)`.
    Edges(Vec<(usize, bool)>),
}

/// Orders `face`'s loops with the outer loop first (falling back to the loop of largest
/// 3D bounding-box diagonal when the STEP face omits `FACE_OUTER_BOUND`) and splices
/// loops that touch at a shared vertex so no two loops on the same face share a vertex
/// (SMLib Face Test 39).
pub(crate) fn normalized_face_loops(body: &Body, face: &Face) -> Vec<NormalizedLoop> {
    let mut ordered: Vec<&step_import::model::Loop> = Vec::new();
    if face.loops.iter().any(|l| l.outer) {
        ordered.extend(face.loops.iter().filter(|l| l.outer));
        ordered.extend(face.loops.iter().filter(|l| !l.outer));
    } else if face.loops.len() > 1 {
        let mut best_idx = 0usize;
        let mut best_diag_sq = -1.0_f64;
        for (i, l) in face.loops.iter().enumerate() {
            let diag_sq = loop_bbox_diag_sq(body, l);
            if diag_sq > best_diag_sq {
                best_diag_sq = diag_sq;
                best_idx = i;
            }
        }
        if let Some(outer) = face.loops.get(best_idx) {
            ordered.push(outer);
        }
        for (i, l) in face.loops.iter().enumerate() {
            if i != best_idx {
                ordered.push(l);
            }
        }
    } else {
        ordered.extend(face.loops.iter());
    }

    let mut slots: Vec<Option<NormalizedLoop>> = ordered
        .into_iter()
        .map(|l| {
            Some(match &l.bound {
                Bound::Vertex(v) => NormalizedLoop::Vertex(*v),
                Bound::Edges(uses) => {
                    let oriented_uses: Vec<(usize, bool)> = if l.orientation {
                        uses.iter().map(|u| (u.edge, u.forward)).collect()
                    } else {
                        uses.iter().rev().map(|u| (u.edge, !u.forward)).collect()
                    };
                    NormalizedLoop::Edges(oriented_uses)
                }
            })
        })
        .collect();

    loop {
        let mut found = None;
        'search: for (i, slot_i) in slots.iter().enumerate() {
            let Some(NormalizedLoop::Edges(uses_i)) = slot_i else {
                continue;
            };
            for (j_off, slot_j) in slots.iter().skip(i + 1).enumerate() {
                let j = i + 1 + j_off;
                let Some(NormalizedLoop::Edges(uses_j)) = slot_j else {
                    continue;
                };
                for (ei, &(edge_i, fwd_i)) in uses_i.iter().enumerate() {
                    let Some(e_i) = body.edges.get(edge_i) else {
                        continue;
                    };
                    let vi = if fwd_i { e_i.start } else { e_i.end };
                    for (ej, &(edge_j, fwd_j)) in uses_j.iter().enumerate() {
                        let Some(e_j) = body.edges.get(edge_j) else {
                            continue;
                        };
                        let vj = if fwd_j { e_j.start } else { e_j.end };
                        if vi == vj {
                            found = Some((i, j, ei, ej));
                            break 'search;
                        }
                    }
                }
            }
        }
        let Some((i, j, ei, ej)) = found else {
            break;
        };
        let taken_j = slots.get_mut(j).and_then(Option::take);
        let taken_i = slots.get_mut(i).and_then(Option::take);
        if let (Some(NormalizedLoop::Edges(uses_i)), Some(NormalizedLoop::Edges(uses_j))) =
            (taken_i, taken_j)
        {
            let mut merged = Vec::new();
            merged.extend(uses_i.iter().copied().take(ei));
            merged.extend(uses_j.iter().copied().skip(ej));
            merged.extend(uses_j.iter().copied().take(ej));
            merged.extend(uses_i.iter().copied().skip(ei));
            if let Some(slot) = slots.get_mut(i) {
                *slot = Some(NormalizedLoop::Edges(merged));
            }
        } else {
            break;
        }
    }

    slots.into_iter().flatten().collect()
}

/// Squared 3D bounding-box diagonal of `l`, sampling interior points along each edge so
/// closed single-edge loops (circles, ellipses, closed B-splines) report their true span.
fn loop_bbox_diag_sq(body: &Body, l: &step_import::model::Loop) -> f64 {
    let mut min = Point3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut max = Point3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    let mut update = |p: [f64; 3]| {
        min.x = min.x.min(p[0]);
        min.y = min.y.min(p[1]);
        min.z = min.z.min(p[2]);
        max.x = max.x.max(p[0]);
        max.y = max.y.max(p[1]);
        max.z = max.z.max(p[2]);
    };
    match &l.bound {
        Bound::Vertex(v) => {
            if let Some(vtx) = body.vertices.get(*v) {
                update(vtx.point);
            }
        }
        Bound::Edges(uses) => {
            for u in uses {
                let Some(e) = body.edges.get(u.edge) else {
                    continue;
                };
                if let Some(vtx) = body.vertices.get(e.start) {
                    update(vtx.point);
                }
                if let Some(vtx) = body.vertices.get(e.end) {
                    update(vtx.point);
                }
                for sample in 0_u8..=8 {
                    let s = f64::from(sample) * 0.125;
                    if let Some(pt) = step_import::script::edge_point(e, s) {
                        update(pt);
                    }
                }
            }
        }
    }
    if min.x.is_finite() && max.x.is_finite() {
        let dx = max.x - min.x;
        let dy = max.y - min.y;
        let dz = max.z - min.z;
        dx * dx + dy * dy + dz * dz
    } else {
        0.0
    }
}

/// The calls of one face; returns the face made.
fn face_calls(
    k: &mut Kernel,
    c: &mut BrepConstructor,
    body: &Body,
    shell: &Shell,
    face: &Face,
    void: bool,
    edges: &mut [Option<EdgeId>],
    vertices: &mut [Option<VertexId>],
) -> Result<FaceId, String> {
    let vertex_point = |v: usize| body.vertices.get(v).map_or([0.0; 3], |x| x.point);
    // the face normal (after the shell's orientation) points away from the material
    let normal_is_surface_normal = face.same_sense == shell.orientation;
    let mut surface = crate::make_surface(k, &face.surface, None).map_err(call("face surface"))?;
    if let Some(d) = step_import::script::step_uv_domain(body, face) {
        surface
            .adjust_step_uv_domain(k, &Extent2d::from_bounds(d[0], d[1], d[2], d[3]))
            .map_err(call("AdjustSTEPUVDomain"))?;
    }
    let made = c
        .start_face(orient(normal_is_surface_normal != void))
        .map_err(call("StartFace"))?;
    let domain = surface
        .natural_uv_domain_with_nurb(k)
        .map_err(call("GetNaturalUVDomain"))?;
    c.set_face_surface(surface, &domain)
        .map_err(call("SetFaceSurface"))?;

    for l in normalized_face_loops(body, face) {
        let uses = match l {
            NormalizedLoop::Vertex(v) => {
                c.start_and_end_single_vertex_loop(point(vertex_point(v)))
                    .map_err(call("StartAndEndSingleVertexLoop"))?;
                continue;
            }
            NormalizedLoop::Edges(uses) => uses,
        };
        // the loop runs counter-clockwise about the face's own normal
        c.start_loop(orient(face.same_sense))
            .map_err(call("StartLoop"))?;
        for (ei, forward) in uses {
            let Some(e) = body.edges.get(ei) else {
                continue;
            };
            let o = orient(forward == e.same_sense);
            if let Some(existing) = edges.get(ei).copied().flatten() {
                c.start_edge(k, o, Some(existing), None, None)
                    .map_err(call("StartEdge"))?;
                c.end_edge().map_err(call("EndEdge"))?;
                continue;
            }
            // the edge's vertices in its curve's direction, made at first use
            let (a, b) = if e.same_sense {
                (e.start, e.end)
            } else {
                (e.end, e.start)
            };
            let mut ends = [None, None];
            for (slot, v) in ends.iter_mut().zip([a, b]) {
                if let Some(id) = vertices.get_mut(v) {
                    if id.is_none() {
                        *id = Some(
                            c.start_vertex_of_loop(point(vertex_point(v)))
                                .map_err(call("StartVertexOfLoop"))?,
                        );
                    }
                    *slot = *id;
                }
            }
            let [start, end] = ends;
            let id = c
                .start_edge(k, o, None, start, end)
                .map_err(call("StartEdge"))?;
            if let Some(slot) = edges.get_mut(ei) {
                *slot = Some(id);
            }
            let (lo, hi) = if e.same_sense {
                (e.t_start, e.t_end)
            } else {
                (e.t_end, e.t_start)
            };
            let curve = crate::make_curve(k, &e.curve, Some((lo.min(hi), lo.max(hi))))
                .map_err(call("edge curve"))?;
            c.set_edge_curve(k, Some(curve), false)
                .map_err(call("SetEdgeCurve"))?;
            c.set_limits(
                k,
                &VertexNode::from_point(point(vertex_point(a))),
                &VertexNode::from_point(point(vertex_point(b))),
            )
            .map_err(call("SetLimits"))?;
            c.end_edge().map_err(call("EndEdge"))?;
        }
        c.end_loop(k).map_err(call("EndLoop"))?;
    }
    c.end_face().map_err(call("EndFace"))?;
    Ok(made)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{build_body, BuildMode, BuildOptions};
    use solidcore::geom::attribute::{AttributeStore, ObjectKey};
    use solidcore::mass::integrator::properties::{brep_volume, PropertyAttributes};
    use solidcore::topology::brep::validation::{check_brep, CheckLevel};
    use step_import::model::{Curve, Edge, Frame, Loop, OrientedEdge, Surface, Vertex};

    fn xy_frame(origin: [f64; 3], flip_z: bool) -> Frame {
        if flip_z {
            Frame {
                origin,
                x: [1.0, 0.0, 0.0],
                y: [0.0, -1.0, 0.0],
                z: [0.0, 0.0, -1.0],
            }
        } else {
            Frame {
                origin,
                x: [1.0, 0.0, 0.0],
                y: [0.0, 1.0, 0.0],
                z: [0.0, 0.0, 1.0],
            }
        }
    }

    fn line_edge(step_id: u64, vertices: &[Vertex], start: usize, end: usize) -> Edge {
        let p0 = vertices[start].point;
        let p1 = vertices[end].point;
        let dx = p1[0] - p0[0];
        let dy = p1[1] - p0[1];
        let dz = p1[2] - p0[2];
        let len = (dx * dx + dy * dy + dz * dz).sqrt();
        Edge {
            step_id,
            start,
            end,
            curve: Curve::Line {
                origin: p0,
                dir: [dx / len, dy / len, dz / len],
                magnitude: 1.0,
            },
            same_sense: true,
            t_start: 0.0,
            t_end: len,
        }
    }

    fn measure_volume(brep: &Brep) -> f64 {
        let k = Kernel::default();
        let store = AttributeStore::default();
        let attrs = PropertyAttributes {
            store: &store,
            key_of: &|_| ObjectKey(0),
        };
        let mut copy = brep.clone_standing_for_self().unwrap();
        brep_volume(&k, &mut copy, 1e-7, &attrs).unwrap()
    }

    fn assert_level2_valid(brep: &Brep) {
        let mut k = Kernel::default();
        let mut copy = brep.clone_standing_for_self().unwrap();
        let res = check_brep(&mut k, &mut copy, CheckLevel::Level2).unwrap();
        assert!(
            res.pointer_reports.is_empty(),
            "pointer reports: {:?}",
            res.pointer_reports
        );
        assert_eq!(
            res.assert_valid,
            Some(true),
            "assert_valid reports: {:?}",
            res.reports.reports()
        );
        assert!(
            res.reports.reports().is_empty(),
            "reports: {:?}",
            res.reports.reports()
        );
    }

    #[test]
    fn missing_face_outer_bound_orders_closed_circle_outer_loop_first() {
        // Annular cylinder z in [0, 10], R_out = 20, R_in = 8.
        // Both planar caps omit FACE_OUTER_BOUND (`outer: false` on every loop) and list
        // the small inner circle before the closed single-edge outer circle (`start == end`).
        let r_out = 20.0;
        let r_in = 8.0;
        let h = 10.0;
        let vertices = vec![
            Vertex {
                step_id: 1,
                point: [r_in, 0.0, 0.0],
            },
            Vertex {
                step_id: 2,
                point: [r_out, 0.0, 0.0],
            },
            Vertex {
                step_id: 3,
                point: [r_in, 0.0, h],
            },
            Vertex {
                step_id: 4,
                point: [r_out, 0.0, h],
            },
        ];
        let edges = vec![
            Edge {
                step_id: 10,
                start: 0,
                end: 0,
                curve: Curve::Circle {
                    frame: xy_frame([0.0, 0.0, 0.0], false),
                    radius: r_in,
                },
                same_sense: true,
                t_start: 0.0,
                t_end: 360.0,
            },
            Edge {
                step_id: 11,
                start: 1,
                end: 1,
                curve: Curve::Circle {
                    frame: xy_frame([0.0, 0.0, 0.0], false),
                    radius: r_out,
                },
                same_sense: true,
                t_start: 0.0,
                t_end: 360.0,
            },
            Edge {
                step_id: 12,
                start: 2,
                end: 2,
                curve: Curve::Circle {
                    frame: xy_frame([0.0, 0.0, h], false),
                    radius: r_in,
                },
                same_sense: true,
                t_start: 0.0,
                t_end: 360.0,
            },
            Edge {
                step_id: 13,
                start: 3,
                end: 3,
                curve: Curve::Circle {
                    frame: xy_frame([0.0, 0.0, h], false),
                    radius: r_out,
                },
                same_sense: true,
                t_start: 0.0,
                t_end: 360.0,
            },
            line_edge(14, &vertices, 0, 2),
            line_edge(15, &vertices, 1, 3),
        ];
        let faces = vec![
            // Bottom cap at z = 0 (outward normal -Z): inner loop listed first, outer: false.
            Face {
                step_id: 100,
                surface: Surface::Plane {
                    frame: xy_frame([0.0, 0.0, 0.0], true),
                },
                same_sense: true,
                loops: vec![
                    Loop {
                        outer: false,
                        orientation: true,
                        bound: Bound::Edges(vec![OrientedEdge {
                            edge: 0,
                            forward: true,
                        }]),
                    },
                    Loop {
                        outer: false,
                        orientation: true,
                        bound: Bound::Edges(vec![OrientedEdge {
                            edge: 1,
                            forward: false,
                        }]),
                    },
                ],
            },
            // Top cap at z = h (outward normal +Z): inner loop listed first, outer: false.
            Face {
                step_id: 101,
                surface: Surface::Plane {
                    frame: xy_frame([0.0, 0.0, h], false),
                },
                same_sense: true,
                loops: vec![
                    Loop {
                        outer: false,
                        orientation: true,
                        bound: Bound::Edges(vec![OrientedEdge {
                            edge: 2,
                            forward: false,
                        }]),
                    },
                    Loop {
                        outer: false,
                        orientation: true,
                        bound: Bound::Edges(vec![OrientedEdge {
                            edge: 3,
                            forward: true,
                        }]),
                    },
                ],
            },
            // Outer cylinder wall (outward normal +R) with seam edge 5.
            Face {
                step_id: 102,
                surface: Surface::Cylinder {
                    frame: xy_frame([0.0, 0.0, 0.0], false),
                    radius: r_out,
                },
                same_sense: true,
                loops: vec![Loop {
                    outer: false,
                    orientation: true,
                    bound: Bound::Edges(vec![
                        OrientedEdge {
                            edge: 1,
                            forward: true,
                        },
                        OrientedEdge {
                            edge: 5,
                            forward: true,
                        },
                        OrientedEdge {
                            edge: 3,
                            forward: false,
                        },
                        OrientedEdge {
                            edge: 5,
                            forward: false,
                        },
                    ]),
                }],
            },
            // Inner bore wall (outward normal -R, same_sense: false) with seam edge 4.
            Face {
                step_id: 103,
                surface: Surface::Cylinder {
                    frame: xy_frame([0.0, 0.0, 0.0], false),
                    radius: r_in,
                },
                same_sense: false,
                loops: vec![Loop {
                    outer: false,
                    orientation: true,
                    bound: Bound::Edges(vec![
                        OrientedEdge {
                            edge: 0,
                            forward: false,
                        },
                        OrientedEdge {
                            edge: 4,
                            forward: true,
                        },
                        OrientedEdge {
                            edge: 2,
                            forward: true,
                        },
                        OrientedEdge {
                            edge: 4,
                            forward: false,
                        },
                    ]),
                }],
            },
        ];
        let body = Body {
            step_id: 1,
            name: None,
            kind: BodyKind::ManifoldSolid,
            unit_mm: 1.0,
            uncertainty_mm: Some(1e-6),
            shells: vec![Shell {
                step_id: 200,
                closed: true,
                orientation: true,
                faces: vec![0, 1, 2, 3],
            }],
            faces,
            edges,
            vertices,
        };

        let norm = normalized_face_loops(&body, &body.faces[1]);
        assert_eq!(norm.len(), 2);
        match &norm[0] {
            NormalizedLoop::Edges(uses) => assert_eq!(uses, &[(3, true)]),
            NormalizedLoop::Vertex(_) => panic!("expected edge loop"),
        }

        let expected_vol = core::f64::consts::PI * (r_out * r_out - r_in * r_in) * h;
        for mode in [BuildMode::Default, BuildMode::StitchOnly] {
            let built = build_body(
                &body,
                BuildOptions {
                    mode,
                    ..BuildOptions::default()
                },
            )
            .unwrap();
            assert!(built.face_failures.is_empty());
            assert_level2_valid(&built.brep);
            let vol = measure_volume(&built.brep);
            assert!(
                (vol - expected_vol).abs() < 1e-2,
                "mode {mode:?}: got {vol}, expected {expected_vol}"
            );
        }
    }

    #[test]
    fn touching_loops_sharing_vertex_splice_without_face_test_39() {
        // Stepped 2-manifold solid:
        // - Base [0, 20] x [0, 20] x [0, 5] with a triangular through-hole
        //   (10, 10) - (6, 16) - (14, 16) for z in [0, 5].
        // - Boss [0, 20] x [0, 10] x [5, 10] on y in [0, 10].
        // On the z = 5 step face [0, 20] x [10, 20], the hole's apex (10, 10, 5) lies on
        // the outer boundary y = 10, so the STEP outer loop and inner hole loop share a
        // vertex while every 3D edge of the solid remains strictly 2-manifold.
        let pts = [
            // 0..4: z = 0 outer square
            [0.0, 0.0, 0.0],
            [20.0, 0.0, 0.0],
            [20.0, 20.0, 0.0],
            [0.0, 20.0, 0.0],
            // 4..7: z = 0 triangular hole
            [10.0, 10.0, 0.0],
            [6.0, 16.0, 0.0],
            [14.0, 16.0, 0.0],
            // 7..14: z = 5 vertices (note vertex 12 = (10, 10, 5) shared on step face!)
            [0.0, 0.0, 5.0],   // 7
            [20.0, 0.0, 5.0],  // 8
            [20.0, 10.0, 5.0], // 9
            [20.0, 20.0, 5.0], // 10
            [0.0, 20.0, 5.0],  // 11
            [0.0, 10.0, 5.0],  // 12
            [10.0, 10.0, 5.0], // 13 (shared vertex!)
            [6.0, 16.0, 5.0],  // 14
            [14.0, 16.0, 5.0], // 15
            // 16..20: z = 10 boss top
            [0.0, 0.0, 10.0],   // 16
            [20.0, 0.0, 10.0],  // 17
            [20.0, 10.0, 10.0], // 18
            [0.0, 10.0, 10.0],  // 19
        ];
        let vertices: Vec<Vertex> = pts
            .iter()
            .enumerate()
            .map(|(i, &point)| Vertex {
                step_id: (i + 1) as u64,
                point,
            })
            .collect();

        let edge_pairs = [
            // 0..4: z = 0 outer CCW (+Z)
            (0, 1), // 0
            (1, 2), // 1
            (2, 3), // 2
            (3, 0), // 3
            // 4..7: z = 0 hole CCW (+Z): 4->6->5->4
            (4, 6), // 4
            (6, 5), // 5
            (5, 4), // 6
            // 7..13: z = 5 perimeter segments
            (7, 8),   // 7: y=0
            (8, 9),   // 8: x=20, y in [0,10]
            (9, 10),  // 9: x=20, y in [10,20]
            (10, 11), // 10: y=20
            (11, 12), // 11: x=0, y in [10,20]
            (12, 7),  // 12: x=0, y in [0,10]
            // 13..15: z = 5 step split line y=10 via shared vertex 13=(10,10,5)
            (12, 13), // 13: (0,10,5)->(10,10,5)
            (13, 9),  // 14: (10,10,5)->(20,10,5)
            // 15..18: z = 5 hole CW (+Z): 13->14->15->13
            (13, 14), // 15: (10,10,5)->(6,16,5)
            (14, 15), // 16: (6,16,5)->(14,16,5)
            (15, 13), // 17: (14,16,5)->(10,10,5)
            // 18..22: base verticals
            (0, 7),  // 18
            (1, 8),  // 19
            (2, 10), // 20
            (3, 11), // 21
            // 22..25: hole verticals
            (4, 13), // 22
            (5, 14), // 23
            (6, 15), // 24
            // 25..29: boss verticals
            (7, 16),  // 25
            (8, 17),  // 26
            (9, 18),  // 27
            (12, 19), // 28
            // 29..33: z = 10 boss top CCW (+Z)
            (16, 17), // 29
            (17, 18), // 30
            (18, 19), // 31
            (19, 16), // 32
        ];
        let edges: Vec<Edge> = edge_pairs
            .iter()
            .enumerate()
            .map(|(i, &(s, e))| line_edge((100 + i) as u64, &vertices, s, e))
            .collect();

        let oe = |edge: usize, forward: bool| OrientedEdge { edge, forward };
        let lp = |outer: bool, uses: Vec<OrientedEdge>| Loop {
            outer,
            orientation: true,
            bound: Bound::Edges(uses),
        };
        let plane = |origin: [f64; 3], z: [f64; 3], x: [f64; 3]| {
            let y = [
                z[1] * x[2] - z[2] * x[1],
                z[2] * x[0] - z[0] * x[2],
                z[0] * x[1] - z[1] * x[0],
            ];
            Surface::Plane {
                frame: Frame { origin, x, y, z },
            }
        };

        let faces = vec![
            // 0: Bottom z = 0 (normal -Z): outer CW from +Z = 0<-3<-2<-1, hole CCW from +Z = 4->6->5
            Face {
                step_id: 200,
                surface: plane([0.0, 0.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]),
                same_sense: true,
                loops: vec![
                    lp(
                        true,
                        vec![oe(3, false), oe(2, false), oe(1, false), oe(0, false)],
                    ),
                    lp(false, vec![oe(4, true), oe(5, true), oe(6, true)]),
                ],
            },
            // 1: Step face z = 5 (normal +Z): TWO loops sharing vertex 13!
            Face {
                step_id: 201,
                surface: plane([0.0, 0.0, 5.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
                same_sense: true,
                loops: vec![
                    lp(
                        true,
                        vec![
                            oe(13, true),
                            oe(14, true),
                            oe(9, true),
                            oe(10, true),
                            oe(11, true),
                        ],
                    ),
                    lp(false, vec![oe(15, true), oe(16, true), oe(17, true)]),
                ],
            },
            // 2: Boss top z = 10 (normal +Z)
            Face {
                step_id: 202,
                surface: plane([0.0, 0.0, 10.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
                same_sense: true,
                loops: vec![lp(
                    true,
                    vec![oe(29, true), oe(30, true), oe(31, true), oe(32, true)],
                )],
            },
            // 3: Front y = 0, z in [0, 10] (normal -Y, x = +X => y_axis = -Z; CCW about -Y: (0,0,0)->(20,0,0)->(20,0,5)->(20,0,10)->(0,0,10)->(0,0,5)->(0,0,0) is +X then +Z => normal (+X)x(+Z) = -Y!)
            Face {
                step_id: 203,
                surface: plane([0.0, 0.0, 0.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]),
                same_sense: true,
                loops: vec![lp(
                    true,
                    vec![
                        oe(0, true),
                        oe(19, true),
                        oe(26, true),
                        oe(29, false),
                        oe(25, false),
                        oe(18, false),
                    ],
                )],
            },
            // 4: Back y = 20, z in [0, 5] (normal +Y: (-X)x(+Z) = +Y)
            Face {
                step_id: 204,
                surface: plane([0.0, 20.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]),
                same_sense: true,
                loops: vec![lp(
                    true,
                    vec![oe(2, true), oe(21, true), oe(10, false), oe(20, false)],
                )],
            },
            // 5: Boss back wall y = 10, z in [5, 10] (normal +Y: (-X)x(+Z) = +Y)
            Face {
                step_id: 205,
                surface: plane([0.0, 10.0, 5.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]),
                same_sense: true,
                loops: vec![lp(
                    true,
                    vec![
                        oe(14, false),
                        oe(13, false),
                        oe(28, true),
                        oe(31, false),
                        oe(27, false),
                    ],
                )],
            },
            // 6: Right x = 20 (normal +X: (+Y)x(+Z) = +X)
            Face {
                step_id: 206,
                surface: plane([20.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
                same_sense: true,
                loops: vec![lp(
                    true,
                    vec![
                        oe(1, true),
                        oe(20, true),
                        oe(9, false),
                        oe(27, true),
                        oe(30, false),
                        oe(26, false),
                        oe(19, false),
                    ],
                )],
            },
            // 7: Left x = 0 (normal -X: (-Y)x(+Z) = -X)
            Face {
                step_id: 207,
                surface: plane([0.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]),
                same_sense: true,
                loops: vec![lp(
                    true,
                    vec![
                        oe(3, true),
                        oe(18, true),
                        oe(25, true),
                        oe(32, false),
                        oe(28, false),
                        oe(11, false),
                        oe(21, false),
                    ],
                )],
            },
            // 8: Hole wall 4->13->15->6 (bottom edge 4=(10,10)->6=(14,16), d=(4,6,0); inward to hole is (-6,4,0): (+d)x(+Z) = (6,-4,0) points into solid? Wait: hole interior is around (10,14), edge (10,10)->(14,16) has left side (-6,4) toward (10,14), so outward from solid into hole is (-6,4,0) = (+Z)x(d)! So loop (-d) then (+Z) has (-d)x(+Z) = (-6,4,0)!)
            Face {
                step_id: 208,
                surface: plane(
                    [10.0, 10.0, 0.0],
                    [-6.0 / 52.0_f64.sqrt(), 4.0 / 52.0_f64.sqrt(), 0.0],
                    [-4.0 / 52.0_f64.sqrt(), -6.0 / 52.0_f64.sqrt(), 0.0],
                ),
                same_sense: true,
                loops: vec![lp(
                    true,
                    vec![oe(4, false), oe(22, true), oe(17, false), oe(24, false)],
                )],
            },
            // 9: Hole wall 6=(14,16)->5=(6,16), d=(-8,0,0), outward from solid into hole is (0,-1,0) = (-d)x(+Z) where -d=(+8,0,0), (+8,0,0)x(0,0,1) = (0,-8,0)!
            Face {
                step_id: 209,
                surface: plane([14.0, 16.0, 0.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]),
                same_sense: true,
                loops: vec![lp(
                    true,
                    vec![oe(5, false), oe(24, true), oe(16, false), oe(23, false)],
                )],
            },
            // 10: Hole wall 5=(6,16)->4=(10,10), d=(4,-6,0), hole interior is right of d? Left of (6,16)->(10,10) is (6,4) toward (10,14)! So (-d)x(+Z) = (-4,6,0)x(0,0,1) = (6,4,0) points into hole!
            Face {
                step_id: 210,
                surface: plane(
                    [6.0, 16.0, 0.0],
                    [6.0 / 52.0_f64.sqrt(), 4.0 / 52.0_f64.sqrt(), 0.0],
                    [-4.0 / 52.0_f64.sqrt(), 6.0 / 52.0_f64.sqrt(), 0.0],
                ),
                same_sense: true,
                loops: vec![lp(
                    true,
                    vec![oe(6, false), oe(23, true), oe(15, false), oe(22, false)],
                )],
            },
        ];

        let body = Body {
            step_id: 2,
            name: None,
            kind: BodyKind::ManifoldSolid,
            unit_mm: 1.0,
            uncertainty_mm: Some(1e-6),
            shells: vec![Shell {
                step_id: 300,
                closed: true,
                orientation: true,
                faces: (0..faces.len()).collect(),
            }],
            faces,
            edges,
            vertices,
        };

        // Step face (index 1) had 2 loops sharing vertex 13; splicing merges them into 1 loop of 8 edges.
        let spliced = normalized_face_loops(&body, &body.faces[1]);
        assert_eq!(spliced.len(), 1);
        match &spliced[0] {
            NormalizedLoop::Edges(uses) => assert_eq!(uses.len(), 8),
            NormalizedLoop::Vertex(_) => panic!("expected merged edge loop"),
        }

        let expected_vol = 20.0 * 20.0 * 5.0 + 20.0 * 10.0 * 5.0 - 0.5 * 8.0 * 6.0 * 5.0; // 2880.0
        let built = build_body(&body, BuildOptions::default()).unwrap();
        assert!(built.face_failures.is_empty());
        assert_level2_valid(&built.brep);
        let vol = measure_volume(&built.brep);
        assert!(
            (vol - expected_vol).abs() < 1e-6,
            "got {vol}, expected {expected_vol}"
        );
    }
}
