//! The blends: rounds and chamfers, blends of each kind with one size per edge, corner
//! setbacks (closed, or left open for the caller to fill), the variable round and the
//! unequal chamfer.
//!
//! Source: tangency-labs/smlib `capi/src/smlib_ops.cpp` (Apache-2.0, Copyright SolidCraft
//! contributors), replay/1.0.1 at `d51a723` (the file's latest commit there). Selection
//! by points, `outer_faceuse` and `finish` are the crate's (`pick`, `session`).
//!
//! | Rust | C API |
//! |------|-------|
//! | [`Smlib::fillet`] | `smlib_fillet` |
//! | [`Smlib::chamfer`] | `smlib_chamfer` |
//! | [`Smlib::blend`] | `smlib_blend` |
//! | [`Smlib::blend_setback`] | `smlib_blend_setback` |
//! | [`Smlib::blend_setback_open`] | `smlib_blend_setback_open` |
//! | [`Smlib::fillet_variable`] | `smlib_fillet_variable` |
//! | [`Smlib::chamfer_sides`] | `smlib_chamfer_sides` |
//! | `blend_call` (private) | `blend_call` |
//! | `round_edges`, `chamfer_offsets`, `chamfer_edges` (private) | the same |
//! | `chamfer_radius`, `sides_offsets`, `face_frame` (private) | the same |
//! | `blend_failure`, `setbacks_fit`, `open_corners` (private) | the same |
//! | `law_ok`, `law_pinches`, `law_along`, `law_corner` (private) | `law_ok`, `law_pinches`, `law_along`; `law_corner` is new |
//!
//! Deviations:
//! - **Kernel.** The calls run on the API's kernel; the C API runs the fork's SMLib, whose
//!   fixes solidcore carries as opt-in `Kernel::fixes` flags (the SolidCraft preset).
//! - **New faces.** The C API marks every face before the operation with an attribute that
//!   SMLib copies through splits and merges, and counts a face without one as made by it
//!   (`made_faces`: none made is `SMLIB_ERR_NO_EFFECT`). Here a face is made by the
//!   operation when the body had no face of its id before: an arena id is never reused.
//!   A face SMLib only split also gets a new id, so a split alone counts as an effect; a
//!   blend that succeeds always makes its blend faces, so the status is the same.
//! - **Radius corners.** The fork's `SmFilletExecutive::SplitFilletFacesAtRadiusCorners`
//!   (tangency-labs/smlib 6eecf04 and 5fbffc9) is not in solidcore yet. A round without a
//!   radius law has no radius corner, so `round_edges` skips the step. A variable round
//!   whose law has a corner inside the edge (an interior knot repeated `degree` times or
//!   more, where its slope may jump) is refused with `SMLIB_ERR_UNSUPPORTED` rather than
//!   built unsplit, which would be another body than the C API's.
//! - **Not exact.** The fork's `SM_FILERR_NOT_EXACT` comes only from its cliff-overflow
//!   blends (tangency-labs/smlib 8ba4667, tangency-labs/solidcore#960), which solidcore does
//!   not build: `blend_failure` maps `NO_ROOM` and `CORNER_UNSUPPORTED` only.
//! - **Ease.** A setback's ease is 0 (`SetEaseDist` is not in solidcore); the eased entry
//!   points are not ported.
//! - **Corner patches.** The faces of a set-back corner's patch are read from the session
//!   (its `SM_AI_FILLET_CORNER_FACE` marks) before the stitch, as face ids, and must still
//!   be faces of the body after it (else `SMLIB_ERR_UNSUPPORTED`); the C API reads the marks
//!   after the stitch and then removes them. The marks live in the session and go with it.
//! - **Why.** The fillet executive's error comment (`g_blend_why`, which the C API dumps as
//!   the call's `why` input) is not kept: the status and [`Smlib::last_error`] say why.
//! - **Outputs.** `out_failed_point` and `out_failed_corner` are read after the call with
//!   [`Smlib::last_failed_point`] and [`Smlib::last_failed_corner`]; the holes of
//!   [`Smlib::blend_setback_open`] are returned with the body (no `SMLIB_ERR_BUFFER`).
//! - **Recursion.** None.

use std::collections::BTreeSet;

use solidcore::blend::fillet_executive::session::FilletSession;
use solidcore::blend::fillet_geom::p1::{FilletCornerType, FilletErrorType};
use solidcore::blend::fillet_solver::p1::FilletSurfaceGeneratorType;
use solidcore::blend::fillet_standard_solver::p1::FilletLaw;
use solidcore::geom::extent1d::Extent1d;
use solidcore::mass::integrator::properties::edge_length;
use solidcore::math::{Point3, Vec3, PI};
use solidcore::nurbs::b_spline_curve::p1::{
    BSplineCurve, BSplineCurveForm, CanonicalData, KnotType,
};
use solidcore::nurbs::b_spline_surface_n_lib::p2::OrientType;
use solidcore::nurbs::knots::Side;
use solidcore::query::api_queries::brep_is_manifold_solid;
use solidcore::topology::brep::p1;
use solidcore::topology::brep::p2::ArenaNotifyState;
use solidcore::topology::brep::p3::{delete_face, DeleteFaceOptions};
use solidcore::topology::brep::p5::{stitch_and_orient, StitchAndOrientOptions};
use solidcore::topology::edgeuse::p1 as eu1;
use solidcore::topology::edgeuse::p3::{evaluate_binormal_uv, BinormalOptions};
use solidcore::topology::topology_traverser::Marks;
use solidcore::topology::{edge, face, vertex, Brep, EdgeId, EdgeuseId, FaceId, VertexId};
use solidcore::{EntityRef, Kernel};

use crate::pick::{outer_faceuse, pick_edges, pick_vertices};
use crate::session::{finish, finite, num, vec, Fail, Smlib};
use crate::status::{Result, Status};

/// `SMLIB_BLEND_ROUND`: a circular cross-section at constant radius (as [`Smlib::fillet`]).
pub const BLEND_ROUND: i32 = 0;
/// `SMLIB_BLEND_CURVATURE`: a curvature-continuous (G2) blend-curve cross-section between
/// the contact lines the round of that radius has.
pub const BLEND_CURVATURE: i32 = 1;
/// `SMLIB_BLEND_CHORD`: a circular cross-section whose chord (width from contact line to
/// contact line) is constant; the size is that width.
pub const BLEND_CHORD: i32 = 2;
/// `SMLIB_BLEND_CHAMFER`: an equal-distance chamfer (as [`Smlib::chamfer`]); the size is
/// the setback.
pub const BLEND_CHAMFER: i32 = 3;

/// `MAX_POINTS`: the most points one call takes.
const MAX_POINTS: usize = 4096;

/// A corner set back: its vertex and the setback along each of its edges (the ease is 0).
#[derive(Debug, Clone, Copy)]
struct SetBack {
    vertex: VertexId,
    distance: f64,
}

/// One edge round an opened corner's hole: the face beside it, and whether that face's
/// boundary, run counter-clockwise about its outward normal, runs the edge from its start
/// vertex to its end.
#[derive(Debug, Clone, Copy)]
struct HoleEdge {
    edge: EdgeId,
    face: FaceId,
    forward: bool,
}

/// How a round fails: a status, or a setback corner SMLib does not set back (its index in
/// the setbacks, for `out_failed_corner`).
enum RoundFail {
    Fail(Fail),
    Refused(usize),
}

impl From<Fail> for RoundFail {
    fn from(f: Fail) -> Self {
        Self::Fail(f)
    }
}

impl From<solidcore::SolidError> for RoundFail {
    fn from(e: solidcore::SolidError) -> Self {
        Self::Fail(Fail::Solid(e))
    }
}

impl From<Status> for RoundFail {
    fn from(s: Status) -> Self {
        Self::Fail(Fail::Status(s))
    }
}

/// `size_ok`.
fn size_ok(x: f64) -> bool {
    x.is_finite() && x > 0.0 && x <= 1e6
}

/// `tolerance_ok`.
fn tolerance_ok(t: f64) -> bool {
    t.is_finite() && t > 0.0 && t < 1e6
}

/// `points_ok`: 1 to [`MAX_POINTS`] finite points (0 too with `allow_none`).
fn points_ok(points: &[[f64; 3]], allow_none: bool) -> bool {
    if points.len() > MAX_POINTS || (points.is_empty() && !allow_none) {
        return false;
    }
    points.iter().all(|p| finite(p))
}

/// The points, three numbers each, for the log line.
fn pts(points: &[[f64; 3]]) -> String {
    let flat: Vec<f64> = points.iter().flatten().copied().collect();
    if flat.is_empty() {
        "[]".to_string()
    } else {
        vec(&flat)
    }
}

/// The binormal at the middle of edge use `eu`'s edge (`EvaluateBinormal(t, FALSE, ...)`),
/// with the surface parameters there.
fn binormal_at_middle(
    k: &Kernel,
    b: &mut Brep,
    eu: EdgeuseId,
) -> solidcore::Result<Option<(Vec3, solidcore::math::Point2)>> {
    let Some(e) = eu1::edge(b, eu)? else {
        return Ok(None);
    };
    let t = edge::interval(b, e)?.evaluate(0.5);
    match evaluate_binormal_uv(k, b, eu, t, false, &BinormalOptions::default()) {
        Ok(r) => Ok(Some((r.eval.binormal, r.eval.uv))),
        Err(e) if e.is_unported() => Err(e),
        Err(_) => Ok(None),
    }
}

/// `chamfer_radius`: the radius of the ball whose tangent points on the two faces of `e`
/// lie `distance` back from the edge at its midpoint, `distance * tan(a / 2)` with `a` the
/// angle between the directions into the faces.
fn chamfer_radius(
    k: &Kernel,
    b: &mut Brep,
    e: EdgeId,
    distance: f64,
) -> solidcore::Result<Option<f64>> {
    let Some(u0) = edge::blend_edgeuse(k, b, e)? else {
        return Ok(None);
    };
    let Some(u1) = eu1::radial(b, u0)? else {
        return Ok(None);
    };
    let (Some((b0, _)), Some((b1, _))) =
        (binormal_at_middle(k, b, u0)?, binormal_at_middle(k, b, u1)?)
    else {
        return Ok(None);
    };
    let c = b0.dot(b1).clamp(-1.0, 1.0);
    let radius = distance * ((1.0 - c) / (1.0 + c)).sqrt();
    Ok((radius.is_finite() && radius > 0.0).then_some(radius))
}

/// `blend_failure`: the status a blend that failed with `f` reports, by the fillet
/// executive's error code: `NO_ROOM` is `SMLIB_ERR_NO_ROOM` and `CORNER_UNSUPPORTED`
/// `SMLIB_ERR_CORNER_UNSUPPORTED`; any other failure is `f` (`NOT_EXACT`: see the module's
/// deviations).
fn blend_failure(session: &FilletSession<'_>, f: solidcore::SolidError) -> Fail {
    let code = session.executive().error.code;
    if code == FilletErrorType::NoRoom as i32 {
        Fail::Status(Status::ERR_NO_ROOM)
    } else if code == FilletErrorType::CornerUnsupported as i32 {
        Fail::Status(Status::ERR_CORNER_UNSUPPORTED)
    } else {
        Fail::Solid(f)
    }
}

/// The faces `b` has now and had not in `before`, after its stitch: `made_faces`
/// (`SMLIB_ERR_NO_EFFECT` when there are none; see the module's deviations).
fn made_faces(b: &Brep, before: &BTreeSet<FaceId>) -> core::result::Result<(), Fail> {
    if b.face_list().iter().any(|f| !before.contains(f)) {
        Ok(())
    } else {
        Err(Fail::Status(Status::ERR_NO_EFFECT))
    }
}

/// `chamfer_offsets`: chamfers `edges` of `b` in place, edge `i` with its blend edge use's
/// surface offset `first[i]` and the other's `second[i]`: a constant-radius solver with the
/// second radius and a linear cross-section per edge (the fork's `SetFilletParameters` with
/// `pdOptSecondRadius`, approximation tolerance 4 times the edge's, 20 and 10 degrees),
/// bevelled corners, then the stitch.
fn chamfer_offsets(
    k: &mut Kernel,
    b: &mut Brep,
    edges: &[EdgeId],
    first: &[f64],
    second: &[f64],
) -> core::result::Result<(), Fail> {
    let before: BTreeSet<FaceId> = b.face_list().iter().copied().collect();
    let mut tolerances = Vec::new();
    for &e in edges {
        tolerances.push(edge::tolerance(b, e)?);
    }
    {
        let mut session = FilletSession::new(k, b)?;
        let generator = session.add_surface_generator(
            solidcore::blend::fillet_solver::p1::FilletSurfaceGenerator::Linear,
        )?;
        for (((&e, &r0), &r1), &tol) in edges.iter().zip(first).zip(second).zip(&tolerances) {
            let Some(eu) = session.blend_edgeuse(k, e)? else {
                // SetFilletParameters: an edge that cannot be filleted
                return Err(Fail::Status(Status::SM_ERR_INVALID_INPUT));
            };
            let solver = session
                .new_constant_radius_solver_with_second(
                    k,
                    eu,
                    r0,
                    r1,
                    4.0 * tol,
                    20.0 * (PI / 180.0),
                    10.0 * (PI / 180.0),
                )
                .map_err(|e| blend_failure(&session, e))?;
            session
                .load_solver(k, solver)
                .map_err(|e| blend_failure(&session, e))?;
            session.set_solver_generator(solver, generator)?;
        }
        session
            .create_corners(k, None)
            .map_err(|e| blend_failure(&session, e))?;
        for corner in session.fillet_corners()? {
            session.set_corner_bevel(corner, true)?;
        }
        session
            .do_filleting(k, None)
            .map_err(|e| blend_failure(&session, e))?;
    }
    stitch_and_orient(k, b, &StitchAndOrientOptions::default())?;
    // A chamfer that made no face did nothing, whatever SMLib's status said.
    made_faces(b, &before)
}

/// `chamfer_edges`: chamfers `edges` of `b` in place, edge `i` `distances[i]` back along
/// each face.
fn chamfer_edges(
    k: &mut Kernel,
    b: &mut Brep,
    edges: &[EdgeId],
    distances: &[f64],
) -> core::result::Result<(), Fail> {
    let mut radii = Vec::new();
    for (&e, &d) in edges.iter().zip(distances) {
        let Some(r) = chamfer_radius(k, b, e, d)? else {
            return Err(Fail::Status(Status::SM_ERR_INVALID_INPUT));
        };
        radii.push(r);
    }
    chamfer_offsets(k, b, edges, &radii, &radii)
}

/// `face_frame`: the outward unit normal, at the middle of `u`'s edge, of the face `u`
/// bounds, and the unit direction into that face square to the edge.
fn face_frame(k: &Kernel, b: &mut Brep, u: EdgeuseId) -> solidcore::Result<Option<(Vec3, Vec3)>> {
    let Some(f) = eu1::face(b, u)? else {
        return Ok(None);
    };
    let Some(fu) = outer_faceuse(b, f)? else {
        return Ok(None);
    };
    let Some((into, uv)) = binormal_at_middle(k, b, u)? else {
        return Ok(None);
    };
    let Some(surface) = face::p2::surface(b, f)? else {
        return Ok(None);
    };
    let normal = match surface.evaluate_normal(k, uv, Side::Left, Side::Left) {
        Ok(n) => n,
        Err(e) if e.is_unported() => return Err(e),
        Err(_) => return Ok(None),
    };
    let len = normal.length();
    if !(len > 0.0) || !len.is_finite() {
        return Ok(None);
    }
    let sign = if eu1::faceuse_orientation(b, fu)? == OrientType::Same {
        1.0
    } else {
        -1.0
    };
    Ok(Some((normal * (sign / len), into)))
}

/// `sides_offsets`: the surface offsets of `e` (its blend edge use's surface first) at which
/// `chamfer_offsets` cuts a chamfer `d1` back along the face whose outward normal at the
/// edge's midpoint leans most along `up`, and `second` back along the other (`angle`: the
/// angle in radians the chamfer makes with the first face). With `c` the cosine between the
/// directions into the faces, setbacks `s0`, `s1` lie at offsets `(s1 - s0 c) / sin` and
/// `(s0 - s1 c) / sin`; an angle `phi` to the first face gives
/// `s1 = s0 sin phi / sin(theta + phi)`.
fn sides_offsets(
    k: &Kernel,
    b: &mut Brep,
    e: EdgeId,
    d1: f64,
    second: f64,
    angle: bool,
    up: Vec3,
) -> solidcore::Result<Option<(f64, f64)>> {
    let Some(u0) = edge::blend_edgeuse(k, b, e)? else {
        return Ok(None);
    };
    let Some(u1) = eu1::radial(b, u0)? else {
        return Ok(None);
    };
    let (Some((n0, b0)), Some((n1, b1))) = (face_frame(k, b, u0)?, face_frame(k, b, u1)?) else {
        return Ok(None);
    };
    let c = b0.dot(b1).clamp(-1.0, 1.0);
    let sn = (1.0 - c * c).max(0.0).sqrt();
    if !(sn > 1e-6) {
        return Ok(None);
    }
    let mut s0 = d1;
    let mut s1 = second;
    if angle {
        let d = (c.acos() + second).sin();
        if !(d > 1e-9) {
            return Ok(None);
        }
        s1 = d1 * second.sin() / d;
    }
    if n1.dot(up) > n0.dot(up) {
        core::mem::swap(&mut s0, &mut s1);
    }
    let r0 = (s1 - s0 * c) / sn;
    let r1 = (s0 - s1 * c) / sn;
    Ok((r0.is_finite() && r1.is_finite() && r0 > 0.0 && r1 > 0.0).then_some((r0, r1)))
}

/// `round_edges`: rounds `edges` of `b` in place, edge `i` at `sizes[i]`, as
/// `SmApiFilletEdgesPerEdge` does with one cross-section and solver for all:
/// [`BLEND_ROUND`] circular at constant radius, [`BLEND_CURVATURE`] the blend curve at
/// thirds (continuity 2) between the round's contact lines, [`BLEND_CHORD`] circular at
/// constant width. `law` makes the radius vary (one edge, [`BLEND_ROUND`] only). Each
/// setback corner (all its edges rounded, all of one convexity) is set back between making
/// the corners and filleting; with `patches`, the faces of each one's patch afterwards.
fn round_edges(
    k: &mut Kernel,
    b: &mut Brep,
    edges: &[EdgeId],
    sizes: &[f64],
    kind: i32,
    law: Option<FilletLaw>,
    setbacks: &[SetBack],
    patches: bool,
) -> core::result::Result<Vec<Vec<FaceId>>, RoundFail> {
    let (section, chord, continuity) = match kind {
        BLEND_ROUND => (FilletSurfaceGeneratorType::Circular, false, 1),
        BLEND_CURVATURE => (FilletSurfaceGeneratorType::BlendCurve, false, 2),
        BLEND_CHORD => (FilletSurfaceGeneratorType::Circular, true, 1),
        _ => return Err(Status::ERR_ARG.into()),
    };
    if law.is_some() && (edges.len() != 1 || kind != BLEND_ROUND) {
        return Err(Status::ERR_ARG.into());
    }
    let before: BTreeSet<FaceId> = b.face_list().iter().copied().collect();
    let mut found = Vec::new();
    {
        let mut session = FilletSession::new(k, b)?;
        // The curvature-continuous section is the quintic whose inner control points lie a
        // third and two thirds of the way from each contact to its tangents' corner.
        session.set_blend_curve_at_thirds(kind == BLEND_CURVATURE);
        let prepared = match law {
            Some(law) => {
                let size = sizes.first().copied().unwrap_or(0.0);
                session.prepare_variable_radius_edges(k, edges, size, law, 1, continuity)
            }
            None if chord => session.prepare_per_edge_constant_distance(
                k,
                edges,
                sizes,
                &[section],
                continuity,
                1.0,
            ),
            None => session.prepare_per_edge_constant_radius(
                k,
                edges,
                sizes,
                &[section],
                continuity,
                1.0,
            ),
        };
        prepared.map_err(|e| blend_failure(&session, e))?;
        session
            .create_corners(k, None)
            .map_err(|e| blend_failure(&session, e))?;
        let corners = session.fillet_corners()?;
        let mut set_back = Vec::new();
        for (j, sb) in setbacks.iter().enumerate() {
            let mut corner = None;
            for &c in &corners {
                if session.corner(c)?.filleted_vertex() == Some(sb.vertex) {
                    corner = Some(c);
                    break;
                }
            }
            let Some(corner) = corner else {
                return Err(Status::ERR_ARG.into());
            };
            // SMLib sets back an N x N corner whose edges all have one convexity; a mixed
            // one's sections come out kinked and the result no solid.
            if session.corner(corner)?.corner_type() != FilletCornerType::NByNConvex {
                return Err(RoundFail::Refused(j));
            }
            session.set_corner_setback(corner, sb.distance)?;
            set_back.push(corner);
        }
        session.set_mark_corner_patches(patches);
        session
            .do_filleting(k, None)
            .map_err(|e| blend_failure(&session, e))?;
        // (SplitFilletFacesAtRadiusCorners: see the module's deviations)
        if patches {
            let marked = session.corner_patch_faces()?;
            for corner in &set_back {
                let index = corners.iter().position(|c| c == corner);
                let patch: Vec<FaceId> = marked
                    .iter()
                    .filter(|(i, _)| Some(*i) == index)
                    .map(|&(_, f)| f)
                    .collect();
                found.push(patch);
            }
        }
    }
    stitch_and_orient(k, b, &StitchAndOrientOptions::default())?;
    // A round that made no face did nothing, whatever SMLib's status said.
    made_faces(b, &before)?;
    for patch in &found {
        if patch.is_empty() || patch.iter().any(|&f| !b.is_live_entity(EntityRef::Face(f))) {
            return Err(Status::ERR_UNSUPPORTED.into());
        }
    }
    Ok(found)
}

/// `setbacks_fit`: each setback shorter than every edge at its corner, and two on the ends of
/// one edge together shorter than it.
fn setbacks_fit(k: &Kernel, b: &Brep, setbacks: &[SetBack]) -> solidcore::Result<bool> {
    for sb in setbacks {
        for e in vertex::edges(b, sb.vertex, None)? {
            let Ok(length) = edge_length(k, b, e, 1e-6) else {
                return Ok(false);
            };
            let far = edge::other_vertex(b, e, Some(sb.vertex))?;
            let mut other = 0.0;
            for o in setbacks {
                if Some(o.vertex) == far && far != Some(sb.vertex) {
                    other = o.distance;
                }
            }
            if sb.distance + other >= length {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

/// `open_corners`: removes each patch's faces from the solid `b`, with the edges and vertices
/// only they used, and gives the edges round each hole as the faces beside it run them, each
/// starting where the one before ends. Each must lie between a face of its patch and one
/// face outside every patch, and close one loop, else `SMLIB_ERR_UNSUPPORTED`.
fn open_corners(
    k: &mut Kernel,
    b: &mut Brep,
    patches: &[Vec<FaceId>],
) -> core::result::Result<Vec<Vec<HoleEdge>>, Fail> {
    let unsupported = || Fail::Status(Status::ERR_UNSUPPORTED);
    let any: BTreeSet<FaceId> = patches.iter().flatten().copied().collect();
    let mut holes = Vec::new();
    for patch in patches {
        let mine: BTreeSet<FaceId> = patch.iter().copied().collect();
        let mut seen = BTreeSet::new();
        let mut rim = Vec::new();
        for &f in patch {
            for e in face::p4::edges(b, f, None, false, None)? {
                if !seen.insert(e) {
                    continue;
                }
                let beside = edge::faces(b, e, None)?;
                let others: Vec<FaceId> = beside
                    .iter()
                    .copied()
                    .filter(|g| !mine.contains(g))
                    .collect();
                let Some(&other) = others.first() else {
                    continue;
                };
                if others.len() != 1 || beside.len() != 2 || any.contains(&other) {
                    return Err(unsupported());
                }
                let fu = outer_faceuse(b, other)?.ok_or_else(unsupported)?;
                let mut found = None;
                for u in edge::edgeuses(b, e, None, None)? {
                    if eu1::faceuse(b, u)? == Some(fu) {
                        found = Some(u);
                    }
                }
                let u = found.ok_or_else(unsupported)?;
                rim.push(HoleEdge {
                    edge: e,
                    face: other,
                    forward: eu1::orientation(b, u)? == OrientType::Same,
                });
            }
        }
        let start = |h: &HoleEdge| {
            if h.forward {
                edge::start_vertex(b, h.edge)
            } else {
                edge::end_vertex(b, h.edge)
            }
        };
        let end = |h: &HoleEdge| {
            if h.forward {
                edge::end_vertex(b, h.edge)
            } else {
                edge::start_vertex(b, h.edge)
            }
        };
        let Some(&head) = rim.first() else {
            return Err(unsupported());
        };
        let mut lp = vec![head];
        let mut used = vec![false; rim.len()];
        if let Some(u) = used.first_mut() {
            *u = true;
        }
        let first_start = start(&head)?;
        loop {
            let last = *lp.last().ok_or_else(unsupported)?;
            let last_end = end(&last)?;
            if last_end == first_start {
                break;
            }
            let mut next = None;
            for (i, h) in rim.iter().enumerate() {
                if !used.get(i).copied().unwrap_or(true) && start(h)? == last_end {
                    if next.is_some() {
                        return Err(unsupported());
                    }
                    next = Some(i);
                }
            }
            let i = next.ok_or_else(unsupported)?;
            if let Some(u) = used.get_mut(i) {
                *u = true;
            }
            lp.push(*rim.get(i).ok_or_else(unsupported)?);
        }
        if lp.len() != rim.len() {
            return Err(unsupported());
        }
        holes.push(lp);
    }
    // A face is removed only with the B-rep open for editing.
    let editing = b.editing_enabled();
    b.set_editing_enabled(true);
    let mut state = ArenaNotifyState::new(b);
    let mut result = Ok(());
    for &f in &any {
        let mut env = state.env(None);
        if let Err(e) = delete_face(
            k,
            b,
            &mut env,
            &mut Marks::default(),
            f,
            DeleteFaceOptions {
                delete_edges_vertices: true,
                ..DeleteFaceOptions::default()
            },
        ) {
            result = Err(Fail::Solid(e));
            break;
        }
    }
    b.set_editing_enabled(editing);
    result.map(|()| holes)
}

/// `law_ok`: a clamped B-spline of `degree` (1 to 3) from 0 to 1, 2 to 100 control radii,
/// knots rising, each end repeated `degree + 1` times, radii finite in `[0, 1e6]`.
fn law_ok(degree: i32, knots: &[f64], radii: &[f64]) -> bool {
    let Ok(d) = usize::try_from(degree) else {
        return false;
    };
    let n = radii.len();
    if !(1..=3).contains(&d) || !(2..=100).contains(&n) || n <= d || knots.len() != n + d + 1 {
        return false;
    }
    if !finite(knots) || !finite(radii) {
        return false;
    }
    for (i, &kv) in knots.iter().enumerate() {
        let low_end = i <= d;
        let high_end = i >= n;
        let previous = i.checked_sub(1).and_then(|j| knots.get(j)).copied();
        if (low_end && kv != 0.0) || (high_end && kv != 1.0) || previous.is_some_and(|p| kv < p) {
            return false;
        }
    }
    radii.iter().all(|&r| (0.0..=1e6).contains(&r))
}

/// `law_pinches`: whether the law's radius is zero somewhere inside (0, 1).
fn law_pinches(degree: usize, knots: &[f64], radii: &[f64]) -> bool {
    let n = radii.len();
    let all_zero = |from: isize, to: isize| -> bool {
        let lo = usize::try_from(from.max(0)).unwrap_or(0);
        let mut j = lo;
        while (j as isize) <= to && j < n {
            if radii.get(j).copied().unwrap_or(0.0) != 0.0 {
                return false;
            }
            j += 1;
        }
        from <= to
    };
    let d = degree as isize;
    let n_knots = n + degree + 1;
    let mut i = degree;
    while i + 1 < n_knots && i < n {
        let (a, c) = (knots.get(i).copied(), knots.get(i + 1).copied());
        if matches!((a, c), (Some(a), Some(c)) if a < c) && all_zero(i as isize - d, i as isize) {
            return true;
        }
        i += 1;
    }
    let mut s = degree + 1;
    while s < n {
        let mut m = 1;
        while s + m < n && knots.get(s + m) == knots.get(s) {
            m += 1;
        }
        if all_zero(s as isize + m as isize - 1 - d, s as isize - 1) {
            return true;
        }
        s += m;
    }
    false
}

/// Whether the law may have a corner inside (0, 1): an interior knot repeated `degree`
/// times or more, where its slope can jump and the fork splits the round
/// (`SplitFilletFacesAtRadiusCorners`; see the module's deviations).
fn law_corner(degree: usize, knots: &[f64]) -> bool {
    let mut i = 0;
    while i < knots.len() {
        let kv = knots.get(i).copied().unwrap_or(0.0);
        let mut m = 1;
        while knots.get(i + m) == Some(&kv) {
            m += 1;
        }
        if kv > 0.0 && kv < 1.0 && m >= degree {
            return true;
        }
        i += m;
    }
    false
}

/// `law_along`: the radius law along `e` as a B-spline whose x is the radius and whose y is
/// the parameter itself (control points at the Greville abscissae), reversed when `start`
/// is nearer the end of the edge's interval than its start; its map is 0 to 1.
fn law_along(
    k: &Kernel,
    b: &Brep,
    e: EdgeId,
    start: [f64; 3],
    degree: usize,
    knots: &[f64],
    radii: &[f64],
) -> solidcore::Result<Option<FilletLaw>> {
    let Some(curve) = edge::curve(b, e)? else {
        return Ok(None);
    };
    let ivl = edge::interval(b, e)?;
    let (Ok(lo), Ok(hi)) = (
        curve.evaluate_point(k, ivl.min()),
        curve.evaluate_point(k, ivl.max()),
    ) else {
        return Ok(None);
    };
    let at = Point3::new(start[0], start[1], start[2]);
    let reverse = (hi - at).length() < (lo - at).length();
    let n_knots = knots.len();
    let kv: Vec<f64> = (0..n_knots)
        .map(|i| {
            if reverse {
                1.0 - knots.get(n_knots - 1 - i).copied().unwrap_or(0.0)
            } else {
                knots.get(i).copied().unwrap_or(0.0)
            }
        })
        .collect();
    let n = radii.len();
    let mut cps = Vec::new();
    for i in 0..n {
        let mut g = 0.0;
        for j in 1..=degree {
            g += kv.get(i + j).copied().unwrap_or(0.0);
        }
        let r = if reverse {
            radii.get(n - 1 - i)
        } else {
            radii.get(i)
        };
        cps.push(Point3::new(
            r.copied().unwrap_or(0.0),
            g / degree as f64,
            0.0,
        ));
    }
    let mut unique: Vec<f64> = Vec::new();
    let mut mults: Vec<usize> = Vec::new();
    for &x in &kv {
        match (unique.last(), mults.last_mut()) {
            (Some(&u), Some(m)) if u == x => *m += 1,
            _ => {
                unique.push(x);
                mults.push(1);
            }
        }
    }
    let Ok(c) = BSplineCurve::create_canonical(&CanonicalData {
        dim: 3,
        degree,
        control_points: &cps,
        form: BSplineCurveForm::default(),
        knot_multiplicities: &mults,
        knots: &unique,
        knot_type: KnotType::default(),
        weights: None,
        trim: None,
    }) else {
        return Ok(None);
    };
    Ok(Some(FilletLaw::bspline(
        c,
        &Extent1d::from_min_max(0.0, 1.0),
    )))
}

/// The arguments of one blend call (`blend_call`'s).
struct BlendCall<'a> {
    name: &'static str,
    kind: i32,
    points: &'a [[f64; 3]],
    sizes: &'a [f64],
    tolerance: f64,
    corners: &'a [[f64; 3]],
    setbacks: &'a [f64],
    open: bool,
}

impl Smlib {
    /// `blend_call`: the work of the blend entry points. Returns the body and, when `open`,
    /// the holes as `smlib_blend_setback_open` describes them.
    fn blend_call(&mut self, brep: &Brep, c: &BlendCall<'_>) -> Result<(Brep, Vec<i64>)> {
        let kind_ok = matches!(
            c.kind,
            BLEND_ROUND | BLEND_CURVATURE | BLEND_CHORD | BLEND_CHAMFER
        );
        if !kind_ok
            || !points_ok(c.points, false)
            || !tolerance_ok(c.tolerance)
            || c.sizes.len() != c.points.len()
            || !points_ok(c.corners, true)
            || c.setbacks.len() != c.corners.len()
        {
            return Err(Status::ERR_ARG);
        }
        if !c.sizes.iter().all(|&s| size_ok(s)) || !c.setbacks.iter().all(|&s| size_ok(s)) {
            return Err(Status::ERR_ARG);
        }
        // A chamfer's corners are bevelled, which builds no setback.
        if !c.corners.is_empty() && c.kind == BLEND_CHAMFER {
            return Err(Status::ERR_UNSUPPORTED);
        }
        let line = || {
            let mut args = format!(
                "{}(brep, kind={}, points={}, sizes={}, tolerance={}",
                c.name,
                c.kind,
                pts(c.points),
                vec(c.sizes),
                num(c.tolerance)
            );
            if !c.corners.is_empty() {
                args += &format!(", corners={}, setbacks={}", pts(c.corners), vec(c.setbacks));
            }
            args + ")"
        };
        self.run(line, &[("brep", brep)], |api| {
            let mut w = api.copy(brep)?;
            let edges = pick_edges(api.kernel(), &w, c.points, c.tolerance);
            let edges = edges.map_err(|e| {
                if let Some(i) = e.index() {
                    api.set_failed_point(i);
                }
                Fail::from(e)
            })?;
            let (picked, first) = (edges.items, edges.first);
            let each: Vec<f64> = first
                .iter()
                .map(|&i| c.sizes.get(i).copied().unwrap_or(0.0))
                .collect();
            // Each corner point picks a vertex; a vertex picked twice takes the first
            // point's setback.
            let corners = pick_vertices(&w, c.corners, c.tolerance).map_err(|e| {
                if let Some(i) = e.index() {
                    api.set_failed_corner(i);
                }
                Fail::from(e)
            })?;
            let (at, first_corner) = (corners.items, corners.first);
            let back: Vec<SetBack> = at
                .iter()
                .zip(&first_corner)
                .map(|(&vertex, &j)| SetBack {
                    vertex,
                    distance: c.setbacks.get(j).copied().unwrap_or(0.0),
                })
                .collect();
            if !setbacks_fit(api.kernel(), &w, &back)? {
                return Err(Fail::Status(Status::ERR_NO_ROOM));
            }
            // The corners' places, read now: filleting deletes the vertices.
            let mut at_points = Vec::new();
            for &v in &at {
                at_points.push(vertex::point(&w, v)?);
            }
            let solid_in = brep_is_manifold_solid(&w)?;
            let k = api.kernel_mut();
            let patches = if c.kind == BLEND_CHAMFER {
                chamfer_edges(k, &mut w, &picked, &each)?;
                Vec::new()
            } else {
                match round_edges(k, &mut w, &picked, &each, c.kind, None, &back, c.open) {
                    Ok(p) => p,
                    Err(RoundFail::Fail(f)) => return Err(f),
                    // A corner SMLib does not set back: the first point that picked it.
                    Err(RoundFail::Refused(j)) => {
                        if let Some(&i) = first_corner.get(j) {
                            api.set_failed_corner(i);
                        }
                        return Err(Fail::Status(Status::ERR_UNSUPPORTED));
                    }
                }
            };
            if !c.open {
                return finish(w, solid_in, false).map(|b| (b, Vec::new()));
            }
            // Each corner opened: its patch removed, the edges round its hole in the
            // result's orders, once per corner point.
            if !brep_is_manifold_solid(&w)? {
                return Err(Fail::Status(Status::ERR_NOT_SOLID));
            }
            let rims = open_corners(api.kernel_mut(), &mut w, &patches)?;
            let all_edges = p1::edges(&w, None)?;
            let all_faces = p1::faces(&w, None)?;
            let mut values = Vec::new();
            for p in c.corners {
                let p = Point3::new(p[0], p[1], p[2]);
                let mut which = None;
                let mut best = f64::INFINITY;
                for (k, q) in at_points.iter().enumerate() {
                    let d = (*q - p).length();
                    if d < best {
                        best = d;
                        which = Some(k);
                    }
                }
                let rim = which
                    .and_then(|w| rims.get(w))
                    .ok_or(Fail::Status(Status::ERR_ARG))?;
                values.push(rim.len() as i64);
                for h in rim {
                    let e = all_edges.iter().position(|&x| x == h.edge);
                    let f = all_faces.iter().position(|&x| x == h.face);
                    let (Some(e), Some(f)) = (e, f) else {
                        return Err(Fail::Status(Status::ERR_UNSUPPORTED));
                    };
                    values.extend([e as i64, f as i64, i64::from(h.forward)]);
                }
            }
            finish(w, false, false).map(|b| (b, values))
        })
    }

    /// `smlib_fillet`: a new body, `brep` with the edges `points` select rounded at
    /// `radius` (a circular cross-section at constant radius). A point that selects no edge
    /// or several is `SMLIB_ERR_PICK_NONE` or `SMLIB_ERR_PICK_AMBIGUOUS`, its index in
    /// [`Smlib::last_failed_point`].
    pub fn fillet(
        &mut self,
        brep: &Brep,
        points: &[[f64; 3]],
        tolerance: f64,
        radius: f64,
    ) -> Result<Brep> {
        let sizes = vec![radius; points.len()];
        self.blend_call(
            brep,
            &BlendCall {
                name: "smlib_fillet",
                kind: BLEND_ROUND,
                points,
                sizes: &sizes,
                tolerance,
                corners: &[],
                setbacks: &[],
                open: false,
            },
        )
        .map(|(b, _)| b)
    }

    /// `smlib_chamfer`: a new body, `brep` with the selected edges cut by an equal-distance
    /// chamfer `distance` back from the edge along each face (measured at each edge's
    /// midpoint), with bevelled corners.
    pub fn chamfer(
        &mut self,
        brep: &Brep,
        points: &[[f64; 3]],
        tolerance: f64,
        distance: f64,
    ) -> Result<Brep> {
        let sizes = vec![distance; points.len()];
        self.blend_call(
            brep,
            &BlendCall {
                name: "smlib_chamfer",
                kind: BLEND_CHAMFER,
                points,
                sizes: &sizes,
                tolerance,
                corners: &[],
                setbacks: &[],
                open: false,
            },
        )
        .map(|(b, _)| b)
    }

    /// `smlib_blend`: a new body, `brep` with the selected edges blended as `kind`
    /// (`BLEND_*`), each at its point's size (`sizes`, one per point, each > 0); an edge
    /// reached from several points takes the first point's size.
    pub fn blend(
        &mut self,
        brep: &Brep,
        kind: i32,
        points: &[[f64; 3]],
        sizes: &[f64],
        tolerance: f64,
    ) -> Result<Brep> {
        self.blend_call(
            brep,
            &BlendCall {
                name: "smlib_blend",
                kind,
                points,
                sizes,
                tolerance,
                corners: &[],
                setbacks: &[],
                open: false,
            },
        )
        .map(|(b, _)| b)
    }

    /// `smlib_blend_setback`: as [`Smlib::blend`], with corner setbacks: each point of
    /// `corners` (at least one) picks a vertex where every edge is blended, of one
    /// convexity, and moves where that corner's patch meets each blend `setbacks[j]` along
    /// the edge from the vertex. A corner SMLib does not set back is
    /// `SMLIB_ERR_UNSUPPORTED` ([`Smlib::last_failed_corner`] names its first point); a
    /// setback as long as an edge there is `SMLIB_ERR_NO_ROOM`.
    pub fn blend_setback(
        &mut self,
        brep: &Brep,
        kind: i32,
        points: &[[f64; 3]],
        sizes: &[f64],
        corners: &[[f64; 3]],
        setbacks: &[f64],
        tolerance: f64,
    ) -> Result<Brep> {
        if corners.is_empty() {
            return Err(Status::ERR_ARG);
        }
        self.blend_call(
            brep,
            &BlendCall {
                name: "smlib_blend_setback",
                kind,
                points,
                sizes,
                tolerance,
                corners,
                setbacks,
                open: false,
            },
        )
        .map(|(b, _)| b)
    }

    /// `smlib_blend_setback_open`: as [`Smlib::blend_setback`], with each set-back corner's
    /// patch removed: a sheet whose only free edges are those round the holes, and the
    /// holes: for each point of `corners` in order, `n` then `n` triples (edge, face,
    /// forward), the edges round that corner's hole as the faces beside it run them
    /// (counter-clockwise about their outward normal), each starting where the one before
    /// ends; `edge` and `face` indices in the result's edge and face orders, `forward` 1
    /// where the face runs the edge from its start vertex to its end. A patch that is not a
    /// disc of faces of its own is `SMLIB_ERR_UNSUPPORTED`.
    pub fn blend_setback_open(
        &mut self,
        brep: &Brep,
        kind: i32,
        points: &[[f64; 3]],
        sizes: &[f64],
        corners: &[[f64; 3]],
        setbacks: &[f64],
        tolerance: f64,
    ) -> Result<(Brep, Vec<i64>)> {
        if corners.is_empty() {
            return Err(Status::ERR_ARG);
        }
        self.blend_call(
            brep,
            &BlendCall {
                name: "smlib_blend_setback_open",
                kind,
                points,
                sizes,
                tolerance,
                corners,
                setbacks,
                open: true,
            },
        )
    }

    /// `smlib_fillet_variable`: a new body, `brep` with the one edge `point` selects
    /// rounded by a circular cross-section whose radius follows a clamped B-spline of
    /// `degree` (1 to 3) over the position along the edge, 0 at its end nearest `start`:
    /// `knots` (`radii.len() + degree + 1`, rising from 0 to 1, each end repeated
    /// `degree + 1` times) and 2 to 100 control `radii`, each >= 0. A law whose radius is
    /// zero inside the edge is `SMLIB_ERR_UNSUPPORTED`, and so (for now) is one with a
    /// corner inside the edge (see the module's deviations).
    pub fn fillet_variable(
        &mut self,
        brep: &Brep,
        point: [f64; 3],
        tolerance: f64,
        start: [f64; 3],
        degree: i32,
        knots: &[f64],
        radii: &[f64],
    ) -> Result<Brep> {
        if !finite(&point)
            || !finite(&start)
            || !tolerance_ok(tolerance)
            || !law_ok(degree, knots, radii)
        {
            return Err(Status::ERR_ARG);
        }
        let d = usize::try_from(degree).map_err(|_| Status::ERR_ARG)?;
        if law_pinches(d, knots, radii) || law_corner(d, knots) {
            return Err(Status::ERR_UNSUPPORTED);
        }
        let line = || {
            format!(
                "smlib_fillet_variable(brep, point={}, tolerance={}, start={}, degree={degree}, knots={}, radii={})",
                vec(&point),
                num(tolerance),
                vec(&start),
                vec(knots),
                vec(radii)
            )
        };
        self.run(line, &[("brep", brep)], |api| {
            let mut w = api.copy(brep)?;
            // (no index is written: the call has one point)
            let picked = pick_edges(api.kernel(), &w, &[point], tolerance)?;
            let Some(&e) = picked.items.first() else {
                return Err(Fail::Status(Status::ERR_PICK_NONE));
            };
            let Some(law) = law_along(api.kernel(), &w, e, start, d, knots, radii)? else {
                return Err(Fail::Status(Status::SM_ERR_INVALID_INPUT));
            };
            let solid_in = brep_is_manifold_solid(&w)?;
            // The size is the start radius.
            let first = radii.first().copied().unwrap_or(0.0);
            let last = radii.last().copied().unwrap_or(0.0);
            let size = [first.max(last)];
            match round_edges(
                api.kernel_mut(),
                &mut w,
                &[e],
                &size,
                BLEND_ROUND,
                Some(law),
                &[],
                false,
            ) {
                Ok(_) => {}
                Err(RoundFail::Fail(f)) => return Err(f),
                Err(RoundFail::Refused(_)) => return Err(Fail::Status(Status::ERR_UNSUPPORTED)),
            }
            finish(w, solid_in, false)
        })
    }

    /// `smlib_chamfer_sides`: a new body, `brep` with the selected edges cut by an unequal
    /// chamfer with bevelled corners: `distance` back from the edge along the face whose
    /// outward normal at the edge's midpoint leans most along `first` (non-zero), and along
    /// the other `second` back, or (`second_is_angle` 1) at the angle `second` (radians, in
    /// (0, pi)) the chamfer makes with the first face. Setbacks are measured at each edge's
    /// midpoint.
    pub fn chamfer_sides(
        &mut self,
        brep: &Brep,
        points: &[[f64; 3]],
        tolerance: f64,
        distance: f64,
        second: f64,
        second_is_angle: i32,
        first: [f64; 3],
    ) -> Result<Brep> {
        let angle = second_is_angle == 1;
        let second_ok = if angle {
            second.is_finite() && second > 1e-6 && second < PI - 1e-6
        } else {
            size_ok(second)
        };
        if !points_ok(points, false)
            || !tolerance_ok(tolerance)
            || !size_ok(distance)
            || (second_is_angle != 0 && second_is_angle != 1)
            || !second_ok
            || !finite(&first)
        {
            return Err(Status::ERR_ARG);
        }
        let up = Vec3::new(first[0], first[1], first[2]);
        if !(up.length() > 0.0) {
            return Err(Status::ERR_ARG);
        }
        let line = || {
            format!(
                "smlib_chamfer_sides(brep, points={}, tolerance={}, distance={}, second={}, angle={second_is_angle}, first={})",
                pts(points),
                num(tolerance),
                num(distance),
                num(second),
                vec(&first)
            )
        };
        self.run(line, &[("brep", brep)], |api| {
            let mut w = api.copy(brep)?;
            let picked = pick_edges(api.kernel(), &w, points, tolerance);
            let picked = picked
                .map_err(|e| {
                    if let Some(i) = e.index() {
                        api.set_failed_point(i);
                    }
                    Fail::from(e)
                })?
                .items;
            let solid_in = brep_is_manifold_solid(&w)?;
            let mut r0 = Vec::new();
            let mut r1 = Vec::new();
            for &e in &picked {
                let Some((a, b)) =
                    sides_offsets(api.kernel(), &mut w, e, distance, second, angle, up)?
                else {
                    return Err(Fail::Status(Status::SM_ERR_INVALID_INPUT));
                };
                r0.push(a);
                r1.push(b);
            }
            chamfer_offsets(api.kernel_mut(), &mut w, &picked, &r0, &r1)?;
            finish(w, solid_in, false)
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn law_checks_follow_the_c_api() {
        assert!(law_ok(1, &[0.0, 0.0, 1.0, 1.0], &[3.0, 1.0]));
        assert!(!law_ok(1, &[0.0, 0.5, 1.0, 1.0], &[3.0, 1.0]));
        assert!(!law_ok(4, &[0.0, 0.0, 1.0, 1.0], &[3.0, 1.0]));
        assert!(!law_ok(1, &[0.0, 0.0, 1.0], &[3.0, 1.0]));
        // zero only at an end is no pinch; zero inside is
        assert!(!law_pinches(1, &[0.0, 0.0, 1.0, 1.0], &[0.0, 1.0]));
        assert!(law_pinches(1, &[0.0, 0.0, 0.5, 1.0, 1.0], &[1.0, 0.0, 1.0]));
        // a degree-1 law with an interior knot can turn there; a quadratic with simple
        // knots cannot
        assert!(law_corner(1, &[0.0, 0.0, 0.5, 1.0, 1.0]));
        assert!(!law_corner(2, &[0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0]));
        assert!(!law_corner(1, &[0.0, 0.0, 1.0, 1.0]));
    }
}
