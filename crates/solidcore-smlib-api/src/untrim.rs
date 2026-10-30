//! Untrim Face: a face restored towards its surface. Picked inner loops (holes) are
//! taken away with `SmFace::DeleteInnerLoop` on a copy where they are lamina; the outer
//! boundary goes back to a parameter box of the surface (`SmBrep::CreateFaceFromSurface`,
//! a face bounded by the box's own iso-parameter edges), and the holes kept are rebuilt
//! in that face from their own edges' curves and intervals (`SmBrep::MakeVertexLoop` and
//! `MakeEdgeInFace`, as `SmBrep::CreateFacesWithCurves` inserts loops: closing a loop
//! splits off the face inside it, which is deleted). The new face and the rebuilt edges
//! take the attributes of the face and edges they stand for, so the history tags carry
//! through.
//!
//! Source: tangency-labs/smlib `capi/src/smlib_untrim.cpp` (Apache-2.0, Copyright
//! SolidCraft contributors) at fork commit f1a39ad (replay/1.0.1).
//!
//! | Rust | C API (SMLib call) |
//! |------|-------|
//! | [`Smlib::untrim_face`] | `smlib_untrim_face` |
//! | `face_at`, `loop_distance`, `copy_attributes`, `analytic`, `span_of`, `Piece`, `pieces_of`, `is_disc`, `rebuild_hole` (private) | the same |
//!
//! Deviations: the hole points are a slice of points (`None` for the C API's negative
//! count: every hole), `outer` and `in_place` are `bool`s. `SmAttribute::MakeCopy` is
//! [`Attribute::make_copy`]. The face's own copy (`CopyFaces` through `SmMerge`) notifies
//! through one arena-only store for both bodies, and runs on a deep copy of the body (the
//! C API's reads the body itself). The holes' edge curves, intervals and attributes are
//! copied out before the new face is built (the C API keeps pointers to the old edges).

use solidcore::geom::attribute::Attribute;
use solidcore::geom::extent2d::Extent2d;
use solidcore::geom::object::ObjectRef;
use solidcore::geom::tol::USE_DEFAULT;
use solidcore::geom::vector3d::Point3;
use solidcore::nurbs::b_spline_surface_n_lib::p2::OrientType;
use solidcore::surface::surface::SurfaceKind;
use solidcore::topology::brep::p1 as bp1;
use solidcore::topology::brep::p2::ArenaNotifyState;
use solidcore::topology::brep::p3::{delete_face, new_brep, DeleteFaceOptions};
use solidcore::topology::brep::p5::{make_edge_in_face, EdgeInFace};
use solidcore::topology::brep::p7::copy_faces;
use solidcore::topology::brep::p8::{create_face_from_surface, FaceBuildEnv};
use solidcore::topology::brep::p9::make_vertex_loop;
use solidcore::topology::edge::{self, NotifyEnv};
use solidcore::topology::edgeuse::{p1 as eu1, p3 as eu3};
use solidcore::topology::face::p1::{delete_inner_loop, flip_face_orientation};
use solidcore::topology::face::p2 as face2;
use solidcore::topology::loopuse;
use solidcore::topology::topology_traverser::Marks;
use solidcore::topology::{Brep, EdgeId, FaceId, LoopId};
use solidcore::{EntityRef, Kernel};

use crate::pick::face_distance;
use crate::session::{finite, num, take, vec, Fail, Smlib};
use crate::status::{Result, Status};

/// The most hole points.
const MAX_HOLE_POINTS: usize = 10_000;

fn empty() -> Fail {
    Fail::Status(Status::ERR_EMPTY)
}

fn unsupported() -> Fail {
    Fail::Status(Status::ERR_UNSUPPORTED)
}

/// The face of `b` nearest `p`, within `tol`.
fn face_at(k: &mut Kernel, b: &mut Brep, p: Point3, tol: f64) -> solidcore::Result<Option<FaceId>> {
    let mut best = None;
    let mut d0 = f64::INFINITY;
    for f in bp1::faces(b, None)? {
        let d = face_distance(k, b, f, p).map_or(f64::INFINITY, |(d, _)| d);
        if d.is_finite() && d < d0 {
            d0 = d;
            best = Some(f);
        }
    }
    Ok(best.filter(|_| d0 <= tol))
}

/// The edges of loop `l` (`SmLoop::GetEdges`: its first loop use's).
fn loop_edges(b: &Brep, l: LoopId) -> solidcore::Result<Vec<EdgeId>> {
    let (lu, _) = loopuse::loop_loopuses(b, l)?;
    Ok(loopuse::edges(b, lu, None)?.into_iter().flatten().collect())
}

fn loop_distance(k: &Kernel, b: &Brep, l: LoopId, p: Point3) -> f64 {
    let Ok(edges) = loop_edges(b, l) else {
        return f64::INFINITY;
    };
    edges
        .into_iter()
        .map(|e| {
            solidcore::query::api_queries::edge_closest_point(k, b, e, p)
                .ok()
                .map(|c| c.closest.distance)
                .filter(|d| d.is_finite())
                .unwrap_or(f64::INFINITY)
        })
        .fold(f64::INFINITY, f64::min)
}

/// Copies of the attributes of `from` in `b` (`MakeCopy`), for another object.
fn attribute_copies(b: &Brep, from: EntityRef) -> Vec<Attribute> {
    b.entity_attributes(from)
        .into_iter()
        .map(|(_, a)| a.make_copy())
        .collect()
}

fn add_attributes(b: &mut Brep, to: EntityRef, copies: Vec<Attribute>) -> solidcore::Result<()> {
    for a in copies {
        b.add_attribute(to, a)?;
    }
    Ok(())
}

/// Whether the surface is one of SMLib's analytic kinds (each is held as a B-spline
/// subclass): it has no parameter box of its own, so a face goes back to the box of its
/// own span.
fn analytic(kind: &SurfaceKind) -> bool {
    !matches!(kind, SurfaceKind::BSpline(_))
}

/// The box of the face's outer loop in the surface's parameters: the union of the tight
/// boxes of its edge uses' UV trim curves (made where the face has none attached yet).
fn span_of(k: &Kernel, b: &mut Brep, outer: LoopId) -> solidcore::Result<Option<Extent2d>> {
    let (lu, _) = loopuse::loop_loopuses(b, outer)?;
    let eus = eu1::loopuse_edgeuses(b, lu)?;
    let mut any = false;
    let (mut lo, mut hi) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
    for eu in eus {
        let fetched = eu3::get_or_create_uv_trim_curve(k, b, eu, false, false, true)?;
        if fetched.status.is_err() {
            return Ok(None);
        }
        let Some(uv) = b.edgeuse(eu)?.uv_trim_curve.as_ref() else {
            return Ok(None);
        };
        let Ok(e) = uv.calculate_tight_bounding_box(k, &uv.natural_interval()) else {
            return Ok(None);
        };
        let (m, x) = (e.min(), e.max());
        lo = [lo[0].min(m.x), lo[1].min(m.y)];
        hi = [hi[0].max(x.x), hi[1].max(x.y)];
        any = true;
    }
    let ok = any && hi[0] > lo[0] && hi[1] > lo[1] && (lo[0] + lo[1] + hi[0] + hi[1]).is_finite();
    if !ok {
        return Ok(None);
    }
    Ok(Some(Extent2d::from_bounds(lo[0], lo[1], hi[0], hi[1])))
}

/// One edge of a hole to rebuild: (a copy of) its edge's curve and interval, how the loop
/// runs along it, its ends as the loop runs, and copies of its attributes.
struct Piece {
    curve: solidcore::curve::curve::Curve,
    interval: solidcore::geom::extent1d::Extent1d,
    orient: OrientType,
    from: Point3,
    to: Point3,
    attributes: Vec<Attribute>,
}

fn pieces_of(k: &Kernel, b: &Brep, l: LoopId) -> solidcore::Result<Option<Vec<Piece>>> {
    let (lu, _) = loopuse::loop_loopuses(b, l)?;
    let mut out = Vec::new();
    for eu in eu1::loopuse_edgeuses(b, lu)? {
        let Some(e) = eu1::edge(b, eu)? else {
            return Ok(None);
        };
        let data = b.edge(e)?;
        let Some(curve) = data.curve.as_ref() else {
            return Ok(None);
        };
        let (curve, interval) = (curve.copy(k)?, data.interval);
        let vs = edge::vertices(b, e, None)?;
        let (Some(Some(v0)), v1) = (vs.first().copied(), vs.get(1).copied()) else {
            return Ok(None);
        };
        let v1 = v1.flatten().unwrap_or(v0);
        let same = eu1::orientation(b, eu)? == OrientType::Same;
        let (p0, p1) = (b.vertex(v0)?.point, b.vertex(v1)?.point);
        out.push(Piece {
            curve,
            interval,
            attributes: attribute_copies(b, EntityRef::Edge(e)),
            orient: if same {
                OrientType::Same
            } else {
                OrientType::Opposite
            },
            from: if same { p0 } else { p1 },
            to: if same { p1 } else { p0 },
        });
    }
    Ok((!out.is_empty()).then_some(out))
}

/// Whether `f` is bounded by one loop of exactly the edges `made`.
fn is_disc(b: &Brep, f: FaceId, made: &[EdgeId]) -> solidcore::Result<bool> {
    let loops = face2::loops(b, f, false, None)?;
    let &[l] = loops.as_slice() else {
        return Ok(false);
    };
    let edges = loop_edges(b, l)?;
    Ok(edges.len() == made.len() && edges.iter().all(|e| made.contains(e)))
}

/// The hole `pieces` rebuilt in `face` of `r`: vertices and edges made in the face, the
/// face the closed loop splits off inside it deleted. Returns the face that keeps the
/// rest. Making edges in a face is an edit, which SMLib allows only in edit mode.
fn rebuild_hole(
    k: &mut Kernel,
    r: &mut Brep,
    env: &mut NotifyEnv<'_>,
    face: FaceId,
    pieces: &[Piece],
) -> core::result::Result<FaceId, Fail> {
    let editing = r.editing_enabled();
    r.set_editing_enabled(true);
    let result = rebuild_hole_in(k, r, env, face, pieces);
    r.set_editing_enabled(editing);
    result
}

fn rebuild_hole_in(
    k: &mut Kernel,
    r: &mut Brep,
    env: &mut NotifyEnv<'_>,
    face: FaceId,
    pieces: &[Piece],
) -> core::result::Result<FaceId, Fail> {
    let tol = r.tolerance();
    let start = pieces.first().ok_or_else(empty)?;
    let (_, first) = make_vertex_loop(r, face, start.from, OrientType::Opposite, env)?;
    let mut at = first;
    let mut made = Vec::new();
    let mut split = None;
    for (i, p) in pieces.iter().enumerate() {
        let last = i + 1 == pieces.len();
        let next = if last {
            first
        } else {
            make_vertex_loop(r, face, p.to, OrientType::Opposite, env)?.1
        };
        let mut out = EdgeInFace::default();
        let editing = r.editing_enabled();
        make_edge_in_face(
            k,
            r,
            editing,
            Some(face),
            Some(at),
            Some(next),
            p.curve.clone(),
            None,
            p.interval,
            p.orient,
            tol,
            &mut out,
            None,
        )?;
        let new_edge = out.edge.ok_or_else(empty)?;
        add_attributes(r, EntityRef::Edge(new_edge), p.attributes.clone())?;
        made.push(new_edge);
        if let Some(nf) = out.face {
            if !last || split.is_some() {
                return Err(Fail::Status(Status::ERR_EXCEPTION));
            }
            split = Some(nf);
        }
        at = next;
    }
    let split = split.ok_or(Fail::Status(Status::ERR_EXCEPTION))?;
    // One of the two faces is the disc inside the hole: it goes, its edges stay on the
    // other.
    let (disc, rest) = if is_disc(r, split, &made)? {
        (split, face)
    } else if is_disc(r, face, &made)? {
        (face, split)
    } else {
        return Err(Fail::Status(Status::ERR_EXCEPTION));
    };
    let mut marks = Marks::default();
    let options = DeleteFaceOptions {
        delete_edges_vertices: false,
        do_region_nesting: true,
        region_to_keep: None,
    };
    delete_face(k, r, env, &mut marks, disc, options)?;
    Ok(rest)
}

/// A face on `surface` over `domain` in `b` (`CreateFaceFromSurface`).
fn face_from_surface(
    k: &mut Kernel,
    b: &mut Brep,
    env: &mut NotifyEnv<'_>,
    surface: solidcore::surface::surface::Surface,
    domain: &Extent2d,
) -> solidcore::Result<FaceId> {
    let key_of: &dyn Fn(EntityRef) -> ObjectRef = env.key_of;
    let this_face = |f: FaceId| key_of(EntityRef::Face(f));
    let mut marks = Marks::default();
    let mut fenv = FaceBuildEnv {
        env,
        marks: &mut marks,
        this_face: &this_face,
    };
    create_face_from_surface(k, b, Some(surface), domain, &mut fenv)
}

impl Smlib {
    /// `smlib_untrim_face`: a new body, the face of `body` nearest `point` (within `tol`)
    /// with holes taken away (`holes`: the hole nearest each point, within `tol`; `None`:
    /// all of them) and, with `outer`, its outer boundary put back to its surface's
    /// parameter box (an analytic surface's: the span of the face), the holes not taken
    /// away kept. With `in_place` the result is the whole body with the face replaced;
    /// otherwise the one face alone.
    ///
    /// # Errors
    ///
    /// `SMLIB_ERR_ARG` for a non-finite point, a tolerance not above 0, more than 10000
    /// hole points, or nothing to do (no holes and not `outer`); `SMLIB_ERR_PICK_NONE` when
    /// no face or hole lies within `tol`; `SMLIB_ERR_UNSUPPORTED` for a face without loops
    /// or surface, a hole another face uses, or a span that cannot be measured;
    /// `SMLIB_ERR_EMPTY` / `SMLIB_ERR_EXCEPTION` where SMLib makes nothing or an
    /// unexpected split; the kernel's status otherwise.
    pub fn untrim_face(
        &mut self,
        body: &Brep,
        point: [f64; 3],
        tol: f64,
        holes: Option<&[[f64; 3]]>,
        outer: bool,
        in_place: bool,
    ) -> Result<Brep> {
        let n_holes = holes.map_or(0, <[[f64; 3]]>::len);
        let ok = finite(&point)
            && tol.is_finite()
            && tol > 0.0
            && n_holes <= MAX_HOLE_POINTS
            && holes.is_none_or(|h| h.iter().all(|p| finite(p)))
            && !(holes.is_some() && n_holes == 0 && !outer);
        if !ok {
            return Err(Status::ERR_ARG);
        }
        self.run(
            || {
                let count = holes.map_or_else(|| "-1".to_string(), |h| h.len().to_string());
                format!(
                    "smlib_untrim_face(body, {}, {}, holes={count}, outer={}, in_place={})",
                    vec(&point),
                    num(tol),
                    i32::from(outer),
                    i32::from(in_place)
                )
            },
            &[("body", body)],
            |api| untrim(api, body, point, tol, holes, outer, in_place),
        )
    }
}

#[allow(clippy::too_many_lines)] // one capi routine
fn untrim(
    api: &mut Smlib,
    body: &Brep,
    point: [f64; 3],
    tol: f64,
    holes: Option<&[[f64; 3]]>,
    outer: bool,
    in_place: bool,
) -> core::result::Result<Brep, Fail> {
    let at = Point3::new(point[0], point[1], point[2]);
    // The face to untrim, in a body of its own (its hole edges lamina there) or in a copy
    // of the whole body.
    let mut whole = api.copy(body)?;
    let k = api.kernel_mut();
    let (mut w, f) = if in_place {
        let f = face_at(k, &mut whole, at, tol)?;
        (whole, f)
    } else {
        let src = face_at(k, &mut whole, at, tol)?.ok_or(Fail::Status(Status::ERR_PICK_NONE))?;
        let mut w = new_brep(k, USE_DEFAULT, USE_DEFAULT)?;
        w.set_tolerance(whole.tolerance());
        // Through SmMerge, which copies the topology with its attributes.
        let mut state = ArenaNotifyState::new(&whole);
        state.seed(&w);
        let mut env = state.env(None);
        let key_of: &dyn Fn(EntityRef) -> ObjectRef = env.key_of;
        let this_face = |x: FaceId| key_of(EntityRef::Face(x));
        let (copied, _) = copy_faces(k, &mut whole, &[src], &mut w, true, &this_face, &mut env)?;
        let &[f] = copied.as_slice() else {
            return Err(empty());
        };
        (w, Some(f))
    };
    let f = f.ok_or(Fail::Status(Status::ERR_PICK_NONE))?;
    let loops = face2::loops(&w, f, false, None)?;
    let Some(&outer_loop) = loops.first() else {
        return Err(unsupported());
    };
    let inner = loops.get(1..).unwrap_or_default();
    // The holes taken away: all of them, or the one nearest each point.
    let mut removed: Vec<LoopId> = Vec::new();
    match holes {
        None => removed.extend_from_slice(inner),
        Some(points) => {
            for p in points {
                let p = Point3::new(p[0], p[1], p[2]);
                let mut best = None;
                let mut d0 = f64::INFINITY;
                for &l in inner {
                    let d = loop_distance(k, &w, l, p);
                    if d < d0 {
                        d0 = d;
                        best = Some(l);
                    }
                }
                let best = best
                    .filter(|_| d0 <= tol)
                    .ok_or(Fail::Status(Status::ERR_PICK_NONE))?;
                if !removed.contains(&best) {
                    removed.push(best);
                }
            }
        }
    }
    if !outer {
        let mut state = ArenaNotifyState::new(&w);
        let mut env = state.env(None);
        for &l in &removed {
            // A hole whose edges another face of the body uses is no hole of a sheet.
            if delete_inner_loop(k, &mut w, &mut env, f, l).is_err() {
                return Err(unsupported());
            }
        }
        return take(w);
    }
    // The outer boundary: the surface's parameter box (a free-form surface's own, an
    // analytic one's the span of the face).
    let surf = face2::surface(&w, f)?.ok_or_else(unsupported)?.clone();
    let mut bx = surf.natural_uv_domain_with_nurb(k)?;
    if analytic(surf.kind_with_nurb(k)?) {
        bx = span_of(k, &mut w, outer_loop)?.ok_or_else(unsupported)?;
    }
    let mut kept = Vec::new();
    for &l in inner {
        if removed.contains(&l) {
            continue;
        }
        kept.push(pieces_of(k, &w, l)?.ok_or_else(unsupported)?);
    }
    let copy = surf.copy(k)?;
    let face_copies = attribute_copies(&w, EntityRef::Face(f));
    let was = match face2::upward_faceuse(&w, f)? {
        Some(fu) => Some(eu1::faceuse_orientation(&w, fu)?),
        None => None,
    };
    if in_place {
        let mut state = ArenaNotifyState::new(&w);
        let mut env = state.env(None);
        new_face(k, &mut w, &mut env, copy, &bx, &kept, face_copies, was)?;
        // The old face goes, with the edges and vertices only it used: an edit of the
        // brep, which SMLib allows only in edit mode (`SmTemporaryChangeValue` on
        // `m_bEditingEnabled`, restored however the delete ends).
        let mut marks = Marks::default();
        let options = DeleteFaceOptions {
            delete_edges_vertices: true,
            do_region_nesting: true,
            region_to_keep: None,
        };
        let editing = w.editing_enabled();
        w.set_editing_enabled(true);
        let deleted = delete_face(k, &mut w, &mut env, &mut marks, f, options);
        w.set_editing_enabled(editing);
        deleted?;
        return take(w);
    }
    let mut r = new_brep(k, USE_DEFAULT, USE_DEFAULT)?;
    r.set_tolerance(w.tolerance());
    let mut state = ArenaNotifyState::new(&r);
    let mut env = state.env(None);
    new_face(k, &mut r, &mut env, copy, &bx, &kept, face_copies, was)?;
    take(r)
}

/// The face on `surface` over `bx` in `target`, the `kept` holes rebuilt in it, the old
/// face's attributes given to it, and turned to face as the old one did (`was`).
#[allow(clippy::too_many_arguments)] // one step of the capi routine
fn new_face(
    k: &mut Kernel,
    target: &mut Brep,
    env: &mut NotifyEnv<'_>,
    surface: solidcore::surface::surface::Surface,
    bx: &Extent2d,
    kept: &[Vec<Piece>],
    face_copies: Vec<Attribute>,
    was: Option<OrientType>,
) -> core::result::Result<FaceId, Fail> {
    let mut nf = face_from_surface(k, target, env, surface, bx)?;
    for hole in kept {
        nf = rebuild_hole(k, target, env, nf, hole)?;
    }
    add_attributes(target, EntityRef::Face(nf), face_copies)?;
    // The face faces the way the old one did.
    if let (Some(was), Some(now)) = (was, face2::upward_faceuse(target, nf)?) {
        if was != eu1::faceuse_orientation(target, now)? {
            flip_face_orientation(target, nf)?;
        }
    }
    Ok(nf)
}
