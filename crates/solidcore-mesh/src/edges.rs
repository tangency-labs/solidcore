//! Edge polylines: every B-rep edge divided once, adaptively, and shared by its faces.

use std::collections::{BTreeMap, BTreeSet};

use solidcore::curve::curve::Curve;
use solidcore::kernel::Kernel;
use solidcore::math::{Point2, Point3, Vec3};
use solidcore::surface::surface::Surface;
use solidcore::topology::brep::p1 as bl;
use solidcore::topology::edge as ed;
use solidcore::topology::edgeuse::p1 as eu1;
use solidcore::topology::Brep;
use solidcore::{EdgeId, EntityRef, FaceId, Result, VertexId};

use crate::geom::{angle, curve_eval, invalid, segment_distance, surface_normal};
use crate::mesh::{EdgePolyline, Failure};
use crate::options::{EdgeParams, FaceParams, TessellateOptions};

/// Most times one segment is halved.
pub const MAX_EDGE_DEPTH: u32 = 30;
/// Most points of one edge polyline; more is an error.
pub const MAX_EDGE_POINTS: usize = 1 << 20;
/// Edge segments smaller than this fraction of the chord tolerance are not split.
const MIN_SPLIT_FRACTION: f64 = 1e-3;
/// Most pieces one [`SegmentSplit`] divides a segment into.
pub const MAX_SPLIT_PIECES: u32 = 64;

/// A boundary segment a face asks to have divided ([`EdgeSet::split_segments`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SegmentSplit {
    /// The mesh vertices of two neighbouring points of one edge polyline (sorted).
    pub segment: [u32; 2],
    /// The equal parameter pieces it should become (at most [`MAX_SPLIT_PIECES`]).
    pub pieces: u32,
}

/// The edge polylines of a B-rep and the mesh vertices they hold: the shared input of
/// every face ([`crate::tessellate_face`]). Its positions are the mesh's first vertices:
/// the B-rep vertices, then the edges' interior points.
#[derive(Debug, Clone, Default)]
pub struct EdgeSet {
    pub(crate) positions: Vec<Point3>,
    pub(crate) vertex_index: BTreeMap<VertexId, u32>,
    pub(crate) edges: Vec<EdgePolyline>,
    pub(crate) index: BTreeMap<EdgeId, usize>,
    pub(crate) failures: Vec<Failure>,
    pub(crate) model_size: f64,
    /// Each edge's segment limit at the last division.
    pub(crate) max_segments: BTreeMap<EdgeId, u32>,
}

impl EdgeSet {
    /// The mesh vertex positions so far.
    #[must_use]
    pub fn positions(&self) -> &[Point3] {
        &self.positions
    }

    /// The edge polylines, in the B-rep's edge order.
    #[must_use]
    pub fn polylines(&self) -> &[EdgePolyline] {
        &self.edges
    }

    /// The polyline of edge `e`, if it was tessellated.
    #[must_use]
    pub fn polyline(&self, e: EdgeId) -> Option<&EdgePolyline> {
        self.index.get(&e).and_then(|i| self.edges.get(*i))
    }

    /// Edges that could not be divided.
    #[must_use]
    pub fn failures(&self) -> &[Failure] {
        &self.failures
    }

    /// The model size relative tolerances resolve against (bounding-box diagonal).
    #[must_use]
    pub fn model_size(&self) -> f64 {
        self.model_size
    }

    pub(crate) fn position(&self, i: u32) -> Result<Point3> {
        self.positions
            .get(i as usize)
            .copied()
            .ok_or_else(|| invalid("mesh vertex out of range"))
    }

    pub(crate) fn push(&mut self, p: Point3) -> Result<u32> {
        push_position(&mut self.positions, p)
    }

    /// Divides (or, after an earlier level, further divides) every edge of `b` by
    /// `opts`. A finer level only adds points: every earlier point stays, with its id.
    pub(crate) fn refine(&mut self, k: &Kernel, b: &Brep, opts: &TessellateOptions) -> Result<()> {
        let size = self.model_size;
        for e in bl::edges(b, None)? {
            let faces = edge_faces(b, e)?;
            let params: Vec<FaceParams> = faces
                .iter()
                .map(|f| opts.face_params(f.face, size))
                .collect();
            let p = opts.edge_params(e, size, &params);
            self.max_segments.insert(e, p.max_segments);
            let slot = self.index.get(&e).copied();
            let result = match slot.and_then(|i| self.edges.get(i)).cloned() {
                Some(poly) => self.refine_edge(k, b, &faces, &p, poly),
                None => self.first_edge(k, b, e, &faces, &p),
            };
            match (result, slot) {
                (Ok(poly), Some(i)) => {
                    if let Some(s) = self.edges.get_mut(i) {
                        *s = poly;
                    }
                }
                (Ok(poly), None) => {
                    self.index.insert(e, self.edges.len());
                    self.edges.push(poly);
                }
                (Err(error), _) => {
                    let entity = EntityRef::Edge(e);
                    if !self.failures.iter().any(|f| f.entity == entity) {
                        self.failures.push(Failure { entity, error });
                    }
                }
            }
        }
        Ok(())
    }

    /// Splits each polyline segment of `segments` (a pair of neighbouring points of one
    /// edge polyline, in either order; the most pieces asked when one is listed twice)
    /// into equal parameter pieces, keeping every other point; returns the edges
    /// changed, in id order. Pairs that are no polyline segment are left alone, and an
    /// edge stops gaining points at its segment limit.
    ///
    /// # Errors
    ///
    /// The curve evaluation's errors.
    pub fn split_segments(
        &mut self,
        k: &Kernel,
        b: &Brep,
        segments: &[SegmentSplit],
    ) -> Result<Vec<EdgeId>> {
        let mut wanted: BTreeMap<[u32; 2], u32> = BTreeMap::new();
        for r in segments {
            let [x, y] = r.segment;
            let key = if x < y { [x, y] } else { [y, x] };
            let n = wanted.entry(key).or_insert(0);
            *n = (*n).max(r.pieces.min(MAX_SPLIT_PIECES));
        }
        // (edge index, segment index, pieces) of every wanted segment
        let mut at: BTreeMap<usize, Vec<(usize, u32)>> = BTreeMap::new();
        for (ei, poly) in self.edges.iter().enumerate() {
            if poly.degenerate {
                continue;
            }
            for (si, w) in poly.vertices.windows(2).enumerate() {
                if let [x, y] = w {
                    let key = if x < y { [*x, *y] } else { [*y, *x] };
                    if let Some(n) = wanted.get(&key) {
                        at.entry(ei).or_default().push((si, *n));
                    }
                }
            }
        }
        let mut changed = Vec::new();
        for (ei, mut segs) in at {
            let Some(poly) = self.edges.get(ei).cloned() else {
                continue;
            };
            let rec = b.edge(poly.edge)?;
            let curve = rec
                .curve
                .as_ref()
                .ok_or_else(|| invalid("edge without curve"))?;
            let mut poly = poly;
            let limit = self
                .max_segments
                .get(&poly.edge)
                .copied()
                .unwrap_or(u32::MAX) as usize;
            // back to front, so earlier indices stay valid
            segs.sort_unstable_by(|a, c| c.cmp(a));
            let before = poly.params.len();
            for (si, n) in segs {
                let (Some(&t0), Some(&t1)) = (poly.params.get(si), poly.params.get(si + 1)) else {
                    continue;
                };
                // the new points in order after point si
                for i in 1..n {
                    let tm = t0 + (t1 - t0) * f64::from(i) / f64::from(n);
                    let prev = poly.params.get(si + i as usize - 1).copied().unwrap_or(t0);
                    if !(tm > prev && tm < t1)
                        || poly.params.len() >= MAX_EDGE_POINTS
                        || poly.params.len() > limit
                    {
                        break;
                    }
                    let (p, _) = curve_eval(k, curve, tm, false)?;
                    let v = self.push(p)?;
                    poly.params.insert(si + i as usize, tm);
                    poly.vertices.insert(si + i as usize, v);
                }
            }
            if poly.params.len() == before {
                continue;
            }
            changed.push(poly.edge);
            if let Some(slot) = self.edges.get_mut(ei) {
                *slot = poly;
            }
        }
        changed.sort_unstable();
        Ok(changed)
    }

    fn vertex_slot(&mut self, b: &Brep, v: VertexId) -> Result<u32> {
        if let Some(i) = self.vertex_index.get(&v) {
            return Ok(*i);
        }
        let i = self.push(b.vertex(v)?.point)?;
        self.vertex_index.insert(v, i);
        Ok(i)
    }

    /// The first division of edge `e`.
    fn first_edge(
        &mut self,
        k: &Kernel,
        b: &Brep,
        e: EdgeId,
        faces: &[EdgeFace<'_>],
        p: &EdgeParams,
    ) -> Result<EdgePolyline> {
        let (v0, v1) = ed::start_end_vertices(b, e)?;
        let (Some(v0), Some(v1)) = (v0, v1) else {
            return Err(invalid("edge without vertices"));
        };
        let i0 = self.vertex_slot(b, v0)?;
        let i1 = self.vertex_slot(b, v1)?;
        let rec = b.edge(e)?;
        let ivl = rec.interval;
        let (t0, t1) = (ivl.min(), ivl.max());
        let degenerate = EdgePolyline {
            edge: e,
            params: vec![t0, t1],
            vertices: vec![i0, i1],
            degenerate: true,
        };
        let Some(curve) = rec.curve.as_ref() else {
            return if v0 == v1 {
                Ok(degenerate)
            } else {
                Err(invalid("edge without curve"))
            };
        };
        if !(t0.is_finite() && t1.is_finite() && t0 < t1) {
            return if v0 == v1 {
                Ok(degenerate)
            } else {
                Err(invalid("edge with an empty interval"))
            };
        }
        // a closed edge whose curve stays at its vertex is a pole
        if v0 == v1 {
            let pv = b.vertex(v0)?.point;
            let tol = p.chord.min(rec.tolerance.max(0.0) * 10.0).max(1e-12);
            let mut collapsed = true;
            for i in 0..=16 {
                let t = ivl.evaluate(f64::from(i) / 16.0);
                if curve_eval(k, curve, t, i == 16)?.0.distance(pv) > tol {
                    collapsed = false;
                    break;
                }
            }
            if collapsed {
                return Ok(degenerate);
            }
        }
        let n0 = p.min_segments.max(if v0 == v1 { 2 } else { 1 });
        let mut params = Vec::new();
        let mut vertices = Vec::new();
        for i in 0..=n0 {
            params.push(ivl.evaluate(f64::from(i) / f64::from(n0)));
            vertices.push(if i == 0 {
                i0
            } else if i == n0 {
                i1
            } else {
                u32::MAX
            });
        }
        // the inner starting points become vertices as they are kept
        for (t, v) in params.iter().zip(vertices.iter_mut()) {
            if *v == u32::MAX {
                let (pt, _) = curve_eval(k, curve, *t, false)?;
                *v = self.push(pt)?;
            }
        }
        let poly = EdgePolyline {
            edge: e,
            params,
            vertices,
            degenerate: false,
        };
        self.split(k, curve, (t0, t1), faces, p, &poly)
    }

    /// A further division of an edge already divided.
    fn refine_edge(
        &mut self,
        k: &Kernel,
        b: &Brep,
        faces: &[EdgeFace<'_>],
        p: &EdgeParams,
        poly: EdgePolyline,
    ) -> Result<EdgePolyline> {
        if poly.degenerate {
            return Ok(poly);
        }
        let rec = b.edge(poly.edge)?;
        let curve = rec
            .curve
            .as_ref()
            .ok_or_else(|| invalid("edge without curve"))?;
        let ivl = rec.interval;
        self.split(k, curve, (ivl.min(), ivl.max()), faces, p, &poly)
    }

    /// Splits every segment of `poly` until it meets `p`, keeping every point.
    fn split(
        &mut self,
        k: &Kernel,
        curve: &Curve,
        interval: (f64, f64),
        faces: &[EdgeFace<'_>],
        p: &EdgeParams,
        poly: &EdgePolyline,
    ) -> Result<EdgePolyline> {
        let n = poly.params.len();
        let mut params = Vec::with_capacity(n);
        let mut vertices = Vec::with_capacity(n);
        let mut segments = n.saturating_sub(1);
        let last = n.saturating_sub(1);
        // the existing points, with their tangents
        let mut ends = Vec::with_capacity(n);
        for (i, (t, v)) in poly.params.iter().zip(&poly.vertices).enumerate() {
            let (_, d) = curve_eval(k, curve, *t, i == last)?;
            ends.push((*t, self.position(*v)?, d, *v));
        }
        let Some(&(t_first, _, _, v_first)) = ends.first() else {
            return Err(invalid("edge polyline without points"));
        };
        params.push(t_first);
        vertices.push(v_first);
        let mut stack: Vec<Segment> = Vec::new();
        for w in ends.windows(2) {
            let [(ta, pa, da, _), (tb, pb, db, vb)] = w else {
                continue;
            };
            stack.push(Segment {
                t0: *ta,
                t1: *tb,
                p0: *pa,
                p1: *pb,
                d0: *da,
                d1: *db,
                v1: Some(*vb),
                mid: None,
                depth: 0,
            });
            // process this original segment fully before the next, in curve order
            while let Some(s) = stack.pop() {
                let (tm, mid) = if let Some(m) = s.mid {
                    (m.0, (m.1, m.2))
                } else {
                    let tm = (s.t0 + s.t1) / 2.0;
                    let (pm, dm) = curve_eval(k, curve, tm, false)?;
                    (tm, (pm, dm))
                };
                // a segment far below the chord tolerance is left alone (a tangent or
                // normal still turning there is a singular point's, as at a fillet's
                // collapsed corner)
                let size =
                    s.p0.distance(s.p1)
                        .max(s.p0.distance(mid.0))
                        .max(s.p1.distance(mid.0));
                let can = s.depth < MAX_EDGE_DEPTH
                    && segments < p.max_segments as usize
                    && s.p0.distance(s.p1) >= 2.0 * p.min_len
                    && size >= MIN_SPLIT_FRACTION * p.chord;
                let quarters = if can {
                    Some(Self::needs_split(k, curve, &s, mid.0, interval, faces, p)?)
                } else {
                    None
                };
                if let Some(Some((q1, q3))) = quarters {
                    segments += 1;
                    stack.push(Segment {
                        t0: tm,
                        t1: s.t1,
                        p0: mid.0,
                        p1: s.p1,
                        d0: mid.1,
                        d1: s.d1,
                        v1: s.v1,
                        mid: Some(q3),
                        depth: s.depth + 1,
                    });
                    stack.push(Segment {
                        t0: s.t0,
                        t1: tm,
                        p0: s.p0,
                        p1: mid.0,
                        d0: s.d0,
                        d1: mid.1,
                        v1: None,
                        mid: Some(q1),
                        depth: s.depth + 1,
                    });
                    continue;
                }
                let v = match s.v1 {
                    Some(v) => v,
                    None => self.push(s.p1)?,
                };
                params.push(s.t1);
                vertices.push(v);
                if params.len() > MAX_EDGE_POINTS {
                    return Err(invalid("edge needs too many points"));
                }
            }
        }
        Ok(EdgePolyline {
            edge: poly.edge,
            params,
            vertices,
            degenerate: false,
        })
    }

    /// When segment `s` (midpoint `pm`) is outside the tolerances, its two quarter
    /// points (the midpoints of its halves); `None` when it is within them. The chord
    /// is checked at the quarter points and the midpoint; the curve's tangent and every
    /// face's normal must not turn more than the angle tolerance.
    #[allow(clippy::type_complexity)]
    fn needs_split(
        k: &Kernel,
        curve: &Curve,
        s: &Segment,
        pm: Point3,
        interval: (f64, f64),
        faces: &[EdgeFace<'_>],
        p: &EdgeParams,
    ) -> Result<Option<((f64, Point3, Vec3), (f64, Point3, Vec3))>> {
        let t1q = s.t0 + (s.t1 - s.t0) * 0.25;
        let t3q = s.t0 + (s.t1 - s.t0) * 0.75;
        let (p1q, d1q) = curve_eval(k, curve, t1q, false)?;
        let (p3q, d3q) = curve_eval(k, curve, t3q, false)?;
        let quarters = Some(((t1q, p1q, d1q), (t3q, p3q, d3q)));
        if let Some(l) = p.max_len {
            if s.p0.distance(s.p1) > l {
                return Ok(quarters);
            }
        }
        if angle(s.d0, s.d1) > p.angle {
            return Ok(quarters);
        }
        for q in [p1q, pm, p3q] {
            if segment_distance(q, s.p0, s.p1) > p.chord {
                return Ok(quarters);
            }
        }
        let len = interval.1 - interval.0;
        for face in faces {
            let normal = |t: f64| {
                let ti = face.trim.natural_interval();
                let q = face
                    .trim
                    .evaluate_point(k, ti.evaluate((t - interval.0) / len))
                    .ok()?;
                surface_normal(k, face.surface, Point2::new(q.x, q.y))
            };
            if let (Some(a), Some(c)) = (normal(s.t0), normal(s.t1)) {
                if angle(a, c) > p.angle {
                    return Ok(quarters);
                }
            }
        }
        Ok(None)
    }
}

pub(crate) fn push_position(positions: &mut Vec<Point3>, p: Point3) -> Result<u32> {
    let i =
        u32::try_from(positions.len()).map_err(|_| invalid("more than u32::MAX mesh vertices"))?;
    positions.push(p);
    Ok(i)
}

/// One pending piece of an edge's division; its midpoint is cached when a parent
/// evaluated it as a quarter point.
#[derive(Debug, Clone, Copy)]
struct Segment {
    t0: f64,
    t1: f64,
    p0: Point3,
    p1: Point3,
    d0: Vec3,
    d1: Vec3,
    /// The mesh vertex at `t1`, when it exists already.
    v1: Option<u32>,
    mid: Option<(f64, Point3, Vec3)>,
    depth: u32,
}

/// A face along an edge: its surface and the edge use's UV trim curve.
pub(crate) struct EdgeFace<'a> {
    pub face: FaceId,
    pub surface: &'a Surface,
    pub trim: &'a Curve,
}

/// The faces along edge `e` with a surface and a UV trim curve for it, once each.
pub(crate) fn edge_faces(b: &Brep, e: EdgeId) -> Result<Vec<EdgeFace<'_>>> {
    let mut seen = BTreeSet::new();
    let mut out = Vec::new();
    for eu in &b.edge(e)?.edgeuses {
        let Some(f) = eu1::face(b, *eu)? else {
            continue;
        };
        let (Some(trim), _) = eu1::uv_trim_curve(b, *eu)? else {
            continue;
        };
        let Some(surface) = b.face(f)?.surface.as_ref() else {
            continue;
        };
        if seen.insert(f) {
            out.push(EdgeFace {
                face: f,
                surface,
                trim,
            });
        }
    }
    Ok(out)
}
