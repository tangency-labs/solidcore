//! One face: its boundary loops in UV from the shared edge polylines, a constrained
//! Delaunay triangulation in scaled UV, grid points and refinement, kept between levels
//! of detail so each level only inserts points.

use std::collections::{BTreeMap, BTreeSet, VecDeque};

use spade::handles::{FixedFaceHandle, FixedUndirectedEdgeHandle, FixedVertexHandle, InnerTag};
use spade::{
    ConstrainedDelaunayTriangulation, HasPosition, PositionInTriangulation, Triangulation,
};

use solidcore::curve::curve::Curve;
use solidcore::geom::extent1d::Extent1d;
use solidcore::geom::extent2d::Extent2d;
use solidcore::kernel::Kernel;
use solidcore::math::{Point2, Point3, Vec3, INFINITE_PARAMETER};
use solidcore::nurbs::b_spline_surface_n_lib::p2::OrientType;
use solidcore::surface::surface::Surface;
use solidcore::topology::edgeuse::p1 as eu1;
use solidcore::topology::{Brep, EdgeuseId, FaceuseId};
use solidcore::{FaceId, Result};

use crate::backend::{BackendFace, FaceBackend, FaceMeshView};
use crate::edges::{EdgeSet, SegmentSplit, MAX_SPLIT_PIECES};
use crate::geom::{
    angle, closes, foot, invalid, invert_point, power_of_two, segment_distance, surface_normal,
    surface_point, surface_scales,
};
use crate::grid::{bounded, grid_lines, MAX_CELL_ASPECT, SAMPLE_GRID};
use crate::mesh::EdgePolyline;
use crate::nurbs::Nurbs;
use crate::options::FaceParams;

/// A boundary point's trim-curve UV is polished when its surface point is farther than
/// this fraction of the chord tolerance from the edge point.
const POLISH_FRACTION: f64 = 0.1;
/// A refinement point keeps this fraction of its triangle's (or edge's) size from the
/// boundary.
const CLEARANCE_FRACTION: f64 = 0.1;
/// Triangles smaller than this fraction of the chord tolerance are not split.
const MIN_SPLIT_FRACTION: f64 = 1e-3;
/// A triangle whose doubled area is below this fraction of its longest edge squared is
/// a sliver: its normal is not compared with the surface's.
const SLIVER: f64 = 0.02;
/// The B-spline fast path is used when it agrees with the surface's own evaluation
/// within this fraction of the face's size.
const FAST_PATH_AGREEMENT: f64 = 1e-12;
/// A finer level inserts its grid when it has this many times the cells of the last
/// grid inserted.
const REGRID_FACTOR: usize = 4;
/// Refinement points closer than this fraction of the grid step to a point already
/// there are dropped.
const MIN_POINT_GAP: f64 = 1e-6;
/// A triangle across a fold of its surface is split while larger than this many
/// chord tolerances.
const FOLD_SIZE: f64 = 4.0;
/// Most UV samples of one degenerate (pole) segment.
const MAX_POLE_SEGMENTS: u32 = 1024;
/// Most generations of refinement points below a boundary or grid point: a criterion
/// still failing there (a singular point, a kink in the surface) stops refining.
pub const MAX_LOCAL_DEPTH: u8 = 48;
/// Largest loop a face with no area in UV is ear-clipped from, in points.
const MAX_EAR_LOOP: usize = 256;
/// An interior point closer than this fraction of the face's scaled UV extent to a
/// point already there is not inserted (spade's rounded geometry fails on such pairs).
const NEAR_DUPLICATE: f64 = 1e-12;
/// Two boundary segments crossing in UV are split at the crossing only when the sine of
/// their angle is at least this: the rounded crossing point of a grazing pair can land
/// on the wrong side of an edge and corrupt the triangulation, so such a face fails.
const MIN_CROSSING_SINE: f64 = 1e-3;
/// ... and when the crossing is at least this fraction of each segment from its ends.
const CROSSING_MARGIN: f64 = 1e-6;

/// One point of a face's boundary loop in UV.
#[derive(Debug, Clone, Copy)]
struct LoopPoint {
    uv: Point2,
    xyz: Point3,
    vertex: u32,
    /// The segment to the next point collapses to this point's vertex (a pole).
    pole_next: bool,
}

/// A boundary point as inserted: UV, position and mesh vertex.
type BoundaryPoint = (Point2, Point3, u32);

/// A face's loops and, for each, the direction it goes once around a period in, when
/// it does ([`FaceState::face_loops`]).
type FaceLoops = (Vec<Vec<LoopPoint>>, Vec<Option<bool>>);

/// How a loop meets itself where it closes, or two edge uses meet: at one UV point, at
/// a pole, or one period away on a closed surface.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Junction {
    /// Within the tolerance in UV.
    Same,
    /// Far apart in UV, but the surface stays at the vertex between them.
    Pole,
    /// One whole number of periods away: the shift `(du, dv)` that closes the gap.
    Period(f64, f64),
}

/// The plane a periodic face bounded by loops that each go once around the period (no
/// seam edge) is triangulated in: the periodic coordinate `a` turns about the origin
/// (one period one turn), the other coordinate `b` runs out along the radius. The face
/// then is an annulus there, and no seam is needed.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Annulus {
    /// The periodic coordinate is u (else v).
    along_u: bool,
    /// Its period.
    period: f64,
    /// The `a` at angle zero.
    ref_a: f64,
    /// The `b` at the inner radius.
    ref_b: f64,
    /// The span of `b` over the loops (the inner radius is at least this, scaled, so
    /// the outer radius is at most twice the inner).
    span_b: f64,
}

impl Annulus {
    /// `(a, b)` of `uv` and their scale factors.
    fn split(&self, uv: Point2, scale: (f64, f64)) -> (f64, f64, f64, f64) {
        if self.along_u {
            (uv.x, uv.y, scale.0, scale.1)
        } else {
            (uv.y, uv.x, scale.1, scale.0)
        }
    }

    /// The inner radius: one period around at the speed of `a`, at least the scaled span
    /// of `b`.
    fn inner(&self, sa: f64, sb: f64) -> f64 {
        (sa * self.period / core::f64::consts::TAU).max(self.span_b * sb)
    }

    /// The turn of one unit of `a`: clockwise for u, counter-clockwise for v, so the
    /// map keeps the orientation of (u, v).
    fn turn(&self) -> f64 {
        let t = core::f64::consts::TAU / self.period;
        if self.along_u {
            -t
        } else {
            t
        }
    }

    /// The plane point of `uv` with the triangulation's `scale`.
    fn plane_point(&self, uv: Point2, scale: (f64, f64)) -> (f64, f64) {
        let (a, b, sa, sb) = self.split(uv, scale);
        let r = self.inner(sa, sb) + (b - self.ref_b) * sb;
        let th = (a - self.ref_a) * self.turn();
        (r * th.cos(), r * th.sin())
    }

    /// The UV of plane point `(x, y)` (its `a` within half a period of
    /// [`Annulus::ref_a`]).
    fn plane_uv(&self, x: f64, y: f64, scale: (f64, f64)) -> Point2 {
        let (sa, sb) = if self.along_u {
            (scale.0, scale.1)
        } else {
            (scale.1, scale.0)
        };
        let r = (x * x + y * y).sqrt();
        let a = self.ref_a + y.atan2(x) / self.turn();
        let b = self.ref_b + (r - self.inner(sa, sb)) / sb;
        if self.along_u {
            Point2::new(a, b)
        } else {
            Point2::new(b, a)
        }
    }
}

/// A vertex of a face's triangulation.
#[derive(Debug, Clone, Copy)]
struct FaceVertex {
    uv: Point2,
    xyz: Point3,
    /// The shared mesh vertex: a boundary point's from the start, an interior point's
    /// once it is output.
    vertex: Option<u32>,
    /// Unit outward normal, when the surface has one there.
    normal: Option<Vec3>,
    /// On a boundary loop (never removed by refinement).
    boundary: bool,
    /// Still in the triangulation.
    alive: bool,
    /// How far the point is off the surface at its UV: the B-rep's own edge/face gap
    /// for a boundary point, `0` for an interior one. A triangle may be off the surface
    /// by its corners' gap on top of the chord tolerance.
    gap: f64,
    /// Refinement generation: `0` for boundary and grid points, one more than the
    /// deepest corner of the triangle (or edge) a refinement point splits.
    depth: u8,
}

/// The triangulation's vertex: its position in scaled UV and its [`FaceVertex`] index.
#[derive(Debug, Clone, Copy)]
struct CdtVertex {
    pos: spade::Point2<f64>,
    id: usize,
}

impl HasPosition for CdtVertex {
    type Scalar = f64;

    fn position(&self) -> spade::Point2<f64> {
        self.pos
    }
}

/// The triangulation, located through a hierarchy (logarithmic point location,
/// deterministic: a vertex's levels follow from its index). Spade's
/// `add_constraint_and_split` appends its split vertices without telling the hierarchy
/// and its `remove` then corrupts it: vertices are never removed here (a new
/// triangulation is made instead), and after a split the hierarchy is not used
/// ([`FaceState::stale_hints`]) until the triangulation is made again.
type Cdt =
    ConstrainedDelaunayTriangulation<CdtVertex, (), (), (), spade::HierarchyHintGenerator<f64>>;

/// The mesh of one face, ready to join a [`crate::Mesh`] ([`crate::assemble`]).
#[derive(Debug, Clone, PartialEq)]
#[allow(clippy::struct_excessive_bools)] // independent flags
pub struct FaceMesh {
    /// The face.
    pub face: FaceId,
    /// Each local vertex's (u, v).
    pub uvs: Vec<Point2>,
    /// Each local vertex's position (a boundary vertex's is the edge polyline's point).
    pub positions: Vec<Point3>,
    /// Each local vertex's unit outward normal (zero where the surface has none).
    pub normals: Vec<Vec3>,
    /// Each local vertex's shared mesh vertex: boundary vertices have one, interior
    /// vertices get one when assembled.
    pub shared: Vec<Option<u32>>,
    /// Triangles over local vertices, counter-clockwise seen from outside.
    pub triangles: Vec<[u32; 3]>,
    /// Each local vertex's distance from the surface at its UV: the B-rep's own
    /// edge/face gap at a boundary vertex, `0` inside.
    pub gaps: Vec<f64>,
    /// Edge uses mapped into UV by point inversion.
    pub inverted_edgeuses: Vec<EdgeuseId>,
    /// Refinement stopped at a bound with work left.
    pub unconverged: bool,
    /// The boundary crossed itself in UV and was split there.
    pub split_boundary: bool,
    /// The face had no area in UV and was ear-clipped in 3D.
    pub ear_clipped: bool,
    /// Largest measured deviation of the triangles from the surface (when measured).
    pub deviation: Option<f64>,
    /// Triangles whose measured deviation exceeds the chord tolerance plus the largest
    /// gap of their corners.
    pub triangles_over: usize,
    /// The backend that placed the face's grid lines (`None`: the built-in grid).
    pub backend: Option<&'static str>,
    /// A backend took the face and failed; the built-in grid was used.
    pub backend_failed: bool,
    /// The backend's own measure of the face's deviation (when measured).
    pub backend_deviation: Option<f64>,
}

/// A face being meshed, kept between levels.
#[allow(clippy::struct_excessive_bools)] // independent state flags
pub(crate) struct FaceState {
    face: FaceId,
    /// The UV the triangulation is centred on.
    origin: Point2,
    /// UV scale factors of the triangulation (powers of two): the speed scaling with
    /// the metric of the grid's cells (see [`FaceState::refine`]).
    scale: (f64, f64),
    /// Powers of two near the mean `|Su|`, `|Sv|`.
    speed: (f64, f64),
    /// `+1` when the outward side is the surface normal's.
    sense: f64,
    /// Triangles counter-clockwise in UV face inward and are wound the other way.
    flip: bool,
    /// The surface's natural domain (UV searches stay in it).
    domain: Extent2d,
    /// Periods of a surface closed in u or v (over its natural domain).
    period: (Option<f64>, Option<f64>),
    nurbs: Option<Nurbs>,
    /// The model size (bounds the fast path's disagreement).
    size: f64,
    /// With [`FaceState::nurbs`] checked against the surface: the sign turning its
    /// `Su x Sv` outward. Points and normals then come from the B-spline directly.
    fast: Option<f64>,
    cdt: Cdt,
    vertices: Vec<FaceVertex>,
    /// The boundary loops as inserted (pole runs expanded).
    boundary: Vec<Vec<(Point2, u32)>>,
    /// Samples of each pole run, by its first point (vertex, u bits, v bits).
    pole_segments: BTreeMap<(u32, u64, u64), u32>,
    inverted: Vec<EdgeuseId>,
    unconverged: bool,
    split_boundary: bool,
    /// Ear-clipped triangles of a face with no area in UV.
    ears: Vec<[usize; 3]>,
    /// Triangles found within the tolerances at the current level, by sorted vertex ids.
    good: BTreeSet<[usize; 3]>,
    /// Interior edges found within the tolerances at the current level.
    good_edges: BTreeSet<[usize; 2]>,
    /// The tolerances [`FaceState::good`] and [`FaceState::good_edges`] were found
    /// for (a vertex id keeps its position, so they hold until those change).
    good_for: Option<FaceParams>,
    /// Boundary segments too long for a bad triangle on or next to them to be refined
    /// from inside: the edge polylines should be split there.
    requests: Vec<SegmentSplit>,
    /// The grid lines `(us, vs)` of the current level.
    lines: (Vec<f64>, Vec<f64>),
    /// The edges of the face's loops.
    edges_used: BTreeSet<solidcore::EdgeId>,
    /// A refinement point was refused at [`MAX_LOCAL_DEPTH`] this level.
    depth_capped: bool,
    /// Cells of the last grid inserted.
    grid_cells: usize,
    /// The backend that placed the face's grid lines at this level.
    backend: Option<&'static str>,
    /// The UV box of the loops at the last refinement.
    bbox: Extent2d,
    /// A backend took the face and failed; the triangulation meshed it instead.
    backend_failed: bool,
    /// The face's loops each go once around a period with no seam edge: it is
    /// triangulated in that annulus.
    annulus: Option<Annulus>,
    /// A split constraint added vertices the triangulation's hierarchy does not know:
    /// points are located and inserted from [`FaceState::last_vertex`] instead until
    /// the triangulation is made again.
    stale_hints: bool,
    /// The vertex last inserted (the hint while the hierarchy is stale).
    last_vertex: Option<FixedVertexHandle>,
    /// The face's UV domain is unbounded (an infinite surface's parameter range): the
    /// triangulation's origin and scale are taken from the loops on the first level.
    rescale: bool,
}

/// The outward side of face `f`: `+1` when it is the surface normal's side. A solid's
/// outward side faces its void region: the face use whose region is void looks out
/// along its orientation. A sheet (no or both sides void) takes its upward face use's.
pub(crate) fn outward_sense(b: &Brep, f: FaceId) -> Result<f64> {
    let up = b
        .face(f)?
        .faceuse
        .ok_or_else(|| invalid("face without face use"))?;
    let mate = b.faceuse(up)?.mate;
    let sign = |fu: FaceuseId| -> Result<f64> {
        Ok(if b.faceuse(fu)?.orientation == OrientType::Opposite {
            -1.0
        } else {
            1.0
        })
    };
    let void = |fu: FaceuseId| -> Result<Option<bool>> {
        let Some(shell) = b.faceuse(fu)?.shell else {
            return Ok(None);
        };
        let Some(region) = b.shell(shell)?.region else {
            return Ok(None);
        };
        Ok(Some(b.region(region)?.is_void() == 1))
    };
    if let Some(mate) = mate {
        match (void(up)?, void(mate)?) {
            (Some(false), Some(true)) => return sign(mate),
            (Some(true), Some(false)) => return sign(up),
            _ => {}
        }
    }
    sign(up)
}

/// The periods of surface `s` with natural domain `natural`, in u and v: where the
/// surface closes across `probe` in that direction ([`closes`] at `1e-9` of the model
/// size), its natural width; else, where it closes a whole turn (`2π`, `360`) or two,
/// four or eight natural widths on (a surface of revolution whose natural domain is
/// only part of the turn: its faces may lie past it), that.
fn periods(
    k: &Kernel,
    s: &Surface,
    natural: &Extent2d,
    probe: &Extent2d,
    model_size: f64,
) -> (Option<f64>, Option<f64>) {
    let tol = 1e-9 * model_size.max(1e-300);
    let along = |u: bool| {
        let width = if u {
            natural.u_max() - natural.u_min()
        } else {
            natural.v_max() - natural.v_min()
        };
        if !(width > 0.0 && width < INFINITE_PARAMETER) {
            return None;
        }
        let candidates = [
            width,
            core::f64::consts::TAU,
            360.0,
            2.0 * width,
            4.0 * width,
            8.0 * width,
        ];
        candidates.into_iter().find(|&t| {
            if t < width {
                return false;
            }
            let span = if u {
                Extent2d::from_bounds(
                    probe.u_min(),
                    probe.v_min(),
                    probe.u_min() + t,
                    probe.v_max(),
                )
            } else {
                Extent2d::from_bounds(
                    probe.u_min(),
                    probe.v_min(),
                    probe.u_max(),
                    probe.v_min() + t,
                )
            };
            closes(k, s, &span, u, tol)
        })
    };
    (along(true), along(false))
}

/// Whether a bound of `d` is at the unbounded-parameter sentinel (or beyond, or not a
/// number).
fn unbounded(d: &Extent2d) -> bool {
    [d.u_min(), d.u_max(), d.v_min(), d.v_max()]
        .iter()
        .any(|x| x.is_nan() || x.abs() >= INFINITE_PARAMETER)
}

fn surface_of(b: &Brep, f: FaceId) -> Result<&Surface> {
    b.face(f)?
        .surface
        .as_ref()
        .ok_or_else(|| invalid("face without surface"))
}

fn same_uv(a: Point2, b: Point2) -> bool {
    a.x.to_bits() == b.x.to_bits() && a.y.to_bits() == b.y.to_bits()
}

impl FaceState {
    /// A face not meshed yet.
    pub(crate) fn new(k: &Kernel, b: &Brep, f: FaceId, model_size: f64) -> Result<Self> {
        let surface = surface_of(b, f)?;
        let sense = outward_sense(b, f)?;
        // (a lazy analytic face reads the domain of the NURB its MakeNurb builds)
        let natural = surface.natural_uv_domain_with_nurb(k)?;
        let mut domain = b.face(f)?.uv_domain;
        if !(domain.u_min() < domain.u_max() && domain.v_min() < domain.v_max()) {
            domain = natural;
        }
        // an unbounded domain (an infinite cone or cylinder) says nothing about the
        // face's size: sampled far out it would set the scale from points thousands of
        // times farther than the face
        let rescale = unbounded(&domain) || unbounded(&natural);
        // (until then: the scale near the parameter origin in an unbounded direction)
        let bound = |x: f64, near: f64| {
            if x.abs() < INFINITE_PARAMETER {
                x
            } else {
                near
            }
        };
        let sampled = Extent2d::from_bounds(
            bound(domain.u_min(), -1.0),
            bound(domain.v_min(), -1.0),
            bound(domain.u_max(), 1.0),
            bound(domain.v_max(), 1.0),
        );
        let scale = surface_scales(k, surface, &sampled);
        // across the natural domain where it is bounded, near the parameter origin where
        // not (its far ends say nothing; the loops' box refines it, see
        // FaceState::rescale_to)
        let probe = Extent2d::from_bounds(
            bound(natural.u_min(), -1.0),
            bound(natural.v_min(), -1.0),
            bound(natural.u_max(), 1.0),
            bound(natural.v_max(), 1.0),
        );
        let period = periods(k, surface, &natural, &probe, model_size);
        Ok(Self {
            face: f,
            origin: domain.mid(),
            scale,
            speed: scale,
            sense,
            flip: false,
            domain: natural,
            period,
            nurbs: Nurbs::from_surface(k, surface),
            size: model_size,
            fast: None,
            cdt: Cdt::new(),
            vertices: Vec::new(),
            boundary: Vec::new(),
            pole_segments: BTreeMap::new(),
            inverted: Vec::new(),
            unconverged: false,
            split_boundary: false,
            ears: Vec::new(),
            good: BTreeSet::new(),
            good_edges: BTreeSet::new(),
            good_for: None,
            requests: Vec::new(),
            lines: (Vec::new(), Vec::new()),
            edges_used: BTreeSet::new(),
            depth_capped: false,
            grid_cells: 0,
            backend: None,
            backend_failed: false,
            bbox: Extent2d::new(),
            annulus: None,
            stale_hints: false,
            last_vertex: None,
            rescale,
        })
    }

    fn scaled(&self, uv: Point2) -> spade::Point2<f64> {
        if let Some(an) = &self.annulus {
            let (x, y) = an.plane_point(uv, self.scale);
            return spade::mitigate_underflow(spade::Point2::new(x, y));
        }
        spade::mitigate_underflow(spade::Point2::new(
            (uv.x - self.origin.x) * self.scale.0,
            (uv.y - self.origin.y) * self.scale.1,
        ))
    }

    fn unscaled(&self, p: spade::Point2<f64>) -> Point2 {
        if let Some(an) = &self.annulus {
            return an.plane_uv(p.x, p.y, self.scale);
        }
        Point2::new(
            self.origin.x + p.x / self.scale.0,
            self.origin.y + p.y / self.scale.1,
        )
    }

    /// The UV at barycentric weights `w` of the UVs `uv`, blended in the plane of the
    /// triangulation (on an annulus a triangle may straddle the period, where blending
    /// the UVs themselves would land on the far side).
    fn plane_blend(&self, uv: [Point2; 3], w: [f64; 3]) -> Point2 {
        let [p, q, r] = uv.map(|x| self.scaled(x));
        self.unscaled(spade::Point2::new(
            p.x * w[0] + q.x * w[1] + r.x * w[2],
            p.y * w[0] + q.y * w[1] + r.y * w[2],
        ))
    }

    /// The surface point at `uv` (from the B-spline directly when checked).
    fn point(&self, k: &Kernel, s: &Surface, uv: Point2) -> Result<Point3> {
        if let (Some(nb), Some(_)) = (&self.nurbs, self.fast) {
            if let Some(p) = nb.point(uv.x, uv.y).filter(|p| crate::geom::finite3(*p)) {
                return Ok(p);
            }
        }
        surface_point(k, s, uv)
    }

    /// Checks the B-spline fast path against the surface's own evaluation at a sample
    /// grid of `bbox` and turns it on when they agree (and `Su x Sv` gives the normal's
    /// side).
    fn check_fast_path(&mut self, k: &Kernel, s: &Surface, bbox: &Extent2d, size: f64) {
        self.fast = None;
        let Some(nb) = &self.nurbs else {
            return;
        };
        for j in 0..=4 {
            for i in 0..=4 {
                let uv = Point2::new(
                    bbox.u_min() + (bbox.u_max() - bbox.u_min()) * f64::from(i) / 4.0,
                    bbox.v_min() + (bbox.v_max() - bbox.v_min()) * f64::from(j) / 4.0,
                );
                let (Some((p, _, _)), Ok(q)) = (nb.eval(uv.x, uv.y), s.evaluate_point(k, uv))
                else {
                    return;
                };
                let d = p.distance(q);
                if d > FAST_PATH_AGREEMENT * size || d.is_nan() {
                    return;
                }
            }
        }
        self.fast = self.cross_sign(k, s).map(|c| c * self.sense);
    }

    fn outward(&self, k: &Kernel, s: &Surface, uv: Point2) -> Option<Vec3> {
        if let (Some(nb), Some(sign)) = (&self.nurbs, self.fast) {
            if let Some((_, su, sv)) = nb.eval(uv.x, uv.y) {
                if let Some(n) = unit_cross(su, sv) {
                    return Some(n * sign);
                }
            }
        }
        surface_normal(k, s, uv).map(|n| n * self.sense)
    }

    fn vertex(&self, h: FixedVertexHandle) -> Result<&FaceVertex> {
        let id = self.cdt.vertex(h).data().id;
        self.vertices
            .get(id)
            .ok_or_else(|| invalid("triangulation vertex out of range"))
    }

    fn id(&self, h: FixedVertexHandle) -> usize {
        self.cdt.vertex(h).data().id
    }

    /// Adds `v` to the vertex table and the triangulation (its position must be free).
    fn insert(&mut self, v: FaceVertex) -> Result<FixedVertexHandle> {
        let id = self.vertices.len();
        self.vertices.push(v);
        self.insert_existing(id)?
            .ok_or_else(|| invalid("two points at one UV position"))
    }

    /// The triangulation's `locate` of a finite position (`NoTriangulation` for any
    /// other: spade asserts on NaN), from [`FaceState::last_vertex`] while the hierarchy
    /// is stale.
    fn locate(&self, pos: spade::Point2<f64>) -> PositionInTriangulation {
        if !(pos.x.is_finite() && pos.y.is_finite()) {
            return PositionInTriangulation::NoTriangulation;
        }
        if self.stale_hints {
            let hint = self
                .last_vertex
                .unwrap_or_else(|| FixedVertexHandle::from_index(0));
            self.cdt.locate_with_hint(pos, hint)
        } else {
            self.cdt.locate(pos)
        }
    }

    /// Whether `pos` is on, or within [`NEAR_DUPLICATE`] of the face's scaled UV extent
    /// of, a vertex of the triangulation around it.
    fn near_duplicate(&self, pos: spade::Point2<f64>) -> bool {
        let b = self.bbox;
        let extent = ((b.u_max() - b.u_min()) * self.scale.0)
            .abs()
            .max(((b.v_max() - b.v_min()) * self.scale.1).abs());
        let gap = if extent.is_finite() {
            NEAR_DUPLICATE * extent
        } else {
            0.0
        };
        let near: Vec<spade::Point2<f64>> = match self.locate(pos) {
            PositionInTriangulation::OnVertex(_) => return true,
            PositionInTriangulation::OnFace(f) => self.cdt.face(f).positions().to_vec(),
            PositionInTriangulation::OnEdge(e) => self.cdt.directed_edge(e).positions().to_vec(),
            _ => Vec::new(),
        };
        near.iter()
            .any(|q| (q.x - pos.x).abs() <= gap && (q.y - pos.y).abs() <= gap)
    }

    /// Puts vertex `id` of the table into the triangulation; `None` (nothing inserted)
    /// when a vertex is already there or all but there ([`FaceState::near_duplicate`]).
    fn insert_existing(&mut self, id: usize) -> Result<Option<FixedVertexHandle>> {
        let uv = self
            .vertices
            .get(id)
            .ok_or_else(|| invalid("vertex out of range"))?
            .uv;
        let pos = self.scaled(uv);
        let vertex = CdtVertex { pos, id };
        spade::validate_vertex(&vertex).map_err(|_| invalid("UV point cannot be triangulated"))?;
        if self.near_duplicate(pos) {
            return Ok(None);
        }
        let inserted = if self.stale_hints {
            let hint = self
                .last_vertex
                .unwrap_or_else(|| FixedVertexHandle::from_index(0));
            self.cdt.insert_with_hint(vertex, hint)
        } else {
            self.cdt.insert(vertex)
        };
        let h = inserted.map_err(|_| invalid("UV point cannot be triangulated"))?;
        self.last_vertex = Some(h);
        if let Some(v) = self.vertices.get_mut(id) {
            v.alive = true;
        }
        Ok(Some(h))
    }

    /// Inserts an interior point at `uv` with its surface point and normal (`None`:
    /// a point is already there, nothing inserted).
    fn insert_interior(
        &mut self,
        k: &Kernel,
        s: &Surface,
        uv: Point2,
        xyz: Option<(Point3, Option<Vec3>)>,
        depth: u8,
    ) -> Result<Option<FixedVertexHandle>> {
        if self.near_duplicate(self.scaled(uv)) {
            return Ok(None);
        }
        let (xyz, normal) = match xyz {
            Some(x) => x,
            None => (self.point(k, s, uv)?, self.outward(k, s, uv)),
        };
        let id = self.vertices.len();
        self.vertices.push(FaceVertex {
            uv,
            xyz,
            vertex: None,
            normal,
            boundary: false,
            alive: false,
            gap: 0.0,
            depth,
        });
        let h = self.insert_existing(id)?;
        if h.is_none() {
            // not in the triangulation: drop it from the table again
            self.vertices.truncate(id);
        }
        Ok(h)
    }

    /// The handle of the boundary point at `uv` on mesh vertex `vertex`, inserting it
    /// when it is new.
    fn boundary_handle(
        &mut self,
        k: &Kernel,
        s: &Surface,
        uv: Point2,
        xyz: Point3,
        vertex: u32,
    ) -> Result<FixedVertexHandle> {
        if let PositionInTriangulation::OnVertex(h) = self.locate(self.scaled(uv)) {
            let v = self.vertex(h)?;
            return if v.vertex == Some(vertex) && v.boundary {
                Ok(h)
            } else {
                Err(invalid("two boundary points at one UV position"))
            };
        }
        let normal = self.outward(k, s, uv);
        let gap = surface_gap(k, s, uv, xyz);
        self.insert(FaceVertex {
            uv,
            xyz,
            vertex: Some(vertex),
            normal,
            boundary: true,
            alive: true,
            gap,
            depth: 0,
        })
    }

    /// Joins two boundary vertices by a constraint; where it would cross a constraint
    /// already there, both are split at the crossing (the face then is closed there only
    /// geometrically).
    fn constrain(
        &mut self,
        k: &Kernel,
        s: &Surface,
        a: FixedVertexHandle,
        c: FixedVertexHandle,
    ) -> Result<()> {
        if a == c || self.cdt.exists_constraint(a, c) {
            return Ok(());
        }
        if !self.cdt.try_add_constraint(a, c).is_empty() || self.cdt.exists_constraint(a, c) {
            return Ok(());
        }
        // the constraints in the way: each must cross cleanly, else spade's split at a
        // rounded crossing point can corrupt the triangulation (and panic)
        let (pa, pc) = (self.cdt.vertex(a).position(), self.cdt.vertex(c).position());
        let mut crossings: Vec<[spade::Point2<f64>; 2]> = Vec::new();
        for e in self.cdt.get_conflicting_edges_between_vertices(a, c) {
            crossings.push(e.positions());
        }
        if crossings.is_empty() {
            return Err(invalid("boundary constraint cannot be added"));
        }
        if !crossings
            .iter()
            .all(|[p, q]| clean_crossing(pa, pc, *p, *q))
        {
            return Err(invalid(
                "boundary crosses itself in UV at a grazing angle or near a point",
            ));
        }
        let edges = self.cdt.add_constraint_and_split(a, c, |pos| CdtVertex {
            pos,
            id: usize::MAX,
        });
        if edges.is_empty() {
            return Err(invalid("boundary constraint cannot be added"));
        }
        self.split_boundary = true;
        self.stale_hints = true;
        // the split points become boundary vertices on the surface
        for e in edges {
            let to = self.cdt.directed_edge(e).to().fix();
            if self.cdt.vertex(to).data().id == usize::MAX {
                let pos = self.cdt.vertex(to).position();
                let uv = self.unscaled(pos);
                let xyz = surface_point(k, s, uv)?;
                let normal = self.outward(k, s, uv);
                let id = self.vertices.len();
                self.vertices.push(FaceVertex {
                    uv,
                    xyz,
                    vertex: None,
                    normal,
                    boundary: true,
                    alive: true,
                    gap: 0.0,
                    depth: 0,
                });
                self.cdt.vertex_data_mut(to).id = id;
            }
        }
        Ok(())
    }

    /// Meshes (or, after an earlier level, refines) the face to `p` against the edge
    /// polylines of `edges`. Only inserts points into what an earlier level made, except
    /// interior points the refined boundary leaves outside or against it.
    pub(crate) fn refine(
        &mut self,
        k: &Kernel,
        b: &Brep,
        edges: &EdgeSet,
        p: &FaceParams,
        backend: Option<&dyn FaceBackend>,
    ) -> Result<()> {
        let s = surface_of(b, self.face)?;
        if self.good_for != Some(*p) {
            self.good.clear();
            self.good_edges.clear();
            self.good_for = Some(*p);
        }
        self.unconverged = false;
        self.depth_capped = false;
        self.backend_failed = false;
        self.requests.clear();
        let (mut loops, mut wraps) = self.face_loops(k, b, s, edges, p, None)?;
        if self.rescale {
            // an unbounded domain: the origin, scale and periods of the loops' own
            // box, and the loops mapped again with them
            if self.rescale_to(k, s, &loops) {
                (loops, wraps) = self.face_loops(k, b, s, edges, p, None)?;
            }
        }
        if loops.is_empty() {
            return Err(invalid("face without boundary loops"));
        }
        let mut bbox = Extent2d::new();
        for q in loops.iter().flatten() {
            bbox.add_point(q.uv);
        }
        if !(bbox.u_min() <= bbox.u_max() && bbox.v_min() <= bbox.v_max()) {
            return Err(invalid("face loops have no UV points"));
        }
        let first = self.boundary.is_empty();
        let (mut us, mut vs) = if bbox.u_min() < bbox.u_max() && bbox.v_min() < bbox.v_max() {
            if first {
                // the B-spline fast path, before anything is evaluated
                self.check_fast_path(k, s, &bbox, self.size);
            }
            // a B-spline (on the fast path) is sampled once per knot span at least
            let samples =
                self.nurbs
                    .as_ref()
                    .filter(|_| self.fast.is_some())
                    .map_or(SAMPLE_GRID, |nb| {
                        nb.spans_in(true, bbox.u_min(), bbox.u_max())
                            .max(nb.spans_in(false, bbox.v_min(), bbox.v_max()))
                    });
            let eval = |uv: Point2| Ok((self.point(k, s, uv)?, self.outward(k, s, uv)));
            grid_lines(&eval, &bbox, self.speed, p, samples)?
        } else {
            (Vec::new(), Vec::new())
        };
        self.bbox = bbox;
        // a backend's lines, where denser than the sampled ones (those keep the angle
        // tolerance)
        self.backend = None;
        if let Some(be) =
            backend.filter(|_| bbox.u_min() < bbox.u_max() && bbox.v_min() < bbox.v_max())
        {
            let face = BackendFace {
                face: self.face,
                surface: s,
                bbox,
                chord: p.chord,
                angle: p.angle,
            };
            match be.grid_lines(k, &face) {
                Some(Ok((bu, bv))) => {
                    self.backend = Some(be.name());
                    if bu.len() > us.len() {
                        us = bu;
                    }
                    if bv.len() > vs.len() {
                        vs = bv;
                    }
                }
                // the sampled grid for a face the backend failed on
                Some(Err(_)) => self.backend_failed = true,
                None => {}
            }
        }
        // whichever source they came from, the lines are bounded like the built-in ones
        let (us, vs) = bounded(us, vs, &bbox, p.max_triangles);
        let step =
            |lines: &[f64], lo: f64, hi: f64, sc: f64| (hi - lo) * sc / (lines.len() + 1) as f64;
        let has_pole = loops.iter().flatten().any(|q| q.pole_next);
        // loops each once around a period (no seam edge): the face is an annulus
        let annulus = self.annulus_of(&wraps, &bbox, has_pole)?;
        let new_annulus = annulus.map(|a| a.along_u) != self.annulus.map(|a| a.along_u);
        if new_annulus {
            self.annulus = annulus;
        }
        let mut metric = self.speed;
        if !has_pole {
            // the triangulation's metric: UV scaled so the grid cells the tolerances
            // ask for are about square (at most MAX_CELL_ASPECT to one), so Delaunay
            // triangles are long where the surface is straight (along a ruling) and
            // refinement splits across it. A face with a pole keeps the speed scaling
            // (near a pole the grid's cells say nothing about the surface).
            let hu = step(&us, bbox.u_min(), bbox.u_max(), self.speed.0);
            let hv = step(&vs, bbox.v_min(), bbox.v_max(), self.speed.1);
            let h = hu.min(hv);
            if h > 0.0 && h.is_finite() {
                let (hu, hv) = (hu.min(MAX_CELL_ASPECT * h), hv.min(MAX_CELL_ASPECT * h));
                metric = (
                    power_of_two(self.speed.0 * h / hu),
                    power_of_two(self.speed.1 * h / hv),
                );
            }
        }
        // a finer level whose cells have another shape triangulates its points (all
        // kept) again in the new metric
        let ratio = |m: (f64, f64)| m.0 / m.1;
        let remetric = !first && {
            let r = ratio(metric) / ratio(self.scale);
            !(0.5..=2.0).contains(&r)
        };
        if first || remetric {
            self.scale = metric;
        }
        let scaled_step = step(&us, bbox.u_min(), bbox.u_max(), self.scale.0).min(step(
            &vs,
            bbox.v_min(),
            bbox.v_max(),
            self.scale.1,
        ));

        // the boundary, pole runs expanded
        let expanded = self.expand_poles(&loops, scaled_step);

        if first || remetric || new_annulus || !self.update_boundary(k, s, &loops, &expanded)? {
            self.rebuild(k, s, &loops, &expanded)?;
        }
        self.boundary = expanded
            .iter()
            .map(|l| l.iter().map(|q| (q.0, q.2)).collect())
            .collect();
        self.flip = self.uv_ccw_faces_inward(k, s, &bbox)?;
        if !first {
            self.remove_stale_interior(k, s)?;
        }
        // split constraints left the hierarchy behind: a new triangulation of the same
        // points and constraints (no crossings now) has a whole one
        if self.stale_hints {
            self.retriangulate(k, s, &BTreeSet::new())?;
        }

        // grid points, away from the boundary and the points already there; a finer
        // level adds its grid only when it is much denser than the last one inserted
        // (else it refines the earlier triangles: a grid on top of a similar one would
        // only crowd it)
        let cells = (us.len() + 1) * (vs.len() + 1);
        if first || cells >= REGRID_FACTOR * self.grid_cells {
            self.insert_grid(k, s, &us, &vs, scaled_step, p.max_triangles)?;
            self.grid_cells = cells;
        }
        self.lines = (us, vs);

        // refinement
        let mut converged = false;
        let mut changed: Option<Vec<FixedVertexHandle>> = None;
        for _ in 0..p.max_passes {
            if self.cdt.num_inner_faces() >= p.max_triangles {
                break;
            }
            let candidates = self.refine_candidates(k, s, p, changed.as_deref())?;
            if candidates.is_empty() {
                converged = true;
                break;
            }
            let mut inserted = Vec::new();
            for (uv, depth) in candidates {
                if self.cdt.num_inner_faces() >= p.max_triangles {
                    break;
                }
                let pos = self.scaled(uv);
                let near = match self.locate(pos) {
                    PositionInTriangulation::OnFace(f) => self.cdt.face(f).vertices().to_vec(),
                    PositionInTriangulation::OnEdge(e)
                        if !self.cdt.is_constraint_edge(e.as_undirected()) =>
                    {
                        self.cdt.directed_edge(e).vertices().to_vec()
                    }
                    _ => continue,
                };
                // a point (all but) on one already there (two candidates of one pass
                // can meet) would make a sliver
                let gap = MIN_POINT_GAP * scaled_step;
                if near.iter().any(|v| {
                    let q = v.position();
                    (q.x - pos.x).abs() <= gap && (q.y - pos.y).abs() <= gap
                }) {
                    continue;
                }
                if let Some(h) = self.insert_interior(k, s, uv, None, depth)? {
                    inserted.push(h);
                }
            }
            let added = inserted.len();
            changed = Some(inserted);
            if added == 0 {
                converged = true;
                break;
            }
        }
        self.unconverged = !converged || self.depth_capped;

        // a face with no area in UV: clip ears off its one loop in 3D
        self.ears.clear();
        let inside = inside_faces(&self.cdt);
        if !self.cdt.inner_faces().any(|f| is_inside(&inside, f.fix())) {
            self.ear_clip(k, s)?;
        }
        Ok(())
    }

    /// The domains mapping edge uses into UV searches ([`MapContext`]): Newton feet
    /// within the natural domain widened by a period on each side of a closed
    /// direction, and a point with no seed first sampled over the natural domain with
    /// an unbounded direction replaced by the loops' box (once known, else `[-1, 1]`).
    fn map_domains(&self) -> (Extent2d, Extent2d) {
        let n = self.domain;
        let widen = |lo: f64, hi: f64, p: Option<f64>| match p {
            Some(p) if p > 0.0 && p.is_finite() => (lo - p, hi + p),
            _ => (lo, hi),
        };
        let (u0, u1) = widen(n.u_min(), n.u_max(), self.period.0);
        let (v0, v1) = widen(n.v_min(), n.v_max(), self.period.1);
        let search = Extent2d::from_bounds(u0, v0, u1, v1);
        let bb = self.bbox;
        let known = bb.u_min() < bb.u_max() && bb.v_min() < bb.v_max();
        let pick = |x: f64, from_box: f64, fallback: f64| {
            if x.abs() < INFINITE_PARAMETER {
                x
            } else if known {
                from_box
            } else {
                fallback
            }
        };
        let grid = Extent2d::from_bounds(
            pick(n.u_min(), bb.u_min(), -1.0),
            pick(n.v_min(), bb.v_min(), -1.0),
            pick(n.u_max(), bb.u_max(), 1.0),
            pick(n.v_max(), bb.v_max(), 1.0),
        );
        (search, grid)
    }

    /// On an unbounded domain ([`FaceState::rescale`]): takes the triangulation's origin
    /// and scale from the box of `loops`, and the periods from the surface over it (the
    /// natural domain's far ends say nothing about the face). `true` when it did (the
    /// caller maps the loops again).
    fn rescale_to(&mut self, k: &Kernel, s: &Surface, loops: &[Vec<LoopPoint>]) -> bool {
        self.rescale = false;
        let mut bb = Extent2d::new();
        for q in loops.iter().flatten() {
            bb.add_point(q.uv);
        }
        if !(bb.u_min() < bb.u_max() && bb.v_min() < bb.v_max()) {
            return false;
        }
        self.origin = bb.mid();
        self.speed = surface_scales(k, s, &bb);
        self.scale = self.speed;
        self.bbox = bb;
        // the natural domain where bounded, the loops' box where not
        let n = self.domain;
        let pick = |x: f64, y: f64| if x.abs() < INFINITE_PARAMETER { x } else { y };
        let probe = Extent2d::from_bounds(
            pick(n.u_min(), bb.u_min()),
            pick(n.v_min(), bb.v_min()),
            pick(n.u_max(), bb.u_max()),
            pick(n.v_max(), bb.v_max()),
        );
        self.period = periods(k, s, &n, &probe, self.size);
        true
    }

    /// A text report of how the face's loops map into UV, for diagnosing a face that
    /// does not mesh: the triangulation's origin, scale and periods, then the trace of
    /// [`Self::face_loops`] itself (the same start, seeds, end hints and joins as the
    /// mesh): per edge use how it was mapped, its ends and how it meets the previous
    /// one, and how each loop closes; last the outcome. Takes the scale of the loops'
    /// box on an unbounded domain, as the first level does.
    pub(crate) fn loop_report(
        &mut self,
        k: &Kernel,
        b: &Brep,
        edges: &EdgeSet,
        p: &FaceParams,
    ) -> Result<String> {
        use core::fmt::Write as _;
        let s = surface_of(b, self.face)?;
        if self.rescale {
            if let Ok((loops, _)) = self.face_loops(k, b, s, edges, p, None) {
                self.rescale_to(k, s, &loops);
            }
            self.rescale = false;
        }
        let mut out = String::new();
        let _ = writeln!(
            out,
            "face {} origin {:?} scale {:?} period {:?} natural {:?} chord {}",
            self.face.index(),
            self.origin,
            self.scale,
            self.period,
            self.domain,
            p.chord
        );
        let mut trace = Vec::new();
        let result = self.face_loops(k, b, s, edges, p, Some(&mut trace));
        for line in &trace {
            let _ = writeln!(out, "{line}");
        }
        match result {
            Ok((loops, wraps)) => {
                let _ = writeln!(out, " ok: {} loops, wraps {wraps:?}", loops.len());
            }
            Err(e) => {
                let _ = writeln!(out, " error: {e:?}");
            }
        }
        Ok(out)
    }

    /// The annulus a face whose loops `wraps` go around a period is triangulated in:
    /// `None` when no loop does; an error unless at least two do, all around the same
    /// closed direction, with no pole.
    fn annulus_of(
        &self,
        wraps: &[Option<bool>],
        bbox: &Extent2d,
        has_pole: bool,
    ) -> Result<Option<Annulus>> {
        let around: Vec<bool> = wraps.iter().flatten().copied().collect();
        let Some(&along_u) = around.first() else {
            return Ok(None);
        };
        let period = if along_u {
            self.period.0
        } else {
            self.period.1
        };
        let fail = || invalid("loop does not close in UV (periodic face without seam edge)");
        let Some(period) = period.filter(|p| *p > 0.0 && p.is_finite()) else {
            return Err(fail());
        };
        if has_pole || around.len() < 2 || around.iter().any(|d| *d != along_u) {
            return Err(fail());
        }
        // a later level keeps the map (its points stay where they were inserted)
        if let Some(an) = self.annulus.filter(|a| a.along_u == along_u) {
            return Ok(Some(an));
        }
        let (a_lo, a_hi, b_lo, b_hi) = if along_u {
            (bbox.u_min(), bbox.u_max(), bbox.v_min(), bbox.v_max())
        } else {
            (bbox.v_min(), bbox.v_max(), bbox.u_min(), bbox.u_max())
        };
        Ok(Some(Annulus {
            along_u,
            period,
            ref_a: (a_lo + a_hi) / 2.0,
            ref_b: b_lo,
            span_b: b_hi - b_lo,
        }))
    }

    /// The boundary loops of the face in UV, outer first, each closed implicitly, and
    /// for each whether it goes once around a period instead of closing (`Some(true)`
    /// around u, `Some(false)` around v).
    /// With `trace`, a line per edge use (how it was mapped, its ends, how it meets the
    /// previous one) and per loop closure is pushed as it is made ([`Self::loop_report`]).
    fn face_loops(
        &mut self,
        k: &Kernel,
        b: &Brep,
        s: &Surface,
        edges: &EdgeSet,
        p: &FaceParams,
        mut trace: Option<&mut Vec<String>>,
    ) -> Result<FaceLoops> {
        let domain = self.domain;
        let (search, grid) = self.map_domains();
        let scale = self.scale;
        // the domain's size, a direction it is unbounded in left out (else the face's
        // own box, else unbounded: no jump is too long)
        let extent = |lo: f64, hi: f64, sc: f64| {
            if lo.abs() >= INFINITE_PARAMETER || hi.abs() >= INFINITE_PARAMETER {
                0.0
            } else {
                (hi - lo) * sc
            }
        };
        let mut span = extent(domain.u_min(), domain.u_max(), scale.0).max(extent(
            domain.v_min(),
            domain.v_max(),
            scale.1,
        ));
        if span.is_nan() || span <= 0.0 {
            let bb = self.bbox;
            span = ((bb.u_max() - bb.u_min()) * scale.0).max((bb.v_max() - bb.v_min()) * scale.1);
        }
        if span.is_nan() || span <= 0.0 {
            span = f64::INFINITY;
        }
        let polish_above = POLISH_FRACTION * p.chord;
        let fu = b
            .face(self.face)?
            .faceuse
            .ok_or_else(|| invalid("face without face use"))?;
        let mut loops = Vec::new();
        let mut wraps = Vec::new();
        let mut inverted = Vec::new();
        let ctx = MapContext {
            k,
            s,
            domain: &search,
            grid: &grid,
            edges,
            polish_above,
            period: self.period,
        };
        for (li, lu) in eu1::faceuse_loopuses(b, fu, false, None)?
            .into_iter()
            .enumerate()
        {
            if let Some(t) = trace.as_deref_mut() {
                t.push(format!(" loop {li} (span {span})"));
            }
            let mut pts: Vec<LoopPoint> = Vec::new();
            // how far the loop has moved the trim curves' own UVs by whole periods (it
            // stays continuous across a seam it goes around)
            let mut offset = (0.0, 0.0);
            let eus = eu1::loopuse_edgeuses(b, lu)?;
            // every edge use mapped by its UV trim curve where that follows the edge
            let mut uses = Vec::with_capacity(eus.len());
            for eu in &eus {
                let e = eu1::edge(b, *eu)?.ok_or_else(|| invalid("edge use without edge"))?;
                self.edges_used.insert(e);
                let poly = edges
                    .polyline(e)
                    .ok_or_else(|| invalid("face uses an edge that failed"))?;
                let ivl = b.edge(e)?.interval;
                let reverse = eu1::orientation(b, *eu)? == OrientType::Opposite;
                let by_trim = match eu1::uv_trim_curve(b, *eu)?.0 {
                    Some(c) => ctx
                        .edgeuse_points(poly, &ivl, Some(c), reverse, None, None)
                        .ok(),
                    None => None,
                };
                uses.push((*eu, poly, ivl, reverse, by_trim));
            }
            // start at an edge use so mapped, so edge uses mapped by point inversion
            // continue from a known side of any seam
            let n = uses.len();
            let start = uses.iter().position(|u| u.4.is_some()).unwrap_or(0);
            for step in 0..n {
                let i = (start + step) % n;
                // the next edge use's start, when its trim curve maps it: where an
                // inverted edge use ends (its Newton feet are seeded towards it, so a
                // badly conditioned stretch near a pole does not wander off)
                // (taken as it is, moved only as the loop has moved the trim curves: the
                // trim curve puts it on the face's own side of a seam, which moving it
                // by periods towards the seed could only lose, as a seed half a period
                // away rounds either way)
                let end = match uses.get((i + 1) % n) {
                    Some((.., Some(next))) => next
                        .first()
                        .map(|q| (Point2::new(q.uv.x + offset.0, q.uv.y + offset.1), true)),
                    _ => None,
                };
                // else where the use's end vertex is found on the surface (the
                // seed's side of a seam is taken when it is used)
                let end = match (end, uses.get(i)) {
                    (Some(e), _) => Some(e),
                    (None, Some((_, poly, _, reverse, None))) => {
                        let last = if *reverse {
                            poly.vertices.first()
                        } else {
                            poly.vertices.last()
                        };
                        if let Some(v) = last {
                            invert_point(k, s, edges.position(*v)?, &grid, &search)
                                .ok()
                                .map(|uv| (uv, false))
                        } else {
                            None
                        }
                    }
                    _ => None,
                };
                let Some((eu, poly, ivl, reverse, by_trim)) = uses.get_mut(i) else {
                    continue;
                };
                let trimmed = by_trim.is_some();
                let seq = if let Some(seq) = by_trim.take() {
                    seq
                } else {
                    inverted.push(*eu);
                    let seed = pts.last().map(|q| q.uv);
                    ctx.edgeuse_points(poly, ivl, None, *reverse, seed, end)?
                };
                if let Some(t) = trace.as_deref_mut() {
                    let meets = match (pts.last(), seq.first()) {
                        (Some(l), Some(f)) => format!(
                            "{:?}",
                            junction(k, s, p.chord, *l, *f, span, scale, self.period)
                        ),
                        _ => "start".to_string(),
                    };
                    let ends = (seq.first(), seq.last());
                    t.push(format!(
                        "  eu {} {} {} n {} end-hint {end:?} first {:?} last {:?} meets {meets}",
                        eu.index(),
                        if *reverse { "rev" } else { "fwd" },
                        if trimmed { "trim" } else { "inverted" },
                        seq.len(),
                        ends.0.map(|q| (q.vertex, q.uv)),
                        ends.1.map(|q| (q.vertex, q.uv)),
                    ));
                }
                let moved = join(&mut pts, seq, k, s, p.chord, span, scale, self.period)?;
                if trimmed {
                    offset = moved;
                }
            }
            // close the loop: at its start, at a pole, or a period on (a loop once
            // around a closed surface with no seam edge)
            let mut wrap = None;
            if pts.len() > 1 {
                if let (Some(first), Some(last)) = (pts.first().copied(), pts.last().copied()) {
                    if first.vertex != last.vertex {
                        return Err(invalid("loop does not close"));
                    }
                    let closes = junction(k, s, p.chord, last, first, span, scale, self.period);
                    if let Some(t) = trace.as_deref_mut() {
                        t.push(format!("  closes {closes:?}"));
                    }
                    match closes? {
                        Junction::Pole => {
                            if let Some(l) = pts.last_mut() {
                                l.pole_next = true;
                            }
                        }
                        Junction::Same => {
                            pts.pop();
                        }
                        Junction::Period(du, dv) => {
                            // once around exactly (a loop wound more often, or around both
                            // ways, bounds no annulus)
                            let once = |d: f64, p: Option<f64>| {
                                d == 0.0 || p.is_some_and(|p| ((d.abs() / p) - 1.0).abs() < 1e-6)
                            };
                            if (du != 0.0) == (dv != 0.0)
                                || !once(du, self.period.0)
                                || !once(dv, self.period.1)
                            {
                                return Err(invalid(
                                    "loop does not close in UV (periodic face without seam edge)",
                                ));
                            }
                            pts.pop();
                            wrap = Some(du != 0.0);
                        }
                    }
                }
            }
            if pts.len() >= 2 {
                loops.push(pts);
                wraps.push(wrap);
            }
        }
        for eu in inverted {
            if !self.inverted.contains(&eu) {
                self.inverted.push(eu);
            }
        }
        Ok((loops, wraps))
    }

    /// The loops with every pole run sampled at about the grid spacing (a power of two
    /// of samples, never fewer than an earlier level's, so levels nest): `(uv, xyz,
    /// vertex)` per point.
    fn expand_poles(
        &mut self,
        loops: &[Vec<LoopPoint>],
        scaled_step: f64,
    ) -> Vec<Vec<(Point2, Point3, u32)>> {
        let mut out = Vec::with_capacity(loops.len());
        for lp in loops {
            let count = lp.len();
            let mut run = Vec::with_capacity(count);
            for (i, q) in lp.iter().enumerate() {
                run.push((q.uv, q.xyz, q.vertex));
                if !q.pole_next {
                    continue;
                }
                let Some(next) = lp.get((i + 1) % count) else {
                    continue;
                };
                let (a, c) = (self.scaled(q.uv), self.scaled(next.uv));
                let len = ((c.x - a.x).powi(2) + (c.y - a.y).powi(2)).sqrt();
                let need = if scaled_step > 0.0 && scaled_step.is_finite() {
                    ((len / scaled_step).ceil() as u32).clamp(1, MAX_POLE_SEGMENTS)
                } else {
                    1
                };
                let key = (q.vertex, q.uv.x.to_bits(), q.uv.y.to_bits());
                let segs = need
                    .next_power_of_two()
                    .max(self.pole_segments.get(&key).copied().unwrap_or(1));
                self.pole_segments.insert(key, segs);
                for j in 1..segs {
                    let w = f64::from(j) / f64::from(segs);
                    run.push((
                        Point2::new(
                            q.uv.x + (next.uv.x - q.uv.x) * w,
                            q.uv.y + (next.uv.y - q.uv.y) * w,
                        ),
                        q.xyz,
                        q.vertex,
                    ));
                }
            }
            out.push(run);
        }
        out
    }

    /// A new triangulation of the boundary; interior points of earlier levels that are
    /// still clear of it are put back.
    fn rebuild(
        &mut self,
        k: &Kernel,
        s: &Surface,
        _loops: &[Vec<LoopPoint>],
        expanded: &[Vec<(Point2, Point3, u32)>],
    ) -> Result<()> {
        let old: Vec<usize> = (0..self.vertices.len())
            .filter(|i| {
                self.vertices
                    .get(*i)
                    .is_some_and(|v| v.alive && !v.boundary)
            })
            .collect();
        self.cdt = Cdt::new();
        for v in &mut self.vertices {
            v.alive = false;
        }
        self.split_boundary = false;
        // every boundary point once (one UV position holds one mesh vertex), and the
        // loop segments between them; bulk loading is fast and stable
        let mut index: BTreeMap<(u64, u64), usize> = BTreeMap::new();
        let mut points: Vec<CdtVertex> = Vec::new();
        let mut segments: Vec<[usize; 2]> = Vec::new();
        for lp in expanded {
            let mut ids: Vec<usize> = Vec::with_capacity(lp.len());
            for (uv, xyz, vertex) in lp {
                // by position in the triangulation (on an annulus, UVs a period apart
                // are one point)
                let pos = self.scaled(*uv);
                let key = (pos.x.to_bits(), pos.y.to_bits());
                let i = if let Some(i) = index.get(&key) {
                    let id = points.get(*i).map_or(usize::MAX, |p| p.id);
                    if self.vertices.get(id).and_then(|v| v.vertex) != Some(*vertex) {
                        return Err(invalid("two boundary points at one UV position"));
                    }
                    *i
                } else {
                    let normal = self.outward(k, s, *uv);
                    let id = self.vertices.len();
                    self.vertices.push(FaceVertex {
                        uv: *uv,
                        xyz: *xyz,
                        vertex: Some(*vertex),
                        normal,
                        boundary: true,
                        alive: true,
                        gap: surface_gap(k, s, *uv, *xyz),
                        depth: 0,
                    });
                    spade::validate_vertex(&CdtVertex { pos, id })
                        .map_err(|_| invalid("UV point cannot be triangulated"))?;
                    index.insert(key, points.len());
                    points.push(CdtVertex { pos, id });
                    points.len() - 1
                };
                if ids.last() != Some(&i) {
                    ids.push(i);
                }
            }
            if ids.len() > 1 && ids.first() == ids.last() {
                ids.pop();
            }
            let n = ids.len();
            for i in 0..n {
                if let (Some(a), Some(c)) = (ids.get(i), ids.get((i + 1) % n)) {
                    if a != c {
                        segments.push([*a, *c]);
                    }
                }
            }
        }
        self.bulk_load(k, s, points, segments)?;
        if !old.is_empty() {
            // which points are inside, all decided on the boundary's triangulation (an
            // insertion makes faces the flood fill has not seen)
            let inside = inside_faces(&self.cdt);
            let keep: Vec<usize> = old
                .into_iter()
                .filter(|id| {
                    self.vertices
                        .get(*id)
                        .is_some_and(|v| self.clear_of_boundary(&inside, v.uv, 0.0))
                })
                .collect();
            for id in keep {
                self.insert_existing(id)?;
            }
        }
        Ok(())
    }

    /// A new triangulation of the boundary `points` and its `segments` (indices into
    /// `points`); a segment crossing another is added with splits.
    fn bulk_load(
        &mut self,
        k: &Kernel,
        s: &Surface,
        points: Vec<CdtVertex>,
        segments: Vec<[usize; 2]>,
    ) -> Result<()> {
        let ids: Vec<usize> = points.iter().map(|p| p.id).collect();
        let mut conflicts = Vec::new();
        self.cdt = Cdt::try_bulk_load_cdt(points, segments, |e| conflicts.push(e))
            .map_err(|_| invalid("UV boundary cannot be triangulated"))?;
        self.stale_hints = false;
        self.last_vertex = None;
        for [a, c] in conflicts {
            let (Some(ia), Some(ic)) = (ids.get(a), ids.get(c)) else {
                continue;
            };
            let (Some(ha), Some(hc)) = (self.handle_of(*ia), self.handle_of(*ic)) else {
                continue;
            };
            self.constrain(k, s, ha, hc)?;
        }
        Ok(())
    }

    /// A new triangulation of the points still in it but those in `drop`, with its
    /// constraint edges: how points leave the triangulation (spade's vertex removal
    /// is not used: it has corrupted the triangulation's hint and point location).
    fn retriangulate(&mut self, k: &Kernel, s: &Surface, drop: &BTreeSet<usize>) -> Result<()> {
        let mut index: BTreeMap<usize, usize> = BTreeMap::new();
        let mut points: Vec<CdtVertex> = Vec::new();
        let mut interior: Vec<usize> = Vec::new();
        for v in self.cdt.vertices() {
            let id = v.data().id;
            if drop.contains(&id) {
                continue;
            }
            if self.vertices.get(id).is_some_and(|x| x.boundary) {
                index.insert(id, points.len());
                points.push(*v.data());
            } else {
                interior.push(id);
            }
        }
        let mut segments: Vec<[usize; 2]> = Vec::new();
        for e in self.cdt.undirected_edges() {
            if !self.cdt.is_constraint_edge(e.fix()) {
                continue;
            }
            let [x, y] = e.vertices().map(|v| v.data().id);
            if let (Some(i), Some(j)) = (index.get(&x), index.get(&y)) {
                segments.push([*i, *j]);
            }
        }
        for id in drop {
            if let Some(v) = self.vertices.get_mut(*id) {
                v.alive = false;
            }
        }
        for id in &interior {
            if let Some(v) = self.vertices.get_mut(*id) {
                v.alive = false;
            }
        }
        self.bulk_load(k, s, points, segments)?;
        for id in interior {
            self.insert_existing(id)?;
        }
        Ok(())
    }

    /// Inserts the points the edges gained since the last level into the boundary,
    /// replacing each constraint they split. `false` when the boundary changed other
    /// than by insertion (the caller rebuilds).
    fn update_boundary(
        &mut self,
        k: &Kernel,
        s: &Surface,
        _loops: &[Vec<LoopPoint>],
        expanded: &[Vec<(Point2, Point3, u32)>],
    ) -> Result<bool> {
        if expanded.len() != self.boundary.len() {
            return Ok(false);
        }
        // per loop: for every old segment, the new points inside it
        let mut plans: Vec<Vec<(usize, Vec<BoundaryPoint>)>> = Vec::new();
        for (old, new) in self.boundary.iter().zip(expanded) {
            let mut plan = Vec::new();
            let mut j = 0;
            let mut pending: Vec<BoundaryPoint> = Vec::new();
            for (i, (ouv, ov)) in old.iter().enumerate() {
                loop {
                    let Some(q) = new.get(j) else {
                        return Ok(false);
                    };
                    j += 1;
                    if q.2 == *ov && same_uv(q.0, *ouv) {
                        break;
                    }
                    if i == 0 {
                        return Ok(false);
                    }
                    pending.push(*q);
                }
                if i > 0 && !pending.is_empty() {
                    plan.push((i - 1, std::mem::take(&mut pending)));
                }
            }
            // after the last old point: inside the closing segment
            pending.extend(new.iter().skip(j).copied());
            if !pending.is_empty() {
                plan.push((old.len().saturating_sub(1), pending));
            }
            plans.push(plan);
        }
        for (li, plan) in plans.into_iter().enumerate() {
            for (si, points) in plan {
                let Some(old) = self.boundary.get(li) else {
                    return Ok(false);
                };
                let n = old.len();
                let (Some(&(auv, av)), Some(&(cuv, cv))) = (old.get(si), old.get((si + 1) % n))
                else {
                    return Ok(false);
                };
                let ha = self.find_boundary(auv, av)?;
                let hc = self.find_boundary(cuv, cv)?;
                let (Some(ha), Some(hc)) = (ha, hc) else {
                    return Ok(false);
                };
                let Some(e) = self
                    .cdt
                    .get_edge_from_neighbors(ha, hc)
                    .map(|e| e.fix().as_undirected())
                else {
                    return Ok(false);
                };
                if !self.cdt.is_constraint_edge(e) {
                    return Ok(false);
                }
                self.cdt.remove_constraint_edge(e);
                let mut chain = vec![(auv, av)];
                for (uv, xyz, v) in &points {
                    self.boundary_handle(k, s, *uv, *xyz, *v)?;
                    chain.push((*uv, *v));
                }
                chain.push((cuv, cv));
                for w in chain.windows(2) {
                    let [(xa, va), (xc, vc)] = w else {
                        continue;
                    };
                    let (Some(h0), Some(h1)) =
                        (self.find_boundary(*xa, *va)?, self.find_boundary(*xc, *vc)?)
                    else {
                        return Ok(false);
                    };
                    self.constrain(k, s, h0, h1)?;
                }
            }
        }
        Ok(true)
    }

    /// The handle of vertex `id` of the table, by its position.
    fn handle_of(&self, id: usize) -> Option<FixedVertexHandle> {
        let uv = self.vertices.get(id)?.uv;
        match self.locate(self.scaled(uv)) {
            PositionInTriangulation::OnVertex(h) => Some(h),
            _ => None,
        }
    }

    fn find_boundary(&self, uv: Point2, vertex: u32) -> Result<Option<FixedVertexHandle>> {
        match self.locate(self.scaled(uv)) {
            PositionInTriangulation::OnVertex(h) if self.vertex(h)?.vertex == Some(vertex) => {
                Ok(Some(h))
            }
            _ => Ok(None),
        }
    }

    /// Removes interior points the refined boundary leaves outside, or closer to it
    /// than a tenth of their shortest edge.
    fn remove_stale_interior(&mut self, k: &Kernel, s: &Surface) -> Result<()> {
        let inside = inside_faces(&self.cdt);
        let mut stale = BTreeSet::new();
        for v in self.cdt.vertices() {
            let id = v.data().id;
            if self.vertices.get(id).is_none_or(|x| x.boundary) {
                continue;
            }
            let mut out = false;
            let mut shortest = f64::INFINITY;
            let mut near = f64::INFINITY;
            let pos = v.position();
            for e in v.out_edges() {
                match e.face().as_inner() {
                    Some(f) if is_inside(&inside, f.fix()) => {
                        near = near.min(constraint_clearance(&self.cdt, f.fix(), pos));
                    }
                    _ => out = true,
                }
                let q = e.to().position();
                shortest = shortest.min(((q.x - pos.x).powi(2) + (q.y - pos.y).powi(2)).sqrt());
            }
            if out || near < CLEARANCE_FRACTION * shortest {
                stale.insert(id);
            }
        }
        if stale.is_empty() {
            return Ok(());
        }
        self.retriangulate(k, s, &stale)
    }

    /// Inserts the grid points inside the loops that are clear of the boundary (0.3 of
    /// the step) and of the points already there (half the step), until the triangulation has
    /// `max_triangles` triangles.
    fn insert_grid(
        &mut self,
        k: &Kernel,
        s: &Surface,
        us: &[f64],
        vs: &[f64],
        scaled_step: f64,
        max_triangles: usize,
    ) -> Result<()> {
        if us.is_empty() && vs.is_empty() {
            return Ok(());
        }
        let inside = inside_faces(&self.cdt);
        let clearance = 0.3 * scaled_step;
        let spacing = 0.5 * scaled_step;
        // the four corners of a grid cell are cocircular, which sends every Delaunay
        // test of a regular grid down spade's exact-arithmetic path; a fixed pattern of
        // shifts of 1e-5 of the line spacing breaks the ties (and costs nothing in
        // accuracy: the points are evaluated where they end up)
        let b = self.bbox;
        let du = 1e-5 * (b.u_max() - b.u_min()) / (us.len() + 1) as f64;
        let dv = 1e-5 * (b.v_max() - b.v_min()) / (vs.len() + 1) as f64;
        let shift = |i: usize, j: usize, a: usize, c: usize| {
            ((i.wrapping_mul(a).wrapping_add(j.wrapping_mul(c))) % 7) as f64 / 7.0 - 0.5
        };
        let at = |i: usize, j: usize| {
            let (u, v) = (us.get(i).copied()?, vs.get(j).copied()?);
            Some(Point2::new(
                u + du * shift(i, j, 3, 5),
                v + dv * shift(i, j, 5, 3),
            ))
        };
        let mut wanted = Vec::new();
        for j in 0..vs.len() {
            for i in 0..us.len() {
                let Some(uv) = at(i, j) else {
                    continue;
                };
                let pos = self.scaled(uv);
                let PositionInTriangulation::OnFace(fh) = self.locate(pos) else {
                    continue;
                };
                if !is_inside(&inside, fh) || constraint_clearance(&self.cdt, fh, pos) < clearance {
                    continue;
                }
                let near = self.cdt.face(fh).vertices().iter().any(|w| {
                    let q = w.position();
                    ((q.x - pos.x).powi(2) + (q.y - pos.y).powi(2)).sqrt() < spacing
                });
                if !near {
                    wanted.push((i, j));
                }
            }
        }
        for (i, j) in wanted {
            // the face's triangle budget holds whatever the lines asked for
            if self.cdt.num_inner_faces() >= max_triangles {
                break;
            }
            let Some(uv) = at(i, j) else {
                continue;
            };
            self.insert_interior(k, s, uv, None, 0)?;
        }
        Ok(())
    }

    /// The sign of `(Su x Sv) . EvaluateNormal` at the face's origin, so a fast-path
    /// normal `Su x Sv` can be turned like the surface normal.
    fn cross_sign(&self, k: &Kernel, s: &Surface) -> Option<f64> {
        let uv = self.origin;
        let (_, su, sv) = s.evaluate_1st_derivatives(k, uv, false, false).ok()?;
        let n = surface_normal(k, s, uv)?;
        let d = su.cross(sv).dot(n);
        (d != 0.0 && d.is_finite()).then(|| d.signum())
    }

    /// Whether a triangle counter-clockwise in UV faces inward: the sign of
    /// `(Su x Sv) . outward` where `Su x Sv` is largest on a sample grid of `bbox`.
    fn uv_ccw_faces_inward(&self, k: &Kernel, s: &Surface, bbox: &Extent2d) -> Result<bool> {
        let mut best: Option<(f64, f64)> = None;
        for j in 0..=4 {
            for i in 0..=4 {
                let uv = Point2::new(
                    bbox.u_min() + (bbox.u_max() - bbox.u_min()) * (0.1 + 0.8 * f64::from(i) / 4.0),
                    bbox.v_min() + (bbox.v_max() - bbox.v_min()) * (0.1 + 0.8 * f64::from(j) / 4.0),
                );
                let Ok((_, su, sv)) = s.evaluate_1st_derivatives(k, uv, false, false) else {
                    continue;
                };
                let Some(n) = surface_normal(k, s, uv) else {
                    continue;
                };
                let c = su.cross(sv);
                let l = c.length();
                if l.is_finite() && best.is_none_or(|(bl, _)| l > bl) {
                    best = Some((l, c.dot(n) * self.sense));
                }
            }
        }
        match best {
            Some((l, d)) if l > 0.0 => Ok(d < 0.0),
            _ => Err(invalid("surface has no normal on the face")),
        }
    }

    /// Whether `uv` lies inside the loops, off every vertex and boundary edge, and at
    /// least `need` (scaled UV) from the boundary edges near it.
    fn clear_of_boundary(&self, inside: &[bool], uv: Point2, need: f64) -> bool {
        let pos = self.scaled(uv);
        let faces = match self.locate(pos) {
            PositionInTriangulation::OnFace(f) => vec![f],
            PositionInTriangulation::OnEdge(e)
                if !self.cdt.is_constraint_edge(e.as_undirected()) =>
            {
                let d = self.cdt.directed_edge(e);
                [d.face().as_inner(), d.rev().face().as_inner()]
                    .into_iter()
                    .flatten()
                    .map(|f| f.fix())
                    .collect()
            }
            _ => return false,
        };
        !faces.is_empty()
            && faces
                .into_iter()
                .all(|f| is_inside(inside, f) && constraint_clearance(&self.cdt, f, pos) >= need)
    }

    /// Asks for boundary segment `r.segment` to be split (into the most pieces asked).
    fn request(&mut self, r: SegmentSplit) {
        match self.requests.iter_mut().find(|q| q.segment == r.segment) {
            Some(q) => q.pieces = q.pieces.max(r.pieces),
            None => self.requests.push(r),
        }
    }

    /// The request to split boundary segment `e` (between two mesh vertices) into one
    /// piece more than the grid lines it crosses, at least two.
    fn split_request(&self, e: FixedUndirectedEdgeHandle) -> Option<SegmentSplit> {
        let [va, vb] = self.cdt.undirected_edge(e).vertices();
        let a = self.vertices.get(va.data().id)?;
        let b = self.vertices.get(vb.data().id)?;
        let (x, y) = (a.vertex?, b.vertex?);
        if x == y {
            return None;
        }
        let crossed = |lines: &[f64], p: f64, q: f64| {
            let (lo, hi) = if p < q { (p, q) } else { (q, p) };
            lines.iter().filter(|t| **t > lo && **t < hi).count()
        };
        let n = crossed(&self.lines.0, a.uv.x, b.uv.x)
            .max(crossed(&self.lines.1, a.uv.y, b.uv.y))
            .saturating_add(1);
        let pieces = u32::try_from(n)
            .unwrap_or(MAX_SPLIT_PIECES)
            .clamp(2, MAX_SPLIT_PIECES);
        Some(SegmentSplit {
            segment: if x < y { [x, y] } else { [y, x] },
            pieces,
        })
    }

    /// The split request for the boundary segment nearest `uv` (inside the loops) when
    /// it is longer than `len` (scaled UV).
    fn blocking_segment(&self, inside: &[bool], uv: Point2, len: f64) -> Option<SegmentSplit> {
        let pos = self.scaled(uv);
        let f = match self.locate(pos) {
            PositionInTriangulation::OnFace(f) => f,
            PositionInTriangulation::OnEdge(e) => {
                let d = self.cdt.directed_edge(e);
                [d.face().as_inner(), d.rev().face().as_inner()]
                    .into_iter()
                    .flatten()
                    .next()?
                    .fix()
            }
            _ => return None,
        };
        if !is_inside(inside, f) {
            return None;
        }
        let (_, e) = nearest_constraint(&self.cdt, f, pos);
        let e = e?;
        let [x, y] = self.cdt.undirected_edge(e).positions();
        let seg = ((y.x - x.x).powi(2) + (y.y - x.y).powi(2)).sqrt();
        if seg > len {
            self.split_request(e)
        } else {
            None
        }
    }

    /// A point a quarter of directed edge `e`'s length beyond its midpoint, on the side
    /// away from `e`'s face, when it falls inside the loops clear of the boundary.
    fn beyond_chord(
        &self,
        inside: &[bool],
        e: spade::handles::DirectedEdgeHandle<'_, CdtVertex, (), spade::CdtEdge<()>, ()>,
    ) -> Option<Point2> {
        let [p, q] = e.positions();
        let (dx, dy) = (q.x - p.x, q.y - p.y);
        let len = (dx * dx + dy * dy).sqrt();
        if len <= 0.0 || len.is_nan() {
            return None;
        }
        // the face is on the edge's left: go right
        let pos = spade::Point2::new((p.x + q.x) / 2.0 + dy * 0.25, (p.y + q.y) / 2.0 - dx * 0.25);
        let PositionInTriangulation::OnFace(f) = self.locate(pos) else {
            return None;
        };
        if !is_inside(inside, f) || constraint_clearance(&self.cdt, f, pos) < 0.05 * len {
            return None;
        }
        Some(self.unscaled(pos))
    }

    /// The UV points the next refinement pass inserts: the midpoint of every interior
    /// edge off the chord (or too long), and for every triangle off the surface,
    /// flipped or turned too far from it, the midpoint of its longest edge in scaled UV
    /// (its centroid when that edge is on the boundary). Every point keeps clear of the
    /// boundary, so no sliver along it is made.
    fn refine_candidates(
        &mut self,
        k: &Kernel,
        s: &Surface,
        p: &FaceParams,
        changed: Option<&[FixedVertexHandle]>,
    ) -> Result<Vec<(Point2, u8)>> {
        let inside = inside_faces(&self.cdt);
        // every triangle and edge, or (after a pass that only inserted points) those
        // around the new points: every other one was checked with the same neighbours
        let (faces, edges): (
            Vec<FixedFaceHandle<InnerTag>>,
            Vec<FixedUndirectedEdgeHandle>,
        ) = match changed {
            None => (
                self.cdt.inner_faces().map(|f| f.fix()).collect(),
                self.cdt.undirected_edges().map(|e| e.fix()).collect(),
            ),
            Some(vs) => {
                let mut fs = BTreeSet::new();
                let mut es = BTreeSet::new();
                for v in vs {
                    for e in self.cdt.vertex(*v).out_edges() {
                        es.insert(e.fix().as_undirected());
                        if let Some(f) = e.face().as_inner() {
                            fs.insert(f.fix());
                        }
                    }
                }
                (fs.into_iter().collect(), es.into_iter().collect())
            }
        };
        let tol = p.chord;
        let mut split_edges: Vec<FixedUndirectedEdgeHandle> = Vec::new();
        let mut marked: BTreeSet<usize> = BTreeSet::new();
        let mut centroids: Vec<(Point2, f64, Option<SegmentSplit>, u8)> = Vec::new();
        let mut newly_good = Vec::new();
        for fh in faces {
            if !is_inside(&inside, fh) {
                continue;
            }
            let face = self.cdt.face(fh);
            let [va, vb, vc] = face.vertices();
            let ids = {
                let mut x = [va.data().id, vb.data().id, vc.data().id];
                x.sort_unstable();
                x
            };
            if self.good.contains(&ids) {
                continue;
            }
            let a = self.vertex(va.fix())?;
            let b = self.vertex(vb.fix())?;
            let c = self.vertex(vc.fix())?;
            // a triangle collapsed onto a pole is dropped on output
            let same = |x: &FaceVertex, y: &FaceVertex| x.vertex.is_some() && x.vertex == y.vertex;
            if same(a, b) || same(b, c) || same(c, a) {
                continue;
            }
            let uv = if self.annulus.is_some() {
                self.plane_blend([a.uv, b.uv, c.uv], [1.0 / 3.0; 3])
            } else {
                Point2::new(
                    (a.uv.x + b.uv.x + c.uv.x) / 3.0,
                    (a.uv.y + b.uv.y + c.uv.y) / 3.0,
                )
            };
            let mut n = (b.xyz - a.xyz).cross(c.xyz - a.xyz);
            if self.flip {
                n = -n;
            }
            let nl = n.length();
            let size = a
                .xyz
                .distance(b.xyz)
                .max(b.xyz.distance(c.xyz))
                .max(c.xyz.distance(a.xyz));
            let mut bad = p.max_len.is_some_and(|l| size > l);
            if !bad && nl > 0.0 && nl.is_finite() {
                // the surface at the centroid and between it and each corner (where
                // a twisted patch leaves a flat triangle) against the triangle's plane
                for w in [
                    [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
                    [2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0],
                    [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0],
                    [1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0],
                ] {
                    let q = if self.annulus.is_some() {
                        self.plane_blend([a.uv, b.uv, c.uv], w)
                    } else {
                        Point2::new(
                            a.uv.x * w[0] + b.uv.x * w[1] + c.uv.x * w[2],
                            a.uv.y * w[0] + b.uv.y * w[1] + c.uv.y * w[2],
                        )
                    };
                    let pt = self.point(k, s, q)?;
                    if ((pt - a.xyz).dot(n) / nl).abs() > tol + a.gap.max(b.gap).max(c.gap) {
                        bad = true;
                        break;
                    }
                }
                // the surface folds over inside the triangle (its normal reverses
                // between corners): no triangle there can follow it, so no flip or
                // angle test
                let normals: Vec<Vec3> = [a, b, c].iter().filter_map(|x| x.normal).collect();
                let folded = normals
                    .iter()
                    .enumerate()
                    .any(|(i, x)| normals.iter().skip(i + 1).any(|y| x.dot(*y) < 0.0));
                if folded {
                    // ... and its size: a fold is resolved down to a band of small
                    // triangles along it
                    bad = bad || size > FOLD_SIZE * tol;
                } else if bad {
                } else if let Some(o) = self.outward(k, s, uv) {
                    // flipped, or turned from the surface by more than the angle
                    // tolerance at the centroid or an interior corner (a boundary
                    // corner's normal may be a singular point's)
                    let facet = n / nl;
                    let interior = [a, b, c]
                        .into_iter()
                        .filter(|x| !x.boundary)
                        .filter_map(|x| x.normal);
                    // a sliver's facet normal is no measure of anything (a fan at a
                    // pole is all slivers): only the flip counts there
                    let sliver = nl < SLIVER * size * size;
                    bad = o.dot(n) < 0.0
                        || (!sliver
                            && core::iter::once(o)
                                .chain(interior)
                                .any(|q| angle(facet, q) > p.angle));
                }
            }
            // a triangle far below the tolerance or the length floor is left alone
            if !bad || size < MIN_SPLIT_FRACTION * tol || size < 2.0 * p.min_len {
                if !bad {
                    newly_good.push(ids);
                }
                continue;
            }
            let mut longest = None;
            let mut longest_len = -1.0;
            for e in face.adjacent_edges() {
                let [x, y] = e.positions();
                let l = (y.x - x.x).powi(2) + (y.y - x.y).powi(2);
                if l > longest_len {
                    longest_len = l;
                    longest = Some(e);
                }
            }
            let Some(e) = longest else {
                continue;
            };
            let ue = e.fix().as_undirected();
            let depth = a.depth.max(b.depth).max(c.depth).saturating_add(1);
            if self.cdt.is_constraint_edge(ue) {
                // too close to that boundary segment, the segment itself is split
                let request = self.split_request(ue);
                centroids.push((uv, CLEARANCE_FRACTION * longest_len.sqrt(), request, depth));
            } else if [a, b, c].iter().all(|x| x.boundary) {
                // all corners on the boundary (an ear between boundary points):
                // splitting the chord would put a point next to the boundary, so a
                // point goes beyond it, into the face across it
                if let Some(q) = self.beyond_chord(&inside, e) {
                    centroids.push((q, 0.0, None, depth));
                }
            } else if marked.insert(ue.index()) {
                split_edges.push(ue);
            }
        }
        let mut good_edges = Vec::new();
        for eh in edges {
            let edge = self.cdt.undirected_edge(eh);
            if self.cdt.is_constraint_edge(eh) || marked.contains(&eh.index()) {
                continue;
            }
            let [va, vb] = edge.vertices();
            let key = {
                let (x, y) = (va.data().id, vb.data().id);
                if x < y {
                    [x, y]
                } else {
                    [y, x]
                }
            };
            if self.good_edges.contains(&key) {
                continue;
            }
            let d = edge.as_directed();
            let in_face =
                |f: Option<FixedFaceHandle<InnerTag>>| f.is_some_and(|f| is_inside(&inside, f));
            if !(in_face(d.face().as_inner().map(|f| f.fix()))
                || in_face(d.rev().face().as_inner().map(|f| f.fix())))
            {
                continue;
            }
            let a = self.vertex(va.fix())?;
            let b = self.vertex(vb.fix())?;
            if a.vertex.is_some() && a.vertex == b.vertex {
                continue;
            }
            let len = a.xyz.distance(b.xyz);
            if len < 2.0 * p.min_len {
                continue;
            }
            let uv = if self.annulus.is_some() {
                self.plane_blend([a.uv, b.uv, b.uv], [0.5, 0.5, 0.0])
            } else {
                Point2::new((a.uv.x + b.uv.x) / 2.0, (a.uv.y + b.uv.y) / 2.0)
            };
            let mut split = p.max_len.is_some_and(|l| len > l);
            if !split {
                let pt = self.point(k, s, uv)?;
                split = segment_distance(pt, a.xyz, b.xyz) > tol + a.gap.max(b.gap);
            }
            if split && marked.insert(edge.fix().index()) {
                split_edges.push(edge.fix());
            } else if !split {
                good_edges.push(key);
            }
        }
        self.good.extend(newly_good);
        self.good_edges.extend(good_edges);
        let mut out = Vec::with_capacity(split_edges.len() + centroids.len());
        for e in split_edges {
            let edge = self.cdt.undirected_edge(e);
            let [x, y] = edge.positions();
            let len = ((y.x - x.x).powi(2) + (y.y - x.y).powi(2)).sqrt();
            let [va, vb] = edge.vertices();
            let a = self.vertex(va.fix())?;
            let b = self.vertex(vb.fix())?;
            let uv = if self.annulus.is_some() {
                self.plane_blend([a.uv, b.uv, b.uv], [0.5, 0.5, 0.0])
            } else {
                Point2::new((a.uv.x + b.uv.x) / 2.0, (a.uv.y + b.uv.y) / 2.0)
            };
            let depth = a.depth.max(b.depth).saturating_add(1);
            if depth > MAX_LOCAL_DEPTH {
                self.depth_capped = true;
            } else if self.clear_of_boundary(&inside, uv, CLEARANCE_FRACTION * len) {
                out.push((uv, depth));
            } else if let Some(r) = self.blocking_segment(&inside, uv, len) {
                // a boundary segment longer than the edge keeps the split off it (a fan
                // from a long segment's ends): that segment is split
                self.request(r);
            }
        }
        for (uv, need, request, depth) in centroids {
            if depth > MAX_LOCAL_DEPTH {
                self.depth_capped = true;
            } else if self.clear_of_boundary(&inside, uv, need) {
                out.push((uv, depth));
            } else if let Some(r) = request {
                self.request(r);
            }
        }
        Ok(out)
    }

    /// Fills a face with no triangle in UV from its one loop, by clipping ears in the
    /// loop's best-fit plane (no point added, so its neighbours still meet it).
    fn ear_clip(&mut self, k: &Kernel, s: &Surface) -> Result<()> {
        let [lp] = self.boundary.as_slice() else {
            return Ok(());
        };
        if lp.len() < 3 || lp.len() > MAX_EAR_LOOP {
            return Ok(());
        }
        let mut ids = Vec::with_capacity(lp.len());
        for (uv, v) in lp {
            let Some(h) = self.find_boundary(*uv, *v)? else {
                return Ok(());
            };
            let id = self.id(h);
            if ids.last() != Some(&id) {
                ids.push(id);
            }
        }
        let pts: Vec<Point3> = ids
            .iter()
            .filter_map(|i| self.vertices.get(*i).map(|v| v.xyz))
            .collect();
        // Newell normal of the loop
        let mut nrm = Vec3::new(0.0, 0.0, 0.0);
        for (i, a) in pts.iter().enumerate() {
            let Some(b) = pts.get((i + 1) % pts.len()) else {
                continue;
            };
            nrm += Vec3::new(
                (a.y - b.y) * (a.z + b.z),
                (a.z - b.z) * (a.x + b.x),
                (a.x - b.x) * (a.y + b.y),
            );
        }
        let Ok(nrm) = nrm.unitize() else {
            return Ok(());
        };
        let tris = ear_clip_polygon(&pts, nrm);
        if tris.is_empty() {
            return Ok(());
        }
        // wind the triangles outward
        let outward = ids
            .first()
            .and_then(|i| self.vertices.get(*i))
            .and_then(|v| self.outward(k, s, v.uv));
        let reverse = outward.is_some_and(|o| o.dot(nrm) < 0.0);
        self.ears = tris
            .into_iter()
            .filter_map(|[a, b, c]| {
                let (a, b, c) = (*ids.get(a)?, *ids.get(b)?, *ids.get(c)?);
                Some(if reverse { [a, c, b] } else { [a, b, c] })
            })
            .collect();
        Ok(())
    }

    /// The face's current triangles over its vertex table, wound outward (corners on
    /// one pole vertex dropped).
    fn triangles(&self) -> Vec<[usize; 3]> {
        if !self.ears.is_empty() {
            return self.ears.clone();
        }
        let inside = inside_faces(&self.cdt);
        let mut out = Vec::new();
        for face in self.cdt.inner_faces() {
            if !is_inside(&inside, face.fix()) {
                continue;
            }
            let [a, b, c] = face.vertices().map(|v| v.data().id);
            let shared = |i: usize| self.vertices.get(i).and_then(|v| v.vertex);
            let same = |x: usize, y: usize| shared(x).is_some() && shared(x) == shared(y);
            if same(a, b) || same(b, c) || same(c, a) {
                continue;
            }
            out.push(if self.flip { [a, c, b] } else { [a, b, c] });
        }
        out
    }

    /// The face's mesh now, over the vertices its triangles use.
    pub(crate) fn snapshot(
        &self,
        k: &Kernel,
        b: &Brep,
        chord: Option<f64>,
        backend: Option<&dyn FaceBackend>,
    ) -> Result<FaceMesh> {
        let tris = self.triangles();
        let mut local: BTreeMap<usize, u32> = BTreeMap::new();
        let mut m = FaceMesh {
            face: self.face,
            inverted_edgeuses: self.inverted.clone(),
            unconverged: self.unconverged,
            split_boundary: self.split_boundary,
            ear_clipped: !self.ears.is_empty(),
            uvs: Vec::new(),
            positions: Vec::new(),
            normals: Vec::new(),
            shared: Vec::new(),
            triangles: Vec::new(),
            gaps: Vec::new(),
            deviation: None,
            triangles_over: 0,
            backend: self.backend,
            backend_failed: self.backend_failed,
            backend_deviation: None,
        };
        for t in &tris {
            let mut corners = [0u32; 3];
            for (slot, id) in corners.iter_mut().zip(t) {
                if let Some(l) = local.get(id) {
                    *slot = *l;
                    continue;
                }
                let v = self
                    .vertices
                    .get(*id)
                    .ok_or_else(|| invalid("vertex out of range"))?;
                let l = u32::try_from(m.uvs.len())
                    .map_err(|_| invalid("face has too many vertices"))?;
                m.uvs.push(v.uv);
                m.positions.push(v.xyz);
                m.normals.push(v.normal.unwrap_or(Vec3::new(0.0, 0.0, 0.0)));
                m.shared.push(v.vertex);
                m.gaps.push(v.gap);
                local.insert(*id, l);
                *slot = l;
            }
            m.triangles.push(corners);
        }
        if let Some(chord) = chord {
            let s = surface_of(b, self.face)?;
            let (search, _) = self.map_domains();
            let (dev, over) =
                crate::deviation::face_mesh_deviation(k, s, &search, self.period, &m, chord);
            m.deviation = Some(dev);
            m.triangles_over = over;
            // the backend's own measure of a face whose grid it placed
            if let (Some(be), Some(_)) = (backend, self.backend) {
                let face = BackendFace {
                    face: self.face,
                    surface: s,
                    bbox: self.bbox,
                    chord,
                    angle: 0.0,
                };
                let view = FaceMeshView {
                    uvs: &m.uvs,
                    positions: &m.positions,
                    triangles: &m.triangles,
                };
                m.backend_deviation = be.deviation(k, &face, &view);
            }
        }
        Ok(m)
    }

    /// The face.
    pub(crate) fn face(&self) -> FaceId {
        self.face
    }

    /// Boundary segments (sorted mesh vertex pairs, with the pieces each should become)
    /// that the last refinement could not resolve from inside the face: splitting them in the edge polylines
    /// ([`EdgeSet::split_segments`]) and refining again lets it.
    pub(crate) fn requests(&self) -> &[SegmentSplit] {
        &self.requests
    }

    /// Whether edge `e` bounds the face.
    pub(crate) fn uses_edge(&self, e: solidcore::EdgeId) -> bool {
        self.edges_used.contains(&e)
    }

    /// Records the shared mesh vertices an assembly gave the snapshot's interior
    /// vertices (matched by UV), so later levels keep them.
    pub(crate) fn adopt(&mut self, m: &FaceMesh) {
        let mut by_uv: BTreeMap<(u64, u64), u32> = BTreeMap::new();
        for (uv, sh) in m.uvs.iter().zip(&m.shared) {
            if let Some(v) = sh {
                by_uv.insert((uv.x.to_bits(), uv.y.to_bits()), *v);
            }
        }
        for v in &mut self.vertices {
            if v.vertex.is_none() {
                if let Some(id) = by_uv.get(&(v.uv.x.to_bits(), v.uv.y.to_bits())) {
                    v.vertex = Some(*id);
                }
            }
        }
    }
}

/// What mapping an edge use into UV needs.
struct MapContext<'a> {
    k: &'a Kernel,
    s: &'a Surface,
    /// Where Newton feet may go: the natural domain, a period wider on each side in a
    /// closed direction (the loop is unwrapped across it, not stopped at it).
    domain: &'a Extent2d,
    /// Where a point with no seed is first looked for ([`invert_point`]): the natural
    /// domain where bounded.
    grid: &'a Extent2d,
    edges: &'a EdgeSet,
    polish_above: f64,
    period: (Option<f64>, Option<f64>),
}

impl MapContext<'_> {
    /// `uv` moved by whole periods to within half a period of `near`.
    fn unwrap(&self, uv: Point2, near: Point2) -> Point2 {
        let shift = |x: f64, c: f64, p: Option<f64>| match p {
            Some(p) if p > 0.0 => x - ((x - c) / p).round() * p,
            _ => x,
        };
        Point2::new(
            shift(uv.x, near.x, self.period.0),
            shift(uv.y, near.y, self.period.1),
        )
    }

    /// The points of one edge use in loop order (`reverse`: against the edge) mapped
    /// into UV: by the UV trim curve at the same normalised edge parameter when given
    /// (polished towards the edge point only when its surface point is off it), else
    /// by point inversion continuing from the previous point (`seed` for the first),
    /// kept within half a period of it on a closed surface. With `end` (where the edge
    /// use is known to end) and `seed`, each inversion starts from the point that far
    /// along the line between them instead: `end` as given when it is exact (a trim
    /// curve's point, `true`), else moved by whole periods to the seed's side.
    fn edgeuse_points(
        &self,
        poly: &EdgePolyline,
        ivl: &Extent1d,
        trim: Option<&Curve>,
        reverse: bool,
        seed: Option<Point2>,
        end: Option<(Point2, bool)>,
    ) -> Result<Vec<LoopPoint>> {
        // the end: as given when exact, else a whole number of periods moved to the
        // side of the start
        // (not for a closed edge: its end is its start, a whole period on or not, so
        // only continuity from the start finds which way round it goes)
        let closed = poly.vertices.len() > 1 && poly.vertices.first() == poly.vertices.last();
        let towards = match (seed, end) {
            (Some(a), Some((z, exact))) if !closed => {
                Some((a, if exact { z } else { self.unwrap(z, a) }))
            }
            _ => None,
        };
        let (k, s) = (self.k, self.s);
        let (t0, t1) = (ivl.min(), ivl.max());
        let len = t1 - t0;
        let n = poly.params.len();
        let mut out = Vec::with_capacity(n);
        let mut prev = seed;
        // a UV trim curve must follow its edge: at its ends (an imported curve may run
        // against the edge, or over another interval than the edge's) and half way
        // (a closed edge's curve run backwards meets the right ends)
        let check = if trim.is_some() && !poly.degenerate {
            let mut length = 0.0;
            for w in poly.vertices.windows(2) {
                if let [a, c] = w {
                    length += self.edges.position(*a)?.distance(self.edges.position(*c)?);
                }
            }
            Some((self.polish_above.max(0.05 * length), 0.25 * length))
        } else {
            None
        };
        for pos in 0..n {
            let i = if reverse { n - 1 - pos } else { pos };
            let (Some(t), Some(v)) = (poly.params.get(i), poly.vertices.get(i)) else {
                continue;
            };
            let xyz = self.edges.position(*v)?;
            let sn = if len > 0.0 {
                ((t - t0) / len).clamp(0.0, 1.0)
            } else if i == 0 {
                0.0
            } else {
                1.0
            };
            let uv = if let Some(c) = trim {
                let ti = c.natural_interval();
                let q = c.evaluate_point(k, ti.evaluate(sn))?;
                if !(q.x.is_finite() && q.y.is_finite()) {
                    return Err(invalid("UV trim curve is not finite"));
                }
                let uv = Point2::new(q.x, q.y);
                let off = surface_point(k, s, uv)?.distance(xyz);
                if let Some((at_end, inside)) = check {
                    let limit = if pos == 0 || pos + 1 == n {
                        at_end
                    } else if pos == n / 2 {
                        inside
                    } else {
                        f64::INFINITY
                    };
                    if off > limit {
                        return Err(invalid("UV trim curve does not follow its edge"));
                    }
                }
                let uv = if poly.degenerate || off <= self.polish_above {
                    uv
                } else {
                    foot(k, s, uv, xyz, self.domain).0
                };
                // continuous across a seam: a curve the kernel keeps within the domain
                // may jump a period between two points
                prev.map_or(uv, |q| self.unwrap(uv, q))
            } else {
                // along the use: the edge parameter's fraction in loop order
                let along = if reverse { 1.0 - sn } else { sn };
                let guess = towards.map(|(a, z)| {
                    Point2::new(a.x + (z.x - a.x) * along, a.y + (z.y - a.y) * along)
                });
                // a foot from the guess and one from the previous point: the guess's
                // unless the other is clearly closer (an edge use the long way round)
                let polish = |seed: Point2| {
                    let (uv, d) = foot(k, s, seed, xyz, self.domain);
                    (self.unwrap(uv, seed), d.unwrap_or(f64::INFINITY))
                };
                match (guess.map(polish), prev.map(polish)) {
                    (Some((g, dg)), Some((q, dq))) => {
                        if dg <= dq + self.polish_above {
                            g
                        } else {
                            q
                        }
                    }
                    (Some((uv, _)), None) | (None, Some((uv, _))) => uv,
                    (None, None) => invert_point(k, s, xyz, self.grid, self.domain)?,
                }
            };
            prev = Some(uv);
            out.push(LoopPoint {
                uv,
                xyz,
                vertex: *v,
                pole_next: poly.degenerate && pos == 0,
            });
        }
        Ok(out)
    }
}

/// `su x sv` as a unit vector, `None` where it is too small against the tangents to
/// give a direction (a pole).
fn unit_cross(su: Vec3, sv: Vec3) -> Option<Vec3> {
    let c = su.cross(sv);
    let l = c.length();
    let m = su.length_squared().max(sv.length_squared());
    (l.is_finite() && l > 1e-10 * m).then(|| c / l)
}

/// How far `xyz` is from the surface point at `uv` (`0` when it does not evaluate).
fn surface_gap(k: &Kernel, s: &Surface, uv: Point2, xyz: Point3) -> f64 {
    s.evaluate_point(k, uv).map_or(0.0, |p| p.distance(xyz))
}

/// Distance in scaled UV between `a` and `c`.
fn scaled_gap(a: Point2, c: Point2, scale: (f64, f64)) -> f64 {
    let (du, dv) = ((a.x - c.x) * scale.0, (a.y - c.y) * scale.1);
    (du * du + dv * dv).sqrt()
}

/// How the UV jump from `a` to `b` (one mesh vertex) closes: a pole when they are far
/// apart in UV but the surface stays at the vertex between them; a whole number of
/// periods on a closed surface (`period`) when moving `b` by that brings it to `a`. Any
/// other jump across most of the domain is a loop that does not close.
fn junction(
    k: &Kernel,
    s: &Surface,
    chord: f64,
    a: LoopPoint,
    b: LoopPoint,
    span: f64,
    scale: (f64, f64),
    period: (Option<f64>, Option<f64>),
) -> Result<Junction> {
    let close = if span.is_finite() {
        chord.max(1e-9 * span)
    } else {
        chord
    };
    let gap = scaled_gap(a.uv, b.uv, scale);
    if gap <= close {
        return Ok(Junction::Same);
    }
    let mid = Point2::new((a.uv.x + b.uv.x) / 2.0, (a.uv.y + b.uv.y) / 2.0);
    if surface_point(k, s, mid)?.distance(a.xyz) <= chord {
        return Ok(Junction::Pole);
    }
    let near = |x: f64, c: f64, p: Option<f64>| match p {
        Some(p) if p > 0.0 => x - ((x - c) / p).round() * p,
        _ => x,
    };
    let w = Point2::new(
        near(b.uv.x, a.uv.x, period.0),
        near(b.uv.y, a.uv.y, period.1),
    );
    if !same_uv(w, b.uv) && scaled_gap(a.uv, w, scale) <= close {
        return Ok(Junction::Period(w.x - b.uv.x, w.y - b.uv.y));
    }
    if gap > 0.1 * span {
        return Err(invalid(
            "loop does not close in UV (periodic face without seam edge)",
        ));
    }
    Ok(Junction::Same)
}

/// Appends an edge use's points to a loop, merging the shared junction point (the edge
/// use moved by whole periods where it starts a period away). Returns that move.
fn join(
    pts: &mut Vec<LoopPoint>,
    mut seq: Vec<LoopPoint>,
    k: &Kernel,
    s: &Surface,
    chord: f64,
    span: f64,
    scale: (f64, f64),
    period: (Option<f64>, Option<f64>),
) -> Result<(f64, f64)> {
    let mut moved = (0.0, 0.0);
    if let (Some(last), Some(first)) = (pts.last().copied(), seq.first().copied()) {
        if last.vertex != first.vertex {
            return Err(invalid("loop edge uses do not connect"));
        }
        match junction(k, s, chord, last, first, span, scale, period)? {
            Junction::Pole => {
                if let Some(l) = pts.last_mut() {
                    l.pole_next = true;
                }
            }
            j => {
                if let Junction::Period(du, dv) = j {
                    for q in &mut seq {
                        q.uv = Point2::new(q.uv.x + du, q.uv.y + dv);
                    }
                    moved = (du, dv);
                }
                pts.pop();
                if let Some(f) = seq.first_mut() {
                    f.pole_next |= last.pole_next;
                }
            }
        }
    }
    pts.extend(seq);
    Ok(moved)
}

/// Flood fill from the outer face, flipping at every constraint edge: `true` for faces
/// inside the loops, indexed by face handle.
fn inside_faces(cdt: &Cdt) -> Vec<bool> {
    let n = cdt.num_all_faces();
    let mut label: Vec<Option<bool>> = vec![None; n];
    let outer = cdt.outer_face().fix();
    if let Some(l) = label.get_mut(outer.index()) {
        *l = Some(false);
    }
    let mut queue = VecDeque::from([outer]);
    let limit = cdt.num_directed_edges();
    while let Some(fh) = queue.pop_front() {
        let inside = label.get(fh.index()).copied().flatten().unwrap_or(false);
        let Some(start) = cdt.face(fh).adjacent_edge() else {
            continue;
        };
        let mut e = start;
        for _ in 0..limit {
            let across = e.rev().face().fix();
            let constraint = cdt.is_constraint_edge(e.fix().as_undirected());
            if let Some(slot) = label.get_mut(across.index()) {
                if slot.is_none() {
                    *slot = Some(inside ^ constraint);
                    queue.push_back(across);
                }
            }
            e = e.next();
            if e == start {
                break;
            }
        }
    }
    label.into_iter().map(|l| l.unwrap_or(false)).collect()
}

fn is_inside(inside: &[bool], f: FixedFaceHandle<InnerTag>) -> bool {
    inside.get(f.index()).copied().unwrap_or(false)
}

/// Whether segment `a`-`c` crosses segment `p`-`q` at an angle whose sine is at least
/// [`MIN_CROSSING_SINE`] and at least [`CROSSING_MARGIN`] of each from its ends: a
/// crossing whose rounded point stays on the right side of both.
fn clean_crossing(
    a: spade::Point2<f64>,
    c: spade::Point2<f64>,
    p: spade::Point2<f64>,
    q: spade::Point2<f64>,
) -> bool {
    let cross = |x: (f64, f64), y: (f64, f64)| x.0 * y.1 - x.1 * y.0;
    let d1 = (c.x - a.x, c.y - a.y);
    let d2 = (q.x - p.x, q.y - p.y);
    let l1 = (d1.0 * d1.0 + d1.1 * d1.1).sqrt();
    let l2 = (d2.0 * d2.0 + d2.1 * d2.1).sqrt();
    let den = cross(d1, d2);
    if !den.is_finite() || den.abs() < MIN_CROSSING_SINE * l1 * l2 {
        return false;
    }
    let ap = (p.x - a.x, p.y - a.y);
    let t = cross(ap, d2) / den;
    let u = cross(ap, d1) / den;
    let inner = |x: f64| (CROSSING_MARGIN..=1.0 - CROSSING_MARGIN).contains(&x);
    inner(t) && inner(u)
}

/// Distance in the plane from `p` to segment `a`-`b`.
fn segment_distance_2d(p: spade::Point2<f64>, a: spade::Point2<f64>, b: spade::Point2<f64>) -> f64 {
    let (dx, dy) = (b.x - a.x, b.y - a.y);
    let l2 = dx * dx + dy * dy;
    let s = if l2 > 0.0 {
        (((p.x - a.x) * dx + (p.y - a.y) * dy) / l2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let (qx, qy) = (a.x + dx * s - p.x, a.y + dy * s - p.y);
    (qx * qx + qy * qy).sqrt()
}

/// The distance from `p` (in face `f`) to the nearest constraint edge of `f` or of its
/// neighbours.
fn constraint_clearance(cdt: &Cdt, f: FixedFaceHandle<InnerTag>, p: spade::Point2<f64>) -> f64 {
    nearest_constraint(cdt, f, p).0
}

/// The nearest constraint edge to `p` of face `f` and its neighbours, and its distance
/// (infinite and `None` when there is none).
fn nearest_constraint(
    cdt: &Cdt,
    f: FixedFaceHandle<InnerTag>,
    p: spade::Point2<f64>,
) -> (f64, Option<FixedUndirectedEdgeHandle>) {
    let mut best = (f64::INFINITY, None);
    for e in cdt.face(f).adjacent_edges() {
        let mut check =
            |x: spade::handles::DirectedEdgeHandle<'_, CdtVertex, (), spade::CdtEdge<()>, ()>| {
                let u = x.fix().as_undirected();
                if cdt.is_constraint_edge(u) {
                    let [a, b] = x.positions();
                    let d = segment_distance_2d(p, a, b);
                    if d < best.0 {
                        best = (d, Some(u));
                    }
                }
            };
        check(e);
        if let Some(other) = e.rev().face().as_inner() {
            for x in other.adjacent_edges() {
                check(x);
            }
        }
    }
    best
}

/// Ear clipping of the polygon `pts` (no repeated closing point) projected on the plane
/// of unit normal `n`: triangles counter-clockwise about `n`, empty when it fails.
fn ear_clip_polygon(pts: &[Point3], n: Vec3) -> Vec<[usize; 3]> {
    // a frame of the plane
    // (`n` itself is the first axis it returns)
    let Ok([_, x, y]) = n.make_unit_ortho_vectors(None) else {
        return Vec::new();
    };
    let (x, y) = if x.cross(y).dot(n) < 0.0 {
        (y, x)
    } else {
        (x, y)
    };
    let p2: Vec<(f64, f64)> = pts.iter().map(|p| (p.dot(x), p.dot(y))).collect();
    let at = |i: usize| p2.get(i).copied().unwrap_or((0.0, 0.0));
    let cross = |a: (f64, f64), b: (f64, f64), c: (f64, f64)| {
        (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0)
    };
    let mut idx: Vec<usize> = (0..pts.len()).collect();
    let mut out = Vec::new();
    let mut guard = 0;
    while idx.len() > 3 && guard < pts.len() * pts.len() {
        guard += 1;
        let m = idx.len();
        let mut clipped = false;
        for i in 0..m {
            let (Some(&ia), Some(&ib), Some(&ic)) =
                (idx.get((i + m - 1) % m), idx.get(i), idx.get((i + 1) % m))
            else {
                continue;
            };
            let (a, b, c) = (at(ia), at(ib), at(ic));
            if cross(a, b, c) <= 0.0 {
                continue;
            }
            let blocked = idx.iter().any(|&j| {
                j != ia && j != ib && j != ic && {
                    let q = at(j);
                    cross(a, b, q) >= 0.0 && cross(b, c, q) >= 0.0 && cross(c, a, q) >= 0.0
                }
            });
            if blocked {
                continue;
            }
            out.push([ia, ib, ic]);
            idx.remove(i);
            clipped = true;
            break;
        }
        if !clipped {
            return Vec::new();
        }
    }
    if let [a, b, c] = idx.as_slice() {
        if cross(at(*a), at(*b), at(*c)) > 0.0 {
            out.push([*a, *b, *c]);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn area(pts: &[Point3], tris: &[[usize; 3]]) -> f64 {
        tris.iter()
            .map(|t| {
                (pts[t[1]] - pts[t[0]])
                    .cross(pts[t[2]] - pts[t[0]])
                    .length()
                    / 2.0
            })
            .sum()
    }

    #[test]
    fn ear_clipping_covers_the_polygon() {
        let z = Vec3::new(0.0, 0.0, 1.0);
        // an L of area 3, counter-clockwise about +z
        let l = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(2.0, 1.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(1.0, 2.0, 0.0),
            Point3::new(0.0, 2.0, 0.0),
        ];
        let t = ear_clip_polygon(&l, z);
        assert_eq!(t.len(), 4);
        assert!((area(&l, &t) - 3.0).abs() < 1e-12);
        // every triangle counter-clockwise about the normal
        for [a, b, c] in &t {
            assert!((l[*b] - l[*a]).cross(l[*c] - l[*a]).dot(z) > 0.0);
        }
        // seen from the other side the polygon is clockwise: no ear, nothing clipped
        assert!(ear_clip_polygon(&l, -z).is_empty());
    }

    #[test]
    fn clean_crossing_rejects_grazing_and_end_crossings() {
        let p = |x: f64, y: f64| spade::Point2::new(x, y);
        // a plain X
        assert!(clean_crossing(
            p(0.0, 0.0),
            p(2.0, 2.0),
            p(0.0, 2.0),
            p(2.0, 0.0)
        ));
        // nearly parallel: rejected
        assert!(!clean_crossing(
            p(0.0, 0.0),
            p(2.0, 0.0),
            p(0.0, 1e-5),
            p(2.0, -1e-5)
        ));
        // through an end point: rejected
        assert!(!clean_crossing(
            p(0.0, 0.0),
            p(2.0, 0.0),
            p(1.0, 0.0),
            p(1.0, 1.0)
        ));
        // apart: rejected
        assert!(!clean_crossing(
            p(0.0, 0.0),
            p(1.0, 0.0),
            p(2.0, -1.0),
            p(2.0, 1.0)
        ));
    }

    /// A split constraint, then many locates and inserts: the hint generator must
    /// follow the split vertices (spade's hierarchy did not and indexed out of
    /// bounds).
    #[test]
    fn split_constraints_keep_point_location_valid() {
        let mut cdt = Cdt::new();
        let ins = |x: f64, y: f64, id: usize, cdt: &mut Cdt| {
            cdt.insert(CdtVertex {
                pos: spade::Point2::new(x, y),
                id,
            })
            .unwrap()
        };
        let a = ins(0.0, 0.0, 0, &mut cdt);
        let b = ins(4.0, 4.0, 1, &mut cdt);
        let c = ins(0.0, 4.0, 2, &mut cdt);
        let d = ins(4.0, 0.0, 3, &mut cdt);
        assert!(!cdt.try_add_constraint(a, b).is_empty());
        let split = cdt.add_constraint_and_split(c, d, |pos| CdtVertex {
            pos,
            id: usize::MAX,
        });
        assert!(!split.is_empty());
        for i in 0..64 {
            let x = 0.1 + 3.8 * f64::from(i % 8) / 7.0;
            let y = 0.1 + 3.8 * f64::from(i / 8) / 7.0;
            let _ = cdt.locate(spade::Point2::new(x, y));
            ins(x + 1e-3, y + 2e-3, 10 + i as usize, &mut cdt);
        }
        assert_eq!(cdt.num_vertices(), 5 + 64);
    }

    fn cylinder_annulus() -> Annulus {
        Annulus {
            along_u: true,
            period: core::f64::consts::TAU,
            ref_a: core::f64::consts::PI,
            ref_b: 0.0,
            span_b: 1.0,
        }
    }

    #[test]
    fn annulus_map_round_trips_and_keeps_orientation() {
        for an in [
            cylinder_annulus(),
            Annulus {
                along_u: false,
                ..cylinder_annulus()
            },
        ] {
            let scale = (2.0, 0.5);
            for (u, v) in [(0.3, 0.2), (3.0, 0.9), (6.0, 0.5), (1.0, 0.0)] {
                let uv = if an.along_u {
                    Point2::new(u, v)
                } else {
                    Point2::new(v, u)
                };
                let (x, y) = an.plane_point(uv, scale);
                let back = an.plane_uv(x, y, scale);
                let (a0, a1) = if an.along_u {
                    (uv.x, back.x)
                } else {
                    (uv.y, back.y)
                };
                let turns = (a1 - a0) / an.period;
                assert!((turns - turns.round()).abs() < 1e-12, "{uv:?} {back:?}");
                let (b0, b1) = if an.along_u {
                    (uv.y, back.y)
                } else {
                    (uv.x, back.x)
                };
                assert!((b0 - b1).abs() < 1e-12, "{uv:?} {back:?}");
            }
            // a small triangle counter-clockwise in UV stays counter-clockwise
            let [p, q, r] = [(1.0, 0.5), (1.1, 0.5), (1.0, 0.6)]
                .map(|(x, y)| an.plane_point(Point2::new(x, y), scale));
            let area = (q.0 - p.0) * (r.1 - p.1) - (q.1 - p.1) * (r.0 - p.0);
            assert!(area > 0.0, "{:?}", an.along_u);
        }
    }

    /// A cylinder face bounded by two circles and no seam edge: both loops go once
    /// around u. In the annulus plane they are two concentric rings and the inside is
    /// the band between them, every triangle of it within the face's v range.
    #[test]
    fn seamless_cylinder_band_is_an_annulus() {
        let an = cylinder_annulus();
        let scale = (1.0, 1.0);
        let n = 32;
        let mut cdt = Cdt::new();
        let mut rings = Vec::new();
        for (ring, v) in [0.0, 1.0].into_iter().enumerate() {
            let mut hs = Vec::new();
            for i in 0..n {
                let u = core::f64::consts::TAU * f64::from(i) / f64::from(n);
                let (x, y) = an.plane_point(Point2::new(u, v), scale);
                hs.push(
                    cdt.insert(CdtVertex {
                        pos: spade::Point2::new(x, y),
                        id: ring * 100 + i as usize,
                    })
                    .unwrap(),
                );
            }
            rings.push(hs);
        }
        for hs in &rings {
            for i in 0..hs.len() {
                assert!(!cdt
                    .try_add_constraint(hs[i], hs[(i + 1) % hs.len()])
                    .is_empty());
            }
        }
        let inside = inside_faces(&cdt);
        let mut count = 0;
        for f in cdt.inner_faces() {
            if !is_inside(&inside, f.fix()) {
                continue;
            }
            count += 1;
            let [p, q, r] = f.positions();
            let c = an.plane_uv((p.x + q.x + r.x) / 3.0, (p.y + q.y + r.y) / 3.0, scale);
            assert!(c.y > 0.0 && c.y < 1.0, "{c:?}");
        }
        // a band of 32 + 32 boundary points: 64 triangles
        assert_eq!(count, 2 * n as usize);
    }

    /// A face of a unit box, as a [`FaceState`] with its UV box set (the inputs below
    /// are put into its triangulation directly). `None` where the box cannot be built.
    fn box_face() -> Option<(Kernel, Brep, FaceId)> {
        let mut k = Kernel::default();
        let b = solidcore::primitive::api_primitives::create_box(
            &mut k,
            Point3::new(0.0, 0.0, 0.0),
            1.0,
            1.0,
            1.0,
        )
        .ok()?;
        let f = b.faces().ids().next()?;
        Some((k, b, f))
    }

    fn boundary_point(uv: Point2, vertex: u32) -> FaceVertex {
        FaceVertex {
            uv,
            xyz: Point3::new(uv.x, uv.y, 0.0),
            vertex: Some(vertex),
            normal: None,
            boundary: true,
            alive: true,
            gap: 0.0,
            depth: 0,
        }
    }

    /// Grazing boundary crossings and near-duplicate points fail the face (or are
    /// skipped) instead of reaching spade's rounded splits, which panicked.
    #[test]
    fn grazing_and_near_duplicate_inputs_do_not_panic() {
        let Some((k, b, f)) = box_face() else {
            return;
        };
        let s = surface_of(&b, f).unwrap();
        let mut st = FaceState::new(&k, &b, f, 1.0).unwrap();
        st.bbox = Extent2d::from_bounds(0.0, 0.0, 1.0, 1.0);
        st.origin = Point2::new(0.0, 0.0);
        st.scale = (1.0, 1.0);
        // a segment, then one crossing it at a grazing angle
        let a = st.insert(boundary_point(Point2::new(0.0, 0.5), 0)).unwrap();
        let c = st.insert(boundary_point(Point2::new(1.0, 0.5), 1)).unwrap();
        st.constrain(&k, s, a, c).unwrap();
        let p = st
            .insert(boundary_point(Point2::new(0.0, 0.5 + 1e-7), 2))
            .unwrap();
        let q = st
            .insert(boundary_point(Point2::new(1.0, 0.5 - 1e-7), 3))
            .unwrap();
        assert!(st.constrain(&k, s, p, q).is_err());
        // a clean crossing is split, and the triangulation stays usable
        let x = st.insert(boundary_point(Point2::new(0.5, 0.0), 4)).unwrap();
        let y = st.insert(boundary_point(Point2::new(0.5, 1.0), 5)).unwrap();
        st.constrain(&k, s, x, y).unwrap();
        assert!(st.stale_hints && st.split_boundary);
        // near-duplicate interior points are skipped, distinct ones go in
        let uv = Point2::new(0.3, 0.3);
        assert!(st.insert_interior(&k, s, uv, None, 0).unwrap().is_some());
        assert!(st.insert_interior(&k, s, uv, None, 0).unwrap().is_none());
        let nudged = Point2::new(0.3 + 1e-15, 0.3);
        assert!(st
            .insert_interior(&k, s, nudged, None, 0)
            .unwrap()
            .is_none());
        assert!(st
            .insert_interior(&k, s, Point2::new(0.6, 0.6), None, 0)
            .unwrap()
            .is_some());
        // non-finite positions are refused, not located
        assert!(matches!(
            st.locate(spade::Point2::new(f64::NAN, 0.0)),
            PositionInTriangulation::NoTriangulation
        ));
        assert!(st
            .insert_interior(&k, s, Point2::new(f64::NAN, 0.0), None, 0)
            .is_err());
        // a new triangulation of the same points has a whole hierarchy again
        st.retriangulate(&k, s, &BTreeSet::new()).unwrap();
        assert!(!st.stale_hints);
        assert!(matches!(
            st.locate(st.scaled(Point2::new(0.6, 0.6))),
            PositionInTriangulation::OnVertex(_)
        ));
    }

    #[test]
    fn flood_fill_marks_the_inside_of_a_square_with_a_hole() {
        let mut cdt = Cdt::new();
        let ring = |pts: &[(f64, f64)], base: usize, cdt: &mut Cdt| {
            let hs: Vec<_> = pts
                .iter()
                .enumerate()
                .map(|(i, (x, y))| {
                    cdt.insert(CdtVertex {
                        pos: spade::Point2::new(*x, *y),
                        id: base + i,
                    })
                    .unwrap()
                })
                .collect();
            for i in 0..hs.len() {
                cdt.add_constraint(hs[i], hs[(i + 1) % hs.len()]);
            }
        };
        ring(
            &[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)],
            0,
            &mut cdt,
        );
        ring(
            &[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0)],
            4,
            &mut cdt,
        );
        let inside = inside_faces(&cdt);
        let mut area_in = 0.0;
        for f in cdt.inner_faces() {
            if is_inside(&inside, f.fix()) {
                area_in += f.area();
            }
        }
        assert!((area_in - 12.0).abs() < 1e-12, "{area_in}");
    }
}
