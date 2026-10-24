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
use solidcore::math::{Point2, Point3, Vec3};
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
/// deterministic: a vertex's levels follow from its index).
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
        let natural = surface.natural_uv_domain();
        let mut domain = b.face(f)?.uv_domain;
        if !(domain.u_min() < domain.u_max() && domain.v_min() < domain.v_max()) {
            domain = natural;
        }
        let scale = surface_scales(k, surface, &domain);
        let tol = 1e-9 * model_size.max(1e-300);
        let period = (
            closes(k, surface, &natural, true, tol).then(|| natural.u_max() - natural.u_min()),
            closes(k, surface, &natural, false, tol).then(|| natural.v_max() - natural.v_min()),
        );
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
        })
    }

    fn scaled(&self, uv: Point2) -> spade::Point2<f64> {
        spade::mitigate_underflow(spade::Point2::new(
            (uv.x - self.origin.x) * self.scale.0,
            (uv.y - self.origin.y) * self.scale.1,
        ))
    }

    fn unscaled(&self, p: spade::Point2<f64>) -> Point2 {
        Point2::new(
            self.origin.x + p.x / self.scale.0,
            self.origin.y + p.y / self.scale.1,
        )
    }

    /// The surface point at `uv` (from the B-spline directly when checked).
    fn point(&self, k: &Kernel, s: &Surface, uv: Point2) -> Result<Point3> {
        if let (Some(nb), Some(_)) = (&self.nurbs, self.fast) {
            if let Some((p, _, _)) = nb.eval(uv.x, uv.y).filter(|e| crate::geom::finite3(e.0)) {
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
        self.insert_existing(id)
    }

    /// Puts vertex `id` of the table into the triangulation.
    fn insert_existing(&mut self, id: usize) -> Result<FixedVertexHandle> {
        let uv = self
            .vertices
            .get(id)
            .ok_or_else(|| invalid("vertex out of range"))?
            .uv;
        let pos = self.scaled(uv);
        let h = self
            .cdt
            .insert(CdtVertex { pos, id })
            .map_err(|_| invalid("UV point cannot be triangulated"))?;
        if let Some(v) = self.vertices.get_mut(id) {
            v.alive = true;
        }
        Ok(h)
    }

    /// Inserts an interior point at `uv` with its surface point and normal.
    fn insert_interior(
        &mut self,
        k: &Kernel,
        s: &Surface,
        uv: Point2,
        xyz: Option<(Point3, Option<Vec3>)>,
        depth: u8,
    ) -> Result<FixedVertexHandle> {
        let (xyz, normal) = match xyz {
            Some(x) => x,
            None => (self.point(k, s, uv)?, self.outward(k, s, uv)),
        };
        self.insert(FaceVertex {
            uv,
            xyz,
            vertex: None,
            normal,
            boundary: false,
            alive: true,
            gap: 0.0,
            depth,
        })
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
        if let PositionInTriangulation::OnVertex(h) = self.cdt.locate(self.scaled(uv)) {
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
        let edges = self.cdt.add_constraint_and_split(a, c, |pos| CdtVertex {
            pos,
            id: usize::MAX,
        });
        if edges.is_empty() {
            return Err(invalid("boundary constraint cannot be added"));
        }
        self.split_boundary = true;
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
        let loops = self.face_loops(k, b, s, edges, p)?;
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

        if first || remetric || !self.update_boundary(k, s, &loops, &expanded)? {
            self.rebuild(k, s, &loops, &expanded)?;
        }
        self.boundary = expanded
            .iter()
            .map(|l| l.iter().map(|q| (q.0, q.2)).collect())
            .collect();
        self.flip = self.uv_ccw_faces_inward(k, s, &bbox)?;
        if !first {
            self.remove_stale_interior();
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
                let near = match self.cdt.locate(pos) {
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
                inserted.push(self.insert_interior(k, s, uv, None, depth)?);
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

    /// The boundary loops of the face in UV, outer first, each closed implicitly.
    fn face_loops(
        &mut self,
        k: &Kernel,
        b: &Brep,
        s: &Surface,
        edges: &EdgeSet,
        p: &FaceParams,
    ) -> Result<Vec<Vec<LoopPoint>>> {
        let domain = self.domain;
        let scale = self.scale;
        let span = ((domain.u_max() - domain.u_min()) * scale.0)
            .max((domain.v_max() - domain.v_min()) * scale.1);
        let polish_above = POLISH_FRACTION * p.chord;
        let fu = b
            .face(self.face)?
            .faceuse
            .ok_or_else(|| invalid("face without face use"))?;
        let mut loops = Vec::new();
        let mut inverted = Vec::new();
        for lu in eu1::faceuse_loopuses(b, fu, false, None)? {
            let mut pts: Vec<LoopPoint> = Vec::new();
            let eus = eu1::loopuse_edgeuses(b, lu)?;
            // start at an edge use with a UV trim curve, so edge uses mapped by point
            // inversion continue from a known side of any seam
            let mut start = 0;
            for (i, eu) in eus.iter().enumerate() {
                if eu1::uv_trim_curve(b, *eu)?.0.is_some() {
                    start = i;
                    break;
                }
            }
            for eu in eus.iter().cycle().skip(start).take(eus.len()) {
                let eu = *eu;
                let e = eu1::edge(b, eu)?.ok_or_else(|| invalid("edge use without edge"))?;
                self.edges_used.insert(e);
                let poly = edges
                    .polyline(e)
                    .ok_or_else(|| invalid("face uses an edge that failed"))?;
                let ivl = b.edge(e)?.interval;
                let trim = eu1::uv_trim_curve(b, eu)?.0;
                let reverse = eu1::orientation(b, eu)? == OrientType::Opposite;
                let seed = pts.last().map(|q| q.uv);
                let ctx = MapContext {
                    k,
                    s,
                    domain: &domain,
                    edges,
                    polish_above,
                    period: self.period,
                };
                let mapped = match trim {
                    Some(c) => ctx.edgeuse_points(poly, &ivl, Some(c), reverse, seed).ok(),
                    None => None,
                };
                let seq = if let Some(seq) = mapped {
                    seq
                } else {
                    inverted.push(eu);
                    ctx.edgeuse_points(poly, &ivl, None, reverse, seed)?
                };
                join(&mut pts, seq, k, s, p.chord, span, scale)?;
            }
            // close the loop
            if pts.len() > 1 {
                if let (Some(first), Some(last)) = (pts.first().copied(), pts.last().copied()) {
                    if first.vertex != last.vertex {
                        return Err(invalid("loop does not close"));
                    }
                    if junction_is_pole(k, s, p.chord, last, first, span, scale)? {
                        if let Some(l) = pts.last_mut() {
                            l.pole_next = true;
                        }
                    } else {
                        pts.pop();
                    }
                }
            }
            if pts.len() >= 2 {
                loops.push(pts);
            }
        }
        for eu in inverted {
            if !self.inverted.contains(&eu) {
                self.inverted.push(eu);
            }
        }
        Ok(loops)
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
                let key = (uv.x.to_bits(), uv.y.to_bits());
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
                    let pos = self.scaled(*uv);
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
        let ids: Vec<usize> = points.iter().map(|p| p.id).collect();
        let mut conflicts = Vec::new();
        self.cdt = Cdt::try_bulk_load_cdt(points, segments, |e| conflicts.push(e))
            .map_err(|_| invalid("UV boundary cannot be triangulated"))?;
        // a segment crossing another is added with splits
        for [a, c] in conflicts {
            let (Some(ia), Some(ic)) = (ids.get(a), ids.get(c)) else {
                continue;
            };
            let (Some(ha), Some(hc)) = (self.handle_of(*ia), self.handle_of(*ic)) else {
                continue;
            };
            self.constrain(k, s, ha, hc)?;
        }
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
        match self.cdt.locate(self.scaled(uv)) {
            PositionInTriangulation::OnVertex(h) => Some(h),
            _ => None,
        }
    }

    fn find_boundary(&self, uv: Point2, vertex: u32) -> Result<Option<FixedVertexHandle>> {
        match self.cdt.locate(self.scaled(uv)) {
            PositionInTriangulation::OnVertex(h) if self.vertex(h)?.vertex == Some(vertex) => {
                Ok(Some(h))
            }
            _ => Ok(None),
        }
    }

    /// Removes interior points the refined boundary leaves outside, or closer to it
    /// than a tenth of their shortest edge.
    fn remove_stale_interior(&mut self) {
        let inside = inside_faces(&self.cdt);
        let mut stale = Vec::new();
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
                stale.push(pos);
            }
        }
        for pos in stale {
            if let PositionInTriangulation::OnVertex(h) = self.cdt.locate(pos) {
                let id = self.id(h);
                self.cdt.remove(h);
                if let Some(v) = self.vertices.get_mut(id) {
                    v.alive = false;
                }
            }
        }
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
        let mut wanted = Vec::new();
        for (j, v) in vs.iter().enumerate() {
            for (i, u) in us.iter().enumerate() {
                let uv = Point2::new(*u, *v);
                let pos = self.scaled(uv);
                let PositionInTriangulation::OnFace(fh) = self.cdt.locate(pos) else {
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
            let (Some(u), Some(v)) = (us.get(i), vs.get(j)) else {
                continue;
            };
            self.insert_interior(k, s, Point2::new(*u, *v), None, 0)?;
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
        let faces = match self.cdt.locate(pos) {
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
        let f = match self.cdt.locate(pos) {
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
        let PositionInTriangulation::OnFace(f) = self.cdt.locate(pos) else {
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
            let uv = Point2::new(
                (a.uv.x + b.uv.x + c.uv.x) / 3.0,
                (a.uv.y + b.uv.y + c.uv.y) / 3.0,
            );
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
                    let q = Point2::new(
                        a.uv.x * w[0] + b.uv.x * w[1] + c.uv.x * w[2],
                        a.uv.y * w[0] + b.uv.y * w[1] + c.uv.y * w[2],
                    );
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
            let uv = Point2::new((a.uv.x + b.uv.x) / 2.0, (a.uv.y + b.uv.y) / 2.0);
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
            let uv = Point2::new((a.uv.x + b.uv.x) / 2.0, (a.uv.y + b.uv.y) / 2.0);
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
            let (dev, over) = crate::deviation::face_mesh_deviation(k, s, &self.domain, &m, chord);
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
    domain: &'a Extent2d,
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
    /// kept within half a period of it on a closed surface.
    fn edgeuse_points(
        &self,
        poly: &EdgePolyline,
        ivl: &Extent1d,
        trim: Option<&Curve>,
        reverse: bool,
        seed: Option<Point2>,
    ) -> Result<Vec<LoopPoint>> {
        let (k, s) = (self.k, self.s);
        let (t0, t1) = (ivl.min(), ivl.max());
        let len = t1 - t0;
        let n = poly.params.len();
        let mut out = Vec::with_capacity(n);
        let mut prev = seed;
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
            let uv = match trim {
                Some(c) => {
                    let ti = c.natural_interval();
                    let q = c.evaluate_point(k, ti.evaluate(sn))?;
                    if !(q.x.is_finite() && q.y.is_finite()) {
                        return Err(invalid("UV trim curve is not finite"));
                    }
                    let uv = Point2::new(q.x, q.y);
                    if poly.degenerate
                        || surface_point(k, s, uv)?.distance(xyz) <= self.polish_above
                    {
                        uv
                    } else {
                        foot(k, s, uv, xyz, self.domain).0
                    }
                }
                None => match prev {
                    Some(seed) => {
                        let uv = foot(k, s, seed, xyz, self.domain).0;
                        self.unwrap(uv, seed)
                    }
                    None => invert_point(k, s, xyz, self.domain)?,
                },
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

/// Whether the UV jump from `a` to `b` (one mesh vertex) is a pole: far apart in UV but
/// the surface stays at the vertex between them. A jump across most of the domain where
/// the surface does not stay is a loop that only closes across a seam.
fn junction_is_pole(
    k: &Kernel,
    s: &Surface,
    chord: f64,
    a: LoopPoint,
    b: LoopPoint,
    span: f64,
    scale: (f64, f64),
) -> Result<bool> {
    let gap = scaled_gap(a.uv, b.uv, scale);
    if gap <= chord.max(1e-9 * span) {
        return Ok(false);
    }
    let mid = Point2::new((a.uv.x + b.uv.x) / 2.0, (a.uv.y + b.uv.y) / 2.0);
    if surface_point(k, s, mid)?.distance(a.xyz) <= chord {
        return Ok(true);
    }
    if gap > 0.1 * span {
        return Err(invalid(
            "loop does not close in UV (periodic face without seam edge)",
        ));
    }
    Ok(false)
}

/// Appends an edge use's points to a loop, merging the shared junction point.
fn join(
    pts: &mut Vec<LoopPoint>,
    mut seq: Vec<LoopPoint>,
    k: &Kernel,
    s: &Surface,
    chord: f64,
    span: f64,
    scale: (f64, f64),
) -> Result<()> {
    if let (Some(last), Some(first)) = (pts.last().copied(), seq.first().copied()) {
        if last.vertex != first.vertex {
            return Err(invalid("loop edge uses do not connect"));
        }
        if junction_is_pole(k, s, chord, last, first, span, scale)? {
            if let Some(l) = pts.last_mut() {
                l.pole_next = true;
            }
        } else {
            pts.pop();
            if let Some(f) = seq.first_mut() {
                f.pole_next |= last.pole_next;
            }
        }
    }
    pts.extend(seq);
    Ok(())
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
