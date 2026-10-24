//! Watertight, configurable triangle meshes of solidcore B-reps.
//!
//! This crate is new code, not a port: it sits beside solidcore's ported SMLib
//! tessellator (`solidcore::query::tess`, kept for C++ parity) and is judged against the
//! true geometry instead ([`deviation`]). The design follows the usual trimmed-surface
//! scheme (as in truck's `truck-meshalgo`, Apache-2.0, and the robustness fixes
//! SolidCraft made to it) and the knot-span grid of the nurbs-grid-mesh design note
//! (Apache-2.0); no code is copied. The constrained Delaunay triangulation is the
//! [`spade`] crate (MIT or Apache-2.0).
//!
//! # Pipeline
//!
//! 1. **Edges** ([`EdgeSet`]). Every edge is divided once, adaptively in its curve
//!    parameter (chord at the quarter points and midpoint, tangent turning, the normal
//!    turning of every face along it, segment length, segment counts), from its start
//!    vertex to its end vertex. Every face reuses that polyline, so neighbouring faces
//!    share their boundary vertices exactly: the mesh is watertight by construction.
//! 2. **Loops in UV.** Each edge use maps the shared polyline into its face's UV with
//!    its UV trim curve at the same normalised edge parameter (as upstream
//!    `SmEdgeuse::NormalizedEvaluate`), polished by Gauss-Newton only where the trim's
//!    surface point is off the edge point. An edge use without a UV trim curve falls
//!    back to point inversion, continuing from the previous point and kept within half
//!    a period of it on a closed surface ([`Mesh::inverted_edgeuses`]).
//! 3. **Seams and poles.** A seam is a pair of loop edges at the two sides of the
//!    period and needs no special case. A pole (a degenerate edge, or a UV jump between
//!    edge uses where the surface stays at one point) is sampled in UV at the grid
//!    spacing, every sample on the one pole vertex; triangles collapsing there are
//!    dropped.
//! 4. **Interior.** Grid lines from the chord, angle and length needs of a sample grid
//!    (per knot span from a second-derivative bound on B-splines), with cells of
//!    bounded aspect, give Steiner points inside the loops and clear of them. The
//!    triangulation runs in UV scaled by powers of two near the surface's mean speeds.
//!    Inside triangles are found by flooding from the outer face across constraints.
//! 5. **Refinement.** Bad triangles (off the surface, flipped, or turned past the angle
//!    tolerance) get the midpoint of their longest edge, or their centroid when that
//!    edge is on the boundary; interior edges off the chord or too long get their
//!    midpoint. No point goes onto or next to the boundary, which stays the shared
//!    polyline.
//! 6. **Boundary rounds.** A bad triangle on a boundary segment too long to refine
//!    from inside asks for that segment to be split; the edge polyline gains its
//!    midpoint and the faces along it refine again (inserting points only).
//! 7. **Fallbacks** (SolidCraft's hard cases). A boundary crossing itself in UV is
//!    split where it crosses ([`Mesh::split_boundary_faces`]); a face with no area in
//!    UV is ear-clipped from its loop in 3D ([`Mesh::ear_clipped_faces`]). Where the
//!    B-rep's edge curve lies off its face's surface (an edge/face gap), triangles at
//!    that boundary may be off the surface by the gap on top of the chord tolerance
//!    ([`MeshStats::boundary_gap`]).
//!
//! # Levels of detail
//!
//! [`tessellate_lods`] (or a [`Tessellator`] fed finer and finer options) refines one
//! state: each level inserts points into the previous level's edge polylines and face
//! triangulations, so vertex ids are shared and every level is watertight.
//!
//! # Parallel use
//!
//! The library runs no threads (wasm). A caller can divide the edges once
//! ([`tessellate_edge_set`]) and mesh faces independently, one [`FaceTessellation`]
//! per face, which only reads the edge set, on its own threads; the steps
//! [`FaceTessellation`] lists give exactly what [`tessellate`] returns.
//!
//! Deviations: not a port. Bounds (no recursion anywhere): edge segments split at most
//! [`MAX_EDGE_DEPTH`] times and an edge has at most [`MAX_EDGE_POINTS`] points (more is
//! an error); at most [`MAX_GRID_LINES`] grid lines per direction; refinement points
//! at most [`MAX_LOCAL_DEPTH`] generations below a boundary or grid point; at most
//! [`MAX_BOUNDARY_ROUNDS`] rounds of boundary splits; refinement stops at
//! [`TessellateOptions::max_refine_passes`] passes or
//! [`TessellateOptions::max_triangles_per_face`] triangles (a face that stops at any of
//! these with work left is listed in [`Mesh::unconverged_faces`]). Output is
//! deterministic: every collection is visited in id or insertion order, and the
//! triangulation (with its hierarchy point location) is deterministic.

#![forbid(unsafe_code)]

pub mod backend;
mod deviation;
mod edges;
mod face;
mod geom;
mod grid;
mod mesh;
mod nurbs;
mod options;

pub use deviation::{deviation, DeviationReport};
pub use edges::{EdgeSet, SegmentSplit, MAX_EDGE_DEPTH, MAX_EDGE_POINTS, MAX_SPLIT_PIECES};
pub use face::{FaceMesh, MAX_LOCAL_DEPTH};
pub use grid::MAX_GRID_LINES;
pub use mesh::{EdgePolyline, EdgeReport, Failure, Mesh, MeshStats};
pub use options::{
    EdgeHook, EdgeOptions, EdgeOverride, FaceHook, FaceOverride, OutputOptions, Quality,
    TessellateOptions, Tolerance,
};

use solidcore::kernel::Kernel;
use solidcore::math::Vec3;
use solidcore::topology::brep::p1 as bl;
use solidcore::topology::Brep;
use solidcore::{EdgeId, EntityRef, FaceId, Result, SolidError};

use crate::edges::push_position;
use crate::face::FaceState;
use crate::geom::invalid;

/// Tessellates every face of `b` into one mesh within `opts`.
///
/// Edges and faces that cannot be meshed are skipped and listed in [`Mesh::failures`]
/// (the mesh then has holes there).
///
/// # Errors
///
/// [`SolidError::InvalidInput`] for options out of range; [`SolidError::NotFound`] for
/// a B-rep whose lists name missing entities.
pub fn tessellate(k: &Kernel, b: &Brep, opts: &TessellateOptions) -> Result<Mesh> {
    Tessellator::new(k, b)?.level(opts)
}

/// Nested levels of detail, one mesh per entry of `levels` (coarsest first): each
/// level refines the previous one's edge polylines and face triangulations by inserting
/// points, so the levels share vertex ids and each is watertight. A level asking for
/// less than the one before stays as fine as that one. Each level's
/// [`MeshStats::max_deviation`] is measured.
///
/// # Errors
///
/// As [`tessellate`].
pub fn tessellate_lods(k: &Kernel, b: &Brep, levels: &[TessellateOptions]) -> Result<Vec<Mesh>> {
    let mut t = Tessellator::new(k, b)?;
    levels
        .iter()
        .map(|o| {
            let mut o = o.clone();
            o.measure_deviation = true;
            t.level(&o)
        })
        .collect()
}

/// Only the edge polylines of `b` (a wireframe; no face is meshed): the mesh has no
/// triangles, and its positions are the polylines' points.
///
/// # Errors
///
/// As [`tessellate`].
pub fn tessellate_edges(k: &Kernel, b: &Brep, opts: &TessellateOptions) -> Result<Mesh> {
    Tessellator::new(k, b)?.edge_level(opts)
}

/// The edge polylines of `b`, the shared input of [`tessellate_face`].
///
/// # Errors
///
/// As [`tessellate`].
pub fn tessellate_edge_set(k: &Kernel, b: &Brep, opts: &TessellateOptions) -> Result<EdgeSet> {
    check_options(opts)?;
    let mut edges = EdgeSet {
        model_size: model_size(k, b)?,
        ..EdgeSet::default()
    };
    edges.refine(k, b, opts)?;
    Ok(edges)
}

/// The mesh of face `f` against the edge polylines `edges` (made with the same
/// options), in one call: [`FaceTessellation::new`], one [`FaceTessellation::refine`]
/// and [`FaceTessellation::mesh`], without the boundary rounds of [`tessellate`].
///
/// # Errors
///
/// What stops the face (its UV loops, a surface that does not evaluate, ...).
pub fn tessellate_face(
    k: &Kernel,
    b: &Brep,
    f: FaceId,
    edges: &EdgeSet,
    opts: &TessellateOptions,
) -> Result<FaceMesh> {
    let mut t = FaceTessellation::new(k, b, f, edges)?;
    t.refine(k, b, edges, opts)?;
    t.mesh(k, b, edges, opts)
}

/// Most rounds of splitting boundary segments faces ask for
/// ([`FaceTessellation::boundary_requests`]) and refining the faces along them.
pub const MAX_BOUNDARY_ROUNDS: usize = 4;

/// One face's tessellation, kept so it can be refined again: after its edges gained
/// points (a finer level of detail, or [`EdgeSet::split_segments`]) a refinement only
/// inserts points. Faces are independent given the edge set, so a caller may hold one
/// per face and refine them on several threads; [`tessellate`] is
///
/// 1. [`tessellate_edge_set`];
/// 2. per face, [`FaceTessellation::new`] and [`FaceTessellation::refine`];
/// 3. up to [`MAX_BOUNDARY_ROUNDS`] times: every face's
///    [`FaceTessellation::boundary_requests`] (sorted) into
///    [`EdgeSet::split_segments`], then [`FaceTessellation::refine`] again for the
///    faces that [use](FaceTessellation::uses_edge) a changed edge;
/// 4. per face, [`FaceTessellation::mesh`], joined in face order by [`assemble`].
pub struct FaceTessellation {
    state: FaceState,
}

impl FaceTessellation {
    /// A tessellation of face `f` with nothing meshed yet.
    ///
    /// # Errors
    ///
    /// [`SolidError::InvalidInput`] for a face without surface or face use.
    pub fn new(k: &Kernel, b: &Brep, f: FaceId, edges: &EdgeSet) -> Result<Self> {
        Ok(Self {
            state: FaceState::new(k, b, f, edges.model_size)?,
        })
    }

    /// Meshes the face to `opts` against `edges`, or refines what an earlier call made.
    ///
    /// # Errors
    ///
    /// What stops the face (its UV loops, a surface that does not evaluate, ...).
    pub fn refine(
        &mut self,
        k: &Kernel,
        b: &Brep,
        edges: &EdgeSet,
        opts: &TessellateOptions,
    ) -> Result<()> {
        check_options(opts)?;
        let p = opts.face_params(self.state.face(), edges.model_size);
        self.state.refine(k, b, edges, &p, opts.backend.as_deref())
    }

    /// Boundary segments (sorted mesh vertex pairs, with the pieces each should become)
    /// the last refinement could not resolve from inside the face; see
    /// [`EdgeSet::split_segments`].
    #[must_use]
    pub fn boundary_requests(&self) -> &[SegmentSplit] {
        self.state.requests()
    }

    /// Whether edge `e` bounds the face.
    #[must_use]
    pub fn uses_edge(&self, e: EdgeId) -> bool {
        self.state.uses_edge(e)
    }

    /// The face's mesh now (its deviation measured when `opts` asks for it).
    ///
    /// # Errors
    ///
    /// [`SolidError::InvalidInput`] for a face too large to number.
    pub fn mesh(
        &self,
        k: &Kernel,
        b: &Brep,
        edges: &EdgeSet,
        opts: &TessellateOptions,
    ) -> Result<FaceMesh> {
        let p = opts.face_params(self.state.face(), edges.model_size);
        self.state.snapshot(
            k,
            b,
            opts.measure_deviation.then_some(p.chord),
            opts.backend.as_deref(),
        )
    }

    /// Keeps the mesh vertices [`assemble`] gave this face's interior vertices, so a
    /// later refinement's meshes reuse them.
    pub fn adopt(&mut self, m: &FaceMesh) {
        self.state.adopt(m);
    }
}

/// Joins face meshes (in the order given; give them in the B-rep's face order to match
/// [`tessellate`]) and the edge polylines into one mesh; a face's `Err` becomes a
/// [`Mesh::failures`] entry. Edge deviation is measured when `opts` asks for it and
/// `k` and `b` are given.
pub fn assemble(
    edges: EdgeSet,
    faces: Vec<(FaceId, Result<FaceMesh>)>,
    opts: &TessellateOptions,
    measure: Option<(&Kernel, &Brep)>,
) -> Mesh {
    let mut edges = edges;
    let mut mesh = Mesh::default();
    mesh.failures.clone_from(&edges.failures);
    for (f, r) in faces {
        match r {
            Ok(mut fm) => append_face(&mut mesh, &mut edges.positions, &mut fm, opts),
            Err(error) => mesh.failures.push(Failure {
                entity: EntityRef::Face(f),
                error,
            }),
        }
    }
    finish(&mut mesh, &edges, opts, measure);
    mesh
}

fn check_options(opts: &TessellateOptions) -> Result<()> {
    let positive = |x: f64| x.is_finite() && x > 0.0;
    let chord = match opts.chord_tolerance {
        Tolerance::Absolute(x) | Tolerance::Relative(x) => x,
    };
    if !positive(chord) {
        return Err(invalid("chord tolerance must be finite and positive"));
    }
    if !(opts.angle_tolerance > 0.0 && opts.angle_tolerance <= core::f64::consts::PI) {
        return Err(invalid("angle tolerance must be in (0, pi]"));
    }
    let lengths = [opts.max_edge_length, opts.edge.max_segment_length];
    if lengths.iter().flatten().any(|l| !positive(*l)) {
        return Err(invalid("maximum lengths must be finite and positive"));
    }
    if !(opts.min_edge_length.is_finite() && opts.min_edge_length >= 0.0) {
        return Err(invalid(
            "minimum edge length must be finite and not negative",
        ));
    }
    Ok(())
}

/// The model size relative tolerances resolve against: the diagonal of the box of the
/// B-rep's vertices and 17 points along each edge curve. An estimate of the
/// bounding-box diagonal (low for a face that bulges past its edges, such as a whole
/// sphere's) that is cheap where the tight box query is not; 1 when it is empty.
fn model_size(k: &Kernel, b: &Brep) -> Result<f64> {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    let mut add = |p: solidcore::math::Point3| {
        if !(p.x.is_finite() && p.y.is_finite() && p.z.is_finite()) {
            return;
        }
        for ((l, h), x) in lo.iter_mut().zip(hi.iter_mut()).zip([p.x, p.y, p.z]) {
            *l = l.min(x);
            *h = h.max(x);
        }
    };
    for v in bl::vertices(b, None)? {
        add(b.vertex(v)?.point);
    }
    for e in bl::edges(b, None)? {
        let rec = b.edge(e)?;
        if let Some(c) = rec.curve.as_ref() {
            for i in 0..=16 {
                if let Ok(p) = c.evaluate_point(k, rec.interval.evaluate(f64::from(i) / 16.0)) {
                    add(p);
                }
            }
        }
    }
    let d = ((hi[0] - lo[0]).powi(2) + (hi[1] - lo[1]).powi(2) + (hi[2] - lo[2]).powi(2)).sqrt();
    Ok(if d.is_finite() && d > 0.0 { d } else { 1.0 })
}

/// Adds one face mesh: its interior vertices get mesh vertices (recorded in
/// `fm.shared`), its triangles and per-corner data join the mesh.
fn append_face(
    mesh: &mut Mesh,
    positions: &mut Vec<solidcore::math::Point3>,
    fm: &mut FaceMesh,
    opts: &TessellateOptions,
) {
    for (sh, p) in fm.shared.iter_mut().zip(&fm.positions) {
        if sh.is_none() {
            if let Ok(i) = push_position(positions, *p) {
                *sh = Some(i);
            }
        }
    }
    for t in &fm.triangles {
        let get = |i: u32| fm.shared.get(i as usize).copied().flatten();
        let (Some(a), Some(b), Some(c)) = (get(t[0]), get(t[1]), get(t[2])) else {
            continue;
        };
        if a == b || b == c || c == a {
            continue;
        }
        mesh.triangles.push([a, b, c]);
        if opts.output.normals {
            let mut n = t.map(|i| {
                fm.normals
                    .get(i as usize)
                    .copied()
                    .unwrap_or(Vec3::new(0.0, 0.0, 0.0))
            });
            // a corner without surface normal, or whose surface normal points against
            // the triangle (a singular point or crease at that corner), takes the
            // triangle's own
            let p = t.map(|i| {
                fm.positions
                    .get(i as usize)
                    .copied()
                    .unwrap_or(Vec3::new(0.0, 0.0, 0.0))
            });
            let g = (p[1] - p[0]).cross(p[2] - p[0]);
            for x in &mut n {
                if x.length_squared() == 0.0 || x.dot(g) <= 0.0 {
                    *x = g.unitize_or_z();
                }
            }
            mesh.normals.push(n);
        }
        if opts.output.uvs {
            mesh.uvs.push(t.map(|i| {
                fm.uvs
                    .get(i as usize)
                    .copied()
                    .unwrap_or(solidcore::math::Point2::new(0.0, 0.0))
            }));
        }
        if opts.output.face_ids {
            mesh.triangle_faces.push(fm.face);
        }
    }
    mesh.inverted_edgeuses
        .extend(fm.inverted_edgeuses.iter().copied());
    if fm.unconverged {
        mesh.unconverged_faces.push(fm.face);
    }
    if fm.split_boundary {
        mesh.split_boundary_faces.push(fm.face);
    }
    if fm.backend.is_some() {
        mesh.backend_faces.push(fm.face);
    }
    if fm.backend_failed {
        mesh.backend_failed_faces.push(fm.face);
    }
    if let Some(d) = fm.backend_deviation {
        mesh.stats.backend_deviation =
            Some(mesh.stats.backend_deviation.map_or(d, |x: f64| x.max(d)));
    }
    if fm.ear_clipped {
        mesh.ear_clipped_faces.push(fm.face);
    }
    if fm.deviation.is_some() {
        let g = fm.gaps.iter().fold(0.0_f64, |x, y| x.max(*y));
        mesh.stats.boundary_gap = Some(mesh.stats.boundary_gap.map_or(g, |x: f64| x.max(g)));
    }
    if let Some(d) = fm.deviation {
        mesh.stats.face_deviation = Some(mesh.stats.face_deviation.map_or(d, |x: f64| x.max(d)));
        mesh.stats.triangles_over += fm.triangles_over;
    }
}

/// Positions, edge polylines and statistics once every face is in.
fn finish(
    mesh: &mut Mesh,
    edges: &EdgeSet,
    opts: &TessellateOptions,
    measure: Option<(&Kernel, &Brep)>,
) {
    mesh.positions.clone_from(&edges.positions);
    if opts.output.edge_polylines {
        mesh.edges.clone_from(&edges.edges);
    }
    mesh.stats.model_size = edges.model_size;
    mesh.stats.chord_tolerance = opts.chord_tolerance.resolve(edges.model_size);
    if let (true, Some((k, b))) = (opts.measure_deviation, measure) {
        let mut worst = 0.0_f64;
        for e in &edges.edges {
            if let Ok(d) = deviation::edge_deviation(k, b, e, &edges.positions) {
                worst = worst.max(d);
            }
        }
        mesh.stats.edge_deviation = Some(worst);
        mesh.stats.face_deviation = Some(mesh.stats.face_deviation.unwrap_or(0.0));
        mesh.stats.max_deviation = Some(worst.max(mesh.stats.face_deviation.unwrap_or(0.0)));
    }
}

/// One face's place in a [`Tessellator`].
struct FaceSlot {
    face: FaceId,
    state: Option<FaceState>,
    failure: Option<SolidError>,
}

/// A tessellation kept between levels of detail: each [`Tessellator::level`] refines
/// what the previous one made, inserting points only.
pub struct Tessellator<'a> {
    k: &'a Kernel,
    b: &'a Brep,
    edges: EdgeSet,
    faces: Vec<FaceSlot>,
}

impl<'a> Tessellator<'a> {
    /// A tessellation of `b` with nothing meshed yet.
    ///
    /// # Errors
    ///
    /// [`SolidError::NotFound`] for a B-rep whose lists name missing entities.
    pub fn new(k: &'a Kernel, b: &'a Brep) -> Result<Self> {
        let faces = bl::faces(b, None)?
            .into_iter()
            .map(|face| FaceSlot {
                face,
                state: None,
                failure: None,
            })
            .collect();
        Ok(Self {
            k,
            b,
            edges: EdgeSet {
                model_size: model_size(k, b)?,
                ..EdgeSet::default()
            },
            faces,
        })
    }

    /// The edge polylines so far.
    #[must_use]
    pub fn edges(&self) -> &EdgeSet {
        &self.edges
    }

    /// Refines the edges only, to `opts`, and returns them as a mesh without triangles.
    ///
    /// # Errors
    ///
    /// [`SolidError::InvalidInput`] for options out of range.
    pub fn edge_level(&mut self, opts: &TessellateOptions) -> Result<Mesh> {
        check_options(opts)?;
        self.edges.refine(self.k, self.b, opts)?;
        let mut mesh = Mesh::default();
        mesh.failures.clone_from(&self.edges.failures);
        let mut o = opts.clone();
        o.output.edge_polylines = true;
        finish(&mut mesh, &self.edges, &o, Some((self.k, self.b)));
        Ok(mesh)
    }

    /// The next level: edges and faces refined to `opts` (inserting points only) and
    /// the whole mesh returned.
    ///
    /// # Errors
    ///
    /// [`SolidError::InvalidInput`] for options out of range.
    pub fn level(&mut self, opts: &TessellateOptions) -> Result<Mesh> {
        check_options(opts)?;
        let (k, b) = (self.k, self.b);
        self.edges.refine(k, b, opts)?;
        let size = self.edges.model_size;
        let refine = |slot: &mut FaceSlot, edges: &EdgeSet| {
            if slot.failure.is_some() {
                return;
            }
            let p = opts.face_params(slot.face, size);
            let run = |slot: &mut FaceSlot| -> Result<()> {
                if slot.state.is_none() {
                    slot.state = Some(FaceState::new(k, b, slot.face, size)?);
                }
                slot.state
                    .as_mut()
                    .ok_or_else(|| invalid("face state missing"))?
                    .refine(k, b, edges, &p, opts.backend.as_deref())
            };
            if let Err(error) = run(slot) {
                slot.state = None;
                slot.failure = Some(error);
            }
        };
        for slot in &mut self.faces {
            refine(slot, &self.edges);
        }
        // boundary segments faces could not refine from inside are split in the edge
        // polylines, and the faces along them refined again
        for _ in 0..MAX_BOUNDARY_ROUNDS {
            let mut requests: Vec<SegmentSplit> = self
                .faces
                .iter()
                .filter_map(|s| s.state.as_ref())
                .flat_map(|s| s.requests().iter().copied())
                .collect();
            requests.sort_unstable();
            requests.dedup();
            if requests.is_empty() {
                break;
            }
            let changed = self.edges.split_segments(k, b, &requests)?;
            if changed.is_empty() {
                break;
            }
            for slot in &mut self.faces {
                let uses = slot
                    .state
                    .as_ref()
                    .is_some_and(|s| changed.iter().any(|e| s.uses_edge(*e)));
                if uses {
                    refine(slot, &self.edges);
                }
            }
        }
        let mut mesh = Mesh::default();
        mesh.failures.clone_from(&self.edges.failures);
        for slot in &mut self.faces {
            if let Some(error) = &slot.failure {
                mesh.failures.push(Failure {
                    entity: EntityRef::Face(slot.face),
                    error: error.clone(),
                });
                continue;
            }
            let p = opts.face_params(slot.face, size);
            let Some(st) = slot.state.as_mut() else {
                continue;
            };
            match st.snapshot(
                k,
                b,
                opts.measure_deviation.then_some(p.chord),
                opts.backend.as_deref(),
            ) {
                Ok(mut fm) => {
                    append_face(&mut mesh, &mut self.edges.positions, &mut fm, opts);
                    st.adopt(&fm);
                }
                Err(error) => {
                    mesh.failures.push(Failure {
                        entity: EntityRef::Face(slot.face),
                        error: error.clone(),
                    });
                    slot.state = None;
                    slot.failure = Some(error);
                }
            }
        }
        finish(&mut mesh, &self.edges, opts, Some((k, b)));
        Ok(mesh)
    }
}
