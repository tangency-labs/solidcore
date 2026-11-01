//! The output: [`Mesh`], its edge polylines and the checks run on it.

use std::collections::BTreeMap;

use solidcore::math::{Point2, Point3, Vec3};
use solidcore::topology::EdgeuseId;
use solidcore::{EdgeId, EntityRef, FaceId, SolidError};

/// The shared polyline of one B-rep edge.
#[derive(Debug, Clone, PartialEq)]
pub struct EdgePolyline {
    /// The edge.
    pub edge: EdgeId,
    /// Curve parameters of the points, increasing over the edge's interval.
    pub params: Vec<f64>,
    /// Mesh vertex of each point (the first and last are the edge's vertices).
    pub vertices: Vec<u32>,
    /// The edge collapses to its vertex (a pole): two points on one mesh vertex.
    pub degenerate: bool,
}

/// An entity that could not be meshed, and why.
#[derive(Debug, Clone, PartialEq)]
pub struct Failure {
    /// The edge or face.
    pub entity: EntityRef,
    /// What stopped it.
    pub error: SolidError,
}

/// What one tessellation run measured and how hard it worked.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct MeshStats {
    /// The chord tolerance asked for, in model units (the global one; overrides aside).
    pub chord_tolerance: f64,
    /// The model size relative tolerances were resolved against (the diagonal of the
    /// box of the B-rep's vertices and edge samples).
    pub model_size: f64,
    /// Largest measured distance between the mesh and the true geometry: triangles
    /// against their face's surface (seven samples per triangle, each to its Newton foot
    /// on the surface) and edge segments against their curve. `None` unless
    /// [`crate::TessellateOptions::measure_deviation`] was set.
    pub max_deviation: Option<f64>,
    /// Of [`MeshStats::max_deviation`], the part from faces.
    pub face_deviation: Option<f64>,
    /// Of [`MeshStats::max_deviation`], the part from edge segments.
    pub edge_deviation: Option<f64>,
    /// Largest deviation of the faces whose grid a backend placed, by the backend's own
    /// measure (when measured).
    pub backend_deviation: Option<f64>,
    /// Triangles whose measured deviation exceeds their face's chord tolerance plus the
    /// gap of their corners from the surface (the B-rep's own edge/face gaps).
    pub triangles_over: usize,
    /// Largest distance of a face's boundary vertex (an edge polyline point) from the
    /// face's surface: the B-rep's own edge/face gap, which no mesh sharing the edge
    /// polylines can close. Measured with the deviation.
    pub boundary_gap: Option<f64>,
}

/// A triangle mesh of a B-rep. Vertex ids are stable across the levels of
/// [`crate::tessellate_lods`]: a finer level's positions start with the coarser
/// level's (some of which a finer level may no longer use).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Mesh {
    /// Vertex positions.
    pub positions: Vec<Point3>,
    /// Triangles as vertex indices, counter-clockwise seen from the outward side.
    pub triangles: Vec<[u32; 3]>,
    /// Unit outward surface normal at each triangle's corners (parallel to
    /// [`Mesh::triangles`]; empty when not asked for). A corner takes its own face's
    /// normal, so a vertex on a sharp edge has one normal per face; where the surface
    /// has none, or one pointing against the triangle (a singular point), the corner
    /// takes the triangle's normal.
    pub normals: Vec<[Vec3; 3]>,
    /// The face's (u, v) at each triangle's corners (parallel to [`Mesh::triangles`];
    /// empty when not asked for).
    pub uvs: Vec<[Point2; 3]>,
    /// The B-rep face of each triangle (parallel to [`Mesh::triangles`]; empty when not
    /// asked for).
    pub triangle_faces: Vec<FaceId>,
    /// The polyline of every tessellated edge, in the B-rep's edge order (empty when not
    /// asked for).
    pub edges: Vec<EdgePolyline>,
    /// Edges and faces left out, in the order met.
    pub failures: Vec<Failure>,
    /// Edge uses mapped into UV by point inversion (no usable UV trim curve).
    pub inverted_edgeuses: Vec<EdgeuseId>,
    /// Faces whose refinement stopped at a bound (passes or triangles) with work left.
    pub unconverged_faces: Vec<FaceId>,
    /// Faces whose boundary crossed itself in UV and was split where it crossed: the
    /// split points are not on the neighbouring face, so the mesh is closed there only
    /// geometrically.
    pub split_boundary_faces: Vec<FaceId>,
    /// Faces with no area in UV that were filled by clipping ears off their boundary
    /// loop in 3D.
    pub ear_clipped_faces: Vec<FaceId>,
    /// Faces whose grid lines a [`crate::backend::FaceBackend`] placed.
    pub backend_faces: Vec<FaceId>,
    /// Faces a backend took and failed on; the built-in grid was used.
    pub backend_failed_faces: Vec<FaceId>,
    /// Measurements and effort.
    pub stats: MeshStats,
}

/// How the triangle edges of a [`Mesh`] are shared ([`Mesh::edge_report`]).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct EdgeReport {
    /// Distinct undirected triangle edges.
    pub edges: usize,
    /// Edges of exactly one triangle (the boundary of an open mesh).
    pub open: Vec<[u32; 2]>,
    /// Edges of more than two triangles.
    pub non_manifold: Vec<[u32; 2]>,
    /// Edges of two triangles that run the same way along them (inconsistent winding).
    pub inconsistent: Vec<[u32; 2]>,
}

impl EdgeReport {
    /// Every edge has two triangles running opposite ways: a closed, consistently
    /// oriented mesh.
    #[must_use]
    pub fn is_watertight(&self) -> bool {
        self.open.is_empty() && self.non_manifold.is_empty() && self.inconsistent.is_empty()
    }
}

impl Mesh {
    /// The triangle's corner positions, `None` for an index out of range.
    #[must_use]
    pub fn triangle_points(&self, t: [u32; 3]) -> Option<[Point3; 3]> {
        Some([
            *self.positions.get(t[0] as usize)?,
            *self.positions.get(t[1] as usize)?,
            *self.positions.get(t[2] as usize)?,
        ])
    }

    /// Total triangle area.
    #[must_use]
    pub fn area(&self) -> f64 {
        self.triangles
            .iter()
            .filter_map(|t| self.triangle_points(*t))
            .map(|[a, b, c]| (b - a).cross(c - a).length() / 2.0)
            .sum()
    }

    /// Signed enclosed volume (divergence theorem; positive for a closed mesh wound
    /// outward).
    #[must_use]
    pub fn volume(&self) -> f64 {
        self.triangles
            .iter()
            .filter_map(|t| self.triangle_points(*t))
            .map(|[a, b, c]| a.dot(b.cross(c)) / 6.0)
            .sum()
    }

    /// Mesh vertices used by at least one triangle or edge polyline.
    #[must_use]
    pub fn used_vertices(&self) -> usize {
        let mut used = vec![false; self.positions.len()];
        let ids = self
            .triangles
            .iter()
            .flatten()
            .chain(self.edges.iter().flat_map(|e| e.vertices.iter()));
        for i in ids {
            if let Some(u) = used.get_mut(*i as usize) {
                *u = true;
            }
        }
        used.into_iter().filter(|u| *u).count()
    }

    /// How the triangle edges are shared: open, non-manifold and inconsistently wound
    /// edges, each as its vertex pair in ascending order.
    #[must_use]
    pub fn edge_report(&self) -> EdgeReport {
        // (count, sum of directions: +1 for low -> high, -1 for high -> low)
        let mut uses: BTreeMap<[u32; 2], (usize, i64)> = BTreeMap::new();
        for t in &self.triangles {
            for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
                let (key, dir) = if a < b { ([a, b], 1) } else { ([b, a], -1) };
                let e = uses.entry(key).or_insert((0, 0));
                e.0 += 1;
                e.1 += dir;
            }
        }
        let mut r = EdgeReport {
            edges: uses.len(),
            ..EdgeReport::default()
        };
        for (key, (count, dir)) in uses {
            match count {
                1 => r.open.push(key),
                2 if dir != 0 => r.inconsistent.push(key),
                2 => {}
                _ => r.non_manifold.push(key),
            }
        }
        r
    }
}
