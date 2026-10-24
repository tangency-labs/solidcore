//! Measuring a mesh against the true geometry: each triangle sampled at seven points,
//! each point's distance to its Newton foot on the face's surface (a surface point, so
//! the distance bounds the true one from above), and each edge segment against its
//! curve. (The sampling follows the nurbs-grid-mesh design note's `deviation`.)

use solidcore::geom::extent2d::Extent2d;
use solidcore::kernel::Kernel;
use solidcore::math::{Point2, Point3};
use solidcore::surface::surface::Surface;
use solidcore::topology::Brep;
use solidcore::Result;

use crate::face::FaceMesh;
use crate::geom::{curve_eval, foot, invalid, segment_distance, surface_normal};
use crate::mesh::{EdgePolyline, Mesh};

/// Barycentric weights of the seven samples: centroid, edge midpoints, and three points
/// between the centroid and the corners.
const SAMPLES: [[f64; 3]; 7] = [
    [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
    [0.5, 0.5, 0.0],
    [0.0, 0.5, 0.5],
    [0.5, 0.0, 0.5],
    [2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0],
    [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0],
    [1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0],
];

/// The measured deviation of a whole mesh ([`deviation`]).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct DeviationReport {
    /// Largest distance over every sample of every triangle and edge segment.
    pub max: f64,
    /// Largest distance of a triangle sample from its surface.
    pub face_max: f64,
    /// Largest distance of an edge segment from its curve.
    pub edge_max: f64,
    /// Mean of the triangles' largest sample distances.
    pub mean: f64,
    /// Triangles whose largest sample distance exceeds `tolerance`.
    pub over: usize,
    /// The triangle with the largest distance.
    pub worst: Option<usize>,
    /// Triangles wound against their face's outward side: the facet normal opposes the
    /// outward surface normal at the foot of the triangle's centroid.
    pub flipped: usize,
    /// The flipped triangles' indices, in order.
    pub flipped_triangles: Vec<usize>,
}

/// The largest sample distance of one triangle (UVs `uv`, corners `p`) from surface `s`.
pub(crate) fn triangle_deviation(
    k: &Kernel,
    s: &Surface,
    domain: &Extent2d,
    uv: [Point2; 3],
    p: [Point3; 3],
) -> f64 {
    let mut worst = 0.0_f64;
    for w in SAMPLES {
        let q = p[0] * w[0] + p[1] * w[1] + p[2] * w[2];
        let start = Point2::new(
            uv[0].x * w[0] + uv[1].x * w[1] + uv[2].x * w[2],
            uv[0].y * w[0] + uv[1].y * w[1] + uv[2].y * w[2],
        );
        if let (_, Some(d)) = foot(k, s, start, q, domain) {
            worst = worst.max(d);
        }
    }
    worst
}

/// The largest triangle deviation of a face mesh and how many triangles exceed `chord`
/// plus their corners' largest gap from the surface.
pub(crate) fn face_mesh_deviation(
    k: &Kernel,
    s: &Surface,
    domain: &Extent2d,
    m: &FaceMesh,
    chord: f64,
) -> (f64, usize) {
    let mut max = 0.0_f64;
    let mut over = 0;
    for t in &m.triangles {
        let get = |i: u32| Some((*m.uvs.get(i as usize)?, *m.positions.get(i as usize)?));
        let (Some(a), Some(b), Some(c)) = (get(t[0]), get(t[1]), get(t[2])) else {
            continue;
        };
        let d = triangle_deviation(k, s, domain, [a.0, b.0, c.0], [a.1, b.1, c.1]);
        let gap = t
            .iter()
            .filter_map(|i| m.gaps.get(*i as usize))
            .fold(0.0_f64, |x, y| x.max(*y));
        if d > chord + gap {
            over += 1;
        }
        max = max.max(d);
    }
    (max, over)
}

/// The largest distance of an edge polyline's segments from the edge curve (sampled
/// at each segment's quarter points and midpoint).
pub(crate) fn edge_deviation(
    k: &Kernel,
    b: &Brep,
    poly: &EdgePolyline,
    positions: &[Point3],
) -> Result<f64> {
    if poly.degenerate {
        return Ok(0.0);
    }
    let rec = b.edge(poly.edge)?;
    let curve = rec
        .curve
        .as_ref()
        .ok_or_else(|| invalid("edge without curve"))?;
    let mut worst = 0.0_f64;
    for (ts, vs) in poly.params.windows(2).zip(poly.vertices.windows(2)) {
        let ([t0, t1], [v0, v1]) = (ts, vs) else {
            continue;
        };
        let (Some(a), Some(c)) = (positions.get(*v0 as usize), positions.get(*v1 as usize)) else {
            continue;
        };
        for f in [0.25, 0.5, 0.75] {
            let (q, _) = curve_eval(k, curve, t0 + (t1 - t0) * f, false)?;
            worst = worst.max(segment_distance(q, *a, *c));
        }
    }
    Ok(worst)
}

/// Measures `mesh` against `b`'s surfaces and curves: its triangles need
/// [`Mesh::uvs`] and [`Mesh::triangle_faces`], its edges [`Mesh::edges`]. `tolerance`
/// counts the triangles over it.
///
/// # Errors
///
/// [`solidcore::SolidError::InvalidInput`] for a mesh without UVs or face ids;
/// [`solidcore::SolidError::NotFound`] for faces or edges not in `b`.
pub fn deviation(k: &Kernel, b: &Brep, mesh: &Mesh, tolerance: f64) -> Result<DeviationReport> {
    if mesh.uvs.len() != mesh.triangles.len() || mesh.triangle_faces.len() != mesh.triangles.len() {
        return Err(invalid("deviation needs the mesh's UVs and face ids"));
    }
    let mut r = DeviationReport::default();
    let mut sum = 0.0;
    for (i, ((t, uv), f)) in mesh
        .triangles
        .iter()
        .zip(&mesh.uvs)
        .zip(&mesh.triangle_faces)
        .enumerate()
    {
        let s = b
            .face(*f)?
            .surface
            .as_ref()
            .ok_or_else(|| invalid("face without surface"))?;
        let Some(p) = mesh.triangle_points(*t) else {
            continue;
        };
        let domain = s.natural_uv_domain();
        let d = triangle_deviation(k, s, &domain, *uv, p);
        sum += d;
        let facet = (p[1] - p[0]).cross(p[2] - p[0]);
        let centroid = (p[0] + p[1] + p[2]) / 3.0;
        let start = Point2::new(
            (uv[0].x + uv[1].x + uv[2].x) / 3.0,
            (uv[0].y + uv[1].y + uv[2].y) / 3.0,
        );
        let (at, _) = foot(k, s, start, centroid, &domain);
        if let Some(n) = surface_normal(k, s, at) {
            if n.dot(facet) * crate::face::outward_sense(b, *f)? < 0.0 {
                r.flipped += 1;
                r.flipped_triangles.push(i);
            }
        }
        if d > tolerance {
            r.over += 1;
        }
        if d > r.face_max || r.worst.is_none() {
            r.face_max = r.face_max.max(d);
            r.worst = Some(i);
        }
    }
    if !mesh.triangles.is_empty() {
        r.mean = sum / mesh.triangles.len() as f64;
    }
    for e in &mesh.edges {
        r.edge_max = r.edge_max.max(edge_deviation(k, b, e, &mesh.positions)?);
    }
    r.max = r.face_max.max(r.edge_max);
    Ok(r)
}
