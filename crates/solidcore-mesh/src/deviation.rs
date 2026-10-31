//! Measuring a mesh against the true geometry: each triangle sampled at seven points,
//! each point's distance to its Newton foot on the face's surface (a surface point, so
//! the distance bounds the true one from above), and each edge segment against its
//! curve. (The sampling follows the nurbs-grid-mesh design note's `deviation`.)

use std::collections::BTreeMap;

use solidcore::geom::extent2d::Extent2d;
use solidcore::kernel::Kernel;
use solidcore::math::{Point2, Point3};
use solidcore::surface::surface::Surface;
use solidcore::topology::Brep;
use solidcore::{FaceId, Result};

use crate::face::FaceMesh;
use crate::geom::{closes, curve_eval, foot, invalid, segment_distance, surface_normal};
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

/// `uv` with its second and third corners moved by whole periods (`period`, of a
/// surface closed in u or v) to within half a period of the first, when that makes the
/// triangle smaller in UV: a triangle of a face with no seam edge may straddle the
/// period, and blending its UVs as given would start the samples on the far side.
pub(crate) fn unwrap_triangle(uv: [Point2; 3], period: (Option<f64>, Option<f64>)) -> [Point2; 3] {
    let [a, b, c] = uv;
    let near = |x: f64, to: f64, p: Option<f64>| match p {
        Some(p) if p > 0.0 && p.is_finite() => x - ((x - to) / p).round() * p,
        _ => x,
    };
    let span = |xs: [f64; 3]| {
        let lo = xs[0].min(xs[1]).min(xs[2]);
        let hi = xs[0].max(xs[1]).max(xs[2]);
        hi - lo
    };
    let mut out = uv;
    for (axis, p) in [(0, period.0), (1, period.1)] {
        let get = |q: Point2| if axis == 0 { q.x } else { q.y };
        let given = [get(a), get(b), get(c)];
        let moved = [
            given[0],
            near(given[1], given[0], p),
            near(given[2], given[0], p),
        ];
        if span(moved) < span(given) {
            for (q, x) in out.iter_mut().zip(moved) {
                if axis == 0 {
                    q.x = x;
                } else {
                    q.y = x;
                }
            }
        }
    }
    out
}

/// The largest sample distance of one triangle (UVs `uv`, corners `p`) from surface `s`
/// (closed with `period` in u or v, see [`unwrap_triangle`]).
pub(crate) fn triangle_deviation(
    k: &Kernel,
    s: &Surface,
    domain: &Extent2d,
    period: (Option<f64>, Option<f64>),
    uv: [Point2; 3],
    p: [Point3; 3],
) -> f64 {
    let uv = unwrap_triangle(uv, period);
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
    period: (Option<f64>, Option<f64>),
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
        let d = triangle_deviation(k, s, domain, period, [a.0, b.0, c.0], [a.1, b.1, c.1]);
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
    // the surfaces' periods (a closed surface's face may have no seam edge), by face
    let mut lo = Point3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut hi = Point3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for q in &mesh.positions {
        lo = Point3::new(lo.x.min(q.x), lo.y.min(q.y), lo.z.min(q.z));
        hi = Point3::new(hi.x.max(q.x), hi.y.max(q.y), hi.z.max(q.z));
    }
    let size = if mesh.positions.is_empty() {
        0.0
    } else {
        lo.distance(hi)
    };
    let tol = 1e-9 * size.max(1e-300);
    let mut periods: BTreeMap<FaceId, (Option<f64>, Option<f64>)> = BTreeMap::new();
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
        // the domain upstream's GetNaturalUVDomain gives: a lazy analytic face reads the
        // NURB its MakeNurb builds (made once, kept on the surface)
        let domain = s.natural_uv_domain_with_nurb(k)?;
        let period = *periods.entry(*f).or_insert_with(|| {
            let along = |u: bool| {
                let width = if u {
                    domain.u_max() - domain.u_min()
                } else {
                    domain.v_max() - domain.v_min()
                };
                closes(k, s, &domain, u, tol).then_some(width)
            };
            (along(true), along(false))
        });
        let uv = &unwrap_triangle(*uv, period);
        let d = triangle_deviation(k, s, &domain, period, *uv, p);
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_triangle_straddling_the_period_is_unwrapped() {
        let tau = core::f64::consts::TAU;
        let p = (Some(tau), None);
        let t = [
            Point2::new(0.05, 0.2),
            Point2::new(tau - 0.05, 0.2),
            Point2::new(0.0, 0.4),
        ];
        let u = unwrap_triangle(t, p);
        assert!((u[1].x - (-0.05)).abs() < 1e-12, "{u:?}");
        assert_eq!(u[2].x, 0.0);
        // a triangle that is already the smallest stays as it is
        let wide = [
            Point2::new(0.0, 0.0),
            Point2::new(tau / 3.0, 0.0),
            Point2::new(2.0 * tau / 3.0, 0.0),
        ];
        assert_eq!(unwrap_triangle(wide, p), wide);
        // no period: unchanged
        assert_eq!(unwrap_triangle(t, (None, None)), t);
    }
}
