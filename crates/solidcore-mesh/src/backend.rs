//! Face backends: other code that places a face's interior points.
//!
//! The B-rep layer divides the edges, maps each face's loops into UV (seams, poles and
//! all) and triangulates each face with a constrained Delaunay triangulation whose
//! constraints are the shared edge polylines, so the mesh is watertight by construction.
//! A [`FaceBackend`] may supply the face's interior grid lines (where the Steiner points
//! go) in place of the built-in sampled grid, and may measure the finished face. It never
//! triangulates: the triangulation stays the only one, and refinement still checks every
//! triangle against the surface.
//!
//! With the `nurbs-grid-mesh` feature, [`NurbsGridBackend`] places the lines of dense
//! B-spline faces per knot span with the nurbs-grid-mesh crate.

use std::fmt;

use solidcore::geom::extent2d::Extent2d;
use solidcore::kernel::Kernel;
use solidcore::math::{Point2, Point3};
use solidcore::surface::surface::Surface;
use solidcore::{FaceId, Result};

/// A face offered to a backend.
#[derive(Debug, Clone, Copy)]
pub struct BackendFace<'a> {
    /// The face.
    pub face: FaceId,
    /// Its surface.
    pub surface: &'a Surface,
    /// The UV box of its boundary loops.
    pub bbox: Extent2d,
    /// Chord tolerance (model units).
    pub chord: f64,
    /// Angle tolerance (radians); the built-in grid's lines stay where they are denser.
    pub angle: f64,
}

/// A finished face mesh, as a backend measures it: per-vertex UVs and positions and
/// triangles over them.
#[derive(Debug, Clone, Copy)]
pub struct FaceMeshView<'a> {
    /// Each vertex's (u, v).
    pub uvs: &'a [Point2],
    /// Each vertex's position.
    pub positions: &'a [Point3],
    /// Triangles over the vertices.
    pub triangles: &'a [[u32; 3]],
}

/// Supplies a face's interior grid lines and, optionally, measures its mesh.
pub trait FaceBackend: Send + Sync {
    /// A short name for reports.
    fn name(&self) -> &'static str;

    /// The face's interior grid lines `(us, vs)`, increasing and strictly inside
    /// `face.bbox`; `None` to decline the face (the built-in grid then). Whatever it
    /// returns is bounded like the built-in grid: sorted, clipped to the box and thinned
    /// to [`crate::MAX_GRID_LINES`] per direction and the face's triangle budget.
    ///
    /// # Errors
    ///
    /// What stopped the backend on a face it accepted.
    fn grid_lines(
        &self,
        k: &Kernel,
        face: &BackendFace<'_>,
    ) -> Option<Result<(Vec<f64>, Vec<f64>)>>;

    /// The largest distance of the mesh from the face's surface by the backend's own
    /// measure, `None` when it has none for this face.
    fn deviation(
        &self,
        _k: &Kernel,
        _face: &BackendFace<'_>,
        _mesh: &FaceMeshView<'_>,
    ) -> Option<f64> {
        None
    }
}

impl fmt::Debug for dyn FaceBackend {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "FaceBackend({})", self.name())
    }
}

#[cfg(feature = "nurbs-grid-mesh")]
pub use grid::NurbsGridBackend;

#[cfg(feature = "nurbs-grid-mesh")]
mod grid {
    use nurbs_grid_mesh::{deviation, grid_lines, Dir, MeshError, NurbsSurface, TrimmedMesh};
    use solidcore::kernel::Kernel;
    use solidcore::math::Point3;
    use solidcore::{Result, SolidError};

    use super::{BackendFace, FaceBackend, FaceMeshView};
    use crate::grid::MAX_GRID_LINES;
    use crate::nurbs::canonical_arrays;

    const OP: &str = "solidcore_mesh::NurbsGridBackend";

    /// Grid lines of dense B-spline faces from nurbs-grid-mesh (per knot span, from the
    /// second-derivative control-point bound), and its Newton-foot deviation: faces whose
    /// surface evaluates through a B-spline with at least [`NurbsGridBackend::min_poles`]
    /// control points.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct NurbsGridBackend {
        /// Fewest control points (u count times v count) of a face it takes.
        pub min_poles: usize,
        /// Most grid lines per knot span (fewer when the spans times this would pass
        /// [`crate::MAX_GRID_LINES`]).
        pub max_lines_per_span: usize,
        /// The lines are placed for this fraction of the chord tolerance (the
        /// triangulation's diagonals sag like the cells, so a margin saves refinement).
        pub chord_fraction: f64,
    }

    impl Default for NurbsGridBackend {
        /// Faces of 1 000 control points or more, 64 lines per span, 0.7 of the chord.
        fn default() -> Self {
            Self {
                min_poles: 1000,
                max_lines_per_span: 64,
                chord_fraction: 0.7,
            }
        }
    }

    fn err(e: &MeshError) -> SolidError {
        let reason = match e {
            MeshError::BadSurface(r) | MeshError::BadOptions(r) => r,
            MeshError::NoLoops => "no boundary loop",
            MeshError::LoopInsideCell => "a loop inside one grid cell",
            MeshError::OpenCell { .. } => "loop pieces in a grid cell did not close",
            MeshError::TooLarge => "too many vertices",
            MeshError::Internal => "nurbs-grid-mesh internal error",
        };
        SolidError::Unsupported {
            operation: OP,
            reason,
        }
    }

    impl NurbsGridBackend {
        /// The face's surface as nurbs-grid-mesh's, when it is dense enough and the
        /// arrays evaluate like the kernel's surface, with its knot span counts in u
        /// and v.
        fn surface(
            &self,
            k: &Kernel,
            face: &BackendFace<'_>,
        ) -> Option<(NurbsSurface, [usize; 2])> {
            let c = canonical_arrays(k, face.surface)?;
            let spans = |kv: &[f64]| {
                kv.windows(2)
                    .filter(|w| matches!(w, [a, b] if b > a))
                    .count()
            };
            let spans = [spans(&c.uk), spans(&c.vk)];
            if c.nu.saturating_mul(c.nv) < self.min_poles {
                return None;
            }
            let ns = NurbsSurface::new(
                c.p,
                c.q,
                c.uk,
                c.vk,
                c.nu,
                c.nv,
                c.poles,
                c.weights.as_deref(),
            )
            .ok()?;
            let b = face.bbox;
            for j in 0..=2 {
                for i in 0..=2 {
                    let uv = solidcore::math::Point2::new(
                        b.u_min() + (b.u_max() - b.u_min()) * f64::from(i) / 2.0,
                        b.v_min() + (b.v_max() - b.v_min()) * f64::from(j) / 2.0,
                    );
                    let mine = ns.evaluate(uv.x, uv.y).ok()?.point;
                    let theirs = face.surface.evaluate_point(k, uv).ok()?;
                    let d = Point3::new(mine[0], mine[1], mine[2]).distance(theirs);
                    if d > 1e-9 * (1.0 + theirs.length()) {
                        return None;
                    }
                }
            }
            Some((ns, spans))
        }
    }

    impl FaceBackend for NurbsGridBackend {
        fn name(&self) -> &'static str {
            "nurbs-grid-mesh"
        }

        fn grid_lines(
            &self,
            k: &Kernel,
            face: &BackendFace<'_>,
        ) -> Option<Result<(Vec<f64>, Vec<f64>)>> {
            let (ns, spans) = self.surface(k, face)?;
            let tol = self.chord_fraction * face.chord;
            let b = face.bbox;
            let lines = |dir: Dir, spans: usize, lo: f64, hi: f64| -> Result<Vec<f64>> {
                // lines per span times spans stays within MAX_GRID_LINES
                let per_span = self
                    .max_lines_per_span
                    .min(MAX_GRID_LINES / spans.max(1))
                    .max(1);
                let g = grid_lines(&ns, dir, tol, hi - lo, per_span).map_err(|e| err(&e))?;
                Ok(g.lines.into_iter().filter(|t| *t > lo && *t < hi).collect())
            };
            Some(
                lines(Dir::U, spans[0], b.u_min(), b.u_max())
                    .and_then(|us| Ok((us, lines(Dir::V, spans[1], b.v_min(), b.v_max())?))),
            )
        }

        fn deviation(
            &self,
            k: &Kernel,
            face: &BackendFace<'_>,
            mesh: &FaceMeshView<'_>,
        ) -> Option<f64> {
            let (ns, _) = self.surface(k, face)?;
            let m = TrimmedMesh {
                positions: mesh.positions.iter().map(|p| [p.x, p.y, p.z]).collect(),
                uvs: mesh.uvs.iter().map(|q| [q.x, q.y]).collect(),
                triangles: mesh.triangles.to_vec(),
                grid: [0, 0],
                capped_spans: 0,
                refined_lines: 0,
                dropped_degenerate: 0,
            };
            deviation(&ns, &m, face.chord).ok().map(|d| d.max)
        }
    }
}
