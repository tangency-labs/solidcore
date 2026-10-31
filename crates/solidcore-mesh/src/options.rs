//! What a tessellation aims for: [`TessellateOptions`], its presets and per-entity
//! overrides.

use std::collections::BTreeMap;
use std::fmt;
use std::sync::Arc;

use solidcore::{EdgeId, FaceId};

use crate::backend::FaceBackend;

/// A length given outright or as a fraction of the model size (the diagonal of the
/// box of the B-rep's vertices and edges, sampled).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Tolerance {
    /// Model units.
    Absolute(f64),
    /// A fraction of the bounding-box diagonal.
    Relative(f64),
}

impl Tolerance {
    /// The length in model units for a model of size `model_size`.
    #[must_use]
    pub fn resolve(self, model_size: f64) -> f64 {
        match self {
            Self::Absolute(t) => t,
            Self::Relative(r) => r * model_size,
        }
    }
}

/// Ready-made option sets, coarsest first.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Quality {
    /// Chord 1/500 of the model size, 30 degrees: previews and picking.
    Draft,
    /// Chord 1/2 000 of the model size, 15 degrees: display (the default).
    Normal,
    /// Chord 1/10 000 of the model size, 8 degrees: export and measurement.
    Fine,
}

/// How each edge is divided, beside the chord and angle tolerances it shares with the
/// faces.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EdgeOptions {
    /// Fewest segments of a curved or straight edge (at least 1).
    pub min_segments: u32,
    /// Most segments of an edge; splitting stops there (earlier segments in curve
    /// order are split first).
    pub max_segments: u32,
    /// Longest segment (model units), if bounded; [`TessellateOptions::max_edge_length`]
    /// also applies.
    pub max_segment_length: Option<f64>,
}

impl Default for EdgeOptions {
    /// One segment at least, 100 000 at most, no length bound.
    fn default() -> Self {
        Self {
            min_segments: 1,
            max_segments: 100_000,
            max_segment_length: None,
        }
    }
}

/// What the [`crate::Mesh`] carries beside positions and triangles.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(clippy::struct_excessive_bools)] // independent toggles
pub struct OutputOptions {
    /// Unit outward surface normals at every triangle corner ([`crate::Mesh::normals`]).
    pub normals: bool,
    /// The face's (u, v) at every triangle corner ([`crate::Mesh::uvs`]).
    pub uvs: bool,
    /// The B-rep face of every triangle ([`crate::Mesh::triangle_faces`]).
    pub face_ids: bool,
    /// The shared edge polylines ([`crate::Mesh::edges`]).
    pub edge_polylines: bool,
}

impl Default for OutputOptions {
    /// Everything on.
    fn default() -> Self {
        Self {
            normals: true,
            uvs: true,
            face_ids: true,
            edge_polylines: true,
        }
    }
}

/// Tolerances of one face that replace the global ones; `None` keeps the global value.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct FaceOverride {
    /// Chord tolerance (model units).
    pub chord_tolerance: Option<f64>,
    /// Angle tolerance (radians).
    pub angle_tolerance: Option<f64>,
    /// Longest interior triangle edge (model units).
    pub max_edge_length: Option<f64>,
    /// Most triangles of the face.
    pub max_triangles: Option<usize>,
}

/// Division settings of one edge that replace the derived ones; `None` keeps the
/// derived value.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct EdgeOverride {
    /// Chord tolerance (model units).
    pub chord_tolerance: Option<f64>,
    /// Angle tolerance (radians).
    pub angle_tolerance: Option<f64>,
    /// Longest segment (model units).
    pub max_segment_length: Option<f64>,
    /// Fewest segments.
    pub min_segments: Option<u32>,
    /// Most segments.
    pub max_segments: Option<u32>,
}

/// A per-face override hook: called once per face and level, `None` for no override.
pub type FaceHook = Arc<dyn Fn(FaceId) -> Option<FaceOverride> + Send + Sync>;
/// A per-edge override hook: called once per edge and level, `None` for no override.
pub type EdgeHook = Arc<dyn Fn(EdgeId) -> Option<EdgeOverride> + Send + Sync>;

/// What [`crate::tessellate`] aims for. Every option has a documented default
/// ([`Quality::Normal`]); the output is deterministic for any options.
#[derive(Clone)]
pub struct TessellateOptions {
    /// Largest distance between the mesh and the surface or curve it stands for.
    pub chord_tolerance: Tolerance,
    /// Largest angle (radians) a triangle's normal may turn from the surface normal at
    /// its centroid and interior corners, and that an edge segment's tangent or a
    /// face's normal along it may turn.
    pub angle_tolerance: f64,
    /// Longest triangle edge and edge segment (model units), if bounded.
    pub max_edge_length: Option<f64>,
    /// Shortest length refinement still splits (model units): no segment or triangle
    /// edge is split below twice this. `0` for no floor.
    pub min_edge_length: f64,
    /// Per-edge division settings.
    pub edge: EdgeOptions,
    /// Most triangles per face; refinement of a face stops near it.
    pub max_triangles_per_face: usize,
    /// Most refinement passes per face and level.
    pub max_refine_passes: u32,
    /// Which per-vertex and per-triangle data the mesh carries.
    pub output: OutputOptions,
    /// Measure the finished mesh against the true surfaces and curves
    /// ([`crate::MeshStats::max_deviation`]). Costs a few surface evaluations per
    /// triangle.
    pub measure_deviation: bool,
    /// Per-face tolerance overrides (they take precedence over [`Self::face_hook`]).
    pub face_overrides: BTreeMap<FaceId, FaceOverride>,
    /// Per-edge division overrides (they take precedence over [`Self::edge_hook`]).
    pub edge_overrides: BTreeMap<EdgeId, EdgeOverride>,
    /// Per-face override hook, consulted for faces without a map entry.
    pub face_hook: Option<FaceHook>,
    /// Per-edge override hook, consulted for edges without a map entry.
    pub edge_hook: Option<EdgeHook>,
    /// Supplies faces' interior grid lines (it may decline a face; the built-in sampled
    /// grid then). The triangulation and refinement stay the same. `None`: the built-in
    /// grid for every face.
    pub backend: Option<Arc<dyn FaceBackend>>,
}

impl fmt::Debug for TessellateOptions {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("TessellateOptions")
            .field("chord_tolerance", &self.chord_tolerance)
            .field("angle_tolerance", &self.angle_tolerance)
            .field("max_edge_length", &self.max_edge_length)
            .field("min_edge_length", &self.min_edge_length)
            .field("edge", &self.edge)
            .field("max_triangles_per_face", &self.max_triangles_per_face)
            .field("max_refine_passes", &self.max_refine_passes)
            .field("output", &self.output)
            .field("measure_deviation", &self.measure_deviation)
            .field("face_overrides", &self.face_overrides)
            .field("edge_overrides", &self.edge_overrides)
            .field("face_hook", &self.face_hook.as_ref().map(|_| "Fn"))
            .field("edge_hook", &self.edge_hook.as_ref().map(|_| "Fn"))
            .field("backend", &self.backend.as_ref().map(|b| b.name()))
            .finish()
    }
}

impl Default for TessellateOptions {
    /// [`Quality::Normal`].
    fn default() -> Self {
        Self::preset(Quality::Normal)
    }
}

impl TessellateOptions {
    /// The options of a quality preset: chord tolerance relative to the model size and
    /// angle tolerance as listed on [`Quality`]; no length bounds, at most 250 000
    /// triangles per face, 32 refinement passes, every output on, no deviation
    /// measurement, no overrides, no backend.
    #[must_use]
    pub fn preset(q: Quality) -> Self {
        let (chord, angle) = match q {
            Quality::Draft => (2e-3, 30.0_f64),
            Quality::Normal => (5e-4, 15.0),
            Quality::Fine => (1e-4, 8.0),
        };
        Self {
            chord_tolerance: Tolerance::Relative(chord),
            angle_tolerance: angle.to_radians(),
            max_edge_length: None,
            min_edge_length: 0.0,
            edge: EdgeOptions::default(),
            max_triangles_per_face: 250_000,
            max_refine_passes: 32,
            output: OutputOptions::default(),
            measure_deviation: false,
            face_overrides: BTreeMap::new(),
            edge_overrides: BTreeMap::new(),
            face_hook: None,
            edge_hook: None,
            backend: None,
        }
    }

    /// These options with an absolute chord tolerance.
    #[must_use]
    pub fn with_chord(mut self, chord: f64) -> Self {
        self.chord_tolerance = Tolerance::Absolute(chord);
        self
    }

    /// These options with a chord tolerance relative to the model size.
    #[must_use]
    pub fn with_relative_chord(mut self, fraction: f64) -> Self {
        self.chord_tolerance = Tolerance::Relative(fraction);
        self
    }

    /// The override of face `f`: its map entry, else the hook's answer.
    pub(crate) fn face_override(&self, f: FaceId) -> Option<FaceOverride> {
        self.face_overrides
            .get(&f)
            .copied()
            .or_else(|| self.face_hook.as_ref().and_then(|h| h(f)))
    }

    /// The override of edge `e`: its map entry, else the hook's answer.
    pub(crate) fn edge_override(&self, e: EdgeId) -> Option<EdgeOverride> {
        self.edge_overrides
            .get(&e)
            .copied()
            .or_else(|| self.edge_hook.as_ref().and_then(|h| h(e)))
    }
}

/// The tolerances one face is meshed to, every option resolved.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct FaceParams {
    pub chord: f64,
    pub angle: f64,
    pub max_len: Option<f64>,
    pub min_len: f64,
    pub max_triangles: usize,
    pub max_passes: u32,
}

/// The settings one edge is divided by, every option resolved.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct EdgeParams {
    pub chord: f64,
    pub angle: f64,
    pub max_len: Option<f64>,
    pub min_len: f64,
    pub min_segments: u32,
    pub max_segments: u32,
}

fn min_opt(a: Option<f64>, b: Option<f64>) -> Option<f64> {
    match (a, b) {
        (Some(x), Some(y)) => Some(x.min(y)),
        (x, None) => x,
        (None, y) => y,
    }
}

impl TessellateOptions {
    /// The resolved tolerances of face `f` for a model of size `model_size`.
    pub(crate) fn face_params(&self, f: FaceId, model_size: f64) -> FaceParams {
        let o = self.face_override(f).unwrap_or_default();
        FaceParams {
            chord: o
                .chord_tolerance
                .unwrap_or_else(|| self.chord_tolerance.resolve(model_size)),
            angle: o.angle_tolerance.unwrap_or(self.angle_tolerance),
            max_len: o.max_edge_length.or(self.max_edge_length),
            min_len: self.min_edge_length,
            max_triangles: o.max_triangles.unwrap_or(self.max_triangles_per_face),
            max_passes: self.max_refine_passes,
        }
    }

    /// The resolved settings of edge `e` along faces with resolved `faces` tolerances:
    /// the tightest of the global and the faces' chord, angle and length, then the
    /// edge's own override.
    pub(crate) fn edge_params(
        &self,
        e: EdgeId,
        model_size: f64,
        faces: &[FaceParams],
    ) -> EdgeParams {
        let mut chord = self.chord_tolerance.resolve(model_size);
        let mut angle = self.angle_tolerance;
        let mut max_len = min_opt(self.max_edge_length, self.edge.max_segment_length);
        for f in faces {
            chord = chord.min(f.chord);
            angle = angle.min(f.angle);
            max_len = min_opt(max_len, f.max_len);
        }
        let o = self.edge_override(e).unwrap_or_default();
        EdgeParams {
            chord: o.chord_tolerance.unwrap_or(chord),
            angle: o.angle_tolerance.unwrap_or(angle),
            max_len: o.max_segment_length.or(max_len),
            min_len: self.min_edge_length,
            min_segments: o.min_segments.unwrap_or(self.edge.min_segments).max(1),
            max_segments: o
                .max_segments
                .unwrap_or(self.edge.max_segments)
                .max(o.min_segments.unwrap_or(self.edge.min_segments).max(1)),
        }
    }
}
