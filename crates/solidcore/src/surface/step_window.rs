//! (not upstream) A face's STEP parameter window from its boundary: the box in the
//! surface's STEP (u, v) parameters that the face's bounding edges span on it.
//!
//! Upstream: none (solidcore's STEP reader, #2578). The idea is SolidCraft's reader
//! (`crates/step/src/topology/swept.rs` `span`, `range`, `trimmed_range`,
//! `widest_gap`, `common_angle`, `crates/step/src/topology.rs` `Builder::extent` and
//! `face_extent`, and `crates/geom/src/brep/form.rs` `extent`, at a5b71ca2; MIT OR
//! Apache-2.0).
//!
//! | Rust                 | SolidCraft                                                  |
//! |----------------------|-------------------------------------------------------------|
//! | [`BoundaryEdge`]     | `form::extent`'s `(&Curve, [f64; 2], bool)` loop entries    |
//!
//! Deviations: not a port.

use crate::curve::curve::Curve;
use crate::geom::extent1d::Extent1d;

/// One edge of a face's boundary loop, as the loop uses it: the edge's curve, its
/// STEP interval on that curve and the direction the loop runs along it.
///
/// A face's loops are lists of these in loop order, each edge running on from where
/// the one before it ends. The loops are oriented as they bound the face: the face
/// lies on their left about the surface's STEP normal `Su x Sv` (a file's
/// `FACE_BOUND.orientation` applied, and every loop reversed when the face's
/// `same_sense` is false).
#[derive(Debug, Clone, Copy)]
pub struct BoundaryEdge<'a> {
    /// The edge's curve (for a `CIRCLE` or an `ELLIPSE` the full conic), whose STEP
    /// parameters [`Self::interval`] is in (those of [`Curve::evaluate_step_point`]).
    pub curve: &'a Curve,
    /// The edge's STEP interval on [`Self::curve`], lower parameter first: the
    /// [`crate::curve::step_edge::StepEdgeInterval::interval`] of
    /// [`crate::curve::step_edge::edge_interval`].
    pub interval: Extent1d,
    /// The loop runs along the curve's parameter, from the interval's minimum to its
    /// maximum; when false, from the maximum down to the minimum.
    pub forward: bool,
}
