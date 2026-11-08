//! (not upstream) A body's soundness: its check at `SM_LEVEL_2`, and for a closed body,
//! whether it is a manifold solid of positive volume. The local operation measures its
//! input's [`standing`] (the check with the body's own owner naming, and the manifold
//! test) and its result's, and the result's [`volume`] only to test it is positive
//! ([`crate::topology::brep::local_operation::local_operation`]); the goldens soundness
//! probe (`goldens::soundness`) measures result bodies with [`soundness`] (the geometry's
//! recorded owners).

use std::collections::BTreeMap;

use crate::blend::fillet_section::interval::Ivl;
use crate::error::{EntityRef, Result};
use crate::geom::attribute::ObjectKey;
use crate::geom::object::ObjectRef;
use crate::heal::stitch::p2::ClosedGap;
use crate::kernel::Kernel;
use crate::surface::surface_cache::p3::definition_bound;
use crate::topology::Brep;

use super::assert_valid::{check_brep_with_keys, BrepCheck, CheckLevel};

/// The relative accuracy the volume is integrated at: its sign is a verdict on the body,
/// which a closed solid's whole volume holds away from zero.
const VOLUME_ACCURACY: f64 = 1e-6;

/// What [`soundness`] measured of one body.
#[derive(Debug, Clone)]
pub struct Soundness {
    /// The check at [`CheckLevel::Level2`], its geometry's owners taken as the geometry
    /// records them.
    pub check: BrepCheck,
    /// Whether the body is a manifold solid ([`crate::topology::brep::is_manifold_solid`]).
    pub manifold: bool,
    /// Its volume, when it is a manifold solid.
    pub volume: Option<f64>,
}

impl Soundness {
    /// Whether the check is clean and, when the body is a manifold solid, its volume is
    /// positive (a NaN volume is not).
    #[must_use]
    pub fn sound(&self) -> bool {
        self.check.is_valid() && (!self.manifold || self.volume.is_some_and(|v| v > 0.0))
    }
}

/// What [`standing`] measured of one body: its check and whether it is a manifold solid,
/// not its volume.
#[derive(Debug, Clone)]
pub struct Standing {
    /// The check at [`CheckLevel::Level2`], its geometry's owners the identities the
    /// body's naming scheme names its entities by.
    pub check: BrepCheck,
    /// Whether the body is a manifold solid ([`crate::topology::brep::is_manifold_solid`]).
    pub manifold: bool,
    /// What the body's data lets an entity hold: the largest tolerance it holds (its own,
    /// an edge's or a vertex's) and the enclosed gap of the vertices its tolerances join.
    /// An edge or vertex of a local operation's result that the operation kept is held to
    /// the input's tolerance of it; after the merge the kept ones are not told apart from
    /// the rebuilt, so every closed-form entity is held to this one
    /// ([`tolerance_beyond_bound`]).
    pub tolerance_cap: f64,
}

/// Measures `b`: [`check_brep_with_keys`] at [`CheckLevel::Level2`], each geometry's owner
/// the identity the geometry records (an edge's curve's, a face's surface's, an edge use's
/// UV trim curve's; an edge use without a trim the identity the body names it by, which the
/// trim the check makes is owned under), then the manifold test and, for a manifold solid,
/// the volume. What is checked is the body's geometry, not a caller's naming of it (a
/// runner names a body's entities by identities of its own): the owner tests cannot fail.
/// [`standing`] checks the owners.
///
/// The check makes the UV trim curves and caches its evaluations need on `b`: measure a
/// private copy ([`Brep::clone_standing_for_self`]) to leave a body as it was.
///
/// # Errors
///
/// Those of the check, the manifold test and the volume.
pub fn soundness(k: &mut Kernel, b: &mut Brep) -> Result<Soundness> {
    let mut owners = BTreeMap::new();
    for (id, e) in b.edges().iter() {
        if let Some(o) = e.curve.as_ref().and_then(crate::curve::curve::Curve::owner) {
            owners.insert(EntityRef::Edge(id), o);
        }
    }
    for (id, f) in b.faces().iter() {
        if let Some(o) = f
            .surface
            .as_ref()
            .and_then(crate::surface::surface::Surface::owner)
        {
            owners.insert(EntityRef::Face(id), o);
        }
    }
    for (id, u) in b.edgeuses().iter() {
        let o = match u.uv_trim_curve.as_ref() {
            Some(trim) => trim.owner(),
            None => Some(b.owner_identity(EntityRef::Edgeuse(id))),
        };
        if let Some(o) = o {
            owners.insert(EntityRef::Edgeuse(id), o);
        }
    }
    let keys = |e: EntityRef| {
        owners.get(&e).copied().unwrap_or(ObjectRef {
            type_code: 0,
            key: ObjectKey(0),
        })
    };
    let check = check_brep_with_keys(k, b, CheckLevel::Level2, &keys)?;
    let manifold = crate::topology::brep::p2::is_manifold_solid(b)?;
    let volume = if manifold { Some(volume(k, b)?) } else { None };
    Ok(Soundness {
        check,
        manifold,
        volume,
    })
}

/// Measures `b`'s standing: [`check_brep_with_keys`] at [`CheckLevel::Level2`], each
/// geometry's owner the identity `b`'s naming scheme names the entity by
/// ([`Brep::owner_identity`], as the operations that bind geometry set it), so a curve or
/// surface owned under another entity's identity fails the owner tests; then the manifold
/// test. Its volume is not measured ([`volume`]).
///
/// As for [`soundness`], the check edits what it measures: measure a private copy.
///
/// # Errors
///
/// Those of the check and the manifold test.
pub fn standing(k: &mut Kernel, b: &mut Brep) -> Result<Standing> {
    // the B-rep is borrowed mutably by the check: its naming and recorded keys are copied
    let naming = b.naming();
    let recorded = b.geometry_owner_keys().clone();
    let keys =
        |e: EntityRef| crate::topology::owner_identity_of(naming, recorded.get(&e).copied(), e);
    let check = check_brep_with_keys(k, b, CheckLevel::Level2, &keys)?;
    let manifold = crate::topology::brep::p2::is_manifold_solid(b)?;
    let tolerance_cap = input_cap(b)?;
    Ok(Standing {
        check,
        manifold,
        tolerance_cap,
    })
}

/// What the body's own data lets an entity hold: the largest of the body's tolerance and
/// its edges' and vertices', the enclosed gaps between its edges' ends and their vertices
/// ([`edge_end_gap`]), and the enclosed gap of the vertices its own tolerances already join.
///
/// Two vertices closer than the sum of their tolerances are one point by the body's own
/// certificate (`docs/validity-and-tolerances.md` section 5.1, invariant 6); merging them
/// grows the survivor's tolerance to cover that gap, which is the imported-data case of
/// section 8 (growth to an enclosed gap, the data's own), not a deviation of the
/// operation's construction. The gap is bounded above with its rounding
/// ([`crate::heal::stitch::bounds::points`]), the sum of the tolerances rounded up, and a
/// pair is joined when the bound is at most the sum. Vertices are swept along x so only
/// pairs within the largest sum are compared.
fn input_cap(b: &Brep) -> Result<f64> {
    let mut largest = b.tolerance();
    for (id, _) in b.edges().iter() {
        largest = largest.max(crate::topology::edge::tolerance(b, id)?);
        largest = largest.max(edge_end_gap(b, id)?.unwrap_or(0.0));
    }
    let mut vertices = Vec::new();
    let mut widest = 0.0_f64;
    for (id, v) in b.vertices().iter() {
        let tolerance = crate::topology::vertex::tolerance(b, id)?;
        largest = largest.max(tolerance);
        widest = widest.max(tolerance);
        vertices.push((v.point, tolerance));
    }
    vertices.sort_by(|a, c| a.0.x.total_cmp(&c.0.x));
    let reach = (Ivl::from(widest) + Ivl::from(widest)).hi;
    let mut joined = 0.0_f64;
    for (i, &(p, tp)) in vertices.iter().enumerate() {
        for &(q, tq) in vertices.iter().skip(i + 1) {
            if q.x - p.x > reach {
                break;
            }
            let gap = crate::heal::stitch::bounds::points(p, q)?;
            if gap <= (Ivl::from(tp) + Ivl::from(tq)).hi {
                joined = joined.max(gap);
            }
        }
    }
    Ok(largest.max(joined))
}

/// The enclosed gap between the ends of edge `e`'s curve and its vertices: the larger of
/// the distances at its two ends, the vertices taken to the ends they are nearer (an upper
/// bound with its rounding, [`crate::heal::stitch::bounds::point_at_curve`]). The gap is a
/// fact of the geometry the body holds, whatever the tolerances say
/// (`docs/validity-and-tolerances.md` invariant 6: an edge's ends are within the sum of its
/// and the vertex's tolerance): a tolerance raised to cover it is raised to a proved gap,
/// section 8's enclosed-gap case. `None` for an edge without a curve or a vertex, and for a
/// curve kind without an enclosure (only lines, ellipses and B-splines have one): its gap
/// is not certified, and never replaced by a sampled one.
///
/// # Errors
///
/// Those of the topology reads and the enclosure other than an unsupported curve kind.
fn edge_end_gap(b: &Brep, e: crate::topology::EdgeId) -> Result<Option<f64>> {
    use crate::heal::stitch::bounds::point_at_curve;
    let Some(curve) = b.edge(e)?.curve.as_ref() else {
        return Ok(None);
    };
    let (Some(start), Some(end)) = (
        crate::topology::edge::start_vertex(b, e)?,
        crate::topology::edge::end_vertex(b, e)?,
    ) else {
        return Ok(None);
    };
    let interval = crate::topology::edge::interval(b, e)?;
    let (low, high) = (interval.min(), interval.max());
    let (a, z) = (b.vertex(start)?.point, b.vertex(end)?.point);
    let gaps = (|| {
        let straight = point_at_curve(a, curve, low)?.max(point_at_curve(z, curve, high)?);
        let swapped = point_at_curve(z, curve, low)?.max(point_at_curve(a, curve, high)?);
        Ok(straight.min(swapped))
    })();
    match gaps {
        Ok(gap) => Ok(Some(gap)),
        Err(crate::SolidError::Unsupported { .. }) => Ok(None),
        Err(error) => Err(error),
    }
}

/// What [`tolerance_beyond_bound`] found: an edge or a vertex whose tolerance is above the
/// bound its certificate allows.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ToleranceBeyondBound {
    /// The edge or vertex.
    pub entity: EntityRef,
    /// Its tolerance.
    pub tolerance: f64,
    /// The bound.
    pub bound: f64,
}

/// The bound on the tolerance of edge `e`, `Some` only where a certificate gives one.
///
/// A constructed edge's deviation is delta + max `eps_i` (`docs/validity-and-tolerances.md`
/// sections 5.2 and 8), its tolerance that bound raised to the body's. Here delta is 0
/// where the stored curve is a closed form (`Curve::is_analytic`: a line, a circle or
/// ellipse, a parabola, a hyperbola, which are the intersection itself), and `eps_i` is
/// the certified bound of each face's definition ([`definition_bound`], 0 for a surface
/// that is its net). The bound is the largest of those, `cap` (the body's tolerance, or
/// the input's for an edge the operation kept) with the enclosed gap of the curve's ends to
/// its vertices ([`edge_end_gap`]) added, rounded up: a gap is covered by the sum of the
/// tolerances it lies between (`docs/validity-and-tolerances.md` section 5.2, Edge test 13
/// and invariant 6), so an entity may hold the body's tolerance and the gap, and no more. A traced curve (a B-spline from a
/// surface-surface intersection) has no proved delta (`docs/ssi.md` sections 6.1 and 9,
/// tangency-labs/solidcore#3211): `None`, the edge is not claimed certified.
fn edge_bound(b: &Brep, e: crate::topology::EdgeId, cap: f64) -> Result<Option<f64>> {
    let Some(curve) = b.edge(e)?.curve.as_ref() else {
        return Ok(None);
    };
    if !curve.is_analytic() {
        return Ok(None);
    }
    // its ends' enclosed gap to its vertices (a closed-form kind without an enclosure is
    // not certified)
    let Some(ends) = edge_end_gap(b, e)? else {
        return Ok(None);
    };
    let mut bound = (Ivl::from(cap) + Ivl::from(ends)).hi;
    for f in crate::topology::edge::faces(b, e, None)? {
        let Some(surface) = crate::topology::face::p2::surface(b, f)? else {
            return Ok(None);
        };
        bound = bound.max(definition_bound(surface));
    }
    Ok(Some(bound))
}

/// The first edge or vertex of `b` (in arena order) whose tolerance is above the bound a
/// certificate gives it, if any: an edge on a closed-form curve is bounded by
/// [`edge_bound`], a vertex all of whose edges are by the largest of theirs. Entities
/// whose bound is not certified (an edge on a traced curve, and the vertices it ends at)
/// are not judged here: that part of the decision waits for
/// tangency-labs/solidcore#3211.
///
/// `input` is the standing of the body the operation edited: its
/// [`Standing::tolerance_cap`] is what an entity the operation kept may hold. `gaps` are the
/// gaps the operation's stitch closed ([`ClosedGap`], each bound enclosed with its
/// rounding): a closure of a gap the cap already covers (`gap <= cap`) needs no growth, and
/// what it applied is the rounding of computing the tolerance, which the bound allows,
/// summed upward. A closure of a larger gap allows nothing: the growth it applied is beyond
/// the bound, and the entity that carries it is found here. Stitches inside the merge are
/// not reported, and their growth is held to the cap alone.
///
/// # Errors
///
/// Those of the topology reads.
pub fn tolerance_beyond_bound(
    b: &Brep,
    input: &Standing,
    gaps: &[ClosedGap],
) -> Result<Option<ToleranceBeyondBound>> {
    let base = input.tolerance_cap.max(b.tolerance());
    let mut allowance = Ivl::from(0.0);
    for gap in gaps.iter().filter(|g| g.gap <= base) {
        allowance = allowance + Ivl::from(gap.growth);
    }
    let cap = (Ivl::from(base) + allowance).hi;
    for (id, _) in b.edges().iter() {
        if let Some(bound) = edge_bound(b, id, cap)? {
            let tolerance = crate::topology::edge::tolerance(b, id)?;
            if tolerance > bound {
                return Ok(Some(ToleranceBeyondBound {
                    entity: EntityRef::Edge(id),
                    tolerance,
                    bound,
                }));
            }
        }
    }
    for (id, _) in b.vertices().iter() {
        let mut bound: Option<f64> = None;
        let mut certified = true;
        for e in crate::topology::vertex::edges(b, id, None)? {
            match edge_bound(b, e, cap)? {
                Some(edge) => bound = Some(bound.map_or(edge, |b: f64| b.max(edge))),
                None => certified = false,
            }
        }
        if let (true, Some(bound)) = (certified, bound) {
            let tolerance = crate::topology::vertex::tolerance(b, id)?;
            if tolerance > bound {
                return Ok(Some(ToleranceBeyondBound {
                    entity: EntityRef::Vertex(id),
                    tolerance,
                    bound,
                }));
            }
        }
    }
    Ok(None)
}

/// The volume of the manifold solid `b`, integrated at [`VOLUME_ACCURACY`].
///
/// # Errors
///
/// Those of [`crate::mass::integrator::properties::brep_volume_default`].
pub fn volume(k: &mut Kernel, b: &mut Brep) -> Result<f64> {
    crate::mass::integrator::properties::brep_volume_default(k, b, VOLUME_ACCURACY)
}
