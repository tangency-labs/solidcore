//! How far past its face an offset face's surface must reach for the offsets of its
//! neighbours to meet its own: the margin the offset executive extends each surface by before
//! it intersects the offsets (not upstream: upstream extends every surface by five offset
//! distances).
//!
//! At a point `p` of an edge between a face `F` (unit normal `nF`, binormal `bF` into `F`,
//! offset `dF`) and its neighbour `G` (unit normal `nG`, offset `dG`), `c = nF·nG`, the
//! offsets' tangent planes `(x − p)·nF = dF` and `(x − p)·nG = dG` meet on the line through
//! `p + a nF + β nG`, `a + β c = dF`, `a c + β = dG`. Its foot on `F`'s tangent plane lies
//! past the edge (along the outward binormal `m = −bF`) by the edge's signed reach
//!
//! ```text
//! r = β (nG·m),   β = (dG − c dF) / (1 − c²)   ( = d / (1 + c) where dF = dG = d ),
//! ```
//!
//! `d·cot(α / 2)` for equal offsets, `α` the sector angle of the solid between the faces:
//! `d` across a box's edge, `0` across a tangent edge, negative where the offsets meet inside
//! `F`. At a vertex where `F`'s edge uses meet with reaches `r1`, `r2` and outward binormals
//! `m1`, `m2`, the offset corner's foot `x` on the tangent plane solves `x·m1 = r1`,
//! `x·m2 = r2`, at `|x|` from the vertex.
//!
//! Each foot is within `|r|` (`|x|` at a corner) of a point of the face, so within that of the
//! face's surface's box, however far the meeting runs along the face: an inward meeting that
//! leaves the face across another of its edges (a hole's offset growing past a side) is held
//! as well as an outward one. The reach of a face is the largest of these, sampled along each
//! edge (its sector angle varies) and at each corner. Faces whose offsets never meet, a knife
//! edge's (`c = −1` for equal offsets) or tangent faces offset apart (`c = 1`, `dF ≠ dG`),
//! add nothing: no extension makes them meet.
//!
//! An extension meets its distance at the middle of each side, and a reach is needed at every
//! station. [`shortest_across`] bounds the length across a B-spline's sides from below, for
//! the executive to ask again where the rows narrow along a side (#2908; the measurements and
//! the choice of asking again are `MarkusReichart`'s on #2926). Analytic faces are left as
//! they were.

use crate::error::offset::OffsetFailure;
use crate::error::Failure;
use crate::error::{Result, SolidError};
use crate::geom::extent2d::Extent2d;
use crate::geom::vector3d::{sqrt, Vec3};
use crate::kernel::Kernel;
use crate::math::{Point2, EFF_ZERO};
use crate::nurbs::b_spline_surface::p1::SurfParamType;
use crate::nurbs::b_spline_surface::p2::IsoCurveOptions;
use crate::nurbs::b_spline_surface_n_lib::p2::OrientType;
use crate::nurbs::knots::Side;
use crate::surface::extend::{ExtendedSides, SideExtension};
use crate::surface::implicit::{self, AnalyticSupport, Implicit};
use crate::surface::side::BoxSide;
use crate::surface::surface::{Surface, SurfaceKind};
use crate::topology::edgeuse::{p1 as eu, p3 as eu3};
use crate::topology::{edge, loopuse, Brep, EdgeuseId, FaceId, FaceuseId};

const OP: &str = "offset::offset_executive::reach::extension_reach";
const LOCAL_OP: &str = "offset::offset_executive::reach::local_extension_reach";

/// The number of intervals each edge is sampled at for its sector angle.
const EDGE_SAMPLES: u32 = 8;

/// One end of an edge use, as the corner beside it reads it.
#[derive(Debug, Clone, Copy)]
struct End {
    /// The signed reach past the edge.
    reach: f64,
    /// The outward binormal (away from the face, perpendicular to the edge).
    outward: Vec3,
}

/// The signed reach past the edge of a face offset by `own` whose neighbour is offset by
/// `other`, from their unit normals and the face's outward binormal (module docs); zero where
/// the offsets do not meet: faces folded back on each other (equal offsets) or tangent to
/// each other (different offsets: their offsets run apart), each within `tangent_sine`, the
/// sine of the tangency tolerance.
fn signed_reach(
    n_own: Vec3,
    n_other: Vec3,
    outward: Vec3,
    own: f64,
    other: f64,
    tangent_sine: f64,
) -> f64 {
    let c = n_own.dot(n_other);
    let sine_squared = (tangent_sine * tangent_sine).max(EFF_ZERO);
    #[allow(clippy::float_cmp)] // the same offset, not a computed value
    let beta = if own == other {
        // 1 + c = 1 - cos(θ) ≈ θ² / 2 for normals θ off opposite
        let denominator = 1.0 + c;
        if denominator <= 0.5 * sine_squared {
            return 0.0;
        }
        own / denominator
    } else {
        // 1 - c² = sin²(θ) for normals θ apart
        let denominator = 1.0 - c * c;
        if denominator <= sine_squared {
            return 0.0;
        }
        (other - c * own) / denominator
    };
    beta * n_other.dot(outward)
}

/// The signed reach past `eu`'s edge at `parameter` (module docs), with `eu`'s outward
/// binormal there, and whether the edge has a neighbour (zero reach without).
fn edge_reach(
    k: &Kernel,
    b: &mut Brep,
    eu: EdgeuseId,
    parameter: f64,
    offset: f64,
    unmoved: &[FaceId],
    tangent_sine: f64,
) -> Result<(End, bool)> {
    let options = eu3::BinormalOptions::default();
    let own = eu3::evaluate_binormal(k, b, eu, parameter, false, &options)?;
    let outward = -own.binormal;
    // A seam's radial is on the face itself (its offset continues across), a lamina
    // edge's is its mate on the face's other side (its offset is closed by a side wall).
    let own_face = eu::face(b, eu)?;
    let neighbour = match eu::radial(b, eu)? {
        Some(r) if r != eu && eu::face(b, r)? != own_face => Some(r),
        _ => None,
    };
    let Some(radial) = neighbour else {
        let end = End {
            reach: 0.0,
            outward,
        };
        return Ok((end, false));
    };
    let other_offset = match eu::face(b, radial)? {
        Some(f) if unmoved.contains(&f) => 0.0,
        _ => offset,
    };
    let other = eu3::evaluate_binormal(k, b, radial, parameter, false, &options)?;
    let n_own = own.face_normal.unitize()?;
    let n_other = other.face_normal.unitize()?;
    let end = End {
        reach: signed_reach(n_own, n_other, outward, offset, other_offset, tangent_sine),
        outward,
    };
    Ok((end, true))
}

/// The distance from a vertex of the foot of its offset corner on the face's tangent plane
/// (module docs), from the reaches and outward binormals of the edge uses meeting there.
/// Where the edges run on tangent to each other (within `tangent_sine`, the sine of the
/// tangency tolerance) the corner is no corner: the offset edges run on into each other, and
/// the larger reach stands (the tangent lines' crossing, far off for nearly parallel lines,
/// is not where curved edges' offsets meet).
fn corner_reach(into: End, out: End, tangent_sine: f64) -> f64 {
    let g = into.outward.dot(out.outward);
    let det = 1.0 - g * g;
    if det <= (tangent_sine * tangent_sine).max(EFF_ZERO) {
        return into.reach.abs().max(out.reach.abs());
    }
    let a = (into.reach - g * out.reach) / det;
    let c = (out.reach - g * into.reach) / det;
    sqrt((a * a + c * c + 2.0 * a * c * g).max(0.0))
}

/// How far past its face the surface of `faceuse` must be extended for the offsets of its
/// neighbours to meet its own, the face offset by `offset` (signed, along the face use's
/// normal) and each neighbour by the same, except the `unmoved` faces, which stay where they
/// are (module docs): the largest distance from the face of the feet of the meetings, sampled
/// along each edge at its sector angle and at each corner (edges meeting within
/// `tangency_tol_radians` of tangent make no corner); `None` for a face with no
/// neighbour (every edge a lamina edge or a seam), whose offset meets no other.
///
/// # Errors
///
/// [`SolidError::InvalidInput`] for an offset that is not finite or a vanishing face normal;
/// those of [`eu3::evaluate_binormal`]; [`SolidError::NotFound`] for a dangling id.
pub fn extension_reach(
    k: &Kernel,
    b: &mut Brep,
    faceuse: FaceuseId,
    offset: f64,
    unmoved: &[FaceId],
    tangency_tol_radians: f64,
) -> Result<Option<f64>> {
    if !offset.is_finite() {
        return Err(SolidError::invalid_input(OP, "the offset is not finite"));
    }
    let tangent_sine = tangency_tol_radians.abs().min(crate::math::PI / 2.0).sin();
    let mut reach = 0.0_f64;
    let mut met = false;
    for lu in eu::faceuse_loopuses(b, faceuse, false, None)? {
        let mut ends: Vec<(End, End)> = Vec::new();
        for edgeuse in loopuse::edgeuses(b, lu, None, None)? {
            let Some(e) = eu::edge(b, edgeuse)? else {
                continue;
            };
            let interval = edge::interval(b, e)?;
            let (lo, hi) = (interval.min(), interval.max());
            let mut first = None;
            let mut last = None;
            for i in 0..=EDGE_SAMPLES {
                let parameter = lo + (hi - lo) * f64::from(i) / f64::from(EDGE_SAMPLES);
                let (end, neighbour) =
                    edge_reach(k, b, edgeuse, parameter, offset, unmoved, tangent_sine)?;
                met |= neighbour;
                reach = reach.max(end.reach.abs());
                if i == 0 {
                    first = Some(end);
                }
                last = Some(end);
            }
            let (Some(first), Some(last)) = (first, last) else {
                continue;
            };
            // The edge use runs from its start to its end vertex along the loop.
            ends.push(match eu::orientation(b, edgeuse)? {
                OrientType::Opposite => (last, first),
                _ => (first, last),
            });
        }
        // Each edge use's end meets the next one's start (the loop closes).
        let next = ends.iter().cycle().skip(1);
        for (&(_, into), &(out, _)) in ends.iter().zip(next) {
            reach = reach.max(corner_reach(into, out, tangent_sine));
        }
    }
    Ok(met.then_some(reach))
}

/// A lower bound of the distance from each extended side of `base`'s B-spline to the same
/// side of `extended`, over the side's stations on `base`'s box: `None` where `base` is no
/// B-spline or the offset of one, a side's two rows differ in weights, or the bound is not positive. Each row is a
/// curve in the same parameter and basis, with the same weights, so the difference at a
/// station is a convex combination of the differences of their control points, and its length
/// is at least the least of those along their mean direction. Past the box's corners the
/// rows go on narrowing, to nothing at a ruled face's apex; the margin's one offset past the
/// reach is what stands for them.
///
/// # Errors
///
/// Those of the domains and the iso curves.
pub(super) fn shortest_across(
    k: &Kernel,
    base: &Surface,
    extended: &ExtendedSides,
) -> Result<Option<f64>> {
    // an offset surface is the offset of its base extended: the rows are its base's
    let (near_surface, far_surface) = match (base.kind(), extended.surface.kind()) {
        (SurfaceKind::OffsetSurface(a), SurfaceKind::OffsetSurface(b)) => {
            (a.base_surface(), b.base_surface())
        }
        _ => (base, &extended.surface),
    };
    let (SurfaceKind::BSpline(_), SurfaceKind::BSpline(net)) =
        (near_surface.kind(), far_surface.kind())
    else {
        return Ok(None);
    };
    let from = near_surface.natural_uv_domain_with_nurb(k)?;
    let to = far_surface.natural_uv_domain_with_nurb(k)?;
    let mut shortest = f64::INFINITY;
    for (&side, done) in BoxSide::ALL.iter().zip(&extended.sides) {
        if !matches!(done, SideExtension::Extended(_)) {
            continue;
        }
        let along = if side.is_u() {
            SurfParamType::U
        } else {
            SurfParamType::V
        };
        // the side's stations on the source's box, the row wherever it lies
        let box_of_row = if side.is_u() {
            Extent2d::from_points(
                Point2::new(to.min().x, from.min().y),
                Point2::new(to.max().x, from.max().y),
            )
        } else {
            Extent2d::from_points(
                Point2::new(from.min().x, to.min().y),
                Point2::new(from.max().x, to.max().y),
            )
        };
        let row = |at: f64| {
            let options = IsoCurveOptions {
                domain: Some(box_of_row),
                uv_iso_curve: false,
            };
            net.create_iso_parametric_curve(k, along, at, options)?
                .curve
                .canonical()
        };
        let (near, far) = (row(side.at(&from))?, row(side.at(&to))?);
        if near.control_points.len() != far.control_points.len() || near.weights != far.weights {
            return Ok(None);
        }
        let mut mean = Vec3::ZERO;
        for (a, b) in near.control_points.iter().zip(&far.control_points) {
            mean += *b - *a;
        }
        let Ok(direction) = mean.unitize() else {
            return Ok(None);
        };
        for (a, b) in near.control_points.iter().zip(&far.control_points) {
            shortest = shortest.min((*b - *a).dot(direction));
        }
    }
    Ok((shortest.is_finite() && shortest > 0.0).then_some(shortest))
}

/// A face's surface near a point, as a plane: a point on it and its unit normal.
#[derive(Debug, Clone, Copy)]
struct Tangent {
    point: Vec3,
    normal: Vec3,
    /// Whether the face moves to this surface (else it is the face's own, at the point).
    moved: bool,
}

/// The tangent plane of `surface` at its point nearest `p`: an infinite plane's own (no
/// solve), any other surface's at its closest point.
fn tangent_near(k: &Kernel, surface: &Surface, p: Vec3) -> Result<Tangent> {
    if let SurfaceKind::Plane(plane) = surface.kind() {
        let normal = plane.position().z_axis().unitize()?;
        let origin = plane.position().origin();
        return Ok(Tangent {
            point: p - normal * (p - origin).dot(normal),
            normal,
            moved: true,
        });
    }
    let closest = crate::query::api_queries::surface_closest_point(k, surface, p)?;
    let uv = closest
        .uvs
        .first()
        .copied()
        .ok_or_else(|| SolidError::sm_err(LOCAL_OP))?;
    let normal = surface
        .evaluate_normal(k, uv, Side::Left, Side::Left)?
        .unitize()?;
    Ok(Tangent {
        point: closest.closest.point,
        normal,
        moved: true,
    })
}

/// The point nearest `p` where the planes `a` and `b` meet: on their line of meeting, the
/// point `p + α a.n + β b.n` (`α + β c = h_a`, `α c + β = h_b`, `c = a.n·b.n`, `h` each
/// plane's height over `p`). Planes within `tangent_sine` of parallel meet only where they
/// are one plane (within `tolerance` apart; then at `p`'s foot on it), else not at all.
fn meeting_near(
    p: Vec3,
    a: Tangent,
    b: Tangent,
    tangent_sine: f64,
    tolerance: f64,
) -> Option<Vec3> {
    let c = a.normal.dot(b.normal);
    let h_a = (a.point - p).dot(a.normal);
    let h_b = (b.point - p).dot(b.normal);
    let det = 1.0 - c * c;
    if det <= (tangent_sine * tangent_sine).max(EFF_ZERO) {
        // one plane where the other's point lies on it, else apart
        return ((b.point - a.point).dot(a.normal).abs() <= tolerance).then(|| p + a.normal * h_a);
    }
    let alpha = (h_a - c * h_b) / det;
    let beta = (h_b - c * h_a) / det;
    Some(p + a.normal * alpha + b.normal * beta)
}

/// The point where the planes `a`, `b` and `c` meet, `None` where they meet in no one point
/// (two of them within `tangent_sine` of parallel: their triple product is below it).
fn corner_of(a: Tangent, b: Tangent, c: Tangent, tangent_sine: f64) -> Option<Vec3> {
    let bc = b.normal.cross(c.normal);
    let det = a.normal.dot(bc);
    if det.abs() <= tangent_sine.max(EFF_ZERO) {
        return None;
    }
    let (ha, hb, hc) = (
        a.point.dot(a.normal),
        b.point.dot(b.normal),
        c.point.dot(c.normal),
    );
    Some((bc * ha + c.normal.cross(a.normal) * hb + a.normal.cross(b.normal) * hc) * (1.0 / det))
}

/// An exact class a local operation's meeting is found on exactly: a plane, or a cone
/// (a cylinder for `tan_angle = 0`) or a sphere, as [`Implicit`] holds them.
#[derive(Debug, Clone, Copy)]
enum Exact {
    Plane { origin: Vec3, normal: Vec3 },
    Quadric(Implicit),
}

/// `surface`'s exact class: its analytic label (a plane, cylinder, cone or sphere) with its
/// NURBS form proved on that label's implicit equation ([`implicit::lies_on`], the owner's
/// rule for an exact class at bound 0); `None` for any other kind, a torus (its sections
/// are quartics), or a form not proved on its label.
fn exact_class(k: &Kernel, surface: &Surface) -> Result<Option<Exact>> {
    let Some(support) = AnalyticSupport::of(surface) else {
        return Ok(None);
    };
    if matches!(support, AnalyticSupport::Torus { .. }) {
        return Ok(None);
    }
    let Some(form) = crate::blend::fillet_section::station::as_b_spline(k, surface)? else {
        return Ok(None);
    };
    let equation = support.implicit();
    if !implicit::lies_on(&form, &equation)? {
        return Ok(None);
    }
    Ok(Some(match equation {
        Implicit::Plane { origin, normal } => Exact::Plane {
            origin,
            normal: normal.unitize()?,
        },
        Implicit::Cone {
            origin,
            axis,
            radius,
            tan_angle,
        } => Exact::Quadric(Implicit::Cone {
            origin,
            axis: axis.unitize()?,
            radius,
            tan_angle,
        }),
        other => Exact::Quadric(other),
    }))
}

/// Where two surfaces meet in an edge's normal section, found on their exact equations
/// ([`exact_meeting`]).
#[derive(Debug, Clone, Copy, PartialEq)]
enum Meeting {
    /// They meet at this point, the meeting nearest the edge point.
    At(Vec3),
    /// They do not meet in the section, nor come within the B-rep's tolerance of it.
    Apart,
    /// Not both of exact classes with one a plane: not found exactly.
    NotExact,
}

/// The coefficients `(A, B, C)` of `F(x0 + τ d) = A τ² + B τ + C` for the quadric's
/// implicit function `F` (a cone's `|P v|² − (r + t z)²`, `v = x − origin`, `z = v·a`, `P`
/// the projection off the axis; a sphere's `|v|² − r²`), with the cone's radius
/// `r + t z` at `x0` and its rate along `d` (zero for a sphere).
fn along_line(q: &Implicit, x0: Vec3, d: Vec3) -> Option<[f64; 5]> {
    match *q {
        Implicit::Cone {
            origin,
            axis,
            radius,
            tan_angle,
        } => {
            let v = x0 - origin;
            let (z0, zd) = (v.dot(axis), d.dot(axis));
            let (p0, pd) = (v - axis * z0, d - axis * zd);
            let (r0, rd) = (radius + tan_angle * z0, tan_angle * zd);
            Some([
                pd.dot(pd) - rd * rd,
                2.0 * (p0.dot(pd) - r0 * rd),
                p0.dot(p0) - r0 * r0,
                r0,
                rd,
            ])
        }
        Implicit::Sphere { centre, radius } => {
            let v = x0 - centre;
            Some([
                d.dot(d),
                2.0 * v.dot(d),
                v.dot(v) - radius * radius,
                0.0,
                0.0,
            ])
        }
        _ => None,
    }
}

/// The distance from `x` to the quadric `q`: a cone's (a cylinder's) from its generator
/// in the meridian half-plane, `|ρ − (r + t z)| / √(1 + t²)`; a sphere's `||v| − r|`.
fn distance_to(q: &Implicit, x: Vec3) -> Option<f64> {
    match *q {
        Implicit::Cone {
            origin,
            axis,
            radius,
            tan_angle,
        } => {
            let v = x - origin;
            let z = v.dot(axis);
            let rho = (v - axis * z).length();
            Some((rho - radius - tan_angle * z).abs() / sqrt(1.0 + tan_angle * tan_angle))
        }
        Implicit::Sphere { centre, radius } => Some(((x - centre).length() - radius).abs()),
        _ => None,
    }
}

/// Where the plane `(origin, normal)` meets the quadric `q` in the normal section through
/// `p` spanned by the unit, perpendicular `m` and `n` (module docs of
/// [`local_extension_reach`]): the plane's trace there is the line `x0 + τ d` (`x0` the
/// foot of `p` on it), the quadric's equation along it `A τ² + B τ + C = 0`, and the meeting
/// the root nearest `p` (smallest `|τ|`; on a cone, on its nappe `r + t z ≥ 0`). Without a
/// real root the line passes the quadric; where it passes within `tolerance` of it
/// (measured at `τ = −B / 2A`, the least of `F` along it) they touch there, else they are
/// apart. A plane whose trace is not a line (the section within rounding of parallel to
/// it) is [`Meeting::NotExact`].
fn plane_meets_quadric(
    p: Vec3,
    m: Vec3,
    n: Vec3,
    (origin, normal): (Vec3, Vec3),
    q: &Implicit,
    tolerance: f64,
) -> Meeting {
    let (nm, nn) = (normal.dot(m), normal.dot(n));
    let nu2 = nm * nm + nn * nn;
    if nu2 <= EFF_ZERO {
        return Meeting::NotExact;
    }
    let h = (origin - p).dot(normal);
    let x0 = p + (m * nm + n * nn) * (h / nu2);
    let d = (m * nn - n * nm) * (1.0 / sqrt(nu2));
    let Some([a, b, c, r0, rd]) = along_line(q, x0, d) else {
        return Meeting::NotExact;
    };
    let on_nappe = |tau: f64| r0 + rd * tau >= 0.0;
    let scale = a.abs().max(b.abs()).max(c.abs());
    let mut roots: Vec<f64> = Vec::new();
    if a.abs() <= EFF_ZERO {
        // the line runs along the quadric's ruling: one root, or none
        if b.abs() <= EFF_ZERO * scale.max(1.0) {
            return Meeting::NotExact;
        }
        roots.push(-c / b);
    } else {
        let disc = b * b - 4.0 * a * c;
        if disc < 0.0 {
            let tau = -b / (2.0 * a);
            let x = x0 + d * tau;
            let near = distance_to(q, x).is_some_and(|gap| gap <= tolerance) && on_nappe(tau);
            return if near { Meeting::At(x) } else { Meeting::Apart };
        }
        // the root of larger size from the sum, the other from the product (no cancellation)
        let w = -0.5 * (b + b.signum() * sqrt(disc));
        // (`w` vanishes only where `B` and the discriminant do, so `C` too: a double root
        // at the foot)
        if w.abs() > 0.0 {
            roots.push(w / a);
            roots.push(c / w);
        } else {
            roots.push(0.0);
        }
    }
    roots
        .into_iter()
        .filter(|&tau| tau.is_finite() && on_nappe(tau))
        .min_by(|x, y| x.abs().total_cmp(&y.abs()))
        .map_or(Meeting::Apart, |tau| Meeting::At(x0 + d * tau))
}

/// Where the new surfaces `own` and `other` of an edge's two faces meet in its normal
/// section at `p` (spanned by the unit, perpendicular `m` and `n`), when that can be found
/// on their exact equations: a plane and a plane (their line, the point of it nearest `p`;
/// planes within `tangent_sine` of parallel are one within `tolerance`, else apart, as
/// [`meeting_near`] has them), a plane and a cone, cylinder or sphere
/// ([`plane_meets_quadric`]). Two curved quadrics (their sections meet in a quartic) or any
/// other kind are [`Meeting::NotExact`].
#[allow(clippy::too_many_arguments)] // the section's frame, both surfaces and both tolerances
fn exact_meeting(
    p: Vec3,
    m: Vec3,
    n: Vec3,
    (own, own_tangent): (Option<Exact>, Tangent),
    (other, other_tangent): (Option<Exact>, Tangent),
    tangent_sine: f64,
    tolerance: f64,
) -> Meeting {
    match (own, other) {
        (Some(Exact::Plane { .. }), Some(Exact::Plane { .. })) => {
            meeting_near(p, own_tangent, other_tangent, tangent_sine, tolerance)
                .map_or(Meeting::Apart, Meeting::At)
        }
        (Some(Exact::Plane { origin, normal }), Some(Exact::Quadric(q)))
        | (Some(Exact::Quadric(q)), Some(Exact::Plane { origin, normal })) => {
            plane_meets_quadric(p, m, n, (origin, normal), &q, tolerance)
        }
        _ => Meeting::NotExact,
    }
}

/// One end of an edge use for [`local_extension_reach`]: the edge point and the face's plane
/// there.
#[derive(Debug, Clone, Copy)]
struct LocalEnd {
    point: Vec3,
    own: Tangent,
}

/// The surface a face takes in a local operation: its own (the face stays) or the one it
/// moves to.
fn local_tangent(
    k: &Kernel,
    face: Option<FaceId>,
    moved: &[(FaceId, &Surface)],
    p: Vec3,
    normal: Vec3,
) -> Result<Tangent> {
    match face.and_then(|f| moved.iter().find(|(m, _)| *m == f)) {
        Some((_, surface)) => tangent_near(k, surface, p),
        None => Ok(Tangent {
            point: p,
            normal: normal.unitize()?,
            moved: false,
        }),
    }
}

/// How far past its face the surface of `faceuse` must be extended in a local operation for
/// the surfaces of its neighbours to meet its own, and how far the surfaces moved there
/// (module docs of [`extension_reach`], with each face's offset replaced by the surface it
/// takes): the faces in `moved` take the surface given with them, every other keeps its
/// own (the copy given with it in `kept`, else the face's surface). The move is the largest
/// distance from an edge point of a moved surface's nearest point. Where every surface
/// stays, both are zero. `None` for a face with no neighbour (every edge a lamina edge or a
/// seam). `tolerance` is the B-rep's: surfaces closer than it are one.
///
/// # The reach of an edge
///
/// The new edge between the face and its neighbour crosses the edge's normal section at
/// each sample `p` (the plane through `p` spanned by the outward binormal `m` and the face
/// normal `n`) where the traces of the two new surfaces there meet; the reach is the
/// largest distance from `p` of those meetings, at each sample (the sector varies along
/// the edge) of each edge where a surface beside it moves.
///
/// Where the two surfaces are of exact classes ([`implicit::lies_on`] of their labels) and
/// one is a plane, the meeting is found on their equations ([`exact_meeting`]): two planes
/// meet on their line, a plane's trace `x0 + τ d` meets a cylinder's, cone's or sphere's
/// where `A τ² + B τ + C = 0`, the root nearest `p`. This is exact to rounding at any angle
/// between the surfaces, so it holds in both limits. Across an edge where they cross at a
/// sector angle `α`, a move `δ` of one is met `δ / sin α` along the other's trace, to first
/// order (the tangent planes' meeting). Across a tangent edge the planes no longer meet
/// near (`sin α → 0`), and the meeting is set by how the traces curve apart: the traces
/// at heights `δ + ½ k₁ s²` and `½ k₂ s²` over the common tangent (normal curvatures `k₁`,
/// `k₂` across the edge) meet where `½ |k₁ − k₂| s² = δ`, `s = √(2 δ / |k₁ − k₂|)` (the
/// second-order contact of two curves; a round of radius `r` beside a plane moved `δ`
/// into it is met `√(2 r δ − δ²)` along the plane from the old tangent line, `√(2 r δ)`
/// (1 for `r = 5`, `δ = 0.1`) from the old edge, where the planes give nothing). With both terms, `δ + θ s + ½ κ s² = 0`, the quadratic
/// is the plane–quadric one, solved exactly.
///
/// Where they have no meeting in the section ([`Meeting::Apart`]: tangent planes apart,
/// equal curvature, or a curved trace bending away from the plane's), no extension of either
/// makes them meet near the edge: [`OffsetFailure::FacesApart`], unless `step_walls`
/// (`Kernel::fixes.local_operation_step_walls`) closes the gap with a wall, when the edge
/// adds nothing.
///
/// Any other pair (two curved quadrics, a torus, a B-spline or a blend) is taken as its
/// tangent planes (a moved surface's at its point nearest the edge): their meeting is the
/// first-order reach ([`meeting_near`]), off the true meeting by the bend over the reach,
/// as the offset's reach is ([`extension_reach`]). Where those planes are tangent (within
/// `tangency_tol_radians`) and apart (beyond `tolerance`), the meeting is set by the
/// surfaces' relative curvature alone, and bounding `s` from above takes a certified
/// *lower* bound on `|k₁ − k₂|` over the region (a bound from a net bounds the curvature
/// from above, which bounds the meeting from below only): none is had here, so the reach
/// is not bounded: [`OffsetFailure::ReachUncertified`], unless `step_walls`.
///
/// At each corner, the reach is the largest distance of the points where the face's plane
/// and two others of the faces at its vertex meet (where three meet in no one point, the
/// edges' meetings stand), for every two of them, not only its edge neighbours there: where
/// four or more faces meet and one across the vertex moves, the vertex splits and the
/// face's new edges run to the split points, which are among these meetings (first order,
/// as the tangent planes are). Every pair is taken, so the reach bounds the split points
/// from above (a meeting that is no vertex of the result only extends the face further).
/// A face at the vertex moved adds its move there.
///
/// # Errors
///
/// [`OffsetFailure::FacesApart`] and [`OffsetFailure::ReachUncertified`] as above;
/// [`SolidError::InvalidInput`] for a vanishing face normal; those of
/// [`eu3::evaluate_binormal`], of a moved surface's closest point and normal, and of the
/// exact classes' proof ([`implicit::lies_on`]); [`SolidError::NotFound`] for a dangling
/// id.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)] // the operation's surfaces and flags
pub fn local_extension_reach(
    k: &Kernel,
    b: &mut Brep,
    faceuse: FaceuseId,
    moved: &[(FaceId, &Surface)],
    kept: &[(FaceId, &Surface)],
    tangency_tol_radians: f64,
    step_walls: bool,
) -> Result<Option<(f64, f64)>> {
    let tangent_sine = tangency_tol_radians.abs().min(crate::math::PI / 2.0).sin();
    let tolerance = b.tolerance();
    let options = eu3::BinormalOptions::default();
    let (mut reach, mut shift) = (0.0_f64, 0.0_f64);
    let mut met = false;
    // each face's exact class, proved once
    let mut classes: Vec<(FaceId, Option<Exact>)> = Vec::new();
    let mut class_of = |b: &Brep, face: FaceId| -> Result<Option<Exact>> {
        if let Some((_, c)) = classes.iter().find(|(f, _)| *f == face) {
            return Ok(*c);
        }
        let given = moved
            .iter()
            .chain(kept)
            .find(|(f, _)| *f == face)
            .map(|(_, s)| *s);
        let c = match given {
            Some(s) => exact_class(k, s)?,
            None => match crate::topology::face::p2::surface(b, face)? {
                Some(s) => exact_class(k, s)?,
                None => None,
            },
        };
        classes.push((face, c));
        Ok(c)
    };
    let refuse =
        |face: Option<FaceId>, why: fn(FaceId) -> OffsetFailure| -> Result<Option<(f64, f64)>> {
            let face = face.ok_or_else(|| SolidError::null_pointer(LOCAL_OP))?;
            Err(Failure::new(LOCAL_OP, why(face)).into_error())
        };
    for lu in eu::faceuse_loopuses(b, faceuse, false, None)? {
        let mut ends: Vec<(Option<crate::topology::VertexId>, LocalEnd, LocalEnd)> = Vec::new();
        for edgeuse in loopuse::edgeuses(b, lu, None, None)? {
            let Some(e) = eu::edge(b, edgeuse)? else {
                continue;
            };
            let interval = edge::interval(b, e)?;
            let (lo, hi) = (interval.min(), interval.max());
            let own_face = eu::face(b, edgeuse)?;
            let radial = match eu::radial(b, edgeuse)? {
                Some(r) if r != edgeuse && eu::face(b, r)? != own_face => Some(r),
                _ => None,
            };
            let mut first = None;
            let mut last = None;
            for i in 0..=EDGE_SAMPLES {
                let parameter = lo + (hi - lo) * f64::from(i) / f64::from(EDGE_SAMPLES);
                let at = eu3::evaluate_binormal(k, b, edgeuse, parameter, false, &options)?;
                let own = local_tangent(k, own_face, moved, at.point, at.face_normal)?;
                shift = shift.max((own.point - at.point).length());
                if let Some(r) = radial {
                    let theirs = eu3::evaluate_binormal(k, b, r, parameter, false, &options)?;
                    let face = eu::face(b, r)?;
                    let other = local_tangent(k, face, moved, at.point, theirs.face_normal)?;
                    shift = shift.max((other.point - at.point).length());
                    met = true;
                    // where both surfaces stay they meet at the edge: nothing to reach
                    if own.moved || other.moved {
                        let m = (-at.binormal).unitize()?;
                        let n = at.face_normal.unitize()?;
                        let n = (n - m * n.dot(m)).unitize()?;
                        let own_class = match own_face {
                            Some(f) => class_of(b, f)?,
                            None => None,
                        };
                        let other_class = match face {
                            Some(f) => class_of(b, f)?,
                            None => None,
                        };
                        let meeting = exact_meeting(
                            at.point,
                            m,
                            n,
                            (own_class, own),
                            (other_class, other),
                            tangent_sine,
                            tolerance,
                        );
                        let x = match meeting {
                            Meeting::At(x) => Some(x),
                            Meeting::Apart if step_walls => None,
                            Meeting::Apart => {
                                return refuse(own_face, |face| OffsetFailure::FacesApart { face });
                            }
                            Meeting::NotExact => {
                                match meeting_near(at.point, own, other, tangent_sine, tolerance) {
                                    Some(x) => Some(x),
                                    None if step_walls => None,
                                    None => {
                                        return refuse(own_face, |face| {
                                            OffsetFailure::ReachUncertified { face }
                                        });
                                    }
                                }
                            }
                        };
                        if let Some(x) = x {
                            reach = reach.max((x - at.point).length());
                        }
                    }
                }
                let end = LocalEnd {
                    point: at.point,
                    own,
                };
                if i == 0 {
                    first = Some(end);
                }
                last = Some(end);
            }
            let (Some(first), Some(last)) = (first, last) else {
                continue;
            };
            let start = match eu::vertexuse(b, edgeuse)? {
                Some(vu) => crate::topology::vertexuse::vertex(b, vu)?,
                None => None,
            };
            let (into, out) = match eu::orientation(b, edgeuse)? {
                OrientType::Opposite => (last, first),
                _ => (first, last),
            };
            ends.push((start, into, out));
        }
        // Each edge use's end meets the next one's start (the loop closes), at that one's
        // start vertex: the corner there is met by every face around the vertex.
        let next = ends.iter().cycle().skip(1);
        for (&(_, _, into), &(vertex, _, _)) in ends.iter().zip(next) {
            let Some(vertex) = vertex else {
                continue;
            };
            let around = faces_around(k, b, vertex, b.faceuse(faceuse)?.face, moved)?;
            // a face that stays is where it is (its edge point differs from this one by
            // rounding only)
            for (_, t) in around.iter().filter(|(_, t)| t.moved) {
                shift = shift.max((t.point - into.point).length());
            }
            // a corner of faces that all stay stays: it adds nothing (its planes' meeting
            // computed again would differ from the vertex by rounding only)
            if !(into.own.moved || around.iter().any(|(_, t)| t.moved)) {
                continue;
            }
            for (i, &(_, a)) in around.iter().enumerate() {
                for &(_, c) in around.iter().skip(i + 1) {
                    if let Some(x) = corner_of(into.own, a, c, tangent_sine) {
                        reach = reach.max((x - into.point).length());
                    }
                }
            }
        }
    }
    Ok(met.then_some((reach, shift)))
}

/// Every face at `vertex` other than `own` (each once), with the plane it takes there in a
/// local operation ([`local_tangent`], its normal at the vertex from one of its edge uses
/// starting there): its edge neighbours at the vertex and, where four or more faces meet,
/// those across the vertex from it, which share no edge with it.
///
/// # Errors
///
/// Those of [`eu3::evaluate_binormal`] and [`local_tangent`]; [`SolidError::NotFound`] for
/// a dangling id.
fn faces_around(
    k: &Kernel,
    b: &mut Brep,
    vertex: crate::topology::VertexId,
    own: Option<FaceId>,
    moved: &[(FaceId, &Surface)],
) -> Result<Vec<(FaceId, Tangent)>> {
    use crate::topology::{vertex as vx, vertexuse as vu};
    let options = eu3::BinormalOptions::default();
    let mut around: Vec<(FaceId, Tangent)> = Vec::new();
    for use_ in vx::vertexuses(b, vertex)? {
        if !vu::is_edge_vertexuse(b, use_)? {
            continue;
        }
        let Some(edgeuse) = vu::edgeuse(b, use_)? else {
            continue;
        };
        let face = eu::face(b, edgeuse)?;
        let Some(f) = face.filter(|f| Some(*f) != own) else {
            continue;
        };
        if around.iter().any(|(g, _)| *g == f) {
            continue;
        }
        let parameter = vu::compute_t_point(b, use_)?;
        let at = eu3::evaluate_binormal(k, b, edgeuse, parameter, false, &options)?;
        around.push((f, local_tangent(k, face, moved, at.point, at.face_normal)?));
    }
    Ok(around)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Point3;
    use crate::topology::brep::p1 as bp1;

    /// The basin's ruled side: rows 46 long at its rim and 20 at its floor.
    #[test]
    fn a_narrowing_side_is_short_of_its_distance_where_it_narrows() {
        use crate::nurbs::b_spline_curve::KnotType;
        use crate::nurbs::b_spline_surface::{
            BSplineSurface, BSplineSurfaceForm, CanonicalSurfaceData,
        };
        use crate::surface::extend::Continuity;
        let k = Kernel::default();
        let base = Surface::from(
            BSplineSurface::create_canonical(&CanonicalSurfaceData {
                u_degree: 1,
                v_degree: 1,
                control_points: &[
                    Point3::new(20.0, 4.0, 30.0),
                    Point3::new(20.0, 50.0, 30.0),
                    Point3::new(30.0, 20.0, 10.0),
                    Point3::new(30.0, 40.0, 10.0),
                ],
                form: BSplineSurfaceForm::Unspecified,
                u_multiplicities: &[2, 2],
                v_multiplicities: &[2, 2],
                u_knots: &[0.0, 1.0],
                v_knots: &[0.0, 1.0],
                knot_type: KnotType::Unspecified,
                weights: None,
                uv_domain: None,
            })
            .unwrap(),
        );
        let made = base
            .extended_sides_keeping_parameters(&k, &BoxSide::ALL, Continuity::Ruled, 11.24)
            .unwrap();
        // the rows are 33 at the side's middle, where the extension meets its distance, and
        // 20 at the floor, so the v sides reach 11.24 * 20 / 33 at the floor's end
        let shortest = shortest_across(&k, &base, &made).unwrap().unwrap();
        assert!((shortest - 11.24 * 20.0 / 33.0).abs() <= 1e-6, "{shortest}");
    }

    fn box_faceuses(k: &mut Kernel, inset: bool) -> (Brep, Vec<(FaceId, FaceuseId)>) {
        let brep = crate::primitive::api_primitives::create_box(
            k,
            Point3::new(0.0, 0.0, 0.0),
            1.0,
            2.0,
            3.0,
        )
        .unwrap();
        let faceuses = bp1::faces(&brep, None)
            .unwrap()
            .into_iter()
            .map(|f| {
                let fu = super::super::p2::select_offset_faceuse(&brep, f, false, inset, OP);
                (f, fu.unwrap())
            })
            .collect();
        (brep, faceuses)
    }

    /// Across a box's edges (sector angle a right angle) the offsets meet `d` past each
    /// edge, and at its corners `d` past both: `d·√2` from the vertex. Offset inward they meet
    /// as far inside: the same distance from the face.
    #[test]
    fn a_box_offset_either_way_reaches_its_corners_offset() {
        for inset in [false, true] {
            let mut k = Kernel::default();
            let (mut brep, faceuses) = box_faceuses(&mut k, inset);
            for (_, fu) in faceuses {
                let reach = extension_reach(&k, &mut brep, fu, 0.25, &[], 1e-3)
                    .unwrap()
                    .unwrap();
                assert!((reach - 0.25 * sqrt(2.0)).abs() < 1e-12, "{reach}");
            }
        }
    }

    /// A box shelled open at one face: that face stays where it is, and each wall's offset
    /// meets it at the wall's edge (a right angle, `β = 0`). A wall's corner on the open face
    /// lies `d` from its vertex, and its other corners `d·√2`.
    #[test]
    fn an_unmoved_face_is_met_where_it_stands() {
        let mut k = Kernel::default();
        let (mut brep, faceuses) = box_faceuses(&mut k, true);
        let (open, _) = faceuses[0];
        for &(_, fu) in faceuses.iter().filter(|(f, _)| *f != open) {
            let reach = extension_reach(&k, &mut brep, fu, 0.25, &[open], 1e-3)
                .unwrap()
                .unwrap();
            assert!((reach - 0.25 * sqrt(2.0)).abs() < 1e-12, "{reach}");
        }
    }

    /// The signed reach: `d·cot(α/2)` for equal offsets (a right angle `d`, a tangent edge
    /// `0`), `0` for a knife edge (no meeting); an unmoved neighbour at a right angle is met at
    /// the edge, at 60° between normals `dF·c·sin / (1 − c²)` inside.
    #[test]
    fn signed_reach_is_the_foot_of_the_offsets_meeting() {
        let up = Vec3::new(0.0, 0.0, 1.0);
        let out = Vec3::X;
        assert!((signed_reach(up, Vec3::X, out, 0.5, 0.5, 1e-3) - 0.5).abs() < 1e-15);
        assert!((signed_reach(up, -Vec3::X, out, 0.5, 0.5, 1e-3) + 0.5).abs() < 1e-15);
        assert_eq!(signed_reach(up, up, out, 0.5, 0.5, 1e-3), 0.0);
        assert_eq!(signed_reach(up, -up, out, 0.5, 0.5, 1e-3), 0.0);
        assert_eq!(signed_reach(up, up, out, 0.5, 0.0, 1e-3), 0.0);
        assert!(signed_reach(up, Vec3::X, out, 0.5, 0.0, 1e-3).abs() < 1e-15);
        let half = sqrt(3.0) / 2.0;
        let tilted = Vec3::new(half, 0.0, 0.5);
        let expected = -0.5 * 0.5 * half / (1.0 - 0.25);
        assert!((signed_reach(up, tilted, out, 0.5, 0.0, 1e-3) - expected).abs() < 1e-15);
        // within the tangency tolerance of tangent, offsets 0.5 apart never meet near
        let nearly = Vec3::new(1e-4, 0.0, (1.0 - 1e-8_f64).sqrt());
        assert_eq!(signed_reach(up, nearly, out, 0.5, 0.0, 1e-3), 0.0);
        assert!(signed_reach(up, nearly, out, 0.5, 0.0, 1e-5).abs() > 1e3);
    }

    /// The corner's foot: perpendicular edges each reached `1` and `2` past put the corner
    /// at `(1, 2)`; edges at 60° reached `d` past both put it on the bisector, `d / cos 30°`
    /// from the vertex, as reached `d` inside; parallel edges, or edges within the tangency
    /// tolerance of tangent, take the larger reach.
    #[test]
    fn corner_reach_is_the_foot_of_the_offset_corner() {
        let end = |reach: f64, outward: Vec3| End { reach, outward };
        let square = corner_reach(end(1.0, Vec3::X), end(2.0, Vec3::Y), 1e-3);
        assert!((square - sqrt(5.0)).abs() < 1e-12);
        let half = sqrt(3.0) / 2.0;
        let sixty = Vec3::new(0.5, half, 0.0);
        let acute = corner_reach(end(1.0, Vec3::X), end(1.0, sixty), 1e-3);
        assert!((acute - 1.0 / half).abs() < 1e-12, "{acute}");
        let inside = corner_reach(end(-1.0, Vec3::X), end(-1.0, sixty), 1e-3);
        assert!((inside - 1.0 / half).abs() < 1e-12, "{inside}");
        let parallel = corner_reach(end(1.0, Vec3::X), end(-3.0, Vec3::X), 1e-3);
        assert_eq!(parallel, 3.0);
        // edges 1e-4 off tangent, within a tolerance of 1e-3: no corner, not 1e4 off
        let nearly = Vec3::new((1.0 - 1e-8_f64).sqrt(), 1e-4, 0.0);
        let tangent = corner_reach(end(1.0, Vec3::X), end(2.0, nearly), 1e-3);
        assert_eq!(tangent, 2.0);
    }

    /// A moved input is refused: an offset that is not finite.
    #[test]
    fn a_non_finite_offset_is_refused() {
        let mut k = Kernel::default();
        let (mut brep, faceuses) = box_faceuses(&mut k, false);
        assert!(extension_reach(&k, &mut brep, faceuses[0].1, f64::NAN, &[], 1e-3).is_err());
    }

    /// The box's face whose plane has unit normal `n` (outward), with its upward face use.
    fn box_face(brep: &Brep, n: Vec3) -> (FaceId, FaceuseId) {
        for f in bp1::faces(brep, None).unwrap() {
            let surface = brep.face(f).unwrap().surface.as_ref().unwrap();
            let SurfaceKind::Plane(plane) = surface.kind() else {
                continue;
            };
            let z = plane.position().z_axis();
            let fu = super::super::p2::select_offset_faceuse(brep, f, false, false, OP).unwrap();
            let up = brep.faceuse(fu).unwrap().orientation != OrientType::Opposite;
            let outward = if up { z } else { -z };
            if (outward - n).length() < 1e-12 {
                return (f, fu);
            }
        }
        panic!("no face of normal {n:?}");
    }

    /// The box [0, 1] x [0, 2] x [0, 3] with its top (z = 3) moved to `top`: each face's
    /// (reach, move).
    fn moved_top(top: &Surface) -> [(f64, f64); 4] {
        let k = Kernel::default();
        let mut kk = Kernel::default();
        let (mut brep, _) = box_faceuses(&mut kk, false);
        let (top_face, top_fu) = box_face(&brep, Vec3::Z);
        let moved = [(top_face, top)];
        let mut at = |n: Vec3| {
            let fu = if n == Vec3::Z {
                top_fu
            } else {
                box_face(&brep, n).1
            };
            local_extension_reach(&k, &mut brep, fu, &moved, &[], 1e-3, false)
                .unwrap()
                .unwrap()
        };
        [at(Vec3::Z), at(-Vec3::Y), at(-Vec3::X), at(-Vec3::Z)]
    }

    /// The top raised by `h`: the sides meet it `h` above their top edges (and corners), the
    /// top meets the sides `h` above its old edges, both moved `h`; the bottom meets
    /// sides that stay: nothing to reach, nothing moved.
    #[test]
    fn a_raised_top_is_met_as_far_as_it_rose() {
        let h = 0.5;
        let top = Surface::from(crate::surface::plane::Plane::from_normal(
            Point3::new(0.0, 0.0, 3.0 + h),
            Vec3::Z,
        ));
        let [top, front, left, bottom] = moved_top(&top);
        for (reach, shift) in [top, front, left] {
            assert!((reach - h).abs() < 1e-12, "{reach}");
            assert!((shift - h).abs() < 1e-12, "{shift}");
        }
        assert_eq!(bottom, (0.0, 0.0));
    }

    /// The top tilted by `θ` about the line `y = 1, z = 3`: `z = 3 + (y - 1) tan θ`. The
    /// front (y = 0) meets it `tan θ` below its top edge, and moves `sin θ` there (the
    /// edge's distance from the tilted plane); the left side (x = 0) meets it at its
    /// corners `tan θ` off (along its top edge only `|y - 1| sin θ`); the top's own meetings
    /// are `tan θ` from its old edges at the front and back.
    #[test]
    fn a_tilted_top_is_met_by_its_tangent_at_the_corners() {
        let theta: f64 = 0.3;
        let (s, c) = (theta.sin(), theta.cos());
        let top = Surface::from(crate::surface::plane::Plane::from_normal(
            Point3::new(0.5, 1.0, 3.0),
            Vec3::new(0.0, -s, c),
        ));
        let [top, front, left, bottom] = moved_top(&top);
        let t = theta.tan();
        for (name, (reach, shift)) in [("top", top), ("front", front), ("left", left)] {
            assert!((reach - t).abs() < 1e-12, "{name} reach {reach}");
            assert!((shift - s).abs() < 1e-12, "{name} move {shift}");
        }
        assert_eq!(bottom, (0.0, 0.0));
    }

    /// The pyramid on the square [-1, 1]² with its apex at (0, 0, 1), its side y = -1..0
    /// (`z = 1 + y`) moved in by `δ` along its normal: `z = 1 + y - a`, `a = δ √2`. Four
    /// sides meet at the apex, and the side across it (`z = 1 - y`) shares no edge with the
    /// moved one: its edge neighbours (`z = 1 ∓ x`) and the base stay. The moved side cuts
    /// the apex off, and the side across meets it on a new edge from
    /// (±a/2, a/2, 1 - a/2), `a √3 / 2 = δ √(3/2)` from the apex; the move there is `δ`
    /// (the apex's distance from the moved plane). With only its edge neighbours read, the
    /// side across had neither.
    #[test]
    fn a_side_across_a_four_face_vertex_from_a_moved_one_reaches_its_split() {
        let k = Kernel::default();
        let mut kk = Kernel::default();
        let mut brep =
            crate::primitive::api_primitives::create_pyramid(&mut kk, Point3::ZERO, 2.0, 1.0)
                .unwrap();
        let s = sqrt(0.5);
        let delta = 0.1;
        let (moved_face, _) = box_face(&brep, Vec3::new(0.0, -s, s));
        let (_, across) = box_face(&brep, Vec3::new(0.0, s, s));
        let n = Vec3::new(0.0, -s, s);
        let plane = Surface::from(crate::surface::plane::Plane::from_normal(
            Point3::new(0.0, 0.0, 1.0) - n * delta,
            n,
        ));
        let moved = [(moved_face, &plane)];
        let (reach, shift) = local_extension_reach(&k, &mut brep, across, &moved, &[], 1e-3, false)
            .unwrap()
            .unwrap();
        assert!((reach - delta * sqrt(1.5)).abs() < 1e-12, "{reach}");
        assert!((shift - delta).abs() < 1e-12, "{shift}");
    }

    /// Nothing moved: every face meets its neighbours at its edges.
    #[test]
    fn a_body_nothing_moves_in_needs_no_extension() {
        let k = Kernel::default();
        let mut kk = Kernel::default();
        let (mut brep, faceuses) = box_faceuses(&mut kk, false);
        for (_, fu) in faceuses {
            let found = local_extension_reach(&k, &mut brep, fu, &[], &[], 1e-3, false).unwrap();
            assert_eq!(found, Some((0.0, 0.0)));
        }
    }

    /// A 20 mm cube with its vertical edge at x = y = 20 rounded r 5 (a cylinder about
    /// (15, 15)), and its faces: the round, the face on x = 20 and the bottom.
    fn rounded_cube(k: &mut Kernel) -> (Brep, FaceId, FaceId, FaceId) {
        use crate::blend::fillet_executive::api::fillet_edges;
        use crate::topology::vertex;
        let mut b = crate::primitive::api_primitives::create_box(k, Point3::ZERO, 20.0, 20.0, 20.0)
            .unwrap();
        let e = b
            .edges()
            .ids()
            .find(|&e| {
                let (s, t) = edge::start_end_vertices(&b, e).unwrap();
                [s, t].into_iter().flatten().all(|v| {
                    let p = vertex::point(&b, v).unwrap();
                    p.x == 20.0 && p.y == 20.0
                })
            })
            .unwrap();
        fillet_edges(k, &mut b, &[e], 5.0, 1, 1, 1.0).unwrap();
        let mut round = None;
        let mut side = None;
        let mut bottom = None;
        for f in bp1::faces(&b, None).unwrap() {
            let face = b.face(f).unwrap();
            let surface = face.surface.as_ref().unwrap();
            let p = surface
                .evaluate_point(k, face.uv_domain.evaluate(0.5, 0.5))
                .unwrap();
            match surface.kind() {
                SurfaceKind::Cone(_) => round = Some(f),
                SurfaceKind::Plane(_) if (p.x - 20.0).abs() < 1e-9 => side = Some(f),
                SurfaceKind::Plane(_) if p.z.abs() < 1e-9 => bottom = Some(f),
                _ => {}
            }
        }
        (b, round.unwrap(), side.unwrap(), bottom.unwrap())
    }

    fn faceuse_of(b: &Brep, f: FaceId) -> FaceuseId {
        super::super::p2::select_offset_faceuse(b, f, false, false, OP).unwrap()
    }

    /// Review finding 6: the round's side face (x = 20, tangent to the round along
    /// y = 15) moved into the body by 0.1 meets the round `√(2 r δ − δ²) = √0.99` along it
    /// from the old tangent line (`½ s² / r ≈ δ`), `√(2 r δ) = 1` from the old edge's
    /// points, where the tangent planes, parallel, gave nothing;
    /// the side's other edges, across right angles, reach 0.1, and so does the bottom where
    /// it meets the moved side. The round, which stays, is reached as far at that edge.
    #[test]
    fn a_plane_moved_into_a_round_is_reached_where_it_meets_it() {
        let mut k = Kernel::default();
        let (mut b, round, side, bottom) = rounded_cube(&mut k);
        let moved_side = Surface::from(crate::surface::plane::Plane::from_normal(
            Point3::new(19.9, 0.0, 0.0),
            Vec3::X,
        ));
        let moved = [(side, &moved_side)];
        let want = sqrt(2.0 * 5.0 * 0.1);
        for (f, expected) in [(side, want), (round, want), (bottom, 0.1)] {
            let fu = faceuse_of(&b, f);
            let (reach, shift) = local_extension_reach(&k, &mut b, fu, &moved, &[], 1e-3, false)
                .unwrap()
                .unwrap();
            assert!(
                (reach - expected).abs() < 1e-9,
                "{f:?}: {reach} vs {expected}"
            );
            assert!((shift - 0.1).abs() < 1e-9, "{f:?}: {shift}");
        }
    }

    /// The round given the coaxial cylinder r 5.1 (moved out 0.1): the side faces, which
    /// stay, meet it `√(5.1² − 5²) ≈ 1.005` past their tangent edges, the bottom 0.1
    /// past its arc.
    #[test]
    fn a_round_moved_out_is_reached_by_its_tangent_planes() {
        let mut k = Kernel::default();
        let (mut b, round, side, bottom) = rounded_cube(&mut k);
        let wider = Surface::from(
            crate::primitive::api_primitives::create_cylinder_surface(
                &k,
                Point3::new(15.0, 15.0, -1.0),
                5.1,
                22.0,
            )
            .unwrap(),
        );
        let moved = [(round, &wider)];
        let want = sqrt(5.1 * 5.1 - 25.0);
        for (f, expected) in [(side, want), (round, want), (bottom, 0.1)] {
            let fu = faceuse_of(&b, f);
            let (reach, _) = local_extension_reach(&k, &mut b, fu, &moved, &[], 1e-3, false)
                .unwrap()
                .unwrap();
            assert!(
                (reach - expected).abs() < 1e-9,
                "{f:?}: {reach} vs {expected}"
            );
        }
    }

    /// The side moved out by 0.1, off the round: the plane x = 20.1 passes the cylinder,
    /// whose trace bends away from it: no extension makes them meet, refused as apart;
    /// with the step walls asked for, the edge adds nothing (the wall closes the gap) and
    /// the side's other edges reach 0.1.
    #[test]
    fn a_plane_moved_off_a_round_is_apart() {
        let mut k = Kernel::default();
        let (mut b, _, side, _) = rounded_cube(&mut k);
        let off = Surface::from(crate::surface::plane::Plane::from_normal(
            Point3::new(20.1, 0.0, 0.0),
            Vec3::X,
        ));
        let moved = [(side, &off)];
        let fu = faceuse_of(&b, side);
        let e = local_extension_reach(&k, &mut b, fu, &moved, &[], 1e-3, false).unwrap_err();
        assert!(
            matches!(
                e.failure().map(|f| &f.reason),
                Some(crate::error::Reason::Offset(OffsetFailure::FacesApart { face })) if *face == side
            ),
            "{e}"
        );
        let (reach, _) = local_extension_reach(&k, &mut b, fu, &moved, &[], 1e-3, true)
            .unwrap()
            .unwrap();
        assert!((reach - 0.1).abs() < 1e-9, "{reach}");
    }

    /// Nothing moved on the rounded cube: every face meets its neighbours at its edges,
    /// tangent ones too.
    #[test]
    fn a_rounded_body_nothing_moves_in_needs_no_extension() {
        let mut k = Kernel::default();
        let (mut b, _, _, _) = rounded_cube(&mut k);
        for f in bp1::faces(&b, None).unwrap() {
            let fu = faceuse_of(&b, f);
            let found = local_extension_reach(&k, &mut b, fu, &[], &[], 1e-3, false).unwrap();
            assert_eq!(found, Some((0.0, 0.0)));
        }
    }

    /// The round kept as its bare B-spline net (no exact class): beside the side moved in,
    /// tangent and apart, the meeting is set by the curvature, which is not bounded for
    /// that kind: refused as uncertified (with the step walls asked for, it adds nothing).
    #[test]
    fn a_tangent_neighbour_of_no_exact_class_is_uncertified() {
        let mut k = Kernel::default();
        let (mut b, round, side, _) = rounded_cube(&mut k);
        let net = Surface::from(
            crate::blend::fillet_section::station::as_b_spline(
                &k,
                b.face(round).unwrap().surface.as_ref().unwrap(),
            )
            .unwrap()
            .unwrap()
            .into_owned(),
        );
        let moved_side = Surface::from(crate::surface::plane::Plane::from_normal(
            Point3::new(19.9, 0.0, 0.0),
            Vec3::X,
        ));
        let moved = [(side, &moved_side)];
        let kept = [(round, &net)];
        let fu = faceuse_of(&b, side);
        let e = local_extension_reach(&k, &mut b, fu, &moved, &kept, 1e-3, false).unwrap_err();
        assert!(
            matches!(
                e.failure().map(|f| &f.reason),
                Some(crate::error::Reason::Offset(OffsetFailure::ReachUncertified { face })) if *face == side
            ),
            "{e}"
        );
        assert!(local_extension_reach(&k, &mut b, fu, &moved, &kept, 1e-3, true).is_ok());
    }

    /// The exact meetings: a plane's trace crossing a unit circle (a cylinder's section
    /// across its axis) at height `h` meets it `√(1 − h²)` along, the nearer root; a trace
    /// passing it is apart, unless within the tolerance; a plane trace and a sphere's great
    /// circle the same.
    #[test]
    fn a_plane_meets_a_quadric_at_its_nearest_root() {
        let cylinder = Implicit::Cone {
            origin: Vec3::ZERO,
            axis: Vec3::Z,
            radius: 1.0,
            tan_angle: 0.0,
        };
        let sphere = Implicit::Sphere {
            centre: Vec3::ZERO,
            radius: 1.0,
        };
        // the section z = 0 at p = (0, 1, 0), m = x, n = y: the plane y = 1 - d
        let (p, m, n) = (Vec3::new(0.0, 1.0, 0.0), Vec3::X, Vec3::Y);
        for q in [cylinder, sphere] {
            for d in [0.1, 0.5, 1e-6] {
                let plane = (Vec3::new(0.0, 1.0 - d, 0.0), Vec3::Y);
                let Meeting::At(x) = plane_meets_quadric(p, m, n, plane, &q, 1e-9) else {
                    panic!("{d}");
                };
                let want = sqrt(1.0 - (1.0 - d) * (1.0 - d));
                assert!((x.x.abs() - want).abs() < 1e-12, "{d}: {x:?}");
                assert!((x.y - (1.0 - d)).abs() < 1e-15);
            }
            let past = (Vec3::new(0.0, 1.1, 0.0), Vec3::Y);
            assert_eq!(plane_meets_quadric(p, m, n, past, &q, 1e-9), Meeting::Apart);
            let touching = (Vec3::new(0.0, 1.0 + 1e-10, 0.0), Vec3::Y);
            assert!(matches!(
                plane_meets_quadric(p, m, n, touching, &q, 1e-9),
                Meeting::At(_)
            ));
            // across: the plane x = 0.2 meets the circle at y = √0.96, not at its far side
            let across = (Vec3::new(0.2, 0.0, 0.0), Vec3::X);
            let Meeting::At(x) = plane_meets_quadric(p, m, n, across, &q, 1e-9) else {
                panic!("across");
            };
            assert!((x.y - sqrt(0.96)).abs() < 1e-12, "{x:?}");
        }
    }
}
