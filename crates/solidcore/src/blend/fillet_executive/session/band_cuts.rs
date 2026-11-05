//! Whether a convex blend would take all of the part between its rails somewhere along its
//! edge, cutting it into pieces, and the largest size that keeps it whole.
//!
//! Upstream: not in SMLib at the pinned commit. SolidCraft's fork (tangency-labs/smlib#79,
//! merged at 5e8c9c1) adds `SmFilletSolver::BandCutsThePart` with `sm_EdgesCrossingPlane`
//! and `sm_KeepsMaterialBehind` (`SmFilletSolver.cpp`, `SmFilletSolver.h`) in commits
//! 862971f, e568e45 (the section intersected with every face, not sampled) and aac4288 (a
//! run of cutting sections reaching an end of the edge is no cut). They run only with the
//! opt-in `Kernel::fixes.fillet_gn_section_generator`.
//!
//! | Rust                          | Upstream (fork)                          |
//! |-------------------------------|------------------------------------------|
//! | [`band_cuts_the_part`]        | `SmFilletSolver::BandCutsThePart`        |
//! | [`BandCut`]                   | its outputs `rbCuts`, `rdRadius`, `rdSetback` |
//! | `Station::crossings_of` (private) | `sm_EdgesCrossingPlane`, one edge   |
//! | `keeps_material_behind` (private) | `sm_KeepsMaterialBehind`             |
//! | `keeps` (private)             | the `Keeps` lambda of `BandCutsThePart`  |
//!
//! Not the fork's: [`band_breaks_through_the_part`] (WR-707, #2299), the same sections
//! judged for breaking through the part into a hollow behind a wall, with
//! `fixes.fillet_refuse_rounds_breaking_through`; and [`BandCheck`] (#2648), what a check
//! keeps across its solvers and sections.
//!
//! Deviations: the fork's out-params are returned (`None` for no cut). Each section is
//! judged only against the edges and faces within its reach (#2648, [`BandCheck`]): the
//! answer is the fork's (an edge or face out of reach has no crossing that counts), but an
//! error the fork would meet on an out-of-reach edge or face is not met. A station's plane
//! crossings are found edge by edge as its sections ask for them, not all up front; each
//! face's trimmed cache is made once per check, not once per section; the body's cache (the
//! trees a classification solves through) is made once per check and again when a
//! tolerance of the body differs from the ones it was made with, not once per
//! classification; the check polls the
//! kernel's cancel token. The solutions of
//! `SmSurface::GlobalCurveIntersect` carry no object handles (the fork reads only their
//! parameters). Every loop is bounded by fixed counts (32 steps, 40 bisections, halvings
//! down to a millionth of the radius); there is no recursion.

// the fork's NaN-safe tests (`!(a > b)`) are kept as written
#![allow(clippy::neg_cmp_op_on_partial_ord)]
use crate::blend::fillet_executive::session::ball_regular::{
    ball_law, ball_section_at, rolls_a_ball, set_ball_radius,
};
use crate::blend::fillet_executive::session::SolverObject;
use crate::blend::fillet_solver::p1::FilletSurfaceGeneratorType;
use crate::curve::curve::Curve;
use crate::curve::curve::CurveKind;
use crate::curve::ellipse::p2::EllipseOptions;
use crate::curve::ellipse::Ellipse;
use crate::curve::line::Line;
use crate::error::EntityRef;
use crate::error::{Result, SolidError};
use crate::geom::axis2_placement::Axis2Placement;
use crate::geom::extent1d::Extent1d;
use crate::geom::extent2d::Extent2d;
use crate::geom::extent3d::Extent3d;
use crate::geom::solution_array::{ObjectKind, SolutionObject, SolutionType};
use crate::geom::tol::{brep_zone_tol_3d, SmallTopoType};
use crate::intersect::point_class::{PointClassType, PointClassification};
use crate::intersect::surface_intersector::p1::{SrfIndex, TsectPnt};
use crate::kernel::Kernel;
use crate::math::{rad_to_deg, Point2, Point3, Vec3, EFF_ZERO};
use crate::nurbs::b_spline_curve::BoundingBoxRequest;
use crate::query::brep::{BrepFaceServices, BrepQuery, FaceCaches};
use crate::surface::plane::p1::Plane;
use crate::surface::surface::{Surface, SurfaceKind};
use crate::topology::brep_cache::CacheObject;
use crate::topology::edge;
use crate::topology::edgeuse::p1 as eu;
use crate::topology::{face, vertex, Brep, EdgeId, FaceId};
use std::collections::{BTreeMap, BTreeSet};

const OP: &str = "blend::fillet_executive::session::band_cuts";

/// `lSteps`: the sections along the edge.
const STEPS: usize = 32;

/// `lBisections`: the bisection steps on the ball's radius.
const BISECTIONS: usize = 40;

/// `smos_Max` of two doubles.
fn smos_max(a: f64, b: f64) -> f64 {
    if a > b {
        a
    } else {
        b
    }
}

/// Where a convex blend would cut the part ([`band_cuts_the_part`]): the largest ball
/// radius that keeps it whole, least over the cutting sections (`rdRadius`), and that
/// radius's setback from the edge to the first rail (`rdSetback`, a chamfer's size).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BandCut {
    /// `rdRadius`.
    pub radius: f64,
    /// `rdSetback`.
    pub setback: f64,
}

/// The operands of the plane's intersections: an edge's curve and the plane, named as
/// upstream's solutions name them (the caches they build own their trees).
#[derive(Debug, Clone, Copy, PartialEq)]
enum Operand {
    /// The edge's curve, with its natural interval.
    Curve(Extent1d),
    /// The plane, with its natural domain.
    Plane(Extent2d),
}

impl SolutionObject for Operand {
    fn kind(&self) -> ObjectKind {
        match self {
            Self::Curve(_) => ObjectKind::Curve,
            Self::Plane(_) => ObjectKind::Surface,
        }
    }

    fn natural_interval(&self) -> Option<Extent1d> {
        match self {
            Self::Curve(ivl) => Some(*ivl),
            Self::Plane(_) => None,
        }
    }

    fn natural_uv_domain(&self) -> Option<Extent2d> {
        match self {
            Self::Plane(domain) => Some(*domain),
            Self::Curve(_) => None,
        }
    }
}

/// What one check of a body's convex blends keeps from section to section and from solver
/// to solver ([`band_cuts_the_part`] and [`band_breaks_through_the_part`] over every solver
/// of a session, #2648). The body does not change during a check: the trimmed-surface
/// caches and UV trims its solves make on it change no answer. So each face's trimmed
/// cache is made once ([`FaceCaches`]), the body's cache once ([`BrepQuery`]), and each
/// face's and edge's box once.
///
/// The body's cache holds the trees of its vertices, edges and faces, keyed by the body's
/// naming scheme and its proven geometry owners, each boxed with its tolerance of the
/// moment the cache was made (`make_brep_cache`). What the check's solves edit on the body
/// (a missing UV trim curve made on an edge use for a face's trimmed cache or a ray's
/// preparation, the edge's tolerance grown to fit it and its vertices' grown with it, and
/// the vertex uses' gap caches reset) adds, removes and renames nothing the trees hold;
/// but a tolerance grown after the cache was made is not in its boxes, and a solve's
/// candidate test reads the boxes (`SmTopologySolver::MayContainAnswer`), where a cache
/// made for that classification would hold it. So the check keeps, with the cache, every
/// vertex's, edge's and face's tolerance as the cache was made ([`tolerances`]) and before
/// each classification compares the body's to them: on any difference it makes the cache
/// again. Each classification is then that of a cache made for it.
///
/// The boxes keep each section's tests to the geometry within its reach, by construction
/// giving the answers of testing everything: a curve and a surface (or a plane) whose
/// boxes lie farther apart than the intersections' tolerance have no intersection, and a
/// plane crossing an edge lies on the edge, inside its box.
#[derive(Default)]
pub struct BandCheck {
    /// The faces' trimmed-surface caches made so far.
    caches: FaceCaches,
    /// The body's query with its cache, made at the first classification, and the
    /// tolerances its boxes were made with ([`tolerances`]).
    query: Option<(BrepQuery<'static>, Vec<u64>)>,
    /// Every face in the body's order, with the box of its surface over the domain its
    /// solves use (`None`: unbounded, or a kind without a certain box; always tested).
    faces: Option<Vec<(FaceId, Option<Extent3d>)>>,
    /// Every edge in the body's order, with the box of its curve over its interval
    /// (`None`: always tested).
    edges: Option<Vec<(EdgeId, Option<Extent3d>)>>,
}

impl core::fmt::Debug for BandCheck {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("BandCheck")
            .field("caches", &self.caches.len())
            .field("query", &self.query.is_some())
            .finish_non_exhaustive()
    }
}

impl BandCheck {
    /// A check that has made nothing yet.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Every face of `body` with its box, made the first time they are asked for.
    fn faces(&mut self, k: &Kernel, body: &Brep) -> Result<&[(FaceId, Option<Extent3d>)]> {
        if self.faces.is_none() {
            let mut out = Vec::new();
            for f in crate::topology::brep::p1::faces(body, None)? {
                let bound = match face::p2::surface(body, f)? {
                    Some(surface) => surface_bound(k, surface),
                    None => None,
                };
                out.push((f, bound));
            }
            self.faces = Some(out);
        }
        Ok(self.faces.as_deref().unwrap_or(&[]))
    }

    /// Every edge of `body` with its box, made the first time they are asked for.
    fn edges(&mut self, k: &Kernel, body: &Brep) -> Result<&[(EdgeId, Option<Extent3d>)]> {
        if self.edges.is_none() {
            let mut out = Vec::new();
            for e in crate::topology::brep::p1::edges(body, None)? {
                let bound = match edge::evaluation_curve(body, e)? {
                    Some(curve) => curve_bound(k, &curve, &edge::interval(body, e)?),
                    None => None,
                };
                out.push((e, bound));
            }
            self.edges = Some(out);
        }
        Ok(self.edges.as_deref().unwrap_or(&[]))
    }

    /// Solves `surface` (face `f`'s, over `domain`) against `section` through the face's
    /// trimmed cache, kept from solve to solve.
    #[allow(clippy::too_many_arguments)]
    fn face_solve(
        &mut self,
        k: &Kernel,
        body: &mut Brep,
        f: FaceId,
        surface: &Surface,
        domain: &Extent2d,
        section: &Curve,
        interval: &Extent1d,
        tolerance: f64,
    ) -> Result<crate::geom::solution_array::SolutionArray<CacheObject>> {
        k.work_counters.add_band_face_solve();
        let before = self.caches.len();
        let mut services =
            BrepFaceServices::with_caches(&mut *body, core::mem::take(&mut self.caches));
        // the face's surface solves through its face (its trimmed cache), as upstream's
        // `GlobalCurveIntersect` reaches the owner's cache
        let solutions = surface.global_curve_intersect_with_context::<CacheObject>(
            k,
            [None, Some(CacheObject::Surface(f))],
            None,
            &mut services,
            domain,
            section,
            interval,
            tolerance,
        );
        self.caches = services.into_caches();
        for _ in before..self.caches.len() {
            k.work_counters.add_band_face_cache();
        }
        solutions
    }

    /// `point`'s classification against `body` (`SmBrep::Point3DClassify`, with the boundary
    /// intersections) through the body's cache, made at the first classification and again
    /// when a tolerance of the body differs from the ones the cache was made with.
    fn classify(
        &mut self,
        k: &Kernel,
        body: &mut Brep,
        point: Point3,
        zone: f64,
    ) -> Result<PointClassification<EntityRef>> {
        let now = tolerances(body)?;
        if self
            .query
            .as_ref()
            .is_some_and(|(_, made_with)| *made_with != now)
        {
            self.query = None;
        }
        let query = if let Some((query, _)) = self.query.as_ref() {
            query
        } else {
            let query = BrepQuery::of_brep(body);
            query.cache(k, body)?;
            k.work_counters.add_band_brep_cache();
            &self.query.insert((query, now)).0
        };
        query.point_3d_classify(k, body, point, zone, true)
    }
}

/// The tolerance of every vertex, edge and face of `body`, as bits, in the order the body's
/// cache boxes them: what [`BandCheck`] compares before a classification to know whether
/// the cache's boxes still carry the body's tolerances (bits, so a NaN compares equal to
/// itself and does not make the cache again at every classification).
fn tolerances(body: &Brep) -> Result<Vec<u64>> {
    let mut out = Vec::new();
    for v in crate::topology::brep::p1::vertices(body, None)? {
        out.push(body.vertex(v)?.tolerance.to_bits());
    }
    for e in crate::topology::brep::p1::edges(body, None)? {
        out.push(body.edge(e)?.tolerance.to_bits());
    }
    for f in crate::topology::brep::p1::faces(body, None)? {
        out.push(body.face(f)?.tolerance.to_bits());
    }
    Ok(out)
}

/// A box that bounds something: finite and not inside out.
fn usable(b: Extent3d) -> Option<Extent3d> {
    let finite = [b.min, b.max]
        .iter()
        .all(|p| p.x.is_finite() && p.y.is_finite() && p.z.is_finite());
    (finite && !b.has_negative_volume()).then_some(b)
}

/// The box of `surface` over its natural domain (the domain the check's solves use), for
/// the kinds whose box is the hull of a control net that holds the surface (the B-spline
/// override of `CalculateBoundingBox`, unexpanded); `None` for the others or a box that
/// cannot be made: such a face is always tested.
fn surface_bound(k: &Kernel, surface: &Surface) -> Option<Extent3d> {
    if !matches!(
        surface.kind(),
        SurfaceKind::BSpline(_)
            | SurfaceKind::Plane(_)
            | SurfaceKind::Cone(_)
            | SurfaceKind::Sphere(_)
            | SurfaceKind::Torus(_)
            | SurfaceKind::ProceduralBlend(_)
            | SurfaceKind::SurfOfRevolution(_)
            | SurfaceKind::SurfOfExtrusion(_)
    ) {
        return None;
    }
    let surface = surface.clone_for_evaluation().ok()?;
    let domain = surface.natural_uv_domain_with_nurb(k).ok()?;
    let boxes = surface
        .calculate_bounding_box(k, &domain, [true, false, false], false)
        .ok()?;
    usable(boxes.normal_box?)
}

/// The box of `curve` over `interval`, for the kinds whose box is the hull of a control
/// polygon that holds the curve (a B-spline, conic or line); `None` for the others or a box
/// that cannot be made.
fn curve_bound(k: &Kernel, curve: &Curve, interval: &Extent1d) -> Option<Extent3d> {
    if !matches!(
        curve.kind(),
        CurveKind::BSpline(_)
            | CurveKind::Ellipse(_)
            | CurveKind::Line(_)
            | CurveKind::Hyperbola(_)
            | CurveKind::Parabola(_)
    ) {
        return None;
    }
    let req = BoundingBoxRequest {
        normal: true,
        pseudo: false,
        polar: false,
    };
    usable(
        curve
            .calculate_bounding_box(k, interval, req, true)
            .ok()?
            .normal?,
    )
}

/// The distance from `p` to the box `b` (0 inside it).
fn distance_to_box(b: &Extent3d, p: Point3) -> f64 {
    let gap = |v: f64, lo: f64, hi: f64| {
        if v < lo {
            lo - v
        } else if v > hi {
            v - hi
        } else {
            0.0
        }
    };
    Vec3::new(
        gap(p.x, b.min.x, b.max.x),
        gap(p.y, b.min.y, b.max.y),
        gap(p.z, b.min.z, b.max.z),
    )
    .length()
}

/// A station's plane, through `origin` square to the unit `normal`, and the body's edges
/// crossing it (`sm_EdgesCrossingPlane`), each edge's crossings found the first time a
/// section asks for them.
struct Station {
    /// The plane.
    plane: Surface,
    /// Its natural domain.
    domain: Extent2d,
    /// The edges whose boxes reach the plane (within twice the body's tolerance, the
    /// intersections' tolerance), in the body's order, with their boxes.
    near: Vec<(EdgeId, Option<Extent3d>)>,
    /// The crossings found so far, by edge.
    crossings: BTreeMap<EdgeId, Vec<Point3>>,
}

impl Station {
    /// The plane of the station at `origin` with the unit `normal` (the edge's tangent).
    fn new(
        k: &Kernel,
        check: &mut BandCheck,
        body: &Brep,
        origin: Point3,
        normal: Vec3,
    ) -> Result<Self> {
        let x = if normal.x.abs() < 0.9 {
            Vec3::new(1.0, 0.0, 0.0).cross(normal)
        } else {
            Vec3::new(0.0, 1.0, 0.0).cross(normal)
        };
        let x = x.unitize()?;
        let y = normal.cross(x);
        let plane = Surface::from(Plane::create_canonical(&Axis2Placement::new(origin, x, y))?);
        let domain = plane.natural_uv_domain_with_nurb(k)?;
        let reach = 2.0 * body.tolerance();
        let near = check
            .edges(k, body)?
            .iter()
            .filter(|(_, bound)| {
                bound.is_none_or(|b| !(b.distance_to_plane(origin, normal).abs() > reach))
            })
            .copied()
            .collect();
        Ok(Self {
            plane,
            domain,
            near,
            crossings: BTreeMap::new(),
        })
    }

    /// `sm_EdgesCrossingPlane` for one edge: the points where its curve crosses the plane,
    /// the curve intersected with the plane; an edge lying in the plane gives its ends there.
    fn crossings_of(&mut self, k: &Kernel, body: &Brep, e: EdgeId) -> Result<&[Point3]> {
        if !self.crossings.contains_key(&e) {
            let mut crossings = Vec::new();
            if let Some(curve) = edge::evaluation_curve(body, e)? {
                k.work_counters.add_band_edge_solve();
                let interval = edge::interval(body, e)?;
                let solutions = self.plane.global_curve_intersect(
                    k,
                    [
                        Operand::Curve(curve.natural_interval()),
                        Operand::Plane(self.domain),
                    ],
                    &self.domain,
                    &curve,
                    &interval,
                    body.tolerance(),
                )?;
                for sol in solutions.as_slice() {
                    let t = sol
                        .start
                        .parameters
                        .first()
                        .copied()
                        .ok_or_else(|| SolidError::null_pointer(OP))?;
                    crossings.push(curve.evaluate_point(k, t)?);
                    if sol.solution_type == SolutionType::RangeOfValues {
                        let t = sol
                            .end
                            .parameters
                            .first()
                            .copied()
                            .ok_or_else(|| SolidError::null_pointer(OP))?;
                        crossings.push(curve.evaluate_point(k, t)?);
                    }
                }
            }
            self.crossings.insert(e, crossings);
        }
        Ok(self.crossings.get(&e).map_or(&[][..], Vec::as_slice))
    }
}

/// The box of a section over its interval; `None` (every face tested) where it cannot be
/// made.
fn section_bound(k: &Kernel, section: &Curve, interval: &Extent1d) -> Option<Extent3d> {
    curve_bound(k, section, interval)
}

/// `sm_KeepsMaterialBehind`: whether the part keeps material behind the section of a
/// convex blend at a station, the region between the section and the ball's centre
/// `centre` (the ball's sector between the rails for a round, the triangle of the centre
/// and the rails for a chamfer's `chord`).
///
/// Kept where a body edge crosses the station's plane inside that region (`crossings`),
/// or the section (as an exact curve) crosses a face's surface inside the face. With
/// neither, one point of the section classified against the body says whether it lies
/// inside the part; outside, the blend takes all of the part between its rails there.
///
/// Only the geometry within reach of the section is tested (#2648), by construction with
/// the answer of testing all of it: a crossing counts only within the ball's radius of its
/// centre (a round's) or on the station's plane (a chamfer's), so an edge whose box lies
/// farther from the centre, or from the plane, has none that counts; and a face whose
/// surface's box lies farther from the section's box than twice the intersections'
/// tolerance has no intersection with it.
#[allow(clippy::too_many_arguments)]
fn keeps_material_behind(
    k: &Kernel,
    body: &mut Brep,
    check: &mut BandCheck,
    station: &mut Station,
    centre: Point3,
    rails: [Point3; 2],
    chord: bool,
) -> Result<bool> {
    let [rail1, rail2] = rails;
    let a = rail1 - centre;
    let b = rail2 - centre;
    let r2 = smos_max(a.length_squared(), b.length_squared());
    let (aa, ab, bb) = (a.dot(a), a.dot(b), b.dot(b));
    let det = aa * bb - ab * ab;
    if !(det > EFF_ZERO * aa * bb) {
        return Ok(true); // rails opposite or together: no region to judge
    }
    let tol = 1.0e-9 * r2.sqrt();
    // a round's crossing counts within `r2.sqrt() + tol` of the centre: the edges whose boxes
    // lie farther (by more than the body's tolerance, as a crossing lies on its curve) have
    // none that counts; a chamfer's counts anywhere on the plane over its triangle
    let reach = (!chord).then(|| r2.sqrt() + tol + body.tolerance());
    let mut i = 0;
    while let Some(&(e, bound)) = station.near.get(i) {
        i += 1;
        if let (Some(reach), Some(b)) = (reach, bound) {
            if distance_to_box(&b, centre) > reach {
                continue;
            }
        }
        k.check_cancelled(OP)?;
        for &crossing in station.crossings_of(k, body, e)? {
            let x = crossing - centre;
            let (xa, xb) = (x.dot(a), x.dot(b));
            let alpha = (xa * bb - xb * ab) / det;
            let beta = (xb * aa - xa * ab) / det;
            if alpha < -1.0e-12 || beta < -1.0e-12 {
                continue;
            }
            let inside = if chord {
                alpha + beta <= 1.0 + 1.0e-12
            } else {
                x.length() <= r2.sqrt() + tol
            };
            if inside {
                return Ok(true);
            }
        }
    }

    // the section as an exact curve: the chord, or the circle's arc from the first rail to the
    // second the short way
    let (section, middle) = if chord {
        let line = Line::create_line_segment(3, rail1, rail2, None)?;
        (Curve::from(line), (rail1 + rail2) * 0.5)
    } else {
        let rho = aa.sqrt();
        let e1 = a / rho;
        let e2 = (b - a * (ab / aa)).unitize()?;
        let phi = b.dot(e2).atan2(b.dot(e1));
        let circle = Ellipse::new_circle(
            k,
            centre,
            e1,
            e2,
            Extent1d::from_min_max(0.0, rad_to_deg(phi)),
            rho,
            &EllipseOptions {
                dim: 3,
                ..EllipseOptions::default()
            },
        )?;
        let middle = centre + e1 * (rho * (0.5 * phi).cos()) + e2 * (rho * (0.5 * phi).sin());
        (Curve::from(circle), middle)
    };
    let interval = section.natural_interval();
    let ends = body.tolerance();
    let zone = brep_zone_tol_3d(k, body, SmallTopoType::default());
    let reach = section_bound(k, &section, &interval);
    let mut i = 0;
    while let Some(&(f, bound)) = check.faces(k, body)?.get(i) {
        i += 1;
        if let (Some(s), Some(b)) = (reach, bound) {
            if s.are_disjoint(&b, 2.0 * ends) {
                continue; // no point of the section within the tolerance of the surface
            }
        }
        k.check_cancelled(OP)?;
        // the crossings inside the face's surface's domain, judged after the borrow ends
        let mut on_face: Vec<Point2> = Vec::new();
        {
            let Some(surface) = face::p2::surface(body, f)? else {
                continue;
            };
            let surface = surface.clone_for_evaluation()?;
            let domain = surface.natural_uv_domain_with_nurb(k)?;
            let solutions =
                check.face_solve(k, body, f, &surface, &domain, &section, &interval, ends)?;
            for sol in solutions.as_slice() {
                if sol.solution_type == SolutionType::RangeOfValues {
                    return Ok(true); // the section runs along the face's surface
                }
                let params = &sol.start.parameters;
                let (Some(&t), Some(&u), Some(&v)) = (params.first(), params.get(1), params.get(2))
                else {
                    return Err(SolidError::null_pointer(OP));
                };
                let at = section.evaluate_point(k, t)?;
                if at.distance(rail1) <= ends || at.distance(rail2) <= ends {
                    continue; // the section's own ends, on the blended faces
                }
                on_face.push(Point2::new(u, v));
            }
        }
        for uv in on_face {
            let mut class = PointClassification::<EntityRef>::new(zone);
            crate::query::face::p2::point_classify(
                k,
                body,
                f,
                uv,
                &mut class,
                crate::query::face::p2::PointClassifyOptions {
                    src_zone_tol_3d: zone,
                    do_boundary_intersections: true,
                    use_uv_classification: true,
                    ..crate::query::face::p2::PointClassifyOptions::default()
                },
            )?;
            if class.point_class() == PointClassType::Face {
                return Ok(true); // the section crosses the part's boundary inside a face
            }
        }
    }
    // no crossing: the section between its rails is wholly in the part or wholly outside it
    let class = check.classify(k, body, middle, zone)?;
    if class.point_class() != PointClassType::Region {
        return Ok(true);
    }
    if let Some(EntityRef::Region(region)) = class.object() {
        if body.region(region)?.is_void() == 0 {
            return Ok(true);
        }
    }
    Ok(false)
}

/// One section judged by `keeps`.
#[derive(Debug, Clone, Copy, Default)]
struct Section {
    /// `rbJudged`: the section was found.
    found: bool,
    /// `rbKeeps`: material is kept behind it.
    keeps: bool,
    /// `rdSetbackHere`: the distance from the edge to the first rail.
    setback: f64,
}

/// The `Keeps` lambda of `BandCutsThePart`: at the edge's `t` with the ball's radius
/// `trial` (0: the solver's own, restored to `own` after), whether material is kept behind
/// the section; `near` is the guess, then this section if found.
#[allow(clippy::too_many_arguments)]
fn keeps(
    k: &Kernel,
    target: &mut Brep,
    solver: &mut SolverObject,
    check: &mut BandCheck,
    edge_curve: &Curve,
    chord: bool,
    own: f64,
    t: f64,
    trial: f64,
    judge: Judge<'_>,
    near: &mut Option<TsectPnt>,
) -> Result<Section> {
    k.check_cancelled(OP)?;
    let mut out = Section {
        found: false,
        keeps: true,
        setback: 0.0,
    };
    if trial > 0.0 && !set_ball_radius(solver, trial)? {
        return Ok(out);
    }
    let evaluated = (|| -> Result<Option<(TsectPnt, [Point3; 3])>> {
        let Some(pnt) = ball_section_at(k, target, solver, edge_curve, t, near.as_ref())? else {
            return Ok(None);
        };
        let uvs = [pnt.uv_pos(SrfIndex::First), pnt.uv_pos(SrfIndex::Second)];
        let base = solver.base_mut();
        let mut at = [Point3::UNDEF; 3];
        {
            let slot = base
                .surface_mut(0)?
                .ok_or_else(|| SolidError::null_pointer(OP))?;
            at[0] = slot.surface.surface.evaluate_point_extending(k, uvs[0])?;
        }
        for (rail, uv) in uvs.into_iter().enumerate() {
            let slot = base
                .surface(rail)?
                .ok_or_else(|| SolidError::null_pointer(OP))?;
            let p = slot.surface.surface.base_surface().evaluate_point(k, uv)?;
            if let Some(dst) = at.get_mut(rail + 1) {
                *dst = p;
            }
        }
        // the rails past a base surface's domain lie on its extension, as the offset evaluates it
        for (rail, uv) in uvs.into_iter().enumerate() {
            let slot = base
                .surface(rail)?
                .ok_or_else(|| SolidError::null_pointer(OP))?;
            let offset = &slot.surface.surface;
            if let Some(extended) = offset.extended_base_surface() {
                if !offset
                    .base_surface()
                    .natural_uv_domain_with_nurb(k)?
                    .contains_point(uv, 0.0)
                    && extended
                        .natural_uv_domain_with_nurb(k)?
                        .contains_point(uv, 0.0)
                {
                    let p = extended.evaluate_point(k, uv)?;
                    if let Some(dst) = at.get_mut(rail + 1) {
                        *dst = p;
                    }
                }
            }
        }
        Ok(Some((pnt, at)))
    })();
    if trial > 0.0 {
        // the solver's own radius back, whatever the solve did
        let _ = set_ball_radius(solver, own);
    }
    let Some((pnt, [centre, rail1, rail2])) = evaluated? else {
        return Ok(out);
    };
    *near = Some(pnt);
    out.found = true;
    let p = edge_curve.evaluate_point(k, t)?;
    out.setback = p.distance(rail1);
    k.work_counters.add_band_section();
    out.keeps = match judge {
        Judge::KeepsMaterial(station) => {
            keeps_material_behind(k, target, check, station, centre, [rail1, rail2], chord)?
        }
        Judge::StaysInside(faces) => {
            !section_breaks_through(k, target, check, centre, [rail1, rail2], chord, faces)?
        }
    };
    Ok(out)
}

/// What `keeps` judges at a section: material kept behind it (`sm_KeepsMaterialBehind`,
/// with the body's edges crossing the station's plane), or (WR-707, not the fork's) the
/// section not breaking through the part (`section_breaks_through`).
enum Judge<'a> {
    /// `sm_KeepsMaterialBehind` with the station's plane and its crossings.
    KeepsMaterial(&'a mut Station),
    /// `section_breaks_through` false, with these faces.
    StaysInside(&'a ThroughFaces),
}

/// The faces `section_breaks_through` tells apart (WR-707, not the fork's).
struct ThroughFaces {
    /// The faces at the blended edge's ends, which the corners trim: not judged.
    skip: BTreeSet<FaceId>,
    /// The blended faces and every face sharing an edge with one: the faces the blend
    /// reaches over the blended faces' boundaries, which it is trimmed to.
    reached: BTreeSet<FaceId>,
}

/// Not the fork's (WR-707, #2299): whether the section of a convex blend at a station,
/// between its rails (the ball's arc, or the chamfer's chord, as an exact curve), breaks
/// through the part's boundary: it leaves the part through a face and comes back into it
/// through another (or the same) face, so the part between those crossings is a hollow
/// behind a wall that the blend opens (a round on a wall thinner than it reaches into the
/// hollow behind, #2299).
///
/// The section is split at every point where it meets the body's boundary (crossings
/// within the body's tolerance of its own ends are its rails), and each piece between two
/// such points is classified by its middle. The section breaks through where a piece
/// lies outside the part's material, one of its ends crosses a face inside the face (not
/// at an edge), neither end is on a face at the blended edge's ends (which the corners
/// trim), and one end is on a face that shares no edge with the blended faces: an inner
/// face of the part, bounding material behind a wall, that the blend does not reach over
/// the blended faces' boundaries; and the material beside the piece runs to a rail on one
/// side (WR-707 addendum, 2026-11-03, #2328). That material is the wall, part of what the
/// blend rounds. A piece between islands of material, with void on either side of each
/// along the section (a thin rib across a slot), is not behind a wall: the blend runs past
/// the islands and is trimmed to their faces.
///
/// So a piece outside the part that reaches a rail is not a break-through: there the
/// blend runs off its face (over a cliff, into a hole or an opening) and is trimmed to the
/// faces it runs onto. Nor is a piece through a hole or an opening in a blended face: its
/// walls bound that face. The construction builds both to their closed-form volumes
/// (SolidCraft rig cases the fork builds). It builds a blind pocket or a cavity under a
/// blended face wrong, and those are refused (WR-707 addendum, 2026-11-02).
fn section_breaks_through(
    k: &Kernel,
    body: &mut Brep,
    check: &mut BandCheck,
    centre: Point3,
    rails: [Point3; 2],
    chord: bool,
    faces: &ThroughFaces,
) -> Result<bool> {
    let [rail1, rail2] = rails;
    let a = rail1 - centre;
    let b = rail2 - centre;
    let (aa, ab, bb) = (a.dot(a), a.dot(b), b.dot(b));
    let det = aa * bb - ab * ab;
    if !(det > EFF_ZERO * aa * bb) {
        return Ok(false); // rails opposite or together: no section to judge
    }
    let section = if chord {
        Curve::from(Line::create_line_segment(3, rail1, rail2, None)?)
    } else {
        let rho = aa.sqrt();
        let e1 = a / rho;
        let e2 = (b - a * (ab / aa)).unitize()?;
        let phi = b.dot(e2).atan2(b.dot(e1));
        Curve::from(Ellipse::new_circle(
            k,
            centre,
            e1,
            e2,
            Extent1d::from_min_max(0.0, rad_to_deg(phi)),
            rho,
            &EllipseOptions {
                dim: 3,
                ..EllipseOptions::default()
            },
        )?)
    };
    let interval = section.natural_interval();
    let ends = body.tolerance();
    let zone = brep_zone_tol_3d(k, body, SmallTopoType::default());
    // where the section meets the boundary of face `f`: its parameter, whether it crosses
    // `f` inside it (not on its boundary), whether `f` is judged and whether the blend
    // reaches `f`
    let mut meets: Vec<Meet> = Vec::new();
    // (as `keeps_material_behind`: a face whose surface's box lies apart from the section's
    // has no meet)
    let reach = section_bound(k, &section, &interval);
    let mut i = 0;
    while let Some(&(f, bound)) = check.faces(k, body)?.get(i) {
        i += 1;
        if let (Some(s), Some(b)) = (reach, bound) {
            if s.are_disjoint(&b, 2.0 * ends) {
                continue;
            }
        }
        k.check_cancelled(OP)?;
        let mut on_face: Vec<(f64, Point2)> = Vec::new();
        {
            let Some(surface) = face::p2::surface(body, f)? else {
                continue;
            };
            let surface = surface.clone_for_evaluation()?;
            let domain = surface.natural_uv_domain_with_nurb(k)?;
            let solutions =
                check.face_solve(k, body, f, &surface, &domain, &section, &interval, ends)?;
            for sol in solutions.as_slice() {
                let params = &sol.start.parameters;
                let (Some(&t), Some(&u), Some(&v)) = (params.first(), params.get(1), params.get(2))
                else {
                    return Err(SolidError::null_pointer(OP));
                };
                let at = section.evaluate_point(k, t)?;
                if at.distance(rail1) <= ends || at.distance(rail2) <= ends {
                    continue; // the section's own ends, on the blended faces
                }
                on_face.push((t, Point2::new(u, v)));
            }
        }
        for (t, uv) in on_face {
            let mut class = PointClassification::<EntityRef>::new(zone);
            crate::query::face::p2::point_classify(
                k,
                body,
                f,
                uv,
                &mut class,
                crate::query::face::p2::PointClassifyOptions {
                    src_zone_tol_3d: zone,
                    do_boundary_intersections: true,
                    use_uv_classification: true,
                    ..crate::query::face::p2::PointClassifyOptions::default()
                },
            )?;
            let inside = match class.point_class() {
                PointClassType::Unknown => continue, // outside the face
                PointClassType::Face => true,
                _ => false, // on the face's boundary
            };
            meets.push(Meet {
                t,
                inside,
                judged: !faces.skip.contains(&f),
                reached: faces.reached.contains(&f),
            });
        }
    }
    meets.sort_by(|x, y| x.t.total_cmp(&y.t));
    // the bounds of the section's pieces: its rails and the meets between them; piece `i`
    // runs from `at[i]` to `at[i + 1]`, so meet `j` ends piece `j` and starts piece `j + 1`
    let mut at = Vec::with_capacity(meets.len() + 2);
    at.push(interval.min());
    at.extend(meets.iter().map(|m| m.t));
    at.push(interval.max());
    let mut pieces = Pieces {
        section: &section,
        at: &at,
        classes: vec![None; at.len() - 1],
        zone,
    };
    for (j, pair) in meets.windows(2).enumerate() {
        let [m0, m1] = pair else {
            continue;
        };
        if !(m0.judged && m1.judged) || !(m0.inside || m1.inside) || (m0.reached && m1.reached) {
            continue;
        }
        let piece = j + 1;
        if pieces.class(k, body, check, piece)? != Piece::Void {
            continue;
        }
        // the wall: the material beside the piece runs to a rail
        if pieces.material_to_rail(k, body, check, (0..piece).rev())?
            || pieces.material_to_rail(k, body, check, piece + 1..at.len() - 1)?
        {
            return Ok(true); // a hollow behind a wall: the section breaks through
        }
    }
    Ok(false)
}

/// The pieces of a section between the points where it meets the part's boundary
/// (`section_breaks_through`), each classified once, by its middle.
struct Pieces<'a> {
    /// The section.
    section: &'a Curve,
    /// The pieces' bounds, in the section's parameter.
    at: &'a [f64],
    /// Each piece's class, once classified.
    classes: Vec<Option<Piece>>,
    /// The body's zone tolerance.
    zone: f64,
}

impl Pieces<'_> {
    /// The class of piece `i`.
    fn class(
        &mut self,
        k: &Kernel,
        body: &mut Brep,
        check: &mut BandCheck,
        i: usize,
    ) -> Result<Piece> {
        if let Some(Some(known)) = self.classes.get(i) {
            return Ok(*known);
        }
        let (Some(&t0), Some(&t1)) = (self.at.get(i), self.at.get(i + 1)) else {
            return Err(SolidError::null_pointer(OP));
        };
        let class = if !(t1 - t0 > 0.0) {
            Piece::Point
        } else if point_in_material(
            k,
            body,
            check,
            self.section.evaluate_point(k, 0.5 * (t0 + t1))?,
            self.zone,
        )? {
            Piece::Material
        } else {
            Piece::Void
        };
        if let Some(slot) = self.classes.get_mut(i) {
            *slot = Some(class);
        }
        Ok(class)
    }

    /// Whether the pieces `run`, taken from a void piece out to a rail, all lie in the
    /// material (a piece of no length does not interrupt them): the material beside the
    /// void piece is the part the blend rounds, reaching the blended face at its rail, not
    /// an island with void on both sides that the blend runs past and is trimmed to.
    fn material_to_rail(
        &mut self,
        k: &Kernel,
        body: &mut Brep,
        check: &mut BandCheck,
        run: impl Iterator<Item = usize>,
    ) -> Result<bool> {
        for i in run {
            if self.class(k, body, check, i)? == Piece::Void {
                return Ok(false);
            }
        }
        Ok(true)
    }
}

/// What a piece of a section between two meets is (`section_breaks_through`).
#[derive(Clone, Copy, PartialEq, Eq)]
enum Piece {
    /// In the part's material.
    Material,
    /// Outside it.
    Void,
    /// Of no length: one point, met on two faces.
    Point,
}

/// Where a section meets the boundary of a face (`section_breaks_through`).
#[derive(Clone, Copy)]
struct Meet {
    /// The section's parameter there.
    t: f64,
    /// Inside the face, not on its boundary.
    inside: bool,
    /// The face is judged: not one at the blended edge's ends.
    judged: bool,
    /// The face is one the blend reaches ([`ThroughFaces::reached`]).
    reached: bool,
}

/// Whether `point` lies in the part's material, inside a region that is not void (the
/// final classification of `sm_KeepsMaterialBehind`, with a point on the boundary, within
/// the zone tolerance of a face, edge or vertex, not counted in: a piece of a section
/// between two crossings of faces that lies within the tolerance of the boundary is too
/// shallow to keep material).
fn point_in_material(
    k: &Kernel,
    body: &mut Brep,
    check: &mut BandCheck,
    point: Point3,
    zone: f64,
) -> Result<bool> {
    let class = check.classify(k, body, point, zone)?;
    if class.point_class() != PointClassType::Region {
        return Ok(false);
    }
    if let Some(EntityRef::Region(region)) = class.object() {
        return Ok(body.region(region)?.is_void() == 0);
    }
    Ok(false)
}

/// `SmFilletSolver::BandCutsThePart`: whether this convex blend (one that removes
/// material) would take all of the part between its rails somewhere along its edge,
/// cutting it into pieces (a round or chamfer through the bridge over a slot); and where it
/// would, the largest size that keeps the part whole. `None` where it does not cut.
///
/// At the middles of 32 equal steps of the edge (as [`super::ball_regular::ball_is_regular`]),
/// with the body's edges crossing each section's plane found once. A section where no
/// material lies behind the blend's section is a cut where sections either side of its run
/// along the edge keep material. The limit is then found at each cutting section by
/// bisection on the ball's radius, least over the sections. A constant ball only.
///
/// `target` is the body blended; it holds the solver's edge use.
///
/// # Errors
///
/// Those of the section solves, the intersections, the classifications and the
/// evaluations; `SM_ERR_NULL_POINTER` for an edge without a curve.
pub fn band_cuts_the_part(
    k: &Kernel,
    target: &mut Brep,
    solver: &mut SolverObject,
    check: &mut BandCheck,
) -> Result<Option<BandCut>> {
    let base = solver.base();
    if base.surface(0)?.is_none() || base.surface(1)?.is_none() {
        return Ok(None);
    }
    let Some(edgeuse) = base.edgeuse(0)? else {
        return Ok(None);
    };
    let edge_id = eu::edge(target, edgeuse)?.ok_or_else(|| SolidError::null_pointer(OP))?;
    let curve = edge::evaluation_curve(target, edge_id)?
        .ok_or_else(|| SolidError::null_pointer(OP))?
        .clone();
    let ivl = edge::interval(target, edge_id)?;
    let (rolls, on_edge) = rolls_a_ball(solver);
    if !rolls || on_edge {
        return Ok(None); // a constant ball only
    }
    let probe = TsectPnt::new();
    let own = ball_law(k, solver, &probe)?.r[0];
    let chord = solver.base().fsg.as_ref().is_some_and(|fsg| {
        fsg.fillet_surface_generator_type() == FilletSurfaceGeneratorType::Linear
    });

    // each section: its parameter, station (its plane and crossings) and section found near,
    // and whether it keeps material behind it (0), takes all of the part between the rails
    // (1), or is not judged (2)
    let mut ts: Vec<f64> = Vec::new();
    let mut stations: Vec<Option<Station>> = Vec::new();
    let mut nears: Vec<Option<TsectPnt>> = Vec::new();
    let mut states: Vec<u8> = Vec::new();
    let mut near: Option<TsectPnt> = None;
    for j in 0..STEPS {
        let t = ivl.evaluate((j as f64 + 0.5) / STEPS as f64);
        ts.push(t);
        let pv = curve.evaluate(k, t, 1, true, false)?;
        let (Some(&point), Some(&tangent)) = (pv.first(), pv.get(1)) else {
            return Err(SolidError::null_pointer(OP));
        };
        if !(tangent.length() > EFF_ZERO) {
            stations.push(None);
            nears.push(None);
            states.push(2);
            continue;
        }
        let normal = tangent.unitize()?;
        let mut station = Station::new(k, check, target, point, normal)?;
        let judged = keeps(
            k,
            target,
            solver,
            check,
            &curve,
            chord,
            own,
            t,
            0.0,
            Judge::KeepsMaterial(&mut station),
            &mut near,
        )?;
        stations.push(Some(station));
        nears.push(near.clone());
        states.push(if judged.found {
            u8::from(!judged.keeps)
        } else {
            2
        });
    }

    // a run of sections that take all of the part cuts it in two only where sections either side
    // of the run along the edge keep material: a run that reaches an end of the edge (a tip, where
    // the part runs out thinner than the band) blends that end away and leaves the part whole
    let state = |i: usize| states.get(i).copied().unwrap_or(2);
    let mut cut: Option<BandCut> = None;
    let mut j = 0;
    while j < STEPS {
        if state(j) != 1 {
            j += 1;
            continue;
        }
        let first = j;
        let mut last = j;
        while last + 1 < STEPS && state(last + 1) != 0 {
            last += 1;
        }
        j = last + 1;
        let kept_before = (0..first).any(|i| state(i) == 0);
        let kept_after = (last + 1..STEPS).any(|i| state(i) == 0);
        if !kept_before || !kept_after {
            continue;
        }
        for m in first..=last {
            if state(m) != 1 {
                continue;
            }
            let t = ts.get(m).copied().unwrap_or(0.0);
            // (a section taking all of the part was judged at its station)
            let Some(station) = stations.get_mut(m).and_then(Option::as_mut) else {
                continue;
            };
            // the largest radius that keeps material behind its section here
            let mut at = nears.get(m).cloned().flatten();
            let radius = own;
            let (mut lo, mut hi, mut lo_setback) = (0.0_f64, radius, 0.0);
            let mut tr = 0.5 * radius;
            while tr > 1.0e-6 * radius {
                let judged = keeps(
                    k,
                    target,
                    solver,
                    check,
                    &curve,
                    chord,
                    own,
                    t,
                    tr,
                    Judge::KeepsMaterial(station),
                    &mut at,
                )?;
                if judged.found && judged.keeps {
                    lo = tr;
                    lo_setback = judged.setback;
                    break;
                }
                hi = tr;
                tr *= 0.5;
            }
            let mut i = 0;
            while i < BISECTIONS && lo > 0.0 {
                let mid = 0.5 * (lo + hi);
                let judged = keeps(
                    k,
                    target,
                    solver,
                    check,
                    &curve,
                    chord,
                    own,
                    t,
                    mid,
                    Judge::KeepsMaterial(station),
                    &mut at,
                )?;
                if judged.found && judged.keeps {
                    lo = mid;
                    lo_setback = judged.setback;
                } else {
                    hi = mid;
                }
                i += 1;
            }
            if cut.is_none_or(|c| lo < c.radius) {
                cut = Some(BandCut {
                    radius: lo,
                    setback: lo_setback,
                });
            }
        }
    }
    Ok(cut)
}

/// Not the fork's (WR-707, #2299; `fixes.fillet_refuse_rounds_breaking_through` with
/// `fixes.fillet_gn_section_generator`): whether this convex blend's section breaks
/// through the part somewhere along its edge (it leaves the part and comes back into it
/// through an inner face, one sharing no edge with the blended faces: a hollow behind a
/// wall, `section_breaks_through`), and where it would, the largest size that keeps every
/// section from breaking through. At the same 32 sections as [`band_cuts_the_part`]; faces
/// at the blended edge's ends (the corners trim the blend there) are not judged. The limit at
/// each crossing section is found as the cut's is (halving from the ball's radius, then
/// 40 bisections), least over the sections. A constant ball only. `None` where no
/// section breaks through.
///
/// # Errors
///
/// As [`band_cuts_the_part`].
pub fn band_breaks_through_the_part(
    k: &Kernel,
    target: &mut Brep,
    solver: &mut SolverObject,
    check: &mut BandCheck,
) -> Result<Option<BandCut>> {
    let base = solver.base();
    if base.surface(0)?.is_none() || base.surface(1)?.is_none() {
        return Ok(None);
    }
    let Some(edgeuse) = base.edgeuse(0)? else {
        return Ok(None);
    };
    let edge_id = eu::edge(target, edgeuse)?.ok_or_else(|| SolidError::null_pointer(OP))?;
    let curve =
        edge::evaluation_curve(target, edge_id)?.ok_or_else(|| SolidError::null_pointer(OP))?;
    let ivl = edge::interval(target, edge_id)?;
    let (rolls, on_edge) = rolls_a_ball(solver);
    if !rolls || on_edge {
        return Ok(None); // a constant ball only
    }
    let probe = TsectPnt::new();
    let own = ball_law(k, solver, &probe)?.r[0];
    let chord = solver.base().fsg.as_ref().is_some_and(|fsg| {
        fsg.fillet_surface_generator_type() == FilletSurfaceGeneratorType::Linear
    });
    // the faces at the edge's ends
    let mut skip = BTreeSet::new();
    let v0 = edge::vertex(target, edge_id)?;
    for v in [v0, edge::other_vertex(target, edge_id, v0)?]
        .into_iter()
        .flatten()
    {
        skip.extend(vertex::faces(target, v)?);
    }
    // the blended faces and the faces sharing an edge with them
    let mut reached = BTreeSet::new();
    for f in edge::faces(target, edge_id, None)? {
        reached.insert(f);
        for e in face::p4::edges(target, f, None, false, None)? {
            reached.extend(edge::faces(target, e, None)?);
        }
    }
    let faces = ThroughFaces { skip, reached };
    let mut out: Option<BandCut> = None;
    let mut near: Option<TsectPnt> = None;
    for j in 0..STEPS {
        let t = ivl.evaluate((j as f64 + 0.5) / STEPS as f64);
        let judged = keeps(
            k,
            target,
            solver,
            check,
            &curve,
            chord,
            own,
            t,
            0.0,
            Judge::StaysInside(&faces),
            &mut near,
        )?;
        if !judged.found || judged.keeps {
            continue;
        }
        // the largest radius whose section here stays inside the part
        let mut at = near.clone();
        let (mut lo, mut hi, mut lo_setback) = (0.0_f64, own, 0.0);
        let mut tr = 0.5 * own;
        while tr > 1.0e-6 * own {
            let s = keeps(
                k,
                target,
                solver,
                check,
                &curve,
                chord,
                own,
                t,
                tr,
                Judge::StaysInside(&faces),
                &mut at,
            )?;
            if s.found && s.keeps {
                lo = tr;
                lo_setback = s.setback;
                break;
            }
            hi = tr;
            tr *= 0.5;
        }
        let mut i = 0;
        while i < BISECTIONS && lo > 0.0 {
            let mid = 0.5 * (lo + hi);
            let s = keeps(
                k,
                target,
                solver,
                check,
                &curve,
                chord,
                own,
                t,
                mid,
                Judge::StaysInside(&faces),
                &mut at,
            )?;
            if s.found && s.keeps {
                lo = mid;
                lo_setback = s.setback;
            } else {
                hi = mid;
            }
            i += 1;
        }
        if out.is_none_or(|c| lo < c.radius) {
            out = Some(BandCut {
                radius: lo,
                setback: lo_setback,
            });
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]
    use super::*;
    use crate::kernel::CancelToken;
    use crate::primitive::api_primitives::create_box;

    /// A 10 mm cube and the section of a 2 mm round of its top front edge at its middle: the
    /// station's plane `x = 5`, the ball's centre 2 in from the top and the front.
    fn cube_section(k: &mut Kernel) -> (Brep, Point3, [Point3; 2]) {
        let b = create_box(k, Point3::ZERO, 10.0, 10.0, 10.0).unwrap();
        let centre = Point3::new(5.0, 2.0, 8.0);
        (
            b,
            centre,
            [Point3::new(5.0, 2.0, 10.0), Point3::new(5.0, 0.0, 8.0)],
        )
    }

    /// The section keeps the cube behind it, judged against the faces and edges within its
    /// reach only (#2648): the top and the front (the section meets them at its rails), no
    /// edge (the four crossing the plane lie farther than the ball's radius from its centre),
    /// and each face's cache, the body's cache and the body's box made once over two
    /// judgements.
    #[test]
    fn a_section_is_judged_against_the_geometry_within_its_reach() {
        let mut k = Kernel::default();
        let (mut b, centre, rails) = cube_section(&mut k);
        let mut check = BandCheck::new();
        let mut station = Station::new(
            &k,
            &mut check,
            &b,
            Point3::new(5.0, 0.0, 10.0),
            Vec3::new(1.0, 0.0, 0.0),
        )
        .unwrap();
        assert_eq!(station.near.len(), 4, "the edges along x reach the plane");
        k.work_counters.reset();
        for _ in 0..2 {
            assert!(keeps_material_behind(
                &k,
                &mut b,
                &mut check,
                &mut station,
                centre,
                rails,
                false
            )
            .unwrap());
        }
        assert_eq!(k.work_counters.band_edge_solves(), 0);
        assert_eq!(k.work_counters.band_face_solves(), 4);
        assert_eq!(k.work_counters.band_face_caches(), 2);
        assert_eq!(k.work_counters.band_brep_caches(), 1);
        assert_eq!(k.work_counters.band_brep_boxes(), 1);
    }

    /// A tolerance of the body grown between two judgements (as a UV trim curve made for a
    /// solve grows its edge's and vertices') makes the body's cache again for the second,
    /// whose boxes carry the tolerances, and the body's box again with it: the answer is
    /// that of a cache made for it.
    #[test]
    fn a_tolerance_grown_between_judgements_makes_the_body_s_cache_again() {
        let mut k = Kernel::default();
        let (mut b, centre, rails) = cube_section(&mut k);
        let mut check = BandCheck::new();
        let mut station = Station::new(
            &k,
            &mut check,
            &b,
            Point3::new(5.0, 0.0, 10.0),
            Vec3::new(1.0, 0.0, 0.0),
        )
        .unwrap();
        k.work_counters.reset();
        let e = b.edges().ids().next().unwrap();
        for grown in [None, Some(2.0 * b.edge(e).unwrap().tolerance), None] {
            if let Some(tolerance) = grown {
                edge::set_tolerance(&mut b, e, tolerance, false, true).unwrap();
            }
            assert!(keeps_material_behind(
                &k,
                &mut b,
                &mut check,
                &mut station,
                centre,
                rails,
                false
            )
            .unwrap());
        }
        assert_eq!(k.work_counters.band_brep_caches(), 2);
        assert_eq!(k.work_counters.band_brep_boxes(), 2);
    }

    /// A token fired before a judgement stops it at the band check's own poll, with the
    /// check's operation (#2648: the phase polled no token before).
    #[test]
    fn a_fired_token_stops_a_judgement() {
        let mut k = Kernel::default();
        let (mut b, centre, rails) = cube_section(&mut k);
        let token = CancelToken::new();
        k.set_cancel_token(token.clone());
        token.cancel();
        let mut check = BandCheck::new();
        let mut station = Station::new(
            &k,
            &mut check,
            &b,
            Point3::new(5.0, 0.0, 10.0),
            Vec3::new(1.0, 0.0, 0.0),
        )
        .unwrap();
        for chord in [false, true] {
            let judged =
                keeps_material_behind(&k, &mut b, &mut check, &mut station, centre, rails, chord);
            assert_eq!(
                judged,
                Err(SolidError::Cancelled { operation: OP }),
                "chord {chord}"
            );
        }
        let faces = ThroughFaces {
            skip: BTreeSet::new(),
            reached: BTreeSet::new(),
        };
        let judged = section_breaks_through(&k, &mut b, &mut check, centre, rails, false, &faces);
        assert_eq!(judged, Err(SolidError::Cancelled { operation: OP }));
    }
}
