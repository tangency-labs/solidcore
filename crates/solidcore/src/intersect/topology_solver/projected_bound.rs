//! Whole-domain bounds for a projected point/surface maximum.
//!
//! Not upstream. Candidate generation is deliberately independent of this check:
//! every knot rectangle must exclude a value larger than the candidate plus retention
//! slack and the arithmetic enclosure widths of individual evaluated points.
//! The mean-value enclosure uses the existing outward-rounded NURBS partials.
//! Deviations: no recursion; at most 16,384 boxes and 48 bisections along a path.

use crate::blend::fillet_section::certificate::{knot_cuts, nurbs_partials_over_piece};
use crate::blend::fillet_section::interval::Ivl;
use crate::blend::fillet_section::series::V3;
use crate::blend::fillet_section::station::Nurbs;
use crate::error::{Result, SolidError};
use crate::geom::extent2d::Extent2d;
use crate::kernel::Kernel;
use crate::math::Vec2;
use crate::math::Vec3;

const OP: &str = "intersect::topology_solver::projected_bound";
const MAX_BOXES: usize = 16_384;
const MAX_DEPTH: usize = 48;

fn refused(reason: &'static str) -> SolidError {
    SolidError::Unsupported {
        operation: OP,
        reason,
    }
}

/// Certify the scalar maximum for the single-answer point/face operation.
/// A previous answer on another object supplies a feasible lower bound; a new
/// answer on this face must additionally pass the exact trim and point checks.
pub(super) fn certify<O: crate::geom::solution_array::SolutionObject>(
    k: &Kernel,
    objects: &dyn super::p1::TopologyObjects<O>,
    face: O,
    point: Vec3,
    solver: &crate::geom::global_solver::GlobalSolver<'_, O>,
) -> Result<()> {
    let (net, domain, normal, candidate_lower, rounding) =
        checked_candidate(k, objects, face, point, solver)?;
    let upper =
        (Ivl::from(candidate_lower) + Ivl::from(solver.tolerance_3d) + Ivl::from(rounding)).lo;
    if !maximum_with_rounding_is_bounded(k, &net, &domain, point, normal, upper, true)? {
        return Err(refused(
            "whole-face projected maximum could not be certified within tolerance",
        ));
    }
    Ok(())
}

/// Check a candidate before it can supply a feasible bound for tree pruning.
/// The whole-face upper bound may be checked against a later, better candidate.
pub(super) fn certify_candidate<O: crate::geom::solution_array::SolutionObject>(
    k: &Kernel,
    objects: &dyn super::p1::TopologyObjects<O>,
    face: O,
    point: Vec3,
    solver: &crate::geom::global_solver::GlobalSolver<'_, O>,
) -> Result<()> {
    checked_candidate(k, objects, face, point, solver).map(|_| ())
}

fn checked_candidate<O: crate::geom::solution_array::SolutionObject>(
    k: &Kernel,
    objects: &dyn super::p1::TopologyObjects<O>,
    face: O,
    point: Vec3,
    solver: &crate::geom::global_solver::GlobalSolver<'_, O>,
) -> Result<(Nurbs, Extent2d, Vec3, f64, f64)> {
    use crate::blend::fillet_section::station::as_b_spline;
    use crate::surface::surface::SurfaceKind;
    let surface = objects.surface(objects.face_surface(face)?)?;
    if surface.procedural_blend().is_some() || surface.procedural_sweep().is_some() {
        return Err(refused(
            "projected maximum needs an enclosure of the surface definition",
        ));
    }
    match surface.kind_with_nurb(k)? {
        SurfaceKind::BSpline(_) => {}
        SurfaceKind::Plane(_)
        | SurfaceKind::Cone(_)
        | SurfaceKind::Sphere(_)
        | SurfaceKind::Torus(_) => {
            crate::surface::implicit::AnalyticSupport::check(k, surface)?;
        }
        _ => return Err(refused("projected maximum has no exact surface enclosure")),
    }
    let bs = as_b_spline(k, surface)?.ok_or_else(|| refused("missing surface net"))?;
    let net = Nurbs::of(&bs)?;
    let domain = objects.face_uv_domain(face)?;
    let normal = *solver
        .opt_vectors
        .first()
        .ok_or_else(|| refused("missing projection vector"))?;
    let solution = solver
        .solutions
        .as_slice()
        .first()
        .ok_or_else(|| refused("projected maximum has no feasible candidate"))?;
    let reported = solution.start.solution_value;
    if !reported.is_finite()
        || !solver.tolerance_3d.is_finite()
        || solver.tolerance_3d < 0.0
        || ![point.x, point.y, point.z, normal.x, normal.y, normal.z]
            .iter()
            .all(|x| x.is_finite())
    {
        return Err(refused("projected maximum has no finite error budget"));
    }
    let mut rounding = 0.0;
    // Other objects supply a previously validated feasible incumbent. Without
    // this face's point enclosure, grant it no additional arithmetic allowance.
    let mut candidate_lower = reported;
    if solution.objects.contains(&Some(face)) {
        if solution.solution_type != crate::geom::solution_array::SolutionType::SingleValue {
            return Err(refused(
                "projected maximum requires a single feasible point",
            ));
        }
        let (b, f) = objects
            .face_entity(face)
            .ok_or_else(|| refused("missing face boundary"))?;
        if !is_exact_rectangle(k, b, f, &domain)? {
            return Err(refused(
                "candidate feasibility needs an exact rectangular trim boundary",
            ));
        }
        let uv = Vec2::new(solution.start.param(0)?, solution.start.param(1)?);
        if !uv.x.is_finite() || !uv.y.is_finite() || !domain.contains_point(uv, 0.0) {
            return Err(refused("projected candidate is outside the face domain"));
        }
        let value = value_at(&net, uv, point, normal)?;
        candidate_lower = value.lo;
        rounding = (Ivl::from(value.hi) - Ivl::from(value.lo)).hi;
        let budget = (Ivl::from(solver.tolerance_3d) + Ivl::from(rounding)).lo;
        if !(value - Ivl::from(reported)).mag().is_finite()
            || !rounding.is_finite()
            || (value - Ivl::from(reported)).mag() > budget
        {
            return Err(refused(
                "projected candidate value is not enclosed within tolerance",
            ));
        }
    }
    Ok((net, domain, normal, candidate_lower, rounding))
}

#[derive(Clone, Copy)]
struct Cell {
    u: Ivl,
    v: Ivl,
    ku: usize,
    kv: usize,
    depth: usize,
}

fn dot(a: V3<Ivl>, b: V3<Ivl>) -> Ivl {
    (a.x * b.x + a.y * b.y) + a.z * b.z
}

fn cross(a: V3<Ivl>, b: V3<Ivl>) -> V3<Ivl> {
    V3 {
        x: a.y * b.z - a.z * b.y,
        y: a.z * b.x - a.x * b.z,
        z: a.x * b.y - a.y * b.x,
    }
}

fn projected(v: V3<Ivl>, normal: V3<Ivl>) -> V3<Ivl> {
    cross(cross(normal, v), normal)
}

fn relative(v: V3<Ivl>, point: Vec3) -> V3<Ivl> {
    V3 {
        x: v.x - Ivl::from(point.x),
        y: v.y - Ivl::from(point.y),
        z: v.z - Ivl::from(point.z),
    }
}

fn square(x: Ivl) -> Ivl {
    Ivl::checked(x.mig() * x.mig(), x.mag() * x.mag())
}

fn norm_squared(v: V3<Ivl>) -> Ivl {
    (square(v.x) + square(v.y)) + square(v.z)
}

/// Enclose the metric at a feasible parameter, including arithmetic rounding.
pub(super) fn value_at(net: &Nurbs, uv: Vec2, point: Vec3, normal: Vec3) -> Result<Ivl> {
    let p = crate::blend::fillet_section::certificate::nurbs_partials_over_box(
        net, uv.x, uv.x, uv.y, uv.y, 0, 0,
    )?;
    Ok(norm_squared(projected(
        relative(p.at(0, 0), point),
        V3::from_vec3(normal),
    ))
    .sqrt())
}

/// Check the stored boundary, without trusting the tolerance-based rectangular flag.
pub(super) fn is_exact_rectangle(
    k: &Kernel,
    b: &crate::topology::Brep,
    f: crate::topology::FaceId,
    domain: &Extent2d,
) -> Result<bool> {
    Ok(exact_rectangle_boundary(k, b, f, domain)?.is_some())
}

/// Return the exact stored boundary curves only after proving the rectangle.
pub(super) fn exact_rectangle_boundary(
    k: &Kernel,
    b: &crate::topology::Brep,
    f: crate::topology::FaceId,
    domain: &Extent2d,
) -> Result<Option<Vec<crate::curve::curve::Curve>>> {
    if !crate::topology::trimming_tools::face_uv_curves_present(b, f)? {
        return Ok(None);
    }
    let boundary = crate::topology::face::p1::create_oriented_trim_boundary(
        k,
        &mut b.clone_standing_for_self()?,
        f,
        true,
        false,
        None,
        &|entity| b.owner_identity(entity),
    )?;
    Ok(boundary_is_rectangle(&boundary, domain).then_some(boundary.curves))
}

fn boundary_is_rectangle(
    boundary: &crate::topology::face::p1::OrientedTrimBoundary,
    domain: &Extent2d,
) -> bool {
    if boundary.loops != [4] || boundary.curves.len() != 4 {
        return false;
    }
    let corners = [
        domain.min,
        Vec2::new(domain.max.x, domain.min.y),
        domain.max,
        Vec2::new(domain.min.x, domain.max.y),
    ];
    let mut edges = 0_u8;
    let mut first = None;
    let mut previous = None;
    for curve in &boundary.curves {
        let Some((a, b)) = rectangle_segment_ends(curve, &corners) else {
            return false;
        };
        let Some(i) = corners.iter().position(|p| *p == a) else {
            return false;
        };
        let Some(j) = corners.iter().position(|p| *p == b) else {
            return false;
        };
        let edge = if (i + 1) % 4 == j {
            i
        } else if (j + 1) % 4 == i {
            j
        } else {
            return false;
        };
        let bit = 1_u8 << edge;
        if edges & bit != 0 || previous.is_some_and(|p| p != a) {
            return false;
        }
        edges |= bit;
        first.get_or_insert(a);
        previous = Some(b);
    }
    edges == 15 && previous == first
}

/// Recognize affine trim definitions, rather than their sampled endpoints.
fn rectangle_segment_ends(
    curve: &crate::curve::curve::Curve,
    corners: &[Vec2; 4],
) -> Option<(Vec2, Vec2)> {
    use crate::curve::curve::p1::dyadic_four_products_zero;
    use crate::curve::curve::CurveKind;
    use crate::nurbs::math::NOW;
    if let CurveKind::Line(line) = curve.kind() {
        let natural = line.natural_interval();
        let analytic = line.step_interval();
        if !line.is_bounded()
            || !natural.length().is_finite()
            || natural.length() <= 0.0
            || !analytic.length().is_finite()
            || analytic.length() <= 0.0
            || crate::geom::tol::is_zero(natural.length())
        {
            return None;
        }
        let p = line.line_point();
        let v = line.line_vector();
        let scale = line.line_scale();
        let component_matches = |origin, direction, parameter, target| {
            let scaled = scale * direction;
            dyadic_four_products_zero([
                (scale, direction, false),
                (scaled, 1.0, true),
                (0.0, 0.0, false),
                (0.0, 0.0, false),
            ]) && dyadic_four_products_zero([
                (origin, 1.0, false),
                (parameter, scaled, false),
                (target, 1.0, true),
                (0.0, 0.0, false),
            ])
        };
        let endpoint = |t| {
            corners.iter().copied().find(|c| {
                component_matches(p.x, v.x, t, c.x) && component_matches(p.y, v.y, t, c.y)
            })
        };
        let (a, b) = if line.inside_out() {
            (analytic.max(), analytic.min())
        } else {
            (analytic.min(), analytic.max())
        };
        return Some((endpoint(a)?, endpoint(b)?));
    }
    let CurveKind::BSpline(bs) = curve.kind() else {
        return None;
    };
    let net = bs.gw_nurb()?;
    let [a, b] = net.cpts() else { return None };
    let [t0, t1, t2, t3] = net.knots() else {
        return None;
    };
    if net.degree() != 1
        || t0 != t1
        || t2 != t3
        || t1 >= t2
        || curve.natural_interval().min() != *t1
        || curve.natural_interval().max() != *t2
        || ![a.w, b.w].iter().all(|w| *w == NOW || *w == 1.0)
    {
        return None;
    }
    Some((Vec2::new(a.x, a.y), Vec2::new(b.x, b.y)))
}

/// A conservative distance range between two projected boxes, in the same
/// double-cross metric as the queries. In particular, N need not be unit length.
pub(super) fn projected_box_distance(
    a: &crate::geom::extent3d::Extent3d,
    b: &crate::geom::extent3d::Extent3d,
    normal: Vec3,
) -> Ivl {
    if [a, b].iter().any(|extent| {
        let (lo, hi) = (extent.min(), extent.max());
        ![lo.x, lo.y, lo.z, hi.x, hi.y, hi.z]
            .iter()
            .all(|x| x.is_finite())
            || lo.x > hi.x
            || lo.y > hi.y
            || lo.z > hi.z
    }) {
        return Ivl::whole();
    }
    let delta = V3 {
        x: Ivl::new(b.min().x, b.max().x) - Ivl::new(a.min().x, a.max().x),
        y: Ivl::new(b.min().y, b.max().y) - Ivl::new(a.min().y, a.max().y),
        z: Ivl::new(b.min().z, b.max().z) - Ivl::new(a.min().z, a.max().z),
    };
    norm_squared(projected(delta, V3::from_vec3(normal))).sqrt()
}

fn bound(net: &Nurbs, cell: Cell, point: Vec3, normal: V3<Ivl>) -> Result<Ivl> {
    let p = nurbs_partials_over_piece(
        net, cell.ku, cell.kv, cell.u.lo, cell.u.hi, cell.v.lo, cell.v.hi, 1,
    )?;
    let q = projected(relative(p.at(0, 0), point), normal);
    let direct = norm_squared(q);
    let (u, v) = (cell.u.mid(), cell.v.mid());
    let middle = nurbs_partials_over_piece(net, cell.ku, cell.kv, u, u, v, v, 0)?;
    let at_middle = norm_squared(projected(relative(middle.at(0, 0), point), normal));
    let du = Ivl::from(2.0) * dot(q, projected(p.at(1, 0), normal));
    let dv = Ivl::from(2.0) * dot(q, projected(p.at(0, 1), normal));
    let centred = at_middle + du * (cell.u - Ivl::from(u)) + dv * (cell.v - Ivl::from(v));
    Ok(Ivl::new(
        direct.lo.max(centred.lo),
        direct.hi.min(centred.hi),
    ))
}

/// Prove an upper bound on the projected distance over the supplied rectangle.
/// This alone does NOT certify candidate feasibility, trimming, or an analytic
/// definition's agreement with this net; callers must establish those separately.
#[cfg(test)]
pub(super) fn maximum_is_bounded(
    k: &Kernel,
    net: &Nurbs,
    domain: &Extent2d,
    point: Vec3,
    normal: Vec3,
    upper: f64,
) -> Result<bool> {
    maximum_with_rounding_is_bounded(k, net, domain, point, normal, upper, false)
}

/// Point evaluations contribute arithmetic enclosure width, never the width of
/// a parameter box. This keeps zero retention slack distinct from exact arithmetic.
fn maximum_with_rounding_is_bounded(
    k: &Kernel,
    net: &Nurbs,
    domain: &Extent2d,
    point: Vec3,
    normal: Vec3,
    upper: f64,
    allow_point_rounding: bool,
) -> Result<bool> {
    if !upper.is_finite()
        || upper < 0.0
        || ![domain.min.x, domain.min.y, domain.max.x, domain.max.y]
            .iter()
            .all(|x| x.is_finite())
        || domain.min.x >= domain.max.x
        || domain.min.y >= domain.max.y
    {
        return Ok(false);
    }
    let us = knot_cuts(&net.u, net.p, net.nu, domain.min.x, domain.max.x, 0)?;
    let vs = knot_cuts(&net.v, net.q, net.nv, domain.min.y, domain.max.y, 0)?;
    let normal_vector = normal;
    let normal = V3::from_vec3(normal);
    let mut visited = 0;
    for &(ua, ub, ku) in &us {
        for &(va, vb, kv) in &vs {
            let mut pending = vec![Cell {
                u: Ivl::new(ua, ub),
                v: Ivl::new(va, vb),
                ku,
                kv,
                depth: 0,
            }];
            while let Some(cell) = pending.pop() {
                k.check_cancelled(OP)?;
                visited += 1;
                if visited > MAX_BOXES {
                    return Ok(false);
                }
                let cell_upper = if allow_point_rounding {
                    let value = value_at(
                        net,
                        Vec2::new(cell.u.mid(), cell.v.mid()),
                        point,
                        normal_vector,
                    )?;
                    let width = (Ivl::from(value.hi) - Ivl::from(value.lo)).hi;
                    if !width.is_finite() {
                        return Ok(false);
                    }
                    (Ivl::from(upper) + Ivl::from(width)).lo
                } else {
                    upper
                };
                match bound(net, cell, point, normal) {
                    Ok(value) if value.lo <= value.hi && value.sqrt().hi <= cell_upper => continue,
                    Ok(value) if value.sqrt().lo > cell_upper => return Ok(false),
                    Err(e @ SolidError::Cancelled { .. }) => return Err(e),
                    // A denominator not separated from zero may become provable
                    // after subdivision. An unresolved box never counts as passed.
                    _ => {}
                }
                if cell.depth == MAX_DEPTH {
                    return Ok(false);
                }
                let split_u =
                    (cell.u.hi - cell.u.lo) / (ub - ua) >= (cell.v.hi - cell.v.lo) / (vb - va);
                let range = if split_u { cell.u } else { cell.v };
                let middle = range.mid();
                if !(middle > range.lo && middle < range.hi) {
                    return Ok(false);
                }
                let mut left = Cell {
                    depth: cell.depth + 1,
                    ..cell
                };
                let mut right = left;
                if split_u {
                    left.u.hi = middle;
                    right.u.lo = middle;
                } else {
                    left.v.hi = middle;
                    right.v.lo = middle;
                }
                pending.push(right);
                pending.push(left);
            }
        }
    }
    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hidden_peak() -> Nurbs {
        // S(u,v)=(4u(1-u),4v(1-v),0): all four corners project to
        // the origin, but the interior reaches (1,1,0).
        let mut pw = Vec::new();
        for x in [0.0, 2.0, 0.0] {
            for y in [0.0, 2.0, 0.0] {
                pw.extend([x, y, 0.0, 1.0]);
            }
        }
        Nurbs {
            p: 2,
            q: 2,
            nu: 3,
            nv: 3,
            u: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            v: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            pw,
        }
    }

    #[test]
    fn zero_retention_uses_only_point_evaluation_rounding() {
        let k = Kernel::default();
        let net = hidden_peak();
        let domain = Extent2d::from_bounds(0.0, 0.0, 1.0, 1.0);
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let value = value_at(&net, Vec2::new(0.5, 0.5), Vec3::ZERO, normal).unwrap();
        let width = (Ivl::from(value.hi) - Ivl::from(value.lo)).hi;
        let upper = (Ivl::from(value.lo) + Ivl::from(width)).lo;
        assert!(maximum_with_rounding_is_bounded(
            &k,
            &net,
            &domain,
            Vec3::ZERO,
            normal,
            upper,
            true,
        )
        .unwrap());
        assert!(!maximum_with_rounding_is_bounded(
            &k,
            &net,
            &domain,
            Vec3::ZERO,
            normal,
            1.0,
            true,
        )
        .unwrap());
    }

    #[test]
    fn analytic_trim_proof_uses_exact_members_not_a_stale_net_or_rounded_endpoints() {
        use crate::curve::curve::Curve;
        use crate::curve::line::Line;
        use crate::geom::extent1d::Extent1d;
        use crate::nurbs::b_spline_curve::BSplineCurve;
        let y = 1e16;
        let corners = [
            Vec2::new(0.0, y),
            Vec2::new(1.0, y),
            Vec2::new(1.0, y + 2.0),
            Vec2::new(0.0, y + 2.0),
        ];
        let base =
            BSplineCurve::create_line_segment(2, Vec3::new(0.0, y, 0.0), Vec3::new(1.0, y, 0.0))
                .unwrap();
        let line = |dy| {
            Curve::from(Line::from_parts(
                base.clone(),
                Vec3::new(0.0, y, 0.0),
                Vec3::new(1.0, dy, 0.0),
                Extent1d::from_min_max(0.0, 1.0),
                1.0,
                false,
            ))
        };
        assert_eq!(
            rectangle_segment_ends(&line(0.0), &corners),
            Some((corners[0], corners[1]))
        );
        let tilted = line(0.25);
        assert_eq!(tilted.evaluate_point(&Kernel::default(), 1.0).unwrap().y, y);
        assert!(rectangle_segment_ends(&tilted, &corners).is_none());
    }

    #[test]
    fn projected_box_bounds_preserve_the_double_cross_scale() {
        use crate::geom::extent3d::Extent3d;
        let a = Extent3d::from_point(Vec3::new(0.2, 0.3, 2.0));
        let b = Extent3d::from_points(Vec3::new(-0.5, -0.5, -1.5), Vec3::new(0.5, 0.5, -0.5));
        for scale in [0.001, 1.0, 1000.0] {
            let n = Vec3::new(0.5, 0.0, 0.5) * scale;
            let range = projected_box_distance(&a, &b, n);
            for x in [-0.5, 0.5] {
                for y in [-0.5, 0.5] {
                    for z in [-1.5, -0.5] {
                        let value = n.cross(Vec3::new(x, y, z) - a.min()).cross(n).length();
                        assert!(range.lo <= value && value <= range.hi);
                    }
                }
            }
        }
    }

    #[test]
    fn almost_rectangular_trims_do_not_prove_feasibility() {
        use crate::curve::curve::Curve;
        use crate::nurbs::b_spline_curve::BSplineCurve;
        use crate::topology::face::p1::OrientedTrimBoundary;
        let domain = Extent2d::from_bounds(0.0, 0.0, 1.0, 1.0);
        let boundary = |inset| {
            let corners = [
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(1.0 - inset, 1.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
            ];
            OrientedTrimBoundary {
                loops: vec![4],
                curves: corners
                    .iter()
                    .zip(corners.iter().cycle().skip(1))
                    .map(|(&a, &b)| {
                        Curve::from(BSplineCurve::create_line_segment(2, a, b).unwrap())
                    })
                    .collect(),
            }
        };
        assert!(boundary_is_rectangle(&boundary(0.0), &domain));
        assert!(!boundary_is_rectangle(&boundary(1e-12), &domain));
        let mut disconnected = boundary(0.0);
        disconnected.curves.swap(1, 2);
        assert!(!boundary_is_rectangle(&disconnected, &domain));
        let mut duplicate = boundary(0.0);
        duplicate.curves[3] = duplicate.curves[2].clone();
        assert!(!boundary_is_rectangle(&duplicate, &domain));
    }

    #[test]
    fn a_missed_interior_maximum_cannot_be_certified() {
        let k = Kernel::default();
        let net = hidden_peak();
        let domain = Extent2d::from_bounds(0.0, 0.0, 1.0, 1.0);
        let normal = Vec3::new(0.0, 0.0, 1.0);
        assert!(!maximum_is_bounded(&k, &net, &domain, Vec3::ZERO, normal, 1.0).unwrap());
        assert!(
            maximum_is_bounded(&k, &net, &domain, Vec3::ZERO, normal, 2.0_f64.sqrt() + 1e-6,)
                .unwrap()
        );
    }

    #[test]
    fn non_unit_projection_is_bounded_in_the_requested_metric() {
        let k = Kernel::default();
        let net = hidden_peak();
        let domain = Extent2d::from_bounds(0.0, 0.0, 1.0, 1.0);
        let normal = Vec3::new(0.0, 0.0, 2.0);
        assert!(!maximum_is_bounded(&k, &net, &domain, Vec3::ZERO, normal, 2.0).unwrap());
        assert!(maximum_is_bounded(
            &k,
            &net,
            &domain,
            Vec3::ZERO,
            normal,
            4.0 * 2.0_f64.sqrt() + 1e-6,
        )
        .unwrap());
    }

    #[test]
    fn unresolved_rounding_and_invalid_domains_refuse() {
        let k = Kernel::default();
        let net = hidden_peak();
        let domain = Extent2d::from_bounds(0.0, 0.0, 1.0, 1.0);
        assert!(!maximum_is_bounded(&k, &net, &domain, Vec3::ZERO, Vec3::ZERO, 0.0).unwrap());
        let invalid = Extent2d::from_bounds(0.0, 0.0, f64::INFINITY, 1.0);
        assert!(!maximum_is_bounded(
            &k,
            &net,
            &invalid,
            Vec3::ZERO,
            Vec3::new(0.0, 0.0, 1.0),
            2.0,
        )
        .unwrap());
    }

    #[test]
    fn cancellation_is_not_reported_as_an_unproved_bound() {
        let token = crate::kernel::CancelToken::default();
        token.cancel();
        let k = Kernel {
            cancel: Some(token),
            ..Kernel::default()
        };
        let result = maximum_is_bounded(
            &k,
            &hidden_peak(),
            &Extent2d::from_bounds(0.0, 0.0, 1.0, 1.0),
            Vec3::ZERO,
            Vec3::new(0.0, 0.0, 1.0),
            2.0,
        );
        assert!(matches!(result, Err(SolidError::Cancelled { .. })));
    }
}
