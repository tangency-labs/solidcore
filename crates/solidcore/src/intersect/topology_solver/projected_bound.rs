//! Whole-domain bounds for a projected point/surface maximum.
//!
//! Not upstream. Candidate generation is deliberately independent of this check:
//! every knot rectangle must exclude a value larger than the candidate plus tolerance.
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
        if !(value - Ivl::from(reported)).mag().is_finite()
            || (value - Ivl::from(reported)).mag() > solver.tolerance_3d
        {
            return Err(refused(
                "projected candidate value is not enclosed within tolerance",
            ));
        }
    }
    let upper = (Ivl::from(reported) + Ivl::from(solver.tolerance_3d)).lo;
    if !maximum_is_bounded(k, &net, &domain, point, normal, upper)? {
        return Err(refused(
            "whole-face projected maximum could not be certified within tolerance",
        ));
    }
    Ok(())
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
    if !crate::topology::trimming_tools::face_uv_curves_present(b, f)? {
        return Ok(false);
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
    Ok(boundary_is_rectangle(&boundary, domain))
}

fn boundary_is_rectangle(
    boundary: &crate::topology::face::p1::OrientedTrimBoundary,
    domain: &Extent2d,
) -> bool {
    use crate::curve::curve::CurveKind;
    use crate::nurbs::math::NOW;
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
        let CurveKind::BSpline(bs) = curve.kind() else {
            return false;
        };
        let Some(net) = bs.gw_nurb() else {
            return false;
        };
        let [a, b] = net.cpts() else { return false };
        let [t0, t1, t2, t3] = net.knots() else {
            return false;
        };
        if net.degree() != 1
            || t0 != t1
            || t2 != t3
            || t1 >= t2
            || curve.natural_interval().min() != *t1
            || curve.natural_interval().max() != *t2
            || ![a.w, b.w].iter().all(|w| *w == NOW || *w == 1.0)
        {
            return false;
        }
        let (a, b) = (Vec2::new(a.x, a.y), Vec2::new(b.x, b.y));
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
pub(super) fn maximum_is_bounded(
    k: &Kernel,
    net: &Nurbs,
    domain: &Extent2d,
    point: Vec3,
    normal: Vec3,
    upper: f64,
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
                match bound(net, cell, point, normal) {
                    Ok(value) if value.lo <= value.hi && value.sqrt().hi <= upper => continue,
                    Ok(value) if value.sqrt().lo > upper => return Ok(false),
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
