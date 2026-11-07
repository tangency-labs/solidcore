//! Native projected maxima for an affine curve against a rectangular face.
//!
//! Not upstream. Norm convexity reduces every continuous affine span to its
//! endpoints. Native UV refinement proposes candidates; the shared homogeneous
//! Bernstein residual certifies the whole face for every exact endpoint.
//! Deviations: at most 128 defining curve spans and 8,192 UV seeds; native
//! refinement and residual subdivision retain their existing finite budgets.

use super::p1::TopologyObjects;
use super::{projected_bound, projected_residual};
use crate::blend::fillet_section::certificate::knot_cuts;
use crate::blend::fillet_section::interval::Ivl;
use crate::blend::fillet_section::station::Nurbs;
use crate::curve::curve::Curve;
use crate::error::{Result, SolidError};
use crate::geom::extent2d::Extent2d;
use crate::geom::solution_array::SolutionObject;
use crate::math::{Vec2, Vec3};
use crate::Kernel;

const OP: &str = "intersect::topology_solver::projected_affine_face";
pub(super) fn refused(reason: &'static str) -> SolidError {
    SolidError::Unsupported {
        operation: OP,
        reason,
    }
}

pub(super) struct Candidate {
    pub parameters: [f64; 3],
    pub reported: f64,
    pub enclosure: Ivl,
}

/// Construct on the native definitions without treating a fitted silhouette as
/// either operand. The returned candidate is not accepted until `certify` runs.
pub(super) fn construct<O: SolutionObject>(
    k: &Kernel,
    objects: &dyn TopologyObjects<O>,
    face: O,
    curve: &Curve,
    endpoints: &[(f64, [Ivl; 3])],
    normal: Vec3,
) -> Result<(Candidate, Nurbs, Extent2d)> {
    if ![normal.x, normal.y, normal.z].iter().all(|x| x.is_finite()) {
        return Err(refused("nonfinite projection vector"));
    }
    let net = projected_bound::face_net(k, objects, face)?;
    let domain = objects.face_uv_domain(face)?;
    let (brep, face_id) = objects
        .face_entity(face)
        .ok_or_else(|| refused("missing native face boundary"))?;
    if !projected_bound::is_exact_rectangle(k, brep, face_id, &domain)? {
        return Err(refused(
            "affine pair maximum requires an exact rectangular face domain",
        ));
    }
    let surface = objects.surface(objects.face_surface(face)?)?;
    let us = knot_cuts(&net.u, net.p, net.nu, domain.min().x, domain.max().x, 0)?;
    let vs = knot_cuts(&net.v, net.q, net.nv, domain.min().y, domain.max().y, 0)?;
    let mut best: Option<Candidate> = None;
    let mut seeds = 0usize;
    for &(t, definition) in endpoints {
        k.check_cancelled(OP)?;
        let point = curve.evaluate_point(k, t)?;
        for &(a, b, _) in &us {
            for &(c, d, _) in &vs {
                for u in [a, a + (b - a) * 0.5, b] {
                    for v in [c, c + (d - c) * 0.5, d] {
                        seeds += 1;
                        if seeds > 8192 {
                            return Err(refused(
                                "native affine pair candidate seed budget exhausted",
                            ));
                        }
                        let uv = projected_bound::refine_uv(
                            k,
                            surface,
                            &domain,
                            point,
                            normal,
                            Vec2::new(u, v),
                        )?;
                        let reported = normal
                            .cross(curve.evaluate_point(k, t)? - surface.evaluate_point(k, uv)?)
                            .cross(normal)
                            .length();
                        let enclosure =
                            projected_bound::value_at_definition(&net, uv, definition, normal)?;
                        if !(reported.is_finite()
                            && enclosure.lo.is_finite()
                            && enclosure.hi.is_finite()
                            && enclosure.lo <= reported
                            && reported <= enclosure.hi
                            && domain.contains_point(uv, 0.0))
                        {
                            return Err(refused(
                                "native pair report is outside its own definition enclosure",
                            ));
                        }
                        if best.as_ref().is_none_or(|old| reported > old.reported) {
                            best = Some(Candidate {
                                parameters: [t, uv.x, uv.y],
                                reported,
                                enclosure,
                            });
                        }
                    }
                }
            }
        }
    }
    Ok((
        best.ok_or_else(|| refused("affine pair has no native candidate"))?,
        net,
        domain,
    ))
}

/// Prove every admissible curve/face pair fits the retained candidate's own bound.
pub(super) fn certify(
    k: &Kernel,
    net: &Nurbs,
    domain: &Extent2d,
    endpoints: &[(f64, [Ivl; 3])],
    normal: Vec3,
    upper: f64,
) -> Result<()> {
    for &(_, point) in endpoints {
        if !projected_residual::maximum_is_bounded_from_definition(
            k, net, domain, point, normal, upper,
        )? {
            return Err(refused(
                "native affine curve/face maximum exceeds the retained pair enclosure",
            ));
        }
    }
    Ok(())
}

/// Reconstruct the actual retained pair; no discarded candidate contributes its
/// arithmetic uncertainty. Unsupported owners explicitly refuse.
pub(super) fn retained_bound<O: SolutionObject>(
    k: &Kernel,
    objects: &dyn TopologyObjects<O>,
    solution: &crate::geom::solution_array::Solution<O>,
    query_owner: O,
    query_curve: &Curve,
    normal: Vec3,
) -> Result<Ivl> {
    use crate::geom::solution_array::{ObjectKind, SolutionType};
    if solution.solution_type != SolutionType::SingleValue || solution.num_objects != 2 {
        return Err(refused("retained affine pair must name two native points"));
    }
    let mut parameter = 0usize;
    let mut positions = Vec::new();
    for handle in solution.objects.iter().take(2) {
        let handle = handle.ok_or_else(|| refused("retained pair has an unnamed operand"))?;
        let position = if handle == query_owner
            || matches!(handle.kind(), ObjectKind::Curve | ObjectKind::Edge)
        {
            let curve = if handle == query_owner {
                query_curve
            } else if handle.kind() == ObjectKind::Edge {
                objects.curve(
                    objects
                        .edge_curve(handle)?
                        .ok_or_else(|| refused("retained edge has no defining curve"))?,
                )?
            } else {
                objects.curve(handle)?
            };
            let t = solution.start.param(parameter)?;
            parameter += 1;
            if handle.kind() == ObjectKind::Edge {
                let domain = objects.edge_interval(handle)?;
                if t < domain.min() || t > domain.max() {
                    return Err(refused("retained edge point is outside its domain"));
                }
            }
            crate::curve::curve::p1::affine_point_definition_enclosure(curve, t)
                .ok_or_else(|| refused("retained pair curve lacks an exact point enclosure"))?
        } else if handle.kind() == ObjectKind::Face {
            let uv = Vec2::new(
                solution.start.param(parameter)?,
                solution.start.param(parameter + 1)?,
            );
            parameter += 2;
            let domain = objects.face_uv_domain(handle)?;
            let (b, f) = objects
                .face_entity(handle)
                .ok_or_else(|| refused("retained face lacks a boundary"))?;
            if !domain.contains_point(uv, 0.0)
                || !projected_bound::is_exact_rectangle(k, b, f, &domain)?
            {
                return Err(refused("retained face point lacks exact domain membership"));
            }
            projected_bound::native_point_at(&projected_bound::face_net(k, objects, handle)?, uv)?
        } else if handle.kind() == ObjectKind::Vertex {
            let point = objects.vertex_point(handle)?;
            [point.x, point.y, point.z].map(Ivl::from)
        } else {
            return Err(refused(
                "retained pair operand has no native point certificate",
            ));
        };
        positions.push(position);
    }
    let [a, b] = positions.as_slice() else {
        return Err(refused("retained pair has incomplete operands"));
    };
    if parameter != solution.num_variables {
        return Err(refused("retained pair parameter layout is ambiguous"));
    }
    let bound = projected_bound::pair_metric(*a, *b, normal);
    if !(bound.lo.is_finite()
        && bound.hi.is_finite()
        && bound.lo <= solution.start.solution_value
        && solution.start.solution_value <= bound.hi)
    {
        return Err(refused(
            "retained pair report is outside its own native enclosure",
        ));
    }
    Ok(bound)
}

/// Independently recover prior-domain coverage; attainment of an earlier fitted
/// answer is not a certificate of the domain it came from.
pub(super) fn certify_previous_domain<O: SolutionObject>(
    k: &Kernel,
    objects: &dyn TopologyObjects<O>,
    solution: &crate::geom::solution_array::Solution<O>,
    query_owner: O,
    query_curve: &Curve,
    normal: Vec3,
    upper: f64,
) -> Result<()> {
    use crate::geom::solution_array::ObjectKind;
    enum Domain {
        Points(Vec<[Ivl; 3]>),
        Face(Nurbs, Extent2d),
    }
    let mut domains = Vec::new();
    for handle in solution.objects.iter().take(2) {
        let handle = handle.ok_or_else(|| refused("previous domain has an unnamed owner"))?;
        let domain = if handle == query_owner
            || matches!(handle.kind(), ObjectKind::Curve | ObjectKind::Edge)
        {
            let curve = if handle == query_owner {
                query_curve
            } else if handle.kind() == ObjectKind::Edge {
                objects.curve(
                    objects
                        .edge_curve(handle)?
                        .ok_or_else(|| refused("previous edge has no curve"))?,
                )?
            } else {
                objects.curve(handle)?
            };
            let interval = if handle.kind() == ObjectKind::Edge {
                objects.edge_interval(handle)?
            } else {
                curve.natural_interval()
            };
            let endpoints = crate::curve::curve::p1::affine_curve_endpoints(curve, interval)?
                .ok_or_else(|| refused("previous curve domain lacks affine coverage provenance"))?;
            Domain::Points(endpoints.into_iter().map(|(_, point)| point).collect())
        } else if handle.kind() == ObjectKind::Vertex {
            let point = objects.vertex_point(handle)?;
            Domain::Points(vec![[point.x, point.y, point.z].map(Ivl::from)])
        } else if handle.kind() == ObjectKind::Face {
            // Bounding the native UV rectangle is conservative even if the old
            // face has additional trims; its returned point was checked separately.
            Domain::Face(
                projected_bound::face_net(k, objects, handle)?,
                objects.face_uv_domain(handle)?,
            )
        } else {
            return Err(refused(
                "previous operand has no domain coverage certificate",
            ));
        };
        domains.push(domain);
    }
    match domains.as_slice() {
        [Domain::Points(a), Domain::Points(b)] => {
            for a in a {
                for b in b {
                    k.check_cancelled(OP)?;
                    let bound = projected_bound::pair_metric(*a, *b, normal);
                    if !(bound.hi.is_finite() && bound.hi <= upper) {
                        return Err(refused(
                            "previous affine domain exceeds final retained enclosure",
                        ));
                    }
                }
            }
        }
        [Domain::Points(points), Domain::Face(net, domain)]
        | [Domain::Face(net, domain), Domain::Points(points)] => {
            for point in points {
                if !projected_residual::maximum_is_bounded_from_definition(
                    k, net, domain, *point, normal, upper,
                )? {
                    return Err(refused(
                        "previous native face domain exceeds final retained enclosure",
                    ));
                }
            }
        }
        _ => {
            return Err(refused(
                "previous pair has no supported whole-domain certificate",
            ))
        }
    }
    Ok(())
}
