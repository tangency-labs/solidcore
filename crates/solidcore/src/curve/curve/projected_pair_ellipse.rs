//! Native-parameter edge certificates for finite ellipse-pair enumeration.
//!
//! Not upstream. A strict edge maximum has a zero first derivative and negative
//! second derivative. Interval exclusion covers the entire edge; strict endpoint
//! signs and a negative derivative enclose exactly one stationary maximum.
//! Deviations: experimental, test-only construction; at most 4096 cells, depth
//! 48 and 64 interval-Newton contractions. Endpoint/multiple-root ambiguity
//! refuses. No production dispatch or All completeness claim yet.

use super::projected_pair::{dot, norm_squared, project, sub, Definition, Piece};
use crate::blend::fillet_section::certificate::curve_derivs_over;
use crate::blend::fillet_section::interval::Ivl;
use crate::blend::fillet_section::series::V3;
use crate::error::{Result, SolidError};
use crate::kernel::Kernel;
use crate::math::Vec3;

const OP: &str = "curve::curve::projected_pair_ellipse";
const MAX_CELLS: usize = 4096;
const MAX_DEPTH: usize = 48;
const NEWTON_STEPS: usize = 64;

fn refused(reason: &'static str) -> SolidError {
    SolidError::Unsupported {
        operation: OP,
        reason,
    }
}

fn finite(x: Ivl) -> bool {
    x.lo.is_finite() && x.hi.is_finite() && x.lo <= x.hi
}

/// One unique strict stationary maximum, in the original parameter coordinates.
/// This parameter enclosure is geometric uncertainty, never point rounding.
#[derive(Clone, Copy, Debug)]
struct EdgeMaximum {
    parameter: Ivl,
    original_bracket: Ivl,
}

/// Half the first and second derivatives of the squared projected edge distance.
fn edge_jet(
    definition: &Definition<'_>,
    piece: Piece,
    fixed: V3<Ivl>,
    normal: Vec3,
    parameter: Ivl,
) -> Result<(Ivl, Ivl)> {
    let Definition::Nurbs(net) = definition else {
        return Err(refused(
            "ellipse edge needs its certified native rational definition",
        ));
    };
    let jet = curve_derivs_over(net, parameter, 2, piece.knot)?;
    let point = *jet
        .first()
        .ok_or_else(|| refused("missing edge point enclosure"))?;
    let first = *jet
        .get(1)
        .ok_or_else(|| refused("missing edge tangent enclosure"))?;
    let second = *jet
        .get(2)
        .ok_or_else(|| refused("missing edge curvature enclosure"))?;
    let n = V3::from_vec3(normal);
    let delta = project(sub(point, fixed), n);
    let first = project(first, n);
    let second = project(second, n);
    Ok((dot(delta, first), dot(first, first) + dot(delta, second)))
}

/// Objective enclosure on a root bracket. Its geometric remainder stays in the
/// bound: only a separately evaluated point interval may supply rounding slack.
fn edge_value_over(
    definition: &Definition<'_>,
    piece: Piece,
    fixed: V3<Ivl>,
    normal: Vec3,
    range: Ivl,
) -> Result<Ivl> {
    let midpoint = range.mid();
    let (point, _) = definition.jet(Ivl::from(midpoint), piece.knot)?;
    let at_midpoint = norm_squared(project(sub(point, fixed), V3::from_vec3(normal)));
    let (first, _) = edge_jet(definition, piece, fixed, normal, Ivl::from(midpoint))?;
    let (_, second) = edge_jet(definition, piece, fixed, normal, range)?;
    let delta = range - Ivl::from(midpoint);
    let squared_delta = Ivl::checked(delta.mig() * delta.mig(), delta.mag() * delta.mag());
    let squared = at_midpoint + Ivl::from(2.0) * first * delta + second * squared_delta;
    let value = squared.sqrt();
    if !finite(value) {
        return Err(refused("root objective has no finite Taylor enclosure"));
    }
    Ok(value)
}

/// Contract a bracket whose strict sign change already proves existence.
/// Every interval Newton intersection retains its unique root; a contraction
/// stopping at arithmetic resolution keeps the remaining geometric enclosure.
fn contract(mut range: Ivl, evaluate: &impl Fn(Ivl) -> Result<(Ivl, Ivl)>) -> Result<Ivl> {
    for _ in 0..NEWTON_STEPS {
        let midpoint = range.mid();
        let (value, _) = evaluate(Ivl::from(midpoint))?;
        let (_, derivative) = evaluate(range)?;
        if !(finite(value) && finite(derivative) && derivative.hi < 0.0) {
            return Err(refused(
                "stationary maximum lost its strict derivative certificate",
            ));
        }
        let newton = Ivl::from(midpoint) - value / derivative;
        let next = Ivl::new(range.lo.max(newton.lo), range.hi.min(newton.hi));
        if !finite(next) {
            return Err(refused("interval Newton lost a previously certified root"));
        }
        if next == range {
            return Ok(range);
        }
        range = next;
    }
    Ok(range)
}

/// Isolate every strict interior edge maximum, or refuse an unresolved region.
/// Endpoint roots belong to the corner construction but are not silently dropped:
/// until an exact endpoint identity can be proved they cause explicit refusal.
fn isolate(
    kernel: &Kernel,
    initial: Ivl,
    evaluate: impl Fn(Ivl) -> Result<(Ivl, Ivl)>,
) -> Result<Vec<EdgeMaximum>> {
    if !(finite(initial) && initial.lo < initial.hi) {
        return Err(refused(
            "edge parameter domain is not finite and increasing",
        ));
    }
    let mut todo = vec![(initial, 0_usize)];
    let mut maxima = Vec::new();
    let mut visited = 0_usize;
    while let Some((range, depth)) = todo.pop() {
        kernel.check_cancelled(OP)?;
        visited += 1;
        if visited > MAX_CELLS {
            return Err(refused("edge stationary-root cell budget exhausted"));
        }
        let (value, derivative) = evaluate(range)?;
        if !(finite(value) && finite(derivative)) {
            return Err(refused("edge stationary function has no finite enclosure"));
        }
        // Nonzero first derivative excludes stationarity. Positive second
        // derivative excludes a local maximum even if a minimum root is present.
        if value.lo > 0.0 || value.hi < 0.0 || derivative.lo > 0.0 {
            continue;
        }
        if derivative.hi < 0.0 {
            let (left, _) = evaluate(Ivl::from(range.lo))?;
            let (right, _) = evaluate(Ivl::from(range.hi))?;
            if !(finite(left) && finite(right)) {
                return Err(refused("edge endpoint derivative has no finite enclosure"));
            }
            if left.hi < 0.0 || right.lo > 0.0 {
                continue;
            }
            if left.lo > 0.0 && right.hi < 0.0 {
                maxima.push(EdgeMaximum {
                    parameter: contract(range, &evaluate)?,
                    original_bracket: range,
                });
                continue;
            }
        }
        let mid = range.mid();
        if depth >= MAX_DEPTH || !(range.lo < mid && mid < range.hi) {
            return Err(refused(
                "edge stationary root or endpoint identity is unresolved",
            ));
        }
        todo.push((Ivl::new(mid, range.hi), depth + 1));
        todo.push((Ivl::new(range.lo, mid), depth + 1));
    }
    Ok(maxima)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curve::curve::Curve;
    use crate::curve::ellipse::{Ellipse, EllipseOptions};
    use crate::geom::extent1d::Extent1d;

    #[test]
    fn edge_isolation_proves_unique_maximum_and_excludes_minimum() {
        let roots = isolate(&Kernel::default(), Ivl::new(-2.0, 2.0), |x| {
            Ok((Ivl::from(0.3) - x, Ivl::from(-1.0)))
        })
        .unwrap();
        assert_eq!(roots.len(), 1);
        assert!(roots[0].parameter.lo <= 0.3 && roots[0].parameter.hi >= 0.3);
        assert!(roots[0].parameter.hi - roots[0].parameter.lo < 1e-14);
        assert_eq!(roots[0].original_bracket, Ivl::new(-2.0, 2.0));
        assert!(isolate(&Kernel::default(), Ivl::new(-2.0, 2.0), |x| {
            Ok((x, Ivl::from(1.0)))
        })
        .unwrap()
        .is_empty());
    }

    #[test]
    fn edge_isolation_refuses_unresolved_endpoint_and_continuum() {
        assert!(isolate(&Kernel::default(), Ivl::new(0.0, 1.0), |x| {
            Ok((Ivl::from(0.0) - x, Ivl::from(-1.0)))
        })
        .is_err());
        assert!(isolate(&Kernel::default(), Ivl::new(0.0, 1.0), |_| {
            Ok((Ivl::from(0.0), Ivl::from(0.0)))
        })
        .is_err());
    }

    #[test]
    fn original_ellipse_edge_encloses_the_stronger_native_parameter_witness() {
        let kernel = Kernel::default();
        let curve = Curve::from(
            Ellipse::new(
                &kernel,
                Vec3::new(1.9892188469923955, -3.246996295465655, -3.4858186076942532),
                Vec3::X,
                Vec3::Y,
                Extent1d::from_min_max(-90.0, 101.59684508222117),
                4.723327346620286,
                2.505285780586296,
                &EllipseOptions::default(),
            )
            .unwrap(),
        );
        let definition = Definition::of(&curve, curve.natural_interval()).unwrap();
        let (fixed, _) = definition.jet(Ivl::from(1.0), None).unwrap();
        let mut roots = Vec::new();
        for piece in definition.pieces(curve.natural_interval()).unwrap() {
            roots.extend(
                isolate(&kernel, piece.range, |t| {
                    edge_jet(&definition, piece, fixed, Vec3::Z, t)
                })
                .unwrap(),
            );
        }
        assert_eq!(roots.len(), 1);
        let root = roots[0].parameter;
        assert!(root.lo < 0.381209626 && root.hi > 0.381209624);
        assert!(root.hi - root.lo < 1e-10);
        // This tests the original native parameter map, not STEP degrees.
        let (point, _) = definition.jet(root, None).unwrap();
        let distance = norm_squared(project(sub(point, fixed), V3::from_vec3(Vec3::Z))).sqrt();
        assert!(distance.lo <= 6.327472493603755 && distance.hi >= 6.327472493603755);
        let piece = definition
            .pieces(curve.natural_interval())
            .unwrap()
            .into_iter()
            .find(|piece| piece.range.lo <= root.lo && root.hi <= piece.range.hi)
            .unwrap();
        let objective = edge_value_over(&definition, piece, fixed, Vec3::Z, root).unwrap();
        assert!(objective.lo <= 6.327472493603755 && objective.hi >= 6.327472493603755);
        assert!(objective.hi - objective.lo < 2e-11);
    }
}
