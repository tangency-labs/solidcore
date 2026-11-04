//! Intersections of screw-invariant walls reduce to two-dimensional section roots.
//!
//! A root generates one helix. Its first pcurve is a constant-profile line; the
//! second is the exact fractional-linear change between rational arc parameters.
//! No spatial intersection trace or sampled root certificate is used.
//! Deviations: new construction. Root subdivision is iterative, at most 16,384
//! boxes and depth 48; windings and pcurve spans have fixed resource budgets.

mod cylinder;

use std::f64::consts::TAU;

use crate::blend::fillet_section::interval::Ivl;
use crate::curve::crv_on_surf::CrvOnSurf;
use crate::curve::curve::Curve;
use crate::error::{Result, SolidError};
use crate::geom::extent2d::Extent2d;
use crate::kernel::Kernel;
use crate::math::{Point2, Point3, Vec3};
use crate::nurbs::b_spline_curve::{BSplineCurve, BSplineCurveForm, CanonicalData, KnotType};
use crate::surface::procedural_sweep::screw_offset::{
    IntersectionIsolation, ScrewOffsetSection, SectionIntersection,
};
use crate::surface::procedural_sweep::ProceduralSweep;
use crate::surface::surface::p2::{SurfaceIntersection, TsectCurveType};
use crate::surface::surface::Surface;

const OP: &str = "offset::procedural_screw::intersect";
const MAX_BOXES: usize = 16_384;
const MAX_DEPTH: usize = 48;
const MAX_WINDINGS: i64 = 12_000;

fn bad(reason: &'static str) -> SolidError {
    SolidError::invalid_input(OP, reason)
}

struct Wall<'a> {
    surface: &'a Surface,
    definition: &'a ProceduralSweep,
    section: ScrewOffsetSection,
    source: [f64; 2],
    profile: Ivl,
    angle: Ivl,
    face_domain: Extent2d,
}

impl<'a> Wall<'a> {
    fn new(surface: &'a Surface, domain: &Extent2d) -> Result<Option<Self>> {
        let Some(sweep) = surface.procedural_sweep() else {
            return Ok(None);
        };
        let definition = sweep.definition();
        let Some((section, source)) = definition.screw_offset_section()? else {
            return Ok(None);
        };
        let a = definition.definition_parameters(domain.min());
        let b = definition.definition_parameters(domain.max());
        let box_ = Extent2d::from_points(a, b);
        let [lo, hi] = source;
        let profile = (Ivl::new(box_.u_min(), box_.u_max()) - Ivl::from(lo))
            / (Ivl::from(hi) - Ivl::from(lo));
        let angle = Ivl::new(
            definition.turn_at(box_.v_min()).0,
            definition.turn_at(box_.v_max()).0,
        );
        if ![profile.lo, profile.hi, angle.lo, angle.hi]
            .iter()
            .all(|v| v.is_finite())
            || profile.lo >= profile.hi
            || angle.lo >= angle.hi
        {
            return Err(bad("invalid screw intersection domain"));
        }
        Ok(Some(Self {
            surface,
            definition,
            section,
            source,
            profile,
            angle,
            face_domain: *domain,
        }))
    }

    fn profile_parameter(&self, x: f64) -> f64 {
        let [lo, hi] = self.source;
        lo + (hi - lo) * x
    }
}

/// Returns `None` only for a different geometric class, before any construction.
/// Once classified, a failed certificate is an error, never a spatial-SSI retry.
pub(crate) fn intersect(
    k: &Kernel,
    first: &Surface,
    first_domain: &Extent2d,
    second: &Surface,
    second_domain: &Extent2d,
    tolerance: f64,
) -> Result<Option<SurfaceIntersection>> {
    if first.procedural_sweep().is_some()
        && matches!(second.kind(), crate::surface::surface::SurfaceKind::Cone(_))
    {
        return cylinder::intersect(k, first, first_domain, second, second_domain, tolerance);
    }
    if second.procedural_sweep().is_some()
        && matches!(first.kind(), crate::surface::surface::SurfaceKind::Cone(_))
    {
        let mut result =
            cylinder::intersect(k, second, second_domain, first, first_domain, tolerance)?;
        if let Some(out) = &mut result {
            std::mem::swap(&mut out.surface1_uv_curves, &mut out.surface2_uv_curves);
        }
        return Ok(result);
    }
    let (Some(sa), Some(sb)) = (first.procedural_sweep(), second.procedural_sweep()) else {
        return Ok(None);
    };
    if !sa.definition().is_normal_offset() && !sb.definition().is_normal_offset() {
        return Ok(None);
    }
    if !sa.definition().is_screw_invariant() || !sb.definition().is_screw_invariant() {
        return Ok(None);
    }
    let (Some(ma), Some(mb)) = (sa.definition().screw(), sb.definition().screw()) else {
        return Ok(None);
    };
    if ma.origin() != mb.origin() || ma.axis() != mb.axis() || ma.pitch() != mb.pitch() {
        return Ok(None);
    }
    let (Some(a), Some(b)) = (
        Wall::new(first, first_domain)?,
        Wall::new(second, second_domain)?,
    ) else {
        return Err(bad(
            "common-screw intersection requires invariant affine source sections",
        ));
    };
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(bad(
            "a screw intersection needs a positive finite tolerance",
        ));
    }
    let reference = a.section.evaluate(a.profile.mid())?.0 - ma.origin();
    let radial = reference - ma.axis() * reference.dot(ma.axis());
    if a.section == b.section {
        return coincident(k, &a, &b, radial, tolerance).map(Some);
    }
    let bounds_a = a
        .section
        .meridian_bounds(a.profile.lo, a.profile.hi, radial)?;
    let bounds_b = b
        .section
        .meridian_bounds(b.profile.lo, b.profile.hi, radial)?;
    let pitch = ma.pitch();
    let winding = if pitch == 0.0 {
        Ivl::from(0.0)
    } else {
        (bounds_a.height - bounds_b.height) / Ivl::from(pitch)
    };
    if !winding.lo.is_finite()
        || !winding.hi.is_finite()
        || winding.lo < -(MAX_WINDINGS as f64)
        || winding.hi > MAX_WINDINGS as f64
    {
        return Err(bad("screw intersection exceeds its winding budget"));
    }
    let mut out = SurfaceIntersection::default();
    let mut remaining = MAX_BOXES;
    for n in (winding.lo.ceil() as i64)..=(winding.hi.floor() as i64) {
        let shift = Ivl::from(pitch) * Ivl::from(n as f64);
        let roots = roots(
            [&a.section, &b.section],
            [a.profile, b.profile],
            radial,
            shift,
            tolerance / 64.0,
            &mut remaining,
        )?;
        for root in roots {
            let orbit = isolated_orbit(&a, &b, radial, n, root)?;
            append_helix(k, &a, &b, radial, orbit, tolerance, &mut out)?;
        }
    }
    Ok(Some(out))
}

#[derive(Clone, Copy)]
struct Orbit {
    parameters: [f64; 2],
    winding: i64,
    deviation: f64,
    kind: TsectCurveType,
}

fn isolated_orbit(
    a: &Wall<'_>,
    b: &Wall<'_>,
    radial: Vec3,
    winding: i64,
    root: SectionIntersection,
) -> Result<Orbit> {
    let box_a = a
        .section
        .meridian_bounds(root.first.lo, root.first.hi, radial)?;
    let box_b = b
        .section
        .meridian_bounds(root.second.lo, root.second.hi, radial)?;
    let dr = box_a.radius - box_b.radius;
    let dz = box_a.height
        - box_b.height
        - Ivl::from(winding as f64)
            * Ivl::from(
                a.definition
                    .screw()
                    .ok_or_else(|| bad("missing screw motion"))?
                    .pitch(),
            );
    let deviation = (Ivl::from(dr.mag()) * Ivl::from(dr.mag())
        + Ivl::from(dz.mag()) * Ivl::from(dz.mag()))
    .sqrt()
    .hi;
    Ok(Orbit {
        parameters: [root.first.mid(), root.second.mid()],
        winding,
        deviation,
        kind: TsectCurveType::Crossing,
    })
}

// Identical retained sections prove coincidence directly. Boundaries of the
// overlap rectangle are exact section curves and helices; no root is invented.
fn coincident(
    k: &Kernel,
    a: &Wall<'_>,
    b: &Wall<'_>,
    radial: Vec3,
    tolerance: f64,
) -> Result<SurfaceIntersection> {
    let rounding = (Ivl::from(a.definition.floor()?) + Ivl::from(b.definition.floor()?)).hi;
    if !rounding.is_finite() || rounding > tolerance {
        return Err(bad(
            "coincident screw boundaries exceed the requested numerical tolerance",
        ));
    }
    let lo = a.profile.lo.max(b.profile.lo);
    let hi = a.profile.hi.min(b.profile.hi);
    let first = a.angle.lo.max(b.angle.lo);
    let last = a.angle.hi.min(b.angle.hi);
    let mut out = SurfaceIntersection::default();
    if lo > hi || first > last {
        return Ok(out);
    }
    if first < last {
        for x in [lo, hi] {
            append_helix(
                k,
                a,
                b,
                radial,
                Orbit {
                    parameters: [x, x],
                    winding: 0,
                    deviation: 0.0,
                    kind: TsectCurveType::RegionBoundary,
                },
                tolerance,
                &mut out,
            )?;
            if lo == hi {
                break;
            }
        }
    }
    if lo < hi {
        for angle in [first, last] {
            let first_uv = section_pcurve(a, lo, hi, angle)?;
            let second_uv = section_pcurve(b, lo, hi, angle)?;
            append_pair(
                k,
                a,
                first_uv,
                second_uv,
                TsectCurveType::RegionBoundary,
                rounding,
                &mut out,
            )?;
            if first == last {
                break;
            }
        }
    }
    Ok(out)
}

fn section_pcurve(wall: &Wall<'_>, from: f64, to: f64, angle: f64) -> Result<Curve> {
    let t = wall.definition.parameter_at_angle(angle)?;
    let points = [from, to]
        .into_iter()
        .map(|x| {
            let p = wall
                .definition
                .face_parameters(Point2::new(wall.profile_parameter(x), t))?;
            Ok(Point3::new(p.x, p.y, 0.0))
        })
        .collect::<Result<Vec<_>>>()?;
    uv_curve(&points, None, &[0.0, 1.0], &[2, 2])
}

fn append_pair(
    k: &Kernel,
    a: &Wall<'_>,
    first_uv: Curve,
    second_uv: Curve,
    kind: TsectCurveType,
    deviation: f64,
    out: &mut SurfaceIntersection,
) -> Result<()> {
    let curve = Curve::from(CrvOnSurf::new(
        k,
        first_uv.copy(k)?,
        a.surface.copy(k)?,
        Some(a.face_domain),
        3,
    )?);
    out.curves_3d.push(curve);
    out.surface1_uv_curves.push(Some(first_uv));
    out.surface2_uv_curves.push(Some(second_uv));
    out.curve_types.push(kind);
    out.deviations.push(deviation);
    Ok(())
}

fn roots(
    sections: [&ScrewOffsetSection; 2],
    domains: [Ivl; 2],
    radial: Vec3,
    shift: Ivl,
    tolerance: f64,
    remaining: &mut usize,
) -> Result<Vec<SectionIntersection>> {
    let [a, b] = sections;
    let [first, second] = domains;
    let mut pending = vec![(first, second, 0_usize)];
    let mut result = Vec::new();
    while let Some((x, y, depth)) = pending.pop() {
        *remaining = remaining
            .checked_sub(1)
            .ok_or_else(|| bad("screw section root isolation exhausted its box budget"))?;
        match a.isolate_intersection_shifted(b, radial, x, y, shift, tolerance)? {
            IntersectionIsolation::NoRoot => {}
            IntersectionIsolation::Unique(root) => {
                if !result.iter().any(|old: &SectionIntersection| {
                    old.first.lo <= root.first.hi
                        && root.first.lo <= old.first.hi
                        && old.second.lo <= root.second.hi
                        && root.second.lo <= old.second.hi
                }) {
                    result.push(root);
                }
            }
            IntersectionIsolation::Undecided { first, second } => {
                if depth >= MAX_DEPTH {
                    return Err(bad("a screw section is tangent, coincident or not isolated within the depth budget"));
                }
                if first.hi - first.lo >= second.hi - second.lo {
                    let m = first.mid();
                    if m <= first.lo || m >= first.hi {
                        return Err(bad("screw root interval cannot be subdivided"));
                    }
                    pending.push((Ivl::new(m, first.hi), second, depth + 1));
                    pending.push((Ivl::new(first.lo, m), second, depth + 1));
                } else {
                    let m = second.mid();
                    if m <= second.lo || m >= second.hi {
                        return Err(bad("screw root interval cannot be subdivided"));
                    }
                    pending.push((first, Ivl::new(m, second.hi), depth + 1));
                    pending.push((first, Ivl::new(second.lo, m), depth + 1));
                }
            }
        }
    }
    Ok(result)
}

fn uv_curve(
    points: &[Point3],
    weights: Option<&[f64]>,
    knots: &[f64],
    multiplicities: &[usize],
) -> Result<Curve> {
    // Adjacent rational pieces may have differently scaled homogeneous endpoints.
    // Match those scales before sharing the endpoint. Otherwise the canonical
    // constructor's removal of a full knot averages unequal weights and changes
    // the fractional-linear parameter map even though the Euclidean joins agree.
    let mut joined_points = Vec::new();
    let mut joined_weights = Vec::new();
    let mut joined_multiplicities = Vec::new();
    let join = knots.len() > 2 && multiplicities.iter().all(|&m| m == 2);
    if join {
        if points.len() != 2 * (knots.len() - 1) {
            return Err(bad("incompatible piecewise linear pcurve coefficients"));
        }
        for (i, pair) in points.chunks_exact(2).enumerate() {
            let [start, end] = *pair else { continue };
            let weight = |j: usize| -> Result<f64> {
                weights
                    .map_or(Some(1.0_f64), |w| w.get(j).copied())
                    .filter(|w| w.is_finite() && *w > 0.0)
                    .ok_or_else(|| bad("invalid paired pcurve weight"))
            };
            let wa = weight(2 * i)?;
            let wb = weight(2 * i + 1)?;
            let scale = joined_weights.last().copied().unwrap_or(wa) / wa;
            let end_weight = wb * scale;
            if !end_weight.is_finite() || end_weight <= 0.0 {
                return Err(bad("paired pcurve homogeneous scales overflowed"));
            }
            if i == 0 {
                joined_points.push(start);
                joined_weights.push(wa);
            }
            joined_points.push(end);
            joined_weights.push(end_weight);
        }
        joined_multiplicities = (0..knots.len())
            .map(|i| if i == 0 || i + 1 == knots.len() { 2 } else { 1 })
            .collect();
    }
    let (points, weights, multiplicities) = if join {
        (
            joined_points.as_slice(),
            Some(joined_weights.as_slice()),
            joined_multiplicities.as_slice(),
        )
    } else {
        (points, weights, multiplicities)
    };
    Ok(Curve::from(BSplineCurve::create_canonical(
        &CanonicalData {
            dim: 2,
            degree: 1,
            control_points: points,
            form: BSplineCurveForm::Unspecified,
            knot_multiplicities: multiplicities,
            knots,
            knot_type: KnotType::Unspecified,
            weights,
            trim: None,
        },
    )?))
}

fn span(definition: &ProceduralSweep, angle: f64) -> Result<(f64, f64)> {
    let stations = definition.stations();
    let i = stations
        .partition_point(|&v| v <= angle)
        .saturating_sub(1)
        .min(stations.len().saturating_sub(2));
    match (stations.get(i), stations.get(i + 1)) {
        (Some(&a), Some(&b)) if a < b => Ok((a, b)),
        _ => Err(bad("missing screw arc span")),
    }
}

#[allow(clippy::too_many_arguments)]
fn append_helix(
    k: &Kernel,
    a: &Wall<'_>,
    b: &Wall<'_>,
    radial: Vec3,
    orbit: Orbit,
    tolerance: f64,
    out: &mut SurfaceIntersection,
) -> Result<()> {
    let [x, y] = orbit.parameters;
    let qa = a.section.meridian(x, radial)?;
    let qb = b.section.meridian(y, radial)?;
    let shift = qa.phase - qb.phase + orbit.winding as f64 * TAU;
    let lo = a.angle.lo.max(b.angle.lo - shift);
    let hi = a.angle.hi.min(b.angle.hi - shift);
    if lo >= hi {
        return Ok(());
    }
    let mut angles = vec![lo];
    angles.extend(
        a.definition
            .stations()
            .iter()
            .copied()
            .filter(|&v| v > lo && v < hi),
    );
    angles.extend(
        b.definition
            .stations()
            .iter()
            .map(|&v| v - shift)
            .filter(|&v| v > lo && v < hi),
    );
    angles.push(hi);
    angles.sort_by(f64::total_cmp);
    angles.dedup();
    if angles.len() > 24_002 {
        return Err(bad("screw helix exceeds its arc budget"));
    }
    let wa = a.profile_parameter(x);
    let wb = b.profile_parameter(y);
    let mut points = Vec::new();
    let mut weights = Vec::new();
    let mut knots = Vec::new();
    for pair in angles.windows(2) {
        let [from, to] = *pair else { continue };
        let mid = from + (to - from) * 0.5;
        let (a0, a1) = span(a.definition, mid)?;
        let (b0, b1) = span(b.definition, mid + shift)?;
        let am = (a0 + a1) * 0.5;
        let bm = (b0 + b1) * 0.5;
        let at = ((a1 - a0) * 0.25).tan();
        let bt = ((b1 - b0) * 0.25).tan();
        let q = ((am + shift - bm) * 0.5).tan();
        for angle in [from, to] {
            let t = a.definition.parameter_at_angle(angle)?;
            let z = at * (2.0 * (t - a0) / (a1 - a0) - 1.0);
            let denominator = 1.0 - q * z;
            let other_t = bm + (b1 - b0) / (2.0 * bt) * (z + q) / denominator;
            if !denominator.is_finite() || denominator <= 0.0 || !other_t.is_finite() {
                return Err(bad("singular screw arc parameter correspondence"));
            }
            let uv = b.definition.face_parameters(Point2::new(wb, other_t))?;
            points.push(Point3::new(uv.x, uv.y, 0.0));
            weights.push(denominator);
            if angle == from {
                knots.push(t);
            }
        }
    }
    let begin = a.definition.parameter_at_angle(lo)?;
    let end = a.definition.parameter_at_angle(hi)?;
    knots.push(end);
    let multiplicities = vec![2; knots.len()];
    let second_uv = uv_curve(&points, Some(&weights), &knots, &multiplicities)?;
    let first_points = [begin, end]
        .map(|t| a.definition.face_parameters(Point2::new(wa, t)))
        .into_iter()
        .collect::<Result<Vec<_>>>()?
        .into_iter()
        .map(|p| Point3::new(p.x, p.y, 0.0))
        .collect::<Vec<_>>();
    let first_uv = uv_curve(&first_points, None, &[begin, end], &[2, 2])?;
    let deviation = (Ivl::from(orbit.deviation)
        + Ivl::from(a.definition.floor()?)
        + Ivl::from(b.definition.floor()?))
    .hi;
    if !deviation.is_finite() || deviation > tolerance {
        return Err(bad(
            "the isolated screw helix does not meet the requested tolerance",
        ));
    }
    append_pair(k, a, first_uv, second_uv, orbit.kind, deviation, out)
}
