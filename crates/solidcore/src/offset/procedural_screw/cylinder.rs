//! A coaxial cylinder is the vertical line of the same meridian section.
//!
//! The 3D edge remains a helix on the exact screw wall. Its cylinder angular
//! pcurve is a fractional-linear change of rational circle parameters. Height is
//! interpolated as a function of that circle parameter, with the analytic
//! `sup|theta''| * delta_u^2 / 8` bound. No spatial intersection trace is used.
//! Deviations: new construction; at most 24,000 pcurve spans, no recursion.

use super::{
    append_pair, bad, roots, span, uv_curve, Extent2d, Ivl, Kernel, Point2, Point3, Result,
    ScrewOffsetSection, Surface, SurfaceIntersection, TsectCurveType, Vec3, Wall, MAX_BOXES,
    MAX_WINDINGS, TAU,
};
use crate::surface::cone::Cone;
use crate::surface::surface::SurfaceKind;

const MAX_SPANS: usize = 24_000;

struct CircleArc {
    parameter: [f64; 2],
    angle: [f64; 2],
    endpoint_scale: [f64; 2],
    tau: f64,
    second_angle_bound: f64,
}

impl CircleArc {
    // For endpoint homogeneous scales a²,b², a rational quadratic circle has
    // tan((theta-mid)/2) = tau*((a+b)s-a)/(a+(b-a)s), s=(u-u0)/(u1-u0).
    fn parameter_at(&self, theta: f64) -> Result<f64> {
        let [a, b] = self.endpoint_scale;
        let [t0, t1] = self.angle;
        let q = ((theta - (t0 + t1) * 0.5) * 0.5).tan();
        let denominator = self.tau * (a + b) - q * (b - a);
        let s = a * (q + self.tau) / denominator;
        let [u0, u1] = self.parameter;
        let u = u0 + (u1 - u0) * s;
        if !u.is_finite() || denominator == 0.0 {
            return Err(bad("singular cylinder circle parameterization"));
        }
        Ok(u)
    }

    fn angle_at(&self, u: f64) -> f64 {
        let [u0, u1] = self.parameter;
        let [a, b] = self.endpoint_scale;
        let [t0, t1] = self.angle;
        let s = (u - u0) / (u1 - u0);
        let q = self.tau * ((a + b) * s - a) / (a + (b - a) * s);
        (t0 + t1) * 0.5 + 2.0 * q.atan()
    }
}

struct Cylinder<'a> {
    cone: &'a Cone,
    arcs: Vec<CircleArc>,
    angle: Ivl,
    height: Ivl,
    // Cylinder native linear parameter = a + b * cylinder STEP height.
    linear_map: [f64; 2],
    sense: f64,
}

impl<'a> Cylinder<'a> {
    fn new(k: &Kernel, cone: &'a Cone, domain: &Extent2d, axis: Vec3) -> Result<Self> {
        let frame = cone.canonical();
        let sense = if frame.origin.z_axis() == axis {
            1.0
        } else {
            -1.0
        };
        let native = |u: f64, v: f64| {
            if frame.swap_uv {
                Point2::new(v, u)
            } else {
                Point2::new(u, v)
            }
        };
        let ([u0, u1], [v0, v1]) = if frame.swap_uv {
            (
                [domain.v_min(), domain.v_max()],
                [domain.u_min(), domain.u_max()],
            )
        } else {
            (
                [domain.u_min(), domain.u_max()],
                [domain.v_min(), domain.v_max()],
            )
        };
        let p0 = cone.convert_uv_from_nurbs_to_step(k, native(u0, v0))?;
        let p1 = cone.convert_uv_from_nurbs_to_step(k, native(u1, v1))?;
        let angle = Ivl::new(p0.x.min(p1.x).to_radians(), p0.x.max(p1.x).to_radians());
        let height = Ivl::new(p0.y.min(p1.y), p0.y.max(p1.y));
        let b = (v1 - v0) / (p1.y - p0.y);
        let a = v0 - b * p0.y;
        if ![angle.lo, angle.hi, height.lo, height.hi, a, b]
            .iter()
            .all(|v| v.is_finite())
            || angle.lo >= angle.hi
            || height.lo >= height.hi
            || b == 0.0
        {
            return Err(bad("invalid cylinder intersection domain"));
        }
        let circle = cone
            .rev()
            .polar_curve()
            .ok_or_else(|| bad("cylinder has no rational angular parameterization"))?
            .base()
            .canonical()?;
        if circle.degree != 2
            || circle.knots.len() < 2
            || circle.knots.len() > MAX_SPANS + 1
            || circle.control_points.len() != 2 * circle.knots.len() - 1
            || circle.weights.len() != circle.control_points.len()
            || circle
                .knot_multiplicities
                .iter()
                .enumerate()
                .any(|(i, &m)| {
                    m != if i == 0 || i + 1 == circle.knots.len() {
                        3
                    } else {
                        2
                    }
                })
        {
            return Err(bad(
                "cylinder angular parameterization must have quadratic Bezier circle spans",
            ));
        }
        let mut arcs = Vec::new();
        for (i, pair) in circle.knots.windows(2).enumerate() {
            let [from, to] = *pair else { continue };
            let theta0 = cone
                .convert_uv_from_nurbs_to_step(k, native(from, v0))?
                .x
                .to_radians();
            let theta1 = cone
                .convert_uv_from_nurbs_to_step(k, native(to, v0))?
                .x
                .to_radians();
            let weight = |j| {
                circle
                    .weights
                    .get(j)
                    .copied()
                    .ok_or_else(|| bad("missing circle weight"))
            };
            let w0 = weight(2 * i)?;
            let w1 = weight(2 * i + 2)?;
            if ![from, to, theta0, theta1, w0, w1]
                .iter()
                .all(|v| v.is_finite())
                || from >= to
                || w0 <= 0.0
                || w1 <= 0.0
                || theta0 == theta1
                || (theta1 - theta0).abs() >= std::f64::consts::PI
            {
                return Err(bad("invalid cylinder rational circle span"));
            }
            let ea = Ivl::from(w0).sqrt();
            let eb = Ivl::from(w1).sqrt();
            let tau = ((theta1 - theta0) * 0.25).tan();
            let t = Ivl::checked(tau.abs(), tau.abs());
            let low = Ivl::from(ea.lo.min(eb.lo));
            let difference = Ivl::from((eb - ea).mag());
            let width = Ivl::from(to) - Ivl::from(from);
            let second = (Ivl::from(8.0) * t * ea * eb * difference / (low * low * low)
                + Ivl::from(16.0) * t * t * t * ea * ea * eb * eb / (low * low * low * low))
                / (width * width);
            if !second.hi.is_finite() || low.lo <= 0.0 {
                return Err(bad("cylinder angular derivative bound is not finite"));
            }
            arcs.push(CircleArc {
                parameter: [from, to],
                angle: [theta0, theta1],
                endpoint_scale: [w0.sqrt(), w1.sqrt()],
                tau,
                second_angle_bound: second.hi,
            });
        }
        Ok(Self {
            cone,
            arcs,
            angle,
            height,
            linear_map: [a, b],
            sense,
        })
    }

    fn uv(&self, angular: f64, height: f64) -> Point3 {
        let [a, b] = self.linear_map;
        let linear = a + b * height;
        if self.cone.rev().swap_uv() {
            Point3::new(linear, angular, 0.0)
        } else {
            Point3::new(angular, linear, 0.0)
        }
    }
}

pub(super) fn intersect(
    k: &Kernel,
    surface: &Surface,
    domain: &Extent2d,
    cylinder_surface: &Surface,
    cylinder_domain: &Extent2d,
    tolerance: f64,
) -> Result<Option<SurfaceIntersection>> {
    let Some(sweep) = surface.procedural_sweep() else {
        return Ok(None);
    };
    if !sweep.definition().is_normal_offset() {
        return Ok(None);
    }
    let Some(motion) = sweep.definition().screw() else {
        return Ok(None);
    };
    let SurfaceKind::Cone(cone) = cylinder_surface.kind_with_nurb(k)? else {
        return Ok(None);
    };
    let frame = cone.canonical();
    if frame.semi_angle_deg != 0.0
        || (frame.origin.z_axis() != motion.axis() && frame.origin.z_axis() != -motion.axis())
        || (frame.origin.origin() - motion.origin()).cross(motion.axis()) != Vec3::ZERO
    {
        return Ok(None);
    }
    let Some(wall) = Wall::new(surface, domain)? else {
        return Err(bad(
            "a coaxial cylinder requires an invariant affine screw section",
        ));
    };
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(bad(
            "cylinder intersection needs a positive finite tolerance",
        ));
    }
    let cylinder = Cylinder::new(k, cone, cylinder_domain, motion.axis())?;
    let arm = wall.section.evaluate(wall.profile.mid())?.0 - motion.origin();
    let radial = arm - motion.axis() * arm.dot(motion.axis());
    let radial = radial / radial.length();
    let bounds = wall
        .section
        .meridian_bounds(wall.profile.lo, wall.profile.hi, radial)?;
    // A vertical section wider than the wall's entire height enclosure makes its
    // height parameter auxiliary; the isolated root is the radius equality.
    let margin = (bounds.height.hi - bounds.height.lo).max(bounds.radius.hi);
    let lo = bounds.height.lo - margin;
    let hi = bounds.height.hi + margin;
    let radius = frame.radius;
    let axis = motion.axis();
    let origin = motion.origin();
    let vertical = ScrewOffsetSection::new(
        origin + radial * radius + axis * lo,
        origin + radial * radius + axis * hi,
        motion,
        0.0,
    )?;
    let mut remaining = MAX_BOXES;
    let roots = roots(
        [&wall.section, &vertical],
        [wall.profile, Ivl::new(0.0, 1.0)],
        radial,
        Ivl::from(0.0),
        tolerance / 64.0,
        &mut remaining,
    )?;
    let mut out = SurfaceIntersection::default();
    let mut span_budget = MAX_SPANS;
    for root in roots {
        let x = root.first.mid();
        let point = wall.section.evaluate(x)?.0;
        let meridian = wall.section.meridian(x, radial)?;
        let reference = radial
            .dot(frame.origin.y_axis())
            .atan2(radial.dot(frame.origin.x_axis()));
        let phase = cylinder.sense * meridian.phase + reference;
        let base_height = (point - frame.origin.origin()).dot(frame.origin.z_axis());
        let rise = motion.pitch() / TAU;
        let slope = cylinder.sense * rise;
        let heights = [
            (cylinder.height.lo - base_height) / slope,
            (cylinder.height.hi - base_height) / slope,
        ];
        let theta = Ivl::new(
            wall.angle.lo.max(heights[0].min(heights[1])),
            wall.angle.hi.min(heights[0].max(heights[1])),
        );
        if theta.lo >= theta.hi {
            continue;
        }
        let physical = Ivl::from(cylinder.sense) * theta + Ivl::from(phase);
        let winds = (physical - cylinder.angle) / Ivl::from(TAU);
        if !winds.lo.is_finite()
            || !winds.hi.is_finite()
            || winds.lo < -(MAX_WINDINGS as f64)
            || winds.hi > MAX_WINDINGS as f64
        {
            return Err(bad("cylinder intersection exceeds its winding budget"));
        }
        let enclosure = wall
            .section
            .meridian_bounds(root.first.lo, root.first.hi, radial)?;
        let deviation = (enclosure.radius - Ivl::from(radius)).mag();
        for n in (winds.lo.ceil() as i64)..=(winds.hi.floor() as i64) {
            let shift = phase - n as f64 * TAU;
            let limits = [
                (cylinder.angle.lo - shift) / cylinder.sense,
                (cylinder.angle.hi - shift) / cylinder.sense,
            ];
            let lo = theta.lo.max(limits[0].min(limits[1]));
            let hi = theta.hi.min(limits[0].max(limits[1]));
            if lo < hi {
                append(
                    k,
                    &wall,
                    &cylinder,
                    x,
                    [lo, hi],
                    shift,
                    base_height,
                    rise,
                    deviation,
                    tolerance,
                    &mut span_budget,
                    &mut out,
                )?;
            }
        }
    }
    Ok(Some(out))
}

#[allow(clippy::too_many_arguments)]
fn append(
    k: &Kernel,
    wall: &Wall<'_>,
    cylinder: &Cylinder<'_>,
    x: f64,
    angles: [f64; 2],
    shift: f64,
    base_height: f64,
    rise: f64,
    root_error: f64,
    tolerance: f64,
    budget: &mut usize,
    out: &mut SurfaceIntersection,
) -> Result<()> {
    let [lo, hi] = angles;
    let mut cuts = vec![lo, hi];
    cuts.extend(
        wall.definition
            .stations()
            .iter()
            .copied()
            .filter(|&t| t > lo && t < hi),
    );
    for arc in &cylinder.arcs {
        for t in arc.angle {
            let theta = (t - shift) / cylinder.sense;
            if theta > lo && theta < hi {
                cuts.push(theta);
            }
        }
    }
    cuts.sort_by(f64::total_cmp);
    cuts.dedup();
    let mut points = Vec::new();
    let mut weights = Vec::new();
    let mut knots = Vec::new();
    let mut fit_bound = 0.0_f64;
    for interval in cuts.windows(2) {
        let [from, to] = *interval else { continue };
        let mid = (from + to) * 0.5;
        let beta_mid = cylinder.sense * mid + shift;
        let arc = cylinder
            .arcs
            .iter()
            .find(|arc| {
                let [a, b] = arc.angle;
                beta_mid >= a.min(b) && beta_mid <= a.max(b)
            })
            .ok_or_else(|| bad("cylinder angle is outside its rational circle spans"))?;
        let first_u = arc.parameter_at(cylinder.sense * from + shift)?;
        let last_u = arc.parameter_at(cylinder.sense * to + shift)?;
        let width = Ivl::from(last_u) - Ivl::from(first_u);
        let error = Ivl::from(rise.abs()) * Ivl::from(arc.second_angle_bound) * width * width
            / Ivl::from(8.0);
        let mut pieces = 1_usize;
        while (error / Ivl::from(pieces as f64) / Ivl::from(pieces as f64)).hi > tolerance * 0.25 {
            pieces = pieces
                .checked_mul(2)
                .filter(|&n| n <= MAX_SPANS)
                .ok_or_else(|| bad("cylinder pcurve cannot meet its bounded fit budget"))?;
        }
        *budget = budget
            .checked_sub(pieces)
            .ok_or_else(|| bad("cylinder pcurve span budget exceeded"))?;
        fit_bound = fit_bound.max((error / Ivl::from(pieces as f64) / Ivl::from(pieces as f64)).hi);
        let (a0, a1) = span(wall.definition, mid)?;
        let am = (a0 + a1) * 0.5;
        let at = ((a1 - a0) * 0.25).tan();
        let [b0, b1] = arc.angle;
        let q0 = ((cylinder.sense * am + shift - (b0 + b1) * 0.5) * 0.5).tan();
        let [ea, eb] = arc.endpoint_scale;
        for i in 0..pieces {
            for j in [i, i + 1] {
                let u = first_u + (last_u - first_u) * j as f64 / pieces as f64;
                let theta = if j == 0 {
                    from
                } else if j == pieces {
                    to
                } else {
                    (arc.angle_at(u) - shift) / cylinder.sense
                };
                let t = wall.definition.parameter_at_angle(theta)?;
                let z = cylinder.sense * at * (2.0 * (t - a0) / (a1 - a0) - 1.0);
                let qn = z + q0;
                let qd = 1.0 - q0 * z;
                let denominator = (arc.tau * (ea + eb) * qd - (eb - ea) * qn) * arc.tau.signum();
                if !denominator.is_finite() || denominator <= 0.0 {
                    return Err(bad("singular cylinder paired parameter map"));
                }
                points.push(cylinder.uv(u, base_height + cylinder.sense * rise * theta));
                weights.push(denominator);
                if j == i {
                    knots.push(t);
                }
            }
        }
    }
    let begin = wall.definition.parameter_at_angle(lo)?;
    let end = wall.definition.parameter_at_angle(hi)?;
    knots.push(end);
    let second_uv = uv_curve(&points, Some(&weights), &knots, &vec![2; knots.len()])?;
    let w = wall.profile_parameter(x);
    let first_points = [begin, end]
        .into_iter()
        .map(|t| {
            let p = wall.definition.face_parameters(Point2::new(w, t))?;
            Ok(Point3::new(p.x, p.y, 0.0))
        })
        .collect::<Result<Vec<_>>>()?;
    let first_uv = uv_curve(&first_points, None, &[begin, end], &[2, 2])?;
    let deviation = (Ivl::from(root_error)
        + Ivl::from(fit_bound)
        + Ivl::from(2.0) * Ivl::from(wall.definition.floor()?))
    .hi;
    if !deviation.is_finite() || deviation > tolerance {
        return Err(bad(
            "cylinder helix pcurves are not certified within the requested tolerance",
        ));
    }
    append_pair(
        k,
        wall,
        first_uv,
        second_uv,
        TsectCurveType::Crossing,
        deviation,
        out,
    )
}
