//! Opt-in split of the intersection curves at singular crossings
//! (`Kernel::fixes.ssi_split_curves_at_singular_points`, WR-TBD, tangency-labs/solidcore#393).
//!
//! Upstream: none. Where two intersection branches cross, the surface normals are
//! parallel, `TraceCurve` marches through the point without ending its curve, and
//! `ComputeSingularityPoint` only runs on start points, so `DoIntersection` returns the
//! complete point set with a curve running through the crossing. This pass runs on its
//! final curve list (part 5, [`super::p5`]) when the flag is on, and is not a port.
//!
//! | Rust                                  | Upstream |
//! |---------------------------------------|----------|
//! | [`split_curves_at_singular_points`]   | none (opt-in fix) |
//!
//! Deviations:
//! - **Not upstream.** Off by default; the flag-off path never reaches this file.
//! - **Heuristic.** A crossing is a local minimum of the sine of the angle between the
//!   surface normals along a crossing or near tangent curve, found on [`SAMPLES`]
//!   uniform samples of the curve parameter and refined by golden section. It splits
//!   the curve (and both UV curves) there when the sine is below [`SINGULAR_SINE`], it
//!   is above [`RISE_SINE`] at both neighbouring samples (so a contact of high order,
//!   where the normals stay nearly parallel, is not cut into slivers) and the point is
//!   more than ten times the approximation tolerance from either end. Two crossings
//!   closer than a sample spacing on one curve are not told apart.
//! - **Failures.** A curve whose normals cannot be evaluated, or whose pieces cannot be
//!   built, is left whole; only [`SolidError::Unported`] propagates.

use crate::curve::curve::Curve;
use crate::error::{ignored, Result, SolidError};
use crate::geom::extent1d::Extent1d;
use crate::geom::vector2d::Point2;
use crate::geom::vector3d::Vec3;
use crate::intersect::adv_surface_intersector::InteriorCurves;
use crate::kernel::Kernel;
use crate::math::EFF_ZERO;
use crate::nurbs::b_spline_curve_n_lib::TrimMode;
use crate::nurbs::knots::Side;
use crate::surface::surface::p1::evaluate_1st_derivatives;
use crate::surface::surface::p2::TsectCurveType;
use crate::surface::surface::Surface;

/// Samples (spans) of the curve parameter searched for crossings.
const SAMPLES: usize = 200;

/// The sine of the angle between the normals below which a minimum is a crossing.
const SINGULAR_SINE: f64 = 1e-6;

/// The sine both neighbouring samples must exceed: the angle rises again on both sides.
const RISE_SINE: f64 = 1e-4;

/// Golden section steps refining a minimum (the bracket shrinks to about `1e-13`).
const REFINE_STEPS: usize = 60;

/// A split point is this many approximation tolerances from either end of its curve.
const END_CLEARANCE: f64 = 10.0;

/// The sine of the angle between the unit normals of the two surfaces at the curve
/// parameter `t`, from the UV curves; `None` where a normal is zero.
fn normals_sine(
    k: &Kernel,
    surfaces: [&Surface; 2],
    uv: [&Curve; 2],
    t: f64,
) -> Result<Option<f64>> {
    let mut n = [Vec3::ZERO; 2];
    for ((surface, curve), n) in surfaces.into_iter().zip(uv).zip(&mut n) {
        let p = curve.evaluate_point(k, t)?;
        let (_, du, dv) =
            evaluate_1st_derivatives(k, surface, Point2::new(p.x, p.y), Side::Left, Side::Left)?;
        *n = du.cross(dv);
    }
    let [n1, n2] = n;
    let (l1, l2) = (n1.length(), n2.length());
    if l1 <= EFF_ZERO || l2 <= EFF_ZERO {
        return Ok(None);
    }
    Ok(Some(n1.cross(n2).length() / (l1 * l2)))
}

/// The minimum of `f` on `[a, b]` by golden section: `(parameter, value)`, or `None`
/// where `f` is undefined.
fn golden_minimum(
    f: &mut impl FnMut(f64) -> Result<Option<f64>>,
    mut a: f64,
    mut b: f64,
) -> Result<Option<(f64, f64)>> {
    let r = 0.5 * (5.0_f64.sqrt() - 1.0);
    let mut x1 = b - r * (b - a);
    let mut x2 = a + r * (b - a);
    let (Some(mut f1), Some(mut f2)) = (f(x1)?, f(x2)?) else {
        return Ok(None);
    };
    for _ in 0..REFINE_STEPS {
        if f1 <= f2 {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = b - r * (b - a);
            let Some(v) = f(x1)? else { return Ok(None) };
            f1 = v;
        } else {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + r * (b - a);
            let Some(v) = f(x2)? else { return Ok(None) };
            f2 = v;
        }
    }
    Ok(Some(if f1 <= f2 { (x1, f1) } else { (x2, f2) }))
}

/// The parameters at which the crossing curve `curve` (UV curves `uv`) passes through a
/// singular crossing, increasing; `None` where its normals cannot be evaluated.
fn crossing_parameters(
    k: &Kernel,
    surfaces: [&Surface; 2],
    curve: &Curve,
    uv: [&Curve; 2],
    tol: f64,
) -> Result<Option<Vec<f64>>> {
    let ivl = curve.natural_interval();
    let (t0, t1) = (ivl.min(), ivl.max());
    let at = |i: usize| t0 + (t1 - t0) * i as f64 / SAMPLES as f64;
    let mut f = |t: f64| -> Result<Option<f64>> {
        Ok(ignored(normals_sine(k, surfaces, uv, t))?.flatten())
    };

    // a zero normal (a surface singularity) is never a minimum
    let mut sines = Vec::with_capacity(SAMPLES + 1);
    for i in 0..=SAMPLES {
        let Some(v) = f(at(i))? else {
            return Ok(None);
        };
        sines.push(v);
    }

    let ends = [curve.evaluate_point(k, t0)?, curve.evaluate_point(k, t1)?];
    let mut found: Vec<f64> = Vec::new();
    for (i, w) in sines.windows(3).enumerate() {
        let [before, here, after] = *w else { continue };
        // a local minimum (the first sample of a plateau), rising again on both sides
        if !(here <= before && here < after && before > RISE_SINE && after > RISE_SINE) {
            continue;
        }
        let Some((t, sine)) = golden_minimum(&mut f, at(i), at(i + 2))? else {
            continue;
        };
        if sine >= SINGULAR_SINE {
            continue;
        }
        let Some(p) = ignored(curve.evaluate_point(k, t))? else {
            continue;
        };
        if ends.iter().all(|&e| !p.closer_than(END_CLEARANCE * tol, e)) {
            found.push(t);
        }
    }
    Ok(Some(found))
}

/// A copy of `curve` trimmed to `[lo, hi]`.
fn piece(k: &Kernel, curve: &Curve, lo: f64, hi: f64) -> Result<Curve> {
    let mut p = curve.copy(k)?;
    p.trim(k, &mut Extent1d::from_min_max(lo, hi), TrimMode::Notify)?;
    Ok(p)
}

/// Splits `curve` at `params` (increasing, inside its interval) into pieces.
fn pieces(k: &Kernel, curve: &Curve, params: &[f64]) -> Result<Vec<Curve>> {
    let ivl = curve.natural_interval();
    let mut bounds = vec![ivl.min()];
    bounds.extend_from_slice(params);
    bounds.push(ivl.max());
    bounds
        .windows(2)
        .map(|w| match *w {
            [lo, hi] => piece(k, curve, lo, hi),
            _ => Err(SolidError::sm_err(OP_PIECES)),
        })
        .collect()
}

const OP_PIECES: &str = "intersect::surface_intersector::split_singular::pieces";

/// The pieces of a 3D curve and its two UV curves split at its crossings, or `None`
/// where it has none or cannot be split.
fn split_one(
    k: &Kernel,
    surfaces: [&Surface; 2],
    tol: f64,
    curve: &Curve,
    uv: [&Curve; 2],
) -> Result<Option<[Vec<Curve>; 3]>> {
    let Some(params) = crossing_parameters(k, surfaces, curve, uv, tol)? else {
        return Ok(None);
    };
    if params.is_empty() {
        return Ok(None);
    }
    let build = || -> Result<[Vec<Curve>; 3]> {
        Ok([
            pieces(k, curve, &params)?,
            pieces(k, uv[0], &params)?,
            pieces(k, uv[1], &params)?,
        ])
    };
    ignored(build())
}

/// Splits the crossing and near tangent curves of `set` (3D curves with their UV curves
/// on `surfaces`) at the points where the surface normals are parallel and the curve
/// passes through, as described in the module header. The five arrays stay parallel: a
/// split curve's type and deviation are repeated for each piece.
///
/// # Errors
///
/// [`SolidError::Unported`] from the evaluations; any other failure leaves the curve it
/// concerns whole.
pub fn split_curves_at_singular_points(
    k: &Kernel,
    surfaces: [&Surface; 2],
    tol: f64,
    set: &mut InteriorCurves,
) -> Result<()> {
    let curves = core::mem::take(&mut set.curves_3d);
    let uv1 = core::mem::take(&mut set.surface1_uv_curves);
    let uv2 = core::mem::take(&mut set.surface2_uv_curves);
    let types = core::mem::take(&mut set.curve_types);
    let deviations = core::mem::take(&mut set.deviations);
    let all = curves
        .into_iter()
        .zip(uv1)
        .zip(uv2)
        .zip(types)
        .zip(deviations);
    for ((((curve, u1), u2), ty), dev) in all {
        let crossing = matches!(ty, TsectCurveType::Crossing | TsectCurveType::NearTangent);
        let split = if crossing {
            split_one(k, surfaces, tol, &curve, [&u1, &u2])?
        } else {
            None
        };
        let [c, a, b] = split.unwrap_or_else(|| [vec![curve], vec![u1], vec![u2]]);
        set.curve_types.extend(core::iter::repeat_n(ty, c.len()));
        set.deviations.extend(core::iter::repeat_n(dev, c.len()));
        set.curves_3d.extend(c);
        set.surface1_uv_curves.extend(a);
        set.surface2_uv_curves.extend(b);
    }
    Ok(())
}
