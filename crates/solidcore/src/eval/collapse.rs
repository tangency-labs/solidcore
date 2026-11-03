//! Which parameter direction a surface is collapsed along at a point: the one classifier of
//! poles ([`uv::degenerate`](super::uv::degenerate), [`uv::degenerate_within`](super::uv::degenerate_within),
//! [`SurfaceEval::pole`](super::SurfaceEval::pole), [`SurfaceEval::pole_within`](super::SurfaceEval::pole_within)).
//!
//! A surface is collapsed along `d` at `(u, v)` when its iso-line through the point along `d`
//! is one point: moving along `d` does not move the point, to any order. That is decided on the
//! iso-line itself, exactly, never from derivatives at the point (which vanish at an isolated
//! singular point too, and leave a pole of high order undecided):
//!
//! - A rational B-spline's iso-line along `v` at `u0` is the B-spline whose control points are
//!   the columns' points at `u0` (de Boor along `u`), and a B-spline with positive weights is
//!   one point exactly when its control points coincide: this holds at any degree. The points
//!   are computed with a certified bound on their rounding ([`Row`]), relative to the net's size,
//!   so "coincide" means "coincide to rounding", at any scale and offset.
//! - A revolution's iso-line along its angle is a circle of the profile point's distance from
//!   the axis; along its profile it is the turned profile, one point when the profile's control
//!   points coincide. The affine map after the turn scales a movement by at most its
//!   Frobenius norm.
//! - The analytic forms ([`SurfaceEval`](super::SurfaceEval)) give their circles' radii in closed form.
//!
//! Within `reach` (a boundary point that lies on a pole to within a tolerance): the iso-line
//! lies in the convex hull of its control points, so it moves by no more than twice their
//! spread from the first, plus rounding; a circle moves by its diameter. When that certified
//! movement is within `reach` the direction is collapsed within reach. Never `speed × span`,
//! which near an apex (where the speed grows from zero) reads both directions as collapsed.
//!
//! What cannot be decided is refused with a typed error: non-finite parameters, reach or
//! control points, a parameter off the surface, a malformed net, a direction past v.

#![allow(clippy::neg_cmp_op_on_partial_ord)] // NaN-refusing tests

use super::spline::{onto_knots, EvalCurve, EvalSurface};
use super::vec::VecExt;
use super::EvalError;
use crate::math::Vec3;

/// The unit roundoff of `f64`.
const UNIT: f64 = f64::EPSILON * 0.5;

/// Higham's `γ_n = n u / (1 − n u)`: the relative error bound of `n` rounded operations.
pub(crate) fn gamma(n: usize) -> f64 {
    let nu = n as f64 * UNIT;
    nu / (1.0 - nu)
}

/// The directions a surface is collapsed along at a point.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Collapse {
    /// Moving along either direction moves the point.
    Neither,
    /// Moving along u does not move the point (the iso-line along u is one point).
    U,
    /// Moving along v does not move the point.
    V,
    /// Neither direction moves it: a corner where two collapsed sides meet.
    Both,
}

impl Collapse {
    /// From whether each direction is collapsed.
    #[must_use]
    pub fn of([u, v]: [bool; 2]) -> Collapse {
        match (u, v) {
            (false, false) => Collapse::Neither,
            (true, false) => Collapse::U,
            (false, true) => Collapse::V,
            (true, true) => Collapse::Both,
        }
    }

    /// Whether each direction (u, v) is collapsed.
    #[must_use]
    pub fn flags(self) -> [bool; 2] {
        [
            matches!(self, Collapse::U | Collapse::Both),
            matches!(self, Collapse::V | Collapse::Both),
        ]
    }

    /// Whether direction `d` (0: u, 1: v) is collapsed; [`EvalError::Invalid`] past v.
    pub fn along(self, d: usize) -> Result<bool, EvalError> {
        self.flags()
            .get(d)
            .copied()
            .ok_or_else(|| EvalError::Invalid(format!("no parameter direction {d}")))
    }
}

/// An iso-line's control points as seen from the first: their largest distance from it
/// (`spread`), and a certified bound on that distance's rounding (`err`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Row {
    spread: f64,
    err: f64,
}

impl Row {
    /// One point to rounding, or one that moves by no more than `reach` (certified: the
    /// iso-line lies in its control points' hull, within twice their spread of any of its
    /// points).
    fn collapsed(self, reach: f64, scale: f64) -> bool {
        self.spread <= self.err || 2.0 * scale * (self.spread + self.err) <= reach
    }
}

fn unevaluable(what: &str) -> EvalError {
    EvalError::Unevaluable(format!("pole classification: {what}"))
}

/// The point of the homogeneous B-spline of degree `p` on `knots` over `ctrl` at `t`, by de
/// Boor's algorithm; `None` off the knot range or on malformed data.
fn de_boor(p: usize, knots: &[f64], ctrl: &[[f64; 4]], t: f64) -> Option<[f64; 4]> {
    let n = ctrl.len();
    let t = onto_knots(knots, p, n, t)?;
    // The span: the last k in [p, n − 1] with knots[k] <= t.
    let last = n.checked_sub(1)?;
    let mut k = p;
    while k < last && knots.get(k.checked_add(1)?).is_some_and(|x| *x <= t) {
        k = k.checked_add(1)?;
    }
    let first = k.checked_sub(p)?;
    let mut d: Vec<[f64; 4]> = ctrl.get(first..=k)?.to_vec();
    for r in 1..=p {
        for j in (r..=p).rev() {
            let i = first.checked_add(j)?;
            let lo = *knots.get(i)?;
            let hi = *knots.get(i.checked_add(p)?.checked_add(1)?.checked_sub(r)?)?;
            let span = hi - lo;
            if !(span > 0.0) {
                return None;
            }
            let a = (t - lo) / span;
            let prev = *d.get(j.checked_sub(1)?)?;
            let cur = d.get_mut(j)?;
            for (c, q) in cur.iter_mut().zip(prev) {
                *c = (1.0 - a) * q + a * *c;
            }
        }
    }
    d.get(p).copied()
}

/// The sizes a rounding bound is taken relative to, over homogeneous points `ctrl`: the
/// largest Euclidean point, the largest homogeneous coordinate vector, the largest and the
/// smallest weight. `None` when a coordinate is not finite or a weight is not positive.
fn sizes(ctrl: &[[f64; 4]]) -> Option<(f64, f64, f64, f64)> {
    let (mut m, mut hx, mut wmax, mut wmin) = (0.0f64, 0.0f64, 0.0f64, f64::INFINITY);
    for h in ctrl {
        if !h.iter().all(|x| x.is_finite()) || !(h[3] > 0.0) {
            return None;
        }
        let x = Vec3::new(h[0], h[1], h[2]);
        hx = hx.max(x.len());
        m = m.max((x * (1.0 / h[3])).len());
        wmax = wmax.max(h[3]);
        wmin = wmin.min(h[3]);
    }
    Some((m, hx, wmax, wmin))
}

/// The iso-line control points `pts` (homogeneous, each a de Boor point of degree `p`, or the
/// net's own points at `p = 0`) of a net whose sizes are `net` ([`sizes`]), as a [`Row`].
fn row(pts: &[[f64; 4]], p: usize, net: (f64, f64, f64, f64)) -> Option<Row> {
    let (m, hx, wmax, wmin) = net;
    // de Boor: p levels of convex combinations, each with its coefficient rounded (two
    // subtractions and a division) and its combination (two products, a subtraction and a
    // sum): within γ_(12p) of the largest coordinate per level, summed over the levels.
    let g = gamma(12usize.checked_mul(p)?.checked_add(4)?);
    let floor = wmin - g * wmax;
    if !(floor > 0.0) {
        return None;
    }
    // A Euclidean point x / w: the numerator's and the weight's errors over the weight, plus
    // the division's own rounding.
    let point = g * (3f64.sqrt() * hx + m * wmax) / floor + 2.0 * UNIT * m;
    let euclid = |h: &[f64; 4]| -> Option<Vec3> {
        let x = Vec3::new(h[0], h[1], h[2]) * (1.0 / h[3]);
        (h[3] > 0.0 && x.is_finite()).then_some(x)
    };
    let first = euclid(pts.first()?)?;
    let mut spread = 0.0f64;
    for h in pts {
        spread = spread.max(euclid(h)?.dist(first));
    }
    // Both points' errors, and the distance's own (a difference, three squares, two sums, a
    // square root) on points within the net's size.
    let err = 2.0 * point + gamma(7) * 2.0 * m;
    (spread.is_finite() && err.is_finite()).then_some(Row { spread, err })
}

/// The Frobenius norm of the linear part of a row-major affine map: a bound on how far it
/// moves a vector of unit length.
fn frobenius(m: &[[f64; 4]; 3]) -> f64 {
    m.iter()
        .flat_map(|r| r.iter().take(3))
        .map(|x| x * x)
        .sum::<f64>()
        .sqrt()
        * (1.0 + gamma(12))
}

/// A profile point's distance from the axis through `origin` along the unit `axis`, with a
/// certified bound on its rounding; the profile is evaluated at `s` ([`de_boor`]).
fn radius(c: &EvalCurve, s: f64, origin: Vec3, axis: Vec3) -> Option<Row> {
    let net = sizes(c.ctrl())?;
    let h = de_boor(c.degree(), c.knots(), c.ctrl(), s)?;
    let pt = row(&[h], c.degree(), net)?;
    // The point's own error (its row's `err` is twice it, plus a distance's).
    let point_err = 0.5 * pt.err;
    let x = Vec3::new(h[0], h[1], h[2]) * (1.0 / h[3]);
    let d = x - origin;
    let r = (d - axis * d.dot(axis)).len();
    // The difference, a dot product, a scaled difference and a length, on |x| + |origin|;
    // the axis is unit to rounding.
    let err = point_err + gamma(16) * (x.len() + origin.len());
    (r.is_finite() && err.is_finite()).then_some(Row { spread: r, err })
}

/// Which directions `s` is collapsed along at `uv`, exactly (to rounding), or within `reach`
/// mm (module doc). Refused with a typed error where that cannot be decided.
pub fn collapse(s: &EvalSurface, uv: [f64; 2], reach: f64) -> Result<Collapse, EvalError> {
    if !(uv[0].is_finite() && uv[1].is_finite()) {
        return Err(unevaluable(&format!(
            "parameters ({}, {}) not finite",
            uv[0], uv[1]
        )));
    }
    if !(reach >= 0.0 && reach.is_finite()) {
        return Err(EvalError::Invalid(format!(
            "pole classification: reach {reach} is not a finite length"
        )));
    }
    match s {
        // An unbounded parameter moves the point without end unless its axis is zero.
        EvalSurface::Plane { du, dv, .. } => {
            if !(du.is_finite() && dv.is_finite()) {
                return Err(unevaluable("a plane's axes are not finite"));
            }
            Ok(Collapse::of([du.len() == 0.0, dv.len() == 0.0]))
        }
        EvalSurface::Nurbs(n) => {
            let [u, v] = s.wrap(uv);
            let [pu, pv] = n.degree();
            let [nu, nv] = n.count();
            let [ku, kv] = n.knots();
            let ctrl = n.ctrl();
            let net = sizes(ctrl).ok_or_else(|| {
                unevaluable("a control point is not finite or a weight not positive")
            })?;
            let at = |i: usize, j: usize| -> Option<[f64; 4]> {
                ctrl.get(j.checked_mul(nu)?.checked_add(i)?).copied()
            };
            // Along v at u: each column j's point at u (de Boor along u).
            let along_v: Option<Vec<[f64; 4]>> = (0..nv)
                .map(|j| {
                    let col: Option<Vec<[f64; 4]>> = (0..nu).map(|i| at(i, j)).collect();
                    de_boor(pu, ku, &col?, u)
                })
                .collect();
            // Along u at v: each row i's point at v (de Boor along v).
            let along_u: Option<Vec<[f64; 4]>> = (0..nu)
                .map(|i| {
                    let line: Option<Vec<[f64; 4]>> = (0..nv).map(|j| at(i, j)).collect();
                    de_boor(pv, kv, &line?, v)
                })
                .collect();
            let off = || unevaluable(&format!("({}, {}) is off the spline", uv[0], uv[1]));
            let rv = row(&along_v.ok_or_else(off)?, pu, net).ok_or_else(off)?;
            let ru = row(&along_u.ok_or_else(off)?, pv, net).ok_or_else(off)?;
            Ok(Collapse::of([
                ru.collapsed(reach, 1.0),
                rv.collapsed(reach, 1.0),
            ]))
        }
        EvalSurface::Revolved {
            profile,
            origin,
            axis,
            transform,
            swap,
        } => {
            let [a, b] = s.wrap(uv);
            let prof = if *swap { b } else { a };
            let scale = frobenius(transform);
            if !(scale.is_finite() && scale > 0.0) {
                return Err(unevaluable("the revolution's map is not finite or is zero"));
            }
            let net = sizes(profile.ctrl()).ok_or_else(|| {
                unevaluable("a profile control point is not finite or a weight not positive")
            })?;
            // Along the angle: the circle of the profile point's distance from the axis; it
            // moves by its diameter.
            let circle = radius(profile, prof, *origin, *axis).ok_or_else(|| {
                unevaluable(&format!("({}, {}) is off the profile", uv[0], uv[1]))
            })?;
            // A circle moves by its diameter, twice its radius: a row of spread r.
            let turned = circle.collapsed(reach, scale);
            // Along the profile: the turned profile, one point when its control points coincide.
            let line = row(profile.ctrl(), 0, net)
                .ok_or_else(|| unevaluable("the profile's control points"))?;
            let along = line.collapsed(reach, scale);
            let (pu, pv) = if *swap {
                (turned, along)
            } else {
                (along, turned)
            };
            Ok(Collapse::of([pu, pv]))
        }
    }
}

/// [`collapse`] for a circle of radius `r` (in closed form, rounded within `err`) along one
/// direction: one point to rounding, or one whose diameter is within `reach`.
pub(crate) fn circle(r: f64, err: f64, reach: f64) -> Result<bool, EvalError> {
    if !(r.is_finite() && err.is_finite()) {
        return Err(unevaluable("a radius is not finite"));
    }
    let r = r.abs();
    Ok(r <= err || 2.0 * (r + err) <= reach)
}

/// The unit roundoff, for the closed forms' own bounds.
pub(crate) const ROUNDOFF: f64 = UNIT;

#[cfg(test)]
mod tests;
