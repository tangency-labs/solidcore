//! The normal at a pole approached along a parameter direction.
//!
//! At a pole ([`super::collapse`] says a direction is collapsed there) `S_u × S_v` vanishes and
//! the normal has no value of its own, but along a line of parameters `(u, v) + t (a, b)`,
//! `t > 0`, the normal of the points it passes through has a limit. Expanding both first
//! derivatives in `t`,
//!
//! `S_u(t) = Σ_m t^m / m! · D^m S_u`, `S_v(t) = Σ_n t^n / n! · D^n S_v`,
//! with `D = a ∂/∂u + b ∂/∂v`, so `D^m S_u = Σ_r C(m, r) a^r b^(m−r) S_(u^(r+1) v^(m−r))`,
//!
//! gives `S_u × S_v = Σ_k t^k N_k`, `N_k = Σ_(m+n=k) (D^m S_u × D^n S_v) / (m! n!)`. For
//! small `t > 0` the normal is along the first `N_k` that does not vanish: that is its limit
//! along the line, exactly (no step is taken). The partials go to the fourth order
//! ([`MAX_PARTIAL_ORDER`]), so `k` runs to three; a pole whose `N_0` to `N_3` all vanish is of
//! higher order and is refused ([`EvalError::Degenerate`]).
//!
//! There is one pole decision, the classifier's: where it says direction `d` is collapsed, the
//! iso-line along `d` is one point, so every pure derivative `∂_d^n S` (n ≥ 1) is identically
//! zero there, and those partials are set to exact zero ([`zero_collapsed`]) before the
//! expansion. A coefficient then vanishes only when it is within its own error: each partial
//! carries the rounding of its evaluation ([`Scale`]: `γ` times the net's size and the basis
//! derivatives' scale along each direction to the partial's orders, `2p / Δ` per order for a
//! spline, the knot span `Δ` holding the point), propagated through the terms of `N_k`:
//! `|N_k| ≤ Σ_m (|D^m S_u| E_n + E_m |D^n S_v| + E_m E_n) / (m! n!) + γ Σ_m |D^m S_u| |D^n S_v| / (m! n!)`.
//! No absolute or cross-order threshold: under parameters scaled by any factor each order of
//! partials, each `N_k` and each bound scale alike, so the order and the limit are the same, at
//! any offset (the caller takes the partials with the pole point at the origin) and size.
use super::spline::{Partials, MAX_PARTIAL_ORDER};
use super::vec::VecExt;
use super::EvalError;
use crate::math::Vec3;

/// The highest order `k` of the expansion the partials reach.
const MAX_ORDER: usize = MAX_PARTIAL_ORDER - 1;

/// The rounding of a coefficient against its terms' magnitudes: the binomial sums, the powers,
/// the cross products and the sum over the terms, generously.
fn cancellation() -> f64 {
    super::collapse::gamma(64)
}

/// The scale of a surface's partials at a point, for their rounding: the net's size (its
/// largest point from the pole, times its weights' spread) and, along u and v, the factor one
/// more derivative can grow them by (`2p / Δ` for a spline's basis, 1 for an angle).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Scale {
    /// The net's size, mm.
    pub size: f64,
    /// Per direction (u, v), the growth of one derivative.
    pub per: [f64; 2],
}

impl Scale {
    /// A bound on the rounding of the partial `∂_u^i ∂_v^j S`: the rational quotient's Leibniz
    /// sums over basis derivatives of these scales, generously.
    fn err(self, i: usize, j: usize) -> f64 {
        super::collapse::gamma(1024)
            * self.size
            * self.per[0].powi(i as i32)
            * self.per[1].powi(j as i32)
    }
}

/// The partials `d` with every pure derivative along a collapsed direction (`collapsed[0]`: u,
/// `[1]`: v) set to exact zero, as the classifier's verdict makes them.
pub fn zero_collapsed(d: &mut Partials, collapsed: [bool; 2]) {
    for n in 1..=MAX_PARTIAL_ORDER {
        if collapsed[0] {
            if let Some(x) = d.get_mut(n).and_then(|row| row.get_mut(0)) {
                *x = Vec3::ZERO;
            }
        }
        if collapsed[1] {
            if let Some(x) = d.get_mut(0).and_then(|row| row.get_mut(n)) {
                *x = Vec3::ZERO;
            }
        }
    }
}

fn binomial(n: usize, k: usize) -> f64 {
    (1..=k).map(|i| (n + 1 - i) as f64 / i as f64).product()
}

fn factorial(n: usize) -> f64 {
    (1..=n).map(|k| k as f64).product()
}

/// `D^m` of the first derivative along `first` (0: `S_u`, 1: `S_v`), from the partials.
fn directional(d: &Partials, first: usize, m: usize, [a, b]: [f64; 2]) -> Option<Vec3> {
    let mut x = Vec3::ZERO;
    for r in 0..=m {
        let (i, j) = if first == 0 {
            (r + 1, m - r)
        } else {
            (r, m - r + 1)
        };
        let s = *d.get(i)?.get(j)?;
        x += s * (binomial(m, r) * a.powi(r as i32) * b.powi((m - r) as i32));
    }
    Some(x)
}

/// The unit normal `S_u × S_v` approached at the point of `d` along the parameter direction
/// `dir` (from the point into the surface's domain), and its order `k`: the first coefficient
/// `N_k` of the expansion that does not vanish by cancellation, normalised (module doc). The
/// partials are the classifier's ([`zero_collapsed`]). [`EvalError::Degenerate`] when every
/// coefficient up to `k = 3` vanishes (a pole of higher order), or `dir` is zero or not finite;
/// [`EvalError::Unevaluable`] when a partial is not finite.
pub fn limit_along(d: &Partials, dir: [f64; 2], scale: Scale) -> Result<(Vec3, usize), EvalError> {
    let degenerate = |why: &str| {
        EvalError::Degenerate(format!("no normal along ({}, {}): {why}", dir[0], dir[1]))
    };
    if !(dir.iter().all(|x| x.is_finite()) && (dir[0] != 0.0 || dir[1] != 0.0)) {
        return Err(degenerate("not a direction"));
    }
    if !d.iter().flatten().all(|p| p.is_finite()) {
        return Err(EvalError::Unevaluable("a partial is not finite".into()));
    }
    let err = |i: usize, j: usize| scale.err(i, j);
    // `D^m` of S_u (first = 0) or S_v (1) and its error bound.
    let directional_err = |first: usize, m: usize| -> f64 {
        (0..=m)
            .map(|r| {
                let (i, j) = if first == 0 {
                    (r + 1, m - r)
                } else {
                    (r, m - r + 1)
                };
                binomial(m, r)
                    * dir[0].abs().powi(r as i32)
                    * dir[1].abs().powi((m - r) as i32)
                    * err(i, j)
            })
            .sum()
    };
    for k in 0..=MAX_ORDER {
        let (mut n, mut bound) = (Vec3::ZERO, 0.0f64);
        for m in 0..=k {
            let su = directional(d, 0, m, dir).ok_or_else(|| degenerate("partials missing"))?;
            let sv = directional(d, 1, k - m, dir).ok_or_else(|| degenerate("partials missing"))?;
            let (eu, ev) = (directional_err(0, m), directional_err(1, k - m));
            let f = 1.0 / (factorial(m) * factorial(k - m));
            n += su.cross(sv) * f;
            bound +=
                (su.len() * ev + eu * sv.len() + eu * ev + cancellation() * su.len() * sv.len())
                    * f;
        }
        let len = n.len();
        if len > bound {
            return Ok((n * (1.0 / len), k));
        }
    }
    Err(degenerate("a pole of order above three"))
}

#[cfg(test)]
mod tests;
