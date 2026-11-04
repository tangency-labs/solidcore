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
//! "Vanishes" is relative to the surface's own scale there: `|N_k|` within `1e-9` of the
//! square of its largest partial up to the fourth order, so the test is the same at any offset
//! and size (a translated pole's first derivatives are rounding, not zero).

use super::spline::{Partials, MAX_PARTIAL_ORDER};
use super::vec::VecExt;
use super::EvalError;
use crate::math::Vec3;

/// The highest order `k` of the expansion the partials reach.
const MAX_ORDER: usize = MAX_PARTIAL_ORDER - 1;

/// How small, against the square of the largest partial, a coefficient is to vanish.
const VANISHES: f64 = 1e-9;

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
/// `dir` (from the point into the surface's domain): the first coefficient `N_k` of the
/// expansion that does not vanish, normalised (module doc). [`EvalError::Degenerate`] when
/// every coefficient up to `k = 3` vanishes (a pole of higher order), or `dir` is zero or not
/// finite; [`EvalError::Unevaluable`] when a partial is not finite.
pub fn normal_along(d: &Partials, dir: [f64; 2]) -> Result<Vec3, EvalError> {
    let degenerate = |why: &str| {
        EvalError::Degenerate(format!("no normal along ({}, {}): {why}", dir[0], dir[1]))
    };
    if !(dir.iter().all(|x| x.is_finite()) && (dir[0] != 0.0 || dir[1] != 0.0)) {
        return Err(degenerate("not a direction"));
    }
    let mut scale = 0.0f64;
    for p in d.iter().flatten().skip(1) {
        if !p.is_finite() {
            return Err(EvalError::Unevaluable("a partial is not finite".into()));
        }
        scale = scale.max(p.len());
    }
    let floor = VANISHES * scale * scale;
    for k in 0..=MAX_ORDER {
        let mut n = Vec3::ZERO;
        for m in 0..=k {
            let su = directional(d, 0, m, dir).ok_or_else(|| degenerate("partials missing"))?;
            let sv = directional(d, 1, k - m, dir).ok_or_else(|| degenerate("partials missing"))?;
            n += su.cross(sv) * (1.0 / (factorial(m) * factorial(k - m)));
        }
        if n.len() > floor {
            return n
                .normalized()
                .ok_or_else(|| degenerate("a coefficient with no direction"));
        }
    }
    Err(degenerate("a pole of order above three"))
}

#[cfg(test)]
mod tests;
