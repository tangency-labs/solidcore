//! A B-spline curve's speed `|C'|` on one knot span, as an integrand (not upstream): whether
//! the span is a point, exactly, and the bound on the rounding error of evaluating the
//! speed there ([`BSplineCurve::speed_rounding_bound`]), which a length integral's step
//! acceptance is floored at ([`crate::mass::integrator::OdeOptions::rounding`]).

use super::BSplineCurve;
use crate::error::{Result, SolidError};
use crate::geom::extent1d::Extent1d;
use crate::mass::integrator::gamma;
use crate::nurbs::math::{CPoint, NOW, NOZ};

const OP: &str = "nurbs::b_spline_curve::speed";

/// The span's control points and the length of its knot span.
struct Span<'a> {
    points: &'a [CPoint],
    length: f64,
    degree: usize,
}

impl BSplineCurve {
    /// The control points `P_(j-p) ..= P_j` of the knot span `[t_j, t_(j+1))` holding the
    /// middle of `interval` (which lies within one knot span), and that span's length.
    fn span_of(&self, interval: &Extent1d) -> Result<Span<'_>> {
        let cur = self.nurb_ref()?;
        let degree = cur.degree();
        let knots = cur.knots();
        let points = cur.cpts();
        let mid = 0.5 * (interval.min() + interval.max());
        let missing = || SolidError::invalid_input(OP, "interval outside the knot spans");
        // the last knot at or below the middle, at most the last control point's index
        let j = knots
            .iter()
            .rposition(|&t| t <= mid)
            .ok_or_else(missing)?
            .min(points.len().checked_sub(1).ok_or_else(missing)?);
        let (Some(&lo), Some(&hi)) = (knots.get(j), knots.get(j + 1)) else {
            return Err(missing());
        };
        let first = j.checked_sub(degree).ok_or_else(missing)?;
        Ok(Span {
            points: points.get(first..=j).ok_or_else(missing)?,
            length: hi - lo,
            degree,
        })
    }

    /// Whether the curve is a single point on the knot span holding `interval`: all of
    /// that span's control points are the same Euclidean point, exactly. For a rational
    /// curve (homogeneous `(x w, y w, z w, w)`) that is `h_i w_0 = h_0 w_i` for every
    /// coordinate, compared as exact products (each `a b` as `fl(a b)` plus its exact
    /// remainder, which is unique for a real value): every weight then multiplies the
    /// same point, and `C = sum R_i P_i` is that point whatever the weights.
    ///
    /// # Errors
    ///
    /// `SM_ERR` without a NURB; `InvalidInput` when `interval` is outside the knot
    /// spans.
    pub fn is_point_on(&self, interval: &Extent1d) -> Result<bool> {
        let span = self.span_of(interval)?;
        let Some(&p0) = span.points.first() else {
            return Ok(false);
        };
        let coords = |p: CPoint| [p.x, p.y, if p.z == NOZ { 0.0 } else { p.z }];
        Ok(span.points.iter().all(|&p| {
            if p.w == NOW || p0.w == NOW {
                return p.w == p0.w && coords(p) == coords(p0);
            }
            coords(p)
                .into_iter()
                .zip(coords(p0))
                .all(|(h, h0)| exact_product(h, p0.w) == exact_product(h0, p.w))
        }))
    }

    /// A bound on the absolute rounding error of the computed speed `|C'(t)|` for `t` in
    /// the knot span `[t_j, t_(j+1))` (length `D`) holding `interval`: with `p` the
    /// degree, `M` the largest Euclidean control point norm of the span and `rho` the
    /// ratio of its largest to its smallest weight,
    ///
    /// `delta = gamma_(7p+10) S M`, `S = 2p / D` (polynomial), `S = 4 p rho / D`
    /// (rational).
    ///
    /// The evaluator forms `C' = sum_i R'_i P_i` over the span's `p + 1` control points,
    /// `R'_i` the derivatives of the (rational) basis.
    /// - `S` bounds `sum_i |R'_i|`: `N'_(i,p) = p (N_(i,p-1) / (t_(i+p) - t_i) -
    ///   N_(i+1,p-1) / (t_(i+p+1) - t_(i+1)))`; each denominator of a function nonzero
    ///   on the span is a support holding the span, so at least `D`, and the `N_(i,p-1)`
    ///   sum to 1: `sum |N'_i| <= 2p / D`. Rationally `R'_i = (w_i N'_i - R_i w') / w`
    ///   with `w = sum N_i w_i >= w_min`, `|w'| <= w_max 2p / D` and `sum R_i = 1`:
    ///   `sum |R'_i| <= 2 (w_max / w_min) 2p / D`.
    /// - Each computed `R'_i` is within `gamma_k` of the sum of the absolute values of its
    ///   terms (Higham, "Accuracy and Stability of Numerical Algorithms", 2nd ed., §3.1):
    ///   the basis triangle's `p` levels each take two knot differences, a division, a
    ///   product and a sum (`5p`), the derivative's coefficient a knot difference, a
    ///   division and two products (4), the weights' sums `p + 1` and the rational
    ///   quotient 4 more; the sum over the `p + 1` points adds `p + 1`:
    ///   `k = 5p + 4 + (p + 1) + 4 + (p + 1) = 7p + 10`.
    /// - So each coordinate of the computed `C'` is within `gamma_k S max|P_i|`, and its
    ///   norm within `gamma_k S M` of the exact `|C'|` (the norm is 1-Lipschitz; its own
    ///   roundings are relative to the sample, counted by the integrator's estimate,
    ///   [`crate::mass::integrator::OdeOptions::rounding`]).
    ///
    /// On a curve whose control points are all equal `C' = 0` exactly and the computed
    /// speed is noise of this size: no length integral can resolve an accuracy below
    /// it.
    ///
    /// # Errors
    ///
    /// `SM_ERR` without a NURB; `InvalidInput` when `interval` is outside the knot
    /// spans, or for a rational span with a weight that is not positive (the bound
    /// needs `w >= w_min > 0`).
    pub fn speed_rounding_bound(&self, interval: &Extent1d) -> Result<f64> {
        let span = self.span_of(interval)?;
        let largest = span
            .points
            .iter()
            .map(|p| p.to_euclid().length())
            .fold(0.0_f64, f64::max);
        let p = span.degree as f64;
        let rational = span.points.iter().any(|q| q.w != NOW);
        let sum_of_basis_derivatives = if rational {
            let (lo, hi) = span
                .points
                .iter()
                .fold((f64::INFINITY, 0.0_f64), |(lo, hi), q| {
                    (lo.min(q.w), hi.max(q.w))
                });
            if lo.is_nan() || lo <= 0.0 {
                return Err(SolidError::invalid_input(OP, "a weight is not positive"));
            }
            4.0 * p * (hi / lo) / span.length
        } else {
            2.0 * p / span.length
        };
        let n = u32::try_from(span.degree)
            .ok()
            .and_then(|d| d.checked_mul(7))
            .and_then(|d| d.checked_add(10))
            .ok_or_else(|| SolidError::invalid_input(OP, "degree out of range"))?;
        let bound = gamma(n) * sum_of_basis_derivatives * largest;
        if bound.is_finite() {
            Ok(bound)
        } else {
            Err(SolidError::invalid_input(OP, "degenerate knot span"))
        }
    }
}

/// `a b` exactly, as `fl(a b)` and the remainder `a b - fl(a b)` (exact by a fused
/// multiply-add, barring underflow).
fn exact_product(a: f64, b: f64) -> (f64, f64) {
    let p = a * b;
    (p, a.mul_add(b, -p))
}
