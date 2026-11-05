//! First-order affine forms for the G^n section's certificate (not upstream, #2690): a
//! quantity over a span as `c + sum_i g_i e_i + [-r, r]`, the `e_i` in `[-1, 1]`: the
//! feet's four residuals off their predicted line (`e_0..e_3`) and the span's parameter
//! (`e_4`, `s = m + rad e_4`), in rounded floating point: every operation adds its own
//! rounding error to `r`, so the form holds the exact result for every value of the `e_i`.
//!
//! The certificate's system over a box is a sum of terms that nearly cancel where a
//! ball's offset nearly collapses (`S + r N`, the face's radius of curvature near r): plain
//! intervals add the terms' widths, the form keeps their shared first-order dependence on
//! the box and cancels it, so what is left of a width is the system's own variation and a
//! remainder of the box's square. Linear operations are exact in the form (to rounding);
//! a product puts the product of the two deviations in `r`; `1 / x` and `sqrt x` are their
//! Taylor polynomials of degree one at the middle with the second derivative's bound over
//! the form's interval in `r`.

use core::ops::{Add, Div, Mul, Sub};

use super::interval::{cmax, round_down, round_up, Ivl};
use super::series::GnScalar;

/// The noise symbols of a form: the feet's four residuals and the span's parameter.
pub const NOISE: usize = 5;

/// The span's parameter's noise symbol.
pub const S_NOISE: usize = 4;

/// `c + sum_i g_i e_i + [-r, r]`, `e_i` in `[-1, 1]`; `r` infinite: undefined (the whole
/// line).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Aff {
    /// The middle.
    pub c: f64,
    /// The coefficients of the noise symbols.
    pub g: [f64; NOISE],
    /// The remainder's radius (rounding included).
    pub r: f64,
}

/// `a + b`, rounded up.
#[inline]
fn up_add(a: f64, b: f64) -> f64 {
    round_up(a + b)
}

/// `a b` (non-negative), rounded up.
#[inline]
fn up_mul(a: f64, b: f64) -> f64 {
    round_up(a * b)
}

impl Aff {
    /// The whole line.
    pub const WHOLE: Self = Self {
        c: 0.0,
        g: [0.0; NOISE],
        r: f64::INFINITY,
    };

    /// The number `d`.
    #[must_use]
    pub const fn point(d: f64) -> Self {
        Self {
            c: d,
            g: [0.0; NOISE],
            r: 0.0,
        }
    }

    /// An interval, with no dependence on the noise.
    #[must_use]
    pub fn of_ivl(x: Ivl) -> Self {
        let m = x.mid();
        let r = round_up(cmax(x.hi - m, m - x.lo));
        Self {
            c: m,
            g: [0.0; NOISE],
            r,
        }
        .checked()
    }

    /// The mean value form `mid + sum_j s_j d_j`: `mid` the value at the middle, each
    /// `s_j` enclosing a slope over the region and `d_j` the deviation along it from the
    /// middle, as a form: the slope's middle times the deviation, and its half width
    /// times the deviation's largest size in `r`.
    #[must_use]
    pub fn mean_value(mid: Ivl, terms: &[(Ivl, Self)]) -> Self {
        let mut a = Self::of_ivl(mid);
        for &(s, d) in terms {
            let sm = s.mid();
            let sr = round_up(cmax(s.hi - sm, sm - s.lo));
            a = a + Self::point(sm) * d;
            a.r = up_add(a.r, up_mul(sr, d.radius()));
        }
        a.checked()
    }

    /// The whole line where a part is not finite.
    #[must_use]
    fn checked(self) -> Self {
        let finite = self.c.is_finite() && self.g.iter().all(|g| g.is_finite()) && !self.r.is_nan();
        if finite {
            self
        } else {
            Self::WHOLE
        }
    }

    /// `sum_i |g_i| + r`, rounded up: the largest deviation from the middle.
    #[must_use]
    pub fn radius(&self) -> f64 {
        let mut s = self.r;
        for g in &self.g {
            s = up_add(s, g.abs());
        }
        s
    }

    /// The interval the form ranges over.
    #[must_use]
    pub fn ivl(&self) -> Ivl {
        let s = self.radius();
        if !(s.is_finite() && self.c.is_finite()) {
            return Ivl::whole();
        }
        Ivl::new(round_down(self.c - s), round_up(self.c + s))
    }

    /// The form plus `[-b, b]`.
    #[must_use]
    pub fn widened(mut self, b: f64) -> Self {
        self.r = up_add(self.r, b.abs());
        self.checked()
    }

    /// `f(self)` for f monotone and convex or concave over the form's interval `[a, b]`:
    /// the min-range linearization (de Figueiredo and Stolfi), `f(x) = alpha x + d(x)` with
    /// `alpha` f's slope at the end where it is least steep, so `d` is monotone and lies
    /// between `d(a)` and `d(b)`: the form `alpha self` plus that band, whose range is no
    /// wider than f's own over `[a, b]` (plain intervals'), and which keeps the form's
    /// dependence on the noise.
    fn min_range(self, alpha: f64, d_a: f64, d_b: f64) -> Self {
        let (lo, hi) = if d_a <= d_b { (d_a, d_b) } else { (d_b, d_a) };
        // the band's ends each within a few roundings (they are a few operations each), and
        // alpha within a rounding of f's slope, so d is monotone but within that slope's
        // error over the interval
        let x = self.ivl();
        let slack =
            8.0 * U * cmax(lo.abs(), hi.abs()) + 4.0 * U * alpha.abs() * (x.hi - x.lo) + TINY;
        let (lo, hi) = (lo - slack, hi + slack);
        let zeta = 0.5 * (lo + hi);
        let delta = round_up(cmax(hi - zeta, zeta - lo));
        let mut out = Self::point(alpha) * self + Self::point(zeta);
        out.r = up_add(out.r, delta);
        out.checked()
    }

    /// `1 / self`; the whole line where the form's interval holds 0.
    #[must_use]
    pub fn recip(self) -> Self {
        let x = self.ivl();
        if x.lo > 0.0 {
            // convex, decreasing: least steep at b; d(x) = 1 / x - x / b^2 decreases
            let (a, b) = (x.lo, x.hi);
            let alpha = -1.0 / (b * b);
            return self.min_range(alpha, 1.0 / a - alpha * a, 1.0 / b - alpha * b);
        }
        if x.hi < 0.0 {
            return Self::point(0.0) - (Self::point(0.0) - self).recip();
        }
        Self::WHOLE
    }
}

impl GnScalar for Aff {
    fn of(d: f64) -> Self {
        Self::point(d)
    }

    fn is_positive(self) -> bool {
        self.ivl().lo > 0.0
    }

    /// `sqrt self`; the whole line unless the form's interval is positive.
    fn sqrt(self) -> Self {
        let x = self.ivl();
        let positive = x.lo > 0.0;
        if !positive {
            return Self::WHOLE;
        }
        // concave, increasing: least steep at b; d(x) = sqrt x - x / (2 sqrt b) increases
        let (a, b) = (x.lo, x.hi);
        let (ra, rb) = (
            crate::geom::vector3d::sqrt(a),
            crate::geom::vector3d::sqrt(b),
        );
        let alpha = 0.5 / rb;
        self.min_range(alpha, ra - alpha * a, rb - alpha * b)
    }
}

/// Twice the unit roundoff: each operation's result is within `U |result|` of exact (as
/// [`round_up`]'s).
const U: f64 = 2.3e-16;

/// Below the normal range, each operation's absolute error.
const TINY: f64 = 1.0e-300;

/// A remainder summed in floating point from at most a dozen non-negative terms (each
/// already a bound), each addition and product within `U` of exact: moved up past their
/// rounding, and past `ops` operations' errors below the normal range.
#[inline]
fn settle(r: f64, ops: f64) -> f64 {
    r * (1.0 + 32.0 * U) + ops * TINY
}

impl Add for Aff {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        let c = self.c + b.c;
        let mut g = [0.0; NOISE];
        let mut size = c.abs();
        for ((o, &x), &y) in g.iter_mut().zip(&self.g).zip(&b.g) {
            *o = x + y;
            size += o.abs();
        }
        // each sum's rounding, and the operands' remainders
        let r = settle(self.r + b.r + U * size, (NOISE + 1) as f64);
        Self { c, g, r }.checked()
    }
}

impl Sub for Aff {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        let c = self.c - b.c;
        let mut g = [0.0; NOISE];
        let mut size = c.abs();
        for ((o, &x), &y) in g.iter_mut().zip(&self.g).zip(&b.g) {
            *o = x - y;
            size += o.abs();
        }
        let r = settle(self.r + b.r + U * size, (NOISE + 1) as f64);
        Self { c, g, r }.checked()
    }
}

impl Mul for Aff {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        let c = self.c * b.c;
        let mut g = [0.0; NOISE];
        let (mut la, mut lb, mut lc) = (0.0_f64, 0.0_f64, 0.0_f64);
        for ((o, &x), &y) in g.iter_mut().zip(&self.g).zip(&b.g) {
            *o = self.c * y + b.c * x;
            la += x.abs();
            lb += y.abs();
            lc += o.abs();
        }
        let (ca, cb) = (self.c.abs(), b.c.abs());
        // rounding: the middle's product, each coefficient's two products and their sum
        let rounding = U * (c.abs() + ca * lb + cb * la + lc);
        // the deviations' product, and each middle times the other's remainder
        let rest = ca * b.r + cb * self.r + (la + self.r) * (lb + b.r);
        let r = settle(rounding + rest, (3 * NOISE + 1) as f64);
        Self { c, g, r }.checked()
    }
}

impl Div for Aff {
    type Output = Self;
    // a quotient is the product with the reciprocal's form
    #[allow(clippy::suspicious_arithmetic_impl)]
    fn div(self, b: Self) -> Self {
        self * b.recip()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `c + rad e_i`.
    fn noise(i: usize, c: f64, rad: f64) -> Aff {
        let mut g = [0.0; NOISE];
        g[i] = rad;
        Aff { c, g, r: 0.0 }
    }

    /// Every point of the noise's cube: the form at it lies in the form's interval.
    fn holds(a: &Aff, exact: impl Fn([f64; NOISE]) -> f64) {
        let iv = a.ivl();
        for k in 0..243_u32 {
            let mut e = [0.0; NOISE];
            let mut q = k;
            for x in &mut e {
                *x = f64::from(q % 3) - 1.0;
                q /= 3;
            }
            let v = exact(e);
            assert!(iv.lo <= v && v <= iv.hi, "{v} not in {iv:?}");
        }
    }

    /// The operations hold their exact results over the box, and cancel what plain
    /// intervals cannot: `x - x` is 0, `(1 + x)(1 - x) + x x` is 1, to rounding and the
    /// product's remainder.
    #[test]
    fn affine_forms_enclose_and_cancel() {
        let x = noise(0, 2.0, 0.1);
        let y = noise(1, 0.0, 0.5);
        let ex = |e: [f64; NOISE]| (2.0 + 0.1 * e[0], 0.5 * e[1]);
        let z = (x - x).ivl();
        assert!(z.lo <= 0.0 && z.hi >= 0.0 && z.hi - z.lo < 1e-14, "{z:?}");
        let s = x * y + x / (y + Aff::point(3.0));
        holds(&s, |e| {
            let (a, b) = ex(e);
            a * b + a / (b + 3.0)
        });
        let q = (x + y).sqrt();
        holds(&q, |e| {
            let (a, b) = ex(e);
            (a + b).sqrt()
        });
        let one = Aff::point(1.0);
        let w = ((one + x) * (one - x) + x * x).ivl();
        // the products' remainders, 2 rad^2, are all that is left
        assert!(w.lo <= 1.0 && w.hi >= 1.0 && w.hi - w.lo < 0.05, "{w:?}");
        let plain = {
            let xi = Ivl::new(1.9, 2.1);
            let o = Ivl::from(1.0);
            (o + xi) * (o - xi) + xi * xi
        };
        assert!(plain.hi - plain.lo > 1.0, "{plain:?}");
        assert_eq!(noise(0, 0.0, 1.0).recip(), Aff::WHOLE);
        assert_eq!(noise(0, 0.0, 1.0).sqrt(), Aff::WHOLE);
        // the mean value form of x^2 about 2 over x in [1.9, 2.1]: 4 + [3.8, 4.2] (x - 2)
        let mv = Aff::mean_value(Ivl::from(4.0), &[(Ivl::new(3.8, 4.2), noise(0, 0.0, 0.1))]);
        holds(&mv, |e| (2.0 + 0.1 * e[0]) * (2.0 + 0.1 * e[0]));
        assert_eq!(Aff::of_ivl(Ivl::whole()).ivl(), Ivl::whole());
    }
}
