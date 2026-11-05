//! The input contract of the curve blocks the curve runners generate: B-spline, Hermite
//! and line data as upstream's constructors take them (ported from the C++ reference
//! wrappers' `bspline` / `line` checks, `crates/solidcore-ffi/src/units/curve_*.rs`). The
//! generators redraw until these hold, so their cases stay on the working path; hostile
//! data is covered by the fixture and fuzz files.

use solidcore::nurbs::math::{NOW, NOZ};

/// The size bounds of a unit's curve blocks.
#[derive(Debug, Clone, Copy)]
pub struct Limits {
    /// The largest last control point index `n`.
    pub max_index: usize,
    /// The largest degree `p`.
    pub max_degree: usize,
    /// The largest magnitude of a coordinate, knot or scalar.
    pub max_coord: f64,
}

/// Whether `v` is finite and at most `max` in magnitude.
#[must_use]
pub fn finite(v: f64, max: f64) -> bool {
    v.is_finite() && v.abs() <= max
}

/// Whether every coordinate of `v` is [`finite`].
#[must_use]
pub fn vec3(v: &[f64; 3], max: f64) -> bool {
    v.iter().all(|&c| finite(c, max))
}

/// The squared length of `v`.
#[must_use]
pub fn len_sq(v: [f64; 3]) -> f64 {
    v.iter().map(|c| c * c).sum()
}

/// Whether `v` is a [`vec3`] of non-zero length.
#[must_use]
pub fn nonzero3(v: &[f64; 3], max: f64) -> bool {
    vec3(v, max) && len_sq(*v) >= 1e-12
}

/// Whether `x` and `y` are non-zero and far from parallel (`|x × y|² ≥ 1e-6 |x|² |y|²`).
#[must_use]
pub fn frame(x: &[f64; 3], y: &[f64; 3], max: f64) -> bool {
    let cross = [
        x[1] * y[2] - x[2] * y[1],
        x[2] * y[0] - x[0] * y[2],
        x[0] * y[1] - x[1] * y[0],
    ];
    nonzero3(x, max) && nonzero3(y, max) && len_sq(cross) >= 1e-6 * len_sq(*x) * len_sq(*y)
}

/// The Euclidean point of a stored control point (`z = NOZ` is 2D, `w = NOW` polynomial).
#[must_use]
pub fn euclid([x, y, z, w]: [f64; 4]) -> [f64; 3] {
    let z = if z == NOZ { 0.0 } else { z };
    if w == NOW {
        [x, y, z]
    } else {
        [x / w, y / w, z / w]
    }
}

/// What a B-spline block that keeps the contract tells the rest of a case.
#[derive(Debug, Clone)]
pub struct BSpline {
    /// The first knot.
    pub a: f64,
    /// The last knot.
    pub b: f64,
    /// The distinct knots, in order.
    pub knots: Vec<f64>,
    /// The Euclidean first control point.
    pub start: [f64; 3],
    /// The Euclidean last control point.
    pub end: [f64; 3],
}

/// A B-spline of control points `cpts` (stored `[x, y, z, w]`), degree `p` and knots `u`, if
/// it keeps the contract: `1 ≤ p ≤ n` within the limits and `n + p + 2` knots, one
/// dimension and rationality, finite control points with weights above `1e-3` and
/// Euclidean points within range, clamped non-decreasing knots within range over a
/// parameter range of at least `1e-6` relative, and no interior knot of multiplicity above
/// `p`.
#[must_use]
pub fn bspline(cpts: &[[f64; 4]], p: usize, u: &[f64], lim: &Limits) -> Option<BSpline> {
    let n = cpts.len().checked_sub(1)?;
    let m = u.len().checked_sub(1)?;
    if n > lim.max_index || p > lim.max_degree || m > lim.max_index + lim.max_degree + 1 {
        return None;
    }
    if p == 0 || p > n || m != n + p + 1 {
        return None;
    }
    let [_, _, z0, w0] = *cpts.first()?;
    for &[x, y, z, w] in cpts {
        if (z == NOZ) != (z0 == NOZ) || (w == NOW) != (w0 == NOW) {
            return None;
        }
        if !x.is_finite()
            || !y.is_finite()
            || (z != NOZ && !z.is_finite())
            || (w != NOW && !w.is_finite())
        {
            return None;
        }
        if w != NOW && w <= 1e-3 {
            return None;
        }
        if !euclid([x, y, z, w])
            .iter()
            .all(|&c| finite(c, lim.max_coord))
        {
            return None;
        }
    }
    if !u.iter().all(|&k| finite(k, lim.max_coord)) {
        return None;
    }
    let (a, b) = clamped(u, p)?;
    let mut knots: Vec<f64> = Vec::new();
    for &k in u {
        if knots.last() != Some(&k) {
            knots.push(k);
        }
    }
    Some(BSpline {
        a,
        b,
        knots,
        start: euclid(*cpts.first()?),
        end: euclid(*cpts.last()?),
    })
}

/// The parameter range `(a, b)` of knots `u` of degree `p`, if they are non-decreasing and
/// clamped (`p + 1` equal knots at each end and no more), span at least `1e-6` relative, and
/// have no interior knot of multiplicity above `p`.
#[must_use]
pub fn clamped(u: &[f64], p: usize) -> Option<(f64, f64)> {
    let m = u.len().checked_sub(1)?;
    if u.windows(2).any(|w| w.first() > w.last()) || m < 2 * p + 1 {
        return None;
    }
    let (&a, &b) = (u.first()?, u.last()?);
    let head = u.get(..=p)?;
    let tail = u.get(m - p..)?;
    let inner = u.get(p + 1..m - p)?;
    if a >= b || head.iter().any(|&k| k != a) || tail.iter().any(|&k| k != b) {
        return None;
    }
    if b - a < 1e-6 * a.abs().max(b.abs()).max(1.0) {
        return None;
    }
    if inner.iter().any(|&k| k == a || k == b) {
        return None;
    }
    let mut run = 1;
    for w in inner.windows(2) {
        run = if w.first() == w.last() { run + 1 } else { 1 };
        if run > p {
            return None;
        }
    }
    Some((a, b))
}

/// The domain `[min, max]` of a line through `p` along `v` over `domain` with `scale`, if it
/// keeps the contract: finite data, a non-zero direction, a domain of at least `1e-6` and a
/// scale above `1e-6`.
#[must_use]
pub fn line(
    p: &[f64; 3],
    v: &[f64; 3],
    domain: [f64; 2],
    scale: f64,
    max: f64,
) -> Option<(f64, f64)> {
    let [d0, d1] = domain;
    let ok = vec3(p, max)
        && nonzero3(v, max)
        && finite(d0, max)
        && finite(d1, max)
        && (d1 - d0).abs() >= 1e-6
        && finite(scale, max)
        && scale > 1e-6;
    ok.then(|| (d0.min(d1), d0.max(d1)))
}

/// Whether the four vectors of a Hermite segment are [`vec3`]s.
#[must_use]
pub fn hermite(vs: [&[f64; 3]; 4], max: f64) -> bool {
    vs.iter().all(|v| vec3(v, max))
}

/// Whether one of at most `max_transforms` transforms applied to a curve keeps the
/// contract: an origin within `max / max_transforms`, two non-parallel axes ([`frame`]) and,
/// when given, scales in `[1e-3, 1]`.
#[must_use]
pub fn transform(
    o: &[f64; 3],
    x: &[f64; 3],
    y: &[f64; 3],
    scale: Option<[f64; 3]>,
    max: f64,
    max_transforms: u64,
) -> bool {
    let s = scale.unwrap_or([1.0; 3]);
    vec3(o, max)
        && o.iter().all(|c| c.abs() <= max / max_transforms as f64)
        && frame(x, y, max)
        && vec3(&s, max)
        && (scale.is_none() || s.iter().all(|c| (1e-3..=1.0).contains(&c.abs())))
}

/// Whether a bounded plane `SmPlane(o, x, y, scale, dom)` keeps the contract: finite data,
/// axes of squared length in `[1e-6, 1e6]` far from parallel, scales in `[1e-3, 1e3]`, and
/// an ordered domain of sides at least `1e-6` within `1e4`.
#[must_use]
pub fn plane(
    o: &[f64; 3],
    x: &[f64; 3],
    y: &[f64; 3],
    scale: [f64; 2],
    dom: [f64; 4],
    max: f64,
) -> bool {
    let [lx, ly] = [len_sq(*x), len_sq(*y)];
    vec3(o, max)
        && vec3(x, max)
        && vec3(y, max)
        && (1e-6..=1e6).contains(&lx)
        && (1e-6..=1e6).contains(&ly)
        && frame(x, y, max)
        && scale
            .iter()
            .all(|&s| finite(s, max) && (1e-3..=1e3).contains(&s))
        && dom.iter().all(|&v| finite(v, max) && v.abs() <= 1e4)
        && dom[2] - dom[0] >= 1e-6
        && dom[3] - dom[1] >= 1e-6
}
