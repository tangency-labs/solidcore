//! The G^n section's certificate: a proved bound on how far a Hermite span between two
//! stations is from the section of the true jets, from rounded interval arithmetic (the
//! feet over the span enclosed by Krawczyk's operator, their Taylor coefficients order by
//! order, and the Hermite remainder `h^4 / 384 max |c''''|`).
//!
//! For the rational section the blend's distance from its definition (the section of the
//! true jets `c_ik(s)` at every s) is at most `max_i sum_k |c_ik - C_ik|`, `C_ik` the jets'
//! Hermite curves: f is in `[0, 1]` and `|u - 1| <= 1` in `T_i`. On a span of length h a
//! cubic Hermite interpolant of a function's values and derivatives at both ends is off it
//! by at most `h^4 / 384` times the largest fourth derivative. That bound comes from
//! rounded interval arithmetic: the feet over the span are enclosed by Krawczyk's operator,
//! which proves a unique root of G in its box for every s in the span when it maps the box
//! into its interior and the preconditioned Jacobian's norm is below one; the feet's Taylor
//! coefficients, to order 4 in s and n in tau, are then enclosed order by order through the
//! same preconditioned system.
//!
//! Fork: `source/SMLib/src/SmFilletSection.cpp` (tangency-labs/smlib@11ea2d9), the
//! certificate (lines 825-1430); not in the pinned upstream.
//!
//! | Rust                              | Fork (`SmFilletSection.cpp`)            |
//! |-----------------------------------|-----------------------------------------|
//! | [`NurbsCurve`], [`NurbsCurve::of`] | `sm_NurbsCurve`, `sm_NurbsCurveOf`     |
//! | [`curve_derivs_over`]             | `sm_CurveDerivsOver`                    |
//! | [`nurbs_partials_over_piece`]     | `sm_NurbsPartialsOverPiece`             |
//! | [`leader_series_over`]            | `sm_LeaderSeriesOver`                   |
//! | [`GnLaw`]                         | `sm_GnLaw` (0d43b4c)                    |
//! | [`law_series`]                    | `sm_LawSeries` (0d43b4c)                |
//! | [`law_series_over`]               | `sm_LawSeriesOver` (0d43b4c)            |
//! | [`inverse4`]                      | `sm_Inverse4`                           |
//! | [`Enclosure`], [`gn_enclose`]     | `sm_GnEnclose` (its `rF`, `pX`, `pDS`)  |
//! | [`SpanAt`], [`gn_enclose_over`], [`face_partials_over`] | `sm_GnSpanAt`, `sm_GnEncloseOver`, `sm_FacePartialsOver` (c5b6bb7) |
//! | [`gn_fourth_over`]                | `sm_GnFourthOver`                       |
//! | [`gn_certify_span`]               | `sm_GnCertifySpan`                      |
//! | `gn_try!` (private)               | `SM_GN_TRY`                             |
//!
//! Deviations:
//! - **Out-params.** Every function returns its results. `sm_GnEnclose`'s `rbProved` FALSE
//!   is `Ok(None)`; its optional `pX` and `pDS` are always returned (computing them changes
//!   nothing else). `sm_CurveDerivsOver`'s `lPiece = -1` is `None`.
//!   `sm_NurbsPartialsOverPiece` returns fresh partials instead of overwriting the low
//!   orders of the caller's (the orders it leaves are never read).
//! - **`SM_GN_TRY`.** A failed evaluation inside the certificate proves nothing, as in the
//!   fork: `gn_enclose` returns `Ok(None)`, `gn_fourth_over` its unproved (infinite)
//!   bounds, `gn_certify_span` the bound it has so far. Only a truncation past the series'
//!   fixed size (`n > SM_GN_MAXE`, which the fork leaves to its callers) is an `Err`.
//! - **Bounds.** Every loop is the fork's, with its fixed counts (8 Krawczyk steps, 200
//!   contraction steps, levels 0..=6); there is no recursion.
//! - **Cancellation (#2136, #2371).** `gn_certify_span` polls the kernel's cancel token
//!   once per level, per piece's end station and per piece's bound; a cancel passes through
//!   `SM_GN_TRY` and the unsolved-piece fallback as `Cancelled` instead of proving
//!   nothing. A token that never fires changes nothing.
//! - **The uv curves' term (#2690).** Not the fork's: a face's point on a foot's uv curve
//!   is off the true foot by at most `|S_u| |du| + |S_v| |dv|`, each uv error `du`, `dv`
//!   bounded by its own parameter's fourth derivative, where the fork takes
//!   `(|S_u| + |S_v|) max(|du|, |dv|)`. The fork's form is never smaller, and on a face
//!   whose `u` (say) moves little while `|S_u|` is large it pairs that large partial with
//!   the other parameter's error and refines spans the tolerance does not need.
//!   `gn_certify_span` certifies a level by this bound but steers its levels (the least
//!   value it stops on, the halving it goes on with) by the fork's form, so a span takes no
//!   more levels or pieces than the fork's, and is certified at the same level or sooner;
//!   steered by the tighter bound, a round on a round took more spans and solves.
//! - **Indices.** Knots and nets are read by checked indices (a knot or a coordinate
//!   outside reads NaN, which shows in the result); the fork's indices stay inside.

use super::affine::{Aff, NOISE, S_NOISE};
use super::definition_box::{definition_partials_over_box, GnDefFace, GnDefs};
use super::form::{binom, gn_form_normalized, GnSectionFormT};
use super::interval::{cmax, cmin, round_down, round_up, Ivl};
use super::series::{GnScalar, Ser, Trunc, VSer, MAX_E, MAX_K, MAX_S, TAYLOR, V3};
use super::station::{
    ders_basis, gn_evaluate, gn_evaluate_vn, gn_jacobian, gn_solve_station_counted, gn_top_jets,
    knot, smax, solve4, span, GnStation, GnSystem, Nurbs, Partials,
};
use super::tay::{
    gn_half_angle, gn_invert, gn_target, tay_sqrt, GnAreaCells, Ring, Tay, TayScalar, TAY,
};
use crate::error::{Result, SolidError};
use crate::kernel::Kernel;
use crate::nurbs::b_spline_curve::BSplineCurve;
use crate::surface::surface::Surface;

const OP: &str = "blend::fillet_section::certificate";
const OP_CERTIFY_SPAN: &str = "blend::fillet_section::certificate::gn_certify_span";

/// `SM_GN_TRY(a)`: a failed evaluation proves nothing, the enclosing function returns
/// `Ok($none)` (the fork's `return SM_SUCCESS` with its outputs as they stand).
macro_rules! gn_try {
    ($e:expr, $none:expr) => {
        match $e {
            Ok(v) => v,
            // a cancel (#2136) is never a failed evaluation: it ends the operation
            Err(e @ SolidError::Cancelled { .. }) => return Err(e),
            Err(_) => return Ok($none),
        }
    };
}

/// `sm_NurbsCurve`: a NURBS curve's data (the leader), as [`Nurbs`].
#[derive(Clone, Debug, PartialEq)]
pub struct NurbsCurve {
    /// The degree.
    pub p: usize,
    /// The number of control points.
    pub n: usize,
    /// The full knot vector.
    pub u: Vec<f64>,
    /// The homogeneous control points `(x w, y w, z w, w)`.
    pub pw: Vec<f64>,
}

impl NurbsCurve {
    /// `sm_NurbsCurveOf`: the data of a B-spline curve.
    ///
    /// # Errors
    ///
    /// `SM_ERR` when the curve has no NURB, its degree, knots and control points do not
    /// agree, or a weight is not positive.
    pub fn of(c: &BSplineCurve) -> Result<Self> {
        let (knots, mult) = c.knots(None)?;
        let (cps, ws) = c.control_polygon()?;
        let p = c.degree()?;
        let n = cps.len();
        let mut u = Vec::new();
        for (&k, &m) in knots.iter().zip(&mult) {
            for _ in 0..m {
                u.push(k);
            }
        }
        if p < 1 || u.len() != n + p + 1 {
            // "SmGnSectionFSG: the centre line's NURBS data do not agree"
            return Err(SolidError::sm_err(OP));
        }
        let rational = ws.len() == cps.len();
        let mut pw = Vec::new();
        for (i, cp) in cps.iter().enumerate() {
            let w = if rational {
                ws.get(i).copied().unwrap_or(f64::NAN)
            } else {
                1.0
            };
            let positive = w > 0.0;
            if !positive {
                return Err(SolidError::sm_err(OP));
            }
            pw.extend_from_slice(&[cp.x * w, cp.y * w, cp.z * w, w]);
        }
        Ok(Self { p, n, u, pw })
    }
}

/// `crPw[i]` (NaN outside).
#[inline]
fn pw_at(pw: &[f64], i: usize) -> f64 {
    pw.get(i).copied().unwrap_or(f64::NAN)
}

/// `rM[i][j]` of a nested vector (0 outside; the fork's indices stay inside).
#[inline]
fn at2<T: Copy>(m: &[Vec<T>], i: usize, j: usize, zero: T) -> T {
    m.get(i).and_then(|r| r.get(j)).copied().unwrap_or(zero)
}

/// `v[i]` of a vector (0 outside).
#[inline]
fn at1<T: Copy>(v: &[T], i: usize, zero: T) -> T {
    v.get(i).copied().unwrap_or(zero)
}

const IZ: Ivl = Ivl::new(0.0, 0.0);
const VZ: V3<Ivl> = V3 {
    x: IZ,
    y: IZ,
    z: IZ,
};

/// `sm_CurveDerivsOver`: the derivatives `D[0..=ld]` of the curve's polynomial piece
/// `piece` (the span of the knots it starts at; `None`: each piece `s` covers) for every
/// parameter in `s`, enclosed (hulled over the pieces).
///
/// Not upstream (#2690): on each piece the quotient is taken about the curve's point at the
/// interval's middle, as [`nurbs_partials_over_piece`] does, so the enclosures do not
/// widen with the curve's distance from the origin.
///
/// # Errors
///
/// `SM_ERR` when no piece covers `s`.
pub fn curve_derivs_over(
    cn: &NurbsCurve,
    s: Ivl,
    ld: usize,
    piece: Option<usize>,
) -> Result<Vec<V3<Ivl>>> {
    let k0 = piece.unwrap_or_else(|| span(&cn.u, cn.p, cn.n, s.lo));
    let k1 = piece.unwrap_or_else(|| span(&cn.u, cn.p, cn.n, s.hi));
    let fixed = piece.is_some();
    let p = cn.p;
    let mut first = true;
    let mut rd = vec![VZ; ld + 1];
    for k in k0..=k1 {
        let ki = k as isize;
        let has_span = knot(&cn.u, ki) < knot(&cn.u, ki + 1);
        if !has_span {
            continue; // a knot repeated: no span
        }
        let a = if k == p || fixed {
            s.lo
        } else {
            cmax(s.lo, knot(&cn.u, ki))
        };
        let b = if k + 1 == cn.n || fixed {
            s.hi
        } else {
            cmin(s.hi, knot(&cn.u, ki + 1))
        };
        // an interval that ends on a knot is the span it lies in, not also the next one's
        // end point
        if a > b || (a == b && s.lo < s.hi) {
            continue;
        }
        let cv = piece_derivs_centred(cn, k, a, b, ld)?;
        for (r, c) in rd.iter_mut().zip(&cv) {
            if first {
                *r = *c;
            } else {
                r.x = r.x.hull(c.x);
                r.y = r.y.hull(c.y);
                r.z = r.z.hull(c.z);
            }
        }
        first = false;
    }
    if first {
        return Err(SolidError::sm_err(OP));
    }
    Ok(rd)
}

/// The derivatives `D[0..=ld]` of the curve's polynomial piece k over `[a, b]`, enclosed
/// directly: the homogeneous curve's Taylor polynomials at the middle m (it is of degree p
/// there), from exact point derivatives (rounding only), then the quotient (the NURBS
/// book's A4.2) in intervals, about the point at m (as [`nurbs_partials_direct`]). A point
/// `a = b` is enclosed to rounding.
fn piece_derivs_direct(cn: &NurbsCurve, k: usize, a: f64, b: f64, ld: usize) -> Vec<V3<Ivl>> {
    let p = cn.p;
    // the homogeneous curve is a polynomial of degree p on the span: its derivatives
    // over [a, b] are its Taylor polynomials at the middle m, from exact point
    // derivatives (rounding only)
    let m = 0.5 * (a + b);
    let nb = ders_basis(&cn.u, k, p, Ivl::from(m), p);
    let mut am = vec![VZ; p + 1];
    let mut wm = vec![IZ; p + 1];
    for d in 0..=p {
        let (mut x, mut y, mut z, mut w) = (IZ, IZ, IZ, IZ);
        for j in 0..=p {
            // k >= p (sm_Span)
            let li = 4 * (k + j - p);
            let nj = at2(&nb, d, j, IZ);
            x = x + nj * Ivl::from(pw_at(&cn.pw, li));
            y = y + nj * Ivl::from(pw_at(&cn.pw, li + 1));
            z = z + nj * Ivl::from(pw_at(&cn.pw, li + 2));
            w = w + nj * Ivl::from(pw_at(&cn.pw, li + 3));
        }
        if let (Some(ad), Some(wd)) = (am.get_mut(d), wm.get_mut(d)) {
            *ad = V3 { x, y, z };
            *wd = w;
        }
    }
    // the quotient about c, the point at the middle (rounded; any c is exact)
    let c = {
        let (a0, w0) = (at1(&am, 0, VZ), at1(&wm, 0, IZ).mid());
        let c = [a0.x.mid() / w0, a0.y.mid() / w0, a0.z.mid() / w0];
        if c.iter().all(|x| x.is_finite()) {
            c
        } else {
            [0.0; 3]
        }
    };
    let [cx, cy, cz] = c.map(Ivl::from);
    for (ad, &wd) in am.iter_mut().zip(&wm) {
        ad.x = ad.x - cx * wd;
        ad.y = ad.y - cy * wd;
        ad.z = ad.z - cz * wd;
    }
    let dt = Ivl::new(a, b) - Ivl::from(m);
    let mut av = vec![VZ; ld + 1];
    let mut wv = vec![IZ; ld + 1];
    for d in 0..=ld {
        let (mut x, mut y, mut z, mut w, mut pwr) = (IZ, IZ, IZ, IZ, Ivl::from(1.0));
        let mut f = 1.0_f64;
        for j in d..=p {
            if j > d {
                pwr = pwr * dt;
                f *= (j - d) as f64;
            }
            let c = pwr * (Ivl::from(1.0) / Ivl::from(f));
            let aj = at1(&am, j, VZ);
            x = x + aj.x * c;
            y = y + aj.y * c;
            z = z + aj.z * c;
            w = w + at1(&wm, j, IZ) * c;
        }
        if let (Some(ad), Some(wd)) = (av.get_mut(d), wv.get_mut(d)) {
            *ad = V3 { x, y, z };
            *wd = w;
        }
    }
    let w0 = at1(&wv, 0, IZ);
    let mut cv = vec![VZ; ld + 1];
    for d in 0..=ld {
        let mut v = at1(&av, d, VZ);
        let mut c = 1.0_f64;
        for i in 1..=d {
            c = c * (d - i + 1) as f64 / i as f64;
            let wi = at1(&wv, i, IZ);
            let cdi = at1(&cv, d - i, VZ);
            v.x = v.x - Ivl::from(c) * wi * cdi.x;
            v.y = v.y - Ivl::from(c) * wi * cdi.y;
            v.z = v.z - Ivl::from(c) * wi * cdi.z;
        }
        if let Some(cd) = cv.get_mut(d) {
            *cd = V3 {
                x: v.x / w0,
                y: v.y / w0,
                z: v.z / w0,
            };
        }
    }
    if let Some(c0) = cv.first_mut() {
        c0.x = c0.x + cx;
        c0.y = c0.y + cy;
        c0.z = c0.z + cz;
    }
    cv
}

/// [`piece_derivs_direct`] in the centred (mean value) form: `D_d(m) + D_(d+1)([a, b])
/// (s - m)`, the derivatives at the middle enclosed on that point and the next one's over
/// the interval, within the direct enclosure. The direct quotient widens with the interval
/// by the numerator's size times w's spread, which the derivatives' own ranges need not
/// (a line's tangent is constant, its rational quotient's is not): the centred form's
/// excess over the range is of the interval's square.
///
/// # Errors
///
/// `SM_ERR` where the two enclosures do not meet (they hold the same derivatives: only a
/// curve that is not finite there).
fn piece_derivs_centred(
    cn: &NurbsCurve,
    k: usize,
    a: f64,
    b: f64,
    ld: usize,
) -> Result<Vec<V3<Ivl>>> {
    if a >= b {
        return Ok(piece_derivs_direct(cn, k, a, b, ld));
    }
    let wide = piece_derivs_direct(cn, k, a, b, ld + 1);
    let m = 0.5 * (a + b);
    let mid = piece_derivs_direct(cn, k, m, m, ld);
    let dt = Ivl::new(a, b) - Ivl::from(m);
    let mut out = vec![VZ; ld + 1];
    for (d, o) in out.iter_mut().enumerate() {
        let (cm, cn1, cw) = (at1(&mid, d, VZ), at1(&wide, d + 1, VZ), at1(&wide, d, VZ));
        let comps = [
            (cm.x, cn1.x, cw.x),
            (cm.y, cn1.y, cw.y),
            (cm.z, cn1.z, cw.z),
        ];
        let mut got = [IZ; 3];
        for (g, (m, n1, w)) in got.iter_mut().zip(comps) {
            *g = centred_within(m + n1 * dt, w)?;
        }
        let [x, y, z] = got;
        *o = V3 { x, y, z };
    }
    Ok(out)
}

/// The centred form `c` within the direct enclosure `w` (both hold the value: they meet).
///
/// # Errors
///
/// `SM_ERR` where they do not meet.
fn centred_within(c: Ivl, w: Ivl) -> Result<Ivl> {
    let x = Ivl::new(cmax(c.lo, w.lo), cmin(c.hi, w.hi));
    let meet = x.lo <= x.hi;
    if !meet {
        return Err(SolidError::sm_err(OP));
    }
    Ok(x)
}

/// `sm_NurbsPartialsOverPiece`: the partials `P[a][b]`, `a + b <= lk`, of one polynomial
/// piece (spans `ku`, `kv`) over `[ua, ub] x [va, vb]`, enclosed.
///
/// Not upstream (#2690): the fork takes them directly, the homogeneous surface's Taylor
/// polynomials at the box's middle (exact) and the rational quotient (the NURBS book's
/// A4.4) in intervals. Interval arithmetic does not see that the numerator `A` and the
/// weight w vary together, so each quotient widens by about the numerator's size times w's
/// relative spread over the box: `|P| |dw| / w` for the point, measured from the origin
/// however near the box the face lies, and `|P_ab| |dw| / w` for each partial, which its
/// range need not have (a cylinder's partial along its axis is constant; its quotient's
/// enclosure is not). The Jacobian of the blend's system over a feet's box then comes out
/// thousands of times wider than its range, Krawczyk's test fails on spans far shorter
/// than it needs, and the certificate subdivides to make up for it. Here:
///
/// - the quotient is taken about the point at the box's middle, c: `P - c = (A - c w) /
///   w`, an identity for any c, with `A - c w` the homogeneous surface moved by -c, so the
///   point's excess is `|P - c| |dw| / w`, of the box's size, wherever the face lies;
/// - each partial is in the centred (mean value) form, `P_ab(m) + P_(a+1)b(box) du +
///   P_a(b+1)(box) dv`, the partials at the middle on that point and the next order's over
///   the box, within the direct enclosure: its excess over the range is of the box's
///   square.
///
/// # Errors
///
/// `SM_ERR` where the weights are not proved positive over the box ("a face's weights are
/// not positive over a box").
#[allow(clippy::too_many_arguments)]
pub fn nurbs_partials_over_piece(
    cn: &Nurbs,
    ku: usize,
    kv: usize,
    ua: f64,
    ub: f64,
    va: f64,
    vb: f64,
    lk: usize,
) -> Result<Partials<Ivl>> {
    let point = ua >= ub && va >= vb;
    // the centred form takes the next order over the box (a series holds MAX_K at most)
    if point || lk >= MAX_K {
        return nurbs_partials_direct(cn, ku, kv, ua, ub, va, vb, lk);
    }
    let wide = nurbs_partials_direct(cn, ku, kv, ua, ub, va, vb, lk + 1)?;
    let (um, vm) = (0.5 * (ua + ub), 0.5 * (va + vb));
    let mid = nurbs_partials_direct(cn, ku, kv, um, um, vm, vm, lk)?;
    let du = Ivl::new(ua, ub) - Ivl::from(um);
    let dv = Ivl::new(va, vb) - Ivl::from(vm);
    let mut out = Partials::default();
    for a in 0..=lk {
        for b in 0..=(lk - a) {
            let (cm, cu, cv, cw) = (
                mid.at(a, b),
                wide.at(a + 1, b),
                wide.at(a, b + 1),
                wide.at(a, b),
            );
            let comps = [
                (cm.x, cu.x, cv.x, cw.x),
                (cm.y, cu.y, cv.y, cw.y),
                (cm.z, cu.z, cv.z, cw.z),
            ];
            let mut got = [IZ; 3];
            for (g, (m, u, v, w)) in got.iter_mut().zip(comps) {
                *g = centred_within(m + u * du + v * dv, w)?;
            }
            let [x, y, z] = got;
            out.set(a, b, V3 { x, y, z });
        }
    }
    Ok(out)
}

/// [`nurbs_partials_over_piece`]'s direct enclosure: the homogeneous surface's Taylor
/// polynomials at the box's middle, from exact point derivatives (it is of degree p, q
/// there), then the rational quotient's (the NURBS book's A4.4) in intervals, about the
/// point at the middle. A point box is enclosed to rounding.
///
/// # Errors
///
/// As [`nurbs_partials_over_piece`].
#[allow(clippy::too_many_arguments)]
fn nurbs_partials_direct(
    cn: &Nurbs,
    ku: usize,
    kv: usize,
    ua: f64,
    ub: f64,
    va: f64,
    vb: f64,
    lk: usize,
) -> Result<Partials<Ivl>> {
    let um = 0.5 * (ua + ub);
    let vm = 0.5 * (va + vb);
    let nu = ders_basis(&cn.u, ku, cn.p, Ivl::from(um), cn.p);
    let nv = ders_basis(&cn.v, kv, cn.q, Ivl::from(vm), cn.q);
    let (la, lb) = (cn.p, cn.q);
    // homogeneous point derivatives Hm[a][b] (x, y, z, w) at the middle
    let mut hm = vec![vec![[IZ; 4]; lb + 1]; la + 1];
    for (a, row) in hm.iter_mut().enumerate() {
        for (b, hab) in row.iter_mut().enumerate() {
            let mut h = [IZ; 4];
            for i in 0..=cn.p {
                for j in 0..=cn.q {
                    // ku >= p and kv >= q (sm_Span)
                    let li = 4 * ((ku + i - cn.p) * cn.nv + (kv + j - cn.q));
                    let c = at2(&nu, a, i, IZ) * at2(&nv, b, j, IZ);
                    for (e, he) in h.iter_mut().enumerate() {
                        *he = *he + c * Ivl::from(pw_at(&cn.pw, li + e));
                    }
                }
            }
            *hab = h;
        }
    }
    // the quotient about c, the point at the middle (rounded; any c is exact)
    let c = {
        let [x, y, z, w] = at2(&hm, 0, 0, [IZ; 4]);
        [x.mid() / w.mid(), y.mid() / w.mid(), z.mid() / w.mid()]
    };
    let c_finite = c.iter().all(|x| x.is_finite());
    let c = if c_finite { c } else { [0.0; 3] };
    for h in hm.iter_mut().flatten() {
        let w = h[3];
        for (he, &ce) in h.iter_mut().zip(&c) {
            *he = *he - Ivl::from(ce) * w;
        }
    }
    let du = Ivl::new(ua, ub) - Ivl::from(um);
    let dv = Ivl::new(va, vb) - Ivl::from(vm);
    let mut pu = vec![IZ; la + 1];
    let mut pv = vec![IZ; lb + 1];
    if let (Some(u0), Some(v0)) = (pu.first_mut(), pv.first_mut()) {
        *u0 = Ivl::from(1.0);
        *v0 = Ivl::from(1.0);
    }
    for i in 1..=la {
        let x = at1(&pu, i - 1, IZ) * du * (Ivl::from(1.0) / Ivl::from(i as f64));
        if let Some(pi) = pu.get_mut(i) {
            *pi = x;
        }
    }
    for j in 1..=lb {
        let x = at1(&pv, j - 1, IZ) * dv * (Ivl::from(1.0) / Ivl::from(j as f64));
        if let Some(pj) = pv.get_mut(j) {
            *pj = x;
        }
    }
    let mut av = vec![vec![VZ; lk + 1]; lk + 1];
    let mut wv = vec![vec![IZ; lk + 1]; lk + 1];
    for a in 0..=lk {
        for b in 0..=(lk - a) {
            let mut h = [IZ; 4];
            for i in a..=la {
                for j in b..=lb {
                    let c = at1(&pu, i - a, IZ) * at1(&pv, j - b, IZ);
                    let hij = at2(&hm, i, j, [IZ; 4]);
                    for (he, &x) in h.iter_mut().zip(&hij) {
                        *he = *he + x * c;
                    }
                }
            }
            if let (Some(ab), Some(wab)) = (
                av.get_mut(a).and_then(|r| r.get_mut(b)),
                wv.get_mut(a).and_then(|r| r.get_mut(b)),
            ) {
                let [x, y, z, w] = h;
                *ab = V3 { x, y, z };
                *wab = w;
            }
        }
    }
    let w00 = at2(&wv, 0, 0, IZ);
    let positive = w00.lo > 0.0;
    if !positive {
        // "SmGnSectionFSG: a face's weights are not positive over a box"
        return Err(SolidError::sm_err(OP));
    }
    let mut out = Partials::default();
    for a in 0..=lk {
        for b in 0..=(lk - a) {
            let mut v3 = at2(&av, a, b, VZ);
            let mut cka = 1.0_f64;
            for i in 0..=a {
                if i > 0 {
                    cka = cka * (a - i + 1) as f64 / i as f64;
                }
                let mut cbj = 1.0_f64;
                for j in 0..=b {
                    if j > 0 {
                        cbj = cbj * (b - j + 1) as f64 / j as f64;
                    }
                    if i == 0 && j == 0 {
                        continue;
                    }
                    let c = Ivl::from(cka * cbj) * at2(&wv, i, j, IZ);
                    let pij = out.at(a - i, b - j);
                    v3.x = v3.x - c * pij.x;
                    v3.y = v3.y - c * pij.y;
                    v3.z = v3.z - c * pij.z;
                }
            }
            out.set(
                a,
                b,
                V3 {
                    x: v3.x / w00,
                    y: v3.y / w00,
                    z: v3.z / w00,
                },
            );
        }
    }
    let [cx, cy, cz] = c;
    let p00 = out.at(0, 0);
    out.set(
        0,
        0,
        V3 {
            x: p00.x + Ivl::from(cx),
            y: p00.y + Ivl::from(cy),
            z: p00.z + Ivl::from(cz),
        },
    );
    Ok(out)
}

/// `sm_LeaderSeriesOver`: the leader's point and unit tangent over `s` as series in sigma
/// to `tr.s`, enclosed (the span lies in one piece of the leader: the stations hold its
/// knots). `piece`: the leader's piece whose polynomial is taken over `s` (`lPiece`,
/// c5b6bb7; `None`, the fork's -1: the one the middle of `s` lies in).
///
/// # Errors
///
/// Those of [`curve_derivs_over`]; `SM_ERR` where the tangent is not proved non-zero.
pub fn leader_series_over(
    leader: &NurbsCurve,
    s: Ivl,
    tr: Trunc,
    piece: Option<usize>,
) -> Result<(VSer<Ivl>, VSer<Ivl>)> {
    let sigma = Trunc::new(tr.s(), 0)?;
    let piece = piece.unwrap_or_else(|| span(&leader.u, leader.p, leader.n, s.mid()));
    let d = curve_derivs_over(leader, s, tr.s() + 1, Some(piece))?;
    let mut l = VSer::constant(VZ);
    let mut d1 = l;
    let mut fact = 1.0_f64;
    for j in 0..=tr.s() {
        if j > 0 {
            fact *= j as f64;
        }
        let f = Ivl::from(1.0) / Ivl::from(fact);
        let dj = at1(&d, j, VZ);
        let dj1 = at1(&d, j + 1, VZ);
        let cl = [dj.x, dj.y, dj.z];
        let cd = [dj1.x, dj1.y, dj1.z];
        for ((la, da), (&x, &y)) in l.a.iter_mut().zip(d1.a.iter_mut()).zip(cl.iter().zip(&cd)) {
            la.set(j, 0, x * f);
            da.set(j, 0, y * f);
        }
    }
    let inv = d1.dot(&d1, sigma).inv_sqrt(sigma)?;
    let t = d1.times(&inv, sigma);
    Ok((l, t))
}

/// [`leader_series_over`] over an `s` that may cross the leader's knots (not upstream): `s`
/// cut at the knots strictly inside it, each part's series on its own piece, their hull
/// coefficient by coefficient. Each coefficient encloses its derivative at every point of
/// its part (at a knot, from either side), and the feet's series use them point by point
/// (Krawczyk's test is over the feet, for each s), so the hull is an enclosure whatever
/// the leader's continuity: where it has a corner in its parameter the two sides' values
/// widen the hull, and the bounds with it. Only a leader that is not continuous at a
/// crossed knot, or an `s` past its knots, is refused.
///
/// # Errors
///
/// [`SolidError::Unsupported`] for such a knot or `s`; those of [`leader_series_over`].
pub fn leader_series_piecewise(
    leader: &NurbsCurve,
    s: Ivl,
    tr: Trunc,
) -> Result<(VSer<Ivl>, VSer<Ivl>)> {
    let hull = |a: &mut VSer<Ivl>, b: &VSer<Ivl>| {
        for (sa, sb) in a.a.iter_mut().zip(&b.a) {
            for (ra, rb) in sa.c.iter_mut().zip(&sb.c) {
                for (x, y) in ra.iter_mut().zip(rb) {
                    *x = x.hull(*y);
                }
            }
        }
    };
    let mut out: Option<(VSer<Ivl>, VSer<Ivl>)> = None;
    for (a, b, k) in knot_cuts(&leader.u, leader.p, leader.n, s.lo, s.hi, 0)? {
        let (l, t) = leader_series_over(leader, Ivl::new(a, b), tr, Some(k))?;
        out = Some(match out {
            None => (l, t),
            Some((mut hl, mut ht)) => {
                hull(&mut hl, &l);
                hull(&mut ht, &t);
                (hl, ht)
            }
        });
    }
    out.ok_or_else(|| SolidError::sm_err(OP))
}

/// `sm_GnLaw` (0d43b4c): a ball whose radius varies, its law r(s) over the leader's
/// parameter (the x of a B-spline curve), with the curve's NURBS data.
#[derive(Clone, Debug)]
pub struct GnLaw<'a> {
    /// `pCurve`.
    pub curve: &'a BSplineCurve,
    /// `sN`.
    pub n: NurbsCurve,
}

impl<'a> GnLaw<'a> {
    /// The law of `curve` (`sLaw.pCurve = pLaw; sm_NurbsCurveOf(*pLaw, sLaw.sN)`).
    ///
    /// # Errors
    ///
    /// Those of [`NurbsCurve::of`].
    pub fn of(curve: &'a BSplineCurve) -> Result<Self> {
        Ok(Self {
            curve,
            n: NurbsCurve::of(curve)?,
        })
    }
}

/// `sm_LawSeries` (0d43b4c): r at `s` to order `order` in sigma (its derivatives past a
/// polynomial law's degree zero), `None` without a law; `before`: the law's piece before
/// `s` at one of its knots.
///
/// # Errors
///
/// Those of the law's evaluation; `SM_ERR` where the radius is not positive ("the radius
/// law is not positive").
pub fn law_series(
    law: Option<&GnLaw<'_>>,
    s: f64,
    order: usize,
    before: bool,
) -> Result<Option<Ser<f64>>> {
    let Some(law) = law else {
        return Ok(None);
    };
    let mut r = Ser::constant(0.0_f64);
    let mut derivs = order;
    if !law.curve.is_rational()? {
        derivs = derivs.min(law.curve.degree()?);
    }
    if derivs > MAX_S {
        return Err(SolidError::sm_err(OP));
    }
    let d = law.curve.evaluate(s, derivs, !before, false)?;
    let mut fact = 1.0_f64;
    for i in 0..=derivs {
        if i > 0 {
            fact *= i as f64;
        }
        let di = d.get(i).ok_or_else(|| SolidError::sm_err(OP))?;
        r.set(i, 0, di.x / fact);
    }
    let positive = r.at(0, 0) > 0.0;
    if !positive {
        // "SmGnSectionFSG: the radius law is not positive"
        return Err(SolidError::sm_err(OP));
    }
    Ok(Some(r))
}

/// `sm_LawSeriesOver` (0d43b4c): r over the span `s` (in one piece of the law: the
/// stations hold its knots) to order `order`, enclosed.
///
/// # Errors
///
/// Those of [`curve_derivs_over`]; `SM_ERR` where r is not proved positive.
pub fn law_series_over(law: &GnLaw<'_>, s: Ivl, order: usize) -> Result<Ser<Ivl>> {
    let mut r = Ser::constant(IZ);
    let piece = span(&law.n.u, law.n.p, law.n.n, s.mid());
    let d = curve_derivs_over(&law.n, s, order, Some(piece))?;
    let mut fact = 1.0_f64;
    for i in 0..=order {
        if i > 0 {
            fact *= i as f64;
        }
        r.set(i, 0, at1(&d, i, VZ).x * (Ivl::from(1.0) / Ivl::from(fact)));
    }
    if !r.at(0, 0).is_positive() {
        return Err(SolidError::sm_err(OP));
    }
    Ok(r)
}

/// `sm_Inverse4`: the inverse of a 4 x 4 matrix (the preconditioner; any approximation
/// will do), column by column with [`solve4`].
///
/// # Errors
///
/// Those of [`solve4`] (a singular matrix).
pub fn inverse4(a: &[[f64; 4]; 4]) -> Result<[[f64; 4]; 4]> {
    let mut y = [[0.0_f64; 4]; 4];
    for j in 0..4 {
        let mut e = [0.0_f64; 4];
        if let Some(ej) = e.get_mut(j) {
            *ej = 1.0;
        }
        let col = solve4(a, &e)?;
        for (row, &c) in y.iter_mut().zip(&col) {
            if let Some(x) = row.get_mut(j) {
                *x = c;
            }
        }
    }
    Ok(y)
}

/// The proved enclosure of a span (`sm_GnEnclose`'s outputs).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Enclosure {
    /// The feet's series over the span (`rF`).
    pub f: [VSer<Ivl>; 2],
    /// The feet's parameters' series (`pX`).
    pub x: [Ser<Ivl>; 4],
    /// Per face, `[|S_u|, |S_v|]` over the feet's box (the fork's `pDS` is their sum):
    /// `|S_u| |du| + |S_v| |dv|` bounds `|dS|` for a uv step `(du, dv)` in the box (#2690).
    pub ds: [[f64; 2]; 2],
    /// The same, least (the fork's `pDSLo`, their sum): no smaller anywhere in the box.
    pub ds_lo: [[f64; 2]; 2],
    /// The point V's series over the span (`pV`, 7102bc5).
    pub v: VSer<Ivl>,
}

fn ivl_at(x: &[Ivl; 4], i: usize) -> Ivl {
    x.get(i).copied().unwrap_or(IZ)
}

fn f_at(x: &[f64; 4], i: usize) -> f64 {
    x.get(i).copied().unwrap_or(f64::NAN)
}

/// The system over the partials `p` (`sSys` with `apP` at `asP`).
fn system<'a, T: GnScalar>(
    p: &'a [Partials<T>; 2],
    leader: &VSer<T>,
    tangent: &VSer<T>,
    r: [f64; 2],
    side: [f64; 2],
    law: Option<&'a Ser<T>>,
) -> GnSystem<'a, T> {
    let [p0, p1] = p;
    GnSystem {
        p: [p0, p1],
        leader: *leader,
        tangent: *tangent,
        r,
        side,
        hold: false,
        law,
    }
}

/// A face's partials over `[ua, ub] x [va, vb]`, enclosed (`sm_FacePartialsOver`,
/// c5b6bb7): a NURBS face's piece `(ku, kv)`, continued; a blend by its definition's on
/// the piece of its net along s that `(ku, kv)` is in
/// ([`super::definition_box::definition_partials_over_box`]).
///
/// # Errors
///
/// Those of [`nurbs_partials_over_piece`] and
/// [`super::definition_box::definition_partials_over_box`].
#[allow(clippy::too_many_arguments)]
pub fn face_partials_over(
    k: &Kernel,
    cn: &Nurbs,
    def: Option<&mut GnDefFace<'_>>,
    ku: usize,
    kv: usize,
    ua: f64,
    ub: f64,
    va: f64,
    vb: f64,
    lk: usize,
) -> Result<Partials<Ivl>> {
    if let Some(d) = def {
        return definition_partials_over_box(k, d, cn, ku, kv, ua, ub, va, vb, lk);
    }
    nurbs_partials_over_piece(cn, ku, kv, ua, ub, va, vb, lk)
}

/// The parts of `[lo, hi]` between the knots strictly inside it, each with the span whose
/// closed interval holds it: a knot crossed where the direction is continuous to less than
/// `continuity` (degree minus multiplicity) is refused, and so is an interval past the
/// knots' range (the face or curve is not defined there). A point `lo = hi` on an inner
/// knot lies in both spans' closed intervals: it is returned twice, once with each, so a
/// hull over the parts holds both one-sided values there (#2570).
pub(crate) fn knot_cuts(
    knots: &[f64],
    p: usize,
    n: usize,
    lo: f64,
    hi: f64,
    continuity: usize,
) -> Result<Vec<(f64, f64, usize)>> {
    let refused = |reason| SolidError::Unsupported {
        operation: OP,
        reason,
    };
    let (first, last) = (knot(knots, p as isize), knot(knots, n as isize));
    if !(lo >= first && hi <= last && lo <= hi) {
        return Err(refused("a box past the face's knots"));
    }
    let after = span(knots, p, n, lo);
    if lo == hi && after > p && knot(knots, after as isize) == lo {
        // the span that ends at the knot: the last one before it that is not empty
        let mut before = after - 1;
        while before > p && knot(knots, before as isize) >= knot(knots, before as isize + 1) {
            before -= 1;
        }
        return Ok(vec![(lo, hi, before), (lo, hi, after)]);
    }
    let mut at = vec![lo];
    let mut prev = f64::NAN;
    for &x in knots {
        if x > lo && x < hi && x != prev {
            let m = knots.iter().filter(|&&y| y == x).count();
            if p.checked_sub(m).is_none_or(|c| c < continuity) {
                return Err(refused(
                    "a box crosses a knot line where the face or curve is not continuous enough",
                ));
            }
            at.push(x);
        }
        prev = x;
    }
    at.push(hi);
    Ok(at
        .iter()
        .zip(at.iter().skip(1))
        .map(|(&a, &b)| (a, b, span(knots, p, n, 0.5 * (a + b))))
        .collect())
}

/// Whether `[lo, hi]` lies on more than one piece of `knots` (it crosses an inner knot, or
/// is a point on one), or past their range.
pub(crate) fn crosses_a_knot(knots: &[f64], p: usize, n: usize, lo: f64, hi: f64) -> bool {
    knot_cuts(knots, p, n, lo, hi, 0).map_or(true, |c| c.len() > 1)
}

/// The continuity a face must have at a knot line a feet's box crosses
/// ([`nurbs_partials_over_box`]): Krawczyk's mean-value form needs the feet's system (the
/// faces' points and first partials) continuous, with its Jacobian (their second partials)
/// enclosed wherever it exists; at C^2 the system is C^1 across the line.
pub const FEET_BOX_CONTINUITY: usize = 2;

/// A NURBS face's partials `P[a][b]`, `a + b <= lk`, over `[ua, ub] x [va, vb]` wherever
/// it lies on the face's knots, enclosed (not upstream): the box cut at the knot lines
/// strictly inside it, each part enclosed on its own polynomial piece
/// ([`nurbs_partials_over_piece`] over the part, its closed span), and the enclosures'
/// hull. Every point of the box lies in some closed part, where that piece is the face, so
/// the hull holds, at every point, each partial the face has there (on a knot line, from
/// either side): an enclosure point by point, whatever the continuity. A use that also
/// integrates over the box (a mean-value form) needs the face continuous to some order
/// across the lines it crosses: a crossed line where the face's continuity (degree minus
/// the knot's multiplicity) is below `continuity` is refused, as is a box past the knots'
/// range. No station invariant is assumed.
///
/// # Errors
///
/// [`SolidError::Unsupported`] for a crossed knot line of too low a continuity or a box
/// past the knots; those of [`nurbs_partials_over_piece`].
#[allow(clippy::too_many_arguments)]
pub fn nurbs_partials_over_box(
    cn: &Nurbs,
    ua: f64,
    ub: f64,
    va: f64,
    vb: f64,
    lk: usize,
    continuity: usize,
) -> Result<Partials<Ivl>> {
    let us = knot_cuts(&cn.u, cn.p, cn.nu, ua, ub, continuity)?;
    let vs = knot_cuts(&cn.v, cn.q, cn.nv, va, vb, continuity)?;
    let mut out: Option<Partials<Ivl>> = None;
    for &(a0, a1, ku) in &us {
        for &(b0, b1, kv) in &vs {
            let p = nurbs_partials_over_piece(cn, ku, kv, a0, a1, b0, b1, lk)?;
            out = Some(match out {
                None => p,
                Some(h) => partials_hull(h, &p, lk),
            });
        }
    }
    out.ok_or_else(|| SolidError::sm_err(OP))
}

/// The hull of two enclosures of partials, `a + b <= lk`.
pub(crate) fn partials_hull(mut h: Partials<Ivl>, p: &Partials<Ivl>, lk: usize) -> Partials<Ivl> {
    for a in 0..=lk {
        for b in 0..=(lk - a) {
            let (x, y) = (h.at(a, b), p.at(a, b));
            let v = V3 {
                x: x.x.hull(y.x),
                y: x.y.hull(y.y),
                z: x.z.hull(y.z),
            };
            h.set(a, b, v);
        }
    }
    h
}

/// A span the feet are enclosed over (`sm_GnSpanAt`, c5b6bb7): S, the feet's guesses at
/// its ends (their hull, widened, is the first box), each face's polynomial piece over it
/// and the leader's (`None`: the one the middle of S lies in).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpanAt {
    /// `sS`.
    pub s: Ivl,
    /// `axA`: the feet's guesses at S's start.
    pub xa: [f64; 4],
    /// `axB`: and at its end.
    pub xb: [f64; 4],
    /// `aku`: each face's piece in u.
    pub ku: [usize; 2],
    /// `akv`: and in v.
    pub kv: [usize; 2],
    /// `lLeader`.
    pub leader: Option<usize>,
    /// Not upstream: a NURBS face's partials over a box are enclosed piece by piece
    /// ([`nurbs_partials_over_box`]), so the feet may cross knot lines (`ku`, `kv` unused for
    /// them), and the leader's series over S on the piece `leader` names (the caller cut S
    /// at the leader's knots: a point S on a knot, the side it names) or, `None`, piece by
    /// piece ([`leader_series_piecewise`]); `false`: on the pieces `ku`, `kv` and `leader`,
    /// continued (upstream's, which relies on the station invariant).
    pub piecewise: bool,
}

/// `sm_GnEnclose`: the feet's series over `S = [a.s, b.s]` (a point when they are one
/// station), to order `order_s` in s and `n` in tau, enclosed; `None` where Krawczyk's test
/// fails on S (or an evaluation fails: nothing is proved). `defs`: the faces that are
/// blends by their definitions, enclosed over the boxes (c5b6bb7).
///
/// # Errors
///
/// `InvalidInput` when `order_s` or `n` exceeds the series' size.
#[allow(clippy::too_many_arguments)]
pub fn gn_enclose(
    k: &Kernel,
    faces_n: [&Nurbs; 2],
    leader_n: &NurbsCurve,
    a: &GnStation,
    b: &GnStation,
    r: [f64; 2],
    side: [f64; 2],
    n: usize,
    order_s: usize,
    law: Option<&GnLaw<'_>>,
    defs: Option<&mut GnDefs<'_>>,
) -> Result<Option<Enclosure>> {
    let ax = |c: usize| a.x.get(c).map_or(f64::NAN, |x| x.at(0, 0));
    let bx = |c: usize| b.x.get(c).map_or(f64::NAN, |x| x.at(0, 0));
    let mut at = SpanAt {
        s: Ivl::new(a.s, b.s),
        xa: [ax(0), ax(1), ax(2), ax(3)],
        xb: [bx(0), bx(1), bx(2), bx(3)],
        ku: [0; 2],
        kv: [0; 2],
        leader: None,
        piecewise: false,
    };
    // each face's polynomial piece over the span: the one its middle lies in. A station
    // stands where a foot crosses a knot line, so every foot over the span lies in that
    // piece, and the piece's polynomial, continued over the rest of a box, is the face
    // wherever a foot can be: the root Krawczyk's test proves unique in a box is the foot
    // (a definition's piece of its net is one between two of its stations, its polynomial
    // pieces the ones the faces under it have there)
    for (f, face) in faces_n.iter().enumerate() {
        let su = span(&face.u, face.p, face.nu, 0.5 * (ax(2 * f) + bx(2 * f)));
        let sv = span(
            &face.v,
            face.q,
            face.nv,
            0.5 * (ax(2 * f + 1) + bx(2 * f + 1)),
        );
        if let (Some(u), Some(v)) = (at.ku.get_mut(f), at.kv.get_mut(f)) {
            *u = su;
            *v = sv;
        }
    }
    gn_enclose_over(k, faces_n, leader_n, &at, r, side, n, order_s, law, defs)
}

/// `sm_GnEncloseOver` (c5b6bb7): the feet's series over `at.s`, to order `order_s` in s and
/// `n` in tau, enclosed; `None` where Krawczyk's test fails on S (or an evaluation fails:
/// nothing is proved).
///
/// Not upstream (#2690): the meet of two enclosures, the affine one
/// ([`gn_enclose_affine`]) and the fork's in plain intervals ([`gn_enclose_plain`]). Each
/// holds, for every s in S, the feet's unique root and its Taylor coefficients, so their
/// intersection, coefficient by coefficient, does too; it is proved where either is. The
/// affine forms cancel the system's first-order dependence on the feet and on s and win
/// by orders of magnitude where the ball's offset nearly collapses; plain intervals win
/// where a face's partials have no form (a definition's, or a box across a knot line),
/// and on long spans where the forms' second-order remainders outgrow the box: the meet
/// is never wider than the fork's.
///
/// # Errors
///
/// `InvalidInput` when `order_s` or `n` exceeds the series' size.
#[allow(clippy::too_many_arguments)]
pub fn gn_enclose_over(
    k: &Kernel,
    faces_n: [&Nurbs; 2],
    leader_n: &NurbsCurve,
    at: &SpanAt,
    r: [f64; 2],
    side: [f64; 2],
    n: usize,
    order_s: usize,
    law: Option<&GnLaw<'_>>,
    defs: Option<&mut GnDefs<'_>>,
) -> Result<Option<Enclosure>> {
    crate::profile::profile_scope!(k, "blend::certificate::gn_enclose_over");
    let mut defs = defs;
    let affine = gn_enclose_affine(
        k,
        faces_n,
        leader_n,
        at,
        r,
        side,
        n,
        order_s,
        law,
        defs.as_deref_mut(),
    )?;
    let plain = gn_enclose_plain(k, faces_n, leader_n, at, r, side, n, order_s, law, defs)?;
    Ok(match (affine, plain) {
        (Some(a), Some(p)) => Some(gn_try!(meet(&a, &p), None)),
        (a, p) => a.or(p),
    })
}

/// The meet of two enclosures of the same span, coefficient by coefficient.
///
/// # Errors
///
/// `SM_ERR` where two do not meet (both hold the true value: only a defect).
fn meet(a: &Enclosure, b: &Enclosure) -> Result<Enclosure> {
    let ivl = |x: Ivl, y: Ivl| -> Result<Ivl> {
        let m = Ivl::new(cmax(x.lo, y.lo), cmin(x.hi, y.hi));
        let meets = m.lo <= m.hi;
        if !meets {
            return Err(SolidError::sm_err(OP));
        }
        Ok(m)
    };
    let ser = |x: &Ser<Ivl>, y: &Ser<Ivl>| -> Result<Ser<Ivl>> {
        let mut out = *x;
        for (orow, yrow) in out.c.iter_mut().zip(&y.c) {
            for (o, &yv) in orow.iter_mut().zip(yrow) {
                *o = ivl(*o, yv)?;
            }
        }
        Ok(out)
    };
    let vser = |x: &VSer<Ivl>, y: &VSer<Ivl>| -> Result<VSer<Ivl>> {
        let [x0, x1, x2] = &x.a;
        let [y0, y1, y2] = &y.a;
        Ok(VSer {
            a: [ser(x0, y0)?, ser(x1, y1)?, ser(x2, y2)?],
        })
    };
    let [af0, af1] = &a.f;
    let [bf0, bf1] = &b.f;
    let [ax0, ax1, ax2, ax3] = &a.x;
    let [bx0, bx1, bx2, bx3] = &b.x;
    let mut ds = a.ds;
    let mut ds_lo = a.ds_lo;
    for ((d, l), (bd, bl)) in ds
        .iter_mut()
        .zip(ds_lo.iter_mut())
        .zip(b.ds.iter().zip(&b.ds_lo))
    {
        for ((dv, lv), (&bdv, &blv)) in d.iter_mut().zip(l.iter_mut()).zip(bd.iter().zip(bl)) {
            *dv = cmin(*dv, bdv);
            *lv = cmax(*lv, blv);
        }
    }
    Ok(Enclosure {
        f: [vser(af0, bf0)?, vser(af1, bf1)?],
        x: [
            ser(ax0, bx0)?,
            ser(ax1, bx1)?,
            ser(ax2, bx2)?,
            ser(ax3, bx3)?,
        ],
        ds,
        ds_lo,
        v: vser(&a.v, &b.v)?,
    })
}

/// The fork's `sm_GnEncloseOver` (c5b6bb7), in plain intervals: the feet over S in one
/// box about the guesses' middle, Krawczyk's operator on it, and the feet's Taylor
/// coefficients order by order (each `-Y R + C x`, `|x| <= |Y R| / (1 - q)`).
/// [`gn_enclose_over`] meets its enclosures with the affine ones (#2690).
///
/// # Errors
///
/// `InvalidInput` when `order_s` or `n` exceeds the series' size.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn gn_enclose_plain(
    k: &Kernel,
    faces_n: [&Nurbs; 2],
    leader_n: &NurbsCurve,
    at: &SpanAt,
    r: [f64; 2],
    side: [f64; 2],
    n: usize,
    order_s: usize,
    law: Option<&GnLaw<'_>>,
    defs: Option<&mut GnDefs<'_>>,
) -> Result<Option<Enclosure>> {
    let s = at.s;
    let tr = Trunc::new(order_s, n)?;
    // a ball whose radius varies: its law over the span (0d43b4c)
    let law_r = match law {
        Some(l) => Some(gn_try!(law_series_over(l, s, order_s), None)),
        None => None,
    };
    let law_r = law_r.as_ref();
    let (leader, tangent) = if at.piecewise && at.leader.is_none() {
        gn_try!(leader_series_piecewise(leader_n, s, tr), None)
    } else {
        gn_try!(leader_series_over(leader_n, s, tr, at.leader), None)
    };
    let ax = |c: usize| at.xa.get(c).copied().unwrap_or(f64::NAN);
    let bx = |c: usize| at.xb.get(c).copied().unwrap_or(f64::NAN);
    let [ku0, ku1] = at.ku;
    let [kv0, kv1] = at.kv;
    let mut defs = defs;
    let mut partials_over = |xl: [f64; 4], xh: [f64; 4], lk: usize| -> Result<[Partials<Ivl>; 2]> {
        let [n0, n1] = faces_n;
        let [l0, l1, l2, l3] = xl;
        let [h0, h1, h2, h3] = xh;
        let (d0, d1) = match defs.as_deref_mut() {
            Some(d) => {
                let [a0, a1] = &mut d.faces;
                (a0.as_mut(), a1.as_mut())
            }
            None => (None, None),
        };
        if at.piecewise && d0.is_none() && d1.is_none() {
            return Ok([
                nurbs_partials_over_box(n0, l0, h0, l1, h1, lk, FEET_BOX_CONTINUITY)?,
                nurbs_partials_over_box(n1, l2, h2, l3, h3, lk, FEET_BOX_CONTINUITY)?,
            ]);
        }
        Ok([
            face_partials_over(k, n0, d0, ku0, kv0, l0, h0, l1, h1, lk)?,
            face_partials_over(k, n1, d1, ku1, kv1, l2, h2, l3, h3, lk)?,
        ])
    };

    // the feet's box: Krawczyk's operator from the guesses' hull
    let mut xm = [0.0_f64; 4];
    let mut bx_ = [IZ; 4];
    for (c, (m, xc)) in xm.iter_mut().zip(bx_.iter_mut()).enumerate() {
        let (lo, hi) = (ax(c), bx(c));
        *m = 0.5 * (lo + hi);
        let w = 0.5 * (hi - lo).abs() + 1.0e-9 * (1.0 + m.abs());
        *xc = Ivl::new(*m - 2.0 * w, *m + 2.0 * w);
    }
    let mut y = [[0.0_f64; 4]; 4];
    let mut cm = [[IZ; 4]; 4];
    let mut q = f64::INFINITY;
    let mut proved = false;
    let point = Trunc::new(0, 0)?;
    // G at the box's centre (which the iterations keep) for every s in the span
    let sg = {
        let pc = gn_try!(partials_over(xm, xm, 2), None);
        let sxc = xm.map(|x| Ser::constant(Ivl::from(x)));
        let (sg, _) = gn_try!(
            gn_evaluate(&system(&pc, &leader, &tangent, r, side, law_r), &sxc, point),
            None
        );
        sg
    };
    let mut it = 0;
    while it < 8 && !proved {
        it += 1;
        // the Jacobian over the box
        let pb = gn_try!(partials_over(bx_.map(|x| x.lo), bx_.map(|x| x.hi), 2), None);
        let sxb = bx_.map(Ser::constant);
        let jac = gn_try!(
            gn_jacobian(&system(&pb, &leader, &tangent, r, side, law_r), &sxb),
            None
        );
        let jm = jac.map(|row| row.map(Ivl::mid));
        y = gn_try!(inverse4(&jm), None);
        q = 0.0;
        for (i, (yi, ci)) in y.iter().zip(cm.iter_mut()).enumerate() {
            let mut row = 0.0_f64;
            for (j, cij) in ci.iter_mut().enumerate() {
                let mut c = Ivl::from(if i == j { 1.0 } else { 0.0 });
                for (&yik, jk) in yi.iter().zip(&jac) {
                    c = c - Ivl::from(yik) * jk.get(j).copied().unwrap_or(IZ);
                }
                *cij = c;
                row = round_up(row + c.mag());
            }
            q = cmax(q, row);
        }
        let mut kk = [IZ; 4];
        let mut inside = q < 1.0;
        for (i, ((yi, ci), ki)) in y.iter().zip(&cm).zip(kk.iter_mut()).enumerate() {
            let mut k = Ivl::from(f_at(&xm, i));
            for (&yij, gj) in yi.iter().zip(&sg) {
                k = k - Ivl::from(yij) * gj.at(0, 0);
            }
            for (j, &cij) in ci.iter().enumerate() {
                k = k + cij * (ivl_at(&bx_, j) - Ivl::from(f_at(&xm, j)));
            }
            *ki = k;
            let xi = ivl_at(&bx_, i);
            if !(k.lo > xi.lo && k.hi < xi.hi) {
                inside = false;
            }
        }
        if inside {
            for (xi, ki) in bx_.iter_mut().zip(&kk) {
                *xi = Ivl::new(cmax(xi.lo, ki.lo), cmin(xi.hi, ki.hi));
            }
            proved = true;
        } else {
            let contracts = q < 1.0;
            if !contracts {
                return Ok(None); // the span is too long for its Jacobian to be one preconditioner's
            }
            // the next box from the contraction: K's radii are at most d + |C| w for a box
            // of radii w about the centre, d those of Y G, so the box of radii
            // w = (I - |C|)^-1 d holds its K; a half more, to clear rounding and the change
            // of C with the box
            let mut d = [0.0_f64; 4];
            for (di, yi) in d.iter_mut().zip(&y) {
                let mut v = IZ;
                for (&yij, gj) in yi.iter().zip(&sg) {
                    v = v + Ivl::from(yij) * gj.at(0, 0);
                }
                *di = (v - Ivl::from(v.mid())).mag() + v.mid().abs();
            }
            let mut w = d;
            for _ in 0..200 {
                let mut wn = [0.0_f64; 4];
                for ((wni, &di), ci) in wn.iter_mut().zip(&d).zip(&cm) {
                    *wni = di;
                    for (cij, &wj) in ci.iter().zip(&w) {
                        *wni += cij.mag() * wj;
                    }
                }
                w = wn;
            }
            for ((xi, &wi), &mi) in bx_.iter_mut().zip(&w).zip(&xm) {
                let rr = 1.5 * wi + 1.0e-12 * (1.0 + mi.abs());
                if !rr.is_finite() {
                    return Ok(None);
                }
                *xi = Ivl::new(mi - rr, mi + rr);
            }
        }
    }
    if !proved {
        return Ok(None);
    }

    // the feet's Taylor coefficients over the span, order by order: x = -Y R + C x,
    // |x| <= rho
    let lk = tr.s() + tr.e() + 1;
    let pk = gn_try!(
        partials_over(bx_.map(|x| x.lo), bx_.map(|x| x.hi), lk),
        None
    );
    let sysk = system(&pk, &leader, &tangent, r, side, law_r);
    let mut sx = bx_.map(Ser::constant);
    for tot in 1..=(tr.s() + tr.e()) {
        let (sg, _) = gn_try!(gn_evaluate(&sysk, &sx, tr), None);
        for i in 0..=cmin_usize(tr.s(), tot) {
            let k = tot - i;
            // the section's jets of order n are the face curve's (gn_top_jets, 3bc54a4)
            if k > tr.e() || (n >= 1 && k >= n) {
                continue;
            }
            let mut yr = [IZ; 4];
            let mut rho = 0.0_f64;
            for (yrr, yrow) in yr.iter_mut().zip(&y) {
                let mut v = IZ;
                for (&yrc, gc) in yrow.iter().zip(&sg) {
                    v = v - Ivl::from(yrc) * gc.at(i, k);
                }
                *yrr = v;
                rho = cmax(rho, v.mag());
            }
            // 1 - q rounded down: the quotient rounded up is then at least the exact one
            let gap = round_down(1.0 - q);
            if !gap.is_finite() || gap <= 0.0 {
                return Ok(None);
            }
            rho = round_up(rho / gap);
            if !rho.is_finite() {
                return Ok(None);
            }
            for ((sxr, &yrr), crow) in sx.iter_mut().zip(&yr).zip(&cm) {
                let mut v = yrr;
                for &crc in crow {
                    v = v + crc * Ivl::new(-rho, rho);
                }
                sxr.set(i, k, v);
            }
        }
    }
    let (_, mut f, v, nrm) = gn_try!(gn_evaluate_vn(&sysk, &sx, tr), None);
    gn_top_jets(n, tr, &nrm, &mut f);
    // the faces' first partials over the box: |S_u| |du| + |S_v| |dv| bounds |dS| for a uv
    // step (du, dv) (#2690: kept apart, where the fork sums them)
    let mut ds = [[0.0_f64; 2]; 2];
    let mut ds_lo = [[0.0_f64; 2]; 2];
    for ((dsf, dsl), pf) in ds.iter_mut().zip(ds_lo.iter_mut()).zip(&pk) {
        let su = pf.at(1, 0);
        let sv = pf.at(0, 1);
        let du = crate::geom::vector3d::sqrt(
            su.x.mag() * su.x.mag() + su.y.mag() * su.y.mag() + su.z.mag() * su.z.mag(),
        );
        let dv = crate::geom::vector3d::sqrt(
            sv.x.mag() * sv.x.mag() + sv.y.mag() * sv.y.mag() + sv.z.mag() * sv.z.mag(),
        );
        *dsf = [du, dv];
        let lu = crate::geom::vector3d::sqrt(
            su.x.mig() * su.x.mig() + su.y.mig() * su.y.mig() + su.z.mig() * su.z.mig(),
        );
        let lv = crate::geom::vector3d::sqrt(
            sv.x.mig() * sv.x.mig() + sv.y.mig() * sv.y.mig() + sv.z.mig() * sv.z.mig(),
        );
        *dsl = [lu, lv];
    }
    Ok(Some(Enclosure {
        f,
        x: sx,
        ds,
        ds_lo,
        v,
    }))
}

/// The feet's series over `at.s` in affine forms ([`gn_enclose_over`]'s first enclosure).
///
/// Not upstream (#2690): the fork encloses the feet over S in one box and evaluates the
/// system over it in plain intervals. The feet are a curve in s; a box treats their four
/// parameters as free, and interval arithmetic adds the widths of the system's terms,
/// which nearly cancel where the ball's offset nearly collapses (`S + r N`): the
/// Jacobian's enclosure came out hundreds to thousands of times its range, and each
/// Taylor coefficient, solved through the preconditioner from the ones before, widened by
/// that factor again order by order (the sixth by 10^5 to 10^6 on #2406's plate), so
/// Krawczyk's test failed on all but short spans and the certificate's bound, made of the
/// top coefficient, subdivided for what was its own excess. Here every quantity is a
/// first-order affine form ([`Aff`]) in the span's parameter, `s = m + rad e_4`, and the
/// feet's residuals off their predicted line, `x = x^(s) + E`, `x^(s)` the line through
/// the guesses at S's ends and E a box `|e_q| <= w_q`:
///
/// - the leader, its tangent and a radius law are their mean value forms in s
///   ([`series_affine`]), the faces' partials their mean value forms over the region the
///   feet range over ([`face_partials_affine`]);
/// - Krawczyk's operator is taken for each s about `x^(s)`, with the preconditioner
///   `Y(s) = Y0 - Y0 J1 Y0 (s - m) / rad`, the inverse of the Jacobian's middle to first
///   order in s, so `C = I - Y J` is of the span's square and the box's size:
///   `-Y G(x^(s)) + C E` inside E proves one root in `x^(s) + E` for every s;
/// - each Taylor coefficient is the preconditioned solution's Neumann series to its
///   second term, `x = a + C a + [-q^2 |a| / (1 - q), ...]`, `a = -Y R`, in the same forms:
///   exact for every s but the tail.
///
/// A span on one point is the same with no s.
///
/// # Errors
///
/// `InvalidInput` when `order_s` or `n` exceeds the series' size.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn gn_enclose_affine(
    k: &Kernel,
    faces_n: [&Nurbs; 2],
    leader_n: &NurbsCurve,
    at: &SpanAt,
    r: [f64; 2],
    side: [f64; 2],
    n: usize,
    order_s: usize,
    law: Option<&GnLaw<'_>>,
    defs: Option<&mut GnDefs<'_>>,
) -> Result<Option<Enclosure>> {
    let s = at.s;
    let tr = Trunc::new(order_s, n)?;
    let sm = s.mid();
    let rs = round_up(cmax(s.hi - sm, sm - s.lo));
    let wider = Trunc::new(cmin_usize(order_s + 1, MAX_S), 0)?;
    // a ball whose radius varies: its law over the span (0d43b4c), in s
    let law_r = match law {
        Some(l) => {
            let wide = gn_try!(law_series_over(l, s, wider.s()), None);
            let mid = gn_try!(law_series_over(l, Ivl::from(sm), order_s), None);
            Some(series_affine(&mid, &wide, order_s, rs))
        }
        None => None,
    };
    let law_r = law_r.as_ref();
    let (leader, tangent) = if at.piecewise && at.leader.is_none() {
        // across the leader's knots: the pieces' hull, with no dependence on s
        let (l, t) = gn_try!(leader_series_piecewise(leader_n, s, tr), None);
        (vser_aff(&l), vser_aff(&t))
    } else {
        let piece = at
            .leader
            .or_else(|| Some(span(&leader_n.u, leader_n.p, leader_n.n, sm)));
        let (lw, tw) = gn_try!(leader_series_over(leader_n, s, wider, piece), None);
        let (lm, tm) = gn_try!(leader_series_over(leader_n, Ivl::from(sm), tr, piece), None);
        (
            vseries_affine(&lm, &lw, order_s, rs),
            vseries_affine(&tm, &tw, order_s, rs),
        )
    };
    let ax = |c: usize| at.xa.get(c).copied().unwrap_or(f64::NAN);
    let bx = |c: usize| at.xb.get(c).copied().unwrap_or(f64::NAN);
    // the feet's predicted line x^(s) = xm + g (s - m) / rad, through the guesses
    let xm: [f64; 4] = core::array::from_fn(|c| 0.5 * (ax(c) + bx(c)));
    // a face whose partials have no form (a definition's, or a box across a knot line
    // taken piecewise) moves along the line only in intervals: there the feet are a box
    // about the middle, as the fork's
    let any_def = defs
        .as_ref()
        .is_some_and(|d| d.faces.iter().any(Option::is_some));
    let tube = rs > 0.0 && !at.piecewise && !any_def;
    let slope: [f64; 4] = core::array::from_fn(|c| if tube { 0.5 * (bx(c) - ax(c)) } else { 0.0 });
    let feet = |w: &[f64; 4]| -> [Aff; 4] {
        core::array::from_fn(|q| {
            let mut g = [0.0; NOISE];
            if let Some(gs) = g.get_mut(S_NOISE) {
                *gs = f_at(&slope, q);
            }
            if let Some(gq) = g.get_mut(q) {
                *gq = f_at(w, q);
            }
            Aff {
                c: f_at(&xm, q),
                g,
                r: 0.0,
            }
        })
    };
    let [ku0, ku1] = at.ku;
    let [kv0, kv1] = at.kv;
    let mut defs = defs;
    let mut partials_over = |x: &[Aff; 4], lk: usize| -> Result<[Partials<Aff>; 2]> {
        let [n0, n1] = faces_n;
        let [x0, x1, x2, x3] = *x;
        let (d0, d1) = match defs.as_deref_mut() {
            Some(d) => {
                let [a0, a1] = &mut d.faces;
                (a0.as_mut(), a1.as_mut())
            }
            None => (None, None),
        };
        let piecewise = at.piecewise && d0.is_none() && d1.is_none();
        Ok([
            face_partials_affine(k, n0, d0, piecewise, [ku0, kv0], [x0, x1], lk)?,
            face_partials_affine(k, n1, d1, piecewise, [ku1, kv1], [x2, x3], lk)?,
        ])
    };
    let consts = |x: &[Aff; 4]| -> [Ser<Aff>; 4] { x.map(Ser::constant) };
    // the preconditioner Y(s) = Y0 + Y1 e_4 from a Jacobian's middle and its slope in s
    let precondition = |jac: &[[Aff; 4]; 4]| -> Result<[[Aff; 4]; 4]> {
        let j0 = jac.map(|row| row.map(|x| x.c));
        let j1 = jac.map(|row| row.map(|x| x.g.get(S_NOISE).copied().unwrap_or(0.0)));
        let y0 = inverse4(&j0)?;
        let at = |m: &[[f64; 4]; 4], i: usize, j: usize| {
            m.get(i).and_then(|r| r.get(j)).copied().unwrap_or(f64::NAN)
        };
        let mul = |a: &[[f64; 4]; 4], b: &[[f64; 4]; 4]| -> [[f64; 4]; 4] {
            core::array::from_fn(|i| {
                core::array::from_fn(|j| (0..4).map(|l| at(a, i, l) * at(b, l, j)).sum::<f64>())
            })
        };
        let y1 = mul(&mul(&y0, &j1), &y0);
        Ok(core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                let mut g = [0.0; NOISE];
                if let Some(gs) = g.get_mut(S_NOISE) {
                    *gs = -at(&y1, i, j);
                }
                Aff {
                    c: at(&y0, i, j),
                    g,
                    r: 0.0,
                }
            })
        }))
    };
    // C = I - Y J and its largest row sum over every s and residual
    let contraction = |y: &[[Aff; 4]; 4], jac: &[[Aff; 4]; 4]| -> ([[Aff; 4]; 4], f64) {
        let mut c = [[Aff::point(0.0); 4]; 4];
        let mut q = 0.0_f64;
        for (i, (yi, ci)) in y.iter().zip(c.iter_mut()).enumerate() {
            let mut row = 0.0_f64;
            for (j, cij) in ci.iter_mut().enumerate() {
                let mut v = Aff::point(if i == j { 1.0 } else { 0.0 });
                for (&yil, jl) in yi.iter().zip(jac) {
                    v = v - yil * jl.get(j).copied().unwrap_or(Aff::WHOLE);
                }
                *cij = v;
                row = round_up(row + v.ivl().mag());
            }
            q = cmax(q, row);
        }
        (c, q)
    };
    let times = |m: &[[Aff; 4]; 4], v: &[Aff; 4]| -> [Aff; 4] {
        core::array::from_fn(|i| {
            let mut a = Aff::point(0.0);
            if let Some(mi) = m.get(i) {
                for (&mij, &vj) in mi.iter().zip(v) {
                    a = a + mij * vj;
                }
            }
            a
        })
    };

    // G along the predicted line, for every s in the span
    let point = Trunc::new(0, 0)?;
    let line = feet(&[0.0; 4]);
    let gl = {
        let pc = gn_try!(partials_over(&line, 2), None);
        let (g, _) = gn_try!(
            gn_evaluate(
                &system(&pc, &leader, &tangent, r, side, law_r),
                &consts(&line),
                point
            ),
            None
        );
        g.map(|x| x.at(0, 0))
    };
    // Krawczyk's operator about the line: the residual E from the guesses' spread
    // (about the line the guesses lie on it; about the middle, as the fork's box, twice
    // their spread, which they lie inside)
    let mut w: [f64; 4] = core::array::from_fn(|c| {
        let spread = 0.5 * (bx(c) - ax(c)).abs() + 1.0e-9 * (1.0 + f_at(&xm, c).abs());
        if tube {
            spread
        } else {
            2.0 * spread
        }
    });
    let mut proved: Option<[[Aff; 4]; 4]> = None;
    let mut it = 0;
    while it < 8 && proved.is_none() {
        it += 1;
        let x = feet(&w);
        let pb = gn_try!(partials_over(&x, 2), None);
        let jac = gn_try!(
            gn_jacobian(&system(&pb, &leader, &tangent, r, side, law_r), &consts(&x)),
            None
        );
        let y = gn_try!(precondition(&jac), None);
        let (c, q) = contraction(&y, &jac);
        let d = times(&y, &gl).map(|v| Aff::point(0.0) - v);
        let e: [Aff; 4] = core::array::from_fn(|j| {
            let mut g = [0.0; NOISE];
            if let Some(gj) = g.get_mut(j) {
                *gj = f_at(&w, j);
            }
            Aff { c: 0.0, g, r: 0.0 }
        });
        let ce = times(&c, &e);
        let mut inside = q < 1.0;
        let mut kw = [0.0_f64; 4];
        for (q_, ((di, ci), kq)) in d.iter().zip(&ce).zip(kw.iter_mut()).enumerate() {
            let kk = (*di + *ci).ivl();
            *kq = cmax(kk.lo.abs(), kk.hi.abs());
            let wq = f_at(&w, q_);
            if !(kk.lo > -wq && kk.hi < wq) {
                inside = false;
            }
        }
        if inside {
            // the root's residual lies in K: its radius bounds it
            w = kw;
            proved = Some(y);
        } else {
            let contracts = q < 1.0;
            if !contracts {
                return Ok(None); // the span is too long for its Jacobian to be one preconditioner's
            }
            // the next residual box from the contraction: K's radii are at most d + |C| w,
            // so the box of radii w = (I - |C|)^-1 d holds its K; a half more, to clear
            // rounding and the change of C with the box
            let dm = d.map(|v| v.ivl().mag());
            let cmag = c.map(|row| row.map(|v| v.ivl().mag()));
            let mut wn = dm;
            for _ in 0..200 {
                let mut nx = [0.0_f64; 4];
                for ((nxi, &di), ci) in nx.iter_mut().zip(&dm).zip(&cmag) {
                    *nxi = di;
                    for (cij, &wj) in ci.iter().zip(&wn) {
                        *nxi += cij * wj;
                    }
                }
                wn = nx;
            }
            for ((wq, &wi), &mi) in w.iter_mut().zip(&wn).zip(&xm) {
                let rr = 1.5 * wi + 1.0e-12 * (1.0 + mi.abs());
                if !rr.is_finite() {
                    return Ok(None);
                }
                *wq = rr;
            }
        }
    }
    let Some(y) = proved else {
        return Ok(None);
    };

    // the feet's Taylor coefficients over the span, order by order: J x = -R, so
    // x = a + C x, a = -Y R, C = I - Y J, whose Neumann series to its second term holds x
    // with its tail, |C^2 a + C^3 a + ...| <= q^2 |a| / (1 - q)
    // the partials the series take, and the Jacobian's second ones
    let lk = cmax_usize(tr.s() + tr.e() + 1, 2);
    let x0 = feet(&w);
    let pk = gn_try!(partials_over(&x0, lk), None);
    let sysk = system(&pk, &leader, &tangent, r, side, law_r);
    let mut sx = consts(&x0);
    let jk = gn_try!(gn_jacobian(&sysk, &sx), None);
    let (ck, qk) = contraction(&y, &jk);
    // 1 - q rounded down: the quotient rounded up is then at least the exact one
    let gap = round_down(1.0 - qk);
    if !gap.is_finite() || gap <= 0.0 {
        return Ok(None);
    }
    for tot in 1..=(tr.s() + tr.e()) {
        // the coefficients of total order tot take G's to that order in each variable
        // only: a truncated series' coefficients below its truncation are the full one's
        let to = Trunc::new(cmin_usize(tr.s(), tot), cmin_usize(tr.e(), tot))?;
        let (sg, _) = gn_try!(gn_evaluate(&sysk, &sx, to), None);
        for i in 0..=cmin_usize(tr.s(), tot) {
            let k = tot - i;
            // the section's jets of order n are the face curve's (gn_top_jets, 3bc54a4)
            if k > tr.e() || (n >= 1 && k >= n) {
                continue;
            }
            let rk: [Aff; 4] =
                core::array::from_fn(|q| sg.get(q).map_or(Aff::WHOLE, |g| g.at(i, k)));
            let a = times(&y, &rk).map(|v| Aff::point(0.0) - v);
            let amax = a.iter().fold(0.0_f64, |m, v| cmax(m, v.ivl().mag()));
            let tail = round_up(round_up(qk * qk) * amax / gap);
            if !tail.is_finite() {
                return Ok(None);
            }
            let ca = times(&ck, &a);
            for ((sxq, &aq), &caq) in sx.iter_mut().zip(&a).zip(&ca) {
                sxq.set(i, k, (aq + caq).widened(tail));
            }
        }
    }
    let (_, mut f, v, nrm) = gn_try!(gn_evaluate_vn(&sysk, &sx, tr), None);
    gn_top_jets(n, tr, &nrm, &mut f);
    // the faces' first partials over the region: |S_u| |du| + |S_v| |dv| bounds |dS| for
    // a uv step (du, dv) (#2690: kept apart, where the fork sums them)
    let mut ds = [[0.0_f64; 2]; 2];
    let mut ds_lo = [[0.0_f64; 2]; 2];
    for ((dsf, dsl), pf) in ds.iter_mut().zip(ds_lo.iter_mut()).zip(&pk) {
        let su = v3_ivl(pf.at(1, 0));
        let sv = v3_ivl(pf.at(0, 1));
        let du = crate::geom::vector3d::sqrt(
            su.x.mag() * su.x.mag() + su.y.mag() * su.y.mag() + su.z.mag() * su.z.mag(),
        );
        let dv = crate::geom::vector3d::sqrt(
            sv.x.mag() * sv.x.mag() + sv.y.mag() * sv.y.mag() + sv.z.mag() * sv.z.mag(),
        );
        *dsf = [du, dv];
        let lu = crate::geom::vector3d::sqrt(
            su.x.mig() * su.x.mig() + su.y.mig() * su.y.mig() + su.z.mig() * su.z.mig(),
        );
        let lv = crate::geom::vector3d::sqrt(
            sv.x.mig() * sv.x.mig() + sv.y.mig() * sv.y.mig() + sv.z.mig() * sv.z.mig(),
        );
        *dsl = [lu, lv];
    }
    let [f0, f1] = &f;
    Ok(Some(Enclosure {
        f: [vser_ivl(f0), vser_ivl(f1)],
        x: sx.map(|x| ser_ivl(&x)),
        ds,
        ds_lo,
        v: vser_ivl(&v),
    }))
}

/// A series in sigma over a span, `s = m + rad e_4`, as affine forms in s: each
/// coefficient `c_j` in its mean value form `c_j(m) + (j + 1) c_(j+1)(S) (s - m)` (the
/// coefficients of a Taylor series at s have `d c_j / ds = (j + 1) c_(j+1)`), from the
/// series at the middle (`mid`) and over the span (`wide`, one order more where the series
/// holds it; the top coefficient past that is its enclosure over the span).
fn series_affine(mid: &Ser<Ivl>, wide: &Ser<Ivl>, order: usize, rad: f64) -> Ser<Aff> {
    let mut out = Ser::constant(Aff::point(0.0));
    let mut ds = [0.0; NOISE];
    if let Some(g) = ds.get_mut(S_NOISE) {
        *g = rad;
    }
    let dev = Aff {
        c: 0.0,
        g: ds,
        r: 0.0,
    };
    for j in 0..=order {
        let c = if j < MAX_S {
            let slope = Ivl::from((j + 1) as f64) * wide.at(j + 1, 0);
            Aff::mean_value(mid.at(j, 0), &[(slope, dev)])
        } else {
            Aff::of_ivl(wide.at(j, 0))
        };
        out.set(j, 0, c);
    }
    out
}

/// [`series_affine`] for a vector series.
fn vseries_affine(mid: &VSer<Ivl>, wide: &VSer<Ivl>, order: usize, rad: f64) -> VSer<Aff> {
    let [m0, m1, m2] = &mid.a;
    let [w0, w1, w2] = &wide.a;
    VSer {
        a: [
            series_affine(m0, w0, order, rad),
            series_affine(m1, w1, order, rad),
            series_affine(m2, w2, order, rad),
        ],
    }
}

/// A series of intervals as affine forms with no dependence on the noise.
fn ser_aff(s: &Ser<Ivl>) -> Ser<Aff> {
    Ser {
        c: s.c.map(|row| row.map(Aff::of_ivl)),
    }
}

/// [`ser_aff`] for a vector series.
fn vser_aff(v: &VSer<Ivl>) -> VSer<Aff> {
    VSer {
        a: v.a.map(|s| ser_aff(&s)),
    }
}

/// The intervals a series of affine forms ranges over.
fn ser_ivl(s: &Ser<Aff>) -> Ser<Ivl> {
    Ser {
        c: s.c.map(|row| row.map(|x| x.ivl())),
    }
}

/// [`ser_ivl`] for a vector series.
fn vser_ivl(v: &VSer<Aff>) -> VSer<Ivl> {
    VSer {
        a: v.a.map(|s| ser_ivl(&s)),
    }
}

/// The intervals a vector of affine forms ranges over.
fn v3_ivl(v: V3<Aff>) -> V3<Ivl> {
    V3 {
        x: v.x.ivl(),
        y: v.y.ivl(),
        z: v.z.ivl(),
    }
}

/// Partials enclosed in intervals as affine forms with no dependence on the noise.
fn partials_aff(p: &Partials<Ivl>, lk: usize) -> Partials<Aff> {
    let mut out = Partials::default();
    for a in 0..=lk {
        for b in 0..=(lk - a) {
            let x = p.at(a, b);
            out.set(
                a,
                b,
                V3 {
                    x: Aff::of_ivl(x.x),
                    y: Aff::of_ivl(x.y),
                    z: Aff::of_ivl(x.z),
                },
            );
        }
    }
    out
}

/// A face's partials `P[a][b]`, `a + b <= lk`, at its foot `(u, v)`, given as affine forms
/// (`uv`), as affine forms: on one polynomial piece, their mean value forms `P_ab(m) +
/// P_(a+1)b(B) (u - m_u) + P_a(b+1)(B) (v - m_v)` about the forms' middles m, the partials
/// at m on that point and the next order's over the box B the forms range over
/// ([`nurbs_partials_over_piece`]). A face that is a blend by its definition, and a box
/// across a knot line of a face taken piecewise (`piecewise`,
/// [`nurbs_partials_over_box`]), are their interval enclosures over B, with no dependence
/// on the noise.
///
/// # Errors
///
/// Those of the enclosures.
pub(crate) fn face_partials_affine(
    k: &Kernel,
    cn: &Nurbs,
    def: Option<&mut GnDefFace<'_>>,
    piecewise: bool,
    piece: [usize; 2],
    uv: [Aff; 2],
    lk: usize,
) -> Result<Partials<Aff>> {
    let [ku, kv] = piece;
    let [fu, fv] = uv;
    let (bu, bv) = (fu.ivl(), fv.ivl());
    if let Some(d) = def {
        let p = definition_partials_over_box(k, d, cn, ku, kv, bu.lo, bu.hi, bv.lo, bv.hi, lk)?;
        return Ok(partials_aff(&p, lk));
    }
    let (ku, kv) = if piecewise {
        let us = knot_cuts(&cn.u, cn.p, cn.nu, bu.lo, bu.hi, FEET_BOX_CONTINUITY)?;
        let vs = knot_cuts(&cn.v, cn.q, cn.nv, bv.lo, bv.hi, FEET_BOX_CONTINUITY)?;
        if let ([(_, _, a)], [(_, _, b)]) = (us.as_slice(), vs.as_slice()) {
            (*a, *b)
        } else {
            let p =
                nurbs_partials_over_box(cn, bu.lo, bu.hi, bv.lo, bv.hi, lk, FEET_BOX_CONTINUITY)?;
            return Ok(partials_aff(&p, lk));
        }
    } else {
        (ku, kv)
    };
    let mid = nurbs_partials_direct(cn, ku, kv, fu.c, fu.c, fv.c, fv.c, lk)?;
    let point = bu.lo >= bu.hi && bv.lo >= bv.hi;
    if point {
        return Ok(partials_aff(&mid, lk));
    }
    if lk >= MAX_K {
        let p = nurbs_partials_over_piece(cn, ku, kv, bu.lo, bu.hi, bv.lo, bv.hi, lk)?;
        return Ok(partials_aff(&p, lk));
    }
    let wide = nurbs_partials_over_piece(cn, ku, kv, bu.lo, bu.hi, bv.lo, bv.hi, lk + 1)?;
    let du = fu - Aff::point(fu.c);
    let dv = fv - Aff::point(fv.c);
    let mut out = Partials::default();
    for a in 0..=lk {
        for b in 0..=(lk - a) {
            let (m, pu, pv) = (mid.at(a, b), wide.at(a + 1, b), wide.at(a, b + 1));
            let form = |m: Ivl, su: Ivl, sv: Ivl| Aff::mean_value(m, &[(su, du), (sv, dv)]);
            out.set(
                a,
                b,
                V3 {
                    x: form(m.x, pu.x, pv.x),
                    y: form(m.y, pu.y, pv.y),
                    z: form(m.z, pu.z, pv.z),
                },
            );
        }
    }
    Ok(out)
}

/// The larger of two orders.
#[inline]
fn cmax_usize(a: usize, b: usize) -> usize {
    if a < b {
        b
    } else {
        a
    }
}

/// `std::min` of the loop bound `i <= tr.s && i <= lTot`.
#[inline]
fn cmin_usize(a: usize, b: usize) -> usize {
    if b < a {
        b
    } else {
        a
    }
}

/// The fourth s-derivatives' bounds of a span: `m4[f][k]` per face and jet, `uv4[f]` per
/// face's uv curves (infinite where not proved), then the same bounds' least values at the
/// span's middle (`pM4Lo`, `pUV4Lo`): no piece holding the middle bounds them lower. Then
/// the uv curves' bound and least value as the fork takes them,
/// `(|S_u| + |S_v|) max(|u''''|, |v''''|)` (#2690: `gn_certify_span` steers its levels by
/// these, and certifies by `uv4`).
pub type Fourth = (
    [[f64; MAX_E + 1]; 2],
    [f64; 2],
    [[f64; MAX_E + 1]; 2],
    [f64; 2],
    [f64; 2],
    [f64; 2],
);

/// The Taylor-form bound `sum_{j<M} (4+j)! |c[4+j](m)| (h/2)^j / j!` plus the remainder
/// `(4+M)! |c[4+M](box)| (h/2)^M / M!`, `M = SM_GN_TAYLOR`, of one coefficient row.
fn taylor_bound(mid: impl Fn(usize) -> f64, rem: f64, h: f64) -> f64 {
    let mut m = 0.0_f64;
    let mut fact = 24.0_f64;
    let mut pow = 1.0_f64;
    let mut jf = 1.0_f64;
    for j in 0..TAYLOR {
        if j > 0 {
            fact *= (4 + j) as f64;
            pow *= 0.5 * h;
            jf *= j as f64;
        }
        m += fact * mid(4 + j) * pow / jf;
    }
    fact *= (4 + TAYLOR) as f64;
    pow *= 0.5 * h;
    jf *= TAYLOR as f64;
    m += fact * rem * pow / jf;
    m
}

/// `sm_GnFourthOver`: the bound on the jets' fourth s-derivatives over `[a.s, b.s]`, per
/// face and jet (`m4[f][k]`, infinite where not proved), and on the feet's uv curves',
/// each parameter's times the face's partial in it over the span's box,
/// `|S_u| |u''''| + |S_v| |v''''|` (`uv4[f]`; the fork's is
/// `(|S_u| + |S_v|) max(|u''''|, |v''''|)`, #2690), in their Taylor form at the
/// interval's middle m: `c''''(s)` is within `sum_{j<M} c^(4+j)(m) (s - m)^j / j!` plus the
/// remainder `c^(4+M)(interval) (h/2)^M / M!`, `M = SM_GN_TAYLOR`, the derivatives at the
/// middle enclosed on that one point and the remainder's over the interval.
///
/// # Errors
///
/// `InvalidInput` when `n` exceeds `SM_GN_MAXE`.
#[allow(clippy::too_many_arguments)]
pub fn gn_fourth_over(
    k: &Kernel,
    faces_n: [&Nurbs; 2],
    leader_n: &NurbsCurve,
    faces: [&Surface; 2],
    center: &BSplineCurve,
    a: &GnStation,
    b: &GnStation,
    r: [f64; 2],
    side: [f64; 2],
    n: usize,
) -> Result<Fourth> {
    let (fourth, _) = gn_fourth_over_series(
        k, faces_n, leader_n, faces, center, a, b, r, side, n, None, 0, None,
    )?;
    Ok(fourth)
}

/// The series `sm_GnFourthOver` hands a varying thumbweight's row terms (its `pF`, `pV`,
/// `pFm`, `pVm`, `pdMid`, 7102bc5): the feet's and V's over the span to order
/// `4 + SM_GN_TAYLOR` in s, at its middle to one less, and the middle.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PieceSeries {
    /// `pF`: the feet over the span.
    pub f: [VSer<Ivl>; 2],
    /// `pV`: V over the span.
    pub v: VSer<Ivl>,
    /// `pFm`: the feet at the middle.
    pub fm: [VSer<Ivl>; 2],
    /// `pVm`: V at the middle.
    pub vm: VSer<Ivl>,
    /// `pdMid`: the middle.
    pub mid: f64,
}

/// [`gn_fourth_over`] with the series of [`PieceSeries`] (`None` where they are not all
/// proved, as the fork's out-params are then not all set), and the fork's `pLaw` and
/// `lExtra` (0d43b4c): a ball whose radius varies, and the orders in s the series take
/// past `4 + SM_GN_TAYLOR` (the touch circle's m'); and its `pDefs` (c5b6bb7): the faces
/// that are blends by their definitions, enclosed over the boxes.
///
/// # Errors
///
/// `InvalidInput` when `n` exceeds `SM_GN_MAXE`.
#[allow(clippy::too_many_arguments)]
pub fn gn_fourth_over_series(
    k: &Kernel,
    faces_n: [&Nurbs; 2],
    leader_n: &NurbsCurve,
    faces: [&Surface; 2],
    center: &BSplineCurve,
    a: &GnStation,
    b: &GnStation,
    r: [f64; 2],
    side: [f64; 2],
    n: usize,
    law: Option<&GnLaw<'_>>,
    extra: usize,
    defs: Option<&mut GnDefs<'_>>,
) -> Result<(Fourth, Option<PieceSeries>)> {
    let mut defs = defs;
    let none: (Fourth, Option<PieceSeries>) = (
        (
            [[f64::INFINITY; MAX_E + 1]; 2],
            [f64::INFINITY; 2],
            [[0.0; MAX_E + 1]; 2],
            [0.0; 2],
            [f64::INFINITY; 2],
            [0.0; 2],
        ),
        None,
    );
    let h = b.s - a.s;
    let Some(e5) = gn_enclose(
        k,
        faces_n,
        leader_n,
        a,
        b,
        r,
        side,
        n,
        4 + TAYLOR + extra,
        law,
        defs.as_deref_mut(),
    )?
    else {
        return Ok(none);
    };
    // the middle's station, and the series there
    let mid = 0.5 * (a.s + b.s);
    let mut seed = [0.0_f64; 4];
    for (c, sc) in seed.iter_mut().enumerate() {
        let ac = a.x.get(c).map_or(f64::NAN, |x| x.at(0, 0));
        let bc = b.x.get(c).map_or(f64::NAN, |x| x.at(0, 0));
        *sc = 0.5 * (ac + bc);
    }
    let sm = gn_try!(center.evaluate_point(mid), none);
    let point = Trunc::new(0, 0)?;
    let tol = 1.0e-9 * smax(1.0, sm.length());
    let mid_r = gn_try!(law_series(law, mid, 0, false), none);
    let smid = gn_try!(
        gn_solve_station_counted(
            k,
            faces,
            center,
            mid,
            seed,
            r,
            side,
            point,
            tol,
            mid_r.as_ref()
        )
        .0,
        none
    );
    let Some(e4) = gn_enclose(
        k,
        faces_n,
        leader_n,
        &smid,
        &smid,
        r,
        side,
        n,
        4 + TAYLOR - 1 + extra,
        law,
        defs,
    )?
    else {
        return Ok(none);
    };
    let series = PieceSeries {
        f: e5.f,
        v: e5.v,
        fm: e4.f,
        vm: e4.v,
        mid,
    };
    let (mut m4, mut uv4, mut m4_lo, mut uv4_lo, mut uv4_fork, mut uv4_lo_fork) = none.0;
    // the uv curves' fourth derivatives, in the same Taylor form, each times the face's
    // partial in its parameter over the span's box (#2690)
    for (f, ((uvf, uvk), dsf)) in uv4
        .iter_mut()
        .zip(uv4_fork.iter_mut())
        .zip(&e5.ds)
        .enumerate()
    {
        let (mut sum, mut dmax) = (0.0_f64, 0.0_f64);
        for (c, &dsc) in ((2 * f)..(2 * f + 2)).zip(dsf) {
            let (Some(x4), Some(x5)) = (e4.x.get(c), e5.x.get(c)) else {
                return Ok((none.0, Some(series)));
            };
            let m = taylor_bound(|i| x4.at(i, 0).mag(), x5.at(4 + TAYLOR, 0).mag(), h);
            sum += m * dsc;
            dmax = cmax(dmax, m);
        }
        *uvf = sum;
        // (the fork's: the largest fourth derivative times |S_u| + |S_v|)
        *uvk = dmax * dsf.iter().sum::<f64>();
    }
    // the same, least: the middle's 4! c[4] alone, each times the face's least partial
    for (f, ((uvl, uvk), dsl)) in uv4_lo
        .iter_mut()
        .zip(uv4_lo_fork.iter_mut())
        .zip(&e4.ds_lo)
        .enumerate()
    {
        let (mut sum, mut dlo) = (0.0_f64, 0.0_f64);
        for (c, &dsc) in ((2 * f)..(2 * f + 2)).zip(dsl) {
            let Some(x4) = e4.x.get(c) else {
                return Ok((none.0, Some(series)));
            };
            sum += 24.0 * x4.at(4, 0).mig() * dsc;
            dlo = cmax(dlo, 24.0 * x4.at(4, 0).mig());
        }
        *uvl = sum;
        *uvk = dlo * dsl.iter().sum::<f64>();
    }
    for ((m4f, f4), f5) in m4.iter_mut().zip(&e4.f).zip(&e5.f) {
        for (kk, m4fk) in m4f.iter_mut().enumerate().take(n + 1) {
            let mut sq = 0.0_f64;
            for (s4, s5) in f4.a.iter().zip(&f5.a) {
                // c^(4+j)(m) = (4+j)! F[4+j]: the Taylor polynomial's bound over
                // |s - m| <= h/2, and the remainder
                let m = taylor_bound(|i| s4.at(i, kk).mag(), s5.at(4 + TAYLOR, kk).mag(), h);
                sq += m * m;
            }
            *m4fk = crate::geom::vector3d::sqrt(sq);
        }
    }
    for (m4lf, f4) in m4_lo.iter_mut().zip(&e4.f) {
        for (kk, m4lfk) in m4lf.iter_mut().enumerate().take(n + 1) {
            let mut sq = 0.0_f64;
            for s4 in &f4.a {
                let m = 24.0 * s4.at(4, kk).mig();
                sq += m * m;
            }
            *m4lfk = crate::geom::vector3d::sqrt(sq);
        }
    }
    Ok((
        (m4, uv4, m4_lo, uv4_lo, uv4_fork, uv4_lo_fork),
        Some(series),
    ))
}

/// `sm_GnVarying`: a section whose thumbweight varies along s (the equal area one,
/// resolved per station): its order, `r1 r2`, and I's enclosures on their grid.
#[derive(Clone, Debug, PartialEq)]
pub struct GnVarying {
    /// `n`.
    pub n: usize,
    /// `r1`, `r2` (the fork's `dRR` is their product, rounded; the certificate takes it
    /// in intervals from the radii).
    pub r: [f64; 2],
    /// `pCells`: I's enclosures, computed once per cell.
    pub cells: GnAreaCells,
}

/// `sm_GnRowTerms`: the largest `|(sum_ik dg_ikj c'_ik)|` (`x`) and `|W_j|` (`y`),
/// the least `W_j` (`w_lo`) per row, and the largest `|S - O|` (`r`).
#[derive(Clone, Debug, PartialEq)]
pub struct GnRowTerms {
    /// `X`.
    pub x: Vec<f64>,
    /// `Y`.
    pub y: Vec<f64>,
    /// `Wlo`.
    pub w_lo: Vec<f64>,
    /// `R`.
    pub r: f64,
}

/// `sm_Centred`: the coefficients `0..=SM_GN_TAY` of a series over a piece in their
/// Taylor form at its middle m: `c_i(s) = sum_{j < K - i} C(i + j, i) c_(i+j)(m) d^j +
/// C(K, i) c_K(xi) d^(K-i)`, `d = s - m`, from the middle's coefficients `mid[0..K-1]` and
/// the piece's `top` at order K.
#[must_use]
pub fn centred(mid: &[Ivl], top: Ivl, kk: usize, d: Ivl) -> Tay<Ivl> {
    let mut r = Tay::<Ivl>::zero();
    for i in 0..=TAY {
        let mut sum = Ivl::from(0.0);
        let mut pow = Ivl::from(1.0);
        let mut j = 0;
        while i + j < kk {
            if j > 0 {
                pow = pow * d;
            }
            let m = mid.get(i + j).copied().unwrap_or(IZ);
            sum = sum + Ivl::from(binom(i + j, i)) * m * pow;
            j += 1;
        }
        pow = pow * d;
        r.set(i, sum + Ivl::from(binom(kk, i)) * top * pow);
    }
    r
}

/// A vector series' coordinate q's coefficient `(i, k)`.
fn vat(v: &VSer<Ivl>, q: usize, i: usize, k: usize) -> Ivl {
    v.a.get(q).map_or(IZ, |c| c.at(i, k))
}

/// `sm_GnRowTermsOver`: the row terms over one piece, folded into `rt` (largest X, Y, R;
/// least W); `false` where w's enclosure or a term is not proved. `f`, `centre`: the
/// feet's and V's series over the piece (to order `K = 4 + SM_GN_TAYLOR` in s); `fm`,
/// `cm`: at its middle (to K - 1); `d`: s less the middle over the piece.
#[allow(clippy::too_many_arguments)]
pub fn gn_row_terms_over(
    k: &Kernel,
    vary: &mut GnVarying,
    f: &[VSer<Ivl>; 2],
    centre: &VSer<Ivl>,
    fm: &[VSer<Ivl>; 2],
    cm: &VSer<Ivl>,
    d: Ivl,
    guess: f64,
    bar: &GnSectionFormT<Ivl>,
    o: [f64; 3],
    rt: &mut GnRowTerms,
) -> bool {
    const K: usize = 4 + TAYLOR;
    let n = vary.n;
    let [f0, f1] = f;
    let [fm0, fm1] = fm;
    // x = (F1 - V).(F2 - V) / (r1 r2): at the middle to order K - 1, over the piece at
    // order K
    let mut x_mid = [IZ; K];
    let mut x_top = IZ;
    for q in 0..3 {
        for i in 0..=K {
            for l in 0..=i {
                if i < K {
                    let t = (vat(fm0, q, l, 0) - vat(cm, q, l, 0))
                        * (vat(fm1, q, i - l, 0) - vat(cm, q, i - l, 0));
                    if let Some(xi) = x_mid.get_mut(i) {
                        *xi = *xi + t;
                    }
                } else {
                    x_top = x_top
                        + (vat(f0, q, l, 0) - vat(centre, q, l, 0))
                            * (vat(f1, q, i - l, 0) - vat(centre, q, i - l, 0));
                }
            }
        }
    }
    let [r1, r2] = vary.r;
    let inv_rr = Ivl::from(1.0) / (Ivl::from(r1) * Ivl::from(r2));
    for xi in &mut x_mid {
        *xi = *xi * inv_rr;
    }
    x_top = x_top * inv_rr;
    let u = gn_half_angle(&centred(&x_mid, x_top, K, d));
    let tg = gn_target(&u);
    let phi = vary.cells.at(k, guess) - tg.at(0);
    let mut w0 = IZ;
    let mut root = false;
    let mut rho = 1.0e-12 * (1.0 + guess.abs());
    let mut it = 0;
    while it < 40 && !root {
        let w = Ivl::new(guess - rho, guess + rho);
        let below_one = w.hi < 1.0;
        if !below_one {
            break;
        }
        let d1 = vary.cells.over(k, w).at(1);
        if d1.hi < 0.0 {
            let nn = Ivl::from(guess) - phi / d1;
            if nn.lo > w.lo && nn.hi < w.hi {
                w0 = nn;
                root = true;
            }
        }
        it += 1;
        rho *= 4.0;
    }
    if !root {
        return false;
    }
    let over_w0 = vary.cells.over(k, w0);
    let sw = gn_invert(&tg, &over_w0, w0);
    let form = gn_form_normalized(n, sw);
    // the jets c'_ik as series in s over the piece
    let mut c = [[[Tay::<Ivl>::zero(); 3]; MAX_E + 1]; 2];
    let mut r_max = 0.0_f64;
    for (cf, (ff, ffm)) in c.iter_mut().zip(f.iter().zip(fm)) {
        let mut sum = 0.0_f64;
        for (kk, ck) in cf.iter_mut().enumerate().take(n + 1) {
            let mut sq = 0.0_f64;
            for (q, ckq) in ck.iter_mut().enumerate() {
                let mut ac_mid = [IZ; K];
                for (i, a) in ac_mid.iter_mut().enumerate() {
                    *a = vat(ffm, q, i, kk);
                }
                *ckq = centred(&ac_mid, vat(ff, q, K, kk), K, d);
                if kk == 0 {
                    let oq = o.get(q).copied().unwrap_or(f64::NAN);
                    ckq.set(0, ckq.at(0) - Ivl::from(oq));
                }
                sq += ckq.at(0).mag() * ckq.at(0).mag();
            }
            sum += crate::geom::vector3d::sqrt(sq);
        }
        r_max = cmax(r_max, sum);
    }
    for (j, wj) in form.w.iter().enumerate() {
        let mut sq = 0.0_f64;
        for q in 0..3 {
            // the coefficient of sigma^4 of sum_ik dg_ikj c'_ik
            let mut s4 = IZ;
            for ((gf, bf), cf) in form.g.iter().zip(&bar.g).zip(&c) {
                for ((gk, bk), ck) in gf.iter().zip(bf).zip(cf).take(n + 1) {
                    let (Some(&g), Some(&b), Some(ckq)) = (gk.get(j), bk.get(j), ck.get(q)) else {
                        continue;
                    };
                    let mut dg = g;
                    dg.set(0, dg.at(0) - b);
                    for i in 0..=TAY {
                        s4 = s4 + dg.at(i) * ckq.at(TAY - i);
                    }
                }
            }
            sq += s4.mag() * s4.mag();
        }
        let dx = 24.0 * crate::geom::vector3d::sqrt(sq);
        let dy = 24.0 * wj.at(TAY).mag();
        let dwlo = wj.at(0).lo;
        if !(dx.is_finite() && dy.is_finite() && dwlo > 0.0) {
            return false;
        }
        let (Some(xj), Some(yj), Some(lj)) = (rt.x.get_mut(j), rt.y.get_mut(j), rt.w_lo.get_mut(j))
        else {
            return false;
        };
        *xj = cmax(*xj, dx);
        *yj = cmax(*yj, dy);
        *lj = cmin(*lj, dwlo);
    }
    if !r_max.is_finite() {
        return false;
    }
    rt.r = cmax(rt.r, r_max);
    true
}

/// The arc's homogeneous rows `H[j][q]` and weights `W[j]` ([`arc_rows`]).
pub type ArcRows<S> = ([[Tay<S>; 3]; 3], [Tay<S>; 3]);

/// `sm_ArcRows` (cc08236): the ball's circle as rows: the rational quadratic arc from the
/// rail `rho1` to `rho2` about the ball's centre `m`, its control points `rho1, P1, rho2`
/// and weights `1, w, 1`, `w = cos(theta / 2)` (theta the angle between the radii
/// `b_i = rho_i - m`), `P1 = m + (mid - m) / w^2` where the rails' tangents meet; as
/// homogeneous rows `H_j = W_j (P_j - O)` about a fixed point `o`. The series are in s,
/// over a ring S (`f64`, or intervals). `None` where a radius is not positive or the arc
/// is a half turn or more (no single arc).
#[must_use]
pub fn arc_rows<S: TayScalar + GnScalar>(
    rho1: &[Tay<S>; 3],
    rho2: &[Tay<S>; 3],
    m: &[Tay<S>; 3],
    o: [f64; 3],
) -> Option<ArcRows<S>> {
    let mut b11 = Tay::<S>::zero();
    let mut b22 = Tay::<S>::zero();
    let mut b12 = Tay::<S>::zero();
    for ((&r1, &r2), &mq) in rho1.iter().zip(rho2).zip(m) {
        let b1 = r1 - mq;
        let b2 = r2 - mq;
        b11 = b11 + b1 * b1;
        b22 = b22 + b2 * b2;
        b12 = b12 + b1 * b2;
    }
    if !b11.at(0).is_positive() || !b22.at(0).is_positive() {
        return None;
    }
    // cos^2(theta / 2) = (1 + cos theta) / 2
    let half = (Tay::<S>::lift(1.0) + b12 / tay_sqrt(&(b11 * b22))) * Tay::<S>::lift(0.5);
    if !half.at(0).is_positive() {
        return None; // a half turn or more: no single arc
    }
    let w = tay_sqrt(&half);
    let inv = Tay::<S>::lift(1.0) / half;
    let ws = [Tay::<S>::lift(1.0), w, Tay::<S>::lift(1.0)];
    let mut h = [[Tay::<S>::zero(); 3]; 3];
    for (q, (((&r1, &r2), &mq), &oq)) in rho1.iter().zip(rho2).zip(m).zip(&o).enumerate() {
        let so = Tay::<S>::lift(oq);
        let p1 = mq + ((r1 + r2) * Tay::<S>::lift(0.5) - mq) * inv;
        let [h0, h1, h2] = &mut h;
        if let (Some(a), Some(b), Some(c)) = (h0.get_mut(q), h1.get_mut(q), h2.get_mut(q)) {
            *a = r1 - so;
            *b = w * (p1 - so);
            *c = r2 - so;
        }
    }
    Some((h, ws))
}

/// `sm_ArcRowTermsOver` (cc08236): the arc's row terms over one piece
/// ([`gn_row_terms_over`]'s, every row varying), folded into `rt`; `false` where they are
/// not proved. `f`, `centre`: the rails' and the centre's series over the piece (to order
/// `K = 4 + SM_GN_TAYLOR`); `fm`, `cm`: at its middle (to K - 1); `d`: s less the middle
/// over the piece. `touch` (0d43b4c): the ball's radius varies, and the circle is the one
/// it touches the envelope along, about `m + ((rho1 - m) . m' / m' . m') m'` (`centre` and
/// `cm` one order more, for m').
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn arc_row_terms_over(
    f: &[VSer<Ivl>; 2],
    centre: &VSer<Ivl>,
    fm: &[VSer<Ivl>; 2],
    cm: &VSer<Ivl>,
    d: Ivl,
    o: [f64; 3],
    rt: &mut GnRowTerms,
    touch: bool,
) -> bool {
    const K: usize = 4 + TAYLOR;
    let mut rho = [[Tay::<Ivl>::zero(); 3]; 2];
    let mut mc = [Tay::<Ivl>::zero(); 3];
    for q in 0..3 {
        let mut ac_mid = [IZ; K];
        for (rf, (ff, ffm)) in rho.iter_mut().zip(f.iter().zip(fm)) {
            for (i, a) in ac_mid.iter_mut().enumerate() {
                *a = vat(ffm, q, i, 0);
            }
            if let Some(x) = rf.get_mut(q) {
                *x = centred(&ac_mid, vat(ff, q, K, 0), K, d);
            }
        }
        for (i, a) in ac_mid.iter_mut().enumerate() {
            *a = vat(cm, q, i, 0);
        }
        if let Some(x) = mc.get_mut(q) {
            *x = centred(&ac_mid, vat(centre, q, K, 0), K, d);
        }
    }
    if touch {
        // m' from the centre's coefficients one order up: (i + 1) c_(i+1)
        let mut dm = [Tay::<Ivl>::zero(); 3];
        for (q, dq) in dm.iter_mut().enumerate() {
            let mut ac_mid = [IZ; K];
            for (i, a) in ac_mid.iter_mut().enumerate() {
                *a = vat(cm, q, i + 1, 0) * Ivl::from((i + 1) as f64);
            }
            *dq = centred(
                &ac_mid,
                vat(centre, q, K + 1, 0) * Ivl::from((K + 1) as f64),
                K,
                d,
            );
        }
        let mut dd = Tay::<Ivl>::zero();
        let mut bd = Tay::<Ivl>::zero();
        for ((&dq, &mq), &r0q) in dm.iter().zip(&mc).zip(&rho[0]) {
            dd = dd + dq * dq;
            bd = bd + (r0q - mq) * dq;
        }
        if !dd.at(0).is_positive() {
            return false;
        }
        let kk = bd / dd;
        for (mq, &dq) in mc.iter_mut().zip(&dm) {
            *mq = *mq + kk * dq;
        }
    }
    let [rho0, rho1] = &rho;
    let Some((h, w)) = arc_rows(rho0, rho1, &mc, o) else {
        return false;
    };
    // |S - O| <= |m - O| + the ball's radius over the piece
    let mut d_mo = 0.0_f64;
    let mut d_rad = 0.0_f64;
    for ((mq, &oq), r0q) in mc.iter().zip(&o).zip(rho0) {
        let dd = (mq.at(0) - Ivl::from(oq)).mag();
        d_mo += dd * dd;
        let b = (r0q.at(0) - mq.at(0)).mag();
        d_rad += b * b;
    }
    let d_r = crate::geom::vector3d::sqrt(d_mo) + crate::geom::vector3d::sqrt(d_rad);
    for (j, (hj, wj)) in h.iter().zip(&w).enumerate() {
        let mut sq = 0.0_f64;
        for hq in hj {
            sq += hq.at(TAY).mag() * hq.at(TAY).mag();
        }
        let dx = 24.0 * crate::geom::vector3d::sqrt(sq);
        let dy = 24.0 * wj.at(TAY).mag();
        let dwlo = wj.at(0).lo;
        if !(dx.is_finite() && dy.is_finite() && dwlo > 0.0) {
            return false;
        }
        let (Some(xj), Some(yj), Some(lj)) = (rt.x.get_mut(j), rt.y.get_mut(j), rt.w_lo.get_mut(j))
        else {
            return false;
        };
        *xj = cmax(*xj, dx);
        *yj = cmax(*yj, dy);
        *lj = cmin(*lj, dwlo);
    }
    if !d_r.is_finite() {
        return false;
    }
    rt.r = cmax(rt.r, d_r);
    true
}

/// `sm_GnCertifySpan`: the certified bound on the Hermite span between two stations, or
/// infinity where the span is not proved. The Hermite remainder is `h^4 / 384` times the
/// jets' largest fourth s-derivative over the span, bounded on its pieces
/// ([`gn_fourth_over`] on each), so that a long span's bound follows its derivatives, not
/// the width of one enclosure over it.
///
/// Not upstream (#2690): the fork cuts the span in `2^L` equal pieces, L raised from 0 (to
/// 6 at most) until every piece is proved and, over the target, while that halves the
/// bound, enclosing every piece of every level afresh: each level's work repeats the last
/// one's ends and refines where the bound is already met as much as where it is not. Here
/// the span's pieces are cut as the bound needs (Moore and Skelboe's refinement of the
/// piece that holds the maximum): a piece that is not proved is cut in two, and over the
/// target the pieces whose own terms keep the bound over it are (all of them where none
/// does on its own), each piece enclosed once and each station solved once, to `2^-6` of
/// the span at most, under the same stops (proved within the target, the least terms over
/// it, or a cut that does not halve the bound). The feet's uv curves are Hermite interpolants of the same stations, and the face at
/// a uv curve's point is off the face's true foot by at most `|S_u, S_v|` times the uv
/// error, `h^4 / 384 |x''''|`: that bound counts too (an edge's curve against its face).
///
/// Every level's bound is at least the same sum with each term's least value at the
/// pieces' middles (every piece's enclosure holds its middle's value): once that is over
/// the target no finer level certifies the span, and the bound so far is returned for the
/// caller to halve the span.
///
/// # Errors
///
/// `InvalidInput` when `n` exceeds `SM_GN_MAXE`.
#[allow(clippy::too_many_arguments)]
pub fn gn_certify_span(
    k: &Kernel,
    faces_n: [&Nurbs; 2],
    leader_n: &NurbsCurve,
    faces: [&Surface; 2],
    center: &BSplineCurve,
    a: &GnStation,
    b: &GnStation,
    r: [f64; 2],
    side: [f64; 2],
    n: usize,
    target: f64,
    mut vary: Option<&mut GnVarying>,
    law: Option<&GnLaw<'_>>,
    mut defs: Option<&mut GnDefs<'_>>,
) -> Result<f64> {
    crate::profile::profile_scope!(k, "blend::certificate::gn_certify_span");
    // the series' orders are checked once (the fork leaves them to its callers)
    Trunc::new(4 + TAYLOR, n)?;
    let mut bound_out = f64::INFINITY;
    let h = b.s - a.s;
    let mut last = f64::INFINITY;
    // the varying thumbweight's fixed one wb and point O: the first station's
    let mut bar = GnSectionFormT::<Ivl> {
        w: Vec::new(),
        g: [Vec::new(), Vec::new()],
    };
    let mut o = [0.0_f64; 3];
    // the ball's circle (n = 1, cc08236): its arc's rows all vary along s, about the first
    // station's centre
    let arc = n == 1;
    let rows_vary = vary.is_some() || arc;
    let l_jets = if arc { 0 } else { n }; // the jets graded on their own: the arc's rails only
    if let Some(v) = vary.as_deref() {
        if !(a.w.is_finite() && b.w.is_finite()) {
            return Err(SolidError::sm_err(OP));
        }
        bar = gn_form_normalized(v.n, Ivl::from(a.w));
        for (q, oq) in o.iter_mut().enumerate() {
            *oq = a.f[0].a.get(q).map_or(f64::NAN, |c| c.at(0, 0));
        }
    }
    if arc {
        bar.w = vec![Ivl::from(0.0); 3];
        for (q, oq) in o.iter_mut().enumerate() {
            *oq = a.v.a.get(q).map_or(f64::NAN, |c| c.at(0, 0));
        }
    }
    // the pieces the span is cut in, each enclosed once: refined where its own terms keep
    // the bound over the target, never deeper than 2^-6 of the span (the fork's levels)
    let mut pieces = vec![CertPiece {
        a: *a,
        b: *b,
        depth: 0,
        terms: None,
    }];
    let ld = bar.w.len();
    // whether the last cut was of some of the pieces only
    let mut partial = false;
    loop {
        k.check_cancelled(OP_CERTIFY_SPAN)?;
        // enclose the pieces not enclosed yet; one that is not proved is cut in two
        let mut i = 0;
        while i < pieces.len() {
            k.check_cancelled(OP_CERTIFY_SPAN)?;
            let Some(p) = pieces.get(i) else {
                break;
            };
            if p.terms.is_some() {
                i += 1;
                continue;
            }
            let (ea, eb, depth) = (p.a, p.b, p.depth);
            k.work_counters.add_certificate_piece(u64::from(depth));
            let ((pm, puv, pm_lo, _, puv_fork, puv_lo_fork), series) = gn_try!(
                gn_fourth_over_series(
                    k,
                    faces_n,
                    leader_n,
                    faces,
                    center,
                    &ea,
                    &eb,
                    r,
                    side,
                    l_jets,
                    law,
                    // the touch circle of a varying ball takes the centre's
                    // s-derivative: one order more
                    usize::from(arc && law.is_some()),
                    defs.as_deref_mut(),
                ),
                bound_out
            );
            let mut proved = pm
                .iter()
                .all(|pmf| pmf.iter().take(l_jets + 1).all(|&x| x < f64::INFINITY))
                && puv.iter().all(|&x| x < f64::INFINITY);
            let mut rows = GnRowTerms {
                x: vec![0.0; ld],
                y: vec![0.0; ld],
                w_lo: vec![f64::INFINITY; ld],
                r: 0.0,
            };
            if proved {
                if let Some(v) = vary.as_deref_mut() {
                    // the root's guess at the piece's middle: the stations' thumbweights,
                    // linearly
                    proved = match series {
                        Some(ps) => {
                            let t = (ps.mid - a.s) / h;
                            let d = Ivl::new(ea.s, eb.s) - Ivl::from(ps.mid);
                            let guess = (1.0 - t) * a.w + t * b.w;
                            gn_row_terms_over(
                                k, v, &ps.f, &ps.v, &ps.fm, &ps.vm, d, guess, &bar, o, &mut rows,
                            )
                        }
                        None => false,
                    };
                }
            }
            if proved && arc {
                proved = match series {
                    Some(ps) => {
                        let d = Ivl::new(ea.s, eb.s) - Ivl::from(ps.mid);
                        arc_row_terms_over(
                            &ps.f,
                            &ps.v,
                            &ps.fm,
                            &ps.vm,
                            d,
                            o,
                            &mut rows,
                            law.is_some(),
                        )
                    }
                    None => false,
                };
            }
            if proved {
                if let Some(p) = pieces.get_mut(i) {
                    p.terms = Some(PieceTerms {
                        m4: pm,
                        uv4: puv,
                        uv4_fork: puv_fork,
                        m4_lo: pm_lo,
                        uv4_lo_fork: puv_lo_fork,
                        rows,
                    });
                }
                i += 1;
                continue;
            }
            // not proved: cut in two at a station solved in the middle, while there is room
            if depth >= CERT_DEPTH {
                return Ok(bound_out);
            }
            let Some(mid) =
                cert_middle(k, faces, center, a, b, h, 0.5 * (ea.s + eb.s), r, side, law)?
            else {
                return Ok(bound_out);
            };
            if let Some(p) = pieces.get_mut(i) {
                p.b = mid;
                p.depth = depth + 1;
            }
            pieces.insert(
                i + 1,
                CertPiece {
                    a: mid,
                    b: eb,
                    depth: depth + 1,
                    terms: None,
                },
            );
        }
        // the span's bound over its pieces, and each piece's own
        let terms: Vec<&PieceTerms> = pieces.iter().filter_map(|p| p.terms.as_ref()).collect();
        // steered by the fork's uv term, certified by each parameter's own (#2690)
        let bound = span_bound(&terms, h, &bar, l_jets, rows_vary, arc, true);
        // a bound that is not proved (varying rows whose remainder is not yet below their
        // least weight) is no stop: every piece is cut, as the fork's next level
        let finite = bound.is_finite();
        let mut halved = false;
        if finite {
            let own = span_bound(&terms, h, &bar, l_jets, rows_vary, arc, false);
            bound_out = cmin(bound_out, round_up(own * (1.0 + 1.0e-12)));
            if bound_out <= target {
                break;
            }
            // the bound no finer cut goes under: over the target, no cut certifies the
            // span (the jets' least terms only where the bound takes the jets' sum, not
            // with varying rows, whose own least terms are not ported: the fork's f916a8d)
            let k4 = h * h * h * h / 384.0;
            let mut lower = 0.0_f64;
            for t in &terms {
                for (m4lf, &uvl) in t.m4_lo.iter().zip(&t.uv4_lo_fork) {
                    let mut sum = 0.0_f64;
                    for &m in m4lf.iter().take(l_jets + 1) {
                        sum += k4 * m;
                    }
                    lower = cmax(lower, cmax(if rows_vary { 0.0 } else { sum }, k4 * uvl));
                }
            }
            if round_down(lower * (1.0 - 1.0e-12)) > target {
                break;
            }
            // once proved, finer pieces only where the last cut halved the bound; a cut of
            // some of them that did not is followed by one of all of them, as the fork's
            // next level
            halved = bound <= 0.5 * last;
            if !halved && !partial {
                break;
            }
            if halved {
                last = bound;
            }
        }
        // cut the pieces whose own terms keep the bound over the target (all of them where
        // none does on its own, or the last cut of some did not halve the bound), each at a
        // station solved in its middle
        let alone: Vec<f64> = pieces
            .iter()
            .map(|p| {
                p.terms.as_ref().map_or(f64::INFINITY, |t| {
                    span_bound(&[t], h, &bar, l_jets, rows_vary, arc, true)
                })
            })
            .collect();
        let any = finite && halved && alone.iter().any(|&x| x > target);
        partial = any && alone.iter().any(|&x| x <= target);
        // a piece at the deepest cut stays as it is; with none left to cut, no finer cut
        // certifies the span
        let cuttable = pieces
            .iter()
            .zip(&alone)
            .any(|(p, &x)| p.depth < CERT_DEPTH && !(any && x <= target));
        if !cuttable {
            break;
        }
        let mut next = Vec::with_capacity(2 * pieces.len());
        for (p, &x) in pieces.into_iter().zip(&alone) {
            if (any && x <= target) || p.depth >= CERT_DEPTH {
                next.push(p);
                continue;
            }
            k.check_cancelled(OP_CERTIFY_SPAN)?;
            let Some(mid) = cert_middle(
                k,
                faces,
                center,
                a,
                b,
                h,
                0.5 * (p.a.s + p.b.s),
                r,
                side,
                law,
            )?
            else {
                return Ok(bound_out);
            };
            next.push(CertPiece {
                a: p.a,
                b: mid,
                depth: p.depth + 1,
                terms: None,
            });
            next.push(CertPiece {
                a: mid,
                b: p.b,
                depth: p.depth + 1,
                terms: None,
            });
        }
        pieces = next;
    }
    Ok(bound_out)
}

/// The deepest cut of a span the certificate takes: its pieces are at least `2^-6` of it
/// (the fork's levels `0..=6`).
const CERT_DEPTH: u32 = 6;

/// One piece of a span [`gn_certify_span`] bounds, between two stations, at `depth` cuts
/// from the span, and its terms once enclosed.
struct CertPiece {
    a: GnStation,
    b: GnStation,
    depth: u32,
    terms: Option<PieceTerms>,
}

/// A piece's terms: the jets' and the uv curves' fourth derivatives' bounds and their least
/// values ([`Fourth`]), and its varying rows' terms.
struct PieceTerms {
    m4: [[f64; MAX_E + 1]; 2],
    /// The uv curves' term per parameter, which certifies (#2690).
    uv4: [f64; 2],
    /// The fork's, which steers the cuts, with its least value.
    uv4_fork: [f64; 2],
    m4_lo: [[f64; MAX_E + 1]; 2],
    uv4_lo_fork: [f64; 2],
    rows: GnRowTerms,
}

/// The station at `s` inside the span `a`, `b` (of length `h`), seeded from the span's
/// end stations linearly; `None` where it does not solve (no piece is cut there).
///
/// # Errors
///
/// A cancel (#2136); those of the radius law.
#[allow(clippy::too_many_arguments)]
fn cert_middle(
    k: &Kernel,
    faces: [&Surface; 2],
    center: &BSplineCurve,
    a: &GnStation,
    b: &GnStation,
    h: f64,
    s: f64,
    r: [f64; 2],
    side: [f64; 2],
    law: Option<&GnLaw<'_>>,
) -> Result<Option<GnStation>> {
    let t = (s - a.s) / h;
    let mut seed = [0.0_f64; 4];
    for (c, sc) in seed.iter_mut().enumerate() {
        let ac = a.x.get(c).map_or(f64::NAN, |x| x.at(0, 0));
        let bc = b.x.get(c).map_or(f64::NAN, |x| x.at(0, 0));
        *sc = (1.0 - t) * ac + t * bc;
    }
    let sm = gn_try!(center.evaluate_point(s), None);
    let point = Trunc::new(0, 0)?;
    let tol = 1.0e-9 * smax(1.0, sm.length());
    let end_r = gn_try!(law_series(law, s, 0, false), None);
    match gn_solve_station_counted(
        k,
        faces,
        center,
        s,
        seed,
        r,
        side,
        point,
        tol,
        end_r.as_ref(),
    )
    .0
    {
        Ok(st) => Ok(Some(st)),
        // a cancel (#2136) is not an unsolved piece: it ends the operation
        Err(e @ SolidError::Cancelled { .. }) => Err(e),
        Err(_) => Ok(None),
    }
}

/// The Hermite span's bound from its pieces' terms (`h^4 / 384` times the largest fourth
/// derivatives, and with varying rows their quotient), infinite where it is not proved.
fn span_bound(
    terms: &[&PieceTerms],
    h: f64,
    bar: &GnSectionFormT<Ivl>,
    l_jets: usize,
    rows_vary: bool,
    arc: bool,
    fork_uv: bool,
) -> f64 {
    let mut m4 = [[0.0_f64; MAX_E + 1]; 2];
    let mut uv4 = [0.0_f64; 2];
    let ld = bar.w.len();
    let mut rows = GnRowTerms {
        x: vec![0.0; ld],
        y: vec![0.0; ld],
        w_lo: vec![f64::INFINITY; ld],
        r: 0.0,
    };
    for t in terms {
        for (m4f, pmf) in m4.iter_mut().zip(&t.m4) {
            for (m, &p) in m4f.iter_mut().zip(pmf).take(l_jets + 1) {
                *m = cmax(*m, p);
            }
        }
        let tu = if fork_uv { &t.uv4_fork } else { &t.uv4 };
        for (u, &p) in uv4.iter_mut().zip(tu) {
            *u = cmax(*u, p);
        }
        for (x, &p) in rows.x.iter_mut().zip(&t.rows.x) {
            *x = cmax(*x, p);
        }
        for (y, &p) in rows.y.iter_mut().zip(&t.rows.y) {
            *y = cmax(*y, p);
        }
        for (l, &p) in rows.w_lo.iter_mut().zip(&t.rows.w_lo) {
            *l = cmin(*l, p);
        }
        rows.r = cmax(rows.r, t.rows.r);
    }
    let mut bound = 0.0_f64;
    let mut da = 0.0_f64;
    for (m4f, &uvf) in m4.iter().zip(&uv4) {
        let mut sum = 0.0_f64;
        for &m4fk in m4f.iter().take(l_jets + 1) {
            sum += h * h * h * h / 384.0 * m4fk;
        }
        da = cmax(da, sum);
        bound = cmax(
            bound,
            cmax(
                if rows_vary { 0.0 } else { sum },
                h * h * h * h / 384.0 * uvf,
            ),
        );
    }
    // the arc: no fixed section, every row's remainder (the varying rows' bound with
    // A = 0)
    if arc {
        da = 0.0;
    }
    if rows_vary {
        let kh = h * h * h * h / 384.0;
        for (j, wj) in bar.w.iter().enumerate() {
            let (Some(&lj), Some(&yj), Some(&xj)) =
                (rows.w_lo.get(j), rows.y.get(j), rows.x.get(j))
            else {
                break;
            };
            let den = lj - kh * yj;
            let num = wj.hi * da + kh * xj + kh * yj * rows.r;
            let positive = den > 0.0;
            if !positive {
                return f64::INFINITY;
            }
            bound = cmax(bound, num / den);
        }
    }
    bound
}

#[cfg(test)]
#[allow(clippy::needless_range_loop)]
mod tests {
    use super::super::station::{as_b_spline, nurbs_partials, vec3};
    use super::*;
    use crate::geom::axis2_placement::Axis2Placement;
    use crate::geom::extent2d::Extent2d;
    use crate::math::{Vec2, Vec3};
    use crate::nurbs::b_spline_curve::NurbCircleParam;
    use crate::surface::cone::Cone;
    use crate::surface::plane::Plane;

    fn plane(origin: Vec3, x_axis: Vec3, y_axis: Vec3) -> Surface {
        Surface::from(Plane::new(
            origin,
            x_axis,
            y_axis,
            Vec2::new(1.0, 1.0),
            &Extent2d::from_bounds(-10.0, -10.0, 10.0, 10.0),
        ))
    }

    fn nurbs(k: &Kernel, s: &Surface) -> Nurbs {
        Nurbs::of(&as_b_spline(k, s).unwrap().unwrap()).unwrap()
    }

    fn station(
        k: &Kernel,
        faces: [&Surface; 2],
        leader: &BSplineCurve,
        s: f64,
        seed: [f64; 4],
        r: [f64; 2],
        side: [f64; 2],
    ) -> GnStation {
        let tr = Trunc::new(0, 0).unwrap();
        super::super::station::gn_solve_station(k, faces, leader, s, seed, r, side, tr, 1e-12)
            .unwrap()
    }

    fn contains(i: Ivl, x: f64) -> bool {
        i.lo <= x && x <= i.hi
    }

    /// A cubic in u with an interior knot of multiplicity `m` at u = 1 (knots 0 x4, 1 xm,
    /// 2 x4), linear in v, its control points' z not a polynomial across the knot.
    fn cubic_with_knot(m: usize) -> Nurbs {
        let mut u = vec![0.0; 4];
        u.extend(std::iter::repeat_n(1.0, m));
        u.extend([2.0; 4]);
        let nu = u.len() - 4;
        let mut pw = Vec::new();
        let zs = [0.0, 1.0, -0.5, 2.0, 0.3, -1.2];
        for (i, z) in (0..nu).zip(zs.iter().cycle()) {
            for y in [0.0, 1.0] {
                pw.extend([i as f64, y, z + 0.25 * y, 1.0]);
            }
        }
        Nurbs {
            p: 3,
            q: 1,
            nu,
            nv: 2,
            u,
            v: vec![0.0, 0.0, 1.0, 1.0],
            pw,
        }
    }

    /// A box across a knot line: the piecewise enclosure holds every partial of the face on
    /// both sides (the third u-partial jumps at the knot), which one piece continued over the
    /// box does not; a knot of too low a continuity, and a box past the knots, are refused.
    #[test]
    fn a_box_across_a_knot_line_is_enclosed_piece_by_piece() {
        let n = cubic_with_knot(1);
        let (ua, ub, va, vb) = (0.9, 1.1, 0.2, 0.3);
        let e = nurbs_partials_over_box(&n, ua, ub, va, vb, 3, FEET_BOX_CONTINUITY).unwrap();
        let left = nurbs_partials_over_piece(&n, 3, 1, ua, ub, va, vb, 3).unwrap();
        let mut continued_misses = false;
        for &u in &[0.9, 0.95, 1.0, 1.05, 1.1] {
            for &v in &[0.2, 0.25, 0.3] {
                // at the knot itself, either side's piece
                let spans = if u == 1.0 {
                    vec![3, 4]
                } else {
                    vec![span(&n.u, 3, n.nu, u)]
                };
                for ku in spans {
                    let x = nurbs_partials::<f64>(&n, ku, 1, u, v, 3);
                    for a in 0..=3 {
                        for b in 0..=(3 - a) {
                            let (p, q) = (x.at(a, b), e.at(a, b));
                            for (c, i) in [(p.x, q.x), (p.y, q.y), (p.z, q.z)] {
                                let slack = 1e-9 * (1.0 + c.abs());
                                assert!(
                                    i.lo - slack <= c && c <= i.hi + slack,
                                    "P[{a}][{b}] at ({u}, {v}): {c} not in {i:?}"
                                );
                            }
                        }
                    }
                    let c = x.at(3, 0).z;
                    continued_misses |= !contains(left.at(3, 0).z, c);
                }
            }
        }
        assert!(
            continued_misses,
            "the left piece continued is not the face past the knot"
        );
        // a point on the knot: both sides' partials (the third u-partial jumps there)
        let at = nurbs_partials_over_box(&n, 1.0, 1.0, 0.25, 0.25, 3, FEET_BOX_CONTINUITY).unwrap();
        let sides = [3, 4].map(|ku| nurbs_partials::<f64>(&n, ku, 1, 1.0, 0.25, 3).at(3, 0).z);
        assert!(sides[0] != sides[1], "{sides:?}");
        for c in sides {
            assert!(contains(at.at(3, 0).z, c), "{c} not in {:?}", at.at(3, 0).z);
        }
        // a double knot (C^1) crossed by a feet's box is refused; at C^1 asked, enclosed
        let d = cubic_with_knot(2);
        assert!(matches!(
            nurbs_partials_over_box(&d, ua, ub, va, vb, 3, FEET_BOX_CONTINUITY),
            Err(SolidError::Unsupported { .. })
        ));
        assert!(nurbs_partials_over_box(&d, ua, ub, va, vb, 3, 1).is_ok());
        // past the knots
        assert!(matches!(
            nurbs_partials_over_box(&n, 1.9, 2.1, va, vb, 2, FEET_BOX_CONTINUITY),
            Err(SolidError::Unsupported { .. })
        ));
    }

    /// Two planes at a square corner, the ball of radius r on the edge along y (the
    /// station test's closed form): the feet (tau r, L s, 0) and (0, L s, tau r) are linear
    /// in s and tau, so every fourth s-derivative is 0 and the span's bound is 0 to the
    /// enclosures' rounding. Checked against that closed form: the enclosed series hold
    /// the exact coefficients (`du1/deps = r`, `dv1/ds = L`, the box holds the feet over
    /// the span), every `m4`, `uv4` is below 1e-9, and the certified bound is in
    /// `[0, 1e-9]`.
    #[test]
    fn two_plane_span_bound_is_zero() {
        let k = Kernel::default();
        let f1 = plane(Vec3::ZERO, Vec3::X, Vec3::Y); // z = 0, (u, v) = (x, y)
        let f2 = plane(Vec3::ZERO, Vec3::Y, Vec3::Z); // x = 0, (u, v) = (y, z)
        let (r, len) = (2.0, 5.0);
        let leader =
            BSplineCurve::create_line_segment(3, Vec3::new(r, 0.0, r), Vec3::new(r, len, r))
                .unwrap();
        let ln = NurbsCurve::of(&leader).unwrap();
        assert_eq!((ln.p, ln.n, ln.u.len()), (1, 2, 4));
        let (n1, n2) = (nurbs(&k, &f1), nurbs(&k, &f2));
        let faces = [&f1, &f2];
        let (rr, side) = ([r, r], [1.0, 1.0]);
        let (sa, sb) = (0.2, 0.6);
        let a = station(&k, faces, &leader, sa, [r, len * sa, len * sa, r], rr, side);
        let b = station(&k, faces, &leader, sb, [r, len * sb, len * sb, r], rr, side);
        let n = 2;
        let e = gn_enclose(
            &k,
            [&n1, &n2],
            &ln,
            &a,
            &b,
            rr,
            side,
            n,
            4 + TAYLOR,
            None,
            None,
        )
        .unwrap()
        .expect("proved");
        assert!(contains(e.x[0].c[0][0], r) && contains(e.x[3].c[0][0], r));
        for q in [1, 2] {
            let x = e.x[q].c[0][0];
            assert!(x.lo <= len * sa && len * sb <= x.hi, "{q} {x:?}");
            assert!(contains(e.x[q].c[1][0], len), "{q}");
        }
        assert!(contains(e.x[0].c[0][1], r) && contains(e.x[3].c[0][1], r));
        // |S_u| and |S_v| of a unit plane are 1
        for dsf in e.ds {
            assert!(dsf.iter().all(|d| (d - 1.0).abs() < 1e-12), "{dsf:?}");
        }
        let (m4, uv4, ..) =
            gn_fourth_over(&k, [&n1, &n2], &ln, faces, &leader, &a, &b, rr, side, n).unwrap();
        for f in 0..2 {
            assert!(uv4[f] >= 0.0 && uv4[f] < 1e-9, "{uv4:?}");
            for kk in 0..=n {
                assert!(m4[f][kk] >= 0.0 && m4[f][kk] < 1e-9, "{m4:?}");
            }
        }
        let bound = gn_certify_span(
            &k,
            [&n1, &n2],
            &ln,
            faces,
            &leader,
            &a,
            &b,
            rr,
            side,
            n,
            1e-6,
            None,
            None,
            None,
        )
        .unwrap();
        assert!((0.0..1e-9).contains(&bound), "{bound}");
    }

    /// The (u, v) of the face's net nearest to `p`, on a grid (a seed for Newton).
    fn nearest_uv(n: &Nurbs, p: Vec3) -> (f64, f64) {
        let (u0, u1) = (n.u[0], n.u[n.u.len() - 1]);
        let (v0, v1) = (n.v[0], n.v[n.v.len() - 1]);
        let mut best = (f64::INFINITY, u0, v0);
        for i in 0..=400 {
            for j in 0..=40 {
                let u = u0 + (u1 - u0) * f64::from(i) / 400.0;
                let v = v0 + (v1 - v0) * f64::from(j) / 40.0;
                let ku = span(&n.u, n.p, n.nu, u);
                let kv = span(&n.v, n.q, n.nv, v);
                let d = (vec3(nurbs_partials(n, ku, kv, u, v, 0).at(0, 0)) - p).length();
                if d < best.0 {
                    best = (d, u, v);
                }
            }
        }
        (best.1, best.2)
    }

    /// A ball of radius r inside a cylinder of radius R about z, on the plane z = 0, r near
    /// R (#2690): the cylinder's offset by r, of radius R - r = R / 100, nearly collapses;
    /// the feet's system `S + r N` is a sum of terms a hundred times its size and its
    /// Jacobian a hundredth of theirs. The centre line is the circle of radius R - r at
    /// height r and the plane's foot is its point below, so the plane's (u, v) are the
    /// circle's, `L_xy`, exactly. Krawczyk's test is proved over spans of 5e-3, 1e-3 and
    /// 2e-4 of the leader (in plain intervals, over a box, only from 2e-4 down), every
    /// enclosed Taylor coefficient of the plane's foot holds the exact one (`L^(i) / i!`,
    /// sampled) over the span, and the fourth's excess over that range falls at least
    /// twentyfold from 1e-3 to 2e-4, as the square of the span (the first-order excess of a
    /// box falls fivefold).
    #[test]
    fn a_nearly_collapsed_offset_s_feet_are_enclosed_to_their_range() {
        let k = Kernel::default();
        let (big_r, r) = (3.0, 2.97);
        let fp = plane(Vec3::ZERO, Vec3::X, Vec3::Y);
        let mut cyl = Cone::new_cylinder_bounded(
            Vec3::new(0.0, 0.0, -2.0),
            Vec3::X,
            Vec3::Y,
            big_r,
            0.0,
            360.0,
            10.0,
            false,
            false,
            false,
        );
        cyl.make_nurb(&k).unwrap();
        let fc = Surface::from(cyl);
        let frame = Axis2Placement::new(Vec3::new(0.0, 0.0, r), Vec3::X, Vec3::Y);
        let leader = BSplineCurve::create_circle_segment(
            3,
            &frame,
            big_r - r,
            10.0,
            80.0,
            NurbCircleParam::Quadratic,
        )
        .unwrap();
        let ln = NurbsCurve::of(&leader).unwrap();
        let (np, nc) = (nurbs(&k, &fp), nurbs(&k, &fc));
        let faces = [&fp, &fc];
        let lpt = |s: f64| leader.evaluate_point(s).unwrap();
        let foot_c = |s: f64| {
            let l = lpt(s);
            let q = big_r / (big_r - r);
            Vec3::new(l.x * q, l.y * q, r)
        };
        let (cu, cv) = nearest_uv(&nc, foot_c(0.45));
        let pc = nurbs_partials(
            &nc,
            span(&nc.u, nc.p, nc.nu, cu),
            span(&nc.v, nc.q, nc.nv, cv),
            cu,
            cv,
            1,
        );
        // the ball is inside: the offset runs against the outward direction
        let nrm = vec3(pc.at(1, 0)).cross(vec3(pc.at(0, 1)));
        let out = foot_c(0.45) - Vec3::new(0.0, 0.0, r);
        let side_c = if nrm.dot(out) > 0.0 { -1.0 } else { 1.0 };
        let (rr, side) = ([r, r], [1.0, side_c]);
        let st = |s: f64| {
            let l = lpt(s);
            let (u, v) = nearest_uv(&nc, foot_c(s));
            station(&k, faces, &leader, s, [l.x, l.y, u, v], rr, side)
        };
        let sa = 0.4;
        let a = st(sa);
        assert!(
            (vec3(V3 {
                x: a.f[1].a[0].c[0][0],
                y: a.f[1].a[1].c[0][0],
                z: a.f[1].a[2].c[0][0],
            }) - foot_c(sa))
            .length()
                < 1e-9
        );
        let mut excess = Vec::new();
        for h in [5e-3, 1e-3, 2e-4] {
            let sb = sa + h;
            let e = gn_enclose(
                &k,
                [&np, &nc],
                &ln,
                &a,
                &st(sb),
                rr,
                side,
                2,
                4 + TAYLOR,
                None,
                None,
            )
            .unwrap()
            .expect("Krawczyk's test is proved over the span");
            // the exact coefficients' ranges over the span: L^(i) / i!, sampled
            let mut fact = 1.0_f64;
            let mut e4 = 0.0_f64;
            for i in 0..=4 {
                if i > 0 {
                    fact *= i as f64;
                }
                for (c, x) in [(0, e.x[0].c[i][0]), (1, e.x[1].c[i][0])] {
                    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
                    for j in 0..=200 {
                        let sj = sa + h * f64::from(j) / 200.0;
                        let d = leader.evaluate(sj, i, true, false).unwrap()[i];
                        let v = [d.x, d.y][c] / fact;
                        lo = lo.min(v);
                        hi = hi.max(v);
                    }
                    let slack = 1e-9 * (1.0 + lo.abs().max(hi.abs()));
                    assert!(x.lo <= lo + slack && hi - slack <= x.hi, "c{i} {c}: {x:?}");
                    if i == 4 {
                        e4 = e4.max((x.hi - x.lo) - (hi - lo));
                    }
                }
            }
            excess.push(e4);
        }
        assert!(
            excess[2] <= excess[1] / 20.0,
            "the fourth coefficient's excess: {excess:?}"
        );
    }

    /// A cylinder of radius R about z standing on the plane z = 0, the ball of radius r
    /// outside it: the centre line L(s) is the circle of radius R + r at height r (a
    /// rational quadratic arc), and the feet are `F1 = (L_xy (R + tau r) / (R + r), 0)` on
    /// the plane and `F2 = (L_xy R / (R + r), tau r)` on the cylinder. So the jets' fourth
    /// s-derivatives are, exactly, `|c10''''| = |L''''|`, `|c11''''| = r / (R + r) |L''''|`,
    /// `|c20''''| = R / (R + r) |L''''|` and 0 for the rest, and the plane's uv curve is
    /// `L_xy` (`uv4 >= |L_x''''| + |L_y''''|`, the plane's partials of length 1). Checked
    /// against that closed form, with
    /// L'''' from the leader's own evaluator: at the span's middle (one point, a tight
    /// enclosure) `4! c[4][k]` of every jet and of the plane's (u, v) holds the closed
    /// form's value; over the span every `m4`, `uv4` is finite and at least the closed
    /// form's maximum (sampled at 201 points); the certified bound is finite, below the
    /// target, at least `h^4 / 384` times the closed-form jets' sum, and level 0's (the
    /// pieces' maximum). The span's bounds are sound but loose (`m4[0][0]` about 36 times
    /// the closed form at h = 0.02: the series enclosed over an interval of s widen order
    /// by order, and the Taylor form's remainder carries that). A span five times as long
    /// fails Krawczyk's test in one enclosure and is proved by the pieces, its bound again
    /// finite, below the target and above the closed form's.
    /// A rational face's partials over a box are enclosed as tightly wherever the face
    /// lies (#2690): the cylinder r 3 at the origin and the same cylinder moved 1000 away,
    /// over the same box of their nets, have enclosures of the same widths, which hold the
    /// partials' sampled values and shrink with the box. The quotient `P = A / w` taken in
    /// intervals on the homogeneous coordinates as they stand widens with the point's
    /// distance from the origin, `|P| |dw| / w`, however small the box.
    #[test]
    fn a_rational_face_s_box_enclosure_does_not_widen_with_its_distance() {
        let k = Kernel::default();
        let cylinder = |at: Vec3| {
            let mut c = Cone::new_cylinder_bounded(
                at,
                Vec3::X,
                Vec3::Y,
                3.0,
                0.0,
                360.0,
                10.0,
                false,
                false,
                false,
            );
            c.make_nurb(&k).unwrap();
            nurbs(&k, &Surface::from(c))
        };
        let near = cylinder(Vec3::new(0.0, 0.0, -2.0));
        let far = cylinder(Vec3::new(1000.0, 500.0, -2.0));
        let lk = 3;
        // a box inside one piece of the net, off its middle
        let ku = span(
            &near.u,
            near.p,
            near.nu,
            0.5 * (near.u[near.p] + near.u[near.p + 1]),
        );
        let kv = span(
            &near.v,
            near.q,
            near.nv,
            0.5 * (near.v[near.q] + near.v[near.q + 1]),
        );
        let (u0, u1) = (near.u[ku], near.u[ku + 1]);
        let (v0, v1) = (near.v[kv], near.v[kv + 1]);
        let (um, vm) = (u0 + 0.37 * (u1 - u0), v0 + 0.61 * (v1 - v0));
        let mut widths = Vec::new();
        for delta in [1e-2, 1e-3] {
            let (du, dv) = (delta * (u1 - u0), delta * (v1 - v0));
            let (ua, ub, va, vb) = (um - du, um + du, vm - dv, vm + dv);
            let en = nurbs_partials_over_piece(&near, ku, kv, ua, ub, va, vb, lk).unwrap();
            let ef = nurbs_partials_over_piece(&far, ku, kv, ua, ub, va, vb, lk).unwrap();
            for a in 0..=lk {
                for b in 0..=(lk - a) {
                    // the partial's range over the box, sampled on a grid
                    let mut lo = [f64::INFINITY; 3];
                    let mut hi = [f64::NEG_INFINITY; 3];
                    for i in 0..=8 {
                        for j in 0..=8 {
                            let u = ua + (ub - ua) * f64::from(i) / 8.0;
                            let v = va + (vb - va) * f64::from(j) / 8.0;
                            let x = nurbs_partials::<f64>(&near, ku, kv, u, v, lk).at(a, b);
                            for (c, xc) in [x.x, x.y, x.z].into_iter().enumerate() {
                                lo[c] = lo[c].min(xc);
                                hi[c] = hi[c].max(xc);
                            }
                        }
                    }
                    let (pn, pf) = (en.at(a, b), ef.at(a, b));
                    let comps = [(pn.x, pf.x), (pn.y, pf.y), (pn.z, pf.z)];
                    for (c, (n, f)) in comps.into_iter().enumerate() {
                        let (wn, wf) = (n.hi - n.lo, f.hi - f.lo);
                        // rounding: the far face's coordinates are 1000 times as large
                        let floor = 1e-9 * (1.0 + hi[c].abs().max(lo[c].abs()));
                        assert!(
                            wf <= 2.0 * wn + floor,
                            "P[{a}][{b}].{c} at {delta}: far {wf:e}, near {wn:e}"
                        );
                        assert!(n.lo <= lo[c] && hi[c] <= n.hi, "P[{a}][{b}].{c}");
                        widths.push(wn);
                    }
                }
            }
        }
        // a tenth of the box, at most a fifth of the width (first order, and rounding)
        let (wide, narrow) = widths.split_at(widths.len() / 2);
        for (w1, w2) in wide.iter().zip(narrow) {
            assert!(*w2 <= 0.2 * w1 + 1e-12, "{w2:e} against {w1:e}");
        }
        // the same for a rational curve: the circle r 3, at the origin and 1000 away
        let circle = |at: Vec3| {
            let frame = Axis2Placement::new(at, Vec3::X, Vec3::Y);
            let c = BSplineCurve::create_circle_segment(
                3,
                &frame,
                3.0,
                10.0,
                80.0,
                NurbCircleParam::Quadratic,
            )
            .unwrap();
            NurbsCurve::of(&c).unwrap()
        };
        let (cn, cf) = (circle(Vec3::ZERO), circle(Vec3::new(1000.0, 500.0, 0.0)));
        let piece = span(&cn.u, cn.p, cn.n, 0.5 * (cn.u[cn.p] + cn.u[cn.p + 1]));
        let (s0, s1) = (cn.u[piece], cn.u[piece + 1]);
        let sm = s0 + 0.37 * (s1 - s0);
        for delta in [1e-2, 1e-3] {
            let ds = delta * (s1 - s0);
            let s = Ivl::new(sm - ds, sm + ds);
            let dn = curve_derivs_over(&cn, s, lk, Some(piece)).unwrap();
            let df = curve_derivs_over(&cf, s, lk, Some(piece)).unwrap();
            for (d, (n, f)) in dn.iter().zip(&df).enumerate() {
                for (c, (n, f)) in [(n.x, f.x), (n.y, f.y), (n.z, f.z)].into_iter().enumerate() {
                    let (wn, wf) = (n.hi - n.lo, f.hi - f.lo);
                    let floor = 1e-9 * (1.0 + f.lo.abs().max(f.hi.abs()));
                    assert!(
                        wf <= 2.0 * wn + floor,
                        "D[{d}].{c} at {delta}: far {wf:e}, near {wn:e}"
                    );
                }
            }
        }
    }

    /// A rational face's partials, and a rational curve's derivatives, over a box are
    /// enclosed to within the box's square of their ranges (#2690): from a box of 1e-2 of a
    /// piece to one of 1e-3, each enclosure's excess over its sampled range falls at least
    /// thirtyfold (the square's hundredfold, less the sampling's). The direct quotient's
    /// excess, the partial's size times w's spread, falls only tenfold: the cylinder's
    /// partial along its axis is constant, its direct enclosure 7e-3 wide at 1e-2.
    #[test]
    fn a_rational_box_enclosure_exceeds_its_range_by_the_box_s_square() {
        let k = Kernel::default();
        let mut c = Cone::new_cylinder_bounded(
            Vec3::new(0.0, 0.0, -2.0),
            Vec3::X,
            Vec3::Y,
            3.0,
            0.0,
            360.0,
            10.0,
            false,
            false,
            false,
        );
        c.make_nurb(&k).unwrap();
        let n = nurbs(&k, &Surface::from(c));
        let lk = 3;
        let ku = span(&n.u, n.p, n.nu, 0.5 * (n.u[n.p] + n.u[n.p + 1]));
        let kv = span(&n.v, n.q, n.nv, 0.5 * (n.v[n.q] + n.v[n.q + 1]));
        let (u0, u1) = (n.u[ku], n.u[ku + 1]);
        let (v0, v1) = (n.v[kv], n.v[kv + 1]);
        let (um, vm) = (u0 + 0.37 * (u1 - u0), v0 + 0.61 * (v1 - v0));
        let mut excess = Vec::new();
        for delta in [1e-2, 1e-3] {
            let (du, dv) = (delta * (u1 - u0), delta * (v1 - v0));
            let (ua, ub, va, vb) = (um - du, um + du, vm - dv, vm + dv);
            let e = nurbs_partials_over_piece(&n, ku, kv, ua, ub, va, vb, lk).unwrap();
            for a in 0..=lk {
                for b in 0..=(lk - a) {
                    let mut lo = [f64::INFINITY; 3];
                    let mut hi = [f64::NEG_INFINITY; 3];
                    for i in 0..=16 {
                        for j in 0..=16 {
                            let u = ua + (ub - ua) * f64::from(i) / 16.0;
                            let v = va + (vb - va) * f64::from(j) / 16.0;
                            let x = nurbs_partials::<f64>(&n, ku, kv, u, v, lk).at(a, b);
                            for (c, xc) in [x.x, x.y, x.z].into_iter().enumerate() {
                                lo[c] = lo[c].min(xc);
                                hi[c] = hi[c].max(xc);
                            }
                        }
                    }
                    let p = e.at(a, b);
                    for (c, i) in [p.x, p.y, p.z].into_iter().enumerate() {
                        assert!(i.lo <= lo[c] && hi[c] <= i.hi, "P[{a}][{b}].{c}");
                        excess.push(((a, b, c), (i.hi - i.lo) - (hi[c] - lo[c])));
                    }
                }
            }
        }
        let (wide, narrow) = excess.split_at(excess.len() / 2);
        for ((at, e1), (_, e2)) in wide.iter().zip(narrow) {
            assert!(*e2 <= e1 / 30.0 + 1e-12, "P{at:?}: {e2:e} against {e1:e}");
        }
        // the circle r 3 likewise, its derivatives over an interval of the parameter
        let frame = Axis2Placement::new(Vec3::ZERO, Vec3::X, Vec3::Y);
        let circle = BSplineCurve::create_circle_segment(
            3,
            &frame,
            3.0,
            10.0,
            80.0,
            NurbCircleParam::Quadratic,
        )
        .unwrap();
        let cn = NurbsCurve::of(&circle).unwrap();
        let piece = span(&cn.u, cn.p, cn.n, 0.5 * (cn.u[cn.p] + cn.u[cn.p + 1]));
        let (s0, s1) = (cn.u[piece], cn.u[piece + 1]);
        let sm = s0 + 0.37 * (s1 - s0);
        let mut excess = Vec::new();
        for delta in [1e-2, 1e-3] {
            let ds = delta * (s1 - s0);
            let d = curve_derivs_over(&cn, Ivl::new(sm - ds, sm + ds), lk, Some(piece)).unwrap();
            for (o, dd) in d.iter().enumerate() {
                let mut lo = [f64::INFINITY; 3];
                let mut hi = [f64::NEG_INFINITY; 3];
                for i in 0..=64 {
                    let s = sm - ds + 2.0 * ds * f64::from(i) / 64.0;
                    let x = circle.evaluate(s, lk, true, false).unwrap()[o];
                    for (c, xc) in [x.x, x.y, x.z].into_iter().enumerate() {
                        lo[c] = lo[c].min(xc);
                        hi[c] = hi[c].max(xc);
                    }
                }
                for (c, i) in [dd.x, dd.y, dd.z].into_iter().enumerate() {
                    let slack = 1e-9 * (1.0 + lo[c].abs().max(hi[c].abs()));
                    assert!(i.lo <= lo[c] + slack && hi[c] - slack <= i.hi, "D[{o}].{c}");
                    excess.push(((o, c), (i.hi - i.lo) - (hi[c] - lo[c])));
                }
            }
        }
        let (wide, narrow) = excess.split_at(excess.len() / 2);
        for ((at, e1), (_, e2)) in wide.iter().zip(narrow) {
            assert!(*e2 <= e1 / 30.0 + 1e-9, "D{at:?}: {e2:e} against {e1:e}");
        }
    }

    #[test]
    fn cylinder_plane_span_bound_is_the_closed_form() {
        let k = Kernel::default();
        let (big_r, r) = (3.0, 1.0);
        let fp = plane(Vec3::ZERO, Vec3::X, Vec3::Y);
        let mut cyl = Cone::new_cylinder_bounded(
            Vec3::new(0.0, 0.0, -2.0),
            Vec3::X,
            Vec3::Y,
            big_r,
            0.0,
            360.0,
            10.0,
            false,
            false,
            false,
        );
        cyl.make_nurb(&k).unwrap();
        let fc = Surface::from(cyl);
        let frame = Axis2Placement::new(Vec3::new(0.0, 0.0, r), Vec3::X, Vec3::Y);
        let leader = BSplineCurve::create_circle_segment(
            3,
            &frame,
            big_r + r,
            10.0,
            80.0,
            NurbCircleParam::Quadratic,
        )
        .unwrap();
        let ln = NurbsCurve::of(&leader).unwrap();
        let (np, nc) = (nurbs(&k, &fp), nurbs(&k, &fc));
        let faces = [&fp, &fc];
        // the sides: the plane's normal is +z, towards the ball; the cylinder's by its net
        let lpt = |s: f64| leader.evaluate_point(s).unwrap();
        let foot_c = |s: f64| {
            let l = lpt(s);
            let q = big_r / (big_r + r);
            Vec3::new(l.x * q, l.y * q, r)
        };
        let (cu, cv) = nearest_uv(&nc, foot_c(0.45));
        let pc = nurbs_partials(
            &nc,
            span(&nc.u, nc.p, nc.nu, cu),
            span(&nc.v, nc.q, nc.nv, cv),
            cu,
            cv,
            1,
        );
        let nrm = vec3(pc.at(1, 0)).cross(vec3(pc.at(0, 1)));
        let out = foot_c(0.45) - Vec3::new(0.0, 0.0, r);
        let side_c = if nrm.dot(out) > 0.0 { 1.0 } else { -1.0 };
        let (rr, side) = ([r, r], [1.0, side_c]);
        let st = |s: f64| {
            let l = lpt(s);
            let (u, v) = nearest_uv(&nc, foot_c(s));
            station(&k, faces, &leader, s, [l.x, l.y, u, v], rr, side)
        };
        let (sa, sb) = (0.4, 0.42);
        let (a, b) = (st(sa), st(sb));
        // the stations are the closed form's
        let fa = vec3(V3 {
            x: a.f[1].a[0].c[0][0],
            y: a.f[1].a[1].c[0][0],
            z: a.f[1].a[2].c[0][0],
        });
        assert!((fa - foot_c(sa)).length() < 1e-9);
        assert!((a.x[0].c[0][0] - lpt(sa).x).abs() < 1e-9);

        // L'''' over the span, sampled
        let (mut l4, mut l4x, mut l4y) = (0.0_f64, 0.0_f64, 0.0_f64);
        for i in 0..=200 {
            let s = sa + (sb - sa) * f64::from(i) / 200.0;
            let d = leader.evaluate(s, 4, true, false).unwrap();
            let d4 = d[4];
            assert!(d4.z.abs() < 1e-9);
            l4 = l4.max(Vec3::new(d4.x, d4.y, 0.0).length());
            l4x = l4x.max(d4.x.abs());
            l4y = l4y.max(d4.y.abs());
        }
        assert!(l4 > 1.0, "{l4}");
        let n = 2;
        // at the span's middle the jets' fourth derivatives are enclosed on one point: they
        // hold the closed form's, `(4! c[4][k])`
        let sm = 0.5 * (sa + sb);
        let ms = st(sm);
        let em = gn_enclose(
            &k,
            [&np, &nc],
            &ln,
            &ms,
            &ms,
            rr,
            side,
            n,
            4 + TAYLOR - 1,
            None,
            None,
        )
        .unwrap()
        .expect("proved at a point");
        let dm = leader.evaluate(sm, 4, true, false).unwrap()[4];
        let holds = |i: Ivl, x: f64| {
            let j = Ivl::from(24.0) * i;
            j.lo - 1e-9 <= x && x <= j.hi + 1e-9 && j.hi - j.lo < 1e-2
        };
        let (qp, qc) = (r / (big_r + r), big_r / (big_r + r));
        for (c, l) in [(0, dm.x), (1, dm.y), (2, 0.0)] {
            assert!(holds(em.f[0].a[c].c[4][0], l), "c10 {c}");
            assert!(holds(em.f[0].a[c].c[4][1], qp * l), "c11 {c}");
            assert!(holds(em.f[1].a[c].c[4][0], qc * l), "c20 {c}");
            assert!(holds(em.f[1].a[c].c[4][1], 0.0), "c21 {c}");
            assert!(holds(em.f[0].a[c].c[4][2], 0.0) && holds(em.f[1].a[c].c[4][2], 0.0));
        }
        // the plane's (u, v) is the foot's (x, y)
        assert!(holds(em.x[0].c[4][0], dm.x) && holds(em.x[1].c[4][0], dm.y));
        let (m4, uv4, ..) =
            gn_fourth_over(&k, [&np, &nc], &ln, faces, &leader, &a, &b, rr, side, n).unwrap();
        let want = [
            [l4, r / (big_r + r) * l4, 0.0],
            [big_r / (big_r + r) * l4, 0.0, 0.0],
        ];
        for f in 0..2 {
            for kk in 0..=n {
                let (got, w) = (m4[f][kk], want[f][kk]);
                assert!(got.is_finite() && w <= got, "{f} {kk} {got} {w}");
            }
        }
        let wuv = l4x + l4y;
        assert!(wuv <= uv4[0] && uv4[0].is_finite(), "{} {wuv}", uv4[0]);
        assert!(uv4[1].is_finite());

        let h = sb - sa;
        let target = 1e-3;
        let bound = gn_certify_span(
            &k,
            [&np, &nc],
            &ln,
            faces,
            &leader,
            &a,
            &b,
            rr,
            side,
            n,
            target,
            None,
            None,
            None,
        )
        .unwrap();
        let jets = h * h * h * h / 384.0 * (want[0][0] + want[0][1]);
        assert!(bound.is_finite() && bound < target, "{bound}");
        assert!(
            bound >= jets && bound >= h * h * h * h / 384.0 * wuv,
            "{bound} {jets}"
        );
        // below the target at level 0: the bound is that level's, the pieces' maximum
        let lvl0 = h * h * h * h / 384.0
            * (m4[0][0] + m4[0][1] + m4[0][2])
                .max(m4[1][0] + m4[1][1] + m4[1][2])
                .max(uv4[0])
                .max(uv4[1]);
        assert!((bound - lvl0).abs() <= 1e-9 * lvl0, "{bound} {lvl0}");

        // a span five times as long: one enclosure proves Krawczyk's test over it (the
        // system in affine forms, #2690; in plain intervals it failed and only the pieces
        // proved it), and the bound still holds the closed form's
        let (lb, hl) = (st(0.5), 0.1);
        let mut l4l = 0.0_f64;
        for i in 0..=500 {
            let d = leader
                .evaluate(sa + hl * f64::from(i) / 500.0, 4, true, false)
                .unwrap();
            l4l = l4l.max(Vec3::new(d[4].x, d[4].y, 0.0).length());
        }
        let one = gn_enclose(
            &k,
            [&np, &nc],
            &ln,
            &a,
            &lb,
            rr,
            side,
            n,
            4 + TAYLOR,
            None,
            None,
        )
        .unwrap();
        assert!(one.is_some());
        let long = gn_certify_span(
            &k,
            [&np, &nc],
            &ln,
            faces,
            &leader,
            &a,
            &lb,
            rr,
            side,
            n,
            target,
            None,
            None,
            None,
        )
        .unwrap();
        let jl = hl * hl * hl * hl / 384.0 * (1.0 + r / (big_r + r)) * l4l;
        assert!(
            long.is_finite() && long < target && long >= jl,
            "{long} {jl}"
        );
    }
}
