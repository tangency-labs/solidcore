//! A blend by its definition re-certified (#2388): its net B proved within the
//! definition's bound ε of S again, as a blend rebuilt from a record must be before it is
//! an operand. A net that is not proved within ε, or is not a net of its definition's
//! stations, is refused.
//!
//! Not upstream. The fork's `SmBallPatch::IsOnSphere` is the same check for the Ball's
//! sphere patch, which an exact algebraic identity settles; a rolling-ball round has
//! none, so the round's own certificate is run again.
//!
//! The net is the construction's (`Build`): along s, on each span between two stations, the
//! cubic Hermite interpolant of the section's rows at its stations; across, the section's
//! exact rational form. The stations are the seeds' knots (the seeds are the feet's
//! Hermite curves on them, and are never cut). So for each piece of B:
//!
//! - its span's stations are solved again from the definition, and the span's Hermite
//!   rows H built from them ([`NetForm::span`], the construction's own);
//! - `|B - H|` is bounded on the piece by its control nets: with both rational Bezier
//!   patches of the same degrees, Bernstein bases `N`, points `P_i`, `Q_i` and positive
//!   weights `w_i`, `v_i`, `B - H = sum_ij N_i N_j [w_i v_j (P_i - Q_i) + (Q_i - c)
//!   (w_i v_j - v_i w_j)] / (W V)` for any point c, so
//!   `|B - H| <= max |P_i - Q_i| + max |Q_i - c| max |w_i v_j - v_i w_j| / (min w min v)`;
//! - `|H - S|` on the span is the round's own certificate ([`gn_certify_span`]), asked for
//!   what ε leaves after that.
//!
//! The bound is their sum, at most ε, or the net is refused. `|B - H|` is taken as at
//! least the net's precision ([`net_precision`]): S is computed to it, so differences
//! below it are not the net's, and every bound a blend stores includes it. A net cut or reparametrized
//! after it was built (`SplitAt`, `Reparameterize`, `Reverse`, `SwapUV`) is a piece of the
//! same Hermite net: H is cut to the piece, in homogeneous coordinates, before the nets
//! are compared.
//!
//! Deviations: none (no upstream counterpart). No recursion.

use core::cmp::Ordering;

use super::certificate::{gn_certify_span, GnLaw, NurbsCurve};
use super::section::{
    certificate_faces, crosses_knot, inner_knots, param, varying, NetForm, StationSystem,
    GN_THUMBWEIGHT_EQUAL_AREA,
};
use super::series::Trunc;
use super::station::GnStation;
use crate::error::{Result, SolidError};
use crate::kernel::Kernel;
use crate::math::Vec3;
use crate::nurbs::b_spline_surface::BSplineSurface;
use crate::nurbs::math::NOW;
use crate::nurbs::srf_tool::p1::decompose_to_bez;
use crate::surface::procedural_blend::ProceduralBlend;

const OP: &str = "blend::fillet_section::recertify";

/// A homogeneous control point `(w P, w)`.
type Hpt = (Vec3, f64);

/// A Bezier patch over the definition's (s, t): its homogeneous control points
/// `[i along s][j along t]` and its box.
struct Piece {
    net: Vec<Vec<Hpt>>,
    s: [f64; 2],
    t: [f64; 2],
}

/// The largest rounding of `a + b x` (its two operations, each within half an ulp, and
/// the operands' size).
fn map_rounding(a: f64, b: f64, x: f64) -> f64 {
    4.0 * f64::EPSILON * (a.abs() + (b * x).abs())
}

/// One Bezier curve of homogeneous points cut to `[lo, hi]` of its `[0, 1]` (de Casteljau,
/// at hi, then the left part at lo / hi).
fn sub_curve(pts: &[Hpt], lo: f64, hi: f64) -> Vec<Hpt> {
    let split = |p: &[Hpt], tau: f64| -> (Vec<Hpt>, Vec<Hpt>) {
        let mut level: Vec<Hpt> = p.to_vec();
        let mut left = Vec::new();
        let mut right = Vec::new();
        while let (Some(&first), Some(&last)) = (level.first(), level.last()) {
            left.push(first);
            right.push(last);
            level = level
                .windows(2)
                .map(|w| match w {
                    [(pa, wa), (pb, wb)] => {
                        (*pa * (1.0 - tau) + *pb * tau, wa * (1.0 - tau) + wb * tau)
                    }
                    _ => (Vec3::ZERO, 0.0),
                })
                .collect();
        }
        right.reverse();
        (left, right)
    };
    let (left, _) = if hi < 1.0 {
        split(pts, hi)
    } else {
        (pts.to_vec(), Vec::new())
    };
    if lo > 0.0 && hi > 0.0 {
        split(&left, lo / hi).1
    } else {
        left
    }
}

/// `|B - H|` over a piece, from the two patches' control nets (module docs).
///
/// # Errors
///
/// `InvalidInput` where the nets differ in size or a weight is not positive.
pub(crate) fn net_distance(b: &[Vec<Hpt>], h: &[Vec<Hpt>]) -> Result<f64> {
    let flat = |n: &[Vec<Hpt>]| -> Vec<Hpt> { n.iter().flatten().copied().collect() };
    let (pb, ph) = (flat(b), flat(h));
    if pb.len() != ph.len() || pb.is_empty() {
        return Err(SolidError::invalid_input(
            OP,
            "the net is not of its definition's degrees",
        ));
    }
    let positive = |p: &[Hpt]| p.iter().all(|&(_, w)| w > 0.0 && w.is_finite());
    if !positive(&pb) || !positive(&ph) {
        return Err(SolidError::invalid_input(
            OP,
            "a weight of the net is not positive",
        ));
    }
    let euclid = |(p, w): Hpt| p * (1.0 / w);
    let c = ph.first().map_or(Vec3::ZERO, |&q| euclid(q));
    let mut d_points = 0.0_f64;
    let mut reach = 0.0_f64;
    for (&x, &y) in pb.iter().zip(&ph) {
        d_points = d_points.max((euclid(x) - euclid(y)).length());
        reach = reach.max((euclid(y) - c).length());
    }
    let mut cross = 0.0_f64;
    for (i, &(_, wi)) in pb.iter().enumerate() {
        let vi = ph.get(i).map_or(f64::NAN, |q| q.1);
        for (j, &(_, wj)) in pb.iter().enumerate() {
            let vj = ph.get(j).map_or(f64::NAN, |q| q.1);
            cross = cross.max((wi * vj - vi * wj).abs());
        }
    }
    let w_min = pb.iter().map(|q| q.1).fold(f64::INFINITY, f64::min);
    let v_min = ph.iter().map(|q| q.1).fold(f64::INFINITY, f64::min);
    let d = d_points + reach * cross / (w_min * v_min);
    if d.is_finite() {
        Ok(d)
    } else {
        Err(SolidError::invalid_input(OP, "the net is not finite"))
    }
}

/// The precision of a net and of S where it lies: 16 units in the last place of its
/// largest coordinate. Below it two computations of S (its stations solved from other
/// starts, its rows summed, the net divided out of them) need not agree, so no bound is
/// certified finer: a blend's bound is at least its certificate and this, and the
/// comparison of a net with its rows takes differences up to it as this.
///
/// # Errors
///
/// `SM_ERR` without NURB data.
pub(crate) fn net_precision(net: &BSplineSurface) -> Result<f64> {
    let (_, _, stored) = net.control_points()?;
    let mut size = 0.0_f64;
    for c in &stored {
        let w = if c.w == NOW { 1.0 } else { c.w };
        size = size.max(Vec3::new(c.x, c.y, c.z).length() / w);
    }
    Ok(16.0 * f64::EPSILON * size)
}

/// A span's Hermite rows ([`NetForm::span`]) as a grid `[i along s][j across]`.
fn span_grid(homog: &[Vec3], w: &[f64], ld: usize) -> Vec<Vec<Hpt>> {
    (0..4)
        .map(|i| {
            (0..ld)
                .map(|j| {
                    let at = i * ld + j;
                    (
                        homog.get(at).copied().unwrap_or(Vec3::ZERO),
                        w.get(at).copied().unwrap_or(f64::NAN),
                    )
                })
                .collect()
        })
        .collect()
}

/// How far the net `Build` stores is from the Hermite rows it was made of, span by span
/// (`spans`, each [`NetForm::span`]'s): its control points are the rows divided by their
/// weights and its weights scaled, each rounded, so the net is off the rows by those
/// roundings, which the definition's bound takes in.
///
/// # Errors
///
/// `SM_ERR` where the net is not the spans' (its size); those of [`net_distance`].
pub(crate) fn stored_net_distance(
    net: &BSplineSurface,
    spans: &[(Vec<Vec3>, Vec<f64>)],
    ld: usize,
) -> Result<f64> {
    let (nu, nv, stored) = net.control_points()?;
    if nv != ld || nu != 3 * spans.len() + 1 {
        return Err(SolidError::sm_err(OP));
    }
    let mut d = 0.0_f64;
    for (k, (h, w)) in spans.iter().enumerate() {
        let b: Vec<Vec<Hpt>> = (0..4)
            .map(|i| {
                (0..ld)
                    .map(|j| {
                        stored
                            .get((3 * k + i) * ld + j)
                            .map_or((Vec3::ZERO, f64::NAN), |c| {
                                (Vec3::new(c.x, c.y, c.z), if c.w == NOW { 1.0 } else { c.w })
                            })
                    })
                    .collect()
            })
            .collect();
        d = d.max(net_distance(&b, &span_grid(h, w, ld))?);
    }
    Ok(d)
}

/// The net's Bezier pieces over the definition's (s, t), each net ordered by increasing s
/// and t.
fn pieces(pb: &ProceduralBlend, ld: usize) -> Result<Vec<Piece>> {
    let not_its = || SolidError::invalid_input(OP, "the net is not of its definition's degrees");
    let d = pb.definition();
    let nurb = pb.net().gw_nurb().ok_or_else(|| SolidError::sm_err(OP))?;
    let patches = decompose_to_bez(nurb)?;
    let mut out = Vec::new();
    for row in &patches {
        for patch in row {
            let (p, q) = patch.degrees();
            // (p, q) along (u, v); s is u unless the map is swapped
            let (ds, dt) = if d.swapped { (q, p) } else { (p, q) };
            if ds != 3 || dt + 1 != ld {
                return Err(not_its());
            }
            let (Some(&u0), Some(&u1), Some(&v0), Some(&v1)) = (
                patch.u_knots().first(),
                patch.u_knots().last(),
                patch.v_knots().first(),
                patch.v_knots().last(),
            ) else {
                return Err(SolidError::sm_err(OP));
            };
            let (ps, pt) = if d.swapped {
                ([v0, v1], [u0, u1])
            } else {
                ([u0, u1], [v0, v1])
            };
            let ([a0, a1], [b0, b1]) = (d.map_a, d.map_b);
            let mut s = [a0 + b0 * ps[0], a0 + b0 * ps[1]];
            let mut t = [a1 + b1 * pt[0], a1 + b1 * pt[1]];
            let flip_s = d.map_b[0] < 0.0;
            let flip_t = d.map_b[1] < 0.0;
            if flip_s {
                s.swap(0, 1);
            }
            if flip_t {
                t.swap(0, 1);
            }
            // [i along s][j along t]
            let mut net = vec![vec![(Vec3::ZERO, 0.0); dt + 1]; ds + 1];
            for (iu, prow) in patch.rows().enumerate() {
                for (iv, c) in prow.iter().enumerate() {
                    let (is, it) = if d.swapped { (iv, iu) } else { (iu, iv) };
                    let is = if flip_s { ds - is } else { is };
                    let it = if flip_t { dt - it } else { it };
                    // (stored weighted when rational: homogeneous either way)
                    let w = if c.w == NOW { 1.0 } else { c.w };
                    if let Some(x) = net.get_mut(is).and_then(|r| r.get_mut(it)) {
                        *x = (Vec3::new(c.x, c.y, c.z), w);
                    }
                }
            }
            out.push(Piece { net, s, t });
        }
    }
    Ok(out)
}

/// [`ProceduralBlend::recertify`]: the net's certified distance from the definition, at
/// most the definition's bound.
///
/// # Errors
///
/// `InvalidInput` where the net is not proved within the bound, or is not a net of its
/// definition's stations (its degrees, a piece across a station or past the seeds, the
/// seeds' stations differing, a foot crossing a face's knot line inside a span);
/// `SM_ERR` for a support that is not a NURBS surface or is itself a blend by its
/// definition (its definition is not enclosed yet); those of the station solves and the
/// certificate.
#[allow(clippy::too_many_lines)]
pub fn recertify(k: &Kernel, pb: &ProceduralBlend) -> Result<f64> {
    let d = pb.definition();
    let bound = d.bound;
    let n = d.order as usize;
    let equal_area = n >= 2 && d.thumbweight == GN_THUMBWEIGHT_EQUAL_AREA;
    let net_form = NetForm {
        n,
        thumbweight: d.thumbweight,
        equal_area,
        law: d.law.is_some(),
    };
    let ld = net_form.width();
    let [s0, s1] = pb.supports();
    let faces = [s0, s1];
    let ([n0, n1], mut defs) = certificate_faces(k, faces, n, d.law.is_some())?;
    let any_defs = defs.faces.iter().any(Option::is_some);
    let face_n = [&n0, &n1];
    let leader_n = NurbsCurve::of(&d.leader)?;
    let law = match &d.law {
        Some(c) => Some(GnLaw::of(c)?),
        None => None,
    };

    // the stations: the seeds' knots
    let [seed0, seed1] = &d.seeds;
    let sigma = seed0.knots(None)?.0;
    if seed1.knots(None)?.0 != sigma || sigma.len() < 2 {
        return Err(SolidError::invalid_input(
            OP,
            "the seeds are not on one set of stations",
        ));
    }
    let r = match &law {
        Some(l) => {
            let first = sigma.first().copied().unwrap_or(f64::NAN);
            let r0 = l.curve.evaluate_point(first)?.x;
            [r0, r0]
        }
        None => d.offsets,
    };
    let tr = Trunc::new(
        if law.is_some() { 2 } else { 1 },
        if n >= 2 { n } else { 0 },
    )?;
    let system = StationSystem {
        faces,
        leader: &d.leader,
        r,
        side: d.sides,
        tr,
        n,
        equal_area,
        law: law.as_ref(),
    };

    // each piece's span
    let pieces = pieces(pb, ld)?;
    let precision = net_precision(pb.net())?;
    let mut span_of = Vec::new();
    let mut used = vec![false; sigma.len()];
    for piece in &pieces {
        let mid = 0.5 * (piece.s[0] + piece.s[1]);
        let span = sigma
            .windows(2)
            .position(|w| matches!(w, [a, b] if *a <= mid && mid <= *b))
            .ok_or_else(|| SolidError::invalid_input(OP, "a piece of the net is past its seeds"))?;
        let (Some(&a), Some(&b)) = (sigma.get(span), sigma.get(span + 1)) else {
            return Err(SolidError::sm_err(OP));
        };
        let room_s = |x: f64| map_rounding(d.map_a[0], d.map_b[0], (x - d.map_a[0]) / d.map_b[0]);
        let room_t = |x: f64| map_rounding(d.map_a[1], d.map_b[1], (x - d.map_a[1]) / d.map_b[1]);
        let inside = piece.s[0] >= a - room_s(a)
            && piece.s[1] <= b + room_s(b)
            && piece.t[0] >= -room_t(0.0)
            && piece.t[1] <= 1.0 + room_t(1.0);
        if !inside {
            return Err(SolidError::invalid_input(
                OP,
                "a piece of the net runs across a station of its definition",
            ));
        }
        span_of.push(span);
        for i in [span, span + 1] {
            if let Some(u) = used.get_mut(i) {
                *u = true;
            }
        }
    }

    // the stations the pieces need, on the pieces after them and (but the first) before
    let mut stations: Vec<Option<GnStation>> = vec![None; sigma.len()];
    for (j, (&s, &needed)) in sigma.iter().zip(&used).enumerate() {
        if !needed {
            continue;
        }
        let mut seed = [0.0_f64; 4];
        let mut after = [0_i32; 4];
        let mut before = [0_i32; 4];
        for (f, sd) in [seed0, seed1].into_iter().enumerate() {
            let p = sd.evaluate_point(s)?;
            let right = sd.evaluate(s, 1, false, false)?;
            let left = sd.evaluate(s, 1, true, false)?;
            let (Some(dr), Some(dl)) = (right.get(1), left.get(1)) else {
                return Err(SolidError::sm_err(OP));
            };
            for (c, (x, (r1, l1))) in [(p.x, (dr.x, dl.x)), (p.y, (dr.y, dl.y))]
                .into_iter()
                .enumerate()
            {
                let q = 2 * f + c;
                if let (Some(sq), Some(aq), Some(bq)) =
                    (seed.get_mut(q), after.get_mut(q), before.get_mut(q))
                {
                    *sq = x;
                    *aq = sign(r1);
                    *bq = -sign(l1);
                }
            }
        }
        let mut st = system.solve_after(k, s, seed, Some(after))?;
        st.before = before;
        if j > 0 {
            system.solve_before(k, &mut st)?;
        }
        if let Some(slot) = stations.get_mut(j) {
            *slot = Some(st);
        }
    }
    let solved: Vec<GnStation> = stations.iter().flatten().copied().collect();
    let mut vary = varying(n, r[0] * r[1], &solved, equal_area);

    // span by span: the nets compared, then the certificate for what the bound leaves
    let mut achieved = 0.0_f64;
    let spans: Vec<usize> = {
        let mut s = span_of.clone();
        s.sort_unstable();
        s.dedup();
        s
    };
    for span in spans {
        let (Some(Some(a)), Some(Some(b))) = (stations.get(span), stations.get(span + 1)) else {
            return Err(SolidError::sm_err(OP));
        };
        // a foot crossing a face's knot line inside the span: the construction puts a
        // station there, so this net is not its
        for (f, nf) in face_n.iter().enumerate() {
            for c in 0..2 {
                let q = 2 * f + c;
                let (xa, xb) = (param(a, q, 0, 0), param(b, q, 0, 0));
                if inner_knots(nf, c)
                    .iter()
                    .any(|&dk| crosses_knot(xa, xb, dk))
                {
                    return Err(SolidError::invalid_input(
                        OP,
                        "a foot crosses a face's knot line between two stations",
                    ));
                }
            }
        }
        let (ra, rb) = (net_form.station_rows(a)?, net_form.station_rows(b)?);
        let (homog, w) = net_form.span(a, &ra, b, &rb);
        let h = b.s - a.s;
        let mut delta = 0.0_f64;
        for (piece, _) in pieces.iter().zip(&span_of).filter(|(_, &sp)| sp == span) {
            // H on the piece: cut along s, then across
            let mut cut: Vec<Vec<Hpt>> = Vec::new();
            let lo = ((piece.s[0] - a.s) / h).max(0.0);
            let hi = ((piece.s[1] - a.s) / h).min(1.0);
            let grid = span_grid(&homog, &w, ld);
            let mut columns: Vec<Vec<Hpt>> = Vec::new();
            for j in 0..ld {
                let col: Vec<Hpt> = grid
                    .iter()
                    .map(|row| row.get(j).copied().unwrap_or((Vec3::ZERO, f64::NAN)))
                    .collect();
                columns.push(sub_curve(&col, lo, hi));
            }
            for i in 0..4 {
                let row: Vec<Hpt> = columns
                    .iter()
                    .map(|c| c.get(i).copied().unwrap_or((Vec3::ZERO, f64::NAN)))
                    .collect();
                cut.push(sub_curve(&row, piece.t[0].max(0.0), piece.t[1].min(1.0)));
            }
            delta = delta.max(net_distance(&piece.net, &cut)?);
        }
        // (differences below the net's precision are not the net's)
        let delta = delta.max(precision);
        let target = bound - delta;
        if target.partial_cmp(&0.0) != Some(Ordering::Greater) {
            return Err(SolidError::invalid_input(
                OP,
                "the net is not within its bound of its definition",
            ));
        }
        let certified = gn_certify_span(
            k,
            face_n,
            &leader_n,
            faces,
            &d.leader,
            a,
            b,
            r,
            d.sides,
            n,
            target,
            if equal_area { Some(&mut vary) } else { None },
            law.as_ref(),
            if any_defs { Some(&mut defs) } else { None },
        )?;
        let total = certified + delta;
        if !matches!(
            total.partial_cmp(&bound),
            Some(Ordering::Less | Ordering::Equal)
        ) {
            return Err(SolidError::invalid_input(
                OP,
                "the net is not within its bound of its definition",
            ));
        }
        achieved = achieved.max(total);
    }
    Ok(achieved)
}

/// -1, 0 or +1.
fn sign(x: f64) -> i32 {
    if x > 0.0 {
        1
    } else if x < 0.0 {
        -1
    } else {
        0
    }
}
