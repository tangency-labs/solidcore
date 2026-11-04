//! The general-path motion of a [`ProceduralSweep`](super::ProceduralSweep) (#2572, slice 2,
//! not upstream): the profile carried along a path with a frame, and the certified distance
//! of the net `sweep_scale` fits to it.
//!
//! **The definition.** `S(w, t) = T(t) + x(t) C_x(w) + y(t) C_y(w) + z(t) C_z(w)`, `T` the
//! path, `C` the profile in the frame's own coordinates and `t` the path's parameter. The
//! frame is the projection normal method of `N_CreateBVectors` in closed form: `x = T' /
//! |T'|`, and from the last station `j` of the frame chain (a parameter `t_j <= t`, with the
//! unit normal `z_j` the chain holds there)
//!
//! `z(t) = (z_j - (z_j . x) x) / |z_j - (z_j . x) x|`, `y = z x x`.
//!
//! At each station this is the net's own frame (the chain is that rule applied from one
//! station to the next); between them it is the frame the rule gives from the last, so `S`
//! is continuous and its sections are everywhere perpendicular to the path. On a straight
//! `x` is constant and so is the frame. The chain is the net's, not a curve's: a frame taken
//! from a curve or a closed path has no definition.
//!
//! **The map.** The net's parameter across the sweep is the skin's chord length, not the
//! path's `t`. The definition carries nodes `(v_k, t_k)`, the net's parameter at each
//! placed section (found by a foot point search, whose result is data: the bound below holds
//! for whatever map is stored), and reads `t` by interpolating them linearly.
//!
//! **The bound.** Every section of the net is the profile's control points moved by one
//! motion, in the same rational basis `R_i(u)` (non-negative, summing to 1) for every `v`,
//! so with `P_i(v)` the net's `i`-th control point as a curve of `v` and `S_i(v) = T + R C_i`
//! the same point of `S`,
//! `B(u, v) - S(u, v) = sum_i R_i(u) (P_i(v) - S_i(v))`: a convex combination, and
//! `|B - S| <= max_i max_v |P_i(v) - S_i(v)|`. Each `f_i = P_i - S_i` is bounded over a cell
//! of `v` about its middle `c`, half width `h`, by its Taylor expansion,
//! `|f_i| <= |f_i(c)| + |f_i'(c)| h + (M_P + M_S) h^2 / 2`, `M_P` the largest `|P_i''|` in the
//! cell (its second derivative's control points, a B-spline's convex hull) and `M_S` that of
//! `S_i''`, from bounds of `|T''|`, `|T'''|`, `|T'|` over the path's span (control polygons of
//! its derivatives) and of the frame's derivatives (below). Cells never cross a knot of the
//! net, of the path, a node or a station, where `P_i''` or `S_i''` jump, and are split until
//! the second order term is small or the budget [`MAX_CELLS`] is spent (the bound is then
//! looser, never wrong).
//!
//! With `a = |T''| / |T'|_min` and `b = |T'''| / |T'|_min + 3 a^2`, `|x'| <= a`, `|x''| <= b`;
//! for `z = u / |u|`, `u = z_j - (z_j . x) x`, `|u| >= n = cos(alpha)` for the angle
//! `alpha` the path's tangent has turned from station `j` (at most the station's span times
//! the largest `a`), `|u'| <= a (1 + sin alpha)`, `|u''| <= b (1 + sin alpha) + 2 a^2`,
//! `|z'| <= 2 |u'| / n`, `|z''| <= 2 |u''| / n + 6 |u'|^2 / n^2`, and for `y = z x x`
//! `|y''| <= |z''| + 2 |z'| a + b`.
//!
//! The floor adds a forward error bound of evaluating both in floating point, `gamma_K`
//! times the sizes of the numbers, and a slack of `2^-30` relative on the cell bounds (they
//! are computed in doubles).
//!
//! Deviations: none (no upstream counterpart). No recursion (an explicit stack of cells).

use crate::error::{Result, SolidError};
use crate::geom::axis2_placement::Axis2Placement;
use crate::geom::vector3d::Vec3;
use crate::nurbs::b_spline_curve::BSplineCurve;
use crate::nurbs::b_spline_surface::BSplineSurface;
use crate::nurbs::math::{gamma, MODEL_TOL};
use crate::nurbs::srf_adv::p6::SweepTrace;

use crate::blend::fillet_section::interval::round_up;

const OP_CREATE: &str = "surface::procedural_sweep::PathMotion::create";
const OP_EVAL: &str = "surface::procedural_sweep::PathMotion::frame";

/// The most cells the bound is split into.
pub const MAX_CELLS: usize = 400_000;

/// The deepest a cell is halved.
const MAX_DEPTH: u32 = 40;

/// The most the tangent may turn from a station (radians): `cos(0.6) > 0.82`.
const MAX_TURN: f64 = 0.6;

/// Slack on the doubles' cell bounds, relative.
const SLACK: f64 = 1.0 / 1_073_741_824.0;

fn bad(reason: &'static str) -> SolidError {
    SolidError::invalid_input(OP_CREATE, reason)
}

/// A B-spline curve of 3D points, unrounded evaluation and its derivatives: the net's
/// control columns and the path's derivatives.
#[derive(Debug, Clone)]
struct Spline {
    degree: usize,
    knots: Vec<f64>,
    ctrl: Vec<Vec3>,
}

impl Spline {
    /// The spline through `ctrl` over the full `knots` (`ctrl.len() + degree + 1` of them,
    /// ascending).
    fn new(degree: usize, knots: Vec<f64>, ctrl: Vec<Vec3>) -> Result<Self> {
        let ok = !ctrl.is_empty()
            && knots.len() == ctrl.len() + degree + 1
            && knots.windows(2).all(|w| matches!(w, [a, b] if a <= b))
            && knots.iter().all(|k| k.is_finite())
            && ctrl
                .iter()
                .all(|p| p.x.is_finite() && p.y.is_finite() && p.z.is_finite());
        if ok {
            Ok(Self {
                degree,
                knots,
                ctrl,
            })
        } else {
            Err(bad(
                "a spline's knots and control points are finite and fit",
            ))
        }
    }

    /// The derivative, a spline of one degree less (of degree 0 and zero for a constant).
    fn derivative(&self) -> Self {
        let p = self.degree;
        if p == 0 {
            return Self {
                degree: 0,
                knots: self.knots.clone(),
                ctrl: vec![Vec3::ZERO; self.ctrl.len()],
            };
        }
        let factor = f64::from(u32::try_from(p).unwrap_or(u32::MAX));
        let ctrl = self
            .ctrl
            .windows(2)
            .enumerate()
            .map(|(i, w)| match w {
                [a, b] => {
                    let lo = self.knots.get(i + 1).copied().unwrap_or(0.0);
                    let hi = self.knots.get(i + p + 1).copied().unwrap_or(0.0);
                    if hi > lo {
                        (*b - *a) * (factor / (hi - lo))
                    } else {
                        Vec3::ZERO
                    }
                }
                _ => Vec3::ZERO,
            })
            .collect();
        let knots = self
            .knots
            .get(1..self.knots.len().saturating_sub(1))
            .map(<[f64]>::to_vec)
            .unwrap_or_default();
        Self {
            degree: p - 1,
            knots,
            ctrl,
        }
    }

    /// The knot span `m`, `degree <= m < ctrl.len()`, holding `t` (clamped to the domain).
    fn span(&self, t: f64) -> usize {
        let p = self.degree;
        let n = self.ctrl.len();
        let m = self
            .knots
            .partition_point(|&k| k <= t)
            .saturating_sub(1)
            .max(p);
        m.min(n.saturating_sub(1))
    }

    /// The point at `t` (de Boor).
    fn eval(&self, t: f64) -> Vec3 {
        let p = self.degree;
        let m = self.span(t);
        let mut d: Vec<Vec3> = (0..=p)
            .map(|j| {
                (m + j)
                    .checked_sub(p)
                    .and_then(|i| self.ctrl.get(i))
                    .copied()
                    .unwrap_or(Vec3::ZERO)
            })
            .collect();
        for r in 1..=p {
            for j in (r..=p).rev() {
                let i = j + m - p;
                let lo = self.knots.get(i).copied().unwrap_or(0.0);
                let hi = self.knots.get(i + p + 1 - r).copied().unwrap_or(lo);
                let alpha = if hi > lo { (t - lo) / (hi - lo) } else { 0.0 };
                let prev = d.get(j - 1).copied().unwrap_or(Vec3::ZERO);
                let cur = d.get(j).copied().unwrap_or(Vec3::ZERO);
                if let Some(slot) = d.get_mut(j) {
                    *slot = prev * (1.0 - alpha) + cur * alpha;
                }
            }
        }
        d.get(p).copied().unwrap_or(Vec3::ZERO)
    }

    /// The control points that carry the span holding `t`.
    fn local(&self, t: f64) -> &[Vec3] {
        let m = self.span(t);
        self.ctrl
            .get(m.saturating_sub(self.degree)..=m)
            .unwrap_or(&[])
    }

    /// An upper bound of `|S(x)|` over the span holding `t` (convex hull).
    fn largest(&self, t: f64) -> f64 {
        self.local(t).iter().fold(0.0_f64, |m, p| m.max(p.length()))
    }

    /// The domain.
    fn domain(&self) -> (f64, f64) {
        let p = self.degree;
        let lo = self.knots.get(p).copied().unwrap_or(0.0);
        let hi = self.knots.get(self.ctrl.len()).copied().unwrap_or(lo);
        (lo, hi)
    }
}

/// The full knot vector of unique knots and multiplicities.
fn full_knots(unique: &[f64], mult: &[usize]) -> Vec<f64> {
    unique
        .iter()
        .zip(mult)
        .flat_map(|(&k, &m)| std::iter::repeat_n(k, m))
        .collect()
}

/// A frame of the definition at a parameter, with its first derivatives in `t`.
#[derive(Debug, Clone, Copy)]
pub struct Frame {
    /// The path's point.
    pub point: Vec3,
    /// The path's first derivative.
    pub tangent: Vec3,
    /// The frame's axes `x`, `y`, `z`.
    pub axes: [Vec3; 3],
    /// Their derivatives in `t`.
    pub rates: [Vec3; 3],
}

impl Frame {
    /// The frame's image of the profile's `c`: `T + x c_x + y c_y + z c_z`.
    pub fn place(&self, c: Vec3) -> Vec3 {
        let [x, y, z] = self.axes;
        self.point + x * c.x + y * c.y + z * c.z
    }

    /// The rate of that image in `t`.
    pub fn rate(&self, c: Vec3) -> Vec3 {
        let [x, y, z] = self.rates;
        self.tangent + x * c.x + y * c.y + z * c.z
    }

    /// The image's rotation of a profile vector `c`.
    pub fn turn(&self, c: Vec3) -> Vec3 {
        let [x, y, z] = self.axes;
        x * c.x + y * c.y + z * c.z
    }

    /// The rate of the rotation of `c`.
    pub fn turn_rate(&self, c: Vec3) -> Vec3 {
        let [x, y, z] = self.rates;
        x * c.x + y * c.y + z * c.z
    }
}

/// A general path and the frame it carries the profile in (module docs).
#[derive(Debug, Clone)]
pub struct PathMotion {
    path: BSplineCurve,
    /// The frame chain: parameter `t_j` and the unit normal `z_j`, ascending.
    chain: Vec<(f64, Vec3)>,
    /// The net's parameter across the sweep and the path's, at each section, both strictly
    /// ascending.
    nodes: Vec<(f64, f64)>,
}

impl PathMotion {
    /// The path.
    #[must_use]
    pub const fn path(&self) -> &BSplineCurve {
        &self.path
    }

    /// The frame chain: the parameter of each station and the unit normal there.
    #[must_use]
    pub fn chain(&self) -> &[(f64, Vec3)] {
        &self.chain
    }

    /// The net's parameter and the path's at each placed section, both ascending: the path's
    /// `t` at a net parameter in between is interpolated linearly.
    #[must_use]
    pub fn nodes(&self) -> &[(f64, f64)] {
        &self.nodes
    }

    /// The path's `t` at the net's parameter `v` (clamped to the nodes) and `dt / dv`.
    #[must_use]
    pub fn path_parameter(&self, v: f64) -> (f64, f64) {
        let (Some(&(v0, t0)), Some(&(v1, t1))) = (self.nodes.first(), self.nodes.last()) else {
            return (v, 1.0);
        };
        if !(v > v0) {
            return (t0, 1.0);
        }
        if !(v < v1) {
            return (t1, 1.0);
        }
        let i = self
            .nodes
            .partition_point(|&(x, _)| x <= v)
            .saturating_sub(1);
        match (self.nodes.get(i), self.nodes.get(i + 1)) {
            (Some(&(va, ta)), Some(&(vb, tb))) if vb > va => {
                let m = (tb - ta) / (vb - va);
                (ta + m * (v - va), m)
            }
            _ => (v, 1.0),
        }
    }

    /// The frame at the path's parameter `t` (clamped to the path's domain).
    ///
    /// # Errors
    ///
    /// `InvalidInput` where the path's tangent vanishes or turned too far from its station;
    /// those of evaluating the path.
    pub fn frame(&self, t: f64) -> Result<Frame> {
        let d = self.path.evaluate(t, 2, true, false)?;
        let (Some(&point), Some(&t1), Some(&t2)) = (d.first(), d.get(1), d.get(2)) else {
            return Err(SolidError::invalid_input(
                OP_EVAL,
                "the path has no derivatives",
            ));
        };
        let speed = t1.length();
        if !(speed.is_finite() && speed > 0.0) {
            return Err(SolidError::invalid_input(
                OP_EVAL,
                "the path's tangent vanishes",
            ));
        }
        let x = t1 * (1.0 / speed);
        let dx = (t2 - x * t2.dot(x)) * (1.0 / speed);
        let j = self
            .chain
            .partition_point(|&(tj, _)| tj <= t)
            .saturating_sub(1);
        let zj = self.chain.get(j).map_or(Vec3::Z, |c| c.1);
        let s = zj.dot(x);
        let u = zj - x * s;
        let du = -(x * zj.dot(dx)) - dx * s;
        let n = u.length();
        if !(n.is_finite() && n > 1e-6) {
            return Err(SolidError::invalid_input(
                OP_EVAL,
                "the frame's normal vanishes: the path turned too far from its station",
            ));
        }
        let z = u * (1.0 / n);
        let dz = (du - z * z.dot(du)) * (1.0 / n);
        let y = z.cross(x);
        let dy = dz.cross(x) + z.cross(dx);
        Ok(Frame {
            point,
            tangent: t1,
            axes: [x, y, z],
            rates: [dx, dy, dz],
        })
    }

    /// Moves the path with its net: `rotate_n_move` after `scale` (a uniform one); the
    /// profile, in the frame's own coordinates, takes the scale only (the caller's).
    ///
    /// # Errors
    ///
    /// Those of the path's transform.
    pub(crate) fn transform(
        &mut self,
        rotate_n_move: &Axis2Placement,
        scale: Option<Vec3>,
    ) -> Result<()> {
        self.path.transform(rotate_n_move, scale)?;
        for (_, z) in &mut self.chain {
            *z = rotate_n_move.transform_vector(*z);
        }
        Ok(())
    }

    /// How far the path reaches from the origin.
    pub(crate) fn reach(&self) -> Result<f64> {
        let c = self.path.canonical()?;
        Ok(c.control_points
            .iter()
            .fold(0.0_f64, |m, p| m.max(p.length())))
    }

    /// The motion of a net `net` of the profile `profile` along `path`, as the sweep placed
    /// the sections in `trace`; with its certified bound `|B - S| <= bound` and whether the
    /// net's `(u, v)` are the sweep's `(v, u)`.
    ///
    /// # Errors
    ///
    /// `InvalidInput` for a net that is not the profile's sections moved along the path (its
    /// control points per section, its weights or its knots differ), a path that is rational
    /// or has a vanishing tangent, a frame chain that turns too far between stations, a
    /// section whose net parameter is not found or a bound that is not finite.
    pub(crate) fn create(
        profile: &BSplineCurve,
        path: &BSplineCurve,
        trace: &SweepTrace,
        net: &BSplineSurface,
    ) -> Result<(Self, f64, bool)> {
        let pc = path.canonical()?;
        if pc.dim != 3 {
            return Err(bad("a path is a 3D curve"));
        }
        let closed = match (pc.control_points.first(), pc.control_points.last()) {
            (Some(a), Some(b)) => (*a - *b).length() <= 1e-9 * (1.0 + a.length()),
            _ => true,
        };
        if closed {
            return Err(bad(
                "a closed path's frame chain is corrected, not the rule's",
            ));
        }
        let mut path = path.copy();
        path.set_out_of_bounds_enabled(false);
        // the frame chain, from the rule
        let mut chain: Vec<(f64, Vec3)> = Vec::with_capacity(trace.chain.len());
        for &t in &trace.chain {
            let tangent = path
                .evaluate(t, 1, true, false)?
                .get(1)
                .copied()
                .unwrap_or(Vec3::ZERO);
            let l = tangent.length();
            if !(l.is_finite() && l > 0.0) {
                return Err(bad("the path's tangent vanishes at a station"));
            }
            let x = tangent * (1.0 / l);
            let z = match chain.last() {
                None => trace.start_normal,
                Some(&(_, z)) => {
                    let u = z - x * z.dot(x);
                    let n = u.length();
                    if !(n.is_finite() && n > 0.0) {
                        return Err(bad("the frame chain's normal vanishes"));
                    }
                    u * (1.0 / n)
                }
            };
            chain.push((t, z));
        }
        if chain.len() < 2
            || chain
                .windows(2)
                .any(|w| matches!(w, [a, b] if !(a.0 < b.0)))
        {
            return Err(bad("a frame chain has at least two ascending stations"));
        }
        let mut motion = Self {
            path,
            chain,
            nodes: Vec::new(),
        };
        let columns = Columns::read(profile, net)?;
        motion.nodes = motion.find_nodes(&columns, &trace.pieces)?;
        let bound = motion.certify(&columns)?;
        Ok((motion, bound, columns.swapped))
    }

    /// The net's parameter at each section: the skin's own. Each piece of the path is
    /// skinned with parameters in `[0, 1]` that are the averaged normalised chord lengths of
    /// the sections' control points (`interp_params`, in homogeneous coordinates), and the
    /// pieces are joined with their parameters shifted to follow one another: the `k`-th
    /// piece's are `k + tl`. The map is data (the bound holds for any), so where the net's
    /// own is not these (a wall whose domain does not end at the piece count) the sweep has
    /// no definition rather than a loose one.
    fn find_nodes(&self, cols: &Columns, pieces: &[Vec<f64>]) -> Result<Vec<(f64, f64)>> {
        let (v0, v1) = cols.domain;
        let count = f64::from(u32::try_from(pieces.len()).unwrap_or(u32::MAX));
        if pieces.is_empty() || !(v0 == 0.0 && (v1 - count).abs() <= 1e-9 * count) {
            return Err(bad(
                "a net's parameter across the sweep is not its pieces' skins",
            ));
        }
        let mut nodes: Vec<(f64, f64)> = Vec::new();
        for (k, piece) in pieces.iter().enumerate() {
            let offset = f64::from(u32::try_from(k).unwrap_or(u32::MAX));
            let tl = self.chord_parameters(cols, piece)?;
            for (&t, &p) in piece.iter().zip(&tl) {
                let v = offset + p;
                match nodes.last() {
                    Some(&(pv, pt)) if v == pv && t == pt => {}
                    Some(&(pv, pt)) if !(v > pv && t > pt) => {
                        return Err(bad("the net's sections are not in the path's order"));
                    }
                    _ => nodes.push((v, t)),
                }
            }
        }
        Ok(nodes)
    }

    /// `interp_params` of the sections of one piece: the sections' control points in
    /// homogeneous coordinates, the cumulative chord lengths normalised per control point
    /// and averaged over those with a length above `MODEL_TOL`.
    fn chord_parameters(&self, cols: &Columns, piece: &[f64]) -> Result<Vec<f64>> {
        let n = piece.len();
        if n < 2 {
            return Err(bad("a piece has at least two sections"));
        }
        let mut sections: Vec<Vec<[f64; 4]>> = Vec::with_capacity(n);
        for &t in piece {
            let frame = self.frame(t)?;
            sections.push(
                cols.profile
                    .iter()
                    .zip(&cols.weights)
                    .map(|(c, &w)| {
                        let p = frame.place(*c);
                        if cols.rational {
                            [w * p.x, w * p.y, w * p.z, w]
                        } else {
                            [p.x, p.y, p.z, 0.0]
                        }
                    })
                    .collect(),
            );
        }
        let mut t = vec![0.0_f64; n];
        let mut used = 0.0_f64;
        for i in 0..cols.profile.len() {
            let mut cum = vec![0.0_f64; n];
            let mut sum = 0.0;
            for k in 1..n {
                let (Some(a), Some(b)) = (
                    sections.get(k - 1).and_then(|s| s.get(i)),
                    sections.get(k).and_then(|s| s.get(i)),
                ) else {
                    return Err(bad("a section lacks a control point"));
                };
                let d = a
                    .iter()
                    .zip(b)
                    .map(|(x, y)| (x - y) * (x - y))
                    .sum::<f64>()
                    .sqrt();
                sum += d;
                if let Some(slot) = cum.get_mut(k) {
                    *slot = sum;
                }
            }
            if sum > MODEL_TOL {
                used += 1.0;
                for (slot, c) in t.iter_mut().zip(&cum) {
                    *slot += c / sum;
                }
            }
        }
        if !(used > 0.0) {
            return Err(bad("a piece's sections are apart"));
        }
        for slot in &mut t {
            *slot /= used;
        }
        if let Some(last) = t.last_mut() {
            *last = 1.0;
        }
        Ok(t)
    }

    /// The certified bound (module docs).
    fn certify(&self, cols: &Columns) -> Result<f64> {
        let bounds = PathBounds::new(&self.path)?;
        let (v0, v1) = cols.domain;
        let to_v = |t: f64| -> Option<f64> {
            let i = self
                .nodes
                .partition_point(|&(_, x)| x <= t)
                .saturating_sub(1)
                .min(self.nodes.len().saturating_sub(2));
            match (self.nodes.get(i), self.nodes.get(i + 1)) {
                (Some(&(va, ta)), Some(&(vb, tb))) if tb > ta => {
                    Some(va + (vb - va) * (t - ta) / (tb - ta))
                }
                _ => None,
            }
        };
        // cells never cross a knot of the net, a node, a station or a knot of the path
        let mut cuts: Vec<f64> = cols.unique_knots.clone();
        cuts.extend(self.nodes.iter().map(|n| n.0));
        cuts.extend(self.chain.iter().filter_map(|c| to_v(c.0)));
        cuts.extend(bounds.unique_knots.iter().filter_map(|&k| to_v(k)));
        cuts.retain(|&v| v >= v0 && v <= v1);
        cuts.push(v0);
        cuts.push(v1);
        cuts.sort_by(f64::total_cmp);
        cuts.dedup_by(|a, b| (*a - *b).abs() <= 1e-14 * (1.0 + b.abs()));
        // the sizes the floor and the goal are relative to
        let reach = cols.profile.iter().fold(0.0_f64, |m, c| m.max(c.length()));
        let big = cols
            .splines
            .iter()
            .flat_map(|s| s.ctrl.iter())
            .fold(0.0_f64, |m, p| m.max(p.length()))
            .max(self.reach()? + reach);
        let goal = 1e-7 * (1.0 + big);
        let firsts: Vec<Spline> = cols.splines.iter().map(Spline::derivative).collect();
        let seconds: Vec<Spline> = firsts.iter().map(Spline::derivative).collect();
        let mut worst = 0.0_f64;
        let mut cells = 0_usize;
        for w in cuts.windows(2) {
            let &[a, b] = w else { continue };
            if !(b > a) {
                continue;
            }
            let mut stack: Vec<(f64, f64, u32)> = vec![(a, b, 0)];
            while let Some((lo, hi, depth)) = stack.pop() {
                cells += 1;
                let c = 0.5 * (lo + hi);
                let h = 0.5 * (hi - lo);
                let (t, m) = self.path_parameter(c);
                let dt = m.abs() * h;
                let frame = self.frame(t)?;
                // the path's derivative bounds over the cell's span, and the tangent's lower
                // bound there
                let (_, d2, d3) = bounds.derivative_bounds(t);
                let smin = frame.tangent.length() - d2 * dt;
                if !(smin > 0.0) {
                    return Err(bad("a path's tangent is not bounded away from zero"));
                }
                // the angle the tangent has turned from the station's, at most
                let j = self
                    .chain
                    .partition_point(|&(tj, _)| tj <= t)
                    .saturating_sub(1);
                let tj = self.chain.get(j).map_or(t, |c| c.0);
                let xj = self.path.evaluate(tj, 1, true, false)?;
                let xj = xj.get(1).copied().unwrap_or(Vec3::ZERO);
                let xj = xj * (1.0 / xj.length());
                let a1 = d2 / smin;
                let turned = frame.axes[0].dot(xj).clamp(-1.0, 1.0).acos() + a1 * dt;
                let can_split = depth < MAX_DEPTH && cells + stack.len() < MAX_CELLS;
                if !(turned <= MAX_TURN) {
                    if can_split {
                        stack.push((c, hi, depth + 1));
                        stack.push((lo, c, depth + 1));
                        continue;
                    }
                    return Err(bad("the path turns too far between frame stations"));
                }
                let sin_a = turned.sin();
                let n = turned.cos();
                let a2 = d3 / smin + 3.0 * a1 * a1;
                let u1 = a1 * (1.0 + sin_a);
                let u2 = a2 * (1.0 + sin_a) + 2.0 * a1 * a1;
                let z1 = 2.0 * u1 / n;
                let z2 = 2.0 * u2 / n + 6.0 * u1 * u1 / (n * n);
                let y2 = z2 + 2.0 * z1 * a1 + a2;
                let mut cell = 0.0_f64;
                let mut split = false;
                for (((p, p1), p2), cc) in cols
                    .splines
                    .iter()
                    .zip(&firsts)
                    .zip(&seconds)
                    .zip(&cols.profile)
                {
                    let f = p.eval(c) - frame.place(*cc);
                    let df = p1.eval(c) - frame.rate(*cc) * m;
                    let mp = p2.largest(c);
                    let ms = m * m * (d2 + cc.x.abs() * a2 + cc.y.abs() * y2 + cc.z.abs() * z2);
                    let rem = 0.5 * (mp + ms) * h * h;
                    if rem > goal && can_split {
                        split = true;
                        break;
                    }
                    cell = cell.max(f.length() + df.length() * h + rem);
                }
                if split {
                    stack.push((c, hi, depth + 1));
                    stack.push((lo, c, depth + 1));
                } else {
                    worst = worst.max(cell);
                }
            }
        }
        let floor = gamma(96) * (big + reach) * 4.0;
        let bound = round_up(worst * (1.0 + SLACK) + floor);
        if bound.is_finite() && bound > 0.0 {
            Ok(bound)
        } else {
            Err(bad(
                "the net's distance from the definition is not certified finite",
            ))
        }
    }
}

/// The path's derivative bounds over a span (module docs): `T = N / w` for the homogeneous
/// numerator `N = w P` and the weight `w`, so
/// `|T'| <= (|N'| + |T| |w'|) / w`, `|T''| <= (|N''| + 2 D1 |w'| + |T| |w''|) / w` and
/// `|T'''| <= (|N'''| + 3 D2 |w'| + 3 D1 |w''| + |T| |w'''|) / w`, each of `N^(k)`, `w^(k)`
/// at most the largest of its control points over the span (convex hull), `w` at least the
/// smallest weight, `|T|` at most the largest `|P|`.
struct PathBounds {
    n: [Spline; 4],
    w: [Spline; 4],
    p: Spline,
    unique_knots: Vec<f64>,
}

impl PathBounds {
    fn new(path: &BSplineCurve) -> Result<Self> {
        let pc = path.canonical()?;
        let knots = full_knots(&pc.knots, &pc.knot_multiplicities);
        let weight = |i: usize| pc.weights.get(i).copied().unwrap_or(1.0);
        if (0..pc.control_points.len()).any(|i| !(weight(i) > 0.0 && weight(i).is_finite())) {
            return Err(bad("a path's weights are positive"));
        }
        let numerator: Vec<Vec3> = pc
            .control_points
            .iter()
            .enumerate()
            .map(|(i, &p)| p * weight(i))
            .collect();
        let weights: Vec<Vec3> = (0..pc.control_points.len())
            .map(|i| Vec3::new(weight(i), 0.0, 0.0))
            .collect();
        let chain = |s: Spline| {
            let s1 = s.derivative();
            let s2 = s1.derivative();
            let s3 = s2.derivative();
            [s, s1, s2, s3]
        };
        Ok(Self {
            n: chain(Spline::new(pc.degree, knots.clone(), numerator)?),
            w: chain(Spline::new(pc.degree, knots.clone(), weights)?),
            p: Spline::new(pc.degree, knots, pc.control_points)?,
            unique_knots: pc.knots,
        })
    }

    /// Upper bounds of `|T'|`, `|T''|` and `|T'''|` over the span holding `t`.
    fn derivative_bounds(&self, t: f64) -> (f64, f64, f64) {
        let wmin = self.w.first().map_or(0.0, |w| {
            w.local(t).iter().fold(f64::INFINITY, |m, p| m.min(p.x))
        });
        let tmax = self.p.largest(t);
        let nk = |k: usize| self.n.get(k).map_or(f64::INFINITY, |s| s.largest(t));
        let wk = |k: usize| self.w.get(k).map_or(f64::INFINITY, |s| s.largest(t));
        if !(wmin > 0.0) {
            return (f64::INFINITY, f64::INFINITY, f64::INFINITY);
        }
        let d1 = (nk(1) + tmax * wk(1)) / wmin;
        let d2 = (nk(2) + 2.0 * d1 * wk(1) + tmax * wk(2)) / wmin;
        let d3 = (nk(3) + 3.0 * d2 * wk(1) + 3.0 * d1 * wk(2) + tmax * wk(3)) / wmin;
        (d1, d2, d3)
    }
}

/// The net's control points grouped as the sweep's columns: one spline in `v` per control
/// point of the profile.
struct Columns {
    splines: Vec<Spline>,
    /// The profile's control points, in the frame's coordinates.
    profile: Vec<Vec3>,
    /// Their weights (1 for a polynomial profile).
    weights: Vec<f64>,
    /// Whether the profile is rational.
    rational: bool,
    unique_knots: Vec<f64>,
    domain: (f64, f64),
    swapped: bool,
}

impl Columns {
    /// Reads `net` as the profile's control points swept along one parameter, in the same
    /// rational basis (module docs).
    fn read(profile: &BSplineCurve, net: &BSplineSurface) -> Result<Self> {
        let pc = profile.canonical()?;
        if pc.dim != 3 {
            return Err(bad("a profile is a 3D curve"));
        }
        let sc = net.canonical()?;
        let grid = net.control_point_net()?;
        let (nu, nv) = (grid.u_count, grid.v_count);
        let sections = pc.control_points.len();
        let same = |a: &[f64], b: &[f64]| {
            a.len() == b.len()
                && a.iter()
                    .zip(b)
                    .all(|(x, y)| (x - y).abs() <= 1e-12 * (1.0 + x.abs()))
        };
        let (swapped, along_knots, along_mult, along_degree, along) = if nu == sections
            && sc.u_degree == pc.degree
            && same(&sc.u_knots, &pc.knots)
            && sc.u_multiplicities == pc.knot_multiplicities
        {
            (false, &sc.v_knots, &sc.v_multiplicities, sc.v_degree, nv)
        } else if nv == sections
            && sc.v_degree == pc.degree
            && same(&sc.v_knots, &pc.knots)
            && sc.v_multiplicities == pc.knot_multiplicities
        {
            (true, &sc.u_knots, &sc.u_multiplicities, sc.u_degree, nu)
        } else {
            return Err(bad("a net's sections are not the profile's control points"));
        };
        let at = |i: usize, j: usize| {
            let idx = if swapped { j * nv + i } else { i * nv + j };
            (
                grid.points.get(idx).copied(),
                grid.weights.get(idx).copied(),
            )
        };
        let pw = |i: usize| pc.weights.get(i).copied().unwrap_or(1.0);
        let full = full_knots(along_knots, along_mult);
        let mut splines = Vec::with_capacity(sections);
        let mut ratio: Option<f64> = None;
        for i in 0..sections {
            let mut ctrl = Vec::with_capacity(along);
            for j in 0..along {
                let (Some(p), w) = at(i, j) else {
                    return Err(bad("a net's control points are missing"));
                };
                // the basis in v is the same for every column: a weight per column, within
                // the profile's up to one common factor
                let w = w.unwrap_or(1.0);
                let r = w / pw(i);
                let r0 = *ratio.get_or_insert(r);
                if !((r - r0).abs() <= 1e-12 * r0.abs()) {
                    return Err(bad("a net's weights are not the profile's"));
                }
                ctrl.push(p);
            }
            splines.push(Spline::new(along_degree, full.clone(), ctrl)?);
        }
        let domain = splines
            .first()
            .map(Spline::domain)
            .ok_or_else(|| bad("a net has sections"))?;
        Ok(Self {
            splines,
            weights: (0..pc.control_points.len()).map(pw).collect(),
            rational: !pc.weights.is_empty(),
            profile: pc.control_points,
            unique_knots: along_knots.clone(),
            domain,
            swapped,
        })
    }
}
