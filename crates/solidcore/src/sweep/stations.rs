//! A sweep's stations from its laws: where along its path, and in what frame, the profile is
//! placed, for a path of lines, circular arcs and helices, with the profile's orientation,
//! twist, end scale and up to two guide rails.
//!
//! Upstream: none (not in SMLib). The laws are SolidCraft's (`crates/kernel/src/sweep_path.rs`
//! at 52afe21, `SweepControl`, `controlled_stations`, `rail_frame`), moved into the kernel so
//! no application computes a frame.
//!
//! | Rust              | Upstream |
//! |-------------------|----------|
//! | [`Path`], [`Piece`] | (none) |
//! | [`Laws`], [`Orientation`], [`RailFit`] | (none; SolidCraft's `SweepControl`) |
//! | [`Station`]       | (none; SolidCraft's `Station`) |
//! | [`stations`]      | (none; SolidCraft's `controlled_stations`) |
//!
//! **A station** is a point `at` of the path and a linear map `m` (its columns the images of
//! the axes): a point `q` of the profile at the first station is placed at
//! `at + m (q - at_0)`, `at_0` the first station's point.
//!
//! **The frame (exact).** Along the path the profile is carried by the path's
//! rotation-minimizing frame (parallel transport), which never turns about the tangent. The
//! map from the first station's frame to a station's is a rotation `R`, independent of the
//! frame the first station starts in, and in closed form on every piece:
//! - a **line**: the identity;
//! - a **circular arc** about the unit axis `n`: the rotation about `n` by the angle run, as
//!   the arc's tangent turns (a planar curve's rotation-minimizing frame is its Frenet frame's
//!   turn in its plane);
//! - a **helix** of radius `a` and rise `b` per radian about `n`: the screw's rotation about
//!   `n` by the angle run `phi`, then about the new tangent by `-tau s = -b phi / sqrt(a^2 +
//!   b^2)`. The Frenet frame of a helix is carried by the screw; the rotation-minimizing normal
//!   `U = cos psi N + sin psi B` keeps `U'` along the tangent where `psi' = -tau` (Frenet:
//!   `N' = -kappa T + tau B`), `tau = b / (a^2 + b^2)` and `s = phi sqrt(a^2 + b^2)` (Wang,
//!   Juettler, Zheng and Liu, *Computation of rotation minimizing frames*, ACM TOG 27, 2008, 2;
//!   Guggenheimer, *Computing frames along a trajectory*, CAGD 6, 1989);
//! - a **corner** between pieces, tangents `t_0` to `t_1`: the least rotation taking one to the
//!   other, about `t_0 x t_1` by the angle between them (the limit of the double reflection
//!   across a vanishing segment); none where the pieces meet G1.
//!
//! Each is a rigid rotation computed by Rodrigues' formula, so a frame is orthonormal to
//! rounding: its defect, the largest of `| |c_i|^2 - 1 |` and `|c_i . c_j|` over its columns,
//! is `O(k eps)` for `k` composed rotations.
//!
//! **The laws (SolidCraft's conventions, carried over, not derived here).**
//! - `twist` turns the profile about its normal (turned to point the way the path leaves it)
//!   by `twist u` radians, right-handed, and `end_scale` scales it about the path by
//!   `1 + (end_scale - 1) u`, `u` the fraction of the path's arc length run.
//! - A station at a corner (its fraction exactly one of [`Path::corner_fractions`], which
//!   decides it: a fraction is placed on a piece by comparing it with the pieces' start
//!   fractions, start over total, exactly) is the mitre: the rotation by half
//!   the corner's angle, stretched by `1 / cos(theta / 2)` along `t_1 - t_0`.
//! - [`Orientation::Parallel`] keeps the profile as drawn (the identity), and needs the path
//!   to cross the profile's plane by 5 degrees or more throughout;
//!   [`Orientation::Direction`] keeps the direction's projection square to the tangent, and
//!   needs the path to stay more than `1e-3` (about 0.06 degrees) away from it.
//! - Guide rails take the default orientation, no twist and scale 1, and a path without
//!   corners. The station's plane, square to the path, meets each rail once; with one rail the
//!   profile turns about the path to follow the line from the path to the rail and is scaled
//!   (or stretched across, [`RailFit::Stretch`]) by its length against the start's; with two,
//!   it moves with the first, turns with the line between them and scales with their distance.
//!
//! **Rails.** A rail is a kernel curve with a NURBS form (a line, an arc or circle, a
//! B-spline; a polyline is a degree-1 B-spline). The station's plane meets it at the roots of
//! the plane's signed distance along it, `(X(t) - p W(t)) . t_path` in its homogeneous form,
//! isolated on its Bezier spans by subdividing their Bernstein coefficients
//! ([`crate::query::tangent_edge`]'s `bernstein_roots`). A plane the rail meets other than
//! exactly once (not at all, or more than once) is refused, naming it.
//!
//! **Refused.** A path piece that is not a line, arc or helix: a spline path's
//! rotation-minimizing frame needs a certified enclosure (tangency-labs/solidcore#2572).
//!
//! Deviations: no upstream counterpart. At most `MAX_PIECES` pieces, `MAX_STATIONS`
//! stations, `MAX_TURNS` turns of an arc or helix; the root isolation's bounds are
//! `bernstein_roots`'.

use crate::curve::curve::Curve;
use crate::error::{Result, SolidError};
use crate::kernel::Kernel;
use crate::math::{Point3, Vec3};
use crate::nurbs::crv_basic::Curve as NurbCurve;
use crate::nurbs::math::{CPoint, NOW};
use crate::query::tangent_edge::{bernstein_roots, MAX_DEPTH, MAX_SUBDIVISIONS};

const OP: &str = "sweep::stations";

/// Most path pieces.
pub const MAX_PIECES: usize = 1000;
/// Most stations asked at once.
pub const MAX_STATIONS: usize = 100_000;
/// Most turns of one arc or helix piece.
pub const MAX_TURNS: f64 = 1000.0;
/// The most twist accepted: four turns (radians; SolidCraft's).
pub const MAX_TWIST: f64 = 4.0 * std::f64::consts::TAU;
/// The smallest end scale accepted (SolidCraft's).
pub const MIN_SCALE: f64 = 0.01;
/// The largest end scale accepted (SolidCraft's).
pub const MAX_SCALE: f64 = 100.0;

fn invalid(reason: &'static str) -> SolidError {
    SolidError::invalid_input(OP, reason)
}

fn unsupported(reason: &'static str) -> SolidError {
    SolidError::Unsupported {
        operation: OP,
        reason,
    }
}

/// `v` turned about the unit `axis` by `angle` (Rodrigues).
fn rotate(axis: Vec3, angle: f64, v: Vec3) -> Vec3 {
    let (s, c) = angle.sin_cos();
    v * c + axis.cross(v) * s + axis * (axis.dot(v) * (1.0 - c))
}

/// A linear map by its columns.
type Map = [Vec3; 3];

const IDENTITY: Map = [Vec3::X, Vec3::Y, Vec3::Z];

fn apply(m: &Map, v: Vec3) -> Vec3 {
    let [a, b, c] = *m;
    a * v.x + b * v.y + c * v.z
}

/// `a` after `b`.
fn compose(a: &Map, b: &Map) -> Map {
    b.map(|col| apply(a, col))
}

fn rotation(axis: Vec3, angle: f64) -> Map {
    IDENTITY.map(|e| rotate(axis, angle, e))
}

fn finite(v: Vec3) -> bool {
    v.x.is_finite() && v.y.is_finite() && v.z.is_finite()
}

fn unit(v: Vec3, why: &'static str) -> Result<Vec3> {
    v.unitize().map_err(|_| invalid(why))
}

/// One piece of a sweep's path.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Piece {
    /// From `a` to `b`.
    Line {
        /// Its start.
        a: Point3,
        /// Its end.
        b: Point3,
    },
    /// From `start`, `angle` radians (positive) right-handed about `axis` (unit) through
    /// `centre`.
    Arc {
        /// A point of the axis.
        centre: Point3,
        /// The axis (unit).
        axis: Vec3,
        /// Its start, off the axis.
        start: Point3,
        /// The angle run (radians, positive).
        angle: f64,
    },
    /// From `start`, `angle` radians (positive) right-handed about `axis` (unit) through
    /// `centre`, rising `rise` along the axis per radian (negative: falling).
    Helix {
        /// A point of the axis.
        centre: Point3,
        /// The axis (unit).
        axis: Vec3,
        /// Its start, off the axis.
        start: Point3,
        /// The angle run (radians, positive).
        angle: f64,
        /// The rise per radian.
        rise: f64,
    },
}

impl Piece {
    /// The piece's radius vector at its start (square to the axis), its unit axis and rise.
    fn screw(&self) -> Option<(Point3, Vec3, Vec3, f64, f64)> {
        let (centre, axis, start, angle, rise) = match *self {
            Piece::Line { .. } => return None,
            Piece::Arc {
                centre,
                axis,
                start,
                angle,
            } => (centre, axis, start, angle, 0.0),
            Piece::Helix {
                centre,
                axis,
                start,
                angle,
                rise,
            } => (centre, axis, start, angle, rise),
        };
        let d = start - centre;
        let radial = d - axis * axis.dot(d);
        let foot = start - radial;
        Some((foot, axis, radial, angle, rise))
    }

    fn check(&self) -> Result<()> {
        match *self {
            Piece::Line { a, b } => {
                let ends_finite = finite(a) && finite(b);
                if !ends_finite || a.distance(b) <= 1e-9 * (1.0 + a.length()) {
                    return Err(invalid("a line piece with no length or not finite"));
                }
            }
            Piece::Arc { axis, angle, .. } | Piece::Helix { axis, angle, .. } => {
                let unit_axis = (axis.length() - 1.0).abs() <= 1e-12;
                if !unit_axis {
                    return Err(invalid("a piece's axis that is not a unit vector"));
                }
                let turns = angle / std::f64::consts::TAU;
                if !(angle > 0.0 && turns <= MAX_TURNS) {
                    return Err(invalid(
                        "a piece's angle not positive or past the most turns",
                    ));
                }
                let Some((_, _, radial, _, rise)) = self.screw() else {
                    return Err(invalid("a piece's axis"));
                };
                let size = radial.length();
                if !(size > 0.0 && size.is_finite() && rise.is_finite()) {
                    return Err(invalid("a piece's start on its axis or not finite"));
                }
            }
        }
        Ok(())
    }

    /// The piece's arc length.
    fn length(&self) -> f64 {
        if let Piece::Line { a, b } = *self {
            a.distance(b)
        } else {
            self.screw().map_or(0.0, |(_, _, radial, angle, rise)| {
                angle * radial.length().hypot(rise)
            })
        }
    }

    /// The point, unit tangent and the rotation-minimizing transport from the piece's start, at
    /// the fraction `f` of the piece's length.
    fn at(&self, f: f64) -> Result<(Point3, Vec3, Map)> {
        if let Piece::Line { a, b } = *self {
            let t = unit(b - a, "a line piece with no length")?;
            return Ok((a + (b - a) * f, t, IDENTITY));
        }
        let Some((foot, axis, radial, angle, rise)) = self.screw() else {
            return Err(invalid("a piece's axis"));
        };
        let phi = angle * f;
        let r = rotate(axis, phi, radial);
        let p = foot + r + axis * (rise * phi);
        let t = unit(axis.cross(r) + axis * rise, "a piece's tangent")?;
        let screw = rotation(axis, phi);
        let a = radial.length();
        let psi = -rise * phi / a.hypot(rise);
        let m = if rise == 0.0 {
            screw
        } else {
            compose(&rotation(t, psi), &screw)
        };
        Ok((p, t, m))
    }

    /// The unit tangents at the piece's start and end.
    fn ends(&self) -> Result<(Point3, Vec3, Point3, Vec3)> {
        let (p0, t0, _) = self.at(0.0)?;
        let (p1, t1, _) = self.at(1.0)?;
        Ok((p0, t0, p1, t1))
    }
}

/// A sweep's path: pieces end to end, each starting where the last ended.
#[derive(Debug, Clone, PartialEq)]
pub struct Path {
    /// The pieces, in order.
    pub pieces: Vec<Piece>,
}

/// A corner between two pieces: the transport across it and the mitre.
#[derive(Clone, Copy)]
struct Corner {
    /// The least rotation from the arriving tangent to the leaving one.
    turn: Map,
    /// Half of it.
    half: Map,
    /// The mitre's stretch direction and factor (SolidCraft's).
    mitre: Option<(Vec3, f64)>,
}

/// A path measured: each piece's start and length along the whole, and the corner before it.
struct Measured {
    pieces: Vec<(Piece, f64, f64, Corner)>,
    total: f64,
}

impl Measured {
    /// The fraction of the path's length at which piece `j` starts: its start / the total, the
    /// one expression [`Path::joints`] and [`Path::corner_fractions`] give the caller, so a
    /// station asked at that value is compared with it exactly.
    fn fraction(&self, j: usize) -> Option<f64> {
        self.pieces
            .get(j)
            .map(|(_, start, _, _)| start / self.total)
    }
}

/// Below this angle (radians) two tangents meet G1: no corner.
const G1: f64 = 1e-12;

fn corner(t0: Vec3, t1: Vec3) -> Result<Corner> {
    let angle = t0.cross(t1).length().atan2(t0.dot(t1));
    if angle <= G1 {
        return Ok(Corner {
            turn: IDENTITY,
            half: IDENTITY,
            mitre: None,
        });
    }
    if angle >= std::f64::consts::PI - 1e-9 {
        return Err(invalid("a path that doubles back at a corner"));
    }
    let axis = unit(t0.cross(t1), "a corner's axis")?;
    let stretch = unit(t1 - t0, "a corner's mitre")?;
    Ok(Corner {
        turn: rotation(axis, angle),
        half: rotation(axis, 0.5 * angle),
        mitre: Some((stretch, 1.0 / (0.5 * angle).cos())),
    })
}

impl Path {
    fn measure(&self) -> Result<Measured> {
        if self.pieces.is_empty() || self.pieces.len() > MAX_PIECES {
            return Err(invalid("a path of 1 to MAX_PIECES pieces"));
        }
        let mut out = Vec::new();
        let mut run = 0.0;
        let mut last: Option<(Point3, Vec3)> = None;
        for piece in &self.pieces {
            piece.check()?;
            let (p0, t0, p1, t1) = piece.ends()?;
            let c = match last {
                None => corner(t0, t0)?,
                Some((end, t_end)) => {
                    if end.distance(p0) > 1e-9 * (1.0 + end.length()) {
                        return Err(invalid("a path whose pieces do not meet end to start"));
                    }
                    corner(t_end, t0)?
                }
            };
            let len = piece.length();
            out.push((*piece, run, len, c));
            run += len;
            last = Some((p1, t1));
        }
        if !(run > 0.0 && run.is_finite()) {
            return Err(invalid("a path with no length"));
        }
        Ok(Measured {
            pieces: out,
            total: run,
        })
    }

    /// The fractions of the path's arc length at which its pieces meet.
    ///
    /// # Errors
    ///
    /// As [`stations`] for the path.
    pub fn joints(&self) -> Result<Vec<f64>> {
        let m = self.measure()?;
        Ok((1..m.pieces.len()).filter_map(|j| m.fraction(j)).collect())
    }

    /// The fractions of the path's arc length at its corners (joints where the pieces do not
    /// meet G1): a station asked at exactly one of these values is that corner's mitre.
    ///
    /// # Errors
    ///
    /// As [`stations`] for the path.
    pub fn corner_fractions(&self) -> Result<Vec<f64>> {
        let m = self.measure()?;
        Ok((1..m.pieces.len())
            .filter(|&j| {
                m.pieces
                    .get(j)
                    .is_some_and(|(_, _, _, c)| c.mitre.is_some())
            })
            .filter_map(|j| m.fraction(j))
            .collect())
    }
}

/// How the profile turns as it follows the path.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub enum Orientation {
    /// Square to the path, carried by its rotation-minimizing frame.
    #[default]
    Perpendicular,
    /// As drawn: the profile only moves along the path.
    Parallel,
    /// Square to the path, its turn about the path fixed by a direction: the direction's
    /// projection square to the tangent.
    Direction(Vec3),
}

/// How guide rails size the profile.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum RailFit {
    /// Evenly, by the distance across.
    #[default]
    Scale,
    /// Only across to the rail (between the rails).
    Stretch,
}

/// A sweep's laws (SolidCraft's conventions; module docs). The default is a plain sweep.
#[derive(Debug, Clone)]
pub struct Laws<'a> {
    /// How the profile turns.
    pub orientation: Orientation,
    /// The twist over the path (radians).
    pub twist: f64,
    /// The profile's scale at the end.
    pub end_scale: f64,
    /// Up to two guide rails.
    pub rails: Vec<&'a Curve>,
    /// How the rails size the profile.
    pub rail_fit: RailFit,
}

impl Default for Laws<'_> {
    fn default() -> Self {
        Laws {
            orientation: Orientation::Perpendicular,
            twist: 0.0,
            end_scale: 1.0,
            rails: Vec::new(),
            rail_fit: RailFit::Scale,
        }
    }
}

impl Laws<'_> {
    fn check(&self) -> Result<()> {
        if !(self.twist.is_finite() && self.twist.abs() <= MAX_TWIST) {
            return Err(invalid("a twist past four turns"));
        }
        if !(MIN_SCALE..=MAX_SCALE).contains(&self.end_scale) {
            return Err(invalid("an end scale outside 0.01 to 100"));
        }
        if let Orientation::Direction(d) = self.orientation {
            if !(finite(d) && d.length() > 1e-9 && d.length() < 1e9) {
                return Err(invalid(
                    "an orientation direction not finite or of no length",
                ));
            }
        }
        if self.rails.len() > 2 {
            return Err(invalid("more than two guide rails"));
        }
        let defaults = self.orientation == Orientation::Perpendicular
            && self.twist == 0.0
            && self.end_scale == 1.0;
        if !self.rails.is_empty() && !defaults {
            return Err(invalid(
                "guide rails with an orientation, twist or end scale: the rails set them",
            ));
        }
        Ok(())
    }
}

/// A station: a point of the path and the linear map placing the profile there (module
/// docs).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Station {
    /// The point.
    pub at: Point3,
    /// The map's columns, the images of the axes.
    pub m: [Vec3; 3],
}

/// Where along the measured path the fraction `u` is: the point, unit tangent, transport from
/// the path's start, and the corner when `u` is exactly at one.
fn locate(m: &Measured, u: f64) -> Result<(Point3, Vec3, Map, Option<Corner>)> {
    // which piece, and whether a corner, is decided on the fraction itself, compared with
    // each piece's start fraction exactly (no window: [`Path::corner_fractions`])
    let mut transport = IDENTITY;
    let count = m.pieces.len();
    for (i, (piece, start, len, c)) in m.pieces.iter().enumerate() {
        let before = compose(&c.turn, &transport);
        if i > 0 && c.mitre.is_some() && m.fraction(i) == Some(u) {
            let (p, t0_piece, _) = piece.at(0.0)?;
            let half = compose(&c.half, &transport);
            let t_prev = apply(&transpose_rotation(&c.turn), t0_piece);
            let bisector = unit(t_prev + t0_piece, "a corner's bisector")?;
            return Ok((p, bisector, half, Some(*c)));
        }
        let last = i + 1 == count;
        let before_next = m.fraction(i + 1).is_some_and(|next| u < next);
        if before_next || last {
            let s = u * m.total;
            let f = if *len > 0.0 {
                ((s - start) / len).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let (p, t, local) = piece.at(f)?;
            return Ok((p, t, compose(&local, &before), None));
        }
        let (_, _, whole) = piece.at(1.0)?;
        transport = compose(&whole, &before);
    }
    Err(invalid("a path without pieces"))
}

/// The inverse of a rotation: its transpose.
fn transpose_rotation(r: &Map) -> Map {
    let [a, b, c] = *r;
    [
        Vec3::new(a.x, b.x, c.x),
        Vec3::new(a.y, b.y, c.y),
        Vec3::new(a.z, b.z, c.z),
    ]
}

/// The roots in its parameter of the signed distance of `rail` from the plane through `p`
/// square to the unit `t`; refused unless exactly one.
fn rail_crossing(k: &Kernel, rail: &Curve, p: Point3, t: Vec3) -> Result<Point3> {
    let base = rail
        .b_spline_base()
        .ok_or_else(|| unsupported("a guide rail with no NURBS form"))?;
    let r = base.degree()?;
    let (unique, mult) = base.knots(None)?;
    let mut knots = Vec::new();
    for (&x, &n) in unique.iter().zip(&mult) {
        knots.extend(std::iter::repeat_n(x, n));
    }
    let (points, weights) = base.control_polygon()?;
    if !weights.is_empty() && weights.len() != points.len() {
        return Err(invalid("a rail's weights do not match its control points"));
    }
    let mut pw = Vec::new();
    for (i, q) in points.iter().enumerate() {
        if weights.is_empty() {
            pw.push(CPoint::new(q.x, q.y, q.z, NOW));
            continue;
        }
        let Some(&w) = weights.get(i) else {
            return Err(invalid("a rail's weights do not match its control points"));
        };
        let positive = w > 0.0;
        if !positive {
            return Err(invalid("a rail's weight that is not positive"));
        }
        pw.push(CPoint::new(q.x * w, q.y * w, q.z * w, w));
    }
    let nurb = NurbCurve::new(pw, r, knots)?;
    let mut roots: Vec<f64> = Vec::new();
    for span in nurb.decompose_bez()? {
        let (Some(&a), Some(&b)) = (span.knots().first(), span.knots().last()) else {
            continue;
        };
        let w = |c: &CPoint| if c.w == NOW { 1.0 } else { c.w };
        let dist = |c: &CPoint| {
            (c.x - p.x * w(c)) * t.x + (c.y - p.y * w(c)) * t.y + (c.z - p.z * w(c)) * t.z
        };
        let scale = span
            .cpts()
            .iter()
            .map(|c| c.x.abs() + c.y.abs() + c.z.abs() + p.length() * w(c))
            .fold(0.0, f64::max);
        let floor = 64.0 * f64::EPSILON * scale;
        let c: Vec<f64> = span
            .cpts()
            .iter()
            .map(|c| {
                let d = dist(c);
                if d.abs() <= floor {
                    0.0
                } else {
                    d
                }
            })
            .collect();
        if c.iter().all(|x| *x == 0.0) {
            return Err(invalid("a guide rail lying in a station's plane"));
        }
        let mut found = Vec::new();
        let mut work = MAX_SUBDIVISIONS;
        bernstein_roots(&c, 0.0, 1.0, MAX_DEPTH, &mut work, &mut found);
        for s in found {
            let x = a + (b - a) * s;
            // a root at a span's end is the next span's start: one crossing
            if !roots
                .iter()
                .any(|y| (y - x).abs() <= 1e-12 * (1.0 + x.abs()))
            {
                roots.push(x);
            }
        }
    }
    match roots.as_slice() {
        [x] => rail.evaluate_point(k, *x),
        [] => Err(invalid(
            "a guide rail that does not reach a station's plane: it must run the path's length",
        )),
        _ => Err(invalid(
            "a guide rail that crosses a station's plane more than once",
        )),
    }
}

/// How guide rails turn and size a section (SolidCraft's `rail_frame`): `along` (unit) turns
/// onto `along_to`, `across` (unit, square to `along`) onto `across_to` sized by `k`, and the
/// third way follows, sized by `k` ([`RailFit::Scale`]) or kept ([`RailFit::Stretch`]).
fn rail_frame(
    along: Vec3,
    across: Vec3,
    along_to: Vec3,
    across_to: Vec3,
    k: f64,
    fit: RailFit,
) -> Map {
    let k_other = match fit {
        RailFit::Scale => k,
        RailFit::Stretch => 1.0,
    };
    let (w0, w1) = (along.cross(across), along_to.cross(across_to));
    let map = |v: Vec3| {
        along_to * v.dot(along) + across_to * (k * v.dot(across)) + w1 * (k_other * v.dot(w0))
    };
    IDENTITY.map(map)
}

/// The stations of a sweep along `path` at the fractions `at` of its arc length, placing a
/// profile whose plane's normal is `profile_normal`, by `laws` (module docs).
///
/// # Errors
///
/// [`SolidError::InvalidInput`] for a path or laws out of range (pieces not meeting, a corner
/// doubling back, a parallel profile the path does not cross, a path along the orientation
/// direction, a rail not crossing a station's plane exactly once, rails along a path with
/// corners, rails meeting each other or the path), fractions outside `[0, 1]` or past
/// `MAX_STATIONS`; [`SolidError::Unsupported`] for a rail with no NURBS form; those of the
/// rails' evaluations.
pub fn stations(
    k: &Kernel,
    path: &Path,
    profile_normal: Vec3,
    laws: &Laws<'_>,
    at: &[f64],
) -> Result<Vec<Station>> {
    laws.check()?;
    if at.is_empty() || at.len() > MAX_STATIONS {
        return Err(invalid("1 to MAX_STATIONS stations"));
    }
    if at.iter().any(|u| !(0.0..=1.0).contains(u)) {
        return Err(invalid("a station's fraction outside [0, 1]"));
    }
    let measured = path.measure()?;
    let n = unit(profile_normal, "a profile normal of no length")?;
    let (p0, t0, _, _) = locate(&measured, 0.0)?;
    let axis = if n.dot(t0) < 0.0 { n * -1.0 } else { n };
    if !laws.rails.is_empty() {
        return railed(k, &measured, laws, at, p0, t0);
    }
    let mut out = Vec::new();
    for &u in at {
        let (p, t, transport, at_corner) = locate(&measured, u)?;
        let turn = match laws.orientation {
            Orientation::Perpendicular => transport,
            Orientation::Parallel => {
                if t.dot(axis) < 5.0_f64.to_radians().sin() {
                    return Err(invalid(
                        "a parallel profile the path does not cross by 5 degrees or more",
                    ));
                }
                IDENTITY
            }
            Orientation::Direction(d) => {
                let frame = |t: Vec3| -> Result<Map> {
                    let d = unit(d, "an orientation direction of no length")?;
                    let across = d - t * t.dot(d);
                    if across.length() < 1e-3 {
                        return Err(invalid("a path running along the orientation direction"));
                    }
                    let r = unit(across, "a path running along the orientation direction")?;
                    Ok([t, r, t.cross(r)])
                };
                let (f0, fi) = (frame(t0)?, frame(t)?);
                compose(&fi, &transpose_rotation(&f0))
            }
        };
        let (sin, cos) = (laws.twist * u).sin_cos();
        let scale = 1.0 + (laws.end_scale - 1.0) * u;
        let shape = |v: Vec3| {
            let along = axis * v.dot(axis);
            let r = (v - along) * scale;
            along + r * cos + axis.cross(r) * sin
        };
        let mut m = IDENTITY.map(|e| apply(&turn, shape(e)));
        if laws.orientation != Orientation::Parallel {
            if let Some(Corner {
                mitre: Some((b, s)),
                ..
            }) = at_corner
            {
                for c in &mut m {
                    *c += b * ((s - 1.0) * c.dot(b));
                }
            }
        }
        out.push(Station { at: p, m });
    }
    Ok(out)
}

/// The stations shaped by guide rails (module docs).
fn railed(
    k: &Kernel,
    measured: &Measured,
    laws: &Laws<'_>,
    at: &[f64],
    s0: Point3,
    t0: Vec3,
) -> Result<Vec<Station>> {
    if measured.pieces.iter().any(|(_, _, _, c)| c.mitre.is_some()) {
        return Err(unsupported("guide rails along a path with corners"));
    }
    let near = |d: Vec3| d.length() <= 1e-9 * (1.0 + measured.total);
    let across_at = |u: f64| -> Result<(Point3, Vec3, Point3, Vec3)> {
        let (p, t, _, _) = locate(measured, u)?;
        let mut hits = Vec::new();
        for rail in &laws.rails {
            hits.push(rail_crossing(k, rail, p, t)?);
        }
        let (anchor, v) = match hits.as_slice() {
            [r] => (p, *r - p),
            [a, b] => (*a, *b - *a),
            _ => return Err(invalid("one or two guide rails")),
        };
        Ok((p, t, anchor, v - t * v.dot(t)))
    };
    let (_, _, a0, d0) = across_at(0.0)?;
    if near(d0) {
        return Err(invalid(
            "a guide rail starting on the path, or two starting at one point",
        ));
    }
    let u0 = unit(d0, "a guide rail's start")?;
    let mut out = Vec::new();
    for &u in at {
        let (_, t, anchor, d) = across_at(u)?;
        if near(d) {
            return Err(invalid(
                "a guide rail meeting the path, or two meeting: the profile shrinks to nothing",
            ));
        }
        let ud = unit(d, "a guide rail's direction across")?;
        let m = rail_frame(t0, u0, t, ud, d.length() / d0.length(), laws.rail_fit);
        // the anchor's start lands on the anchor
        out.push(Station {
            at: anchor - apply(&m, a0 - s0),
            m,
        });
    }
    Ok(out)
}
