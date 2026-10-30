//! The sweeps: extrusion (straight or with draft), revolution, sweep along a path, loft
//! through sections or to a point, the sweep through stations of a profile that twists,
//! scales or follows rails, and planar faces from loops (the screw sweeps follow with
//! `solidcore::sweep::screw`).
//! Profiles, paths and sections cross as curve lists (`kinds`, `reals`: lines, circular
//! arcs and NURBS); each call builds its own curves from them.
//!
//! Source: tangency-labs/smlib `capi/src/smlib_sweeps.cpp` (Apache-2.0, Copyright SolidCraft
//! contributors), replay/1.0.1 at `a3485b96c` (the piecewise sweep along lines and arcs), and
//! the `SMLIB_CURVE_*`, `SMLIB_DRAFT_CORNER_*`, `SMLIB_SWEEP_*`, `SMLIB_ENDS_*` and
//! `SMLIB_SWEEP_MAX_*` values of `capi/include/smlib.h`.
//!
//! | Rust | C API (`SM_API` / SMLib call) |
//! |------|-------|
//! | [`Smlib::extrude`] | `smlib_extrude` (`SmApiCreateLinearSweep`, `SmApiCreateTaperExtrude`) |
//! | [`Smlib::revolve`] | `smlib_revolve` (`SmApiCreateRotationalSweep`) |
//! | [`Smlib::sweep`], [`Smlib::sweep_placed`] | `smlib_sweep`, `smlib_sweep_placed` (`SmApiCreateCurveSweep`, or the linear / rotational sweep) |
//! | [`Smlib::loft`] | `smlib_loft` (`SmApiCreateSkinPrimitive`) |
//! | [`Smlib::loft_to_point`] | `smlib_loft_to_point` (`CreatePlanarFaceWith3DCurves`, `StitchAndOrient`) |
//! | [`Smlib::sweep_stations`] | `smlib_sweep_stations` (`SmApiCreateSkinPrimitive`, `SmApiStitchIntoSolid`) |
//! | [`Smlib::planar_faces`] | `smlib_planar_faces` (`SmApiCreatePlanarFaces`) |
//! | `parse`, `line`, `arc`, `nurbs`, `reals_of` (private) | `parse`, `line`, `arc`, `nurbs`, `reals_of` |
//! | `one_loop`, `straight`, `circular`, `joined` (private) | the same (`finish`: `crate::session`'s) |
//! | `Rigid`, `Segment`, `segmented`, `carried`, `piecewise` (private) | the same |
//! | `put_motion`, `turned`, `maps_ok`, `placed` (private) | the same |
//!
//! Deviations:
//! - **Arguments.** A curve list is two slices (`kinds`, `reals`), a station map list a
//!   slice of `[f64; 12]`, a radial law a slice of `[turn, offset]` pairs; their lengths are
//!   the C counts. No out-params: `smlib_sweep_placed`'s `placed` comes with the body, and
//!   `smlib_extrude`'s `out_corner` is [`Smlib::last_failed_corner_point`].
//! - **Logging.** The C API's replay dumps of curve lists (`dump_curves`) are not written:
//!   the dump directory holds the call's line and its B-rep inputs only (these calls have
//!   none).

use std::f64::consts::{PI, TAU};

use solidcore::curve::curve::p3::compute_crvs_normal;
use solidcore::curve::curve::{Curve, CurveKind};
use solidcore::curve::ellipse::Ellipse;
use solidcore::curve::line::Line;
use solidcore::geom::axis2_placement::Axis2Placement;
use solidcore::geom::extent1d::Extent1d;
use solidcore::geom::tol::USE_DEFAULT;
use solidcore::heal::stitch::p1::api_stitch_into_solid;
use solidcore::math::{Point3, Vec3};
use solidcore::nurbs::b_spline_curve::{BSplineCurve, BSplineCurveForm, CanonicalData, KnotType};
use solidcore::primitive::api_trim_ops::create_planar_faces;
use solidcore::primitive::primitive_creation::SweepOptions;
use solidcore::query::brep::BrepIdentity;
use solidcore::sweep::api_brep::{
    create_curve_sweep, create_curve_sweep_with_placements, create_linear_sweep,
    create_rotational_sweep, create_taper_extrude_with_corner, OffsetCornerType,
    TaperCornerRefusal, TaperExtrude,
};
use solidcore::sweep::api_surfaces::create_skin_primitive;
use solidcore::topology::brep::p2::{caller_cache_key, ArenaNotifyState};
use solidcore::topology::brep::p3::new_brep;
use solidcore::topology::brep::p4::stitch_faces;
use solidcore::topology::brep::p5::{stitch_and_orient_with_identity, StitchAndOrientOptions};
use solidcore::topology::brep::p7::copy_faces;
use solidcore::topology::brep::p8::{
    create_planar_face_with_3d_curves, create_planar_faces_with_3d_curves, FaceBuildEnv,
};
use solidcore::topology::loopuse::notify_ref;
use solidcore::topology::topology_traverser::Marks;
use solidcore::topology::{Brep, FaceId};
use solidcore::{EntityRef, Kernel};

use crate::session::{finish, num, vec, Fail, Smlib};
use crate::status::{Result, Status};

/// `SMLIB_CURVE_LINE`: 6 reals, start then end.
pub const CURVE_LINE: i64 = 1;
/// `SMLIB_CURVE_ARC`: 12 reals, centre, unit axis, unit reference direction (square to
/// the axis), radius, start and end angles (radians, right-handed about the axis from the
/// reference).
pub const CURVE_ARC: i64 = 2;
/// `SMLIB_CURVE_NURBS`: kinds `3, degree, count`, then `4 count` reals (x, y, z, weight)
/// and `count + degree + 1` clamped knots.
pub const CURVE_NURBS: i64 = 3;
/// `SMLIB_CURVE_MAX_CURVES`.
pub const CURVE_MAX_CURVES: usize = 16384;
/// `SMLIB_CURVE_MAX_DEGREE`.
pub const CURVE_MAX_DEGREE: i64 = 25;
/// `SMLIB_CURVE_MAX_POINTS`.
pub const CURVE_MAX_POINTS: i64 = 100_000;
/// `SMLIB_FLAT_MAX_LEN`: the longest real array a call takes.
pub const FLAT_MAX_LEN: usize = 1 << 27;
/// `SMLIB_DRAFT_CORNER_CURVATURE`: each wall continues as itself to its neighbour.
pub const DRAFT_CORNER_CURVATURE: i32 = 0;
/// `SMLIB_DRAFT_CORNER_TANGENT`: lines' and arcs' walls continue as their tangent planes.
pub const DRAFT_CORNER_TANGENT: i32 = 1;
/// `SMLIB_DRAFT_CORNER_ROUND`: an arc about the corner.
pub const DRAFT_CORNER_ROUND: i32 = 2;
/// `SMLIB_SWEEP_FOLLOW`: the profile turns with the path.
pub const SWEEP_FOLLOW: i32 = 0;
/// `SMLIB_SWEEP_TRANSLATE`: the profile only moves along the path.
pub const SWEEP_TRANSLATE: i32 = 1;
/// `SMLIB_ENDS_OPEN`.
pub const ENDS_OPEN: i32 = 0;
/// `SMLIB_ENDS_CAPPED`.
pub const ENDS_CAPPED: i32 = 1;
/// `SMLIB_ENDS_JOINED`: the last section lies on the first; the walls' ends are glued to
/// their starts.
pub const ENDS_JOINED: i32 = 2;
/// `SMLIB_SWEEP_MAX_STATIONS`.
pub const SWEEP_MAX_STATIONS: usize = 4096;
/// `SMLIB_SWEEP_MAX_STATION_CURVES`.
pub const SWEEP_MAX_STATION_CURVES: f64 = 400_000.0;

/// `MAX_COORD`: a coordinate no larger than this (mm, 1000 km).
const MAX_COORD: f64 = 1e9;

type R<T> = core::result::Result<T, Fail>;

const fn arg() -> Fail {
    Fail::Status(Status::ERR_ARG)
}

const fn sm_err() -> Fail {
    Fail::Status(Status::SM_ERR)
}

fn coord_ok(p: &[f64]) -> bool {
    p.iter().all(|x| x.is_finite() && x.abs() <= MAX_COORD)
}

/// `a > b`, false when either is NaN: the C API tests a range as `!(a > b)`, so that a NaN
/// fails it, and reads here as `!gt(a, b)`.
fn gt(a: f64, b: f64) -> bool {
    a > b
}

/// `lo < x < hi`, false for a NaN.
fn open(x: f64, lo: f64, hi: f64) -> bool {
    x > lo && x < hi
}

/// `lo < x <= hi`, false for a NaN.
fn half_open(x: f64, lo: f64, hi: f64) -> bool {
    x > lo && x <= hi
}

fn finite(p: &[f64]) -> bool {
    p.iter().all(|x| x.is_finite())
}

fn length(v: [f64; 3]) -> f64 {
    let [x, y, z] = v;
    (x * x + y * y + z * z).sqrt()
}

/// `unit`: a unit vector from `v`, or `None` when it has no usable length.
fn unit(v: [f64; 3]) -> Option<Vec3> {
    if !finite(&v) {
        return None;
    }
    let n = length(v);
    if !gt(n, 1e-12) || !n.is_finite() {
        return None;
    }
    let [x, y, z] = v;
    Some(Vec3::new(x / n, y / n, z / n))
}

fn point(v: [f64; 3]) -> Point3 {
    let [x, y, z] = v;
    Point3::new(x, y, z)
}

/// `N` numbers of `r` from `at`, as a fixed array.
fn fixed<const N: usize>(r: &[f64], at: usize) -> R<[f64; N]> {
    let end = at.checked_add(N).ok_or_else(arg)?;
    r.get(at..end)
        .and_then(|s| <[f64; N]>::try_from(s).ok())
        .ok_or_else(arg)
}

/// `reals_of`: the number of reals curve kind `kind` takes (NURBS: with `degree` and
/// `count`), or `None`.
fn reals_of(kind: i64, degree: i64, count: i64) -> Option<usize> {
    match kind {
        CURVE_LINE => Some(6),
        CURVE_ARC => Some(12),
        CURVE_NURBS => usize::try_from(4 * count + count + degree + 1).ok(),
        _ => None,
    }
}

/// `line`.
fn line(r: [f64; 6]) -> R<Curve> {
    if !coord_ok(&r) {
        return Err(arg());
    }
    let [ax, ay, az, bx, by, bz] = r;
    let a = Point3::new(ax, ay, az);
    let b = Point3::new(bx, by, bz);
    if !gt(a.distance(b), 1e-9) {
        return Err(arg());
    }
    Ok(Curve::from(Line::create_line_segment(3, a, b, None)?))
}

/// `arc`: a circular arc, its circle's parameter starting at the arc's start (SMLib takes
/// intervals in degrees within `[-360, 360]`).
fn arc(k: &Kernel, r: [f64; 12]) -> R<Curve> {
    let [cx, cy, cz, nx, ny, nz, rx, ry, rz, radius, start, end] = r;
    if !coord_ok(&[cx, cy, cz]) || !finite(&[nx, ny, nz, rx, ry, rz, radius, start, end]) {
        return Err(arg());
    }
    let (Some(axis), Some(reference)) = (unit([nx, ny, nz]), unit([rx, ry, rz])) else {
        return Err(arg());
    };
    if axis.dot(reference).abs() > 1e-9 {
        return Err(arg());
    }
    let sweep = end - start;
    if !open(radius, 1e-9, MAX_COORD)
        || start.abs() > 1e6
        || !half_open(sweep, 1e-9, TAU * (1.0 + 1e-12))
    {
        return Err(arg());
    }
    let y = axis.cross(reference);
    let x0 = reference * start.cos() + y * start.sin();
    let y0 = axis.cross(x0);
    let frame = Axis2Placement::new(Point3::new(cx, cy, cz), x0, y0);
    let degrees = Extent1d::from_min_max(0.0, f64::min(360.0, sweep * 180.0 / PI));
    let c = Ellipse::create_canonical_circle(k, &frame, radius, Some(&degrees), None)?;
    Ok(Curve::from(c))
}

/// `nurbs`: `count` control points (x, y, z, weight), then `count + degree + 1` knots,
/// clamped (end multiplicity `degree + 1`, inner ones at most `degree`), non-decreasing.
fn nurbs(degree: i64, count: i64, r: &[f64]) -> R<Curve> {
    if !(1..=CURVE_MAX_DEGREE).contains(&degree)
        || !(degree + 1..=CURVE_MAX_POINTS).contains(&count)
    {
        return Err(arg());
    }
    let n = usize::try_from(count).map_err(|_| arg())?;
    let d = usize::try_from(degree).map_err(|_| arg())?;
    let mut points = Vec::new();
    let mut weights = Vec::new();
    let mut rational = false;
    for i in 0..n {
        let [x, y, z, w] = fixed::<4>(r, 4 * i)?;
        if !coord_ok(&[x, y, z]) || !w.is_finite() || !open(w, 1e-12, 1e12) {
            return Err(arg());
        }
        points.push(Point3::new(x, y, z));
        weights.push(w);
        rational = rational || w != 1.0;
    }
    let knots = r.get(4 * n..4 * n + n + d + 1).ok_or_else(arg)?;
    if !finite(knots) {
        return Err(arg());
    }
    let mut unique: Vec<f64> = Vec::new();
    let mut mult: Vec<usize> = Vec::new();
    let mut prev: Option<f64> = None;
    for &knot in knots {
        match prev {
            Some(p) if knot < p => return Err(arg()),
            Some(p) if knot == p => {
                if let Some(m) = mult.last_mut() {
                    *m += 1;
                }
            }
            _ => {
                unique.push(knot);
                mult.push(1);
            }
        }
        prev = Some(knot);
    }
    let (Some(&first), Some(&last), Some(&m0), Some(&m1)) =
        (unique.first(), unique.last(), mult.first(), mult.last())
    else {
        return Err(arg());
    };
    if unique.len() < 2 || m0 != d + 1 || m1 != d + 1 || !gt(last - first, 1e-12) {
        return Err(arg());
    }
    let inner = mult.get(1..mult.len() - 1).unwrap_or(&[]);
    if inner.iter().any(|&m| m > d) {
        return Err(arg());
    }
    let c = BSplineCurve::create_canonical(&CanonicalData {
        dim: 3,
        degree: d,
        control_points: &points,
        form: BSplineCurveForm::Unspecified,
        knot_multiplicities: &mult,
        knots: &unique,
        knot_type: KnotType::Unspecified,
        weights: rational.then_some(weights.as_slice()),
        trim: None,
    })?;
    Ok(Curve::from(c))
}

/// `parse`: the curves of a curve list (`kinds`, `reals`: lines, circular arcs and NURBS,
/// as `SMLIB_CURVE_*` describe them), built in order. Both slices must be used up
/// exactly; any malformed entry (unknown kind, short or extra reals, a non-finite or
/// out-of-range number, a degenerate line or arc, unclamped or decreasing knots) is
/// `SMLIB_ERR_ARG`. The crate's one curve-list decoder (F3's sweeps, F5's projections).
pub(crate) fn parse(k: &Kernel, kinds: &[i64], r: &[f64]) -> R<Vec<Curve>> {
    if kinds.is_empty() || kinds.len() > 3 * CURVE_MAX_CURVES || r.len() > FLAT_MAX_LEN {
        return Err(arg());
    }
    let mut out = Vec::new();
    let mut i = 0usize;
    let mut at = 0usize;
    while let Some(&kind) = kinds.get(i) {
        i += 1;
        if out.len() >= CURVE_MAX_CURVES {
            return Err(arg());
        }
        let (mut degree, mut count) = (0, 0);
        if kind == CURVE_NURBS {
            let (Some(&d), Some(&n)) = (kinds.get(i), kinds.get(i + 1)) else {
                return Err(arg());
            };
            i += 2;
            degree = d;
            count = n;
            if !(1..=CURVE_MAX_DEGREE).contains(&degree) || !(2..=CURVE_MAX_POINTS).contains(&count)
            {
                return Err(arg());
            }
        }
        let need = reals_of(kind, degree, count).ok_or_else(arg)?;
        let end = at
            .checked_add(need)
            .filter(|&e| e <= r.len())
            .ok_or_else(arg)?;
        let c = match kind {
            CURVE_LINE => line(fixed::<6>(r, at)?)?,
            CURVE_ARC => arc(k, fixed::<12>(r, at)?)?,
            _ => nurbs(degree, count, r.get(at..end).ok_or_else(arg)?)?,
        };
        out.push(c);
        at = end;
    }
    if at == r.len() {
        Ok(out)
    } else {
        Err(arg())
    }
}

/// `list`: a curve list for the log line (`curves[kinds=n, reals=m]`).
pub(crate) fn list(kinds: &[i64], r: &[f64]) -> String {
    format!("curves[kinds={}, reals={}]", kinds.len(), r.len())
}

/// `GetEnds`.
fn ends(k: &Kernel, c: &Curve) -> R<(Point3, Point3)> {
    let e = c.ends(k)?;
    Ok((e.start, e.end))
}

/// `one_loop`: whether the curves, in order, run end to end round one closed loop (within
/// `tol`).
fn one_loop(k: &Kernel, c: &[Curve], tol: f64) -> R<bool> {
    let mut first = Point3::ZERO;
    let mut prev = Point3::ZERO;
    for (i, x) in c.iter().enumerate() {
        let (a, b) = ends(k, x)?;
        if i == 0 {
            first = a;
        } else if a.distance(prev) > tol {
            return Ok(false);
        }
        prev = b;
    }
    Ok(!c.is_empty() && prev.distance(first) <= tol)
}

/// `straight`: whether the path is lines end to end along one direction: then that unit
/// direction and the path's length.
fn straight(k: &Kernel, path: &[Curve]) -> R<Option<(Vec3, f64)>> {
    let mut distance = 0.0;
    let mut dir = Vec3::ZERO;
    let mut prev = Point3::ZERO;
    for (i, c) in path.iter().enumerate() {
        if !matches!(c.kind(), CurveKind::Line(_)) {
            return Ok(None);
        }
        let (a, b) = ends(k, c)?;
        if i > 0 && a.distance(prev) > 1e-9 {
            return Ok(None);
        }
        prev = b;
        let mut d = b - a;
        let l = d.length();
        if !gt(l, 1e-12) {
            return Ok(None);
        }
        d = d / l;
        if i == 0 {
            dir = d;
        } else if d.dot(dir) < 1.0 - 1e-12 {
            return Ok(None);
        }
        distance += l;
    }
    Ok((!path.is_empty()).then_some((dir, distance)))
}

/// An arc's numbers from a path's curve list (`reals`, 12 from `at`).
fn arc_reals(r: &[f64], at: usize) -> R<[f64; 12]> {
    fixed::<12>(r, at)
}

/// `circular`: whether the path is arcs end to end about one axis: then the centre, the
/// unit axis they turn about counter-clockwise and the total angle (radians).
fn circular(
    k: &Kernel,
    path: &[Curve],
    kinds: &[i64],
    r: &[f64],
) -> R<Option<(Point3, Vec3, f64)>> {
    if kinds.len() != path.len() || path.is_empty() {
        return Ok(None);
    }
    let mut angle = 0.0;
    let mut centre = Point3::ZERO;
    let mut axis = Vec3::ZERO;
    let mut radius = 0.0;
    let mut prev = Point3::ZERO;
    for (i, (&kind, c)) in kinds.iter().zip(path).enumerate() {
        if kind != CURVE_ARC {
            return Ok(None);
        }
        let [cx, cy, cz, nx, ny, nz, _, _, _, rad, a0, a1] = arc_reals(r, 12 * i)?;
        let centre_i = Point3::new(cx, cy, cz);
        let mut axis_i = Vec3::new(nx, ny, nz);
        axis_i = axis_i / axis_i.length();
        let (a, b) = ends(k, c)?;
        if i == 0 {
            centre = centre_i;
            axis = axis_i;
            radius = rad;
        } else if centre_i.distance(centre) > 1e-9
            || axis_i.dot(axis) < 1.0 - 1e-12
            || (rad - radius).abs() > 1e-9
            || a.distance(prev) > 1e-9
        {
            return Ok(None);
        }
        prev = b;
        angle += a1 - a0;
    }
    Ok((angle <= TAU * (1.0 + 1e-12)).then_some((centre, axis, angle)))
}

/// `SM_CAST_PTR(SmBSplineCurve, c)`: the B-spline part of the kinds derived from
/// `SmBSplineCurve`.
fn bspline_part(c: &Curve) -> Option<&BSplineCurve> {
    match c.kind() {
        CurveKind::BSpline(b) => Some(b),
        CurveKind::Line(l) => Some(l.base()),
        CurveKind::Ellipse(e) => Some(e.base()),
        CurveKind::Hyperbola(h) => Some(h.base()),
        CurveKind::Parabola(p) => Some(p.base()),
        _ => None,
    }
}

/// `joined`: a path of end-to-end curves as one B-spline (`CreateCurveSweep` takes a single
/// curve), from plain B-spline copies (lines and circles keep an analytic
/// parameterisation).
fn joined(path: &[Curve]) -> R<Curve> {
    let mut parts = Vec::new();
    for c in path {
        parts.push(bspline_part(c).ok_or_else(arg)?.copy());
    }
    if parts.len() == 1 {
        return parts.pop().map(Curve::from).ok_or_else(arg);
    }
    Ok(Curve::from(BSplineCurve::create_by_joining(
        &mut parts, None,
    )?))
}

/// `put_motion`: writes motion `k` (0 first, 1 middle, 2 last) of `smlib_sweep_placed`: the
/// columns, then the translation.
fn put_motion(placed: &mut [f64; 36], k: usize, x: Vec3, y: Vec3, z: Vec3, t: Point3) {
    let m = [x.x, x.y, x.z, y.x, y.y, y.z, z.x, z.y, z.z, t.x, t.y, t.z];
    if let Some(slot) = placed.get_mut(12 * k..12 * k + 12) {
        slot.copy_from_slice(&m);
    }
}

/// `turned`: `v` turned `angle` right-handed about the unit `k` (Rodrigues).
fn turned(v: Vec3, k: Vec3, angle: f64) -> Vec3 {
    v * angle.cos() + k.cross(v) * angle.sin() + k * (k.dot(v) * (1.0 - angle.cos()))
}

/// `Rigid`: a rigid motion `p -> R p + t`, `R` by its columns.
#[derive(Clone, Copy)]
struct Rigid {
    x: Vec3,
    y: Vec3,
    z: Vec3,
    t: Vec3,
}

impl Default for Rigid {
    fn default() -> Self {
        Self {
            x: Vec3::new(1.0, 0.0, 0.0),
            y: Vec3::new(0.0, 1.0, 0.0),
            z: Vec3::new(0.0, 0.0, 1.0),
            t: Vec3::ZERO,
        }
    }
}

impl Rigid {
    /// This motion, then a turn of `angle` right-handed about the unit `k` through `c`.
    fn turned_by(self, c: Vec3, k: Vec3, angle: f64) -> Self {
        Self {
            x: turned(self.x, k, angle),
            y: turned(self.y, k, angle),
            z: turned(self.z, k, angle),
            t: c + turned(self.t - c, k, angle),
        }
    }

    /// This motion, then a move by `d`.
    fn moved_by(self, d: Vec3) -> Self {
        Self {
            t: self.t + d,
            ..self
        }
    }

    /// As the placement `SmCurve::Transform` takes: the world frame's image.
    fn placement(self) -> Axis2Placement {
        Axis2Placement::new(self.t, self.x, self.y)
    }
}

/// `Segment`: one segment of a path of lines and arcs: a line from `a` along the unit
/// `dir` for `length`, or an arc turning `angle` (radians) right-handed about the unit `k`
/// through `c`, `length` along.
#[derive(Clone, Copy, Default)]
struct Segment {
    arc: bool,
    dir: Vec3,
    c: Vec3,
    k: Vec3,
    angle: f64,
    length: f64,
}

impl Segment {
    /// The motion `m` carried `s` of the way along (`0 <= s <= length`).
    fn carry(&self, m: Rigid, s: f64) -> Rigid {
        if self.arc {
            m.turned_by(self.c, self.k, self.angle * s / self.length)
        } else {
            m.moved_by(self.dir * s)
        }
    }
}

/// `segmented`: whether the path is two or more lines and arcs end to end, each leaving in
/// the direction the one before arrives in, and not closed: then its segments.
fn segmented(k: &Kernel, path: &[Curve], kinds: &[i64], r: &[f64]) -> R<Option<Vec<Segment>>> {
    if kinds.len() != path.len() || kinds.len() < 2 {
        return Ok(None);
    }
    let mut out = Vec::new();
    let mut at = 0usize;
    let mut before = Vec3::ZERO;
    let mut first = Point3::ZERO;
    let mut prev = Point3::ZERO;
    for (i, (&kind, c)) in kinds.iter().zip(path).enumerate() {
        let mut g = Segment::default();
        let (a, b) = ends(k, c)?;
        let (leaving, arriving);
        if kind == CURVE_LINE {
            g.length = (b - a).length();
            if !gt(g.length, 1e-12) {
                return Ok(None);
            }
            g.dir = (b - a) / g.length;
            leaving = g.dir;
            arriving = g.dir;
            at += 6;
        } else if kind == CURVE_ARC {
            let [cx, cy, cz, nx, ny, nz, _, _, _, rad, a0, a1] = arc_reals(r, at)?;
            g.arc = true;
            g.c = Vec3::new(cx, cy, cz);
            g.k = Vec3::new(nx, ny, nz);
            g.k = g.k / g.k.length();
            g.angle = a1 - a0;
            g.length = rad * g.angle;
            let l = g.k.cross(a - g.c);
            let e = g.k.cross(b - g.c);
            leaving = l / l.length();
            arriving = e / e.length();
            at += 12;
        } else {
            return Ok(None);
        }
        if i == 0 {
            first = a;
        } else if a.distance(prev) > 1e-9 || leaving.dot(before) < 1.0 - 1e-12 {
            return Ok(None);
        }
        prev = b;
        before = arriving;
        out.push(g);
    }
    if first.distance(prev) <= 1e-6 {
        return Ok(None);
    }
    Ok(Some(out))
}

/// `carried`: the profile's curves carried by `m`, as copies.
fn carried(k: &Kernel, profile: &[Curve], m: Rigid) -> R<Vec<Curve>> {
    let placement = m.placement();
    let mut out = Vec::new();
    for src in profile {
        let mut copy = src.copy(k).map_err(|_| sm_err())?;
        copy.transform(k, &placement, None).map_err(|_| sm_err())?;
        out.push(copy);
    }
    Ok(out)
}

/// Runs `body` with a face-construction environment over an arena-only notification store
/// (faces named by `loopuse::notify_ref`), as `SmApiCreatePlanarFaces` builds its faces.
fn with_face_env<T>(
    k: &mut Kernel,
    b: &mut Brep,
    body: impl FnOnce(&mut Kernel, &mut Brep, &mut FaceBuildEnv<'_, '_>) -> solidcore::Result<T>,
) -> solidcore::Result<T> {
    let mut state = ArenaNotifyState::new(b);
    let mut env = state.env(None);
    let mut marks = Marks::default();
    let this_face = |f: FaceId| notify_ref(EntityRef::Face(f));
    let mut fenv = FaceBuildEnv {
        env: &mut env,
        marks: &mut marks,
        this_face: &this_face,
    };
    body(k, b, &mut fenv)
}

/// `StitchAndOrient` of a B-rep whose faces were made under `loopuse::notify_ref`, its
/// spatial cache keyed by those names.
fn stitch_and_orient(k: &mut Kernel, b: &mut Brep) -> solidcore::Result<()> {
    let key_of: &dyn Fn(EntityRef) -> solidcore::geom::object::ObjectRef = &notify_ref;
    let keys = |object| caller_cache_key(key_of, object);
    stitch_and_orient_with_identity(
        k,
        b,
        &StitchAndOrientOptions::default(),
        Some(BrepIdentity {
            keys: &keys,
            ..BrepIdentity::default()
        }),
    )
}

/// `piecewise`: the profile swept along a path of lines and arcs, turning with it: each
/// segment swept on its own, exactly (a line by the linear sweep, an arc by the rotational
/// sweep), the profile carried from each to the next, and the walls sewn along the
/// sections they share.
fn piecewise(
    k: &mut Kernel,
    profile: &[Curve],
    segs: &[Segment],
    cap: bool,
    placed: Option<&mut [f64; 36]>,
) -> R<Brep> {
    let mut brep = new_brep(k, USE_DEFAULT, USE_DEFAULT)?;
    let mut m = Rigid::default();
    let mut starts = Vec::new();
    let mut total = 0.0;
    for g in segs {
        starts.push(m);
        total += g.length;
        let section = carried(k, profile, m)?;
        let refs: Vec<Option<&Curve>> = section.iter().map(Some).collect();
        let mut piece = if g.arc {
            create_rotational_sweep(k, &refs, g.c, g.k, g.angle * 180.0 / PI, false)?
        } else {
            create_linear_sweep(k, &refs, g.dir, g.length, false)?
        };
        let faces: Vec<FaceId> = piece.face_list().to_vec();
        let mut state = ArenaNotifyState::new(&brep);
        let mut env = state.env(None);
        let this_face = |f: FaceId| notify_ref(EntityRef::Face(f));
        let (_, all) = copy_faces(
            k, &mut piece, &faces, &mut brep, false, &this_face, &mut env,
        )?;
        if !all {
            return Err(sm_err());
        }
        m = g.carry(m, g.length);
    }
    let first = starts.first().copied().unwrap_or_default();
    if cap {
        // The caps: the profile where the path starts and where it ends.
        for at in [first, m] {
            let section = carried(k, profile, at)?;
            with_face_env(k, &mut brep, |k, b, fenv| {
                let tol = b.tolerance();
                let region = b.infinite_region();
                create_planar_faces_with_3d_curves(k, b, region, &section, tol, fenv)
            })?;
        }
    }
    stitch_and_orient(k, &mut brep)?;
    if let Some(placed) = placed {
        // The motions to the first section, the one half-way along and the last.
        let mut half = 0.5 * total;
        let mut middle = m;
        for (i, (g, start)) in segs.iter().zip(&starts).enumerate() {
            if half <= g.length || i + 1 == segs.len() {
                middle = g.carry(*start, half.min(g.length));
                break;
            }
            half -= g.length;
        }
        for (j, motion) in [first, middle, m].into_iter().enumerate() {
            put_motion(placed, j, motion.x, motion.y, motion.z, motion.t);
        }
    }
    finish(brep, cap, false)
}

/// `maps_ok`: each station map finite, its matrix of bounded size and far from singular, its
/// translation within `MAX_COORD`.
fn maps_ok(maps: &[[f64; 12]]) -> bool {
    maps.iter().all(|m| {
        let [m0, m1, m2, m3, m4, m5, m6, m7, m8, t0, t1, t2] = *m;
        if !finite(m) || !coord_ok(&[t0, t1, t2]) {
            return false;
        }
        let c0 = length([m0, m1, m2]);
        let c1 = length([m3, m4, m5]);
        let c2 = length([m6, m7, m8]);
        if ![c0, c1, c2].iter().all(|&c| c < 1e6) {
            return false;
        }
        let det = m0 * (m4 * m8 - m7 * m5) - m3 * (m1 * m8 - m7 * m2) + m6 * (m1 * m5 - m4 * m2);
        det.abs() > 1e-9 * c0 * c1 * c2
    })
}

/// `placed`: a copy of `c` as a plain B-spline with every control point carried by the
/// affine map `m` (`x -> M x + t`), the weights unchanged.
fn placed(c: &BSplineCurve, m: &[f64; 12]) -> R<Curve> {
    let canonical = c.canonical().map_err(|_| sm_err())?;
    if canonical.dim != 3 {
        return Err(sm_err());
    }
    let [m0, m1, m2, m3, m4, m5, m6, m7, m8, t0, t1, t2] = *m;
    let points: Vec<Point3> = canonical
        .control_points
        .iter()
        .map(|p| {
            Point3::new(
                m0 * p.x + m3 * p.y + m6 * p.z + t0,
                m1 * p.x + m4 * p.y + m7 * p.z + t1,
                m2 * p.x + m5 * p.y + m8 * p.z + t2,
            )
        })
        .collect();
    let r = BSplineCurve::create_canonical(&CanonicalData {
        dim: 3,
        degree: canonical.degree,
        control_points: &points,
        form: BSplineCurveForm::Unspecified,
        knot_multiplicities: &canonical.knot_multiplicities,
        knots: &canonical.knots,
        knot_type: KnotType::Unspecified,
        weights: (!canonical.weights.is_empty()).then_some(canonical.weights.as_slice()),
        trim: None,
    })?;
    Ok(Curve::from(r))
}

impl Smlib {
    /// `smlib_extrude`: the planar profile extruded `distance` (> 0) along `direction`,
    /// square to the profile's plane. `draft` (radians, `|draft| < pi/2`) leans the walls:
    /// 0 straight, positive inward, negative outward; a drafted profile is one closed loop,
    /// its curves in order, and `direction` must be its normal (either way round). `corner`
    /// (one of the `DRAFT_CORNER_*`, read only with a draft) closes its sharp corners; a
    /// corner it cannot close is `SMLIB_ERR_CORNER_APART`, `SMLIB_ERR_WALL_CLOSES` or
    /// `SMLIB_ERR_NOT_EXACT`, and [`Smlib::last_failed_corner_point`] gives that corner's
    /// point on the profile (`out_corner`).
    /// With `cap` 1 the ends are capped and the result must be a closed manifold solid.
    pub fn extrude(
        &mut self,
        kinds: &[i64],
        reals: &[f64],
        direction: [f64; 3],
        distance: f64,
        draft: f64,
        corner: i32,
        cap: i32,
    ) -> Result<Brep> {
        if !distance.is_finite()
            || !open(distance, 1e-9, MAX_COORD)
            || !draft.is_finite()
            || draft.abs() >= PI / 2.0 - 1e-6
            || (cap != 0 && cap != 1)
            || !matches!(
                corner,
                DRAFT_CORNER_CURVATURE | DRAFT_CORNER_TANGENT | DRAFT_CORNER_ROUND
            )
        {
            return Err(Status::ERR_ARG);
        }
        self.run(
            || {
                let args = format!(
                    "{}, direction={}, distance={}, draft={}",
                    list(kinds, reals),
                    vec(&direction),
                    num(distance),
                    num(draft)
                );
                format!("smlib_extrude({args}, corner={corner}, cap={cap})")
            },
            &[],
            |api| {
                let dir = unit(direction).ok_or_else(arg)?;
                let profile = parse(api.kernel(), kinds, reals)?;
                let refs: Vec<Option<&Curve>> = profile.iter().map(Some).collect();
                let r = if draft == 0.0 {
                    create_linear_sweep(api.kernel_mut(), &refs, dir, distance, cap != 0)?
                } else {
                    // SMLib takes the draft as the walls' angle to the profile's normal, in
                    // degrees, positive leaning in; the height's sign picks the side, so it
                    // is matched to `direction`. The draft extrusion takes a single chain.
                    if !one_loop(api.kernel(), &profile, 1e-6)? {
                        return Err(arg());
                    }
                    let normal = compute_crvs_normal(api.kernel(), &profile).map_err(|_| arg())?;
                    let along = normal.dot(dir);
                    if along.abs() < 1.0 - 1e-9 {
                        // A draft extrusion runs square to the profile.
                        return Err(arg());
                    }
                    let height = if along > 0.0 { distance } else { -distance };
                    let angle = draft * 180.0 / PI;
                    let rule = match corner {
                        DRAFT_CORNER_TANGENT => OffsetCornerType::LinearExtension,
                        DRAFT_CORNER_ROUND => OffsetCornerType::FilletCorner,
                        _ => OffsetCornerType::CurveExtension,
                    };
                    let ends = if cap != 0 { 3 } else { 0 };
                    match create_taper_extrude_with_corner(
                        api.kernel_mut(),
                        &refs,
                        height,
                        angle,
                        ends,
                        rule,
                    )? {
                        TaperExtrude::Solid(b) => b,
                        TaperExtrude::Refused(refusal) => {
                            let p = refusal.corner;
                            api.set_failed_corner_point([p.x, p.y, p.z]);
                            return Err(Fail::Status(match refusal.reason {
                                TaperCornerRefusal::Apart => Status::ERR_CORNER_APART,
                                TaperCornerRefusal::WallCloses => Status::ERR_WALL_CLOSES,
                                TaperCornerRefusal::NotExact => Status::ERR_NOT_EXACT,
                            }));
                        }
                    }
                };
                finish(r, cap != 0, false)
            },
        )
    }

    /// `smlib_revolve`: the planar profile turned `angle` radians (`0 < angle <= 2 pi`)
    /// right-handed about the axis through `origin` along `axis`, which lies in the
    /// profile's plane; capped with `cap` 1.
    pub fn revolve(
        &mut self,
        kinds: &[i64],
        reals: &[f64],
        origin: [f64; 3],
        axis: [f64; 3],
        angle: f64,
        cap: i32,
    ) -> Result<Brep> {
        if !coord_ok(&origin)
            || !angle.is_finite()
            || !half_open(angle, 1e-9, TAU * (1.0 + 1e-12))
            || (cap != 0 && cap != 1)
        {
            return Err(Status::ERR_ARG);
        }
        self.run(
            || {
                format!(
                    "smlib_revolve({}, origin={}, axis={}, angle={}, cap={cap})",
                    list(kinds, reals),
                    vec(&origin),
                    vec(&axis),
                    num(angle)
                )
            },
            &[],
            |api| {
                let dir = unit(axis).ok_or_else(arg)?;
                let profile = parse(api.kernel(), kinds, reals)?;
                let refs: Vec<Option<&Curve>> = profile.iter().map(Some).collect();
                let degrees = f64::min(360.0, angle * 180.0 / PI);
                let r = create_rotational_sweep(
                    api.kernel_mut(),
                    &refs,
                    point(origin),
                    dir,
                    degrees,
                    cap != 0,
                )?;
                finish(r, cap != 0, false)
            },
        )
    }

    /// `smlib_sweep`: [`Smlib::sweep_placed`] without the placements.
    pub fn sweep(
        &mut self,
        kinds: &[i64],
        reals: &[f64],
        path_kinds: &[i64],
        path_reals: &[f64],
        frame: i32,
        cap: i32,
    ) -> Result<Brep> {
        self.sweep_with(kinds, reals, path_kinds, path_reals, frame, cap, false)
            .map(|(b, _)| b)
    }

    /// `smlib_sweep_placed`: the profile swept along a path of end-to-end curves, carried
    /// with the path's start (`frame`: [`SWEEP_FOLLOW`] or [`SWEEP_TRANSLATE`]). Along a
    /// straight path the sweep is an extrusion, turning with a path of arcs about one axis
    /// a revolution, and turning with a path of lines and arcs one exact sweep per piece;
    /// otherwise the walls are B-spline approximations of the swept surface. Round a closed
    /// path the walls close on each other (no caps). With the body come three rigid
    /// motions, 12 numbers each (a 3 x 3 matrix by columns, then a translation), carrying
    /// the profile onto the sweep's first, middle and last sections (`placed`).
    pub fn sweep_placed(
        &mut self,
        kinds: &[i64],
        reals: &[f64],
        path_kinds: &[i64],
        path_reals: &[f64],
        frame: i32,
        cap: i32,
    ) -> Result<(Brep, [f64; 36])> {
        self.sweep_with(kinds, reals, path_kinds, path_reals, frame, cap, true)
            .map(|(b, placed)| (b, placed.unwrap_or([0.0; 36])))
    }

    /// `smlib_sweep_placed` with `placed` given (`with_placements`) or `NULL`.
    fn sweep_with(
        &mut self,
        kinds: &[i64],
        reals: &[f64],
        path_kinds: &[i64],
        path_reals: &[f64],
        frame: i32,
        cap: i32,
        with_placements: bool,
    ) -> Result<(Brep, Option<[f64; 36]>)> {
        if (frame != SWEEP_FOLLOW && frame != SWEEP_TRANSLATE) || (cap != 0 && cap != 1) {
            return Err(Status::ERR_ARG);
        }
        self.run(
            || {
                format!(
                    "smlib_sweep({}, path {}, frame={frame}, cap={cap})",
                    list(kinds, reals),
                    list(path_kinds, path_reals)
                )
            },
            &[],
            |api| {
                let mut placed = with_placements.then_some([0.0; 36]);
                let profile = parse(api.kernel(), kinds, reals)?;
                let path = parse(api.kernel(), path_kinds, path_reals)?;
                let refs: Vec<Option<&Curve>> = profile.iter().map(Some).collect();
                // Along a straight path the sweep is an extrusion: exact walls.
                if let Some((dir, distance)) = straight(api.kernel(), &path)? {
                    let r = create_linear_sweep(api.kernel_mut(), &refs, dir, distance, cap != 0)?;
                    if let Some(p) = placed.as_mut() {
                        for k in 0..3 {
                            let t = dir * (distance * 0.5 * k as f64);
                            put_motion(
                                p,
                                k,
                                Vec3::new(1.0, 0.0, 0.0),
                                Vec3::new(0.0, 1.0, 0.0),
                                Vec3::new(0.0, 0.0, 1.0),
                                t,
                            );
                        }
                    }
                    return finish(r, cap != 0, false).map(|b| (b, placed));
                }
                // Turning with a path of arcs about one axis, the profile turns about it.
                if frame == SWEEP_FOLLOW {
                    if let Some((centre, axis, angle)) =
                        circular(api.kernel(), &path, path_kinds, path_reals)?
                    {
                        let degrees = f64::min(360.0, angle * 180.0 / PI);
                        let r = create_rotational_sweep(
                            api.kernel_mut(),
                            &refs,
                            centre,
                            axis,
                            degrees,
                            cap != 0,
                        )?;
                        if let Some(p) = placed.as_mut() {
                            let n = axis.length();
                            if !gt(n, 0.0) {
                                return Err(arg());
                            }
                            let k = axis / n;
                            let turn = degrees * PI / 180.0;
                            for j in 0..3 {
                                let a = turn * 0.5 * j as f64;
                                let t = centre - turned(centre, k, a);
                                put_motion(
                                    p,
                                    j,
                                    turned(Vec3::new(1.0, 0.0, 0.0), k, a),
                                    turned(Vec3::new(0.0, 1.0, 0.0), k, a),
                                    turned(Vec3::new(0.0, 0.0, 1.0), k, a),
                                    t,
                                );
                            }
                        }
                        return finish(r, cap != 0, false).map(|b| (b, placed));
                    }
                    // Along lines and arcs, turning with the path, the sweep is exact piece
                    // by piece.
                    if let Some(segs) = segmented(api.kernel(), &path, path_kinds, path_reals)? {
                        let b = piecewise(
                            api.kernel_mut(),
                            &profile,
                            &segs,
                            cap != 0,
                            placed.as_mut(),
                        )?;
                        return Ok((b, placed));
                    }
                }
                let along = joined(&path)?;
                // Round a closed path the profile comes back onto itself: no caps.
                let (first, _) = ends(api.kernel(), path.first().ok_or_else(arg)?)?;
                let (_, last) = ends(api.kernel(), path.last().ok_or_else(arg)?)?;
                let closed = first.distance(last) <= 1e-6;
                let options = SweepOptions {
                    cap_ends: cap != 0 && !closed,
                    translational_sweep: frame == SWEEP_TRANSLATE,
                    profile_at_path_start: true,
                    ..SweepOptions::default()
                };
                let mut r = if let Some(p) = placed.as_mut() {
                    let s = create_curve_sweep_with_placements(
                        api.kernel_mut(),
                        &refs,
                        Some(&along),
                        None,
                        None,
                        Some(&options),
                    )?;
                    for (k, m) in s.placements.iter().enumerate() {
                        put_motion(p, k, m.x_axis(), m.y_axis(), m.z_axis(), m.origin());
                    }
                    s.sweep.brep
                } else {
                    create_curve_sweep(
                        api.kernel_mut(),
                        &refs,
                        Some(&along),
                        None,
                        None,
                        Some(&options),
                    )?
                    .brep
                };
                if closed {
                    let tol = r.tolerance();
                    let s = stitch_faces(api.kernel_mut(), &mut r, tol, None, false)?;
                    if s.lamina_edges != 0 {
                        return Err(Fail::Status(Status::ERR_NOT_SOLID));
                    }
                }
                finish(r, cap != 0, false).map(|b| (b, placed))
            },
        )
    }

    /// `smlib_loft`: a loft (skin) through the sections given one after the other in the
    /// curve list, `section_sizes[i]` curves for section `i` (at least two sections); a
    /// section with holes holds its loops one after the other, its outline first. `degree`
    /// (1 ruled, 2, 3) is the walls' degree across the sections; `cap` 1 closes the ends.
    pub fn loft(
        &mut self,
        kinds: &[i64],
        reals: &[f64],
        section_sizes: &[i64],
        degree: i32,
        cap: i32,
    ) -> Result<Brep> {
        if section_sizes.len() < 2
            || section_sizes.len() > CURVE_MAX_CURVES
            || !(1..=3).contains(&degree)
            || (cap != 0 && cap != 1)
        {
            return Err(Status::ERR_ARG);
        }
        self.run(
            || {
                format!(
                    "smlib_loft({}, sections={}, degree={degree}, cap={cap})",
                    list(kinds, reals),
                    section_sizes.len()
                )
            },
            &[],
            |api| {
                let all = parse(api.kernel(), kinds, reals)?;
                let mut counts = Vec::new();
                let mut total = 0usize;
                for &size in section_sizes {
                    let n = usize::try_from(size)
                        .ok()
                        .filter(|n| (1..=CURVE_MAX_CURVES).contains(n))
                        .ok_or_else(arg)?;
                    total += n;
                    counts.push(n);
                }
                if total != all.len() {
                    return Err(arg());
                }
                let refs: Vec<Option<&Curve>> = all.iter().map(Some).collect();
                let r = create_skin_primitive(
                    api.kernel_mut(),
                    &counts,
                    &refs,
                    degree as usize,
                    cap != 0,
                )?
                .brep;
                finish(r, cap != 0, false)
            },
        )
    }

    /// `smlib_loft_to_point`: the closed loop of lines `kinds`/`reals` (at least three, in
    /// order) lofted to `apex`, off the loop's plane: a pyramid of planar faces.
    pub fn loft_to_point(&mut self, kinds: &[i64], reals: &[f64], apex: [f64; 3]) -> Result<Brep> {
        if !coord_ok(&apex) {
            return Err(Status::ERR_ARG);
        }
        self.run(
            || {
                format!(
                    "smlib_loft_to_point({}, apex={})",
                    list(kinds, reals),
                    vec(&apex)
                )
            },
            &[],
            |api| {
                let section = parse(api.kernel(), kinds, reals)?;
                if section.len() < 3 || !one_loop(api.kernel(), &section, 1e-6)? {
                    return Err(arg());
                }
                // The base's corners, each line's start.
                let mut corners = Vec::new();
                for x in &section {
                    if !matches!(x.kind(), CurveKind::Line(_)) {
                        return Err(arg());
                    }
                    corners.push(ends(api.kernel(), x)?.0);
                }
                let normal = compute_crvs_normal(api.kernel(), &section).map_err(|_| arg())?;
                let top = point(apex);
                let first = corners.first().copied().ok_or_else(arg)?;
                let size = corners
                    .iter()
                    .fold(0.0, |s: f64, q| s.max(q.distance(first)));
                let up = top - first;
                if !gt(up.dot(normal).abs(), 1e-9 * size.max(1.0)) {
                    return Err(arg());
                }
                let k = api.kernel_mut();
                let mut brep = new_brep(k, USE_DEFAULT, USE_DEFAULT)?;
                let region = brep.infinite_region();
                let tol = 1e-6;
                // Each face is made of new lines it takes over.
                let face = |k: &mut Kernel, brep: &mut Brep, pts: &[Point3]| -> R<()> {
                    let mut edges = Vec::new();
                    for (i, &p) in pts.iter().enumerate() {
                        let q = pts.get((i + 1) % pts.len()).copied().ok_or_else(arg)?;
                        edges.push(Curve::from(Line::create_line_segment(3, p, q, None)?));
                    }
                    with_face_env(k, brep, |k, b, fenv| {
                        create_planar_face_with_3d_curves(k, b, region, &edges, tol, fenv)
                    })?;
                    Ok(())
                };
                face(k, &mut brep, &corners)?;
                for (i, &c) in corners.iter().enumerate() {
                    let next = corners
                        .get((i + 1) % corners.len())
                        .copied()
                        .ok_or_else(arg)?;
                    face(k, &mut brep, &[c, next, top])?;
                }
                stitch_and_orient(k, &mut brep)?;
                finish(brep, true, false)
            },
        )
    }

    /// `smlib_sweep_stations`: the planar profile carried through the stations of `maps`
    /// (2 to [`SWEEP_MAX_STATIONS`]), each an affine map: a 3 x 3 matrix by columns, then a
    /// translation. Each section is the profile moved by its map exactly (its control points
    /// mapped); the walls are skinned through the sections, cubic across them (lower with
    /// fewer than four). `ends`: [`ENDS_OPEN`], [`ENDS_CAPPED`] or [`ENDS_JOINED`].
    pub fn sweep_stations(
        &mut self,
        kinds: &[i64],
        reals: &[f64],
        maps: &[[f64; 12]],
        ends: i32,
    ) -> Result<Brep> {
        if maps.len() < 2
            || maps.len() > SWEEP_MAX_STATIONS
            || !matches!(ends, ENDS_OPEN | ENDS_CAPPED | ENDS_JOINED)
        {
            return Err(Status::ERR_ARG);
        }
        self.run(
            || {
                format!(
                    "smlib_sweep_stations({}, stations={}, ends={ends})",
                    list(kinds, reals),
                    maps.len()
                )
            },
            &[],
            |api| {
                if !maps_ok(maps) {
                    return Err(arg());
                }
                let profile = parse(api.kernel(), kinds, reals)?;
                if profile.len() as f64 * maps.len() as f64 > SWEEP_MAX_STATION_CURVES {
                    return Err(arg());
                }
                // Every section is the profile carried by its station's map: the sections
                // share the profile's degrees and knots, which the skin interpolates.
                let mut sections = Vec::new();
                let mut counts = Vec::new();
                for m in maps {
                    for x in &profile {
                        let b = bspline_part(x).ok_or_else(arg)?;
                        sections.push(placed(b, m)?);
                    }
                    counts.push(profile.len());
                }
                let degree = (maps.len() - 1).min(3);
                let refs: Vec<Option<&Curve>> = sections.iter().map(Some).collect();
                let mut r = create_skin_primitive(
                    api.kernel_mut(),
                    &counts,
                    &refs,
                    degree,
                    ends == ENDS_CAPPED,
                )?
                .brep;
                if ends == ENDS_JOINED {
                    // The last section lies on the first: the walls' ends are glued to their
                    // starts.
                    let s = api_stitch_into_solid(api.kernel_mut(), Some(&mut r))?;
                    if !s.produces_a_solid {
                        return Err(Fail::Status(Status::ERR_NOT_SOLID));
                    }
                }
                finish(r, ends != ENDS_OPEN, false)
            },
        )
    }

    /// `smlib_planar_faces`: a sheet of planar faces bounded by the closed loops of the
    /// coplanar curves (any order and orientation, nested to any level).
    pub fn planar_faces(&mut self, kinds: &[i64], reals: &[f64]) -> Result<Brep> {
        self.run(
            || format!("smlib_planar_faces({})", list(kinds, reals)),
            &[],
            |api| {
                let loops = parse(api.kernel(), kinds, reals)?;
                let r = create_planar_faces(api.kernel_mut(), &loops)?;
                finish(r, false, false)
            },
        )
    }
}
