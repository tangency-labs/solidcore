//! (not upstream) A closed surface's seam turned into the widest gap of its face's
//! boundary: the surface placed again about its own axis so that its angle parameter
//! starts where the face does not reach, or on the face's own seam.
//!
//! Upstream: none (solidcore's STEP reader, #2578). The idea is SolidCraft's reader
//! (`crates/step/src/topology/swept.rs` `widest_gap`, `common_angle`, `span`,
//! `crates/step/src/topology/frames.rs` `about_axis`, `quarters_turned`,
//! `turn_frames`, and `crates/geom/src/brep/form.rs` `tube_turned`, at a5b71ca2;
//! MIT OR Apache-2.0) and the frame step-brep's `seam_frame` and `face_seam_point`
//! choose for SMLib's import (left as they are for #2585). The computations are the
//! kernel's own: the surfaces' STEP angles, [`Axis2Placement`] and
//! [`Surface::transform`].
//!
//! | Rust                      | SolidCraft                                                   |
//! |---------------------------|--------------------------------------------------------------|
//! | [`seam_into_widest_gap`]  | `swept::span` (a closed conic's start), `frames::turn_frames` (the turn), `form::tube_turned` |
//! | [`PlacedSeam`]            | the turned frame and `Span::Turn { start }`                  |
//! | [`ClosedSeam`]            | `swept::widest_gap`'s `(start, width)`, `common_angle`        |
//!
//! Deviations: not a port. See [`seam_into_widest_gap`] for the rules. SolidCraft and
//! step-brep take the boundary's angles from sampled points (step-brep 17 an edge, for
//! cylinders and cones only); here they come from each edge's extrema, for every
//! closed kind. SolidCraft starts a closed profile where the face starts (the end of
//! the gap); the seam here goes in the gap's middle, as SMLib's import puts it, so
//! the face keeps half the gap on either side of it.

use crate::curve::curve::CurveKind;
use crate::curve::ellipse::Ellipse;
use crate::error::{Result, SolidError};
use crate::geom::axis2_placement::Axis2Placement;
use crate::geom::extent1d::Extent1d;
use crate::geom::vector3d::{acos, sqrt};
use crate::kernel::Kernel;
use crate::math::{rad_to_deg, Point2, Point3, Vec2, PI};
use crate::surface::step_window::BoundaryEdge;
use crate::surface::surface::{Surface, SurfaceKind};

const OP: &str = "surface::step_seam::seam_into_widest_gap";

/// One turn, in degrees (the STEP angles' period).
const TURN: f64 = 360.0;
/// The longest step between the key points of a conic edge, in its STEP angle: a
/// quarter turn, so that a step of a circle round the tube turns a quarter of the tube
/// at most and a step of a monotone arc less than a full turn.
const STEP: f64 = 90.0;
/// The longest conic edge taken, in its STEP angle: two turns (the kernel keeps a
/// conic's STEP interval within `[-360, 360]`).
const MAX_SPAN: f64 = 2.0 * TURN;

/// Where the seam of one closed direction of a surface goes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClosedSeam {
    /// The seam's angle in degrees, in `[0, 360)`, in the input surface's STEP
    /// parameter of that direction; `None` when the seam stays where it is (the
    /// boundary covers the full turn and no seam point is given).
    pub angle: Option<f64>,
    /// The widest gap of the boundary's angles in that direction, in degrees: its
    /// minimum in `[0, 360)` is where the gap starts, its length is the gap's width
    /// (its maximum may pass 360). `None` when the boundary covers the full turn.
    pub gap: Option<Extent1d>,
}

/// A closed surface with its seam in the widest gap of a face's boundary.
#[derive(Debug, Clone)]
pub struct PlacedSeam {
    /// The input surface turned about its own axis so that its STEP u (the angle
    /// about the axis) is 0 at [`Self::u`]'s angle: the same geometry, its point at
    /// `(u, v)` the input's at `(u + angle, v)`. The input itself (a copy) when the
    /// seam stays where it is.
    pub surface: Surface,
    /// The seam about the axis (STEP u).
    pub u: ClosedSeam,
    /// For a torus closed round its tube, the seam round the tube (STEP v): the
    /// face's v lies within `[angle - 360, angle]`. `None` for the other surfaces.
    pub v: Option<ClosedSeam>,
}

/// `surface` placed so that its seam lies in the widest gap of the face bounded by
/// `loops`, or through `seam_point` where the boundary covers the full turn; `tol` is
/// the face's 3D tolerance.
///
/// **Closed surfaces.** A cone or cylinder, sphere, torus or surface of revolution
/// whose STEP u domain is a full turn: u is the angle in degrees about the Z axis of
/// the surface's position, from its X axis, or for a surface of revolution from its
/// profile (which the kernel keeps on the position's X/Z plane, on either side of the
/// axis; a profile off such a plane is not taken). A torus whose STEP v domain is a
/// full turn is closed round its tube too: v is the tube angle in degrees from the
/// outer equator.
///
/// **The boundary's angles.** Each edge of `loops` (their order and direction do not
/// matter) covers an arc of angles about the axis, found exactly from the edge's
/// extrema: a line from its ends; a circle or an ellipse from its ends and the
/// parameters where its angle about the axis turns back, the roots in closed form of
/// `(P(t) - O) x P'(t) . Z = 0`, unwound along it step by step (a monotone arc that
/// surrounds the axis by its turning sense). Points within `tol` of the axis (a
/// cone's apex, a sphere's poles) have no angle and are left out; an edge that lies
/// within `tol` of a plane through the axis (a ruling, a meridian, a profile) covers
/// the angle of its points on each side of the axis. Round a torus's tube, a circle
/// about the axis covers its one tube angle and a circle round the tube the arc it
/// runs over.
///
/// **The seam.** It goes in the middle of the widest gap between the covered arcs
/// (the first of equally wide ones, from the smallest angle). The surface is turned
/// about its own axis by that angle ([`Surface::transform`] with the rigid motion
/// that takes its position onto the turned one, so every kind keeps its own data):
/// the new X axis is `cos a X + sin a Y`, with the sine and cosine taken in degrees
/// after an exact reduction to the octant, so a seam on a quarter line gets exact
/// axes. When the boundary covers the full turn (the widest gap is no wider than
/// `tol` along the face's widest circle, or there is no boundary at all), the seam
/// stays where it is, unless `seam_point` is given: then the new X axis is exactly
/// the unit radial direction of `seam_point`, which normally is a vertex of the
/// face's seam edge. A seam point is not used where the boundary leaves a gap.
///
/// **Round a torus's tube.** STEP measures the tube angle from the outer equator,
/// and the kernel's torus keeps that (its evaluation has no tube start), so the seam
/// in v is no placement: [`PlacedSeam::v`] gives its angle, and the face's v lies
/// within `[angle - 360, angle]` (for [`crate::surface::torus::Torus::adjust_step_uv_domain`]).
/// A sphere's v (latitude) and the v of the other surfaces are not closed.
///
/// # Errors
///
/// [`SolidError::InvalidInput`] for a `tol` that is negative or not finite, a
/// `seam_point` that is not finite or (where it is used) within `tol` of the axis,
/// or of a torus's tube circle, an edge interval that is not finite or (on a conic)
/// longer than two turns, a surface that is not one of the closed kinds or is not
/// closed about its axis, a surface of revolution whose profile's middle point is
/// within `tol` of the axis. [`SolidError::Unsupported`] for an edge that is not a
/// line, circle or ellipse, round a torus's tube for an edge that is not a circle
/// about the axis or round the tube, and for a surface of revolution whose profile is
/// not on a plane through its axis. Those of [`crate::curve::curve::Curve::evaluate_step_point`] and
/// [`Surface::transform`].
pub fn seam_into_widest_gap(
    k: &Kernel,
    surface: &Surface,
    loops: &[Vec<BoundaryEdge<'_>>],
    seam_point: Option<Point3>,
    tol: f64,
) -> Result<PlacedSeam> {
    if !(tol.is_finite() && tol >= 0.0) {
        return Err(SolidError::invalid_input(
            OP,
            "tolerance negative or not finite",
        ));
    }
    if seam_point.is_some_and(|p| !(p.x.is_finite() && p.y.is_finite() && p.z.is_finite())) {
        return Err(SolidError::invalid_input(OP, "seam point not finite"));
    }
    let closed = closed(k, surface, tol)?;
    let axis = &closed.frame;

    let mut about = Cover::default();
    let mut round = Cover::default();
    for edge in loops.iter().flatten() {
        if !(edge.interval.min().is_finite() && edge.interval.max().is_finite()) {
            return Err(SolidError::invalid_input(OP, "edge interval not finite"));
        }
        about_axis(k, edge, axis, tol, &mut about)?;
        if let Some(major) = closed.tube {
            round_tube(k, edge, axis, major, tol, &mut round)?;
        }
    }

    // About the axis, with the seam point's radial direction where it is used.
    let (u, through) = match widest_gap(&about, tol) {
        Some(gap) => (gap_seam(gap), None),
        None => match seam_point {
            None => (
                ClosedSeam {
                    angle: None,
                    gap: None,
                },
                None,
            ),
            Some(p) => {
                let q = in_plane(axis, p);
                if q.length() <= tol {
                    return Err(SolidError::invalid_input(OP, "seam point on the axis"));
                }
                let direction = (axis.x_axis() * q.x + axis.y_axis() * q.y).unitize()?;
                let direction = direction * closed.side;
                (
                    ClosedSeam {
                        angle: Some(angle_of(q)),
                        gap: None,
                    },
                    Some(direction),
                )
            }
        },
    };

    let v = match closed.tube {
        None => None,
        Some(major) => Some(match widest_gap(&round, tol) {
            Some(gap) => gap_seam(gap),
            None => ClosedSeam {
                angle: match seam_point {
                    None => None,
                    Some(p) => {
                        let g = in_tube(axis, major, p);
                        if g.length() <= tol {
                            return Err(SolidError::invalid_input(
                                OP,
                                "seam point on the centre circle of the torus's tube",
                            ));
                        }
                        Some(angle_of(g))
                    }
                },
                gap: None,
            },
        }),
    };

    let mut placed = surface.clone();
    if let Some(angle) = u.angle {
        let x = through.unwrap_or_else(|| {
            let (sin, cos) = sin_cos_deg(angle);
            (axis.x_axis() * cos + axis.y_axis() * sin) * closed.side
        });
        let position = &closed.position;
        let turned = Axis2Placement::new(position.origin(), x, position.z_axis().cross(x));
        placed.transform(k, &Axis2Placement::from_to(position, &turned), None)?;
    }
    Ok(PlacedSeam {
        surface: placed,
        u,
        v,
    })
}

/// A closed surface's position, the frame its STEP u is measured in (the position's
/// origin and Z axis, u = 0 along its X axis, or for a surface of revolution whose
/// profile lies on the other side of the axis, along -X: `side` -1), and for a torus
/// closed round its tube, its major radius.
struct Closed {
    position: Axis2Placement,
    frame: Axis2Placement,
    side: f64,
    tube: Option<f64>,
}

fn closed(k: &Kernel, surface: &Surface, tol: f64) -> Result<Closed> {
    let (rev, major) = match surface.kind() {
        SurfaceKind::Cone(c) => (c.rev(), None),
        SurfaceKind::Sphere(s) => (s.base(), None),
        SurfaceKind::Torus(t) => (t.base(), Some(t.major_radius())),
        SurfaceKind::SurfOfRevolution(r) => (r.as_ref(), None),
        _ => {
            return Err(SolidError::invalid_input(
                OP,
                "not a cone, cylinder, sphere, torus or surface of revolution",
            ))
        }
    };
    let domain = rev.anal_uv_domain;
    if domain.u_interval().length() < TURN {
        return Err(SolidError::invalid_input(
            OP,
            "the surface is not closed about its axis",
        ));
    }
    let position = rev.position;
    let side = if let SurfaceKind::SurfOfRevolution(r) = surface.kind() {
        // u = 0 is the profile, which lies on the position's X/Z plane, on one side.
        if !r.planar_generator() {
            return Err(SolidError::Unsupported {
                operation: OP,
                reason: "a surface of revolution whose profile is not on a plane through its axis",
            });
        }
        let middle = surface.evaluate_step_point(k, Point2::new(0.0, domain.v_interval().mid()))?;
        let x = (middle - position.origin()).dot(position.x_axis());
        if x.abs() <= tol {
            return Err(SolidError::invalid_input(
                OP,
                "the profile's middle point is on the axis",
            ));
        }
        x.signum()
    } else {
        1.0
    };
    let frame = Axis2Placement::new(
        position.origin(),
        position.x_axis() * side,
        position.y_axis() * side,
    );
    Ok(Closed {
        position,
        frame,
        side,
        tube: major.filter(|_| domain.v_interval().length() >= TURN),
    })
}

/// The arcs of angles (degrees, unwound: `min <= max`) a boundary covers, and the
/// largest distance of their points from the axis (or from the tube's centre).
#[derive(Default)]
struct Cover {
    arcs: Vec<Extent1d>,
    reach: f64,
}

impl Cover {
    /// The angle of `q` (a point not within `tol` of the centre), alone.
    fn point(&mut self, q: Vec2) {
        let a = angle_of(q);
        self.arcs.push(Extent1d::from_min_max(a, a));
        self.reach = self.reach.max(q.length());
    }

    /// The arc swept by the key points `qs` (each with its parameter), in order:
    /// each point at its own angle, unwound by `turn(previous, next)`, the turn from
    /// one key point to the next. A point within `tol` of the centre ends the arc
    /// before it (it has no angle).
    fn sweep(
        &mut self,
        qs: &[(f64, Vec2)],
        tol: f64,
        turn: impl Fn((f64, Vec2), (f64, Vec2)) -> f64,
    ) {
        let mut run: Option<(f64, f64, f64)> = None;
        let mut previous: Option<(f64, Vec2)> = None;
        for &(t, q) in qs {
            let r = q.length();
            if r <= tol {
                self.end(run.take());
                previous = None;
                continue;
            }
            let a = angle_of(q);
            run = Some(match (run, previous) {
                (Some((at, lo, hi)), Some(prev)) => {
                    let estimate = at + turn(prev, (t, q));
                    let unwound = a + TURN * ((estimate - a) / TURN).round();
                    (unwound, lo.min(unwound), hi.max(unwound))
                }
                _ => (a, a, a),
            });
            previous = Some((t, q));
            self.reach = self.reach.max(r);
        }
        self.end(run);
    }

    fn end(&mut self, run: Option<(f64, f64, f64)>) {
        if let Some((_, lo, hi)) = run {
            self.arcs.push(Extent1d::from_min_max(lo, hi));
        }
    }
}

/// The covered arcs of an edge about the axis of `axis`.
fn about_axis(
    k: &Kernel,
    edge: &BoundaryEdge<'_>,
    axis: &Axis2Placement,
    tol: f64,
    cover: &mut Cover,
) -> Result<()> {
    let (lo, hi) = (edge.interval.min(), edge.interval.max());
    let at = |t: f64| -> Result<Vec2> { Ok(in_plane(axis, edge.curve.evaluate_step_point(k, t)?)) };
    match edge.curve.kind() {
        CurveKind::Line(_) => {
            let (q0, q1) = (at(lo)?, at(hi)?);
            let chord = (q1 - q0).length();
            if chord <= tol || cross(q0, q1).abs() <= tol * chord {
                // Along a ruling or in a plane through the axis: one angle each side.
                for q in [q0, q1] {
                    if q.length() > tol {
                        cover.point(q);
                    }
                }
            } else {
                // Off the axis, a segment turns less than half a turn about it.
                cover.sweep(&[(lo, q0), (hi, q1)], tol, |a, b| signed_angle(a.1, b.1));
            }
            Ok(())
        }
        CurveKind::Ellipse(e) => {
            if hi - lo > MAX_SPAN {
                return Err(SolidError::invalid_input(
                    OP,
                    "a conic edge longer than two turns",
                ));
            }
            let (c, a, b) = conic_in_plane(axis, e);
            match through_axis(c, a, b, tol) {
                Through::Near => Ok(()),
                Through::Plane(w) => {
                    // In a plane through the axis: the points farthest each side of it,
                    // at the ends or where `w . P(t)`, a sinusoid in t, turns.
                    let turning = rad_to_deg(w.dot(b).atan2(w.dot(a)));
                    let mut ts = vec![lo, hi];
                    for base in [turning, turning + 0.5 * TURN] {
                        push_periodic(&mut ts, base, lo, hi);
                    }
                    let mut far: Option<(f64, Vec2)> = None;
                    let mut near: Option<(f64, Vec2)> = None;
                    for t in ts {
                        let q = at(t)?;
                        let f = w.dot(q);
                        if far.is_none_or(|(g, _)| f > g) {
                            far = Some((f, q));
                        }
                        if near.is_none_or(|(g, _)| f < g) {
                            near = Some((f, q));
                        }
                    }
                    if let Some((_, q)) = far.filter(|(f, _)| *f > tol) {
                        cover.point(q);
                    }
                    if let Some((_, q)) = near.filter(|(f, _)| *f < -tol) {
                        cover.point(q);
                    }
                    Ok(())
                }
                Through::No => {
                    // The turn about the axis along `c + a cos t + b sin t` is
                    // `h(t) = D + P cos t + Q sin t` (per radian, over |P|^2).
                    let d = cross(a, b);
                    let (p, q) = (cross(c, b), cross(a, c));
                    let rho = sqrt(p * p + q * q);
                    let surrounds = d.abs() > rho;
                    let mut ts = steps(lo, hi);
                    if !surrounds && rho > 0.0 {
                        // Where the angle turns back: cos(t - phi) = -D / rho.
                        let phi = rad_to_deg(q.atan2(p));
                        let delta = rad_to_deg(acos(-d / rho));
                        for base in [phi - delta, phi + delta] {
                            push_periodic(&mut ts, base, lo, hi);
                        }
                        ts.sort_by(f64::total_cmp);
                        ts.dedup();
                    }
                    let mut qs = Vec::new();
                    for &t in &ts {
                        qs.push((t, at(t)?));
                    }
                    if surrounds {
                        // The arc surrounds the axis: its angle turns one way, with the
                        // sign of D. The affine map taking a, b to the unit axes takes it
                        // to the unit circle about m (|m| < 1) and keeps half turns, and
                        // there a quarter step turns at most 225 degrees: whether a step
                        // turns past half a turn is read there, its angle here.
                        let sense = d.signum();
                        let m = Vec2::new(p / d, q / d);
                        let image = |t: f64| {
                            let (sin, cos) = sin_cos_deg(t);
                            Vec2::new(m.x + cos, m.y + sin)
                        };
                        cover.sweep(&qs, tol, |(ta, qa), (tb, qb)| {
                            let (ia, ib) = (image(ta), image(tb));
                            let past_half = cross(ia, ib) < 0.0 && ia.dot(ib) < 0.0;
                            let turn = sense * signed_angle(qa, qb);
                            let turn = if past_half {
                                if turn < 0.25 * TURN {
                                    turn + TURN
                                } else {
                                    turn
                                }
                            } else if turn < -0.25 * TURN {
                                turn + TURN
                            } else {
                                turn
                            };
                            sense * turn
                        });
                    } else {
                        // The axis lies outside it: it turns less than half a turn.
                        cover.sweep(&qs, tol, |a, b| signed_angle(a.1, b.1));
                    }
                    Ok(())
                }
            }
        }
        _ => Err(SolidError::Unsupported {
            operation: OP,
            reason: "the angles about the axis of an edge that is not a line, circle or ellipse",
        }),
    }
}

/// The covered arcs of an edge round a torus's tube (tube angles).
fn round_tube(
    k: &Kernel,
    edge: &BoundaryEdge<'_>,
    axis: &Axis2Placement,
    major: f64,
    tol: f64,
    cover: &mut Cover,
) -> Result<()> {
    let CurveKind::Ellipse(e) = edge.curve.kind() else {
        return Err(SolidError::Unsupported {
            operation: OP,
            reason: "the tube angles of an edge that is not a circle",
        });
    };
    if !e.is_circle() {
        return Err(SolidError::Unsupported {
            operation: OP,
            reason: "the tube angles of an edge that is not a circle",
        });
    }
    let (lo, hi) = (edge.interval.min(), edge.interval.max());
    let at = |t: f64| -> Result<Vec2> {
        Ok(in_tube(axis, major, edge.curve.evaluate_step_point(k, t)?))
    };
    let position = e.position();
    let (c, a, b) = conic_in_plane(axis, e);
    let tilt = position.z_axis().cross(axis.z_axis()).length() * e.x_radius();
    if c.length() <= tol && tilt <= tol {
        // A circle about the axis: one tube angle.
        let g = at(lo)?;
        if g.length() > tol {
            cover.point(g);
        }
        return Ok(());
    }
    if !matches!(through_axis(c, a, b, tol), Through::Plane(_)) {
        return Err(SolidError::Unsupported {
            operation: OP,
            reason: "the tube angles of a circle neither about the axis nor round the tube",
        });
    }
    // A circle round the tube: its tube angle runs with its own, a quarter a step.
    let mut qs = Vec::new();
    for t in steps(lo, hi) {
        qs.push((t, at(t)?));
    }
    cover.sweep(&qs, tol, |a, b| signed_angle(a.1, b.1));
    Ok(())
}

/// How a conic `c + a cos t + b sin t` (in the plane normal to the axis) lies to the
/// axis.
enum Through {
    /// Within `tol` of the axis everywhere: no angle.
    Near,
    /// Within `tol` of a line through the axis, along the unit `w`.
    Plane(Vec2),
    /// Neither.
    No,
}

fn through_axis(c: Vec2, a: Vec2, b: Vec2, tol: f64) -> Through {
    let longest = [c, a, b]
        .into_iter()
        .max_by(|x, y| x.length().total_cmp(&y.length()))
        .unwrap_or(c);
    let length = longest.length();
    // (each point is within |c| + |a| + |b| of the axis)
    if c.length() + a.length() + b.length() <= tol {
        return Through::Near;
    }
    let w = Vec2::new(longest.x / length, longest.y / length);
    let normal = Vec2::new(-w.y, w.x);
    let (na, nb) = (normal.dot(a), normal.dot(b));
    if normal.dot(c).abs() + sqrt(na * na + nb * nb) <= tol {
        Through::Plane(w)
    } else {
        Through::No
    }
}

/// A conic's centre and axis vectors (radii applied) in the plane normal to `axis`.
fn conic_in_plane(axis: &Axis2Placement, e: &Ellipse) -> (Vec2, Vec2, Vec2) {
    let position = e.position();
    (
        in_plane(axis, position.origin()),
        along_plane(axis, position.x_axis() * e.x_radius()),
        along_plane(axis, position.y_axis() * e.y_radius()),
    )
}

/// A point's coordinates on the X and Y axes of `axis`, from its origin.
fn in_plane(axis: &Axis2Placement, p: Point3) -> Vec2 {
    along_plane(axis, p - axis.origin())
}

/// A vector's coordinates on the X and Y axes of `axis`.
fn along_plane(axis: &Axis2Placement, v: Point3) -> Vec2 {
    Vec2::new(v.dot(axis.x_axis()), v.dot(axis.y_axis()))
}

/// A point's place in its torus's tube: its distance from the axis less the major
/// radius, and its height (the tube angle is its angle).
fn in_tube(axis: &Axis2Placement, major: f64, p: Point3) -> Vec2 {
    let d = p - axis.origin();
    let q = along_plane(axis, d);
    Vec2::new(q.length() - major, d.dot(axis.z_axis()))
}

fn cross(a: Vec2, b: Vec2) -> f64 {
    a.x * b.y - a.y * b.x
}

/// The turn from `a` to `b` about the origin, in degrees in `(-180, 180]`.
fn signed_angle(a: Vec2, b: Vec2) -> f64 {
    rad_to_deg(cross(a, b).atan2(a.dot(b)))
}

/// `a` in `[0, 360)`.
fn wrap(a: f64) -> f64 {
    let r = a.rem_euclid(TURN);
    if r >= TURN {
        0.0
    } else {
        r
    }
}

/// The angle of `q` (not zero) in degrees in `[0, 360)`: the arc tangent of the
/// smaller coordinate over the larger, placed in its octant, so the angle of a point
/// on an axis is an exact multiple of 90.
fn angle_of(q: Vec2) -> f64 {
    let (ax, ay) = (q.x.abs(), q.y.abs());
    let a = if ay <= ax {
        rad_to_deg((ay / ax).atan())
    } else {
        0.25 * TURN - rad_to_deg((ax / ay).atan())
    };
    wrap(if q.x >= 0.0 {
        if q.y >= 0.0 {
            a
        } else {
            TURN - a
        }
    } else if q.y >= 0.0 {
        0.5 * TURN - a
    } else {
        0.5 * TURN + a
    })
}

/// `(sin a, cos a)` of `a` in degrees: `a` reduced exactly to within 45 degrees of a
/// quarter line, so a quarter line gets exact values.
fn sin_cos_deg(a: f64) -> (f64, f64) {
    let a = wrap(a);
    let quarter = (a / STEP).round();
    // (exact: a and the quarter line are within a factor of two, or a < 45)
    let r = a - STEP * quarter;
    let (s, c) = (r * (PI / 180.0)).sin_cos();
    if quarter == 1.0 {
        (c, -s)
    } else if quarter == 2.0 {
        (-s, -c)
    } else if quarter == 3.0 {
        (-c, s)
    } else {
        (s, c)
    }
}

/// `lo`, `hi` and the parameters between them a quarter turn apart at most.
fn steps(lo: f64, hi: f64) -> Vec<f64> {
    let n = ((hi - lo) / STEP).ceil().max(1.0);
    let mut ts = Vec::new();
    let mut i = 0.0;
    while i < n {
        ts.push(lo + (hi - lo) * (i / n));
        i += 1.0;
    }
    ts.push(hi);
    ts
}

/// Adds `base + 360 j` for every `j` that puts it strictly inside `(lo, hi)`.
fn push_periodic(ts: &mut Vec<f64>, base: f64, lo: f64, hi: f64) {
    if !base.is_finite() {
        return;
    }
    let mut t = base + TURN * ((lo - base) / TURN).floor();
    // (at most three turns: the interval is at most two long)
    let mut count = 0;
    while t < hi && count < 4 {
        if t > lo {
            ts.push(t);
        }
        t += TURN;
        count += 1;
    }
}

/// The widest gap between the arcs of `cover`, or `None` when they leave none wider
/// than `tol` along the farthest of their points (or there are none). The arcs are
/// laid out over three turns from the smallest start; there, on the second turn, the
/// covered set is exact, and the gaps starting on it are each gap once.
fn widest_gap(cover: &Cover, tol: f64) -> Option<Extent1d> {
    if cover.arcs.iter().any(|a| a.length() >= TURN) {
        return None;
    }
    let mut spans: Vec<(f64, f64)> = Vec::new();
    for arc in &cover.arcs {
        let start = wrap(arc.min());
        for turn in [0.0, TURN, 2.0 * TURN] {
            spans.push((start + turn, start + turn + arc.length()));
        }
    }
    spans.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.total_cmp(&y.1)));
    let first = spans.first()?.0;
    let mut blocks: Vec<(f64, f64)> = Vec::new();
    for (lo, hi) in spans {
        match blocks.last_mut() {
            Some(last) if lo <= last.1 => last.1 = last.1.max(hi),
            _ => blocks.push((lo, hi)),
        }
    }
    let mut best: Option<(f64, f64)> = None;
    for pair in blocks.windows(2) {
        if let [before, after] = pair {
            let start = before.1;
            let width = after.0 - start;
            if start >= first + TURN
                && start < first + 2.0 * TURN
                && best.is_none_or(|(_, w)| width > w)
            {
                best = Some((start, width));
            }
        }
    }
    let (start, width) = best?;
    if width * (PI / 180.0) * cover.reach <= tol {
        return None;
    }
    let start = wrap(start - TURN);
    Some(Extent1d::from_min_max(start, start + width))
}

/// The seam in the middle of `gap`.
fn gap_seam(gap: Extent1d) -> ClosedSeam {
    ClosedSeam {
        angle: Some(wrap(gap.min() + 0.5 * gap.length())),
        gap: Some(gap),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curve::curve::Curve;
    use crate::curve::ellipse::Ellipse;
    use crate::curve::line::Line;
    use crate::geom::extent2d::Extent2d;
    use crate::math::{Point2, Vec3};
    use crate::nurbs::b_spline_curve::{BSplineCurve, BSplineCurveForm, CanonicalData, KnotType};
    use crate::surface::cone::Cone;
    use crate::surface::plane::Plane;
    use crate::surface::sphere::Sphere;
    use crate::surface::surf_of_revolution::SurfOfRevolution;
    use crate::surface::surface::SurfaceKind;
    use crate::surface::torus::{Torus, TorusOptions};
    use std::f64::consts::FRAC_1_SQRT_2;

    /// The face tolerance of the cases (a usual STEP distance uncertainty).
    const TOL: f64 = 1e-6;

    fn k() -> Kernel {
        Kernel::default()
    }

    fn p(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z)
    }

    fn v3(x: f64, y: f64, z: f64) -> Vec3 {
        Vec3::new(x, y, z)
    }

    /// The circle of `radius` about `centre` in the plane of `x` and `y`, the full
    /// curve (its edges take their intervals from [`edge`]).
    fn circle(k: &Kernel, centre: Point3, x: Vec3, y: Vec3, radius: f64) -> Curve {
        Curve::from(
            Ellipse::create_canonical_circle(
                k,
                &Axis2Placement::new(centre, x, y),
                radius,
                None,
                None,
            )
            .unwrap(),
        )
    }

    /// The line from `a` to `b` over STEP parameters `[0, 1]`.
    fn line(a: Point3, b: Point3) -> Curve {
        Curve::from(Line::create_canonical(a, b - a))
    }

    fn edge(curve: &Curve, lo: f64, hi: f64) -> BoundaryEdge<'_> {
        BoundaryEdge {
            curve,
            interval: Extent1d::from_min_max(lo, hi),
            forward: true,
        }
    }

    fn position(s: &Surface) -> Axis2Placement {
        match s.kind() {
            SurfaceKind::Cone(c) => *c.rev().position(),
            SurfaceKind::Sphere(s) => *s.base().position(),
            SurfaceKind::Torus(t) => *t.base().position(),
            SurfaceKind::SurfOfRevolution(r) => *r.position(),
            _ => panic!("not a closed surface"),
        }
    }

    /// The placed surface has the input's points: its point at `(u, v)` is the input's
    /// at `(u + turn, v)`, for every eighth of a turn (and between) and each `v`. Both
    /// are evaluated by their own `cos`/`sin` of the angle in radians (each within an
    /// ulp) from placements that agree to rounding, so the points agree within a few
    /// ulps of the coordinates' size `scale`: the bound is 8 ulps of it.
    #[track_caller]
    fn assert_same_points(
        k: &Kernel,
        old: &Surface,
        new: &Surface,
        turn: f64,
        vs: &[f64],
        scale: f64,
    ) {
        let bound = 8.0 * f64::EPSILON * scale;
        for i in 0..16 {
            let u = 22.5 * f64::from(i) + 3.25;
            for &v in vs {
                let a = new.evaluate_step_point(k, Point2::new(u, v)).unwrap();
                let b = old
                    .evaluate_step_point(k, Point2::new((u + turn).rem_euclid(360.0), v))
                    .unwrap();
                assert!(
                    a.distance(b) <= bound,
                    "at ({u}, {v}): {a:?} is not {b:?} (distance {}, bound {bound})",
                    a.distance(b)
                );
            }
        }
    }

    /// `a` equals the closed form `b` to rounding (two ulps of 360 degrees).
    #[track_caller]
    fn assert_angle(a: f64, b: f64) {
        assert!((a - b).abs() <= 1e-13, "{a} is not {b}");
    }

    /// The vector `a` is the closed form `b` to rounding (two ulps of 1).
    #[track_caller]
    fn assert_unit(a: Vec3, b: Vec3) {
        assert!(a.distance(b) <= 2.0 * f64::EPSILON, "{a:?} is not {b:?}");
    }

    fn cylinder(k: &Kernel, frame: &Axis2Placement, r: f64) -> Surface {
        Surface::from(Cone::create_canonical_cylinder(k, frame, r).unwrap())
    }

    /// The face of a cylinder of radius 2 about Z between the heights 0 and 3 and the
    /// angles 0 and 180 degrees (its corner arcs split at 90): boundary points at 0,
    /// 90 and 180 degrees.
    struct HalfCylinder {
        arcs: [Curve; 2],
        lines: [Curve; 2],
    }

    impl HalfCylinder {
        fn new(k: &Kernel, frame: &Axis2Placement) -> Self {
            let (x, y, z) = (frame.x_axis(), frame.y_axis(), frame.z_axis());
            let o = frame.origin();
            Self {
                arcs: [circle(k, o, x, y, 2.0), circle(k, o + 3.0 * z, x, y, 2.0)],
                lines: [
                    line(o + 2.0 * x, o + 2.0 * x + 3.0 * z),
                    line(o - 2.0 * x, o - 2.0 * x + 3.0 * z),
                ],
            }
        }

        fn loops(&self) -> Vec<Vec<BoundaryEdge<'_>>> {
            let [bottom, top] = &self.arcs;
            let [at0, at180] = &self.lines;
            vec![vec![
                edge(bottom, 0.0, 90.0),
                edge(bottom, 90.0, 180.0),
                edge(at180, 0.0, 1.0),
                edge(top, 90.0, 180.0),
                edge(top, 0.0, 90.0),
                edge(at0, 0.0, 1.0),
            ]]
        }
    }

    #[test]
    fn cylinder_seam_into_its_widest_gap() {
        let k = k();
        let frame = Axis2Placement::IDENTITY;
        let s = cylinder(&k, &frame, 2.0);
        let face = HalfCylinder::new(&k, &frame);
        let placed = seam_into_widest_gap(&k, &s, &face.loops(), None, TOL).unwrap();

        // The boundary covers 0..180 degrees: the gap is 180..360, its middle 270.
        assert_eq!(placed.u.angle, Some(270.0));
        let gap = placed.u.gap.unwrap();
        assert_eq!((gap.min(), gap.length()), (180.0, 180.0));
        assert!(placed.v.is_none());

        // The new X axis is exactly the direction at 270 degrees.
        let pos = position(&placed.surface);
        assert_eq!(pos.x_axis(), v3(0.0, -1.0, 0.0));
        assert_eq!(pos.y_axis(), v3(1.0, 0.0, 0.0));
        assert_eq!(pos.z_axis(), v3(0.0, 0.0, 1.0));
        assert_eq!(pos.origin(), p(0.0, 0.0, 0.0));
        assert_same_points(&k, &s, &placed.surface, 270.0, &[-1.0, 0.0, 1.5, 3.0], 4.0);
    }

    #[test]
    fn cylinder_seam_in_a_general_frame() {
        // The same face on a cylinder whose axis is (4/5, -3/5, 0) through (1, 2, 3).
        let k = k();
        let frame = Axis2Placement::new(p(1.0, 2.0, 3.0), v3(0.6, 0.8, 0.0), v3(0.0, 0.0, 1.0));
        let s = cylinder(&k, &frame, 2.0);
        let face = HalfCylinder::new(&k, &frame);
        let placed = seam_into_widest_gap(&k, &s, &face.loops(), None, TOL).unwrap();

        assert_angle(placed.u.angle.unwrap(), 270.0);
        let pos = position(&placed.surface);
        // At 270 degrees: -Y of the frame.
        assert_unit(pos.x_axis(), v3(0.0, 0.0, -1.0));
        assert_unit(pos.z_axis(), frame.z_axis());
        assert_same_points(&k, &s, &placed.surface, 270.0, &[-1.0, 0.0, 1.5, 3.0], 8.0);
    }

    #[test]
    fn cylinder_full_band_keeps_its_seam() {
        let k = k();
        let frame = Axis2Placement::IDENTITY;
        let s = cylinder(&k, &frame, 2.0);
        let bottom = circle(&k, p(0.0, 0.0, 0.0), Vec3::X, Vec3::Y, 2.0);
        let top = circle(&k, p(0.0, 0.0, 3.0), Vec3::X, Vec3::Y, 2.0);
        let loops = vec![
            vec![edge(&bottom, 0.0, 360.0)],
            vec![edge(&top, 0.0, 360.0)],
        ];
        let placed = seam_into_widest_gap(&k, &s, &loops, None, TOL).unwrap();

        assert_eq!(
            placed.u,
            ClosedSeam {
                angle: None,
                gap: None
            }
        );
        assert_eq!(position(&placed.surface), frame);
        assert_same_points(&k, &s, &placed.surface, 0.0, &[0.0, 3.0], 4.0);

        // The circles as two half edges each: their ends meet to rounding, still the
        // full turn.
        let halves = vec![
            vec![edge(&bottom, 0.0, 180.0), edge(&bottom, 180.0, 360.0)],
            vec![edge(&top, 0.0, 180.0), edge(&top, 180.0, 360.0)],
        ];
        let placed = seam_into_widest_gap(&k, &s, &halves, None, TOL).unwrap();
        assert_eq!(
            placed.u,
            ClosedSeam {
                angle: None,
                gap: None
            }
        );
    }

    #[test]
    fn cylinder_full_band_seam_through_its_seam_point() {
        // A band of radius 5 with a seam vertex at (3, 4, 1): the X axis is exactly
        // its radial direction (3/5, 4/5, 0).
        let k = k();
        let frame = Axis2Placement::IDENTITY;
        let s = cylinder(&k, &frame, 5.0);
        let bottom = circle(&k, p(0.0, 0.0, 0.0), Vec3::X, Vec3::Y, 5.0);
        let top = circle(&k, p(0.0, 0.0, 3.0), Vec3::X, Vec3::Y, 5.0);
        let loops = vec![
            vec![edge(&bottom, 0.0, 360.0)],
            vec![edge(&top, 0.0, 360.0)],
        ];
        let placed = seam_into_widest_gap(&k, &s, &loops, Some(p(3.0, 4.0, 1.0)), TOL).unwrap();

        let angle = placed.u.angle.unwrap();
        assert_angle(angle, 4.0_f64.atan2(3.0).to_degrees());
        assert!(placed.u.gap.is_none());
        let pos = position(&placed.surface);
        assert_eq!(pos.x_axis(), v3(3.0 / 5.0, 4.0 / 5.0, 0.0));
        assert_unit(pos.y_axis(), v3(-4.0 / 5.0, 3.0 / 5.0, 0.0));
        assert_same_points(&k, &s, &placed.surface, angle, &[0.0, 3.0], 8.0);
    }

    #[test]
    fn cylinder_oblique_ellipse_across_the_seam() {
        // The cylinder of radius 1 cut by the plane z = x: the ellipse (cos t, sin t,
        // cos t) over -90..90 degrees (across the seam), down the lines at -90 and
        // 90 to the arc at z = -2. The boundary covers -90..90: the seam goes to 180.
        let k = k();
        let s = cylinder(&k, &Axis2Placement::IDENTITY, 1.0);
        let ellipse = Curve::from(
            Ellipse::create_canonical(
                &k,
                &Axis2Placement::new(p(0.0, 0.0, 0.0), v3(1.0, 0.0, 1.0), Vec3::Y),
                2.0_f64.sqrt(),
                1.0,
                None,
                None,
            )
            .unwrap(),
        );
        let bottom = circle(&k, p(0.0, 0.0, -2.0), Vec3::X, Vec3::Y, 1.0);
        let up = line(p(0.0, -1.0, -2.0), p(0.0, -1.0, 0.0));
        let down = line(p(0.0, 1.0, 0.0), p(0.0, 1.0, -2.0));
        let loops = vec![vec![
            edge(&bottom, -90.0, 90.0),
            edge(&down, 0.0, 1.0),
            edge(&ellipse, -90.0, 90.0),
            edge(&up, 0.0, 1.0),
        ]];
        let placed = seam_into_widest_gap(&k, &s, &loops, None, TOL).unwrap();

        assert_angle(placed.u.angle.unwrap(), 180.0);
        let gap = placed.u.gap.unwrap();
        assert_angle(gap.min(), 90.0);
        assert_angle(gap.length(), 180.0);
        assert_unit(position(&placed.surface).x_axis(), v3(-1.0, 0.0, 0.0));
        assert_same_points(
            &k,
            &s,
            &placed.surface,
            placed.u.angle.unwrap(),
            &[-2.0, 0.0, 1.0],
            4.0,
        );
    }

    #[test]
    fn cone_seam_with_its_apex_on_the_boundary() {
        // A cone of semi-angle 30 degrees and radius 1 at height 0: its apex at
        // height -sqrt(3). The face runs from the apex along the lines at 90 and 180
        // degrees to the arc 90..180 at height 0. The apex has no angle; the gap is
        // 180..450, its middle 315.
        let k = k();
        let s = Surface::from(
            Cone::create_canonical_with_semi_angle(
                &k,
                &Axis2Placement::IDENTITY,
                1.0,
                30.0,
                None,
                false,
                false,
            )
            .unwrap(),
        );
        let apex = p(0.0, 0.0, -(3.0_f64.sqrt()));
        let arc = circle(&k, p(0.0, 0.0, 0.0), Vec3::X, Vec3::Y, 1.0);
        let at90 = line(apex, p(0.0, 1.0, 0.0));
        let at180 = line(p(-1.0, 0.0, 0.0), apex);
        let loops = vec![vec![
            edge(&at90, 0.0, 1.0),
            edge(&arc, 90.0, 180.0),
            edge(&at180, 0.0, 1.0),
        ]];
        let placed = seam_into_widest_gap(&k, &s, &loops, None, TOL).unwrap();

        assert_angle(placed.u.angle.unwrap(), 315.0);
        let gap = placed.u.gap.unwrap();
        assert_angle(gap.min(), 180.0);
        assert_angle(gap.length(), 270.0);
        assert_unit(
            position(&placed.surface).x_axis(),
            v3(FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0),
        );
        assert_same_points(&k, &s, &placed.surface, 315.0, &[-1.0, 0.0, 2.0], 4.0);
    }

    fn sphere(k: &Kernel, r: f64) -> Surface {
        Surface::from(Sphere::create_canonical(k, &Axis2Placement::IDENTITY, r).unwrap())
    }

    #[test]
    fn sphere_lune_seam_into_its_longitude_gap() {
        // The lune of a sphere of radius 2 between the meridians at 0 and 90 degrees,
        // each from pole to pole: the poles have no longitude, the gap is 90..360.
        let k = k();
        let s = sphere(&k, 2.0);
        let o = p(0.0, 0.0, 0.0);
        let at0 = circle(&k, o, Vec3::X, Vec3::Z, 2.0);
        let at90 = circle(&k, o, Vec3::Y, Vec3::Z, 2.0);
        let loops = vec![vec![edge(&at0, -90.0, 90.0), edge(&at90, -90.0, 90.0)]];
        let placed = seam_into_widest_gap(&k, &s, &loops, None, TOL).unwrap();

        assert_angle(placed.u.angle.unwrap(), 225.0);
        let gap = placed.u.gap.unwrap();
        assert_angle(gap.min(), 90.0);
        assert_angle(gap.length(), 270.0);
        assert!(placed.v.is_none());
        assert_unit(
            position(&placed.surface).x_axis(),
            v3(-FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0),
        );
        assert_same_points(
            &k,
            &s,
            &placed.surface,
            225.0,
            &[-90.0, -45.0, 0.0, 30.0, 90.0],
            2.0,
        );
    }

    #[test]
    fn sphere_cap_at_the_pole() {
        // The cap above latitude 60 of a sphere of radius 2: its circle (radius 1 at
        // height sqrt(3)) covers every longitude, so the seam stays, or goes through
        // the seam point (0, 1, sqrt(3)) at 90 degrees.
        let k = k();
        let s = sphere(&k, 2.0);
        let h = 3.0_f64.sqrt();
        let rim = circle(&k, p(0.0, 0.0, h), Vec3::X, Vec3::Y, 1.0);
        let loops = vec![vec![edge(&rim, 0.0, 360.0)]];

        let kept = seam_into_widest_gap(&k, &s, &loops, None, TOL).unwrap();
        assert_eq!(
            kept.u,
            ClosedSeam {
                angle: None,
                gap: None
            }
        );
        assert_eq!(position(&kept.surface), Axis2Placement::IDENTITY);

        let through = seam_into_widest_gap(&k, &s, &loops, Some(p(0.0, 1.0, h)), TOL).unwrap();
        assert_eq!(through.u.angle, Some(90.0));
        assert_eq!(position(&through.surface).x_axis(), v3(0.0, 1.0, 0.0));
        assert_same_points(&k, &s, &through.surface, 90.0, &[60.0, 75.0, 90.0], 2.0);
    }

    fn torus(k: &Kernel) -> Surface {
        Surface::from(
            Torus::create_canonical(
                k,
                &Axis2Placement::IDENTITY,
                3.0,
                1.0,
                TorusOptions::default(),
                None,
            )
            .unwrap(),
        )
    }

    #[test]
    fn torus_seam_about_its_axis() {
        // Half a torus (major radius 3, minor 1): the tube circles at 0 and 180
        // degrees and the outer equator 0..180 (its seam round the tube, used
        // twice). About the axis the gap is 180..360; round the tube the boundary
        // covers the full turn, so that seam stays.
        let k = k();
        let s = torus(&k);
        let tube0 = circle(&k, p(3.0, 0.0, 0.0), Vec3::X, Vec3::Z, 1.0);
        let tube180 = circle(&k, p(-3.0, 0.0, 0.0), -Vec3::X, Vec3::Z, 1.0);
        let equator = circle(&k, p(0.0, 0.0, 0.0), Vec3::X, Vec3::Y, 4.0);
        let loops = vec![vec![
            edge(&equator, 0.0, 180.0),
            edge(&tube180, 0.0, 360.0),
            edge(&equator, 0.0, 180.0),
            edge(&tube0, 0.0, 360.0),
        ]];
        let placed = seam_into_widest_gap(&k, &s, &loops, None, TOL).unwrap();

        assert_eq!(placed.u.angle, Some(270.0));
        assert_eq!(position(&placed.surface).x_axis(), v3(0.0, -1.0, 0.0));
        assert_eq!(
            placed.v,
            Some(ClosedSeam {
                angle: None,
                gap: None
            })
        );
        assert_same_points(
            &k,
            &s,
            &placed.surface,
            270.0,
            &[0.0, 45.0, 180.0, 300.0],
            4.0,
        );
    }

    #[test]
    fn torus_seam_round_its_tube() {
        // The outer band of the torus between the tube angles -30 and 30 degrees: the
        // circles of radius 3 + cos 30 at heights -/+ 1/2 cover every angle about the
        // axis (that seam stays); round the tube the gap is 30..330, its middle 180
        // (the inner equator).
        let k = k();
        let s = torus(&k);
        let r = 3.0 + 3.0_f64.sqrt() / 2.0;
        let low = circle(&k, p(0.0, 0.0, -0.5), Vec3::X, Vec3::Y, r);
        let high = circle(&k, p(0.0, 0.0, 0.5), Vec3::X, Vec3::Y, r);
        let loops = vec![vec![edge(&low, 0.0, 360.0)], vec![edge(&high, 0.0, 360.0)]];
        let placed = seam_into_widest_gap(&k, &s, &loops, None, TOL).unwrap();

        assert_eq!(
            placed.u,
            ClosedSeam {
                angle: None,
                gap: None
            }
        );
        let v = placed.v.unwrap();
        assert_angle(v.angle.unwrap(), 180.0);
        let gap = v.gap.unwrap();
        assert_angle(gap.min(), 30.0);
        assert_angle(gap.length(), 300.0);
        assert_eq!(position(&placed.surface), Axis2Placement::IDENTITY);
    }

    #[test]
    fn torus_v_seam_through_its_seam_point() {
        // The torus's whole tube between the tube circles at 0 and 90 degrees about
        // the axis, with its seam vertex on the top circle (tube angle 90): round the
        // tube the seam goes through it.
        let k = k();
        let s = torus(&k);
        let tube0 = circle(&k, p(3.0, 0.0, 0.0), Vec3::X, Vec3::Z, 1.0);
        let tube90 = circle(&k, p(0.0, 3.0, 0.0), Vec3::Y, Vec3::Z, 1.0);
        let top = circle(&k, p(0.0, 0.0, 1.0), Vec3::X, Vec3::Y, 3.0);
        let loops = vec![vec![
            edge(&top, 0.0, 90.0),
            edge(&tube90, 0.0, 360.0),
            edge(&top, 0.0, 90.0),
            edge(&tube0, 0.0, 360.0),
        ]];
        let placed = seam_into_widest_gap(&k, &s, &loops, Some(p(3.0, 0.0, 1.0)), TOL).unwrap();

        assert_angle(placed.u.angle.unwrap(), 225.0);
        assert_eq!(
            placed.v,
            Some(ClosedSeam {
                angle: Some(90.0),
                gap: None
            })
        );
    }

    #[test]
    fn revolution_of_a_line_profile() {
        // The line from (2, 0, 0) to (3, 0, 2) about Z, the face between its copies
        // at 0 and 180 degrees and the arcs at heights 0 and 2: the gap is 180..360.
        let k = k();
        let profile = line(p(2.0, 0.0, 0.0), p(3.0, 0.0, 2.0));
        let s = Surface::from(
            SurfOfRevolution::create_canonical(
                &k,
                profile.clone(),
                p(0.0, 0.0, 0.0),
                Vec3::Z,
                false,
            )
            .unwrap(),
        );
        let at180 = line(p(-2.0, 0.0, 0.0), p(-3.0, 0.0, 2.0));
        let bottom = circle(&k, p(0.0, 0.0, 0.0), Vec3::X, Vec3::Y, 2.0);
        let top = circle(&k, p(0.0, 0.0, 2.0), Vec3::X, Vec3::Y, 3.0);
        let loops = vec![vec![
            edge(&bottom, 0.0, 180.0),
            edge(&at180, 0.0, 1.0),
            edge(&top, 0.0, 180.0),
            edge(&profile, 0.0, 1.0),
        ]];
        let placed = seam_into_widest_gap(&k, &s, &loops, None, TOL).unwrap();

        assert_eq!(placed.u.angle, Some(270.0));
        // The kernel places this revolution with X = (-1, 0, 0), its profile on -X: the
        // position turned by 270 degrees about Z has X = (0, 1, 0).
        assert_eq!(position(&s).x_axis(), v3(-1.0, 0.0, 0.0));
        assert_eq!(position(&placed.surface).x_axis(), v3(0.0, 1.0, 0.0));
        assert!(placed.v.is_none());
        assert_same_points(&k, &s, &placed.surface, 270.0, &[0.0, 0.5, 1.0], 4.0);
    }

    #[test]
    fn revolution_of_an_arc_profile() {
        // The half circle of radius 1 about (3, 0, 0) in the XZ plane, from -90 to 90
        // degrees, about Z; the face between its copies at 0 and 90 degrees and the
        // arcs of radius 3 at heights -1 and 1: the gap is 90..360, its middle 225.
        let k = k();
        let profile = circle(&k, p(3.0, 0.0, 0.0), Vec3::X, Vec3::Z, 1.0);
        let half = Curve::from(
            Ellipse::create_canonical_circle(
                &k,
                &Axis2Placement::new(p(3.0, 0.0, 0.0), Vec3::X, Vec3::Z),
                1.0,
                Some(&Extent1d::from_min_max(-90.0, 90.0)),
                None,
            )
            .unwrap(),
        );
        let s = Surface::from(
            SurfOfRevolution::create_canonical(&k, half, p(0.0, 0.0, 0.0), Vec3::Z, false).unwrap(),
        );
        let at90 = circle(&k, p(0.0, 3.0, 0.0), Vec3::Y, Vec3::Z, 1.0);
        let low = circle(&k, p(0.0, 0.0, -1.0), Vec3::X, Vec3::Y, 3.0);
        let high = circle(&k, p(0.0, 0.0, 1.0), Vec3::X, Vec3::Y, 3.0);
        let loops = vec![vec![
            edge(&low, 0.0, 90.0),
            edge(&at90, -90.0, 90.0),
            edge(&high, 0.0, 90.0),
            edge(&profile, -90.0, 90.0),
        ]];
        let placed = seam_into_widest_gap(&k, &s, &loops, None, TOL).unwrap();

        assert_angle(placed.u.angle.unwrap(), 225.0);
        assert_unit(
            position(&placed.surface).x_axis(),
            v3(-FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0),
        );
        assert_same_points(
            &k,
            &s,
            &placed.surface,
            225.0,
            &[-90.0, 0.0, 45.0, 90.0],
            4.0,
        );
    }

    #[test]
    fn a_face_without_boundary_keeps_its_seam() {
        // A whole sphere: no edges, every angle.
        let k = k();
        let s = sphere(&k, 1.0);
        let placed = seam_into_widest_gap(&k, &s, &[], None, TOL).unwrap();
        assert_eq!(
            placed.u,
            ClosedSeam {
                angle: None,
                gap: None
            }
        );
        assert_eq!(position(&placed.surface), Axis2Placement::IDENTITY);
    }

    #[test]
    fn invalid_inputs_are_refused() {
        let k = k();
        let frame = Axis2Placement::IDENTITY;
        let s = cylinder(&k, &frame, 2.0);
        let face = HalfCylinder::new(&k, &frame);
        let loops = face.loops();
        let invalid = |r: Result<PlacedSeam>| {
            assert!(
                matches!(r, Err(SolidError::InvalidInput { operation: OP, .. })),
                "{r:?}"
            );
        };

        // a tolerance negative or not finite
        invalid(seam_into_widest_gap(&k, &s, &loops, None, -1.0));
        invalid(seam_into_widest_gap(&k, &s, &loops, None, f64::NAN));
        // a seam point not finite, or on the axis where the band needs it
        let bottom = circle(&k, p(0.0, 0.0, 0.0), Vec3::X, Vec3::Y, 2.0);
        let band = vec![vec![edge(&bottom, 0.0, 360.0)]];
        invalid(seam_into_widest_gap(
            &k,
            &s,
            &band,
            Some(p(f64::INFINITY, 0.0, 0.0)),
            TOL,
        ));
        invalid(seam_into_widest_gap(
            &k,
            &s,
            &band,
            Some(p(0.0, 0.0, 1.0)),
            TOL,
        ));
        // an edge interval not finite
        let open = vec![vec![edge(&bottom, 0.0, f64::INFINITY)]];
        invalid(seam_into_widest_gap(&k, &s, &open, None, TOL));
        // a surface that is not closed
        let plane = Surface::from(Plane::create_canonical(&frame).unwrap());
        invalid(seam_into_widest_gap(&k, &plane, &loops, None, TOL));
        // a cone over less than a full turn
        let mut part = Cone::create_canonical_cylinder(&k, &frame, 2.0).unwrap();
        part.adjust_step_uv_domain(&k, &Extent2d::from_bounds(0.0, 0.0, 180.0, 3.0))
            .unwrap();
        invalid(seam_into_widest_gap(
            &k,
            &Surface::from(part),
            &loops,
            None,
            TOL,
        ));
    }

    #[test]
    fn edges_of_other_kinds_are_unsupported() {
        // A degree-2 spline edge on the cylinder: its angles are not computed here.
        let k = k();
        let s = cylinder(&k, &Axis2Placement::IDENTITY, 1.0);
        let points = [p(1.0, 0.0, 0.0), p(1.0, 1.0, 0.5), p(0.0, 1.0, 1.0)];
        let spline = Curve::from(
            BSplineCurve::create_canonical(&CanonicalData {
                dim: 3,
                degree: 2,
                control_points: &points,
                form: BSplineCurveForm::Unspecified,
                knot_multiplicities: &[3, 3],
                knots: &[0.0, 1.0],
                knot_type: KnotType::Unspecified,
                weights: Some(&[1.0, FRAC_1_SQRT_2, 1.0]),
                trim: None,
            })
            .unwrap(),
        );
        let loops = vec![vec![edge(&spline, 0.0, 1.0)]];
        let r = seam_into_widest_gap(&k, &s, &loops, None, TOL);
        assert!(
            matches!(r, Err(SolidError::Unsupported { operation: OP, .. })),
            "{r:?}"
        );
    }
}
