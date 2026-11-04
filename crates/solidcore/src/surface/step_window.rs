//! (not upstream) A face's STEP parameter window from its boundary: the box in the
//! surface's STEP (u, v) parameters that the face's bounding edges span on it, each
//! edge's extremes found exactly, unwound along the surface's periodic directions.
//!
//! Upstream: none (solidcore's STEP reader, #2578). The idea is SolidCraft's reader
//! (`crates/step/src/topology/swept.rs` `span`, `range`, `trimmed_range`,
//! `widest_gap`, `common_angle`, `crates/step/src/topology.rs` `Builder::extent` and
//! `face_extent`, and `crates/geom/src/brep/form.rs` `extent`, at a5b71ca2; MIT OR
//! Apache-2.0), which samples each edge; here each edge's extremes are its
//! stationary points. The computations are the kernel's own: the STEP
//! parameterisations of [`Surface::evaluate_step_point`] (in closed form for planes,
//! cones, cylinders, spheres and tori, and the angle of a surface of revolution),
//! [`Surface::step_inversion`] and [`Surface::evaluate_step`] for the other kinds and
//! a revolution's profile parameter, [`Curve::evaluate`] and the line's and
//! ellipse's STEP evaluators for the edges, and [`LocalSolve1d::solve_by_bisection`]
//! for the stationary points without a closed form.
//!
//! | Rust                 | SolidCraft                                                  |
//! |----------------------|-------------------------------------------------------------|
//! | [`face_window`]      | `form::extent`, `Builder::extent`, `face_extent`; `swept::span`, `range`, `trimmed_range`, `widest_gap`, `common_angle` (a swept surface's profile range is its window's v) |
//! | [`BoundaryEdge`]     | `form::extent`'s `(&Curve, [f64; 2], bool)` loop entries    |
//! | [`WindowMargin`]     | `swept::MARGIN` (and step-brep's 5 %): ops' choice now      |
//!
//! Deviations: not a port. Recursion: halving a piece of an edge to decide its turn
//! round a period stops at [`MAX_SPLIT_DEPTH`] levels.

use crate::curve::curve::Curve;
use crate::error::{Result, SolidError};
use crate::geom::extent1d::Extent1d;
use crate::geom::extent2d::Extent2d;
use crate::kernel::Kernel;
use crate::math::Point3;
use crate::surface::surface::Surface;

const OP: &str = "surface::step_window::face_window";

/// One edge of a face's boundary loop, as the loop uses it: the edge's curve, its
/// STEP interval on that curve and the direction the loop runs along it.
///
/// A face's loops are lists of these in loop order, each edge running on from where
/// the one before it ends. The loops are oriented as they bound the face: the face
/// lies on their left about the surface's STEP normal `Su x Sv` (a file's
/// `FACE_BOUND.orientation` applied, and every loop reversed when the face's
/// `same_sense` is false).
#[derive(Debug, Clone, Copy)]
pub struct BoundaryEdge<'a> {
    /// The edge's curve (for a `CIRCLE` or an `ELLIPSE` the full conic), whose STEP
    /// parameters [`Self::interval`] is in (those of [`Curve::evaluate_step_point`]).
    pub curve: &'a Curve,
    /// The edge's STEP interval on [`Self::curve`], lower parameter first: the
    /// [`crate::curve::step_edge::StepEdgeInterval::interval`] of
    /// [`crate::curve::step_edge::edge_interval`].
    pub interval: Extent1d,
    /// The loop runs along the curve's parameter, from the interval's minimum to its
    /// maximum; when false, from the maximum down to the minimum.
    pub forward: bool,
}

/// The margin [`face_window`] adds on each side of the face's box, in the surface's
/// STEP units per direction (degrees for an angle, lengths or the profile's own
/// parameter otherwise). Ops chooses it; the default is none.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct WindowMargin {
    /// Added below the box's u minimum and above its maximum.
    pub u: f64,
    /// Added below the box's v minimum and above its maximum.
    pub v: f64,
}

/// The STEP (u, v) window of a face on `surface`: the box its boundary spans, plus
/// `margin`, as [`Surface::adjust_step_uv_domain`] takes it.
///
/// `loops` are the face's edge loops ([`BoundaryEdge`]), the outer loop first;
/// `vertices` the points of its vertex loops (an apex or a pole); `tol` the distance
/// within which the boundary lies on the surface.
///
/// **Exact.** Between its stationary points an edge's u and v are monotone, so an
/// edge's extremes are at its ends or where `du/dt` or `dv/dt` vanishes. Those are
/// found in closed form for lines and conics (`ELLIPSE`, `CIRCLE`) on a plane (the
/// extremes of `(P - O) . X` and `(P - O) . Y`), on a cone or cylinder (the height)
/// and on a sphere (the latitude), and for the angle about the axis of every surface
/// of revolution (`((P - O) x P') . Z = 0`, a sine equation on a conic); a circle
/// about the axis is an iso line, a line keeps its angle, and a curve in a plane
/// through the axis keeps its angle and its profile parameter is monotone. The rest
/// (other curve kinds; a torus's tube angle and a revolution's profile parameter off
/// those cases; B-spline and other surfaces) are the sign changes of the parameter's
/// rate along the edge, bracketed on the edge's knot spans (`2 (degree + 1)`
/// brackets each) and solved by [`LocalSolve1d::solve_by_bisection`].
///
/// **Units and parameterisations** are the kernel's STEP ones
/// ([`Surface::evaluate_step_point`]): a plane's `O + u sx X + v sy Y`, a cone's or
/// cylinder's angle in degrees and height along the axis, a sphere's longitude and
/// latitude in degrees, a torus's angle about the axis and round the tube in degrees,
/// a surface of revolution's angle in degrees and its profile's STEP parameter, and
/// the natural parameters of the other kinds.
///
/// **Periodic directions** (the angle of every surface of revolution, a torus's tube
/// angle and the profile angle of a revolved conic) are unwound along each loop, and
/// each later loop is moved by whole periods to where the box over it and the loops
/// before it is smallest (into the period of the first). A direction that a loop goes
/// round, or that the boundary spans by a full turn, is exactly one period. A face
/// whose loops go round the angle once more one way than the other contains what lies
/// on their left: its box reaches it (a sphere's pole, a cone's apex, the end of a
/// revolution's profile on the axis; a torus's whole tube). At a point on the axis
/// the angle is undefined: the loop passes round it on the face's side, from the angle
/// it arrives along to the one it leaves along. A face without edge loops on a sphere
/// or a torus is the whole surface.
///
/// **The result's angles** are moved by whole periods to start in `[0, 360)`, and
/// then one period down if they end above 360, within the kernel's angular limits
/// `-360 <= min <= max <= 360` (as [`crate::curve::step_edge::edge_interval`] keeps an
/// arc): a band from 300 degrees over the seam to 60 is `-60..60`, a full turn from 90
/// degrees `-270..90`.
///
/// **The margin** is added on each side of every direction that is not a full turn,
/// never beyond a sphere's poles, a cone's apex or, for the other kinds, the
/// surface's maximum analytic domain ([`Surface::max_analytic_domain`]) in its
/// directions that are not periodic; an angle widened to a full turn or more is one
/// period.
///
/// # Errors
///
/// [`SolidError::InvalidInput`] for a `tol` or a margin that is negative or not
/// finite, an empty loop, an edge interval that is not finite or ends below its start,
/// an edge or vertex farther than `tol` from a plane, cone, cylinder, sphere or torus,
/// no edge loop on a surface that is not closed, a face whose loops wind round the
/// angle more than once more one way than the other, or that reaches a side its
/// surface does not close (a cylinder bounded by one circle).
/// [`SolidError::Unsupported`] for a loop through a point on the axis along the axis.
/// Those of [`Surface::step_inversion`], [`Surface::evaluate_step`],
/// [`Surface::max_analytic_domain`] and the curves' evaluations.
pub fn face_window(
    k: &Kernel,
    surface: &Surface,
    loops: &[Vec<BoundaryEdge<'_>>],
    vertices: &[Point3],
    margin: WindowMargin,
    tol: f64,
) -> Result<Extent2d> {
    let _ = (k, surface, loops, vertices, margin, tol);
    Err(SolidError::Unsupported {
        operation: OP,
        reason: "not yet",
    })
}

#[cfg(test)]
#[allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]
mod tests {
    use super::*;
    use crate::curve::curve::CurveKind;
    use crate::curve::ellipse::Ellipse;
    use crate::curve::line::Line;
    use crate::geom::axis2_placement::Axis2Placement;
    use crate::math::Vec3;
    use crate::nurbs::b_spline_curve::{BSplineCurve, BSplineCurveForm, CanonicalData, KnotType};
    use crate::surface::cone::Cone;
    use crate::surface::plane::Plane;
    use crate::surface::sphere::Sphere;
    use crate::surface::surf_of_revolution::SurfOfRevolution;
    use crate::surface::torus::{Torus, TorusOptions};

    const TOL: f64 = 1e-9;

    fn k() -> Kernel {
        Kernel::default()
    }

    fn p(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z)
    }

    fn rad(deg: f64) -> f64 {
        deg.to_radians()
    }

    /// `a` equals the closed form `b` to rounding.
    #[track_caller]
    fn assert_close(a: f64, b: f64) {
        assert!(
            (a - b).abs() <= 1e-12 * (1.0 + b.abs()),
            "{a} is not {b} (difference {})",
            a - b
        );
    }

    #[track_caller]
    fn assert_window(w: &Extent2d, u: [f64; 2], v: [f64; 2]) {
        assert_close(w.u_min(), u[0]);
        assert_close(w.u_max(), u[1]);
        assert_close(w.v_min(), v[0]);
        assert_close(w.v_max(), v[1]);
    }

    /// The circle of radius `r` about `center`, its 0 degree point along `x` and its
    /// 90 degree point along `y`.
    fn circle(k: &Kernel, center: Point3, x: Vec3, y: Vec3, r: f64) -> Curve {
        let frame = Axis2Placement::new(center, x, y);
        Curve::from(Ellipse::create_canonical_circle(k, &frame, r, None, None).unwrap())
    }

    /// The circle of radius `r` about Z at height `z`.
    fn level(k: &Kernel, z: f64, r: f64) -> Curve {
        circle(k, p(0.0, 0.0, z), Vec3::X, Vec3::Y, r)
    }

    /// The point at `deg` degrees on the circle of radius `r` about Z at height `z`.
    fn round(r: f64, deg: f64, z: f64) -> Point3 {
        p(r * rad(deg).cos(), r * rad(deg).sin(), z)
    }

    /// The segment from `a` to `b`, its STEP parameter 0 at `a` and 1 at `b`.
    fn segment(a: Point3, b: Point3) -> Curve {
        Curve::from(Line::create_line_segment(3, a, b, None).unwrap())
    }

    fn edge(curve: &Curve, lo: f64, hi: f64, forward: bool) -> BoundaryEdge<'_> {
        BoundaryEdge {
            curve,
            interval: Extent1d::from_min_max(lo, hi),
            forward,
        }
    }

    /// A segment used from its start to its end.
    fn along(curve: &Curve) -> BoundaryEdge<'_> {
        edge(curve, 0.0, 1.0, true)
    }

    fn window(k: &Kernel, s: &Surface, loops: &[Vec<BoundaryEdge<'_>>]) -> Extent2d {
        face_window(k, s, loops, &[], WindowMargin::default(), TOL).unwrap()
    }

    fn cylinder(k: &Kernel) -> Surface {
        Surface::from(Cone::create_canonical_cylinder(k, &Axis2Placement::IDENTITY, 1.0).unwrap())
    }

    /// The edges of a band of the cylinder r = 1 between z = 0 and z = 3 and the angles
    /// `a` and `b`.
    struct Band {
        bottom: Curve,
        top: Curve,
        up: Curve,
        down: Curve,
        a: f64,
        b: f64,
    }

    impl Band {
        fn new(k: &Kernel, a: f64, b: f64) -> Self {
            Self {
                bottom: level(k, 0.0, 1.0),
                top: level(k, 3.0, 1.0),
                up: segment(round(1.0, b, 0.0), round(1.0, b, 3.0)),
                down: segment(round(1.0, a, 3.0), round(1.0, a, 0.0)),
                a,
                b,
            }
        }

        /// Counterclockwise in (u, v): along the bottom, up, back along the top, down.
        fn outer(&self) -> Vec<BoundaryEdge<'_>> {
            vec![
                edge(&self.bottom, self.a, self.b, true),
                along(&self.up),
                edge(&self.top, self.a, self.b, false),
                along(&self.down),
            ]
        }
    }

    /// Cylinder r = 1, the band z in [0, 3] over a quarter: u = [0, 90], v = [0, 3].
    #[test]
    fn a_cylinder_band_over_a_quarter() {
        let k = k();
        let band = Band::new(&k, 0.0, 90.0);
        let w = window(&k, &cylinder(&k), &[band.outer()]);
        assert_window(&w, [0.0, 90.0], [0.0, 3.0]);
    }

    /// The band from 300 degrees over the seam to 60, its arcs as #2609 gives them
    /// (-60..60) or split at the seam (300..360 and 0..60): either way u = [-60, 60],
    /// the kernel's convention for an angle over the seam.
    #[test]
    fn an_arc_across_the_seam() {
        let k = k();
        let s = cylinder(&k);
        let band = Band::new(&k, -60.0, 60.0);
        assert_window(&window(&k, &s, &[band.outer()]), [-60.0, 60.0], [0.0, 3.0]);

        let split = vec![
            edge(&band.bottom, 300.0, 360.0, true),
            edge(&band.bottom, 0.0, 60.0, true),
            along(&band.up),
            edge(&band.top, 0.0, 60.0, false),
            edge(&band.top, 300.0, 360.0, false),
            along(&band.down),
        ];
        assert_window(&window(&k, &s, &[split]), [-60.0, 60.0], [0.0, 3.0]);
    }

    /// Two full circles bound the whole band: u is exactly one period. From 0 degrees
    /// it is 0..360; from 90 degrees (#2609's -270..90 for a full turn from there) it
    /// starts at 90, one period down: -270..90.
    #[test]
    fn a_full_circle_is_one_period() {
        let k = k();
        let s = cylinder(&k);
        let (bottom, top) = (level(&k, 0.0, 1.0), level(&k, 3.0, 1.0));
        let full = |lo: f64, hi: f64| {
            vec![
                vec![edge(&bottom, lo, hi, true)],
                vec![edge(&top, lo, hi, false)],
            ]
        };
        assert_window(&window(&k, &s, &full(0.0, 360.0)), [0.0, 360.0], [0.0, 3.0]);
        assert_window(
            &window(&k, &s, &full(-270.0, 90.0)),
            [-270.0, 90.0],
            [0.0, 3.0],
        );
    }

    /// A rectangle on a plane whose frame is moved and turned: its box exactly.
    #[test]
    fn a_rectangle_on_a_plane() {
        let k = k();
        let (o, x, y) = (p(1.0, -2.0, 0.5), Vec3::Y, Vec3::Z);
        let s = Surface::from(Plane::create_canonical(&Axis2Placement::new(o, x, y)).unwrap());
        let at = |u: f64, v: f64| o + u * x + v * y;
        let sides = [
            segment(at(1.0, 2.0), at(4.0, 2.0)),
            segment(at(4.0, 2.0), at(4.0, 5.0)),
            segment(at(4.0, 5.0), at(1.0, 5.0)),
            segment(at(1.0, 5.0), at(1.0, 2.0)),
        ];
        let l: Vec<_> = sides.iter().map(along).collect();
        assert_window(&window(&k, &s, &[l]), [1.0, 4.0], [2.0, 5.0]);
    }

    /// An ellipse (a = 3, b = 1) turned by 30 degrees in a plane tilted by 45 degrees,
    /// centred at (2, 1) in the plane's frame: its box is the centre plus or minus
    /// sqrt(a^2 cos^2 30 + b^2 sin^2 30) = sqrt(7) in u and
    /// sqrt(a^2 sin^2 30 + b^2 cos^2 30) = sqrt(3) in v. On the arc -45..135 closed by
    /// its chord, the u and v maxima are inside the arc (at -10.9 and 30 degrees) and
    /// the u minimum is at an end.
    #[test]
    fn a_turned_ellipse_on_a_tilted_plane() {
        let k = k();
        let h = 0.5f64.sqrt();
        let (x, y) = (Vec3::new(h, 0.0, h), Vec3::Y);
        let s =
            Surface::from(Plane::create_canonical(&Axis2Placement::new(Vec3::ZERO, x, y)).unwrap());
        let (c30, s30) = (rad(30.0).cos(), rad(30.0).sin());
        let (ex, ey) = (c30 * x + s30 * y, -s30 * x + c30 * y);
        let centre = 2.0 * x + y;
        let frame = Axis2Placement::new(centre, ex, ey);
        let e = Curve::from(Ellipse::create_canonical(&k, &frame, 3.0, 1.0, None, None).unwrap());
        let (r7, r3) = (7f64.sqrt(), 3f64.sqrt());
        let w = window(&k, &s, &[vec![edge(&e, 0.0, 360.0, true)]]);
        assert_window(&w, [2.0 - r7, 2.0 + r7], [1.0 - r3, 1.0 + r3]);

        let at = |t: f64| centre + 3.0 * rad(t).cos() * ex + rad(t).sin() * ey;
        let chord = segment(at(135.0), at(-45.0));
        let w = window(&k, &s, &[vec![edge(&e, -45.0, 135.0, true), along(&chord)]]);
        let u = |t: f64| at(t).dot(x);
        let v = |t: f64| at(t).dot(y);
        assert_window(
            &w,
            [u(135.0).min(u(-45.0)), 2.0 + r7],
            [v(135.0).min(v(-45.0)), 1.0 + r3],
        );
    }

    fn sphere(k: &Kernel) -> Surface {
        Surface::from(Sphere::create_canonical(k, &Axis2Placement::IDENTITY, 2.0).unwrap())
    }

    /// Sphere r = 2: the cap above latitude 30, bounded by its rim alone (running up u,
    /// the pole on its left), reaches the pole: u one period, v = [30, 90]. The quarter
    /// of it between the meridians at 0 and 90 degrees, up one to the pole and down the
    /// other: u = [0, 90], v = [30, 90]. The cap below latitude -30, its rim running
    /// down u: v = [-90, -30].
    #[test]
    fn a_spherical_cap_to_the_pole() {
        let k = k();
        let s = sphere(&k);
        let (z, r) = (2.0 * rad(30.0).sin(), 2.0 * rad(30.0).cos());
        let rim = level(&k, z, r);
        let w = window(&k, &s, &[vec![edge(&rim, 0.0, 360.0, true)]]);
        assert_window(&w, [0.0, 360.0], [30.0, 90.0]);

        let east = circle(&k, Vec3::ZERO, Vec3::Y, Vec3::Z, 2.0);
        let west = circle(&k, Vec3::ZERO, Vec3::X, Vec3::Z, 2.0);
        let quarter = vec![
            edge(&rim, 0.0, 90.0, true),
            edge(&east, 30.0, 90.0, true),
            edge(&west, 30.0, 90.0, false),
        ];
        assert_window(&window(&k, &s, &[quarter]), [0.0, 90.0], [30.0, 90.0]);

        let south = level(&k, -z, r);
        let w = window(&k, &s, &[vec![edge(&south, 0.0, 360.0, false)]]);
        assert_window(&w, [0.0, 360.0], [-90.0, -30.0]);
    }

    /// A great circle of the sphere tilted by 20 degrees about Y: its latitude reaches
    /// +-20 inside the edge (at 0 and 180 degrees on it) and its longitude goes round
    /// once, so the face on its left reaches the north pole: v = [-20, 90].
    #[test]
    fn a_tilted_great_circle_on_a_sphere() {
        let k = k();
        let tilt = Vec3::new(rad(20.0).cos(), 0.0, rad(20.0).sin());
        let g = circle(&k, Vec3::ZERO, tilt, Vec3::Y, 2.0);
        let w = window(&k, &sphere(&k), &[vec![edge(&g, 0.0, 360.0, true)]]);
        assert_window(&w, [0.0, 360.0], [-20.0, 90.0]);
    }

    fn frustum(k: &Kernel) -> Surface {
        Surface::from(
            Cone::create_canonical(k, &Axis2Placement::IDENTITY, 2.0, 1.0, 2.0, false, false)
                .unwrap(),
        )
    }

    /// The cone from radius 2 at z = 0 to radius 1 at z = 2 (apex at z = 4): the band
    /// between its rims is u one period, v = [0, 2]; the part between the generators
    /// at 0 and 120 degrees is u = [0, 120]; bounded by its lower rim alone (running up
    /// u), the face reaches the apex: v = [0, 4].
    #[test]
    fn a_cone_frustum() {
        let k = k();
        let s = frustum(&k);
        let (bottom, top) = (level(&k, 0.0, 2.0), level(&k, 2.0, 1.0));
        let w = window(
            &k,
            &s,
            &[
                vec![edge(&bottom, 0.0, 360.0, true)],
                vec![edge(&top, 0.0, 360.0, false)],
            ],
        );
        assert_window(&w, [0.0, 360.0], [0.0, 2.0]);

        let up = segment(round(2.0, 120.0, 0.0), round(1.0, 120.0, 2.0));
        let down = segment(round(1.0, 0.0, 2.0), round(2.0, 0.0, 0.0));
        let part = vec![
            edge(&bottom, 0.0, 120.0, true),
            along(&up),
            edge(&top, 0.0, 120.0, false),
            along(&down),
        ];
        assert_window(&window(&k, &s, &[part]), [0.0, 120.0], [0.0, 2.0]);

        let w = window(&k, &s, &[vec![edge(&bottom, 0.0, 360.0, true)]]);
        assert_window(&w, [0.0, 360.0], [0.0, 4.0]);
    }

    fn torus(k: &Kernel) -> Surface {
        let opts = TorusOptions::default();
        Surface::from(
            Torus::create_canonical(k, &Axis2Placement::IDENTITY, 3.0, 1.0, opts, None).unwrap(),
        )
    }

    /// The tube circle of the torus (R = 3, r = 1) at `deg` degrees about Z, its angle
    /// the torus's v.
    fn tube(k: &Kernel, deg: f64) -> Curve {
        let radial = Vec3::new(rad(deg).cos(), rad(deg).sin(), 0.0);
        circle(k, 3.0 * radial, radial, Vec3::Z, 1.0)
    }

    /// Torus R = 3, r = 1: the quarter between 0 and 90 degrees over the upper half of
    /// the tube (outer equator over the top to the inner one): u = [0, 90],
    /// v = [0, 180]; the outer half of the tube (v from -90 over the seam to 90):
    /// v = [-90, 90].
    #[test]
    fn a_torus_section() {
        let k = k();
        let s = torus(&k);
        let (outer, inner) = (level(&k, 0.0, 4.0), level(&k, 0.0, 2.0));
        let (t0, t90) = (tube(&k, 0.0), tube(&k, 90.0));
        let upper = vec![
            edge(&outer, 0.0, 90.0, true),
            edge(&t90, 0.0, 180.0, true),
            edge(&inner, 0.0, 90.0, false),
            edge(&t0, 0.0, 180.0, false),
        ];
        assert_window(&window(&k, &s, &[upper]), [0.0, 90.0], [0.0, 180.0]);

        let (low, high) = (level(&k, -1.0, 3.0), level(&k, 1.0, 3.0));
        let outside = vec![
            edge(&low, 0.0, 90.0, true),
            edge(&t90, -90.0, 90.0, true),
            edge(&high, 0.0, 90.0, false),
            edge(&t0, -90.0, 90.0, false),
        ];
        assert_window(&window(&k, &s, &[outside]), [0.0, 90.0], [-90.0, 90.0]);
    }

    /// The profile segment from (1, 0, 0) to (2, 0, 1) turned about Z: the band between
    /// its end circles over a quarter is u = [0, 90] and v the profile's [0, 1].
    #[test]
    fn a_revolution_of_a_line() {
        let k = k();
        let profile = segment(p(1.0, 0.0, 0.0), p(2.0, 0.0, 1.0));
        let s = Surface::from(
            SurfOfRevolution::create_canonical(&k, profile, Vec3::ZERO, Vec3::Z, false).unwrap(),
        );
        let (bottom, top) = (level(&k, 0.0, 1.0), level(&k, 1.0, 2.0));
        let up = segment(round(1.0, 90.0, 0.0), round(2.0, 90.0, 1.0));
        let down = segment(round(2.0, 0.0, 1.0), round(1.0, 0.0, 0.0));
        let l = vec![
            edge(&bottom, 0.0, 90.0, true),
            along(&up),
            edge(&top, 0.0, 90.0, false),
            along(&down),
        ];
        assert_window(&window(&k, &s, &[l]), [0.0, 90.0], [0.0, 1.0]);
    }

    /// A band of the cylinder from 150 to 270 degrees (across the angle's branch at 180)
    /// with a hole between 200 and 220 at z in [1, 2]: the hole's angles come out
    /// -160..-140 on their own and are moved into the outer loop's period, so the
    /// window stays u = [150, 270]. The hole alone is 200..220.
    #[test]
    fn a_hole_in_the_period_of_the_outer_loop() {
        let k = k();
        let s = cylinder(&k);
        let band = Band::new(&k, 150.0, 270.0);
        let (lower, upper) = (level(&k, 1.0, 1.0), level(&k, 2.0, 1.0));
        let left = segment(round(1.0, 200.0, 1.0), round(1.0, 200.0, 2.0));
        let right = segment(round(1.0, 220.0, 2.0), round(1.0, 220.0, 1.0));
        // clockwise in (u, v): the face on its left is outside it
        let hole = vec![
            along(&left),
            edge(&upper, 200.0, 220.0, true),
            along(&right),
            edge(&lower, 200.0, 220.0, false),
        ];
        let w = window(&k, &s, &[band.outer(), hole.clone()]);
        assert_window(&w, [150.0, 270.0], [0.0, 3.0]);
        assert_window(&window(&k, &s, &[hole]), [200.0, 220.0], [1.0, 2.0]);
    }

    /// A margin widens every direction that is not a full turn by exactly itself, but
    /// never past a pole.
    #[test]
    fn a_margin_widens_by_exactly_itself() {
        let k = k();
        let band = Band::new(&k, 0.0, 90.0);
        let m = WindowMargin { u: 5.0, v: 0.5 };
        let w = face_window(&k, &cylinder(&k), &[band.outer()], &[], m, TOL).unwrap();
        assert_window(&w, [-5.0, 95.0], [-0.5, 3.5]);

        let rim = level(&k, 2.0 * rad(30.0).sin(), 2.0 * rad(30.0).cos());
        let m = WindowMargin { u: 5.0, v: 10.0 };
        let cap = [vec![edge(&rim, 0.0, 360.0, true)]];
        let w = face_window(&k, &sphere(&k), &cap, &[], m, TOL).unwrap();
        assert_window(&w, [0.0, 360.0], [20.0, 90.0]);
    }

    /// A vertex loop bounds the directions that are not periodic: a cone's apex.
    #[test]
    fn a_vertex_loop_at_the_apex() {
        let k = k();
        let bottom = level(&k, 0.0, 2.0);
        let rim = [vec![edge(&bottom, 0.0, 360.0, true)]];
        let apex = [p(0.0, 0.0, 4.0)];
        let w = face_window(&k, &frustum(&k), &rim, &apex, WindowMargin::default(), TOL).unwrap();
        assert_window(&w, [0.0, 360.0], [0.0, 4.0]);
    }

    /// Faces without edge loops: the whole sphere and torus; a plane has no whole.
    #[test]
    fn faces_without_edge_loops() {
        let k = k();
        assert_window(&window(&k, &sphere(&k), &[]), [0.0, 360.0], [-90.0, 90.0]);
        assert_window(&window(&k, &torus(&k), &[]), [0.0, 360.0], [0.0, 360.0]);
        let plane = Surface::from(Plane::create_canonical(&Axis2Placement::IDENTITY).unwrap());
        assert!(matches!(
            face_window(&k, &plane, &[], &[], WindowMargin::default(), TOL),
            Err(SolidError::InvalidInput { .. })
        ));
    }

    /// The quadratic Bezier from (0, 0) over (1, 2) to (2, 0) on the XY plane, closed
    /// by its chord: y = 4 t (1 - t) peaks at 1 inside the edge, found as the root of
    /// its rate (a spline has no closed form here): u = [0, 2], v = [0, 1].
    #[test]
    fn a_spline_arc_on_a_plane() {
        let k = k();
        let s = Surface::from(Plane::create_canonical(&Axis2Placement::IDENTITY).unwrap());
        let arc = Curve::from(
            BSplineCurve::create_canonical(&CanonicalData {
                dim: 3,
                degree: 2,
                control_points: &[p(0.0, 0.0, 0.0), p(1.0, 2.0, 0.0), p(2.0, 0.0, 0.0)],
                form: BSplineCurveForm::PolylineForm,
                knot_multiplicities: &[3, 3],
                knots: &[0.0, 1.0],
                knot_type: KnotType::Unspecified,
                weights: None,
                trim: None,
            })
            .unwrap(),
        );
        let chord = segment(p(2.0, 0.0, 0.0), p(0.0, 0.0, 0.0));
        let w = window(&k, &s, &[vec![edge(&arc, 0.0, 1.0, true), along(&chord)]]);
        assert_window(&w, [0.0, 2.0], [0.0, 1.0]);
    }

    /// `circle` as its rational B-spline (the kernel's NURB of it), over its natural
    /// interval.
    fn as_spline(c: &Curve) -> (Curve, Extent1d) {
        let CurveKind::Ellipse(e) = c.kind() else {
            panic!("not a conic");
        };
        let b = Curve::from(e.base().clone());
        let ivl = b.natural_interval();
        (b, ivl)
    }

    /// The cylinder band, a torus section and the tilted great circle of the sphere
    /// with their circles as B-splines: the angles unwind along the splines' natural
    /// parameters, and the great circle's latitude extremes (+-20) are the roots of its
    /// rate.
    #[test]
    fn spline_circles_on_a_cylinder_and_a_sphere() {
        let k = k();
        let (bottom, bi) = as_spline(&level(&k, 0.0, 1.0));
        let (top, ti) = as_spline(&level(&k, 3.0, 1.0));
        let spline = |c, i: Extent1d, forward| edge(c, i.min(), i.max(), forward);
        let band = [
            vec![spline(&bottom, bi, true)],
            vec![spline(&top, ti, false)],
        ];
        assert_window(&window(&k, &cylinder(&k), &band), [0.0, 360.0], [0.0, 3.0]);

        // the upper quarter of `a_torus_section`, its tube circles' halves 0..180 as
        // the splines' 0..0.5 (their knots at the quarter turns)
        let (outer, inner) = (level(&k, 0.0, 4.0), level(&k, 0.0, 2.0));
        let ((t0, _), (t90, _)) = (as_spline(&tube(&k, 0.0)), as_spline(&tube(&k, 90.0)));
        let upper = vec![
            edge(&outer, 0.0, 90.0, true),
            edge(&t90, 0.0, 0.5, true),
            edge(&inner, 0.0, 90.0, false),
            edge(&t0, 0.0, 0.5, false),
        ];
        assert_window(&window(&k, &torus(&k), &[upper]), [0.0, 90.0], [0.0, 180.0]);

        let tilt = Vec3::new(rad(20.0).cos(), 0.0, rad(20.0).sin());
        let (g, gi) = as_spline(&circle(&k, Vec3::ZERO, tilt, Vec3::Y, 2.0));
        let w = window(&k, &sphere(&k), &[vec![spline(&g, gi, true)]]);
        assert_window(&w, [0.0, 360.0], [-20.0, 90.0]);
    }

    #[test]
    fn invalid_inputs_are_refused() {
        let k = k();
        let s = cylinder(&k);
        let (bottom, top) = (level(&k, 0.0, 1.0), level(&k, 3.0, 1.0));
        let band = vec![
            vec![edge(&bottom, 0.0, 360.0, true)],
            vec![edge(&top, 0.0, 360.0, false)],
        ];
        let refused = |loops: &[Vec<BoundaryEdge<'_>>], m: WindowMargin, tol: f64| {
            matches!(
                face_window(&k, &s, loops, &[], m, tol),
                Err(SolidError::InvalidInput { .. })
            )
        };
        let none = WindowMargin::default();
        assert!(refused(&band, none, -1.0));
        assert!(refused(&band, none, f64::NAN));
        assert!(refused(&band, WindowMargin { u: -1.0, v: 0.0 }, TOL));
        assert!(refused(
            &band,
            WindowMargin {
                u: 0.0,
                v: f64::INFINITY
            },
            TOL
        ));
        assert!(refused(&[vec![]], none, TOL));
        let not_finite = BoundaryEdge {
            curve: &bottom,
            interval: Extent1d::from_raw(0.0, f64::NAN),
            forward: true,
        };
        assert!(refused(&[vec![not_finite]], none, TOL));
        let reversed = BoundaryEdge {
            curve: &bottom,
            interval: Extent1d::from_raw(90.0, 0.0),
            forward: true,
        };
        assert!(refused(&[vec![reversed]], none, TOL));
        // a circle off the cylinder
        let off = level(&k, 0.0, 1.5);
        assert!(refused(&[vec![edge(&off, 0.0, 360.0, true)]], none, TOL));
        // one circle does not bound a cylinder
        assert!(refused(&[vec![edge(&bottom, 0.0, 360.0, true)]], none, TOL));
        // both circles round the same way
        let twice = vec![
            vec![edge(&bottom, 0.0, 360.0, true)],
            vec![edge(&top, 0.0, 360.0, true)],
        ];
        assert!(refused(&twice, none, TOL));
    }
}
