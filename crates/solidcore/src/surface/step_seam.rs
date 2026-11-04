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
//! Deviations: not a port. See [`seam_into_widest_gap`] for the rules.

use crate::error::{Result, SolidError};
use crate::geom::axis2_placement::Axis2Placement;
use crate::geom::extent1d::Extent1d;
use crate::kernel::Kernel;
use crate::math::Point3;
use crate::surface::step_window::BoundaryEdge;
use crate::surface::surface::Surface;

const OP: &str = "surface::step_seam::seam_into_widest_gap";

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
/// # Errors
///
/// [`SolidError::InvalidInput`] for invalid input; [`SolidError::Unsupported`] for
/// edges whose angles are not computed.
pub fn seam_into_widest_gap(
    k: &Kernel,
    surface: &Surface,
    loops: &[Vec<BoundaryEdge<'_>>],
    seam_point: Option<Point3>,
    tol: f64,
) -> Result<PlacedSeam> {
    let _ = (k, surface, loops, seam_point, tol, Axis2Placement::IDENTITY);
    Err(SolidError::invalid_input(OP, "not implemented yet"))
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
        assert_eq!(position(&placed.surface).x_axis(), v3(0.0, -1.0, 0.0));
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
