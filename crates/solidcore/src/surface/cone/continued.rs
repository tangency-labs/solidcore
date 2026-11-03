//! A cone or cylinder continued naturally past one side of a face's parameter box: the
//! analytic itself over a wider STEP domain, so a sheet extended past a side keeps its kind.
//!
//! Not upstream: SolidCraft's composition (`cone_continued`), moved into the kernel so the
//! modelling layer computes no geometry (docs/layers.md §1). It is an exact class at bound 0
//! (owner rule 13564): the result is a [`Cone`] with the source's placement, radius and
//! semi-angle (its label), and its NURBS form is proved on that one implicit surface
//! ([`AnalyticSupport::check`]); a continuation is refused unless both hold.
//!
//! Deviations: none (no upstream counterpart). No recursion.

use crate::error::{Result, SolidError};
use crate::geom::extent2d::Extent2d;
use crate::kernel::Kernel;
use crate::surface::cone::Cone;
use crate::surface::implicit::AnalyticSupport;
use crate::surface::side::BoxSide;
use crate::surface::surface::Surface;

const OP: &str = "surface::cone::Cone::continued";

/// A full turn, degrees.
const TURN: f64 = 360.0;

impl Cone {
    /// The cone continued naturally `distance` past side `side` of the face box `dom` (in the
    /// cone's NURBS parameters): its STEP domain widened at that side, by the angle
    /// `distance / r` (degrees) where the side runs along the axis (`r` the radius at the
    /// side's middle), or by the height `distance cos(semi-angle)` where it runs round it
    /// (so the slant length added is `distance`). The rest of the domain is `dom`'s.
    ///
    /// Which STEP coordinate the side holds constant follows from the parametrisation (the
    /// NURBS u is the height when the cone's uv is swapped), not from a measurement. A
    /// sweep that comes out over a turn by no more than the rounding of the sum is a full
    /// turn.
    ///
    /// # Errors
    ///
    /// `InvalidInput` for a distance not finite and positive, a side along the axis at a
    /// radius not positive, a sweep proved over a turn, a height past the apex (from
    /// [`Cone::adjust_step_uv_domain`]), or a result whose label or NURBS form is not the
    /// source's surface; those of the parameter conversions and of
    /// [`AnalyticSupport::check`].
    pub fn continued(
        &self,
        k: &Kernel,
        dom: &Extent2d,
        side: BoxSide,
        distance: f64,
    ) -> Result<Self> {
        if !(distance.is_finite() && distance > 0.0) {
            return Err(SolidError::invalid_input(
                OP,
                "the distance is not finite and positive",
            ));
        }
        let step = |uv| self.convert_uv_from_nurbs_to_step(k, uv);
        let (lo, hi, mid) = (
            step(dom.min())?,
            step(dom.max())?,
            step(side.point(dom, 0.5))?,
        );
        let mut angle = [lo.x.min(hi.x), lo.x.max(hi.x)];
        let mut height = [lo.y.min(hi.y), lo.y.max(hi.y)];
        // the side holds the NURBS u constant for a u side; that is the height when swapped
        let in_height = side.is_u() == self.canonical().swap_uv;
        let (row, at, by) = if in_height {
            let by = distance * self.semi_angle_deg().to_radians().cos();
            (&mut height, mid.y, by)
        } else {
            let r = self.radius(mid.y);
            if r.is_nan() || r <= 0.0 {
                return Err(SolidError::invalid_input(
                    OP,
                    "the side runs along the axis at a radius not positive",
                ));
            }
            (&mut angle, mid.x, (distance / r).to_degrees())
        };
        // the side is the row's end nearer its middle's coordinate
        let high = (at - row[1]).abs() <= (at - row[0]).abs();
        if high {
            row[1] += by;
        } else {
            row[0] -= by;
        }
        if !in_height {
            // the sweep's sum rounds by at most 2 ulps of its largest term
            let span = angle[1] - angle[0];
            let rounding = 2.0 * f64::EPSILON * angle[0].abs().max(angle[1].abs()).max(TURN);
            if span - TURN > rounding {
                return Err(SolidError::invalid_input(OP, "the sweep is over a turn"));
            }
            if span > TURN {
                // a full turn, to rounding: the end grown is the one moved
                if high {
                    angle[1] = angle[0] + TURN;
                } else {
                    angle[0] = angle[1] - TURN;
                }
            }
        }
        let wider = Extent2d::from_bounds(angle[0], height[0], angle[1], height[1]);
        let mut cone = self.clone();
        cone.adjust_step_uv_domain(k, &wider)?;
        // bound 0: the source's label, and the result's form on its implicit surface
        let was = AnalyticSupport::of(&Surface::from(self.clone()));
        let now = AnalyticSupport::check(k, &Surface::from(cone.clone()))?;
        if was != Some(now) {
            return Err(SolidError::invalid_input(
                OP,
                "the continuation is not the source's surface",
            ));
        }
        Ok(cone)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::axis2_placement::Axis2Placement;
    use crate::geom::vector2d::Point2;
    use crate::geom::vector3d::Vec3;

    fn k() -> Kernel {
        Kernel::default()
    }

    /// A cone over the STEP box `[a0, a1] x [h0, h1]`, with its NURBS box.
    fn cone(bot: f64, top: f64, h: f64, a: [f64; 2], hs: [f64; 2], swap: bool) -> (Cone, Extent2d) {
        let k = k();
        let pl = Axis2Placement::new(Vec3::new(1.0, -2.0, 0.5), Vec3::X, Vec3::Y);
        let mut c = Cone::create_canonical(&k, &pl, bot, top, h, swap, false).unwrap();
        c.adjust_step_uv_domain(&k, &Extent2d::from_bounds(a[0], hs[0], a[1], hs[1]))
            .unwrap();
        let dom = Surface::from(c.clone())
            .natural_uv_domain_with_nurb(&k)
            .unwrap();
        (c, dom)
    }

    fn step_box(c: &Cone) -> Extent2d {
        Surface::from(c.clone())
            .step_uv_domain_with_nurb(&k())
            .unwrap()
    }

    /// The side of `dom` that holds the STEP angle (`in_height` false) or height constant
    /// at its high or low end.
    fn side_of(c: &Cone, dom: &Extent2d, in_height: bool, high: bool) -> BoxSide {
        let k = k();
        BoxSide::ALL
            .into_iter()
            .find(|s| {
                let a = c
                    .convert_uv_from_nurbs_to_step(&k, s.point(dom, 0.0))
                    .unwrap();
                let z = c
                    .convert_uv_from_nurbs_to_step(&k, s.point(dom, 1.0))
                    .unwrap();
                let held = if in_height { (a.y, z.y) } else { (a.x, z.x) };
                let sb = step_box(c);
                let end = if in_height {
                    if high {
                        sb.max().y
                    } else {
                        sb.min().y
                    }
                } else if high {
                    sb.max().x
                } else {
                    sb.min().x
                };
                (held.0 - end).abs() < 1e-9 && (held.1 - end).abs() < 1e-9
            })
            .unwrap()
    }

    #[test]
    fn a_cylinder_continued_round_its_axis_by_the_arc() {
        // radius 2 over a quarter turn: pi more of arc is a half turn
        for swap in [false, true] {
            let (c, dom) = cone(2.0, 2.0, 3.0, [0.0, 90.0], [0.0, 3.0], swap);
            let s = side_of(&c, &dom, false, true);
            let d = c.continued(&k(), &dom, s, std::f64::consts::PI).unwrap();
            let b = step_box(&d);
            assert!((b.min().x - 0.0).abs() < 1e-12 && (b.max().x - 180.0).abs() < 1e-12);
            assert!((b.min().y - 0.0).abs() < 1e-12 && (b.max().y - 3.0).abs() < 1e-12);
            assert!(d.is_cylinder());
        }
    }

    #[test]
    fn a_cone_continued_along_its_axis_by_the_slant_length() {
        // radius 1 to 2 over height 3: tan(semi) = 1/3
        let (c, dom) = cone(1.0, 2.0, 3.0, [0.0, 90.0], [0.0, 3.0], false);
        let s = side_of(&c, &dom, true, true);
        let d = c.continued(&k(), &dom, s, 2.0).unwrap();
        let b = step_box(&d);
        let dh = 2.0 * 3.0 / 10f64.sqrt();
        assert!((b.max().y - (3.0 + dh)).abs() < 1e-12);
        // the slant from the old top to the new along a generator is the distance
        let at = |h: f64| d.evaluate_step_point(Point2::new(45.0, h)).unwrap();
        assert!((at(3.0 + dh).distance(at(3.0)) - 2.0).abs() < 1e-12);
        // the radius grows as the source's: r(h) = 1 + h / 3
        assert!((d.radius(3.0 + dh) - (1.0 + (3.0 + dh) / 3.0)).abs() < 1e-12);
    }

    #[test]
    fn the_low_end_grows_down() {
        let (c, dom) = cone(1.0, 2.0, 3.0, [10.0, 100.0], [1.0, 3.0], true);
        let s = side_of(&c, &dom, true, false);
        let d = c.continued(&k(), &dom, s, 0.5).unwrap();
        let b = step_box(&d);
        assert!((b.min().y - (1.0 - 0.5 * 3.0 / 10f64.sqrt())).abs() < 1e-12);
        assert!((b.max().y - 3.0).abs() < 1e-12);
        let s = side_of(&c, &dom, false, false);
        let d = c.continued(&k(), &dom, s, 0.5).unwrap();
        let b = step_box(&d);
        // r at the side's middle (height 2) is 5/3
        assert!((b.min().x - (10.0 - (0.5 * 3.0 / 5.0f64).to_degrees())).abs() < 1e-12);
    }

    #[test]
    fn past_the_apex_or_over_a_turn_is_refused() {
        let (c, dom) = cone(1.0, 2.0, 3.0, [0.0, 90.0], [0.0, 3.0], false);
        // the apex is 3 below the origin: a slant of 4 reaches past it
        let s = side_of(&c, &dom, true, false);
        assert!(c.continued(&k(), &dom, s, 4.0).is_err());
        let (c, dom) = cone(1.0, 1.0, 3.0, [0.0, 270.0], [0.0, 3.0], false);
        let s = side_of(&c, &dom, false, true);
        assert!(c.continued(&k(), &dom, s, 2.0).is_err());
        // exactly to a full turn is one
        let d = c
            .continued(&k(), &dom, s, std::f64::consts::FRAC_PI_2)
            .unwrap();
        let b = step_box(&d);
        assert_eq!(b.max().x - b.min().x, TURN);
        assert!(c.continued(&k(), &dom, s, 0.0).is_err());
        assert!(c.continued(&k(), &dom, s, f64::NAN).is_err());
    }

    #[test]
    fn a_moved_form_is_refused() {
        // the NURBS form moved off the label: no longer an exact class
        let (mut c, dom) = cone(1.0, 2.0, 3.0, [0.0, 90.0], [0.0, 3.0], false);
        let s = side_of(&c, &dom, true, true);
        let d = c.continued(&k(), &dom, s, 1.0).unwrap();
        let mut moved = d.clone();
        let shift = Axis2Placement::new(Vec3::new(0.0, 0.0, 1e-6), Vec3::X, Vec3::Y);
        moved.rev.base.transform(&shift, None).unwrap();
        assert!(AnalyticSupport::check(&k(), &Surface::from(moved)).is_err());
        // a source whose label differs from its form gives no continuation
        c.base_radius += 1e-6;
        assert!(c.continued(&k(), &dom, s, 1.0).is_err());
    }
}
