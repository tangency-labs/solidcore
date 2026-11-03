//! An analytic surface continued naturally past one side of a face's parameter box: the
//! same plane, cone, sphere or torus over a wider STEP domain, so a sheet extended past a side
//! keeps its kind and is never re-expressed as a B-spline.
//!
//! Not upstream: SolidCraft's composition (`cone_continued`, generalised), moved into the
//! kernel so the modelling layer computes no geometry (docs/layers.md §1). Each result is an
//! exact class at bound 0 (owner rule 13564): the source's label (placement, radii, angle)
//! unchanged, and the result's NURBS form proved on that one implicit surface
//! ([`AnalyticSupport::check`]). A continuation is refused unless both hold.
//!
//! Deviations: none (no upstream counterpart). No recursion.

use crate::error::{Result, SolidError};
use crate::geom::extent2d::Extent2d;
use crate::geom::vector2d::Point2;
use crate::geom::vector3d::{Point3, Vec3};
use crate::kernel::Kernel;
use crate::surface::implicit::AnalyticSupport;
use crate::surface::side::BoxSide;
use crate::surface::surface::{Surface, SurfaceKind};

const OP: &str = "surface::continued::Surface::continued";

/// A full turn, degrees.
const TURN: f64 = 360.0;

/// Widens `row` at its end nearer `at` by `by`; whether that was its high end.
fn widen(row: &mut [f64; 2], at: f64, by: f64) -> bool {
    let high = (at - row[1]).abs() <= (at - row[0]).abs();
    if high {
        row[1] += by;
    } else {
        row[0] -= by;
    }
    high
}

/// Holds `row` within `[lo, hi]`: refused when proved past them (beyond the rounding of the
/// sums that made it, 2 ulps of its largest term), set to them when within that rounding.
fn within(row: &mut [f64; 2], lo: f64, hi: f64, what: &'static str) -> Result<()> {
    let rounding = 2.0 * f64::EPSILON * row[0].abs().max(row[1].abs()).max(hi.abs()).max(lo.abs());
    if row[0] < lo - rounding || row[1] > hi + rounding {
        return Err(SolidError::invalid_input(OP, what));
    }
    row[0] = row[0].max(lo);
    row[1] = row[1].min(hi);
    Ok(())
}

/// Holds a sweep `row` grown at its `high` end (or low) to at most a turn: refused when proved
/// over one, set to exactly one when within the rounding of the sum.
fn at_most_a_turn(row: &mut [f64; 2], high: bool, what: &'static str) -> Result<()> {
    let span = row[1] - row[0];
    let rounding = 2.0 * f64::EPSILON * row[0].abs().max(row[1].abs()).max(TURN);
    if span - TURN > rounding {
        return Err(SolidError::invalid_input(OP, what));
    }
    if span > TURN {
        if high {
            row[1] = row[0] + TURN;
        } else {
            row[0] = row[1] - TURN;
        }
    }
    Ok(())
}

impl Surface {
    /// The surface continued naturally `distance` past side `side` of the face box `dom` (in
    /// its NURBS parameters), as the same analytic surface over a wider STEP domain: the
    /// length added across the side, measured at the side's middle, is `distance`.
    ///
    /// - **Plane:** `distance / s` more of the parameter the side holds, `s` its scale.
    /// - **Cone, cylinder:** [`crate::surface::cone::Cone::continued`].
    /// - **Sphere, torus:** round the axis, the angle `distance / rho` (degrees), `rho` the
    ///   distance from the axis of the side's middle; along the generator, `distance / r`, `r`
    ///   its circle's radius. A sphere's latitude stays within the poles, every angle within a
    ///   turn (both to the rounding of the sum, then held to them exactly).
    ///
    /// Which STEP coordinate the side holds follows from the parametrisation (a revolution's
    /// NURBS u is its sweep unless swapped; a plane's NURBS uv are its STEP uv), not from a
    /// measurement.
    ///
    /// # Errors
    ///
    /// `InvalidInput` for a source not proved an exact class ([`AnalyticSupport::check`]:
    /// a B-spline's extension is not this, and a form moved off its label is refused);
    /// `Unsupported` for an analytic kind with no continuation here; `InvalidInput` for a distance not finite and positive, a side at the axis, a
    /// domain proved past a pole or over a turn, or a result not the source's surface (label
    /// and form); those of the conversions, the domain adjustments and
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
        // the source an exact class first: a form moved off its label is refused, not remade
        AnalyticSupport::check(k, self).map_err(|e| match e {
            SolidError::InvalidInput { reason, .. } => SolidError::invalid_input(OP, reason),
            other => other,
        })?;
        let grown = match self.kind() {
            SurfaceKind::Cone(c) => Self::from(c.continued(k, dom, side, distance)?),
            SurfaceKind::Plane(p) => {
                let mut p = p.clone();
                let (s, d) = (p.uv_scale(), p.step_uv_domain());
                let mid = side.point(dom, 0.5);
                let mut u = [d.min().x, d.max().x];
                let mut v = [d.min().y, d.max().y];
                if side.is_u() {
                    widen(&mut u, mid.x, distance / s.x);
                } else {
                    widen(&mut v, mid.y, distance / s.y);
                }
                p.adjust_step_uv_domain(&Extent2d::from_bounds(u[0], v[0], u[1], v[1]))?;
                Self::from(p)
            }
            SurfaceKind::Sphere(sp) => {
                let mut sp = sp.as_ref().clone();
                let (o, z) = (sp.base.position.origin(), sp.base.position.z_axis());
                let swap = sp.base.swap_uv;
                let radius = sp.radius();
                let step = |uv| sp.convert_uv_from_nurbs_to_step(k, uv);
                let rows = revolution_rows(&step, dom, side, swap)?;
                let along = side.is_u() == swap;
                let rho = if along {
                    radius
                } else {
                    axis_distance(sp.evaluate_step_point(rows.mid)?, o, z)
                };
                let (mut a, mut b) = (rows.sweep, rows.generator);
                let row = if along { &mut b } else { &mut a };
                let high = widen(row, rows.at, rate(rho, distance)?);
                if along {
                    within(&mut b, -90.0, 90.0, "the latitude is past a pole")?;
                } else {
                    at_most_a_turn(&mut a, high, "the sweep is over a turn")?;
                }
                sp.adjust_step_uv_domain(k, &Extent2d::from_bounds(a[0], b[0], a[1], b[1]))?;
                Self::from(sp)
            }
            SurfaceKind::Torus(t) => {
                let mut t = t.as_ref().clone();
                let (o, z) = (t.base.position.origin(), t.base.position.z_axis());
                let swap = t.base.swap_uv;
                let minor = t.minor_radius();
                let step = |uv| t.convert_uv_from_nurbs_to_step(k, uv);
                let rows = revolution_rows(&step, dom, side, swap)?;
                let along = side.is_u() == swap;
                let rho = if along {
                    minor
                } else {
                    axis_distance(t.evaluate_step_point(rows.mid)?, o, z)
                };
                let (mut a, mut b) = (rows.sweep, rows.generator);
                let row = if along { &mut b } else { &mut a };
                let high = widen(row, rows.at, rate(rho, distance)?);
                if along {
                    at_most_a_turn(&mut b, high, "the tube angle is over a turn")?;
                } else {
                    at_most_a_turn(&mut a, high, "the sweep is over a turn")?;
                }
                t.adjust_step_uv_domain(k, &Extent2d::from_bounds(a[0], b[0], a[1], b[1]))?;
                Self::from(t)
            }
            _ => {
                return Err(SolidError::Unsupported {
                    operation: OP,
                    reason: "the surface has no analytic continuation",
                })
            }
        };
        // bound 0: the source's label (past a plane's domain), the result's form on it
        let was = AnalyticSupport::of(self);
        let now = AnalyticSupport::check(k, &grown)?;
        let same = match (was, now) {
            (
                Some(AnalyticSupport::Plane {
                    placement: a,
                    uv_scale: sa,
                    ..
                }),
                AnalyticSupport::Plane {
                    placement: b,
                    uv_scale: sb,
                    ..
                },
            ) => a == b && sa == sb,
            (w, n) => w == Some(n),
        };
        if !same {
            return Err(SolidError::invalid_input(
                OP,
                "the continuation is not the source's surface",
            ));
        }
        Ok(grown)
    }
}

/// A revolution's STEP box rows at a face box: the sweep (STEP u) and generator (STEP v)
/// ranges, the side's middle in STEP, and the held coordinate's value there.
struct Rows {
    sweep: [f64; 2],
    generator: [f64; 2],
    mid: Point2,
    at: f64,
}

fn revolution_rows(
    step: &dyn Fn(Point2) -> Result<Point2>,
    dom: &Extent2d,
    side: BoxSide,
    swap: bool,
) -> Result<Rows> {
    let (lo, hi, mid) = (
        step(dom.min())?,
        step(dom.max())?,
        step(side.point(dom, 0.5))?,
    );
    // the side holds the NURBS u for a u side: the generator coordinate when swapped
    let along = side.is_u() == swap;
    Ok(Rows {
        sweep: [lo.x.min(hi.x), lo.x.max(hi.x)],
        generator: [lo.y.min(hi.y), lo.y.max(hi.y)],
        mid,
        at: if along { mid.y } else { mid.x },
    })
}

/// The distance of `p` from the axis through `o` along the unit `z`.
fn axis_distance(p: Point3, o: Point3, z: Vec3) -> f64 {
    let r = p - o;
    (r - z * r.dot(z)).length()
}

/// The angle, degrees, of an arc of length `distance` on a circle of radius `rho`.
fn rate(rho: f64, distance: f64) -> Result<f64> {
    if rho.is_nan() || rho <= 0.0 {
        return Err(SolidError::invalid_input(
            OP,
            "the side's middle is on the axis",
        ));
    }
    Ok((distance / rho).to_degrees())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::vector2d::Vec2;
    use crate::surface::plane::Plane;
    use crate::surface::sphere::{Sphere, SphereOptions};
    use crate::surface::torus::{Torus, TorusOptions};
    use std::f64::consts::{FRAC_PI_2, PI};

    fn k() -> Kernel {
        Kernel::default()
    }

    fn sphere(lon: [f64; 2], lat: [f64; 2]) -> Surface {
        let opts = SphereOptions {
            swap_uv: false,
            inside_out: false,
            circle_nurb: None,
        };
        let anal = Extent2d::from_bounds(lon[0], lat[0], lon[1], lat[1]);
        Surface::from(
            Sphere::new(
                &k(),
                Vec3::new(0.0, 0.0, 1.0),
                Vec3::X,
                Vec3::Y,
                anal,
                2.0,
                opts,
            )
            .unwrap(),
        )
    }

    fn sphere_as(lon: [f64; 2], lat: [f64; 2], swap_uv: bool, inside_out: bool) -> Surface {
        let opts = SphereOptions {
            swap_uv,
            inside_out,
            circle_nurb: None,
        };
        let anal = Extent2d::from_bounds(lon[0], lat[0], lon[1], lat[1]);
        Surface::from(
            Sphere::new(&k(), Vec3::new(0.0, 0.0, 1.0), Vec3::X, Vec3::Y, anal, 2.0, opts)
                .unwrap(),
        )
    }

    fn torus_as(a: [f64; 2], b: [f64; 2], swap_uv: bool, inside_out: bool) -> Surface {
        let anal = Extent2d::from_bounds(a[0], b[0], a[1], b[1]);
        let opts = TorusOptions {
            swap_uv,
            inside_out,
            ..TorusOptions::default()
        };
        Surface::from(
            Torus::new(&k(), Vec3::ZERO, Vec3::X, Vec3::Y, anal, 5.0, 1.0, opts).unwrap(),
        )
    }

    /// Swapped and inside-out spheres and tori continue as the plain ones do: the same STEP
    /// growth, the kind kept (REVIEW2's missing tests).
    #[test]
    fn swapped_and_inside_out_spheres_and_tori_grow_alike() {
        for (swap, inside) in [(false, true), (true, false), (true, true)] {
            let s = sphere_as([0.0, 90.0], [0.0, 30.0], swap, inside);
            let (dom, side) = side_at(&s, true, 30.0);
            let g = s.continued(&k(), &dom, side, 2.0 * PI / 3.0).unwrap();
            assert!(close(steps(&g).max().y, 90.0), "{swap} {inside}: {:?}", steps(&g));
            assert!(matches!(g.kind(), SurfaceKind::Sphere(_)));
            let (dom, side) = side_at(&s, false, 90.0);
            let g = s.continued(&k(), &dom, side, 1.0).unwrap();
            let want = 90.0 + (1.0 / (2.0 * 15f64.to_radians().cos())).to_degrees();
            assert!(close(steps(&g).max().x, want), "{swap} {inside}: {:?}", steps(&g));
            let t = torus_as([0.0, 90.0], [0.0, 90.0], swap, inside);
            let (dom, side) = side_at(&t, true, 90.0);
            let g = t.continued(&k(), &dom, side, FRAC_PI_2).unwrap();
            assert!(close(steps(&g).max().y, 180.0), "{swap} {inside}: {:?}", steps(&g));
            assert!(matches!(g.kind(), SurfaceKind::Torus(_)));
        }
    }

    fn torus(a: [f64; 2], b: [f64; 2]) -> Surface {
        let anal = Extent2d::from_bounds(a[0], b[0], a[1], b[1]);
        Surface::from(
            Torus::new(
                &k(),
                Vec3::ZERO,
                Vec3::X,
                Vec3::Y,
                anal,
                5.0,
                1.0,
                TorusOptions::default(),
            )
            .unwrap(),
        )
    }

    fn steps(s: &Surface) -> Extent2d {
        s.step_uv_domain_with_nurb(&k()).unwrap()
    }

    /// The side whose STEP coordinate (`in_v`: the generator's) is `value` along it.
    fn side_at(s: &Surface, in_v: bool, value: f64) -> (Extent2d, BoxSide) {
        let k = k();
        let dom = s.natural_uv_domain_with_nurb(&k).unwrap();
        let step = |uv| match s.kind() {
            SurfaceKind::Sphere(sp) => sp.convert_uv_from_nurbs_to_step(&k, uv).unwrap(),
            SurfaceKind::Torus(t) => t.convert_uv_from_nurbs_to_step(&k, uv).unwrap(),
            _ => uv,
        };
        let side = BoxSide::ALL
            .into_iter()
            .find(|sd| {
                [0.0, 1.0].iter().all(|&t| {
                    let p = step(sd.point(&dom, t));
                    ((if in_v { p.y } else { p.x }) - value).abs() < 1e-9
                })
            })
            .unwrap();
        (dom, side)
    }

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-11
    }

    #[test]
    fn a_plane_grows_by_the_distance_over_its_scale() {
        let p = Surface::from(Plane::new(
            Point3::ZERO,
            Vec3::X,
            Vec3::Y,
            Vec2::new(2.0, 0.5),
            &Extent2d::from_bounds(0.0, 0.0, 1.0, 4.0),
        ));
        let dom = p.natural_uv_domain_with_nurb(&k()).unwrap();
        let g = p.continued(&k(), &dom, BoxSide::UMax, 3.0).unwrap();
        let d = g.step_uv_domain();
        assert!(close(d.max().x, 2.5) && close(d.min().x, 0.0) && close(d.max().y, 4.0));
        let g = p.continued(&k(), &dom, BoxSide::VMin, 1.0).unwrap();
        assert!(close(g.step_uv_domain().min().y, -2.0));
        assert!(matches!(g.kind(), SurfaceKind::Plane(_)));
    }

    #[test]
    fn a_sphere_grows_round_and_toward_its_pole() {
        let s = sphere([0.0, 90.0], [0.0, 30.0]);
        // round, at the side's middle latitude 15: rho = 2 cos 15
        let (dom, side) = side_at(&s, false, 90.0);
        let g = s.continued(&k(), &dom, side, 1.0).unwrap();
        let want = 90.0 + (1.0 / (2.0 * 15f64.to_radians().cos())).to_degrees();
        assert!(close(steps(&g).max().x, want), "{:?}", steps(&g));
        assert!(matches!(g.kind(), SurfaceKind::Sphere(_)));
        // toward the pole by an arc of pi/3 (60 degrees of radius 2): to the pole exactly
        let (dom, side) = side_at(&s, true, 30.0);
        let g = s.continued(&k(), &dom, side, 2.0 * PI / 3.0).unwrap();
        assert_eq!(steps(&g).max().y, 90.0);
        assert!(s
            .continued(&k(), &dom, side, 2.0 * PI / 3.0 + 1e-6)
            .is_err());
    }

    #[test]
    fn a_torus_grows_round_its_axis_and_its_tube() {
        let t = torus([0.0, 90.0], [0.0, 90.0]);
        // round the tube: r = 1, a quarter more is pi/2
        let (dom, side) = side_at(&t, true, 90.0);
        let g = t.continued(&k(), &dom, side, FRAC_PI_2).unwrap();
        assert!(close(steps(&g).max().y, 180.0), "{:?}", steps(&g));
        // round the axis at the side's middle (tube angle 45): rho = 5 + cos 45
        let (dom, side) = side_at(&t, false, 90.0);
        let g = t.continued(&k(), &dom, side, 1.0).unwrap();
        let want = 90.0 + (1.0 / (5.0 + 45f64.to_radians().cos())).to_degrees();
        assert!(close(steps(&g).max().x, want), "{:?}", steps(&g));
        assert!(matches!(g.kind(), SurfaceKind::Torus(_)));
        assert!(t.continued(&k(), &dom, side, 100.0).is_err());
    }

    #[test]
    fn no_analytic_continuation_for_a_spline_and_a_moved_form_is_refused() {
        let s = sphere([0.0, 90.0], [0.0, 30.0]);
        let (dom, side) = side_at(&s, false, 90.0);
        let form = crate::blend::fillet_section::station::as_b_spline(&k(), &s)
            .unwrap()
            .unwrap()
            .into_owned();
        assert!(Surface::from(form)
            .continued(&k(), &dom, side, 1.0)
            .is_err());
        let SurfaceKind::Sphere(sp) = s.kind() else {
            unreachable!()
        };
        let mut moved = sp.as_ref().clone();
        let shift = crate::geom::axis2_placement::Axis2Placement::new(
            Vec3::new(0.0, 0.0, 1e-6),
            Vec3::X,
            Vec3::Y,
        );
        moved.base.base.transform(&shift, None).unwrap();
        assert!(Surface::from(moved)
            .continued(&k(), &dom, side, 1.0)
            .is_err());
        assert!(s.continued(&k(), &dom, side, -1.0).is_err());
    }
}
