//! The intersection of two exact surfaces over their parameter domains, as the pieces of
//! each exact curve inside both domains with their UV curves.
//!
//! Not upstream. An analytic special case of [`Surface::global_surface_intersect`] (a
//! sphere or cone with a plane, two spheres, ...) gives its exact line or conic made from
//! the analytic surfaces' numbers, clipped at most to an analytic domain, and no UV curves
//! (upstream's `AssociateOutputs`). Upstream's offset never meets them there: it intersects
//! fitted B-spline extensions, whose marched curves are bounded by the domains and carry UV
//! curves. A surface extended exactly keeps its analytic kind, so the special case answers,
//! and a curve may run across the other surface's seam or past its domain.
//! [`SurfaceIntersection::clip_exact_curves`] gives such a curve as the marched intersection
//! gives its boundary curves (`AddBoundaryCurve`): dropped and trimmed onto each surface over
//! its domain (the drop trace, which breaks the curve at the domain's boundaries and seams),
//! each overlap of the two surfaces' pieces one piece of the exact curve with its two UV
//! curves. Each piece is proved on both analytic surfaces ([`curve_lies_on`]) and each UV
//! curve's drop is within the tolerance; anything else is refused with its typed error.
//!
//! Deviations: none (no upstream counterpart). No recursion.

use crate::curve::curve::{Curve, CurveKind};
use crate::error::{Result, SolidError};
use crate::geom::extent1d::Extent1d;
use crate::geom::extent2d::Extent2d;
use crate::kernel::Kernel;
use crate::math::EFF_ZERO;
use crate::nurbs::b_spline_curve_n_lib::TrimMode;
use crate::surface::implicit::{curve_lies_on, AnalyticSupport, Implicit};
use crate::surface::surface::p4::NoObject;
use crate::surface::surface::{Surface, SurfaceIntersection};

/// The operation [`SurfaceIntersection::clip_exact_curves`] names in its errors.
pub const OP: &str = "surface::surface::SurfaceIntersection::clip_exact_curves";

impl SurfaceIntersection {
    /// Each curve of this intersection of `surface` over `uv_domain` with `other` over
    /// `other_uv_domain` that came without its UV curves (an analytic special case's exact
    /// curve) replaced, in place, by its pieces inside both domains, each with its UV
    /// curves on both surfaces (module docs); its type and its largest drop distance are
    /// each piece's. A curve inside neither domain gives no piece. Curves with both UV
    /// curves and degenerate curves (a touching point) are kept as they are.
    ///
    /// # Errors
    ///
    /// [`SolidError::Unsupported`] for a curve without UV curves whose surfaces are not
    /// both analytic, a piece not proved on both, or a UV curve whose drop is farther than
    /// `approx_tol_3d` from the curve; those of the drops, the trims and the proof.
    #[allow(clippy::too_many_arguments)]
    pub fn clip_exact_curves(
        &mut self,
        k: &Kernel,
        surface: &Surface,
        uv_domain: &Extent2d,
        other: &Surface,
        other_uv_domain: &Extent2d,
        approx_tol_3d: f64,
    ) -> Result<()> {
        let old = core::mem::take(self);
        let mut uv1 = old.surface1_uv_curves.into_iter();
        let mut uv2 = old.surface2_uv_curves.into_iter();
        let mut types = old.curve_types.into_iter();
        let mut deviations = old.deviations.into_iter();
        for curve in old.curves_3d {
            let (a, b) = (uv1.next().flatten(), uv2.next().flatten());
            let (Some(curve_type), Some(deviation)) = (types.next(), deviations.next()) else {
                return Err(SolidError::port_invariant(
                    OP,
                    "fewer curve types or deviations than curves",
                ));
            };
            let keep = (a.is_some() && b.is_some()) || curve.is_degenerate(k, EFF_ZERO, None)?;
            if keep {
                self.curves_3d.push(curve);
                self.surface1_uv_curves.push(a);
                self.surface2_uv_curves.push(b);
                self.curve_types.push(curve_type);
                self.deviations.push(deviation);
                continue;
            }
            let implicits = [implicit_of(surface)?, implicit_of(other)?];
            let interval = curve.natural_interval();
            let on1 = surface.drop_and_trim_curve::<NoObject>(
                k,
                uv_domain,
                &curve,
                &interval,
                approx_tol_3d,
                true,
                true,
                None,
            )?;
            let on2 = other.drop_and_trim_curve::<NoObject>(
                k,
                other_uv_domain,
                &curve,
                &interval,
                approx_tol_3d,
                true,
                true,
                None,
            )?;
            if on1.max_drop_to_surf > approx_tol_3d || on2.max_drop_to_surf > approx_tol_3d {
                return Err(SolidError::Unsupported {
                    operation: OP,
                    reason: "a UV curve of an exact intersection curve is farther than the \
                             tolerance from it",
                });
            }
            let drop = on1.max_drop_to_surf.max(on2.max_drop_to_surf);
            for p1 in &on1.uv_curves {
                for p2 in &on2.uv_curves {
                    let (i1, i2) = (p1.natural_interval(), p2.natural_interval());
                    if i1.are_disjoint(&i2, 0.0) {
                        continue;
                    }
                    let mut ivl = i1.intersect(&i2)?;
                    if ivl.length() <= 0.0 {
                        continue;
                    }
                    let mut piece = curve.copy(k)?;
                    let mut piece1 = p1.copy(k)?;
                    let mut piece2 = p2.copy(k)?;
                    // (each trim may snap the interval to a knot: the next gets the snapped one)
                    for c in [&mut piece, &mut piece1, &mut piece2] {
                        trim_to(k, c, &mut ivl)?;
                    }
                    for implicit in &implicits {
                        if !on_implicit(&piece, implicit)? {
                            return Err(SolidError::Unsupported {
                                operation: OP,
                                reason: "a piece of an exact intersection curve is not proved \
                                         on both analytic surfaces",
                            });
                        }
                    }
                    self.curves_3d.push(piece);
                    self.surface1_uv_curves.push(Some(piece1));
                    self.surface2_uv_curves.push(Some(piece2));
                    self.curve_types.push(curve_type);
                    self.deviations.push(drop);
                }
            }
        }
        Ok(())
    }
}

/// The implicit equation of an analytic `surface`.
fn implicit_of(surface: &Surface) -> Result<Implicit> {
    AnalyticSupport::of(surface)
        .map(|s| s.implicit())
        .ok_or(SolidError::Unsupported {
            operation: OP,
            reason: "an intersection curve without UV curves on a surface that is not analytic",
        })
}

/// `curve` trimmed to `ivl` (which the trim may snap).
fn trim_to(k: &Kernel, curve: &mut Curve, ivl: &mut Extent1d) -> Result<()> {
    curve.trim(k, ivl, TrimMode::Notify)
}

/// Whether `curve`, a B-spline kind, lies on `implicit` ([`curve_lies_on`]).
fn on_implicit(curve: &Curve, implicit: &Implicit) -> Result<bool> {
    match curve.kind() {
        CurveKind::BSpline(b) => curve_lies_on(b, implicit),
        CurveKind::Line(l) => curve_lies_on(&*l.base_with_nurb()?, implicit),
        CurveKind::Ellipse(e) => curve_lies_on(e.base(), implicit),
        CurveKind::Hyperbola(h) => curve_lies_on(h.base(), implicit),
        CurveKind::Parabola(p) => curve_lies_on(p.base(), implicit),
        _ => Err(SolidError::Unsupported {
            operation: OP,
            reason: "an exact intersection curve that is not of a B-spline kind",
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::axis2_placement::Axis2Placement;
    use crate::geom::vector2d::Vec2;
    use crate::geom::vector3d::Vec3;
    use crate::surface::plane::Plane;
    use crate::surface::sphere::Sphere;

    /// A whole sphere of radius 5 about the origin (its seam the half plane `y = 0`,
    /// `x > 0`) and a plane through `(-0.5, 0, 0)` tilted about `y`, which cuts it in a
    /// circle round the sphere's axis (the axis meets the plane at `z = -5/3`): across its
    /// seam once, at `z` about `-4.9`.
    fn sphere_and_plane(k: &Kernel) -> (Surface, Surface) {
        let frame = Axis2Placement::new(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        );
        let sphere = Surface::from(Sphere::create_canonical(k, &frame, 5.0).expect("sphere"));
        let plane = Surface::from(
            Plane::create_canonical(&Axis2Placement::new(
                Vec3::new(-0.5, 0.0, 0.0),
                Vec3::new(-0.3, 0.0, 1.0),
                Vec3::new(0.0, 1.0, 0.0),
            ))
            .expect("plane"),
        );
        (sphere, plane)
    }

    fn intersect(
        k: &Kernel,
        s: &Surface,
        o: &Surface,
    ) -> (SurfaceIntersection, Extent2d, Extent2d) {
        let (d1, d2) = (
            s.natural_uv_domain_with_nurb(k).expect("domain"),
            o.natural_uv_domain_with_nurb(k).expect("domain"),
        );
        let ssi = s
            .global_surface_intersect(k, &d1, o, &d2, [true, true], Some(1e-6), None)
            .expect("ssi");
        (ssi, d1, d2)
    }

    #[test]
    fn a_circle_across_the_seam_is_its_pieces_with_their_uv_curves() {
        let k = Kernel::default();
        let (sphere, plane) = sphere_and_plane(&k);
        let (mut ssi, d1, d2) = intersect(&k, &sphere, &plane);
        // the special case: one whole exact circle, no UV curves
        assert_eq!(ssi.curves_3d.len(), 1);
        assert!(ssi.surface1_uv_curves[0].is_none() && ssi.surface2_uv_curves[0].is_none());
        let whole = ssi.curves_3d[0].natural_interval();
        ssi.clip_exact_curves(&k, &sphere, &d1, &plane, &d2, 1e-6)
            .expect("clip");
        // broken at the sphere's seam: two pieces covering the circle
        assert_eq!(ssi.curves_3d.len(), 2);
        let mut covered = 0.0;
        for (i, piece) in ssi.curves_3d.iter().enumerate() {
            let ivl = piece.natural_interval();
            covered += ivl.length();
            let uv1 = ssi.surface1_uv_curves[i]
                .as_ref()
                .expect("uv on the sphere");
            let uv2 = ssi.surface2_uv_curves[i].as_ref().expect("uv on the plane");
            for j in 0..=8 {
                let t = ivl.evaluate(f64::from(j) / 8.0);
                let p = piece.evaluate_point(&k, t).expect("point");
                for (uv, s) in [(uv1, &sphere), (uv2, &plane)] {
                    let q = uv.evaluate_point(&k, t).expect("uv");
                    let on = s.evaluate_point(&k, Vec2::new(q.x, q.y)).expect("on");
                    assert!((on - p).length() <= 1e-6, "{} off", (on - p).length());
                }
            }
            // each piece proved on both analytic surfaces
            for s in [&sphere, &plane] {
                assert!(on_implicit(piece, &implicit_of(s).unwrap()).unwrap());
            }
        }
        assert!((covered - whole.length()).abs() <= 1e-12 * whole.length());
    }

    #[test]
    fn a_circle_moved_off_its_sphere_is_refused() {
        let k = Kernel::default();
        let (sphere, plane) = sphere_and_plane(&k);
        let (ssi, _, _) = intersect(&k, &sphere, &plane);
        let circle = &ssi.curves_3d[0];
        let on = implicit_of(&sphere).unwrap();
        assert!(on_implicit(circle, &on).unwrap());
        // the same sphere 1e-9 larger: the circle is not on it
        let Implicit::Sphere { centre, radius } = on else {
            panic!("a sphere")
        };
        let moved = Implicit::Sphere {
            centre,
            radius: radius + 1e-9,
        };
        assert!(!on_implicit(circle, &moved).unwrap());
    }

    #[test]
    fn a_curve_on_a_surface_that_is_not_analytic_is_refused() {
        let k = Kernel::default();
        let (sphere, plane) = sphere_and_plane(&k);
        let (mut ssi, d1, d2) = intersect(&k, &sphere, &plane);
        // the plane's own B-spline net: the same points, no analytic kind
        let crate::surface::surface::SurfaceKind::Plane(p) = plane.kind() else {
            panic!("a plane")
        };
        let nurbs = Surface::from(p.base().clone());
        let err = ssi
            .clip_exact_curves(&k, &sphere, &d1, &nurbs, &d2, 1e-6)
            .unwrap_err();
        assert!(matches!(err, SolidError::Unsupported { operation: OP, .. }));
    }
}
