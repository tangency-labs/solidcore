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
//! curves. A closed curve that crosses a seam or a domain boundary is first turned to start
//! at a crossing, so that its pieces end only at crossings and each is given once.
//!
//! What is checked, and how: each piece's control points are on both analytic surfaces to
//! the rounding of evaluating their implicit equations ([`curve_lies_on`]; not covering the
//! rounding of the knot insertion or the trim that made the piece). Each UV curve's drop
//! distance is the drop trace's *sampled* largest distance (`DropAndTrimCurve`'s
//! `rdMaxDropToSurf`), within the tolerance: it is not certified between the samples.
//! Anything else is refused with its typed error.
//!
//! Deviations: none (no upstream counterpart). No recursion.

use crate::curve::curve::{Curve, CurveKind};
use crate::curve::ellipse::{Ellipse, EllipseOptions};
use crate::error::{Result, SolidError};
use crate::geom::extent1d::Extent1d;
use crate::geom::extent2d::Extent2d;
use crate::kernel::Kernel;
use crate::math::{EFF_ZERO, EFF_ZERO_SQRT};
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
    /// each piece's. Each piece is given once, and ends only where the curve crosses a
    /// seam or a domain boundary of either surface: a closed curve crossing one is first
    /// turned to start at such a crossing ([`turned_to`]), so its own parameter origin
    /// splits nothing; a closed curve crossing none is one closed piece. A curve inside
    /// neither domain gives no piece. Curves with both UV curves and degenerate curves (a
    /// touching point) are kept as they are.
    ///
    /// # Errors
    ///
    /// [`SolidError::Unsupported`] for a curve without UV curves whose surfaces are not
    /// both analytic, a piece not on both ([`curve_lies_on`]), a UV curve whose sampled
    /// drop distance is more than `approx_tol_3d`, or a closed curve to be turned that is
    /// not an ellipse; [`SolidError::PortInvariant`] for a drop whose pieces overlap or a
    /// turned ellipse that does not start at the crossing; those of the drops, the trims
    /// and the check.
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
        let pair = Pair {
            surface,
            uv_domain,
            other,
            other_uv_domain,
            approx_tol_3d,
        };
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
            let mut drops = pair.drops(k, &curve)?;
            let mut curve = curve;
            if let Some(t) = drops.first_inner_break(&curve.natural_interval()) {
                if is_closed_curve(k, &curve, approx_tol_3d)? {
                    curve = turned_to(k, &curve, t, approx_tol_3d)?;
                    drops = pair.drops(k, &curve)?;
                }
            }
            for (mut ivl, p1, p2) in drops.overlaps() {
                let mut piece = curve.copy(k)?;
                let mut piece1 = p1.copy(k)?;
                let mut piece2 = p2.copy(k)?;
                // (each trim may snap the interval to a knot: the next gets the snapped one)
                for c in [&mut piece, &mut piece1, &mut piece2] {
                    c.trim(k, &mut ivl, TrimMode::Notify)?;
                }
                for implicit in &implicits {
                    if !on_implicit(&piece, implicit)? {
                        return Err(SolidError::Unsupported {
                            operation: OP,
                            reason: "a piece of an exact intersection curve is not on both \
                                     analytic surfaces",
                        });
                    }
                }
                self.curves_3d.push(piece);
                self.surface1_uv_curves.push(Some(piece1));
                self.surface2_uv_curves.push(Some(piece2));
                self.curve_types.push(curve_type);
                self.deviations.push(drops.drop);
            }
        }
        Ok(())
    }
}

/// The two surfaces of [`SurfaceIntersection::clip_exact_curves`] over their domains.
struct Pair<'a> {
    surface: &'a Surface,
    uv_domain: &'a Extent2d,
    other: &'a Surface,
    other_uv_domain: &'a Extent2d,
    approx_tol_3d: f64,
}

/// A curve's UV pieces on each surface of a [`Pair`], each list in parameter order and
/// without overlaps, and their largest (sampled) drop distance.
struct Drops {
    on1: Vec<Curve>,
    on2: Vec<Curve>,
    drop: f64,
}

impl Pair<'_> {
    /// `curve` dropped and trimmed onto both surfaces over their domains.
    fn drops(&self, k: &Kernel, curve: &Curve) -> Result<Drops> {
        let (on1, d1) = self.drop_onto(k, self.surface, self.uv_domain, curve)?;
        let (on2, d2) = self.drop_onto(k, self.other, self.other_uv_domain, curve)?;
        Ok(Drops {
            on1,
            on2,
            drop: d1.max(d2),
        })
    }

    /// `curve`'s UV pieces on `surface` over `domain`, in parameter order, and their
    /// largest sampled drop distance.
    fn drop_onto(
        &self,
        k: &Kernel,
        surface: &Surface,
        domain: &Extent2d,
        curve: &Curve,
    ) -> Result<(Vec<Curve>, f64)> {
        let on = surface.drop_and_trim_curve::<NoObject>(
            k,
            domain,
            curve,
            &curve.natural_interval(),
            self.approx_tol_3d,
            true,
            true,
            None,
        )?;
        if on.max_drop_to_surf > self.approx_tol_3d {
            return Err(SolidError::Unsupported {
                operation: OP,
                reason: "a UV curve of an exact intersection curve is farther than the \
                         tolerance from it",
            });
        }
        let mut pieces = on.uv_curves;
        pieces.sort_by(|a, b| {
            a.natural_interval()
                .min()
                .total_cmp(&b.natural_interval().min())
        });
        // the trace breaks the curve at the crossings: consecutive pieces at most meet
        let mut ends = pieces.iter().map(Curve::natural_interval);
        if let Some(mut before) = ends.next() {
            for ivl in ends {
                if ivl.min() < before.max() {
                    return Err(SolidError::port_invariant(
                        OP,
                        "the drop of an exact intersection curve gave overlapping pieces",
                    ));
                }
                before = ivl;
            }
        }
        Ok((pieces, on.max_drop_to_surf))
    }
}

impl Drops {
    /// The first end of a piece on either surface inside `interval`: where the curve
    /// crosses a seam or a domain boundary other than at its own ends. An end nearer an
    /// end of `interval` than the drop's own break distance (`SM_EFF_ZERO_SQRT / 100`,
    /// below which `FindBoundaryStartPoints` takes a crossing for the end) is the end.
    fn first_inner_break(&self, interval: &Extent1d) -> Option<f64> {
        let near = EFF_ZERO_SQRT / 100.0;
        self.on1
            .iter()
            .chain(&self.on2)
            .flat_map(|c| {
                let ivl = c.natural_interval();
                [ivl.min(), ivl.max()]
            })
            .filter(|&t| t > interval.min() + near && t < interval.max() - near)
            .min_by(f64::total_cmp)
    }

    /// Each overlap of a piece on the first surface with one on the second, once, in
    /// parameter order: a sweep over both ordered lists.
    fn overlaps(&self) -> Vec<(Extent1d, &Curve, &Curve)> {
        let mut out = Vec::new();
        let (mut i, mut j) = (0, 0);
        while let (Some(p1), Some(p2)) = (self.on1.get(i), self.on2.get(j)) {
            let (i1, i2) = (p1.natural_interval(), p2.natural_interval());
            let (lo, hi) = (i1.min().max(i2.min()), i1.max().min(i2.max()));
            if hi > lo {
                out.push((Extent1d::from_min_max(lo, hi), p1, p2));
            }
            // the piece that ends first meets nothing further on the other list
            if i1.max() <= i2.max() {
                i += 1;
            } else {
                j += 1;
            }
        }
        out
    }
}

/// Whether `curve` is closed: an ellipse all the way round, or another curve whose ends
/// are within `tol` ([`Curve::is_closed`]).
fn is_closed_curve(k: &Kernel, curve: &Curve, tol: f64) -> Result<bool> {
    match curve.kind() {
        // the full turn as the ellipse's own largest domain takes it
        CurveKind::Ellipse(e) => Ok(e.step_interval().length() > 360.0 - EFF_ZERO),
        _ => curve.is_closed(k, &curve.natural_interval(), tol),
    }
}

/// The closed `curve` (an ellipse all the way round) turned to start at its parameter
/// `t`: the same ellipse (placement, radii, class, sense) over the full turn from the
/// angle of `t`, so its points are the same to rounding and its parameter origin is at
/// `t`.
///
/// # Errors
///
/// [`SolidError::Unsupported`] for a closed curve that is not an ellipse;
/// [`SolidError::PortInvariant`] for a turned ellipse without NURB or not starting
/// within `tol` of the curve's point at `t`; those of the angle and the constructor.
fn turned_to(k: &Kernel, curve: &Curve, t: f64, tol: f64) -> Result<Curve> {
    let CurveKind::Ellipse(e) = curve.kind() else {
        return Err(SolidError::Unsupported {
            operation: OP,
            reason: "a closed exact intersection curve crossing a seam or domain boundary \
                     that is not an ellipse",
        });
    };
    let angle = e.convert_t_from_nurbs_to_step(k, t)?;
    // an inside-out ellipse runs from its largest angle down
    let domain = if e.inside_out() {
        Extent1d::from_min_max(angle - 360.0, angle)
    } else {
        Extent1d::from_min_max(angle, angle + 360.0)
    };
    let opts = EllipseOptions {
        dim: e.base().dim(),
        nurb: None,
        inside_out: e.inside_out(),
    };
    let at = e.position();
    let (o, x, y) = (at.origin(), at.x_axis(), at.y_axis());
    let turned = if e.is_circle_class() {
        Ellipse::new_circle(k, o, x, y, domain, e.x_radius(), &opts)?
    } else {
        Ellipse::new(k, o, x, y, domain, e.x_radius(), e.y_radius(), &opts)?
    };
    if turned.base().nurb.is_none() {
        return Err(SolidError::port_invariant(
            OP,
            "a turned exact intersection ellipse has no NURB",
        ));
    }
    let turned = Curve::from(turned);
    let start = turned.evaluate_point(k, turned.natural_interval().min())?;
    if start.distance(curve.evaluate_point(k, t)?) > tol {
        return Err(SolidError::port_invariant(
            OP,
            "a turned exact intersection ellipse does not start at the crossing",
        ));
    }
    Ok(turned)
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

    /// A whole sphere of radius 5 about the origin, its seam the half plane `y = 0`,
    /// `x > 0` (`u` 0 and 1).
    fn sphere(k: &Kernel) -> Surface {
        let frame = Axis2Placement::new(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        );
        Surface::from(Sphere::create_canonical(k, &frame, 5.0).expect("sphere"))
    }

    fn plane(origin: Vec3, x: Vec3, y: Vec3) -> Surface {
        Surface::from(Plane::create_canonical(&Axis2Placement::new(origin, x, y)).expect("plane"))
    }

    /// The sphere and a plane through `(-0.5, 0, 0)` tilted about `y`, which cuts it in a
    /// circle round the sphere's axis (the axis meets the plane at `z = -5/3`) that starts
    /// off the seam and crosses it once, at `z` about `-4.9`.
    fn sphere_and_plane(k: &Kernel) -> (Surface, Surface) {
        let tilted = plane(
            Vec3::new(-0.5, 0.0, 0.0),
            Vec3::new(-0.3, 0.0, 1.0),
            Vec3::new(0.0, 1.0, 0.0),
        );
        (sphere(k), tilted)
    }

    /// The plane `z = 2`: its circle on the sphere starts on the seam.
    fn plane_z2() -> Surface {
        plane(
            Vec3::new(0.0, 0.0, 2.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        )
    }

    fn intersect(
        k: &Kernel,
        s: &Surface,
        o: &Surface,
    ) -> (SurfaceIntersection, Extent2d, Extent2d) {
        let d2 = o.natural_uv_domain_with_nurb(k).expect("domain");
        intersect_over(k, s, o, d2)
    }

    fn intersect_over(
        k: &Kernel,
        s: &Surface,
        o: &Surface,
        d2: Extent2d,
    ) -> (SurfaceIntersection, Extent2d, Extent2d) {
        let d1 = s.natural_uv_domain_with_nurb(k).expect("domain");
        let ssi = s
            .global_surface_intersect(k, &d1, o, &d2, [true, true], Some(1e-6), None)
            .expect("ssi");
        (ssi, d1, d2)
    }

    /// Whether `uv` is on a side of `d`.
    fn on_side(uv: Vec2, d: &Extent2d) -> bool {
        let near = |a: f64, b: f64| (a - b).abs() <= 1e-9 * (1.0 + b.abs());
        near(uv.x, d.u_min())
            || near(uv.x, d.u_max())
            || near(uv.y, d.v_min())
            || near(uv.y, d.v_max())
    }

    /// The pieces of the special case's one exact circle of `sphere` with `other` over
    /// `d2`, checked: each piece's UV curves are on the piece (sampled) and the piece is
    /// on both surfaces ([`curve_lies_on`]); no two pieces overlap; each end of a piece
    /// is a crossing (on a side of either domain) unless `closed_free`, when the circle
    /// crosses nothing and is one closed piece. The pieces' 3D ends, in order.
    fn checked_pieces(
        k: &Kernel,
        sphere: &Surface,
        other: &Surface,
        d2: Extent2d,
        closed_free: bool,
    ) -> Vec<(Vec3, Vec3)> {
        let (mut ssi, d1, d2) = intersect_over(k, sphere, other, d2);
        // the special case: one whole exact circle, no UV curves
        assert_eq!(ssi.curves_3d.len(), 1);
        assert!(ssi.surface1_uv_curves[0].is_none() && ssi.surface2_uv_curves[0].is_none());
        ssi.clip_exact_curves(k, sphere, &d1, other, &d2, 1e-6)
            .expect("clip");
        let mut ends = Vec::new();
        let mut intervals = Vec::new();
        for (i, piece) in ssi.curves_3d.iter().enumerate() {
            let ivl = piece.natural_interval();
            let uv1 = ssi.surface1_uv_curves[i]
                .as_ref()
                .expect("uv on the sphere");
            let uv2 = ssi.surface2_uv_curves[i].as_ref().expect("uv on the other");
            for j in 0..=8 {
                let t = ivl.evaluate(f64::from(j) / 8.0);
                let p = piece.evaluate_point(k, t).expect("point");
                for (uv, s) in [(uv1, sphere), (uv2, other)] {
                    let q = uv.evaluate_point(k, t).expect("uv");
                    let on = s.evaluate_point(k, Vec2::new(q.x, q.y)).expect("on");
                    assert!((on - p).length() <= 1e-6, "{} off", (on - p).length());
                }
            }
            for s in [sphere, other] {
                assert!(on_implicit(piece, &implicit_of(s).unwrap()).unwrap());
            }
            for t in [ivl.min(), ivl.max()] {
                let (a, b) = (
                    uv1.evaluate_point(k, t).expect("uv"),
                    uv2.evaluate_point(k, t).expect("uv"),
                );
                let crossing =
                    on_side(Vec2::new(a.x, a.y), &d1) || on_side(Vec2::new(b.x, b.y), &d2);
                assert!(
                    crossing || closed_free,
                    "piece {i} ends at {t}, not at a crossing"
                );
            }
            intervals.push(ivl);
            ends.push((
                piece.evaluate_point(k, ivl.min()).expect("start"),
                piece.evaluate_point(k, ivl.max()).expect("end"),
            ));
        }
        for (i, a) in intervals.iter().enumerate() {
            for b in intervals.iter().skip(i + 1) {
                assert!(a.max() <= b.min() || b.max() <= a.min(), "pieces overlap");
            }
        }
        ends
    }

    fn same(a: Vec3, b: Vec3) -> bool {
        (a - b).length() <= 1e-12 * 5.0
    }

    #[test]
    fn a_circle_across_the_seam_once_is_one_closed_piece_from_the_seam() {
        let k = Kernel::default();
        let (sphere, tilted) = sphere_and_plane(&k);
        let d2 = tilted.natural_uv_domain_with_nurb(&k).expect("domain");
        let ends = checked_pieces(&k, &sphere, &tilted, d2, false);
        // one seam-to-seam piece, all the way round (not split at the circle's own start)
        assert_eq!(ends.len(), 1);
        let (start, end) = ends[0];
        assert!(same(start, end));
        assert!(start.y.abs() <= 1e-9 && start.x > 0.0 && start.z < -4.8);
    }

    #[test]
    fn a_circle_starting_on_the_seam_is_one_piece_given_once() {
        let k = Kernel::default();
        let sphere = sphere(&k);
        let h = plane_z2();
        let d2 = h.natural_uv_domain_with_nurb(&k).expect("domain");
        let ends = checked_pieces(&k, &sphere, &h, d2, false);
        assert_eq!(ends.len(), 1);
        let (start, end) = ends[0];
        assert!(same(start, end));
        assert!(same(start, Vec3::new(21.0_f64.sqrt(), 0.0, 2.0)));
    }

    #[test]
    fn a_circle_starting_on_the_seam_over_a_quadrant_is_its_quarter_once() {
        let k = Kernel::default();
        let sphere = sphere(&k);
        let h = plane_z2();
        let quadrant = Extent2d::from_points(Vec2::new(0.0, 0.0), Vec2::new(10.0, 10.0));
        let ends = checked_pieces(&k, &sphere, &h, quadrant, false);
        assert_eq!(ends.len(), 1);
        let r = 21.0_f64.sqrt();
        let (start, end) = ends[0];
        assert!(same(start, Vec3::new(r, 0.0, 2.0)), "{start:?}");
        assert!(same(end, Vec3::new(0.0, r, 2.0)), "{end:?}");
    }

    #[test]
    fn a_circle_across_the_seam_twice_is_two_pieces_between_the_crossings() {
        let k = Kernel::default();
        let sphere = sphere(&k);
        // the plane x = 1: a circle off the axis, across the seam at z = +-sqrt(24)
        let v = plane(
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        );
        let d2 = v.natural_uv_domain_with_nurb(&k).expect("domain");
        let ends = checked_pieces(&k, &sphere, &v, d2, false);
        assert_eq!(ends.len(), 2);
        let z = 24.0_f64.sqrt();
        let (top, bottom) = (Vec3::new(1.0, 0.0, z), Vec3::new(1.0, 0.0, -z));
        // consecutive: the one's end the other's start, each a seam crossing
        assert!(same(ends[0].1, ends[1].0) && same(ends[1].1, ends[0].0));
        for (a, b) in ends {
            assert!((same(a, top) && same(b, bottom)) || (same(a, bottom) && same(b, top)));
        }
    }

    #[test]
    fn a_circle_crossing_no_seam_is_one_closed_piece() {
        let k = Kernel::default();
        let sphere = sphere(&k);
        // the plane y = 1: a circle on the side of the sphere away from its seam
        let w = plane(
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(1.0, 0.0, 0.0),
        );
        let d2 = w.natural_uv_domain_with_nurb(&k).expect("domain");
        let ends = checked_pieces(&k, &sphere, &w, d2, true);
        assert_eq!(ends.len(), 1);
        assert!(same(ends[0].0, ends[0].1));
    }

    /// The drop itself: a circle starting on the sphere's seam, on a `v` iso line, drops
    /// to that iso line once (the seam's other point is the same point of it).
    #[test]
    fn the_drop_of_an_iso_circle_starting_on_the_seam_gives_it_once() {
        let k = Kernel::default();
        let sphere = sphere(&k);
        let h = plane_z2();
        let (ssi, d1, _) = intersect(&k, &sphere, &h);
        let circle = &ssi.curves_3d[0];
        let on = sphere
            .drop_and_trim_curve::<NoObject>(
                &k,
                &d1,
                circle,
                &circle.natural_interval(),
                1e-6,
                true,
                true,
                None,
            )
            .expect("drop");
        assert_eq!(on.uv_curves.len(), 1);
    }

    /// A full circular edge (a hemisphere's base circle `z = 2` when its base plane moves
    /// there) is covered by one piece: trimming that piece to the closed edge's whole span
    /// gives the circle, so the edge is mated.
    #[test]
    fn a_full_circular_edge_is_mateable_on_its_one_piece() {
        let k = Kernel::default();
        let sphere = sphere(&k);
        let h = plane_z2();
        let (mut ssi, d1, d2) = intersect(&k, &sphere, &h);
        let whole = ssi.curves_3d[0].natural_interval();
        ssi.clip_exact_curves(&k, &sphere, &d1, &h, &d2, 1e-6)
            .expect("clip");
        assert_eq!(ssi.curves_3d.len(), 1);
        let mut piece = ssi.curves_3d[0].copy(&k).expect("copy");
        let mut span = piece.natural_interval();
        assert!(span == whole, "{span:?}");
        piece
            .trim(&k, &mut span, TrimMode::Notify)
            .expect("trimmed to the closed edge");
        let r = 21.0_f64.sqrt();
        for j in 0..=16 {
            let p = piece
                .evaluate_point(&k, span.evaluate(f64::from(j) / 16.0))
                .expect("p");
            assert!(((p - Vec3::new(0.0, 0.0, 2.0)).length() - r).abs() <= 1e-12 * 5.0);
            assert!((p.z - 2.0).abs() <= 1e-12 * 5.0);
        }
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
