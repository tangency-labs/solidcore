//! The setback fill at order 0 with no setback: a skew polygon with straight sides, closed
//! by the setback split's interior map (not the fork's; tangency-labs/solidcore#2345, the
//! owner's answer 3: chamfers meeting at a vertex whose rails' crossings are not in one
//! plane).
//!
//! The blends stop at their crossings, so the hole is the polygon `X_0 .. X_(n-1)` of their
//! end sections, each a straight line, and nothing of it lies on the faces: there are no
//! setback patches, and each side is its edge's inner curve. The interior map F of the
//! regular n-gon (vertices `B_m = (cos 2 pi m / n, sin 2 pi m / n)`, as [`fill_hole`]'s) is
//! held to the side `X_m`-`X_(m+1)` along its edge m from `B_m` to `B_(m+1)`, linearly, and
//! is otherwise the one of least bending energy over the n quadrilaterals the domain splits
//! into at its centre ([`fill_hole`]'s energy and solves), at the least degree whose
//! conditions can be met. It meets each blend to G0 exactly: along the side it is the side.
//! At four sides it is the bilinear map of the corners (the one quadratic linear along all
//! four sides of a square).
//!
//! The corner face is one face: F is one polynomial, so over the n-gon's bounding square
//! `[-1, 1]^2` it is one tensor Bezier patch, its parameters the domain's own, and the face
//! is that patch trimmed by the n-gon, whose sides are straight in the parameters and map
//! onto the polygon's sides. The n quadrilaterals' patches are kept too: their control
//! points bound the face (the patch over the square runs past the n-gon).
//!
//! | Rust               | Role                                                  |
//! |--------------------|-------------------------------------------------------|
//! | [`PolygonFill`]    | the fill: its patch over the square, its parts, its degree |
//! | [`fill_polygon`]   | the fill of a skew polygon at order 0                  |
//!
//! Deviations: not a port. No recursion; the degree is tried from 1 to 9.
//!
//! [`fill_hole`]: super::fill_hole

use super::{
    affine_solutions, el, el_mut, el_ref, first, least_squares, monomials, over_patch, second,
    smax, Affine, Net, Refusal, Tensor, GW, GX,
};
use crate::error::Result as SolidResult;
use crate::geom::vector2d::Point2;
use crate::geom::vector3d::sqrt;
use crate::math::{Point3, PI};
use crate::nurbs::b_spline_curve::KnotType;
use crate::nurbs::b_spline_surface::{BSplineSurface, BSplineSurfaceForm, CanonicalSurfaceData};

/// The most sides a polygon is filled for (a bound on input-driven work).
const MAX_SIDES: usize = 64;

/// A skew polygon's order-0 fill ([`fill_polygon`]).
#[derive(Debug, Clone)]
pub struct PolygonFill {
    /// F over the n-gon's bounding square `[-1, 1]^2`, its parameters the domain's: corner
    /// m at `B_m` ([`Self::domain`]), side m the straight line from `B_m` to `B_(m+1)`. One
    /// tensor Bezier patch of degree `(degree, degree)`.
    pub surface: BSplineSurface,
    /// The n quadrilaterals of the n-gon at its centre, F over each (the setback split's
    /// interior parts, `(s, t)`: `(0, 0)` the centre, `(1, 1)` `B_q`): the face's own
    /// pieces, whose control points bound it.
    pub parts: Vec<BSplineSurface>,
    /// The domain's vertices `B_m`.
    pub domain: Vec<Point2>,
    /// F's degree.
    pub degree: usize,
}

/// The fill of the skew polygon `corners` (three or more, in order round it; its sides
/// straight) at order 0: the interior map of least bending energy that is each side along
/// its edge of the regular n-gon, at the least degree from 1 that meets that to `tol`
/// (module docs).
///
/// # Errors
///
/// [`Refusal`] with `SM_ERR_INVALID_INPUT` for fewer than three or more than 64 corners or a
/// tolerance not positive; with `SM_ERR` for conditions no degree up to 9 meets, or a fill
/// that folds. A failing call's error, without a reason.
#[allow(clippy::too_many_lines)] // one construction
pub fn fill_polygon(corners: &[Point3], tol: f64) -> Result<PolygonFill, Refusal> {
    let n = corners.len();
    if !(3..=MAX_SIDES).contains(&n) || !(tol > 0.0) {
        return Err(Refusal::invalid(
            "SmBlendSetback: a polygon's fill needs three to 64 corners and a positive tolerance",
        ));
    }
    let nf = n as f64;
    let bs: Vec<Point2> = (0..n)
        .map(|m| {
            Point2::new(
                (2.0 * PI * m as f64 / nf).cos(),
                (2.0 * PI * m as f64 / nf).sin(),
            )
        })
        .collect();

    for d in 1..=9 {
        let mono = monomials(d);
        let mc = mono.len();
        let unknowns = 3 * mc;

        // each side's Bernstein coefficients per monomial (degree d along the edge)
        let mut edge_basis: Vec<Vec<Vec<f64>>> = Vec::with_capacity(n);
        for m in 0..n {
            let b0 = el(&bs, m)?;
            let b1 = el(&bs, (m + 1) % n)?;
            let mut eb = Vec::with_capacity(mc);
            for &(a, b) in &mono {
                eb.push(super::along_line(a, b, &[b0.x, b1.x], &[b0.y, b1.y], d)?);
            }
            edge_basis.push(eb);
        }
        // each interior part's net per monomial, degree (d, d): corners (0, 0) the centre,
        // (1, 0) the middle of edge q, (1, 1) B_q, (0, 1) the middle of edge q - 1
        let mut part_basis: Vec<Vec<Tensor>> = Vec::with_capacity(n);
        for q in 0..n {
            let bq = el(&bs, q)?;
            let bn = el(&bs, (q + 1) % n)?;
            let bp = el(&bs, (q + n - 1) % n)?;
            let mq = Point2::new((bq.x + bn.x) * 0.5, (bq.y + bn.y) * 0.5);
            let mp = Point2::new((bq.x + bp.x) * 0.5, (bq.y + bp.y) * 0.5);
            let tx = Tensor {
                p: 1,
                q: 1,
                c: vec![0.0, mp.x, mq.x, bq.x],
            };
            let ty = Tensor {
                p: 1,
                q: 1,
                c: vec![0.0, mp.y, mq.y, bq.y],
            };
            let mut parts = Vec::with_capacity(mc);
            for &(a, b) in &mono {
                parts.push(over_patch(a, b, &tx, &ty, d)?);
            }
            part_basis.push(parts);
        }
        let nets = |xs: &[f64]| -> SolidResult<Vec<Net>> {
            let mut out = Vec::with_capacity(n);
            for parts in &part_basis {
                out.push(net_of(xs, parts, mc, d)?);
            }
            Ok(out)
        };
        // the conditions: F's coefficients along each edge less the side's, raised to d
        let conditions = |xs: &[f64]| -> SolidResult<Vec<f64>> {
            let mut out = Vec::with_capacity(n * (d + 1) * 3);
            for (m, eb) in edge_basis.iter().enumerate() {
                let p0 = *el_ref(corners, m)?;
                let p1 = *el_ref(corners, (m + 1) % n)?;
                for j in 0..=d {
                    let t = j as f64 / d as f64;
                    let side = p0 + (p1 - p0) * t;
                    for (c, want) in [side.x, side.y, side.z].into_iter().enumerate() {
                        let mut s = -want;
                        for (jj, ebj) in eb.iter().enumerate() {
                            let x = el(xs, c * mc + jj)?;
                            if x != 0.0 {
                                s += x * el(ebj, j)?;
                            }
                        }
                        out.push(s);
                    }
                }
            }
            Ok(out)
        };
        // the bending energy's terms, each squared and summed
        let energy = |parts: &[Net]| -> SolidResult<Vec<f64>> {
            let mut out = Vec::new();
            for net in parts {
                for (&gx1, &gw1) in GX.iter().zip(&GW) {
                    for (&gx2, &gw2) in GX.iter().zip(&GW) {
                        let [ruu, ruv, rvv] = second(net, gx1, gx2)?;
                        let w = sqrt(gw1 * gw2);
                        for (uu, (uv, vv)) in [ruu.x, ruu.y, ruu.z]
                            .into_iter()
                            .zip([ruv.x, ruv.y, ruv.z].into_iter().zip([rvv.x, rvv.y, rvv.z]))
                        {
                            out.push(w * uu);
                            out.push(w * sqrt(2.0) * uv);
                            out.push(w * vv);
                        }
                    }
                }
            }
            Ok(out)
        };

        // both affine in x: their matrices column by column
        let x0 = vec![0.0; unknowns];
        let c0 = conditions(&x0)?;
        let e0 = energy(&nets(&x0)?)?;
        let (lc, le) = (c0.len(), e0.len());
        let mut ac = vec![0.0; lc * unknowns];
        let mut ae = vec![0.0; le * unknowns];
        for u in 0..unknowns {
            let mut xs = vec![0.0; unknowns];
            *el_mut(&mut xs, u)? = 1.0;
            let cu = conditions(&xs)?;
            let eu = energy(&nets(&xs)?)?;
            for i in 0..lc {
                *el_mut(&mut ac, i * unknowns + u)? = el(&cu, i)? - el(&c0, i)?;
            }
            for i in 0..le {
                *el_mut(&mut ae, i * unknowns + u)? = el(&eu, i)? - el(&e0, i)?;
            }
        }
        let mut rhs = vec![0.0; lc];
        let mut size = 0.0;
        for (ri, &ci) in rhs.iter_mut().zip(&c0) {
            *ri = -ci;
            size = smax(size, ci.abs());
        }
        let Some(Affine { x0: xp, z }) =
            affine_solutions(ac, rhs, lc, unknowns, 1.0e-9 * (1.0 + size))?
        else {
            continue;
        };
        // the energy over x = xp + Z y; with no freedom left the conditions alone are F
        let mut xs = xp;
        if !z.is_empty() {
            let free = z.len();
            let mut mm = vec![0.0; le * free];
            let mut rr = vec![0.0; le];
            for i in 0..le {
                let mut s = el(&e0, i)?;
                for (u, &xu) in xs.iter().enumerate() {
                    s += el(&ae, i * unknowns + u)? * xu;
                }
                *el_mut(&mut rr, i)? = -s;
                for (f, zf) in z.iter().enumerate() {
                    let mut t = 0.0;
                    for (u, &zu) in zf.iter().enumerate() {
                        t += el(&ae, i * unknowns + u)? * zu;
                    }
                    *el_mut(&mut mm, i * free + f)? = t;
                }
            }
            let Some(y) = least_squares(mm, rr, le, free, 1.0e-12)? else {
                continue;
            };
            for (zf, &yf) in z.iter().zip(&y) {
                for (xu, &zu) in xs.iter_mut().zip(zf) {
                    *xu += zu * yf;
                }
            }
        }
        let mut worst = 0.0;
        for &c in &conditions(&xs)? {
            worst = smax(worst, c.abs());
        }
        if !(worst <= 1.0e-3 * tol) {
            continue;
        }

        // no fold: every part's normal keeps to one side, and to the centre part's side
        let part_nets = nets(&xs)?;
        let (_, su, sv) = first(el_ref(&part_nets, 0)?, 0.0, 0.0)?;
        let at_centre = su.cross(sv);
        for net in &part_nets {
            for a in 0..=8 {
                for b in 0..=8 {
                    let (_, su, sv) = first(net, f64::from(a) / 8.0, f64::from(b) / 8.0)?;
                    if !(su.cross(sv).dot(at_centre) > 0.0) {
                        return Err(Refusal::sm_err("SmBlendSetback: the polygon's fill folds"));
                    }
                }
            }
        }

        // the patch over the square, and the parts
        let square = |c: [f64; 4]| Tensor {
            p: 1,
            q: 1,
            c: c.to_vec(),
        };
        let (sx, sy) = (
            square([-1.0, -1.0, 1.0, 1.0]),
            square([-1.0, 1.0, -1.0, 1.0]),
        );
        let mut over = Vec::with_capacity(mc);
        for &(a, b) in &mono {
            over.push(over_patch(a, b, &sx, &sy, d)?);
        }
        let surface = bezier(&net_of(&xs, &over, mc, d)?, [-1.0, 1.0])?;
        let mut parts = Vec::with_capacity(n);
        for net in &part_nets {
            parts.push(bezier(net, [0.0, 1.0])?);
        }
        return Ok(PolygonFill {
            surface,
            parts,
            domain: bs,
            degree: d,
        });
    }
    Err(Refusal::sm_err(
        "SmBlendSetback: no interior map up to degree 9 meets the polygon's sides",
    ))
}

/// F's net over a patch: the coefficients `xs[c mc + j]` of the monomials' nets `basis`.
fn net_of(xs: &[f64], basis: &[Tensor], mc: usize, d: usize) -> SolidResult<Net> {
    let mut net = Net {
        p: d,
        q: d,
        pts: vec![[0.0, 0.0, 0.0, 1.0]; (d + 1) * (d + 1)],
    };
    for c in 0..3 {
        for (j, part) in basis.iter().enumerate() {
            let x = el(xs, c * mc + j)?;
            if x == 0.0 {
                continue;
            }
            for (h, &pc) in net.pts.iter_mut().zip(&part.c) {
                *el_mut(h, c)? += x * pc;
            }
        }
    }
    Ok(net)
}

/// The polynomial Bezier patch of a net (its weights 1) over `span x span`.
fn bezier(net: &Net, span: [f64; 2]) -> SolidResult<BSplineSurface> {
    let pts: Vec<Point3> = net
        .pts
        .iter()
        .map(|h| Point3::new(h[0], h[1], h[2]))
        .collect();
    let um = [net.p + 1, net.p + 1];
    let vm = [net.q + 1, net.q + 1];
    BSplineSurface::create_canonical(&CanonicalSurfaceData {
        u_degree: net.p,
        v_degree: net.q,
        control_points: &pts,
        form: BSplineSurfaceForm::Unspecified,
        u_multiplicities: &um,
        v_multiplicities: &vm,
        u_knots: &span,
        v_knots: &span,
        knot_type: KnotType::Unspecified,
        weights: None,
        uv_domain: None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A skew quadrilateral's fill is the bilinear map of its corners: at the domain's
    /// centre their mean, and along each side the side, at degree 2.
    #[test]
    fn a_skew_quadrilateral_is_filled_by_its_bilinear_map() {
        let xs = [
            Point3::new(1.5, 0.5, 8.5),
            Point3::new(1.0, 2.0, 8.0),
            Point3::new(-2.5, -0.5, 7.5),
            Point3::new(0.0, -2.0, 8.0),
        ];
        let fill = fill_polygon(&xs, 1e-9).unwrap();
        assert_eq!(fill.degree, 2);
        assert_eq!(fill.parts.len(), 4);
        let centre = (xs[0] + xs[1] + xs[2] + xs[3]) * 0.25;
        let at = fill.surface.evaluate_point(Point2::new(0.0, 0.0)).unwrap();
        assert!(at.distance(centre) < 1e-13, "{at:?}");
        // the bilinear map p = (x + y + 1) / 2, q = (y - x + 1) / 2 of the diamond
        // B_0 = (1, 0) .. B_3 = (0, -1): X_3 (1 - p)(1 - q) + X_0 p (1 - q) + X_1 p q
        // + X_2 (1 - p) q
        for (x, y) in [(0.25, 0.25), (-0.5, 0.1), (0.0, -0.75), (0.3, -0.6)] {
            let (p, q) = ((x + y + 1.0) / 2.0, (y - x + 1.0) / 2.0);
            let want = xs[3] * ((1.0 - p) * (1.0 - q))
                + xs[0] * (p * (1.0 - q))
                + xs[1] * (p * q)
                + xs[2] * ((1.0 - p) * q);
            let at = fill.surface.evaluate_point(Point2::new(x, y)).unwrap();
            assert!(at.distance(want) < 1e-13, "{x} {y}: {at:?} {want:?}");
        }
    }

    /// A planar regular hexagon's fill is its plane: every part's control points on it.
    #[test]
    fn a_planar_hexagon_is_filled_flat() {
        let xs: Vec<Point3> = (0..6)
            .map(|m| {
                let a = f64::from(m) * PI / 3.0;
                Point3::new(2.0 * a.cos(), 2.0 * a.sin(), 5.0)
            })
            .collect();
        let fill = fill_polygon(&xs, 1e-9).unwrap();
        for part in &fill.parts {
            for p in &part.control_point_net().unwrap().points {
                assert!((p.z - 5.0).abs() < 1e-12, "{p:?}");
            }
        }
    }

    #[test]
    fn a_polygon_needs_three_corners() {
        assert!(fill_polygon(&[Point3::ZERO, Point3::new(1.0, 0.0, 0.0)], 1e-9).is_err());
    }
}
