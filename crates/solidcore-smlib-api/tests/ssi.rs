//! Surface/surface intersection: the surface constructors, the three methods, the
//! curves' B-splines checked against both surfaces, touching points, and hostile
//! arguments.
//!
//! Source: tangency-labs/smlib `capi/tests/ssi.cpp` (Apache-2.0, Copyright SolidCraft
//! contributors) at fork commit c2f1788 (replay/1.0.1): `plane_cylinder`, `plane_sphere`,
//! `nurbs_paraboloids`, `sphere_touching_plane`, `partly_coincident_with_tangent_line`,
//! `sphere_cone_top_circle`, `sunk_cylinder`, `hostile`. The buffer-size checks have no
//! counterpart (the accessors return values); `SSI_BASIC` is `SMLIB_ERR_UNSUPPORTED` in
//! this port.

#![allow(clippy::unwrap_used, clippy::panic, clippy::indexing_slicing)] // test helpers

use solidcore_smlib_api::{
    Smlib, SsiNurbs, SsiResult, Status, Surface, SSI_ADVANCED, SSI_BASIC, SSI_CROSSING,
    SSI_CURVE_3D, SSI_CURVE_UV_A, SSI_CURVE_UV_B, SSI_GLOBAL, SSI_HAS_UV_A, SSI_HAS_UV_B,
    SSI_TOUCHING,
};

const X: [f64; 3] = [1.0, 0.0, 0.0];
const Y: [f64; 3] = [0.0, 1.0, 0.0];
const O: [f64; 3] = [0.0, 0.0, 0.0];

/// The curve's point at `t` (de Boor on the homogeneous points), `dim` coordinates.
fn point(c: &SsiNurbs, dim: usize, t: f64) -> Vec<f64> {
    let p = usize::try_from(c.degree).unwrap();
    let n = c.weights.len();
    let mut k = p;
    while k + 1 < n && c.knots[k + 1] <= t {
        k += 1;
    }
    let mut d = vec![0.0; (p + 1) * 4];
    for j in 0..=p {
        let q = k - p + j;
        for e in 0..dim {
            d[4 * j + e] = c.points[dim * q + e] * c.weights[q];
        }
        d[4 * j + 3] = c.weights[q];
    }
    for rr in 1..=p {
        for j in (rr..=p).rev() {
            let (lo, hi) = (c.knots[k - p + j], c.knots[k + 1 + j - rr]);
            let a = if hi > lo { (t - lo) / (hi - lo) } else { 0.0 };
            for e in 0..4 {
                d[4 * j + e] = (1.0 - a) * d[4 * (j - 1) + e] + a * d[4 * j + e];
            }
        }
    }
    (0..dim).map(|e| d[4 * p + e] / d[4 * p + 3]).collect()
}

/// `n` points spread evenly over the B-spline's knot range (none for a missing curve).
#[allow(clippy::cast_precision_loss)] // small counts
fn sample(c: Option<&SsiNurbs>, dim: usize, n: usize) -> Vec<Vec<f64>> {
    let Some(c) = c else {
        return vec![vec![f64::NAN; dim]; n];
    };
    let p = usize::try_from(c.degree).unwrap();
    let (t0, t1) = (c.knots[p], c.knots[c.weights.len()]);
    (0..n)
        .map(|k| point(c, dim, t0 + (t1 - t0) * k as f64 / (n - 1) as f64))
        .collect()
}

struct Samples {
    flags: i32,
    xyz: Vec<Vec<f64>>,
    uv_a: Vec<Vec<f64>>,
    uv_b: Vec<Vec<f64>>,
}

fn curve(r: &SsiResult, i: usize) -> Samples {
    let info = r.curve_info(i).unwrap();
    let c3 = r.curve_nurbs(i, SSI_CURVE_3D).unwrap();
    let ua = r.curve_nurbs(i, SSI_CURVE_UV_A);
    let ub = r.curve_nurbs(i, SSI_CURVE_UV_B);
    assert_eq!(ua.is_ok(), info.flags & SSI_HAS_UV_A != 0);
    assert_eq!(ub.is_ok(), info.flags & SSI_HAS_UV_B != 0);
    for e in [&ua, &ub] {
        if let Err(s) = e {
            assert_eq!(*s, Status::ERR_EMPTY);
        }
    }
    let n = 201;
    Samples {
        flags: info.flags,
        xyz: sample(Some(c3), 3, n),
        uv_a: sample(ua.ok(), 2, n),
        uv_b: sample(ub.ok(), 2, n),
    }
}

fn dist(p: &[f64], q: &[f64]) -> f64 {
    ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt()
}

/// Every sample of every curve: `on` (both surfaces' implicit equations) within `tol`, and
/// where the method returned UV curves the surfaces at those UVs within `tol` too.
fn check_on_both(
    api: &mut Smlib,
    r: &SsiResult,
    a: &Surface,
    b: &Surface,
    tol: f64,
    on: impl Fn(&[f64]) -> f64,
) {
    let mut worst: f64 = 0.0;
    for i in 0..r.curve_count() {
        let s = curve(r, i);
        for (k, p) in s.xyz.iter().enumerate() {
            assert!(p.iter().all(|x| x.is_finite()));
            worst = worst.max(on(p));
            if s.flags & SSI_HAS_UV_A != 0 {
                let q = api.surface_evaluate(a, s.uv_a[k][0], s.uv_a[k][1]).unwrap();
                worst = worst.max(dist(p, &q));
            } else {
                assert!(s.uv_a[k][0].is_nan());
            }
            if s.flags & SSI_HAS_UV_B != 0 {
                let q = api.surface_evaluate(b, s.uv_b[k][0], s.uv_b[k][1]).unwrap();
                worst = worst.max(dist(p, &q));
            }
        }
    }
    assert!(worst < tol, "worst {worst} >= {tol}");
}

fn type_of(r: &SsiResult, i: usize) -> i32 {
    r.curve_info(i).unwrap().curve_type
}

/// `c1_plane_cyl_ellipse`: a tilted plane through a cylinder's axis, one closed ellipse.
#[test]
fn plane_cylinder() {
    let mut api = Smlib::default();
    let py = [0.0, 0.894_427_190_999_915_9, 0.447_213_595_499_957_9];
    let plane = api
        .surface_plane(O, X, py, [-3.0, 3.0], [-3.0, 3.0])
        .unwrap();
    let cyl = api
        .surface_cylinder(O, X, Y, 1.0, [0.0, 360.0], [-2.0, 2.0])
        .unwrap();
    let p = api.surface_evaluate(&plane, 1.0, 2.0).unwrap();
    assert!(
        (p[0] - 1.0).abs() < 1e-12
            && (p[1] - 2.0 * py[1]).abs() < 1e-12
            && (p[2] - 2.0 * py[2]).abs() < 1e-12
    );
    assert_eq!(Smlib::surface_domain(&plane), [-3.0, 3.0, -3.0, 3.0]);
    // On the plane: n . p = 0 with n = X x py; on the cylinder: x^2 + y^2 = 1.
    let nrm = [0.0, -py[2], py[1]];
    let on = |q: &[f64]| {
        (q[0] * nrm[0] + q[1] * nrm[1] + q[2] * nrm[2])
            .abs()
            .max((q[0].hypot(q[1]) - 1.0).abs())
    };
    for method in [SSI_GLOBAL, SSI_ADVANCED] {
        let r = api.ssi(&plane, &cyl, method, 1e-6).unwrap();
        assert_eq!(r.curve_count(), 1);
        assert_eq!(type_of(&r, 0), SSI_CROSSING);
        assert_eq!(r.point_count(), 0);
        check_on_both(&mut api, &r, &plane, &cyl, 1e-5, on);
        let s = curve(&r, 0);
        if method == SSI_ADVANCED {
            assert_eq!(s.flags, SSI_HAS_UV_A | SSI_HAS_UV_B);
        } else {
            // The analytic special case returns no UV curves.
            assert_eq!(s.flags, 0);
        }
    }
}

/// `c1_plane_sphere_circle`: z = 1 through a sphere of radius 2, a circle of radius sqrt(3).
#[test]
fn plane_sphere() {
    let mut api = Smlib::default();
    let plane = api
        .surface_plane([0.0, 0.0, 1.0], X, Y, [-3.0, 3.0], [-3.0, 3.0])
        .unwrap();
    let sphere = api
        .surface_sphere(O, X, Y, 2.0, [0.0, 360.0], [-90.0, 90.0])
        .unwrap();
    let r = api.ssi(&plane, &sphere, SSI_GLOBAL, 0.0).unwrap();
    assert_eq!(r.curve_count(), 1);
    assert_eq!(type_of(&r, 0), SSI_CROSSING);
    check_on_both(&mut api, &r, &plane, &sphere, 1e-6, |q| {
        (q[2] - 1.0).abs().max((dist(q, &O) - 2.0).abs())
    });
}

/// `c1n_parab_parab_circle`: z = x^2 + y^2 against z = 2 - x^2 - y^2 over [-2, 2]^2, both as
/// biquadratic Beziers; they cross in the unit circle at z = 1.
#[test]
fn nurbs_paraboloids() {
    let mut api = Smlib::default();
    let k = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
    let up = [
        -2.0, -2.0, 8.0, -2.0, 0.0, 0.0, -2.0, 2.0, 8.0, 0.0, -2.0, 0.0, 0.0, 0.0, -8.0, 0.0, 2.0,
        0.0, 2.0, -2.0, 8.0, 2.0, 0.0, 0.0, 2.0, 2.0, 8.0,
    ];
    let down = [
        -2.0, -2.0, -6.0, -2.0, 0.0, 2.0, -2.0, 2.0, -6.0, 0.0, -2.0, 2.0, 0.0, 0.0, 10.0, 0.0,
        2.0, 2.0, 2.0, -2.0, -6.0, 2.0, 0.0, 2.0, 2.0, 2.0, -6.0,
    ];
    let w = [1.0; 9];
    let a = api.surface_nurbs(2, 2, 3, 3, &k, &k, &up, None).unwrap();
    let b = api
        .surface_nurbs(2, 2, 3, 3, &k, &k, &down, Some(&w))
        .unwrap();
    let p = api.surface_evaluate(&a, 0.5, 0.5).unwrap();
    assert!(dist(&p, &O) < 1e-12);
    for method in [SSI_GLOBAL, SSI_ADVANCED] {
        let r = api.ssi(&a, &b, method, 1e-6).unwrap();
        assert_eq!(r.curve_count(), 1);
        assert_eq!(type_of(&r, 0), SSI_CROSSING);
        check_on_both(&mut api, &r, &a, &b, 1e-5, |q| {
            (q[2] - 1.0).abs().max((q[0].hypot(q[1]) - 1.0).abs())
        });
        assert_eq!(curve(&r, 0).flags, SSI_HAS_UV_A | SSI_HAS_UV_B);
    }
    // the basic intersector's own FindInteriorCurves is not in the port
    assert_eq!(
        api.ssi(&a, &b, SSI_BASIC, 1e-6).unwrap_err(),
        Status::ERR_UNSUPPORTED
    );
}

/// `c3_sphere_on_plane`: a unit sphere resting on z = 0 touches it at the origin.
#[test]
fn sphere_touching_plane() {
    let mut api = Smlib::default();
    let sphere = api
        .surface_sphere([0.0, 0.0, 1.0], X, Y, 1.0, [0.0, 360.0], [-90.0, 90.0])
        .unwrap();
    let plane = api
        .surface_plane(O, X, Y, [-3.0, 3.0], [-3.0, 3.0])
        .unwrap();
    let r = api.ssi(&sphere, &plane, SSI_GLOBAL, 1e-6).unwrap();
    assert_eq!(r.curve_count(), 1);
    assert_eq!(type_of(&r, 0), SSI_TOUCHING);
    assert_eq!(r.point_count(), 1);
    assert!(dist(&r.point(0).unwrap().xyz, &O) < 1e-9);
    assert_eq!(r.point(1).unwrap_err(), Status::ERR_ARG);
}

/// A plane z = 0 against z = f(x) over x in [0, 4], y in [0, 1], with f = 0 for x <= 1 and
/// f = s (x - 1)^2 (x - 3)^2 beyond: the surfaces coincide over x in [0, 1] and touch along
/// the line x = 3 without crossing. `mirror` maps x to 4 - x. The tangent line must come
/// out whichever contact the tangent-curve search meets first. (The fork's tangent-curve
/// tracer: `ssi_tangent_curve_tracer` here.)
fn partly_coincident_with_tangent_line(mirror: bool) {
    let mut api = Smlib::default();
    api.kernel_mut().fixes.ssi_tangent_curve_tracer = true;
    let s = 0.25;
    let ku = [
        0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4.0, 4.0, 4.0, 4.0, 4.0,
    ];
    let kv = [0.0, 0.0, 1.0, 1.0];
    let gx = [0.0, 0.25, 0.5, 0.75, 1.75, 2.5, 3.25, 4.0];
    let gz = [0.0, 0.0, 0.0, 0.0, 0.0, 6.0 * s, -9.0 * s, 9.0 * s];
    let mut pts = Vec::new();
    for i in 0..8 {
        for j in 0..2 {
            pts.extend([
                if mirror { 4.0 - gx[i] } else { gx[i] },
                f64::from(j),
                gz[i],
            ]);
        }
    }
    let xt = if mirror { 1.0 } else { 3.0 };
    let plane = api
        .surface_plane(O, X, Y, [-1.0, 5.0], [-1.0, 2.0])
        .unwrap();
    let bump = api.surface_nurbs(4, 1, 8, 2, &ku, &kv, &pts, None).unwrap();
    let q = api.surface_evaluate(&bump, 2.0, 0.5).unwrap();
    assert!((q[0] - 2.0).abs() < 1e-12 && (q[2] - s).abs() < 1e-12);
    let q = api.surface_evaluate(&bump, 3.0, 0.5).unwrap();
    assert!((q[0] - xt).abs() < 1e-12 && q[2].abs() < 1e-12);
    let q = api.surface_evaluate(&bump, 0.5, 0.5).unwrap();
    assert!(q[2].abs() < 1e-15);
    for method in [SSI_GLOBAL, SSI_ADVANCED] {
        let r = api.ssi(&plane, &bump, method, 1e-6).unwrap();
        // The tangent line: a curve all along x = xt, z = 0, from y = 0 to y = 1.
        let found = (0..r.curve_count()).any(|i| {
            let c = curve(&r, i);
            let off = c
                .xyz
                .iter()
                .map(|p| (p[0] - xt).abs().max(p[2].abs()))
                .fold(0.0, f64::max);
            let ylo = c.xyz.iter().map(|p| p[1]).fold(f64::INFINITY, f64::min);
            let yhi = c.xyz.iter().map(|p| p[1]).fold(f64::NEG_INFINITY, f64::max);
            off < 1e-4 && ylo < 1e-4 && yhi > 1.0 - 1e-4
        });
        assert!(
            found,
            "mirror={mirror}, method {method}: no curve along x = {xt}"
        );
    }
}

#[test]
fn partly_coincident_with_tangent_line_both_ends() {
    partly_coincident_with_tangent_line(false);
    partly_coincident_with_tangent_line(true);
}

/// A sphere centred on a cone's axis crosses the cone's straight line just beyond its top
/// end (by 2 tol along the line): that crossing is the cone's top circle, on both surfaces
/// to within tol; the line's other crossing gives a second circle lower down.
#[test]
fn sphere_cone_top_circle() {
    let mut api = Smlib::default();
    let (tol, half) = (1e-3, 70.0_f64);
    let rad = half.to_radians();
    let (ta, sa, ca) = (rad.tan(), rad.sin(), rad.cos());
    let c = [0.0, 0.0, 6.0];
    let p = [ta + 2.0 * tol * sa, 0.0, 1.0 + 2.0 * tol * ca];
    let radius = dist(&p, &c);
    let sphere = api
        .surface_sphere(c, X, Y, radius, [0.0, 360.0], [-90.0, 90.0])
        .unwrap();
    let cone = api.surface_cone(O, X, Y, half, [0.0, 1.0]).unwrap();
    assert!((dist(&[ta, 0.0, 1.0], &c) - radius).abs() < tol);
    let r = api.ssi(&sphere, &cone, SSI_GLOBAL, tol).unwrap();
    assert_eq!(r.curve_count(), 2);
    let mut at_top = false;
    for i in 0..r.curve_count() {
        assert_eq!(type_of(&r, i), SSI_CROSSING);
        let s = curve(&r, i);
        let zlo = s.xyz.iter().map(|p| p[2]).fold(f64::INFINITY, f64::min);
        let zhi = s.xyz.iter().map(|p| p[2]).fold(f64::NEG_INFINITY, f64::max);
        at_top = at_top || ((zlo - 1.0).abs() < tol && (zhi - 1.0).abs() < tol);
    }
    assert!(at_top);
    check_on_both(&mut api, &r, &sphere, &cone, 2.0 * tol, |q| {
        (dist(q, &c) - radius)
            .abs()
            .max((q[0].hypot(q[1]) * ca - q[2] * sa).abs())
    });
}

/// `c7_sunk_cyl_*`: a unit cylinder along -x over x in [-2, 2], its axis at z = cos a, sunk
/// 1 - cos a (5e-7 for a = 1e-3, 5e-13 for 1e-6) into the plane z = 0, which it crosses at
/// the angle a in the two lines y = +-sin a. At approximation tolerance 1e-6 both lines come
/// out, each where it is: one line in the middle would be off them by sin a, up to a
/// thousand times the tolerance. (Bare surfaces bound no face, so no zone tolerance makes
/// the plane tangent.)
#[test]
fn sunk_cylinder() {
    let mut api = Smlib::default();
    let (u, angle, height, z) = ([-3.0, 3.0], [0.0, 360.0], [0.0, 4.0], [0.0, 0.0, 1.0]);
    for a in [1e-3, 1e-4, 1e-5, 1e-6_f64] {
        let (axis, half) = ([2.0, 0.0, a.cos()], a.sin());
        let plane = api.surface_plane(O, X, Y, u, u).unwrap();
        let cyl = api
            .surface_cylinder(axis, z, Y, 1.0, angle, height)
            .unwrap();
        let r = api.ssi(&plane, &cyl, SSI_GLOBAL, 1e-6).unwrap();
        // A line along x from -2 to 2 at z = 0 on each side, within a tenth of sin a of
        // y = +-sin a.
        let mut sides = [0; 2];
        for i in 0..r.curve_count() {
            let s = curve(&r, i);
            let (mut xlo, mut xhi, mut ylo, mut yhi, mut zoff) =
                (1e300_f64, -1e300_f64, 1e300_f64, -1e300_f64, 0.0_f64);
            for p in &s.xyz {
                xlo = xlo.min(p[0]);
                xhi = xhi.max(p[0]);
                ylo = ylo.min(p[1]);
                yhi = yhi.max(p[1]);
                zoff = zoff.max(p[2].abs());
            }
            let whole = xlo < -2.0 + 1e-6 && xhi > 2.0 - 1e-6 && zoff < 1e-9;
            for (side, n) in sides.iter_mut().enumerate() {
                let y = if side == 0 { -half } else { half };
                let near = (ylo - y).abs() < 0.1 * half && (yhi - y).abs() < 0.1 * half;
                *n += i32::from(whole && near);
            }
        }
        assert!(
            r.curve_count() == 2 && sides == [1, 1],
            "sunk_cylinder(a={a}): {} curve(s), {} at y = -sin a, {} at y = +sin a",
            r.curve_count(),
            sides[0],
            sides[1]
        );
    }
}

#[test]
fn hostile() {
    let mut api = Smlib::default();
    let arg = Err(Status::ERR_ARG);
    let u = [-1.0, 1.0];
    let (nan3, inf3, zero) = ([f64::NAN, 0.0, 0.0], [0.0, f64::INFINITY, 0.0], [0.0; 3]);
    let huge3 = [1e300, 0.0, 0.0];
    let (ang, wide, h, rev_h) = ([0.0, 360.0], [0.0, 400.0], [0.0, 1.0], [1.0, 0.0]);
    let (lat_bad, neg_h) = ([-100.0, 90.0], [-1.0, 1.0]);
    assert_eq!(api.surface_plane(nan3, X, Y, u, u).map(|_| ()), arg);
    assert_eq!(api.surface_plane(O, inf3, Y, u, u).map(|_| ()), arg);
    assert_eq!(api.surface_plane(O, zero, Y, u, u).map(|_| ()), arg);
    assert_eq!(api.surface_plane(O, X, X, u, u).map(|_| ()), arg);
    assert_eq!(api.surface_plane(huge3, X, Y, u, u).map(|_| ()), arg);
    assert_eq!(api.surface_plane(O, X, Y, [1.0, -1.0], u).map(|_| ()), arg);
    assert_eq!(api.surface_cylinder(O, X, Y, 0.0, ang, h).map(|_| ()), arg);
    assert_eq!(api.surface_cylinder(O, X, Y, -1.0, ang, h).map(|_| ()), arg);
    assert_eq!(
        api.surface_cylinder(O, X, Y, f64::NAN, ang, h).map(|_| ()),
        arg
    );
    assert_eq!(api.surface_cylinder(O, X, Y, 1.0, wide, h).map(|_| ()), arg);
    assert_eq!(
        api.surface_cylinder(O, X, Y, 1.0, ang, rev_h).map(|_| ()),
        arg
    );
    assert_eq!(api.surface_cone(O, X, Y, 0.0, h).map(|_| ()), arg);
    assert_eq!(api.surface_cone(O, X, Y, 90.0, h).map(|_| ()), arg);
    assert_eq!(api.surface_cone(O, X, Y, 30.0, neg_h).map(|_| ()), arg);
    assert_eq!(
        api.surface_sphere(O, X, Y, 1.0, ang, lat_bad).map(|_| ()),
        arg
    );
    assert_eq!(api.surface_sphere(O, X, zero, 1.0, ang, h).map(|_| ()), arg);
    assert_eq!(
        api.surface_torus(O, X, Y, 1.0, 2.0, ang, ang).map(|_| ()),
        arg
    );
    assert_eq!(
        api.surface_torus(O, X, Y, 3.0, 0.0, ang, ang).map(|_| ()),
        arg
    );

    let k = [0.0, 0.0, 1.0, 1.0];
    let (k_desc, k_flat, k_nan) = ([0.0, 1.0, 0.0, 1.0], [0.0; 4], [0.0, 0.0, f64::NAN, 1.0]);
    let pts = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0];
    let mut pts_nan = pts;
    pts_nan[4] = f64::NAN;
    let (w, w_zero, w_neg) = ([1.0; 4], [1.0, 0.0, 1.0, 1.0], [1.0, -1.0, 1.0, 1.0]);
    assert!(api
        .surface_nurbs(1, 1, 2, 2, &k, &k, &pts, Some(&w))
        .is_ok());
    assert_eq!(
        api.surface_nurbs(1, 1, 2, 2, &k_desc, &k, &pts, Some(&w))
            .map(|_| ()),
        arg
    );
    assert_eq!(
        api.surface_nurbs(1, 1, 2, 2, &k_flat, &k, &pts, Some(&w))
            .map(|_| ()),
        arg
    );
    assert_eq!(
        api.surface_nurbs(1, 1, 2, 2, &k_nan, &k, &pts, Some(&w))
            .map(|_| ()),
        arg
    );
    assert_eq!(
        api.surface_nurbs(1, 1, 2, 2, &k[..3], &k, &pts, Some(&w))
            .map(|_| ()),
        arg
    );
    assert_eq!(
        api.surface_nurbs(1, 1, 2, 2, &k, &k, &pts_nan, Some(&w))
            .map(|_| ()),
        arg
    );
    assert_eq!(
        api.surface_nurbs(1, 1, 2, 2, &k, &k, &pts, Some(&w_zero))
            .map(|_| ()),
        arg
    );
    assert_eq!(
        api.surface_nurbs(1, 1, 2, 2, &k, &k, &pts, Some(&w_neg))
            .map(|_| ()),
        arg
    );
    assert_eq!(
        api.surface_nurbs(0, 1, 2, 2, &k, &k, &pts, Some(&w))
            .map(|_| ()),
        arg
    );
    assert_eq!(
        api.surface_nurbs(2, 1, 2, 2, &k, &k, &pts, Some(&w))
            .map(|_| ()),
        arg
    );
    // Huge counts are rejected before any array is read.
    assert_eq!(
        api.surface_nurbs(1, 1, 2_000_000, 2, &k, &k, &pts, Some(&w))
            .map(|_| ()),
        arg
    );
    assert_eq!(
        api.surface_nurbs(1, 1, 100_000, 100_000, &k, &k, &pts, Some(&w))
            .map(|_| ()),
        arg
    );
    assert_eq!(
        api.surface_nurbs(40, 1, 41, 2, &k, &k, &pts, Some(&w))
            .map(|_| ()),
        arg
    );

    let a = api.surface_plane(O, X, Y, u, u).unwrap();
    assert_eq!(api.ssi(&a, &a, 3, 0.0).map(|_| ()), arg);
    assert_eq!(api.ssi(&a, &a, -1, 0.0).map(|_| ()), arg);
    assert_eq!(api.ssi(&a, &a, 0, f64::NAN).map(|_| ()), arg);
    assert_eq!(api.surface_evaluate(&a, f64::NAN, 0.0).map(|_| ()), arg);
}
