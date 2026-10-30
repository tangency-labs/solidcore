//! Surface/surface intersection: plain surface descriptions (built as the SSI corpus
//! runner builds them from a case file), SMLib's three SSI entry points, and the curves
//! they return, as SMLib holds them.
//!
//! Source: tangency-labs/smlib `capi/src/smlib_ssi.cpp` (Apache-2.0, Copyright SolidCraft
//! contributors) at fork commit ef62892 (replay/1.0.1).
//!
//! | Rust | C API (SMLib call) |
//! |------|-------|
//! | [`Surface`] | `smlib_surface` |
//! | [`Smlib::surface_plane`] | `smlib_surface_plane` (`SmPlane`) |
//! | [`Smlib::surface_cylinder`] | `smlib_surface_cylinder` (`SmCylinder`) |
//! | [`Smlib::surface_cone`] | `smlib_surface_cone` (`SmCone::CreateCanonical`) |
//! | [`Smlib::surface_sphere`] | `smlib_surface_sphere` (`SmSphere`) |
//! | [`Smlib::surface_torus`] | `smlib_surface_torus` (`SmTorus`) |
//! | [`Smlib::surface_nurbs`] | `smlib_surface_nurbs` (`SmBSplineSurface`) |
//! | [`Smlib::surface_domain`] | `smlib_surface_domain` |
//! | [`Smlib::surface_evaluate`] | `smlib_surface_evaluate` (`EvaluatePoint`) |
//! | [`Smlib::ssi`] | `smlib_ssi` (`GlobalSurfaceIntersect`, `SmAdvSurfaceIntersector`) |
//! | [`SsiResult`] and its methods | `smlib_ssi_result`, `smlib_ssi_curve_count`, `smlib_ssi_curve_info`, `smlib_ssi_curve_nurbs`, `smlib_ssi_point_count`, `smlib_ssi_point` |
//! | `nurbs`, `start`, `bounded`, `axes`, `range`, `length`, `knots`, `clamped`, `unique_knots` (private) | the same |
//!
//! Deviations: no `smlib_surface_free` / `smlib_ssi_result_free` (values are dropped); the
//! accessors return values instead of filling caller buffers, so `SMLIB_ERR_BUFFER` does
//! not arise; no case-file or surface dumps (the call's line is still logged).
//! [`SSI_BASIC`] (`SmSurfaceIntersector::DoIntersection` with the basic intersector's
//! `FindInteriorCurves`) is `SMLIB_ERR_UNSUPPORTED`: the port runs `DoIntersection` only
//! through the advanced intersector. A curve whose B-spline cannot make its NURB is
//! returned as no curve, as the C API returns one whose NURB pointer is null.

use solidcore::curve::curve::{Curve as KCurve, CurveKind};
use solidcore::geom::axis2_placement::Axis2Placement;
use solidcore::geom::extent2d::Extent2d;
use solidcore::geom::tol::approx_tol3d;
use solidcore::geom::vector2d::{Point2, Vec2};
use solidcore::geom::vector3d::{Point3, Vec3};
use solidcore::intersect::surface_intersector::p5::do_intersection;
use solidcore::nurbs::b_spline_curve::{BSplineCurve, KnotType};
use solidcore::nurbs::b_spline_surface::{
    BSplineSurface, BSplineSurfaceForm, CanonicalSurfaceData,
};
use solidcore::nurbs::math::{CPoint, NOW};
use solidcore::nurbs::srf_basic::Surface as NurbSurface;
use solidcore::surface::cone::Cone;
use solidcore::surface::plane::Plane;
use solidcore::surface::sphere::{Sphere, SphereOptions};
use solidcore::surface::surface::p2::{SurfaceIntersection, TsectCurveType};
use solidcore::surface::surface::Surface as KSurface;
use solidcore::surface::torus::{Torus, TorusOptions};
use solidcore::Kernel;

use crate::session::{finite, num, vec, Fail, Smlib};
use crate::status::{Result, Status};

/// `SMLIB_SURFACE_MAX_COORD`: the largest coordinate, length or parameter a surface
/// description may use.
pub const SURFACE_MAX_COORD: f64 = 1e8;
/// `SMLIB_SURFACE_MAX_DEGREE`.
pub const SURFACE_MAX_DEGREE: i32 = 25;
/// `SMLIB_SURFACE_MAX_POLES`: the most control points in either direction and in all.
pub const SURFACE_MAX_POLES: i32 = 1_000_000;

/// `SMLIB_SSI_GLOBAL`: `SmSurface::GlobalSurfaceIntersect` (analytic special cases, else
/// [`SSI_ADVANCED`]).
pub const SSI_GLOBAL: i32 = 0;
/// `SMLIB_SSI_ADVANCED`: `SmAdvSurfaceIntersector` (the general NURBS path).
pub const SSI_ADVANCED: i32 = 1;
/// `SMLIB_SSI_BASIC`: `SmSurfaceIntersector` (the basic intersector; not in this port).
pub const SSI_BASIC: i32 = 2;

/// `SMLIB_SSI_NONE`: the method reported no type for this curve.
pub const SSI_NONE: i32 = -1;
/// `SMLIB_SSI_UNKNOWN`.
pub const SSI_UNKNOWN: i32 = 0;
/// `SMLIB_SSI_TOUCHING`: a degenerate curve, one point where the surfaces touch.
pub const SSI_TOUCHING: i32 = 1;
/// `SMLIB_SSI_CROSSING`: a transversal curve.
pub const SSI_CROSSING: i32 = 2;
/// `SMLIB_SSI_TANGENT`: a curve along which the surfaces are tangent.
pub const SSI_TANGENT: i32 = 3;
/// `SMLIB_SSI_COINCIDENT`: where the surfaces butt up against each other.
pub const SSI_COINCIDENT: i32 = 4;
/// `SMLIB_SSI_NEAR_TANGENT`: a curve with a small angle of intersection.
pub const SSI_NEAR_TANGENT: i32 = 5;
/// `SMLIB_SSI_REGION_BOUNDARY`: bounds a region where the surfaces coincide.
pub const SSI_REGION_BOUNDARY: i32 = 6;

/// `SMLIB_SSI_HAS_UV_A`: the method returned the curve's UV curve on `a`.
pub const SSI_HAS_UV_A: i32 = 1;
/// `SMLIB_SSI_HAS_UV_B`: the method returned the curve's UV curve on `b`.
pub const SSI_HAS_UV_B: i32 = 2;

/// `SMLIB_SSI_CURVE_3D`.
pub const SSI_CURVE_3D: i32 = 0;
/// `SMLIB_SSI_CURVE_UV_A`: on `a`, in `a`'s UV.
pub const SSI_CURVE_UV_A: i32 = 1;
/// `SMLIB_SSI_CURVE_UV_B`: on `b`, in `b`'s UV.
pub const SSI_CURVE_UV_B: i32 = 2;

/// `smlib_surface`: a surface with the UV domain [`Smlib::ssi`] intersects within.
#[derive(Debug, Clone)]
pub struct Surface {
    surface: KSurface,
    domain: Extent2d,
}

impl Surface {
    /// The kernel surface.
    #[must_use]
    pub const fn surface(&self) -> &KSurface {
        &self.surface
    }

    /// The UV domain.
    pub(crate) const fn domain(&self) -> &Extent2d {
        &self.domain
    }
}

/// A curve as SMLib holds it: its (rational) B-spline. `dim` 0: none (the method
/// returned no such curve, or SMLib's curve is not a B-spline).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SsiNurbs {
    /// 3 for a 3D curve, 2 for a UV curve, 0 for none.
    pub dim: i32,
    /// The degree.
    pub degree: i32,
    /// All the knots.
    pub knots: Vec<f64>,
    /// The Cartesian control points, `dim` numbers each.
    pub points: Vec<f64>,
    /// One weight per control point.
    pub weights: Vec<f64>,
}

/// What [`SsiResult::curve_info`] gives (`smlib_ssi_curve_info`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SsiCurveInfo {
    /// Its `SSI_*` type ([`SSI_NONE`] when the method gave none).
    pub curve_type: i32,
    /// Its `SSI_HAS_UV_*` flags.
    pub flags: i32,
    /// SMLib's deviation estimate (NaN when the method gave none).
    pub deviation: f64,
    /// The 3D curve's natural interval.
    pub interval: [f64; 2],
}

#[derive(Debug, Clone)]
struct SsiCurve {
    info: SsiCurveInfo,
    c3: SsiNurbs,
    uv_a: SsiNurbs,
    uv_b: SsiNurbs,
}

/// A touching point (`smlib_ssi_point`): where on both surfaces, NaN where not known.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SsiPoint {
    /// The point.
    pub xyz: [f64; 3],
    /// Its UV on `a`.
    pub uv_a: [f64; 2],
    /// Its UV on `b`.
    pub uv_b: [f64; 2],
}

/// `smlib_ssi_result`: the curves and touching points of one intersection.
#[derive(Debug, Clone, Default)]
pub struct SsiResult {
    curves: Vec<SsiCurve>,
    points: Vec<SsiPoint>,
}

impl SsiResult {
    /// `smlib_ssi_curve_count`.
    #[must_use]
    pub fn curve_count(&self) -> usize {
        self.curves.len()
    }

    /// `smlib_ssi_curve_info`: curve `i`'s type, flags, deviation and interval.
    ///
    /// # Errors
    ///
    /// `SMLIB_ERR_ARG` for no curve `i`.
    pub fn curve_info(&self, i: usize) -> Result<SsiCurveInfo> {
        self.curves.get(i).map(|k| k.info).ok_or(Status::ERR_ARG)
    }

    /// `smlib_ssi_curve_nurbs`: curve `i` (`which`: [`SSI_CURVE_3D`], [`SSI_CURVE_UV_A`]
    /// or [`SSI_CURVE_UV_B`]) as SMLib's B-spline.
    ///
    /// # Errors
    ///
    /// `SMLIB_ERR_ARG` for no curve `i` or another `which`; `SMLIB_ERR_EMPTY` for a UV
    /// curve the method did not return; `SMLIB_ERR_UNSUPPORTED` for a curve it returned
    /// that is not a B-spline.
    pub fn curve_nurbs(&self, i: usize, which: i32) -> Result<&SsiNurbs> {
        let k = self.curves.get(i).ok_or(Status::ERR_ARG)?;
        let c = match which {
            SSI_CURVE_3D => &k.c3,
            SSI_CURVE_UV_A => &k.uv_a,
            SSI_CURVE_UV_B => &k.uv_b,
            _ => return Err(Status::ERR_ARG),
        };
        if c.dim == 0 {
            let returned = which == SSI_CURVE_3D
                || k.info.flags
                    & (if which == SSI_CURVE_UV_A {
                        SSI_HAS_UV_A
                    } else {
                        SSI_HAS_UV_B
                    })
                    != 0;
            return Err(if returned {
                Status::ERR_UNSUPPORTED
            } else {
                Status::ERR_EMPTY
            });
        }
        Ok(c)
    }

    /// `smlib_ssi_point_count`: the touching points (one per [`SSI_TOUCHING`] curve).
    #[must_use]
    pub fn point_count(&self) -> usize {
        self.points.len()
    }

    /// `smlib_ssi_point`: touching point `i`.
    ///
    /// # Errors
    ///
    /// `SMLIB_ERR_ARG` for no point `i`.
    pub fn point(&self, i: usize) -> Result<SsiPoint> {
        self.points.get(i).copied().ok_or(Status::ERR_ARG)
    }
}

/// `v` finite and every value within `limit`.
fn bounded(v: &[f64], limit: f64) -> bool {
    finite(v) && v.iter().all(|x| x.abs() <= limit)
}

/// Two finite, non-zero axes that are not (nearly) parallel.
fn axes(x: &[f64; 3], y: &[f64; 3]) -> bool {
    if !bounded(x, SURFACE_MAX_COORD) || !bounded(y, SURFACE_MAX_COORD) {
        return false;
    }
    let lx = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
    let ly = (y[0] * y[0] + y[1] * y[1] + y[2] * y[2]).sqrt();
    let positive = lx > 1e-12 && ly > 1e-12;
    if !positive {
        return false;
    }
    let c = [
        x[1] * y[2] - x[2] * y[1],
        x[2] * y[0] - x[0] * y[2],
        x[0] * y[1] - x[1] * y[0],
    ];
    (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt() > 1e-9 * lx * ly
}

/// An increasing interval `[r[0], r[1]]` inside `[lo, hi]` spanning at most `span`.
fn range(r: &[f64; 2], lo: f64, hi: f64, span: f64) -> bool {
    bounded(r, SURFACE_MAX_COORD) && r[0] < r[1] && r[0] >= lo && r[1] <= hi && r[1] - r[0] <= span
}

fn length(x: f64) -> bool {
    x.is_finite() && x > 1e-9 && x <= SURFACE_MAX_COORD
}

fn v3(v: &[f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}

/// The runner's `unit_cross`.
fn unit_cross(a: Vec3, b: Vec3) -> Vec3 {
    let c = Vec3::new(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    );
    let l = (c.x * c.x + c.y * c.y + c.z * c.z).sqrt();
    Vec3::new(c.x / l, c.y / l, c.z / l)
}

/// `construct`'s tail: the surface's natural UV domain, which must be finite and not
/// empty (`SM_ERR_INVALID_INPUT`).
fn finish(k: &Kernel, surface: KSurface) -> core::result::Result<Surface, Fail> {
    let domain = surface.natural_uv_domain_with_nurb(k)?;
    let d = [
        domain.min().x,
        domain.max().x,
        domain.min().y,
        domain.max().y,
    ];
    let ordered = d[0] < d[1] && d[2] < d[3];
    if !finite(&d) || !ordered {
        return Err(Fail::Status(Status::SM_ERR_INVALID_INPUT));
    }
    Ok(Surface { surface, domain })
}

/// A knot vector of `n + p + 1` finite, non-decreasing values, with a non-empty domain
/// `[k[p], k[n]]`, no interior knot repeated more than `p` times and no end knot more
/// than `p + 1`.
fn knots(k: &[f64], n: usize, p: usize) -> bool {
    if k.len() != n + p + 1 || !bounded(k, SURFACE_MAX_COORD) {
        return false;
    }
    let mut run = 1;
    for w in k.windows(2) {
        let &[a, b] = w else {
            return false;
        };
        if b < a {
            return false;
        }
        run = if b == a { run + 1 } else { 1 };
        if run > p + 1 {
            return false;
        }
    }
    let (Some(&lo), Some(&hi)) = (k.get(p), k.get(n)) else {
        return false;
    };
    // Interior multiplicity: a knot strictly inside (k[p], k[n]) at most p times.
    run = 1;
    for w in k.windows(2) {
        let &[a, b] = w else {
            return false;
        };
        run = if b == a { run + 1 } else { 1 };
        if run > p && b > lo && b < hi {
            return false;
        }
    }
    lo < hi
}

/// A knot vector [`knots`] accepted is clamped: its first and its last value each `p + 1`
/// times.
fn clamped(k: &[f64], p: usize) -> bool {
    let n = k.len();
    k.first() == k.get(p) && k.last() == n.checked_sub(1 + p).and_then(|i| k.get(i))
}

/// Unique values and their multiplicities, the form `SmBSplineSurface::CreateCanonical`
/// takes.
fn unique_knots(k: &[f64]) -> (Vec<f64>, Vec<usize>) {
    let mut values: Vec<f64> = Vec::new();
    let mut mults: Vec<usize> = Vec::new();
    for &x in k {
        match (values.last(), mults.last_mut()) {
            (Some(&last), Some(m)) if last == x => *m += 1,
            _ => {
                values.push(x);
                mults.push(1);
            }
        }
    }
    (values, mults)
}

/// The `SmBSplineCurve` (base) of a curve kind derived from it, `None` for the others.
const fn b_spline_base(c: &KCurve) -> Option<&BSplineCurve> {
    match c.kind() {
        CurveKind::BSpline(b) => Some(b),
        CurveKind::Line(l) => Some(l.base()),
        CurveKind::Ellipse(e) => Some(e.base()),
        CurveKind::Hyperbola(h) => Some(h.base()),
        CurveKind::Parabola(p) => Some(p.base()),
        _ => None,
    }
}

/// `c` as its B-spline, exactly: SMLib's lines, circles and ellipses are B-spline curves
/// too, and its intersection curves B-splines it fitted within the approximation
/// tolerance. Nothing (dim 0) for a curve that is not a B-spline: it is never
/// approximated here.
fn nurbs(c: Option<&KCurve>, dim: usize) -> SsiNurbs {
    let mut out = SsiNurbs::default();
    let Some(base) = c.and_then(b_spline_base) else {
        return out;
    };
    // GetOrCreateGwNurbPointer through a const cast: made on a copy here
    let mut b = base.clone_exact();
    let Ok(n) = b.get_or_create_gw_nurb() else {
        return out;
    };
    let (Ok(d), Ok(degree)) = (i32::try_from(dim), i32::try_from(n.degree())) else {
        return out;
    };
    out.dim = d;
    out.degree = degree;
    out.knots = n.knots().to_vec();
    for p in n.cpts() {
        let w = if p.w == NOW { 1.0 } else { p.w };
        let q = [p.x / w, p.y / w, p.z / w];
        out.points.extend(q.iter().take(dim));
        out.weights.push(w);
    }
    out
}

/// `c`'s point at the start of its natural interval, NaN where there is no curve or it
/// cannot be evaluated.
fn start<const D: usize>(k: &Kernel, c: Option<&KCurve>) -> [f64; D] {
    let mut out = [f64::NAN; D];
    if let Some(c) = c {
        if let Ok(p) = c.evaluate_point(k, c.natural_interval().min()) {
            for (o, v) in out.iter_mut().zip([p.x, p.y, p.z]) {
                *o = v;
            }
        }
    }
    out
}

impl Smlib {
    /// `smlib_surface_plane`: the plane through `origin` spanned by `xaxis` and `yaxis`
    /// over `u` x `v` (an `SmPlane` with unit parameter scale).
    ///
    /// # Errors
    ///
    /// `SMLIB_ERR_ARG` for a coordinate beyond [`SURFACE_MAX_COORD`], zero or parallel axes,
    /// or an empty or reversed range; a kernel status otherwise.
    pub fn surface_plane(
        &mut self,
        origin: [f64; 3],
        xaxis: [f64; 3],
        yaxis: [f64; 3],
        u: [f64; 2],
        v: [f64; 2],
    ) -> Result<Surface> {
        let m = SURFACE_MAX_COORD;
        if !bounded(&origin, m)
            || !axes(&xaxis, &yaxis)
            || !range(&u, -m, m, 2.0 * m)
            || !range(&v, -m, m, 2.0 * m)
        {
            return Err(Status::ERR_ARG);
        }
        self.run(
            || {
                format!(
                    "smlib_surface_plane(origin={}, xaxis={}, yaxis={}, u={}, v={})",
                    vec(&origin),
                    vec(&xaxis),
                    vec(&yaxis),
                    vec(&u),
                    vec(&v)
                )
            },
            &[],
            |api| {
                let dom = Extent2d::from_bounds(u[0], v[0], u[1], v[1]);
                let plane = Plane::new(
                    v3(&origin),
                    v3(&xaxis),
                    v3(&yaxis),
                    Vec2::new(1.0, 1.0),
                    &dom,
                );
                finish(api.kernel(), KSurface::from(plane))
            },
        )
    }

    /// `smlib_surface_cylinder`: the cylinder of `radius` about the axis `xaxis` x `yaxis`
    /// through `origin`, from `angle_deg[0]` to `angle_deg[1]` round it and from
    /// `height[0]` to `height[1]` along it.
    ///
    /// # Errors
    ///
    /// `SMLIB_ERR_ARG` for a bad coordinate, axes, radius or range (an angle range within
    /// -360..360 and at most 360); a kernel status otherwise.
    pub fn surface_cylinder(
        &mut self,
        origin: [f64; 3],
        xaxis: [f64; 3],
        yaxis: [f64; 3],
        radius: f64,
        angle_deg: [f64; 2],
        height: [f64; 2],
    ) -> Result<Surface> {
        let m = SURFACE_MAX_COORD;
        if !bounded(&origin, m)
            || !axes(&xaxis, &yaxis)
            || !length(radius)
            || !range(&angle_deg, -360.0, 360.0, 360.0)
            || !range(&height, -m, m, m)
        {
            return Err(Status::ERR_ARG);
        }
        self.run(
            || {
                format!(
                    "smlib_surface_cylinder(origin={}, xaxis={}, yaxis={}, radius={}, angle={}, height={})",
                    vec(&origin),
                    vec(&xaxis),
                    vec(&yaxis),
                    num(radius),
                    vec(&angle_deg),
                    vec(&height)
                )
            },
            &[],
            |api| {
                let (x, y) = (v3(&xaxis), v3(&yaxis));
                let z = unit_cross(x, y);
                let o: Point3 = v3(&origin) + z * height[0];
                let cylinder = Cone::new_cylinder_bounded(
                    o,
                    x,
                    y,
                    radius,
                    angle_deg[0],
                    angle_deg[1],
                    height[1] - height[0],
                    false,
                    false,
                    false,
                );
                finish(api.kernel(), KSurface::from(cylinder))
            },
        )
    }

    /// `smlib_surface_cone`: the cone with its apex at `apex` about the axis `xaxis` x
    /// `yaxis`, of half angle `half_angle_deg` (0.01 to 89), between `height[0]` and
    /// `height[1]` from the apex.
    ///
    /// # Errors
    ///
    /// `SMLIB_ERR_ARG` for a bad apex, axes, angle or heights (0 <= h0, h1 - h0 > 1e-9); a
    /// kernel status (`SmCone::CreateCanonical`) otherwise.
    pub fn surface_cone(
        &mut self,
        apex: [f64; 3],
        xaxis: [f64; 3],
        yaxis: [f64; 3],
        half_angle_deg: f64,
        height: [f64; 2],
    ) -> Result<Surface> {
        let m = SURFACE_MAX_COORD;
        let angle_ok = half_angle_deg.is_finite() && (0.01..=89.0).contains(&half_angle_deg);
        let height_ok = range(&height, 0.0, m, m) && height[1] - height[0] > 1e-9;
        if !bounded(&apex, m) || !axes(&xaxis, &yaxis) || !angle_ok || !height_ok {
            return Err(Status::ERR_ARG);
        }
        self.run(
            || {
                format!(
                    "smlib_surface_cone(apex={}, xaxis={}, yaxis={}, half_angle_deg={}, height={})",
                    vec(&apex),
                    vec(&xaxis),
                    vec(&yaxis),
                    num(half_angle_deg),
                    vec(&height)
                )
            },
            &[],
            |api| {
                let (x, y) = (v3(&xaxis), v3(&yaxis));
                let z = unit_cross(x, y);
                let [h0, h1] = height;
                let t = (half_angle_deg * core::f64::consts::PI / 180.0).tan();
                let o: Point3 = v3(&apex) + z * h0;
                let pl = Axis2Placement::new(o, x, y);
                let cone = Cone::create_canonical(
                    api.kernel(),
                    &pl,
                    h0 * t,
                    h1 * t,
                    h1 - h0,
                    false,
                    false,
                )?;
                finish(api.kernel(), KSurface::from(cone))
            },
        )
    }

    /// `smlib_surface_sphere`: the sphere of `radius` about `center`, its poles on the axis
    /// `xaxis` x `yaxis`, between longitudes `lon_deg` and latitudes `lat_deg`.
    ///
    /// # Errors
    ///
    /// `SMLIB_ERR_ARG` for a bad center, axes, radius or range (longitude within -360..360,
    /// latitude within -90..90); a kernel status otherwise.
    pub fn surface_sphere(
        &mut self,
        center: [f64; 3],
        xaxis: [f64; 3],
        yaxis: [f64; 3],
        radius: f64,
        lon_deg: [f64; 2],
        lat_deg: [f64; 2],
    ) -> Result<Surface> {
        if !bounded(&center, SURFACE_MAX_COORD)
            || !axes(&xaxis, &yaxis)
            || !length(radius)
            || !range(&lon_deg, -360.0, 360.0, 360.0)
            || !range(&lat_deg, -90.0, 90.0, 180.0)
        {
            return Err(Status::ERR_ARG);
        }
        self.run(
            || {
                format!(
                    "smlib_surface_sphere(center={}, xaxis={}, yaxis={}, radius={}, lon={}, lat={})",
                    vec(&center),
                    vec(&xaxis),
                    vec(&yaxis),
                    num(radius),
                    vec(&lon_deg),
                    vec(&lat_deg)
                )
            },
            &[],
            |api| {
                let anal = Extent2d::from_bounds(lon_deg[0], lat_deg[0], lon_deg[1], lat_deg[1]);
                let sphere = Sphere::new(
                    api.kernel(),
                    v3(&center),
                    v3(&xaxis),
                    v3(&yaxis),
                    anal,
                    radius,
                    SphereOptions::default(),
                )?;
                finish(api.kernel(), KSurface::from(sphere))
            },
        )
    }

    /// `smlib_surface_torus`: the torus about `center` and the axis `xaxis` x `yaxis`, of
    /// radii `major_radius` > `minor_radius`, over `major_deg` round the axis and
    /// `minor_deg` round the tube.
    ///
    /// # Errors
    ///
    /// `SMLIB_ERR_ARG` for a bad center, axes, radii or range; a kernel status otherwise.
    pub fn surface_torus(
        &mut self,
        center: [f64; 3],
        xaxis: [f64; 3],
        yaxis: [f64; 3],
        major_radius: f64,
        minor_radius: f64,
        major_deg: [f64; 2],
        minor_deg: [f64; 2],
    ) -> Result<Surface> {
        if !bounded(&center, SURFACE_MAX_COORD)
            || !axes(&xaxis, &yaxis)
            || !length(major_radius)
            || !length(minor_radius)
            || minor_radius >= major_radius
            || !range(&major_deg, -360.0, 360.0, 360.0)
            || !range(&minor_deg, -360.0, 360.0, 360.0)
        {
            return Err(Status::ERR_ARG);
        }
        self.run(
            || {
                format!(
                    "smlib_surface_torus(center={}, xaxis={}, yaxis={}, R={}, r={}, major={}, minor={})",
                    vec(&center),
                    vec(&xaxis),
                    vec(&yaxis),
                    num(major_radius),
                    num(minor_radius),
                    vec(&major_deg),
                    vec(&minor_deg)
                )
            },
            &[],
            |api| {
                let anal = Extent2d::from_bounds(major_deg[0], minor_deg[0], major_deg[1], minor_deg[1]);
                let torus = Torus::new(
                    api.kernel(),
                    v3(&center),
                    v3(&xaxis),
                    v3(&yaxis),
                    anal,
                    major_radius,
                    minor_radius,
                    TorusOptions::default(),
                )?;
                finish(api.kernel(), KSurface::from(torus))
            },
        )
    }

    /// `smlib_surface_nurbs`: the (rational) B-spline surface of degrees `degree_u`,
    /// `degree_v` (1 to [`SURFACE_MAX_DEGREE`]) with `nu` x `nv` control points
    /// (`points`, `[x, y, z]` per point, u index slowest), their `weights` (each 1e-6 to
    /// 1e6; `None` is all 1) and the knot vectors (`nu + degree_u + 1` and
    /// `nv + degree_v + 1` values). An unclamped knot vector gives the surface over its
    /// domain `[k[p], k[n]]`, clamped.
    ///
    /// # Errors
    ///
    /// `SMLIB_ERR_ARG` for a bad degree, count (more than [`SURFACE_MAX_POLES`]), knot
    /// vector, point or weight; a kernel status otherwise.
    #[allow(clippy::too_many_arguments)] // the C API's arguments
    pub fn surface_nurbs(
        &mut self,
        degree_u: i32,
        degree_v: i32,
        nu: i32,
        nv: i32,
        knots_u: &[f64],
        knots_v: &[f64],
        points: &[f64],
        weights: Option<&[f64]>,
    ) -> Result<Surface> {
        let deg_ok = |d: i32| (1..=SURFACE_MAX_DEGREE).contains(&d);
        let count_ok = |n: i32, d: i32| n > d && n <= SURFACE_MAX_POLES;
        if !deg_ok(degree_u)
            || !deg_ok(degree_v)
            || !count_ok(nu, degree_u)
            || !count_ok(nv, degree_v)
            || i64::from(nu) * i64::from(nv) > i64::from(SURFACE_MAX_POLES)
        {
            return Err(Status::ERR_ARG);
        }
        // (all positive, and the product is within SURFACE_MAX_POLES)
        let to = |x: i32| usize::try_from(x).map_err(|_| Status::ERR_ARG);
        let (pu, pv, cu, cv) = (to(degree_u)?, to(degree_v)?, to(nu)?, to(nv)?);
        if !knots(knots_u, cu, pu) || !knots(knots_v, cv, pv) {
            return Err(Status::ERR_ARG);
        }
        let n = cu * cv;
        if points.len() != 3 * n || !bounded(points, SURFACE_MAX_COORD) {
            return Err(Status::ERR_ARG);
        }
        if let Some(w) = weights {
            if w.len() != n
                || !w
                    .iter()
                    .all(|&x| x.is_finite() && (1e-6..=1e6).contains(&x))
            {
                return Err(Status::ERR_ARG);
            }
        }
        let weight = |i: usize| weights.and_then(|w| w.get(i).copied()).unwrap_or(1.0);
        self.run(
            || format!("smlib_surface_nurbs(degree=[{degree_u},{degree_v}], nu={nu}, nv={nv})"),
            &[],
            |api| {
                let surface = if !clamped(knots_u, pu) || !clamped(knots_v, pv) {
                    // NLib's surfaces are clamped (sm_DecomposeSrf, which every surface
                    // cache runs, reads past the control net otherwise). SMLib's
                    // component constructor clamps the surface to its domain [k[p], k[n]]
                    // (N_Iges128Srf), the same surface there.
                    let cpts: Vec<Point3> = points
                        .as_chunks::<3>()
                        .0
                        .iter()
                        .map(|&[x, y, z]| Point3::new(x, y, z))
                        .collect();
                    let wts: Vec<f64> = (0..n).map(weight).collect();
                    let (ku, mu) = unique_knots(knots_u);
                    let (kv, mv) = unique_knots(knots_v);
                    let bs = BSplineSurface::create_canonical(&CanonicalSurfaceData {
                        u_degree: pu,
                        v_degree: pv,
                        control_points: &cpts,
                        form: BSplineSurfaceForm::Unspecified,
                        u_multiplicities: &mu,
                        v_multiplicities: &mv,
                        u_knots: &ku,
                        v_knots: &kv,
                        knot_type: KnotType::Unspecified,
                        weights: Some(&wts),
                        uv_domain: None,
                    })?;
                    KSurface::from(bs)
                } else {
                    // What the runner does: NLib's surface filled directly (homogeneous
                    // control points), sidestepping SmBSplineSurface's component
                    // constructor.
                    let mut net: Vec<Vec<CPoint>> = Vec::new();
                    for (i, row) in points.chunks_exact(3 * cv).enumerate() {
                        let mut r = Vec::new();
                        for (j, &[x, y, z]) in row.as_chunks::<3>().0.iter().enumerate() {
                            let w = weight(i * cv + j);
                            r.push(CPoint::new(x * w, y * w, z * w, w));
                        }
                        net.push(r);
                    }
                    let nurb = NurbSurface::new(net, pu, pv, knots_u.to_vec(), knots_v.to_vec())?;
                    KSurface::from(BSplineSurface::from_nurb(nurb, false))
                };
                finish(api.kernel(), surface)
            },
        )
    }

    /// `smlib_surface_domain`: the surface's UV domain as `[u0, u1, v0, v1]`, what
    /// [`Smlib::ssi`] intersects within.
    #[must_use]
    pub fn surface_domain(s: &Surface) -> [f64; 4] {
        let d = s.domain;
        [d.min().x, d.max().x, d.min().y, d.max().y]
    }

    /// `smlib_surface_evaluate`: the point of `s` at `(u, v)`.
    ///
    /// # Errors
    ///
    /// `SMLIB_ERR_ARG` for a non-finite parameter; the kernel's status, or `SM_ERR` for a
    /// non-finite point.
    pub fn surface_evaluate(&mut self, s: &Surface, u: f64, v: f64) -> Result<[f64; 3]> {
        if !u.is_finite() || !v.is_finite() {
            return Err(Status::ERR_ARG);
        }
        self.run(
            || format!("smlib_surface_evaluate(u={}, v={})", num(u), num(v)),
            &[],
            |api| {
                let q = s.surface.evaluate_point(api.kernel(), Point2::new(u, v))?;
                let r = [q.x, q.y, q.z];
                if !finite(&r) {
                    return Err(Fail::Status(Status::SM_ERR));
                }
                Ok(r)
            },
        )
    }

    /// `smlib_ssi`: the intersection of `a` and `b` over their domains by `method`
    /// ([`SSI_GLOBAL`] or [`SSI_ADVANCED`]; [`SSI_BASIC`] is not in this port) within
    /// `approx_tol` (0 or less: SMLib's default), surface edges used on both.
    ///
    /// # Errors
    ///
    /// `SMLIB_ERR_ARG` for another method or a non-finite tolerance;
    /// `SMLIB_ERR_UNSUPPORTED` for [`SSI_BASIC`]; the kernel's status otherwise.
    pub fn ssi(
        &mut self,
        a: &Surface,
        b: &Surface,
        method: i32,
        approx_tol: f64,
    ) -> Result<SsiResult> {
        if !(SSI_GLOBAL..=SSI_BASIC).contains(&method) || !approx_tol.is_finite() {
            return Err(Status::ERR_ARG);
        }
        self.run(
            || format!("smlib_ssi(method={method}, approx_tol={})", num(approx_tol)),
            &[],
            |api| {
                let k = api.kernel();
                let tol = (approx_tol > 0.0).then(|| approx_tol3d(k, approx_tol));
                let edges = [true, true];
                let x: SurfaceIntersection = match method {
                    SSI_GLOBAL => a.surface.global_surface_intersect(
                        k, &a.domain, &b.surface, &b.domain, edges, tol, None,
                    )?,
                    SSI_ADVANCED => do_intersection(
                        k, &a.surface, &a.domain, &b.surface, &b.domain, edges, tol, None,
                    )?,
                    _ => return Err(Fail::Status(Status::ERR_UNSUPPORTED)),
                };
                Ok(result(k, &x))
            },
        )
    }
}

/// The C API's copy of an intersection: each curve's type, flags, deviation, interval
/// and B-splines, and a touching point per [`SSI_TOUCHING`] curve.
fn result(k: &Kernel, x: &SurfaceIntersection) -> SsiResult {
    let mut r = SsiResult::default();
    for (i, c3) in x.curves_3d.iter().enumerate() {
        let ua = x.surface1_uv_curves.get(i).and_then(Option::as_ref);
        let ub = x.surface2_uv_curves.get(i).and_then(Option::as_ref);
        let curve_type = x.curve_types.get(i).map_or(SSI_NONE, |t| *t as i32);
        let ivl = c3.natural_interval();
        let info = SsiCurveInfo {
            curve_type,
            flags: (if ua.is_some() { SSI_HAS_UV_A } else { 0 })
                | (if ub.is_some() { SSI_HAS_UV_B } else { 0 }),
            deviation: x.deviations.get(i).copied().unwrap_or(f64::NAN),
            interval: [ivl.min(), ivl.max()],
        };
        if curve_type == TsectCurveType::Touching as i32 {
            r.points.push(SsiPoint {
                xyz: start::<3>(k, Some(c3)),
                uv_a: start::<2>(k, ua),
                uv_b: start::<2>(k, ub),
            });
        }
        r.curves.push(SsiCurve {
            info,
            c3: nurbs(Some(c3), 3),
            uv_a: nurbs(ua, 2),
            uv_b: nurbs(ub, 2),
        });
    }
    r
}
