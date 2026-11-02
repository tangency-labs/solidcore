//! Whether a NURBS surface lies on an analytic surface: the analytic surface's implicit
//! equation composed with the NURBS in homogeneous coordinates is a polynomial in Bernstein
//! form on each Bézier patch, which vanishes identically exactly when the patch lies on the
//! surface. Its coefficients are checked against their rounding (#2388: a support kept as
//! its analytic kind beside its NURBS form, the two refused when they disagree).
//!
//! Not upstream; the fork's `SmBallPatch::IsOnSphere` makes the sphere's check this way.
//!
//! Deviations: none (no upstream counterpart). No recursion.

use crate::error::{Result, SolidError};
use crate::geom::axis2_placement::Axis2Placement;
use crate::geom::extent2d::Extent2d;
use crate::kernel::Kernel;
use crate::math::{Vec2, Vec3, PI};
use crate::nurbs::b_spline_surface::BSplineSurface;
use crate::nurbs::math::NOW;
use crate::nurbs::srf_tool::p1::decompose_to_bez;
use crate::surface::cone::Cone;
use crate::surface::plane::Plane;
use crate::surface::sphere::Sphere;
use crate::surface::surface::{Surface, SurfaceKind};
use crate::surface::torus::{Torus, TorusOptions};

const OP: &str = "surface::implicit::lies_on";

/// An analytic surface by its implicit equation's data.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Implicit {
    /// The plane through `origin` with unit `normal`: `n . X = 0`.
    Plane {
        /// A point of it.
        origin: Vec3,
        /// Its unit normal.
        normal: Vec3,
    },
    /// The cone (a cylinder for `tan_angle = 0`) about the unit `axis` through `origin`,
    /// of radius `radius + tan_angle (X . a)`: `|X - (X . a) a|^2 = (r + t X . a)^2`.
    Cone {
        /// A point of the axis, where the radius is `radius`.
        origin: Vec3,
        /// The unit axis.
        axis: Vec3,
        /// The radius at `origin`.
        radius: f64,
        /// The tangent of the semi-angle.
        tan_angle: f64,
    },
    /// The sphere about `centre`: `|X|^2 = r^2`.
    Sphere {
        /// Its centre.
        centre: Vec3,
        /// Its radius.
        radius: f64,
    },
    /// The torus about the unit `axis` through `centre`:
    /// `(|X|^2 + R^2 - r^2)^2 = 4 R^2 (|X|^2 - (X . a)^2)`.
    Torus {
        /// Its centre.
        centre: Vec3,
        /// The unit axis.
        axis: Vec3,
        /// The major radius R.
        major: f64,
        /// The minor radius r.
        minor: f64,
    },
}

/// A polynomial in Bernstein form over a patch, degrees `(p, q)`, `c[i (q + 1) + j]`.
#[derive(Debug, Clone, PartialEq)]
struct BPoly {
    p: usize,
    q: usize,
    c: Vec<f64>,
}

/// `C(n, k)` as a double.
fn binom(n: usize, k: usize) -> f64 {
    let mut r = 1.0_f64;
    for i in 0..k {
        #[allow(clippy::cast_precision_loss)]
        let f = (n - i) as f64 / (i + 1) as f64;
        r *= f;
    }
    r
}

impl BPoly {
    fn at(&self, i: usize, j: usize) -> f64 {
        self.c.get(i * (self.q + 1) + j).copied().unwrap_or(0.0)
    }

    /// The product, of degrees `(p + p', q + q')`.
    fn mul(&self, o: &Self) -> Self {
        let (p, q) = (self.p + o.p, self.q + o.q);
        let mut c = vec![0.0; (p + 1) * (q + 1)];
        for i in 0..=self.p {
            for j in 0..=self.q {
                for k in 0..=o.p {
                    for l in 0..=o.q {
                        let f = binom(self.p, i) * binom(o.p, k) / binom(p, i + k)
                            * (binom(self.q, j) * binom(o.q, l) / binom(q, j + l));
                        if let Some(x) = c.get_mut((i + k) * (q + 1) + j + l) {
                            *x += f * self.at(i, j) * o.at(k, l);
                        }
                    }
                }
            }
        }
        Self { p, q, c }
    }

    /// `a self + b o`, of the same degrees.
    fn combine(&self, a: f64, o: &Self, b: f64) -> Self {
        Self {
            p: self.p,
            q: self.q,
            c: self
                .c
                .iter()
                .zip(&o.c)
                .map(|(x, y)| a * x + b * y)
                .collect(),
        }
    }

    fn max_abs(&self) -> f64 {
        self.c.iter().fold(0.0_f64, |m, x| m.max(x.abs()))
    }
}

/// `gamma_k = k u / (1 - k u)`.
fn gamma(k: usize) -> f64 {
    #[allow(clippy::cast_precision_loss)]
    let ku = k as f64 * f64::EPSILON / 2.0;
    ku / (1.0 - ku)
}

/// Whether `form` lies on `a`: on every Bézier patch, every Bernstein coefficient of the
/// implicit equation in the patch's homogeneous coordinates `(X W, W)` (`X` from the
/// surface's origin) is within the rounding of computing it, `gamma_K` times the size of
/// its terms, `K` four times the operations of the products of the equation's degree d
/// (`K = 4 d (p + 1) (q + 1) + 16`), the stored points a few ulps each.
///
/// # Errors
///
/// `SM_ERR` without NURB data; those of the Bézier decomposition.
pub fn lies_on(form: &BSplineSurface, a: &Implicit) -> Result<bool> {
    let nurb = form.gw_nurb().ok_or_else(|| SolidError::sm_err(OP))?;
    let origin = match *a {
        Implicit::Plane { origin, .. } | Implicit::Cone { origin, .. } => origin,
        Implicit::Sphere { centre, .. } | Implicit::Torus { centre, .. } => centre,
    };
    for row in decompose_to_bez(nurb)? {
        for patch in row {
            let (p, q) = patch.degrees();
            let mut x = [Vec::new(), Vec::new(), Vec::new()];
            let mut w = Vec::new();
            for prow in patch.rows() {
                for c in prow {
                    // (stored weighted when rational: H, homogeneous either way)
                    let wi = if c.w == NOW { 1.0 } else { c.w };
                    // X W = H - O W
                    let xw = Vec3::new(c.x, c.y, c.z) - origin * wi;
                    x[0].push(xw.x);
                    x[1].push(xw.y);
                    x[2].push(xw.z);
                    w.push(wi);
                }
            }
            let poly = |c: Vec<f64>| BPoly { p, q, c };
            let [x0, x1, x2] = x;
            let (xs, ws) = ([poly(x0), poly(x1), poly(x2)], poly(w));
            let size = xs.iter().map(BPoly::max_abs).fold(0.0_f64, f64::max);
            let wmax = ws.max_abs();
            let dot = |v: Vec3| xs[0].combine(v.x, &xs[1], v.y).combine(1.0, &xs[2], v.z);
            let norm2 = || {
                xs[0]
                    .mul(&xs[0])
                    .combine(1.0, &xs[1].mul(&xs[1]), 1.0)
                    .combine(1.0, &xs[2].mul(&xs[2]), 1.0)
            };
            let w2 = ws.mul(&ws);
            let (f, scale, degree) = match *a {
                Implicit::Plane { normal, .. } => (dot(normal), size, 1),
                Implicit::Cone {
                    axis,
                    radius,
                    tan_angle,
                    ..
                } => {
                    let along = dot(axis);
                    let r = ws.combine(radius, &along, tan_angle);
                    let f = norm2().combine(1.0, &along.mul(&along), -1.0).combine(
                        1.0,
                        &r.mul(&r),
                        -1.0,
                    );
                    let s = size + (radius + tan_angle.abs() * size) * wmax;
                    (f, s * s, 2)
                }
                Implicit::Sphere { radius, .. } => {
                    let f = norm2().combine(1.0, &w2, -radius * radius);
                    let s = size + radius * wmax;
                    (f, s * s, 2)
                }
                Implicit::Torus {
                    axis, major, minor, ..
                } => {
                    let n2 = norm2();
                    let along = dot(axis);
                    let g = n2.combine(1.0, &w2, major * major - minor * minor);
                    let radial = n2.combine(1.0, &along.mul(&along), -1.0);
                    let f = g
                        .mul(&g)
                        .combine(1.0, &radial.mul(&w2), -4.0 * major * major);
                    let s = size + (major + minor) * wmax;
                    (f, s * s * s * s, 4)
                }
            };
            let tol = gamma(4 * degree * (p + 1) * (q + 1) + 16) * scale;
            let within = |v: &f64| {
                matches!(
                    v.abs().partial_cmp(&tol),
                    Some(core::cmp::Ordering::Less | core::cmp::Ordering::Equal)
                )
            };
            if !f.c.iter().all(within) {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

const OP_BUILD: &str = "surface::implicit::AnalyticSupport::build";

/// A support's analytic kind and parameters, beside its NURBS form (#2388): what a
/// record keeps of a plane, cylinder, cone, sphere or torus support.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AnalyticSupport {
    /// A plane `O + u su X + v sv Y` over its domain.
    Plane {
        /// Its placement.
        placement: Axis2Placement,
        /// `(su, sv)`.
        uv_scale: Vec2,
        /// The domain.
        domain: Extent2d,
    },
    /// A cylinder about the placement's z axis.
    Cylinder {
        /// Its placement.
        placement: Axis2Placement,
        /// The radius.
        radius: f64,
    },
    /// A cone about the placement's z axis.
    Cone {
        /// Its placement.
        placement: Axis2Placement,
        /// The radius at the origin.
        base_radius: f64,
        /// The semi-angle, degrees.
        semi_angle_deg: f64,
    },
    /// A sphere about the placement's origin.
    Sphere {
        /// Its placement.
        placement: Axis2Placement,
        /// The radius.
        radius: f64,
    },
    /// A torus about the placement's z axis.
    Torus {
        /// Its placement.
        placement: Axis2Placement,
        /// The major radius.
        major_radius: f64,
        /// The minor radius.
        minor_radius: f64,
    },
}

impl AnalyticSupport {
    /// The analytic kind of a surface, if it is a plane, cylinder, cone, sphere or torus.
    #[must_use]
    pub fn of(s: &Surface) -> Option<Self> {
        Some(match s.kind() {
            SurfaceKind::Plane(p) => Self::Plane {
                placement: *p.position(),
                uv_scale: p.uv_scale(),
                domain: p.step_uv_domain(),
            },
            SurfaceKind::Cone(c) if c.cylinder => Self::Cylinder {
                placement: c.rev.position,
                radius: c.base_radius,
            },
            SurfaceKind::Cone(c) => Self::Cone {
                placement: c.rev.position,
                base_radius: c.base_radius,
                semi_angle_deg: c.semi_angle_deg,
            },
            SurfaceKind::Sphere(sp) => Self::Sphere {
                placement: sp.base.position,
                radius: sp.radius,
            },
            SurfaceKind::Torus(t) => Self::Torus {
                placement: t.base.position,
                major_radius: t.major_radius,
                minor_radius: t.minor_radius,
            },
            _ => return None,
        })
    }

    /// Its implicit equation's data.
    fn implicit(&self) -> Implicit {
        let axis = |p: &Axis2Placement| p.x_axis().cross(p.y_axis());
        match self {
            Self::Plane { placement, .. } => Implicit::Plane {
                origin: placement.origin(),
                normal: axis(placement),
            },
            Self::Cylinder { placement, radius } => Implicit::Cone {
                origin: placement.origin(),
                axis: axis(placement),
                radius: *radius,
                tan_angle: 0.0,
            },
            Self::Cone {
                placement,
                base_radius,
                semi_angle_deg,
            } => Implicit::Cone {
                origin: placement.origin(),
                axis: axis(placement),
                radius: *base_radius,
                tan_angle: (semi_angle_deg * PI / 180.0).tan(),
            },
            Self::Sphere { placement, radius } => Implicit::Sphere {
                centre: placement.origin(),
                radius: *radius,
            },
            Self::Torus {
                placement,
                major_radius,
                minor_radius,
            } => Implicit::Torus {
                centre: placement.origin(),
                axis: axis(placement),
                major: *major_radius,
                minor: *minor_radius,
            },
        }
    }

    /// The analytic surface on the NURB `form`: refused unless `form` lies on it
    /// ([`lies_on`]) or its placement's axes are not orthonormal. The surface keeps its placement as given and takes `form` as its NURB,
    /// so it evaluates as `form` wherever its NURB is read.
    ///
    /// # Errors
    ///
    /// `InvalidInput` where `form` is off the analytic surface; those of the canonical
    /// constructors and of [`lies_on`].
    pub fn build(&self, k: &Kernel, form: BSplineSurface) -> Result<Surface> {
        // the placement's axes as stored: unit and square to rounding (8 ulps)
        let p = match self {
            Self::Plane { placement, .. }
            | Self::Cylinder { placement, .. }
            | Self::Cone { placement, .. }
            | Self::Sphere { placement, .. }
            | Self::Torus { placement, .. } => placement,
        };
        let ulps = 8.0 * f64::EPSILON;
        let unit = |v: Vec3| (v.length() - 1.0).abs() <= ulps;
        let square = p.x_axis().dot(p.y_axis()).abs() <= ulps;
        if !(unit(p.x_axis()) && unit(p.y_axis()) && square) {
            return Err(SolidError::invalid_input(
                OP_BUILD,
                "a support's placement is not orthonormal",
            ));
        }
        if !lies_on(&form, &self.implicit())? {
            return Err(SolidError::invalid_input(
                OP_BUILD,
                "a support's analytic kind disagrees with its NURBS form",
            ));
        }
        Ok(match *self {
            Self::Plane {
                placement,
                uv_scale,
                domain,
            } => Surface::from(Plane::from_stored(placement, uv_scale, domain, form)),
            Self::Cylinder { placement, radius } => {
                let mut c = Cone::create_canonical_cylinder(k, &placement, radius)?;
                c.rev.position = placement;
                c.rev.base = form;
                Surface::from(c)
            }
            Self::Cone {
                placement,
                base_radius,
                semi_angle_deg,
            } => {
                let mut c = Cone::create_canonical_with_semi_angle(
                    k,
                    &placement,
                    base_radius,
                    semi_angle_deg,
                    None,
                    false,
                    false,
                )?;
                c.rev.position = placement;
                c.rev.base = form;
                Surface::from(c)
            }
            Self::Sphere { placement, radius } => {
                let mut sp = Sphere::create_canonical(k, &placement, radius)?;
                sp.base.position = placement;
                sp.base.base = form;
                Surface::from(sp)
            }
            Self::Torus {
                placement,
                major_radius,
                minor_radius,
            } => {
                let mut t = Torus::create_canonical(
                    k,
                    &placement,
                    major_radius,
                    minor_radius,
                    TorusOptions {
                        swap_uv: false,
                        inside_out: false,
                        circle_nurb: None,
                    },
                    None,
                )?;
                t.base.position = placement;
                t.base.base = form;
                Surface::from(t)
            }
        })
    }
}

const OP_DERIVE: &str = "surface::implicit::AnalyticSupport::derive";

/// The analytic kind of a support, without its parameters (#2388: a record's label beside
/// the support's NURBS form, which is the support's one geometric truth).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SupportKind {
    /// A plane.
    Plane,
    /// A cylinder.
    Cylinder,
    /// A cone.
    Cone,
    /// A sphere.
    Sphere,
    /// A torus.
    Torus,
}

/// A NURBS surface's control points, Euclidean, `[i along u][j along v]`, with its
/// degrees.
fn euclidean_net(form: &BSplineSurface) -> Result<(Vec<Vec<Vec3>>, usize, usize)> {
    use crate::nurbs::b_spline_surface::SurfParamType;
    let (nu, nv, cps) = form.control_points()?;
    let mut net = vec![Vec::new(); nu];
    for (i, row) in net.iter_mut().enumerate() {
        for j in 0..nv {
            let c = cps
                .get(i * nv + j)
                .ok_or_else(|| SolidError::sm_err(OP_DERIVE))?;
            let w = if c.w == NOW { 1.0 } else { c.w };
            row.push(Vec3::new(c.x, c.y, c.z) * (1.0 / w));
        }
    }
    Ok((
        net,
        form.degree(SurfParamType::U),
        form.degree(SurfParamType::V),
    ))
}

/// The circle of a rational quadratic arc `p0, p1, p2` (tangent to `p0 p1` at `p0` and
/// to `p1 p2` at `p2`): its centre, radius and unit normal; the centre solves the two
/// tangency conditions in the arc's plane. `None` for a degenerate arc.
fn arc_circle(p0: Vec3, p1: Vec3, p2: Vec3) -> Option<(Vec3, f64, Vec3)> {
    let (a, b) = (p1 - p0, p1 - p2);
    let n = a.cross(p2 - p1);
    // rows a, b, n; right-hand side a.p0, b.p2, n.p0 (Cramer)
    let det = a.dot(b.cross(n));
    let scale = a.length() * b.length() * n.length();
    // (a pick of the arcs that span a plane, not a tolerance on the result)
    if det.abs().partial_cmp(&(1e-12 * scale)) != Some(core::cmp::Ordering::Greater) {
        return None;
    }
    let (ra, rb, rn) = (a.dot(p0), b.dot(p2), n.dot(p0));
    let centre = (b.cross(n) * ra + n.cross(a) * rb + a.cross(b) * rn) * (1.0 / det);
    Some((centre, (p0 - centre).length(), n * (1.0 / n.length())))
}

/// The first non-degenerate arc's circle of the row `k` of the net along `u` (`along_u`)
/// or `v`.
fn row_circle(net: &[Vec<Vec3>], along_u: bool, k: usize) -> Option<(Vec3, f64, Vec3)> {
    let at = |i: usize| -> Option<Vec3> {
        if along_u {
            net.get(i)?.get(k).copied()
        } else {
            net.get(k)?.get(i).copied()
        }
    };
    arc_circle(at(0)?, at(1)?, at(2)?)
}

/// The Greville abscissae of a NURBS surface along u (`along_u`) or v: each control
/// point's average of its `degree` inner knots.
fn greville(form: &BSplineSurface, along_u: bool) -> Result<Vec<f64>> {
    use crate::nurbs::b_spline_surface::SurfParamType;
    let dir = if along_u {
        SurfParamType::U
    } else {
        SurfParamType::V
    };
    let p = form.degree(dir);
    let (knots, mults) = form.knots(dir, None)?;
    let mut full = Vec::new();
    for (&k, &m) in knots.iter().zip(&mults) {
        full.extend(std::iter::repeat_n(k, m));
    }
    let n = full.len().saturating_sub(p + 1);
    #[allow(clippy::cast_precision_loss)]
    let pf = p as f64;
    Ok((0..n)
        .map(|i| full.iter().skip(i + 1).take(p).sum::<f64>() / pf)
        .collect())
}

/// A unit vector square to `a`.
fn square_to(a: Vec3) -> Vec3 {
    let t = if a.x.abs() < 0.9 { Vec3::X } else { Vec3::Y };
    let x = t - a * t.dot(a);
    x * (1.0 / x.length())
}

impl AnalyticSupport {
    /// The analytic surface of kind `kind` a NURBS form is, derived from the form: a
    /// plane from its corners (its scales and domain from the knots); a surface of
    /// revolution from its circles, the rational quadratic arcs of its rows. Two rows of
    /// the revolving direction are circles about one axis (parallel normals, centres on
    /// the axis); a cylinder's radius is theirs, a cone's grows between them, a sphere's
    /// and a torus' centre and radii come from a generating arc. Exact to rounding: a
    /// rational NURBS represents these surfaces exactly. [`AnalyticSupport::build`] then
    /// checks the form against the derived surface.
    ///
    /// # Errors
    ///
    /// `InvalidInput` where the form is not of that kind (no circles about one axis, a
    /// degenerate corner); those of the control net.
    pub fn derive(kind: SupportKind, form: &BSplineSurface) -> Result<Self> {
        let not_kind =
            || SolidError::invalid_input(OP_DERIVE, "a support's NURBS form is not of its kind");
        let (net, p, q) = euclidean_net(form)?;
        let (nu, nv) = (net.len(), net.first().map_or(0, Vec::len));
        let corner = |i: usize, j: usize| net.get(i).and_then(|r| r.get(j)).copied();
        if kind == SupportKind::Plane {
            let (Some(p00), Some(pu), Some(pv)) = (
                corner(0, 0),
                corner(nu.saturating_sub(1), 0),
                corner(0, nv.saturating_sub(1)),
            ) else {
                return Err(not_kind());
            };
            let d = form.natural_uv_domain();
            let (du, dv) = (d.max().x - d.min().x, d.max().y - d.min().y);
            let (eu, ev) = (pu - p00, pv - p00);
            let x = eu * (1.0 / eu.length());
            let ey = ev - x * ev.dot(x);
            let y = ey * (1.0 / ey.length());
            let (su, sv) = (eu.length() / du, ev.dot(y) / dv);
            let origin = p00 - x * (d.min().x * su) - y * (d.min().y * sv);
            let ok = [x, y, origin]
                .iter()
                .all(|v| v.x.is_finite() && v.y.is_finite() && v.z.is_finite())
                && su.is_finite()
                && sv.is_finite();
            if !ok {
                return Err(not_kind());
            }
            // the form is the plane's parametrisation: of equal weights (or none), each
            // control point the plane's affine point at its Greville parameters, as a
            // polynomial form is exactly where it reproduces the affine map; otherwise
            // its (u, v) is not the plane's, and it is refused
            let (_, _, cps) = form.control_points()?;
            let w0 = cps.first().map_or(NOW, |c| c.w);
            if cps.iter().any(|c| c.w != w0) {
                return Err(not_kind());
            }
            let (gu, gv) = (greville(form, true)?, greville(form, false)?);
            let size = net
                .iter()
                .flatten()
                .fold(0.0_f64, |a, p| a.max(p.length()))
                .max(origin.length());
            let tol = gamma(8 + p + q) * size;
            for (i, row) in net.iter().enumerate() {
                for (j, pt) in row.iter().enumerate() {
                    let (Some(&u), Some(&v)) = (gu.get(i), gv.get(j)) else {
                        return Err(not_kind());
                    };
                    let affine = origin + x * (u * su) + y * (v * sv);
                    let within = matches!(
                        (*pt - affine).length().partial_cmp(&tol),
                        Some(core::cmp::Ordering::Less | core::cmp::Ordering::Equal)
                    );
                    if !within {
                        return Err(not_kind());
                    }
                }
            }
            return Ok(Self::Plane {
                placement: Axis2Placement::from_stored(origin, x, y),
                uv_scale: Vec2::new(su, sv),
                domain: d,
            });
        }
        // the revolving direction: rows of quadratic arcs about one axis
        let mut found = None;
        for along_u in [true, false] {
            let (deg, rows) = if along_u { (p, nv) } else { (q, nu) };
            if deg != 2 {
                continue;
            }
            let circles: Vec<_> = (0..rows)
                .filter_map(|k| row_circle(&net, along_u, k))
                .collect();
            let (Some(&(ca, ra, na)), Some(&(cb, rb, nb))) = (circles.first(), circles.last())
            else {
                continue;
            };
            let size = ra.max(rb).max((cb - ca).length());
            // (these pick the revolving direction; the form is checked against the
            // derived surface by `build`)
            let parallel = na.cross(nb).length() <= 1e-9;
            let on_axis = (cb - ca).cross(na).length() <= 1e-9 * size.max(1.0);
            if parallel && on_axis {
                found = Some((along_u, ca, ra, cb, rb, na));
                break;
            }
        }
        let Some((along_u, ca, ra, cb, rb, axis)) = found else {
            return Err(not_kind());
        };
        let placement =
            |o: Vec3| Axis2Placement::from_stored(o, square_to(axis), axis.cross(square_to(axis)));
        // a generating arc: the first column across the revolving direction
        let generator = if along_u {
            row_circle(&net, false, 0)
        } else {
            row_circle(&net, true, 0)
        };
        Ok(match kind {
            SupportKind::Plane => return Err(not_kind()),
            SupportKind::Cylinder => Self::Cylinder {
                placement: placement(ca),
                radius: ra,
            },
            SupportKind::Cone => {
                let h = (cb - ca).dot(axis);
                if h.abs().partial_cmp(&0.0) != Some(core::cmp::Ordering::Greater) {
                    return Err(not_kind());
                }
                Self::Cone {
                    placement: placement(ca),
                    base_radius: ra,
                    semi_angle_deg: ((rb - ra) / h).atan() * 180.0 / PI,
                }
            }
            SupportKind::Sphere => {
                let Some((c, r, _)) = generator else {
                    return Err(not_kind());
                };
                Self::Sphere {
                    placement: placement(c),
                    radius: r,
                }
            }
            SupportKind::Torus => {
                let Some((g, minor, _)) = generator else {
                    return Err(not_kind());
                };
                let t = (g - ca).dot(axis);
                let centre = ca + axis * t;
                Self::Torus {
                    placement: placement(centre),
                    major_radius: (g - centre).length(),
                    minor_radius: minor,
                }
            }
        })
    }

    /// Its kind, without its parameters.
    #[must_use]
    pub const fn kind(&self) -> SupportKind {
        match self {
            Self::Plane { .. } => SupportKind::Plane,
            Self::Cylinder { .. } => SupportKind::Cylinder,
            Self::Cone { .. } => SupportKind::Cone,
            Self::Sphere { .. } => SupportKind::Sphere,
            Self::Torus { .. } => SupportKind::Torus,
        }
    }
}
