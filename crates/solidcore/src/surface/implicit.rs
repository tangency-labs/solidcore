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
