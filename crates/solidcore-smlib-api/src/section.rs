//! Curves on a face: a face's section by a plane, handed back as a curve list (NURBS
//! only) with each curve's interval.
//!
//! Source: tangency-labs/smlib `capi/src/smlib_sweeps.cpp` (`smlib_face_section`,
//! `smlib_curves`, `put_curve`, `smlib_curves_size`, `smlib_curves_read`; Apache-2.0,
//! Copyright SolidCraft contributors) at fork commit a3485b9 (replay/1.0.1).
//!
//! | Rust | C API (SMLib call) |
//! |------|-------|
//! | [`Curves`] | `smlib_curves` |
//! | [`Smlib::face_section`] | `smlib_face_section` (`SmBrep::CreatePlanarSectionCurves`) |
//! | [`Curves::size`], [`Curves::kinds`], [`Curves::reals`], [`Curves::intervals`] | `smlib_curves_size`, `smlib_curves_read` |
//! | `put_curve` (private) | `put_curve` |
//!
//! Deviations: no `smlib_curves_free` (the value is dropped); the arrays are read as
//! slices instead of copied into caller buffers of the exact sizes. The section runs on a
//! deep copy of the body (the C API's runs on the body itself, whose lazy UV trims and
//! caches it may add). `smlib_project_curves` follows with the curve-list parser.

use solidcore::curve::curve::{Curve as KCurve, CurveKind};
use solidcore::geom::vector3d::{Point3, Vec3};
use solidcore::nurbs::b_spline_curve::BSplineCurve;
use solidcore::nurbs::math::NOW;
use solidcore::topology::brep::planar_section::{
    create_planar_section_curves, PlanarSectionOptions, SectionOutputs,
};
use solidcore::topology::Brep;

use crate::pick::pick_faces;
use crate::session::{finite, num, vec, Fail, Smlib};
use crate::status::{Result, Status};

/// `SMLIB_CURVE_NURBS`: a curve-list entry's kind.
const CURVE_NURBS: i64 = 3;
/// `SMLIB_CURVE_MAX_DEGREE`.
const CURVE_MAX_DEGREE: usize = 25;
/// `SMLIB_CURVE_MAX_POINTS`.
const CURVE_MAX_POINTS: usize = 100_000;

/// `smlib_curves`: curves as a curve list (each `[SMLIB_CURVE_NURBS, degree, count]` in
/// the kinds, its control points `x, y, z, weight` then all its knots in the reals) and
/// each curve's interval.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Curves {
    kinds: Vec<i64>,
    reals: Vec<f64>,
    intervals: Vec<f64>,
}

impl Curves {
    /// `smlib_curves_size`: the lengths of the kinds and the reals, and the number of
    /// curves.
    #[must_use]
    pub fn size(&self) -> (usize, usize, usize) {
        (self.kinds.len(), self.reals.len(), self.intervals.len() / 2)
    }

    /// The kinds (`smlib_curves_read`).
    #[must_use]
    pub fn kinds(&self) -> &[i64] {
        &self.kinds
    }

    /// The reals (`smlib_curves_read`).
    #[must_use]
    pub fn reals(&self) -> &[f64] {
        &self.reals
    }

    /// Each curve's interval, `[min, max]` per curve (`smlib_curves_read`).
    #[must_use]
    pub fn intervals(&self) -> &[f64] {
        &self.intervals
    }
}

/// The `SmBSplineCurve` (`GetRootCurve`) of a curve kind derived from it.
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

/// `put_curve`: `c` appended to `out` as the NURBS it is: control points `(x, y, z,
/// weight)`, all knots, and its interval. `SMLIB_ERR_UNSUPPORTED` for a curve that is no
/// B-spline or whose knots are not clamped (the curve list takes clamped ones).
fn put_curve(c: &KCurve, out: &mut Curves) -> core::result::Result<(), Fail> {
    let unsupported = || Fail::Status(Status::ERR_UNSUPPORTED);
    let mut b = b_spline_base(c).ok_or_else(unsupported)?.clone_exact();
    let n = b.get_or_create_gw_nurb().map_err(|_| unsupported())?;
    let degree = n.degree();
    let cpts = n.cpts();
    let u = n.knots();
    let (count, knots) = (cpts.len(), u.len());
    if !(1..=CURVE_MAX_DEGREE).contains(&degree)
        || count < degree + 1
        || count > CURVE_MAX_POINTS
        || knots != count + degree + 1
    {
        return Err(unsupported());
    }
    let (Some(&first), Some(&last)) = (u.first(), u.last()) else {
        return Err(unsupported());
    };
    for i in 1..=degree {
        if u.get(i) != Some(&first) || u.get(knots - 1 - i) != Some(&last) {
            return Err(unsupported());
        }
    }
    let to = |x: usize| i64::try_from(x).map_err(|_| unsupported());
    out.kinds.extend([CURVE_NURBS, to(degree)?, to(count)?]);
    for q in cpts {
        let w = if q.w == NOW { 1.0 } else { q.w };
        out.reals.extend([q.x / w, q.y / w, q.z / w, w]);
    }
    out.reals.extend_from_slice(u);
    let ivl = c.natural_interval();
    out.intervals.extend([ivl.min(), ivl.max()]);
    Ok(())
}

impl Smlib {
    /// `smlib_face_section`: the curves where the plane through `plane_point` with normal
    /// `plane_normal` cuts the face of `brep` picked by `face_point` (within
    /// `pick_tolerance`), within `tolerance`, as NURBS.
    ///
    /// # Errors
    ///
    /// `SMLIB_ERR_ARG` for a non-finite point or normal, a zero normal, a pick tolerance
    /// not in (0, 1e6) or a tolerance not in [1e-9, 1e3); `SMLIB_ERR_PICK_NONE` /
    /// `SMLIB_ERR_PICK_AMBIGUOUS` for the pick; `SMLIB_ERR_EMPTY` when the plane misses the
    /// face; `SMLIB_ERR_UNSUPPORTED` for a section curve that is no clamped B-spline; the
    /// kernel's status otherwise.
    pub fn face_section(
        &mut self,
        brep: &Brep,
        face_point: [f64; 3],
        pick_tolerance: f64,
        plane_point: [f64; 3],
        plane_normal: [f64; 3],
        tolerance: f64,
    ) -> Result<Curves> {
        let ok = finite(&face_point)
            && finite(&plane_point)
            && finite(&plane_normal)
            && pick_tolerance.is_finite()
            && pick_tolerance > 0.0
            && pick_tolerance < 1e6
            && tolerance.is_finite()
            && (1e-9..1e3).contains(&tolerance);
        if !ok {
            return Err(Status::ERR_ARG);
        }
        let n = Vec3::new(plane_normal[0], plane_normal[1], plane_normal[2]);
        let long = n.length() > 1e-12;
        if !long {
            return Err(Status::ERR_ARG);
        }
        let n = n.unitize().map_err(|_| Status::ERR_ARG)?;
        self.run(
            || {
                format!(
                    "smlib_face_section(brep, face_point={}, pick_tolerance={}, plane_point={}, plane_normal={}, tolerance={})",
                    vec(&face_point),
                    num(pick_tolerance),
                    vec(&plane_point),
                    vec(&plane_normal),
                    num(tolerance)
                )
            },
            &[("brep", brep)],
            |api| {
                let mut b = api.copy(brep)?;
                let k = api.kernel_mut();
                let picked = pick_faces(k, &mut b, &[face_point], pick_tolerance)?;
                let face = *picked.items.first().ok_or(Fail::Status(Status::ERR_PICK_NONE))?;
                let options = PlanarSectionOptions {
                    approx_tol_3d: Some(tolerance),
                    ang_tol_rad: None,
                    outputs: SectionOutputs {
                        curves_3d: true,
                        uv_curves: true,
                        faces: true,
                    },
                    force_uv_classify: false,
                };
                let at = Point3::new(plane_point[0], plane_point[1], plane_point[2]);
                let section = create_planar_section_curves(k, &mut b, at, n, &options)?;
                let made = section.curves_3d.unwrap_or_default();
                let faces = section.faces.unwrap_or_default();
                if faces.len() != made.len() {
                    return Err(Fail::Status(Status::ERR_UNSUPPORTED));
                }
                let mut result = Curves::default();
                for (c, f) in made.iter().zip(&faces) {
                    if *f == face {
                        put_curve(c, &mut result)?;
                    }
                }
                if result.intervals.is_empty() {
                    return Err(Fail::Status(Status::ERR_EMPTY));
                }
                Ok(result)
            },
        )
    }
}
