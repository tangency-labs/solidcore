//! Split Face and Split Body.
//!
//! [`Smlib::split_faces`] cuts chosen faces where tool surfaces meet them, the way SMLib's
//! own plane cutter (`SmPlaneCutter::MergeIntersection`) cuts a subset of faces: the
//! face's surface is intersected with the tool over the face's UV domain (the planar
//! section for a plane reaching over the face, `GlobalSurfaceIntersect` otherwise), and
//! the untrimmed curves go to `SmBrep::MergeCurvesOnSurface` with the tool as the
//! intersecting surface, which trims them to the face, cuts the edges they end on and
//! splits the face. Unlike the plane cutter, a failed merge fails the call.
//!
//! [`Smlib::split_body`] slices a solid with a one-face sheet on the tool surface
//! (`SM_BO_SLICE`, which keeps the solid on the side the sheet's normal points away from
//! and closes it with the sheet's piece inside it), once with the sheet and once with it
//! flipped.
//!
//! Source: tangency-labs/smlib `capi/src/smlib_split.cpp` (Apache-2.0, Copyright
//! SolidCraft contributors) at fork commit 3e03bd2 (replay/1.0.1).
//!
//! | Rust | C API (SMLib call) |
//! |------|-------|
//! | [`Smlib::split_faces`] | `smlib_split_faces` (`MergeCurvesOnSurface`) |
//! | [`Smlib::split_body`] | `smlib_split_body` (`SmApiBooleanWithOptions`, `SM_BO_SLICE`) |
//! | `plane_covers`, `cut_face` (private) | the same |
//!
//! Deviations: the faces are given by index into the body's face list and the tools as
//! [`Surface`] values; the counts of new faces per chosen face come back with the body
//! (the C API's `made`). The new objects are named by an arena-only notification store.

use std::collections::BTreeSet;

use solidcore::boolean::api_brep::boolean_with_options;
use solidcore::boolean::merge::MergeOptions;
use solidcore::geom::a_object::BooleanOperation;
use solidcore::geom::extent3d::Extent3d;
use solidcore::geom::object::ObjectRef;
use solidcore::geom::tol::{approx_tol3d, USE_DEFAULT};
use solidcore::geom::vector2d::Point2;
use solidcore::math::EFF_ZERO;
use solidcore::query::api_queries::brep_is_manifold_solid;
use solidcore::surface::plane::Plane;
use solidcore::surface::surface::p2::TsectCurveType;
use solidcore::topology::brep::p1 as bp1;
use solidcore::topology::brep::p10::merge_curves_on_surface_through_surface;
use solidcore::topology::brep::p2::ArenaNotifyState;
use solidcore::topology::brep::p3::new_brep;
use solidcore::topology::brep::p8::{create_face_from_surface, FaceBuildEnv};
use solidcore::topology::brep::planar_section::face_planar_section_curves;
use solidcore::topology::edge::NotifyEnv;
use solidcore::topology::face::p2 as face2;
use solidcore::topology::loopuse::notify_ref;
use solidcore::topology::topology_traverser::Marks;
use solidcore::topology::{Brep, EdgeId, FaceId};
use solidcore::{EntityRef, Kernel};

use crate::session::{take, Fail, Smlib};
use crate::ssi::Surface;
use crate::status::{Result, Status};

/// `SMLIB_SPLIT_MAX_FACES`.
pub const SPLIT_MAX_FACES: usize = 100_000;
/// `SMLIB_SPLIT_MAX_TOOLS`.
pub const SPLIT_MAX_TOOLS: usize = 4096;

/// The C API's `list` for log lines: at most 16 entries.
fn list(v: &[usize]) -> String {
    let parts: Vec<String> = v.iter().take(16).map(ToString::to_string).collect();
    let more = if v.len() > 16 { ",..." } else { "" };
    format!("[{}{more}]", parts.join(","))
}

/// Whether the plane `tool` reaches over all of `bx` (each corner's foot on the plane
/// lies in the tool's parameter box): then its section with a face in the box is the
/// infinite plane's.
fn plane_covers(k: &Kernel, tool: &Surface, bx: &Extent3d, tol: f64) -> bool {
    let [u0, u1, v0, v1] = Smlib::surface_domain(tool);
    let s = tool.surface();
    let (Ok(o), Ok(pu), Ok(pv)) = (
        s.evaluate_point(k, Point2::new(u0, v0)),
        s.evaluate_point(k, Point2::new(u1, v0)),
        s.evaluate_point(k, Point2::new(u0, v1)),
    ) else {
        return false;
    };
    // The plane is affine in (u, v): P = o + (u - u0) X + (v - v0) Y.
    let x = (pu - o) / (u1 - u0);
    let y = (pv - o) / (v1 - v0);
    let (a, b, c) = (x.dot(x), x.dot(y), y.dot(y));
    let det = a * c - b * b;
    let nonzero = det.abs() > 0.0;
    if !nonzero {
        return false;
    }
    let (m, big) = (bx.min(), bx.max());
    for corner in 0..8 {
        let q = solidcore::geom::vector3d::Point3::new(
            if corner & 1 != 0 { big.x } else { m.x },
            if corner & 2 != 0 { big.y } else { m.y },
            if corner & 4 != 0 { big.z } else { m.z },
        );
        let r = q - o;
        let (g0, g1) = (x.dot(r), y.dot(r));
        let u = u0 + (c * g0 - b * g1) / det;
        let v = v0 + (a * g1 - b * g0) / det;
        let (su, sv) = (tol / a.sqrt(), tol / c.sqrt());
        let inside = u >= u0 - su && u <= u1 + su && v >= v0 - sv && v <= v1 + sv;
        if !inside {
            return false;
        }
    }
    true
}

/// `cut_face`: cuts `face` where `tool` meets it; the faces it was split into besides
/// itself.
fn cut_face(
    k: &mut Kernel,
    b: &mut Brep,
    env: &mut NotifyEnv<'_>,
    curve_ref: &dyn Fn(EdgeId, i64) -> ObjectRef,
    face: FaceId,
    tool: &Surface,
    tol: f64,
) -> core::result::Result<Vec<FaceId>, Fail> {
    let surface = face2::surface(b, face)?
        .ok_or(Fail::Status(Status::SM_ERR))?
        .clone_for_evaluation()?;
    let domain = face2::uv_domain(b, face)?;
    let approx = approx_tol3d(k, approx_tol3d(k, tol) / 4.0);
    let bx = solidcore::query::face::p1::calculate_bounding_box(k, &mut *b, face, false)?;
    // A plane reaching over the face sections it exactly (SmPlaneCutter's way); a plane
    // strip, like any other tool, is intersected over its own domain.
    let planar = if tool.surface().is_kind_of(Plane::TYPE)? {
        let p = tool.surface().is_planar(k, EFF_ZERO)?;
        (p.is_planar && plane_covers(k, tool, &bx, tol)).then_some(p)
    } else {
        None
    };
    let (mut c3, uv1) = if let Some(p) = planar {
        // As SmPlaneCutter: a face far from the plane is not sectioned.
        if bx.distance_to_plane(p.point, p.normal).abs() > tol {
            return Ok(Vec::new());
        }
        let section = face_planar_section_curves(
            k,
            b,
            face,
            &surface,
            &domain,
            p.point,
            p.normal,
            Some(approx),
            None,
        )?;
        let uv = (section.uv_curves.len() == section.curves_3d.len()).then_some(section.uv_curves);
        (section.curves_3d, uv)
    } else {
        let x = surface.global_surface_intersect(
            k,
            &domain,
            tool.surface(),
            tool.domain(),
            [true, true],
            Some(approx),
            None,
        )?;
        // A single point of contact cuts nothing.
        let mut c3 = Vec::new();
        let mut uv1 = Vec::new();
        for (i, c) in x.curves_3d.into_iter().enumerate() {
            if x.curve_types.get(i) == Some(&TsectCurveType::Touching) {
                continue;
            }
            c3.push(c);
            uv1.push(x.surface1_uv_curves.get(i).cloned().flatten());
        }
        let all_uv = uv1.iter().all(Option::is_some);
        (c3, all_uv.then_some(uv1))
    };
    if c3.is_empty() {
        return Ok(Vec::new());
    }
    let (made, _on) = merge_curves_on_surface_through_surface(
        k,
        b,
        env,
        curve_ref,
        face,
        tol,
        &mut c3,
        uv1.as_deref(),
        Some(tool.surface().clone()),
    )?;
    Ok(made)
}

impl Smlib {
    /// `smlib_split_faces`: a new body, `body` with each face `faces[i]` (indices into its
    /// face list, distinct) cut where each of `tools` in turn meets it, the pieces kept as
    /// faces; and per chosen face the number of new faces it was split into. `body` is
    /// left unchanged.
    ///
    /// # Errors
    ///
    /// `SMLIB_ERR_ARG` for 0 or more than [`SPLIT_MAX_FACES`] faces, 0 or more than
    /// [`SPLIT_MAX_TOOLS`] tools, or a face index out of range or repeated;
    /// `SMLIB_ERR_NO_CUT` when no tool cuts any chosen face; the kernel's status otherwise
    /// (a failed merge fails the call).
    pub fn split_faces(
        &mut self,
        body: &Brep,
        faces: &[usize],
        tools: &[&Surface],
    ) -> Result<(Brep, Vec<i64>)> {
        if faces.is_empty()
            || faces.len() > SPLIT_MAX_FACES
            || tools.is_empty()
            || tools.len() > SPLIT_MAX_TOOLS
        {
            return Err(Status::ERR_ARG);
        }
        self.run(
            || {
                format!(
                    "smlib_split_faces(body, faces={}, n_tools={})",
                    list(faces),
                    tools.len()
                )
            },
            &[("body", body)],
            |api| {
                let count = bp1::num_faces(body);
                let mut seen = BTreeSet::new();
                for &f in faces {
                    if f >= count || !seen.insert(f) {
                        return Err(Fail::Status(Status::ERR_ARG));
                    }
                }
                let mut a = api.copy(body)?;
                let all = bp1::faces(&a, None)?;
                if all.len() != count {
                    return Err(Fail::Status(Status::SM_ERR));
                }
                // Each chosen face and, as tools cut it, its pieces.
                let mut pieces: Vec<Vec<FaceId>> = Vec::new();
                for &f in faces {
                    pieces.push(vec![*all.get(f).ok_or(Fail::Status(Status::ERR_ARG))?]);
                }
                let mut made = vec![0_i64; faces.len()];
                let editing = a.editing_enabled();
                a.set_editing_enabled(true);
                let k = api.kernel_mut();
                let doing_boolean = k.context.doing_boolean();
                k.context.set_doing_boolean(true);
                let tol = a.tolerance();
                let mut state = ArenaNotifyState::new(&a);
                let mut env = state.env(None);
                let curve_ref = |e: EdgeId, _type_code: i64| notify_ref(EntityRef::Edge(e));
                let mut cut = || -> core::result::Result<i64, Fail> {
                    let mut total = 0_i64;
                    for tool in tools {
                        for (piece, n) in pieces.iter_mut().zip(made.iter_mut()) {
                            let now = piece.clone();
                            for f in now {
                                let split =
                                    cut_face(k, &mut a, &mut env, &curve_ref, f, tool, tol)?;
                                let added = i64::try_from(split.len())
                                    .map_err(|_| Fail::Status(Status::SM_ERR))?;
                                piece.extend(split);
                                *n += added;
                                total += added;
                            }
                        }
                    }
                    Ok(total)
                };
                let total = cut();
                k.context.set_doing_boolean(doing_boolean);
                a.set_editing_enabled(editing);
                if total? == 0 {
                    return Err(Fail::Status(Status::ERR_NO_CUT));
                }
                Ok((take(a)?, made))
            },
        )
    }

    /// `smlib_split_body`: `body` (a manifold solid) cut by the surface `tool` into the part
    /// behind it and the part in front of it (the side its normal points to), each closed
    /// with the tool's piece inside the body. `body` is left unchanged.
    ///
    /// # Errors
    ///
    /// `SMLIB_ERR_NOT_SOLID` for a body that is not a manifold solid; `SMLIB_ERR_NO_CUT`
    /// when either side is empty; the kernel's status otherwise.
    pub fn split_body(&mut self, body: &Brep, tool: &Surface) -> Result<(Brep, Brep)> {
        self.run(
            || "smlib_split_body(body)".to_string(),
            &[("body", body)],
            |api| {
                if !brep_is_manifold_solid(body)? {
                    return Err(Fail::Status(Status::ERR_NOT_SOLID));
                }
                let back = api.split_side(body, tool, true)?;
                let front = api.split_side(body, tool, false)?;
                match (back, front) {
                    (Some(b), Some(f)) => Ok((b, f)),
                    _ => Err(Fail::Status(Status::ERR_NO_CUT)),
                }
            },
        )
    }

    /// One side: the body sliced by the tool's sheet, its normal turned away from it
    /// (`flip`: the surface turned over, its normal reversed; SMLib's slice reads the side
    /// from the surface, not the face's use). `None` when nothing is left.
    fn split_side(
        &mut self,
        body: &Brep,
        tool: &Surface,
        flip: bool,
    ) -> core::result::Result<Option<Brep>, Fail> {
        let a = self.copy(body)?;
        let k = self.kernel_mut();
        let mut surface = tool.surface().clone();
        if flip {
            surface.swap_uv(k)?;
        }
        let mut sheet = new_brep(k, USE_DEFAULT, USE_DEFAULT)?;
        sheet.set_tolerance(a.tolerance());
        let domain = surface.natural_uv_domain_with_nurb(k)?;
        {
            let mut state = ArenaNotifyState::new(&sheet);
            let mut env = state.env(None);
            let key_of: &dyn Fn(EntityRef) -> ObjectRef = env.key_of;
            let this_face = |f: FaceId| key_of(EntityRef::Face(f));
            let mut marks = Marks::default();
            let mut fenv = FaceBuildEnv {
                env: &mut env,
                marks: &mut marks,
                this_face: &this_face,
            };
            create_face_from_surface(k, &mut sheet, Some(surface), &domain, &mut fenv)?;
        }
        // The sheet is consumed by the slice.
        let options = MergeOptions::new(false, false, false, false, false);
        let r = boolean_with_options(k, a, sheet, BooleanOperation::Slice, options)?.result;
        if bp1::num_faces(&r) == 0 {
            return Ok(None);
        }
        Ok(Some(r))
    }
}
