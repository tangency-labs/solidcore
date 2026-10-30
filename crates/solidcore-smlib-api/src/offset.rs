//! The offsets: shells, thickened sheets, drafted faces, offset sheets and merged faces.
//!
//! Source: tangency-labs/smlib `capi/src/smlib_ops.cpp` (`smlib_shell`, `smlib_draft`,
//! `smlib_thicken`, `smlib_merge_faces`, `points_ok`, `tolerance_ok`, `pts`, `FaceSet`,
//! `leaves_cavity`, `outward_normal`, `drafted_cone`, `draft_faces`, `moved_faces`), branch
//! replay/1.0.1 at `d51a723`, and `capi/src/smlib_sheets.cpp` (`smlib_offset_sheet`) at
//! `b018013` (Apache-2.0, Copyright SolidCraft contributors).
//!
//! | Rust | C API (`SM_API` / SMLib call) |
//! |------|-------|
//! | [`Smlib::shell`] | `smlib_shell` (`SmApiShellBrepFull`) |
//! | [`Smlib::thicken`] | `smlib_thicken` (`SmApiOffsetBrepFull`, `SmApiShellBrepFull`) |
//! | [`Smlib::draft`] | `smlib_draft` (`SmBrep::LocalOperation`) |
//! | [`Smlib::offset_sheet`] | `smlib_offset_sheet` (`SmApiOffsetBrepFull`) |
//! | [`Smlib::merge_faces`] | `smlib_merge_faces` (`SmBrep::RemoveTopologicalEdgesAndVertices`) |
//! | `MAX_POINTS` (private) | `MAX_POINTS` |
//! | `points_ok`, `tolerance_ok`, `pts` (private) | `points_ok`, `tolerance_ok`, `pts` |
//! | `FaceSet`, `leaves_cavity` (private) | `FaceSet`, `leaves_cavity` |
//! | `outward_normal`, `drafted_cone` (private) | `outward_normal`, `drafted_cone` |
//! | `draft_faces`, `moved_faces` (private) | `draft_faces`, `moved_faces` |
//!
//! Deviations:
//! - **Selections** are `[x, y, z]` points; a point that picks no face (edge) or cannot
//!   decide is the call's status, its index [`Smlib::last_failed_point`] (the C API's
//!   `out_failed_point`).
//! - **Draft's no-effect test.** The C API tags the tilted surfaces it hands
//!   `LocalOperation` with an attribute and reports `SMLIB_ERR_NO_EFFECT` when no face of
//!   the result carries one. Surfaces carry no attributes here: a face counts as moved
//!   when the middle of its surface's domain lies on one of the surfaces handed over
//!   (within the body's tolerance): on a tilted plane, or on a drafted cone (its radius
//!   at that height along the axis).
//! - **Inversion.** `outward_normal`'s `STEPInversion` counts as answered when it
//!   succeeds.
//! - **Sheet side walls.** The fork's shell joins consecutive side walls of a thickened
//!   sheet that lie on one plane; here that is the opt-in
//!   `ShellOptions::merge_sheet_side_walls`, which the C API's calls set.
//! - **Merge status.** A failed `RemoveTopologicalEdgesAndVertices` is its SMLib status
//!   unchanged, as the C API returns it (not through `smlib_shim::status`).

use solidcore::geom::axis2_placement::Axis2Placement;
use solidcore::intersect::point_class::PointClassType;
use solidcore::math::{Point3, Vec3};
use solidcore::nurbs::b_spline_curve::ControlPointForm;
use solidcore::nurbs::b_spline_surface::SurfParamType;
use solidcore::nurbs::b_spline_surface_n_lib::p2::OrientType;
use solidcore::nurbs::knots::Side;
use solidcore::offset::offset_executive::p2::{
    offset_brep_simple, shell_brep_simple, OffsetOptions, ShellOptions,
};
use solidcore::query::api_queries::{
    brep_is_manifold_solid, face_bounding_box, face_classify_uv, face_internal_point,
};
use solidcore::surface::cone::Cone;
use solidcore::surface::surface::{Surface, SurfaceKind};
use solidcore::topology::brep::local_operation::local_operation_simple;
use solidcore::topology::brep::p5::remove_topological_edges_and_vertices;
use solidcore::topology::EdgeId;
use solidcore::topology::{Brep, FaceId};
use solidcore::Kernel;

use crate::pick::{face_distance, outer_faceuse, pick_edges, pick_faces};
use crate::session::{finish, finite, num, take, vec, Fail, Smlib};
use crate::status::{Result, Status};

/// The most selection points a call takes (the C API's `MAX_POINTS`).
const MAX_POINTS: usize = 4096;

/// `points_ok`: at most `MAX_POINTS` finite points, at least one unless `allow_none`.
fn points_ok(points: &[[f64; 3]], allow_none: bool) -> bool {
    if points.len() > MAX_POINTS || (points.is_empty() && !allow_none) {
        return false;
    }
    points.iter().all(|p| finite(p))
}

/// `tolerance_ok`.
fn tolerance_ok(t: f64) -> bool {
    t.is_finite() && t > 0.0 && t < 1e6
}

/// `pts`: the points for the log line.
fn pts(points: &[[f64; 3]]) -> String {
    if points.is_empty() {
        return "[]".to_string();
    }
    let flat: Vec<f64> = points.iter().flatten().copied().collect();
    vec(&flat)
}

/// `std::max`: `b` when `a` is less, else `a`.
fn max2(a: f64, b: f64) -> f64 {
    if a < b {
        b
    } else {
        a
    }
}

/// The shell options of the C API's `SmApiShellBrepFull(.., TRUE, TRUE, TRUE, ..)`.
fn shell_options(thickness: f64) -> ShellOptions {
    ShellOptions {
        offset: OffsetOptions {
            distance: thickness,
            extend_convex_edges: true,
            merge_results: true,
        },
        create_offset_solid: true,
        merge_sheet_side_walls: true,
    }
}

/// The offset options of the C API's `SmApiOffsetBrepFull(.., TRUE, TRUE, ..)`.
fn offset_options(distance: f64) -> OffsetOptions {
    OffsetOptions {
        distance,
        extend_convex_edges: true,
        merge_results: true,
    }
}

/// Faces with their loose bounding boxes, for nearest-face queries.
struct FaceSet {
    faces: Vec<FaceId>,
    boxes: Vec<(Point3, Point3)>,
}

impl FaceSet {
    const fn new() -> Self {
        Self {
            faces: Vec::new(),
            boxes: Vec::new(),
        }
    }

    fn add(&mut self, k: &Kernel, b: &Brep, f: FaceId) -> bool {
        let Ok(e) = face_bounding_box(k, b, f, false) else {
            return false;
        };
        self.faces.push(f);
        self.boxes.push((e.min(), e.max()));
        true
    }

    /// The distance from `p` to the nearest face, with the nearest point; infinite when
    /// none answers. A face whose box is farther than the nearest face so far is skipped.
    fn nearest(&self, k: &mut Kernel, b: &mut Brep, p: Point3) -> (f64, Option<Point3>) {
        let mut best = f64::INFINITY;
        let mut at = None;
        for (&f, &(lo, hi)) in self.faces.iter().zip(&self.boxes) {
            let mut bx = 0.0;
            for (l, h, x) in [(lo.x, hi.x, p.x), (lo.y, hi.y, p.y), (lo.z, hi.z, p.z)] {
                let out = max2(l - x, max2(0.0, x - h));
                bx += out * out;
            }
            if f64::sqrt(bx) >= best {
                continue;
            }
            if let Ok((d, q)) = face_distance(k, b, f, p) {
                if d.is_finite() && d < best {
                    best = d;
                    at = Some(q);
                }
            }
        }
        (best, at)
    }
}

/// `leaves_cavity`: whether an inward offset of `wall` (> 0) from every face of the solid
/// `b` but `open` leaves a cavity: a point inside `b` farther than `wall`, and a margin,
/// from every face that is offset. Certified, never guessed: a point counts only as the
/// centre of a ball inside `b`, reached from a face sample by stepping inward along the
/// normal, then from ball to ball. SMLib's offset fails where no cavity remains
/// (tangency-labs/solidcore#380), so a wall this finds no room for is refused first.
#[allow(clippy::neg_cmp_op_on_partial_ord)] // the C API's NaN-safe tests
fn leaves_cavity(k: &mut Kernel, b: &mut Brep, open: &[FaceId], wall: f64) -> bool {
    const GRID: i32 = 5;
    const MAX_SAMPLES: usize = 2048;
    const STEPS: i32 = 48;
    let Ok(list) = solidcore::topology::brep::p1::faces(b, None) else {
        return false;
    };
    let mut all = FaceSet::new();
    let mut walls = FaceSet::new();
    for &f in &list {
        if !all.add(k, b, f) {
            return false;
        }
        if !open.contains(&f) && !walls.add(k, b, f) {
            return false;
        }
    }
    if walls.faces.is_empty() {
        return true;
    }
    let tol = b.tolerance();
    let margin = max2(10.0 * tol, 1e-3 * wall);
    let need = wall + margin;
    // Face samples: a grid over each offset face's domain, the points inside its trims,
    // each with the inward unit normal there.
    let mut samples: Vec<(Point3, Vec3)> = Vec::new();
    for &f in &walls.faces {
        let Ok(Some(fu)) = outer_faceuse(b, f) else {
            continue;
        };
        let (Ok(face), Ok(use_)) = (b.face(f), b.faceuse(fu)) else {
            continue;
        };
        if face.surface.is_none() {
            continue;
        }
        let inward = if use_.orientation == OrientType::Same {
            -1.0
        } else {
            1.0
        };
        let dom = face.uv_domain;
        for i in 0..GRID {
            if samples.len() >= MAX_SAMPLES {
                break;
            }
            for j in 0..GRID {
                if samples.len() >= MAX_SAMPLES {
                    break;
                }
                let uv = dom.evaluate(
                    (f64::from(i) + 0.5) / f64::from(GRID),
                    (f64::from(j) + 0.5) / f64::from(GRID),
                );
                if !matches!(face_classify_uv(k, b, f, uv), Ok(PointClassType::Face)) {
                    continue;
                }
                let Some(surf) = b.face(f).ok().and_then(|x| x.surface.as_ref()) else {
                    continue;
                };
                let (Ok(at), Ok(n)) = (
                    surf.evaluate_point(k, uv),
                    surf.evaluate_normal(k, uv, Side::Left, Side::Left),
                ) else {
                    continue;
                };
                let len = n.length();
                if !(len > 0.0) || !len.is_finite() {
                    continue;
                }
                samples.push((at, n * (inward / len)));
            }
        }
    }
    // First the direct step of `need` from each sample.
    for &(s, n) in &samples {
        if let Some((c, _)) = ball(k, b, &all, s, n, need) {
            if walls.nearest(k, b, c).0 >= need {
                return true;
            }
        }
    }
    // Then climbing, from a small ball at each sample: step away from the nearest offset
    // face by most of the ball's radius (staying inside the ball, so inside the body)
    // while the distance from the offset faces grows.
    for &(s, n) in &samples {
        let Some((mut c, mut r)) = ball(k, b, &all, s, n, max2(margin, 1e3 * tol)) else {
            continue;
        };
        for _ in 0..STEPS {
            let (d, q) = walls.nearest(k, b, c);
            if d >= need {
                return true;
            }
            let Some(q) = q else {
                break;
            };
            if !(d > 0.0) || !d.is_finite() {
                break;
            }
            let next = c + (c - q) * (0.9 * r / d);
            let (rn, _) = all.nearest(k, b, next);
            if !(rn > 0.0) || !rn.is_finite() || walls.nearest(k, b, next).0 <= d * (1.0 + 1e-9) {
                break;
            }
            c = next;
            r = rn;
        }
    }
    false
}

/// A step of `h` inward from the sample `s` (inward normal `n`) ends at a ball's centre
/// when no face of `all` is nearer to it than `h` (the step crossed none, and the ball
/// touches only the sample's face): the centre and the ball's radius.
fn ball(
    k: &mut Kernel,
    b: &mut Brep,
    all: &FaceSet,
    s: Point3,
    n: Vec3,
    h: f64,
) -> Option<(Point3, f64)> {
    let c = s + n * h;
    let (r, _) = all.nearest(k, b, c);
    (r.is_finite() && r >= h * (1.0 - 1e-9)).then_some((c, r))
}

/// `outward_normal`: the outward normal of face `f` of `b` at `at` (a point on it): its
/// surface's normal there, turned with the outer faceuse. `None` where the point cannot be
/// placed on the surface.
fn outward_normal(k: &Kernel, b: &Brep, f: FaceId, at: Point3) -> Option<Vec3> {
    let surf = b.face(f).ok()?.surface.as_ref()?;
    let fu = outer_faceuse(b, f).ok()??;
    let (step, _) = surf
        .step_inversion(k, &surf.step_uv_domain(), at, 1e6, None)
        .ok()?;
    let nurb = surf.convert_uv_from_step_to_nurbs(k, step).ok()?;
    let mut n = surf.evaluate_normal(k, nurb, Side::Left, Side::Left).ok()?;
    if b.faceuse(fu).ok()?.orientation != OrientType::Same {
        n = -n;
    }
    Some(n)
}

/// `drafted_cone`: the cone a face on the cone or cylinder `cyl` about the pull (its
/// outward normal `n` at `inside`) becomes drafted by `angle` radians about its circle on
/// the neutral plane (through `origin`, normal `normal`), so that it leans in going along
/// `pull`, which runs along the axis: its generators turned by `angle` in their planes
/// through the axis, about that circle. With `R0` the radius on the neutral plane and `t`
/// the slope of the radius going along the pull (0 for a cylinder), the radius `h` along
/// the pull from that circle is `R0 + h tan(atan(t) - s angle)`, `s = 1` where the face
/// looks away from the axis (a pin), `-1` where it looks towards it (a hole). It is the
/// surface's own NURBS with each control point moved radially by the new radius over the
/// old at its height, so it keeps the surface's parameters and the face's UV trims stay
/// where they were; then the cone found in it (`CopyAndAddAnalytics`). `None` when the
/// pull does not run along the axis; `SMLIB_ERR_UNSUPPORTED` when the cone would close
/// within the face's span.
#[allow(clippy::neg_cmp_op_on_partial_ord)] // the C API's NaN-safe tests
fn drafted_cone(
    k: &Kernel,
    cyl: &Cone,
    n: Vec3,
    inside: Point3,
    origin: Point3,
    normal: Vec3,
    pull: Vec3,
    angle: f64,
) -> core::result::Result<Option<Surface>, Status> {
    let pos = cyl.rev().position();
    let Ok(z) = pos.z_axis().unitize() else {
        return Ok(None);
    };
    if (z.dot(pull).abs() - 1.0).abs() > 1e-9 {
        return Ok(None);
    }
    let axis = pos.canonical().0;
    let along = z.dot(normal);
    if along.abs() < 1e-12 {
        return Err(Status::ERR_UNSUPPORTED);
    }
    // The hinge circle's centre: the axis on the neutral plane, at the cone's height v0.
    let v0 = (origin - axis).dot(normal) / along;
    let centre = axis + z * v0;
    let radial = (inside - centre) - pull * (inside - centre).dot(pull);
    let s = if radial.dot(n) >= 0.0 { 1.0 } else { -1.0 };
    // The cone's radius is linear in its height v along z; along the pull it changes by
    // `slope` per unit, and drafted by `turned`.
    let r0 = cyl.radius(v0);
    let slope = (cyl.radius(1.0) - cyl.radius(0.0)) * z.dot(pull);
    let turned = (slope.atan() - s * angle).tan();
    if !(r0 > 0.0) || !turned.is_finite() {
        return Err(Status::ERR_UNSUPPORTED);
    }
    let mut net = cyl.rev().base().clone();
    let nu = net.number_control_points(SurfParamType::U);
    let nv = net.number_control_points(SurfParamType::V);
    for i in 0..nu {
        for j in 0..nv {
            let Ok((p, w)) = net.control_point(ControlPointForm::EuclideanRational, i, j) else {
                return Err(Status::ERR_UNSUPPORTED);
            };
            let h = (p - centre).dot(pull);
            let was = cyl.radius(v0 + h * z.dot(pull));
            let now = r0 + h * turned;
            if !(was > 0.0) || !(now > 0.0) {
                return Err(Status::ERR_UNSUPPORTED);
            }
            let foot = centre + pull * h;
            if net
                .set_control_point(
                    ControlPointForm::EuclideanRational,
                    i,
                    j,
                    foot + (p - foot) * (now / was),
                    w,
                )
                .is_err()
            {
                return Err(Status::ERR_UNSUPPORTED);
            }
        }
    }
    let Ok(cone) = Surface::from(net).copy_and_add_analytics(k) else {
        return Err(Status::ERR_UNSUPPORTED);
    };
    Ok(Some(cone))
}

/// `draft_faces`: tilts each face of `picked` by `angle` radians about its line on the
/// neutral plane (through `origin`, normal `normal`), so that it leans in going along
/// `pull`, then lets the local operation extend or trim the faces around it to meet the
/// tilted surfaces. A plane must be parallel to the pull; the hinge is where a line from
/// a point inside the face, along the pull, meets the neutral plane. A cylinder or cone
/// must run about the pull, and turns about its circle on the neutral plane
/// (`drafted_cone`).
#[allow(clippy::neg_cmp_op_on_partial_ord)] // the C API's NaN-safe test
fn draft_faces(
    k: &mut Kernel,
    b: &mut Brep,
    picked: &[FaceId],
    origin: Point3,
    normal: Vec3,
    pull: Vec3,
    angle: f64,
) -> core::result::Result<(), Fail> {
    let mut faces = Vec::new();
    let mut surfaces: Vec<Surface> = Vec::new();
    for &f in picked {
        let fu = outer_faceuse(b, f)?;
        let surf = b.face(f).ok().and_then(|x| x.surface.as_ref());
        let (Some(surf), Some(fu)) = (surf, fu) else {
            return Err(Status::SM_ERR_INVALID_INPUT.into());
        };
        // A cylinder or cone about the pull drafted becomes a cone (drafted_cone).
        if let SurfaceKind::Cone(cyl) = surf.kind() {
            let cyl = cyl.as_ref().clone();
            let Ok(inside) = face_internal_point(k, b, f) else {
                return Err(Status::SM_ERR_INVALID_INPUT.into());
            };
            let Some(n) = outward_normal(k, b, f, inside) else {
                return Err(Status::SM_ERR_INVALID_INPUT.into());
            };
            match drafted_cone(k, &cyl, n, inside, origin, normal, pull, angle) {
                Ok(Some(cone)) => {
                    faces.push(f);
                    surfaces.push(cone);
                    continue;
                }
                Ok(None) => return Err(Status::ERR_UNSUPPORTED.into()),
                Err(st) => return Err(st.into()),
            }
        }
        let SurfaceKind::Plane(plane) = surf.kind() else {
            return Err(Status::ERR_UNSUPPORTED.into());
        };
        // The face's outward normal: the plane's, turned with the outer faceuse.
        let mut n = plane.position().z_axis();
        if b.faceuse(fu)?.orientation != OrientType::Same {
            n = -n;
        }
        if n.dot(pull).abs() > 1e-6 {
            return Err(Status::ERR_UNSUPPORTED.into());
        }
        let Ok(inside) = face_internal_point(k, b, f) else {
            return Err(Status::SM_ERR_INVALID_INPUT.into());
        };
        let along = pull.dot(normal);
        if along.abs() < 1e-12 {
            return Err(Status::ERR_ARG.into());
        }
        let hinge = inside + pull * ((origin - inside).dot(normal) / along);
        // Turning n about n x pull by `angle` gives n cos(angle) + pull sin(angle).
        let mut axis = n.cross(pull);
        let len = axis.length();
        if !(len > 0.5) {
            return Err(Status::ERR_UNSUPPORTED.into());
        }
        axis /= len;
        let surf = b
            .face(f)
            .ok()
            .and_then(|x| x.surface.as_ref())
            .ok_or(Fail::Status(Status::SM_ERR_INVALID_INPUT))?;
        let Ok(mut tilted) = surf.copy(k) else {
            return Err(Status::SM_ERR_INVALID_INPUT.into());
        };
        let mut turn = Axis2Placement::IDENTITY;
        turn.rotate_about_axis_at_point(angle, hinge, axis);
        if tilted.transform(k, &turn, None).is_err() {
            return Err(Status::SM_ERR_INVALID_INPUT.into());
        }
        faces.push(f);
        surfaces.push(tilted);
    }
    let targets = surfaces.clone();
    local_operation_simple(k, b, &faces, surfaces)?;
    moved_faces(k, b, &targets)
}

/// Whether `p` lies on `target` within `tol` (planes and cones, the surfaces a draft
/// hands over).
fn on_surface(target: &Surface, p: Point3, tol: f64) -> bool {
    match target.kind() {
        SurfaceKind::Plane(plane) => {
            let o = plane.position().canonical().0;
            (p - o).dot(plane.position().z_axis()).abs() <= tol
        }
        SurfaceKind::Cone(cone) => {
            let pos = cone.rev().position();
            let Ok(z) = pos.z_axis().unitize() else {
                return false;
            };
            let o = pos.canonical().0;
            let h = (p - o).dot(z);
            let r = (p - o - z * h).length();
            (r - cone.radius(h)).abs() <= tol
        }
        _ => false,
    }
}

/// `moved_faces`: after a local operation reported success, some face of `b` must lie on
/// one of the `targets` it was handed, or it moved none (`SMLIB_ERR_NO_EFFECT`).
fn moved_faces(k: &Kernel, b: &Brep, targets: &[Surface]) -> core::result::Result<(), Fail> {
    let tol = b.tolerance();
    for f in solidcore::topology::brep::p1::faces(b, None)? {
        let face = b.face(f)?;
        let Some(surf) = face.surface.as_ref() else {
            continue;
        };
        let Ok(p) = surf.evaluate_point(k, face.uv_domain.evaluate(0.5, 0.5)) else {
            continue;
        };
        if targets.iter().any(|t| on_surface(t, p, tol)) {
            return Ok(());
        }
    }
    Err(Status::ERR_NO_EFFECT.into())
}

/// A unit vector from `v`, or `None` when it has no usable length (`smlib_draft`'s
/// `unit`).
#[allow(clippy::neg_cmp_op_on_partial_ord)] // the C API's NaN-safe test
fn unit(v: [f64; 3]) -> Option<Vec3> {
    if !finite(&v) {
        return None;
    }
    let u = Vec3::new(v[0], v[1], v[2]);
    let len = u.length();
    if !(len > 1e-12) || !len.is_finite() {
        return None;
    }
    Some(u / len)
}

impl Smlib {
    /// `smlib_shell`: `brep` hollowed to walls `thickness` thick (negative: inward), open
    /// at the faces the `points` pick (none: a closed cavity), each point within
    /// `tolerance` of its face. An inward wall that leaves no cavity is
    /// `SMLIB_ERR_NO_ROOM`, refused before the offset runs.
    pub fn shell(
        &mut self,
        brep: &Brep,
        points: &[[f64; 3]],
        tolerance: f64,
        thickness: f64,
    ) -> Result<Brep> {
        if !points_ok(points, true)
            || !tolerance_ok(tolerance)
            || !thickness.is_finite()
            || thickness == 0.0
            || thickness.abs() > 1e6
        {
            return Err(Status::ERR_ARG);
        }
        self.run(
            || {
                format!(
                    "smlib_shell(brep, points={}, tolerance={}, thickness={})",
                    pts(points),
                    num(tolerance),
                    num(thickness)
                )
            },
            &[("brep", brep)],
            |api| {
                let mut w = api.copy(brep)?;
                let picked = picked_faces(api, &mut w, points, tolerance)?;
                let solid_in = brep_is_manifold_solid(&w)?;
                // An inward wall with no room left inside it corrupts SMLib's working body
                // (tangency-labs/solidcore#380): refuse it before the offset sees it.
                if thickness < 0.0
                    && solid_in
                    && !leaves_cavity(api.kernel_mut(), &mut w, &picked, -thickness)
                {
                    return Err(Status::ERR_NO_ROOM.into());
                }
                let r = shell_brep_simple(api.kernel_mut(), &w, shell_options(thickness), &picked)?;
                finish(r, solid_in, false)
            },
        )
    }

    /// `smlib_thicken`: a sheet thickened into a solid `thickness` thick along its normals
    /// (negative: against them); `symmetric`: half on each side. A solid is
    /// `SMLIB_ERR_ARG`.
    pub fn thicken(&mut self, brep: &Brep, thickness: f64, symmetric: bool) -> Result<Brep> {
        if !thickness.is_finite() || thickness == 0.0 || thickness.abs() > 1e6 {
            return Err(Status::ERR_ARG);
        }
        self.run(
            || {
                format!(
                    "smlib_thicken(brep, thickness={}, symmetric={})",
                    num(thickness),
                    i32::from(symmetric)
                )
            },
            &[("brep", brep)],
            |api| {
                let mut w = api.copy(brep)?;
                if brep_is_manifold_solid(&w)? {
                    return Err(Status::ERR_ARG.into());
                }
                // Symmetric: the sheet moved back by half first, then thickened by the whole.
                if symmetric {
                    w = offset_brep_simple(
                        api.kernel_mut(),
                        &mut w,
                        offset_options(-0.5 * thickness),
                    )?;
                }
                let r = shell_brep_simple(api.kernel_mut(), &w, shell_options(thickness), &[])?;
                finish(r, true, false)
            },
        )
    }

    /// `smlib_draft`: each face the `points` pick (planes parallel to `pull`) tilted
    /// `angle` radians (`0 < |angle| < 89.94` degrees) about its line on the neutral plane
    /// through `neutral_origin` with normal `neutral_normal`, leaning in going along
    /// `pull`; the faces around them extended or trimmed to meet. A curved or
    /// non-parallel face is `SMLIB_ERR_UNSUPPORTED`; the body must be a solid.
    pub fn draft(
        &mut self,
        brep: &Brep,
        points: &[[f64; 3]],
        tolerance: f64,
        neutral_origin: [f64; 3],
        neutral_normal: [f64; 3],
        pull: [f64; 3],
        angle: f64,
    ) -> Result<Brep> {
        let (Some(normal), Some(along)) = (unit(neutral_normal), unit(pull)) else {
            return Err(Status::ERR_ARG);
        };
        if !points_ok(points, false)
            || !tolerance_ok(tolerance)
            || !finite(&neutral_origin)
            || neutral_origin.iter().any(|x| x.abs() > 1e9)
            || !angle.is_finite()
            || angle == 0.0
            || angle.abs() >= 1.5697963267948966
        {
            return Err(Status::ERR_ARG);
        }
        let origin = Point3::new(neutral_origin[0], neutral_origin[1], neutral_origin[2]);
        self.run(
            || {
                format!(
                    "smlib_draft(brep, points={}, tolerance={}, neutral={}/{}, pull={}, angle={})",
                    pts(points),
                    num(tolerance),
                    vec(&neutral_origin),
                    vec(&neutral_normal),
                    vec(&pull),
                    num(angle)
                )
            },
            &[("brep", brep)],
            |api| {
                let mut w = api.copy(brep)?;
                if !brep_is_manifold_solid(&w)? {
                    return Err(Status::ERR_ARG.into());
                }
                let picked = picked_faces(api, &mut w, points, tolerance)?;
                draft_faces(
                    api.kernel_mut(),
                    &mut w,
                    &picked,
                    origin,
                    normal,
                    along,
                    angle,
                )?;
                finish(w, true, false)
            },
        )
    }

    /// `smlib_offset_sheet`: every face of a sheet, or of faces that close on themselves
    /// (a sphere's one face), moved `distance` along its normal (`SmApiOffsetBrepFull`).
    pub fn offset_sheet(&mut self, sheet: &Brep, distance: f64) -> Result<Brep> {
        if !distance.is_finite() || distance == 0.0 || distance.abs() > 1e6 {
            return Err(Status::ERR_ARG);
        }
        self.run(
            || format!("smlib_offset_sheet(sheet, distance={})", num(distance)),
            &[("sheet", sheet)],
            |api| {
                let mut w = api.copy(sheet)?;
                let moved = offset_brep_simple(api.kernel_mut(), &mut w, offset_options(distance))?;
                take(moved)
            },
        )
    }

    /// `smlib_merge_faces`: the body with every edge and vertex that separates nothing
    /// removed (faces on one surface merged), but the edges the `keep` points pick, each
    /// within `tolerance`.
    pub fn merge_faces(&mut self, brep: &Brep, keep: &[[f64; 3]], tolerance: f64) -> Result<Brep> {
        if !points_ok(keep, true) || !tolerance_ok(tolerance) {
            return Err(Status::ERR_ARG);
        }
        self.run(
            || {
                format!(
                    "smlib_merge_faces(brep, keep={}, tolerance={})",
                    pts(keep),
                    num(tolerance)
                )
            },
            &[("brep", brep)],
            |api| {
                let mut w = api.copy(brep)?;
                let kept = if keep.is_empty() {
                    Vec::new()
                } else {
                    picked_edges(api, &w, keep, tolerance)?
                };
                let solid_in = brep_is_manifold_solid(&w)?;
                let keep_edges = (!kept.is_empty()).then_some(kept.as_slice());
                if let Err(e) = remove_topological_edges_and_vertices(
                    api.kernel_mut(),
                    &mut w,
                    keep_edges,
                    None,
                ) {
                    // the SMLib status as is (`static_cast<int32_t>(r)`)
                    return Err(match e.sm_status() {
                        Some(code) => Fail::Status(Status::from_code(code)),
                        None => Fail::Solid(e),
                    });
                }
                finish(w, solid_in, false)
            },
        )
    }
}

/// `pick_faces` on the working copy `b`, the failed point's index recorded.
fn picked_faces(
    api: &mut Smlib,
    b: &mut Brep,
    points: &[[f64; 3]],
    tolerance: f64,
) -> core::result::Result<Vec<FaceId>, Fail> {
    match pick_faces(api.kernel_mut(), b, points, tolerance) {
        Ok(picked) => Ok(picked.items),
        Err(fail) => {
            if let Some(i) = fail.index() {
                api.set_failed_point(i);
            }
            Err(fail.into())
        }
    }
}

/// `pick_edges` on `b`, the failed point's index recorded.
fn picked_edges(
    api: &mut Smlib,
    b: &Brep,
    points: &[[f64; 3]],
    tolerance: f64,
) -> core::result::Result<Vec<EdgeId>, Fail> {
    match pick_edges(api.kernel(), b, points, tolerance) {
        Ok(picked) => Ok(picked.items),
        Err(fail) => {
            if let Some(i) = fail.index() {
                api.set_failed_point(i);
            }
            Err(fail.into())
        }
    }
}
