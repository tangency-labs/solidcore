//! Extend Sheet: a one-face sheet bounded by its surface's parameter lines, continued past
//! one of those sides by SMLib's extension of the surface
//! (`SmBSplineSurface::CreateExtendedSurface`: its end spans continued, or a G1 or G2
//! run-off). SMLib holds every surface, its analytic ones included, as a B-spline, so one
//! path serves all: the extended B-spline continues a plane, a circle or any conic section
//! on itself, and `CopyAndAddAnalytics` gives it back its analytic kind where it is one.
//! The new face, and each of its sides, takes the attributes of the face and side it
//! continues (history tags).
//!
//! Source: tangency-labs/smlib `capi/src/smlib_extend.cpp` (Apache-2.0, Copyright
//! SolidCraft contributors) at fork commit b3dd0e5 (replay/1.0.1).
//!
//! | Rust | C API (SMLib call) |
//! |------|-------|
//! | [`Smlib::extend_sheet`] | `smlib_extend_sheet` (`CreateExtendedSurface`, `CopyAndAddAnalytics`) |
//! | `Side`, `SIDES`, `on_side`, `derivatives_at` / `point_at` (`eval`), `across`, `cone_continued`, `same_patch`, `edge_distance`, `nearest_edge`, `is_boundary`, `side_edge`, `copy_attributes` (private) | the same |
//!
//! Deviations: `SmAttribute::MakeCopy` is [`Attribute::make_copy`] onto the new body's
//! object.

use solidcore::geom::attribute::Attribute;
use solidcore::geom::extent2d::Extent2d;
use solidcore::geom::geom_utility::p2::ContinuityType;
use solidcore::geom::object::ObjectRef;
use solidcore::geom::tol::USE_DEFAULT;
use solidcore::geom::vector2d::Point2;
use solidcore::geom::vector3d::Vec3;
use solidcore::nurbs::b_spline_surface::{BSplineSurface, SurfParamType};
use solidcore::nurbs::b_spline_surface_n_lib::p2::EvaluateFlags;
use solidcore::nurbs::knots::Side as KnotSide;
use solidcore::query::api_queries::edge_closest_point;
use solidcore::surface::surface::{Surface as KSurface, SurfaceKind};
use solidcore::topology::brep::p1 as bp1;
use solidcore::topology::brep::p2::ArenaNotifyState;
use solidcore::topology::brep::p3::new_brep;
use solidcore::topology::brep::p8::{create_face_from_surface, FaceBuildEnv};
use solidcore::topology::edgeuse::p1::faceuse_orientation;
use solidcore::topology::face::p1::flip_face_orientation;
use solidcore::topology::face::{p2 as face2, p4 as face4};
use solidcore::topology::topology_traverser::Marks;
use solidcore::topology::{Brep, EdgeId, FaceId};
use solidcore::{EntityRef, Kernel};

use crate::session::{finite, num, take, vec, Fail, Smlib};
use crate::status::{Result, Status};

/// `SMLIB_EXTEND_NATURAL`: the end spans continued (`SM_CT_CINFINITY`).
pub const EXTEND_NATURAL: i32 = 0;
/// `SMLIB_EXTEND_TANGENT`: a G1 run-off (`SM_CT_G1`).
pub const EXTEND_TANGENT: i32 = 1;
/// `SMLIB_EXTEND_CURVATURE`: a G2 run-off (`SM_CT_G1_G2`).
pub const EXTEND_CURVATURE: i32 = 2;

#[derive(Debug, Clone, Copy)]
struct Side {
    which: SurfParamType,
    /// The side is a u = const line (the face grows in u).
    u: bool,
    /// At the box's max.
    high: bool,
}

const SIDES: [Side; 4] = [
    Side {
        which: SurfParamType::UMin,
        u: true,
        high: false,
    },
    Side {
        which: SurfParamType::UMax,
        u: true,
        high: true,
    },
    Side {
        which: SurfParamType::VMin,
        u: false,
        high: false,
    },
    Side {
        which: SurfParamType::VMax,
        u: false,
        high: true,
    },
];

/// The point of side `s` of `dom` at `t` in `[0, 1]` along it, in UV.
fn on_side(dom: &Extent2d, s: Side, t: f64) -> Point2 {
    let (lo, hi) = (dom.min(), dom.max());
    if s.u {
        return Point2::new(if s.high { hi.x } else { lo.x }, lo.y + t * (hi.y - lo.y));
    }
    Point2::new(lo.x + t * (hi.x - lo.x), if s.high { hi.y } else { lo.y })
}

/// `eval`, `surf->Evaluate(uv, 1, 1, TRUE, TRUE, TRUE, d)`: the point and its first
/// derivatives in u and v, `None` unless all three are finite.
fn derivatives_at(k: &Kernel, surf: &KSurface, uv: Point2) -> Option<(Vec3, Vec3, Vec3)> {
    let flags = EvaluateFlags {
        u_side: KnotSide::Left,
        v_side: KnotSide::Left,
        only_upper_half: true,
        non_zero_tangents: true,
    };
    let d = surf.evaluate(k, uv, 1, 1, flags).ok()?;
    let p = *d.first()?.first()?;
    let dv = *d.first()?.get(1)?;
    let du = *d.get(1)?.first()?;
    let ok = p.length().is_finite() && dv.length().is_finite() && du.length().is_finite();
    ok.then_some((p, du, dv))
}

/// `eval` for the point only.
fn point_at(k: &Kernel, surf: &KSurface, uv: Point2) -> Option<Vec3> {
    derivatives_at(k, surf, uv).map(|(p, _, _)| p)
}

/// The 16-point Gauss-Legendre abscissae (positive half).
const GAUSS_X: [f64; 8] = [
    0.095_012_509_837_637_4,
    0.281_603_550_779_258_9,
    0.458_016_777_657_227_4,
    0.617_876_244_402_643_8,
    0.755_404_408_355_003,
    0.865_631_202_387_831_8,
    0.944_575_023_073_232_6,
    0.989_400_934_991_649_9,
];
/// Their weights.
const GAUSS_W: [f64; 8] = [
    0.189_450_610_455_068_5,
    0.182_603_415_044_923_6,
    0.169_156_519_395_002_5,
    0.149_595_988_816_576_7,
    0.124_628_971_255_533_9,
    0.095_158_511_682_492_8,
    0.062_253_523_938_647_9,
    0.027_152_459_411_754_1,
];

/// `across`: how far the surface `ext` (the face's surface extended past side `s` of its
/// box `dom`, on the same parameters) runs past that side, measured across it: the length
/// of the parameter line through the side's middle from the side to `ext`'s own bound, by
/// 16-point Gauss-Legendre quadrature on 64 pieces. Negative when it can't be evaluated.
fn across(k: &Kernel, ext: &KSurface, dom: &Extent2d, s: Side) -> f64 {
    const PIECES: u32 = 64;
    let Ok(bx) = ext.natural_uv_domain_with_nurb(k) else {
        return -1.0;
    };
    let mid = on_side(dom, s, 0.5);
    let from = if s.u { mid.x } else { mid.y };
    let to = match (s.u, s.high) {
        (true, true) => bx.max().x,
        (true, false) => bx.min().x,
        (false, true) => bx.max().y,
        (false, false) => bx.min().y,
    };
    let h = (to - from) / f64::from(PIECES);
    let mut length = 0.0;
    for piece in 0..PIECES {
        let c = from + (f64::from(piece) + 0.5) * h;
        // (abscissae -X[0..8] then X[0..8], as the C API's `i < 8 ? -X[i] : X[i - 8]`)
        for (x, w) in GAUSS_X
            .iter()
            .zip(GAUSS_W)
            .map(|(&x, w)| (-x, w))
            .chain(GAUSS_X.iter().zip(GAUSS_W).map(|(&x, w)| (x, w)))
        {
            let t = c + 0.5 * h * x;
            let uv = if s.u {
                Point2::new(t, mid.y)
            } else {
                Point2::new(mid.x, t)
            };
            let Some((_, du, dv)) = derivatives_at(k, ext, uv) else {
                return -1.0;
            };
            length += 0.5 * h.abs() * w * (if s.u { du } else { dv }).length();
        }
    }
    length
}

/// `cone_continued`: a cone's or cylinder's face continued naturally `distance` past side
/// `s` of its box `dom`, as the analytic itself (the cone over a wider STEP domain: degrees
/// round the axis, height along it): round the axis by the arc of `distance` at the side's
/// middle, along it by `distance` on the slant. `None` where the cone can't be widened so
/// (more than a turn, or past the apex).
fn cone_continued(
    k: &Kernel,
    surf: &KSurface,
    dom: &Extent2d,
    s: Side,
    distance: f64,
) -> Option<KSurface> {
    let SurfaceKind::Cone(cone) = surf.kind_with_nurb(k).ok()? else {
        return None;
    };
    let mut cone = cone.as_ref().clone();
    let a = cone
        .convert_uv_from_nurbs_to_step(k, on_side(dom, s, 0.0))
        .ok()?;
    let b = cone
        .convert_uv_from_nurbs_to_step(k, on_side(dom, s, 1.0))
        .ok()?;
    let m = cone
        .convert_uv_from_nurbs_to_step(k, on_side(dom, s, 0.5))
        .ok()?;
    let lo = cone.convert_uv_from_nurbs_to_step(k, dom.min()).ok()?;
    let hi = cone.convert_uv_from_nurbs_to_step(k, dom.max()).ok()?;
    let semi = cone.canonical().semi_angle_deg;
    let mut bx = [
        [lo.x.min(hi.x), lo.x.max(hi.x)],
        [lo.y.min(hi.y), lo.y.max(hi.y)],
    ];
    let eps = 1e-9 * (1.0 + m.x.abs() + m.y.abs());
    // The STEP coordinate constant along the side is the one the face grows in.
    let along_u = if (a.x - b.x).abs() <= eps {
        true
    } else if (a.y - b.y).abs() <= eps {
        false
    } else {
        return None;
    };
    let at = if along_u { m.x } else { m.y };
    let row = if along_u { &mut bx[0] } else { &mut bx[1] };
    let high = (at - row[1]).abs() <= (at - row[0]).abs();
    let by = if along_u {
        let r = cone.radius(m.y);
        let positive = r > 0.0;
        if !positive {
            return None;
        }
        distance / r * 180.0 / core::f64::consts::PI
    } else {
        distance * (semi * core::f64::consts::PI / 180.0).cos()
    };
    if high {
        row[1] += by;
    } else {
        row[0] -= by;
    }
    if bx[0][1] - bx[0][0] > 360.0 + 1e-9 {
        return None;
    }
    let wider = Extent2d::from_bounds(bx[0][0], bx[1][0], bx[0][1], bx[1][1]);
    cone.adjust_step_uv_domain(k, &wider).ok()?;
    Some(KSurface::from(cone))
}

/// `same_patch`: whether surfaces `p` and `q` span the same patch: their boxes' corners
/// the same points, in any order, within `tol`.
fn same_patch(k: &Kernel, p: &KSurface, q: &KSurface, tol: f64) -> bool {
    let (Ok(bp), Ok(bq)) = (
        p.natural_uv_domain_with_nurb(k),
        q.natural_uv_domain_with_nurb(k),
    ) else {
        return false;
    };
    let corners = |b: &Extent2d| {
        [
            b.min(),
            Point2::new(b.max().x, b.min().y),
            b.max(),
            Point2::new(b.min().x, b.max().y),
        ]
    };
    for u in corners(&bp) {
        let Some(x) = point_at(k, p, u) else {
            return false;
        };
        let found = corners(&bq)
            .iter()
            .any(|&v| point_at(k, q, v).is_some_and(|y| (x - y).length() <= tol));
        if !found {
            return false;
        }
    }
    true
}

fn edge_distance(k: &Kernel, b: &Brep, e: EdgeId, p: Vec3) -> f64 {
    match edge_closest_point(k, b, e, p) {
        Ok(c) if c.closest.distance.is_finite() => c.closest.distance,
        _ => f64::INFINITY,
    }
}

fn nearest_edge(k: &Kernel, b: &Brep, edges: &[EdgeId], p: Vec3) -> (f64, Option<EdgeId>) {
    let mut best = f64::INFINITY;
    let mut which = None;
    for &e in edges {
        let d = edge_distance(k, b, e, p);
        if d < best {
            best = d;
            which = Some(e);
        }
    }
    (best, which)
}

/// The samples along a side.
const SAMPLES: u32 = 9;

fn sample_t(i: u32) -> f64 {
    f64::from(i) / f64::from(SAMPLES)
}

/// Whether side `s` of the face's box is a side of the face: its samples lie on the
/// face's edges (`Some(true)` for a side that collapses to a point, a pole; `None` when
/// not a side).
fn is_boundary(
    k: &Kernel,
    b: &Brep,
    surf: &KSurface,
    dom: &Extent2d,
    s: Side,
    edges: &[EdgeId],
    tol: f64,
) -> Option<bool> {
    let a = point_at(k, surf, on_side(dom, s, 0.0))?;
    let e = point_at(k, surf, on_side(dom, s, 1.0))?;
    let m = point_at(k, surf, on_side(dom, s, 0.5))?;
    if (a - e).length() <= tol && (a - m).length() <= tol {
        return Some(true);
    }
    for i in 0..=SAMPLES {
        let p = point_at(k, surf, on_side(dom, s, sample_t(i)))?;
        if nearest_edge(k, b, edges, p).0 > tol {
            return None;
        }
    }
    Some(false)
}

/// The edge of `edges` that side `s` of `dom` runs along (nearest the side's middle), or
/// `None` for a side collapsed to a point.
fn side_edge(
    k: &Kernel,
    b: &Brep,
    surf: &KSurface,
    dom: &Extent2d,
    s: Side,
    edges: &[EdgeId],
    on: f64,
) -> Option<EdgeId> {
    let a = point_at(k, surf, on_side(dom, s, 0.0))?;
    let e = point_at(k, surf, on_side(dom, s, 1.0))?;
    let m = point_at(k, surf, on_side(dom, s, 0.5))?;
    if (a - e).length() <= on && (a - m).length() <= on {
        return None;
    }
    let (d, which) = nearest_edge(k, b, edges, m);
    if d <= on {
        which
    } else {
        None
    }
}

/// `to` (of `rb`) takes a copy of each attribute of `from` (of `b`), as SMLib's copies
/// pass them on.
fn copy_attributes(
    b: &Brep,
    from: EntityRef,
    rb: &mut Brep,
    to: EntityRef,
) -> solidcore::Result<()> {
    let copies: Vec<Attribute> = b
        .entity_attributes(from)
        .into_iter()
        .map(|(_, a)| a.make_copy())
        .collect();
    for a in copies {
        rb.add_attribute(to, a)?;
    }
    Ok(())
}

/// `SM_CAST_PTR(SmBSplineSurface, pSurface)` of a surface kind with its NURB made.
fn b_spline_part(kind: &SurfaceKind) -> Option<&BSplineSurface> {
    match kind {
        SurfaceKind::BSpline(b) => Some(b),
        SurfaceKind::Plane(p) => Some(p.base()),
        SurfaceKind::SurfOfRevolution(s) => Some(s.base()),
        SurfaceKind::SurfOfExtrusion(s) => Some(s.base()),
        SurfaceKind::Cone(c) => Some(c.rev().base()),
        SurfaceKind::Sphere(s) => Some(s.base().base()),
        SurfaceKind::Torus(t) => Some(t.base().base()),
        _ => None,
    }
}

impl Smlib {
    /// `smlib_extend_sheet`: a new one-face sheet, `sheet`'s face (all of its surface,
    /// bounded by its parameter lines) continued by `distance` past the side whose edge
    /// lies within `tol` of `point`, by `mode` ([`EXTEND_NATURAL`], [`EXTEND_TANGENT`] or
    /// [`EXTEND_CURVATURE`]). The new face and its sides keep the attributes of the face
    /// and sides they continue, and the face faces as the old one did.
    ///
    /// # Errors
    ///
    /// `SMLIB_ERR_ARG` for a non-finite point, a tolerance or distance not above 0, or
    /// another mode; `SMLIB_ERR_UNSUPPORTED` for a sheet of more than one face, a face with
    /// holes, one that is not all of its surface or not bounded by its parameter lines;
    /// `SMLIB_ERR_PICK_NONE` when no edge lies within `tol` of `point`; `SMLIB_ERR_EMPTY`
    /// for an empty result; the kernel's status otherwise.
    pub fn extend_sheet(
        &mut self,
        sheet: &Brep,
        point: [f64; 3],
        tol: f64,
        distance: f64,
        mode: i32,
    ) -> Result<Brep> {
        let ok = finite(&point)
            && tol.is_finite()
            && tol > 0.0
            && distance.is_finite()
            && distance > 0.0
            && (EXTEND_NATURAL..=EXTEND_CURVATURE).contains(&mode);
        if !ok {
            return Err(Status::ERR_ARG);
        }
        self.run(
            || {
                format!(
                    "smlib_extend_sheet(sheet, {}, {}, {}, {mode})",
                    vec(&point),
                    num(tol),
                    num(distance)
                )
            },
            &[("sheet", sheet)],
            |api| extend(api.kernel_mut(), sheet, point, tol, distance, mode),
        )
    }
}

fn empty() -> Fail {
    Fail::Status(Status::ERR_EMPTY)
}

const fn unsupported() -> Fail {
    Fail::Status(Status::ERR_UNSUPPORTED)
}

fn extend(
    k: &mut Kernel,
    b: &Brep,
    point: [f64; 3],
    tol: f64,
    distance: f64,
    mode: i32,
) -> core::result::Result<Brep, Fail> {
    let faces = bp1::faces(b, None)?;
    let &[f] = faces.as_slice() else {
        return Err(unsupported());
    };
    let Some(surf) = face2::surface(b, f)? else {
        return Err(unsupported());
    };
    if face2::loops(b, f, false, None)?.len() != 1 {
        return Err(unsupported());
    }
    let edges = face4::edges(b, f, None, false, None)?;
    let at = Vec3::new(point[0], point[1], point[2]);
    let (d, picked) = nearest_edge(k, b, &edges, at);
    let picked = match picked {
        Some(e) if d <= tol => e,
        _ => return Err(Fail::Status(Status::ERR_PICK_NONE)),
    };
    let model_tol = b.tolerance().max(1e-9);
    let on = 10.0 * model_tol;
    let dom = face2::uv_domain(b, f)?;
    // Every side of the box is a side of the face (the face is its whole box), and the
    // picked edge is all of one side.
    let mut grow = None;
    for s in SIDES {
        let Some(collapsed) = is_boundary(k, b, surf, &dom, s, &edges, on) else {
            return Err(unsupported());
        };
        if collapsed {
            continue;
        }
        let all = (0..=SAMPLES).all(|i| {
            point_at(k, surf, on_side(&dom, s, sample_t(i)))
                .is_some_and(|p| edge_distance(k, b, picked, p) <= on)
        });
        if all {
            grow = Some(s);
        }
    }
    let Some(grow) = grow else {
        return Err(unsupported());
    };
    // SMLib extends the whole surface; the face must be all of it.
    let natural = surf.natural_uv_domain_with_nurb(k)?;
    let (nmin, nmax, dmin, dmax) = (natural.min(), natural.max(), dom.min(), dom.max());
    let ptol = 1e-9 * 1.0_f64.max((nmax.x - nmin.x).max(nmax.y - nmin.y));
    let whole = (nmin.x - dmin.x).abs() <= ptol
        && (nmin.y - dmin.y).abs() <= ptol
        && (nmax.x - dmax.x).abs() <= ptol
        && (nmax.y - dmax.y).abs() <= ptol;
    if !whole {
        return Err(unsupported());
    }
    let ct = match mode {
        EXTEND_TANGENT => ContinuityType::G1,
        EXTEND_CURVATURE => ContinuityType::G1G2,
        _ => ContinuityType::CInfinity,
    };
    // CreateExtendedSurface takes its distance as an estimate (a natural parabola falls
    // 2.5% short): the distance asked of it is corrected until the extension measures
    // `distance` across the side's middle, keeping the surface's parameters on the face.
    let base = b_spline_part(surf.kind_with_nurb(k)?).ok_or_else(unsupported)?;
    let mut extended: Option<KSurface> = None;
    let mut last_err = None;
    let mut ask = distance;
    for _round in 0..12 {
        let tried = match base.create_extended_surface_dir(k, grow.which, ask, ct) {
            Ok(t) => KSurface::from(t),
            Err(e) => {
                last_err = Some(e);
                break;
            }
        };
        let ext = extended.insert(tried);
        // The face's side stays where it was on the extended surface's parameters.
        let mid = on_side(&dom, grow, 0.5);
        let (Some(was), Some(now)) = (point_at(k, surf, mid), point_at(k, ext, mid)) else {
            break;
        };
        if (was - now).length() > model_tol {
            break;
        }
        let got = across(k, ext, &dom, grow);
        let measured = got > 0.0;
        if !measured || (got - distance).abs() <= 1e-9 * distance {
            break;
        }
        ask *= distance / got;
    }
    let Some(extended) = extended else {
        return Err(last_err.map_or_else(empty, Fail::from));
    };
    // A cone or cylinder continued naturally stays one, as the analytic widened (its patch
    // the extended B-spline's), else SMLib gives it back its analytic kind where it can.
    let continued = if mode == EXTEND_NATURAL {
        cone_continued(k, surf, &dom, grow, distance)
    } else {
        None
    }
    .filter(|g| same_patch(k, g, &extended, 10.0 * model_tol));
    let grown = match continued {
        Some(g) => g,
        None => extended.copy_and_add_analytics(k)?,
    };
    let bx = grown.natural_uv_domain_with_nurb(k)?;
    let mut r = new_brep(k, USE_DEFAULT, USE_DEFAULT)?;
    r.set_tolerance(b.tolerance());
    let nf = {
        let mut state = ArenaNotifyState::new(&r);
        let mut env = state.env(None);
        let key_of: &dyn Fn(EntityRef) -> ObjectRef = env.key_of;
        let this_face = |x: FaceId| key_of(EntityRef::Face(x));
        let mut marks = Marks::default();
        let mut fenv = FaceBuildEnv {
            env: &mut env,
            marks: &mut marks,
            this_face: &this_face,
        };
        create_face_from_surface(k, &mut r, Some(grown.clone()), &bx, &mut fenv)?
    };
    // The face and each side stand for the face and side they continue.
    copy_attributes(b, EntityRef::Face(f), &mut r, EntityRef::Face(nf))?;
    let made = face4::edges(&r, nf, None, false, None)?;
    let new_on = 10.0 * r.tolerance().max(1e-9);
    for side in SIDES {
        let from = side_edge(k, b, surf, &dom, side, &edges, on);
        let to = side_edge(k, &r, &grown, &bx, side, &made, new_on);
        if let (Some(from), Some(to)) = (from, to) {
            copy_attributes(b, EntityRef::Edge(from), &mut r, EntityRef::Edge(to))?;
        }
    }
    // The new face faces the way the old one did.
    let was = face2::upward_faceuse(b, f)?;
    let now = face2::upward_faceuse(&r, nf)?;
    if let (Some(was), Some(now)) = (was, now) {
        if faceuse_orientation(b, was)? != faceuse_orientation(&r, now)? {
            flip_face_orientation(&mut r, nf)?;
        }
    }
    take(r)
}
