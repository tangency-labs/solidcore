//! A procedural blend to its record and back (#2388): [`to_record`] writes the definition
//! and the net as plain data, [`from_record`] makes the blend from it exactly, and
//! [`rebuild`] also re-certifies the net against the definition
//! ([`ProceduralBlend::recertify`]), refusing a net that is off it. This is the record ↔
//! kernel conversion `docs/layers.md` gives `solidcore-ops`: it copies numbers and computes
//! nothing, and moves there when that crate exists.
//!
//! Not upstream. The fork carries a definition through its database only
//! (`SmProceduralBlend::WriteToDB`, [`ProceduralBlend::read_from_db`]).
//!
//! Deviations: none (no upstream counterpart). Recursion: nested supports are walked at
//! most [`BLEND_MAX_DEPTH`] deep, with an explicit depth refused past it.

use super::record::{
    BlendDefinitionRecord, NurbsCurveRecord, NurbsSurfaceRecord, ParameterMapRecord,
    ProceduralBlendRecord, RadiusRecord, SupportKindRecord, SupportRecord,
};
use super::{BlendDefinitionInput, ProceduralBlend, BLEND_MAX_DEPTH};
use crate::blend::fillet_section::section::GnSectionFsg;
use crate::blend::fillet_section::station::as_b_spline;
use crate::error::{EntityRef, Result, SolidError, SM_ERR_INVALID_INPUT};
use crate::kernel::Kernel;
use crate::nurbs::b_spline_curve::BSplineCurve;
use crate::nurbs::b_spline_surface::BSplineSurface;
use crate::nurbs::crv_basic::Curve as NurbCurve;
use crate::nurbs::knots::MAX_DEGREE;
use crate::nurbs::math::{CPoint, NOW, NOZ};
use crate::nurbs::srf_basic::Surface as NurbSurface;
use crate::surface::implicit::{AnalyticSupport, SupportKind};
use crate::surface::surface::{Surface, SurfaceKind};
use crate::topology::FaceId;

const OP_TO: &str = "surface::procedural_blend::convert::to_record";
const OP_FROM: &str = "surface::procedural_blend::convert::from_record";
const OP_FACE: &str = "surface::procedural_blend::convert::rebuild_for_face";

/// The record of a procedural blend: its definition (supports by value, nested blends as
/// their own records) and its net, every number as the kernel stores it.
///
/// # Errors
///
/// `Unsupported` for a support that is not a NURBS surface (an offset, a curve-bounded or a
/// STEP surface: the definition's stations evaluate those by their own evaluators, which a
/// record does not carry); `SM_ERR` for a curve or surface without NURB data; those of
/// making an analytic support's NURBS form.
pub fn to_record(k: &Kernel, pb: &ProceduralBlend) -> Result<ProceduralBlendRecord> {
    to_record_at(k, pb, 1)
}

fn to_record_at(k: &Kernel, pb: &ProceduralBlend, depth: u32) -> Result<ProceduralBlendRecord> {
    if depth > BLEND_MAX_DEPTH {
        return Err(SolidError::invalid_input(
            OP_TO,
            "blends nest at most 3 definitions deep",
        ));
    }
    let d = pb.definition();
    let [s0, s1] = pb.supports();
    let supports = [support_record(k, s0, depth)?, support_record(k, s1, depth)?];
    let radius = match &d.law {
        Some(law) => RadiusRecord::Law(curve_record(law)?),
        None => RadiusRecord::Constant { offsets: d.offsets },
    };
    let [seed0, seed1] = &d.seeds;
    Ok(ProceduralBlendRecord {
        definition: BlendDefinitionRecord {
            supports,
            sides: d.sides,
            radius,
            order: d.order,
            thumbweight: d.thumbweight,
            leader: curve_record(&d.leader)?,
            seeds: [curve_record(seed0)?, curve_record(seed1)?],
            bound: d.bound,
            map: ParameterMapRecord {
                a: d.map_a,
                b: d.map_b,
                swapped: d.swapped,
            },
        },
        net: surface_record(pb.net())?,
    })
}

/// A support's record: a nested blend's own record, or the NURBS form its stations are
/// solved on ([`as_b_spline`]: an analytic support is evaluated as its net).
fn support_record(k: &Kernel, s: &Surface, depth: u32) -> Result<SupportRecord> {
    if let SurfaceKind::ProceduralBlend(inner) = s.kind() {
        return Ok(SupportRecord::Blend(Box::new(to_record_at(
            k,
            inner,
            depth + 1,
        )?)));
    }
    let Some(b) = as_b_spline(k, s)? else {
        return Err(SolidError::Unsupported {
            operation: OP_TO,
            reason: "a support that is not a NURBS surface has no record",
        });
    };
    let kind = match AnalyticSupport::of(s).map(|a| a.kind()) {
        None => SupportKindRecord::Free,
        Some(SupportKind::Plane) => SupportKindRecord::Plane,
        Some(SupportKind::Cylinder) => SupportKindRecord::Cylinder,
        Some(SupportKind::Cone) => SupportKindRecord::Cone,
        Some(SupportKind::Sphere) => SupportKindRecord::Sphere,
        Some(SupportKind::Torus) => SupportKindRecord::Torus,
    };
    Ok(SupportRecord::Nurbs {
        kind,
        form: surface_record(&b)?,
    })
}

/// The kernel's kind a record's label names.
const fn kind_of(k: SupportKindRecord) -> Option<SupportKind> {
    match k {
        SupportKindRecord::Free => None,
        SupportKindRecord::Plane => Some(SupportKind::Plane),
        SupportKindRecord::Cylinder => Some(SupportKind::Cylinder),
        SupportKindRecord::Cone => Some(SupportKind::Cone),
        SupportKindRecord::Sphere => Some(SupportKind::Sphere),
        SupportKindRecord::Torus => Some(SupportKind::Torus),
    }
}

fn to_u32(x: usize) -> Result<u32> {
    u32::try_from(x).map_err(|_| SolidError::invalid_input(OP_TO, "a count past u32"))
}

/// A NURBS curve's stored numbers.
fn curve_record(c: &BSplineCurve) -> Result<NurbsCurveRecord> {
    let nurb = c.nurb.as_ref().ok_or_else(|| SolidError::sm_err(OP_TO))?;
    let rational = nurb.is_rational();
    let points_3d = nurb.is_3d();
    let cpts = nurb.cpts();
    Ok(NurbsCurveRecord {
        dim: u32::try_from(c.dim).map_err(|_| SolidError::sm_err(OP_TO))?,
        points_3d,
        degree: to_u32(nurb.degree())?,
        knots: nurb.knots().to_vec(),
        points: cpts
            .iter()
            .map(|p| [p.x, p.y, if points_3d { p.z } else { 0.0 }])
            .collect(),
        weights: rational.then(|| cpts.iter().map(|p| p.w).collect()),
    })
}

/// A NURBS surface's stored numbers.
fn surface_record(s: &BSplineSurface) -> Result<NurbsSurfaceRecord> {
    let nurb = s.gw_nurb().ok_or_else(|| SolidError::sm_err(OP_TO))?;
    let rational = nurb.is_rational();
    let mut points = Vec::new();
    let mut weights = Vec::new();
    for row in nurb.rows() {
        for p in row {
            points.push([p.x, p.y, p.z]);
            weights.push(p.w);
        }
    }
    Ok(NurbsSurfaceRecord {
        u_degree: to_u32(nurb.p())?,
        v_degree: to_u32(nurb.q())?,
        u_count: to_u32(nurb.n().saturating_add(1))?,
        v_count: to_u32(nurb.m().saturating_add(1))?,
        u_knots: nurb.u_knots().to_vec(),
        v_knots: nurb.v_knots().to_vec(),
        points,
        weights: rational.then_some(weights),
    })
}

/// The blend a record holds, exactly: its supports, curves and net as their stored
/// numbers, the definition made by [`ProceduralBlend::create`] (which refuses an order,
/// side, radius, thumbweight or bound out of range) with the record's map. The net is
/// taken as it is, not checked against the definition: [`rebuild`] does that.
///
/// # Errors
///
/// `InvalidInput` for a record its structural check refuses
/// ([`ProceduralBlendRecord::validate`]), a degree past the kernel's, knots that decrease,
/// a weight that is not positive, or a singular map; those of
/// [`ProceduralBlend::create`].
pub fn from_record(k: &Kernel, rec: &ProceduralBlendRecord) -> Result<ProceduralBlend> {
    rec.validate()
        .map_err(|e| SolidError::invalid_input(OP_FROM, e.reason))?;
    from_record_at(k, rec, 1)
}

/// [`from_record`] and the net re-certified against the definition
/// ([`ProceduralBlend::recertify`]): the blend a record holds, refused unless its net is
/// proved within the record's bound of its definition.
///
/// # Errors
///
/// Those of [`from_record`] and [`ProceduralBlend::recertify`].
pub fn rebuild(k: &Kernel, rec: &ProceduralBlendRecord) -> Result<ProceduralBlend> {
    let pb = from_record(k, rec)?;
    pb.recertify(k)?;
    Ok(pb)
}

/// [`rebuild`] for the face `face` of a body: a refusal names the face
/// (`SM_ERR_INVALID_INPUT` with [`EntityRef::Face`]), as a body read back from its record
/// reports it.
///
/// # Errors
///
/// Those of [`rebuild`], naming `face`.
pub fn rebuild_for_face(
    k: &Kernel,
    face: FaceId,
    rec: &ProceduralBlendRecord,
) -> Result<ProceduralBlend> {
    rebuild(k, rec).map_err(|_| SolidError::Upstream {
        operation: OP_FACE,
        code: SM_ERR_INVALID_INPUT,
        entity: Some(EntityRef::Face(face)),
    })
}

fn from_record_at(k: &Kernel, rec: &ProceduralBlendRecord, depth: u32) -> Result<ProceduralBlend> {
    if depth > BLEND_MAX_DEPTH {
        return Err(SolidError::invalid_input(
            OP_FROM,
            "blends nest at most 3 definitions deep",
        ));
    }
    let d = &rec.definition;
    let mut supports = Vec::new();
    for (i, s) in d.supports.iter().enumerate() {
        supports.push(match s {
            SupportRecord::Nurbs { kind, form } => {
                let form = surface_of(form)?;
                match kind_of(*kind) {
                    None => Surface::from(form),
                    // the analytic surface derived from the form, the form checked on it
                    Some(kind) => AnalyticSupport::derive(kind, &form)
                        .and_then(|a| a.build(k, form))
                        .map_err(|e| match e {
                            SolidError::InvalidInput { .. } => SolidError::invalid_input(
                                OP_FROM,
                                if i == 0 {
                                    "support 1's NURBS form is not of its kind"
                                } else {
                                    "support 2's NURBS form is not of its kind"
                                },
                            ),
                            other => other,
                        })?,
                }
            }
            SupportRecord::Blend(b) => Surface::from(from_record_at(k, b, depth + 1)?),
        });
    }
    let (law, offsets) = match &d.radius {
        RadiusRecord::Constant { offsets } => (None, *offsets),
        RadiusRecord::Law(law) => (Some(curve_of(law)?), [0.0, 0.0]),
    };
    let leader = curve_of(&d.leader)?;
    let [seed0, seed1] = &d.seeds;
    let seeds = [curve_of(seed0)?, curve_of(seed1)?];
    let net = surface_of(&rec.net)?;
    let (Some(s0), Some(s1)) = (supports.first(), supports.get(1)) else {
        return Err(SolidError::sm_err(OP_FROM));
    };
    if !(d.map.b[0] != 0.0 && d.map.b[1] != 0.0) {
        return Err(SolidError::invalid_input(
            OP_FROM,
            "its parametrization is singular",
        ));
    }
    let input = BlendDefinitionInput {
        supports: [s0, s1],
        sides: d.sides,
        offsets,
        law: law.as_ref(),
        order: d.order,
        thumbweight: d.thumbweight,
        leader: &leader,
        seeds: [&seeds[0], &seeds[1]],
        bound: d.bound,
    };
    let mut pb = ProceduralBlend::create(k, &net, &input, |n, w| {
        GnSectionFsg::is_valid(n as usize, w)
    })?;
    pb.def.map_a = d.map.a;
    pb.def.map_b = d.map.b;
    pb.def.swapped = d.map.swapped;
    Ok(pb)
}

/// `degree` as a kernel degree: 1 to `NL_MAXDEG` (the NURBS evaluators size their stack
/// arrays by it, as the database readers check).
fn degree_of(degree: u32) -> Result<usize> {
    let p = degree as usize;
    if !(1..=MAX_DEGREE).contains(&p) {
        return Err(SolidError::invalid_input(
            OP_FROM,
            "a degree is 1 to the kernel's highest",
        ));
    }
    Ok(p)
}

fn non_decreasing(knots: &[f64]) -> Result<()> {
    if knots.windows(2).all(|w| matches!(w, [a, b] if a <= b)) {
        Ok(())
    } else {
        Err(SolidError::invalid_input(
            OP_FROM,
            "a knot vector decreases",
        ))
    }
}

/// The stored point `i` of a record (weighted when rational).
fn stored_point(
    points: &[[f64; 3]],
    weights: Option<&Vec<f64>>,
    i: usize,
    z: bool,
) -> Result<CPoint> {
    let [x, y, pz] = *points.get(i).ok_or_else(|| SolidError::sm_err(OP_FROM))?;
    let w = match weights {
        Some(ws) => {
            let w = *ws.get(i).ok_or_else(|| SolidError::sm_err(OP_FROM))?;
            if w.partial_cmp(&0.0) != Some(core::cmp::Ordering::Greater) {
                return Err(SolidError::invalid_input(
                    OP_FROM,
                    "a weight is not positive",
                ));
            }
            w
        }
        None => NOW,
    };
    Ok(CPoint::new(x, y, if z { pz } else { NOZ }, w))
}

/// The curve a record holds, its NURB as stored (no knot fix-up).
fn curve_of(r: &NurbsCurveRecord) -> Result<BSplineCurve> {
    let p = degree_of(r.degree)?;
    non_decreasing(&r.knots)?;
    let mut pw = Vec::new();
    for i in 0..r.points.len() {
        pw.push(stored_point(&r.points, r.weights.as_ref(), i, r.points_3d)?);
    }
    let nurb = NurbCurve::new(pw, p, r.knots.clone())?;
    let dim = u64::from(r.dim);
    Ok(BSplineCurve {
        nurb: Some(nurb),
        ..BSplineCurve::new(dim)
    })
}

/// The surface a record holds, its NURB as stored (no knot fix-up).
fn surface_of(r: &NurbsSurfaceRecord) -> Result<BSplineSurface> {
    let p = degree_of(r.u_degree)?;
    let q = degree_of(r.v_degree)?;
    non_decreasing(&r.u_knots)?;
    non_decreasing(&r.v_knots)?;
    let (nu, nv) = (r.u_count as usize, r.v_count as usize);
    let mut net = Vec::new();
    for i in 0..nu {
        let mut row = Vec::new();
        for j in 0..nv {
            let at = i
                .checked_mul(nv)
                .and_then(|x| x.checked_add(j))
                .ok_or_else(|| SolidError::sm_err(OP_FROM))?;
            row.push(stored_point(&r.points, r.weights.as_ref(), at, true)?);
        }
        net.push(row);
    }
    let nurb = NurbSurface::new(net, p, q, r.u_knots.clone(), r.v_knots.clone())?;
    Ok(BSplineSurface {
        nurb: Some(nurb).into(),
        ..BSplineSurface::new()
    })
}
