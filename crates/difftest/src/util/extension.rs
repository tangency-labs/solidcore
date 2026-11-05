//! The exact extension (`Surface::extended_sides`, `BSplineCurve::extended_by_distance`) in
//! the terms of the C++ extension the runners compare it with (#2730): the C++ side still
//! runs `CreateExtendedSurface` / `N_SrfExtendByDist` / `N_CrvExtendByDist`, whose least-squares
//! fits solidcore no longer has, so these cases move (their records name them).

use solidcore::geom::geom_utility::p2::ContinuityType;
use solidcore::kernel::Kernel;
use solidcore::nurbs::b_spline_curve::BSplineCurve;
use solidcore::nurbs::b_spline_surface::{BSplineSurface, SurfParamType};
use solidcore::nurbs::crv_adv::p1::{Continuity as NlibContinuity, CurveEnd};
use solidcore::nurbs::crv_basic::Curve as NurbsCurve;
use solidcore::nurbs::crv_geom::ToleranceKind;
use solidcore::nurbs::srf_basic::{Direction, Surface as Net};
use solidcore::surface::extend::Continuity;
use solidcore::surface::side::BoxSide;
use solidcore::surface::surface::{Surface, SurfaceKind};
use solidcore::{Result, SolidError};

const OP: &str = "difftest::util::extension";

/// The continuity a caller asking for `c` gets: Natural for `SM_CT_CINFINITY`, Curvature for
/// `SM_CT_G1_G2`, Ruled for `SM_CT_G1R`, Tangent for `SM_CT_G1` and every other type (as
/// upstream's `_ => NL_G1`).
#[must_use]
pub fn continuity(c: ContinuityType) -> Continuity {
    match c {
        ContinuityType::CInfinity => Continuity::Natural,
        ContinuityType::G1G2 => Continuity::Curvature,
        ContinuityType::G1R => Continuity::Ruled,
        _ => Continuity::Tangent,
    }
}

/// The NLib continuity flag (`NL_G1 = 2`, `NL_G2 = 3`, `NL_CMAX = 4`, `NL_G1R = 6`) as the
/// exact extension's; `None` for another flag.
#[must_use]
pub fn continuity_of_flag(flag: i32) -> Option<Continuity> {
    match flag {
        2 => Some(Continuity::Tangent),
        6 => Some(Continuity::Ruled),
        3 => Some(Continuity::Curvature),
        4 => Some(Continuity::Natural),
        _ => None,
    }
}

/// The sides of a parameter box `p` names.
#[must_use]
pub fn sides(p: SurfParamType) -> &'static [BoxSide] {
    match p {
        SurfParamType::U => &[BoxSide::UMin, BoxSide::UMax],
        SurfParamType::V => &[BoxSide::VMin, BoxSide::VMax],
        SurfParamType::UMin => &[BoxSide::UMin],
        SurfParamType::UMax => &[BoxSide::UMax],
        SurfParamType::VMin => &[BoxSide::VMin],
        SurfParamType::VMax => &[BoxSide::VMax],
        SurfParamType::Both => &BoxSide::ALL,
        _ => &[],
    }
}

/// `s` extended `dist` past the sides `dir` names (every side for `None`).
///
/// # Errors
///
/// Those of `Surface::extended_sides`.
pub fn surface(
    k: &Kernel,
    s: &Surface,
    dir: Option<SurfParamType>,
    dist: f64,
    cont: ContinuityType,
) -> Result<Surface> {
    let asked = dir.map_or(&BoxSide::ALL[..], sides);
    Ok(s.extended_sides(k, asked, continuity(cont), dist)?.surface)
}

/// [`surface`] of a B-spline surface, as a B-spline surface (a plane widened as its NURBS
/// form).
///
/// # Errors
///
/// Those of [`surface`] and of making the NURB.
pub fn b_spline(
    k: &Kernel,
    s: &BSplineSurface,
    dir: Option<SurfParamType>,
    dist: f64,
    cont: ContinuityType,
) -> Result<BSplineSurface> {
    let e = surface(k, &Surface::from(s.copy()), dir, dist, cont)?;
    as_b_spline(k, &e)
}

/// The NURBS form of `s` as a B-spline surface.
///
/// # Errors
///
/// Those of making the NURB; `SM_ERR` for a surface without one.
pub fn as_b_spline(k: &Kernel, s: &Surface) -> Result<BSplineSurface> {
    if let SurfaceKind::BSpline(b) = s.kind() {
        return Ok(b.copy());
    }
    solidcore::blend::fillet_section::station::as_b_spline(k, s)?
        .map(std::borrow::Cow::into_owned)
        .ok_or_else(|| SolidError::sm_err(OP))
}

/// A NURBS net extended `dist` past the end `end` of direction `dir` (NLib's direction flag:
/// every value but `U` is v), as `N_SrfExtendByDist` is asked.
///
/// # Errors
///
/// Those of `Surface::extended_by_distance`; `InvalidInput` for a continuity flag it has no
/// counterpart of.
pub fn net(
    k: &Kernel,
    n: &Net,
    dist: f64,
    dir: Direction,
    end: CurveEnd,
    flag: i32,
) -> Result<Net> {
    let how = continuity_of_flag(flag)
        .ok_or_else(|| SolidError::invalid_input(OP, "the continuity flag"))?;
    let side = match (dir == Direction::U, end) {
        (true, CurveEnd::Start) => BoxSide::UMin,
        (true, CurveEnd::End) => BoxSide::UMax,
        (false, CurveEnd::Start) => BoxSide::VMin,
        (false, CurveEnd::End) => BoxSide::VMax,
    };
    let s = Surface::from(BSplineSurface::from_gw_nurb(n));
    let dom = s.natural_uv_domain_with_nurb(k)?;
    let e = s.extended_by_distance(k, &dom, side, how, dist)?;
    as_b_spline(k, &e.surface)?
        .gw_nurb()
        .cloned()
        .ok_or_else(|| SolidError::sm_err(OP))
}

/// An NLib curve extended `dist` past `end` (absolute, or `relative` to the curve's arc
/// length), as `N_CrvExtendByDist` is asked with `cont`.
///
/// # Errors
///
/// Those of `BSplineCurve::extended_by_distance` and of the arc length.
pub fn nlib_curve(
    c: &NurbsCurve,
    dist: f64,
    relative: bool,
    end: CurveEnd,
    cont: NlibContinuity,
) -> Result<NurbsCurve> {
    let how = match cont {
        NlibContinuity::G1 => Continuity::Tangent,
        NlibContinuity::G2 => Continuity::Curvature,
        NlibContinuity::Cmax => Continuity::Natural,
    };
    let dist = if relative {
        let (a, b) = c.param_bounds();
        dist * c.arc_length(a, b, 1e-12, ToleranceKind::Relative)?
    } else {
        dist
    };
    let e = BSplineCurve::from_gw_nurb(3, c).extended_by_distance(end, how, dist)?;
    e.curve
        .gw_nurb()
        .cloned()
        .ok_or_else(|| SolidError::sm_err(OP))
}
