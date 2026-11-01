//! The G^n section generator: the blend surface of the order-n section between the rails,
//! on stations at the knots of the rails and the centre line (and where a foot crosses a
//! face's knot line), certified span by span, and the section's order of contact with a
//! face at a rail.
//!
//! Fork: `SmGnSectionFSG` (`SM_GN_MAX_ORDER`, `SM_GN_DEFAULT_THUMBWEIGHT`, the constructor,
//! `IsValid`, `CreateSurface`, `CreateSurfaceFromCurves`, `Build`, `ContactOrderAt`,
//! `GetOrder`, `GetThumbweight`, `GetCertifiedBound`, `GetCertifiedSpans`) in
//! `source/SMLib/inc/SmFilletSolver.h` and `source/SMLib/src/SmFilletSection.cpp`
//! (tangency-labs/smlib@11ea2d9; the equal area thumbweight resolved at every station
//! 7102bc5, any thumbweight below 1 b9d6a72); see [`super`] for the table.

use super::certificate::{gn_certify_span, GnVarying, NurbsCurve};
use super::form::{gn_form, gn_form_normalized, hermite_knots, hermite_points};
use super::series::Trunc;
use super::station::{gn_solve_station, smax, surface_partials, vec3, GnStation, Nurbs};
use super::tay::{
    acos_clamped, gn_area_sum, gn_half_angle, gn_invert, gn_target, gn_wedge_cos, GnAreaCells,
    Ring, Tay,
};
use crate::error::{Result, SolidError};
use crate::kernel::Kernel;
use crate::math::{Point2, Point3, Vec3, BIG_DOUBLE, PI};
use crate::nurbs::b_spline_curve::{BSplineCurve, BSplineCurveForm, CanonicalData, KnotType};
use crate::nurbs::b_spline_surface::p1::CanonicalSurfaceData;
use crate::nurbs::b_spline_surface::{BSplineSurface, BSplineSurfaceForm};
use crate::nurbs::knots::Side;
use crate::surface::surface::{Surface, SurfaceKind};

const OP: &str = "blend::fillet_section::GnSectionFsg";

/// `SM_GN_MAX_ORDER`: the highest order of a G^n section.
pub const GN_MAX_ORDER: usize = 4;
/// `SM_GN_THUMBWEIGHT_EQUAL_AREA` (the fork's 6148773 = d01d7b0's 1947a45; 2, outside the
/// thumbweights below 1, since b9d6a72): a value outside the thumbweights' range that a
/// generator resolves at every station, from that station's wedge, to the w at which its
/// section, on the same rails, takes off the area the ball's circle does
/// ([`GnSectionFsg::equal_area_thumbweight`]; 7102bc5): w varies along the blend.
pub const GN_THUMBWEIGHT_EQUAL_AREA: f64 = 2.0;
/// `SM_GN_DEFAULT_THUMBWEIGHT`: the equal area one, [`GN_THUMBWEIGHT_EQUAL_AREA`] (0, the
/// section with no thumbweight, before 6148773).
pub const GN_DEFAULT_THUMBWEIGHT: f64 = GN_THUMBWEIGHT_EQUAL_AREA;

/// `SmGnSectionFSG`: a blend's cross-section of contact order n with its faces along the
/// rails, `n = 2..=GN_MAX_ORDER`, with thumbweight w. It needs a ball of one radius
/// (constant offsets).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GnSectionFsg {
    /// `m_lOrder`: n.
    pub order: usize,
    /// `m_dThumbweight`: w, `w < 1`, or [`GN_THUMBWEIGHT_EQUAL_AREA`] (resolved per
    /// station, not here).
    pub thumbweight: f64,
    /// `m_dCertifiedBound`: the last surface built: its certified distance from the blend
    /// (-1: none yet).
    pub certified_bound: f64,
    /// `m_lCertifiedSpans`: and its number of spans along the rails.
    pub certified_spans: usize,
}

/// The faces' surfaces under a generator's rails and their uv curves.
#[derive(Clone, Copy, Debug)]
pub struct GnFaces<'a> {
    /// The faces' surfaces (`pSurf1`, `pSurf2`).
    pub surfaces: [&'a Surface; 2],
    /// The rails' uv curves on them (`pUV1`, `pUV2`).
    pub uv: [&'a BSplineCurve; 2],
}

/// The surface under an offset surface (`GetBaseSurface`), else the surface itself.
fn base_of(s: &Surface) -> &Surface {
    match s.kind() {
        SurfaceKind::OffsetSurface(o) => o.base_surface(),
        _ => s,
    }
}

/// `GetNaturalInterval` as `(min, max)`.
fn interval(c: &BSplineCurve) -> (f64, f64) {
    let i = c.natural_interval();
    (i.min(), i.max())
}

/// `EvaluateNormal(uv, TRUE, TRUE, sN)` then `Unitize` (both `SER`).
fn unit_normal(k: &Kernel, s: &Surface, uv: Point2) -> Result<Vec3> {
    s.evaluate_normal(k, uv, Side::Left, Side::Left)?.unitize()
}

/// A station's foot on face f at sigma^i eps^e.
fn foot(st: &GnStation, f: usize, i: usize, e: usize) -> Vec3 {
    let [a, b] = &st.f;
    let fs = if f == 0 { a } else { b };
    let [x, y, z] = &fs.a;
    Vec3::new(x.at(i, e), y.at(i, e), z.at(i, e))
}

/// A station's foot parameter q (u1, v1, u2, v2) at sigma^i eps^e.
fn param(st: &GnStation, q: usize, i: usize, e: usize) -> f64 {
    st.x.get(q).map_or(f64::NAN, |x| x.at(i, e))
}

/// `smos_Min(a, b)`: `(a < b) ? a : b`.
#[inline]
fn smin(a: f64, b: f64) -> f64 {
    if a < b {
        a
    } else {
        b
    }
}

/// `sm_GnEqualAreaW`: the thumbweight at which the order-n section takes off what the
/// ball's circle does, for half the angle u between the ball's radii
/// ([`GnSectionFsg::equal_area_thumbweight`]): the root of `I(w) = target(u)`.
fn gn_equal_area_w(order: usize, u: f64) -> Result<f64> {
    const OP_EA: &str = "blend::fillet_section::GnSectionFsg::equal_area_thumbweight";
    if !((2..=GN_MAX_ORDER).contains(&order) && u > 0.0 && u < 0.5 * PI) {
        return Err(SolidError::invalid_input(OP_EA, "order or wedge angle"));
    }
    let target = gn_target(&Tay::<f64>::lift(u)).at(0);
    let mut lo = -1.0_f64;
    let mut hi = 1.0 - 1.0e-9;
    let mut held = target >= gn_area_sum::<f64>(order, hi);
    while held && gn_area_sum::<f64>(order, lo) < target {
        hi = lo;
        lo = 2.0 * lo - 1.0;
        held = lo >= -1.0e3;
    }
    if !held {
        return Err(SolidError::invalid_input(
            OP_EA,
            "no thumbweight gives this wedge's section the ball's circle's area",
        ));
    }
    let mut it = 0;
    while it < 200 && hi - lo > 1.0e-15 * (1.0 - lo) {
        let mid = 0.5 * (lo + hi);
        if gn_area_sum::<f64>(order, mid) > target {
            lo = mid;
        } else {
            hi = mid;
        }
        it += 1;
    }
    Ok(0.5 * (lo + hi))
}

/// `sm_GnStationThumbweight`: a station's equal area thumbweight and its s-derivative: x
/// from the station's feet and V to first order in s, w the root at its u,
/// `dw/ds = target' / I'(w)`.
fn gn_station_thumbweight(n: usize, rr: f64, st: &mut GnStation) -> Result<()> {
    let u = gn_half_angle(&gn_wedge_cos::<f64>(&st.f, &st.v, rr, 1));
    let w = gn_equal_area_w(n, u.at(0))?;
    let d = gn_area_sum(n, Tay::of_terms(w, 1.0));
    let falls = d.at(1) < 0.0;
    if !falls {
        // "the section's area does not fall with its thumbweight"
        return Err(SolidError::sm_err(OP));
    }
    let sw = gn_invert(&gn_target(&u), &d, w);
    st.w = w;
    st.wds = sw.at(1);
    Ok(())
}

impl GnSectionFsg {
    /// `SmGnSectionFSG(lOrder, dThumbweight)`.
    #[must_use]
    pub const fn new(order: usize, thumbweight: f64) -> Self {
        Self {
            order,
            thumbweight,
            certified_bound: -1.0,
            certified_spans: 0,
        }
    }

    /// `IsValid`: whether an order and a thumbweight make a G^n section:
    /// `2 <= n <= GN_MAX_ORDER`, and `w < 1` (the section's weights `1 - w + w t` and
    /// `1 - w t` stay positive on [0, 1] and its reparametrizations monotone; b9d6a72) or
    /// the equal area default.
    #[must_use]
    pub fn is_valid(order: usize, thumbweight: f64) -> bool {
        (2..=GN_MAX_ORDER).contains(&order)
            && (thumbweight == GN_THUMBWEIGHT_EQUAL_AREA
                || (thumbweight.is_finite() && thumbweight < 1.0))
    }

    /// `EqualAreaThumbweight`: the w at which the order-n section between two planes
    /// meeting at the angle `alpha` (radians, the angle of the wedge the ball sits in,
    /// `0 < alpha < pi`) takes off, on the same rails, the area the ball's circle does:
    /// `r^2 (cot(alpha/2) - (pi - alpha)/2)`.
    ///
    /// With `e1`, `e2` the rails' points from the edge (`|e_i| = d = r cot(alpha/2)`), the
    /// section takes off `(1/2) d^2 sin(alpha) I(w)`, `I` the integral of `a b' - a' b`
    /// (`sm_GnAreaSum`), so `I(w) = (cot(alpha/2) - (pi - alpha)/2) / (cot(alpha/2)
    /// cos^2(alpha/2))`: with `u = (pi - alpha)/2`, half the angle between the ball's radii
    /// to its contacts, `(sin u - u cos u) / sin^3 u`. `I` falls as w rises (to 0 as w
    /// nears 1) and nears the chord's 1 as w falls, so the root is bisected on a bracket
    /// whose lower end steps down from -1 until it holds it (b9d6a72): a sharp wedge's w is
    /// well below that.
    ///
    /// # Errors
    ///
    /// `SM_ERR_INVALID_INPUT` for an order past 2..=4, an angle outside (0, pi), or a
    /// target the bracket does not reach by w = -1000 (a wedge so open or so shut that the
    /// circle takes off next to nothing, or all but the chord's area).
    pub fn equal_area_thumbweight(order: usize, alpha: f64) -> Result<f64> {
        if !(alpha > 0.0 && alpha < PI) {
            return Err(SolidError::invalid_input(
                "blend::fillet_section::GnSectionFsg::equal_area_thumbweight",
                "wedge angle",
            ));
        }
        gn_equal_area_w(order, 0.5 * (PI - alpha))
    }

    /// `ThumbweightAt`: the thumbweight at a station whose ball's radii to its contacts
    /// meet at the angle whose cosine is `cos_radii`: the one given, or the equal area one
    /// of the wedge `alpha = pi - (that angle)`. Exact for two planes; elsewhere the
    /// section of the station's own wedge (7102bc5).
    ///
    /// # Errors
    ///
    /// Those of [`GnSectionFsg::equal_area_thumbweight`].
    pub fn thumbweight_at(&self, cos_radii: f64) -> Result<f64> {
        if self.thumbweight != GN_THUMBWEIGHT_EQUAL_AREA {
            return Ok(self.thumbweight);
        }
        gn_equal_area_w(self.order, 0.5 * acos_clamped(cos_radii))
    }

    /// `CreateSurface`: the blend surface between the solver's rails. The solver's faces
    /// are its offset surfaces; the section is built on the faces under them (an offset
    /// surface's base, else the surface itself). The rails and their uv curves are rebuilt
    /// on the section's stations from the same exact feet, so that the surface's
    /// boundaries are the rails.
    ///
    /// # Errors
    ///
    /// Those of `Build`.
    pub fn create_surface(
        &mut self,
        k: &Kernel,
        approx_tol_3d: f64,
        center_line: &BSplineCurve,
        rails: [&mut BSplineCurve; 2],
        surfaces: [&Surface; 2],
        uvs: [&mut BSplineCurve; 2],
    ) -> Result<BSplineSurface> {
        let [s1, s2] = surfaces;
        let faces = [base_of(s1), base_of(s2)];
        self.build(k, center_line, rails, approx_tol_3d, faces, uvs, true)
    }

    /// `CreateSurfaceFromCurves`: the blend surface between two rails on the faces'
    /// surfaces, about the centre line. The rails and uv curves are not changed.
    ///
    /// # Errors
    ///
    /// Those of `Build`.
    pub fn create_surface_from_curves(
        &mut self,
        k: &Kernel,
        center_line: &BSplineCurve,
        rails: [&mut BSplineCurve; 2],
        approx_tol_3d: f64,
        faces: GnFaces<'_>,
    ) -> Result<BSplineSurface> {
        let [u1, u2] = faces.uv;
        let mut uv1 = u1.clone();
        let mut uv2 = u2.clone();
        self.build(
            k,
            center_line,
            rails,
            approx_tol_3d,
            faces.surfaces,
            [&mut uv1, &mut uv2],
            false,
        )
    }

    /// `Build`: the G^n blend surface between the rails, its stations the knots of the
    /// rails and the centre line. `replace_rails`: the rails and uv curves become the
    /// jets' own order-0 curves (the surface's boundaries). The surface's u is the rails'
    /// parameter, v the section's t, from the first rail (v = 0) to the second (v = 1).
    #[allow(clippy::too_many_arguments)]
    fn build(
        &mut self,
        k: &Kernel,
        center: &BSplineCurve,
        rails: [&mut BSplineCurve; 2],
        tol: f64,
        faces: [&Surface; 2],
        uvs: [&mut BSplineCurve; 2],
        replace_rails: bool,
    ) -> Result<BSplineSurface> {
        let n = self.order;
        if !Self::is_valid(n, self.thumbweight) {
            return Err(SolidError::invalid_input(
                OP,
                "the order must be 2 to 4 and the thumbweight below 1",
            ));
        }
        let [rail1, rail2] = rails;
        let (c_min, c_max) = interval(center);
        let scale = smax(1.0, c_max - c_min);
        // the equal area thumbweight is resolved per station, from its own wedge
        let equal_area = self.thumbweight == GN_THUMBWEIGHT_EQUAL_AREA;

        // the stations: the unique knots of the rails and the centre line
        let mut ss: Vec<f64> = Vec::new();
        for c in [center, &*rail1, &*rail2] {
            let (lo, hi) = interval(c);
            if (lo - c_min).abs() > 1.0e-9 * scale || (hi - c_max).abs() > 1.0e-9 * scale {
                // "the rails and the centre line run over different intervals"
                return Err(SolidError::sm_err(OP));
            }
            ss.extend(c.knots(None)?.0);
        }
        ss.sort_by(f64::total_cmp);
        {
            let mut unique: Vec<f64> = Vec::new();
            for &s in &ss {
                if unique.last().is_none_or(|&b| s - b > 1.0e-9 * scale) {
                    unique.push(s);
                }
            }
            if let Some(f) = unique.first_mut() {
                *f = c_min;
            }
            if let Some(l) = unique.last_mut() {
                *l = c_max;
            }
            ss = unique;
        }
        let (Some(&s_first), true) = (ss.first(), ss.len() >= 2) else {
            // "the rails have no length"
            return Err(SolidError::sm_err(OP));
        };

        // the offsets and the sides, at the first station: the centre line there is the
        // ball's centre
        let mut r = [0.0_f64; 2];
        let mut side = [0.0_f64; 2];
        {
            let m = center.evaluate_point(s_first)?;
            let rl: [&BSplineCurve; 2] = [&*rail1, &*rail2];
            let rows = r
                .iter_mut()
                .zip(side.iter_mut())
                .zip(rl)
                .zip(uvs.iter())
                .zip(faces);
            for ((((rf, sf), rail), uv), face) in rows {
                let rho = rail.evaluate_point(s_first)?;
                let w = uv.evaluate_point(s_first)?;
                let nrm = unit_normal(k, face, Point2::new(w.x, w.y))?;
                let off = m - rho;
                *rf = off.length();
                *sf = if off.dot(nrm) >= 0.0 { 1.0 } else { -1.0 };
                let has_radius = *rf > tol;
                if !has_radius {
                    // "the ball has no radius at the rails' start"
                    return Err(SolidError::sm_err(OP));
                }
            }
        }

        // every station's feet and their series: values and s-derivatives (sigma order
        // 1), jets to order n
        let tr = Trunc::new(1, n)?;
        let mut stations: Vec<GnStation> = Vec::new();
        for &s in &ss {
            let w0 = uvs[0].evaluate_point(s)?;
            let w1 = uvs[1].evaluate_point(s)?;
            let seed = [w0.x, w0.y, w1.x, w1.y];
            // a constant ball: the station's centre on the centre line at the same offsets
            let m = center.evaluate_point(s)?;
            let solve_tol = 1.0e-9 * smax(1.0, m.length());
            let st = gn_solve_station(k, faces, center, s, seed, r, side, tr, solve_tol)?;
            for (f, &rf) in r.iter().enumerate() {
                let rho = foot(&st, f, 0, 0);
                if (m.distance(rho) - rf).abs() > smax(tol, 1.0e-6 * rf) {
                    return Err(SolidError::invalid_input(
                        OP,
                        "a section of order 2 and up is built for a ball of one radius only",
                    ));
                }
            }
            stations.push(st);
        }

        // the certificate
        let (stations, bound, spans) =
            self.certify(k, faces, center, stations, r, side, tr, tol, equal_area)?;
        self.certified_bound = bound;
        self.certified_spans = spans;
        let ss: Vec<f64> = stations.iter().map(|st| st.s).collect();

        let ld = 3 * n + 2; // the section's control points across

        // the jets' Hermite curves along s, all on the stations' knots
        let (u_knots, u_mult) = hermite_knots(&ss);
        let mut jets: [Vec<Vec<Vec3>>; 2] = [Vec::new(), Vec::new()];
        for (f, jf) in jets.iter_mut().enumerate() {
            for kk in 0..=n {
                let v: Vec<Vec3> = stations.iter().map(|st| foot(st, f, 0, kk)).collect();
                let d: Vec<Vec3> = stations.iter().map(|st| foot(st, f, 1, kk)).collect();
                jf.push(hermite_points(&ss, &v, &d));
            }
        }
        let jet = |f: usize, kk: usize, i: usize| -> Vec3 {
            jets.get(f)
                .and_then(|j| j.get(kk))
                .and_then(|j| j.get(i))
                .copied()
                .unwrap_or(Vec3::new(f64::NAN, f64::NAN, f64::NAN))
        };
        let lnu = jets[0].first().map_or(0, Vec::len);

        // the homogeneous control points W P and weights W, [i along s][j across]
        let mut homog = vec![Vec3::new(0.0, 0.0, 0.0); lnu * ld];
        let mut wv = vec![0.0_f64; lnu * ld];
        if equal_area {
            // each station's own form (over W_0): every row W_j P_j and weight W_j a
            // Hermite curve of its values and s-derivatives there (the form's in w(s) and
            // the jets')
            let mut hv: Vec<Vec<Vec3>> = vec![Vec::new(); ld];
            let mut hd: Vec<Vec<Vec3>> = vec![Vec::new(); ld];
            let mut wvv: Vec<Vec<Vec3>> = vec![Vec::new(); ld];
            let mut wdd: Vec<Vec<Vec3>> = vec![Vec::new(); ld];
            for st in &stations {
                let sf = gn_form_normalized(n, Tay::of_terms(st.w, st.wds));
                for j in 0..ld {
                    let mut sv = Vec3::new(0.0, 0.0, 0.0);
                    let mut sd = Vec3::new(0.0, 0.0, 0.0);
                    for (f, gf) in sf.g.iter().enumerate() {
                        for (kk, gk) in gf.iter().enumerate().take(n + 1) {
                            let c = foot(st, f, 0, kk);
                            let c1 = foot(st, f, 1, kk);
                            let g = gk.get(j).copied().unwrap_or_else(|| Tay::lift(f64::NAN));
                            sv += c * g.at(0);
                            sd = sd + c * g.at(1) + c1 * g.at(0);
                        }
                    }
                    let wj = sf.w.get(j).copied().unwrap_or_else(|| Tay::lift(f64::NAN));
                    if let (Some(a), Some(b), Some(c), Some(d)) =
                        (hv.get_mut(j), hd.get_mut(j), wvv.get_mut(j), wdd.get_mut(j))
                    {
                        a.push(sv);
                        b.push(sd);
                        c.push(Vec3::new(wj.at(0), 0.0, 0.0));
                        d.push(Vec3::new(wj.at(1), 0.0, 0.0));
                    }
                }
            }
            for j in 0..ld {
                let (Some(a), Some(b), Some(c), Some(d)) =
                    (hv.get(j), hd.get(j), wvv.get(j), wdd.get(j))
                else {
                    return Err(SolidError::sm_err(OP));
                };
                let hp = hermite_points(&ss, a, b);
                let wp = hermite_points(&ss, c, d);
                for i in 0..lnu {
                    if let (Some(h), Some(w), Some(&p), Some(q)) = (
                        homog.get_mut(i * ld + j),
                        wv.get_mut(i * ld + j),
                        hp.get(i),
                        wp.get(i),
                    ) {
                        *h = p;
                        *w = q.x;
                    }
                }
            }
        } else {
            // one section's form: its rows are linear in the jets, so the jets' curves
            // carry them
            let form = gn_form(n, self.thumbweight);
            for i in 0..lnu {
                for j in 0..ld {
                    let mut p = Vec3::new(0.0, 0.0, 0.0);
                    for (f, gf) in form.g.iter().enumerate() {
                        for (kk, g) in gf.iter().enumerate().take(n + 1) {
                            p += jet(f, kk, i) * g.get(j).copied().unwrap_or(f64::NAN);
                        }
                    }
                    if let (Some(h), Some(w)) = (homog.get_mut(i * ld + j), wv.get_mut(i * ld + j))
                    {
                        *h = p;
                        *w = form.w.get(j).copied().unwrap_or(f64::NAN);
                    }
                }
            }
        }
        let mut w_min = BIG_DOUBLE;
        let mut w_max = 0.0_f64;
        for &w in &wv {
            w_min = if w_min < w { w_min } else { w }; // smos_Min
            w_max = smax(w_max, w);
        }
        let positive = w_min > 0.0;
        if !positive {
            // "a weight of the section's rows is not positive"
            return Err(SolidError::sm_err(OP));
        }
        let w_scale = 1.0 / crate::geom::vector3d::sqrt(w_min * w_max);
        let mut cps: Vec<Point3> = Vec::new();
        let mut wts: Vec<f64> = Vec::new();
        for (&h, &w) in homog.iter().zip(&wv) {
            cps.push(h * (1.0 / w));
            wts.push(w * w_scale);
        }
        let surface = BSplineSurface::create_canonical(&CanonicalSurfaceData {
            u_degree: 3,
            v_degree: ld - 1,
            control_points: &cps,
            form: BSplineSurfaceForm::Unspecified,
            u_multiplicities: &u_mult,
            v_multiplicities: &[ld, ld],
            u_knots: &u_knots,
            v_knots: &[0.0, 1.0],
            knot_type: KnotType::Unspecified,
            weights: Some(&wts),
            uv_domain: None,
        })?;

        if replace_rails {
            for (f, (rail, uv)) in [rail1, rail2].into_iter().zip(uvs).enumerate() {
                let rail_cps: Vec<Point3> = (0..lnu).map(|i| jet(f, 0, i)).collect();
                rail.set_canonical(&CanonicalData {
                    dim: 3,
                    degree: 3,
                    control_points: &rail_cps,
                    form: BSplineCurveForm::Unspecified,
                    knot_multiplicities: &u_mult,
                    knots: &u_knots,
                    knot_type: KnotType::Unspecified,
                    weights: None,
                    trim: None,
                })?;
                let v: Vec<Vec3> = stations
                    .iter()
                    .map(|st| Vec3::new(param(st, 2 * f, 0, 0), param(st, 2 * f + 1, 0, 0), 0.0))
                    .collect();
                let d: Vec<Vec3> = stations
                    .iter()
                    .map(|st| Vec3::new(param(st, 2 * f, 1, 0), param(st, 2 * f + 1, 1, 0), 0.0))
                    .collect();
                let uv_cps = hermite_points(&ss, &v, &d);
                let dim = uv.dim;
                uv.set_canonical(&CanonicalData {
                    dim,
                    degree: 3,
                    control_points: &uv_cps,
                    form: BSplineCurveForm::Unspecified,
                    knot_multiplicities: &u_mult,
                    knots: &u_knots,
                    knot_type: KnotType::Unspecified,
                    weights: None,
                    trim: None,
                })?;
            }
        }
        Ok(surface)
    }

    /// The certificate block of `Build`: a station wherever a foot crosses a knot line of
    /// its face, then each span's bound (`sm_GnCertifySpan`) at most a tenth of the
    /// tolerance, a span over it halved (a station solved at its middle) until none is.
    /// Returns the stations, the largest bound and the number of spans.
    #[allow(clippy::too_many_arguments)]
    fn certify(
        &self,
        k: &Kernel,
        faces: [&Surface; 2],
        center: &BSplineCurve,
        mut stations: Vec<GnStation>,
        r: [f64; 2],
        side: [f64; 2],
        tr: Trunc,
        tol: f64,
        equal_area: bool,
    ) -> Result<(Vec<GnStation>, f64, usize)> {
        let mut nurbs = Vec::new();
        for f in faces {
            let Some(bs) = super::station::as_b_spline(k, f)? else {
                // "a face that is not a NURBS surface: its blend of order 2 and up cannot
                // be certified"
                return Err(SolidError::sm_err(OP));
            };
            nurbs.push(Nurbs::of(&bs)?);
        }
        let (Some(n0), Some(n1)) = (nurbs.first(), nurbs.get(1)) else {
            return Err(SolidError::sm_err(OP));
        };
        let face_n = [n0, n1];
        let leader = NurbsCurve::of(center)?;
        let target = 0.1 * tol;
        let solve_tol = |s: f64| -> Result<f64> {
            let m = center.evaluate_point(s)?;
            Ok(1.0e-9 * smax(1.0, m.length()))
        };
        // a station wherever a foot crosses a knot line of its face: the face is only as
        // smooth there as the knot allows, and a span's enclosures are taken on one side
        for (f, nf) in face_n.iter().enumerate() {
            for c in 0..2 {
                let (kv, deg) = if c == 0 { (&nf.u, nf.p) } else { (&nf.v, nf.q) };
                let mut inner: Vec<f64> = Vec::new();
                let mut m = deg + 1;
                while m + deg + 1 < kv.len() {
                    let km = kv.get(m).copied().unwrap_or(f64::NAN);
                    if inner.last().is_none_or(|&b| km > b) {
                        inner.push(km);
                    }
                    m += 1;
                }
                let mut next: Vec<GnStation> = Vec::new();
                for pair in stations.windows(2) {
                    let [a, b] = pair else { continue };
                    next.push(*a);
                    for &dk in &inner {
                        let q = 2 * f + c;
                        let xa = param(a, q, 0, 0);
                        let xb = param(b, q, 0, 0);
                        let margin = 1.0e-9 * (1.0 + dk.abs());
                        let crosses = (xa < dk - margin && xb > dk + margin)
                            || (xa > dk + margin && xb < dk - margin);
                        if !crosses {
                            continue;
                        }
                        // the parameter where that foot is on the knot line: the secant,
                        // kept in its bracket
                        let (mut s0, mut s1, mut x0, mut x1) = (a.s, b.s, xa - dk, xb - dk);
                        let mut at: Option<GnStation> = None;
                        for _ in 0..60 {
                            let mut ds = s0 - x0 * (s1 - s0) / (x1 - x0);
                            let inside = ds > s0 && ds < s1;
                            if !inside {
                                ds = 0.5 * (s0 + s1);
                            }
                            let t = (ds - a.s) / (b.s - a.s);
                            let mut seed = [0.0_f64; 4];
                            for (q2, sq) in seed.iter_mut().enumerate() {
                                *sq = (1.0 - t) * param(a, q2, 0, 0) + t * param(b, q2, 0, 0);
                            }
                            let st = gn_solve_station(
                                k,
                                faces,
                                center,
                                ds,
                                seed,
                                r,
                                side,
                                tr,
                                solve_tol(ds)?,
                            )?;
                            let x = param(&st, q, 0, 0) - dk;
                            at = Some(st);
                            if x.abs() <= 1.0e-12 * (1.0 + dk.abs()) {
                                break;
                            }
                            if (x < 0.0) == (x0 < 0.0) {
                                s0 = ds;
                                x0 = x;
                            } else {
                                s1 = ds;
                                x1 = x;
                            }
                        }
                        if let Some(st) = at {
                            next.push(st);
                        }
                    }
                }
                if let Some(last) = stations.last() {
                    next.push(*last);
                }
                next.sort_by(|a, b| a.s.total_cmp(&b.s));
                stations = next;
            }
        }
        // the equal area thumbweight at every station, and I's enclosures on a grid over
        // their range (a little wider: the stations put in later lie between these)
        let n = self.order;
        let rr = r[0] * r[1];
        let mut vary = GnVarying {
            n,
            rr,
            cells: GnAreaCells::default(),
        };
        if equal_area {
            let mut w_lo = f64::INFINITY;
            let mut w_hi = f64::NEG_INFINITY;
            for st in &mut stations {
                gn_station_thumbweight(n, rr, st)?;
                w_lo = smin(w_lo, st.w);
                w_hi = smax(w_hi, st.w);
            }
            let pad = 0.05 * (w_hi - w_lo) + 1.0e-6 * (1.0 + smax(w_lo.abs(), w_hi.abs()));
            vary.cells
                .init(n, w_lo - pad, smin(w_hi + pad, 0.5 * (1.0 + w_hi)), 64);
        }
        let mut bounds = vec![-1.0_f64; stations.len().saturating_sub(1)];
        loop {
            let mut all = true;
            for (j, bj) in bounds.iter_mut().enumerate() {
                if *bj < 0.0 {
                    let (Some(a), Some(b)) = (stations.get(j), stations.get(j + 1)) else {
                        return Err(SolidError::sm_err(OP));
                    };
                    *bj = gn_certify_span(
                        k,
                        face_n,
                        &leader,
                        faces,
                        center,
                        a,
                        b,
                        r,
                        side,
                        self.order,
                        target,
                        if equal_area { Some(&mut vary) } else { None },
                    )?;
                }
                let within = *bj <= target;
                if !within {
                    all = false;
                }
            }
            if all {
                break;
            }
            if stations.len() > 4096 {
                // "the blend could not be certified within its tolerance"
                return Err(SolidError::sm_err(OP));
            }
            let mut next: Vec<GnStation> = Vec::new();
            let mut next_bounds: Vec<f64> = Vec::new();
            for (pair, &bj) in stations.windows(2).zip(&bounds) {
                let [a, b] = pair else { continue };
                next.push(*a);
                if bj <= target {
                    next_bounds.push(bj);
                    continue;
                }
                let mid = 0.5 * (a.s + b.s);
                let inside = mid > a.s && mid < b.s;
                if !inside {
                    return Err(SolidError::sm_err(OP));
                }
                let mut seed = [0.0_f64; 4];
                for (c, sc) in seed.iter_mut().enumerate() {
                    *sc = 0.5 * (param(a, c, 0, 0) + param(b, c, 0, 0));
                }
                let mut st =
                    gn_solve_station(k, faces, center, mid, seed, r, side, tr, solve_tol(mid)?)?;
                if equal_area {
                    gn_station_thumbweight(n, rr, &mut st)?;
                }
                next.push(st);
                next_bounds.push(-1.0);
                next_bounds.push(-1.0);
            }
            if let Some(last) = stations.last() {
                next.push(*last);
            }
            stations = next;
            bounds = next_bounds;
        }
        let mut bound = 0.0_f64;
        for &b in &bounds {
            bound = smax(bound, b);
        }
        // (the fork's SM_ERR_MESSAGE report of the spans and the bound is not logged)
        let spans = bounds.len();
        Ok((stations, bound, spans))
    }

    /// `ContactOrderAt`: the order of contact k of this section with face `rail` (0 or 1)
    /// at the station of the centre line's parameter `s`: the power at which the section
    /// stands off the face, `c s^k`. The coefficient of `t^(n+1)` in the standoff along the
    /// face's normal `N_i` is `[(T_o(0) - rho_i) - c_i,n+1 (-1/(1 - w))^(n+1)] . N_i`:
    /// `k = n + 1` where it is not zero, `n + 2` where it is (two planes).
    ///
    /// # Errors
    ///
    /// `SM_ERR_INVALID_INPUT` for a rail past 1 or an invalid generator; those of the
    /// station solve and the evaluations.
    #[allow(clippy::too_many_arguments)]
    pub fn contact_order_at(
        &self,
        k: &Kernel,
        faces: [&Surface; 2],
        center: &BSplineCurve,
        s: f64,
        seed: [f64; 4],
        rail_points: [Point3; 2],
        rail: usize,
    ) -> Result<usize> {
        if rail > 1 || !Self::is_valid(self.order, self.thumbweight) {
            return Err(SolidError::invalid_input(OP, "ContactOrderAt"));
        }
        let n = self.order;
        let m = center.evaluate_point(s)?;
        let mut r = [0.0_f64; 2];
        let mut side = [0.0_f64; 2];
        let rows = r
            .iter_mut()
            .zip(side.iter_mut())
            .zip(faces)
            .zip(rail_points)
            .zip(seed.chunks_exact(2));
        for ((((rf, sf), face), rho), uv) in rows {
            let (Some(&u), Some(&v)) = (uv.first(), uv.get(1)) else {
                return Err(SolidError::sm_err(OP));
            };
            let nrm = unit_normal(k, face, Point2::new(u, v))?;
            let off = m - rho;
            *rf = off.length();
            *sf = if off.dot(nrm) >= 0.0 { 1.0 } else { -1.0 };
        }
        let tr = Trunc::new(0, n + 1)?;
        let at = gn_solve_station(
            k,
            faces,
            center,
            s,
            seed,
            r,
            side,
            tr,
            1.0e-9 * smax(1.0, m.length()),
        )?;
        // the thumbweight of the station's own wedge: the ball's radii to its contacts from V
        let v = Vec3::new(at.v.a[0].at(0, 0), at.v.a[1].at(0, 0), at.v.a[2].at(0, 0));
        let w = self
            .thumbweight_at((foot(&at, 0, 0, 0) - v).dot(foot(&at, 1, 0, 0) - v) / (r[0] * r[1]))?;
        let (i, o) = (rail, 1 - rail);
        let mut t_o = Vec3::new(0.0, 0.0, 0.0);
        for kk in 0..=n {
            t_o += foot(&at, o, 0, kk) * if kk % 2 == 0 { 1.0 } else { -1.0 };
        }
        let p = surface_partials(
            k,
            if i == 0 { faces[0] } else { faces[1] },
            Point2::new(param(&at, 2 * i, 0, 0), param(&at, 2 * i + 1, 0, 0)),
            1,
        )?;
        let ni = vec3(p.at(1, 0)).cross(vec3(p.at(0, 1))).unitize()?;
        let scale = (-1.0 / (1.0 - w)).powf((n + 1) as f64);
        let coef = ((t_o - foot(&at, i, 0, 0)) - foot(&at, i, 0, n + 1) * scale).dot(ni);
        let size = foot(&at, 0, 0, 0).distance(foot(&at, 1, 0, 0));
        Ok(if coef.abs() > 1.0e-9 * smax(1.0, size) {
            n + 1
        } else {
            n + 2
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::extent2d::Extent2d;
    use crate::math::Vec2;
    use crate::surface::plane::Plane;

    fn plane(x_axis: Vec3, y_axis: Vec3) -> Surface {
        Surface::from(Plane::new(
            Vec3::ZERO,
            x_axis,
            y_axis,
            Vec2::new(1.0, 1.0),
            &Extent2d::from_bounds(-10.0, -10.0, 10.0, 10.0),
        ))
    }

    /// `EqualAreaThumbweight` for orders 2 to 4 over wedges from 0.3 to 3 radians: the
    /// values fork d01d7b0 gives (built locally, `-ffp-contract=off`), bit for bit; the
    /// wedges too shut for any thumbweight are refused, as there. At a square corner the
    /// section takes off what the circle does, `(1 - pi/4) r^2` (closed form).
    /// The fork's `EqualAreaThumbweight` at 7102bc5 (b9d6a72's bracket below -1, the
    /// closed-form integrand and the `sin`/`cos` target; built locally, `SMLib` and `SM_API`
    /// only, and run on these inputs): every value bit for bit. Before b9d6a72 the sharp
    /// wedges (0.05, 0.3, 1.0) had none.
    #[test]
    fn equal_area_thumbweights_are_the_fork_s() {
        let al = [
            0.05,
            0.3,
            1.0,
            core::f64::consts::FRAC_PI_2,
            2.0,
            2.5,
            3.0,
            3.1,
        ];
        let fork: [[f64; 8]; 3] = [
            [
                -40.158432198687606,
                -5.9593132523273784,
                -1.2200674332497643,
                -0.52521240458755813,
                -0.29072725674815347,
                -0.1505880430432826,
                -0.097163726970223796,
                -0.094769962349426001,
            ],
            [
                -41.536688806991449,
                -6.1966199386802323,
                -1.299175030152564,
                -0.58094928730219575,
                -0.3385247245324452,
                -0.19361792885821072,
                -0.13837023281846145,
                -0.13589468192029808,
            ],
            [
                -42.664956195284432,
                -6.3873279899518245,
                -1.3609708400917411,
                -0.62397886337619157,
                -0.37522338561877655,
                -0.22652960634299257,
                -0.16983693311526854,
                -0.16729661577416338,
            ],
        ];
        for (n, row) in (2..=4).zip(fork) {
            for (&alpha, want) in al.iter().zip(row) {
                let got = GnSectionFsg::equal_area_thumbweight(n, alpha).unwrap();
                assert_eq!(got.to_bits(), want.to_bits(), "n {n} alpha {alpha}: {got}");
            }
        }
        for alpha in al {
            assert!(GnSectionFsg::equal_area_thumbweight(1, alpha).is_err());
            assert!(GnSectionFsg::equal_area_thumbweight(5, alpha).is_err());
        }
        assert!(GnSectionFsg::equal_area_thumbweight(2, 0.0).is_err());
        assert!(GnSectionFsg::equal_area_thumbweight(2, PI).is_err());

        // a square corner: half the integral is the area over r^2, 1 - pi/4 (closed form)
        for n in 2..=GN_MAX_ORDER {
            let w = GnSectionFsg::equal_area_thumbweight(n, core::f64::consts::FRAC_PI_2).unwrap();
            let area = 0.5 * gn_area_sum(n, w);
            assert!((area - (1.0 - PI / 4.0)).abs() < 1e-12, "n {n}: {area}");
        }

        // at a station: the radii at a right angle (cos 0) give the square corner's w; a
        // thumbweight given is kept
        let g = GnSectionFsg::new(2, GN_DEFAULT_THUMBWEIGHT);
        assert_eq!(
            g.thumbweight_at(0.0).unwrap().to_bits(),
            (-0.52521240458755813_f64).to_bits()
        );
        assert_eq!(g.thumbweight, GN_THUMBWEIGHT_EQUAL_AREA);
        assert_eq!(GnSectionFsg::new(2, 0.3).thumbweight_at(0.7).unwrap(), 0.3);
        // radii that coincide or oppose: no wedge
        assert!(g.thumbweight_at(1.0).is_err());
        assert!(g.thumbweight_at(-1.0).is_err());
    }

    #[test]
    fn valid_orders_and_thumbweights() {
        assert!(GnSectionFsg::is_valid(2, 0.0));
        assert!(GnSectionFsg::is_valid(4, 0.999));
        assert!(!GnSectionFsg::is_valid(1, 0.0));
        assert!(!GnSectionFsg::is_valid(5, 0.0));
        assert!(!GnSectionFsg::is_valid(2, 1.0));
        // the fork's b9d6a72 (its IsValid at 7102bc5, run): any finite w below 1, and the
        // equal area default 2
        assert!(GnSectionFsg::is_valid(2, -0.1));
        assert!(GnSectionFsg::is_valid(2, -1.0));
        assert!(GnSectionFsg::is_valid(2, -5.0));
        assert!(GnSectionFsg::is_valid(2, -1e300));
        assert!(GnSectionFsg::is_valid(2, -2.0));
        assert!(GnSectionFsg::is_valid(2, GN_THUMBWEIGHT_EQUAL_AREA));
        assert_eq!(GN_THUMBWEIGHT_EQUAL_AREA, 2.0);
        assert!(!GnSectionFsg::is_valid(2, f64::NEG_INFINITY));
        assert!(!GnSectionFsg::is_valid(2, f64::NAN));
        let g = GnSectionFsg::new(3, 0.25);
        assert_eq!(g.certified_bound, -1.0);
        assert_eq!(g.certified_spans, 0);
    }

    /// Two planes: the feet are linear in tau, so `c_i,n+1 = 0` and `T_o(0)` is the edge's
    /// point, on face i: the `t^(n+1)` coefficient vanishes and the contact is of order
    /// n + 2 on both faces (the fork's `ContactOrderAt` note), at a square corner and at
    /// 60 degrees, for every order and a thumbweight.
    #[test]
    fn two_planes_meet_the_section_at_order_n_plus_2() {
        let k = Kernel::default();
        let r = 1.5;
        let (c, s60) = (0.5_f64, 3.0_f64.sqrt() / 2.0);
        // face 1: z = 0 (normal +z); face 2: the plane through the y axis at angle a to it
        for (dir2, center) in [
            (Vec3::Z, Vec3::new(r, 0.0, r)),
            // the second face along (c, 0, s60): the ball's centre on the bisector at
            // distance r from both
            (
                Vec3::new(c, 0.0, s60),
                Vec3::new(r / (30.0_f64.to_radians()).tan(), 0.0, r),
            ),
        ] {
            let f1 = plane(Vec3::X, Vec3::Y);
            let f2 = plane(Vec3::Y, dir2);
            let leader =
                BSplineCurve::create_line_segment(3, center, center + Vec3::new(0.0, 4.0, 0.0))
                    .unwrap();
            let s = 0.5;
            let m = center + Vec3::new(0.0, 2.0, 0.0);
            let rho1 = Vec3::new(m.x, m.y, 0.0);
            // the foot on face 2 along its direction, at the centre's projection
            let t2 = m.x * dir2.x + m.z * dir2.z;
            let rho2 = Vec3::new(t2 * dir2.x, m.y, t2 * dir2.z);
            assert!(((m - rho1).length() - r).abs() < 1e-12);
            assert!(((m - rho2).length() - r).abs() < 1e-12);
            let seed = [rho1.x, rho1.y, m.y, t2];
            for n in 2..=GN_MAX_ORDER {
                for w in [0.0, 0.3] {
                    let g = GnSectionFsg::new(n, w);
                    for rail in 0..2 {
                        let order = g
                            .contact_order_at(&k, [&f1, &f2], &leader, s, seed, [rho1, rho2], rail)
                            .unwrap();
                        assert_eq!(order, n + 2, "n {n} w {w} rail {rail}");
                    }
                }
            }
        }
    }
}
