//! A B-spline curve extended past one end by an arc length, exactly:
//! [`BSplineCurve::extended_by_distance`].
//!
//! Not upstream: it replaces the port of `N_CrvExtendByDist`, which estimated the distance and
//! fitted its extension. A curve's net is one row, so this is the surface extension's row
//! construction ([`Surface::extended_by_distance`](crate::surface::surface::Surface)) on that
//! row, with the same distance solve: Natural continues the end span's polynomial (on the
//! homogeneous points, so a rational curve stays on its conic), Tangent mirrors the end piece,
//! Curvature mirrors it with its second point from the C2 condition. The source's points and
//! knots are kept bit for bit; the arc length past the end is measured as a certified
//! enclosure (`curve_length`), to `1e-9` of the distance; a weight not positive or a length
//! not measurable to that stop is refused ([`SolidError::ExtendRefused`]) naming a distance
//! that works.
//!
//! Deviations: none (no upstream counterpart). No recursion.

use super::{not_exact, op_err, solve, Continuity, Grown, Local, Made};
use crate::error::{Result, SolidError};
use crate::nurbs::b_spline_curve::BSplineCurve;
use crate::nurbs::crv_adv::p1::CurveEnd;
use crate::nurbs::crv_basic::Curve as Net;
use crate::nurbs::math::CPoint;
use crate::surface::side::{curve_length, Across};

const OP: &str = "surface::extend::BSplineCurve::extended_by_distance";

/// An extended curve and its certified bounds.
#[derive(Debug, Clone)]
pub struct ExtendedCurve {
    /// The extension: the source's net and knots, bit for bit, then the extension's.
    pub curve: BSplineCurve,
    /// How far the part past the end is from the construction's definition: a certified bound
    /// on the distance, from the rounding of making its net.
    pub extension: f64,
    /// The arc length from the source's end to the result's end: its certified enclosure and
    /// the enclosure's parts.
    pub across: Across,
}

impl BSplineCurve {
    /// The curve continued `distance` (arc length) past its end `end` as `how` says, exactly
    /// (module docs): the result, its construction bound and its arc length's enclosure,
    /// within `max(1e-9 distance, the enclosure's width)` of `distance`.
    ///
    /// # Errors
    ///
    /// - `InvalidInput`: a distance not finite and positive; a net of degree 0, too small to
    ///   extend, or not clamped at the end;
    /// - `SM_ERR` for a curve without a net;
    /// - [`SolidError::ExtendRefused`]: a weight would not be positive, or the length at the
    ///   distance not measurable to the stop, with a distance that works when one was found;
    /// - `NotConverged`: the distance not bracketed or not met.
    pub fn extended_by_distance(
        &self,
        end: CurveEnd,
        how: Continuity,
        distance: f64,
    ) -> Result<ExtendedCurve> {
        if !(distance.is_finite() && distance > 0.0) {
            return Err(SolidError::invalid_input(
                OP,
                "the distance is not finite and positive",
            ));
        }
        let net = self.nurb.as_ref().ok_or_else(|| SolidError::sm_err(OP))?;
        let high = end == CurveEnd::End;
        let row = net.cpts().iter().map(super::homogeneous).collect();
        let source = Local::from_rows(
            vec![row],
            net.knots().to_vec(),
            net.degree(),
            net.is_rational(),
            high,
        )?;
        let joint = if high { source.end()? } else { -source.end()? };
        let made = |e: f64| -> Result<Option<Made<Net>>> {
            let grown = match source.extended(how, e)? {
                Grown::Made(g) => g,
                Grown::WeightNotPositive => return Ok(None),
                Grown::OverBudget => return Ok(Some((None, f64::INFINITY))),
            };
            // a net that overflowed is no extension to measure: unproved, not refused
            if !grown.is_finite() {
                return Ok(Some((None, f64::INFINITY)));
            }
            let bound = grown.bound();
            Ok(Some((Some(grown.curve(high)?), bound)))
        };
        // a length past the measure's leaf budget is no measurement, as a leaf not converged
        let across = |m: &Option<Net>, rel: f64| {
            let length = m.as_ref().map(|c| {
                let (lo, hi) = c.param_bounds();
                if high {
                    curve_length(c, joint, hi, rel)
                } else {
                    curve_length(c, lo, joint, rel)
                }
            });
            match length {
                None | Some(Err(SolidError::NotConverged { .. })) => Ok(None),
                Some(Ok(a)) => Ok(Some(a)),
                Some(Err(e)) => Err(e),
            }
        };
        let ((curve, extension), across) = solve(OP, &made, &across, &source, distance)?;
        let curve = curve.ok_or_else(op_err)?;
        // a bound that is not finite certifies nothing of the extension: not exact
        if !extension.is_finite() {
            return Err(not_exact(OP, None));
        }
        if !keeps(net, &curve, high) {
            return Err(SolidError::port_invariant(
                OP,
                "the extension moved the source's net",
            ));
        }
        Ok(ExtendedCurve {
            curve: Self {
                nurb: Some(curve),
                knot_type: self.knot_type,
                ..Self::new(self.dim)
            },
            extension,
            across,
        })
    }
}

/// Whether `out` holds `net`'s points and knots, bit for bit, at its kept end (`high`: the
/// extension is past the knots' high end): every knot but the extended end's last copy.
fn keeps(net: &Net, out: &Net, high: bool) -> bool {
    let (ka, kb) = (net.knots(), out.knots());
    let (pa, pb) = (net.cpts(), out.cpts());
    let (Some(keep), Some(shift), Some(moved)) = (
        ka.len().checked_sub(1),
        kb.len().checked_sub(ka.len()),
        pb.len().checked_sub(pa.len()),
    ) else {
        return false;
    };
    let same = |x: &CPoint, y: &CPoint| x.x == y.x && x.y == y.y && x.z == y.z && x.w == y.w;
    if high {
        ka.iter().take(keep).eq(kb.iter().take(keep)) && pa.iter().zip(pb).all(|(x, y)| same(x, y))
    } else {
        ka.iter().skip(1).eq(kb.iter().skip(shift + 1))
            && pa
                .iter()
                .zip(pb.iter().skip(moved))
                .all(|(x, y)| same(x, y))
    }
}

#[cfg(test)]
mod tests;
