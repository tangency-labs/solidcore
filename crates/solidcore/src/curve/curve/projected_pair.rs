//! Whole-domain construction of a single projected curve-pair maximum.
//!
//! Not upstream. The branch-and-bound invariant follows Neumaier (2004),
//! sections 11 and 12: only a feasible point supplies a lower bound, and every
//! parameter rectangle must have an upper bound within the requested tolerance.
//! Deviations: iterative subdivision; at most 32,768 boxes, depth 48 and degree
//! 18. Unsupported definitions and unresolved bounds return a typed refusal.

use super::p1::{affine_definition, AffineDefinition};
use super::{Curve, CurveKind};
use crate::blend::fillet_section::certificate::{curve_derivs_over, NurbsCurve};
use crate::blend::fillet_section::interval::Ivl;
use crate::blend::fillet_section::series::V3;
use crate::curve::hermite_curve::HermiteCurve;
use crate::error::{Result, SolidError};
use crate::geom::extent1d::Extent1d;
use crate::geom::solution_array::{Solution, SolutionArray, SolutionObject, SolutionType};
use crate::kernel::Kernel;
use crate::math::Vec3;
use crate::nurbs::math::{NOW, NOZ};

const OP: &str = "curve::curve::projected_pair";
const MAX_BOXES: usize = 32_768;
const MAX_DEPTH: usize = 48;
// The Taylor helper forms factorials in f64: through 18! they are exact integers.
const MAX_DEGREE: usize = 18;

fn refused(reason: &'static str) -> SolidError {
    SolidError::Unsupported {
        operation: OP,
        reason,
    }
}

pub(super) enum Definition<'a> {
    Nurbs(NurbsCurve),
    Line(AffineDefinition),
    Hermite(&'a HermiteCurve),
}

#[derive(Clone, Copy)]
pub(super) struct Piece {
    pub(super) range: Ivl,
    pub(super) knot: Option<usize>,
}

fn vector([x, y, z]: [Ivl; 3]) -> V3<Ivl> {
    V3 { x, y, z }
}

impl<'a> Definition<'a> {
    pub(super) fn of(curve: &'a Curve, domain: Extent1d) -> Result<Self> {
        match curve.kind() {
            CurveKind::Line(_) => affine_definition(curve, domain)
                .map(Self::Line)
                .ok_or_else(|| refused("line definition cannot be enclosed")),
            CurveKind::Hermite(h) => Ok(Self::Hermite(h)),
            CurveKind::BSpline(bs) => Self::nurbs(bs),
            CurveKind::Ellipse(e) => {
                // Validate size/degree before the identity checker decomposes the net.
                let definition = Self::nurbs(e.base())?;
                if !crate::curve::from_parameters::ellipse_form_is_exact(e)? {
                    return Err(refused(
                        "ellipse form does not certify its analytic support",
                    ));
                }
                Ok(definition)
            }
            // An analytic/procedural kind's cached form is not an identity proof.
            _ => Err(refused(
                "projected pair maximum needs an exact curve-definition enclosure",
            )),
        }
    }

    fn nurbs(bs: &crate::nurbs::b_spline_curve::BSplineCurve) -> Result<Self> {
        let raw = bs
            .gw_nurb()
            .ok_or_else(|| refused("missing curve definition"))?;
        let p = raw.degree();
        let n = raw.cpts().len();
        if p == 0
            || p > MAX_DEGREE
            || n <= p
            || n > MAX_BOXES
            || n.checked_add(p).and_then(|x| x.checked_add(1)) != Some(raw.knots().len())
            || !matches!(bs.dim(), 2 | 3)
            || raw.knots().iter().any(|x| !x.is_finite())
            || raw.knots().windows(2).any(|x| x.first() > x.last())
        {
            return Err(refused("curve definition exceeds the enclosure contract"));
        }
        // These controls are already homogeneous. Do not divide by w
        // and multiply back: that would certify a different rounded net.
        let mut pw = Vec::new();
        for cp in raw.cpts() {
            let w = if cp.w == NOW { 1.0 } else { cp.w };
            let z = if cp.z == NOZ { 0.0 } else { cp.z };
            if !(w > 0.0 && w.is_finite() && cp.x.is_finite() && cp.y.is_finite() && z.is_finite())
            {
                return Err(refused(
                    "curve definition has no finite positive-weight enclosure",
                ));
            }
            pw.extend([cp.x, cp.y, z, w]);
        }
        Ok(Self::Nurbs(NurbsCurve {
            p,
            n,
            u: raw.knots().to_vec(),
            pw,
        }))
    }

    pub(super) fn pieces(&self, domain: Extent1d) -> Result<Vec<Piece>> {
        let range = Ivl::new(domain.min(), domain.max());
        let Self::Nurbs(net) = self else {
            return Ok(vec![Piece { range, knot: None }]);
        };
        let mut pieces = Vec::new();
        for (knot, ends) in net.u.windows(2).enumerate().skip(net.p).take(net.n - net.p) {
            let (&lo, &hi) = ends
                .first()
                .zip(ends.last())
                .ok_or_else(|| refused("missing knot span"))?;
            let a = lo.max(range.lo);
            let b = hi.min(range.hi);
            if lo < hi && a <= b {
                pieces.push(Piece {
                    range: Ivl::new(a, b),
                    knot: Some(knot),
                });
            }
            if pieces.len() > MAX_BOXES {
                return Err(refused(
                    "curve has too many knot spans for the enclosure budget",
                ));
            }
        }
        if pieces.is_empty() {
            return Err(refused("curve interval has no knot piece"));
        }
        // The public B-spline evaluator returns first/last controls at natural
        // endpoints, including nonclamped nets. Enclose these isolated values
        // separately from the one-sided polynomial pieces.
        for (endpoint, knots) in [
            (net.u.get(net.p), net.u.get(..=net.p)),
            (net.u.get(net.n), net.u.get(net.n..)),
        ] {
            if let (Some(&endpoint), Some(knots)) = (endpoint, knots) {
                if range.lo <= endpoint
                    && endpoint <= range.hi
                    && !knots.iter().all(|&knot| knot == endpoint)
                {
                    pieces.push(Piece {
                        range: Ivl::from(endpoint),
                        knot: None,
                    });
                }
            }
        }
        Ok(pieces)
    }

    pub(super) fn jet(&self, range: Ivl, knot: Option<usize>) -> Result<(V3<Ivl>, V3<Ivl>)> {
        match self {
            Self::Nurbs(net) => {
                if range.lo == range.hi && knot.is_none() {
                    let control = if net.u.get(net.p) == Some(&range.lo) {
                        net.pw.get(..4)
                    } else if net.u.get(net.n) == Some(&range.lo) {
                        net.pw.get(net.pw.len().saturating_sub(4)..)
                    } else {
                        None
                    };
                    if let Some([x, y, z, w]) = control {
                        let weight = Ivl::from(*w);
                        return Ok((
                            V3 {
                                x: Ivl::from(*x) / weight,
                                y: Ivl::from(*y) / weight,
                                z: Ivl::from(*z) / weight,
                            },
                            V3::from_vec3(Vec3::new(0.0, 0.0, 0.0)),
                        ));
                    }
                }
                let jet = curve_derivs_over(net, range, 1, knot)?;
                Ok((
                    *jet.first()
                        .ok_or_else(|| refused("missing curve value enclosure"))?,
                    *jet.get(1)
                        .ok_or_else(|| refused("missing curve derivative enclosure"))?,
                ))
            }
            Self::Line(line) => {
                let (point, derivative) = line.jet_over(range);
                Ok((vector(point), vector(derivative)))
            }
            Self::Hermite(h) => {
                let (c, d) = h.power_basis();
                let a = V3::from_vec3(h.start_point());
                let b = V3::from_vec3(h.start_derivative());
                let c = V3::from_vec3(c);
                let d = V3::from_vec3(d);
                Ok((
                    add(
                        a,
                        scale(add(b, scale(add(c, scale(d, range)), range)), range),
                    ),
                    add(
                        b,
                        scale(
                            add(scale(c, Ivl::from(2.0)), scale(d, Ivl::from(3.0) * range)),
                            range,
                        ),
                    ),
                ))
            }
        }
    }
}

fn add(a: V3<Ivl>, b: V3<Ivl>) -> V3<Ivl> {
    V3 {
        x: a.x + b.x,
        y: a.y + b.y,
        z: a.z + b.z,
    }
}

pub(super) fn sub(a: V3<Ivl>, b: V3<Ivl>) -> V3<Ivl> {
    V3 {
        x: a.x - b.x,
        y: a.y - b.y,
        z: a.z - b.z,
    }
}

fn scale(a: V3<Ivl>, b: Ivl) -> V3<Ivl> {
    V3 {
        x: a.x * b,
        y: a.y * b,
        z: a.z * b,
    }
}

fn cross(a: V3<Ivl>, b: V3<Ivl>) -> V3<Ivl> {
    V3 {
        x: a.y * b.z - a.z * b.y,
        y: a.z * b.x - a.x * b.z,
        z: a.x * b.y - a.y * b.x,
    }
}

pub(super) fn project(v: V3<Ivl>, normal: V3<Ivl>) -> V3<Ivl> {
    cross(cross(normal, v), normal)
}

pub(super) fn dot(a: V3<Ivl>, b: V3<Ivl>) -> Ivl {
    (a.x * b.x + a.y * b.y) + a.z * b.z
}

fn square(x: Ivl) -> Ivl {
    Ivl::checked(x.mig() * x.mig(), x.mag() * x.mag())
}

pub(super) fn norm_squared(v: V3<Ivl>) -> Ivl {
    (square(v.x) + square(v.y)) + square(v.z)
}

#[derive(Clone, Copy)]
pub(super) struct Cell {
    pub(super) a: Piece,
    pub(super) b: Piece,
    pub(super) depth: usize,
}

pub(super) struct Query<'a> {
    pub(super) kernel: &'a Kernel,
    pub(super) curves: [&'a Curve; 2],
    pub(super) definitions: [Definition<'a>; 2],
    pub(super) normal: Vec3,
    pub(super) tolerance: f64,
}

#[derive(Clone, Copy)]
pub(super) struct Candidate {
    pub(super) parameters: [f64; 2],
    pub(super) reported: f64,
    pub(super) lower: f64,
    pub(super) upper: f64,
}

impl Query<'_> {
    pub(super) fn candidate(&self, s: f64, t: f64) -> Result<Candidate> {
        let [a, b] = &self.definitions;
        // Point evaluation uses the definition's own knot-side convention,
        // independently of the left/right box whose endpoint nominated it.
        let (pa, _) = a.jet(Ivl::from(s), None)?;
        let (pb, _) = b.jet(Ivl::from(t), None)?;
        let value = norm_squared(project(sub(pa, pb), V3::from_vec3(self.normal))).sqrt();
        let [ca, cb] = self.curves;
        let delta = ca.evaluate_point(self.kernel, s)? - cb.evaluate_point(self.kernel, t)?;
        let reported = self.normal.cross(delta).cross(self.normal).length();
        if !(reported.is_finite()
            && value.lo.is_finite()
            && value.hi.is_finite()
            && (value - Ivl::from(reported)).mag() <= self.tolerance)
        {
            return Err(refused(
                "candidate value cannot be enclosed within the requested tolerance",
            ));
        }
        Ok(Candidate {
            parameters: [s, t],
            reported,
            lower: value.lo.min(reported),
            upper: value.hi.max(reported),
        })
    }

    pub(super) fn bound(&self, cell: Cell) -> Result<Ivl> {
        let [a, b] = &self.definitions;
        let n = V3::from_vec3(self.normal);
        let (pa, da) = a.jet(cell.a.range, cell.a.knot)?;
        let (pb, db) = b.jet(cell.b.range, cell.b.knot)?;
        let d = project(sub(pa, pb), n);
        let direct = norm_squared(d);
        let s = cell.a.range.mid();
        let t = cell.b.range.mid();
        let (ca, _) = a.jet(Ivl::from(s), cell.a.knot)?;
        let (cb, _) = b.jet(Ivl::from(t), cell.b.knot)?;
        let centered = norm_squared(project(sub(ca, cb), n))
            + Ivl::from(2.0) * dot(d, project(da, n)) * (cell.a.range - Ivl::from(s))
            - Ivl::from(2.0) * dot(d, project(db, n)) * (cell.b.range - Ivl::from(t));
        Ok(Ivl::new(
            direct.lo.max(centered.lo),
            direct.hi.min(centered.hi),
        ))
    }
}

pub(super) fn bounded_domain(curve: &Curve, interval: Extent1d) -> Result<Extent1d> {
    let natural = curve.natural_interval();
    let lo = natural.min().max(interval.min());
    let hi = natural.max().min(interval.max());
    if !(lo.is_finite()
        && hi.is_finite()
        && lo <= hi
        && interval.min().is_finite()
        && interval.max().is_finite())
    {
        return Err(refused(
            "projected maximum needs finite overlapping parameter domains",
        ));
    }
    Ok(Extent1d::from_min_max(lo, hi))
}

/// Construct a single maximum directly from feasible candidates and complete bounds.
pub(super) fn solve<O: SolutionObject>(
    kernel: &Kernel,
    curves: [&Curve; 2],
    intervals: [Extent1d; 2],
    normal: Vec3,
    tolerance: f64,
    handles: [Option<O>; 2],
    target: Option<f64>,
) -> Result<SolutionArray<O>> {
    kernel.check_cancelled(OP)?;
    if !(tolerance.is_finite()
        && tolerance >= 0.0
        && normal.x.is_finite()
        && normal.y.is_finite()
        && normal.z.is_finite())
    {
        return Err(refused("projected maximum has no finite error budget"));
    }
    let [ca, cb] = curves;
    let [ia, ib] = intervals;
    let da = bounded_domain(ca, ia)?;
    let db = bounded_domain(cb, ib)?;
    let a = Definition::of(ca, da)?;
    let b = Definition::of(cb, db)?;
    let pieces_a = a.pieces(da)?;
    let pieces_b = b.pieces(db)?;
    let query = Query {
        kernel,
        curves,
        definitions: [a, b],
        normal,
        tolerance,
    };
    let mut best = query.candidate(da.min(), db.min())?;
    let mut visited = 0;
    let mut global_upper = best.upper;
    for a in pieces_a {
        for &b in &pieces_b {
            let mut pending = vec![Cell { a, b, depth: 0 }];
            while let Some(cell) = pending.pop() {
                kernel.check_cancelled(OP)?;
                visited += 1;
                if visited > MAX_BOXES {
                    return Err(refused("projected maximum exhausted its enclosure budget"));
                }
                for s in [cell.a.range.lo, cell.a.range.mid(), cell.a.range.hi] {
                    for t in [cell.b.range.lo, cell.b.range.mid(), cell.b.range.hi] {
                        let candidate = query.candidate(s, t)?;
                        if candidate.lower > best.lower {
                            best = candidate;
                        }
                    }
                }
                let bound = query.bound(cell)?;
                let threshold = (Ivl::from(best.lower) + Ivl::from(tolerance)).lo;
                if bound.lo <= bound.hi && bound.sqrt().hi <= threshold {
                    global_upper = global_upper.max(bound.sqrt().hi);
                    continue;
                }
                if cell.depth == MAX_DEPTH {
                    return Err(refused(
                        "projected maximum reached its subdivision depth bound",
                    ));
                }
                let split_a = (cell.a.range.hi - cell.a.range.lo)
                    / da.length().max(f64::MIN_POSITIVE)
                    >= (cell.b.range.hi - cell.b.range.lo) / db.length().max(f64::MIN_POSITIVE);
                let range = if split_a { cell.a.range } else { cell.b.range };
                let mid = range.mid();
                if !(range.lo < mid && mid < range.hi) {
                    return Err(refused(
                        "no representable subdivision closes the maximum bound",
                    ));
                }
                let mut left = Cell {
                    depth: cell.depth + 1,
                    ..cell
                };
                let mut right = left;
                if split_a {
                    left.a.range.hi = mid;
                    right.a.range.lo = mid;
                } else {
                    left.b.range.hi = mid;
                    right.b.range.lo = mid;
                }
                pending.push(right);
                pending.push(left);
            }
        }
    }
    if let Some(target) = target {
        if !target.is_finite() {
            return Err(refused("projected maximum target must be finite"));
        }
        let cutoff = Ivl::from(target.abs()) - Ivl::from(tolerance);
        if best.lower < cutoff.hi {
            if global_upper < cutoff.lo {
                return Ok(SolutionArray::new());
            }
            return Err(refused("maximum enclosure straddles the incumbent cutoff"));
        }
    }
    let mut solution = Solution::new();
    solution.num_objects = 2;
    for (slot, handle) in solution.objects.iter_mut().zip(handles) {
        *slot = handle;
    }
    solution.num_variables = 2;
    solution.solution_type = SolutionType::SingleValue;
    let [s, t] = best.parameters;
    *solution.start.param_mut(0)? = s;
    *solution.start.param_mut(1)? = t;
    solution.start.solution_value = best.reported;
    let mut solutions = SolutionArray::new();
    solutions.add(&solution);
    Ok(solutions)
}
