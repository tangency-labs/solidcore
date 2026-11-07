//! `SmCurve`, part 1: the global point solver (`GlobalPointSolve`) with the local
//! solve of its curve/point global solver, the Hermite approximation of a curve by a
//! piecewise cubic B-spline (`ApproximateCurve`), tessellation by bisection, the
//! classification of a curve/curve range intersection, and the refinement of range
//! intersections.
//!
//! Upstream: `source/SMLib/src/SmCurve.cpp`, `source/SMLib/inc/SmCurve.h`.
//!
//! | Rust                                              | Upstream                                         |
//! |---------------------------------------------------|--------------------------------------------------|
//! | [`refine_cc_solutions`]                           | `sm_RefineCCSolutions` (static)                  |
//! | [`test_span_accuracy`], [`SpanAccuracy`]          | `sm_TestSpanAccuracy` (static)                   |
//! | [`Curve::approximate_curve`], [`Curve::approximate_curve_base`], [`CurveApproximation`], [`ApproximateCurveOptions`] | `ApproximateCurve` (virtual: the `SmCrvOnSurf` override), `SmCurve::ApproximateCurve` |
//! | [`Curve::classify_intersecting_range_solution`], [`TangentPoint`] | `SmCurve::ClassifyIntersectingRangeSolution` |
//! | `CpGlobalSolver::local_solve_branch` ([`GlobalSolve::local_solve`] of [`CpGlobalSolver`]) | `SmCPGlobalSolver::LocalSolve` |
//! | [`Curve::global_point_solve`], [`GlobalPointSolveOptions`] | `SmCurve::GlobalPointSolve`             |
//! | `Curve::global_point_solve_cached` (crate)        | `SmCurve::GlobalPointSolve` (caller's cache)     |
//! | [`Curve::tessellate_by_bisection`]                | `SmCurve::TessellateByBisection`                 |
//! | [`Curve::tessellate`], [`Curve::tessellate_base`] | `Tessellate` (virtual), `SmCurve::Tessellate`    |
//!
//! Deviations:
//! - **Projected affine extrema.** Exact affine spans use a closed-form local
//!   fraction and an outward-rounded objective certificate, including raw tiny
//!   knot spans. Interior full-multiplicity upper knots are isolated right-valued
//!   candidates beside a left span with its upper limit excluded; query bounds are intersected
//!   exactly. Certificates account for their own rounding separately from
//!   retention slack. Projected extrema refuse uncovered cache domains and
//!   propagate solve errors rather than publishing uncertified partial results.
//!   Native joins are omitted only after defining continuity/collinearity or
//!   one-sided objective-sign proofs, never from an approximate G1 label alone.
//! - **Projected endpoint extrema.** Eligible finite boundaries are candidates
//!   independently of the stationary-point search. Finding a stationary
//!   point does not establish that it beats the interval endpoints.
//! - **Out-params.** Outputs are returned: [`SpanAccuracy`], [`CurveApproximation`]
//!   (whose `achieved_tol` is `None` where upstream leaves `rdAchievedTol` unwritten,
//!   the degenerate curve), [`TangentPoint`] (written only for a range classified as
//!   degenerate), a [`Tessellation`] (all three arrays; upstream fills the ones it is
//!   given, each independently, and fails with `SM_ERR_INVALID_INPUT` when none is),
//!   and a [`SolutionArray`] for `GlobalPointSolve`.
//! - **Objects.** Solutions name the caller's handles: `GlobalPointSolve` takes the
//!   handle of the curve (named by the line drop and the tree solutions, upstream
//!   `this`) and of its analytic replacement ([`GlobalPointSolveOptions`]; upstream
//!   the cache's `SmLine` / `SmCircle`). `sm_RefineCCSolutions` resolves solution
//!   objects to curves through a caller-supplied lookup (upstream casts the object
//!   pointers; a failing cast is `SM_ERR_NULL_POINTER`).
//! - **Caches.** `GlobalPointSolve` builds the curve's cache with
//!   [`make_curve_cache`] (co-ported from unit `curve.curve_p5`) instead of looking it
//!   up in the cache store (see there); a failed build is `SM_ERR_NULL_POINTER`, as
//!   upstream's `NER` on the missing cache. The cache tree's owner object is the
//!   curve's handle, and [`CpGlobalSolver`] carries the curve
//!   ([`CpGlobalSolver::curve`]) that upstream reads back from the tree's cache.
//! - **Persistent solution.** `SmCPGlobalSolver::LocalSolve` reuses one
//!   `SmSolution` across its guesses: [`Curve::local_point_solve`] returns a fresh
//!   solution, merged into the persistent one exactly as upstream writes it (type,
//!   value and first parameter when an answer is found; the whole solution for the
//!   ray operations).
//! - **Ignored statuses.** Where upstream ignores a status (the span end
//!   evaluations and the closed-curve continuity of `LocalSolve`, `SolveIt` in
//!   non-projected `GlobalPointSolve`, `GetKnots` in `TessellateByBisection`, the local intersect
//!   of `sm_RefineCCSolutions`, the local solves and evaluations of
//!   `ClassifyIntersectingRangeSolution`, `Copy` and `EditParameterization` in
//!   `ApproximateCurve`, `AngleBetween`) the error is ignored and upstream's value
//!   kept (the previous one, `0` for a failed angle, [`Vec3::UNDEF`] for a
//!   never-written `SmVector3d`), except [`SolidError::Unported`], which is always
//!   propagated. `sm_GetProjDist` returns a failing status as the distance
//!   upstream; that status is the solution value here too.
//! - **Undefined behaviour.** An `ApproximateCurve` that builds no curve (upstream
//!   dereferences the null curve) and a failing `Copy` of a B-spline in the
//!   just-copy branch (upstream returns a null curve with success) are errors.
//!   Option vectors an operation reads without a `NULL` check are `InvalidInput`
//!   when missing.
//! - **Approximation work.** The Hermite approximation keeps upstream's depth-first
//!   span order, but stops with `SM_ERR` after 8192 subdivisions or when the native
//!   midpoint formula cannot produce a finite, distinct parameter. The subdivision
//!   count is total across break intervals: upstream's natural-interval / 1500
//!   stopping rule permits at most 2047 subdivisions for one in-domain interval.
//!   Across a partition of the natural interval, the total stays below roughly
//!   2 x 1500 plus the break count, rather than 2047 per break, leaving margin under
//!   the 8192 cap. Far-outside breaks can otherwise require trillions of
//!   spans. The bound limits the pending stack and generated polygon growth, not
//!   the supplied breaks (including duplicates), which are still visited in order.
//! - **Output allocation.** The separate [`MAX_APPROX_SPANS`] appended-span bound
//!   is retained and returns `InvalidInput` before unbounded control-polygon growth
//!   (upstream fuzz `curve.curve_p1` fuzz01255 aborted on a 6 GiB allocation).
//! - `SM_ASSERT`, `WARN` and the `SM_DEBUG_CODE` drawing blocks are not ported;
//!   `TessellateByBisection`'s unused `GetDegree` (its only consumer,
//!   `sm_FindTestParams`, ignores it) is not called.
//! - **Callees of unit `curve.curve_p5`** (`SmCurve::LocalPropertyAnalysis`, the
//!   evaluators behind `LocalCurveSolve` / `LocalCurveIntersect`, `SmCurve::IsArc`)
//!   are wired; [`SolidError::Unported`] is propagated from what they reach.
//! - **Statuses passed on (WR-663).** With the opt-in
//!   `Kernel::fixes.query_statuses_pass_on` the statuses SMLib drops are passed on, as the
//!   fork (tangency-labs/smlib 8311f8e): `sm_RefineCCSolutions` (its local intersections;
//!   its `Inversion`s not, tangency-labs/smlib 8ea6a20), the degenerate `ApproximateCurve`'s `EditParameterization`,
//!   `ClassifyIntersectingRangeSolution`, `SmCPGlobalSolver::LocalSolve` and
//!   `GlobalPointSolve`'s `SolveIt`.

use crate::curve::curve::knots_ignoring_status;
use crate::curve::curve::p2::{
    add_to_curve, evaluate_continuity, find_test_params, flush_crv, output, pop_stack,
    proj_signed_angle, push_point, stack_look, test_step, CpGlobalSolver, Curve, CurveKind,
    CurvePropertyType, TessStepResult, BSPLINE_CURVE_TYPE, CRV_ON_SURF_TYPE,
};
use crate::curve::curve::p3::{eval_n, local_curve_intersect, PointSolveOptions};
use crate::curve::curve::p4::proj_dist;
use crate::curve::curve::p5::{local_property_analysis, make_curve_cache};
use crate::curve::curve_cache::CurveCache;
use crate::curve::hermite_curve::HermiteCurve;
use crate::curve::line::TessellationRequest;
use crate::error::{ignore_or_pass_on, ignored, ignored_or_pass_on, passed_on, Result, SolidError};
use crate::geom::extent1d::Extent1d;
use crate::geom::geom_utility::p1::segment_point_angle_deg_to_perp;
use crate::geom::geom_utility::p2::ContinuityType;
use crate::geom::global_solver::{
    Branch, GlobalSolve, ObjectZoneTol, SolutionPairType, SolutionRequested, SolveRequest,
    SolverOperation,
};
use crate::geom::solution_array::{
    OperandObject, Solution, SolutionArray, SolutionObject, SolutionType, SortKey,
};
use crate::geom::tol::{are_same, is_zero, is_zero_to_tol};
use crate::geom::tree::{NodeData, Tree};
use crate::geom::vector3d::{Point3, Vec3};
use crate::kernel::Kernel;
use crate::math::{
    max2, BIG_DOUBLE, EFF_ZERO, EFF_ZERO_RAD, EFF_ZERO_SQ, EFF_ZERO_SQRT, PI, UNDEF_DOUBLE,
};
use crate::nurbs::b_spline_curve::{BSplineCurve, Tessellation};

const OP_REFINE: &str = "curve::curve::p1::refine_cc_solutions";
const OP_APPROX: &str = "curve::curve::p1::Curve::approximate_curve";

/// Most Hermite spans [`Curve::approximate_curve`] appends before it gives up (not
/// upstream). Upstream grows its control polygon without limit, so a curve whose
/// spans pass the accuracy test only when tiny (coordinates near `1e10` and beyond,
/// e.g. a break at `1e300`) runs until memory is exhausted. Upstream fixtures need
/// hundreds of spans at most; 2^20 spans are about 75 MB of control points.
pub const MAX_APPROX_SPANS: usize = 1 << 20;
const OP_LOCAL: &str = "curve::curve::p1::CpGlobalSolver::local_solve_branch";
const OP_GLOBAL: &str = "curve::curve::p1::Curve::global_point_solve";

/// `SM_CONTINUITY_ANGLE` (`SM_ANG_TOL_DEG`, degrees): `EvaluateContinuity`'s default.
const CONTINUITY_ANGLE: f64 = 0.5;

/// Maximum adaptive Hermite subdivisions in one approximation (see deviations).
const MAX_APPROXIMATION_SUBDIVISIONS: usize = 8192;

/// Entry `i` of a three-slot sample array (`d` past its end).
fn at3<T: Copy>(a: &[T; 3], i: usize, d: T) -> T {
    a.get(i).copied().unwrap_or(d)
}

/// Sets entry `i` of a three-slot sample array.
fn put3<T>(a: &mut [T; 3], i: usize, v: T) {
    if let Some(x) = a.get_mut(i) {
        *x = v;
    }
}

/// The first `N` derivatives of an evaluation (missing entries stay `UNDEF`).
fn first_n<const N: usize>(v: &[Vec3]) -> [Vec3; N] {
    let mut out = [Vec3::UNDEF; N];
    for (d, s) in out.iter_mut().zip(v) {
        *d = *s;
    }
    out
}

/// The `SmStatus` value of an error a `double` function returns as its value
/// (`sm_GetProjDist`'s `SER`).
fn status_value(e: &SolidError) -> f64 {
    match e {
        SolidError::Upstream { code, .. } => f64::from(*code),
        SolidError::InvalidInput { .. } | SolidError::PortInvariant { .. } => 1007.0,
        _ => 1001.0,
    }
}

// ---- sm_RefineCCSolutions -----------------------------------------------------------

/// `sm_RefineCCSolutions`: polishes the ends of the curve/curve range solutions by a
/// local intersection from each end (when the end deviates by more than `1e-6`),
/// keeping a result that deviates less and moves both parameters by less than
/// `0.002` of the range (as upstream, the end's move is measured from the range's start,
/// so a forward range's end is never refined, unless
/// `Kernel::fixes.refine_cc_range_ends` is set: WR-383). `curve_of` resolves a solution
/// object to its curve.
///
/// # Errors
///
/// `SM_ERR` for a range solution without two variables; `SM_ERR_NULL_POINTER` for
/// an object that is not a curve; those of the local intersection.
pub fn refine_cc_solutions<'c, O, F>(
    k: &Kernel,
    solutions: &mut SolutionArray<O>,
    curve_of: F,
) -> Result<()>
where
    O: SolutionObject,
    F: Fn(Option<O>) -> Option<&'c Curve>,
{
    let desired_precision = 1.0e-10; // for curve/curve intersection
    let shift_limit = 0.002; // Maximum relative parameter shift allowed.
    let pass_on = k.fixes.query_statuses_pass_on;

    for ii in 0..solutions.len() {
        let Some(sol) = solutions.get_mut(ii) else {
            break;
        };

        // Check only range solutions.
        if sol.solution_type != SolutionType::RangeOfValues {
            continue;
        }

        if sol.num_variables != 2 {
            return Err(SolidError::sm_err(OP_REFINE));
        }

        let crv1 = curve_of(sol.objects[0]).ok_or(SolidError::null_pointer(OP_REFINE))?;
        let crv2 = curve_of(sol.objects[1]).ok_or(SolidError::null_pointer(OP_REFINE))?;

        let natural_dom1 = crv1.natural_interval();
        let natural_dom2 = crv2.natural_interval();

        // Set up the intervals.
        let mut sol_ivl1 = Extent1d::new();
        sol_ivl1.add_value(sol.start.parameters[0]);
        sol_ivl1.add_value(sol.end.parameters[0]);
        let reverse1 = sol.start.parameters[0] > sol.end.parameters[0];

        let mut sol_ivl2 = Extent1d::new();
        sol_ivl2.add_value(sol.start.parameters[1]);
        sol_ivl2.add_value(sol.end.parameters[1]);
        let reverse2 = sol.start.parameters[1] > sol.end.parameters[1];

        for at_end in [false, true] {
            let end = if at_end { sol.end } else { sol.start };
            if end.solution_value <= 0.000001 {
                continue;
            }
            let guess1 = end.parameters[0];
            let guess2 = end.parameters[1];

            // (status ignored; passed on with the opt-in fix, tangency-labs/smlib 8311f8e)
            let found = ignored_or_pass_on(
                local_curve_intersect(
                    k,
                    crv1,
                    &natural_dom1,
                    crv2,
                    &natural_dom2,
                    desired_precision,
                    guess1,
                    guess2,
                ),
                pass_on,
            )?
            .flatten();

            if let Some(x) = found {
                if x.deviation < end.solution_value {
                    // Do not allow big shifts.  [140915b]
                    let inv1 = sol_ivl1.inversion(x.my_param);
                    let inv2 = sol_ivl2.inversion(x.other_param);
                    // The fork (tangency-labs/smlib 8ea6a20): a refined parameter outside the
                    // range is an answer, which the shift limit judges. Only a zero-length
                    // range gives no fraction, and then the end stays where it is.
                    if pass_on && (inv1.parameter.is_none() || inv2.parameter.is_none()) {
                        continue;
                    }
                    let mut frac1 = inv1.parameter.unwrap_or(UNDEF_DOUBLE);
                    let mut frac2 = inv2.parameter.unwrap_or(UNDEF_DOUBLE);
                    if reverse1 {
                        frac1 = 1.0 - frac1;
                    }
                    if reverse2 {
                        frac2 = 1.0 - frac2;
                    }

                    // Not upstream (WR-383, tangency-labs/solidcore#527): with the
                    // opt-in, an end is judged by its distance from 1 (the fraction is
                    // measured from the solution's start); upstream tests |frac| at both
                    if at_end && k.fixes.refine_cc_range_ends {
                        frac1 = 1.0 - frac1;
                        frac2 = 1.0 - frac2;
                    }
                    if frac1.abs() < shift_limit && frac2.abs() < shift_limit {
                        let dst = if at_end { &mut sol.end } else { &mut sol.start };
                        dst.solution_value = x.deviation;
                        dst.parameters[0] = x.my_param;
                        dst.parameters[1] = x.other_param;
                    }
                }
            }
        }
    }

    Ok(())
}

// ---- sm_TestSpanAccuracy ------------------------------------------------------------

/// Result of [`test_span_accuracy`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpanAccuracy {
    /// `rbSpanPassesTest`: the Hermite span is within the tolerance (or very short).
    pub passes: bool,
    /// `rdMaxDistToCrv`: the largest sampled distance (up to the first failing
    /// sample when the span fails).
    pub max_dist: f64,
}

/// `sm_TestSpanAccuracy`: whether the Hermite span between `last` and `curr`
/// (position, first derivative at the parameters `last_t`, `curr_t`) lies within
/// `approx_tol` of `curve` at 3 samples (7 when the span is forced to pass: shorter
/// than `SM_EFF_ZERO` scaled or than 1/1500 of the curve's interval). With
/// `match_parameterization` the distance is taken at equal relative parameters,
/// else to the closest Hermite point (a failed drop counts as `2 approx_tol`).
///
/// # Errors
///
/// Those of the curve evaluations; [`SolidError::Unported`] from the drop.
pub fn test_span_accuracy(
    k: &Kernel,
    curve: &Curve,
    last_t: f64,
    last_pv: [Point3; 2],
    curr_t: f64,
    curr_pv: [Point3; 2],
    approx_tol: f64,
    match_parameterization: bool,
) -> Result<SpanAccuracy> {
    // init output
    let mut out = SpanAccuracy {
        passes: false,
        max_dist: 0.0,
    };

    // locals
    let mut force_pass = false;
    let this_span = curr_t - last_t;
    let t_scale = 1.0 + curr_t.abs();
    let ivl = curve.natural_interval();
    let min_step = 1500.0;
    let mut num_test_points: u32 = 3;

    // when span is getting very small - just pass it
    if this_span < EFF_ZERO * t_scale {
        force_pass = true;
    }

    // Don't let the interval get too small
    if this_span < (ivl.length()) / min_step {
        force_pass = true;
    }

    // hermite cuve between interval end points
    let hermite = HermiteCurve::new(
        last_pv[0],
        last_pv[1] * this_span,
        curr_pv[0],
        curr_pv[1] * this_span,
        3,
    );
    let hermite_curve = Curve::from(hermite);

    if force_pass {
        // get a better error value if we're going to keep this step.
        num_test_points = 7;
    }

    for ii in 0..num_test_points {
        let normalized_param = (f64::from(ii) + 1.0) / (f64::from(num_test_points) + 1.0);
        let crv_param = last_t + normalized_param * this_span;

        // Evaluate the curve.  Then, if we're matching parameterization,
        // evaluate the Hermite curve at the same (relative) parameter,
        // else find the closest point on the Hermite curve.
        let crv_pnt = curve.evaluate_point(k, crv_param)?;

        let approx_dist = if match_parameterization {
            // Approx dist = DistanceBetween(Curve(param), Hermite(param))
            let herm_pnt = hermite.evaluate_point(normalized_param);
            herm_pnt.distance(crv_pnt)
        } else {
            // Approx dist = DistanceBetween(Curve(param), NearestPointOnHermite)
            let hermite_domain = Extent1d::from_min_max(0.0, 1.0);
            let r = hermite_curve.local_point_solve::<()>(
                k,
                &hermite_domain,
                SolverOperation::Minimize,
                crv_pnt,
                &PointSolveOptions::default(),
                normalized_param,
            );
            // save nearest distance (when drop point failed - load a big distance to
            // quit testing)
            match ignored(r)? {
                Some((true, sol)) => sol.start.solution_value,
                _ => 2.0 * approx_tol,
            }
        };

        // When ApproxDist is too big
        if approx_dist > approx_tol && !force_pass {
            // set output and quit testing
            out.passes = false;
            return Ok(out);
        }
        out.max_dist = max2(out.max_dist, approx_dist);
    }

    // arrive here when HermiteCurve is within ApproxTol of Curve at sample points or
    // when bForcePass == TRUE
    out.passes = true;
    Ok(out)
}

// ---- ApproximateCurve ---------------------------------------------------------------

/// The flags of [`Curve::approximate_curve`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ApproximateCurveOptions {
    /// `bOptCreateAnalytics` (default `true`): make an exact arc when the
    /// approximation lies on a circle.
    pub create_analytics: bool,
    /// `bOptMatchParameterization` (default `false`): measure the error at equal
    /// parameters instead of to the closest point.
    pub match_parameterization: bool,
    /// `bJustCopyBSplines` (default `false`): copy B-spline curves instead of
    /// approximating them.
    pub just_copy_bsplines: bool,
}

impl Default for ApproximateCurveOptions {
    fn default() -> Self {
        Self {
            create_analytics: true,
            match_parameterization: false,
            just_copy_bsplines: false,
        }
    }
}

/// Result of [`Curve::approximate_curve`].
#[derive(Debug, Clone)]
pub struct CurveApproximation {
    /// `rpNewBSplineCurve`: the approximation (a B-spline or a copy of this curve).
    pub curve: Curve,
    /// `rdAchievedTol`: the largest sampled distance; `None` where upstream does not
    /// write it (a degenerate curve).
    pub achieved_tol: Option<f64>,
}

impl Curve {
    /// `ApproximateCurve`: a piecewise cubic B-spline through the curve at the
    /// parameters `breaks` (at least 2), each break interval bisected until its
    /// Hermite span passes [`test_span_accuracy`] with `approx_tol_3d`, flushed by
    /// [`flush_crv`] (an exact arc when possible with `create_analytics`). A
    /// degenerate curve becomes a degenerate B-spline over `[breaks[0], last]`;
    /// B-splines are copied with `just_copy_bsplines` (tolerance `0`) and also when
    /// the approximation would have more control points.
    ///
    /// # Errors
    ///
    /// `SM_ERR` for fewer than 2 breaks, exhausted subdivision work or a midpoint
    /// that cannot progress; those of the evaluations, the degenerate curve,
    /// [`flush_crv`] and `Copy`; `InvalidInput` when nothing is built.
    pub fn approximate_curve(
        &self,
        k: &Kernel,
        breaks: &[f64],
        approx_tol_3d: f64,
        opts: &ApproximateCurveOptions,
    ) -> Result<CurveApproximation> {
        match self.kind() {
            // the SmCrvOnSurf override (its cache update is lost on this copy)
            CurveKind::CrvOnSurf(c) => {
                let mut c = c.clone_exact_shared();
                let (curve, gap) =
                    c.approximate_curve(k, breaks, approx_tol_3d, opts.match_parameterization)?;
                Ok(CurveApproximation {
                    curve: Curve::from(curve),
                    achieved_tol: Some(gap),
                })
            }
            // the SmCrvOnSurf override on the derived object (its Evaluate approximates)
            CurveKind::TangentField(c) => {
                let members = c.members();
                let mut c = c.base().clone();
                let (curve, gap) = c.approximate_curve_as(
                    k,
                    breaks,
                    approx_tol_3d,
                    opts.match_parameterization,
                    Some(members),
                )?;
                Ok(CurveApproximation {
                    curve: Curve::from(curve),
                    achieved_tol: Some(gap),
                })
            }
            // the SmCompositeCurve override (breaks and flags unused)
            CurveKind::Composite(c) => {
                let a = c.approximate_curve(k, approx_tol_3d)?;
                let curve = a.curve.ok_or_else(|| {
                    SolidError::invalid_input(OP_APPROX, "composite approximation is not valid")
                })?;
                Ok(CurveApproximation {
                    curve: Curve::from(curve),
                    achieved_tol: Some(a.max_gap),
                })
            }
            _ => self.approximate_curve_base(k, breaks, approx_tol_3d, opts),
        }
    }

    /// `SmCurve::ApproximateCurve` (the base-class body of
    /// [`Curve::approximate_curve`]).
    ///
    /// # Errors
    ///
    /// As [`Curve::approximate_curve`].
    pub fn approximate_curve_base(
        &self,
        k: &Kernel,
        breaks: &[f64],
        approx_tol_3d: f64,
        opts: &ApproximateCurveOptions,
    ) -> Result<CurveApproximation> {
        // low work - when asked, just copy BSplines
        if opts.just_copy_bsplines && self.is_kind_of(BSPLINE_CURVE_TYPE) {
            // (upstream ignores the status of Copy)
            let curve = self.copy(k)?;
            return Ok(CurveApproximation {
                curve,
                achieved_tol: Some(0.0),
            });
        }

        // need 2 breakpoints to define the curve interval
        let (Some(&first), Some(&last)) = (breaks.first(), breaks.last()) else {
            return Err(SolidError::sm_err(OP_APPROX));
        };
        if breaks.len() < 2 {
            return Err(SolidError::sm_err(OP_APPROX));
        }

        // low work - degenerate curves
        if self.is_degenerate(k, EFF_ZERO, None)? {
            let pnt = self.evaluate_point(k, first)?;
            let degenerate = BSplineCurve::create_degenerate_curve(self.dim(), pnt)?;
            let mut curve = Curve::from(degenerate);
            let ivl = Extent1d::from_min_max(first, last);
            // (status ignored; passed on with the opt-in fix, tangency-labs/smlib 8311f8e)
            ignore_or_pass_on(
                curve.edit_parameterization(k, &ivl),
                k.fixes.query_statuses_pass_on,
            )?;
            return Ok(CurveApproximation {
                curve,
                achieved_tol: None,
            });
        }

        // locals
        let mut stack: Vec<f64> = Vec::with_capacity(256 * 4 * 2);
        let mut ctrl_poly: Vec<Point3> = Vec::with_capacity(256);
        let mut knots: Vec<f64> = Vec::with_capacity(256);
        let mut max_deviation = 0.0;
        let mut subdivisions = 0;

        // for every interval defined by input BreakParams array
        for w in breaks.windows(2) {
            let (Some(&a), Some(&b)) = (w.first(), w.get(1)) else {
                continue;
            };
            let mut last_t = a;
            let mut curr_t = b;

            // interval start position and tangent (Last)
            let mut last_pv: [Vec3; 2] = eval_n(k, self, last_t, 1, true, true)?; // NonZeroTangent values
            {
                // walk to a nonZero Tangent evaluation
                let mut step = 0.0001;
                while last_pv[1].length_squared() < EFF_ZERO_SQ * 100000.0 {
                    let tmp = last_pv[0];
                    last_pv = eval_n(k, self, last_t + step * (curr_t - last_t), 1, true, true)?;
                    last_pv[0] = tmp;
                    step *= 2.0;
                    if step > 0.5 {
                        break;
                    }
                }
            }

            // interval end position and tangent (Curr)
            let mut curr_pv: [Vec3; 2] = eval_n(k, self, curr_t, 1, false, true)?; // NonZeroTangent values
            {
                let mut step = 0.0001;
                while curr_pv[1].length_squared() < EFF_ZERO_SQ * 100000.0 {
                    let tmp = curr_pv[0];
                    curr_pv = eval_n(k, self, curr_t - step * (curr_t - last_t), 1, true, true)?;
                    curr_pv[0] = tmp;
                    step *= 2.0;
                    if step > 0.5 {
                        break;
                    }
                }
            }

            // Place Curr Point [param, posXYZ, tanXYZ] on sStack
            push_point(curr_t, curr_pv, &mut stack);

            // For each span subdivide until our approximating curve is accurate enough.
            while !stack.is_empty() {
                // get last point's param, pos, and tan from stack - leave stack unmodified
                (curr_t, curr_pv) = stack_look(&stack)?;

                // Test if HermiteSpan is within ApproxTol of thisCurve interval
                let acc = test_span_accuracy(
                    k,
                    self,
                    last_t,
                    last_pv,
                    curr_t,
                    curr_pv,
                    approx_tol_3d,
                    opts.match_parameterization,
                )?;

                // When Span passes test (or If step size gets real small just accept it)
                if acc.passes || (last_t - curr_t).abs() < EFF_ZERO {
                    // save max deviation
                    max_deviation = max2(max_deviation, acc.max_dist);

                    // (not upstream) bound the growing curve: see MAX_APPROX_SPANS
                    if knots.len() > MAX_APPROX_SPANS {
                        return Err(SolidError::invalid_input(
                            OP_APPROX,
                            "approximation needs more than MAX_APPROX_SPANS spans",
                        ));
                    }

                    // Append Hermite interval to end of growing piecewise cubic bezier
                    add_to_curve(last_pv, curr_pv, last_t, curr_t, &mut ctrl_poly, &mut knots)?;

                    // for next iteration - move Curr param, pos, and tan to Last
                    last_t = curr_t;
                    last_pv = curr_pv;

                    // for next iteration - pop stack
                    pop_stack(&mut stack)?;
                } else {
                    if subdivisions >= MAX_APPROXIMATION_SUBDIVISIONS {
                        return Err(SolidError::sm_err(OP_APPROX));
                    }
                    subdivisions += 1;

                    // get param mid point and associated curve pos and tan
                    let mid_t = (last_t + curr_t) / 2.0;
                    if !mid_t.is_finite() || mid_t == last_t || mid_t == curr_t {
                        return Err(SolidError::sm_err(OP_APPROX));
                    }
                    let mid_pv: [Vec3; 2] = eval_n(k, self, mid_t, 1, true, true)?; // NonZeroTangent values
                    curr_t = mid_t;
                    curr_pv = mid_pv;

                    // push midPoint onto stach
                    push_point(curr_t, curr_pv, &mut stack);
                }
            }
        }

        // create pBSC curve from control polygon
        let bsc = flush_crv(
            approx_tol_3d,
            &mut ctrl_poly,
            &mut knots,
            opts.create_analytics,
            k.fixes.query_statuses_pass_on,
        )?
        .ok_or(SolidError::invalid_input(
            OP_APPROX,
            "no approximating curve built",
        ))?;
        let new_count = bsc.number_control_points()?;
        let mut out = CurveApproximation {
            curve: Curve::from(bsc),
            achieved_tol: Some(max_deviation),
        };

        // return a warning when input curve is better than the approximation
        if self.is_kind_of(BSPLINE_CURVE_TYPE)
            && !self.is_kind_of(CRV_ON_SURF_TYPE)
            && new_count > self.number_control_points()?
        {
            out.curve = self.copy(k)?;
        }

        Ok(out)
    }
}

// ---- ClassifyIntersectingRangeSolution ------------------------------------------------

/// The tangent intersection point of
/// [`Curve::classify_intersecting_range_solution`] (the fields of `rSolPoint` it
/// writes).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TangentPoint {
    /// `m_dSolutionValue`: the deviation at the better range end.
    pub solution_value: f64,
    /// `m_adParameters[0..2]`: the parameters of that end on this and the other curve.
    pub parameters: [f64; 2],
}

impl Curve {
    /// `ClassifyIntersectingRangeSolution`: whether the range intersection `sol` with
    /// `other` is degenerate (a tangent point, not a coincident stretch), and then the
    /// point (`None`: `rSolPoint` left unmodified). A point solution is degenerate.
    /// Otherwise, if a local minimum of the gap (well below the better end's
    /// deviation) lies inside the range, the gap is sampled from the better end at
    /// geometrically growing steps: a non-zero, monotonic, faster than linearly
    /// growing gap (or a sample that does not drop onto the other curve) is a
    /// tangent point; a gap that stays zero over a short range also is.
    ///
    /// # Errors
    ///
    /// Those of the local solves.
    pub fn classify_intersecting_range_solution<O: SolutionObject>(
        &self,
        k: &Kernel,
        other: &Self,
        sol: &Solution<O>,
    ) -> Result<(bool, Option<TangentPoint>)> {
        // low work - Solution already a point solution
        if sol.solution_type == SolutionType::SingleValue {
            return Ok((true, None));
        }

        // locals to run from best solution end to other end
        let this_curve = self;
        let ivl = self.natural_interval();
        let o_ivl = other.natural_interval();
        let (start_dev, end_dev, start_par, delta_par) =
            if sol.start.solution_value < sol.end.solution_value {
                (
                    sol.start.solution_value,
                    sol.end.solution_value,
                    [sol.start.parameters[0], sol.start.parameters[1]],
                    [
                        sol.end.parameters[0] - sol.start.parameters[0],
                        sol.end.parameters[1] - sol.start.parameters[1],
                    ],
                )
            } else {
                (
                    sol.end.solution_value,
                    sol.start.solution_value,
                    [sol.end.parameters[0], sol.end.parameters[1]],
                    [
                        sol.start.parameters[0] - sol.end.parameters[0],
                        sol.start.parameters[1] - sol.end.parameters[1],
                    ],
                )
            };

        // Do this only if it looks like a relative minimum on the interior.  [140915b]
        let mut this_dom = Extent1d::from_value(sol.start.parameters[0]);
        this_dom.add_value(sol.end.parameters[0]);
        let mut othr_dom = Extent1d::from_value(sol.start.parameters[1]);
        othr_dom.add_value(sol.end.parameters[1]);
        // (status ignored; the fork passes these statuses on, tangency-labs/smlib 8311f8e)
        let pass_on = k.fixes.query_statuses_pass_on;
        let solved = ignored_or_pass_on(
            self.local_curve_solve::<O>(
                k,
                &this_dom,
                other,
                &othr_dom,
                SolverOperation::Minimize,
                1.0,
                Some(end_dev),
                None,
                [this_dom.evaluate(0.5), othr_dom.evaluate(0.5)],
                [None, None],
            ),
            pass_on,
        )?;
        let Some((true, min_sol)) = solved else {
            return Ok((false, None)); // Did not find a minimum.
        };

        if min_sol.start.solution_value > start_dev / 10.0 {
            return Ok((false, None)); // Not much of a minimum.
        }

        // get Solution segment start points and gapVector
        let pnt0 = ignored_or_pass_on(this_curve.evaluate_point(k, start_par[0]), pass_on)?
            .unwrap_or(Vec3::UNDEF);
        let pnt1 = ignored_or_pass_on(other.evaluate_point(k, start_par[1]), pass_on)?
            .unwrap_or(Vec3::UNDEF);
        let xsect_gap = pnt1 - pnt0;

        // scale zero to current intersection point location
        let scaled_zero = EFF_ZERO * (1.0 + pnt0.max_dimension());

        // locals - gap function properties
        let mut curve_point = [Vec3::UNDEF; 3];
        let mut o_curve_point = [Vec3::UNDEF; 3];
        let mut gap = [Vec3::UNDEF; 3];
        let mut gap_length = [0.0; 3];

        // init classification booleans
        let mut b_zero = true;
        let mut b_mono = true;
        let mut b_gdiv = true;

        // init first SamplePoint array entries with the xSect StartPoint
        curve_point[0] = pnt0;
        o_curve_point[0] = pnt1;
        gap[0] = Vec3::ZERO; // set gap to exactly zero - see trick later on.
        gap_length[0] = 0.0;

        // init walk parameter (a little smaller than a fraction of a power of two)
        let mut delta = 1.0 / 16.0001;

        // init solution classification
        let mut is_degenerate = false;

        // until solution segment is classified
        let mut done = false;
        let mut ii: usize = 1;
        while !done && delta < 1.0 {
            // get the gapFunction sample array index
            let tgt: usize = if ii > 1 { 2 } else { 1 };

            // get the samplePoint Parameter Values
            let t1 = start_par[0] + delta * delta_par[0];
            let t2 = start_par[1] + delta * delta_par[1];

            // get ThisCurve sample point
            if let Some(p) = ignored_or_pass_on(this_curve.evaluate_point(k, t1), pass_on)? {
                put3(&mut curve_point, tgt, p);
            }

            // Find OtherCurve nearest point to sCurvePoint
            let dropped = ignored_or_pass_on(
                other.local_point_solve::<O>(
                    k,
                    &o_ivl,
                    SolverOperation::Normalize,
                    at3(&curve_point, tgt, Vec3::UNDEF),
                    &PointSolveOptions::default(),
                    t2,
                ),
                pass_on,
            )?;

            if let Some((true, s)) = dropped {
                // evaluate the nearest OtherCurve point
                let o_drop = s.start.parameters[0];
                if let Some(p) = ignored_or_pass_on(other.evaluate_point(k, o_drop), pass_on)? {
                    put3(&mut o_curve_point, tgt, p);
                }

                // get the Curve/Curve gapVector (minus the start gap)
                let g = (at3(&o_curve_point, tgt, Vec3::UNDEF)
                    - at3(&curve_point, tgt, Vec3::UNDEF))
                    - xsect_gap;
                put3(&mut gap, tgt, g);
                put3(&mut gap_length, tgt, g.length());
            } else {
                // Found a point in the intersection range on this curve which does
                // not drop along a normal to the other curve: a point solution.
                is_degenerate = true;

                // get out of this classification loop
                done = true;
                ii += 1;
                delta *= 2.0;
                continue;
            }

            // once we have 3 sample points
            if ii >= 2 {
                let tgt1 = (tgt + 3 - 1) % 3;
                let tgt2 = (tgt + 3 - 2) % 3;
                let len = at3(&gap_length, tgt, 0.0);
                let len1 = at3(&gap_length, tgt1, 0.0);
                let len2 = at3(&gap_length, tgt2, 0.0);

                // 1st: is xSectCurve gapSize zero (coincident curves)
                if len > scaled_zero {
                    b_zero = false;
                }

                // 2nd: is xSectCurve gap function monotonic
                if at3(&gap, tgt, Vec3::UNDEF).dot(at3(&gap, tgt1, Vec3::UNDEF)) <= 0.0 {
                    b_mono = false;
                }

                // 3rd: is xSectCurve gap function geometrically diverging
                if len == 0.0 || len < len1 + 1.0001 * (len1 - len2) {
                    b_gdiv = false;
                }

                // (4th, the linear divergence test, sets nothing upstream)

                // now classify this XSectCurve GapFunction.
                if !b_zero && b_mono && b_gdiv {
                    is_degenerate = true;
                    done = true;
                } else if b_zero && b_mono && b_gdiv {
                    // try another point
                    done = false;

                    // move current guess to last guess
                    curve_point[1] = curve_point[2];
                    o_curve_point[1] = o_curve_point[2];
                    gap[1] = gap[2];
                    gap_length[1] = gap_length[2];
                } else {
                    // this side of the curve can have a coincident region
                    is_degenerate = false;
                    done = true;
                }
            }
            ii += 1;
            delta *= 2.0;
        }

        // check to see if interval was classified
        if !done && b_gdiv {
            // range diverged but never got larger than Scaled Zero:
            // let's try throwing away short intervals
            is_degenerate = delta_par[0].abs() < ivl.length() / 10.0
                && delta_par[1].abs() < o_ivl.length() / 10.0;
        }

        // when soluton is classified as degenerate
        if is_degenerate {
            return Ok((
                true,
                Some(TangentPoint {
                    solution_value: start_dev,
                    parameters: start_par,
                }),
            ));
        }
        Ok((false, None))
    }
}

// ---- SmCPGlobalSolver::LocalSolve -------------------------------------------------------

/// Merges a [`Curve::local_point_solve`] result into the persistent solution of
/// `LocalSolve` the way upstream writes its `rSolution`.
fn merge_point_solution<O: Copy>(
    dst: &mut Solution<O>,
    found: bool,
    src: &Solution<O>,
    op: SolverOperation,
) {
    if !found {
        return;
    }
    if op == SolverOperation::Rayfire || op == SolverOperation::Signed3dDirectedMinimize {
        *dst = *src;
        return;
    }
    dst.solution_type = src.solution_type;
    dst.start.solution_value = src.start.solution_value;
    dst.start.parameters[0] = src.start.parameters[0];
}

// Affine definitions are solved in a local fraction, not by driving a stationary
// residual to zero in a raw knot interval with only a few thousand f64 values.
use crate::blend::fillet_section::interval::Ivl;

/// Report that no representable parameter was certified to its own objective enclosure.
fn affine_refusal() -> SolidError {
    SolidError::Upstream {
        operation: "curve::curve::p1::projected_affine_certificate",
        code: crate::error::SM_ERR_NOT_WITHIN_TOLERANCE,
        entity: None,
    }
}

// MSRV 1.85: f64::next_up/next_down were stabilized in Rust 1.86.
/// Return the next representable parameter toward positive infinity on Rust 1.85.
fn next_parameter_up(value: f64) -> f64 {
    if value.is_nan() || value == f64::INFINITY {
        value
    } else if value == 0.0 {
        f64::from_bits(1)
    } else if value > 0.0 {
        f64::from_bits(value.to_bits() + 1)
    } else {
        f64::from_bits(value.to_bits() - 1)
    }
}

/// Return the next representable parameter toward negative infinity on Rust 1.85.
fn next_parameter_down(value: f64) -> f64 {
    -next_parameter_up(-value)
}

/// Enclose each coordinate of an exactly specified binary64 point.
fn iv_point(p: Vec3) -> [Ivl; 3] {
    [p.x, p.y, p.z].map(Ivl::from)
}
/// Enclose coordinate-wise subtraction without discarding rounding error.
fn iv_sub(a: [Ivl; 3], b: [Ivl; 3]) -> [Ivl; 3] {
    let [ax, ay, az] = a;
    let [bx, by, bz] = b;
    [ax - bx, ay - by, az - bz]
}
/// Enclose the cross product using the original component operation order.
fn iv_cross(a: [Ivl; 3], b: [Ivl; 3]) -> [Ivl; 3] {
    let [ax, ay, az] = a;
    let [bx, by, bz] = b;
    [ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx]
}
/// Enclose the dot product with a fixed accumulation order.
fn iv_dot(a: [Ivl; 3], b: [Ivl; 3]) -> Ivl {
    a.into_iter()
        .zip(b)
        .fold(Ivl::from(0.0), |s, (x, y)| s + x * y)
}
/// Enclose Euclidean length, including intervals that cross zero.
fn iv_length(a: [Ivl; 3]) -> Ivl {
    a.into_iter()
        .fold(Ivl::from(0.0), |s, v| {
            s + Ivl::checked(v.mig() * v.mig(), v.mag() * v.mag())
        })
        .sqrt()
}
/// Enclose affine interpolation at an interval-valued fraction.
fn iv_lerp(a: [Ivl; 3], b: [Ivl; 3], t: Ivl) -> [Ivl; 3] {
    let [ax, ay, az] = a;
    let [bx, by, bz] = b;
    [ax + (bx - ax) * t, ay + (by - ay) * t, az + (bz - az) * t]
}

/// Tight directed bounds for affine definition proofs. Each IEEE binary64
/// basic operation is rounded once, then stepped one representable value outward;
/// signed zero and subnormal results are covered without a fixed absolute floor.
fn carrier_outward(lo: f64, hi: f64) -> Ivl {
    if !lo.is_finite() || !hi.is_finite() || lo > hi {
        Ivl::whole()
    } else {
        Ivl::new(next_parameter_down(lo), next_parameter_up(hi))
    }
}
/// Whether a carrier operand is an ordered, finite enclosure.
fn carrier_finite(a: Ivl) -> bool {
    a.lo.is_finite() && a.hi.is_finite() && a.lo <= a.hi
}
/// Enclose a sum, preserving exact addition of zero.
fn carrier_add(a: Ivl, b: Ivl) -> Ivl {
    if !carrier_finite(a) || !carrier_finite(b) {
        return Ivl::whole();
    }
    if b.lo == 0.0 && b.hi == 0.0 {
        a
    } else if a.lo == 0.0 && a.hi == 0.0 {
        b
    } else {
        carrier_outward(a.lo + b.lo, a.hi + b.hi)
    }
}
/// Enclose a difference, preserving exact subtraction of identical point inputs.
fn carrier_sub(a: Ivl, b: Ivl) -> Ivl {
    if !carrier_finite(a) || !carrier_finite(b) {
        return Ivl::whole();
    }
    if a.lo == a.hi && a == b {
        Ivl::from(0.0)
    } else if b.lo == 0.0 && b.hi == 0.0 {
        a
    } else {
        carrier_outward(a.lo - b.hi, a.hi - b.lo)
    }
}
/// Enclose a product; a nonfinite intermediate cannot establish a certificate.
fn carrier_mul(a: Ivl, b: Ivl) -> Ivl {
    if !carrier_finite(a) || !carrier_finite(b) {
        return Ivl::whole();
    }
    let products = [a.lo * b.lo, a.lo * b.hi, a.hi * b.lo, a.hi * b.hi];
    if products.iter().any(|value| !value.is_finite()) {
        return Ivl::whole();
    }
    if (a.lo == 0.0 && a.hi == 0.0) || (b.lo == 0.0 && b.hi == 0.0) {
        return Ivl::from(0.0);
    }
    if a.lo == 1.0 && a.hi == 1.0 {
        return b;
    }
    if b.lo == 1.0 && b.hi == 1.0 {
        return a;
    }
    carrier_outward(
        products.into_iter().fold(f64::INFINITY, f64::min),
        products.into_iter().fold(f64::NEG_INFINITY, f64::max),
    )
}
/// Enclose division directly at its four corners; a zero denominator refuses.
fn carrier_div(a: Ivl, b: Ivl) -> Ivl {
    if !carrier_finite(a) || !carrier_finite(b) {
        return Ivl::whole();
    }
    if !(b.lo > 0.0 || b.hi < 0.0) {
        return Ivl::whole();
    }
    if b.lo == 1.0 && b.hi == 1.0 {
        return a;
    }
    let values = [a.lo / b.lo, a.lo / b.hi, a.hi / b.lo, a.hi / b.hi];
    if values.iter().any(|value| !value.is_finite()) {
        return Ivl::whole();
    }
    carrier_outward(
        values.into_iter().fold(f64::INFINITY, f64::min),
        values.into_iter().fold(f64::NEG_INFINITY, f64::max),
    )
}
/// Dot product with outward rounding after each multiply and addition.
fn carrier_dot(a: [Ivl; 3], b: [Ivl; 3]) -> Ivl {
    a.into_iter().zip(b).fold(Ivl::from(0.0), |sum, (a, b)| {
        carrier_add(sum, carrier_mul(a, b))
    })
}
/// Cross product with outward rounding after each multiply and subtraction.
fn carrier_cross(a: [Ivl; 3], b: [Ivl; 3]) -> [Ivl; 3] {
    let [ax, ay, az] = a;
    let [bx, by, bz] = b;
    [
        carrier_sub(carrier_mul(ay, bz), carrier_mul(az, by)),
        carrier_sub(carrier_mul(az, bx), carrier_mul(ax, bz)),
        carrier_sub(carrier_mul(ax, by), carrier_mul(ay, bx)),
    ]
}
/// Nonnegative length bound; IEEE sqrt is correctly rounded before stepping out.
fn carrier_length(a: [Ivl; 3]) -> Ivl {
    let squared = a.into_iter().fold(Ivl::from(0.0), |sum, x| {
        let square = carrier_outward(x.mig() * x.mig(), x.mag() * x.mag());
        carrier_add(sum, square)
    });
    let bounds = carrier_outward(squared.lo.max(0.0).sqrt(), squared.hi.sqrt());
    Ivl::new(bounds.lo.max(0.0), bounds.hi)
}

// An exact zero predicate for four products of finite binary64 inputs. A
// product has at most 106 significand bits and exponent in [-2148, 1942].
// Four same-sign products need at most 4198 bits above 2^-2148; 66 u64 words
// supply 4224 bits. This fixed bound is independent of caller magnitudes.
/// Test a signed sum of four binary64 products exactly using bounded integer storage.
pub(crate) fn dyadic_four_products_zero(terms: [(f64, f64, bool); 4]) -> bool {
    dyadic_products_zero(&terms)
}

/// At most eight products need 4199 bits; the same 4224-bit storage suffices.
fn dyadic_products_zero(terms: &[(f64, f64, bool)]) -> bool {
    /// Split a finite binary64 value into its sign, integer significand and exponent.
    fn parts(value: f64) -> Option<(bool, u64, i32)> {
        if !value.is_finite() {
            return None;
        }
        let bits = value.to_bits();
        // The mask is at most 2047, so this cast is exact. Finite inputs exclude
        // 2047; the resulting power-of-two exponent is in [-1074, 971].
        let exponent = ((bits >> 52) & 0x7ff) as i32;
        let fraction = bits & ((1_u64 << 52) - 1);
        Some((
            bits >> 63 != 0,
            if exponent == 0 {
                fraction
            } else {
                fraction | (1_u64 << 52)
            },
            if exponent == 0 {
                -1074
            } else {
                exponent - 1075
            },
        ))
    }
    /// Accumulate one shifted word with checked access to the fixed-capacity sum.
    fn add_word(words: &mut [u64; 66], word: u64, shift: usize) -> bool {
        let mut index = shift / 64;
        // Constant nonzero divisor; remainder is 0..=63, below u128's 128 bits.
        let mut carry = u128::from(word) << (shift % 64);
        while carry != 0 {
            let Some(slot) = words.get_mut(index) else {
                return false;
            };
            // Two u64 values sum to at most 2^65-2. The cast deliberately keeps
            // the low word; sum >> 64 carries its single high bit forward.
            let sum = u128::from(*slot) + (carry & u128::from(u64::MAX));
            *slot = sum as u64;
            carry = (carry >> 64) + (sum >> 64);
            index += 1;
        }
        true
    }
    if terms.len() > 8 {
        return false;
    }
    let mut positive = [0_u64; 66];
    let mut negative = [0_u64; 66];
    for &(a, b, negate) in terms {
        let Some((sa, ma, ea)) = parts(a) else {
            return false;
        };
        let Some((sb, mb, eb)) = parts(b) else {
            return false;
        };
        // Each significand is below 2^53, so multiplication fits in 106 bits.
        // ea+eb+2148 is in [0,4090], fitting both i32 and wasm32 usize. The
        // second word's shift is at most 4154, within the 4224-bit accumulator.
        let product = u128::from(ma) * u128::from(mb);
        let shift = (ea + eb + 2148) as usize;
        let words = if sa ^ sb ^ negate {
            &mut negative
        } else {
            &mut positive
        };
        if !add_word(words, product as u64, shift)
            || !add_word(words, (product >> 64) as u64, shift + 64)
        {
            return false;
        }
    }
    positive == negative
}

/// Omit an affine boundary only when adjoining definitions are exactly smooth.
/// Cache G1 is an angular tolerance and cannot establish this identity.
fn exact_affine_join(curve: &Curve, t: f64) -> bool {
    if matches!(curve.kind(), CurveKind::Line(_)) {
        return true;
    }
    exact_native_spline_join(curve, t, true)
}

/// Establish smoothness from native knot multiplicity or exact endpoint chords.
/// At a C0 Bezier join, tangent directions use the adjacent defining controls;
/// equal positive local weights cancel without evaluating or rounding endpoints.
fn exact_native_spline_join(curve: &Curve, t: f64, require_linear: bool) -> bool {
    let CurveKind::BSpline(spline) = curve.kind() else {
        return false;
    };
    let Some(nurb) = spline.gw_nurb() else {
        return false;
    };
    let degree = nurb.degree();
    if degree == 0 || (require_linear && degree != 1) {
        return false;
    }
    let Some(order) = degree.checked_add(1) else {
        return false;
    };
    let natural = spline.natural_interval();
    if t > natural.min() && t < natural.max() {
        let multiplicity = nurb.knots().partition_point(|knot| *knot <= t)
            - nurb.knots().partition_point(|knot| *knot < t);
        // Interior non-knots and knots of multiplicity < degree have at least
        // C1 continuity by the defining spline space, not the cache angle test.
        if multiplicity < degree {
            let Some(first) = nurb.cpts().first() else {
                return false;
            };
            return nurb.cpts().iter().all(|point| {
                point.x.is_finite()
                    && point.y.is_finite()
                    && point.z.is_finite()
                    && if first.w == crate::nurbs::math::NOW {
                        point.w == crate::nurbs::math::NOW
                    } else {
                        point.w.is_finite() && point.w > 0.0
                    }
            });
        }
    } else if !nurb
        .knots()
        .iter()
        .take(order)
        .all(|knot| *knot == natural.min())
        || !nurb
            .knots()
            .iter()
            .rev()
            .take(order)
            .all(|knot| *knot == natural.max())
    {
        return false;
    }
    let left_t = if t == natural.min() { natural.max() } else { t };
    let right_t = if t == natural.max() { natural.min() } else { t };
    let span = |parameter: f64, left: bool| {
        let after = nurb.knots().partition_point(|knot| {
            if left {
                *knot < parameter
            } else {
                *knot <= parameter
            }
        });
        let i = after.checked_sub(1)?;
        let (&lo, &hi) = nurb.knots().get(i).zip(nurb.knots().get(i + 1))?;
        let contains = if left {
            lo < parameter && parameter <= hi
        } else {
            lo <= parameter && parameter < hi
        };
        (contains && hi > lo && i > 0).then_some(i)
    };
    let Some((li, ri)) = span(left_t, true).zip(span(right_t, false)) else {
        return false;
    };
    if li == ri && left_t == right_t {
        // An artificial subdivision inside the same affine span.
        return true;
    }
    let controls = |i: usize| {
        let (a, b) = nurb.cpts().get(i - 1).zip(nurb.cpts().get(i))?;
        let weight = if a.w == crate::nurbs::math::NOW && b.w == a.w {
            1.0
        } else if a.w == b.w && a.w > 0.0 && a.w.is_finite() {
            a.w
        } else {
            return None;
        };
        let coordinates = |p: &crate::nurbs::math::CPoint| {
            [
                p.x,
                p.y,
                if p.z == crate::nurbs::math::NOZ {
                    0.0
                } else {
                    p.z
                },
            ]
        };
        Some((coordinates(a), coordinates(b), weight))
    };
    let Some(right_control) = ri
        .checked_sub(degree)
        .and_then(|index| index.checked_add(1))
    else {
        return false;
    };
    let Some(((a, b, w), (c, d, v))) = controls(li).zip(controls(right_control)) else {
        return false;
    };
    // Homogeneous positions must agree even at repeated knots or a closed seam.
    if !b
        .into_iter()
        .zip(c)
        .all(|(b, c)| dyadic_products_zero(&[(b, v, false), (c, w, true)]))
    {
        return false;
    }
    let [ax, ay, az] = a;
    let [bx, by, bz] = b;
    let [cx, cy, cz] = c;
    let [dx, dy, dz] = d;
    let parallel = [
        (ax, ay, bx, by, cx, cy, dx, dy),
        (ay, az, by, bz, cy, cz, dy, dz),
        (az, ax, bz, bx, cz, cx, dz, dx),
    ]
    .into_iter()
    .all(|(ai, aj, bi, bj, ci, cj, di, dj)| {
        dyadic_products_zero(&[
            (bi, dj, false),
            (bi, cj, true),
            (ai, dj, true),
            (ai, cj, false),
            (bj, di, true),
            (bj, ci, false),
            (aj, di, false),
            (aj, ci, true),
        ])
    });
    // Parallel but reversed chords form a cusp, not a smooth join. Comparisons
    // determine difference signs exactly without subtracting or normalizing.
    parallel
        && a.into_iter()
            .zip(b)
            .zip(c.into_iter().zip(d))
            .any(|((a, b), (c, d))| (b > a && d > c) || (b < a && d < c))
}

/// Require defining smoothness for native spline and plane-UV joins.
/// A plane is an affine map of UV: it preserves an exact UV tangent identity.
/// In particular, native UV non-knots remain artificial subdivisions rather
/// than becoming extrema merely because the curve is wrapped by a surface.
fn omit_projected_boundary(curve: &Curve, t: f64, continuity: ContinuityType) -> bool {
    if (continuity as i32) < (ContinuityType::G1 as i32) {
        return false;
    }
    match curve.kind() {
        CurveKind::BSpline(_) => exact_native_spline_join(curve, t, false),
        CurveKind::CrvOnSurf(on_surface) => {
            if on_surface.surface().procedural_definition().is_some() {
                return false;
            }
            let crate::surface::surface::p2::SurfaceKind::Plane(plane) =
                on_surface.surface().kind()
            else {
                return false;
            };
            let uv = on_surface.uv_curve();
            let CurveKind::BSpline(spline) = uv.kind() else {
                return false;
            };
            let Some(net) = spline.gw_nurb() else {
                return false;
            };
            // Plane evaluation clamps outside its analytic domain. Prove the
            // defining positive-weight UV hull stays strictly inside, so that
            // clamping cannot introduce another, hidden nonsmooth join.
            let domain = plane.step_uv_domain();
            let inside = net.cpts().iter().all(|control| {
                let weight = if control.w == crate::nurbs::math::NOW {
                    1.0
                } else {
                    control.w
                };
                if !weight.is_finite() || weight <= 0.0 {
                    return false;
                }
                let x = Ivl::from(control.x) / Ivl::from(weight);
                let y = Ivl::from(control.y) / Ivl::from(weight);
                x.lo > domain.min().x
                    && x.hi < domain.max().x
                    && y.lo > domain.min().y
                    && y.hi < domain.max().y
            });
            inside && exact_native_spline_join(uv, t, false)
        }
        // Other curve kinds retain their existing boundary policy; this proof
        // does not establish smoothness of their wrapped definitions.
        _ => true,
    }
}

/// Exclude a continuous native C0 knot when a one-sided derivative proves that
/// the projected squared distance moves away from the requested extremum there.
fn native_join_is_not_extremum(
    curve: &Curve,
    t: f64,
    point: Vec3,
    normal: Vec3,
    maximize: bool,
) -> bool {
    let CurveKind::BSpline(spline) = curve.kind() else {
        return false;
    };
    let Some(nurb) = spline.gw_nurb() else {
        return false;
    };
    let degree = nurb.degree();
    let Some(order) = degree.checked_add(1) else {
        return false;
    };
    let first = nurb.knots().partition_point(|knot| *knot < t);
    let after = nurb.knots().partition_point(|knot| *knot <= t);
    let natural = spline.natural_interval();
    if degree == 0 || after - first < degree || t <= natural.min() || t >= natural.max() {
        return false;
    }
    let Some(left) = first.checked_sub(1) else {
        return false;
    };
    let Some(right) = after.checked_sub(order) else {
        return false;
    };
    let Some(a) = left.checked_sub(1).and_then(|index| nurb.cpts().get(index)) else {
        return false;
    };
    let Some(b) = nurb.cpts().get(left) else {
        return false;
    };
    let Some(c) = nurb.cpts().get(right) else {
        return false;
    };
    let Some(d) = right
        .checked_add(1)
        .and_then(|index| nurb.cpts().get(index))
    else {
        return false;
    };
    let data = |p: &crate::nurbs::math::CPoint| {
        let w = if p.w == crate::nurbs::math::NOW {
            1.0
        } else {
            p.w
        };
        let xyz = [
            p.x,
            p.y,
            if p.z == crate::nurbs::math::NOZ {
                0.0
            } else {
                p.z
            },
        ];
        (w.is_finite() && w > 0.0 && xyz.iter().all(|x| x.is_finite())).then_some((xyz, w))
    };
    let Some(((a, wa), (b, wb))) = data(a).zip(data(b)) else {
        return false;
    };
    let Some(((c, wc), (d, wd))) = data(c).zip(data(d)) else {
        return false;
    };
    if !b
        .into_iter()
        .zip(c)
        .all(|(b, c)| dyadic_products_zero(&[(b, wc, false), (c, wb, true)]))
    {
        return false;
    }
    let a = a.map(|x| Ivl::from(x) / Ivl::from(wa));
    let b = b.map(|x| Ivl::from(x) / Ivl::from(wb));
    let c = c.map(|x| Ivl::from(x) / Ivl::from(wc));
    let d = d.map(|x| Ivl::from(x) / Ivl::from(wd));
    let project = |v| iv_cross(iv_cross(iv_point(normal), v), iv_point(normal));
    let residual = project(iv_sub(b, iv_point(point)));
    let left = iv_dot(residual, project(iv_sub(b, a)));
    let right = iv_dot(residual, project(iv_sub(d, c)));
    if ![left.lo, left.hi, right.lo, right.hi]
        .into_iter()
        .all(f64::is_finite)
    {
        return false;
    }
    if maximize {
        left.hi < 0.0 || right.lo > 0.0
    } else {
        left.lo > 0.0 || right.hi < 0.0
    }
}

/// Test exact parallelism before rounding the difference of the two input points.
fn exact_parallel_difference(a: [f64; 3], b: [f64; 3], normal: Vec3) -> bool {
    let [ax, ay, az] = a;
    let [bx, by, bz] = b;
    // Expand N cross (B-A) without first rounding B-A. Constant positive
    // homogeneous weights and a line's scalar speed cancel from this identity.
    [
        (normal.x, normal.y, ax, ay, bx, by),
        (normal.y, normal.z, ay, az, by, bz),
        (normal.z, normal.x, az, ax, bz, bx),
    ]
    .into_iter()
    .all(|(ni, nj, ai, aj, bi, bj)| {
        dyadic_four_products_zero([
            (ni, bj, false),
            (ni, aj, true),
            (nj, bi, true),
            (nj, ai, false),
        ])
    })
}

/// Exact factors of an affine direction, before endpoint evaluation rounds them.
enum AffineDirection {
    Controls([f64; 3], [f64; 3]),
    MappedLine {
        vector: Vec3,
        scale: f64,
        parameter_ends: (f64, f64),
    },
}

pub(super) struct AffineDefinition {
    start: [Ivl; 3],
    end: [Ivl; 3],
    interval: Extent1d,
    direction: AffineDirection,
    isolated_end: Option<[Ivl; 3]>,
}
impl AffineDefinition {
    /// Enclose one continuous affine piece and its derivative without refitting it.
    pub(super) fn jet_over(&self, t: Ivl) -> ([Ivl; 3], [Ivl; 3]) {
        if self.interval.min() == self.interval.max() {
            // The restricted definition is a single point; its extension here is constant.
            return (self.start, [Ivl::from(0.0); 3]);
        }
        let span = Ivl::from(self.interval.max()) - Ivl::from(self.interval.min());
        let fraction = (t - Ivl::from(self.interval.min())) / span;
        let derivative = iv_sub(self.end, self.start).map(|x| x / span);
        (iv_lerp(self.start, self.end, fraction), derivative)
    }

    /// A directed enclosure from defining factors, independent of the wider
    /// arithmetic enclosure used to validate a reported native evaluation.
    /// At a discontinuity this evaluates this span's left limit.
    fn tight_point(&self, curve: &Curve, t: f64) -> Option<[Ivl; 3]> {
        match self.direction {
            AffineDirection::Controls(a, b) => {
                let CurveKind::BSpline(spline) = curve.kind() else {
                    return None;
                };
                let net = spline.gw_nurb()?;
                let index = net
                    .knots()
                    .partition_point(|knot| *knot <= self.interval.min())
                    .checked_sub(1)?;
                let weight = net.cpts().get(index)?.w;
                let weight = Ivl::from(if weight == crate::nurbs::math::NOW {
                    1.0
                } else {
                    weight
                });
                let fraction = carrier_div(
                    carrier_sub(Ivl::from(t), Ivl::from(self.interval.min())),
                    carrier_sub(
                        Ivl::from(self.interval.max()),
                        Ivl::from(self.interval.min()),
                    ),
                );
                let coordinate = |a, b| {
                    // Endpoints are defining controls, with no interpolation error.
                    let value = if t == self.interval.min() {
                        Ivl::from(a)
                    } else if t == self.interval.max() {
                        Ivl::from(b)
                    } else {
                        carrier_add(
                            Ivl::from(a),
                            carrier_mul(carrier_sub(Ivl::from(b), Ivl::from(a)), fraction),
                        )
                    };
                    carrier_div(value, weight)
                };
                let [ax, ay, az] = a;
                let [bx, by, bz] = b;
                Some([coordinate(ax, bx), coordinate(ay, by), coordinate(az, bz)])
            }
            AffineDirection::MappedLine { .. } => {
                let CurveKind::Line(line) = curve.kind() else {
                    return None;
                };
                let natural = line.natural_interval();
                let analytic = line.step_interval();
                let mut step = Ivl::from(t);
                if line.is_bounded() && !is_zero(natural.length()) {
                    let (a, b) = if line.inside_out() {
                        (natural.max(), natural.min())
                    } else {
                        (natural.min(), natural.max())
                    };
                    let fraction = carrier_div(
                        carrier_sub(step, Ivl::from(a)),
                        carrier_sub(Ivl::from(b), Ivl::from(a)),
                    );
                    step = carrier_add(
                        Ivl::from(analytic.min()),
                        carrier_mul(
                            fraction,
                            carrier_sub(Ivl::from(analytic.max()), Ivl::from(analytic.min())),
                        ),
                    );
                } else if !line.is_bounded() && line.inside_out() {
                    step = carrier_sub(Ivl::from(0.0), step);
                }
                let scale = carrier_mul(step, Ivl::from(line.line_scale()));
                let p = line.line_point();
                let v = line.line_vector();
                Some([
                    carrier_add(Ivl::from(p.x), carrier_mul(scale, Ivl::from(v.x))),
                    carrier_add(Ivl::from(p.y), carrier_mul(scale, Ivl::from(v.y))),
                    carrier_add(Ivl::from(p.z), carrier_mul(scale, Ivl::from(v.z))),
                ])
            }
        }
    }

    /// Enclose the mathematical endpoint difference from defining factors before
    /// rounded world-space endpoint evaluations can erase a short direction.
    fn difference(&self, curve: &Curve, interval: Extent1d) -> Option<[Ivl; 3]> {
        let width = Ivl::from(interval.max()) - Ivl::from(interval.min());
        let (direction, factor) = match self.direction {
            AffineDirection::Controls(a, b) => {
                let CurveKind::BSpline(spline) = curve.kind() else {
                    return None;
                };
                let nurb = spline.gw_nurb()?;
                let index = nurb
                    .knots()
                    .partition_point(|knot| *knot <= self.interval.min())
                    .checked_sub(1)?;
                let control = nurb.cpts().get(index)?;
                let weight = if control.w == crate::nurbs::math::NOW {
                    1.0
                } else {
                    control.w
                };
                (
                    iv_sub(b.map(Ivl::from), a.map(Ivl::from)),
                    width
                        / (Ivl::from(self.interval.max()) - Ivl::from(self.interval.min()))
                        / Ivl::from(weight),
                )
            }
            AffineDirection::MappedLine { vector, scale, .. } => {
                let CurveKind::Line(line) = curve.kind() else {
                    return None;
                };
                let natural = line.natural_interval();
                let mut factor = width * Ivl::from(scale);
                let bounded_map = line.is_bounded() && !is_zero(natural.length());
                if bounded_map {
                    let analytic = line.step_interval();
                    factor = factor * (Ivl::from(analytic.max()) - Ivl::from(analytic.min()))
                        / (Ivl::from(natural.max()) - Ivl::from(natural.min()));
                }
                if line.inside_out() && (bounded_map || !line.is_bounded()) {
                    factor = Ivl::from(0.0) - factor;
                }
                (iv_point(vector), factor)
            }
        };
        Some(direction.map(|component| component * factor))
    }

    /// Lower-bound distance to the infinite projected carrier without a rounded
    /// stationary parameter: |N|^2 |(A-Q).(N x D)| / |N x D|.
    /// The raw direction's nonzero scalar factors cancel. This bound remains
    /// valid on every bounded subinterval of the carrier.
    fn carrier_minimum(&self, curve: &Curve, point: Vec3, normal: Vec3) -> Option<Ivl> {
        let (anchor, direction) = match self.direction {
            AffineDirection::Controls(a, b) => {
                let CurveKind::BSpline(spline) = curve.kind() else {
                    return None;
                };
                let nurb = spline.gw_nurb()?;
                let index = nurb
                    .knots()
                    .partition_point(|knot| *knot <= self.interval.min())
                    .checked_sub(1)?;
                let weight = nurb.cpts().get(index)?.w;
                let weight = Ivl::from(if weight == crate::nurbs::math::NOW {
                    1.0
                } else {
                    weight
                });
                let [ax, ay, az] = a.map(Ivl::from);
                let [bx, by, bz] = b.map(Ivl::from);
                (
                    a.map(|x| carrier_div(Ivl::from(x), weight)),
                    [
                        carrier_sub(bx, ax),
                        carrier_sub(by, ay),
                        carrier_sub(bz, az),
                    ],
                )
            }
            AffineDirection::MappedLine { vector, .. } => {
                let CurveKind::Line(line) = curve.kind() else {
                    return None;
                };
                (iv_point(line.line_point()), iv_point(vector))
            }
        };
        let n = iv_point(normal);
        let transverse = carrier_cross(n, direction);
        let denominator = carrier_length(transverse);
        if denominator.lo <= 0.0 || !denominator.hi.is_finite() {
            return None;
        }
        let [ax, ay, az] = anchor;
        let offset = [
            carrier_sub(ax, Ivl::from(point.x)),
            carrier_sub(ay, Ivl::from(point.y)),
            carrier_sub(az, Ivl::from(point.z)),
        ];
        let numerator = carrier_dot(offset, transverse);
        let bound = carrier_div(
            carrier_mul(
                Ivl::new(numerator.mig(), numerator.mag()),
                carrier_dot(n, n),
            ),
            denominator,
        );
        (bound.lo.is_finite() && bound.hi.is_finite()).then_some(bound)
    }

    /// Enclose the affine point at the caller's raw parameter.
    fn point(&self, t: f64) -> [Ivl; 3] {
        let fraction = (Ivl::from(t) - Ivl::from(self.interval.min()))
            / (Ivl::from(self.interval.max()) - Ivl::from(self.interval.min()));
        iv_lerp(self.start, self.end, fraction)
    }
    /// Recognize an exactly constant projected objective without a tolerance test.
    fn constant_projection(&self, normal: Vec3) -> bool {
        match self.direction {
            AffineDirection::Controls(a, b) => exact_parallel_difference(a, b, normal),
            AffineDirection::MappedLine {
                vector,
                scale,
                parameter_ends: (a, b),
            } => {
                // The mapped endpoint difference factors as scale * (b-a) * V,
                // times a nonzero query/natural-span ratio for bounded maps.
                // Test its exact zero factors before multiplying or rounding
                // endpoint positions: neither underflow nor lost low bits may
                // turn a transverse direction into a constant projection.
                scale == 0.0
                    || a == b
                    || exact_parallel_difference([0.0; 3], [vector.x, vector.y, vector.z], normal)
            }
        }
    }
}

/// Enclose an analytic line or an equal-weight linear span, preserving isolated upper-knot values.
pub(super) fn affine_definition(curve: &Curve, interval: Extent1d) -> Option<AffineDefinition> {
    if let CurveKind::Line(line) = curve.kind() {
        let point = iv_point(line.line_point());
        let vector = iv_point(line.line_vector());
        let natural = line.natural_interval();
        let analytic = line.step_interval();
        let bounded_map = line.is_bounded() && !is_zero(natural.length());
        let endpoint = |t: f64| {
            let mut step = Ivl::from(t);
            if bounded_map {
                let fraction = if line.inside_out() {
                    (step - Ivl::from(natural.max()))
                        / (Ivl::from(natural.min()) - Ivl::from(natural.max()))
                } else {
                    (step - Ivl::from(natural.min()))
                        / (Ivl::from(natural.max()) - Ivl::from(natural.min()))
                };
                step = Ivl::from(analytic.min())
                    + fraction * (Ivl::from(analytic.max()) - Ivl::from(analytic.min()));
            } else if !line.is_bounded() && line.inside_out() {
                step = Ivl::from(0.0) - step;
            }
            let [x, y, z] = point;
            let [dx, dy, dz] = vector;
            let scale = step * Ivl::from(line.line_scale());
            [x + scale * dx, y + scale * dy, z + scale * dz]
        };
        let mapped_domain = if bounded_map { analytic } else { interval };
        return Some(AffineDefinition {
            start: endpoint(interval.min()),
            end: endpoint(interval.max()),
            interval,
            isolated_end: None,
            direction: AffineDirection::MappedLine {
                vector: line.line_vector(),
                scale: line.line_scale(),
                parameter_ends: (mapped_domain.min(), mapped_domain.max()),
            },
        });
    }
    let CurveKind::BSpline(spline) = curve.kind() else {
        return None;
    };
    let nurb = spline.gw_nurb()?;
    if nurb.degree() != 1 {
        return None;
    }
    // A right-sided lookup at the lower boundary avoids rounding a tiny span's
    // midpoint onto a knot. A leaf crossing a knot is not certified affine here.
    let index = nurb
        .knots()
        .partition_point(|t| *t <= interval.min())
        .saturating_sub(1);
    let (&lo, &hi) = nurb.knots().get(index).zip(nurb.knots().get(index + 1))?;
    if hi <= lo || interval.min() < lo || interval.max() > hi || index == 0 {
        return None;
    }
    let (a, b) = nurb.cpts().get(index - 1).zip(nurb.cpts().get(index))?;
    let nonrational = a.w == crate::nurbs::math::NOW && b.w == crate::nurbs::math::NOW;
    if !(nonrational || a.w == b.w && a.w > 0.0 && a.w.is_finite()) {
        return None;
    }
    let coords = |p: &crate::nurbs::math::CPoint| {
        let z = if p.z == crate::nurbs::math::NOZ {
            0.0
        } else {
            p.z
        };
        [p.x, p.y, z]
    };
    let ac = coords(a);
    let bc = coords(b);
    let w = Ivl::from(if nonrational { 1.0 } else { a.w });
    // Side::Left names the span's left boundary, so evaluate_point at an
    // interior full-multiplicity knot returns the NEXT span's first control.
    // The natural upper endpoint is special and retains the last left value.
    let isolated_end = if interval.max() == hi
        && hi < spline.natural_interval().max()
        && nurb.knots().get(index + 2) == Some(&hi)
    {
        let right = nurb.cpts().get(index + 1)?;
        let weight = if right.w == crate::nurbs::math::NOW {
            1.0
        } else if right.w > 0.0 && right.w.is_finite() {
            right.w
        } else {
            return None;
        };
        Some(coords(right).map(|x| Ivl::from(x) / Ivl::from(weight)))
    } else {
        None
    };
    Some(AffineDefinition {
        start: ac.map(|x| Ivl::from(x) / w),
        end: bc.map(|x| Ivl::from(x) / w),
        interval: Extent1d::from_min_max(lo, hi),
        direction: AffineDirection::Controls(ac, bc),
        isolated_end,
    })
}

/// Enclose the projected metric at an actual native point of an analytic line
/// or an equal-weight linear spline span. Native linear-spline knots and natural
/// endpoints use the evaluator's actual control value, including a right-valued
/// discontinuity. Unsupported curve definitions or parameters return `None`.
pub(crate) fn affine_projected_point_enclosure(
    curve: &Curve,
    t: f64,
    point: Vec3,
    normal: Vec3,
) -> Option<Ivl> {
    let natural = curve.natural_interval();
    if !t.is_finite() || t < natural.min() || t > natural.max() {
        return None;
    }
    let position = match curve.kind() {
        CurveKind::Line(_) => affine_definition(curve, Extent1d::from_min_max(t, t))?.start,
        CurveKind::BSpline(spline) => {
            let nurb = spline.gw_nurb()?;
            if nurb.degree() != 1 {
                return None;
            }
            let control = if t == natural.min() {
                nurb.cpts().first()
            } else if t == natural.max() {
                nurb.cpts().last()
            } else {
                let index = nurb
                    .knots()
                    .partition_point(|knot| *knot <= t)
                    .checked_sub(1)?;
                if nurb.knots().get(index) == Some(&t) {
                    index
                        .checked_sub(1)
                        .and_then(|index| nurb.cpts().get(index))
                } else {
                    None
                }
            };
            if let Some(control) = control {
                let weight = if control.w == crate::nurbs::math::NOW {
                    1.0
                } else {
                    control.w
                };
                if !weight.is_finite() || weight <= 0.0 {
                    return None;
                }
                [
                    control.x,
                    control.y,
                    if control.z == crate::nurbs::math::NOZ {
                        0.0
                    } else {
                        control.z
                    },
                ]
                .map(|x| Ivl::from(x) / Ivl::from(weight))
            } else {
                affine_definition(curve, Extent1d::from_min_max(t, t))?.point(t)
            }
        }
        _ => return None,
    };
    let bound = iv_length(iv_cross(
        iv_cross(iv_point(normal), iv_sub(position, iv_point(point))),
        iv_point(normal),
    ));
    (bound.lo.is_finite() && bound.hi.is_finite()).then_some(bound)
}

/// Certificates belong to the actual leaf candidates, not their scalar ordering.
/// The selected structural cohort traverses every defining knot span; omitted
/// smooth boundaries have their extrema owned by adjoining leaves. A proved
/// nonretained minimum leaf lies above a feasible incumbent's upper enclosure,
/// so cannot lower the global optimum below the recorded minimum bound.
#[derive(Debug, Clone)]
pub(super) struct ProjectedAffineSelection<O> {
    candidates: Vec<(Solution<O>, Ivl)>,
    maximum_upper: f64,
    minimum_lower: f64,
    span_lower: f64,
    pub(super) pending_spans: Option<Vec<Extent1d>>,
}

/// Select before solving: every active native piece must have the existing
/// exact affine definition. Mixed rational/non-affine definitions retain their
/// established construction; they are never omitted from a selected proof.
fn complete_affine_definition(curve: &Curve, interval: Extent1d) -> Option<Vec<Extent1d>> {
    if !interval.min().is_finite()
        || !interval.max().is_finite()
        || interval.min() >= interval.max()
    {
        return None;
    }
    let mut spans = Vec::new();
    match curve.kind() {
        CurveKind::Line(_) => {
            affine_definition(curve, interval)?;
            spans.push(interval);
        }
        CurveKind::BSpline(spline) => {
            let net = spline.gw_nurb()?;
            if net.degree() != 1 {
                return None;
            }
            for pair in net.knots().windows(2) {
                let &[a, b] = pair else {
                    return None;
                };
                let a = a.max(interval.min());
                let b = b.min(interval.max());
                if a < b {
                    let span = Extent1d::from_min_max(a, b);
                    affine_definition(curve, span)?;
                    spans.push(span);
                }
            }
        }
        _ => return None,
    }
    (spans.first()?.min() == interval.min() && spans.last()?.max() == interval.max())
        .then_some(spans)
}

impl<O: SolutionObject> ProjectedAffineSelection<O> {
    fn finish(
        &self,
        curve: &Curve,
        point: Vec3,
        normal: Vec3,
        maximize: bool,
        requested: SolutionRequested,
        target: Option<f64>,
        retention: f64,
        solutions: &mut SolutionArray<O>,
        mut identical: impl FnMut(&Solution<O>, &Solution<O>) -> Result<bool>,
    ) -> Result<()> {
        if self.pending_spans.is_some() {
            // An empty cache traversal is not a proof about the definition.
            return Err(affine_refusal());
        }
        let limit = if maximize {
            self.maximum_upper
        } else {
            self.minimum_lower
        };
        if !limit.is_finite() {
            return Err(affine_refusal());
        }
        let cutoff = target.map(|target| {
            if maximize {
                carrier_sub(Ivl::from(target.abs()), Ivl::from(retention))
            } else {
                carrier_add(Ivl::from(target), Ivl::from(retention))
            }
        });
        let Some(current) = solutions.as_slice().first() else {
            // Scalar filtering can discard every rounded report even though
            // the actual optimum qualifies. Empty success requires a proof of
            // exclusion; overlap with the external band explicitly refuses.
            return if cutoff.is_some_and(|cutoff| {
                carrier_finite(cutoff)
                    && if maximize {
                        limit < cutoff.lo
                    } else {
                        limit > cutoff.hi
                    }
            }) {
                Ok(())
            } else {
                Err(affine_refusal())
            };
        };
        let certifies = |solution: &Solution<O>| {
            let value = solution.start.solution_value;
            let certifies_bound = |own: Ivl| {
                value >= own.lo
                    && value <= own.hi
                    && if maximize {
                        own.hi >= limit
                    } else {
                        own.lo <= limit
                    }
            };
            // Reuse only this exact returned native point's original proof.
            // A different algebraic evaluation can be tighter without making
            // the original own-point enclosure invalid. Never borrow another
            // parameter's width, a left limit, or a whole-range hull.
            if self.candidates.iter().any(|(candidate, own)| {
                candidate.start.parameters[0] == solution.start.parameters[0]
                    && candidate.start.solution_value == value
                    && certifies_bound(*own)
            }) {
                return true;
            }
            let Some(own) = affine_projected_point_enclosure(
                curve,
                solution.start.parameters[0],
                point,
                normal,
            ) else {
                return false;
            };
            certifies_bound(own)
        };
        if certifies(current) {
            if requested == SolutionRequested::All {
                let feasible = self
                    .candidates
                    .iter()
                    .map(|(_, bound)| if maximize { bound.lo } else { bound.hi })
                    .fold(
                        if maximize {
                            f64::NEG_INFINITY
                        } else {
                            f64::INFINITY
                        },
                        |a, b| if maximize { a.max(b) } else { a.min(b) },
                    );
                let band = if maximize {
                    carrier_sub(Ivl::from(feasible), Ivl::from(retention))
                } else {
                    carrier_add(Ivl::from(feasible), Ivl::from(retention))
                };
                if !carrier_finite(band) {
                    return Err(affine_refusal());
                }
                for (candidate, bound) in &self.candidates {
                    let outside = |cutoff: Ivl| {
                        carrier_finite(cutoff)
                            && if maximize {
                                bound.hi < cutoff.lo
                            } else {
                                bound.lo > cutoff.hi
                            }
                    };
                    if outside(band) || cutoff.is_some_and(outside) {
                        continue;
                    }
                    let mut represented = false;
                    for retained in solutions.as_slice() {
                        if identical(candidate, retained)?
                            || affine_range_represents(retained, candidate)
                        {
                            represented = true;
                            break;
                        }
                    }
                    if !represented {
                        return Err(affine_refusal());
                    }
                }
            }
            return Ok(());
        }
        // All's existing canonical retention and range merging is preserved.
        // An ambiguous reported ordering cannot silently become a certified set.
        if requested != SolutionRequested::Single {
            return Err(affine_refusal());
        }
        let Some((replacement, _)) = self
            .candidates
            .iter()
            .find(|(solution, _)| certifies(solution))
        else {
            return Err(affine_refusal());
        };
        if let Some(cutoff) = cutoff {
            // Replacing a rounded winner must not resurrect a report excluded
            // by the caller's external incumbent. Overlapping cutoff arithmetic
            // refuses instead of silently changing retention semantics.
            let value = replacement.start.solution_value;
            if !carrier_finite(cutoff)
                || if maximize {
                    value < cutoff.hi
                } else {
                    value > cutoff.lo
                }
            {
                return Err(affine_refusal());
            }
        }
        *solutions = SolutionArray::default();
        solutions.add(replacement);
        Ok(())
    }
}

/// The canonical accumulator absorbs points and subranges into constant ranges.
/// This exact containment check introduces no new geometric dedup tolerance.
fn affine_range_represents<O: SolutionObject>(
    range: &Solution<O>,
    candidate: &Solution<O>,
) -> bool {
    if range.solution_type != SolutionType::RangeOfValues
        || range.num_variables != 1
        || candidate.num_variables != 1
        || range.num_objects != candidate.num_objects
        || range.objects != candidate.objects
    {
        return false;
    }
    let a = range.start.parameters[0].min(range.end.parameters[0]);
    let b = range.start.parameters[0].max(range.end.parameters[0]);
    let t = candidate.start.parameters[0];
    a <= t
        && t <= b
        && (candidate.solution_type != SolutionType::RangeOfValues
            || (a <= candidate.end.parameters[0] && candidate.end.parameters[0] <= b))
}

impl<'t, O: SolutionObject + 't> CpGlobalSolver<'t, O> {
    fn record_affine_candidate(&mut self, solution: &Solution<O>, bound: Ivl) {
        if let Some(selection) = &mut self.projected_affine_selection {
            selection.candidates.push((*solution, bound));
            selection.minimum_lower = selection
                .minimum_lower
                .min(bound.lo.max(selection.span_lower));
        }
    }
    /// `AddSortedSolution` by value (with parameter) or by first parameter.
    fn add_solution(&mut self, solution: &Solution<O>, by_value: bool) -> Result<()> {
        let key = if by_value {
            SortKey::ByValueWithParameter
        } else {
            SortKey::ByFirstParameter
        };
        self.add_sorted_solution(solution, key, true)
    }

    /// Exact affine geometry dispatch; `false` means this leaf is not affine.
    fn projected_affine(
        &mut self,
        k: &Kernel,
        curve: &Curve,
        interval: Extent1d,
        conts: [ContinuityType; 2],
        point: Vec3,
        branch: &Branch,
        owner: Option<O>,
    ) -> Result<bool> {
        let Some(definition) = affine_definition(curve, interval) else {
            if self.projected_affine_selection.is_some() {
                return Err(affine_refusal());
            }
            return Ok(false);
        };
        let normal = self.opt(0)?;
        let constant = definition.constant_projection(normal);
        if let Some(selection) = &mut self.projected_affine_selection {
            // A norm of an affine vector is convex: both defining endpoints
            // bound the whole span, even when an endpoint is not retained.
            // An isolated right-valued knot is an additional native point.
            let project = |position: [Ivl; 3]| {
                let [x, y, z] = position;
                let offset = [
                    carrier_sub(x, Ivl::from(point.x)),
                    carrier_sub(y, Ivl::from(point.y)),
                    carrier_sub(z, Ivl::from(point.z)),
                ];
                carrier_cross(carrier_cross(iv_point(normal), offset), iv_point(normal))
            };
            let mut hull: Option<[Ivl; 3]> = None;
            let mut include = |position| -> Result<()> {
                let projected = project(position);
                let bound = carrier_length(projected);
                if !carrier_finite(bound) {
                    return Err(affine_refusal());
                }
                selection.maximum_upper = selection.maximum_upper.max(bound.hi);
                hull = Some(match hull {
                    None => projected,
                    Some(previous) => {
                        let union = |a: Ivl, b: Ivl| Ivl::new(a.lo.min(b.lo), a.hi.max(b.hi));
                        let [ax, ay, az] = previous;
                        let [bx, by, bz] = projected;
                        [union(ax, bx), union(ay, by), union(az, bz)]
                    }
                });
                Ok(())
            };
            // Exact parallelism proves that the projected vector is constant.
            // Re-evaluating a remote endpoint would lose that shared-factor
            // identity to interval dependency without adding any coverage.
            let upper = if constant {
                interval.min()
            } else {
                interval.max()
            };
            for t in [interval.min(), upper] {
                include(
                    definition
                        .tight_point(curve, t)
                        .ok_or_else(affine_refusal)?,
                )?;
            }
            if definition.isolated_end.is_some() {
                let t = interval.max();
                // Only the right knot point belongs to this query; the next
                // span need not itself have equal weights or be affine.
                let CurveKind::BSpline(spline) = curve.kind() else {
                    return Err(affine_refusal());
                };
                let net = spline.gw_nurb().ok_or_else(affine_refusal)?;
                let index = net
                    .knots()
                    .partition_point(|knot| *knot <= t)
                    .checked_sub(2)
                    .ok_or_else(affine_refusal)?;
                let control = net.cpts().get(index).ok_or_else(affine_refusal)?;
                let weight = Ivl::from(if control.w == crate::nurbs::math::NOW {
                    1.0
                } else {
                    control.w
                });
                let z = if control.z == crate::nurbs::math::NOZ {
                    0.0
                } else {
                    control.z
                };
                include(
                    [control.x, control.y, z].map(|value| carrier_div(Ivl::from(value), weight)),
                )?;
            }
            // The affine vector lies inside its endpoint hull. Its norm's
            // lower bound can tighten the local stationary-point proof,
            // particularly when the objective is almost constant.
            selection.span_lower = carrier_length(hull.ok_or_else(affine_refusal)?).lo;
        }
        let projected = |p| {
            iv_cross(
                iv_cross(iv_point(normal), iv_sub(p, iv_point(point))),
                iv_point(normal),
            )
        };
        let start = projected(definition.point(interval.min()));
        let end = projected(definition.point(interval.max()));
        let Some(difference) = definition.difference(curve, interval) else {
            return Err(affine_refusal());
        };
        let delta = iv_cross(iv_cross(iv_point(normal), difference), iv_point(normal));
        let left_end = if definition.isolated_end.is_some() {
            next_parameter_down(interval.max())
        } else {
            interval.max()
        };
        let maximize = self.base.operation == SolverOperation::ProjectedMaximize;
        let operation = self.base.operation;
        let tolerance = self.base.tolerance_3d;
        if !tolerance.is_finite() || tolerance < 0.0 {
            return Err(affine_refusal());
        }
        let make_solution = |t: f64, value: f64| {
            let mut solution = Solution::new();
            solution.solution_type = SolutionType::SingleValue;
            solution.nodes[0] = branch[0];
            solution.nodes[1] = branch[1];
            solution.objects[0] = owner;
            solution.num_objects = 1;
            solution.num_variables = 1;
            solution.start.parameters[0] = t;
            solution.start.solution_value = value;
            solution
        };
        let metric = |t: f64| -> Result<(f64, Ivl)> {
            let value = proj_dist(k, curve, t, point, normal, true, Vec3::ZERO, operation)?;
            let position = if t == interval.max() {
                definition
                    .isolated_end
                    .unwrap_or_else(|| definition.point(t))
            } else {
                definition.point(t)
            };
            let bounds = iv_length(projected(position));
            if !value.is_finite()
                || !bounds.lo.is_finite()
                || !bounds.hi.is_finite()
                || value < bounds.lo
                || value > bounds.hi
            {
                return Err(affine_refusal());
            }
            // The enclosure certifies optimization error, not exact binary
            // representation of a generally irrational reported distance.
            Ok((value, bounds))
        };
        if constant {
            // The exact parallel identity proves zero optimization error for
            // every parameter, including a zero requested distance tolerance.
            let (a, a_bound) = metric(interval.min())?;
            let (b, _) = metric(left_end)?;
            if definition.isolated_end.is_some() {
                let (value, bound) = metric(interval.max())?;
                let isolated = make_solution(interval.max(), value);
                self.record_affine_candidate(&isolated, bound);
                self.add_solution(&isolated, true)?;
            }
            let mut solution = make_solution(interval.min(), a);
            if interval.min() < left_end {
                solution.solution_type = SolutionType::RangeOfValues;
            }
            solution.end = solution.start;
            solution.end.parameters[0] = left_end;
            solution.end.solution_value = b;
            self.record_affine_candidate(&solution, a_bound);
            self.add_solution(&solution, true)?;
            return Ok(true);
        }
        let length = iv_length(delta);
        if !length.hi.is_finite() {
            return Err(affine_refusal());
        }
        let foot = if length.lo <= 0.0 {
            // Division by an enclosure containing zero cannot locate the foot.
            // Enclose the whole span instead; the candidate-gap check below
            // still has to certify the requested accuracy, including a distinct
            // upper-knot value. This is not an exact constant-projection range.
            Ivl::new(0.0, 1.0)
        } else {
            let length_sq = length * length;
            (Ivl::from(0.0) - iv_dot(start, delta)) / length_sq
        };
        if !foot.lo.is_finite() || !foot.hi.is_finite() {
            return Err(affine_refusal());
        }
        let fraction = Ivl::new(foot.lo.clamp(0.0, 1.0), foot.hi.clamp(0.0, 1.0));
        let smooth = [
            (conts[0] as i32) >= (ContinuityType::G1 as i32)
                && (exact_affine_join(curve, interval.min())
                    || native_join_is_not_extremum(curve, interval.min(), point, normal, maximize)),
            (conts[1] as i32) >= (ContinuityType::G1 as i32)
                && (exact_affine_join(curve, interval.max())
                    || native_join_is_not_extremum(curve, interval.max(), point, normal, maximize)),
        ];
        if !maximize
            && definition.isolated_end.is_none()
            && ((foot.hi <= 0.0 && smooth[0]) || (foot.lo >= 1.0 && smooth[1]))
        {
            return Ok(true);
        }
        let mut optimum = if foot.hi <= 0.0 {
            iv_length(start)
        } else if foot.lo >= 1.0 {
            iv_length(end)
        } else {
            iv_length(iv_lerp(start, end, fraction))
        };
        if let Some(carrier) = definition.carrier_minimum(curve, point, normal) {
            optimum.lo = optimum.lo.max(carrier.lo);
        }
        if let Some(position) = definition.isolated_end {
            let isolated = iv_length(projected(position));
            optimum.lo = optimum.lo.min(isolated.lo);
            optimum.hi = optimum.hi.min(isolated.hi);
        }
        if let Some(selection) = &self.projected_affine_selection {
            optimum.lo = optimum.lo.max(selection.span_lower);
        }
        if !maximize {
            // A proved nonretained leaf needs no representable local minimizer.
            // Re-enclose the incumbent from its actual defining span: no report
            // value or another evaluation's width is used as solve accuracy.
            if let Some(incumbent) = self.base.solutions.as_slice().first() {
                let t = incumbent.start.parameters[0];
                if let Some(bound) = affine_projected_point_enclosure(curve, t, point, normal) {
                    let cutoff = Ivl::from(bound.hi) + Ivl::from(tolerance);
                    if cutoff.hi.is_finite() && optimum.lo > cutoff.hi {
                        return Ok(true);
                    }
                }
            }
        }
        let mut parameters = Vec::new();
        for (t, omit) in [(interval.min(), smooth[0]), (interval.max(), smooth[1])] {
            if !omit {
                parameters.push(t);
            }
        }
        if definition.isolated_end.is_some() && !parameters.contains(&left_end) {
            parameters.push(left_end);
        }
        if !maximize && fraction.mid() > 0.0 && fraction.mid() < 1.0 {
            let t = interval.evaluate(fraction.mid());
            for candidate in [next_parameter_down(t), t, next_parameter_up(t)] {
                if candidate > interval.min() && candidate < interval.max() {
                    parameters.push(candidate);
                }
            }
        }
        // Smooth joins belong to their adjoining stationary search. They need
        // not become artificial endpoint answers of this leaf.
        if parameters.is_empty() {
            // A smooth affine maximum belongs to adjoining leaves. An interior
            // minimum, however, cannot disappear merely because no binary64
            // parameter lies inside this span.
            return if maximize {
                Ok(true)
            } else {
                Err(affine_refusal())
            };
        }
        let mut candidates = Vec::new();
        for t in parameters {
            let (value, bound) = metric(t)?;
            candidates.push((t, value, bound));
        }
        candidates.sort_by(|a, b| {
            if maximize {
                b.2.hi.total_cmp(&a.2.hi).then_with(|| b.1.total_cmp(&a.1))
            } else {
                a.2.lo.total_cmp(&b.2.lo).then_with(|| a.1.total_cmp(&b.1))
            }
        });
        let Some(&(best_t, best_value, best_bound)) = candidates.first() else {
            return Ok(true);
        };
        let proven_optimal = if definition.isolated_end.is_some() {
            false
        } else if maximize {
            candidates
                .iter()
                .all(|(t, _, bound)| *t == best_t || bound.hi <= best_bound.lo)
        } else {
            (best_t == interval.min() && foot.hi < 0.0)
                || (best_t == interval.max() && foot.lo > 1.0)
        };
        if !proven_optimal {
            let certified_to_own_enclosure = if maximize {
                // Only eligible boundaries belong to this leaf's maximum search.
                let mut upper = candidates
                    .iter()
                    .map(|(_, _, b)| b.hi)
                    .fold(f64::NEG_INFINITY, f64::max);
                if definition.isolated_end.is_some() {
                    upper = upper.max(iv_length(end).hi);
                }
                upper <= best_bound.hi
            } else {
                optimum.lo >= best_bound.lo
            };
            // Retention slack only filters solutions in add_sorted_solution.
            // Neither that slack nor a different evaluation's interval width
            // may make this representable candidate an accurate construction.
            if !certified_to_own_enclosure {
                return Err(affine_refusal());
            }
        }
        let best = make_solution(best_t, best_value);
        self.record_affine_candidate(&best, best_bound);
        self.add_solution(&best, true)?;
        if self.base.solution_requested == SolutionRequested::Single {
            return Ok(true);
        }
        // Preserve eligible endpoint candidates for All mode, but do not emit
        // adjacent floating-point representatives as distinct stationary roots.
        for (t, value, bound) in candidates.into_iter().skip(1) {
            if t == interval.min() || t == interval.max() || t == left_end {
                let candidate = make_solution(t, value);
                self.record_affine_candidate(&candidate, bound);
                self.add_solution(&candidate, true)?;
            }
        }
        Ok(true)
    }

    /// Option vector `i` (`m_cpOptVectors[i]`, read without a `NULL` check).
    fn opt(&self, i: usize) -> Result<Vec3> {
        self.base
            .opt_vectors
            .get(i)
            .copied()
            .ok_or(SolidError::invalid_input(OP_LOCAL, "missing option vector"))
    }

    /// `SmCPGlobalSolver::LocalSolve` on the point node `branch[0]` and curve leaf
    /// `branch[1]` (always returns `rbNeedsMoreSubdivision` false): the span interval
    /// is limited to the solver's (span ends at a closed curve's seam take the seam
    /// continuity), spans out of reach of the best answer or tolerance are culled,
    /// a degenerate span is a point solution, then each guess (the closest point of
    /// the span chord, its middle and ends) is solved locally (by
    /// `LocalPropertyAnalysis` for the directed, projected and pivot operations, else
    /// [`Curve::local_point_solve`]); boundary solutions that continue into the
    /// neighbouring span are skipped, solutions equal to their guess are collected
    /// and three of them make the span a range solution.
    ///
    /// # Errors
    ///
    /// `SM_ERR_NULL_POINTER` for a node without data; `InvalidInput` for a missing
    /// branch node, curve or option vector or a node that is not a curve leaf; those
    /// of `Unitize`, the signed angle, `AddSortedSolution` and `LocalPropertyAnalysis`.
    pub fn local_solve_branch(&mut self, k: &Kernel, branch: &Branch) -> Result<bool> {
        // The affine construction follows defining knot pieces directly. A cache
        // may omit geometrically degenerate pieces, whose projected distance can
        // nevertheless be the global extremum. Cache topology is not coverage.
        if let Some(spans) = self
            .projected_affine_selection
            .as_mut()
            .and_then(|selection| selection.pending_spans.take())
        {
            let curve = self.curve.ok_or_else(affine_refusal)?;
            let point = self.base.branch_node(branch, 0)?.bbox.min();
            let tree = self.base.tree(1)?;
            let owner = tree.owner_object()?;
            let mut definition_branch = *branch;
            definition_branch[1] = tree.top_node();
            for span in spans {
                k.check_cancelled(OP_GLOBAL)?;
                let conts = [
                    if span.min() == self.interval.min() {
                        ContinuityType::Discontinuous
                    } else {
                        ContinuityType::G1
                    },
                    if span.max() == self.interval.max() {
                        ContinuityType::Discontinuous
                    } else {
                        ContinuityType::G1
                    },
                ];
                self.projected_affine(k, curve, span, conts, point, &definition_branch, owner)?;
            }
            return Ok(false);
        }
        // init outputs
        let needs_more_subdivision = false;
        // the fork passes on the statuses below (tangency-labs/smlib 8311f8e)
        let pass_on = k.fixes.query_statuses_pass_on;

        // locals
        let point_node = self.base.branch_node(branch, 0)?;
        let curve_node = self.base.branch_node(branch, 1)?;
        if matches!(curve_node.data, NodeData::None) {
            return Err(SolidError::null_pointer(OP_LOCAL));
        }
        let curve_tree = self.base.tree(1)?;
        let bez_span = *curve_tree.curve_span(curve_node)?;
        let mut bez_ivl = *bez_span.interval();
        let test_pnt = point_node.bbox.min();
        let curve = self
            .curve
            .ok_or(SolidError::invalid_input(OP_LOCAL, "solver without curve"))?;
        let curve_obj = curve_tree.owner_object()?;
        let op = self.base.operation;
        let requested = self.base.solution_requested;
        let tol_3d = self.base.tolerance_3d;
        let opt_vectors = self.base.opt_vectors.clone();
        let interval = self.interval;

        // Do interval tests
        let mut conts = bez_span.aux.conts;

        // no work - disjoint leafNode and problem intervals
        if bez_ivl.are_disjoint(&interval, 0.0) {
            return Ok(needs_more_subdivision);
        }

        // when BezierInterval extends beyond TotalProblemInterval - mark this range
        // ends as discontinuous
        if bez_ivl.min() <= interval.min() {
            conts[0] = ContinuityType::Discontinuous;
        }
        if bez_ivl.max() >= interval.max() {
            conts[1] = ContinuityType::Discontinuous;
        }

        // trim sBezIvl to m_crInterval
        let projected_extremum = matches!(
            op,
            SolverOperation::ProjectedMinimize | SolverOperation::ProjectedMaximize
        );
        let containment_tolerance = if projected_extremum { 0.0 } else { EFF_ZERO };
        if !bez_ivl.is_contained_by(&interval, containment_tolerance) {
            bez_ivl = bez_ivl.intersect(&interval)?;
            // no work - zero length trimmed interval
            if if projected_extremum {
                bez_ivl.max() == bez_ivl.min()
            } else {
                are_same(bez_ivl.max(), bez_ivl.min())
            } {
                return Ok(needs_more_subdivision);
            }
        }

        // adjust continuity for boundary endSpans on closed curves
        let parent_ivl = curve.natural_interval();
        if (is_zero(bez_ivl.min() - parent_ivl.min()) || is_zero(bez_ivl.max() - parent_ivl.max()))
            && (is_zero(interval.min() - parent_ivl.min())
                && is_zero(interval.max() - parent_ivl.max()))
            && curve.is_closed(k, &parent_ivl, 0.0)?
        {
            // get pCurve1's end continuity (status ignored; passed on with the opt-in
            // fix, tangency-labs/smlib 8311f8e)
            let mut continuity = ContinuityType::Discontinuous;
            if let Some(c) = ignored_or_pass_on(
                evaluate_continuity(
                    k,
                    curve,
                    parent_ivl.min(),
                    false,
                    curve,
                    parent_ivl.max(),
                    true,
                    CONTINUITY_ANGLE,
                ),
                pass_on,
            )? {
                continuity = c.continuity;
            }
            // set span endContinuities as appropriate
            if is_zero(bez_ivl.min() - parent_ivl.min()) {
                conts[0] = continuity;
            }
            if is_zero(bez_ivl.max() - parent_ivl.max()) {
                conts[1] = continuity;
            }
        }

        if matches!(
            op,
            SolverOperation::ProjectedMinimize | SolverOperation::ProjectedMaximize
        ) && self.projected_affine(k, curve, bez_ivl, conts, test_pnt, branch, curve_obj)?
        {
            return Ok(needs_more_subdivision);
        }

        // cull no hope minimization and SingleSolution normal situations
        if op == SolverOperation::Minimize
            || (op == SolverOperation::Normalize && requested == SolutionRequested::Single)
        {
            // quit when point/pseudoBox distance is greater than tolerance to current
            // best answer
            let mut dist = bez_span.pseudo_box().distance_to_point(test_pnt);
            dist -= tol_3d;
            if dist < 0.0 {
                dist = 0.0;
            }
            if dist * dist > self.base.best_answer_so_far_sq {
                return Ok(needs_more_subdivision);
            }
        }

        // cull no hope intersection cases
        if op == SolverOperation::Intersect {
            let dist = bez_span.pseudo_box().distance_to_point(test_pnt);
            if dist > tol_3d {
                return Ok(needs_more_subdivision);
            }
        }

        let max_turn = bez_span.aux.max_turning_angle_deg;
        let do_line_test = max_turn < 90.0
            && (op == SolverOperation::Normalize
                || ((op == SolverOperation::Minimize || op == SolverOperation::Maximize)
                    && ((conts[0] as i32) > (ContinuityType::C0 as i32)
                        && (conts[1] as i32) > (ContinuityType::C0 as i32))));

        // locals
        let mut guesses: Vec<f64> = Vec::with_capacity(6);

        // get BezierCurve EndPoints and the vector between them (statuses ignored)
        let mut start_pt = [Vec3::UNDEF; 2];
        let mut end_pt = [Vec3::UNDEF; 2];
        if let Some(v) =
            ignored_or_pass_on(curve.evaluate(k, bez_ivl.min(), 1, true, true), pass_on)?
        {
            start_pt = first_n(&v);
        }
        if let Some(v) =
            ignored_or_pass_on(curve.evaluate(k, bez_ivl.max(), 1, true, true), pass_on)?
        {
            end_pt = first_n(&v);
        }

        let vec = end_pt[0] - start_pt[0];

        // negate sStartPt tangent so that both end Point tangent vectors point outside
        start_pt[1] = -start_pt[1];

        let by_value_degenerate = op == SolverOperation::Minimize
            || op == SolverOperation::Maximize
            || (op == SolverOperation::Normalize && requested == SolutionRequested::Single);

        // when bezierCurve segment is degenerate (zeroLength)
        if vec.length_squared() < EFF_ZERO_SQ {
            // when looking for min/max/intersect solutions
            if op == SolverOperation::Minimize
                || op == SolverOperation::Maximize
                || op == SolverOperation::Intersect
                || (op == SolverOperation::Normalize && requested == SolutionRequested::Single)
            {
                // get point/point distance
                let dist = end_pt[0].distance(test_pnt);

                // cull out of tolerance intersect cases
                if op == SolverOperation::Intersect && dist > tol_3d {
                    return Ok(needs_more_subdivision);
                }

                // save this as a solution
                let mut solution: Solution<O> = Solution::new();
                solution.solution_type = SolutionType::SingleValue;
                solution.nodes[0] = branch[0];
                solution.nodes[1] = branch[1];
                solution.num_variables = 1;
                solution.start.solution_value = dist;
                solution.start.parameters[0] = bez_ivl.min();
                solution.num_objects = 1;
                solution.objects[0] = curve_obj;

                // add solution to sSolution array in proper order
                self.add_solution(&solution, by_value_degenerate)?;
            }

            return Ok(needs_more_subdivision);
        }

        // arrive here when bezier segment is not zeroLength
        // get closest point on endPoint lineSegment to testPoint
        let (ang_deg_from_perp, t) =
            segment_point_angle_deg_to_perp(start_pt[0], end_pt[0], test_pnt);

        // when doing the line test
        if do_line_test
            && ang_deg_from_perp * 0.99 - EFF_ZERO > max_turn
            && !test_pnt.closer_than(tol_3d, start_pt[0])
            && !test_pnt.closer_than(tol_3d, end_pt[0])
        {
            return Ok(needs_more_subdivision); // no Normal point possible on this segment
        }

        // for all cases but SM_SO_AT_DISTANCE
        if op == SolverOperation::AtDistance {
            // add end and mid Point locations as guesses
            guesses.push(bez_ivl.evaluate(0.5));
            guesses.push(bez_ivl.evaluate(0.0));
            guesses.push(bez_ivl.evaluate(1.0));
        } else {
            // add closest linePoint to guess array and unique end and mid point positions
            guesses.push(bez_ivl.evaluate(t));
            if !is_zero_to_tol(t - 0.5, EFF_ZERO_SQRT) {
                guesses.push(bez_ivl.evaluate(0.5));
            }
            if !is_zero_to_tol(t - 0.0, EFF_ZERO_SQRT) {
                guesses.push(bez_ivl.evaluate(0.0));
            }
            if !is_zero_to_tol(t - 1.0, EFF_ZERO_SQRT) {
                guesses.push(bez_ivl.evaluate(1.0));
            }
        }

        // use minimization to find point/curve intersection points
        let local_op = if op == SolverOperation::Intersect {
            SolverOperation::Minimize
        } else {
            op
        };

        // init solution search locals
        let mut solution: Solution<O> = Solution::new();
        let mut coin_solutions: Vec<Solution<O>> = Vec::with_capacity(6);
        let mut found_solution: bool;
        let mut last_found_solution = false;
        let mut solution_value = BIG_DOUBLE;
        let mut last_deviation = 0.0;
        let mut last_parameter = 0.0;

        // when intersecting - only try the 1st guess
        let num_guesses = if op == SolverOperation::Intersect {
            1
        } else {
            guesses.len()
        };

        let by_value_main = matches!(
            op,
            SolverOperation::Minimize
                | SolverOperation::Maximize
                | SolverOperation::Rayfire
                | SolverOperation::SignedDirectedMinimize
                | SolverOperation::Signed3dDirectedMinimize
                | SolverOperation::DirectedMinimize
                | SolverOperation::DirectedMaximize
                | SolverOperation::ProjectedMinimize
                | SolverOperation::ProjectedMaximize
                | SolverOperation::ProjectedTangentThroughPoint
        ) || (op == SolverOperation::Normalize
            && requested == SolutionRequested::Single);

        let opt_slice: Option<&[Vec3]> = if opt_vectors.is_empty() {
            None
        } else {
            Some(&opt_vectors)
        };

        let projected_extremum = matches!(
            op,
            SolverOperation::ProjectedMinimize | SolverOperation::ProjectedMaximize
        );
        if projected_extremum {
            // Constrained extrema include finite boundaries as well as stationary
            // points. Enumerate them independently: Newton can find an interior
            // minimum while a requested maximum lies at an endpoint (or vice versa).
            for (parameter, continuity) in [(bez_ivl.min(), conts[0]), (bez_ivl.max(), conts[1])] {
                if omit_projected_boundary(curve, parameter, continuity)
                    || ((continuity as i32) >= (ContinuityType::G1 as i32)
                        && native_join_is_not_extremum(
                            curve,
                            parameter,
                            test_pnt,
                            self.opt(0)?,
                            op == SolverOperation::ProjectedMaximize,
                        ))
                {
                    continue;
                }
                let mut endpoint = Solution::new();
                endpoint.solution_type = SolutionType::SingleValue;
                endpoint.nodes[0] = branch[0];
                endpoint.nodes[1] = branch[1];
                endpoint.num_variables = 1;
                endpoint.num_objects = 1;
                endpoint.objects[0] = curve_obj;
                endpoint.start.parameters[0] = parameter;
                // A failed evaluation is not a distance or an excluded candidate.
                endpoint.start.solution_value = proj_dist(
                    k,
                    curve,
                    parameter,
                    test_pnt,
                    self.opt(0)?,
                    true,
                    Vec3::ZERO,
                    op,
                )?;
                self.add_solution(&endpoint, true)?;
            }
        }

        // for every guess point
        for i in 0..num_guesses {
            let guess = guesses.get(i).copied().unwrap_or(0.0);

            // for projected, directed, and signed pivot cases
            if matches!(
                op,
                SolverOperation::SignedDirectedMinimize
                    | SolverOperation::DirectedMinimize
                    | SolverOperation::DirectedMaximize
                    | SolverOperation::ProjectedMinimize
                    | SolverOperation::ProjectedMaximize
                    | SolverOperation::SignedPivotMinimize
                    | SolverOperation::ProjectedTangentThroughPoint
            ) {
                // locals for LocalPropertyAnalysis call
                let mut curve_prop = CurvePropertyType::PerpendicularToVector;
                let mut vectors = [Vec3::UNDEF; 2];
                let mut direction = Vec3::ZERO;
                let mut property_value = 0.0;

                match op {
                    SolverOperation::DirectedMinimize
                    | SolverOperation::DirectedMaximize
                    | SolverOperation::SignedDirectedMinimize => {
                        curve_prop = CurvePropertyType::PlaneIntersection;
                        vectors[0] = self.opt(1)?.cross(self.opt(0)?);
                        vectors[0] = vectors[0].unitize()?;
                        direction = self.opt(1)?;
                        property_value = -vectors[0].dot(test_pnt);
                    }
                    SolverOperation::SignedPivotMinimize => {
                        curve_prop = CurvePropertyType::CylinderIntersection;
                        vectors[0] = self.opt(1)?; // Pivot
                        vectors[1] = self.opt(0)?; // Normal of projection plane
                        direction = self.opt(2)?; // 'Side' of Pivot pt
                        let v = vectors[1].cross(test_pnt - vectors[0]).cross(vectors[1]);
                        property_value = v.length(); // Radius of an infinite cylinder
                    }
                    SolverOperation::ProjectedMinimize | SolverOperation::ProjectedMaximize => {
                        curve_prop = if op == SolverOperation::ProjectedMaximize {
                            CurvePropertyType::ProjectedPointMaximize
                        } else {
                            CurvePropertyType::ProjectedPointMinimize
                        };
                        vectors[0] = test_pnt;
                        vectors[1] = self.opt(0)?;
                    }
                    SolverOperation::ProjectedTangentThroughPoint => {
                        curve_prop = CurvePropertyType::ProjectedTangentThroughPoint;
                        vectors[0] = test_pnt;
                        vectors[1] = self.opt(0)?;
                    }
                    _ => {}
                }

                // execute NewtonRaphson serarch to find curve param value
                // that sets requested CurveProperty to zero
                match local_property_analysis(
                    k,
                    curve,
                    &bez_ivl,
                    curve_prop,
                    guess,
                    Some(property_value),
                    &vectors,
                    &mut solution,
                ) {
                    Ok(f) => found_solution = f,
                    Err(e) if e.is_unported_or_cancelled() => return Err(e),
                    Err(_) => continue, // If error try another guess
                }

                // A missed stationary point does not rule out another guess.
                // Boundary candidates have already been evaluated independently.
                if !found_solution && projected_extremum {
                    continue;
                }

                if op == SolverOperation::SignedPivotMinimize {
                    // set SolutionValue = CCW angle from sTestPnt to sSolPt
                    let sol_pt = curve.evaluate_point(k, solution.start.parameters[0])?;
                    solution.start.solution_value =
                        proj_signed_angle(test_pnt, sol_pt, vectors[0], vectors[1], direction)?;
                } else if op == SolverOperation::ProjectedTangentThroughPoint {
                    // Ok: The return function value is not the distance between the points.
                } else {
                    // set SolutionValue = appropriate signed testPoint/CurvePoint distance
                    solution.start.solution_value = match proj_dist(
                        k,
                        curve,
                        solution.start.parameters[0],
                        test_pnt,
                        self.opt(0)?,
                        true,
                        direction,
                        op,
                    ) {
                        Ok(d) => d,
                        Err(e) if e.is_unported_or_cancelled() => return Err(e),
                        Err(e) => status_value(&e),
                    };
                }
            } else {
                // else try SmCurve::LocalPointSolve for current guess
                let opts = PointSolveOptions {
                    distance_tolerance: Some(tol_3d),
                    target_distance: Some(self.base.at_distance_sq),
                    vectors: opt_slice,
                };
                match curve.local_point_solve::<O>(k, &bez_ivl, local_op, test_pnt, &opts, guess) {
                    Ok((f, s)) => {
                        found_solution = f;
                        merge_point_solution(&mut solution, f, &s, local_op);
                    }
                    Err(e) if e.is_unported_or_cancelled() => return Err(e),
                    Err(_) => continue, // If have an error move on to next guess
                }
            }

            // quit on 1st guesses that fail (except for SM_SO_NORMALIZE)
            if i == 0 && !found_solution {
                if op == SolverOperation::Normalize {
                    continue;
                }
                return Ok(needs_more_subdivision);
            }

            // skip same answers from different guesses
            if last_found_solution
                && found_solution
                && is_zero(last_deviation - solution.start.solution_value)
                && is_zero(last_parameter - solution.start.parameters[0])
            {
                continue;
            }

            // remember last unique solution
            if found_solution {
                last_found_solution = true;
                last_deviation = solution.start.solution_value;
                last_parameter = solution.start.parameters[0];
            }

            // set outputs
            if i == 0 {
                solution_value = solution.start.solution_value;
            }

            // quit for out-of-tolerance intersection solutions
            if op == SolverOperation::Intersect && solution.start.solution_value > tol_3d {
                return Ok(needs_more_subdivision);
            }

            // move on to next guess points for nonZero RayFires
            if op == SolverOperation::Rayfire {
                let scaled_zero = EFF_ZERO * (1.0 + test_pnt.max_dimension());
                if solution.start.solution_value < -scaled_zero {
                    continue;
                }
            }

            // now Cull Boundary Pt Solutions when real solution will be in next segment
            let mut bndry: Option<(Point3, Vec3)> = None;
            if is_zero(solution.start.parameters[0] - bez_ivl.min())
                && (conts[0] as i32) >= (ContinuityType::G1 as i32)
                && (!projected_extremum || omit_projected_boundary(curve, bez_ivl.min(), conts[0]))
            {
                bndry = Some((start_pt[0], start_pt[1])); // been negated: points outside
            }
            if is_zero(solution.start.parameters[0] - bez_ivl.max())
                && (conts[1] as i32) >= (ContinuityType::G1 as i32)
                && (!projected_extremum || omit_projected_boundary(curve, bez_ivl.max(), conts[1]))
            {
                bndry = Some((end_pt[0], end_pt[1])); // points outside interval
            }

            // when solution is a tolerant boundary solution to a non SM_SO_NORMALIZE problem
            if let Some((bndry_pt, bndry_vec)) = bndry {
                if local_op != SolverOperation::Normalize
                    && !is_zero_to_tol(solution.start.solution_value, EFF_ZERO * 100.0)
                {
                    // skip tolerant answers for problems that should solve exactly
                    if local_op == SolverOperation::AtDistance {
                        continue;
                    }

                    // Solution/BoundaryPt vector
                    let sol_vec = test_pnt - bndry_pt;
                    let mut ang_from_perp_rad = 0.0;
                    let mut outside_dir = false;

                    if matches!(
                        local_op,
                        SolverOperation::Minimize
                            | SolverOperation::Maximize
                            | SolverOperation::Intersect
                            | SolverOperation::Rayfire
                    ) {
                        if sol_vec.length_squared() > EFF_ZERO_SQ {
                            // (status ignored, 0 on failure; passed on with the opt-in
                            // fix, tangency-labs/smlib 8311f8e)
                            ang_from_perp_rad =
                                passed_on(bndry_vec.angle_between(sol_vec), pass_on)?
                                    .unwrap_or(0.0);
                        }
                        ang_from_perp_rad = (ang_from_perp_rad - PI / 2.0).abs();
                        let dot = bndry_vec.dot(sol_vec);
                        outside_dir = dot > 0.0; // TRUE = sSolVec points outside interval
                    } else if local_op == SolverOperation::Signed3dDirectedMinimize {
                        ang_from_perp_rad =
                            passed_on(self.opt(0)?.angle_between(sol_vec), pass_on)?.unwrap_or(0.0);
                    }

                    // skip solutions with tolerant angles that will get better in the
                    // next interval
                    if !is_zero_to_tol(ang_from_perp_rad, EFF_ZERO_RAD) && outside_dir {
                        continue;
                    }
                }
            }

            // If we get an answer back which is the same as the guess point
            // then we might have an interval which matches not just a single point
            if is_zero(solution.start.parameters[0] - guess)
                && is_zero(solution.start.solution_value - solution_value)
            {
                coin_solutions.push(solution);
                continue;
            }

            if !found_solution {
                continue;
            }

            // save the solution
            solution.nodes[0] = branch[0];
            solution.nodes[1] = branch[1];
            solution.num_variables = 1;
            solution.num_objects = 1;
            solution.objects[0] = curve_obj;

            // add solution to sorted sSolution array
            self.add_solution(&solution, by_value_main)?;
        }

        let by_value_coin = matches!(
            op,
            SolverOperation::Minimize
                | SolverOperation::Maximize
                | SolverOperation::SignedDirectedMinimize
                | SolverOperation::Signed3dDirectedMinimize
                | SolverOperation::DirectedMinimize
                | SolverOperation::DirectedMaximize
                | SolverOperation::ProjectedMinimize
                | SolverOperation::ProjectedMaximize
        ) || (op == SolverOperation::Normalize
            && requested == SolutionRequested::Single);

        // If there is less than 3 which are coincident then just put them on to the list.
        // If there are 3 then assume that it is a coincident segment.
        if coin_solutions.len() < 3 {
            for coin in &mut coin_solutions {
                coin.solution_type = SolutionType::SingleValue;
                coin.nodes[0] = branch[0];
                coin.nodes[1] = branch[1];
                coin.num_variables = 1;
                coin.num_objects = 1;
                coin.objects[0] = curve_obj;
            }
            for coin in &coin_solutions {
                self.add_solution(coin, by_value_coin)?;
            }
        } else {
            // Coincident segment found
            let mut num_done = 0;
            let mut done_first = false;
            let mut done_last = false;
            let mut ivl_solution: Solution<O> = Solution::new();
            ivl_solution.solution_type = SolutionType::RangeOfValues;

            for coin in &coin_solutions {
                if is_zero(coin.start.parameters[0] - bez_ivl.min()) {
                    num_done += 1;
                    done_first = true;
                    ivl_solution.start = coin.start;
                }
                if is_zero(coin.start.parameters[0] - bez_ivl.max()) {
                    num_done += 1;
                    done_last = true;
                    ivl_solution.end = coin.start;
                }
            }

            if num_done == 2 && done_first && done_last {
                ivl_solution.nodes[0] = branch[0];
                ivl_solution.nodes[1] = branch[1];
                ivl_solution.num_variables = 1;
                ivl_solution.solution_type = SolutionType::RangeOfValues;
                ivl_solution.num_objects = 1;
                ivl_solution.objects[0] = curve_obj;
                self.add_solution(&ivl_solution, by_value_coin)?;
            }
        }

        Ok(needs_more_subdivision)
    }
}

// ---- GlobalPointSolve -----------------------------------------------------------------

/// The optional inputs of [`Curve::global_point_solve`].
#[derive(Debug, Clone, Copy)]
pub struct GlobalPointSolveOptions<'a, O> {
    /// `cpdOptTargetDistance`: the target distance of
    /// [`SolverOperation::AtDistance`] (required there), or the distance limit of the
    /// minimising and maximising operations.
    pub target_distance: Option<f64>,
    /// `cpOptVectors`: the vectors of the ray, directed and projected operations.
    pub opt_vectors: Option<&'a [Vec3]>,
    /// The handle of this curve, named by its solutions. With this and
    /// [`GlobalPointSolveOptions::analytic_handle`] both `None` the solutions name no
    /// object, but are merged as upstream's (which name the curve).
    pub handle: Option<O>,
    /// The handle of the curve cache's analytic replacement (an `SmLine` the cache
    /// makes of a straight B-spline), named by its fast drop's solutions.
    pub analytic_handle: Option<O>,
    /// `SmTol::GetZoneTol3d(const SmObject*)` of the handles, which a
    /// [`SolverOperation::Rayfire`] solve reads (`GetObjectTolerance`); `None` when
    /// the caller cannot resolve them (that solve then returns
    /// [`SolidError::Unported`]).
    pub zone_tol_3d: Option<ObjectZoneTol<'a, O>>,
}

impl<O> Default for GlobalPointSolveOptions<'_, O> {
    fn default() -> Self {
        Self {
            target_distance: None,
            opt_vectors: None,
            handle: None,
            analytic_handle: None,
            zone_tol_3d: None,
        }
    }
}

/// `DropPointFast(interval, op, point, NULL, tol, target, requested, ..)` of `curve`
/// naming `handle` (base: `SM_ERR`; a line drops the point onto itself, an ellipse
/// handles circles).
fn drop_point_fast<O: Copy>(
    k: &Kernel,
    curve: &Curve,
    interval: &Extent1d,
    op: SolverOperation,
    point: Point3,
    tol: f64,
    target_distance: Option<f64>,
    requested: SolutionRequested,
    handle: Option<O>,
) -> Result<Vec<Solution<O>>> {
    match curve.kind() {
        CurveKind::Line(l) => l.drop_point_fast(interval, op, point, tol, target_distance, handle),
        CurveKind::Ellipse(e) => e.drop_point_fast(
            k,
            interval,
            op,
            point,
            None,
            tol,
            target_distance,
            requested,
            handle,
        ),
        _ => Err(SolidError::sm_err("curve::curve::p1::drop_point_fast")),
    }
}

impl Curve {
    /// `GlobalPointSolve`: the points of the curve over `interval` solving `op` for
    /// `test_point` (closest / farthest, normal feet, intersections within
    /// `dist_tol`, points at a distance, rays, directed and projected extrema). For
    /// [`SolverOperation::Minimize`] and [`SolverOperation::Intersect`] the line drop
    /// is tried first, then the drop onto the cache's analytic replacement (polished
    /// by [`Curve::local_point_solve`]); otherwise the point is solved against every
    /// span of the curve's cache tree by the curve/point global solver.
    ///
    /// # Errors
    ///
    /// `SM_ERR_NULL_POINTER` when the cache cannot be built or
    /// [`SolverOperation::AtDistance`] has no target distance; those of the polishing
    /// local solve; [`SolidError::Unported`] from unported callees. Projected
    /// minimum/maximum queries propagate tree-solve errors, including
    /// `SM_ERR_NOT_WITHIN_TOLERANCE` when an affine optimum cannot be certified
    /// at a representable parameter or the cache does not cover the requested
    /// interval. Other operations retain the legacy policy
    /// of ignoring ordinary tree-solve failures unless status propagation is enabled.
    pub fn global_point_solve<O: SolutionObject>(
        &self,
        k: &Kernel,
        interval: &Extent1d,
        op: SolverOperation,
        test_point: Point3,
        dist_tol: f64,
        requested: SolutionRequested,
        opts: &GlobalPointSolveOptions<'_, O>,
    ) -> Result<SolutionArray<O>> {
        // no handles: upstream's solutions still name this curve (whose natural
        // interval IdenticalSolutions reads); the caller gets them without objects
        if opts.handle.is_none() && opts.analytic_handle.is_none() {
            let ivl = self.natural_interval();
            let inner = GlobalPointSolveOptions {
                target_distance: opts.target_distance,
                opt_vectors: opts.opt_vectors,
                handle: Some(OperandObject::Curve(0, ivl)),
                analytic_handle: Some(OperandObject::Curve(2, ivl)),
                zone_tol_3d: None,
            };
            return Ok(self
                .global_point_solve_cached(
                    k, interval, op, test_point, dist_tol, requested, &inner, &mut None,
                )?
                .without_objects());
        }
        self.global_point_solve_cached(
            k, interval, op, test_point, dist_tol, requested, opts, &mut None,
        )
    }

    /// [`Curve::global_point_solve`] for a caller that names no handles, keeping the
    /// curve's cache in `cache` across calls as [`Curve::global_point_solve_cached`]
    /// does. The cache names the stand-in handle of the unnamed solve (the curve with its
    /// natural interval), so the slot must belong to this curve, unchanged since the
    /// cache was built.
    ///
    /// # Errors
    ///
    /// As [`Curve::global_point_solve`].
    pub(crate) fn global_point_solve_unnamed_cached<O: SolutionObject>(
        &self,
        k: &Kernel,
        interval: &Extent1d,
        op: SolverOperation,
        test_point: Point3,
        dist_tol: f64,
        requested: SolutionRequested,
        target_distance: Option<f64>,
        cache: &mut Option<CurveCache<OperandObject>>,
    ) -> Result<SolutionArray<O>> {
        let ivl = self.natural_interval();
        let inner = GlobalPointSolveOptions {
            target_distance,
            opt_vectors: None,
            handle: Some(OperandObject::Curve(0, ivl)),
            analytic_handle: Some(OperandObject::Curve(2, ivl)),
            zone_tol_3d: None,
        };
        Ok(self
            .global_point_solve_cached(
                k, interval, op, test_point, dist_tol, requested, &inner, cache,
            )?
            .without_objects())
    }

    /// [`Curve::global_point_solve`] keeping the curve's cache in `cache` across calls,
    /// as upstream keeps it on the curve (`GetOrCreateObjectCache`): the cache is built
    /// where the solve first needs it (with the same errors) and reused by later calls
    /// given the same slot. The slot must belong to this curve, unchanged since the
    /// cache was built, and to the same `opts.handle` (the cache names it in its
    /// solutions); the build is deterministic, so a reused cache holds the same spans,
    /// tree and analytic replacement as a fresh one.
    ///
    /// # Errors
    ///
    /// As [`Curve::global_point_solve`].
    pub(crate) fn global_point_solve_cached<O: SolutionObject>(
        &self,
        k: &Kernel,
        interval: &Extent1d,
        op: SolverOperation,
        test_point: Point3,
        dist_tol: f64,
        requested: SolutionRequested,
        opts: &GlobalPointSolveOptions<'_, O>,
        cache: &mut Option<CurveCache<O>>,
    ) -> Result<SolutionArray<O>> {
        let target = opts.target_distance;

        // Try quick solves when MINIMIZING and INTERSECTING (without using curve cache)
        if op == SolverOperation::Minimize || op == SolverOperation::Intersect {
            if let Some(sols) = ignored(drop_point_fast(
                k,
                self,
                interval,
                op,
                test_point,
                dist_tol,
                target,
                requested,
                opts.handle,
            ))? {
                let mut out = SolutionArray::new();
                for s in &sols {
                    out.add(s);
                }
                return Ok(out);
            }
        }

        // get target curve's cache
        let mut solutions: SolutionArray<O> = SolutionArray::new();
        if cache.is_none() {
            *cache = Some(match make_curve_cache::<O>(k, self, opts.handle) {
                Ok(c) => c,
                Err(e) if e.is_unported_or_cancelled() => return Err(e),
                Err(_) => return Err(SolidError::null_pointer(OP_GLOBAL)),
            });
        }
        let cache = cache.as_mut().ok_or(SolidError::null_pointer(OP_GLOBAL))?;
        if matches!(
            op,
            SolverOperation::ProjectedMinimize | SolverOperation::ProjectedMaximize
        ) && !interval.is_contained_by(cache.interval(), 0.0)
        {
            // A result over only the cached intersection cannot certify the
            // caller's whole request. Do not silently discard its remainder.
            return Err(SolidError::Upstream {
                operation: OP_GLOBAL,
                code: crate::error::SM_ERR_NOT_WITHIN_TOLERANCE,
                entity: None,
            });
        }

        // Analytical Curve Optimization
        if op == SolverOperation::Minimize || op == SolverOperation::Intersect {
            if let Some(analytic) = cache.analytic_curve(k)? {
                // (the cache's analytic curve is this curve when it is analytic)
                let analytic_handle = if self.is_analytic() {
                    opts.handle
                } else {
                    opts.analytic_handle
                };
                if let Some(sols) = ignored(drop_point_fast(
                    k,
                    analytic,
                    interval,
                    op,
                    test_point,
                    dist_tol,
                    target,
                    requested,
                    analytic_handle,
                ))? {
                    for s in &sols {
                        solutions.add(s);
                    }

                    // when there are solutions that need tightening
                    if let Some(first) = solutions.get(0).copied() {
                        let anal_value = first.start.parameters[0];
                        let anal_solution_value = first.start.solution_value;

                        // when solution could be tighter
                        if anal_solution_value > EFF_ZERO {
                            let (found, sol) = self.local_point_solve::<O>(
                                k,
                                interval,
                                op,
                                test_point,
                                &PointSolveOptions::default(),
                                anal_value,
                            )?;

                            // when LocalPointSolve converged to a better answer
                            if found && sol.start.solution_value < anal_solution_value {
                                if let Some(slot) = solutions.get_mut(0) {
                                    *slot = sol;
                                }
                            }
                        }
                    }

                    // all done
                    return Ok(solutions);
                }
            }
        }

        // Execute a Point/CurveLeafNode solution for every LeafNode in the target
        // curve's cache.
        let point_tree: Tree<O> = Tree::from_point(test_point, 0.0, false);
        let curve_tree = cache.tree().ok_or(SolidError::null_pointer(OP_GLOBAL))?;

        // Construct and init global solver
        let mut gs = CpGlobalSolver::new(k, &point_tree, curve_tree, *interval).with_curve(self);
        if matches!(
            op,
            SolverOperation::ProjectedMinimize | SolverOperation::ProjectedMaximize
        ) && matches!(
            requested,
            SolutionRequested::Single | SolutionRequested::All
        ) {
            if let Some(spans) = complete_affine_definition(self, *interval) {
                gs.projected_affine_selection = Some(ProjectedAffineSelection {
                    candidates: Vec::new(),
                    maximum_upper: f64::NEG_INFINITY,
                    minimum_lower: f64::INFINITY,
                    span_lower: 0.0,
                    pending_spans: Some(spans),
                });
            }
        }
        gs.base.fixes = k.fixes;
        gs.zone_tol_3d = opts.zone_tol_3d;

        // set atDistance value - when appropriate
        if op == SolverOperation::AtDistance {
            let d = target.ok_or(SolidError::null_pointer(OP_GLOBAL))?;
            gs.base.set_at_distance(d);
        }

        // Set dBestAnswerSoFarSq for various cases
        let best_answer_so_far_sq = match op {
            SolverOperation::Rayfire
            | SolverOperation::Minimize
            | SolverOperation::DirectedMinimize
            | SolverOperation::ProjectedMinimize
            | SolverOperation::SignedDirectedMinimize
            | SolverOperation::Signed3dDirectedMinimize => match target {
                Some(d) => {
                    if d >= 0.0 {
                        d * d
                    } else {
                        -d * d
                    }
                }
                None => BIG_DOUBLE,
            },
            SolverOperation::Maximize
            | SolverOperation::DirectedMaximize
            | SolverOperation::ProjectedMaximize => match target {
                Some(d) => d * d,
                None => 0.0,
            },
            _ => BIG_DOUBLE,
        };

        // execute the spatial decomposition tree leafNode/Point solutions (status
        // ignored; passed on with the opt-in fix, tangency-labs/smlib 8311f8e)
        let request = SolveRequest {
            operation: op,
            solution_requested: requested,
            tolerance_3d: dist_tol,
            best_answer_so_far_sq,
            opt_vectors: opts.opt_vectors.unwrap_or(&[]),
        };
        let solved = gs.solve_it(&request, &mut solutions);
        if matches!(
            op,
            SolverOperation::ProjectedMinimize | SolverOperation::ProjectedMaximize
        ) {
            solved?;
        } else {
            ignore_or_pass_on(solved, k.fixes.query_statuses_pass_on)?;
        }

        if let Some(selection) = &gs.projected_affine_selection {
            selection.finish(
                self,
                test_point,
                gs.opt(0)?,
                op == SolverOperation::ProjectedMaximize,
                requested,
                target,
                dist_tol,
                &mut solutions,
                |a, b| Ok(gs.identical_solutions(a, b)?.0 == SolutionPairType::Identical),
            )?;
        }

        Ok(solutions)
    }
}

// ---- TessellateByBisection / Tessellate ---------------------------------------------

impl Curve {
    /// `TessellateByBisection`: sample points of the curve over `interval` whose
    /// segments pass the chord height ratio, tangent angle (degrees) and maximum 3D
    /// distance tests ([`test_step`]; a zero tolerance is ignored). A degenerate curve
    /// gives its start point (zero tangent); a line, when no maximum distance is
    /// asked for, its two ends (with the start tangent). Otherwise steps grow by
    /// doubling and shrink by halving from the whole interval, each accepted or
    /// growing step being re-tested at three points per knot span it covers.
    ///
    /// # Errors
    ///
    /// Those of the evaluations and of the degeneracy and line tests.
    pub fn tessellate_by_bisection(
        &self,
        k: &Kernel,
        interval: &Extent1d,
        chord_height_tol: f64,
        angle_tol_deg: f64,
        max_3d_distance: f64,
    ) -> Result<Tessellation> {
        let mut out = Tessellation::default();

        // locals
        let max_3d_distance_squared = max_3d_distance * max_3d_distance;

        // SmCurve locals - knots (status ignored)
        let (knots, _) = knots_ignoring_status(self.knots_with_status(k, None))?;

        // get this SmCurve's start position and derivative
        let mut cur_pv: [Vec3; 2] = eval_n(k, self, interval.min(), 1, true, true)?; // nonZeroTangent values

        // when SmCurve's start derivative is zero
        if cur_pv[1].length_squared() < EFF_ZERO_SQ {
            // get SmCurve's pos and der values near beginning
            let tmp: [Vec3; 2] = eval_n(k, self, interval.evaluate(0.001), 1, true, true)?;

            // if der near end is zero - make one up to continue processing
            if tmp[1].length_squared() < EFF_ZERO_SQ {
                cur_pv[1] = Vec3::new(1.0, 1.0, 1.0); // Recover by adding arbitrary vector here
            } else {
                cur_pv[1] = tmp[1];
            }
        }

        // check for degenerate or linear curves
        if self.is_degenerate(k, EFF_ZERO, None)? {
            // degenerate curve is a point - load one tessellation point and exit
            out.parameters.push(interval.min());
            out.points.push(cur_pv[0]);
            out.tangents.push(Vec3::ZERO);
            return Ok(out);
        }

        // If linear, and we're not concerned with segment length, return the whole line.
        let tol = if chord_height_tol > EFF_ZERO {
            chord_height_tol / 100.0
        } else {
            EFF_ZERO * 100.0
        };
        if max_3d_distance < EFF_ZERO && self.is_line(k, 5, tol)?.is_some() {
            // load two tessellation points and exit
            out.parameters.push(interval.min());
            out.parameters.push(interval.max());
            out.points.push(cur_pv[0]);
            cur_pv[0] = self.evaluate_point(k, interval.max())?;
            out.points.push(cur_pv[0]);
            out.tangents.push(cur_pv[1]);
            out.tangents.push(cur_pv[1]);
            return Ok(out);
        }

        // Get step size estimates
        let min_step_size = interval.length() / 100000.0;

        // load start param and value into output arrays
        output(
            cur_pv,
            interval.min(),
            Some(&mut out.parameters),
            Some(&mut out.points),
            Some(&mut out.tangents),
        );

        // begin at the curve start with step size set to interval size
        let mut cur_t = interval.min();
        let mut param_step_size = interval.length();
        let scaled_zero = EFF_ZERO * (1.0 + interval.max_dimension());

        // locals
        let mut done = false;
        let mut last_step_size = 0.0;
        let mut num_segments: u64 = 0;

        // search for next tessellation point
        while !done {
            num_segments += 1;

            let mut next_pv: [Vec3; 2];
            let mut step_direction: i64 = 0;
            let mut tests_made: u64 = 0;
            loop {
                // count the number of iterations to find the next tessellation point
                tests_made += 1;

                // limit the next step size to the interval's max value
                if cur_t + param_step_size > interval.max() {
                    param_step_size = interval.max() - cur_t;
                }

                // make sure that our last step is not too close to the end of the curve
                if cur_t + param_step_size < interval.max() - scaled_zero
                    && cur_t + 2.0 * param_step_size > interval.max() + scaled_zero
                {
                    param_step_size = interval.max() - cur_t;
                }

                next_pv = eval_n(k, self, cur_t + param_step_size, 1, true, true)?; // nonZeroTangent values
                if next_pv[1].length_squared() < EFF_ZERO_SQ {
                    let tmp: [Vec3; 2] =
                        eval_n(k, self, cur_t + param_step_size * 0.99, 1, true, true)?;
                    if tmp[1].length_squared() < EFF_ZERO_SQ {
                        next_pv[1] = cur_pv[1];
                    } else {
                        next_pv[1] = tmp[1];
                    }
                }

                if param_step_size > min_step_size - EFF_ZERO {
                    let mut result = test_step(
                        param_step_size,
                        cur_pv,
                        next_pv,
                        chord_height_tol,
                        angle_tol_deg,
                        max_3d_distance_squared,
                    );

                    // If we don't have a lot of segments force some additional testing.
                    if result == TessStepResult::Accept || result == TessStepResult::Increase {
                        // generate 3 sample points per knot interval
                        let params = find_test_params(cur_t, cur_t + param_step_size, &knots);
                        for &tmp_t in &params {
                            let next_pv2: [Vec3; 2] = eval_n(k, self, tmp_t, 1, true, true)?; // nonZeroTangent values
                            if next_pv2[1].length_squared() < EFF_ZERO_SQ {
                                continue;
                            }
                            let result2 = test_step(
                                tmp_t - cur_t,
                                cur_pv,
                                next_pv2,
                                chord_height_tol,
                                angle_tol_deg,
                                max_3d_distance_squared,
                            );
                            if result2 == TessStepResult::Cut {
                                result = result2;
                                break;
                            }
                        }
                    }

                    if tests_made > 100 {
                        break; // Just accept what we have so far
                    }

                    if result == TessStepResult::Increase
                        && num_segments == 1
                        && are_same(param_step_size, interval.length())
                    {
                        let test_tol = EFF_ZERO * (1.0 + cur_pv[0].max_dimension());
                        if next_pv[0].distance(cur_pv[0]) < test_tol {
                            result = TessStepResult::Cut;
                        }
                    }
                    if result == TessStepResult::Accept {
                        break;
                    }
                    if result == TessStepResult::Cut {
                        if step_direction == 1 {
                            param_step_size = last_step_size;
                            next_pv = eval_n(k, self, cur_t + param_step_size, 1, true, true)?;
                            break;
                        }
                        param_step_size /= 2.0;

                        step_direction = -1;
                        if param_step_size < min_step_size {
                            param_step_size = min_step_size;
                            next_pv = eval_n(k, self, cur_t + param_step_size, 1, true, true)?;
                            break;
                        }
                        continue;
                    } else if result == TessStepResult::Increase {
                        step_direction = 1;
                        last_step_size = param_step_size;

                        // If the current step goes effectively all the way to the end,
                        // accept it: evaluate there and break the while loop.
                        if cur_t + param_step_size > interval.max() - scaled_zero {
                            param_step_size = interval.max() - cur_t;
                            next_pv = eval_n(k, self, interval.max(), 1, true, true)?;
                            break;
                        }

                        param_step_size *= 2.0; // Double the step.

                        // If the doubled step goes past end, set it to the end
                        if cur_t + param_step_size > interval.max() + EFF_ZERO {
                            param_step_size = interval.max() - cur_t;
                        }

                        continue;
                    }
                } else if param_step_size < min_step_size
                    && param_step_size < interval.max() - cur_t
                {
                    param_step_size = min_step_size;
                    next_pv = eval_n(k, self, cur_t + param_step_size, 1, true, true)?;
                    break;
                }

                // If we made it to here no adjustment was required
                break;
            }

            // If made it to here then the point was accepted.
            cur_t += param_step_size;
            output(
                next_pv,
                cur_t,
                Some(&mut out.parameters),
                Some(&mut out.points),
                Some(&mut out.tangents),
            );
            cur_pv = next_pv;
            if are_same(cur_t, interval.max()) {
                done = true;
            }
            let est_step = param_step_size;
            let max_step = interval.max() - cur_t; // Max step possible
            param_step_size = max_step.min(est_step);
        }

        Ok(out)
    }

    /// `Tessellate` (virtual): the line's even parameter steps (its tolerances
    /// ignored), the B-spline override (the same body as the base), or the base body
    /// [`Curve::tessellate_base`].
    ///
    /// # Errors
    ///
    /// Those of [`Curve::tessellate_base`] and of the B-spline override.
    pub fn tessellate(
        &self,
        k: &Kernel,
        interval: &Extent1d,
        chord_height_tol: f64,
        angle_tol_deg: f64,
        min_segments: u64,
    ) -> Result<Tessellation> {
        match self.kind() {
            CurveKind::Line(l) => l.tessellate(
                interval,
                min_segments,
                TessellationRequest {
                    parameters: true,
                    points: true,
                    tangents: true,
                },
            ),
            CurveKind::BSpline(c) => c.tessellate(
                k,
                interval,
                chord_height_tol,
                angle_tol_deg,
                min_segments as usize,
            ),
            CurveKind::Hyperbola(c) => c.base().tessellate(
                k,
                interval,
                chord_height_tol,
                angle_tol_deg,
                min_segments as usize,
            ),
            CurveKind::Ellipse(c) => c.base().tessellate(
                k,
                interval,
                chord_height_tol,
                angle_tol_deg,
                min_segments as usize,
            ),
            CurveKind::Parabola(c) => c.base().tessellate(
                k,
                interval,
                chord_height_tol,
                angle_tol_deg,
                min_segments as usize,
            ),
            CurveKind::Hermite(_)
            | CurveKind::Conic(_)
            | CurveKind::Projected(_)
            | CurveKind::Iso(_)
            | CurveKind::Offset(_)
            | CurveKind::CrvInVolume(_)
            | CurveKind::CrvOnSurf(_)
            | CurveKind::TangentField(_)
            | CurveKind::Composite(_) => {
                self.tessellate_base(k, interval, chord_height_tol, angle_tol_deg, min_segments)
            }
        }
    }

    /// `SmCurve::Tessellate` (the base-class body): [`Curve::tessellate_by_bisection`]
    /// with the maximum segment length `ApproximateLength(5) / min_segments` (none
    /// for at most one segment).
    ///
    /// # Errors
    ///
    /// Those of [`Curve::tessellate_by_bisection`] and the length.
    pub fn tessellate_base(
        &self,
        k: &Kernel,
        interval: &Extent1d,
        chord_height_tol: f64,
        angle_tol_deg: f64,
        min_segments: u64,
    ) -> Result<Tessellation> {
        // init MaxLength control parameter
        let max_length = if min_segments > 1 {
            self.approximate_length(k, interval, 5, false)?.0 / (min_segments as f64 * 1.0)
        } else {
            0.0
        };

        // pass the call along
        self.tessellate_by_bisection(k, interval, chord_height_tol, angle_tol_deg, max_length)
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn complete_affine_selection_excludes_mixed_weight_spans() {
        use crate::nurbs::math::CPoint;
        for (last_weight, supported) in [(1.0, true), (2.0, false)] {
            let net = crate::nurbs::crv_basic::Curve::new(
                vec![
                    CPoint::new(0.0, 0.0, 0.0, 1.0),
                    CPoint::new(1.0, 0.0, 0.0, 1.0),
                    CPoint::new(2.0, 1.0, 0.0, last_weight),
                ],
                1,
                vec![0.0, 0.0, 1.0, 2.0, 2.0],
            )
            .unwrap();
            let curve = Curve::from(BSplineCurve::from_nurb(net, 3, true));
            assert_eq!(
                complete_affine_definition(&curve, Extent1d::from_min_max(0.0, 2.0)).is_some(),
                supported
            );
            for outside in [
                Extent1d::from_min_max(-1.0, 2.0),
                Extent1d::from_min_max(0.0, 3.0),
            ] {
                assert!(complete_affine_definition(&curve, outside).is_none());
            }
        }
    }

    #[test]
    fn projected_plane_uv_boundary_uses_defining_join_not_cache_angle() {
        use crate::curve::crv_on_surf::CrvOnSurf;
        use crate::geom::extent2d::Extent2d;
        use crate::math::Vec2;
        use crate::nurbs::math::{CPoint, NOW};
        use crate::surface::{plane::Plane, surface::p2::Surface};
        let k = Kernel::default();
        for (height, domain_y, kept_definition, smooth) in [
            (0.999999, 2.0, false, false),
            (1.0, 2.0, false, true),
            (1.0, 1.0, false, false),
            (1.0, 2.0, true, false),
        ] {
            let uv = Curve::from(BSplineCurve::from_nurb(
                crate::nurbs::crv_basic::Curve::new(
                    vec![
                        CPoint::new(-0.001, height, 0.0, NOW),
                        CPoint::new(0.0, 1.0, 0.0, NOW),
                        CPoint::new(0.001, height, 0.0, NOW),
                    ],
                    1,
                    vec![0.0, 0.0, 1.0, 2.0, 2.0],
                )
                .unwrap(),
                2,
                false,
            ));
            let domain = Extent2d::from_points(Vec2::new(-2.0, -2.0), Vec2::new(2.0, domain_y));
            let mut surface = Surface::from(Plane::new(
                Vec3::ZERO,
                Vec3::X,
                Vec3::Y,
                Vec2::new(1.0, 1.0),
                &domain,
            ));
            if kept_definition {
                use crate::surface::procedural_blend::corner::{CornerDefinition, EndDefinition};
                surface.set_corner_definition(CornerDefinition::End(EndDefinition {
                    support: Box::new(surface.clone()),
                    boundary: Vec::new(),
                }));
            }
            let curve = Curve::from(CrvOnSurf::new(&k, uv, surface, None, 2).unwrap());
            assert_eq!(
                omit_projected_boundary(&curve, 1.0, ContinuityType::G1),
                smooth
            );
            if domain_y == 2.0 && !kept_definition {
                // The same underlying affine span on either side of an
                // artificial subdivision remains exactly smooth.
                assert!(omit_projected_boundary(&curve, 0.5, ContinuityType::G1));
                assert!(omit_projected_boundary(&curve, 1.5, ContinuityType::G1));
            }
        }
    }

    #[test]
    fn affine_point_enclosure_uses_actual_right_value_at_a_jump() {
        use crate::nurbs::math::CPoint;
        let net = crate::nurbs::crv_basic::Curve::new(
            vec![
                CPoint::new(100.0, 0.0, 0.0, 1.0),
                CPoint::new(1.0, 0.0, 0.0, 1.0),
                CPoint::new(20.0, 0.0, 0.0, 2.0),
                CPoint::new(10.0, 1.0, 0.0, 1.0),
            ],
            1,
            vec![0.0, 0.0, 1.0, 1.0, 2.0, 2.0],
        )
        .unwrap();
        let curve = Curve::from(BSplineCurve::from_nurb(net, 3, true));
        let at = affine_projected_point_enclosure(&curve, 1.0, Vec3::ZERO, Vec3::Z).unwrap();
        assert!(at.lo <= 10.0 && at.hi >= 10.0 && at.lo > 9.0);
        let left =
            affine_projected_point_enclosure(&curve, next_parameter_down(1.0), Vec3::ZERO, Vec3::Z)
                .unwrap();
        assert!(left.lo < 1.000001 && left.hi < 1.000001);
        let end = affine_projected_point_enclosure(&curve, 2.0, Vec3::ZERO, Vec3::Z).unwrap();
        assert!(end.lo <= 101.0_f64.sqrt() && end.hi >= 101.0_f64.sqrt());
        assert!(affine_projected_point_enclosure(&curve, 3.0, Vec3::ZERO, Vec3::Z).is_none());
        let line = Curve::from(crate::curve::line::Line::infinite(
            Vec3::new(1.0, 2.0, 3.0),
            Vec3::X * 2.0,
            3,
            false,
        ));
        let metric = affine_projected_point_enclosure(&line, 0.5, Vec3::ZERO, Vec3::Z).unwrap();
        assert!(metric.lo <= 8.0_f64.sqrt() && metric.hi >= 8.0_f64.sqrt());
    }
    #[test]
    fn carrier_directed_arithmetic_covers_binary64_boundaries() {
        let third = carrier_div(Ivl::from(1.0), Ivl::from(3.0));
        // FMA retains the tiny signed differences that separately rounded
        // endpoint*3 can round back to exactly one.
        assert!(third.lo.mul_add(3.0, -1.0) < 0.0 && third.hi.mul_add(3.0, -1.0) > 0.0);
        let tiny = f64::from_bits(1);
        let half = carrier_mul(Ivl::from(tiny), Ivl::from(0.5));
        assert!(half.lo <= 0.0 && half.hi >= tiny);
        let zero = carrier_sub(Ivl::from(-0.0), Ivl::from(0.0));
        assert_eq!(zero, Ivl::from(0.0));
        let overflow = carrier_mul(Ivl::from(f64::MAX), Ivl::from(2.0));
        assert_eq!(overflow, Ivl::whole());
        assert_eq!(
            carrier_div(Ivl::from(1.0), Ivl::new(-0.0, 0.0)),
            Ivl::whole()
        );
        assert_eq!(
            carrier_sub(Ivl::from(f64::INFINITY), Ivl::from(f64::INFINITY)),
            Ivl::whole()
        );
        let product = carrier_mul(
            Ivl::from(next_parameter_up(1.0)),
            Ivl::from(next_parameter_down(1.0)),
        );
        assert!(product.lo < 1.0 && product.hi > 1.0);
        let length = carrier_length([Ivl::from(3.0), Ivl::from(4.0), Ivl::from(0.0)]);
        assert!(length.lo <= 5.0 && length.hi >= 5.0);
    }
    #[test]
    fn cubic_c0_join_proof_distinguishes_straight_and_small_kink() {
        for height in [0.0, 1.0 / 1024.0] {
            let points = (-3_i32..=3)
                .map(|x| {
                    crate::nurbs::math::CPoint::new(
                        f64::from(x),
                        f64::from(x.abs()) * height,
                        0.0,
                        crate::nurbs::math::NOW,
                    )
                })
                .collect();
            let net = crate::nurbs::crv_basic::Curve::new(
                points,
                3,
                vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0],
            )
            .unwrap();
            let curve = Curve::from(BSplineCurve::from_nurb(net, 3, false));
            assert_eq!(exact_native_spline_join(&curve, 1.0, false), height == 0.0);
            assert!(exact_native_spline_join(&curve, 0.5, false));
        }
    }

    #[test]
    fn affine_join_requires_exact_position_direction_and_orientation() {
        use crate::nurbs::math::CPoint;
        let join = |points: Vec<CPoint>, repeated: bool| {
            let knots = if repeated {
                vec![0.0, 0.0, 1.0, 1.0, 2.0, 2.0]
            } else {
                vec![0.0, 0.0, 1.0, 2.0, 2.0]
            };
            let net = crate::nurbs::crv_basic::Curve::new(points, 1, knots).unwrap();
            Curve::from(BSplineCurve::from_nurb(net, 3, false))
        };
        let p = |x, y, w| CPoint::new(x, y, 0.0, w);
        for scale in [f64::from_bits(1), 1.0, 2.0_f64.powi(500)] {
            let straight = join(
                vec![p(-scale, 0.0, 1.0), p(0.0, 0.0, 1.0), p(scale, 0.0, 1.0)],
                false,
            );
            assert!(exact_affine_join(&straight, 1.0));
            let kink = join(
                vec![
                    p(-scale, scale, 1.0),
                    p(0.0, 0.0, 1.0),
                    p(scale, scale, 1.0),
                ],
                false,
            );
            assert!(!exact_affine_join(&kink, 1.0));
        }
        let cusp = join(
            vec![p(1.0, 0.0, 1.0), p(0.0, 0.0, 1.0), p(1.0, 0.0, 1.0)],
            false,
        );
        assert!(!exact_affine_join(&cusp, 1.0));
        let weighted = join(
            vec![
                p(0.0, 0.0, 1.0),
                p(1.0, 0.0, 1.0),
                p(2.0, 0.0, 2.0),
                p(4.0, 0.0, 2.0),
            ],
            true,
        );
        assert!(exact_affine_join(&weighted, 1.0));
        let gap = join(
            vec![
                p(0.0, 0.0, 1.0),
                p(1.0, 0.0, 1.0),
                p(next_parameter_up(2.0), 0.0, 2.0),
                p(4.0, 0.0, 2.0),
            ],
            true,
        );
        assert!(!exact_affine_join(&gap, 1.0));
        // Rounded direction subtraction loses the middle control's tiny offset.
        let tiny = 2.0_f64.powi(-60);
        let rounded = join(
            vec![p(-1.0, -1.0, 1.0), p(0.0, tiny, 1.0), p(1.0, 1.0, 1.0)],
            false,
        );
        assert!(!exact_affine_join(&rounded, 1.0));
    }

    #[test]
    fn affine_parallel_identity_is_orientation_and_scale_independent() {
        for scale in [
            f64::from_bits(1),
            2.0_f64.powi(-500),
            1.0,
            2.0_f64.powi(500),
            f64::MAX,
        ] {
            assert!(exact_parallel_difference(
                [0.0; 3],
                [scale, -scale, scale],
                Vec3::new(1.0, -1.0, 1.0)
            ));
            assert!(exact_parallel_difference(
                [-scale, scale, -scale],
                [scale, -scale, scale],
                Vec3::new(1.0, -1.0, 1.0)
            ));
        }
        assert!(exact_parallel_difference(
            [3.0, -7.0, 9.0],
            [5.0, -4.0, 14.0],
            Vec3::new(-4.0, -6.0, -10.0)
        ));
        assert!(exact_parallel_difference(
            [0.0; 3],
            [1.0, 1.0, 0.0],
            Vec3::new(1.0, 1.0, 0.0)
        ));
    }

    #[test]
    fn affine_parallel_identity_rejects_near_parallel_and_nonfinite() {
        assert!(!exact_parallel_difference(
            [0.0; 3],
            [1.0, next_parameter_up(1.0), 0.0],
            Vec3::new(1.0, 1.0, 0.0)
        ));
        // Ordinary products underflow here; the exact determinant is nonzero.
        let tiny = f64::from_bits(1);
        assert!(!exact_parallel_difference(
            [0.0; 3],
            [tiny, tiny * 2.0, 0.0],
            Vec3::new(tiny, tiny, 0.0)
        ));
        // Ordinary products overflow; exact cancellation still establishes parallelism.
        assert!(exact_parallel_difference(
            [0.0; 3],
            [f64::MAX, f64::MAX, 0.0],
            Vec3::new(f64::MAX, f64::MAX, 0.0)
        ));
        assert!(!exact_parallel_difference(
            [0.0; 3],
            [f64::MAX, next_parameter_down(f64::MAX), 0.0],
            Vec3::new(f64::MAX, f64::MAX, 0.0)
        ));
        assert!(!exact_parallel_difference(
            [0.0; 3],
            [f64::INFINITY, 1.0, 0.0],
            Vec3::Z
        ));
    }
    use super::*;

    #[test]
    fn projected_query_refuses_incomplete_cached_domain() {
        use crate::curve::curve_cache::CurveCacheOptions;
        use crate::topology::ObjectHandle;
        let k = Kernel::default();
        let curve = Curve::from(BSplineCurve::create_line_segment(3, Vec3::ZERO, Vec3::X).unwrap());
        for op in [
            SolverOperation::ProjectedMinimize,
            SolverOperation::ProjectedMaximize,
        ] {
            let result = curve.global_point_solve::<ObjectHandle>(
                &k,
                &Extent1d::from_min_max(0.0, 2.0),
                op,
                Vec3::X * 2.0,
                1e-9,
                SolutionRequested::Single,
                &GlobalPointSolveOptions {
                    opt_vectors: Some(&[Vec3::Z]),
                    ..Default::default()
                },
            );
            assert!(
                matches!(result, Err(e) if e.kind() == crate::error::ErrorKind::NotWithinTolerance)
            );
        }
        let mut cache = Some(CurveCache::<ObjectHandle>::new(
            curve.clone(),
            Extent1d::from_min_max(0.0, 0.5),
            CurveCacheOptions::default(),
        ));
        for op in [
            SolverOperation::ProjectedMinimize,
            SolverOperation::ProjectedMaximize,
        ] {
            let result = curve.global_point_solve_cached(
                &k,
                &Extent1d::from_min_max(0.0, 1.0),
                op,
                Vec3::X,
                1e-9,
                SolutionRequested::Single,
                &GlobalPointSolveOptions {
                    opt_vectors: Some(&[Vec3::Z]),
                    ..Default::default()
                },
                &mut cache,
            );
            assert!(
                matches!(result, Err(e) if e.kind() == crate::error::ErrorKind::NotWithinTolerance)
            );
        }
    }

    #[test]
    fn smooth_adjacent_float_span_cannot_drop_an_interior_minimum() {
        use crate::topology::ObjectHandle;
        let k = Kernel::default();
        let lo = 1.0;
        let hi = next_parameter_up(lo);
        let spline = crate::nurbs::crv_basic::Curve::new(
            vec![
                crate::nurbs::math::CPoint::new(-1.0, 1.0, 0.0, 1.0),
                crate::nurbs::math::CPoint::new(1.0, 1.0, 0.0, 1.0),
            ],
            1,
            vec![lo, lo, hi, hi],
        )
        .unwrap();
        // Preserve the raw adjacent knots: an owned NURB normalizes them.
        let spline = BSplineCurve::from_nurb(spline, 3, true);
        assert_eq!(spline.gw_nurb().unwrap().knots(), &[lo, lo, hi, hi]);
        let curve = Curve::from(spline);
        let point_tree = Tree::<ObjectHandle>::new(None);
        let curve_tree = Tree::<ObjectHandle>::new(None);
        let interval = Extent1d::from_min_max(lo, hi);
        let mut solver = CpGlobalSolver::new(&k, &point_tree, &curve_tree, interval);
        solver.base.operation = SolverOperation::ProjectedMinimize;
        solver.base.tolerance_3d = 1e-12;
        solver.base.opt_vectors = vec![Vec3::Z];
        let result = solver.projected_affine(
            &k,
            &curve,
            interval,
            [ContinuityType::G1; 2],
            Vec3::ZERO,
            &[None; crate::geom::global_solver::MAX_SOLVER_TREES],
            None,
        );
        assert!(
            matches!(result, Err(e) if e.kind() == crate::error::ErrorKind::NotWithinTolerance)
        );
    }

    /// tangency-labs/solidcore#527 (WR-383): a forward range whose end sits 0.001 off
    /// the crossing; the local intersection moves the end by 0.001 of the range, which
    /// upstream's test (|frac| near 0) refuses and the opt-in (|1 - frac|) accepts.
    #[test]
    fn issue_527_refine_cc_range_end() {
        use crate::geom::attribute::ObjectKey;
        use crate::geom::solution_array::{ObjectKind, Solution, SolutionType};
        use crate::topology::ObjectHandle;
        let line = |a: [f64; 3], b: [f64; 3]| {
            let p = |c: [f64; 3]| Vec3::new(c[0], c[1], c[2]);
            Curve::from(BSplineCurve::create_line_segment(3, p(a), p(b)).unwrap())
        };
        // crv1 along x; crv2 crosses it at x = 0.999, at t2 = 1 / 1.001
        let crv1 = line([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let crv2 = line([0.999, 1.0, 0.0], [0.999, -0.001, 0.0]);
        let h = |key| ObjectHandle {
            kind: ObjectKind::Curve,
            key: ObjectKey(key),
        };
        let run = |k: &Kernel| {
            let mut sol = Solution::<ObjectHandle>::new();
            sol.solution_type = SolutionType::RangeOfValues;
            sol.num_variables = 2;
            sol.num_objects = 2;
            sol.objects[0] = Some(h(1));
            sol.objects[1] = Some(h(2));
            sol.start.parameters[0] = 0.0;
            sol.start.parameters[1] = 0.0;
            sol.start.solution_value = 0.0;
            sol.end.parameters[0] = 1.0;
            sol.end.parameters[1] = 1.0;
            sol.end.solution_value = 0.001;
            let mut sols = SolutionArray::with_capacity(1);
            sols.try_set_size(1).unwrap();
            sols.set_at(0, &sol).unwrap();
            refine_cc_solutions(k, &mut sols, |o| match o.map(|o| o.key) {
                Some(ObjectKey(1)) => Some(&crv1),
                Some(ObjectKey(2)) => Some(&crv2),
                _ => None,
            })
            .unwrap();
            sols.get(0).unwrap().end
        };
        let mut k = Kernel::default();
        let end = run(&k);
        assert_eq!((end.parameters[0], end.solution_value), (1.0, 0.001));
        k.fixes.refine_cc_range_ends = true;
        let end = run(&k);
        assert!((end.parameters[0] - 0.999).abs() < 1e-9, "{end:?}");
        assert!(end.solution_value < 1e-9, "{end:?}");
    }

    fn cubic() -> Curve {
        Curve::from(HermiteCurve::new(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(3.0, 3.0, 0.0),
            Vec3::new(4.0, 0.0, 0.0),
            Vec3::new(3.0, -3.0, 0.0),
            3,
        ))
    }

    #[test]
    fn tessellation_covers_the_interval() {
        let k = &Kernel::default();
        let c = cubic();
        let ivl = Extent1d::from_min_max(0.0, 1.0);
        let t = c.tessellate_by_bisection(k, &ivl, 0.01, 5.0, 0.0).unwrap();
        assert!(t.parameters.len() > 2);
        assert_eq!(t.parameters.first().copied(), Some(0.0));
        assert!((t.parameters.last().copied().unwrap() - 1.0).abs() < 1e-12);
        assert!(t.parameters.windows(2).all(|w| w[0] < w[1]));
        assert_eq!(t.points.len(), t.parameters.len());
    }

    #[test]
    fn approximation_is_within_tolerance() {
        let k = &Kernel::default();
        let c = cubic();
        let r = c
            .approximate_curve(
                &Kernel::default(),
                &[0.0, 1.0],
                1e-3,
                &ApproximateCurveOptions::default(),
            )
            .unwrap();
        assert!(r.achieved_tol.unwrap() <= 1e-3);
        let p = r.curve.evaluate_point(k, 0.5).unwrap();
        assert!(p.distance(c.evaluate_point(k, 0.5).unwrap()) < 1e-2);
    }

    #[test]
    fn approximation_bounds_outside_domain_refinement_without_rejecting_tolerance() {
        let k = Kernel::default();
        let c = Curve::from(BSplineCurve::create_line_segment(3, Vec3::ZERO, Vec3::X).unwrap());
        let opts = ApproximateCurveOptions {
            create_analytics: false,
            match_parameterization: true,
            just_copy_bsplines: false,
        };
        // A negative tolerance still follows native subdivision and force-pass
        // semantics over the natural interval; it is not replaced or rejected.
        let natural = c
            .approximate_curve_base(&k, &[0.0, 1.0], -1e300, &opts)
            .unwrap();
        assert_eq!(natural.curve.natural_interval(), c.natural_interval());
        assert_eq!(
            natural.curve.evaluate_point(&k, 0.5).unwrap(),
            Vec3::X * 0.5
        );

        // The same stopping rule over an enormous clamped interval would need
        // trillions of spans, despite this curve having only two control points.
        assert!(matches!(
            c.approximate_curve_base(&k, &[0.0, 1e100], -1e300, &opts),
            Err(SolidError::Upstream { code: 1001, .. })
        ));
    }

    #[test]
    fn approximation_fails_when_native_midpoint_overflows() {
        let c = Curve::from(BSplineCurve::create_line_segment(3, Vec3::ZERO, Vec3::X).unwrap());
        let opts = ApproximateCurveOptions {
            create_analytics: false,
            match_parameterization: true,
            just_copy_bsplines: false,
        };
        assert!(matches!(
            c.approximate_curve_base(&Kernel::default(), &[f64::MAX * 0.75, f64::MAX], 0.0, &opts,),
            Err(SolidError::Upstream { code: 1001, .. })
        ));
    }

    #[test]
    fn global_point_solve_finds_closest_point() {
        let k = Kernel::default();
        let c = cubic();
        let ivl = Extent1d::from_min_max(0.0, 1.0);
        let sols = c
            .global_point_solve::<Obj>(
                &k,
                &ivl,
                SolverOperation::Maximize,
                Vec3::new(2.0, -5.0, 0.0),
                1e-6,
                SolutionRequested::Single,
                &GlobalPointSolveOptions::default(),
            )
            .unwrap();
        assert_eq!(sols.len(), 1);
    }

    #[test]
    fn point_solution_range_is_degenerate() {
        let k = &Kernel::default();
        let c = cubic();
        let mut s: Solution<Obj> = Solution::new();
        s.solution_type = SolutionType::SingleValue;
        assert_eq!(
            c.classify_intersecting_range_solution(k, &c, &s).unwrap(),
            (true, None)
        );
    }

    #[derive(Debug, Clone, Copy, PartialEq)]
    struct Obj;
    impl SolutionObject for Obj {
        fn kind(&self) -> crate::geom::solution_array::ObjectKind {
            crate::geom::solution_array::ObjectKind::Other
        }
    }

    #[derive(Debug, Clone, Copy, PartialEq)]
    struct CurveId(usize);
    impl SolutionObject for CurveId {
        fn kind(&self) -> crate::geom::solution_array::ObjectKind {
            crate::geom::solution_array::ObjectKind::Other
        }
    }

    fn range_solution(start: [f64; 2], end: [f64; 2]) -> Solution<CurveId> {
        let mut sol = Solution::new();
        sol.solution_type = SolutionType::RangeOfValues;
        sol.num_variables = 2;
        sol.num_objects = 2;
        sol.objects[0] = Some(CurveId(0));
        sol.objects[1] = Some(CurveId(1));
        sol.start.parameters[..2].copy_from_slice(&start);
        sol.end.parameters[..2].copy_from_slice(&end);
        sol
    }

    #[test]
    fn range_refinement_polishes_both_ends_in_each_parameter_direction() {
        let mut k = Kernel::default();
        k.fixes.refine_cc_range_ends = true;
        // The line and quadratic bow meet exactly at x=0 and x=1.
        for reverse1 in [false, true] {
            for reverse2 in [false, true] {
                let line = Curve::from(HermiteCurve::new(
                    if reverse1 { Vec3::X } else { Vec3::ZERO },
                    if reverse1 { -Vec3::X } else { Vec3::X },
                    if reverse1 { Vec3::ZERO } else { Vec3::X },
                    if reverse1 { -Vec3::X } else { Vec3::X },
                    3,
                ));
                let bow = Curve::from(HermiteCurve::new(
                    if reverse2 { Vec3::X } else { Vec3::ZERO },
                    Vec3::new(if reverse2 { -1.0 } else { 1.0 }, 0.1, 0.0),
                    if reverse2 { Vec3::ZERO } else { Vec3::X },
                    Vec3::new(if reverse2 { -1.0 } else { 1.0 }, -0.1, 0.0),
                    3,
                ));
                let parameter = |x: f64, reverse| if reverse { 1.0 - x } else { x };
                let mut sol = range_solution(
                    [parameter(0.0001, reverse1), parameter(0.0001, reverse2)],
                    [parameter(0.9999, reverse1), parameter(0.9999, reverse2)],
                );
                for end in [&mut sol.start, &mut sol.end] {
                    end.solution_value = line
                        .evaluate_point(&k, end.parameters[0])
                        .unwrap()
                        .distance(bow.evaluate_point(&k, end.parameters[1]).unwrap());
                    assert!(end.solution_value > 1e-6);
                }
                let mut solutions = SolutionArray::new();
                solutions.add(&sol);
                let curves = [&line, &bow];
                refine_cc_solutions(&k, &mut solutions, |id| {
                    id.and_then(|CurveId(i)| curves.get(i).copied())
                })
                .unwrap();
                let refined = solutions.get(0).unwrap();
                for (end, x) in [(refined.start, 0.0), (refined.end, 1.0)] {
                    assert!(
                        end.solution_value < 1e-10,
                        "{reverse1:?}, {reverse2:?}: {end:?}"
                    );
                    assert!((end.parameters[0] - parameter(x, reverse1)).abs() < 1e-8);
                    assert!((end.parameters[1] - parameter(x, reverse2)).abs() < 1e-8);
                }
            }
        }
    }

    #[test]
    fn range_refinement_rejects_an_end_jump_to_the_start() {
        let mut k = Kernel::default();
        k.fixes.refine_cc_range_ends = true;
        let line = Curve::from(HermiteCurve::new(Vec3::ZERO, Vec3::X, Vec3::X, Vec3::X, 3));
        let direction = Vec3::new(1.0, 0.1, 0.0);
        let crossing = Curve::from(HermiteCurve::new(
            Vec3::ZERO,
            direction,
            direction,
            direction,
            3,
        ));
        // The only intersection is at zero. The end solve converges there too,
        // but accepting it would collapse almost the entire candidate range.
        let mut sol = range_solution([0.0001; 2], [0.9; 2]);
        sol.start.solution_value = 0.00001;
        sol.end.solution_value = 0.09;
        let mut solutions = SolutionArray::new();
        solutions.add(&sol);
        let curves = [&line, &crossing];
        refine_cc_solutions(&k, &mut solutions, |id| {
            id.and_then(|CurveId(i)| curves.get(i).copied())
        })
        .unwrap();
        let refined = solutions.get(0).unwrap();
        assert!(refined.start.solution_value < 1e-10);
        assert_eq!(refined.end, sol.end);
    }

    #[test]
    fn span_accuracy_of_exact_hermite_passes() {
        let k = &Kernel::default();
        let c = cubic();
        let a = c.evaluate(k, 0.0, 1, true, false).unwrap();
        let b = c.evaluate(k, 1.0, 1, false, false).unwrap();
        let r =
            test_span_accuracy(k, &c, 0.0, [a[0], a[1]], 1.0, [b[0], b[1]], 1e-9, true).unwrap();
        assert!(r.passes);
        assert!(r.max_dist < 1e-12);
    }
}
