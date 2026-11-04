//! Regression test for tangency-labs/solidcore#2574: the radius test of
//! `SmBSplineCurve::AssertValid` (13 in 2D, 14 in 3D, read with
//! `Kernel::fixes.bspline_check_min_radius`) reported a cusp on a trim that is a straight
//! line, because its `SM_CP_MINIMIZE_RADIUS_OF_CURVATURE` search returned a radius of 0
//! where its objective cannot be evaluated.
//!
//! The curve is the mapped UV trim of a face of `issue_1945_fork_1207_step4`'s result
//! built without `oscillating_skips_failed_guesses` (its `run(false)`), read out of the
//! check: a 2D cubic with 13 control points of weight 1, all at y = 2 sqrt 2, x rising
//! unevenly from 5.6e-15 to 4.70465, over four spans joined `C0` (interior knots of
//! multiplicity 3 = p).
//!
//! Closed form: every control point has the same y, so the curve is the segment of the
//! line y = 2 sqrt 2 from x = 5.6e-15 to x = 4.70465. A straight line has no curvature,
//! so C' x C'' vanishes and the radius of curvature is infinite at every parameter: no
//! point of this trim is a cusp. Measured over 41 parameters, y stays within one unit in
//! the last place of 2 sqrt 2, the speed stays between 33.9571 and 35.2679 (never zero),
//! and |C' x C''| / |C'|^3 stays below 4.7e-15, a radius above 2.1e14 — eleven orders of
//! magnitude above the radius tolerance 3.46062814695851e-3 the check compares against.
//!
//! The search used to return a solution of radius 0 at t = 0.00583189 anyway:
//! `CurvePropertyEfo::evaluate` zeroed `out.f` and then failed where |C' x C''| is zero,
//! which is every point of a straight piece; `LocalSolve1d::eval_into` kept that 0 as the
//! iteration value; and `local_property_analysis` read it back as the solution value of a
//! minimisation that found nothing. The check compared 0 against its radius tolerance and
//! reported the cusp. The three sites now answer the closed form, take no value from a
//! failed evaluation, and report no solution, so the search finds no finite minimum and
//! both checks pass. The last test guards the other side of that: where the search does
//! stop at a value it evaluated, the check must still see it, so a cubic with a tight but
//! regular bend next to a domain end is still reported.
//! Written with Claude Code (Anthropic).
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    missing_docs
)] // test helpers
use solidcore::curve::curve::p2::CurvePropertyType;
use solidcore::curve::curve::p7::PropertyAnalysisOptions;
use solidcore::curve::curve::Curve;
use solidcore::geom::assert_array::{AssertArray, ObjectKey};
use solidcore::geom::solution_array::{ObjectKind, SolutionObject};
use solidcore::kernel::{Kernel, KernelOptions};
use solidcore::nurbs::b_spline_curve::BSplineCurve;
use solidcore::nurbs::math::{CPoint, NOZ};

/// The y every control point sits at (2 sqrt 2).
const Y: f64 = 2.8284271247461916;

/// The control points' x, as the check saw them.
const XS: [f64; 13] = [
    5.560586550410673e-15,
    0.38469966347659956,
    0.7712034706794334,
    1.1595370454427245,
    1.5478706202060155,
    1.9380339617640256,
    2.330053181489324,
    2.7220724012146218,
    3.1159474983184556,
    3.5117050833605332,
    3.9074626684026104,
    4.3051027405713445,
    4.704652421599879,
];

/// The four spans' ends; the inner three carry multiplicity 3, the outer two 4.
const BREAKS: [f64; 5] = [
    0.0,
    0.033986983156039985,
    0.06797396631207998,
    0.10196094946811998,
    0.13594793262415997,
];

/// The solution separation the check passes for this trim (`mapped_size / 20`).
const SEP: f64 = 0.23523262107999368;

/// The radius tolerance the check compares against for this trim.
const RADIUS_TOL: f64 = 3.460628146958512e-3;

/// Parameters sampled by [`straight_line`].
const SAMPLES: usize = 40;

/// The offset of the tight bend's middle control points ([`tight_bend`]).
const BEND_H: f64 = 1e-4;

/// The radius tolerance the check compares [`tight_bend`] against.
const BEND_TOL: f64 = 1.4142135623730951e-4;

/// The smallest radius of curvature of [`tight_bend`], over 20001 parameters.
const BEND_MIN_RADIUS: f64 = 5.815568702310391e-6;

/// The trim as a `dim`-dimensional cubic (`dim == 3` puts it in the plane z = 0).
fn trim(dim: u64) -> BSplineCurve {
    let z = if dim == 2 { NOZ } else { 0.0 };
    let cps: Vec<CPoint> = XS.iter().map(|&x| CPoint::new(x, Y, z, 1.0)).collect();
    let mut knots: Vec<f64> = Vec::with_capacity(17);
    for (i, &b) in BREAKS.iter().enumerate() {
        let mult = if i == 0 || i == BREAKS.len() - 1 {
            4
        } else {
            3
        };
        knots.extend(vec![b; mult]);
    }
    BSplineCurve::from_data(dim, cps.len() - 1, knots.len() - 1, 3, &cps, &knots).unwrap()
}

/// Over [`SAMPLES`] + 1 parameters of `curve`: the largest distance of a point from the
/// line y = [`Y`], the smallest and largest speed, and the largest curvature
/// |C' x C''| / |C'|^3.
fn straight_line(curve: &BSplineCurve) -> (f64, f64, f64, f64) {
    let (mut off, mut slow, mut fast, mut kappa) = (0.0_f64, f64::MAX, 0.0_f64, 0.0_f64);
    for i in 0..=SAMPLES {
        let t = BREAKS[0] + (BREAKS[BREAKS.len() - 1] - BREAKS[0]) * i as f64 / SAMPLES as f64;
        let d = curve.evaluate(t, 3, false, false).unwrap();
        let speed = d[1].length();
        off = off.max((d[0].y - Y).abs());
        slow = slow.min(speed);
        fast = fast.max(speed);
        kappa = kappa.max(d[1].cross(d[2]).length() / (speed * speed * speed));
    }
    (off, slow, fast, kappa)
}

/// The trim's handle in the solutions of its own property analysis.
#[derive(Debug, Clone, Copy, PartialEq)]
struct TrimCurve;

impl SolutionObject for TrimCurve {
    fn kind(&self) -> ObjectKind {
        ObjectKind::Curve
    }
}

/// A kernel whose only opt-in fix is the one the radius test needs: it reads the search's
/// minimum radius (WR-382). The search itself takes no value from a failed evaluation by
/// default since #2574, so `oscillating_skips_failed_guesses` stays off here.
fn check_kernel() -> Kernel {
    let mut k = Kernel::new(&KernelOptions::default());
    k.fixes.bspline_check_min_radius = true;
    k
}

/// The radius test's reports of `curve`'s `AssertValid`: the failed test's index, value
/// and tolerance, with `None` for a clean check.
fn radius_report(k: &Kernel, curve: &BSplineCurve) -> (bool, Option<(u64, f64, f64)>) {
    let mut list = AssertArray::new();
    let ok = curve
        .assert_valid(k, Some(&mut list), false, ObjectKey(1))
        .unwrap();
    let failed = list
        .reports()
        .iter()
        .find(|r| !r.ok)
        .map(|r| (r.test_index, r.value, r.tol3d));
    (ok, failed)
}

/// A cubic Bezier with a tight but regular bend at one end: control points (0,0),
/// (h,0), (h,h), (1,1), reversed for `at_end`. In closed form C'(0) = (3h, 0) and
/// C''(0) = (-6h, 6h), so the speed never vanishes and the radius of curvature at the
/// bend's end is 27h^3 / 18h^2 = 1.5h.
fn tight_bend(h: f64, at_end: bool) -> BSplineCurve {
    let mut pts = [[0.0, 0.0], [h, 0.0], [h, h], [1.0, 1.0]];
    if at_end {
        pts.reverse();
    }
    let cps: Vec<CPoint> = pts
        .iter()
        .map(|q| CPoint::new(q[0], q[1], 0.0, 1.0))
        .collect();
    let knots = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
    BSplineCurve::from_data(3, cps.len() - 1, knots.len() - 1, 3, &cps, &knots).unwrap()
}

/// The smallest radius of curvature |C'|^3 / |C' x C''| of `curve` over 20001
/// parameters, with the parameter it sits at, and the smallest speed.
fn smallest_radius(curve: &BSplineCurve) -> (f64, f64, f64) {
    let (mut best, mut at, mut slow) = (f64::MAX, f64::NAN, f64::MAX);
    for i in 0..=20000 {
        let t = f64::from(i) / 20000.0;
        let d = curve.evaluate(t, 3, false, false).unwrap();
        let speed = d[1].length();
        slow = slow.min(speed);
        let cross = d[1].cross(d[2]).length();
        if cross > 0.0 && speed > 0.0 {
            let r = speed * speed * speed / cross;
            if r < best {
                best = r;
                at = t;
            }
        }
    }
    (best, at, slow)
}

/// **The trim is a straight line, so its radius of curvature is infinite.** Every
/// control point shares one y, so the curve is a segment of the line y = 2 sqrt 2
/// traversed forward at a speed that never vanishes. Its curvature is therefore zero, the
/// radius everywhere is far above the tolerance the check compares against, and this
/// holds independently of the search below.
#[test]
fn the_trim_is_a_straight_line_without_curvature() {
    let (off, slow, fast, kappa) = straight_line(&trim(2));
    assert!(
        off <= 1e-15,
        "the trim leaves the line y = 2 sqrt 2 by {off:e}"
    );
    assert!(slow > 33.9 && fast < 35.3, "speed [{slow},{fast}]");
    assert!(
        kappa < 1e-14,
        "curvature {kappa:e}: the trim is not straight"
    );
    assert!(
        1.0 / kappa > 1e6 * RADIUS_TOL,
        "radius {:e} against the tolerance {RADIUS_TOL:e}",
        1.0 / kappa
    );
    assert!(XS.windows(2).all(|w| w[0] < w[1]), "x is not increasing");
}

/// **The minimum-radius search finds no finite minimum on the straight trim.** Its
/// objective is the radius of curvature, which is unbounded at every parameter of a
/// straight piece, so there is no minimum to report and the search returns no solution.
#[test]
fn the_min_radius_search_finds_no_finite_minimum() {
    let curve = trim(2);
    let k = check_kernel();
    let opts: PropertyAnalysisOptions<'_, TrimCurve> = PropertyAnalysisOptions::default();
    let found = Curve::from(curve.clone())
        .global_property_analysis(
            &k,
            &curve.natural_interval(),
            CurvePropertyType::MinimizeRadiusOfCurvature,
            SEP,
            &opts,
        )
        .unwrap();
    let sols = found.solutions.as_slice();
    assert!(
        sols.is_empty(),
        "{} solutions on a curve without curvature, first value {:e} at t = {:e}",
        sols.len(),
        sols[0].start.solution_value,
        sols[0].start.parameters[0]
    );
}

/// **The check passes on the straight 2D trim (test 13).** `AssertValid` reads the
/// search's first solution as the minimum radius; with no solution there is no radius
/// below the tolerance, and no test fails.
#[test]
fn the_check_passes_on_the_straight_2d_trim() {
    assert_eq!(radius_report(&check_kernel(), &trim(2)), (true, None));
}

/// **The same straight curve in 3D passes test 14.** The fix is in the search, not in the
/// 2D path: the trim put in the plane z = 0 checks clean as well.
#[test]
fn the_check_passes_on_the_straight_3d_curve() {
    assert_eq!(radius_report(&check_kernel(), &trim(3)), (true, None));
}

/// **A real small radius next to a domain end is still reported.** The fix must not turn
/// into a blanket "no solution": where the search does stop at a value it evaluated, that
/// value is the minimum and the check has to see it. This cubic Bezier bends tightly at
/// one end while its speed never vanishes, so in closed form its radius of curvature
/// there is 1.5 h = 1.5e-4, and its smallest radius over the curve is 5.8156e-6 just
/// inside the domain — far below the radius tolerance 1.4142e-4. Test 14 reports it with
/// the bend at either end of the domain.
#[test]
fn the_check_still_reports_a_tight_bend_next_to_a_domain_end() {
    for at_end in [false, true] {
        let curve = tight_bend(BEND_H, at_end);
        let bend_t = if at_end { 1.0 } else { 0.0 };

        // closed form at the bend's end: a radius of 1.5 h at a speed of 3 h
        let d = curve.evaluate(bend_t, 2, false, false).unwrap();
        let radius = d[1].length().powi(3) / d[1].cross(d[2]).length();
        assert!(
            (radius - 1.5 * BEND_H).abs() <= 1e-12 * radius,
            "radius {radius:e} at the bend, not 1.5 h"
        );

        let (min_radius, at, slow) = smallest_radius(&curve);
        assert!(slow > 0.0, "the speed vanishes at some parameter");
        assert!(
            (min_radius - BEND_MIN_RADIUS).abs() <= 1e-9 * BEND_MIN_RADIUS,
            "smallest radius {min_radius:e}"
        );
        assert!(
            min_radius < BEND_TOL,
            "{min_radius:e} is above the tolerance"
        );
        assert!(
            (at - bend_t).abs() < 0.01,
            "the smallest radius sits at t = {at:e}, not next to the bend"
        );

        // the check reports it, and the value it reports is that minimum
        let (ok, failed) = radius_report(&check_kernel(), &curve);
        assert!(!ok, "the tight bend checks clean (at_end = {at_end})");
        let (index, value, tol) = failed.expect("a failed test");
        assert_eq!(index, 14, "failed test {index}, not the 3D radius test");
        assert!(
            (tol - BEND_TOL).abs() <= 1e-15 * tol,
            "radius tolerance {tol:e}"
        );
        assert!(
            value < tol,
            "reported radius {value:e} is not below {tol:e}"
        );
        assert!(
            (value - min_radius).abs() <= 1e-4 * min_radius,
            "reported radius {value:e} is not the minimum {min_radius:e}"
        );
    }
}
