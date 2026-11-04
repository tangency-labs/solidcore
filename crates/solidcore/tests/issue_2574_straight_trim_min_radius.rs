//! Regression test for tangency-labs/solidcore#2574: the radius test of
//! `SmBSplineCurve::AssertValid` (13 in 2D, 14 in 3D, read with
//! `Kernel::fixes.bspline_check_min_radius`) reports a cusp on a trim that is a straight
//! line, because its `SM_CP_MINIMIZE_RADIUS_OF_CURVATURE` search returns a radius of 0
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
//! and |C' x C''| / |C'|^3 stays below 4.7e-15, a radius above 2.1e14.
//!
//! The search returns a solution of radius 0 at t = 0.00583189 anyway:
//! `CurvePropertyEfo::evaluate` zeroes `out.f` and then fails where |C' x C''| is zero,
//! which is every point of a straight piece; `LocalSolve1d::eval_into` keeps that 0 as
//! the iteration value; and `local_property_analysis` reads it back as the solution value
//! of a minimisation that found nothing. The check compares 0 against its radius
//! tolerance 3.46062814695851e-3 and reports the cusp. With `oscillating_skips_failed_guesses`
//! (WR-253, tangency-labs/solidcore#481) the failed evaluations are not taken as values
//! and the same curve checks clean.
//!
//! The three tests named `_today` hold that defect and turn round when #2574 is fixed.
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

/// The parameter the search reports its radius of 0 at.
const CUSP_T: f64 = 5.831894412731113e-3;

/// Parameters sampled by [`straight_line`].
const SAMPLES: usize = 40;

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

/// A kernel whose only opt-in fixes are the two the radius test needs: it reads the
/// search's minimum radius (WR-382), and takes or skips failed evaluations as values.
fn check_kernel(skip_failed_guesses: bool) -> Kernel {
    let mut k = Kernel::new(&KernelOptions::default());
    k.fixes.bspline_check_min_radius = true;
    k.fixes.oscillating_skips_failed_guesses = skip_failed_guesses;
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

/// **The trim is a straight line, so its radius of curvature is infinite.** Every
/// control point shares one y, so the curve is a segment of the line y = 2 sqrt 2
/// traversed forward at a speed that never vanishes. Its curvature is therefore zero,
/// and this holds independently of the defect below.
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
    assert!(XS.windows(2).all(|w| w[0] < w[1]), "x is not increasing");
}

/// **The minimum-radius search returns a radius of 0 on the straight trim.** The search
/// reports one solution, at t = 0.00583189 with the value 0, although the trim has no
/// curvature anywhere.
///
/// Right would be no solution at all (no curvature, so no finite radius to minimise).
/// When the search no longer takes a failed evaluation as a value, this test becomes
/// that check.
#[test]
fn the_min_radius_search_returns_a_zero_radius_today() {
    let curve = trim(2);
    let k = check_kernel(false);
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
    assert_eq!(
        sols.len(),
        1,
        "the search returns {} solutions now: turn this test into the no-solution check",
        sols.len()
    );
    let sol = &sols[0].start;
    assert_eq!(
        sol.solution_value, 0.0,
        "solution value {:e}",
        sol.solution_value
    );
    assert!(
        (sol.parameters[0] - CUSP_T).abs() < 1e-9,
        "radius 0 reported at t = {:e}",
        sol.parameters[0]
    );
}

/// **The check reports the straight 2D trim as a cusp (test 13).** `AssertValid` reads
/// the search's first solution as the minimum radius, compares 0 against its radius
/// tolerance 3.46e-3 and fails.
///
/// Right would be a clean check. When the search no longer returns that 0, this test
/// becomes that check.
#[test]
fn the_check_reports_the_straight_2d_trim_as_a_cusp_today() {
    let (ok, failed) = radius_report(&check_kernel(false), &trim(2));
    assert!(
        !ok,
        "the straight trim checks clean now: turn this test round"
    );
    let (index, value, tol) = failed.expect("a failed test");
    assert_eq!(index, 13, "failed test {index}, not the 2D radius test");
    assert_eq!(value, 0.0, "reported radius {value:e}");
    assert!((tol - RADIUS_TOL).abs() < 1e-15, "radius tolerance {tol:e}");
}

/// **The same straight curve in 3D is reported as a cusp by test 14.** The defect is in
/// the search, not in the 2D path: the trim put in the plane z = 0 fails the 3D radius
/// test with the same value and tolerance.
///
/// Right would be a clean check, as for the 2D trim above.
#[test]
fn the_check_reports_the_straight_3d_curve_as_a_cusp_today() {
    let (ok, failed) = radius_report(&check_kernel(false), &trim(3));
    assert!(
        !ok,
        "the straight curve checks clean now: turn this test round"
    );
    let (index, value, tol) = failed.expect("a failed test");
    assert_eq!(index, 14, "failed test {index}, not the 3D radius test");
    assert_eq!(value, 0.0, "reported radius {value:e}");
    assert!((tol - RADIUS_TOL).abs() < 1e-15, "radius tolerance {tol:e}");
}

/// **Skipping failed guesses leaves the straight trim without a radius.** With
/// `oscillating_skips_failed_guesses` the search finds no solution on the straight trim,
/// and both the 2D and the 3D check pass. This is the result the fix of #2574 has to
/// reach without the flag.
#[test]
fn skipping_failed_guesses_leaves_the_straight_trim_uncurved() {
    let curve = trim(2);
    let k = check_kernel(true);
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
    assert!(
        found.solutions.as_slice().is_empty(),
        "{} solutions on a curve without curvature",
        found.solutions.as_slice().len()
    );
    assert_eq!(radius_report(&k, &curve), (true, None));
    assert_eq!(radius_report(&k, &trim(3)), (true, None));
}
