//! Red regressions for the original All requests in #3060.
//! Compile with rustc --test, linking solidcore and serde_json from the same build.
//! These are diagnostic tests, not a claim that All certification is implemented.

use serde_json::Value;
use solidcore::curve::curve::{p9::CurveSolveOptions, Curve};
use solidcore::curve::ellipse::{Ellipse, EllipseOptions};
use solidcore::curve::line::{BoundedLineOptions, Line};
use solidcore::geom::extent1d::Extent1d;
use solidcore::geom::global_solver::{SolutionRequested, SolverOperation};
use solidcore::geom::solution_array::{ObjectKind, SolutionObject, SolutionType};
use solidcore::kernel::Kernel;
use solidcore::math::Vec3;
use solidcore::nurbs::b_spline_curve::BSplineCurve;
use solidcore::nurbs::nurbs_crv::allocate_nurb_curve;

#[derive(Clone, Copy, PartialEq, Debug)]
struct Handle;
impl SolutionObject for Handle {
    fn kind(&self) -> ObjectKind {
        ObjectKind::Curve
    }
}

fn number(v: &Value) -> f64 {
    v.as_f64().unwrap()
}
fn vector(v: &Value) -> Vec3 {
    Vec3::new(number(&v[0]), number(&v[1]), number(&v[2]))
}
fn interval(v: &Value) -> Extent1d {
    Extent1d::from_min_max(number(&v[0]), number(&v[1]))
}
fn curve(k: &Kernel, v: &Value) -> Curve {
    if let Some(v) = v.get("line") {
        return Curve::from(Line::bounded(
            vector(&v["p"]),
            vector(&v["v"]),
            &interval(&v["domain"]),
            &BoundedLineOptions {
                scale: number(&v["scale"]),
                inside_out: v["inside_out"].as_bool().unwrap(),
                ..Default::default()
            },
        ));
    }
    if let Some(v) = v.get("ellipse") {
        return Curve::from(
            Ellipse::new(
                k,
                vector(&v["c"]),
                vector(&v["x"]),
                vector(&v["y"]),
                interval(&v["domain"]),
                number(&v["radii"][0]),
                number(&v["radii"][1]),
                &EllipseOptions::default(),
            )
            .unwrap(),
        );
    }
    let v = &v["bspline"];
    let points = v["cpts"].as_array().unwrap();
    let knots = v["u"].as_array().unwrap();
    let mut net = allocate_nurb_curve(
        points.len() - 1,
        v["p"].as_u64().unwrap() as usize,
        knots.len() - 1,
    )
    .unwrap();
    for (dst, src) in net.cpts_mut().iter_mut().zip(points) {
        dst.x = number(&src[0]);
        dst.y = number(&src[1]);
        dst.z = number(&src[2]);
        dst.w = number(&src[3]);
    }
    for (dst, src) in net.knots_mut().iter_mut().zip(knots) {
        *dst = number(src);
    }
    let mut spline = BSplineCurve::from_nurb(net, v["dim"].as_u64().unwrap(), false);
    spline.set_out_of_bounds_enabled(v["oob"].as_bool().unwrap());
    Curve::from(spline)
}

fn check(name: &str, witness: Option<[f64; 2]>, required_corners: &[[f64; 2]]) {
    let cases: Value =
        serde_json::from_str(include_str!("projected_pair_all_cases_3060.json")).unwrap();
    let input = &cases
        .as_array()
        .unwrap()
        .iter()
        .find(|case| case["name"] == name)
        .unwrap()["inputs"];
    assert_eq!(input["requested"], 1);
    let k = Kernel::default();
    let a = curve(&k, &input["c1"]);
    let own_b = (!input["c2"].is_null()).then(|| curve(&k, &input["c2"]));
    let b = own_b.as_ref().unwrap_or(&a);
    let ia = interval(&input["ivl1"]);
    let ib = interval(&input["ivl2"]);
    let tol = number(&input["tol"]);
    let n = vector(&input["vecs"][0]);
    let witness = witness.unwrap_or_else(|| {
        // Ellipse fixture's input domain is STEP degrees, but its solver
        // interval and returned parameters are NURBS coordinates.
        let v = &input["c1"]["ellipse"];
        let e = Ellipse::new(
            &k,
            vector(&v["c"]),
            vector(&v["x"]),
            vector(&v["y"]),
            interval(&v["domain"]),
            number(&v["radii"][0]),
            number(&v["radii"][1]),
            &EllipseOptions::default(),
        )
        .unwrap();
        [
            e.convert_t_from_step_to_nurbs(&k, 0.0).unwrap(),
            e.convert_t_from_step_to_nurbs(&k, 90.0).unwrap(),
        ]
    });
    for (t, requested, c) in [(witness[0], ia, &a), (witness[1], ib, b)] {
        assert!(requested.contains_value(t, 0.0));
        assert!(c.natural_interval().contains_value(t, 0.0));
    }
    let measure = |s, t| {
        n.cross(a.evaluate_point(&k, s).unwrap() - b.evaluate_point(&k, t).unwrap())
            .cross(n)
            .length()
    };
    let lower_bound = measure(witness[0], witness[1]);
    let options = CurveSolveOptions::<Handle> {
        opt_vectors: Some(&[n]),
        ..Default::default()
    };
    let answers = a
        .global_curve_solve(
            &k,
            &ia,
            b,
            &ib,
            SolverOperation::ProjectedMaximize,
            tol,
            &options,
            SolutionRequested::All,
        )
        .unwrap();
    assert!(!answers.is_empty(), "{name}: no answers");
    for answer in answers.as_slice() {
        assert_eq!(answer.solution_type, SolutionType::SingleValue);
        let [s, t] = [answer.start.parameters[0], answer.start.parameters[1]];
        assert!(ia.contains_value(s, 0.0) && ib.contains_value(t, 0.0));
        assert!(a.natural_interval().contains_value(s, 0.0));
        assert!(b.natural_interval().contains_value(t, 0.0));
        let reported = answer.start.solution_value;
        assert!((reported - measure(s, t)).abs() <= tol);
        assert!(
            reported + tol >= lower_bound,
            "{name}: reported={reported:.17e}, admissible={lower_bound:.17e}, tolerance={tol}"
        );
    }
    // These affine cases have distinct constrained corner maxima within the
    // caller's tolerance of the best; All must not collapse them to Single.
    let mut matched = std::collections::BTreeSet::new();
    for wanted in required_corners {
        // Bound the objective perturbation by the sum of the projected
        // endpoint displacements (triangle inequality), in caller units.
        // Matching is injective so one answer cannot stand for two maxima.
        let found = answers.as_slice().iter().enumerate().find(|(i, answer)| {
            if matched.contains(i) {
                return false;
            }
            let displacement = |c: &Curve, actual, expected| {
                n.cross(
                    c.evaluate_point(&k, actual).unwrap() - c.evaluate_point(&k, expected).unwrap(),
                )
                .cross(n)
                .length()
            };
            displacement(&a, answer.start.parameters[0], wanted[0])
                + displacement(b, answer.start.parameters[1], wanted[1])
                <= tol
        });
        let (index, _) =
            found.unwrap_or_else(|| panic!("{name}: missing distinct maximum near {wanted:?}"));
        matched.insert(index);
    }
}

#[test]
fn gen0472_ellipse_all_must_dominate_a_feasible_converted_witness() {
    check("gen0472", None, &[]);
}
#[test]
fn gen1839_rational_all_must_dominate_the_opposite_corner() {
    check(
        "gen1839",
        Some([0.43561427937496744, 0.03474215986040488]),
        &[],
    );
}
#[test]
fn gen2050_all_retains_both_tolerance_equivalent_corner_maxima() {
    let lo = 0.30573599521598394;
    let hi = 0.42689846142600174;
    check("gen2050", Some([hi, lo]), &[[hi, lo], [lo, hi]]);
}
#[test]
fn gen2187_all_must_dominate_the_larger_corner() {
    check(
        "gen2187",
        Some([3.254294670614625, -3.304387257364806]),
        &[],
    );
}
#[test]
fn gen2585_all_retains_both_tolerance_equivalent_corner_maxima() {
    check(
        "gen2585",
        Some([0.6795950883084139, 1.0610357492278748]),
        &[
            [0.6795950883084139, 1.0610357492278748],
            [0.15819072816225468, 3.977500125051104],
        ],
    );
}
