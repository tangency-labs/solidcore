//! Endpoint candidates must not hide better interior projected extrema (#2840).
//!
//! Independent oracle for these rational quadratic fixtures: construct homogeneous
//! Bernstein polynomials directly, then certify the stationary polynomial's sign
//! using Bernstein coefficients and bisect its unique root. No kernel evaluator,
//! optimizer, or C++ output supplies the oracle. Certificates use ordinary f64
//! arithmetic with large sign margins, not formally outward-rounded intervals.
#![allow(clippy::unwrap_used, clippy::indexing_slicing)]

use solidcore::curve::curve::{p1::GlobalPointSolveOptions, Curve};
use solidcore::geom::extent1d::Extent1d;
use solidcore::geom::global_solver::{SolutionRequested, SolverOperation};
use solidcore::geom::solution_array::{ObjectKind, SolutionObject, SolutionType};
use solidcore::math::Vec3;
use solidcore::nurbs::b_spline_curve::p1::BSplineCurve;
use solidcore::nurbs::crv_basic::Curve as NurbCurve;
use solidcore::nurbs::math::CPoint;
use solidcore::Kernel;

#[derive(Clone, Copy, PartialEq)]
struct Object;
impl SolutionObject for Object {
    fn kind(&self) -> ObjectKind {
        ObjectKind::Other
    }
}

fn eval(p: &[f64], t: f64) -> f64 {
    p.iter()
        .rev()
        .fold(0.0, |value, coefficient| value * t + coefficient)
}

fn choose(n: usize, k: usize) -> f64 {
    (0..k).fold(1.0, |v, i| v * (n - i) as f64 / (i + 1) as f64)
}

fn bernstein(power: &[f64]) -> Vec<f64> {
    let degree = power.len() - 1;
    (0..=degree)
        .map(|i| {
            (0..=i)
                .map(|j| power[j] * choose(i, j) / choose(degree, j))
                .sum()
        })
        .collect()
}

struct Span {
    bounds: [f64; 2],
    homogeneous: [[f64; 3]; 4],
}

impl Span {
    fn new(control: &[[f64; 4]], original: [f64; 2], bounds: [f64; 2]) -> Self {
        assert_eq!(control.len(), 3);
        assert!(control.iter().all(|p| p[3] > 0.0));
        let alpha = (bounds[0] - original[0]) / (original[1] - original[0]);
        let beta = (bounds[1] - bounds[0]) / (original[1] - original[0]);
        let homogeneous = std::array::from_fn(|axis| {
            let a = control[0][axis];
            let b = 2.0 * (control[1][axis] - a);
            let c = a - 2.0 * control[1][axis] + control[2][axis];
            [
                a + b * alpha + c * alpha * alpha,
                beta * (b + 2.0 * c * alpha),
                c * beta * beta,
            ]
        });
        Self {
            bounds,
            homogeneous,
        }
    }

    fn point(&self, s: f64) -> Vec3 {
        let w = eval(&self.homogeneous[3], s);
        Vec3::new(
            eval(&self.homogeneous[0], s) / w,
            eval(&self.homogeneous[1], s) / w,
            eval(&self.homogeneous[2], s) / w,
        )
    }

    fn distance(&self, s: f64, target: Vec3, normal: Vec3) -> f64 {
        normal.cross(target - self.point(s)).cross(normal).length()
    }

    fn stationary(&self, target: Vec3, normal: Vec3) -> [f64; 5] {
        let w = self.homogeneous[3];
        let projected: [Vec3; 3] = std::array::from_fn(|i| {
            let xyz = Vec3::new(
                self.homogeneous[0][i],
                self.homogeneous[1][i],
                self.homogeneous[2][i],
            );
            normal.cross(target * w[i] - xyz).cross(normal)
        });
        let mut numerator = [0.0; 5];
        for i in 0..3 {
            for j in 0..3 {
                numerator[i + j] += projected[i].dot(projected[j]);
            }
        }
        // Distance squared is N/W²; its derivative has numerator H=N'W-2NW'.
        let mut h = [0.0; 6];
        for i in 1..5 {
            for j in 0..3 {
                h[i - 1 + j] += i as f64 * numerator[i] * w[j];
            }
        }
        for i in 0..5 {
            for j in 1..3 {
                h[i + j - 1] -= 2.0 * numerator[i] * j as f64 * w[j];
            }
        }
        let magnitude = h.iter().fold(1.0_f64, |a, b| a.max(b.abs()));
        assert!(
            h[5].abs() <= 1e-12 * magnitude,
            "degree-five terms cancel algebraically"
        );
        [h[0], h[1], h[2], h[3], h[4]]
    }

    fn extremum(&self, target: Vec3, normal: Vec3, maximize: bool, monotone: bool) -> (f64, f64) {
        let h = self.stationary(target, normal);
        let s = if monotone {
            // The first span of gen0002 increases all the way to the join.
            assert!(bernstein(&h).iter().all(|&v| v > 1.0));
            assert!(maximize);
            1.0
        } else {
            let derivative = [h[1], 2.0 * h[2], 3.0 * h[3], 4.0 * h[4]];
            let direction = if maximize { -1.0 } else { 1.0 };
            assert!(
                bernstein(&derivative).iter().all(|&v| direction * v > 0.1),
                "H must be strictly monotone"
            );
            assert!(direction * eval(&h, 0.0) < -0.01 && direction * eval(&h, 1.0) > 0.01);
            let (mut lo, mut hi) = (0.0, 1.0);
            for _ in 0..64 {
                let mid = (lo + hi) / 2.0;
                if direction * eval(&h, mid) < 0.0 {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            (lo + hi) / 2.0
        };
        (
            self.bounds[0] + s * (self.bounds[1] - self.bounds[0]),
            self.distance(s, target, normal),
        )
    }
}

#[allow(clippy::too_many_arguments)]
fn check(
    control: &[[f64; 4]],
    knots: &[f64],
    bounds: [f64; 2],
    target: Vec3,
    normal: Vec3,
    tolerance: f64,
    maximize: bool,
    spans: &[Span],
) {
    let cp = control
        .iter()
        .map(|p| CPoint {
            x: p[0],
            y: p[1],
            z: p[2],
            w: p[3],
        })
        .collect();
    let curve = Curve::from(BSplineCurve::from_nurb(
        NurbCurve::new(cp, 2, knots.to_vec()).unwrap(),
        3,
        false,
    ));
    let candidates: Vec<_> = spans
        .iter()
        .enumerate()
        .map(|(i, s)| s.extremum(target, normal, maximize, spans.len() == 2 && i == 0))
        .collect();
    let expected = candidates
        .iter()
        .copied()
        .reduce(|a, b| {
            if (maximize && b.1 > a.1) || (!maximize && b.1 < a.1) {
                b
            } else {
                a
            }
        })
        .unwrap();
    assert!(expected.0 > bounds[0] && expected.0 < bounds[1]);
    let k = Kernel::default();
    for request in [SolutionRequested::Single, SolutionRequested::All] {
        let answers = curve
            .global_point_solve::<Object>(
                &k,
                &Extent1d::from_min_max(bounds[0], bounds[1]),
                if maximize {
                    SolverOperation::ProjectedMaximize
                } else {
                    SolverOperation::ProjectedMinimize
                },
                target,
                tolerance,
                request,
                &GlobalPointSolveOptions {
                    opt_vectors: Some(&[normal]),
                    ..Default::default()
                },
            )
            .unwrap();
        assert!(!answers.is_empty());
        if request == SolutionRequested::Single {
            assert_eq!(answers.len(), 1);
        }
        for answer in answers.as_slice() {
            assert_eq!(answer.solution_type, SolutionType::SingleValue);
            let t = answer.start.parameters[0];
            let span = spans
                .iter()
                .find(|s| t >= s.bounds[0] && t <= s.bounds[1])
                .unwrap();
            let actual = span.distance(
                (t - span.bounds[0]) / (span.bounds[1] - span.bounds[0]),
                target,
                normal,
            );
            assert!(
                (actual - answer.start.solution_value).abs() <= 1e-12_f64.max(1e-9 * actual.abs()),
                "reported metric disagrees with rational polynomial"
            );
            // Preserve the caller's original tolerance for global optimality.
            assert!((actual-expected.1).abs()<=tolerance,
                "endpoint displaced interior optimum: t={t}, distance={actual}, oracle={expected:?}, tol={tolerance}");
        }
    }
}

// Original generated inputs, retained exactly so the native comparison and the
// independent mathematical regressions exercise the same geometry.

#[test]
fn gen0002_interior_extremum_survives_endpoint_candidates() {
    let control = [
        [
            5.477596627102355,
            -4.556269464099508,
            1.114995725656613,
            1.0,
        ],
        [
            9.076471354595526,
            -2.132022809946599,
            1.2093799310872764,
            0.710088106679831,
        ],
        [11.18257991594583, 3.841957587213183, 4.293909950576297, 1.0],
        [
            6.804762646024233,
            7.58827958804344,
            4.8887288430295435,
            0.710088106679831,
        ],
        [2.251369955507407, 9.248272788188828, 6.340319531047773, 1.0],
        [
            -3.6074205877400116,
            5.545897440403772,
            4.115642140064191,
            0.710088106679831,
        ],
        [
            -3.604555041788431,
            0.9414150425449053,
            3.1959906041602584,
            1.0,
        ],
        [
            -1.5116827423535617,
            -4.208922190082522,
            0.42322769408518396,
            0.710088106679831,
        ],
        [
            5.227686992962869,
            -4.605290016641698,
            1.096440365002791,
            1.0,
        ],
    ];
    let knots = [
        0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
    ];
    let bounds = [0.6256821477105545, 0.9635580665549801];
    let spans = [
        Span::new(&control[4..7], [0.5, 0.75], [bounds[0], 0.75]),
        Span::new(&control[6..9], [0.75, 1.0], [0.75, bounds[1]]),
    ];
    check(
        &control,
        &knots,
        bounds,
        Vec3::new(3.864483291304881, 2.34600166204466, 3.727657628352193),
        Vec3::new(0.6184107265151202, -0.5901631565932342, -0.695866229492025),
        0.05,
        true,
        &spans,
    );
}

#[test]
fn gen2419_interior_extremum_survives_endpoint_candidates() {
    let control = [
        [
            4.68246458913153,
            0.2683108077599766,
            0.7857380281272487,
            1.0,
        ],
        [
            4.405156176361579,
            0.5572137913733559,
            0.660330035392789,
            0.8275822849855569,
        ],
        [5.18329866623959, 1.417855376481916, 0.820266301512733, 1.0],
        [
            4.1740561315767195,
            1.789570193122566,
            0.6973456848123297,
            0.8275822849855569,
        ],
        [
            4.2999559828848914,
            2.308063461512357,
            0.8470050195809443,
            1.0,
        ],
        [
            2.94307866372967,
            2.030654673617574,
            0.7045870141857388,
            0.8275822849855569,
        ],
        [
            3.145828375026399,
            1.8168898669671454,
            0.832251893331091,
            1.0,
        ],
    ];
    let knots = [
        0.0,
        0.0,
        0.0,
        0.3333333333333333,
        0.3333333333333333,
        0.6666666666666666,
        0.6666666666666666,
        1.0,
        1.0,
        1.0,
    ];
    let bounds = [0.18190734190117286, 0.27563504714764536];
    let spans = [Span::new(&control[0..3], [0.0, 1.0 / 3.0], bounds)];
    check(
        &control,
        &knots,
        bounds,
        Vec3::new(7.695548237243017, 1.1786104807110709, 0.9369883409089445),
        Vec3::new(-0.2870092592131892, -0.84405099161452, -0.2736338163951917),
        0.001,
        false,
        &spans,
    );
}
