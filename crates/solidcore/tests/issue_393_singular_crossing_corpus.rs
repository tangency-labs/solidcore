//! Issue #393 (WR-77): the singular-crossing corpus against its closed forms.
//!
//! The same surface pairs as `fixtures/surface.surface_p2/ssi_singular_crossings.json`
//! (`harness/fixturegen/ssi_singular_crossings.py`), rebuilt here: exact rational
//! cylinders, cones, spheres and tori, an exact bilinear saddle and planes, whose
//! intersection branches cross at known nodes. Each pair runs through
//! `GlobalSurfaceIntersect` at the ssi corpus tolerance (1e-6) with the default kernel
//! (upstream behaviour, which difftest holds at parity with the C++ in the `sx_*` cases),
//! and the result is measured against the closed form: the share of each exact curve a
//! returned curve covers, how far the returned curves are from the exact ones, and at
//! each node how many returned curves end there and how many pass through it.
//!
//! #393 is that upstream does not end the curves at a node where the surface normals are
//! parallel (the marcher steps over it): the point set is complete, the topology is not.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::too_many_lines
)] // test helpers

use solidcore::curve::curve::Curve;
use solidcore::geom::tol::approx_tol3d;
use solidcore::geom::vector3d::Vec3;
use solidcore::kernel::Kernel;
use solidcore::nurbs::b_spline_surface::BSplineSurface;
use solidcore::nurbs::math::CPoint;
use solidcore::nurbs::srf_basic::Surface as NurbSurface;
use solidcore::surface::surface::p2::{SurfaceIntersection, TsectCurveType};
use solidcore::surface::surface::Surface;

const W: f64 = std::f64::consts::FRAC_1_SQRT_2;
/// The 9-point unit circle `(x, y, w)` from `+x`, counter-clockwise.
const CIRCLE: [(f64, f64, f64); 9] = [
    (1.0, 0.0, 1.0),
    (1.0, 1.0, W),
    (0.0, 1.0, 1.0),
    (-1.0, 1.0, W),
    (-1.0, 0.0, 1.0),
    (-1.0, -1.0, W),
    (0.0, -1.0, 1.0),
    (1.0, -1.0, W),
    (1.0, 0.0, 1.0),
];
const CIRCLE_KNOTS: [f64; 12] = [
    0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
];
/// The 5-point half circle `(radial, axial, w)` from `-axis` to `+axis`.
const HALF: [(f64, f64, f64); 5] = [
    (0.0, -1.0, 1.0),
    (1.0, -1.0, W),
    (1.0, 0.0, 1.0),
    (1.0, 1.0, W),
    (0.0, 1.0, 1.0),
];
const HALF_KNOTS: [f64; 8] = [0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0];
/// Where each rational circle starts, off the frame axes (degrees).
const TURN_DEG: f64 = 20.0;

const X: Vec3 = Vec3::new(1.0, 0.0, 0.0);
const Y: Vec3 = Vec3::new(0.0, 1.0, 0.0);
const Z: Vec3 = Vec3::new(0.0, 0.0, 1.0);

fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}

fn hom(p: Vec3, w: f64) -> CPoint {
    CPoint::new(p.x * w, p.y * w, p.z * w, w)
}

fn bss(net: Vec<Vec<CPoint>>, p: usize, q: usize, u: Vec<f64>, kv: Vec<f64>) -> Surface {
    let nurb = NurbSurface::new(net, p, q, u, kv).unwrap();
    Surface::from(BSplineSurface::from_nurb(nurb, false))
}

fn clamped(n: usize) -> Vec<f64> {
    let mut k = vec![0.0; n + 1];
    k.extend(vec![1.0; n + 1]);
    k
}

/// The bilinear patch of the plane `origin + u du + v dv`, `[u0, u1] x [v0, v1]`.
fn patch(origin: Vec3, du: Vec3, dv: Vec3, (u0, u1): (f64, f64), (v0, v1): (f64, f64)) -> Surface {
    let net = [u0, u1]
        .iter()
        .map(|&u| {
            [v0, v1]
                .iter()
                .map(|&w| hom(origin + du * u + dv * w, 1.0))
                .collect()
        })
        .collect();
    bss(net, 1, 1, clamped(1), clamped(1))
}

/// `profile` (`(radial, axial, w)`, degree `q`) revolved about `axis` through `center`,
/// the circle starting at `e1` turned by [`TURN_DEG`] towards `e2`.
fn revolved(
    center: Vec3,
    (e1, e2, axis): (Vec3, Vec3, Vec3),
    profile: &[(f64, f64, f64)],
    q: usize,
    knots: &[f64],
) -> Surface {
    let (s, c) = TURN_DEG.to_radians().sin_cos();
    let (e1, e2) = (e1 * c + e2 * s, e1 * -s + e2 * c);
    let net = CIRCLE
        .iter()
        .map(|&(cx, cy, wu)| {
            profile
                .iter()
                .map(|&(rho, z, wv)| {
                    hom(
                        center + e1 * (rho * cx) + e2 * (rho * cy) + axis * z,
                        wu * wv,
                    )
                })
                .collect()
        })
        .collect();
    bss(net, 2, q, CIRCLE_KNOTS.to_vec(), knots.to_vec())
}

/// `center + r (cos e1 + sin e2) + h axis`, `h` in `[0, len]`.
fn cylinder(center: Vec3, frame: (Vec3, Vec3, Vec3), r: f64, len: f64) -> Surface {
    revolved(
        center,
        frame,
        &[(r, 0.0, 1.0), (r, len, 1.0)],
        1,
        &clamped(1),
    )
}

/// The double cone `rho = slope |z|` about the frame's axis through the origin, `z` in
/// `[-1, 1]`: its apex is inside the surface.
fn double_cone(frame: (Vec3, Vec3, Vec3), slope: f64) -> Surface {
    revolved(
        Vec3::ZERO,
        frame,
        &[(-slope, -1.0, 1.0), (slope, 1.0, 1.0)],
        1,
        &clamped(1),
    )
}

fn sphere(frame: (Vec3, Vec3, Vec3), r: f64) -> Surface {
    let prof: Vec<_> = HALF.iter().map(|&(a, b, w)| (r * a, r * b, w)).collect();
    revolved(Vec3::ZERO, frame, &prof, 2, &HALF_KNOTS)
}

fn torus(big: f64, small: f64) -> Surface {
    let prof: Vec<_> = CIRCLE
        .iter()
        .map(|&(a, b, w)| (big + small * a, small * b, w))
        .collect();
    revolved(Vec3::ZERO, (X, Y, Z), &prof, 2, &CIRCLE_KNOTS)
}

/// `z = x y` over `[x0, x1] x [y0, y1]` (bilinear: its corners are its net).
fn saddle((x0, x1): (f64, f64), (y0, y1): (f64, f64)) -> Surface {
    let net = [x0, x1]
        .iter()
        .map(|&x| [y0, y1].iter().map(|&y| hom(v(x, y, x * y), 1.0)).collect())
        .collect();
    bss(net, 1, 1, clamped(1), clamped(1))
}

/// A node of the exact intersection, where branches cross.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum NodeKind {
    /// The surface normals are parallel there (the surfaces touch): #393's case.
    Parallel,
    /// An apex: one or both surfaces have no normal there.
    Apex,
}

/// The exact intersection: polylines (sampled finely) and the nodes.
struct Exact {
    curves: Vec<Vec<Vec3>>,
    nodes: Vec<(Vec3, NodeKind)>,
}

/// `f(t)` for `n + 1` values of `t` evenly over `[a, b]`.
fn sampled(a: f64, b: f64, n: u32, f: impl Fn(f64) -> Vec3) -> Vec<Vec3> {
    (0..=n)
        .map(|i| f(a + (b - a) * f64::from(i) / f64::from(n)))
        .collect()
}

/// One corpus pair (the fixture's `sx_<name>`).
struct Pair {
    name: &'static str,
    a: Surface,
    b: Surface,
    exact: Exact,
}

fn corpus() -> Vec<Pair> {
    use std::f64::consts::{PI, TAU};
    let s3 = 3.0_f64.sqrt();
    let b60 = v(s3 / 2.0, 0.0, 0.5);
    let b60_perp = v(-0.5, 0.0, s3 / 2.0);
    let d30 = v(0.0, s3 / 2.0, 0.5);
    let (s20, c20) = 20.0_f64.to_radians().sin_cos();
    let poles = v(s20, c20, 0.0);
    // the curves of cylinder A (axis z, unit) cut by the plane z = m x
    let on_cyl_z = |m: f64| sampled(0.0, TAU, 1440, move |t| v(t.cos(), t.sin(), m * t.cos()));
    let line = |d: Vec3, t: f64| sampled(-t, t, 400, move |s| d * s);
    vec![
        Pair {
            name: "steinmetz_perpendicular",
            a: cylinder(v(0.0, 0.0, -2.0), (X, Y, Z), 1.0, 4.0),
            b: cylinder(v(-2.0, 0.0, 0.0), (Y, Z, X), 1.0, 4.0),
            exact: Exact {
                curves: vec![on_cyl_z(1.0), on_cyl_z(-1.0)],
                nodes: vec![
                    (v(0.0, 1.0, 0.0), NodeKind::Parallel),
                    (v(0.0, -1.0, 0.0), NodeKind::Parallel),
                ],
            },
        },
        Pair {
            name: "steinmetz_oblique",
            a: cylinder(v(0.0, 0.0, -2.0), (X, Y, Z), 1.0, 4.0),
            b: cylinder(b60 * -2.0, (Y, b60_perp, b60), 1.0, 4.0),
            exact: Exact {
                // p.(a - b) = 0: z = x sin60 / (1 - cos60); p.(a + b) = 0: z = -x sin60 / (1 + cos60)
                curves: vec![on_cyl_z(s3), on_cyl_z(-s3 / 3.0)],
                nodes: vec![
                    (v(0.0, 1.0, 0.0), NodeKind::Parallel),
                    (v(0.0, -1.0, 0.0), NodeKind::Parallel),
                ],
            },
        },
        Pair {
            name: "torus_bitangent_plane",
            a: torus(2.0, 1.0),
            b: patch(Vec3::ZERO, X, d30, (-3.5, 3.5), (-3.5, 3.5)),
            exact: Exact {
                curves: [1.0, -1.0]
                    .iter()
                    .map(|&cx| {
                        sampled(0.0, TAU, 1440, move |t| {
                            v(cx, 0.0, 0.0) + X * (2.0 * t.cos()) + d30 * (2.0 * t.sin())
                        })
                    })
                    .collect(),
                nodes: vec![
                    (v(0.0, 1.5, s3 / 2.0), NodeKind::Parallel),
                    (v(0.0, -1.5, -s3 / 2.0), NodeKind::Parallel),
                ],
            },
        },
        Pair {
            name: "saddle_xy_plane",
            a: saddle((-1.0, 1.5), (-0.7, 1.0)),
            b: patch(Vec3::ZERO, X, Y, (-2.0, 2.0), (-2.0, 2.0)),
            exact: Exact {
                curves: vec![
                    sampled(-1.0, 1.5, 400, |x| v(x, 0.0, 0.0)),
                    sampled(-0.7, 1.0, 400, |y| v(0.0, y, 0.0)),
                ],
                nodes: vec![(Vec3::ZERO, NodeKind::Parallel)],
            },
        },
        Pair {
            name: "viviani_sphere_cylinder",
            a: sphere((Z, v(c20, -s20, 0.0), poles), 2.0),
            b: cylinder(v(1.0, 0.0, -2.5), (X, Y, Z), 1.0, 5.0),
            exact: Exact {
                curves: vec![sampled(0.0, 4.0 * PI, 2880, |t| {
                    v(1.0 + t.cos(), t.sin(), 2.0 * (t / 2.0).sin())
                })],
                nodes: vec![(v(2.0, 0.0, 0.0), NodeKind::Parallel)],
            },
        },
        Pair {
            name: "cone_plane_apex",
            a: double_cone((X, Y, Z), 1.0),
            b: patch(Vec3::ZERO, Y, v(0.5, 0.0, 1.0), (-1.5, 1.5), (-1.5, 1.5)),
            exact: Exact {
                curves: vec![
                    line(v(0.5, s3 / 2.0, 1.0), 1.0),
                    line(v(0.5, -s3 / 2.0, 1.0), 1.0),
                ],
                nodes: vec![(Vec3::ZERO, NodeKind::Apex)],
            },
        },
        Pair {
            name: "cones_shared_apex",
            a: double_cone((X, Y, Z), 1.0),
            b: double_cone((Y, Z, X), s3),
            exact: Exact {
                curves: [(1.0, 1.0), (1.0, -1.0), (-1.0, 1.0), (-1.0, -1.0)]
                    .iter()
                    .map(|&(sy, sz)| line(v(1.0, sy, sz * 2.0_f64.sqrt()), W))
                    .collect(),
                nodes: vec![(Vec3::ZERO, NodeKind::Apex)],
            },
        },
    ]
}

/// How far a sample may be from the exact curve (or the result) and still count.
const NEAR: f64 = 1e-3;
/// How close a sample of a returned curve must come to a node for the minimum there to
/// be refined (the curves are sampled at 4001 parameters).
const CANDIDATE: f64 = 1e-2;
/// How close a returned curve must come to a node to reach it (the corpus tolerance is
/// 1e-6).
const AT_NODE: f64 = 1e-5;

fn seg_dist(p: Vec3, a: Vec3, b: Vec3) -> f64 {
    let d = b - a;
    let l2 = d.dot(d);
    let t = if l2 > 0.0 {
        ((p - a).dot(d) / l2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    p.distance(a + d * t)
}

fn poly_dist(p: Vec3, poly: &[Vec3]) -> f64 {
    poly.windows(2)
        .map(|w| seg_dist(p, w[0], w[1]))
        .fold(f64::INFINITY, f64::min)
}

/// What the returned curves do at one node.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
struct AtNode {
    /// Curve ends within [`AT_NODE`] of the node.
    ends: usize,
    /// Passes of a curve through the node (within [`AT_NODE`]) away from its ends.
    through: usize,
}

/// What a result is, measured against the exact intersection.
#[derive(Debug)]
struct Measure {
    /// The status error, if the intersection failed.
    error: Option<String>,
    /// The curve types, in order.
    types: Vec<TsectCurveType>,
    /// Per exact curve: the share of its samples a returned non-point curve covers.
    covered: Vec<f64>,
    /// The largest distance from a sample of a returned curve to the exact curves.
    off: f64,
    /// Per node: what the returned curves do there.
    nodes: Vec<AtNode>,
    /// The largest distance between a returned 3D curve and either surface evaluated on
    /// its UV curve, at the same parameter (the three curves share one parameterisation).
    uv_off: f64,
}

/// The distance from `node` to `cv` at the local minimum bracketed by `[a, b]`
/// (golden section).
fn refined_min(k: &Kernel, cv: &Curve, node: Vec3, mut a: f64, mut b: f64) -> f64 {
    let g = 0.5 * (5.0_f64.sqrt() - 1.0);
    let f = |t: f64| cv.evaluate_point(k, t).unwrap().distance(node);
    let (mut c, mut d) = (b - g * (b - a), a + g * (b - a));
    let (mut fc, mut fd) = (f(c), f(d));
    for _ in 0..100 {
        if fc < fd {
            (b, d, fd) = (d, c, fc);
            c = b - g * (b - a);
            fc = f(c);
        } else {
            (a, c, fc) = (c, d, fd);
            d = a + g * (b - a);
            fd = f(d);
        }
    }
    f(0.5 * (a + b))
}

/// The ends of `cv` at `node` and its passes through `node` away from its ends.
fn at_node(k: &Kernel, cv: &Curve, node: Vec3) -> AtNode {
    let ivl = cv.natural_interval();
    let n = 4000;
    let ts: Vec<f64> = (0..=n)
        .map(|j| ivl.evaluate(f64::from(j) / f64::from(n)))
        .collect();
    let ds: Vec<f64> = ts
        .iter()
        .map(|&t| cv.evaluate_point(k, t).unwrap().distance(node))
        .collect();
    let mut r = AtNode {
        ends: usize::from(ds[0] < AT_NODE) + usize::from(ds[n as usize] < AT_NODE),
        through: 0,
    };
    // interior local minima, away from the ends (a hundredth of the curve)
    let margin = (n / 100) as usize;
    for i in margin.max(1)..(n as usize - margin.max(1)) {
        if ds[i] < CANDIDATE
            && ds[i] <= ds[i - 1]
            && ds[i] < ds[i + 1]
            && refined_min(k, cv, node, ts[i - 1], ts[i + 1]) < AT_NODE
        {
            r.through += 1;
        }
    }
    r
}

fn measure(
    r: &Result<SurfaceIntersection, solidcore::SolidError>,
    surfaces: [&Surface; 2],
    e: &Exact,
) -> Measure {
    let k = Kernel::default();
    let x = match r {
        Ok(x) => x,
        Err(err) => {
            return Measure {
                error: Some(err.to_string()),
                types: vec![],
                covered: vec![0.0; e.curves.len()],
                off: 0.0,
                nodes: vec![AtNode::default(); e.nodes.len()],
                uv_off: 0.0,
            }
        }
    };
    let mut polys: Vec<Vec<Vec3>> = Vec::new();
    let mut nodes = vec![AtNode::default(); e.nodes.len()];
    let mut uv_off: f64 = 0.0;
    for (i, cv) in x.curves_3d.iter().enumerate() {
        let ivl = cv.natural_interval();
        for (s, uv) in surfaces
            .iter()
            .zip([&x.surface1_uv_curves[i], &x.surface2_uv_curves[i]])
        {
            let uv = uv.as_ref().unwrap();
            for j in 0..=100 {
                let t = ivl.evaluate(f64::from(j) / 100.0);
                let q = uv.evaluate_point(&k, t).unwrap();
                let on = s
                    .evaluate_point(&k, solidcore::Point2::new(q.x, q.y))
                    .unwrap();
                uv_off = uv_off.max(on.distance(cv.evaluate_point(&k, t).unwrap()));
            }
        }
    }
    for (cv, t) in x.curves_3d.iter().zip(&x.curve_types) {
        let ivl = cv.natural_interval();
        polys.push(
            (0..=400)
                .map(|j| {
                    cv.evaluate_point(&k, ivl.evaluate(f64::from(j) / 400.0))
                        .unwrap()
                })
                .collect(),
        );
        if *t != TsectCurveType::Touching {
            for (acc, &(node, _)) in nodes.iter_mut().zip(&e.nodes) {
                let a = at_node(&k, cv, node);
                acc.ends += a.ends;
                acc.through += a.through;
            }
        }
    }
    let off = polys
        .iter()
        .flatten()
        .map(|&p| {
            e.curves
                .iter()
                .map(|c| poly_dist(p, c))
                .fold(f64::INFINITY, f64::min)
        })
        .fold(0.0, f64::max);
    let covered = e
        .curves
        .iter()
        .map(|c| {
            let hit = c
                .iter()
                .filter(|&&p| polys.iter().any(|poly| poly_dist(p, poly) < NEAR))
                .count();
            hit as f64 / c.len() as f64
        })
        .collect();
    Measure {
        error: None,
        types: x.curve_types.clone(),
        covered,
        off,
        nodes,
        uv_off,
    }
}

/// `GlobalSurfaceIntersect` over the natural domains at the corpus tolerance (1e-6),
/// with the default kernel (upstream behaviour).
fn intersect(a: &Surface, b: &Surface) -> Result<SurfaceIntersection, solidcore::SolidError> {
    let k = Kernel::default();
    let d1 = a.natural_uv_domain_with_nurb(&k)?;
    let d2 = b.natural_uv_domain_with_nurb(&k)?;
    let tol = approx_tol3d(&k, 1e-6);
    a.global_surface_intersect(&k, &d1, b, &d2, [true, true], Some(tol), None)
}

fn run(p: &Pair, swapped: bool) -> Measure {
    let (a, b) = if swapped { (&p.b, &p.a) } else { (&p.a, &p.b) };
    measure(&intersect(a, b), [a, b], &p.exact)
}

/// Prints the whole corpus, both operand orders (`--nocapture`).
#[test]
#[ignore = "survey: prints the corpus table"]
fn issue_393_survey() {
    for p in corpus() {
        for swapped in [false, true] {
            let m = run(&p, swapped);
            eprintln!(
                "{:<26} {:<8} err={:?} n={} covered={:?} off={:.1e} uv_off={:.1e} nodes={:?} types={:?}",
                p.name,
                if swapped { "swapped" } else { "" },
                m.error,
                m.types.len(),
                m.covered
                    .iter()
                    .map(|c| (c * 1000.0).round() / 10.0)
                    .collect::<Vec<_>>(),
                m.off,
                m.uv_off,
                m.nodes
                    .iter()
                    .map(|n| (n.ends, n.through))
                    .collect::<Vec<_>>(),
                m.types
            );
        }
    }
}

/// Upstream finds the whole point set of every pair, in both operand orders: each exact
/// curve is covered end to end, and every returned curve lies on the exact curves (within
/// the sampling of the closed forms).
#[test]
fn issue_393_upstream_point_set_complete() {
    for p in corpus() {
        for swapped in [false, true] {
            let m = run(&p, swapped);
            assert!(m.error.is_none(), "{} {swapped}: {m:?}", p.name);
            assert!(
                m.covered.iter().all(|&c| c >= 0.995),
                "{} {swapped}: {m:?}",
                p.name
            );
            assert!(m.off < 1e-5, "{} {swapped}: {m:?}", p.name);
            assert!(m.uv_off < 1e-5, "{} {swapped}: {m:?}", p.name);
            assert!(
                m.types.iter().all(|&t| t == TsectCurveType::Crossing),
                "{} {swapped}: {m:?}",
                p.name
            );
        }
    }
}

/// What upstream (the default kernel, at parity with the C++ in difftest's `sx_*` cases)
/// does at each node, as `(curve ends there, curves passing through)`:
///
/// - where the normals are parallel (Steinmetz, Villarceau, Viviani), two curves pass
///   through the node and none ends there: #393. The saddle on a plane is split at its
///   node in one operand order and passed through in the other.
/// - at a cone's apex (no normal) a plane's two lines end there (four ends); two cones
///   sharing their apex run all four lines through it.
#[test]
fn issue_393_upstream_runs_curves_through_the_nodes() {
    let through = AtNode {
        ends: 0,
        through: 2,
    };
    let split = AtNode {
        ends: 4,
        through: 0,
    };
    let expected: [(&str, [Vec<AtNode>; 2]); 7] = [
        (
            "steinmetz_perpendicular",
            [vec![through, through], vec![through, through]],
        ),
        (
            "steinmetz_oblique",
            [vec![through, through], vec![through, through]],
        ),
        (
            "torus_bitangent_plane",
            [vec![through, through], vec![through, through]],
        ),
        ("saddle_xy_plane", [vec![split], vec![through]]),
        ("viviani_sphere_cylinder", [vec![through], vec![through]]),
        ("cone_plane_apex", [vec![split], vec![split]]),
        (
            "cones_shared_apex",
            [
                vec![AtNode {
                    ends: 0,
                    through: 4,
                }],
                vec![AtNode {
                    ends: 0,
                    through: 4,
                }],
            ],
        ),
    ];
    let pairs = corpus();
    assert_eq!(pairs.len(), expected.len());
    for (name, [fwd, swp]) in expected {
        let p = pairs.iter().find(|p| p.name == name).unwrap();
        assert_eq!(run(p, false).nodes, fwd, "{name}");
        assert_eq!(run(p, true).nodes, swp, "{name} swapped");
    }
}

/// An implicit equation, zero on its surface.
type Implicit = fn(Vec3) -> f64;

/// The closed forms the fixture's surfaces are built on: each surface's control net
/// reproduces its implicit equation at sampled parameters.
#[test]
fn issue_393_surfaces_are_the_closed_forms() {
    let k = Kernel::default();
    let cyl_z: Implicit = |p| p.x * p.x + p.y * p.y - 1.0;
    let checks: [(&str, Implicit, Implicit); 7] = [
        ("steinmetz_perpendicular", cyl_z, |p| {
            p.y * p.y + p.z * p.z - 1.0
        }),
        ("steinmetz_oblique", cyl_z, |p| {
            let b = v(3.0_f64.sqrt() / 2.0, 0.0, 0.5);
            p.dot(p) - p.dot(b) * p.dot(b) - 1.0
        }),
        (
            "torus_bitangent_plane",
            |p| {
                let q = (p.x * p.x + p.y * p.y).sqrt() - 2.0;
                q * q + p.z * p.z - 1.0
            },
            |p| p.z * 3.0_f64.sqrt() - p.y,
        ),
        ("saddle_xy_plane", |p| p.z - p.x * p.y, |p| p.z),
        (
            "viviani_sphere_cylinder",
            |p| p.dot(p) - 4.0,
            |p| (p.x - 1.0) * (p.x - 1.0) + p.y * p.y - 1.0,
        ),
        (
            "cone_plane_apex",
            |p| p.x * p.x + p.y * p.y - p.z * p.z,
            |p| p.x - 0.5 * p.z,
        ),
        (
            "cones_shared_apex",
            |p| p.x * p.x + p.y * p.y - p.z * p.z,
            |p| p.y * p.y + p.z * p.z - 3.0 * p.x * p.x,
        ),
    ];
    for p in corpus() {
        let (_, ea, eb) = checks.iter().find(|c| c.0 == p.name).unwrap();
        for (s, e) in [(&p.a, ea), (&p.b, eb)] {
            let d = s.natural_uv_domain_with_nurb(&k).unwrap();
            for i in 0..=10 {
                for j in 0..=10 {
                    let uv = d.evaluate(f64::from(i) / 10.0, f64::from(j) / 10.0);
                    let x = s.evaluate_point(&k, uv).unwrap();
                    assert!(e(x).abs() < 1e-12, "{} {x:?}", p.name);
                }
            }
        }
    }
}
