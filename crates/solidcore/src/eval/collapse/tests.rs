use super::*;
use crate::eval::uv::{degenerate, degenerate_within};

const MID: f64 = std::f64::consts::FRAC_1_SQRT_2;

/// The spline of `degree` whose net is `rows[i][j]` (i along u), Euclidean points with weights,
/// on clamped knots with one span along each direction.
fn spline(degree: [usize; 2], rows: &[Vec<([f64; 3], f64)>]) -> EvalSurface {
    let (nu, nv) = (rows.len(), rows[0].len());
    let clamped = |p: usize, n: usize| {
        let mut k = vec![0.0; p + 1];
        let inner = n - p - 1;
        k.extend((1..=inner).map(|i| i as f64 / (inner + 1) as f64));
        k.extend(vec![1.0; p + 1]);
        k
    };
    let mut ctrl = Vec::new();
    for j in 0..nv {
        for row in rows {
            let ([x, y, z], w) = row[j];
            ctrl.push([w * x, w * y, w * z, w]);
        }
    }
    EvalSurface::nurbs(
        degree,
        [clamped(degree[0], nu), clamped(degree[1], nv)],
        [nu, nv],
        ctrl,
    )
}

/// The quarter turn about z of the profile `profile` ((radius, height) control points with
/// weights, degree one less than their count), scaled by `scale` and moved by `offset`, as the
/// record holds a revolution: u the profile, v the quarter arc's rational parameter.
fn turned(profile: &[([f64; 2], f64)], scale: f64, offset: [f64; 3]) -> EvalSurface {
    let arc = [([1.0, 0.0], 1.0), ([1.0, 1.0], MID), ([0.0, 1.0], 1.0)];
    let rows: Vec<Vec<([f64; 3], f64)>> = profile
        .iter()
        .map(|([r, z], wp)| {
            arc.iter()
                .map(|([ax, ay], wa)| {
                    let p = [
                        offset[0] + scale * r * ax,
                        offset[1] + scale * r * ay,
                        offset[2] + scale * z,
                    ];
                    (p, wp * wa)
                })
                .collect()
        })
        .collect();
    spline([profile.len() - 1, 2], &rows)
}

/// The cone of the comparison's face 3: the straight profile (s², s²) of a quadratic Bézier
/// with its first two control points at the apex, so `S_u = S_v = 0` there.
fn cone(scale: f64, offset: [f64; 3]) -> EvalSurface {
    turned(
        &[([0.0, 0.0], 1.0), ([0.0, 0.0], 1.0), ([1.0, 1.0], 1.0)],
        scale,
        offset,
    )
}

/// The apex of a cone whose first derivatives both vanish there is collapsed along its row (v,
/// the arc) and not along u, wherever the cone lies and at any scale, at every v, exactly and
/// within a reach; off the apex nothing is collapsed. The speeds alone read both directions as
/// collapsed there, and an apex off the origin (its derivatives then rounding, not zero) fell
/// back to that.
#[test]
fn an_apex_of_zero_speed_is_collapsed_along_its_row_anywhere_and_at_any_scale() {
    let cases = [
        (10.0, [0.0, 0.0, 0.0]),
        (10.0, [10.0, 0.0, 0.0]),
        (10.0, [1000.0, 1000.0, 1000.0]),
        (1e-2, [0.01, 0.0, 0.0]),
        (1e4, [1e4, -2e4, 3e4]),
    ];
    for (scale, offset) in cases {
        let s = cone(scale, offset);
        let reach = 1e-7 * scale;
        for v in [0.0, 0.3, 0.77, 1.0] {
            let (_, su, sv) = s.derivs(0.0, v).unwrap();
            let size = scale + offset.iter().map(|x| x.abs()).fold(0.0, f64::max);
            assert!(
                su.len() <= 1e-12 * size && sv.len() <= 1e-12 * size,
                "{scale} {offset:?} {v}"
            );
            for r in [0.0, reach] {
                assert_eq!(
                    collapse(&s, [0.0, v], r).unwrap(),
                    Collapse::V,
                    "{scale} {offset:?} v {v} reach {r}"
                );
            }
            assert_eq!(collapse(&s, [0.5, v], 0.0).unwrap(), Collapse::Neither);
            assert_eq!(collapse(&s, [0.5, v], reach).unwrap(), Collapse::Neither);
            // Beside the apex, within the reach: its circle (radius scale · u²) is within it,
            // the meridian is not. A speed × span test read both as collapsed here.
            assert_eq!(collapse(&s, [1e-9, v], reach).unwrap(), Collapse::V);
            assert!(degenerate(&s, [0.0, v], 1).unwrap() && !degenerate(&s, [0.0, v], 0).unwrap());
            assert!(degenerate_within(&s, [0.0, v], 1, reach).unwrap());
            assert!(!degenerate_within(&s, [0.0, v], 0, reach).unwrap());
        }
    }
    // A rational profile (weights 1, 0.5, 1) off the origin.
    let s = turned(
        &[([0.0, 0.0], 1.0), ([0.0, 0.0], 0.5), ([1.0, 1.0], 1.0)],
        10.0,
        [10.0, 0.0, 0.0],
    );
    for v in [0.0, 0.3, 1.0] {
        assert_eq!(collapse(&s, [0.0, v], 1e-6).unwrap(), Collapse::V);
    }
}

/// An apex of degree 5 (S ~ s⁵, every derivative to the fourth zero along both directions) is
/// collapsed along its row; a row of degree 5 whose first five points coincide but whose sixth
/// does not is not collapsed, though its derivatives to the fourth vanish at v = 0.
#[test]
fn degree_five_apexes_and_rows_are_decided_on_the_net() {
    let mut profile = vec![([0.0, 0.0], 1.0); 5];
    profile.push(([1.0, 1.0], 1.0));
    let s = turned(&profile, 10.0, [10.0, 20.0, 30.0]);
    assert_eq!(collapse(&s, [0.0, 0.3], 0.0).unwrap(), Collapse::V);
    assert_eq!(collapse(&s, [0.5, 0.3], 0.0).unwrap(), Collapse::Neither);

    let o = [0.0, 0.0, 0.0];
    let generic = |i: usize, j: usize| [i as f64, j as f64, (i * j) as f64 + 1.0];
    let rows: Vec<Vec<([f64; 3], f64)>> = (0..3)
        .map(|i| {
            (0..6)
                .map(|j| {
                    let p = match (i, j) {
                        (0, 5) => [0.0, 0.0, 1.0],
                        (0, _) | (1, 0) => o,
                        _ => generic(i, j),
                    };
                    (p, 1.0)
                })
                .collect()
        })
        .collect();
    let s = spline([2, 5], &rows);
    let p = s.partials(0.0, 0.0).unwrap();
    assert!(
        (1..=4).all(|k| p[0][k].len() == 0.0),
        "the row's derivatives to the fourth vanish"
    );
    assert_ne!(
        s.derivs(0.0, 1.0).unwrap().0,
        s.derivs(0.0, 0.0).unwrap().0,
        "the row moves"
    );
    for v in [0.0, 0.5] {
        assert_eq!(
            collapse(&s, [0.0, v], 0.0).unwrap(),
            Collapse::Neither,
            "{v}"
        );
    }
}

/// An isolated singular point (P00 = P10 = P01: `S_u = S_v = 0` at the corner, `S_uu`, `S_vv` not)
/// is no pole; a doubly degenerate corner (rows u = 0 and v = 0 both at one point) is collapsed
/// along both directions there and along each side's own direction elsewhere on it.
#[test]
fn a_singular_point_is_no_pole_and_a_doubly_degenerate_corner_is_both() {
    let o = [5.0, 6.0, 7.0];
    let generic = |i: usize, j: usize| [5.0 + i as f64, 6.0 + j as f64, 7.0 + (i * i + j) as f64];
    let net = |at_o: &dyn Fn(usize, usize) -> bool| -> EvalSurface {
        let rows: Vec<Vec<([f64; 3], f64)>> = (0..3)
            .map(|i| {
                (0..3)
                    .map(|j| (if at_o(i, j) { o } else { generic(i, j) }, 1.0))
                    .collect()
            })
            .collect();
        spline([2, 2], &rows)
    };
    let singular = net(&|i, j| (i, j) == (0, 0) || (i, j) == (1, 0) || (i, j) == (0, 1));
    let (_, su, sv) = singular.derivs(0.0, 0.0).unwrap();
    assert!(su.len() == 0.0 && sv.len() == 0.0);
    assert_eq!(
        collapse(&singular, [0.0, 0.0], 0.0).unwrap(),
        Collapse::Neither
    );

    let corner = net(&|i, j| i == 0 || j == 0);
    assert_eq!(collapse(&corner, [0.0, 0.0], 0.0).unwrap(), Collapse::Both);
    assert_eq!(collapse(&corner, [0.0, 0.5], 0.0).unwrap(), Collapse::V);
    assert_eq!(collapse(&corner, [0.5, 0.0], 0.0).unwrap(), Collapse::U);
    assert_eq!(
        collapse(&corner, [0.5, 0.5], 0.0).unwrap(),
        Collapse::Neither
    );
}

/// A sphere as a revolution (two rational quarter arcs from pole to pole), off the origin and at
/// scales 1e-3 and 1e3: collapsed along the angle at both poles, nowhere else.
#[test]
fn a_sphere_of_revolution_is_collapsed_at_its_poles_only() {
    for (r, o) in [
        (2.0, [0.0, 0.0, 0.0]),
        (2e-3, [1.0, -1.0, 2.0]),
        (2e3, [3e3, 0.0, -4e3]),
    ] {
        let pts = [[0.0, -r], [r, -r], [r, 0.0], [r, r], [0.0, r]];
        let ws = [1.0, MID, 1.0, MID, 1.0];
        let ctrl = pts
            .iter()
            .zip(ws)
            .map(|([x, z], w)| [w * (o[0] + x), w * o[1], w * (o[2] + z), w])
            .collect();
        let s = EvalSurface::Revolved {
            profile: EvalCurve::new(2, vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0], ctrl),
            origin: Vec3::new(o[0], o[1], o[2]),
            axis: Vec3::Z,
            transform: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
            ],
            swap: false,
        };
        for a in [0.0, 1.0, 4.0] {
            assert_eq!(collapse(&s, [0.0, a], 0.0).unwrap(), Collapse::V, "{r}");
            assert_eq!(collapse(&s, [1.0, a], 0.0).unwrap(), Collapse::V, "{r}");
            assert_eq!(
                collapse(&s, [0.25, a], 0.0).unwrap(),
                Collapse::Neither,
                "{r}"
            );
            assert_eq!(
                collapse(&s, [0.25, a], 1e-3 * r).unwrap(),
                Collapse::Neither,
                "{r}"
            );
        }
    }
}

/// What cannot be decided is refused, never guessed: parameters or a reach that are not finite,
/// a net with a point that is not finite, a parameter off the spline, a direction past v.
#[test]
fn what_cannot_be_decided_is_refused() {
    let s = cone(10.0, [1.0, 2.0, 3.0]);
    assert!(collapse(&s, [f64::NAN, 0.3], 0.0).is_err());
    assert!(collapse(&s, [0.0, f64::INFINITY], 0.0).is_err());
    assert!(collapse(&s, [0.0, 0.3], f64::NAN).is_err());
    assert!(collapse(&s, [0.0, 0.3], -1.0).is_err());
    assert!(collapse(&s, [3.0, 0.3], 0.0).is_err(), "off the spline");
    assert!(degenerate(&s, [0.0, 0.3], 2).is_err());
    let mut rows = vec![vec![([0.0, 0.0, 0.0], 1.0); 3]; 3];
    rows[2][2] = ([f64::NAN, 0.0, 0.0], 1.0);
    let bad = spline([2, 2], &rows);
    assert!(collapse(&bad, [0.0, 0.0], 0.0).is_err());
}
