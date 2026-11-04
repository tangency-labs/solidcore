use super::*;
use crate::eval::spline::EvalSurface;

const L: f64 = 10.0;
const MID: f64 = std::f64::consts::FRAC_1_SQRT_2;

/// The quarter turn about z of a quadratic or cubic profile `(radius, height)` from the apex,
/// moved by `offset`, as the record holds a revolution: a rational spline, u the profile, v the
/// quarter arc's rational parameter.
fn turned(profile: &[[f64; 2]], offset: [f64; 3]) -> EvalSurface {
    let arc = [([1.0, 0.0], 1.0), ([1.0, 1.0], MID), ([0.0, 1.0], 1.0)];
    let mut ctrl = Vec::new();
    for ([ax, ay], w) in arc {
        for [r, z] in profile {
            let p = [offset[0] + r * ax, offset[1] + r * ay, offset[2] + z];
            ctrl.push([w * p[0], w * p[1], w * p[2], w]);
        }
    }
    let p = profile.len() - 1;
    let mut ku = vec![0.0; p + 1];
    ku.extend(vec![1.0; p + 1]);
    EvalSurface::nurbs(
        [p, 2],
        [ku, vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]],
        [p + 1, 3],
        ctrl,
    )
}

/// θ of the rational quarter arc at v.
fn angle(v: f64) -> f64 {
    let (b0, b1, b2) = ((1.0 - v) * (1.0 - v), 2.0 * v * (1.0 - v) * MID, v * v);
    (b1 + b2).atan2(b0 + b1)
}

fn close(a: Vec3, b: Vec3) -> bool {
    a.dist(b) <= 1e-13
}

/// The comparison's three faces, each at a higher-order pole where the evaluator has no normal,
/// and their limits in closed form along every direction into the face:
/// face 1, the profile r = z²/L tangent to the axis (order 2): `−e_r(θ)`;
/// face 2, the patch L (u³, u³ v, u) with its row u = 0 collapsed (order 3): (−1, 0, 0);
/// face 3, the 45° cone of profile (L s², L s²) (order 3): (−cos θ, −sin θ, 1) / √2,
/// face 3 also with its apex at (1000, −500, 200).
#[test]
fn the_three_faces_have_their_closed_form_limits() {
    let face1 = turned(&[[0.0, 0.0], [0.0, 0.5 * L], [L, L]], [0.0; 3]);
    let cone = |o| turned(&[[0.0, 0.0], [0.0, 0.0], [L, L]], o);
    let cube = [0.0, 0.0, 0.0, 1.0];
    let line = [0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0];
    let mut ctrl = Vec::new();
    for j in 0..2 {
        for i in 0..4 {
            ctrl.push([L * cube[i], L * cube[i] * f64::from(j), L * line[i], 1.0]);
        }
    }
    let face2 = EvalSurface::nurbs(
        [3, 1],
        [
            vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
        ],
        [4, 2],
        ctrl,
    );
    for v in [0.0, 0.3, 0.77, 1.0] {
        let th = angle(v);
        let cases = [
            (&face1, Vec3::new(-th.cos(), -th.sin(), 0.0)),
            (&face2, Vec3::new(-1.0, 0.0, 0.0)),
            (&cone([0.0; 3]), Vec3::new(-th.cos(), -th.sin(), 1.0) * MID),
            (
                &cone([1000.0, -500.0, 200.0]),
                Vec3::new(-th.cos(), -th.sin(), 1.0) * MID,
            ),
        ];
        for (k, (s, want)) in cases.iter().enumerate() {
            assert!(
                matches!(s.eval(0.0, v), Err(EvalError::Degenerate(_))),
                "case {k} at v = {v}: no normal of its own"
            );
            for dir in [[1.0, 0.0], [1.0, 0.3], [0.5, -0.2]] {
                let n = s.normal_along(0.0, v, dir).unwrap();
                assert!(close(n, *want), "case {k} v {v} {dir:?}: {n:?} vs {want:?}");
            }
        }
    }
}

/// Off a pole the normal along a direction is the normal itself.
#[test]
fn off_a_pole_it_is_the_normal() {
    let s = turned(&[[0.0, 0.0], [0.0, 0.0], [L, L]], [3.0, 2.0, 1.0]);
    let (_, n) = s.eval(0.5, 0.3).unwrap();
    assert!(close(s.normal_along(0.5, 0.3, [0.2, 1.0]).unwrap(), n));
}

/// A cone whose profile (L s³, L s³) starts with its first three control points at the apex:
/// `S_u × S_v ~ s⁵`, so `N_0` to `N_3` vanish and the pole is of order five. Refused with
/// `Degenerate`, which is what a mesher's fallback keys on; so is a direction that is none.
#[test]
fn a_pole_of_order_above_three_is_refused() {
    let s = turned(
        &[[0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [L, L]],
        [7.0, 0.0, -3.0],
    );
    for v in [0.0, 0.4] {
        assert!(matches!(
            s.normal_along(0.0, v, [1.0, 0.0]),
            Err(EvalError::Degenerate(_))
        ));
    }
    let face = turned(&[[0.0, 0.0], [0.0, 0.0], [L, L]], [0.0; 3]);
    assert!(matches!(
        face.normal_along(0.0, 0.4, [0.0, 0.0]),
        Err(EvalError::Degenerate(_))
    ));
    assert!(matches!(
        face.normal_along(0.0, 0.4, [f64::NAN, 1.0]),
        Err(EvalError::Degenerate(_))
    ));
}
