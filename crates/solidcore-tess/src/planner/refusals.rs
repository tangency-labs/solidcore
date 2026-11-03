//! Faces on a surface the kernel's evaluator refuses on one circle: the planner leaves such a
//! face out whole and names the point, wherever it meets it (a node, a middle the sag check
//! needs, a node of a seam cut, or a point a size rule probes between the nodes it places).
//! Nothing is substituted for the point, and no point is skipped.

use solidcore::eval::{EvalCurve, EvalSurface};
use solidcore_convert as kernel;
use solidcore_record::Vec3;

use crate::vec::VecExt as _;

use super::size::Rules;
use super::{NoValue, plan_or_fill, seam, uv};

/// An exact cylinder of radius 3 about z, its profile `u` in [0, 1] running from z = 0 to
/// z = 40 (`z = 40 u`), `v` the angle: a rational cubic whose weight has a double root at
/// `hole` (0.5 or 0.25), where both its numerator and weight vanish. The evaluator has no point
/// on the circle u = `hole` (0 / 0) and the cylinder's own everywhere else. The weights are
/// exact in binary, so the root is exactly zero there.
pub(crate) fn cylinder_without(hole: f64) -> EvalSurface {
    // Homogeneous control points (w x, w y, w z, w): the weight is (1 - u / hole)² scaled, and
    // w z is 40 u times it, in the cubic Bernstein basis.
    let ctrl = if hole == 0.5 {
        // Weight 3 (1 - 2u)²: Bernstein (3, -1, -1, 3); w z = 40 u · weight: (0, 40, -80, 120).
        vec![
            [9.0, 0.0, 0.0, 3.0],
            [-3.0, 0.0, 40.0, -1.0],
            [-3.0, 0.0, -80.0, -1.0],
            [9.0, 0.0, 120.0, 3.0],
        ]
    } else {
        assert_eq!(hole, 0.25);
        // Weight 3 (1 - 4u)²: Bernstein (3, -5, 3, 27); w z = 40 u · weight: (0, 40, -240, 1080).
        vec![
            [9.0, 0.0, 0.0, 3.0],
            [-15.0, 0.0, 40.0, -5.0],
            [9.0, 0.0, -240.0, 3.0],
            [81.0, 0.0, 1080.0, 27.0],
        ]
    };
    let s = EvalSurface::Revolved {
        profile: EvalCurve::new(3, vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0], ctrl),
        origin: solidcore::Vec3::ZERO,
        axis: solidcore::Vec3::Z,
        transform: [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
        ],
        swap: false,
    };
    assert!(
        s.derivs(hole, 0.3).is_none(),
        "no point on the circle u = {hole}"
    );
    for u in [0.0, 0.125, 0.375, 0.5, 0.625, 0.75, 1.0]
        .into_iter()
        .filter(|u| *u != hole)
    {
        let p = point(&s, [u, 0.3]);
        let want = Vec3::new(3.0 * 0.3f64.cos(), 3.0 * 0.3f64.sin(), 40.0 * u);
        assert!(p.dist(want) < 1e-9, "the cylinder at u = {u}: {p:?}");
    }
    s
}

fn point(s: &EvalSurface, q: [f64; 2]) -> Vec3 {
    kernel::record_vec3(s.derivs(q[0], q[1]).unwrap().0)
}

/// A quarter turn of that cylinder, u in [0, 1] and v in [0, pi/2], as a loop counter-clockwise
/// in (u, v): its straight sides (u = 0 to 1 at v = 0 and v = pi/2) one segment each, its arcs
/// (u = 0 and u = 1) `n` segments.
fn quarter(s: &EvalSurface, n: u32) -> seam::Loop {
    let quarter = std::f64::consts::FRAC_PI_2;
    let at = |k: u32| quarter * f64::from(k) / f64::from(n);
    let mut uv: Vec<[f64; 2]> = vec![[0.0, 0.0]];
    uv.extend((0..=n).map(|k| [1.0, at(k)]));
    // Down the arc at u = 0, from v = pi/2 to just above the start.
    uv.extend((1..=n).rev().map(|k| [0.0, at(k)]));
    let pts = uv.iter().map(|q| point(s, *q)).collect();
    seam::Loop { uv, pts }
}

/// The face on [`cylinder_without`]`(0.5)`: whatever layout meshes it puts nodes or middles on
/// the circle u = 0.5. The face is refused with the point, not meshed around it, and not filled
/// leniently: no layout meshes a surface that is not there.
#[test]
fn a_face_whose_surface_refuses_an_inner_point_is_refused_whole() {
    let s = cylinder_without(0.5);
    let loops = [quarter(&s, 8)];
    let rules = Rules::new(0.01, 40.0, None);
    let r = plan_or_fill(
        0,
        "revolution",
        &s,
        &loops,
        &rules,
        (1.0, [40.0, 3.0]),
        0.0,
        None,
    );
    match r {
        Err(e) => assert!(
            e.to_string()
                .contains("the surface does not evaluate at [0.5"),
            "{e}"
        ),
        Ok((fp, stats)) => {
            panic!(
                "meshed with {} triangles, lenient {:?}",
                fp.triangles.len(),
                stats.lenient
            )
        }
    }
}

/// On [`cylinder_without`]`(0.25)` the straight direction needs no line inside [0, 1] (the
/// profile is straight), so no node lies on u = 0.25, and the sag check looks at middles
/// (u = 0.5) and centroids (u = 1/3, 2/3) only. The size rules' probes between the lines (the
/// quarter points of an interval, and the iso-lines across it) do reach u = 0.25: there the
/// surface refuses, and the division names the point instead of taking the interval as fine.
#[test]
fn the_division_names_a_point_between_its_lines_the_surface_does_not_evaluate_at() {
    let s = cylinder_without(0.25);
    let rules = Rules::new(0.01, 40.0, None);
    let quarter = (0.0, std::f64::consts::FRAC_PI_2);
    assert!(matches!(
        uv::divide(&s, 0, (0.0, 1.0), quarter, &rules, 50.0, 1000),
        Err(NoValue([u, _])) if u == 0.25
    ));
    assert!(matches!(
        uv::divide(&s, 1, quarter, (0.0, 1.0), &rules, 50.0, 1000),
        Err(NoValue([u, _])) if u == 0.25
    ));
    assert_eq!(
        uv::fine_between(&s, &rules, 50.0, [0.0, 0.1], [0.5, 0.1]),
        Err(NoValue([0.25, 0.1]))
    );
    // Off the circle the rules hold: the profile is straight, and 40 mm is within 50.
    assert_eq!(
        uv::fine_between(&s, &rules, 50.0, [0.5, 0.1], [1.0, 0.1]),
        Ok(true)
    );
}

/// The face on [`cylinder_without`]`(0.25)`, its arcs fine enough (16 segments on r = 3) that no
/// side across the face sags past the tolerance: its nodes stay off u = 0.25 and so do its
/// triangles' middles (see above). Taking a point the surface does not evaluate at as fine, the
/// planner meshed it as a strip of 32 triangles straight over the circle the surface does not
/// have. The size rules probe u = 0.25, and the face is refused there.
#[test]
fn a_face_whose_surface_refuses_a_point_between_its_nodes_is_refused_whole() {
    let s = cylinder_without(0.25);
    let loops = [quarter(&s, 16)];
    let rules = Rules::new(0.01, 40.0, None);
    let r = plan_or_fill(
        0,
        "revolution",
        &s,
        &loops,
        &rules,
        (50.0, [40.0, 3.0]),
        0.0,
        None,
    );
    match r {
        Err(e) => assert!(
            e.to_string()
                .contains("the surface does not evaluate at [0.25"),
            "{e}"
        ),
        Ok((fp, stats)) => {
            panic!(
                "meshed with {} triangles, lenient {:?}",
                fp.triangles.len(),
                stats.lenient
            )
        }
    }
}

/// The Delaunay sag check evaluates every triangle's corners and middle: a surface that refuses
/// a middle stops it with that point, instead of taking the triangle as on the surface.
#[test]
fn the_delaunay_sag_check_names_a_middle_the_surface_does_not_evaluate_at() {
    use super::cdt::{Fill, Sag, triangulate};
    let b = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    // A plane with no point inside the open square.
    let eval = |u: f64, v: f64| {
        let inside = u > 0.0 && u < 1.0 && v > 0.0 && v < 1.0;
        if inside {
            Err(NoValue([u, v]))
        } else {
            Ok((Vec3::new(u, v, 0.0), Some(Vec3::Z)))
        }
    };
    let sag = Sag {
        eval: &eval,
        tol: 1e-3,
        h: 10.0,
        turn: None,
        scale: [1.0, 1.0],
    };
    let fill = Fill {
        seeds: &[],
        refine: Some((sag, 100)),
        shape: false,
        lenient: false,
    };
    match triangulate(&b, std::slice::from_ref(&(0..4)), fill) {
        Err(NoValue([u, v])) => assert!(u > 0.0 && u < 1.0 && v > 0.0 && v < 1.0, "{u} {v}"),
        other => panic!("{other:?}"),
    }
}
