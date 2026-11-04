//! #2652 (solidcraft#11848; the generic blend, hub #2243): a round whose ball is as wide as
//! a face it rolls on, where that face is a B-spline. The ball's offset of that face
//! collapses onto its focal curve, which the other face's offset crosses at one point: the
//! ball's centre does not move along the edge and the round is a piece of that one ball
//! (#2403's ball triangle). #2403 read that event only off a cylinder or torus; it is now
//! read from the geometry where the round is built (`fillet_exact`'s
//! `ball_on_face_along_arc`): the centre is the same at every point of the round and the
//! ball lies on the face along its contact arc, within the solver's tolerance.
//!
//! The part: a 10 x 10 x 6 block (z from -6 to 0) whose upright edge at (10, 10) is a
//! quarter round of radius 2 about (8, 8), built as a plain rational B-spline wall (cubic in
//! z over three spans). Its bottom row's middle point is pushed out 0.5 along the diagonal,
//! so the wall is no cylinder (no natural quadric is read off it), while its top two spans
//! (z from -4 to 0) are exactly the quarter cylinder. Its top outline rounded R2 runs over
//! the wall's top arc with the ball's centre fixed at (8, 8, -2). Before the fix that round
//! was traced through the fixed centre and failed (1001); now it is the ball's piece and the
//! part is the same as the box whose upright edge was rounded R2 (an exact cylinder) and
//! then its top outline, to 1e-12 of the volume removed.
//!
//! The twisted disc of solidcraft#11848 meets this event over its side rounds' top edges:
//! see `plan/findings/general--twisted-disc-side-round-fit-and-vertex-gaps.md` for what remains
//! there.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers

use solidcore::blend::fillet_executive::api::{create_face_from_surface, fillet_edges};
use solidcore::curve::curve::Curve;
use solidcore::error::EntityRef;
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::geom::tol;
use solidcore::kernel::UpstreamFixes;
use solidcore::mass::integrator::properties::{brep_volume, PropertyAttributes};
use solidcore::math::{Point2, Point3};
use solidcore::nurbs::b_spline_curve::{BSplineCurve, BSplineCurveForm, CanonicalData, KnotType};
use solidcore::nurbs::b_spline_surface::{
    BSplineSurface, BSplineSurfaceForm, CanonicalSurfaceData,
};
use solidcore::primitive::api_primitives::create_box;
use solidcore::surface::surface::Surface;
use solidcore::topology::brep::p2::is_manifold_solid;
use solidcore::topology::brep::p5::{stitch_and_orient, StitchAndOrientOptions};
use solidcore::topology::brep::p8::{create_planar_faces_with_3d_curves, FaceBuildEnv};
use solidcore::topology::brep::validation::{check_brep, CheckLevel};
use solidcore::topology::{edge, vertex, EdgeId, FaceId};
use solidcore::{Brep, Kernel, SolidError};

/// SolidCraft's fix list (`tests/kernel_cases/fixes.json`), `fillet_blend_corner` as given.
fn kernel(blend_corner: bool) -> Kernel {
    let names: Vec<String> = serde_json::from_str(include_str!("kernel_cases/fixes.json")).unwrap();
    let mut f = UpstreamFixes::default();
    for name in &names {
        if let Some(flag) = f.flag_mut(name) {
            *flag = true;
        }
    }
    assert!(f.fillet_blend_corner);
    f.fillet_blend_corner = blend_corner;
    Kernel {
        fixes: f,
        ..Kernel::default()
    }
}

/// The block's rows in z: cubic over three spans, at their Greville abscissae, z = -6 + 6 g.
const ROWS: [f64; 6] = [0.0, 1.0 / 9.0, 1.0 / 3.0, 2.0 / 3.0, 8.0 / 9.0, 1.0];

/// How far the bottom row's middle point of the round wall is pushed out.
const BULGE: f64 = 0.5;

/// The control points of `p` on row `j` (at `z`): the round's middle point pushed out on the
/// bottom row.
fn row(p: &Piece, j: usize, z: f64) -> Vec<Point3> {
    p.points
        .iter()
        .enumerate()
        .map(|(i, q)| {
            let d = if j == 0 && p.degree == 2 && i == 1 {
                BULGE
            } else {
                0.0
            };
            Point3::new(q.x + d, q.y + d, z)
        })
        .collect()
}

/// The wall over `p`: `p` along u, z along v.
fn wall(p: &Piece) -> Surface {
    let mut pts = Vec::new();
    let mut weights = Vec::new();
    let rows: Vec<Vec<Point3>> = ROWS
        .iter()
        .enumerate()
        .map(|(j, g)| row(p, j, -6.0 + 6.0 * g))
        .collect();
    for i in 0..p.points.len() {
        for r in &rows {
            pts.push(r[i]);
            weights.push(p.weights.as_ref().map_or(1.0, |w| w[i]));
        }
    }
    Surface::from(
        BSplineSurface::create_canonical(&CanonicalSurfaceData {
            u_degree: p.degree,
            v_degree: 3,
            control_points: &pts,
            form: BSplineSurfaceForm::Unspecified,
            u_multiplicities: &p.mults,
            v_multiplicities: &[4, 1, 1, 4],
            u_knots: &p.knots,
            v_knots: &[0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0],
            knot_type: KnotType::Unspecified,
            weights: p.weights.as_ref().map(|_| weights.as_slice()),
            uv_domain: None,
        })
        .unwrap(),
    )
}

/// The block's outline, counter-clockwise: the quarter round about (8, 8) at its corner.
fn outline() -> Vec<Piece> {
    let p = |x: f64, y: f64| Point3::new(x, y, 0.0);
    vec![
        line(p(0.0, 0.0), p(10.0, 0.0)),
        line(p(10.0, 0.0), p(10.0, 8.0)),
        arc(
            p(8.0, 8.0),
            2.0,
            0.0,
            std::f64::consts::FRAC_PI_2,
            p(10.0, 8.0),
            p(8.0, 10.0),
        ),
        line(p(8.0, 10.0), p(0.0, 10.0)),
        line(p(0.0, 10.0), p(0.0, 0.0)),
    ]
}

/// The block with its B-spline round wall.
fn block(k: &mut Kernel) -> Brep {
    let pieces = outline();
    let mut b =
        solidcore::topology::brep::p3::new_brep(k, tol::USE_DEFAULT, tol::USE_DEFAULT).unwrap();
    for p in &pieces {
        create_face_from_surface(k, &mut b, &wall(p)).unwrap();
    }
    for (j, z) in [(0, -6.0), (ROWS.len() - 1, 0.0)] {
        let cap: Vec<Curve> = pieces
            .iter()
            .map(|p| {
                let pts = row(p, j, z);
                curve_of(
                    &Piece {
                        degree: p.degree,
                        points: pts,
                        weights: p.weights.clone(),
                        knots: p.knots.clone(),
                        mults: p.mults.clone(),
                    },
                    |q| q,
                )
            })
            .collect();
        with_face_env(k, &mut b, |k, b, fenv| {
            let tol = b.tolerance();
            let region = b.infinite_region();
            create_planar_faces_with_3d_curves(k, b, region, &cap, tol, fenv).unwrap();
        });
    }
    stitch_and_orient(k, &mut b, &StitchAndOrientOptions::default()).unwrap();
    b.refresh_cached_topology();
    assert!(is_manifold_solid(&b).unwrap());
    b
}

/// The box with the same outline: its upright edge at (10, 10) rounded R2, an exact
/// cylinder.
fn reference(k: &mut Kernel) -> Brep {
    let mut b = create_box(k, Point3::new(0.0, 0.0, -6.0), 10.0, 10.0, 6.0).unwrap();
    let upright: Vec<EdgeId> = b
        .edges()
        .ids()
        .filter(|&e| {
            ends(&b, e)
                .iter()
                .all(|p| (p.x - 10.0).abs() < 1e-9 && (p.y - 10.0).abs() < 1e-9)
        })
        .collect();
    assert_eq!(upright.len(), 1);
    fillet_edges(k, &mut b, &upright, 2.0, 1, 1, 1.0).unwrap();
    b.refresh_cached_topology();
    b
}

fn top_outline(b: &Brep) -> Vec<EdgeId> {
    b.edges()
        .ids()
        .filter(|&e| ends(b, e).iter().all(|p| p.z.abs() < 1e-9))
        .collect()
}

fn volume(k: &Kernel, b: &mut Brep) -> f64 {
    let store = AttributeStore::default();
    let attrs = PropertyAttributes {
        store: &store,
        key_of: &|_| ObjectKey(0),
    };
    brep_volume(k, b, 1e-11, &attrs).unwrap()
}

/// The volume the top outline's round R2 removes, and the part after it.
fn round_top(k: &mut Kernel, mut b: Brep) -> solidcore::Result<(f64, Brep)> {
    let before = volume(k, &mut b);
    let top = top_outline(&b);
    assert_eq!(top.len(), 5);
    fillet_edges(k, &mut b, &top, 2.0, 1, 1, 1.0)?;
    let after = volume(k, &mut b);
    Ok((before - after, b))
}

/// The faces every point of whose (natural) surface is on the ball about `c` of radius 2.
fn faces_on_the_ball(k: &Kernel, b: &Brep, c: Point3) -> Vec<FaceId> {
    b.faces()
        .ids()
        .filter(|&f| {
            let s = b.face(f).unwrap().surface.as_ref().unwrap();
            let d = s.natural_uv_domain_with_nurb(k).unwrap();
            (0..=4).all(|i| {
                (0..=4).all(|j| {
                    let uv = Point2::new(
                        d.min().x + (d.max().x - d.min().x) * f64::from(i) / 4.0,
                        d.min().y + (d.max().y - d.min().y) * f64::from(j) / 4.0,
                    );
                    (s.evaluate_point(k, uv).unwrap().distance(c) - 2.0).abs() <= 1e-9
                })
            })
        })
        .collect()
}

#[test]
fn the_round_over_a_bspline_round_of_its_radius_is_the_ball_s_piece() {
    let mut k = kernel(true);
    let b = block(&mut k);
    // the round wall is a B-spline, the ball's radius across its top arc
    let (removed, mut b) = round_top(&mut k, b).unwrap();
    let lines = check_brep(&mut k, &mut b, CheckLevel::Level2)
        .unwrap()
        .lines();
    assert!(lines.is_empty(), "{lines:?}");
    assert!(is_manifold_solid(&b).unwrap());
    // the round over the arc: one face on the ball about the fixed centre
    let ball = faces_on_the_ball(&k, &b, Point3::new(8.0, 8.0, -2.0));
    assert_eq!(ball.len(), 1, "{ball:?}");
    // the same part as the box whose upright round is an exact cylinder
    let mut kr = kernel(true);
    let r = reference(&mut kr);
    let (want, _) = round_top(&mut kr, r).unwrap();
    assert!(
        (removed - want).abs() <= 1e-12 * want,
        "removed {removed}, the cylinder's {want}"
    );
}

/// Without `fillet_blend_corner` (#2403's flag, which gates the event) the round is traced
/// through the fixed centre and fails, as it did before the fix.
#[test]
fn without_the_blend_corner_the_round_over_the_bspline_fails() {
    let mut k = kernel(false);
    let b = block(&mut k);
    match round_top(&mut k, b) {
        Err(SolidError::Upstream { code: 1001, .. }) => {}
        other => panic!("expected 1001, got {:?}", other.map(|(v, _)| v)),
    }
}

// ---- helpers (as `issue_2648_twisted_disc_top_round.rs`'s) -----------------------------

/// A profile curve in the XY plane as a (rational) B-spline: its control points, weights
/// (`None`: 1), knots with multiplicities, and degree.
struct Piece {
    degree: usize,
    points: Vec<Point3>,
    weights: Option<Vec<f64>>,
    knots: Vec<f64>,
    mults: Vec<usize>,
}

fn line(a: Point3, b: Point3) -> Piece {
    Piece {
        degree: 1,
        points: vec![a, b],
        weights: None,
        knots: vec![0.0, 1.0],
        mults: vec![2, 2],
    }
}

/// The circular arc about `c` of radius `r` from angle `a0` to `a1` (radians, either
/// way), in segments of at most 90 degrees, from `start` to `end` (its exact ends).
fn arc(c: Point3, r: f64, a0: f64, a1: f64, start: Point3, end: Point3) -> Piece {
    let n = ((a1 - a0).abs() / std::f64::consts::FRAC_PI_2)
        .ceil()
        .max(1.0) as usize;
    let d = (a1 - a0) / n as f64;
    let w = (0.5 * d).cos();
    let at = |a: f64, rr: f64| Point3::new(c.x + rr * a.cos(), c.y + rr * a.sin(), 0.0);
    let mut points = vec![start];
    let mut weights = vec![1.0];
    let mut knots = vec![0.0];
    let mut mults = vec![3];
    for i in 0..n {
        let mid = a0 + d * (i as f64 + 0.5);
        points.push(at(mid, r / w));
        weights.push(w);
        if i + 1 == n {
            points.push(end);
        } else {
            points.push(at(a0 + d * (i as f64 + 1.0), r));
        }
        weights.push(1.0);
        knots.push((i + 1) as f64);
        mults.push(if i + 1 == n { 3 } else { 2 });
    }
    Piece {
        degree: 2,
        points,
        weights: Some(weights),
        knots,
        mults,
    }
}

/// The curve of `p`'s control points under `map`.
fn curve_of(p: &Piece, map: impl Fn(Point3) -> Point3) -> Curve {
    let pts: Vec<Point3> = p.points.iter().map(|&q| map(q)).collect();
    Curve::from(
        BSplineCurve::create_canonical(&CanonicalData {
            dim: 3,
            degree: p.degree,
            control_points: &pts,
            form: BSplineCurveForm::Unspecified,
            knot_multiplicities: &p.mults,
            knots: &p.knots,
            knot_type: KnotType::Unspecified,
            weights: p.weights.as_deref(),
            trim: None,
        })
        .unwrap(),
    )
}

fn with_face_env<T>(
    k: &mut Kernel,
    b: &mut Brep,
    body: impl FnOnce(&mut Kernel, &mut Brep, &mut FaceBuildEnv<'_, '_>) -> T,
) -> T {
    let mut state = solidcore::topology::brep::p2::ArenaNotifyState::new(b);
    let mut env = state.env(None);
    let mut marks = solidcore::topology::topology_traverser::Marks::default();
    let this_face = |f: FaceId| solidcore::topology::loopuse::notify_ref(EntityRef::Face(f));
    let mut fenv = FaceBuildEnv {
        env: &mut env,
        marks: &mut marks,
        this_face: &this_face,
    };
    body(k, b, &mut fenv)
}

fn ends(b: &Brep, e: EdgeId) -> [Point3; 2] {
    let v0 = edge::vertex(b, e).unwrap().unwrap();
    let v1 = edge::other_vertex(b, e, Some(v0)).unwrap().unwrap();
    [vertex::point(b, v0).unwrap(), vertex::point(b, v1).unwrap()]
}
