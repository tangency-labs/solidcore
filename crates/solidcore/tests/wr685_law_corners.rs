//! WR-685 law-corner addendum (#2273; the fork's tangency-labs/smlib 7babe17, 58d03ef and
//! 62fcc82 in 5e8c9c1): a variable round whose radius law has a corner is the union of
//! its balls there. Each side is its own exact piece, ending on the corner ball's touch
//! circle of its own side through the same two rail points. Where the law dips, the
//! circle after lies ahead of the one before and the ball's sphere between them is a face
//! of its own (`sphere_zone`); where it peaks, the sides cross.
//!
//! The corner's geometry is ported (`solidcore::blend::law_corners`) and checked here
//! against closed forms and the fork's build. The round itself is not built yet: the
//! fork keeps the pieces from its piecewise skin between the law's corners and splits the
//! faces there (tangency-labs/smlib@6eecf04 and 5fbffc9), which solidcore does not have
//! (solidcore#1403, #1828). The fork's capi cases are kept below, ignored until then.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use solidcore::blend::fillet_executive::session::FilletSession;
use solidcore::blend::fillet_solver::p1::FilletSolverType;
use solidcore::blend::fillet_standard_solver::p1::{CircularCrossSectionFsg, FilletLaw};
use solidcore::blend::law_corners::{
    curve_between, find_law_corner_pieces, has_built_law_corners, sphere_zone, ArcApproximation,
    LawCornerKind, LawCornerPieces,
};
use solidcore::curve::curve::Curve;
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::geom::extent1d::Extent1d;
use solidcore::geom::vector2d::Point2;
use solidcore::kernel::UpstreamFixes;
use solidcore::mass::integrator::properties::{brep_volume, PropertyAttributes};
use solidcore::math::{Point3, Vec3};
use solidcore::nurbs::b_spline_curve::p1::{
    BSplineCurve, BSplineCurveForm, CanonicalData, KnotType,
};
use solidcore::nurbs::b_spline_surface::{BSplineSurface, SurfParamType};
use solidcore::primitive::api_primitives::create_box;
use solidcore::surface::surface::Surface;
use solidcore::topology::{edge, EdgeId};
use solidcore::{Brep, Kernel};

// ---- the corner ball's zone at a dip ------------------------------------------------------

/// The dip corner of the law 4 -> 2 -> 4 on the 20 cube's edge y = z = 20, at x = 10: the
/// corner ball's centre and radius, and the two rail points.
const CENTRE: Point3 = Point3::new(10.0, 18.0, 18.0);
const ONE: Point3 = Point3::new(10.0, 18.0, 20.0);
const TWO: Point3 = Point3::new(10.0, 20.0, 18.0);

/// The middles of the two sides' touch arcs at the corner, as the fork's scratch driver
/// computed them: the side before lies in the plane `(p - CENTRE).(1, 0.2, 0.2) = 0.4`
/// (the envelope's `<p - c, c'> = -r r'`, r falling 0.2 a unit), the side after in
/// `(p - CENTRE).(1, -0.2, -0.2) = -0.4`.
const ON_FIRST: Point3 = Point3::new(9.8362146258571865, 19.409463435357033, 19.409463435357033);
const ON_SECOND: Point3 = Point3::new(10.163785374142813, 19.409463435357033, 19.409463435357033);

/// The fork's `sm_SphereZone` for those points (5e8c9c1 built locally, the function
/// copied into a scratch driver as it is static there): the control net row by row in u,
/// Euclidean, and its weights.
const FORK_NET: [(f64, f64, f64, f64); 9] = [
    (10.0, 20.0, 18.0, 8.0),
    (9.5999999999999996, 20.0, 20.0, 1.9611613513818402),
    (10.0, 18.0, 20.0, 1.0),
    (10.0, 20.0, 18.0, 7.8446454055273609),
    (10.0, 20.0, 20.0, 2.0),
    (10.0, 18.0, 20.0, 0.98058067569092011),
    (10.0, 20.0, 18.0, 8.0),
    (10.4, 20.0, 20.0, 1.9611613513818402),
    (10.0, 18.0, 20.0, 1.0),
];

#[test]
fn a_dips_zone_is_the_fork_s_patch_bit_for_bit() {
    let zone = sphere_zone(ONE, TWO, ON_FIRST, ON_SECOND).unwrap();
    let net = zone.control_point_net().unwrap();
    assert_eq!((net.u_count, net.v_count), (3, 3));
    for (i, &(x, y, z, w)) in FORK_NET.iter().enumerate() {
        let p = net.points[i];
        assert_eq!(
            [
                p.x.to_bits(),
                p.y.to_bits(),
                p.z.to_bits(),
                net.weights[i].to_bits()
            ],
            [x.to_bits(), y.to_bits(), z.to_bits(), w.to_bits()],
            "control point {i}: {p:?} w {} against the fork's ({x}, {y}, {z}) w {w}",
            net.weights[i]
        );
    }
}

#[test]
fn a_dips_zone_lies_on_the_corner_ball_between_its_touch_circles() {
    let zone = sphere_zone(ONE, TWO, ON_FIRST, ON_SECOND).unwrap();
    let before = Vec3::new(1.0, 0.2, 0.2);
    let after = Vec3::new(1.0, -0.2, -0.2);
    for i in 0..=8 {
        for j in 0..=8 {
            let (u, v) = (f64::from(i) / 8.0, f64::from(j) / 8.0);
            let x = zone.evaluate_point(Point2::new(u, v)).unwrap();
            // on the ball
            assert!(
                ((x - CENTRE).length() - 2.0).abs() < 1e-13,
                "{u} {v}: {x:?}"
            );
            // the u ends are the two touch circles
            if i == 0 {
                assert!(((x - CENTRE).dot(before) - 0.4).abs() < 1e-13, "{v}: {x:?}");
            }
            if i == 8 {
                assert!(((x - CENTRE).dot(after) + 0.4).abs() < 1e-13, "{v}: {x:?}");
            }
            // the v ends are poles at the rail points
            if j == 0 {
                assert!(x.distance(TWO) < 1e-13, "{u}: {x:?}");
            }
            if j == 8 {
                assert!(x.distance(ONE) < 1e-13, "{u}: {x:?}");
            }
            // between the circles: ahead of the first, behind the second along the edge
            assert!(
                x.x >= ON_FIRST.x - 1e-12 && x.x <= ON_SECOND.x + 1e-12,
                "{x:?}"
            );
        }
    }
    // the ends' arcs run through the arcs' middles given
    for (u, on) in [(0.0, ON_FIRST), (1.0, ON_SECOND)] {
        let near = (0..=1000)
            .map(|j| {
                let v = f64::from(j) / 1000.0;
                zone.evaluate_point(Point2::new(u, v)).unwrap().distance(on)
            })
            .fold(f64::INFINITY, f64::min);
        assert!(near < 1e-2, "{u}: {near}");
    }
}

#[test]
fn a_zone_needs_two_rail_points_and_a_wedge_under_a_half_turn() {
    // one point twice
    assert!(sphere_zone(ONE, ONE, ON_FIRST, ON_SECOND).is_err());
    // the two arcs of one circle between the rail points: inverted, opposite rays from the
    // image of the second point, a half turn
    let n = Vec3::new(1.0, 0.2, 0.2);
    let circle_centre = CENTRE + n * (0.4 / n.dot(n));
    let other_arc = circle_centre * 2.0 - ON_FIRST;
    assert!(sphere_zone(ONE, TWO, ON_FIRST, other_arc).is_err());
    // the given pair is well inside a half turn
    assert!(sphere_zone(ONE, TWO, ON_FIRST, ON_SECOND).is_ok());
}

#[test]
fn a_crossing_curve_runs_from_the_edge_s_start_to_its_end() {
    let k = Kernel::default();
    let (a, b) = (Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 2.0, 3.0));
    let line = Curve::from(BSplineCurve::create_line_segment(3, a, b).unwrap());
    let same = curve_between(&k, &line, a, b, 1e-6)
        .unwrap()
        .expect("forward");
    let ivl = same.natural_interval();
    assert_eq!(same.evaluate_point(&k, ivl.min()).unwrap(), a);
    let back = curve_between(&k, &line, b, a, 1e-6)
        .unwrap()
        .expect("reversed");
    let ivl = back.natural_interval();
    assert!(back.evaluate_point(&k, ivl.min()).unwrap().distance(b) < 1e-15);
    assert!(back.evaluate_point(&k, ivl.max()).unwrap().distance(a) < 1e-15);
    // ends within the tolerance only
    let off = Point3::new(1.0, 2.0, 3.0 + 1e-3);
    assert!(curve_between(&k, &line, a, off, 1e-6).unwrap().is_none());
    assert!(curve_between(&k, &line, a, off, 1e-2).unwrap().is_some());
}

// ---- each corner's kind, from the round's pieces ------------------------------------------

/// The degree-1 curve through `points` at the parameters `at`.
fn polyline(points: &[Point3], at: &[f64]) -> BSplineCurve {
    let mut mults = vec![1; at.len()];
    mults[0] = 2;
    *mults.last_mut().unwrap() = 2;
    BSplineCurve::create_canonical(&CanonicalData {
        dim: 3,
        degree: 1,
        control_points: points,
        form: BSplineCurveForm::default(),
        knot_multiplicities: &mults,
        knots: at,
        knot_type: KnotType::default(),
        weights: None,
        trim: None,
    })
    .unwrap()
}

/// A round on the 20 cube's edge y = z = 20 whose law is linear between the radii `r` at
/// x = 0, 10 and 20: the centre line and the rails on y = 20 and z = 20, each parameter
/// x, whole and per piece.
fn law_curves(r: [f64; 3]) -> ([BSplineCurve; 3], Vec<[BSplineCurve; 3]>) {
    let x = [0.0, 10.0, 20.0];
    let centre = |i: usize| Point3::new(x[i], 20.0 - r[i], 20.0 - r[i]);
    let rail1 = |i: usize| Point3::new(x[i], 20.0, 20.0 - r[i]);
    let rail2 = |i: usize| Point3::new(x[i], 20.0 - r[i], 20.0);
    let whole = |f: &dyn Fn(usize) -> Point3| polyline(&[f(0), f(1), f(2)], &x);
    let piece = |f: &dyn Fn(usize) -> Point3, i: usize| polyline(&[f(i), f(i + 1)], &x[i..i + 2]);
    let pieces = (0..2)
        .map(|i| [piece(&centre, i), piece(&rail1, i), piece(&rail2, i)])
        .collect();
    ([whole(&centre), whole(&rail1), whole(&rail2)], pieces)
}

const APPROX_TOL: f64 = 1e-5;

/// The fork's section: exact arcs (its `SmArcSection` default).
const EXACT: ArcApproximation = ArcApproximation {
    approximate_arcs: false,
    tol_percent_radius: 0.05,
};

/// The round's pieces skinned on the balls' touch circles (the G^n section flag on, a
/// variable-radius solver), joined as the fork's `CircleSurface` joins them, and each
/// corner classified.
fn law_corner_pieces(k: &Kernel, r: [f64; 3]) -> (LawCornerPieces, BSplineSurface) {
    let ([center_line, rail1, _], pieces) = law_curves(r);
    let generator = CircularCrossSectionFsg {
        approximate_arcs: false,
        ..CircularCrossSectionFsg::default()
    };
    let mut surfaces = Vec::new();
    let mut joined: Option<BSplineSurface> = None;
    let mut joined_ends = Vec::new();
    for [mut c, mut r1, mut r2] in pieces {
        let mut piece = generator
            .create_surface(
                k,
                APPROX_TOL,
                FilletSolverType::VariableRadius,
                &mut c,
                &mut r1,
                &mut r2,
            )
            .unwrap();
        surfaces.push(piece.clone());
        joined = Some(match joined {
            None => piece,
            Some(mut j) => j.join_surface(&mut piece, SurfParamType::U, 0.0).unwrap(),
        });
        joined_ends.push(joined.as_ref().unwrap().natural_uv_domain().max().x);
    }
    let kept = LawCornerPieces::classify(
        &center_line,
        &rail1,
        &[0.0, 10.0, 20.0],
        surfaces,
        joined_ends,
        APPROX_TOL,
        EXACT,
    )
    .unwrap();
    (kept, joined.unwrap())
}

fn gn_kernel() -> Kernel {
    let mut k = Kernel::default();
    k.fixes.fillet_gn_section_generator = true;
    k
}

#[test]
fn a_law_that_dips_has_the_corner_ball_s_zone_between_its_sides() {
    let k = gn_kernel();
    let (kept, _) = law_corner_pieces(&k, [4.0, 2.0, 4.0]);
    assert_eq!(kept.kinds, vec![LawCornerKind::Dip]);
    assert_eq!(kept.centres, vec![CENTRE]);
    assert_eq!(kept.radii, vec![2.0]);
    assert_eq!(kept.tol_3d, APPROX_TOL);
    assert!(kept.has_dip_or_peak());
    // the sides end on the touch circles the zone runs between
    let before = &kept.pieces[0];
    let after = &kept.pieces[1];
    for (side, at_max, plane, rhs) in [
        (before, true, Vec3::new(1.0, 0.2, 0.2), 0.4),
        (after, false, Vec3::new(1.0, -0.2, -0.2), -0.4),
    ] {
        let d = side.natural_uv_domain();
        let u = if at_max { d.max().x } else { d.min().x };
        for j in 0..=8 {
            let v = d.min().y + (d.max().y - d.min().y) * f64::from(j) / 8.0;
            let x = side.evaluate_point(Point2::new(u, v)).unwrap();
            assert!(((x - CENTRE).length() - 2.0).abs() < 1e-9, "{x:?}");
            assert!(((x - CENTRE).dot(plane) - rhs).abs() < 1e-9, "{x:?}");
        }
    }
}

#[test]
fn a_law_that_peaks_has_sides_that_run_past_each_other() {
    let k = gn_kernel();
    let (kept, _) = law_corner_pieces(&k, [2.0, 4.0, 2.0]);
    assert_eq!(kept.kinds, vec![LawCornerKind::Peak]);
    assert_eq!(kept.centres, vec![Point3::new(10.0, 16.0, 16.0)]);
    assert_eq!(kept.radii, vec![4.0]);
}

#[test]
fn a_constant_radius_corner_matches() {
    let k = gn_kernel();
    let (kept, _) = law_corner_pieces(&k, [3.0, 3.0, 3.0]);
    assert_eq!(kept.kinds, vec![LawCornerKind::Matching]);
    assert!(!kept.has_dip_or_peak());
    assert!(!has_built_law_corners(std::slice::from_ref(&kept)));
}

#[test]
fn approximated_arcs_widen_the_pieces_accuracy() {
    let k = gn_kernel();
    let (kept, _) = law_corner_pieces(&k, [4.0, 2.0, 4.0]);
    let ([center_line, rail1, _], _) = law_curves([4.0, 2.0, 4.0]);
    let approx = LawCornerPieces::classify(
        &center_line,
        &rail1,
        &[0.0, 10.0, 20.0],
        kept.pieces.clone(),
        kept.joined_ends.clone(),
        APPROX_TOL,
        ArcApproximation {
            approximate_arcs: true,
            tol_percent_radius: 0.05,
        },
    )
    .unwrap();
    // 0.05 percent of the corner ball's radius 2
    assert_eq!(approx.tol_3d, 0.01 * 0.05 * 2.0);
}

#[test]
fn a_rounds_pieces_are_found_by_its_corners_once() {
    let k = gn_kernel();
    let (dip, joined) = law_corner_pieces(&k, [4.0, 2.0, 4.0]);
    let (peak, _) = law_corner_pieces(&k, [2.0, 4.0, 2.0]);
    let mut kept = vec![peak, dip];
    assert!(has_built_law_corners(&kept));
    let joined = Surface::from(joined);
    // the joined surface turns where both laws have their corner: the first not taken
    assert_eq!(
        find_law_corner_pieces(&k, &mut kept, &joined).unwrap(),
        Some(0)
    );
    assert!(kept[0].taken && !kept[1].taken);
    assert_eq!(
        find_law_corner_pieces(&k, &mut kept, &joined).unwrap(),
        Some(1)
    );
    assert_eq!(
        find_law_corner_pieces(&k, &mut kept, &joined).unwrap(),
        None
    );
    // a surface that does not turn there has none
    let mut fresh = vec![law_corner_pieces(&k, [4.0, 2.0, 4.0]).0];
    let smooth = Surface::from(fresh[0].pieces[0].clone());
    assert_eq!(
        find_law_corner_pieces(&k, &mut fresh, &smooth).unwrap(),
        None
    );
    assert!(!fresh[0].taken);
}

// ---- the fork's capi cases (not built yet) -------------------------------------------------
//
// At the time of the port (solidcore main f5a031c9c), every fix on, the three cases are
// refused at the law's corner: the trace's rails run smoothly across it (no through point
// there, 6eecf04), and the fold check (`sm_BallsFold`) reads that as balls folding "where
// its radius is 2" (the dip; 2.5 for the first peak), status 1001. With the pinned
// upstream (no fix) the dip removes 35.1156 against the union of the balls' 39.5421, the
// first peak 19.1640 against 21.7999 (seven faces each), and the dips case builds an
// invalid body.

fn volume(k: &Kernel, b: &mut Brep) -> f64 {
    let store = AttributeStore::default();
    let attrs = PropertyAttributes {
        store: &store,
        key_of: &|_| ObjectKey(0),
    };
    brep_volume(k, b, 1e-8, &attrs).unwrap()
}

/// The edge from `p` to `q`, and whether it runs from `q` to `p`.
fn edge_between(k: &Kernel, b: &Brep, p: Point3, q: Point3) -> (EdgeId, bool) {
    b.edges()
        .ids()
        .find_map(|e| {
            let ivl = edge::interval(b, e).unwrap();
            let c = edge::curve(b, e).unwrap().unwrap();
            let s = c.evaluate_point(k, ivl.min()).unwrap();
            let t = c.evaluate_point(k, ivl.max()).unwrap();
            if s.distance(p) < 1e-9 && t.distance(q) < 1e-9 {
                Some((e, false))
            } else if s.distance(q) < 1e-9 && t.distance(p) < 1e-9 {
                Some((e, true))
            } else {
                None
            }
        })
        .expect("edge")
}

/// The fork capi's `law_along`: the B-spline law of `degree` with the full `knots` and
/// control radii `radii`, the radius as x and the parameter as y (control points at the
/// Greville abscissae), reversed when the edge runs the other way.
fn spline_law(degree: usize, knots: &[f64], radii: &[f64], reverse: bool) -> FilletLaw {
    let n = knots.len();
    let kv: Vec<f64> = (0..n)
        .map(|i| {
            if reverse {
                1.0 - knots[n - 1 - i]
            } else {
                knots[i]
            }
        })
        .collect();
    let control: Vec<Vec3> = (0..radii.len())
        .map(|i| {
            let g: f64 = (1..=degree).map(|j| kv[i + j]).sum::<f64>() / degree as f64;
            let r = radii[if reverse { radii.len() - 1 - i } else { i }];
            Vec3::new(r, g, 0.0)
        })
        .collect();
    let mut unique: Vec<f64> = Vec::new();
    let mut mults: Vec<usize> = Vec::new();
    for &t in &kv {
        if unique.last() == Some(&t) {
            *mults.last_mut().unwrap() += 1;
        } else {
            unique.push(t);
            mults.push(1);
        }
    }
    let curve = BSplineCurve::create_canonical(&CanonicalData {
        dim: 3,
        degree,
        control_points: &control,
        form: BSplineCurveForm::default(),
        knot_multiplicities: &mults,
        knots: &unique,
        knot_type: KnotType::default(),
        weights: None,
        trim: None,
    })
    .unwrap();
    FilletLaw::bspline(curve, &Extent1d::from_min_max(0.0, 1.0))
}

/// The box `size` after the law rounds its edge from `from` to `to` (the law's start at
/// `from`): the volume removed and the face count, as `smlib_fillet_variable` runs it,
/// every fix on (the fork's own behaviour).
fn law_round(
    size: [f64; 3],
    from: Point3,
    to: Point3,
    degree: usize,
    knots: &[f64],
    radii: &[f64],
) -> (f64, usize) {
    let mut k = Kernel::default();
    k.fixes.clone_from(&UpstreamFixes::all());
    let mut b = create_box(&mut k, Vec3::ZERO, size[0], size[1], size[2]).unwrap();
    let (e, reversed) = edge_between(&k, &b, from, to);
    let law = spline_law(degree, knots, radii, reversed);
    let base = radii[0].max(radii[radii.len() - 1]);
    {
        let mut session = FilletSession::new(&mut k, &mut b).unwrap();
        session
            .prepare_variable_radius_edges(&mut k, &[e], base, law, 1, 1)
            .unwrap();
        session.do_fillet(&mut k).unwrap();
    }
    let faces = b.faces().ids().count();
    (size[0] * size[1] * size[2] - volume(&k, &mut b), faces)
}

/// The fork's 5e8c9c1 build, through its C API (a scratch driver running the capi cases):
/// the volume removed at `smlib_brep_volume`'s 1e-7.
const BLOCKED: &str = "the law's corners are not trace through points and the round is not \
    skinned piece by piece or split there (tangency-labs/smlib@6eecf04, 5fbffc9; solidcore#1403, #1828)";

#[test]
#[ignore = "blocked: law corners need the fork's 6eecf04 (solidcore#1403, #1828)"]
fn law_peak() {
    // 2 -> w -> 2 on the 10 x 20 x 30 box's top edge x = 10, z = 30, from y = 0: the
    // block's six faces and the two sides' rounds, meeting where they cross
    for (w, want, fork) in [
        (2.5, 21.799882, 21.799877917178),
        (3.0, 27.093683, 27.093677260694),
        (3.2, 29.389966, 29.389911194516),
        (4.0, 39.541034, 39.541083296321),
    ] {
        let (removed, faces) = law_round(
            [10.0, 20.0, 30.0],
            Point3::new(10.0, 0.0, 30.0),
            Point3::new(10.0, 20.0, 30.0),
            1,
            &[0.0, 0.0, 0.5, 1.0, 1.0],
            &[2.0, w, 2.0],
        );
        eprintln!(
            "peak {w}: removed {removed}, want {want} (fork {fork}), {faces} faces; {BLOCKED}"
        );
        assert!((removed - want).abs() <= 1e-4 * want, "{removed}");
        assert_eq!(faces, 8);
    }
}

#[test]
#[ignore = "blocked: law corners need the fork's 6eecf04 (solidcore#1403, #1828)"]
fn law_dip() {
    // 4 -> 2 -> 4 on the 20 cube's edge y = z = 20 from either end, and the same law as a
    // quadratic with a corner at 0.5: the cube's six faces, the two sides' rounds and the
    // corner ball's zone between them (the fork: 39.542055909305 linear, 39.542055908618
    // quadratic)
    let (a, b) = (Point3::new(0.0, 20.0, 20.0), Point3::new(20.0, 20.0, 20.0));
    for (from, to, degree, knots, radii) in [
        (
            a,
            b,
            1,
            &[0.0, 0.0, 0.5, 1.0, 1.0][..],
            &[4.0, 2.0, 4.0][..],
        ),
        (
            b,
            a,
            1,
            &[0.0, 0.0, 0.5, 1.0, 1.0][..],
            &[4.0, 2.0, 4.0][..],
        ),
        (
            a,
            b,
            2,
            &[0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0][..],
            &[4.0, 3.0, 2.0, 3.0, 4.0][..],
        ),
    ] {
        let (removed, faces) = law_round([20.0; 3], from, to, degree, knots, radii);
        eprintln!("dip: removed {removed}, want 39.542073, {faces} faces; {BLOCKED}");
        assert!((removed - 39.542073).abs() <= 1e-4 * 39.542073, "{removed}");
        assert_eq!(faces, 9);
    }
}

#[test]
#[ignore = "blocked: law corners need the fork's 6eecf04 (solidcore#1403, #1828)"]
fn law_dips() {
    // 4, 2, 4, 2 at 0, 0.25, 0.75, 1 on the 20 cube's edge y = z = 20: a dip, then a peak;
    // the cube's six faces, three rounds and the zone at the dip (the fork:
    // 38.818141402826)
    let (removed, faces) = law_round(
        [20.0; 3],
        Point3::new(0.0, 20.0, 20.0),
        Point3::new(20.0, 20.0, 20.0),
        1,
        &[0.0, 0.0, 0.25, 0.75, 1.0, 1.0],
        &[4.0, 2.0, 4.0, 2.0],
    );
    eprintln!("dips: removed {removed}, want 38.818277, {faces} faces; {BLOCKED}");
    assert!((removed - 38.818277).abs() <= 1e-4 * 38.818277, "{removed}");
    assert_eq!(faces, 10);
}
