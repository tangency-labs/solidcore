//! A general-path sweep is a `ProceduralSweep` definition with a certified bound (#2572,
//! slice 2; not upstream): the profile carried along the path in the frame the sweep's
//! projection normal method gives, in closed form, beside the net `sweep_scale` fits and the
//! net's certified distance from it, `|B - S| <= bound`.
//!
//! - A circle along a line, an arc, or a line and then an arc: the definition's points are
//!   the exact cylinder, torus or both, at off-net parameters, and the net lies within the
//!   bound of them.
//! - A circle along a 3D spline of straights and R2 rounds (a battery pack's balance wire,
//!   a 1.5 mm circle swept along a fitted spline): the wall reports its definition and its
//!   bound; the definition's points are the circle about the path's point at the net's
//!   parameter, and the net lies within the bound of them at parameters off its knots.
//! - The definition follows its net's edits.
//! - Every wall face is a definition, or of a kind `data/sweep_definition_ledger/` lists.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers
use solidcore::curve::curve::Curve;
use solidcore::curve::ellipse::Ellipse;
use solidcore::curve::line::Line;
use solidcore::geom::axis2_placement::Axis2Placement;
use solidcore::geom::extent1d::Extent1d;
use solidcore::math::{Point2, Point3, Vec3};
use solidcore::nurbs::b_spline_curve::BSplineCurve;
use solidcore::nurbs::b_spline_surface::SurfParamType;
use solidcore::nurbs::crv_fit::p1::ParamMethod;
use solidcore::primitive::api_curves::create_arc;
use solidcore::primitive::primitive_creation::SweepOptions;
use solidcore::surface::procedural_sweep::{ProceduralSweep, SweepRef};
use solidcore::surface::surface::{Surface, SurfaceKind};
use solidcore::sweep::api_brep::{create_curve_sweep, create_pipe_sweep};
use solidcore::topology::brep::p1;
use solidcore::topology::face;
use solidcore::{Brep, Kernel};

use std::collections::{BTreeMap, BTreeSet};

/// The balance wire's fit points: straights joined by R2 rounds, in mm.
const FIT: &[(f64, f64, f64)] = &[
    (128.855, 22.2, 63.309),
    (130.937, 22.2, 63.309),
    (133.018, 22.2, 63.309),
    (135.1, 22.2, 63.309),
    (135.711, 22.2, 63.253),
    (136.211, 22.2, 63.087),
    (136.6, 22.2, 62.809),
    (136.878, 22.2, 62.42),
    (137.044, 22.2, 61.92),
    (137.1, 22.2, 61.309),
    (137.1, 22.2, 59.1),
    (137.1, 22.144, 58.489),
    (137.1, 21.978, 57.989),
    (137.1, 21.7, 57.6),
    (137.1, 21.311, 57.322),
    (137.1, 20.811, 57.156),
    (137.1, 20.2, 57.1),
    (137.1, 17.247, 57.1),
    (137.1, 14.295, 57.1),
    (137.1, 11.342, 57.1),
    (137.1, 8.389, 57.1),
    (137.1, 5.437, 57.1),
    (137.1, 2.484, 57.1),
    (137.1, -0.468, 57.1),
    (137.1, -3.421, 57.1),
    (137.1, -6.374, 57.1),
    (137.1, -9.326, 57.1),
    (137.1, -12.279, 57.1),
    (137.1, -15.232, 57.1),
    (137.1, -18.184, 57.1),
    (137.1, -21.137, 57.1),
    (137.1, -24.089, 57.1),
    (137.1, -27.042, 57.1),
    (137.1, -29.995, 57.1),
    (137.1, -32.947, 57.1),
    (137.1, -35.9, 57.1),
    (137.1, -36.511, 57.044),
    (137.1, -37.011, 56.878),
    (137.1, -37.4, 56.6),
    (137.1, -37.678, 56.211),
    (137.1, -37.844, 55.711),
    (137.1, -37.9, 55.1),
    (137.1, -37.9, 52.55),
    (137.156, -37.9, 51.939),
    (137.322, -37.9, 51.439),
    (137.6, -37.9, 51.05),
    (137.989, -37.9, 50.772),
    (138.489, -37.9, 50.606),
    (139.1, -37.9, 50.55),
    (142.055, -37.9, 50.55),
    (145.011, -37.9, 50.55),
    (147.966, -37.9, 50.55),
    (150.922, -37.9, 50.55),
    (153.877, -37.9, 50.55),
    (156.833, -37.9, 50.55),
    (159.788, -37.9, 50.55),
    (162.744, -37.9, 50.55),
    (165.7, -37.9, 50.55),
    (168.655, -37.9, 50.55),
    (168.67, -37.9, 48.519),
    (168.685, -37.9, 46.489),
];

/// The pipe's radius: the wire's 1.5 mm circle.
const RADIUS: f64 = 0.75;

fn spline_path() -> Curve {
    let pts: Vec<Point3> = FIT.iter().map(|&(x, y, z)| Point3::new(x, y, z)).collect();
    Curve::from(
        BSplineCurve::interpolate_points(&pts, None, 3, None, None, true, ParamMethod::ChordLength)
            .unwrap(),
    )
}

fn line_path() -> Curve {
    Curve::from(
        BSplineCurve::create_line_segment(
            3,
            Point3::new(1.0, 2.0, 3.0),
            Point3::new(4.0, 6.0, 15.0),
        )
        .unwrap(),
    )
}

fn arc_path() -> Curve {
    Curve::from(create_arc(Point3::new(0.0, 0.0, 0.0), 10.0, 0.0, 90.0).unwrap())
}

/// The options `create_pipe_sweep` sweeps with: the circle centred on the path, square to it.
fn pipe_options() -> SweepOptions {
    SweepOptions {
        move_profile_center_to_path: true,
        orient_profile_perpendicular_to_path: true,
        ..SweepOptions::default()
    }
}

/// The disc square to the start of the line-then-arc path, and the path: ten along x, then
/// a quarter arc of radius 10 turning to +y.
fn line_then_arc(k: &Kernel) -> (Curve, Curve) {
    let disc = Axis2Placement::new(
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    );
    let circle =
        Curve::from(Ellipse::create_canonical_circle(k, &disc, RADIUS, None, None).unwrap());
    let line =
        Line::create_line_segment(3, Vec3::new(0.0, 0.0, 0.0), Vec3::new(10.0, 0.0, 0.0), None)
            .unwrap();
    let arc_frame = Axis2Placement::new(
        Vec3::new(10.0, 10.0, 0.0),
        Vec3::new(0.0, -1.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
    );
    let arc = Ellipse::create_canonical_circle(
        k,
        &arc_frame,
        10.0,
        Some(&Extent1d::from_min_max(0.0, 90.0)),
        None,
    )
    .unwrap();
    let path = solidcore::curve::composite_curve::CompositeCurve::build_composites_from_curves(
        k,
        vec![Curve::from(line), Curve::from(arc)],
        1e-6,
        0.0,
    )
    .unwrap()
    .composites
    .into_iter()
    .next()
    .unwrap()
    .make_composite_nurb()
    .unwrap();
    (circle, Curve::from(path))
}

fn walls(b: &Brep) -> Vec<Surface> {
    p1::faces(b, None)
        .unwrap()
        .into_iter()
        .map(|f| face::p2::surface(b, f).unwrap().unwrap().clone())
        .collect()
}

/// The farthest the net strays from the definition over a grid of the face's `(u, v)` off
/// the net's knots, the net's own parameter domain (a long axis sampled densely), and the
/// definition's.
fn deviation(s: &SweepRef<'_>, along: usize) -> f64 {
    let dom = s.net().natural_uv_domain();
    let (u, v) = (dom.u_interval(), dom.v_interval());
    let (nu, nv) = if s.definition().definition_map().2 {
        (along, 24)
    } else {
        (24, along)
    };
    let mut far = 0.0_f64;
    for i in 0..nu {
        for j in 0..nv {
            let uv = Point2::new(
                u.min() + (u.max() - u.min()) * (f64::from(i as u32) + 0.37) / f64::from(nu as u32),
                v.min() + (v.max() - v.min()) * (f64::from(j as u32) + 0.61) / f64::from(nv as u32),
            );
            let gap = (s.net().evaluate_point(uv).unwrap()
                - s.definition().evaluate_definition_point(uv).unwrap())
            .length();
            far = far.max(gap);
        }
    }
    far
}

/// Every wall of `b` is a definition of a path motion; returns the first.
fn only_wall(b: &Brep) -> Surface {
    let w = walls(b);
    assert_eq!(w.len(), 1);
    assert!(
        w[0].procedural_sweep().is_some(),
        "the wall keeps its definition"
    );
    w.into_iter().next().unwrap()
}

/// The definition's point at `uv` is `radius` from the path's point at the net's
/// parameter `v`: the circle about the path.
fn on_tube(d: &ProceduralSweep, uv: Point2) -> f64 {
    let m = d.path().unwrap();
    let wt = d.definition_parameters(uv);
    let (t, _) = m.path_parameter(wt.y);
    let on_path = m.path().evaluate_point(t).unwrap();
    (d.evaluate_definition_point(uv).unwrap() - on_path).length()
}

/// A circle along a line is the exact cylinder: the definition's points are `RADIUS` from
/// the axis and between the ends, the bound is the rounding alone, and the net is within it.
#[test]
fn a_circle_along_a_line_is_the_exact_cylinder() {
    let mut k = Kernel::default();
    let b = create_pipe_sweep(&mut k, RADIUS, Some(&line_path()), false).unwrap();
    let wall = only_wall(&b);
    let sweep = wall.procedural_sweep().unwrap();
    assert!(sweep.bound() < 1e-9, "bound {}", sweep.bound());
    let (a, e) = (Vec3::new(1.0, 2.0, 3.0), Vec3::new(4.0, 6.0, 15.0));
    let axis = (e - a) * (1.0 / (e - a).length());
    let dom = sweep.net().natural_uv_domain();
    for i in 0..37 {
        for j in 0..29 {
            let uv = Point2::new(
                dom.u_interval().min()
                    + (dom.u_interval().max() - dom.u_interval().min()) * (f64::from(i) + 0.3)
                        / 37.0,
                dom.v_interval().min()
                    + (dom.v_interval().max() - dom.v_interval().min()) * (f64::from(j) + 0.4)
                        / 29.0,
            );
            let p = sweep.definition().evaluate_definition_point(uv).unwrap();
            let rel = p - a;
            let along = rel.dot(axis);
            let radial = (rel - axis * along).length();
            assert!((radial - RADIUS).abs() <= 1e-12, "{radial}");
            assert!(
                along > -1e-12 && along < (e - a).length() + 1e-12,
                "{along}"
            );
        }
    }
    assert!(deviation(&sweep, 29) <= sweep.bound());
}

/// A circle along an arc is the exact torus: the definition's points are `RADIUS` from the
/// arc's circle, the net is within the certified bound, and the bound is no looser than a
/// few times the net's true distance.
#[test]
fn a_circle_along_an_arc_is_the_exact_torus() {
    let mut k = Kernel::default();
    let b = create_pipe_sweep(&mut k, RADIUS, Some(&arc_path()), false).unwrap();
    let wall = only_wall(&b);
    let sweep = wall.procedural_sweep().unwrap();
    let dom = sweep.net().natural_uv_domain();
    for i in 0..61 {
        for j in 0..23 {
            let uv = Point2::new(
                dom.u_interval().min()
                    + (dom.u_interval().max() - dom.u_interval().min()) * (f64::from(i) + 0.3)
                        / 61.0,
                dom.v_interval().min()
                    + (dom.v_interval().max() - dom.v_interval().min()) * (f64::from(j) + 0.4)
                        / 23.0,
            );
            let p = sweep.definition().evaluate_definition_point(uv).unwrap();
            let ring = p.x.hypot(p.y);
            let minor = (ring - 10.0).hypot(p.z);
            assert!((minor - RADIUS).abs() <= 1e-12, "{minor}");
            assert!(p.x >= -1e-12 && p.y >= -1e-12, "{p:?}");
        }
    }
    let far = deviation(&sweep, 2001);
    assert!(far <= sweep.bound(), "{far} past {}", sweep.bound());
    assert!(sweep.bound() <= 4.0 * far, "{} loose, {far}", sweep.bound());
}

/// A circle along a straight and then an arc: two pieces joined, one wall, whose definition
/// carries the piece's parameters.
#[test]
fn a_circle_along_a_line_then_an_arc_has_one_definition() {
    let mut k = Kernel::default();
    let (circle, path) = line_then_arc(&k);
    let sweep = create_curve_sweep(
        &mut k,
        &[Some(&circle)],
        Some(&path),
        None,
        None,
        Some(&pipe_options()),
    )
    .unwrap();
    let all = walls(&sweep.brep);
    assert!(!all.is_empty());
    for wall in &all {
        let s = wall.procedural_sweep().expect("a definition");
        let far = deviation(&s, 3001);
        assert!(far <= s.bound(), "{far} past {}", s.bound());
        let dom = s.net().natural_uv_domain();
        for i in 0..40 {
            for j in 0..12 {
                let uv = Point2::new(
                    dom.u_interval().min()
                        + (dom.u_interval().max() - dom.u_interval().min()) * (f64::from(i) + 0.3)
                            / 40.0,
                    dom.v_interval().min()
                        + (dom.v_interval().max() - dom.v_interval().min()) * (f64::from(j) + 0.4)
                            / 12.0,
                );
                let gap = on_tube(s.definition(), uv);
                assert!((gap - RADIUS).abs() <= 1e-9, "{gap}");
            }
        }
    }
}

/// The balance wire: a circle along a fitted 3D spline of straights and R2 rounds reports
/// its definition and its bound; the definition's points are the circle of `RADIUS` about
/// the path's point, and the net lies within the bound of them off its knots.
#[test]
fn a_circle_sweeps_along_a_fitted_3d_spline_with_straights_and_rounds() {
    let mut k = Kernel::default();
    let b = create_pipe_sweep(&mut k, RADIUS, Some(&spline_path()), false).unwrap();
    let wall = only_wall(&b);
    let sweep = wall.procedural_sweep().unwrap();
    let bound = sweep.bound();
    assert!(bound.is_finite() && bound > 0.0);
    let d = sweep.definition();
    assert!(d.path().is_some() && d.screw().is_none());
    assert!(d.path().unwrap().chain().len() > 60);
    let dom = sweep.net().natural_uv_domain();
    for i in 0..30 {
        for j in 0..500 {
            let uv = Point2::new(
                dom.u_interval().min()
                    + (dom.u_interval().max() - dom.u_interval().min()) * (f64::from(i) + 0.3)
                        / 30.0,
                dom.v_interval().min()
                    + (dom.v_interval().max() - dom.v_interval().min()) * (f64::from(j) + 0.4)
                        / 500.0,
            );
            let gap = on_tube(d, uv);
            assert!((gap - RADIUS).abs() <= 1e-9, "{gap}");
        }
    }
    let far = deviation(&sweep, 6001);
    assert!(far <= bound, "{far} past {bound}");
}

/// The definition follows its net: a rigid motion moves both, a reversal, a swap and a
/// reparametrisation move its map, a piece keeps it.
#[test]
fn the_path_definition_follows_its_nets_edits() {
    let mut k = Kernel::default();
    let b = create_pipe_sweep(&mut k, RADIUS, Some(&arc_path()), false).unwrap();
    let wall = only_wall(&b);
    let base = deviation(&wall.procedural_sweep().unwrap(), 801);
    assert!(base <= wall.procedural_sweep().unwrap().bound());
    let place = Axis2Placement::new(
        Point3::new(10.0, -3.0, 7.0),
        Vec3::new(0.6, 0.8, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    );
    let mut moved = wall.copy(&k).unwrap();
    moved.transform(&k, &place, None).unwrap();
    let s = moved.procedural_sweep().expect("moved with its net");
    let far = deviation(&s, 801);
    assert!(far <= s.bound(), "{far} past {}", s.bound());
    assert!((far - base).abs() <= 1e-6 * (1.0 + base));
    for edit in 0..4 {
        let mut e = wall.copy(&k).unwrap();
        match edit {
            0 => e.reverse(&k, SurfParamType::U).unwrap(),
            1 => e.reverse(&k, SurfParamType::V).unwrap(),
            2 => e.swap_uv(&k).unwrap(),
            _ => {
                let dom = solidcore::geom::extent2d::Extent2d::from_bounds(-1.0, 3.0, 2.0, 9.0);
                e.reparameterize(&k, &dom).unwrap();
            }
        }
        let s = e.procedural_sweep().expect("edited with its net");
        let far = deviation(&s, 801);
        assert!(far <= s.bound(), "edit {edit}: {far} past {}", s.bound());
    }
    // a scale that is not uniform leaves no definition, and the net moves without it
    let mut squashed = wall.copy(&k).unwrap();
    squashed
        .transform(
            &k,
            &Axis2Placement::default(),
            Some(Vec3::new(1.0, 1.0, 2.0)),
        )
        .unwrap();
    assert!(squashed.procedural_sweep().is_none());
    assert!(matches!(squashed.kind(), SurfaceKind::BSpline(_)));
}

/// The ledger: each kind of wall face that is not a definition yet, with its issue (the
/// first `#n` of its file); only the `path-` kinds are this slice's.
fn path_ledger() -> BTreeMap<String, u32> {
    let dir = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/data/sweep_definition_ledger"
    );
    let mut out = BTreeMap::new();
    for entry in std::fs::read_dir(dir).unwrap() {
        let path = entry.unwrap().path();
        let stem = path.file_stem().unwrap().to_string_lossy().into_owned();
        if !stem.starts_with("path-") {
            continue;
        }
        let text = std::fs::read_to_string(&path).unwrap();
        let issue = text
            .split('#')
            .nth(1)
            .and_then(|t| {
                t.chars()
                    .take_while(char::is_ascii_digit)
                    .collect::<String>()
                    .parse()
                    .ok()
            })
            .unwrap_or_else(|| panic!("{stem}: the ledger entry names no issue"));
        out.insert(stem, issue);
    }
    out
}

/// Every wall face of every path sweep case is a definition, or of a kind the ledger lists;
/// and every listed kind still appears.
#[test]
fn issue_2572_every_path_wall_is_a_definition_or_in_the_ledger() {
    let mut k = Kernel::default();
    let (circle, path) = line_then_arc(&k);
    let translational = SweepOptions {
        translational_sweep: true,
        ..pipe_options()
    };
    let reference = Curve::from(
        Line::create_line_segment(3, Vec3::new(0.0, 0.0, 0.0), Vec3::new(10.0, 0.0, 0.0), None)
            .unwrap(),
    );
    let scale = Curve::from(
        Line::create_line_segment(3, Vec3::new(0.0, 1.0, 0.0), Vec3::new(10.0, 2.0, 0.0), None)
            .unwrap(),
    );
    let closed = Curve::from(create_arc(Point3::new(0.0, 0.0, 0.0), 10.0, 0.0, 360.0).unwrap());
    // (the kind its bare walls are listed as, the body)
    let cases: Vec<(&str, Brep)> = vec![
        (
            "path-b-spline",
            create_pipe_sweep(&mut k, RADIUS, Some(&line_path()), false).unwrap(),
        ),
        (
            "path-b-spline",
            create_pipe_sweep(&mut k, RADIUS, Some(&arc_path()), false).unwrap(),
        ),
        (
            "path-b-spline",
            create_pipe_sweep(&mut k, RADIUS, Some(&spline_path()), false).unwrap(),
        ),
        (
            "path-b-spline",
            create_curve_sweep(
                &mut k,
                &[Some(&circle)],
                Some(&path),
                None,
                None,
                Some(&pipe_options()),
            )
            .unwrap()
            .brep,
        ),
        (
            "path-translational",
            create_curve_sweep(
                &mut k,
                &[Some(&circle)],
                Some(&path),
                None,
                None,
                Some(&translational),
            )
            .unwrap()
            .brep,
        ),
        (
            "path-scaled",
            create_curve_sweep(
                &mut k,
                &[Some(&circle)],
                Some(&reference),
                Some(&reference),
                Some(&scale),
                None,
            )
            .unwrap()
            .brep,
        ),
        (
            "path-closed",
            create_pipe_sweep(&mut k, RADIUS, Some(&closed), false).unwrap(),
        ),
    ];
    let ledger = path_ledger();
    let mut bare: BTreeSet<String> = BTreeSet::new();
    let mut definitions = 0;
    for (kind, b) in &cases {
        for s in walls(b) {
            if s.procedural_definition().is_some() {
                definitions += 1;
                continue;
            }
            if matches!(s.kind(), SurfaceKind::Plane(_)) {
                continue;
            }
            bare.insert((*kind).to_owned());
        }
    }
    assert!(definitions >= 4, "{definitions}");
    let unlisted: Vec<&String> = bare.iter().filter(|kd| !ledger.contains_key(*kd)).collect();
    assert!(
        unlisted.is_empty(),
        "bare wall faces of kinds the ledger does not list: {unlisted:?}"
    );
    let fixed: Vec<&String> = ledger.keys().filter(|kd| !bare.contains(*kd)).collect();
    assert!(
        fixed.is_empty(),
        "listed kinds no case makes bare any more (fixed? delete their ledger files): {fixed:?}"
    );
}
