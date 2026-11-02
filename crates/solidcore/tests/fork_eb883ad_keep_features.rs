//! The SolidCraft fork's tangency-labs/smlib@eb883ad (WR-703, #2273): a boss or a hole inside a
//! round's band stays, the round trimmed round it in place; the cases are the fork's ddbb782.
//!
//! The 40 x 30 x 20 block's top front edge (y = 0, z = 20) is rounded; the feature, a cylinder of
//! radius 2 about the upright x = 20, y = 4, stands 2 to 6 in from the front on the top. A round
//! of radius r takes the band y < r off the top and z > 20 - r off the front, so at r 8 the
//! footprint lies inside the band, at r 1.9 the band stops short of it.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)] // test helpers

use solidcore::blend::fillet_executive::api::fillet_edges;
use solidcore::boolean::api_brep::{boolean_difference, boolean_union};
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::kernel::{KernelOptions, UpstreamFixes};
use solidcore::mass::integrator::properties::{brep_volume, PropertyAttributes};
use solidcore::math::{Point3, PI};
use solidcore::primitive::api_primitives::{create_box, create_cone, create_cylinder};
use solidcore::query::api_queries::edge_closest_point;
use solidcore::topology::brep::validation::{check_brep, CheckLevel};
use solidcore::topology::{edge, EdgeId};
use solidcore::{Brep, Kernel, Result};

fn kernel(keep: bool) -> Kernel {
    let mut k = Kernel::default();
    k.fixes.fillet_keep_features = keep;
    k
}

/// The block with a boss 5 tall on its top, or a hole 5 deep in it.
fn featured_block(k: &mut Kernel, boss: bool) -> Brep {
    let block = create_box(k, Point3::new(0.0, 0.0, 0.0), 40.0, 30.0, 20.0).unwrap();
    if boss {
        let c = create_cylinder(k, Point3::new(20.0, 4.0, 20.0), 2.0, 5.0).unwrap();
        boolean_union(k, block, c).unwrap()
    } else {
        let c = create_cylinder(k, Point3::new(20.0, 4.0, 15.0), 2.0, 6.0).unwrap();
        boolean_difference(k, block, c).unwrap()
    }
}

/// The edge nearest the point, within 0.1.
fn pick(k: &Kernel, b: &Brep, q: Point3) -> EdgeId {
    let mut best = (0.1f64, None);
    for &e in b.edge_list() {
        if let Ok(c) = edge_closest_point(k, b, e, q) {
            if c.closest.distance <= best.0 {
                best = (c.closest.distance, Some(e));
            }
        }
    }
    best.1.expect("an edge there")
}

fn volume_of(k: &Kernel, b: &mut Brep) -> f64 {
    let store = AttributeStore::default();
    let attrs = PropertyAttributes {
        store: &store,
        key_of: &|_| ObjectKey(0),
    };
    brep_volume(k, b, 1e-7, &attrs).unwrap()
}

/// The material the round r 8 takes from under the footprint, down from the top to the round's
/// surface z = 12 + sqrt(64 - (8 - y)^2): over the disc, y = 4 + 2 sin t, its chord 4 cos t,
/// Int 8 cos^2 t (8 - sqrt(64 - (4 - 2 sin t)^2)) dt over t in [-pi/2, pi/2], by Simpson's rule.
fn column_under_footprint() -> f64 {
    let spans = 4000;
    let h = PI / f64::from(spans);
    let f = |t: f64| {
        let (c, d) = (t.cos(), 4.0 - 2.0 * t.sin());
        8.0 * c * c * (8.0 - (64.0 - d * d).sqrt())
    };
    let mut sum = f(-0.5 * PI) + f(0.5 * PI);
    for i in 1..spans {
        sum += f(-0.5 * PI + h * f64::from(i)) * if i % 2 == 1 { 4.0 } else { 2.0 };
    }
    sum * h / 3.0
}

/// The round's cross-section area taken off, for its 40 long edge.
fn spandrel(r: f64) -> f64 {
    (1.0 - PI / 4.0) * r * r
}

/// The block less its round r, with the boss's 20 pi on it or the hole's 20 pi out of it.
fn featured_rounded(boss: bool, r: f64) -> f64 {
    24000.0 - 40.0 * spandrel(r) + if boss { 20.0 * PI } else { -20.0 * PI }
}

/// How far round the feature's axis the edges lying on both its wall and the round r 8's surface
/// run, in radians: 2 pi where the wall comes down to the round all round its footprint.
fn wall_on_round(k: &Kernel, b: &Brep) -> f64 {
    let tol = 100.0 * b.tolerance();
    let mut turned = 0.0;
    for &e in b.edge_list() {
        let Some(curve) = edge::curve(b, e).unwrap() else {
            continue;
        };
        let ivl = edge::interval(b, e).unwrap();
        let (mut span, mut last) = (0.0, 0.0);
        let mut on = true;
        for s in 0..=256 {
            let p = curve
                .evaluate_point(k, ivl.evaluate(f64::from(s) / 256.0))
                .unwrap();
            on = ((p.x - 20.0).hypot(p.y - 4.0) - 2.0).abs() < tol
                && ((p.y - 8.0).hypot(p.z - 12.0) - 8.0).abs() < tol;
            if !on {
                break;
            }
            let a = (p.y - 4.0).atan2(p.x - 20.0);
            if s > 0 {
                let d: f64 = a - last;
                span += d - (2.0 * PI) * (d / (2.0 * PI)).round();
            }
            last = a;
        }
        if on {
            turned += f64::abs(span);
        }
    }
    turned
}

/// One sound solid of volume `want` within `rel` with `faces` faces; with `reach`, the wall on
/// the round all round its footprint.
fn assert_feature_in_band(
    boss: bool,
    keep: bool,
    r: f64,
    want: f64,
    rel: f64,
    faces: usize,
    reach: bool,
) {
    let mut k = kernel(keep);
    let mut b = featured_block(&mut k, boss);
    let e = pick(&k, &b, Point3::new(20.0, 0.0, 20.0));
    fillet_edges(&mut k, &mut b, &[e], r, 1, 1, 1.0).unwrap();
    let v = volume_of(&k, &mut b);
    let turned = wall_on_round(&k, &b);
    println!(
        "boss {boss} keep {keep} r {r}: volume {v:.9}, want {want:.9}; {} faces, want {faces}; wall on the round {turned:.9}",
        b.face_list().len()
    );
    assert!((v - want).abs() <= rel * want, "volume {v}, want {want}");
    assert_eq!(b.face_list().len(), faces);
    if reach {
        assert!(
            (turned - 2.0 * PI).abs() < 1e-9,
            "wall on the round {turned}"
        );
    } else {
        assert!(turned.abs() < 1e-12, "wall on the round {turned}");
    }
    assert!(solidcore::topology::brep::p2::is_manifold_solid(&b).unwrap());
    let lines = check_brep(&mut k, &mut b, CheckLevel::Level2)
        .unwrap()
        .lines();
    assert!(lines.is_empty(), "{lines:?}");
}

/// To 1e-7, as the fork's own case at 5e8c9c1 (its 1e-9 at ddbb782 went to 1e-7 in the same PR):
/// the crease is fitted within a tenth of the fillet's tolerance, which moves the column's volume
/// by about ten times that.
#[test]
fn the_boss_in_the_band_is_kept() {
    let want = featured_rounded(true, 8.0) + column_under_footprint();
    assert_feature_in_band(true, true, 8.0, want, 1e-7, 9, true);
}

#[test]
fn the_hole_in_the_band_is_kept() {
    let want = featured_rounded(false, 8.0) + column_under_footprint();
    assert_feature_in_band(false, true, 8.0, want, 1e-7, 9, true);
}

#[test]
fn with_keep_off_the_round_takes_the_boss() {
    assert_feature_in_band(
        true,
        false,
        8.0,
        24000.0 - 40.0 * spandrel(8.0),
        1e-9,
        7,
        false,
    );
}

#[test]
fn with_keep_off_the_round_takes_the_hole() {
    assert_feature_in_band(
        false,
        false,
        8.0,
        24000.0 - 40.0 * spandrel(8.0),
        1e-9,
        7,
        false,
    );
}

#[test]
fn a_band_short_of_the_feature_leaves_it_either_way() {
    for keep in [false, true] {
        for boss in [true, false] {
            assert_feature_in_band(boss, keep, 1.9, featured_rounded(boss, 1.9), 1e-9, 9, false);
        }
    }
}

// ---- Beyond the fork's own cases: every flag on, the figures from the fork's kernel at eb883ad ----

/// The block with a feature on or in its top.
fn block_with(k: &mut Kernel, feature: fn(&mut Kernel) -> Result<Brep>, boss: bool) -> Brep {
    let block = create_box(k, Point3::new(0.0, 0.0, 0.0), 40.0, 30.0, 20.0).unwrap();
    let cutter = feature(k).unwrap();
    if boss {
        boolean_union(k, block, cutter).unwrap()
    } else {
        boolean_difference(k, block, cutter).unwrap()
    }
}

/// The block's edge rounded 8, picked away from the feature: its volume and faces.
fn rounded_with(
    keep: bool,
    feature: fn(&mut Kernel) -> Result<Brep>,
    boss: bool,
) -> Result<(f64, usize)> {
    let mut k = Kernel::new(&KernelOptions::default());
    k.fixes = UpstreamFixes::all();
    k.fixes.fillet_keep_features = keep;
    let mut b = block_with(&mut k, feature, boss);
    let e = pick(&k, &b, Point3::new(5.0, 0.0, 20.0));
    fillet_edges(&mut k, &mut b, &[e], 8.0, 1, 1, 1.0)?;
    assert!(solidcore::topology::brep::p2::is_manifold_solid(&b).unwrap());
    let lines = check_brep(&mut k, &mut b, CheckLevel::Level2)
        .unwrap()
        .lines();
    assert!(lines.is_empty(), "{lines:?}");
    Ok((volume_of(&k, &mut b), b.face_list().len()))
}

fn assert_kept_as_the_fork(
    feature: fn(&mut Kernel) -> Result<Brep>,
    boss: bool,
    volume: f64,
    faces: usize,
) {
    let (v, n) = rounded_with(true, feature, boss).unwrap();
    assert!(
        (v - volume).abs() <= 1e-9 * volume,
        "volume {v}, the fork's {volume}"
    );
    assert_eq!(n, faces);
    let (v, n) = rounded_with(false, feature, boss).unwrap();
    let taken = 24000.0 - 40.0 * spandrel(8.0);
    assert!(
        (v - taken).abs() <= 1e-9 * taken,
        "volume {v}, want {taken}"
    );
    assert_eq!(n, 7);
}

fn prism_boss(k: &mut Kernel) -> Result<Brep> {
    create_box(k, Point3::new(18.0, 2.0, 20.0), 4.0, 4.0, 5.0)
}

fn prism_pocket(k: &mut Kernel) -> Result<Brep> {
    create_box(k, Point3::new(18.0, 2.0, 15.0), 4.0, 4.0, 6.0)
}

fn tapered_boss(k: &mut Kernel) -> Result<Brep> {
    create_cone(k, Point3::new(20.0, 4.0, 20.0), 2.0, 1.0, 5.0)
}

fn tapered_pit(k: &mut Kernel) -> Result<Brep> {
    create_cone(k, Point3::new(20.0, 4.0, 15.0), 1.0, 2.0, 5.5)
}

/// Planar walls: each crease is an exact section, and the volume is the fork's in 17 digits.
#[test]
fn a_prism_s_planar_walls_are_kept_as_the_fork_keeps_them() {
    assert_kept_as_the_fork(prism_boss, true, 23549.896259768986, 12);
    assert_kept_as_the_fork(prism_pocket, false, 23389.896259952176, 12);
}

/// Walls that lean: as the fork, kept, with no refusal for the boss's.
#[test]
fn a_cone_s_walls_are_kept_as_the_fork_keeps_them() {
    assert_kept_as_the_fork(tapered_boss, true, 23501.985588925883, 11);
    assert_kept_as_the_fork(tapered_pit, false, 23427.540889420579, 9);
}

/// A boss astride the round's contact line: refused with the flag on or off, as the fork.
#[test]
fn a_boss_astride_the_contact_line_is_refused_either_way() {
    fn astride(k: &mut Kernel) -> Result<Brep> {
        create_box(k, Point3::new(18.0, 6.0, 20.0), 4.0, 4.0, 5.0)
    }
    for keep in [true, false] {
        assert!(rounded_with(keep, astride, true).is_err(), "keep {keep}");
    }
}

/// A pocket astride the contact line stands either way: the merge keeps it, as the fork.
#[test]
fn a_pocket_astride_the_contact_line_stands_either_way() {
    fn astride(k: &mut Kernel) -> Result<Brep> {
        create_box(k, Point3::new(18.0, 6.0, 15.0), 4.0, 4.0, 6.0)
    }
    for keep in [true, false] {
        let (v, n) = rounded_with(keep, astride, false).unwrap();
        let fork = 23371.292359325398;
        assert!((v - fork).abs() <= 1e-9 * fork, "keep {keep}: volume {v}");
        assert_eq!(n, 12);
    }
}
