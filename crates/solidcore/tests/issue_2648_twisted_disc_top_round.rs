//! #2648 (solidcraft#11848): a 2 mm round of a twisted slotted disc's top outline spent
//! most of its time in `band_cuts`, the check that a convex blend keeps the part whole
//! (`fixes.fillet_gn_section_generator`): every section of every blended edge intersected
//! every face (building the face's trimmed cache again each time) and every edge.
//!
//! - The part, rebuilt from the issue's script (the corpus entry
//!   `kernel_cases/fillets/11848-twisted-slotted-disc-top-outline`): a disc of radius 30
//!   with 2 V-slots, swept up 15 with a 30 degree twist (cubic walls over the entry's 28
//!   spans), its 4 twisted side edges rounded 2 mm. It has the entry's 14 faces, 36 edges and
//!   24 vertices; the slot-bottom vertices on the top, exact images of the profile's, are the
//!   entry's to 1e-9, and the side rounds' ends lie within 0.02 of the entry's (SolidCraft's
//!   own sweep and round are not this test's).
//! - The check judges each section against the geometry within its reach and makes each
//!   face's cache once per check: counted on the WR-696 bridge (its refusal and fit size
//!   unchanged, `wr696_band_cuts.rs`).
//! - The top outline's round: a probe (`--ignored`, `SLOTS=n`), timed in the pull request.
//!
//! Fix set: the corpus's `fixes.json` (SolidCraft's `fixes()`; the entry's own list adds
//! only `fillet_through_extended_rail_end`, a flag removed in #2500).
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::approx_constant // the entry's vertices, as printed
)] // test helpers

use solidcore::blend::fillet_executive::api::{create_face_from_surface, fillet_edges};
use solidcore::curve::curve::Curve;
use solidcore::error::EntityRef;
use solidcore::geom::tol;
use solidcore::kernel::UpstreamFixes;
use solidcore::math::Point3;
use solidcore::nurbs::b_spline_curve::{BSplineCurve, BSplineCurveForm, CanonicalData, KnotType};
use solidcore::nurbs::b_spline_surface::{
    BSplineSurface, BSplineSurfaceForm, CanonicalSurfaceData,
};
use solidcore::surface::surface::Surface;
use solidcore::topology::brep::p5::{stitch_and_orient, StitchAndOrientOptions};
use solidcore::topology::brep::p8::{create_planar_faces_with_3d_curves, FaceBuildEnv};
use solidcore::topology::{edge, vertex, EdgeId, FaceId};
use solidcore::{Brep, Kernel};

/// The disc's radius, height and twist (degrees over the height).
const RADIUS: f64 = 30.0;
const HEIGHT: f64 = 15.0;
const TWIST: f64 = 30.0;

/// The corpus's fix list (`tests/kernel_cases/fixes.json`, SolidCraft's `fixes()`).
fn solidcraft_kernel() -> Kernel {
    let names: Vec<String> = serde_json::from_str(include_str!("kernel_cases/fixes.json")).unwrap();
    let mut f = UpstreamFixes::default();
    for name in &names {
        if let Some(flag) = f.flag_mut(name) {
            *flag = true;
        }
    }
    assert!(f.fillet_gn_section_generator);
    Kernel {
        fixes: f,
        ..Kernel::default()
    }
}

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

fn rotate(p: Point3, a: f64) -> Point3 {
    let (s, c) = a.sin_cos();
    Point3::new(c * p.x - s * p.y, s * p.x + c * p.y, p.z)
}

/// The profile of a disc of radius 30 with `n` V-slots (walls at 45 degrees to the slot's
/// axis, a 2 mm arc at the bottom 2.17 mm into the disc), the first slot's axis along +y.
fn profile(n: usize) -> Vec<Piece> {
    let s2 = std::f64::consts::SQRT_2;
    let centre_y = 27.828_427_124_746_188;
    let a = Point3::new(-4.639_136_501_002_611_5, 29.639_136_501_002_61, 0.0);
    let b = Point3::new(-s2, centre_y - s2, 0.0);
    let c = Point3::new(s2, centre_y - s2, 0.0);
    let d = Point3::new(4.639_136_501_002_611_5, 29.639_136_501_002_61, 0.0);
    let step = std::f64::consts::TAU / n as f64;
    let slots: Vec<[Point3; 4]> = (0..n)
        .map(|k| {
            let t = step * k as f64;
            [rotate(a, t), rotate(b, t), rotate(c, t), rotate(d, t)]
        })
        .collect();
    let q = std::f64::consts::FRAC_PI_4;
    let mut pieces = Vec::new();
    for k in 0..n {
        let t = step * k as f64;
        let [a, b, c, d] = slots[k];
        let centre = rotate(Point3::new(0.0, centre_y, 0.0), t);
        pieces.push(line(a, b));
        pieces.push(arc(centre, 2.0, t - 3.0 * q, t - q, b, c));
        pieces.push(line(c, d));
        // the outer arc from this slot's right wall to the next slot's left wall (clockwise)
        let next = slots[(k + n - 1) % n][0];
        let ad = d.y.atan2(d.x);
        let an = next.y.atan2(next.x);
        let mut an = an;
        while an > ad {
            an -= std::f64::consts::TAU;
        }
        pieces.push(arc(Point3::ZERO, RADIUS, ad, an, d, next));
    }
    pieces
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

/// The twist's v knots, as SolidCraft's sweep makes them (the corpus entry's walls): cubic,
/// 28 spans, the interior knots `k / 30` for `k = 2 .. 28`.
fn twist_knots() -> (Vec<f64>, Vec<usize>) {
    let mut knots = vec![0.0];
    let mut mults = vec![4];
    for k in 2..=28 {
        knots.push(f64::from(k) / 30.0);
        mults.push(1);
    }
    knots.push(1.0);
    mults.push(4);
    (knots, mults)
}

/// The rows of the twist in v: each row the profile turned by the twist at its Greville
/// abscissa and lifted to its height (z exactly linear in v).
fn twist_rows() -> Vec<(f64, f64)> {
    let (knots, mults) = twist_knots();
    let mut flat = Vec::new();
    for (&k, &m) in knots.iter().zip(&mults) {
        flat.extend(std::iter::repeat_n(k, m));
    }
    (0..flat.len() - 4)
        .map(|j| {
            let g = (flat[j + 1] + flat[j + 2] + flat[j + 3]) / 3.0;
            (TWIST.to_radians() * g, HEIGHT * g)
        })
        .collect()
}

fn place(p: Point3, (angle, z): (f64, f64)) -> Point3 {
    let q = rotate(p, angle);
    Point3::new(q.x, q.y, z)
}

fn wall(p: &Piece) -> Surface {
    let rows = twist_rows();
    let (vk, vm) = twist_knots();
    // u along the profile, v up: control points u-major
    let mut pts = Vec::new();
    let mut weights = Vec::new();
    for (i, &q) in p.points.iter().enumerate() {
        for &row in &rows {
            pts.push(place(q, row));
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
            v_multiplicities: &vm,
            u_knots: &p.knots,
            v_knots: &vk,
            knot_type: KnotType::Unspecified,
            weights: p.weights.as_ref().map(|_| weights.as_slice()),
            uv_domain: None,
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

/// The disc with `n` slots swept up 15 with a 30 degree twist.
fn twisted_disc(k: &mut Kernel, n: usize) -> Brep {
    let pieces = profile(n);
    let mut b =
        solidcore::topology::brep::p3::new_brep(k, tol::USE_DEFAULT, tol::USE_DEFAULT).unwrap();
    for p in &pieces {
        create_face_from_surface(k, &mut b, &wall(p)).unwrap();
    }
    let rows = twist_rows();
    for row in [rows[0], rows[rows.len() - 1]] {
        let cap: Vec<Curve> = pieces
            .iter()
            .map(|p| curve_of(p, |q| place(q, row)))
            .collect();
        with_face_env(k, &mut b, |k, b, fenv| {
            let tol = b.tolerance();
            let region = b.infinite_region();
            create_planar_faces_with_3d_curves(k, b, region, &cap, tol, fenv).unwrap();
        });
    }
    stitch_and_orient(k, &mut b, &StitchAndOrientOptions::default()).unwrap();
    b.refresh_cached_topology();
    assert!(solidcore::topology::brep::p2::is_manifold_solid(&b).unwrap());
    b
}

fn ends(b: &Brep, e: EdgeId) -> [Point3; 2] {
    let v0 = edge::vertex(b, e).unwrap().unwrap();
    let v1 = edge::other_vertex(b, e, Some(v0)).unwrap().unwrap();
    [vertex::point(b, v0).unwrap(), vertex::point(b, v1).unwrap()]
}

/// The side edges where the outer walls meet the slots' walls: the edges running from the
/// bottom to the top whose bottom end is at the disc's radius.
fn side_edges(b: &Brep) -> Vec<EdgeId> {
    b.edges()
        .ids()
        .filter(|&e| {
            let [p, q] = ends(b, e);
            let (lo, hi) = if p.z < q.z { (p, q) } else { (q, p) };
            lo.z.abs() < 1e-6
                && (hi.z - HEIGHT).abs() < 1e-6
                && ((lo.x * lo.x + lo.y * lo.y).sqrt() - RADIUS).abs() < 1e-3
        })
        .collect()
}

/// The top outline: every edge with both ends on the top.
fn top_edges(b: &Brep) -> Vec<EdgeId> {
    b.edges()
        .ids()
        .filter(|&e| ends(b, e).iter().all(|p| (p.z - HEIGHT).abs() < 1e-6))
        .collect()
}

fn side_rounded_disc(k: &mut Kernel, n: usize) -> Brep {
    let mut b = twisted_disc(k, n);
    let sides = side_edges(&b);
    assert_eq!(sides.len(), 2 * n);
    fillet_edges(k, &mut b, &sides, 2.0, 1, 1, 1.0).unwrap();
    b.refresh_cached_topology();
    b
}

#[test]
fn the_rebuilt_disc_is_the_corpus_entry_s() {
    let mut k = solidcraft_kernel();
    let b = side_rounded_disc(&mut k, 2);
    assert_eq!(b.faces().ids().count(), 14);
    assert_eq!(b.edges().ids().count(), 36);
    let points: Vec<Point3> = b
        .vertices()
        .ids()
        .map(|v| vertex::point(&b, v).unwrap())
        .collect();
    assert_eq!(points.len(), 24);
    // the entry's vertices (`11848-twisted-slotted-disc-top-outline.bodies.json.gz`, b7)
    let entry: [[f64; 3]; 24] = [
        [-6.421243, 29.304738, 0.0],
        [-6.399435, -29.309508, 0.0],
        [-3.71353, 28.71353, 0.0],
        [-3.690963, -28.690963, 0.0],
        [-1.414214, -26.414214, 0.0],
        [-1.414214, 26.414214, 0.0],
        [1.414214, -26.414214, 0.0],
        [1.414214, 26.414214, 0.0],
        [3.690963, 28.690963, 0.0],
        [3.71353, -28.71353, 0.0],
        [6.399435, 29.309508, 0.0],
        [6.421243, -29.304738, 0.0],
        [-20.196827, 22.183061, 15.0],
        [-17.54195, 23.001621, 15.0],
        [-14.431852, 22.168273, 15.0],
        [-11.982362, 23.582487, 15.0],
        [-11.140754, 26.723411, 15.0],
        [-9.09141, 28.589268, 15.0],
        [9.09141, -28.589268, 15.0],
        [11.140754, -26.723411, 15.0],
        [11.982362, -23.582487, 15.0],
        [14.431852, -22.168273, 15.0],
        [17.54195, -23.001621, 15.0],
        [20.196827, -22.183061, 15.0],
    ];
    for e in entry {
        let e = Point3::new(e[0], e[1], e[2]);
        let near = points
            .iter()
            .map(|p| p.distance(e))
            .fold(f64::INFINITY, f64::min);
        // the slot bottoms' vertices (|x| < 2 at the bottom, their images turned 30 degrees
        // on the top) are the profile's exactly (the entry's to its 6 printed decimals)
        let exact = (e.z == 0.0 && e.x.abs() < 2.0)
            || (e.z == HEIGHT && e.x.abs() < 15.0 && e.x.abs() > 11.5);
        let within = if exact { 1e-6 } else { 0.02 };
        assert!(near <= within, "entry vertex {e:?}: nearest {near}");
    }
}

/// The WR-696 bridge (`wr696_band_cuts.rs`): a round of 9 would take all of the bridge, refused
/// naming 7.828427 as before. The check judges its sections against the geometry within
/// their reach and makes each face's cache once: the counts of its work.
#[test]
fn the_band_check_tests_only_what_its_sections_reach() {
    use solidcore::blend::fillet_executive::session::FilletSession;
    use solidcore::boolean::api_brep::boolean_difference;
    use solidcore::primitive::api_primitives::create_box;
    let mut k = Kernel::default();
    let mut fixes = UpstreamFixes::all();
    fixes.fillet_round_ends_on_outline = false;
    k.fixes = fixes;
    let block = create_box(&mut k, Point3::ZERO, 20.0, 4.0, 10.0).unwrap();
    let slot = create_box(&mut k, Point3::new(8.0, -1.0, -1.0), 4.0, 6.0, 10.0).unwrap();
    let mut b = boolean_difference(&mut k, block, slot).unwrap();
    let faces = b.faces().ids().count() as u64;
    let edges = b.edges().ids().count() as u64;
    let e = b
        .edges()
        .ids()
        .find(|&e| {
            let [p, q] = ends(&b, e);
            p.y.abs() < 1e-9
                && q.y.abs() < 1e-9
                && (p.z - 10.0).abs() < 1e-9
                && (q.z - 10.0).abs() < 1e-9
        })
        .unwrap();
    k.work_counters.reset();
    let mut session = FilletSession::new(&mut k, &mut b).unwrap();
    session
        .prepare_constant_radius_edges(&mut k, &[e], 9.0, 1, 2, 0.0)
        .unwrap();
    assert!(session.do_fillet(&mut k).is_err());
    assert_eq!(
        session.error_info().comment,
        "Fillet Error: the round would take all of the part between its rails along the edge and cut it into pieces: a radius below 7.828427 keeps it whole"
    );
    let w = &k.work_counters;
    let sections = w.band_sections();
    eprintln!(
        "sections {sections}, face solves {}, face caches {}, edge solves {}; {faces} faces, {edges} edges",
        w.band_face_solves(),
        w.band_face_caches(),
        w.band_edge_solves()
    );
    // before #2648 each of the 278 sections intersected all 10 faces, making each face's
    // cache again (2780 of each), and each of the 32 stations all 24 edges (768)
    assert_eq!((sections, faces, edges), (278, 10, 24));
    assert_eq!(w.band_face_solves(), 510);
    assert_eq!(w.band_face_caches(), 4);
    assert_eq!(w.band_edge_solves(), 32);
}

/// The top outline's round, timed (`SLOTS=n`, default 2).
#[test]
#[ignore = "timing probe"]
fn probe_top_outline_round() {
    let n: usize = std::env::var("SLOTS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(2);
    let mut k = solidcraft_kernel();
    let mut b = side_rounded_disc(&mut k, n);
    let top = top_edges(&b);
    k.work_counters.reset();
    let t = std::time::Instant::now();
    let r = fillet_edges(&mut k, &mut b, &top, 2.0, 1, 1, 1.0);
    let w = &k.work_counters;
    eprintln!(
        "top round of {} edges: {:?} -> {:?}; band sections {}, face solves {}, face caches {}, edge solves {}",
        top.len(),
        t.elapsed(),
        r.as_ref().err(),
        w.band_sections(),
        w.band_face_solves(),
        w.band_face_caches(),
        w.band_edge_solves()
    );
}
