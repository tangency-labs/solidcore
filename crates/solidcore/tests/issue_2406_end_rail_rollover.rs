//! #2406 (from tangency-labs/smlib#98): a round ending at a smooth vertex while its rail
//! has already rolled over its side face's tangent edge onto the next face.
//!
//! The part is SolidCraft #851's plate with its concave R21 corner wall (the fork's capi
//! `arc_wall_plate`): the plate x 0..41 1/6, y -21..0, z 0..10 and the corner
//! x 41 1/6..62 1/6 less the cylinder r 21 about the axis A = (62 1/6, -21). Its top front
//! edge (y = -21, z = 10) is rounded R2, then the middle edge of the arc top edge's tangent
//! chain, between that round and the arc wall, is rounded alone at r (the fork's capi
//! `scenario_middle_edge_of_an_arc_tangent_chain_r0_25`). The middle edge ends at two smooth
//! vertices: at the bottom (41 1/6, -21, 8) it runs on into the wall's upright edge, at the
//! top V = (62 1/6 - sqrt(437), -19, 10) into the arc's top edge. Each End cuts the round on
//! the plane through its vertex square to its edge: z = 8 at the bottom, the plane P
//! through V and the wall's axis at the top.
//!
//! The rail on the R2 round reaches the round's tangent edge T with the top face
//! (y = -19, z = 10) before P: a ball touching the round on T at x has its centre at
//! (x, -19, 10 - r), 21 + r from A, so (x - 62 1/6)^2 + 4 = (21 + r)^2, at the angle
//! asin(2 / (21 + r)) about A, short of V's asin(2 / 21) for every r > 0. Past T the ball
//! rolls on the top face: from there to P the round is the one between the top face and
//! the wall (centres on the circle of radius 21 + r about A at z = 10 - r), its rail on the
//! top face the circle of radius 21 + r at z = 10, and the cap's edge on the top face is
//! the top face's line in P from that rail to V. The fork (18e96f1) refuses it,
//! `SMLIB_ERR_NO_ROOM` (`CalcEndCapEdges`: "a rail rolling over into an end is not built").
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers
use solidcore::blend::fillet_executive::session::FilletSession;
use solidcore::blend::fillet_solver::p1::FilletSurfaceGeneratorType;
use solidcore::boolean::api_brep::{boolean_difference, boolean_union};
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::geom::global_solver::SolverOperation;
use solidcore::kernel::UpstreamFixes;
use solidcore::mass::integrator::properties::{brep_volume, PropertyAttributes};
use solidcore::math::{Point3, Vec3, BIG_DOUBLE};
use solidcore::nurbs::knots::Side;
use solidcore::primitive::api_primitives::{create_box, create_cylinder};
use solidcore::query::api_queries::{get_edges, get_faces};
use solidcore::surface::surface::p4::NoObject;
use solidcore::surface::surface::p5::DropPointOptions;
use solidcore::topology::brep::p5::{stitch_and_orient, StitchAndOrientOptions};
use solidcore::topology::brep::validation::{check_brep, CheckLevel};
use solidcore::topology::face::p2 as face2;
use solidcore::topology::{edge, vertex, EdgeId, FaceId};
use solidcore::{Brep, Kernel};

/// Every fix but `geometric_g2_continuity`, with or without `fillet_blend_corner`.
fn all_fixes(blend: bool) -> Kernel {
    let mut k = Kernel::default();
    let mut f = UpstreamFixes::default();
    for &name in UpstreamFixes::NAMES {
        if name != "geometric_g2_continuity" && (blend || name != "fillet_blend_corner") {
            *f.flag_mut(name).unwrap() = true;
        }
    }
    k.fixes = f;
    k
}

fn pick(k: &Kernel, b: &Brep, at: [f64; 3]) -> EdgeId {
    let p = Point3::new(at[0], at[1], at[2]);
    let mut best: Option<(EdgeId, f64)> = None;
    for e in get_edges(b).unwrap() {
        let Some(c) = edge::curve(b, e).unwrap() else {
            continue;
        };
        let ivl = edge::interval(b, e).unwrap();
        let d = c
            .drop_point(
                k,
                &ivl,
                p,
                None,
                BIG_DOUBLE,
                None,
                SolverOperation::Minimize,
            )
            .unwrap();
        if d.success && d.distance < 0.1 && best.is_none_or(|(_, bd)| d.distance < bd) {
            best = Some((e, d.distance));
        }
    }
    best.unwrap().0
}

/// The edge at `at` rounded `r` alone (a circular round), then the body stitched.
fn round(k: &mut Kernel, b: &mut Brep, at: [f64; 3], r: f64) -> Result<(), String> {
    let e = pick(k, b, at);
    {
        let mut s = FilletSession::new(k, b).unwrap();
        s.prepare_per_edge_constant_radius(
            k,
            &[e],
            &[r],
            &[FilletSurfaceGeneratorType::Circular],
            1,
            1.0,
        )
        .unwrap();
        if let Err(err) = s.do_fillet(k) {
            let info = s.error_info();
            return Err(format!(
                "{err:?}: {:?} {} (fit {:?})",
                info.code,
                info.comment,
                s.fit_size()
            ));
        }
    }
    stitch_and_orient(k, b, &StitchAndOrientOptions::default()).unwrap();
    Ok(())
}

fn volume(k: &Kernel, b: &mut Brep) -> f64 {
    let store = AttributeStore::default();
    let attrs = PropertyAttributes {
        store: &store,
        key_of: &|_| ObjectKey(0),
    };
    brep_volume(k, b, 1e-12, &attrs)
        .or_else(|_| brep_volume(k, b, 1e-9, &attrs))
        .unwrap()
}

fn clean(k: &mut Kernel, b: &mut Brep) {
    let lines = check_brep(k, b, CheckLevel::Level2).unwrap().lines();
    assert!(lines.is_empty(), "{lines:?}");
}

/// The plate's width, x of the arc wall's axis less 21.
const W: f64 = 41.0 + 1.0 / 6.0;
/// The arc wall's axis (x, y), vertical.
const AXIS: [f64; 2] = [W + 21.0, -21.0];
/// The arc wall's radius and the front round's.
const WALL: f64 = 21.0;
const FRONT_R: f64 = 2.0;
/// A point of the plate's top front edge, and of the middle edge (the fork's capi picks).
const FRONT: [f64; 3] = [20.583333333333332, -21.0, 10.0];
const MIDDLE: [f64; 3] = [41.17486688363082, -20.413218067968561, 9.415188112149071];

/// SolidCraft #851's plate with its concave R21 corner wall.
fn arc_wall_plate(k: &mut Kernel) -> Brep {
    let corner0 = create_box(k, Point3::new(W, -21.0, 0.0), 21.0, 21.0, 10.0).unwrap();
    let disc = create_cylinder(k, Point3::new(W + 21.0, -21.0, 0.0), WALL, 10.0).unwrap();
    let corner = boolean_difference(k, corner0, disc).unwrap();
    let plate = create_box(k, Point3::new(0.0, -21.0, 0.0), W, 21.0, 10.0).unwrap();
    boolean_union(k, plate, corner).unwrap()
}

/// The plate with its top front edge rounded R2.
fn front_rounded(k: &mut Kernel) -> Brep {
    let mut b = arc_wall_plate(k);
    round(k, &mut b, FRONT, FRONT_R).unwrap();
    b
}

/// The distance of `p` from the arc wall's axis.
fn rho(p: Point3) -> f64 {
    let (dx, dy) = (p.x - AXIS[0], p.y - AXIS[1]);
    (dx * dx + dy * dy).sqrt()
}

/// `p` less the front round's axis (y = -19, z = 8, along x), in (y, z).
fn off_front_axis(p: Point3) -> (f64, f64) {
    (p.y + 21.0 - FRONT_R, p.z - 10.0 + FRONT_R)
}

/// The point at angle `phi` about the arc wall's axis (from -x towards +y), `rho` from it,
/// at height `z`.
fn about_axis(phi: f64, rho: f64, z: f64) -> Point3 {
    Point3::new(AXIS[0] - rho * phi.cos(), AXIS[1] + rho * phi.sin(), z)
}

/// The closed forms of the middle edge rounded `r`, at its top End.
struct TopEnd {
    /// V's angle about the axis, asin(2 / 21): P is the half-plane at this angle.
    phi_v: f64,
    /// Where the rail on the R2 round meets its tangent edge T with the top face.
    q: Point3,
    /// The rolled rail's end on P (on the top face), and the wall rail's end on P.
    rail_end: Point3,
    wall_end: Point3,
    /// The ball centres at Q and on P.
    centre_q: Point3,
    centre_p: Point3,
    /// The vertex V.
    v: Point3,
}

impl TopEnd {
    fn new(r: f64) -> Self {
        let phi_v = (2.0 / WALL).asin();
        let phi_q = (2.0 / (WALL + r)).asin();
        Self {
            phi_v,
            q: about_axis(phi_q, WALL + r, 10.0),
            rail_end: about_axis(phi_v, WALL + r, 10.0),
            wall_end: about_axis(phi_v, WALL, 10.0 - r),
            centre_q: about_axis(phi_q, WALL + r, 10.0 - r),
            centre_p: about_axis(phi_v, WALL + r, 10.0 - r),
            v: about_axis(phi_v, WALL, 10.0),
        }
    }

    /// The signed distance of `p` from P (positive past it, away from the bottom).
    fn past_p(&self, p: Point3) -> f64 {
        let n = Vec3::new(self.phi_v.sin(), self.phi_v.cos(), 0.0);
        (p - self.v).dot(n)
    }
}

/// `n + 1` points along edge `e`'s curve.
fn samples(k: &Kernel, b: &Brep, e: EdgeId, n: usize) -> Vec<Point3> {
    let c = edge::curve(b, e).unwrap().unwrap();
    let ivl = edge::interval(b, e).unwrap();
    (0..=n)
        .map(|i| {
            let t = ivl.min() + (ivl.max() - ivl.min()) * (i as f64) / (n as f64);
            c.evaluate(k, t, 0, true, false).unwrap()[0]
        })
        .collect()
}

/// Edge `e`'s end points.
fn ends(b: &Brep, e: EdgeId) -> [Point3; 2] {
    let v0 = edge::vertex(b, e).unwrap().unwrap();
    let v1 = edge::other_vertex(b, e, Some(v0)).unwrap().unwrap();
    [vertex::point(b, v0).unwrap(), vertex::point(b, v1).unwrap()]
}

/// The edges every sample of which `on` takes within `tol`, with their worst value.
fn edges_where(k: &Kernel, b: &Brep, tol: f64, on: impl Fn(Point3) -> f64) -> Vec<(EdgeId, f64)> {
    let mut out = Vec::new();
    for e in get_edges(b).unwrap() {
        if edge::curve(b, e).unwrap().is_none() {
            continue;
        }
        let worst = samples(k, b, e, 16)
            .into_iter()
            .map(|p| on(p).abs())
            .fold(0.0, f64::max);
        if worst <= tol {
            out.push((e, worst));
        }
    }
    out
}

/// Face `f`'s unit normal (the surface's, unsigned) at its point nearest `p`.
fn normal_at(k: &Kernel, b: &Brep, f: FaceId, p: Point3) -> Vec3 {
    let s = face2::surface(b, f).unwrap().unwrap();
    let dom = s.natural_uv_domain_with_nurb(k).unwrap();
    let d = s
        .drop_point::<NoObject>(k, p, &dom, &DropPointOptions::default())
        .unwrap();
    assert!(d.success, "no foot of {p:?}");
    s.evaluate_normal(k, d.result_uv, Side::Right, Side::Right)
        .unwrap()
}

/// The largest angle (radians) between the normals of edge `e`'s two faces along it.
fn crease(k: &Kernel, b: &Brep, e: EdgeId) -> f64 {
    let fs = edge::faces(b, e, None).unwrap();
    assert_eq!(fs.len(), 2, "edge {e:?} is not between two faces");
    let mut worst = 0.0f64;
    for p in samples(k, b, e, 16) {
        let n0 = normal_at(k, b, fs[0], p);
        let n1 = normal_at(k, b, fs[1], p);
        worst = worst.max(n0.dot(n1).abs().min(1.0).acos());
    }
    worst
}

/// The one edge whose samples all satisfy `on` within `tol`, and its worst value.
fn the_edge(
    k: &Kernel,
    b: &Brep,
    what: &str,
    tol: f64,
    on: impl Fn(Point3) -> f64,
) -> (EdgeId, f64) {
    let found = edges_where(k, b, tol, on);
    assert_eq!(found.len(), 1, "{what}: {found:?}");
    found[0]
}

/// Some end of edge `e` at `want`.
fn assert_ends_at(b: &Brep, e: EdgeId, what: &str, want: Point3, tol: f64) {
    let [a, z] = ends(b, e);
    let d = a.distance(want).min(z.distance(want));
    assert!(
        d <= tol,
        "{what}: ends {a:?}, {z:?}, want {want:?} ({d:e} off)"
    );
}

/// Curves with a closed form the kernel computes exactly (the rolled rail on the top face,
/// the cap's edges, the rollover seam and the End's vertices): off it by rounding only.
const EXACT_TOL: f64 = 1e-9;
/// The end rib is the traced round's section by P, off the ball's circle by up to 3.7e-6.
const SECTION_TOL: f64 = 1e-5;
/// The rail on the R2 round is traced again from the bottom End to Q at the solver's
/// approximation tolerance, 1.6e-4 on this part (`approx_tol_3d`, as the fork's
/// `ReCalcFilletGeom`): it lies off the round by up to 1.5e-4, and its ball off the wall
/// by up to 1.2e-5 (r 0.25).
const TRACE_TOL: f64 = 1.6e-4;
/// Two faces meet G1 where their normals agree within this angle (radians): the seam's
/// faces are traced, 8.2e-5 apart at most.
const G1_TOL: f64 = 1e-4;
/// The round and the R2 round meet along the rail traced again within this angle: the
/// re-traced rail is off the round by up to 1.5e-4 all along it, and the faces' normals
/// there by up to 1.7e-3 (r 0.25). The trace's own angular tolerance is the solver's 20
/// degrees.
const TRACE_G1_TOL: f64 = 2e-3;

/// Checks the middle edge rounded `r` at its top End: the rail rolls over T onto the top
/// face, the round from there to P is the top face's and the wall's, and the cap's edges
/// lie on the faces' sections by P.
fn check_top_end(k: &Kernel, b: &Brep, r: f64) {
    let t = TopEnd::new(r);

    // the rolled rail: on the top face (z = 10) on the circle of radius 21 + r about the
    // axis, from T (at Q) to P, the top face G1 to the round on it
    let (rail, off) = the_edge(k, b, "the rolled rail", EXACT_TOL, |p| {
        if t.past_p(p) > EXACT_TOL || p.distance(t.q) > 4.0 * r {
            1.0
        } else {
            (p.z - 10.0).abs().max((rho(p) - (WALL + r)).abs())
        }
    });
    eprintln!("r {r}: the rolled rail off its closed form by {off:e}");
    assert_ends_at(b, rail, "the rolled rail at T", t.q, EXACT_TOL);
    assert_ends_at(b, rail, "the rolled rail on P", t.rail_end, EXACT_TOL);
    let c = crease(k, b, rail);
    eprintln!("r {r}: the rolled rail's crease {c:e}");
    assert!(c <= EXACT_TOL, "the rolled rail is not G1: {c:e}");

    // the cap's edge on the top face: the top face's line in P from the rail's end to V
    let seg = t.v - t.rail_end;
    let (cap_top, off) = the_edge(k, b, "the cap's edge on the top face", EXACT_TOL, |p| {
        let s = ((p - t.rail_end).dot(seg) / seg.dot(seg)).clamp(0.0, 1.0);
        p.distance(t.rail_end + seg * s)
    });
    eprintln!("r {r}: the cap's edge on the top face off by {off:e}");
    assert_ends_at(
        b,
        cap_top,
        "the top cap edge at the rail",
        t.rail_end,
        EXACT_TOL,
    );
    assert_ends_at(b, cap_top, "the top cap edge at V", t.v, EXACT_TOL);

    // the cap's edge on the wall: the wall's line in P (vertical) from the wall rail to V
    let (cap_wall, _) = the_edge(k, b, "the cap's edge on the wall", EXACT_TOL, |p| {
        p.distance(Point3::new(t.v.x, t.v.y, p.z.clamp(10.0 - r, 10.0)))
    });
    assert_ends_at(
        b,
        cap_wall,
        "the wall cap edge at its rail",
        t.wall_end,
        EXACT_TOL,
    );

    // the end rib: the ball's circle on P from the wall's rail to the top face's
    let (rib, off) = the_edge(k, b, "the end rib", SECTION_TOL, |p| {
        (p.distance(t.centre_p) - r).abs().max(t.past_p(p).abs())
    });
    eprintln!("r {r}: the end rib off its circle by {off:e}");
    assert_ends_at(b, rib, "the end rib on the top face", t.rail_end, EXACT_TOL);
    assert_ends_at(b, rib, "the end rib on the wall", t.wall_end, EXACT_TOL);

    // the seam where the rail rolls over: the ball's circle in the section plane at Q (the
    // half-plane through the axis and Q), between the round on the R2 round and the round
    // on the top face, G1
    let nq = {
        let (dx, dy) = (t.centre_q.x - AXIS[0], t.centre_q.y - AXIS[1]);
        Vec3::new(dy, -dx, 0.0) * (1.0 / (dx * dx + dy * dy).sqrt())
    };
    let (seam, off) = the_edge(k, b, "the rollover seam", EXACT_TOL, |p| {
        (p.distance(t.centre_q) - r)
            .abs()
            .max((p - t.centre_q).dot(nq).abs())
    });
    eprintln!("r {r}: the rollover seam off its circle by {off:e}");
    assert_ends_at(b, seam, "the seam at T", t.q, EXACT_TOL);
    let c = crease(k, b, seam);
    eprintln!("r {r}: the rollover seam's crease {c:e}");
    assert!(c <= G1_TOL, "the rollover seam is not G1: {c:e}");

    // the rail on the R2 round runs up to T at Q, the ball on it touching the wall: each
    // point is on the round, and the ball's centre there (r in from the round) is 21 + r
    // from the wall's axis
    let on_round = |p: Point3| {
        let (dy, dz) = off_front_axis(p);
        ((dy * dy + dz * dz).sqrt() - FRONT_R).abs()
    };
    let mut round_rail = None;
    for e in get_edges(b).unwrap() {
        if edge::curve(b, e).unwrap().is_none() || e == seam {
            continue;
        }
        let [a, z] = ends(b, e);
        if a.distance(t.q).min(z.distance(t.q)) > EXACT_TOL {
            continue;
        }
        let pts = samples(k, b, e, 16);
        if pts.iter().all(|&p| on_round(p) <= TRACE_TOL) && pts.iter().any(|p| p.z < 10.0 - r) {
            round_rail = Some(e);
        }
    }
    let round_rail = round_rail.expect("the rail on the R2 round, ending at Q");
    let mut worst = (0.0_f64, t.q);
    for p in samples(k, b, round_rail, 32) {
        let (dy, dz) = off_front_axis(p);
        let l = (dy * dy + dz * dz).sqrt();
        let centre = Point3::new(p.x, p.y - r * dy / l, p.z - r * dz / l);
        let miss = (rho(centre) - (WALL + r)).abs();
        if miss > worst.0 {
            worst = (miss, p);
        }
    }
    let (miss, p) = worst;
    eprintln!("r {r}: the R2 round's rail's ball misses the wall by {miss:e} at {p:?}");
    assert!(
        miss <= TRACE_TOL,
        "the R2 round's rail at {p:?}: its ball misses the wall by {miss:e}"
    );
    let c = crease(k, b, round_rail);
    eprintln!("r {r}: the R2 round's rail's crease {c:e}");
    assert!(c <= TRACE_G1_TOL, "the R2 round's rail is not G1: {c:e}");
}

/// The middle edge rounded 0.25 builds: its rail rolls over T onto the top face before the
/// top End, and the cap's edge on the top face is the top face's line in P. The fork
/// (18e96f1) refuses it, `SMLIB_ERR_NO_ROOM` with no size.
#[test]
fn the_middle_edge_r0_25_rolls_over_onto_the_top_face_before_its_end() {
    let mut k = all_fixes(true);
    // (built twice, not cloned: a cloned body's surfaces have no owning faces)
    let mut b = front_rounded(&mut k);
    let mut before = front_rounded(&mut k);
    round(&mut k, &mut b, MIDDLE, 0.25).unwrap();
    clean(&mut k, &mut b);
    let removed = volume(&k, &mut before) - volume(&k, &mut b);
    eprintln!(
        "r 0.25: removed {removed}, faces {}",
        get_faces(&b).unwrap().len()
    );
    // the round takes less than its ball's corner, r^2 (1 - pi/4), along the edge (a
    // quarter of the R2 round's circle and the rolled stretch)
    let corner = 0.0625 * (1.0 - std::f64::consts::FRAC_PI_4);
    assert!(
        removed > 0.0 && removed < corner * (std::f64::consts::PI + 0.1),
        "{removed}"
    );
    check_top_end(&k, &b, 0.25);
}

/// At r 1 the rolled stretch is four times as long (0.096): it builds the same way.
#[test]
fn the_middle_edge_r1_rolls_over_onto_the_top_face_before_its_end() {
    let mut k = all_fixes(true);
    let mut b = front_rounded(&mut k);
    round(&mut k, &mut b, MIDDLE, 1.0).unwrap();
    clean(&mut k, &mut b);
    check_top_end(&k, &b, 1.0);
}

/// Without `fillet_blend_corner` the middle edge's round is the pinned one's: refused.
#[test]
fn without_the_blend_corner_the_middle_edge_is_refused() {
    let mut k = all_fixes(false);
    let mut b = front_rounded(&mut k);
    assert!(round(&mut k, &mut b, MIDDLE, 0.25).is_err());
}
