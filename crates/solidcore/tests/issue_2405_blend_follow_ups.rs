//! #2405 (from tangency-labs/solidcraft#10850 and tangency-labs/smlib#97): follow-ups of
//! the generic blend's corner work, each checked against a closed form.
//!
//! - The exact round classes' cylinder (`CalcLineSpineFilletGeom`) is built on its exact
//!   spine, the line the offsets meet in, and keeps that frame. The fork's 1.06e-12 on the
//!   oblique Ball came from the solved stations it took the centres from (the far one
//!   2e-12 off the spine here), and the analytic recognition of the patch
//!   (`CopyAndAddAnalytics` → `IsNurbSurfaceCone`, the end circles' `IsArc`) then refitted
//!   the frame from the net (7e-15 in the axis's direction, 2e-13 at C).
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers
use solidcore::blend::fillet_executive::session::FilletSession;
use solidcore::blend::fillet_solver::p1::FilletSurfaceGeneratorType;
use solidcore::boolean::api_brep::boolean_difference;
use solidcore::curve::curve::Curve;
use solidcore::geom::global_solver::SolverOperation;
use solidcore::kernel::UpstreamFixes;
use solidcore::math::{Point3, Vec3, BIG_DOUBLE};
use solidcore::primitive::api_curves::create_line_segment;
use solidcore::primitive::api_primitives::create_box;
use solidcore::query::api_queries::{get_edges, get_faces};
use solidcore::surface::surface::p2::SurfaceKind;
use solidcore::sweep::api_brep::create_linear_sweep;
use solidcore::topology::brep::p5::{stitch_and_orient, StitchAndOrientOptions};
use solidcore::topology::face::p2 as face2;
use solidcore::topology::{edge, EdgeId};
use solidcore::{Brep, Kernel};

/// Every fix but `geometric_g2_continuity`, with or without `fillet_exact_round_classes`.
fn all_fixes(exact_classes: bool) -> Kernel {
    let mut k = Kernel::default();
    let mut f = UpstreamFixes::default();
    for &name in UpstreamFixes::NAMES {
        if name != "geometric_g2_continuity"
            && (exact_classes || name != "fillet_exact_round_classes")
        {
            *f.flag_mut(name).unwrap() = true;
        }
    }
    k.fixes = f;
    k
}

fn unit(v: Vec3) -> Vec3 {
    v / v.length()
}

/// The edge of `b` nearest `at` (within 0.1).
fn pick(k: &Kernel, b: &Brep, at: Vec3) -> EdgeId {
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
                at,
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

/// Rounds the edges at `at` at radius `r` together.
fn round(k: &mut Kernel, b: &mut Brep, at: &[Vec3], r: f64) {
    let edges: Vec<EdgeId> = at.iter().map(|&p| pick(k, b, p)).collect();
    let xs = vec![FilletSurfaceGeneratorType::Circular; edges.len()];
    let radii = vec![r; edges.len()];
    {
        let mut s = FilletSession::new(k, b).unwrap();
        s.prepare_per_edge_constant_radius(k, &edges, &radii, &xs, 1, 1.0)
            .unwrap();
        s.do_fillet(k).unwrap();
    }
    stitch_and_orient(k, b, &StitchAndOrientOptions::default()).unwrap();
}

/// The prism swept `height` along `dir` from the closed polygon `pts`.
fn prism(k: &mut Kernel, pts: &[Point3], dir: Vec3, height: f64) -> Brep {
    let sides: Vec<Curve> = (0..pts.len())
        .map(|i| Curve::from(create_line_segment(pts[i], pts[(i + 1) % pts.len()]).unwrap()))
        .collect();
    let refs: Vec<Option<&Curve>> = sides.iter().map(Some).collect();
    create_linear_sweep(k, &refs, dir, height, true).unwrap()
}

/// The box 400 wide on the far side (along `n`) of the plane through `at` with normal `n`.
fn half_space(k: &mut Kernel, at: Vec3, n: Vec3) -> Brep {
    let u = unit(n.cross(if n.x.abs() < 0.9 { Vec3::X } else { Vec3::Y }));
    let v = n.cross(u);
    let h = 200.0;
    let c = [
        at + (u + v) * h,
        at + (v - u) * h,
        at - (u + v) * h,
        at + (u - v) * h,
    ];
    prism(k, &c, n, 400.0)
}

/// The fork's capi `ball oblique_corner` (as `wr697_ball_corner`): three planes through
/// the origin, none square to another, with outward normals `n`, their edges `e` 20 long
/// and cut square at their far ends; the three edges at the origin rounded r 2. Returns
/// the rounded body, `n`, `e` and the corner ball's centre C (`n_i . C = -r`).
fn oblique_corner(k: &mut Kernel, rad: f64) -> (Brep, [Vec3; 3], [Vec3; 3], Vec3) {
    let len = 20.0;
    let n = [
        unit(Vec3::new(1.0, 0.2, 0.0)),
        unit(Vec3::new(0.0, 1.0, 0.2)),
        unit(Vec3::new(0.2, 0.0, 1.0)),
    ];
    let mut e = [Vec3::ZERO; 3];
    for i in 0..3 {
        e[i] = unit(n[i].cross(n[(i + 1) % 3]));
        if e[i].dot(n[(i + 2) % 3]) > 0.0 {
            e[i] = -e[i];
        }
    }
    let det = n[0].dot(n[1].cross(n[2]));
    let c = (n[1].cross(n[2]) + n[2].cross(n[0]) + n[0].cross(n[1])) * (-rad / det);
    let mut body = create_box(k, Point3::new(-40.0, -40.0, -40.0), 80.0, 80.0, 80.0).unwrap();
    for i in 0..3 {
        let hs = half_space(k, Vec3::ZERO, n[i]);
        body = boolean_difference(k, body, hs).unwrap();
        let hs = half_space(k, e[i] * len, e[i]);
        body = boolean_difference(k, body, hs).unwrap();
    }
    let points: Vec<Vec3> = e.iter().map(|&ei| ei * (len / 2.0)).collect();
    round(k, &mut body, &points, rad);
    (body, n, e, c)
}

/// The three rounds' cones of the oblique corner, each with the edge it rounds: the cone
/// whose axis is nearest that edge's direction.
fn round_cones(b: &Brep, e: &[Vec3; 3]) -> Vec<(usize, solidcore::surface::cone::Cone)> {
    let mut out = Vec::new();
    for f in get_faces(b).unwrap() {
        let s = face2::surface(b, f).unwrap().unwrap();
        if let SurfaceKind::Cone(cone) = s.kind() {
            let z = cone.rev().position().z_axis();
            let i = (0..3)
                .min_by(|&a, &bb| {
                    z.cross(e[a])
                        .length()
                        .partial_cmp(&z.cross(e[bb]).length())
                        .unwrap()
                })
                .unwrap();
            out.push((i, (**cone).clone()));
        }
    }
    out
}

/// Item 1: each round of the oblique corner is the exact cylinder of its spine. The ball
/// of radius r touching planes i and i + 1 rolls along their edge with its centre on the
/// offsets' line through C along `e_i = n_i x n_(i+1) / |n_i x n_(i+1)|`, so the round's
/// axis is that line and its radius r.
#[test]
fn the_exact_rounds_cylinder_is_on_its_exact_spine() {
    let mut k = all_fixes(true);
    let rad = 2.0;
    let (b, _n, e, c) = oblique_corner(&mut k, rad);
    let cones = round_cones(&b, &e);
    assert_eq!(cones.len(), 3, "three rounds, each a cone");
    let mut worst = (0.0f64, 0.0f64, 0.0f64);
    for (i, cone) in &cones {
        let frame = cone.rev().position();
        let z = frame.z_axis();
        // the axis's distance from C, its direction's from e_i, the radius's from r
        let off = (c - frame.origin()).cross(z).length();
        let dir = z.cross(e[*i]).length();
        let r = (cone.bot_radius() - rad).abs();
        eprintln!(
            "round {i}: axis {off:.3e} from C, {dir:.3e} from e; radius {r:.3e} off; semi-angle {}",
            cone.semi_angle_deg()
        );
        assert_eq!(cone.semi_angle_deg(), 0.0);
        worst = (worst.0.max(off), worst.1.max(dir), worst.2.max(r));
    }
    // to rounding: the direction to 1e-15, the axis's distance from C and the radius to
    // 1e-15 of the part's size (the frame's origin is the spine's far end, 28 from C,
    // whose coordinates round at 3.6e-15)
    let size = 48.66;
    assert!(
        worst.0 <= 1e-15 * size && worst.1 <= 1e-15 && worst.2 <= 1e-15 * size,
        "{worst:?}"
    );
}
