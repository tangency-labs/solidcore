//! #2933 (WR-445): a corner setback is measured from the set-back vertex on every edge, with
//! the default kernel. A box 60 x 40 x 20 with r5 rounds on the three edges at (60, 40, 20)
//! and that corner set back 8, as the goldens' `real-box-high-corner-setback-eight` builds it.
//! Each round ends on its section 8 from the vertex, so the six section ends lie on the
//! setback planes x = 52, y = 32 and z = 12, and the body is convex: its volume lies in the
//! bounds the kernel case `ladder/07-fillets/013-corner-setback` derives in closed form.
//! The corner's patch meets each round tangentially, and every edge lies within the
//! kernel's 4e-5 of its faces.
#![allow(clippy::unwrap_used, clippy::indexing_slicing)] // test helpers

use solidcore::blend::fillet_executive::session::FilletSession;
use solidcore::mass::integrator::properties::brep_volume_default;
use solidcore::math::Point3;
use solidcore::nurbs::knots::Side;
use solidcore::primitive::api_primitives::create_box;
use solidcore::query::api_queries::{
    brep_bounding_box, brep_is_manifold_solid, surface_closest_point,
};
use solidcore::surface::surface::{Surface, SurfaceKind};
use solidcore::topology::brep::p5::{stitch_and_orient, StitchAndOrientOptions};
use solidcore::topology::face::p2 as face;
use solidcore::topology::{edge, vertex, Brep};
use solidcore::Kernel;

const CORNER: [f64; 3] = [60.0, 40.0, 20.0];

/// The box with its far corner set back 8, built with the default kernel.
fn build() -> (Kernel, Brep) {
    let mut k = Kernel::default();
    let mut b = create_box(&mut k, Point3::new(0.0, 0.0, 0.0), 60.0, 40.0, 20.0).unwrap();
    let corner = Point3::new(CORNER[0], CORNER[1], CORNER[2]);
    let v = b
        .vertices()
        .ids()
        .find(|&v| vertex::point(&b, v).unwrap() == corner)
        .unwrap();
    let edges: Vec<_> = b
        .edges()
        .ids()
        .filter(|&e| {
            let (s, t) = edge::start_end_vertices(&b, e).unwrap();
            s == Some(v) || t == Some(v)
        })
        .collect();
    assert_eq!(edges.len(), 3);
    {
        let mut s = FilletSession::new(&mut k, &mut b).unwrap();
        s.prepare_constant_radius_edges(&mut k, &edges, 5.0, 1, 1, 1.0)
            .unwrap();
        // set back (`AddSetBackVertex`), so the corner is the N x N corner, not the Ball
        s.add_set_back_vertex(v);
        s.create_corners(&k, None).unwrap();
        let corners: Vec<_> = s
            .fillet_corners()
            .unwrap()
            .into_iter()
            .filter(|&id| s.corner(id).unwrap().filleted_vertex() == Some(v))
            .collect();
        assert_eq!(corners.len(), 1);
        s.set_corner_setback(corners[0], 8.0).unwrap();
        s.do_filleting(&mut k, None).unwrap();
    }
    stitch_and_orient(&mut k, &mut b, &StitchAndOrientOptions::default()).unwrap();
    (k, b)
}

#[test]
fn each_round_ends_eight_from_the_vertex() {
    let (_, b) = build();
    // the rounds' sections at the setback planes: x = 52 (the x edge), y = 32, z = 12
    let ends = [
        [52.0, 35.0, 20.0],
        [52.0, 40.0, 15.0],
        [55.0, 32.0, 20.0],
        [60.0, 32.0, 15.0],
        [55.0, 40.0, 12.0],
        [60.0, 35.0, 12.0],
    ];
    let points: Vec<Point3> = b
        .vertices()
        .ids()
        .map(|v| vertex::point(&b, v).unwrap())
        .collect();
    for end in ends {
        let want = Point3::new(end[0], end[1], end[2]);
        let nearest = points
            .iter()
            .map(|p| p.distance(want))
            .fold(f64::INFINITY, f64::min);
        assert!(nearest < 1e-6, "no vertex at {end:?}: nearest {nearest}");
    }
}

#[test]
fn the_set_back_box_is_convex_and_inside_its_bounds() {
    let (k, mut b) = build();
    assert!(brep_is_manifold_solid(&b).unwrap());
    let x = brep_bounding_box(&k, &b, true).unwrap();
    for (i, limit) in CORNER.iter().enumerate() {
        assert!(x.max_at(i) <= limit + 1e-6, "{i}: {}", x.max_at(i));
    }
    // 013-corner-setback's closed-form bounds: the hull of the sections in the corner cube
    // below, the three rounds' prisms above
    let volume = brep_volume_default(&k, &mut b, 1e-6).unwrap();
    assert!(
        (47_271.975_7..=47_384.893_5).contains(&volume),
        "volume {volume}"
    );
}

/// The largest angle (degrees) between the normals of an edge's two faces, sampled along
/// it, over the edges where a face is not a plane and the faces meet within 45 degrees (the
/// rounds' and the patch's seams).
fn max_seam_angle(k: &Kernel, b: &Brep) -> f64 {
    let mut worst: f64 = 0.0;
    for e in b.edges().ids() {
        let faces = edge::faces(b, e, None).unwrap();
        let [f1, f2] = faces[..] else { continue };
        let (s1, s2) = (
            face::surface(b, f1).unwrap().unwrap(),
            face::surface(b, f2).unwrap().unwrap(),
        );
        if matches!(s1.kind(), SurfaceKind::Plane(_)) && matches!(s2.kind(), SurfaceKind::Plane(_))
        {
            continue;
        }
        let c = edge::curve(b, e).unwrap().unwrap();
        let iv = edge::interval(b, e).unwrap();
        let mut along: f64 = 0.0;
        for i in 0..=20 {
            let p = c
                .evaluate_point(k, iv.min() + iv.length() * f64::from(i) / 20.0)
                .unwrap();
            let normal = |s: &Surface| {
                let uv = surface_closest_point(k, s, p).unwrap().uvs[0];
                s.evaluate_normal(k, uv, Side::Left, Side::Left).unwrap()
            };
            let (n1, n2) = (normal(s1), normal(s2));
            along = along.max(n1.cross(n2).length().atan2(n1.dot(n2).abs()).to_degrees());
        }
        // a round ending square on a box face far from the corner is no seam
        if along < 45.0 {
            worst = worst.max(along);
        }
    }
    worst
}

#[test]
fn the_patch_meets_the_rounds_tangentially_within_tolerance() {
    let (k, mut b) = build();
    b.refresh_cached_topology();
    assert!(solidcore::topology::integrity::check(&b).is_ok());
    let seam = max_seam_angle(&k, &b);
    assert!(seam < 0.01, "seam {seam} degrees");
    let tol = b
        .edges()
        .ids()
        .map(|e| edge::tolerance(&b, e).unwrap())
        .fold(0.0, f64::max);
    assert!(tol <= 4.0e-5 * (1.0 + 1e-9), "edge tolerance {tol}");
}
