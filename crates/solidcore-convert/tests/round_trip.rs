//! Every record kind goes to the kernel as the kernel kind of the same fields (checked one way
//! against literals, so a swap made in both directions is caught), and back unchanged, bit for
//! bit; a malformed record is carried over for the kernel to refuse.

use solidcore::eval::{self as k, CurveEval, SurfaceEval};
use solidcore::prism::profile as kp;
use solidcore_convert::{
    curve, curve2, frame, loop2, plane, record_curve, record_curve2, record_frame, record_list,
    record_loop2, record_pair, record_plane, record_region2, record_seg2, record_surface,
    record_table, record_triple, record_vec2, record_vec3, region2, seg2, surface, vec2, vec3,
};
use solidcore_record::brep::{Curve, Curve2, Frame, Surface};
use solidcore_record::{Loop2, Plane, Region2, Seg2, Vec2, Vec3, MAX_BEZIER_CTRL};

const O: [f64; 3] = [1.5, -2.25, 3.125];
const X: [f64; 3] = [0.6, 0.8, 0.0];
const Y: [f64; 3] = [0.0, 0.0, 1.0];

fn f() -> Frame {
    Frame {
        origin: O,
        x_axis: X,
        y_axis: Y,
    }
}

fn kf() -> k::Frame {
    k::Frame {
        origin: O,
        x_axis: X,
        y_axis: Y,
    }
}

const M: [[f64; 4]; 4] = [
    [1.0, 2.0, 3.0, 4.0],
    [5.0, 6.0, 7.0, 8.0],
    [9.0, 10.0, 11.0, 12.0],
    [13.0, 14.0, 15.0, 16.0],
];

/// Each record curve with the kernel curve it must become.
fn curves() -> Vec<(Curve, k::Curve)> {
    let pts = vec![[0.0, 0.5, 1.0], [1.0, 1.5, 2.0], [-1.0, 0.25, 7.0]];
    let hom = vec![[0.0, 0.5, 1.0, 1.0], [2.0, 3.0, 4.0, 2.0]];
    vec![
        (
            Curve::Line {
                start: [0.1, 0.2, 0.3],
                end: [1.0, 2.0, 3.0],
            },
            k::Curve::Line {
                start: [0.1, 0.2, 0.3],
                end: [1.0, 2.0, 3.0],
            },
        ),
        (
            Curve::Circle {
                frame: f(),
                radius: 2.5,
            },
            k::Curve::Circle {
                frame: kf(),
                radius: 2.5,
            },
        ),
        (
            Curve::Ellipse {
                frame: f(),
                radii: [3.0, 1.25],
            },
            k::Curve::Ellipse {
                frame: kf(),
                radii: [3.0, 1.25],
            },
        ),
        (
            Curve::Polyline {
                points: pts.clone(),
            },
            k::Curve::Polyline {
                points: pts.clone(),
            },
        ),
        (
            Curve::Spline {
                knots: vec![0.0, 0.0, 0.5, 1.0, 1.0],
                points: pts.clone(),
            },
            k::Curve::Spline {
                knots: vec![0.0, 0.0, 0.5, 1.0, 1.0],
                points: pts,
            },
        ),
        (
            Curve::RationalSpline {
                knots: vec![0.0, 0.0, 1.0, 1.0],
                homogeneous: hom.clone(),
            },
            k::Curve::RationalSpline {
                knots: vec![0.0, 0.0, 1.0, 1.0],
                homogeneous: hom,
            },
        ),
    ]
}

/// Each record surface with the kernel surface it must become: every u/v pair distinct.
fn surfaces() -> Vec<(Surface, k::Surface)> {
    let (uk, vk) = (vec![0.0, 0.0, 1.0, 1.0], vec![-1.0, -1.0, 0.0, 2.0, 2.0]);
    // Two rows along u of three points along v: a transposition changes the shape.
    let pts = vec![
        vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.1], [0.0, 2.0, 0.3]],
        vec![[1.0, 0.0, 0.5], [1.0, 1.0, 0.7], [1.0, 2.0, 1.1]],
    ];
    let hom = vec![
        vec![
            [0.0, 0.0, 0.0, 1.0],
            [0.0, 0.5, 0.25, 0.5],
            [0.0, 2.0, 0.3, 1.0],
        ],
        vec![
            [2.0, 0.0, 1.0, 2.0],
            [1.0, 1.0, 0.7, 1.0],
            [0.75, 1.5, 0.5, 0.75],
        ],
    ];
    let line = (
        Curve::Line {
            start: [2.0, 0.0, 0.0],
            end: [3.0, 0.0, 1.0],
        },
        k::Curve::Line {
            start: [2.0, 0.0, 0.0],
            end: [3.0, 0.0, 1.0],
        },
    );
    vec![
        (
            Surface::Plane {
                origin: [0.0, 1.0, 2.0],
                u_axis: [1.0, 0.0, 0.0],
                v_axis: [0.0, 1.0, 0.0],
            },
            k::Surface::Plane {
                origin: [0.0, 1.0, 2.0],
                u_axis: [1.0, 0.0, 0.0],
                v_axis: [0.0, 1.0, 0.0],
            },
        ),
        (
            Surface::PlaneThrough {
                origin: [0.0; 3],
                u_point: [4.0, 0.0, 0.0],
                v_point: [0.0, 5.0, 0.0],
            },
            k::Surface::PlaneThrough {
                origin: [0.0; 3],
                u_point: [4.0, 0.0, 0.0],
                v_point: [0.0, 5.0, 0.0],
            },
        ),
        (
            Surface::Cylinder {
                frame: f(),
                radius: 2.0,
            },
            k::Surface::Cylinder {
                frame: kf(),
                radius: 2.0,
            },
        ),
        (
            Surface::Cone {
                frame: f(),
                radius_at_origin: 2.0,
                half_angle: 0.3,
            },
            k::Surface::Cone {
                frame: kf(),
                radius_at_origin: 2.0,
                half_angle: 0.3,
            },
        ),
        (
            Surface::Sphere {
                frame: f(),
                radius: 2.0,
            },
            k::Surface::Sphere {
                frame: kf(),
                radius: 2.0,
            },
        ),
        (
            Surface::Torus {
                frame: f(),
                major_radius: 5.0,
                minor_radius: 1.0,
                tube_start: 0.5,
            },
            k::Surface::Torus {
                frame: kf(),
                major_radius: 5.0,
                minor_radius: 1.0,
                tube_start: 0.5,
            },
        ),
        (
            Surface::BallPatch {
                center: [1.0, 2.0, 3.0],
                radius: 2.0,
            },
            k::Surface::BallPatch {
                center: [1.0, 2.0, 3.0],
                radius: 2.0,
            },
        ),
        (
            Surface::Spline {
                u_knots: uk.clone(),
                v_knots: vk.clone(),
                points: pts.clone(),
            },
            k::Surface::Spline {
                u_knots: uk.clone(),
                v_knots: vk.clone(),
                points: pts,
            },
        ),
        (
            Surface::RationalSpline {
                u_knots: uk.clone(),
                v_knots: vk.clone(),
                homogeneous: hom.clone(),
            },
            k::Surface::RationalSpline {
                u_knots: uk,
                v_knots: vk,
                homogeneous: hom,
            },
        ),
        (
            Surface::PlacedRevolution {
                profile: Box::new(line.0.clone()),
                origin: [0.5, 0.0, 0.0],
                axis: [0.0, 0.0, 1.0],
                matrix: M,
                forward: false,
            },
            k::Surface::PlacedRevolution {
                profile: Box::new(line.1.clone()),
                origin: [0.5, 0.0, 0.0],
                axis: [0.0, 0.0, 1.0],
                matrix: M,
                forward: false,
            },
        ),
        (
            Surface::PlacedRevolution {
                profile: Box::new(line.0),
                origin: [0.0; 3],
                axis: [1.0, 0.0, 0.0],
                matrix: M,
                forward: true,
            },
            k::Surface::PlacedRevolution {
                profile: Box::new(line.1),
                origin: [0.0; 3],
                axis: [1.0, 0.0, 0.0],
                matrix: M,
                forward: true,
            },
        ),
    ]
}

#[test]
fn every_curve_kind_becomes_the_kernel_kind_and_comes_back() {
    let mut all = curves();
    let s = surfaces();
    let (sa, ka) = s[0].clone();
    let (sb, kb) = s[4].clone();
    let (la, lk) = all[0].clone();
    all.push((
        Curve::SurfaceIntersection {
            surfaces: [Box::new(sa), Box::new(sb)],
            leader: Box::new(la),
            tolerance: 1e-6,
        },
        k::Curve::SurfaceIntersection {
            surfaces: [Box::new(ka), Box::new(kb)],
            leader: Box::new(lk),
            tolerance: 1e-6,
        },
    ));
    for (r, want) in all {
        assert_eq!(curve(&r), want);
        assert_eq!(record_curve(curve(&r)), r);
    }
}

#[test]
fn every_surface_kind_becomes_the_kernel_kind_and_comes_back() {
    for (r, want) in surfaces() {
        assert_eq!(surface(&r), want);
        assert_eq!(record_surface(surface(&r)), r);
    }
}

#[test]
fn pcurves_frames_and_planes_become_the_kernel_kinds_and_come_back() {
    let (kn, pts, hom) = (
        vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        vec![[0.0, 0.5], [1.0, 2.0], [3.0, -1.0]],
        vec![[1.0, 0.0, 1.0], [0.5, 0.25, 0.5], [0.0, 1.0, 1.0]],
    );
    for (r, want) in [
        (
            Curve2::Spline {
                knots: kn.clone(),
                points: pts.clone(),
            },
            k::Curve2::Spline {
                knots: kn.clone(),
                points: pts,
            },
        ),
        (
            Curve2::RationalSpline {
                knots: kn.clone(),
                homogeneous: hom.clone(),
            },
            k::Curve2::RationalSpline {
                knots: kn,
                homogeneous: hom,
            },
        ),
    ] {
        assert_eq!(curve2(&r), want);
        assert_eq!(record_curve2(curve2(&r)), r);
    }
    assert_eq!(frame(&f()), kf());
    assert_eq!(record_frame(&frame(&f())), f());
    let p = Plane {
        origin: Vec3::new(1.0, 2.0, 3.0),
        x: Vec3::new(0.0, 1.0, 0.0),
        y: Vec3::new(0.0, 0.0, 1.0),
    };
    assert_eq!(
        plane(&p),
        k::Frame {
            origin: [1.0, 2.0, 3.0],
            x_axis: [0.0, 1.0, 0.0],
            y_axis: [0.0, 0.0, 1.0]
        }
    );
    assert_eq!(record_plane(&plane(&p)), p);
}

/// Copies keep every bit: a negative zero, a NaN payload, an infinity, a subnormal.
#[test]
fn vectors_keep_every_bit() {
    let odd = [
        -0.0,
        f64::from_bits(0x7ff8_0000_dead_beef),
        f64::NEG_INFINITY,
        f64::from_bits(1),
    ];
    for x in odd {
        for v in [
            Vec3::new(x, 1.0, -2.0),
            Vec3::new(3.0, x, 0.25),
            Vec3::new(0.5, -7.0, x),
        ] {
            let kv = vec3(v);
            assert_eq!(
                [kv.x.to_bits(), kv.y.to_bits(), kv.z.to_bits()],
                [v.x.to_bits(), v.y.to_bits(), v.z.to_bits()]
            );
            let back = record_vec3(kv);
            assert_eq!(
                [back.x.to_bits(), back.y.to_bits(), back.z.to_bits()],
                [v.x.to_bits(), v.y.to_bits(), v.z.to_bits()]
            );
        }
    }
    let (a, b, c) = (
        solidcore::Vec3::new(1.0, 2.0, 3.0),
        solidcore::Vec3::new(-0.0, 5.0, 6.0),
        solidcore::Vec3::new(7.0, 8.0, 9.0),
    );
    assert_eq!(
        record_pair((a, b)),
        (Vec3::new(1.0, 2.0, 3.0), Vec3::new(-0.0, 5.0, 6.0))
    );
    assert_eq!(
        record_triple((a, b, c)),
        (
            Vec3::new(1.0, 2.0, 3.0),
            Vec3::new(-0.0, 5.0, 6.0),
            Vec3::new(7.0, 8.0, 9.0)
        )
    );
    assert_eq!(
        record_list(&[c, a]),
        vec![Vec3::new(7.0, 8.0, 9.0), Vec3::new(1.0, 2.0, 3.0)]
    );
    let mut d: k::Partials =
        [[solidcore::Vec3::ZERO; k::MAX_PARTIAL_ORDER + 1]; k::MAX_PARTIAL_ORDER + 1];
    d[1][2] = a;
    d[3][0] = c;
    let t = record_table(&d);
    assert_eq!(
        (t[1][2], t[2][1], t[3][0], t[0][3]),
        (
            Vec3::new(1.0, 2.0, 3.0),
            Vec3::ZERO,
            Vec3::new(7.0, 8.0, 9.0),
            Vec3::ZERO
        )
    );
}

/// Every record kind has its kernel kind, so the conversion never refuses: a malformed record
/// (knots that do not match its points, a NaN) is carried over as it is, and the kernel's
/// evaluator is what refuses it.
#[test]
fn a_malformed_record_is_carried_over_for_the_kernel_to_refuse() {
    let bad_spline = Surface::Spline {
        u_knots: vec![0.0, 1.0],
        v_knots: vec![0.0, 0.0, 1.0, 1.0],
        points: vec![vec![[0.0; 3]; 2]; 2],
    };
    assert_eq!(record_surface(surface(&bad_spline)), bad_spline);
    assert!(SurfaceEval::new(&surface(&bad_spline)).is_err());
    let bad_curve = Curve::Spline {
        knots: vec![0.0, 1.0],
        points: vec![[0.0; 3]; 3],
    };
    assert_eq!(record_curve(curve(&bad_curve)), bad_curve);
    assert!(CurveEval::new(&curve(&bad_curve)).is_err());
    let nan = Surface::Sphere {
        frame: Frame {
            origin: [f64::NAN, 0.0, 0.0],
            ..f()
        },
        radius: 1.0,
    };
    let Surface::Sphere { frame: g, .. } = record_surface(surface(&nan)) else {
        panic!("a sphere")
    };
    assert!(g.origin[0].is_nan());
}

/// Every planar profile segment kind goes to the kernel's of the same fields (checked against
/// literals) and back unchanged; loops and regions keep their segments in order.
#[test]
fn planar_profiles_become_the_kernel_kinds_and_come_back() {
    let (a, b, c, d) = (
        Vec2::new(1.0, 2.0),
        Vec2::new(-3.5, 0.25),
        Vec2::new(0.0, -0.0),
        Vec2::new(7.0, 1e-300),
    );
    let ka = |v: Vec2| solidcore::math::Vec2::new(v.x, v.y);
    let mut p = [Vec2::new(0.0, 0.0); MAX_BEZIER_CTRL];
    p[0] = a;
    p[1] = b;
    p[2] = c;
    let segs = [
        (Seg2::Line { a, b }, kp::Seg2::Line { a: ka(a), b: ka(b) }),
        (
            Seg2::Arc {
                center: c,
                radius: 2.5,
                start: -0.5,
                sweep: -3.0,
            },
            kp::Seg2::Arc {
                center: ka(c),
                radius: 2.5,
                start: -0.5,
                sweep: -3.0,
            },
        ),
        (
            Seg2::Cubic {
                p0: a,
                p1: b,
                p2: c,
                p3: d,
            },
            kp::Seg2::Cubic {
                p0: ka(a),
                p1: ka(b),
                p2: ka(c),
                p3: ka(d),
            },
        ),
        (
            Seg2::Conic {
                a,
                apex: b,
                b: c,
                w: 0.75,
            },
            kp::Seg2::Conic {
                a: ka(a),
                apex: ka(b),
                b: ka(c),
                w: 0.75,
            },
        ),
        (
            Seg2::Bezier { n: 2, p },
            kp::Seg2::Bezier { n: 2, p: p.map(ka) },
        ),
    ];
    for (r, want) in &segs {
        assert_eq!(&seg2(r), want);
        assert_eq!(&record_seg2(&seg2(r)), r);
    }
    assert_eq!(record_vec2(vec2(d)), d);
    let outer = Loop2 {
        segs: segs.iter().map(|(r, _)| *r).collect(),
    };
    let hole = Loop2 {
        segs: vec![segs[0].0, segs[1].0],
    };
    assert_eq!(
        loop2(&outer).segs,
        segs.iter().map(|(_, k)| *k).collect::<Vec<_>>()
    );
    assert_eq!(record_loop2(&loop2(&outer)), outer);
    let region = Region2 {
        outer,
        holes: vec![hole.clone(), hole],
    };
    let k = region2(&region);
    assert_eq!(k.holes.len(), 2);
    assert_eq!(record_region2(&k), region);
}

// ---- #2357: a ball patch is made from the record only once it is certified on its sphere ----

#[allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers
mod ball_patch {
    use solidcore::blend::fillet_executive::session::FilletSession;
    use solidcore::blend::fillet_solver::p1::FilletSurfaceGeneratorType;
    use solidcore::curve::curve::p2::CurveKind;
    use solidcore::kernel::{KernelOptions, UpstreamFixes};
    use solidcore::math::Point3;
    use solidcore::nurbs::b_spline_curve::BSplineCurve;
    use solidcore::nurbs::math::{NOW, NOZ};
    use solidcore::primitive::api_primitives::create_box;
    use solidcore::query::api_queries::{get_edges, get_faces};
    use solidcore::surface::ball_patch::BallPatch;
    use solidcore::surface::surface::p2::SurfaceKind;
    use solidcore::topology::brep::p5::{stitch_and_orient, StitchAndOrientOptions};
    use solidcore::topology::edge;
    use solidcore::topology::face::p2 as face2;
    use solidcore::topology::topology_traverser::face_edges;
    use solidcore::Kernel;
    use solidcore_convert::ball_patch;
    use solidcore_record::brep::{
        Brep, Curve, Curve2, Edge, Face, Loop, OrientedEdge, Surface, VERSION,
    };

    /// A 40 x 30 x 60 box, every edge rounded r 1.5: eight Ball corners.
    fn rounded_box() -> (Kernel, solidcore::Brep) {
        let mut k = Kernel::new(&KernelOptions::default());
        k.fixes = UpstreamFixes::all();
        let mut b = create_box(&mut k, Point3::new(0.0, 0.0, 0.0), 40.0, 30.0, 60.0).unwrap();
        let edges = get_edges(&b).unwrap();
        let xs = vec![FilletSurfaceGeneratorType::Circular; edges.len()];
        {
            let mut s = FilletSession::new(&mut k, &mut b).unwrap();
            s.prepare_per_edge_constant_radius(&mut k, &edges, &[1.5; 12], &xs, 1, 1.0)
                .unwrap();
            s.do_fillet(&mut k).unwrap();
        }
        stitch_and_orient(&mut k, &mut b, &StitchAndOrientOptions::default()).unwrap();
        (k, b)
    }

    /// A kernel NURBS curve as the record's: its full knot vector and homogeneous points.
    fn record_nurbs_curve(c: &BSplineCurve) -> Curve {
        let (unique, mults) = c.knots(None).unwrap();
        let knots = unique
            .iter()
            .zip(&mults)
            .flat_map(|(&k, &m)| std::iter::repeat_n(k, m))
            .collect();
        let homogeneous = c
            .control_points()
            .unwrap()
            .iter()
            .map(|p| {
                let w = if p.w == NOW { 1.0 } else { p.w };
                [p.x, p.y, if p.z == NOZ { 0.0 } else { p.z }, w]
            })
            .collect();
        Curve::RationalSpline { knots, homogeneous }
    }

    fn record_curve(c: &solidcore::curve::curve::Curve) -> Curve {
        match c.kind() {
            CurveKind::BSpline(b) => record_nurbs_curve(b),
            CurveKind::Ellipse(e) => record_nurbs_curve(e.base()),
            CurveKind::Line(l) => record_nurbs_curve(l.base()),
            other => panic!("a loop curve of kind {other:?}"),
        }
    }

    /// The record of the box's ball patches: each face its `BallPatch` label over its net, its
    /// loop the face's edges with their curves (the rest of the record is not read), and the
    /// kernel's patches in the same order.
    fn ball_record(b: &solidcore::Brep) -> (Brep, Vec<BallPatch>) {
        let mut rec = Brep {
            version: VERSION,
            tolerance: 1e-6,
            regions: Vec::new(),
            shells: Vec::new(),
            faces: Vec::new(),
            edges: Vec::new(),
            vertices: Vec::new(),
        };
        let mut patches = Vec::new();
        for f in get_faces(b).unwrap() {
            let s = face2::surface(b, f).unwrap().unwrap();
            let SurfaceKind::BallPatch(bp) = s.kind() else {
                continue;
            };
            let n = bp.net().gw_nurb().unwrap();
            let w = |c: &solidcore::nurbs::math::CPoint| if c.w == NOW { 1.0 } else { c.w };
            let net = Surface::RationalSpline {
                u_knots: n.u_knots().to_vec(),
                v_knots: n.v_knots().to_vec(),
                homogeneous: n
                    .rows()
                    .map(|r| r.iter().map(|c| [c.x, c.y, c.z, w(c)]).collect())
                    .collect(),
            };
            let mut edges = Vec::new();
            for e in face_edges(b, f).unwrap() {
                let curve = record_curve(edge::curve(b, e).unwrap().unwrap());
                let index = u32::try_from(rec.edges.len()).unwrap();
                rec.edges.push(Edge {
                    curve,
                    nurbs: None,
                    interval: [0.0, 1.0],
                    vertices: [0, 0],
                    tolerance: 1e-6,
                    uses: Vec::new(),
                });
                edges.push(OrientedEdge {
                    edge: index,
                    forward: true,
                    pcurve: Curve2::Spline {
                        knots: vec![0.0, 0.0, 1.0, 1.0],
                        points: vec![[0.0, 0.0], [1.0, 1.0]],
                    },
                });
            }
            let c = bp.centre();
            rec.faces.push(Face {
                surface: Surface::BallPatch {
                    center: [c.x, c.y, c.z],
                    radius: bp.radius(),
                },
                nurbs: Some(net),
                domain: [[0.0, 1.0], [0.0, 1.0]],
                loops: vec![Loop::Edges { edges }],
            });
            patches.push(bp.copy());
        }
        (rec, patches)
    }

    fn same_bits(a: &BallPatch, b: &BallPatch) -> bool {
        let bits = |p: Point3| [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()];
        let (na, nb) = (
            a.net().control_point_net().unwrap(),
            b.net().control_point_net().unwrap(),
        );
        bits(a.centre()) == bits(b.centre())
            && a.radius().to_bits() == b.radius().to_bits()
            && na
                .points
                .iter()
                .map(|&p| bits(p))
                .eq(nb.points.iter().map(|&p| bits(p)))
            && na
                .weights
                .iter()
                .map(|w| w.to_bits())
                .eq(nb.weights.iter().map(|w| w.to_bits()))
    }

    /// The eight corners (each centre r in from its three faces, R = 1.5) come back from the
    /// record as certified ball patches, bit for bit the kernel's.
    #[test]
    fn the_rounded_box_s_corners_come_back_certified() {
        let (_k, b) = rounded_box();
        let (rec, patches) = ball_record(&b);
        assert_eq!(rec.faces.len(), 8);
        for (i, want) in patches.iter().enumerate() {
            let got = ball_patch(&rec, i).unwrap();
            assert!(same_bits(&got, want), "face {i}");
            let c = got.centre();
            for (x, side) in [(c.x, 40.0), (c.y, 30.0), (c.z, 60.0)] {
                let near = (x - 1.5).abs().min((x - (side - 1.5)).abs());
                assert!(near <= 1e-12, "face {i}: {c:?}");
            }
            assert_eq!(got.radius(), 1.5);
        }
    }

    /// A sphere 1e-9 larger, a NaN centre, a loop curve moved 1e-7 are refused, naming the face
    /// and, for the sphere, the coefficient off it against its bound.
    #[test]
    fn a_ball_patch_off_its_sphere_is_refused() {
        let (_k, b) = rounded_box();
        let (rec, _) = ball_record(&b);
        for i in 0..rec.faces.len() {
            let mut off = rec.clone();
            if let Surface::BallPatch { radius, .. } = &mut off.faces[i].surface {
                *radius += 1e-9;
            }
            let why = ball_patch(&off, i).unwrap_err().to_string();
            assert!(
                why.contains(&format!("face {i}")) && why.contains("bound"),
                "{why}"
            );
            let mut off = rec.clone();
            if let Surface::BallPatch { center, .. } = &mut off.faces[i].surface {
                center[0] = f64::NAN;
            }
            assert!(ball_patch(&off, i).is_err());
        }
        let mut off = rec.clone();
        let Loop::Edges { edges } = &off.faces[0].loops[0] else {
            panic!("an edge loop");
        };
        let e = edges[0].edge as usize;
        if let Curve::RationalSpline { homogeneous, .. } = &mut off.edges[e].curve {
            for h in homogeneous.iter_mut() {
                let w = h[3];
                for x in &mut h[..3] {
                    *x += 1e-7 * w;
                }
            }
        }
        let why = ball_patch(&off, 0).unwrap_err().to_string();
        assert!(why.contains("loop edge"), "{why}");
    }

    /// Far from the origin the certificate is as tight: the corners moved to about
    /// (1e5, -7e4, 3e4) are accepted, and a sphere 1e-9 larger there is still refused.
    #[test]
    fn a_ball_patch_far_from_the_origin_is_accepted() {
        let (_k, b) = rounded_box();
        let (mut rec, _) = ball_record(&b);
        let t = [1e5, -7e4, 3e4];
        let shift = |h: &mut [f64; 4]| {
            let w = h[3];
            for (x, d) in h.iter_mut().zip(t) {
                *x += d * w;
            }
        };
        for f in &mut rec.faces {
            if let Surface::BallPatch { center, .. } = &mut f.surface {
                for (c, d) in center.iter_mut().zip(t) {
                    *c += d;
                }
            }
            if let Some(Surface::RationalSpline { homogeneous, .. }) = &mut f.nurbs {
                homogeneous.iter_mut().flatten().for_each(shift);
            }
        }
        for e in &mut rec.edges {
            if let Curve::RationalSpline { homogeneous, .. } = &mut e.curve {
                homogeneous.iter_mut().for_each(shift);
            }
        }
        for i in 0..rec.faces.len() {
            ball_patch(&rec, i).unwrap_or_else(|e| panic!("face {i}: {e}"));
            let mut off = rec.clone();
            if let Surface::BallPatch { radius, .. } = &mut off.faces[i].surface {
                *radius += 1e-9;
            }
            assert!(ball_patch(&off, i).is_err(), "face {i}");
        }
    }
}
