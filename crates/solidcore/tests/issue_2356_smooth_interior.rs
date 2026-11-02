//! Deterministic measurements for #2356; these are not interval certificates.
//!
//! The top is z = 0 before x = 5, z = 0.008 (x - 5)^3 afterwards:
//! cubic C2 at the simple knot, with a jump in its third derivative. A ball
//! between it and y = 0 crosses that knot exactly at leader parameter s = 0.5.
//! The analytic control replaces the top with z = 0. No fitted B-rep is needed.
#![allow(clippy::unwrap_used, clippy::panic, clippy::indexing_slicing)]

use solidcore::blend::fillet_section::section::{GnSectionFsg, GN_THUMBWEIGHT_EQUAL_AREA};
use solidcore::blend::fillet_section::smooth_interior::correction_kernels;
use solidcore::kernel::Kernel;
use solidcore::math::{Point2, Point3, Vec3};
use solidcore::nurbs::b_spline_curve::{BSplineCurve, KnotType};
use solidcore::nurbs::b_spline_surface::{
    BSplineSurface, BSplineSurfaceForm, CanonicalSurfaceData, SurfParamType,
};
use solidcore::nurbs::b_spline_surface_n_lib::p2::EvaluateFlags;
use solidcore::nurbs::knots::Side;
use solidcore::surface::procedural_blend::{BlendDefinitionInput, ProceduralBlend};
use solidcore::surface::surface::{Surface, SurfaceKind};

fn line(dim: u64, a: Point3, b: Point3) -> BSplineCurve {
    BSplineCurve::create_line_segment(dim, a, b).unwrap()
}

fn fixture(curved: bool) -> (Kernel, ProceduralBlend) {
    let mut k = Kernel::default();
    k.fixes.fillet_gn_section_generator = true;
    k.fixes.fillet_round_as_definition = true;
    // Sided definition evaluation is provided by #2322, not #2274 alone.
    k.fixes.fillet_round_on_round_as_definition = true;
    let wall = Surface::from(
        BSplineSurface::create_bilinear_surface(
            Point3::ZERO,
            Point3::new(10.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 10.0),
            Point3::new(10.0, 0.0, 10.0),
        )
        .unwrap(),
    );
    let mut cp = Vec::new();
    for (i, x) in [0.0, 10.0 / 6.0, 5.0, 50.0 / 6.0, 10.0]
        .into_iter()
        .enumerate()
    {
        for y in [0.0, 10.0] {
            cp.push(Point3::new(x, y, if curved && i == 4 { 1.0 } else { 0.0 }));
        }
    }
    let top = Surface::from(
        BSplineSurface::create_canonical(&CanonicalSurfaceData {
            u_degree: 3,
            v_degree: 1,
            control_points: &cp,
            form: BSplineSurfaceForm::Unspecified,
            u_multiplicities: &[4, 1, 4],
            v_multiplicities: &[2, 2],
            u_knots: &[0.0, 0.5, 1.0],
            v_knots: &[0.0, 1.0],
            knot_type: KnotType::Unspecified,
            weights: None,
            uv_domain: None,
        })
        .unwrap(),
    );
    let leader = line(3, Point3::new(0.0, 2.0, 2.0), Point3::new(10.0, 2.0, 2.0));
    let seed = line(2, Point3::new(0.0, 0.2, 0.0), Point3::new(1.0, 0.2, 0.0));
    let dummy_net =
        BSplineSurface::create_bilinear_surface(Point3::ZERO, Vec3::X, Vec3::Y, Vec3::X + Vec3::Y)
            .unwrap();
    let pb = ProceduralBlend::create(
        &k,
        &dummy_net,
        &BlendDefinitionInput {
            supports: [&wall, &top],
            sides: [-1.0, 1.0],
            offsets: [2.0, 2.0],
            law: None,
            order: 2,
            thumbweight: GN_THUMBWEIGHT_EQUAL_AREA,
            leader: &leader,
            seeds: [&seed, &seed],
            bound: 1e-5,
        },
        |n, w| GnSectionFsg::is_valid(n as usize, w),
    )
    .unwrap();
    (k, pb)
}

fn normal(d: &[Vec<Vec3>]) -> Vec3 {
    let n = d[1][0].cross(d[0][1]);
    assert!(n.length() > 1e-8, "singular sample");
    n / n.length()
}

fn normal_jump(a: &[Vec<Vec3>], b: &[Vec<Vec3>]) -> f64 {
    // atan2 retains accuracy near zero; acos(dot) loses the 1e-9 gate.
    let (a, b) = (normal(a), normal(b));
    a.cross(b).length().atan2(a.dot(b))
}

fn s_curvature(d: &[Vec<Vec3>]) -> f64 {
    normal(d).dot(d[2][0]) / d[1][0].dot(d[1][0])
}

#[test]
fn corrected_station_jets_cancel_measured_jumps_including_equal_area_rows() {
    let (k, pb) = fixture(true);
    let lambda = 0.15;
    for j in 1..32 {
        let t = f64::from(j) / 32.0;
        let eps = lambda * 4.0 * t * (1.0 - t);
        let eps_t = lambda * 4.0 * (1.0 - 2.0 * t);
        let mut a = pb
            .evaluate_definition_sided(&k, 0.5, t, 2, 2, true)
            .unwrap();
        let mut b = pb
            .evaluate_definition_sided(&k, 0.5, t, 2, 2, false)
            .unwrap();
        // Jumps of the complete rational section, including w(s)'s series.
        let delta = b[1][0] - a[1][0];
        let delta_t = b[1][1] - a[1][1];
        let delta_tt = b[1][2] - a[1][2];
        let delta2 = b[2][0] - a[2][0];
        let delta2_t = b[2][1] - a[2][1];
        for (d, before) in [(&mut a, true), (&mut b, false)] {
            let (g, g2) = correction_kernels(0.0, eps, before).unwrap();
            d[0][0] = d[0][0] - delta * g.value - delta2 * g2.value;
            d[1][0] = d[1][0] - delta * g.first - delta2 * g2.first;
            d[2][0] = d[2][0] - delta * g.second - delta2 * g2.second;
            // G(0,eps)=-3 eps/16; G2(0,eps)=0 identically in t.
            d[0][1] += (delta_t * eps + delta * eps_t) * (3.0 / 16.0);
            d[1][1] =
                d[1][1] - delta_t * g.first - delta2_t * g2.first + delta2 * (3.0 * eps_t / 32.0);
            d[0][2] +=
                (delta_tt * eps + delta_t * (2.0 * eps_t) - delta * (8.0 * lambda)) * (3.0 / 16.0);
        }
        assert!((a[0][0] - b[0][0]).length() <= 1e-9);
        assert!((a[1][0] - b[1][0]).length() <= 1e-9);
        assert!((a[2][0] - b[2][0]).length() <= 1e-9);
        assert!((a[1][1] - b[1][1]).length() <= 1e-9);
        assert!((a[0][2] - b[0][2]).length() <= 1e-9);
        assert!(normal_jump(&a, &b) <= 1e-9);
        assert!((s_curvature(&a) - s_curvature(&b)).abs() <= 1e-9);
    }
}

#[test]
fn built_net_reproduces_the_crease_at_a_support_crossing() {
    for curved in [false, true] {
        let (k, pb) = fixture(curved);
        let mut rail0 = line(3, Point3::new(0.0, 0.0, 2.0), Point3::new(10.0, 0.0, 2.0));
        let mut rail1 = line(3, Point3::new(0.0, 2.0, 0.0), Point3::new(10.0, 2.0, 1.0));
        let mut uv0 = pb.seed(0).unwrap().clone();
        let mut uv1 = pb.seed(1).unwrap().clone();
        let mut generator = GnSectionFsg::new(2, GN_THUMBWEIGHT_EQUAL_AREA);
        let built = generator
            .build_with_radius(
                &k,
                pb.leader(),
                [&mut rail0, &mut rail1],
                1e-5,
                [pb.support(0).unwrap(), pb.support(1).unwrap()],
                [&mut uv0, &mut uv1],
                true,
                2.0,
            )
            .unwrap();
        let SurfaceKind::ProceduralBlend(built) = built.kind() else {
            panic!("expected definition with its built net")
        };
        let net = built.net();
        let (mut knots, _) = net.knots(SurfParamType::U, None).unwrap();
        let crossing = if curved {
            *knots
                .iter()
                .min_by(|a, b| ((**a - 0.5).abs()).total_cmp(&(**b - 0.5).abs()))
                .unwrap()
        } else {
            0.5
        };
        if !curved && !knots.contains(&0.5) {
            knots.push(0.5);
        }
        assert!((crossing - 0.5).abs() <= 1e-9, "missing support crossing");
        let mut max_jump = 0.0_f64;
        let mut measured_lines = 0;
        for &s in knots.iter().filter(|&&s| s > 0.0 && s < 1.0) {
            measured_lines += 1;
            for j in 0..=32 {
                let t = f64::from(j) / 32.0;
                let a = net
                    .evaluate_simple(Point2::new(s, t), 2, 2, true, true, false)
                    .unwrap();
                let b = net
                    .evaluate_simple(Point2::new(s, t), 2, 2, false, true, false)
                    .unwrap();
                assert!((a[0][0] - b[0][0]).length() <= 1e-9);
                let jump = normal_jump(&a, &b);
                max_jump = max_jump.max(jump);
                if !curved || (s - crossing).abs() > 1e-9 {
                    assert!(jump <= 1e-9, "unexpected crease at s={s}, t={t}: {jump:e}");
                }
                if (s - crossing).abs() <= 1e-9 {
                    for (d, before) in [(&a, false), (&b, true)] {
                        // NLib's FromLeft=true selects the span starting at the knot.
                        let definition = built
                            .evaluate_definition_sided(&k, s, t, 1, 1, before)
                            .unwrap();
                        assert!(
                            (d[1][0] - definition[1][0]).length() <= 1e-9,
                            "net must reproduce its definition's sided s-derivative"
                        );
                    }
                }
            }
        }
        println!("built net curved={curved}; crossing={crossing:e}; max normal jump={max_jump:e} rad; spans={}; measured lines={measured_lines}", generator.certified_spans);
        if curved {
            assert!(
                max_jump > 1e-6,
                "known net crease disappeared: enable the 1e-9 acceptance gate"
            );
        }
    }
}

#[test]
fn analytic_control_has_no_one_sided_crease() {
    let (k, pb) = fixture(false);
    for j in 0..=32 {
        let t = f64::from(j) / 32.0;
        let a = pb
            .evaluate_definition_sided(&k, 0.5, t, 2, 2, true)
            .unwrap();
        let b = pb
            .evaluate_definition_sided(&k, 0.5, t, 2, 2, false)
            .unwrap();
        assert!(normal_jump(&a, &b) <= 1e-9);
        assert!((s_curvature(&a) - s_curvature(&b)).abs() <= 1e-9);
    }
}

#[test]
fn c2_support_reproduces_definition_crease_without_losing_rail_contact() {
    let (k, pb) = fixture(true);
    let mut max_jump = 0.0_f64;
    let mut max_curvature_jump = 0.0_f64;
    for j in 0..=32 {
        let t = f64::from(j) / 32.0;
        let a = pb
            .evaluate_definition_sided(&k, 0.5, t, 2, 2, true)
            .unwrap();
        let b = pb
            .evaluate_definition_sided(&k, 0.5, t, 2, 2, false)
            .unwrap();
        assert!((a[0][0] - b[0][0]).length() <= 1e-9);
        max_jump = max_jump.max(normal_jump(&a, &b));
        max_curvature_jump = max_curvature_jump.max((s_curvature(&a) - s_curvature(&b)).abs());
        if j == 0 || j == 32 {
            let n = if j == 0 { -Vec3::Y } else { Vec3::Z };
            for d in [&a, &b] {
                assert!((d[0][0].x - 5.0).abs() <= 1e-9);
                assert!(d[0][0].dot(n).abs() <= 1e-9, "rail position");
                assert!(d[0][1].dot(n).abs() <= 1e-9, "rail tangent");
                assert!(d[0][2].dot(n).abs() <= 1e-9, "rail normal curvature");
            }
            assert!(normal_jump(&a, &b) <= 1e-9);
        }
    }
    println!(
        "C2 crossing: max normal jump={max_jump:e} rad; s-curvature jump={max_curvature_jump:e}"
    );
    // A known-defect regression, deliberately inverted until production correction
    // and its certificate land. It must stop passing when the defect is fixed.
    assert!(
        max_jump > 1e-6,
        "known crease disappeared: turn this into the 1e-9 acceptance gate"
    );
    assert!(max_curvature_jump > 1e-6);
}

#[test]
fn support_is_c2_but_not_c3_at_the_crossing() {
    let (k, pb) = fixture(true);
    let top = pb.support(1).unwrap();
    let mut sides = Vec::new();
    for side in [Side::Left, Side::Right] {
        sides.push(
            top.evaluate_simple(
                &k,
                Point2::new(0.5, 0.2),
                3,
                0,
                EvaluateFlags {
                    u_side: side,
                    ..EvaluateFlags::default()
                },
            )
            .unwrap(),
        );
    }
    for i in 0..=2 {
        assert!((sides[0][i][0] - sides[1][i][0]).length() <= 1e-9);
    }
    assert!((sides[0][3][0].z - sides[1][3][0].z).abs() > 1.0);
}
