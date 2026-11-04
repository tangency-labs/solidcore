//! #2048: `check_brep` of a body with a face on an `OffsetSurface` returned
//! `Unported { upstream: "SmFace handle of a surface's face", .. }`.
//! `SmOffsetSurface::AssertValid` test 2 runs `GlobalSurfaceIntersect` of the surface
//! with itself, and upstream's `smsurf_GetSurfaceCache` reaches the surface's face
//! through the surface (`GetFace()`); the port's face check now hands the face-aware
//! cache context to that solve.
//!
//! The input: a quarter cylinder (radius 10, length 20) as a rational B-spline face,
//! offset by 1 to a quarter cylinder of radius 11 on one `OffsetSurface` face, of area
//! `11 * pi / 2 * 20`.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::indexing_slicing)] // test helpers

use solidcore::geom::object::ObjectRef;
use solidcore::geom::tol::USE_DEFAULT;
use solidcore::mass::integrator::properties::brep_area;
use solidcore::nurbs::b_spline_surface::p1::BSplineSurface;
use solidcore::nurbs::knots::Side;
use solidcore::nurbs::math::CPoint;
use solidcore::offset::offset_executive::p2::{offset_brep_simple, OffsetOptions};
use solidcore::surface::offset_surface::OffsetSurface;
use solidcore::surface::surface::{Surface, SurfaceKind};
use solidcore::topology::brep::local_operation::local_operation_simple;
use solidcore::topology::brep::p2::ArenaNotifyState;
use solidcore::topology::brep::p3::new_brep;
use solidcore::topology::brep::p8::{create_face_from_surface, FaceBuildEnv};
use solidcore::topology::brep::validation::{check_brep, CheckLevel};
use solidcore::topology::loopuse::notify_ref;
use solidcore::topology::topology_traverser::Marks;
use solidcore::topology::{Brep, FaceId};
use solidcore::{EntityRef, Kernel};

/// The issue's quarter cylinder: u degree 2 (rational quarter circle), v degree 1.
fn quarter_cylinder(k: &mut Kernel) -> Brep {
    let c = std::f64::consts::FRAC_1_SQRT_2;
    let pw: Vec<Vec<CPoint>> = [(10.0, 0.0, 1.0), (10.0, 10.0, c), (0.0, 10.0, 1.0)]
        .iter()
        .map(|&(x, y, w)| {
            vec![
                CPoint::new(x * w, y * w, 0.0, w),
                CPoint::new(x * w, y * w, 20.0 * w, w),
            ]
        })
        .collect();
    let s = Surface::from(
        BSplineSurface::from_data(
            2,
            1,
            2,
            1,
            5,
            3,
            &pw,
            &[0., 0., 0., 1., 1., 1.],
            &[0., 0., 1., 1.],
        )
        .unwrap(),
    );
    face_on_surface(k, s)
}

fn face_on_surface(k: &mut Kernel, s: Surface) -> Brep {
    let d = s.natural_uv_domain_with_nurb(k).unwrap();
    face_on_surface_domain(k, s, d)
}

fn face_on_surface_domain(
    k: &mut Kernel,
    s: Surface,
    d: solidcore::geom::extent2d::Extent2d,
) -> Brep {
    let mut b = new_brep(k, USE_DEFAULT, USE_DEFAULT).unwrap();
    {
        let this_face = |f: FaceId| -> ObjectRef { notify_ref(EntityRef::Face(f)) };
        let mut state = ArenaNotifyState::new(&b);
        let mut env = state.env(None);
        let mut marks = Marks::default();
        let mut fenv = FaceBuildEnv {
            env: &mut env,
            marks: &mut marks,
            this_face: &this_face,
        };
        create_face_from_surface(k, &mut b, Some(s), &d, &mut fenv).unwrap();
    }
    b.refresh_cached_topology();
    b
}

#[test]
fn check_brep_answers_on_an_offset_surface_face() {
    let mut k = Kernel::default();
    let mut b = quarter_cylinder(&mut k);
    assert_eq!(
        check_brep(&mut k, &mut b, CheckLevel::Level2)
            .unwrap()
            .assert_valid,
        Some(true)
    );

    let mut r = offset_brep_simple(
        &mut k,
        &mut b,
        OffsetOptions {
            distance: 1.0,
            extend_convex_edges: true,
            merge_results: true,
        },
    )
    .unwrap();
    r.refresh_cached_topology();

    let faces: Vec<FaceId> = r.faces().ids().collect();
    assert_eq!(faces.len(), 1);
    let surface = r.face(faces[0]).unwrap().surface.as_ref().unwrap();
    assert!(matches!(surface.kind(), SurfaceKind::OffsetSurface(_)));
    let area = brep_area(&k, &mut r, 1e-6).unwrap();
    let expected = 11.0 * std::f64::consts::FRAC_PI_2 * 20.0;
    assert!((area - expected).abs() < 1e-3, "area {area} vs {expected}");

    let c = check_brep(&mut k, &mut r, CheckLevel::Level2).expect("check_brep answers");
    assert_eq!(c.assert_valid, Some(true), "{:?}", c.reports);
    assert_eq!(c.failures().count(), 0);
}

#[test]
fn extended_offset_discards_only_unused_support_continuations() {
    for distance in [-0.1, 0.1] {
        let mut k = Kernel {
            fixes: solidcore::kernel::UpstreamFixes::all(),
            ..Kernel::default()
        };
        // WR-736's normal-offset corner image must not consume LocalOperation's
        // bounding-box extension distance as the replacement's displacement.
        assert!(k.fixes.offset_lamina_corner_at_vertex_image);
        let base = rising_eight_turns();
        let mut b = non_iso_trimmed_patch(&mut k, base.clone());
        assert!(check_brep(&mut k, &mut b, CheckLevel::Level2)
            .unwrap()
            .is_valid());
        let faces: Vec<_> = b.faces().ids().collect();
        let counts = (b.faces().len(), b.edges().len(), b.vertices().len());
        // The retained boundary extends beyond this replacement's original
        // support, so preserving its paired UV curves cannot restore that
        // support. Radial rulings are linear: their continuation is exact.
        let mut target = base.clone();
        let mut target_domain = solidcore::geom::extent2d::Extent2d::from_points(
            solidcore::math::Point2::new(0.0, 0.25),
            solidcore::math::Point2::new(32.0, 0.75),
        );
        target.trim_with_domain(&k, &mut target_domain).unwrap();
        let offset = Surface::from(OffsetSurface::new(&k, distance, target, true, 0).unwrap());
        local_operation_simple(&mut k, &mut b, &faces, vec![offset]).unwrap();
        assert_eq!(
            counts,
            (b.faces().len(), b.edges().len(), b.vertices().len())
        );
        let check = check_brep(&mut k, &mut b, CheckLevel::Level2).unwrap();
        assert!(check.is_valid(), "{distance}: {:?}", check.lines());
        let face_id = b.faces().ids().next().unwrap();
        let trimmed = solidcore::topology::face::p1::calculate_uv_domain_from_uv_trim_curves(
            &k, &mut b, face_id, true, false, false,
        )
        .unwrap()
        .domain;
        let surface = b.face(face_id).unwrap().surface.as_ref().unwrap();
        assert!(matches!(surface.kind(), SurfaceKind::OffsetSurface(_)));
        let retained = surface.natural_uv_domain_with_nurb(&k).unwrap();
        // ShrinkGeometry retains nearby existing bounds within one thousandth
        // of the trim-domain diagonal. It discards the large unused extension.
        assert!(
            retained.are_equal(&trimmed, trimmed.size().length() / 1000.0),
            "{retained:?} / {trimmed:?}"
        );
        for u in [1.0, 8.0, 16.0, 24.0, 31.0] {
            for v in [0.4, 0.5, 0.6] {
                let uv = solidcore::math::Point2::new(u, v);
                let expected = base.evaluate_point(&k, uv).unwrap()
                    + base
                        .evaluate_normal(&k, uv, Side::Left, Side::Left)
                        .unwrap()
                        * distance;
                assert!(surface.evaluate_point(&k, uv).unwrap().distance(expected) < 1e-8);
            }
        }
    }
}

#[test]
fn local_operation_matches_normals_on_the_same_cylinder_quadrant() {
    for reversed_parameter in [false, true] {
        let mut k = Kernel {
            fixes: solidcore::kernel::UpstreamFixes::all(),
            ..Kernel::default()
        };
        let base = Surface::from(
            solidcore::primitive::api_primitives::create_cylinder_surface(
                &k,
                solidcore::math::Vec3::ZERO,
                10.0,
                20.0,
            )
            .unwrap(),
        );
        let natural = base.natural_uv_domain_with_nurb(&k).unwrap();
        let quadrant = solidcore::geom::extent2d::Extent2d::from_points(
            natural.evaluate(0.75, 0.0),
            natural.max(),
        );
        let reference = base.evaluate_point(&k, quadrant.mid()).unwrap();
        let mut b = face_on_surface_domain(&mut k, base.clone(), quadrant);
        assert!(check_brep(&mut k, &mut b, CheckLevel::Level2)
            .unwrap()
            .is_valid());
        let faces: Vec<_> = b.faces().ids().collect();
        let mut target = base.create_offset_surface(&k, -0.1, 1e-6).unwrap().unwrap();
        if reversed_parameter {
            target
                .reverse(&k, solidcore::nurbs::b_spline_surface::SurfParamType::U)
                .unwrap();
        }
        local_operation_simple(&mut k, &mut b, &faces, vec![target]).unwrap();
        assert_eq!(
            (b.faces().len(), b.edges().len(), b.vertices().len()),
            (1, 4, 4)
        );
        assert!(check_brep(&mut k, &mut b, CheckLevel::Level2)
            .unwrap()
            .is_valid());
        let face = b.face(b.faces().ids().next().unwrap()).unwrap();
        let surface = face.surface.as_ref().unwrap();
        let uv = face.uv_domain.mid();
        let point = surface.evaluate_point(&k, uv).unwrap();
        let normal = surface
            .evaluate_normal(&k, uv, Side::Left, Side::Left)
            .unwrap();
        let outward = point.x * normal.x + point.y * normal.y;
        assert!(outward > 9.8, "{point:?} / {normal:?}");
        let expected =
            solidcore::math::Vec3::new(reference.x * 0.99, reference.y * 0.99, reference.z);
        assert!(point.distance(expected) < 1e-8, "{point:?} / {expected:?}");
        let area = brep_area(&k, &mut b, 1e-6).unwrap();
        let expected_area = 9.9 * std::f64::consts::FRAC_PI_2 * 20.0;
        assert!((area - expected_area).abs() < 1e-3, "{area}");
    }
}

/// Eight turns of an exact rational ruled patch. Four non-isoparametric UV
/// segments bound a sheet; no screw sweep or imported geometry is involved.
fn rising_eight_turns() -> Surface {
    let quarters = 32usize;
    let mut points = Vec::new();
    for i in 0..=2 * quarters {
        let a = i as f64 * std::f64::consts::FRAC_PI_4;
        let w = if i % 2 == 0 {
            1.0
        } else {
            std::f64::consts::FRAC_1_SQRT_2
        };
        let z = i as f64 * 2.5 / 8.0;
        points.push(vec![
            CPoint::new(
                8.646835306586816 * a.cos(),
                8.646835306586816 * a.sin(),
                z * w,
                w,
            ),
            CPoint::new(10.0 * a.cos(), 10.0 * a.sin(), (z - 0.78125) * w, w),
        ]);
    }
    let mut knots = vec![0.0; 3];
    for i in 1..quarters {
        knots.extend([i as f64, i as f64]);
    }
    knots.extend([quarters as f64; 3]);
    Surface::from(
        BSplineSurface::from_data(
            2 * quarters,
            1,
            2,
            1,
            2 * quarters + 3,
            3,
            &points,
            &knots,
            &[0., 0., 1., 1.],
        )
        .unwrap(),
    )
}
fn non_iso_trimmed_patch(k: &mut Kernel, s: Surface) -> Brep {
    use solidcore::curve::curve::Curve;
    use solidcore::math::Point3;
    use solidcore::nurbs::b_spline_curve::BSplineCurve;
    use solidcore::nurbs::b_spline_surface_n_lib::p2::OrientType;
    let domain = s.natural_uv_domain_with_nurb(k).unwrap();
    let mut b = new_brep(k, USE_DEFAULT, USE_DEFAULT).unwrap();
    let corners = [
        Point3::new(0.1, 0.2, 0.),
        Point3::new(31.9, 0.3, 0.),
        Point3::new(31.8, 0.8, 0.),
        Point3::new(0.2, 0.7, 0.),
    ];
    let uv: Vec<_> = (0..4)
        .map(|i| {
            Curve::from(
                BSplineCurve::create_line_segment(2, corners[i], corners[(i + 1) % 4]).unwrap(),
            )
        })
        .collect();
    let input = solidcore::topology::brep::p7::FaceCurvesInput {
        loops: vec![4],
        curves_3d: None,
        uv_curves: Some(uv),
        orientations: vec![OrientType::Same; 4],
        loop_points: vec![],
    };
    let mut attrs = solidcore::geom::attribute::AttributeStore::default();
    let mut records = solidcore::topology::topology::NotifyRecords::new();
    let mut stores = solidcore::topology::edge::NotificationStores {
        records: &mut records,
        attributes: &mut attrs,
        key_of: &notify_ref,
    };
    solidcore::topology::brep::p7::make_face_with_curves_with_kernel_context(
        k,
        &mut b,
        None,
        input,
        s,
        &domain,
        OrientType::Same,
        Some(true),
        &|f| notify_ref(EntityRef::Face(f)),
        &mut stores,
    )
    .unwrap();
    b.refresh_cached_topology();
    b
}

#[test]
fn extended_offset_preserves_paired_non_iso_parameter_curves() {
    let mut k = Kernel {
        fixes: solidcore::kernel::UpstreamFixes::all(),
        ..Kernel::default()
    };
    // Keep the uniform-offset lamina correction enabled: local replacements
    // must still construct their corners from their supplied supports.
    assert!(k.fixes.offset_lamina_corner_at_vertex_image);
    let base = rising_eight_turns();
    let mut b = non_iso_trimmed_patch(&mut k, base.clone());
    assert!(check_brep(&mut k, &mut b, CheckLevel::Level2)
        .unwrap()
        .is_valid());
    let faces: Vec<_> = b.faces().ids().collect();
    let counts = (b.faces().len(), b.edges().len(), b.vertices().len());
    let distance = -0.1;
    let offset = Surface::from(OffsetSurface::new(&k, distance, base.clone(), true, 0).unwrap());
    let result = local_operation_simple(&mut k, &mut b, &faces, vec![offset]);
    result.unwrap();
    assert_eq!(
        counts,
        (b.faces().len(), b.edges().len(), b.vertices().len())
    );
    let check = check_brep(&mut k, &mut b, CheckLevel::Level2).unwrap();
    assert!(check.is_valid(), "{:?}", check.lines());
    assert!(b
        .edges()
        .iter()
        .all(|(_, e)| e.tolerance <= 20.0 * b.tolerance()));
    let surface = b
        .face(b.faces().ids().next().unwrap())
        .unwrap()
        .surface
        .as_ref()
        .unwrap();
    for u in [1.0, 8.0, 16.0, 24.0, 31.0] {
        for v in [0.4, 0.5, 0.6] {
            let uv = solidcore::math::Point2::new(u, v);
            let expected = base.evaluate_point(&k, uv).unwrap()
                + base
                    .evaluate_normal(&k, uv, Side::Left, Side::Left)
                    .unwrap()
                    * distance;
            assert!(surface.evaluate_point(&k, uv).unwrap().distance(expected) < 1e-8);
        }
    }
}
