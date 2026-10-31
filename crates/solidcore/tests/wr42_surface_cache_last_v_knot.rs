//! WR-42: a surface cache must include its sole internal C0 v knot when opted in.
#![allow(clippy::unwrap_used, clippy::indexing_slicing)]

use solidcore::geom::extent2d::Extent2d;
use solidcore::geom::geom_utility::p2::ContinuityType;
use solidcore::nurbs::b_spline_surface::BSplineSurface;
use solidcore::nurbs::math::{CPoint, NOW};
use solidcore::nurbs::srf_basic::Surface as NurbSurface;
use solidcore::surface::surface::Surface;
use solidcore::surface::surface_cache::p1::{SurfaceCache, SurfaceCacheOptions};
use solidcore::Kernel;

#[test]
fn sole_internal_c0_v_knot_is_opt_in() {
    // Two bilinear patches meet at v = 0.5 with different v tangents.
    let net = [0.0, 1.0]
        .into_iter()
        .map(|u| {
            vec![
                CPoint::new(u, 0.0, 0.0, NOW),
                CPoint::new(u, 0.5, 0.0, NOW),
                CPoint::new(u, 1.0, 1.0, NOW),
            ]
        })
        .collect();
    let nurb = NurbSurface::new(
        net,
        1,
        1,
        vec![0.0, 0.0, 1.0, 1.0],
        vec![0.0, 0.0, 0.5, 1.0, 1.0],
    )
    .unwrap();
    let surface = Surface::from(BSplineSurface::from_nurb(nurb, false));
    let domain = Extent2d::from_bounds(0.0, 0.0, 1.0, 1.0);
    let mut cache = SurfaceCache::<()>::new(surface, domain, SurfaceCacheOptions::default());
    let mut k = Kernel::default();
    cache.restore_continuities(&k).unwrap();
    assert_eq!(cache.v_continuities().len(), 3);
    assert_eq!(cache.v_continuities()[1], ContinuityType::C0);
    assert!(!k.fixes.surface_cache_last_v_knot);
    assert_eq!(
        cache.surface_continuity(&k, &domain).unwrap(),
        ContinuityType::CInfinity
    );

    k.fixes.surface_cache_last_v_knot = true;
    assert_eq!(
        cache.surface_continuity(&k, &domain).unwrap(),
        ContinuityType::C0
    );
    for (v0, v1) in [(0.0, 0.25), (0.75, 1.0)] {
        let away = Extent2d::from_bounds(0.0, v0, 1.0, v1);
        assert_eq!(
            cache.surface_continuity(&k, &away).unwrap(),
            ContinuityType::CInfinity
        );
    }
    k.fixes.surface_cache_last_v_knot = false;
    assert_eq!(
        cache.surface_continuity(&k, &domain).unwrap(),
        ContinuityType::CInfinity
    );
}

#[test]
fn tracer_checks_the_last_internal_v_entry_when_opted_in() {
    use solidcore::intersect::surface_tracer::{SurfaceTrace, SurfaceTracer};

    // The tracer normally only reports a C0 diagnostic and continues. An oversized
    // cached array makes the extra visit observable through the query's bounds
    // error, without changing that production diagnostic into a new failure.
    let mut k = Kernel::default();
    let surface = Surface::from(BSplineSurface::from_nurb(
        NurbSurface::new(
            vec![
                vec![
                    CPoint::new(0.0, 0.0, 0.0, NOW),
                    CPoint::new(0.0, 1.0, 0.0, NOW),
                ],
                vec![
                    CPoint::new(1.0, 0.0, 0.0, NOW),
                    CPoint::new(1.0, 1.0, 0.0, NOW),
                ],
            ],
            1,
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
        )
        .unwrap(),
        false,
    ));
    let domain = Extent2d::from_bounds(0.0, 0.0, 1.0, 1.0);
    let mut cache = SurfaceCache::<solidcore::surface::surface::p4::NoObject>::new(
        surface.clone(),
        domain,
        SurfaceCacheOptions::default(),
    );
    cache.set_continuities(Vec::new(), vec![ContinuityType::CInfinity; 4]);
    let mut tracer = SurfaceTracer::new(&k, &surface, &domain).unwrap();
    let legacy_error = tracer.do_trace(&k, &mut cache, None, None).unwrap_err();
    assert!(!legacy_error
        .to_string()
        .contains("continuity array longer than the knots"));
    k.fixes.surface_cache_last_v_knot = true;
    let error = tracer.do_trace(&k, &mut cache, None, None).unwrap_err();
    assert!(error
        .to_string()
        .contains("continuity array longer than the knots"));
}
