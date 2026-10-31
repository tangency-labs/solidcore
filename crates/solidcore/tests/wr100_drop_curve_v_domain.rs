//! WR-100: `DropCurve` must distinguish natural U and V domain degeneracy.
#![allow(clippy::unwrap_used)]
use solidcore::curve::curve::Curve;
use solidcore::geom::extent1d::Extent1d;
use solidcore::geom::extent2d::Extent2d;
use solidcore::intersect::surface_drop_lib::p1::DropFailRecord;
use solidcore::intersect::surface_drop_lib::p3::SM_ERR_DEGENERATE_SURFACE;
use solidcore::kernel::Kernel;
use solidcore::nurbs::b_spline_curve::BSplineCurve;
use solidcore::nurbs::b_spline_surface::BSplineSurface;
use solidcore::nurbs::crv_basic::Curve as NurbCurve;
use solidcore::nurbs::math::{CPoint, NOW};
use solidcore::nurbs::srf_basic::Surface as NurbSurface;
use solidcore::surface::surface::p2::{CreateUvTrimCurvePathType, DropCurveFail};
use solidcore::surface::surface::p4::NoObject;
use solidcore::surface::surface::Surface;

fn run(k: &Kernel, domain: Extent2d) -> (bool, DropCurveFail<NoObject>) {
    let net = vec![
        vec![
            CPoint::new(0.0, 0.0, 0.0, NOW),
            CPoint::new(0.0, 1.0, 0.0, NOW),
        ],
        vec![
            CPoint::new(1.0, 0.0, 0.0, NOW),
            CPoint::new(1.0, 1.0, 0.0, NOW),
        ],
    ];
    let nurb = NurbSurface::new(
        net,
        1,
        1,
        vec![
            domain.min().x,
            domain.min().x,
            domain.max().x,
            domain.max().x,
        ],
        vec![
            domain.min().y,
            domain.min().y,
            domain.max().y,
            domain.max().y,
        ],
    )
    .unwrap();
    // Borrowed NURBs preserve the caller's knot domain without knot fixup.
    let surface = Surface::from(BSplineSurface::from_nurb(nurb, true));
    let nurb = NurbCurve::new(
        vec![CPoint::new(0.5, 0.0, 0.5, NOW); 2],
        1,
        vec![0.0, 0.0, 1.0, 1.0],
    )
    .unwrap();
    let curve = Curve::from(BSplineCurve::from_nurb(nurb, 3, false));
    let mut fail = DropCurveFail::default();
    let result = surface.drop_curve_base(
        k,
        &domain,
        &curve,
        &Extent1d::from_min_max(0.0, 1.0),
        1e-4,
        true,
        Some(DropFailRecord {
            fail: &mut fail,
            surface: None,
            curve: None,
        }),
    );
    (result.is_err(), fail)
}

#[test]
fn opt_in_detects_v_degenerate_natural_domain() {
    let mut k = Kernel::default();
    k.fixes.drop_curve_v_domain_degeneracy = true;
    let (failed, fail) = run(&k, Extent2d::from_bounds(0.0, 0.0, 1.0, 0.0));
    assert!(failed);
    assert_eq!(fail.fail_status, SM_ERR_DEGENERATE_SURFACE);
    assert_eq!(
        fail.create_path_type,
        CreateUvTrimCurvePathType::DegenSurfaceFail
    );
    assert_eq!((fail.degen_surf_u, fail.degen_surf_v), (0, 1));
}

#[test]
fn u_degeneracy_retains_default_flags_and_opt_in_distinguishes_v() {
    let mut k = Kernel::default();
    assert!(!k.fixes.drop_curve_v_domain_degeneracy);
    let domain = Extent2d::from_bounds(0.0, 0.0, 0.0, 1.0);
    let (failed, fail) = run(&k, domain);
    assert!(failed);
    assert_eq!(fail.fail_status, SM_ERR_DEGENERATE_SURFACE);
    assert_eq!((fail.degen_surf_u, fail.degen_surf_v), (1, 1));
    k.fixes.drop_curve_v_domain_degeneracy = true;
    let (failed, fail) = run(&k, domain);
    assert!(failed);
    assert_eq!(fail.fail_status, SM_ERR_DEGENERATE_SURFACE);
    assert_eq!((fail.degen_surf_u, fail.degen_surf_v), (1, 0));
}

#[test]
fn default_preserves_upstream_v_degenerate_success() {
    let k = Kernel::default();
    assert!(!k.fixes.drop_curve_v_domain_degeneracy);
    let (failed, fail) = run(&k, Extent2d::from_bounds(0.0, 0.0, 1.0, 0.0));
    assert!(!failed);
    assert_ne!(
        fail.create_path_type,
        CreateUvTrimCurvePathType::DegenSurfaceFail
    );
    assert_eq!((fail.degen_surf_u, fail.degen_surf_v), (0, 0));
}
