//! tangency-labs/solidcore#391 (WR-80): `GlobalSurfaceIntersect` fails with
//! `SM_ERR_NULL_POINTER` when a boundary start point lies on an intersection curve that
//! was already found. `SmSurfaceIntersector::TraceCurve` then returns `SM_SUCCESS` with an
//! empty point list (and leaves `m_bCurveIsClosed` at the previous trace's FALSE), and
//! `FindInteriorCurves`/`DoIntersection` call `ReverseCurveDirection` on the empty list.
//!
//! The surfaces are two single-span Bezier patches (degrees 2 and 4), built here from
//! their formulas, that touch to the 4th order along the parabola `y = 0`, `z = 0.2 x^2`
//! (`x = 2s - 1`, `y = t - 1/2`, `s` and `t` in `[0, 1]`):
//!
//! - `A`: `z = 0.5 y^3 - 0.5 y + 0.2 x^2 + y^4`
//! - `B`: `z = 0.5 y^3 - 0.5 y + 0.2 x^2`
//!
//! Their difference is `y^4`, so the boundary start points on `x = +-1` come in clusters
//! about 1e-3 apart and a later point of a cluster lies on the curve its neighbour traced.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers
use solidcore::error::SM_ERR_NULL_POINTER;
use solidcore::geom::tol::approx_tol3d;
use solidcore::geom::vector2d::Point2;
use solidcore::kernel::Kernel;
use solidcore::nurbs::b_spline_surface::BSplineSurface;
use solidcore::nurbs::math::CPoint;
use solidcore::nurbs::srf_basic::Surface as NurbSurface;
use solidcore::surface::surface::p2::SurfaceIntersection;
use solidcore::surface::surface::Surface;
use solidcore::SolidError;

/// `n choose k` as a float.
fn binom(n: usize, k: usize) -> f64 {
    (0..k).fold(1.0, |acc, i| acc * (n - i) as f64 / (i + 1) as f64)
}

/// Bernstein coefficients of the polynomial with ascending monomial coefficients `mono`.
fn bernstein(mono: &[f64]) -> Vec<f64> {
    let n = mono.len() - 1;
    (0..=n)
        .map(|j| {
            (0..=j)
                .map(|k| binom(j, k) / binom(n, k) * mono[k])
                .sum::<f64>()
        })
        .collect()
}

/// `z = 0.5 y^3 - 0.5 y + 0.2 x^2` (+ `y^4` when `with_quartic`) as a degree (2, 4) patch.
fn patch(with_quartic: bool) -> Surface {
    // 0.2 (2s - 1)^2 in s
    let g = bernstein(&[0.2, -0.8, 0.8]);
    // 0.5 w^3 - 0.5 w (+ w^4) with w = t - 1/2, expanded in t
    let f = if with_quartic {
        bernstein(&[0.25, -0.625, 0.75, -1.5, 1.0])
    } else {
        bernstein(&[0.1875, -0.125, -0.75, 0.5, 0.0])
    };
    let xs = [-1.0, 0.0, 1.0]; // x = 2s - 1, degree-elevated to 2
    let ys = [-0.5, -0.25, 0.0, 0.25, 0.5]; // y = t - 1/2, degree-elevated to 4
    let net = (0..3)
        .map(|i| {
            (0..5)
                .map(|j| CPoint::new(xs[i], ys[j], g[i] + f[j], 1.0))
                .collect()
        })
        .collect();
    let nurb = NurbSurface::new(
        net,
        2,
        4,
        vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    )
    .unwrap();
    Surface::from(BSplineSurface::from_nurb(nurb, false))
}

#[test]
fn quartic_contact_patches_have_the_stated_shape() {
    let k = Kernel::default();
    for (with_quartic, a) in [(true, patch(true)), (false, patch(false))] {
        for s in [0.0, 0.25, 0.5, 0.8, 1.0] {
            for t in [0.0, 0.3, 0.5, 0.9, 1.0] {
                let p = a.evaluate_point(&k, Point2::new(s, t)).unwrap();
                let (x, w) = (2.0 * s - 1.0, t - 0.5);
                let mut z = 0.5 * w * w * w - 0.5 * w + 0.2 * x * x;
                if with_quartic {
                    z += w * w * w * w;
                }
                assert!((p.x - x).abs() < 1e-12, "x at ({s}, {t})");
                assert!((p.y - w).abs() < 1e-12, "y at ({s}, {t})");
                assert!((p.z - z).abs() < 1e-12, "z at ({s}, {t}): {} vs {z}", p.z);
            }
        }
    }
}

/// `GlobalSurfaceIntersect` of the two patches at the corpus' approximation tolerance.
fn intersect(k: &Kernel) -> solidcore::Result<SurfaceIntersection> {
    let (a, b) = (patch(true), patch(false));
    let (da, db) = (a.natural_uv_domain(), b.natural_uv_domain());
    a.global_surface_intersect(
        k,
        &da,
        &b,
        &db,
        [true, true],
        Some(approx_tol3d(k, 1e-6)),
        None,
    )
}

/// With `Kernel::fixes.ssi_empty_trace_not_reversed` the start points that lie on a found
/// curve add nothing and the intersection goes on. This fix alone does not find the
/// tangent curve (that needs the cluster fixes of #392/#394), so only the failure is
/// asserted gone.
#[test]
fn issue_391_empty_trace_is_not_reversed_with_the_fix() {
    let mut k = Kernel::default();
    k.fixes.ssi_empty_trace_not_reversed = true;
    let x = intersect(&k).unwrap();
    // whatever comes back is consistent: one entry per curve in every array
    let n = x.curves_3d.len();
    assert_eq!(x.surface1_uv_curves.len(), n);
    assert_eq!(x.surface2_uv_curves.len(), n);
    assert_eq!(x.curve_types.len(), n);
    assert_eq!(x.deviations.len(), n);
}

/// As SMLib (the default): the reversal of the empty trace is `SM_ERR_NULL_POINTER`.
/// The operation is checked too, so another `SM_ERR_NULL_POINTER` of the intersector
/// cannot keep this test green while the fixed path becomes unreachable.
#[test]
fn issue_391_empty_trace_fails_by_default() {
    let k = Kernel::default();
    let r = intersect(&k);
    match r {
        Err(SolidError::Upstream {
            code, operation, ..
        }) => {
            assert_eq!(code, SM_ERR_NULL_POINTER);
            assert_eq!(
                operation,
                "intersect::surface_intersector::p1::SurfaceIntersector::reverse_curve_direction"
            );
        }
        other => panic!("expected SM_ERR_NULL_POINTER, got {:?}", other.map(|_| ())),
    }
}
