//! Issue #2954 / review follow-up: curved-seam shell-tangent handling.
//!
//! Verifies that when shelling with open faces adjacent to tangent curved seams,
//! rails derived from boundary curves track the curved boundary within certified tolerance.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers

#[path = "support/offset_fixture.rs"]
mod offset_fixture;

use offset_fixture::{planar_faces_at, shell};
use solidcore::primitive::api_primitives::create_cylinder;
use solidcore::query::api_queries::brep_is_manifold_solid;
use solidcore::topology::brep::p1;
use solidcore::{Kernel, Point3};

#[test]
fn test_curved_seam_shell_tangent_rail_to_boundary() -> solidcore::Result<()> {
    let mut k = Kernel::default();

    // A cylinder r=10, h=20 shelled inward 1.0 with top cap open.
    // The cylinder's rim is a circular curved seam adjacent to the open face and the cylindrical wall.
    let cyl = create_cylinder(&mut k, Point3::new(0.0, 0.0, 0.0), 10.0, 20.0)?;
    let top = planar_faces_at(&k, &cyl, 2, 20.0);
    assert_eq!(top.len(), 1, "top planar face");

    let shelled = shell(&mut k, cyl, 1.0, &top)?;
    assert!(
        brep_is_manifold_solid(&shelled)?,
        "shelled cylinder must be a manifold solid"
    );

    // Verify all curved edges on the shelled body:
    // The rim and bottom curves must lie on their respective analytic boundaries (radius 10.0 or 9.0)
    // within the edge's certified tolerance.
    let mut checked_curved_edges = 0;
    for e in p1::edges(&shelled, None)? {
        let edge = shelled.edge(e)?;
        if let Some(ref curve) = edge.curve {
            let ivl = edge.interval;
            let tol = edge.tolerance.max(1e-5);
            let p_start = curve.evaluate_point(&k, ivl.min())?;
            let p_mid = curve.evaluate_point(&k, ivl.evaluate(0.5))?;
            let p_end = curve.evaluate_point(&k, ivl.max())?;

            let r_start = (p_start.x * p_start.x + p_start.y * p_start.y).sqrt();
            let r_mid = (p_mid.x * p_mid.x + p_mid.y * p_mid.y).sqrt();
            let r_end = (p_end.x * p_end.x + p_end.y * p_end.y).sqrt();

            if (r_start - 10.0).abs() < 1e-3 || (r_start - 9.0).abs() < 1e-3 {
                checked_curved_edges += 1;
                let expected_r = if (r_start - 10.0).abs() < 1e-3 {
                    10.0
                } else {
                    9.0
                };
                assert!(
                    (r_start - expected_r).abs() <= tol,
                    "curve start radius {r_start} must match boundary {expected_r} within {tol}"
                );
                assert!(
                    (r_mid - expected_r).abs() <= tol,
                    "curve midpoint radius {r_mid} must match boundary {expected_r} within {tol}"
                );
                assert!(
                    (r_end - expected_r).abs() <= tol,
                    "curve end radius {r_end} must match boundary {expected_r} within {tol}"
                );
            }
        }
    }

    assert!(
        checked_curved_edges >= 2,
        "must have verified at least 2 curved boundary edges (found {checked_curved_edges})"
    );
    Ok(())
}
