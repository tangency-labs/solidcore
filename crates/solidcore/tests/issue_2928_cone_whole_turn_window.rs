//! #2928: `Cone::adjust_step_uv_domain` (`AdjustSTEPUVDomain`) on a whole-turn u window whose
//! ends sit a rounding error outside ±180° (the NIST `ctc_02` AP203 face #1152: u
//! −180.0000000000005 … 179.99999999999955°). The window is inside [−360°, 360°], but the
//! sweep circle was built a turn too high ([≈180°, ≈540°]) and `create_circle_segment`
//! refused it. A whole-turn window is accepted wherever its ends fall.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers
use solidcore::geom::{axis2_placement::Axis2Placement, extent2d::Extent2d};
use solidcore::{kernel::Kernel, math::Vec3, surface::cone::Cone};

const V_MAX: f64 = 39.31893628625494;

fn frame() -> Axis2Placement {
    Axis2Placement::new(
        Vec3::new(-2.4980018054066e-13, 820.0, -140.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(-1.0, 0.0, 0.0),
    )
}

#[test]
fn a_whole_turn_window_a_rounding_error_past_180_is_accepted() {
    let k = Kernel::default();
    for semi_angle in [3.0, 2.999999996572] {
        for (u0, u1) in [
            (-180.0000000000005, 179.99999999999955),
            (-180.0, 180.0),
            (0.0, 360.0),
            (-179.9999999999995, 180.0000000000005),
        ] {
            let mut cone = Cone::create_canonical_with_semi_angle(
                &k,
                &frame(),
                100.0,
                semi_angle,
                None,
                false,
                false,
            )
            .unwrap();
            let r = cone.adjust_step_uv_domain(&k, &Extent2d::from_bounds(u0, 0.0, u1, V_MAX));
            assert!(r.is_ok(), "semi-angle {semi_angle}, u {u0}..{u1}: {r:?}");
        }
    }
}
