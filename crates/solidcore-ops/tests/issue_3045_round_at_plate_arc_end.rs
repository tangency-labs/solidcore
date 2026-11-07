//! #3045: a round on the edge where a plate's straight side meets its arc end, from about 0.43 of
//! the arc's radius up. The round is the cylinder of radius R tangent to the side and inside the
//! arc's cylinder; it removes 8 times the area between the corner and that circle.
//!
//! The failing configuration needs three of the fixes together: the arc face a cylinder
//! (`linear_sweep_arc_is_cylinder`) and its offsets continued as themselves
//! (`offset_extension_from_base`, `offset_extension_as_itself`). The continued cylinder is a
//! whole circle, so the rail's end on the side edge has a second solution outside the face.

#![allow(clippy::disallowed_methods)] // closed-form oracle, not ops geometry
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)] // test

mod common;

use std::f64::consts::FRAC_PI_2;

use common::edge_at;
use solidcore_ops::blend::{blend, Blend, Corners};
use solidcore_ops::sweep::{extrude, Curve, DraftCorner};
use solidcore_ops::{session, Body, CheckLevel};

const CX: f64 = 93.75;
const CY: f64 = 15.0;
const RC: f64 = 16.25;
const THICKNESS: f64 = 8.0;

/// The plate: 100 x 30 with the +x end an arc of radius 16.25 about (93.75, 15), 8 thick.
fn plate() -> Body {
    let half = (CY / RC).asin();
    let line = |start: [f64; 3], end: [f64; 3]| Curve::Line { start, end };
    let profile = vec![
        line([0.0, 0.0, 0.0], [100.0, 0.0, 0.0]),
        Curve::Arc {
            center: [CX, CY, 0.0],
            axis: [0.0, 0.0, 1.0],
            reference: [1.0, 0.0, 0.0],
            radius: RC,
            start: -half,
            end: half,
        },
        line([100.0, 30.0, 0.0], [0.0, 30.0, 0.0]),
        line([0.0, 30.0, 0.0], [0.0, 0.0, 0.0]),
    ];
    extrude(
        &profile,
        [0.0, 0.0, 1.0],
        THICKNESS,
        0.0,
        DraftCorner::default(),
        true,
    )
    .unwrap()
}

/// The plate's volume: the 100 x 30 rectangle and the arc's segment past x = 100.
fn plate_volume() -> f64 {
    let turn = 2.0 * (CY / RC).asin();
    THICKNESS * (3000.0 + RC * RC / 2.0 * (turn - turn.sin()))
}

/// The area the round of radius `r` removes: between the corner (100, 0), the side and the arc,
/// and the round's circle (centre `(CX + sqrt((RC - r)^2 - (CY - r)^2), r)`, tangent to the side
/// and inside the arc's circle). Green's theorem round the loop corner, foot on the side,
/// contact on the arc: each circular arc contributes `1/2 (r^2 (b - a) + cx r (sin b - sin a) +
/// cy r (cos a - cos b))`, the side nothing.
fn removed_area(r: f64) -> f64 {
    let dx = ((RC - r).powi(2) - (CY - r).powi(2)).sqrt();
    let (fx, fy) = (CX + dx, r);
    let contact = (r - CY).atan2(dx);
    let corner = (-CY).atan2(100.0 - CX);
    let arc = |(cx, cy): (f64, f64), rad: f64, a: f64, b: f64| {
        0.5 * (rad * rad * (b - a)
            + cx * rad * (b.sin() - a.sin())
            + cy * rad * (a.cos() - b.cos()))
    };
    (arc((fx, fy), r, -FRAC_PI_2, contact) + arc((CX, CY), RC, contact, corner)).abs()
}

/// The round of radius `r` on the plate's edge at the corner: its volume against the closed form.
fn check_round(r: f64, flags: &str) {
    let b = plate();
    let e = edge_at(&b, [100.0, 0.0, THICKNESS / 2.0]);
    blend(&b, Blend::Round, &[e], &[r], &[], Corners::Filled)
        .unwrap_or_else(|err| panic!("{flags}, R {r}: {err:?}"));
    let report = b.check(CheckLevel::All).unwrap();
    assert!(report.is_valid(), "{flags}, R {r}: {report}");
    let volume = b.volume(1e-9).unwrap();
    let expected = plate_volume() - THICKNESS * removed_area(r);
    assert!(
        (volume - expected).abs() < 1e-6,
        "{flags}, R {r}: volume {volume}, closed form {expected}"
    );
}

/// From R 6.5 (0.40 of the arc's radius) up the round failed with status 1001 under the three
/// fixes; with none of them it always built. Both build to the closed form for every R.
#[test]
fn the_round_builds_for_every_radius() {
    const RADII: [f64; 9] = [3.0, 5.0, 6.0, 6.5, 7.0, 7.5, 8.0, 10.0, 12.0];
    // A single test: the fixes are one process-wide setting.
    for (flags, names) in [
        ("default flags", &[][..]),
        (
            "the three fixes",
            &[
                "linear_sweep_arc_is_cylinder",
                "offset_extension_from_base",
                "offset_extension_as_itself",
            ][..],
        ),
    ] {
        session::set_fixes(names.iter().copied()).unwrap();
        for r in RADII {
            check_round(r, flags);
        }
    }
    session::set_fixes([]).unwrap();
}
