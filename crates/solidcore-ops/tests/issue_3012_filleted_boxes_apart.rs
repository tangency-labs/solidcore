//! #3012: a 40 x 30 x 10 box with its four vertical edges rounded R2, less the same shape
//! (an extruded rounded rectangle) shifted along x by dx, failed in both orders: the
//! boolean read the two fillets as coincident faces at their mated tangent edge, though
//! they part by up to dx along the surface, and the region of a face found coincident
//! there was the wrong side's. The result is one body of the sliver, dx * 30 * 10 mm3
//! less its two tapered ends. Its own test binary: the fix set is the process kernel's.

#![allow(
    clippy::unwrap_used,
    clippy::indexing_slicing,
    clippy::disallowed_methods
)]

mod common;

use common::{edge_at, volume};
use solidcore_ops::boolean::boolean;
use solidcore_ops::sweep::{extrude, Curve, DraftCorner};
use solidcore_ops::{blend, make_box, session, Blend, Body, BooleanOp, CheckLevel, Corners};
use std::f64::consts::{FRAC_PI_2, PI};

fn rounded_box() -> Body {
    let b = make_box([0.0; 3], [40.0, 30.0, 10.0]).unwrap();
    let edges: Vec<usize> = [[0.0, 0.0], [40.0, 0.0], [0.0, 30.0], [40.0, 30.0]]
        .iter()
        .map(|p| edge_at(&b, [p[0], p[1], 5.0]))
        .collect();
    blend(&b, Blend::Round, &edges, &[2.0], &[], Corners::Filled).unwrap();
    b
}

/// The same outline as a rounded rectangle, shifted by `dx`, extruded.
fn extruded(dx: f64) -> Body {
    let (x0, x1, y1, r) = (dx, 40.0 + dx, 30.0, 2.0);
    let line = |s: [f64; 2], e: [f64; 2]| Curve::Line {
        start: [s[0], s[1], 0.0],
        end: [e[0], e[1], 0.0],
    };
    let arc = |c: [f64; 2], start: f64| Curve::Arc {
        center: [c[0], c[1], 0.0],
        axis: [0.0, 0.0, 1.0],
        reference: [1.0, 0.0, 0.0],
        radius: r,
        start,
        end: start + FRAC_PI_2,
    };
    let outline = [
        line([x0 + r, 0.0], [x1 - r, 0.0]),
        arc([x1 - r, r], -FRAC_PI_2),
        line([x1, r], [x1, y1 - r]),
        arc([x1 - r, y1 - r], 0.0),
        line([x1 - r, y1], [x0 + r, y1]),
        arc([x0 + r, y1 - r], FRAC_PI_2),
        line([x0, y1 - r], [x0, r]),
        arc([x0 + r, r], PI),
    ];
    extrude(
        &outline,
        [0.0, 0.0, 1.0],
        10.0,
        0.0,
        DraftCorner::Curvature,
        true,
    )
    .unwrap()
}

#[test]
fn issue_3012_cut_of_filleted_boxes_a_sliver_apart() {
    session::set_fixes(solidcore::kernel::UpstreamFixes::NAMES.iter().copied()).unwrap();
    for dx in [1e-4, 1e-3] {
        for reversed in [false, true] {
            let (a, b) = (rounded_box(), extruded(dx));
            let (x, y) = if reversed { (b, a) } else { (a, b) };
            let r = boolean(x, y, BooleanOp::Subtract)
                .unwrap_or_else(|e| panic!("dx {dx} reversed {reversed}: {e}"))
                .unwrap();
            assert!(r.check(CheckLevel::All).unwrap().is_valid(), "dx {dx}");
            // the sliver's volume dx * 300, less the taper at its ends: 1e-3 of it at most
            let expected = dx * 300.0;
            let v = volume(&r);
            assert!(
                (v - expected).abs() < 1e-3 * expected,
                "dx {dx} reversed {reversed}: {v} vs {expected}"
            );
        }
    }
}
