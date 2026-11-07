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

use common::{edge_at, kernel_copy, volume};
use solidcore::surface::surface::SurfaceKind;
use solidcore::topology::face::p2 as face;
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
    // the fillets of both bodies are cones, not surfaces of extrusion
    session::set_fixes(["linear_sweep_arc_is_cylinder", "classify_edge_interval"]).unwrap();
    for dx in [1e-4, 1e-3, 1e-2] {
        for reversed in [false, true] {
            let (a, b) = (rounded_box(), extruded(dx));
            let (x, y) = if reversed { (b, a) } else { (a, b) };
            let r = boolean(x, y, BooleanOp::Subtract)
                .unwrap_or_else(|e| panic!("dx {dx} reversed {reversed}: {e}"))
                .unwrap();
            let report = r.check(CheckLevel::All).unwrap();
            assert!(report.is_valid(), "dx {dx} reversed {reversed}: {report}");
            assert_eq!(r.solid_count().unwrap(), 1, "dx {dx} reversed {reversed}");
            // Every face of the sliver: the two flat sides (26 x 10 each) and the four fillets
            // (a quarter turn of R2 by 10 each), plus the thin ends. A fillet face the sector
            // skipped and a neighbouring edge put in the wrong region would add or drop 31.4.
            let faces = 2.0 * 26.0 * 10.0 + 4.0 * FRAC_PI_2 * 2.0 * 10.0;
            let area = r.area(1e-8).unwrap();
            assert!(
                (area - faces).abs() < 100.0 * dx,
                "dx {dx} reversed {reversed}: area {area} vs {faces}"
            );
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

/// The straight-edge fillets of both operands are cones of semi-angle 0 (`SurfaceKind::Cone`,
/// `is_cylinder`), the kind `faces_coincide` gates: a cylinder has no kind of its own.
#[test]
fn issue_3012_the_fillets_are_cylinders_of_the_cone_kind() {
    session::set_fixes(["linear_sweep_arc_is_cylinder", "classify_edge_interval"]).unwrap();
    for body in [rounded_box(), extruded(1e-4)] {
        let (_, brep) = kernel_copy(&body);
        let round: Vec<_> = brep
            .face_list()
            .iter()
            .filter_map(
                |&f| match face::surface(&brep, f).unwrap().unwrap().kind() {
                    SurfaceKind::Plane(_) => None,
                    SurfaceKind::Cone(c) => Some(c.is_cylinder()),
                    other => panic!("a fillet face of kind {other:?}"),
                },
            )
            .collect();
        assert_eq!(round, [true; 4]);
    }
}
