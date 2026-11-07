//! R10 on a 5 x 5 bar consumes both original supports: its new seams cross, not touch.
//!
//! The section is the square intersected with the circle centred at (y,z)=(10,-5),
//! radius 10. Its two surviving long cylinder/plane seams have outward-normal dot
//! product -1/2 (120 degrees), not the nonexistent original supports' G1 contact.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::disallowed_methods
)] // Independent closed-form oracle and corpus test data, not ops geometry.

use std::io::Read;
use std::path::{Path, PathBuf};

use serde_json::Value;
use solidcore_ops::blend::{blend, Blend, Corners};
use solidcore_ops::{session, Body, CheckLevel};
use solidcore_record::brep::Brep;

const CASE: &str = "4782-fillet-bar-edge-wider-than-both-faces";

fn corpus() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../solidcore/tests/kernel_cases")
}

fn entry() -> Value {
    serde_json::from_str(
        &std::fs::read_to_string(corpus().join(format!("fillets/{CASE}.json"))).unwrap(),
    )
    .unwrap()
}

#[test]
fn the_case_grades_the_existing_overflow_not_absent_support_seams() {
    let entry = entry();
    assert!(
        entry["expected"].get("continuity").is_none(),
        "the original y=0 and z=5 supports are consumed; their G1 seams do not exist"
    );
    assert_eq!(entry["expected"]["volume_mm3"], 630.2934873);
    assert_eq!(entry["expected"]["tolerance"]["rel"], 1e-6);
    let known: Value = serde_json::from_str(
        &std::fs::read_to_string(corpus().join("known_failing.json")).unwrap(),
    )
    .unwrap();
    assert!(!known
        .as_array()
        .unwrap()
        .iter()
        .any(|row| { row["case"].as_str() == Some(&format!("fillets/{CASE}")) }));
}

#[test]
fn the_oversized_round_keeps_its_closed_form_and_two_crossing_seams() {
    let root = corpus();
    let fixes: Vec<String> =
        serde_json::from_str(&std::fs::read_to_string(root.join("fixes.json")).unwrap()).unwrap();
    session::set_fixes(fixes.iter().map(String::as_str)).unwrap();
    let entry = entry();
    let path = root
        .join("fillets")
        .join(entry["records"].as_str().unwrap());
    let mut text = String::new();
    flate2::read::GzDecoder::new(std::fs::File::open(path).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    let records: Value = serde_json::from_str(&text).unwrap();
    let record: Brep = serde_json::from_value(records["b2"].clone()).unwrap();
    let (input, order) = Body::from_record(&record, true).unwrap();
    let before = input.to_record().unwrap().0;
    let (body, _) = input.copy_with_order(&order).unwrap();
    blend(&body, Blend::Round, &[1], &[10.0], &[], Corners::Filled).unwrap();

    assert_eq!(body.counts().unwrap(), [5, 9, 6]);
    assert!(body.is_manifold_solid().unwrap());
    let check = body.check(CheckLevel::All).unwrap();
    assert!(check.is_valid(), "{check}");
    let volume = body.volume(1e-9).unwrap();
    let exact = 2000.0 * std::f64::consts::PI / 3.0 - 400.0 * 75.0_f64.sqrt() + 2000.0;
    assert!(
        (volume - exact).abs() <= 1e-9 * exact,
        "{volume} != {exact}"
    );

    let (result, order) = body.to_record().unwrap();
    let mut crossings = 0;
    let mut sharp_planar = 0;
    let mut long_edges = 0;
    for (i, edge) in result.edges.iter().enumerate() {
        let a = result.vertices[edge.vertices[0] as usize].point;
        let b = result.vertices[edge.vertices[1] as usize].point;
        if ((a[0] - b[0]).abs() - 80.0).abs() > 1e-9 {
            continue;
        }
        long_edges += 1;
        assert_eq!(edge.uses.len(), 2);
        let faces: Vec<_> = edge.uses.iter().map(|u| u.face as usize).collect();
        let normals: Vec<_> = faces
            .iter()
            .map(|&f| body.face_normal_at_edge(&order, i, f).unwrap())
            .collect();
        let dot: f64 = normals[0]
            .value
            .iter()
            .zip(normals[1].value)
            .map(|(x, y)| x * y)
            .sum();
        if (dot + 0.5).abs() <= 1e-9 {
            crossings += 1;
            let yz = [a[1], a[2]];
            let lo = [10.0 - 75.0_f64.sqrt(), 0.0];
            let hi = [5.0, 75.0_f64.sqrt() - 5.0];
            assert!([lo, hi]
                .iter()
                .any(|want| { yz.iter().zip(want).all(|(x, y)| (x - y).abs() <= 1e-9) }));
        } else {
            assert!(dot.abs() <= 1e-9, "unexpected long seam normal dot: {dot}");
            sharp_planar += 1;
        }
        for &face in &faces {
            let samples = body.edge_face_gap_samples(&order, i, face).unwrap();
            assert!(!samples.is_empty());
            for sample in samples {
                assert!(sample.distance.value <= edge.tolerance + sample.distance.bound);
            }
        }
    }
    assert_eq!((long_edges, crossings, sharp_planar), (3, 2, 1));
    assert_eq!(input.to_record().unwrap().0, before, "operand changed");
    eprintln!("#2948 volume={volume:.17}; exact={exact:.17}; topology=5/9/6; long seams=3, 120-degree crossings=2");
}
