//! tangency-labs/solidcore#2539: SMLib's check refused exact results where a plane nearly
//! holds a cylinder's axis.
//!
//! SolidCraft's enclosure part cuts its case body `periph_case_R`, and in #2539's reduction
//! `e3` a cone, with `t24`: a box clipped by a slanted half-space and notched by a cylinder of
//! radius 8 whose axis lies along the clipped box's slanted face, to 4.5e-7 rad. The face meets
//! the notch in an ellipse of radii 1.8e7 and 8 over 2.3e-4 degrees, a line to 6e-11 over its
//! length. Its analytic form, `centre + rx cos(t) X + ry sin(t) Y`, carries rounding of about
//! `rx` times an ulp of `t` (4e-9 here). `SmEllipse`'s tests 1 and 2 hold the analytic form to
//! its NURB within `100 SM_EFF_ZERO (1 + |p|)`, 1.7e-9 to 3.8e-9 at the part's size, so the
//! results, exact to rounding, were refused. The tests now allow the rounding the evaluation
//! carries, `gamma_11 max(|rx|, |ry|) (1 + |t|)` (WR-741).
//!
//! The bodies are SolidCraft's records, not rebuilt here: upstream's verdict on this ellipse
//! turns on its last bits, so the same operands built in a test can pass or fail it. On main
//! each test fails on `SmEllipse` alone.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers

use solidcore::boolean::api_brep::boolean;
use solidcore::geom::a_object::BooleanOperation;
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::kernel::UpstreamFixes;
use solidcore::mass::integrator::properties::{brep_volume, PropertyAttributes};
use solidcore::topology::api_files::read_brep_from_file;
use solidcore::topology::brep::p10::TopologyEnv;
use solidcore::topology::brep::p2::is_manifold_solid;
use solidcore::topology::brep::validation::{check_brep, CheckLevel};
use solidcore::topology::loopuse;
use solidcore::topology::topology::NotifyRecords;
use solidcore::{Brep, Kernel};

/// SolidCraft's profile, every fix on but `fillet_refuse_rounds_breaking_through` (its
/// `HELD_BACK`).
fn solidcraft() -> Kernel {
    let mut k = Kernel {
        fixes: UpstreamFixes::all(),
        ..Kernel::default()
    };
    k.fixes.fillet_refuse_rounds_breaking_through = false;
    k
}

/// `b` checks clean at Level 2.
fn checks_clean(k: &mut Kernel, b: &mut Brep) {
    let lines = check_brep(k, b, CheckLevel::Level2).unwrap().lines();
    assert!(lines.is_empty(), "{lines:?}");
}

/// A body as SolidCraft records it (binary SMB), read under its naming.
fn read_notify_named(k: &mut Kernel, bytes: &[u8]) -> Brep {
    let mut attributes = AttributeStore::default();
    let mut records = NotifyRecords::new();
    let mut env = TopologyEnv {
        all_attributes: &[],
        attributes: &mut attributes,
        key_of: &loopuse::notify_ref,
        records: &mut records,
        next_key: 1 << 50,
    };
    let mut b = read_brep_from_file(k, Some(bytes.to_vec()), false, false, &mut env).unwrap();
    b.refresh_cached_topology();
    b
}

/// `t24` and `e3`'s cone as SolidCraft main 4bdffef06 (on solidcore 3df8a06b5) passes them to
/// the cut `e3`, written where its check refused the cut. `t24` is SolidCraft's 198 x 105.5 x
/// 71.5 box minus the half-space `t25` and minus the notch `t26`, a rectangle 8 wide turned
/// about its side on the clipped box's slanted face. The cone is the trapezoid (0, -55),
/// (24, -55), (19, 16.5), (0, 16.5) turned about world Z through (0, -5.8).
const TOOL: &[u8] = include_bytes!("data/issue_2539_e3_tool.smb");
const CONE: &[u8] = include_bytes!("data/issue_2539_e3_cone.smb");
/// `periph_case_R` before its last cut (SolidCraft's record, as above) cut down by solidcore to
/// its part below z = -14, where `t24` takes its corner: 18 faces, written by
/// `write_brep_to_file`. The whole part's record is 338 KB.
const PART: &[u8] = include_bytes!("data/issue_2539_periph_case_lower.smb");

fn volume(k: &Kernel, b: &mut Brep) -> f64 {
    let store = AttributeStore::default();
    let attrs = PropertyAttributes {
        store: &store,
        key_of: &|_| ObjectKey(0),
    };
    brep_volume(k, b, 1e-9, &attrs).unwrap()
}

/// `t24` itself: its slanted face meets the notch in the near-line ellipse, which checks
/// clean (upstream's tolerance refused it: `SmEllipse` test 2, 3.97e-9 against 3.81e-9).
#[test]
fn issue_2539_the_notched_tool_checks_clean() {
    let mut k = solidcraft();
    let mut tool = read_notify_named(&mut k, TOOL);
    assert!(is_manifold_solid(&tool).unwrap());
    checks_clean(&mut k, &mut tool);
}

/// `e3`: the cone minus `t24` is one closed solid of 8 faces, the cone's 104300.353 less
/// what the tool takes, and checks clean (upstream's tolerance refused it on the ellipse it
/// keeps from `t24`, as above).
#[test]
fn issue_2539_e3_cone_minus_the_notched_tool() {
    let mut k = solidcraft();
    let cone = read_notify_named(&mut k, CONE);
    let tool = read_notify_named(&mut k, TOOL);
    let mut r = boolean(&mut k, cone, tool, BooleanOperation::Difference)
        .unwrap()
        .result;
    assert!(is_manifold_solid(&r).unwrap());
    assert_eq!(r.face_list().len(), 8);
    let v = volume(&k, &mut r);
    assert!((v - 100_808.551_793).abs() < 1e-5, "{v}");
    checks_clean(&mut k, &mut r);
}

/// `periph_case_R`'s lower part minus `t24` is one closed solid that loses exactly what it
/// shares with `t24`, and checks clean (upstream's tolerance refused it on the ellipse the cut
/// trims, `SmEllipse` tests 1 and 2 at 4.64e-9 against 1.7e-9, as it refused the whole part's
/// cut, 72170.402 to 70319.173).
#[test]
fn issue_2539_periph_case_minus_the_notched_tool() {
    let mut k = solidcraft();
    let (part, tool) = (
        read_notify_named(&mut k, PART),
        read_notify_named(&mut k, TOOL),
    );
    let mut common = boolean(&mut k, part, tool, BooleanOperation::Intersection)
        .unwrap()
        .result;
    let shared = volume(&k, &mut common);
    let mut part = read_notify_named(&mut k, PART);
    let before = volume(&k, &mut part);
    let tool = read_notify_named(&mut k, TOOL);
    let mut r = boolean(&mut k, part, tool, BooleanOperation::Difference)
        .unwrap()
        .result;
    assert!(is_manifold_solid(&r).unwrap());
    let after = volume(&k, &mut r);
    // the cut takes `t24`'s corner of the slab, 1532.499 measured, whatever the intersection
    // says: two disjoint classifications cannot pass
    assert!(before - after > 1500.0, "{before} - {after}");
    assert!(
        (before - after - shared).abs() < 1e-4,
        "{before} - {after} against {shared}"
    );
    checks_clean(&mut k, &mut r);
}
