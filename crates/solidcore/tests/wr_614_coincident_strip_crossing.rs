//! tangency-labs/solidcore#1794 (WR-614): two 4 mm cubes, B turned atan(0.1) (5.71 deg)
//! about z, whose faces y = 20 (A) and F1 (B's y' = 0) cross at x = 1.9999 and stay
//! within 1e-5 of each other across the strip x in [1.9998, 2]: B's corner edge lies
//! 1e-5 below A's face, A's corner edge 1e-5 off F1 (the default tolerances; at 5e-6 the
//! union happens to succeed). The crossing is fixed only to
//! tol / sin, the whole strip, yet becomes an edge whose end vertices are picked per
//! body and never related, and the merge fails (`SmMerge::MakeFaceuse`, `SM_ERR`). The
//! opt-in `Kernel::fixes.boolean_coincident_strip_drops_crossing` drops that crossing: the
//! union then succeeds with the topology of the union a hair apart (d = 5e-6), which
//! upstream already makes. (That result keeps every face use in the infinite region's
//! shell, with no solid region, in both engines: a separate matter.)
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers
use solidcore::boolean::api_brep::boolean_union;
use solidcore::error::SM_ERR;
use solidcore::math::{Point3, Vec3};
use solidcore::primitive::api_primitives::create_box;
use solidcore::topology::brep::p1;
use solidcore::topology::proto_topology::p1::api_transform::{
    api_rotate, api_translate, TransformObject,
};
use solidcore::{Brep, Kernel, SolidError};

const SIDE: f64 = 4.0;
const SLOPE: f64 = 0.1;
/// How far each bounding edge lies off the other cube's face.
const D: f64 = 1e-5;

/// A: [-2, 2] x [16, 20] x [0, 4]. B: a 4 mm cube turned about z so that its face F1
/// crosses y = 20 at x = 2 - d / SLOPE, with its corner edge at x = 2 - 2 d / SLOPE,
/// d below y = 20.
fn operands(k: &mut Kernel, d: f64) -> (Brep, Brep) {
    let a = create_box(k, Point3::new(-2.0, 16.0, 0.0), SIDE, SIDE, SIDE).unwrap();
    let mut b = create_box(k, Point3::new(0.0, 0.0, 0.0), SIDE, SIDE, SIDE).unwrap();
    let o = Point3::new(0.0, 0.0, 0.0);
    let deg = SLOPE.atan().to_degrees();
    api_rotate(
        k,
        Some(TransformObject::Brep(&mut b)),
        o,
        Vec3::new(0.0, 0.0, 1.0),
        deg,
    )
    .unwrap();
    let half = d / SLOPE;
    let shift = Vec3::new((2.0 - half) - half, 20.0 - d, 0.0);
    api_translate(k, Some(TransformObject::Brep(&mut b)), shift).unwrap();
    (a, b)
}

/// Face, edge and vertex counts and the regions' void flags.
fn shape(b: &Brep) -> (usize, usize, usize, Vec<i32>) {
    let voids = p1::regions(b, None)
        .unwrap()
        .into_iter()
        .map(|r| b.region(r).unwrap().is_void)
        .collect();
    (
        p1::faces(b, None).unwrap().len(),
        p1::edges(b, None).unwrap().len(),
        p1::vertices(b, None).unwrap().len(),
        voids,
    )
}

#[test]
fn wr_614_coincident_strip_drops_crossing_is_off_by_default() {
    assert!(
        !Kernel::default()
            .fixes
            .boolean_coincident_strip_drops_crossing
    );
}

#[test]
fn wr_614_union_fails_upstream() {
    let mut k = Kernel::default();
    let (a, b) = operands(&mut k, D);
    match boolean_union(&mut k, a, b) {
        Err(SolidError::Upstream { code, .. }) => assert_eq!(code, SM_ERR),
        other => panic!("upstream's union fails with SM_ERR, got {other:?}"),
    }
}

#[test]
fn wr_614_union_with_the_fix() {
    // upstream's union of the same cubes 5e-6 apart, which it makes
    let mut k = Kernel::default();
    let (a, b) = operands(&mut k, 5e-6);
    let want = shape(&boolean_union(&mut k, a, b).unwrap());
    assert_eq!(want.0, 12);

    let mut k = Kernel::default();
    k.fixes.boolean_coincident_strip_drops_crossing = true;
    let (a, b) = operands(&mut k, D);
    let u = boolean_union(&mut k, a, b).unwrap();
    assert_eq!(shape(&u), want);
}

#[test]
fn wr_614_fix_leaves_the_union_apart_alone() {
    let mut k = Kernel::default();
    let (a, b) = operands(&mut k, 5e-6);
    let want = shape(&boolean_union(&mut k, a, b).unwrap());
    let mut k = Kernel::default();
    k.fixes.boolean_coincident_strip_drops_crossing = true;
    let (a, b) = operands(&mut k, 5e-6);
    assert_eq!(shape(&boolean_union(&mut k, a, b).unwrap()), want);
}
