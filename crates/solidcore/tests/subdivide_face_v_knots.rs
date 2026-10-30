//! `SmBrep::SubdivideFaceAtDiscontinuities` takes the U knot count for its V pass
//! (`SmBrep.cpp:18483`, WR-589). The fixtures are unit-square sheets on degree-1
//! B-spline surfaces (`fixtures/topology.brep_p10/subdivide_*.smb`), subdivided at
//! every interior knot (`SM_CT_CINFINITY`):
//! - one U knot and two V knots: upstream cuts at the U knot and the first V knot only
//!   (4 faces); `subdivide_face_counts_v_knots` cuts at both V knots (6 faces);
//! - no U knot and one V knot: upstream makes no cut (1 face); the fix cuts it (2);
//! - one U knot and no V knot: not compared with upstream ([`SolidError::InvalidInput`]
//!   here); the fix cuts at the U knot (2 faces).
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers

use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::geom::geom_utility::p2::ContinuityType;
use solidcore::geom::object::ObjectRef;
use solidcore::topology::brep::p1::faces;
use solidcore::topology::brep::p10::{subdivide_face_at_discontinuities, TopologyEnv};
use solidcore::topology::brep::p7::read_from_file;
use solidcore::topology::database_io::FileType;
use solidcore::topology::edge::NotifyEnv;
use solidcore::topology::topology::NotifyRecords;
use solidcore::topology::{EdgeId, FaceId};
use solidcore::{Brep, EntityRef, Kernel, SolidError};

/// A distinct identity per entity kind and index.
fn key_of(e: EntityRef) -> ObjectRef {
    let (kind, index) = match e {
        EntityRef::Brep(x) => (0, x.index()),
        EntityRef::Region(x) => (1, x.index()),
        EntityRef::Shell(x) => (2, x.index()),
        EntityRef::Faceuse(x) => (3, x.index()),
        EntityRef::Face(x) => (4, x.index()),
        EntityRef::Loopuse(x) => (5, x.index()),
        EntityRef::Loop(x) => (6, x.index()),
        EntityRef::Edgeuse(x) => (7, x.index()),
        EntityRef::Edge(x) => (8, x.index()),
        EntityRef::Vertexuse(x) => (9, x.index()),
        EntityRef::Vertex(x) => (10, x.index()),
    };
    ObjectRef {
        type_code: 0,
        key: ObjectKey((1u64 << 40) * (kind + 1) + u64::from(index)),
    }
}

/// Reads fixture `name`, subdivides its face at every interior knot; the result and the
/// number of faces afterwards.
fn subdivide(name: &str, counts_v: bool) -> (Result<Vec<FaceId>, SolidError>, usize) {
    let mut k = Kernel::default();
    k.fixes.subdivide_face_counts_v_knots = counts_v;
    let mut b = Brep::new_with_infinite_region().unwrap();
    let mut records = NotifyRecords::new();
    let mut attributes = AttributeStore::default();
    let mut env = TopologyEnv {
        all_attributes: &[],
        attributes: &mut attributes,
        key_of: &key_of,
        records: &mut records,
        next_key: 1 << 50,
    };
    read_from_file(
        &mut k,
        &mut b,
        fixture(name),
        FileType::Ascii,
        false,
        &mut env,
    )
    .unwrap();
    let face = faces(&b, None).unwrap()[0];
    let curve_ref = |e: EdgeId, _type_code: i64| key_of(EntityRef::Edge(e));
    let mut notify = NotifyEnv {
        records: &mut records,
        attributes: &mut attributes,
        context: None,
        key_of: &key_of,
    };
    let r = subdivide_face_at_discontinuities(
        &mut k,
        &mut b,
        &mut notify,
        &curve_ref,
        face,
        ContinuityType::CInfinity,
    );
    let n = faces(&b, None).unwrap().len();
    (r, n)
}

/// A crafted fixture of the `topology.brep_p10` difftest.
fn fixture(name: &str) -> Vec<u8> {
    let root =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/topology.brep_p10");
    std::fs::read(root.join(name)).unwrap()
}

#[test]
fn v_pass_takes_the_u_count_upstream() {
    let (r, n) = subdivide("subdivide_u1v2.smb", false);
    r.unwrap();
    assert_eq!(n, 4);
    let (r, n) = subdivide("subdivide_u0v1.smb", false);
    r.unwrap();
    assert_eq!(n, 1);
    let (r, _) = subdivide("subdivide_u1v0.smb", false);
    assert!(matches!(r, Err(SolidError::InvalidInput { .. })), "{r:?}");
}

#[test]
fn fix_cuts_at_every_v_knot() {
    let (r, n) = subdivide("subdivide_u1v2.smb", true);
    r.unwrap();
    assert_eq!(n, 6);
    let (r, n) = subdivide("subdivide_u0v1.smb", true);
    r.unwrap();
    assert_eq!(n, 2);
    let (r, n) = subdivide("subdivide_u1v0.smb", true);
    r.unwrap();
    assert_eq!(n, 2);
}
