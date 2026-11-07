//! Projected face classification must preserve certificate refusals.
use solidcore::curve::curve::Curve;
use solidcore::geom::{
    attribute::{AttributeStore, ObjectKey},
    object::ObjectRef,
};
use solidcore::nurbs::{b_spline_curve::BSplineCurve, crv_basic::Curve as Net, math::CPoint};
use solidcore::topology::brep::p10::{
    make_topology_from_data, project_and_trim, MakeTopologyOptions, ProjectAndTrimOptions,
    ProjectEnv, TopologyEnv, TrimType,
};
use solidcore::topology::{brep_data::BrepData, database_io::FileType};
use solidcore::{math::Vec3, Kernel};
use solidcore::{
    topology::{brep_cache::CacheObject, edge::NotifyEnv, topology::NotifyRecords, Brep, BrepId},
    EntityRef,
};
fn key_of(e: EntityRef) -> ObjectRef {
    let (t, k, i) = match e {
        EntityRef::Brep(x) => (16011, 0, x.index()),
        EntityRef::Region(x) => (16012, 1, x.index()),
        EntityRef::Shell(x) => (16013, 2, x.index()),
        EntityRef::Faceuse(x) => (16022, 3, x.index()),
        EntityRef::Face(x) => (16002, 4, x.index()),
        EntityRef::Loopuse(x) => (16004, 5, x.index()),
        EntityRef::Loop(x) => (16003, 6, x.index()),
        EntityRef::Edgeuse(x) => (16005, 7, x.index()),
        EntityRef::Edge(x) => (16016, 8, x.index()),
        EntityRef::Vertexuse(x) => (16006, 9, x.index()),
        EntityRef::Vertex(x) => (16021, 10, x.index()),
    };
    ObjectRef {
        type_code: t,
        key: ObjectKey((1u64 << 40) * (k + 1) + u64::from(i)),
    }
}
#[test]
fn projected_side_refusals_do_not_become_successful_face_deletions() {
    let base =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/topology.brep_p10");
    for file in ["raw_vote_rectangle.smb", "raw_vote_midpoint_lens.smb"] {
        for gap in [1.5e-12, 1.9e-12, 1e-8] {
            for trim in [0, 1] {
                let mut k = Kernel::default();
                let mut data = BrepData::read_from_file(
                    &k,
                    std::fs::read(base.join(file)).unwrap(),
                    FileType::Ascii,
                )
                .unwrap();
                let mut b = Brep::new_with_infinite_region().unwrap();
                let mut store = AttributeStore::new();
                let mut records = NotifyRecords::new();
                let mut env = TopologyEnv {
                    all_attributes: &[],
                    attributes: &mut store,
                    key_of: &key_of,
                    records: &mut records,
                    next_key: 1u64 << 50,
                };
                make_topology_from_data(
                    &mut k,
                    &mut b,
                    &mut data,
                    &mut env,
                    MakeTopologyOptions {
                        make_uv_trim_curves: false,
                        add_analytics: false,
                        force_same_tol: false,
                        check_uv_trim_curves: false,
                        always_heal: false,
                        healer_enabled: false,
                    },
                )
                .unwrap();
                let keys = |o: CacheObject| match o {
                    CacheObject::Edge(e) | CacheObject::Curve(e) => key_of(EntityRef::Edge(e)),
                    CacheObject::Face(f) | CacheObject::Surface(f) => key_of(EntityRef::Face(f)),
                    CacheObject::Vertex(v) => key_of(EntityRef::Vertex(v)),
                    CacheObject::Other(r) => r,
                    _ => key_of(EntityRef::Brep(BrepId::from_index(0))),
                };
                let owner_of = |_: ObjectRef| None;
                let penv = ProjectEnv {
                    brep_id: BrepId::from_index(0),
                    keys: &keys,
                    owner_of: &owner_of,
                };
                let mut notify = NotifyEnv {
                    records: &mut records,
                    attributes: &mut store,
                    context: None,
                    key_of: &key_of,
                };
                let curve = Curve::from(BSplineCurve::from_nurb(
                    Net::new(
                        vec![
                            CPoint::new(-2., 0., 0., 1.),
                            CPoint::new(-1., 0., 0., 1.),
                            CPoint::new(1., 1., 0., 1.),
                            CPoint::new(2., 1., 0., 1.),
                        ],
                        1,
                        vec![0., 0., 1., 1. + gap, 2., 2.],
                    )
                    .unwrap(),
                    3,
                    true,
                ));
                let result = project_and_trim(
                    &mut k,
                    &mut b,
                    &penv,
                    &mut notify,
                    &curve,
                    Vec3::Z,
                    Vec3::new(-1.5, -1., 0.),
                    if trim == 0 {
                        TrimType::KeepPoint
                    } else {
                        TrimType::DeletePoint
                    },
                    &ProjectAndTrimOptions {
                        subset_faces: Some(&[]),
                        ..Default::default()
                    },
                );
                let faces = solidcore::topology::brep::p1::faces(&b, None).unwrap();
                if gap < 1e-8 {
                    let error =
                        result.expect_err("a refused side projection cannot classify a face");
                    assert_eq!(
                        error.kind(),
                        solidcore::error::ErrorKind::NotWithinTolerance,
                        "{file} gap={gap} trim={trim}: {error:?}"
                    );
                    // Empty projection subset makes this a classification-only operation:
                    // no face can be deleted before all required votes have been classified.
                    assert_eq!(faces.len(), 1, "{file} gap={gap} trim={trim}");
                } else {
                    result.unwrap();
                    assert_eq!(
                        faces.len(),
                        usize::from(trim == 0),
                        "{file} gap={gap} trim={trim}"
                    );
                }
            }
        }
    }
}
