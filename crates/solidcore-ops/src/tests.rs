//! The skeleton every family shares: primitives, copy and transform, counts and volume, `.smb`
//! bytes, the check, the result contract and the provenance hook.

use solidcore::kernel::UpstreamFixes;
use solidcore::topology::face::p1::flip_face_orientation;
use solidcore::{Brep, EntityRef};

use crate::body::{finish, Expect};
use crate::provenance::{self, Provenance, FIRST_ID};
use crate::session::{self, run, with_body};
use crate::{make_box, make_cone, make_cylinder, make_sphere, make_torus, Body, CheckLevel, Error};

/// Every test here runs with every fix on, as an application with all of them named would.
fn all_fixes() {
    session::set_fixes(UpstreamFixes::NAMES.iter().copied()).unwrap();
}

#[test]
fn a_box_measures_exactly() {
    all_fixes();
    let b = make_box([1.0, -2.0, 3.0], [10.0, 20.0, 30.0]).unwrap();
    assert_eq!(b.counts().unwrap(), [6, 12, 8]);
    let v = b.volume(1e-8).unwrap();
    assert!((v - 6000.0).abs() < 6000.0 * 1e-9, "{v}");
}

/// The accuracy goes to the kernel as given: one it refuses is its status.
#[test]
fn a_volume_accuracy_is_the_kernels_to_refuse() {
    all_fixes();
    let ball = make_sphere([1.0, 2.0, 3.0], 5.0).unwrap();
    let at = ball.volume(1e-8).unwrap();
    assert!(
        (at - 4.0 / 3.0 * std::f64::consts::PI * 125.0).abs() < 1e-6,
        "{at}"
    );
    assert_eq!(
        ball.volume(f64::INFINITY).err(),
        Some(Error::Status(crate::error::INVALID_INPUT))
    );
}

#[test]
fn hostile_numbers_are_refused() {
    all_fixes();
    assert_eq!(
        make_box([f64::NAN, 0.0, 0.0], [1.0; 3]).err(),
        Some(Error::Argument)
    );
    assert_eq!(
        make_box([0.0; 3], [1.0, -1.0, 1.0]).err(),
        Some(Error::Argument)
    );
    assert_eq!(
        make_cylinder([0.0; 3], f64::INFINITY, 1.0).err(),
        Some(Error::Argument)
    );
    assert_eq!(
        make_cone([0.0; 3], 2.0, 2.0, 1.0).err(),
        Some(Error::Argument)
    );
    assert_eq!(make_torus([0.0; 3], 1.0, 2.0).err(), Some(Error::Argument));
    let b = make_box([0.0; 3], [1.0; 3]).unwrap();
    let mut m = IDENTITY;
    m[3][0] = f64::NAN;
    assert_eq!(b.transform(&m).err(), Some(Error::Argument));
    let mut m = IDENTITY;
    m[0][3] = 0.5;
    assert_eq!(
        b.transform(&m).err(),
        Some(Error::Argument),
        "a projective fourth column"
    );
    let mut m = IDENTITY;
    m[3][3] = 2.0;
    assert_eq!(
        b.transform(&m).err(),
        Some(Error::Argument),
        "a projective fourth column"
    );
    assert_eq!(b.counts().unwrap(), [6, 12, 8]);
}

#[test]
fn every_primitive_is_valid_at_every_level() {
    all_fixes();
    let o = [1.0, 2.0, 3.0];
    let bodies = [
        ("box", make_box(o, [2.0, 3.0, 4.0]).unwrap()),
        ("cylinder", make_cylinder(o, 15.0, 50.0).unwrap()),
        ("sphere", make_sphere(o, 12.0).unwrap()),
        ("frustum", make_cone(o, 10.0, 4.0, 20.0).unwrap()),
        ("cone", make_cone(o, 10.0, 0.0, 20.0).unwrap()),
        ("torus", make_torus(o, 22.0, 4.0).unwrap()),
    ];
    for (name, b) in &bodies {
        for level in [CheckLevel::Basic, CheckLevel::Geometry, CheckLevel::All] {
            let r = b.check(level).unwrap();
            assert!(r.is_valid(), "{name}, {level:?}: {r}");
            assert_eq!(r.text, "");
        }
    }
}

/// A vertex whose tolerance is below the body's fails the check with its numbered line; a
/// vertex pointing at no vertexuse fails `ValidatePointers`, and `AssertValid` is then skipped.
#[test]
fn a_failed_check_reports_its_lines() {
    all_fixes();
    let b = make_box([0.0; 3], [1.0; 3]).unwrap();
    with_body("loosen", &b, |_, br| {
        let v = br.vertex_list()[0];
        br.vertex_mut(v)?.tolerance = 0.5;
        Ok(())
    })
    .unwrap();
    let r = b.check(CheckLevel::Geometry).unwrap();
    assert!(!r.is_valid());
    assert_eq!(
        r.failures,
        r.text.lines().filter(|l| l.contains("::Test #")).count(),
        "{r}"
    );
    assert!(!r.failed_tests().is_empty(), "{r}");
    let c = make_box([0.0; 3], [1.0; 3]).unwrap();
    with_body("unhook", &c, |_, br| {
        let v = br.vertex_list()[0];
        br.vertex_mut(v)?.vertexuses.clear();
        Ok(())
    })
    .unwrap();
    let r = c.check(CheckLevel::Basic).unwrap();
    assert!(!r.is_valid());
    assert!(
        r.text
            .contains("AssertValid skipped: ValidatePointers failed"),
        "{r}"
    );
}

#[test]
fn smb_bytes_round_trip_valid() {
    all_fixes();
    let b = make_torus([0.0; 3], 10.0, 3.0).unwrap();
    for ascii in [false, true] {
        let back = Body::from_smb(b.to_smb(ascii).unwrap(), ascii).unwrap();
        assert_eq!(back.counts().unwrap(), b.counts().unwrap());
        let r = back.check(CheckLevel::All).unwrap();
        assert!(r.is_valid(), "ascii {ascii}: {r}");
    }
}

const IDENTITY: [[f64; 4]; 4] = [
    [1.0, 0.0, 0.0, 0.0],
    [0.0, 1.0, 0.0, 0.0],
    [0.0, 0.0, 1.0, 0.0],
    [0.0, 0.0, 0.0, 1.0],
];

fn first_point(b: &Body) -> [f64; 3] {
    let br = b.lock();
    let p = br.vertex(br.vertex_list()[0]).unwrap().point;
    [p.x, p.y, p.z]
}

/// A translation, a rotation, a mirror and a non-uniform scale each move the body they are given
/// as their closed forms say: the volume times the determinant's size, a valid solid after.
#[test]
fn a_transform_moves_the_body_it_is_given() {
    all_fixes();
    let cases: [(&str, [[f64; 4]; 4], f64); 4] = [
        (
            "translation",
            [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [5.0, 6.0, 7.0, 1.0],
            ],
            1.0,
        ),
        (
            "quarter turn about z",
            [
                [0.0, 1.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0; 4].map(|_| 0.0),
            ],
            1.0,
        ),
        (
            "mirror in x",
            [
                [-1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
            1.0,
        ),
        (
            "scale 2, 3, 0.5",
            [
                [2.0, 0.0, 0.0, 0.0],
                [0.0, 3.0, 0.0, 0.0],
                [0.0, 0.0, 0.5, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
            3.0,
        ),
    ];
    for (name, mut m, det) in cases {
        m[3][3] = 1.0;
        let b = make_box([1.0, 1.0, 1.0], [1.0, 2.0, 3.0]).unwrap();
        let p = first_point(&b);
        b.transform(&m).unwrap();
        let q = first_point(&b);
        let image: Vec<f64> = (0..3)
            .map(|c| p[0] * m[0][c] + p[1] * m[1][c] + p[2] * m[2][c] + m[3][c])
            .collect();
        assert_eq!(q.to_vec(), image, "{name}");
        let v = b.volume(1e-8).unwrap();
        assert!((v - 6.0 * det).abs() < 1e-8 * 6.0 * det, "{name}: {v}");
        assert!(b.check(CheckLevel::All).unwrap().is_valid(), "{name}");
    }
}

#[derive(Clone, Debug, PartialEq)]
struct Origin(Vec<usize>);

impl Provenance for Origin {
    const ID: u64 = FIRST_ID + 41;
    fn split(&self, child: Option<&Self>) -> Option<Self> {
        Some(child.map_or_else(|| self.clone(), |c| self.joined(c)))
    }
    fn merge(&self, result: Option<&Self>) -> Option<Self> {
        Some(result.map_or_else(|| self.clone(), |r| self.joined(r)))
    }
    fn coincident(&self, onto: Option<&Self>) -> Option<Self> {
        Some(onto.map_or_else(|| self.clone(), |o| self.joined(o)))
    }
}

impl Origin {
    fn joined(&self, other: &Self) -> Self {
        let mut v = self.0.clone();
        v.extend(&other.0);
        v.sort_unstable();
        v.dedup();
        Origin(v)
    }
}

/// Another kind of data under the same ID.
#[derive(Clone, Debug, PartialEq)]
struct Impostor;

impl Provenance for Impostor {
    const ID: u64 = FIRST_ID + 41;
    fn split(&self, _: Option<&Self>) -> Option<Self> {
        None
    }
    fn merge(&self, _: Option<&Self>) -> Option<Self> {
        None
    }
    fn coincident(&self, _: Option<&Self>) -> Option<Self> {
        None
    }
}

#[test]
fn one_provenance_type_per_id() {
    provenance::take::<Origin>().unwrap();
    provenance::take::<Origin>().unwrap();
    assert_eq!(provenance::take::<Impostor>(), Err(Error::Argument));
}

/// Each face of `b` carries its index, each edge 100 + its index.
fn tag(b: &mut Brep, base: usize) {
    provenance::take::<Origin>().unwrap();
    for (i, f) in b.face_list().to_vec().into_iter().enumerate() {
        provenance::attach(b, EntityRef::Face(f), Origin(vec![base + i])).unwrap();
    }
    for (i, e) in b.edge_list().to_vec().into_iter().enumerate() {
        provenance::attach(b, EntityRef::Edge(e), Origin(vec![base + 100 + i])).unwrap();
    }
}

/// A copy keeps every attribute its source carries, in the same order (solidcore#1823's
/// review: a copy made under a fresh owner identity lost them).
#[test]
fn a_copy_keeps_the_provenance_on_its_faces_and_edges() {
    all_fixes();
    let b = make_box([0.0; 3], [1.0; 3]).unwrap();
    {
        let mut k = b.lock();
        tag(&mut k, 0);
        // Attaching again replaces.
        let f0 = k.face_list()[0];
        provenance::attach(&mut k, EntityRef::Face(f0), Origin(vec![0])).unwrap();
    }
    let c = b.copy().unwrap();
    let (src, copy) = (b.lock(), c.lock());
    assert_eq!((copy.face_list().len(), copy.edge_list().len()), (6, 12));
    for (i, (&f, &g)) in src.face_list().iter().zip(copy.face_list()).enumerate() {
        assert_eq!(
            provenance::read::<Origin>(&copy, EntityRef::Face(g)),
            Some(Origin(vec![i])),
            "face {i}"
        );
        assert_eq!(
            provenance::read::<Origin>(&src, EntityRef::Face(f)),
            Some(Origin(vec![i]))
        );
    }
    for (i, &e) in copy.edge_list().iter().enumerate() {
        assert_eq!(
            provenance::read::<Origin>(&copy, EntityRef::Edge(e)),
            Some(Origin(vec![100 + i])),
            "edge {i}"
        );
    }
}

/// A slab cut out of a box's middle splits each side face in two: both pieces carry the side's
/// data. A box joined to its neighbour along a face merges their coplanar sides: the merged face
/// carries both. Faces the boolean finds lying on one another join theirs onto the survivor.
#[test]
fn provenance_goes_through_splits_and_merges() {
    use solidcore::boolean::api_brep::boolean_difference;
    use solidcore::primitive::api_primitives::create_box;
    all_fixes();
    let split = run("split", |k| {
        let mut a = create_box(k, solidcore::Vec3::new(0.0, 0.0, 0.0), 3.0, 3.0, 3.0)?;
        tag(&mut a, 0);
        let slab = create_box(k, solidcore::Vec3::new(-1.0, -1.0, 1.0), 5.0, 5.0, 1.0)?;
        boolean_difference(k, a, slab)
    })
    .unwrap();
    let seen: Vec<Option<Origin>> = split
        .face_list()
        .iter()
        .map(|&f| provenance::read::<Origin>(&split, EntityRef::Face(f)))
        .collect();
    // the box's faces: its top and bottom whole, each side in two pieces; the slab's two faces
    // (the cut's walls) carry nothing, the slab carrying nothing
    assert_eq!(seen.len(), 2 + 2 * 4 + 2, "{seen:?}");
    assert_eq!(seen.iter().filter(|o| o.is_none()).count(), 2, "{seen:?}");
    let mut counts = [0usize; 6];
    for o in seen.iter().flatten() {
        assert_eq!(o.0.len(), 1, "{o:?}");
        counts[o.0[0]] += 1;
    }
    let mut sorted = counts;
    sorted.sort_unstable();
    assert_eq!(
        sorted,
        [1, 1, 2, 2, 2, 2],
        "pieces per face of the box: {counts:?}"
    );
}

/// The merge and coincidence rules, as the kernel's attribute store runs them: the survivor of
/// two merged entities carries both's data, and an entity found lying on another hands its data
/// on to that one. (The kernel's boolean merges coplanar faces and collinear edges without a
/// merge event, so no body edit reaches them yet: the data of the one removed is dropped.)
#[test]
fn provenance_goes_through_merges_and_coincidences() {
    use crate::provenance::Carried;
    use solidcore::geom::attribute::{Attribute, AttributeStore, Behavior, ObjectKey};
    let mut store = AttributeStore::default();
    let (a, b, c, d) = (ObjectKey(1), ObjectKey(2), ObjectKey(3), ObjectKey(4));
    let on = |store: &mut AttributeStore, x: ObjectKey, o: Origin| {
        let h = store
            .insert(Attribute::application(
                Origin::ID,
                Behavior::Copy,
                Carried(o),
            ))
            .unwrap();
        store.add_attribute(x, h).unwrap();
        h
    };
    let read = |store: &AttributeStore, x: ObjectKey| {
        let h = store.find_attribute(x, Origin::ID)?;
        store
            .get(h)?
            .application_value::<Carried<Origin>>()
            .map(|c| c.0.clone())
    };
    let ha = on(&mut store, a, Origin(vec![1]));
    on(&mut store, b, Origin(vec![2]));
    store.merge(ha, a, b, b).unwrap();
    assert_eq!(read(&store, b), Some(Origin(vec![1, 2])));
    let hc = on(&mut store, c, Origin(vec![3]));
    on(&mut store, d, Origin(vec![4]));
    store.coincident(hc, c, d).unwrap();
    assert_eq!(read(&store, d), Some(Origin(vec![3, 4])));
}

/// A solid input must give a solid: an empty result is `Empty`, a box turned inside out is
/// `NotSolid`, a box is a solid. What is expected comes from the input.
#[test]
fn the_result_contract() {
    all_fixes();
    let out = run("contract", |k| {
        let mut input = solidcore::primitive::api_primitives::create_box(
            k,
            solidcore::Vec3::new(0.0, 0.0, 0.0),
            1.0,
            1.0,
            1.0,
        )?;
        let expect = Expect::of(k, &mut input)?;
        let empty = finish(k, Brep::new_with_infinite_region()?, expect)?.err();
        let mut inside_out = solidcore::primitive::api_primitives::create_box(
            k,
            solidcore::Vec3::new(0.0, 0.0, 0.0),
            1.0,
            1.0,
            1.0,
        )?;
        for f in inside_out.face_list().to_vec() {
            flip_face_orientation(&mut inside_out, f)?;
        }
        let flipped = finish(k, inside_out, expect)?.err();
        let same = finish(k, input, expect)?.is_ok();
        let any = solidcore::primitive::api_primitives::create_box(
            k,
            solidcore::Vec3::new(0.0, 0.0, 0.0),
            1.0,
            1.0,
            1.0,
        )?;
        let sheet_ok = finish(k, any, Expect::FACES)?.is_ok();
        Ok(Some((empty, flipped, same, sheet_ok)))
    })
    .unwrap();
    assert_eq!(
        out,
        Some((Some(Error::Empty), Some(Error::NotSolid), true, true))
    );
}
