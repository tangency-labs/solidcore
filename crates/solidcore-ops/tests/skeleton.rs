//! The skeleton every family shares: the session, primitives, copy and transform, counts and
//! volume, `.smb` bytes, the check, the result contract, the index and the provenance hook.

use std::sync::{Arc, Mutex};

use solidcore::{Brep, EntityRef, SolidError};
use solidcore_ops::provenance::{self, FIRST_ID, Provenance};
use solidcore_ops::session::{self, Failure};
use solidcore_ops::{Body, CheckLevel, Error, Expect, Index, finish, make_box, make_cone, make_cylinder, make_sphere, make_torus};

/// Every test here runs with every fix on, as an application with all of them named would.
#[allow(clippy::unwrap_used)]
fn all_fixes() {
    session::set_fixes(solidcore::kernel::UpstreamFixes::NAMES.iter().copied()).unwrap();
}

#[test]
fn a_box_measures_exactly() {
    all_fixes();
    let b = make_box([1.0, -2.0, 3.0], [10.0, 20.0, 30.0]).unwrap();
    assert_eq!(b.counts().unwrap(), [6, 12, 8]);
    let v = b.volume(1e-8).unwrap();
    assert!((v - 6000.0).abs() < 6000.0 * 1e-9, "{v}");
}

/// An accuracy finer than 1e-8 or not positive is 1e-8: a sphere's volume at 0, -1 and 1e-15 is
/// its volume at 1e-8, on its closed form.
#[test]
fn a_volume_accuracy_out_of_range_is_clamped() {
    all_fixes();
    let ball = make_sphere([1.0, 2.0, 3.0], 5.0).unwrap();
    let at = ball.volume(1e-8).unwrap();
    assert!((at - 4.0 / 3.0 * std::f64::consts::PI * 125.0).abs() < 1e-6, "{at}");
    for acc in [0.0, -1.0, 1e-15] {
        assert_eq!(ball.volume(acc).ok(), Some(at), "accuracy {acc}");
    }
}

#[test]
fn hostile_numbers_are_refused() {
    all_fixes();
    assert_eq!(make_box([f64::NAN, 0.0, 0.0], [1.0; 3]).err(), Some(Error::Argument));
    assert_eq!(make_box([0.0; 3], [1.0, -1.0, 1.0]).err(), Some(Error::Argument));
    assert_eq!(make_cylinder([0.0; 3], f64::INFINITY, 1.0).err(), Some(Error::Argument));
    assert_eq!(make_cone([0.0; 3], 2.0, 2.0, 1.0).err(), Some(Error::Argument));
    assert_eq!(make_torus([0.0; 3], 1.0, 2.0).err(), Some(Error::Argument));
    let b = make_box([0.0; 3], [1.0; 3]).unwrap();
    assert_eq!(b.volume(f64::INFINITY).err(), Some(Error::Argument));
    let mut m = [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0; 4]];
    m[3][0] = f64::NAN;
    assert_eq!(b.transform(&m).err(), Some(Error::Argument));
}

#[test]
fn a_panic_is_an_error_and_the_kernel_is_fresh_after() {
    all_fixes();
    assert!(!session::in_call());
    let r: solidcore_ops::Result<()> = session::run("test", |_| panic!("boom"));
    assert!(matches!(r, Err(Error::Internal(ref m)) if m.contains("boom")), "{r:?}");
    assert!(matches!(session::crash_test(), Err(Error::Internal(ref m)) if m.contains("crash test")));
    assert!(make_box([0.0; 3], [1.0; 3]).is_ok());
    assert!(!session::in_call());
    let inside = session::run("inside", |_| Ok(session::in_call() && session::current_call().as_deref() == Some("inside"))).unwrap();
    assert!(inside);
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

/// A translation moves the body and keeps its volume; the transform changes the body it is given.
#[test]
fn a_transform_moves_the_body_it_is_given() {
    all_fixes();
    let b = make_box([0.0; 3], [1.0, 2.0, 3.0]).unwrap();
    let m = [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [5.0, 6.0, 7.0, 1.0]];
    b.transform(&m).unwrap();
    let first = b.lock().vertex_list()[0];
    let p = b.lock().vertex(first).unwrap().point;
    assert!(p.x >= 5.0 && p.y >= 6.0 && p.z >= 7.0, "{p:?}");
    assert!((b.volume(1e-8).unwrap() - 6.0).abs() < 1e-9);
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

/// A copy keeps every attribute its source carries, in the same record order (solidcore#1823's
/// review: a copy made under a fresh owner identity lost them).
#[test]
fn a_copy_keeps_the_provenance_on_its_faces_and_edges() {
    all_fixes();
    let b = make_box([0.0; 3], [1.0; 3]).unwrap();
    {
        let mut k = b.lock();
        let faces = k.face_list().to_vec();
        for (i, f) in faces.into_iter().enumerate() {
            provenance::attach(&mut k, EntityRef::Face(f), Origin(vec![i])).unwrap();
        }
        let edges = k.edge_list().to_vec();
        for (i, e) in edges.into_iter().enumerate() {
            provenance::attach(&mut k, EntityRef::Edge(e), Origin(vec![100 + i])).unwrap();
        }
        // Attaching again replaces.
        let f0 = k.face_list()[0];
        provenance::attach(&mut k, EntityRef::Face(f0), Origin(vec![0])).unwrap();
    }
    let c = b.copy().unwrap();
    let (src, copy) = (b.lock(), c.lock());
    let (si, ci) = (Index::of_kernel_order(&src), Index::of_kernel_order(&copy));
    assert_eq!((ci.faces.len(), ci.edges.len()), (6, 12));
    for i in 0..6 {
        let got: Option<Origin> = provenance::read(&copy, EntityRef::Face(ci.face(i).unwrap()));
        assert_eq!(got, Some(Origin(vec![i])), "face {i}");
        assert_eq!(provenance::read::<Origin>(&src, EntityRef::Face(si.face(i).unwrap())), Some(Origin(vec![i])));
    }
    for i in 0..12 {
        let got: Option<Origin> = provenance::read(&copy, EntityRef::Edge(ci.edge(i).unwrap()));
        assert_eq!(got, Some(Origin(vec![100 + i])), "edge {i}");
    }
    assert_eq!(ci.face(6).err(), Some(Error::Argument));
}

#[test]
fn the_result_contract() {
    all_fixes();
    let empty = Brep::new_with_infinite_region().unwrap();
    assert_eq!(finish(empty, Expect::Faces).unwrap().err(), Some(Error::Empty));
    let b = make_box([0.0; 3], [1.0; 3]).unwrap().into_inner();
    assert!(finish(b, Expect::Solid).unwrap().is_ok());
}

/// A failed call hands its name, its reason and its inputs as `.smb` bytes to the callback the
/// application installed; nothing is written without one.
#[test]
fn a_failed_call_reaches_the_failure_hook_with_its_inputs() {
    type Seen = Vec<(String, String, Vec<(&'static str, Vec<u8>)>)>;
    all_fixes();
    let seen: Arc<Mutex<Seen>> = Arc::default();
    let s = Arc::clone(&seen);
    session::set_failure_hook(Some(Box::new(move |f: &Failure<'_>| {
        if f.call == "refused on purpose" {
            let inputs = f.inputs.iter().map(|(w, b)| (*w, b.clone().unwrap())).collect();
            s.lock().unwrap().push((f.call.to_owned(), f.why.to_owned(), inputs));
        }
    })));
    let b = make_box([0.0; 3], [1.0; 3]).unwrap();
    let r: solidcore_ops::Result<()> = session::with_body("refused on purpose", &b, |_, _| Err(SolidError::invalid_input("test", "refused")));
    assert!(matches!(r, Err(Error::Kernel(_))), "{r:?}");
    session::set_failure_hook(None);
    let seen = seen.lock().unwrap();
    assert_eq!(seen.len(), 1);
    let (call, why, inputs) = &seen[0];
    assert_eq!(call, "refused on purpose");
    assert!(why.contains("refused"), "{why}");
    assert_eq!(inputs.len(), 1);
    assert_eq!(inputs[0].0, "brep");
    let back = Body::from_smb(inputs[0].1.clone(), false).unwrap();
    assert_eq!(back.counts().unwrap(), [6, 12, 8]);
}
