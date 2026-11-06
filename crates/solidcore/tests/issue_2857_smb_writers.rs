//! #2857: a body read from an `.smb` file can hold kinds no solidcore operation makes
//! (hyperbolas, composite and offset curves, STEP surfaces, surfaces in volumes), and
//! writing it back must give the record it was read from.
//!
//! Each record of the reader units (`curve.composite_curve_p2`, `curve.offset_curve`,
//! `surface.offset_surface`, `surface.surf_of_revolution_p3`, `surface.srf_in_volume`)
//! is frozen under `tests/data/issue_2857_smb_writers/<unit>.json`. A record upstream's
//! `WriteToDB` wrote (`upstream: true`: the fixtures' read cases filled by the C++
//! reference, and generated cases at seed 2857) is read and written back in its own
//! format, which must give the same bytes in binary, and in ASCII the same text with
//! every double as the file writer writes the record's value (`same_text`). Every
//! record, the hand-written fixtures too,
//! is then written in binary and read again: the second binary write is the first byte
//! for byte, and both objects evaluate to the same bits on a grid of their domain.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::float_cmp
)] // test helpers

use serde_json::Value;
use solidcore::curve::curve::p2::{Curve, COMPOSITE_CURVE_TYPE, OFFSET_CURVE_TYPE};
use solidcore::geom::vector2d::Point2;
use solidcore::surface::surface::p2::{
    OFFSET_SURFACE_TYPE, SRF_IN_VOLUME_TYPE, SURF_OF_REVOLUTION_TYPE,
};
use solidcore::surface::surface::Surface;
use solidcore::topology::database_io::{fixed16, DatabaseReader, DatabaseWriter, FileType};
use solidcore::topology::geometry_write::{write_curve, write_surface};
use solidcore::Kernel;

/// One frozen record.
struct Record {
    name: String,
    dim: u64,
    binary: bool,
    version: u64,
    upstream: bool,
    bytes: Vec<u8>,
}

/// A curve or a surface read from a record.
enum Object {
    Curve(Curve),
    Surface(Surface),
}

fn records(unit: &str) -> Vec<Record> {
    let path = format!(
        "{}/tests/data/issue_2857_smb_writers/{unit}.json",
        env!("CARGO_MANIFEST_DIR")
    );
    let v: Value = serde_json::from_str(&std::fs::read_to_string(&path).unwrap()).unwrap();
    v.as_array()
        .unwrap()
        .iter()
        .map(|r| {
            let hex = r["hex"].as_str().unwrap();
            let bytes = (0..hex.len())
                .step_by(2)
                .map(|i| u8::from_str_radix(&hex[i..i + 2], 16).unwrap())
                .collect();
            Record {
                name: r["name"].as_str().unwrap().to_string(),
                dim: r["dim"].as_u64().unwrap_or(3),
                binary: r["binary"].as_bool().unwrap(),
                version: r["version"].as_u64().unwrap(),
                upstream: r["upstream"].as_bool().unwrap(),
                bytes,
            }
        })
        .collect()
}

const fn file_type(binary: bool) -> FileType {
    if binary {
        FileType::Binary
    } else {
        FileType::Ascii
    }
}

fn read(
    k: &Kernel,
    type_code: i64,
    curve: bool,
    r: &Record,
    bytes: Vec<u8>,
    binary: bool,
) -> Object {
    let mut db = DatabaseReader::new();
    db.open_for_read(bytes, file_type(binary));
    if curve {
        let c = Curve::read_from_db(k, type_code, &mut db, r.dim, r.version);
        Object::Curve(c.unwrap_or_else(|e| panic!("{}: read: {e:?}", r.name)))
    } else {
        let s = Surface::read_from_db(k, type_code, &mut db, r.version);
        Object::Surface(s.unwrap_or_else(|e| panic!("{}: read: {e:?}", r.name)))
    }
}

fn write(o: &Object, r: &Record, binary: bool) -> Vec<u8> {
    let mut db = DatabaseWriter::open_for_write(file_type(binary));
    let written = match o {
        Object::Curve(c) => write_curve(c, &mut db, r.version),
        Object::Surface(s) => write_surface(s, &mut db, r.version),
    };
    written.unwrap_or_else(|e| panic!("{}: write: {e:?}", r.name));
    db.into_bytes()
}

/// The object's points on a 5 x 5 grid (5 parameters for a curve) of its natural
/// domain, as bits; an evaluation error as its text.
fn grid(k: &Kernel, o: &Object) -> Vec<String> {
    let at = |lo: f64, hi: f64, i: usize| lo + (hi - lo) * (i as f64) / 4.0;
    let bits = |p: solidcore::Result<solidcore::Point3>| match p {
        Ok(p) => format!(
            "{:x} {:x} {:x}",
            p.x.to_bits(),
            p.y.to_bits(),
            p.z.to_bits()
        ),
        Err(e) => format!("{e:?}"),
    };
    match o {
        Object::Curve(c) => {
            let ivl = c.natural_interval();
            (0..5)
                .map(|i| bits(c.evaluate_point(k, at(ivl.min(), ivl.max(), i))))
                .collect()
        }
        Object::Surface(s) => {
            let d = match s.natural_uv_domain_with_nurb(k) {
                Ok(d) => d,
                Err(e) => return vec![format!("{e:?}")],
            };
            let mut out = Vec::new();
            for i in 0..5 {
                for j in 0..5 {
                    let uv = Point2::new(at(d.u_min(), d.u_max(), i), at(d.v_min(), d.v_max(), j));
                    out.push(bits(s.evaluate_point(k, uv)));
                }
            }
            out
        }
    }
}

/// An ASCII record written back as `written`, against the C++ reference's `record`:
/// the same lines, labels and integers, and every double written as the file writer
/// writes the record's value ([`fixed16`]: the reference's shim writes doubles with 17
/// significant digits instead, so its integral doubles have no decimals).
fn same_text(name: &str, written: &[u8], record: &[u8]) {
    let (w, r) = (
        String::from_utf8_lossy(written),
        String::from_utf8_lossy(record),
    );
    let (wl, rl): (Vec<&str>, Vec<&str>) = (w.lines().collect(), r.lines().collect());
    assert_eq!(wl.len(), rl.len(), "{name}: other line count");
    for (i, (a, b)) in wl.iter().zip(&rl).enumerate() {
        let (at, bt): (Vec<&str>, Vec<&str>) = (a.split(' ').collect(), b.split(' ').collect());
        // a volume's label follows its value without a space ("0.5BendVolume ...")
        let split = |t: &str| {
            let k = t
                .find(|c: char| c.is_ascii_alphabetic() && c != 'e' && c != 'E')
                .unwrap_or(t.len());
            (t[..k].to_string(), t[k..].to_string())
        };
        let same = at.len() == bt.len()
            && at.iter().zip(&bt).all(|(x, y)| {
                let ((xn, xl), (yn, yl)) = (split(x), split(y));
                xl == yl && (xn == yn || yn.parse::<f64>().is_ok_and(|v| xn == fixed16(v)))
            });
        assert!(same, "{name}: line {i} written as\n{a:?}, record\n{b:?}");
    }
}

/// Every record of `unit` (read as `type_code`) through the checks of the module docs.
fn round_trip(unit: &str, type_code: i64, curve: bool) {
    let k = Kernel::default();
    let records = records(unit);
    assert!(!records.is_empty(), "{unit}: no records");
    for r in &records {
        let o = read(&k, type_code, curve, r, r.bytes.clone(), r.binary);
        if r.upstream {
            let back = write(&o, r, r.binary);
            if r.binary {
                assert!(back == r.bytes, "{}: written back to other bytes", r.name);
            } else {
                same_text(&r.name, &back, &r.bytes);
            }
        }
        let b1 = write(&o, r, true);
        let o2 = read(&k, type_code, curve, r, b1.clone(), true);
        assert!(
            write(&o2, r, true) == b1,
            "{}: binary write not stable",
            r.name
        );
        assert_eq!(grid(&k, &o), grid(&k, &o2), "{}: other geometry", r.name);
    }
}

#[test]
fn composite_curves_write_back_their_records() {
    round_trip("curve.composite_curve_p2", COMPOSITE_CURVE_TYPE, true);
}

#[test]
fn offset_curves_write_back_their_records() {
    round_trip("curve.offset_curve", OFFSET_CURVE_TYPE, true);
}

#[test]
fn offset_surfaces_write_back_their_records() {
    round_trip("surface.offset_surface", OFFSET_SURFACE_TYPE, false);
}

#[test]
fn surfaces_of_revolution_write_back_their_records() {
    round_trip(
        "surface.surf_of_revolution_p3",
        SURF_OF_REVOLUTION_TYPE,
        false,
    );
}

#[test]
fn surfaces_in_volumes_write_back_their_records() {
    round_trip("surface.srf_in_volume", SRF_IN_VOLUME_TYPE, false);
}
