//! The Part 21 parser over a STEP corpus: every file parses, the units of its representation
//! contexts, and its Part 42 curves and surfaces.
//!
//! Ignored by default; the corpus is not in the repository. Run it with the directory of
//! `.stp` / `.step` files (read in place, recursively) in `SOLIDCORE_STEP_CORPUS`:
//!
//! ```sh
//! SOLIDCORE_STEP_CORPUS=/path/to/corpus cargo test --release -p solidcore-record \
//!     --test step_p21_corpus -- --ignored --nocapture
//! ```
//!
//! Each file prints one line: its entity count, duplicate ids and schema, or its error. Files
//! are read as bytes and decoded lossily (some writers put Latin-1 bytes in strings), as a
//! caller holding the file's bytes would. The test fails when any file does not parse. The
//! STEP files committed to the repository are parsed by `repository_files_parse`, which always
//! runs.
//!
//! `corpus_units` (ignored too, same variable and command with `corpus_units` as the filter)
//! reads every geometric representation context with assigned units in each file and prints
//! one line per file and a table of how often each length and plane angle unit occurs, and
//! how many contexts state an uncertainty. It fails when a file does not parse. Files with no
//! context whose length unit can be read are listed and counted, not failed: a file written
//! with Part 21 short names (`SUNT`, `LNGUNT`, `GMRPCN`, …) or one without geometry has none.
//!
//! `corpus_geometry` (ignored too, same variable, `corpus_geometry` as the filter) reads every
//! curve and surface entity of each file with `step::read::Reader`, in the units of the file's
//! first geometric context: curves in model space with `curve`, the curves of a parameter space
//! through their `PCURVE` with `pcurve`, surfaces with `surface`. It prints one line per file
//! and, per group and entity kind, how many were read and how many refused, with the reasons.
//! It also lists the entities where a check the record no longer makes (a direction shorter than
//! 1e-12, unit factors equal within 1e-12) would have decided otherwise. It fails when a file
//! does not parse; refusals are counted, not failed.
//!
//! `corpus_topology` (ignored too, same variable, `corpus_topology` as the filter) maps every
//! solid item and loose face of each file with `step::topology::read_solids` and checks it
//! against the file's own census, counted from its entities without the reader as #2593 counts
//! it: one body per `MANIFOLD_SOLID_BREP`, `BREP_WITH_VOIDS`, `FACETED_BREP` and
//! `SHELL_BASED_SURFACE_MODEL`, its faces the lengths of its shells' face lists. Every body is
//! read or refused with its reason; a body read has as many faces as its shells list; and no file
//! maps more distinct faces, edges or vertices than it has `ADVANCED_FACE`/`FACE_SURFACE`,
//! `EDGE_CURVE` and `VERTEX_POINT` entities. It prints one line per file, the totals, the reasons
//! for refusals and the files whose distinct counts differ from the file's. It fails when a file
//! does not parse or a check fails; refusals are counted, not failed.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use solidcore_record::step::p21;
use solidcore_record::step::units::{self, UnitKind};

fn files(dir: &Path, out: &mut Vec<PathBuf>) {
    let Ok(entries) = std::fs::read_dir(dir) else {
        return;
    };
    let mut paths = Vec::new();
    for entry in entries.flatten() {
        paths.push(entry.path());
    }
    paths.sort();
    for path in paths {
        if path.is_dir() {
            files(&path, out);
        } else if path
            .extension()
            .and_then(std::ffi::OsStr::to_str)
            .is_some_and(|e| matches!(e.to_ascii_lowercase().as_str(), "stp" | "step"))
        {
            out.push(path);
        }
    }
}

/// Parse each file; print one line per file and a total. Returns the failures.
fn parse_all(root: &Path, paths: &[PathBuf]) -> Vec<String> {
    let mut failures = Vec::new();
    let mut parsed = 0usize;
    let mut entities = 0usize;
    for path in paths {
        let name = path
            .strip_prefix(root)
            .unwrap_or(path)
            .display()
            .to_string();
        let bytes = match std::fs::read(path) {
            Ok(bytes) => bytes,
            Err(e) => {
                println!("{name}: unreadable: {e}");
                failures.push(format!("{name}: unreadable: {e}"));
                continue;
            }
        };
        match p21::parse(&String::from_utf8_lossy(&bytes)) {
            Ok(ex) => {
                parsed += 1;
                entities += ex.entities.len();
                println!(
                    "{name}: {} entities, {} duplicate ids, schema {:?}",
                    ex.entities.len(),
                    ex.duplicates.len(),
                    ex.schemas
                );
            }
            Err(e) => {
                println!("{name}: ERROR {e}");
                failures.push(format!("{name}: {e}"));
            }
        }
    }
    println!(
        "parsed {parsed} of {} files, {entities} entities",
        paths.len()
    );
    failures
}

#[test]
#[ignore = "needs a STEP corpus in SOLIDCORE_STEP_CORPUS"]
fn corpus_parses() {
    let Some(dir) = std::env::var_os("SOLIDCORE_STEP_CORPUS") else {
        eprintln!("SOLIDCORE_STEP_CORPUS is not set");
        return;
    };
    let root = Path::new(&dir);
    let mut paths = Vec::new();
    files(root, &mut paths);
    assert!(!paths.is_empty(), "no STEP files under {dir:?}");
    let failures = parse_all(root, &paths);
    assert!(
        failures.is_empty(),
        "files that did not parse:\n{}",
        failures.join("\n")
    );
}

/// A length factor's name (millimetres per unit).
fn length_name(f: f64) -> String {
    match f {
        1.0 => "mm".into(),
        25.4 => "inch".into(),
        1000.0 => "m".into(),
        10.0 => "cm".into(),
        304.8 => "foot".into(),
        _ => format!("{f} mm"),
    }
}

/// A plane angle factor's name (radians per unit). A degree is the number the file writes:
/// `degree` when it is π/180 to 12 digits or more, else the written factor, so a file that
/// rounds π/180 shorter (`0.0174532925`) shows in the table as what it is.
fn angle_name(f: f64) -> String {
    if f == 1.0 {
        "rad".into()
    } else if (f - std::f64::consts::PI / 180.0).abs() <= 1e-14 {
        "degree".into()
    } else {
        format!("{f} rad")
    }
}

#[test]
#[ignore = "needs a STEP corpus in SOLIDCORE_STEP_CORPUS"]
fn corpus_units() {
    let Some(dir) = std::env::var_os("SOLIDCORE_STEP_CORPUS") else {
        eprintln!("SOLIDCORE_STEP_CORPUS is not set");
        return;
    };
    let root = Path::new(&dir);
    let mut paths = Vec::new();
    files(root, &mut paths);
    assert!(!paths.is_empty(), "no STEP files under {dir:?}");
    let mut failures = Vec::new();
    let mut unread = Vec::new();
    let mut lengths: BTreeMap<String, usize> = BTreeMap::new();
    let mut angles: BTreeMap<String, usize> = BTreeMap::new();
    let (mut contexts, mut with_uncertainty, mut files_read) = (0usize, 0usize, 0usize);
    for path in &paths {
        let name = path
            .strip_prefix(root)
            .unwrap_or(path)
            .display()
            .to_string();
        let ex = match std::fs::read(path)
            .map_err(|e| e.to_string())
            .and_then(|b| p21::parse(&String::from_utf8_lossy(&b)))
        {
            Ok(ex) => ex,
            Err(e) => {
                println!("{name}: ERROR {e}");
                failures.push(format!("{name}: {e}"));
                continue;
            }
        };
        let mut here = Vec::new();
        let mut readable = false;
        for (&id, e) in &ex.entities {
            if !(e.has("GEOMETRIC_REPRESENTATION_CONTEXT") && e.has("GLOBAL_UNIT_ASSIGNED_CONTEXT"))
            {
                continue;
            }
            // Whether the context names a length unit that reads as a factor, rather than
            // falling back to millimetres.
            readable |= e
                .record("GLOBAL_UNIT_ASSIGNED_CONTEXT")
                .and_then(<[p21::Param]>::first)
                .and_then(p21::Param::as_list)
                .is_some_and(|l| {
                    l.iter().filter_map(p21::Param::as_ref_id).any(|u| {
                        matches!(units::unit_factor(&ex, u),
                            Some((UnitKind::Length, f)) if f.is_finite() && f > 0.0)
                    })
                });
            let c = units::context(&ex, id);
            contexts += 1;
            with_uncertainty += usize::from(c.uncertainty.is_some());
            let (l, a) = (length_name(c.units.len), angle_name(c.units.ang));
            *lengths.entry(l.clone()).or_default() += 1;
            *angles.entry(a.clone()).or_default() += 1;
            here.push(format!("#{id} {l}/{a} uncertainty {:?}", c.uncertainty));
        }
        println!("{name}: {} contexts: {}", here.len(), here.join(", "));
        if readable {
            files_read += 1;
        } else {
            println!("{name}: no context with a readable length unit");
            unread.push(name);
        }
    }
    println!(
        "{files_read} of {} files have a context with a readable length unit; {contexts} contexts, {with_uncertainty} with an uncertainty",
        paths.len()
    );
    println!("length units (contexts): {lengths:?}");
    println!("plane angle units (contexts): {angles:?}");
    println!(
        "no context with a readable length unit: {} files {unread:?}",
        unread.len()
    );
    assert!(
        failures.is_empty(),
        "files that did not parse:\n{}",
        failures.join("\n")
    );
}

#[test]
fn repository_files_parse() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let paths = [
        "crates/difftest/src/units/topology_brep_constructor_tetra.stp",
        "crates/step-brep/tests/fixtures/discord-2026-10-30-fillet.step",
    ]
    .map(|p| root.join(p));
    let failures = parse_all(&root, &paths);
    assert!(
        failures.is_empty(),
        "files that did not parse:\n{}",
        failures.join("\n")
    );
}

/// Curve entity names (simple instances) the corpus count looks at: every one the reader maps,
/// and the other ISO 10303-42 curves a file may hold, which it refuses.
const CURVES: &[&str] = &[
    "B_SPLINE_CURVE",
    "B_SPLINE_CURVE_WITH_KNOTS",
    "BEZIER_CURVE",
    "BOUNDED_SURFACE_CURVE",
    "CIRCLE",
    "CLOTHOID",
    "COMPOSITE_CURVE",
    "COMPOSITE_CURVE_ON_SURFACE",
    "CURVE_REPLICA",
    "ELLIPSE",
    "HYPERBOLA",
    "INTERSECTION_CURVE",
    "LINE",
    "OFFSET_CURVE_2D",
    "OFFSET_CURVE_3D",
    "PARABOLA",
    "POLYLINE",
    "QUASI_UNIFORM_CURVE",
    "SEAM_CURVE",
    "SURFACE_CURVE",
    "TRIMMED_CURVE",
    "UNIFORM_CURVE",
];

/// Surface entity names, likewise.
const SURFACES: &[&str] = &[
    "B_SPLINE_SURFACE",
    "B_SPLINE_SURFACE_WITH_KNOTS",
    "BEZIER_SURFACE",
    "CONICAL_SURFACE",
    "CURVE_BOUNDED_SURFACE",
    "CYLINDRICAL_SURFACE",
    "DEGENERATE_TOROIDAL_SURFACE",
    "OFFSET_SURFACE",
    "ORIENTED_SURFACE",
    "PLANE",
    "QUASI_UNIFORM_SURFACE",
    "RECTANGULAR_COMPOSITE_SURFACE",
    "RECTANGULAR_TRIMMED_SURFACE",
    "SPHERICAL_SURFACE",
    "SURFACE_OF_LINEAR_EXTRUSION",
    "SURFACE_OF_REVOLUTION",
    "SURFACE_REPLICA",
    "TOROIDAL_SURFACE",
    "UNIFORM_SURFACE",
];

/// Whether `e` is a curve or a surface for the count, and its kind: the entity name, or for a
/// complex instance `(B_SPLINE_CURVE …)` with `rational` when it is.
fn geometry_kind(e: &p21::Entity) -> Option<(bool, String)> {
    if e.is_complex() {
        let rational = if e.has("RATIONAL_B_SPLINE_CURVE") || e.has("RATIONAL_B_SPLINE_SURFACE") {
            " rational"
        } else {
            ""
        };
        let is_curve =
            e.has("B_SPLINE_CURVE") || e.records.iter().any(|r| CURVES.contains(&r.name.as_str()));
        let is_surface = e.has("B_SPLINE_SURFACE")
            || e.records
                .iter()
                .any(|r| SURFACES.contains(&r.name.as_str()));
        let names: Vec<&str> = e
            .records
            .iter()
            .map(|r| r.name.as_str())
            .filter(|n| CURVES.contains(n) || SURFACES.contains(n))
            .collect();
        let kind = format!("({}){rational}", names.join(" "));
        return match (is_curve, is_surface) {
            (true, false) => Some((true, kind)),
            (false, true) => Some((false, kind)),
            _ => None,
        };
    }
    let n = e.name();
    if CURVES.contains(&n) {
        Some((true, n.to_string()))
    } else if SURFACES.contains(&n) {
        Some((false, n.to_string()))
    } else {
        None
    }
}

/// An error message with its entity ids replaced by `#n`, so equal reasons count together.
fn reason(e: &str) -> String {
    let mut out = String::new();
    let mut chars = e.chars().peekable();
    while let Some(c) = chars.next() {
        out.push(c);
        if c == '#' && chars.peek().is_some_and(char::is_ascii_digit) {
            out.push('n');
            while chars.peek().is_some_and(char::is_ascii_digit) {
                chars.next();
            }
        }
    }
    out
}

/// The curves of a file that lie in a parameter space: the items of its
/// `DEFINITIONAL_REPRESENTATION`s and the curves they refer to (a 2D trimmed curve's basis …).
fn curves_2d(ex: &p21::Exchange) -> std::collections::BTreeSet<u64> {
    fn refs(p: &p21::Param, out: &mut Vec<u64>) {
        match p {
            p21::Param::Ref(r) => out.push(*r),
            p21::Param::List(l) | p21::Param::Typed(_, l) => l.iter().for_each(|x| refs(x, out)),
            _ => {}
        }
    }
    let mut stack = Vec::new();
    for e in ex.entities.values() {
        if e.name() == "DEFINITIONAL_REPRESENTATION"
            && let Some(items) = e.params().get(1)
        {
            refs(items, &mut stack);
        }
    }
    let mut seen = std::collections::BTreeSet::new();
    while let Some(id) = stack.pop() {
        let Some(e) = ex.get(id) else { continue };
        if !matches!(geometry_kind(e), Some((true, _))) || !seen.insert(id) {
            continue;
        }
        for r in &e.records {
            r.params.iter().for_each(|p| refs(p, &mut stack));
        }
    }
    seen
}

/// Read or refused, by reason.
#[derive(Default)]
struct Tally {
    read: usize,
    refused: BTreeMap<String, usize>,
}

impl Tally {
    fn add<T>(&mut self, r: &Result<T, String>) {
        match r {
            Ok(_) => self.read += 1,
            Err(e) => *self.refused.entry(reason(e)).or_default() += 1,
        }
    }
    fn refused(&self) -> usize {
        self.refused.values().sum()
    }
}

#[test]
#[ignore = "needs a STEP corpus in SOLIDCORE_STEP_CORPUS"]
fn corpus_geometry() {
    use solidcore_record::step::read::Reader;
    let Some(dir) = std::env::var_os("SOLIDCORE_STEP_CORPUS") else {
        eprintln!("SOLIDCORE_STEP_CORPUS is not set");
        return;
    };
    let root = Path::new(&dir);
    let mut paths = Vec::new();
    files(root, &mut paths);
    assert!(!paths.is_empty(), "no STEP files under {dir:?}");
    let mut failures = Vec::new();
    // (group, kind) -> tally, groups "curve", "surface", "pcurve" (by its 2D curve's kind).
    let mut kinds: BTreeMap<(&str, String), Tally> = BTreeMap::new();
    let mut mixed_units = Vec::new();
    let mut removed = Vec::new();
    for path in &paths {
        let name = path
            .strip_prefix(root)
            .unwrap_or(path)
            .display()
            .to_string();
        let ex = match std::fs::read(path)
            .map_err(|e| e.to_string())
            .and_then(|b| p21::parse(&String::from_utf8_lossy(&b)))
        {
            Ok(ex) => ex,
            Err(e) => {
                println!("{name}: ERROR {e}");
                failures.push(format!("{name}: {e}"));
                continue;
            }
        };
        // The units of the file's first geometric context (the topology reader, which comes
        // later, reads each representation in its own context).
        let contexts: Vec<units::Units> = ex
            .entities
            .iter()
            .filter(|(_, e)| {
                e.has("GEOMETRIC_REPRESENTATION_CONTEXT") && e.has("GLOBAL_UNIT_ASSIGNED_CONTEXT")
            })
            .map(|(&id, _)| units::context(&ex, id).units)
            .collect();
        let u = contexts.first().copied().unwrap_or_default();
        if contexts.iter().any(|c| *c != u) {
            mixed_units.push(name.clone());
        }
        // Where a removed check would have decided otherwise: SolidCraft refused a direction
        // whose length is at most 1e-12 (here: only one that is exactly zero), and took the two
        // unit factors of a conic's parameter space as equal within 1e-12 relative (here: only
        // when equal).
        for (&id, e) in &ex.entities {
            if e.name() == "DIRECTION"
                && let Some(r) = e.params().get(1).and_then(p21::Param::as_list)
            {
                let sq: f64 = r.iter().filter_map(p21::Param::as_f64).map(|x| x * x).sum();
                if sq > 0.0 && sq <= 1e-24 {
                    removed.push(format!("{name}: #{id} a direction of length at most 1e-12"));
                }
            }
        }
        if u.len != u.ang && (u.len - u.ang).abs() <= 1e-12 * u.len.abs().max(u.ang.abs()) {
            removed.push(format!(
                "{name}: length and angle factors equal within 1e-12"
            ));
        }
        let flat = curves_2d(&ex);
        let mut reader = Reader::new(&ex, u);
        let mut here = (Tally::default(), Tally::default(), Tally::default());
        for (&id, e) in &ex.entities {
            if e.name() == "PCURVE" {
                let kind = e
                    .params()
                    .get(2)
                    .and_then(p21::Param::as_ref_id)
                    .and_then(|r| ex.get(r))
                    .and_then(|r| r.params().get(1)?.as_list()?.first()?.as_ref_id())
                    .and_then(|c| ex.get(c))
                    .and_then(geometry_kind)
                    .map_or_else(|| "?".to_string(), |k| k.1);
                let r = reader.pcurve(id);
                here.2.add(&r);
                kinds.entry(("pcurve", kind)).or_default().add(&r);
                continue;
            }
            let Some((is_curve, kind)) = geometry_kind(e) else {
                continue;
            };
            if is_curve {
                if flat.contains(&id) {
                    continue;
                }
                let r = reader.curve(id);
                here.0.add(&r);
                kinds.entry(("curve", kind)).or_default().add(&r);
            } else {
                let r = reader.surface(id);
                here.1.add(&r);
                kinds.entry(("surface", kind)).or_default().add(&r);
            }
        }
        println!(
            "{name}: curves {} read, {} refused; surfaces {} read, {} refused; pcurves {} read, {} refused",
            here.0.read,
            here.0.refused(),
            here.1.read,
            here.1.refused(),
            here.2.read,
            here.2.refused()
        );
    }
    println!("by kind (group, kind: read / refused, with reasons):");
    let mut totals: BTreeMap<&str, (usize, usize)> = BTreeMap::new();
    for ((group, kind), t) in &kinds {
        println!("  {group} {kind}: {} / {}", t.read, t.refused());
        for (why, n) in &t.refused {
            println!("      {n} x {why}");
        }
        let g = totals.entry(group).or_default();
        g.0 += t.read;
        g.1 += t.refused();
    }
    println!("totals (read / refused): {totals:?}");
    println!(
        "files whose contexts have different units (read in the first): {} {mixed_units:?}",
        mixed_units.len()
    );
    println!(
        "entities a removed check would have decided otherwise: {} {removed:?}",
        removed.len()
    );
    assert!(
        failures.is_empty(),
        "files that did not parse:\n{}",
        failures.join("\n")
    );
}

/// The solid items #2593 counts as bodies.
const BODIES: [&str; 4] = [
    "MANIFOLD_SOLID_BREP",
    "BREP_WITH_VOIDS",
    "FACETED_BREP",
    "SHELL_BASED_SURFACE_MODEL",
];

/// A body's faces as its entities list them, without the topology reader: the length of each
/// shell's face list (an oriented shell's inner shell), summed over the shells it holds.
fn census_faces(ex: &p21::Exchange, body: &p21::Entity) -> usize {
    let p = body.params();
    let one = |i: usize| p.get(i).and_then(p21::Param::as_ref_id).into_iter();
    let many = |i: usize| {
        p.get(i)
            .and_then(p21::Param::as_list)
            .unwrap_or(&[])
            .iter()
            .filter_map(p21::Param::as_ref_id)
    };
    let shells: Vec<u64> = match body.name() {
        "BREP_WITH_VOIDS" => one(1).chain(many(2)).collect(),
        "SHELL_BASED_SURFACE_MODEL" => many(1).collect(),
        _ => one(1).collect(),
    };
    let mut n = 0;
    for s in shells {
        let Some(mut e) = ex.get(s) else { continue };
        if e.name().starts_with("ORIENTED_")
            && let Some(inner) = e
                .params()
                .get(2)
                .and_then(p21::Param::as_ref_id)
                .and_then(|i| ex.get(i))
        {
            e = inner;
        }
        n += e
            .params()
            .get(1)
            .and_then(p21::Param::as_list)
            .map_or(0, <[p21::Param]>::len);
    }
    n
}

#[test]
#[ignore = "needs a STEP corpus in SOLIDCORE_STEP_CORPUS"]
fn corpus_topology() {
    use solidcore_record::step::topology::{self, SolidKind};
    use std::collections::BTreeSet;
    let Some(dir) = std::env::var_os("SOLIDCORE_STEP_CORPUS") else {
        eprintln!("SOLIDCORE_STEP_CORPUS is not set");
        return;
    };
    let root = Path::new(&dir);
    let mut paths = Vec::new();
    files(root, &mut paths);
    assert!(!paths.is_empty(), "no STEP files under {dir:?}");
    let mut failures = Vec::new();
    let mut refused: BTreeMap<String, usize> = BTreeMap::new();
    let mut loose_refused: BTreeMap<String, usize> = BTreeMap::new();
    // Census bodies, bodies read, bodies refused, loose faces read and refused; listed and
    // mapped faces of the bodies read; mapped and file faces, edges and vertices.
    let mut total = [0usize; 13];
    let mut open_solids = Vec::new();
    let mut differ = Vec::new();
    for path in &paths {
        let name = path
            .strip_prefix(root)
            .unwrap_or(path)
            .display()
            .to_string();
        let ex = match std::fs::read(path)
            .map_err(|e| e.to_string())
            .and_then(|b| p21::parse(&String::from_utf8_lossy(&b)))
        {
            Ok(ex) => ex,
            Err(e) => {
                println!("{name}: ERROR {e}");
                failures.push(format!("{name}: {e}"));
                continue;
            }
        };
        let r = topology::read_solids(&ex);
        let count = |names: &[&str]| {
            ex.entities
                .values()
                .filter(|e| !e.is_complex() && names.contains(&e.name()))
                .count()
        };
        let bodies: Vec<u64> = ex
            .entities
            .iter()
            .filter(|(_, e)| !e.is_complex() && BODIES.contains(&e.name()))
            .map(|(id, _)| *id)
            .collect();
        let read: BTreeMap<u64, &topology::StepSolid> =
            r.solids.iter().map(|s| (s.item, s)).collect();
        let failed: BTreeSet<u64> = r.failures.iter().map(|f| f.item).collect();
        let (mut body_read, mut body_refused, mut census, mut mapped) = (0, 0, 0, 0);
        for &b in &bodies {
            let pieces = topology::items(&ex, b);
            if pieces.iter().all(|p| read.contains_key(p)) {
                body_read += 1;
                let c = census_faces(&ex, ex.get(b).unwrap());
                let m: usize = pieces.iter().map(|p| read[p].faces.len()).sum();
                census += c;
                mapped += m;
                if c != m {
                    failures.push(format!("{name}: body #{b}: {m} faces mapped, {c} listed"));
                }
            } else {
                body_refused += 1;
                assert!(
                    pieces.iter().any(|p| failed.contains(p)),
                    "{name}: body #{b} neither read nor refused"
                );
            }
        }
        let loose: BTreeSet<u64> = topology::loose_faces(&ex).into_iter().collect();
        let mut loose_failed = 0;
        for f in &r.failures {
            let tally = if loose.contains(&f.item) {
                loose_failed += 1;
                &mut loose_refused
            } else {
                &mut refused
            };
            *tally.entry(reason(&f.message)).or_default() += 1;
            println!("{name}: #{} refused: {}", f.item, f.message);
        }
        let loose_read = r
            .solids
            .iter()
            .filter(|s| s.kind == SolidKind::LooseFace)
            .count();
        for s in &r.solids {
            if s.kind != SolidKind::LooseFace
                && s.shells
                    .iter()
                    .chain(&s.voids)
                    .any(|sh| sh.declared_closed && !sh.closed)
            {
                open_solids.push(format!("{name}: #{}", s.item));
            }
        }
        let mut faces = BTreeSet::new();
        let mut edges = BTreeSet::new();
        let mut vertices = BTreeSet::new();
        for s in &r.solids {
            faces.extend(s.faces.iter().map(|f| f.entity));
            edges.extend(s.edges.iter().map(|e| e.entity));
            vertices.extend(s.vertices.iter().map(|v| v.entity));
        }
        let file = [
            count(&["ADVANCED_FACE", "FACE_SURFACE"]),
            count(&["EDGE_CURVE"]),
            count(&["VERTEX_POINT"]),
        ];
        let got = [faces.len(), edges.len(), vertices.len()];
        for (k, what) in ["faces", "edges", "vertices"].iter().enumerate() {
            if got[k] > file[k] {
                failures.push(format!(
                    "{name}: {} {what} mapped, {} in the file",
                    got[k], file[k]
                ));
            }
        }
        if got != file {
            differ.push(format!(
                "{name}: faces {}/{}, edges {}/{}, vertices {}/{} (mapped/file); {} bodies refused",
                got[0], file[0], got[1], file[1], got[2], file[2], body_refused
            ));
        }
        println!(
            "{name}: bodies {} (read {body_read}, refused {body_refused}), loose faces {loose_read} read \
             {loose_failed} refused; faces {mapped}/{census} (mapped/listed, bodies read); \
             faces {}/{}, edges {}/{}, vertices {}/{} (mapped/file)",
            bodies.len(),
            got[0],
            file[0],
            got[1],
            file[1],
            got[2],
            file[2]
        );
        let add = [
            bodies.len(),
            body_read,
            body_refused,
            loose_read,
            loose_failed,
            census,
            mapped,
            got[0],
            file[0],
            got[1],
            file[1],
            got[2],
            file[2],
        ];
        for (t, a) in total.iter_mut().zip(add) {
            *t += a;
        }
    }
    println!("files: {}", paths.len());
    println!(
        "bodies: {} in the census, {} read, {} refused; loose faces: {} read, {} refused",
        total[0], total[1], total[2], total[3], total[4]
    );
    println!(
        "faces of the bodies read: {} mapped, {} listed",
        total[6], total[5]
    );
    println!(
        "distinct entities mapped / in the file: faces {}/{}, edges {}/{}, vertices {}/{}",
        total[7], total[8], total[9], total[10], total[11], total[12]
    );
    println!(
        "bodies read with a shell the file calls closed whose edges are not all used twice: {}",
        open_solids.len()
    );
    for o in &open_solids {
        println!("  {o}");
    }
    println!("bodies refused, by reason: {refused:?}");
    println!("loose faces refused, by reason: {loose_refused:?}");
    println!("files whose mapped faces, edges or vertices differ from the file's count:");
    for d in &differ {
        println!("  {d}");
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
