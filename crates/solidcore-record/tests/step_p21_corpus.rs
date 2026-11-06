//! The Part 21 parser over a STEP corpus: every file parses, and the units of its representation
//! contexts.
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
//! reads every geometric representation context in each file (`GEOMETRIC_REPRESENTATION_CONTEXT`,
//! or `GMRPCN` in Part 21 short names; not a `PARAMETRIC_REPRESENTATION_CONTEXT`, which has no
//! units) and prints one line per file, a table of how often each
//! length and plane angle unit occurs (its name and factor), how many contexts state a distance
//! uncertainty, and each context whose length unit is unreadable, with why. It fails when a file
//! does not parse; unreadable units are listed and counted, not failed (a file written with
//! Part 21 short names has them, and the importer decides).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use solidcore_record::step::p21;
use solidcore_record::step::units::{self, Unit, Unreadable};

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

/// An error as text, whichever type the parser's error is.
fn text(e: &impl std::fmt::Display) -> String {
    e.to_string()
}

/// A unit as a table key: its name and factor, or why it is unreadable.
fn key(u: &Result<Unit, Unreadable>) -> String {
    match u {
        Ok(u) => format!("{} {}", u.name, u.factor),
        Err(e) => format!("unreadable {:?}", e.why),
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
    let mut unreadable = Vec::new();
    let mut without = Vec::new();
    let mut lengths: BTreeMap<String, usize> = BTreeMap::new();
    let mut angles: BTreeMap<String, usize> = BTreeMap::new();
    let (mut contexts, mut with_uncertainty) = (0usize, 0usize);
    for path in &paths {
        let name = path
            .strip_prefix(root)
            .unwrap_or(path)
            .display()
            .to_string();
        let ex = match std::fs::read(path)
            .map_err(|e| text(&e))
            .and_then(|b| p21::parse(&String::from_utf8_lossy(&b)).map_err(|e| text(&e)))
        {
            Ok(ex) => ex,
            Err(e) => {
                println!("{name}: ERROR {e}");
                failures.push(format!("{name}: {e}"));
                continue;
            }
        };
        let mut here = Vec::new();
        let mut bad = 0usize;
        for (&id, e) in &ex.entities {
            // A parametric context (the (u, v) space of a pcurve) has no units.
            if !(e.has("GEOMETRIC_REPRESENTATION_CONTEXT") || e.has("GMRPCN"))
                || e.has("PARAMETRIC_REPRESENTATION_CONTEXT")
            {
                continue;
            }
            let c = units::context(&ex, id);
            contexts += 1;
            with_uncertainty += usize::from(matches!(c.uncertainty, Some(Ok(_))));
            let (l, a) = (key(&c.units.len), key(&c.units.ang));
            *lengths.entry(l.clone()).or_default() += 1;
            *angles.entry(a.clone()).or_default() += 1;
            if let Err(u) = c.units.len {
                bad += 1;
                println!(
                    "{name}: #{id} length unit unreadable at #{}: {:?}",
                    u.entity, u.why
                );
            }
            here.push(format!("#{id} {l} / {a} uncertainty {:?}", c.uncertainty));
        }
        println!("{name}: {} contexts: {}", here.len(), here.join(", "));
        if here.is_empty() {
            without.push(name);
        } else if bad > 0 {
            unreadable.push(format!("{name} ({bad} of {})", here.len()));
        }
    }
    println!(
        "{} files; {contexts} contexts, {with_uncertainty} with a distance uncertainty",
        paths.len()
    );
    println!("length units (contexts): {lengths:?}");
    println!("plane angle units (contexts): {angles:?}");
    println!(
        "files with an unreadable length unit: {} {unreadable:?}",
        unreadable.len()
    );
    println!(
        "files with no geometric representation context: {} {without:?}",
        without.len()
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
