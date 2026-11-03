//! Shared corpus helpers (`#[path]`-included): the exporter's own geometric validation
//! volumes of a STEP file, read through `step_import`'s Part 21 parser and converted to
//! mm³, and the file's body and face census, so that expected values come from the file
//! and never from the reader under test.
//!
//! A validation volume is written as the CAx-IF Recommended Practices for Geometric and
//! Assembly Validation Properties (release 4.2, §4 and its summary table) describe it:
//! a `PROPERTY_DEFINITION` named `'geometric validation property'`, linked by a
//! `PROPERTY_DEFINITION_REPRESENTATION` to a `REPRESENTATION` whose items include a
//! `MEASURE_REPRESENTATION_ITEM` with a `VOLUME_MEASURE` value and a volume unit. The
//! unit is a `DERIVED_UNIT` (or `VOLUME_UNIT`) of `DERIVED_UNIT_ELEMENT`s whose length
//! units, an `SI_UNIT` or a `CONVERSION_BASED_UNIT` such as `'INCH'` (25.4 mm), are
//! resolved through the file (ISO 10303-41); the volume in mm³ is the written value
//! times the product of each element's mm-per-unit raised to its exponent.
#![allow(
    clippy::expect_used,
    clippy::indexing_slicing,
    clippy::panic,
    dead_code
)] // test helpers; each test uses a subset

use std::path::{Path, PathBuf};

use step_import::part21::{self, Instance, Param, Record, StepFile};
use step_import::Limits;

/// One geometric validation volume of a file.
#[derive(Debug, Clone, PartialEq)]
pub struct ValidationVolume {
    /// The `PROPERTY_DEFINITION`'s `#id`.
    pub property: u64,
    /// What the property is defined on (its `definition` attribute's `#id`).
    pub definition: u64,
    /// The property's description, e.g. `'volume of NIST_CTC_05_ASME1_RD'` (may be empty).
    pub description: String,
    /// The value as written, in the file's volume unit.
    pub written: f64,
    /// The value in mm³.
    pub mm3: f64,
}

/// The corpus directory, or `None` when `SOLIDCORE_STEP_CORPUS` is unset.
pub fn corpus_dir() -> Option<PathBuf> {
    std::env::var_os("SOLIDCORE_STEP_CORPUS").map(PathBuf::from)
}

/// The file named `name` anywhere under the corpus `dir` (searched recursively, in
/// sorted order, so the STEP Tools and NIST directories can sit side by side or be
/// given on their own), read in place.
pub fn corpus_file(dir: &Path, name: &str) -> Vec<u8> {
    let path = find(dir, name).unwrap_or_else(|| {
        panic!(
            "{name} is not under {}; fetch the corpus with scripts/fetch-steptools-corpus.sh \
             and scripts/fetch-nist-pmi-corpus.sh",
            dir.display()
        )
    });
    std::fs::read(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

/// The first path named `name` under `dir`, depth first in sorted order.
fn find(dir: &Path, name: &str) -> Option<PathBuf> {
    let mut paths: Vec<PathBuf> = std::fs::read_dir(dir)
        .ok()?
        .flatten()
        .map(|entry| entry.path())
        .collect();
    paths.sort();
    for path in paths {
        if path.is_dir() {
            if let Some(found) = find(&path, name) {
                return Some(found);
            }
        } else if path.file_name().and_then(std::ffi::OsStr::to_str) == Some(name) {
            return Some(path);
        }
    }
    None
}

/// `bytes` parsed by `step_import`'s Part 21 parser.
pub fn parse(bytes: &[u8]) -> StepFile {
    part21::parse(bytes, &Limits::default()).expect("the file parses")
}

/// Every geometric validation volume of `file`, in file order, converted to mm³.
///
/// # Errors
///
/// When a volume's unit cannot be resolved to lengths.
pub fn validation_volumes(file: &StepFile) -> Result<Vec<ValidationVolume>, String> {
    let mut out = Vec::new();
    for instance in file.entities.values() {
        let Some(link) = instance.record("PROPERTY_DEFINITION_REPRESENTATION") else {
            continue;
        };
        let (Some(property), Some(representation)) = (ref_at(link, 0), ref_at(link, 1)) else {
            continue;
        };
        let Some(definition) = record(file, property, "PROPERTY_DEFINITION") else {
            continue;
        };
        if str_at(definition, 0) != Some("geometric validation property") {
            continue;
        }
        let Some(representation) = record(file, representation, "REPRESENTATION") else {
            continue;
        };
        let items = representation
            .params
            .get(1)
            .and_then(Param::as_list)
            .unwrap_or_default();
        for item in items.iter().filter_map(Param::as_ref_id) {
            let Some(measure) = record(file, item, "MEASURE_REPRESENTATION_ITEM") else {
                continue;
            };
            let Some(Param::Typed(kind, value)) = measure.params.get(1) else {
                continue;
            };
            if kind != "VOLUME_MEASURE" {
                continue;
            }
            let written = value
                .as_f64()
                .ok_or_else(|| format!("#{item}: VOLUME_MEASURE is not a number"))?;
            let unit = ref_at(measure, 2).ok_or_else(|| format!("#{item}: no unit"))?;
            let factor = volume_mm3(file, unit).map_err(|error| format!("#{item}: {error}"))?;
            out.push(ValidationVolume {
                property,
                definition: ref_at(definition, 2).unwrap_or_default(),
                description: str_at(definition, 1).unwrap_or_default().to_string(),
                written,
                mm3: written * factor,
            });
        }
    }
    Ok(out)
}

/// The one validation volume of `file` in mm³: `None` where the file has none, else the
/// value every validation volume of the file gives (a file often writes the same value
/// twice, on the part and on its shape aspect).
///
/// # Errors
///
/// When a unit cannot be resolved, or the file's validation volumes disagree.
pub fn validation_volume(file: &StepFile) -> Result<Option<f64>, String> {
    let volumes = validation_volumes(file)?;
    let Some(first) = volumes.first() else {
        return Ok(None);
    };
    for other in &volumes {
        if (other.mm3 - first.mm3).abs() > 1e-12 * first.mm3.abs().max(other.mm3.abs()) {
            return Err(format!(
                "validation volumes disagree: #{} {} mm³, #{} {} mm³",
                first.property, first.mm3, other.property, other.mm3
            ));
        }
    }
    Ok(Some(first.mm3))
}

/// The body entities `step_import` maps, as the file declares them: (`#id`, entity,
/// faces of its shells). Faces are counted over every shell the body holds (a solid's
/// outer shell and voids, a surface model's shells), so they are the file's, not the
/// reader's.
pub fn body_census(file: &StepFile) -> Vec<(u64, String, usize)> {
    const BODIES: [&str; 4] = [
        "MANIFOLD_SOLID_BREP",
        "BREP_WITH_VOIDS",
        "FACETED_BREP",
        "SHELL_BASED_SURFACE_MODEL",
    ];
    let mut out = Vec::new();
    for (&id, instance) in &file.entities {
        let Instance::Simple(body) = instance else {
            continue;
        };
        if !BODIES.contains(&body.name.as_str()) {
            continue;
        }
        let mut shells: Vec<u64> = Vec::new();
        for param in body.params.iter().skip(1) {
            match param {
                Param::Ref(shell) => shells.push(*shell),
                Param::List(list) => shells.extend(list.iter().filter_map(Param::as_ref_id)),
                _ => {}
            }
        }
        let faces = shells.iter().map(|&shell| shell_faces(file, shell)).sum();
        out.push((id, body.name.clone(), faces));
    }
    out
}

/// The number of faces of the shell `#id` (through an `ORIENTED_CLOSED_SHELL` or
/// `ORIENTED_OPEN_SHELL` to the shell it orients).
fn shell_faces(file: &StepFile, id: u64) -> usize {
    let Some(Instance::Simple(shell)) = file.entities.get(&id) else {
        return 0;
    };
    match shell.name.as_str() {
        "CLOSED_SHELL" | "OPEN_SHELL" => shell
            .params
            .get(1)
            .and_then(Param::as_list)
            .map_or(0, <[Param]>::len),
        "ORIENTED_CLOSED_SHELL" | "ORIENTED_OPEN_SHELL" => {
            ref_at(shell, 2).map_or(0, |inner| shell_faces(file, inner))
        }
        _ => 0,
    }
}

/// mm³ per unit of the volume unit `#id`.
fn volume_mm3(file: &StepFile, id: u64) -> Result<f64, String> {
    let instance = file
        .entities
        .get(&id)
        .ok_or_else(|| format!("unit #{id} missing"))?;
    let derived = instance
        .record("DERIVED_UNIT")
        .or_else(|| instance.record("VOLUME_UNIT"))
        .ok_or_else(|| {
            format!(
                "unit #{id} ({}) is not a derived unit",
                instance.type_name()
            )
        })?;
    let elements = derived
        .params
        .first()
        .and_then(Param::as_list)
        .ok_or_else(|| format!("unit #{id}: no elements"))?;
    let mut factor = 1.0;
    let mut dimension = 0.0;
    for element in elements.iter().filter_map(Param::as_ref_id) {
        let element = record(file, element, "DERIVED_UNIT_ELEMENT")
            .ok_or_else(|| format!("unit #{id}: element is not a DERIVED_UNIT_ELEMENT"))?;
        let unit = ref_at(element, 0).ok_or_else(|| format!("unit #{id}: element unit"))?;
        let exponent = element
            .params
            .get(1)
            .and_then(Param::as_f64)
            .ok_or_else(|| format!("unit #{id}: element exponent"))?;
        factor *= length_mm(file, unit, 0)?.powf(exponent);
        dimension += exponent;
    }
    if dimension != 3.0 {
        return Err(format!(
            "unit #{id} has length dimension {dimension}, not 3"
        ));
    }
    Ok(factor)
}

/// mm per unit of the length unit `#id`: an `SI_UNIT` of metres with its prefix, or a
/// `CONVERSION_BASED_UNIT` whose `MEASURE_WITH_UNIT` gives it in another length unit.
fn length_mm(file: &StepFile, id: u64, depth: usize) -> Result<f64, String> {
    if depth > 8 {
        return Err(format!("length unit #{id}: conversion chain too deep"));
    }
    let instance = file
        .entities
        .get(&id)
        .ok_or_else(|| format!("length unit #{id} missing"))?;
    if instance.record("LENGTH_UNIT").is_none() {
        return Err(format!(
            "unit #{id} ({}) is not a length unit",
            instance.type_name()
        ));
    }
    if let Some(si) = instance.record("SI_UNIT") {
        if si.params.get(1).and_then(Param::as_enum) != Some("METRE") {
            return Err(format!("length unit #{id}: SI unit is not METRE"));
        }
        let prefix = match si.params.first() {
            Some(Param::Null) | None => 1.0,
            Some(prefix) => match prefix.as_enum() {
                Some("MILLI") => 1e-3,
                Some("CENTI") => 1e-2,
                Some("DECI") => 1e-1,
                Some("KILO") => 1e3,
                Some("MICRO") => 1e-6,
                Some("NANO") => 1e-9,
                other => return Err(format!("length unit #{id}: prefix {other:?}")),
            },
        };
        return Ok(prefix * 1000.0);
    }
    if let Some(conversion) = instance.record("CONVERSION_BASED_UNIT") {
        let measure = ref_at(conversion, 1)
            .and_then(|m| file.entities.get(&m))
            .ok_or_else(|| format!("length unit #{id}: no conversion factor"))?;
        let measure = measure
            .record("LENGTH_MEASURE_WITH_UNIT")
            .or_else(|| measure.record("MEASURE_WITH_UNIT"))
            .ok_or_else(|| format!("length unit #{id}: conversion is not a measure"))?;
        let value = measure
            .params
            .first()
            .and_then(Param::as_f64)
            .ok_or_else(|| format!("length unit #{id}: conversion value"))?;
        let unit = ref_at(measure, 1).ok_or_else(|| format!("length unit #{id}: base unit"))?;
        return Ok(value * length_mm(file, unit, depth + 1)?);
    }
    Err(format!(
        "length unit #{id} ({}) is neither SI nor conversion based",
        instance.type_name()
    ))
}

/// The record `name` of the instance `#id`.
fn record<'a>(file: &'a StepFile, id: u64, name: &str) -> Option<&'a Record> {
    file.entities.get(&id).and_then(|i| i.record(name))
}

/// The reference at parameter `index` of `record`.
fn ref_at(record: &Record, index: usize) -> Option<u64> {
    record.params.get(index).and_then(Param::as_ref_id)
}

/// The string at parameter `index` of `record`.
fn str_at(record: &Record, index: usize) -> Option<&str> {
    match record.params.get(index) {
        Some(Param::Str(s)) => Some(s),
        _ => None,
    }
}

/// `|a - b| / max(|a|, |b|)`.
pub fn relative(a: f64, b: f64) -> f64 {
    (a - b).abs() / a.abs().max(b.abs()).max(f64::MIN_POSITIVE)
}
