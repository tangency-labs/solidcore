//! The cases of #2580, pinned against today's reader so that the port into the reader of
//! #2578 has an acceptance test: faces that wind round a torus or sphere, multi-winding
//! loops, loopless faces and edges across a closed curve, in nine corpus files.
//!
//! Expected values are the file's own, never the reader's output:
//!
//! - **Volumes** where the file carries its exporter's geometric validation volume
//!   (read through `step_import`'s parser and converted from the file's unit, see
//!   [`validation`]), or where the geometry has a closed form (`ftc_11_rb`, below).
//! - **Counts and closedness** everywhere, from the file's census: one body per
//!   `MANIFOLD_SOLID_BREP` / `BREP_WITH_VOIDS` / `FACETED_BREP` /
//!   `SHELL_BASED_SURFACE_MODEL`, and the faces of the shells each body holds
//!   ([`validation::body_census`]). Every solid is expected closed, a surface model open.
//!
//! | file (corpus path) | bodies | faces | closed | volume expected (mm³) |
//! |---|---|---|---|---|
//! | `nist_ftc_06_asme1_rd` (AP203 geometry only) | 1 | 144 | 1 | none in the file |
//! | `nist_stc_06_asme1_ap242-e3` | 1 | 144 | 1 | the file writes `VOLUME_MEASURE(0.0)` |
//! | `nist_ftc_11_asme1_rb` (AP203 geometry only) | 1 | 6 | 1 | 5129.386136 (closed form) |
//! | `nist_ftc_11_asme1_ap242-e2` | 1 | 42 | 1 | none in the file |
//! | `turbine` (STEP Tools AP203) | 1 | 1241 | 1 | none in the file |
//! | `nist_ctc_05_asme1_rd` (AP203 geometry only) | 2 | 209 | 1 | 774.6912251592 in³ = 12 694 914.687 |
//! | `nist_ctc_05_asme1_ap203` (AP203 with PMI) | 1 | 156 | 1 | 774.691593579675 in³ = 12 694 920.724 |
//! | `nist_ctc_05_asme1_ap242-e1` | 1 | 156 | 1 | none in the file |
//! | `nist_stc_10_asme1_ap242-e2` | 1 | 256 | 1 | 188 312.244641244 |
//!
//! `ctc_05_rd`'s second body is its `SHELL_BASED_SURFACE_MODEL` of one face, an open
//! sheet. `stc_06`'s validation properties are all zero (volume, area, and a centroid at
//! the origin): an exporter placeholder, not a measurement, so no volume is asserted for
//! it. The issue lists `ctc_05` twice; the corpus has three exports of it, and all three
//! are pinned (the two AP203 ones on their validation volumes).
//!
//! **Tolerance.** [`VALIDATION_TOLERANCE`] is the smallest round relative tolerance the
//! three validation volumes meet with today's reader (measured 1.1e-6, 1.2e-7 and 3.4e-6);
//! it is not derived from the files. The CAx-IF Recommended Practices for Geometric and
//! Assembly Validation Properties (release 4.2, §4.11) leave the threshold to the
//! business process and count a volume within 1 % as a successful exchange, which is far
//! too loose to catch a wrong face here. It is an open question to the maintainers.
//!
//! **Known failing today.** `ftc_11_rb` is a washer of six faces, two of them tori whose
//! faces run over the top of the tube; today's reader closes it at 3179.68 mm³, 38 % below
//! its closed form ([`ftc_11_rb_volume_is_wrong_known_failing`]).
//!
//! The corpus tests are ignored by default; the corpus is not in the repository. Fetch it
//! with `scripts/fetch-steptools-corpus.sh` and `scripts/fetch-nist-pmi-corpus.sh` and run
//! them with a directory holding both in `SOLIDCORE_STEP_CORPUS` (searched recursively):
//!
//! ```sh
//! SOLIDCORE_STEP_CORPUS=corpus cargo test --release -p step-brep \
//!     --test issue_2580_winding_faces -- --ignored --nocapture
//! ```

#![allow(clippy::expect_used, clippy::indexing_slicing, clippy::panic)] // test helpers

#[path = "support/validation.rs"]
mod validation;

use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::kernel::Kernel;
use solidcore::mass::integrator::properties::{brep_volume, PropertyAttributes};
use solidcore::topology::Brep;
use step_brep::{import_with_report, is_closed_solid, ImportReport};
use step_import::Limits;

use validation::{corpus_dir, corpus_file, relative};

/// The relative tolerance of a measured volume against the file's validation volume
/// (see the module docs: chosen from today's measurements, an open question).
const VALIDATION_TOLERANCE: f64 = 5e-6;

/// The relative tolerance of a measured volume against a closed form: the integrator's
/// accuracy is 1e-9, so this is left at 1e-7 for the trims of the curved faces.
const CLOSED_FORM_TOLERANCE: f64 = 1e-7;

/// What a file's volume is checked against.
#[derive(Debug, Clone, Copy)]
enum Volume {
    /// The file's geometric validation volume.
    Validation,
    /// A closed form of the file's geometry, in mm³.
    ClosedForm(fn() -> f64),
    /// Nothing: counts and closedness only.
    None,
}

/// One corpus file: (file name, bodies, faces, closed solids, volume).
type Case = (&'static str, usize, usize, usize, Volume);

/// The files of #2580 with their census (see the module docs).
const CASES: [Case; 9] = [
    ("nist_ftc_06_asme1_rd.stp", 1, 144, 1, Volume::None),
    ("nist_stc_06_asme1_ap242-e3.stp", 1, 144, 1, Volume::None),
    (
        "nist_ftc_11_asme1_rb.stp",
        1,
        6,
        1,
        Volume::ClosedForm(ftc_11_rb_closed_form), // known failing
    ),
    ("nist_ftc_11_asme1_ap242-e2.stp", 1, 42, 1, Volume::None),
    ("turbine.stp", 1, 1241, 1, Volume::None),
    ("nist_ctc_05_asme1_rd.stp", 2, 209, 1, Volume::Validation),
    ("nist_ctc_05_asme1_ap203.stp", 1, 156, 1, Volume::Validation),
    ("nist_ctc_05_asme1_ap242-e1.stp", 1, 156, 1, Volume::None),
    (
        "nist_stc_10_asme1_ap242-e2.stp",
        1,
        256,
        1,
        Volume::Validation,
    ),
];

/// The volume of `brep`, or `None` where the integrator gives none.
fn volume(brep: &Brep) -> Option<f64> {
    let kernel = Kernel::default();
    let store = AttributeStore::default();
    let attrs = PropertyAttributes {
        store: &store,
        key_of: &|_| ObjectKey(0),
    };
    let mut copy = brep.clone();
    brep_volume(&kernel, &mut copy, 1e-9, &attrs).ok()
}

/// (bodies, faces over all bodies, closed bodies) of `report`.
fn counts(report: &ImportReport) -> (usize, usize, usize) {
    (
        report.bodies.len(),
        report.bodies.iter().map(|b| b.face_list().len()).sum(),
        report.bodies.iter().filter(|b| is_closed_solid(b)).count(),
    )
}

/// The summed volume of the closed bodies of `report`.
fn closed_volume(report: &ImportReport) -> f64 {
    report
        .bodies
        .iter()
        .filter(|b| is_closed_solid(b))
        .map(|b| volume(b).expect("a closed body has a volume"))
        .sum()
}

/// `name` read from the corpus and imported by today's reader.
fn import(dir: &std::path::Path, name: &str) -> (Vec<u8>, ImportReport) {
    let bytes = corpus_file(dir, name);
    let report = import_with_report(&bytes, &Limits::default())
        .unwrap_or_else(|error| panic!("{name}: import failed: {}", error.message));
    (bytes, report)
}

/// `∫ (x + c) √(a² − x²) dx`: the r-moment of a circular arc of radius `a` about a centre
/// at radius `c`, as a function of `x = r − c` (the antiderivative, closed form).
fn arc_moment(x: f64, a: f64, c: f64) -> f64 {
    let s = (a * a - x * x).max(0.0).sqrt();
    -(s * s * s) / 3.0 + c * (x * s / 2.0 + a * a / 2.0 * (x / a).asin())
}

/// The volume of `ftc_11_rb` in closed form, from the file's own geometry: a washer about
/// the z axis, faces `#95` (cylinder r 16, z −1.5…0), `#135` (cylinder r 31.5, z
/// −1.5…0), `#115` (plane z −1.5), `#163` (plane z 0.425, r 18.9385…28.5615) and two tori
/// on z = 0 of tube radius 1.5, `#187` about radius 17.5 and `#175` about radius 30, each
/// running from its cylinder over the top of its tube down to the plane at z = 0.425.
///
/// In the (r, z) half-plane the section is the rectangle 16…31.5 × −1.5…0 under the
/// curve z(r) = √(1.5² − (r − 17.5)²) for r in 16…18.9385, 0.425 on 18.9385…28.5615, and
/// √(1.5² − (r − 30)²) on 28.5615…31.5, where 18.9385 = 17.5 + √(1.5² − 0.425²) and
/// 28.5615 = 30 − √(1.5² − 0.425²) are the plane's circles `#158` and `#147`. By Pappus,
/// V = 2π ∫∫ r dr dz over the section.
fn ftc_11_rb_closed_form() -> f64 {
    let (tube, top, base) = (1.5_f64, 0.425_f64, 1.5_f64);
    let (r_in, r_out) = (16.0_f64, 31.5_f64);
    let (c_in, c_out) = (17.5_f64, 30.0_f64);
    let reach = (tube * tube - top * top).sqrt();
    let (flat_in, flat_out) = (c_in + reach, c_out - reach);
    let rectangle = base * (r_out * r_out - r_in * r_in) / 2.0;
    let inner = arc_moment(reach, tube, c_in) - arc_moment(-tube, tube, c_in);
    let flat = top * (flat_out * flat_out - flat_in * flat_in) / 2.0;
    let outer = arc_moment(tube, tube, c_out) - arc_moment(-reach, tube, c_out);
    2.0 * std::f64::consts::PI * (rectangle + inner + flat + outer)
}

/// The file's validation volume of a case, in mm³.
fn validation_volume_of(name: &str, bytes: &[u8]) -> Option<f64> {
    validation::validation_volume(&validation::parse(bytes))
        .unwrap_or_else(|error| panic!("{name}: {error}"))
}

/// One inch-unit file: a `DERIVED_UNIT` of the cubed `CONVERSION_BASED_UNIT('INCH')`.
const INCH_VOLUME: &str = "#1=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.));\
#2=LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(25.4),#1);\
#3=(CONVERSION_BASED_UNIT('INCH',#2)LENGTH_UNIT()NAMED_UNIT(*));\
#4=DERIVED_UNIT_ELEMENT(#3,3.0);#5=VOLUME_UNIT((#4));\
#6=MEASURE_REPRESENTATION_ITEM('volume measure',VOLUME_MEASURE(2.0),#5);\
#7=PRODUCT_DEFINITION_SHAPE('','',$);\
#8=PROPERTY_DEFINITION('geometric validation property','volume of part',#7);\
#9=REPRESENTATION('volume',(#6),$);#10=PROPERTY_DEFINITION_REPRESENTATION(#8,#9);\
#11=PROPERTY_DEFINITION('pmi validation property','',#7);\
#12=MEASURE_REPRESENTATION_ITEM('volume measure',VOLUME_MEASURE(9.0),#5);\
#13=REPRESENTATION('volume',(#12),$);#14=PROPERTY_DEFINITION_REPRESENTATION(#11,#13);";

/// `data` as a Part 21 file.
fn part21(data: &str) -> Vec<u8> {
    format!(
        "ISO-10303-21;HEADER;FILE_SCHEMA(('CONFIG_CONTROL_DESIGN'));ENDSEC;\
         DATA;{data}ENDSEC;END-ISO-10303-21;"
    )
    .into_bytes()
}

/// The validation-volume helper reads a geometric validation property only, and converts
/// it from the file's unit: 2 in³ is 2 × 25.4³ = 32 774.128 mm³. The `'pmi validation
/// property'` beside it is not a geometric one.
#[test]
fn validation_volume_converts_the_files_unit() {
    let file = validation::parse(&part21(INCH_VOLUME));
    let volumes = validation::validation_volumes(&file).expect("units resolve");
    assert_eq!(volumes.len(), 1, "{volumes:?}");
    assert_eq!(volumes[0].property, 8);
    assert_eq!(volumes[0].description, "volume of part");
    assert_eq!(volumes[0].written, 2.0);
    let expected = 2.0 * 25.4 * 25.4 * 25.4;
    assert!(
        relative(volumes[0].mm3, expected) <= 1e-15,
        "{} mm³, expected {expected}",
        volumes[0].mm3
    );
    assert_eq!(
        validation::validation_volume(&file).expect("one value"),
        Some(volumes[0].mm3)
    );
}

/// A file whose geometric validation volumes disagree is reported, not averaged.
#[test]
fn disagreeing_validation_volumes_are_an_error() {
    let disagreeing = INCH_VOLUME.replace(
        "'pmi validation property'",
        "'geometric validation property'",
    );
    let file = validation::parse(&part21(&disagreeing));
    assert_eq!(
        validation::validation_volumes(&file).map(|v| v.len()),
        Ok(2)
    );
    assert!(validation::validation_volume(&file).is_err());
}

/// The closed form of `ftc_11_rb` at its derived parameters (see [`ftc_11_rb_closed_form`]).
#[test]
fn ftc_11_rb_closed_form_is_the_washer() {
    let v = ftc_11_rb_closed_form();
    assert!((v - 5129.386136422674).abs() <= 1e-9 * v, "{v}");
}

/// Each file imports with the bodies, faces and closed solids of its census, and no body
/// is skipped.
#[test]
#[ignore = "needs SOLIDCORE_STEP_CORPUS"]
fn winding_face_files_match_their_census() {
    let Some(dir) = corpus_dir() else {
        eprintln!("SOLIDCORE_STEP_CORPUS is not set");
        return;
    };
    let mut failures = Vec::new();
    for (name, bodies, faces, closed, _) in CASES {
        let (bytes, report) = import(&dir, name);
        let census = validation::body_census(&validation::parse(&bytes));
        let census = (
            census.len(),
            census.iter().map(|(_, _, faces)| faces).sum::<usize>(),
        );
        if census != (bodies, faces) {
            failures.push(format!(
                "{name}: census (bodies, faces) {census:?}, table ({bodies}, {faces})"
            ));
        }
        let got = counts(&report);
        println!("{name}: (bodies, faces, closed) {got:?}");
        if got != (bodies, faces, closed) {
            failures.push(format!(
                "{name}: (bodies, faces, closed) {got:?}, expected ({bodies}, {faces}, {closed})"
            ));
        }
        if !report.skipped.is_empty() {
            failures.push(format!("{name}: skipped {:?}", report.skipped));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// The files with a geometric validation volume close at it, within
/// [`VALIDATION_TOLERANCE`]; the files without one say so, and `stc_06` writes its zero
/// placeholder.
#[test]
#[ignore = "needs SOLIDCORE_STEP_CORPUS"]
fn winding_face_files_match_their_validation_volume() {
    let Some(dir) = corpus_dir() else {
        eprintln!("SOLIDCORE_STEP_CORPUS is not set");
        return;
    };
    let mut failures = Vec::new();
    for (name, _, _, _, expected) in CASES {
        let (bytes, report) = import(&dir, name);
        let file_volume = validation_volume_of(name, &bytes);
        let measured = closed_volume(&report);
        match (expected, file_volume) {
            (Volume::Validation, Some(file_volume)) => {
                let deviation = relative(measured, file_volume);
                println!(
                    "{name}: volume {measured:.6} mm³, validation {file_volume:.6} mm³, \
                     {deviation:.2e} relative"
                );
                if deviation > VALIDATION_TOLERANCE {
                    failures.push(format!(
                        "{name}: volume {measured} mm³, validation {file_volume} mm³ \
                         ({deviation:.2e} relative)"
                    ));
                }
            }
            (Volume::Validation, None) => {
                failures.push(format!("{name}: no validation volume in the file"));
            }
            (Volume::None, Some(placeholder)) if name.starts_with("nist_stc_06") => {
                println!("{name}: volume {measured:.6} mm³, validation placeholder {placeholder}");
                if placeholder != 0.0 {
                    failures.push(format!(
                        "{name}: validation volume {placeholder}, the zero placeholder expected; \
                         assert it as a validation volume"
                    ));
                }
            }
            (Volume::None | Volume::ClosedForm(_), Some(file_volume)) => {
                failures.push(format!(
                    "{name}: carries a validation volume {file_volume} mm³; assert it"
                ));
            }
            (Volume::None | Volume::ClosedForm(_), None) => {
                println!("{name}: volume {measured:.6} mm³, no validation volume");
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// `ftc_11_rb`'s washer (see [`ftc_11_rb_closed_form`]) is one closed solid of six faces
/// today, but of the wrong volume: 3179.68 mm³ against 5129.386 mm³ in closed form (38 %
/// low). Its tori run over the top of their tubes, the faces #2580 is about.
///
/// Fails when the reader gets the volume right, which is the point: then move the file
/// to a closed-form volume test with [`CLOSED_FORM_TOLERANCE`].
#[test]
#[ignore = "needs SOLIDCORE_STEP_CORPUS"]
fn ftc_11_rb_volume_is_wrong_known_failing() {
    let Some(dir) = corpus_dir() else {
        eprintln!("SOLIDCORE_STEP_CORPUS is not set");
        return;
    };
    let name = "nist_ftc_11_asme1_rb.stp";
    let closed_form = CASES
        .iter()
        .find_map(|(n, _, _, _, v)| match v {
            Volume::ClosedForm(f) if *n == name => Some(f()),
            _ => None,
        })
        .expect("a closed-form case");
    let (_, report) = import(&dir, name);
    let measured = closed_volume(&report);
    let deviation = relative(measured, closed_form);
    println!(
        "{name}: volume {measured:.6} mm³, closed form {closed_form:.6} mm³, \
         {deviation:.2e} relative"
    );
    assert!(
        deviation.is_finite(),
        "{name}: no measurable volume ({measured} mm³), so the known failure is not the one pinned"
    );
    assert!(
        deviation > CLOSED_FORM_TOLERANCE,
        "{name} now closes at its closed form ({measured} against {closed_form}); move it \
         to a closed-form volume test"
    );
}

/// A volume that is not finite never counts as close to its expected value.
#[test]
fn relative_deviation_of_a_non_finite_volume_is_infinite() {
    assert_eq!(relative(f64::NAN, 1.0), f64::INFINITY);
    assert_eq!(relative(1.0, f64::INFINITY), f64::INFINITY);
    assert_eq!(relative(3.0, 4.0), 0.25);
}
