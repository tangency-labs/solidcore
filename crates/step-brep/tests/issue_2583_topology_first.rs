//! The cases of #2583, pinned against today's reader so that the port into the build of
//! #2578 has an acceptance test: a body built from the file's own topology
//! ([`step_brep::build_body`], `BrepConstructor`) where building each face and stitching
//! is not the way, and what the stitch-only build ([`BuildMode::StitchOnly`]) gives on
//! the same body.
//!
//! Expected values are the file's own, never the reader's output: the exporter's
//! geometric validation volume (read through `step_import`'s parser, see [`validation`],
//! shared with #2580's test), a closed form, and counts from the file's census.
//!
//! | file | bodies | faces | closed | volume expected (mm³) |
//! |---|---|---|---|---|
//! | `fixtures/rod_doubled_edge.step` (SolidCraft) | 1 | 7 | 1 | 74 900.620862 (closed form) |
//! | `AP203 geometry only/nist_ctc_02_asme1_rc` | 2 | 664 | 1 | 47 101 527.47498 (validation) |
//! | `AP203 with PMI/nist_ctc_02_asme1_ap203` | 1 | 487 | 1 | none in the file |
//! | `nist_ctc_02_asme1_ap242-e2` | 3 | 637 | 1 | 47 106 788.0321802 (validation) |
//!
//! The extra bodies are one-face `SHELL_BASED_SURFACE_MODEL`s, open sheets.
//!
//! **`rod_doubled_edge`** comes from SolidCraft's `tests/step-regressions/` (MIT OR
//! Apache-2.0; written by SolidCraft from the public commands of SolidCraft #4598): a
//! Ø40 × 60 rod with a Ø8 cross bore, whose bore rim was exported as a B-spline of 1 399
//! control points running back and forth. Its volume in closed form is
//! [`rod_closed_form`].
//!
//! **Known failing today.** The topology build of `ctc_02_ap203` leaves out face `#1152`
//! (its UV domain fails: "end angle outside [-360, 360]"), so the body is not closed and
//! the stitched build is kept ([`ctc_02_ap203_topology_build_leaves_out_a_face_known_failing`]).
//!
//! **Tolerances.** [`VALIDATION_TOLERANCE`] is #2580's (#2593): 5e-6 relative, chosen
//! from measurements, an open question to the maintainers. [`ROD_TOLERANCE`] likewise:
//! the rod's faces lie on exact surfaces, but its bore rim is the exporter's approximating
//! B-spline, so the trims only approximate the exact intersection; today's builds are
//! 7.8e-6 (topology first) and 7.0e-7 (stitch only) off, and the tolerance is the
//! smallest round value both meet.
//!
//! The corpus tests are ignored by default; the corpus is not in the repository. Fetch it
//! with `scripts/fetch-nist-pmi-corpus.sh` and run them with its directory (or one above
//! it) in `SOLIDCORE_STEP_CORPUS` (searched recursively):
//!
//! ```sh
//! SOLIDCORE_STEP_CORPUS=corpus cargo test --release -p step-brep \
//!     --test issue_2583_topology_first -- --include-ignored --nocapture
//! ```

#![allow(clippy::expect_used, clippy::indexing_slicing, clippy::panic)] // test helpers

#[path = "support/validation.rs"]
mod validation;

use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::kernel::Kernel;
use solidcore::mass::integrator::properties::{brep_volume, PropertyAttributes};
use solidcore::topology::Brep;
use step_brep::{
    build_body, build_from_topology, is_closed_solid, BuildMode, BuildOptions, BuildPath, BuiltBody,
};
use step_import::Limits;

use validation::{corpus_dir, corpus_file, relative};

/// The relative tolerance against a file's validation volume (#2593, an open question).
const VALIDATION_TOLERANCE: f64 = 5e-6;

/// The relative tolerance of the rod against its closed form (see the module docs).
const ROD_TOLERANCE: f64 = 1e-5;

/// SolidCraft's `tests/step-regressions/rod_doubled_edge.step`.
const ROD: &[u8] = include_bytes!("fixtures/rod_doubled_edge.step");

/// The volume of `rod_doubled_edge` in closed form, from the file's geometry.
///
/// The rod is the cylinder `#1713` of radius 20 about the z axis through (100, 0)
/// (its side is split into that cylinder and the rational quadratic B-spline `#1688` on
/// the same circle), from the plane z = 0 (`#1730`) to z = 60 (`#1658`): π · 20² · 60. The
/// bore is a cylinder of radius 4 about the line y = 0, z = 30 along x (the rational
/// quadratic B-spline faces `#32` and `#1202`, whose control points lie at distance 4),
/// from its flat bottom, the plane x = 110 (`#1649`), out through the rod's side at
/// x = 100 + √(400 − y²). Its volume is the integral of that depth over the bore's disk:
///
/// `∫∫ (√(400 − y²) − 10) dy dz` over `y² + (z − 30)² ≤ 16`
/// `= ∫_{−π/2}^{π/2} 32 cos²θ √(400 − 16 sin²θ) dθ − 10 · 16π` (y = 4 sin θ),
///
/// an elliptic integral, evaluated here by the trapezoidal rule over a full period of its
/// smooth periodic integrand, which converges geometrically: 32 points already agree
/// with 128 to the last bit.
fn rod_closed_form() -> f64 {
    use std::f64::consts::PI;
    let n = 128_u32;
    let step = 2.0 * PI / f64::from(n);
    let wall: f64 = (0..n)
        .map(|k| {
            let t = step * f64::from(k);
            16.0 * t.cos() * t.cos() * (400.0 - 16.0 * t.sin() * t.sin()).sqrt()
        })
        .sum::<f64>()
        * step;
    let bore = wall - 10.0 * 16.0 * PI;
    PI * 400.0 * 60.0 - bore
}

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

/// Every body of `bytes` built in `mode`.
fn build_all(bytes: &[u8], mode: BuildMode) -> Vec<BuiltBody> {
    let model = step_import::import(bytes, &Limits::default()).expect("imports");
    assert!(model.skipped.is_empty(), "skipped: {:?}", model.skipped);
    let options = BuildOptions {
        mode,
        ..BuildOptions::default()
    };
    model
        .bodies
        .iter()
        .map(|body| build_body(body, options).expect("builds"))
        .collect()
}

/// (bodies, faces, closed bodies, summed volume of the closed bodies) of `built`.
fn summary(built: &[BuiltBody]) -> (usize, usize, usize, f64) {
    let closed: Vec<&BuiltBody> = built.iter().filter(|b| is_closed_solid(&b.brep)).collect();
    (
        built.len(),
        built.iter().map(|b| b.brep.face_list().len()).sum(),
        closed.len(),
        closed
            .iter()
            .map(|b| volume(&b.brep).expect("a closed body has a volume"))
            .sum(),
    )
}

/// The census (bodies, faces) of `bytes` ([`validation::body_census`]).
fn census(bytes: &[u8]) -> (usize, usize) {
    let census = validation::body_census(&validation::parse(bytes));
    (census.len(), census.iter().map(|(_, _, faces)| faces).sum())
}

#[test]
fn rod_closed_form_is_the_bored_rod() {
    let v = rod_closed_form();
    assert!((v - 74900.62086234294).abs() <= 1e-9 * v, "{v}");
}

/// `rod_doubled_edge` builds from its topology into one closed solid of seven faces at
/// its closed-form volume.
#[test]
fn rod_doubled_edge_builds_topology_first_at_its_closed_form() {
    assert_eq!(census(ROD), (1, 7));
    let built = build_all(ROD, BuildMode::Default);
    let (bodies, faces, closed, measured) = summary(&built);
    assert_eq!((bodies, faces, closed), (1, 7, 1));
    assert_eq!(built[0].path, BuildPath::Topology);
    assert!(
        built[0].face_failures.is_empty(),
        "{:?}",
        built[0].face_failures
    );
    let expected = rod_closed_form();
    println!(
        "rod topology first: {measured:.6} mm³, closed form {expected:.6}, {:.2e} relative",
        relative(measured, expected)
    );
    assert!(
        relative(measured, expected) <= ROD_TOLERANCE,
        "volume {measured}, closed form {expected}"
    );
}

/// What the stitch-only build gives on `rod_doubled_edge` today: one closed solid of seven
/// faces at its closed-form volume, too.
#[test]
fn rod_doubled_edge_stitch_only() {
    let built = build_all(ROD, BuildMode::StitchOnly);
    let (bodies, faces, closed, measured) = summary(&built);
    assert_eq!((bodies, faces, closed), (1, 7, 1));
    assert_eq!(built[0].path, BuildPath::Stitch);
    let expected = rod_closed_form();
    println!(
        "rod stitch only: {measured:.6} mm³, closed form {expected:.6}, {:.2e} relative",
        relative(measured, expected)
    );
    assert!(
        relative(measured, expected) <= ROD_TOLERANCE,
        "volume {measured}, closed form {expected}"
    );
}

/// `ctc_02`'s two exports with a validation volume close at it, built topology first
/// (the default) and stitched only, with the bodies and faces of their census.
#[test]
#[ignore = "needs SOLIDCORE_STEP_CORPUS"]
fn ctc_02_closes_at_its_validation_volume_both_ways() {
    // (file, bodies, faces, path of the solid's default build)
    const CASES: [(&str, usize, usize, BuildPath); 2] = [
        ("nist_ctc_02_asme1_rc.stp", 2, 664, BuildPath::Topology),
        (
            "nist_ctc_02_asme1_ap242-e2.stp",
            3,
            637,
            BuildPath::Topology,
        ),
    ];
    let Some(dir) = corpus_dir() else {
        eprintln!("SOLIDCORE_STEP_CORPUS is not set");
        return;
    };
    let mut failures = Vec::new();
    for (name, bodies, faces, path) in CASES {
        let bytes = corpus_file(&dir, name);
        if census(&bytes) != (bodies, faces) {
            failures.push(format!("{name}: census {:?}", census(&bytes)));
        }
        let expected = validation::validation_volume(&validation::parse(&bytes))
            .unwrap_or_else(|error| panic!("{name}: {error}"))
            .unwrap_or_else(|| panic!("{name}: no validation volume"));
        for mode in [BuildMode::Default, BuildMode::StitchOnly] {
            let built = build_all(&bytes, mode);
            let (got_bodies, got_faces, closed, measured) = summary(&built);
            let deviation = relative(measured, expected);
            let solid_path = built
                .iter()
                .find(|b| is_closed_solid(&b.brep))
                .map(|b| b.path);
            println!(
                "{name} {mode:?}: ({got_bodies}, {got_faces}, {closed}) solid {solid_path:?}, \
                 {measured:.6} mm³, validation {expected:.6}, {deviation:.2e} relative"
            );
            if (got_bodies, got_faces, closed) != (bodies, faces, 1) {
                failures.push(format!(
                    "{name} {mode:?}: (bodies, faces, closed) ({got_bodies}, {got_faces}, \
                     {closed}), expected ({bodies}, {faces}, 1)"
                ));
            }
            let want_path = if mode == BuildMode::StitchOnly {
                BuildPath::Stitch
            } else {
                path
            };
            if solid_path != Some(want_path) {
                failures.push(format!("{name} {mode:?}: solid built {solid_path:?}"));
            }
            if deviation > VALIDATION_TOLERANCE {
                failures.push(format!(
                    "{name} {mode:?}: volume {measured}, validation {expected} \
                     ({deviation:.2e} relative)"
                ));
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// `ctc_02`'s AP203 export with PMI carries no validation volume: it closes as one solid
/// of its 487 faces, by the stitched build (see the known-failing test below).
#[test]
#[ignore = "needs SOLIDCORE_STEP_CORPUS"]
fn ctc_02_ap203_closes_by_stitching() {
    let Some(dir) = corpus_dir() else {
        eprintln!("SOLIDCORE_STEP_CORPUS is not set");
        return;
    };
    let name = "nist_ctc_02_asme1_ap203.stp";
    let bytes = corpus_file(&dir, name);
    assert_eq!(census(&bytes), (1, 487));
    assert_eq!(
        validation::validation_volume(&validation::parse(&bytes)).expect("units"),
        None
    );
    for mode in [BuildMode::Default, BuildMode::StitchOnly] {
        let built = build_all(&bytes, mode);
        let (bodies, faces, closed, measured) = summary(&built);
        println!("{name} {mode:?}: ({bodies}, {faces}, {closed}), {measured:.6} mm³");
        assert_eq!((bodies, faces, closed), (1, 487, 1), "{mode:?}");
        assert_eq!(built[0].path, BuildPath::Stitch, "{mode:?}");
    }
}

/// The topology build of `ctc_02_ap203` leaves out face `#1152`, whose UV domain fails
/// ("end angle outside [-360, 360]"), so it is not closed and [`build_body`] keeps the
/// stitched build. Fails once the topology build closes the body: then expect
/// [`BuildPath::Topology`] in [`ctc_02_ap203_closes_by_stitching`].
#[test]
#[ignore = "needs SOLIDCORE_STEP_CORPUS"]
fn ctc_02_ap203_topology_build_leaves_out_a_face_known_failing() {
    let Some(dir) = corpus_dir() else {
        eprintln!("SOLIDCORE_STEP_CORPUS is not set");
        return;
    };
    let bytes = corpus_file(&dir, "nist_ctc_02_asme1_ap203.stp");
    let model = step_import::import(&bytes, &Limits::default()).expect("imports");
    let built = build_from_topology(&model.bodies[0], BuildOptions::default()).expect("builds");
    let failed: Vec<u64> = built.face_failures.iter().map(|f| f.face).collect();
    println!(
        "topology build: {} faces, closed {}, failures {:?}",
        built.brep.face_list().len(),
        is_closed_solid(&built.brep),
        built.face_failures
    );
    assert_eq!(
        (
            failed,
            built.brep.face_list().len(),
            is_closed_solid(&built.brep)
        ),
        (vec![1152], 486, false),
        "the topology build changed; if it closes now, expect BuildPath::Topology"
    );
}
